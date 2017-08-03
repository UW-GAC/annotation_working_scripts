library(tidyverse)
library(stringr)
library(digest) # to generate md5 hashes as unique keys

# ----------------------------------------------------------------------------------------
# Define some variables
# ----------------------------------------------------------------------------------------

# path to input file
# this is a 10,005-line example input file (first 10k from
# /projects/topmed/downloaded_data/WGSA_annotation/freeze4/freezes_2a_3a_4_annot/freezes_2a_3a_4.chr16.snp.annotated.general20170422.gz
# plus last 5 from data/field_parsing_challenges.txt)
#test_file <- "data/6_test.gz"
#test_file <- "data/10k_test.gz"
test_file <- "/projects/topmed/downloaded_data/WGSA_annotation/freeze4/freezes_2a_3a_4_annot/freezes_2a_3a_4.chr1.snp.annotated.general20170422.gz"

# path to output file
out_file <- "/projects/topmed/variant_annotation/annotation_parsing/chr_1.csv"

# chunk size
chunk_size <- 10000

# path to mysql config file
#cnf_path <- "/home/staff/bheavner/.mysql-topmed-annotation.cnf"

# target database
#database <- "mysql_topmed_annotation_full_test"

# target table for writing
# (schema def in ~/projects/variant_work/db_build/pilot_test_schema.sql)
#db_table <- "pilot_16_test"
#db_table <- "pilot_15_test"
#db_table <- "pilot_14_test"
#db_table <- "pilot_13_test"
#db_table <- "pilot_12_test"
#db_table <- "pilot_11_test"
#db_table <- "pilot_10_test"

# fields for database import
desired_columns = c(
  "`#chr`",
  "pos",
  "ref",
  "alt",
  "rs_dbSNP147",
  "CADD_phred",
  "ENCODE_TFBS_score",
  "ENCODE_TFBS_cells",
  "ENCODE_Dnase_score",
  "ENCODE_Dnase_cells",
  "FANTOM5_enhancer_permissive",
  "FANTOM5_enhancer_robust",
  "ANNOVAR_ensembl_precedent_consequence",
  "ANNOVAR_ensembl_precedent_gene",
  "unique_variant",
  "VEP_ensembl_Consequence",
  "VEP_ensembl_Transcript_ID",
  "VEP_ensembl_Gene_Name",
  "VEP_ensembl_Gene_ID",
  "VEP_ensembl_Protein_ID",
  "VEP_ensembl_CCDS",
  "VEP_ensembl_SWISSPROT",
  "VEP_ensembl_Codon_Change_or_Distance",
  "VEP_ensembl_Amino_Acid_Change",
  "VEP_ensembl_HGVSc",
  "VEP_ensembl_HGVSp",
  "VEP_ensembl_cDNA_position",
  "VEP_ensembl_CDS_position",
  "VEP_ensembl_Protein_position",
  "VEP_ensembl_Exon_or_Intron_Rank",
  "VEP_ensembl_STRAND",
  "VEP_ensembl_CANONICAL",
  "VEP_ensembl_LoF",
  "VEP_ensembl_LoF_filter",
  "VEP_ensembl_LoF_flags",
  "VEP_ensembl_LoF_info"
)

# |-separated fields that will need to be parsed prior to import
to_split <-
  c(
    "VEP_ensembl_Consequence",
    "VEP_ensembl_Transcript_ID",
    "VEP_ensembl_Gene_Name",
    "VEP_ensembl_Gene_ID",
    "VEP_ensembl_Protein_ID",
    "VEP_ensembl_CCDS",
    "VEP_ensembl_SWISSPROT",
    "VEP_ensembl_Codon_Change_or_Distance",
    "VEP_ensembl_Amino_Acid_Change",
    "VEP_ensembl_HGVSc",
    "VEP_ensembl_HGVSp",
    "VEP_ensembl_cDNA_position",
    "VEP_ensembl_CDS_position",
    "VEP_ensembl_Protein_position",
    "VEP_ensembl_Exon_or_Intron_Rank",
    "VEP_ensembl_STRAND",
    "VEP_ensembl_CANONICAL",
    "VEP_ensembl_LoF",
    "VEP_ensembl_LoF_filter",
    "VEP_ensembl_LoF_flags",
    "VEP_ensembl_LoF_info"
  )

# ----------------------------------------------------------------------------------------
# helper functions
# ----------------------------------------------------------------------------------------

has_header <- function(file_chunk){
  any(str_sub(raw_line, 1, 4) == "#chr")
}


# ----------------------------------------------------------------------------------------
# File parsing and DB inport
# ----------------------------------------------------------------------------------------
# get connection to file for reading 
# (see http://stackoverflow.com/questions/12626637/reading-a-text-file-in-r-line-by-line)
readfile_con = gzfile(test_file, "r")

index = 0L

readr.show_progress = FALSE

while (TRUE) {
  # read a chunk
  raw_line = suppressWarnings(readLines(readfile_con, n = chunk_size))
  
  # check for header and grab if in chunk
  if (has_header(raw_line)) {
    raw_header <- raw_line[str_sub(raw_line, 1, 4) == "#chr"]
    all_fields <-
      read_tsv(paste0(raw_line, collapse = "\n"),
               #paste all lines to a single string
               col_types = cols(.default = col_character()))
    header_flag <- TRUE
  } else {
    all_fields <-
      read_tsv(paste0(raw_header, "\n", # if the header isn't in this chunk, add it
                      paste0(raw_line, collapse = "\n")),
               col_types = cols(.default = col_character())) # 1 obs of 359 variables
    header_flag <- FALSE
  }
  
  if (dim(all_fields)[1] == 0) { break } # end iteration if all_fields has 0 observations (to avoid dplyr errors relating to passing empty tibbles)
  
  # pick out the desired columns for further operation
  selected_columns <- all_fields %>%
    select_(.dots = desired_columns) %>% # select fields of interest
    select(chr = `#chr`, everything()) %>% # get rid of the # in the chr column name
    mutate(wgsa_version = "WGSA065") # add wgsa version
  
  # so now we've got a tibble with 9 observations of 37 variables.
  # Next, need to parse the VEP_ variables to unpack them -> does this do it by line?
  expanded <- selected_columns %>%
    separate_rows_(to_split, sep = "\\|")
  
  # add a hash of each line as a unique key
  # e.g. digest(paste(data.frame(letters[1:10], letters[11:20])[1,], collapse = ""), algo = "md5", serialize = FALSE)
  
  # first, combine columns by row for hashing
  lines <- expanded %>% unite(foo, everything())
  
  # add hash of each string and save resulting tibble
  parsed_lines <-
    mutate(expanded, hash = map_chr(lines$foo, function(x)
      digest(
        x, algo = "md5", serialize = FALSE
      ))) %>%
    distinct() # hopefully a minor pre-filter
  
  #TODO
  # sanity check columns
  
  # write tibble to tsv file
  if (header_flag) {
    write_tsv(parsed_lines, out_file)
  } else {
    write_tsv(parsed_lines, out_file, append = TRUE)
  }

  index <- index + 1L
  print(paste0("Chunks: ", index, " Lines: <= ", chunk_size * index, 
               " Records in current import: ", dim(parsed_lines)[1]))
}
close(readfile_con)
