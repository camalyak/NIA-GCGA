# Author: Kayla Mac
# Date: 8/3/26
# Dataset: Plasmidsaurus NIA-GCGA
# Purpose: format counts file for further analysis


pacman::p_load(readr,
               tidyverse)


# read in counts file
tbl <- read_tsv("data/4ZJQV2-expression-matrix.tsv")


# get columns we want and write as a new file
new_tbl <- tbl |>
  select(-c("gene_id", "gene_biotype"))
