# Author: Kayla Mac
# Date: 8/3/26
# Dataset: Plasmidsaurus NIA-GCGA
# Purpose: format counts file for further analysis


pacman::p_load(readr,
               tidyverse)


# read in counts file
tbl <- read_tsv("data/4ZJQV2-expression-matrix.tsv")


# get columns we want and write as a new file
tbl_norm <- tbl |>
  select(-c("gene_id", "gene_biotype")) |>
  select(gene_name, ends_with("_cpm"))

write_csv(tbl_norm, "data/normalized_matrix.csv")


tbl_ct <- tbl |>
  select(-c("gene_id", "gene_biotype")) |>
  select(gene_name, ends_with("_count"))

write_csv(tbl_ct, "data/counts_matrix.csv")


