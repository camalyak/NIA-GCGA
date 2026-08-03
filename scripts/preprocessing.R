# Author: Kayla Mac
# Date: 8/3/26
# Dataset: Plasmidsaurus NIA-GCGA
# Purpose: format counts file for further analysis


pacman::p_load(readr,
               tidyverse)


# read in counts file
tbl <- read_tsv("data/4ZJQV2-expression-matrix.tsv")


# split by cpm and counts
# normalized matrix
cpm_cols <- tbl |> 
  select(ends_with("_cpm")) |> 
  names()

cpm_cols_sorted <- cpm_cols[order(as.numeric(str_extract(cpm_cols, "(?<=_)\\d+(?=_cpm)")))]

tbl_norm <- tbl |>
  select(gene_name, all_of(cpm_cols_sorted))

write_csv(tbl_norm, "data/normalized_matrix.csv")


# raw counts
ct_cols <- tbl |>
  select(ends_with("_count")) |>
  names() 

ct_cols_sorted <- ct_cols[order(as.numeric(str_extract(ct_cols, "(?<=_)\\d+(?=_count)")))]

tbl_ct <- tbl |>
  select(gene_name, all_of(ct_cols_sorted))

# for DESeq2, need to convert from dbl to int
tbl_ct <- tbl_ct |> 
  mutate(across(ends_with("_count"), ~ as.integer(round(.x))))

write_csv(tbl_ct, "data/counts_matrix.csv")
