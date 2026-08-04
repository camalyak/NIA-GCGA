# Author: Kayla Mac
# Date: 8/3/26
# Dataset: Plasmidsaurus NIA-GCGA
# Purpose: format counts file for further analysis


pacman::p_load(readr,
               tidyverse,
               org.Mm.eg.db,
               AnnotationDbi)


# read in counts file
tbl <- read_tsv("data/4ZJQV2-expression-matrix.tsv")
symbols <- tbl$gene_name


# add GeneID and make rowname = GeneID
gene_ids <- mapIds(org.Mm.eg.db,
                   keys = symbols,
                   column = "ENTREZID",
                   keytype = "SYMBOL",
                   multiVals = "first")

tbl$GeneID <- gene_ids



# raw counts
ct_cols <- tbl |>
  dplyr::select(ends_with("_count")) |>
  names() 

ct_cols_sorted <- ct_cols[order(as.numeric(str_extract(ct_cols, "(?<=_)\\d+(?=_count)")))]

tbl_ct <- tbl |>
  dplyr::select(GeneID, all_of(ct_cols_sorted))


# for DESeq2, need to convert from dbl to int
tbl_ct <- tbl_ct |>
  filter(!is.na(GeneID)) |>
  group_by(GeneID) |>
  summarise(across(everything(), sum)) |>
  column_to_rownames("GeneID") |>
  as.matrix() |>
  round()

storage.mode(tbl_ct) <- "integer"


tbl_ct |>
  as.data.frame() |>
  tibble::rownames_to_column("GeneID") |>
  write.csv("counts_matrix.csv", row.names = FALSE)
