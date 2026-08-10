# Author: Kayla Mac
# Date: 8/10/26
# Dataset: Plasmidsaurus NIA-GCGA
# Purpose: plot genes from TPM matrix

pacman::p_load(ComplexHeatmap,
               tidyr,
               dplyr)


# load tpm
tpm <- read_csv("data/tpm_matrix.csv")


tpm_mat <- tpm |>
  tibble::column_to_rownames("...1") |>  # adjust to your actual gene ID column name
  as.matrix()


# set up group names and labels
gs <- c(10, 10, 10, 10)
gn <- c("9_mo_S", "9_mo_A", "18_mo_S", "18_mo_A")
group_labels <- rep(gn, times = gs)


# labeling for groups in heatmap
group_colors <- setNames(rainbow(length(gn)), gn)
top_anno <- HeatmapAnnotation(
  Group = factor(group_labels, levels = gn),
  col = list(Group = group_colors),
  show_annotation_name = FALSE
)


tpm_zscore <- t(scale(t(tpm_mat)))


group_means <- sapply(gn, function(g) {
  rowMeans(tpm_zscore[, group_labels == g], na.rm = TRUE)
})

gene_between_var <- apply(group_means, 1, var, na.rm = TRUE)

top20_genes <- names(sort(gene_between_var, decreasing = TRUE))[1:20]
mat_top20_genes <- tpm_zscore[top20_genes, ]


p_var <- Heatmap(mat_top20_genes,
        name = "Z-score",
        top_annotation = top_anno,
        column_title = "Top 20 Genes by Between-Group Variance",
        border = TRUE, 
        show_row_names = TRUE, 
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        row_names_gp = gpar(fontsize = 6), 
        column_names_gp = gpar(fontsize = 8))

png("images/tpm_htmp_var.png", width = 12, height = 8, units = "in", res = 300)
draw(p_var)
dev.off()


