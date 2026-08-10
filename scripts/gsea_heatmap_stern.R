# Author: Kayla Mac
# Date: 8/4/26
# Dataset: Plasmidsaurus NIA-GCGA
# Purpose: heatmap of GSEA processes where fdr <= 0.05 

pacman::p_load(ComplexHeatmap,
               tidyr,
               dplyr)


# load normalized matrix
ssgsea <- read_csv("data/ssGSEA_scores.csv")


# set up group names and labels
gs <- c(10, 10, 10, 10)
gn <- c("9_mo_S", "9_mo_A", "18_mo_S", "18_mo_A")
group_labels <- rep(gn, times = gs)


# convert to wide format
ssgsea_wide <- ssgsea |>
  dplyr::select(Name, Term, NES) |>
  pivot_wider(names_from = Name, values_from = NES)


# put samples in ascending order
sample_order <- ssgsea_wide |>
  dplyr::select(-Term) |>
  names() |>
  (\(x) x[order(as.numeric(gsub("\\D", "", x)))])()


# labeling for groups in heatmap
group_colors <- setNames(rainbow(length(gn)), gn)
top_anno <- HeatmapAnnotation(
  Group = factor(group_labels, levels = gn),
  col = list(Group = group_colors),
  show_annotation_name = FALSE
)


# convert to matrix
ssgsea_mat <- ssgsea_wide |>
  column_to_rownames("Term") |>
  dplyr::select(all_of(sample_order)) |>
  as.matrix()


# get z-score
ssgsea_zscore <- t(scale(t(ssgsea_mat)))


# get mean of each group for each pathway
group_means <- sapply(gn, function(g) {
  rowMeans(ssgsea_zscore[, group_labels == g], na.rm = TRUE)
})


# calculate variance between groups means
pathway_between_var <- apply(group_means, 1, var, na.rm = TRUE)


# get top 20 pathways
top20_pathways <- names(sort(pathway_between_var, decreasing = TRUE))[1:20]
mat_top20 <- ssgsea_zscore[top20_pathways, ]


# shorten the pathway string
wrapped_names <- sapply(rownames(mat_top20), function(x) {
  paste(strwrap(x, width = 40), collapse = "\n")
})
rownames(mat_top20) <- wrapped_names


p_var <- Heatmap(mat_top20,
                 name = "NES",
                 top_annotation = top_anno,
                 column_title = "Top 20 GO Pathways by Variance",
                 border = TRUE, 
                 show_row_names = TRUE, 
                 cluster_columns = FALSE,
                 cluster_rows = TRUE,
                 row_names_gp = gpar(fontsize = 6), 
                 column_names_gp = gpar(fontsize = 8))

png("images/ssgsea_htmp_var.png", width = 12, height = 8, units = "in", res = 300)
draw(p_var)
dev.off()






