# Author: Kayla Mac
# Date: 8/10/26
# Dataset: Plasmidsaurus NIA-GCGA
# Purpose: plot samples in PCA


tpm <- read_csv("data/tpm_matrix.csv")

tpm_mat <- tpm |>
  tibble::column_to_rownames("...1") |>
  as.matrix()

tpm_log <- log2(tpm_mat + 1)

gene_var <- apply(tpm_log, 1, var)
top_var_genes <- names(sort(gene_var, decreasing = TRUE))[1:3000]

pca_input <- t(tpm_log[top_var_genes, ])
pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)


# set up group names and labels
gs <- c(10, 10, 10, 10)
gn <- c("9_mo_S", "9_mo_A", "18_mo_S", "18_mo_A")
group_labels <- rep(gn, times = gs)


pca_df <- as.data.frame(pca_result$x) |>
  tibble::rownames_to_column("Sample") |>
  mutate(Group = group_labels)


p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("purple", "cyan", "green", "red")) +
  theme_minimal() +
  labs(title = "PCA of TPM Expression",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"))

ggsave("images/pca.png", p_pca, width = 7, height = 6, units = "in", dpi = 300)

