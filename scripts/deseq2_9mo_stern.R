# Author: Kayla Mac
# Date: 8/3/26
# Dataset: Plasmidsaurus NIA-GCGA
# Purpose: differential gene analysis using DESeq2 - 9 month saline v agonist


pacman::p_load(DESeq2,
               ggrepel,
               apeglm)


# load counts
path <- ("data/counts_matrix.csv")
counts <- as.matrix(data.table::fread(path, header = T, colClasses = "integer"), rownames = "GeneID")


# load annotations
apath <- ("data/annot_mouse.csv")
annot <- data.table::fread(apath, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)
rownames(annot) <- annot$GeneID


# sample selection: 
sm <- "11111111110000000000XXXXXXXXXXXXXXXXXXXX"
sml <- strsplit(sm, split="")[[1]]


# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
counts <- counts[ ,sel]



# group membership for samples
gs <- factor(sml)
groups <- make.names(c("9_mo_A","9_mo_S"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(counts))


# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( counts >= 10 ) >= min(table(sml))
counts <- counts[keep, ]


ds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~Group)

ds <- DESeq(ds, test = "Wald", sfType = "poscount")


# extract results for top genes table
r <- results(ds, contrast = c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")


# save dge file as csv
r |>
  as.data.frame() |>
  tibble::rownames_to_column("GeneID") |>
  merge(annot[, c("GeneID", "Symbol", "Description")], by = "GeneID", sort = FALSE) |>
  write.csv("data/deseq_9mo_stern.csv", row.names = FALSE)

##########################################################

# plot dge volcano plot

make_volcano_fc <- function(result_table, expt_name, title, p_cutoff = 0.05, log2FoldChange_cutoff = 1) {
  # Standardize column names
  result_table$gene <- result_table$Symbol
  result_table$negLogP <- -log10(result_table$padj)
  result_table$enrichment <- "Not Significant"
  
  # Calculate counts of significant genes
  up_count <- sum(result_table$log2FoldChange >= log2FoldChange_cutoff & result_table$padj < p_cutoff, na.rm = TRUE)
  down_count <- sum(result_table$log2FoldChange <= -log2FoldChange_cutoff & result_table$padj < p_cutoff, na.rm = TRUE)
  
  # Create clear labels with dynamic gene counts
  label_expt <- paste0("Upregulated in ", expt_name, " (n = ", up_count, ")")
  label_ctrl <- paste0("Downregulated in ", expt_name, " (n = ", down_count, ")")
  
  # Positive log2FoldChange = Experimental > Control (Right side)
  result_table$enrichment[result_table$log2FoldChange >= log2FoldChange_cutoff & result_table$padj < p_cutoff] <- label_expt
  # Negative log2FoldChange = Control > Experimental (Left side)
  result_table$enrichment[result_table$log2FoldChange <= -log2FoldChange_cutoff & result_table$padj < p_cutoff] <- label_ctrl
  
  valid_genes <- result_table[!is.na(result_table$gene), ]
  
  expt_subset <- valid_genes[valid_genes$enrichment == label_expt, ]
  ctrl_subset <- valid_genes[valid_genes$enrichment == label_ctrl, ]
  
  top_expt <- head(expt_subset[order(expt_subset$log2FoldChange, decreasing = TRUE), ], 10)
  top_ctrl <- head(ctrl_subset[order(ctrl_subset$log2FoldChange, decreasing = FALSE), ], 10)
  
  # Bind datasets cleanly
  plot_labels <- rbind(top_expt, top_ctrl)
  
  # Plotting
  p <- ggplot(result_table, aes(x = log2FoldChange, y = negLogP, color = enrichment)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = setNames(
      c("#ac0f0f", "#2248b9", "grey80"),
      c(label_expt, label_ctrl, "Not Significant")
    ), na.translate = FALSE) + 
    geom_vline(xintercept = c(-log2FoldChange_cutoff, log2FoldChange_cutoff), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", alpha = 0.5) +
    theme_minimal() +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = NULL
    ) +
    theme(panel.border = element_rect(color = "cornsilk4", 
                                      fill = NA, 
                                      linewidth = 1),
          legend.position = "top")
  
  # Safe check for row count 
  if (!is.null(plot_labels) && nrow(plot_labels) > 0) {
    p <- p + geom_text_repel(data = plot_labels, aes(label = gene), size = 3, max.overlaps = 10)
  }
  
  return(p)
}


dge <- read_csv("data/deseq_9mo_stern.csv")
png("images/9_mo_volcano.png", width = 7, height = 6, units = "in", res = 300)
make_volcano_fc(dge, "9 mo agonist", "Up / Downregulation in Agonist Relative to Saline (9 months)")
dev.off()












