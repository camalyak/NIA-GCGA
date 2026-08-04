# Author: Kayla Mac
# Date: 8/3/26
# Dataset: Plasmidsaurus NIA-GCGA
# Purpose: differential gene analysis using DESeq2


pacman::p_load(DESeq2)


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
groups <- make.names(c("9_mo_S","9_mo_A"))
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
  write.csv("data/DESeq_9mo_stern.csv", row.names = FALSE)



