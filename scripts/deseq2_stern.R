# Author: Kayla Mac
# Date: 8/3/26
# Dataset: Plasmidsaurus NIA-GCGA
# Purpose: differential gene analysis using DESeq2


pacman::p_load(DESeq2)


# load counts
path <- ("data/counts_matrix.csv")
counts <- as.matrix(data.table::fread(path, header = T, colClasses = "integer"), rownames = "gene_name")


# load annotations
apath <- ("data/annot_mouse.csv")
annot <- data.table::fread(apath, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)
rownames(annot) <- annot$GeneID
