# Author: Kayla Mac
# Date: 8/4/26
# Dataset: Plasmidsaurus NIA-GCGA
# Purpose: dotplots of GSEA processes where fdr <= 0.05 

pacman::p_load(openxlsx,
               tidyverse)


make_dotplot <- function(gsea_results, title="GSEA Dotplot", top_n=20, fdr_cutoff=0.05, save_path=NULL) {
  
  # Filter and calculate Ratio
  filtered_results <- gsea_results |>
    filter(fdr <= fdr_cutoff) |>
    separate(gene_ratio, into = c("hits", "total"), sep = "/", remove = FALSE, convert = TRUE) |>
    mutate(
      ratio = hits / total,
      minus_log10_fdr = -log10(ifelse(fdr == 0, 1e-10, fdr)),
      wrapped_term = stringr::str_wrap(Term, width = 35)
    ) |>
    arrange(fdr) |>
    head(top_n)
  
  p <- ggplot(filtered_results, aes(x = nes, 
                                    y = reorder(wrapped_term, nes), 
                                    size = ratio,        # Dot size based on Ratio
                                    color = minus_log10_fdr)) +
    geom_point() +
    scale_color_gradientn(colors = c("blue", "red"), name = "Adj. p-val\n(-log10)") +
    labs(title = title, 
         subtitle =  "Top 20 Pathways with FDR <= 0.05",
         x = "Normalized Enrichment Score (NES)", 
         y = NULL,
         size = "Enrichment\nRatio") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          panel.border = element_rect(color = "cornsilk4", fill = NA, linewidth = 1),
          plot.subtitle = element_text(color = "cornsilk4"))
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 7.25, height = 12, dpi = 300)
  }
  
  return(p)
}


# chow v hfd non-tumor control
gsea_file <- read.xlsx("data/9mo_stern_gsea.xlsx")

make_dotplot(gsea_file, 
             title = "Pathway Enrichment in 9 Month Agonist\nRelative to Saline", 
             top_n = 20, 
             fdr_cutoff = 0.05, 
             save_path = "images/9mo_stern_gsea.png")



gsea_file <- read.xlsx("data/18mo_stern_gsea.xlsx")

make_dotplot(gsea_file, 
             title = "Pathway Enrichment in 18 Month Agonist\nRelative to Saline", 
             top_n = 20, 
             fdr_cutoff = 0.05, 
             save_path = "images/18mo_stern_gsea.png")


gsea_file <- read.xlsx("data/saline_stern_gsea.xlsx")

make_dotplot(gsea_file, 
             title = "Pathway Enrichment in 18 Month Saline\nRelative to 9 Month", 
             top_n = 20, 
             fdr_cutoff = 0.05, 
             save_path = "images/saline_stern_gsea.png")


gsea_file <- read.xlsx("data/agonist_stern_gsea.xlsx")

make_dotplot(gsea_file, 
             title = "Pathway Enrichment in 18 Month Agonist\nRelative to 9 Month", 
             top_n = 20, 
             fdr_cutoff = 0.05, 
             save_path = "images/agonist_stern_gsea.png")
