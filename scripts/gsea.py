# Create virtual environment
# python -m venv .venv

# Activate it
# .\.venv\Scripts\Activate.ps1

# Install packages
# pip install gseapy

# Run your script
# python gsea.py

# Check if gseapy is installed
import sys
try:
    from gseapy import Msigdb
except ImportError:
    print("gseapy not installed. Run: pip install gseapy")
    sys.exit(1)


# Import necessary libraries
import matplotlib
import matplotlib.pyplot as plt
from gseapy import Msigdb
from gseapy import heatmap
from gseapy import dotplot
import pandas as pd
import gseapy as gp
import numpy as np
import seaborn as sns
import textwrap


# Function to perform GSEA analysis
def gsea_analysis(dge_file):
    df = dge_file

    df_clean = df.groupby('Symbol')[['pvalue', 'log2FoldChange']].mean().reset_index()

    df_clean["Rank"] = np.sign(df_clean.log2FoldChange) * -np.log10(df_clean.pvalue)

    ranking = df_clean[['Symbol', 'Rank']].sort_values("Rank", ascending=False)

    pre_res = gp.prerank(rnk=ranking, gene_sets='GO_Biological_Process_2025',
                          permutation_num=1000, weight=1, seed=6)

    out = []
    for term in list(pre_res.results):
        out.append([term,
                    pre_res.results[term]['pval'],
                    pre_res.results[term]['fdr'],
                    pre_res.results[term]['nes'],
                    pre_res.results[term]['es'],
                    pre_res.results[term]['tag %'],
                    pre_res.results[term]['lead_genes']])

    out_df = pd.DataFrame(out, columns=['Term', 'p_value', 'fdr', 'nes', 'es', 'gene_ratio', 'gene'])

    return out_df


dge_file = pd.read_csv("C:\\Users\\kaylamac\\Documents\\Code\\NIA-GCGA\\data\\deseq_9mo_stern.csv")
gsea_file = gsea_analysis(dge_file)
gsea_file.to_excel("C:\\Users\\kaylamac\\Documents\\Code\\NIA-GCGA\\data\\9mo_stern_gsea.xlsx", index = False)


dge_file = pd.read_csv("C:\\Users\\kaylamac\\Documents\\Code\\NIA-GCGA\\data\\deseq_18mo_stern.csv")
gsea_file = gsea_analysis(dge_file)
gsea_file.to_excel("C:\\Users\\kaylamac\\Documents\\Code\\NIA-GCGA\\data\\18mo_stern_gsea.xlsx", index = False)


dge_file = pd.read_csv("C:\\Users\\kaylamac\\Documents\\Code\\NIA-GCGA\\data\\deseq_saline_stern.csv")
gsea_file = gsea_analysis(dge_file)
gsea_file.to_excel("C:\\Users\\kaylamac\\Documents\\Code\\NIA-GCGA\\data\\saline_stern_gsea.xlsx", index = False)


dge_file = pd.read_csv("C:\\Users\\kaylamac\\Documents\\Code\\NIA-GCGA\\data\\deseq_agonist_stern.csv")
gsea_file = gsea_analysis(dge_file)
gsea_file.to_excel("C:\\Users\\kaylamac\\Documents\\Code\\NIA-GCGA\\data\\agonist_stern_gsea.xlsx", index = False)