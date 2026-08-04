import gseapy as gp
import pandas as pd
import numpy as np

tpm_df = pd.read_csv("data/tpm_matrix.csv", index_col=0)
tpm_log = np.log2(tpm_df + 1)
tpm_log = tpm_log.dropna()

ss = gp.ssgsea(data=tpm_log, 
               gene_sets='GO_Biological_Process_2025', 
               outdir='ssGSEA_GO_Results', 
               sample_norm_method='rank', 
               no_plot=False)


ss.res2d.to_csv("ssGSEA_scores.csv")

