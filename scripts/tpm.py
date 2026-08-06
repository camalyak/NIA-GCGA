# Import necessary libraries
import requests
import json
import pandas as pd
import numpy as np

def get_mouse_gene_lengths(gene_list):
    server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/mouse"
    
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    lengths = {}
    
    # Process in chunks of 1000
    for i in range(0, len(gene_list), 1000):
        print(f"Processing genes {i} to {i+1000}...")
        chunk = gene_list[i:i+1000]
        data = {"symbols": chunk}
        
        response = requests.post(server + ext, headers=headers, json=data)
        
        if response.ok:
            results = response.json()
            for symbol, info in results.items():
                # Ensembl returns 'None' if a gene symbol isn't found
                if info and 'start' in info and 'end' in info:
                    lengths[symbol] = abs(info['end'] - info['start'])
        else:
            print(f"Error fetching chunk starting at {i}: {response.status_code}")
            
    return lengths


df = pd.read_csv("data/counts_matrix_tpm.csv", index_col=0)
lengths_dict = get_mouse_gene_lengths(df.index.dropna().tolist())
gene_lengths = pd.Series(lengths_dict)

# Align and calculate TPM
common_genes = df.index.intersection(gene_lengths.index)
df_filtered = df.loc[common_genes]
lengths_kb = gene_lengths.loc[common_genes] / 1000.0

rpk = df_filtered.div(lengths_kb, axis=0)
tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6

tpm.to_csv("tpm_matrix.csv")
print("TPM matrix saved to tpm_matrix.csv")