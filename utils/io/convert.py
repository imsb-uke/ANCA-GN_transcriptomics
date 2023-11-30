import os
import pandas as pd
import scanpy as sc
import argparse
from rpy2 import robjects

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

data = sc.read(args.filename)
data.X = data.raw.to_adata()[:, data.var_names].X
counts = data.to_df()
metadata = data.obs

metadata_tmp, counts_tmp = args.filename.split(".h5ad")[0]+"_metadata.csv", args.filename.split(".h5ad")[0]+"_counts.csv"
print(f"Saving tmp files to {metadata_tmp} and {counts_tmp}.")

counts.to_csv(counts_tmp)
metadata.to_csv(metadata_tmp)