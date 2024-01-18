import scanpy as sc
from scipy import io
import os
import shutil
import gzip
import os

def adata2raw(adata, save_dir):
    os.makedirs(save_dir, exist_ok=True)
    with open(f"{save_dir}/barcodes.tsv", "w") as f:
        for item in adata.obs_names:
            f.write(item + "\n")

    with open(f"{save_dir}/features.tsv", "w") as f:
        for item in ["\t".join([x, x, "Gene Expression"]) for x in adata.var_names]:
            f.write(item + "\n")
    io.mmwrite(f"{save_dir}/matrix", adata.layers["counts"].T)
    adata.obs.to_csv(f"{save_dir}/metadata.csv")


def gzip_file(file_path, delete_original=False):
    with open(file_path, "rb") as f_in:
        with gzip.open(f"{file_path}.gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    if delete_original:
        os.unlink(file_path)