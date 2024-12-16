

import os
import tempfile
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scvi


import torch
from scvi.model.utils import mde
scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

torch.set_float32_matmul_precision("high")
save_dir = "temp"

adata = anndata.read_h5ad(filename=os.path.join("Python_Objects", "adata_preprocessed.h5ad"))  

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample")

scvi_model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)

scvi.settings.dl_num_workers = 64

scvi_model.train()

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()

SCVI_MDE_KEY = "X_scVI_mde"
adata.obsm[SCVI_MDE_KEY] = mde(adata.obsm[SCVI_LATENT_KEY])

adata.write_h5ad(filename=os.path.join("Python_Objects", "adata_scvi.h5ad"))

