import typing as tp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

import anndata
import scanpy as sc
from imc.types import DataFrame, Path
from imc.graphics import rasterize_scanpy

# remove batch using combat
# refer to https://github.com/ElementoLab/imc/imc/ops/clustering.py

z_score = True
remove_batch = True





# merge all the csv file

# convert to AnnaData




# Scaling/Normalization
print("Performing data scaling/normalization.")
z_score_cap: float = 3.0
if z_score:
    _ads = list()
    for roi_name in a.obs["roi"].unique():
        a2 = a[a.obs["roi"] == roi_name, :].copy()
        sc.pp.scale(a2, max_value=z_score_cap)
        a2.X[a2.X < -z_score_cap] = -z_score_cap
        # print(a2.X.min(), a2.X.max())
        _ads.append(a2)
    a = anndata.concat(_ads)
    sc.pp.scale(a)
if remove_batch:
    sc.pp.combat(a, batch_variable)
    sc.pp.scale(a)



















