import gzip
import os
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, coo_matrix

total_lasso = pd.read_csv("E:\BGI_Paper\E14_16\lasso\E14-16h_a_bin1\E14-16h_a_S04_bin1.gem", sep="\t")
print(total_lasso)
nucleus_lasso = pd.read_csv("E:\BGI_Paper\E14_16\lasso\E14-16h_a_bin1_nucleus_all\E14-16h_a_S04_bin1.gem", sep="\t")

with gzip.open("E:\BGI_Paper\E14_16\segmentation\E14-16h_a_S04_geneall\E14-16h_a_S04.npz.gz", "r") as f:
    mtx = coo_matrix(np.load(f))
    x = pd.Series(mtx.row) + np.min(nucleus_lasso["x"])
    y = pd.Series(mtx.col) + np.min(nucleus_lasso["y"])
    value = pd.Series(mtx.data)
    cells = pd.concat([x, y, value], axis=1)
    cells.columns = ["x", "y", "cell"]
    print(cells)

total_cells = pd.merge(total_lasso, cells, on=["x", "y"], how="inner")
print(total_cells)