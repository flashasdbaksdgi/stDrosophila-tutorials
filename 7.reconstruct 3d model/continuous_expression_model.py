import anndata as ad
import numpy as np

from dynamo.vectorfield.scVectorField import SparseVFC

# load data
adata = ad.read("E:\BGI_Paper\Data\B_Anno_H5ad\sct_analysis_anno_h5ad\E16-18_a_SCT_anno.h5ad")

genes = ['RpL38', 'RpL39', 'RpL41', 'RpL9', 'RpS12', 'RpS15', 'RpS20']

adata = adata[:, genes].copy()
adata.obs["x"] = adata.obs["x"] - adata.obs["x"].min()
adata.obs["y"] = adata.obs["y"] - adata.obs["y"].min()
adata.obs["z"] = adata.obs["z"] - 4.9
adata.obs[["x", "y", "z"]] = adata.obs[["x", "y", "z"]].round(2)
adata.obsm["spatial"] = adata.obs[["x", "y", "z"]].values

# SparseVFC Parameters
## X (ndarray (dimension: n_obs x n_features)) – Original data.
X = adata.obsm["spatial"]
## V (ndarray (dimension: n_obs x n_features)) – Velocities of cells in the same order and dimension of X.
V = adata[:, genes].X - adata[:, genes].X.mean(0)
## Grid (ndarray) – The function that returns diffusion matrix which can be dependent on the variables (for example, genes)
grid_num = 50
min_vec, max_vec = (
    X.min(0),
    X.max(0),
)
min_vec = min_vec - 0.01 * np.abs(max_vec - min_vec)
max_vec = max_vec + 0.01 * np.abs(max_vec - min_vec)
Grid_list = np.meshgrid(*[np.linspace(i, j, grid_num) for i, j in zip(min_vec, max_vec)])
Grid = np.array([i.flatten() for i in Grid_list]).T

# Apply sparseVFC (vector field consensus) algorithm to learn a functional form of the vector field from random samples with outlier on the entire space robustly and efficiently.
lambda_ = 0.02  # 1e-10
lstsq_method = "scipy"
res = SparseVFC(X, V, Grid, lstsq_method=lstsq_method, M=1000, lambda_=lambda_)

import matplotlib.pyplot as plt
import dynamo as dyn
from dynamo.vectorfield import vector_field_function

three_d_func = lambda x: vector_field_function(x, res)
three_d_adata = adata.copy()
three_d_adata[:, genes].X = three_d_func(X) + adata[:, genes].X.mean(0)

"""
dyn.pl.space(three_d_adata[adata.obs['z'] == 0.0], genes=[genes[0]])
plt.show()
dyn.pl.space(adata[adata.obs['z'] == 0.0], genes=[genes[0]])
"""

import stDrosophila as sd

mesh = sd.pl.build_three_d_model(adata=three_d_adata, groupby=None, group_show="all",
                                 gene_show=genes[0], voxelize=True, voxel_size=[0.5, 0.5, 0.7])
sd.pl.easy_three_d_plot(mesh=mesh, scalar="genes", save="exmaple.png")