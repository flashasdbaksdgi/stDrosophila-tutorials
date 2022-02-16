import anndata as ad
import numpy as np
import pandas as pd
import pyacvd
import pyvista as pv
import stDrosophila as sd

from anndata import AnnData
from pandas import DataFrame
from pyvista import PolyData
from typing import Tuple

# Example data
file = '/media/yao/Elements SE/BGI_Paper/Data/B_Anno_H5ad/sct_analysis_anno_h5ad/E16-18_a_SCT_anno.h5ad'
adata = ad.read(file)
adata.obs["x"] = adata.obs["x"] - adata.obs["x"].min()
adata.obs["y"] = adata.obs["y"] - adata.obs["y"].min()
adata.obs["z"] = adata.obs["z"] - 4.9
adata.obs[["x", "y", "z"]] = adata.obs[["x", "y", "z"]].round(2)
adata.obsm["spatial"] = adata.obs[["x", "y", "z"]].values
arr = adata.obs[["x", "y", "z"]].values

clusters = ["fat body", "CNS"]
clusters_color = ["#F56867", "#FEB915"]
sub_meshes = []
for cluster, color in zip(clusters, clusters_color):
    sub_adata = adata[adata.obs["anno"] == cluster, :].copy()
    sub_adata = sd.tl.om_kde(sub_adata, percent=0.2)
    sub_adata = sd.tl.om_EllipticEnvelope(sub_adata, percent=0.2)
    sub_mesh = sd.tl.build_three_d_model(adata=sub_adata, groupby='anno', group_show=cluster, group_cmap=[color],
                                             mask_alpha=0, surf_alpha=0.5, surf_color=color,
                                             n_surf=10000, voxel_smooth=200, voxelize=True, voxel_size=[0.5, 0.5, 0.7])
    sub_meshes.append(sub_mesh)


mesh = sd.tl.build_three_d_model(adata=adata, groupby='anno', group_show='all', surf_alpha= 0.2, group_amap=0,
                                 n_surf=10000, voxel_smooth=200, voxelize=True, voxel_size=[0.5, 0.5, 0.7])


p = pv.Plotter()
p.add_mesh(mesh, rgba=True, scalars="groups_rgba")
for sub_mesh in sub_meshes:
    p.add_mesh(sub_mesh, rgba=True, scalars="groups_rgba")
p.show(screenshot="ex.png")
