import anndata as ad
import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib as mpl
from pandas import DataFrame
from three_d_plots import smoothing_mesh, build_three_d_model, easy_three_d_plot, three_d_slicing


adata = ad.read("/media/yao/Elements SE/BGI_Paper/Data/B_Anno_H5ad/sct_analysis_anno_h5ad/E16-18_a_SCT_anno.h5ad")
adata.obsm["spatial"] = adata.obs[["x", "y", "z"]].values
coordsby = "spatial"

voxel_model = build_three_d_model(adata=adata, coordsby=coordsby,  n_surf=10000, groupby="anno", gene_show="Adh",voxelize=False, mask_alpha=0.05, group_show="fat body")
p = pv.Plotter()
p.add_mesh(voxel_model, rgba=True, scalars="genes_rgba")
p.show()


"""
sphere_elv = surface.elevation(low_point=None, high_point=None)
sphere_elv.plot(smooth_shading=True)
print(sphere_elv["Elevation"])

"""
"""
plane = pv.Plane()
_ = surface.compute_implicit_distance(plane, inplace=True)
dist = surface['implicit_distance']
type(dist)

pl = pv.Plotter()
_ = pl.add_mesh(surface, scalars='implicit_distance', cmap='bwr')
_ = pl.add_mesh(plane, color='w', style='wireframe')
pl.show()"""
