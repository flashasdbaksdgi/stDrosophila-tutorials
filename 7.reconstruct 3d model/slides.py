import anndata as ad
import numpy as np
from stDrosophila.pl import build_three_d_model, easy_three_d_plot, three_d_slicing

from dynamo.vectorfield import vector_field_function
from dynamo.vectorfield.scVectorField import SparseVFC

# Example data
file = 'E:\BGI_Paper\Data\B_Anno_H5ad\sct_analysis_anno_h5ad\E16-18_a_SCT_anno.h5ad'
adata = ad.read(file)
genes = ['RpL38', 'RpL39', 'RpL41', 'RpL9', 'RpS12', 'RpS15', 'RpS20']
adata = adata[:, genes].copy()
adata.obs["x"] = adata.obs["x"] - adata.obs["x"].min()
adata.obs["y"] = adata.obs["y"] - adata.obs["y"].min()
adata.obs["z"] = adata.obs["z"] - 4.9
adata.obs[["x", "y", "z"]] = adata.obs[["x", "y", "z"]].round(2)
adata.obsm["spatial"] = adata.obs[["x", "y", "z"]].values
"""
# plot all groups and one gene(RpL38) exp
mesh1 = build_three_d_model(adata=adata, groupby='anno', group_show='all', gene_show='RpL38',
                            mask_alpha=0, n_surf=20000, voxel_smooth=500,
                            voxelize=True, voxel_size=[0.5, 0.5, 0.7])
easy_three_d_plot(mesh=mesh1, scalar="groups", save="exmaple1.png")
easy_three_d_plot(mesh=mesh1, scalar="genes", save="exmaple2.png")
"""
# plot one groups and one gene(RpL38) exp
mesh2 = build_three_d_model(adata=adata, groupby='anno', group_show='fat body', gene_show='RpL38',
                            mask_alpha=0, n_surf=20000, voxel_smooth=500,
                            voxelize=True, voxel_size=[0.5, 0.5, 0.7])
easy_three_d_plot(mesh=mesh2, scalar="groups", save="exmaple3.png")
easy_three_d_plot(mesh=mesh2, scalar="genes", save="exmaple4.png")

# continuous gene exp
X = adata.obsm["spatial"]
V = adata[:, genes].X - adata[:, genes].X.mean(0)
grid_num = 50
min_vec, max_vec = (
    X.min(0),
    X.max(0),
)
min_vec = min_vec - 0.01 * np.abs(max_vec - min_vec)
max_vec = max_vec + 0.01 * np.abs(max_vec - min_vec)
Grid_list = np.meshgrid(*[np.linspace(i, j, grid_num) for i, j in zip(min_vec, max_vec)])
Grid = np.array([i.flatten() for i in Grid_list]).T
lambda_ = 0.02  # 1e-10
lstsq_method = "scipy"
res = SparseVFC(X, V, Grid, lstsq_method=lstsq_method, M=1000, lambda_=lambda_)

three_d_func = lambda x: vector_field_function(x, res)
three_d_adata = adata.copy()
three_d_adata[:, genes].X = three_d_func(X) + adata[:, genes].X.mean(0)

mesh3 = build_three_d_model(adata=three_d_adata, groupby=None, group_show="all",gene_show=genes[0],
                            mask_alpha=0, n_surf=20000, voxel_smooth=500,
                            voxelize=True, voxel_size=[0.5, 0.5, 0.7])
easy_three_d_plot(mesh=mesh3, scalar="genes", save="exmaple5.png")

# slicing
slices_x = three_d_slicing(mesh3, axis='x', n_slices=7)
slices_y = three_d_slicing(mesh3, axis='y', n_slices=7)
slices_z = three_d_slicing(mesh3, axis='z', n_slices=5)
slices_o = three_d_slicing(mesh3, n_slices='orthogonal')
import pyvista as pv
p = pv.Plotter(shape=(2,2))
p.subplot(0, 0)
p.add_mesh(slices_x, scalars="genes_rgba", rgba=True)
p.background_color = "white"
p.camera_position = "iso"
p.subplot(0, 1)
p.add_mesh(slices_y, scalars="genes_rgba", rgba=True)
p.background_color = "white"
p.camera_position = "iso"
p.subplot(1, 0)
p.add_mesh(slices_z, scalars="genes_rgba", rgba=True)
p.background_color = "white"
p.camera_position = "iso"
p.subplot(1, 1)
p.add_mesh(slices_o, scalars="genes_rgba", rgba=True)
p.background_color = "white"
p.camera_position = "iso"
p.show(screenshot=r"exmaple6.png")

"""
# compute vol
volume_size_all = sp.tl.compute_volume(mesh=mesh3, group_show="all")
volume_size_tissue = sp.tl.compute_volume(mesh=mesh3, group_show="fat body")
"""


