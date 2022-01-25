
import os
import stDrosophila as sd
import pandas as pd
import numpy as np

slices = ["E14-16h_a_S01","E14-16h_a_S02","E14-16h_a_S04","E14-16h_a_S05","E14-16h_a_S06",
          "E14-16h_a_S07","E14-16h_a_S08","E14-16h_a_S09","E14-16h_a_S10","E14-16h_a_S12",
          "E14-16h_a_S13","E14-16h_a_S14","E14-16h_a_S15","E14-16h_a_S16"]

adata_list = []
for slice in slices:
    lasso_file = f"/media/yao/Elements SE/BGI_Paper/E14_16/lasso/E14-16h_a_bin1_cells/{slice}_bin1.gem.gz"
    label_mtx_file = f"/media/yao/Elements SE/BGI_Paper/E14_16/segmentation/{slice}_geneall/{slice}.npz.gz"
    lasso = sd.io.read_lasso(lasso_file)
    adata = sd.io.lasso2adata(data=lasso, slice=slice, label_path=label_mtx_file, z=int(slice[-2:]), z_gap=7, cellbin=False)
    print(adata.obs)
    adata_list.append(adata)

integrate_adata = adata_list[0].concatenate(adata_list[1:], batch_key="slice", batch_categories=slices, join='outer')
integrate_adata.obs["cluster"] = "no_cluster"
integrate_adata.obsm["spatial"] = integrate_adata.obs[["x", "y", "z"]].values
mesh, surface = sd.pl.build_three_d_model(adata=integrate_adata, groupby="cluster", group_show="all", group_cmap="rainbow", voxelize=True, voxel_size=[1,1,1])
sd.pl.easy_three_d_plot(mesh, None, scalar="groups", surface_color="gainsboro",surface_opacity=0.5,
                        save=r"points_all_groups.png")
# along x-axis
slices_x = sd.pl.three_d_slicing(mesh, axis='x', n_slices=7)
# along y-axis
slices_y = sd.pl.three_d_slicing(mesh, axis='y', n_slices=7)
# along z-axis
slices_z = sd.pl.three_d_slicing(mesh, axis='z', n_slices=5)
# orthogonal slicing
slices_o = sd.pl.three_d_slicing(mesh, n_slices='orthogonal')
import pyvista as pv
# display
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
p.show(screenshot=r"slicing.png")