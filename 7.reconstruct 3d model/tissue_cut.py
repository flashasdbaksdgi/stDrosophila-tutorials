import anndata as ad
import stDrosophila as sd

# Example data
file = '/media/yao/Elements SE/BGI_Paper/Data/B_Anno_H5ad/sct_analysis_anno_h5ad/E16-18_a_SCT_anno.h5ad'
adata = ad.read(file)
adata.obs["x"] = adata.obs["x"] - adata.obs["x"].min()
adata.obs["y"] = adata.obs["y"] - adata.obs["y"].min()
adata.obs["z"] = adata.obs["z"] - 4.9
adata.obs[["x", "y", "z"]] = adata.obs[["x", "y", "z"]].round(2)
adata.obsm["spatial"] = adata.obs[["x", "y", "z"]].values
adata = adata[adata.obs["anno"] == "fat body", :]

adata = sd.tl.om_kde(adata, percent=0.2)
adata = sd.tl.om_EllipticEnvelope(adata, percent=0.05)

"""
mesh1 = sd.tl.build_three_d_model(adata=adata, groupby=None, n_surf=100000, voxel_smooth=200,
                            voxelize=True, voxel_size=[0.5, 0.5, 0.7])
sd.pl.easy_three_d_plot(mesh=mesh1, scalar="groups", save="exmaple1.png")
"""

