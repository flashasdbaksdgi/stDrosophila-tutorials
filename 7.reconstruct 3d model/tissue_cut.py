import anndata as ad
import numpy as np
import pandas as pd
from sklearn.covariance import EllipticEnvelope
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
from stDrosophila.pl import build_three_d_model, easy_three_d_plot, three_d_slicing

from dynamo.vectorfield import vector_field_function
from dynamo.vectorfield.scVectorField import SparseVFC


def coords_KDE(
    coords: np.array,
    kernel: str = "gaussian",
    bandwidth: float = 1.0
):
    """Kernel density estimation based on sklearn"""

    return KernelDensity(kernel=kernel, bandwidth=bandwidth).fit(coords).score_samples(coords)


def anomaly_detection(
    coords,
    contamination=0.05
):
    """Outlier handling"""

    return EllipticEnvelope(contamination=contamination).fit(coords).predict(coords)


def coords_filter_kde(
    adata: ad.AnnData,
    coordsby: str = "spatial",
    outlier: float = 0.2,
    kernel: str = "gaussian",
    bandwidth: float = 1.0
):
    float_type = np.float64

    # Kernel density estimation
    coords = pd.DataFrame(adata.obsm[coordsby].astype(float_type))
    adata.obs["coords_kde"] = coords_KDE(coords=coords, kernel=kernel, bandwidth=bandwidth)

    # Filter outliers
    percent_outlier = f"{int(outlier*100)}%"
    CV = adata.obs["coords_kde"].describe(percentiles=[outlier])[percent_outlier]

    return adata[adata.obs["coords_kde"] > CV, :]


def coords_filter_ee(
    adata: ad.AnnData,
    coordsby: str = "spatial",
    outlier: float = 0.2,
):
    float_type = np.float64

    coords = pd.DataFrame(adata.obsm[coordsby].astype(float_type))
    adata.obs["outlier"] = anomaly_detection(coords=coords, contamination=outlier)

    return adata[adata.obs["outlier"] != -1, :]

# Example data
file = 'E:\BGI_Paper\Data\B_Anno_H5ad\sct_analysis_anno_h5ad\E16-18_a_SCT_anno.h5ad'
adata = ad.read(file)
adata.obs["x"] = adata.obs["x"] - adata.obs["x"].min()
adata.obs["y"] = adata.obs["y"] - adata.obs["y"].min()
adata.obs["z"] = adata.obs["z"] - 4.9
adata.obs[["x", "y", "z"]] = adata.obs[["x", "y", "z"]].round(2)
adata.obsm["spatial"] = adata.obs[["x", "y", "z"]].values
adata = adata[adata.obs["anno"] == "fat body", :]

adata = coords_filter_kde(adata, outlier=0.2)
adata = coords_filter_ee(adata, outlier=0.2)
"""fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(new_adata.obs["x"], new_adata.obs["y"], new_adata.obs["z"], c=new_adata.obs["coords_kde"])
plt.tight_layout()
plt.show()
"""
"""
mesh1 = build_three_d_model(adata=adata, groupby=None, n_surf=100000, voxel_smooth=200,
                            voxelize=True, voxel_size=[0.5, 0.5, 0.7])
easy_three_d_plot(mesh=mesh1, scalar="groups", save="exmaple1.png")"""

import pyvista as pv
import pyacvd
cloud = pv.PolyData(adata.obsm["spatial"])
#surf = cloud.reconstruct_surface()
surf = cloud.delaunay_3d().extract_geometry()
surf.subdivide(nsub=3, subfilter="loop", inplace=True)
clustered = pyacvd.Clustering(surf)
clustered.cluster(10000)
uniform_surf = clustered.create_mesh()
density = uniform_surf.length / 200
surface = pv.voxelize(uniform_surf, density=density, check_surface=False)

surface.plot(show_grid=True)

