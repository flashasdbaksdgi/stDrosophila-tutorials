## Enter raw anndata data(slices) and output path
folder = r"/media/yao/Elements SE/BGI_Paper/ssDNA_test/A1/A1_sample/center_right/result/adata"
opath = r"/media/yao/Elements SE/BGI_Paper/ssDNA_test/A1/A1_sample/center_right/result/align"
slice_col="slice"
numItermax = 200
numItermaxEmd = 100000
# ----------------------------------------------------------------------------------------------------------------------
import os
import anndata as ad
import stDrosophila as sd
import torch
files = [os.path.join(root, filename) for root, dirs1, files in os.walk(folder) for filename in files]
files.sort()
slices = [ad.read(file) for file in files]
align_slices = sd.tl.slice_alignment(slices, alpha=0.1, numItermax=numItermax, numItermaxEmd=numItermaxEmd,
                                     device=torch.device("cuda:0"), verbose=True)
if not os.path.exists(opath):
    os.mkdir(opath)
for slice in align_slices:
    subSave = os.path.join(opath, f"{slice.obs[slice_col][0]}.h5ad")
    slice.write_h5ad(subSave)

opath = "A1_center_right_align.png"
sd.pl.spatial_plot(adata=align_slices, cluster_col=None, save=opath, slice_col=slice_col, spot_size=5)

