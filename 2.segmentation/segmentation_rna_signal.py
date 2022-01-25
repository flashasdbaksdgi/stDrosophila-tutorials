import os
import gzip
import cv2
import matplotlib.pyplot as plt
import numpy as np
import skimage
import sklearn
import spateo as st
from scipy import ndimage

slice = "E14-16h_a_S14"

plt.rcParams['image.interpolation'] = 'none'
save = f"/media/yao/Elements SE/BGI_Paper/E14_16/segmentation/{slice}_geneall"
if not os.path.exists(save):
    os.mkdir(save)
total = st.io.read_bgi_agg(f"/media/yao/Elements SE/BGI_Paper/E14_16/lasso/E14-16h_a_bin1/{slice}_bin1.gem", binsize=1)
nucleus = st.io.read_bgi_agg(f"/media/yao/Elements SE/BGI_Paper/E14_16/lasso/E14-16h_a_bin1_nucleus_all/{slice}_bin1.gem",
                             binsize=1)

total_mtx = total[0]
nucleus_mtx = nucleus[0]
fig, axes = plt.subplots(ncols=2, figsize=(6, 3), tight_layout=True)
axes[0].imshow(total_mtx.A, vmin=0, vmax=5)
axes[0].set_title('total')
axes[1].imshow(nucleus_mtx.A, vmin=0, vmax=5)
axes[1].set_title('nucleus')
plt.savefig(os.path.join(save, "raw_img.png"), dpi=100)
plt.close(fig)

scores = st.pp.segmentation.icell.score_pixels(
    nucleus_mtx,
    k=5,
    method='EM+BP',
    em_kwargs=dict(downsample=100000, seed=2022),
    bp_kwargs=dict(n_threads=8, k=3, square=False, p=2 / 3, q=1 / 3),
)

thresholds = skimage.filters.threshold_multiotsu(scores, classes=3)
est_nuclei_mask = st.pp.segmentation.utils.apply_threshold(scores, 7, thresholds[0])

fig, axes = plt.subplots(ncols=2, figsize=(6, 3), tight_layout=True)
axes[0].imshow(scores)
axes[0].set_title('scores')
axes[1].imshow(est_nuclei_mask)
axes[1].set_title('nuclei segmentation')
plt.savefig(os.path.join(save, "mid_img.png"), dpi=100)
plt.close(fig)

est_marker_mask = st.pp.segmentation.utils.safe_erode(
    scores, 3, square=False, min_area=100, n_iter=10, float_k=5, float_threshold=thresholds[1]
)
est_watershed = st.pp.segmentation.label.watershed(
    nucleus_mtx.A, est_nuclei_mask, est_marker_mask, 9
)

fig, axes = plt.subplots(ncols=2, figsize=(6, 3), tight_layout=True)
axes[0].imshow(est_nuclei_mask)
axes[0].set_title('markers')

axes[1].imshow(skimage.color.label2rgb(est_watershed, bg_label=0))
axes[1].set_title('final segmentation')
plt.savefig(os.path.join(save, "segmentation.png"), dpi=100)
plt.close(fig)

with gzip.open(os.path.join(save, f"{slice}.npz.gz"), "w") as f:
    np.save(f, est_watershed)
