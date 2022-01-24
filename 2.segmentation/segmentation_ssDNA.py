
import gzip

import cv2
import matplotlib.pyplot as plt
import numpy as np
import skimage
import sklearn
import spateo as st
from scipy import ndimage

plt.rcParams['image.interpolation'] = 'none'

tif = skimage.io.imread('/media/yao/Elements SE/BGI_Paper/mouse_diaphragm/ssDNA_image/FP200000514BR_C1.tif')

fig, axes = plt.subplots(ncols=1, figsize=(3, 3), tight_layout=True)
axes.imshow(tif)
axes.set_title('ssDNA')
plt.show()
plt.close(fig)

nuclei_mask = st.pp.segmentation.icell.mask_nuclei_from_stain(
    tif, otsu_classes=3, otsu_index=0, local_k=45, mk=7
)
fig, axes = plt.subplots(ncols=2, figsize=(6, 3), tight_layout=True)
axes[0].imshow(tif)
axes[0].set_title('nuclei')
axes[1].imshow(nuclei_mask)
axes[1].set_title('nuclei segmentation')
plt.show()
plt.close(fig)

marker_mask = st.pp.segmentation.utils.safe_erode(
    nuclei_mask, 3, square=False, min_area=100, n_iter=10
)
watershed = st.pp.segmentation.label.watershed(
    tif, nuclei_mask, marker_mask, 5
)

fig, axes = plt.subplots(ncols=2, figsize=(6, 3), tight_layout=True)
axes[0].imshow(nuclei_mask)
axes[0].imshow(marker_mask, alpha=0.5)
axes[0].set_title('markers')

axes[1].imshow(skimage.color.label2rgb(watershed, bg_label=0))
axes[1].set_title('final segmentation')
plt.show()
plt.close(fig)