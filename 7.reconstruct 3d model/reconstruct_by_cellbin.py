

import stDrosophila as sd

lasso = sd.io.read_lasso("E:\BGI_Paper\E14_16\lasso\E14-16h_a_bin1_cells\E14-16h_a_S05_bin1.gem.gz")
adata = sd.io.lasso2adata(data=lasso, slice="E14-16h_a_S05", physical_coords=False)
print(adata)
print(adata.obs)



