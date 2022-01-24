import gzip
import os
import numpy as np
import pandas as pd
import stDrosophila as sd

slice = "E14-16h_a_S14"
total_file = f"E:\BGI_Paper\E14_16\lasso\E14-16h_a_bin1\{slice}_bin1.gem"
nucleus_file = f"E:\BGI_Paper\E14_16\lasso\E14-16h_a_bin1_nucleus_all\{slice}_bin1.gem"
cells_file = f"E:\BGI_Paper\E14_16\segmentation\{slice}_geneall\{slice}.npz.gz"
save = f"E:\BGI_Paper\E14_16\lasso\E14-16h_a_bin1_cells\{slice}_bin1.gem.gz"
total_cells = sd.tl.mapping2lasso(total_file=total_file, nucleus_file=nucleus_file, cells_file=cells_file, save=save)
print(total_cells)

