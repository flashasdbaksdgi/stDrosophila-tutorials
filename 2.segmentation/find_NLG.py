
import stDrosophila as sd

import os
import platform

folder = "/media/yao/Elements SE/BGI_Paper/E14_16/lasso/E14-16h_a_bin1"
opath = "/media/yao/Elements SE/BGI_Paper/E14_16/lasso/E14-16h_a_bin1_nucleus_all"
# ----------------------------------------------------------------------------------------------------------------------
files = [os.path.join(root, filename) for root, dirs, files in os.walk(folder) for filename in files]
files.sort()
for file in files[14:]:
    filename = file.split("\\")[-1] if platform.system() == "Windows" else file.split("/")[-1]
    print(filename)
    nucleus_lasso = sd.tl.find_nuclear_genes(file=file, save=os.path.join(opath, filename), gene_num="all")
    print(nucleus_lasso)


