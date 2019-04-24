import numpy as np
from math import ceil
from glob import glob
import os

ratio = ___KEEP_ACF_FRACTION___
input_dir = "cor_NH____FIRST_DAT_FILE___-___LAST_DAT_FILE___"
output_dir = "cor_NH____FIRST_DAT_FILE___-___LAST_DAT_FILE____diluted"
scaling = ___SCALING___
stride = ___DILUTE_STRIDE___
remove_first_point = True

files = sorted(glob(input_dir+"/*.cor"))
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for fn in files:
    data = np.loadtxt(fn)
    length = len(data[:,0])
    nlim = length * ratio
    grid = []
    if not remove_first_point:
        grid.append(0)
    # if (scaling==1.0):
    #     tau = stride
    # else:
    #     tau = 1.0
    tau = 1.0
    while tau <= nlim:
        grid.append(int(tau))
        if scaling == 1.0:
            tau += stride
        else:
            tau = ceil(tau * scaling)
    np.savetxt(output_dir+fn[len(input_dir):], data[grid, :], fmt="%14.6e")
