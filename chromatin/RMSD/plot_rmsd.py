import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from glob import glob


plt.rcParams.update({'font.size': 12})
colors = ["blue", "red", "green", "magenta"]


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


files = sorted(glob('rmsd_*.csv'))
plt.figure(figsize=cm2inch(16, 12), dpi=300)

for i, fn in enumerate(files):
    data = np.loadtxt(fn, skiprows=1, delimiter=",")
    x = data[:, 0] / 1000
    y = data[:, 1]
    plt.plot(x, y, color=colors[i], linewidth=1.0)

plt.xlabel("Time, ns")
plt.ylabel(r"RMSD, ${\rm\AA}$")

plt.savefig("RMSD.png", bbox_inches="tight")
plt.savefig("RMSD.pdf", bbox_inches="tight")
