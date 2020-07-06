import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json


#
#  setup trajectory parameters
#
with open("input.json", "r") as f:
    pars = json.load(f)
rmsd_fn = pars["rmsd_fnout"]
axis_limits = pars["axis_limits"]
fig_fn = pars["fig_fn"]


#
#  read RMSD data
#
rmsd = np.loadtxt(rmsd_fn, skiprows=1, delimiter=",")
with open(rmsd_fn, "r") as f:
    header = f.readline()
lab = header.split(",")


#
#  plot parameters
#
plt.rcParams.update({'font.size': 12})
colors = ["red", "blue", "green"]
plt.figure(figsize=(16/2.54, 12/2.54))


#
#  plot
#
plt.plot(rmsd[:, 0]/1000, rmsd[:, 1], color=colors[0], linewidth=1.0, label=lab[1][:-11])
plt.plot(rmsd[:, 0]/1000, rmsd[:, 2], color=colors[1], linewidth=1.0, label=lab[2][:-11])
plt.plot(rmsd[:, 0]/1000, rmsd[:, 3], color=colors[2], linewidth=1.0, label=lab[3][:-11])

plt.xlabel("Time, ns")
plt.ylabel(r"RMSD, ${\rm\AA}$")
plt.axis(axis_limits)

plt.legend()

plt.savefig(fig_fn + ".pdf", bbox_inches="tight")
plt.savefig(fig_fn + ".png", bbox_inches="tight")
