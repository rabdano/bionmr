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
first_dat_file = int(pars["first_dat_file"])
last_dat_file = int(pars["last_dat_file"])
stride = int(pars["stride"])
reference_pdb = pars["reference_pdb"]
rmsd_fn = pars["rmsd_fnout"]
len_inner_turn = pars["len_inner_turn"]
panel_labels = pars["panel_labels"]
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
colors = ["blue", "red", "black", "green", "magenta"]
fig, axs = plt.subplots(1, 2, figsize=(32/2.54, 12/2.54), dpi=300)


#
#  plot
#
axs[0].plot(rmsd[:, 0]/1000, rmsd[:, 1], color=colors[0], linewidth=1.0, label=lab[1][:-4])
axs[0].plot(rmsd[:, 0]/1000, rmsd[:, 2], color=colors[1], linewidth=1.0, label=lab[2][:-4])
axs[0].plot(rmsd[:, 0]/1000, rmsd[:, 3], color=colors[2], linewidth=1.0, label=lab[3][:-4])
axs[1].plot(rmsd[:, 0]/1000, rmsd[:, 4], color=colors[3], linewidth=1.0, label=lab[4][:-4])
axs[1].plot(rmsd[:, 0]/1000, rmsd[:, 5], color=colors[4], linewidth=1.0, label=lab[5][:-4])

axs[0].set(xlabel="Time, ns", ylabel=r"RMSD, ${\rm\AA}$")
axs[1].set(xlabel="Time, ns", ylabel=r"RMSD, ${\rm\AA}$")
axs[0].set_xlim(axis_limits[:2])
axs[1].set_xlim(axis_limits[:2])
axs[0].set_ylim(axis_limits[2:])
axs[1].set_ylim(axis_limits[2:])
axs[0].legend()
axs[1].legend()
axs[0].text(.03, .92, panel_labels[0],
            fontsize=16,
            horizontalalignment='left',
            transform=axs[0].transAxes)
axs[1].text(.03, .92, panel_labels[1],
            fontsize=16,
            horizontalalignment='left',
            transform=axs[1].transAxes)

plt.savefig(fig_fn + ".png", bbox_inches="tight")
plt.savefig(fig_fn + ".pdf", bbox_inches="tight")
