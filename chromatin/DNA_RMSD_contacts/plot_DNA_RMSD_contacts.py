import matplotlib.pyplot as plt
import numpy as np
import json

# setup parameters
with open("input.json", "r") as f:
    pars = json.load(f)
first_dat_file = int(pars["first_dat_file"])
last_dat_file = int(pars["last_dat_file"])
stride = int(pars["stride"])
rmsd_fnout = pars["rmsd_fnout"]
contacts_1_fnout = pars["contacts_1_fnout"]
contacts_2_fnout = pars["contacts_2_fnout"]
fig_fn = pars["fig_fn"]


# load data
position, rmsd = np.genfromtxt(rmsd_fnout,
                               delimiter=",",
                               skip_header=1,
                               usecols=(0, 1),
                               unpack=True)
contacts_1 = np.genfromtxt(contacts_1_fnout,
                           delimiter=",",
                           skip_header=1,
                           usecols=1,
                           unpack=True)
contacts_2 = np.genfromtxt(contacts_2_fnout,
                           delimiter=",",
                           skip_header=1,
                           usecols=1,
                           unpack=True)

# plot
fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'hspace': 0.05})
axs[0].bar(position, contacts_1, color="b")
axs[0].bar(position, contacts_2, color="g")
axs[0].set(ylabel="Contacts", ylim=[0, np.max([np.max(contacts_1), np.max(contacts_2)]) * 1.05])
axs[0].text(0.03, 0.9, "A", fontsize=12, ha="center", va="center", transform=axs[0].transAxes)
axs[1].grid()
axs[1].plot(position, rmsd, marker="o", ms=4, linewidth=0.5, color="r")
axs[1].set(ylabel=r"RMSD, $\rm\AA$", ylim=[0, 12.5], yticks=range(13))
axs[1].text(0.03, 0.9, "B", fontsize=12, ha="center", va="center", transform=axs[1].transAxes)
for ax in axs:
    ax.set(xlabel="Position", xlim=[-73, 73], xticks=range(-70, 74, 10))
    ax.label_outer()
plt.savefig(fig_fn + ".pdf", bbox_inches="tight")
plt.savefig(fig_fn + ".png", bbox_inches="tight", dpi=300)
