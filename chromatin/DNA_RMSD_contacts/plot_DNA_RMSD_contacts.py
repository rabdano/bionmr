import matplotlib.pyplot as plt
import numpy as np
import json

# setup parameters
with open("input.json", "r") as f:
    pars = json.load(f)
first_dat_file = int(pars["first_dat_file"])
last_dat_file = int(pars["last_dat_file"])
stride = int(pars["stride"])
rmsf_fnout = pars["rmsf_fnout"]
contacts_fnout = pars["contacts_fnout"]
fig_fn = pars["fig_fn"]


# load data
position, rmsf = np.genfromtxt(rmsf_fnout,
                               delimiter=",",
                               skip_header=1,
                               usecols=(0, 1),
                               unpack=True)
contacts = np.genfromtxt(contacts_fnout,
                         delimiter=",",
                         skip_header=1,
                         usecols=1,
                         unpack=True)

# plot
fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'hspace': 0.05})
axs[0].bar(position, contacts)
axs[0].set(ylabel="Contacts", ylim=[0, np.max(contacts) * 1.05])
axs[0].text(0.03, 0.9, "A", fontsize=12, ha="center", va="center", transform=axs[0].transAxes)
axs[1].grid()
axs[1].plot(position, rmsf, marker="o", ms=4, linewidth=0.5)
axs[1].set(ylabel=r"RMSF, $\rm\AA$", ylim=[0, 5.999])
axs[1].text(0.03, 0.9, "B", fontsize=12, ha="center", va="center", transform=axs[1].transAxes)
for ax in axs:
    ax.set(xlabel="Position", xlim=[-73, 73], xticks=range(-70, 74, 10))
    ax.label_outer()
plt.savefig(fig_fn + ".pdf", bbox_inches="tight")
plt.savefig(fig_fn + ".png", bbox_inches="tight", dpi=300)
