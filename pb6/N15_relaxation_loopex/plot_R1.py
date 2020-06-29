import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np

plot_out_name = 'R1_full'

# read data
with open("R1_results.txt", "r") as f:
    lines = f.readlines()
n_resids = len(lines)
resids = np.zeros(n_resids, dtype=int)
resnames = []
R1 = np.zeros(n_resids)
for i, line in enumerate(lines):
    resids[i] = int(line.split()[0].split('-')[0])
    resnames.append(line.split()[0].split('-')[1])
    R1[i] = float(line[:-1].split()[1])
# print(resids, resnames, R1)

plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(16.0/2.54, 12.0/2.54))

body_resids = resids[resids < 376]
body_R1 = R1[resids < 376]
ig_like_resids = resids[resids >= 376]
ig_like_R1 = R1[resids >= 376]


# Experimetnal data
# pb6_exp_resids = np.array([3, 7, 9, 10, 11, 14, 15])
# pb6_exp_R1 = np.array([1.26765865, 1.291394278, 1.166972683, 1.213193526, 1.17605744, 1.228871478, 1.201527367])

# plot
fig, ax = plt.subplots()
ax.plot(body_resids, body_R1, marker="D", ms=3, markeredgecolor="b", markerfacecolor="b",
        linewidth=1.0, color="b", label="tube body")
ax.plot(ig_like_resids, ig_like_R1, marker="D", ms=3, markeredgecolor="g", markerfacecolor="g",
        linewidth=1.0, color="g", label="Ig-like")
# ax.plot(pb6_exp_resids, pb6_exp_R1, marker="o", ms=7, markeredgecolor="r", markerfacecolor="r",
#         linestyle='None', label="Experiment")

plt.ylabel(r'${\rm R_{1},\ s^{-1}}$')
plt.xlabel('Residue')
plt.axis([0, 465, 0, 1.6])
plt.xticks([1] + list(range(50, 464, 50)))
plt.yticks(np.arange(0, 1.61, 0.2))

# Put a legend above current axis
ax.legend(loc='lower left', bbox_to_anchor=(0.01, 1.01), ncol=2,
          borderaxespad=0, frameon=False, numpoints=1)

plt.savefig(plot_out_name + ".png", dpi=300, bbox_inches='tight')
plt.savefig(plot_out_name + ".pdf", bbox_inches='tight')
