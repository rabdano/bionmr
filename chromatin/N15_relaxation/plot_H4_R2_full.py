import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np

plot_out_name = 'R2.png'

H4_seq = 'SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG'

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

# read data
with open("R2_results.txt", "r") as f:
    lines = f.readlines()
n_resids = len(lines)
resids = np.zeros(n_resids, dtype=int)
resnames = []
R2 = np.zeros(n_resids)
for i, line in enumerate(lines):
    resids[i] = int(line.split()[0].split('-')[0])
    resnames.append(line.split()[0].split('-')[1])
    R2[i] = float(line[:-1].split()[1])
# print(resids,resnames,R2)

plt.rcParams.update({'font.size': 14})
plt.figure(figsize=cm2inch(16.0, 12.0))

H4_resids = resids[:100] - 135
H4_R2_1 = R2[:100]
H4_R2_2 = R2[100:]
H4_R2 = 0.5 * (H4_R2_1 + H4_R2_2)
H4_R2_err = np.abs(H4_R2 - H4_R2_1)

# Experimetnal data
H4_exp_resids = np.array([3, 7, 9, 10, 11, 14, 15])
H4_exp_R2 = np.array([16.39998754, 12.36114589, 11.62683573, 12.52521179, 10.97878061, 9.914836342, 20.27790368])

# plot
fig, ax = plt.subplots()
ax.plot(H4_resids, H4_R2_1, marker="D", ms=7, markeredgecolor="b", markerfacecolor="b",
        linewidth=2.0, color="b", label="H4-1")
# ax.plot([], [], marker="D", ms=7, markeredgecolor="b", markerfacecolor="b",
#         linewidth=2.0, color="b", label="H4-1 [136-159]")
ax.plot(H4_resids, H4_R2_2, marker="D", ms=7, markeredgecolor="g", markerfacecolor="g",
        linewidth=2.0, color="g", label="H4-2")
# ax.plot([], [], marker="D", ms=7, markeredgecolor="g", markerfacecolor="g",
#         linewidth=2.0, color="g", label="H4-2 [623-646]")
ax.plot(H4_exp_resids, H4_exp_R2, marker="o", ms=10, markeredgecolor="r", markerfacecolor="r",
        linestyle='None', label="Experiment")

plt.ylabel(r'${\rm R_{2},\ s^{-1}}$')
plt.xlabel('Residue')
plt.axis([0, 103, 0, 350])
plt.xticks(np.arange(1, 102, 1))
plt.yticks(np.arange(0, 351, 50))
n_labels = len(ax.get_xticklabels())
labels = []
for i in range(0, n_labels):
    if (((i+1) % 10) == 0) | (i == 0):
        dig = str(i+1)
    else:
        dig = ''
    # labels.append('{}\n{}'.format(dig, H4_seq[i]))
    labels.append('{}'.format(dig))
ax.set_xticklabels(labels)
# Put a legend above current axis
ax.legend(loc='lower left', bbox_to_anchor=(0.01, 1.01), ncol=3,
          borderaxespad=0, frameon=False, numpoints=1)






axins = inset_axes(ax, width="80%", height="80%",
                   bbox_to_anchor=(.3, .1, .6, .5),
                   bbox_transform=ax.transAxes, loc=1)

axins.plot(H4_resids, H4_R2_1, marker="D", ms=7, markeredgecolor="b", markerfacecolor="b",
        linewidth=2.0, color="b", label="H4-1")
axins.plot(H4_resids, H4_R2_2, marker="D", ms=7, markeredgecolor="g", markerfacecolor="g",
        linewidth=2.0, color="g", label="H4-2")
axins.plot(H4_exp_resids, H4_exp_R2, marker="o", ms=10, markeredgecolor="r", markerfacecolor="r",
        linestyle='None', label="Experiment")

axins.axis([0, 16, 0, 40])

axins.set_xticks(np.arange(1, 16, 1))
axins.set_yticks(np.arange(0, 41, 10))
n_labels = len(axins.get_xticklabels())
labels = []
for i in range(0, n_labels):
    if (((i+1) % 5) == 0) | (i == 0):
        dig = str(i+1)
    else:
        dig = ''
    labels.append('{}\n{}'.format(dig, H4_seq[i]))
    # labels.append('{}'.format(dig))
axins.set_xticklabels(labels)





plt.savefig(plot_out_name, dpi=300, bbox_inches='tight')

