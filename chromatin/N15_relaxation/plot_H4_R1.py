import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plot_out_name = 'R1.png'

H4_seq = 'SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG'

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

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
# print(resids,resnames,R2)

plt.rcParams.update({'font.size': 14})
plt.figure(figsize=cm2inch(16.0, 12.0))

H4_resids = resids[:100] - 135
H4_R1_1 = R1[:100]
H4_R1_2 = R1[100:]
H4_R1 = 0.5 * (H4_R1_1 + H4_R1_2)
H4_R1_err = np.abs(H4_R1 - H4_R1_1)

# Experimetnal data
H4_exp_resids = np.array([3, 7, 9, 10, 11, 14, 15])
H4_exp_R1 = np.array([1.26765865, 1.291394278, 1.166972683, 1.213193526, 1.17605744, 1.228871478, 1.201527367])

# plot
fig, ax = plt.subplots()
ax.plot(H4_resids, H4_R1_1, marker="D", ms=7, markeredgecolor="b", markerfacecolor="b",
        linewidth=2.0, color="b", label="H4-1")
# ax.plot([], [], marker="D", ms=7, markeredgecolor="b", markerfacecolor="b",
#         linewidth=2.0, color="b", label="H4-1 [136-159]")
ax.plot(H4_resids, H4_R1_2, marker="D", ms=7, markeredgecolor="g", markerfacecolor="g",
        linewidth=2.0, color="g", label="H4-2")
# ax.plot([], [], marker="D", ms=7, markeredgecolor="g", markerfacecolor="g",
#         linewidth=2.0, color="g", label="H4-2 [623-646]")
ax.plot(H4_exp_resids, H4_exp_R1, marker="o", ms=10, markeredgecolor="r", markerfacecolor="r",
        linestyle='None', label="Experiment")

plt.ylabel(r'${\rm R_{1},\ s^{-1}}$')
plt.axis([0, 26, 0, 2])
plt.xticks(np.arange(1, 25, 1))
plt.yticks(np.arange(0, 2.1, 0.2))
n_labels = len(ax.get_xticklabels())
labels = []
for i in range(0, n_labels):
    if (((i+1) % 5) == 0) | (i == 0):
        dig = str(i+1)
    else:
        dig = ''
    labels.append('{}\n{}'.format(dig, H4_seq[i]))
ax.set_xticklabels(labels)
# Put a legend above current axis
ax.legend(loc='lower left', bbox_to_anchor=(0.01, 1.01), ncol=3,
          borderaxespad=0, frameon=False, numpoints=1)
plt.savefig(plot_out_name, dpi=300, bbox_inches='tight')

