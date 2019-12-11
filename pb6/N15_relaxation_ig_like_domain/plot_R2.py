import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plot_out_name = 'R2_full.png'

pb6_seq = 'VTAITVKSAGNVTTLNRSATLQMSVEVTPSSARNKEVTWAITAGDAATINATGLLRADASKTGAVTVEATAKDGSGVKGTKVITVTAGG'

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
# print(resids, resnames, R2)

plt.rcParams.update({'font.size': 14})
plt.figure(figsize=cm2inch(16.0, 12.0))

pb6_resids = resids
pb6_R2 = R2

# Experimetnal data
# pb6_exp_resids = np.array([3, 7, 9, 10, 11, 14, 15])
# pb6_exp_R2 = np.array([1.26765865, 1.291394278, 1.166972683, 1.213193526, 1.17605744, 1.228871478, 1.201527367])

# plot
fig, ax = plt.subplots()
ax.plot(pb6_resids, pb6_R2, marker="D", ms=5, markeredgecolor="b", markerfacecolor="b",
        linewidth=2.0, color="b", label="pb6")
# ax.plot(pb6_exp_resids, pb6_exp_R2, marker="o", ms=7, markeredgecolor="r", markerfacecolor="r",
#         linestyle='None', label="Experiment")

plt.ylabel(r'${\rm R_{2},\ s^{-1}}$')
plt.xlabel('Residue')
plt.axis([0, 90, 0, 20])
plt.xticks([1] + list(range(10, 93, 10)))
plt.yticks(np.arange(0, 20.01, 5.0))

# Put a legend above current axis
ax.legend(loc='lower left', bbox_to_anchor=(0.01, 1.01), ncol=1,
          borderaxespad=0, frameon=False, numpoints=1)

plt.savefig(plot_out_name, dpi=300, bbox_inches='tight')
