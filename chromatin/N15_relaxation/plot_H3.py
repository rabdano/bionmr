import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plot_out_name = 'H3_ff14SB_TIP3P_0M_NaCl_NPT.png'
label_H3 = "H3, ff14SB, TIP3P, 0 M NaCl, NPT"

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

plt.figure(figsize=cm2inch(16.0, 12.0))

H3_resids = resids[0:39]
H3_R2_1 = R2[0:39]
H3_R2_2 = R2[62:101]
H3_R2 = 0.5 * (H3_R2_1 + H3_R2_2)
H3_R2_err = np.abs(H3_R2 - H3_R2_1)

# plot
fig, ax = plt.subplots()
ax.errorbar(H3_resids, H3_R2, yerr=H3_R2_err, marker="D", ms=5, markeredgecolor="r", markerfacecolor="r",
    ecolor="k", color="k", linewidth=1.0, label=label_H3)
plt.xlabel(r'${\rm Residue}$')
plt.ylabel(r'${\rm R_{2},\ s^{-1}}$')
plt.axis([0,45,0,165])
plt.yticks(np.arange(0,221,20))
# Put a legend above current axis
ax.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2, 
            borderaxespad=0, frameon=False, numpoints=1)
plt.grid(True)
plt.savefig(plot_out_name, dpi=300, bbox_inches='tight')

# plt.close()
# pp.close()