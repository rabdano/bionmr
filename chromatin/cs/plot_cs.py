import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plot_out_name = 'H4'
RC_path = '/home/seva/scripts/git/bionmr/chromatin/cs' + '/H4_Tamiola.cs'
H4_seq = 'SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG'

plt.rcParams.update({'font.size': 14})

# read data
with open("report.txt", "r") as f:
    lines = f.readlines()
n_resids = 102
resids = np.zeros(n_resids * 2, dtype=int)
resnames = []
labels = ["H", "N", "CA", "CB"]
cs = np.zeros((n_resids * 2, 4))  # H, N, CA, CB (C and HA shifts are ignored)
first = 0
second = 0
for i, line in enumerate(lines):
    if "First tail" in line:
        first = i
    if "Second tail" in line:
        second = i
data = lines[first + 2:first + 2 + n_resids] + lines[second + 2:second + 2 + n_resids]

# get data from final_report.txt
for i, line in enumerate(data):
    resids[i] = int(line.split()[0])
    resnames.append(line.split()[1])
    cs[i, 0] = float(line.split()[2])
    cs[i, 1] = float(line.split()[3])
    cs[i, 2] = float(line.split()[4])
    cs[i, 3] = float(line.split()[5])

# Experimetnal data
H4_exp_resids = np.array([3, 7, 9, 10, 11, 13, 14, 15])
H4_exp_cs = np.array([[8.500, 120.775, 56.725, 31.226],
                      [8.338, 108.846, 45.665, np.nan],
                      [8.545, 109.823, 45.686, np.nan],
                      [8.206, 121.442, 55.551, 42.715],
                      [8.517, 109.533, 45.662, np.nan],
                      [8.592, 110.159, 45.770, np.nan],
                      [8.302, 108.902, 45.451, np.nan],
                      [8.230, 123.793, 52.755, 19.764]])


# Random coil shifts
def get_RC(rId, aName, path_to_RC=RC_path):
    if aName == "H":
        aName = "HN"
    f = open(path_to_RC, "r")
    lines = [line.rstrip("\n") for line in f]
    f.close()
    # print(rId, aName, lines)
    RC = 0.0
    for line in lines:
        s = line.split()
        # print(s,rId, int(s[0].strip()), aName, str(s[2].strip()))
        if (rId == int(s[0].strip())) & (aName == str(s[3].strip())):
            RC = float(s[4].strip())
            break
    return RC


rc_cs = np.zeros((n_resids, 4))  # H, N, CA, CB (C and HA shifts are ignored)
for i in range(n_resids):
    for j, l in enumerate(labels):
        rc_cs[i, j] = get_RC(i + 1, l)

# Pick data to H4-tail only
H4_resids = list(range(1, 25))
d_mins = [6.5, 100, 40, 10]
d_maxs = [10, 135, 70, 70]
d_step = [0.5, 5.0, 5.0, 5.0]
d_mins_diff = [-1, -5, -3, -5]
d_maxs_diff = [1, 5, 3, 5]
d_step_diff = [0.2, 1.0, 1.0, 1.0]

# Plot

for i, l in enumerate(labels):

    # Calculate RMSD between SHIFTS and Experiment

    sd_1 = np.zeros(H4_exp_resids.shape)
    sd_2 = np.zeros(H4_exp_resids.shape)

    for idx, rId in enumerate(H4_exp_resids):
        cs_exp = H4_exp_cs[idx, i]
        cs_1 = cs[rId - 1, i]
        cs_2 = cs[rId + 102 - 1, i]
        sd_1[idx] = np.sqrt((cs_exp - cs_1) ** 2)
        sd_2[idx] = np.sqrt((cs_exp - cs_2) ** 2)

    rmsd_1 = np.nanmean(sd_1)
    rmsd_2 = np.nanmean(sd_2)

    # Plot absolute values
    plt.figure(figsize=(8, 6))

    fig, ax = plt.subplots()
    ax.plot(H4_resids, cs[:24, i], marker="D", ms=5, markeredgecolor="b", markerfacecolor="b",
            linewidth=1.0, color="b", label="H4-1 136-159")
    ax.plot(H4_resids, cs[102:102 + 24, i], marker="D", ms=5, markeredgecolor="g", markerfacecolor="g",
            linewidth=1.0, color="g", label="H4-2 623-646")
    ax.plot(H4_exp_resids, H4_exp_cs[:, i], marker="o", ms=5, markeredgecolor="r", markerfacecolor="r",
            linestyle='None', label="Experimental")
    ax.plot(H4_resids, rc_cs[:24, i], marker="_", ms=5, markeredgecolor="k", markerfacecolor="k",
            linestyle='None', label="RC")

    # plt.xlabel(r'${\rm Residue}$')

    plt.xticks(np.arange(1, 25, 1))
    n_labels = len(ax.get_xticklabels())
    labels = []
    for j in range(0, n_labels):
        if (((j + 1) % 5) == 0) | (j == 0):
            dig = str(j + 1)
        else:
            dig = ''
        labels.append('{}\n{}'.format(dig, H4_seq[j]))
    ax.set_xticklabels(labels)

    plt.ylabel('\u03B4' + l + " , ppm")
    plt.axis([0, 26, d_mins[i], d_maxs[i]])
    plt.yticks(np.arange(d_mins[i], d_maxs[i] + 0.1, d_step[i]))
    # Put a legend above current axis
    ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=2,
              borderaxespad=0, frameon=False, numpoints=1)
    # plt.grid(True)

    plt.text(0.98, 0.98, "RMSD: %.2f, %.2f" % (rmsd_1, rmsd_2),
             horizontalalignment='right',
             verticalalignment='top',
             transform=ax.transAxes)

    plt.savefig(plot_out_name + "_" + l + ".png", dpi=300, bbox_inches='tight')

    plt.close()

    # Plot difference with random coil
    plt.figure(figsize=(8, 6))

    fig, ax = plt.subplots()
    ax.plot(H4_resids, cs[:24, i] - rc_cs[:24, i], marker="D", ms=7, markeredgecolor="b", markerfacecolor="b",
            linewidth=2.0, color="b", label="H4-1 [136-159]")
    ax.plot(H4_resids, cs[102:102 + 24, i] - rc_cs[:24, i], marker="D", ms=7, markeredgecolor="g", markerfacecolor="g",
            linewidth=2.0, color="g", label="H4-2 [623-646]")
    ax.plot(H4_exp_resids, H4_exp_cs[:, i] - rc_cs[H4_exp_resids - 1, i], marker="o", ms=10, markeredgecolor="r",
            markerfacecolor="r",
            linestyle='None', label="Experimental")

    # plt.xlabel('Residue')

    plt.xticks(np.arange(1, 25, 1))
    n_labels = len(ax.get_xticklabels())
    labels = []
    for j in range(0, n_labels):
        if (((j + 1) % 5) == 0) | (j == 0):
            dig = str(j + 1)
        else:
            dig = ''
        labels.append('{}\n{}'.format(dig, H4_seq[j]))
    ax.set_xticklabels(labels)

    plt.ylabel('\u03B4' + l + ' - \u03B4' + l + r'$_{\rm RC}$, ppm')
    plt.axis([0, 26, d_mins_diff[i], d_maxs_diff[i]])
    plt.yticks(np.arange(d_mins_diff[i], d_maxs_diff[i] + 0.1, d_step_diff[i]))
    # Put a legend above current axis
    ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=2,
              borderaxespad=0, frameon=False, numpoints=1)
    # plt.grid(True)

    plt.text(0.98, 0.98, "RMSD: %.2f, %.2f" % (rmsd_1, rmsd_2),
             horizontalalignment='right',
             verticalalignment='top',
             transform=ax.transAxes)

    plt.savefig(plot_out_name + "_" + l + "_diff.png", dpi=300, bbox_inches='tight')

    plt.close()
