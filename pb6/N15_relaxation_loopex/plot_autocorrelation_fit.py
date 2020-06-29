import matplotlib as mpl
mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from correlation_functions import Corfun
import json

seq = "  "
with open("seq.fasta", "r") as f:
    for line in f:
        if line[0] != ">":
            seq += line.rstrip("\n")

#
#  setup trajectory parameters
#
with open("input.json", "r") as f:
    pars = json.load(f)
first_dat_file = int(pars["first_dat_file"])
last_dat_file = int(pars["last_dat_file"])
n_exp = pars["n_exp"]
step = pars["step"]
fit_out = pars["fit_out"]
plot_out_name = pars["autocor_fig_name"]
axis_limits = pars["autocor_axis_limits"]

# read data
cor_align = "cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file)
cor = "cor_NH_%d-%d_tumbling_on" % (first_dat_file, last_dat_file)
files_align = sorted(glob(cor_align+'/*.cor'))
files = sorted(glob(cor+'/*.cor'))

assert len(files) == len(files_align), (len(files), len(files_align))

# read fitted amps and taus
with open('autocor_fit_results.txt', 'r') as af:
    acfs = af.readlines()

# create file for plot
plt.rcParams.update({'font.size': 10})
pp = PdfPages(plot_out_name)
plt.figure(figsize=(20.0/2.54, 15.0/2.54))

for file, file_align, a in zip(files, files_align, acfs):
    print(file)

    # get rId and rName
    rId = int(a.split()[0].split('-')[0])
    rName = a.split()[0].split('-')[1]

    # get amps and taus
    amps = [float(x) for x in a.split()[1:1 + n_exp]]
    taus = [float(x) / step for x in a.split()[1 + n_exp:1 + 2 * n_exp]]

    # sort amps and taus by tau
    amps_sorted = [a for t, a in sorted(zip(taus, amps))]
    taus_sorted = [t for t, a in sorted(zip(taus, amps))]
    amps = amps_sorted
    taus = taus_sorted

    # get data
    data = np.genfromtxt(file)
    data_align = np.genfromtxt(file_align)
    parms = amps + taus
    data_fit = Corfun.multiexponent(data[:, 0], *parms)

    # plot
    plt.clf()
    plt.plot(data_align[:, 0] / 1000, data_align[:, 1], 'g-', linewidth=2)
    plt.plot(data[:, 0] / 1000, data[:, 1], 'b-', linewidth=3)
    plt.plot(data[:, 0] / 1000, data_fit, 'r-', linewidth=2)
    plt.xlabel('Time, ns')
    plt.ylabel('Correlation function')
    plt.title(f"{seq[rId - 1]}{rId:d}", fontsize=18)

    txt = "amp / tau[ns]\n"
    i = 0
    for a, t in zip(amps, taus):
        i += 1
        if a < 0.01:
            continue
        if i == len(amps):
            txt += "{:.2f} / {:.4f}".format(np.around(a, decimals=2), np.around(t / 1e3, decimals=4))
        else:
            txt += "{:.2f} / {:.4f}\n".format(np.around(a, decimals=2), np.around(t / 1e3, decimals=4))
    t = plt.text(x=0.65,
                 y=0.75,
                 s=txt,
                 horizontalalignment='left',
                 verticalalignment='top',
                 bbox=dict(facecolor='white'),
                 fontsize=10,
                 transform=plt.gca().transAxes)
    t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))

    plt.axis(axis_limits)
    plt.legend(['%s%d tumbling off' % (rName, rId), '%s%d tumbling on' % (rName, rId), '%s%d fit' % (rName, rId)],
               loc="upper center")
    plt.grid(True)
    pp.savefig()

plt.close()
pp.close()
