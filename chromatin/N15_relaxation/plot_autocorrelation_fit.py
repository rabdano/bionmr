import matplotlib as mpl
mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from correlation_functions import Corfun

n_exp = ___N_EXP___
plot_out_name = 'autocor_fit_results.pdf'
step = 1e-12

# read data
cor = 'cor_NH____FIRST_DAT_FILE___-___LAST_DAT_FILE____tumbling_on'
files = sorted(glob(cor+'/*.cor'))

# read fitted amps and taus
with open('autocor_fit_results.txt', 'r') as af:
    acfs = af.readlines()

# create file for plot
pp = PdfPages(plot_out_name)
plt.figure(figsize=(16.0/2.54, 12.0/2.54))

for file, a in zip(files, acfs):
    print(file)

    # get rId and rName
    rId = int(a.split()[0].split('-')[0])
    rName = a.split()[0].split('-')[1]
    if rId < 300:
        rId = rId - 135
    else:
        rId = rId - 622

    # get amps and taus
    amps = [float(x) for x in a.split()[1:1+n_exp]]
    taus = [float(x)/step for x in a.split()[1+n_exp:1+2*n_exp]]

    # get data
    data = np.genfromtxt(file)
    parms = amps + taus
    data_fit = Corfun.multiexponent(data[:, 0], *parms)

    # plot
    plt.clf()
    plt.plot(data[:, 0], data[:, 1], 'b-', linewidth=0.5)
    plt.plot(data[:, 0], data_fit, 'r-', linewidth=0.5)
    plt.xlabel('Time, ps')
    # axis_limits = [-10, data[-1, 0]+10, -0.2, 1.01]
    axis_limits = [-1000, 101000, 0, 1.01]
    txt = "amp / tau[ns]\n"
    i = 0
    for a, t in zip(amps, taus):
        i += 1
        if i == len(amps):
            txt += "{:.2f} / {:.3f}".format(a, t/1e3)
        else:
            txt += "{:.2f} / {:.3f}\n".format(a, t/1e3)
    plt.text(x=0.7,
             y=0.7,
             s=txt,
             horizontalalignment='left',
             verticalalignment='top',
             bbox=dict(facecolor='white'),
             fontsize=12,
             transform=plt.gca().transAxes)
    plt.axis(axis_limits)
    plt.legend(['%s%d data' % (rName, rId), '%s%d fit' % (rName, rId)])
    plt.grid(True)
    pp.savefig()

plt.close()
pp.close()
