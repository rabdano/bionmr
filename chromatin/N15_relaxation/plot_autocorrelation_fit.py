import matplotlib as mpl
mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from correlation_functions import Corfun

n_exp = 3
plot_out_name = 'autocor_fit_results.pdf'


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


# read data
cor = 'cor_NH_51-300_diluted'
files = sorted(glob(cor+'/*.cor'))

# read fitted amps and taus
with open('autocor_fit_results.txt','r') as af:
    acfs = af.readlines()

# create file for plot
pp = PdfPages(plot_out_name)
plt.figure(figsize=cm2inch(16.0, 12.0))

for file, a in zip(files, acfs):
    print(file)

    # get rId and rName
    rId = int(a.split()[0].split('-')[0])
    rName = a.split()[0].split('-')[1]

    # get amps and taus
    amps = [float(x) for x in a.split()[1:1+n_exp]]
    taus = [float(x) for x in a.split()[1+n_exp:1+2*n_exp]]

    # get data
    data = np.genfromtxt(file)
    parms = amps + taus
    data_fit = Corfun.multiexponent(data[:, 0], *parms)

    # plot
    plt.clf()
    plt.plot(data[:, 0], data[:, 1], 'b-', linewidth=0.5)
    plt.plot(data[:, 0], data_fit, 'r-', linewidth=0.5)
    plt.xlabel('Time, ps')
    axis_limits = [-10, data[-1, 0]+10, -0.2, 1.01]
    plt.axis(axis_limits)
    plt.legend(['%s%d data' % (rName, rId), '%s%d fit' % (rName, rId)])
    plt.grid(True)
    pp.savefig()

plt.close()
pp.close()
