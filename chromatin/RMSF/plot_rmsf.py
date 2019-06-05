import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import os
from matplotlib.backends.backend_pdf import PdfPages


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


files = sorted(glob('RMSF_chain_*.csv'))
chains = [f[-5] for f in files]

for fn, c in zip(files, chains):
    plt.figure(figsize=cm2inch(16.0, 12.0))
    fig, ax = plt.subplots()

    data = np.genfromtxt(fn, delimiter=',', unpack=True, skip_header=1, dtype=None)
    resids = [d[0] for d in data]
    resnames = [d[1].decode('UTF-8') for d in data]
    rmsf = [d[2] for d in data]

    plt.step(range(len(rmsf)), rmsf, where="mid")

    plt.ylabel("RMSF, $\AA$")
    plt.grid(color="#CCCCCC", lw=0.1)


    def to_label(a):
        from Bio.PDB.Polypeptide import three_to_one
        if a == 'HID':
            a = 'HIS'
        return "%s" % (three_to_one(a))

    plt.xticks(range(len(rmsf)),
               [to_label(a) for a in resnames],
               rotation=0, fontsize="x-small")

    plt.ylim((0,25))

    plt.savefig('RMSF_chain_' + c + '.png')




#
#
# nve = np.genfromtxt('NVE/cor_NH_n300_s1/' + f)
# npt_g_2 = np.genfromtxt('NPT_gamma_ln_2/cor_NH_n245_s1/' + f)
# npt_g_001 = np.genfromtxt('NPT_gamma_ln_0.01/cor_NH_n245_s1/' + f)
#
# ax.plot(nve[:,0], nve[:,1], linewidth=0.5, color="k", label="NVE")
# ax.plot(npt_g_2[:,0], npt_g_2[:,1], linewidth=0.5, color="b", label="NPT_gamma_ln_2")
# ax.plot(npt_g_001[:,0], npt_g_001[:,1], linewidth=0.5, color="r", label="NPT_gamma_ln_0.01")
#
# plt.xlim([-100, 2000])
# plt.ylim([-0.1, 1.1])
# plt.xlabel('Time, ps')
# plt.ylabel('ACF NH')
#
# ax.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=3,
#         borderaxespad=0, frameon=False, numpoints=1)
# plt.grid(True)
#
# #         plt.savefig('figures/' + f + '.png')
# pdf.savefig()
# plt.close('all')