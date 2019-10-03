import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import os
from matplotlib.backends.backend_pdf import PdfPages


plt.rcParams.update({'font.size': 14})
H4_seq = 'SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG'


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
    # plt.grid(color="#CCCCCC", lw=0.1)


    def to_label(a):
        from Bio.PDB.Polypeptide import three_to_one
        if (a == 'HID') | (a == 'HIP') | (a == 'HIE'):
            a = 'HIS'
        return "%s" % (three_to_one(a))

    plt.xticks(range(len(rmsf)),
               [to_label(a) for a in resnames],
               rotation=0, fontsize="x-small")

    plt.ylim((0, 25))

    plt.savefig('RMSF_chain_' + c + '.png')




plt.figure(figsize=cm2inch(16.0, 12.0))
fig, ax = plt.subplots()

data1 = np.genfromtxt('RMSF_chain_B.csv', delimiter=',', unpack=True, skip_header=1, dtype=None)
data2 = np.genfromtxt('RMSF_chain_F.csv', delimiter=',', unpack=True, skip_header=1, dtype=None)
resids = [d[0] for d in data]
resnames = [d[1].decode('UTF-8') for d in data]
rmsf1 = [d[2] for d in data1]
rmsf2 = [d[2] for d in data2]

plt.plot(range(len(rmsf1)), rmsf1, marker='D', ms=7, markeredgecolor='b', markerfacecolor='b',
         linewidth=2.0, color='b', label='H4-1')
plt.plot(range(len(rmsf2)), rmsf2, marker='D', ms=7, markeredgecolor='g', markerfacecolor='g',
        linewidth=2.0, color='g', label='H4-2')

plt.ylabel("RMSF, $\AA$")
# plt.grid(color="#CCCCCC", lw=0.1)


def to_label(a):
    from Bio.PDB.Polypeptide import three_to_one
    if (a == 'HID') | (a == 'HIP') | (a == 'HIE'):
        a = 'HIS'
    return "%s" % (three_to_one(a))


plt.xticks(range(len(rmsf)))
plt.yticks(range(0, 21, 2))

n_labels = len(ax.get_xticklabels())
labels = []
for i in range(0, n_labels):
    if (((i+1) % 5) == 0) | (i == 0):
        dig = str(i+1)
    else:
        dig = ''
    labels.append('{}\n{}'.format(dig, H4_seq[i]))
ax.set_xticklabels(labels)

plt.ylim((0, 20))

plt.legend()

plt.savefig('RMSF_chain_BF.png', dpi=300, bbox_inches='tight')
