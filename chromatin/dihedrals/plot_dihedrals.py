import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import sys

stride = 100


H4_seq = 'SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG'


def correct_numbering(residues_of_interest, renumber_map=None):
    if not renumber_map:
        renumber_map = [rid for rid in range(1, 136)] + \
                       [rid for rid in range(1, 103)] + \
                       [rid for rid in range(1, 129)] + \
                       [rid for rid in range(1, 123)] + \
                       [rid for rid in range(1, 136)] + \
                       [rid for rid in range(1, 103)] + \
                       [rid for rid in range(1, 129)] + \
                       [rid for rid in range(1, 123)] + \
                       [rid for rid in range(1, 148)] + \
                       [rid for rid in range(1, 148)]

    residues_of_interest = [renumber_map[resid - 1] for resid in residues_of_interest]

    return residues_of_interest


phi_df = pd.read_csv('phi.csv')
psi_df = pd.read_csv('psi.csv')

header = phi_df.columns.values.tolist()
residues_of_interest = [int(x) for x in header]
residues_of_interest = correct_numbering(residues_of_interest)
n = len(phi_df.index)

with PdfPages('phi_figures.pdf') as pdf:
    for i, resid in enumerate(header):
        sys.stdout.write('Plotting phi figure %d of %d...\r' % (i, len(header)))
        plt.figure(figsize=(16, 9))
        x = np.linspace(0, n, n) / 1000 * stride
        y = phi_df[resid]
        plt.scatter(x, y, marker='o', s=4, alpha=0.7, lw=0, color="b",
                    label='{}{}'.format(H4_seq[residues_of_interest[i] - 1], residues_of_interest[i]))
        plt.xlabel('Time, ns')
        plt.ylabel(r'$\varphi$')
        plt.ylim([-185.0, 185.0])
        plt.yticks(range(-180, 181, 20))
        plt.legend()
        pdf.savefig(bbox_inches='tight')
        plt.close()
sys.stdout.write('phi - done!' + ' ' * 20 + '\n')

with PdfPages('psi_figures.pdf') as pdf:
    for i, resid in enumerate(header):
        sys.stdout.write('Plotting psi figure %d of %d...\r' % (i, len(header)))
        plt.figure(figsize=(16, 9))
        x = np.linspace(0, n, n) / 1000 * stride
        y = psi_df[resid]
        plt.scatter(x, y, marker='o', s=4, alpha=0.7, lw=0, color="b",
                    label='{}{}'.format(H4_seq[residues_of_interest[i] - 1], residues_of_interest[i]))
        plt.xlabel('Time, ns')
        plt.ylabel(r'$\psi$')
        plt.ylim([-185.0, 185.0])
        plt.yticks(range(-180, 181, 20))
        plt.legend()
        pdf.savefig(bbox_inches='tight')
        plt.close()
sys.stdout.write('psi - done!' + ' ' * 20 + '\n')
