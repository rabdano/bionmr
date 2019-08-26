import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

H4_seq = 'SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG'

with open('extract_contact_map_Ca.py', 'r') as f:
    for line in f:
        if 'fout_1' in line:
            h4_map_1_fn = line.split("'")[-2]
        if 'fout_1' in line:
            h4_map_2_fn = line.split("'")[-2]
            break

h4_map_1 = np.flip(np.loadtxt(h4_map_1_fn), axis=0)
h4_map_2 = np.flip(np.loadtxt(h4_map_2_fn), axis=0)
n_residues = len(h4_map_1)

fig, ax = plt.subplots(figsize=(8, 8))
cmap = plt.get_cmap('jet', np.ceil(np.max(np.max(h4_map_1))/5))

im = ax.imshow(h4_map_1, cmap=cmap, aspect='equal', vmin=0.0, vmax=5.0*np.ceil(np.max(np.max(h4_map_1))/5))

plt.xticks(range(n_residues))
plt.yticks(range(n_residues))
n_labels = len(ax.get_xticklabels())
x_labels = []
y_labels = []
for i in range(0, n_labels):
    if (((i+1) % 5) == 0) | (i == 0):
        dig = str(i+1) + ' '
    else:
        dig = ''
    x_labels.append('{}\n{}'.format(dig, H4_seq[i]))
    y_labels.append('{}{}'.format(dig, H4_seq[i]))
y_labels.reverse()

ax.set_xticklabels(x_labels)
ax.set_yticklabels(y_labels)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
cbar.ax.set_ylabel('Distance / $\AA$')

plt.savefig('H4-1_contact_map.png', bbox_inches='tight')
plt.clf()
plt.close(fig)


fig, ax = plt.subplots(figsize=(8, 8))
cmap = plt.get_cmap('jet', np.ceil(np.max(np.max(h4_map_2))/5))

im = ax.imshow(h4_map_2, cmap=cmap, aspect='equal', vmin=0.0, vmax=5.0*np.ceil(np.max(np.max(h4_map_1))/5))

plt.xticks(range(n_residues))
plt.yticks(range(n_residues))
n_labels = len(ax.get_xticklabels())
x_labels = []
y_labels = []
for i in range(0, n_labels):
    if (((i+1) % 5) == 0) | (i == 0):
        dig = str(i+1) + ' '
    else:
        dig = ''
    x_labels.append('{}\n{}'.format(dig, H4_seq[i]))
    y_labels.append('{}{}'.format(dig, H4_seq[i]))
y_labels.reverse()

ax.set_xticklabels(x_labels)
ax.set_yticklabels(y_labels)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
cbar.ax.set_ylabel('Distance / $\AA$')

plt.savefig('H4-2_contact_map.png', bbox_inches='tight')
plt.clf()
plt.close(fig)
