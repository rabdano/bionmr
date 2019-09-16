import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 14})
H4_seq = 'SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG'

h4_1_adherence = np.loadtxt('adherence_avg_H4-1_KDtree.txt')
h4_2_adherence = np.loadtxt('adherence_avg_H4-2_KDtree.txt')

plt.figure(figsize=(12, 9))
fig, ax = plt.subplots()

plt.plot(range(len(h4_1_adherence)), h4_1_adherence, marker='D', ms=7, markeredgecolor='b', markerfacecolor='b',
         linewidth=2.0, color='b', label='H4-1 [136-159]')
plt.plot(range(len(h4_2_adherence)), h4_2_adherence, marker='D', ms=7, markeredgecolor='g', markerfacecolor='g',
        linewidth=2.0, color='g', label='H4-2 [623-646]')

plt.ylabel('dist(Ca, nucleosome core), $\AA$')

plt.xticks(range(len(h4_1_adherence)))

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
plt.savefig('adherence.png', dpi=300, bbox_inches='tight')
