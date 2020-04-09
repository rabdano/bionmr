import matplotlib as mpl
from bionmr_utils.md import *
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import json
from traces_func_lib import *
from timeit import default_timer as timer

# User input:
sb_trace = 'sb_trace.dat'  # file with data about hydrogen bonds
hb_input_json = 'hb_input.json'  # file with description of hydrogen bond extraction parameters
hb_trace = 'hb_trace.dat'  # file with data about salt bridges
hb_avg_win = 50  # window for averaging of data
sb_avg_win = 50
hb_stride = 100
sb_stride = 100
size_of_figure = (21.0/2/2.54, 29.7/2.54*3)

h4_1 = list(range(136, 160))
h4_2 = list(range(623, 647))
rid_of_interest = h4_1 + h4_2
H4_seq = 'SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG'
ref = 'ref.pdb'
ref_frame = PdbFile(ref).get_frame()


t2 = timer()


# read HB parameters
with open(hb_input_json, 'r') as f:
    pars = json.load(f)
hb_trj_filename_first = int(pars['trj_filename_first'])
hb_trj_filename_last = int(pars['trj_filename_last'])


# extract HB header
with open(hb_trace, 'r') as f:
    header = f.readline()
hbs = hb_split_header(header)


# read HB data
# hb_data = np.genfromtxt(hb_trace, skip_header=1)  # slow version
n_frames = int(1000 * (hb_trj_filename_last - hb_trj_filename_first + 1) / hb_stride)
hb_data = np.zeros((n_frames, len(hbs)), dtype=int)
with open(hb_trace, 'r') as f:
    for i, line in enumerate(f):
        if i > 0:
            hb_data[i-1, :] = np.fromstring(line.rstrip('\n'), dtype=int, sep=' ')


# average HB data
hb_avg_data = average_data(hb_data, hb_avg_win)


# extract SB header
with open(sb_trace, 'r') as f:
    header = f.readline()
sbs = sb_split_header(header, ref_frame)


# read SB data
# sb_data = np.genfromtxt(sb_trace, skip_header=1)  # slow version
n_frames = int(1000 * (hb_trj_filename_last - hb_trj_filename_first + 1) / sb_stride)
sb_data = np.zeros((n_frames, len(sbs)), dtype=int)
with open(sb_trace, 'r') as f:
    for i, line in enumerate(f):
        if i > 0:
            sb_data[i-1, :] = np.fromstring(line.rstrip('\n'), dtype=int, sep=' ')


# average SB data
sb_avg_data = average_data(sb_data, sb_avg_win)


t1, t2 = t2, timer()
print("Read data:", t2-t1, "s")


# prune data for not H4 tails
hbs, hb_avg_data = prune_not_h4_contacts(hbs, hb_avg_data, rid_of_interest)
sbs, sb_avg_data = prune_not_h4_contacts(sbs, sb_avg_data, rid_of_interest)


# correct numbering
hbs = correct_numbering(hbs)
sbs = correct_numbering(sbs)


# prune intra-tail interactions
hbs, hb_avg_data = prune_intra_tail_contacts(hbs, hb_avg_data, len(h4_1))
sbs, sb_avg_data = prune_intra_tail_contacts(sbs, sb_avg_data, len(h4_1))


# prune traces with no color on map
hbs, hb_avg_data = prune_empty_traces(hbs, hb_avg_data)
sbs, sb_avg_data = prune_empty_traces(sbs, sb_avg_data)


# join SB data for same residues
sbs, sb_avg_data = joint_sb_data_for_same_residue(sbs, sb_avg_data)


# remove traces for Na+ ions
hbs, hb_avg_data = prune_sodium_contacts(hbs, hb_avg_data)
sbs, sb_avg_data = prune_sodium_contacts(sbs, sb_avg_data)


t1, t2 = t2, timer()
print("Prune redundant traces:", t2-t1, "s")


#
#  H4-1
#


# build cmap from HB and SB
cmap, labels, numbers_of_lines = build_cmap(hbs, sbs, hb_avg_data, sb_avg_data, "B", len(h4_1))
rows = np.sum(numbers_of_lines)


t1, t2 = t2, timer()
print("Build cmap from HB and SB:", t2-t1, "s")


# create pixel maps
hb_pixel_map, sb_pixel_map, tcks = build_pixel_map(cmap, labels, len(hb_avg_data[:, 0]), numbers_of_lines)


t1, t2 = t2, timer()
print("Create pixel_maps:", t2-t1, "s")


# plot
fig, ax = plt.subplots(figsize=size_of_figure, dpi=300)

hb_cmap = ListedColormap(['#b3cde0', '#6497b1', '#005b96', '#03396c'])
hb_cmap.set_under('k', alpha=0)
sb_cmap = ListedColormap(['#ff7b7b', '#ff5252', '#ff0000', '#a70000'])
sb_cmap.set_under('k', alpha=0)
im1 = ax.imshow(hb_pixel_map.T, cmap=hb_cmap, aspect='auto', interpolation='none', clim=[0.2, 1.0],
                extent=(hb_trj_filename_first-1, hb_trj_filename_last, 0, rows))
im2 = ax.imshow(sb_pixel_map.T, cmap=sb_cmap, aspect='auto', interpolation='none', clim=[0.2, 1.0],
                extent=(hb_trj_filename_first-1, hb_trj_filename_last, 0, rows))
for resid in range(len(h4_1)+1):
    y = np.sum(numbers_of_lines) - np.sum(numbers_of_lines[:resid])
    plt.plot([hb_trj_filename_first-1, hb_trj_filename_last], [y]*2, color='k', linewidth=0.5)

# setup ticks and labels
tcks_pos = []
tcks_seq = []
for i in range(len(h4_1)):
    tcks_pos += [0.5 * np.flip(numbers_of_lines)[i] + np.sum(np.flip(numbers_of_lines)[:i])]
    tcks_seq += [H4_seq[i] + str(i+1)]

tcks_seq.reverse()
tcks.reverse()
ax.set_yticks(tcks_pos)
ax.set_yticklabels(tcks_seq)
ax_right = ax.twinx()
ax_right.set_yticks(np.linspace(0.5, rows + 0.5, rows + 1))
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(axis='both', which='major', labelsize=6)
ax_right.set_yticklabels(tcks)
ax_right.yaxis.tick_right()

ax.set_xlabel('Time, ns')

# save figure
plt.savefig('H4-1_nonuniform.pdf', bbox_inches='tight')
plt.savefig('H4-1_nonuniform.png', bbox_inches='tight')

t1, t2 = t2, timer()
print("Plotting:", t2-t1)


#
# H4-2
#


# build cmap from HB and SB
cmap, labels, numbers_of_lines = build_cmap(hbs, sbs, hb_avg_data, sb_avg_data, "F", len(h4_2))
rows = np.sum(numbers_of_lines)


t1, t2 = t2, timer()
print("Build cmap from HB and SB:", t2-t1, "s")


# create pixel maps
hb_pixel_map, sb_pixel_map, tcks = build_pixel_map(cmap, labels, len(hb_avg_data[:, 0]), numbers_of_lines)


t1, t2 = t2, timer()
print("Create pixel_maps:", t2-t1, "s")


# plot
fig, ax = plt.subplots(figsize=size_of_figure, dpi=300)

hb_cmap = ListedColormap(['#b3cde0', '#6497b1', '#005b96', '#03396c'])
hb_cmap.set_under('k', alpha=0)
sb_cmap = ListedColormap(['#ff7b7b', '#ff5252', '#ff0000', '#a70000'])
sb_cmap.set_under('k', alpha=0)
im1 = ax.imshow(hb_pixel_map.T, cmap=hb_cmap, aspect='auto', interpolation='none', clim=[0.2, 1.0],
                extent=(hb_trj_filename_first-1, hb_trj_filename_last, 0, rows))
im2 = ax.imshow(sb_pixel_map.T, cmap=sb_cmap, aspect='auto', interpolation='none', clim=[0.2, 1.0],
                extent=(hb_trj_filename_first-1, hb_trj_filename_last, 0, rows))
for resid in range(len(h4_2)+1):
    y = np.sum(numbers_of_lines) - np.sum(numbers_of_lines[:resid])
    plt.plot([hb_trj_filename_first-1, hb_trj_filename_last], [y]*2, color='k', linewidth=0.5)

# setup ticks and labels
tcks_pos = []
tcks_seq = []
for i in range(len(h4_2)):
    tcks_pos += [0.5 * np.flip(numbers_of_lines)[i] + np.sum(np.flip(numbers_of_lines)[:i])]
    tcks_seq += [H4_seq[i] + str(i+1)]

tcks_seq.reverse()
tcks.reverse()
ax.set_yticks(tcks_pos)
ax.set_yticklabels(tcks_seq)
ax_right = ax.twinx()
ax_right.set_yticks(np.linspace(0.5, rows + 0.5, rows + 1))
ax_right.set_ylim(ax.get_ylim())
ax_right.tick_params(axis='both', which='major', labelsize=6)
ax_right.set_yticklabels(tcks)
ax_right.yaxis.tick_right()

ax.set_xlabel('Time, ns')

# save figure
plt.savefig('H4-2_nonuniform.pdf', bbox_inches='tight')
plt.savefig('H4-2_nonuniform.png', bbox_inches='tight')

t1, t2 = t2, timer()
print("Plotting:", t2-t1)
