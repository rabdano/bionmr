import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json
from matplotlib.backends.backend_pdf import PdfPages

# User input:
sb_trace = 'sb_trace.dat'  # file with data about hydrogen bonds
avg_win = 500  # window for averaging of data
step = 25  # traces per page in figures
h4_1 = list(range(136, 160))
h4_2 = list(range(623, 647))
rid_of_interest = h4_1 + h4_2
H4_seq = 'SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG'

def split_header(header):
    result = []
    header = header.rstrip('\n').strip()
    sbs = header.split(' ')
    for sb in sbs:
        donor, acceptor = tuple(sb.split('--'))
        donor_rName = donor.split('~')[0]
        donor_rId = int(donor.split('~')[1])
        acceptor_rName = acceptor.split('~')[0]
        acceptor_rId = int(acceptor.split('~')[1])
        result.append([[donor_rName, donor_rId],
                       [acceptor_rName, acceptor_rId]])
    return result


def average_data(data, window=100):
    # average data over windows
    # e.g. if in window hydrogen bond is present > 50% then value for window is 1
    n, m = data.shape
    nw = int(n/window)
    result = np.zeros((nw, m))
    for start in range(0, len(data), window):
        result[int(start/window), :] = np.mean(data[start:start+window, :], axis=0)
    return result


def correct_numbering(sbs, renumber_map=None):
    if not renumber_map:
        renumber_map = {'A': [rid for rid in range(1, 136)],
                        'B': [rid for rid in range(1, 103)],
                        'C': [rid for rid in range(1, 129)],
                        'D': [rid for rid in range(1, 123)],
                        'E': [rid for rid in range(1, 136)],
                        'F': [rid for rid in range(1, 103)],
                        'G': [rid for rid in range(1, 129)],
                        'H': [rid for rid in range(1, 123)],
                        'I': [rid for rid in range(1, 148)],
                        'J': [rid for rid in range(1, 148)]}

    absolute_map = {}
    renumber_array = []
    last_resid = 0

    for k in sorted(renumber_map):
        v = renumber_map[k]
        absolute_map[k] = [rid for rid in range(last_resid + 1, last_resid + len(v) + 1)]
        renumber_array += v
        last_resid = len(renumber_array)

    for sb in sbs:
        i1 = sb[0][1]
        i2 = sb[1][1]
        for k, v in absolute_map.items():
            if i1 in v:
                sb[0].insert(0, k)
            if i2 in v:
                sb[1].insert(0, k)
        sb[0][2] = renumber_array[sb[0][2] - 1]
        sb[1][2] = renumber_array[sb[1][2] - 1]

    return sbs


# read data and parameters
with open(sb_trace, 'r') as f:
    header = f.readline()
data = np.genfromtxt(sb_trace, skip_header=1)
with open('input.json', 'r') as f:
    pars = json.load(f)
trj_filename_first = max(int(pars['trj_filename_first']), int(pars['start_frame_of_interest']/1000)+1)
trj_filename_last = min(int(pars['trj_filename_last']), int(pars['end_frame_of_interest']/1000+1))
stride = int(pars['stride'])


# extract header and average data
sbs = split_header(header)
avg_data = average_data(data, avg_win)


# plot figures for all HBs

cmap = plt.get_cmap('Greys', 5)

with PdfPages('SB_figures.pdf') as pdf:
    for sb_idx in range(0, len(sbs), step):
        plt.figure(figsize=(8, 8), dpi=96)

        # get indexes
        i1 = sb_idx
        i2 = sb_idx + step
        if i2 >= len(sbs):
            i2 = len(sbs)

        # get data
        cur_data = avg_data[:, i1:i2]

        # plot image
        ax = plt.gca()
        im = ax.imshow(cur_data.T, cmap=cmap, aspect='auto', extent=(trj_filename_first-1, trj_filename_last, -0.5, step-0.5))

        # setup ticks and labels
        tcks = ['%s' % sbs[i] for i in range(i1, i2)]
        tcks.reverse()
        plt.yticks(range(step), tcks)

        plt.xlabel('Time, ns')

        # save figure
        pdf.savefig(bbox_inches='tight')
        plt.close()


# plot figures for H4 only

# prune data for not H4 tails
to_del = []
for i, sb in enumerate(sbs):
    if not (sb[0][1] in rid_of_interest) and not (sb[1][1] in rid_of_interest):
        to_del.append(i)
sbs = correct_numbering(sbs)
sbs = [[sb, i] for i, sb in enumerate(sbs) if i not in to_del]
avg_data = np.delete(avg_data, to_del, axis=1)
data = np.delete(data, to_del, axis=1)

# prune non protein-DNA SB data
to_del = []
for i, sb in enumerate([s[0] for s in sbs]):
    if not (((sb[0][0] in ['B', 'F']) and (sb[1][0] in ['I', 'J'])) or
            ((sb[0][0] in ['I', 'J']) and (sb[1][0] in ['B', 'F']))):
        to_del.append(i)
sbs = [[sb, i] for i, sb in enumerate([s[0] for s in sbs]) if i not in to_del]
sb_avg_data = np.delete(avg_data, to_del, axis=1)
sb_data = np.delete(data, to_del, axis=1)


sbs_1 = [sb for sb in sbs if ((sb[0][0][0] == 'B') | (sb[0][1][0] == 'B'))]
sbs_2 = [sb for sb in sbs if ((sb[0][0][0] == 'F') | (sb[0][1][0] == 'F'))]
avg_data_1 = avg_data[:, :len(sbs_1)]
avg_data_2 = avg_data[:, len(sbs_1):]

h4_1_cur = [s[0][0][2] for s in sbs_1]
h4_2_cur = [s[0][0][2] for s in sbs_2]


# H4-1
if len(sbs_1) > 0:
    labels_1 = []
    h4_1_traces = []
    counter = 0
    for resid in range(1, len(h4_1)+1):
        resname = H4_seq[resid-1]
        if resid in h4_1_cur:
            for sb in sbs_1:
                if sb[0][0][2] == resid:
                    labels_1.append(resname + '-{} '.format(resid) + '({})'.format(sb[0][0][0]) +
                                    ' - {}'.format(sb[0][1][1]) + '-{}'.format(sb[0][1][2]) + ' ({})'.format(sb[0][1][0]))
                    h4_1_traces.append(avg_data_1[:, counter])
                    counter += 1
        else:
            labels_1.append('{}-{} ({})'.format(resname, resid, 'B'))
            h4_1_traces.append(np.zeros(avg_data_1[:, 0].shape))

    with PdfPages('SB_figures_H4-1.pdf') as pdf:
        fig, ax = plt.subplots(figsize=(8, 8), dpi=96)

        # get data
        cur_data = np.array(h4_1_traces)

        # setup ticks and labels
        tcks = labels_1
        tcks.reverse()
        plt.yticks(range(len(labels_1)), tcks)
        ax.yaxis.tick_right()
        plt.xlabel('Time, ns')

        # plot image
        ax = plt.gca()
        im = ax.imshow(cur_data, cmap=cmap, aspect='auto', extent=(trj_filename_first-1, trj_filename_last, -0.5, len(labels_1)-0.5))

        # save figure
        pdf.savefig(bbox_inches='tight')
        plt.close()


# H4-2
if len(sbs_2) > 0:
    labels_2 = []
    h4_2_traces = []
    counter = 0
    for resid in range(1, len(h4_2)+1):
        resname = H4_seq[resid-1]
        if resid in h4_2_cur:
            for sb in sbs_2:
                if sb[0][0][2] == resid:
                    labels_2.append(resname + '-{} '.format(resid) + '({})'.format(sb[0][0][0]) +
                                    ' - {}'.format(sb[0][1][1]) + '-{}'.format(sb[0][1][2]) + ' ({})'.format(sb[0][1][0]))
                    h4_2_traces.append(avg_data_2[:, counter])
                    counter += 1
        else:
            labels_2.append('{}-{} ({})'.format(resname, resid, 'F'))
            h4_2_traces.append(np.zeros(avg_data_2[:, 0].shape))

    with PdfPages('SB_figures_H4-2.pdf') as pdf:
        fig, ax = plt.subplots(figsize=(8, 8), dpi=96)

        # get data
        cur_data = np.array(h4_2_traces)

        # setup ticks and labels
        tcks = labels_2
        tcks.reverse()
        plt.yticks(range(len(labels_2)), tcks)
        ax.yaxis.tick_right()
        plt.xlabel('Time, ns')

        # plot image
        ax = plt.gca()
        im = ax.imshow(cur_data, cmap=cmap, aspect='auto', extent=(trj_filename_first-1, trj_filename_last, -0.5, len(labels_2)-0.5))

        # save figure
        pdf.savefig(bbox_inches='tight')
        plt.close()


# plot traces of interaction partners
# H4-1
# get unique resids in H4-1 that has SB
has_sb_1 = sorted(list(set([sb[0][0][2] for sb in sbs if sb[0][0][0] == 'B'])))

# group sb by common interaction partners from H4-1
sbs_grouped = []
for resid in has_sb_1:
    p = [sb[1] for sb in sbs if ((sb[0][0][0] == 'B') & (sb[0][0][2] == resid))]
    sbs_grouped.append(p)

# get avg_data traces for each group
traces = []
for i, resid in enumerate(has_sb_1):
    idxs = []
    for j, sb in enumerate(sbs):
        if sb[1] in sbs_grouped[i]:
            idxs.append(j)
    traces.append(avg_data[:, idxs])

# calculate functions for plotting
fn = []
for trace, group in zip(traces, sbs_grouped):
    vals = np.zeros(len(trace), dtype=int)
    for i, t in enumerate(trace):
        # fill vals: 0 - no SB; 1 - SB with first partner; 2 - SB with second partner ...
        if np.sum(t) == 0.0:
            vals[i] = 0
            continue
        for j, s in enumerate(t):
            if s > 0.5:
                vals[i] = j + 1
                # what if sb is formed simultaneously to two partners? Take first only since it is a rare case
                continue
    fn.append(vals)

# plot traces
if len(fn) > 0:
    x = np.linspace(trj_filename_first, trj_filename_last, len(fn[0]))
    with PdfPages('SB_figures_H4_traces_H4-1.pdf') as pdf:
        for i, y in enumerate(fn):
            plt.figure(figsize=(10, 1))
            plt.scatter(x, y, marker='o', linewidths=0, alpha=0.7)
            plt.xlabel('Time, ns')
            vals = range(max(y) + 1)
            plt.yticks(vals, ['no SB'] + [str(sb) for sb in sbs if sb[1] in sbs_grouped[i]])
            pdf.savefig(bbox_inches='tight')
            plt.close()






# H4-2
# get unique resids in H4-2 that has SB
has_sb_2 = sorted(list(set([sb[0][0][2] for sb in sbs if sb[0][0][0] == 'F'])))

# group sb by common interaction partners from H4-1
sbs_grouped = []
for resid in has_sb_2:
    p = [sb[1] for sb in sbs if ((sb[0][0][0] == 'F') & (sb[0][0][2] == resid))]
    sbs_grouped.append(p)

# get avg_data traces for each group
traces = []
for i, resid in enumerate(has_sb_2):
    idxs = []
    for j, sb in enumerate(sbs):
        if sb[1] in sbs_grouped[i]:
            idxs.append(j)
    traces.append(avg_data[:, idxs])

# calculate functions for plotting
fn = []
for trace, group in zip(traces, sbs_grouped):
    vals = np.zeros(len(trace), dtype=int)
    for i, t in enumerate(trace):
        # fill vals: 0 - no SB; 1 - SB with first partner; 2 - SB with second partner ...
        if np.sum(t) == 0.0:
            vals[i] = 0
            continue
        for j, s in enumerate(t):
            if s > 0.5:
                vals[i] = j + 1
                # what if sb is formed simultaneously to two partners? Take first only since it is a rare case
                continue
    fn.append(vals)

# plot traces
if len(fn) > 0:
    x = np.linspace(trj_filename_first, trj_filename_last, len(fn[0]))
    with PdfPages('SB_figures_H4_traces_H4-2.pdf') as pdf:
        for i, y in enumerate(fn):
            plt.figure(figsize=(10, 1))
            plt.scatter(x, y, marker='o', linewidths=0, alpha=0.7)
            plt.xlabel('Time, ns')
            vals = range(max(y) + 1)
            plt.yticks(vals, ['no SB'] + [str(sb) for sb in sbs if sb[1] in sbs_grouped[i]])
            pdf.savefig(bbox_inches='tight')
            plt.close()
