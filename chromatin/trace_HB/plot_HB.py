import matplotlib.pyplot as plt
import numpy as np
import json
from matplotlib.backends.backend_pdf import PdfPages

# User input:
hb_trace = 'hb_trace.dat'  # file with data about hydrogen bonds
avg_win = 100  # window for averaging of data
step = 25  # traces per page in figures
rid_of_interest = list(range(136, 160)) + list(range(623, 647))

def split_header(header):
    # THR~3::N--DC~1133::OP1 THR~3::OG1--DC~1133::OP1
    result = []
    header = header.rstrip('\n').strip()
    hbs = header.split(' ')
    for hb in hbs:
        donor, acceptor = tuple(hb.split('--'))
        donor_r = donor.split('::')[0]
        donor_aName = donor.split('::')[1]
        donor_rName = donor_r.split('~')[0]
        donor_rId = int(donor_r.split('~')[1])
        acceptor_r = acceptor.split('::')[0]
        acceptor_aName = acceptor.split('::')[1]
        acceptor_rName = acceptor_r.split('~')[0]
        acceptor_rId = int(acceptor_r.split('~')[1])
        result.append([[donor_rName, donor_rId, donor_aName],
                       [acceptor_rName, acceptor_rId, acceptor_aName]])
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


def correct_numbering(hbs, renumber_map=None):
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

    for hb in hbs:
        i1 = hb[0][1]
        i2 = hb[1][1]
        for k, v in absolute_map.items():
            if i1 in v:
                hb[0].insert(0, k)
            if i2 in v:
                hb[1].insert(0, k)
        hb[0][2] = renumber_array[hb[0][2] - 1]
        hb[1][2] = renumber_array[hb[1][2] - 1]

    return hbs


# read data and parameters
with open(hb_trace, 'r') as f:
    header = f.readline()
data = np.genfromtxt(hb_trace, skip_header=1)
with open('input.json', 'r') as f:
    pars = json.load(f)
trj_filename_first = int(pars['trj_filename_first'])
trj_filename_last = int(pars['trj_filename_last'])
stride = int(pars['stride'])


# extract header and average data
hbs = split_header(header)
avg_data = average_data(data, avg_win)


# plot figures for all HBs

cmap = plt.get_cmap('Greys', 5)

with PdfPages('HB_figures.pdf') as pdf:
    for hb_idx in range(0, len(hbs), step):
        plt.figure(figsize=(8, 8), dpi=96)

        # get indexes
        i1 = hb_idx
        i2 = hb_idx + step
        if i2 >= len(hbs):
            i2 = len(hbs)

        # get data
        cur_data = avg_data[:, i1:i2]

        # plot image
        ax = plt.gca()
        im = ax.imshow(cur_data.T, cmap=cmap, aspect='auto', extent=(trj_filename_first-1, trj_filename_last, -0.5, step-0.5))

        # setup ticks and labels
        tcks = ['%s' % hbs[i] for i in range(i1, i2)]
        tcks.reverse()
        plt.yticks(range(step), tcks)

        plt.xlabel('Time, ns')

        # save figure
        pdf.savefig(bbox_inches='tight')
        plt.close()


# plot figures for H4 only

# prune data for not H4 tails
to_del = []
for i, hb in enumerate(hbs):
    if not (hb[0][1] in rid_of_interest) and not (hb[1][1] in rid_of_interest):
        to_del.append(i)
hbs = correct_numbering(hbs)
hbs = [[hb, i] for i, hb in enumerate(hbs) if i not in to_del]
avg_data = np.delete(avg_data, to_del, axis=1)
data = np.delete(data, to_del, axis=1)

with PdfPages('HB_figures_H4.pdf') as pdf:
    plt.figure(figsize=(8, 8), dpi=96)

    # get data
    cur_data = avg_data

    # setup ticks and labels
    tcks = ['%s' % hbs[i] for i in range(len(hbs))]
    tcks.reverse()
    plt.yticks(range(len(hbs)), tcks)
    plt.xlabel('Time, ns')

    # plot image
    ax = plt.gca()
    im = ax.imshow(cur_data.T, cmap=cmap, aspect='auto', extent=(trj_filename_first-1, trj_filename_last, -0.5, len(hbs)-0.5))

    # save figure
    pdf.savefig(bbox_inches='tight')
    plt.close()





# plot traces of interaction partners

# H4-1
# get unique resids in H4-1 that has HB
has_hb_1 = sorted(list(set([hb[0][0][2] for hb in hbs if hb[0][0][0] == 'B'])))

# group hb by common interaction partners from H4-1
hbs_grouped = []
for resid in has_hb_1:
    p = [hb[1] for hb in hbs if ((hb[0][0][0] == 'B') & (hb[0][0][2] == resid))]
    hbs_grouped.append(p)

# get avg_data traces for each group
traces = []
for i, resid in enumerate(has_hb_1):
    idxs = []
    for j, hb in enumerate(hbs):
        if hb[1] in hbs_grouped[i]:
            idxs.append(j)
    traces.append(avg_data[:, idxs])

# calculate functions for plotting
fn = []
for trace, group in zip(traces, hbs_grouped):
    vals = np.zeros(len(trace), dtype=int)
    for i, t in enumerate(trace):
        # fill vals: 0 - no HB; 1 - HB with first partner; 2 - HB with second partner ...
        if np.sum(t) == 0.0:
            vals[i] = 0
            continue
        for j, s in enumerate(t):
            if s > 0.5:
                vals[i] = j + 1
                # what if hb is formed simultaneously to two partners? Take first only since it is a rare case
                continue
    fn.append(vals)

# plot traces
x = np.linspace(trj_filename_first, trj_filename_last, len(fn[0]))
with PdfPages('HB_figures_H4_traces_H4-1.pdf') as pdf:
    for i, y in enumerate(fn):
        plt.figure(figsize=(10, 1))
        plt.scatter(x, y, marker='o', linewidths=0, alpha=0.7)
        plt.xlabel('Time, ns')
        vals = range(max(y) + 1)
        plt.yticks(vals, ['no HB'] + [str(hb) for hb in hbs if hb[1] in hbs_grouped[i]])
        pdf.savefig(bbox_inches='tight')
        plt.close()






# H4-2
# get unique resids in H4-2 that has HB
has_hb_2 = sorted(list(set([hb[0][0][2] for hb in hbs if hb[0][0][0] == 'F'])))

# group hb by common interaction partners from H4-1
hbs_grouped = []
for resid in has_hb_2:
    p = [hb[1] for hb in hbs if ((hb[0][0][0] == 'F') & (hb[0][0][2] == resid))]
    hbs_grouped.append(p)

# get avg_data traces for each group
traces = []
for i, resid in enumerate(has_hb_2):
    idxs = []
    for j, hb in enumerate(hbs):
        if hb[1] in hbs_grouped[i]:
            idxs.append(j)
    traces.append(avg_data[:, idxs])

# calculate functions for plotting
fn = []
for trace, group in zip(traces, hbs_grouped):
    vals = np.zeros(len(trace), dtype=int)
    for i, t in enumerate(trace):
        # fill vals: 0 - no HB; 1 - HB with first partner; 2 - HB with second partner ...
        if np.sum(t) == 0.0:
            vals[i] = 0
            continue
        for j, s in enumerate(t):
            if s > 0.5:
                vals[i] = j + 1
                # what if hb is formed simultaneously to two partners? Take first only since it is a rare case
                continue
    fn.append(vals)

# plot traces
x = np.linspace(trj_filename_first, trj_filename_last, len(fn[0]))
with PdfPages('HB_figures_H4_traces_H4-2.pdf') as pdf:
    for i, y in enumerate(fn):
        plt.figure(figsize=(10, 1))
        plt.scatter(x, y, marker='o', linewidths=0, alpha=0.7)
        plt.xlabel('Time, ns')
        vals = range(max(y) + 1)
        plt.yticks(vals, ['no HB'] + [str(hb) for hb in hbs if hb[1] in hbs_grouped[i]])
        pdf.savefig(bbox_inches='tight')
        plt.close()
