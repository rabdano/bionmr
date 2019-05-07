import matplotlib.pyplot as plt
import numpy as np
import json
from matplotlib.backends.backend_pdf import PdfPages

# User input:
sb_trace = 'sb_trace.dat'  # file with data about hydrogen bonds
avg_win = 100  # window for averaging of data
step = 25  # traces per page in figures
rid_of_interest = list(range(136, 160)) + list(range(623, 647))

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
        result[int(start/window), :] = np.around(np.mean(data[start:start+window, :], axis=0))
    return result


# read data and parameters
with open(sb_trace, 'r') as f:
    header = f.readline()
data = np.genfromtxt(sb_trace, skip_header=1)
with open('input.json', 'r') as f:
    pars = json.load(f)
trj_filename_first = pars['trj_filename_first']
trj_filename_last = pars['trj_filename_last']
stride = pars['stride']


# extract header and average data
sbs = split_header(header)
avg_data = average_data(data, avg_win)


# plot figures for all HBs

cmap = plt.get_cmap('binary', 2)

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
        im = ax.imshow(cur_data.T, cmap=cmap, aspect='auto')

        # setup ticks and labels
        plt.yticks(range(step), ['%s' % sbs[i] for i in range(i1, i2)])
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
sbs = [sb for i, sb in enumerate(sbs) if i not in to_del]
avg_data = np.delete(avg_data, to_del, axis=1)

with PdfPages('SB_figures_H4.pdf') as pdf:
    plt.figure(figsize=(8, 8), dpi=96)

    # get data
    cur_data = avg_data

    # setup ticks and labels
    plt.yticks(range(len(sbs)), ['%s' % sbs[i] for i in range(len(sbs))])
    plt.xlabel('Time, ns')

    # plot image
    ax = plt.gca()
    im = ax.imshow(cur_data.T, cmap=cmap, aspect='auto')

    # save figure
    pdf.savefig(bbox_inches='tight')
    plt.close()
