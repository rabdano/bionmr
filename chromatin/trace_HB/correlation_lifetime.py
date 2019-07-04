import json
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys


hb_trace = 'hb_trace.dat'  # file with data about hydrogen bonds
remove_first_acf_point = True


def autocorr(x):
    # return real part of autocorrelation function
    f = np.fft.fft(np.pad(x, len(x), mode='constant')[len(x):])
    result = np.fft.ifft(f * np.conj(f))
    result = result[:len(x)]
    result /= np.linspace(len(x), 1, len(x))
    return np.real(result)


def exp_decay(x, a, tau, c):
    return a * np.exp(-x/tau) + c


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
header = header.rstrip('\n').strip()
hbs = split_header(header)
hbs = correct_numbering(hbs)
data = np.genfromtxt(hb_trace, skip_header=1).T
with open('input.json', 'r') as f:
    pars = json.load(f)
trj_filename_first = pars['trj_filename_first']
trj_filename_last = pars['trj_filename_last']

# calculate time axis vector
time = np.linspace(0, (trj_filename_last-trj_filename_first+1)*1000, len(data[0]))

# calculate ACF for each HB
acf = np.zeros(data.shape)
acf_fit = np.zeros((len(data), 2))

for i, trace in enumerate(data):
    sys.stdout.write('Processing hydrogen bond #{} of {}\r'.format(i+1, len(data)))

    trace_acf = autocorr(trace)

    # remove first point if necessary
    if remove_first_acf_point:
        x = time[1:]
        y = trace_acf[1:]
    else:
        x = time
        y = trace_acf

    # cut acf for fitting (0.8)
    x = x[:int(np.floor(len(x) * 0.8))]
    y = y[:int(np.floor(len(y) * 0.8))]

    popt, pcov = curve_fit(exp_decay, x, y, p0=(np.max(y), 1000, 0),
                           bounds=([np.max(y)/10, np.min(x)/10, -1e-15], [np.max(y)*10, np.max(x)*10, np.max(y)]),
                           xtol=1e-4, max_nfev=10000)

    acf_fit[i, 0] = popt[1]/1000
    acf_fit[i, 1] = np.abs(popt[2])/(popt[0] + popt[2])
    acf[i] = trace_acf

    plt.figure(figsize=(8, 6))
    plt.plot(x, y)
    x_fit = np.linspace(np.min(time), np.max(time), 1000)
    plt.plot(x_fit, exp_decay(x_fit, *popt))
    plt.xlabel('Time, ps')
    plt.ylabel('Correlation')
    plt.savefig('Figures/{:04d}.png'.format(i), bbox_inches='tight')
    plt.close()
sys.stdout.write('\n')

with open('result.csv', 'w') as f:
    h = 'bond_id;'
    h += 'donor_chain;donor_resname;donor_resid;donor_atom;'
    h += 'acceptor_chain;acceptor_resname;acceptor_resid;acceptor_atom;'
    h += 'tau, ns; S2\n'
    f.write(h)
    i = 1
    fmt = '{:d}; {};{};{:d};{}; {};{};{:d};{}; {:.2f};{:.3f}\n'
    for hb, r in zip(hbs, acf_fit):
        f.write(fmt.format(i, hb[0][0], hb[0][1], hb[0][2], hb[0][3],
                              hb[1][0], hb[1][1], hb[1][2], hb[1][3],
                              r[0], r[1]))
        i += 1
# np.savetxt('taus.csv', acf_fit, fmt='%.6f;%.6f', header='tau, ns; S2')
# np.savetxt('acf.txt', acf.T, fmt='%.6e', header=header)
