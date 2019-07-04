import json
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys


sb_trace = 'sb_trace.dat'  # file with data about hydrogen bonds
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
header = header.rstrip('\n').strip()
sbs = split_header(header)
sbs = correct_numbering(sbs)
data = np.genfromtxt(sb_trace, skip_header=1).T
with open('input.json', 'r') as f:
    pars = json.load(f)
trj_filename_first = pars['trj_filename_first']
trj_filename_last = pars['trj_filename_last']

# calculate time axis vector
time = np.linspace(0, (trj_filename_last-trj_filename_first+1)*1000, len(data[0]))

# calculate ACF for each SB
acf = np.zeros(data.shape)
acf_fit = np.zeros((len(data), 2))

for i, trace in enumerate(data):
    sys.stdout.write('Processing salt bridge #{} of {}\r'.format(i+1, len(data)))

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
    plt.plot(x, y, label=str(sbs[i]))
    x_fit = np.linspace(np.min(time), np.max(time), 1000)
    plt.plot(x_fit, exp_decay(x_fit, *popt), label=r'$\tau$={:.2f}ns; S$^2$={:.3f}'.format(acf_fit[i, 0], acf_fit[i, 1]))
    plt.xlabel('Time, ps')
    plt.ylabel('Correlation')
    plt.legend()
    plt.savefig('Figures/{:04d}.png'.format(i), bbox_inches='tight')
    plt.close()
sys.stdout.write('\n')

with open('result.csv', 'w') as f:
    h = 'bridge_id;'
    h += 'donor_chain;donor_resname;donor_resid;'
    h += 'acceptor_chain;acceptor_resname;acceptor_resid;'
    h += 'tau, ns; S2\n'
    f.write(h)
    i = 1
    fmt = '{:d}; {};{};{:d}; {};{};{:d}; {:.2f};{:.3f}\n'
    for hb, r in zip(sbs, acf_fit):
        f.write(fmt.format(i, hb[0][0], hb[0][1], hb[0][2],
                              hb[1][0], hb[1][1], hb[1][2],
                              r[0], r[1]))
        i += 1
# np.savetxt('taus.csv', acf_fit, fmt='%.6f;%.6f', header='tau, ns; S2')
# np.savetxt('acf.txt', acf.T, fmt='%.6e', header=header)
