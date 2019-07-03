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


# read data and parameters
with open(hb_trace, 'r') as f:
    header = f.readline()
data = np.genfromtxt(hb_trace, skip_header=1).T
print(data.shape)
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

np.savetxt('taus.csv', acf_fit, fmt='%.6f;%.6f', header='tau, ns; S2')
np.savetxt('acf.txt', acf.T, fmt='%.6e', header=header)
