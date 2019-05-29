import nmrglue as ng
import numpy as np
from scipy.optimize import curve_fit, fmin
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
from glob import glob


def lorentzian(x, amp, x0, w):
    return amp * (1 / (1 + ((x0 - x) / (w / 2))**2))


def T2_decay(x, T2):
    return np.exp(-x/T2)


def T1_decay(x, T1, M0=-2.0):
    return (1 + M0 * np.exp(-x/T1))


def phase_correction(data, fn, p0=0.0, p1=0.0):
    if not callable(fn):
        fn = {
            'peak_minima': ng.proc_autophase._ps_peak_minima_score,
            'acme': ng.proc_autophase._ps_acme_score,
        }[fn]

    opt = [p0, p1]
    opt = fmin(fn, x0=opt, args=(data, ), disp=False)

    phasedspc = ng.proc_autophase.ps(data, p0=opt[0], p1=opt[1])

    return opt


def process_T1(path_to_dataset, fig_out='T1.pdf'):
    # Load data
    dic, data = ng.bruker.read(path_to_dataset)
    sw_hz = dic['acqus']['SW_h']
    o1 = dic['acqus']['O1']
    delays = np.genfromtxt(path_to_dataset + '/vdlist')

    # Processing parameters
    n_points = 2 ** 14
    cut_first = 70
    lb = 3.0  # Hz

    # Fourier transform
    data = ng.proc_base.em(data, lb=lb / sw_hz)
    data = ng.proc_base.zf(data, n_points - data.shape[1])
    data = ng.proc_base.fft(data[:, cut_first:])

    # Phase and baseline correction
    # get phase for last spectrum
    ph0, ph1 = phase_correction(data[-1], 'acme')
    # apply phase correction to all spectra
    data = ng.proc_autophase.ps(data, p0=ph0, p1=ph1)
    nl = [i for i in range(0, int((n_points - cut_first) / 4))] + \
         [i for i in range(3 * int((n_points - cut_first) / 4), (n_points - cut_first))]
    data = ng.proc_bl.base(data.real, nl=nl)

    # fit
    hz = np.linspace(o1 - sw_hz / 2, o1 + sw_hz / 2, data.shape[1])
    pos = np.argmin(data[0])
    left = pos - int(data.shape[1] / 40)
    right = pos + int(data.shape[1] / 40)

    relax_data = np.zeros((data.shape[0],))

    pp = PdfPages('spec' + fig_out + '_T1.pdf')
    plt.figure()
    for i, spec in enumerate(data):
        x = hz[left:right]
        y = spec[left:right]
        plt.clf()
        plt.plot(x, y, label=delays[i])
        plt.xlabel('Frequency, Hz')
        plt.legend()
        if np.abs(np.max(y)) > np.abs(np.min(y)):
            amp0 = np.max(y)
        else:
            amp0 = np.min(y)
        popt, pcov = curve_fit(lorentzian, x, y, p0=(amp0, hz[pos], 2), maxfev=10000,
                               bounds=([-1e10, hz[pos]-1, 1], [1e10, hz[pos]+1, 10]),
                               xtol=1e-4)
        relax_data[i] = popt[0] * np.pi * popt[2]
        pp.savefig()
    pp.close()
    relax_data = relax_data / np.max(relax_data)

    # Plot relaxation curve and fit
    popt, pcov = curve_fit(T1_decay, delays, relax_data, p0=0.5)
    xfit = np.linspace(np.min(delays), np.max(delays), 1000)
    yfit = T1_decay(xfit, popt[0])
    residual = relax_data - T1_decay(delays, popt[0])

    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    ax0 = plt.subplot(gs[0])
    ax0.scatter(delays, relax_data, c='k')
    plt.plot(xfit, yfit, 'r')
    plt.grid()
    ax0.set_ylabel(r'M$_z$')
    ax0.set_title('R$_1$=%.3fs$^{-1}$ T$_1$=%.3fs' % (1 / popt[0], popt[0]))
    ax1 = plt.subplot(gs[1])
    ax1.set_xlabel('Delay, s')
    ax1.set_ylabel('Residual')
    ax1.scatter(delays, residual)
    ax1.set_ylim([-np.abs(np.max(residual)) * 1.2, np.abs(np.max(residual)) * 1.2])
    plt.grid()

    plt.tight_layout()
    plt.savefig(fig_out + '_T1.pdf', bbox_inches='tight')

    return 1 / popt[0], popt[0]


def process_T2(path_to_dataset, fig_out='T2.pdf'):
    # Load data
    dic, data = ng.bruker.read(path_to_dataset)
    sw_hz = dic['acqus']['SW_h']
    o1 = dic['acqus']['O1']
    d20 = dic['acqus']['D'][20]
    p2 = dic['acqus']['P'][2]
    delays = np.genfromtxt(path_to_dataset + '/vclist')
    delays = delays * (2 * d20 + p2 * 1e-6)

    # Processing parameters
    n_points = 2 ** 14
    cut_first = 70
    lb = 3.0  # Hz

    # Fourier transform
    data = ng.proc_base.em(data, lb=lb / sw_hz)
    data = ng.proc_base.zf(data, n_points - data.shape[1])
    data = ng.proc_base.fft(data[:, cut_first:])

    # Phase and baseline correction
    # get phase for first spectrum
    ph0, ph1 = phase_correction(data[0], 'acme')
    # apply phase correction to all spectra
    data = ng.proc_autophase.ps(data, p0=ph0, p1=ph1)
    nl = [i for i in range(0, int((n_points - cut_first) / 4))] + \
         [i for i in range(3 * int((n_points - cut_first) / 4), (n_points - cut_first))]
    data = ng.proc_bl.base(data.real, nl=nl)

    # fit
    hz = np.linspace(o1 - sw_hz / 2, o1 + sw_hz / 2, data.shape[1])
    pos = np.argmax(data[0])
    left = pos - int(data.shape[1] / 40)
    right = pos + int(data.shape[1] / 40)

    relax_data = np.zeros((data.shape[0],))

    pp = PdfPages('spec' + fig_out + '_T2.pdf')
    plt.figure()
    for i, spec in enumerate(data):
        x = hz[left:right]
        y = spec[left:right]
        plt.clf()
        plt.plot(x, y, label=delays[i])
        plt.xlabel('Frequency, Hz')
        plt.legend()
        if np.abs(np.max(y)) > np.abs(np.min(y)):
            amp0 = np.max(y)
        else:
            amp0 = np.min(y)
        popt, pcov = curve_fit(lorentzian, x, y, p0=(amp0, hz[pos], 2), maxfev=10000,
                               bounds=([-1e10, hz[pos]-1, 1], [1e10, hz[pos]+1, 10]),
                               xtol=1e-4)
        relax_data[i] = popt[0] * np.pi * popt[2]
        pp.savefig()
    pp.close()
    relax_data = relax_data / np.max(relax_data)

    # Plot relaxation curve and fit
    popt, pcov = curve_fit(T2_decay, delays, relax_data, p0=0.5)
    xfit = np.linspace(np.min(delays), np.max(delays), 1000)
    yfit = T2_decay(xfit, popt[0])
    residual = relax_data - T2_decay(delays, popt[0])

    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    ax0 = plt.subplot(gs[0])
    ax0.scatter(delays, relax_data, c='k')
    plt.plot(xfit, yfit, 'r')
    plt.grid()
    ax0.set_ylabel(r'M$_z$')
    ax0.set_title('R$_2$=%.3fs$^{-1}$ T$_2$=%.3fs' % (1 / popt[0], popt[0]))
    ax1 = plt.subplot(gs[1])
    ax1.set_xlabel('Delay, s')
    ax1.set_ylabel('Residual')
    ax1.scatter(delays, residual)
    ax1.set_ylim([-np.abs(np.max(residual)) * 1.2, np.abs(np.max(residual)) * 1.2])
    plt.grid()

    plt.tight_layout()
    plt.savefig(fig_out + '_T2.pdf', bbox_inches='tight')

    return 1 / popt[0], popt[0]


# datasets
T1_datasets = glob('../RSO2h_T1/*')
T2_datasets = glob('../RSO2h_T2/*')
T1_names = [d.split('\\')[-1] for d in T1_datasets]
T2_names = [d.split('\\')[-1] for d in T2_datasets]

# process and write results
with open('results_T1.txt', 'w') as f:
    f.write('Dataset, R1, T1\n')
    for dataset, name in zip(T1_datasets, T1_names):
        R1, T1 = process_T1(dataset, fig_out=name)
        f.write(dataset + ',%.3f,%.3f\n' % (R1, T1))
with open('results_T2.txt', 'w') as f:
    f.write('Dataset, R2, T2\n')
    for dataset, name in zip(T2_datasets, T2_names):
        R2, T2 = process_T2(dataset, fig_out=name)
        f.write(dataset + ',%.3f,%.3f\n' % (R2, T2))



