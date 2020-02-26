from correlation_functions import Corfun
import numpy as np
import os

n_exp = ___N_EXP___
repetitions = ___REP___
step = 1e-12
NMR_freq = 850e6
tumbling_tau_exp = 165.1e-9  # 165.1 ns was obtained by HYDRONMR for 3LZ0 trajectory (100 mM NaCl in 100% H2O)
fit_out = "autocor_fit_results.txt"
R1_out = "R1_results.txt"
R2_out = "R2_results.txt"

fun_aligned = Corfun("cor_NH____FIRST_DAT_FILE___-___LAST_DAT_FILE____diluted")

# Multiply ACFs for aligned frames by tumbling exponential and save
time_in_steps = fun_aligned.data[0][:, 0]
exp_tau_exp = np.exp(-time_in_steps/(tumbling_tau_exp / step))
for i, acf in enumerate(fun_aligned.data):
    fun_aligned.data[i][:, 1] = np.multiply(acf[:, 1], exp_tau_exp)
tumbling_on_dir = "cor_NH____FIRST_DAT_FILE___-___LAST_DAT_FILE____tumbling_on"
if not os.path.exists(tumbling_on_dir):
    os.makedirs(tumbling_on_dir)
fun_aligned.write(tumbling_on_dir)

# Process ACF
fun = fun_aligned

up_mask = 1e-9 * (___LAST_DAT_FILE___ - ___FIRST_DAT_FILE___ + 1) / step
lb = [0.0]*2*n_exp
ub = [1.0]*n_exp + [up_mask]*n_exp

# Fit autocorrelation functions for aligned trajectory
sp = np.zeros(n_exp*2)
sp[:n_exp] = np.array([1.0/n_exp]*n_exp)
sp[n_exp:] = np.logspace(np.log10(1e-12), np.log10(tumbling_tau_exp), n_exp)/step
fun.fit(n_exp=n_exp, rep=repetitions, step=step, sp=sp, lb=lb, ub=ub, fout=fit_out)

# Calculate relaxation rates
fun.calc_R1(opfreq=NMR_freq, X='N15', fnout=R1_out)
fun.calc_R2(opfreq=NMR_freq, X='N15', fnout=R2_out)
