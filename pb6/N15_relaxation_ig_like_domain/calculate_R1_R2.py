from correlation_functions import Corfun
import numpy as np

n_exp = 6
repetitions = 1
step = 1e-12
NMR_freq = 800e6
tumbling_tau_exp = 3.8e-9  # 3.8 ns is mean tumbling time obtained from 1 exp fitting of not aligned trajectory ACFs
expected_tau_r_range = (1e-12, 4e-9)
fit_out = "autocor_fit_results.txt"
R1_out = "R1_results.txt"
R2_out = "R2_results.txt"

fun = Corfun("cor_NH_1-1000_diluted")

up_mask = 1e15
lb = [0.0]*2*n_exp
ub = [1.0]*n_exp + [up_mask]*n_exp

# Obtain autocorrelation functions for aligned trajectory
sp = np.zeros(n_exp*2)
sp[:n_exp] = np.array([1.0/n_exp]*n_exp)
sp[n_exp:] = np.logspace(np.log10(expected_tau_r_range[0]), np.log10(expected_tau_r_range[1]), n_exp)/step

fun.fit(n_exp=n_exp, rep=repetitions, step=step, sp=sp, lb=lb, ub=ub, fout=fit_out)

print(fun.pars[:][n_exp:n_exp*2+1])

for i in range(len(fun.pars)):
    for j in range(n_exp, n_exp*2):
        fun.pars[i][j] = (fun.pars[i][j] * tumbling_tau_exp) / (fun.pars[i][j] + tumbling_tau_exp)

fun.calc_R1(opfreq=NMR_freq, X='N15', fnout=R1_out)
fun.calc_R2(opfreq=NMR_freq, X='N15', fnout=R2_out)
