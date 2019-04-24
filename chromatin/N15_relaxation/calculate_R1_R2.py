from correlation_functions import Corfun
import numpy as np

n_exp = ___N_EXP___
repetitions = ___REP___
step = 1e-12
NMR_freq = 500e6
tumbling_tau_exp = 185.9 # 176.6 ns was obtained by HYDRONMR for mononucleosome
fit_out = "autocor_fit_results.txt"
R1_out = "R1_results.txt"
R2_out = "R2_results.txt"

fun = Corfun("cor_NH____FIRST_DAT_FILE___-___LAST_DAT_FILE____diluted")

up_mask = 1e15
lb = [0.0]*2*n_exp
ub = [1.0]*n_exp + [up_mask]*n_exp

# Obtain autocorrelation functions for aligned trajectory
sp = np.zeros(n_exp*2)
sp[:n_exp] = np.array([1.0/n_exp]*n_exp)
sp[n_exp:] = np.logspace(np.log10(1e-12),np.log10(tumbling_tau_exp*1e-9),n_exp)/step

fun.fit(n_exp=n_exp, rep=repetitions, sp=sp, lb=lb, ub=ub, fout=fit_out)

for i in range(len(fun.pars)):
    for j in range(n_exp, n_exp*2):
        fun.pars[i][j] = (fun.pars[i][j] * tumbling_tau_exp*1000/(step*1e12)) \
                       / (fun.pars[i][j] + tumbling_tau_exp*1000/(step*1e12))

fun.calc_R1(step, opfreq=NMR_freq, X='N15', fnout=R1_out)
fun.calc_R2(step, opfreq=NMR_freq, X='N15', fnout=R2_out)
