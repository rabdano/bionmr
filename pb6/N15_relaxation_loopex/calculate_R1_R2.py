from correlation_functions import Corfun
import numpy as np
import os
import json


#
#  setup trajectory parameters
#
with open("input.json", "r") as f:
    pars = json.load(f)
first_dat_file = int(pars["first_dat_file"])
last_dat_file = int(pars["last_dat_file"])
n_exp = pars["n_exp"]
repetitions = pars["repetitions"]
step = pars["step"]
NMR_freq = pars["NMR_freq"]
tumbling_tau_exp = pars["tumbling_tau_exp"]
fit_out = pars["fit_out"]
R1_out = pars["R1_out"]
R2_out = pars["R2_out"]

fun_aligned = Corfun("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file))

# Multiply ACFs for aligned frames by tumbling exponential and save
time_in_steps = fun_aligned.data[0][:, 0]
exp_tau_exp = np.exp(-time_in_steps/(tumbling_tau_exp / step))
for i, acf in enumerate(fun_aligned.data):
    fun_aligned.data[i][:, 1] = np.multiply(acf[:, 1], exp_tau_exp)
tumbling_on_dir = "cor_NH_%d-%d_tumbling_on" % (first_dat_file, last_dat_file)
if not os.path.exists(tumbling_on_dir):
    os.makedirs(tumbling_on_dir)
fun_aligned.write(tumbling_on_dir)

# Process ACF
fun = fun_aligned

dn_mask = 0
up_mask = 1e15 #1e-9 * (last_dat_file - first_dat_file + 1) / step
lb = [0.0]*n_exp + [dn_mask]*n_exp
ub = [1.0]*n_exp + [up_mask]*n_exp

# Fit autocorrelation functions for aligned trajectory
sp = np.zeros(n_exp*2)
sp[:n_exp] = np.array([1.0/n_exp]*n_exp)
sp[n_exp:] = np.logspace(np.log10(step), np.log10(tumbling_tau_exp), n_exp)/step
fun.fit(n_exp=n_exp, rep=repetitions, step=step, sp=sp, lb=lb, ub=ub, fout=fit_out)

# Calculate relaxation rates
fun.calc_R1(opfreq=NMR_freq, X='N15', fnout=R1_out)
fun.calc_R2(opfreq=NMR_freq, X='N15', fnout=R2_out)
