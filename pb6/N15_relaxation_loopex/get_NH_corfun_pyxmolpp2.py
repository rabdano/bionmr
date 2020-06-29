import numpy as np
from math import ceil
from bionmr_utils.md import *
import os
import sys
from correlation_functions import cor
from get_secondary_structure_residues import *
from tqdm import tqdm
import json


#
#  setup trajectory parameters
#
with open("input.json", "r") as f:
    pars = json.load(f)
path_to_traj = pars["path_to_traj"]
first_dat_file = int(pars["first_dat_file"])
last_dat_file = int(pars["last_dat_file"])
n_steps = last_dat_file - first_dat_file + 1
stride = int(pars["stride"])
cut_autocorr_function = n_steps
fft_acf = pars["fft_acf"]
scaling = float(pars["scaling"])
ratio = float(pars["ratio"])
remove_first_point = pars["remove_first_point"]
align = pars["align"]
ig_like_domain_start=pars["ig_like_domain_start"]


#
#  read trajectory
#
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/0_prepare/ref.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file,
                          subdir="6_run",
                          filetype="nc")

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" % (len(traj),
                                                                                len(traj[0].asChains),
                                                                                len(traj[0].asResidues),
                                                                                len(traj[0].asAtoms)))
print("Using run%05d.nc - run%05d.nc" % (first_dat_file, last_dat_file))


#
#  get secondary structure information:
#  get residue ids for secondary structure from DSSP using Biopython and PDB
#
protein_chains = "A"
ss_residues_by_chain = []
incr = 0
for c in protein_chains:
    ss = get_secondary_structure_residues(c, pdb_code="5NGJ")
    ss_residues_by_chain.append(set([(x+incr) for x in ss]))
    incr += len(ref.asResidues.filter(cName == c))
# print("Secondary structure residues by chain")
# print(ss_residues_by_chain)
ss_residues = []
for chain_r in ss_residues_by_chain:
    for r in chain_r:
        if r < ig_like_domain_start:
            ss_residues.append(r)
print("Secondary structure residues combined")
print(ss_residues)
ss_residues = set(ss_residues)


#
#  set atom selections & predicates
#
align_pred = (aName == "CA") & (rId.is_in(ss_residues))
residues_of_interest = set(list(range(1, 465)))


#
#  create folders and open files
#
resids = []
resnames = []
for at in ref.asAtoms.filter((aName == "H") & (rId.is_in(residues_of_interest))):
    resids.append(at.rId.serial)
    resnames.append(at.rName.str)
resi = tuple(zip(resids, resnames))
print("Autocorrelation functions will be calculated for following residues:")
print(resi)

if not os.path.exists("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file)):
    os.makedirs("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file))
# vectors - list of VectXYZ
vectors = np.zeros((int(traj.size / stride), len(resi), 3))
skip_frames = set()


#
#  run through trajectory and calculate vectors
#
for frame in tqdm(traj[::stride], disable=False):
    if frame.index == 0:
        frame_ats = frame.asAtoms
        # align frame by ref
        if align:
            ref_ats = ref.asAtoms.filter(align_pred)
            align_ats = frame.asAtoms.filter(align_pred)
            frame_ats.transform(align_ats.alignment_to(ref_ats))
        # NH vectors
        N, H = [], []
        for resid in resids:
            N.append(frame.asAtoms.filter((rId == resid) & (aName == "N"))[0])
            H.append(frame.asAtoms.filter((rId == resid) & (aName == "H"))[0])

    # align by sec.str. Ca for tube domain
    if align:
        frame_ats.transform(align_ats.alignment_to(ref_ats))

    # calculate vectors
    for i, N_at, H_at in zip(range(len(resids)), N, H):
        vec = H_at.r - N_at.r
        if vec.len() < 0.5:
            skip_frames.add(frame.index)
        vectors[int(frame.index / stride), i] = vec.to_np


if len(skip_frames) > 0:
    fft_acf = False
    np.savetxt('skip_frames.txt', np.array(list(skip_frames)), delimiter=',', header='')


#
#  calculate autocorrelation functions
#
steps = np.arange(int(cut_autocorr_function * 1000 / stride))
grid = []
nlim = ratio * len(vectors[:, 0, 0])
if not remove_first_point:
    grid.append(0)
tau = 1.0
while tau <= nlim:
    grid.append(int(tau))
    tau = ceil(tau * scaling)

print("Calculating autocorrelation functions...")
for i, rid, rname in zip(range(len(resids)), resids, resnames):
    sys.stdout.write("Residue %d of %d\r" % (i+1, len(resids)))
    if fft_acf:
        ac = np.array(calc_autocorr_order_2(VectorXYZ.from_numpy(vectors[:, i, :]),
                                            limit=int(cut_autocorr_function*1000/stride)))
        np.savetxt("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file)+"/%04d_%s.cor" % (rid, rname),
                   np.vstack((steps[grid], ac[grid])).T, fmt="%14.6e")
    else:
        ac = cor(vectors[i], grid, skip_frames)
        np.savetxt("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file) + "/%04d_%s.cor" % (rid, rname),
                   np.vstack((steps[grid], ac)).T, fmt="%14.6e")


sys.stdout.write("\n")
print("Done!")
