import numpy as np
from math import ceil
from bionmr_utils.md import *
import os
import sys
from correlation_functions import cor


# setup parameters.
path_to_traj = "../.."
first_dat_file = 1
last_dat_file = 1000
n_steps = last_dat_file - first_dat_file + 1
stride = 1
cut_autocorr_function = n_steps  # ns
fft_acf = True
scaling = 1.005
ratio = 0.8
remove_first_point = True
align = True
align_pred = (aName == "CA")

residues_of_interest = set(list(range(1, 90)))

# read trajectory
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/5_run/run00001.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file)

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" % (len(traj),
                                                                                len(traj[0].asChains),
                                                                                len(traj[0].asResidues),
                                                                                len(traj[0].asAtoms)))
print("Using run%05d.dat - run%05d.dat" % (first_dat_file, last_dat_file))

# create folders and open files
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
vectors = []
skip_frames = set()

first_time = True

# run through trajectory and calculate vectors
print("Processing frames...")
for frame in traj[::stride]:
    sys.stdout.write("Frame %d of %d\r" % (frame.index+1, traj.size/stride))
    
    if first_time:
        # align frame by ref
        if align:
            ref_ats = ref.asAtoms.filter(align_pred)
            align_ats = frame.asAtoms.filter(align_pred)
        alignment = calc_alignment(ref_ats.toCoords, align_ats.toCoords)
        frame.asAtoms.transform(alignment)
        # NH vectors
        N, H = [], []
        for resid in resids:
            N.append(frame.asAtoms.filter((rId == resid) & (aName == "N"))[0])
            H.append(frame.asAtoms.filter((rId == resid) & (aName == "H"))[0])
            vectors.append(VectorXYZ([(H[-1].r - N[-1].r)]))
        first_time = False

    else:
        # align frame by ref
        alignment = calc_alignment(ref_ats.toCoords, align_ats.toCoords)
        frame.asAtoms.transform(alignment)
        for i, N_at, H_at in zip(range(len(resids)), N, H):
            vec = H_at.r - N_at.r
            if vec.len() < 0.5:
                skip_frames.add(frame.index)
            vectors[i].append(H_at.r - N_at.r)
sys.stdout.write("\n")


if len(skip_frames) > 0:
    fft_acf = False
    np.savetxt('skip_frames.txt', np.array(list(skip_frames)), delimiter=',', header='')


# calculate autocorrelation functions
steps = np.arange(int(cut_autocorr_function*1000/stride))
grid = []
nlim = ratio * len(vectors[0])
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
        ac = np.array(calc_autocorr_order_2(vectors[i], limit=int(cut_autocorr_function*1000/stride)))
        np.savetxt("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file)+"/%04d_%s.cor" % (rid, rname),
                   np.vstack((steps[grid], ac[grid])).T, fmt="%14.6e")
    else:
        ac = cor(vectors[i], grid, skip_frames)
        np.savetxt("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file) + "/%04d_%s.cor" % (rid, rname),
                   np.vstack((steps[grid], ac)).T, fmt="%14.6e")


sys.stdout.write("\n")
print("Done!")
