import numpy as np
from math import ceil
from bionmr_utils.md import *
import os
import sys
from correlation_functions import cor, gets2
from glob import glob


#
#  Setup parameters
#
path_to_traj = "../.."
first_dat_file = 0
last_dat_file = 1000
n_chunk = 16  # 29 residues per ring
n_steps = last_dat_file - first_dat_file + 1
stride = 1
cut_autocorr_function = n_steps  # ns
fft_acf = True
scaling = 1.005
ratio = 0.8
remove_first_point = True
align = True
align_pred = (aName == "CA")
n_res = 464
n_chains = 9


#
#  Residues ranges for alignment: Secondary structure ranges
#
ssca = [(93, 99), (186, 189), (193, 200), (271, 274), (279, 287), (359, 368), (384, 386), (9, 14), (25, 28),
        (34, 38), (68, 75), (79, 81), (86, 88), (112, 114), (119, 121), (135, 140), (143, 148), (154, 158),
        (166, 173), (175, 179), (209, 214), (221, 234), (260, 267), (294, 299), (309, 314), (317, 319),
        (333, 338), (350, 354), (378, 383), (389, 391), (395, 403), (412, 418), (422, 424), (429, 432),
        (438, 446), (453, 461)]
ssca_no_ig_like = [(93, 99), (186, 189), (193, 200), (271, 274), (279, 287), (359, 368), (384, 386), (9, 14), (25, 28),
                   (34, 38), (68, 75), (79, 81), (86, 88), (112, 114), (119, 121), (135, 140), (143, 148), (154, 158),
                   (166, 173), (175, 179), (209, 214), (221, 234), (260, 267), (294, 299), (309, 314), (317, 319),
                   (333, 338), (350, 354)]
all_rings_ssca_ranges = [(a + i * 464, b + i * 464) for (a, b) in ssca for i in [0, 1, 2, 3, 4, 5, 6, 7, 8]]
all_rings_ssca_no_ig_like_ranges = [(a + i * 464, b + i * 464) for (a, b) in ssca_no_ig_like for i in [0, 1, 2, 3, 4, 5, 6, 7, 8]]

# pick only CA of tube excluding ig-like domain
all_rings_ssca = set()
for ran in all_rings_ssca_no_ig_like_ranges:
    all_rings_ssca = all_rings_ssca.union(set([i for i in range(ran[0], ran[1]+1)]))

align_pred = (rId.is_in(all_rings_ssca) & (aName == "CA"))
# residues_of_interest = set(range(1, 464*9+1))
# residues_of_interest = set(list(range(1, 20)))
# residues_of_interest = [[c for c in range(a + i * 464, b + i * 464 + 1)] for (a, b) in [(1, 5)] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8]]
# residues_of_interest = set([item for sublist in residues_of_interest for item in sublist])
residues_of_interest = []
for i in range(n_chunk):
    chunk = []
    for j in range(n_chains):
        idx1 = 464 * j + int(n_res/n_chunk) * i + 1
        idx2 = 464 * j + int(n_res/n_chunk) * (i + 1)
        chunk += list(range(idx1, idx2 + 1, 1))
    residues_of_interest.append(chunk)


#
# Read trajectory
#
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/0_prepare/protein_with_chains.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file,
                          subdir="8_unwrap",
                          filetype="nc")

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" % (len(traj),
                                                                                len(traj[0].asChains),
                                                                                len(traj[0].asResidues),
                                                                                len(traj[0].asAtoms)))
print("Using run%05d.nc - run%05d.nc" % (first_dat_file, last_dat_file))


#
#  Initialize folders
#
if not os.path.exists("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file)):
    os.makedirs("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file))
if not os.path.exists("cor_NH_%d-%d_averaged" % (first_dat_file, last_dat_file)):
    os.makedirs("cor_NH_%d-%d_averaged" % (first_dat_file, last_dat_file))



#
#  Process trajectory by chunks
#
skip_frames = []
s2s = np.zeros((n_res, 1))

for c, chunk in enumerate(residues_of_interest):
    print("\n" + ("=" * 80))
    print("Processing frames for chunk %d/%d..." % (c + 1, n_chunk))

    hydrogen_pred = ((aName == "H") & (rId.is_in(set(chunk))))
    resids = []
    resnames = []
    vectors = []
    chain_names = [chain.cName.str for chain in ref.asChains]
    for cn in chain_names:
        resids_ = []
        resnames_ = []
        for at in ref.asAtoms.filter(cName == cn).filter(hydrogen_pred):
            resids_.append(at.rId.serial)
            resnames_.append(at.rName.str)
        resids.append(resids_)
        resnames.append(resnames_)
        # vectors are shaped [n_chains, n_resids_in_chain, n_frames, xyz]
        vectors.append(np.zeros((len(ref.asAtoms.filter(cName == cn).filter(hydrogen_pred)), int(traj.size / stride), 3), dtype="float"))

    #
    #  Run through trajectory and calculate vectors
    #
    s = "Autocorrelation functions will be calculated for following residues:\n"
    for chain_resnames, chain_resids in zip(resnames, resids):
        for resname, resid in zip(chain_resnames, chain_resids):
            s += "{:s}-{:d} ".format(resname, resid)
    print(s)

    for frame in traj[::stride]:
        sys.stdout.write("Frame %d of %d\r" % (int(frame.index / stride) + 1, int(traj.size / stride)))

        # create selections on first frame
        if frame.index == 0:
            ref_align_ats = ref.asAtoms.filter(align_pred)
            frame_align_ats = frame.asAtoms.filter(align_pred)
            frame_ats = frame.asAtoms.filter(aName.is_in({"H", "N"}))

            # NH atom selections
            set_resids = set([item for resids_ in resids for item in resids_])
            Ns = frame_ats.filter(rId.is_in(set_resids)).filter(aName == "N")
            Hs = frame_ats.filter(rId.is_in(set_resids)).filter(aName == "H")

        frame_ats.transform(frame_align_ats.alignment_to(ref_align_ats))

        # calculate vectors
        for N, H in zip(Ns, Hs):
            assert N.rId == H.rId, "N-H mismatch"
            assert H.cName.str in chain_names, "H not in chain"
            j = chain_names.index(H.cName.str)
            assert H.rId.serial in resids[j], "H resid not in list of resids"
            i = resids[j].index(H.rId.serial)
            vec = H.r - N.r
            if vec.len() > 0.5:
                vectors[j][i, int(frame.index / stride), :] = vec.to_np
            else:
                skip_frames.append([frame.index, vec.len()])

    sys.stdout.write("\n")



    #
    #  Check vectors
    #
    print("Checking extracted vectors...")
    flag = True
    for i, vec in enumerate(vectors):
        for j, res_vec in enumerate(vec):
            for k, r in enumerate(res_vec):
                if np.isnan(r[0]):
                    print("NAN!", i, j, k)
                    flag = False
                elif r[0] == 0.0:
                    print("ZERO!", i, j, k)
                    flag = False
    if flag:
        print("Vectors are fine.")
    else:
        print("Some problems found.")

    #
    #  Save indexes of corrupted frames
    #
    print("Checking trajectory for corrupted frames...")
    flag = True
    if len(skip_frames) > 0:
        fft_acf = False
        flag = False
        print("WARNING: Some corrupted frames found!")
        np.savetxt('skip_frames.txt', np.array(skip_frames), delimiter=',', fmt=("%d", "%.3f"))
    if flag:
        print("Frames are fine.")
    else:
        print("Some problems found.")


    #
    #  Calculate autocorrelation functions
    #
    print("Calculating autocorrelation functions...")
    steps = np.arange(int(cut_autocorr_function * 1000 / stride))
    grid = []
    nlim = ratio * len(vectors[0][0])
    if not remove_first_point:
        grid.append(0)
    tau = 1.0
    while tau <= nlim:
        grid.append(int(tau))
        tau = ceil(tau * scaling)

    # loop through chains
    for j, chain_resnames, chain_resids in zip(range(len(vectors)), resnames, resids):
        # loop through residues in chain
        for i, rid, rname in zip(range(len(chain_resids)), chain_resids, chain_resnames):
            sys.stdout.write("{:s}-{:d} ".format(rname, rid))
            if fft_acf:
                ac = np.array(calc_autocorr_order_2(VectorXYZ.from_numpy(vectors[j][i, :, :]), limit=int(cut_autocorr_function*1000/stride)))
                np.savetxt("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file)+"/%04d_%s.cor" % (rid, rname),
                           np.vstack((steps[grid], ac[grid])).T, fmt="%14.6e")
            else:
                ac = cor(VectorXYZ.from_numpy(vectors[j][i, :, :]), grid, skip_frames)
                np.savetxt("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file) + "/%04d_%s.cor" % (rid, rname),
                           np.vstack((steps[grid], ac)).T, fmt="%14.6e")

    sys.stdout.write("\n")


    #
    #  Calculate order parameters
    #
    print("Calculating order parameters...")
    if not os.path.exists("vec_traces_for_s2"):
        os.makedirs("vec_traces_for_s2")

    n_frames = int(traj.size / stride)
    # loop through residues of first chain
    for i, resid, resname in zip(range(len(resids[0])), resids[0], resnames[0]):
        sys.stdout.write("%s-%d " % (resname, resid))

        # concatenate NH vectors for 9 chains
        long_vec = np.zeros((n_frames * n_chains, 3))
        for j in range(len(resids)):
            long_vec[n_frames * j:n_frames * j + n_frames, :] = vectors[j][i, :, :]

        # calculate order parameter from concatenated NH vectors
        s2 = gets2(long_vec)
        s2s[resid - 1] = s2

        # save vectors for s2
        np.savetxt("vec_traces_for_s2/%03d-%s.txt" % (resid, resname), long_vec, fmt="%7.3f")

    sys.stdout.write("\n")

#
#  Write order parameters to file
#
one_chain_resids = range(1, n_res + 1)
one_chain_resnames = [r.rName.str for r in ref.asChains[0].asResidues]
with open("s2.txt", "w") as f:
    for resid, resname, s2 in zip(one_chain_resids, one_chain_resnames, s2s):
        f.write("%04d-%3s %.3f\n" % (resid, resname, s2))


#
#  Average autocorrelation functions
#
print("Averaging autocorrelation functions...")

acf_files = sorted(glob("cor_NH_%d-%d_diluted" % (first_dat_file, last_dat_file)+"/*.cor"))

resids = []
resnames = []

for acf_file in acf_files:
    bn = os.path.basename(acf_file)
    name = bn.rstrip(".cor")
    resids.append(int(name.split("_")[0]))
    resnames.append(name.split("_")[1])

for i in range(n_res):
    fs = []
    for acf_file, res in zip(acf_files, resids):
        if (((res - 1) % n_res) + 1) == i:
            fs.append(acf_file)
    if len(fs) == 9:
        acf_cum = np.genfromtxt(fs[0], usecols=(0, 1))
        for k in range(1, 9, 1):
            acf_cum += np.genfromtxt(fs[k], usecols=(0, 1))
        np.savetxt("cor_NH_%d-%d_averaged" % (first_dat_file, last_dat_file) + "/" + os.path.basename(fs[0]),
                   acf_cum / 9,
                   fmt="%14.6e")
    else:
        print(fs)
        continue


print("Done!")
