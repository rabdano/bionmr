from bionmr_utils.md import *
from get_secondary_structure_residues import *
from tqdm import tqdm
import numpy as np
from scipy import spatial
import json

#
#  setup parameters
#
with open("input.json", "r") as f:
    pars = json.load(f)
first_dat_file = int(pars["first_dat_file"])
last_dat_file = int(pars["last_dat_file"])
stride = int(pars["stride"])
cutoff = float(pars["cutoff"])
reference_pdb = pars["reference_pdb"]
rmsd_fnout = pars["rmsd_fnout"]
rmsd_contacted_fnout = pars["rmsd_contacted_fnout"]
contacts_1_fnout = pars["contacts_1_fnout"]
contacts_2_fnout = pars["contacts_2_fnout"]
n_steps = last_dat_file - first_dat_file + 1


#
#  set atom selections
#
protein_and_dna_chains = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
residues_of_interest_1 = set(range(136, 160))
residues_of_interest_2 = set(range(623, 647))
probe = ((cName == "I") & (aName == "P"))
dna_align_pred = aName.is_in({"N1", "N9"})
heavy = lambda a: (a.name.str[0] != 'H')
h4_1_pred = rId.is_in(residues_of_interest_1)
h4_2_pred = rId.is_in(residues_of_interest_2)


#
#  load trajectory
#
path_to_traj = "../.."
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/1_build/ref.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file)

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" %
      (len(traj), len(traj[0].asChains), len(traj[0].asResidues), len(traj[0].asAtoms)))
print("Using run%05d.dat - run%05d.dat" % (first_dat_file, last_dat_file))


#
#  get secondary structure information:
#  get residue ids for secondary structure from DSSP using Biopython and PDB
#
protein_chains = "ABCDEFGH"
ss_residues_by_chain = []
incr = 0
for c in protein_chains:
    ss = get_secondary_structure_residues(c, pdb_code="1KX5")
    ss_residues_by_chain.append(set([(x+incr) for x in ss]))
    incr += len(ref.asResidues.filter(cName == c))
# print("Secondary structure residues by chain")
# print(ss_residues_by_chain)
ss_residues = [r for chain_r in ss_residues_by_chain for r in chain_r]
# print("Secondary structure residues combined")
# print(ss_residues)
ss_residues = set(ss_residues)


#
#  get residue ids for DNA pairs
#
dna = ref.asAtoms.filter(aName == "N1")
dna_resids = []
for at in dna:
    dna_resids.append(int(at.rId.serial))
n_dna_resids = len(dna_resids)
n_dna_bp = int(n_dna_resids / 2)
idx_01 = int(dna_resids[0] + n_dna_bp // 2)
idx_02 = int(dna_resids[-1] - n_dna_bp // 2)
dna_bp = [(r1, r2) for r1, r2 in zip(range(idx_01 - int(n_dna_bp//2), idx_01 + int(n_dna_bp//2) + 1, 1),
                                     range(idx_02 + int(n_dna_bp//2), idx_02 - int(n_dna_bp//2) - 1, -1))]
print("Nucleotides:", n_dna_resids)
print("Base pairs:", n_dna_bp)
print("Middle index I:", idx_01)
print("Middle index J:", idx_02)
print("Base pairs:\n", [(r1, r2) for r1, r2 in [dna_bp[0], dna_bp[-1]]])


#
#  reference atoms for alignment
#  use 3LZ0 structure with reconstructed tails and propka protonation
#
reference = PdbFile(reference_pdb).get_frame()
ref_ats = reference.asAtoms
ref_align_ca_ss = ref_ats.filter((aName == "CA") & (rId.is_in(ss_residues)))
ref_align_dna = ref_ats.filter(dna_align_pred)
ref_align_bp = [ref_ats.filter(dna_align_pred & rId.is_in({r1, r2})) for r1, r2 in dna_bp]


#
#  create data variables
#
bp_coordinates = np.zeros((n_dna_bp, int(traj.size/stride), 3), dtype=float)
bp_ref_coordinates = np.zeros((n_dna_bp, 3), dtype=float)
contacts_1 = np.zeros(n_dna_bp, dtype=int)
contacts_2 = np.zeros(n_dna_bp, dtype=int)
rmsd = np.zeros(n_dna_bp, dtype=float)
bp_contacted_1 = np.zeros((n_dna_bp, int(traj.size/stride)), dtype=int)
bp_contacted_2 = np.zeros((n_dna_bp, int(traj.size/stride)), dtype=int)



#
#  Handle PBC:
#
def get_XST(path_to_traj):
    time, V = np.genfromtxt(path_to_traj + "/5_run/summary.VOLUME", unpack=True)

    inpcrd = path_to_traj + "/1_build/box.inpcrd"

    with open(inpcrd) as f:
        content = f.readlines()

    content = content[-1].split()[:6]

    a, b, c, alpha, beta, gamma = [float(x) for x in content]

    # https://en.wikipedia.org/wiki/Parallelepiped
    V0 = a * b * c * np.sqrt(
        1.0 + 2.0 * np.cos(alpha / 180 * np.pi) * np.cos(beta / 180 * np.pi) * np.cos(gamma / 180 * np.pi) - np.cos(
            alpha / 180 * np.pi) ** 2 - np.cos(beta / 180 * np.pi) ** 2 - np.cos(gamma / 180 * np.pi) ** 2)

    # XST file format
    # http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2003-2004/0234.html

    m_v1 = np.array([a, 0.0, 0.0])
    m_v2 = np.array([b * np.cos(gamma / 180 * np.pi), b * np.sin(gamma / 180 * np.pi), 0])
    m_v3 = np.zeros(3)

    m_v3[0] = c * np.cos(beta / 180 * np.pi)
    m_v3[1] = b / m_v2[1] * c * np.cos(alpha / 180 * np.pi) - m_v2[0] / m_v2[1] * m_v3[0]
    m_v3[2] = np.sqrt(c * c - m_v3[0] * m_v3[0] - m_v3[1] * m_v3[1])

    xst = np.zeros((len(V), 13))
    origin = np.zeros(3)

    for i, f_V in enumerate(V):
        scaling_factor = np.cbrt(f_V / V0)
        # print(scaling_factor)
        xst[i, 0] = time[i]
        xst[i, 1:10] = np.array(
            [m_v1[0], m_v1[1], m_v1[2], m_v2[0], m_v2[1], m_v2[2], m_v3[0], m_v3[1], m_v3[2]]) * scaling_factor
        xst[i, 10:13] = [origin[0], origin[1], origin[2]]

    return xst


# get PBC from XST file and inpcrd file
pbc = get_XST(path_to_traj)
v1 = VectorXYZ.from_numpy(pbc[:, 1:4])
v2 = VectorXYZ.from_numpy(pbc[:, 4:7])
v3 = VectorXYZ.from_numpy(pbc[:, 7:10])
assert len(v1) >= traj.size, "PBC data is not enough"
scaling_factors = [v1[i].len()/v1[0].len() for i in range(len(v1))]
lat_vec = LatticeVectors(v1[0], v2[0], v3[0])
shift_finder = BestShiftFinder(lat_vec)
inpcrd = path_to_traj + "/1_build/box.inpcrd"
with open(inpcrd) as f:
    content = f.readlines()
content = content[-1].split()[:6]
a, b, c, alpha, beta, gamma = [float(x) for x in content]

#
#  run through trajectory
#
for frame in tqdm(traj[::stride], disable=False):
    if frame.index == 0:
        # create variables for selections from atoms of frame
        ref_pbc_ats = ref.asAtoms
        frame_ats = frame.asAtoms

        ref_probe = ref_pbc_ats.filter(probe)
        frame_probe = frame_ats.filter(probe)

        frame_align_ca_ss = frame_ats.filter((aName == "CA") & (rId.is_in(ss_residues)))
        frame_align_bp = [frame_ats.filter(dna_align_pred & rId.is_in({r1, r2})) for r1, r2 in dna_bp]
        frame_contact_bp = [frame_ats.filter(rId.is_in({r1, r2})).filter(heavy) for r1, r2 in dna_bp]

        ats_ref = []
        ats_frame = []
        for cid in sorted(protein_and_dna_chains):
            ats_ref.append(ref_pbc_ats.filter(cName == cid))
            ats_frame.append(frame_ats.filter(cName == cid))

        h4_1 = frame.asAtoms.filter(h4_1_pred).filter(heavy)
        h4_2 = frame.asAtoms.filter(h4_2_pred).filter(heavy)

        # align reference by first frame
        # put reference bp coordinates to bp_ref_coordinates
        ref_pbc_ats.transform(ref_probe.alignment_to(frame_probe))
        for i, ref_bp in enumerate(ref_align_bp):
            assert len(ref_bp) == 3, len(ref_bp)
            bp_ref_coordinates[i, :] = calc_geom_center(ref_bp.toCoords).to_np

    # unwrapping
    # align all reference atoms by one of DNA chains
    ref_pbc_ats.transform(ref_probe.alignment_to(frame_probe))

    # cycle through chains and check first atom displacement
    shift_finder.scale_lattice_by(scaling_factors[int(frame.index / stride)])
    for i, cid in enumerate(sorted(protein_and_dna_chains)):
        dif, shift = shift_finder.find_best_shift(ats_ref[i].geom_center(),
                                                  ats_frame[i].geom_center())
        ats_frame[i].transform(Translation3d(shift))
    shift_finder.scale_lattice_by(1.0 / scaling_factors[int(frame.index / stride)])

    # calculate contacts
    # create trees
    h4_1_coords = h4_1.toCoords.to_numpy()
    h4_2_coords = h4_2.toCoords.to_numpy()
    bp_trees = []
    for bp_ats in frame_contact_bp:
        bp_trees.append(spatial.cKDTree(bp_ats.toCoords.to_numpy()))

    # count contacts and save indexes for RMSD calculation (pick bp)
    for i, bp_tree in enumerate(bp_trees):
        # H4-1
        dist, _ = bp_tree.query(h4_1_coords, k=1, eps=1e-5, n_jobs=-1, distance_upper_bound=cutoff)
        assert len(dist.shape) == 1
        if np.sum(np.isfinite(dist)) > 0:
            contacts_1[i] += 1
            bp_contacted_1[i, int(frame.index / stride)] = 1
        # H4-2
        dist, _ = bp_tree.query(h4_2_coords, k=1, eps=1e-5, n_jobs=-1, distance_upper_bound=cutoff)
        assert len(dist.shape) == 1
        if np.sum(np.isfinite(dist)) > 0:
            contacts_2[i] += 1
            bp_contacted_2[i, int(frame.index / stride)] = 1

    # align nucleosome in frame by reference 3LZ0 structure using sec.str. CA atoms
    frame_ats.transform(frame_align_ca_ss.alignment_to(ref_align_ca_ss))

    # calculate bp center of mass coordinates
    for i, frame_bp in enumerate(frame_align_bp):
        assert len(frame_bp) == 3, len(frame_bp)
        bp_coordinates[i, int(frame.index / stride), :] = calc_geom_center(frame_bp.toCoords).to_np


#
#  calculate rmsd
#
cum_sd_contacted_1 = 0
cum_sd_contacted_2 = 0
cum_sd_contacted = 0
n_contacted_1 = 0
n_contacted_2 = 0
n_contacted = 0
for i in range(n_dna_bp):
    cum_sd = 0
    for j in range(int(traj.size/stride)):
        sd = (np.linalg.norm(bp_coordinates[i, j, :] - bp_ref_coordinates[i, :]))**2
        assert np.sqrt(sd) < (a/2), ("sd =", np.sqrt(sd), "- translation error (?)", "frame.index =", j, "bp =", i)
        cum_sd += sd
        if bp_contacted_1[i, j] != 0:
            cum_sd_contacted_1 += sd
            n_contacted_1 += 1
        if bp_contacted_2[i, j] != 0:
            cum_sd_contacted_2 += sd
            n_contacted_2 += 1
        if (bp_contacted_1[i, j] != 0) or (bp_contacted_2[i, j] != 0):
            cum_sd_contacted += sd
            n_contacted += 1
    rmsd[i] = np.sqrt(cum_sd / int(traj.size/stride))


#
#  write RMSD contacted to file
#
with open(rmsd_contacted_fnout, 'w') as f:
    f.write("Mean RMSD for bp contacted by H4-1: %.5f\n" % np.sqrt(cum_sd_contacted_1 / n_contacted_1))
    f.write("Mean RMSD for bp contacted by H4-2: %.5f\n" % np.sqrt(cum_sd_contacted_2 / n_contacted_2))
    f.write("Mean RMSD for bp contacted by H4-1 or H4-2: %.5f\n" % np.sqrt(cum_sd_contacted / n_contacted))


#
#  write RMSD to file
#
position = np.arange(-(n_dna_bp//2), n_dna_bp//2 + 1, 1)
np.savetxt(
    fname=rmsd_fnout,
    X=np.vstack((position, rmsd)).T,
    fmt="%.5f",
    header="position, rmsd[A]",
    delimiter=","
)

# write contacts to file
position = np.arange(-(n_dna_bp//2), n_dna_bp//2 + 1, 1)
np.savetxt(
    fname=contacts_1_fnout,
    X=np.vstack((position, contacts_1)).T,
    fmt="%d",
    header="position, contacts",
    delimiter=","
)
np.savetxt(
    fname=contacts_2_fnout,
    X=np.vstack((position, contacts_2)).T,
    fmt="%d",
    header="position, contacts",
    delimiter=","
)
