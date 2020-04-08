import numpy as np
from scipy import spatial
from bionmr_utils.md import *
from tqdm import tqdm
import json


# setup parameters
with open("input.json", "r") as f:
    pars = json.load(f)
first_dat_file = int(pars["first_dat_file"])
last_dat_file = int(pars["last_dat_file"])
stride = int(pars["stride"])
cutoff = float(pars["cutoff"])
contacts_fnout = pars["contacts_fnout"]
n_steps = last_dat_file - first_dat_file + 1

protein_and_dna_chains = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
residues_of_interest_1 = set(range(136, 160))
residues_of_interest_2 = set(range(623, 647))


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


# load trajectory
path_to_traj = "../.."
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/1_build/ref.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file)

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" %
      (len(traj), len(traj[0].asChains), len(traj[0].asResidues), len(traj[0].asAtoms)))
print("Using run%05d.dat - run%05d.dat" % (first_dat_file, last_dat_file))


# get residue ids for DNA pairs
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


# create cumulative variables
contacts = np.zeros(n_dna_bp)


# create predicates
heavy = lambda a: not (a.name.str[0] == 'H')
h4_1_pred = rId.is_in(residues_of_interest_1)
h4_2_pred = rId.is_in(residues_of_interest_2)
probe = ((cName == "I") & (aName == "P"))


# get PBC from XST file and inpcrd file
pbc = get_XST(path_to_traj)
v1 = VectorXYZ.from_numpy(pbc[:, 1:4])
v2 = VectorXYZ.from_numpy(pbc[:, 4:7])
v3 = VectorXYZ.from_numpy(pbc[:, 7:10])
scaling_factors = [v1[i].len()/v1[0].len() for i in range(len(v1))]
lat_vec = LatticeVectors(v1[0], v2[0], v3[0])
shift_finder = BestShiftFinder(lat_vec)
inpcrd = path_to_traj + "/1_build/box.inpcrd"
with open(inpcrd) as f:
    content = f.readlines()
content = content[-1].split()[:6]
a, b, c, alpha, beta, gamma = [float(x) for x in content]


# run through trajectory and calculate number of contacts
print("Processing frames...")
for frame in tqdm(traj[::stride]):
    if frame.index == 0:
        ref_ats = ref.asAtoms
        frame_ats = frame.asAtoms

        # unwrpapping
        ref_probe = ref_ats.filter(probe)
        frame_probe = frame_ats.filter(probe)
        alignment = calc_alignment(frame_ats.filter(probe).toCoords, ref_ats.filter(probe).toCoords)
        ref_ats.transform(alignment)

        ats_ref = []
        ats_frame = []
        for cid in sorted(protein_and_dna_chains):
            ats_ref.append(ref_ats.filter(cName == cid))
            ats_frame.append(frame_ats.filter(cName == cid))

        # frame ats for calculating contacts
        h4 = frame.asAtoms.filter(h4_1_pred | h4_2_pred).filter(heavy)
        frame_bp = [frame_ats.filter(rId.is_in({r1, r2})).filter(heavy) for r1, r2 in dna_bp]

    # unwrapping
    # align all reference atoms by one of DNA chains
    ref_ats.transform(ref_probe.alignment_to(frame_probe))

    # cycle through chains and check first atom displacement
    shift_finder.scale_lattice_by(scaling_factors[frame.index])
    for i, cid in enumerate(sorted(protein_and_dna_chains)):
        dif, shift = shift_finder.find_best_shift(ats_ref[i].geom_center(),
                                                  ats_frame[i].geom_center())
        ats_frame[i].transform(Translation3d(shift))
    shift_finder.scale_lattice_by(1.0 / scaling_factors[frame.index])

    # calculate contacts
    # create trees
    h4_coords = h4.toCoords.to_numpy()
    bp_trees = []
    for bp_ats in frame_bp:
        bp_trees.append(spatial.cKDTree(bp_ats.toCoords.to_numpy()))

    # count contacts
    for i, bp_tree in enumerate(bp_trees):
        dist, indexes = bp_tree.query(h4_coords, k=1, eps=1e-5, n_jobs=-1, distance_upper_bound=cutoff)
        assert len(dist.shape) == 1
        if np.sum(np.isfinite(dist)) > 0:
            contacts[i] += 1

# write contacts to file
position = np.arange(-(n_dna_bp//2), n_dna_bp//2 + 1, 1)
np.savetxt(
    fname=contacts_fnout,
    X=np.vstack((position, contacts)).T,
    fmt="%d",
    header="position, contacts",
    delimiter=","
)
