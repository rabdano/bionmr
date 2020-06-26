from bionmr_utils.md import *
from get_secondary_structure_residues import *
from tqdm import tqdm
import numpy as np
import json


#
#  setup trajectory parameters
#
with open("input.json", "r") as f:
    pars = json.load(f)
path_to_traj = pars["path_to_traj"]
first_dat_file = int(pars["first_dat_file"])
last_dat_file = int(pars["last_dat_file"])
stride = int(pars["stride"])
reference_pdb = pars["reference_pdb"]
rmsf_fnout = pars["rmsf_fnout"]


#
#  set atom selections
#
protein_and_dna_chains = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
probe = ((cName == "I") & (aName == "P"))
residues_of_interest = set(list(range(136, 160)) + list(range(623, 647)))
chain_letters = 'BF'


#
#  load trajectory
#
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
#  reference atoms for alignment
#  use 3LZ0 structure with reconstructed tails and propka protonation
#
reference = PdbFile(reference_pdb).get_frame()
ref_ats = reference.asAtoms
ref_align_ca_ss = ref_ats.filter((aName == "CA") & (rId.is_in(ss_residues)))
ca_ref_by_chain = [ref_ats.filter((aName == 'CA') & (cName == ChainName(c)) & (rId.is_in(residues_of_interest))) for c in chain_letters]


#
#  initialize average coordinates with (0,0,0)
#
avg_coords = [ca_ref.toCoords.transform(UniformScale3d(0)) for ca_ref in ca_ref_by_chain]


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
#  run through trajectory to calculate avg. coordinates
#
for frame in tqdm(traj[::stride], disable=False):
    if frame.index == 0:
        # create variables for selections from atoms of frame
        ref_pbc_ats = ref.asAtoms
        frame_ats = frame.asAtoms

        ref_probe = ref_pbc_ats.filter(probe)
        frame_probe = frame_ats.filter(probe)

        frame_align_ca_ss = frame_ats.filter((aName == "CA") & (rId.is_in(ss_residues)))

        ats_ref = []
        ats_frame = []
        for cid in sorted(protein_and_dna_chains):
            ats_ref.append(ref_pbc_ats.filter(cName == cid))
            ats_frame.append(frame_ats.filter(cName == cid))

        ca_frame_by_chain = [
            frame_ats.filter((aName == 'CA') & (cName == ChainName(c)) & (rId.is_in(residues_of_interest)))
            for c in chain_letters]

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

    # align nucleosome in frame by reference 3LZ0 structure using sec.str. CA atoms
    frame_ats.transform(frame_align_ca_ss.alignment_to(ref_align_ca_ss))

    # calculate average coordinates

    for i, ca_frame in enumerate(ca_frame_by_chain):
        # sum all coordinates
        for j, a in enumerate(ca_frame):
            avg_coords[i][j] += a.r

# 'divide' by number of frames
for crds in avg_coords:
    crds.transform(UniformScale3d(1.0 / traj.size * stride))


#
#  align to average coordinates and calculate RMSF for each CA
#
print('align to average coordinates and calculate RMSF for each CA (per residue)')
rmsf = [np.zeros((len(x),)) for x in ca_ref_by_chain]
for frame in tqdm(traj[::stride], disable=False):
    if frame.index == 0:
        # create variables for selections from atoms of frame
        ref_pbc_ats = ref.asAtoms
        frame_ats = frame.asAtoms

        ref_probe = ref_pbc_ats.filter(probe)
        frame_probe = frame_ats.filter(probe)

        frame_align_ca_ss = frame_ats.filter((aName == "CA") & (rId.is_in(ss_residues)))

        ats_ref = []
        ats_frame = []
        for cid in sorted(protein_and_dna_chains):
            ats_ref.append(ref_pbc_ats.filter(cName == cid))
            ats_frame.append(frame_ats.filter(cName == cid))

        ca_frame_by_chain = [
            frame_ats.filter((aName == 'CA') & (cName == ChainName(c)) & (rId.is_in(residues_of_interest)))
            for c in chain_letters]

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

    # align nucleosome in frame by reference 3LZ0 structure using sec.str. CA atoms
    frame_ats.transform(frame_align_ca_ss.alignment_to(ref_align_ca_ss))

    for i, ca_frame in enumerate(ca_frame_by_chain):
        # calculate RMSF
        for j, a in enumerate(ca_frame):
            rmsf[i][j] += (a.r - avg_coords[i][j]).len2()

rmsf = [np.sqrt(rmsf_i / traj.size * stride) for rmsf_i in rmsf]


# write RMSF to file
for i, rmsf_i in enumerate(rmsf):
    with open('RMSF_chain_' + chain_letters[i] + '.csv', 'w') as f:
        f.write('rId,rName,RMSF\n')
        for j, at in enumerate(ca_ref_by_chain[i]):
            f.write('%d,%s,%.8f\n' % (at.rId.serial, at.rName.str, rmsf_i[j]))
