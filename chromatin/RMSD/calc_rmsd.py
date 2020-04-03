from bionmr_utils.md import *
from get_secondary_structure_residues import *
from tqdm import tqdm
import numpy as np
from timeit import default_timer as timer


# setup trajectory parameters
first_dat_file = 51
last_dat_file = 1000
stride = 1000  # ps
# reference PDB
reference_pdb = "/home/seva/chromatin/5_solution_Widom_601/pdb/Amber/1_propka/01531-propka.pdb"
# probe atoms
protein_and_dna_chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
probe = ((cName == 'I') & (aName == 'P'))
dna_align_pred = ((aName == "N1") | (aName == "N9"))


def get_XST(path_to_traj):
    time, V = np.genfromtxt(path_to_traj + '/5_run/summary.VOLUME', unpack=True)

    inpcrd = path_to_traj + '/1_build/box.inpcrd'

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

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" % (
    len(traj), len(traj[0].asChains), len(traj[0].asResidues), len(traj[0].asAtoms)))
print("Using run%05d.dat - run%05d.dat" % (first_dat_file, last_dat_file))


# get residue ids for secondary structure from DSSP using Biopython and PDB
protein_chains = 'ABCDEFGH'
ss_residues_by_chain = []
incr = 0
for c in protein_chains:
    ss = get_secondary_structure_residues(c, pdb_code='1KX5')
    ss_residues_by_chain.append(set([(x+incr) for x in ss]))
    incr += len(ref.asResidues.filter(cName == c))
# print("Secondary structure residues by chain")
# print(ss_residues_by_chain)
ss_residues = [r for chain_r in ss_residues_by_chain for r in chain_r]
# print("Secondary structure residues combined")
# print(ss_residues)
ss_residues = set(ss_residues)


# reference atoms for alignment
reference = PdbFile(reference_pdb).get_frame()
ref_ats = reference.asAtoms
ref_align_protein = reference.asAtoms.filter(cName.is_in({"A", "B", "C", "D", "E", "F", "G", "H"}))
ref_align_ca = reference.asAtoms.filter(aName == "CA")
ref_align_ca_ss = reference.asAtoms.filter((aName == "CA") & (rId.is_in(ss_residues)))
ref_align_dna = reference.asAtoms.filter(dna_align_pred)

rmsd_ref_align_protein = np.zeros(int(traj.size / stride))
rmsd_ref_align_ca = np.zeros(int(traj.size / stride))
rmsd_ref_align_ca_ss = np.zeros(int(traj.size / stride))
rmsd_ref_align_dna = np.zeros(int(traj.size / stride))


# get PBC from XST file and inpcrd file
pbc = get_XST(path_to_traj)
v1 = VectorXYZ.from_numpy(pbc[:, 1:4])
v2 = VectorXYZ.from_numpy(pbc[:, 4:7])
v3 = VectorXYZ.from_numpy(pbc[:, 7:10])
scaling_factors = [v1[i].len()/v1[0].len() for i in range(len(v1))]
lat_vec = LatticeVectors(v1[0], v2[0], v3[0])
shift_finder = BestShiftFinder(lat_vec)
inpcrd = path_to_traj + '/1_build/box.inpcrd'
with open(inpcrd) as f:
    content = f.readlines()
content = content[-1].split()[:6]
a, b, c, alpha, beta, gamma = [float(x) for x in content]

t2 = timer()
# run through trajectory
for frame in tqdm(traj[::stride], disable=False):
    if frame.index == 0:
        frame_ats = frame.asAtoms

        ref_probe = ref_ats.filter(probe)
        frame_probe = frame_ats.filter(probe)

        frame_align_protein = frame_ats.filter(cName.is_in({"A", "B", "C", "D", "E", "F", "G", "H"}))
        frame_align_ca = frame_ats.filter(aName == "CA")
        frame_align_ca_ss = frame_ats.filter((aName == "CA") & (rId.is_in(ss_residues)))
        frame_align_dna = frame_ats.filter(dna_align_pred)

        # align reference by first frame nucleic P
        alignment = calc_alignment(frame_ats.filter(probe).toCoords, ref_ats.filter(probe).toCoords)
        ref_ats.transform(alignment)

        ats_ref = []
        ats_frame = []
        for cid in sorted(protein_and_dna_chains):
            ats_ref.append(ref_ats.filter(cName == cid))
            ats_frame.append(frame_ats.filter(cName == cid))

    # align all reference atoms by one of DNA chains
    ref_ats.transform(ref_probe.alignment_to(frame_probe))

    # cycle through chains and check first atom displacement
    shift_finder.scale_lattice_by(scaling_factors[frame.index])
    for i, cid in enumerate(sorted(protein_and_dna_chains)):
        dif, shift = shift_finder.find_best_shift(ats_ref[i].geom_center(),
                                                  ats_frame[i].geom_center())
        ats_frame[i].transform(Translation3d(shift))
    shift_finder.scale_lattice_by(1.0/scaling_factors[frame.index])

    # calculate RMSD
    # frame_align_protein.align_to(ref_align_protein)
    # rmsd_ref_align_protein[int(frame.index / stride)] = frame_align_protein.rmsd(ref_align_protein)
    alignment = calc_alignment(ref_align_protein.toCoords, frame_align_protein.toCoords)
    rmsd_ref_align_protein[int(frame.index / stride)] = calc_rmsd(ref_align_protein.toCoords, frame_align_protein.toCoords, alignment)

    # frame_align_ca.align_to(ref_align_ca)
    # rmsd_ref_align_ca[int(frame.index / stride)] = frame_align_ca.rmsd(ref_align_ca)
    alignment = calc_alignment(ref_align_ca.toCoords, frame_align_ca.toCoords)
    rmsd_ref_align_ca[int(frame.index / stride)] = calc_rmsd(ref_align_ca.toCoords, frame_align_ca.toCoords, alignment)

    # frame_align_ca_ss.align_to(ref_align_ca_ss)
    # rmsd_ref_align_ca_ss[int(frame.index / stride)] = frame_align_ca_ss.rmsd(ref_align_ca_ss)
    alignment = calc_alignment(ref_align_ca_ss.toCoords, frame_align_ca_ss.toCoords)
    rmsd_ref_align_ca_ss[int(frame.index / stride)] = calc_rmsd(ref_align_ca_ss.toCoords, frame_align_ca_ss.toCoords, alignment)

    # frame_align_dna.align_to(ref_align_dna)
    # rmsd_ref_align_dna[int(frame.index / stride)] = frame_align_dna.rmsd(ref_align_dna)
    alignment = calc_alignment(ref_align_dna.toCoords, frame_align_dna.toCoords)
    rmsd_ref_align_dna[int(frame.index / stride)] = calc_rmsd(ref_align_dna.toCoords, frame_align_dna.toCoords, alignment)


# write RMSD to file
time = np.linspace((first_dat_file - 1) * 1000 + stride, last_dat_file * 1000, int(traj.size / stride))
np.savetxt(
    fname="rmsd_ref_align_protein.csv",
    X=np.vstack((time, rmsd_ref_align_protein)).T,
    fmt="%.5f",
    header="time[ps], rmsd_ref_align_protein[A]",
    delimiter=","
)
np.savetxt(
    fname="rmsd_ref_align_ca.csv",
    X=np.vstack((time, rmsd_ref_align_ca)).T,
    fmt="%.5f",
    header="time[ps], rmsd_ref_align_ca[A]",
    delimiter=","
)
np.savetxt(
    fname="rmsd_ref_align_ca_ss.csv",
    X=np.vstack((time, rmsd_ref_align_ca_ss)).T,
    fmt="%.5f",
    header="time[ps], rmsd_ref_align_ca_ss[A]",
    delimiter=","
)
np.savetxt(
    fname="rmsd_ref_align_dna.csv",
    X=np.vstack((time, rmsd_ref_align_dna)).T,
    fmt="%.5f",
    header="time[ps], rmsd_ref_align_dna[A]",
    delimiter=","
)
