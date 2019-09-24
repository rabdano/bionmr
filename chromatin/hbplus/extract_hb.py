import numpy as np
from bionmr_utils.md import *
from glob import glob
import errno
import os
import sys
from subprocess import call
from tqdm import tqdm
import time

DEVNULL = open(os.devnull, 'wb')

# setup parameters.
path_to_traj = "../.."
first_dat_file = 51
last_dat_file = 800
n_steps = last_dat_file - first_dat_file + 1
stride = 100

# probe atoms
probe = ((cName == 'J') & (aName == 'P'))


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


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


def parse_hbplus_output(hbplus_out):
    #============================================================================
    #                          Table I: *.hb2 format
    #
    # 01-13 Donor Atom, including . . .
    #
    # 01    Chain ID (defaults to '-')
    # 02-05 Residue Number
    # 06    Insertion Code (defaults to '-')
    # 07-09 Amino Acid Three Letter Code
    # 10-13 Atom Type Four Letter Code
    #
    # 15-27 Acceptor Atom, same format as Donor atom
    # 28-32 Donor - Acceptor distance, in Angstroms
    # 34-35 Atom Categories - M(ain-chain), S(ide-chain) or H(etatm) - of D & A
    # 37-39 Gap between donor and acceptor groups, in amino acids
    #       (-1 if not applicable)
    # 41-45 Distance between the CA atoms of the donor and acceptor residues
    #       (-1 if one of the two atoms is in a hetatm)
    # 47-51 Angle formed by the Donor and Acceptor at the hydrogen, in degrees.
    #       (-1 if the hydrogen is not defined)
    # 53-57 Distance between the hydrogen and the Acceptor, in Angstroms
    #       (-1 if the hydrogen is not defined)
    # 59-63 The smaller angle at the Acceptor formed by the hydrogen and an
    #       acceptor antecedent (-1 if the hydrogen, or the acceptor antecedent,
    #       is not defined)
    # 65-69 The smaller angle at the Acceptor formed by the donor and an acceptor
    #       antecedent (-1 if not applicable)
    # 71-75 Count of hydrogen bonds
    #============================================================================

    header = 8
    f = open(hbplus_out, 'r')
    hbs = set()
    for i, line in enumerate(f):
        if i >= header:
            d_cId = line[0].strip()
            d_rId = int(line[1:5].strip())
            d_rName = line[6:9].strip()
            d_aName = line[9:13].strip()
            a_cId = line[14].strip()
            a_rId = int(line[15:19].strip())
            a_rName = line[20:23].strip()
            a_aName = line[23:27].strip()
            hbs.add('{}~{}::{}--{}~{}::{}'.format(d_rName, d_rId, d_aName, a_rName, a_rId, a_aName))
    f.close()
    return hbs


# create output dir
output_dir = 'pdb_unwrapped'
mkdir_p(output_dir)

# read trajectory
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/1_build/ref.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file)

chain_ids = [chain.name.str for chain in ref.asChains]

for i, chain in enumerate(ref.asAtoms.asChains):
    print(i, chain.cName, len(chain.asResidues))

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

hydrogen_bonds = list()

start = time.time()
# run through trajectory
for frame in tqdm(traj[::stride]):
    # sys.stdout.write(' '*40 + '\r' + 'Frames: {:d}/{:d}:'.format(int((frame.index/stride)+1), int(traj.size/stride)))
    if frame.index == 0:
        ref_ats = ref.asAtoms
        frame_ats = frame.asAtoms
        ref_align_atoms = ref_ats.filter(probe)
        frame_align_atoms = frame_ats.filter(probe)
        # align reference by first frame nucleic P
        ref_ats.transform(ref_align_atoms.alignment_to(frame_align_atoms))
        ats_ref = []
        ats_frame = []
        for cid in sorted(chain_ids):
            ats_ref.append(ref_ats.filter(cName == cid))
            ats_frame.append(frame_ats.filter(cName == cid))


    # align all atoms by one of DNA chains
    ref_ats.transform(ref_align_atoms.alignment_to(frame_align_atoms))

    # cycle through chains and check first atom displacement
    shift_finder.scale_lattice_by(scaling_factors[frame.index])

    # sys.stdout.write('\t\t')
    for i, cid in enumerate(sorted(chain_ids)):
        # sys.stdout.write(cid + ' ')
        # sys.stdout.flush()

        dif, shift = shift_finder.find_best_shift(ats_ref[i].geom_center(),
                                                  ats_frame[i].geom_center())

        ats_frame[i].transform(Translation3d(shift))

    shift_finder.scale_lattice_by(1.0/scaling_factors[i])

    # write updated pdb
    fn = '{:07d}.pdb'.format(frame.index)
    with open(output_dir + '/' + fn, 'w') as f:
        a = v1[i].len()
        b = v2[i].len()
        c = v3[i].len()
        fmt = 'CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%10s%4s\n'
        f.write(fmt % (a, b, c, alpha, beta, gamma, ' ' * 10, ' ' * 4))
        frame_ats[0].frame.to_pdb(f)
sys.stdout.write('\n')

files = sorted(glob(output_dir + '/*.pdb'))
for i, file in tqdm(enumerate(files)):
    # sys.stdout.write(' '*40 + '\r' + 'HBPLUS: {}/{}\r'.format(i, len(files)))
    fn = file[-11:]
    cmd = ['/home/seva/bin/hbplus/hbplus']
    cmd += [fn]
    call(cmd, cwd=output_dir, stdout=DEVNULL)

    frame_hbs = parse_hbplus_output(file[:-4] + '.hb2')
    hydrogen_bonds.append(frame_hbs)
sys.stdout.write('\n')
end = time.time()
print('Elapsed time: ', end - start)


unique_hydrogen_bonds = set()
for frame_hbs in hydrogen_bonds:
    unique_hydrogen_bonds |= frame_hbs

unique_hydrogen_bonds = sorted(list(unique_hydrogen_bonds), key=lambda hb: int(hb.split('::')[0].split('~')[1]))

hb_trace = np.zeros((len(hydrogen_bonds), len(unique_hydrogen_bonds)), dtype=int)

for i, unique_hb in enumerate(unique_hydrogen_bonds):
    for j, frame_hbs in enumerate(hydrogen_bonds):
        if unique_hb in frame_hbs:
            hb_trace[j, i] = 1

np.savetxt('hb_trace.dat', hb_trace, fmt='%d', header=' '.join(unique_hydrogen_bonds), comments='')

call(['rm', '-rf', 'output_dir'])
