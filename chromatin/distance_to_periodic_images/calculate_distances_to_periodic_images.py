import numpy as np
from scipy.spatial import cKDTree
from bionmr_utils.md import *
from tqdm import tqdm


# setup parameters.
path_to_traj = "."
fnout = 'distances.txt'
first_dat_file = 1
last_dat_file = 1
n_steps = last_dat_file - first_dat_file + 1
stride = 1

residues_of_interest_1 = set(range(136, 160))
residues_of_interest_2 = set(range(623, 647))


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


# read trajectory
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb="ref.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file,
                          subdir=".",
                          filetype="nc")

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" %
      (len(traj), len(traj[0].asChains), len(traj[0].asResidues), len(traj[0].asAtoms)))
print("Using run%05d.dat - run%05d.dat" % (first_dat_file, last_dat_file))


# open files and cumulative variables
fn = open(fnout, 'w')
fn.write("min_dist_H4-1, min_dis_H4-2\n")


# create predicates
h4_1_pred = rId.is_in(residues_of_interest_1)
h4_2_pred = rId.is_in(residues_of_interest_2)
heavy = lambda a: not (a.name.str[0] == 'H')


# get PBC from XST file and inpcrd file
pbc = get_XST("../..")
v1 = VectorXYZ.from_numpy(pbc[:, 1:4])
v2 = VectorXYZ.from_numpy(pbc[:, 4:7])
v3 = VectorXYZ.from_numpy(pbc[:, 7:10])
scaling_factors = [v1[i].len()/v1[0].len() for i in range(len(v1))]
lat_vec = LatticeVectors(v1[0], v2[0], v3[0])
shift_finder = BestShiftFinder(lat_vec)
inpcrd = "../../1_build/box.inpcrd"
with open(inpcrd) as f:
    content = f.readlines()
content = content[-1].split()[:6]
a, b, c, alpha, beta, gamma = [float(x) for x in content]
distance_upper_bound = np.max(np.array([a/2, b/2, c/2]))


# run through trajectory and calculate distances
print("Processing frames...")
for frame in tqdm(traj[::stride], disable=False):
    if frame.index == 0:
        h4_1 = frame.asResidues.filter(h4_1_pred).asAtoms.filter(heavy)
        h4_2 = frame.asResidues.filter(h4_2_pred).asAtoms.filter(heavy)
        im = frame.asAtoms.filter(heavy)

    # initialize coordinates
    h4_1_coords_np = h4_1.toCoords.to_numpy()
    h4_2_coords_np = h4_2.toCoords.to_numpy()
    n = len(im)
    all_im_coords_np = np.zeros((n * 26, 3))

    # build arrays of coordinates for periodic images
    lat_vec.scale_by(scaling_factors[frame.index])
    counter = 0
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                if not ((i == 0) & (j == 0) & (k == 0)):
                    # get shift and shift_back
                    shift = lat_vec.get_shift(i, j, k)
                    shift_back = lat_vec.get_shift(-i, -j, -k)
                    # apply shift
                    im.transform(Translation3d(shift))
                    # store coordinates
                    im_coords_np = im.toCoords.to_numpy()
                    all_im_coords_np[counter * n:counter * n + n, :] = im_coords_np
                    counter += 1
                    # apply shift_back
                    im.transform(Translation3d(shift_back))
    lat_vec.scale_by(1.0 / scaling_factors[frame.index])

    # build binary trees
    im_tree = cKDTree(all_im_coords_np)

    # query trees and get minimal distances
    dist1, idx1 = im_tree.query(h4_1_coords_np,
                                k=1,
                                eps=1.0,
                                distance_upper_bound=distance_upper_bound,
                                n_jobs=-1)
    h4_1_min_dist = np.min(dist1)
    dist2, idx2 = im_tree.query(h4_2_coords_np,
                                k=1,
                                eps=1.0,
                                distance_upper_bound=distance_upper_bound,
                                n_jobs=-1)
    h4_2_min_dist = np.min(dist2)
    fn.write("%.2f, %.2f\n" % (h4_1_min_dist, h4_2_min_dist))


fn.close()
