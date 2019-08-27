import numpy as np
from scipy import spatial
from bionmr_utils.md import *
from tqdm import tqdm

H4_seq = 'SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG'

# setup parameters.
path_to_traj = "../.."
fnout_avg_1 = 'adherence_avg_H4-1_KDtree.txt'
fnout_avg_2 = 'adherence_avg_H4-2_KDtree.txt'
fnout_1 = 'adherence_H4-1_KDtree.txt'
fnout_2 = 'adherence_H4-2_KDtree.txt'
first_dat_file = 51
last_dat_file = 800
n_steps = last_dat_file - first_dat_file + 1
stride = 100

residues_of_interest_1 = set(range(136, 160))
residues_of_interest_2 = set(range(623, 647))

# read trajectory
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/5_run/run00001.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file)

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" %
      (len(traj), len(traj[0].asChains), len(traj[0].asResidues), len(traj[0].asAtoms)))
print("Using run%05d.dat - run%05d.dat" % (first_dat_file, last_dat_file))

# open files and cumulative variables
f_1 = open(fnout_1, 'w')
f_2 = open(fnout_2, 'w')
f_1.write(','.join([H4_seq[i] + str(i+1) for i in range(len(residues_of_interest_1))]) + '\n')
f_2.write(','.join([H4_seq[i] + str(i+1) for i in range(len(residues_of_interest_2))]) + '\n')
min_dists_1 = np.zeros((len(residues_of_interest_1), int(traj.size / stride)))
min_dists_2 = np.zeros((len(residues_of_interest_2), int(traj.size / stride)))

# create predicates
h4_1_ca_pred = rId.is_in(residues_of_interest_1) & (aName == 'CA')
h4_2_ca_pred = rId.is_in(residues_of_interest_2) & (aName == 'CA')

nucleosome_core_pred = ~ (rId.is_in(residues_of_interest_1) | rId.is_in(residues_of_interest_2))
heavy = lambda a: not (a.name.str[0] == 'H')

# run through trajectory and calculate distances
print("Processing frames...")
for frame in tqdm(traj[::stride]):
    if frame.index == 0:
        Ca_1 = frame.asAtoms.filter(h4_1_ca_pred)
        Ca_2 = frame.asAtoms.filter(h4_2_ca_pred)
        nucleosome_core = frame.asAtoms.filter(heavy).filter(nucleosome_core_pred)

    coords_Ca_1 = Ca_1.toCoords
    coords_Ca_2 = Ca_2.toCoords
    coords_nucleosome_core = nucleosome_core.toCoords

    tree_nucleosome_core = spatial.cKDTree(coords_nucleosome_core.to_numpy())

    for i, Ca in enumerate(Ca_1):
        min_dist, idx = tree_nucleosome_core.query([[Ca.r.x, Ca.r.y, Ca.r.z]], k=1, eps=1e-5)
        f_1.write('{:.2f},'.format(min_dist[0]))
        min_dists_1[i, int(frame.index / stride)] = min_dist[0]
    f_1.write('\n')

    for i, Ca in enumerate(Ca_2):
        min_dist, idx = tree_nucleosome_core.query([[Ca.r.x, Ca.r.y, Ca.r.z]], k=1, eps=1e-5)
        f_2.write('{:.2f},'.format(min_dist[0]))
        min_dists_2[i, int(frame.index / stride)] = min_dist[0]
    f_2.write('\n')

f_1.close()
f_2.close()

np.savetxt(fnout_avg_1, np.mean(min_dists_1, axis=1), fmt='%.3f')
np.savetxt(fnout_avg_2, np.mean(min_dists_2, axis=1), fmt='%.3f')
