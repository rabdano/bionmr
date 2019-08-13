from bionmr_utils.md import *
from tqdm import tqdm
import sys


# setup parameters
path_to_traj = "../.."
first_dat_file = 1
last_dat_file = 600
n_steps = last_dat_file - first_dat_file + 1
stride = 100
phi_out_file = 'phi.csv'
psi_out_file = 'psi.csv'

# residues_of_interest = sorted(list(range(1, 45)) + list(range(136, 160)) + list(range(488, 532)) + list(range(623, 647)))
residues_of_interest = sorted(list(range(136, 160)) + list(range(623, 647)))


# read trajectory
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/5_run/run00001.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file)

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" % (len(traj), len(traj[0].asChains), len(traj[0].asResidues), len(traj[0].asAtoms)))
print("Using run%05d.dat - run%05d.dat" % (first_dat_file, last_dat_file))


# open files
fphi = open(phi_out_file, 'w')
fpsi = open(psi_out_file, 'w')

header = [str(x) for x in residues_of_interest]
fphi.write(','.join([str(x) for x in header]) + '\n')
fpsi.write(','.join([str(x) for x in header]) + '\n')

# run through trajectory and calculate dihedrals
print("Processing frames...")
for frame in tqdm(traj[::stride]):
    if frame.index == 0:
        residues = []
        for resid in residues_of_interest:
            residues.append(frame.asResidues.filter(rId == resid)[0])
        phi_objects = []
        psi_objects = []
        for residue in residues:
            phi = TorsionAngleFactory.phi(residue)
            psi = TorsionAngleFactory.psi(residue)
            phi_objects.append(phi)
            psi_objects.append(psi)

    phi_angles = []
    psi_angles = []
    for phi_object, psi_object in zip(phi_objects, psi_objects):
        if phi_object:
            phi_angles.append(phi_object.value().degrees)
        else:
            phi_angles.append(0.0)
        if psi_object:
            psi_angles.append(psi_object.value().degrees)
        else:
            psi_angles.append(0.0)

    fphi.write(','.join(['%.2f' % x for x in phi_angles]) + '\n')
    fpsi.write(','.join(['%.2f' % x for x in psi_angles]) + '\n')
