from bionmr_utils.md import *
from get_secondary_structure_residues import *
from tqdm import tqdm
import numpy as np

# setup trajectory parameters
n_dat_files = 500  # ns
stride = 1  # ps
residues_of_interest = set(list(range(1, 45)) + list(range(136, 160)) + list(range(488, 532)) + list(range(623, 647)))
chain_letters = 'ABEF'

# load trajectory
path_to_traj = "../.."
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/5_run/run00001.pdb",
                          stride=1,
                          first=1,
                          last=n_dat_files)

# get residue ids for secondary structure from DSSP using Biopython and PDB
protein_chains = 'ABCDEFGH'
ss_residues = []
incr = 0
for c in protein_chains:
    ss = get_secondary_structure_residues(c, pdb_code='1KX5')
    ss_residues.append(set([(x+incr) for x in ss]))
    incr += len(ref.asResidues.filter(cName == c))
del ss_residues[7]
del ss_residues[6]
del ss_residues[3]
del ss_residues[2]

# reference atoms for alignment
# DNA_ref = ref.asAtoms.filter(cName.is_in({'I', 'J'}) & (aName == 'P'))
align_ref = [ref.asAtoms.filter((aName == 'CA') & (rId.is_in(ss_res))) for ss_res in ss_residues]


# atoms for RMSF calculation
ca_ref = [ref.asAtoms.filter((aName == 'CA') & (cName == ChainName(c)) & (rId.is_in(residues_of_interest))) for c in chain_letters]


# initialize average coordinates with (0,0,0)
avg_coords = [c.toCoords.transform(UniformScale3d(0)) for c in ca_ref]


# calculate average coordinates across N frames
print('calculate average coordinates across N frames')
skipped_frames = set()
first = True
for frame in tqdm(traj[::stride]):
    if first:
        # DNA = frame.asAtoms.filter(cName.is_in({'I', 'J'}) & (aName == 'P'))
        check_atom_1 = frame.asAtoms[0]
        check_atom_2 = frame.asAtoms[1]
        align = [frame.asAtoms.filter((aName == 'CA') & (rId.is_in(ss_res))) for ss_res in ss_residues]
        ca = [frame.asAtoms.filter((aName == 'CA') & (cName == ChainName(c)) & (rId.is_in(residues_of_interest))) for c in chain_letters]
        first = False

    # # align to DNA
    # alignment = DNA.alignment_to(DNA_ref)
    # frame.asAtoms.transform(alignment)

    # skip corrupted frames
    d = distance(check_atom_1, check_atom_2)
    if d < 0.5:
        skipped_frames.add(frame.index)
        continue

    for i, chain in enumerate(ca):

        # align to secondary structured CA
        alignment = align[i].alignment_to(align_ref[i])
        frame.asAtoms.transform(alignment)

        # sum all coordinates
        for j, a in enumerate(chain):
            avg_coords[i][j] += a.r

# 'divide' by number of frames
print('Frames skipped:\n')
print(skipped_frames)
for crds in avg_coords:
    crds.transform(UniformScale3d(1.0 / (traj.size - len(skipped_frames)) * stride))


# align to average coordinates and calculate RMSF for each CA
print('align to average coordinates and calculate RMSF for each CA (per residue)')
rmsf = [np.zeros((len(x),)) for x in ca_ref]
first = True
for frame in tqdm(traj[::stride]):
    if first:
        # DNA = frame.asAtoms.filter(cName.is_in({'I', 'J'}) & (aName == 'P'))
        align = [frame.asAtoms.filter((aName == 'CA') & (rId.is_in(ss_res))) for ss_res in ss_residues]
        ca = [frame.asAtoms.filter((aName == 'CA') & (cName == ChainName(c)) & (rId.is_in(residues_of_interest))) for c in chain_letters]
        first = False

    # # align to DNA
    # alignment = DNA.alignment_to(DNA_ref)
    # frame.asAtoms.transform(alignment)

    # skip corrupted frames
    if frame.index in skipped_frames:
        continue

    for i, chain in enumerate(ca):
        # # align to average coordinates
        # chain.align_to(avg_coords[i])

        # align to secondary structured CA
        alignment = align[i].alignment_to(align_ref[i])
        frame.asAtoms.transform(alignment)

        # calculate RMSF
        for j, a in enumerate(chain):
            rmsf[i][j] += (a.r - avg_coords[i][j]).len2()

rmsf = [np.sqrt(rmsf_i / (traj.size - len(skipped_frames)) * stride) for rmsf_i in rmsf]


# write RMSF to file
for i, rmsf_i in enumerate(rmsf):
    with open('RMSF_chain_' + chain_letters[i] + '.csv', 'w') as f:
        f.write('rId,rName,RMSF\n')
        for j, at in enumerate(ca_ref[i]):
            f.write('%d,%s,%.8f\n' % (at.rId.serial, at.rName.str, rmsf_i[j]))
