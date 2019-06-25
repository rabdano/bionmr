import numpy as np
from glob import glob


path_to_pdb = 'pdb_shiftx2'

# common variables
n_residues = 974
h4_1 = [136, 237]
h4_2 = [623, 724]
resids = []
resnames = []


# define AtomNMR class for convenience
class AtomNMR:
    def __init__(self, line):
        s = line.split(',')
        self.cName = str(s[0])
        self.rId = int(s[1])
        self.rName = str(s[2])
        self.aName = str(s[3])
        self.cs = float(s[4])


files = sorted(glob(path_to_pdb + '/run?????.cs'))

H_shifts = np.zeros((len(files), n_residues))
N_shifts = np.zeros((len(files), n_residues))
CA_shifts = np.zeros((len(files), n_residues))
CB_shifts = np.zeros((len(files), n_residues))
C_shifts = np.zeros((len(files), n_residues))
HA_shifts = np.zeros((len(files), n_residues))

for i, file in enumerate(files):
    # read file contents
    fid = open(file, 'r')
    lines = [line.rstrip('\n') for line in fid]
    del lines[0]
    fid.close()

    # create atoms with AtomNMR
    atoms = []
    for line in lines:
        atoms.append(AtomNMR(line))

    # fill in tables of chemical shifts
    current_resid = 0
    rid = -1
    for at in atoms:
        # update residue id
        if current_resid != at.rId:
            current_resid = at.rId
            rid += 1
            # fill resnames and resids variables
            if i == 0:
                resids.append(current_resid)
                resnames.append(at.rName)
        # put chemical shift to corresponding table
        if at.aName == 'H':
            H_shifts[i, rid] = at.cs
        if at.aName == 'N':
            N_shifts[i, rid] = at.cs
        if at.aName == 'CA':
            CA_shifts[i, rid] = at.cs
        if at.aName == 'CB':
            CB_shifts[i, rid] = at.cs
        if at.aName == 'C':
            C_shifts[i, rid] = at.cs
        if at.aName == 'HA':
            HA_shifts[i, rid] = at.cs

str1 = '%3s%6s%10s%10s%10s%10s%10s%10s\n' % ('rId', ' rName', '  mean_H_1', '  mean_N_1', ' mean_CA_1', ' mean_CB_1', '  mean_C_1', ' mean_HA_1')
str2 = '%3s%6s%10s%10s%10s%10s%10s%10s\n' % ('rId', ' rName', '  mean_H_2', '  mean_N_2', ' mean_CA_2', ' mean_CB_2', '  mean_C_2', ' mean_HA_2')
str3 = '%3s%6s%10s%10s%10s%10s%10s%10s\n' % ('rId', ' rName', '    mean_H', '    mean_N', '   mean_CA', '   mean_CB', '    mean_C', '   mean_HA')

h4_1_indexes = [resids.index(h4_1[0]), resids.index(h4_1[1])+1]
h4_2_indexes = [resids.index(h4_2[0]), resids.index(h4_2[1])+1]

h4_1_resids = resids[h4_1_indexes[0]: h4_1_indexes[1]]
h4_1_resnames = resnames[h4_1_indexes[0]: h4_1_indexes[1]]

h4_2_resids = resids[h4_2_indexes[0]: h4_2_indexes[1]]
h4_2_resnames = resnames[h4_2_indexes[0]: h4_2_indexes[1]]

for rId1, rId2, rName1, rName2 in zip(h4_1_resids, h4_2_resids, h4_1_resnames, h4_2_resnames):
    i1 = resids.index(rId1)
    i2 = resids.index(rId2)
    H1, H2 = H_shifts[:, i1], H_shifts[:, i2]
    N1, N2 = N_shifts[:, i1], N_shifts[:, i2]
    CA1, CA2 = CA_shifts[:, i1], CA_shifts[:, i2]
    CB1, CB2 = CB_shifts[:, i1], CB_shifts[:, i2]
    C1, C2 = C_shifts[:, i1], C_shifts[:, i2]
    HA1, HA2 = HA_shifts[:, i1], HA_shifts[:, i2]

    mean_H_1 = np.mean(H1[np.nonzero(H1)])
    mean_H_2 = np.mean(H2[np.nonzero(H2)])
    mean_N_1 = np.mean(N1[np.nonzero(N1)])
    mean_N_2 = np.mean(N2[np.nonzero(N2)])

    mean_CA_1 = np.mean(CA1[np.nonzero(CA1)])
    mean_CA_2 = np.mean(CA2[np.nonzero(CA2)])
    mean_CB_1 = np.mean(CB1[np.nonzero(CB1)])
    mean_CB_2 = np.mean(CB2[np.nonzero(CB2)])

    mean_C_1 = np.mean(C1[np.nonzero(C1)])
    mean_C_2 = np.mean(C2[np.nonzero(C2)])
    mean_HA_1 = np.mean(HA1[np.nonzero(HA1)])
    mean_HA_2 = np.mean(HA2[np.nonzero(HA2)])
    
    str1 += '%3d%6s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n' % \
        (rId1, rName1, mean_H_1, mean_N_1, mean_CA_1, mean_CB_1, mean_C_1, mean_HA_1)
    str2 += '%3d%6s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n' % \
        (rId2, rName2, mean_H_2, mean_N_2, mean_CA_2, mean_CB_2, mean_C_2, mean_HA_2)
    str3 += '%3d%6s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n' % \
        (rId1-135, rName1,
         0.5*(mean_H_1 + mean_H_2), 0.5*(mean_N_1 + mean_N_2),
         0.5*(mean_CA_1 + mean_CA_2), 0.5*(mean_CB_1 + mean_CB_2),
         0.5*(mean_C_1 + mean_C_2), 0.5*(mean_HA_1 + mean_HA_2))

with open('shiftx2_report.txt', 'w') as f:
    f.write('First tail:\n')
    f.write(str1)
    f.write('Second tail:\n')
    f.write(str2)
    f.write('Mean for two tails:\n')
    f.write(str3)
