import numpy as np
from glob import glob
import struct

# common variables
path_to_pdb = 'pdb_shifts'
n_residues = 204
n_tail_residues = int(n_residues/2)
residue_map = list(range(136, 136+102)) + list(range(623, 623+102))
qdb_header = 2
resnames = []


# define AtomNMR class for convenience
class AtomNMR:
    def __init__(self, line, emp_qdb='emp'):
        if emp_qdb == 'emp':
            s = line.split()
            # .emp description
            # s[0] - basename
            # s[1] - rName
            # s[2] - rId
            # s[3] - aName
            # s[4] - RingCur
            # s[5] - El
            # s[6] - P_anis
            # s[7] - Const
            # s[8] - RC
            # s[9] - pred
            self.basename = str(s[0])
            self.rName = str(s[1])
            self.rId = int(s[2])
            self.aName = str(s[3])
            self.RingCur = float(s[4])
            self.El = float(s[5])
            self.P_anis = float(s[6])
            self.Const = float(s[7])
            self.RC = float(s[8])
            self.pred = float(s[9])
            self.cs = float(s[9])
        else:
            # .qdb description
            # PDB cf residue atom     bb_p  bb_s  bb_f  chi_p chi_s  HB-D  HB-I  REF         pred
            # s[0] - basename
            # s[1] - cf
            # s[2] - rName
            # s[3] - rId
            # s[4] - aName
            # s[5] - bb_p
            # s[6] - bb_s
            # s[7] - bb_f
            # s[8] - chi_p
            # s[9] - chi_s
            # s[10] - HB-D
            # s[11] - HB-I
            # s[12] - REF
            # s[13] - pred
            # fmt = '%-5s%-2s%-5s%-4d%-3s    %6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%8.2f    %8.2f'
            # [0:5]
            # [5:7]
            # [7:12]
            # [12:16]
            # [16:19]
            # [23:29]
            # [29:35]
            # [35:41]
            # [41:47]
            # [47:53]
            # [53:59]
            # [59:64]
            # [64:72]
            # [76:84]
            self.basename = str(line[0:5].strip())
            self.cf = str(line[5:7].strip())
            self.rName = str(line[7:12].strip())
            self.rId = int(line[12:16].strip())
            self.aName = str(line[16:19].strip())
            self.bb_p = float(line[23:29].strip())
            self.bb_s = float(line[29:35].strip())
            self.bb_f = float(line[36:41].strip())
            self.chi_p = float(line[42:47].strip())
            self.chi_s = float(line[48:53].strip())
            self.HB_D = float(line[54:59].strip())
            self.HB_I = float(line[60:65].strip())
            self.REF = float(line[66:73].strip())
            self.cs = float(line[78:85].strip())


basenames = [x[:-4] for x in sorted(glob(path_to_pdb + '/run?????.pdb'))]

# create empty arrays for shifts storage
H_shifts = np.zeros((len(basenames), n_residues))
N_shifts = np.zeros((len(basenames), n_residues))
CA_shifts = np.zeros((len(basenames), n_residues))
CB_shifts = np.zeros((len(basenames), n_residues))
C_shifts = np.zeros((len(basenames), n_residues))
HA_shifts = np.zeros((len(basenames), n_residues))

# extract chemical shifts
for i, file in enumerate(basenames):
    # print(i,file)
    # read file contents
    emp_all = open(file + '.emp', 'r')
    emp_protein = open(file + '_p.emp', 'r')
    qdb_all = open(file + '_cn.qdb', 'r')
    qdb_protein = open(file + "_p_cn.qdb","r")
    emp_all_lines = [line.rstrip('\n') for line in emp_all]
    emp_protein_lines = [line.rstrip('\n') for line in emp_protein]
    qdb_all_lines = [line.rstrip('\n') for line in qdb_all]
    qdb_protein_lines = [line.rstrip("\n") for line in qdb_protein]

    emp_all_header = 0
    for j, line in enumerate(emp_all_lines):
        if '=====' in line:
            emp_all_header = j

    emp_protein_header = 0
    for j, line in enumerate(emp_protein_lines):
        if '=====' in line:
            emp_protein_header = j

    del emp_all_lines[0:emp_all_header + 1]
    del emp_protein_lines[0:emp_protein_header + 1]
    del qdb_all_lines[0:qdb_header]
    del qdb_protein_lines[0:qdb_header]
    emp_all.close()
    emp_protein.close()
    qdb_all.close()
    qdb_protein.close()

    # create atoms: one line - one atom
    atoms_all = []
    atoms_protein = []
    atoms_cn_all = []
    atoms_cn_protein = []
    for line in emp_all_lines:
        atoms_all.append(AtomNMR(line, emp_qdb="emp"))
    for line in emp_protein_lines:
        atoms_protein.append(AtomNMR(line, emp_qdb="emp"))
    for line in qdb_all_lines:
        atoms_cn_all.append(AtomNMR(line, emp_qdb="qdb"))
    for line in qdb_protein_lines:
        atoms_cn_protein.append(AtomNMR(line, emp_qdb="qdb"))

    # fill in tables
    # proton chemical shifts
    for a_all, a_protein in zip(atoms_all, atoms_protein):
        # print(a_all.rId, a_protein.rId)
        assert residue_map.index(a_all.rId) == residue_map.index(a_protein.rId)
        index = residue_map.index(a_all.rId)

        # combine 5 contributions for protons, replace random coil shift
        a_all.cs = a_all.RC  # Random coil (constant term)
        a_all.cs += a_all.RingCur  # Ring current (DNA and protein contibution)
        a_all.cs += a_protein.El  # Electorstatic (only protein contibution!!!)
        a_all.cs += a_all.P_anis  # Anisotropy (depends on geometry)
        a_all.cs += a_protein.Const  # Constant term (atom type specific)

        a_all.cs = a_all.cs - a_protein.cs

        # put chemical shift to corresponding table
        if a_all.aName == 'H':
            H_shifts[i, index] = a_all.cs
        if a_all.aName == 'HA':
            HA_shifts[i, index] = a_all.cs

    # carbon and nitrogen
    for a_all, a_protein in zip(atoms_cn_all, atoms_cn_protein):
        if a_all.rId in residue_map:
            assert residue_map.index(a_all.rId) == residue_map.index(a_protein.rId)
            index = residue_map.index(a_all.rId)

            # combine contributions, replace random coil shift
            a_all.cs = a_all.REF  # reference chemical shift
            a_all.cs += a_all.bb_p  # preceding backbone effect
            a_all.cs += a_all.bb_s  # self backbone effect
            a_all.cs += a_all.bb_f  # following backbone effect
            a_all.cs += a_all.chi_p  # preceding chi effect
            a_all.cs += a_all.chi_s  # self chi effect
            a_all.cs += a_all.HB_D  # direct HB effect
            a_all.cs += a_all.HB_I  # indirect HB effect

            a_all.cs = a_all.cs - a_protein.cs

            # put chemical shift to corresponding table
            if a_all.aName == "N":
                N_shifts[i, index] = a_all.cs
            if a_all.aName == "Ca":
                CA_shifts[i, index] = a_all.cs
            if a_all.aName == "Cb":
                CB_shifts[i, index] = a_all.cs
            if a_all.aName == "CO":
                C_shifts[i, index] = a_all.cs

    if i == 0:
        for rId in residue_map:
            for at in (atoms_all + atoms_cn_all):
                if at.rId == rId:
                    resnames.append(at.rName)
                    break


# for shift in H_shifts:
#     print(shift)
# TODO
# 1. mean for residue across all frames
# 2. *historgram for residue for each shift
# 3. mean for two tails
str1 = '%3s%6s%10s%10s%10s%10s%10s%10s\n' % ('rId', ' rName', '  mean_H_1', '  mean_N_1', ' mean_CA_1', ' mean_CB_1', '  mean_C_1', ' mean_HA_1')
str2 = '%3s%6s%10s%10s%10s%10s%10s%10s\n' % ('rId', ' rName', '  mean_H_2', '  mean_N_2', ' mean_CA_2', ' mean_CB_2', '  mean_C_2', ' mean_HA_2')
str3 = '%3s%6s%10s%10s%10s%10s%10s%10s\n' % ('rId', ' rName', '    mean_H', '    mean_N', '   mean_CA', '   mean_CB', '    mean_C', '   mean_HA')
for i, rId_1, rId_2, rName in zip(range(n_tail_residues), residue_map[:n_tail_residues], residue_map[n_tail_residues:], resnames[:n_tail_residues]):
    H1, H2 = H_shifts[:, i], H_shifts[:, i+n_tail_residues]
    N1, N2 = N_shifts[:, i], N_shifts[:, i+n_tail_residues]
    CA1, CA2 = CA_shifts[:, i], CA_shifts[:, i+n_tail_residues]
    CB1, CB2 = CB_shifts[:, i], CB_shifts[:, i+n_tail_residues]
    C1, C2 = C_shifts[:, i], C_shifts[:, i+n_tail_residues]
    HA1, HA2 = HA_shifts[:, i], HA_shifts[:, i+n_tail_residues]

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
        (rId_1, rName, mean_H_1, mean_N_1, mean_CA_1, mean_CB_1, mean_C_1, mean_HA_1)
    str2 += '%3d%6s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n' % \
        (rId_2, rName, mean_H_2, mean_N_2, mean_CA_2, mean_CB_2, mean_C_2, mean_HA_2)
    str3 += '%3d%6s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n' % \
        (rId_1-135, rName,
         0.5*(mean_H_1 + mean_H_2), 0.5*(mean_N_1 + mean_N_2),
         0.5*(mean_CA_1 + mean_CA_2), 0.5*(mean_CB_1 + mean_CB_2),
         0.5*(mean_C_1 + mean_C_2), 0.5*(mean_HA_1 + mean_HA_2))

with open('shifts_diff.txt', 'w') as f:
    f.write('First tail:\n')
    f.write(str1)
    f.write('Second tail:\n')
    f.write(str2)
    f.write('Mean for two tails:\n')
    f.write(str3)
