from Bio.PDB import *
from Bio.PDB.Polypeptide import three_to_one

sequences = []
chains = []

p = PDBParser(QUIET=True)
structure = p.get_structure('X', 'pdb_unwrapped/run00001.pdb')
for model in structure:
    for chain in model:
        chains.append(chain.get_id())
        chain_seq = ''
        for residue in chain:
            rName_3 = residue.get_resname()
            if (rName_3 == 'HID') or (rName_3 == 'HIE') or (rName_3 == 'HIP'):
                rName_3 = 'HIS'
            if (rName_3 == 'GLH'):
                rName_3 = 'GLU'
            if (rName_3 == 'ASH'):
                rName_3 = 'ASP'
            if (rName_3 == 'LYN'):
                rName_3 = 'LYS'
            if ('na' in rName_3.lower()) or ('cl' in rName_3.lower()):
                continue
            if rName_3.strip()[0] == 'D':
                chain_seq += rName_3.strip()[1]
            else:
                chain_seq += three_to_one(rName_3)
        sequences.append(chain_seq)

for c, s in zip(chains, sequences):
    print('{} : {}'.format(c, s))
