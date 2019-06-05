from Bio.PDB import *
from Bio.PDB.DSSP import dssp_dict_from_pdb_file


def get_secondary_structure_residues(chain, pdb_code='1KX5'):
    p = PDBList()
    fn = p.retrieve_pdb_file(pdb_code=pdb_code, file_format='pdb', overwrite=False)

    dssp_dict = dssp_dict_from_pdb_file(fn)[0]

    residues = []

    for k in dssp_dict.keys():
        cName = k[0]
        rId = k[1][1]
        DSSP = dssp_dict[k][1]
        if not (DSSP in 'TS-'):
            if cName == chain:
                residues.append(rId)

    return residues
