from bionmr_utils.md import *
from glob import glob
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm
import pyxmolpp2

cutoff = 8.0  # cutoff N-P[DNA] distance
protein_resnames = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                    "PRO", "SER", "THR", "TRP", "TYR", "VAL"}

files = glob("pdbs/????.pdb")
fn_out_lys_po4 = "lys-po4.csv"
fn_out_nh3_po4 = "nh3-po4.csv"
fn_out_arg_po4 = "arg-po4.csv"

lys_n_p, lys_ce_n_p, lys_n_op1, lys_n_op2, lys_n_op12, lys_n_o3, lys_n_o5 = [], [], [], [], [], [], []
nh3_n_p, nh3_ca_n_p, nh3_n_op1, nh3_n_op2, nh3_n_op12, nh3_n_o3, nh3_n_o5 = [], [], [], [], [], [], []
arg_n12_p, arg_ne_op12, arg_n1_op1, arg_n1_op2, arg_n1_p, arg_n2_op1, arg_n2_op2, arg_n2_p = [], [], [], [], [], [], [], []
arg_ne_n1_p, arg_ne_n2_p, arg_ne_cz_p, arg_cz_n1_p, arg_cz_n2_p = [], [], [], [], []
arg_pairs = []


for file in tqdm(files, disable=False):
    frame = PdbFile(file).get_frame()

    #
    #  LYS - PO4
    #
    ps = frame.asAtoms.filter(aName == "P")
    ps_crd = ps.toCoords
    ps_tree = cKDTree(ps_crd.to_numpy())

    lysns = frame.asAtoms.filter((rName == "LYS") & (aName == "NZ"))
    lysns_crd = lysns.toCoords

    dist, idx = ps_tree.query(lysns_crd.to_numpy(), k=1, eps=1e-5, distance_upper_bound=cutoff)

    for n_i, p_i, d in zip(range(len(lysns)), idx, dist):
        if np.isfinite(d):
            # do not consider 5'-terminal nucleotide
            if (p_i > len(ps)) | ((p_i - 1) < 0):
                continue

            # get atom coordinates
            n_crd = lysns_crd[int(n_i)]
            p_crd = ps_crd[int(p_i)]
            try:
                op1_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "OP1")[0].r
                op2_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "OP2")[0].r
                o3_crd = ps[int(p_i) - 1].residue.asAtoms.filter(aName == "O3'")[0].r
                o5_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "O5'")[0].r
                ce_crd = lysns[int(n_i)].residue.asAtoms.filter(aName == "CE")[0].r
            except:
                continue

            # store distances and angles
            lys_n_p.append(distance(n_crd, p_crd))
            lys_n_op1.append(distance(n_crd, op1_crd))
            lys_n_op2.append(distance(n_crd, op2_crd))
            lys_n_op12.append(min(distance(n_crd, op1_crd), distance(n_crd, op2_crd)))
            lys_n_o3.append(distance(n_crd, o3_crd))
            lys_n_o5.append(distance(n_crd, o5_crd))
            lys_ce_n_p.append(calc_angle(ce_crd, n_crd, p_crd).degrees)

    #
    #  NH3
    #
    residues = []
    for chain in frame.asChains:
        residues.append(chain.asResidues[0])
    ats = ResidueSelection(residues).asAtoms.filter((aName == "N") & (rName.is_in(protein_resnames)))
    nh3ns = AtomSelection(ats)
    nh3ns_crd = nh3ns.toCoords

    dist, idx = ps_tree.query(nh3ns_crd.to_numpy(), k=1, eps=1e-5, distance_upper_bound=cutoff)

    for n_i, p_i, d in zip(range(len(nh3ns)), idx, dist):
        if np.isfinite(d):
            # do not consider 5'-terminal nucleotide
            if (p_i > len(ps)) | ((p_i - 1) < 0):
                continue

            # get atom coordinates
            n_crd = nh3ns_crd[int(n_i)]
            p_crd = ps_crd[int(p_i)]
            try:
                op1_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "OP1")[0].r
                op2_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "OP2")[0].r
                o3_crd = ps[int(p_i) - 1].residue.asAtoms.filter(aName == "O3'")[0].r
                o5_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "O5'")[0].r
                ca_crd = nh3ns[int(n_i)].residue.asAtoms.filter(aName == "CA")[0].r
            except:
                continue

            # store distances and angles
            nh3_n_p.append(distance(n_crd, p_crd))
            nh3_n_op1.append(distance(n_crd, op1_crd))
            nh3_n_op2.append(distance(n_crd, op2_crd))
            nh3_n_op12.append(min(distance(n_crd, op1_crd), distance(n_crd, op2_crd)))
            nh3_n_o3.append(distance(n_crd, o3_crd))
            nh3_n_o5.append(distance(n_crd, o5_crd))
            nh3_ca_n_p.append(calc_angle(ca_crd, n_crd, p_crd).degrees)

    #
    #  ARG
    #
    argns = frame.asAtoms.filter((rName == "ARG") & (aName.is_in({"NH1", "NH2"})))
    argns_crd = argns.toCoords

    dist, idx = ps_tree.query(argns_crd.to_numpy(), k=1, eps=1e-5, distance_upper_bound=cutoff)

    for n_i, p_i, d in zip(range(len(argns)), idx, dist):
        if np.isfinite(d):
            # do not consider 5'-terminal nucleotide
            if (p_i > len(ps)) | ((p_i - 1) < 0):
                continue

            # get atom coordinates
            p_crd = ps_crd[int(p_i)]
            if n_i % 2 == 0:
                n1_i = n_i
                n2_i = n_i + 1
                n1_crd = argns_crd[n1_i]
                n2_crd = argns_crd[n2_i]
            else:
                n1_i = n_i - 1
                n2_i = n_i
                n1_crd = argns_crd[n1_i]
                n2_crd = argns_crd[n2_i]
            try:
                op1_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "OP1")[0].r
                op2_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "OP2")[0].r
                o3_crd = ps[int(p_i) - 1].residue.asAtoms.filter(aName == "O3'")[0].r
                o5_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "O5'")[0].r
                cz_crd = argns[n1_i].residue.asAtoms.filter(aName == "CZ")[0].r
                ne_crd = argns[n1_i].residue.asAtoms.filter(aName == "NE")[0].r
            except:
                continue

            # store distances and angles
            arg_n12_p.append(min(distance(n1_crd, p_crd), distance(n2_crd, p_crd)))
            arg_ne_op12.append(min(distance(ne_crd, op1_crd), distance(ne_crd, op2_crd)))
            arg_n1_op1.append(distance(n1_crd, op1_crd))
            arg_n1_op2.append(distance(n1_crd, op2_crd))
            arg_n1_p.append(distance(n1_crd, p_crd))
            arg_n2_op1.append(distance(n2_crd, op1_crd))
            arg_n2_op2.append(distance(n2_crd, op2_crd))
            arg_n2_p.append(distance(n2_crd, p_crd))
            arg_ne_n1_p.append(calc_angle(ne_crd, n1_crd, p_crd).degrees)
            arg_ne_n2_p.append(calc_angle(ne_crd, n2_crd, p_crd).degrees)
            arg_ne_cz_p.append(calc_angle(ne_crd, cz_crd, p_crd).degrees)
            arg_cz_n1_p.append(calc_angle(cz_crd, n1_crd, p_crd).degrees)
            arg_cz_n2_p.append(calc_angle(cz_crd, n2_crd, p_crd).degrees)

            arg_pairs.append(f"{argns[n1_i].cName.str}{argns[n1_i].rId.serial}")

# save distances and angles to files
results_lys = np.array([lys_n_p, lys_ce_n_p, lys_n_op1, lys_n_op2, lys_n_op12, lys_n_o3, lys_n_o5]).T
np.savetxt(fn_out_lys_po4, results_lys, fmt="%.5f", delimiter=",",
           header="lys_n_p, lys_ce_n_p, lys_n_op1, lys_n_op2, lys_n_op12, lys_n_o3, lys_n_o5")

results_nh3 = np.array([nh3_n_p, nh3_ca_n_p, nh3_n_op1, nh3_n_op2, nh3_n_op12, nh3_n_o3, nh3_n_o5]).T
np.savetxt(fn_out_nh3_po4, results_nh3, fmt="%.5f", delimiter=",",
           header="nh3_n_p, nh3_ca_n_p, nh3_n_op1, nh3_n_op2, nh3_n_op12, nh3_n_o3, nh3_n_o5")

results_arg = np.array(
    [arg_n12_p, arg_ne_op12, arg_n1_op1, arg_n1_op2, arg_n1_p, arg_n2_op1, arg_n2_op2, arg_n2_p, arg_ne_n1_p,
     arg_ne_n2_p, arg_ne_cz_p, arg_cz_n1_p, arg_cz_n2_p]).T
# remove duplicate information for ARG
to_del = []
for i, pair in enumerate(arg_pairs[:-1]):
    if pair == arg_pairs[i + 1]:
        to_del.append(i + 1)
results_arg = np.delete(results_arg, to_del, axis=0)
np.savetxt(fn_out_arg_po4, results_arg, fmt="%.5f", delimiter=",",
           header=("arg_n12_p, arg_ne_op12, arg_n1_op1, arg_n1_op2, arg_n1_p, arg_n2_op1, arg_n2_op2, arg_n2_p," +
                   " arg_ne_n1_p, arg_ne_n2_p, arg_ne_cz_p, arg_cz_n1_p, arg_cz_n2_p"))
