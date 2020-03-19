from bionmr_utils.md import *
from glob import glob
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm

cutoff = 8.0  # cutoff N-P[DNA] distance
protein_resnames = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                    "PRO", "SER", "THR", "TRP", "TYR", "VAL"}

files = glob("pdbs/????.pdb")
fn_out_lys_po4 = "lys-po4.csv"
fn_out_nh3_po4 = "nh3-po4.csv"
fn_out_arg_po4 = "arg-po4.csv"

lys_n_op12, lys_ce_n_p, lys_n_p, lys_n_op1, lys_n_op2, lys_n_o3, lys_n_o5 = [], [], [], [], [], [], []
nh3_n_op12, nh3_ca_n_p, nh3_n_p, nh3_n_op1, nh3_n_op2, nh3_n_o3, nh3_n_o5 = [], [], [], [], [], [], []
arg_n12_p, arg_ne_op12, arg_n12_op12, arg_n1_op1, arg_n1_op2, arg_n1_p, arg_n2_op1, arg_n2_op2, arg_n2_p, arg_ne_n1_p, \
    arg_ne_n2_p, arg_ne_cz_p, arg_cz_n1_p, arg_cz_n2_p = [], [], [], [], [], [], [], [], [], [], [], [], [], []
lys_n, lys_ce, lys_p, lys_op1, lys_op2, lys_o3, lys_o5 = [], [], [], [], [], [], []
nh3_n, nh3_ca, nh3_p, nh3_op1, nh3_op2, nh3_o3, nh3_o5 = [], [], [], [], [], [], []
arg_n1, arg_n2, arg_cz, arg_ne, arg_p, arg_op1, arg_op2, arg_o3, arg_o5 = [], [], [], [], [], [], [], [], []
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

            # get atoms
            n = lysns[n_i]
            p = ps[int(p_i)]
            try:
                op1 = ps[int(p_i)].residue.asAtoms.filter(aName == "OP1")[0]
                op2 = ps[int(p_i)].residue.asAtoms.filter(aName == "OP2")[0]
                o3 = ps[int(p_i) - 1].residue.asAtoms.filter(aName == "O3'")[0]
                o5 = ps[int(p_i)].residue.asAtoms.filter(aName == "O5'")[0]
                ce = lysns[n_i].residue.asAtoms.filter(aName == "CE")[0]
            except:
                continue

            # store distances and angles
            lys_n_p.append(distance(n.r, p.r))
            lys_n_op1.append(distance(n.r, op1.r))
            lys_n_op2.append(distance(n.r, op2.r))
            lys_n_op12.append(min(distance(n.r, op1.r), distance(n.r, op2.r)))
            lys_n_o3.append(distance(n.r, o3.r))
            lys_n_o5.append(distance(n.r, o5.r))
            lys_ce_n_p.append(calc_angle(ce.r, n.r, p.r).degrees)

            # store atom IDs
            lys_n.append(n.aId)
            lys_ce.append(ce.aId)
            lys_p.append(p.aId)
            lys_op1.append(op1.aId)
            lys_op2.append(op2.aId)
            lys_o3.append(o3.aId)
            lys_o5.append(o5.aId)

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
            n = nh3ns[n_i]
            p = ps[int(p_i)]
            try:
                op1 = ps[int(p_i)].residue.asAtoms.filter(aName == "OP1")[0]
                op2 = ps[int(p_i)].residue.asAtoms.filter(aName == "OP2")[0]
                o3 = ps[int(p_i) - 1].residue.asAtoms.filter(aName == "O3'")[0]
                o5 = ps[int(p_i)].residue.asAtoms.filter(aName == "O5'")[0]
                ca = nh3ns[int(n_i)].residue.asAtoms.filter(aName == "CA")[0]
            except:
                continue

            # store distances and angles
            nh3_n_p.append(distance(n.r, p.r))
            nh3_n_op1.append(distance(n.r, op1.r))
            nh3_n_op2.append(distance(n.r, op2.r))
            nh3_n_op12.append(min(distance(n.r, op1.r), distance(n.r, op2.r)))
            nh3_n_o3.append(distance(n.r, o3.r))
            nh3_n_o5.append(distance(n.r, o5.r))
            nh3_ca_n_p.append(calc_angle(ca.r, n.r, p.r).degrees)

            # store atom IDs
            nh3_n.append(n.aId)
            nh3_ca.append(ca.aId)
            nh3_p.append(p.aId)
            nh3_op1.append(op1.aId)
            nh3_op2.append(op2.aId)
            nh3_o3.append(o3.aId)
            nh3_o5.append(o5.aId)

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
            p = ps[int(p_i)]

            try:
                n1 = argns[n_i].residue.asAtoms.filter(aName == "NH1")[0]
                n2 = argns[n_i].residue.asAtoms.filter(aName == "NH2")[0]
                op1 = ps[int(p_i)].residue.asAtoms.filter(aName == "OP1")[0]
                op2 = ps[int(p_i)].residue.asAtoms.filter(aName == "OP2")[0]
                o3 = ps[int(p_i) - 1].residue.asAtoms.filter(aName == "O3'")[0]
                o5 = ps[int(p_i)].residue.asAtoms.filter(aName == "O5'")[0]
                cz = argns[n_i].residue.asAtoms.filter(aName == "CZ")[0]
                ne = argns[n_i].residue.asAtoms.filter(aName == "NE")[0]
            except:
                continue

            # store distances and angles
            arg_n12_p.append(min(distance(n1.r, p.r), distance(n2.r, p.r)))
            arg_ne_op12.append(min(distance(ne.r, op1.r), distance(ne.r, op2.r)))
            arg_n12_op12.append(min(distance(n1.r, op1.r),
                                    distance(n1.r, op2.r),
                                    distance(n2.r, op1.r),
                                    distance(n2.r, op2.r)))
            arg_n1_op1.append(distance(n1.r, op1.r))
            arg_n1_op2.append(distance(n1.r, op2.r))
            arg_n1_p.append(distance(n1.r, p.r))
            arg_n2_op1.append(distance(n2.r, op1.r))
            arg_n2_op2.append(distance(n2.r, op2.r))
            arg_n2_p.append(distance(n2.r, p.r))
            arg_ne_n1_p.append(calc_angle(ne.r, n1.r, p.r).degrees)
            arg_ne_n2_p.append(calc_angle(ne.r, n2.r, p.r).degrees)
            arg_ne_cz_p.append(calc_angle(ne.r, cz.r, p.r).degrees)
            arg_cz_n1_p.append(calc_angle(cz.r, n1.r, p.r).degrees)
            arg_cz_n2_p.append(calc_angle(cz.r, n2.r, p.r).degrees)

            # store atom IDs
            arg_n1.append(n1.aId)
            arg_n2.append(n2.aId)
            arg_cz.append(cz.aId)
            arg_ne.append(ne.aId)
            arg_p.append(p.aId)
            arg_op1.append(op1.aId)
            arg_op2.append(op2.aId)
            arg_o3.append(o3.aId)
            arg_o5.append(o5.aId)

            arg_pairs.append(f"{n1.aId}-{n2.aId}-{p.aId}")


# save distances and angles to files
# LYS
results_lys = np.array([lys_n_op12, lys_ce_n_p, lys_n_p, lys_n_op1, lys_n_op2, lys_n_o3, lys_n_o5,
                        lys_n, lys_ce, lys_p, lys_op1, lys_op2, lys_o3, lys_o5]).T
header = "lys_n_op12, lys_ce_n_p, lys_n_p, lys_n_op1, lys_n_op2, lys_n_o3, lys_n_o5, " + \
         "lys_n, lys_ce, lys_p, lys_op1, lys_op2, lys_o3, lys_o5"
fmt = ["%.5f"] * 7 + ["%d"] * 7
np.savetxt(fn_out_lys_po4, results_lys, fmt=fmt, delimiter=",", header=header)

# NH3
results_nh3 = np.array([nh3_n_op12, nh3_ca_n_p, nh3_n_p, nh3_n_op1, nh3_n_op2, nh3_n_o3, nh3_n_o5,
                        nh3_n, nh3_ca, nh3_p, nh3_op1, nh3_op2, nh3_o3, nh3_o5]).T
header = "nh3_n_op12, nh3_ca_n_p, nh3_n_p, nh3_n_op1, nh3_n_op2, nh3_n_o3, nh3_n_o5, " + \
         "nh3_n, nh3_ca, nh3_p, nh3_op1, nh3_op2, nh3_o3, nh3_o5"
fmt = ["%.5f"] * 7 + ["%d"] * 7
np.savetxt(fn_out_nh3_po4, results_nh3, fmt=fmt, delimiter=",", header=header)

# ARG
results_arg = np.array(
    [arg_n12_p, arg_ne_op12, arg_n12_op12, arg_n1_op1, arg_n1_op2, arg_n1_p, arg_n2_op1, arg_n2_op2, arg_n2_p,
     arg_ne_n1_p, arg_ne_n2_p, arg_ne_cz_p, arg_cz_n1_p, arg_cz_n2_p,
     arg_n1, arg_n2, arg_cz, arg_ne, arg_p, arg_op1, arg_op2, arg_o3, arg_o5]).T
header = "arg_n12_p, arg_ne_op12, arg_n12_op12, arg_n1_op1, arg_n1_op2, arg_n1_p, arg_n2_op1, arg_n2_op2, arg_n2_p, " + \
         "arg_ne_n1_p, arg_ne_n2_p, arg_ne_cz_p, arg_cz_n1_p, arg_cz_n2_p, " + \
         "arg_n1, arg_n2, arg_cz, arg_ne, arg_p, arg_op1, arg_op2, arg_o3, arg_o5"
fmt = ["%.5f"] * 14 + ["%d"] * 9

# remove duplicate information for ARG
to_del = []
for i, pair in enumerate(arg_pairs[:-1]):
    if pair == arg_pairs[i + 1]:
        to_del.append(i + 1)
results_arg = np.delete(results_arg, to_del, axis=0)
np.savetxt(fn_out_arg_po4, results_arg, fmt=fmt, delimiter=",", header=header)
