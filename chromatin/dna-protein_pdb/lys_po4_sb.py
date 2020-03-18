from bionmr_utils.md import *
from glob import glob
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm
import pyxmolpp2

cutoff = 8.0  # cutoff [LYS]N-P[DNA] distance

files = glob("pdbs/????.pdb")
# files = glob("pdbs/1kx5.pdb")
fn_out = "results.csv"

n_p = []
n_op1 = []
n_op2 = []
n_op = []
n_o3 = []
n_o5 = []
ce_n_p = []

for file in tqdm(files):
    frame = PdbFile(file).get_frame()

    ps = frame.asAtoms.filter(aName == "P")
    ps_crd = ps.toCoords
    ps_tree = cKDTree(ps_crd.to_numpy())

    lysns = frame.asAtoms.filter((rName == "LYS") & (aName == "NZ"))
    lysns_crd = lysns.toCoords
    # lysns_tree = cKDTree(lysns_crd.to_numpy())

    dist, idx = ps_tree.query(lysns_crd.to_numpy(), k=1, eps=1e-5, distance_upper_bound=cutoff)

    # print(lysns_crd.to_numpy())
    # print(ps.toCoords.to_numpy())
    # print(len(lysns_crd.to_numpy()), len(ps.toCoords.to_numpy()))
    # print(dist)
    # print(idx)
    # print(len(ps), len(dist), len(idx))

    # find interactions
    for n_i, p_i, d in zip(range(len(lysns)), idx, dist):
        if np.isfinite(d):
            try:
                n_crd = lysns_crd[int(n_i)]
                p_crd = ps_crd[int(p_i)]
                op1_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "OP1")[0].r
                op2_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "OP2")[0].r
                o3_crd = ps[int(p_i) - 1].residue.asAtoms.filter(aName == "O3'")[0].r
                o5_crd = ps[int(p_i)].residue.asAtoms.filter(aName == "O5'")[0].r
                ce_crd = lysns[int(n_i)].residue.asAtoms.filter(aName == "CE")[0].r

                # if calc_angle(ce_crd, n_crd, p_crd).degrees < 80:
                #     continue

                n_p.append(distance(n_crd, p_crd))
                n_op1.append(distance(n_crd, op1_crd))
                n_op2.append(distance(n_crd, op2_crd))
                n_op.append(min(distance(n_crd, op1_crd), distance(n_crd, op2_crd)))
                n_o3.append(distance(n_crd, o3_crd))
                n_o5.append(distance(n_crd, o5_crd))
                ce_n_p.append(calc_angle(ce_crd, n_crd, p_crd).degrees)

                # # save geometries
                # # structure, a_cName, a_rName, a_rId, d_cName, d_rName, d_rId
                # a_cName = ps[int(p_i)].cName.str
                # a_rName = ps[int(p_i)].rName.str
                # a_rId = ps[int(p_i)].rId.serial
                # d_cName = lysns[int(n_i)].cName.str
                # d_rName = lysns[int(n_i)].rName.str
                # d_rId = lysns[int(n_i)].rId.serial
                # donor = lysns[int(n_i)].residue.asAtoms
                # acceptor = ps[int(p_i)].residue.asAtoms
                # acceptor_m1 = ps[int(p_i)].residue.chain.asResidues.filter(rId.is_in({ps[int(p_i)].rId.serial, ps[int(p_i)].rId.serial-1})).asAtoms.filter(aName == "O3'")
                # with open("structures/%s-%s-%s-%d-%s-%s-%d.pdb" % (file.split("/")[-1].rstrip(".pdb"),
                #                                         d_cName,
                #                                         d_rName,
                #                                         d_rId,
                #                                         a_cName,
                #                                         a_rName,
                #                                         a_rId), "w") as f:
                #     donor.to_pdb(f)
                #     acceptor.to_pdb(f)
                #     acceptor_m1.to_pdb(f)
            except pyxmolpp2.polymer.OutOfRangeAtomSelection:
                pass


# print(n_p)
# print(n_op1)
# print(n_op2)
# print(n_o3)
# print(n_o5)
# print(ce_n_p)

results = np.array([n_p, n_op1, n_op2, n_op, n_o3, n_o5, ce_n_p]).T
np.savetxt(fn_out, results, fmt="%.5f", delimiter=",", header="n_p,n_op1,n_op2,n_op,n_o3,n_o5,ce_n_p")
