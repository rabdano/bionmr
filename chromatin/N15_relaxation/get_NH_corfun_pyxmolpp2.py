import numpy as np
from pyxmolpp2.polymer import *
from pyxmolpp2.pdb import *
from pyxmolpp2.geometry import *
from pyxmolpp2_lib_sr import *
import os
import sys
import pickle

altered_records = AlteredPdbRecords(StandardPdbRecords.instance())
altered_records.alter_record(RecordName("ATOM"),FieldName("serial"),[7,12])


# setup parameters.
path_to_traj = "../.."
n_steps = ___N_DAT_FILES___
stride = ___STRIDE___
residues_of_interest = set(list(range(1,45)) + list(range(136,160)) + list(range(488,532)) + list(range(623,647)))
cut_autocorr_function = n_steps #ns
HN_mask = "___HN_MASK___"
align_dna = True

# read trajectory
ref = PdbFile(path_to_traj + "/5_run/run00001.pdb", altered_records).get_frame()
traj = make_trajectory(path_to_traj, ref, first=(___FIRST_DAT_FILE___ - 1), limit=(n_steps + ___FIRST_DAT_FILE___ - 1))
print("Trajectory contain %d frames, %d residues / %d atoms in each." % (len(traj), len(traj[0].asResidues), len(traj[0].asAtoms)))

# create folders and open files
resids = []
resnames = []
for at in ref.asAtoms.filter((aName == HN_mask) & (rId.is_in(residues_of_interest))):
    resids.append(at.rId.serial)
    resnames.append(at.rName.str)
resi = tuple(zip(resids, resnames))
print("Autocorrelation functions will be calculated for following residues:")
print(resi)

if not os.path.exists("cor_NH_n%d_s%d" % (n_steps, stride)): os.makedirs("cor_NH_n%d_s%d" % (n_steps, stride))
# vectors - list of VectXYZ
vectors = []

first_time = True
frame_id = 0
is_nucleic = lambda a: a.name.str == "P"

# run through trajectory and calculate vectors
print("Processing frames...")
for frame in traj[::stride]:
    sys.stdout.write("Frame %d of %d\r" % (frame_id+1,traj.size/stride))
    
    if first_time:
        s1, s2  = frame.asAtoms.filter(is_nucleic), ref.asAtoms.filter(is_nucleic)
        moved_atoms = frame.asAtoms.filter((aName == "N") | (aName == HN_mask) | (aName == "P"))
        assert len(s1) == len(s2)
        N, H = [], []
        for resid in resids:
            N.append(frame.asAtoms.filter((rId == resid) & (aName == "N"))[0])
            H.append(frame.asAtoms.filter((rId == resid) & (aName == HN_mask))[0])
            vectors.append(VectorXYZ([(H[-1].r - N[-1].r)]))
        first_time = False

    else:
        # align frame by ref
        alignment = calc_alignment(s2.toCoords, s1.toCoords)
        moved_atoms.transform(alignment)
        # for a in moved_atoms:
        #     a.r = alignment.transform(a.r)
        for i, N_at, H_at in zip(range(len(resids)), N, H):
            vectors[i].append(H_at.r - N_at.r)
    frame_id += 1
sys.stdout.write("\n")

# # save vectors
# print("Saving vectors...")
# pickle.dump([ x.to_np for x in vectors ], open( "vectors.pickle", "wb" ))
# print("vectors.pickle created")

# calculate autocorrelation functions
print("Calculating autocorrelation functions...")
for i, rid, rname in zip(range(len(resids)), resids, resnames):
    sys.stdout.write("Residue %d of %d\r" % (i+1,len(resids)))
    ac = np.array(calc_autocorr_order_2(vectors[i], limit=int(cut_autocorr_function*1000/stride)))
    steps = np.arange(int(cut_autocorr_function*1000/stride))
    np.savetxt("cor_NH_n%d_s%d"%(n_steps,stride)+"/%04d_%s.cor" % (rid, rname), np.vstack((steps,ac)).T,fmt="%14.6e")
sys.stdout.write("\n")
print("Done!")








    