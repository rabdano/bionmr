from bionmr_utils.md import *
from get_secondary_structure_residues import *
from tqdm import tqdm
import numpy as np
import json


#
#  setup trajectory parameters
#
with open("input.json", "r") as f:
    pars = json.load(f)
path_to_traj = pars["path_to_traj"]
first_dat_file = int(pars["first_dat_file"])
last_dat_file = int(pars["last_dat_file"])
stride = int(pars["stride"])
reference_pdb = pars["reference_pdb"]
rmsd_fnout = pars["rmsd_fnout"]
ig_like_domain_start = pars["ig_like_domain_start"]


#
#  read trajectory
#
traj, ref = traj_from_dir(path=path_to_traj,
                          reference_pdb=path_to_traj + "/0_prepare/ref.pdb",
                          stride=1,
                          first=first_dat_file,
                          last=last_dat_file,
                          subdir="6_run",
                          filetype="nc")

print("Trajectory contain %d frames with %d chains / %d residues / %d atoms" % (len(traj),
                                                                                len(traj[0].asChains),
                                                                                len(traj[0].asResidues),
                                                                                len(traj[0].asAtoms)))
print("Using run%05d.nc - run%05d.nc" % (first_dat_file, last_dat_file))


#
#  get secondary structure information:
#  get residue ids for secondary structure from DSSP using Biopython and PDB
#
protein_chains = "A"
ss_residues_by_chain = []
incr = 0
for c in protein_chains:
    ss = get_secondary_structure_residues(c, pdb_code="5NGJ")
    ss_residues_by_chain.append(set([(x+incr) for x in ss]))
    incr += len(ref.asResidues.filter(cName == c))
# print("Secondary structure residues by chain")
# print(ss_residues_by_chain)
ss_residues = []
for chain_r in ss_residues_by_chain:
    for r in chain_r:
        ss_residues.append(r)
print("Secondary structure residues combined")
print(ss_residues)
ss_residues = set(ss_residues)


#
#  set predicates
#
sec_str_ca = (aName == "CA") & (rId.is_in(ss_residues))
sec_str_ca_tube = (aName == "CA") & (rId.is_in(ss_residues)) & (rId < ig_like_domain_start)
sec_str_ca_ig_like = (aName == "CA") & (rId.is_in(ss_residues)) & (rId >= ig_like_domain_start)


#
#  reference atoms for alignment
#
reference = PdbFile(reference_pdb).get_frame()
ref_ats = reference.asAtoms


#
#  create data variables
#
rmsd_ca_ss = np.zeros(int(traj.size / stride))
rmsd_ca_ss_tube = np.zeros(int(traj.size / stride))
rmsd_ca_ss_ig_like = np.zeros(int(traj.size / stride))


#
#  run through trajectory
#
for frame in tqdm(traj[::stride], disable=False):
    if frame.index == 0:
        # create variables for selections from atoms of frame
        frame_ats = frame.asAtoms

        ref_align_ca_ss = ref_ats.filter(sec_str_ca)
        ref_align_ca_ss_tube = ref_ats.filter(sec_str_ca_tube)
        ref_align_ca_ss_ig_like = ref_ats.filter(sec_str_ca_ig_like)

        frame_align_ca_ss = frame_ats.filter(sec_str_ca)
        frame_align_ca_ss_tube = frame_ats.filter(sec_str_ca_tube)
        frame_align_ca_ss_ig_like = frame_ats.filter(sec_str_ca_ig_like)

    # rmsd for all sec.str. CA atoms
    rmsd_ca_ss[int(frame.index / stride)] = calc_rmsd(ref_align_ca_ss.toCoords,
                                                      frame_align_ca_ss.toCoords,
                                                      frame_align_ca_ss.alignment_to(ref_align_ca_ss))

    # rmsd for sec.str. CA atoms of tube body domain
    rmsd_ca_ss_tube[int(frame.index / stride)] = calc_rmsd(ref_align_ca_ss_tube.toCoords,
                                                           frame_align_ca_ss_tube.toCoords,
                                                           frame_align_ca_ss_tube.alignment_to(ref_align_ca_ss_tube))

    # rmsd for sec.str. CA atoms of Ig-like domain
    rmsd_ca_ss_ig_like[int(frame.index / stride)] = calc_rmsd(ref_align_ca_ss_ig_like.toCoords,
                                                              frame_align_ca_ss_ig_like.toCoords,
                                                              frame_align_ca_ss_ig_like.alignment_to(ref_align_ca_ss_ig_like))

#
#  write RMSD to file
#
header = ["time [ps]",
          r"sec.str. C$\rm\alpha$ [$\rm\AA$]",
          r"sec.str. C$\rm\alpha$ tube body [$\rm\AA$]",
          r"sec.str. C$\rm\alpha$ ig-like [$\rm\AA$]"]
time = np.linspace((first_dat_file - 1) * 1000 + stride, last_dat_file * 1000, int(traj.size / stride))
rmsd = np.vstack((time,
                  rmsd_ca_ss,
                  rmsd_ca_ss_tube,
                  rmsd_ca_ss_ig_like)).T
np.savetxt(
    fname=rmsd_fnout,
    X=rmsd,
    fmt="%.5f",
    header=",".join(header),
    delimiter=","
)
