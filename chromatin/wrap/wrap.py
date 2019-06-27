from bionmr_utils.md import *
from glob import glob
import numpy as np
import errno
import os
import sys

altered_records = AlteredPdbRecords(StandardPdbRecords.instance())
altered_records.alter_record(PdbRecordName("ATOM"), PdbFieldName("serial"), [7, 12])


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def split_chains(frame: Frame, cids="ABCDEFGHIJ", c_ranges=None):
    if not c_ranges:
        c_ranges = [set(range(1, 136)),
                    set(range(136, 238)),
                    set(range(238, 366)),
                    set(range(366, 488)),
                    set(range(488, 623)),
                    set(range(623, 725)),
                    set(range(725, 853)),
                    set(range(853, 975)),
                    set(range(975, 1122)),
                    set(range(1122, 1269))]
    frm = Frame(0)
    for i, cid in enumerate(cids):
        chain = frm.emplace(ChainName(cid))
        for res in frame.asResidues.filter(rId.is_in(c_ranges[i])):
            chain.emplace(res)
    return frm


# create output dir
output_dir = '../pdb_unwrapped/'
mkdir_p('../pdb_unwrapped/')

# get reference structure with good arrangement of chains
ref = PdbFile('../1_build/box.pdb', altered_records).get_frame()
ref_frm = split_chains(ref)
ref_ats = ref_frm.asAtoms

chain_ids = [chain.name.str for chain in ref_frm.asChains]

for i, chain in enumerate(ref_ats.asChains):
    print(i, chain.cName, len(chain.asResidues))

# get PBC from XST file and inpcrd file
pbc = np.loadtxt('../5_run/traj.xst')
v1 = VectorXYZ.from_numpy(pbc[:, 1:4])
v2 = VectorXYZ.from_numpy(pbc[:, 4:7])
v3 = VectorXYZ.from_numpy(pbc[:, 7:10])
scaling_factors = [v1[i].len()/v1[0].len() for i in range(len(v1))]
lat_vec = LatticeVectors(v1[0], v2[0], v3[0])
shift_finder = BestShiftFinder(lat_vec)
inpcrd = '../1_build/box.inpcrd'
with open(inpcrd) as f:
    content = f.readlines()
content = content[-1].split()[:6]
a, b, c, alpha, beta, gamma = [float(x) for x in content]

# get frames from pdb files
pdbs = sorted(glob('run?????.pdb'))

# probe atoms
probe = ((cName == 'J') & (aName == 'P'))

for i, pdb in enumerate(pdbs):
    # get frame atoms
    frame = PdbFile(pdb, altered_records).get_frame()
    frame = split_chains(frame)
    frame_ats = frame.asAtoms

    # align reference by first frame nucleic P
    if i == 0:
        ref_ats.transform(ref_ats.filter(probe).alignment_to(frame_ats.filter(probe)))

        with open(output_dir + 'ref.pdb', 'w') as f:
            fmt = 'CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%10s%4s\n'
            f.write(fmt % (a, b, c, alpha, beta, gamma, ' ' * 10, ' ' * 4))
            ref_ats[0].frame.to_pdb(f)

    # align all atoms by one of DNA chains
    ref_ats.transform(ref_ats.filter(probe).alignment_to(frame_ats.filter(probe)))

    # cycle through chains and check first atom displacement
    shift_finder.scale_lattice_by(scaling_factors[i])

    sys.stdout.write(pdb + '\t:\t')
    for cid in sorted(chain_ids):
        at_ref = ref_ats.filter(cName == cid)
        at_frame = frame_ats.filter(cName == cid)
        sys.stdout.write(cid + ' ')
        sys.stdout.flush()

        dif, shift = shift_finder.find_best_shift(at_ref.filter(cName == cid).geom_center(),
                                                  at_frame.filter(cName == cid).geom_center())

        at_frame.transform(Translation3d(shift))

    sys.stdout.write('\n')
    shift_finder.scale_lattice_by(1.0/scaling_factors[i])

    # write updated pdb
    with open(output_dir + pdb, 'w') as f:
        a = v1[i].len()
        b = v2[i].len()
        c = v3[i].len()
        fmt = 'CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%10s%4s\n'
        f.write(fmt % (a, b, c, alpha, beta, gamma, ' ' * 10, ' ' * 4))
        frame_ats[0].frame.to_pdb(f)
