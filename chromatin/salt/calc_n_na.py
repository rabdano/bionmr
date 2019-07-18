from subprocess import call
from glob import glob
import os

path_to_pdb = '../../pdb'
path_to_pdb_unwrapped = '../../pdb_unwrapped'
tmp_file = 'tmp.pdb'
n_na_file = 'n_na.txt'
vmd_script_file = 'vmd_script.tcl'

pdbs = sorted(glob(path_to_pdb + '/run?????.pdb'))
pdbs_unwrapped = sorted(glob(path_to_pdb_unwrapped + '/run?????.pdb'))

assert len(pdbs) == len(pdbs_unwrapped), 'Different number of pdbs and pdbs_unwrapped.'

vmd_script = """set filename %s
mol new $filename type {pdb} first 0 last -1 step 1 waitfor all
set sel [atomselect top \"name 'Na+' and within 5 of (nucleic and index 348879 to 358224)\"]
set n_na [$sel num]
set fo [open %s a]
puts $fo "%s $n_na"
close $fo
exit
"""


cryst_records = []

for fn in pdbs_unwrapped:
    f = open(fn, 'r')
    cryst_records.append(f.readline())
    f.close()

call(['rm', '-rf', n_na_file])

for i, fn in enumerate(pdbs):
    f = open(fn, 'r')

    # save temporary pdb file
    nf = open(tmp_file, 'w')
    nf.write(cryst_records[i])
    for line in f:
        nf.write(line)

    # run PropPDB
    cmd = ['/opt/amber16'+'/bin/PropPDB', '-p', tmp_file]
    cmd += ['-o', 'x8_' + tmp_file, '-ix', '3', '-iy', '3', '-iz', '3']
    call(cmd)

    script = vmd_script % ('x8_' + tmp_file, n_na_file, os.path.basename(fn))
    with open(vmd_script_file, 'w') as f:
        f.write(script)

    cmd = ['vmd', '-e', vmd_script_file]

    call(cmd)
    call(['rm', '-rf', tmp_file, 'x8_' + tmp_file, vmd_script_file])
