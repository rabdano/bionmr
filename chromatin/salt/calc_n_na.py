from subprocess import call, Popen, PIPE
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
set sel [atomselect top \"name 'Na+' and within 5 of (nucleic and index %d to %d)\"]
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

    # get indexes for DNA in 14th copy
    process = Popen(['grep', "HO5' DA5 I", 'x8_' + tmp_file], stdout=PIPE)
    out, err = process.communicate()
    lines = out.split('\n')
    line = lines[13]
    index1 = int(line.split()[1]) - 1

    process = Popen(['grep', "HO3' DT3 J", 'x8_' + tmp_file], stdout=PIPE)
    out, err = process.communicate()
    lines = out.split('\n')
    line = lines[13]
    index2 = int(line.split()[1]) + 1

    print(index1, index2)

    script = vmd_script % ('x8_' + tmp_file, index1, index2, n_na_file, os.path.basename(fn))
    with open(vmd_script_file, 'w') as f:
        f.write(script)

    cmd = ['vmd', '-e', vmd_script_file]

    call(cmd)
    call(['rm', '-rf', tmp_file, 'x8_' + tmp_file, vmd_script_file])
