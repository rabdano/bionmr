from subprocess import call
from glob import glob
import os
import numpy as np

hydropro = "/home/seva/bin/HYDROPRO/hydropro10/hydropro10-lnx.exe"
pdb_glob = "../../pdb_unwrapped/run?????.pdb"
first = 51
last = 1000
stride = 1

input_template = """%s              !Name of molecule
%s              !Name for output file 
%s              !Strucutural (PDB) file       
2               !Type of calculation
4.8,            !AER, radius of primary elements
-1,             !NSIG
25.,            !T (temperature, centigrade)
0.008475,        !ETA (Viscosity), Malvern Zetasizer Solvent Builder: 5 mM Tris buffer with pH 7, 100 mM NaCl, 0.5 mM EDTA and 0.1 mM MgCl2 @ 25C
198182.7,       !RM (Molecular weight)
0.702,          !Partial specific volume, cm3/g
1.0,            !Solvent density, g/cm3
0               !Number of values of Q
-1              !Number of intervals
0,              !Number of trials for MC calculation of covolume
1               !IDIF=1 (yes) for full diffusion tensors
"""

f = open("hydropro.dat", "w+")
for file in sorted(glob(pdb_glob))[first-1:last:stride]:
    bn = os.path.basename(file)
    cp_cmd = ["cp", file, bn]
    print(cp_cmd)
    call(cp_cmd)
    mol_name = bn[:-4]
    out_name = bn[:-4]
    pdb_file = bn

    inp_file = input_template % (mol_name, out_name, pdb_file)

    f.write(inp_file)

f.write("*               !End of file")
f.close()

call([hydropro])

results = sorted(glob("run?????-res.txt"))
taus = np.zeros(len(results), dtype="float")

for i, res in enumerate(results):
    with open(res, "r") as f:
        for line in f:
            if "Harmonic mean (correlation) time" in line:
                tau = float(line.split(":")[1].strip().split()[0])
                taus[i] = tau
                print(res,tau*1e9,"ns")
                break

print("Mean tau =", "%.1f" % (np.mean(taus) * 1e9), "ns")
print("STDEV =", "%.1f" % (np.std(taus) * 1e9), "ns")

with open('results.txt', 'w') as f:
    f.write("Mean tau = %.1f ns\n" % (np.mean(taus) * 1e9))
    f.write(("STDEV = %.1f ns\n" % (np.std(taus) * 1e9)))

call(['rm', '-rf', './*.pdb'])
