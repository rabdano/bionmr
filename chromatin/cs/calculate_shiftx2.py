from multiprocessing import Pool, cpu_count, current_process
import subprocess
from glob import glob
import os


path_to_pdb = 'pdb_shiftx2'

def shiftx2(file):
    """
    ## SHIFTX2 Ver 1.10A (Released Sep 12, 2016) ##
    python /home/seva/bin/shiftx2-1.13/shiftx2-linux/shiftx2.py -h
    Usage: usage shiftx2.py -i infile.pdb [options]

    Options:
      -h, --help            show this help message and exit
      -i INFILE, --infile=INFILE
                            Input file
      -o OUTFILE, --outfile=OUTFILE
                            Output file
      -c CHAINID, --chain=CHAINID
                            the chain ID to predict, by defaut it is the first
                            chain to be predicted
      -f OUTFORMAT, --outformat=OUTFORMAT
                            [CSV|BMRB|TABULAR|NEF]
      -a ATOMS, --atoms=ATOMS
                            [ALL|BACKBONE|SIDECHAIN]
      -p PH, --ph=PH        PH value, default PH=5.0
      -t TEMP, --temperature=TEMP
                            Temperature, default Temp=298K
      -b BATCH, --batch=BATCH
                            run shiftx2 in batch mode for PDBs matching the given
                            regular expression. For example: -b '*.pdb' predicts
                            all PDB files in the current directory (chain
                            selection is disabled in batch mode)
      -m, --nmr             run shiftx2 in NMR mode to calculate the averaged
                            shifts for multiple models
      -g OBSCS, --graph=OBSCS
                            given an observed chemical shift file in BMRB format,
                            generate correlation graphs for multiple NMR model
                            predictions. this option can only be used in NMR mode.
      -u, --multichain      toggle option to predict multi-chain PDB as a single
                            long chain. SHIFTX2 will concatenate all chains in a
                            single input PDB file and make prediction for the
                            concatenated chain.
      -d, --deuterated      toggle option for deuterated proteins
      -v, --verbose         toggle option to be verbose and print debug output
                            messages
      -e, --explicit        toggle option to keep intermediate SHIFTX+ (.sxp
                            files) and SHIFTY+ (.shifty files) results
      -x, --shiftx1         toggle option to run SHIFTX 1.0 instead of SHIFTX+
      -n, --noshifty        toggle option to exclude SHIFTY+ results
      -z TEMPFOLDER, --temp=TEMPFOLDER
                            Specify which folder to store temporary files
      -r, --nmr-mc          run shiftx2 in NMR-Multichain mode for PDB with
                            multiple models, each of which contains multiple
                            chains.
      -k, --phosphorylated  toggle option to calculate phosphorylated amino
                            acids's chemical shifts
    """
    curr_proc = current_process()
    wd = curr_proc.name
    fn = os.path.basename(file)
    subprocess.call(["mkdir", "-p", wd])
    subprocess.call(["cp", file, wd])
    
    shiftx2_cmd = ["python", "/home/seva/bin/shiftx2-1.13/shiftx2-linux/shiftx2.py"]
    shiftx2_cmd += ["-i", wd + "/" + fn]
    # shiftx2_cmd += ["-o", file[:-4]+".cs"]
    shiftx2_cmd += ["-f", "CSV"]
    shiftx2_cmd += ["-a", "BACKBONE"]
    shiftx2_cmd += ["-p", "7.0"]
    shiftx2_cmd += ["-u"]
    shiftx2_cmd += ["-z", curr_proc.name]
    print(" ".join(shiftx2_cmd))
    subprocess.call(shiftx2_cmd)
    
    subprocess.call(["cp", wd + "/" + fn[:-4] + ".pdb.cs", file[:-4]+".cs"])
    subprocess.call(["rm", "-rf", wd])

files = sorted(glob(path_to_pdb + "/run?????.pdb"))

ncpu = cpu_count()
pool = Pool(ncpu)

pool.map(shiftx2, files)

