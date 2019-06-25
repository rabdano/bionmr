from multiprocessing import Pool, cpu_count, current_process
import subprocess
from glob import glob
import os

path_to_shifts = '/home/seva/bin/shifts-5.5/bin/shifts'
path_to_pdb = 'pdb_shifts'


def shifts(file):
    """
    $ /home/seva/bin/shifts-5.5/bin/shifts -h
    shifts [Options] [<expression>] <basename>

    Arguments
       <basename>
                  The base name of the system you wish to analyze.
                  The <basename>.pdb file is used as the starting
                  structure. <basename>.obs must be present if the
                  -readobs flag (below) is used. The predicted shifts
                  are printed to <basename>.emp.
       <expression>
                  The NAB atom expression used to select atoms to
                  compute chemical shifts for. "::H*" will select
                  all hydrogens.

    Options:
       -csa       Compute chemical shift anisotropy tensors based on
                  ring and peptide group susceptibilities.*
       -nocoil    Do not attempt to compare computed shifts to 'random
                  coil' reference values. Only valid without -qdb
       -reslib    Use an Amber residue library; good for nucleic acids;
                  needed for electrostatic effects in non-proteins
       -noreduce  Don't run reduce: assume all protons in input file are OK
       -readobs   Read observed shifts (<basename>.obs) and compare
                  with predicted shifts.
       -details   Print out individual contributions to each shift.
       -sander    Write out a sander input file based on calculated
                  shifts.
       -qdb       Use a 'quantum database' approach for 13C and 15N
                  shifts in proteins.
       -refine    Try to find side-chain angles that improve qdb fit
       -HN        Use a module that predicts amide proton shifts.

    *The CSA approach is not yet fully parametrized.

    This program computes chemical shifts for proteins and nucleic
    acids based on empirical formulas. See the shifts manual for
    more details and citations.
    """
    curr_proc = current_process()
    wd = curr_proc.name
    basename = os.path.basename(file)[:-4]
    subprocess.call(["mkdir", "-p", wd])
    subprocess.call(["cp " + file[:-4] + '*.pdb ' + wd], shell=True)

    # PROTONS
    # calculate for protein + dna
    shifts_cmd = [path_to_shifts]
    shifts_cmd += ['-noreduce']
    shifts_cmd += ['2,6::H,HA']
    shifts_cmd += [basename]
    subprocess.call(shifts_cmd, cwd=wd)

    # calculate for protein only
    shifts_cmd = [path_to_shifts]
    shifts_cmd += ['-noreduce']
    shifts_cmd += ['2,6::H,HA']
    shifts_cmd += [basename + '_p']
    subprocess.call(shifts_cmd, cwd=wd)

    # CARBON and NITROGEN
    # calculate for protein + dna
    shifts_cmd = [path_to_shifts]
    shifts_cmd += ['-noreduce']
    shifts_cmd += ['-qdb']
    shifts_cmd += [basename + '_cn']
    subprocess.call(shifts_cmd, cwd=wd)

    # calculate for protein only
    shifts_cmd = [path_to_shifts]
    shifts_cmd += ['-noreduce']
    shifts_cmd += ['-qdb']
    shifts_cmd += [basename + '_p_cn']
    subprocess.call(shifts_cmd, cwd=wd)

    subprocess.call(["cp " + wd + "/" + "*.emp " + path_to_pdb], shell=True)
    subprocess.call(["cp " + wd + "/" + "*.qdb " + path_to_pdb], shell=True)
    subprocess.call(["rm", "-rf", wd])


files = sorted(glob(path_to_pdb + "/run?????.pdb"))

ncpu = cpu_count()
pool = Pool(ncpu)

pool.map(shifts, files)
