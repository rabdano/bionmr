import sys
from parmed.amber import *

# Run parameters
print('Run parameters:')
print('=' * 80)
try:
    f = open('5_run/run00001.in')
    for line in f:
        sys.stdout.write(line)
    sys.stdout.write('\n')
    f.close()
except:
    print('There are no .in files in 5_run.')
print('=' * 80)
print('\n' * 2)


# Number of molecules
parm = AmberParm('1_build/box.mod.prmtop')
residue_names = [res.name for res in parm.residues]

sodium_names = {'Na+', 'SOD'}
chlorine_names = {'Cl-', 'CLA'}
water_names = {'T3P', 'WAT', 'HOH', 'SOL', 'SPC', 'T4D', 'T4E'}

n_sodium = 0
for sn in sodium_names:
    for res_name in residue_names:
        if sn in res_name:
            n_sodium += 1

n_chlorine = 0
for cl in chlorine_names:
    for res_name in residue_names:
        if cl in res_name:
            n_chlorine += 1


n_water = 0
for wat in water_names:
    for res_name in residue_names:
        if wat in res_name:
            n_water += 1

print('Number of molecules, residues and atoms:')
print('=' * 80)
print('Atoms: %d' % len(parm.atoms))
print('Residues: %d' % len(residue_names))
print('Water: %d' % n_water)
print('Na+: %d' % n_sodium)
print('Cl-: %d' % n_chlorine)
print('=' * 80)
print('\n' * 2)


# Concentrations
print('Concentrations:')
print('=' * 80)

Avogadro_number = 6.022140857e23
volume_of_H2O_molecule = 2.99150757642e-26  # L
water_volume = n_water * volume_of_H2O_molecule

histones_concentration = 2 / Avogadro_number / water_volume
DNA_concentration = 1 / Avogadro_number / water_volume
sodium_conentration = n_sodium / Avogadro_number / water_volume
chlorine_conentration = n_chlorine / Avogadro_number / water_volume
print('Histones concentration (monomers): %.3f mM' % (histones_concentration * 1000))
print('DNA concentration: %.3f mM' % (DNA_concentration * 1000))
print('Na+ concentration: %.3f mM' % (sodium_conentration * 1000))
print('Cl- concentration: %.3f mM' % (chlorine_conentration * 1000))
# print('Na+ concentration: %f mM' % (55.5 * n_sodium / n_water * 18.01528 / 22.989769))
# print('Cl- concentration: %f mM' % (55.5 * n_chlorine / n_water * 18.01528 / 35.453))
print('=' * 80)
print('\n' * 2)