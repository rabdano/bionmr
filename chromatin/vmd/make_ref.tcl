mol new box.pdb type pdb first 0 last -1 step 1 waitfor all

set A [atomselect top "residue 0 to 134"]
set B [atomselect top "residue 135 to 236"]
set C [atomselect top "residue 237 to 364"]
set D [atomselect top "residue 365 to 486"]
set E [atomselect top "residue 487 to 621"]
set F [atomselect top "residue 622 to 723"]
set G [atomselect top "residue 724 to 851"]
set H [atomselect top "residue 852 to 973"]
set I [atomselect top "residue 974 to 1120"]
set J [atomselect top "residue 1121 to 1267"]

$A set chain A
$B set chain B
$C set chain C
$D set chain D
$E set chain E
$F set chain F
$G set chain G
$H set chain H
$I set chain I
$J set chain J

set sel [atomselect top "not (chain A B C D E F G H I J)"]
$sel set chain X

set sel [atomselect top all]
$sel writepdb ref.pdb
exit