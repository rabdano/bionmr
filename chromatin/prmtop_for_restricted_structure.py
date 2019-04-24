from parmed.amber import *

p = AmberParm('1_build/box.mod.prmtop')
p.strip('!:1-1268')
p.remake_parm()

p.write_parm('1_build/restricted_structure.parm7')
