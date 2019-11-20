# split chains for pb6 tube made from 9 rings
set n 464
set tot_n [expr $n * 9]

set chain_names "ABCDEFGHI"

set i 0
for {set res 0} {$res < $tot_n} {incr res $n} {
    set r1 $res
    set r2 [expr $res + $n - 1]
    set sel [atomselect top "residue $r1 to $r2"]
    set name [string range $chain_names $i $i]
    $sel set chain $name
    set chain_resid 1
    for {set r $r1} {$r <= $r2} {incr r} {
        set r_ats [atomselect top "residue $r"]
        $r_ats set resid $chain_resid
        incr chain_resid
    }
    unset sel
    incr i
}
