# Save possible contacts

set cut 3.5
set fp [open "contacts.txt" w+]
set n_frames [molinfo top get numframes]
set sel_string "((hydrophobic or aromatic) and sidechain) and (exwithin $cut of resname CML) and noh"


for { set frame 0 } { $frame < $n_frames } { incr frame } {
    animate goto $frame
    set c [atomselect top $sel_string]
    puts $fp "$frame [lsort -unique [$c get resid]]"
}

close $fp
