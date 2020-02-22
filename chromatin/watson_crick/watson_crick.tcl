# this script is memory-consuming because of VMD memory leaks
# use system with at least 32 GB of RAM

set n_bp 145
set n_prot_res 974
set n_res 1264
set n_frames [molinfo top get numframes]

set outfile [open "report_full.txt" w]
set outfile2 [open "report.txt" w]

for { set frame 0 } { $frame < $n_frames } { incr frame } {
    animate goto $frame
    if { $frame == 0 } {
        for { set bp 0 } { $bp < $n_bp } { incr bp } {
            set res_id [expr $bp + $n_prot_res]
            set res_cmp_id [expr $n_res - $bp - 1]
            set res [atomselect top "residue $res_id"]
            set res_cmp [atomselect top "residue $res_cmp_id"]
            set res_name [lindex [$res get resname] 0]
            set res_cmp_name [lindex [$res_cmp get resname] 0]
            puts -nonewline $outfile "$res_name$res_id-$res_cmp_name$res_cmp_id,"
        }
    }
    puts -nonewline $outfile "\n"

    set d_cum 0.0

    for { set bp 0 } { $bp < $n_bp } { incr bp } {
        set res_id [expr $bp + $n_prot_res]
        set res_cmp_id [expr $n_res - $bp - 1]
        set res [atomselect top "residue $res_id"]
        set res_cmp [atomselect top "residue $res_cmp_id"]
        set res_name [lindex [$res get resname] 0]
        set res_cmp_name [lindex [$res_cmp get resname] 0]

        if { $res_name == "DA" || $res_name == "DA3" || $res_name == "DA5" } {
            set N1 [atomselect top "residue $res_id and name N1"]
            set N6 [atomselect top "residue $res_id and name N6"]
            set N3 [atomselect top "residue $res_cmp_id and name N3"]
            set O4 [atomselect top "residue $res_cmp_id and name O4"]
            set N1_idx [lindex [$N1 get index] 0]
            set N6_idx [lindex [$N6 get index] 0]
            set N3_idx [lindex [$N3 get index] 0]
            set O4_idx [lindex [$O4 get index] 0]
            set d1 [measure bond "$N1_idx $N3_idx"]
            set d2 [measure bond "$N6_idx $O4_idx"]
            set mean_d [expr ($d1 + $d2) / 2]
            $N1 delete
            $N6 delete
            $N3 delete
            $O4 delete

        } elseif { $res_name == "DT" || $res_name == "DT3" || $res_name == "DT5" } {
            set N1 [atomselect top "residue $res_cmp_id and name N1"]
            set N6 [atomselect top "residue $res_cmp_id and name N6"]
            set N3 [atomselect top "residue $res_id and name N3"]
            set O4 [atomselect top "residue $res_id and name O4"]
            set N1_idx [lindex [$N1 get index] 0]
            set N6_idx [lindex [$N6 get index] 0]
            set N3_idx [lindex [$N3 get index] 0]
            set O4_idx [lindex [$O4 get index] 0]
            set d1 [measure bond "$N1_idx $N3_idx"]
            set d2 [measure bond "$N6_idx $O4_idx"]
            set mean_d [expr ($d1 + $d2) / 2]
            $N1 delete
            $N6 delete
            $N3 delete
            $O4 delete
        } elseif { $res_name == "DG" || $res_name == "DG3" || $res_name == "DG5" } {
            set O6 [atomselect top "residue $res_id and name O6"]
            set N1 [atomselect top "residue $res_id and name N1"]
            set N2 [atomselect top "residue $res_id and name N2"]
            set N4 [atomselect top "residue $res_cmp_id and name N4"]
            set N3 [atomselect top "residue $res_cmp_id and name N3"]
            set O2 [atomselect top "residue $res_cmp_id and name O2"]
            set O6_idx [lindex [$O6 get index] 0]
            set N1_idx [lindex [$N1 get index] 0]
            set N2_idx [lindex [$N2 get index] 0]
            set N4_idx [lindex [$N4 get index] 0]
            set N3_idx [lindex [$N3 get index] 0]
            set O2_idx [lindex [$O2 get index] 0]
            set d1 [measure bond "$O6_idx $N4_idx"]
            set d2 [measure bond "$N1_idx $N3_idx"]
            set d3 [measure bond "$N2_idx $O2_idx"]
            set mean_d [expr ($d1 + $d2 + $d3) / 3]
            $O6 delete
            $N1 delete
            $N2 delete
            $N4 delete
            $N3 delete
            $O2 delete
        } elseif { $res_name == "DC" || $res_name == "DC3" || $res_name == "DC5" } {
            set O6 [atomselect top "residue $res_cmp_id and name O6"]
            set N1 [atomselect top "residue $res_cmp_id and name N1"]
            set N2 [atomselect top "residue $res_cmp_id and name N2"]
            set N4 [atomselect top "residue $res_id and name N4"]
            set N3 [atomselect top "residue $res_id and name N3"]
            set O2 [atomselect top "residue $res_id and name O2"]
            set O6_idx [lindex [$O6 get index] 0]
            set N1_idx [lindex [$N1 get index] 0]
            set N2_idx [lindex [$N2 get index] 0]
            set N4_idx [lindex [$N4 get index] 0]
            set N3_idx [lindex [$N3 get index] 0]
            set O2_idx [lindex [$O2 get index] 0]
            set d1 [measure bond "$O6_idx $N4_idx"]
            set d2 [measure bond "$N1_idx $N3_idx"]
            set d3 [measure bond "$N2_idx $O2_idx"]
            set mean_d [expr ($d1 + $d2 + $d3) / 3]
            $O6 delete
            $N1 delete
            $N2 delete
            $N4 delete
            $N3 delete
            $O2 delete
        }
        unset -nocomplain N1 N6 N3 O4 O6 N2 N4 O2 d1 d2 d3
        puts -nonewline $outfile "$mean_d,"
        set d_cum [expr $d_cum + $mean_d]
    }
    $res delete
    $res_cmp delete
    unset -nocomplain res res_cmp res_id res_cmp_id res_name res_cmp_name
    set mean_frame_d [expr $d_cum / $n_bp]
    puts $outfile2 "$frame - $mean_frame_d"
    puts "$frame - $mean_frame_d"
}

close $outfile
close $outfile2
