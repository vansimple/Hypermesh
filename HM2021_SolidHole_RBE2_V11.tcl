# HM2021_SolidHole_RBE2_V11.tcl -- V11 (direct cross-component pairs; no temporary spiders) -- HyperMesh 2021, OptiStruct / Nastran
# ASCII source avoids Windows Tcl encoding problems. All lengths use model units.
# Usage: source {C:/your/path/HM2021_SolidHole_RBE2_V11.tcl}
# Select COMPONENTS, not partial solid elements. Hidden elements in selected
# components are included, so the extracted solid skin is complete.
# Supported: linear TET4 / PYRAMID5 / PENTA6 / HEX8, conforming mesh.
# Higher-order solids cause that component to be skipped (not truncated).
# Closed near-planar sharp-edged profiles: round, ellipse, polygon; >=4 nodes.
# Diameter filter uses equivalent area diameter D=2*sqrt(projected area/pi).
# Supports near-prismatic walls; strong taper, fillets or branched rims may fail.
# V4: one vertex-adjacent solid layer around the entire cylindrical wall.
# Include inner/outer nodes at all depths; NOT only the cylindrical skin.
# Exclude blind-floor interior nodes; retain the wall/floor shared perimeter.
# Expand from wall nodes to incident solids ONCE, never recursively.
# Filter interior disk nodes and nodes beyond the axial extent of the hole wall.
# A valid mouth must adjoin outward annular faces. Inward disk faces are rejected.
# Through-hole mouths merge by connected cylindrical wall; blind floor is NOT a mouth.
# Master position = arithmetic mean of final unique dependent-node coordinates.
# Master nodes remain required by RBE2; only their temporary display is cleared.
# Components split through a hole, chamfers, fillets and intersecting holes
# are outside this conservative recognizer's scope. Inspect the solver deck.
# Existing rigid/RBE3 node overlap skips the new spider. Remove old V1/V2/V3/V4 RBE2
# ELEMENTS before regeneration; do not delete original structural mesh nodes.
# Save your .hm model first. Marks 1/2 are used. Original mesh is never deleted.

namespace eval ::SHR11 {
    variable sourcePath [info script]
    variable diagnosticOnly 0
    variable verbose 0
    variable pairmode 1
    variable alignOnly 1
    variable axisTol 0.2
    variable angleTol 5.0
    variable diamTol 1.0
    variable gapTol 2.0
    variable dmin 3.0
    variable dmax 30.0
    variable tol 0.04
    variable angle 55.0
    variable dofs 123456
    variable through 0
    variable washer 1
    variable status "Set diameter range, then select components."
    variable busy 0
    variable made 0
    variable skipped 0
    variable xyz; variable en; variable ne
    variable fn; variable fc; variable normal; variable owner; variable ef
    variable reserved
}
proc ::SHR11::add {a b} {lassign $a x y z; lassign $b u v w; list [expr {$x+$u}] [expr {$y+$v}] [expr {$z+$w}]}
proc ::SHR11::sub {a b} {lassign $a x y z; lassign $b u v w; list [expr {$x-$u}] [expr {$y-$v}] [expr {$z-$w}]}
proc ::SHR11::mul {a s} {lassign $a x y z; list [expr {$x*$s}] [expr {$y*$s}] [expr {$z*$s}]}
proc ::SHR11::dot {a b} {lassign $a x y z; lassign $b u v w; expr {$x*$u+$y*$v+$z*$w}}
proc ::SHR11::cross {a b} {lassign $a x y z; lassign $b u v w; list [expr {$y*$w-$z*$v}] [expr {$z*$u-$x*$w}] [expr {$x*$v-$y*$u}]}
proc ::SHR11::norm {a} {expr {sqrt([dot $a $a])}}
proc ::SHR11::unit {a} {set s [norm $a]; if {$s<1.e-20} {error "Degenerate vector"}; mul $a [expr {1.0/$s}]}
proc ::SHR11::key {a b} {lsort -integer [list $a $b]}
proc ::SHR11::edges {ns} {
    set out {}; set prev [lindex $ns end]
    foreach n $ns {lappend out [key $prev $n]; set prev $n}; return $out
}
proc ::SHR11::mean {ns} {
    variable xyz; set c {0.0 0.0 0.0}
    foreach n $ns {set c [add $c $xyz($n)]}
    mul $c [expr {1.0/[llength $ns]}]
}
proc ::SHR11::templates {cfg} {
    switch -- $cfg {
        204 {return {{0 1 2} {0 3 1} {1 3 2} {2 3 0}}}
        205 {return {{0 1 2 3} {0 4 1} {1 4 2} {2 4 3} {3 4 0}}}
        206 {return {{0 1 2} {3 5 4} {0 3 4 1} {1 4 5 2} {2 5 3 0}}}
        208 {return {{0 1 2 3} {4 7 6 5} {0 4 5 1} {1 5 6 2} {2 6 7 3} {3 7 4 0}}}
    }
    return {}
}
# Scalar fields only. Query IDs with each batch so returned ordering is respected.
# Older/build-specific unsupported batch queries fall back to scalar reads.
proc ::SHR11::readFields {entity ids fields} {
    set result {}; set count [llength $ids]
    for {set start 0} {$start<$count} {incr start 2000} {
        set chunk [lrange $ids $start [expr {$start+1999}]]
        set rc [catch {
            *createmark $entity 1 {*}$chunk
            set order [hm_getvalue $entity mark=1 dataname=id]
            if {[lsort -integer $order] ne [lsort -integer $chunk]} {error "Batch ID mismatch"}
            set columns {}
            foreach field $fields {
                set values [hm_getvalue $entity mark=1 dataname=$field]
                if {[llength $values]!=[llength $order]} {error "Batch field length mismatch"}
                lappend columns $values
            }
        }]
        if {$rc} {
            foreach id $chunk {
                set row {}; foreach field $fields {lappend row [hm_getvalue $entity id=$id dataname=$field]}
                dict set result $id $row
            }
        } else {
            set i 0
            foreach id $order {
                set row {}; foreach col $columns {lappend row [lindex $col $i]}
                dict set result $id $row; incr i
            }
        }
    }
    return $result
}
proc ::SHR11::loadComponent {cid} {
    variable xyz; variable en; variable ne
    variable fn; variable fc; variable normal; variable owner; variable ef
    foreach a {xyz en ne fn fc normal owner ef} {array unset $a}
    *createmark elems 1 "by collector id" $cid
    set ids [hm_getmark elems 1]; set solids {}
    set cfgs [readFields elems $ids {config}]
    foreach e $ids {
        set cfg [lindex [dict get $cfgs $e] 0]
        if {$cfg in {210 213 215 220}} {error "Component $cid contains second-order solids; skipped."}
        if {[llength [templates $cfg]]} {lappend solids [list $e $cfg]}
    }
    array set faces {}; array set counts {}; array set fe {}
    foreach item $solids {
        lassign $item e cfg
        set ns [hm_getvalue elems id=$e dataname=nodes]
        if {[llength [lsort -unique $ns]] != [llength $ns]} {error "Degenerate connectivity at element $e"}
        set en($e) $ns
        foreach n $ns {
            lappend ne($n) $e
        }
        foreach ix [templates $cfg] {
            set f {}; foreach i $ix {lappend f [lindex $ns $i]}
            set k [lsort -integer $f]
            if {![info exists counts($k)]} {set counts($k) 0; set faces($k) $f; set fe($k) $e}
            incr counts($k)
        }
    }
    dict for {n point} [readFields nodes [array names ne] {x y z}] {set xyz($n) $point}
    array set centers {}
    set fid 0
    foreach k [array names counts] {
        if {$counts($k)>2} {error "Non-manifold solid face in component $cid"}
        if {$counts($k)!=1} {continue}
        set f $faces($k); set c [mean $f]; set e $fe($k)
        # Sum cross products for a more stable quad normal.
        set v {0 0 0}; set p $xyz([lindex $f end])
        foreach n $f {set q $xyz($n); set v [add $v [cross [sub $p $c] [sub $q $c]]]; set p $q}
        set v [unit $v]
        if {![info exists centers($e)]} {set centers($e) [mean $en($e)]}
        if {[dot $v [sub $c $centers($e)]]<0} {set v [mul $v -1]}
        incr fid; set fn($fid) $f; set fc($fid) $c; set normal($fid) $v; set owner($fid) $e
        foreach edge [edges $f] {lappend ef($edge) $fid}
    }
    return [llength $solids]
}
proc ::SHR11::loops {} {
    variable ef; variable normal; variable fc; variable angle
    set ct [expr {cos($angle*acos(-1)/180.0)}]; array set adj {}
    foreach edge [array names ef] {
        if {[llength $ef($edge)]!=2} {continue}
        lassign $ef($edge) a b
        if {[dot $normal($a) $normal($b)]>$ct} {continue}
        # Reject concave solid creases along polygonal hole walls; keep mouth edges.
        set delta [sub $fc($b) $fc($a)]
        set eps [expr {1.e-8*max(1.e-12,[norm $delta])}]
        if {[dot $normal($a) $delta]>$eps && [dot $normal($b) [mul $delta -1]]>$eps} {continue}
        lassign $edge n m; lappend adj($n) $m; lappend adj($m) $n
    }
    array set seen {}; set result {}
    foreach start [lsort -integer [array names adj]] {
        if {[info exists seen($start)]} {continue}
        set queue [list $start]; set seen($start) 1; set group {}; set good 1
        for {set i 0} {$i<[llength $queue]} {incr i} {
            set n [lindex $queue $i]; lappend group $n
            if {[llength $adj($n)]!=2} {set good 0}
            foreach m $adj($n) {if {![info exists seen($m)]} {set seen($m) 1; lappend queue $m}}
        }
        if {!$good || [llength $group]<4} {continue}
        set ordered {}; set prev -1; set n $start
        while {1} {
            lappend ordered $n; lassign $adj($n) a b
            set next [expr {$a==$prev ? $b : $a}]; set prev $n; set n $next
            if {$n==$start} {break}
            if {[llength $ordered]>[llength $group]} {set good 0; break}
        }
        if {$good} {lappend result $ordered}
    }
    return $result
}
# Historical procedure name retained; now returns a planar closed-profile model.
# Fields: area centroid, plane normal, equivalent radius, rim IDs, u, v, polygon2D.
proc ::SHR11::circle {ns} {
    variable xyz; variable tol; variable dmin; variable dmax
    set base [mean $ns]; set normal {0 0 0}
    set prev [sub $xyz([lindex $ns end]) $base]
    foreach n $ns {
        set cur [sub $xyz($n) $base]
        set normal [add $normal [cross $prev $cur]]; set prev $cur
    }
    if {[norm $normal]<1.e-20} {return {}}
    set axis [unit $normal]
    set u [unit [sub $xyz([lindex $ns 0]) $base]]
    set v [unit [cross $axis $u]]; set u [unit [cross $v $axis]]
    set poly {}
    foreach n $ns {set d [sub $xyz($n) $base]; lappend poly [list [dot $d $u] [dot $d $v]]}
    set area2 0.0; set cx 0.0; set cy 0.0; lassign [lindex $poly end] px py
    foreach pt $poly {
        lassign $pt x y; set q [expr {$px*$y-$x*$py}]
        set area2 [expr {$area2+$q}]
        set cx [expr {$cx+($px+$x)*$q}]; set cy [expr {$cy+($py+$y)*$q}]
        set px $x; set py $y
    }
    if {$area2<=1.e-20} {return {}}
    set r [expr {sqrt($area2/(2*acos(-1)))}]
    if {2*$r<$dmin*(1-1.e-10) || 2*$r>$dmax*(1+1.e-10)} {return {}}
    set cx [expr {$cx/(3*$area2)}]; set cy [expr {$cy/(3*$area2)}]
    set center [add $base [add [mul $u $cx] [mul $v $cy]]]
    set shifted {}; foreach pt $poly {lassign $pt x y; lappend shifted [list [expr {$x-$cx}] [expr {$y-$cy}]]}
    foreach n $ns {if {abs([dot [sub $xyz($n) $center] $axis])>$tol*$r} {return {}}}
    # Reject self-crossing profiles (proper crossings of nonadjacent edges).
    set count [llength $shifted]
    for {set i 0} {$i<$count} {incr i} {
        set a [lindex $shifted $i]; set b [lindex $shifted [expr {($i+1)%$count}]]
        for {set j [expr {$i+2}]} {$j<$count} {incr j} {
            if {$i==0 && $j==$count-1} {continue}
            set c [lindex $shifted $j]; set d [lindex $shifted [expr {($j+1)%$count}]]
            if {[orient2 $a $b $c]*[orient2 $a $b $d]<0 && [orient2 $c $d $a]*[orient2 $c $d $b]<0} {return {}}
        }
    }
    return [list $center $axis $r $ns $u $v $shifted]
}
proc ::SHR11::orient2 {a b c} {
    lassign $a ax ay; lassign $b bx by; lassign $c cx cy
    expr {($bx-$ax)*($cy-$ay)-($by-$ay)*($cx-$ax)}
}
# Returns signed distance to the actual projected polygon: negative inside.
# Also returns the inward unit normal of the nearest edge, in global coordinates.
# A one-profile cache bounds memory and eliminates repeated node/face projections.
proc ::SHR11::activateProfile {profile} {
    variable activeProfile; variable segments; variable pointCache
    variable wallCache; variable annulusCache
    if {[info exists activeProfile] && $activeProfile eq $profile} {return}
    set activeProfile $profile
    array unset pointCache; array unset wallCache; array unset annulusCache
    lassign $profile center axis r ns u v poly
    set segments {}; lassign [lindex $poly end] ax ay
    foreach pt $poly {
        lassign $pt bx by
        set dx [expr {$bx-$ax}]; set dy [expr {$by-$ay}]; set l2 [expr {$dx*$dx+$dy*$dy}]
        set inward {0 0 0}
        if {$l2>1.e-24} {set inward [mul [add [mul $u [expr {-$dy}]] [mul $v $dx]] [expr {1./sqrt($l2)}]]}
        lappend segments [list $ax $ay $bx $by $dx $dy $l2 $inward]
        set ax $bx; set ay $by
    }
}
proc ::SHR11::profileDistance {point profile} {
    variable segments; variable pointCache
    activateProfile $profile
    if {[info exists pointCache($point)]} {return $pointCache($point)}
    lassign $profile center axis r ns u v poly
    set q [sub $point $center]; set x [dot $q $u]; set y [dot $q $v]
    set inside 0; set best Inf; set inward {0 0 0}
    foreach seg $segments {
        lassign $seg ax ay bx by dx dy l2 edgeNormal
        if {$l2>1.e-24} {
            set t [expr {max(0.,min(1.,(($x-$ax)*$dx+($y-$ay)*$dy)/$l2))}]
            set ex [expr {$x-$ax-$t*$dx}]; set ey [expr {$y-$ay-$t*$dy}]
            set d2 [expr {$ex*$ex+$ey*$ey}]
            if {$d2<$best} {set best $d2; set inward $edgeNormal}
        }
        if {($ay>$y)!=($by>$y)} {
            if {$x<($bx-$ax)*($y-$ay)/($by-$ay)+$ax} {set inside [expr {!$inside}]}
        }
    }
    set distance [expr {sqrt($best)}]
    if {$inside} {set distance [expr {-$distance}]}
    if {[array size pointCache]>16384} {array unset pointCache}
    set pointCache($point) [list $distance $inward]
    return $pointCache($point)
}
proc ::SHR11::isWall {f profile} {
    variable wallCache
    activateProfile $profile
    if {![info exists wallCache($f)]} {set wallCache($f) [isWallRaw $f $profile]}
    return $wallCache($f)
}
proc ::SHR11::isWallRaw {f profile} {
    variable fn; variable xyz; variable fc; variable normal; variable tol
    lassign $profile c a r
    if {abs([dot $normal($f) $a])>0.35} {return 0}
    foreach n $fn($f) {
        if {abs([lindex [profileDistance $xyz($n) $profile] 0])>$tol*$r} {return 0}
    }
    lassign [profileDistance $fc($f) $profile] dist inward
    expr {[dot $normal($f) $inward]>0.5}
}
proc ::SHR11::wall {circle} {
    variable ef; variable fn
    lassign $circle c a r ns
    set queue {}; array set seen {}
    foreach edge [edges $ns] {
        if {![info exists ef($edge)]} {return {}}
        set found 0
        foreach f $ef($edge) {
            if {[isWall $f $circle]} {
                set found 1
                if {![info exists seen($f)]} {set seen($f) 1; lappend queue $f}
            }
        }
        if {!$found} {return {}}
    }
    for {set i 0} {$i<[llength $queue]} {incr i} {
        set f [lindex $queue $i]
        foreach edge [edges $fn($f)] {
            foreach g $ef($edge) {
                if {[info exists seen($g)]} {continue}
                set seen($g) 1
                if {[isWall $g $circle]} {lappend queue $g}
            }
        }
    }
    return $queue
}
# Is a face part of the outward annulus, in the mouth plane?
# Crucially: a disk/floor face reaches inward, or has no node outside the radius.
proc ::SHR11::isAnnulus {f profile} {
    variable annulusCache
    activateProfile $profile
    if {![info exists annulusCache($f)]} {set annulusCache($f) [isAnnulusRaw $f $profile]}
    return $annulusCache($f)
}
proc ::SHR11::isAnnulusRaw {f circle} {
    variable fn; variable xyz; variable normal; variable tol; variable fc
    lassign $circle c a r
    if {abs([dot $normal($f) $a])<0.94} {return 0}
    set outside 0; set eps [expr {1.e-7*$r}]
    foreach n $fn($f) {
        if {abs([dot [sub $xyz($n) $c] $a])>$tol*$r} {return 0}
        set d [lindex [profileDistance $xyz($n) $circle] 0]
        if {$d<-$eps} {return 0}
        if {$d>$eps} {set outside 1}
    }
    # Also rejects a floor triangulated using perimeter nodes only.
    if {[lindex [profileDistance $fc($f) $circle] 0]<-$eps} {return 0}
    return $outside
}
proc ::SHR11::mouth {circle} {
    variable ef; variable fn; variable nf; variable washer
    set rim [lindex $circle 3]
    # Every rim edge needs a verified outward annulus; bottom edges fail here.
    foreach edge [edges $rim] {
        if {![info exists ef($edge)]} {return {}}
        set ok 0
        foreach f $ef($edge) {if {[isAnnulus $f $circle]} {set ok 1; break}}
        if {!$ok} {return {}}
    }
    set selected {}; foreach n $rim {dict set selected $n 1}
    if {$washer} {
        # Vertex-adjacent surface layer handles triangulated washers too.
        set faceids {}
        foreach n $rim {foreach f $nf($n) {dict set faceids $f 1}}
        foreach f [dict keys $faceids] {
            if {![isAnnulus $f $circle]} {continue}
            foreach n $fn($f) {dict set selected $n 1}
        }
    }
    return [lsort -integer [dict keys $selected]]
}
proc ::SHR11::solidLayer {wallnodes circle} {
    variable xyz; variable ne; variable en; variable washer
    set selected {}; foreach n $wallnodes {dict set selected $n 1}
    if {!$washer} {return [dict keys $selected]}
    lassign $circle c a r
    set zmin Inf; set zmax -Inf; set solids {}
    foreach n $wallnodes {
        set z [dot [sub $xyz($n) $c] $a]
        set zmin [expr {min($zmin,$z)}]; set zmax [expr {max($zmax,$z)}]
        foreach e $ne($n) {dict set solids $e 1}
    }
    # Numerical tolerance only: do not admit a physical layer below a blind floor.
    set eps [expr {max(1.e-10,1.e-7*$r)}]
    set checked $selected
    foreach e [dict keys $solids] {
        foreach n $en($e) {
            if {[dict exists $checked $n]} {continue}
            dict set checked $n 1
            set v [sub $xyz($n) $c]; set z [dot $v $a]
            if {$z<$zmin-$eps || $z>$zmax+$eps} {continue}
            # Wall nodes are always retained, including floor perimeter nodes.
            # Other nodes inside the cylinder (the blind floor) are excluded.
            if {[lindex [profileDistance $xyz($n) $circle] 0]<-$eps} {continue}
            dict set selected $n 1
        }
    }
    return [lsort -integer [dict keys $selected]]
}
proc ::SHR11::recognize {} {
    variable activeProfile; unset -nocomplain activeProfile
    variable fn; variable nf; variable through; variable tol
    array unset nf
    foreach f [array names fn] {foreach n $fn($f) {lappend nf($n) $f}}
    set circles {}; set mouths {}
    foreach loop [loops] {
        set c [circle $loop]
        if {![llength $c]} {continue}
        set nodes [mouth $c]
        if {![llength $nodes]} {continue}
        lappend circles $c; lappend mouths $nodes
    }
    array set done {}; set result {}
    foreach c $circles {
        set ns [lindex $c 3]
        if {[info exists done([lindex $ns 0])]} {continue}
        activateProfile $c
        set wall [wall $c]; if {![llength $wall]} {continue}
        array unset wn; array set wn {}
        foreach f $wall {foreach n $fn($f) {set wn($n) 1}}
        set rims {}; set nodes {}
        # One solid layer at EVERY depth; suppress blind-floor interior/below-floor nodes.
        foreach n [solidLayer [array names wn] $c] {dict set nodes $n 1}
        foreach c2 $circles mn $mouths {
            set all 1
            foreach n [lindex $c2 3] {if {![info exists wn($n)]} {set all 0; break}}
            if {$all} {
                lappend rims $c2
                foreach n $mn {dict set nodes $n 1}
            }
        }
        foreach rim $rims {foreach n [lindex $rim 3] {set done($n) 1}}
        if {$through && [llength $rims]!=2} {continue}
        if {[llength $rims]>2 || ![llength $rims]} {continue}
        # Use the FINAL unique dependent-node set for BOTH blind and through holes.
        # Arithmetic coordinate mean: no mouth-plane placement or axis projection.
        set dependents [lsort -integer [dict keys $nodes]]
        set center [mean $dependents]
        if {[llength $rims]==2} {
            set gap [norm [sub [lindex [lindex $rims 0] 0] [lindex [lindex $rims 1] 0]]]
            if {$gap<$tol*[lindex $c 2]} {continue}
        }
        # Persist geometry before the next component replaces the local mesh cache.
        set mouthcenters {}; foreach rim $rims {lappend mouthcenters [lindex $rim 0]}
        lassign $c ref axis radius
        set lo Inf; set hi -Inf
        foreach n [array names wn] {
            variable xyz
            set z [dot [sub $xyz($n) $ref] $axis]
            set lo [expr {min($lo,$z)}]; set hi [expr {max($hi,$z)}]
        }
        set meta [dict create axis $axis ref $ref mouths $mouthcenters ends [list [add $ref [mul $axis $lo]] [add $ref [mul $axis $hi]]]]
        lappend result [list $center $dependents [expr {2*[lindex $c 2]}] [llength $rims] $meta]
    }
    return $result
}
proc ::SHR11::reserveExisting {} {
    variable reserved; array unset reserved
    *createmark elems 1 all
    set ids [hm_getmark elems 1]
    dict for {e row} [readFields elems $ids {config}] {
        set cfg [lindex $row 0]
        if {$cfg ni {5 55 56}} {continue}
        # Conservatively reserve BOTH independent and dependent nodes.
        foreach n [hm_getvalue elems id=$e dataname=nodes] {set reserved($n) $e}
    }
}
proc ::SHR11::createSpider {center ns cname} {
    variable dofs
    *currentcollector component $cname
    *elementtype 55 1
    set node 0; set elem 0
    set rc [catch {
        lassign $center x y z
        *createnode $x $y $z 0 0 0
        set node [hm_latestentityid nodes]
        *createmark nodes 2 {*}$ns
        *rigidlink $node 2 $dofs
        set elem [hm_latestentityid elems]
        set tn [string toupper [hm_getvalue elems id=$elem dataname=typename]]
        if {$tn ne "RBE2"} {error "Created type is '$tn', not RBE2. Load OptiStruct/Nastran profile."}
        set actual [hm_getvalue elems id=$elem dataname=independentnode]
        if {$actual!=$node} {error "Independent node verification failed"}
        # Clear ONLY this newly created node from the temporary display mark.
        # Never delete the independent GRID belonging to an existing RBE2.
        *createmark nodes 2 $node
        *nodemarkcleartempmark 2
    } msg opts]
    if {$rc} {
        if {$elem>0} {catch {*createmark elems 2 $elem; *deletemark elems 2}}
        if {$node>0} {catch {*createmark nodes 2 $node; *deletemark nodes 2}}
        return -options $opts $msg
    }
    return [list $elem $node]
}
proc ::SHR11::executeSingle {} {
    variable dmin; variable dmax; variable tol; variable angle; variable dofs
    variable reserved; variable made; variable skipped; variable status
    set made 0; set skipped 0
    foreach x [list $dmin $dmax $tol $angle] {if {![string is double -strict $x]} {error "Parameters must be numbers"}}
    if {$dmin<=0 || $dmax<$dmin || $tol<=0 || $tol>0.1 || $angle<46 || $angle>85} {
        error "Use 0 < Dmin <= Dmax; 0 < tolerance <= 0.1; 46 <= angle <= 85."
    }
    if {![regexp {^[1-6]+$} $dofs] || [llength [lsort -unique [split $dofs {}]]]!=[string length $dofs]} {error "Invalid DOF digits"}
    *createmarkpanel comps 1 "Select complete SOLID components; then Proceed"
    set cids [hm_getmark comps 1]
    if {![llength $cids]} {set status "Cancelled: no components selected."; return}
    reserveExisting
    set names {}
    *createmark comps 1 all
    foreach id [hm_getmark comps 1] {dict set names [hm_getvalue comps id=$id dataname=name] 1}
    foreach cid $cids {
        set src [hm_getvalue comps id=$cid dataname=name]
        set status "Reading component $cid: $src"; update idletasks
        puts "SHR: processing $cid ($src)"
        set loadStart [clock milliseconds]
        if {[catch {set cnt [loadComponent $cid]} msg]} {puts "SHR SKIP: $msg"; incr skipped; continue}
        if {!$cnt} {continue}
        set candidates [recognize]; set output ""
        puts "SHR: [llength $candidates] hole candidates"
        foreach hole $candidates {
            lassign $hole center ns diameter mouths
            set conflict 0
            foreach n $ns {if {[info exists reserved($n)]} {set conflict 1; break}}
            if {$conflict} {puts "SHR SKIP: D=$diameter overlaps rigid element $reserved($n) at node $n"; incr skipped; continue}
            if {$output eq ""} {
                set base "RBE2_HOLES_${cid}_${src}"; set output $base; set i 1
                while {[dict exists $names $output]} {set output "${base}_[incr i]"}
                *createentity comps name=$output
                dict set names $output 1
            }
            lassign [createSpider $center $ns $output] eid master
            foreach n $ns {set reserved($n) $eid}; set reserved($master) $eid
            incr made
            puts "SHR CREATED: comp=$output RBE2=$eid independent=$master D=$diameter dependent_nodes=[llength $ns] mouths=$mouths"
        }
    }
    set status "Finished: $made RBE2 created; $skipped conflicts/components skipped. See Tcl console."
    puts "SHR: $status"
}
# Returns true only for adjacent, coaxial holes with facing OPEN mouths.
proc ::SHR11::compatible {ha hb} {
    variable alignOnly
    variable axisTol; variable angleTol; variable diamTol; variable gapTol
    if {[dict get $ha cid]==[dict get $hb cid]} {return 0}
    set a [dict get $ha hole]; set b [dict get $hb hole]
    if {abs([lindex $a 2]-[lindex $b 2])>$diamTol+1.e-9} {return 0}
    set ma [lindex $a 4]; set mb [lindex $b 4]
    set ua [dict get $ma axis]; set ub [dict get $mb axis]
    if {abs([dot $ua $ub])<cos($angleTol*acos(-1)/180.0)} {return 0}
    set ra [dict get $ma ref]; set rb [dict get $mb ref]
    # Both endpoints must stay near the opposite axis (also limits angular drift).
    foreach point [dict get $ma ends] {
        set d [sub $point $rb]
        if {[norm [sub $d [mul $ub [dot $d $ub]]]]>$axisTol+1.e-9} {return 0}
    }
    foreach point [dict get $mb ends] {
        set d [sub $point $ra]
        if {[norm [sub $d [mul $ua [dot $d $ua]]]]>$axisTol+1.e-9} {return 0}
    }
    set aa {}; foreach p [dict get $ma ends] {lappend aa [dot [sub $p $ra] $ua]}
    set bb {}; foreach p [dict get $mb ends] {lappend bb [dot [sub $p $ra] $ua]}
    lassign [lsort -real $aa] alo ahi; lassign [lsort -real $bb] blo bhi
    set eps [expr {1.e-7*max(1.0,[lindex $a 2],[lindex $b 2])}]
    # V7 alignment mode uses axis/diameter and axial separation only.
    # Does not demand a facing mouth or reject overlapping axial intervals.
    if {$alignOnly} {
        set separation [expr {max(0.0,$blo-$ahi,$alo-$bhi)}]
        return [expr {$separation<=$gapTol+$eps}]
    }
    if {$ahi<=$blo+$eps} {set ta $ahi; set tb $blo; set gap [expr {$blo-$ahi}]} elseif {$bhi<=$alo+$eps} {
        set ta $alo; set tb $bhi; set gap [expr {$alo-$bhi}]
    } else {return 0}
    if {$gap>$gapTol+$eps} {return 0}
    set openA 0; set openB 0
    foreach p [dict get $ma mouths] {if {abs([dot [sub $p $ra] $ua]-$ta)<=$eps} {set openA 1}}
    foreach p [dict get $mb mouths] {if {abs([dot [sub $p $ra] $ua]-$tb)<=$eps} {set openB 1}}
    expr {$openA && $openB}
}
# Diagnostic only: measure even candidates rejected by the spatial broadphase.
proc ::SHR11::diagnosePair {ha hb} {
    variable axisTol; variable angleTol; variable diamTol; variable gapTol
    set a [dict get $ha hole]; set b [dict get $hb hole]
    set ma [lindex $a 4]; set mb [lindex $b 4]
    set ua [dict get $ma axis]; set ub [dict get $mb axis]
    set ra [dict get $ma ref]; set rb [dict get $mb ref]
    set theta [expr {acos(min(1.0,abs([dot $ua $ub])))*180.0/acos(-1)}]
    set offset 0.0
    foreach {m ref axis} [list $ma $rb $ub $mb $ra $ua] {
        foreach p [dict get $m ends] {
            set d [sub $p $ref]
            set offset [expr {max($offset,[norm [sub $d [mul $axis [dot $d $axis]]]])}]
        }
    }
    set aa {}; set bb {}
    foreach p [dict get $ma ends] {lappend aa [dot [sub $p $ra] $ua]}
    foreach p [dict get $mb ends] {lappend bb [dot [sub $p $ra] $ua]}
    lassign [lsort -real $aa] alo ahi; lassign [lsort -real $bb] blo bhi
    set gap [expr {max(0.0,$blo-$ahi,$alo-$bhi)}]
    set dd [expr {abs([lindex $a 2]-[lindex $b 2])}]
    set eps [expr {1.e-7*max(1.0,[lindex $a 2],[lindex $b 2])}]
    set reasons {}
    if {$offset>$axisTol+1.e-9} {lappend reasons AXIS_OFFSET}
    if {$theta>$angleTol+1.e-9} {lappend reasons AXIS_ANGLE}
    if {$dd>$diamTol+1.e-9} {lappend reasons DIAMETER}
    if {$gap>$gapTol+$eps} {lappend reasons AXIAL_GAP}
    if {![llength $reasons]} {
        if {[compatible $ha $hb]} {set reasons PASS} else {set reasons FACING_MOUTH_OR_OVERLAP}
    }
    return [format "axisOffset=%.6g(limit=%.6g) angle=%.6g(limit=%.6g) diameterDiff=%.6g(limit=%.6g) axialGap=%.6g(limit=%.6g) result=%s refA=%s refB=%s" $offset $axisTol $theta $angleTol $dd $diamTol $gap $gapTol $reasons $ra $rb]
}
proc ::SHR11::diagnoseUnmatched {holes indices} {
    # Cap diagnostics: nearest three opposing hole-axis midpoints for 12 holes.
    # Nearest spatial neighbour is NOT asserted to be the intended bolt partner.
    set centers {}
    foreach h $holes {
        lassign [dict get [lindex [dict get $h hole] 4] ends] a b
        lappend centers [mul [add $a $b] 0.5]
    }
    foreach i [lrange $indices 0 11] {
        set h [lindex $holes $i]; set choices {}; set j -1
        foreach other $holes {
            incr j
            if {[dict get $h cid]==[dict get $other cid]} {continue}
            lappend choices [list [norm [sub [lindex $centers $i] [lindex $centers $j]]] $j]
        }
        foreach choice [lrange [lsort -real -index 0 $choices] 0 2] {
            set j [lindex $choice 1]; set other [lindex $holes $j]
            puts "SHR DIAG hole=$i comp=[dict get $h cid] neighbour=$j comp=[dict get $other cid]: [diagnosePair $h $other]"
        }
    }
    if {[llength $indices]>12} {puts "SHR DIAG: output limited to first 12 unmatched/ambiguous holes."}
}
proc ::SHR11::pairHoles {holes} {
    set count [llength $holes]; array set adj {}
    for {set i 0} {$i<$count} {incr i} {set adj($i) {}}
    variable gapTol; variable axisTol
    set boxes {}; set maxD 1.0
    for {set i 0} {$i<$count} {incr i} {
        set h [lindex $holes $i]; set hole [dict get $h hole]
        set maxD [expr {max($maxD,[lindex $hole 2])}]
        lassign [dict get [lindex $hole 4] ends] a b
        lassign $a ax ay az; lassign $b bx by bz
        lappend boxes [list [expr {min($ax,$bx)}] [expr {max($ax,$bx)}] [expr {min($ay,$by)}] [expr {max($ay,$by)}] [expr {min($az,$bz)}] [expr {max($az,$bz)}] $i]
    }
    set boxes [lsort -real -index 0 $boxes]
    set pad [expr {$gapTol+$axisTol+1.e-7*$maxD+1.e-8}]
    set comparisons 0
    for {set k 0} {$k<$count} {incr k} {
        lassign [lindex $boxes $k] xmin xmax ymin ymax zmin zmax i
        for {set l [expr {$k+1}]} {$l<$count} {incr l} {
            lassign [lindex $boxes $l] x0 x1 y0 y1 z0 z1 j
            if {$x0>$xmax+$pad} {break}
            if {$y0>$ymax+$pad || $ymin>$y1+$pad || $z0>$zmax+$pad || $zmin>$z1+$pad} {continue}
            # Keep original argument order and final ordering exactly as before.
            set first [expr {min($i,$j)}]; set second [expr {max($i,$j)}]
            if {[dict get [lindex $holes $first] cid]==[dict get [lindex $holes $second] cid]} {continue}
            incr comparisons
            if {[compatible [lindex $holes $first] [lindex $holes $second]]} {lappend adj($i) $j; lappend adj($j) $i}
        }
    }
    for {set i 0} {$i<$count} {incr i} {set adj($i) [lsort -integer $adj($i)]}
    puts "SHR PERF: $count holes; $comparisons cross-component detailed pair checks (all unordered pairs: [expr {$count*($count-1)/2}])"
    set pairs {}; set unused 0; set unmatched {}
    for {set i 0} {$i<$count} {incr i} {
        if {[llength $adj($i)]==1} {
            set j [lindex $adj($i) 0]
            if {[llength $adj($j)]==1} {if {$i<$j} {lappend pairs [list $i $j]}; continue}
        }
        incr unused
        lappend unmatched $i
        set h [lindex $holes $i]
        if {$unused<=12} {puts "SHR SKIP hole $i in component [dict get $h cid]: candidates=$adj($i)"}
    }
    if {$unused} {diagnoseUnmatched $holes $unmatched}
    return [list $pairs $unused]
}
proc ::SHR11::executePairs {} {
    variable dmin; variable dmax; variable tol; variable angle; variable dofs
    variable axisTol; variable angleTol; variable diamTol; variable gapTol
    variable reserved; variable made; variable skipped; variable status; variable xyz
    set made 0; set skipped 0
    foreach x [list $dmin $dmax $tol $angle $axisTol $angleTol $diamTol $gapTol] {
        if {![string is double -strict $x]} {error "Parameters must be numeric"}
    }
    if {$dmin<=0 || $dmax<$dmin || $tol<=0 || $tol>0.1 || $angle<46 || $angle>85 || $axisTol<0 || $angleTol<0 || $angleTol>15 || $diamTol<0 || $gapTol<0} {error "Invalid diameter, tolerance, angle or gap range"}
    if {![regexp {^[1-6]+$} $dofs] || [llength [lsort -unique [split $dofs {}]]]!=[string length $dofs]} {error "Invalid DOF digits"}
    *createmarkpanel comps 1 "Select TWO OR MORE complete solid components; then Proceed"
    set cids [hm_getmark comps 1]
    if {[llength $cids]<2} {set status "Select at least two solid components for pair mode."; return}
    reserveExisting
    set names {}; set holes {}; set coordinates {}
    *createmark comps 1 all
    foreach id [hm_getmark comps 1] {dict set names [hm_getvalue comps id=$id dataname=name] 1}
    foreach cid $cids {
        set src [hm_getvalue comps id=$cid dataname=name]
        set status "Recognizing holes in $cid ($src)"; update idletasks
        set loadStart [clock milliseconds]
        if {[catch {set cnt [loadComponent $cid]} msg]} {puts "SHR SKIP: $msg"; incr skipped; continue}
        if {!$cnt} {continue}
        set recognizeStart [clock milliseconds]
        set found [recognize]
        puts "SHR PERF component=$cid load=[expr {($recognizeStart-$loadStart)/1000.0}]s recognize=[expr {([clock milliseconds]-$recognizeStart)/1000.0}]s"
        puts "SHR: component $cid, [llength $found] holes"
        # Candidate geometry and node IDs only: NO master, RBE2 or component writes.
        foreach hole $found {
            set ns [lindex $hole 1]; set conflict 0
            foreach n $ns {if {[info exists reserved($n)]} {set conflict 1; break}}
            if {$conflict} {puts "SHR SKIP candidate: existing rigid $reserved($n) at node $n"; incr skipped; continue}
            lappend holes [dict create cid $cid name $src hole $hole]
            foreach n $ns {dict set coordinates $n $xyz($n)}
        }
    }
    lassign [pairHoles $holes] pairs unused
    variable diagnosticOnly
    if {$diagnosticOnly} {
        set status "DIAGNOSTIC ONLY: [llength $pairs] pairs; $unused unmatched/ambiguous; no entities created."
        puts "SHR: $status"; return
    }
    set outputs {}
    foreach pair $pairs {
        lassign $pair ia ib
        set a [lindex $holes $ia]; set b [lindex $holes $ib]
        set key [lsort -integer [list [dict get $a cid] [dict get $b cid]]]
        set ns [lsort -integer -unique [concat [lindex [dict get $a hole] 1] [lindex [dict get $b hole] 1]]]
        set conflict 0
        foreach n $ns {if {[info exists reserved($n)]} {set conflict 1; break}}
        if {$conflict} {puts "SHR SKIP pair $key: node $n already in rigid $reserved($n)"; incr skipped; continue}
        # All unique dependents from BOTH components, with equal node weight.
        set center {0 0 0}
        foreach n $ns {set center [add $center [dict get $coordinates $n]]}
        set center [mul $center [expr {1.0/[llength $ns]}]]
        if {![dict exists $outputs $key]} {
            lassign $key ca cb
            set base "RBE2_PAIR_${ca}_${cb}"; set out $base; set suffix 1
            while {[dict exists $names $out]} {set out "${base}_[incr suffix]"}
            *createentity comps name=$out
            dict set names $out 1; dict set outputs $key $out
        }
        set out [dict get $outputs $key]
        lassign [createSpider $center $ns $out] eid master
        foreach n $ns {set reserved($n) $eid}; set reserved($master) $eid
        incr made
        puts "SHR CREATED: pair=$key comp=$out RBE2=$eid independent=$master dependents=[llength $ns]"
    }
    set status "Finished: $made merged RBE2; $unused unmatched/ambiguous candidates ignored (no elements created); $skipped conflicts/components skipped."
    puts "SHR: $status"
}
proc ::SHR11::execute {} {
    variable pairmode; variable sourcePath
    puts "SHR V11 DIRECT: source=$sourcePath"
    set start [clock milliseconds]
    set rc [catch {if {$pairmode} {executePairs} else {executeSingle}} msg opts]
    puts "SHR PERF: total elapsed [format %.3f [expr {([clock milliseconds]-$start)/1000.0}]] seconds"
    if {$rc} {return -options $opts $msg}
}
proc ::SHR11::run {} {
    variable busy; variable status; variable made
    if {$busy} {return}; set busy 1
    .shr.go configure -state disabled
    set rc [catch {execute} msg opts]
    set busy 0; .shr.go configure -state normal
    if {$rc} {
        set status "Stopped after $made RBE2: $msg"
        puts [dict get $opts -errorinfo]
        tk_messageBox -icon error -title "Solid hole RBE2" -message $status
    } else {tk_messageBox -icon info -title "Solid hole RBE2" -message $status}
}
proc ::SHR11::ui {} {
    if {[winfo exists .shr]} {destroy .shr}
    toplevel .shr; wm title .shr "HM2021 Hole RBE2 V11 - direct connection holes"
    set row 0
    foreach {v label} {dmin "Minimum equivalent-area diameter" dmax "Maximum equivalent-area diameter" tol "Profile / plane tolerance (relative, 0.04)" angle "Feature angle, degrees (55)" dofs "RBE2 dependent DOFs" axisTol "Pair: maximum axis offset (model units)" angleTol "Pair: maximum axis angle (degrees)" diamTol "Pair: maximum diameter difference" gapTol "Pair: maximum axial mouth gap"} {
        label .shr.l$v -text $label -anchor w
        entry .shr.e$v -textvariable ::SHR11::$v -width 16
        grid .shr.l$v -row $row -column 0 -sticky w -padx 10 -pady 4
        grid .shr.e$v -row $row -column 1 -padx 10 -pady 4; incr row
    }
    checkbutton .shr.pair -text "DIRECT PAIRS: match holes first, create final RBE2 only" -variable ::SHR11::pairmode
    grid .shr.pair -row $row -columnspan 2 -sticky w -padx 10; incr row
    checkbutton .shr.diag -text "DIAGNOSTIC ONLY (pair mode): no entities created" -variable ::SHR11::diagnosticOnly
    grid .shr.diag -row $row -columnspan 2 -sticky w -padx 10; incr row
    checkbutton .shr.align -text "Match by hole alignment (ignore facing-mouth test)" -variable ::SHR11::alignOnly
    grid .shr.align -row $row -columnspan 2 -sticky w -padx 10; incr row
    checkbutton .shr.through -text "Require through holes (two mouths -> one RBE2)" -variable ::SHR11::through
    grid .shr.through -row $row -columnspan 2 -sticky w -padx 10; incr row
    checkbutton .shr.washer -text "Include ONE solid layer along full depth (inner + outer nodes)" -variable ::SHR11::washer
    grid .shr.washer -row $row -columnspan 2 -sticky w -padx 10; incr row
    label .shr.note -justify left -text "OptiStruct/Nastran profile required. Save model first.\nSelect whole components, including their hidden mesh.\nLinear solids; closed near-planar profiles, at least 4 rim nodes.\nOne solid layer at all depths; blind-floor interior/below-floor nodes excluded.\nMaster temporary marks cleared; required GRID is retained.\nDelete old RBE2 elements on target holes before rebuilding.\nPair mode: ambiguous/unmatched holes skipped; no single-hole fallback.\nTwo-part pairs only. No temporary single-hole RBE2 are created.\nExisting rigid-node overlap skips the entire hole."
    grid .shr.note -row $row -columnspan 2 -padx 10 -pady 8; incr row
    button .shr.go -text "Select components and CREATE RBE2" -command ::SHR11::run
    grid .shr.go -row $row -columnspan 2 -pady 8; incr row
    label .shr.status -textvariable ::SHR11::status -wraplength 480 -justify left
    grid .shr.status -row $row -columnspan 2 -padx 10 -pady 8
    wm protocol .shr WM_DELETE_WINDOW {if {!$::SHR11::busy} {destroy .shr}}
}
puts "SHR V11 DIRECT loaded: $::SHR11::sourcePath"
if {![info exists ::SHR11_NO_UI]} {::SHR11::ui}
