##################################################################
#  TCL Script to visualize CG structures in VMD
#  
#  Version 3.0
##################################################################
#  
# Adapted from vmds
#
##################################################################

# ==================================================================================================
# SETUP THE CG REPRESENTATION
# ==================================================================================================
proc g_cg { args } {

	### Some variable need to be global to the file
	global helix_radius
	global helix_color
	global helix_angle
	global helix_angle_tol
	global CG_bb
	global CG_bb_index
	global CG_bb_num
	global CG_prot_res
	global bond
	
	### Set defaults
	set NAME 		"\[ g_cg \]"
	set molid		top;				# Molecule Id
	set first		0;				# The first frame to process
	set last		-1;				# The last frame to process
	set skip		1;				# The step between 2 frames
	set nframe		0;				# number of frame
	set sel_list		{{name W WF} {name "BH.*" "SH.*"} {not name W WF "BH.*" "SH.*"}};	# Selection list
	set rep_list		{Lines Licorice Lines};							# Representation list
	set color_list		{"ColorId 0" "ColorId 7" "ColorId 10"};					# Color list
	set calpha		"name \"BH.*\"";		# Calpha Selection text
	set other		"not name \"BH.*\" \"SH.*\"";	# Non-Protein Selection text
	set tpr			"topol.tpr";			# GROMACS TPR file
	set helix_radius	3;				# Cylinder representation for helix
	set helix_color		"red";
	set helix_angle		57.0;				# Helix definition parameters
	set helix_angle_tol	5.0;				# Helix definition parameters
	set gdistrib		"/usr/local/gromacs";	# GROMACS distribution
	
	### Print Help message if no args
	if {[llength $args] == 0} {
		puts "Synopsis:"
		puts "========="
		puts "	Rebuild the bonds between th particles of a coarse grained system"
		puts "	Require a a GROMACS topology file and gmxdump"
		puts ""
		puts "Usage:"
		puts "======"
		puts "	g_cg \[options] -tpr a/gromacs/tpr/file.tpr"
		puts ""
		puts "Options:"
		puts "========"
		puts "	-molid		top				Molecule Id"
		puts "	-tpr		topo.tpr			GROMACS tolopolgy file corresping to your system"
		puts "	-sel		see below			List of the selections to display"
		puts "	-rep		see below			List of the representations to apply on each selection"
		puts "	-color		see below			List of the color to apply on each selection"
		puts "	-calpha		name \"BH.*\"			Calpha particle (for the helical residues detection)"
		puts "	-hangle		57.0				Angle applied to retrieve helical residues"
		puts "	-htol		5.0				Angle tolerance applied to retrieve helical residue"
		puts "	-hradius	3				Radius of the cylinder for the helix representation"
		puts "	-hcolor		red				Color of the helix representation"
		puts "	-distrib	/usr/export/gromacs-3.3.1	Path to a GROMACS distribution"
		puts ""
		puts "Details of the default representations:"
		puts "======================================="
		puts "A) Default representation"
		puts "	Compounds	water		Protein				Other"
		puts "	Selection	name W WF	name \"BH.*\" \"SH.*\"		not name W WF \"BH.*\" \"SH.*\""
		puts "  Representation	Lines		Licorice			Lines"
		puts "	Color		ColorId 0	ColorId 7			ColorId 10"
		puts ""
		puts "B) Modified the default representation"
		puts "Warning! The lists of selections (-sel), representations (-rep) and colors (-color)"
		puts "-MUST- have the same size!"
		puts ""
		puts "	1. represent only chain A of a protein"
		puts "		g_cg -tpr topol.tpr -sel {chain A and name \"BH.*\" \"SH.*\"} -rep {Licorice} -color {NAME}"
		puts "	2. represent a protein and lipids in 2 different styles"
		puts "		g_cg -tpr topol.tpr -sel {{name \"BH.*\" \"SH.*\"} {resname DOPC}} -rep {CPK Line} -color {Residue Blue}"
		return
	}
	
	### Parse options with one argument
	foreach {i j} $args {
		if { $i=="-tpr"    } { set tpr			$j }
		if { $i=="-sel"    } { set sel_list		$j }
		if { $i=="-rep"    } { set rep_list		$j }
		if { $i=="-color"  } { set color_list		$j }
		if { $i=="-calpha" } { set calpha		$j }
		if { $i=="-other"  } { set other		$j }
		if { $i=="-molid"  } { set molid		$j }
		if { $i=="-hradius"} { set helix_radius		$j }
		if { $i=="-hangle" } { set helix_angle		$j }
		if { $i=="-htol"   } { set helix_angle_tol	$j }
		if { $i=="-hcolor" } { set helix_color		$j }
		if { $i=="-distrib"} { set gdistrib		$j }
		
	}

	
	### Do some checking
	set nsel [llength $sel_list]
	set nrep [llength $rep_list]
	set ncol [llength $color_list]

	if { $nsel != $nrep || $nsel != $ncol || $nrep != $ncol } {
		puts "$NAME -- Lists of selections, representations and colors -MUST- have the same size"
		return
	}

	if {[file exists $tpr]} {
		set f [open "|$gdistrib/bin/gmxdump -s $tpr 2> /dev/null" r]
	} else {
		puts "$NAME -- Cannot open \"$tpr\""
		return
	}

	## Process the .tpr file
	puts "$NAME Processing \"$tpr\"..."
	array unset bond
	while { [gets $f line]>=0 } {
		# Parse line
		regexp {natoms = (\d+)} $line dum N
	 
		if { [regexp {\(BONDS\)\s+(\d+)\s+(\d+)} $line dum n1 n2] } {
			if { [info exists bond($n1)] } {
				lappend bond($n1) $n2
				lappend bond($n2) $n1
			} else {
				set bond($n1) $n2
    				set bond($n2) $n1
			}
 		}
	 
		if { [regexp {\(CONSTR\)\s+(\d+)\s+(\d+)} $line dum n1 n2] } {
			if { [info exists bond($n1)] } {
				lappend bond($n1) $n2
				lappend bond($n2) $n1
			} else {
				set bond($n1) $n2
				set bond($n2) $n1
			}
		} 
	}

	### Create the bond list
	puts "$NAME Create the bond list for $N atoms..."
	set blist {}
	for {set i 0} {$i<$N} {incr i} {
		if { [info exists bond($i)] } {
			lappend blist $bond($i)
		} else {
			lappend blist {}
		}
	}

	set all [atomselect top all]

	### Rebuild the bonds
	puts "$NAME Rebuild bonds..."
	$all setbonds $blist
	

	### Set variables for the helical residue detection
	# Get all ca-ca dihedrals
	set CG_bb [atomselect top "$calpha"]
	set CG_bb_index [$CG_bb get index]
	set CG_bb_num [$CG_bb num]
	set CG_prot_res [$CG_bb get resid]
	
	### Create representation
	# Delete previous reprensentations
	puts "$NAME Create representations..."
	#puts "[molinfo $molid get numreps] representation"
	for { set i 0 } {$i < [molinfo $molid get numreps]} {incr i} { mol delrep 0 $molid }
	mol delrep 0 $molid
	#return

	# Create new ones
	for { set i 0 } { $i < $nsel } { incr i } {
		mol selection	 [lindex $sel_list $i]
		mol rep		 [lindex $rep_list $i]
		mol color	 [lindex $color_list $i]
		mol addrep	 $molid
	}
	
	### Create the cylinder representation for helices
	#puts "Rendering secondary structure for frame 0..."
	#draw_cg_helix 0
	#puts "Registering callback..."
	#trace variable vmd_frame w trace_func
	#puts "Done."
	close $f
	return
}

# ==================================================================================================
# GET_HELIX
# ==================================================================================================
proc get_cg_helix {frame} {
	global CG_bb_index
	global CG_bb_num
	global bond helix_angle helix_angle_tol
 
	set helix {}
	for {set i 0} {$i<$CG_bb_num} {incr i} {
		set atoms [lrange $CG_bb_index $i [expr $i+3]]
		if {[llength $atoms]==4} {
		# See if all 4 atoms are bound together
			set ok 0
			if { [lsearch $bond([lindex $atoms 0]) [lindex $atoms 1]]>=0 } { incr ok }
			if { [lsearch $bond([lindex $atoms 1]) [lindex $atoms 2]]>=0 } { incr ok }
			if { [lsearch $bond([lindex $atoms 2]) [lindex $atoms 3]]>=0 } { incr ok }
			if {$ok == 3} {
				set val [measure dihed "$atoms" frame $frame]
				if {[expr abs($val-$helix_angle)] < $helix_angle_tol} {set helix "$helix $atoms"}
			}
		}
	}

	set helix [lsort -unique -integer $helix]
 
	# Now split the list into the bound sublists
	set helices {}
	set temp {}
	for {set i 0} {$i<[llength $helix]} {incr i} {
		if { [lsearch $bond([lindex $helix $i]) [lindex $helix [expr $i+1]]]>=0 } {
			# bound pair
			lappend temp [lindex $helix $i] [lindex $helix [expr $i+1]]
		} else {
			# non-bound pair found
			lappend helices [lsort -unique $temp]
			set temp {}
		}
	}
 
	return $helices
}

proc draw_cg_helix {frame} {
	global CG_bb_index
	global CG_bb_num
	global helix_radius helix_color

	draw delete all

	set hel [get_cg_helix $frame]
 
	foreach helix $hel {
		if {[llength $helix]>4} {
			# If the helix is long, find two points to determine axis
			set first4 [lrange $helix 0 3]
			set last4  [lrange $helix end-3 end]
			set sel_first [atomselect top "index $first4"]
			set sel_last [atomselect top "index $last4"]
			set a [measure center $sel_first]
			set b [measure center $sel_last]
			set c1 [lindex [$sel_first get "x y z"] 0]
			set c2 [lindex [$sel_last  get "x y z"] 3]
			$sel_first delete
			$sel_last  delete
			# Project first and last points to this line to get ends of cylinder
			set point1 [project2line $a $b $c1]
			set point2 [project2line $a $b $c2]
			draw color $helix_color
			draw cylinder $point1 $point2 radius $helix_radius filled yes resolution 36
		}
	}
}

# ==================================================================================================
# PROJECT2LINE
# ==================================================================================================
proc project2line {a b c} {
	set a1 [lindex $a 0]
	set a2 [lindex $a 1]
	set a3 [lindex $a 2]

	set b1 [lindex $b 0]
	set b2 [lindex $b 1]
	set b3 [lindex $b 2]

	set c1 [lindex $c 0]
	set c2 [lindex $c 1]
	set c3 [lindex $c 2]
 
	set C [expr $c1*$a1-$c1*$b1 + $c2*$a2-$c2*$b2 + $c3*$a3-$c3*$b3]
	set A [expr (($b1-$a1)*($b1-$a1) + ($b2-$a2)*($b2-$a2) + ($b3-$a3)*($b3-$a3))/($b1-$a1)]
	set B [expr ( ($a2*$b1-$a1*$b2)*($b2-$a2) + ($a3*$b1-$a1*$b3)*($b3-$a3)  )/($b1-$a1)]

	set p1 [expr (-$C-$B)/$A ]
	set p2 [expr ($p1*($b2-$a2) + $a2*$b1 - $a1*$b2)/($b1-$a1) ]
	set p3 [expr ($p1*($b3-$a3) + $a3*$b1 - $a1*$b3)/($b1-$a1) ]
 
	return "$p1 $p2 $p3"
}

# ==================================================================================================
# TCL callbacks to capture the change of frame
# ==================================================================================================
proc trace_func {args} {
	global vmd_frame
	#puts "called $vmd_frame([molinfo top])"
	draw_cg_helix $vmd_frame([molinfo top])
}
