#!/bin/sh
# the next line restarts using wish \
exec wish "$0" "$@"
#########################################################
wm title . "DISCRETA Solvers"

########################################################
#  Initial text
########################################################
set intro "\n\n \
The algorithms in DISCRETA for\n \
solving the diophantine linear systems \n \
were programmed by:\n\n \

Brendan McKay      <bdm@cs.anu.edu.au>: mckay
Kari Nurmela       <Kari.Nurmela@hut.fi> : tabudes makeKM
Patric Oestergard  <Patric.Ostergard@hut.fi> : wclique 
Alfred Wassermann  <Alfred.Wassermann@uni-bayreuth.de> : \n\
                   LLL spread ls_det ls_rand\
\n\n"
########################################################
#  Default values
########################################################
set lambda 1

set lllc0 20
set lllbeta 120
set lllp 18
set lllguess 1000

set clique_n 16
set clique_d 4
set clique_w 5

set tabudes_seed 2
set tabudes_verbose 2
set tabudes_maxiter 10000
set tabudes_attr "add 5 tl"
set tabudes_neigh "addone removeone correctrandom"
set tabudes_volcounts ""
set tabudes_neigh_cover "vchangeone empty correctrandom"
set tabudes_neigh_pack "empty vchangeone correctrandom"

set ls_n 2
set ls_maxrounds 50
set ls_maxsol 10
set ls_maxloops 2000000

set makekm_t 5
set makekm_k 6
set makekm_order 14
set makekm_filename ""

#############################
# Set initial filename
#############################
catch {
	set filehandle [open "km_fname" r]
	gets $filehandle filename
	close filehandle
}
#############################
# Get initial lambda value
# from command line
#############################
if {$argc > 0} {
	set lambda [lindex $argv 0]
}
########################################################
#  Colors and fonts and ...
########################################################
set filetypes {
    {{KM Files} 	      {KM*.txt}     }
    {{Gen. Files} 	   {gen*.txt}    }
    {{Text Files}       {.txt}        }
    {{All Files}        *             }
}
set filetypes2 {
    {{Gen. Files} 	   {gen*.txt}    }
    {{KM Files} 	      {KM*.txt}     }
    {{Text Files}       {.txt}        }
    {{All Files}        *             }
}

set logfile "solver.log"
set method ""
set processid -1
set proc_handle -1

set color1 wheat1
set color2 royalblue4
set color3 wheat2
option add *Background $color1
option add *Button.Foreground $color2
option add *Radiobutton.Foreground $color2
option add *Radiobutton.insertBackground $color3

#option add *font -adobe-helvetica-bold-r-normal--20-140-100-100-p-105-iso8859-1
#option add *font -adobe-helvetica-bold-r-normal--14-140-75-75-p-82-iso8859-1
option add *font -adobe-helvetica-bold-r-normal--17-120-100-100-p-92-iso8859-1
set font2 -*-courier-medium-r-*-*-14-*-*-*-*-*-*-*

########################################################
#  The frames partitioning the application window
########################################################
frame .headline         
frame .midpart           
frame .midpart.right    
frame .midpart.type     -relief groove -bd 2 
frame .bottom           
frame .bottom.l
frame .bottom.r

pack .headline -fill x
pack .midpart -anchor w -fill x
pack .midpart.type .midpart.right  -side left -anchor n
pack .bottom -fill x
pack .bottom.l .bottom.r -in .bottom -side left -fill x

########################################################
# The headline
########################################################

label .headline.text -text "DISCRETA solver" 
pack .headline.text -in .headline 

########################################################
# The left side: type of combinatorial object:
#   design, large set, covering, packing
########################################################
set p .midpart.type
button $p.designs   -text "t-designs"   -command {choose_problemtype designs}
button $p.coverings -text "coverings"   -command {choose_problemtype coverings}
button $p.packings  -text "packings"    -command {choose_problemtype packings}
button $p.largesets -text "large sets of t-designs"  -command {choose_problemtype largesets}
button $p.dismiss   -text "exit"        -width 40 \
   -command { if {$processid > 0} { stopit $processid $proc_handle } ; exit}

frame  $p.ffilename 
label  $p.lfilename -text "Filename:" 
entry  $p.filename  -textvariable filename -width 25 -relief sunken -bd 2 \
        -font $font2
button $p.selectfilename  -text "select" \
	-command { 
		set filename [tk_getOpenFile -filetypes $filetypes] 
		catch { set filename [ file tail $filename ] }
	} 
         
frame  $p.flambda   
label  $p.llambda   -text "Lambda:" 
entry  $p.lambda    -textvariable lambda -width 20 -relief sunken -bd 2 
#################################
pack $p.lfilename $p.filename -side left -in $p.ffilename -ipady 2m -ipadx 5 -fill x 
pack $p.selectfilename -in $p.ffilename -ipady 2m -side left 
pack $p.llambda $p.lambda   -side left -in $p.flambda     -ipady 2m -ipadx 10 -fill x 
#################################
pack $p.ffilename $p.flambda $p.designs $p.coverings $p.packings \
     $p.largesets $p.dismiss\
     -in $p -ipady 2m -ipadx 5m -fill x 
########################################################
#  Problem type
########################################################
set problemtype ""
set p .midpart.right.problemtype

frame $p 
pack $p -in .midpart.right -side top -anchor n -fill x

label $p.problemtype -textvariable problemtype 
pack $p.problemtype -in $p -ipady 2m -ipadx 10

########################################################
#
# On start some empty frames are set:
#
# Empty solver frame on startup
########################################################
set p .midpart.right.method
frame $p  -relief groove -bd 2 
pack $p   -in .midpart.right -side top -anchor n -fill x

########################################################
# Empty parameter frame on startup:
########################################################
set p .midpart.right
frame $p.param   
pack $p.param -in $p -fill x

########################################################
# open logfile window
########################################################
set pr .bottom.r
set pl .bottom.l
		
text $pr.text -relief groove -bd 2  -yscrollcommand "$pr.scroll set" \
     -width 72 -height 20  -font $font2
scrollbar $pr.scroll -command "$pr.text yview" 

button $pl.clear -text "Clear" -command "$pr.text delete 1.0 end" 
button $pl.stopit   -text "Stop it"     \
   -command { if {$processid > 0} { stopit $processid	$proc_handle} }
button $pl.showresult   -text "Show solution file"  \
   -command { if {$processid > 0} { showresult } }
button $pl.showkm   -text "Show Kramer Mesner file"  \
   -command { showkmfile $filename }
button $pl.startbutton -text "Start solver" -width 40 -bg $color3 \
	-command { start_solver $method}

pack $pr.scroll -in $pr -side right -fill y
pack $pr.text -in $pr -side right -fill x

pack $pl.clear -side bottom -ipady 2m -ipadx 2m -fill x
pack $pl.stopit -side bottom -ipady 2m -ipadx 2m -fill x
pack $pl.showresult -side bottom -ipady 2m -ipadx 2m -fill x
pack $pl.showkm -side bottom -ipady 2m -ipadx 2m -fill x
pack $pl.startbutton -in $pl -side bottom -ipady 2m -ipadx 5m -fill x

########################################################
#  The intro text
########################################################
set p .bottom.r
$p.text insert end "$intro\n"
$p.text yview moveto 1









########################################################
#
# Here are the procedures to react to
# the buttons
#
########################################################

########################################################
#  Show the available solvers for each problem
########################################################
proc show_solvers type {
	
	set p .midpart.right.method
	catch { destroy $p }
	frame $p     -relief groove -bd 2 
	
	switch $type {
		"designs" {
			set p1 $p.part1
			set p2 $p.part2
			frame $p1 
			frame $p2
			pack $p1 $p2 -in $p -side top -anchor w
			
		 	radiobutton $p1.lllwith -text "LLL with output" -variable method -value lllwith  \
				-command {make_lll_param .midpart.right.param}
			radiobutton $p1.lll     -text "LLL"				  -variable method -value lll  \
				-command {make_lll_param .midpart.right.param}
			radiobutton $p1.spread  -text "Spread"			  -variable method -value spread  \
				-command {make_no_param .midpart.right.param}
			radiobutton $p1.mckay   -text "McKay"			  -variable method -value mckay  \
				-command {make_no_param .midpart.right.param}

			radiobutton $p2.wclique -text "wclique"			  -variable method -value wclique  \
				-command {make_no_param .midpart.right.param}
			radiobutton $p2.wcliquedirect  -text "wclique directly"  -variable method -value wcliquedirect \
				-command {make_wcliquedirect_param .midpart.right.param}
			radiobutton $p2.tabudes -text "tabudes"			  -variable method -value tabudes  \
				-command {make_tabudes_param .midpart.right.param}
				
			radiobutton $p2.makekm -text "makeKM"			  -variable method -value makekm  \
				-command {make_makekm_param .midpart.right.param}
		
			pack $p1.lllwith $p1.lll $p1.spread $p1.mckay -in $p1\
			     -ipady 2m -ipadx 5m -side left -anchor w
			pack $p2.wclique $p2.wcliquedirect $p2.tabudes $p2.makekm -in $p2\
			     -ipady 2m -ipadx 5m -side left -anchor w
		}
		"coverings" { 
			radiobutton $p.tabudes -text "tabudes"			  -variable method -value tabudes \
				-command {make_tabudes_cover_param .midpart.right.param}
			radiobutton $p.makekm -text "makeKM"			  -variable method -value makekm  \
				-command {make_makekm_param .midpart.right.param}
				
			pack $p.tabudes $p.makekm \
			     -in $p -ipady 2m -ipadx 5m -side left -anchor w
		}
		"packings" { 
			set p1 $p.part1
			set p2 $p.part2
			frame $p1 
			frame $p2
			pack $p1 $p2 -in $p -side top -anchor w

			radiobutton $p1.wclique -text "wclique"			  -variable method -value wclique \
				-command {make_no_param .midpart.right.param}
			radiobutton $p1.wcliquedirect  -text "wclique directly"  -variable method -value wcliquedirect \
				-command {make_wcliquedirect_param .midpart.right.param}
			radiobutton $p1.spread  -text "Spread"			  -variable method -value spread \
				-command {make_no_param .midpart.right.param}
			radiobutton $p1.lll  -text "LLL" 			        -variable method -value lllpacking \
				-command {make_lllpacking_param .midpart.right.param}
			radiobutton $p2.tabudes -text "tabudes"			  -variable method -value tabudes \
				-command {make_tabudes_pack_param .midpart.right.param}

			radiobutton $p2.makekm -text "makeKM"			  -variable method -value makekm  \
				-command {make_makekm_param .midpart.right.param}

			pack $p1.wclique $p1.wcliquedirect $p1.spread $p1.lll \
			     -in $p1 -ipady 2m -ipadx 5m -side left -anchor w
			pack $p2.tabudes $p2.makekm \
			     -in $p2 -ipady 2m -ipadx 5m -side left -anchor w
		}
		"largesets" {
			radiobutton $p.lsdet -text "Large sets determ."	 -variable method -value lsdet \
				-command {make_lsdet_param .midpart.right.param}
			radiobutton $p.lsrand -text "Large sets random"	 -variable method -value lsrand \
				-command {make_lsrand_param .midpart.right.param}
			pack $p.lsdet $p.lsrand \
			     -in $p -ipady 2m -ipadx 5m -side left -anchor w
		}
	}
	pack $p -in .midpart.right -side top -anchor n -fill x
}

########################################################
#  The parameter section for LLL
########################################################
proc make_lll_param p {
	set pp .midpart.right
	catch { destroy $p }
	
	frame $p     -relief groove -bd 2 
	frame $p.l      
	frame $p.r   
	pack $p.l $p.r -in $p -side left -fill x
	
	label $p.l.lc0     -text "c0:"  
	entry $p.r.c0      -textvariable lllc0 -width 5 
	label $p.l.lbeta   -text "beta:" 
	entry $p.r.beta    -textvariable lllbeta -width 5 
	label $p.l.lp      -text "p:" 
	entry $p.r.p       -textvariable lllp -width 5 

	pack $p.l.lc0 $p.l.lbeta $p.l.lp -in $p.l -side top -anchor w -ipady 2m -ipadx 10
	pack $p.r.c0 $p.r.beta $p.r.p -in $p.r    -side top -anchor w -ipady 2m -ipadx 10
	pack $p -in .midpart.right -side top -anchor n -fill x
}

########################################################
#  The parameter section for LLL-packings
########################################################
proc make_lllpacking_param p {
	set pp .midpart.right
	catch { destroy $p }
	
	frame $p     -relief groove -bd 2 
	frame $p.l      
	frame $p.r   
	pack $p.l $p.r -in $p -side left -fill x
	
	label $p.l.lguess  -text "Diff. from Opt.:" 
	entry $p.r.guess   -textvariable lllguess -width 5 
	label $p.l.lbeta   -text "beta:" 
	entry $p.r.beta    -textvariable lllbeta -width 5 
	label $p.l.lp      -text "p:" 
	entry $p.r.p       -textvariable lllp -width 5 
	label $p.l.comment  -text "" 
	label $p.r.comment  -text "Diff...:  upper bound of difference to Steiner system." 

	pack $p.l.lguess $p.l.lbeta $p.l.lp $p.l.comment -in $p.l -side top -anchor w -ipady 2m -ipadx 10
	pack $p.r.guess $p.r.beta $p.r.p $p.r.comment -in $p.r    -side top -anchor w -ipady 2m -ipadx 10
	pack $p -in .midpart.right -side top -anchor n -fill x
}

########################################################
#  The parameter section for wcliquedirect
########################################################
proc make_wcliquedirect_param p {
	global font2
	
	set pp .midpart.right
	catch { destroy $p }
	
	frame $p     -relief groove -bd 2 
	frame $p.l      
	frame $p.r   
	pack $p.l $p.r -in $p -side left -fill x
	
	label $p.l.ln  -text "n:" 
	entry $p.r.n   -textvariable clique_n -width 5 
	label $p.l.ld  -text "d:" 
	entry $p.r.d   -textvariable clique_d -width 5 
	label $p.l.lw  -text "w:" 
	entry $p.r.w   -textvariable clique_w -width 5 
	label $p.l.comment  -text "" 
	label $p.r.comment  -text "Use file of generators as input !" 

	label $p.l.lgenfile  -text "Gener. filename:" 
	frame $p.r.f
	entry $p.r.f.genfile  -textvariable genfile -width 40 -font $font2
	button $p.r.f.selectfilename  -text "select" \
		-command { 
			set genfile [tk_getOpenFile -filetypes $filetypes2] 
			catch { set genfile [ file tail $genfile] }
		} 

	pack $p.r.f.genfile $p.r.f.selectfilename -in $p.r.f -side left \
	     -anchor w -ipady 2m -ipadx 5 -fill x

	pack $p.l.ln $p.l.ld $p.l.lw $p.l.lgenfile $p.l.comment -in $p.l -side top -anchor w -ipady 2m -ipadx 10
	pack $p.r.n $p.r.d $p.r.w $p.r.f $p.r.comment -in $p.r    -side top -anchor w -ipady 2m -ipadx 10
	pack $p -in .midpart.right -side top -anchor n -fill x
}

########################################################
#  The parameter section for tabudes for designs
########################################################
proc make_tabudes_param p {
	set pp .midpart.right
	catch { destroy $p }
	
	frame $p     -relief groove -bd 2 
	frame $p.l      
	frame $p.r   
	pack $p.l $p.r -in $p -side left -fill x
	
	label $p.l.lseed -text "Seed:" 
	entry $p.r.seed   -textvariable tabudes_seed -width 5 
	label $p.l.lverbose -text "Verbose:" 
	entry $p.r.verbose   -textvariable tabudes_verbose -width 5 
	label $p.l.lmaxiter -text "Max. Iter.:" 
	entry $p.r.maxiter   -textvariable tabudes_maxiter -width 10 
	label $p.l.lattr -text "Attr.:" 
	entry $p.r.attr   -textvariable tabudes_attr -width 40
	label $p.l.lneigh -text "Neighb.:" 
	entry $p.r.neigh   -textvariable tabudes_neigh -width 40

	pack $p.l.lmaxiter $p.l.lattr $p.l.lneigh $p.l.lseed $p.l.lverbose -in $p.l \
	-side top -anchor w -ipady 2m -ipadx 10
	pack $p.r.maxiter $p.r.attr $p.r.neigh $p.r.seed $p.r.verbose -in $p.r \
	-side top -anchor w -ipady 2m -ipadx 10
	pack $p -in .midpart.right -side top -anchor n -fill x
}

########################################################
#  The parameter section for tabudes for coverings
########################################################
proc make_tabudes_cover_param p {
	set pp .midpart.right
	catch { destroy $p }
	
	frame $p     -relief groove -bd 2 
	frame $p.l      
	frame $p.r   
	pack $p.l $p.r -in $p -side left -fill x
	
	label $p.l.lseed -text "Seed:" 
	entry $p.r.seed   -textvariable tabudes_seed -width 5 
	label $p.l.lverbose -text "Verbose:" 
	entry $p.r.verbose   -textvariable tabudes_verbose -width 5 
	label $p.l.lmaxiter -text "Max. Iter.:" 
	entry $p.r.maxiter   -textvariable tabudes_maxiter -width 10 
	label $p.l.lattr -text "Attr.:" 
	entry $p.r.attr   -textvariable tabudes_attr -width 40
	label $p.l.lneigh -text "Neighb.:" 
	entry $p.r.neigh   -textvariable tabudes_neigh_cover -width 40
	label $p.l.lvolcounts -text "Volcounts:" 
	entry $p.r.volcounts   -textvariable tabudes_volcounts -width 40

	pack $p.l.lvolcounts $p.l.lmaxiter $p.l.lattr $p.l.lneigh $p.l.lseed $p.l.lverbose -in $p.l \
	-side top -anchor w -ipady 2m -ipadx 10
	pack $p.r.volcounts $p.r.maxiter $p.r.attr $p.r.neigh $p.r.seed $p.r.verbose -in $p.r \
	-side top -anchor w -ipady 2m -ipadx 10
	pack $p -in .midpart.right -side top -anchor n -fill x
}

########################################################
#  The parameter section for tabudes for packings
########################################################
proc make_tabudes_pack_param p {
	set pp .midpart.right
	catch { destroy $p }
	
	frame $p     -relief groove -bd 2 
	frame $p.l      
	frame $p.r   
	pack $p.l $p.r -in $p -side left -fill x
	
	label $p.l.lseed -text "Seed:" 
	entry $p.r.seed   -textvariable tabudes_seed -width 5 
	label $p.l.lverbose -text "Verbose:" 
	entry $p.r.verbose   -textvariable tabudes_verbose -width 5 
	label $p.l.lmaxiter -text "Max. Iter.:" 
	entry $p.r.maxiter   -textvariable tabudes_maxiter -width 10 
	label $p.l.lattr -text "Attr.:" 
	entry $p.r.attr   -textvariable tabudes_attr -width 40
	label $p.l.lneigh -text "Neighb.:" 
	entry $p.r.neigh   -textvariable tabudes_neigh_pack -width 40
	label $p.l.lvolcounts -text "Volcounts:" 
	entry $p.r.volcounts   -textvariable tabudes_volcounts -width 40
	label $p.l.lhint -text "" 
	button $p.r.hint   -width 30 -text "Give me a hint for Volcounts" \
		-command { orbithint $filename } 

	pack $p.l.lvolcounts $p.l.lhint $p.l.lmaxiter $p.l.lattr $p.l.lneigh \
			$p.l.lseed $p.l.lverbose -in $p.l \
			-side top -anchor w -ipady 2m -ipadx 10
	pack $p.r.volcounts $p.r.hint  $p.r.maxiter $p.r.attr $p.r.neigh \
			$p.r.seed $p.r.verbose -in $p.r \
			-side top -anchor w -ipady 2m -ipadx 10
	pack $p -in .midpart.right -side top -anchor n -fill x
}

########################################################
#  The parameter section for tabudes for large sets determin.
########################################################
proc make_lsdet_param p {
	set pp .midpart.right
	catch { destroy $p }
	
	frame $p     -relief groove -bd 2 
	frame $p.l      
	frame $p.r   
	pack $p.l $p.r -in $p -side left -fill x
	
	label $p.l.ln  -text "N:" 
	entry $p.r.n   -textvariable ls_n -width 5 
	label $p.l.comment  -text "" 
	label $p.r.comment  -text "LLL parameters in t-designs are used also" 

	pack $p.l.ln $p.l.comment -in $p.l -side top -anchor w -ipady 2m -ipadx 10
	pack $p.r.n $p.r.comment -in $p.r    -side top -anchor w -ipady 2m -ipadx 10
	pack $p -in .midpart.right -side top -anchor n -fill x
}

########################################################
#  The parameter section for tabudes for large sets random.
########################################################
proc make_lsrand_param p {
	set pp .midpart.right
	catch { destroy $p }
	
	frame $p     -relief groove -bd 2 
	frame $p.l      
	frame $p.r   
	pack $p.l $p.r -in $p -side left -fill x
	
	label $p.l.ln  -text "N:" 
	entry $p.r.n   -textvariable ls_n -width 5 
	label $p.l.lmaxround  -text "Max. rounds:" 
	entry $p.r.maxround   -textvariable ls_maxrounds -width 5 
	label $p.l.lmaxsol  -text "Max. solutions:" 
	entry $p.r.maxsol   -textvariable ls_maxsol -width 5 
	label $p.l.lmaxloops  -text "Max. loops:" 
	entry $p.r.maxloops   -textvariable ls_maxloops -width 5 
	label $p.l.comment  -text "" 
	label $p.r.comment  -text "LLL parameters in t-designs are used also" 

	pack $p.l.ln $p.l.lmaxround $p.l.lmaxsol $p.l.lmaxloops $p.l.comment -in $p.l -side top -anchor w -ipady 2m -ipadx 10
	pack $p.r.n $p.r.maxround $p.r.maxsol $p.r.maxloops $p.r.comment -in $p.r    -side top -anchor w -ipady 2m -ipadx 10
	pack $p -in .midpart.right -side top -anchor n -fill x
}

########################################################
#  The parameter section for makeKM
########################################################
proc make_makekm_param p {
	global font2

	set pp .midpart.right
	catch { destroy $p }
	
	frame $p     -relief groove -bd 2 
	frame $p.l      
	frame $p.r   
	pack $p.l $p.r -in $p -side left -fill x
	
	label $p.l.lt     -text "t:"  
	entry $p.r.t      -textvariable makekm_t -width 5 
	label $p.l.lk     -text "k:" 
	entry $p.r.k      -textvariable makekm_k -width 5 
	label $p.l.lorder -text "order:" 
	entry $p.r.order  -textvariable makekm_order -width 10
	label $p.l.lcomment  -text "" 
	label $p.r.comment  -text "Use file of generators as input !" 

	label $p.l.lgenfile  -text "Gener. filename:" 
	frame $p.r.f
	entry $p.r.f.genfile  -textvariable genfile -width 40 -font $font2
	button $p.r.f.selectfilename  -text "select" \
	  -command { 
	  		set genfile [tk_getOpenFile -filetypes $filetypes2] 
			catch { set genfile [ file tail $genfile] }
			set f [open $genfile r]
			set ff 0
			while { [gets $f line] >= 0 } {
				if [regexp "order=" $line] {
					set ff 1
					regsub {^.*order=} $line {} makekm_order
				}
			} 
			if {$ff==0} {
				set makekm_order "?"
			}
	  } 

	pack $p.r.f.genfile $p.r.f.selectfilename -in $p.r.f -side left \
	     -anchor w -ipady 2m -ipadx 5 -fill x

	pack $p.l.lcomment $p.l.lgenfile $p.l.lorder  $p.l.lk $p.l.lt -in $p.l -side top -anchor w -ipady 2m -ipadx 10
	pack $p.r.comment $p.r.f $p.r.order $p.r.k $p.r.t -in $p.r    -side top -anchor w -ipady 2m -ipadx 10
	pack $p -in .midpart.right -side top -anchor n -fill x
}

########################################################
#  No parameters available
########################################################
proc make_no_param p {
	catch { destroy .midpart.right.param }
	catch { destroy .midpart.right.startbutton }
	frame $p         -relief groove -bd 2 
	label $p.lempty   -text "no parameters" 
	pack $p.lempty -in $p -ipady 2m -ipadx 10
	pack $p -in .midpart.right -side top -anchor n -fill x
}

########################################################
#  Choose problem type (left buttons)
########################################################
proc choose_problemtype type {
	global problemtype 
	
	switch $type {
		designs {set problemtype "t-designs"}
		coverings {set problemtype "coverings"}
		packings {set problemtype "packings"}
		largesets {set problemtype "Large sets of designs"}
	}
	set p .midpart.right.problemtype
	catch { destroy  $p.problemtype }
	label $p.problemtype -text $problemtype 
	pack $p.problemtype -in $p -ipady 2m -ipadx 10
	
	set p .midpart.right
	catch { destroy $p.param }
	
	show_solvers $type
	make_no_param $p.param
}

########################################################
#  Launch solver
########################################################
proc start_solver method {
	global lambda filename logfile problemtype genfile

	global lllc0 lllbeta lllp 
	global lllguess

	global clique_n clique_d clique_w

	global tabudes_seed tabudes_verbose tabudes_maxiter tabudes_neigh
	global tabudes_attr tabudes_volcounts 
	global tabudes_neigh_cover tabudes_neigh_pack 
	
	global ls_n ls_maxsol ls_maxrounds ls_maxloops
	
	global makekm_t makekm_k makekm_order makekm_filename
#
# error handling
#	
	if {($filename == "") && ($method != "makekm") && ($method != "wcliquedirect")} { 
		msg "No filename given !"
		return
	}
	if {($lambda == "") && ($method != "makekm") && ($method != "wcliquedirect")} { 
		msg "lambda not specified !"
		return
	}
	if {($lambda < 0) && ($method != "makekm") && ($method != "wcliquedirect")} { 
		msg "lambda has wrong value !"
		return
	}

# Just for LLL"
	set c0 [ expr $lllc0 * $lambda ]
	
	switch $problemtype {
		"t-designs" {
			switch $method {
				lllwith {
					startsolver "discreta_lll $c0 $lambda bkz $lllbeta $lllp $filename"
				}
				lll {
					startsolver "discreta_lll $c0 $lambda silent bkz $lllbeta $lllp $filename"
				}
				spread {
					startsolver "discreta_spread $lambda $filename"
				}
				mckay {
					startsolver "discreta_mckay mckay.log $lambda <$filename"
				}
				wclique {
					if {$lambda != 1} {
						msg "wclique works only for lambda=1 !"
						return
					}
					startsolver "discreta_clique $filename"
				}
				wcliquedirect {
					startsolver "discreta_cliqued $clique_n $clique_d $clique_w $genfile"
				}
				tabudes {
					startsolver "tabudes discretafile $filename algo tabu maxpen 0 \
					discretaoutput solutions verbose $tabudes_verbose seed $tabudes_seed\
					lambda $lambda attr \"$tabudes_attr\" maxiter $tabudes_maxiter \
					neigh \"$tabudes_neigh\""
				}
				makekm {
					regsub {gen\_} $genfile "KM\_s\_" makekm_filename
					regsub {\.txt} $makekm_filename "\_t$makekm_t\_k$makekm_k\.txt" makekm_filename
					startsolver "discreta_makeKM $genfile $makekm_order $makekm_k \
					$makekm_t $makekm_filename"
					catch { set filename [file tail $makekm_filename] }
				}
				default {
					msg "no solver method selected"
				}
			}
		}
		coverings {
			switch $method {
				tabudes {
					startsolver "tabudes discretafile $filename algo tabu maxpen 0 \
					discretaoutput solutions verbose $tabudes_verbose seed $tabudes_seed\
					lambda $lambda attr \"$tabudes_attr\" maxiter $tabudes_maxiter \
					neigh \"$tabudes_neigh_cover\" coverpenalty cover \
					initsol volumes volcounts \"$tabudes_volcounts\""
				}
				makekm {
					regsub {gen\_} $genfile "KM\_s\_" makekm_filename
					regsub {\.txt} $makekm_filename "\_t$makekm_t\_k$makekm_k\.txt" makekm_filename
					startsolver "discreta_makeKM $genfile $makekm_order $makekm_k \
					$makekm_t $makekm_filename"
					catch { set filename [file tail $makekm_filename] }
				}
				default {
					msg "no solver method selected"
				}
			}
		}
		packings {
			switch $method {
				wclique {
					startsolver "discreta_clique $filename"
				}
				wcliquedirect {
					startsolver "discreta_cliqued $clique_n $clique_d $clique_w $genfile"
				}
				spread {
					startsolver "discreta_spreadpacking 1 $filename"
				}
				lllpacking {
					startsolver "discreta_lllpacking_shell $lllguess $lambda $lllbeta $lllp $filename"
				}
				tabudes {
					startsolver "tabudes discretafile $filename algo tabu maxpen 0 \
					discretaoutput solutions verbose $tabudes_verbose seed $tabudes_seed\
					lambda $lambda attr \"$tabudes_attr\" maxiter $tabudes_maxiter \
					neigh \"$tabudes_neigh_pack\" coverpenalty pack \
					initsol volumes volcounts \"$tabudes_volcounts\""
				}
				makekm {
					regsub {gen\_} $genfile "KM\_s\_" makekm_filename
					regsub {\.txt} $makekm_filename "\_t$makekm_t\_k$makekm_k\.txt" makekm_filename
					startsolver "discreta_makeKM $genfile $makekm_order $makekm_k \
					$makekm_t $makekm_filename"
					catch { set filename [file tail $makekm_filename] }
				}
				default {
					msg "no solver method selected"
				}
			}
		}
		"Large sets of designs" {
			switch $method {
				lsdet {
					startsolver "discreta_ls1 $filename $ls_n $lllc0 $lllbeta $lllp"
				}
				lsrand {
					startsolver "discreta_ls2 $filename $ls_n $ls_maxsol $ls_maxloops \
					$ls_maxrounds $lllc0 $lllbeta $lllp"
				}
				default {
					msg "no solver method selected"
				}
			}
#			msg "$problemtype not implemented"
		}
		default {
			msg "no problem type selected"
		}
	}
}

########################################################
#  start a solver
########################################################
proc startsolver text {
	global processid 
	global proc_handle
	
global filename	
	set p .bottom.r
	catch {
		$p.text insert end "Start program:\n$text ...\n"
		puts "Start program\n$text ...\n"

		set proc_handle [ open |$text w] 
		set processid [ pid $proc_handle ]
		$p.text yview moveto 1
	}
}

########################################################
#  message window
########################################################
proc msg text {
	tk_messageBox -message $text -type ok 
}

########################################################
#  show result
########################################################
proc showresult {} {
	set p .bottom.r

	set result [ exec tail solutions ]
	$p.text insert end "$result\n"
	set result [ exec wc solutions ]
	$p.text insert end "$result\n"
	$p.text yview moveto 1
}

########################################################
#  show Kramer Mesner file
########################################################
proc showkmfile filename {
	set p .bottom.r

	set f [ open $filename r]
	while { [gets $f line] >= 0 } {
		$p.text insert end "$line\n"
	}
	close $f
#	$p.text yview moveto 1
}

########################################################
#  show result
########################################################
proc showmsg text {
	set p .bottom.r

	$p.text insert end "$text\n"
	$p.text yview moveto 1
}

########################################################
#  stop the solver program
########################################################
proc stopit {processid proc_handle} {
		if {$proc_handle>0} {
			exec kill -TERM $processid 
	
			if [eof $proc_handle ] {
				catch {close $proc_handle}
				return
			}
			
			puts "Program with procid $processid stopped !\n"
			.bottom.r.text insert end "Program stopped !\n\n"
		   .bottom.r.text yview moveto 1
		}
}

########################################################
#  stop the solver program
########################################################
proc orbithint filename {
	set p .bottom.r.text
	catch {
		set text [ exec orbithint $filename ]
		$p insert end "$text\n"
	   $p yview moveto 1
	}
}
