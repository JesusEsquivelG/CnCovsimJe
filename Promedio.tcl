set tcl_precision 4
set archivo [lindex $argv 0]
set filec [open $archivo RDONLY]
set num_nodx [lindex $argv 1]
set simulacionesx [lindex $argv 2]

scan $num_nodx "%d" num_nod
scan $simulacionesx "%d" simulaciones

#puts $num_nod
#puts $simulaciones

set cantidades_enlaces_procesadas {}
while {[eof $filec]==0} {
    set linea [gets $filec]
    if {[string index $linea 0]!={*} && [string index $linea 0]!={ }} {
    set enlaces [lindex $linea 0]
	if {[lsearch -exact $cantidades_enlaces_procesadas $enlaces]==-1} {
	    set nodos_con_($enlaces) [lindex $linea 1]
	    set cantidades_enlaces_procesadas [lappend cantidades_enlaces_procesadas $enlaces]
	} else {
	    set nodos_con_($enlaces) [lappend nodos_con_($enlaces) [lindex $linea 1]]
	}
    }
}
close $filec
set mayor 0
for {set i 0} {$i<[expr [llength $cantidades_enlaces_procesadas]-1]} {incr i} {
	if {[lindex $cantidades_enlaces_procesadas $i]>=$mayor} {
		set mayor [lindex $cantidades_enlaces_procesadas $i]
	}
}
set suma_prom 0
for {set i -1} {$i<=$mayor} {incr i} {
	    set cant_enlaces $i
	if {[info exists nodos_con_($cant_enlaces)]==1} {
    set suma 0
    for {set j 0} {$j<[llength $nodos_con_($cant_enlaces)]} {incr j} {
    set suma [expr $suma + [lindex $nodos_con_($cant_enlaces) $j]]
    }
    set promedio_de_nodos [expr 1.0*$suma/$simulaciones]
    puts "$cant_enlaces	[expr 1.0*$promedio_de_nodos/$num_nod]"
    set suma_prom [expr $suma_prom+$promedio_de_nodos]
}
}



