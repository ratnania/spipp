##########################################################################
# bspline directory
##########################################################################
set(files 
	banfac.F90
	banslv.F90
	bchfac.F90
	bchslv.F90
	bsplpp.F90
	bsplvb.F90
	bsplvb_vect.F90
	bsplvd.F90
	bsplvd_vect.F90
	bspp2d.F90
	bvalue.F90
	chol1d.F90
	colloc.F90
	colpnt.F90
	cubspl.F90
	cwidth.F90
	difequ.F90
	eqblok.F90
	interv.F90
	knots.F90
	l2appr.F90
	l2err.F90
	l2knts.F90
	l2main.F90
	newnot.F90
	newnotfake.F90
	pppack.F90
	ppvalu.F90
	putit.F90
	round.F90
	setupq.F90
	slvblktext.F90
	spli2d.F90
	splint.F90
	splint2d.F90
	splopt.F90)

set(mods 
	banfac.mod
	banslv.mod
	bchfac.mod
	bchslv.mod
	bsplpp.mod
	bsplvb.mod
	bsplvb_vect.mod
	bsplvd.mod
	bsplvd_vect.mod
	bspp2d.mod
	bvalue.mod
	chol1d.mod
	colloc.mod
	colpnt.mod
	cubspl.mod
	cwidth.mod
	difequ.mod
	eqblok.mod
	interv.mod
	knots.mod
	l2appr.mod
	l2err.mod
	l2knts.mod
	l2main.mod
	newnot.mod
	newnotfake.mod
	pppack.mod
	ppvalu.mod
	putit.mod
	round.mod
	setupq.mod
	slvblktext.mod
	spli2d.mod
	splint.mod
	splint2d.mod
	splopt.mod)
#
add_library(pppack SHARED ${files})
install(TARGETS pppack DESTINATION lib)
##########################################################################
