##########################################################################
# bspline directory
##########################################################################
set(files 
	bspline.F90 )
set(mods 
	bspline.mod  
	bsp.mod  
	crv.mod  
	iga.mod  
	srf.mod  
	vol.mod)
#
add_library(bspline SHARED ${files})
install(TARGETS bspline DESTINATION lib)
##########################################################################
