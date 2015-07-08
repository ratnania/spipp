##########################################################################
set(files 
	${SPIPP_BSPLINE_DIR}/bspline.F90 )
set(mods 
	${CMAKE_BINARY_DIR}/bspline.mod  
	${CMAKE_BINARY_DIR}/bsp.mod  
	${CMAKE_BINARY_DIR}/crv.mod  
	${CMAKE_BINARY_DIR}/iga.mod  
	${CMAKE_BINARY_DIR}/srf.mod  
	${CMAKE_BINARY_DIR}/vol.mod)
#
add_library(spi_bspline SHARED ${files})
install(TARGETS spi_bspline DESTINATION lib)
##########################################################################
