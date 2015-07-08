##########################################################################
set(files 
	${SPIPP_GLOBAL_DIR}/global_def.F90 
	${SPIPP_GLOBAL_DIR}/global.F90 )
set(mods 
	${CMAKE_BINARY_DIR}/spi_global_def.mod  
	${CMAKE_BINARY_DIR}/spi_global.mod)
#
add_library(spi_global SHARED ${files})
install(TARGETS spi_global DESTINATION lib)
##########################################################################

