##########################################################################
set(files 
	  ${SPIPP_QUADRATURES_DIR}/quadratures_def.F90
	  ${SPIPP_QUADRATURES_DIR}/ptpdgl.F90
	  ${SPIPP_QUADRATURES_DIR}/quadratures_line.F90 
	  ${SPIPP_QUADRATURES_DIR}/quadratures.F90 )

set(mods 
          ${CMAKE_BINARY_DIR}/spi_quadratures_def.mod
          ${CMAKE_BINARY_DIR}/spi_quadratures_line.mod
          ${CMAKE_BINARY_DIR}/spi_quadratures.mod)
#
add_library(spi_quadratures ${files})
target_link_libraries(spi_quadratures spi_global)
install(TARGETS spi_quadratures DESTINATION lib)
##########################################################################
