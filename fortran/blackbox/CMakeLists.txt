##########################################################################
set(files 
	  ${SPIPP_BLACKBOX_DIR}/blackbox_def.F90
	  ${SPIPP_BLACKBOX_DIR}/blackbox_bsplines.F90
#	  ${SPIPP_BLACKBOX_DIR}/blackbox_fourier.F90
#	  ${SPIPP_BLACKBOX_DIR}/blackbox_hbezier.F90
	  ${SPIPP_BLACKBOX_DIR}/blackbox.F90 )

set(mods 
          ${CMAKE_BINARY_DIR}/modules/spi_blackbox_def.mod
          ${CMAKE_BINARY_DIR}/modules/spi_blackbox_bsplines.mod
#          ${CMAKE_BINARY_DIR}/modules/spi_blackbox_fourier.mod
#          ${CMAKE_BINARY_DIR}/modules/spi_blackbox_hbezier.mod
          ${CMAKE_BINARY_DIR}/modules/spi_blackbox.mod)

add_library(spi_blackbox ${files})
target_link_libraries(spi_blackbox spi_global spi_quadratures spi_mesh spi_basis)
install(TARGETS spi_blackbox DESTINATION lib)
##########################################################################
