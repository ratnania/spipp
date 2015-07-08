##########################################################################
set(files 
	  ${SPIPP_BLACKBOX_DIR}/blackbox_def.F90
	  ${SPIPP_BLACKBOX_DIR}/blackbox.F90 )

set(mods 
          ${CMAKE_BINARY_DIR}/modules/spi_blackbox_def.mod
          ${CMAKE_BINARY_DIR}/modules/spi_blackbox.mod)

add_library(spi_blackbox SHARED ${files})
target_link_libraries(spi_blackbox spi_global spi_quadratures)
install(TARGETS spi_blackbox DESTINATION lib)
##########################################################################
