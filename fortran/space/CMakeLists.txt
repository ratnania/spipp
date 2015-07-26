##########################################################################
set(files 
	${SPIPP_SPACE_DIR}/space_def.F90
	${SPIPP_SPACE_DIR}/space.F90 )

set(mods 
          ${CMAKE_BINARY_DIR}/modules/spi_space_def.mod
          ${CMAKE_BINARY_DIR}/modules/spi_space.mod)

add_library(spi_space SHARED ${files})
target_link_libraries(spi_space spi_global spi_quadratures spi_mesh spi_basis spi_blackbox spi_numbering)
install(TARGETS spi_space DESTINATION lib)
##########################################################################
