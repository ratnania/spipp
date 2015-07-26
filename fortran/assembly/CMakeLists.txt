##########################################################################
set(files 
	${SPIPP_ASSEMBLY_DIR}/assembly_def.F90
	${SPIPP_ASSEMBLY_DIR}/assembly_1d_bspline.F90
	${SPIPP_ASSEMBLY_DIR}/assembly.F90)

set(mods 
          ${CMAKE_BINARY_DIR}/modules/spi_assembly_def.mod
          ${CMAKE_BINARY_DIR}/modules/spi_assembly_1d_bspline.mod
          ${CMAKE_BINARY_DIR}/modules/spi_assembly.mod)

add_library(spi_assembly SHARED ${files})
target_link_libraries(spi_assembly spi_global spi_greenbox spi_matrix spi_mesh spi_space)
install(TARGETS spi_assembly DESTINATION lib)
##########################################################################
