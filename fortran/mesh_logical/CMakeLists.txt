##########################################################################
set(files 
	${SPIPP_MESH_LOGICAL_DIR}/mesh_logical_def.F90
	${SPIPP_MESH_LOGICAL_DIR}/mesh_logical_bspline.F90
	${SPIPP_MESH_LOGICAL_DIR}/mesh_logical.F90)

set(mods 
          ${CMAKE_BINARY_DIR}/modules/spi_mesh_logical_def.mod
          ${CMAKE_BINARY_DIR}/modules/spi_mesh_logical_bspline.mod
          ${CMAKE_BINARY_DIR}/modules/spi_mesh_logical.mod)

add_library(spi_mesh_logical ${files})
target_link_libraries(spi_mesh_logical spi_global)
install(TARGETS spi_mesh_logical DESTINATION lib)
##########################################################################
