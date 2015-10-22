##########################################################################
set(files 
	${SPIPP_MESH_DIR}/mesh_def.F90
	${SPIPP_MESH_DIR}/mesh_transform.F90 
	${SPIPP_MESH_DIR}/mesh.F90 )

set(mods 
          ${CMAKE_BINARY_DIR}/modules/spi_mesh_def.mod
	  ${CMAKE_BINARY_DIR}/modules/spi_mesh_transform.mod
          ${CMAKE_BINARY_DIR}/modules/spi_mesh.mod)

add_library(spi_mesh ${files})
target_link_libraries(spi_mesh spi_global)
install(TARGETS spi_mesh DESTINATION lib)
##########################################################################
