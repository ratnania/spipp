##########################################################################
set(files 
	${SPIPP_MATRIX_DIR}/matrix_def.F90
	${SPIPP_MATRIX_DIR}/matrix.F90)

set(mods 
          ${CMAKE_BINARY_DIR}/modules/spi_matrix_def.mod
          ${CMAKE_BINARY_DIR}/modules/spi_matrix.mod)

add_library(spi_matrix SHARED ${files})
target_link_libraries(spi_matrix spi_global spi_sparse_matrix ${LAPACK_LIBRARIES})
install(TARGETS spi_matrix DESTINATION lib)
##########################################################################
