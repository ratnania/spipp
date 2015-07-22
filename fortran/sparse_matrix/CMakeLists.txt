##########################################################################
set(files 
	${SPIPP_SPARSE_MATRIX_DIR}/sparse_matrix_def.F90
	${SPIPP_SPARSE_MATRIX_DIR}/sparse_matrix.F90
	${SPIPP_SPARSE_MATRIX_DIR}/sparsekit.F90
	${SPIPP_SPARSE_MATRIX_DIR}/csr_matrix.F90
	${SPIPP_SPARSE_MATRIX_DIR}/bnd_matrix.F90)

set(mods 
	${CMAKE_BINARY_DIR}/modules/SPI_SPARSE_MATRIX.mod
	${CMAKE_BINARY_DIR}/modules/SPI_BND_MATRIX_MODULE.mod
	${CMAKE_BINARY_DIR}/modules/SPI_CSR_MATRIX_MODULE.mod
        ${CMAKE_BINARY_DIR}/modules/spi_utilities.mod)

add_library(spi_matrix_utilities SHARED ${files})
target_link_libraries(spi_matrix_utilities spi_global)
install(TARGETS spi_matrix_utilities DESTINATION lib)
##########################################################################
