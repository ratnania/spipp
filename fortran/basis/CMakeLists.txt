##########################################################################
set(files 
	${SPIPP_BASIS_DIR}/basis_def.F90
        ${SPIPP_BASIS_DIR}/bsplines.F90
        ${SPIPP_BASIS_DIR}/fourier_mods.F90
        ${SPIPP_BASIS_DIR}/fourier.F90
	${SPIPP_BASIS_DIR}/basis.F90 )

set(mods 
          ${CMAKE_BINARY_DIR}/modules/basis_def.mod
          ${CMAKE_BINARY_DIR}/modules/spi_basis_bsplines.mod
          ${CMAKE_BINARY_DIR}/modules/spi_fourier_mods.mod
          ${CMAKE_BINARY_DIR}/modules/spi_basis_fourier.mod
          ${CMAKE_BINARY_DIR}/modules/basis.mod)

add_library(spi_basis SHARED ${files})
target_link_libraries(spi_basis spi_global spi_mesh)
install(TARGETS spi_basis DESTINATION lib)
##########################################################################
