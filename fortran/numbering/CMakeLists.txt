##########################################################################
# numbering directory
##########################################################################
set(files 
	${SPIPP_NUMBERING_DIR}/numbering_def.F90
	${SPIPP_NUMBERING_DIR}/numbering.F90)

set(mods 
          ${CMAKE_BINARY_DIR}/modules/numbering_def.mod
          ${CMAKE_BINARY_DIR}/modules/numbering.mod
)
add_library(spi_numbering ${files})
target_link_libraries(spi_numbering spi_global spi_quadratures spi_mesh)
install(TARGETS spi_numbering DESTINATION lib)
##########################################################################
