##########################################################################
# greenbox directory
##########################################################################
set(files 
	${SPIPP_GREENBOX_DIR}/greenbox_def.F90
	${SPIPP_GREENBOX_DIR}/greenbox.F90)

set(mods 
          ${CMAKE_BINARY_DIR}/modules/greenbox_def.mod
          ${CMAKE_BINARY_DIR}/modules/greenbox.mod
)
add_library(spi_greenbox SHARED ${files})
target_link_libraries(spi_greenbox spi_global spi_quadratures spi_mesh)
install(TARGETS spi_greenbox DESTINATION lib)
##########################################################################
