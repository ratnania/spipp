##########################################################################
set(files 
	${SPIPP_UTILITIES_DIR}/quicksort.F90)

set(mods 
          ${CMAKE_BINARY_DIR}/modules/spi_qsort.mod)

add_library(spi_utilities SHARED ${files})
target_link_libraries(spi_utilities spi_global)
install(TARGETS spi_utilities DESTINATION lib)
##########################################################################
