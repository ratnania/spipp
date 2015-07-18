##########################################################################
# space directory
##########################################################################
set(space_files 
	${JOREK_SPACE_DIR}/space_def.F90 
	${JOREK_SPACE_DIR}/space.F90)

set(space_mods 
          ${CMAKE_BINARY_DIR}/modules/space_def.mod
          ${CMAKE_BINARY_DIR}/modules/space.mod)

IF(APPLE)
  add_library(space STATIC ${space_files})
  target_link_libraries(space  ${SPM_LIBRARIES} jorek_glob jorek_param quadratures mesh basis blackbox greenbox numbering graph boundary_conditions )
ELSE()
  add_library(space SHARED ${space_files})
  target_link_libraries(space  ${SPM_LIBRARIES} jorek_glob jorek_param quadratures mesh basis blackbox greenbox numbering graph boundary_conditions ${MPI_LIBRARY})
ENDIF()
install(TARGETS space DESTINATION lib)
install(FILES ${space_mods} DESTINATION include)
##########################################################################
