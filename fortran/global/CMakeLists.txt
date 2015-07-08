##########################################################################
set(self_files 
	${SPIPP_GLOBAL_DIR}/global_def.F90
	${SPIPP_GLOBAL_DIR}/global.F90 )

set(self_mods 
          ${CMAKE_BINARY_DIR}/modules/global_def.mod
          ${CMAKE_BINARY_DIR}/modules/global.mod)

IF(APPLE)
  add_library(self STATIC ${self_files})
  target_link_libraries(self   ) # TODO add linked libs 
ELSE()
  add_library(self SHARED ${self_files})
  target_link_libraries(self  ) # TODO add linked libs
ENDIF()
install(TARGETS self DESTINATION lib)
install(FILES ${self_mods} DESTINATION include)
##########################################################################
