##########################################################################
add_subdirectory(bspline)
#add_subdirectory(pppack)
IF(BUILD_TESTS)
   add_subdirectory(bspline/tests)
ENDIF(BUILD_TESTS)
##########################################################################
