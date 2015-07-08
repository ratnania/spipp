##########################################################################
add_subdirectory(bspline)
IF(BUILD_TESTS)
   add_subdirectory(bspline/tests)
ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
#add_subdirectory(pppack)
#IF(BUILD_TESTS)
#   add_subdirectory(bspline/tests)
#ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
add_subdirectory(global)
#IF(BUILD_TESTS)
#   add_subdirectory(global/tests)
#ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
add_subdirectory(quadratures)
#IF(BUILD_TESTS)
#   add_subdirectory(quadratures/tests)
#ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
#add_subdirectory(space)
#IF(BUILD_TESTS)
#   add_subdirectory(space/tests)
#ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
#add_subdirectory(blackbox)
#IF(BUILD_TESTS)
#   add_subdirectory(blackbox/tests)
#ENDIF(BUILD_TESTS)
##########################################################################

