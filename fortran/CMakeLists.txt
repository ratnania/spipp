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
##########################################################################

##########################################################################
add_subdirectory(utilities)
#IF(BUILD_TESTS)
#   add_subdirectory(utilities/tests)
#ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
add_subdirectory(quadratures)
IF(BUILD_TESTS)
   add_subdirectory(quadratures/tests)
ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
add_subdirectory(mesh)
IF(BUILD_TESTS)
   add_subdirectory(mesh/tests)
ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
add_subdirectory(basis)
IF(BUILD_TESTS)
   add_subdirectory(basis/tests)
ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
add_subdirectory(blackbox)
IF(BUILD_TESTS)
   add_subdirectory(blackbox/tests)
ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
add_subdirectory(greenbox)
#IF(BUILD_TESTS)
#   add_subdirectory(greenbox/tests)
#ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
add_subdirectory(numbering)
IF(BUILD_TESTS)
   add_subdirectory(numbering/tests)
ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
add_subdirectory(sparse_matrix)
IF(BUILD_TESTS)
   add_subdirectory(sparse_matrix/tests)
ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
add_subdirectory(matrix)
#IF(BUILD_TESTS)
#   add_subdirectory(matrix/tests)
#ENDIF(BUILD_TESTS)
##########################################################################

##########################################################################
add_subdirectory(space)
#IF(BUILD_TESTS)
#   add_subdirectory(space/tests)
#ENDIF(BUILD_TESTS)
##########################################################################

