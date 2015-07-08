##########################################################################
# splines tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(ex1 ex1.F90)
TARGET_LINK_LIBRARIES(ex1 bspline)
add_test(NAME splines_test_1 COMMAND ex1)
# ...
# ... ex2
ADD_EXECUTABLE(ex2 ex2.F90)
TARGET_LINK_LIBRARIES(ex2 bspline)
add_test(NAME splines_test_2 COMMAND ex2)
# ...
# ... ex3
ADD_EXECUTABLE(ex3 ex3.F90)
TARGET_LINK_LIBRARIES(ex3 bspline)
add_test(NAME splines_test_3 COMMAND ex3)
# ...
##########################################################################
