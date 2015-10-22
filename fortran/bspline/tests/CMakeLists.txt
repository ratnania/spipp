##########################################################################
# splines tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(bspline_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(bspline_ex1 spi_bspline)
add_test(NAME splines_test_1 COMMAND bspline_ex1)
# ...
# ... ex2
ADD_EXECUTABLE(bspline_ex2 ex2.F90)
TARGET_LINK_LIBRARIES(bspline_ex2 spi_bspline)
add_test(NAME splines_test_2 COMMAND bspline_ex2)
# ...
# ... ex3
ADD_EXECUTABLE(bspline_ex3 ex3.F90)
TARGET_LINK_LIBRARIES(bspline_ex3 spi_bspline)
add_test(NAME splines_test_3 COMMAND bspline_ex3)
# ...
##########################################################################
