##########################################################################
# space tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(space_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(space_ex1 spi_global spi_quadratures spi_mesh spi_space)
add_test(NAME space_test_1 COMMAND space_ex1)
# ...
##########################################################################
