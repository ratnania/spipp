##########################################################################
# quadratures tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(quadratures_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(quadratures_ex1 spi_global spi_quadratures)
add_test(NAME quadratures_test_1 COMMAND quadratures_ex1)
# ...
##########################################################################
