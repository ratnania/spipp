##########################################################################
# blackbox tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(blackbox_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(blackbox_ex1 spi_global spi_quadratures spi_basis spi_blackbox)
add_test(NAME blackbox_test_1 COMMAND blackbox_ex1)
# ...
##########################################################################
