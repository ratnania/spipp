##########################################################################
# greenbox tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(greenbox_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(greenbox_ex1 spi_global spi_quadratures spi_mesh spi_blackbox spi_greenbox)
add_test(NAME greenbox_test_1 COMMAND greenbox_ex1)
# ...
##########################################################################
