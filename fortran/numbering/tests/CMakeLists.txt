##########################################################################
# numbering tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(numbering_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(numbering_ex1 spi_global spi_quadratures spi_mesh spi_numbering)
add_test(NAME numbering_test_1 COMMAND numbering_ex1)
# ...
##########################################################################
