##########################################################################
# basis tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(basis_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(basis_ex1 spi_global spi_mesh spi_quadratures spi_basis)
add_test(NAME basis_test_1 COMMAND basis_ex1)
# ...

# ... ex2
ADD_EXECUTABLE(basis_ex2 ex2.F90)
TARGET_LINK_LIBRARIES(basis_ex2 spi_global spi_mesh spi_quadratures spi_basis)
add_test(NAME basis_test_2 COMMAND basis_ex2)
# ...

# ... ex3
ADD_EXECUTABLE(basis_ex3 ex3.F90)
TARGET_LINK_LIBRARIES(basis_ex3 spi_global spi_mesh spi_quadratures spi_basis)
add_test(NAME basis_test_3 COMMAND basis_ex3)
# ...
##########################################################################
