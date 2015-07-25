##########################################################################
# matrix tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(matrix_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(matrix_ex1 spi_global spi_quadratures spi_mesh spi_basis spi_numbering spi_blackbox spi_greenbox  spi_sparse_matrix spi_matrix ${LAPACK_LIBRARIES})
add_test(NAME matrix_test_1 COMMAND matrix_ex1)
# ...
##########################################################################
