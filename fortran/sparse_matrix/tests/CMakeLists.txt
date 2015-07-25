##########################################################################
# sparse_matrix tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(sparse_matrix_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(sparse_matrix_ex1 spi_global spi_sparse_matrix ${LAPACK_LIBRARIES})
add_test(NAME sparse_matrix_test_1 COMMAND sparse_matrix_ex1)
# ...
# ... ex2
ADD_EXECUTABLE(sparse_matrix_ex2 ex2.F90)
TARGET_LINK_LIBRARIES(sparse_matrix_ex2 spi_global spi_sparse_matrix  spi_quadratures spi_mesh spi_numbering ${LAPACK_LIBRARIES})
add_test(NAME sparse_matrix_test_2 COMMAND sparse_matrix_ex2)
# ...
##########################################################################
