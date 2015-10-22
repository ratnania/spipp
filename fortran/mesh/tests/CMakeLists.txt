##########################################################################
# mesh tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(mesh_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(mesh_ex1 spi_mesh spi_quadratures spi_global)
add_test(NAME mesh_test_1 COMMAND mesh_ex1)
# ...
##########################################################################
