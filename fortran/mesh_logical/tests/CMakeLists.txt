##########################################################################
# mesh_logical tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(mesh_logical_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(mesh_logical_ex1 spi_mesh_logical spi_quadratures spi_global)
add_test(NAME mesh_logical_test_1 COMMAND mesh_logical_ex1)
# ...

##########################################################################
