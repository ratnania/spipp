##########################################################################
# mesh tests directory
##########################################################################

# ... ex1
ADD_EXECUTABLE(mesh_ex1 ex1.F90)
TARGET_LINK_LIBRARIES(mesh_ex1 spi_mesh spi_mesh_logical spi_quadratures spi_global)
add_test(NAME mesh_test_1 COMMAND mesh_ex1)
# ...

# ... ex2
#ADD_EXECUTABLE(mesh_ex2 ex2.F90)
#TARGET_LINK_LIBRARIES(mesh_ex2 spi_mesh spi_mesh_logical spi_quadratures spi_global)
#add_test(NAME mesh_test_2 COMMAND mesh_ex2)
# ...
##########################################################################
