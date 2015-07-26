!# -*- coding: utf8 -*-
MODULE SPI_ASSEMBLY_DEF 
  USE SPI_GLOBAL_DEF
  IMPLICIT NONE
  
  LOGICAL :: ml_assembly_initialized

  INTEGER, DIMENSION(:)  , POINTER :: mpi_ltog_rhsrows
  INTEGER, DIMENSION(:,:), POINTER :: mpi_ltog_rows
  INTEGER, DIMENSION(:,:), POINTER :: mpi_ltog_cols

END MODULE SPI_ASSEMBLY_DEF
