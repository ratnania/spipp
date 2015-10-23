!# -*- coding: utf8 -*-
MODULE SPI_MESH_LOGICAL_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  IMPLICIT NONE

  ! .........................................................
  TYPE, ABSTRACT, PUBLIC :: DEF_MESH_LOGICAL_ABSTRACT
     INTEGER :: n_elements = 0

     REAL(SPI_RK), DIMENSION (:)  , ALLOCATABLE :: grid
     REAL(SPI_RK), DIMENSION (:,:), ALLOCATABLE :: points

     CLASS(DEF_QUADRATURE_ABSTRACT), POINTER :: ptr_quad => NULL()
  END TYPE DEF_MESH_LOGICAL_ABSTRACT
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_LOGICAL_ABSTRACT) :: DEF_MESH_LOGICAL_BSPLINE
     integer :: n 
     integer :: p
     INTEGER :: type_bc

     REAL(SPI_RK), DIMENSION (:)  , ALLOCATABLE :: knots
  END TYPE DEF_MESH_LOGICAL_BSPLINE
  ! .........................................................

END MODULE SPI_MESH_LOGICAL_DEF
