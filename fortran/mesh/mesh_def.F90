!# -*- coding: utf8 -*-
MODULE SPI_MESH_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_MESH_LOGICAL_DEF
  IMPLICIT NONE

  ! .........................................................
  TYPE, ABSTRACT, PUBLIC :: DEF_MESH_CONTROL_POINTS_ABSTRACT
  END TYPE DEF_MESH_CONTROL_POINTS_ABSTRACT
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_CONTROL_POINTS_ABSTRACT) :: DEF_MESH_CONTROL_POINTS_1D
     REAL(SPI_RK), DIMENSION (:,:), ALLOCATABLE :: coef
  END TYPE DEF_MESH_CONTROL_POINTS_1D
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_CONTROL_POINTS_ABSTRACT) :: DEF_MESH_CONTROL_POINTS_2D
     REAL(SPI_RK), DIMENSION (:,:,:), ALLOCATABLE :: coef
  END TYPE DEF_MESH_CONTROL_POINTS_2D
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_CONTROL_POINTS_ABSTRACT) :: DEF_MESH_CONTROL_POINTS_3D
     REAL(SPI_RK), DIMENSION (:,:,:,:), ALLOCATABLE :: coef
  END TYPE DEF_MESH_CONTROL_POINTS_3D
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC :: DEF_MESH_LOGICAL_BSPLINE_POINTER
     CLASS(DEF_MESH_LOGICAL_BSPLINE), POINTER :: ptr_mesh
  END TYPE DEF_MESH_LOGICAL_BSPLINE_POINTER
  ! .........................................................

  ! .........................................................
  TYPE, ABSTRACT, PUBLIC :: DEF_MESH_BSPLINE_ABSTRACT
     INTEGER :: n_elements = 0
     INTEGER :: n_dim      = 0 
     INTEGER :: d_dim      = 0
     INTEGER :: c_dim      = 0

     TYPE(DEF_MESH_LOGICAL_BSPLINE_POINTER), DIMENSION(:), ALLOCATABLE :: logicals
  END TYPE DEF_MESH_BSPLINE_ABSTRACT
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_BSPLINE_ABSTRACT) :: DEF_MESH_BSPLINE_1D
     ! c_dim= 1 
     TYPE(DEF_MESH_CONTROL_POINTS_1D), DIMENSION(:), ALLOCATABLE :: control_points 
  END TYPE DEF_MESH_BSPLINE_1D
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_BSPLINE_ABSTRACT) :: DEF_MESH_BSPLINE_2D
     ! c_dim= 1 
     TYPE(DEF_MESH_CONTROL_POINTS_2D), DIMENSION(:), ALLOCATABLE :: control_points
  END TYPE DEF_MESH_BSPLINE_2D
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_BSPLINE_ABSTRACT) :: DEF_MESH_BSPLINE_3D
     ! c_dim= 1 
     TYPE(DEF_MESH_CONTROL_POINTS_3D), DIMENSION(:), ALLOCATABLE :: control_points
  END TYPE DEF_MESH_BSPLINE_3D
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_BSPLINE_ABSTRACT) :: DEF_MESH_BSPLINE_1D1D
     ! c_dim= 2 
     TYPE(DEF_MESH_CONTROL_POINTS_1D), DIMENSION(:), ALLOCATABLE :: control_points 
  END TYPE DEF_MESH_BSPLINE_1D1D
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_BSPLINE_ABSTRACT) :: DEF_MESH_1D1D1D_BSPLINE
     ! c_dim= 3 
     TYPE(DEF_MESH_CONTROL_POINTS_1D), DIMENSION(:), ALLOCATABLE :: control_points 
  END TYPE DEF_MESH_1D1D1D_BSPLINE
  ! .........................................................

END MODULE SPI_MESH_DEF
