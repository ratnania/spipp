!# -*- coding: utf8 -*-
MODULE SPI_MESH_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  IMPLICIT NONE

  ! .........................................................
!  TYPE, PUBLIC :: DEF_ELEMENT
!     INTEGER :: oi_nen
!     REAL(KIND=SPI_RK) , DIMENSION(:),   ALLOCATABLE :: opr_knots 
!
!     INTEGER :: oi_n_vtex
!     INTEGER :: oi_n_bnet
!     INTEGER :: oi_n_order
!     INTEGER :: oi_loc_id   ! local number to the current patch
!
!     INTEGER, DIMENSION(:), ALLOCATABLE    :: opi_LocToGlob  
!  END TYPE DEF_ELEMENT
  ! .........................................................

  ! .........................................................
  TYPE, ABSTRACT, PUBLIC :: DEF_MESH_ABSTRACT
     INTEGER :: n_elements

     integer :: n 
     integer :: p
     INTEGER :: n_dim

     REAL(SPI_RK), dimension (:), allocatable :: opr_grid
     REAL(SPI_RK), dimension (:,:), allocatable :: opr_points
     CLASS(DEF_QUADRATURE_ABSTRACT), POINTER :: ptr_quad => NULL()
  END TYPE DEF_MESH_ABSTRACT
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_ABSTRACT) :: DEF_MESH_1D_BSPLINE
     INTEGER :: type_bc

     !> KNOT VECTOR FOR EACH DIRECTION
     real(SPI_RK), dimension (:), ALLOCATABLE :: knots
     REAL(SPI_RK), DIMENSION (:,:), ALLOCATABLE :: control_points
  END TYPE DEF_MESH_1D_BSPLINE
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_ABSTRACT) :: DEF_MESH_2D_BSPLINE
     TYPE(DEF_MESH_1D_BSPLINE) :: mesh_u
     TYPE(DEF_MESH_1D_BSPLINE) :: mesh_v

     REAL(SPI_RK), DIMENSION (:,:,:), ALLOCATABLE :: control_points
  END TYPE DEF_MESH_2D_BSPLINE
  ! .........................................................

END MODULE SPI_MESH_DEF
