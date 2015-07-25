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
     INTEGER :: oi_n_dim    
     INTEGER :: oi_n_elmts
     INTEGER :: oi_n_vtex_per_elmt

     !> NUBMER OF NON VANISHING BASIS 
     INTEGER :: oi_nen
     !> NUMBER NON ZERO ELEMENTS PER DIRECTION
     INTEGER :: oi_nnp
     !> TYPE BOUNDARY CONDITION
     INTEGER :: oi_type_bc

     CLASS(DEF_QUADRATURE_1D), POINTER :: ptr_quad => NULL()
  END TYPE DEF_MESH_ABSTRACT
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_ABSTRACT) :: DEF_MESH_1D
  END TYPE DEF_MESH_1D
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_1D) :: DEF_MESH_1D_HBEZIER
  END TYPE DEF_MESH_1D_HBEZIER
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_1D) :: DEF_MESH_1D_FOURIER
     INTEGER :: oi_n_plane
  END TYPE DEF_MESH_1D_FOURIER
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_1D) :: DEF_MESH_1D_BSPLINE
      integer		:: oi_n 
      integer		:: oi_p	

      !> KNOT VECTOR FOR EACH DIRECTION
      real(SPI_RK), dimension ( : ), pointer	:: opr_knot

  END TYPE DEF_MESH_1D_BSPLINE
  ! .........................................................

END MODULE SPI_MESH_DEF
