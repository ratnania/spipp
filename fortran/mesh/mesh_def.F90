!# -*- coding: utf8 -*-
MODULE SPI_MESH_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF

  ! .........................................................
  TYPE, PUBLIC :: DEF_ELEMENT
     INTEGER :: oi_nen
     REAL(KIND=SPI_RK) , DIMENSION(:),   ALLOCATABLE :: opr_knots 

     INTEGER :: oi_n_vtex
     INTEGER :: oi_n_bnet
     INTEGER :: oi_n_order
     INTEGER :: oi_loc_id   ! local number to the current patch

     INTEGER, DIMENSION(:), ALLOCATABLE    :: opi_LocToGlob  
  END TYPE DEF_ELEMENT
  ! .........................................................

  ! .........................................................
  TYPE, ABSTRACT, PUBLIC :: DEF_MESH_ABSTRACT
     INTEGER :: oi_n_dim    
     INTEGER :: oi_n_Elmts
     INTEGER :: oi_n_vtex_per_elmt

     INTEGER :: oi_nen
  END TYPE DEF_MESH_ABSTRACT
  ! .........................................................

  ! .........................................................
  TYPE, PUBLIC, EXTENDS(DEF_MESH_ABSTRACT) :: DEF_MESH_1D
     INTEGER :: oi_n_plane

     CLASS(DEF_QUADRATURE_1D), POINTER :: ptr_quad => NULL()
  END TYPE DEF_MESH_1D
  ! .........................................................

END MODULE SPI_MESH_DEF
