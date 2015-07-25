!# -*- coding: utf8 -*-
MODULE SPI_BASIS_DEF 
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  USE SPI_MESH_DEF
  IMPLICIT NONE
  
  INTEGER :: n_deriv=SPI_INT_DEFAULT

! ..........................................................        
  TYPE, ABSTRACT, PUBLIC :: DEF_BASIS_ABSTRACT
     INTEGER :: oi_n_order
     INTEGER :: oi_n_max_order
     INTEGER :: oi_nen
     INTEGER :: oi_nderiv
     CHARACTER(LEN=1024) :: dirname 

     CLASS(DEF_QUADRATURE_ABSTRACT), POINTER :: ptr_quad => NULL()
  END TYPE DEF_BASIS_ABSTRACT
! ..........................................................        

! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BASIS_ABSTRACT) :: DEF_BASIS_1D
     REAL(KIND=SPI_RK), DIMENSION(:,:,:)  , ALLOCATABLE :: TestfT_0
     REAL(KIND=SPI_RK), DIMENSION(:,:,:)  , ALLOCATABLE :: TestfT_p, TestfT_pp 
  END TYPE DEF_BASIS_1D
! ..........................................................        

! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BASIS_1D) :: DEF_BASIS_1D_BSPLINE
     CLASS(DEF_MESH_1D_BSPLINE), POINTER :: ptr_mesh => NULL()
  END TYPE DEF_BASIS_1D_BSPLINE
! ..........................................................        

! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BASIS_1D) :: DEF_BASIS_1D_FOURIER
     CLASS(DEF_MESH_1D_FOURIER), POINTER :: ptr_mesh => NULL()
  END TYPE DEF_BASIS_1D_FOURIER
! ..........................................................        

! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BASIS_1D) :: DEF_BASIS_1D_HBEZIER
     CLASS(DEF_MESH_1D_HBEZIER), POINTER :: ptr_mesh => NULL()
  END TYPE DEF_BASIS_1D_HBEZIER
! ..........................................................        

END MODULE SPI_BASIS_DEF
