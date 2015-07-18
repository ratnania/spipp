!# -*- coding: utf8 -*-
MODULE SPI_BASIS_DEF 
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  IMPLICIT NONE
  
  INTEGER :: n_deriv=AnyI

! ..........................................................        
  TYPE, ABSTRACT, PUBLIC :: DEF_BASIS_ABSTRACT
     INTEGER :: oi_n_order
     INTEGER :: oi_n_max_order
     INTEGER :: oi_nen
     INTEGER :: oi_nderiv
     INTEGER :: oi_type
     CHARACTER(LEN=1024) :: dirname 

     CLASS(DEF_QUADRATURE_ABSTRACT), POINTER :: ptr_quad => NULL()
  END TYPE DEF_BASIS_ABSTRACT
! ..........................................................        

! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BASIS_ABSTRACT) :: DEF_BASIS_1D
     REAL(KIND=RK), DIMENSION(:,:,:)  , ALLOCATABLE :: TestfT_0
     REAL(KIND=RK), DIMENSION(:,:,:)  , ALLOCATABLE :: TestfT_p, TestfT_pp 
  END TYPE DEF_BASIS_1D
! ..........................................................        

  TYPE, PUBLIC, EXTENDS(DEF_BASIS_1D) :: DEF_BASIS_1D_BSPLINES
  END TYPE DEF_BASIS_1D_BSPLINES

  TYPE, PUBLIC, EXTENDS(DEF_BASIS_1D) :: DEF_BASIS_1D_FOURIER
  END TYPE DEF_BASIS_1D_FOURIER

END MODULE SPI_BASIS_DEF
