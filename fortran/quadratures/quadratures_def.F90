!# -*- coding: utf8 -*-
MODULE SPI_QUADRATURES_DEF
  USE SPI_GLOBAL_DEF
  IMPLICIT NONE

! ..........................................................        
  TYPE, ABSTRACT, PUBLIC :: DEF_QUADRATURE_ABSTRACT
     INTEGER :: oi_n_points
     INTEGER :: oi_type
     CHARACTER(LEN=1024) :: dirname 

     REAL(KIND=SPI_RK), DIMENSION(:)  , ALLOCATABLE :: weights 
     REAL(KIND=SPI_RK), DIMENSION(:), ALLOCATABLE :: points
  END TYPE DEF_QUADRATURE_ABSTRACT

  TYPE, PUBLIC, EXTENDS(DEF_QUADRATURE_ABSTRACT) :: DEF_QUADRATURE_1D
  END TYPE DEF_QUADRATURE_1D
! ..........................................................        

! ..........................................................        

END MODULE SPI_QUADRATURES_DEF
