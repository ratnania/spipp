!# -*- coding: utf8 -*-
MODULE SPI_QUADRATURES_DEF
  USE SPI_GLOBAL_DEF
  IMPLICIT NONE

! ..........................................................        
  TYPE, ABSTRACT, PUBLIC :: DEF_QUADRATURE_ABSTRACT
     INTEGER :: oi_n_points
     INTEGER :: oi_type
     INTEGER :: oi_dim
     CHARACTER(LEN=1024) :: dirname 

     REAL(KIND=SPI_RK), DIMENSION(:)  , ALLOCATABLE :: opr_weights 
     REAL(KIND=SPI_RK), DIMENSION(:,:), ALLOCATABLE :: opr_points
  END TYPE DEF_QUADRATURE_ABSTRACT

  TYPE, PUBLIC, EXTENDS(DEF_QUADRATURE_ABSTRACT) :: DEF_QUADRATURE_1D
  END TYPE DEF_QUADRATURE_1D
! ..........................................................        

! ..........................................................        

END MODULE SPI_QUADRATURES_DEF
