!# -*- coding: utf8 -*-
MODULE SPI_GREENBOX_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF

  ! ..........................................................        
  TYPE, ABSTRACT, PUBLIC :: DEF_GREENBOX_ABSTRACT
     INTEGER :: oi_nvar
  END TYPE DEF_GREENBOX_ABSTRACT
  ! ..........................................................        

  ! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_GREENBOX_ABSTRACT) :: DEF_GREENBOX_1D
     REAL(KIND=SPI_RK), DIMENSION(:,:), ALLOCATABLE :: VarN_0  

     ! ... Logical Derivatives
     REAL(KIND=SPI_RK), DIMENSION(:,:), ALLOCATABLE :: VarN_s1
     REAL(KIND=SPI_RK), DIMENSION(:,:), ALLOCATABLE :: VarN_s1s1
     ! ...

     ! ... Physical Derivatives
     REAL(KIND=SPI_RK), DIMENSION(:,:), ALLOCATABLE :: VarN_x1
     REAL(KIND=SPI_RK), DIMENSION(:,:), ALLOCATABLE :: VarN_x1x1
     ! ...

     REAL(KIND=SPI_RK), DIMENSION(:,:,:), ALLOCATABLE :: sol_analytical

     CLASS(DEF_QUADRATURE_1D), POINTER :: ptr_quad => NULL()
  END TYPE DEF_GREENBOX_1D
  ! ..........................................................        


END MODULE SPI_GREENBOX_DEF
