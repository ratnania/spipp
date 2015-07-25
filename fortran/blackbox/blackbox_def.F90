!# -*- coding: utf8 -*-
MODULE SPI_BLACKBOX_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  USE SPI_BASIS_DEF
  IMPLICIT NONE

  ! ..........................................................        
  TYPE, PUBLIC, ABSTRACT :: DEF_BLACKBOX_ABSTRACT
     INTEGER :: oi_n_points
     INTEGER :: oi_basis_type

     REAL(KIND=SPI_RK), DIMENSION(:), ALLOCATABLE  :: Vol
     REAL(KIND=SPI_RK), DIMENSION(:), ALLOCATABLE  :: wVol

     ! ...
     REAL(KIND=SPI_RK), DIMENSION(:, :), ALLOCATABLE  :: Xp_0
     ! ...    

     ! ... Position Logical Derivatives
     REAL(KIND=SPI_RK), DIMENSION(:, :), ALLOCATABLE  :: Xp_s1 
     REAL(KIND=SPI_RK), DIMENSION(:, :), ALLOCATABLE  :: Xp_s1s1 
     ! ...

     ! ... Position Physical Derivatives
     REAL(KIND=SPI_RK), DIMENSION(:, :), ALLOCATABLE  :: Xp_x1 
     REAL(KIND=SPI_RK), DIMENSION(:, :), ALLOCATABLE  :: Xp_x1x1 
     ! ...

     ! ...
     REAL(KIND=SPI_RK), DIMENSION(:), ALLOCATABLE :: B_0
     ! ...

     ! ... Basis Logical Derivatives
     REAL(KIND=SPI_RK), DIMENSION(:), ALLOCATABLE :: B_s1
     REAL(KIND=SPI_RK), DIMENSION(:), ALLOCATABLE :: B_s1s1
     ! ...

     ! ... Basis Physical Derivatives
     REAL(KIND=SPI_RK), DIMENSION(:), ALLOCATABLE :: B_x1
     REAL(KIND=SPI_RK), DIMENSION(:), ALLOCATABLE :: B_x1x1
     ! ...
    
     CLASS(DEF_QUADRATURE_1D), POINTER :: ptr_quad => NULL()
  END TYPE DEF_BLACKBOX_ABSTRACT
  ! ..........................................................        

  ! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BLACKBOX_ABSTRACT) :: DEF_BLACKBOX_1D
  END TYPE DEF_BLACKBOX_1D
  ! ..........................................................        

  ! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BLACKBOX_1D) :: DEF_BLACKBOX_1D_BSPLINE
     CLASS(DEF_BASIS_1D_BSPLINE), POINTER :: ptr_basis => NULL()
  END TYPE DEF_BLACKBOX_1D_BSPLINE
  ! ..........................................................        

  ! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BLACKBOX_1D) :: DEF_BLACKBOX_1D_FOURIER
     CLASS(DEF_BASIS_1D_FOURIER), POINTER :: ptr_basis => NULL()
  END TYPE DEF_BLACKBOX_1D_FOURIER
  ! ..........................................................        

  ! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BLACKBOX_1D) :: DEF_BLACKBOX_1D_HBEZIER
     CLASS(DEF_BASIS_1D_HBEZIER), POINTER :: ptr_basis => NULL()
  END TYPE DEF_BLACKBOX_1D_HBEZIER
  ! ..........................................................        

END MODULE SPI_BLACKBOX_DEF
