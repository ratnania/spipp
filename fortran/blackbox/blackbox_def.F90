!# -*- coding: utf8 -*-
MODULE SPI_BLACKBOX_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
!  USE SPI_BASIS_DEF

! ..........................................................        
  TYPE, PUBLIC :: DEF_BLACKBOX_1D
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
    
     ! TODO
!     CLASS(DEF_BASIS_1D), POINTER :: ptr_basis => NULL()
     CLASS(DEF_QUADRATURE_1D), POINTER :: ptr_quad => NULL()
  END TYPE DEF_BLACKBOX_1D
  ! -----------------------------------------------------------------------------

END MODULE SPI_BLACKBOX_DEF
