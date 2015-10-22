!# -*- coding: utf8 -*-
MODULE SPI_BLACKBOX_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  USE SPI_BASIS_DEF
  IMPLICIT NONE

  ! ..........................................................        
  TYPE, PUBLIC, ABSTRACT :: DEF_BLACKBOX_ABSTRACT
     INTEGER :: oi_n_points
     INTEGER :: oi_n_elements
     INTEGER :: n_dim
     INTEGER :: oi_basis_type

     ! ...
     REAL(KIND=SPI_RK), DIMENSION(:, :, :), ALLOCATABLE  :: Xp_0
     ! ...    

     ! ... Position Logical Derivatives
     REAL(KIND=SPI_RK), DIMENSION(:, :, :), ALLOCATABLE  :: Xp_s1 
     REAL(KIND=SPI_RK), DIMENSION(:, :, :), ALLOCATABLE  :: Xp_s1s1 
     ! ...

     ! ... Position Physical Derivatives
     REAL(KIND=SPI_RK), DIMENSION(:, :, :), ALLOCATABLE  :: Xp_x1 
     REAL(KIND=SPI_RK), DIMENSION(:, :, :), ALLOCATABLE  :: Xp_x1x1 
     ! ...

     ! ... Volume and wieghted volume
     REAL(KIND=SPI_RK), DIMENSION(:, :), ALLOCATABLE  :: Jacobians 
     REAL(KIND=SPI_RK), DIMENSION(:, :), ALLOCATABLE  :: Vol
     REAL(KIND=SPI_RK), DIMENSION(:, :), ALLOCATABLE  :: wVol
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
     REAL(SPI_RK), DIMENSION (:,:), ALLOCATABLE :: control_points

     CLASS(DEF_BASIS_1D_BSPLINE), POINTER :: ptr_basis => NULL()
     CLASS(DEF_MESH_1D_BSPLINE), POINTER :: ptr_mesh => NULL()
  END TYPE DEF_BLACKBOX_1D_BSPLINE
  ! ..........................................................        

  ! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BLACKBOX_ABSTRACT) :: DEF_BLACKBOX_2D_BSPLINE
     REAL(SPI_RK), DIMENSION (:, :, :), ALLOCATABLE :: control_points

     CLASS(DEF_BASIS_1D_BSPLINE), POINTER :: ptr_basis_u => NULL()
     CLASS(DEF_MESH_1D_BSPLINE) , POINTER :: ptr_mesh_u  => NULL()

     CLASS(DEF_BASIS_1D_BSPLINE), POINTER :: ptr_basis_v => NULL()
     CLASS(DEF_MESH_1D_BSPLINE) , POINTER :: ptr_mesh_v  => NULL()
  END TYPE DEF_BLACKBOX_2D_BSPLINE
  ! ..........................................................        

  ! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BLACKBOX_ABSTRACT) :: DEF_BLACKBOX_3D_BSPLINE
     REAL(SPI_RK), DIMENSION (:, :, :, :), ALLOCATABLE :: control_points

     CLASS(DEF_BASIS_1D_BSPLINE), POINTER :: ptr_basis_u => NULL()
     CLASS(DEF_MESH_1D_BSPLINE) , POINTER :: ptr_mesh_u  => NULL()

     CLASS(DEF_BASIS_1D_BSPLINE), POINTER :: ptr_basis_v => NULL()
     CLASS(DEF_MESH_1D_BSPLINE) , POINTER :: ptr_mesh_v  => NULL()

     CLASS(DEF_BASIS_1D_BSPLINE), POINTER :: ptr_basis_w => NULL()
     CLASS(DEF_MESH_1D_BSPLINE) , POINTER :: ptr_mesh_w  => NULL()
  END TYPE DEF_BLACKBOX_3D_BSPLINE
  ! ..........................................................        

  ! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BLACKBOX_1D) :: DEF_BLACKBOX_1D_FOURIER
     CLASS(DEF_BASIS_1D_FOURIER), POINTER :: ptr_basis => NULL()
     CLASS(DEF_MESH_1D_FOURIER), POINTER :: ptr_mesh => NULL()
  END TYPE DEF_BLACKBOX_1D_FOURIER
  ! ..........................................................        

  ! ..........................................................        
  TYPE, PUBLIC, EXTENDS(DEF_BLACKBOX_1D) :: DEF_BLACKBOX_1D_HBEZIER
     CLASS(DEF_BASIS_1D_HBEZIER), POINTER :: ptr_basis => NULL()
     CLASS(DEF_MESH_1D_HBEZIER), POINTER :: ptr_mesh => NULL()
  END TYPE DEF_BLACKBOX_1D_HBEZIER
  ! ..........................................................        

END MODULE SPI_BLACKBOX_DEF
