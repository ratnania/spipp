!# -*- coding: utf8 -*-
MODULE SPACE_DEF
  USE MESH_DEF
  USE BLACKBOX_DEF
  USE GREENBOX_DEF
  USE QUADRATURES_DEF 
  USE BASIS_DEF
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 2 

  ! ..................................................
  TYPE, PUBLIC, ABSTRACT :: DEF_SPACE_ABSTRACT
     CHARACTER(LEN = 1024) :: dirname
     CHARACTER(LEN = 1024) :: dirname_other
  END TYPE DEF_SPACE_ABSTRACT
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_SPACE_ABSTRACT) :: DEF_SPACE_2D
     CLASS(DEF_MESH_2D)    , POINTER  :: ptr_mesh => NULL()
     TYPE(DEF_BASIS_2D)     :: basis
  END TYPE DEF_SPACE_2D
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_SPACE_2D) :: DEF_SPACE_QUAD_2D
     TYPE(DEF_QUADRATURE_SQUARE)   :: quadrature
  END TYPE DEF_SPACE_QUAD_2D
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_SPACE_2D) :: DEF_SPACE_TRIANGLE_2D
     TYPE(DEF_QUADRATURE_TRIANGLE) :: quadrature
  END TYPE DEF_SPACE_TRIANGLE_2D
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_SPACE_ABSTRACT) :: DEF_SPACE_2D1D
     CLASS(DEF_MESH_2D)    , POINTER  :: ptr_mesh_2d => NULL()
     CLASS(DEF_MESH_1D)    , POINTER  :: ptr_mesh_1d => NULL()
     TYPE(DEF_BASIS_2D)     :: basis_2d
     TYPE(DEF_BASIS_1D)     :: basis_1d

     TYPE(DEF_QUADRATURE_1D) :: quadrature_1d
  END TYPE DEF_SPACE_2D1D
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_SPACE_2D1D) :: DEF_SPACE_QUAD_2D1D
     TYPE(DEF_QUADRATURE_SQUARE)   :: quadrature_2d
  END TYPE DEF_SPACE_QUAD_2D1D
  ! ..................................................

END MODULE SPACE_DEF
