!# -*- coding: utf8 -*-
MODULE SPI_SPACE_DEF
  USE SPI_MESH_DEF
  USE SPI_NUMBERING_DEF
  USE SPI_QUADRATURES_DEF 
  USE SPI_BASIS_DEF
  USE SPI_BLACKBOX_DEF
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 2 

  ! ..................................................
  TYPE, PUBLIC, ABSTRACT :: DEF_SPACE_ABSTRACT
     CHARACTER(LEN = 1024) :: dirname
  END TYPE DEF_SPACE_ABSTRACT
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_SPACE_ABSTRACT) :: DEF_SPACE_1D
     TYPE(DEF_QUADRATURE_1D) :: oo_quad
  END TYPE DEF_SPACE_1D
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_SPACE_1D) :: DEF_SPACE_1D_BSPLINE
     TYPE(DEF_NUMBERING_1D_BSPLINE) :: oo_numbering
     TYPE(DEF_BASIS_1D_BSPLINE)     :: oo_basis
     TYPE(DEF_BLACKBOX_1D_BSPLINE)  :: oo_bbox

     CLASS(DEF_MESH_1D_BSPLINE)     , POINTER  :: ptr_mesh => NULL()
  END TYPE DEF_SPACE_1D_BSPLINE
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_SPACE_1D) :: DEF_SPACE_1D_FOURIER
     TYPE(DEF_NUMBERING_1D_FOURIER) :: oo_numbering
     TYPE(DEF_BASIS_1D_FOURIER)     :: oo_basis
     TYPE(DEF_BLACKBOX_1D_FOURIER)  :: oo_bbox

     CLASS(DEF_MESH_1D_FOURIER)     , POINTER  :: ptr_mesh => NULL()
  END TYPE DEF_SPACE_1D_FOURIER
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_SPACE_1D) :: DEF_SPACE_1D_HBEZIER
     TYPE(DEF_NUMBERING_1D_HBEZIER) :: oo_numbering
     TYPE(DEF_BASIS_1D_HBEZIER)     :: oo_basis
     TYPE(DEF_BLACKBOX_1D_HBEZIER)  :: oo_bbox

     CLASS(DEF_MESH_1D_HBEZIER)     , POINTER  :: ptr_mesh => NULL()
  END TYPE DEF_SPACE_1D_HBEZIER
  ! ..................................................

END MODULE SPI_SPACE_DEF
