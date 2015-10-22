!# -*- coding: utf8 -*-
MODULE SPI_MESH_TRANSFORM 
  USE SPI_MESH_DEF
  USE SPI_GLOBAL_DEF
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: MESH_TRANSLATE_1D

CONTAINS

  ! .........................................................
  SUBROUTINE MESH_TRANSLATE_1D(self, apr_displ)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE), INTENT(INOUT)  :: self
     REAL(KIND=SPI_RK), DIMENSION(:)       :: apr_displ
     ! LOCAL
     INTEGER :: li_mesh
     INTEGER :: li_err 


  END SUBROUTINE MESH_TRANSLATE_1D
  ! .........................................................

END MODULE SPI_MESH_TRANSFORM
