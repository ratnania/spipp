!# -*- coding: utf8 -*-
MODULE SPI_BASIS
  USE SPI_BASIS_DEF
  USE SPI_BASIS_BSPLINE
  USE SPI_BASIS_FOURIER 
  USE SPI_BASIS_HBEZIER
  USE SPI_MESH_DEF
  USE SPI_GLOBAL_DEF
  IMPLICIT NONE

CONTAINS
   ! ...................................................
   SUBROUTINE CREATE_BASIS(self, ao_mesh, ao_quad, dirname)
   IMPLICIT NONE
     CLASS(DEF_BASIS_ABSTRACT)         , INTENT(INOUT) :: self
     CLASS(DEF_MESH_ABSTRACT)          , INTENT(INOUT) :: ao_mesh
     CLASS(DEF_QUADRATURE_1D), TARGET  , INTENT(INOUT) :: ao_quad
     CHARACTER(LEN = 1024)   , OPTIONAL, INTENT(IN)    :: dirname
     ! LOCAL

     IF (PRESENT(dirname)) THEN
        self % dirname = TRIM(dirname)
     END IF

     ! ...
     SELECT TYPE (self)

     ! ...
     CLASS IS (DEF_BASIS_1D_BSPLINE)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D_BSPLINE)
           self % ptr_quad => ao_quad

           CALL CREATE_BASIS_1D_BSPLINE(self, ao_mesh)
        CLASS DEFAULT
           STOP 'CREATE_BASIS: unexpected type for ao_mesh object!'
        END SELECT
     ! ...

     ! ...
     CLASS IS (DEF_BASIS_1D_FOURIER)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D_FOURIER)
           self % ptr_quad => ao_quad

           CALL CREATE_BASIS_1D_FOURIER(self, ao_mesh)
        CLASS DEFAULT
           STOP 'CREATE_BASIS: unexpected type for ao_mesh object!'
        END SELECT
     ! ...

     ! ...
     CLASS IS (DEF_BASIS_1D_HBEZIER)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D_HBEZIER)
           self % ptr_quad => ao_quad

           CALL CREATE_BASIS_1D_HBEZIER(self, ao_mesh)
        CLASS DEFAULT
           STOP 'CREATE_BASIS: unexpected type for ao_mesh object!'
        END SELECT
     ! ...

     CLASS DEFAULT
        STOP 'CREATE_BASIS: unexpected type for self object!'
     END SELECT
     ! ...

   END SUBROUTINE CREATE_BASIS
   ! ...................................................

   ! ...................................................
   SUBROUTINE FREE_BASIS(self)
   IMPLICIT NONE
     CLASS(DEF_BASIS_ABSTRACT), INTENT(INOUT) :: self
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_BASIS_1D_BSPLINE)
        CALL FREE_BASIS_1D_BSPLINE(self)
     CLASS IS (DEF_BASIS_1D_HBEZIER)
        CALL FREE_BASIS_1D_HBEZIER(self)
     CLASS IS (DEF_BASIS_1D_FOURIER)
        CALL FREE_BASIS_1D_FOURIER(self)
     CLASS DEFAULT
        STOP 'FREE_BASIS: unexpected type for self object!'
     END SELECT
     ! ...
              
   END SUBROUTINE FREE_BASIS
   ! ...................................................

   ! ...................................................
   SUBROUTINE RESET_BASIS(self)
   IMPLICIT NONE
     CLASS(DEF_BASIS_ABSTRACT), INTENT(INOUT) :: self
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_BASIS_1D_BSPLINE)
        CALL RESET_BASIS_1D_BSPLINE(self)
     CLASS IS (DEF_BASIS_1D_HBEZIER)
        CALL RESET_BASIS_1D_HBEZIER(self)
     CLASS IS (DEF_BASIS_1D_FOURIER)
        CALL RESET_BASIS_1D_FOURIER(self)
     CLASS DEFAULT
        STOP 'RESET_BASIS: unexpected type for self object!'
     END SELECT
     ! ...
              
   END SUBROUTINE RESET_BASIS
   ! ...................................................

   ! ...................................................
   SUBROUTINE UPDATE_BASIS(self, apr_points)
   IMPLICIT NONE
     CLASS(DEF_BASIS_ABSTRACT), INTENT(INOUT) :: self
     REAL(SPI_RK), DIMENSION(:) :: apr_points
     ! LOCAL

     ! ...
     SELECT TYPE (self)

     ! ...
     CLASS IS (DEF_BASIS_1D_BSPLINE)
        CALL UPDATE_BASIS_1D_BSPLINE(self, apr_points)
     ! ...

     ! ...
     CLASS IS (DEF_BASIS_1D_FOURIER)
        CALL UPDATE_BASIS_1D_FOURIER(self, apr_points)
     ! ...

     ! ...
     CLASS IS (DEF_BASIS_1D_HBEZIER)
        CALL UPDATE_BASIS_1D_HBEZIER(self, apr_points)
     ! ...

     CLASS DEFAULT
        STOP 'UPDATE_BASIS: unexpected type for self object!'
     END SELECT
     ! ...

   END SUBROUTINE UPDATE_BASIS
   ! ...................................................

END MODULE SPI_BASIS
