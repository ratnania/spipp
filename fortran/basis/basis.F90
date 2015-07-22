!# -*- coding: utf8 -*-
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
!*                                                                           */
!*                  This file is part of the PlaTo program                   */
!*                                                                           */
!* Copyright (C) 2011-2012 B. Nkonga                                         */
!*                                                                           */
!*  PlaTo  is distributed under the terms of the Cecill-B License.           */
!*                                                                           */
!*  You can find a copy of the Cecill-B License at :                         */
!*  http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt                 */
!*                                                                           */
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
MODULE SPI_BASIS
  USE SPI_BASIS_DEF
  USE SPI_BASIS_BSPLINES 
  USE SPI_BASIS_FOURIER 
  USE SPI_MESH_DEF
  USE SPI_GLOBAL_DEF
  IMPLICIT NONE

CONTAINS
   ! ...................................................
   SUBROUTINE CREATE_BASIS(self, ao_mesh, ai_type, dirname)
   IMPLICIT NONE
     CLASS(DEF_BASIS_ABSTRACT)      , INTENT(INOUT) :: self
     CLASS(DEF_MESH_ABSTRACT)       , INTENT(INOUT) :: ao_mesh
     INTEGER                        , INTENT(IN)    :: ai_type
     CHARACTER(LEN = 1024), OPTIONAL, INTENT(IN)    :: dirname
     ! LOCAL

     IF (PRESENT(dirname)) THEN
        self % dirname = TRIM(dirname)
     END IF

     ! ...
     SELECT TYPE (self)

     ! ...
     CLASS IS (DEF_BASIS_1D_BSPLINES)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D)
           CALL CREATE_BASIS_1D_BSPLINES(self, ao_mesh)
        CLASS DEFAULT
           STOP 'CREATE_BASIS: unexpected type for ao_mesh object!'
        END SELECT
     ! ...

     ! ...
     CLASS IS (DEF_BASIS_1D_FOURIER)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D)
           CALL CREATE_BASIS_1D_FOURIER(self, ao_mesh)
        CLASS DEFAULT
           STOP 'CREATE_BASIS: unexpected type for ao_mesh object!'
        END SELECT
     ! ...

     ! ...
     CLASS IS (DEF_BASIS_1D_HBEZIER)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D)
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
   SUBROUTINE RESET_BASIS(self, ao_mesh, ai_elmt_id)
   IMPLICIT NONE
     CLASS(DEF_BASIS_ABSTRACT), INTENT(INOUT) :: self
     CLASS(DEF_MESH_ABSTRACT), INTENT(INOUT) :: ao_mesh
     INTEGER, INTENT(IN)       :: ai_elmt_id
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_BASIS_1D_BSPLINES)
        CALL RESET_BASIS_1D_BSPLINES(self)
     CLASS IS (DEF_BASIS_1D_HBEZIER)
        CALL RESET_BASIS_1D_FOURIER(self)
     CLASS IS (DEF_BASIS_1D_FOURIER)
        CALL RESET_BASIS_1D_HBEZIER(self)
     CLASS DEFAULT
        STOP 'RESET_BASIS: unexpected type for self object!'
     END SELECT
     ! ...
              
   END SUBROUTINE RESET_BASIS
   ! ...................................................

   ! ...................................................
   SUBROUTINE UPDATE_BASIS(self, ao_mesh, ai_elmt_id)
   IMPLICIT NONE
     CLASS(DEF_BASIS_ABSTRACT), INTENT(INOUT) :: self
     CLASS(DEF_MESH_ABSTRACT), INTENT(INOUT) :: ao_mesh
     INTEGER, INTENT(IN)       :: ai_elmt_id
     ! LOCAL

     ! ...
     SELECT TYPE (self)

     ! ...
     CLASS IS (DEF_BASIS_1D_BSPLINES)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D)
           CALL UPDATE_BASIS_1D_BSPLINES(self, ao_mesh, ai_elmt_id)
        CLASS DEFAULT
           STOP 'UPDATE_BASIS: unexpected type for ao_mesh object!'
        END SELECT
     ! ...

     ! ...
     CLASS IS (DEF_BASIS_1D_FOURIER)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D)
           CALL UPDATE_BASIS_1D_FOURIER(self, ao_mesh, ai_elmt_id)
        CLASS DEFAULT
           STOP 'UPDATE_BASIS: unexpected type for ao_mesh object!'
        END SELECT
     ! ...

     ! ...
     CLASS IS (DEF_BASIS_1D_HBEZIER)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D)
           CALL UPDATE_BASIS_1D_HBEZIER(self, ao_mesh, ai_elmt_id)
        CLASS DEFAULT
           STOP 'UPDATE_BASIS: unexpected type for ao_mesh object!'
        END SELECT
     ! ...

     CLASS DEFAULT
        STOP 'UPDATE_BASIS: unexpected type for self object!'
     END SELECT
     ! ...

   END SUBROUTINE UPDATE_BASIS
   ! ...................................................
END MODULE SPI_BASIS
