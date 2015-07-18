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
MODULE FEBasis
  USE TypeDef
  USE BASIS_DEF
  USE MESH_DEF
  USE JOREK_PARAM_DEF
  USE JOREK_PARAM
  USE FEBasis1D_Bezier
  USE FEBasis1D_Fourier
  USE FEBasis2D_Bezier
!  USE FEBasis2D_Splines
  USE FEBasis2D_BoxSplines
  USE FEBasis2D_Bezier_description

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

     self % oi_type = ai_type

     SELECT TYPE (self)
        ! ... 1D CASE
        CLASS IS (DEF_BASIS_1D)
           SELECT TYPE (ao_mesh)
              CLASS IS (DEF_MESH_1D)
                 ! ... toroidal basis
                 SELECT CASE(self % oi_type)
                 ! ... Bezier
                 CASE(INT_TOROIDAL_BASIS_HBEZIER)
                    CALL CREATE_BASIS_1D_HBEZIER(self, ao_mesh)
                 ! ... 
                 ! ... Fourier 
                 CASE(INT_TOROIDAL_BASIS_FOURIER)
                    CALL CREATE_BASIS_1D_FOURIER(self, ao_mesh)
                 ! ... 
                 END SELECT
                 ! ...

              CLASS DEFAULT
                 STOP 'CREATE_BASIS: unexpected type for ao_mesh object!'
           END SELECT
        ! ...

        ! ... 2D CASE
        CLASS IS (DEF_BASIS_2D)
           SELECT TYPE (ao_mesh)
              CLASS IS (DEF_MESH_2D)
                 ! ... poloidal basis
                 SELECT CASE(self % oi_type)
                 ! ... Bezier
                 CASE(INT_POLOIDAL_BASIS_HBEZIER)
                    CALL CREATE_BASIS_2D_HBEZIER(self, ao_mesh)
                 ! ... Box-Splines 
                 CASE(INT_POLOIDAL_BASIS_BOXSPLINES)
                    CALL CREATE_BASIS_2D_BOXSPLINES(self, ao_mesh)
                 ! ... 
                 ! ... generic bezier description 
                 CASE(INT_POLOIDAL_BASIS_BEZIER_DESCRIPTION)
                    CALL CREATE_BASIS_2D_BEZIER_DESCRIPTION(self, ao_mesh)
                 ! ... 
                 END SELECT
              CLASS DEFAULT
                 STOP 'CREATE_BASIS: unexpected type for ao_mesh object!'
           END SELECT
        ! ...

        CLASS DEFAULT
           STOP 'CREATE_BASIS: unexpected type for self object!'
     END SELECT

   END SUBROUTINE CREATE_BASIS
   ! ...................................................

   ! ...................................................
   SUBROUTINE RESET_BASIS(self, ao_mesh, ai_elmt_id)
   IMPLICIT NONE
     CLASS(DEF_BASIS_ABSTRACT), INTENT(INOUT) :: self
     CLASS(DEF_MESH_ABSTRACT), INTENT(INOUT) :: ao_mesh
     INTEGER, INTENT(IN)       :: ai_elmt_id

     SELECT TYPE (self)
        ! ... 1D CASE
!        CLASS IS (DEF_BASIS_1D)
!           SELECT TYPE (ao_mesh)
!              CLASS IS (DEF_MESH_1D)
!                 ! ... toroidal basis
!                 SELECT CASE(self % oi_type)
!                 ! ... Bezier
!                 CASE(INT_TOROIDAL_BASIS_HBEZIER)
!                    CALL RESET_BASIS_1D_HBEZIER(self, ao_mesh, ai_elmt_id)
!                 ! ... 
!                 ! ... Fourier 
!                 CASE(INT_TOROIDAL_BASIS_FOURIER)
!                    CALL RESET_BASIS_1D_FOURIER(self, ao_mesh, ai_elmt_id)
!                 ! ... 
!                 END SELECT
!                 ! ...
!
!              CLASS DEFAULT
!                 STOP 'RESET_BASIS: unexpected type for ao_mesh object!'
!           END SELECT
        ! ...

        ! ... 2D CASE
        CLASS IS (DEF_BASIS_2D)
           SELECT TYPE (ao_mesh)
              CLASS IS (DEF_MESH_2D)
                 ! ... poloidal basis
                 SELECT CASE(self % oi_type)
                 ! ... Bezier
                 CASE(INT_POLOIDAL_BASIS_HBEZIER)
                    CALL RESET_BASIS_2D_HBEZIER(self, ao_mesh, ai_elmt_id)
                 ! ... 
                 ! ... Box-Splines 
                 CASE(INT_POLOIDAL_BASIS_BOXSPLINES)
                    CALL RESET_BASIS_2D_BOXSPLINES(self, ao_mesh, ai_elmt_id)
                 ! ... 
                 ! ... generic bezier description 
                 CASE(INT_POLOIDAL_BASIS_BEZIER_DESCRIPTION)
                    CALL RESET_BASIS_2D_BEZIER_DESCRIPTION(self, ao_mesh, ai_elmt_id)
                 ! ... 
                 END SELECT
              CLASS DEFAULT
                 STOP 'RESET_BASIS: unexpected type for ao_mesh object!'
           END SELECT
        ! ...

        CLASS DEFAULT
           STOP 'RESET_BASIS: unexpected type for self object!'
     END SELECT

   END SUBROUTINE RESET_BASIS
   ! ...................................................

   ! ...................................................
   SUBROUTINE UPDATE_BASIS(self, ao_mesh, ai_elmt_id)
   IMPLICIT NONE
     CLASS(DEF_BASIS_ABSTRACT), INTENT(INOUT) :: self
     CLASS(DEF_MESH_ABSTRACT), INTENT(INOUT) :: ao_mesh
     INTEGER, INTENT(IN)       :: ai_elmt_id

     SELECT TYPE (self)
        ! ... 1D CASE
        CLASS IS (DEF_BASIS_1D)
           SELECT TYPE (ao_mesh)
              CLASS IS (DEF_MESH_1D)
                 ! ... toroidal basis
                 SELECT CASE(self % oi_type)
                 ! ... Bezier
                 CASE(INT_TOROIDAL_BASIS_HBEZIER)
                    CALL UPDATE_BASIS_1D_HBEZIER(self, ao_mesh, ai_elmt_id)
                 ! ... 
                 ! ... Fourier 
                 CASE(INT_TOROIDAL_BASIS_FOURIER)
                    CALL UPDATE_BASIS_1D_FOURIER(self, ao_mesh, ai_elmt_id)
                 ! ... 
                 END SELECT
                 ! ...

              CLASS DEFAULT
                 STOP 'UPDATE_BASIS: unexpected type for ao_mesh object!'
           END SELECT
        ! ...

        ! ... 2D CASE
        CLASS IS (DEF_BASIS_2D)
           SELECT TYPE (ao_mesh)
              CLASS IS (DEF_MESH_2D)
                 ! ... poloidal basis
                 SELECT CASE(self % oi_type)
                 ! ... Bezier
                 CASE(INT_POLOIDAL_BASIS_HBEZIER)
                    CALL UPDATE_BASIS_2D_HBEZIER(self, ao_mesh, ai_elmt_id)
                 ! ... 
                 ! ... Box-Splines 
                 CASE(INT_POLOIDAL_BASIS_BOXSPLINES)
                    CALL UPDATE_BASIS_2D_BOXSPLINES(self, ao_mesh, ai_elmt_id)
                 ! ... 
                 ! ... generic bezier description 
                 CASE(INT_POLOIDAL_BASIS_BEZIER_DESCRIPTION)
                    CALL UPDATE_BASIS_2D_BEZIER_DESCRIPTION(self, ao_mesh, ai_elmt_id)
                 ! ... 
                 END SELECT
              CLASS DEFAULT
                 STOP 'UPDATE_BASIS: unexpected type for ao_mesh object!'
           END SELECT
        ! ...

        CLASS DEFAULT
           STOP 'UPDATE_BASIS: unexpected type for self object!'
     END SELECT

   END SUBROUTINE UPDATE_BASIS
   ! ...................................................
END MODULE FEBasis
