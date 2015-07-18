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
MODULE SPI_BASIS_FOURIER 
  USE SPI_BASIS_DEF
  USE SPI_FOURIER_MODS
  USE SPI_MESH_DEF
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CREATE_BASIS_1D_FOURIER, RESET_BASIS_1D_FOURIER, UPDATE_BASIS_1D_FOURIER

  INTEGER :: n_tor
  INTEGER :: n_plane
  INTEGER :: n_period
  INTEGER :: ierror

CONTAINS

   ! ...................................................
   SUBROUTINE CREATE_BASIS_1D_FOURIER(self, ao_mesh)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D), INTENT(INOUT) :: self
     CLASS(DEF_MESH_1D), INTENT(INOUT) :: ao_mesh

     CALL Init_FE_BasisF1D_Parameters(self, ao_mesh)
     CALL Init_FE_BasisF1D(self, ao_mesh)

   END SUBROUTINE CREATE_BASIS_1D_FOURIER
   ! ...................................................

   ! ...................................................
   SUBROUTINE RESET_BASIS_1D_FOURIER(self, ao_mesh, ai_elmt_id)
    IMPLICIT NONE
     CLASS(DEF_BASIS_1D), INTENT(INOUT) :: self
     CLASS(DEF_MESH_1D), INTENT(INOUT) :: ao_mesh
     INTEGER, INTENT(IN)       :: ai_elmt_id

   END SUBROUTINE RESET_BASIS_1D_FOURIER
   ! ...................................................

   ! ...................................................
   SUBROUTINE UPDATE_BASIS_1D_FOURIER(self, ao_mesh, ai_elmt_id)
    IMPLICIT NONE
     CLASS(DEF_BASIS_1D), INTENT(INOUT) :: self
     CLASS(DEF_MESH_1D), INTENT(INOUT) :: ao_mesh
     INTEGER, INTENT(IN)       :: ai_elmt_id

   END SUBROUTINE UPDATE_BASIS_1D_FOURIER
   ! ...................................................

   ! ...................................................
   SUBROUTINE Init_FE_BasisF1D_Parameters(self, ao_mesh)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D), INTENT(INOUT) :: self
     CLASS(DEF_MESH_1D), INTENT(INOUT) :: ao_mesh
     ! LOCAL
     INTEGER :: li_err

     ao_mesh % oi_n_vtex_per_elmt     = 1
     ao_mesh % oi_n_max_vtex_per_elmt = 1

     CALL JOREK_Param_GETInt(INT_N_TOR_ID,n_tor,li_err)
     CALL JOREK_Param_GETInt(INT_N_PLANE_ID,n_plane,li_err)
     CALL JOREK_Param_GETInt(INT_N_PERIOD_ID,n_period,li_err)

     self % oi_n_order     = n_tor
     self % oi_n_max_order = n_tor

     ao_mesh % oi_n_order             = self % oi_n_order

    END SUBROUTINE Init_FE_BasisF1D_Parameters
   ! ...................................................

   ! ...................................................
   SUBROUTINE Init_FE_BasisF1D(self, ao_mesh)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D), INTENT(INOUT) :: self
     CLASS(DEF_MESH_1D), INTENT(INOUT) :: ao_mesh
    ! LOCAL
    INTEGER       :: iv, jv, ig, il
    REAL(KIND=RK) :: s , t, phi
    INTEGER :: mode(n_tor) 

     ALLOCATE(self % TestfT_0 (ao_mesh % ptr_quad % oi_n_points, self % oi_n_order, ao_mesh % oi_n_vtex_per_elmt))
     ALLOCATE(self % TestfT_p (ao_mesh % ptr_quad % oi_n_points, self % oi_n_order, ao_mesh % oi_n_vtex_per_elmt))
     ALLOCATE(self % TestfT_pp(ao_mesh % ptr_quad % oi_n_points, self % oi_n_order, ao_mesh % oi_n_vtex_per_elmt))

    mode(1)=0
    DO il=2,n_tor
       mode(il)=int(il/2)/n_period
    END DO
    

    DO ig = 1, ao_mesh % ptr_quad % oi_n_points


      CALL Compute_fourier_mods(n_tor,n_period,1,phi,self % TestfT_0(ig,1,1))
      CALL Compute_FD_fourier_mods(n_tor,n_period,1,phi,self % TestfT_p(ig,1,1)) 

        DO il=1,(n_tor-1)/2
  
        phi= ao_mesh % ptr_quad % opr_points(1,ig)
        CALL Compute_fourier_mods(n_tor,n_period,2*il,phi,self % TestfT_0(ig,2*il,1))
        CALL Compute_FD_fourier_mods(n_tor,n_period,2*il,phi,self % TestfT_p(ig,2*il,1)) 
        CALL Compute_SD_fourier_mods(n_tor,n_period,2*il,phi,self % TestfT_pp(ig,2*il,1)) 

        CALL Compute_fourier_mods(n_tor,n_period,2*il+1,phi,self % TestfT_0(ig,2*il+1,1))
        CALL Compute_FD_fourier_mods(n_tor,n_period,2*il+1,phi,self % TestfT_p(ig,2*il+1,1)) 
        CALL Compute_SD_fourier_mods(n_tor,n_period,2*il+1,phi,self % TestfT_pp(ig,2*il+1,1)) 
        
        END DO
    END DO

  END SUBROUTINE Init_FE_BasisF1D


END MODULE SPI_BASIS_FOURIER
