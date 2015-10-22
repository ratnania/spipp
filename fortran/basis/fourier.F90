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
  PUBLIC :: CREATE_BASIS_1D_FOURIER, &
          & RESET_BASIS_1D_FOURIER,  &
          & UPDATE_BASIS_1D_FOURIER, &
          & FREE_BASIS_1D_FOURIER

  INTEGER :: n_tor
  INTEGER :: n_plan
  INTEGER :: n_period
  INTEGER :: ierror

CONTAINS

   ! ...................................................
   SUBROUTINE CREATE_BASIS_1D_FOURIER(self, ao_mesh)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D_FOURIER)       , INTENT(INOUT) :: self
     CLASS(DEF_MESH_1D_FOURIER), TARGET, INTENT(INOUT) :: ao_mesh

     self % ptr_mesh => ao_mesh

     CALL Init_FE_BasisF1D_Parameters(self)
     CALL Init_FE_BasisF1D(self)

   END SUBROUTINE CREATE_BASIS_1D_FOURIER
   ! ...................................................

   ! ...................................................
   SUBROUTINE FREE_BASIS_1D_FOURIER(self)
    IMPLICIT NONE
     CLASS(DEF_BASIS_1D_FOURIER), INTENT(INOUT) :: self

   END SUBROUTINE FREE_BASIS_1D_FOURIER
   ! ...................................................

   ! ...................................................
   SUBROUTINE RESET_BASIS_1D_FOURIER(self)
    IMPLICIT NONE
     CLASS(DEF_BASIS_1D_FOURIER), INTENT(INOUT) :: self

   END SUBROUTINE RESET_BASIS_1D_FOURIER
   ! ...................................................

   ! ...................................................
   SUBROUTINE UPDATE_BASIS_1D_FOURIER(self, apr_points)
    IMPLICIT NONE
     CLASS(DEF_BASIS_1D_FOURIER), INTENT(INOUT) :: self
     REAL(SPI_RK), DIMENSION(:) :: apr_points

   END SUBROUTINE UPDATE_BASIS_1D_FOURIER
   ! ...................................................

   ! ...................................................
   SUBROUTINE Init_FE_BasisF1D_Parameters(self)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D_FOURIER), INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: li_err

     self % ptr_mesh % oi_n_vtex_per_elmt     = 1

     ! ... TODO to set
     n_tor    = SPI_INT_DEFAULT 
     n_plan   = SPI_INT_DEFAULT
     n_period = SPI_INT_DEFAULT
     ! ...

     self % oi_n_order     = n_tor
     self % oi_n_max_order = n_tor

    END SUBROUTINE Init_FE_BasisF1D_Parameters
   ! ...................................................

   ! ...................................................
   SUBROUTINE Init_FE_BasisF1D(self)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D_FOURIER), INTENT(INOUT) :: self
    ! LOCAL
    INTEGER       :: iv, jv, ig, il
    REAL(KIND=SPI_RK) :: s , t, phi
    INTEGER :: mode(n_tor) 

     ALLOCATE(self % TestfT_0 (self % ptr_quad % oi_n_points, self % oi_n_order, self % ptr_mesh % oi_n_vtex_per_elmt))
     ALLOCATE(self % TestfT_p (self % ptr_quad % oi_n_points, self % oi_n_order, self % ptr_mesh % oi_n_vtex_per_elmt))
     ALLOCATE(self % TestfT_pp(self % ptr_quad % oi_n_points, self % oi_n_order, self % ptr_mesh % oi_n_vtex_per_elmt))

    mode(1)=0
    DO il=2,n_tor
       mode(il)=int(il/2)/n_period
    END DO
    

    DO ig = 1, self % ptr_quad % oi_n_points


      CALL Compute_fourier_mods(n_tor,n_period,1,phi,self % TestfT_0(ig,1,1))
      CALL Compute_FD_fourier_mods(n_tor,n_period,1,phi,self % TestfT_p(ig,1,1)) 

        DO il=1,(n_tor-1)/2
  
        phi= self % ptr_quad % opr_points(1,ig)
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
