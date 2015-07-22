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
MODULE SPI_BASIS_BSPLINES 
  USE SPI_BASIS_DEF
  USE SPI_MESH_DEF

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CREATE_BASIS_1D_BSPLINES, RESET_BASIS_1D_BSPLINES, UPDATE_BASIS_1D_BSPLINES

CONTAINS

   ! ...................................................
   SUBROUTINE CREATE_BASIS_1D_BSPLINES(self, ao_mesh)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D_BSPLINES), INTENT(INOUT) :: self
     CLASS(DEF_MESH_1D), INTENT(INOUT) :: ao_mesh

     CALL Init_FE_BasisF1D_Parameters(self, ao_mesh)
     CALL Init_FE_BasisF1D(self, ao_mesh)

   END SUBROUTINE CREATE_BASIS_1D_BSPLINES
   ! ...................................................

   ! ...................................................
   SUBROUTINE RESET_BASIS_1D_BSPLINES(self)
    IMPLICIT NONE
     CLASS(DEF_BASIS_1D_BSPLINES), INTENT(INOUT) :: self

   END SUBROUTINE RESET_BASIS_1D_BSPLINES
   ! ...................................................

   ! ...................................................
   SUBROUTINE UPDATE_BASIS_1D_BSPLINES(self, ao_mesh, ai_elmt_id)
    IMPLICIT NONE
     CLASS(DEF_BASIS_1D_BSPLINES), INTENT(INOUT) :: self
     CLASS(DEF_MESH_1D), INTENT(INOUT) :: ao_mesh
     INTEGER, INTENT(IN)       :: ai_elmt_id

   END SUBROUTINE UPDATE_BASIS_1D_BSPLINES
   ! ...................................................

   ! ...................................................
   SUBROUTINE INITIALIZE_PARAMETERS_BASIS_1D_BSPLINES(self, ao_mesh)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D_BSPLINES), INTENT(INOUT) :: self
     CLASS(DEF_MESH_1D), INTENT(INOUT) :: ao_mesh


   END SUBROUTINE INITIALIZE_PARAMETERS_BASIS_1D_BSPLINES
   ! ...................................................

   ! ...................................................
   SUBROUTINE Init_FE_BasisF1D(self, ao_mesh)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D_BSPLINES), INTENT(INOUT) :: self
     CLASS(DEF_MESH_1D), INTENT(INOUT) :: ao_mesh
     ! LOCAL
     INTEGER       :: iv, jv, ig, il
     REAL(KIND=SPI_RK) :: s , t, phi
     REAL(KIND=SPI_RK), DIMENSION(ao_mesh % oi_n_vtex_per_elmt,self % oi_n_order)   :: BT   ,  BT_p , BT_pp
    
    
     ALLOCATE(self % TestfT_0 (ao_mesh % ptr_quad % oi_n_points, self % oi_n_order, ao_mesh % oi_n_vtex_per_elmt))
     ALLOCATE(self % TestfT_p (ao_mesh % ptr_quad % oi_n_points, self % oi_n_order, ao_mesh % oi_n_vtex_per_elmt))
     ALLOCATE(self % TestfT_pp(ao_mesh % ptr_quad % oi_n_points, self % oi_n_order, ao_mesh % oi_n_vtex_per_elmt))
    
     DO ig = 1, ao_mesh % ptr_quad % oi_n_points
    
        phi = ao_mesh % ptr_quad % opr_points(1,ig)
    
        CALL BasisFunctions1D(phi, ao_mesh % oi_n_vtex_per_elmt, self % oi_n_order, BT, BT_p, BT_pp)
        DO il = 1, self % oi_n_order
           self % TestfT_0  (ig, il, 1:ao_mesh % oi_n_vtex_per_elmt) = BT(1:ao_mesh % oi_n_vtex_per_elmt,il)
           self % TestfT_p  (ig, il, 1:ao_mesh % oi_n_vtex_per_elmt) = BT_p(1:ao_mesh % oi_n_vtex_per_elmt,il)
           self % TestfT_pp (ig, il, 1:ao_mesh % oi_n_vtex_per_elmt) = BT_pp(1:ao_mesh % oi_n_vtex_per_elmt,il)
        END DO
     END DO

   END SUBROUTINE Init_FE_BasisF1D
   ! ...................................................

   ! ...................................................
  SUBROUTINE basisfunctions1D(s, ai_n_vtex_per_elmt_Tor, ai_n_order_Tor,Bt,Bt_s,Bt_ss)

    ! --- Routine parameters
    REAL(KIND=SPI_RK), INTENT(in)  :: s          !< s-coordinate in the element
    INTEGER                    :: ai_n_vtex_per_elmt_Tor
    INTEGER                    :: ai_n_order_Tor
    REAL(KIND=SPI_RK), INTENT(out) :: Bt(ai_n_vtex_per_elmt_Tor,ai_n_order_Tor)     !< Basis functions
    REAL(KIND=SPI_RK), INTENT(out) :: Bt_s(ai_n_vtex_per_elmt_Tor,ai_n_order_Tor)   !< Basis functions derived with respect to s
    REAL(KIND=SPI_RK), INTENT(out) :: Bt_ss(ai_n_vtex_per_elmt_Tor,ai_n_order_Tor)  !< Basis functions derived two times with respect to s

    !---------------------------------------------------------- vertex (1)
    Bt   (1,1)=  2.0*(s**3) - 3.0*(s**2) + 1
    Bt_s (1,1)=  6.0*(s**2) - 6.0*s
    Bt_ss(1,1)= 12.0*s      - 6.0

    Bt   (1,2)=    (s**3) - 2.0*(s**2) + s
    Bt_s (1,2)=3.0*(s**2) - 4.0*s      + 1
    Bt_ss(1,2)=6.0*s      - 4.0


    !---------------------------------------------------------- vertex (2)
    Bt   (2,1)= -  2.0*(s**3) + 3.0*(s**2) 
    Bt_s (2,1)= -  6.0*(s**2) + 6.0*s
    Bt_ss(2,1)= - 12.0*s      + 6.0

    Bt   (2,2)=    (s**3) -     (s**2)
    Bt_s (2,2)=3.0*(s**2) - 2.0*s
    Bt_ss(2,2)=6.0*s      - 2.0


  END SUBROUTINE basisfunctions1D
   ! ...................................................

END MODULE SPI_BASIS_BSPLINES
