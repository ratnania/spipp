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
MODULE SPI_BASIS_BSPLINE 
  USE SPI_BASIS_DEF
  USE SPI_MESH_DEF
  USE BSP

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CREATE_BASIS_1D_BSPLINE, &
          & RESET_BASIS_1D_BSPLINE,  &
          & UPDATE_BASIS_1D_BSPLINE, &
          & FREE_BASIS_1D_BSPLINE

CONTAINS

   ! ...................................................
   SUBROUTINE CREATE_BASIS_1D_BSPLINE(self, ao_mesh)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D_BSPLINE)       , INTENT(INOUT) :: self
     CLASS(DEF_MESH_1D_BSPLINE), TARGET, INTENT(INOUT) :: ao_mesh
     ! LOCAL
     INTEGER :: n_elements
     INTEGER :: n_points

     self % ptr_mesh => ao_mesh

     self % oi_n_order = self % ptr_mesh % oi_p + 1 

     n_elements = ao_mesh % n_elements 
     n_points   = self % ptr_quad % oi_n_points 

     ALLOCATE(self % TestfT_0 (n_elements, n_points, self % oi_n_order))
     ALLOCATE(self % TestfT_p (n_elements, n_points, self % oi_n_order))
     ALLOCATE(self % TestfT_pp(n_elements, n_points, self % oi_n_order))

   END SUBROUTINE CREATE_BASIS_1D_BSPLINE
   ! ...................................................

   ! ...................................................
   SUBROUTINE FREE_BASIS_1D_BSPLINE(self)
    IMPLICIT NONE
     CLASS(DEF_BASIS_1D_BSPLINE), INTENT(INOUT) :: self

   END SUBROUTINE FREE_BASIS_1D_BSPLINE
   ! ...................................................

   ! ...................................................
   SUBROUTINE RESET_BASIS_1D_BSPLINE(self)
    IMPLICIT NONE
     CLASS(DEF_BASIS_1D_BSPLINE), INTENT(INOUT) :: self

   END SUBROUTINE RESET_BASIS_1D_BSPLINE
   ! ...................................................

   ! ...................................................
   SUBROUTINE UPDATE_BASIS_1D_BSPLINE(self)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D_BSPLINE), INTENT(INOUT) :: self
     ! LOCAL
     INTEGER       :: li_ig
     REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1) :: B     
     REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1) :: B_s   
     REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1) :: B_ss  
     INTEGER :: li_n_points
     INTEGER :: i_element
     REAL(SPI_RK) :: lr_x

     DO i_element=1, self % ptr_mesh % n_elements 
!     DO i_element=self % ptr_mesh % n_elements, self % ptr_mesh % n_elements 
        DO li_ig = 1, self % ptr_quad % oi_n_points
           lr_x = self % ptr_mesh % opr_points(i_element,li_ig)
!           print *, lr_x
           CALL EVALUATE_BASIS_1D_BSPLINE(self, lr_x, B, B_s, B_ss)
           print *, lr_x, B

           self % TestfT_0  (i_element, li_ig, :) = B   (:)
           self % TestfT_p  (i_element, li_ig, :) = B_s (:)
           self % TestfT_pp (i_element, li_ig, :) = B_ss(:)
        END DO
     END DO

   END SUBROUTINE UPDATE_BASIS_1D_BSPLINE
   ! ...................................................

   ! ...................................................
  SUBROUTINE EVALUATE_BASIS_1D_BSPLINE(self, s, B, B_s, B_ss)
  IMPLICIT NONE
    CLASS(DEF_BASIS_1D_BSPLINE), INTENT(INOUT) :: self
    !< s-coordinate in the element
    REAL(KIND=SPI_RK), INTENT(IN)  :: s          
    !< Basis functions
    REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1), INTENT(INOUT) :: B     
    !< Basis functions derived with respect to s
    REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1), INTENT(INOUT) :: B_s   
    !< Basis functions derived two times with respect to s
    REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1), INTENT(INOUT) :: B_ss  
    ! LOCAL
    INTEGER, PARAMETER :: N_DERIV = 2
    REAL(SPI_RK), DIMENSION(0:self % ptr_mesh % oi_p, 0:N_DERIV) :: lpr_dbatx
    REAL(SPI_RK), DIMENSION(0:self % ptr_mesh % oi_n + self % ptr_mesh % oi_p) :: lpr_knots
    INTEGER :: li_span

    lpr_knots(0:self % ptr_mesh % oi_n + self % ptr_mesh % oi_p) = &
            & self % ptr_mesh % opr_knots(1:self % ptr_mesh % oi_n + self % ptr_mesh % oi_p + 1)

    li_span = -1
    CALL EvalBasisFunsDers( self % ptr_mesh % oi_p, &
                          & self % ptr_mesh % oi_n + self % ptr_mesh % oi_p, &
                          & lpr_knots, &
                          & s, N_DERIV, li_span, lpr_dbatx)

    B(:)    = lpr_dbatx(:,0)
    B_s(:)  = lpr_dbatx(:,1)
    B_ss(:) = lpr_dbatx(:,2)

  END SUBROUTINE EVALUATE_BASIS_1D_BSPLINE
   ! ...................................................

END MODULE SPI_BASIS_BSPLINE
