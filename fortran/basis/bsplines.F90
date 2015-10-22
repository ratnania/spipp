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

     self % ptr_mesh => ao_mesh
     self % ptr_knot => ao_mesh % opr_knot 

     self % oi_n_order = self % ptr_mesh % oi_p + 1 

     ALLOCATE(self % TestfT_0 (self % ptr_quad % oi_n_points, self % ptr_mesh % oi_p + 1, 1))
     ALLOCATE(self % TestfT_p (self % ptr_quad % oi_n_points, self % ptr_mesh % oi_p + 1, 1))
     ALLOCATE(self % TestfT_pp(self % ptr_quad % oi_n_points, self % ptr_mesh % oi_p + 1, 1))

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
     INTEGER       :: i_point
     INTEGER       :: i_element
     REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1) :: B     
     REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1) :: B_s   
     REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1) :: B_ss  
     INTEGER :: n_points

     n_points = SIZE(apr_points, 1)
    
     DO i_element = 1, self % ptr_mesh % oi_n_elmts 
        DO i_point = 1, n_points
           CALL EVALUATE_BASIS_1D_BSPLINE(self, apr_points(li_ig), B, B_s, B_ss)
           self % TestfT_0  (i_element, :, i_point) = B   (:)
           self % TestfT_p  (i_element, :, i_point) = B_s (:)
           self % TestfT_pp (i_element, :, i_point) = B_ss(:)
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
    INTEGER :: li_span

    CALL FindSpan( self % ptr_mesh % oi_p, &
                 & self % ptr_mesh % oi_n, &
                 & self % ptr_knot       , &
                 & s, li_span)
    
    CALL EvalBasisFunsDers( self % ptr_mesh % oi_p, &
                          & self % ptr_mesh % oi_n, &
                          & self % ptr_knot       , &
                          & s, N_DERIV, li_span, lpr_dbatx)

    B(:)    = lpr_dbatx(:,0)
    B_s(:)  = lpr_dbatx(:,1)
    B_ss(:) = lpr_dbatx(:,2)

  END SUBROUTINE EVALUATE_BASIS_1D_BSPLINE
   ! ...................................................

END MODULE SPI_BASIS_BSPLINE
