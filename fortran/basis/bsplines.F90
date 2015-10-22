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
          & EVALUATE_BASIS_ON_QUADRATURE_POINTS_1D_BSPLINE, &
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

     ALLOCATE(self % B_0 (n_elements, n_points, self % oi_n_order))
     ALLOCATE(self % B_s1 (n_elements, n_points, self % oi_n_order))
     ALLOCATE(self % B_s1s1(n_elements, n_points, self % oi_n_order))

     ALLOCATE ( self % opi_i_begin( n_elements) ) 
     self % opi_i_begin = SPI_INT_DEFAULT

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
   SUBROUTINE EVALUATE_BASIS_ON_QUADRATURE_POINTS_1D_BSPLINE(self)
   IMPLICIT NONE
     CLASS(DEF_BASIS_1D_BSPLINE), INTENT(INOUT) :: self
     ! LOCAL
     INTEGER       :: li_ig
     REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1) :: B     
     REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1) :: B_s   
     REAL(KIND=SPI_RK), DIMENSION(self % ptr_mesh % oi_p + 1) :: B_ss  
     INTEGER :: li_n_points
     INTEGER :: i_element
     INTEGER :: i_span
     REAL(SPI_RK) :: lr_x

     DO i_element=1, self % ptr_mesh % n_elements 
        DO li_ig = 1, self % ptr_quad % oi_n_points
           lr_x = self % ptr_mesh % opr_points(i_element,li_ig)
           CALL EVALUATE_BASIS_1D_BSPLINE(self, lr_x, B, B_s, B_ss, i_span)
!           print *, lr_x, B

           self % B_0    (i_element, li_ig, :) = B   (:)
           self % B_s1   (i_element, li_ig, :) = B_s (:)
           self % B_s1s1 (i_element, li_ig, :) = B_ss(:)

           self % opi_i_begin (i_element) = i_span
        END DO
     END DO

   END SUBROUTINE EVALUATE_BASIS_ON_QUADRATURE_POINTS_1D_BSPLINE
   ! ...................................................

   ! ...................................................
  SUBROUTINE EVALUATE_BASIS_1D_BSPLINE(self, s, B, B_s, B_ss, span)
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
    INTEGER, INTENT(INOUT) :: span
    ! LOCAL
    INTEGER, PARAMETER :: N_DERIV = 2
    REAL(SPI_RK), DIMENSION(0:self % ptr_mesh % oi_p, 0:N_DERIV) :: lpr_dbatx
    REAL(SPI_RK), DIMENSION(0:self % ptr_mesh % oi_n + self % ptr_mesh % oi_p) :: lpr_knots
    INTEGER :: i_element

    lpr_knots(0:self % ptr_mesh % oi_n + self % ptr_mesh % oi_p) = &
            & self % ptr_mesh % opr_knots(1:self % ptr_mesh % oi_n + self % ptr_mesh % oi_p + 1)

    span = -1
    CALL EvalBasisFunsDers( self % ptr_mesh % oi_p, &
                          & self % ptr_mesh % oi_n + self % ptr_mesh % oi_p, &
                          & lpr_knots, &
                          & s, N_DERIV, span, lpr_dbatx)

    B(:)    = lpr_dbatx(:,0)
    B_s(:)  = lpr_dbatx(:,1)
    B_ss(:) = lpr_dbatx(:,2)

  END SUBROUTINE EVALUATE_BASIS_1D_BSPLINE
   ! ...................................................

END MODULE SPI_BASIS_BSPLINE
