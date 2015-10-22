!# -*- coding: utf8 -*-
MODULE SPI_BLACKBOX_BSPLINE
  USE SPI_BLACKBOX_DEF
  USE SPI_QUADRATURES_DEF
  USE SPI_BASIS_DEF
  IMPLICIT NONE

    CONTAINS

  ! ........................................................
  SUBROUTINE CREATE_BLACKBOX_1D_BSPLINE(self, ao_mesh, ao_basis)
  implicit none
     type(DEF_BLACKBOX_1D_BSPLINE)        , INTENT(INOUT) :: self
     TYPE(DEF_MESH_1D_BSPLINE)   , TARGET, INTENT(IN)    :: ao_mesh
     TYPE(DEF_BASIS_1D_BSPLINE)   , TARGET, INTENT(IN)    :: ao_basis
     ! LOCAL VARIABLES
     integer  :: li_ref

     self % ptr_basis => ao_basis
     self % ptr_mesh => ao_mesh

  END SUBROUTINE CREATE_BLACKBOX_1D_BSPLINE
  ! ........................................................

  ! ........................................................
  SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_BSPLINE(self)
  implicit none
     type(DEF_BLACKBOX_1D_BSPLINE), INTENT(INOUT) :: self
     ! LOCAL VARIABLES
     REAL(SPI_RK) :: lr_a
     REAL(SPI_RK) :: lr_b
     INTEGER :: i_element
     INTEGER :: i_point
     INTEGER :: i

     DO i_point = 1, self % ptr_quad %  oi_n_points
        DO i_element = 1, self % ptr_mesh % n_elements 
           lr_a = self % ptr_mesh % opr_grid(i_element)
           lr_b = self % ptr_mesh % opr_grid(i_element+1)

           self % Vol(i_element, i_point) = ABS( lr_b - lr_a )
           self % wVol(i_element, i_point) = self % Vol(i_element, i_point) * self % ptr_quad % opr_weights(i_point)
        END DO
     END DO

     self % Jacobians = 0.0d0
     DO i = 1, self % n_dim 
        DO i_point = 1, self % ptr_quad %  oi_n_points
           DO i_element = 1, self % ptr_mesh % n_elements 
              self % Jacobians(i_element, i_point) = self % Jacobians(i_element, i_point) &
                      & + self % Xp_s1(i_element, i_point, i)**2
           END DO
        END DO
     END DO
     self % Jacobians = SQRT(self % Jacobians)
  END SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_BSPLINE
  ! .........................................................

  ! ........................................................        
  SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_BSPLINE(self)
  IMPLICIT NONE
     CLASS(DEF_BLACKBOX_1D_BSPLINE) :: self
     INTEGER :: ai_i
     ! LOCAL
     INTEGER :: i_element
     INTEGER :: i_basis
     INTEGER :: i_point
     INTEGER :: i_ctrl_pt

     self % Xp_0    = 0.0d0
     self % Xp_s1   = 0.0d0
     self % Xp_s1s1 = 0.0d0

     DO i_point = 1, self % ptr_quad %  oi_n_points
        DO i_basis = 1, self % ptr_mesh % oi_p + 1 
           DO i_element = 1, self % ptr_mesh % n_elements 
              i_ctrl_pt = i_basis + self % ptr_basis % opi_i_begin(i_element) - self % ptr_mesh % oi_p 

              self % Xp_0(i_element, i_point, :) = self % Xp_0(i_element, i_point, :) &
                      & + self % ptr_basis % B_0(i_element, i_point, i_basis) * self % ptr_mesh % control_points(i_ctrl_pt, :)

              self % Xp_s1(i_element, i_point, :) = self % Xp_s1(i_element, i_point, :) &
                      & + self % ptr_basis % B_s1(i_element, i_point, i_basis) * self % ptr_mesh %control_points(i_ctrl_pt, :)

              self % Xp_s1s1(i_element, i_point, :) = self % Xp_s1s1(i_element, i_point, :) &
                      & + self % ptr_basis % B_s1s1(i_element, i_point, i_basis) * self % ptr_mesh %control_points(i_ctrl_pt, :)
           END DO
        END DO
     END DO

     DO i_point = 1, self % ptr_quad %  oi_n_points
        DO i_element = 1, self % ptr_mesh % n_elements 
           self % Xp_x1(i_element, i_point, :) = self % Xp_s1(i_element, i_point, :) /  self % Vol(i_element, i_point) 
           self % Xp_x1x1(i_element, i_point, :) = self % Xp_s1s1(i_element, i_point, :) /  self % Vol(i_element, i_point)**2
        END DO
     END DO

  END SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_BSPLINE
  ! ..........................................................  

END MODULE SPI_BLACKBOX_BSPLINE
