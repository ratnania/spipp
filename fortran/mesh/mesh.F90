!# -*- coding: utf8 -*-
MODULE SPI_MESH 
  USE SPI_MESH_DEF
  USE SPI_MESH_BSPLINE
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  IMPLICIT NONE

CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_MESH(self &
                  & , quad_u &
                  & , quad_v &
                  & , n_u &
                  & , n_v &
                  & , p_u &
                  & , p_v &
                  & , type_bc_u &
                  & , type_bc_v &
                  & , knots_u &
                  & , knots_v &
                  & , control_points_1d &
                  & , control_points_2d &
                  & )
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_MESH_ABSTRACT)       , INTENT(INOUT) :: self 
     CLASS(DEF_QUADRATURE_ABSTRACT), TARGET, OPTIONAL, INTENT(INOUT) :: quad_u 
     CLASS(DEF_QUADRATURE_ABSTRACT), TARGET, OPTIONAL, INTENT(INOUT) :: quad_v 
     INTEGER               , OPTIONAL, INTENT(IN)   :: n_u
     INTEGER               , OPTIONAL, INTENT(IN)   :: n_v
     INTEGER               , OPTIONAL, INTENT(IN)   :: p_u
     INTEGER               , OPTIONAL, INTENT(IN)   :: p_v
     INTEGER               , OPTIONAL, INTENT(IN)   :: type_bc_u
     INTEGER               , OPTIONAL, INTENT(IN)   :: type_bc_v
     real(SPI_RK), dimension (:), OPTIONAL, INTENT(IN) :: knots_u
     real(SPI_RK), dimension (:), OPTIONAL, INTENT(IN) :: knots_v
     real(SPI_RK), dimension (:,:), OPTIONAL, INTENT(IN) :: control_points_1d
     real(SPI_RK), dimension (:,:,:), OPTIONAL, INTENT(IN) :: control_points_2d
     ! LOCAL
     INTEGER :: li_err     

     SELECT TYPE (self)
     CLASS IS (DEF_MESH_1D_BSPLINE)
        IF ((PRESENT(n_u)) .AND. (PRESENT(p_u)) .AND. (PRESENT(quad_u))) THEN
           CALL CREATE_MESH_BSPLINES_1D(self, &
                   & n_u, &
                   & p_u, & 
                   & type_bc=type_bc_u, &
                   & knots=knots_u, &
                   & control_points=control_points_1d&
                   & )
           CALL CREATE_MESH_POINTS_1D(self, quad_u) 
        END IF

     CLASS DEFAULT
        STOP 'CREATE_MESH: unexpected type for self object!'
     END SELECT

  END SUBROUTINE CREATE_MESH
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_POINTS_1D(self, ao_quad)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     CLASS(DEF_QUADRATURE_ABSTRACT), TARGET , INTENT(INOUT) :: ao_quad 
     ! LOCAL
     INTEGER :: i_element
     INTEGER :: li_ig 
     REAL(SPI_RK) :: lr_a
     REAL(SPI_RK) :: lr_b
     REAL(SPI_RK) :: lr_x

     self % ptr_quad => ao_quad

     ALLOCATE(self % opr_points(self % n_elements, self % ptr_quad % oi_n_points))

     DO i_element=1, self % n_elements 
        lr_a = self % opr_grid(i_element)
        lr_b = self % opr_grid(i_element+1)

        DO li_ig = 1, self % ptr_quad % oi_n_points
           lr_x = lr_a + self % ptr_quad % opr_points(li_ig) * (lr_b - lr_a )
           self % opr_points(i_element, li_ig) = lr_x
        END DO
     END DO

  END SUBROUTINE CREATE_MESH_POINTS_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE FREE_MESH(self)
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_MESH_ABSTRACT)       , INTENT(INOUT) :: self 
     ! LOCAL

     DEALLOCATE ( self % opr_grid) 

     SELECT TYPE (self)
     CLASS IS (DEF_MESH_1D_BSPLINE)
        CALL FREE_MESH_BSPLINES_1D(self)

     CLASS DEFAULT
        STOP 'FREE_MESH: unexpected type for self object!'
     END SELECT

  END SUBROUTINE FREE_MESH
  ! .........................................................

  ! .........................................................
END MODULE SPI_MESH 
