!# -*- coding: utf8 -*-
MODULE SPI_MESH_LOGICAL 
  USE SPI_MESH_LOGICAL_DEF
  USE SPI_MESH_LOGICAL_BSPLINE
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  IMPLICIT NONE

CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_MESH_LOGICAL(self, quad, n, p, type_bc, knots)
  IMPLICIT NONE
     CLASS(DEF_MESH_LOGICAL_ABSTRACT)  , INTENT(INOUT) :: self
     CLASS(DEF_QUADRATURE_ABSTRACT), TARGET, INTENT(INOUT) :: quad 
     INTEGER               , OPTIONAL, INTENT(IN)   :: n
     INTEGER               , OPTIONAL, INTENT(IN)   :: p
     INTEGER               , OPTIONAL, INTENT(IN)   :: type_bc
     real(SPI_RK), dimension (:), OPTIONAL, INTENT(IN) :: knots
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_MESH_LOGICAL_BSPLINE)
        IF ((PRESENT(n)) .AND. (PRESENT(p))) THEN
           CALL CREATE_MESH_LOGICAL_BSPLINE(self, &
                   & n, &
                   & p, & 
                   & type_bc=type_bc, &
                   & knots=knots &
                   & )
           CALL SET_MESH_LOGICAL_QUADRATURE(self, quad) 
           CALL INTIALIZE_MESH_LOGICAL_POINTS(self)
        ELSE
           STOP 'CREATE_MESH_LOGICAL: wrong arguments!'
        END IF
     CLASS DEFAULT
        STOP 'CREATE_MESH_LOGICAL: unexpected type for self object!'
     END SELECT
     ! ...
  END SUBROUTINE CREATE_MESH_LOGICAL
  ! .........................................................

  ! .........................................................
  SUBROUTINE FREE_MESH_LOGICAL(self)
  IMPLICIT NONE
     CLASS(DEF_MESH_LOGICAL_ABSTRACT)      , INTENT(INOUT) :: self
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_MESH_LOGICAL_BSPLINE)
        CALL FREE_MESH_LOGICAL_BSPLINE(self)
     CLASS DEFAULT
        STOP 'FREE_MESH_LOGICAL: unexpected type for self object!'
     END SELECT
     ! ...

  END SUBROUTINE FREE_MESH_LOGICAL
  ! .........................................................

  ! .........................................................
  SUBROUTINE SET_MESH_LOGICAL_QUADRATURE(self, quad)
  IMPLICIT NONE
     CLASS(DEF_MESH_LOGICAL_BSPLINE)  , INTENT(INOUT) :: self
     CLASS(DEF_QUADRATURE_ABSTRACT), TARGET, INTENT(INOUT) :: quad 
     ! LOCAL

     self % ptr_quad => quad

  END SUBROUTINE SET_MESH_LOGICAL_QUADRATURE
  ! .........................................................

  ! .........................................................
  SUBROUTINE INTIALIZE_MESH_LOGICAL_POINTS(self)
  IMPLICIT NONE
     CLASS(DEF_MESH_LOGICAL_ABSTRACT)      , INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: i_element
     INTEGER :: li_ig 
     REAL(SPI_RK) :: lr_a
     REAL(SPI_RK) :: lr_b
     REAL(SPI_RK) :: lr_x

     ALLOCATE(self % points(self % n_elements, self % ptr_quad % oi_n_points))

     DO i_element=1, self % n_elements 
        lr_a = self % grid(i_element)
        lr_b = self % grid(i_element+1)

        DO li_ig = 1, self % ptr_quad % oi_n_points
           lr_x = lr_a + self % ptr_quad % points(li_ig) * (lr_b - lr_a )
           self % points(i_element, li_ig) = lr_x
        END DO
     END DO

  END SUBROUTINE INTIALIZE_MESH_LOGICAL_POINTS
  ! .........................................................

END MODULE SPI_MESH_LOGICAL 
