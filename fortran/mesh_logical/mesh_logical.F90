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
END MODULE SPI_MESH_LOGICAL 
