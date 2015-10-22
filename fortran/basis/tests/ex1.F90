!
! File: ex1.f90
!
!
! Usage:
!   > ./ex1 
!
! Authors:
!   Ahmed RATNANI  - ratnaniahmed@gmail.com
!

! ............................................
subroutine test1 ()
USE SPI_GLOBAL_DEF
USE SPI_GLOBAL
USE SPI_MESH_DEF
USE SPI_MESH
USE SPI_QUADRATURES_DEF
USE SPI_QUADRATURES
USE SPI_BASIS_DEF
USE SPI_BASIS
implicit none
   ! LOCAL
   TYPE(DEF_MESH_1D_BSPLINE), TARGET :: lo_mesh
   TYPE(DEF_QUADRATURE_1D), TARGET :: lo_quad
   TYPE(DEF_BASIS_1D_BSPLINE), TARGET :: lo_basis
   ! ... number of internal knots is = N - P - 1
!   INTEGER, PARAMETER :: P = 3
   INTEGER, PARAMETER :: P = 2 
   INTEGER, PARAMETER :: N = P + 1 + 3 
   REAL(SPI_RK), DIMENSION(N+P+1) :: KNOTS
   INTEGER :: li_elmt
   INTEGER :: i, e, i_point

   KNOTS(1:P+1) = 0.0
   DO i =1, N - P - 1 
      KNOTS(i+P+1) = i * 1.0 / (N-P) 
   END DO
   KNOTS(N+1:N+P+1) = 1.0

   CALL CREATE_QUADRATURE(lo_quad, SPI_QUADRATURES_LEGENDRE, P) 
   CALL CREATE_MESH(lo_mesh, quad_u=lo_quad, n_u=N, p_u=P, knots_u=KNOTS) 
   CALL CREATE_BASIS(lo_basis, lo_mesh, lo_quad) 

   PRINT *, ">>> points"
   DO e = 1, lo_mesh % n_elements 
      print *, "========"
      PRINT *, lo_mesh % opr_points(e, :)
      print *, "   --- "
      DO i_point = 1, P+1
        PRINT *, lo_basis % B_0(e, i_point, :)
      END DO
   END DO

   CALL FREE_BASIS(lo_basis) 
   CALL FREE_MESH(lo_mesh) 
   CALL FREE_QUADRATURE(lo_quad) 
end subroutine test1
! ............................................
   
! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
