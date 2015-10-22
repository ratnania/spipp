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
implicit none
   ! LOCAL
   TYPE(DEF_MESH_1D_BSPLINE), TARGET :: lo_mesh
   ! ... number of internal knots is = N - P - 1
   INTEGER, PARAMETER :: N = 8 
   INTEGER, PARAMETER :: P = 3
   REAL(SPI_RK), DIMENSION(N+P+1) :: KNOTS
   TYPE(DEF_QUADRATURE_1D), TARGET :: lo_quad
   INTEGER, PARAMETER :: K = 3 
   INTEGER :: e

   KNOTS(1:P+1) = 0.0
   KNOTS(P+2) = 0.25
   KNOTS(P+3) = 0.5
   KNOTS(P+4) = 0.5
   KNOTS(P+5) = 0.75
   KNOTS(N+1:N+P+1) = 1.0

   CALL CREATE_QUADRATURE(lo_quad, SPI_QUADRATURES_LEGENDRE, K)
   CALL CREATE_MESH(lo_mesh, lo_quad, N, P, apr_knots=KNOTS) 

   PRINT *, ">>> N, P :", lo_mesh % oi_N, lo_mesh % oi_P

   PRINT *, ">>> knots"
   PRINT *, lo_mesh % opr_knots

   PRINT *, ">>> grid"
   PRINT *, lo_mesh % opr_grid

   PRINT *, ">>> control points"
   PRINT *, lo_mesh % opr_control_points

   PRINT *, ">>> points"
   DO e = 1, 4
      PRINT *, lo_mesh % opr_points(e,:)
   END DO

   CALL FREE_MESH(lo_mesh) 
end subroutine test1
! ............................................
   
! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
