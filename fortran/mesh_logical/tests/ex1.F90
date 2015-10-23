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
USE SPI_MESH_LOGICAL_DEF
USE SPI_MESH_LOGICAL
USE SPI_QUADRATURES_DEF
USE SPI_QUADRATURES
implicit none
   ! LOCAL
   TYPE(DEF_MESH_LOGICAL_BSPLINE), TARGET :: mesh_logical
   ! ... number of internal knots is = N - P - 1
   INTEGER, PARAMETER :: N = 8 
   INTEGER, PARAMETER :: P = 3
   REAL(SPI_RK), DIMENSION(N+P+1) :: KNOTS
   TYPE(DEF_QUADRATURE_1D), TARGET :: quad
   INTEGER, PARAMETER :: K = 3 
   INTEGER :: e

   KNOTS(1:P+1) = 0.0
   KNOTS(P+2) = 0.25
   KNOTS(P+3) = 0.5
   KNOTS(P+4) = 0.5
   KNOTS(P+5) = 0.75
   KNOTS(N+1:N+P+1) = 1.0

   CALL CREATE_QUADRATURE(quad, SPI_QUADRATURES_LEGENDRE, K)
   CALL CREATE_MESH_LOGICAL(mesh_logical, quad=quad, n=N, p=P, knots=KNOTS) 

   PRINT *, ">>> N, P :", mesh_logical % n, mesh_logical % p

   PRINT *, ">>> knots"
   PRINT *, mesh_logical % knots

   PRINT *, ">>> grid"
   PRINT *, mesh_logical % grid

!   PRINT *, ">>> points"
!   DO e = 1, 4
!      PRINT *, mesh_logical % points(e,:)
!   END DO

   CALL FREE_MESH_LOGICAL(mesh_logical) 
end subroutine test1
! ............................................
   
! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
