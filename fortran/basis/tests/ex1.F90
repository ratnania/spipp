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
!   INTEGER, PARAMETER :: N = P + 5 
   INTEGER, PARAMETER :: P = 1 
   INTEGER, PARAMETER :: N = P + 4 
   REAL(SPI_RK), DIMENSION(N+P+1) :: KNOTS
   INTEGER :: li_elmt
   INTEGER :: i, e

   KNOTS(1:P+1) = 0.0
   KNOTS(P+2) = 0.25
   KNOTS(P+3) = 0.5
!   KNOTS(P+4) = 0.5
   KNOTS(N)   = 0.75
   KNOTS(N+1:N+P+1) = 1.0

   CALL CREATE_QUADRATURE(lo_quad, SPI_QUADRATURES_LEGENDRE, P) 
   CALL CREATE_MESH(lo_mesh, lo_quad, N, P, apr_knots=KNOTS) 
   CALL CREATE_BASIS(lo_basis, lo_mesh, lo_quad) 

   CALL RESET_BASIS(lo_basis) 


   PRINT *, ">>> points"
   DO e = 1, 4
      PRINT *, lo_mesh % opr_points(e,:)
   END DO


   CALL UPDATE_BASIS(lo_basis)
   
!   PRINT *, SIZE(lo_basis % TestfT_0,1)
!   PRINT *, SIZE(lo_basis % TestfT_0,2)
!   PRINT *, SIZE(lo_basis % TestfT_0,3)
!   DO e = 1, 4
!      PRINT *, ">>> element ", e
!      DO i = 1, P+1 
!         PRINT *, lo_basis % TestfT_0(e,i,:)
!      END DO
!   END DO
!   PRINT *, lo_basis % TestfT_P
!   PRINT *, lo_basis % TestfT_PP

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
   
