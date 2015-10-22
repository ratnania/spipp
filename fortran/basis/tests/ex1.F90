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
   INTEGER, PARAMETER :: N = 5 
   INTEGER, PARAMETER :: P = 3
   INTEGER, PARAMETER :: K = 3 
   INTEGER :: li_elmt
   REAL(SPI_RK), DIMENSION(3) :: lpr_points

   CALL CREATE_MESH(lo_mesh, ai_n=N, ai_p=P, ai_type_bc=SPI_BC_PERIODIC) 
   CALL CREATE_QUADRATURE(lo_quad, SPI_QUADRATURES_LEGENDRE, K) 
   CALL CREATE_BASIS(lo_basis, lo_mesh, lo_quad) 

   CALL RESET_BASIS(lo_basis) 

   lpr_points(1) = 0.2
   lpr_points(2) = 0.4
   lpr_points(3) = 0.6
   CALL UPDATE_BASIS(lo_basis, lpr_points) 
   
   PRINT *, lo_basis % TestfT_0
   PRINT *, lo_basis % TestfT_P
   PRINT *, lo_basis % TestfT_PP

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
   
