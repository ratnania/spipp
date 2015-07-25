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
USE SPI_QUADRATURES_DEF
USE SPI_QUADRATURES
USE SPI_MESH_DEF
USE SPI_MESH
USE SPI_BASIS_DEF
USE SPI_BASIS
implicit none
   ! LOCAL
   TYPE(DEF_MESH_1D_BSPLINE), TARGET :: lo_mesh
   ! ... number of internal knots is = N - P - 1
   INTEGER, PARAMETER :: N = 5 
   INTEGER, PARAMETER :: P = 3

   CALL CREATE_MESH(lo_mesh, ai_n=N, ai_p=P, ai_type_bc=SPI_BC_PERIODIC) 
!   CALL CREATE_MESH(lo_mesh, ai_n=N, ai_p=P, ai_type_bc=SPI_BC_DIRICHLET_HOMOGEN) 

   PRINT *, ">>> knots"
   PRINT *, lo_mesh % opr_knot

   CALL FREE_MESH(lo_mesh) 
end subroutine test1
! ............................................
   
! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
