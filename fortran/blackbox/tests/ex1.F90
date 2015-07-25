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
USE SPI_BLACKBOX_DEF
USE SPI_BLACKBOX
implicit none
   ! LOCAL
   TYPE(DEF_QUADRATURE_1D), TARGET :: lo_quad
   TYPE(DEF_MESH_1D_BSPLINE), TARGET :: lo_mesh
   TYPE(DEF_BASIS_1D_BSPLINE), TARGET :: lo_basis
   TYPE(DEF_BLACKBOX_1D_BSPLINE), TARGET :: lo_bbox
   ! ... number of internal knots is = N - P - 1
   INTEGER, PARAMETER :: N = 5 
   INTEGER, PARAMETER :: P = 3
   INTEGER, PARAMETER :: K = 3 
   INTEGER :: li_elmt
   REAL(SPI_RK) :: lr_a
   REAL(SPI_RK) :: lr_b

   CALL CREATE_QUADRATURE(lo_quad, SPI_QUADRATURES_LEGENDRE, K)
   CALL CREATE_MESH(lo_mesh, lo_quad, ai_n=N, ai_p=P, ai_type_bc=SPI_BC_PERIODIC) 
   CALL CREATE_BASIS(lo_basis, lo_mesh) 
   CALL CREATE_BLACKBOX(lo_bbox, lo_basis, lo_quad)

   PRINT *, ">>> knots"
   PRINT *, lo_mesh % opr_knot

   CALL RESET_BASIS(lo_basis) 

   li_elmt = 1
   CALL UPDATE_BASIS(lo_basis, li_elmt) 

   CALL BLACKBOX_RESET_POSITION(lo_bbox) 

   lr_a = 0.0 ; lr_b = 1.0
   CALL BLACKBOX_COMPUTE_METRIC(lo_bbox, lr_a, lr_b) 

   CALL FREE_QUADRATURE(lo_quad) 
   CALL FREE_MESH(lo_mesh) 
   CALL FREE_BASIS(lo_basis) 
   CALL FREE_BLACKBOX(lo_bbox) 

end subroutine test1
! ............................................
   
! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
