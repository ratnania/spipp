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
USE SPI_SPACE_DEF
USE SPI_SPACE
implicit none
   ! LOCAL
   TYPE(DEF_MESH_1D_BSPLINE), TARGET :: lo_mesh
   TYPE(DEF_SPACE_1D_BSPLINE), TARGET :: lo_space
   integer  :: li_n, li_nel
   real(SPI_RK),dimension(:),pointer :: lpr_x,lpr_y
   integer  :: li_err,li_flag,li_i,li_j
   character(len=14) :: ls_file		
   ! ... number of internal knots is = N - P - 1
   INTEGER, PARAMETER :: N = 5 
   INTEGER, PARAMETER :: P = 3
   INTEGER, PARAMETER :: K = 3 

   ls_file = "ex_2_matrix.mm"
   
   CALL CREATE_MESH(lo_mesh, ai_n=N, ai_p=P, ai_type_bc=SPI_BC_PERIODIC) 
   CALL CREATE_SPACE(lo_space, lo_mesh, SPI_QUADRATURES_LEGENDRE, ai_k=K) 
   
   ! ... free object
   CALL FREE_MESH(lo_mesh) 
   CALL FREE_SPACE(lo_space)
   ! ...
   
end subroutine test1
! ............................................

! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
