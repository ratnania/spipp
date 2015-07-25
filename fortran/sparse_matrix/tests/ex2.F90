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
subroutine test2 ()
USE SPI_GLOBAL_DEF
USE SPI_GLOBAL
USE SPI_QUADRATURES_DEF
USE SPI_QUADRATURES
USE SPI_MESH_DEF
USE SPI_MESH
USE SPI_NUMBERING_DEF
USE SPI_NUMBERING
USE SPI_SPARSE_MATRIX_DEF
USE SPI_SPARSE_MATRIX
implicit none
   ! LOCAL
   TYPE(DEF_QUADRATURE_1D), TARGET :: lo_quad
   TYPE(DEF_MESH_1D_BSPLINE), TARGET :: lo_mesh
   TYPE(DEF_NUMBERING_1D_BSPLINE), TARGET :: lo_numbering
   type(DEF_MATRIX_CSR) :: lo_A
   integer  :: li_n, li_nel
   real(SPI_RK),dimension(:),pointer :: lpr_x,lpr_y
   integer  :: li_err,li_flag,li_i,li_j
   character(len=14) :: ls_file		
   ! ... number of internal knots is = N - P - 1
   INTEGER, PARAMETER :: N = 5 
   INTEGER, PARAMETER :: P = 3
   INTEGER, PARAMETER :: K = 3 

   ls_file = "ex_2_matrix.mm"
   
   CALL CREATE_QUADRATURE(lo_quad, SPI_QUADRATURES_LEGENDRE, K)
   CALL CREATE_MESH(lo_mesh, lo_quad, ai_n=N, ai_p=P, ai_type_bc=SPI_BC_PERIODIC) 
   CALL CREATE_NUMBERING(lo_numbering, lo_mesh) 

   call create_sparse_matrix( lo_A, ai_nR=li_n, ai_nC=li_n, ai_nel=li_nel )

   PRINT *, ">>> LM"
   PRINT *, lo_numbering % opi_LM

   print *, "%%%%%"
   CALL initialize_sparse_matrix_with_LM(lo_A, lo_mesh % oi_n_elmts, &
           & lo_numbering % opi_LM, lo_mesh % oi_nen, &
           & lo_numbering % opi_LM, lo_mesh % oi_nen)
   print *, "%%%%%"
   
   ! ... save the matrix in the MM format
   CALL save_sparse_matrix(lo_A, ls_file, SPI_MATRIX_OUTPUT_FORMAT_MM)  
   ! ...
   
   ! ... free object
   CALL FREE_QUADRATURE(lo_quad) 
   CALL FREE_NUMBERING(lo_numbering) 
   CALL FREE_MESH(lo_mesh) 
   call free_sparse_matrix ( lo_A )
   ! ...
   
end subroutine test2
! ............................................

! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test2()

END PROGRAM Main
! ............................................
   
