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
USE SPI_GLOBAL
USE SPI_SPARSE_MATRIX_DEF
USE SPI_SPARSE_MATRIX
implicit none
   ! LOCAL
   type(DEF_MATRIX_CSR) :: lo_A
   integer  :: li_n, li_nnz
   real(SPI_RK),dimension(:),pointer :: lpr_x,lpr_y
   integer  :: li_err,li_flag,li_i,li_j
   character(len=14) :: ls_file		

   ls_file = "ex_1_matrix.mm"
   
   li_n = 8
   li_nnz = 18
   
   allocate(lpr_x(li_n),stat=li_err)
   if (li_err.ne.0) li_flag=10	
   
   allocate(lpr_y(li_n),stat=li_err)
   if (li_err.ne.0) li_flag=10		
   
   call create_sparse_matrix( lo_A, ai_nR=li_n, ai_nC=li_n, ai_nnz=li_nnz )
   
   lo_A%opi_ia(1) =1
   lo_A%opi_ia(2) =5	
   lo_A%opi_ia(3) =	8	
   lo_A%opi_ia(4) =	10	
   lo_A%opi_ia(5) =	12
   lo_A%opi_ia(6) =	15	
   lo_A%opi_ia(7) =	17
   lo_A%opi_ia(8) =	18
   lo_A%opi_ia(9) =	19
   
   lo_A%opi_ja(1) = 1
   lo_A%opi_ja(2) =	3
   lo_A%opi_ja(3) =	6	
   lo_A%opi_ja(4) =	7	
   lo_A%opi_ja(5) =	2
   lo_A%opi_ja(6) =	3	
   lo_A%opi_ja(7) =	5
   lo_A%opi_ja(8) =	3
   lo_A%opi_ja(9) =8
   lo_A%opi_ja(10) =4	
   lo_A%opi_ja(11) =7		
   lo_A%opi_ja(12) =5		
   lo_A%opi_ja(13) =6	
   lo_A%opi_ja(14) =7		
   lo_A%opi_ja(15) =6	
   lo_A%opi_ja(16) =8
   lo_A%opi_ja(17) =7
   lo_A%opi_ja(18) =8	
   
   lo_A%opr_a(1) =7
   lo_A%opr_a(2) =	1
   lo_A%opr_a(3) =	2	
   lo_A%opr_a(4) =	7	
   lo_A%opr_a(5) =	-4
   lo_A%opr_a(6) =	8	
   lo_A%opr_a(7) =	2
   lo_A%opr_a(8) =	1
   lo_A%opr_a(9) =5
   lo_A%opr_a(10) =7	
   lo_A%opr_a(11) =9		
   lo_A%opr_a(12) =5		
   lo_A%opr_a(13) =-1	
   lo_A%opr_a(14) =5		
   lo_A%opr_a(15) =0	
   lo_A%opr_a(16) =5
   lo_A%opr_a(17) =11
   lo_A%opr_a(18) =5	
   
   do li_i=1,li_n
      lpr_x(li_i) = 1.0*li_i
   enddo
   
   lpr_y(1:li_n) = 0.0	
   
   ! ... save the matrix in the MM format
   CALL save_sparse_matrix(lo_A, ls_file, SPI_MATRIX_OUTPUT_FORMAT_MM)  
   ! ...

   ! ... matrix vector product 
   call mult_sparse_matrix_vector ( lo_A, lpr_x, lpr_y )
   print*,'y=',lpr_y
   ! ...
   
   ! ... free object
   call free_sparse_matrix ( lo_A )
   ! ...
   
   deallocate(lpr_x)
   deallocate(lpr_y)	
end subroutine test1
! ............................................

! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
