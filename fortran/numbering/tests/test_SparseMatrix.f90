!
!	test_SparseMatrix.f90
!	SparseMatrixStr
!
!	Created by ahmed ratnani on 10/04/09.
!	Copyright 2009 __MyCompanyName__. All rights reserved.
!

program test_SparseMatrix
use used_precision
use SparseMatrix_Module
use io_module
implicit none

	call test13 ()
	
contains
	include "tools.f90"
!---------------------------------------------------------------------------------
	subroutine test1 ()
	implicit none
		type(csr_matrix) :: lo_A
		integer  :: li_n, li_nel
		real(wp),dimension(:),pointer :: lpr_x,lpr_y
		integer  :: li_err,li_flag,li_i,li_j
		
		li_n = 8
		li_nel = 18
		
		allocate(lpr_x(li_n),stat=li_err)
			if (li_err.ne.0) li_flag=10	
			
		allocate(lpr_y(li_n),stat=li_err)
			if (li_err.ne.0) li_flag=10		
			
		call create_SparseMatrix ( lo_A, li_n, li_n, li_nel )
		
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
		
			lpr_x(li_i) = 1.0_wp*li_i
			
		enddo
		
		lpr_y(1:li_n) = 0.0_wp	
		
		call printProfil_csrMatrix ( lo_A )
		
		call Mult_MV ( lo_A, lpr_x, lpr_y )
		
		print*,'y=',lpr_y
		
		call free_SparseMatrix ( lo_A )
	
		deallocate(lpr_x)
		deallocate(lpr_y)	
	end subroutine test1
!---------------------------------------------------------------------------------	
	subroutine test2 ()
	implicit none
		type(csr_matrix) :: lo_A
		integer  :: li_n, li_nel
		real(wp),dimension(:),pointer :: lpr_x,lpr_y
		integer  :: li_err,li_flag,li_i,li_j
		
		li_n = 8
		li_nel = 22
		
		allocate(lpr_x(li_n),stat=li_err)
			if (li_err.ne.0) li_flag=10	
			
		allocate(lpr_y(li_n),stat=li_err)
			if (li_err.ne.0) li_flag=10				
			
		call create_SparseMatrix ( lo_A, li_n, li_n, li_nel )
		
		lo_A%opi_ia(1) =	1
		lo_A%opi_ia(2) =	3	
		lo_A%opi_ia(3) =	6	
		lo_A%opi_ia(4) =	9	
		lo_A%opi_ia(5) =	12
		lo_A%opi_ia(6) =	15	
		lo_A%opi_ia(7) =	18
		lo_A%opi_ia(8) =	21
		lo_A%opi_ia(9) =	23
		
		lo_A%opi_ja(1) =	1
		lo_A%opi_ja(2) =	2
		
		lo_A%opi_ja(3) =	1	
		lo_A%opi_ja(4) =	2	
		lo_A%opi_ja(5) =	3
		
		lo_A%opi_ja(6) =	2	
		lo_A%opi_ja(7) =	3
		lo_A%opi_ja(8) =	4
		
		lo_A%opi_ja(9) =	3
		lo_A%opi_ja(10) =	4	
		lo_A%opi_ja(11) =	5		
				
		lo_A%opi_ja(12) =	4		
		lo_A%opi_ja(13) =	5	
		lo_A%opi_ja(14) =	6				
		
		lo_A%opi_ja(15) =	5
		lo_A%opi_ja(16) =	6
		lo_A%opi_ja(17) =	7		
		
		lo_A%opi_ja(18) =	6	
		lo_A%opi_ja(19) =	7
		lo_A%opi_ja(20) =	8
		
		lo_A%opi_ja(21) =	7
		lo_A%opi_ja(22) =	8
						
		lo_A%opr_a(1) =		2.0_wp
		lo_A%opr_a(2) =		1.0_wp
		lo_A%opr_a(3) =		1.0_wp	
		lo_A%opr_a(4) =		2.0_wp	
		lo_A%opr_a(5) =		1.0_wp
		lo_A%opr_a(6) =		1.0_wp	
		lo_A%opr_a(7) =		2.0_wp
		lo_A%opr_a(8) =		1.0_wp
		lo_A%opr_a(9) =		1.0_wp
		lo_A%opr_a(10) =	2.0_wp	
		lo_A%opr_a(11) =	1.0_wp		
		lo_A%opr_a(12) =	1.0_wp		
		lo_A%opr_a(13) =	2.0_wp	
		lo_A%opr_a(14) =	1.0_wp		
		lo_A%opr_a(15) =	1.0_wp	
		lo_A%opr_a(16) =	2.0_wp
		lo_A%opr_a(17) =	1.0_wp
		lo_A%opr_a(18) =	1.0_wp
		lo_A%opr_a(19) =	2.0_wp		
		lo_A%opr_a(20) =	1.0_wp
		lo_A%opr_a(21) =	1.0_wp					
		lo_A%opr_a(22) =	2.0_wp		
				
		lpr_y ( 1 : li_n ) = 1.0_wp	
		
		call printProfil_csrMatrix ( lo_A )		

	PRINT*,'> Solving Linear system with Conjugate Gradiant'
		call Gradient_conj ( lo_A, lpr_y, lpr_x, 600, 0.00000001_wp )

		print*,'x=',lpr_x

	PRINT*,'> Solving Linear system with Pardiso'				
		call solvePardiso2D ( li_n, lo_A , lpr_y, lpr_x )

		print*,'x=',lpr_x
		
		call free_SparseMatrix ( lo_A )
		
		deallocate(lpr_x)
		deallocate(lpr_y)	
	end subroutine test2
!---------------------------------------------------------------------------------	
	subroutine test3 ()
	implicit none
		type(csr_matrix) :: lo_A
		integer  :: li_nR, li_nC, li_nel
		real(wp),dimension(:),pointer :: lpr_x,lpr_y	
		integer  :: li_err,li_flag,li_i,li_j
		
		li_nC = 8
		li_nR = 7		
		li_nel = 22
		
		allocate(lpr_x(li_nC),stat=li_err)
			if (li_err.ne.0) li_flag=10	
			
		allocate(lpr_y(li_nR),stat=li_err)
			if (li_err.ne.0) li_flag=10			
			
		call create_SparseMatrix ( lo_A, li_nR, li_nC, li_nel )
		
		lo_A%opi_ia(1) =	1
		lo_A%opi_ia(2) =	3	
		lo_A%opi_ia(3) =	6	
		lo_A%opi_ia(4) =	9	
		lo_A%opi_ia(5) =	12
		lo_A%opi_ia(6) =	15	
		lo_A%opi_ia(7) =	18
		lo_A%opi_ia(8) =	21
		
		lo_A%opi_ja(1) =	1
		lo_A%opi_ja(2) =	2
		
		lo_A%opi_ja(3) =	1	
		lo_A%opi_ja(4) =	2	
		lo_A%opi_ja(5) =	3
		
		lo_A%opi_ja(6) =	2	
		lo_A%opi_ja(7) =	3
		lo_A%opi_ja(8) =	4
		
		lo_A%opi_ja(9) =	3
		lo_A%opi_ja(10) =	4	
		lo_A%opi_ja(11) =	5		
				
		lo_A%opi_ja(12) =	4		
		lo_A%opi_ja(13) =	5	
		lo_A%opi_ja(14) =	6				
		
		lo_A%opi_ja(15) =	5
		lo_A%opi_ja(16) =	6
		lo_A%opi_ja(17) =	7		
		
		lo_A%opi_ja(18) =	6	
		lo_A%opi_ja(19) =	7
		lo_A%opi_ja(20) =	8
						
		lo_A%opr_a(1) =		2.0_wp
		lo_A%opr_a(2) =		1.0_wp
		lo_A%opr_a(3) =		1.0_wp	
		lo_A%opr_a(4) =		2.0_wp	
		lo_A%opr_a(5) =		1.0_wp
		lo_A%opr_a(6) =		1.0_wp	
		lo_A%opr_a(7) =		2.0_wp
		lo_A%opr_a(8) =		1.0_wp
		lo_A%opr_a(9) =		1.0_wp
		lo_A%opr_a(10) =	2.0_wp	
		lo_A%opr_a(11) =	1.0_wp		
		lo_A%opr_a(12) =	1.0_wp		
		lo_A%opr_a(13) =	2.0_wp	
		lo_A%opr_a(14) =	1.0_wp		
		lo_A%opr_a(15) =	1.0_wp	
		lo_A%opr_a(16) =	2.0_wp
		lo_A%opr_a(17) =	1.0_wp
		lo_A%opr_a(18) =	1.0_wp
		lo_A%opr_a(19) =	2.0_wp		
		lo_A%opr_a(20) =	1.0_wp	
				
		lpr_x ( 1 : li_nC ) = 1.0_wp	
		
		call printProfil_csrMatrix ( lo_A )		
		
		call Mult_MV ( lo_A, lpr_x, lpr_y )
		
		print*,'y=',lpr_y
		
		call free_SparseMatrix ( lo_A )
		
		deallocate(lpr_x)
		deallocate(lpr_y)	
	end subroutine test3	
!---------------------------------------------------------------------------------	
	subroutine test10 ()
	implicit none
		type(csr_matrix) :: lo_A, lo_A_sym
		integer  :: li_nR, li_nC, li_nel
		integer  :: li_err,li_flag,li_i,li_j, li_ndof, li_nen, li_nnz, li_e
		integer, dimension(:,:), pointer :: lpi_LM			
					
		call import_param ( li_nC, li_nR, li_nel, li_nnz, li_nen, li_ndof )
		
		allocate ( lpi_LM ( li_nen, li_nel ) , stat = li_err )
		if ( li_err .ne. 0 ) li_flag = 10	
								
		call create_SparseMatrix ( lo_A_sym, li_nR, li_nC, li_nnz )
		
		call import_data ( lo_A_sym, lpi_LM, li_nel, li_nen, li_ndof )

		call create_SparseMatrix_from_sym ( lo_A, lo_A_sym )

		lo_A%opr_a ( : ) = 1.0_wp
		call printProfil_csrMatrix ( lo_A )			

		call free_SparseMatrix ( lo_A )	
	end subroutine test10	
!---------------------------------------------------------------------------------	
	subroutine test13 ()
	implicit none
		type(csr_matrix) :: lo_Mat
		integer  :: li_nR, li_nC, li_nel
		integer  :: li_err,li_flag,li_i,li_j, li_ndof, li_nen, li_nnz, li_e
		integer, dimension(:,:), pointer :: lpi_LM			

		call import_param ( li_nC, li_nR, li_nel, li_nnz, li_nen, li_ndof )

	! CREATING MATRIX FROM LM ARRAY
		allocate ( lpi_LM ( li_nen, li_nel ) , stat = li_err )
		if ( li_err .ne. 0 ) li_flag = 10	

		call import_LM ( lpi_LM, li_nel, li_nen )
			
		call create_CSR ( lo_Mat, li_nR, li_nC, li_nel, lpi_LM, li_nen, lpi_LM, li_nen )
										
PRINT*,'////////////////////////////////////'	

		lo_Mat%opr_a ( : ) = 1.0_wp	 	

!		call printProfil_csrMatrix ( lo_Mat )	
		
		PRINT*,'lo_Mat%oi_nel =',lo_Mat%oi_nel	
		
		call free_SparseMatrix ( lo_Mat )
		
		deallocate ( lpi_LM )		
		
	end subroutine test13	
!---------------------------------------------------------------------------------	
	subroutine test13_MULT ()
	implicit none
		type(csr_matrix) :: lo_Mat
		integer  :: li_nR, li_nC, li_nel
		integer  :: li_err,li_flag,li_i,li_j, li_ndof, li_nen, li_nnz, li_e
		integer, dimension(:,:,:), pointer :: lpi_LM	
					
		call import_param ( li_nC, li_nR, li_nel, li_nnz, li_nen, li_ndof )

	! CREATING MATRIX FROM LM ARRAY
		allocate ( lpi_LM ( li_ndof, li_nen, li_nel ) , stat = li_err )
		if ( li_err .ne. 0 ) li_flag = 10	

		call import_LM_MULT ( lpi_LM, li_nel, li_nen, li_ndof )
		
		call create_CSR ( lo_Mat, li_nR, li_nC, li_nel, li_ndof, lpi_LM, li_nen, li_ndof, lpi_LM, li_nen )
										
PRINT*,'////////////////////////////////////'	

		lo_Mat%opr_a ( : ) = 1.0_wp	 	

		call printProfil_csrMatrix ( lo_Mat )	
		
		PRINT*,'lo_Mat%oi_nel =',lo_Mat%oi_nel	
		
		call free_SparseMatrix ( lo_Mat )
		
		deallocate ( lpi_LM )
	end subroutine test13_MULT							
end program test_SparseMatrix	