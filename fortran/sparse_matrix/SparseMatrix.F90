!
!	SparseMatrix.f90
!	SparseMatrix
!
!	Created by ahmed ratnani on 02/01/10.
!	Copyright 2010 __MyCompanyName__. All rights reserved.
!
module SparseMatrix_Module
use used_precision
use ISOBOX_Def_Module
use qsort_c_module
implicit none
	
contains
    include "smout.f90"
    include "SparseMatrix2.f90"
    include "SparseMatrix3.f90"    
    include "SparseMatrix4.f90"	
    include "smcg.f90"
    include "smoperations.f90"    
    include "smtools.f90"	
!---------------------------------------------------------------------------------
	subroutine create_SparseMatrix ( this, ai_nR, ai_nC, ai_nel )
	implicit none
		!> param[inout] this : CSR MATRIX STRUCTURE
		type(csr_matrix) :: this		
		!> param[in] ai_nC : NUMBER OF COLUMNS
		integer :: ai_nC
		!> param[in] ai_nR : NUMBER OF ROWS
		integer :: ai_nR
		!> param[in] ai_nel : NUMBER OF NON ZERO ELEMENTS		
		integer :: ai_nel				
		!local var
		integer :: li_err,li_flag
		
		this%ol_use_mm_format = .FALSE.
		
		this%oi_nR   = ai_nR
		this%oi_nC   = ai_nC		
		this%oi_nel  = ai_nel
		
		allocate(this%opi_ia(this%oi_nR+1),stat=li_err)
		if (li_err.ne.0) li_flag=10	

		allocate(this%opi_ja(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=20	
				
		allocate(this%opr_a(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=30			
		
		this%opr_a ( : ) = 0.0_wp
	end subroutine create_SparseMatrix
!---------------------------------------------------------------------------------
! THIS ROUTINE CREATES A NEW FULL MATRIX FROM A COMPRESSED SYMMETRIC MATRIX
	subroutine create_SparseMatrix_from_sym ( this, ao_A )
	implicit none
		!> param[inout] this : CSR MATRIX STRUCTURE	
		type(csr_matrix) :: this
		!> param[in] ao_A : CSR MATRIX STRUCTURE		
		type(csr_matrix) :: ao_A		
		!local var
		integer :: li_err,li_flag
		integer :: li_i, li_k, li_j, li_count, li_i_tmp, li_j_tmp, li_index
		integer, dimension(ao_A%oi_nR+1) :: lpi_occ
					
		lpi_occ ( : ) = 0
		
		! COMPUTING THE NUMBER OF NON ZERO DIAGONAL ELEMENT 
		li_count = 0
		do li_i = 1, ao_A%oi_nR 
			if ( ao_A%opi_ja ( ao_A%opi_ia ( li_i ) ) == li_i ) then
				li_count = li_count + 1
			end if
		end do	
		
		this%ol_use_mm_format = .FALSE.		
		
		this%oi_nel = 2 * ao_A%oi_nel - li_count
		this%oi_nR = ao_A%oi_nR
		this%oi_nC = ao_A%oi_nC	
	
		allocate(this%opi_ia(this%oi_nR+1),stat=li_err)
		if (li_err.ne.0) li_flag=10	

		allocate(this%opi_ja(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=20	
				
		allocate(this%opr_a(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=30	
		
		! COPY THE OLD MATRIX IN THE NEW MATRIX		
		do li_i = 1, ao_A%oi_nR
			lpi_occ ( li_i ) = ao_A%opi_ia ( li_i + 1 ) - ao_A%opi_ia ( li_i )
			do li_i_tmp = 1, li_i - 1
				li_k = ao_A%opi_ia ( li_i_tmp )
				do while ( ( li_k <= ao_A%opi_ia ( li_i_tmp + 1 ) - 1 ) .AND. ( ao_A%opi_ja ( li_k ) <= li_i ) )			
					if ( li_i == ao_A%opi_ja ( li_k ) ) then
						lpi_occ ( li_i ) = lpi_occ ( li_i ) + 1
					end if
					li_k = li_k + 1
				end do
			end do
		end do

		this%opi_ia ( 1 ) = 1
		
		do li_i = 1, this%oi_nR 
		
			this%opi_ia ( li_i + 1 ) = this%opi_ia ( 1 ) + SUM ( lpi_occ ( 1 : li_i ) )
			
		end do

		li_index = 1
		do li_i = 1, ao_A%oi_nR
			do li_i_tmp = 1, li_i - 1
				li_k = ao_A%opi_ia ( li_i_tmp )
				do while ( ( li_k <= ao_A%opi_ia ( li_i_tmp + 1 ) - 1 ) .AND. ( ao_A%opi_ja ( li_k ) <= li_i ) )			
					if ( li_i == ao_A%opi_ja ( li_k ) ) then
						this%opi_ja ( li_index ) = li_i_tmp
						li_index = li_index + 1
					end if
					li_k = li_k + 1
				end do
			end do
			
			this%opi_ja ( li_index : li_index + ao_A%opi_ia ( li_i + 1 ) - 1 - ao_A%opi_ia ( li_i ) ) = &
						ao_A%opi_ja ( ao_A%opi_ia ( li_i ) : ao_A%opi_ia ( li_i + 1 ) - 1 )	
			
			li_index = li_index + ao_A%opi_ia ( li_i + 1 ) - ao_A%opi_ia ( li_i )
		end do
				
	end subroutine create_SparseMatrix_from_sym	
!---------------------------------------------------------------------------------
	subroutine create_SparseMatrix_two_fields ( this, ai_nel, ao_field1, ao_field2 )
	implicit none
		!> param[inout] this : CSR MATRIX STRUCTURE
		type(csr_matrix) :: this	
		integer :: ai_nel	
		type(FIELD) :: ao_field1
		type(FIELD) :: ao_field2		

		call create_SparseMatrix_with_LM    ( this					&
											, ao_field1%oi_sizePB	&
											, ao_field2%oi_sizePB	&
											, ai_nel				&
											, ao_field1%opi_LM		&
											, ao_field1%oi_nen		& 
											, ao_field2%opi_LM		&
											, ao_field2%oi_nen		)

	end subroutine create_SparseMatrix_two_fields	
!---------------------------------------------------------------------------------
	subroutine create_SparseMatrix_two_fields_filter ( this, ai_nel, ao_field1, ao_field2, ai_nfilter, api_filter )
	implicit none
		!> param[inout] this : CSR MATRIX STRUCTURE
		type(csr_matrix) :: this	
		integer :: ai_nel	
		type(FIELD) :: ao_field1
		type(FIELD) :: ao_field2	
		integer :: ai_nfilter	
		integer, dimension ( : ) :: api_filter						

		call create_SparseMatrix_with_LM_filter ( this					&
												, ao_field1%oi_sizePB	&
												, ao_field2%oi_sizePB	&
												, ai_nel				&
												, ao_field1%opi_LM		&
												, ao_field1%oi_nen		& 
												, ao_field2%opi_LM		&
												, ao_field2%oi_nen		&
												, ai_nfilter			&
												, api_filter			)

	end subroutine create_SparseMatrix_two_fields_filter		
!---------------------------------------------------------------------------------
	subroutine create_SparseMatrix_with_LM ( this, ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2 )
	implicit none
		!> param[inout] this : CSR MATRIX STRUCTURE		
		type(csr_matrix) :: this
		!> param[in] ai_nC : NUMBER OF COLUMNS, IT IS THE DIMENSION OF THE 1st SPACE		
		integer :: ai_nC
		!> param[in] ai_nR : NUMBER OF ROWS, IT IS THE DIMENSION OF THE 2nd SPACE		
		integer :: ai_nR
		!> param[in] ai_nel : NUMBER OF ELEMENTS (WITH NON ZERO MEASURE) IN THE PATCH   		
		integer :: ai_nel		
		!> param[in] api_LM_1 : LM ARRAY FOR ROWS
		integer, dimension(:,:), pointer :: api_LM_1
		!> param[in] api_LM_1 : LM ARRAY FOR COLUMNS
		integer, dimension(:,:), pointer :: api_LM_2					
		!> param[in] ai_nen_1 : NUMBER OF NON VANISHING FUNCTIONS PER ELEMENT, IN THE 1st SPACE 
		integer :: ai_nen_1
		!> param[in] ai_nen_2 : NUMBER OF NON VANISHING FUNCTIONS PER ELEMENT, IN THE 2nd SPACE
		integer :: ai_nen_2					
		!local var
		integer :: li_err,li_flag
		integer :: li_nnz
		integer, dimension(:,:), pointer :: lpi_columns	
		integer, dimension(:), pointer :: lpi_occ			

		allocate(lpi_columns(ai_nR, 0:8*ai_nen_2),stat=li_err)
		if (li_err.ne.0) li_flag=1	
		
		allocate(lpi_occ(ai_nR+1),stat=li_err)
		if (li_err.ne.0) li_flag=2	
								
		lpi_columns ( :, : ) = 0
		lpi_occ ( : ) = 0
				
		! COUNTING NON ZERO ELEMENTS

		li_nnz = count_non_zero_elts ( ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, lpi_columns, lpi_occ )
	
		this%ol_use_mm_format = .FALSE.
			
		this%oi_nR   = ai_nR
		this%oi_nC   = ai_nC		
		this%oi_nel  = li_nnz
		
		allocate(this%opi_ia(this%oi_nR+1),stat=li_err)
		if (li_err.ne.0) li_flag=10	

		allocate(this%opi_ja(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=20	
				
		allocate(this%opr_a(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=30		

		call init_SparseMatrix ( this, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, lpi_columns, lpi_occ )					

		this%opr_a ( : ) = 0.0_wp
		
		deallocate(lpi_columns)		
		deallocate(lpi_occ)		
	end subroutine create_SparseMatrix_with_LM	
!---------------------------------------------------------------------------------
	subroutine init_SparseMatrix ( this, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, api_columns, api_occ )
	! _1 FOR ROWS
	! _2 FOR COLUMNS	
	implicit none
		type(csr_matrix) :: this
		integer, dimension(:,:), pointer :: api_LM_1, api_LM_2			
		integer :: ai_nel, ai_nen_1, ai_nen_2
		integer, dimension(:,:), pointer :: api_columns	
		integer, dimension(:), pointer :: api_occ	
		!local var
		integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_index, li_i	, li_size		
		integer :: li_err,li_flag				
		real(wp), dimension(:), pointer :: lpr_tmp							

	! INITIALIZING ia 
		this%opi_ia ( 1 ) = 1
		
		do li_i = 1, this%oi_nR 
		
			this%opi_ia ( li_i + 1 ) = this%opi_ia ( 1 ) + SUM ( api_occ ( 1 : li_i ) )
			
		end do

	! INITIALIZING ja		
		do li_e = 1, ai_nel
		
			do li_b_1 = 1, ai_nen_1

				li_A_1 = api_LM_1 ( li_b_1, li_e )
			
				if ( li_A_1 == 0 ) then
					cycle
				end if		
							
				if ( api_columns ( li_A_1, 0 ) == 0 ) then
					cycle
				end if		
											
				li_size = api_columns ( li_A_1, 0 )

				allocate ( lpr_tmp ( li_size ) , stat = li_err )
				if ( li_err .ne. 0 ) li_flag = 10								
				
				lpr_tmp ( 1 : li_size ) = real( api_columns ( li_A_1, 1 : li_size ) )
				
				call QsortC ( lpr_tmp )
					
				do li_i = 1, li_size

					this%opi_ja ( this%opi_ia ( li_A_1 ) + li_i - 1 ) = int ( lpr_tmp ( li_i ) )																	
					
				end do
				
				api_columns ( li_A_1, 0 ) = 0					
				deallocate ( lpr_tmp )

			end do

		end do	

	end subroutine init_SparseMatrix	
!---------------------------------------------------------------------------------
	integer function count_non_zero_elts ( ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, api_columns, api_occ )
	! _1 FOR ROWS
	! _2 FOR COLUMNS	
	implicit none
		integer :: ai_nR, ai_nC
		integer, dimension(:,:), pointer :: api_LM_1, api_LM_2			
		integer :: ai_nel, ai_nen_1, ai_nen_2
		integer, dimension(:,:), pointer :: api_columns	
		integer, dimension(:), pointer :: api_occ	
		!local var
		integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_i			
		integer :: li_err,li_flag			
		real(wp), dimension(:), pointer :: lpr_tmp	
		integer :: li_result	
		logical :: ll_done
				
	! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED	
		do li_e = 1, ai_nel

			do li_b_1 = 1, ai_nen_1

				li_A_1 = api_LM_1 ( li_b_1, li_e )
				if ( li_A_1 == 0 ) then
					cycle
				end if					
									
				do li_b_2 = 1, ai_nen_2
						
					li_A_2 = api_LM_2 ( li_b_2, li_e )
					if ( li_A_2 == 0 ) then
						cycle
					end if
					
					ll_done = .false.
					! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (li_A_1, li_A_2)
					do li_i = 1 , api_columns ( li_A_1, 0 ) 

						if ( api_columns ( li_A_1, li_i ) /= li_A_2 ) then
							cycle
						end if	

						ll_done = .true.	
						exit				

					end do

					if ( .not. ll_done ) then
						
						api_occ ( li_A_1 ) = api_occ ( li_A_1 ) + 1	
												
						! li_A_1 IS THE ROW NUM, li_A_2 THE COLUMN NUM
						! INITIALIZATION OF THE SPARSE MATRIX	
						api_columns ( li_A_1, 0 ) = api_columns ( li_A_1, 0 ) + 1
						api_columns ( li_A_1, api_columns ( li_A_1, 0 ) ) = li_A_2							
						
					end if

				end do					

			end do

		end do

	! COUNT NON ZERO ELEMENTS		
		li_result = SUM ( api_occ ( 1 : ai_nR ) )
		
		count_non_zero_elts = li_result
	end function count_non_zero_elts	
!---------------------------------------------------------------------------------
	subroutine free_SparseMatrix(this)
	implicit none
		type(csr_matrix) :: this
		
		deallocate(this%opi_ia)
		deallocate(this%opi_ja)
		deallocate(this%opr_a)

		if ( this%ol_use_mm_format ) then
		
			DEALLOCATE ( this%opi_i )
			
		end if
		
	end subroutine free_SparseMatrix	
			
end module SparseMatrix_Module
