!---------------------------------------------------------------------------------
	subroutine create_SparseMatrix_with_LM_filter   ( this, ai_nR, ai_nC, ai_nel	&
													, api_LM_1, ai_nen_1			&
												    , api_LM_2, ai_nen_2			&
													, ai_nfilter , api_filter )
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
		!> param [in] ai_nfilter :: FILTER DIMENSION
		integer :: ai_nfilter	
		!> param[in] api_filter : IS AN ARRAY CONTAINING A SET OF ELEMENTS THAT WE WILL LOOP OVER
		integer, dimension ( : ) :: api_filter			
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

		li_nnz = count_non_zero_elts_filter ( ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1	&
											, api_LM_2, ai_nen_2, lpi_columns, lpi_occ	&
											, ai_nfilter, api_filter )
	
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

		call init_SparseMatrix_filter  ( this, ai_nel, api_LM_1, ai_nen_1			&
										, api_LM_2, ai_nen_2, lpi_columns, lpi_occ	&
										, ai_nfilter, api_filter )					

		this%opr_a ( : ) = 0.0_wp
		
		deallocate(lpi_columns)		
		deallocate(lpi_occ)		
	end subroutine create_SparseMatrix_with_LM_filter	
!---------------------------------------------------------------------------------
	subroutine init_SparseMatrix_filter ( this, ai_nel, api_LM_1, ai_nen_1			&
										, api_LM_2, ai_nen_2, api_columns, api_occ	&
										, ai_nfilter, api_filter )
	! _1 FOR ROWS
	! _2 FOR COLUMNS	
	implicit none
		type(csr_matrix) :: this
		integer, dimension(:,:), pointer :: api_LM_1, api_LM_2			
		integer :: ai_nel, ai_nen_1, ai_nen_2
		integer, dimension(:,:), pointer :: api_columns	
		integer, dimension(:), pointer :: api_occ
		integer :: ai_nfilter
		integer, dimension ( : ) :: api_filter				
		!local var
		integer :: li_e
		integer :: li_b_1
		integer :: li_A_1
		integer :: li_b_2
		integer :: li_A_2		
		integer :: li_index
		integer :: li_i
		integer :: li_size		
		integer :: li_err
		integer :: li_flag				
		real(wp), dimension(:), pointer :: lpr_tmp							

		allocate ( lpr_tmp ( 1000 ) , stat = li_err )
		if ( li_err .ne. 0 ) li_flag = 10								

	! INITIALIZING ia 
		this%opi_ia ( 1 ) = 1
		
		do li_i = 1, this%oi_nR 
		
			this%opi_ia ( li_i + 1 ) = this%opi_ia ( 1 ) + SUM ( api_occ ( 1 : li_i ) )
			
		end do

	! INITIALIZING ja		
		do li_index = 1, ai_nfilter
		
			li_e = api_filter ( li_index )
		
			do li_b_1 = 1, ai_nen_1

				li_A_1 = api_LM_1 ( li_b_1, li_e )
			
				if ( li_A_1 == 0 ) then
					cycle
				end if		
							
				if ( api_columns ( li_A_1, 0 ) == 0 ) then
					cycle
				end if		
											
				li_size = api_columns ( li_A_1, 0 )
				
				lpr_tmp ( 1 : li_size ) = real( api_columns ( li_A_1, 1 : li_size ) )
				
				call QsortC ( lpr_tmp ( 1 : li_size ) )
					
				do li_i = 1, li_size

					this%opi_ja ( this%opi_ia ( li_A_1 ) + li_i - 1 ) = int ( lpr_tmp ( li_i ) )																	
					
				end do
				
				api_columns ( li_A_1, 0 ) = 0					

			end do

		end do	
		
		deallocate ( lpr_tmp )		

	end subroutine init_SparseMatrix_filter	
!---------------------------------------------------------------------------------
	integer function count_non_zero_elts_filter ( ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1	&
												, api_LM_2, ai_nen_2, api_columns, api_occ	&
												, ai_nfilter, api_filter )
	! _1 FOR ROWS
	! _2 FOR COLUMNS	
	implicit none
		integer :: ai_nR, ai_nC
		integer, dimension(:,:), pointer :: api_LM_1, api_LM_2			
		integer :: ai_nel, ai_nen_1, ai_nen_2
		integer, dimension(:,:), pointer :: api_columns	
		integer, dimension(:), pointer :: api_occ	
		integer  :: ai_nfilter
		integer, dimension ( : ) :: api_filter			
		!local var
		integer :: li_e
		integer :: li_b_1
		integer :: li_A_1
		integer :: li_b_2
		integer :: li_A_2		
		integer :: li_index
		integer :: li_i
		integer :: li_size		
		integer :: li_err
		integer :: li_flag	
		real(wp), dimension(:), pointer :: lpr_tmp	
		integer :: li_result	
		logical :: ll_done
				
	! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED	
		do li_index = 1, ai_nfilter
		
			li_e = api_filter ( li_index )

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
		
		count_non_zero_elts_filter = li_result
	end function count_non_zero_elts_filter