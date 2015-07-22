!---------------------------------------------------------------------------------
!
!	THE FOLLOWING ROUTINES ARE THE SAME AS DEFINED BEFORE BUT WE ALLOW THE USE OF LM DEFINED WITH DEGREE OF FREEDOM
!	
!---------------------------------------------------------------------------------
	subroutine create_SparseMatrix_with_LM_MULT ( this, ai_nR, ai_nC, ai_nel, ai_ndof_1, api_LM_1, ai_nen_1, ai_ndof_2, api_LM_2, &
		ai_nen_2 )
	implicit none
		type(csr_matrix) :: this
		integer :: ai_nC, ai_nR
		integer, dimension(:,:,:), pointer :: api_LM_1, api_LM_2			
		integer :: ai_nel, ai_ndof_1, ai_nen_1, ai_nen_2, ai_ndof_2		
		!local var
		integer :: li_err,li_flag
		integer :: li_nnz
		
		! COUNTING NON ZERO ELEMENTS
		li_nnz = count_non_zero_elts_MULT ( ai_nR, ai_nC, ai_nel, ai_ndof_1, api_LM_1, ai_nen_1, ai_ndof_2, api_LM_2, ai_nen_2 )

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
		
		call init_SparseMatrix_MULT ( this, ai_nel, ai_ndof_1, api_LM_1, ai_nen_1, ai_ndof_2, api_LM_2, ai_nen_2 )
		
		this%opr_a ( : ) =0.0_wp
	end subroutine create_SparseMatrix_with_LM_MULT	
!---------------------------------------------------------------------------------
	subroutine init_SparseMatrix_MULT ( this, ai_nel, ai_ndof_1, api_LM_1, ai_nen_1, ai_ndof_2, api_LM_2, ai_nen_2 )
	! _1 FOR ROWS
	! _2 FOR COLUMNS	
	implicit none
		type(csr_matrix) :: this
		integer, dimension(:,:,:), pointer :: api_LM_1, api_LM_2			
		integer :: ai_nel, ai_ndof_1, ai_nen_1, ai_nen_2, ai_ndof_2
		!local var
		integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_i_1, li_i_2, li_index, li_i	, li_size
		integer, dimension(:,:,:,:,:), pointer :: lpi_Marked	
		integer, dimension(this%oi_nR+1) :: lpi_occ		
		integer :: li_err,li_flag			
		integer, dimension(:,:), pointer :: lpi_columns		
		real(wp), dimension(:), pointer :: lpr_tmp			
				
		allocate ( lpi_Marked ( ai_nel, ai_nen_1, ai_nen_2, ai_ndof_1, ai_ndof_2 ) , stat = li_err )
		if ( li_err .ne. 0 ) li_flag = 10		
			
		lpi_Marked ( : , : , : , : , : ) = 0

		lpi_occ ( : ) = 0
		
	! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED	
		do li_e = 1, ai_nel

			do li_b_1 = 1, ai_nen_1

				do li_i_1 = 1, ai_ndof_1

					li_A_1 = api_LM_1 ( li_i_1, li_b_1, li_e )
					if ( li_A_1 == 0 ) then
						cycle
					end if					
										
					do li_b_2 = 1, ai_nen_2

						do li_i_2 = 1, ai_ndof_2
							
							li_A_2 = api_LM_2 ( li_i_2, li_b_2, li_e )
							if ( li_A_2 == 0 ) then
								cycle
							end if
							
							! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (li_A_1, li_A_2)
							if ( done_couple_MULT ( lpi_Marked, li_e, li_A_1, li_A_2, api_LM_1, api_LM_2, &
											   ai_ndof_1, ai_ndof_2, ai_nen_1, ai_nen_2 ) == 0 ) then
								
								! li_A_1 IS THE ROW NUM, li_A_2 THE COLUMN NUM
								! INITIALIZATION OF THE SPARSE MATRIX
								lpi_occ ( li_A_1 ) = lpi_occ ( li_A_1 ) + 1				
								
								lpi_Marked ( li_e, li_b_1, li_b_2, li_i_1, li_i_2 ) = 1							
								
							end if
													
						end do

					end do
					
				end do

			end do

		end do
				
	! WE CREATE lpi_columns, THAT CONTAINS FOR EACH ROW THE COLUMNS INDICES, NOT NECESSARY IN THE ASCENDING WAY.	
	! lpi_columns ( li_A_1, 0 ) CONTAINS THE INDEX OF THE LAST INSERTED COLUMN		
		allocate ( lpi_columns ( this%oi_nR, 0:MAXVAL(lpi_occ(:)) ) , stat = li_err )
		if ( li_err .ne. 0 ) li_flag = 10	
		
		lpi_Marked ( : , : , : , : , : ) = 0
		lpi_columns ( :, : ) = 0
		
	! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED	
		do li_e = 1, ai_nel

			do li_b_1 = 1, ai_nen_1

				do li_i_1 = 1, ai_ndof_1

					li_A_1 = api_LM_1 ( li_i_1, li_b_1, li_e )
					if ( li_A_1 == 0 ) then
						cycle
					end if					
										
					do li_b_2 = 1, ai_nen_2

						do li_i_2 = 1, ai_ndof_2
							
							li_A_2 = api_LM_2 ( li_i_2, li_b_2, li_e )
							if ( li_A_2 == 0 ) then
								cycle
							end if
							
							! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (li_A_1, li_A_2)
							if ( done_couple_MULT ( lpi_Marked, li_e, li_A_1, li_A_2, api_LM_1, api_LM_2, &
											   ai_ndof_1, ai_ndof_2, ai_nen_1, ai_nen_2 ) == 0 ) then
								
								! li_A_1 IS THE ROW NUM, li_A_2 THE COLUMN NUM
								! INITIALIZATION OF THE SPARSE MATRIX	
								lpi_columns ( li_A_1, 0 ) = lpi_columns ( li_A_1, 0 ) + 1
								lpi_columns ( li_A_1, lpi_columns ( li_A_1, 0 ) ) = li_A_2	
								
								lpi_Marked ( li_e, li_b_1, li_b_2, li_i_1, li_i_2 ) = 1							
								
							end if
													
						end do

					end do
					
				end do

			end do

		end do

	! INITIALIZING ia 
		this%opi_ia ( 1 ) = 1
		
		do li_i = 1, this%oi_nR 
		
			this%opi_ia ( li_i + 1 ) = this%opi_ia ( 1 ) + SUM ( lpi_occ ( 1 : li_i ) )
			
		end do

	! INITIALIZING ja		
		do li_e = 1, ai_nel
		
			do li_b_1 = 1, ai_nen_1

				do li_i_1 = 1, ai_ndof_1

					li_A_1 = api_LM_1 ( li_i_1, li_b_1, li_e )
					if ( ( li_A_1 == 0 ) .OR. ( lpi_columns ( li_A_1, 0 ) == 0 ) ) then
						cycle
					end if		
								
					li_size = lpi_columns ( li_A_1, 0 )
					
					allocate ( lpr_tmp ( li_size ) , stat = li_err )
					if ( li_err .ne. 0 ) li_flag = 10								
					
					lpr_tmp ( 1 : li_size ) = real( lpi_columns ( li_A_1, 1 : li_size ) )
					
					call QsortC ( lpr_tmp )
						
					do li_i = 1, li_size

						this%opi_ja ( this%opi_ia ( li_A_1 ) + li_i - 1 ) = int ( lpr_tmp ( li_i ) )																	
						
					end do
					
					lpi_columns ( li_A_1, 0 ) = 0					
					deallocate ( lpr_tmp )
					
				end do

			end do

		end do
						
		deallocate ( lpi_Marked ) 		
		deallocate ( lpi_columns ) 		

	end subroutine init_SparseMatrix_MULT	
!---------------------------------------------------------------------------------
	integer function count_non_zero_elts_MULT ( ai_nR, ai_nC, ai_nel, ai_ndof_1, api_LM_1, ai_nen_1, ai_ndof_2, api_LM_2, ai_nen_2 )
	! _1 FOR ROWS
	! _2 FOR COLUMNS	
	implicit none
		integer :: ai_nR, ai_nC
		integer, dimension(:,:,:), pointer :: api_LM_1, api_LM_2			
		integer :: ai_nel, ai_ndof_1, ai_nen_1, ai_nen_2, ai_ndof_2
		!local var
		integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_i_1, li_i_2, li_index, li_i	, li_size
		integer, dimension(:,:,:,:,:), pointer :: lpi_Marked	
		integer, dimension(ai_nR+1) :: lpi_occ		
		integer :: li_err,li_flag			
		integer, dimension(:,:), pointer :: lpi_columns		
		real(wp), dimension(:), pointer :: lpr_tmp	
		integer :: li_result		
				
		allocate ( lpi_Marked ( ai_nel, ai_nen_1, ai_nen_2, ai_ndof_1, ai_ndof_2 ) , stat = li_err )
		if ( li_err .ne. 0 ) li_flag = 10		
			
		lpi_Marked ( : , : , : , : , : ) = 0

		lpi_occ ( : ) = 0
		
	! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED	
		do li_e = 1, ai_nel

			do li_b_1 = 1, ai_nen_1

				do li_i_1 = 1, ai_ndof_1

					li_A_1 = api_LM_1 ( li_i_1, li_b_1, li_e )
					if ( li_A_1 == 0 ) then
						cycle
					end if					
										
					do li_b_2 = 1, ai_nen_2

						do li_i_2 = 1, ai_ndof_2
							
							li_A_2 = api_LM_2 ( li_i_2, li_b_2, li_e )
							if ( li_A_2 == 0 ) then
								cycle
							end if
							
							! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (li_A_1, li_A_2)
							if ( done_couple_MULT ( lpi_Marked, li_e, li_A_1, li_A_2, api_LM_1, api_LM_2, &
											   ai_ndof_1, ai_ndof_2, ai_nen_1, ai_nen_2 ) == 0 ) then
								
								! li_A_1 IS THE ROW NUM, li_A_2 THE COLUMN NUM
								! INITIALIZATION OF THE SPARSE MATRIX
								lpi_occ ( li_A_1 ) = lpi_occ ( li_A_1 ) + 1				
								
								lpi_Marked ( li_e, li_b_1, li_b_2, li_i_1, li_i_2 ) = 1							
								
							end if
													
						end do

					end do
					
				end do

			end do

		end do
				
	! WE CREATE lpi_columns, THAT CONTAINS FOR EACH ROW THE COLUMNS INDICES, NOT NECESSARY IN THE ASCENDING WAY.	
	! lpi_columns ( li_A_1, 0 ) CONTAINS THE INDEX OF THE LAST INSERTED COLUMN		
		allocate ( lpi_columns ( ai_nR, 0:MAXVAL(lpi_occ(:)) ) , stat = li_err )
		if ( li_err .ne. 0 ) li_flag = 10	
		
		lpi_Marked ( : , : , : , : , : ) = 0
		lpi_columns ( :, : ) = 0
		
	! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED	
		do li_e = 1, ai_nel

			do li_b_1 = 1, ai_nen_1

				do li_i_1 = 1, ai_ndof_1

					li_A_1 = api_LM_1 ( li_i_1, li_b_1, li_e )
					if ( li_A_1 == 0 ) then
						cycle
					end if					
										
					do li_b_2 = 1, ai_nen_2

						do li_i_2 = 1, ai_ndof_2
							
							li_A_2 = api_LM_2 ( li_i_2, li_b_2, li_e )
							if ( li_A_2 == 0 ) then
								cycle
							end if
							
							! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (li_A_1, li_A_2)
							if ( done_couple_MULT ( lpi_Marked, li_e, li_A_1, li_A_2, api_LM_1, api_LM_2, &
											   ai_ndof_1, ai_ndof_2, ai_nen_1, ai_nen_2 ) == 0 ) then
								
								! li_A_1 IS THE ROW NUM, li_A_2 THE COLUMN NUM
								! INITIALIZATION OF THE SPARSE MATRIX	
								lpi_columns ( li_A_1, 0 ) = lpi_columns ( li_A_1, 0 ) + 1
								lpi_columns ( li_A_1, lpi_columns ( li_A_1, 0 ) ) = li_A_2	
								
								lpi_Marked ( li_e, li_b_1, li_b_2, li_i_1, li_i_2 ) = 1							
								
							end if
													
						end do

					end do
					
				end do

			end do

		end do

	! COUNT NON ZERO ELEMENTS		
		li_result = SUM ( lpi_occ ( 1 : ai_nR ) )
						
		deallocate ( lpi_Marked ) 		
		deallocate ( lpi_columns ) 		
		
		count_non_zero_elts_MULT = li_result
	end function count_non_zero_elts_MULT	
!---------------------------------------------------------------------------------	
	integer function done_couple_MULT ( api_Marked, ai_e, ai_A_1, ai_A_2, api_LM_1, api_LM_2, ai_ndof_1, ai_ndof_2, ai_nen_1, &
		ai_nen_2 )
	implicit none
		integer, dimension(:,:,:,:,:) :: api_Marked
		integer, dimension(:,:,:) :: api_LM_1, api_LM_2		
		integer  :: ai_e, ai_A_1, ai_A_2, ai_ndof_1, ai_ndof_2, ai_nen_1, ai_nen_2
		! LOCAL VARIABLES
		integer  :: li_e, li_i_1, li_i_2, li_A_1, li_b_1, li_A_2, li_b_2, li_sum 
		
		li_sum = 0
		 
		do li_e = 1, ai_e
		
			do li_b_2 = 1, ai_nen_2
				
				do li_i_2 = 1, ai_ndof_2
				
					li_A_2 = api_LM_2 ( li_i_2, li_b_2, li_e )

					if ( li_A_2 == ai_A_2 ) then

						do li_b_1 = 1, ai_nen_1
							
							do li_i_1 = 1, ai_ndof_1
							
								li_A_1 = api_LM_1 ( li_i_1, li_b_1, li_e )	
								
								if ( li_A_1 == ai_A_1 ) then																										
									li_sum = li_sum + api_Marked ( li_e, li_b_1, li_b_2, li_i_1, li_i_2 )
								end if
								
							end do
						
						end do

					end if						
				
				end do
			
			end do
					
		end do
		
		if ( li_sum > 0 ) then
			done_couple_MULT = 1
		else
			done_couple_MULT = 0						
		end if
			
		return
		
	end function done_couple_MULT