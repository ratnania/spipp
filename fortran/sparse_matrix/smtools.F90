subroutine remove_zeros ( ao_A, al_resize, ao_M )
implicit none
	type(CSR_MATRIX)  :: ao_A
	logical, optional :: al_resize
	type(CSR_MATRIX), optional  :: ao_M	
	! LOCAL VARIABLES
	integer  :: li_nnz
	integer  :: li_i
	integer  :: li_k
	integer  :: li_j
	integer  :: li_nR	
	integer  :: li_nC
	integer  :: li_nel		
	logical  :: ll_resize	
	logical  :: ll_use_loc
	type(CSR_MATRIX)  :: lo_M

	ll_resize = .FALSE.
	if ( present ( al_resize ) ) then
		ll_resize = al_resize
	end if
	
	ll_use_loc = .TRUE.
	if ( present ( ao_M ) ) then
		ll_use_loc = .FALSE.
	end if
	
	li_nR	= ao_A%oi_nR			
	li_nC	= ao_A%oi_nC			
	li_nel	= ao_A%oi_nel
	
	if ( ll_use_loc ) then
	
		call create_SparseMatrix ( lo_M, li_nR, li_nC, li_nel )

		!... parcours pour copier les elts non nuls dans le nouveau tableau
		!... nbr d elt non nuls
		li_nnz   = 0		
		do li_i = 1 , li_nR
		
			lo_M%opi_ia(li_i) = li_nnz + 1
			
			do li_k = ao_A%opi_ia ( li_i ) , ao_A%opi_ia ( li_i + 1 ) - 1
			
				li_j = ao_A%opi_ja ( li_k )
				
				if ( dabs ( ao_A%opr_a ( li_k ) ) >=  epsilon ) then
				
					li_nnz = li_nnz + 1
					
					lo_M%opi_ja ( li_nnz ) = li_j
					
					lo_M%opr_a ( li_nnz )  =  ao_A%opr_a ( li_k )
					
				end if			
				
			end do
			
		end do	
	
		lo_M%opi_ia ( li_nR + 1 ) = li_nnz + 1
		
	else
			
		!... parcours pour copier les elts non nuls dans le nouveau tableau
		!... nbr d elt non nuls
		li_nnz   = 0		
		do li_i = 1 , li_nR
		
			ao_M%opi_ia(li_i) = li_nnz + 1
			
			do li_k = ao_A%opi_ia ( li_i ) , ao_A%opi_ia ( li_i + 1 ) - 1
			
				li_j = ao_A%opi_ja ( li_k )
				
				if ( dabs ( ao_A%opr_a ( li_k ) ) >=  epsilon ) then
				
					li_nnz = li_nnz + 1
					
					ao_M%opi_ja ( li_nnz ) = li_j
					
					ao_M%opr_a ( li_nnz )  =  ao_A%opr_a ( li_k )
					
				end if			
				
			end do
			
		end do	
	
		ao_M%opi_ia ( li_nR + 1 ) = li_nnz + 1
					
	end if					

	if ( ll_resize ) then
		!... free the old matrix
		call free_SparseMatrix ( ao_A )
		
		!... create a new matrix with the good sizes
		call create_SparseMatrix ( ao_A, li_nR, li_nC, li_nnz )	
	end if

	if ( ll_use_loc ) then
		!... copy data from the temp matrix
		ao_A%oi_nel = li_nnz
		ao_A%opi_ia ( 1 : li_nR + 1 ) = lo_M%opi_ia ( 1 : li_nR + 1 )
		ao_A%opi_ja ( 1 : li_nnz    ) = lo_M%opi_ja ( 1 : li_nnz )
		ao_A%opr_a  ( 1 : li_nnz    ) = lo_M%opr_a  ( 1 : li_nnz )	

		!... free the temp matrix	
		call free_SparseMatrix ( lo_M )	
		
	else
		!... copy data from the temp matrix
		ao_A%oi_nel = li_nnz
		ao_A%opi_ia ( 1 : li_nR + 1 ) = ao_M%opi_ia ( 1 : li_nR + 1 )
		ao_A%opi_ja ( 1 : li_nnz    ) = ao_M%opi_ja ( 1 : li_nnz )
		ao_A%opr_a  ( 1 : li_nnz    ) = ao_M%opr_a  ( 1 : li_nnz )	
	
	end if

end subroutine remove_zeros
!---------------------------------------------------------------------------------
subroutine set_direct_access_i ( ao_A )
implicit none
	type(CSR_MATRIX)  :: ao_A
	! LOCAL VARIABLES
	integer  :: li_i
	integer  :: li_k
	integer  :: li_nR		
	
	ao_A%ol_use_mm_format = .TRUE.
	
	ALLOCATE ( ao_A%opi_i ( ao_A%oi_nel ) )

	li_nR	= ao_A%oi_nR
		
	do li_i = 1 , li_nR
		
		do li_k = ao_A%opi_ia ( li_i ) , ao_A%opi_ia ( li_i + 1 ) - 1
				
			ao_A%opi_i  ( li_k ) = li_i		
			
		end do
		
	end do				
	
end subroutine set_direct_access_i