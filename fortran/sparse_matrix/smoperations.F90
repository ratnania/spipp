!---------------------------------------------------------------------------------
	subroutine Mult_CSR_Matrix_Vector(this,apr_x,apr_y)
	implicit none
		type(csr_matrix) :: this
		real(wp),dimension(:) :: apr_x,apr_y  
		!local var
		integer  :: li_i, li_k_1, li_k_2		
		
		apr_y = 0.0_wp
		do li_i = 1, this%oi_nR
		
			li_k_1 = this%opi_ia ( li_i )
			li_k_2 = this%opi_ia ( li_i + 1 ) - 1			

			apr_y ( li_i ) = DOT_PRODUCT ( this%opr_a ( li_k_1 : li_k_2 ), apr_x ( this%opi_ja ( li_k_1 : li_k_2 ) ) )
			 
			
			
		end do
		
	end subroutine Mult_CSR_Matrix_Vector	
!---------------------------------------------------------------------------------
	subroutine add_to_csr_Matrix ( this, ar_value, ai_A, ai_Aprime )	
	implicit none
		type(csr_matrix) :: this
		real(wp) :: ar_value
		integer  :: ai_A, ai_Aprime
		!local var
		integer :: li_result, li_i, li_j, li_k
		
		li_result = 0
		
		! THE CURRENT LINE IS this%opi_ia(ai_A)
		do li_k = this%opi_ia(ai_A),this%opi_ia(ai_A+1)-1
			li_j = this%opi_ja(li_k)
			if ( li_j == ai_Aprime ) then
				this%opr_a(li_k) = this%opr_a(li_k) + ar_value
				exit
			end if
		end do	
		
	end subroutine add_to_csr_Matrix
!---------------------------------------------------------------------------------
	subroutine add_csr_Matrix ( this, ao_csr )	
	implicit none
		type(csr_matrix) :: this, ao_csr
		
		this%opr_a = this%opr_a + ao_csr%opr_a
		
	end subroutine add_csr_Matrix	
!---------------------------------------------------------------------------------
	subroutine add_to_csr_Matrix_condition ( this, ar_value, ai_A, ai_Aprime, al_flag )	
	implicit none
		type(csr_matrix) :: this
		real(wp) :: ar_value
		integer  :: ai_A, ai_Aprime
		logical  :: al_flag
		!local var
		integer :: li_result, li_i, li_j, li_k
		
		li_result = 0
		
		! THE CURRENT LINE IS this%opi_ia(ai_A)
		do li_k = this%opi_ia(ai_A),this%opi_ia(ai_A+1)-1
			li_j = this%opi_ja(li_k)
			if ( ( li_j == ai_Aprime ) .AND. ( al_flag ) .AND. ( this%opr_a(li_k) == 0.0_wp ) ) then
				this%opr_a(li_k) = this%opr_a(li_k) + ar_value
				exit
			end if
		end do	
		
	end subroutine add_to_csr_Matrix_condition