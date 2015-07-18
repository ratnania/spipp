module bnd_Matrix_module
use used_precision
use ISOBOX_Def_Module
implicit none

contains
!---------------------------------------------------------------------------------
	subroutine create_bnd_matrix ( this, ao_csr )
	implicit none
		type(bnd_matrix) :: this
		type(csr_matrix) :: ao_csr
		! LOCAL VARIABLES		
		integer :: li_flag, li_err, li_ierr, li_info		
		
		if ( ao_csr%oi_nR /= ao_csr%oi_nC ) then
			print*,'Error create_bnd_matrix : Cannot solve the linear system, the matrix is not square'
			stop
		end if

		this%oi_n = ao_csr%oi_nR

		call getbwd ( this%oi_n, ao_csr%opr_a, ao_csr%opi_ja, ao_csr%opi_ia, this%oi_ml, this%oi_mu )

		this%oi_LDA  = 2 * this%oi_ml + this%oi_mu + 1	
		this%oi_nabd = 2 * this%oi_ml + this%oi_mu + 1	
		this%oi_lowd = 2 * this%oi_ml + this%oi_mu + 1	

		allocate ( this%opr_abd ( this%oi_nabd , this%oi_n ),stat=li_err)
		if (li_err.ne.0) li_flag=10	

		allocate ( this%opi_ipiv ( this%oi_n ),stat=li_flag)
		if (li_err.ne.0) li_flag=11	

		this%opr_abd = 0.0_wp
		call csrbnd ( this%oi_n, ao_csr%opr_a, ao_csr%opi_ja, ao_csr%opi_ia, 1		&
					, this%opr_abd, this%oi_nabd, this%oi_lowd, this%oi_ml, this%oi_mu, li_ierr )

		call DGBTRF ( this%oi_n, this%oi_n, this%oi_ml, this%oi_mu		&
					, this%opr_abd, this%oi_LDA, this%opi_ipiv, li_info )

	end subroutine create_bnd_matrix
!---------------------------------------------------------------------------------
	subroutine free_bnd_matrix ( this )
	implicit none
		type(bnd_matrix) :: this
		
		deallocate(this%opr_abd)
		deallocate(this%opi_ipiv)

	end subroutine free_bnd_matrix		
!---------------------------------------------------------------------------------
	subroutine BND_solve_v ( this, apr_B, apr_U )
	implicit none
		type(bnd_matrix) :: this 	
		real(wp), dimension(:) :: apr_U
		real(wp), dimension(:) :: apr_B	  
		! LOCAL VARIABLES
		integer	 :: li_info	 

		apr_U = apr_B
	
		call DGBTRS ( 'N' , this%oi_n, this%oi_ml, this%oi_mu, 1		&
					, this%opr_abd, this%oi_LDA, this%opi_ipiv			&
					, apr_U, this%oi_n, li_info )

	end subroutine BND_solve_v
!---------------------------------------------------------------------------------
	subroutine BND_solve_m ( this, apr_B, apr_U )
	implicit none
		type(bnd_matrix) :: this 	
		real(wp), dimension(:,:) :: apr_U
		real(wp), dimension(:,:) :: apr_B	  
		! LOCAL VARIABLES
		integer	 :: li_info	 
		integer  :: li_m
		
		li_m = size ( apr_B , 2 )

		apr_U = apr_B
	
		call DGBTRS ( 'N' , this%oi_n, this%oi_ml, this%oi_mu, li_m		&
					, this%opr_abd, this%oi_LDA, this%opi_ipiv			&
					, apr_U, this%oi_n, li_info )

	end subroutine BND_solve_m
!---------------------------------------------------------------------------------
!	subroutine inverseMat( ai_sizePB , apr_A , apr_invA )
!	!calcul l inverse de sigma
!	! X et Y sont deux vecteurs 
!	implicit none
!		integer  :: ai_sizePB
!		double precision, dimension(:,:),intent(in) :: apr_A 
!		double precision, dimension(:,:),intent(inout) :: apr_invA		
!		!var local
!		integer :: li_flag,li_err
!		integer	:: li_info	 !util pour la routine DGETRI
!		integer, parameter  :: li_MaxWorkSize = 5000		
!			
!		apr_invA = apr_A
!		
!		call DGBTRS ( 'N' , this%oi_n, this%oi_ml, this%oi_mu, 1		&
!					, this%opr_abd, this%oi_LDA, this%opi_ipiv			&
!					, apr_U, this%oi_n, li_info )		
!
!		call DGBTRF ( this%oi_n, this%oi_n, this%oi_ml, this%oi_mu		&
!					, this%opr_abd, this%oi_LDA, this%opi_ipiv, li_info )
!
!		call DGBTRI( this%oi_n, apr_invA, ai_sizePB, lpi_ipiv, lpd_work, li_MaxWorkSize, li_info )		
!
!	end subroutine inverseMat	
end module bnd_Matrix_module	