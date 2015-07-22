MODULE SPI_BND_MATRIX_MODULE
USE SPI_SPARSE_MATRIX_DEF 
USE SPI_SPARSEKIT
IMPLICIT NONE

CONTAINS
   ! ...................................................
   subroutine create_bnd_matrix ( self, ao_csr )
   implicit none
      type(DEF_MATRIX_BND) :: self
      type(DEF_MATRIX_CSR) :: ao_csr
      ! LOCAL VARIABLES		
      integer :: li_flag, li_err, li_ierr, li_info		
      
      if ( ao_csr%oi_nR /= ao_csr%oi_nC ) then
         print*,'Error create_bnd_matrix : Cannot solve the linear system, the matrix is not square'
         stop
      end if
      
      self%oi_n = ao_csr%oi_nR
      
      call getbwd ( self%oi_n, ao_csr%opr_a, ao_csr%opi_ja, ao_csr%opi_ia, self%oi_ml, self%oi_mu )
      
      self%oi_LDA  = 2 * self%oi_ml + self%oi_mu + 1	
      self%oi_nabd = 2 * self%oi_ml + self%oi_mu + 1	
      self%oi_lowd = 2 * self%oi_ml + self%oi_mu + 1	
      
      allocate ( self%opr_abd ( self%oi_nabd , self%oi_n ),stat=li_err)
      if (li_err.ne.0) li_flag=10	
      
      allocate ( self%opi_ipiv ( self%oi_n ),stat=li_flag)
      if (li_err.ne.0) li_flag=11	
      
      self%opr_abd = 0.0
      call csrbnd ( self%oi_n, ao_csr%opr_a, ao_csr%opi_ja, ao_csr%opi_ia, 1		&
              , self%opr_abd, self%oi_nabd, self%oi_lowd, self%oi_ml, self%oi_mu, li_ierr )
      
      call DGBTRF ( self%oi_n, self%oi_n, self%oi_ml, self%oi_mu		&
              , self%opr_abd, self%oi_LDA, self%opi_ipiv, li_info )
   
   end subroutine create_bnd_matrix
   ! ...................................................

   ! ...................................................
   subroutine free_bnd_matrix ( self )
   implicit none
      type(DEF_MATRIX_BND) :: self
      
      deallocate(self%opr_abd)
      deallocate(self%opi_ipiv)
   
   end subroutine free_bnd_matrix		
   ! ...................................................
   
   ! ...................................................
   subroutine BND_solve_v ( self, apr_B, apr_U )
   implicit none
      type(DEF_MATRIX_BND) :: self 	
      real(SPI_RK), dimension(:) :: apr_U
      real(SPI_RK), dimension(:) :: apr_B	  
      ! LOCAL VARIABLES
      integer	 :: li_info	 

      apr_U = apr_B
      
      call DGBTRS ( 'N' , self%oi_n, self%oi_ml, self%oi_mu, 1		        &
              , self%opr_abd, self%oi_LDA, self%opi_ipiv			&
              , apr_U, self%oi_n, li_info )
   
   end subroutine BND_solve_v
   ! ...................................................
   
   ! ...................................................
   subroutine BND_solve_m ( self, apr_B, apr_U )
   implicit none
      type(DEF_MATRIX_BND) :: self 	
      real(SPI_RK), dimension(:,:) :: apr_U
      real(SPI_RK), dimension(:,:) :: apr_B	  
      ! LOCAL VARIABLES
      integer	 :: li_info	 
      integer  :: li_m
      
      li_m = size ( apr_B , 2 )
      
      apr_U = apr_B
      
      call DGBTRS ( 'N' , self%oi_n, self%oi_ml, self%oi_mu, li_m		&
              , self%opr_abd, self%oi_LDA, self%opi_ipiv			&
              , apr_U, self%oi_n, li_info )
   
   end subroutine BND_solve_m
   ! ...................................................

   ! ...................................................
END MODULE SPI_BND_MATRIX_MODULE
