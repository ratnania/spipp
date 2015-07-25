MODULE SPI_CSR_MATRIX_MODULE 
USE SPI_GLOBAL_DEF 
USE SPI_SPARSE_MATRIX_DEF
USE SPI_QSORT 
IMPLICIT NONE

CONTAINS

   ! ...................................................
   SUBROUTINE create_csr_matrix( self, ai_nR, ai_nC, ai_nel )
   IMPLICIT NONE
      !> param[inout] self : CSR MATRIX STRUCTURE
      type(DEF_MATRIX_CSR) :: self		
      !> param[in] ai_nC : NUMBER OF COLUMNS
      INTEGER :: ai_nC
      !> param[in] ai_nR : NUMBER OF ROWS
      INTEGER :: ai_nR
      !> param[in] ai_nel : NUMBER OF NON ZERO ELEMENTS		
      INTEGER :: ai_nel				
      !local var
      INTEGER :: li_err,li_flag
      
      self%ol_use_mm_format = .FALSE.
      
      self%oi_nR   = ai_nR
      self%oi_nC   = ai_nC		
      self%oi_nel  = ai_nel
      
      ALLOCATE(self%opi_ia(self%oi_nR+1),stat=li_err)
      IF (li_err.ne.0) li_flag=10	
      
      ALLOCATE(self%opi_ja(self%oi_nel),stat=li_err)
      IF (li_err.ne.0) li_flag=20	
   		
      ALLOCATE(self%opr_a(self%oi_nel),stat=li_err)
      IF (li_err.ne.0) li_flag=30			

      self%opr_a ( : ) = 0.0
   END SUBROUTINE create_csr_matrix
   ! ...................................................

   ! ...................................................
   ! THIS ROUTINE CREATES A NEW FULL MATRIX FROM A COMPRESSED SYMMETRIC MATRIX
   subroutine create_csr_matrix_from_symmetric( self, ao_A )
   implicit none
      !> param[inout] self : CSR MATRIX STRUCTURE	
      type(DEF_MATRIX_CSR) :: self
      !> param[in] ao_A : CSR MATRIX STRUCTURE		
      type(DEF_MATRIX_CSR) :: ao_A		
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
      
      self%ol_use_mm_format = .FALSE.		
      
      self%oi_nel = 2 * ao_A%oi_nel - li_count
      self%oi_nR = ao_A%oi_nR
      self%oi_nC = ao_A%oi_nC	
      
      ALLOCATE(self%opi_ia(self%oi_nR+1),stat=li_err)
      if (li_err.ne.0) li_flag=10	
      
      ALLOCATE(self%opi_ja(self%oi_nel),stat=li_err)
      if (li_err.ne.0) li_flag=20	
                      
      ALLOCATE(self%opr_a(self%oi_nel),stat=li_err)
      if (li_err.ne.0) li_flag=30	
      		
      ! COPY THE OLD MATRIX IN THE NEW MATRIX		
      do li_i = 1, ao_A%oi_nR
         lpi_occ ( li_i ) = ao_A%opi_ia ( li_i + 1 ) - ao_A%opi_ia ( li_i )
         do li_i_tmp = 1, li_i - 1
            li_k = ao_A%opi_ia ( li_i_tmp )
            do while ( ( li_k <= ao_A%opi_ia ( li_i_tmp + 1 ) - 1 ) &
                    & .AND. ( ao_A%opi_ja ( li_k ) <= li_i ) )			
               if ( li_i == ao_A%opi_ja ( li_k ) ) then
                  lpi_occ ( li_i ) = lpi_occ ( li_i ) + 1
               end if
               li_k = li_k + 1
            end do
         end do
      end do
      
      self%opi_ia ( 1 ) = 1
      
      do li_i = 1, self%oi_nR 
         self%opi_ia ( li_i + 1 ) = self%opi_ia ( 1 ) + SUM ( lpi_occ ( 1 : li_i ) )
      end do
      
      li_index = 1
      do li_i = 1, ao_A%oi_nR
         do li_i_tmp = 1, li_i - 1
            li_k = ao_A%opi_ia ( li_i_tmp )
            do while ( ( li_k <= ao_A%opi_ia ( li_i_tmp + 1 ) - 1 ) &
                    & .AND. ( ao_A%opi_ja ( li_k ) <= li_i ) )			
                if ( li_i == ao_A%opi_ja ( li_k ) ) then
                  self%opi_ja ( li_index ) = li_i_tmp
                  li_index = li_index + 1
                end if
                li_k = li_k + 1
             end do
          end do
      
          self%opi_ja ( li_index : li_index + ao_A%opi_ia ( li_i + 1 ) - 1 - ao_A%opi_ia ( li_i ) ) = &
          ao_A%opi_ja ( ao_A%opi_ia ( li_i ) : ao_A%opi_ia ( li_i + 1 ) - 1 )	
      
          li_index = li_index + ao_A%opi_ia ( li_i + 1 ) - ao_A%opi_ia ( li_i )
      end do
      
   end subroutine create_csr_matrix_from_symmetric
   ! ...................................................

   ! ...................................................
   subroutine free_csr_matrix(self)
   implicit none
      type(DEF_MATRIX_CSR) :: self
      
      deallocate(self%opi_ia)
      deallocate(self%opi_ja)
      deallocate(self%opr_a)
      
      if ( self%ol_use_mm_format ) then
         DEALLOCATE ( self%opi_i )
      end if
   
   end subroutine free_csr_matrix
   ! ...................................................

   ! ...................................................
   integer function count_non_zero_elts ( ai_nR, ai_nel, &
                                        & api_LM_1, ai_nen_1, &
                                        & api_LM_2, ai_nen_2, &
                                        & api_columns, api_occ )
   ! _1 FOR ROWS
   ! _2 FOR COLUMNS	
   implicit none
      integer :: ai_nR
      integer, dimension(:,:) :: api_LM_1, api_LM_2			
      integer :: ai_nel, ai_nen_1, ai_nen_2
      integer, dimension(:,:), pointer :: api_columns	
      integer, dimension(:), pointer :: api_occ	
      !local var
      integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_i			
      integer :: li_err,li_flag			
      real(SPI_RK), dimension(:), pointer :: lpr_tmp	
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
      print *, ">>>>>>>>>"
      
      ! COUNT NON ZERO ELEMENTS		
      li_result = SUM ( api_occ ( 1 : ai_nR ) )
               print *, "pass"
      
      count_non_zero_elts = li_result
   end function count_non_zero_elts	   
   ! ...................................................

   ! ...................................................
   subroutine initialize_csr_matrix_with_LM ( self, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2 )
   implicit none
      !> param[inout] self : CSR MATRIX STRUCTURE		
      type(DEF_MATRIX_CSR) :: self
      INTEGER :: ai_nel
      !> param[in] api_LM_1 : LM ARRAY FOR ROWS
      integer, dimension(:,:) :: api_LM_1
      !> param[in] api_LM_1 : LM ARRAY FOR COLUMNS
      integer, dimension(:,:) :: api_LM_2					
      !> param[in] ai_nen_1 : NUMBER OF NON VANISHING FUNCTIONS PER ELEMENT, IN THE 1st SPACE 
      integer :: ai_nen_1
      !> param[in] ai_nen_2 : NUMBER OF NON VANISHING FUNCTIONS PER ELEMENT, IN THE 2nd SPACE
      integer :: ai_nen_2					
      !local var
      integer :: li_err,li_flag
      integer :: li_nnz
      integer, dimension(:,:), pointer :: lpi_columns	
      integer, dimension(:), pointer :: lpi_occ			
      
      self % oi_nR =  MAXVAL(api_LM_1)
      self % oi_nC =  MAXVAL(api_LM_2)
      self % oi_nel = ai_nel

      allocate(lpi_columns(self%oi_nR, 0:8*ai_nen_2),stat=li_err)
      if (li_err.ne.0) li_flag=1	
      
      allocate(lpi_occ(self%oi_nR+1),stat=li_err)
      if (li_err.ne.0) li_flag=2	
                              
      lpi_columns ( :, : ) = 0
      lpi_occ ( : ) = 0
      
      ! COUNTING NON ZERO ELEMENTS
      li_nnz = count_non_zero_elts ( self % oi_nR, ai_nel, &
              & api_LM_1, ai_nen_1, &
              & api_LM_2, ai_nen_2, &
              & lpi_columns, lpi_occ )
      print *, "nnz ", li_nnz
      
      call init_SparseMatrix ( self, self%oi_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, lpi_columns, lpi_occ )					
      deallocate(lpi_columns)		
      deallocate(lpi_occ)		
   end subroutine initialize_csr_matrix_with_LM
   ! ...................................................
     
   ! ...................................................
   subroutine init_SparseMatrix ( self, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, api_columns, api_occ )
   ! _1 FOR ROWS
   ! _2 FOR COLUMNS	
   implicit none
      type(DEF_MATRIX_CSR) :: self
      integer, dimension(:,:) :: api_LM_1, api_LM_2			
      integer :: ai_nel, ai_nen_1, ai_nen_2
      integer, dimension(:,:), pointer :: api_columns	
      integer, dimension(:), pointer :: api_occ	
      !local var
      integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_index, li_i	, li_size		
      integer :: li_err,li_flag				
      real(SPI_RK), dimension(:), pointer :: lpr_tmp							
      
      ! INITIALIZING ia 
      self%opi_ia ( 1 ) = 1
      
      do li_i = 1, self%oi_nR 
         self%opi_ia ( li_i + 1 ) = self%opi_ia ( 1 ) + SUM ( api_occ ( 1 : li_i ) )
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
               self%opi_ja ( self%opi_ia ( li_A_1 ) + li_i - 1 ) = int ( lpr_tmp ( li_i ) )	
            end do
      
            api_columns ( li_A_1, 0 ) = 0					
            deallocate ( lpr_tmp )

         end do
      end do	
      
   end subroutine init_SparseMatrix	
   ! ...................................................

   ! ...................................................
   subroutine mult_csr_matrix_vector(self,apr_x,apr_y)
   implicit none
      type(DEF_MATRIX_CSR) :: self
      real(SPI_RK),dimension(:) :: apr_x
      real(SPI_RK),dimension(:) :: apr_y  
      !local var
      integer  :: li_i, li_k_1, li_k_2		
      
      apr_y = 0.0
      do li_i = 1, self%oi_nR
         li_k_1 = self%opi_ia ( li_i )
         li_k_2 = self%opi_ia ( li_i + 1 ) - 1			
         
         apr_y ( li_i ) = DOT_PRODUCT ( self%opr_a ( li_k_1 : li_k_2 ), apr_x ( self%opi_ja ( li_k_1 : li_k_2 ) ) )
      end do
      
   end subroutine mult_csr_matrix_vector
   ! ...................................................

   ! ...................................................
   subroutine add_csr_matrix_value ( self, ar_value, ai_A, ai_Aprime )	
   implicit none
      type(DEF_MATRIX_CSR) :: self
      real(SPI_RK) :: ar_value
      integer  :: ai_A, ai_Aprime
      !local var
      integer :: li_result, li_i, li_j, li_k
      
      li_result = 0
      
      ! THE CURRENT LINE IS self%opi_ia(ai_A)
      do li_k = self%opi_ia(ai_A),self%opi_ia(ai_A+1)-1
         li_j = self%opi_ja(li_k)
         if ( li_j == ai_Aprime ) then
            self%opr_a(li_k) = self%opr_a(li_k) + ar_value
            exit
         end if
      end do	
      
   end subroutine add_csr_matrix_value
   ! ...................................................

   ! ...................................................
   subroutine add_csr_matrix_csr_matrix ( self, ao_csr )	
   implicit none
      type(DEF_MATRIX_CSR) :: self, ao_csr
      
      self%opr_a = self%opr_a + ao_csr%opr_a
      
   end subroutine add_csr_matrix_csr_matrix	
   ! ...................................................
   
   ! ...................................................
   subroutine save_csr_matrix_mm_format(self, as_file)
   implicit none
      type(DEF_MATRIX_CSR) :: self 
      character(len=*), intent(in) :: as_file		
      !local var
      integer  :: li_i,li_j,li_k	
      integer  :: li_file=1
      integer  :: li_ios

      open(unit=li_file, file=TRIM(as_file), iostat=li_ios)
      if (li_ios/=0) STOP "print_nnz_csrMatrix : erreur d'ouverture du fichier "
              
      write(li_file, *)'%%MatrixMarket matrix coordinate real general'
      write(li_file, *)self%oi_nR,',',self%oi_nC,',',self%oi_nel

      do li_i =1,self%oi_nR 
         do li_k = self%opi_ia(li_i),self%opi_ia(li_i+1)-1
            li_j = self%opi_ja(li_k)

            write(li_file, *)li_i,',',li_j,',',self%opr_a(li_k)		
         end do
      end do	

      close(li_file)

   end subroutine save_csr_matrix_mm_format
   ! ...................................................
       
END MODULE SPI_CSR_MATRIX_MODULE
