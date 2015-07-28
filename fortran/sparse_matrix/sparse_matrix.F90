MODULE SPI_SPARSE_MATRIX 
USE SPI_GLOBAL_DEF 
USE SPI_SPARSE_MATRIX_DEF
USE SPI_CSR_MATRIX_MODULE 
USE SPI_BND_MATRIX_MODULE 
IMPLICIT NONE

CONTAINS

   ! ...................................................
   SUBROUTINE create_sparse_matrix ( self, ai_nR, ai_nC, ai_nnz, ao_csr)
   IMPLICIT NONE
      CLASS(DEF_SPARSE_MATRIX_ABSTRACT)          , INTENT(INOUT) :: self
      INTEGER                          , OPTIONAL, INTENT(IN)    :: ai_nR
      INTEGER                          , OPTIONAL, INTENT(IN)    :: ai_nC
      INTEGER                          , OPTIONAL, INTENT(IN)    :: ai_nnz
      CLASS(DEF_MATRIX_CSR)            , OPTIONAL, INTENT(INOUT) :: ao_csr 
      ! LOCAL
      INTEGER :: li_err
      INTEGER :: li_flag
      
      ! ...
      SELECT TYPE (self)
      CLASS IS (DEF_MATRIX_CSR)
         CALL create_csr_matrix( self, ai_nR=ai_nR, ai_nC=ai_nC, ai_nnz=ai_nnz ) 
     
      CLASS IS (DEF_MATRIX_BND)
         IF (.NOT. PRESENT(ao_csr) ) THEN
            PRINT *, "create_sparse_matrix: try to create bnd matrix, missing arguments."
            STOP
         END IF
         CALL create_bnd_matrix(self, ao_csr)
     
      CLASS DEFAULT
         STOP 'create_sparse_matrix: unexpected type for self object!'
      END SELECT
      ! ...

   END SUBROUTINE create_sparse_matrix
   ! ...................................................

   ! ...................................................
   SUBROUTINE free_sparse_matrix ( self )
   IMPLICIT NONE
      CLASS(DEF_SPARSE_MATRIX_ABSTRACT), INTENT(INOUT) :: self
      ! LOCAL
      INTEGER :: li_err
      INTEGER :: li_flag
      
      ! ...
      SELECT TYPE (self)
      CLASS IS (DEF_MATRIX_CSR)
         CALL free_csr_matrix(self) 
     
      CLASS IS (DEF_MATRIX_BND)
         CALL free_bnd_matrix(self)
     
      CLASS DEFAULT
         STOP 'FREE_SPARSE_MATRIX: unexpected type for self object!'
      END SELECT
      ! ...

   END SUBROUTINE free_sparse_matrix
   ! ...................................................   

   ! ...................................................
   SUBROUTINE reset_sparse_matrix ( self )
   IMPLICIT NONE
      CLASS(DEF_SPARSE_MATRIX_ABSTRACT), INTENT(INOUT) :: self
      ! LOCAL
      INTEGER :: li_err
      INTEGER :: li_flag
      
      ! ...
      SELECT TYPE (self)
      CLASS IS (DEF_MATRIX_CSR)
         CALL reset_csr_matrix(self) 
     
      CLASS IS (DEF_MATRIX_BND)
         CALL reset_bnd_matrix(self)
     
      CLASS DEFAULT
         STOP 'RESET_SPARSE_MATRIX: unexpected type for self object!'
      END SELECT
      ! ...

   END SUBROUTINE reset_sparse_matrix
   ! ...................................................   

   ! ...................................................
   subroutine mult_sparse_matrix_vector(self, apr_x, apr_y)
   implicit none
      CLASS(DEF_SPARSE_MATRIX_ABSTRACT), INTENT(INOUT) :: self
      real(SPI_RK),dimension(:) :: apr_x
      real(SPI_RK),dimension(:) :: apr_y  
      ! LOCAL

      ! ...
      SELECT TYPE (self)
      CLASS IS (DEF_MATRIX_CSR)
         CALL mult_csr_matrix_vector(self, apr_x, apr_y) 
     
      CLASS IS (DEF_MATRIX_BND)
         STOP "mult_sparse_matrix_vector: not yet implemented"
     
      CLASS DEFAULT
         STOP 'mult_sparse_matrix_vector: unexpected type for self object!'
      END SELECT
      ! ...
      
   end subroutine mult_sparse_matrix_vector
   ! ...................................................

   ! ...................................................
   subroutine add_sparse_matrix_value ( self, ar_value, ai_A, ai_Aprime )	
   implicit none
      CLASS(DEF_SPARSE_MATRIX_ABSTRACT), INTENT(INOUT) :: self
      real(SPI_RK) :: ar_value
      integer  :: ai_A, ai_Aprime
      ! LOCAL

      ! ...
      SELECT TYPE (self)
      CLASS IS (DEF_MATRIX_CSR)
         CALL add_csr_matrix_value(self, ar_value, ai_A, ai_Aprime) 
     
      CLASS IS (DEF_MATRIX_BND)
         STOP "add_sparse_matrix_value: not yet implemented"
     
      CLASS DEFAULT
         STOP 'add_sparse_matrix_value: unexpected type for self object!'
      END SELECT
      ! ...
      
   end subroutine add_sparse_matrix_value
   ! ...................................................

   ! ...................................................
   subroutine add_sparse_matrix_sparse_matrix ( self, other)	
   implicit none
      CLASS(DEF_SPARSE_MATRIX_ABSTRACT), INTENT(INOUT) :: self
      CLASS(DEF_SPARSE_MATRIX_ABSTRACT), INTENT(INOUT) :: other 
      
      ! ...
      SELECT TYPE (self)
      CLASS IS (DEF_MATRIX_CSR)
         SELECT TYPE (other)
         CLASS IS (DEF_MATRIX_CSR)
            CALL add_csr_matrix_csr_matrix(self, other) 
        
         CLASS IS (DEF_MATRIX_BND)
            STOP 'add_sparse_matrix_sparse_matrix: unexpected type for other object!'
        
         CLASS DEFAULT
            STOP 'add_sparse_matrix_sparse_matrix: unexpected type for self object!'
         END SELECT
     
      CLASS IS (DEF_MATRIX_BND)
         STOP "add_sparse_matrix_sparse_matrix: not yet implemented"
     
      CLASS DEFAULT
         STOP 'add_sparse_matrix_sparse_matrix: unexpected type for self object!'
      END SELECT
      ! ...
      
   end subroutine add_sparse_matrix_sparse_matrix
   ! ...................................................
   
   ! ...................................................
   subroutine save_sparse_matrix(self, as_file, ai_format)
   implicit none
      CLASS(DEF_SPARSE_MATRIX_ABSTRACT), INTENT(INOUT) :: self
      character(len=*), intent(in) :: as_file		
      integer, intent(in) :: ai_format
      ! LOCAL

      ! ...
      SELECT TYPE (self)
      CLASS IS (DEF_MATRIX_CSR)
         IF (ai_format == SPI_MATRIX_OUTPUT_FORMAT_MM) THEN
            CALL save_csr_matrix_mm_format(self, as_file)
         ELSE
            STOP "save_sparse_matrix: format not yet supported"
         END IF
     
      CLASS IS (DEF_MATRIX_BND)
         STOP "save_sparse_matrix: not yet implemented"
     
      CLASS DEFAULT
         STOP 'save_sparse_matrix: unexpected type for self object!'
      END SELECT
      ! ...

   end subroutine save_sparse_matrix
   ! ...................................................

   ! ...................................................
   SUBROUTINE initialize_sparse_matrix_with_LM(self, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2 )
   IMPLICIT NONE
      CLASS(DEF_SPARSE_MATRIX_ABSTRACT), INTENT(INOUT) :: self
      !> param[in] api_LM_1 : LM ARRAY FOR ROWS
      integer, dimension(:,:) :: api_LM_1
      !> param[in] api_LM_1 : LM ARRAY FOR COLUMNS
      integer, dimension(:,:) :: api_LM_2					
      !> param[in] ai_nen_1 : NUMBER OF NON VANISHING FUNCTIONS PER ELEMENT, IN THE 1st SPACE 
      integer :: ai_nen_1
      !> param[in] ai_nen_2 : NUMBER OF NON VANISHING FUNCTIONS PER ELEMENT, IN THE 2nd SPACE
      integer :: ai_nen_2					
      ! LOCAL
      INTEGER :: li_err
      INTEGER :: li_flag
      
      ! ...
      SELECT TYPE (self)
      CLASS IS (DEF_MATRIX_CSR)
         CALL initialize_csr_matrix_with_LM (self, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2 )
     
      CLASS IS (DEF_MATRIX_BND)
         STOP "initialize_sparse_matrix_with_LM: not yet implemented"
     
      CLASS DEFAULT
         STOP 'initialize_sparse_matrix_with_LM: unexpected type for self object!'
      END SELECT
      ! ...

   END SUBROUTINE initialize_sparse_matrix_with_LM
   ! ...................................................   
   
   ! ...................................................

END MODULE SPI_SPARSE_MATRIX
