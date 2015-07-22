MODULE SPI_SPARSE_MATRIX 
USE SPI_GLOBAL_DEF 
USE SPI_SPARSE_MATRIX_DEF
USE SPI_CSR_MATRIX_MODULE 
USE SPI_BND_MATRIX_MODULE 
IMPLICIT NONE

CONTAINS

   ! ...................................................
   SUBROUTINE create_sparse_matrix ( self, ai_nR, ai_nC, ai_nel, ao_csr)
   IMPLICIT NONE
      CLASS(DEF_SPARSE_MATRIX_ABSTRACT)          , INTENT(INOUT) :: self
      INTEGER                          , OPTIONAL, INTENT(IN)    :: ai_nR
      INTEGER                          , OPTIONAL, INTENT(IN)    :: ai_nC
      INTEGER                          , OPTIONAL, INTENT(IN)    :: ai_nel
      CLASS(DEF_MATRIX_CSR)            , OPTIONAL, INTENT(INOUT) :: ao_csr 
      ! LOCAL
      INTEGER :: li_err
      INTEGER :: li_flag
      
      ! ...
      SELECT TYPE (self)
      CLASS IS (DEF_MATRIX_CSR)
         IF ( (.NOT. PRESENT(ai_nR) ) .OR. &
            & (.NOT. PRESENT(ai_nC) ) .OR. &
            & (.NOT. PRESENT(ai_nel) ) ) THEN
            PRINT *, "create_sparse_matrix: try to create csr matrix, missing arguments."
            STOP
         END IF
         CALL create_csr_matrix( self, ai_nR, ai_nC, ai_nel ) 
     
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

END MODULE SPI_SPARSE_MATRIX
      
      
      
      
      
      

      
      
