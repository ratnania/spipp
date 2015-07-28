MODULE SPI_SPARSE_MATRIX_DEF 
USE SPI_GLOBAL_DEF 
IMPLICIT NONE

   ! ...................................................
   TYPE, PUBLIC, ABSTRACT :: DEF_SPARSE_MATRIX_ABSTRACT
      INTEGER  :: oi_nR   !NUMBER OF ROWS
      INTEGER  :: oi_nC   !NUMBER OF COLUMNS
      INTEGER  :: oi_nnz  !NUMBER OF NON ZERO ELTS
      INTEGER  :: oi_n_elmts !NUMBER OF MESH ELEMENTS
   END TYPE DEF_SPARSE_MATRIX_ABSTRACT
   ! ...................................................

   ! ...................................................
   TYPE, PUBLIC, EXTENDS(DEF_SPARSE_MATRIX_ABSTRACT) :: DEF_MATRIX_CSR
      LOGICAL  :: ol_allocated_ia 
      LOGICAL  :: ol_allocated_jaa 
      INTEGER, DIMENSION(:), POINTER  :: opi_ia
      INTEGER, DIMENSION(:), POINTER  :: opi_ja		

      REAL(SPI_RK), DIMENSION(:), POINTER  :: opr_a		

      !................
      logical :: ol_use_mm_format
      INTEGER, DIMENSION(:), POINTER  :: opi_i
      !................		
   END TYPE DEF_MATRIX_CSR
   ! ...................................................

   ! ...................................................
   TYPE, PUBLIC, EXTENDS(DEF_SPARSE_MATRIX_ABSTRACT) :: DEF_MATRIX_BND
      INTEGER  :: oi_n   !NUMBER OF ROWS
      INTEGER  :: oi_lowd   
      INTEGER  :: oi_ml   		
      INTEGER  :: oi_mu		
      INTEGER  :: oi_nabd  	
      INTEGER  :: oi_LDA		

      REAL(SPI_RK), DIMENSION(:,:), POINTER  :: opr_abd		
      INTEGER, DIMENSION(:),POINTER	:: opi_ipiv		
   END TYPE DEF_MATRIX_BND
   ! ...................................................

END MODULE SPI_SPARSE_MATRIX_DEF
