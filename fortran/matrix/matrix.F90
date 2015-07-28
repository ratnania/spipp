!# -*- coding: utf8 -*-
MODULE SPI_MATRIX
  USE SPI_GLOBAL_DEF
  USE SPI_MATRIX_DEF
  USE SPI_SPARSE_MATRIX_DEF
  USE SPI_SPARSE_MATRIX
  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 2 

  ! .............................................
  
  ! .............................................
  PRIVATE
  PUBLIC :: CREATE_MATRIX               &
       ,    GET_NR_MATRIX               &
       ,    GET_NC_MATRIX               &
       ,    RESET_MATRIX_MATRIX         &
       ,    RESET_MATRIX_RHS            &
       ,    RESET_ELEMENT_MATRIX        &
       ,    RESET_ELEMENT_RHS_MATRIX    &
       ,    ASSEMBLYBEGIN_MATRIX        &
       ,    ASSEMBLYEND_MATRIX          &
       ,    GLOBAL_TO_INDEX_MATRIX      &
       ,    FREE_MATRIX
  
CONTAINS
  ! ..................................................

  ! ..................................................
  SUBROUTINE CREATE_MATRIX(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT), INTENT(INOUT) :: self 
     ! LOCAL
     INTEGER :: li_matrix_id
     INTEGER :: li_err
    
     self % oi_matrix_id   = -1
     self % ol_to_assembly = .FALSE.
     self % ol_allocated_locmat = .FALSE. 
     self % oi_npush_locmat = 0

  END SUBROUTINE CREATE_MATRIX
  ! ..................................................

  ! ..................................................
  SUBROUTINE GET_NR_MATRIX(self, ai_value)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT), INTENT(INOUT) :: self 
     INTEGER, INTENT(INOUT)           :: ai_value
     ! LOCAL
     INTEGER :: li_err

  END SUBROUTINE GET_NR_MATRIX
  ! ..................................................

  ! ..................................................
  SUBROUTINE GET_NC_MATRIX(self, ai_value)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT), INTENT(INOUT) :: self 
     INTEGER, INTENT(INOUT)           :: ai_value
     ! LOCAL
     INTEGER :: li_err

  END SUBROUTINE GET_NC_MATRIX
  ! ..................................................

  ! ..........................................................        
  SUBROUTINE SET_NDOF_ROWS_MATRIX(self, ai_value)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     INTEGER, INTENT(IN) :: ai_value
     ! LOCAL

     self % oi_ndof_rows = ai_value 
     
  END SUBROUTINE SET_NDOF_ROWS_MATRIX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE SET_NDOF_COLS_MATRIX(self, ai_value)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     INTEGER, INTENT(IN) :: ai_value
     ! LOCAL

     self % oi_ndof_cols = ai_value 
     
  END SUBROUTINE SET_NDOF_COLS_MATRIX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE SET_LOCALSIZE_ROWS_MATRIX(self, ai_value)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     INTEGER, INTENT(IN) :: ai_value
     ! LOCAL

     self % oi_localsize_rows = ai_value 
     
  END SUBROUTINE SET_LOCALSIZE_ROWS_MATRIX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE SET_LOCALSIZE_COLS_MATRIX(self, ai_value)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     INTEGER, INTENT(IN) :: ai_value
     ! LOCAL

     self % oi_localsize_cols = ai_value 
     
  END SUBROUTINE SET_LOCALSIZE_COLS_MATRIX
  ! ..........................................................  

  ! ..................................................
  SUBROUTINE GLOBAL_TO_INDEX_MATRIX(self, ai_glob, ai_real)
    IMPLICIT NONE
    CLASS(DEF_MATRIX_ABSTRACT) :: self
    INTEGER                   :: ai_glob
    INTEGER                   :: ai_real
    !LOCAL
    INTEGER                   :: li_err

  END SUBROUTINE GLOBAL_TO_INDEX_MATRIX
  ! ..................................................

  ! ..........................................................        
  SUBROUTINE RESET_MATRIX_MATRIX(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL
     
  END SUBROUTINE RESET_MATRIX_MATRIX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE RESET_MATRIX_RHS(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL

  END SUBROUTINE RESET_MATRIX_RHS
  ! ..........................................................  
    
  ! ..........................................................        
  SUBROUTINE RESET_ELEMENT_MATRIX(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL

     self % Matrix_Contribution = 0.0
     
  END SUBROUTINE RESET_ELEMENT_MATRIX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE RESET_ELEMENT_RHS_MATRIX(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL

     self % Rhs_Contribution    = 0.0
     
  END SUBROUTINE RESET_ELEMENT_RHS_MATRIX
  ! ..........................................................  
    
  ! ..........................................................        
  SUBROUTINE RESET_LOCAL_MATRIX(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL

    IF ( .NOT. self % ol_allocated_locmat ) THEN
       RETURN
    END IF

    self % opr_locmat = 0.0
     
  END SUBROUTINE RESET_LOCAL_MATRIX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE ASSEMBLYBEGIN_MATRIX(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL
     INTEGER :: li_ndof_rows 
     INTEGER :: li_ndof_cols
     INTEGER :: li_localsize_rows 
     INTEGER :: li_localsize_cols

     IF (.NOT. self % ol_allocated_locmat) THEN
        li_ndof_rows      = self % oi_ndof_rows 
        li_ndof_cols      = self % oi_ndof_cols 
        li_localsize_rows = self % oi_localsize_rows
        li_localsize_cols = self % oi_localsize_cols
        IF ( li_localsize_rows * li_localsize_cols > 0) THEN
           ALLOCATE(self % opr_locmat(li_localsize_rows*li_ndof_rows, &
                                    & li_localsize_cols*li_ndof_cols))
           self % ol_allocated_locmat = .TRUE.
           self % opr_locmat = 0.0
        END IF
     END IF

  END SUBROUTINE ASSEMBLYBEGIN_MATRIX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE ASSEMBLYEND_MATRIX(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL

  END SUBROUTINE ASSEMBLYEND_MATRIX
  ! ..........................................................  

  ! ..........................................................        
!  SUBROUTINE _MATRIX(self)
!  IMPLICIT NONE
!     CLASS(DEF_MATRIX_ABSTRACT) :: self
!     ! LOCAL
!
!  END SUBROUTINE _MATRIX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE FREE_MATRIX(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL

     DEALLOCATE(self % Matrix_Contribution)
     DEALLOCATE(self % Rhs_Contribution   )

     IF (self % ol_allocated_locmat) THEN
        DEALLOCATE(self % opr_locmat)
        self % ol_allocated_locmat = .FALSE.
     END IF

  END SUBROUTINE FREE_MATRIX
  ! ..........................................................  

  ! .............................................
  SUBROUTINE ASSEMBLYPUSH_MATRIX(self, ROWS, COLS, IERROR)
  ! pushes add value from localdata to the matrix
    !(ID, IERROR)
    IMPLICIT NONE
    CLASS(DEF_MATRIX_ABSTRACT) :: self
    INTEGER, DIMENSION(:,:), POINTER, INTENT(IN)  :: ROWS
    INTEGER, DIMENSION(:,:), POINTER, INTENT(IN)  :: COLS
    INTEGER,                          INTENT(OUT) :: IERROR
    ! LOCAL
    INTEGER :: li_ierr 
    INTEGER :: li_I
    INTEGER :: li_J
    INTEGER :: li_lrow
    INTEGER :: li_lcol
    INTEGER :: li_iloc
    INTEGER :: li_jloc
    INTEGER :: li_idof
    INTEGER :: li_jdof
    INTEGER :: li_loc
    INTEGER :: li_dof
    INTEGER :: li_err
    REAL(KIND=SPI_RK)    :: lr_val
    LOGICAL :: ll_condition
    INTEGER :: li_ndof_rows 
    INTEGER :: li_ndof_cols
    INTEGER :: li_localsize_rows 
    INTEGER :: li_localsize_cols           

    ! ... TODO
    li_localsize_rows = SPI_INT_DEFAULT
    li_localsize_cols = SPI_INT_DEFAULT 
    li_ndof_rows      = SPI_INT_DEFAULT 
    li_ndof_cols      = SPI_INT_DEFAULT 

    DO li_lcol=1,li_localsize_cols
       DO li_jdof=1, li_ndof_cols
          li_jloc = li_jdof + (li_lcol-1) * li_ndof_cols

          DO li_lrow=1,li_localsize_rows
             DO li_idof=1, li_ndof_rows
                li_iloc = li_idof + (li_lrow-1) * li_ndof_rows

                li_I   = ROWS(li_iloc, li_jloc)
                li_J   = COLS(li_iloc, li_jloc)
                lr_val = self % opr_locmat(li_iloc, li_jloc)

                CALL ADD_SPARSE_MATRIX_VALUE(self % oo_csr, lr_val, li_I, li_J) 

             END DO
          END DO
       END DO
    END DO
    ! ... 

    ! ... 
    ! Reset
    self % opr_locmat = 0.0
    self % oi_npush_locmat = self % oi_npush_locmat + 1
    ! ... 

  END SUBROUTINE ASSEMBLYPUSH_MATRIX
  ! .............................................


END MODULE SPI_MATRIX 
