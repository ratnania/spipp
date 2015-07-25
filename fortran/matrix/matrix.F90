!# -*- coding: utf8 -*-
MODULE SPI_MATRIX
  USE SPI_MATRIX_DEF
  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 2 

  ! .............................................
  
  ! .............................................
  PRIVATE
  PUBLIC :: MATRIX_CREATE               &
       ,    GET_NR_MATRIX               &
       ,    GET_NC_MATRIX               &
       ,    MATRIX_RESET_ELEMENT_MATRIX &
       ,    MATRIX_RESET_ELEMENT_RHS    &
       ,    MATRIX_GET_GLOBAL_ROW_INDEX &
       ,    MATRIX_FREE
  
CONTAINS
  ! ..................................................

  ! ..................................................
  SUBROUTINE MATRIX_CREATE(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT), INTENT(INOUT) :: self 
     ! LOCAL
     INTEGER :: li_matrix_id
     INTEGER :: li_err
    
     self % oi_matrix_id   = -1
     self % ol_to_assembly = .FALSE.

  END SUBROUTINE MATRIX_CREATE
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

  ! ..................................................
  SUBROUTINE MATRIX_GET_GLOBAL_ROW_INDEX(self, ai_glob, ai_real)
    IMPLICIT NONE
    CLASS(DEF_MATRIX_ABSTRACT) :: self
    INTEGER                   :: ai_glob
    INTEGER                   :: ai_real
    !LOCAL
    INTEGER                   :: li_err

  END SUBROUTINE MATRIX_GET_GLOBAL_ROW_INDEX
  ! ..................................................
    
  ! ..........................................................        
  SUBROUTINE MATRIX_RESET_ELEMENT_MATRIX(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL

#ifdef DEBUG_TRACE 
    CALL printlog("MATRIX_RESET_ELEMENT_MATRIX: Begin", ai_dtllevel = 2)
#endif     

     self % Matrix_Contribution = 0.0
     
#ifdef DEBUG_TRACE 
CALL printlog("MATRIX_RESET_ELEMENT_MATRIX: End", ai_dtllevel = 2)
#endif

  END SUBROUTINE MATRIX_RESET_ELEMENT_MATRIX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE MATRIX_RESET_ELEMENT_RHS(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL

#ifdef DEBUG_TRACE 
    CALL printlog("MATRIX_RESET_ELEMENT_RHS: Begin", ai_dtllevel = 2)
#endif     

     self % Rhs_Contribution    = 0.0
     
#ifdef DEBUG_TRACE 
CALL printlog("MATRIX_RESET_ELEMENT_RHS: End", ai_dtllevel = 2)
#endif

  END SUBROUTINE MATRIX_RESET_ELEMENT_RHS
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE MATRIX_FREE(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL

#ifdef DEBUG_TRACE 
    CALL printlog("MATRIX_FREE: Begin", ai_dtllevel = 0)
#endif     

     DEALLOCATE(self % Matrix_Contribution)
     DEALLOCATE(self % Rhs_Contribution   )
     
#ifdef DEBUG_TRACE 
CALL printlog("MATRIX_FREE: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE MATRIX_FREE
  ! ..........................................................  
END MODULE SPI_MATRIX 
