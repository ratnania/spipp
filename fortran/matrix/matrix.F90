!# -*- coding: utf8 -*-
MODULE SPI_MATRIX
  USE SPI_MATRIX_DEF
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
  SUBROUTINE ASSEMBLYBEGIN_MATRIX(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT) :: self
     ! LOCAL

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

  END SUBROUTINE FREE_MATRIX
  ! ..........................................................  
END MODULE SPI_MATRIX 
