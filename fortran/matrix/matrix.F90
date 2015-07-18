!# -*- coding: utf8 -*-
MODULE MATRIX
  USE TypeDef        
  USE SPM_DEF
  USE SPM
  USE DIRICHLET_MOD
  USE JOREK_PARAM
  USE NUMBERING
  USE JOREK_GLOB
  USE GRAPH
  USE MATRIX_DEF
  USE FIELD_DEF
  USE FIELD
  USE MATRIX_2D
  USE MATRIX_2D1D
  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 2 

  ! .............................................
  
  ! .............................................
  PRIVATE
  PUBLIC :: MATRIX_NUMBER_UNKNOWN_FIELD & 
       ,    MATRIX_CREATE               &
       ,    GET_NR_MATRIX               &
       ,    GET_NC_MATRIX               &
       ,    MATRIX_APPEND_UNKNOWN_FIELD &
       ,    INITIALIZE_MATRIX           &
       ,    MATRIX_RESET_ELEMENT_MATRIX &
       ,    MATRIX_RESET_ELEMENT_RHS    &
       ,    MATRIX_GET_GLOBAL_ROW_INDEX &
       ,    MATRIX_FREE
  

  interface INITIALIZE_MATRIX
     module procedure  INITIALIZE_MATRIX_2D, INITIALIZE_MATRIX_2D1D
  end interface

CONTAINS
  ! ..................................................

  ! ..................................................
  SUBROUTINE MATRIX_APPEND_UNKNOWN_FIELD(self, ao_field)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT), INTENT(INOUT) :: self 
     CLASS(DEF_FIELD_ABSTRACT), INTENT(IN) :: ao_field
     ! LOCAL
    
#ifdef DEBUG_TRACE 
    CALL printlog("MATRIX_APPEND_UNKNOWN_FIELD: Begin", ai_dtllevel = 0)
#endif

     PRINT *, "INSERT FIELD ", ao_field % oi_id, " into the matrix ", self % oi_matrix_id

     SELECT TYPE (self)
     ! ... 2D case 
     CLASS IS (DEF_MATRIX_2D)
        SELECT TYPE (ao_field)
        ! ... 
        CLASS IS (DEF_FIELD_2D)
           CALL MATRIX_APPEND_UNKNOWN_FIELD_2D(self, ao_field)
        ! ...
        CLASS DEFAULT
           STOP 'MATRIX_APPEND_UNKNOWN_FIELD: unexpected type for ao_field object!'
        END SELECT
     ! ...
     ! ... 2D1D case 
     CLASS IS (DEF_MATRIX_2D1D)
        SELECT TYPE (ao_field)
        ! ... 2D matrix 
        CLASS IS (DEF_FIELD_2D1D)
           CALL MATRIX_APPEND_UNKNOWN_FIELD_2D1D(self, ao_field)
        ! ...
        CLASS DEFAULT
           STOP 'MATRIX_APPEND_UNKNOWN_FIELD: unexpected type for ao_field object!'
        END SELECT
     ! ...
     CLASS DEFAULT
        STOP 'MATRIX_APPEND_UNKNOWN_FIELD: unexpected type for self object!'
     END SELECT

#ifdef DEBUG_TRACE 
     CALL printlog("MATRIX_APPEND_UNKNOWN_FIELD: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE MATRIX_APPEND_UNKNOWN_FIELD
  ! ..................................................
  FUNCTION MATRIX_NUMBER_UNKNOWN_FIELD(self) RESULT (Nedge)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT), INTENT(IN) :: self 
     INTEGER :: Nedge
     ! LOCAL  
#ifdef DEBUG_TRACE 
     CALL printlog("MATRIX_NUMBER_UNKNOWN_FIELD: Begin", ai_dtllevel = 0)
#endif 
     SELECT TYPE (self)
     !2D case
     CLASS IS (DEF_MATRIX_2D)
        Nedge = MATRIX_NUMBER_UNKNOWN_FIELD_2D(self)
     !2D1D case
     CLASS IS (DEF_MATRIX_2D1D)
        Nedge = MATRIX_NUMBER_UNKNOWN_FIELD_2D1D(self)
     !Default
     CLASS DEFAULT
        STOP 'MATRIX_NUMBER_UNKNOWN_FIELD : unexpected type for self object!'
     END SELECT
#ifdef DEBUG_TRACE
     CALL printlog("MATRIX_NUMBER_UNKNOWN_FIELD: End", ai_dtllevel = 0)  
#endif
   END FUNCTION MATRIX_NUMBER_UNKNOWN_FIELD
  ! ..................................................
  SUBROUTINE MATRIX_CREATE(self)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT), INTENT(INOUT) :: self 
     ! LOCAL
     INTEGER :: li_matrix_id
     INTEGER :: li_err
    
#ifdef DEBUG_TRACE 
     CALL printlog("MATRIX_CREATE: Begin", ai_dtllevel = 0)
#endif

     CALL SPM_CREATE_MATRIX(li_matrix_id, li_err)

     self % oi_matrix_id   = li_matrix_id
     self % ol_to_assembly = .FALSE.

#ifdef DEBUG_TRACE 
     CALL printlog("MATRIX_CREATE: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE MATRIX_CREATE
  ! ..................................................

  ! ..................................................
  SUBROUTINE GET_NR_MATRIX(self, ai_value)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT), INTENT(INOUT) :: self 
     INTEGER, INTENT(INOUT)           :: ai_value
     ! LOCAL
     INTEGER :: li_err
    
#ifdef DEBUG_TRACE 
     CALL printlog("GET_NR_MATRIX: Begin", ai_dtllevel = 0)
#endif

     CALL SPM_GetnR(self % oi_matrix_id, ai_value, li_err)

#ifdef DEBUG_TRACE 
     CALL printlog("GET_NR_MATRIX: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE GET_NR_MATRIX
  ! ..................................................

  ! ..................................................
  SUBROUTINE GET_NC_MATRIX(self, ai_value)
  IMPLICIT NONE
     CLASS(DEF_MATRIX_ABSTRACT), INTENT(INOUT) :: self 
     INTEGER, INTENT(INOUT)           :: ai_value
     ! LOCAL
     INTEGER :: li_err
    
#ifdef DEBUG_TRACE 
     CALL printlog("GET_NC_MATRIX: Begin", ai_dtllevel = 0)
#endif

     CALL SPM_GetnC(self % oi_matrix_id, ai_value, li_err)

#ifdef DEBUG_TRACE 
     CALL printlog("GET_NC_MATRIX: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE GET_NC_MATRIX
  ! ..................................................

  ! ..................................................
  ! ..................................................
  SUBROUTINE MATRIX_GET_GLOBAL_ROW_INDEX(self, ai_glob, ai_real)
    IMPLICIT NONE
    CLASS(DEF_MATRIX_ABSTRACT) :: self
    INTEGER                   :: ai_glob
    INTEGER                   :: ai_real
    !LOCAL
    INTEGER                   :: li_err

#ifdef DEBUG_TRACE
    CALL printlog("MATRIX_GET_GLOBAL_ROW_INDEX: Begin", ai_dtllevel = 1)
#endif

    CALL SPM_GLOBAL_ROW_INDEX(self % oi_matrix_id, ai_glob, ai_real, li_err)
    
#ifdef DEBUG_TRACE
    CALL printlog("MATRIX_GET_GLOBAL_ROW_INDEX: Begin", ai_dtllevel = 1)
#endif


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
END MODULE MATRIX 
