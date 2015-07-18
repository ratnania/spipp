!# -*- coding: utf8 -*-
MODULE MATRIX_DEF
  USE SPM_DEF
  USE JOREK_GLOB_DEF
  USE tracelog_module
  USE BLACKBOX_DEF
  USE GREENBOX_DEF
  USE FIELD_DEF
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: INSERT_LIST_MATRICES_2D, &
          & DELETE_LIST_MATRICES_2D, &
          & INSERT_LIST_MATRICES_2D1D, &
          & DELETE_LIST_MATRICES_2D1D

  PUBLIC :: INSERT_LIST_MATRIX_FIELDS_2D  , &
          & DELETE_LIST_MATRIX_FIELDS_2D  , &
          & LENGTH_LIST_MATRIX_FIELDS_2D  , &
          & INSERT_LIST_MATRIX_FIELDS_2D1D, &
          & DELETE_LIST_MATRIX_FIELDS_2D1D, &
          & LENGTH_LIST_MATRIX_FIELDS_2D1D

  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 2 

  ! ..................................................
  TYPE, PUBLIC, ABSTRACT :: DEF_MATRIX_ABSTRACT
     ! ... matrix ID in SPM
     INTEGER :: oi_matrix_id 
     INTEGER :: oi_nvar
     LOGICAL :: ol_to_assembly

     REAL(KIND=RK), DIMENSION(:)  , ALLOCATABLE  :: opr_global_rhs 
     REAL(KIND=RK), DIMENSION(:)  , ALLOCATABLE  :: opr_global_unknown

     REAL(KIND=RK), DIMENSION(:,:), POINTER :: Matrix_Contribution 
     REAL(KIND=RK), DIMENSION(:)  , POINTER :: Rhs_Contribution
  END TYPE DEF_MATRIX_ABSTRACT
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_MATRIX_ABSTRACT) :: DEF_MATRIX_2D
     TYPE(DEF_FIELD_2D_CELL)            ::  oo_field_unknowns

     PROCEDURE(matrix_weak_formulation_2D), POINTER :: ptr_matrix_contribution  => NULL ()
     PROCEDURE(rhs_weak_formulation_2D)   , POINTER :: ptr_rhs_contribution     => NULL ()
     PROCEDURE(Assembly_Diagnostics_2D)   , POINTER :: ptr_assembly_diagnostics => NULL ()
  END TYPE DEF_MATRIX_2D
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_MATRIX_ABSTRACT) :: DEF_MATRIX_2D1D
     TYPE(DEF_FIELD_2D1D_CELL)     ::  oo_field_unknowns

     PROCEDURE(matrix_weak_formulation_2D1D), POINTER :: ptr_matrix_contribution  => NULL ()
     PROCEDURE(rhs_weak_formulation_2D1D)   , POINTER :: ptr_rhs_contribution     => NULL ()
     PROCEDURE(Assembly_Diagnostics_2D1D)   , POINTER :: ptr_assembly_diagnostics => NULL ()
  END TYPE DEF_MATRIX_2D1D
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC :: DEF_MATRIX_2D_CELL
     TYPE(DEF_MATRIX_2D), POINTER  :: ptr_matrix
     TYPE(DEF_MATRIX_2D_CELL), POINTER :: Next=>NULL()
  END TYPE DEF_MATRIX_2D_CELL
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC :: DEF_MATRIX_2D1D_CELL
     TYPE(DEF_MATRIX_2D1D), POINTER  :: ptr_matrix
     TYPE(DEF_MATRIX_2D1D_CELL), POINTER :: Next=>NULL()
  END TYPE DEF_MATRIX_2D1D_CELL
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC :: DEF_MATRIX_FIELD_2D_CELL
     ! ... local field id in the matrix
     INTEGER :: oi_loc_id

     TYPE(DEF_FIELD_2D), POINTER  :: ptr_field
     TYPE(DEF_MATRIX_2D), POINTER  :: ptr_matrix
     TYPE(DEF_MATRIX_FIELD_2D_CELL), POINTER :: Next=>NULL()
  END TYPE DEF_MATRIX_FIELD_2D_CELL
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC :: DEF_MATRIX_FIELD_2D1D_CELL
     ! ... local field id in the matrix
     INTEGER :: oi_loc_id

     TYPE(DEF_FIELD_2D1D), POINTER  :: ptr_field
     TYPE(DEF_MATRIX_2D1D), POINTER  :: ptr_matrix
     TYPE(DEF_MATRIX_FIELD_2D1D_CELL), POINTER :: Next=>NULL()
  END TYPE DEF_MATRIX_FIELD_2D1D_CELL
  ! ..................................................

  ! ..................................................
  INTERFACE
     SUBROUTINE matrix_weak_formulation_2D(self, ao_BBox2Di, ao_BBox2Dj, ao_GBox2D)
       USE BLACKBOX_DEF
       USE GREENBOX_DEF
       IMPORT DEF_MATRIX_2D

       CLASS(DEF_MATRIX_2D)  :: self
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2Di
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2Dj
       TYPE(DEF_GREENBOX_2D) :: ao_GBox2D
     END SUBROUTINE matrix_weak_formulation_2D
  END INTERFACE
  ! ..................................................

  ! ..................................................
  INTERFACE
     SUBROUTINE rhs_weak_formulation_2D(self, ao_BBox2D, ao_GBox2D)
       USE BLACKBOX_DEF
       USE GREENBOX_DEF
       IMPORT DEF_MATRIX_2D

       CLASS(DEF_MATRIX_2D)  :: self
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
       TYPE(DEF_GREENBOX_2D) :: ao_GBox2D
     END SUBROUTINE rhs_weak_formulation_2D
  END INTERFACE
  ! ..................................................

  ! ..................................................
  INTERFACE
     SUBROUTINE matrix_weak_formulation_2D1D(self, ao_BBox2Di, ao_BBox1Di, ao_BBox2Dj, ao_BBox1Dj, ao_GBox2D1D)
       USE BLACKBOX_DEF
       USE GREENBOX_DEF
       IMPORT DEF_MATRIX_2D1D

       CLASS(DEF_MATRIX_2D1D)  :: self
       TYPE(DEF_BLACKBOX_2D)   :: ao_BBox2Di
       TYPE(DEF_BLACKBOX_1D)   :: ao_BBox1Di
       TYPE(DEF_BLACKBOX_2D)   :: ao_BBox2Dj
       TYPE(DEF_BLACKBOX_1D)   :: ao_BBox1Dj
       TYPE(DEF_GREENBOX_2D1D) :: ao_GBox2D1D
     END SUBROUTINE matrix_weak_formulation_2D1D
  END INTERFACE
  ! ..................................................

  ! ..................................................
  INTERFACE
     SUBROUTINE rhs_weak_formulation_2D1D(self, ao_BBox2D, ao_BBox1D, ao_GBox2D1D)
       USE BLACKBOX_DEF
       USE GREENBOX_DEF
       IMPORT DEF_MATRIX_2D1D

       CLASS(DEF_MATRIX_2D1D)  :: self
       TYPE(DEF_BLACKBOX_2D)   :: ao_BBox2D
       TYPE(DEF_BLACKBOX_1D)   :: ao_BBox1D
       TYPE(DEF_GREENBOX_2D1D) :: ao_GBox2D1D
     END SUBROUTINE rhs_weak_formulation_2D1D
  END INTERFACE
  ! ..................................................

  ! ..................................................
  INTERFACE
     SUBROUTINE Assembly_Diagnostics_2D(self, ao_BBox2D, ao_GBox2D, ai_nstep)
       USE BLACKBOX_DEF
       USE GREENBOX_DEF
       IMPORT DEF_MATRIX_2D

       CLASS(DEF_MATRIX_2D)  :: self
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
       TYPE(DEF_GREENBOX_2D) :: ao_GBox2D
       INTEGER :: ai_nstep
     END SUBROUTINE Assembly_Diagnostics_2D
  END INTERFACE
  ! ..................................................

  ! ..................................................
  INTERFACE
     SUBROUTINE Assembly_Diagnostics_2D1D(self, ao_BBox2D, ao_BBox1D, ao_GBox2D1D, ai_nstep)
       USE BLACKBOX_DEF
       USE GREENBOX_DEF
       IMPORT DEF_MATRIX_2D1D

       CLASS(DEF_MATRIX_2D1D)  :: self
       TYPE(DEF_BLACKBOX_2D)   :: ao_BBox2D
       TYPE(DEF_BLACKBOX_1D)   :: ao_BBox1D
       TYPE(DEF_GREENBOX_2D1D) :: ao_GBox2D1D
       INTEGER :: ai_nstep
     END SUBROUTINE Assembly_Diagnostics_2D1D
  END INTERFACE
  ! ..................................................

CONTAINS

  ! .............................................
  SUBROUTINE INSERT_LIST_MATRICES_2D(Head, ptr_matrix)
  IMPLICIT NONE
    TYPE(DEF_MATRIX_2D_CELL),INTENT(INOUT)       :: Head
    TYPE(DEF_MATRIX_2D), TARGET, INTENT(IN)     :: ptr_matrix
    ! LOCAL
    TYPE(DEF_MATRIX_2D_CELL), POINTER            :: TheCell=>NULL()
    TYPE(DEF_MATRIX_2D_CELL), POINTER            :: PrevCell=>NULL()
    TYPE(DEF_MATRIX_2D_CELL), POINTER            :: NewCell
    ! -------------------------------------

#ifdef DEBUG_TRACE 
     CALL printlog("INSERT_LIST_MATRICES_2D: Begin", ai_dtllevel = 0)
#endif

    TheCell  => Head%next

    IF (.NOT. ASSOCIATED(TheCell))  THEN
       ALLOCATE(NewCell); NewCell % ptr_matrix => ptr_matrix
       Head % Next => NewCell
    ELSE
       IF (TheCell % ptr_matrix % oi_matrix_id > ptr_matrix % oi_matrix_id) THEN
          ALLOCATE(NewCell); NewCell % ptr_matrix=> ptr_matrix
          NewCell % Next  => TheCell
          Head % Next => NewCell
       ELSE 
          PrevCell => Head % next
          TheCell  => PrevCell % next

          DO WHILE( ASSOCIATED(TheCell) ) 
             IF (TheCell % ptr_matrix % oi_matrix_id > ptr_matrix % oi_matrix_id) THEN
                EXIT
             END IF
             PrevCell   => PrevCell % Next 
             TheCell    => TheCell % Next 
          END DO

          IF ( PrevCell % ptr_matrix % oi_matrix_id== ptr_matrix % oi_matrix_id)  THEN
          ELSE
             ALLOCATE(NewCell); NewCell % ptr_matrix=> ptr_matrix
             NewCell % Next  => TheCell
             PrevCell % Next => NewCell 
          END IF
       END IF
    END IF

#ifdef DEBUG_TRACE 
     CALL printlog("INSERT_LIST_MATRICES_2D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE INSERT_LIST_MATRICES_2D
  ! .............................................

  ! .............................................
  SUBROUTINE DELETE_LIST_MATRICES_2D(HashTable)
  IMPLICIT NONE
    TYPE(DEF_MATRIX_2D_CELL), INTENT(INOUT) ::  HashTable
    ! LOCAL
    TYPE(DEF_MATRIX_2D_CELL), POINTER         :: TheCell=>NULL()
    TYPE(DEF_MATRIX_2D_CELL), POINTER         :: PrevCell=>NULL()
    INTEGER :: li_start
    INTEGER :: ipos, is, j

#ifdef DEBUG_TRACE 
     CALL printlog("DELETE_LIST_MATRICES_2D: Begin", ai_dtllevel = 0)
#endif

    TheCell => HashTable%Next

    DO WHILE( ASSOCIATED(TheCell) ) 
       PrevCell     => TheCell
       TheCell      => TheCell%Next
       DEALLOCATE( PrevCell )
    END DO

    HashTable % Next => NULL()

#ifdef DEBUG_TRACE 
     CALL printlog("DELETE_LIST_MATRICES_2D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE DELETE_LIST_MATRICES_2D
  ! .............................................

  ! .............................................
  SUBROUTINE INSERT_LIST_MATRICES_2D1D(Head, ptr_matrix)
  IMPLICIT NONE
    TYPE(DEF_MATRIX_2D1D_CELL),INTENT(INOUT)       :: Head
    TYPE(DEF_MATRIX_2D1D), POINTER, INTENT(IN)     :: ptr_matrix
    ! LOCAL
    TYPE(DEF_MATRIX_2D1D_CELL), POINTER            :: TheCell=>NULL()
    TYPE(DEF_MATRIX_2D1D_CELL), POINTER            :: PrevCell=>NULL()
    TYPE(DEF_MATRIX_2D1D_CELL), POINTER            :: NewCell
    ! -------------------------------------

#ifdef DEBUG_TRACE 
     CALL printlog("INSERT_LIST_MATRICES_2D: Begin", ai_dtllevel = 0)
#endif

    TheCell  => Head%next

    IF (.NOT. ASSOCIATED(TheCell))  THEN
       ALLOCATE(NewCell); NewCell % ptr_matrix => ptr_matrix
       Head % Next => NewCell
    ELSE
       IF (TheCell % ptr_matrix % oi_matrix_id > ptr_matrix % oi_matrix_id) THEN
          ALLOCATE(NewCell); NewCell % ptr_matrix=> ptr_matrix
          NewCell % Next  => TheCell
          Head % Next => NewCell
       ELSE 
          PrevCell => Head % next
          TheCell  => PrevCell % next

          DO WHILE( ASSOCIATED(TheCell) ) 
             IF (TheCell % ptr_matrix % oi_matrix_id > ptr_matrix % oi_matrix_id) THEN
                EXIT
             END IF
             PrevCell   => PrevCell % Next 
             TheCell    => TheCell % Next 
          END DO

          IF ( PrevCell % ptr_matrix % oi_matrix_id== ptr_matrix % oi_matrix_id)  THEN
          ELSE
             ALLOCATE(NewCell); NewCell % ptr_matrix=> ptr_matrix
             NewCell % Next  => TheCell
             PrevCell % Next => NewCell 
          END IF
       END IF
    END IF

#ifdef DEBUG_TRACE 
     CALL printlog("INSERT_LIST_MATRICES_2D1D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE INSERT_LIST_MATRICES_2D1D
  ! .............................................

  ! .............................................
  SUBROUTINE DELETE_LIST_MATRICES_2D1D(HashTable)
  IMPLICIT NONE
    TYPE(DEF_MATRIX_2D1D_CELL), INTENT(INOUT) ::  HashTable
    ! LOCAL
    TYPE(DEF_MATRIX_2D1D_CELL), POINTER         :: TheCell=>NULL()
    TYPE(DEF_MATRIX_2D1D_CELL), POINTER         :: PrevCell=>NULL()
    INTEGER :: li_start
    INTEGER :: ipos, is, j

#ifdef DEBUG_TRACE 
     CALL printlog("DELETE_LIST_MATRICES_2D1D: Begin", ai_dtllevel = 0)
#endif

    TheCell => HashTable%Next

    DO WHILE( ASSOCIATED(TheCell) ) 
       PrevCell     => TheCell
       TheCell      => TheCell%Next
       DEALLOCATE( PrevCell )
    END DO

    HashTable % Next => NULL()

#ifdef DEBUG_TRACE 
     CALL printlog("DELETE_LIST_MATRICES_2D1D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE DELETE_LIST_MATRICES_2D1D
  ! .............................................

  ! .............................................
  SUBROUTINE INSERT_LIST_MATRIX_FIELDS_2D(Head, ao_field, ao_matrix)
  IMPLICIT NONE
    TYPE(DEF_MATRIX_FIELD_2D_CELL),INTENT(INOUT)    :: Head
    TYPE(DEF_FIELD_2D) , TARGET   , INTENT(IN)      :: ao_field
    TYPE(DEF_MATRIX_2D), TARGET   , INTENT(IN)      :: ao_matrix
    ! LOCAL
    TYPE(DEF_MATRIX_FIELD_2D_CELL), POINTER            :: TheCell=>NULL()
    TYPE(DEF_MATRIX_FIELD_2D_CELL), POINTER            :: PrevCell=>NULL()
    TYPE(DEF_MATRIX_FIELD_2D_CELL), POINTER            :: field_cell=>NULL()
    TYPE(DEF_MATRIX_FIELD_2D_CELL), POINTER            :: NewCell
    LOGICAL :: ll_field_found
    ! -------------------------------------

#ifdef DEBUG_TRACE 
     CALL printlog("INSERT_LIST_MATRIX_FIELDS_2D: Begin", ai_dtllevel = 0)
#endif

    ll_field_found = .FALSE.

    field_cell => Head % Next
    DO WHILE( ASSOCIATED(field_cell) )
       IF (field_cell % ptr_field % oi_id == ao_field % oi_id) THEN
          ll_field_found = .TRUE.
       END IF

       field_cell => field_cell % Next
    END DO

    IF (.NOT. ll_field_found) THEN
       TheCell  => Head%next
      
       IF( .NOT. ASSOCIATED(TheCell) )  THEN
          ALLOCATE(NewCell); NewCell%ptr_field => ao_field; NewCell%ptr_matrix => ao_matrix
          Head%Next => NewCell
       ELSE
          ALLOCATE(NewCell); NewCell%ptr_field => ao_field; NewCell%ptr_matrix => ao_matrix
          NewCell%Next  => TheCell
          Head%Next => NewCell
       END IF

       NewCell % oi_loc_id = FIND_LIST_FIELDS_2D(ao_matrix % oo_field_unknowns, ao_field)

       IF (NewCell % oi_loc_id < 0) THEN
          PRINT *, "INSERT_LIST_MATRIX_FIELDS_2D: could not find field in the given matrix"
          STOP
       END IF

       !PRINT *,  "field id ", ao_field % oi_id, &
            !   & " matrix id ", ao_matrix % oi_matrix_id, &
             !  & " local id  ", NewCell % oi_loc_id

    END IF

#ifdef DEBUG_TRACE 
     CALL printlog("INSERT_LIST_MATRIX_FIELDS_2D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE INSERT_LIST_MATRIX_FIELDS_2D
  ! .............................................

  ! .............................................
  SUBROUTINE DELETE_LIST_MATRIX_FIELDS_2D(HashTable)
  IMPLICIT NONE
    TYPE(DEF_MATRIX_FIELD_2D_CELL), INTENT(INOUT) ::  HashTable
    ! LOCAL
    TYPE(DEF_MATRIX_FIELD_2D_CELL), POINTER         :: TheCell=>NULL()
    TYPE(DEF_MATRIX_FIELD_2D_CELL), POINTER         :: PrevCell=>NULL()
    INTEGER :: li_start
    INTEGER :: ipos, is, j

#ifdef DEBUG_TRACE 
     CALL printlog("DELETE_LIST_MATRIX_FIELDS_2D: Begin", ai_dtllevel = 0)
#endif

    TheCell => HashTable%Next

    DO WHILE( ASSOCIATED(TheCell) ) 
       PrevCell     => TheCell
       TheCell      => TheCell%Next
       DEALLOCATE( PrevCell )
    END DO

    HashTable % Next => NULL()

#ifdef DEBUG_TRACE 
     CALL printlog("DELETE_LIST_MATRIX_FIELDS_2D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE DELETE_LIST_MATRIX_FIELDS_2D
  ! .............................................

  ! .............................................
  FUNCTION  LENGTH_LIST_MATRIX_FIELDS_2D(HashTable) RESULT(Nedge)
    TYPE(DEF_MATRIX_FIELD_2D_CELL)    ::  HashTable
    TYPE(DEF_MATRIX_FIELD_2D_CELL), POINTER         :: TheCell=>NULL()
    INTEGER                     :: Nedge, is
    ! -------------------------------------

#ifdef DEBUG_TRACE 
     CALL printlog("LENGTH_LIST_MATRIX_FIELDS_2D: Begin", ai_dtllevel = 0)
#endif

    Nedge=0
    TheCell => HashTable%Next
    DO WHILE( ASSOCIATED(TheCell) ) 
       Nedge       = Nedge +1
       TheCell      => TheCell%Next 
    END DO

#ifdef DEBUG_TRACE 
     CALL printlog("LENGTH_LIST_MATRIX_FIELDS_2D: End", ai_dtllevel = 0)
#endif

  END FUNCTION LENGTH_LIST_MATRIX_FIELDS_2D
  ! .............................................

  ! .............................................
  SUBROUTINE INSERT_LIST_MATRIX_FIELDS_2D1D(Head, ao_field, ao_matrix)
  IMPLICIT NONE
    TYPE(DEF_MATRIX_FIELD_2D1D_CELL),INTENT(INOUT)    :: Head
    TYPE(DEF_FIELD_2D1D) , TARGET   , INTENT(IN)      :: ao_field
    TYPE(DEF_MATRIX_2D1D), TARGET   , INTENT(IN)      :: ao_matrix
    ! LOCAL
    TYPE(DEF_MATRIX_FIELD_2D1D_CELL), POINTER            :: TheCell=>NULL()
    TYPE(DEF_MATRIX_FIELD_2D1D_CELL), POINTER            :: PrevCell=>NULL()
    TYPE(DEF_MATRIX_FIELD_2D1D_CELL), POINTER            :: NewCell
    TYPE(DEF_MATRIX_FIELD_2D1D_CELL), POINTER            :: field_cell=>NULL()
    LOGICAL :: ll_field_found
    ! -------------------------------------

#ifdef DEBUG_TRACE 
     CALL printlog("INSERT_LIST_MATRIX_FIELDS_2D1D: Begin", ai_dtllevel = 0)
#endif

    ll_field_found = .FALSE.

    field_cell => Head % Next
    DO WHILE( ASSOCIATED(field_cell) )
       IF (field_cell % ptr_field % oi_id == ao_field % oi_id) THEN
          ll_field_found = .TRUE.
       END IF

       field_cell => field_cell % Next
    END DO

    IF (.NOT. ll_field_found) THEN
       TheCell  => Head%next
      
       IF( .NOT. ASSOCIATED(TheCell) )  THEN
          ALLOCATE(NewCell); NewCell%ptr_field => ao_field; NewCell%ptr_matrix => ao_matrix
          Head%Next => NewCell
       ELSE
          ALLOCATE(NewCell); NewCell%ptr_field => ao_field; NewCell%ptr_matrix => ao_matrix
          NewCell%Next  => TheCell
          Head%Next => NewCell
       END IF

       NewCell % oi_loc_id = FIND_LIST_FIELDS_2D1D(ao_matrix % oo_field_unknowns, ao_field)
       IF (NewCell % oi_loc_id < 0) THEN
          PRINT *, "INSERT_LIST_MATRIX_FIELDS_2D1D: could not find field in the given matrix"
          STOP
       END IF

       PRINT *,  "field id ", ao_field % oi_id, &
               & " matrix id ", ao_matrix % oi_matrix_id, &
               & " local id  ", NewCell % oi_loc_id

    END IF

#ifdef DEBUG_TRACE 
     CALL printlog("INSERT_LIST_MATRIX_FIELDS_2D1D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE INSERT_LIST_MATRIX_FIELDS_2D1D
  ! .............................................

  ! .............................................
  SUBROUTINE DELETE_LIST_MATRIX_FIELDS_2D1D(HashTable)
  IMPLICIT NONE
    TYPE(DEF_MATRIX_FIELD_2D1D_CELL), INTENT(INOUT) ::  HashTable
    ! LOCAL
    TYPE(DEF_MATRIX_FIELD_2D1D_CELL), POINTER         :: TheCell=>NULL()
    TYPE(DEF_MATRIX_FIELD_2D1D_CELL), POINTER         :: PrevCell=>NULL()
    INTEGER :: li_start
    INTEGER :: ipos, is, j

#ifdef DEBUG_TRACE 
     CALL printlog("DELETE_LIST_MATRIX_FIELDS_2D1D: Begin", ai_dtllevel = 0)
#endif

    TheCell => HashTable%Next

    DO WHILE( ASSOCIATED(TheCell) ) 
       PrevCell     => TheCell
       TheCell      => TheCell%Next
       DEALLOCATE( PrevCell )
    END DO

    HashTable % Next => NULL()

#ifdef DEBUG_TRACE 
     CALL printlog("DELETE_LIST_MATRIX_FIELDS_2D1D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE DELETE_LIST_MATRIX_FIELDS_2D1D
  ! .............................................

  ! .............................................
  FUNCTION  LENGTH_LIST_MATRIX_FIELDS_2D1D(HashTable) RESULT(Nedge)
    TYPE(DEF_MATRIX_FIELD_2D1D_CELL)    ::  HashTable
    TYPE(DEF_MATRIX_FIELD_2D1D_CELL), POINTER         :: TheCell=>NULL()
    INTEGER                     :: Nedge, is
    ! -------------------------------------

#ifdef DEBUG_TRACE 
     CALL printlog("LENGTH_LIST_MATRIX_FIELDS_2D1D: Begin", ai_dtllevel = 0)
#endif

    Nedge=0
    TheCell => HashTable%Next
    DO WHILE( ASSOCIATED(TheCell) ) 
       Nedge       = Nedge +1
       TheCell      => TheCell%Next 
    END DO

#ifdef DEBUG_TRACE 
     CALL printlog("LENGTH_LIST_MATRIX_FIELDS_2D1D: End", ai_dtllevel = 0)
#endif

  END FUNCTION LENGTH_LIST_MATRIX_FIELDS_2D1D
  ! .............................................

END MODULE MATRIX_DEF
