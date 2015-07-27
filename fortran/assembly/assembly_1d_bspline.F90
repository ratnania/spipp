!# -*- coding: utf8 -*-
MODULE SPI_ASSEMBLY_1D_BSPLINE 
  USE SPI_GLOBAL_DEF
  USE SPI_ASSEMBLY_DEF
  USE SPI_GREENBOX_DEF
  USE SPI_GREENBOX
  USE SPI_MESH_DEF
  USE SPI_MESH
  USE SPI_MATRIX_DEF
  USE SPI_MATRIX
  USE SPI_SPACE_DEF
  USE SPI_SPACE

  IMPLICIT NONE

CONTAINS
! ----------------------------------------------------------------------
!  SUBROUTINE LOOP_ON_ELMTS_2D_ISOPARAMETRIC_BEZIER(ao_matrix, &
!                  & ao_hashtab_matrix, ao_hashtab_matrix_field, ai_color, &
!                  & ao_mesh, ao_gbox, &
!                  & ao_Basis2Di, ao_Basis2Dj)
  SUBROUTINE LOOP_ON_ELMTS_1D_BSPLINE(ao_matrix, &
                  & ao_mesh, &
                  & ao_trial_space, & 
                  & ao_test_space, & 
                  & ao_gbox)
    IMPLICIT NONE
    TYPE(DEF_MATRIX_1D), TARGET, INTENT(INOUT) :: ao_matrix
    TYPE(DEF_MESH_1D_BSPLINE), INTENT(INOUT)   :: ao_mesh 
    TYPE(DEF_SPACE_1D_BSPLINE), TARGET :: ao_trial_space
    TYPE(DEF_SPACE_1D_BSPLINE), TARGET :: ao_test_space
    TYPE(DEF_GREENBOX_1D)                      :: ao_gbox
    ! LOCAL
    LOGICAL                   :: to_assemble_matrix
    LOGICAL                   :: to_assemble_rhs
    INTEGER                   :: li_elmt_id, is, is1, is2, ipol
    INTEGER                   :: ivar, ipos
    INTEGER                   :: ipol1, iv_pol1, or_pol1
    INTEGER                   :: ipol2, iv_pol2, or_pol2
    INTEGER                   :: li_loc_1, li_loc_2
    INTEGER                   :: ierr
    INTEGER                   :: i_glob
    INTEGER                   :: i_real
    INTEGER                   :: li_nodenbr
    INTEGER                   :: li_ndof
    INTEGER                   :: li_locsize
    INTEGER                   :: li_root
    INTEGER                   :: li_i 
    INTEGER                   :: li_var
    INTEGER                   :: nivar
    INTEGER                   :: njvar
    INTEGER                   :: nivar_V
    INTEGER                   :: li_mesh

    REAL(KIND=SPI_RK)         :: lr_value
    REAL(KIND=SPI_RK)         :: lr_a
    REAL(KIND=SPI_RK)         :: lr_b
    REAL(KIND=SPI_RK)         :: lr_delta

    REAL(KIND=SPI_RK), DIMENSION(2)           :: Xis
    REAL(KIND=SPI_RK), DIMENSION(ao_mesh % oi_n_dim) :: lpr_control_point 

    INTEGER, PARAMETER :: N_DIM = 2
    INTEGER       :: ijg
    INTEGER       :: li_ctrl_pt
    INTEGER, PARAMETER :: li_one = 1
    REAL(SPI_RK), DIMENSION(ao_trial_space % oo_quad % oi_n_points) :: lpr_points
    ! ----------------------------------------------------------------------
    ! ---------     END of Local Declarations  and Specifications ----------
    ! ----------------------------------------------------------------------

    to_assemble_matrix = .TRUE.
    IF ( ASSOCIATED(ao_matrix % ptr_matrix_contribution) .EQV. .FALSE.) THEN
       to_assemble_matrix = .FALSE.
    END IF

    to_assemble_rhs = .TRUE.
    IF ( ASSOCIATED(ao_matrix % ptr_rhs_contribution) .EQV. .FALSE.) THEN
       to_assemble_rhs = .FALSE.
    END IF

    nivar = ao_matrix % oi_nvar 
    njvar = nivar 
    
    ! ... 
    li_locsize = ao_mesh % oi_nen
    ! ...

    ml_assembly_initialized = .FALSE.
    li_ndof = ao_gbox % oi_nvar

    IF (.NOT. ml_assembly_initialized) THEN
       ALLOCATE(mpi_ltog_rhsrows(li_locsize*li_ndof))
       ALLOCATE(mpi_ltog_rows(li_locsize*li_ndof,li_locsize*li_ndof))
       ALLOCATE(mpi_ltog_cols(li_locsize*li_ndof,li_locsize*li_ndof))
       ml_assembly_initialized = .TRUE.
       mpi_ltog_rhsrows  =  SPI_INT_DEFAULT
       mpi_ltog_rows  =  SPI_INT_DEFAULT
       mpi_ltog_cols  =  SPI_INT_DEFAULT
    END IF

    ! ----------------------------------------------------
    ! ----------------------------------------------------
    ! SPM Initialization
    ! ----------------------------------------------------
    IF ( to_assemble_matrix ) THEN
       CALL RESET_MATRIX_MATRIX(ao_matrix)
    END IF

    IF ( to_assemble_rhs ) THEN
       CALL RESET_MATRIX_RHS(ao_matrix)
    END IF

    CALL ASSEMBLYBEGIN_MATRIX(ao_matrix)
    ! ----------------------------------------------------

    DO li_elmt_id = 1, ao_mesh % oi_n_elmts 

       ! set to zero
       ! the position is always computed by the ^ith test function and not the ^jth 
       CALL BLACKBOX_RESET_POSITION(ao_trial_space % oo_bbox) 
       CALL BLACKBOX_RESET_POSITION(ao_test_space % oo_bbox) 

       CALL RESET_GREENBOX(ao_gbox)

       CALL RESET_ELEMENT_MATRIX(ao_matrix)
       CALL RESET_ELEMENT_RHS(ao_matrix)

       CALL RESET_BASIS(ao_trial_space % oo_basis) 
       CALL RESET_BASIS(ao_test_space % oo_basis) 

       ! ------------------------------------------ Loops on i ddl
       ! Loop on poloidal vertices (i)
       ! -------------------------------------------------------------------
       ! TODO
       lr_a = 0.0 
       lr_b = 1.0
       lr_delta = lr_b - lr_a
       CALL UPDATE_LOGICAL_POSITION_BLACKBOX_1D_BSPLINE(ao_trial_space % oo_bbox, lr_a, lr_b)
       DO li_loc_1 = 1, ao_mesh % oi_nen
          ! TODO
          lpr_control_point = 0.0
          CALL UPDATE_POSITION_BLACKBOX_1D_BSPLINE(ao_trial_space % oo_bbox, lpr_control_point) 

          DO li_var = 1, ao_gbox % oi_nvar 
             ! TODO
             lr_value = 1.0
             CALL UPDATE_VARIABLES_GREENBOX_1D(ao_gbox, ao_trial_space % oo_bbox, li_var, lr_value, li_elmt_id) 
          END DO
       END DO
       ! End of Loops 
       ! -------------------------------------------------------------------
      
       CALL Reset_Local_Data()

       ! TODO
       lpr_points = 0.0
       CALL UPDATE_BASIS(ao_trial_space % oo_basis, lpr_points) 
       CALL UPDATE_BASIS(ao_test_space  % oo_basis, lpr_points) 

       CALL COMPUTE_METRIC_BLACKBOX_1D_BSPLINE(ao_trial_space % oo_bbox, lr_a, lr_b)
       CALL COMPUTE_METRIC_BLACKBOX_1D_BSPLINE(ao_test_space  % oo_bbox, lr_a, lr_b)

       CALL BLACKBOX_UPDATE_PHYSICAL_BASIS(ao_trial_space % oo_bbox) 

       ! ... TODO: add here coordinates system update (not sure be needed)
       ! ...

       ! ... TODO: add here pushback variables update       
       ! ...


       ! ... loop over non vanishing basis functions
       DO li_loc_1 = 1, ao_mesh % oi_nen
          ! ... get the local basis
          CALL GET_BASIS_BLACKBOX_1D_BSPLINE(ao_trial_space % oo_bbox, li_loc_1) 
          ! ...

          ! ... update basis in the physical domain
          CALL UPDATE_PHYSICAL_BASIS_BLACKBOX(ao_trial_space % oo_bbox) 
          ! ...

          ! ............................................. 
          !     here RHS bloc-Vector  to be completed 
          ! ............................................. 
          IF ( to_assemble_rhs ) THEN
             ! ... Call the weak formulation for the RHS 
             CALL ao_matrix % ptr_rhs_contribution(ao_gbox)
             ! ...

             IF (ISNAN(ao_matrix % Rhs_Contribution(1))) THEN
                ao_matrix % Rhs_Contribution = 1.0E14 
             END IF
             
!             CALL SPM_ASSEMBLYRHSADDVALUES_LOCAL(ao_matrix % oi_matrix_id, li_loc_1, &
!                     & ao_matrix % Rhs_Contribution, ierr)
             
             ! ... reset RHS contribution array
             CALL MATRIX_RESET_ELEMENT_RHS(ao_matrix)
             ! ...
          END IF
          ! ............................................. 

          IF ( to_assemble_matrix ) THEN
             DO li_loc_2 = 1, ao_mesh % oi_nen
                ! ... get the local basis
                CALL GET_BASIS_BLACKBOX_1D_BSPLINE(ao_trial_space % oo_bbox, li_loc_2) 
                ! ...

                ! ... update basis in the physical domain
                CALL UPDATE_PHYSICAL_BASIS_BLACKBOX(ao_test_space % oo_bbox) 
                ! ...

                ! ... Call the weak formulation for the matrix
                CALL ao_matrix % ptr_matrix_contribution(ao_gbox)
                ! ...
                
                IF (ISNAN(ao_matrix % Matrix_Contribution(1,1))) THEN
                   ao_matrix % Matrix_Contribution = 1.0E14
                END IF
!                
!                CALL SPM_ASSEMBLYADDVALUES_LOCAL(ao_matrix % oi_matrix_id, li_loc_1, li_loc_2, &
!                        & ao_matrix % Matrix_Contribution, ierr)
                
                ! ... reset Matrix contribution array
                CALL MATRIX_RESET_ELEMENT_MATRIX(ao_matrix)
                ! ...
                END DO
             END IF

       END DO
       ! -------------------------------------------------------------------
     
       ! -- UPDATE RHS AND MATRICES FROM LOCAL-RHS
       nivar = ao_matrix % oi_nvar 
       njvar = nivar 

!       CALL GET_LTOG(ao_mesh, ie_Pol, &
!               & mpi_ltog_rows, mpi_ltog_cols, mpi_ltog_rhsrows &
!               , nivar, njvar)
!
!       IF ( to_assemble_matrix ) THEN
!          CALL SPM_ASSEMBLYPUSH(ao_matrix % oi_matrix_id, mpi_ltog_rows, mpi_ltog_cols, ierr)
!       END IF
!       IF ( to_assemble_rhs ) THEN
!          CALL SPM_ASSEMBLYRHSPUSH(ao_matrix % oi_matrix_id, mpi_ltog_rhsrows, ierr)
!       END IF

       mpi_ltog_rhsrows  =  SPI_INT_DEFAULT
       mpi_ltog_rows  =  SPI_INT_DEFAULT
       mpi_ltog_cols  =  SPI_INT_DEFAULT
       ! -- 

    END DO 
    ! -------------------------------------------------------------------
    ! -------------------------------------------------------------------

    ! ----------------------------------------------------
    ! End of the Assembly process
    ! ----------------------------------------------------
    CALL ASSEMBLYEND_MATRIX(ao_matrix)
    ! ----------------------------------------------------

    IF (ml_assembly_initialized) THEN
       DEALLOCATE(mpi_ltog_rhsrows)
       DEALLOCATE(mpi_ltog_rows)
       DEALLOCATE(mpi_ltog_cols)
       ml_assembly_initialized = .FALSE.
    END IF

  CONTAINS

    ! -------------------------------------------------------------------
    SUBROUTINE Reset_Local_Data()
      IMPLICIT NONE
      IF ( to_assemble_matrix ) THEN
         CALL RESET_ELEMENT_MATRIX(ao_matrix)
      END IF
      IF ( to_assemble_rhs ) THEN
         CALL RESET_ELEMENT_RHS_MATRIX(ao_matrix)
      END IF
    END SUBROUTINE Reset_Local_Data
    ! -------------------------------------------------------------------

  END SUBROUTINE LOOP_ON_ELMTS_1D_BSPLINE
    ! -------------------------------------------------------------------
    ! -------------------------------------------------------------------
END MODULE SPI_ASSEMBLY_1D_BSPLINE
