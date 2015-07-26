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
    INTEGER                   :: i_loc1, i_loc2
    INTEGER                   :: ierr
    INTEGER                   :: i_glob
    INTEGER                   :: i_real
    INTEGER                   :: li_nodenbr
    INTEGER                   :: li_ndof
    INTEGER                   :: li_locsize
    INTEGER                   :: li_root
    INTEGER                   :: li_i 
    INTEGER                   :: nivar
    INTEGER                   :: njvar
    INTEGER                   :: ivar_V
    INTEGER                   :: nivar_V
    INTEGER                   :: li_mesh

    REAL(KIND=SPI_RK) :: lr_value

    REAL(KIND=SPI_RK)   , DIMENSION(2)           :: Xis

    INTEGER, PARAMETER :: N_DIM = 2
    INTEGER       :: ijg
    INTEGER       :: li_ctrl_pt
    INTEGER, PARAMETER :: li_one = 1
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
!       DO i_loc1 = 1, lp_elmt % oi_nen
!          Scale2D = 1.0 ! lp_elmt % Scales(or_pol1, iv_pol1)
!          li_ctrl_pt = lp_elmt % opi_vertices_sons(i_loc1)
!
!          CALL getPoloidalBasisFunction(i_loc1, li_one, &
!                                       & Scale2D, &
!                                       & lo_BBox2Di)
!
!          ! R and Z coordinates and derivatives at Gauss points
!          ! ---------------------------------------------------- BEGIN
!          Xis = ao_mesh % TheNodes_bezier % Coor2D(1:N_dim, li_ctrl_pt, li_one)
!!          print *, "i_loc1, li_ctrl_pt",i_loc1, li_ctrl_pt
!!          print *, "control point ", li_ctrl_pt, Xis
!          CALL updatePoloidalPosition(i_loc1, li_one, & 
!               & Scale2D, Xis, &
!               & lp_elmt % coo_transf_A(:,:), &
!               & lp_elmt % coo_transf_b(:), &
!               & lo_BBox2Di)
!          ! ---------------------------------------------------- END
!
!          ! -------------------------------------------------------------------
!          ipos   = lp_elmt % opi_LocToGlob(i_loc1)
!          matrix_field_cell => ao_hashtab_matrix_field % Next
!          DO WHILE( ASSOCIATED(matrix_field_cell) )
!             ao_matrix_V => matrix_field_cell % ao_matrix
!             ptr_field_V  => matrix_field_cell % ptr_field
!
!             ivar_V  = matrix_field_cell % oi_loc_id 
!             nivar_V = ao_matrix_V % oi_nvar
!
!     
!
!             i_glob = Iddl(ivar_V, ipos, nivar)
!             CALL MATRIX_GET_GLOBAL_ROW_INDEX(ao_matrix_V, i_glob, i_real)
!
!	     i_glob = i_real
!	     IF (i_glob > 0) THEN
!                DO ijg = 1, ao_mesh % ptr_quad % oi_n_points
!
!                   lo_BBox2Di % ijg = ijg
! 
!                   lr_value = ao_matrix_V % opr_global_unknown(i_glob)
!                   CALL UPDATE_VARIABLES_GREENBOX_2D(ao_gbox, lo_BBox2Di, ptr_field_V % oi_id, lr_value)
!
!                END DO 
!	     END IF
!          
!             matrix_field_cell => matrix_field_cell % Next
!          END DO
!          ! -------------------------------------------------------------------
!       END DO
!       ! End of Loops 
!       ! -------------------------------------------------------------------
!      
!       CALL Reset_Local_Data()
!
!       CALL UPDATE_BASIS(ao_Basis2Di, ao_mesh, li_elmt_id)
!       CALL UPDATE_BASIS(ao_Basis2Dj, ao_mesh, li_elmt_id)
!        
!       ! -------------------------------------------------------------------
!       !  Put here nonlinear functions of variables computed at gauss points
!       !  Equation of state should be here
!       ! -------------------------------------------------------------------
!       CALL MappingGeom_At_GaussPt(lo_BBox2Di)
!       ! -------------------------------------------------------------------
!
!       ! -------------------------------------------------------------------
!       !  Compute unknowns gradient in the physical domain 
!       ! -------------------------------------------------------------------
!       CALL Physical_Variables_Pol(lo_BBox2Di, ao_gbox)
!       ! -------------------------------------------------------------------
!
!       ! ------------------------------------------ Loops on i ddl
!       ! Loop on poloidal vertices (i)
!       ! -------------------------------------------------------------------
!
!       DO i_loc1 = 1, lp_elmt % oi_nen
!          Scale2Di = 1.0 ! lp_elmt % Scales(or_pol1, iv_pol1)
!
!          ! -------------------------------------------------
!          CALL getPoloidalBasisFunction(i_loc1, li_one, &
!                                       & Scale2Di, &
!                                       & lo_BBox2Di)
!
!          CALL Physical_BasisFunction_Pol( &
!                  & lo_BBox2Di)
!          ! -------------------------------------------------
!          IF ( to_assemble_rhs ) THEN
!             ! -------------------------------------------------------------------
!             ! Loops on Gauss points
!             ! -------------------------------------------------------------------
!             ! TODO remove the loop from here. must be done inside RHS_for_Vi
!             DO ijg = 1,ao_mesh % ptr_quad % oi_n_points
!                lo_BBox2Di % ijg = ijg
!    
!                ! -------------------------------------------------
!                ! ---------------------------------------------
!                ! ---------------------------------------------
!                ! here RHS bloc-Vector  to be completed 
!                ! this RHS is associated  to 
!                !     i == (ipol1, itor1, or_pol1, or_tor1)
!                ! ---------------------------------------------
!                ! ---------------------------------------------
!                CALL ao_matrix % ptr_rhs_contribution(lo_BBox2Di, ao_gbox)
!                IF (ISNAN(ao_matrix % Rhs_Contribution(1))) THEN
!                   ao_matrix % Rhs_Contribution = 1.0E14 
!                END IF
!                
!                CALL SPM_ASSEMBLYRHSADDVALUES_LOCAL(ao_matrix % oi_matrix_id, i_loc1, &
!                        & ao_matrix % Rhs_Contribution, ierr)
!               
!                ! ... reset RHS contribution array
!                CALL MATRIX_RESET_ELEMENT_RHS(ao_matrix)
!                ! ...
!             END DO 
!             ! -------------------------------------------------------------------
!          END IF
!
!          IF ( to_assemble_matrix ) THEN
!             DO i_loc2 = 1, lp_elmt % oi_nen
!                Scale2Dj = 1.0 ! lp_elmt % Scales(or_pol2, iv_pol2)
!                
!                ! -------------------------------------------------
!                CALL getPoloidalBasisFunction(i_loc2, li_one, &
!                                             & Scale2Dj, &
!                                             & lo_BBox2Dj)
!                
!                CALL Physical_BasisFunction_Pol( &
!                        & lo_BBox2Dj)
!                ! -------------------------------------------------
!                
!                ! -------------------------------------------------------------------
!                ! Loops on Gauss points
!                ! -------------------------------------------------------------------
!                ! TODO remove the loop from here. must be done inside Matrix_for_Vi_Vj
!                DO ijg = 1,ao_mesh % ptr_quad % oi_n_points
!                   lo_BBox2Di % ijg = ijg
!               
!                   CALL ao_matrix % ptr_matrix_contribution(lo_BBox2Di, lo_BBox2Dj, ao_gbox)
!                
!                   IF (ISNAN(ao_matrix % Matrix_Contribution(1,1))) THEN
!                      ao_matrix % Matrix_Contribution = 1.0E14
!                    END IF
!                    
!                    CALL SPM_ASSEMBLYADDVALUES_LOCAL(ao_matrix % oi_matrix_id, i_loc1, i_loc2, &
!                            & ao_matrix % Matrix_Contribution, ierr)
!                    
!                    ! ... reset Matrix contribution array
!                    CALL MATRIX_RESET_ELEMENT_MATRIX(ao_matrix)
!                    ! ...
!                 END DO 
!                 ! -------------------------------------------------------------------
!                END DO
!             END IF
!
!       END DO
!       ! -------------------------------------------------------------------
!
!     
!       ! -- UPDATE RHS AND MATRICES FROM LOCAL-RHS
!       nivar = ao_matrix % oi_nvar 
!       njvar = nivar 
!
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
!    ! -------------------------------------------------------------------
!    ! -------------------------------------------------------------------
!
!    ! ----------------------------------------------------
!    ! End of the Assembly process
!    ! ----------------------------------------------------
    CALL ASSEMBLYEND_MATRIX(ao_matrix)
!    ! ----------------------------------------------------
!
!    ! ... get the RHS from SPM
!    !     here we get the global RHS, 
!    !     for the moment ROOT argument is not used
!    !     the local RHS (for the owener proc) can be accessed by SPM_GetLocalRHS
!    IF (to_assemble_rhs) THEN
!       li_root = -1
!       CALL SPM_GetGlobalRHS(ao_matrix % oi_matrix_id, ao_matrix % opr_global_rhs, li_root, ierr)
!    END IF
!    ! ...
!
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
!      IF ( to_assemble_matrix ) THEN
!         CALL SPM_MATRIXRESET_LOCAL(ao_matrix % oi_matrix_id, ierr)
!      END IF
!      IF ( to_assemble_rhs ) THEN
!         CALL SPM_RHSRESET_LOCAL(ao_matrix % oi_matrix_id, ierr)
!      END IF
    END SUBROUTINE Reset_Local_Data
    ! -------------------------------------------------------------------

  END SUBROUTINE LOOP_ON_ELMTS_1D_BSPLINE
    ! -------------------------------------------------------------------
    ! -------------------------------------------------------------------
END MODULE SPI_ASSEMBLY_1D_BSPLINE
