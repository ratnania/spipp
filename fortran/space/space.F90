!# -*- coding: utf8 -*-
MODULE SPACE 
  USE TypeDef        
  USE SPM_DEF
  USE SPM
  USE DIRICHLET_MOD
  USE NUMBERING
  USE JOREK_GLOB
  USE GRAPH
  USE INDICES_DEF
  USE FEBasis
  USE tracelog_module
  USE SPM_DEF
  USE JOREK_PARAM
  USE JOREK_PARAM_DEF
  USE SPACE_DEF
  IMPLICIT NONE

!  PRIVATE

  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 2 

CONTAINS
  ! ..................................................

  ! ..................................................
  SUBROUTINE SPACE_CREATE_QUAD_2D(self, Mesh2D)
  IMPLICIT NONE
     CLASS(DEF_SPACE_QUAD_2D), TARGET, INTENT(INOUT) :: self 
     CLASS(DEF_MESH_2D), TARGET, INTENT(INOUT) :: Mesh2D 
     ! LOCAL
     INTEGER :: li_poloidal_basis
     INTEGER :: li_err
     INTEGER :: li_quadrature_rule
     INTEGER :: li_n1
     INTEGER :: li_n2
     INTEGER, PARAMETER :: N_DIM = 2

#ifdef DEBUG_TRACE 
     CALL printlog("SPACE_CREATE_QUAD_2D: Begin", ai_dtllevel = 0)
#endif
     
     CALL JOREK_Param_GETInt(INT_QUADRATURE_RULE_ID, li_quadrature_rule, li_err)
     CALL JOREK_Param_GETInt(INT_QUADRATURE_N1_ID, li_n1, li_err)
     CALL JOREK_Param_GETInt(INT_QUADRATURE_N2_ID, li_n2, li_err)
     
     CALL CREATE_QUADRATURE(self % quadrature, N_DIM, &
             & li_quadrature_rule, li_n1,li_n2, &
             & dirname=self % dirname)
     
     Mesh2D % ptr_quad => self % quadrature
     self % ptr_mesh => Mesh2D
    
     CALL JOREK_Param_GETInt(INT_POLOIDAL_BASIS_ID, li_poloidal_basis, li_err)
    
     CALL CREATE_BASIS(self % basis, Mesh2D, li_poloidal_basis, dirname=self % dirname)

#ifdef DEBUG_TRACE 
     CALL printlog("SPACE_CREATE_QUAD_2D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE SPACE_CREATE_QUAD_2D
  ! ..................................................

  ! ..................................................
  SUBROUTINE SPACE_CREATE_TRIANGLE_2D(self, Mesh2D)
  IMPLICIT NONE
     CLASS(DEF_SPACE_TRIANGLE_2D), TARGET, INTENT(INOUT) :: self 
     CLASS(DEF_MESH_2D), TARGET, INTENT(INOUT) :: Mesh2D 
     ! LOCAL
     INTEGER :: li_poloidal_basis
     INTEGER :: li_quadrature_rule
     INTEGER :: li_err
     INTEGER :: li_n1
     INTEGER, PARAMETER :: N_DIM = 2

#ifdef DEBUG_TRACE 
     CALL printlog("SPACE_CREATE_TRIANGLE_2D: Begin", ai_dtllevel = 0)
#endif

     CALL JOREK_Param_GETInt(INT_QUADRATURE_RULE_ID, li_quadrature_rule, li_err)
     CALL JOREK_Param_GETInt(INT_QUADRATURE_N1_ID, li_n1, li_err)

     CALL CREATE_QUADRATURE(self % quadrature, N_DIM, &
             & li_quadrature_rule, li_n1, &
             & dirname=self % dirname)

     Mesh2D % ptr_quad => self % quadrature
     self % ptr_mesh => Mesh2D

     CALL JOREK_Param_GETInt(INT_POLOIDAL_BASIS_ID, li_poloidal_basis, li_err)

     CALL CREATE_BASIS(self % basis, Mesh2D, li_poloidal_basis, &
             & dirname=self % dirname)

#ifdef DEBUG_TRACE 
     CALL printlog("SPACE_CREATE_TRIANGLE_2D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE SPACE_CREATE_TRIANGLE_2D
  ! ..................................................

  ! ..................................................
  SUBROUTINE SPACE_CREATE_QUAD_2D1D(self, Mesh2D, Mesh1D)
  IMPLICIT NONE
     CLASS(DEF_SPACE_QUAD_2D1D), TARGET, INTENT(INOUT) :: self 
     CLASS(DEF_MESH_2D), TARGET, INTENT(INOUT) :: Mesh2D 
     CLASS(DEF_MESH_1D), TARGET, INTENT(INOUT) :: Mesh1D 
     ! LOCAL
     INTEGER :: li_poloidal_basis
     INTEGER :: li_toroidal_basis
     INTEGER :: li_err
     INTEGER :: li_n_period
     INTEGER :: li_n_plane
     INTEGER :: li_quadrature_rule_2d
     INTEGER :: li_quadrature_rule_1d
     INTEGER :: li_n1
     INTEGER :: li_n2
     INTEGER :: li_n3
     INTEGER, PARAMETER :: DIM_2D = 2
     INTEGER, PARAMETER :: DIM_1D = 1

#ifdef DEBUG_TRACE 
     CALL printlog("SPACE_CREATE_QUAD_2D1D: Begin", ai_dtllevel = 0)
#endif

     ! ... TODO must have 2 different quad rule param in parameter.txt
     CALL JOREK_Param_GETInt(INT_QUADRATURE_RULE_ID, li_quadrature_rule_2d, li_err)
     li_quadrature_rule_1d = JOREK_QUADRATURES_DEFAULT 
     ! ...
     CALL JOREK_Param_GETInt(INT_QUADRATURE_N1_ID, li_n1, li_err)
     CALL JOREK_Param_GETInt(INT_QUADRATURE_N2_ID, li_n2, li_err)
     CALL JOREK_Param_GETInt(INT_QUADRATURE_N3_ID, li_n3, li_err)

     CALL JOREK_Param_GETInt(INT_N_PLANE_ID,li_n_plane,li_err)
     CALL JOREK_Param_GETInt(INT_N_PERIOD_ID,li_n_period,li_err)

     CALL JOREK_Param_GETInt(INT_POLOIDAL_BASIS_ID, li_poloidal_basis, li_err)
     CALL JOREK_Param_GETInt(INT_TOROIDAL_BASIS_ID, li_toroidal_basis, li_err)
     
     CALL CREATE_QUADRATURE(self % quadrature_2d, DIM_2D, &
             & li_quadrature_rule_2d, li_n1, li_n2, &
             & dirname=self % dirname)

     !    CALL JOREK_Param_GETInt(INT_N_TOR_ID,n_tor,ierror)

     IF(li_toroidal_basis .EQ. INT_TOROIDAL_BASIS_HBEZIER) THEN
        ! ... TODO uncomment when needed
!        CALL CREATE_QUADRATURE(self % quadrature_1d, DIM_1D, &
!                & li_quadrature_rule_1d, li_n3, &
!                & dirname=self % dirname_other)

        CALL CREATE_QUADRATURE(self % quadrature_1d, DIM_1D, &
                & li_quadrature_rule_1d, li_n3)
     ENDIF

     IF(li_toroidal_basis .EQ. INT_TOROIDAL_BASIS_FOURIER) THEN
        ! ... TODO uncomment when needed
!        CALL CREATE_QUADRATURE(self % quadrature_1d, DIM_1D, &
!                & JOREK_QUADRATURES_QUAD_FOURIER, li_n_plane, ai_n_period=li_n_period, &
!                & dirname=self % dirname_other)

        CALL CREATE_QUADRATURE(self % quadrature_1d, DIM_1D, &
                & JOREK_QUADRATURES_QUAD_FOURIER, li_n_plane, ai_n_period=li_n_period)

     ENDIF  
    
     Mesh2D % ptr_quad => self % quadrature_2d
     Mesh1D % ptr_quad => self % quadrature_1d
    
     self % ptr_mesh_2d => Mesh2D
     self % ptr_mesh_1d => Mesh1D

    CALL CREATE_BASIS(self % basis_2d, Mesh2D, li_poloidal_basis, dirname=self % dirname)
    ! ... TODO uncomment when needed
!    CALL CREATE_BASIS(self % basis_1d, Mesh1D, li_toroidal_basis, dirname=self % dirname)

    CALL CREATE_BASIS(self % basis_1d, Mesh1D, li_toroidal_basis)

#ifdef DEBUG_TRACE 
     CALL printlog("SPACE_CREATE_QUAD_2D1D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE SPACE_CREATE_QUAD_2D1D
  ! ..................................................

  ! ..................................................
  SUBROUTINE SPACE_CREATE(self, ao_mesh, dirname, ao_mesh_other, dirname_other)
  IMPLICIT NONE
     CLASS(DEF_SPACE_ABSTRACT), TARGET          , INTENT(INOUT) :: self 
     CLASS(DEF_MESH_ABSTRACT) , TARGET          , INTENT(INOUT) :: ao_mesh
     CHARACTER(LEN = 1024)            , OPTIONAL, INTENT(IN)    :: dirname
     CLASS(DEF_MESH_ABSTRACT) , TARGET, OPTIONAL, INTENT(INOUT) :: ao_mesh_other
     CHARACTER(LEN = 1024)            , OPTIONAL, INTENT(IN)    :: dirname_other
     ! LOCAL

#ifdef DEBUG_TRACE 
     CALL printlog("SPACE_CREATE: Begin", ai_dtllevel = 0)
#endif

     IF (PRESENT(dirname)) THEN
        self % dirname = TRIM(dirname)
     END IF

     IF (PRESENT(dirname_other)) THEN
        self % dirname_other = TRIM(dirname_other)
     END IF

     SELECT TYPE (self)
     ! ... 2D Quad case 
     CLASS IS (DEF_SPACE_QUAD_2D)
        SELECT TYPE (ao_mesh)
        ! ... 
        CLASS IS (DEF_MESH_2D)
           CALL SPACE_CREATE_QUAD_2D(self, ao_mesh)
        ! ...
        CLASS DEFAULT
           STOP 'SPACE_CREATE: unexpected type for ao_mesh object!'
        END SELECT
     ! ...
     ! ... 2D Triangle case 
     CLASS IS (DEF_SPACE_TRIANGLE_2D)
        SELECT TYPE (ao_mesh)
        ! ... 
        CLASS IS (DEF_MESH_2D)
           CALL SPACE_CREATE_TRIANGLE_2D(self, ao_mesh)
        ! ...
        CLASS DEFAULT
           STOP 'SPACE_CREATE: unexpected type for ao_mesh object!'
        END SELECT
     ! ...
     ! ... 2D1D Quad case 
     CLASS IS (DEF_SPACE_QUAD_2D1D)
        SELECT TYPE (ao_mesh)
        ! ... 
        CLASS IS (DEF_MESH_2D)
           IF (PRESENT(ao_mesh_other)) THEN
              SELECT TYPE (ao_mesh_other)
              ! ... mesh 1D 
              CLASS IS (DEF_MESH_1D)
                 CALL SPACE_CREATE_QUAD_2D1D(self, ao_mesh, ao_mesh_other)
             ! ...
             CLASS DEFAULT
                STOP 'SPACE_CREATE: unexpected type for ao_mesh_other object!'
             END SELECT
           ELSE
              STOP 'SPACE_CREATE: unexpected type for ao_space_trial object!'
           END IF

        ! ...
        CLASS DEFAULT
           STOP 'SPACE_CREATE: unexpected type for ao_mesh object!'
        END SELECT
     ! ...
     ! ... TODO uncomment when needed
!     ! ... 2D1D Triangle case 
!     CLASS IS (DEF_SPACE_TRIANGLE_2D1D)
!        SELECT TYPE (ao_mesh)
!        ! ... 
!        CLASS IS (DEF_MESH_2D)
!           IF (PRESENT(ao_mesh_other)) THEN
!              SELECT TYPE (ao_mesh_other)
!              ! ... mesh 1D 
!              CLASS IS (DEF_MESH_1D)
!                 CALL SPACE_CREATE_TRIANGLE_2D1D(self, ao_mesh, ao_mesh_other)
!             ! ...
!             CLASS DEFAULT
!                STOP 'SPACE_CREATE: unexpected type for ao_mesh_other object!'
!             END SELECT
!           ELSE
!              STOP 'SPACE_CREATE: unexpected type for ao_space_trial object!'
!           END IF
!
!        ! ...
!        CLASS DEFAULT
!           STOP 'SPACE_CREATE: unexpected type for ao_mesh object!'
!        END SELECT
!     ! ...
     CLASS DEFAULT
        STOP 'SPACE_CREATE: unexpected type for self object!'
     END SELECT

#ifdef DEBUG_TRACE 
     CALL printlog("SPACE_CREATE: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE SPACE_CREATE
  ! ..................................................

END MODULE SPACE 
