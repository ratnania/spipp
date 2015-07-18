!# -*- coding: utf8 -*-
MODULE SPI_MESH 
  USE SPI_MESH_DEF
  USE SPI_GLOBAL_DEF
!  USE SPI_MESH_HBEZIER_FILE
!  USE SPI_MESH_BSPLINE_FILE
  USE SPI_MESH_TRANSFORM
  IMPLICIT NONE

CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_MESH(self, other, dirname)
  ! ... both other and dirname are optional
  !     other is used for a tensor geometry (2d1d for example)
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_MESH_ABSTRACT)          , INTENT(INOUT) :: self 
     CLASS(DEF_MESH_ABSTRACT), OPTIONAL, INTENT(INOUT) :: other 
     CHARACTER(LEN = 1024)   , OPTIONAL, INTENT(IN)    :: dirname
     ! LOCAL
     CHARACTER(LEN = 1024) :: ls_dirname
     INTEGER :: li_mesh
     INTEGER :: li_err     

#ifdef DEBUG_TRACE 
     CALL printlog("CREATE_MESH: Begin", ai_dtllevel = 0)
#endif

     CALL JOREK_Param_GETInt(INT_TYPEMESH_ID,li_mesh,li_err)

     IF (PRESENT(dirname)) THEN
        ls_dirname = TRIM(dirname)
     ELSE
        CALL GETCWD(ls_dirname) 
     END IF

     SELECT TYPE (self)
     ! ... 2D case 
     CLASS IS (DEF_MESH_1D)
        SELECT CASE(li_mesh)
        CASE(INT_MESH_CAID_HBEZIER) 
           CALL CREATE_MESH_BSPLINES_1D(self, ls_dirname)

        CASE(INT_MESH_BOXSPLINES) 
           CALL CREATE_MESH_FOURIER_1D(self, ls_dirname)

        CASE DEFAULT
           PRINT *, "ERROR in CREATE_MESH: Not a recognized mesh type"
           STOP
        END SELECT
     ! ...
     ! ... TODO 1Dx1D case and other cases
     CLASS DEFAULT
        STOP 'CREATE_MESH: unexpected type for self object!'
     END SELECT

#ifdef DEBUG_TRACE 
     CALL printlog("CREATE_MESH: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE CREATE_MESH
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_BSPLINES_1D(Mesh1D, dirname)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D), INTENT(INOUT)  :: Mesh1D
     CHARACTER(LEN = 1024), INTENT(IN) :: dirname
     ! LOCAL
     INTEGER :: li_err 

#ifdef DEBUG_TRACE 
     CALL printlog("CREATE_MESH_BSPLINES_1D: Begin", ai_dtllevel = 0)
#endif


#ifdef DEBUG_TRACE 
     CALL printlog("CREATE_MESH_BSPLINES_1D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE CREATE_MESH_SPLINES_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_FOURIER_1D(Mesh1D, dirname)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D), INTENT(INOUT)  :: Mesh1D
     CHARACTER(LEN = 1024), INTENT(IN) :: dirname
     ! LOCAL
     INTEGER :: li_err 

#ifdef DEBUG_TRACE 
     CALL printlog("CREATE_MESH_FOURIER_1D: Begin", ai_dtllevel = 0)
#endif


#ifdef DEBUG_TRACE 
     CALL printlog("CREATE_MESH_FOURIER_1D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE CREATE_MESH_FOURIER_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE ToroidalAllocation(Mesh1D, ai_n_nodes_Tor, ai_n_plane_Tor, ai_toroidal_basis)
  IMPLICIT NONE
    TYPE(DEF_MESH_1D)          :: Mesh1D 
    INTEGER, INTENT(IN)        :: ai_n_nodes_Tor
    INTEGER, INTENT(IN)        :: ai_n_plane_Tor
    INTEGER, INTENT(IN)        :: ai_toroidal_basis
    ! LOCAL

#ifdef DEBUG_TRACE 
    CALL printlog("ToroidalAllocation: Begin", ai_dtllevel = 0)
#endif

    Mesh1D % oi_n_nodes = ai_n_nodes_Tor
    Mesh1D % oi_n_plane = ai_n_plane_Tor 

    SELECT CASE(ai_toroidal_basis)
       CASE(INT_TOROIDAL_BASIS_HBEZIER) 
          ALLOCATE(Mesh1D%AngleT  ( ai_n_nodes_Tor+1) )

       CASE(INT_TOROIDAL_BASIS_FOURIER) 
          ALLOCATE(Mesh1D%AngleT  ( ai_n_plane_Tor) )
    END SELECT

#ifdef DEBUG_TRACE 
    CALL printlog("ToroidalAllocation: End", ai_dtllevel = 0)
#endif
  END SUBROUTINE ToroidalAllocation
  ! .........................................................

  ! .........................................................
  SUBROUTINE ToroidalInitialization( Mesh1D, &
                  & ai_assembly_proc, ai_toroidal_basis, &
                  & ai_n_order, ai_n_order_Tor)
    TYPE(DEF_MESH_1D), INTENT(INOUT) :: Mesh1D 
    INTEGER, INTENT(IN)        :: ai_assembly_proc
    INTEGER, INTENT(IN)        :: ai_toroidal_basis
    INTEGER, INTENT(IN)        :: ai_n_order
    INTEGER, INTENT(IN)        :: ai_n_order_Tor
    ! LOCAL
    INTEGER              :: i_order, i_tor, iv
    INTEGER :: li_err 

#ifdef DEBUG_TRACE 
    CALL printlog("ToroidalInitialization: Begin", ai_dtllevel = 0)
#endif

    Mesh1D % oi_n_order = ai_n_order_Tor
!    Mesh1D % oi_ = 
!    Mesh1D % oi_ = 
!    Mesh1D % oi_ = 

    SELECT CASE(ai_assembly_proc)
    CASE(INT_ASSEMBLY_PROCEDURE_1D) 
       CONTINUE
    CASE(INT_ASSEMBLY_PROCEDURE_2D) 
       CONTINUE
    CASE(INT_ASSEMBLY_PROCEDURE_2D1D)
       CALL ToroidalMesh(Mesh1D, ai_toroidal_basis)
    CASE(INT_ASSEMBLY_PROCEDURE_3D)
       CONTINUE
     END SELECT

#ifdef DEBUG_TRACE 
    CALL printlog("ToroidalInitialization: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE ToroidalInitialization
  ! .........................................................

  ! .........................................................
  SUBROUTINE ToroidalMesh(Mesh1D, ai_toroidal_basis)
    TYPE(DEF_MESH_1D)          :: Mesh1D 
    INTEGER          :: iplan
    REAL(KIND=RK)    :: Dphi
    INTEGER, INTENT(IN)        :: ai_toroidal_basis

    ! ----------------------------------
    IF (Mesh1D % oi_n_nodes< 1) THEN
	    PRINT *, "ToroidalMesh: n_nodes_Tor must be >= 2. Jorek will stop immediatly"
            PRINT *, "ToroidalMesh: current value of n_nodes ", Mesh1D % oi_n_nodes
	    STOP
    END IF

    SELECT CASE(ai_toroidal_basis)
       CASE(INT_TOROIDAL_BASIS_HBEZIER) 
          Dphi = 2.0*Pi/REAL(Mesh1D % oi_n_nodes)

          DO iplan = 1, Mesh1D % oi_n_nodes+1
             Mesh1D%AngleT(iplan) = REAL( iplan -1)*Dphi
          END DO

       CASE(INT_TOROIDAL_BASIS_FOURIER) 
          Dphi = 2.0*Pi/REAL(Mesh1D % oi_n_plane)

          DO iplan = 1, Mesh1D % oi_n_plane
             Mesh1D%AngleT(iplan) = REAL( iplan -1)*Dphi
          END DO
    END SELECT

  END SUBROUTINE ToroidalMesh
  ! .........................................................

  ! .........................................................
  SUBROUTINE MESH_TRANSLATE(self, apr_displ)
  IMPLICIT NONE
     CLASS(DEF_MESH_ABSTRACT)          , INTENT(INOUT) :: self 
     REAL(KIND=RK), DIMENSION(:)       :: apr_displ
     ! LOCAL
     INTEGER :: li_mesh
     INTEGER :: li_err     

#ifdef DEBUG_TRACE 
     CALL printlog("MESH_TRANSLATE: Begin", ai_dtllevel = 0)
#endif

     SELECT TYPE (self)
     ! ... 2D case 
     CLASS IS (DEF_MESH_2D)
        CALL MESH_TRANSLATE_2D(self, apr_displ)
     ! ...
     ! ... TODO 1Dx1D case and other cases
     CLASS DEFAULT
        STOP 'MESH_TRANSLATE: unexpected type for self object!'
     END SELECT

#ifdef DEBUG_TRACE 
     CALL printlog("MESH_TRANSLATE: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE MESH_TRANSLATE
  ! .........................................................

  ! .........................................................
END MODULE SPI_MESH 
