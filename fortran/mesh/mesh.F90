!# -*- coding: utf8 -*-
MODULE SPI_MESH 
  USE SPI_MESH_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  IMPLICIT NONE

CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_MESH(self, ai_n, ai_p, ai_type_bc, dirname)
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_MESH_ABSTRACT)       , INTENT(INOUT) :: self 
     INTEGER               , OPTIONAL, INTENT(IN)   :: ai_n
     INTEGER               , OPTIONAL, INTENT(IN)   :: ai_p
     INTEGER               , OPTIONAL, INTENT(IN)   :: ai_type_bc
     CHARACTER(LEN = 1024) , OPTIONAL, INTENT(IN)   :: dirname
     ! LOCAL
     CHARACTER(LEN = 1024) :: ls_dirname
     INTEGER :: li_err     

     IF (PRESENT(dirname)) THEN
        ls_dirname = TRIM(dirname)
     ELSE
        CALL GETCWD(ls_dirname) 
     END IF

     self % oi_n_vtex_per_elmt = 2  

     SELECT TYPE (self)
     CLASS IS (DEF_MESH_1D_BSPLINE)
        CALL CREATE_MESH_BSPLINES_1D(self, &
                & ai_n=ai_n, &
                & ai_p=ai_p, & 
                & ai_type_bc=ai_type_bc, &
                & dirname=dirname)

     CLASS IS (DEF_MESH_1D_FOURIER)
        CALL CREATE_MESH_FOURIER_1D(self, ls_dirname)

     CLASS IS (DEF_MESH_1D_HBEZIER)
        CALL CREATE_MESH_HBEZIER_1D(self, ls_dirname)

     CLASS DEFAULT
        STOP 'CREATE_MESH: unexpected type for self object!'
     END SELECT

  END SUBROUTINE CREATE_MESH
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_BSPLINES_1D(self, ai_n, ai_p, ai_type_bc, dirname)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     INTEGER              , OPTIONAL, INTENT(IN)    :: ai_n
     INTEGER              , OPTIONAL, INTENT(IN)    :: ai_p
     INTEGER              , OPTIONAL, INTENT(IN)    :: ai_type_bc
     CHARACTER(LEN = 1024), OPTIONAL, INTENT(IN)    :: dirname
     ! LOCAL
     INTEGER :: li_err 

     IF (PRESENT(dirname)) THEN
        ! ... TODO
        STOP "not yet implemented"

     ELSEIF ( (PRESENT(ai_n)) .AND. &
            & (PRESENT(ai_p)) .AND. &
            & (PRESENT(ai_type_bc))) THEN

        self % oi_type_bc = ai_type_bc
        self % oi_p       = ai_p
        self % oi_n       = ai_n
       
        if ( ai_type_bc == SPI_BC_PERIODIC ) then
           self % oi_n =  self % oi_n + self % oi_p - 1			
        end if
                                        
        self % oi_n_elmts = self % oi_n - ai_p
        self % oi_nnp     = self % oi_n
        self % oi_nen     = ai_p + 1
        
        ALLOCATE ( self % opr_knot ( self % oi_n + ai_p + 1 ) ) 
        ALLOCATE ( self % opr_control_points( self % oi_n ) ) 
        
        ! INITIALIZING VECTOR KNOTS
        call INIT_MESH_1D_BSPLINE( self, ai_type_bc )		
     ELSE
        STOP "Wrong arguments"
     END IF

  END SUBROUTINE CREATE_MESH_BSPLINES_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_FOURIER_1D(self, dirname)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_FOURIER), INTENT(INOUT)  :: self
     CHARACTER(LEN = 1024), INTENT(IN) :: dirname
     ! LOCAL
     INTEGER :: li_err 

  END SUBROUTINE CREATE_MESH_FOURIER_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_HBEZIER_1D(self, dirname)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_HBEZIER), INTENT(INOUT)  :: self
     CHARACTER(LEN = 1024), INTENT(IN) :: dirname
     ! LOCAL
     INTEGER :: li_err 

  END SUBROUTINE CREATE_MESH_HBEZIER_1D
  ! .........................................................
  
  ! .........................................................
  SUBROUTINE FREE_MESH(self)
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_MESH_ABSTRACT)       , INTENT(INOUT) :: self 
     ! LOCAL

     SELECT TYPE (self)
     CLASS IS (DEF_MESH_1D_BSPLINE)
        CALL FREE_MESH_BSPLINES_1D(self)

     CLASS IS (DEF_MESH_1D_FOURIER)
        CALL FREE_MESH_FOURIER_1D(self)

     CLASS IS (DEF_MESH_1D_HBEZIER)
        CALL FREE_MESH_HBEZIER_1D(self)

     CLASS DEFAULT
        STOP 'FREE_MESH: unexpected type for self object!'
     END SELECT

  END SUBROUTINE FREE_MESH
  ! .........................................................

  ! .........................................................
  SUBROUTINE FREE_MESH_BSPLINES_1D(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: li_err 

     DEALLOCATE ( self % opr_knot ) 

  END SUBROUTINE FREE_MESH_BSPLINES_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE FREE_MESH_FOURIER_1D(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_FOURIER), INTENT(INOUT)  :: self
     ! LOCAL
     INTEGER :: li_err 

  END SUBROUTINE FREE_MESH_FOURIER_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE FREE_MESH_HBEZIER_1D(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_HBEZIER), INTENT(INOUT)  :: self
     ! LOCAL
     INTEGER :: li_err 

  END SUBROUTINE FREE_MESH_HBEZIER_1D
  ! .........................................................

  ! .........................................................
  subroutine INIT_MESH_1D_BSPLINE ( self, ai_type )
  implicit none
     TYPE(DEF_MESH_1D_BSPLINE), INTENT(INOUT) :: self
     integer  :: ai_type
     ! LOCAL VARIABLES
     integer  :: li_i
     integer  :: li_nu ! number of continuity condition for periodic vector knot		
     real(SPI_RK), dimension ( self % oi_n + self % oi_p + 1 ) :: lpr_knot

     if ( self % oi_n < self % oi_p + 1 ) then
        STOP "Error INIT_MESH_1D_BSPLINE: you must have N >= p + 1"
     end if
     
     ! ... knots
     self % opr_knot ( 1 : self % oi_p + 1 ) = 0.0
     do li_i = 1, self % oi_n - self % oi_p - 1
        self % opr_knot ( self % oi_p + 1 + li_i ) = li_i * 1.0 / ( self % oi_n - self % oi_p )
     end do 
     self % opr_knot ( self % oi_n + 1 : self % oi_n + self % oi_p + 1 ) = 1.0
                                     
     if ( ai_type == SPI_BC_PERIODIC ) then
        li_nu = self % oi_p
        
        lpr_knot ( : ) = self % opr_knot ( : )
        call convert_to_periodic_knots ( lpr_Knot ( : ), self % oi_n, self % oi_p, li_nu, self % opr_knot ( : ) )
     end if
     ! ...

     ! ... control points
     self % opr_control_points = 0.0
     do li_i = 1, self % oi_n 
        self % opr_control_points( li_i ) = (li_i - 1) * 1.0 / (self % oi_n - 1)
     end do 
     ! ...
  end subroutine INIT_MESH_1D_BSPLINE
  ! .........................................................

  ! .........................................................
  subroutine convert_to_periodic_knots ( apr_knot, ai_N, ai_P, ai_nu, apr_t )
  !> CONVERT A KNOT VECTOR INTO A PERIODIC KNOT (TO GET PERIODIC B-SPLINES) 
  !> WITH THE POSSIBILITY OF CONTROLLING THE ORDER OF THE CONTINUITY
  !> \todo NOT TESTED
  implicit none
    !> \param[in]
    real(SPI_RK), dimension(:) :: apr_knot
    !> \param[in]
    integer  :: ai_N
    !> \param[in]
    integer  :: ai_P
    !> \param[in]	Nbr OF CONITUITY CONDITIONS
    integer  :: ai_nu	
    !> \param[out] 
    real(SPI_RK), dimension(:) :: apr_t	
    ! LOCAL VARIABLES 
    integer  :: li_i	
    real(SPI_RK) :: lr_L !THE PERIOD
    
    lr_L = apr_knot ( ai_N + ai_P + 1 ) - apr_knot ( 1 )
    
    do li_i = 1 , ai_nu
       apr_t ( li_i ) = apr_knot ( ai_N - ai_nu + li_i ) - lr_L
    end do
    
    do li_i = ai_N + ai_P + 1 - ai_nu + 1 , ai_n + ai_P + 1
       apr_t ( li_i ) = apr_knot ( li_i - ( ai_n + ai_P + 1 - ai_nu + 1 ) + ai_P + 2 ) + lr_L
    end do	
  
  end subroutine convert_to_periodic_knots
  ! .........................................................


  ! .........................................................
!  SUBROUTINE ToroidalAllocation(self, ai_n_nodes_Tor, ai_n_plane_Tor, ai_toroidal_basis)
!  IMPLICIT NONE
!    TYPE(DEF_MESH_1D)          :: self 
!    INTEGER, INTENT(IN)        :: ai_n_nodes_Tor
!    INTEGER, INTENT(IN)        :: ai_n_plane_Tor
!    INTEGER, INTENT(IN)        :: ai_toroidal_basis
!    ! LOCAL
!
!#ifdef DEBUG_TRACE 
!    CALL printlog("ToroidalAllocation: Begin", ai_dtllevel = 0)
!#endif
!
!    self % oi_n_nodes = ai_n_nodes_Tor
!    self % oi_n_plane = ai_n_plane_Tor 
!
!    SELECT CASE(ai_toroidal_basis)
!       CASE(INT_TOROIDAL_BASIS_HBEZIER) 
!          ALLOCATE(self % AngleT  ( ai_n_nodes_Tor+1) )
!
!       CASE(INT_TOROIDAL_BASIS_FOURIER) 
!          ALLOCATE(self % AngleT  ( ai_n_plane_Tor) )
!    END SELECT
!
!#ifdef DEBUG_TRACE 
!    CALL printlog("ToroidalAllocation: End", ai_dtllevel = 0)
!#endif
!  END SUBROUTINE ToroidalAllocation
  ! .........................................................

  ! .........................................................
!  SUBROUTINE ToroidalInitialization( self, &
!                  & ai_assembly_proc, ai_toroidal_basis, &
!                  & ai_n_order, ai_n_order_Tor)
!    TYPE(DEF_MESH_1D), INTENT(INOUT) :: self 
!    INTEGER, INTENT(IN)        :: ai_assembly_proc
!    INTEGER, INTENT(IN)        :: ai_toroidal_basis
!    INTEGER, INTENT(IN)        :: ai_n_order
!    INTEGER, INTENT(IN)        :: ai_n_order_Tor
!    ! LOCAL
!    INTEGER              :: i_order, i_tor, iv
!    INTEGER :: li_err 
!
!#ifdef DEBUG_TRACE 
!    CALL printlog("ToroidalInitialization: Begin", ai_dtllevel = 0)
!#endif
!
!    self % oi_n_order = ai_n_order_Tor
!!    self % oi_ = 
!!    self % oi_ = 
!!    self % oi_ = 
!
!    SELECT CASE(ai_assembly_proc)
!    CASE(INT_ASSEMBLY_PROCEDURE_1D) 
!       CONTINUE
!    CASE(INT_ASSEMBLY_PROCEDURE_2D) 
!       CONTINUE
!    CASE(INT_ASSEMBLY_PROCEDURE_2D1D)
!       CALL ToroidalMesh(self, ai_toroidal_basis)
!    CASE(INT_ASSEMBLY_PROCEDURE_3D)
!       CONTINUE
!     END SELECT
!
!#ifdef DEBUG_TRACE 
!    CALL printlog("ToroidalInitialization: End", ai_dtllevel = 0)
!#endif
!
!  END SUBROUTINE ToroidalInitialization
  ! .........................................................

  ! .........................................................
!  SUBROUTINE ToroidalMesh(self, ai_toroidal_basis)
!    TYPE(DEF_MESH_1D)          :: self 
!    INTEGER          :: iplan
!    REAL(KIND=RK)    :: Dphi
!    INTEGER, INTENT(IN)        :: ai_toroidal_basis
!
!    ! ----------------------------------
!    IF (self % oi_n_nodes< 1) THEN
!	    PRINT *, "ToroidalMesh: n_nodes_Tor must be >= 2. Jorek will stop immediatly"
!            PRINT *, "ToroidalMesh: current value of n_nodes ", self % oi_n_nodes
!	    STOP
!    END IF
!
!    SELECT CASE(ai_toroidal_basis)
!       CASE(INT_TOROIDAL_BASIS_HBEZIER) 
!          Dphi = 2.0*Pi/REAL(self % oi_n_nodes)
!
!          DO iplan = 1, self % oi_n_nodes+1
!             self % AngleT(iplan) = REAL( iplan -1)*Dphi
!          END DO
!
!       CASE(INT_TOROIDAL_BASIS_FOURIER) 
!          Dphi = 2.0*Pi/REAL(self % oi_n_plane)
!
!          DO iplan = 1, self % oi_n_plane
!             self % AngleT(iplan) = REAL( iplan -1)*Dphi
!          END DO
!    END SELECT
!
!  END SUBROUTINE ToroidalMesh
  ! .........................................................

  ! .........................................................
!  SUBROUTINE MESH_TRANSLATE(self, apr_displ)
!  IMPLICIT NONE
!     CLASS(DEF_MESH_ABSTRACT)          , INTENT(INOUT) :: self 
!     REAL(KIND=RK), DIMENSION(:)       :: apr_displ
!     ! LOCAL
!     INTEGER :: li_mesh
!     INTEGER :: li_err     
!
!#ifdef DEBUG_TRACE 
!     CALL printlog("MESH_TRANSLATE: Begin", ai_dtllevel = 0)
!#endif
!
!     SELECT TYPE (self)
!     ! ... 2D case 
!     CLASS IS (DEF_MESH_2D)
!        CALL MESH_TRANSLATE_2D(self, apr_displ)
!     ! ...
!     CLASS DEFAULT
!        STOP 'MESH_TRANSLATE: unexpected type for self object!'
!     END SELECT
!
!#ifdef DEBUG_TRACE 
!     CALL printlog("MESH_TRANSLATE: End", ai_dtllevel = 0)
!#endif
!
!  END SUBROUTINE MESH_TRANSLATE
  ! .........................................................

  ! .........................................................
END MODULE SPI_MESH 
