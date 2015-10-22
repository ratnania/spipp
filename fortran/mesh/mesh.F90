!# -*- coding: utf8 -*-
MODULE SPI_MESH 
  USE SPI_MESH_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  IMPLICIT NONE

CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_MESH(self, ao_quad, ai_n, ai_p, ai_type_bc, apr_knots)
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_MESH_ABSTRACT)       , INTENT(INOUT) :: self 
     CLASS(DEF_QUADRATURE_ABSTRACT), TARGET , INTENT(INOUT) :: ao_quad 
     INTEGER               , INTENT(IN)   :: ai_n
     INTEGER               , INTENT(IN)   :: ai_p
     INTEGER               , OPTIONAL, INTENT(IN)   :: ai_type_bc
     real(SPI_RK), dimension (:), OPTIONAL, INTENT(IN) :: apr_knots
     ! LOCAL
     INTEGER :: li_err     

     self % oi_n_vtex_per_elmt = 2  

     SELECT TYPE (self)
     CLASS IS (DEF_MESH_1D_BSPLINE)
        CALL CREATE_MESH_BSPLINES_1D(self, &
                & ai_n, &
                & ai_p, & 
                & ai_type_bc=ai_type_bc, &
                & apr_knots=apr_knots &
                & )
        CALL CREATE_MESH_POINTS_1D(self, ao_quad) 

!     CLASS IS (DEF_MESH_1D_FOURIER)
!        CALL CREATE_MESH_FOURIER_1D(self)
!        CALL CREATE_MESH_POINTS_1D(self) 
!
!     CLASS IS (DEF_MESH_1D_HBEZIER)
!        CALL CREATE_MESH_HBEZIER_1D(self)
!        CALL CREATE_MESH_POINTS_1D(self) 

     CLASS DEFAULT
        STOP 'CREATE_MESH: unexpected type for self object!'
     END SELECT

  END SUBROUTINE CREATE_MESH
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_POINTS_1D(self, ao_quad)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     CLASS(DEF_QUADRATURE_ABSTRACT), TARGET , INTENT(INOUT) :: ao_quad 
     ! LOCAL
     INTEGER :: i_element
     INTEGER :: li_ig 
     REAL(SPI_RK) :: lr_a
     REAL(SPI_RK) :: lr_b
     REAL(SPI_RK) :: lr_x

     self % ptr_quad => ao_quad

     ALLOCATE(self % opr_points(self % n_elements, self % ptr_quad % oi_n_points))

     DO i_element=1, self % n_elements 
        lr_a = self % opr_grid(i_element)
        lr_b = self % opr_grid(i_element+1)

        DO li_ig = 1, self % ptr_quad % oi_n_points
           lr_x = lr_a + self % ptr_quad % opr_points(li_ig) * (lr_b - lr_a )
           self % opr_points(i_element, li_ig) = lr_x
        END DO
     END DO

  END SUBROUTINE CREATE_MESH_POINTS_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_GRID_BSPLINES_1D(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: n_elements
     INTEGER :: i
     INTEGER :: i_current
     REAL(SPI_RK), DIMENSION(self % oi_n + self % oi_p + 1) :: lpr_grid
     REAL(SPI_RK) :: min_current

     lpr_grid = SPI_INT_DEFAULT * 1.0 

     i_current = 1
     lpr_grid(i_current) = MINVAL(self % opr_knots)
     DO i=1, self % oi_n + self % oi_p
        min_current = MINVAL(self % opr_knots(i : self % oi_n + self % oi_p + 1))
        IF ( min_current > lpr_grid(i_current) ) THEN
                i_current = i_current + 1
                lpr_grid(i_current) = min_current
        END IF
     END DO
     n_elements = i_current - 1 
     self % n_elements = n_elements

     ALLOCATE ( self % opr_grid( n_elements + 1 ) ) 
     self % opr_grid (1:n_elements+1) = lpr_grid(1:n_elements+1)

  END SUBROUTINE CREATE_MESH_GRID_BSPLINES_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE INITIALIZE_MESH_KNOTS_BSPLINES_1D(self, ai_n, ai_p, apr_knots)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     INTEGER, INTENT(IN)    :: ai_n
     INTEGER, INTENT(IN)    :: ai_p
      real(SPI_RK), dimension (:), INTENT(IN):: apr_knots
     ! LOCAL

     self % oi_n = ai_n
     self % oi_p = ai_p

     ALLOCATE ( self % opr_knots( ai_n + ai_p + 1 ) ) 
     self % opr_knots = apr_knots

  END SUBROUTINE INITIALIZE_MESH_KNOTS_BSPLINES_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE INITIALIZE_MESH_KNOTS_DEFAULT_BSPLINES_1D(self, ai_n, ai_p, ai_type_bc)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     INTEGER, INTENT(IN)    :: ai_n
     INTEGER, INTENT(IN)    :: ai_p
     INTEGER, INTENT(IN)    :: ai_type_bc
     ! LOCAL
     INTEGER :: li_err 
     integer  :: li_i
     integer  :: li_nu ! number of continuity condition for periodic vector knot		
     real(SPI_RK), dimension ( ai_n + ai_p + 1 ) :: lpr_knot

     self % oi_type_bc = ai_type_bc
     self % oi_p       = ai_p
     self % oi_n       = ai_n
     
     if ( ai_type_bc == SPI_BC_PERIODIC ) then
        self % oi_n =  self % oi_n + self % oi_p - 1
     end if
                                     
     self % n_elements = self % oi_n - ai_p
     self % oi_nen     = ai_p + 1
     
     ALLOCATE ( self % opr_knots ( self % oi_n + ai_p + 1 ) ) 
     
     ! INITIALIZING VECTOR KNOTS
     if ( self % oi_n < self % oi_p + 1 ) then
        STOP "Error INIT_MESH_1D_BSPLINE: you must have N >= p + 1"
     end if
     
     ! ... knots
     self % opr_knots ( 1 : self % oi_p + 1 ) = 0.0
     do li_i = 1, self % oi_n - self % oi_p - 1
        self % opr_knots ( self % oi_p + 1 + li_i ) = li_i * 1.0 / ( self % oi_n - self % oi_p )
     end do 
     self % opr_knots ( self % oi_n + 1 : self % oi_n + self % oi_p + 1 ) = 1.0
                                     
     if ( ai_type_bc == SPI_BC_PERIODIC ) then
        li_nu = self % oi_p
        
        lpr_knot ( : ) = self % opr_knots ( : )
        call convert_to_periodic_knots ( lpr_Knot ( : ), self % oi_n, self % oi_p, li_nu, self % opr_knots ( : ) )
     end if
     ! ...

  END SUBROUTINE INITIALIZE_MESH_KNOTS_DEFAULT_BSPLINES_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_BSPLINES_1D(self, ai_n, ai_p, ai_type_bc, apr_knots)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     INTEGER              , INTENT(IN)    :: ai_n
     INTEGER              , INTENT(IN)    :: ai_p
     INTEGER              , OPTIONAL, INTENT(IN)    :: ai_type_bc
     real(SPI_RK), dimension (:), OPTIONAL, INTENT(IN) :: apr_knots
     ! LOCAL
     INTEGER :: li_err 

     IF ( (.NOT. PRESENT(apr_knots)) &
            & .AND. (PRESENT(ai_type_bc)) &
            & ) THEN
        CALL INITIALIZE_MESH_KNOTS_DEFAULT_BSPLINES_1D(self, ai_n, ai_p, ai_type_bc) 
     END IF

     IF ( ( PRESENT(apr_knots)) ) THEN
        CALL INITIALIZE_MESH_KNOTS_BSPLINES_1D(self, ai_n, ai_p, apr_knots) 
     END IF

     CALL CREATE_MESH_GRID_BSPLINES_1D(self) 

  END SUBROUTINE CREATE_MESH_BSPLINES_1D
  ! .........................................................

!  ! .........................................................
!  SUBROUTINE CREATE_MESH_FOURIER_1D(self)
!  IMPLICIT NONE
!     TYPE(DEF_MESH_1D_FOURIER), INTENT(INOUT)  :: self
!     ! LOCAL
!     INTEGER :: li_err 
!
!  END SUBROUTINE CREATE_MESH_FOURIER_1D
!  ! .........................................................
!
!  ! .........................................................
!  SUBROUTINE CREATE_MESH_HBEZIER_1D(self)
!  IMPLICIT NONE
!     TYPE(DEF_MESH_1D_HBEZIER), INTENT(INOUT)  :: self
!     ! LOCAL
!     INTEGER :: li_err 
!
!  END SUBROUTINE CREATE_MESH_HBEZIER_1D
!  ! .........................................................
  
  ! .........................................................
  SUBROUTINE FREE_MESH(self)
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_MESH_ABSTRACT)       , INTENT(INOUT) :: self 
     ! LOCAL

     DEALLOCATE ( self % opr_grid) 

     SELECT TYPE (self)
     CLASS IS (DEF_MESH_1D_BSPLINE)
        CALL FREE_MESH_BSPLINES_1D(self)

!     CLASS IS (DEF_MESH_1D_FOURIER)
!        CALL FREE_MESH_FOURIER_1D(self)
!
!     CLASS IS (DEF_MESH_1D_HBEZIER)
!        CALL FREE_MESH_HBEZIER_1D(self)

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

     DEALLOCATE ( self % opr_knots ) 

  END SUBROUTINE FREE_MESH_BSPLINES_1D
  ! .........................................................

!  ! .........................................................
!  SUBROUTINE FREE_MESH_FOURIER_1D(self)
!  IMPLICIT NONE
!     TYPE(DEF_MESH_1D_FOURIER), INTENT(INOUT)  :: self
!     ! LOCAL
!     INTEGER :: li_err 
!
!  END SUBROUTINE FREE_MESH_FOURIER_1D
!  ! .........................................................
!
!  ! .........................................................
!  SUBROUTINE FREE_MESH_HBEZIER_1D(self)
!  IMPLICIT NONE
!     TYPE(DEF_MESH_1D_HBEZIER), INTENT(INOUT)  :: self
!     ! LOCAL
!     INTEGER :: li_err 
!
!  END SUBROUTINE FREE_MESH_HBEZIER_1D
!  ! .........................................................

  ! .........................................................
  subroutine INIT_MESH_1D_BSPLINE ( self, ai_type )
  implicit none
     TYPE(DEF_MESH_1D_BSPLINE), INTENT(INOUT) :: self
     integer  :: ai_type
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
