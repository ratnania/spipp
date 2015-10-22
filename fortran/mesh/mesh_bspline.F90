!# -*- coding: utf8 -*-
MODULE SPI_MESH_BSPLINE 
  USE SPI_MESH_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  IMPLICIT NONE

CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_MESH_GRID_BSPLINES_1D(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: n_elements
     INTEGER :: i
     INTEGER :: i_current
     REAL(SPI_RK), DIMENSION(self % n + self % p + 1) :: lpr_grid
     REAL(SPI_RK) :: min_current

     lpr_grid = SPI_INT_DEFAULT * 1.0 

     i_current = 1
     lpr_grid(i_current) = MINVAL(self % knots)
     DO i=1, self % n + self % p
        min_current = MINVAL(self % knots(i : self % n + self % p + 1))
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

     self % n = ai_n
     self % p = ai_p

     ALLOCATE ( self % knots( ai_n + ai_p + 1 ) ) 
     self % knots = apr_knots

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
     self % p       = ai_p
     self % n       = ai_n
     
     if ( ai_type_bc == SPI_BC_PERIODIC ) then
        self % n =  self % n + self % p - 1
     end if
                                     
     self % n_elements = self % n - ai_p
     
     ALLOCATE ( self % knots ( self % n + ai_p + 1 ) ) 
     
     ! INITIALIZING VECTOR KNOTS
     if ( self % n < self % p + 1 ) then
        STOP "Error INIT_MESH_1D_BSPLINE: you must have N >= p + 1"
     end if
     
     ! ... knots
     self % knots ( 1 : self % p + 1 ) = 0.0
     do li_i = 1, self % n - self % p - 1
        self % knots ( self % p + 1 + li_i ) = li_i * 1.0 / ( self % n - self % p )
     end do 
     self % knots ( self % n + 1 : self % n + self % p + 1 ) = 1.0
                                     
     if ( ai_type_bc == SPI_BC_PERIODIC ) then
        li_nu = self % p
        
        lpr_knot ( : ) = self % knots ( : )
        call convert_to_periodic_knots ( lpr_Knot ( : ), self % n, self % p, li_nu, self % knots ( : ) )
     end if
     ! ...

  END SUBROUTINE INITIALIZE_MESH_KNOTS_DEFAULT_BSPLINES_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE INITIALIZE_MESH_1D_BSPLINE_CONTROL_POINTS(self, control_points_1d)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
      real(SPI_RK), dimension (:,:), INTENT(IN):: control_points_1d
     ! LOCAL
     INTEGER :: n_dim 

     n_dim =  SIZE(control_points_1d, 2) 
     self % n_dim = n_dim

     ALLOCATE ( self % control_points( self % n, n_dim ) ) 

     ! ... control points
     self % control_points = control_points_1d 
  END SUBROUTINE INITIALIZE_MESH_1D_BSPLINE_CONTROL_POINTS
  ! .........................................................

  ! .........................................................
  SUBROUTINE INITIALIZE_MESH_1D_BSPLINE_CONTROL_POINTS_DEFAULT(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: li_i
     INTEGER :: n_dim 

     n_dim = 1 
     self % n_dim = n_dim
     ALLOCATE ( self % control_points( self % n, n_dim ) ) 

     ! ... control points
     self % control_points = 0.0
     do li_i = 1, self % n 
        self % control_points( li_i, 1 ) = (li_i - 1) * 1.0d0 / (self % n - 1)
     end do 
     ! ...
  END SUBROUTINE INITIALIZE_MESH_1D_BSPLINE_CONTROL_POINTS_DEFAULT
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_BSPLINES_1D(self, ai_n, ai_p, ai_type_bc, apr_knots, control_points)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     INTEGER              , INTENT(IN)    :: ai_n
     INTEGER              , INTENT(IN)    :: ai_p
     INTEGER              , OPTIONAL, INTENT(IN)    :: ai_type_bc
     real(SPI_RK), dimension (:), OPTIONAL, INTENT(IN) :: apr_knots
     real(SPI_RK), dimension (:,:), OPTIONAL, INTENT(IN) :: control_points
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

     IF ( ( PRESENT(control_points)) ) THEN
        CALL INITIALIZE_MESH_1D_BSPLINE_CONTROL_POINTS(self, control_points) 
     ELSE
        CALL INITIALIZE_MESH_1D_BSPLINE_CONTROL_POINTS_DEFAULT(self) 
     END IF

     CALL CREATE_MESH_GRID_BSPLINES_1D(self) 

  END SUBROUTINE CREATE_MESH_BSPLINES_1D
  ! .........................................................
  
  ! .........................................................
  SUBROUTINE FREE_MESH_BSPLINES_1D(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_1D_BSPLINE)      , INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: li_err 

     DEALLOCATE ( self % knots ) 

  END SUBROUTINE FREE_MESH_BSPLINES_1D
  ! .........................................................

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
END MODULE SPI_MESH_BSPLINE 
