!# -*- coding: utf8 -*-
MODULE SPI_MESH_LOGICAL_BSPLINE 
  USE SPI_MESH_LOGICAL_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  IMPLICIT NONE

  PRIVATE
  PUBLIC  ::   CREATE_MESH_LOGICAL_BSPLINE &
           & , FREE_MESH_LOGICAL_BSPLINE

CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_MESH_LOGICAL_BSPLINE(self, n, p, type_bc, knots)
  IMPLICIT NONE
     TYPE(DEF_MESH_LOGICAL_BSPLINE)  , INTENT(INOUT) :: self
     INTEGER              , INTENT(IN)    :: n
     INTEGER              , INTENT(IN)    :: p
     INTEGER              , OPTIONAL, INTENT(IN)    :: type_bc
     real(SPI_RK), dimension (:), OPTIONAL, INTENT(IN) :: knots
     ! LOCAL
     INTEGER :: li_err 

     IF ( (.NOT. PRESENT(knots)) &
            & .AND. (PRESENT(type_bc)) &
            & ) THEN
        CALL INITIALIZE_MESH_LOGICAL_BSPLINE_KNOTS_DEFAULT(self, n, p, type_bc) 
     END IF

     IF ( ( PRESENT(knots)) ) THEN
        CALL INITIALIZE_MESH_LOGICAL_BSPLINE_KNOTS(self, n, p, knots) 
     END IF

     CALL CREATE_MESH_LOGICAL_BSPLINE_GRID(self) 

  END SUBROUTINE CREATE_MESH_LOGICAL_BSPLINE
  ! .........................................................

  ! .........................................................
  SUBROUTINE FREE_MESH_LOGICAL_BSPLINE(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_LOGICAL_BSPLINE)      , INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: li_err 

     DEALLOCATE ( self % knots ) 
     DEALLOCATE ( self % grid) 
!     DEALLOCATE ( self % points) 

  END SUBROUTINE FREE_MESH_LOGICAL_BSPLINE
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_LOGICAL_BSPLINE_GRID(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_LOGICAL_BSPLINE)      , INTENT(INOUT) :: self
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

     ALLOCATE ( self % grid( n_elements + 1 ) ) 
     self % grid (1:n_elements+1) = lpr_grid(1:n_elements+1)

  END SUBROUTINE CREATE_MESH_LOGICAL_BSPLINE_GRID
  ! .........................................................

  ! .........................................................
  SUBROUTINE INITIALIZE_MESH_LOGICAL_BSPLINE_KNOTS(self, ai_n, ai_p, knots)
  IMPLICIT NONE
     TYPE(DEF_MESH_LOGICAL_BSPLINE)      , INTENT(INOUT) :: self
     INTEGER, INTENT(IN)    :: ai_n
     INTEGER, INTENT(IN)    :: ai_p
      real(SPI_RK), dimension (:), INTENT(IN):: knots
     ! LOCAL

     self % n = ai_n
     self % p = ai_p

     ALLOCATE ( self % knots( ai_n + ai_p + 1 ) ) 
     self % knots = knots

  END SUBROUTINE INITIALIZE_MESH_LOGICAL_BSPLINE_KNOTS
  ! .........................................................

  ! .........................................................
  SUBROUTINE INITIALIZE_MESH_LOGICAL_BSPLINE_KNOTS_DEFAULT(self, ai_n, ai_p, type_bc)
  IMPLICIT NONE
     TYPE(DEF_MESH_LOGICAL_BSPLINE)      , INTENT(INOUT) :: self
     INTEGER, INTENT(IN)    :: ai_n
     INTEGER, INTENT(IN)    :: ai_p
     INTEGER, INTENT(IN)    :: type_bc
     ! LOCAL
     INTEGER :: li_err 
     integer  :: li_i
     integer  :: li_nu ! number of continuity condition for periodic vector knot		
     real(SPI_RK), dimension ( ai_n + ai_p + 1 ) :: lpr_knot

     self % type_bc = type_bc
     self % p       = ai_p
     self % n       = ai_n
     
     if ( type_bc == SPI_BC_PERIODIC ) then
        self % n =  self % n + self % p - 1
     end if
                                     
     self % n_elements = self % n - ai_p
     
     ALLOCATE ( self % knots ( self % n + ai_p + 1 ) ) 
     
     ! INITIALIZING VECTOR KNOTS
     if ( self % n < self % p + 1 ) then
        STOP "Error INIT_MESH_BSPLINE: you must have N >= p + 1"
     end if
     
     ! ... knots
     self % knots ( 1 : self % p + 1 ) = 0.0
     do li_i = 1, self % n - self % p - 1
        self % knots ( self % p + 1 + li_i ) = li_i * 1.0 / ( self % n - self % p )
     end do 
     self % knots ( self % n + 1 : self % n + self % p + 1 ) = 1.0
                                     
     if ( type_bc == SPI_BC_PERIODIC ) then
        li_nu = self % p
        
        lpr_knot ( : ) = self % knots ( : )
        call convert_to_periodic_knots ( lpr_Knot ( : ), self % n, self % p, li_nu, self % knots ( : ) )
     end if
     ! ...

  END SUBROUTINE INITIALIZE_MESH_LOGICAL_BSPLINE_KNOTS_DEFAULT
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
END MODULE SPI_MESH_LOGICAL_BSPLINE
