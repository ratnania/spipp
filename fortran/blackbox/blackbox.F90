!# -*- coding: utf8 -*-
MODULE SPI_BLACKBOX
  USE SPI_BLACKBOX_DEF
  USE SPI_QUADRATURES_DEF
  USE SPI_BASIS_DEF
  IMPLICIT NONE

    CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_BLACKBOX(self, ao_mesh, ao_basis, ao_quad)
  IMPLICIT NONE
    CLASS(DEF_BLACKBOX_ABSTRACT), INTENT(INOUT) :: self
    CLASS(DEF_MESH_ABSTRACT)    , TARGET        :: ao_mesh
    CLASS(DEF_BASIS_ABSTRACT)   , TARGET        :: ao_basis
    CLASS(DEF_QUADRATURE_1D)    , TARGET        :: ao_quad
    ! LOCAL
    INTEGER, PARAMETER :: N_DIM=1
    INTEGER :: li_err

    self % ptr_quad => ao_quad

    CALL CREATE_BLACKBOX_ABSTRACT(self, ao_mesh % n_elements, ao_quad % oi_n_points)
    
    ! ...
    SELECT TYPE (self)
    CLASS IS (DEF_BLACKBOX_1D_BSPLINE)
       SELECT TYPE (ao_basis)
       CLASS IS (DEF_BASIS_1D_BSPLINE)
         SELECT TYPE (ao_mesh)
         CLASS IS (DEF_MESH_1D_BSPLINE)
            CALL CREATE_BLACKBOX_1D_BSPLINE(self, ao_mesh, ao_basis)
         CLASS DEFAULT
            STOP 'CREATE_BLACKBOX: unexpected type for ao_mesh object!'
         END SELECT
       CLASS DEFAULT
          STOP 'CREATE_BLACKBOX: unexpected type for ao_basis object!'
       END SELECT

    CLASS IS (DEF_BLACKBOX_1D_HBEZIER)
       SELECT TYPE (ao_basis)
       CLASS IS (DEF_BASIS_1D_HBEZIER)
         SELECT TYPE (ao_mesh)
         CLASS IS (DEF_MESH_1D_HBEZIER)
            CALL CREATE_BLACKBOX_1D_HBEZIER(self, ao_mesh, ao_basis)
         CLASS DEFAULT
            STOP 'CREATE_BLACKBOX: unexpected type for ao_mesh object!'
         END SELECT
       CLASS DEFAULT
          STOP 'CREATE_BLACKBOX: unexpected type for ao_basis object!'
       END SELECT

    CLASS IS (DEF_BLACKBOX_1D_FOURIER)
       SELECT TYPE (ao_basis)
       CLASS IS (DEF_BASIS_1D_FOURIER)
         SELECT TYPE (ao_mesh)
         CLASS IS (DEF_MESH_1D_FOURIER)
            CALL CREATE_BLACKBOX_1D_FOURIER(self, ao_mesh, ao_basis)
         CLASS DEFAULT
            STOP 'CREATE_BLACKBOX: unexpected type for ao_mesh object!'
         END SELECT
       CLASS DEFAULT
          STOP 'CREATE_BLACKBOX: unexpected type for ao_basis object!'
       END SELECT
    CLASS DEFAULT
       STOP 'CREATE_BLACKBOX: unexpected type for self object!'
    END SELECT
    ! ...

  END SUBROUTINE CREATE_BLACKBOX
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_BLACKBOX_ABSTRACT(self, ai_n_elements, ai_n_points)
  IMPLICIT NONE
    CLASS(DEF_BLACKBOX_ABSTRACT)             :: self
    INTEGER                    , INTENT(IN) :: ai_n_elements
    INTEGER                    , INTENT(IN) :: ai_n_points
    ! LOCAL
    INTEGER, PARAMETER :: N_DIM=1
    INTEGER :: li_err

    self % oi_n_points   = ai_n_points
    self % oi_n_elements = ai_n_elements 

    ALLOCATE(self%B_0    (ai_n_elements, ai_n_points))

    ALLOCATE(self%B_s1   (ai_n_elements, ai_n_points))
    ALLOCATE(self%B_s1s1 (ai_n_elements, ai_n_points))

    ALLOCATE(self%B_x1   (ai_n_elements, ai_n_points))
    ALLOCATE(self%B_x1x1 (ai_n_elements, ai_n_points))

    ALLOCATE(self%Xp_0    (ai_n_elements, ai_n_points, N_DIM))
                                                            
    ALLOCATE(self%Xp_s1   (ai_n_elements, ai_n_points, N_DIM))
    ALLOCATE(self%Xp_s1s1 (ai_n_elements, ai_n_points, N_DIM))
                                                            
    ALLOCATE(self%Xp_x1   (ai_n_elements, ai_n_points, N_DIM))
    ALLOCATE(self%Xp_x1x1 (ai_n_elements, ai_n_points, N_DIM))

    ALLOCATE(self%Vol     (ai_n_elements, ai_n_points))
    ALLOCATE(self%wVol    (ai_n_elements, ai_n_points))

  END SUBROUTINE CREATE_BLACKBOX_ABSTRACT
  ! .........................................................

  ! ........................................................
  SUBROUTINE CREATE_BLACKBOX_1D_BSPLINE(self, ao_mesh, ao_basis)
  implicit none
     type(DEF_BLACKBOX_1D_BSPLINE)        , INTENT(INOUT) :: self
     TYPE(DEF_MESH_1D_BSPLINE)   , TARGET, INTENT(IN)    :: ao_mesh
     TYPE(DEF_BASIS_1D_BSPLINE)   , TARGET, INTENT(IN)    :: ao_basis
     ! LOCAL VARIABLES
     integer  :: li_ref

     self % ptr_basis => ao_basis
     self % ptr_mesh => ao_mesh

  END SUBROUTINE CREATE_BLACKBOX_1D_BSPLINE
  ! ........................................................

  ! .........................................................
  SUBROUTINE CREATE_BLACKBOX_1D_FOURIER(self, ao_mesh, ao_basis)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_FOURIER)        , INTENT(INOUT)  :: self
     TYPE(DEF_MESH_1D_FOURIER)   , TARGET, INTENT(IN)     :: ao_mesh
     TYPE(DEF_BASIS_1D_FOURIER)   , TARGET, INTENT(IN)     :: ao_basis
     ! LOCAL
     INTEGER :: li_err 

     self % ptr_basis => ao_basis
     self % ptr_mesh => ao_mesh

  END SUBROUTINE CREATE_BLACKBOX_1D_FOURIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_BLACKBOX_1D_HBEZIER(self, ao_mesh, ao_basis)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_HBEZIER)        , INTENT(INOUT)  :: self
     TYPE(DEF_MESH_1D_HBEZIER)   , TARGET, INTENT(IN)     :: ao_mesh
     TYPE(DEF_BASIS_1D_HBEZIER)   , TARGET, INTENT(IN)     :: ao_basis
     ! LOCAL
     INTEGER :: li_err 

     self % ptr_basis => ao_basis
     self % ptr_mesh => ao_mesh

  END SUBROUTINE CREATE_BLACKBOX_1D_HBEZIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE FREE_BLACKBOX(self)
  IMPLICIT NONE
    CLASS(DEF_BLACKBOX_ABSTRACT)   :: self
    ! LOCAL
    INTEGER :: li_err

    CALL FREE_BLACKBOX_ABSTRACT(self)
    
  END SUBROUTINE FREE_BLACKBOX
  ! .........................................................

  ! .........................................................
  SUBROUTINE FREE_BLACKBOX_ABSTRACT(self)
  IMPLICIT NONE
    CLASS(DEF_BLACKBOX_ABSTRACT)  :: self
    ! LOCAL

    DEALLOCATE(self%B_0 )

    DEALLOCATE(self%B_s1 )
    DEALLOCATE(self%B_s1s1)

    DEALLOCATE(self%B_x1 )
    DEALLOCATE(self%B_x1x1)

    DEALLOCATE(self%Xp_0)

    DEALLOCATE(self%Xp_s1)
    DEALLOCATE(self%Xp_s1s1)

    DEALLOCATE(self%Xp_x1)
    DEALLOCATE(self%Xp_x1x1)

    DEALLOCATE(self%Vol)
    DEALLOCATE(self%wVol)

  END SUBROUTINE FREE_BLACKBOX_ABSTRACT
  ! .........................................................

  ! .........................................................
  SUBROUTINE BLACKBOX_RESET_POSITION(self)
  IMPLICIT NONE
    CLASS(DEF_BLACKBOX_ABSTRACT)  :: self
    ! LOCAL

    self%Xp_0  = 0.0 

    self%Xp_s1   = 0.0
    self%Xp_s1s1 = 0.0 
                  
    self%Xp_x1   = 0.0 
    self%Xp_x1x1 = 0.0 

    self%Vol  = 0.0
    self%wVol = 0.0

  END SUBROUTINE BLACKBOX_RESET_POSITION
  ! .........................................................

  ! .........................................................
  SUBROUTINE COMPUTE_METRIC_BLACKBOX(self, ar_a, ar_b, Nutor)
  ! ... TODO subroutine signature to change with the 1D elements object
  IMPLICIT NONE
     CLASS(DEF_BLACKBOX_ABSTRACT)    , INTENT(INOUT) :: self
     REAL(SPI_RK)                   , INTENT(IN)    :: ar_a
     REAL(SPI_RK)                   , INTENT(IN)    :: ar_b
     INTEGER, DIMENSION(2), OPTIONAL, INTENT(INOUT) :: Nutor
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_BLACKBOX_1D_BSPLINE)
        CALL COMPUTE_METRIC_BLACKBOX_1D_BSPLINE(self, ar_a, ar_b)
!     CLASS IS (DEF_BLACKBOX_1D_HBEZIER)
!        CALL COMPUTE_METRIC_BLACKBOX_1D_HBEZIER(self, ar_a, ar_b)
!     CLASS IS (DEF_BLACKBOX_1D_FOURIER)
!        IF (PRESENT(Nutor)) THEN 
!           CALL COMPUTE_METRIC_BLACKBOX_1D_FOURIER(self, ar_a, ar_b, Nutor)
!        ELSE
!           STOP "COMPUTE_METRIC_BLACKBOX: missing argument"
!        END IF
     CLASS DEFAULT
        STOP 'COMPUTE_METRIC_BLACKBOX: unexpected type for self object!'
     END SELECT
     ! ...

  END SUBROUTINE COMPUTE_METRIC_BLACKBOX 
  ! .........................................................

  ! ........................................................
  SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_BSPLINE(self, ar_a, ar_b)
  implicit none
     type(DEF_BLACKBOX_1D_BSPLINE), INTENT(INOUT) :: self
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_a
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_b
     ! LOCAL VARIABLES
     INTEGER :: kg
     INTEGER :: li_err
     ! ... TODO add loop over elements
     INTEGER :: i_element

     DO kg = 1, self % oi_n_points
        self % Vol(i_element, kg) = ABS( ar_b - ar_a )
       
        self % wVol(i_element, kg) = self % Vol(i_element, kg) !* self % ptr_quad % opr_weights(kg)
     END DO

  END SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_BSPLINE
  ! ........................................................

!  ! .........................................................
!  SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_FOURIER(self, ar_a, ar_b, Nutor)
!  IMPLICIT NONE
!     TYPE(DEF_BLACKBOX_1D_FOURIER), INTENT(INOUT) :: self
!     REAL(SPI_RK)                 , INTENT(IN)    :: ar_a
!     REAL(SPI_RK)                 , INTENT(IN)    :: ar_b
!     INTEGER, DIMENSION(2)        , INTENT(INOUT) :: Nutor
!     ! LOCAL
!     INTEGER :: kg
!     INTEGER :: li_err
!
!     DO kg = 1, self % oi_n_points
!        self % Vol(kg) = ABS( ar_b - ar_a )
!       
!        Nutor(1) = 1 
!        Nutor(2) = 1
!        self % Vol(kg) = 1.0
!    
!        self % wVol(kg) = self % Vol(kg) * self % ptr_quad % opr_weights(kg)
!     END DO
!
!  END SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_FOURIER
!  ! .........................................................
!
!  ! .........................................................
!  SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_HBEZIER(self, ar_a, ar_b)
!  IMPLICIT NONE
!     TYPE(DEF_BLACKBOX_1D_HBEZIER), INTENT(INOUT)  :: self
!     REAL(SPI_RK)                 , INTENT(IN)    :: ar_a
!     REAL(SPI_RK)                 , INTENT(IN)    :: ar_b
!     ! LOCAL
!     INTEGER :: kg
!     INTEGER :: li_err
!
!     DO kg = 1, self % oi_n_points
!        self % Vol(kg) = ABS( ar_b - ar_a )
!        self % wVol(kg) = self % Vol(kg) * self % ptr_quad % opr_weights(kg)
!     END DO
!
!  END SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_HBEZIER
!  ! .........................................................

  ! .........................................................
  SUBROUTINE UPDATE_PHYSICAL_BASIS_BLACKBOX(self)
   IMPLICIT NONE
    CLASS(DEF_BLACKBOX_ABSTRACT), INTENT(INOUT) :: self

    self % B_x1   = self % B_s1   /  self % Vol
    self % B_x1x1 = self % B_s1s1 / (self % Vol**2)

  END SUBROUTINE UPDATE_PHYSICAL_BASIS_BLACKBOX 
  ! .........................................................

  ! .........................................................
  SUBROUTINE BLACKBOX_UPDATE_POSITION()
   IMPLICIT NONE

  END SUBROUTINE BLACKBOX_UPDATE_POSITION 
  ! .........................................................

  ! .........................................................
  SUBROUTINE UPDATE_LOGICAL_POSITION_BLACKBOX(self)
  IMPLICIT NONE
     CLASS(DEF_BLACKBOX_ABSTRACT)    , INTENT(INOUT) :: self
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_BLACKBOX_1D_BSPLINE)
        CALL UPDATE_LOGICAL_POSITION_BLACKBOX_1D_BSPLINE(self)
!     CLASS IS (DEF_BLACKBOX_1D_HBEZIER)
!        CALL UPDATE_LOGICAL_POSITION_BLACKBOX_1D_HBEZIER(self)
!     CLASS IS (DEF_BLACKBOX_1D_FOURIER)
!        CALL UPDATE_LOGICAL_POSITION_BLACKBOX_1D_FOURIER(self)
     CLASS DEFAULT
        STOP 'UPDATE_LOGICAL_POSITION_BLACKBOX: unexpected type for self object!'
     END SELECT
     ! ...

  END SUBROUTINE UPDATE_LOGICAL_POSITION_BLACKBOX 
  ! .........................................................

  ! ........................................................
  SUBROUTINE UPDATE_LOGICAL_POSITION_BLACKBOX_1D_BSPLINE(self)
  implicit none
     type(DEF_BLACKBOX_1D_BSPLINE), INTENT(INOUT) :: self
     ! LOCAL VARIABLES
     INTEGER :: kg
     INTEGER :: li_err
     INTEGER :: i_element
     INTEGER :: i
     INTEGER, PARAMETER :: N_DIM = 1 
     REAL(SPI_RK) :: lr_a
     REAL(SPI_RK) :: lr_b

     DO i =1, N_DIM
       self % Xp_0(i_element,:, i) = lr_a + self % ptr_quad % opr_points(:) * (lr_b - lr_a )
     END DO

  END SUBROUTINE UPDATE_LOGICAL_POSITION_BLACKBOX_1D_BSPLINE
  ! ........................................................

  ! .........................................................
  SUBROUTINE GET_BASIS_BLACKBOX(self, ai_loc)
  IMPLICIT NONE
     CLASS(DEF_BLACKBOX_ABSTRACT)    , INTENT(INOUT) :: self
     INTEGER, INTENT(IN) :: ai_loc
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_BLACKBOX_1D_BSPLINE)
        CALL GET_BASIS_BLACKBOX_1D_BSPLINE(self, ai_loc)
!     CLASS IS (DEF_BLACKBOX_1D_HBEZIER)
!        CALL GET_BASIS_BLACKBOX_1D_HBEZIER(self)
!     CLASS IS (DEF_BLACKBOX_1D_FOURIER)
!        CALL GET_BASIS_BLACKBOX_1D_FOURIER(self)
     CLASS DEFAULT
        STOP 'GET_BASIS_BLACKBOX: unexpected type for self object!'
     END SELECT
     ! ...

  END SUBROUTINE GET_BASIS_BLACKBOX 
  ! .........................................................

  ! ..........................................................        
  SUBROUTINE GET_BASIS_BLACKBOX_1D_BSPLINE(self, ai_loc)
  IMPLICIT NONE
     CLASS(DEF_BLACKBOX_1D_BSPLINE) :: self
     INTEGER, INTENT(IN) :: ai_loc
     ! LOCAL

!    self % B_0    (:) = self % ptr_basis % TestFT_0 (:, or_tor1, iv_tor1)
!    self % B_s1   (:) = self % ptr_basis % TestFT_p (:, or_tor1, iv_tor1)
!    self % B_s1s1 (:) = self % ptr_basis % TestFT_pp(:, or_tor1, iv_tor1)

  END SUBROUTINE GET_BASIS_BLACKBOX_1D_BSPLINE
  ! ..........................................................  

  ! .........................................................
  SUBROUTINE GET_BASIS_BLACKBOX_1D_FOURIER(self)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_FOURIER), INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: kg
     INTEGER :: li_err

  END SUBROUTINE GET_BASIS_BLACKBOX_1D_FOURIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE GET_BASIS_BLACKBOX_1D_HBEZIER(self)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_HBEZIER), INTENT(INOUT)  :: self
     ! LOCAL
     INTEGER :: kg
     INTEGER :: li_err

  END SUBROUTINE GET_BASIS_BLACKBOX_1D_HBEZIER
  ! .........................................................



  ! .........................................................
  SUBROUTINE UPDATE_POSITION_BLACKBOX(self, apr_control_point)
  IMPLICIT NONE
     CLASS(DEF_BLACKBOX_ABSTRACT)    , INTENT(INOUT) :: self
     REAL(SPI_RK), DIMENSION(:), OPTIONAL :: apr_control_point
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_BLACKBOX_1D_BSPLINE)
        IF (PRESENT(apr_control_point)) THEN 
           CALL UPDATE_POSITION_BLACKBOX_1D_BSPLINE(self, apr_control_point)
        ELSE
           STOP "UPDATE_POSITION_BLACKBOX: Wrong arguments"
        END IF
     CLASS IS (DEF_BLACKBOX_1D_HBEZIER)
        CALL UPDATE_POSITION_BLACKBOX_1D_HBEZIER(self)
     CLASS IS (DEF_BLACKBOX_1D_FOURIER)
        CALL UPDATE_POSITION_BLACKBOX_1D_FOURIER(self)
     CLASS DEFAULT
        STOP 'UPDATE_POSITION_BLACKBOX: unexpected type for self object!'
     END SELECT
     ! ...

  END SUBROUTINE UPDATE_POSITION_BLACKBOX 
  ! .........................................................

  ! ..........................................................        
  SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_BSPLINE(self, apr_control_point)
  IMPLICIT NONE
     CLASS(DEF_BLACKBOX_1D_BSPLINE) :: self
     REAL(SPI_RK), DIMENSION(:) :: apr_control_point
     INTEGER :: ai_i
     ! LOCAL
     INTEGER :: ijg
     INTEGER :: i_element

     ! .. TODO add loop over elements

    DO ijg = 1, self % ptr_quad % oi_n_points
       self % Xp_0   (i_element, ijg, :) = self % Xp_0   (i_element, ijg, :) &
               & + self % B_0(i_element, ijg)    * apr_control_point(:)
       self % Xp_s1  (i_element, ijg, :) = self % Xp_s1  (i_element, ijg, :) &
               & + self % B_s1(i_element, ijg)   * apr_control_point(:)
       self % Xp_s1s1(i_element, ijg, :) = self % Xp_s1s1(i_element, ijg, :) &
               & + self % B_s1s1(i_element, ijg) * apr_control_point(:)
    END DO

  END SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_BSPLINE
  ! ..........................................................  

  ! .........................................................
  SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_FOURIER(self)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_FOURIER), INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: kg
     INTEGER :: li_err

  END SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_FOURIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_HBEZIER(self)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_HBEZIER), INTENT(INOUT)  :: self
     ! LOCAL
     INTEGER :: kg
     INTEGER :: li_err

  END SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_HBEZIER
  ! .........................................................


END MODULE SPI_BLACKBOX
