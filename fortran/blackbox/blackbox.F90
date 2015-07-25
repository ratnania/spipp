!# -*- coding: utf8 -*-
MODULE SPI_BLACKBOX
  USE SPI_BLACKBOX_DEF
  USE SPI_QUADRATURES_DEF
  USE SPI_BASIS_DEF
  IMPLICIT NONE

    CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_BLACKBOX(self, ao_basis, ao_quad)
  IMPLICIT NONE
    CLASS(DEF_BLACKBOX_ABSTRACT), INTENT(INOUT) :: self
    CLASS(DEF_BASIS_ABSTRACT)   , TARGET        :: ao_basis
    CLASS(DEF_QUADRATURE_1D)    , TARGET        :: ao_quad
    ! LOCAL
    INTEGER, PARAMETER :: N_DIM=1
    INTEGER :: li_err

    self % ptr_quad => ao_quad

    CALL CREATE_BLACKBOX_ABSTRACT(self, ao_quad % oi_n_points)
    
    ! ...
    SELECT TYPE (self)
    CLASS IS (DEF_BLACKBOX_1D_BSPLINE)
       SELECT TYPE (ao_basis)
       CLASS IS (DEF_BASIS_1D_BSPLINE)
          CALL CREATE_BLACKBOX_1D_BSPLINE(self, ao_basis)
       CLASS DEFAULT
          STOP 'CREATE_BLACKBOX: unexpected type for ao_basis object!'
       END SELECT

    CLASS IS (DEF_BLACKBOX_1D_HBEZIER)
       SELECT TYPE (ao_basis)
       CLASS IS (DEF_BASIS_1D_HBEZIER)
          CALL CREATE_BLACKBOX_1D_HBEZIER(self, ao_basis)
       CLASS DEFAULT
          STOP 'CREATE_BLACKBOX: unexpected type for ao_basis object!'
       END SELECT

    CLASS IS (DEF_BLACKBOX_1D_FOURIER)
       SELECT TYPE (ao_basis)
       CLASS IS (DEF_BASIS_1D_FOURIER)
          CALL CREATE_BLACKBOX_1D_FOURIER(self, ao_basis)
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
  SUBROUTINE CREATE_BLACKBOX_ABSTRACT(self, ai_n_points)
  IMPLICIT NONE
    CLASS(DEF_BLACKBOX_ABSTRACT)             :: self
    INTEGER                    , INTENT(IN) :: ai_n_points
    ! LOCAL
    INTEGER, PARAMETER :: N_DIM=1
    INTEGER :: li_err

    self % oi_n_points   = ai_n_points

    ALLOCATE(self%B_0    (ai_n_points))

    ALLOCATE(self%B_s1   (ai_n_points))
    ALLOCATE(self%B_s1s1 (ai_n_points))

    ALLOCATE(self%B_x1   (ai_n_points))
    ALLOCATE(self%B_x1x1 (ai_n_points))

    ALLOCATE(self%Xp_0    (N_DIM, ai_n_points))

    ALLOCATE(self%Xp_s1   (N_DIM, ai_n_points))
    ALLOCATE(self%Xp_s1s1 (N_DIM, ai_n_points))

    ALLOCATE(self%Xp_x1   (N_DIM, ai_n_points))
    ALLOCATE(self%Xp_x1x1 (N_DIM, ai_n_points))

    ALLOCATE(self%Vol     (ai_n_points))
    ALLOCATE(self%wVol    (ai_n_points))

  END SUBROUTINE CREATE_BLACKBOX_ABSTRACT
  ! .........................................................

  ! ........................................................
  SUBROUTINE CREATE_BLACKBOX_1D_BSPLINE(self, ao_basis)
  implicit none
     type(DEF_BLACKBOX_1D_BSPLINE)        , INTENT(INOUT) :: self
     TYPE(DEF_BASIS_1D_BSPLINE)   , TARGET, INTENT(IN)    :: ao_basis
     ! LOCAL VARIABLES
     integer  :: li_ref

     self % ptr_basis => ao_basis

  END SUBROUTINE CREATE_BLACKBOX_1D_BSPLINE
  ! ........................................................

  ! .........................................................
  SUBROUTINE CREATE_BLACKBOX_1D_FOURIER(self, ao_basis)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_FOURIER)        , INTENT(INOUT)  :: self
     TYPE(DEF_BASIS_1D_FOURIER)   , TARGET, INTENT(IN)     :: ao_basis
     ! LOCAL
     INTEGER :: li_err 

     self % ptr_basis => ao_basis

  END SUBROUTINE CREATE_BLACKBOX_1D_FOURIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_BLACKBOX_1D_HBEZIER(self, ao_basis)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_HBEZIER)        , INTENT(INOUT)  :: self
     TYPE(DEF_BASIS_1D_HBEZIER)   , TARGET, INTENT(IN)     :: ao_basis
     ! LOCAL
     INTEGER :: li_err 

     self % ptr_basis => ao_basis

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
  SUBROUTINE BLACKBOX_COMPUTE_METRIC(self, ar_a, ar_b, Nutor)
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
     CLASS IS (DEF_BLACKBOX_1D_HBEZIER)
        CALL COMPUTE_METRIC_BLACKBOX_1D_HBEZIER(self, ar_a, ar_b)
     CLASS IS (DEF_BLACKBOX_1D_FOURIER)
        IF (PRESENT(Nutor)) THEN 
           CALL COMPUTE_METRIC_BLACKBOX_1D_FOURIER(self, ar_a, ar_b, Nutor)
        ELSE
           STOP "BLACKBOX_COMPUTE_METRIC: missing argument"
        END IF
     CLASS DEFAULT
        STOP 'COMPUTE_METRIC_BLACKBOX: unexpected type for self object!'
     END SELECT
     ! ...

  END SUBROUTINE BLACKBOX_COMPUTE_METRIC 
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

     DO kg = 1, self % oi_n_points
        self % Vol(kg) = ABS( ar_b - ar_a )
       
        self % wVol(kg) = self % Vol(kg) !* self % ptr_quad % opr_weights(kg)
     END DO

  END SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_BSPLINE
  ! ........................................................

  ! .........................................................
  SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_FOURIER(self, ar_a, ar_b, Nutor)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_FOURIER), INTENT(INOUT) :: self
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_a
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_b
     INTEGER, DIMENSION(2)        , INTENT(INOUT) :: Nutor
     ! LOCAL
     INTEGER :: kg
     INTEGER :: li_err

     DO kg = 1, self % oi_n_points
        self % Vol(kg) = ABS( ar_b - ar_a )
       
        Nutor(1) = 1 
        Nutor(2) = 1
        self % Vol(kg) = 1.0
    
        self % wVol(kg) = self % Vol(kg) * self % ptr_quad % opr_weights(kg)
     END DO

  END SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_FOURIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_HBEZIER(self, ar_a, ar_b)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_HBEZIER), INTENT(INOUT)  :: self
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_a
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_b
     ! LOCAL
     INTEGER :: kg
     INTEGER :: li_err

     DO kg = 1, self % oi_n_points
        self % Vol(kg) = ABS( ar_b - ar_a )
        self % wVol(kg) = self % Vol(kg) * self % ptr_quad % opr_weights(kg)
     END DO

  END SUBROUTINE COMPUTE_METRIC_BLACKBOX_1D_HBEZIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE BLACKBOX_UPDATE_PHYSICAL_BASIS(self)
   IMPLICIT NONE
    TYPE(DEF_BLACKBOX_1D), INTENT(INOUT) :: self

    self % B_x1   = self % B_s1   /  self % Vol
    self % B_x1x1 = self % B_s1s1 / (self % Vol**2)

  END SUBROUTINE BLACKBOX_UPDATE_PHYSICAL_BASIS 
  ! .........................................................

  ! .........................................................
  SUBROUTINE BLACKBOX_UPDATE_POSITION()
   IMPLICIT NONE

  END SUBROUTINE BLACKBOX_UPDATE_POSITION 
  ! .........................................................

  ! .........................................................
  SUBROUTINE BLACKBOX_COMPUTE_POSITION(self, ar_a, ar_b)
  IMPLICIT NONE
     CLASS(DEF_BLACKBOX_ABSTRACT)    , INTENT(INOUT) :: self
     REAL(SPI_RK)                   , INTENT(IN)    :: ar_a
     REAL(SPI_RK)                   , INTENT(IN)    :: ar_b
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_BLACKBOX_1D_BSPLINE)
        CALL COMPUTE_POSITION_BLACKBOX_1D_BSPLINE(self, ar_a, ar_b)
     CLASS IS (DEF_BLACKBOX_1D_HBEZIER)
        CALL COMPUTE_POSITION_BLACKBOX_1D_HBEZIER(self, ar_a, ar_b)
     CLASS IS (DEF_BLACKBOX_1D_FOURIER)
        CALL COMPUTE_POSITION_BLACKBOX_1D_FOURIER(self, ar_a, ar_b)
     CLASS DEFAULT
        STOP 'COMPUTE_POSITION_BLACKBOX: unexpected type for self object!'
     END SELECT
     ! ...

  END SUBROUTINE BLACKBOX_COMPUTE_POSITION 
  ! .........................................................

  ! ........................................................
  SUBROUTINE COMPUTE_POSITION_BLACKBOX_1D_BSPLINE(self, ar_a, ar_b)
  implicit none
     type(DEF_BLACKBOX_1D_BSPLINE), INTENT(INOUT) :: self
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_a
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_b
     ! LOCAL VARIABLES
     INTEGER :: kg
     INTEGER :: li_err

     self % Xp_0(1,:)   = ar_a + self % ptr_quad % opr_points(1,:) * (ar_b - ar_a )

  END SUBROUTINE COMPUTE_POSITION_BLACKBOX_1D_BSPLINE
  ! ........................................................

  ! .........................................................
  SUBROUTINE COMPUTE_POSITION_BLACKBOX_1D_FOURIER(self, ar_a, ar_b)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_FOURIER), INTENT(INOUT) :: self
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_a
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_b
     ! LOCAL
     INTEGER :: kg
     INTEGER :: li_err

     self % Xp_0(1,:)   =  self % ptr_quad % opr_points(1,:)

  END SUBROUTINE COMPUTE_POSITION_BLACKBOX_1D_FOURIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE COMPUTE_POSITION_BLACKBOX_1D_HBEZIER(self, ar_a, ar_b)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_HBEZIER), INTENT(INOUT)  :: self
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_a
     REAL(SPI_RK)                 , INTENT(IN)    :: ar_b
     ! LOCAL
     INTEGER :: kg
     INTEGER :: li_err

     self % Xp_0(1,:)   = ar_a + self % ptr_quad % opr_points(1,:) * (ar_b - ar_a )

  END SUBROUTINE COMPUTE_POSITION_BLACKBOX_1D_HBEZIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE BLACKBOX_GET_BASIS(self, or_tor1, iv_tor1)
   IMPLICIT NONE
    TYPE(DEF_BLACKBOX_1D), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: or_tor1
    INTEGER, INTENT(IN) :: iv_tor1
    ! LOCAL

!    self % B_0    (:) = self % ptr_basis % TestFT_0 (:, or_tor1, iv_tor1)
!    self % B_s1   (:) = self % ptr_basis % TestFT_p (:, or_tor1, iv_tor1)
!    self % B_s1s1 (:) = self % ptr_basis % TestFT_pp(:, or_tor1, iv_tor1)

  END SUBROUTINE BLACKBOX_GET_BASIS 
  ! .........................................................

END MODULE SPI_BLACKBOX
