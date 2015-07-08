!# -*- coding: utf8 -*-
MODULE SPI_BLACKBOX
  USE SPI_BLACKBOX_DEF
  USE SPI_QUADRATURES_DEF
!  USE SPI_BASIS_DEF

    CONTAINS

  ! .........................................................
  SUBROUTINE BLACKBOX_CREATE(self, ao_quad)
  IMPLICIT NONE
    TYPE(DEF_BLACKBOX_1D)  :: self
    CLASS(DEF_QUADRATURE_ABSTRACT)  :: ao_quad
    ! LOCAL
    INTEGER, PARAMETER :: N_DIM=1
    INTEGER :: li_err

    self % oi_n_points   = ao_quad % oi_n_points
    ! TODO : update subroutine's signature 
    self % oi_basis_type = SPI_BASIS_BSPLINES 

    ALLOCATE(self%B_0    (ao_quad % oi_n_points))

    ALLOCATE(self%B_s1   (ao_quad % oi_n_points))
    ALLOCATE(self%B_s1s1 (ao_quad % oi_n_points))

    ALLOCATE(self%B_x1   (ao_quad % oi_n_points))
    ALLOCATE(self%B_x1x1 (ao_quad % oi_n_points))

    ALLOCATE(self%Xp_0    (N_DIM, ao_quad % oi_n_points))

    ALLOCATE(self%Xp_s1   (N_DIM, ao_quad % oi_n_points))
    ALLOCATE(self%Xp_s1s1 (N_DIM, ao_quad % oi_n_points))

    ALLOCATE(self%Xp_x1   (N_DIM, ao_quad % oi_n_points))
    ALLOCATE(self%Xp_x1x1 (N_DIM, ao_quad % oi_n_points))

    ALLOCATE(self%Vol     (ao_quad % oi_n_points))
    ALLOCATE(self%wVol    (ao_quad % oi_n_points))

  END SUBROUTINE BLACKBOX_CREATE
  ! .........................................................

  ! .........................................................
  SUBROUTINE BLACKBOX_RESET_POSITION(self)
  IMPLICIT NONE
    TYPE(DEF_BLACKBOX_1D)  :: self
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
  SUBROUTINE BLACKBOX_FREE(self)
  IMPLICIT NONE
    TYPE(DEF_BLACKBOX_1D)  :: self
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

  END SUBROUTINE BLACKBOX_FREE
  ! .........................................................

  ! .........................................................
  SUBROUTINE BLACKBOX_COMPUTE_METRIC(self, ar_a, ar_b, Nutor)
  ! ... TODO subroutine signature to change with the 1D elements object
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D)  :: self
     REAL(SPI_RK) :: ar_a
     REAL(SPI_RK) :: ar_b
     INTEGER, DIMENSION(2) :: Nutor
     ! LOCAL
     INTEGER :: kg
     INTEGER :: li_err


    DO kg = 1, self % oi_n_points
       self % Vol(kg) = ABS( ar_b - ar_a )
      
       SELECT CASE(self % oi_basis_type)
          CASE(SPI_BASIS_FOURIER) 
             Nutor(1) = 1 
             Nutor(2) = 1
             self % Vol(kg) = 1.0
       END SELECT

       self % wVol(kg) = self % Vol(kg) * self % ptr_quad % opr_weights(kg)
    END DO
  END SUBROUTINE BLACKBOX_COMPUTE_METRIC 
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
  SUBROUTINE BLACKBOX_GET_POSITION(self, ar_a, ar_b)
   IMPLICIT NONE
    TYPE(DEF_BLACKBOX_1D), INTENT(INOUT) :: self
    REAL(SPI_RK) :: ar_a
    REAL(SPI_RK) :: ar_b
    ! LOCAL
    INTEGER  :: li_err

    SELECT CASE(self % oi_basis_type)
    CASE(SPI_BASIS_BSPLINES) 
       ! TODO to check
       self % Xp_0(1,:)   = ar_a + self % ptr_quad % opr_points(1,:) * (ar_b - ar_a )
    CASE(SPI_BASIS_HBEZIER) 
       self % Xp_0(1,:)   = ar_a + self % ptr_quad % opr_points(1,:) * (ar_b - ar_a )
    CASE(SPI_BASIS_FOURIER) 
       self % Xp_0(1,:)   =  self % ptr_quad % opr_points(1,:)
    END SELECT

  END SUBROUTINE BLACKBOX_GET_POSITION
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
