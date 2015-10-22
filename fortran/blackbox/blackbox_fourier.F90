!# -*- coding: utf8 -*-
MODULE SPI_BLACKBOX_FOURIER
  USE SPI_BLACKBOX_DEF
  USE SPI_QUADRATURES_DEF
  USE SPI_BASIS_DEF
  IMPLICIT NONE

    CONTAINS

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
  SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_FOURIER(self)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_FOURIER), INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: kg
     INTEGER :: li_err

  END SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_FOURIER
  ! .........................................................

END MODULE SPI_BLACKBOX_FOURIER
