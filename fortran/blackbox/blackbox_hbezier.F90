!# -*- coding: utf8 -*-
MODULE SPI_BLACKBOX_HBEZIER
  USE SPI_BLACKBOX_DEF
  USE SPI_QUADRATURES_DEF
  USE SPI_BASIS_DEF
  IMPLICIT NONE

    CONTAINS

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
  SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_HBEZIER(self)
  IMPLICIT NONE
     TYPE(DEF_BLACKBOX_1D_HBEZIER), INTENT(INOUT)  :: self
     ! LOCAL
     INTEGER :: kg
     INTEGER :: li_err

  END SUBROUTINE UPDATE_POSITION_BLACKBOX_1D_HBEZIER
  ! .........................................................


END MODULE SPI_BLACKBOX_HBEZIER
