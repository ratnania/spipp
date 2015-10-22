!# -*- coding: utf8 -*-
MODULE SPI_BLACKBOX
  USE SPI_BLACKBOX_DEF
  USE SPI_BLACKBOX_BSPLINE
!  USE SPI_BLACKBOX_FOURIER
!  USE SPI_BLACKBOX_HBEZIER
  USE SPI_QUADRATURES_DEF
  USE SPI_BASIS_DEF
  IMPLICIT NONE

    CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_BLACKBOX(self, ao_mesh, ao_basis, ao_quad, control_points_1d)
  IMPLICIT NONE
    CLASS(DEF_BLACKBOX_ABSTRACT), INTENT(INOUT) :: self
    CLASS(DEF_MESH_ABSTRACT)    , TARGET        :: ao_mesh
    CLASS(DEF_BASIS_ABSTRACT)   , TARGET        :: ao_basis
    CLASS(DEF_QUADRATURE_1D)    , TARGET        :: ao_quad
    real(SPI_RK), dimension (:,:), OPTIONAL, INTENT(IN) :: control_points_1d
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
            CALL CREATE_BLACKBOX_1D_BSPLINE(self, ao_mesh, ao_basis, control_points_1d=control_points_1d)
         CLASS DEFAULT
            STOP 'CREATE_BLACKBOX: unexpected type for ao_mesh object!'
         END SELECT
       CLASS DEFAULT
          STOP 'CREATE_BLACKBOX: unexpected type for ao_basis object!'
       END SELECT

!    CLASS IS (DEF_BLACKBOX_1D_HBEZIER)
!       SELECT TYPE (ao_basis)
!       CLASS IS (DEF_BASIS_1D_HBEZIER)
!         SELECT TYPE (ao_mesh)
!         CLASS IS (DEF_MESH_1D_HBEZIER)
!            CALL CREATE_BLACKBOX_1D_HBEZIER(self, ao_mesh, ao_basis)
!         CLASS DEFAULT
!            STOP 'CREATE_BLACKBOX: unexpected type for ao_mesh object!'
!         END SELECT
!       CLASS DEFAULT
!          STOP 'CREATE_BLACKBOX: unexpected type for ao_basis object!'
!       END SELECT
!
!    CLASS IS (DEF_BLACKBOX_1D_FOURIER)
!       SELECT TYPE (ao_basis)
!       CLASS IS (DEF_BASIS_1D_FOURIER)
!         SELECT TYPE (ao_mesh)
!         CLASS IS (DEF_MESH_1D_FOURIER)
!            CALL CREATE_BLACKBOX_1D_FOURIER(self, ao_mesh, ao_basis)
!         CLASS DEFAULT
!            STOP 'CREATE_BLACKBOX: unexpected type for ao_mesh object!'
!         END SELECT
!       CLASS DEFAULT
!          STOP 'CREATE_BLACKBOX: unexpected type for ao_basis object!'
!       END SELECT
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

    ALLOCATE(self%Xp_0    (ai_n_elements, ai_n_points, N_DIM))
                                                            
    ALLOCATE(self%Xp_s1   (ai_n_elements, ai_n_points, N_DIM))
    ALLOCATE(self%Xp_s1s1 (ai_n_elements, ai_n_points, N_DIM))
                                                            
    ALLOCATE(self%Xp_x1   (ai_n_elements, ai_n_points, N_DIM))
    ALLOCATE(self%Xp_x1x1 (ai_n_elements, ai_n_points, N_DIM))

    ALLOCATE(self%Jacobians(ai_n_elements, ai_n_points))
    ALLOCATE(self%Vol     (ai_n_elements, ai_n_points))
    ALLOCATE(self%wVol    (ai_n_elements, ai_n_points))

  END SUBROUTINE CREATE_BLACKBOX_ABSTRACT
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
  SUBROUTINE COMPUTE_METRIC_BLACKBOX(self)
  ! ... TODO subroutine signature to change with the 1D elements object
  IMPLICIT NONE
     CLASS(DEF_BLACKBOX_ABSTRACT)    , INTENT(INOUT) :: self
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_BLACKBOX_1D_BSPLINE)
        CALL COMPUTE_METRIC_BLACKBOX_1D_BSPLINE(self)
!     CLASS IS (DEF_BLACKBOX_1D_HBEZIER)
!        CALL COMPUTE_METRIC_BLACKBOX_1D_HBEZIER(self)
!     CLASS IS (DEF_BLACKBOX_1D_FOURIER)
!           CALL COMPUTE_METRIC_BLACKBOX_1D_FOURIER(self)
     CLASS DEFAULT
        STOP 'COMPUTE_METRIC_BLACKBOX: unexpected type for self object!'
     END SELECT
     ! ...

  END SUBROUTINE COMPUTE_METRIC_BLACKBOX 
  ! .........................................................

  ! .........................................................
  SUBROUTINE UPDATE_PHYSICAL_BASIS_BLACKBOX(self)
   IMPLICIT NONE
    CLASS(DEF_BLACKBOX_ABSTRACT), INTENT(INOUT) :: self

!    self % B_x1   = self % B_s1   /  self % Vol
!    self % B_x1x1 = self % B_s1s1 / (self % Vol**2)

  END SUBROUTINE UPDATE_PHYSICAL_BASIS_BLACKBOX 
  ! .........................................................

  ! .........................................................
  SUBROUTINE BLACKBOX_UPDATE_POSITION()
   IMPLICIT NONE

  END SUBROUTINE BLACKBOX_UPDATE_POSITION 
  ! .........................................................

  ! .........................................................
  SUBROUTINE UPDATE_POSITION_BLACKBOX(self)
  IMPLICIT NONE
     CLASS(DEF_BLACKBOX_ABSTRACT)    , INTENT(INOUT) :: self
     ! LOCAL

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_BLACKBOX_1D_BSPLINE)
        CALL UPDATE_POSITION_BLACKBOX_1D_BSPLINE(self)
!     CLASS IS (DEF_BLACKBOX_1D_HBEZIER)
!        CALL UPDATE_POSITION_BLACKBOX_1D_HBEZIER(self)
!     CLASS IS (DEF_BLACKBOX_1D_FOURIER)
!        CALL UPDATE_POSITION_BLACKBOX_1D_FOURIER(self)
     CLASS DEFAULT
        STOP 'UPDATE_POSITION_BLACKBOX: unexpected type for self object!'
     END SELECT
     ! ...

  END SUBROUTINE UPDATE_POSITION_BLACKBOX 
  ! .........................................................

END MODULE SPI_BLACKBOX
