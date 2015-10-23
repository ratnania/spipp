!# -*- coding: utf8 -*-
MODULE SPI_MESH_BSPLINE 
  USE SPI_MESH_DEF
  USE SPI_MESH_LOGICAL_DEF
  USE SPI_MESH_LOGICAL
  USE SPI_GLOBAL_DEF
  USE SPI_QUADRATURES_DEF
  IMPLICIT NONE

CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_MESH_BSPLINE(   self              &
                                & , d_dim             &
                                & , logical_mesh_s1   &
                                & , logical_mesh_s2   &
                                & , logical_mesh_s3   &
                                &  )
  IMPLICIT NONE
     CLASS(DEF_MESH_BSPLINE_ABSTRACT)                , INTENT(INOUT) :: self
     INTEGER                               , OPTIONAL, INTENT(IN)    :: d_dim
     TYPE(DEF_MESH_LOGICAL_BSPLINE), TARGET, OPTIONAL, INTENT(IN)    :: logical_mesh_s1 
     TYPE(DEF_MESH_LOGICAL_BSPLINE), TARGET, OPTIONAL, INTENT(IN)    :: logical_mesh_s2
     TYPE(DEF_MESH_LOGICAL_BSPLINE), TARGET, OPTIONAL, INTENT(IN)    :: logical_mesh_s3
     ! LOCAL
     INTEGER :: n_dim

     ! ... Manifold dimension
     self % d_dim = 1
     IF ( ( PRESENT(d_dim)) ) THEN
        self % d_dim = d_dim 
     END IF
     ! ...

     ! ... parametric dimension
     n_dim = 0
     IF ( ( PRESENT(logical_mesh_s1)) ) THEN
        n_dim = n_dim + 1
     END IF
     IF ( ( PRESENT(logical_mesh_s2)) ) THEN
        n_dim = n_dim + 1
     END IF
     IF ( ( PRESENT(logical_mesh_s3)) ) THEN
        n_dim = n_dim + 1
     END IF
     self % n_dim = n_dim
     ! ...

     ! ... allocate the array of logical meshes and set the pointers
     ALLOCATE(self % logicals(n_dim))

     IF ( ( PRESENT(logical_mesh_s1)) ) THEN
        self % logicals(1) % ptr_mesh => logical_mesh_s1
     END IF
     IF ( ( PRESENT(logical_mesh_s2)) ) THEN
        self % logicals(2) % ptr_mesh => logical_mesh_s2
     END IF
     IF ( ( PRESENT(logical_mesh_s3)) ) THEN
        self % logicals(3) % ptr_mesh => logical_mesh_s3
     END IF
     ! ...

     ! ...
     SELECT TYPE (self)
     CLASS IS (DEF_MESH_BSPLINE_1D)
        CALL CREATE_MESH_BSPLINE_1D(self)
     CLASS IS (DEF_MESH_BSPLINE_2D)
        CALL CREATE_MESH_BSPLINE_2D(self)
     CLASS IS (DEF_MESH_BSPLINE_3D)
        CALL CREATE_MESH_BSPLINE_3D(self)
     CLASS IS (DEF_MESH_BSPLINE_1D1D)
        CALL CREATE_MESH_BSPLINE_1D1D(self)
     CLASS DEFAULT
        STOP 'CREATE_MESH_BSPLINE: unexpected type for self object!'
     END SELECT
     ! ...

  END SUBROUTINE CREATE_MESH_BSPLINE
  ! .........................................................

  ! .........................................................
  SUBROUTINE FREE_MESH_BSPLINE(self)
  IMPLICIT NONE
     CLASS(DEF_MESH_BSPLINE_ABSTRACT)      , INTENT(INOUT) :: self
     ! LOCAL

     DEALLOCATE(self % logicals)

  END SUBROUTINE FREE_MESH_BSPLINE
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_BSPLINE_1D(self)
  IMPLICIT NONE
     CLASS(DEF_MESH_BSPLINE_1D), INTENT(INOUT) :: self
     ! LOCAL

     self % c_dim = 1

     ALLOCATE(self % control_points(self % c_dim))
     ALLOCATE(self % control_points(1) % coef(self % logicals(1) % ptr_mesh % n, self % d_dim))

  END SUBROUTINE CREATE_MESH_BSPLINE_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_BSPLINE_2D(self)
  IMPLICIT NONE
     CLASS(DEF_MESH_BSPLINE_2D), INTENT(INOUT) :: self
     ! LOCAL

     self % c_dim = 1

     ALLOCATE(self % control_points(self % c_dim))
     ALLOCATE(self % control_points(1) % coef( &
                             &   self % logicals(1) % ptr_mesh % n &
                             & , self % logicals(2) % ptr_mesh % n &
                             & , self % d_dim)) 

  END SUBROUTINE CREATE_MESH_BSPLINE_2D
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_BSPLINE_3D(self)
  IMPLICIT NONE
     CLASS(DEF_MESH_BSPLINE_3D), INTENT(INOUT) :: self
     ! LOCAL

     self % c_dim = 1

     ALLOCATE(self % control_points(self % c_dim))
     ALLOCATE(self % control_points(1) % coef( &
                             &   self % logicals(1) % ptr_mesh % n &
                             & , self % logicals(2) % ptr_mesh % n &
                             & , self % logicals(3) % ptr_mesh % n &
                             & , self % d_dim)) 

  END SUBROUTINE CREATE_MESH_BSPLINE_3D
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_MESH_BSPLINE_1D1D(self)
  IMPLICIT NONE
     CLASS(DEF_MESH_BSPLINE_1D1D), INTENT(INOUT) :: self
     ! LOCAL

     self % c_dim = 2 

     ALLOCATE(self % control_points(self % c_dim))
     ALLOCATE(self % control_points(1) % coef(self % logicals(1) % ptr_mesh % n, self % d_dim))
     ALLOCATE(self % control_points(2) % coef(self % logicals(2) % ptr_mesh % n, self % d_dim))

  END SUBROUTINE CREATE_MESH_BSPLINE_1D1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE SET_MESH_BSPLINE_CONTROL_POINTS_1D(self, control_points, i_dimension)
  IMPLICIT NONE
     TYPE(DEF_MESH_BSPLINE_1D)     , INTENT(INOUT) :: self
      REAL(SPI_RK), DIMENSION (:,:), INTENT(IN)    :: control_points
      INTEGER            , OPTIONAL, INTENT(IN)    :: i_dimension 
     ! LOCAL
     INTEGER :: i

     i = 1
     IF (PRESENT(i_dimension)) THEN
        i = i_dimension
     END IF

     ! ... control points
     self % control_points (i) % coef(:,:) = control_points(:,:) 
  END SUBROUTINE SET_MESH_BSPLINE_CONTROL_POINTS_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE SET_MESH_BSPLINE_CONTROL_POINTS_2D(self, control_points, i_dimension)
  IMPLICIT NONE
     TYPE(DEF_MESH_BSPLINE_2D)       , INTENT(INOUT) :: self
      REAL(SPI_RK), DIMENSION (:,:,:), INTENT(IN)    :: control_points
      INTEGER              , OPTIONAL, INTENT(IN)    :: i_dimension 
     ! LOCAL
     INTEGER :: i

     i = 1
     IF (PRESENT(i_dimension)) THEN
        i = i_dimension
     END IF

     ! ... control points
     self % control_points (i) % coef(:,:,:) = control_points(:,:,:) 
  END SUBROUTINE SET_MESH_BSPLINE_CONTROL_POINTS_2D
  ! .........................................................

  ! .........................................................
  SUBROUTINE SET_MESH_BSPLINE_CONTROL_POINTS_3D(self, control_points, i_dimension)
  IMPLICIT NONE
     TYPE(DEF_MESH_BSPLINE_3D)         , INTENT(INOUT) :: self
      REAL(SPI_RK), DIMENSION (:,:,:,:), INTENT(IN)    :: control_points
      INTEGER                , OPTIONAL, INTENT(IN)    :: i_dimension 
     ! LOCAL
     INTEGER :: i

     i = 1
     IF (PRESENT(i_dimension)) THEN
        i = i_dimension
     END IF

     ! ... control points
     self % control_points (i) % coef(:,:,:,:) = control_points(:,:,:,:) 
  END SUBROUTINE SET_MESH_BSPLINE_CONTROL_POINTS_3D
  ! .........................................................

  ! .........................................................
  SUBROUTINE SET_MESH_BSPLINE_CONTROL_POINTS_1D1D(self, control_points, i_dimension)
  IMPLICIT NONE
     TYPE(DEF_MESH_BSPLINE_1D1D)   , INTENT(INOUT) :: self
      REAL(SPI_RK), DIMENSION (:,:), INTENT(IN)    :: control_points
      INTEGER                      , INTENT(IN)    :: i_dimension 
     ! LOCAL

     ! ... control points
     self % control_points (i_dimension) % coef(:,:) = control_points(:,:) 
  END SUBROUTINE SET_MESH_BSPLINE_CONTROL_POINTS_1D1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE INITIALIZE_MESH_CONTROL_POINTS_DEFAULT_1D(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_CONTROL_POINTS_1D), INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: i
     INTEGER :: n_u
     INTEGER :: n_total

     n_u = SIZE(self % coef, 1)

     n_total = (n_u - 1) 

     ! ... control points
     self % coef = 0.0
     do i = 1, n_u 
        self % coef ( i, : ) = (i - 1) * 1.0d0 / n_total
     end do 
     ! ...
  END SUBROUTINE INITIALIZE_MESH_CONTROL_POINTS_DEFAULT_1D
  ! .........................................................

  ! .........................................................
  SUBROUTINE INITIALIZE_MESH_CONTROL_POINTS_DEFAULT_2D(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_CONTROL_POINTS_2D), INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: i
     INTEGER :: j 
     INTEGER :: n_u
     INTEGER :: n_v
     INTEGER :: n_total

     n_u = SIZE(self % coef, 1)
     n_v = SIZE(self % coef, 2)

     n_total = (n_u - 1) * (n_v - 1) 

     ! ... control points
     self % coef = 0.0
     do i = 1, n_u 
        do j = 1, n_v 
           self % coef ( i, j, : ) = (i - 1) * (j - 1) * 1.0d0 / n_total
        end do
     end do 
     ! ...

  END SUBROUTINE INITIALIZE_MESH_CONTROL_POINTS_DEFAULT_2D
  ! .........................................................

  ! .........................................................
  SUBROUTINE INITIALIZE_MESH_CONTROL_POINTS_DEFAULT_3D(self)
  IMPLICIT NONE
     TYPE(DEF_MESH_CONTROL_POINTS_3D), INTENT(INOUT) :: self
     ! LOCAL
     INTEGER :: i
     INTEGER :: j 
     INTEGER :: k 
     INTEGER :: n_u
     INTEGER :: n_v
     INTEGER :: n_w
     INTEGER :: n_total

     n_u = SIZE(self % coef, 1)
     n_v = SIZE(self % coef, 2)
     n_w = SIZE(self % coef, 3)

     n_total = (n_u - 1) * (n_v - 1) * (n_w - 1)

     ! ... control points
     self % coef = 0.0
     do i = 1, n_u 
        do j = 1, n_v 
           do k = 1, n_w 
              self % coef ( i, j, k, : ) = (i - 1) * (j - 1) * (k - 1) * 1.0d0 / n_total
           end do
        end do
     end do 
     ! ...

  END SUBROUTINE INITIALIZE_MESH_CONTROL_POINTS_DEFAULT_3D
  ! .........................................................
!
!  ! .........................................................
!  SUBROUTINE CREATE_MESH_BSPLINES_1D(self, ai_n, ai_p, type_bc, knots, control_points)
!  IMPLICIT NONE
!     TYPE(DEF_MESH_BSPLINE_1D)      , INTENT(INOUT) :: self
!     INTEGER              , INTENT(IN)    :: ai_n
!     INTEGER              , INTENT(IN)    :: ai_p
!     INTEGER              , OPTIONAL, INTENT(IN)    :: type_bc
!     real(SPI_RK), dimension (:), OPTIONAL, INTENT(IN) :: knots
!     real(SPI_RK), dimension (:,:), OPTIONAL, INTENT(IN) :: control_points
!     ! LOCAL
!     INTEGER :: li_err 
!
!     IF ( ( PRESENT(control_points)) ) THEN
!        CALL INITIALIZE_MESH_1D_BSPLINE_CONTROL_POINTS(self, control_points) 
!     ELSE
!        CALL INITIALIZE_MESH_1D_BSPLINE_CONTROL_POINTS_DEFAULT(self) 
!     END IF
!
!  END SUBROUTINE CREATE_MESH_BSPLINES_1D
!  ! .........................................................
  
  ! .........................................................
END MODULE SPI_MESH_BSPLINE 
