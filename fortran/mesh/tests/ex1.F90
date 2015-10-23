!
! File: ex1.f90
!
!
! Usage:
!   > ./ex1 
!
! Authors:
!   Ahmed RATNANI  - ratnaniahmed@gmail.com
!

! ............................................
subroutine test1 ()
USE SPI_GLOBAL_DEF
USE SPI_GLOBAL
USE SPI_MESH_LOGICAL_DEF
USE SPI_MESH_LOGICAL
USE SPI_MESH_DEF
USE SPI_MESH_BSPLINE
USE SPI_QUADRATURES_DEF
USE SPI_QUADRATURES
implicit none
   ! LOCAL
   TYPE(DEF_MESH_BSPLINE_1D), TARGET :: mesh
   ! ... number of internal knots is = N - P - 1
!   INTEGER, PARAMETER :: P = 1 
   INTEGER, PARAMETER :: P =  3 
   INTEGER, PARAMETER :: M =  3 
!   INTEGER, PARAMETER :: M = 127 
   INTEGER, PARAMETER :: N = P + 1 + M 
   REAL(SPI_RK), DIMENSION(N+P+1) :: KNOTS   
   TYPE(DEF_QUADRATURE_1D), TARGET :: quad
   INTEGER, PARAMETER :: K = 3 
   INTEGER, PARAMETER :: D_DIM = 1 
   INTEGER :: e
   INTEGER :: i 
   TYPE(DEF_MESH_LOGICAL_BSPLINE), TARGET :: mesh_logical
   INTEGER, PARAMETER :: N_DIM = 1 
   REAL(SPI_RK), DIMENSION(N, N_DIM) :: control_points_1d 

   KNOTS(1:P+1) = 0.0
   DO i =1, N - P - 1 
      KNOTS(i+P+1) = i * 1.0 / (N-P) 
   END DO
   KNOTS(N+1:N+P+1) = 1.0

   ! ... control points
   control_points_1d = 0.0
   do i = 1, N 
      control_points_1d( i, 1 ) = (i - 1) * 1.0d0 / (N - 1)
   end do 
   ! ...

   CALL CREATE_QUADRATURE(quad, SPI_QUADRATURES_LEGENDRE, K)
   CALL CREATE_MESH_LOGICAL(mesh_logical, quad=quad, n=N, p=P, knots=KNOTS) 
   CALL CREATE_MESH_BSPLINE(   mesh                         &
                           & , d_dim=D_DIM                  &
                           & , logical_mesh_s1=mesh_logical &
                           &  )

   CALL SET_MESH_BSPLINE_CONTROL_POINTS_1D(mesh, control_points_1d)

   PRINT *, ">>> N, P :", mesh_logical % n, mesh_logical % p

   PRINT *, ">>> knots"
   PRINT *, mesh_logical % knots

   PRINT *, ">>> grid"
   PRINT *, mesh_logical % grid

   PRINT *, ">>> points"
   DO e = 1, 4
      PRINT *, mesh_logical % points(e,:)
   END DO

   CALL FREE_MESH_LOGICAL(mesh_logical) 
   CALL FREE_MESH_BSPLINE(mesh) 

end subroutine test1
! ............................................
   
! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
