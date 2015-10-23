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
USE SPI_MESH_DEF
USE SPI_MESH
USE SPI_QUADRATURES_DEF
USE SPI_QUADRATURES
implicit none
   ! LOCAL
   TYPE(DEF_MESH_2D_BSPLINE), TARGET :: mesh
   ! ... number of internal knots is = N - P - 1
   INTEGER, PARAMETER :: P_u =  3 
   INTEGER, PARAMETER :: M_u = 7 
   INTEGER, PARAMETER :: N_u = P_u + 1 + M_u 
   REAL(SPI_RK), DIMENSION(N_u+P_u+1) :: KNOTS_u
   INTEGER, PARAMETER :: P_v =  2 
   INTEGER, PARAMETER :: M_v =  3 
   INTEGER, PARAMETER :: N_v = P_v + 1 + M_v 
   REAL(SPI_RK), DIMENSION(N_v+P_v+1) :: KNOTS_v
   TYPE(DEF_QUADRATURE_1D), TARGET :: quad_u
   TYPE(DEF_QUADRATURE_1D), TARGET :: quad_v
   INTEGER, PARAMETER :: K_u = 3 
   INTEGER, PARAMETER :: K_v = 2 
   INTEGER :: i 
   INTEGER :: j 
   INTEGER, PARAMETER :: N_DIM = 1 
   REAL(SPI_RK), DIMENSION(N_u, N_v, N_DIM) :: control_points_2d 

   KNOTS_u(1:P_u+1) = 0.0
   DO i =1, N_u - P_u - 1 
      KNOTS_u(i+P_u+1) = i * 1.0 / (N_u-P_u) 
   END DO

   KNOTS_v(1:P_v+1) = 0.0
   DO i =1, N_v - P_v - 1 
      KNOTS_v(i+P_v+1) = i * 1.0 / (N_v-P_v) 
   END DO

   ! ... control points
   control_points_2d = 0.0
   do i = 1, N_u 
      do j = 1, N_v 
         control_points_2d( i, j, : ) = (i - 1) * (j - 1) * 1.0d0 / ( (N_u - 1) * (N_v - 1)) 
      end do 
   end do 
   ! ...

   CALL CREATE_QUADRATURE(quad_u, SPI_QUADRATURES_LEGENDRE, K_u)
   CALL CREATE_QUADRATURE(quad_v, SPI_QUADRATURES_LEGENDRE, K_v)
   CALL CREATE_MESH(mesh &
           & , quad_u=quad_u, n_u=N_u, p_u=P_u, knots_u=KNOTS_u &
           & , quad_v=quad_v, n_v=N_v, p_v=P_v, knots_v=KNOTS_v) 

   PRINT *, ">>> N_u, P_u :", mesh % mesh_u % n, mesh % mesh_u % p
   PRINT *, ">>> N_v, P_v :", mesh % mesh_v % n, mesh % mesh_v % p

   CALL FREE_MESH(mesh) 
end subroutine test1
! ............................................
   
! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
