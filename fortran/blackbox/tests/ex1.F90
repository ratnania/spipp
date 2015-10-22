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
USE SPI_QUADRATURES_DEF
USE SPI_QUADRATURES
USE SPI_MESH_DEF
USE SPI_MESH
USE SPI_BASIS_DEF
USE SPI_BASIS
USE SPI_BLACKBOX_DEF
USE SPI_BLACKBOX
implicit none
   ! LOCAL
   TYPE(DEF_QUADRATURE_1D), TARGET :: lo_quad
   TYPE(DEF_MESH_1D_BSPLINE), TARGET :: lo_mesh
   TYPE(DEF_BASIS_1D_BSPLINE), TARGET :: lo_basis
   TYPE(DEF_BLACKBOX_1D_BSPLINE), TARGET :: lo_bbox
   ! ... number of internal knots is = N - P - 1
!   INTEGER, PARAMETER :: P = 1 
   INTEGER, PARAMETER :: P =  3 
!   INTEGER, PARAMETER :: M =  3 
   INTEGER, PARAMETER :: M = 127 
   INTEGER, PARAMETER :: N = P + 1 + M 
   REAL(SPI_RK), DIMENSION(N+P+1) :: KNOTS
   INTEGER :: K = P 
   INTEGER :: li_elmt
   INTEGER :: i, b, e
   REAL(SPI_RK) :: lr_a
   REAL(SPI_RK) :: lr_b
   real :: start, finish
   INTEGER, PARAMETER :: N_DIM = 1 
   REAL(SPI_RK), DIMENSION(N, N_DIM) :: control_points_1d 
   INTEGER, PARAMETER :: n_iterations = 10 
   INTEGER :: i_iteration

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

   CALL CREATE_QUADRATURE(lo_quad, SPI_QUADRATURES_LEGENDRE, K)
   CALL CREATE_MESH(lo_mesh, lo_quad, N, P, apr_knots=KNOTS) 
   CALL CREATE_BASIS(lo_basis, lo_mesh, lo_quad) 
   CALL CREATE_BLACKBOX(lo_bbox, lo_mesh, lo_basis, lo_quad, control_points_1d=control_points_1d)

   ! ...
   call cpu_time(start)
   DO i_iteration=1, n_iterations
     CALL UPDATE_POSITION_BLACKBOX(lo_bbox)
     CALL COMPUTE_METRIC_BLACKBOX(lo_bbox) 
   END DO
   call cpu_time(finish)
   print '("Time to evaluate positions and jacobians on quadrature = ",f6.3," seconds.")',(finish-start)/n_iterations
   ! ...

   CALL FREE_QUADRATURE(lo_quad) 
   CALL FREE_MESH(lo_mesh) 
   CALL FREE_BASIS(lo_basis) 
   CALL FREE_BLACKBOX(lo_bbox) 

end subroutine test1
! ............................................
   
! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
