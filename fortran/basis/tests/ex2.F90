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
USE SPI_BASIS_DEF
USE SPI_BASIS
implicit none
   ! LOCAL
   TYPE(DEF_MESH_1D_BSPLINE), TARGET :: lo_mesh
   TYPE(DEF_QUADRATURE_1D), TARGET :: lo_quad
   TYPE(DEF_BASIS_1D_BSPLINE), TARGET :: lo_basis
   ! ... number of internal knots is = N - P - 1
   INTEGER, PARAMETER :: P = 1 
!   INTEGER, PARAMETER :: N = P + 4 
   INTEGER, PARAMETER :: N = P + 1 + 127 
   REAL(SPI_RK), DIMENSION(N+P+1) :: KNOTS
   INTEGER, PARAMETER :: n_iterations = 10 
   INTEGER :: i_iteration
   INTEGER :: li_elmt
   INTEGER :: i, b, e
   REAL(SPI_RK), DIMENSION(:,:), ALLOCATABLE :: COEFFS_LOC
   REAL(SPI_RK), DIMENSION(:,:), ALLOCATABLE :: s_at_quad
   real :: start, finish

   KNOTS(1:P+1) = 0.0
   DO i =1, N - P - 1 
      KNOTS(i+P+1) = i * 1.0 / (N-P) 
   END DO
   KNOTS(N+1:N+P+1) = 1.0

   CALL CREATE_QUADRATURE(lo_quad, SPI_QUADRATURES_LEGENDRE, P) 
   CALL CREATE_MESH(lo_mesh, lo_quad, N, P, apr_knots=KNOTS) 
   CALL CREATE_BASIS(lo_basis, lo_mesh, lo_quad) 

!   PRINT *, ">>> points"
!   DO e = 1, 4
!      PRINT *, lo_mesh % opr_points(e,:)
!   END DO

   ALLOCATE(COEFFS_LOC(lo_mesh % n_elements, lo_mesh % oi_p + 1))
   ALLOCATE(s_at_quad(lo_mesh % n_elements, lo_quad % oi_n_points))
  
   ! ... TODO set coeff values
   call cpu_time(start)
   DO i_iteration=1, n_iterations
      COEFFS_LOC = 1.0 
      s_at_quad = 0.0
      DO e = 1, lo_mesh % n_elements 
         DO b = 1, lo_mesh % oi_p + 1 
            s_at_quad(e, :) = s_at_quad(e, :) + lo_basis % B_0(e, :, b) * COEFFS_LOC(e, b)
         END DO
      END DO
   END DO
   call cpu_time(finish)
   print '("Time = ",f6.3," seconds.")',(finish-start)/n_iterations

   CALL FREE_BASIS(lo_basis) 
   CALL FREE_MESH(lo_mesh) 
   CALL FREE_QUADRATURE(lo_quad) 
end subroutine test1
! ............................................
   
! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
