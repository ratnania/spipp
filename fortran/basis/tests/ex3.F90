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
   TYPE(DEF_MESH_1D_BSPLINE) , TARGET :: lo_mesh_u
   TYPE(DEF_QUADRATURE_1D)   , TARGET :: lo_quad_u
   TYPE(DEF_BASIS_1D_BSPLINE), TARGET :: lo_basis_u
   TYPE(DEF_MESH_1D_BSPLINE) , TARGET :: lo_mesh_v
   TYPE(DEF_QUADRATURE_1D)   , TARGET :: lo_quad_v
   TYPE(DEF_BASIS_1D_BSPLINE), TARGET :: lo_basis_v
   ! ... number of internal knots is = N - P - 1
!   INTEGER, PARAMETER :: P =  1 
!   INTEGER, PARAMETER :: P =  2 
   INTEGER, PARAMETER :: P =  3 
   INTEGER, PARAMETER :: N = P + 1 + 127 
!   INTEGER, PARAMETER :: N = P + 1 + 1027 
   REAL(SPI_RK), DIMENSION(N+P+1) :: KNOTS
   INTEGER, PARAMETER :: n_iterations = 10 
   INTEGER :: i_iteration
   INTEGER :: li_elmt
   INTEGER :: i, i_u, i_v, b_u, b_v, e, e_u, e_v
   REAL(SPI_RK), DIMENSION(:,:,:,:), ALLOCATABLE :: COEFFS_LOC
   REAL(SPI_RK), DIMENSION(:,:), ALLOCATABLE :: s_at_quad
   REAL(SPI_RK) :: W_at_u 
   REAL(SPI_RK) :: W_at_v 
   REAL(SPI_RK) :: B_at_u 
   REAL(SPI_RK) :: B_at_v 
   REAL(SPI_RK) :: integral  
   real :: start, finish

   KNOTS(1:P+1) = 0.0
   DO i =1, N - P - 1 
      KNOTS(i+P+1) = i * 1.0 / (N-P) 
   END DO
   KNOTS(N+1:N+P+1) = 1.0

   ! ... u
   CALL CREATE_QUADRATURE(lo_quad_u, SPI_QUADRATURES_LEGENDRE, P) 
   CALL CREATE_MESH(lo_mesh_u, quad_u=lo_quad_u, n_u=N, p_u=P, knots_u=KNOTS) 
   CALL CREATE_BASIS(lo_basis_u, lo_mesh_u, lo_quad_u) 
   ! ...

   ! ... v
   CALL CREATE_QUADRATURE(lo_quad_v, SPI_QUADRATURES_LEGENDRE, P) 
   CALL CREATE_MESH(lo_mesh_v, quad_u=lo_quad_u, n_u=N, p_u=P, knots_u=KNOTS) 
   CALL CREATE_BASIS(lo_basis_v, lo_mesh_v, lo_quad_v) 
   ! ...

   ALLOCATE(COEFFS_LOC(lo_mesh_u % n_elements, lo_mesh_v % n_elements, lo_mesh_u % p + 1, lo_mesh_v % p + 1))
   ALLOCATE(s_at_quad(lo_mesh_u % n_elements * lo_mesh_v % n_elements, lo_quad_u % oi_n_points * lo_quad_v % oi_n_points))
  
!   ! ...
!   call cpu_time(start)
!   DO i_iteration=1, n_iterations
!      COEFFS_LOC = 1.0 
!      s_at_quad = 0.0
!      e = 1
!      DO e_u = 1, lo_mesh_u % n_elements 
!         DO e_v = 1, lo_mesh_v % n_elements 
!            DO b_u = 1, lo_mesh_u % p + 1 
!               DO b_v = 1, lo_mesh_v % p + 1 
!                  i = 1
!                  DO i_u = 1, lo_quad_u %  oi_n_points
!                     B_at_u = lo_basis_u % B_0(e_u, i_u, b_u)
!                     DO i_v = 1, lo_quad_u %  oi_n_points
!                        B_at_v = lo_basis_v % B_0(e_v, i_v, b_v)
!                        s_at_quad(e, i) = s_at_quad(e, i) + B_at_u * B_at_v * COEFFS_LOC(e_u, e_v, b_u, b_v)
!                        i = i + 1
!                     END DO
!                  END DO
!               END DO
!            END DO
!            e = e + 1
!         END DO
!      END DO
!   END DO
!   call cpu_time(finish)
!   print '("Time to evaluate surface on quadrature = ",f6.3," seconds.")',(finish-start)/n_iterations
!   ! ...
    
   ! ...
   call cpu_time(start)
   DO i_iteration=1, n_iterations
      COEFFS_LOC = 1.0 
      s_at_quad = 0.0
      DO b_u = 1, lo_mesh_u % p + 1 
         DO b_v = 1, lo_mesh_v % p + 1 
            i = 1
            DO i_u = 1, lo_quad_u %  oi_n_points
               DO i_v = 1, lo_quad_u %  oi_n_points
                  e = 1
                  DO e_u = 1, lo_mesh_u % n_elements 
                     B_at_u = lo_basis_u % B_0(e_u, i_u, b_u)
                     DO e_v = 1, lo_mesh_v % n_elements 
                        B_at_v = lo_basis_v % B_0(e_v, i_v, b_v)
                        s_at_quad(e, i) = s_at_quad(e, i) + B_at_u * B_at_v * COEFFS_LOC(e_u, e_v, b_u, b_v)
                        e = e + 1
                     END DO
                  END DO
                  i = i + 1
               END DO
            END DO
         END DO
      END DO
   END DO
   call cpu_time(finish)
   print '("Time to evaluate surface on quadrature = ",f6.3," seconds.")',(finish-start)/n_iterations
   ! ...

   ! ...
   call cpu_time(start)
   DO i_iteration=1, n_iterations
      COEFFS_LOC = 1.0 
      integral = 0.0
      DO b_u = 1, lo_mesh_u % p + 1 
         DO b_v = 1, lo_mesh_v % p + 1 
            i = 1
            DO i_u = 1, lo_quad_u %  oi_n_points
               W_at_u = lo_quad_u % opr_weights(i_u)
               DO i_v = 1, lo_quad_u %  oi_n_points
                  W_at_v = lo_quad_v % opr_weights(i_v)

                  e = 1
                  DO e_u = 1, lo_mesh_u % n_elements 
                     B_at_u = lo_basis_u % B_0(e_u, i_u, b_u)
                     DO e_v = 1, lo_mesh_v % n_elements 
                        B_at_v = lo_basis_v % B_0(e_v, i_v, b_v)
                        integral = integral + W_at_u * W_at_v * B_at_u * B_at_v * COEFFS_LOC(e_u, e_v, b_u, b_v)
                        e = e + 1
                     END DO
                  END DO
                  i = i + 1
               END DO
            END DO
         END DO
      END DO
   END DO
   call cpu_time(finish)
   print '("Time to compute an integral on quadrature = ",f6.3," seconds.")',(finish-start)/n_iterations
   ! ...

   CALL FREE_BASIS(lo_basis_u) 
   CALL FREE_MESH(lo_mesh_u) 
   CALL FREE_QUADRATURE(lo_quad_u) 
   CALL FREE_BASIS(lo_basis_v) 
   CALL FREE_MESH(lo_mesh_v) 
   CALL FREE_QUADRATURE(lo_quad_v) 
end subroutine test1
! ............................................
   
! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL test1()

END PROGRAM Main
! ............................................
   
