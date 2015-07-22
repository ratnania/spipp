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
SUBROUTINE EvalBasisFunsDers_TEST_1()
USE bsp
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 7 
     INTEGER, PARAMETER :: p = 3 
     REAL*8, DIMENSION(n+p+1) :: U
     INTEGER, PARAMETER :: nx = 5000
     INTEGER, PARAMETER :: nderiv = 2
     REAL*8, DIMENSION(nx) :: X
     REAL*8, DIMENSION(p+1,nderiv+1) :: dbatx
     INTEGER :: i
     INTEGER :: span 

     PRINT *, ">>>> EvalBasisFunsDers_TEST_1: Begin "

     U        = (/0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0/)

     DO i = 1, nx
        X(i) = FLOAT(i) / FLOAT(nx) 
     END DO

     DO i = 1, nx
        CALL FindSpan(p,n,U,X(i),span)
        CALL EvalBasisFunsDers(p,n,U,X(i),nderiv,span,dbatx)
!        PRINT *, ">>>>"
!        PRINT *, X(i)
!        PRINT *, dbatx
     END DO

     PRINT *, ">>>> EvalBasisFunsDers_TEST_1: End "

END SUBROUTINE EvalBasisFunsDers_TEST_1
! ............................................

! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL EvalBasisFunsDers_TEST_1()

END PROGRAM Main
! ............................................
