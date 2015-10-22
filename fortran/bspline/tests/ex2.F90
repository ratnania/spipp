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
SUBROUTINE EvalBasisFuns_TEST_1()
use bspline, FindS => FindSpan
USE bsp
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 6 
     INTEGER, PARAMETER :: p = 2
     REAL*8, DIMENSION(n+p+1) :: U
     INTEGER, PARAMETER :: nx = 5
     REAL*8, DIMENSION(nx) :: X
     REAL*8, DIMENSION(nx, p+1) :: expected
     REAL*8, DIMENSION(p+1) :: batx
     INTEGER :: i
     INTEGER :: span 

     PRINT *, ">>>> EvalBasisFuns_TEST_1: Begin "

     U        = (/0.0,  0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0/)

     ! ...
     DO i = 1, nx
        X(i) = FLOAT(i-1) / FLOAT(nx-1) 
     END DO
     ! ...

     ! ...
     expected(1,:) = (/ 1.0000000000000000, 0.0000000000000000, 0.0000000000000000/)
     expected(2,:) = (/ 0.5000000000000000, 0.5000000000000000, 0.0000000000000000/)
     expected(3,:) = (/ 0.5000000000000000, 0.5000000000000000, 0.0000000000000000/)
     expected(4,:) = (/ 0.5000000000000000, 0.5000000000000000, 0.0000000000000000/)
     expected(5,:) = (/ 0.0000000000000000, 0.0000000000000000, 1.0000000000000000/)
     ! ...

     ! ...
     OPEN(UNIT=12, FILE="bspline_ex2_test_1.txt"&
             & , ACTION="write", STATUS="replace")
     DO i = 1, nx
        span = FindS(n-1,p,X(i),U)
        CALL EvalBasisFuns(p,n,U,X(i),span,Batx)
        print *, MAXVAL(ABS(Batx - expected(i,:)))

        WRITE(12,*) Batx 
     END DO
     CLOSE(12)
     ! ...

     PRINT *, ">>>> EvalBasisFuns_TEST_1: End "

END SUBROUTINE EvalBasisFuns_TEST_1
! ............................................

! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL EvalBasisFuns_TEST_1()

END PROGRAM Main
! ............................................
