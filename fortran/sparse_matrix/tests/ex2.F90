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
USE bsp
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 6 
     INTEGER, PARAMETER :: p = 2
     REAL*8, DIMENSION(n+p+1) :: U
     INTEGER, PARAMETER :: nx = 5
     REAL*8, DIMENSION(nx) :: X
     REAL*8, DIMENSION(p+1) :: batx
     INTEGER :: i
     INTEGER :: span 

     PRINT *, ">>>> EvalBasisFuns_TEST_1: Begin "

     U        = (/0.0,  0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0/)

     DO i = 1, nx
        X(i) = FLOAT(i) / FLOAT(nx) 
     END DO

     DO i = 1, nx
        CALL FindSpan(p,n,U,X(i),span)
        CALL EvalBasisFuns(p,n,U,X(i),span,Batx)
        PRINT *, ">>>>"
        PRINT *, X(i)
        PRINT *, Batx
     END DO

     PRINT *, ">>>> EvalBasisFuns_TEST_1: End "

END SUBROUTINE EvalBasisFuns_TEST_1
! ............................................

! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL EvalBasisFuns_TEST_1()

END PROGRAM Main
! ............................................
