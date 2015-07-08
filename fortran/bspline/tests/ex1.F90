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
SUBROUTINE SPAN_TEST_1()
USE bsp
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 4
     INTEGER, PARAMETER :: p = 2
     REAL*8, DIMENSION(n+p+1) :: U
     INTEGER, PARAMETER :: nx = 5
     REAL*8, DIMENSION(nx) :: X
     INTEGER, DIMENSION(nx) :: expected
     INTEGER :: i
     INTEGER :: span 

     PRINT *, ">>>> SPAN_TEST_1: Begin "

     U        = (/0.0,  0.0, 0.0,  0.5, 1.0, 1.0, 1.0/)
     X        = (/0.0, 0.25, 0.5, 0.75, 1.0/)
     expected = (/  2,    2,   3,    3, 4/)

     DO i = 1, nx
        CALL FindSpan(p,n,U,X(i),span)
        PRINT *, span, expected(i)
     END DO

     PRINT *, ">>>> SPAN_TEST_1: End "
END SUBROUTINE SPAN_TEST_1
! ............................................

! ............................................
SUBROUTINE SPAN_TEST_2()
USE bsp
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 4
     INTEGER, PARAMETER :: p = 2
     REAL*8, DIMENSION(n+p+1) :: U
     INTEGER, PARAMETER :: nx = 5
     REAL*8, DIMENSION(nx) :: X
     INTEGER, DIMENSION(nx) :: expected
     INTEGER :: i
     INTEGER :: span 

     PRINT *, ">>>> SPAN_TEST_2: Begin "

     U        = (/-1.0,0.0,0.0, 0.5, 1.0,1.0, 2.0/)
     X        = (/0.0, 0.25, 0.5, 0.75, 1.0/)
     expected = (/  2,    2,   3,    3, 4/)

     DO i = 1, nx
        CALL FindSpan(p,n,U,X(i),span)
        PRINT *, span, expected(i)
     END DO

     PRINT *, ">>>> SPAN_TEST_2: End "
END SUBROUTINE SPAN_TEST_2
! ............................................

! ............................................
SUBROUTINE MULTIPLICITY_TEST_1()
USE bsp
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 4
     INTEGER, PARAMETER :: p = 2
     REAL*8, DIMENSION(n+p+1) :: U
     INTEGER, DIMENSION(n+p+1) :: expected
     INTEGER :: i
     INTEGER :: span 
     INTEGER :: mult 

     PRINT *, ">>>> MULTIPLICITY_TEST_1: Begin "

     U        = (/0.0,  0.0, 0.0,  0.5, 1.0, 1.0, 1.0/)
     expected = (/  3,    3,   3,    1,   3,   3,   3/)

     DO i = 1, n+p+1
        CALL FindSpanMult(n,p,U,U(i),span,mult)
        PRINT *, mult, expected(i)
     END DO

     PRINT *, ">>>> MULTIPLICITY_TEST_1: End "
END SUBROUTINE MULTIPLICITY_TEST_1
! ............................................



! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL SPAN_TEST_1()
  CALL SPAN_TEST_2()
  CALL MULTIPLICITY_TEST_1()

END PROGRAM Main
! ............................................
