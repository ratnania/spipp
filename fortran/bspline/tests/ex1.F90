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
USE bspline
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 4
     INTEGER, PARAMETER :: p = 2
     REAL(8), DIMENSION(0:n+p) :: U
     INTEGER, PARAMETER :: nx = 5
     REAL(8), DIMENSION(nx) :: X
     INTEGER, DIMENSION(nx) :: expected
     INTEGER :: i
     INTEGER :: span 

     PRINT *, ">>>> SPAN_TEST_1: Begin "

     U        = (/0.0,  0.0, 0.0,  0.5, 1.0, 1.0, 1.0/)
     X        = (/0.0, 0.25, 0.5, 0.75, 1.0/)
     expected = (/  2,    2,   3,    3,   3/)

     DO i = 1, nx
        span = FindSpan(n-1,p,X(i),U)
        PRINT *, "x = ", X(i), "span = ", span, " expected = ", expected(i)
     END DO

     PRINT *, ">>>> SPAN_TEST_1: End "
END SUBROUTINE SPAN_TEST_1
! ............................................

! ............................................
SUBROUTINE SPAN_TEST_2()
USE bspline
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 4
     INTEGER, PARAMETER :: p = 2
     REAL(8), DIMENSION(0:n+p) :: U
     INTEGER, PARAMETER :: nx = 5
     REAL(8), DIMENSION(nx) :: X
     INTEGER, DIMENSION(nx) :: expected
     INTEGER :: i
     INTEGER :: span 

     PRINT *, ">>>> SPAN_TEST_2: Begin "

     U        = (/-1.0,0.0,0.0, 0.5, 1.0,1.0, 2.0/)
     X        = (/0.0, 0.25, 0.5, 0.75, 1.0/)
     expected = (/  2,    2,   3,    3,  3/)

     DO i = 1, nx
        span = FindSpan(n-1,p,X(i),U)
        PRINT *, "x = ", X(i), "span = ", span, " expected = ", expected(i)
     END DO

     PRINT *, ">>>> SPAN_TEST_2: End "
END SUBROUTINE SPAN_TEST_2
! ............................................

! ............................................
SUBROUTINE SPAN_TEST_3()
USE bspline
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 3
     INTEGER, PARAMETER :: p = 2
     REAL(8), DIMENSION(0:n+p) :: U
     INTEGER, PARAMETER :: nx = 3 
     REAL(8), DIMENSION(nx) :: X
     INTEGER, DIMENSION(nx) :: expected
     INTEGER :: i
     INTEGER :: span 

     PRINT *, ">>>> SPAN_TEST_3: Begin "

     U        = (/-2.0, -1.0, 0.0, 1.0, 2.0, 3.0/)
     X        = (/0.0, 0.5, 1.0/)
     expected = (/  2,   2,   2/)

     DO i = 1, nx
        span = FindSpan(n-1,p,X(i),U)
        PRINT *, "x = ", X(i), "span = ", span, " expected = ", expected(i)
     END DO

     PRINT *, ">>>> SPAN_TEST_3: End "
END SUBROUTINE SPAN_TEST_3
! ............................................

! ............................................
SUBROUTINE MULTIPLICITY_TEST_1()
USE bspline
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 6
     INTEGER, PARAMETER :: p = 2
     REAL(8), DIMENSION(0:n+p) :: U
     INTEGER, PARAMETER :: nx = 5 
     REAL(8), DIMENSION(nx) :: X
     INTEGER, DIMENSION(nx) :: expected
     INTEGER :: i
     INTEGER :: span 
     INTEGER :: mult 

     PRINT *, ">>>> MULTIPLICITY_TEST_1: Begin "

     U        = (/0.0,  0.0, 0.0, 0.25, 0.25, 0.5, 1.0, 1.0, 1.0/)
     X        = (/0.0, 0.25, 0.5, 0.75, 1.0/)
     expected = (/  3,    2,   1,    0, 3/)

     DO i = 1, nx
        span = FindSpan(n-1,p,U(i),U)
        mult = FindMult(span,X(i),p,U)
        PRINT *, "u = ", X(i), "mult = ", mult, " expected = ", expected(i)
     END DO

     PRINT *, ">>>> MULTIPLICITY_TEST_1: End "
END SUBROUTINE MULTIPLICITY_TEST_1
! ............................................

! ............................................
SUBROUTINE MULTIPLICITY_TEST_2()
USE bspline
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 6
     INTEGER, PARAMETER :: p = 2
     REAL(8), DIMENSION(0:n+p) :: U
     INTEGER, PARAMETER :: nx = 5 
     REAL(8), DIMENSION(nx) :: X
     INTEGER, DIMENSION(nx) :: expected
     INTEGER :: i
     INTEGER :: span 
     INTEGER :: mult 

     PRINT *, ">>>> MULTIPLICITY_TEST_2: Begin "

     U        = (/-1.0,  0.0, 0.0, 0.25, 0.25, 0.5, 1.0, 1.0, 2.0/)
     X        = (/0.0, 0.25, 0.5, 0.75, 1.0/)
     expected = (/  2,    2,   1,    0, 2/)

     DO i = 1, nx
        span = FindSpan(n-1,p,U(i),U)
        mult = FindMult(span,X(i),p,U)
        PRINT *, "u = ", X(i), "mult = ", mult, " expected = ", expected(i)
     END DO

     PRINT *, ">>>> MULTIPLICITY_TEST_2: End "
END SUBROUTINE MULTIPLICITY_TEST_2
! ............................................

! ............................................
SUBROUTINE MULTIPLICITY_TEST_3()
USE bspline
IMPLICIT NONE
     INTEGER, PARAMETER :: n = 6
     INTEGER, PARAMETER :: p = 2
     REAL(8), DIMENSION(0:n+p) :: U
     INTEGER, PARAMETER :: nx = 5 
     REAL(8), DIMENSION(nx) :: X
     INTEGER, DIMENSION(nx) :: expected
     INTEGER :: i
     INTEGER :: span 
     INTEGER :: mult 

     PRINT *, ">>>> MULTIPLICITY_TEST_3: Begin "

     U        = (/-2.0,  -1.0, 0.0, 0.25, 0.25, 0.5, 1.0, 2.0, 3.0/)
     X        = (/0.0, 0.25, 0.5, 0.75, 1.0/)
     expected = (/  1,    2,   1,    0, 1/)

     DO i = 1, nx
        span = FindSpan(n-1,p,U(i),U)
        mult = FindMult(span,X(i),p,U)
        PRINT *, "u = ", X(i), "mult = ", mult, " expected = ", expected(i)
     END DO

     PRINT *, ">>>> MULTIPLICITY_TEST_3: End "
END SUBROUTINE MULTIPLICITY_TEST_3
! ............................................


! ............................................
PROGRAM Main

  IMPLICIT NONE

  CALL SPAN_TEST_1()
  CALL SPAN_TEST_2()
  CALL SPAN_TEST_3()
  CALL MULTIPLICITY_TEST_1()
  CALL MULTIPLICITY_TEST_2()
  CALL MULTIPLICITY_TEST_3()

END PROGRAM Main
! ............................................
