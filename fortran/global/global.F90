!# -*- coding: utf8 -*-
MODULE SPI_GLOBAL
  USE SPI_GLOBAL_DEF
IMPLICIT NONE

  CONTAINS

  ! .............................................
  SUBROUTINE SPI_GET_ARGUMENTS(argname, argvalue, IERROR)
  IMPLICIT NONE
     CHARACTER(LEN = 1024), INTENT(IN)           :: argname
     character(len=1024), INTENT(INOUT) ::argvalue
     INTEGER,      INTENT(OUT) :: IERROR

     ! LOCAL
     integer::narg,cptArg
     character(len=1024)::locname
     logical::ll_lookForValue=.FALSE.
     logical::fileExist
 
     !Check if arguments are found
     narg = command_argument_count()
     IF (narg .EQ. 0) THEN
        RETURN
     END IF

     !loop across options
     ll_lookForValue=.FALSE. 
     cptArg=1
     DO WHILE ( (.NOT. ll_lookForValue) .AND. (cptArg < narg) ) 
        call get_command_argument(cptArg,locname)
        IF (adjustl(locname) .EQ. TRIM(argname)) THEN
           ll_lookForValue=.TRUE. 
        END IF
        cptArg = cptArg + 1
     END DO

     IF (ll_lookForValue) THEN
        call get_command_argument(cptArg,locname)

        argvalue=adjustl(locname) !assign a value to pedfile
!        inquire(file=argvalue,exist=fileExist)!check if it exist
!        if(.not.fileExist)then
!           write(*,*)'argument value:',argvalue,' not found'
!           stop
!        else
!           write(*,*)'argument value:',TRIM(argvalue),'  found'
!        endif
      END IF

  END SUBROUTINE SPI_GET_ARGUMENTS
  ! .............................................

END MODULE SPI_GLOBAL
