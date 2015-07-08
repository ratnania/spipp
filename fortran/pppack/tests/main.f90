program main
    use tests
    use tracelog_module
    implicit none

    logical, parameter :: ll_stdoutput = .TRUE.
    integer, parameter :: li_detail = 3

    call opentracelog(al_stdoutput = ll_stdoutput, ai_dtllevel = li_detail)

!    call test11()
    call test_bsplvd()
    call test_bsplvd_vect()

!    call test1()
!    call test_bsplvb()
!    call test_bsplvb_vect()

end
