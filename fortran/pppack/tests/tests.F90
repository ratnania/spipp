!     
! File:   tests.F90
! Author: ratnani
!
! Created on November 13, 2011, 18:09 AM
!

module tests
    use pppack
    use tracelog_module

    implicit none

contains
    !---------------------------------------------------------------
    subroutine test1()
        implicit none
        ! LOCAL
        ! LOCAL VARIABLES
        !        integer, parameter :: k = 4
        !        integer, parameter :: N = 7
        
        integer, parameter :: k = 4
        integer, parameter :: N = k
        integer, parameter :: Niter = 200000
        integer, parameter :: dimx = 1000
        real(8), dimension(N + k) :: t
        real(8), dimension(k, dimx) :: B_vectiatx
        real(8), dimension(k, dimx) :: Biatx
        real(8), dimension(k) :: locBiatx
        real(8), dimension(dimx) :: x
        real(8) :: xmin
        real(8) :: xmax
        real(8) :: dx
        real(8) :: error
        integer :: left
        integer :: mflag
        integer :: iter
        integer :: i
        integer :: j
        REAL time_begin, time_end
        real ( kind = 8),  dimension ( k,dimx) :: deltal
        real ( kind = 8),  dimension ( k,dimx) :: deltar
        real ( kind = 8), dimension(dimx) :: saved
        real ( kind = 8), dimension(dimx) :: term

        CALL printlog("===== Begin test1 =====", ai_dtllevel = 0)
        !        call printcputime()

        !        t = (/ 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0 /)
        !        xmin = 0.25
        !        xmax = 0.5

        t(1:k) = 0.0
        t(k + 1:N + k) = 1.0
        xmin = 0.0
        xmax = 1.0

        dx = (xmax - xmin)/float(dimx + 1)
        do i = 1, dimx
            x(i) = xmin + float(i) * dx
        end do

        call interv(t, N + 1, x(1), left, mflag)

        CALL CPU_TIME(time_begin)

        do iter = 1, Niter
            do i = 1, dimx
                CALL bsplvb(t, k, 1, x(i), left, locBiatx(:))
                Biatx(:, i) = locBiatx
            end do
        end do

        CALL CPU_TIME(time_end)
        PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds using bsplvb'

        !        call printcputime()
        ! ****************

        ! ****************
        CALL CPU_TIME(time_begin)

        do iter = 1, Niter            
            CALL bsplvb_vect(t, k, 1, dimx, x, left, B_vectiatx, deltal, deltar, saved, term)
        end do
        
        CALL CPU_TIME(time_end)
        PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds using bsplvb_vect'
        !        call printcputime()
        ! ****************

        error = 0.0
        do i = 1, dimx
            do j = 1, k
                error = error + (B_vectiatx(j, i) - Biatx(j, i))**2
            end do
        end do
        error = SQRT(error) / (dimx * k)
        print *, "error=", error

        CALL printlog("===== End test1 =====", ai_dtllevel = 0)


    end subroutine test1
!    !---------------------------------------------------------------
!    subroutine test1()
!        implicit none
!        ! LOCAL
!        ! LOCAL VARIABLES
!        !        integer, parameter :: k = 4
!        !        integer, parameter :: N = 7
!        
!        integer, parameter :: k = 4
!        integer, parameter :: N = k
!        integer, parameter :: Niter = 200000
!        integer, parameter :: dimx = 6*k
!        real(8), dimension(N + k) :: t
!        real(8), dimension(k, dimx) :: B_vectiatx
!        real(8), dimension(k, dimx) :: Biatx
!        real(8), dimension(k) :: locBiatx
!        real(8), dimension(dimx) :: x
!        real(8) :: xmin
!        real(8) :: xmax
!        real(8) :: dx
!        real(8) :: error
!        integer :: left
!        integer :: mflag
!        integer :: iter
!        integer :: i
!        integer :: j
!        REAL time_begin, time_end
!        real ( kind = 8),  dimension ( k,dimx) :: deltal
!        real ( kind = 8),  dimension ( k,dimx) :: deltar
!        real ( kind = 8), dimension(dimx) :: saved
!        real ( kind = 8), dimension(dimx) :: term
!#ifdef _OPENMP
!        integer, parameter :: CHUNK = 2
!        integer :: OMP_GET_THREAD_NUM
!        integer :: TID
!        integer :: NTHREADS
!        integer :: OMP_GET_NUM_THREADS
!#endif
!
!        CALL printlog("===== Begin test1 =====", ai_dtllevel = 0)
!        !        call printcputime()
!
!        !        t = (/ 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0 /)
!        !        xmin = 0.25
!        !        xmax = 0.5
!
!        t(1:k) = 0.0
!        t(k + 1:N + k) = 1.0
!        xmin = 0.0
!        xmax = 1.0
!
!        dx = (xmax - xmin)/float(dimx + 1)
!        do i = 1, dimx
!            x(i) = xmin + float(i) * dx
!        end do
!
!        call interv(t, N + 1, x(1), left, mflag)
!
!!$OMP PARALLEL SHARED(x, t, left,Biatx) PRIVATE(i,iter,TID,locBiatx,time_end,time_begin)
!#ifdef _OPENMP
!        TID = OMP_GET_THREAD_NUM()
!        IF (TID .EQ. 0) THEN
!            NTHREADS = OMP_GET_NUM_THREADS()
!        PRINT *, 'Number of threads =', NTHREADS
!        END IF
!!        PRINT *, 'Thread',TID,' starting...'
!#endif
!!$OMP SECTIONS
!
!!$OMP SECTION
!        ! ****************
!!        CALL printlog(" using bsplvb ", ai_dtllevel = 1)
!#ifdef _OPENMP
!!        PRINT *, 'Thread',TID, ' using bsplvb'
!#endif
!        CALL CPU_TIME(time_begin)
!!$OMP DO
!        do iter = 1, Niter
!            do i = 1, dimx
!                CALL bsplvb(t, k, 1, x(i), left, locBiatx(:))
!                Biatx(:, i) = locBiatx
!!                print *, "----"
!!                PRINT *, 'Thread',TID, ' with i=', i, ' and x =', x(i)
!!                PRINT *, 'locBiatx =', locBiatx
!!                print *, "----"
!            end do
!        end do
!!$OMP END DO NOWAIT
!
!        CALL CPU_TIME(time_end)
!        PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds using bsplvb'
!
!        !        call printcputime()
!        ! ****************
!
!!$OMP SECTION
!        ! ****************
!!        CALL printlog(" using bsplvb_vect ", ai_dtllevel = 1)
!#ifdef _OPENMP
!!        PRINT *, 'Thread',TID, ' using bsplvb_vect'
!#endif
!        CALL CPU_TIME(time_begin)
!
!!$OMP DO
!        do iter = 1, Niter
!            CALL bsplvb_vect(t, k, 1, dimx, x, left, B_vectiatx, deltal, deltar, saved, term)
!        end do
!!$OMP END DO NOWAIT
!
!        CALL CPU_TIME(time_end)
!        PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds using bsplvb_vect'
!        !        call printcputime()
!        ! ****************
!
!        error = 0.0
!        do i = 1, dimx
!            do j = 1, k
!                error = error + (B_vectiatx(j, i) - Biatx(j, i))**2
!            end do
!        end do
!        error = SQRT(error) / (dimx * k)
!        print *, "error=", error
!!$OMP END SECTIONS NOWAIT
!
!#ifdef _OPENMP
!        IF (TID .EQ. 0) THEN
!            CALL printlog("===== End test1 =====", ai_dtllevel = 0)
!        end if
!
!!        PRINT *, 'Thread',TID,' done.'
!#endif
!!$OMP END PARALLEL
!
!    end subroutine test1
    !---------------------------------------------------------------
    subroutine test_bsplvb_vect()
        implicit none
        ! LOCAL
        ! LOCAL VARIABLES
        !        integer, parameter :: k = 4
        !        integer, parameter :: N = 7
        
        integer, parameter :: k = 4
        integer, parameter :: N = k
        integer, parameter :: Niter = 200000
        integer, parameter :: dimx = 1000
        real(8), dimension(N + k) :: t
        real(8), dimension(k, dimx) :: B_vectiatx
        real(8), dimension(k, dimx) :: Biatx
        real(8), dimension(k) :: locBiatx
        real(8), dimension(dimx) :: x
        real(8) :: xmin
        real(8) :: xmax
        real(8) :: dx
        real(8) :: error
        integer :: left
        integer :: mflag
        integer :: iter
        integer :: i
        integer :: j
        REAL time_begin, time_end
        real ( kind = 8),  dimension ( k,dimx) :: deltal
        real ( kind = 8),  dimension ( k,dimx) :: deltar
        real ( kind = 8), dimension(dimx) :: saved
        real ( kind = 8), dimension(dimx) :: term
#ifdef _OPENMP
        integer, parameter :: CHUNK = 2
        integer :: OMP_GET_THREAD_NUM
        integer :: TID
        integer :: NTHREADS
        integer :: OMP_GET_NUM_THREADS
#endif

        CALL printlog("===== Begin test_bsplvb_vect =====", ai_dtllevel = 0)
        !        call printcputime()

        !        t = (/ 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0 /)
        !        xmin = 0.25
        !        xmax = 0.5

        t(1:k) = 0.0
        t(k + 1:N + k) = 1.0
        xmin = 0.0
        xmax = 1.0

        dx = (xmax - xmin)/float(dimx + 1)
        do i = 1, dimx
            x(i) = xmin + float(i) * dx
        end do

        call interv(t, N + 1, x(1), left, mflag)

!$OMP PARALLEL SHARED(x, t, left,Biatx) PRIVATE(i,iter,TID,locBiatx,time_end,time_begin)
#ifdef _OPENMP
        TID = OMP_GET_THREAD_NUM()
        IF (TID .EQ. 0) THEN
            NTHREADS = OMP_GET_NUM_THREADS()
        PRINT *, 'Number of threads =', NTHREADS
        END IF
!        PRINT *, 'Thread',TID,' starting...'
#endif
!$OMP SECTIONS

!$OMP SECTION
        ! ****************
!        CALL printlog(" using bsplvb_vect ", ai_dtllevel = 1)
#ifdef _OPENMP
!        PRINT *, 'Thread',TID, ' using bsplvb_vect'
#endif
        CALL CPU_TIME(time_begin)

!$OMP DO
        do iter = 1, Niter
            CALL bsplvb_vect(t, k, 1, dimx, x, left, B_vectiatx, deltal, deltar, saved, term)
        end do
!$OMP END DO NOWAIT

        CALL CPU_TIME(time_end)
        PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds using bsplvb_vect'
        !        call printcputime()
        ! ****************
!$OMP END SECTIONS NOWAIT

#ifdef _OPENMP
        IF (TID .EQ. 0) THEN
            CALL printlog("===== End test_bsplvb_vect =====", ai_dtllevel = 0)
        end if

!        PRINT *, 'Thread',TID,' done.'
#endif
!$OMP END PARALLEL

    end subroutine test_bsplvb_vect
    !---------------------------------------------------------------
    subroutine test_bsplvb()
        implicit none
        ! LOCAL
        ! LOCAL VARIABLES
        !        integer, parameter :: k = 4
        !        integer, parameter :: N = 7
        
        integer, parameter :: k = 4
        integer, parameter :: N = k
        integer, parameter :: Niter = 200000
        integer, parameter :: dimx = 1000
        real(8), dimension(N + k) :: t
        real(8), dimension(k, dimx) :: B_vectiatx
        real(8), dimension(k, dimx) :: Biatx
        real(8), dimension(k) :: locBiatx
        real(8), dimension(dimx) :: x
        real(8) :: xmin
        real(8) :: xmax
        real(8) :: dx
        real(8) :: error
        integer :: left
        integer :: mflag
        integer :: iter
        integer :: i
        integer :: j
        REAL time_begin, time_end
        real ( kind = 8),  dimension ( k,dimx) :: deltal
        real ( kind = 8),  dimension ( k,dimx) :: deltar
        real ( kind = 8), dimension(dimx) :: saved
        real ( kind = 8), dimension(dimx) :: term
#ifdef _OPENMP
        integer, parameter :: CHUNK = 2
        integer :: OMP_GET_THREAD_NUM
        integer :: TID
        integer :: NTHREADS
        integer :: OMP_GET_NUM_THREADS
#endif

        CALL printlog("===== Begin test_bsplvb =====", ai_dtllevel = 0)
        !        call printcputime()

        !        t = (/ 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0 /)
        !        xmin = 0.25
        !        xmax = 0.5

        t(1:k) = 0.0
        t(k + 1:N + k) = 1.0
        xmin = 0.0
        xmax = 1.0

        dx = (xmax - xmin)/float(dimx + 1)
        do i = 1, dimx
            x(i) = xmin + float(i) * dx
        end do

        call interv(t, N + 1, x(1), left, mflag)

!$OMP PARALLEL SHARED(x, t, left,Biatx) PRIVATE(i,iter,TID,locBiatx,time_end,time_begin)
#ifdef _OPENMP
        TID = OMP_GET_THREAD_NUM()
        IF (TID .EQ. 0) THEN
            NTHREADS = OMP_GET_NUM_THREADS()
        PRINT *, 'Number of threads =', NTHREADS
        END IF
!        PRINT *, 'Thread',TID,' starting...'
#endif
!$OMP SECTIONS

!$OMP SECTION
        ! ****************
!        CALL printlog(" using bsplvb ", ai_dtllevel = 1)
#ifdef _OPENMP
!        PRINT *, 'Thread',TID, ' using bsplvb'
#endif
        CALL CPU_TIME(time_begin)
!$OMP DO
        do iter = 1, Niter
            do i = 1, dimx
                CALL bsplvb(t, k, 1, x(i), left, locBiatx(:))
                Biatx(:, i) = locBiatx
!                print *, "----"
!                PRINT *, 'Thread',TID, ' with i=', i, ' and x =', x(i)
!                PRINT *, 'locBiatx =', locBiatx
!                print *, "----"
            end do
        end do
!$OMP END DO NOWAIT

        CALL CPU_TIME(time_end)
        PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds using bsplvb'

        !        call printcputime()
        ! ****************
!$OMP END SECTIONS NOWAIT

#ifdef _OPENMP
        IF (TID .EQ. 0) THEN
            CALL printlog("===== End test_bsplvb =====", ai_dtllevel = 0)
        end if

!        PRINT *, 'Thread',TID,' done.'
#endif
!$OMP END PARALLEL

    end subroutine test_bsplvb
    !---------------------------------------------------------------
    subroutine test11()
        implicit none
        ! LOCAL
        ! LOCAL VARIABLES
        !        integer, parameter :: k = 4
        !        integer, parameter :: N = 7

        integer, parameter :: k = 4
        integer, parameter :: N = k
        integer, parameter :: nderiv = 3
        integer, parameter :: Niter = 10000
        integer, parameter :: dimx = 1000
        real(8), dimension(N + k) :: t
        real(8), dimension(k,nderiv, dimx) :: dB_vectiatx
        real(8), dimension(k,nderiv, dimx) :: dBiatx
        real(8), dimension(k,nderiv) :: locdBiatx
        real(8), dimension(dimx) :: x
        real ( kind = 8 ) a(k,k)
        real(8) :: xmin
        real(8) :: xmax
        real(8) :: dx
        real(8) :: error
        integer :: left
        integer :: mflag
        integer :: iter
        integer :: i
        integer :: j
        integer :: d
        REAL time_begin, time_end
        real ( kind = 8),  dimension ( k,dimx) :: deltal
        real ( kind = 8),  dimension ( k,dimx) :: deltar
        real ( kind = 8), dimension(dimx) :: saved
        real ( kind = 8), dimension(dimx) :: term

        CALL printlog("===== Begin test11 =====", ai_dtllevel = 0)
        !        call printcputime()

        !        t = (/ 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0 /)
        !        xmin = 0.25
        !        xmax = 0.5

        t(1:k) = 0.0
        t(k + 1:N + k) = 1.0
        xmin = 0.0
        xmax = 1.0

        dx = (xmax - xmin)/float(dimx + 1)
        do i = 1, dimx
            x(i) = xmin + float(i) * dx
        end do

        call interv(t, N + 1, x(1), left, mflag)

        print *, '****************'
        print *, 'Begin bsplvd'
        print *, '****************'
        CALL CPU_TIME(time_begin)

        do iter = 1, Niter
            do i = 1, dimx
                CALL bsplvd(t, k, x(i), left, a, locdBiatx, nderiv )
                dBiatx(:, :,i) = locdBiatx
            end do
        end do

        CALL CPU_TIME(time_end)
        PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds using bsplvd'

        !        call printcputime()
        ! ****************

        ! ****************
        print *, '****************'
        print *, 'Begin bsplvd_vect'
        print *, '****************'
        CALL CPU_TIME(time_begin)

        do iter = 1, Niter
            CALL bsplvd_vect(t, k, x, dimx, left, a, dB_vectiatx, nderiv, deltal, deltar, saved, term)
        end do

        CALL CPU_TIME(time_end)
        PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds using bsplvd_vect'
        !        call printcputime()
        ! ****************

        error = 0.0
        do i = 1, dimx
            do j = 1, k
                do d = 1, nderiv
                    error = error + (dB_vectiatx(j,d,i) - dBiatx(j,d,i))**2
                end do
            end do
        end do
        error = SQRT(error) / (dimx * k * nderiv)
        print *, "error=", error

!        print *, "dBiatx ="
!        do d = 1, nderiv
!            print *,dBiatx(:,d,:)
!        end do
!
!        print *, "dB_vectiatx ="
!        do d = 1, nderiv
!            print *,dB_vectiatx(:,d,:)
!        end do

        CALL printlog("===== End test11 =====", ai_dtllevel = 0)


    end subroutine test11
    !---------------------------------------------------------------
    subroutine test_bsplvd_vect()
        implicit none
        ! LOCAL
        ! LOCAL VARIABLES
        !        integer, parameter :: k = 4
        !        integer, parameter :: N = 7

        integer, parameter :: k = 4
        integer, parameter :: N = k
        integer, parameter :: nderiv = 2
        integer, parameter :: Niter = 10000
        integer, parameter :: dimx = 1000
        real(8), dimension(N + k) :: t
        real(8), dimension(nderiv,k,dimx) :: dB_vectiatx
        real(8), dimension(dimx) :: x
        real ( kind = 8 ) a(k,k)
        real(8) :: xmin
        real(8) :: xmax
        real(8) :: dx
        real(8) :: error
        integer :: left
        integer :: mflag
        integer :: iter
        integer :: i
        integer :: j
        integer :: d
        REAL time_begin, time_end
        real ( kind = 8),  dimension ( k,dimx) :: deltal
        real ( kind = 8),  dimension ( k,dimx) :: deltar
        real ( kind = 8), dimension(dimx) :: saved
        real ( kind = 8), dimension(dimx) :: term

        print *,"===== Begin test_bsplvd_vect ====="

        t(1:k) = 0.0
        t(k + 1:N + k) = 1.0
        xmin = 0.0
        xmax = 1.0

        dx = (xmax - xmin)/float(dimx + 1)
        do i = 1, dimx
            x(i) = xmin + float(i) * dx
        end do

        call interv(t, N + 1, x(1), left, mflag)

        print *, '****************'
        print *, 'Begin bsplvd_vect'
        print *, '****************'
        CALL CPU_TIME(time_begin)

        do iter = 1, Niter
            CALL bsplvd_vect_opt(t, k, x, dimx, left, a, dB_vectiatx, nderiv, deltal, deltar, saved, term)
        end do

        CALL CPU_TIME(time_end)

        PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds using bsplvd_vect'
        !        call printcputime()
        ! ****************

    end subroutine test_bsplvd_vect
    !---------------------------------------------------------------
    subroutine test_bsplvd()
        implicit none
        ! LOCAL
        ! LOCAL VARIABLES
        !        integer, parameter :: k = 4
        !        integer, parameter :: N = 7

        integer, parameter :: k = 4
        integer, parameter :: N = k
        integer, parameter :: nderiv = 3
        integer, parameter :: Niter = 10000
        integer, parameter :: dimx = 1000
        real(8), dimension(N + k) :: t
        real(8), dimension(k,nderiv, dimx) :: dBiatx
        real(8), dimension(k,nderiv) :: locdBiatx
        real(8), dimension(dimx) :: x
        real ( kind = 8 ) a(k,k)
        real(8) :: xmin
        real(8) :: xmax
        real(8) :: dx
        real(8) :: error
        integer :: left
        integer :: mflag
        integer :: iter
        integer :: i
        integer :: j
        integer :: d
        REAL time_begin, time_end
        real ( kind = 8),  dimension ( k,dimx) :: deltal
        real ( kind = 8),  dimension ( k,dimx) :: deltar
        real ( kind = 8), dimension(dimx) :: saved
        real ( kind = 8), dimension(dimx) :: term

        CALL printlog("===== Begin test_bsplvd =====", ai_dtllevel = 0)
        !        call printcputime()

        !        t = (/ 0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0 /)
        !        xmin = 0.25
        !        xmax = 0.5

        t(1:k) = 0.0
        t(k + 1:N + k) = 1.0
        xmin = 0.0
        xmax = 1.0

        dx = (xmax - xmin)/float(dimx + 1)
        do i = 1, dimx
            x(i) = xmin + float(i) * dx
        end do

        call interv(t, N + 1, x(1), left, mflag)

        print *, '****************'
        print *, 'Begin bsplvd'
        print *, '****************'
        CALL CPU_TIME(time_begin)

        do iter = 1, Niter
            do i = 1, dimx
                CALL bsplvd(t, k, x(i), left, a, locdBiatx, nderiv )
                dBiatx(:, :,i) = locdBiatx
            end do
        end do

        CALL CPU_TIME(time_end)
        PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds using bsplvd'

        !        call printcputime()
        ! ****************

        CALL printlog("===== End test_bsplvd =====", ai_dtllevel = 0)


    end subroutine test_bsplvd

end module tests