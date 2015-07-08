!     
! File:   bsplvb_vect.F90
! Author: root
!
! Created on December 7, 2011, 11:26 AM
!

subroutine bsplvb_vect(t, jhigh, index, dimx, x, left, biatx, deltal, deltar, saved, term)

    !*****************************************************************************80
    !
    !! BSPLVB_VECT evaluates B-splines at a grid point X with a given knot sequence.
    !
    !  see BSPLVB for more details on the algorithm
    !
    implicit none

    integer, parameter :: jmax = 20
    !    integer, parameter :: dimxmax = 100000

    integer :: jhigh
    integer :: dimx
    real(wp) :: t(left + jhigh)
    real(wp), dimension(:) :: x
    real(wp) :: biatx(:,:)
    real(wp), dimension (:,:) :: deltal
    real(wp), dimension (:,:) :: deltar
    real(wp), dimension(:) :: saved
    real(wp), dimension(:) :: term

    integer :: i
    integer :: index
    integer, save :: j = 1
    integer :: left


    if (index == 1) then
        j = 1
        biatx(1, 1:dimx) = 1.0_wp
        if (jhigh <= j) then
            return
        end if
    end if

    if (t(left + 1) <= t(left)) then
        print*, 'x=', x
        write ( *, '(a)') ' '
        write ( *, '(a)') 'BSPLVB - Fatal error!'
        write ( *, '(a)') '  It is required that T(LEFT) < T(LEFT+1).'
        write ( *, '(a,i8)') '  But LEFT = ', left
        write ( *, '(a,g14.6)') '  T(LEFT) =   ', t(left)
        write ( *, '(a,g14.6)') '  T(LEFT+1) = ', t(left + 1)
        stop
    end if

    do

        deltar(j, 1:dimx) = t(left + j) - x(1:dimx)
        deltal(j, 1:dimx) = x(1:dimx) - t(left + 1 - j)

        saved = 0.0_wp

        do i = 1, j
            term(1:dimx) = biatx(i, 1:dimx) / (deltar(i, 1:dimx) + deltal(j + 1 - i, 1:dimx))
            biatx(i, 1:dimx) = saved(1:dimx) + deltar(i, 1:dimx) * term(1:dimx)
            saved(1:dimx) = deltal(j + 1 - i, 1:dimx) * term(1:dimx)
        end do

        biatx(j + 1, 1:dimx) = saved(1:dimx)
        j = j + 1

        if (jhigh <= j) then
            exit
        end if

    end do

    return
end subroutine bsplvb_vect

