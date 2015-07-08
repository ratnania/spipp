!     
! File:   bsplvd_vect.F90
! Author: ratnani
!
! Created on December 12, 2011, 9:55 AM
!

subroutine bsplvd_vect(t, k, x, dimx, left, a, dbiatx, nderiv, deltal, deltar, saved, term)

    !*****************************************************************************80
    !
    !! BSPLVD_VECT calculates the nonvanishing B-splines and derivatives at at a grid point X.
    !
    !    see BSPLVD for more details on the algorithm

    implicit none

    integer k
    integer left
    integer nderiv
    integer dimx

    real(wp) a(:,:)
    real(wp) ta(k,k)
    real(wp) dbiatx(:,:,:)
    real(wp) factor
    real(wp) fkp1mm
    integer i
    integer ideriv
    integer il
    integer j
    integer d
    integer jlow
    integer jp1mid
    integer ldummy
    integer m
    integer mhigh
    real(wp) sum1
    real(wp) t(left + k)
    real(wp), dimension(:) :: x
    real(wp), dimension (:,:) :: deltal
    real(wp), dimension (:,:) :: deltar
    real(wp), dimension(:) :: saved
    real(wp), dimension(:) :: term

    mhigh = max(min(nderiv, k), 1)

    dbiatx = 0.0_wp

    !
    !  MHIGH is usually equal to NDERIV.
    !
    call bsplvb_vect(t, k + 1 - mhigh, 1, dimx, x, left, dbiatx(:,1,:), deltal, deltar, saved, term)

    if (mhigh == 1) then
        return
    end if
    !
    !  The first column of DBIATX always contains the B-spline values
    !  for the current order.  These are stored in column K+1-current
    !  order before BSPLVB is called to put values for the next
    !  higher order on top of it.
    !
    ideriv = mhigh
    do m = 2, mhigh
        jp1mid = 1
        do j = ideriv, k
            dbiatx(j, ideriv,:) = dbiatx(jp1mid, 1,:)
            jp1mid = jp1mid + 1
        end do
        ideriv = ideriv - 1
        call bsplvb_vect(t, k + 1 - ideriv, 2, dimx, x, left, dbiatx(:,1,:), deltal, deltar, saved, term)
    end do
    !
    !  At this point, B(LEFT-K+I, K+1-J)(X) is in DBIATX(I,J) for
    !  I=J,...,K and J=1,...,MHIGH ('=' NDERIV).
    !
    !  In particular, the first column of DBIATX is already in final form.
    !
    !  To obtain corresponding derivatives of B-splines in subsequent columns,
    !  generate their B-representation by differencing, then evaluate at X.
    !
    jlow = 1
    do i = 1, k
        a(jlow:k, i) = 0.0_wp
        jlow = i
        a(i, i) = 1.0_wp
    end do
    !
    !  At this point, A(.,J) contains the B-coefficients for the J-th of the
    !  K B-splines of interest here.
    !
    do m = 2, mhigh

        fkp1mm = real ( k + 1 - m, kind = 8)
        il = left
        i = k
        !
        !  For J = 1,...,K, construct B-coefficients of (M-1)st derivative of
        !  B-splines from those for preceding derivative by differencing
        !  and store again in  A(.,J).  The fact that  A(I,J) = 0 for
        !  I < J is used.
        !
        do ldummy = 1, k + 1 - m

            factor = fkp1mm / (t(il + k + 1 - m) - t(il))
            !
            !  The assumption that T(LEFT) < T(LEFT+1) makes denominator
            !  in FACTOR nonzero.
            !
            a(i, 1:i) = (a(i, 1:i) - a(i - 1, 1:i)) * factor

            il = il - 1
            i = i - 1

        end do
        ta = transpose(a)
        !
        !  For I = 1,...,K, combine B-coefficients A(.,I) with B-spline values
        !  stored in DBIATX(.,M) to get value of (M-1)st derivative of
        !  I-th B-spline (of interest here) at X, and store in DBIATX(I,M).
        !
        !  Storage of this value over the value of a B-spline
        !  of order M there is safe since the remaining B-spline derivatives
        !  of the same order do not use this value due to the fact
        !  that  A(J,I) = 0  for J < I.
        !
        do i = 1, k

            jlow = max(i, m)
!            do d=1, dimx
!                dbiatx(i,m,d) = dot_product ( a(jlow:k,i), dbiatx(jlow:k,m,d) )
!            end do

            dbiatx(i,m,:) = MATMUL ( ta(i,jlow:k), dbiatx(jlow:k,m,:) )

        end do

    end do

    return
end subroutine bsplvd_vect


!
! File:   bsplvd_vect.F90
! Author: ratnani
!
! Created on December 12, 2011, 9:55 AM
!

subroutine bsplvd_vect_opt(t, k, x, dimx, left, a, dbiatx, nderiv, deltal, deltar, saved, term)

    !*****************************************************************************80
    !
    !! BSPLVD_VECT calculates the nonvanishing B-splines and derivatives at at a grid point X.
    !
    !    see BSPLVD for more details on the algorithm

    implicit none

    integer k
    integer left
    integer nderiv
    integer dimx

    real(wp) a(:,:)
    real(wp) ta(k,k)
    real(wp) dbiatx(:,:,:)
    real(wp) factor
    real(wp) fkp1mm
    integer i
    integer ideriv
    integer il
    integer j
    integer d
    integer jlow
    integer jp1mid
    integer ldummy
    integer m
    integer mhigh
    real(wp) sum1
    real(wp) t(left + k)
    real(wp), dimension(:) :: x
    real(wp), dimension (:,:) :: deltal
    real(wp), dimension (:,:) :: deltar
    real(wp), dimension(:) :: saved
    real(wp), dimension(:) :: term

    mhigh = max(min(nderiv, k), 1)

    dbiatx = 0.0_wp

    !
    !  MHIGH is usually equal to NDERIV.
    !
    call bsplvb_vect(t, k + 1 - mhigh, 1, dimx, x, left, dbiatx(1,:,:), deltal, deltar, saved, term)

    if (mhigh == 1) then
        return
    end if
    !
    !  The first column of DBIATX always contains the B-spline values
    !  for the current order.  These are stored in column K+1-current
    !  order before BSPLVB is called to put values for the next
    !  higher order on top of it.
    !
    ideriv = mhigh
    do m = 2, mhigh
        jp1mid = 1
        do j = ideriv, k
            dbiatx(ideriv,j,:) = dbiatx(1,jp1mid,:)
            jp1mid = jp1mid + 1
        end do
        ideriv = ideriv - 1
        call bsplvb_vect(t, k + 1 - ideriv, 2, dimx, x, left, dbiatx(1,:,:), deltal, deltar, saved, term)
    end do
    !
    !  At this point, B(LEFT-K+I, K+1-J)(X) is in DBIATX(I,J) for
    !  I=J,...,K and J=1,...,MHIGH ('=' NDERIV).
    !
    !  In particular, the first column of DBIATX is already in final form.
    !
    !  To obtain corresponding derivatives of B-splines in subsequent columns,
    !  generate their B-representation by differencing, then evaluate at X.
    !
    jlow = 1
    do i = 1, k
        a(jlow:k, i) = 0.0_wp
        jlow = i
        a(i, i) = 1.0_wp
    end do
    !
    !  At this point, A(.,J) contains the B-coefficients for the J-th of the
    !  K B-splines of interest here.
    !
    do m = 2, mhigh

        fkp1mm = real ( k + 1 - m, kind = 8)
        il = left
        i = k
        !
        !  For J = 1,...,K, construct B-coefficients of (M-1)st derivative of
        !  B-splines from those for preceding derivative by differencing
        !  and store again in  A(.,J).  The fact that  A(I,J) = 0 for
        !  I < J is used.
        !
        do ldummy = 1, k + 1 - m

            factor = fkp1mm / (t(il + k + 1 - m) - t(il))
            !
            !  The assumption that T(LEFT) < T(LEFT+1) makes denominator
            !  in FACTOR nonzero.
            !
            a(i, 1:i) = (a(i, 1:i) - a(i - 1, 1:i)) * factor

            il = il - 1
            i = i - 1

        end do
        ta = transpose(a)
        !
        !  For I = 1,...,K, combine B-coefficients A(.,I) with B-spline values
        !  stored in DBIATX(.,M) to get value of (M-1)st derivative of
        !  I-th B-spline (of interest here) at X, and store in DBIATX(I,M).
        !
        !  Storage of this value over the value of a B-spline
        !  of order M there is safe since the remaining B-spline derivatives
        !  of the same order do not use this value due to the fact
        !  that  A(J,I) = 0  for J < I.
        !
        do i = 1, k

            jlow = max(i, m)
!            do d=1, dimx
!                dbiatx(i,m,d) = dot_product ( a(jlow:k,i), dbiatx(jlow:k,m,d) )
!            end do

            dbiatx(m,i,:) = MATMUL ( ta(i,jlow:k), dbiatx(m,jlow:k,:) )

        end do

    end do

    return
end subroutine bsplvd_vect_opt

