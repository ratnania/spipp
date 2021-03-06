subroutine bsplvd ( t, k, x, left, a, dbiatx, nderiv )

!*****************************************************************************80
!
!! BSPLVD calculates the nonvanishing B-splines and derivatives at X.
!
!  Discussion:
!
!    Values at X of all the relevant B-splines of order K:K+1-NDERIV
!    are generated via BSPLVB and stored temporarily in DBIATX.
!
!    Then the B-spline coefficients of the required derivatives
!    of the B-splines of interest are generated by differencing,
!    each from the preceding one of lower order, and combined with
!    the values of B-splines of corresponding order in DBIATX
!    to produce the desired values.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, real(wp) T(LEFT+K), the knot sequence.  It is assumed that
!    T(LEFT) < T(LEFT+1).  Also, the output is correct only if
!    T(LEFT) <= X <= T(LEFT+1).
!
!    Input, integer K, the order of the B-splines to be evaluated.
!
!    Input, real(wp) X, the point at which these values are sought.
!
!    Input, integer LEFT, indicates the left endpoint of the interval of
!    interest.  The K B-splines whose support contains the interval
!    ( T(LEFT), T(LEFT+1) ) are to be considered.
!
!    Workspace, real(wp) A(K,K).
!
!    Output, real(wp) DBIATX(K,NDERIV).  DBIATX(I,M) contains
!    the value of the (M-1)st derivative of the (LEFT-K+I)-th B-spline
!    of order K for knot sequence T, I=M,...,K, M=1,...,NDERIV.
!
!    Input, integer NDERIV, indicates that values of B-splines and their
!    derivatives up to but not including the NDERIV-th are asked for.
!
  implicit none

  integer k
  integer left
  integer nderiv

  real(wp) a(k,k)
  real(wp) dbiatx(k,nderiv)
  real(wp) factor
  real(wp) fkp1mm
  integer i
  integer ideriv
  integer il
  integer j
  integer jlow
  integer jp1mid
  integer ldummy
  integer m
  integer mhigh
  real(wp) sum1
  real(wp) t(left+k)
  real(wp) x

  mhigh = max ( min ( nderiv, k ), 1 )
!
!  MHIGH is usually equal to NDERIV.
!
  call bsplvb ( t, k+1-mhigh, 1, x, left, dbiatx )

  if ( mhigh == 1 ) then
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
      dbiatx(j,ideriv) = dbiatx(jp1mid,1)
      jp1mid = jp1mid + 1
    end do
    ideriv = ideriv - 1
    call bsplvb ( t, k+1-ideriv, 2, x, left, dbiatx )
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
    a(jlow:k,i) = 0.0_wp
    jlow = i
    a(i,i) = 1.0_wp
  end do
!
!  At this point, A(.,J) contains the B-coefficients for the J-th of the
!  K B-splines of interest here.
!
  do m = 2, mhigh

    fkp1mm = real ( k + 1 - m, kind = 8 )
    il = left
    i = k
!
!  For J = 1,...,K, construct B-coefficients of (M-1)st derivative of
!  B-splines from those for preceding derivative by differencing
!  and store again in  A(.,J).  The fact that  A(I,J) = 0 for
!  I < J is used.
!
    do ldummy = 1, k+1-m

      factor = fkp1mm / ( t(il+k+1-m) - t(il) )
!
!  The assumption that T(LEFT) < T(LEFT+1) makes denominator
!  in FACTOR nonzero.
!
      a(i,1:i) = ( a(i,1:i) - a(i-1,1:i) ) * factor

      il = il - 1
      i = i - 1

    end do
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

      jlow = max ( i, m )

      dbiatx(i,m) = dot_product ( a(jlow:k,i), dbiatx(jlow:k,m) )

    end do

  end do

  return
end subroutine bsplvd
