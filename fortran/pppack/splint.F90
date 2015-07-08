subroutine splint ( tau, gtau, t, n, k, q, bcoef, iflag )

!*****************************************************************************80
!
!! SPLINT produces the B-spline coefficients BCOEF of an interpolating spline.
!
!  Discussion:
!
!    The spline is of order K with knots T(1:N+K), and takes on the
!    value GTAU(I) at TAU(I), for I = 1 to N.
!
!    The I-th equation of the linear system
!
!      A * BCOEF = B
!
!    for the B-spline coefficients of the interpolant enforces interpolation
!    at TAU(1:N).
!
!    Hence, B(I) = GTAU(I), for all I, and A is a band matrix with 2*K-1
!    bands, if it is invertible.
!
!    The matrix A is generated row by row and stored, diagonal by diagonal,
!    in the rows of the array Q, with the main diagonal going
!    into row K.  See comments in the program.
!
!    The banded system is then solved by a call to BANFAC, which
!    constructs the triangular factorization for A and stores it again in
!    Q, followed by a call to BANSLV, which then obtains the solution
!    BCOEF by substitution.
!
!    BANFAC does no pivoting, since the total positivity of the matrix
!    A makes this unnecessary.
!
!    The linear system to be solved is (theoretically) invertible if
!    and only if
!      T(I) < TAU(I) < TAU(I+K), for all I.
!    Violation of this condition is certain to lead to IFLAG = 2.
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
!    Input, real(wp) TAU(N), the data point abscissas.  The entries in
!    TAU should be strictly increasing.
!
!    Input, real(wp) GTAU(N), the data ordinates.
!
!    Input, real(wp) T(N+K), the knot sequence.
!
!    Input, integer N, the number of data points.
!
!    Input, integer K, the order of the spline.
!
!    Output, real(wp) Q((2*K-1)*N), the triangular factorization
!    of the coefficient matrix of the linear system for the B-coefficients
!    of the spline interpolant.  The B-coefficients for the interpolant
!    of an additional data set can be obtained without going through all
!    the calculations in this routine, simply by loading HTAU into BCOEF
!    and then executing the call:
!      call banslv ( q, 2*k-1, n, k-1, k-1, bcoef )
!
!    Output, real(wp) BCOEF(N), the B-spline coefficients of
!    the interpolant.
!
!    Output, integer IFLAG, error flag.
!    1, = success.
!    2, = failure.
!
  implicit none

  integer n

  real(wp) bcoef(n)
  real(wp) gtau(n)
  integer i
  integer iflag
  integer ilp1mx
  integer j
  integer jj
  integer k
  integer kpkm2
  integer left
  real(wp) q((2*k-1)*n)
  real(wp) t(n+k)
  real(wp) tau(n)
  real(wp) taui

  kpkm2 = 2 * ( k - 1 )
  left = k
  q(1:(2*k-1)*n) = 0.0_wp
!
!  Loop over I to construct the N interpolation equations.
!
  do i = 1, n

    taui = tau(i)
    ilp1mx = min ( i + k, n + 1 )
!
!  Find LEFT in the closed interval (I,I+K-1) such that
!
!    T(LEFT) <= TAU(I) < T(LEFT+1)
!
!  The matrix is singular if this is not possible.
!
    left = max ( left, i )

    if ( taui < t(left) ) then
      iflag = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINT - Fatal Error!'
      write ( *, '(a)' ) '  The linear system is not invertible!'
      return
    end if

    do while ( t(left+1) <= taui )

      left = left + 1

      if ( left < ilp1mx ) then
        cycle
      end if

      left = left - 1

      if ( t(left+1) < taui ) then
        iflag = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINT - Fatal Error!'
        write ( *, '(a)' ) '  The linear system is not invertible!'
        return
      end if

      exit

    end do
!
!  The I-th equation enforces interpolation at TAUI, hence for all J,
!
!    A(I,J) = B(J,K,T)(TAUI).
!
!  Only the K entries with J = LEFT-K+1,...,LEFT actually might be nonzero.
!
!  These K numbers are returned, in BCOEF (used for temporary storage here),
!  by the following.
!
    call bsplvb ( t, k, 1, taui, left, bcoef )
!
!  We therefore want BCOEF(J) = B(LEFT-K+J)(TAUI) to go into
!  A(I,LEFT-K+J), that is, into Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since
!  A(I+J,J) is to go into Q(I+K,J), for all I, J, if we consider Q
!  as a two-dimensional array, with  2*K-1 rows.  See comments in
!  BANFAC.
!
!  In the present program, we treat Q as an equivalent
!  one-dimensional array, because of fortran restrictions on
!  dimension statements.
!
!  We therefore want  BCOEF(J) to go into the entry of Q with index:
!
!    I -(LEFT+J)+2*K + ((LEFT+J)-K-1)*(2*K-1)
!   = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
!
    jj = i - left + 1 + ( left - k ) * ( k + k - 1 )

    do j = 1, k
      jj = jj + kpkm2
      q(jj) = bcoef(j)
    end do

  end do
!
!  Obtain factorization of A, stored again in Q.
!
  call banfac ( q, k+k-1, n, k-1, k-1, iflag )

  if ( iflag == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINT - Fatal Error!'
    write ( *, '(a)' ) '  The linear system is not invertible!'
    return
  end if
!
!  Solve
!
!    A * BCOEF = GTAU
!
!  by back substitution.
!
  bcoef(1:n) = gtau(1:n)

  call banslv ( q, k+k-1, n, k-1, k-1, bcoef )

  return
end subroutine splint
