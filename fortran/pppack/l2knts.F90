      subroutine l2knts ( break, l, k, t, n )
!  from  * a practical guide to splines *  by c. de boor
!  to be called in main program  l 2 m a i n .
!onverts the breakpoint sequence  b r e a k   into a corresponding knot
!  sequence  t  to allow the repr. of a pp function of order  k  with
!  k-2 continuous derivatives as a spline of order  k  with knot
!  sequence  t . this means that
!  t(1), ..., t(n+k) =  break(1) k times, then break(i), i=2,...,l, each
!                       once, then break(l+1) k times .
!  therefore,  n = k-1 + l.
!
!******  i n p u t  ******
!  k     the order
!  l     the number of polynomial pieces
!  break(1), ...,break(l+1)  the breakpoint sequence
!
!******  o u t p u t  ******
!  t(1),...,t(n+k)   the knot sequence
!  n     the dimension of the corresp. spline space of order  k .
!
      integer k,l,n,   i,km1
      real break(1),t(1)
!     dimension break(l+1),t(n+k)
      km1 = k - 1
      do 5 i=1,km1
    5    t(i) = break(1)
      do 6 i=1,l
    6    t(km1+i) = break(i)
      n = km1 + l
      do 7 i=1,k
    7    t(n+i) = break(l+1)
                                        return
      end
