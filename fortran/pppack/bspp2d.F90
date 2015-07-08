      subroutine bspp2d ( t, bcoef, n, k, m, scrtch, break, coef)
!  from  * a practical guide to splines *  by c. de Boor
!alls  bsplvb
!  this is an extended version of  bsplpp  for use with tensor products --------
!
!onverts the b-representation  t, bcoef(.,j), n, k  of some spline into --------
!  its pp-representation  break, coef(j,.,.), l, k , j=1, ..., m  .     --------
!
!******  i n p u t  ******
!  t.....knot sequence, of length  n+k
!  bcoef(.,j) b-spline coefficient sequence, of length  n ,j=1,...,m    --------
!  n.....length of  bcoef  and  dimension of spline space  spline(k,t)
!  k.....order of the spline
!
!  w a r n i n g  . . .  the restriction   k .le. kmax (= 20)   is impo-
!        sed by the arbitrary dimension statement for  biatx  below, but
!        is  n o w h e r e   c h e c k e d   for.
!
!  m     number of data sets                                            ********
!
!******  w o r k   a r e a  ******
!  scrtch   of size  (k,k,m), needed to contain bcoeffs of a piece of   ********
!        the spline and its  k-1  derivatives   for each of the m sets  --------
!
!******  o u t p u t  ******
!  break.....breakpoint sequence, of length  l+1, contains (in increas-
!        ing order) the distinct points in the sequence  t(k),...,t(n+1)
!  coef(mm,.,.)  array of size (k,n), with  coef(mm,i,j) = (i-1)st der- ********
!        ivative of  mm-th  spline at break(j) from the right, mm=1,.,m ********
!
!******  m e t h o d  ******
!     for each breakpoint interval, the  k  relevant b-coeffs of the
!  spline are found and then differenced repeatedly to get the b-coeffs
!  of all the derivatives of the spline on that interval. the spline and
!  its first  k-1  derivatives are then evaluated at the left end point
!  of that interval, using  bsplvb  repeatedly to obtain the values of
!  all b-splines of the appropriate order at that point.
!
      integer k,m,n,   i,j,jp1,kmax,kmj,left,lsofar,mm
      parameter (kmax = 20)
      real bcoef(n,m),break(n+2-k),coef(m,k,n+1-k),t(n+k) &
                               ,scrtch(k,k,m),biatx(kmax),diff,fkmj,sum
!
      lsofar = 0
      break(1) = t(k)
      do 50 left=k,n
!                                find the next nontrivial knot interval.
         if (t(left+1) .eq. t(left))    go to 50
         lsofar = lsofar + 1
         break(lsofar+1) = t(left+1)
         if (k .gt. 1)                  go to 9
         do 5 mm=1,m
    5       coef(mm,1,lsofar) = bcoef(left,mm)
                                        go to 50
!        store the k b-spline coeff.s relevant to current knot interval
!                             in  scrtch(.,1) .
    9    do 10 i=1,k
            do 10 mm=1,m
   10          scrtch(i,1,mm) = bcoef(left-k+i,mm)
!
!        for j=1,...,k-1, compute the  k-j  b-spline coeff.s relevant to
!        current knot interval for the j-th derivative by differencing
!        those for the (j-1)st derivative, and store in scrtch(.,j+1) .
         do 20 jp1=2,k
            j = jp1 - 1
            kmj = k - j
            fkmj = float(kmj)
            do 20 i=1,kmj
               diff = (t(left+i) - t(left+i - kmj))/fkmj
               if (diff .le. 0.)         go to 20
               do 15 mm=1,m
   15             scrtch(i,jp1,mm) = &
                  (scrtch(i+1,j,mm) - scrtch(i,j,mm))/diff
   20          continue
!
!        for  j = 0, ..., k-1, find the values at  t(left)  of the  j+1
!        b-splines of order  j+1  whose support contains the current
!        knot interval from those of order  j  (in  biatx ), then comb-
!        ine with the b-spline coeff.s (in scrtch(.,k-j) ) found earlier
!        to compute the (k-j-1)st derivative at  t(left)  of the given
!        spline.
!           note. if the repeated calls to  bsplvb  are thought to gene-
!        rate too much overhead, then replace the first call by
!           biatx(1) = 1.
!        and the subsequent call by the statement
!           j = jp1 - 1
!        followed by a direct copy of the lines
!           deltar(j) = t(left+j) - x
!                  ......
!           biatx(j+1) = saved
!        from  bsplvb . deltal(kmax)  and  deltar(kmax)  would have to
!        appear in a dimension statement, of course.
!
         call bsplvb ( t, 1, 1, t(left), left, biatx )
         do 25 mm=1,m
   25       coef(mm,k,lsofar) = scrtch(1,k,mm)
         do 30 jp1=2,k
            call bsplvb ( t, jp1, 2, t(left), left, biatx )
            kmj = k+1 - jp1
            do 30 mm=1,m
               sum = 0.
               do 28 i=1,jp1
   28             sum = biatx(i)*scrtch(i,kmj,mm) + sum
   30          coef(mm,kmj,lsofar) = sum
   50    continue
                                        return
      end
