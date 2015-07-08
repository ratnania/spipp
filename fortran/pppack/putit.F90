      subroutine putit ( t, kpm, left, scrtch, dbiatx, q, nrow, b )
!  from  * a practical guide to splines *  by c. de Boor (7 may 92)
!alls  bsplvd(bsplvb),difequ(*)
!  to be called by  e q b l o k .
!
!  puts together one block of the collocation equation system
!
!******  i n p u t  ******
!  t     knot sequence, of size  left+kpm (at least)
!  kpm   order of spline
!  left  integer indicating interval of interest, viz the interval
!           (t(left), t(left+1))
!  nrow  number of rows in block to be put together
!
!******  w o r k  a r e a  ******
!  scrtch   used in bsplvd, of size (kpm,kpm)
!  dbiatx   used to contain derivatives of b-splines, of size (kpm,m+1)
!           with dbiatx(j,i+1) containing the i-th derivative of the
!           j-th b-spline of interest
!
!******  o u t p u t  ******
!  q  the block, of size (nrow,kpm)
!  b  the corresponding piece of the right side, of size (nrow)
!
!******  m e t h o d  ******
!  the  k  collocation equations for the interval  (t(left),t(left+1))
!  are constructed with the aid of the subroutine  d i f e q u ( 2, .,
!  . ) and interspersed (in order) with the side conditions (if any) in
!  this interval, using  d i f e q u ( 3, ., . )  for the information.
!     the block  q  has  kpm  columns, corresponding to the  kpm  b-
!  splines of order  kpm  which have the interval (t(left),t(left+1))
!  in their support. the block's diagonal is part of the diagonal of the
!  total system. the first equation in this block not overlapped by the
!  preceding block is therefore equation  l o w r o w , with lowrow = 1+
!  number of side conditions in preceding intervals (or blocks).
!
      integer kpm,left,nrow,   i,irow,iside,itermx,j,k,ll,lowrow,m,mode &
                              ,mp1
      real b(1),dbiatx(kpm,1),q(nrow,kpm),scrtch(1),t(1),   dx,rho,sum &
                                                      ,v(20),xm,xside,xx
      common /side/ m, iside, xside(10)
      common /other/ itermx,k,rho(19)
      mp1 = m+1
      do 10 j=1,kpm
         do 10 i=1,nrow
   10       q(i,j) = 0.
      xm = (t(left+1)+t(left))/2.
      dx = (t(left+1)-t(left))/2.
!
      ll = 1
      lowrow = iside
      do 30 irow=lowrow,nrow
         if (ll .gt. k)                 go to 22
         mode = 2
!        next collocation point is ...
         xx = xm + dx*rho(ll)
         ll = ll + 1
!        the corresp.collocation equation is next unless the next side
!        condition occurs at a point at, or to the left of, the next
!        collocation point.
         if (iside .gt. m)              go to 24
         if (xside(iside) .gt. xx)      go to 24
         ll = ll - 1
   22    mode = 3
         xx = xside(iside)
   24    call difequ ( mode, xx, v )
!        the next equation, a collocation equation (mode = 2) or a side
!        condition (mode = 3), reads
!   (*)   (v(m+1)*d**m + v(m)*d**(m-1) +...+ v(1)*d**0)f(xx) = v(m+2)
!        in terms of the info supplied by  difequ . the corresponding
!        equation for the b-coeffs of  f  therefore has the left side of
!        (*), evaluated at each of the  kpm  b-splines having  xx  in
!        their support, as its  kpm  possibly nonzero coefficients.
         call bsplvd ( t, kpm, xx, left, scrtch, dbiatx, mp1 )
         do 26 j=1,kpm
            sum = 0.
            do 25 i=1,mp1
   25          sum = v(i)*dbiatx(j,i) + sum
   26       q(irow,j) = sum
   30    b(irow) = v(m+2)
                                        return
      end
