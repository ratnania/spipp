      subroutine colloc(aleft,aright,lbegin,iorder,ntimes,addbrk,relerr)
!  from  * a practical guide to splines *  by c. de boor
!hapter xv, example. solution of an ode by collocation.
!alls colpnt, difequ(ppvalu(interv)), knots, eqblok(putit(difequ*,
!     bsplvd(bsplvb)))), slvblk(various subprograms), bsplpp(bsplvb*),
!     newnot
!
!******  i n p u t  ******
!  aleft, aright  endpoints of interval of approximation
!  lbegin   initial number of polynomial pieces in the approximation.
!           a uniform breakpoint sequence is chosen.
!  iorder   order of polynomial pieces in the approximation
!  ntimes   number of passes through  n e w n o t  to be made
!  addbrk   the number (possibly fractional) of breaks to be added per
!           pass through newnot. e.g., if addbrk = .33334, then a break-
!           point will be added at every third pass through newnot.
!  relerr   a tolerance. newton iteration is stopped if the difference
!           between the b-coeffs of two successive iterates is no more
!           than  relerr*(absol.largest b-coefficient).
!
!******  p r i n t e d   o u t p u t  ******
!     consists of the pp-representation of the approximate solution,
!     and of the error at selected points.
!
!******  m e t h o d  ******
!  the m-th order ordinary differential equation with  m  side condit-
!  ions, to be specified in subroutine  d i f e q u , is solved approx-
!  imately by collocation.
!  the approximation  f  to the solution  g  is pp of order  k+m  with
!  l  pieces and  m-1 continuous derivatives.  f  is determined by the
!  requirement that it satisfy the d.e. at  k  points per interval (to
!  be specified in  c o l p n t ) and the  m  side conditions.
!     this usually nonlinear system of equations for  f  is solved by
!  newton's method. the resulting linear system for the b-coeffs of an
!  iterate is constructed appropriately in  e q b l o k  and then solved
!  in  s l v b l k , a program designed to solve  a l m o s t  b l o c k
!  d i a g o n a l  linear systems efficiently.
!     there is an opportunity to attempt improvement of the breakpoint
!  sequence (both in number and location) through use of  n e w n o t .
!
          implicit none
          integer j
      integer npiece
      parameter (npiece=100)
      integer iorder,lbegin,ntimes,   i,iflag,ii,integs(3,npiece),iside &
                        ,iter,itermx,k,kmax,kpm,l,lenblk,lnew,m,n,nbloks &
                        ,ndim,ncoef,nncoef,nt
      parameter (ndim=200,kmax=20,ncoef=npiece*kmax,lenblk=ncoef)
      real(8) addbrk,aleft,aright,relerr,   a(ndim),amax,asave(ndim) &
           ,b(ndim),bloks(lenblk),break,coef,dx,err,rho,t(ndim) &
           ,templ(lenblk),temps(ndim),xside
      equivalence (bloks,templ)
      common /approx/ break(npiece), coef(ncoef), l,kpm
      common /side/ m, iside, xside(10)
      common /other/ itermx,k,rho(kmax-1)
!
      kpm = iorder
      if (lbegin*kpm .gt. ncoef)        go to 999
!  *** set the various parameters concerning the particular dif.equ.
!     including a first approx. in case the de is to be solved by
!     iteration ( itermx .gt. 0) .
      call difequ (1, temps(1), temps )
!  *** obtain the  k  collocation points for the standard interval.
      k = kpm - m
      call colpnt ( k, rho )
!  *** the following five statements could be replaced by a read in or-
!     der to obtain a specific (nonuniform) spacing of the breakpnts.
      dx = (aright - aleft)/float(lbegin)
      temps(1) = aleft
      do 4 i=2,lbegin
    4    temps(i) = temps(i-1) + dx
      temps(lbegin+1) = aright
!  *** generate, in knots, the required knots t(1),...,t(n+kpm).
      call knots ( temps, lbegin, kpm, t, n )
      nt = 1
!  *** generate the almost block diagonal coefficient matrix  bloks  and
!     right side  b  from collocation equations and side conditions.
!     then solve via  slvblk , obtaining the b-representation of the ap-
!     proximation in  t , a ,  n  , kpm  .
   10    call eqblok(t,n,kpm,temps,a,bloks,lenblk,integs,nbloks,b)
         call slvblk(bloks,integs,nbloks,b,temps,a,iflag)
         iter = 1
         if (itermx .le. 1)             go to 30
!  *** save b-spline coeff. of current approx. in  asave , then get new
!     approx. and compare with old. if coeff. are more than  relerr
!     apart (relatively) or if no. of iterations is less than  itermx ,
!     continue iterating.
   20       call bsplpp(t,a,n,kpm,templ,break,coef,l)
            do 25 i=1,n
   25          asave(i) = a(i)
            call eqblok(t,n,kpm,temps,a,bloks,lenblk,integs,nbloks,b)
            call slvblk(bloks,integs,nbloks,b,temps,a,iflag)
            err = 0.
            amax = 0.
            do 26 i=1,n
               amax = amax1(amax,abs(a(i)))
   26          err = amax1(err,abs(a(i)-asave(i)))
            if (err .le. relerr*amax)   go to 30
            iter = iter+1
            if (iter .lt. itermx)       go to 20
!  *** iteration (if any) completed. print out approx. based on current
!     breakpoint sequence, then try to improve the sequence.
   30    print 630,kpm,l,n,(break(i),i=2,l)
  630    format(47h approximation from a space of splines of order,i3 &
               ,4h on ,i3,11h intervals,/13h of dimension,i4 &
               ,16h.  breakpoints -/(5e20.10))
         if (itermx .gt. 0)  print 635,iter,itermx
  635    format(6h after,i3,3h of,i3,20h allowed iterations,)
         call bsplpp(t,a,n,kpm,templ,break,coef,l)
         print 637
  637    format(46h the pp representation of the approximation is)
         do 38 i=1,l
            ii = (i-1)*kpm
   38       print 638, break(i),(coef(ii+j),j=1,kpm)
  638    format(f9.3,e13.6,10e11.3)
!  *** the following call is provided here for possible further analysis
!     of the approximation specific to the problem being solved.
!     it is, of course, easily omitted.
         call difequ ( 4, temps(1), temps )
!
         if (nt .gt. ntimes)            return
!  *** from the pp-rep. of the current approx., obtain in  newnot  a new
!     (and possibly better) sequence of breakpoints, adding (on the ave-
!     rage)  a d d b r k  breakpoints per pass through newnot.
         lnew = lbegin + ifix(float(nt)*real(addbrk))
         if (lnew*kpm .gt. ncoef)       go to 999
         call newnot(break,coef,l,kpm,temps,lnew,templ)
         call knots(temps,lnew,kpm,t,n)
         nt = nt + 1
                                        go to 10
  999 nncoef = ncoef
      print 699,nncoef
  699 format(11h **********/23h the assigned dimension,i5 &
            ,25h for  coef  is too small.)
                                        return
      end subroutine colloc
