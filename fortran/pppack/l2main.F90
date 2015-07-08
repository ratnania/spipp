!  main program for least-squares approximation by splines
!  from  * a practical guide to splines *  by c. de Boor (7 may 92)
!alls setdat,l2knts,l2appr(bsplvb,bchfac,bchslv),bsplpp(bsplvb*)
!     ,l2err(ppvalu(interv)),ppvalu*,newnot
!
!  the program, though ostensibly written for l2-approximation, is typ-
!  ical for programs constructing a pp approximation to a function gi-
!  ven in some sense. the subprogram  l 2 a p p r , for instance, could
!  easily be replaced by one carrying out interpolation or some other
!  form of approximation.
!
!******  i n p u t  ******
!  is expected in  s e t d a t  (quo vide), specifying both the data to
!  be approximated and the order and breakpoint sequence of the pp ap-
!  proximating function to be used. further,  s e t d a t  is expected
!  to  t e r m i n a t e  the run (for lack of further input or because
!   i c o u n t  has reached a critical value).
!     the number  n t i m e s  is read in in the main program. it speci
!  fies the number of passes through the knot improvement algorithm in
!  n e w n o t  to be made. also,  a d d b r k  is read in to specify
!  that, on the average, addbrk knots are to be added per pass through
!  newnot. for example,  addbrk = .34  would cause a knot to be added
!  every third pass (as long as  ntimes .lt. 50).
!
!******  p r i n t e d  o u t p u t  ******
!  is governed by the three print control hollerith strings
!  p r b c o  = 'on'  gives printout of b-spline coeffs. of approxim.
!  p r p c o  = 'on'  gives printout of pp repr. of approximation.
!  p r f u n  = 'on'  gives printout of approximation and error at
!                     every data point.
!  the order  k , the number of pieces  l, and the interior breakpoints
!  are always printed out as are (in l2err) the mean, mean square, and
!  maximum errors in the approximation.
!
      integer i,icount,ii,j,k,l,lbegin,lnew,ll,lpkmax,ltkmax,n,nt,ntau &
              ,ntimes,ntmax,on,prbco,prfun,prpco
      parameter (lpkmax=100,ntmax=200,ltkmax=2000)
      real addbrk,bcoef(lpkmax),break,coef,gtau,q(ltkmax),scrtch(ntmax) &
          ,t(ntmax),tau,totalw,weight
      common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
!     real addbrk,bcoef(100),break,coef,gtau,q(2000),scrtch(200)
!    *    ,t(200),tau,totalw,weight
!     common / data / ntau, tau(200),gtau(200),weight(200),totalw
!     common /data/ also occurs in setdat, l2appr and l2err. it is ment-
!     ioned here only because it might otherwise become undefined be-
!     tween calls to those subroutines.
      common /approx/ break(lpkmax),coef(ltkmax),l,k
!     common /approx/ break(100),coef(2000),l,k
!     common /approx/ also occurs in setdat and l2err.
      data on /'ON'/
!
      icount = 0
!        i c o u n t  provides communication with the data-input-and-
!     termination routine  s e t d a t . it is initialized to  0  to
!     signal to setdat when it is being called for the first time. after
!     that, it is up to setdat to use icount for keeping track of the
!     passes through setdat .
!
!     information about the function to be approximated and order and
!     breakpoint sequence of the approximating pp functions is gathered
!     by a
    1 call setdat(icount)
!
!     breakpoints are translated into knots, and the number  n  of
!     b-splines to be used is obtained by a
      call l2knts ( break, l, k, t, n )
!
!     the integer  n t i m e s  and the real  a d d b r k  are requested
!     as well as the print controls  p r b c o ,  p r p c o  and
!     p r f u n .  ntimes  passes  are made through the subroutine new-
!     not, with an increase of  addbrk  knots for every pass .
      print 600
  600 format(' ntimes,addbrk , prbco,prpco,prfun =? (i3,f10.5/3a2)')
      read 500,ntimes,addbrk,prbco,prpco,prfun
  500 format(i3,f10.5/3a2)
!
      lbegin = l
      nt = 0
!        the b-spline coeffs.  b c o e f  of the l2-approx. are obtain-
!        ed by a
   10    call l2appr ( t, n, k, q, scrtch, bcoef )
         if (prbco .eq. on)  print 609, (bcoef(i),i=1,n)
  609    format(//' b-spline coefficients'/(4e20.10))
!
!        convert the b-repr. of the approximation to pp repr.
         call bsplpp ( t, bcoef, n, k, q, break, coef, l )
         print 610, k, l, (break(ll),ll=2,l)
  610    format(//' approximation by splines of order',i3,' on ', &
               i3,' intervals. breakpoints -'/(4e20.10))
         if (prpco .ne. on)             go to 15
         print 611
  611    format(/' pp-representation for approximation')
         do 12 i=1,l
            ii = (i-1)*k
   12       print 613,break(i),(coef(ii+j),j=1,k)
  613    format(f9.3,4e20.10/(11x,4e20.10))
!
!        compute and print out various error norms by a
   15    call l2err ( prfun, scrtch, q )
!
!        if newnot has been applied less than  n t i m e s  times, try
!        it again to obtain, from the current approx. a possibly improv-
!        ed sequence of breakpoints with  addbrk  more breakpoints (on
!        the average) than the current approximation has.
!           if only an increase in breakpoints is wanted, without the
!        adjustment that newnot provides, a fake newnot routine could be
!        used here which merely returns the breakpoints for  l n e w
!        equal intervals .
         if (nt .ge. ntimes)            go to 1
         lnew = lbegin + float(nt)*addbrk
         call newnot (break, coef, l, k, scrtch, lnew, t )
         call l2knts ( scrtch, lnew, k, t, n )
         nt = nt + 1
                                        go to 10
      end
