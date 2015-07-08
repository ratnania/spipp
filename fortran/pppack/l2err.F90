      subroutine l2err ( prfun , ftau , error )
!  from  * a practical guide to splines *  by c. de boor
!  this routine is to be called in the main program  l 2 m a i n .
!alls subprogram  ppvalu(interv)
!  this subroutine computes various errors of the current l2-approxi-
!  mation , whose pp-repr. is contained in common block  approx  ,
!  to the given data contained in common block  data . it prints out
!  the average error  e r r l 1 , the l2-error  e r r l 2,  and the
!  maximum error  e r r m a x .
!
!******  i n p u t  ******
!  prfun  a hollerith string.  if prfun = 'ON', the routine prints out
!          the value of the approximation as well as its error at
!          every data point.
!
!******  o u t p u t  ******
!  ftau(1), ..., ftau(ntau),  with  ftau(i)  the approximation  f at
!          tau(i), all i.
!  error(1), ..., error(ntau),  with  error(i) = scale*(g - f)
!          at tau(i), all i. here,  s c a l e  equals  1. in case
!          prfun .ne. 'ON' , or the abs.error is greater than 100 some-
!          where. otherwise, s c a l e  is such that the maximum of
!          abs(error))  over all  i  lies between  10  and  100. this
!          makes the printed output more illustrative.
!
      integer prfun,   ie,k,l,ll,lpkmax,ltkmax,ntau,ntmax,on
      real ftau(ntau),error(ntau),  break,coef,err,errmax,errl1,errl2 &
                                   ,gtau,scale,tau,totalw,weight
      parameter (lpkmax=100,ntmax=200,ltkmax=2000)
      common / data / ntau, tau(ntmax),gtau(ntmax),weight(ntmax),totalw
      common /approx/ break(lpkmax),coef(ltkmax),l,k
!
      data on /'ON'/
      errl1 = 0.
      errl2 = 0.
      errmax = 0.
      do 10 ll=1,ntau
         ftau(ll) = ppvalu (break, coef, l, k, tau(ll), 0 )
         error(ll) = gtau(ll) - ftau(ll)
         err = abs(error(ll))
         if (errmax .lt. err)   errmax = err
         errl1 = errl1 + err*weight(ll)
   10    errl2 = errl2 + err**2*weight(ll)
      errl1 = errl1/totalw
      errl2 = sqrt(errl2/totalw)
      print 615,errl2,errl1,errmax
  615 format(///' least square error =',e20.6/ &
                ' average error      =',e20.6/ &
                ' maximum error      =',e20.6//)
      if (prfun .ne. on)                return
!     **  scale error curve and print  **
      ie = 0
      scale = 1.
      if (errmax .ge. 10.)              go to 18
      do 17 ie=1,9
         scale = scale*10.
         if (errmax*scale .ge. 10.)     go to 18
   17    continue
   18 do 19 ll=1,ntau
   19    error(ll) = error(ll)*scale
      print 620,ie,(ll,tau(ll),ftau(ll),error(ll),ll=1,ntau)
  620 format(///14x,'approximation and scaled error curve'/7x, &
      'data point',7x,'approximation',3x,'deviation x 10**',i1/ &
      (i4,f16.8,f16.8,f17.6))
                                        return
      end
