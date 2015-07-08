      subroutine difequ ( mode, xx, v )
!  from  * a practical guide to splines *  by c. de boor
!alls  ppvalu(interv)
!  to be called by   c o l l o c ,  p u t i t
!  information about the differential equation is dispensed from here
!
!******  i n p u t  ******
!  mode  an integer indicating the task to be performed.
!        = 1   initialization
!        = 2   evaluate  de  at  xx
!        = 3   specify the next side condition
!        = 4   analyze the approximation
!  xx a point at which information is wanted
!
!******  o u t p u t  ******
!  v  depends on the  mode  . see comments below
!
      integer mode,   i,iside,itermx,k,kmax,kpm,l,m,ncoef,npiece
      parameter (npiece=100, kmax=20, ncoef=npiece*kmax)
      real v(kmax),xx,   break,coef,eps,ep1,ep2,error,factor,rho,solutn &
                      ,s2ovep,un,x,xside
      common /approx/ break(npiece),coef(ncoef),l,kpm
!     common /approx/ break(100),coef(2000),l,kpm
      common /side/ m,iside,xside(10)
      common /other/ itermx,k,rho(kmax-1)
      save eps,factor,s2ovep
!
!  this sample of  difequ  is for the example in chapter xv. it is a
!  nonlinear second order two point boundary value problem.
!
                                        go to (10,20,30,40),mode
!  initialize everything
!  i.e. set the order  m  of the dif.equ., the nondecreasing sequence
!  xside(i),i=1,...,m, of points at which side cond.s are given and
!  anything else necessary.
   10 m = 2
      xside(1) = 0.
      xside(2) = 1.
!  *** print out heading
      print 499
  499 format(' carrier,s nonlinear perturb. problem')
      eps = .5e-2
      print 610, eps
  610 format(' eps ',e20.10)
!  *** set constants used in formula for solution below.
      factor = (sqrt(2.) + sqrt(3.))**2
      s2ovep = sqrt(2./eps)
!  *** initial guess for newton iteration. un(x) = x*x - 1.
      l = 1
      break(1) = 0.
      do 16 i=1,kpm
   16    coef(i) = 0.
      coef(1) = -1.
      coef(3) = 2.
      itermx = 10
                                        return
!
!  provide value of left side coeff.s and right side at  xx .
!  specifically, at  xx  the dif.equ. reads
!        v(m+1)d**m + v(m)d**(m-1) + ... + v(1)d**0  =  v(m+2)
!  in terms of the quantities v(i),i=1,...,m+2, to be computed here.
   20 continue
      v(3) = eps
      v(2) = 0.
      un = ppvalu(break,coef,l,kpm,xx,0)
      v(1) = 2.*un
      v(4) = un**2 + 1.
                                        return
!
!  provide the  m  side conditions. these conditions are of the form
!        v(m+1)d**m + v(m)d**(m-1) + ... + v(1)d**0  =  v(m+2)
!  in terms of the quantities v(i),i=1,...,m+2, to be specified here.
!  note that v(m+1) = 0  for customary side conditions.
   30 v(m+1) = 0.
                                        go to (31,32,39),iside
   31 v(2) = 1.
      v(1) = 0.
      v(4) = 0.
                                        go to 38
   32 v(2) = 0.
      v(1) = 1.
      v(4) = 0.
   38 iside = iside + 1
   39                                   return
!
!  calculate the error near the boundary layer at  1.
   40 continue
      print 640
  640 format(' x, g(x)  and  g(x)-f(x)  at selected points')
      x = .75
      do 41 i=1,9
         ep1 = exp(s2ovep*(1.-x))*factor
         ep2 = exp(s2ovep*(1.+x))*factor
         solutn = 12./(1.+ep1)**2*ep1 + 12./(1.+ep2)**2*ep2 - 1.
         error = solutn - ppvalu(break,coef,l,kpm,x,0)
         print 641,x,solutn,error
  641    format(3e20.10)
   41    x = x + .03125
                                        return
      end
