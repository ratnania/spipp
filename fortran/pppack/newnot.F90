      subroutine newnot ( break, coef, l, k, brknew, lnew, coefg )
!  from  * a practical guide to splines *  by c. de boor
!  returns  lnew+1  knots in  brknew  which are equidistributed on (a,b)
!  = (break(1),break(l+1)) wrto a certain monotone fctn  g  related to
!  the k-th root of the k-th derivative of the pp function  f  whose pp-
!  representation is contained in  break, coef, l, k .
!
!******  i n p u t ******
!  break, coef, l, k.....contains the pp-representation of a certain
!        function  f  of order  k . Specifically,
!        d**(k-1)f(x) = coef(k,i)  for  break(i).le. x .lt.break(i+1)
!  lnew.....number of intervals into which the interval (a,b) is to be
!        sectioned by the new breakpoint sequence  brknew .
!
!******  o u t p u t  ******
!  brknew.....array of length  lnew+1  containing the new breakpoint se-
!        quence
!  coefg.....the coefficient part of the pp-repr.  break, coefg, l, 2
!        for the monotone p.linear function  g  wrto which  brknew  will
!        be equidistributed.
!
!******  optional  p r i n t e d  o u t p u t  ******
!  coefg.....the pp coeffs of  g  are printed out if  iprint  is set
!        .gt. 0  in data statement below.
!
!******  m e t h o d  ******
!     The k-th derivative of the given pp function  f  does not exist
!  (except perhaps as a linear combination of delta functions). Never-
!  theless, we construct a p.constant function  h  with breakpoint se-
!  quence  break  which is approximately proportional to abs(d**k(f)).
!  Specifically, on  (break(i), break(i+1)),
!
!     abs(jump at break(i) of pc)    abs(jump at break(i+1) of pc)
! h = --------------------------  +  ----------------------------
!       break(i+1) - break(i-1)         break(i+2) - break(i)
!
!  with  pc  the p.constant (k-1)st derivative of  f .
!      Then, the p.linear function  g  is constructed as
!
!    g(x)  =  integral of  h(y)**(1/k)  for  y  from  a  to  x
!
!  and its pp coeffs. stored in  coefg .
!     then  brknew  is determined by
!
!        brknew(i)  =  a + g**(-1)((i-1)*step) , i=1,...,lnew+1
!
!  where  step = g(b)/lnew  and  (a,b) = (break(1),break(l+1)) .
!     In the event that  pc = d**(k-1)(f) is constant in  (a,b)  and
!  therefore  h = 0 identically,  brknew  is chosen uniformly spaced.
!
      integer k,l,lnew,   i,iprint,j
      real break(1),brknew(1),coef(k,l),coefg(2,l),   dif,difprv,oneovk &
                                                     ,step,stepi
!     dimension break(l+1), brknew(lnew+1)
!urrent fortran standard makes it impossible to specify the dimension
!  of  break   and  brknew  without the introduction of additional
!  otherwise superfluous arguments.
      data iprint /0/
!
      brknew(1) = break(1)
      brknew(lnew+1) = break(l+1)
!                               if  g  is constant,  brknew  is uniform.
      if (l .le. 1)                     go to 90
!                        construct the continuous p.linear function  g .
      oneovk = 1./float(k)
      coefg(1,1) = 0.
      difprv = abs(coef(k,2) - coef(k,1))/(break(3)-break(1))
      do 10 i=2,l
         dif = abs(coef(k,i) - coef(k,i-1))/(break(i+1) - break(i-1))
         coefg(2,i-1) = (dif + difprv)**oneovk
         coefg(1,i) = coefg(1,i-1)+coefg(2,i-1)*(break(i)-break(i-1))
   10    difprv = dif
      coefg(2,l) = (2.*difprv)**oneovk
!                                                     step  =  g(b)/lnew
      step = (coefg(1,l)+coefg(2,l)*(break(l+1)-break(l)))/float(lnew)
!
      if (iprint .gt. 0) print 600, step,(i,coefg(1,i),coefg(2,i),i=1,l)
  600 format(7h step =,e16.7/(i5,2e16.5))
!                              if  g  is constant,  brknew  is uniform .
      if (step .le. 0.)                 go to 90
!
!      For i=2,...,lnew, construct  brknew(i) = a + g**(-1)(stepi),
!      with  stepi = (i-1)*step .  this requires inversion of the p.lin-
!      ear function  g .  For this,  j  is found so that
!         g(break(j)) .le. stepi .le. g(break(j+1))
!      and then
!         brknew(i)  =  break(j) + (stepi-g(break(j)))/dg(break(j)) .
!      The midpoint is chosen if  dg(break(j)) = 0 .
      j = 1
      do 30 i=2,lnew
         stepi = float(i-1)*step
   21       if (j .eq. l)               go to 27
            if (stepi .le. coefg(1,j+1))go to 27
            j = j + 1
                                        go to 21
   27    if (coefg(2,j) .eq. 0.)        go to 29
            brknew(i) = break(j) + (stepi - coefg(1,j))/coefg(2,j)
                                        go to 30
   29       brknew(i) = (break(j) + break(j+1))/2.
   30    continue
                                        return
!
!                              if  g  is constant,  brknew  is uniform .
   90 step = (break(l+1) - break(1))/float(lnew)
      do 93 i=2,lnew
   93    brknew(i) = break(1) + float(i-1)*step
                                        return
      end
