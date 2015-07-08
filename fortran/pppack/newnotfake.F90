      subroutine newnot ( break, coef, l, k, brknew, lnew, coefg )
!  from  * a practical guide to splines *  by c. de boor
!  this is a fake version of  n e w n o t  , of use in example 3 of
!  chapter xiv .
!  returns  lnew+1  knots in  brknew  which are equidistributed on (a,b)
!  = (break(1),break(l+1)) .
!
      integer k,l,lnew,   i
      real break(1),brknew(1),coef(k,l),coefg(2,l),   step
!******  i n p u t ******
!  break, coef, l, k.....contains the pp-representation of a certain
!        function  f  of order  k . specifically,
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
      brknew(1) = break(1)
      brknew(lnew+1) = break(l+1)
      step = (break(l+1) - break(1))/float(lnew)
      do 93 i=2,lnew
   93    brknew(i) = break(1) + float(i-1)*step
                                        return
      end
