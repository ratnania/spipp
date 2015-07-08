       subroutine eqblok ( t, n, kpm,  work1, work2, &
                       bloks, lenblk, integs, nbloks,  b )
!  from  * a practical guide to splines *  by c. de boor
!alls putit(difequ,bsplvd(bsplvb))
!  to be called in  c o l l o c
!
!******  i n p u t  ******
!  t   the knot sequence, of length n+kpm
!  n   the dimension of the approximating spline space, i.e., the order
!      of the linear system to be constructed.
!  kpm = k+m, the order of the approximating spline
!  lenblk   the maximum length of the array  bloks  as allowed by the
!           dimension statement in  colloc .
!
!******  w o r k   a r e a s  ******
!  work1    used in  putit, of size (kpm,kpm)
!  work2    used in  putit, of size (kpm,m+1)
!
!******  o u t p u t  ******
!  bloks    the coefficient matrix of the linear system, stored in al-
!           most block diagonal form, of size
!              kpm*sum(integs(1,i) , i=1,...,nbloks)
!  integs   an integer array, of size (3,nbloks), describing the block
!           structure.
!           integs(1,i)  =  number of rows in block  i
!           integs(2,i)  =  number of columns in block  i
!           integs(3,i)  =  number of elimination steps which can be
!                       carried out in block  i  before pivoting might
!                       bring in an equation from the next block.
!  nbloks   number of blocks, equals number of polynomial pieces
!  b   the right side of the linear system, stored corresponding to the
!      almost block diagonal form, of size sum(integs(1,i) , i=1,...,
!      nbloks).
!
!******  m e t h o d  ******
!  each breakpoint interval gives rise to a block in the linear system.
!  this block is determined by the  k  colloc.equations in the interval
!  with the side conditions (if any) in the interval interspersed ap-
!  propriately, and involves the  kpm  b-splines having the interval in
!  their support. correspondingly, such a block has  nrow = k + isidel
!  rows, with  isidel = number of side conditions in this and the prev-
!  ious intervals, and  ncol = kpm  columns.
!     further, because the interior knots have multiplicity  k, we can
!  carry out (in slvblk)  k  elimination steps in a block before pivot-
!  ing might involve an equation from the next block. in the last block,
!  of course, all kpm elimination steps will be carried out (in slvblk).
!
!  see the detailed comments in the solveblok package for further in-
!  formation about the almost block diagonal form used here.
      integer integs(3,1),kpm,lenblk,n,nbloks,   i,index,indexb,iside &
                                            ,isidel,itermx,k,left,m,nrow
      real b(1),bloks(1),t(1),work1(1),work2(1),   rho,xside
      common /side/ m, iside, xside(10)
      common /other/ itermx,k,rho(19)
      index = 1
      indexb = 1
      i = 0
      iside = 1
      do 20 left=kpm,n,k
         i = i+1
!        determine integs(.,i)
         integs(2,i) = kpm
         if (left .lt. n)               go to 14
         integs(3,i) = kpm
         isidel = m
                                        go to 16
   14    integs(3,i) = k
!        at this point,  iside-1  gives the number of side conditions
!        incorporated so far. adding to this the side conditions in the
!        current interval gives the number  isidel .
         isidel = iside-1
   15    if (isidel .eq. m)             go to 16
         if (xside(isidel+1) .ge. t(left+1)) &
                                        go to 16
         isidel = isidel+1
                                        go to 15
   16    nrow = k + isidel
         integs(1,i) = nrow
!        the detailed equations for this block are generated and put
!        together in  p u t i t .
         if (lenblk .lt. index+nrow*kpm-1)go to 999
         call putit(t,kpm,left,work1,work2,bloks(index),nrow,b(indexb))
         index = index + nrow*kpm
   20    indexb = indexb + nrow
      nbloks = i
                                        return
  999 print 699,lenblk
  699 format(11h **********/23h the assigned dimension,i5 &
              ,38h for  bloks  in  colloc  is too small.)
                                        stop
      end
