      subroutine cwidth ( w,b,nequ,ncols,integs,nbloks, d, x,iflag )
!  this program is a variation of the theme in the algorithm bandet1
!  by martin and wilkinson (numer.math. 9(1976)279-307). it solves
!  the linear system
!                           a*x  =  b
!  of  nequ  equations in case  a  is almost block diagonal with all
!  blocks having  ncols  columns using no more storage than it takes to
!  store the interesting part of  a . such systems occur in the determ-
!  ination of the b-spline coefficients of a spline approximation.
!
!                           parameters
!  w     on input, a two-dimensional array of size (nequ,ncols) contain-
!        ing the interesting part of the almost block diagonal coeffici-
!        ent matrix  a (see description and example below). the array
!        integs  describes the storage scheme.
!        on output, w  contains the upper triangular factor  u  of the
!        lu factorization of a possibly permuted version of  a . in par-
!        ticular, the determinant of  a  could now be found as
!            iflag*w(1,1)*w(2,1)* ... * w(nequ,1)  .
!  b     on input, the right side of the linear system, of length  nequ.
!        the contents of  b  are changed during execution.
!  nequ  number of equations in system
!  ncols block width, i.e., number of columns in each block.
!  integs  integer array, of size (2,nequ), describing the block struct-
!        ure of  a .
!           integs(1,i) = no. of rows in block i               =  nrow
!           integs(2,i) = no. of elimination steps in block i
!                       = overhang over next block             =  last
!  nbloks  number of blocks
!  d     work array, to contain row sizes . if storage is scarce, the
!        array  x  could be used in the calling sequence for  d .
!  x     on output, contains computed solution (if iflag .ne. 0), of
!        length  nequ .
!  iflag  on output, integer
!         = (-1)**(no.of interchanges during elimination)
!                if  a  is invertible
!         =  0   if  a  is singular
!
!        ------  block structure of  a  ------
!     the interesting part of  a  is taken to consist of  nbloks  con-
!  secutive blocks, with the i-th block made up of  nrowi = integs(1,i)
!  consecutive rows and  ncols  consecutive columns of  a , and with
!  the first  lasti = integs(2,i) columns to the left of the next block.
!  these blocks are stored consecutively in the workarray  w .
!     for example, here is an 11th order matrix and its arrangement in
!  the workarray  w . (the interesting entries of  a  are indicated by
!  their row and column index modulo 10.)
!
!                  ---   a   ---                          ---   w   ---
!
!                     nrow1=3
!          11 12 13 14                                     11 12 13 14
!          21 22 23 24                                     21 22 23 24
!          31 32 33 34      nrow2=2                        31 32 33 34
!   last1=2      43 44 45 46                               43 44 45 46
!                53 54 55 56         nrow3=3               53 54 55 56
!         last2=3         66 67 68 69                      66 67 68 69
!                         76 77 78 79                      76 77 78 79
!                         86 87 88 89   nrow4=1            86 87 88 89
!                  last3=1   97 98 99 90   nrow5=2         97 98 99 90
!                     last4=1   08 09 00 01                08 09 00 01
!                               18 19 10 11                18 19 10 11
!                        last5=4
!
!  for this interpretation of  a  as an almost block diagonal matrix,
!  we have  nbloks = 5 , and the integs array is
!
!                        i=  1   2   3   4   5
!                  k=
!  integs(k,i) =      1      3   2   3   1   2
!                     2      2   3   1   1   4
!
!       --------  method  --------
!     gauss elimination with scaled partial pivoting is used, but mult-
!  ipliers are  n o t  s a v e d  in order to save storage. rather, the
!  right side is operated on during elimination.
!     the two parameters
!                  i p v t e q   and  l a s t e q
!  are used to keep track of the action.  ipvteq is the index of the
!  variable to be eliminated next, from equations  ipvteq+1,...,lasteq,
!  using equation  ipvteq (possibly after an interchange) as the pivot
!  equation. the entries in the pivot column are  a l w a y s  in column
!  1 of  w . this is accomplished by putting the entries in rows
!  ipvteq+1,...,lasteq  revised by the elimination of the  ipvteq-th
!  variable one to the left in  w . in this way, the columns of the
!  equations in a given block (as stored in  w ) will be aligned with
!  those of the next block at the moment when these next equations be-
!  come involved in the elimination process.
!     thus, for the above example, the first elimination steps proceed
!  as follows.
!
!  *11 12 13 14    11 12 13 14    11 12 13 14    11 12 13 14
!  *21 22 23 24   *22 23 24       22 23 24       22 23 24
!  *31 32 33 34   *32 33 34      *33 34          33 34
!   43 44 45 46    43 44 45 46   *43 44 45 46   *44 45 46        etc.
!   53 54 55 56    53 54 55 56   *53 54 55 56   *54 55 56
!   66 67 68 69    66 67 68 69    66 67 68 69    66 67 68 69
!        .              .              .              .
!
!     in all other respects, the procedure is standard, including the
!  scaled partial pivoting.
!
      integer iflag,integs(2,nbloks),ncols,nequ,   i,ii,icount,ipvteq &
                           ,istar,j,lastcl,lasteq,lasti,nexteq,nrowad
      real b(nequ),d(nequ),w(nequ,ncols),x(nequ),   awi1od,colmax &
                                                   ,rowmax,temp
      iflag = 1
      ipvteq = 0
      lasteq = 0
!                                 the i-loop runs over the blocks
      do 50 i=1,nbloks
!
!        the equations for the current block are added to those current-
!        ly involved in the elimination process, by increasing  lasteq
!        by  integs(1,i) after the rowsize of these equations has been
!        recorded in the array  d .
!
         nrowad = integs(1,i)
         do 10 icount=1,nrowad
            nexteq = lasteq + icount
            rowmax = 0.
            do 5 j=1,ncols
    5          rowmax = amax1(rowmax,abs(w(nexteq,j)))
            if (rowmax .eq. 0.)         go to 999
   10       d(nexteq) = rowmax
         lasteq = lasteq + nrowad
!
!        there will be  lasti = integs(2,i)  elimination steps before
!        the equations in the next block become involved. further,
!        l a s t c l  records the number of columns involved in the cur-
!        rent elimination step. it starts equal to  ncols  when a block
!        first becomes involved and then drops by one after each elim-
!        ination step.
!
         lastcl = ncols
         lasti = integs(2,i)
         do 30 icount=1,lasti
            ipvteq = ipvteq + 1
            if (ipvteq .lt. lasteq)     go to 11
            if ( abs(w(ipvteq,1))+d(ipvteq) .gt. d(ipvteq) ) &
                                        go to 50
                                        go to 999
!
!        determine the smallest  i s t a r  in  (ipvteq,lasteq)  for
!        which  abs(w(istar,1))/d(istar)  is as large as possible, and
!        interchange equations  ipvteq  and  istar  in case  ipvteq
!        .lt. istar .
!
   11       colmax = abs(w(ipvteq,1))/d(ipvteq)
            istar = ipvteq
            ipvtp1 = ipvteq + 1
            do 13 ii=ipvtp1,lasteq
               awi1od = abs(w(ii,1))/d(ii)
               if (awi1od .le. colmax)  go to 13
               colmax = awi1od
               istar = ii
   13          continue
            if ( abs(w(istar,1))+d(istar) .eq. d(istar) ) &
                                        go to 999
            if (istar .eq. ipvteq)      go to 16
            iflag = -iflag
            temp = d(istar)
            d(istar) = d(ipvteq)
            d(ipvteq) = temp
            temp = b(istar)
            b(istar) = b(ipvteq)
            b(ipvteq) = temp
            do 14 j=1,lastcl
               temp = w(istar,j)
               w(istar,j) = w(ipvteq,j)
   14          w(ipvteq,j) = temp
!
!        subtract the appropriate multiple of equation  ipvteq  from
!        equations  ipvteq+1,...,lasteq to make the coefficient of the
!        ipvteq-th unknown (presently in column 1 of  w ) zero, but
!        store the new coefficients in  w  one to the left from the old.
!
   16       do 20 ii=ipvtp1,lasteq
               ratio = w(ii,1)/w(ipvteq,1)
               do 18 j=2,lastcl
   18             w(ii,j-1) = w(ii,j) - ratio*w(ipvteq,j)
               w(ii,lastcl) = 0.
   20          b(ii) = b(ii) - ratio*b(ipvteq)
   30       lastcl = lastcl - 1
   50    continue
!
!  at this point,  w  and  b  contain an upper triangular linear system
!  equivalent to the original one, with  w(i,j) containing entry
!  (i, i-1+j ) of the coefficient matrix. solve this system by backsub-
!  stitution, taking into account its block structure.
!
!                          i-loop over the blocks, in reverse order
      i = nbloks
   59    lasti = integs(2,i)
         jmax = ncols - lasti
         do 70 icount=1,lasti
            sum = 0.
            if (jmax .eq. 0)            go to 61
            do 60 j=1,jmax
   60          sum = sum + x(ipvteq+j)*w(ipvteq,j+1)
   61       x(ipvteq) = (b(ipvteq)-sum)/w(ipvteq,1)
            jmax = jmax + 1
   70       ipvteq = ipvteq - 1
         i = i - 1
         if (i .gt. 0)                  go to 59
                                        return
  999 iflag = 0
                                        return
      end
