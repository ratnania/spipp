      subroutine slvblk ( bloks, integs, nbloks, b, ipivot, x, iflag )
!    this program solves  the  linear system  a*x = b  where a is an
!  almost block diagonal matrix.  such almost block diagonal matrices
!  arise naturally in piecewise polynomial interpolation or approx-
!  imation and in finite element methods for two-point boundary value
!  problems.  the plu factorization method is implemented here to take
!  advantage of the special structure of such systems for savings in
!  computing time and storage requirements.
!
!                  parameters
!  bloks   a one-dimenional array, of length
!                   sum( integs(1,i)*integs(2,i) ; i = 1,nbloks )
!          on input, contains the blocks of the almost block diagonal
!          matrix  a  .  the array integs (see below and the example)
!          describes the block structure.
!          on output, contains correspondingly the plu factorization
!          of  a  (if iflag .ne. 0).  certain of the entries into bloks
!          are arbitrary (where the blocks overlap).
!  integs  integer array description of the block structure of  a .
!            integs(1,i) = no. of rows of block i        =  nrow
!            integs(2,i) = no. of colums of block i      =  ncol
!            integs(3,i) = no. of elim. steps in block i =  last
!                          i  = 1,2,...,nbloks
!          the linear system is of order
!                n  =  sum ( integs(3,i) , i=1,...,nbloks ),
!          but the total number of rows in the blocks is
!              nbrows = sum( integs(1,i) ; i = 1,...,nbloks)
!  nbloks  number of blocks
!  b       right side of the linear system, array of length nbrows.
!          certain of the entries are arbitrary, corresponding to
!          rows of the blocks which overlap (see block structure and
!          the example below).
!  ipivot  on output, integer array containing the pivoting sequence
!          used. length is nbrows
!  x       on output, contains the computed solution (if iflag .ne. 0)
!          length is n.
!  iflag   on output, integer
!            = (-1)**(no. of interchanges during factorization)
!                   if  a  is invertible
!            = 0    if  a  is singular
!
!                   auxiliary programs
!  fcblok (bloks,integs,nbloks,ipivot,scrtch,iflag)  factors the matrix
!           a , and is used for this purpose in slvblk. its arguments
!          are as in slvblk, except for
!              scrtch = a work array of length max(integs(1,i)).
!
!  sbblok (bloks,integs,nbloks,ipivot,b,x)  solves the system a*x = b
!          once  a  is factored. this is done automatically by slvblk
!          for one right side b, but subsequent solutions may be
!          obtained for additional b-vectors. the arguments are all
!          as in slvblk.
!
!  dtblok (bloks,integs,nbloks,ipivot,iflag,detsgn,detlog) computes the
!          determinant of  a  once slvblk or fcblok has done the fact-
!          orization.the first five arguments are as in slvblk.
!              detsgn  = sign of the determinant
!              detlog  = natural log of the determinant
!
!             ------ block structure of  a  ------
!  the nbloks blocks are stored consecutively in the array  bloks .
!  the first block has its (1,1)-entry at bloks(1), and, if the i-th
!  block has its (1,1)-entry at bloks(index(i)), then
!         index(i+1) = index(i)  +  nrow(i)*ncol(i) .
!    the blocks are pieced together to give the interesting part of  a
!  as follows.  for i = 1,2,...,nbloks-1, the (1,1)-entry of the next
!  block (the (i+1)st block ) corresponds to the (last+1,last+1)-entry
!  of the current i-th block.  recall last = integs(3,i) and note that
!  this means that
!      a. every block starts on the diagonal of  a .
!      b. the blocks overlap (usually). the rows of the (i+1)st block
!         which are overlapped by the i-th block may be arbitrarily de-
!         fined initially. they are overwritten during elimination.
!    the right side for the equations in the i-th block are stored cor-
!  respondingly as the last entries of a piece of  b  of length  nrow
!  (= integs(1,i)) and following immediately in  b  the corresponding
!  piece for the right side of the preceding block, with the right side
!  for the first block starting at  b(1) . in this, the right side for
!  an equation need only be specified once on input, in the first block
!  in which the equation appears.
!
!             ------ example and test driver ------
!    the test driver for this package contains an example, a linear
!  system of order 11, whose nonzero entries are indicated in the fol-
!  lowing schema by their row and column index modulo 10. next to it
!  are the contents of the  integs  arrray when the matrix is taken to
!  be almost block diagonal with  nbloks = 5, and below it are the five
!  blocks.
!
!                      nrow1 = 3, ncol1 = 4
!           11 12 13 14
!           21 22 23 24   nrow2 = 3, ncol2 = 3
!           31 32 33 34
!  last1 = 2      43 44 45
!                 53 54 55            nrow3 = 3, ncol3 = 4
!        last2 = 3         66 67 68 69   nrow4 = 3, ncol4 = 4
!                          76 77 78 79      nrow5 = 4, ncol5 = 4
!                          86 87 88 89
!                 last3 = 1   97 98 99 90
!                    last4 = 1   08 09 00 01
!                                18 19 10 11
!                       last5 = 4
!
!         actual input to bloks shown by rows of blocks of  a .
!      (the ** items are arbitrary, this storage is used by slvblk)
!
!  11 12 13 14  / ** ** **  / 66 67 68 69  / ** ** ** **  / ** ** ** **
!  21 22 23 24 /  43 44 45 /  76 77 78 79 /  ** ** ** ** /  ** ** ** **
!  31 32 33 34/   53 54 55/   86 87 88 89/   97 98 99 90/   08 09 00 01
!                                                           18 19 10 11
!
!  index = 1      index = 13  index = 22     index = 34     index = 46
!
!         actual right side values with ** for arbitrary values
!  b1 b2 b3 ** b4 b5 b6 b7 b8 ** ** b9 ** ** b10 b11
!
!  (it would have been more efficient to combine block 3 with block 4)
!
      integer integs(3,nbloks),ipivot(1),iflag
      real bloks(1),b(1),x(1)
!     in the call to fcblok,  x  is used for temporary storage.
      call fcblok(bloks,integs,nbloks,ipivot,x,iflag)
      if (iflag .eq. 0)                 return
      call sbblok(bloks,integs,nbloks,ipivot,b,x)
                                        return
      end
      subroutine fcblok ( bloks, integs, nbloks, ipivot, scrtch, iflag )
!alls subroutines  f a c t r b  and  s h i f t b .
!
!   f c b l o k  supervises the plu factorization with pivoting of
!  scaled rows of the almost block diagonal matrix stored in the arrays
!   b l o k s  and  i n t e g s .
!
!   factrb = subprogram which carries out steps 1,...,last of gauss
!            elimination (with pivoting) for an individual block.
!   shiftb = subprogram which shifts the remaining rows to the top of
!            the next block
!
! parameters
!    bloks   an array that initially contains the almost block diagonal
!            matrix  a  to be factored, and on return contains the com-
!            puted factorization of  a .
!    integs  an integer array describing the block structure of  a .
!    nbloks  the number of blocks in  a .
!    ipivot  an integer array of dimension  sum (integs(1,n) ; n=1,
!            ...,nbloks) which, on return, contains the pivoting stra-
!            tegy used.
!    scrtch  work area required, of length  max (integs(1,n) ; n=1,
!            ...,nbloks).
!    iflag   output parameter;
!            = 0  in case matrix was found to be singular.
!            otherwise,
!            = (-1)**(number of row interchanges during factorization)
!
      integer integs(3,nbloks),ipivot(1),iflag, i,index,indexb,indexn, &
              last,ncol,nrow
      real bloks(1),scrtch(1)
      iflag = 1
      indexb = 1
      indexn = 1
      i = 1
!                        loop over the blocks.  i  is loop index
   10    index = indexn
         nrow = integs(1,i)
         ncol = integs(2,i)
         last = integs(3,i)
!        carry out elimination on the i-th block until next block
!        enters, i.e., for columns 1,...,last  of i-th block.
         call factrb(bloks(index),ipivot(indexb),scrtch,nrow,ncol,last, &
                  iflag)
!         check for having reached a singular block or the last block
         if (iflag .eq. 0 .or. i .eq. nbloks) &
                                        return
         i = i+1
         indexn = nrow*ncol + index
!              put the rest of the i-th block onto the next block
         call shiftb(bloks(index),ipivot(indexb),nrow,ncol,last, &
                  bloks(indexn),integs(1,i),integs(2,i))
         indexb = indexb + nrow
                                        go to 10
      end
      subroutine factrb ( w, ipivot, d, nrow, ncol, last, iflag )
!  adapted from p.132 of 'element.numer.analysis' by conte-de boor
!
!  constructs a partial plu factorization, corresponding to steps 1,...,
!   l a s t   in gauss elimination, for the matrix  w  of order
!   ( n r o w ,  n c o l ), using pivoting of scaled rows.
!
!  parameters
!    w       contains the (nrow,ncol) matrix to be partially factored
!            on input, and the partial factorization on output.
!    ipivot  an integer array of length nrow containing a record of the
!            pivoting strategy used; row ipivot(i) is used during the
!            i-th elimination step, i=1,...,last.
!    d       a work array of length nrow used to store row sizes
!            temporarily.
!    nrow    number of rows of w.
!    ncol    number of columns of w.
!    last    number of elimination steps to be carried out.
!    iflag   on output, equals iflag on input times (-1)**(number of
!            row interchanges during the factorization process), in
!            case no zero pivot was encountered.
!            otherwise, iflag = 0 on output.
!
      integer ipivot(nrow),ncol,last,iflag, i,ipivi,ipivk,j,k,kp1
      real w(nrow,ncol),d(nrow), awikdi,colmax,ratio,rowmax
!  initialize ipivot, d
      do 10 i=1,nrow
         ipivot(i) = i
         rowmax = 0.
         do 9 j=1,ncol
    9       rowmax = amax1(rowmax, abs(w(i,j)))
         if (rowmax .eq. 0.)            go to 999
   10    d(i) = rowmax
! gauss elimination with pivoting of scaled rows, loop over k=1,.,last
      k = 1
!        as pivot row for k-th step, pick among the rows not yet used,
!        i.e., from rows ipivot(k),...,ipivot(nrow), the one whose k-th
!        entry (compared to the row size) is largest. then, if this row
!        does not turn out to be row ipivot(k), redefine ipivot(k) ap-
!        propriately and record this interchange by changing the sign
!        of  i f l a g .
   11    ipivk = ipivot(k)
         if (k .eq. nrow)               go to 21
         j = k
         kp1 = k+1
         colmax = abs(w(ipivk,k))/d(ipivk)
!              find the (relatively) largest pivot
         do 15 i=kp1,nrow
            ipivi = ipivot(i)
            awikdi = abs(w(ipivi,k))/d(ipivi)
            if (awikdi .le. colmax)     go to 15
               colmax = awikdi
               j = i
   15       continue
         if (j .eq. k)                  go to 16
         ipivk = ipivot(j)
         ipivot(j) = ipivot(k)
         ipivot(k) = ipivk
         iflag = -iflag
   16    continue
!        if pivot element is too small in absolute value, declare
!        matrix to be noninvertible and quit.
         if (abs(w(ipivk,k))+d(ipivk) .le. d(ipivk)) &
                                        go to 999
!        otherwise, subtract the appropriate multiple of the pivot
!        row from remaining rows, i.e., the rows ipivot(k+1),...,
!        ipivot(nrow), to make k-th entry zero. save the multiplier in
!        its place.
         do 20 i=kp1,nrow
            ipivi = ipivot(i)
            w(ipivi,k) = w(ipivi,k)/w(ipivk,k)
            ratio = -w(ipivi,k)
            do 20 j=kp1,ncol
   20          w(ipivi,j) = ratio*w(ipivk,j) + w(ipivi,j)
         k = kp1
!        check for having reached the next block.
         if (k .le. last)               go to 11
                                        return
!     if  last  .eq. nrow , check now that pivot element in last row
!     is nonzero.
   21 if( abs(w(ipivk,nrow))+d(ipivk) .gt. d(ipivk) ) &
                                        return
!                   singularity flag set
  999 iflag = 0
                                        return
      end
      subroutine shiftb ( ai, ipivot, nrowi, ncoli, last, &
                          ai1, nrowi1, ncoli1 )
!  shifts the rows in current block, ai, not used as pivot rows, if
!  any, i.e., rows ipivot(last+1),...,ipivot(nrowi), onto the first
!  mmax = nrow-last rows of the next block, ai1, with column last+j of
!  ai  going to column j , j=1,...,jmax=ncoli-last. the remaining col-
!  umns of these rows of ai1 are zeroed out.
!
!                             picture
!
!       original situation after         results in a new block i+1
!       last = 2 columns have been       created and ready to be
!       done in factrb (assuming no      factored by next factrb call.
!       interchanges of rows)
!                   1
!              x  x 1x  x  x           x  x  x  x  x
!                   1
!              0  x 1x  x  x           0  x  x  x  x
!  block i          1                       ---------------
!  nrowi = 4   0  0 1x  x  x           0  0 1x  x  x  0  01
!  ncoli = 5        1                       1             1
!  last = 2    0  0 1x  x  x           0  0 1x  x  x  0  01
!  -------------------------------          1             1   new
!                   1x  x  x  x  x          1x  x  x  x  x1  block
!                   1                       1             1   i+1
!  block i+1        1x  x  x  x  x          1x  x  x  x  x1
!  nrowi1= 5        1                       1             1
!  ncoli1= 5        1x  x  x  x  x          1x  x  x  x  x1
!  -------------------------------          1-------------1
!                   1
!
      integer ipivot(nrowi),last, ip,j,jmax,jmaxp1,m,mmax
      real ai(nrowi,ncoli),ai1(nrowi1,ncoli1)
      mmax = nrowi - last
      jmax = ncoli - last
      if (mmax .lt. 1 .or. jmax .lt. 1) return
!              put the remainder of block i into ai1
      do 10 m=1,mmax
         ip = ipivot(last+m)
         do 10 j=1,jmax
   10       ai1(m,j) = ai(ip,last+j)
      if (jmax .eq. ncoli1)             return
!              zero out the upper right corner of ai1
      jmaxp1 = jmax + 1
      do 20 j=jmaxp1,ncoli1
         do 20 m=1,mmax
   20       ai1(m,j) = 0.
                                        return
      end
      subroutine sbblok ( bloks, integs, nbloks, ipivot, b, x )
!alls subroutines  s u b f o r  and  s u b b a k .
!
!  supervises the solution (by forward and backward substitution) of
!  the linear system  a*x = b  for x, with the plu factorization of  a
!  already generated in  f c b l o k .  individual blocks of equations
!  are solved via  s u b f o r  and  s u b b a k .
!
! parameters
!    bloks, integs, nbloks, ipivot    are as on return from fcblok.
!    b       the right side, stored corresponding to the storage of
!            the equations. see comments in  s l v b l k  for details.
!    x       solution vector
!
      integer integs(3,nbloks),ipivot(1), i,index,indexb,indexx,j,last, &
              nbp1,ncol,nrow
      real bloks(1),b(1),x(1)
!
!      forward substitution pass
!
      index = 1
      indexb = 1
      indexx = 1
      do 20 i=1,nbloks
         nrow = integs(1,i)
         last = integs(3,i)
         call subfor(bloks(index),ipivot(indexb),nrow,last,b(indexb), &
                     x(indexx))
         index = nrow*integs(2,i) + index
         indexb = indexb + nrow
   20    indexx = indexx + last
!
!     back substitution pass
!
      nbp1 = nbloks + 1
      do 30 j=1,nbloks
         i = nbp1 - j
         nrow = integs(1,i)
         ncol = integs(2,i)
         last = integs(3,i)
         index = index - nrow*ncol
         indexb = indexb - nrow
         indexx = indexx - last
   30    call subbak(bloks(index),ipivot(indexb),nrow,ncol,last, &
                     x(indexx))
                                        return
      end
      subroutine subfor ( w, ipivot, nrow, last, b, x )
!  carries out the forward pass of substitution for the current block,
!  i.e., the action on the right side corresponding to the elimination
!  carried out in  f a c t r b  for this block.
!     at the end, x(j) contains the right side of the transformed
!  ipivot(j)-th equation in this block, j=1,...,nrow. then, since
!  for i=1,...,nrow-last, b(nrow+i) is going to be used as the right
!  side of equation  i  in the next block (shifted over there from
!  this block during factorization), it is set equal to x(last+i) here.
!
! parameters
!    w, ipivot, nrow, last  are as on return from factrb.
!    b(j)   is expected to contain, on input, the right side of j-th
!           equation for this block, j=1,...,nrow.
!    b(nrow+j)   contains, on output, the appropriately modified right
!           side for equation j in next block, j=1,...,nrow-last.
!    x(j)   contains, on output, the appropriately modified right
!           side of equation ipivot(j) in this block, j=1,...,last (and
!           even for j=last+1,...,nrow).
!
      integer ipivot(nrow), ip,jmax,k
!     dimension b(nrow + nrow-last)
      real w(nrow,last),b(1),x(nrow)
      ip = ipivot(1)
      x(1) = b(ip)
      if (nrow .eq. 1)                  go to 99
      do 15 k=2,nrow
         ip = ipivot(k)
         jmax = amin0(k-1,last)
         sum = 0.
         do 14 j=1,jmax
   14       sum = w(ip,j)*x(j) + sum
   15    x(k) = b(ip) - sum
!
!     transfer modified right sides of equations ipivot(last+1),...,
!     ipivot(nrow) to next block.
      nrowml = nrow - last
      if (nrowml .eq. 0)                go to 99
      lastp1 = last+1
      do 25 k=lastp1,nrow
   25    b(nrowml+k) = x(k)
   99                                   return
      end
      subroutine subbak ( w, ipivot, nrow, ncol, last, x )
!  carries out backsubstitution for current block.
!
! parameters
!    w, ipivot, nrow, ncol, last  are as on return from factrb.
!    x(1),...,x(ncol)  contains, on input, the right side for the
!            equations in this block after backsubstitution has been
!            carried up to but not including equation ipivot(last).
!            means that x(j) contains the right side of equation ipi-
!            vot(j) as modified during elimination, j=1,...,last, while
!            for j .gt. last, x(j) is already a component of the solut-
!            ion vector.
!    x(1),...,x(ncol) contains, on output, the components of the solut-
!            ion corresponding to the present block.
!
      integer ipivot(nrow),last,  ip,j,k,kp1
      real w(nrow,ncol),x(ncol), sum
      k = last
      ip = ipivot(k)
      sum = 0.
      if (k .eq. ncol)                  go to 4
      kp1 = k+1
    2    do 3 j=kp1,ncol
    3       sum = w(ip,j)*x(j) + sum
    4    x(k) = (x(k) - sum)/w(ip,k)
         if (k .eq. 1)                  return
         kp1 = k
         k = k-1
         ip = ipivot(k)
         sum = 0.
                                        go to 2
      end
      subroutine dtblok ( bloks, integs, nbloks, ipivot, iflag, &
                          detsgn, detlog )
!  computes the determinant of an almost block diagonal matrix whose
!  plu factorization has been obtained previously in fcblok.
!  *** the logarithm of the determinant is computed instead of the
!  determinant itself to avoid the danger of overflow or underflow
!  inherent in this calculation.
!
! parameters
!    bloks, integs, nbloks, ipivot, iflag  are as on return from fcblok.
!            in particular, iflag = (-1)**(number of interchanges dur-
!            ing factorization) if successful, otherwise iflag = 0.
!    detsgn  on output, contains the sign of the determinant.
!    detlog  on output, contains the natural logarithm of the determi-
!            nant if determinant is not zero. otherwise contains 0.
!
      integer integs(3,nbloks),ipivot(1),iflag, i,indexp,ip,k,last
      real bloks(1),detsgn,detlog
!
      detsgn = iflag
      detlog = 0.
      if (iflag .eq. 0)                 return
      index = 0
      indexp = 0
      do 2 i=1,nbloks
         nrow = integs(1,i)
         last = integs(3,i)
         do 1 k=1,last
            ip = index + nrow*(k-1) + ipivot(indexp+k)
            detlog = detlog + alog(abs(bloks(ip)))
    1       detsgn = detsgn*sign(1.,bloks(ip))
         index = nrow*integs(2,i) + index
    2    indexp = indexp + nrow
                                        return
      end
