      subroutine bchfac ( w, nbands, nrow, diag )
!  from  * a practical guide to splines *  by c. de boor
!onstructs cholesky factorization
!                     c  =  l * d * l-transpose
!  with l unit lower triangular and d diagonal, for given matrix c of
!  order  n r o w , in case  c  is (symmetric) positive semidefinite
!  and  b a n d e d , having  n b a n d s  diagonals at and below the
!  main diagonal.
!
!******  i n p u t  ******
!  nrow.....is the order of the matrix  c .
!  nbands.....indicates its bandwidth, i.e.,
!          c(i,j) = 0 for i-j .ge. nbands .
!  w.....workarray of size (nbands,nrow)  containing the  nbands  diago-
!        nals in its rows, with the main diagonal in row  1 . precisely,
!        w(i,j)  contains  c(i+j-1,j), i=1,...,nbands, j=1,...,nrow.
!          for example, the interesting entries of a seven diagonal sym-
!        metric matrix  c  of order  9  would be stored in  w  as
!
!
!
!
!
!
!        all other entries of  w  not identified in this way with an en-
!        try of  c  are never referenced .
!  diag.....is a work array of length  nrow .
!
!******  o u t p u t  ******
!  w.....contains the cholesky factorization  c = l*d*l-transp, with
!        w(1,i) containing  1/d(i,i)
!        and  w(i,j)  containing  l(i-1+j,j), i=2,...,nbands.
!
!******  m e t h o d  ******
!   gauss elimination, adapted to the symmetry and bandedness of  c , is
!   used .
!     near zero pivots are handled in a special way. the diagonal ele-
!  ment c(n,n) = w(1,n) is saved initially in  diag(n), all n. at the n-
!  th elimination step, the current pivot element, viz.  w(1,n), is com-
!  pared with its original value, diag(n). if, as the result of prior
!  elimination steps, this element has been reduced by about a word
!  length, (i.e., if w(1,n)+diag(n) .le. diag(n)), then the pivot is de-
!  clared to be zero, and the entire n-th row is declared to be linearly
!  dependent on the preceding rows. this has the effect of producing
!   x(n) = 0  when solving  c*x = b  for  x, regardless of  b. justific-
!  ation for this is as follows. in contemplated applications of this
!  program, the given equations are the normal equations for some least-
!  squares approximation problem, diag(n) = c(n,n) gives the norm-square
!  of the n-th basis function, and, at this point,  w(1,n)  contains the
!  norm-square of the error in the least-squares approximation to the n-
!  th basis function by linear combinations of the first n-1 . having
!  w(1,n)+diag(n) .le. diag(n) signifies that the n-th function is lin-
!  early dependent to machine accuracy on the first n-1 functions, there
!  fore can safely be left out from the basis of approximating functions
!     the solution of a linear system
!                       c*x = b
!   is effected by the succession of the following  t w o  calls:
!     call bchfac ( w, nbands, nrow, diag )       , to get factorization
!     call bchslv ( w, nbands, nrow, b )          , to solve for x.
!
      integer nbands,nrow,   i,imax,j,jmax,n
      real w(nbands,nrow),diag(nrow),   ratio
      if (nrow .gt. 1)                  go to 9
      if (w(1,1) .gt. 0.) w(1,1) = 1./w(1,1)
                                        return
!                                        store diagonal of  c  in  diag.
    9 do 10 n=1,nrow
   10    diag(n) = w(1,n)
!                                                        factorization .
      do 20 n=1,nrow
         if (w(1,n)+diag(n) .gt. diag(n)) go to 15
         do 14 j=1,nbands
   14       w(j,n) = 0.
                                        go to 20
   15    w(1,n) = 1./w(1,n)
         imax = min0(nbands-1,nrow - n)
         if (imax .lt. 1)               go to 20
         jmax = imax
         do 18 i=1,imax
            ratio = w(i+1,n)*w(1,n)
            do 17 j=1,jmax
   17          w(j,n+i) = w(j,n+i) - w(j+i,n)*ratio
            jmax = jmax - 1
   18       w(i+1,n) = ratio
   20    continue
                                        return
      end
