      MODULE SPI_SPARSEKIT
      CONTAINS

! ...................................................
      subroutine getbwd(n,a,ja,ia,ml,mu)
!-----------------------------------------------------------------------
! gets the bandwidth of lower part and upper part of A.
! does not assume that A is sorted.
!-----------------------------------------------------------------------
! on entry:
!----------
! n      = integer = the row dimension of the matrix
! a, ja,
!    ia = matrix in compressed sparse row format.
!
! on return:
!-----------
! ml      = integer. The bandwidth of the strict lower part of A
! mu      = integer. The bandwidth of the strict upper part of A
!
! Notes:
! ===== ml and mu are allowed to be negative or return. This may be
!       useful since it will tell us whether a band is confined
!       in the strict  upper/lower triangular part.
!       indeed the definitions of ml and mu are
!
!       ml = max ( (i-j)  s.t. a(i,j) .ne. 0  )
!       mu = max ( (j-i)  s.t. a(i,j) .ne. 0  )
!----------------------------------------------------------------------c
! Y. Saad, Sep. 21 1989                                                c
!----------------------------------------------------------------------c
      real*8 a(*)
      integer ja(*),ia(n+1),ml,mu,ldist,i,k
      ml = - n
      mu = - n
      do 3 i=1,n
         do 31 k=ia(i),ia(i+1)-1
            ldist = i-ja(k)
            ml = max(ml,ldist)
            mu = max(mu,-ldist)
 31      continue
 3    continue
      return
      end subroutine getbwd
! ...................................................


! ...................................................
      subroutine csrbnd (n,a,ja,ia,job,abd,nabd,lowd,ml,mu,ierr)
      real*8 a(*),abd(nabd,n)
      integer ia(n+1),ja(*)
!-----------------------------------------------------------------------
!   Compressed Sparse Row  to  Banded (Linpack ) format.
!-----------------------------------------------------------------------
! this subroutine converts a general sparse matrix stored in
! compressed sparse row format into the banded format. for the
! banded format,the Linpack conventions are assumed (see below).
!-----------------------------------------------------------------------
! on entry:
!----------
! n      = integer,the actual row dimension of the matrix.
!
! a,
! ja,
! ia    = input matrix stored in compressed sparse row format.
!
! job      = integer. if job=1 then the values of the lower bandwith ml
!         and the upper bandwidth mu are determined internally.
!         otherwise it is assumed that the values of ml and mu
!         are the correct bandwidths on input. See ml and mu below.
!
! nabd  = integer. first dimension of array abd.
!
! lowd  = integer. this should be set to the row number in abd where
!         the lowest diagonal (leftmost) of A is located.
!         lowd should be  ( 1  .le.  lowd  .le. nabd).
!         if it is not known in advance what lowd should be
!         enter lowd = 0 and the default value lowd = ml+mu+1
!         will be chosen. Alternative: call routine getbwd from unary
!         first to detrermione ml and mu then define lowd accordingly.
!         (Note: the banded solvers in linpack use lowd=2*ml+mu+1. )
!
! ml      = integer. equal to the bandwidth of the strict lower part of A
! mu      = integer. equal to the bandwidth of the strict upper part of A
!         thus the total bandwidth of A is ml+mu+1.
!         if ml+mu+1 is found to be larger than lowd then an error
!         flag is raised (unless lowd = 0). see ierr.
!
! note:   ml and mu are assumed to have       the correct bandwidth values
!         as defined above if job is set to zero on entry.
!
! on return:
!-----------
!
! abd   = real array of dimension abd(nabd,n).
!         on return contains the values of the matrix stored in
!         banded form. The j-th column of abd contains the elements
!         of the j-th column of  the original matrix comprised in the
!         band ( i in (j-ml,j+mu) ) with the lowest diagonal at
!         the bottom row (row lowd). See details below for this format.
!
! ml      = integer. equal to the bandwidth of the strict lower part of A
! mu      = integer. equal to the bandwidth of the strict upper part of A
!         if job=1 on entry then these two values are internally computed.
!
! lowd  = integer. row number in abd where the lowest diagonal
!         (leftmost) of A is located on return. In case lowd = 0
!         on return, then it is defined to ml+mu+1 on return and the
!         lowd will contain this value on return. `
!
! ierr  = integer. used for error messages. On return:
!         ierr .eq. 0  :means normal return
!         ierr .eq. -1 : means invalid value for lowd. (either .lt. 0
!         or larger than nabd).
!         ierr .eq. -2 : means that lowd is not large enough and as
!         result the matrix cannot be stored in array abd.
!         lowd should be at least ml+mu+1, where ml and mu are as
!         provided on output.
!
!----------------------------------------------------------------------*
! Additional details on banded format.  (this closely follows the      *
! format used in linpack. may be useful for converting a matrix into   *
! this storage format in order to use the linpack  banded solvers).    *
!----------------------------------------------------------------------*
!             ---  band storage format  for matrix abd ---             *
! uses ml+mu+1 rows of abd(nabd,*) to store the diagonals of           *
! a in rows of abd starting from the lowest (sub)-diagonal  which  is  *
! stored in row number lowd of abd. the minimum number of rows needed  *
! in abd is ml+mu+1, i.e., the minimum value for lowd is ml+mu+1. the  *
! j-th  column  of  abd contains the elements of the j-th column of a, *
! from bottom to top: the element a(j+ml,j) is stored in  position     *
! abd(lowd,j), then a(j+ml-1,j) in position abd(lowd-1,j) and so on.   *
! Generally, the element a(j+k,j) of original matrix a is stored in    *
! position abd(lowd+k-ml,j), for k=ml,ml-1,..,0,-1, -mu.               *
! The first dimension nabd of abd must be .ge. lowd                    *
!                                                                      *
!     example [from linpack ]:   if the original matrix is             *
!                                                                      *
!              11 12 13  0  0  0                                       *
!              21 22 23 24  0  0                                       *
!               0 32 33 34 35  0     original banded matrix            *
!               0  0 43 44 45 46                                       *
!               0  0  0 54 55 56                                       *
!               0  0  0  0 65 66                                       *
!                                                                      *
! then  n = 6, ml = 1, mu = 2. lowd should be .ge. 4 (=ml+mu+1)  and   *
! if lowd = 5 for example, abd  should be:                             *
!                                                                      *
! untouched --> x  x  x  x  x  x                                       *
!               *  * 13 24 35 46                                       *
!               * 12 23 34 45 56    resulting abd matrix in banded     *
!              11 22 33 44 55 66    format                             *
!  row lowd--> 21 32 43 54 65  *                                       *
!                                                                      *
! * = not used                                                         *
!
!
!----------------------------------------------------------------------*
! first determine ml and mu.
!-----------------------------------------------------------------------
      ierr = 0		
!-----------
      if (job .eq. 1) call getbwd(n,a,ja,ia,ml,mu)
      m = ml+mu+1
      if (lowd .eq. 0) lowd = m
      if (m .gt. lowd)  ierr = -2
      if (lowd .gt. nabd .or. lowd .lt. 0) ierr = -1
      if (ierr .lt. 0) return
!------------
      do 15  i=1,m
         ii = lowd -i+1
         do 10 j=1,n
          abd(ii,j) = 0.0d0
 10      continue
 15   continue
!---------------------------------------------------------------------
      mdiag = lowd-ml
      do 30 i=1,n
         do 20 k=ia(i),ia(i+1)-1
            j = ja(k)
            abd(i-j+mdiag,j) = a(k)
 20      continue
 30   continue
      return
      end subroutine csrbnd
! ...................................................

      END MODULE SPI_SPARSEKIT

