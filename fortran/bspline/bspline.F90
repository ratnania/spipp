! -*- f90 -*-

! This file contains two modules
! BSPLINE : basic functions and routines
! BSP     : specific implementations for 1D,2D,3D
! Externaly, the user should call BSP

module bspline
contains

! .......................................................
function FindSpan(n,p,uu,U) result (span)
!        Determine the knot span index 
!        Input: n,p,u,U
!        Return: the knot span index 
  implicit none
  integer(kind=4), intent(in) :: n, p
  real   (kind=8), intent(in) :: uu, U(0:n+p+1)
  integer(kind=4)             :: span
  integer(kind=4) low, high
  if (uu >= U(n+1)) then
     span = n
     return
  end if
  if (uu <= U(p)) then
     span = p
     return
  end if
  low  = p
  high = n+1
  span = (low + high) / 2
  do while (uu < U(span) .or. uu >= U(span+1))
     if (uu < U(span)) then
        high = span
     else
        low  = span
     end if
     span = (low + high) / 2
  end do
end function FindSpan
! .......................................................

! .......................................................
function FindMult(i,uu,p,U) result (mult)
!
!        Input: 
!        Return
  implicit none
  integer(kind=4), intent(in)  :: i, p
  real   (kind=8), intent(in)  :: uu, U(0:i+p+1)
  integer(kind=4)              :: mult
  integer(kind=4) :: j
  mult = 0
  do j = -p, p+1
     if (uu == U(i+j)) mult = mult + 1
  end do
end function FindMult
! .......................................................

! .......................................................
subroutine FindSpanMult(n,p,uu,U,k,s)
!
!        Input: 
!        Return
  implicit none
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: uu, U(0:n+p+1)
  integer(kind=4), intent(out) :: k, s
  k = FindSpan(n,p,uu,U)
  s = FindMult(k,uu,p,U)
end subroutine FindSpanMult
! .......................................................

! .......................................................
subroutine BasisFuns(i,uu,p,U,N)
!          Compute the nonvanishing basis functions 
!          Input: i,u,p,U 
!          Output: N 
  implicit none
  integer(kind=4), intent(in) :: i, p
  real   (kind=8), intent(in) :: uu, U(0:i+p)
  real   (kind=8), intent(out):: N(0:p)
  integer(kind=4) :: j, r
  real   (kind=8) :: left(p), right(p), saved, temp
  N(0) = 1.0
  do j = 1, p
     left(j)  = uu - U(i+1-j)
     right(j) = U(i+j) - uu
     saved = 0.0
     do r = 0, j-1
        temp = N(r) / (right(r+1) + left(j-r))
        N(r) = saved + right(r+1) * temp
        saved = left(j-r) * temp
     end do
     N(j) = saved
  end do
end subroutine BasisFuns
! .......................................................

! .......................................................
subroutine DersBasisFuns(i,uu,p,n,U,ders)
!          Compute nonzero basis functions and their 
!          derivatives. First section is A2.2 (The NURBS Book) modified 
!          to store functions and knot differences. 
!          Input: i,u,p,n,U 
!          Output: ders  
  implicit none
  integer(kind=4), intent(in) :: i, p, n
  real   (kind=8), intent(in) :: uu, U(0:i+p)
  real   (kind=8), intent(out):: ders(0:p,0:n)
  integer(kind=4) :: j, k, r, s1, s2, rk, pk, j1, j2
  real   (kind=8) :: saved, temp, d
  real   (kind=8) :: left(p), right(p)
  real   (kind=8) :: ndu(0:p,0:p), a(0:1,0:p)
  ndu(0,0) = 1.0
  do j = 1, p
     left(j)  = uu - U(i+1-j)
     right(j) = U(i+j) - uu
     saved = 0.0
     do r = 0, j-1
        ndu(j,r) = right(r+1) + left(j-r)
        temp = ndu(r,j-1) / ndu(j,r)
        ndu(r,j) = saved + right(r+1) * temp
        saved = left(j-r) * temp
     end do
     ndu(j,j) = saved
  end do
  ders(:,0) = ndu(:,p)
  do r = 0, p
     s1 = 0; s2 = 1;
     a(0,0) = 1.0
     do k = 1, n
        d = 0.0
        rk = r-k; pk = p-k;
        if (r >= k) then
           a(s2,0) = a(s1,0) / ndu(pk+1,rk)
           d =  a(s2,0) * ndu(rk,pk)
        end if
        if (rk > -1) then
           j1 = 1
        else
           j1 = -rk
        end if
        if (r-1 <= pk) then
           j2 = k-1
        else
           j2 = p-r
        end if
        do j = j1, j2
           a(s2,j) = (a(s1,j) - a(s1,j-1)) / ndu(pk+1,rk+j)
           d =  d + a(s2,j) * ndu(rk+j,pk)
        end do
        if (r <= pk) then
           a(s2,k) = - a(s1,k-1) / ndu(pk+1,r)
           d =  d + a(s2,k) * ndu(r,pk)
        end if
        ders(r,k) = d
        j = s1; s1 = s2; s2 = j;
     end do
  end do
  r = p
  do k = 1, n
     ders(:,k) = ders(:,k) * r
     r = r * (p-k)
  end do
end subroutine DersBasisFuns
! .......................................................

! .......................................................
subroutine CurvePoint(d,n,p,U,Pw,uu,C)
!
!        Input: 
!        Output: 
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: uu
  real   (kind=8), intent(out) :: C(d)
  integer(kind=4) :: j, span
  real   (kind=8) :: basis(0:p)
  span = FindSpan(n,p,uu,U)
  call BasisFuns(span,uu,p,U,basis)
  C = 0.0
  do j = 0, p
     C = C + basis(j)*Pw(:,span-p+j)
  end do
end subroutine CurvePoint
! .......................................................

! .......................................................
subroutine SurfacePoint(d,n,p,U,m,q,V,Pw,uu,vv,S)
!
!        Input: 
!        Output: 
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  integer(kind=4), intent(in)  :: m, q
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: V(0:m+q+1)
  real   (kind=8), intent(in)  :: Pw(d,0:m,0:n)
  real   (kind=8), intent(in)  :: uu, vv
  real   (kind=8), intent(out) :: S(d)
  integer(kind=4) :: uj, vj, uspan, vspan
  real   (kind=8) :: ubasis(0:p), vbasis(0:q)
  uspan = FindSpan(n,p,uu,U)
  call BasisFuns(uspan,uu,p,U,ubasis)
  vspan = FindSpan(m,q,vv,V)
  call BasisFuns(vspan,vv,q,V,vbasis)
  S = 0.0
  do uj = 0, p
     do vj = 0, q
        S = S + ubasis(uj)*vbasis(vj)*Pw(:,vspan-q+vj,uspan-p+uj)
     end do
  end do
end subroutine SurfacePoint
! .......................................................

! .......................................................
subroutine CurvePntByCornerCut(d,n,p,U,Pw,x,Cw)
!
!        Input: 
!        Output: 
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: x
  real   (kind=8), intent(out) :: Cw(d)
  integer(kind=4) :: i, j, k, s, r
  real   (kind=8) :: uu, alpha, Rw(d,0:p)
  if (x <= U(p)) then
     uu = U(p)
     k = p
     s = FindMult(p,uu,p,U)
     if (s >= p) then
        Cw(:) = Pw(:,0)
        return
     end if
  elseif (x >= U(n+1)) then
     uu = U(n+1)
     k = n+1
     s = FindMult(n,uu,p,U)
     if (s >= p) then
        Cw(:) = Pw(:,n)
        return
     end if
  else
     uu = x
     k = FindSpan(n,p,uu,U)
     s = FindMult(k,uu,p,U)
     if (s >= p) then
        Cw(:) = Pw(:,k-p)
        return
     end if
  end if
  r = p-s
  do i = 0, r
     Rw(:,i) = Pw(:,k-p+i)
  end do
  do j = 1, r
     do i = 0, r-j
        alpha = (uu-U(k-p+j+i))/(U(i+k+1)-U(k-p+j+i))
        Rw(:,i) = alpha*Rw(:,i+1)+(1-alpha)*Rw(:,i)
     end do
  end do
  Cw(:) = Rw(:,0)
end subroutine CurvePntByCornerCut
! .......................................................

! .......................................................
subroutine InsertKnot(d,n,p,U,Pw,uu,k,s,r,V,Qw)
!
!        Input: 
!        Output: 
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: uu
  integer(kind=4), intent(in)  :: k, s, r
  real   (kind=8), intent(out) :: V(0:n+p+1+r)
  real   (kind=8), intent(out) :: Qw(d,0:n+r)
  integer(kind=4) :: i, j, idx
  real   (kind=8) :: alpha, Rw(d,0:p)
  ! Load new knot vector
  forall (i = 0:k) V(i) = U(i)
  forall (i = 1:r) V(k+i) = uu
  forall (i = k+1:n+p+1) V(i+r) = U(i)
  ! Save unaltered control points
  forall (i = 0:k-p) Qw(:,i)   = Pw(:,i)
  forall (i = k-s:n) Qw(:,i+r) = Pw(:,i)
  forall (i = 0:p-s) Rw(:,i)   = Pw(:,k-p+i)
  ! Insert the knot r times
  do j = 1, r
     idx = k-p+j
     do i = 0, p-j-s
        alpha = (uu-U(idx+i))/(U(i+k+1)-U(idx+i))
        Rw(:,i) = alpha*Rw(:,i+1)+(1-alpha)*Rw(:,i)
     end do
     Qw(:,idx) = Rw(:,0)
     Qw(:,k+r-j-s) = Rw(:,p-j-s)
  end do
  ! Load remaining control points
  idx = k-p+r
  do i = idx+1, k-s-1
     Qw(:,i) = Rw(:,i-idx)
  end do
end subroutine InsertKnot
! .......................................................

! .......................................................
subroutine RemoveKnot(d,n,p,U,Pw,uu,r,s,num,t,TOL)
!
!        Input: 
!        Output: 
  implicit none
  integer(kind=4), intent(in)    :: d
  integer(kind=4), intent(in)    :: n, p
  real   (kind=8), intent(inout) :: U(0:n+p+1)
  real   (kind=8), intent(inout) :: Pw(d,0:n)
  real   (kind=8), intent(in)    :: uu
  integer(kind=4), intent(in)    :: r, s, num
  integer(kind=4), intent(out)   :: t
  real   (kind=8), intent(in)    :: TOL

  integer(kind=4) :: m,ord,fout,last,first,off
  integer(kind=4) :: i,j,ii,jj,k
  logical         :: remflag
  real   (kind=8) :: temp(d,0:2*p)
  real   (kind=8) :: alfi,alfj

  m = n + p + 1
  ord = p + 1
  fout = (2*r-s-p)/2
  first = r - p
  last  = r - s
  do t = 0,num-1
     off = first - 1
     temp(:,0) = Pw(:,off)
     temp(:,last+1-off) = Pw(:,last+1)
     i = first; ii = 1
     j = last;  jj = last - off
     remflag = .false.
     do while (j-i > t)
        alfi = (uu-U(i))/(U(i+ord+t)-U(i))
        alfj = (uu-U(j-t))/(U(j+ord)-U(j-t))
        temp(:,ii) = (Pw(:,i)-(1.0-alfi)*temp(:,ii-1))/alfi
        temp(:,jj) = (Pw(:,j)-alfj*temp(:,jj+1))/(1.0-alfj)
        i = i + 1; ii = ii + 1
        j = j - 1; jj = jj - 1
     end do
     if (j-i < t) then
        if (Distance(d,temp(:,ii-1),temp(:,jj+1)) <= TOL) then
           remflag = .true.
        end if
     else
        alfi = (uu-U(i))/(U(i+ord+t)-U(i))
        if (Distance(d,Pw(:,i),alfi*temp(:,ii+t+1)+(1-alfi)*temp(:,ii-1)) <= TOL) then
           remflag = .true.
        end if
     end if
     if (remflag .eqv. .false.) then
        exit ! break out of the for loop
     else
        i = first
        j = last
        do while (j-i > t)
           Pw(:,i) = temp(:,i-off)
           Pw(:,j) = temp(:,j-off)
           i = i + 1
           j = j - 1
        end do
     end if
     first = first - 1
     last  = last  + 1
  end do
  if (t == 0) return
  do k = r+1,m
     U(k-t) = U(k)
  end do
  j = fout
  i = j
  do k = 1,t-1
     if (mod(k,2) == 1) then
        i = i + 1
     else
        j = j - 1
     end if
  end do
  do k = i+1,n
     Pw(:,j) = Pw(:,k)
     j = j + 1
  enddo
contains
  function Distance(d,P1,P2) result (dist)
    implicit none
    integer(kind=4), intent(in) :: d
    real   (kind=8), intent(in) :: P1(d),P2(d)
    integer(kind=4) :: i
    real   (kind=8) :: dist
    dist = 0.0
    do i = 1,d
       dist = dist + (P1(i)-P2(i))*(P1(i)-P2(i))
    end do
    dist = sqrt(dist)
  end function Distance
end subroutine RemoveKnot
! .......................................................

! .......................................................
subroutine ClampKnot(d,n,p,U,Pw,l,r)
!
!        Input: 
!        Output: 
  implicit none
  integer(kind=4), intent(in)    :: d
  integer(kind=4), intent(in)    :: n, p
  real   (kind=8), intent(inout) :: U(0:n+p+1)
  real   (kind=8), intent(inout) :: Pw(d,0:n)
  logical(kind=4), intent(in)    :: l, r
  integer(kind=4) :: k, s
  if (l) then ! Clamp at left end
     k = p
     s = FindMult(p,U(p),p,U)
     call KntIns(d,n,p,U,Pw,k,s)
     U(0:p-1) = U(p)
  end if
  if (r) then ! Clamp at right end
     k = n+1
     s = FindMult(n,U(n+1),p,U)
     call KntIns(d,n,p,U,Pw,k,s)
     U(n+2:n+p+1) = U(n+1)
  end if
contains
  subroutine KntIns(d,n,p,U,Pw,k,s)
      implicit none
      integer(kind=4), intent(in)    :: d
      integer(kind=4), intent(in)    :: n, p
      real   (kind=8), intent(in)    :: U(0:n+p+1)
      real   (kind=8), intent(inout) :: Pw(d,0:n)
      integer(kind=4), intent(in)    :: k, s
      integer(kind=4) :: r, i, j, idx
      real   (kind=8) :: uu, alpha, Rw(d,0:p), Qw(d,0:2*p)
      if (s >= p) return
      uu = U(k)
      r = p-s
      Qw(:,0) = Pw(:,k-p)
      Rw(:,0:p-s) = Pw(:,k-p:k-s)
      do j = 1, r
         idx = k-p+j
         do i = 0, p-j-s
            alpha = (uu-U(idx+i))/(U(i+k+1)-U(idx+i))
            Rw(:,i) = alpha*Rw(:,i+1)+(1-alpha)*Rw(:,i)
         end do
         Qw(:,j) = Rw(:,0)
         Qw(:,p-j-s+r) = Rw(:,p-j-s)
      end do
      if (k == p) then ! left end
         Pw(:,0:r-1) = Qw(:,r:r+r-1)
      else             ! right end
         Pw(:,n-r+1:n) = Qw(:,p-r:p-1)
      end if
    end subroutine KntIns
end subroutine ClampKnot
! .......................................................

! .......................................................
subroutine UnclampKnot(d,n,p,U,Pw,l,r)
!
!        Input: 
!        Output: 
  implicit none
  integer(kind=4), intent(in)    :: d
  integer(kind=4), intent(in)    :: n, p
  real   (kind=8), intent(inout) :: U(0:n+p+1)
  real   (kind=8), intent(inout) :: Pw(d,0:n)
  logical(kind=4), intent(in)    :: l, r
  integer(kind=4) :: i, j, k
  real   (kind=8) :: alpha
  if (l) then ! Unclamp at left end
     do i = 0, p-2
        U(p-i-1) = U(p-i) - (U(n-i+1)-U(n-i))
        k = p-1
        do j = i, 0, -1
           alpha = (U(p)-U(k))/(U(p+j+1)-U(k))
           Pw(:,j) = (Pw(:,j)-alpha*Pw(:,j+1))/(1-alpha)
           k = k-1
        end do
     end do
     U(0) = U(1) - (U(n-p+2)-U(n-p+1)) ! Set first knot
  end if
  if (r) then ! Unclamp at right end
     do i = 0, p-2
        U(n+i+2) = U(n+i+1) + (U(p+i+1)-U(p+i))
        do j = i, 0, -1
           alpha = (U(n+1)-U(n-j))/(U(n-j+i+2)-U(n-j))
           Pw(:,n-j) = (Pw(:,n-j)-(1-alpha)*Pw(:,n-j-1))/alpha
        end do
     end do
     U(n+p+1) = U(n+p) + (U(2*p)-U(2*p-1)) ! Set last knot
  end if
end subroutine UnclampKnot
! .......................................................

! .......................................................
subroutine RefineKnotVector(d,n,p,U,Pw,r,X,Ubar,Qw)
!
!        Input: 
!        Output: 
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: r
  real   (kind=8), intent(in)  :: X(0:r)
  real   (kind=8), intent(out) :: Ubar(0:n+r+1+p+1)
  real   (kind=8), intent(out) :: Qw(d,0:n+r+1)
  integer(kind=4) :: m, a, b
  integer(kind=4) :: i, j, k, l
  integer(kind=4) :: idx
  real   (kind=8) :: alpha
  if (r < 0) then
     Ubar = U
     Qw = Pw
     return
  end if
  m = n + p + 1
  a = FindSpan(n,p,X(0),U)
  b = FindSpan(n,p,X(r),U)
  b = b + 1
  forall (j = 0:a-p) Qw(:,j)     = Pw(:,j)
  forall (j = b-1:n) Qw(:,j+r+1) = Pw(:,j)
  forall (j =   0:a) Ubar(j)     = U(j)
  forall (j = b+p:m) Ubar(j+r+1) = U(j)
  i = b + p - 1
  k = b + p + r
  do j = r, 0, -1
     do while (X(j) <= U(i) .and. i > a)
        Qw(:,k-p-1) = Pw(:,i-p-1)
        Ubar(k) = U(i)
        k = k - 1
        i = i - 1
     end do
     Qw(:,k-p-1) = Qw(:,k-p)
     do l = 1, p
        idx = k - p + l
        alpha = Ubar(k+l) - X(j)
        if (abs(alpha) == 0.0) then
           Qw(:,idx-1) = Qw(:,idx)
        else
           alpha = alpha / (Ubar(k+l) - U(i-p+l))
           Qw(:,idx-1) = alpha*Qw(:,idx-1) + (1-alpha)*Qw(:,idx)
        end if
     end do
     Ubar(k) = X(j)
     k = k-1
  end do
end subroutine RefineKnotVector
! .......................................................

! .......................................................
subroutine DegreeElevate(d,n,p,U,Pw,t,nh,Uh,Qw)
!
!        Input: 
!        Output: 
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: t
  integer(kind=4), intent(in)  :: nh
  real   (kind=8), intent(out) :: Uh(0:nh+p+t+1)
  real   (kind=8), intent(out) :: Qw(d,0:nh)

  integer(kind=4) :: i, j, k, kj, tr, a, b
  integer(kind=4) :: m, ph, kind, cind, first, last
  integer(kind=4) :: r, oldr, s, mul, lbz, rbz

  real   (kind=8) :: bezalfs(0:p+t,0:p)
  real   (kind=8) :: bpts(d,0:p), ebpts(d,0:p+t), nextbpts(d,0:p-2)
  real   (kind=8) :: alfs(0:p-2), ua, ub, alf, bet, gam, den
  if (t < 1) then
     Uh = U
     Qw = Pw
     return
  end if
  m = n + p + 1
  ph = p + t
  ! Bezier coefficients
  bezalfs(0,0)  = 1.0
  bezalfs(ph,p) = 1.0
  do i = 1, ph/2
     do j = max(0,i-t), min(p,i)
        bezalfs(i,j) = Bin(p,j)*Bin(t,i-j)*(1.0d+0/Bin(ph,i))
     end do
  end do
  do i = ph/2+1, ph-1
     do j = max(0,i-t), min(p,i)
        bezalfs(i,j) = bezalfs(ph-i,p-j)
     end do
  end do
  kind = ph+1
  cind = 1
  r = -1
  a = p
  b = p+1
  ua = U(a)
  Uh(0:ph) = ua
  Qw(:,0) = Pw(:,0)
  bpts = Pw(:,0:p)
  do while (b < m)
     i = b
     do while (b < m)
        if (U(b) /= U(b+1)) exit
        b = b + 1
     end do
     mul = b - i + 1
     oldr = r
     r = p - mul
     ub = U(b)
     if (oldr > 0) then
        lbz = (oldr+2)/2
     else
        lbz = 1
     end if
     if (r > 0) then
        rbz = ph - (r+1)/2
     else
        rbz = ph
     end if
     ! insert knots
     if (r > 0) then
        do k = p, mul+1, -1
           alfs(k-mul-1) = (ub-ua)/(U(a+k)-ua)
        end do
        do j = 1, r
           s = mul + j
           do k = p, s, -1
              bpts(:,k) = alfs(k-s)  * bpts(:,k) + &
                     (1.0-alfs(k-s)) * bpts(:,k-1)
           end do
           nextbpts(:,r-j) = bpts(:,p)
        end do
     end if
     ! degree elevate
     do i = lbz, ph
        ebpts(:,i) = 0.0
        do j = max(0,i-t), min(p,i)
           ebpts(:,i) = ebpts(:,i) + bezalfs(i,j)*bpts(:,j)
        end do
     end do
     ! remove knots
     if (oldr > 1) then
        first = kind-2
        last = kind
        den = ub-ua
        bet = (ub-Uh(kind-1))/den
        do tr = 1, oldr-1
           i = first
           j = last
           kj = j-kind+1
           do while (j-i > tr)
              if (i < cind) then
                 alf = (ub-Uh(i))/(ua-Uh(i))
                 Qw(:,i) = alf*Qw(:,i) + (1.0-alf)*alf*Qw(:,i-1)
              end if
              if (j >= lbz) then
                 if (j-tr <= kind-ph+oldr) then
                    gam = (ub-Uh(j-tr))/den
                    ebpts(:,kj) = gam*ebpts(:,kj) + (1.0-gam)*ebpts(:,kj+1)
                 else
                    ebpts(:,kj) = bet*ebpts(:,kj) + (1.0-bet)*ebpts(:,kj+1)
                 end if
              end if
              i = i+1
              j = j-1
              kj = kj-1
           end do
           first = first-1
           last = last+1
        end do
     end if
     !
     if (a /= p) then
        do i = 0, ph-oldr-1
           Uh(kind) = ua
           kind = kind+1
        end do
     end if
     do j = lbz, rbz
        Qw(:, cind) = ebpts(:,j)
        cind = cind+1
     end do
     !
     if (b < m) then
        bpts(:,0:r-1) = nextbpts(:,0:r-1)
        bpts(:,r:p) = Pw(:,b-p+r:b)
        a = b
        b = b+1
        ua = ub
     else
        Uh(kind:kind+ph) = ub
     end if
  end do
contains
  pure function Bin(n,k) result (C)
    implicit none
    integer(kind=4), intent(in) :: n, k
    integer(kind=4) :: i, C
    C = 1
    do i = 0, min(k,n-k) - 1
       C = C * (n - i)
       C = C / (i + 1)
    end do
  end function Bin
end subroutine DegreeElevate
! .......................................................

end module bspline

!
! ----------
!

module BSp
contains

! .......................................................
subroutine FindSpan(p,m,U,uu,span)
!
!        Input: 
!        Output: 
  use bspline, FindS => FindSpan
  implicit none
  integer(kind=4), intent(in)  :: p, m
  real   (kind=8), intent(in)  :: U(0:m), uu
  integer(kind=4), intent(out) :: span
  span = FindS(m-(p+1),p,uu,U)
end subroutine FindSpan
! .......................................................

! .......................................................
subroutine FindMult(p,m,U,uu,span,mult)
!
!        Input: 
!        Output: 
  use bspline, FindM => FindMult
  implicit none
  integer(kind=4), intent(in)  :: p, m
  real   (kind=8), intent(in)  :: U(0:m), uu
  integer(kind=4), intent(in)  :: span
  integer(kind=4), intent(out) :: mult
  integer(kind=4) :: k
  if (span >= 0) then
     k = span
  else
     k = FindSpan(m-(p+1),p,uu,U)
  end if
  mult = FindM(k,uu,p,U)
end subroutine FindMult
! .......................................................

! .......................................................
subroutine FindSpanMult(p,m,U,uu,k,s)
!
!        Input: 
!        Output: 
  use bspline, FindSM => FindSpanMult
  implicit none
  integer(kind=4), intent(in)  :: p, m
  real   (kind=8), intent(in)  :: U(0:m), uu
  integer(kind=4), intent(out) :: k, s
  call FindSM(m-(p+1),p,uu,U,k,s)
end subroutine FindSpanMult
! .......................................................

! .......................................................
subroutine EvalBasisFuns(p,m,U,uu,span,N)
!
!        Input: 
!        Output: 
  use bspline
  implicit none
  integer(kind=4), intent(in) :: p, m, span
  real   (kind=8), intent(in) :: U(0:m), uu
  real   (kind=8), intent(out):: N(0:p)
  integer(kind=4) :: i
  if (span >= 0) then
     i = span
  else
     i = FindSpan(m-(p+1),p,uu,U)
  end if
  call BasisFuns(i,uu,p,U,N)
end subroutine EvalBasisFuns
! .......................................................

! .......................................................
subroutine EvalBasisFunsDers(p,m,U,uu,d,span,dN)
!
!        Input: 
!        Output: 
  use bspline
  implicit none
  integer(kind=4), intent(in) :: p, m, d, span
  real   (kind=8), intent(in) :: U(0:m), uu
  real   (kind=8), intent(out):: dN(0:p,0:d)
  integer(kind=4) :: i
  if (span >= 0) then
     i = span
  else
     i = FindSpan(m-(p+1),p,uu,U)
  end if
  call DersBasisFuns(i,uu,p,d,U,dN)
end subroutine EvalBasisFunsDers
! .......................................................

! .......................................................
subroutine SpanIndex(p,m,U,r,I)
!
!        Input: 
!        Output: 
  integer(kind=4), intent(in)  :: p, m
  real   (kind=8), intent(in)  :: U(0:m)
  integer(kind=4), intent(in)  :: r
  integer(kind=4), intent(out) :: I(r)
  integer(kind=4) :: k, s
  s = 1
  do k = p, m-(p+1)
     if (U(k) /= U(k+1)) then
        I(s) = k; s = s + 1
        if (s > r) exit
     end if
  end do
end subroutine SpanIndex
! .......................................................

! .......................................................
subroutine Greville(p,m,U,X)
!
!        Input: 
!        Output: 
  implicit none
  integer(kind=4), intent(in)  :: p, m
  real   (kind=8), intent(in)  :: U(0:m)
  real   (kind=8), intent(out) :: X(0:m-(p+1))
  integer(kind=4) :: i
  do i = 0, m-(p+1)
     X(i) = sum(U(i+1:i+p)) / p
  end do
end subroutine Greville
! .......................................................

! .......................................................
subroutine InsertKnot(d,n,p,U,Pw,uu,r,V,Qw)
!
!        Input: 
!        Output: 
  use bspline, InsKnt => InsertKnot
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: uu
  integer(kind=4), intent(in)  :: r
  real   (kind=8), intent(out) :: V(0:n+p+1+r)
  real   (kind=8), intent(out) :: Qw(d,0:n+r)
  integer(kind=4) :: k, s
  if (r == 0) then
     V = U; Qw = Pw; return
  end if
  call FindSpanMult(n,p,uu,U,k,s)
  call InsKnt(d,n,p,U,Pw,uu,k,s,r,V,Qw)
end subroutine InsertKnot
! .......................................................

! .......................................................
subroutine RemoveKnot(d,n,p,U,Pw,uu,r,t,V,Qw,TOL)
!
!        Input: 
!        Output: 
  use bspline, RemKnt => RemoveKnot
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: uu
  integer(kind=4), intent(in)  :: r
  integer(kind=4), intent(out) :: t
  real   (kind=8), intent(out) :: V(0:n+p+1)
  real   (kind=8), intent(out) :: Qw(d,0:n)
  real   (kind=8), intent(in)  :: TOL
  integer(kind=4) :: k, s
  t = 0
  V = U
  Qw = Pw
  if (r == 0) return
  if (uu <= U(p)) return
  if (uu >= U(n+1)) return
  call FindSpanMult(n,p,uu,U,k,s)
  call RemKnt(d,n,p,V,Qw,uu,k,s,r,t,TOL)
end subroutine RemoveKnot
! .......................................................

! .......................................................
subroutine Clamp(d,n,p,U,Pw,l,r,V,Qw)
!
!        Input: 
!        Output: 
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  logical(kind=4), intent(in)  :: l, r
  real   (kind=8), intent(out) :: V(0:n+p+1)
  real   (kind=8), intent(out) :: Qw(d,0:n)
  V  = U
  Qw = Pw
  call ClampKnot(d,n,p,V,Qw,l,r)
end subroutine Clamp
! .......................................................

! .......................................................
subroutine Unclamp(d,n,p,U,Pw,l,r,V,Qw)
!
!        Input: 
!        Output: 
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  logical(kind=4), intent(in)  :: l, r
  real   (kind=8), intent(out) :: V(0:n+p+1)
  real   (kind=8), intent(out) :: Qw(d,0:n)
  V  = U
  Qw = Pw
  call UnclampKnot(d,n,p,V,Qw,l,r)
end subroutine Unclamp
! .......................................................

! .......................................................
subroutine RefineKnotVector(d,n,p,U,Pw,r,X,Ubar,Qw)
!
!        Input: 
!        Output: 
  use bspline, RefKnt => RefineKnotVector
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: r
  real   (kind=8), intent(in)  :: X(0:r)
  real   (kind=8), intent(out) :: Ubar(0:n+r+1+p+1)
  real   (kind=8), intent(out) :: Qw(d,0:n+r+1)
  call RefKnt(d,n,p,U,Pw,r,X,Ubar,Qw)
end subroutine RefineKnotVector
! .......................................................

! .......................................................
subroutine DegreeElevate(d,n,p,U,Pw,t,nh,Uh,Qw)
!
!        Input: 
!        Output: 
  use bspline, DegElev => DegreeElevate
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: t
  integer(kind=4), intent(in)  :: nh
  real   (kind=8), intent(out) :: Uh(0:nh+p+t+1)
  real   (kind=8), intent(out) :: Qw(d,0:nh)
  call DegElev(d,n,p,U,Pw,t,nh,Uh,Qw)
end subroutine DegreeElevate
! .......................................................

! .......................................................
subroutine Extract(d,n,p,U,Pw,x,Cw)
!
!        Input: 
!        Output: 
  use bspline, CornerCut => CurvePntByCornerCut
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  real   (kind=8), intent(in)  :: x
  real   (kind=8), intent(out) :: Cw(d)
  call CornerCut(d,n,p,U,Pw,x,Cw)
end subroutine Extract
! .......................................................

! .......................................................
subroutine AssembleBasis1(nderiv,N,rationalize,nx,px,U,Ww,rx,X,Basis)
!
!        Input: 
!        Output: 
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: nderiv
  integer(kind=4), intent(in)  :: N   
  integer(kind=4), intent(in)  :: rationalize 
  integer(kind=4), intent(in)  :: nx, px
  real   (kind=8), intent(in)  :: U(0:nx+px+1)
  real   (kind=8), intent(in)  :: Ww(0:nx)
  integer(kind=4), intent(in)  :: rx
  real   (kind=8), intent(in)  :: X(0:rx)
  real   (kind=8), intent(out) :: Basis(0:N,1:(px+1),1:(rx+1)) 
  integer(kind=4) :: i, j, span, o, indX, indB
  real   (kind=8) :: dbasis(0:px,0:nderiv), w(0:nderiv)
  real   (kind=8) :: Rdbasis(0:nderiv)
  real   (kind=8) :: weight 
  !

  do i = 0, rx
     span = FindSpan(nx,px,X(i),U)
     call DersBasisFuns(span,X(i),px,nderiv,U,dbasis)

     o = span - px
     indX = i


     !
     ! compute w = sum wi Ni
     ! and w' = sum wi Ni'
     w  = 0.0
     do j = 0, px
        weight = Ww(o+j)
        w  = w  + weight * dbasis(j,:)
     end do
     ! compute Nurbs
     Rdbasis  = 0.0
     do j = 0, px
        weight = Ww(o+j)
        indB = j
        Rdbasis(0) = dbasis(j,0) / w(0)

        if (nderiv >= 1) then
           if (rationalize==1) then              
              Rdbasis(1) = dbasis(j,1) / w(0)               &
                         - dbasis(j,0) * w(1) / w(0)**2 
           else
              Rdbasis(1) = dbasis(j,1)
           end if

        end if

        if (nderiv>=2) then
           if (rationalize==1) then              
              Rdbasis(2) = dbasis(j,2) / w(0)               &
                         - dbasis(j,0) * w(2) / w(0)**2     &
                         - 2 * Rdbasis(1) * w(1) / w(0) 
           else
              Rdbasis(2) = dbasis(j,2)
           end if
        end if
     
        Basis(0:N, indB+1, indX+1) = weight * Rdbasis(0:N)
     end do

     !
  end do
  !
end subroutine AssembleBasis1
! .......................................................

! .......................................................
subroutine Evaluate1(d,n,p,U,Pw,r,X,Cw)
!
!        Input: 
!        Output: 
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: n, p
  real   (kind=8), intent(in)  :: U(0:n+p+1)
  real   (kind=8), intent(in)  :: Pw(d,0:n)
  integer(kind=4), intent(in)  :: r
  real   (kind=8), intent(in)  :: X(0:r)
  real   (kind=8), intent(out) :: Cw(d,0:r)
  integer(kind=4) :: i, j, span
  real   (kind=8) :: basis(0:p), C(d)
  !
  do i = 0, r
     span = FindSpan(n,p,X(i),U)
     call BasisFuns(span,X(i),p,U,basis)
     !
     C = 0.0
     do j = 0, p
        C = C + basis(j)*Pw(:,span-p+j)
     end do
     Cw(:,i) = C
     !
  end do
  !
end subroutine Evaluate1
! .......................................................

! .......................................................
subroutine EvaluateDeriv1(nderiv,N,rationalize,d,nx,px,U,Pw,rx,X,Cw)
!
!        Input: 
!        Output: 
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: nderiv
  integer(kind=4), intent(in)  :: N   
  integer(kind=4), intent(in)  :: rationalize   
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: nx, px
  real   (kind=8), intent(in)  :: U(0:nx+px+1)
  real   (kind=8), intent(in)  :: Pw(d,0:nx)
  integer(kind=4), intent(in)  :: rx
  real   (kind=8), intent(in)  :: X(0:rx)
  real   (kind=8), intent(out) :: Cw(d,0:rx,0:N)
  integer(kind=4) :: i, j, span, deriv
  real   (kind=8) :: dbasis(0:px,0:nderiv), C(d), w(0:nderiv)
  real   (kind=8) :: Rdbasis(0:px,0:nderiv)
  real   (kind=8) :: basis(0:px)
  !

  do i = 0, rx
     span = FindSpan(nx,px,X(i),U)
     call DersBasisFuns(span,X(i),px,nderiv,U,dbasis)

     !
     ! compute w = sum wi Ni
     ! and w' = sum wi Ni'
     w  = 0.0
     do j = 0, px
        w  = w  + dbasis(j,:)*Pw(d,span-px+j)
     end do
     ! compute Nurbs
     Rdbasis  = 0.0
     Rdbasis(:,0) = dbasis(:,0) / w(0)

     if (nderiv >= 1) then
        if (rationalize==1) then              
           Rdbasis(:,1) = dbasis(:,1) / w(0)               &
                        - dbasis(:,0) * w(1) / w(0)**2 
        else
           Rdbasis(:,1) = dbasis(:,1)
        end if
     end if

     if (nderiv>=2) then
        if (rationalize==1) then              
           Rdbasis(:,2) = dbasis(:,2) / w(0)               &
                        - dbasis(:,0) * w(2) / w(0)**2     &
                        - 2 * Rdbasis(:,1) * w(1) / w(0) 
        else
           Rdbasis(:,2) = dbasis(:,2)
        end if                        
     end if

     do deriv = 0, N
        C = 0.0
        do j = 0, px
           C = C + Rdbasis(j,deriv)*Pw(:,span-px+j)
        end do
        Cw(:,i, deriv) = C
     end do
     !
  end do
  !
end subroutine EvaluateDeriv1
! .......................................................

! .......................................................
subroutine AssembleBasis2(nderiv,N,rationalize,nx,px,Ux,ny,py,Uy,Ww,rx,X,ry,Y,Basis)
!
!        Input: 
!        Output: 
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: nderiv
  integer(kind=4), intent(in)  :: N 
  integer(kind=4), intent(in)  :: rationalize 
  integer(kind=4), intent(in)  :: nx, ny
  integer(kind=4), intent(in)  :: px, py
  integer(kind=4), intent(in)  :: rx, ry
  real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
  real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
  real   (kind=8), intent(in)  :: Ww(0:ny,0:nx)
  real   (kind=8), intent(in)  :: X(0:rx), Y(0:ry)
  real   (kind=8), intent(out) :: Basis(0:N,1:(px+1)*(py+1),1:(rx+1)*(ry+1)) 
  integer(kind=4) :: ix, jx, iy, jy, ox, oy, deriv, indX, indB
  integer(kind=4) :: spanx(0:rx), spany(0:ry)
  real   (kind=8) :: dbasisx(0:px,0:nderiv,0:rx)
  real   (kind=8) :: dbasisy(0:py,0:nderiv,0:ry)
  ! Rdbasis(0) => Rij
  ! Rdbasis(1) => dx Rij
  ! Rdbasis(2) => dy Rij
  real   (kind=8) :: Rdbasis(0:N)  
  real   (kind=8) :: weight 
  real   (kind=8) :: M, Mx, My, Mxy, Mxx, Myy
  real   (kind=8) :: w, wx, wy, wxy, wxx, wyy

  Basis = 0.0

  !
  do ix = 0, rx
     spanx(ix) = FindSpan(nx,px,X(ix),Ux)
     call DersBasisFuns(spanx(ix),X(ix),px,nderiv,Ux,dbasisx(:,0:nderiv,ix))
  end do
  do iy = 0, ry
     spany(iy) = FindSpan(ny,py,Y(iy),Uy)
     call DersBasisFuns(spany(iy),Y(iy),py,nderiv,Uy,dbasisy(:,0:nderiv,iy))
  end do

  !
  ! compute 
  ! w   = sum wij Ni   Nj
  ! wx  = sum wij Ni'  Nj
  ! wy  = sum wij Ni   Nj'
  ! wxx = sum wij Ni'' Nj
  ! wxy = sum wij Ni'  Nj'
  ! wyy = sum wij Ni   Nj''
  do ix = 0, rx
  ox = spanx(ix) - px
  do iy = 0, ry
     oy = spany(iy) - py

     indX = (rx+1) * iy + ix

     ! --- compute w and its Derivatives
     w   = 0.0 ; wx  = 0.0 ; wy  = 0.0
     wxx = 0.0 ; wxy = 0.0 ; wyy = 0.0
     do jx = 0, px
     do jy = 0, py
        weight = Ww(oy+jy,ox+jx)

        M   = dbasisx(jx,0,ix) * dbasisy(jy,0,iy)
        w  = w  + M   * weight

        if (nderiv >= 1) then
           Mx  = dbasisx(jx,1,ix) * dbasisy(jy,0,iy)
           My  = dbasisx(jx,0,ix) * dbasisy(jy,1,iy)
           wx = wx + Mx  * weight
           wy = wy + My  * weight
        end if

        if (nderiv >= 2) then
           Mxx = dbasisx(jx,2,ix) * dbasisy(jy,0,iy)
           Mxy = dbasisx(jx,1,ix) * dbasisy(jy,1,iy)
           Myy = dbasisx(jx,0,ix) * dbasisy(jy,2,iy)

           wxx = wxx + Mxx * weight
           wxy = wxy + Mxy * weight
           wyy = wyy + Myy * weight
        end if 
     end do
     end do
     ! ---

     ! compute Nurbs and their derivatives
        do jx = 0, px
        do jy = 0, py     
           indB = (px+1) * jy + jx

           weight = Ww(oy+jy,ox+jx)
           
           M   = dbasisx(jx,0,ix) * dbasisy(jy,0,iy)
           Rdbasis(0) = M  / w

           if (nderiv >= 1) then
              Mx  = dbasisx(jx,1,ix) * dbasisy(jy,0,iy)
              My  = dbasisx(jx,0,ix) * dbasisy(jy,1,iy)
              if (rationalize==1) then              
                 Rdbasis(1) = Mx / w - M * wx / w**2 
                 Rdbasis(2) = My / w - M * wy / w**2 
              else
                 Rdbasis(1) = Mx
                 Rdbasis(2) = My
              end if              
           end if
        
           if (nderiv >= 2) then
              Mxx = dbasisx(jx,2,ix) * dbasisy(jy,0,iy)
              Mxy = dbasisx(jx,1,ix) * dbasisy(jy,1,iy)
              Myy = dbasisx(jx,0,ix) * dbasisy(jy,2,iy)
              if (rationalize==1) then              
                 Rdbasis(3) = Mxx / w                 &
                            - 2 * Mx * wx / w**2      &
                            - M * wxx / w**2          &
                            + 2 * M * wx**2 / w**3
              
                 Rdbasis(4) = Mxy / w                 &
                            - Mx * wy / w**2          &
                            - My * wx / w**2          &
                            - M * wxy / w**2          &
                            + 2 * M * wx * wy / w**3
              
                 Rdbasis(5) = Myy / w                 &
                            - 2 * My * wy / w**2      &
                            - M * wyy / w**2          &
                            + 2 * M * wy**2 / w**3
              else
                 Rdbasis(3) = Mxx
                 Rdbasis(4) = Mxy
                 Rdbasis(5) = Myy
              end if     
           end if

           Basis(0:N, indB+1, indX+1) = weight * Rdbasis(0:N)

        end do
        end do
     ! ---

  end do
  end do
  !
end subroutine AssembleBasis2
! .......................................................

! .......................................................
subroutine Evaluate2(d,nx,px,Ux,ny,py,Uy,Pw,rx,X,ry,Y,Cw)
!
!        Input: 
!        Output: 
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: nx, ny
  integer(kind=4), intent(in)  :: px, py
  integer(kind=4), intent(in)  :: rx, ry
  real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
  real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
  real   (kind=8), intent(in)  :: Pw(d,0:ny,0:nx)
  real   (kind=8), intent(in)  :: X(0:rx), Y(0:ry)
  real   (kind=8), intent(out) :: Cw(d,0:ry,0:rx)
  integer(kind=4) :: ix, jx, iy, jy, ox, oy
  integer(kind=4) :: spanx(0:rx), spany(0:ry)
  real   (kind=8) :: basisx(0:px,0:rx), basisy(0:py,0:ry)
  real   (kind=8) :: M, C(d)
  !
  do ix = 0, rx
     spanx(ix) = FindSpan(nx,px,X(ix),Ux)
     call BasisFuns(spanx(ix),X(ix),px,Ux,basisx(:,ix))
  end do
  do iy = 0, ry
     spany(iy) = FindSpan(ny,py,Y(iy),Uy)
     call BasisFuns(spany(iy),Y(iy),py,Uy,basisy(:,iy))
  end do
  !
  do ix = 0, rx
     ox = spanx(ix) - px
     do iy = 0, ry
        oy = spany(iy) - py
        ! ---
        C = 0.0
        do jx = 0, px
           do jy = 0, py
              M = basisx(jx,ix) * basisy(jy,iy)
              C = C + M * Pw(:,oy+jy,ox+jx)
           end do
        end do
        Cw(:,iy,ix) = C
        ! ---
     end do
  end do
  !
end subroutine Evaluate2
! .......................................................

! .......................................................
subroutine EvaluateDeriv2(nderiv,N,rationalize,d,nx,px,Ux,ny,py,Uy,Pw,rx,X,ry,Y,Cw)
!
!        Input: 
!        Output: 
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: nderiv
  integer(kind=4), intent(in)  :: N 
  integer(kind=4), intent(in)  :: rationalize 
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: nx, ny
  integer(kind=4), intent(in)  :: px, py
  integer(kind=4), intent(in)  :: rx, ry
  real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
  real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
  real   (kind=8), intent(in)  :: Pw(d,0:ny,0:nx)
  real   (kind=8), intent(in)  :: X(0:rx), Y(0:ry)
  real   (kind=8), intent(out) :: Cw(d,0:ry,0:rx,0:N)
  integer(kind=4) :: ix, jx, iy, jy, ox, oy, deriv
  integer(kind=4) :: spanx(0:rx), spany(0:ry)
  real   (kind=8) :: dbasisx(0:px,0:nderiv,0:rx)
  real   (kind=8) :: dbasisy(0:py,0:nderiv,0:ry)
  ! Rdbasis(0) => Rij
  ! Rdbasis(1) => dx Rij
  ! Rdbasis(2) => dy Rij
  real   (kind=8) :: Rdbasis(0:N)  
  real   (kind=8) :: C(d,0:N)
  real   (kind=8) :: weight 
  real   (kind=8) :: M, Mx, My, Mxy, Mxx, Myy
  real   (kind=8) :: w, wx, wy, wxy, wxx, wyy

  Cw = 0.0

  !
  do ix = 0, rx
     spanx(ix) = FindSpan(nx,px,X(ix),Ux)
     call DersBasisFuns(spanx(ix),X(ix),px,nderiv,Ux,dbasisx(:,0:nderiv,ix))
  end do
  do iy = 0, ry
     spany(iy) = FindSpan(ny,py,Y(iy),Uy)
     call DersBasisFuns(spany(iy),Y(iy),py,nderiv,Uy,dbasisy(:,0:nderiv,iy))
  end do

  !
  ! compute 
  ! w   = sum wij Ni   Nj
  ! wx  = sum wij Ni'  Nj
  ! wy  = sum wij Ni   Nj'
  ! wxx = sum wij Ni'' Nj
  ! wxy = sum wij Ni'  Nj'
  ! wyy = sum wij Ni   Nj''
  do ix = 0, rx
  ox = spanx(ix) - px
  do iy = 0, ry
     oy = spany(iy) - py

     ! --- compute w and its Derivatives
     w   = 0.0 ; wx  = 0.0 ; wy  = 0.0
     wxx = 0.0 ; wxy = 0.0 ; wyy = 0.0
     do jx = 0, px
     do jy = 0, py
        weight = Pw(d,oy+jy,ox+jx)

        M   = dbasisx(jx,0,ix) * dbasisy(jy,0,iy)
        w  = w  + M   * weight

        if (nderiv >= 1) then
           Mx  = dbasisx(jx,1,ix) * dbasisy(jy,0,iy)
           My  = dbasisx(jx,0,ix) * dbasisy(jy,1,iy)
           wx = wx + Mx  * weight
           wy = wy + My  * weight
        end if

        if (nderiv >= 2) then
           Mxx = dbasisx(jx,2,ix) * dbasisy(jy,0,iy)
           Mxy = dbasisx(jx,1,ix) * dbasisy(jy,1,iy)
           Myy = dbasisx(jx,0,ix) * dbasisy(jy,2,iy)

           wxx = wxx + Mxx * weight
           wxy = wxy + Mxy * weight
           wyy = wyy + Myy * weight
        end if 
     end do
     end do
     ! ---

     ! compute Nurbs and their derivatives
        C = 0.0     
        do jx = 0, px
        do jy = 0, py     
           M   = dbasisx(jx,0,ix) * dbasisy(jy,0,iy)
           Rdbasis(0) = M  / w

           if (nderiv >= 1) then
              Mx  = dbasisx(jx,1,ix) * dbasisy(jy,0,iy)
              My  = dbasisx(jx,0,ix) * dbasisy(jy,1,iy)
              if (rationalize==1) then              
                 Rdbasis(1) = Mx / w - M * wx / w**2 
                 Rdbasis(2) = My / w - M * wy / w**2 
              else
                 Rdbasis(1) = Mx
                 Rdbasis(2) = My
              end if
           end if
        
           if (nderiv >= 2) then
              Mxx = dbasisx(jx,2,ix) * dbasisy(jy,0,iy)
              Mxy = dbasisx(jx,1,ix) * dbasisy(jy,1,iy)
              Myy = dbasisx(jx,0,ix) * dbasisy(jy,2,iy)
              if (rationalize==1) then              
                 Rdbasis(3) = Mxx / w                 &
                            - 2 * Mx * wx / w**2      &
                            - M * wxx / w**2          &
                            + 2 * M * wx**2 / w**3
              
                 Rdbasis(4) = Mxy / w                 &
                            - Mx * wy / w**2          &
                            - My * wx / w**2          &
                            - M * wxy / w**2          &
                            + 2 * M * wx * wy / w**3
              
                 Rdbasis(5) = Myy / w                 &
                            - 2 * My * wy / w**2      &
                            - M * wyy / w**2          &
                            + 2 * M * wy**2 / w**3
              else
                 Rdbasis(3) = Mxx
                 Rdbasis(4) = Mxy
                 Rdbasis(5) = Myy
              end if                         
           end if

           do deriv=0,N
              C(:,deriv) = C(:,deriv) + Rdbasis(deriv) * Pw(:,oy+jy,ox+jx)
           end do
        end do
        end do

        Cw(1:d,iy,ix,0:N) = C(1:d,0:N)
     ! ---

  end do
  end do
  !
end subroutine EvaluateDeriv2
! .......................................................

! .......................................................
subroutine Evaluate3(d,nx,px,Ux,ny,py,Uy,nz,pz,Uz,Pw,rx,X,ry,Y,rz,Z,Cw)
!
!        Input: 
!        Output: 
  use bspline
  implicit none
  integer(kind=4), intent(in)  :: d
  integer(kind=4), intent(in)  :: nx, ny, nz
  integer(kind=4), intent(in)  :: px, py, pz
  integer(kind=4), intent(in)  :: rx, ry, rz
  real   (kind=8), intent(in)  :: Ux(0:nx+px+1)
  real   (kind=8), intent(in)  :: Uy(0:ny+py+1)
  real   (kind=8), intent(in)  :: Uz(0:nz+pz+1)
  real   (kind=8), intent(in)  :: Pw(d,0:nz,0:ny,0:nx)
  real   (kind=8), intent(in)  :: X(0:rx), Y(0:ry), Z(0:rz)
  real   (kind=8), intent(out) :: Cw(d,0:rz,0:ry,0:rx)
  integer(kind=4) :: ix, jx, iy, jy, iz, jz, ox, oy, oz
  integer(kind=4) :: spanx(0:rx), spany(0:ry), spanz(0:rz)
  real   (kind=8) :: basisx(0:px,0:rx), basisy(0:py,0:ry), basisz(0:pz,0:rz)
  real   (kind=8) :: M, C(d)
  !
  do ix = 0, rx
     spanx(ix) = FindSpan(nx,px,X(ix),Ux)
     call BasisFuns(spanx(ix),X(ix),px,Ux,basisx(:,ix))
  end do
  do iy = 0, ry
     spany(iy) = FindSpan(ny,py,Y(iy),Uy)
     call BasisFuns(spany(iy),Y(iy),py,Uy,basisy(:,iy))
  end do
  do iz = 0, rz
     spanz(iz) = FindSpan(nz,pz,Z(iz),Uz)
     call BasisFuns(spanz(iz),Z(iz),pz,Uz,basisz(:,iz))
  end do
  !
  do ix = 0, rx
     ox = spanx(ix) - px
     do iy = 0, ry
        oy = spany(iy) - py
        do iz = 0, rz
           oz = spanz(iz) - pz
           ! ---
           C = 0.0
           do jx = 0, px
              do jy = 0, py
                 do jz = 0, pz
                    M = basisx(jx,ix) * basisy(jy,iy) * basisz(jz,iz)
                    C = C + M * Pw(:,oz+jz,oy+jy,ox+jx)
                 end do
              end do
           end do
           Cw(:,iz,iy,ix) = C
           ! ---
        end do
     end do
  end do
  !
end subroutine Evaluate3
! .......................................................


end module BSp

!
! ----------
!
