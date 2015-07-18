!
!	sparseSolver.f90
!	NURBS_FEM
!
!	Created by ahmed ratnani on 26/01/09.
!	Copyright 2009 __MyCompanyName__. All rights reserved.
!


module sparseSolver_module
use used_precision
use SparseMatrix_Module
implicit none
	integer, parameter  :: oi_detail = 0
	type, public :: paramPardiso2D
!..     Internal solver memory pointer for 64-bit architectures
!..     INTEGER*8 pt(64)
!..     Internal solver memory pointer for 32-bit architectures
!..     INTEGER*4 pt(64)
!..     This is OK in both cases
        integer*4, dimension(64) :: opi_pt!
		integer*4, dimension(64) :: opi_iparm
		
        integer :: oi_maxfct !
		integer :: oi_mnum !
		integer :: oi_mtype
		integer :: oi_phase
		integer :: oi_n !
		integer :: oi_nrhs !
		integer :: oi_error
		integer :: oi_msglvl
		integer :: oi_idum
		real*8  :: or_ddum
				
		real*8, dimension(:), pointer  :: opr_a		
		integer, dimension(:), pointer :: opi_ia
        integer, dimension(:), pointer :: opi_ja
		real*8, dimension(:), pointer  :: opr_x	
		real*8, dimension(:), pointer  :: opr_b		
			
	end type paramPardiso2D
	contains
!------------------------------------------------------------------------	
	subroutine initSparseSolver(this, ai_dim, ao_Mat , apr_b, apr_x)
	implicit none
		type(paramPardiso2D) :: this
		type(csr_matrix) :: ao_Mat
		real(wp), dimension(:)  :: apr_b		
		real(wp), dimension(:):: apr_x	
		integer  :: ai_dim		
		!local var
		integer  :: li_i, li_j,li_flag,li_err	,li_k
		integer  :: li_file=1
		integer  :: li_ios			
		
		this%oi_maxfct = 1
		this%oi_mnum   = 1	
		this%oi_nrhs   = 1
		this%oi_n      = ai_dim
		
		do li_i = 1, 64
			this%opi_pt(li_i) = 0
			this%opi_iparm(li_i) = 0
		end do

		this%opi_iparm(1) = 1 ! no solver default
		this%opi_iparm(2) = 2 ! fill-in reordering from METIS
		this%opi_iparm(3) = 2 ! numbers of processors
		this%opi_iparm(4) = 0 ! no iterative-direct algorithm
		this%opi_iparm(5) = 0 ! no user fill-in reducing permutation
		this%opi_iparm(6) = 0 ! solution on the first n compoments of x
		this%opi_iparm(7) = 0 ! not in use
		this%opi_iparm(8) = 5 ! numbers of iterative refinement steps
		this%opi_iparm(9) = 0 ! not in use
		this%opi_iparm(10) = 13 ! perturbe the pivot elements with 1E-13
		this%opi_iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
		this%opi_iparm(12) = 0 ! not in use
		this%opi_iparm(13) = 1 ! maximum weighted matching algorithm is switched-on
		this%opi_iparm(14) = 0 ! Output: number of perturbed pivots
		this%opi_iparm(15) = 0 ! not in use
		this%opi_iparm(16) = 0 ! not in use
		this%opi_iparm(17) = 0 ! not in use
		this%opi_iparm(18) = -1 ! Output: number of nonzeros in the factor LU
		this%opi_iparm(19) = -1 ! Output: Mflops for LU factorization
		this%opi_iparm(20) = 0 ! Output: Numbers of CG Iterations
		
		!on cree l'objet correspondant, mais cette fois pas avec la structure ijaSparseMatrix2D

		!allocation  
		allocate(this%opi_ia(ai_dim+1),stat=li_flag)
		if (li_err.ne.0) li_flag=10		
		allocate(this%opi_ja(ao_Mat%oi_nel),stat=li_flag)
		if (li_err.ne.0) li_flag=20
		allocate(this%opr_a(ao_Mat%oi_nel),stat=li_flag)
		if (li_err.ne.0) li_flag=30
		allocate(this%opr_x(ai_dim),stat=li_flag)
		if (li_err.ne.0) li_flag=40	
		allocate(this%opr_b(ai_dim),stat=li_flag)
		if (li_err.ne.0) li_flag=50	
		
		this%opi_ia ( : ) = ao_Mat%opi_ia ( : )
		this%opi_ja ( : ) = ao_Mat%opi_ja ( : )
		this%opr_a ( : ) = ao_Mat%opr_a ( : )
												
		do li_i = 1, ai_dim
			this%opr_b(li_i) = real(apr_b(li_i))
		end do			
			
	end subroutine initSparseSolver
!------------------------------------------------------------------------	
	subroutine freeThis(this)
	implicit none
		type(paramPardiso2D) :: this
		deallocate(this%opi_ia)	
		deallocate(this%opi_ja)
		deallocate(this%opr_a)				
		deallocate(this%opr_x)		
		deallocate(this%opr_b)		
	end subroutine freeThis
!------------------------------------------------------------------------	
	subroutine solvePardiso2D(ai_dim, ao_Mat , apr_b, apr_x)
	implicit none
		type(csr_matrix) :: ao_Mat
		real(wp), dimension(:)  :: apr_b		
		real(wp), dimension(:):: apr_x	
		integer  :: ai_dim	
		!local var
		type(paramPardiso2D) :: lo_this
		integer :: li_i
		
		!initialisation des tableaux
		call initSparseSolver(lo_this, ai_dim, ao_Mat , apr_b, apr_x)
!
!  .. Setup Pardiso control parameters und initialize the solvers
!     internal adress pointers. This is only necessary for the FIRST
!     call of the PARDISO solver.
!
      lo_this%oi_mtype     = 11  ! unsymmetric matrix symmetric, indefinite
      lo_this%oi_msglvl    = 0       !0-> no screen output ;  1-> with statistical information	
	  lo_this%oi_error = 0

!..   Reordering and Symbolic Factorization, This step also allocates
!     all memory that is necessary for the factorization

      lo_this%oi_phase     = 11     ! only reordering and symbolic factorization
	  
	  CALL pardiso (lo_this%opi_pt, lo_this%oi_maxfct, lo_this%oi_mnum, lo_this%oi_mtype, &
					lo_this%oi_phase, lo_this%oi_n, &
					lo_this%opr_a, lo_this%opi_ia, lo_this%opi_ja, lo_this%oi_idum, lo_this%oi_nrhs, &
					lo_this%opi_iparm, lo_this%oi_msglvl, lo_this%or_ddum, lo_this%or_ddum, lo_this%oi_error)

	  if (lo_this%oi_msglvl == 1) then
		WRITE(*,*) 'Reordering completed ... '
	  end if

      IF (lo_this%oi_error .NE. 0) THEN
        WRITE(*,*) 'The following ERROR was detected: ', lo_this%oi_error
        STOP
      END IF
	if (lo_this%oi_msglvl == 1) then
      WRITE(*,*) 'Number of nonzeros in factors   = ',lo_this%opi_iparm(18)
      WRITE(*,*) 'Number of factorization MFLOPS  = ',lo_this%opi_iparm(19)
	end if  
	  
!.. Factorization.
      lo_this%oi_phase     = 22  ! only factorization
      CALL pardiso (lo_this%opi_pt, lo_this%oi_maxfct, lo_this%oi_mnum, lo_this%oi_mtype, &
					lo_this%oi_phase, lo_this%oi_n, lo_this%opr_a, lo_this%opi_ia, lo_this%opi_ja, &
                    lo_this%oi_idum, lo_this%oi_nrhs, lo_this%opi_iparm, lo_this%oi_msglvl, &
					lo_this%or_ddum, lo_this%or_ddum, lo_this%oi_error)
	if (lo_this%oi_msglvl == 1) then
      WRITE(*,*) 'Factorization completed ... '
	end if
      IF (lo_this%oi_error .NE. 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', lo_this%oi_error
        STOP
      ENDIF

!.. Back substitution and iterative refinement
!      lo_this%opi_iparm(8)  = 2   ! max numbers of iterative refinement steps
      lo_this%oi_phase     = 33  ! only factorization
	  
      CALL pardiso (lo_this%opi_pt, lo_this%oi_maxfct, lo_this%oi_mnum, lo_this%oi_mtype, &
					lo_this%oi_phase, lo_this%oi_n, lo_this%opr_a, lo_this%opi_ia, lo_this%opi_ja, &
                    lo_this%oi_idum, lo_this%oi_nrhs, lo_this%opi_iparm, lo_this%oi_msglvl, &
					lo_this%opr_b, lo_this%opr_x, lo_this%oi_error)

	  if (lo_this%oi_msglvl == 1) then
		WRITE(*,*) 'Solve completed ...  '
	  end if

	  !copie de la solution
      DO li_i = 1, lo_this%oi_n
		apr_x(li_i) = 1.0_wp * lo_this%opr_x(li_i)
      END DO

!.. Termination and release of memory
      lo_this%oi_phase     = -1           ! release internal memory
      CALL pardiso (lo_this%opi_pt, lo_this%oi_maxfct, lo_this%oi_mnum, lo_this%oi_mtype, &
					lo_this%oi_phase, lo_this%oi_n, lo_this%or_ddum, lo_this%oi_idum, lo_this%oi_idum, &
                    lo_this%oi_idum, lo_this%oi_nrhs, lo_this%opi_iparm, lo_this%oi_msglvl, &
					lo_this%or_ddum, lo_this%or_ddum, lo_this%oi_error)

	  call freeThis(lo_this)
	end subroutine solvePardiso2D
end module sparseSolver_module




