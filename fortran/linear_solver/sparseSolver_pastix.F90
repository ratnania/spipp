module sparseSolver_module
use used_precision
use SparseMatrix_Module
implicit none

	integer, parameter  :: mi_maxn		= 100000
	integer, parameter  :: mi_maxnz		= 3000000
	integer, parameter  :: mi_maxnrhs	= 1		

#include "pastix_fortran.h"
!#ifndef FORCE_NOMPI
	include 'mpif.h'
!#endif

!	module pastix_toolbox
		pastix_float_t    avals(mi_maxnz)
		pastix_float_t    b(mi_maxn)
		real*8    dparm(64)

		pastix_int_t   ia(mi_maxn)
		pastix_int_t   ja(mi_maxnz)		
		pastix_int_t   iparm(64)
!	end module pastix_toolbox

		pastix_int_t n
		pastix_int_t nrhs		
		pastix_int_t nnz
		INTEGER   m
		INTEGER   ldb
		INTEGER   i
		INTEGER   ierr
		INTEGER   n_iter
		INTEGER   iter										
		integer   numprocs
		integer   numthreads
		integer   myid				
		INTEGER   required
		INTEGER   provided
		INTEGER   StatInfo				
		!     ATTENTION : l'integer pastix_data doit etre de la taille 
		!     de l'adressage ce qui depend de la machine.
		integer(kind=8) pastix_data

	private 
	public  :: main_Pastix, sparseSolver_FINALIZE, init_sparseSolver_Pastix
contains
	!*******************************************
	!***  Solving subroutine                   *
	!***  Runs PaStiX                          *
	!*******************************************
	subroutine solve( pastix_data	&
					, mpi_comm		&
					, n				&
					, myid			&
					, numthreads	&
					, iter			&
					, n_iter		&
					, ldb			)
	implicit none
		real*8            t0
		real*8            t1
		real*8            t2
		real*8            t3
		real*8            t4
		real*8            t5
		real*8            t6
		real*8            t7
			
		pastix_data_ptr_t pastix_data
		pastix_float_t    b_save(mi_maxn)

		integer mpi_comm
		pastix_int_t perm(mi_maxn)
		pastix_int_t invp(mi_maxn)		
		pastix_int_t n
		pastix_int_t nrhs		

		integer myid
		integer numthreads
		integer iter
		integer n_iter
		integer soliter
		integer n_soliter

		INTEGER   ldb
		
		ldb = mi_maxn
		nrhs= 1

	!***********************************************************************
	!     Initializing PaStiX
	!     
	!     Calling PaStiX with iparm(IPARM_MODIFY_PARAMETER) set to 0 will
	!     fill iparm and dparm with default parameters and return.
	!***********************************************************************

		call cpu_time(t0)

		if(iter==1) then

		iparm(IPARM_MODIFY_PARAMETER) = 0

!		write(*,*) 'FORTRAN : call PastiX'     
		pastix_data = 0;
		nrhs = 1;

		call pastix_fortran(pastix_data ,mpi_comm,n,ia,ja,avals,perm,&
		invp,b,nrhs,iparm,dparm)
		end if  

		call cpu_time(t1)

!		if (myid .eq. 0) then
!		print *,'FORTRAN : PaStiX init in time - ',t1-t0
!		if (iparm(64) .ne. 0) then
!		print *,'FORTRAN : The following ERROR was detected: ',&
!		iparm(64)
!		stop
!		end if
!		endif

	!***********************************************************************
	!     Analysis
	!***********************************************************************

		call cpu_time(t2)

		iparm(IPARM_VERBOSE)			= API_VERBOSE_NOT
		iparm(IPARM_START_TASK)         = API_TASK_ORDERING
		iparm(IPARM_END_TASK)           = API_TASK_ANALYSE  ! Do all preprocessing tasks
		iparm(IPARM_ITERMAX)            = 1               ! refinement : max number of iterations
		iparm(IPARM_FACTORIZATION)      = API_FACT_LU     ! LU, LLT or LDLT
		iparm(IPARM_THREAD_NBR)         = numthreads      ! number of threads
		iparm(IPARM_RHS_MAKING)         = API_RHS_B       ! right hand side (B : use given, 1 : RHS = 1, I : RHS(i) = i)
		iparm(IPARM_LEVEL_OF_FILL)      = 1
		iparm(IPARM_SYM)                = API_SYM_NO

		dparm(DPARM_EPSILON_REFINEMENT) = 1.e-12          ! error level refinement
		dparm(DPARM_EPSILON_MAGN_CTRL)  = 1.e-32          ! pivot threshold


		if(iter==1) then 

		call pastix_fortran(pastix_data, mpi_comm,n,ia,ja,avals, &
		perm,invp,b,nrhs,iparm,dparm)
		call cpu_time(t3)

!		if (myid .eq. 0) then
!		print *,'FORTRAN : Analysis complete in time - ',t3-t2
!		endif

		if (iparm(64) .ne. 0) then
		print *,'FORTRAN : The following ERROR was detected: ',&
		iparm(64)
		stop
		end if

		end if

	!***********************************************************************
	!     Factorization
	!***********************************************************************

		call cpu_time(t4)

		iparm(IPARM_START_TASK) = API_TASK_NUMFACT
		iparm(IPARM_END_TASK)   = API_TASK_NUMFACT

		call pastix_fortran(pastix_data,mpi_comm,n,ia,ja,avals,perm,&
		invp,b,nrhs,iparm,dparm)

		call cpu_time(t5)

!		if (myid .eq. 0) then
!		print *,'FORTRAN : Factorisation complete in time - ',t5-t4
!		endif

		if (iparm(64) .ne. 0) then
		print *,'FORTRAN : The following ERROR was detected: ',iparm(64)
		stop
		end if

!		if (myid .eq. 0) then
!		print *,'FORTRAN : Number of expected nonzeros in factors = ',&
!		iparm(23)
!		print *,&
!		'FORTRAN : Number of expected FLOPS in factorization = ',&
!		dparm(23)
!		endif

	!***********************************************************************
	!     Solve
	!***********************************************************************

		n_soliter = 1
		b_save = b
		call cpu_time(t6)
		do soliter = 1, n_soliter

!		write(*,*) '****************************************'
!		write(*,*) ' Solve Iteration n ', soliter
!		write(*,*) '****************************************'

		!     call cpu_time(t6)

		iparm(IPARM_START_TASK) = API_TASK_SOLVE
		iparm(IPARM_END_TASK)   = API_TASK_REFINE

		b = b_save
		call pastix_fortran(pastix_data, mpi_comm,n,ia,ja,avals,&
		perm,invp,b,nrhs,iparm,dparm)

		!     call cpu_time(t7)

		end do                    !(iteration loop)
		call cpu_time(t7)

!		if (myid .eq. 0) then
!		print *,'FORTRAN : ', n_soliter, &
!		'Solve complete in time - ',t7-t6
!		endif

		if (iparm(64) .ne. 0) then
		print *,'FORTRAN : The following ERROR was detected: ',&
		iparm(IPARM_ERROR_NUMBER)
		stop
		end if


	!***********************************************************************
	!     Fin - cleaning
	!***********************************************************************

		iparm(IPARM_START_TASK) = API_TASK_CLEAN
		iparm(IPARM_END_TASK)   = API_TASK_CLEAN 

		if(iter==n_iter) then 

		call pastix_fortran(pastix_data, mpi_comm,n,ia,ja,avals,perm,&
		invp,b,nrhs,iparm,dparm)
		end if

	end subroutine solve
!----------------------------------------------------------------------------------------------		
	subroutine init_sparseSolver_Pastix ( )
		implicit none
#ifndef FORCE_NOMPI
	include 'mpif.h'
#else
#define MPI_COMM_WORLD 0
#endif

		numthreads=1

!     Initialize MPI environment 
!     call mpi_init(ierr)
#ifndef FORCE_NOMPI
		required=MPI_THREAD_MULTIPLE;
		call MPI_Init_thread(required,provided,StatInfo);

		call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs,ierr)
		call MPI_COMM_RANK(MPI_COMM_WORLD, myid,ierr)
#endif

		!***********************************************************************	  		
		ldb = mi_maxn
		nrhs= 1

!if (myid .ne. 0)  n   = 0

		n_iter = 1

		!***********************************************************************
		! Start iterations
		!***********************************************************************
		pastix_data = 0;
		
	end subroutine init_sparseSolver_Pastix
!----------------------------------------------------------------------------------------------		
	subroutine main_Pastix ( this, apr_rhs, apr_x )
		implicit none
#ifndef FORCE_NOMPI
	include 'mpif.h'
#else
#define MPI_COMM_WORLD 0
#endif
	! -----------------------------
	! INPUTS
		type(csr_matrix) :: this
		real*8, dimension ( : ) :: apr_rhs		
	! -----------------------------
	! OUTPUTS
		real*8, dimension ( : ) :: apr_x		
	! -----------------------------		
!		! LOCAL VARIABLES
!		pastix_int_t n
!		pastix_int_t nrhs		
!		pastix_int_t nnz
!		INTEGER   m
!		INTEGER   ldb
!		INTEGER   i
!		INTEGER   ierr
!		INTEGER   n_iter
!		INTEGER   iter										
!		integer   numprocs
!		integer   numthreads
!		integer   myid				
!		INTEGER   required
!		INTEGER   provided
!		INTEGER   StatInfo				
!		!     ATTENTION : l'integer pastix_data doit etre de la taille 
!		!     de l'adressage ce qui depend de la machine.
!		integer(kind=8) pastix_data
!		!     pour driver mmread ...
!		!     character rep*10
!		!     character field*7
!		!     character symm*19
!		integer, parameter :: li_file = 1 
!		integer  :: li_i, li_ios
!
!		numthreads=1
!
!!     Initialize MPI environment 
!!     call mpi_init(ierr)
!#ifndef FORCE_NOMPI
!		required=MPI_THREAD_MULTIPLE;
!		call MPI_Init_thread(required,provided,StatInfo);
!
!		call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs,ierr)
!		call MPI_COMM_RANK(MPI_COMM_WORLD, myid,ierr)
!#endif
!
!		!***********************************************************************	  		
!		ldb = mi_maxn
!		nrhs= 1
!
!!if (myid .ne. 0)  n   = 0
!
!		n_iter = 1
!
!		!***********************************************************************
!		! Start iterations
!		!***********************************************************************
!		pastix_data = 0;
		do iter = 1, n_iter
         
!			write(*,*) '****************************************'
!			write(*,*) ' Iteration n ', iter
!			write(*,*) '****************************************'

			!***********************************************************************
			! Initialize vectors
			!***********************************************************************

			m     = this%oi_nC
			n     = this%oi_nR	
			nnz   = this%oi_nel	
			ia	  ( 1 : this%oi_nC + 1 ) = this%opi_ia	( 1 : this%oi_nC + 1 )	
			ja    ( 1 : this%oi_nel )    = this%opi_ja	( 1 : this%oi_nel )
			avals ( 1 : this%oi_nel )    = this%opr_a	( 1 : this%oi_nel )
			b	  ( 1 : this%oi_nR )     = apr_rhs		( 1 : this%oi_nR )	
						
			!***********************************************************************
			! Solver call
			!***********************************************************************

			call solve(pastix_data,MPI_COMM_WORLD,n,myid,numthreads,iter,n_iter,ldb)

		end do                    !(iteration loop)

		! Copy the solution 
		apr_x ( 1 : this%oi_nR ) = b ( 1 : this%oi_nR )
		
	end subroutine main_Pastix
!----------------------------------------------------------------------------------------------			
	subroutine sparseSolver_FINALIZE (IERR) 
	implicit none
		integer  :: IERR
	!***********************************************************************
	! End
	!***********************************************************************
#ifndef FORCE_NOMPI
		call MPI_FINALIZE(IERR)   ! clean up MPI
#endif
	end subroutine sparseSolver_FINALIZE

end module sparseSolver_module