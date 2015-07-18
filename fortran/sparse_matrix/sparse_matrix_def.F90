MODULE SPI_SPARSE_MATRIX 
USE SPI_GLOBAL_DEF 
Implicit None

!	integer, parameter  :: BC_NO_BC						= 0
!	integer, parameter  :: BC_XI_PER_ETA_NO_BC			= 1	
!	integer, parameter  :: BC_XI_MAXPER_ETA_NO_BC		= -1		
!	integer, parameter  :: BC_XI_NO_BC_ETA_PER			= 2	
!	integer, parameter  :: BC_XI_NO_BC_ETA_MAXPER		= -2	
!	integer, parameter  :: BC_DIRICHLET_HOMOGEN			= 10 	
!	integer, parameter  :: BC_ALL_PERIODIC				= 3	
!	
!	integer, parameter  :: BC_XI_DIRICHLET_ETA_NO_BC	= 11	
!	integer, parameter  :: BC_XI_NO_BC_ETA_DIRICHLET	= 12
!	integer, parameter  :: BC_XI_MAXPER_ETA_DIRICHLET	= 30	
!	integer, parameter  :: BC_XI_DIRICHLET_ETA_MAXPER	= 31	
!	integer, parameter  :: BC_XI_DIRICHLET_ETA_PER		= 21
!	integer, parameter  :: BC_XI_PER_ETA_DIRICHLET		= 22					
!
!	integer, parameter  :: BC_DIRICHLET_ROT				= 100					
!	integer, parameter  :: BC_XI_NO_BC_ETA_NO_BC_ETA_0_DIRICHLET		= 111						
!******************************************************************************************
!		SparseMatrix
!******************************************************************************************
	type csr_matrix
		integer  :: oi_nR   !NUMBER OF ROWS
		integer  :: oi_nC   !NUMBER OF COLUMNS
		integer  :: oi_nel  !NUMBER OF NON ZERO ELTS
		integer, dimension(:), pointer  :: opi_ia
		integer, dimension(:), pointer  :: opi_ja		
		real(wp), dimension(:), pointer  :: opr_a		
		!................
		logical :: ol_use_mm_format
		integer, dimension(:), pointer  :: opi_i
		!................		
	end type csr_matrix	
				
	type bnd_matrix
		integer  :: oi_n   !NUMBER OF ROWS
		integer  :: oi_lowd   
		integer  :: oi_ml   		
		integer  :: oi_mu		
		integer  :: oi_nabd  	
		integer  :: oi_LDA		
		real(wp), dimension(:,:), pointer  :: opr_abd		
		integer, dimension(:),pointer	:: opi_ipiv		
	end type bnd_matrix
!******************************************************************************************


END MODULE SPI_SPARSE_MATRIX
