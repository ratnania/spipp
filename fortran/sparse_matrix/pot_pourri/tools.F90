subroutine import_param ( ai_nC, ai_nR, ai_nel, ai_nnz, ai_nen, ai_ndof )
implicit none
	integer  :: ai_nC 
	integer  :: ai_nR 
	integer  :: ai_nel
	integer  :: ai_nnz		
	integer  :: ai_nen
	integer  :: ai_ndof
	! LOCAL VARIABLES
	integer :: li_ios, li_flag		
	integer  :: li_param
	character(41)	:: ls_currentDir	
	
	li_param = 1
	ls_currentDir = '/Users/ahmedratnani/Projects/SparseMatrix'	
	
	!OPEN THE DATA FILE param.dat
	open(unit=li_param, file=ls_currentDir//'/data/param.dat',action="read", iostat=li_ios)
	if (li_ios/=0) STOP "erreur d'ouverture du fichier des donnees param.dat"
	
	call lire(li_param,ai_nC,"ai_nC")	
	call lire(li_param,ai_nR,"ai_nR")	
	call lire(li_param,ai_nel,"ai_nel")					
	call lire(li_param,ai_nnz,"ai_nnz")	
	call lire(li_param,ai_nen,"ai_nen")	
	call lire(li_param,ai_ndof,"ai_ndof")			
	
	!CLOSE THE DATA FILE param.dat
	close(li_param)	
end subroutine import_param
!-------------------------------------------------------------------------------------------  
subroutine import_Matrix ( ao_A, ai_nel, ai_nen, ai_ndof )
implicit none
	type(csr_Matrix) :: ao_A
	integer  :: ai_nel, ai_nen, ai_ndof	
	! LOCAL VARIABLES
	integer  :: li_ios, li_flag, li_err		
	integer  :: li_i, li_input, li_e, li_b
	integer  :: li_Matrix
	character(41)	:: ls_currentDir	
	
	li_Matrix = 1	
	ls_currentDir = '/Users/ahmedratnani/Projects/SparseMatrix'	
	
	!OPEN THE DATA FILE Matrix.dat
	open(unit=li_Matrix, file=ls_currentDir//'/data/Matrix.dat',action="read", iostat=li_ios)
	if (li_ios/=0) STOP "erreur d'ouverture du fichier des donnees Matrix.dat"
	
	! READING ia
	do li_i = 1, ao_A%oi_nR+1
		call lire(li_Matrix,li_input,"index")
		ao_A%opi_ia ( li_i ) = li_input
	end do

	! READING ja
	do li_i = 1, ao_A%oi_nel
		call lire(li_Matrix,li_input,"index")
		ao_A%opi_ja ( li_i ) = li_input
	end do
	
	close(li_Matrix)	
	
end subroutine import_Matrix
!-------------------------------------------------------------------------------------------  
subroutine import_LM ( api_LM, ai_nel, ai_nen )
implicit none
	integer, dimension(:,:) :: api_LM
	integer  :: ai_nel, ai_nen	
	! LOCAL VARIABLES
	integer  :: li_ios, li_flag, li_err		
	integer  :: li_input, li_e, li_b
	integer  :: li_LM
	character(41)	:: ls_currentDir	
	
	li_LM = 2	
	ls_currentDir = '/Users/ahmedratnani/Projects/SparseMatrix'	
			
	!OPEN THE DATA FILE LM.dat
	open(unit=li_LM, file=ls_currentDir//'/data/LM.dat',action="read", iostat=li_ios)
	if (li_ios/=0) STOP "erreur d'ouverture du fichier des donnees LM.dat"
	
	! READING LM
	do li_e = 1, ai_nel	
		do li_b = 1, ai_nen
			call lire(li_LM,li_input,"index")
			api_LM ( li_b, li_e ) = li_input
		end do						
	end do
	
	close(li_LM)
	
end subroutine import_LM
!-------------------------------------------------------------------------------------------  
subroutine import_LM_MULT ( api_LM, ai_nel, ai_nen, ai_ndof )
implicit none
	integer, dimension(:,:,:) :: api_LM
	integer  :: ai_nel, ai_nen, ai_ndof	
	! LOCAL VARIABLES
	integer  :: li_ios, li_flag, li_err		
	integer  :: li_i, li_input, li_e, li_b
	integer  :: li_LM
	character(41)	:: ls_currentDir	
	
	li_LM = 2	
	ls_currentDir = '/Users/ahmedratnani/Projects/SparseMatrix'	
			
	!OPEN THE DATA FILE LM.dat
	open(unit=li_LM, file=ls_currentDir//'/data/LM.dat',action="read", iostat=li_ios)
	if (li_ios/=0) STOP "erreur d'ouverture du fichier des donnees LM.dat"
	
	! READING LM
	do li_e = 1, ai_nel	
		do li_b = 1, ai_nen
			do li_i = 1, ai_ndof
				call lire(li_LM,li_input,"index")
				api_LM ( li_i, li_b, li_e ) = li_input
			end do
		end do						
	end do
	
	close(li_LM)
	
end subroutine import_LM_MULT