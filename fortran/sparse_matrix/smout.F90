!*********************************************************************************
!>  This file contains some routines used for output operations, 
!>  on sparsematrix module
!*********************************************************************************

!---------------------------------------------------------------------------------
	subroutine printProfil_csrMatrix(this, as_file, ai_typeprint)
	implicit none
		type(csr_matrix) :: this 
		character(len=*), intent(in) :: as_file		
        integer, optional :: ai_typeprint
		!local var
		real(wp), dimension(:,:), pointer :: lpr_temp
		integer  :: li_err,li_flag
		integer  :: li_i,li_j,li_k	
		integer  :: li_file=1
		integer  :: li_ios
		integer, parameter :: li_MAX = 3000
        integer  :: li_typeprint
        
        ! element per line
        li_typeprint = 0
        if ( present ( ai_typeprint ) ) then 
            li_typeprint = ai_typeprint
        end if
		
		if ( ( this%oi_nR > li_MAX ) .OR. ( this%oi_nC > li_MAX ) ) then
			print*,'error printProfil_csrMatrix: Matrix size must be less than li_MAX'
			return
		end if
		
		allocate(lpr_temp(this%oi_nR,this%oi_nC),stat=li_err)
		if (li_err.ne.0) li_flag=10	
		
		!initialisation
		lpr_temp(1:this%oi_nR,1:this%oi_nC) = 0.0_wp
						
		do li_i =1,this%oi_nR 
			do li_k = this%opi_ia(li_i),this%opi_ia(li_i+1)-1
				li_j = this%opi_ja(li_k)
				lpr_temp(li_i,li_j) = this%opr_a(li_k)		
			end do
		end do	
		
		open(unit=li_file, file='../Runs/'//as_file//'.dat', iostat=li_ios)
		if (li_ios/=0) STOP "erreur d'ouverture du fichier new.dat"
        if ( li_typeprint == 0 ) then
            do li_i =1,this%oi_nR 
                do li_j =1,this%oi_nC 					
                    write(li_file, *)lpr_temp(li_i,li_j)
                end do
            end do
        end if
        if ( li_typeprint == 1 ) then
            do li_i =1,this%oi_nR 
                write(li_file, *)lpr_temp(li_i,1:this%oi_nC)
            end do
        end if
        
		close(li_file)
				
!		call profilMatrix ( lpr_temp, this%oi_nR, this%oi_nC )
!		call printProfilMatrixInFile ( lpr_temp, this%oi_nR, this%oi_nC )
		
		deallocate(lpr_temp)
	end subroutine printProfil_csrMatrix
!---------------------------------------------------------------------------------
	subroutine todense_csrMatrix(this, apr_dense)
	implicit none
		type(csr_matrix) :: this 
		real(wp), dimension(:,:) :: apr_dense		
		!local var
		integer  :: li_i,li_j,li_k	        
		
		!initialisation
		apr_dense ( 1 : this%oi_nR , 1 : this%oi_nC ) = 0.0_wp
						
		do li_i =1,this%oi_nR 
		
			do li_k = this%opi_ia ( li_i ) , this%opi_ia ( li_i + 1 ) - 1
			
				li_j = this%opi_ja(li_k)
				
				apr_dense ( li_i , li_j ) = this%opr_a ( li_k )		
				
			end do
			
		end do	
				
	end subroutine todense_csrMatrix	
!---------------------------------------------------------------------------------
!%%MatrixMarket matrix coordinate real general
!%=================================================================================
!%
!% This ASCII file represents a sparse MxN matrix with L 
!% nonzeros in the following Matrix Market format:
!%
!% +----------------------------------------------+
!% |%%MatrixMarket matrix coordinate real general | <--- header line
!% |%                                             | <--+
!% |% comments                                    |    |-- 0 or more comment lines
!% |%                                             | <--+         
!% |    M  N  L                                   | <--- rows, columns, entries
!% |    I1  J1  A(I1, J1)                         | <--+
!% |    I2  J2  A(I2, J2)                         |    |
!% |    I3  J3  A(I3, J3)                         |    |-- L lines
!% |        . . .                                 |    |
!% |    IL JL  A(IL, JL)                          | <--+
!% +----------------------------------------------+   
!%
!% Indices are 1-based, i.e. A(1,1) is the first element.
!%
!%=================================================================================
	subroutine print_csrMatrix_mmFormat(this, as_file)
	implicit none
		type(csr_matrix) :: this 
		character(len=*), intent(in) :: as_file		
		!local var
		integer  :: li_i,li_j,li_k	
		integer  :: li_file=1
		integer  :: li_ios

		open(unit=li_file, file='../Runs/'//as_file//'.mm', iostat=li_ios)
		if (li_ios/=0) STOP "print_nnz_csrMatrix : erreur d'ouverture du fichier "
						
		write(li_file, *)'%%MatrixMarket matrix coordinate real general'
		write(li_file, *)this%oi_nR,',',this%oi_nC,',',this%oi_nel
		do li_i =1,this%oi_nR 

			do li_k = this%opi_ia(li_i),this%opi_ia(li_i+1)-1

				li_j = this%opi_ja(li_k)
				
				write(li_file, *)li_i,',',li_j,',',this%opr_a(li_k)		

			end do

		end do	
		
		close(li_file)
				
	end subroutine print_csrMatrix_mmFormat
!----------------------------------------------------------! 
	subroutine profilMatrix(apr_A,ai_nR, ai_nC)
	implicit none
		integer  :: ai_nR, ai_nC
		real(wp), dimension(ai_nR, ai_nC),intent(in)  :: apr_A
		!local var
		integer  :: li_i, li_j
		character, dimension(ai_nR, ai_nC)  :: lpr_charA
		
		do li_i = 1, ai_nR
			do li_j = 1, ai_nC	
				if ( dabs(apr_A( li_i , li_j )) < epsilon ) then
					lpr_charA( li_i , li_j ) = '-'
				else
					lpr_charA( li_i , li_j ) = 'x'
				end if
			end do
		end do	
		
		print*,'>>>>>>>profilMatrix : Begin'
		do li_i = 1, ai_nR
			print*,lpr_charA( li_i , 1:ai_nC ) 
		end do			
		print*,'>>>>>>>profilMatrix : End'
				
	end subroutine profilMatrix		
!----------------------------------------------------------------------------------------------	
	subroutine printProfilMatrixInFile(apr_A,ai_nR, ai_nC)
	implicit none
		integer  :: ai_nR, ai_nC
		real(wp), dimension(ai_nR, ai_nC),intent(in)  :: apr_A
		!local var
		integer  :: li_ios, li_ioutput
		integer  :: li_i,li_j		
		
		li_ioutput = 1
		
		open(unit=li_ioutput, file='../Runs/profilMatrix.dat', iostat=li_ios)
		if (li_ios/=0) STOP "erreur d'ouverture du fichier profilMatrix.dat"
				
		do li_i = 1, ai_nR																			
			do li_j = 1, ai_nC
			
				if ( dabs(apr_A( li_i , li_j )) < epsilon ) then
					write(li_ioutput, *)0.0_wp
				else
					write(li_ioutput, *)apr_A( li_i , li_j )
				end if
				
			end do					
			
			write(li_ioutput, *)' '				
		end do	

		close(li_ioutput)
		
	end subroutine printProfilMatrixInFile	