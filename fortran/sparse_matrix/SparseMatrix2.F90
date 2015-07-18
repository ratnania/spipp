!---------------------------------------------------------------------------------
    !> \todo a tester
    !> create the matrix :
    !>  A , 0
    !>  0 , D
    !> INPUTS : A, D 
    !> OUTPUT : this
	subroutine create_SparseMatrixFrom2 ( this, ao_A, ao_D )
	implicit none
		!> param[inout] this : CSR MATRIX STRUCTURE
		type(csr_matrix) :: this		
		!> param[in] ao_A, ao_D : INPUT MATRICES
		type(csr_matrix) :: ao_A
		type(csr_matrix) :: ao_D        		
		!local var
		integer :: li_err,li_flag
        integer  :: li_i
        integer  :: li_j        
        integer  :: li_index
        integer  :: li_k  
        integer  :: li_iloc
        integer  :: li_jloc        
           
		this%ol_use_mm_format = .FALSE.
								     		
		this%oi_nR   = ao_A%oi_nR + ao_D%oi_nR
		this%oi_nC   = ao_A%oi_nC + ao_D%oi_nC		
		this%oi_nel  = ao_A%oi_nel  &
                     + ao_D%oi_nel
		
		allocate(this%opi_ia(this%oi_nR+1),stat=li_err)
		if (li_err.ne.0) li_flag=10	

		allocate(this%opi_ja(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=20	
				
		allocate(this%opr_a(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=30			
		
		this%opr_a ( : ) = 0.0_wp
        
        !remplissage de la matrice, a partir des matrices A,B,C et D        
!print*,'ao_A%oi_nR=',ao_A%oi_nR
        li_index = 1
        do li_i = 1, this%oi_nR
                    
            ! remplissage des termes de A
            if ( li_i <= ao_A%oi_nR ) then 

            !*********************************************************************    
            ! copie des termes de A
            !*********************************************************************            
                ! initialisation de la colonne
                li_k = ao_A%opi_ia ( li_i )
!print*,'begin'
!print*,'li_k=',li_k
                do while ( li_k < ao_A%opi_ia ( li_i + 1 ) )
!print*,'li_indexA=',li_index                
!print*,'li_k=',li_k
                    ! on recupere la colonne de A
                    li_j = ao_A%opi_ja ( li_k )
                
                    ! copie la valeur de Aij
                    this%opi_ja ( li_index ) = li_j 
                    this%opr_a  ( li_index ) = ao_A%opr_a ( li_k )                     
                    
                    ! incrementation de l'indice 
                    li_index = li_index + 1
                    li_k = li_k + 1
                
                end do
            !*********************************************************************    
!print*,'done'                            
            end if

            ! remplissage des termes de D
            if ( li_i >= ao_A%oi_nR + 1 ) then  

                li_iloc = li_i - ao_A%oi_nR            
            
            !*********************************************************************    
            ! copie des termes de D
            !*********************************************************************                            
                ! initialisation de la colonne
                li_k = ao_D%opi_ia ( li_iloc )
            
                do while ( li_k < ao_D%opi_ia ( li_iloc + 1 ) )
!print*,'li_indexD=',li_index                   
                    ! on recupere la colonne de A
                    li_j = ao_D%opi_ja ( li_k )
                    
                    li_jloc = li_j + ao_A%oi_nC                    
                
                    ! copie la valeur de Aij
                    this%opi_ja ( li_index ) = li_jloc 
                    this%opr_a  ( li_index ) = ao_D%opr_a ( li_k )                     
                    
                    ! incrementation de l'indice 
                    li_index = li_index + 1
                    li_k = li_k + 1                    
                
                end do 
                this%opi_ia ( li_i ) = li_jloc                
            !*********************************************************************                                   
!print*,'done'                
            end if
        
        end do
        
        ! mise a jour du dernier terme
        this%opi_ia ( this%oi_nR + 1 ) = li_index
!print*,'this%opi_ia=',this%opi_ia(:)        
        
	end subroutine create_SparseMatrixFrom2
!---------------------------------------------------------------------------------
    !> \todo a tester
    !> create the matrix :
    !>  A , B 
    !>  C , D
    !> INPUTS : A, B, C, D 
    !> OUTPUT : this
    !> We must have :
    !> ao_A%oi_nC = ao_C%oi_nC     
    !> ao_B%oi_nC = ao_D%oi_nC         
    !> ao_A%oi_nR = ao_B%oi_nR     
    !> ao_C%oi_nR = ao_D%oi_nR             
	subroutine create_SparseMatrixFrom4 ( this, ao_A, ao_B, ao_C, ao_D )
	implicit none
		!> param[inout] this : CSR MATRIX STRUCTURE
		type(csr_matrix) :: this		
		!> param[in] ao_A, ao_B, ao_C, ao_D : INPUT MATRICES
		type(csr_matrix) :: ao_A
		type(csr_matrix) :: ao_B
		type(csr_matrix) :: ao_C
		type(csr_matrix) :: ao_D        		
		!local var
		integer :: li_err,li_flag
        integer  :: li_i
        integer  :: li_j        
        integer  :: li_index
        integer  :: li_k  
        integer  :: li_iloc
        integer  :: li_jloc        
                
        if (    ( ao_A%oi_nC /= ao_C%oi_nC ) .OR.   &
                ( ao_B%oi_nC /= ao_D%oi_nC ) .OR.   &
                ( ao_A%oi_nR /= ao_B%oi_nR ) .OR.   &
                ( ao_C%oi_nR /= ao_D%oi_nR )        &
            ) then
                
            print*,'Error create_Sparse: ao_A, ao_B, ao_C, ao_D are not compatibles'
            STOP
            
        end if
		
		this%ol_use_mm_format = .FALSE.		
		
		this%oi_nR   = ao_A%oi_nR + ao_C%oi_nR
		this%oi_nC   = ao_A%oi_nC + ao_B%oi_nC		
		this%oi_nel  = ao_A%oi_nel  &
                     + ao_B%oi_nel  &
                     + ao_C%oi_nel  &
                     + ao_D%oi_nel
		
		allocate(this%opi_ia(this%oi_nR+1),stat=li_err)
		if (li_err.ne.0) li_flag=10	

		allocate(this%opi_ja(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=20	
				
		allocate(this%opr_a(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=30			
		
		this%opr_a ( : ) = 0.0_wp
        
        !remplissage de la matrice, a partir des matrices A,B,C et D        
print*,'ao_A%oi_nR=',ao_A%oi_nR
        li_index = 1
        do li_i = 1, this%oi_nR
                    
            ! remplissage des termes de A, puis B
            if ( li_i <= ao_A%oi_nR ) then 

            !*********************************************************************    
            ! copie des termes de A
            !*********************************************************************            
                ! initialisation de la colonne
                li_k = ao_A%opi_ia ( li_i )
print*,'-------- begin new row'
                do while ( li_k < ao_A%opi_ia ( li_i + 1 ) )
print*,'li_indexA=',li_index                
print*,'li_k=',li_k
                    ! on recupere la colonne de A
                    li_j = ao_A%opi_ja ( li_k )
                
                    ! copie la valeur de Aij
                    this%opi_ja ( li_index ) = li_j 
                    this%opr_a  ( li_index ) = ao_A%opr_a ( li_k )                     
                    
                    ! incrementation de l'indice 
                    li_index = li_index + 1
                    li_k = li_k + 1
                
                end do
            !*********************************************************************    
print*,'-------- done'                            
            !*********************************************************************    
            ! copie des termes de B
            !*********************************************************************                            
                ! initialisation de la colonne
                li_k = ao_B%opi_ia ( li_i )
            
                do while ( li_k < ao_B%opi_ia ( li_i + 1 ) )
print*,'li_indexB=',li_index                   
                    ! on recupere la colonne de A
                    li_j = ao_B%opi_ja ( li_k )

                    li_jloc = li_j + ao_A%oi_nC
                
                    ! copie la valeur de Aij
                    this%opi_ja ( li_index ) = li_jloc 
                    this%opr_a  ( li_index ) = ao_B%opr_a ( li_k )                     
                    
                    ! incrementation de l'indice 
                    li_index = li_index + 1
                    li_k = li_k + 1                    
                
                end do 
                
                this%opi_ia ( li_i ) = li_jloc                
            !*********************************************************************                                   
print*,'done'                
            end if

            ! remplissage des termes de C, puis D
            if ( li_i >= ao_A%oi_nR + 1 ) then  

                li_iloc = li_i - ao_A%oi_nR            
            
            !*********************************************************************    
            ! copie des termes de C
            !*********************************************************************            
                ! initialisation de la colonne
                li_k = ao_C%opi_ia ( li_iloc )
print*,'-------- begin new row'            
                do while ( li_k < ao_C%opi_ia ( li_iloc + 1 ) )
print*,'li_indexC=',li_index                   
                    ! on recupere la colonne de A
                    li_j = ao_C%opi_ja ( li_k )
                
                    ! copie la valeur de Aij
                    this%opi_ja ( li_index ) = li_j 
                    this%opr_a  ( li_index ) = ao_C%opr_a ( li_k )                     
                    
                    ! incrementation de l'indice 
                    li_index = li_index + 1
                    li_k = li_k + 1                    
                
                end do
            !*********************************************************************    
print*,'done'                                
            !*********************************************************************    
            ! copie des termes de D
            !*********************************************************************                            
                ! initialisation de la colonne
                li_k = ao_D%opi_ia ( li_iloc )
            
                do while ( li_k < ao_D%opi_ia ( li_iloc + 1 ) )
print*,'li_indexD=',li_index                   
                    ! on recupere la colonne de A
                    li_j = ao_D%opi_ja ( li_k )
                    
                    li_jloc = li_j + ao_A%oi_nC                    
                
                    ! copie la valeur de Aij
                    this%opi_ja ( li_index ) = li_jloc 
                    this%opr_a  ( li_index ) = ao_D%opr_a ( li_k )                     
                    
                    ! incrementation de l'indice 
                    li_index = li_index + 1
                    li_k = li_k + 1                    
                
                end do 
                this%opi_ia ( li_i ) = li_jloc                
            !*********************************************************************                                   
print*,'-------- done'                
            end if
        
        end do
        
        ! mise a jour du dernier terme
        this%opi_ia ( this%oi_nR + 1 ) = li_index
print*,'this%opi_ia=',this%opi_ia(:)        
        
	end subroutine create_SparseMatrixFrom4
!---------------------------------------------------------------------------------
    !> \todo a tester
    !> create the matrix :
    !>  A11 , A12, A13, A14 
    !>  A21 , A22, A23, A24
    !>  A31 , A32, A33, A34
    !>  A41 , A42, A43, A44        
    !> INPUTS : Aij, i,j in {1,2,3,4}
    !> OUTPUT : this         
    !> to do so, we create 4 matrices:
    !> lo_M11 = [[A11 , A12],[A21 , A22]]
    !> lo_M12 = [[A13 , A14],[A23 , A24]]
    !> lo_M21 = [[A31 , A32],[A41 , A42]]
    !> lo_M22 = [[A33 , A34],[A43 , A44]] 
    !> and then create "this" = [[M11, M12], [M21, M22]]            
	subroutine create_SparseMatrixFrom16 ( this                             &
                                         , ao_A11 , ao_A12, ao_A13, ao_A14  &
                                         , ao_A21 , ao_A22, ao_A23, ao_A24  &
                                         , ao_A31 , ao_A32, ao_A33, ao_A34  &
                                         , ao_A41 , ao_A42, ao_A43, ao_A44  )          
	implicit none
		!> param[inout] this : CSR MATRIX STRUCTURE
		type(csr_matrix) :: this		
        !> param[in] A11 , A12, A13, A14 : INPUT MATRICES
        !> param[in] A21 , A22, A23, A24 : INPUT MATRICES
        !> param[in] A31 , A32, A33, A34 : INPUT MATRICES
        !> param[in] A41 , A42, A43, A44 : INPUT MATRICES        
		type(csr_matrix) ::  ao_A11 , ao_A12, ao_A13, ao_A14
		type(csr_matrix) ::  ao_A21 , ao_A22, ao_A23, ao_A24
		type(csr_matrix) ::  ao_A31 , ao_A32, ao_A33, ao_A34
		type(csr_matrix) ::  ao_A41 , ao_A42, ao_A43, ao_A44        		
		!local var
		type(csr_matrix) ::  lo_M11
		type(csr_matrix) ::  lo_M12        
		type(csr_matrix) ::  lo_M21
		type(csr_matrix) ::  lo_M22        		        
        
        call create_SparseMatrixFrom4 ( lo_M11             &
                                      , ao_A11 , ao_A12    &
                                      , ao_A21 , ao_A22    )
        call printProfil_csrMatrix ( lo_M11 , "lo_M11", ai_typeprint = 1 )    
                                          
        call create_SparseMatrixFrom4 ( lo_M12             &
                                      , ao_A13 , ao_A14    &
                                      , ao_A23 , ao_A24    )
        call printProfil_csrMatrix ( lo_M12 , "lo_M12", ai_typeprint = 1 )    
        
        call create_SparseMatrixFrom4 ( lo_M21             &
                                      , ao_A31 , ao_A32    &
                                      , ao_A41 , ao_A42    )
        call printProfil_csrMatrix ( lo_M21 , "lo_M21", ai_typeprint = 1 )    
                                                                                
        call create_SparseMatrixFrom4 ( lo_M22             &
                                      , ao_A33 , ao_A34    &
                                      , ao_A43 , ao_A44    )
        call printProfil_csrMatrix ( lo_M22 , "lo_M22", ai_typeprint = 1 )    
        
        call create_SparseMatrixFrom4 ( this               &
                                      , lo_M11 , lo_M12    &
                                      , lo_M21 , lo_M22    )                                      
       
        
	end subroutine create_SparseMatrixFrom16    
!---------------------------------------------------------------------------------
	subroutine create_SparseMatrix_Identity ( this, ai_n )
	implicit none
		!> param[inout] this : CSR MATRIX STRUCTURE
		type(csr_matrix) :: this		
		!> param[in] ai_n : NUMBER OF COLUMNS AND ROWS
		integer :: ai_n
		!local var
		integer :: li_err,li_flag
        integer  :: li_i
        integer  :: li_j        
		
		this%ol_use_mm_format = .FALSE.
				
		this%oi_nR   = ai_n
		this%oi_nC   = ai_n		
		this%oi_nel  = ai_n
		
		allocate(this%opi_ia(this%oi_nR+1),stat=li_err)
		if (li_err.ne.0) li_flag=10	

		allocate(this%opi_ja(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=20	
				
		allocate(this%opr_a(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=30			
		
        do li_i = 1, this%oi_nR
        
            this%opi_ia ( li_i ) = li_i
            this%opi_ja ( li_i ) = li_i            
        
        end do
        
        ! mise a jour du dernier terme
        this%opi_ia ( this%oi_nR + 1 ) = this%oi_nR + 1        

		! initialisation des valeurs de la matrice
        this%opr_a ( : ) = 1.0_wp
        
	end subroutine create_SparseMatrix_Identity
!---------------------------------------------------------------------------------
    !> in this routine we must have nR = nC = n
	subroutine create_SparseMatrix_Zero ( this, ai_n )
	implicit none
		!> param[inout] this : CSR MATRIX STRUCTURE
		type(csr_matrix) :: this		
		!> param[in] ai_n : NUMBER OF COLUMNS AND ROWS
		integer :: ai_n
		!local var
		integer :: li_err,li_flag
        integer  :: li_i
        integer  :: li_j        

		this%ol_use_mm_format = .FALSE.
				
		this%oi_nR   = ai_n
		this%oi_nC   = ai_n		
		this%oi_nel  = ai_n
		
		allocate(this%opi_ia(this%oi_nR+1),stat=li_err)
		if (li_err.ne.0) li_flag=10	

		allocate(this%opi_ja(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=20	
				
		allocate(this%opr_a(this%oi_nel),stat=li_err)
		if (li_err.ne.0) li_flag=30			
		
        do li_i = 1, this%oi_nR
        
            this%opi_ia ( li_i ) = li_i
            this%opi_ja ( li_i ) = li_i            
        
        end do
        
        ! mise a jour du dernier terme
        this%opi_ia ( this%oi_nR + 1 ) = this%oi_nR + 1        

		! initialisation des valeurs de la matrice
        this%opr_a ( : ) = 0.0_wp
        
	end subroutine create_SparseMatrix_Zero        