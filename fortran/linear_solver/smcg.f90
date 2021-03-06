	!*******************************************************!
	!                                                       !
	!          METHODE DU GRADIENT CONJUGUE SIMPLE          !
	!                                                       !
	!*******************************************************!

	subroutine Gradient_conj ( this, apr_B,apr_U, ai_maxIter, ar_eps )
	implicit none
		type(csr_matrix) :: this 	
		real(wp), dimension(:) :: apr_U
		real(wp), dimension(:) :: apr_B
		integer  :: ai_maxIter
		real(wp) :: ar_eps
		!local var
		real(wp), dimension(:), pointer :: lpr_Ad
		real(wp), dimension(:), pointer :: lpr_r
		real(wp), dimension(:), pointer :: lpr_d				
		real(wp), dimension(:), pointer :: lpr_Ux			
		real(wp) :: lr_Norm2r1
		real(wp) :: lr_Norm2r0
		real(wp) :: lr_NormInfb
		real(wp) :: lr_NormInfr						
		real(wp) :: lr_ps
		real(wp) :: lr_beta
		real(wp) :: lr_alpha				
		logical  :: ll_continue
		integer  :: li_iter
		integer  :: li_err
		integer  :: li_flag		
		
		if ( this%oi_nR /= this%oi_nC ) then
			PRINT*,'ERROR Gradient_conj: The matrix must be square'
			stop
		end if

		if ( ( dabs ( MAXVAL ( apr_B ) ) < ar_eps ) .AND. ( dabs ( MINVAL ( apr_B ) ) < ar_eps ) ) then
			apr_U = 0.0_wp
			return
		end if

		allocate(lpr_Ad(this%oi_nR),stat=li_err)
		if (li_err.ne.0) li_flag=10	
		allocate(lpr_r(this%oi_nR),stat=li_err)
		if (li_err.ne.0) li_flag=20	
		allocate(lpr_d(this%oi_nR),stat=li_err)
		if (li_err.ne.0) li_flag=30		
		allocate(lpr_Ux(this%oi_nR),stat=li_err)
		if (li_err.ne.0) li_flag=40					
		
		!================!
		! initialisation !
		!================!
		lpr_Ux = apr_U
		li_iter = 0
		call Mult_CSR_Matrix_Vector( this , lpr_Ux , lpr_Ad )
		
		!-------------------!
		! calcul des normes !
		!-------------------!
		lpr_r		= apr_B - lpr_Ad
		lr_Norm2r0  = DOT_PRODUCT( lpr_r , lpr_r )
		lr_NormInfb = maxval( dabs( apr_B ) )

		lpr_d = lpr_r
		!================!
		ll_continue=.true.
		do while(ll_continue)
			li_iter = li_iter + 1
			!--------------------------------------!
			! calcul du ak parametre optimal local !
			!--------------------------------------!

			call Mult_CSR_Matrix_Vector( this , lpr_d , lpr_Ad )
			
			lr_ps = DOT_PRODUCT( lpr_Ad , lpr_d )
			lr_alpha = lr_Norm2r0 / lr_ps
			
			!==================================================!
			! calcul de l'approximation Xk+1 et du residu Rk+1 !
			!==================================================!
			! calcul des composantes residuelles
			!-----------------------------------
			lpr_r = lpr_r - lr_alpha * lpr_Ad
			
			!----------------------------------------!
			! approximations ponctuelles au rang k+1 !
			!----------------------------------------!
			lpr_Ux = lpr_Ux + lr_alpha * lpr_d
			
			!-------------------------------------------------------!
			! (a) extraction de la norme infinie du residu          !
			!     pour le test d'arret                              !
			! (b) extraction de la norme euclidienne du residu rk+1 !
			!-------------------------------------------------------!
			lr_NormInfr = maxval(dabs( lpr_r ))
			lr_Norm2r1 = DOT_PRODUCT( lpr_r , lpr_r )
			
			!==================================================!
			! calcul de la nouvelle direction de descente dk+1 !
			!==================================================!
			lr_beta = lr_Norm2r1 / lr_Norm2r0
			lr_Norm2r0 = lr_Norm2r1
			lpr_d = lpr_r + lr_beta * lpr_d
			
			!-------------------!
			! boucle suivante ? !
			!-------------------!
			ll_continue=( ( lr_NormInfr / lr_NormInfb ) >= ar_eps ) .AND. ( li_iter < ai_maxIter )
		end do
		apr_U = lpr_Ux
		
		if ( li_iter == ai_maxIter ) then
			print*,'Warning Gradient_conj : li_iter == ai_maxIter'
			print*,'Error after CG =',( lr_NormInfr / lr_NormInfb )
		end if

        call printlog ( "Gradient_conj : duree du calcul =", ai_dtllevel = 1 )            

        call printcputime ( )              		
		
		deallocate(lpr_Ad)
		deallocate(lpr_d)
		deallocate(lpr_r)				
		deallocate(lpr_Ux)			
	end subroutine Gradient_conj