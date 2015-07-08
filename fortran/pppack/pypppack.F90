!     
! File:   pypppack.F90
! Author: root
!
! Created on December 7, 2012, 8:37 AM
!

MODULE pypppack
    use used_precision
    use pppack
    IMPLICIT NONE
contains
    !---------------------------------------------------------------
    subroutine pyinterv(ar_x, apr_t, ai_nk, ai_left, ai_mflag)
        !> this routine interpolates a list of 1D points
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: ar_x
        INTEGER, INTENT(IN) :: ai_nk
        REAL(8), DIMENSION(ai_nk), INTENT(IN) :: apr_t
        INTEGER, INTENT(OUT) :: ai_left
        INTEGER, INTENT(OUT) :: ai_mflag

        CALL interv ( apr_t, ai_nk, ar_x, ai_left, ai_mflag )

    end subroutine pyinterv
    !---------------------------------------------------------------
    subroutine pybvalue(ar_x, ai_jderiv, ai_p, apr_t, apr_bcoef, ai_nk, ai_n, ar_val)
        !> this routine interpolates a list of 1D points
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: ar_x
        INTEGER, INTENT(IN) :: ai_jderiv
        INTEGER, INTENT(IN) :: ai_p
        INTEGER, INTENT(IN) :: ai_nk
        REAL(8), DIMENSION(ai_nk), INTENT(IN) :: apr_t
        INTEGER, INTENT(IN) :: ai_n
        REAL(8), DIMENSION(ai_n), INTENT(IN):: apr_bcoef
        REAL(8), INTENT(OUT) :: ar_val
        ! LOCAL
        INTEGER  :: li_k

        li_k =  ai_p + 1
        
        ar_val = bvalue ( apr_t, apr_bcoef, ai_n, li_k, ar_x, ai_jderiv )

    end subroutine pybvalue
    !---------------------------------------------------------------
    subroutine pybsplvb(ar_x, ai_left, apr_t, ai_nk, apr_biatx, ai_k)
        !> this routine interpolates a list of 1D points
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: ar_x
        INTEGER, INTENT(IN) :: ai_left
        INTEGER, INTENT(IN) :: ai_nk
        REAL(8), DIMENSION(ai_nk), INTENT(IN) :: apr_t
        INTEGER, INTENT(IN) :: ai_k
        REAL(8), DIMENSION(ai_k), INTENT(OUT) :: apr_biatx
        ! LOCAL
        INTEGER :: li_index

        CALL bsplvb ( apr_t, ai_k, li_index, ar_x, ai_left, apr_biatx )

    end subroutine pybsplvb
    !---------------------------------------------------------------
    subroutine pybsplvd(ar_x, ai_left, ai_nderiv, apr_t, ai_nk, apr_dbiatx, ai_k, ai_nderivp1)
        !> this routine interpolates a list of 1D points
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: ar_x
        INTEGER, INTENT(IN) :: ai_left
        INTEGER, INTENT(IN) :: ai_nderiv
        INTEGER, INTENT(IN) :: ai_nk
        REAL(8), DIMENSION(ai_nk), INTENT(IN) :: apr_t
        INTEGER, INTENT(IN) :: ai_k
        INTEGER, INTENT(IN) :: ai_nderivp1
        REAL(8), DIMENSION(ai_k, ai_nderivp1), INTENT(OUT) :: apr_dbiatx
        ! LOCAL
        REAL(8), DIMENSION(ai_k,ai_k)  :: lpr_a

        CALL bsplvd ( apr_t, ai_k, ar_x, ai_left, lpr_a, apr_dbiatx, ai_nderiv )

    end subroutine pybsplvd
    !---------------------------------------------------------------
!    subroutine pybsplpp(apr_t, ai_nk, apr_bcoef, ai_n, apr_break, apr_coef, ai_lmax, ai_k, ai_l)
!        !> this routine interpolates a list of 1D points
!        IMPLICIT NONE
!        REAL(8), DIMENSION(ai_nk), INTENT(IN) :: apr_t
!        INTEGER, INTENT(IN) :: ai_nk
!        REAL(8), DIMENSION(ai_n), INTENT(IN)  :: apr_bcoef
!        INTEGER, INTENT(IN) :: ai_n
!        REAL(8), DIMENSION(ai_lmax), INTENT(OUT)  :: apr_break
!        REAL(8), DIMENSION(ai_k,ai_lmax), INTENT(OUT)  :: apr_coef
!        INTEGER, INTENT(IN) :: ai_lmax
!        INTEGER, INTENT(IN) :: ai_k
!        INTEGER, INTENT(OUT) :: ai_l
!        ! LOCAL
!        REAL(8), DIMENSION(ai_k,ai_k)  :: lpr_a
!        REAL(8), DIMENSION(ai_k,ai_k)  :: lpr_scrtch
!
!        CALL bsplpp ( apr_t, apr_bcoef, ai_n, ai_k, lpr_scrtch, apr_break, apr_coef, ai_l )
!
!    end subroutine pybsplpp
    !---------------------------------------------------------------
!    subroutine pycolloc(ar_aleft, ar_aright, ai_lbegin, ai_iorder, ai_ntimes, ar_addbrk, ar_relerr)
!        !> this routine interpolates a list of 1D points
!        IMPLICIT NONE
!        REAL(8), INTENT(IN) :: ar_aleft
!        REAL(8), INTENT(IN) :: ar_aright
!        INTEGER, INTENT(IN) :: ai_lbegin
!        INTEGER, INTENT(IN) :: ai_iorder
!        INTEGER, INTENT(IN) :: ai_ntimes
!        REAL(8), INTENT(IN) :: ar_addbrk
!        REAL(8), INTENT(IN) :: ar_relerr
!
!        CALL colloc(ar_aleft, ar_aright, ai_lbegin, ai_iorder, ai_ntimes, ar_addbrk, ar_relerr)
!
!    end subroutine pycolloc
    !---------------------------------------------------------------
!    subroutine pycubspl(apr_tau, apr_c, ai_n, ai_ibcbeg, ai_ibcend)
!        !> this routine interpolates a list of 1D points
!        IMPLICIT NONE
!        INTEGER, INTENT(IN) :: ai_n
!        INTEGER, INTENT(IN) :: ai_p
!        REAL(8), DIMENSION(ai_nk), INTENT(IN) :: apr_t
!        REAL(8), DIMENSION(ai_npoint), INTENT(IN) :: apr_tau
!        REAL(8), DIMENSION(ai_npoint, ai_dim), INTENT(IN) :: apr_gtau
!        INTEGER, INTENT(IN) :: ai_nk
!        REAL(8), DIMENSION(ai_npoint, ai_dim), INTENT(OUT) :: apr_P
!        INTEGER, INTENT(IN) :: ai_npoint
!        INTEGER, INTENT(IN) :: ai_dim
!
!        CALL cubspl ( apr_tau, apr_c, ai_n, ai_ibcbeg, ai_ibcend )
!
!    end subroutine pycubspl
    !---------------------------------------------------------------
    subroutine pysplint(ai_n, ai_p, apr_t, apr_tau, apr_gtau, ai_nk, apr_P, ai_npoint, ai_dim)
        !> this routine interpolates a list of 1D points
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ai_n
        INTEGER, INTENT(IN) :: ai_p
        INTEGER, INTENT(IN) :: ai_nk
        REAL(8), DIMENSION(ai_nk), INTENT(IN) :: apr_t
        INTEGER, INTENT(IN) :: ai_npoint
        INTEGER, INTENT(IN) :: ai_dim
        REAL(8), DIMENSION(ai_npoint), INTENT(IN) :: apr_tau
        REAL(8), DIMENSION(ai_npoint, ai_dim), INTENT(IN) :: apr_gtau
        REAL(8), DIMENSION(ai_npoint, ai_dim), INTENT(OUT) :: apr_P
        ! LOCAL
        REAL(8), DIMENSION(:), POINTER :: lpr_bcoef
        REAL(8), DIMENSION(:), POINTER :: lpr_q
        INTEGER  :: li_i
        INTEGER  :: li_n
        INTEGER  :: li_k
        INTEGER  :: li_iflag

        li_n =  ai_n
        li_k =  ai_p + 1

        IF (ai_npoint /= li_n) THEN
            PRINT*, 'ERROR pysplint : ai_npoint must be equal to li_n'
            STOP
        END IF

        ALLOCATE(lpr_bcoef(li_n))
        ALLOCATE(lpr_q((2*li_k-1)*li_n))

        DO li_i = 1, ai_dim
            lpr_bcoef = 0.0_8
            lpr_q = 0.0_8
            CALL splint ( apr_tau, apr_gtau(:,li_i), apr_t, li_n, li_k, lpr_q, lpr_bcoef, li_iflag )
            apr_P(1:li_n, li_i) = lpr_bcoef(1:li_n)
        END DO

        DEALLOCATE(lpr_bcoef)
        DEALLOCATE(lpr_q)

    end subroutine pysplint
    !---------------------------------------------------------------
    subroutine pysplint2d(ai_nx, ai_kx, apr_taux, ai_ny, ai_ky, apr_tauy, apr_g, apr_Bcoef, apr_tx, apr_ty, ai_mx, ai_my, ai_nkx, ai_nky)
        !> this routine interpolates a list of 1D points
        IMPLICIT NONE
        ! INPUT
        INTEGER, INTENT(IN) :: ai_nx
        INTEGER, INTENT(IN) :: ai_kx
        REAL(8), DIMENSION (ai_nx), INTENT(IN) :: apr_taux
        INTEGER, INTENT(IN) :: ai_ny
        INTEGER, INTENT(IN) :: ai_ky
        REAL(8), DIMENSION (ai_ny), INTENT(IN) :: apr_tauy
        REAL(8), DIMENSION (ai_nx, ai_ny), INTENT(IN) :: apr_g
        ! OUTPUT
        INTEGER, INTENT(IN) :: ai_mx
        INTEGER, INTENT(IN) :: ai_my
        INTEGER, INTENT(IN) :: ai_nkx
        INTEGER, INTENT(IN) :: ai_nky
        REAL(8), DIMENSION (ai_mx, ai_my), INTENT(OUT) :: apr_Bcoef
        REAL(8), DIMENSION (ai_nkx), INTENT(OUT) :: apr_tx
        REAL(8), DIMENSION (ai_nky), INTENT(OUT) :: apr_ty

        IF (ai_nx /= ai_mx) THEN
            PRINT*, 'ERROR pysplint2d : ai_nx must be equal to ai_mx'
            STOP
        END IF
        IF (ai_ny /= ai_my) THEN
            PRINT*, 'ERROR pysplint2d : ai_ny must be equal to ai_my'
            STOP
        END IF
        IF (ai_nx+ai_kx /= ai_nkx) THEN
            PRINT*, 'ERROR pysplint2d : ai_nx+ai_kx must be equal to ai_nkx'
            STOP
        END IF
        IF (ai_ny+ai_ky /= ai_nky) THEN
            PRINT*, 'ERROR pysplint2d : ai_ny+ai_ky must be equal to ai_nky'
            STOP
        END IF

        CALL splint2d(ai_nx, ai_kx, apr_taux, ai_ny, ai_ky, apr_tauy, apr_g, apr_Bcoef, apr_tx, apr_ty)

    end subroutine pysplint2d    
    !---------------------------------------------------------------
END MODULE pypppack
