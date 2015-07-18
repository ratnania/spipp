!# -*- coding: utf8 -*-
MODULE SPI_FOURIER_MODS
  USE SPI_GLOBAL_DEF
  USE SPI_GLOBAL
  IMPLICIT NONE

CONTAINS
! ..........................................................        
  SUBROUTINE Compute_fourier_mods(n_tor,n_period,i_tor,ar_plane,ar_mod)
  IMPLICIT NONE
    INTEGER :: n_tor, i_tor, n_period
    REAL(KIND=SPI_RK) :: ar_plane, ar_mod
    INTEGER :: mode(n_tor),il 


    mode(1)=0
    DO il=2,n_tor
       mode(il)=int(il/2)/n_period
    END DO
    
    IF(i_tor .eq. 1) THEN
       ar_mod  = 1.d0
    END IF

    IF((mod(i_tor,2) .eq. 0) .and. i_tor > 1) THEN
       ar_mod  = cos(mode(i_tor)*ar_plane)
    END IF

    IF((mod(i_tor,2) .eq. 1) .and. i_tor > 1) THEN
       ar_mod  = sin(mode(i_tor)*ar_plane)
    END IF

  END SUBROUTINE Compute_fourier_mods
! ..........................................................        

! ..........................................................        
  SUBROUTINE Compute_FD_fourier_mods(n_tor,n_period,i_tor,ar_plane,ar_mod)
  IMPLICIT NONE
    INTEGER :: n_tor, i_tor, n_period
    REAL(KIND=SPI_RK) :: ar_plane, ar_mod
    INTEGER :: mode(n_tor),il 

    mode(1)=0
    DO il=2,n_tor
       mode(il)=int(il/2)/n_period
    END DO
    
    IF(i_tor .eq. 1) THEN
       ar_mod  = 0.d0
    END IF

    IF((mod(i_tor,2) .eq. 0) .and. i_tor > 1) THEN
       ar_mod  = - float(mode(i_tor)) * sin(mode(i_tor)*ar_plane)
    END IF

    IF((mod(i_tor,2) .eq. 1) .and. i_tor > 1) THEN
       ar_mod  =  float(mode(i_tor)) * cos(mode(i_tor)*ar_plane)
    END IF

  END SUBROUTINE Compute_FD_fourier_mods
! ..........................................................        

! ..........................................................        
 SUBROUTINE Compute_SD_fourier_mods(n_tor,n_period,i_tor,ar_plane,ar_mod)
  IMPLICIT NONE
    INTEGER :: n_tor, i_tor, n_period
    REAL(KIND=SPI_RK) :: ar_plane, ar_mod
    INTEGER :: mode(n_tor),il 

    mode(1)=0
    DO il=2,n_tor
       mode(il)=int(il/2)/n_period
    END DO
    
    IF(i_tor .eq. 1) THEN
       ar_mod  = 0.d0
    END IF

    IF((mod(i_tor,2) .eq. 0) .and. i_tor > 1) THEN
       ar_mod  = - float(mode(i_tor))**2 * cos(mode(i_tor)*ar_plane)
    END IF

    IF((mod(i_tor,2) .eq. 1) .and. i_tor > 1) THEN
       ar_mod  = - float(mode(i_tor))**2 * sin(mode(i_tor)*ar_plane)
    END IF

  END SUBROUTINE Compute_SD_fourier_mods
! ..........................................................        

END MODULE SPI_FOURIER_MODS
