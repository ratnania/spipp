!# -*- coding: utf8 -*-
MODULE SPI_GLOBAL_DEF
IMPLICIT NONE

    INTEGER, PARAMETER    :: SPI_RK=KIND(1.D0)
    INTEGER, PARAMETER    :: SPI_INT_DEFAULT=-10000000

    ! ... Error Flags
    INTEGER, PARAMETER :: SPI_SUCCESS                         =  1
  
    !!!!!!!!!!!!!!!!!!! ID parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ... QUADRATURE RULES
    INTEGER, PARAMETER :: SPI_QUADRATURES_DEFAULT           = 0

    INTEGER, PARAMETER :: SPI_QUADRATURES_QUAD_LEGENDRE     = 1 
    INTEGER, PARAMETER :: SPI_QUADRATURES_QUAD_LOBATTO      = 2 
    INTEGER, PARAMETER :: SPI_QUADRATURES_QUAD_FOURIER      = 3
    INTEGER, PARAMETER :: SPI_QUADRATURES_QUAD_FILE         = 4

    !!!!!!!!!!!!!!!!! VALUES !!!!!!!!!!!!!!!

    ! ... basis type 
    INTEGER, PARAMETER  :: SPI_BASIS_BSPLINES        = 1 ! B-Splines 
    INTEGER, PARAMETER  :: SPI_BASIS_FOURIER         = 2 ! fourier

END MODULE SPI_GLOBAL_DEF
