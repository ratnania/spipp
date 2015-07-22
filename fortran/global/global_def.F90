!# -*- coding: utf8 -*-
MODULE SPI_GLOBAL_DEF
IMPLICIT NONE

    INTEGER, PARAMETER    :: SPI_RK=KIND(1.D0)
    INTEGER, PARAMETER    :: SPI_INT_DEFAULT=-10000000

    ! ... Error Flags
    INTEGER, PARAMETER :: SPI_SUCCESS                         =  1

    REAL(KIND=SPI_RK)   :: SPI_PI= 4.0D0*ATAN(1.D0)
  
    !!!!!!!!!!!!!!!!!!! ID parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ... QUADRATURE RULES
    INTEGER, PARAMETER :: SPI_QUADRATURES_DEFAULT           = 0

    INTEGER, PARAMETER :: SPI_QUADRATURES_LEGENDRE     = 1 
    INTEGER, PARAMETER :: SPI_QUADRATURES_LOBATTO      = 2 
    INTEGER, PARAMETER :: SPI_QUADRATURES_FOURIER      = 3
    INTEGER, PARAMETER :: SPI_QUADRATURES_FILE         = 4

    !!!!!!!!!!!!!!!!! VALUES !!!!!!!!!!!!!!!

    ! ... basis type 
    INTEGER, PARAMETER  :: SPI_INT_BASIS_BSPLINES        = 1 ! B-Splines 
    INTEGER, PARAMETER  :: SPI_INT_BASIS_HBEZIER         = 2 ! Hermite-Bezier 
    INTEGER, PARAMETER  :: SPI_INT_BASIS_FOURIER         = 3 ! fourier

    ! ... boundary conditions
    INTEGER, PARAMETER  :: SPI_BC_DIRICHLET_HOMOGEN  = 0 
    INTEGER, PARAMETER  :: SPI_BC_PERIODIC           = 1	
    INTEGER, PARAMETER  :: SPI_BC_NO                 = 2


END MODULE SPI_GLOBAL_DEF
