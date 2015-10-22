!# -*- coding: utf8 -*-
MODULE SPI_GLOBAL_DEF
IMPLICIT NONE

    INTEGER, PARAMETER    :: SPI_RK=KIND(1.D0)
    INTEGER, PARAMETER    :: SPI_INT_DEFAULT=-10000000
    INTEGER, PARAMETER    :: SPI_INT_MAX_NNZ=1000

    ! ... Error Flags
    INTEGER, PARAMETER :: SPI_SUCCESS                         =  1

    REAL(KIND=SPI_RK)   :: SPI_PI= 4.0D0*ATAN(1.D0)
  
    !!!!!!!!!!!!!!!!!!! ID parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ... QUADRATURE RULES
    INTEGER, PARAMETER :: SPI_QUADRATURES_DEFAULT        = 0

    INTEGER, PARAMETER :: SPI_QUADRATURES_LEGENDRE       = 1 
    INTEGER, PARAMETER :: SPI_QUADRATURES_LOBATTO        = 2 
    INTEGER, PARAMETER :: SPI_QUADRATURES_FOURIER        = 3
    INTEGER, PARAMETER :: SPI_QUADRATURES_FILE           = 4

    !!!!!!!!!!!!!!!!! VALUES !!!!!!!!!!!!!!!

    ! ... boundary conditions
    INTEGER, PARAMETER  :: SPI_BC_DIRICHLET_HOMOGEN      = 0 
    INTEGER, PARAMETER  :: SPI_BC_PERIODIC               = 1	
    INTEGER, PARAMETER  :: SPI_BC_NONE                   = 2

    ! ... MATRIX FORMAT OUTPUT
    INTEGER, PARAMETER  :: SPI_MATRIX_OUTPUT_FORMAT_MM   = 0


END MODULE SPI_GLOBAL_DEF
