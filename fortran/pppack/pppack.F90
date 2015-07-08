!     
! File:   pppack.F90
! Author: abdelkaderratnani
!
! Created on November 11, 2011, 5:08 PM
!
!
!	splines_tools.F90
!> THIS MODULES CONTAINS ALL FUNCTIONS REALETED TO SPLINES

module pppack
    use used_precision
    implicit none

contains

    include "banfac.F90"
    include "banslv.F90"
    include "bsplvb.F90"
    include "bsplvd.F90"
    include "bvalue.F90"
    include "interv.F90"
    include "splint.F90"
    include "spli2d.F90"
    include "splint2d.F90"
    include "bsplvb_vect.F90"
    include "bsplvd_vect.F90"

!    include "bsplpp.F90"
!    include "colloc.F90"
!    include "cubspl.F90"

end module pppack
