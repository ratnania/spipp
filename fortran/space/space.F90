!# -*- coding: utf8 -*-
MODULE SPI_SPACE 
  USE SPI_GLOBAL_DEF
  USE SPI_GLOBAL
  USE SPI_QUADRATURES_DEF
  USE SPI_QUADRATURES
  USE SPI_MESH_DEF
  USE SPI_MESH
  USE SPI_NUMBERING_DEF
  USE SPI_NUMBERING
  USE SPI_SPACE_DEF
  IMPLICIT NONE

!  PRIVATE

  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 2 

CONTAINS
  ! ..................................................

  ! .........................................................
  SUBROUTINE CREATE_SPACE(self, ao_mesh)
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_SPACE_ABSTRACT)  , INTENT(INOUT) :: self 
     CLASS(DEF_MESH_ABSTRACT)       , INTENT(INOUT) :: ao_mesh 
     ! LOCAL

     SELECT TYPE (self)
     CLASS IS (DEF_SPACE_1D_BSPLINE)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D_BSPLINE)
           CALL CREATE_SPACE_1D_BSPLINE(self, ao_mesh)
        CLASS DEFAULT
           STOP 'CREATE_SPACE: unexpected type for ao_mesh object!'
        END SELECT

     CLASS IS (DEF_SPACE_1D_FOURIER)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D_FOURIER)
           CALL CREATE_SPACE_1D_FOURIER(self, ao_mesh)
        CLASS DEFAULT
           STOP 'CREATE_SPACE: unexpected type for ao_mesh object!'
        END SELECT

     CLASS IS (DEF_SPACE_1D_HBEZIER)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D_HBEZIER)
           CALL CREATE_SPACE_1D_HBEZIER(self, ao_mesh)
        CLASS DEFAULT
           STOP 'CREATE_SPACE: unexpected type for ao_mesh object!'
        END SELECT

     CLASS DEFAULT
        STOP 'CREATE_SPACE: unexpected type for self object!'
     END SELECT

  END SUBROUTINE CREATE_SPACE
  ! .........................................................

  ! ........................................................
  SUBROUTINE CREATE_SPACE_1D_BSPLINE( self, ao_mesh)
  implicit none
     type(DEF_SPACE_1D_BSPLINE),intent(inout) :: self
     TYPE(DEF_MESH_1D_BSPLINE), TARGET, INTENT(IN) :: ao_mesh
     ! LOCAL VARIABLES

  
  END SUBROUTINE CREATE_SPACE_1D_BSPLINE
  ! ........................................................

  ! .........................................................
  SUBROUTINE CREATE_SPACE_1D_FOURIER(self, ao_mesh)
  IMPLICIT NONE
     TYPE(DEF_SPACE_1D_FOURIER), INTENT(INOUT)  :: self
     TYPE(DEF_MESH_1D_FOURIER), TARGET, INTENT(IN) :: ao_mesh
     ! LOCAL
     INTEGER :: li_err 

  END SUBROUTINE CREATE_SPACE_1D_FOURIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_SPACE_1D_HBEZIER(self, ao_mesh)
  IMPLICIT NONE
     TYPE(DEF_SPACE_1D_HBEZIER), INTENT(INOUT)  :: self
     TYPE(DEF_MESH_1D_HBEZIER), TARGET, INTENT(IN) :: ao_mesh
     ! LOCAL
     INTEGER :: li_err 

  END SUBROUTINE CREATE_SPACE_1D_HBEZIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE FREE_SPACE(self)
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_SPACE_ABSTRACT)       , INTENT(INOUT) :: self 
     ! LOCAL

      SELECT TYPE (self)
      CLASS IS (DEF_SPACE_1D_BSPLINE)
         CALL FREE_SPACE_1D_BSPLINE(self)
     
      CLASS IS (DEF_SPACE_1D_FOURIER)
         CALL FREE_SPACE_1D_FOURIER(self)
     
      CLASS IS (DEF_SPACE_1D_HBEZIER)
         CALL FREE_SPACE_1D_HBEZIER(self)
     
      CLASS DEFAULT
         STOP 'FREE_SPACE: unexpected type for self object!'
      END SELECT

  END SUBROUTINE FREE_SPACE
  ! .........................................................
         
   ! ........................................................
   SUBROUTINE FREE_SPACE_1D_BSPLINE( self )
   implicit none
      type(DEF_SPACE_1D_BSPLINE),intent(inout) :: self
      ! LOCAL VARIABLES
      
   END SUBROUTINE FREE_SPACE_1D_BSPLINE
   ! ........................................................

   ! ........................................................
   SUBROUTINE FREE_SPACE_1D_FOURIER( self )
   implicit none
      type(DEF_SPACE_1D_FOURIER),intent(inout) :: self
      ! LOCAL VARIABLES
      
   END SUBROUTINE FREE_SPACE_1D_FOURIER
   ! ........................................................

   ! ........................................................
   SUBROUTINE FREE_SPACE_1D_HBEZIER( self )
   implicit none
      type(DEF_SPACE_1D_HBEZIER),intent(inout) :: self
      ! LOCAL VARIABLES
      
   END SUBROUTINE FREE_SPACE_1D_HBEZIER
   ! ........................................................

END MODULE SPI_SPACE 
