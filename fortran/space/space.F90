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
  USE SPI_BASIS_DEF
  USE SPI_BASIS
  USE SPI_BLACKBOX_DEF
  USE SPI_BLACKBOX
  IMPLICIT NONE

!  PRIVATE

  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 2 

CONTAINS
  ! ..................................................

  ! .........................................................
  SUBROUTINE CREATE_SPACE(self, ao_mesh, ai_type, ai_k)
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_SPACE_ABSTRACT), INTENT(INOUT) :: self 
     CLASS(DEF_MESH_ABSTRACT) , INTENT(INOUT) :: ao_mesh 
     INTEGER                  , INTENT(IN)    :: ai_type
     INTEGER                   , OPTIONAL, INTENT(IN)    :: ai_k
     ! LOCAL

     SELECT TYPE (self)
     CLASS IS (DEF_SPACE_1D_BSPLINE)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D_BSPLINE)
           IF (PRESENT(ai_k)) THEN
              CALL CREATE_SPACE_1D_BSPLINE(self, ao_mesh, ai_type, ai_k)
           ELSE
              STOP "CREATE_SPACE: Missing argument"
           END IF
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
  SUBROUTINE CREATE_SPACE_1D_BSPLINE(self, ao_mesh, ai_type, ai_k)
  implicit none
     type(DEF_SPACE_1D_BSPLINE), intent(inout) :: self
     TYPE(DEF_MESH_1D_BSPLINE) , TARGET, INTENT(INOUT) :: ao_mesh
     INTEGER                   , INTENT(IN)    :: ai_type
     INTEGER                   , INTENT(IN)    :: ai_k
     ! LOCAL VARIABLES

     self % ptr_mesh => ao_mesh

     CALL CREATE_NUMBERING(self % oo_numbering, ao_mesh)
     CALL CREATE_QUADRATURE(self % oo_quad, ai_type, ai_k)
     CALL CREATE_BASIS(self % oo_basis, self % ptr_mesh, self % oo_quad) 
     CALL CREATE_BLACKBOX(self % oo_bbox, self % oo_basis, self % oo_quad)

  END SUBROUTINE CREATE_SPACE_1D_BSPLINE
  ! ........................................................

  ! .........................................................
  SUBROUTINE CREATE_SPACE_1D_FOURIER(self, ao_mesh)
  IMPLICIT NONE
     TYPE(DEF_SPACE_1D_FOURIER), INTENT(INOUT)  :: self
     TYPE(DEF_MESH_1D_FOURIER), TARGET, INTENT(INOUT) :: ao_mesh
     ! LOCAL
     INTEGER :: li_err 

     self % ptr_mesh => ao_mesh

  END SUBROUTINE CREATE_SPACE_1D_FOURIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_SPACE_1D_HBEZIER(self, ao_mesh)
  IMPLICIT NONE
     TYPE(DEF_SPACE_1D_HBEZIER), INTENT(INOUT)  :: self
     TYPE(DEF_MESH_1D_HBEZIER), TARGET, INTENT(INOUT) :: ao_mesh
     ! LOCAL
     INTEGER :: li_err 

     self % ptr_mesh => ao_mesh

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

      CALL FREE_NUMBERING(self % oo_numbering)
      CALL FREE_QUADRATURE(self % oo_quad)
      CALL FREE_BASIS(self % oo_basis)
      CALL FREE_BLACKBOX(self % oo_bbox)
      
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
