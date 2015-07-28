MODULE SPI_NUMBERING_DEF 
USE SPI_GLOBAL_DEF 
USE SPI_MESH_DEF
IMPLICIT NONE

   ! ...................................................
   TYPE, ABSTRACT, PUBLIC :: DEF_NUMBERING_ABSTRACT
      !> DIMENSION OF THE DISCRETE SPACE IN EACH DIRECTION
      integer		:: oi_size
      
      !> LOCAL TO GLOBAL INDEXATION
      integer, dimension( : , : ), allocatable :: opi_IEN
      !> THE GLOBAL INDEX. WE MUST REMOVE DIRICHLET BASIS FUNCTIONS
      integer, dimension( : ), allocatable :: opi_ID	
      !> LOCAL MATRIX, LM ( e , b ) = ID ( IEN ( e , b ) )
      integer, dimension( : , : ), allocatable :: opi_LM

      !> THE LEFT KNOT INDEX CORRESPONDING TO A NON ZERO-MEASURE ELEMENT
      integer, dimension( : ), allocatable :: opi_ELT_INDEX		
      
   END TYPE DEF_NUMBERING_ABSTRACT
   ! ...................................................

   ! .........................................................
   TYPE, PUBLIC, EXTENDS(DEF_NUMBERING_ABSTRACT) :: DEF_NUMBERING_1D_BSPLINE
      CLASS(DEF_MESH_1D_BSPLINE), POINTER :: ptr_mesh => NULL()
   END TYPE DEF_NUMBERING_1D_BSPLINE
   ! .........................................................

   ! .........................................................
   TYPE, PUBLIC, EXTENDS(DEF_NUMBERING_ABSTRACT) :: DEF_NUMBERING_1D_FOURIER
      CLASS(DEF_MESH_1D_FOURIER), POINTER :: ptr_mesh => NULL()
   END TYPE DEF_NUMBERING_1D_FOURIER
   ! .........................................................

   ! .........................................................
   TYPE, PUBLIC, EXTENDS(DEF_NUMBERING_ABSTRACT) :: DEF_NUMBERING_1D_HBEZIER
      CLASS(DEF_MESH_1D_HBEZIER), POINTER :: ptr_mesh => NULL()
   END TYPE DEF_NUMBERING_1D_HBEZIER
   ! .........................................................

CONTAINS

  ! .............................................
  INTEGER FUNCTION SPI_IDDL_1D(ai_ivar, ai_dof, ai_nvar)
    INTEGER :: ai_dof
    INTEGER :: ai_ivar
    INTEGER :: ai_nvar  
    ! LOCAL

    SPI_IDDL_1D = ai_ivar + ai_nvar*(ai_dof-1) 
  END FUNCTION SPI_IDDL_1D
  ! .............................................


END MODULE SPI_NUMBERING_DEF
