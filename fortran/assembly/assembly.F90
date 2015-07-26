!# -*- coding: utf8 -*-
MODULE SPI_ASSEMBLY
!  USE FIELD_DEF
!  USE FIELD
  USE SPI_GLOBAL_DEF
  USE SPI_ASSEMBLY_DEF
  USE SPI_ASSEMBLY_1D_BSPLINE
  USE SPI_GREENBOX_DEF
  USE SPI_GREENBOX

  IMPLICIT NONE

CONTAINS

   ! ...................................................
   SUBROUTINE No_RHS_for_Vi_1D(ao_BBox1D, ao_GBox1D)
   IMPLICIT NONE
      TYPE(DEF_BLACKBOX_1D) :: ao_BBox1D
      TYPE(DEF_GREENBOX_1D)                :: ao_GBox1D

      RETURN
   END SUBROUTINE No_RHS_for_Vi_1D
   ! ...................................................

   ! ...................................................
   SUBROUTINE No_Matrix_for_Vi_Vj_1D(ao_BBox1Di, ao_BBox1Dj, ao_GBox1D)
   IMPLICIT NONE
      TYPE(DEF_BLACKBOX_1D) :: ao_BBox1Di
      TYPE(DEF_BLACKBOX_1D) :: ao_BBox1Dj
      TYPE(DEF_GREENBOX_1D)                :: ao_GBox1D

      RETURN
   END SUBROUTINE No_Matrix_for_Vi_Vj_1D
   ! ...................................................

   ! ...................................................
!   SUBROUTINE LOOP_ON_ELMTS_GENERIC_1D(ptr_matrix, &
!                   & ao_hashtab_matrix, ao_hashtab_matrix_field, ai_color, &
!                   & Mesh1D, ao_GBox1D, &
!                   & ao_Basis1Di, ao_Basis1Dj)
   SUBROUTINE LOOP_ON_ELMTS_GENERIC_1D(ao_matrix, &
                  & ao_mesh, &
                  & ao_trial_space, & 
                  & ao_test_space, & 
                  & ao_gbox)
     IMPLICIT NONE
     TYPE(DEF_MATRIX_1D), TARGET, INTENT(INOUT) :: ao_matrix
     CLASS(DEF_MESH_ABSTRACT), INTENT(INOUT) :: ao_mesh
     CLASS(DEF_SPACE_ABSTRACT), TARGET :: ao_trial_space
     CLASS(DEF_SPACE_ABSTRACT), TARGET :: ao_test_space
     TYPE(DEF_GREENBOX_1D)                      :: ao_gbox

!     TYPE(DEF_MATRIX_1D_CELL)                ::  ao_hashtab_matrix
!     TYPE(DEF_MATRIX_FIELD_1D_CELL)          ::  ao_hashtab_matrix_field
!     INTEGER, INTENT(IN)                   :: ai_color 
!     TYPE(DEF_MESH_1D)                     :: Mesh1D
!     TYPE(DEF_GREENBOX_1D)                 :: ao_GBox1D
!     TYPE(DEF_BASIS_1D) , TARGET           :: ao_Basis1Di
!     TYPE(DEF_BASIS_1D)   , TARGET         :: ao_Basis1Dj
     ! LOCAL
     INTEGER :: ierr
  

     SELECT TYPE (ao_mesh)
     CLASS IS (DEF_MESH_1D_BSPLINE)
        SELECT TYPE (ao_trial_space)
        CLASS IS (DEF_SPACE_1D_BSPLINE)
           SELECT TYPE (ao_test_space)
           CLASS IS (DEF_SPACE_1D_BSPLINE)
              CALL LOOP_ON_ELMTS_1D_BSPLINE(ao_matrix, &
                        & ao_mesh, &
                        & ao_trial_space, & 
                        & ao_test_space, & 
                        & ao_gbox)
           CLASS DEFAULT
              STOP 'LOOP_ON_ELMTS_GENERIC_1D: unexpected type for ao_test_space object!'
           END SELECT
        CLASS DEFAULT
           STOP 'LOOP_ON_ELMTS_GENERIC_1D: unexpected type for ao_trial_space object!'
        END SELECT

!     CLASS IS (DEF_MESH_1D_FOURIER)
!        CALL LOOP_ON_ELMTS_1D_FOURIER(ao_matrix, ao_mesh) 
!
!     CLASS IS (DEF_MESH_1D_HBEZIER)
!        CALL LOOP_ON_ELMTS_1D_HBEZIER(ao_matrix, ao_mesh) 

     CLASS DEFAULT
        STOP 'LOOP_ON_ELMTS_GENERIC_1D: unexpected type for ao_mesh object!'
     END SELECT
  
   END SUBROUTINE LOOP_ON_ELMTS_GENERIC_1D
   ! ...................................................

END MODULE SPI_ASSEMBLY
