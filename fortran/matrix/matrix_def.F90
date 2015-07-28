!# -*- coding: utf8 -*-
MODULE SPI_MATRIX_DEF
  USE SPI_GLOBAL_DEF
  USE SPI_GLOBAL
  USE SPI_QUADRATURES_DEF
  USE SPI_QUADRATURES
  USE SPI_MESH_DEF
  USE SPI_MESH
  USE SPI_NUMBERING_DEF
  USE SPI_NUMBERING
  USE SPI_SPARSE_MATRIX_DEF
  USE SPI_SPARSE_MATRIX
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 2 

  ! ..................................................
  TYPE, PUBLIC, ABSTRACT :: DEF_MATRIX_ABSTRACT
     ! ... matrix ID in SPM
     INTEGER :: oi_matrix_id 
     INTEGER :: oi_nvar
     LOGICAL :: ol_to_assembly

     INTEGER :: oi_ndof_rows 
     INTEGER :: oi_ndof_cols
     INTEGER :: oi_localsize_rows 
     INTEGER :: oi_localsize_cols

     ! ... Local Matrix Data 
     LOGICAL :: ol_allocated_locmat
     INTEGER :: oi_npush_locmat
     REAL(KIND=SPI_RK), DIMENSION(:, :), POINTER :: opr_locmat
     ! ...

     REAL(KIND=SPI_RK), DIMENSION(:)  , ALLOCATABLE  :: opr_global_rhs 
     REAL(KIND=SPI_RK), DIMENSION(:)  , ALLOCATABLE  :: opr_global_unknown

     REAL(KIND=SPI_RK), DIMENSION(:,:), POINTER :: Matrix_Contribution 
     REAL(KIND=SPI_RK), DIMENSION(:)  , POINTER :: Rhs_Contribution

     TYPE(DEF_MATRIX_CSR) :: oo_csr
  END TYPE DEF_MATRIX_ABSTRACT
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC, EXTENDS(DEF_MATRIX_ABSTRACT) :: DEF_MATRIX_1D
     PROCEDURE(matrix_weak_formulation_1D), POINTER :: ptr_matrix_contribution  => NULL ()
     PROCEDURE(rhs_weak_formulation_1D)   , POINTER :: ptr_rhs_contribution     => NULL ()
     PROCEDURE(Assembly_Diagnostics_1D)   , POINTER :: ptr_assembly_diagnostics => NULL ()
  END TYPE DEF_MATRIX_1D
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC :: DEF_MATRIX_1D_CELL
     TYPE(DEF_MATRIX_1D), POINTER  :: ptr_matrix
     TYPE(DEF_MATRIX_1D_CELL), POINTER :: Next=>NULL()
  END TYPE DEF_MATRIX_1D_CELL
  ! ..................................................

  ! ..................................................
  TYPE, PUBLIC :: DEF_MATRIX_FIELD_1D_CELL
     ! ... local field id in the matrix
     INTEGER :: oi_loc_id

     TYPE(DEF_MATRIX_1D), POINTER  :: ptr_matrix
     TYPE(DEF_MATRIX_FIELD_1D_CELL), POINTER :: Next=>NULL()
  END TYPE DEF_MATRIX_FIELD_1D_CELL
  ! ..................................................

  ! ..................................................
  INTERFACE
     SUBROUTINE matrix_weak_formulation_1D(self, ao_gbox)
       USE SPI_GREENBOX_DEF
       IMPORT DEF_MATRIX_1D

       CLASS(DEF_MATRIX_1D)  :: self
       TYPE(DEF_GREENBOX_1D) :: ao_gbox
     END SUBROUTINE matrix_weak_formulation_1D
  END INTERFACE
  ! ..................................................

  ! ..................................................
  INTERFACE
     SUBROUTINE rhs_weak_formulation_1D(self, ao_gbox)
       USE SPI_GREENBOX_DEF
       IMPORT DEF_MATRIX_1D

       CLASS(DEF_MATRIX_1D)  :: self
       TYPE(DEF_GREENBOX_1D) :: ao_gbox
     END SUBROUTINE rhs_weak_formulation_1D
  END INTERFACE
  ! ..................................................

  ! ..................................................
  INTERFACE
     SUBROUTINE Assembly_Diagnostics_1D(self, ao_gbox, ai_nstep)
       USE SPI_GREENBOX_DEF
       IMPORT DEF_MATRIX_1D

       CLASS(DEF_MATRIX_1D)  :: self
       TYPE(DEF_GREENBOX_1D) :: ao_gbox
       INTEGER :: ai_nstep
     END SUBROUTINE Assembly_Diagnostics_1D
  END INTERFACE
  ! ..................................................

END MODULE SPI_MATRIX_DEF
