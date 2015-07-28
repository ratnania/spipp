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
     INTEGER :: oi_nvar
     INTEGER :: oi_n_elmts
     INTEGER :: oi_n_rows 
     INTEGER :: oi_n_cols
     INTEGER :: oi_nen_1
     INTEGER :: oi_nen_2
     LOGICAL :: ol_to_assembly
     LOGICAL :: ol_allocated_csrmatrix

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

     LOGICAL :: ol_allocated_contribution_matrix
     REAL(KIND=SPI_RK), DIMENSION(:,:), POINTER :: Matrix_Contribution 

     LOGICAL :: ol_allocated_contribution_rhs
     REAL(KIND=SPI_RK), DIMENSION(:)  , POINTER :: Rhs_Contribution

     INTEGER, DIMENSION(:, :), POINTER :: ptr_LM_1
     INTEGER, DIMENSION(:, :), POINTER :: ptr_LM_2

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
