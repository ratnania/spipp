!
! File: ex1.f90
!
!
! Usage:
!   > ./ex1 
!
! Authors:
!   Ahmed RATNANI  - ratnaniahmed@gmail.com
!
   
! ............................................
PROGRAM Main
USE SPI_GLOBAL_DEF
USE SPI_GLOBAL
USE SPI_QUADRATURES_DEF
USE SPI_QUADRATURES
USE SPI_MESH_DEF
USE SPI_MESH
USE SPI_BASIS_DEF
USE SPI_BASIS
USE SPI_BLACKBOX_DEF
USE SPI_BLACKBOX
USE SPI_GREENBOX_DEF
USE SPI_GREENBOX
USE SPI_NUMBERING_DEF
USE SPI_NUMBERING
USE SPI_SPARSE_MATRIX_DEF
USE SPI_SPARSE_MATRIX
USE SPI_MATRIX_DEF
USE SPI_MATRIX
IMPLICIT NONE

  CALL test1()

CONTAINS
   ! ..................................................
   FUNCTION ANALYTICAL_RHS(ao_bbox)
   IMPLICIT NONE
      REAL(KIND=SPI_RK) :: ANALYTICAL_RHS
      TYPE(DEF_BLACKBOX_1D) :: ao_bbox
      REAL(KIND=SPI_RK) :: lr_R
      ! LOCAL
      REAL(KIND=SPI_RK) :: lr_k1
      REAL(KIND=SPI_RK) :: lr_k2
      INTEGER       :: ierr
      REAL(KIND=SPI_RK) :: lr_R0
      REAL(KIND=SPI_RK) :: lr_a
      REAL(KIND=SPI_RK) :: lr_acenter
      INTEGER       :: ijg
      INTEGER, PARAMETER :: li_mode_m1 = 1 
   
      lr_R   = ao_bbox%Xp_0(1,ijg)
     
      lr_k1 = 2.0 * SPI_PI * FLOAT(li_mode_m1)
     
      ANALYTICAL_RHS = lr_k1**2 * SIN(lr_k1*lr_R)
   
   END FUNCTION ANALYTICAL_RHS
   ! ..................................................
   
   ! ............................................
   SUBROUTINE  RHS_for_Vi(ptr_matrix, ao_bboxi, ao_gbox)
   IMPLICIT NONE
      CLASS(DEF_MATRIX_1D), POINTER :: ptr_matrix
      TYPE(DEF_BLACKBOX_1D) :: ao_bboxi
      TYPE(DEF_GREENBOX_1D)                :: ao_gbox
      ! LOCAL
      REAL(KIND=SPI_RK) :: f_rhs
      REAL(KIND=SPI_RK) :: contribution
      INTEGER       :: ijg
      REAL(KIND=SPI_RK) :: wVol
      REAL(KIND=SPI_RK) :: Vi_0
      REAL(KIND=SPI_RK) :: Vi_R
      REAL(KIND=SPI_RK) :: Vi_RR
      REAL(KIND=SPI_RK), DIMENSION(:), POINTER :: Rhs_Contribution
   
      Rhs_Contribution => ptr_matrix % Rhs_Contribution
   
      wVol  = ao_bboxi % wVol(ijg)
      Vi_0  = ao_bboxi % B_0(ijg)
      Vi_R  = ao_bboxi % B_x1(ijg)
      Vi_RR = ao_bboxi % B_x1x1(ijg)
   
      ! ... ADD L2 contribution
      f_rhs       = ANALYTICAL_RHS(ao_bboxi)
   
      contribution                = Vi_0 *wVol*f_rhs
      Rhs_Contribution(1)  =  contribution 
      ! ...
   
   END SUBROUTINE RHS_for_Vi
   ! ............................................
   
   ! ............................................
   SUBROUTINE Matrix_for_Vi_Vj(ptr_matrix, ao_bboxi, ao_bboxj, ao_gbox)
   IMPLICIT NONE
      CLASS(DEF_MATRIX_1D), POINTER :: ptr_matrix
      TYPE(DEF_BLACKBOX_1D) :: ao_bboxi
      TYPE(DEF_BLACKBOX_1D) :: ao_bboxj
      TYPE(DEF_GREENBOX_1D)                :: ao_gbox
      ! LOCAL
      REAL(KIND=SPI_RK) :: contribution
      INTEGER       :: ijg
      REAL(KIND=SPI_RK) :: wVol
      REAL(KIND=SPI_RK) :: Vi_0
      REAL(KIND=SPI_RK) :: Vi_R
      REAL(KIND=SPI_RK) :: Vi_RR
      REAL(KIND=SPI_RK) :: Vj_0
      REAL(KIND=SPI_RK) :: Vj_R
      REAL(KIND=SPI_RK) :: Vj_RR
      REAL(KIND=SPI_RK), DIMENSION(:,:), POINTER :: Matrix_Contribution

      Matrix_Contribution => ptr_matrix % Matrix_Contribution
  
      wVol  = ao_bboxi % wVol(ijg)
      Vi_0  = ao_bboxi % B_0(ijg)    ; Vj_0  = ao_bboxj % B_0(ijg)
      Vi_R  = ao_bboxi % B_x1(ijg)   ; Vj_R  = ao_bboxj % B_x1(ijg)
      Vi_RR = ao_bboxi % B_x1x1(ijg) ; Vj_RR = ao_bboxj % B_x1x1(ijg)
   
      ! ... ADD STIFFNESS CONTRIBUTION
      contribution =  Vi_R * Vj_R * wVol
   
      Matrix_Contribution(1, 1) =  contribution 
      ! ...
   
   END SUBROUTINE Matrix_for_Vi_Vj
   ! ............................................

   ! ............................................
   subroutine test1 ()
   implicit none
      ! LOCAL
      TYPE(DEF_QUADRATURE_1D), TARGET :: lo_quad
      TYPE(DEF_MESH_1D_BSPLINE), TARGET :: lo_mesh
      TYPE(DEF_BASIS_1D_BSPLINE), TARGET :: lo_basis
      TYPE(DEF_NUMBERING_1D_BSPLINE), TARGET :: lo_numbering
      TYPE(DEF_BLACKBOX_1D_BSPLINE), TARGET :: lo_bbox
      TYPE(DEF_GREENBOX_1D), TARGET :: lo_gbox
      TYPE(DEF_MATRIX_1D), TARGET :: lo_matrix
      ! ... number of internal knots is = N - P - 1
      INTEGER, PARAMETER :: N = 5 
      INTEGER, PARAMETER :: P = 3
      INTEGER, PARAMETER :: K = 3 
      INTEGER, PARAMETER :: N_VAR = 1 
      INTEGER, PARAMETER :: I_VAR = 1 
      INTEGER :: li_elmt
      INTEGER :: li_i
      INTEGER :: li_n_rows
      INTEGER :: li_n_cols
      INTEGER :: li_glob
      INTEGER :: li_index
      REAL(SPI_RK) :: lr_a
      REAL(SPI_RK) :: lr_b
      REAL(SPI_RK) :: lr_value
   
      CALL CREATE_QUADRATURE(lo_quad, SPI_QUADRATURES_LEGENDRE, K)
      CALL CREATE_MESH(lo_mesh, lo_quad, ai_n=N, ai_p=P, ai_type_bc=SPI_BC_PERIODIC) 
      CALL CREATE_BASIS(lo_basis, lo_mesh) 
      CALL CREATE_NUMBERING(lo_numbering, lo_mesh) 
      CALL CREATE_BLACKBOX(lo_bbox, lo_basis, lo_quad)
      CALL CREATE_GREENBOX(lo_gbox, N_VAR, lo_quad)
      CALL CREATE_MATRIX(lo_matrix)
   
      PRINT *, ">>> knots"
      PRINT *, lo_mesh % opr_knot
   
      CALL GET_NR_MATRIX(lo_matrix, li_n_rows) 
      CALL GET_NR_MATRIX(lo_matrix, li_n_cols) 
   
      CALL GLOBAL_TO_INDEX_MATRIX(lo_matrix, li_glob, li_index) 
   
      lo_matrix % ptr_matrix_contribution => Matrix_for_Vi_Vj
      lo_matrix % ptr_rhs_contribution    => RHS_for_Vi
   
      CALL RESET_BASIS(lo_basis) 
      print *, "%%%%"
      CALL RESET_GREENBOX(lo_gbox) 
      print *, "%%%%"
  
      CALL BLACKBOX_RESET_POSITION(lo_bbox) 
      print *, "%%%%"
      CALL RESET_ELEMENT_MATRIX(lo_matrix) 
      print *, "%%%%"
      CALL RESET_ELEMENT_RHS_MATRIX(lo_matrix) 
      print *, "%%%%"
   
      li_elmt = 1
      CALL UPDATE_BASIS(lo_basis, li_elmt) 
      print *, "%%%%"
   
      lr_a = 0.0 ; lr_b = 1.0
      CALL BLACKBOX_COMPUTE_METRIC(lo_bbox, lr_a, lr_b) 
      print *, "%%%%"
   
      li_i = 1
      lr_value = 1.0
      CALL UPDATE_VARIABLES_GREENBOX_1D(lo_gbox, lo_bbox, I_VAR, lr_value, li_i) 
      print *, "%%%%"

      CALL FREE_QUADRATURE(lo_quad) 
      CALL FREE_MESH(lo_mesh) 
      CALL FREE_BASIS(lo_basis) 
      CALL FREE_NUMBERING(lo_numbering) 
      CALL FREE_BLACKBOX(lo_bbox) 
      CALL FREE_GREENBOX(lo_gbox) 
      CALL FREE_MATRIX(lo_matrix) 

   end subroutine test1
   ! ............................................

END PROGRAM Main
! ............................................
   
