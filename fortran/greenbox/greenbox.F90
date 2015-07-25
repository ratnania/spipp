!# -*- coding: utf8 -*-
MODULE SPI_GREENBOX
  USE SPI_GREENBOX_DEF
  USE SPI_QUADRATURES_DEF
  USE SPI_BLACKBOX_DEF

    CONTAINS

  ! ..........................................................        
  SUBROUTINE CREATE_GREENBOX(self, ai_n_var, ao_quad)
  IMPLICIT NONE
     CLASS(DEF_GREENBOX_ABSTRACT)         :: self
     CLASS(DEF_QUADRATURE_1D)    , TARGET :: ao_quad
     INTEGER ::  ai_n_var
     ! LOCAL
     INTEGER :: ierr
     
     self % oi_nvar = ai_n_var

     SELECT TYPE (self)
     CLASS IS (DEF_GREENBOX_1D)
        CALL CREATE_GREENBOX_1D(self, ao_quad)
     CLASS DEFAULT
        STOP 'CREATE_GREENBOX: unexpected type for self object!'
     END SELECT

  END SUBROUTINE CREATE_GREENBOX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE RESET_GREENBOX(self)
  IMPLICIT NONE
     CLASS(DEF_GREENBOX_ABSTRACT) :: self
     ! LOCAL
     
     SELECT TYPE (self)
     CLASS IS (DEF_GREENBOX_1D)
        CALL RESET_GREENBOX_1D(self)
     CLASS DEFAULT
        STOP 'RESET_GREENBOX: unexpected type for self object!'
     END SELECT

  END SUBROUTINE RESET_GREENBOX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE FREE_GREENBOX(self)
  IMPLICIT NONE
     CLASS(DEF_GREENBOX_ABSTRACT) :: self
     ! LOCAL

     SELECT TYPE (self)
     CLASS IS (DEF_GREENBOX_1D)
        CALL FREE_GREENBOX_1D(self)
     CLASS DEFAULT
        STOP 'FREE_GREENBOX: unexpected type for self object!'
     END SELECT

  END SUBROUTINE FREE_GREENBOX
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE CREATE_GREENBOX_1D(self, ao_quad)
  IMPLICIT NONE
     CLASS(DEF_GREENBOX_1D)           :: self
     CLASS(DEF_QUADRATURE_1D), TARGET :: ao_quad
     ! LOCAL
     INTEGER :: ierr
     INTEGER, PARAMETER :: N_DIM = 1 

     self % ptr_quad => ao_quad

     ALLOCATE(self % VarN_0   (self % oi_nvar, self % ptr_quad % oi_n_points))
     
     ALLOCATE(self % VarN_s1  (self % oi_nvar, self % ptr_quad % oi_n_points))
     ALLOCATE(self % VarN_s1s1(self % oi_nvar, self % ptr_quad % oi_n_points))

     ALLOCATE(self % VarN_x1  (self % oi_nvar, self % ptr_quad % oi_n_points))
     ALLOCATE(self % VarN_x1x1(self % oi_nvar, self % ptr_quad % oi_n_points))

     ALLOCATE(self % Sol_analytical(self % ptr_quad % oi_n_points, self % oi_nvar, N_DIM+1))

  END SUBROUTINE CREATE_GREENBOX_1D
  ! .......................................................... 

  ! ..........................................................        
  SUBROUTINE RESET_GREENBOX_1D(self)
  IMPLICIT NONE
     CLASS(DEF_GREENBOX_1D) :: self
     ! LOCAL

     self % VarN_0  = 0.0

     self % VarN_s1   = 0.0 
     self % VarN_s1s1 = 0.0

     self % VarN_x1   = 0.0 
     self % VarN_x1x1 = 0.0 

  END SUBROUTINE RESET_GREENBOX_1D
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE FREE_GREENBOX_1D(self)
  IMPLICIT NONE
     CLASS(DEF_GREENBOX_1D) :: self
     ! LOCAL

  END SUBROUTINE FREE_GREENBOX_1D
  ! ..........................................................  

  ! ..........................................................        
  SUBROUTINE UPDATE_VARIABLES_GREENBOX_1D(self, ao_bbox, ivar, ar_value, ai_i)
  IMPLICIT NONE
     CLASS(DEF_GREENBOX_1D) :: self
     CLASS(DEF_BLACKBOX_1D) :: ao_bbox
     INTEGER :: ivar
     REAL(SPI_RK) :: ar_value 
     INTEGER :: ai_i
     ! LOCAL

     self % VarN_0    (ivar, ai_i) = self % VarN_0    (ivar, ai_i) + ao_bbox%B_0(ai_i)    * ar_value  
     self % VarN_s1   (ivar, ai_i) = self % VarN_s1   (ivar, ai_i) + ao_bbox%B_s1(ai_i)   * ar_value  
     self % VarN_s1s1 (ivar, ai_i) = self % VarN_s1s1 (ivar, ai_i) + ao_bbox%B_s1s1(ai_i) * ar_value

  END SUBROUTINE UPDATE_VARIABLES_GREENBOX_1D
  ! ..........................................................  

END MODULE SPI_GREENBOX
