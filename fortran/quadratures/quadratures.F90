!# -*- coding: utf8 -*-
MODULE SPI_QUADRATURES
  USE SPI_QUADRATURES_DEF
  USE SPI_QUADRATURES_LINE

    CONTAINS
  ! .............................................
  SUBROUTINE CREATE_QUADRATURE(self, ai_type, ai_n, ai_n_period, dirname)
  IMPLICIT NONE
     CLASS(DEF_QUADRATURE_ABSTRACT) , INTENT(INOUT) :: self
     INTEGER                        , INTENT(IN)    :: ai_type
     INTEGER                        , INTENT(IN)    :: ai_n 
     INTEGER              , OPTIONAL, INTENT(IN)    :: ai_n_period
     CHARACTER(LEN = 1024), OPTIONAL, INTENT(IN)    :: dirname
     ! LOCAL
     INTEGER :: li_n

     li_n = ai_n

     IF (PRESENT(dirname)) THEN
        self % dirname = TRIM(dirname)
     END IF

     self % oi_type     = ai_type 

     SELECT TYPE (self)
     CLASS IS (DEF_QUADRATURE_1D)
        IF (.NOT. PRESENT(ai_n_period)) THEN
           CALL CREATE_QUADRATURE_LINE(self, li_n)
        ELSE
           CALL CREATE_QUADRATURE_LINE(self, li_n, ai_n_period=ai_n_period)
        END IF
     CLASS DEFAULT
        STOP 'CREATE_QUADRATURE: unexpected type for self object!'
     END SELECT

  END SUBROUTINE CREATE_QUADRATURE
  ! .............................................

  ! .............................................
  SUBROUTINE FREE_QUADRATURE(self)
  IMPLICIT NONE
     CLASS(DEF_QUADRATURE_ABSTRACT) :: self

     DEALLOCATE(self % opr_points )
     DEALLOCATE(self % opr_weights)

  END SUBROUTINE FREE_QUADRATURE
  ! .............................................


END MODULE SPI_QUADRATURES
