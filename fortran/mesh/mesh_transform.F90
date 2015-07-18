!# -*- coding: utf8 -*-
MODULE SPI_MESH_TRANSFORM 
  USE SPIMESH_DEF
  USE SPI_GLOBAL_DEF
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: MESH_TRANSLATE_2D

CONTAINS

  ! .........................................................
  SUBROUTINE MESH_TRANSLATE_2D(self, apr_displ)
  IMPLICIT NONE
     TYPE(DEF_MESH_2D), INTENT(INOUT)  :: self
     REAL(KIND=RK), DIMENSION(:)       :: apr_displ
     ! LOCAL
     INTEGER :: li_mesh
     INTEGER :: li_err 

#ifdef DEBUG_TRACE 
     CALL printlog("MESH_TRANSLATE_2D: Begin", ai_dtllevel = 0)
#endif

     CALL JOREK_Param_GETInt(INT_TYPEMESH_ID,li_mesh,li_err)
     SELECT CASE(li_mesh)
     CASE(INT_MESH_CAID_HBEZIER) ! read form CAID files 
        CALL MESH_TRANSLATE_HERMITE_BEZIER_2D(self, apr_displ)
     CASE(INT_MESH_BOXSPLINES) ! hexagonal mesh
        CALL MESH_TRANSLATE_TRIANGLE_2D(self, apr_displ)
     CASE(INT_MESH_BEZIER_DESCRIPTION) ! generic bernstein description 
        CALL MESH_TRANSLATE_BEZIER_2D(self, apr_displ)
     CASE DEFAULT
        PRINT *, "ERROR in MESH_TRANSLATE_2D: Not a recognized type mesh"
        STOP
     END SELECT

#ifdef DEBUG_TRACE 
     CALL printlog("MESH_TRANSLATE_2D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE MESH_TRANSLATE_2D
  ! .........................................................

  ! .........................................................
  SUBROUTINE MESH_TRANSLATE_HERMITE_BEZIER_2D(self, apr_displ)
  IMPLICIT NONE
     TYPE(DEF_MESH_2D), INTENT(INOUT)  :: self
     REAL(KIND=RK), DIMENSION(:)       :: apr_displ
     ! LOCAL
     INTEGER :: li_mesh
     INTEGER :: li_err 

#ifdef DEBUG_TRACE 
     CALL printlog("MESH_TRANSLATE_HERMITE_BEZIER_2D: Begin", ai_dtllevel = 0)
#endif

      self % TheNodes % Coor2D(1,i_D0, :) = self % TheNodes % Coor2D(1, i_D0, :) + apr_displ(1)
     self % TheNodes % Coor2D(2, i_D0, :) = self % TheNodes % Coor2D(2, i_D0, :) + apr_displ(2)

#ifdef DEBUG_TRACE 
     CALL printlog("MESH_TRANSLATE_HERMITE_BEZIER_2D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE MESH_TRANSLATE_HERMITE_BEZIER_2D
  ! .........................................................

  ! .........................................................
  SUBROUTINE MESH_TRANSLATE_BEZIER_2D(self, apr_displ)
  IMPLICIT NONE
     TYPE(DEF_MESH_2D), INTENT(INOUT)  :: self
     REAL(KIND=RK), DIMENSION(:)       :: apr_displ
     ! LOCAL
     INTEGER :: li_mesh
     INTEGER :: li_err 

#ifdef DEBUG_TRACE 
     CALL printlog("MESH_TRANSLATE_BEZIER_2D: Begin", ai_dtllevel = 0)
#endif

     self % TheNodes % Coor2D(1, :, :) = self % TheNodes % Coor2D(1, :, :) + apr_displ(1)
     self % TheNodes % Coor2D(2, :, :) = self % TheNodes % Coor2D(2, :, :) + apr_displ(2)

     self % TheNodes_bezier % Coor2D(1, :, :) = self % TheNodes_bezier % Coor2D(1, :, :) + apr_displ(1)
     self % TheNodes_bezier % Coor2D(2, :, :) = self % TheNodes_bezier % Coor2D(2, :, :) + apr_displ(2)

#ifdef DEBUG_TRACE 
     CALL printlog("MESH_TRANSLATE_BEZIER_2D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE MESH_TRANSLATE_BEZIER_2D
  ! .........................................................

  ! .........................................................
  SUBROUTINE MESH_TRANSLATE_TRIANGLE_2D(self, apr_displ)
  IMPLICIT NONE
     TYPE(DEF_MESH_2D), INTENT(INOUT)  :: self
     REAL(KIND=RK), DIMENSION(:)       :: apr_displ
     ! LOCAL
     INTEGER :: li_mesh
     INTEGER :: li_err 

#ifdef DEBUG_TRACE 
     CALL printlog("MESH_TRANSLATE_TRIANGLE_2D: Begin", ai_dtllevel = 0)
#endif

     PRINT *, "MESH_TRANSLATE_TRIANGLE_2D: TODO"

#ifdef DEBUG_TRACE 
     CALL printlog("MESH_TRANSLATE_TRIANGLE_2D: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE MESH_TRANSLATE_TRIANGLE_2D
  ! .........................................................
END MODULE SPI_MESH_TRANSFORM
