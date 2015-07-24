!# -*- coding: utf8 -*-
MODULE spi_numbering
  USE SPI_GLOBAL_DEF 
  USE SPI_NUMBERING_DEF
  USE SPI_MESH_DEF
  IMPLICIT NONE
  
  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 2 

  CONTAINS

  ! .........................................................
  SUBROUTINE CREATE_NUMBERING(self, ao_mesh)
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_NUMBERING_ABSTRACT)  , INTENT(INOUT) :: self 
     CLASS(DEF_MESH_ABSTRACT)       , INTENT(INOUT) :: ao_mesh 
     ! LOCAL

     SELECT TYPE (self)
     CLASS IS (DEF_NUMBERING_1D_BSPLINE)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D_BSPLINE)
           CALL CREATE_NUMBERING_1D_BSPLINE(self, ao_mesh)
        CLASS DEFAULT
           STOP 'CREATE_NUMBERING: unexpected type for ao_mesh object!'
        END SELECT

     CLASS IS (DEF_NUMBERING_1D_FOURIER)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D_FOURIER)
           CALL CREATE_NUMBERING_1D_FOURIER(self, ao_mesh)
        CLASS DEFAULT
           STOP 'CREATE_NUMBERING: unexpected type for ao_mesh object!'
        END SELECT

     CLASS IS (DEF_NUMBERING_1D_HBEZIER)
        SELECT TYPE (ao_mesh)
        CLASS IS (DEF_MESH_1D_HBEZIER)
           CALL CREATE_NUMBERING_1D_HBEZIER(self, ao_mesh)
        CLASS DEFAULT
           STOP 'CREATE_NUMBERING: unexpected type for ao_mesh object!'
        END SELECT

     CLASS DEFAULT
        STOP 'CREATE_NUMBERING: unexpected type for self object!'
     END SELECT

  END SUBROUTINE CREATE_NUMBERING
  ! .........................................................

  ! ........................................................
  SUBROUTINE CREATE_NUMBERING_1D_BSPLINE( self, ao_mesh)
  implicit none
     type(DEF_NUMBERING_1D_BSPLINE),intent(inout) :: self
     TYPE(DEF_MESH_1D_BSPLINE), TARGET, INTENT(IN) :: ao_mesh
     ! LOCAL VARIABLES
     integer  :: li_ref

     self % ptr_mesh => ao_mesh
     
     ALLOCATE ( self%opi_IEN( ao_mesh%oi_nen, ao_mesh%oi_n_elmts ) )
     ALLOCATE ( self%opi_LM ( ao_mesh%oi_nen, ao_mesh%oi_n_elmts ) )		
     ALLOCATE ( self%opi_ID ( ao_mesh%oi_nnp ) )		 
     
     ALLOCATE ( self%opi_ELT_INDEX ( ao_mesh % oi_n_elmts ) ) 

     ! INITIALIZING THE ID ARRAY		
     call INIT_ID_BSPLINE(self, ao_mesh % oi_type_bc)	
     call INIT_IEN_BSPLINE(self)
     call INIT_LM_BSPLINE(self)
  
  END SUBROUTINE CREATE_NUMBERING_1D_BSPLINE
  ! ........................................................

  ! .........................................................
  SUBROUTINE CREATE_NUMBERING_1D_FOURIER(self, ao_mesh)
  IMPLICIT NONE
     TYPE(DEF_NUMBERING_1D_FOURIER), INTENT(INOUT)  :: self
     TYPE(DEF_MESH_1D_FOURIER), TARGET, INTENT(IN) :: ao_mesh
     ! LOCAL
     INTEGER :: li_err 

  END SUBROUTINE CREATE_NUMBERING_1D_FOURIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE CREATE_NUMBERING_1D_HBEZIER(self, ao_mesh)
  IMPLICIT NONE
     TYPE(DEF_NUMBERING_1D_HBEZIER), INTENT(INOUT)  :: self
     TYPE(DEF_MESH_1D_HBEZIER), TARGET, INTENT(IN) :: ao_mesh
     ! LOCAL
     INTEGER :: li_err 

  END SUBROUTINE CREATE_NUMBERING_1D_HBEZIER
  ! .........................................................

  ! .........................................................
  SUBROUTINE FREE_NUMBERING(self)
  !     dirname is the directory where geometry files are given
  !     if not provided, then dirname is = $PWD
  IMPLICIT NONE
     CLASS(DEF_NUMBERING_ABSTRACT)       , INTENT(INOUT) :: self 
     ! LOCAL

      DEALLOCATE ( self%opi_IEN )
      DEALLOCATE ( self%opi_ID  )
      DEALLOCATE ( self%opi_LM  )		

      SELECT TYPE (self)
      CLASS IS (DEF_NUMBERING_1D_BSPLINE)
         CALL FREE_NUMBERING_1D_BSPLINE(self)
     
      CLASS IS (DEF_NUMBERING_1D_FOURIER)
         CALL FREE_NUMBERING_1D_FOURIER(self)
     
      CLASS IS (DEF_NUMBERING_1D_HBEZIER)
         CALL FREE_NUMBERING_1D_HBEZIER(self)
     
      CLASS DEFAULT
         STOP 'FREE_NUMBERING: unexpected type for self object!'
      END SELECT

  END SUBROUTINE FREE_NUMBERING
  ! .........................................................
         
   ! ........................................................
   SUBROUTINE FREE_NUMBERING_1D_BSPLINE( self )
   implicit none
      type(DEF_NUMBERING_1D_BSPLINE),intent(inout) :: self
      ! LOCAL VARIABLES
      
   END SUBROUTINE FREE_NUMBERING_1D_BSPLINE
   ! ........................................................

   ! ........................................................
   SUBROUTINE FREE_NUMBERING_1D_FOURIER( self )
   implicit none
      type(DEF_NUMBERING_1D_FOURIER),intent(inout) :: self
      ! LOCAL VARIABLES
      
   END SUBROUTINE FREE_NUMBERING_1D_FOURIER
   ! ........................................................

   ! ........................................................
   SUBROUTINE FREE_NUMBERING_1D_HBEZIER( self )
   implicit none
      type(DEF_NUMBERING_1D_HBEZIER),intent(inout) :: self
      ! LOCAL VARIABLES
      
   END SUBROUTINE FREE_NUMBERING_1D_HBEZIER
   ! ........................................................

   ! ........................................................
   SUBROUTINE INIT_ID_BSPLINE( self, ai_bc_type )
   implicit none
      type(DEF_NUMBERING_1D_BSPLINE),intent(inout) :: self
      integer  :: ai_bc_type
      
      self % opi_ID = 0
      
      select case ( ai_bc_type )
      case ( SPI_BC_NONE ) 
         call INIT_ID_BC_NONE_BSPLINE ( self )
      case ( SPI_BC_DIRICHLET_HOMOGEN ) 
         call INIT_ID_BC_DIRICHLET_HOMOGEN_BSPLINE ( self )
      case ( SPI_BC_PERIODIC )
         call INIT_ID_BC_PERIODIC_BSPLINE ( self )
      end select	

   END SUBROUTINE INIT_ID_BSPLINE 
   ! ........................................................
      
   ! ........................................................
   SUBROUTINE INIT_ID_BC_NONE_BSPLINE(self)
   implicit none
      type(DEF_NUMBERING_1D_BSPLINE),intent(inout) :: self
      !LOCAL VARIABLES
      integer  :: li_d
      integer  :: li_i, li_A
      
      li_d = 0
      li_A = 0	
   
      do li_i = 1, self % ptr_mesh % oi_n
         li_A = li_A + 1
         
         li_d = li_d + 1
         
         self%opi_ID ( li_A ) = li_d
      end do
   
      self%oi_size = li_d

   END SUBROUTINE INIT_ID_BC_NONE_BSPLINE 
   ! ........................................................

   ! ........................................................
   SUBROUTINE INIT_ID_BC_DIRICHLET_HOMOGEN_BSPLINE(self)
   implicit none
      type(DEF_NUMBERING_1D_BSPLINE),intent(inout) :: self
      !LOCAL VARIABLES
      integer  :: li_d
      integer  :: li_i, li_A
      
      li_d = 0
      li_A = 0	
      
      do li_i = 1, self % ptr_mesh % oi_n
         li_A = li_A + 1
         
         if ( ( li_i == 1 ) .OR. ( li_i == self % ptr_mesh % oi_n ) ) then
            self%opi_ID ( li_A ) = 0				
         else
            li_d = li_d + 1
            
            self%opi_ID ( li_A ) = li_d
         end if		
      end do
      
      self%oi_size = li_d
      
    END SUBROUTINE INIT_ID_BC_DIRICHLET_HOMOGEN_BSPLINE 
    ! ........................................................

    ! ........................................................
   SUBROUTINE INIT_ID_BC_PERIODIC_BSPLINE(self)
   implicit none
      type(DEF_NUMBERING_1D_BSPLINE),intent(inout) :: self
      !LOCAL VARIABLES
      integer  :: li_d
      integer  :: li_nu		
      integer  :: li_i
      integer  :: li_A		
      integer  :: li_L		
      
      li_A = 0
      li_d = 0		
   
      li_nu = self % ptr_mesh % oi_p

      do li_i = 1, self % ptr_mesh % oi_n			
         li_A = li_A + 1
         
         if ( li_i <= self % ptr_mesh % oi_n - li_nu ) then	
            if ( self%opi_ID ( li_A ) == 0 ) then
               li_d = li_d + 1
               
               self%opi_ID ( li_A ) = li_d
            end if		
            
            if ( ( 1 <= li_i ) .AND. ( li_i <= li_nu ) ) then	
               li_L = self % ptr_mesh % oi_n - li_nu
               
               self%opi_ID  ( li_A + li_L ) = self%opi_ID ( li_A )
            end if
         end if		
      end do

      self%oi_size = li_d

   END SUBROUTINE INIT_ID_BC_PERIODIC_BSPLINE 
   ! ........................................................

   ! ........................................................
   SUBROUTINE INIT_IEN_BSPLINE(self)
   implicit none
      type(DEF_NUMBERING_1D_BSPLINE),intent(inout) :: self
      !LOCAL VARIABLES
      integer  :: li_A
      integer  :: li_B
      integer  :: li_Bloc		
      integer  :: li_e
      integer  :: li_eloc
      integer  :: li_i
      integer  :: li_iloc	
      integer  :: li_mini
      integer  :: li_maxi		
      
      !... computing the max/min for index i, depending on the type of the patch
      li_mini = self % ptr_mesh % oi_p + 1
      li_maxi = self % ptr_mesh % oi_n		
      if ( ( self % ptr_mesh % oi_type_bc /= SPI_BC_PERIODIC )          .AND.	&
           ( self % ptr_mesh % oi_type_bc /= SPI_BC_DIRICHLET_HOMOGEN ) .AND.	&
           ( self % ptr_mesh % oi_type_bc /= SPI_BC_NONE ) ) then

         print*,'init_IEN_singlefem : Not implemented yet'
         STOP
      end if 
      
      li_eloc = 0
      do li_i = li_mini, li_maxi
         ! WE CHECK IF THE ELEMENT HAS ZERO MEASURE
         if ( self % ptr_mesh % opr_knot ( li_i ) /= self % ptr_mesh % opr_knot ( li_i + 1 ) ) then
            li_eloc = li_eloc + 1
            
            self%opi_ELT_INDEX ( li_eloc ) = li_i
         end if	
      end do		
      
      self%opi_IEN = 0	
      do li_e = 1 , self % ptr_mesh % oi_n_elmts
         li_i = self%opi_ELT_INDEX ( li_e )	
         
         li_A = li_i							
         do li_iloc = 0, self % ptr_mesh % oi_p
            li_B = li_A - ( self % ptr_mesh % oi_p - li_iloc )
                
            li_Bloc = li_iloc + 1
                        
            self%opi_IEN ( li_Bloc , li_e ) = li_B 
         end do			
      end do
      
   END SUBROUTINE INIT_IEN_BSPLINE 
   ! ........................................................
   
   ! ........................................................
   SUBROUTINE INIT_LM_BSPLINE(self)
   implicit none
      type(DEF_NUMBERING_1D_BSPLINE),intent(inout) :: self
      !LOCAL VARIABLES
      integer  :: li_e, li_b, li_dof
      
      do li_e = 1 , self % ptr_mesh % oi_n_elmts
         do li_b = 1, self % ptr_mesh % oi_nen
            self%opi_LM ( li_b , li_e ) = self%opi_ID ( self%opi_IEN ( li_b , li_e ) ) 
         end do
      end do
      
   END SUBROUTINE INIT_LM_BSPLINE 
   ! ........................................................
   
   ! ........................................................
  
END MODULE spi_numbering 
