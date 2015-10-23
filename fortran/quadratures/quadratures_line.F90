!# -*- coding: utf8 -*-
MODULE SPI_QUADRATURES_LINE 
  USE SPI_QUADRATURES_DEF

    CONTAINS
  ! .............................................
  SUBROUTINE CREATE_QUADRATURE_LINE(self, ai_n, ai_n_period)
  IMPLICIT NONE
     CLASS(DEF_QUADRATURE_ABSTRACT) :: self
     INTEGER :: ai_n 
     INTEGER, OPTIONAL :: ai_n_period
     ! LOCAL
     REAL(KIND=SPI_RK), DIMENSION(:), ALLOCATABLE :: GL
     REAL(KIND=SPI_RK), DIMENSION(:), ALLOCATABLE :: W
     INTEGER :: li_i
     INTEGER :: li_size

#ifdef DEBUG_TRACE 
     CALL printlog("CREATE_QUADRATURE_LINE: Begin", ai_dtllevel = 0)
#endif

     li_size = ai_n + 1
     ! FOURIER. n_period IS EXPECTED 
     IF (self % oi_type .EQ. SPI_QUADRATURES_FOURIER) THEN
        li_size = ai_n
     END IF

     self % oi_n_points = li_size

     ALLOCATE(self % points (1:self % oi_n_points))
     ALLOCATE(self % weights(self % oi_n_points))

     ALLOCATE(GL(li_size))
     ALLOCATE(W(li_size))

     SELECT CASE(self % oi_type)
        ! ... DEFAULT JOREK QUADRATURES RULE
        CASE(SPI_QUADRATURES_DEFAULT) 
           CALL GaussJorek(ai_n,GL,W)
          
           DO li_i = 1, li_size 
              self % points (li_i) = GL(li_i) 
              self % weights(li_i) = W(li_i)
           END DO

        ! ... LEGENDRE
        CASE(SPI_QUADRATURES_LEGENDRE) 
           CALL GaussLegendre(ai_n,GL,W)
          
           DO li_i = 1, li_size 
              self % points (li_i) = 0.5 * (GL(li_i) + 1.0) 
              self % weights(li_i) = 0.5 * W(li_i)
           END DO

        ! ... LOBATTO
        CASE(SPI_QUADRATURES_LOBATTO) 
           CALL GaussLobatto(ai_n,GL,W)
          
           DO li_i = 1, li_size 
              self % points (li_i) = 0.5 * (GL(li_i) + 1.0) 
              self % weights(li_i) = 0.5 * W(li_i)
           END DO

        ! ... FOURIER. n_period IS EXPECTED 
        CASE(SPI_QUADRATURES_FOURIER) 
           IF (.NOT. PRESENT(ai_n_period)) THEN
              PRINT *, "ai_n_period is expected for Fourier quadrature rule"
              STOP
           END IF
          
           CALL FourierRule(ai_n, ai_n_period, GL,W)
          
           DO li_i = 1, li_size 
              self % points (li_i) = GL(li_i)
              self % weights(li_i) = W(li_i)
           END DO

        ! ... Default
        CASE DEFAULT
           STOP 'CREATE_QUADRATURE_LINE: unexpected type for quadrature!'
     END SELECT

#ifdef DEBUG_TRACE 
     CALL printlog("CREATE_QUADRATURE_LINE: End", ai_dtllevel = 0)
#endif

  END SUBROUTINE CREATE_QUADRATURE_LINE
  ! .............................................

  ! .............................................
SUBROUTINE FourierRule(n_plane, n_period, GL,W) 
    implicit none
    integer :: n_plane
    integer :: n_period
    real(SPI_RK), dimension(n_plane) :: GL, W
    ! LOCAL
    INTEGER :: ig
   
    DO ig=1,n_plane
        GL(ig) = 2.d0*SPI_PI*float(ig-1)/float(n_plane) / float(n_period)
        W(ig)  = 2.d0*SPI_PI/float(n_plane)/ float(n_period)
    END DO
   
END SUBROUTINE FourierRule
  ! .............................................

  ! .............................................
SUBROUTINE GaussJorek(k,GL,W) 
    implicit none
    integer :: k
        real(SPI_RK), dimension(k+1)			:: GL, W
   
   IF ( k == 3 ) THEN
                !---------	 	   	   
                GL(1) = 0.0694318442029735d0
                GL(2) = 0.3300094782075720d0
                GL(3) = 0.6699905217924280d0
                GL(4) = 0.9305681557970265d0
                !---------	 	   			 
                W(1) = 0.173927422568727d0 			
                W(2) = 0.326072577431273d0						
                W(3) = 0.326072577431273d0						
                W(4) = 0.173927422568727d0	
                !---------	 	   								
   ELSEIF ( k == 4 ) THEN
                !---------	 	   	   
                GL(1) = 0.046910077030668d0
                GL(2) = 0.230765344947158d0
                GL(3) = 0.5d0
                GL(4) = 0.769234655052841d0			 
                GL(5) = 0.953089922969332d0
                !---------	 	   			 
                W(1) = 0.118463442528095d0    	
                W(2) = 0.239314335249683d0 				
                W(3) = 0.284444444444444d0    							
                W(4) = 0.239314335249683d0     							
                W(5) = 0.118463442528095d0
                !---------	 	
   ELSEIF ( k == 5 ) THEN
                !---------	 	   	   
                GL(1) = 0.033765242898424d0
                GL(2) = 0.169395306766868d0
                GL(3) = 0.380690406958402d0
                GL(4) = 0.619309593041598d0		 
                GL(5) = 0.830604693233132d0
                GL(6) = 0.966234757101576d0
                !---------	 	   			 
                W(1)  = 0.085662246189585d0    	
                W(2)  = 0.180380786524069d0		
                W(3)  = 0.233956967286345d0						
                W(4)  = 0.233956967286345d0							
                W(5)  = 0.180380786524069d0
                W(6)  = 0.085662246189585d0
                !---------	
   ELSE
        PRINT *, "Quadrature not available for the given degree ", k
        STOP
   END if
   
END SUBROUTINE GaussJorek
  ! .............................................

  ! .............................................
SUBROUTINE GaussLobatto(k,GL,W) 
    implicit none
    integer :: k
        real(SPI_RK), dimension(k+1)			:: GL, W
   
   if ( k == 1) THEN
                !---------	 	   
                GL(1) = -1.0_SPI_RK
                GL(2) = 1.0_SPI_RK	
                !---------	 
                W(1) = 1.0_SPI_RK			
                W(2) = 1.0_SPI_RK
                !---------	 	   			
ELSEif ( k == 2 ) THEN
                !---------	 	   	   
                GL(1) = -1.0_SPI_RK
                GL(2) = 0.0_SPI_RK
                GL(3) = 1.0_SPI_RK		
                !---------	 	   			 
                W(1) = 0.33333333333333_SPI_RK!1_SPI_RK/3_SPI_RK			
                W(2) = 1.33333333333333_SPI_RK!4_SPI_RK/3_SPI_RK						
                W(3) = 0.33333333333333_SPI_RK!1_SPI_RK/3_SPI_RK	
                !---------	 	   								
ELSEif ( k == 3 ) THEN
                !---------	 	   	   
                GL(1) = -1.0_SPI_RK
                GL(2) = -1.0_SPI_RK/sqrt(5.d0)
                GL(3) = 1.0_SPI_RK/sqrt(5.d0)
                GL(4) = 1.0_SPI_RK		
                !---------	 	   			 
                W(1) = 1.0_SPI_RK/6.0_SPI_RK			
                W(2) = 5.0_SPI_RK/6.0_SPI_RK						
                W(3) = 5.0_SPI_RK/6.0_SPI_RK						
                W(4) = 1.0_SPI_RK/6.0_SPI_RK	
                !---------	 	   								
ELSEif ( k == 4 ) THEN
                !---------	 	   	   
                GL(1) = -1.0_SPI_RK
                GL(2) = -sqrt(3.d0/7.d0)*1.0_SPI_RK
                GL(3) = 0.0_SPI_RK		 
                GL(4) = sqrt(3.d0/7.d0)*1.0_SPI_RK				 
                GL(5) = 1.0_SPI_RK		
                !---------	 	   			 
                W(1) = 1.0_SPI_RK/10.0_SPI_RK			
                W(2) = 49.0_SPI_RK/90.0_SPI_RK						
                W(3) = 32.0_SPI_RK/45.0_SPI_RK									
                W(4) = 49.0_SPI_RK/90.0_SPI_RK									
                W(5) = 1.0_SPI_RK/10.0_SPI_RK	
                !---------	 	   											
   ELSE 
                !calcul des points et poids de Gauss Lobatto
                call ptpdgl(k+1,GL,W)
   END if
   
END SUBROUTINE GaussLobatto
  ! .............................................

  ! .............................................

SUBROUTINE GaussLegendre(k,Node,W) 
! k MUST BE EQUAL TO n-1, k = n-1
implicit none
	integer :: k
	real(SPI_RK), dimension(k+1)			:: Node, W
	! LOCAL VARIABLES
	integer  :: li_ios, li_flag, li_err		
	integer  :: i
	real(SPI_RK) :: lr_w, lr_x
	integer, parameter :: li_file_node = 1
	integer, parameter :: li_file_weight = 2	
	
	Node = 0.0_SPI_RK
	W = 0.0_SPI_RK	

	if ( k == 1 ) THEN
		!---------	 	   
		Node(1) = -0.577350269189626_SPI_RK			
		Node(2) = 0.577350269189626_SPI_RK			
		!---------	 
		W(1) = 1.0_SPI_RK		
		W(2) = 1.0_SPI_RK			
		!---------	 	   				
	ELSEif ( k == 2 ) THEN
		!---------	 	   
		Node(1) = -0.774596669241483_SPI_RK	
		Node(2) = 0.0_SPI_RK 
		Node(3) = 0.774596669241483_SPI_RK	
		!---------	 
		W(1) = 0.555555555555556_SPI_RK	
		W(2) = 0.888888888888889_SPI_RK			
		W(3) = 0.555555555555556_SPI_RK	
		!---------	 	   			
	ELSEif ( k == 3 ) THEN
		!---------	 	   	   
		Node(1) = -0.861136311594053_SPI_RK	
		Node(2) = -0.339981043584856_SPI_RK
		Node(3) = 0.339981043584856_SPI_RK
		Node(4) = 0.861136311594053_SPI_RK	
		!---------	
		W(1) = 0.347854845137454_SPI_RK	 	   			 
		W(2) = 0.652145154862546_SPI_RK			
		W(3) = 0.652145154862546_SPI_RK			
		W(4) = 0.347854845137454_SPI_RK	 	   			 
		!---------	 	   								
	ELSEif ( k == 4 ) THEN
		!---------	 	
		Node(1) = -0.906179845938664_SPI_RK			   	   
		Node(2) = -0.538469310105683_SPI_RK
		Node(3) = 0.0_SPI_RK		 
		Node(4) = 0.538469310105683_SPI_RK
		Node(5) = 0.906179845938664_SPI_RK			   	   
		!---------	 
		W(1) = 0.236926885056189_SPI_RK		   			 
		W(2) = 0.478628670499366_SPI_RK
		W(3) = 0.568888888888889_SPI_RK
		W(4) = 0.478628670499366_SPI_RK
		W(5) = 0.236926885056189_SPI_RK	
		!---------		   								
	ELSEif ( k == 5 ) THEN
		!---------	 	
		Node(1) = -0.932469514203152_SPI_RK			   	   
		Node(2) = -0.661209386466265_SPI_RK
		Node(3) = -0.238619186083197_SPI_RK		 
		Node(4) = 0.238619186083197_SPI_RK
		Node(5) = 0.661209386466265_SPI_RK
		Node(6) = 0.932469514203152_SPI_RK	
		!---------	
		W(1) = 0.171324492379170_SPI_RK	 	   			 
		W(2) = 0.360761573048139_SPI_RK
		W(3) = 0.467913934572691_SPI_RK
		W(4) = 0.467913934572691_SPI_RK
		W(5) = 0.360761573048139_SPI_RK
		W(6) = 0.171324492379170_SPI_RK
		!---------	
	ELSEif ( k == 6 ) THEN
		!---------	 	   	   
		Node(1) = -0.949107912342759_SPI_RK	
		Node(2) = -0.741531185500394_SPI_RK
		Node(3) = -0.404845151377397_SPI_RK		 
		Node(4) = 0.0_SPI_RK
		Node(5) = 0.404845151377397_SPI_RK
		Node(6) = 0.741531185500394_SPI_RK
		Node(7) = 0.949107912342759_SPI_RK	
		!---------	
		W(1) = 0.129484966168870_SPI_RK	 	   			 
		W(2) = 0.279705391489277_SPI_RK
		W(3) = 0.381830050505119_SPI_RK
		W(4) = 0.417959183673469_SPI_RK
		W(5) = 0.381830050505119_SPI_RK
		W(6) = 0.279705391489277_SPI_RK
		W(7) = 0.129484966168870_SPI_RK	 	   			 
		!---------		 	   											
	ELSEif ( k == 7 ) THEN
		!---------	 	
		Node(1) = -0.960289856497536_SPI_RK			   	   
		Node(2) = -0.796666477413627_SPI_RK
		Node(3) = -0.525532409916329_SPI_RK		 
		Node(4) = -0.183434642495650_SPI_RK
		Node(5) = 0.183434642495650_SPI_RK
		Node(6) = 0.525532409916329_SPI_RK
		Node(7) = 0.796666477413627_SPI_RK
		Node(8) = 0.960289856497536_SPI_RK	
		!---------	 
		W(1) = 0.101228536290376_SPI_RK		   			 
		W(2) = 0.222381034453374_SPI_RK
		W(3) = 0.313706645877887_SPI_RK
		W(4) = 0.362683783378362_SPI_RK
		W(5) = 0.362683783378362_SPI_RK
		W(6) = 0.313706645877887_SPI_RK
		W(7) = 0.222381034453374_SPI_RK
		W(8) = 0.101228536290376_SPI_RK	
		!---------	
	ELSEif ( k == 8 ) THEN
		!---------	
		Node(1) = -0.9681602395076261_SPI_RK			 	   	   
		Node(2) = -0.8360311073266358_SPI_RK
		Node(3) = -0.6133714327005904_SPI_RK		 
		Node(4) = -0.3242534234038089_SPI_RK
		Node(5) = 0.0_SPI_RK
		Node(6) = 0.3242534234038089_SPI_RK
		Node(7) = 0.6133714327005904_SPI_RK		
		Node(8) = 0.8360311073266358_SPI_RK
		Node(9) = 0.9681602395076261_SPI_RK			
		!---------	 
		W(1) = 0.08127438836157443_SPI_RK
		W(2) = 0.1806481606948574_SPI_RK
		W(3) = 0.2606106964029354_SPI_RK
		W(4) = 0.3123470770400029_SPI_RK
		W(5) = 0.3302393550012598_SPI_RK
		W(6) = 0.3123470770400029_SPI_RK
		W(7) = 0.2606106964029354_SPI_RK		
		W(8) = 0.1806481606948574_SPI_RK
		W(9) = 0.08127438836157443_SPI_RK
		!---------			
	ELSEif ( k == 9 ) THEN
		!---------	
		Node(1) = -0.9739065285171717_SPI_RK			 	   	   
		Node(2) = -0.8650633666889845_SPI_RK
		Node(3) = -0.6794095682990244_SPI_RK		 
		Node(4) = -0.4333953941292472_SPI_RK
		Node(5) = -0.1488743389816312_SPI_RK
		Node(6) = 0.1488743389816312_SPI_RK
		Node(7) = 0.4333953941292472_SPI_RK		
		Node(8) = 0.6794095682990244_SPI_RK
		Node(9) = 0.8650633666889845_SPI_RK			
		Node(10) = 0.9739065285171717_SPI_RK		
		!---------	 
		W(1) = 0.06667134430868804_SPI_RK
		W(2) = 0.1494513491505806_SPI_RK
		W(3) = 0.2190863625159821_SPI_RK
		W(4) = 0.2692667193099965_SPI_RK
		W(5) = 0.2955242247147530_SPI_RK
		W(6) = 0.2955242247147530_SPI_RK
		W(7) = 0.2692667193099965_SPI_RK		
		W(8) = 0.2190863625159821_SPI_RK
		W(9) = 0.1494513491505806_SPI_RK		
		W(10) = 0.06667134430868804_SPI_RK
		!---------
	ELSEif ( k == 14 ) THEN
		!---------	
		Node(1) = -0.9879925180204855_SPI_RK			 	   	   
		Node(2) = -0.9372733924007058_SPI_RK
		Node(3) = -0.8482065834104272_SPI_RK		 
		Node(4) = -0.7244177313601701_SPI_RK
		Node(5) = -0.5709721726085388_SPI_RK
		Node(6) = -0.3941513470775634_SPI_RK
		Node(7) = -0.2011940939974345_SPI_RK		
		Node(8) = 0.0_SPI_RK
		Node(9) = 0.2011940939974345_SPI_RK			
		Node(10) = 0.3941513470775634_SPI_RK		
		Node(11) = 0.5709721726085388_SPI_RK			 	   	   
		Node(12) = 0.7244177313601701_SPI_RK
		Node(13) = 0.8482065834104272_SPI_RK		 
		Node(14) = 0.9372733924007058_SPI_RK
		Node(15) = 0.9879925180204855_SPI_RK		
		!---------	 
		W(1) = 0.3075324199612_SPI_RK * 0.1_SPI_RK		   			 
		W(2) = 0.7036604748811_SPI_RK * 0.1_SPI_RK		   			 
		W(3) = 0.1071592204672_SPI_RK
		W(4) = 0.1395706779262_SPI_RK
		W(5) = 0.1662692058170_SPI_RK
		W(6) = 0.1861610000156_SPI_RK
		W(7) = 0.1984314853271_SPI_RK		
		W(8) = 0.2025782419256_SPI_RK
		W(9) = 0.1984314853271_SPI_RK		
		W(10) = 0.1861610000156_SPI_RK		
		W(11) = 0.1662692058170_SPI_RK		   			 
		W(12) = 0.1395706779262_SPI_RK
		W(13) = 0.1071592204672_SPI_RK
		W(14) = 0.7036604748811_SPI_RK * 0.1_SPI_RK		   			 
		W(15) = 0.3075324199612_SPI_RK * 0.1_SPI_RK		   			 
		!---------	
	ELSEif ( k == 15 ) THEN	
		!---------	 		
		Node(1) =  -0.9894009349916499_SPI_RK    
		Node(2) =  -0.9445750230732326_SPI_RK    
		Node(3) =  -0.8656312023878318_SPI_RK    
		Node(4) =  -0.7554044083550031_SPI_RK    
		Node(5) =  -0.6178762444026438_SPI_RK    
		Node(6) =  -0.4580167776572274_SPI_RK    
		Node(7) =  -0.2816035507792589_SPI_RK    
		Node(8) =  -0.09501250983763745_SPI_RK
		Node(9) =   0.09501250983763745_SPI_RK
		Node(10) =   0.2816035507792589_SPI_RK    
		Node(11) =   0.4580167776572274_SPI_RK    
		Node(12) =   0.6178762444026438_SPI_RK    
		Node(13) =   0.7554044083550031_SPI_RK    
		Node(14) =   0.8656312023878318_SPI_RK    
		Node(15) =   0.9445750230732326_SPI_RK    
		Node(16) =   0.9894009349916499_SPI_RK    
		!---------	 
		W(1) =   0.02715245941175400_SPI_RK
		W(2) =   0.06225352393864783_SPI_RK
		W(3) =   0.09515851168249274_SPI_RK
		W(4) =   0.1246289712555339_SPI_RK    
		W(5) =   0.1495959888165767_SPI_RK    
		W(6) =   0.1691565193950025_SPI_RK    
		W(7) =   0.1826034150449236_SPI_RK    
		W(8) =   0.1894506104550684_SPI_RK    
		W(9) =   0.1894506104550684_SPI_RK    
		W(10) =   0.1826034150449236_SPI_RK    
		W(11) =   0.1691565193950025_SPI_RK    
		W(12) =   0.1495959888165767_SPI_RK    
		W(13) =   0.1246289712555339_SPI_RK    
		W(14) =   0.09515851168249274_SPI_RK
		W(15) =   0.06225352393864783_SPI_RK
		W(16) =   0.02715245941175400_SPI_RK
						
	ELSEif ( k == 19 ) THEN
		!---------	
		Node(1) = 0.1761400713915198_SPI_RK * 0.1_SPI_RK		 	   	   
		Node(2) = 0.4060142980038685_SPI_RK * 0.1_SPI_RK
		Node(3) = 0.6267204833410915_SPI_RK * 0.1_SPI_RK
		Node(4) = 0.8327674157670482_SPI_RK * 0.1_SPI_RK
		Node(5) = 0.1019301198172404_SPI_RK
		Node(6) = 0.1181945319615183_SPI_RK
		Node(7) = 0.1316886384491766_SPI_RK		
		Node(8) = 0.1420961093183820_SPI_RK
		Node(9) = 0.1491729864726038_SPI_RK			
		Node(10) = 0.1527533871307258_SPI_RK	
		Node(11) = 0.1527533871307258_SPI_RK			 	   	   
		Node(12) = 0.1491729864726038_SPI_RK
		Node(13) = 0.1420961093183820_SPI_RK		 
		Node(14) = 0.1316886384491766_SPI_RK
		Node(15) = 0.1181945319615183_SPI_RK
		Node(16) = 0.1019301198172404_SPI_RK
		Node(17) = 0.8327674157670482_SPI_RK		
		Node(18) = 0.6267204833410915_SPI_RK
		Node(19) = 0.4060142980038685_SPI_RK			
		Node(20) = 0.1761400713915198_SPI_RK				
		!---------	 
		W(1) = -0.9931285991850950_SPI_RK		   			 
		W(2) = -0.9639719272779138_SPI_RK
		W(3) = -0.9122344282513258_SPI_RK
		W(4) = -0.8391169718222188_SPI_RK
		W(5) = -0.7463319064601508_SPI_RK
		W(6) = -0.6360536807265150_SPI_RK
		W(7) = -0.5108670019508271_SPI_RK		
		W(8) = -0.3737060887154195_SPI_RK
		W(9) = -0.2277858511416451_SPI_RK		
		W(10) = -0.7652652113349732_SPI_RK * 0.1_SPI_RK		
		W(11) = 0.7652652113349732_SPI_RK * 0.1_SPI_RK		   			 
		W(12) = 0.2277858511416451_SPI_RK
		W(13) = 0.3737060887154195_SPI_RK
		W(14) = 0.5108670019508271_SPI_RK
		W(15) = 0.6360536807265150_SPI_RK
		W(16) = 0.7463319064601508_SPI_RK
		W(17) = 0.8391169718222188_SPI_RK		
		W(18) = 0.9122344282513258_SPI_RK
		W(19) = 0.9639719272779138_SPI_RK		
		W(20) = 0.9931285991850950_SPI_RK		
		!---------	

	ELSEif ( k == 31 ) THEN
	
		!---------	
		Node(1) =    -0.9972638618494816    
		Node(2) =    -0.9856115115452684    
		Node(3) =    -0.9647622555875064    
		Node(4) =    -0.9349060759377397    
		Node(5) =   -0.8963211557660521_SPI_RK    
		Node(6) =   -0.8493676137325700_SPI_RK    
		Node(7) =   -0.7944837959679424_SPI_RK    
		Node(8) =   -0.7321821187402897_SPI_RK    
		Node(9) =   -0.6630442669302152_SPI_RK    
		Node(10) =   -0.5877157572407623_SPI_RK    
		Node(11) =   -0.5068999089322294_SPI_RK    
		Node(12) =   -0.4213512761306353_SPI_RK    
		Node(13) =   -0.3318686022821277_SPI_RK    
		Node(14) =   -0.2392873622521371_SPI_RK    
		Node(15) =   -0.1444719615827965_SPI_RK    
		Node(16) =   -0.04830766568773832_SPI_RK
		Node(17) =    0.04830766568773832_SPI_RK
		Node(18) =    0.1444719615827965_SPI_RK    
		Node(19) =    0.2392873622521371_SPI_RK    
		Node(20) =    0.3318686022821277_SPI_RK    
		Node(21) =    0.4213512761306353_SPI_RK    
		Node(22) =    0.5068999089322294_SPI_RK    
		Node(23) =    0.5877157572407623_SPI_RK    
		Node(24) =    0.6630442669302152_SPI_RK    
		Node(25) =    0.7321821187402897_SPI_RK    
		Node(26) =    0.7944837959679424_SPI_RK    
		Node(27) =    0.8493676137325700_SPI_RK    
		Node(28) =    0.8963211557660521_SPI_RK    
		Node(29) =    0.9349060759377397_SPI_RK    
		Node(30) =    0.9647622555875064_SPI_RK    
		Node(31) =    0.9856115115452684_SPI_RK    
		Node(32) =    0.9972638618494816_SPI_RK    
		!---------	
		W(1) =  0.007018610009470134_SPI_RK
		W(2) =  0.01627439473090564_SPI_RK
		W(3) =  0.02539206530926207_SPI_RK
		W(4) =  0.03427386291302147_SPI_RK
		W(5) =  0.04283589802222671_SPI_RK
		W(6) =  0.05099805926237618_SPI_RK
		W(7) =  0.05868409347853557_SPI_RK
		W(8) =  0.06582222277636185_SPI_RK
		W(9) =  0.07234579410884852_SPI_RK
		W(10) =  0.07819389578707037_SPI_RK
		W(11) =  0.08331192422694682_SPI_RK
		W(12) =  0.08765209300440385_SPI_RK
		W(13) =  0.09117387869576386_SPI_RK
		W(14) =  0.09384439908080454_SPI_RK
		W(15) =  0.09563872007927485_SPI_RK
		W(16) =  0.09654008851472780_SPI_RK
		W(17) =  0.09654008851472780_SPI_RK
		W(18) =  0.09563872007927485_SPI_RK
		W(19) =  0.09384439908080454_SPI_RK
		W(20) =  0.09117387869576386_SPI_RK
		W(21) =  0.08765209300440385_SPI_RK
		W(22) =  0.08331192422694682_SPI_RK
		W(23) =  0.07819389578707037_SPI_RK
		W(24) =  0.07234579410884852_SPI_RK
		W(25) =  0.06582222277636185_SPI_RK
		W(26) =  0.05868409347853557_SPI_RK
		W(27) =  0.05099805926237618_SPI_RK
		W(28) =  0.04283589802222671_SPI_RK
		W(29) =  0.03427386291302147_SPI_RK
		W(30) =  0.02539206530926207_SPI_RK
		W(31) =  0.01627439473090564_SPI_RK
		W(32) =  0.007018610009470134_SPI_RK
		!---------	
	ELSEif ( k == 63 ) THEN	
		!---------			
		Node(1) =    -0.9993050417357722_SPI_RK    
		Node(2) =    -0.9963401167719553_SPI_RK    
		Node(3) =    -0.9910133714767443_SPI_RK    
		Node(4) =    -0.9833362538846260_SPI_RK    
		Node(5) =    -0.9733268277899110_SPI_RK    
		Node(6) =    -0.9610087996520538_SPI_RK    
		Node(7) =    -0.9464113748584028_SPI_RK   
		Node(8) =    -0.9295691721319396_SPI_RK    
		Node(9) =    -0.9105221370785028_SPI_RK    
		Node(10) =    -0.8893154459951141_SPI_RK    
		Node(11) =    -0.8659993981540928_SPI_RK    
		Node(12) =    -0.8406292962525803_SPI_RK    
		Node(13) =    -0.8132653151227975_SPI_RK    
		Node(14) =    -0.7839723589433414_SPI_RK    
		Node(15) =    -0.7528199072605319_SPI_RK    
		Node(16) =    -0.7198818501716108_SPI_RK    
		Node(17) =    -0.6852363130542333_SPI_RK    
		Node(18) =    -0.6489654712546573_SPI_RK    
		Node(19) =    -0.6111553551723933_SPI_RK    
		Node(20) =    -0.5718956462026340_SPI_RK    
		Node(21) =    -0.5312794640198946_SPI_RK    
		Node(22) =    -0.4894031457070530_SPI_RK    
		Node(23) =    -0.4463660172534641_SPI_RK    
		Node(24) =    -0.4022701579639916_SPI_RK    
		Node(25) =    -0.3572201583376681_SPI_RK    
		Node(26) =    -0.3113228719902110_SPI_RK    
		Node(27) =    -0.2646871622087674_SPI_RK    
		Node(28) =    -0.2174236437400071_SPI_RK    
		Node(29) =    -0.1696444204239928_SPI_RK    
		Node(30) =    -0.1214628192961206_SPI_RK    
		Node(31) =    -0.07299312178779904_SPI_RK
		Node(32) =    -0.02435029266342443_SPI_RK			
		!---------	
		W(1) =        0.001783280721696293_SPI_RK
		W(2) =        0.004147033260562447_SPI_RK
		W(3) =        0.006504457968978380_SPI_RK
		W(4) =        0.008846759826363980_SPI_RK
		W(5) =        0.01116813946013109_SPI_RK
		W(6) =        0.01346304789671863_SPI_RK
		W(7) =        0.01572603047602472_SPI_RK
		W(8) =        0.01795171577569737_SPI_RK
		W(9) =        0.02013482315353015_SPI_RK
		W(10) =        0.02227017380838331_SPI_RK
		W(11) =        0.02435270256871090_SPI_RK
		W(12) =        0.02637746971505460_SPI_RK
		W(13) =        0.02833967261425945_SPI_RK
		W(14) =        0.03023465707240246_SPI_RK
		W(15) =        0.03205792835485158_SPI_RK
		W(16) =        0.03380516183714163_SPI_RK
		W(17) =        0.03547221325688238_SPI_RK
		W(18) =        0.03705512854023998_SPI_RK
		W(19) =        0.03855015317861563_SPI_RK
		W(20) =        0.03995374113272038_SPI_RK
		W(21) =        0.04126256324262342_SPI_RK
		W(22) =        0.04247351512365358_SPI_RK
		W(23) =        0.04358372452932344_SPI_RK
		W(24) =        0.04459055816375659_SPI_RK
		W(25) =        0.04549162792741814_SPI_RK
		W(26) =        0.04628479658131443_SPI_RK
		W(27) =        0.04696818281621006_SPI_RK
		W(28) =        0.04754016571483031_SPI_RK
		W(29) =        0.04799938859645831_SPI_RK
		W(30) =        0.04834476223480295_SPI_RK
		W(31) =        0.04857546744150343_SPI_RK
		W(32) =        0.04869095700913969_SPI_RK	
		!---------	
		do i = 1, (k+1)/2
			Node( i + (k+1)/2 ) = - Node( i )
			W( i + (k+1)/2 ) = W( i )			
		END do
		
!	ELSEif ( k == 9 ) THEN
!		!---------	
!		Node(1) = _SPI_RK			 	   	   
!		Node(2) = _SPI_RK
!		Node(3) = _SPI_RK		 
!		Node(4) = _SPI_RK
!		Node(5) = _SPI_RK
!		Node(6) = _SPI_RK
!		Node(7) = _SPI_RK		
!		Node(8) = _SPI_RK
!		Node(9) = _SPI_RK			
!		Node(10) = _SPI_RK		
!		!---------	 
!		W(1) = _SPI_RK		   			 
!		W(2) = _SPI_RK
!		W(3) = _SPI_RK
!		W(4) = _SPI_RK
!		W(5) = _SPI_RK
!		W(6) = _SPI_RK
!		W(7) = _SPI_RK		
!		W(8) = _SPI_RK
!		W(9) = _SPI_RK		
!		W(10) = SPI_RK		
		!---------					
	ELSE 
		print*,'computeLegendre : Error Not implemented Yet for n > 8'
	END if

END SUBROUTINE
  ! .............................................
   

END MODULE SPI_QUADRATURES_LINE 
