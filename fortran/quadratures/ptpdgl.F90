!# -*- coding: utf8 -*-
!
!  ptpdgl.f
!  Maxwell_DG
!
!  Created by Ahmed RATNANI on 24/04/08.
!  Copyright 2008 __MyCompanyName__. All rights reserved.
!

!
      SUBROUTINE PTPDGL  ( M, Y, W )
!     *****************
!
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
	  integer, parameter  :: k_max = 200
      DIMENSION  Y(M) , X(k_max,k_max) , W(M)
      DIMENSION P1(0:k_max,k_max),P2(0:k_max,k_max),PP1(0:k_max,k_max),PP2(0:k_max,k_max),YM(k_max),YP(k_max),C(2)

!
!     CALCUL DES POINTS POUR L'INTERPOLATION DE GAUSS-LOBATTO-LEGENDRE
!     ----------------------------------------------------------------
!
      X(1,3) = - 1.D0
      X(2,3) =   0.D0
      X(3,3) =   1.D0
!
!----NOMBRE DES ITERATIONS ... ;
      IKM = 60
!
      DO 10 N = 2 , M-2
!
         DO 3 K = 1 , N
            YM(K) = X(K,N+1)
            YP(K) = X(K+1,N+1)
   3     CONTINUE
!
         DO 1 J  = 1 , N/2
         DO 1 IK = 1 , IKM

            C(1) = YM(J)
            C(2) = ( YM(J) + YP(J) ) * .5D0
               P1(0,J)  = 1.D0
               P1(1,J)  = C(1)
               PP1(1,J) = 1.D0
               PP1(2,J) = 3.D0 * C(1)

               P2(0,J)  = 1.D0
               P2(1,J)  = C(2)
               PP2(1,J) = 1.D0
               PP2(2,J) = 3.D0 * C(2)
!
               DO 4 I = 2 , N

                  XI  = I
                  X2I = 2*I
                  P1(I,J)    = ( ( X2I - 1.D0 ) * C(1) * P1(I-1,J) - ( XI - 1.D0 ) * P1(I-2,J) ) / XI
                  PP1(I+1,J) = ( ( X2I + 1.D0 ) * ( P1(I,J) + C(1) * PP1(I,J) ) - XI * PP1(I-1,J) ) / (XI + 1.D0)

                  P2(I,J)    = ( ( X2I - 1.D0 ) * C(2) * P2(I-1,J) - ( XI - 1.D0 ) * P2(I-2,J) ) / XI
                  PP2(I+1,J) = ( ( X2I + 1.D0 ) * ( P2(I,J) + C(2) * PP2(I,J) ) - XI * PP2(I-1,J) ) / (XI + 1.D0)

  4             CONTINUE
!
            SMO = PP1(N+1,J) * PP2(N+1,J)
!
            IF ( SMO .LT. 0.D0 )  THEN
                    YP(J) = ( YM(J) + YP(J) ) * .5D0
                ELSE
                    YM(J) = ( YM(J) + YP(J) ) * .5D0
            ENDIF

   1     CONTINUE
!
         DO 5 K = 1 , N
            X(K+1,N+2) = ( YM(K) + YP(K) ) * .5D0
   5     CONTINUE
!
         X(1,N+2) = - 1.D0
!
         IF ( MOD(N,2) .EQ. 1 ) X(N/2+2,N+2) = 0.D0
!
         DO 6 K = N/2+2 , N+2
            X(K,N+2) = - X(N+3-K,N+2)
   6     CONTINUE
!
!
 10   CONTINUE
!
      DO 7 I = 1 , M
           Y(I) = X(I,M)
  7   CONTINUE
!
!
!     CALCUL DES POIDS POUR L'INTERPOLATION DE GAUSS-LOBATTO-LEGENDRE
!     ---------------------------------------------------------------
      XM = M
!
      DO 8 J = 1 , M
!
         P1(0,J) = 1.D0
         P1(1,J) = X(J,M)
!
         DO 9 I = 2 , M-1
            XI = I
            P1(I,J) = ( ( XI + XI - 1.D0 ) * X(J,M) * P1(I-1,J) - ( XI - 1.D0 ) * P1(I-2,J) ) / XI
  9      CONTINUE
!
         IF(( J .NE. 1 ).AND.( J .NE. M ))  THEN
            W(J) = 2.D0/ ( (XM - 1.D0)* XM* P1(M-1,J)*P1(M-1,J))
         ENDIF
!
  8   CONTINUE
!
      W(1) = 2.D0 / ( ( XM - 1.D0 ) * XM )
      W(M) = W(1)
!
 14   FORMAT(1X,5F20.16)
!
      RETURN
      END
	  

 
	
