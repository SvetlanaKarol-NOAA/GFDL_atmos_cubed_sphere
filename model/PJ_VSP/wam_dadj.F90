
      SUBROUTINE WAM_DADJ(P2, DP2, T2, IM,NLEV)     
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NLEV, IM   
      REAL,PARAMETER :: cp= 1003.       
      REAL,PARAMETER :: GAMMA=9.5,RDG=.287./9.5,RDCP=RDG*GAMMA/cp
      REAL,PARAMETER :: P0=100.E2   !100 MB
      INTEGER        :: LOFF,NLAY
        
      REAL,INTENT(INOUT)   :: T2(IM,NLEV) , DP2(IM, NLEV)
      REAL,DIMENSION(NLEV+1),INTENT(IN):: P2(IM, NLEV+1)
      
      REAL,DIMENSION(NLEV)   :: PM,DP, R, Q , P,T
! INTERFACE PRESSURE LEVELS
      INTEGER:: I,J,K,L,N, ik
      
      INTEGER,DIMENSION(NLEV):: NML
      REAL,   DIMENSION(NLEV):: TETA,TPP,PDP      
 
  DO ik=1,IM
  
        P(1:NLEV+1) = P2(ik,1:NLEV+1)
        T(1:NLEV)   = T2(ik,1:NLEV)
        dp(1:NLEV)  = DP2(ik,1:NLEV)
      DO L=NLEV, 2, -1
         IF(P(L) <= P0) THEN
            LOFF= L
            EXIT
         ENDIF
      ENDDO  
      
        NLAY=LOFF

      DO L=1,NLEV
         PM(L)=.5*(P(L)+P(L+1))
      ENDDO

      R(NLEV)=1.
      DO L=NLEV-1,2,1
         R(L)=R(L+1)*(PM(L+1)/PM(L))**RDCP
	 Q(L)=DP(L)/R(L)
      ENDDO

! INITIALIZE FIRST COMBINED LAYER WITH FIRST MODEL LAYER

         I=NLAY
         K=I
         NML(K)=1
         TETA(K)=T(NLAY)*R(NLAY)
         PDP(K)=Q(NLAY)
         TPP(K)=TETA(K)*Q(NLAY)

! SCAN MODEL LAYERS (E.G., GOING UP FROM THE FIRST LAYER ABOVE OFFSET)

         DO L= NLAY, 1, -1

! INITIALIZE NEXT LAYER WITH CURRENT MODEL LAYER

            NML(K-1)=1
            TETA(K-1)=T(L)*R(L)
            PDP(K-1)=Q(L)
            TPP(K-1)=TETA(K-1)*Q(L)

! RECURSIVELY CHECK STABILITY WITH IMMEDIATELY UNDERLYING (COMBINED)
!     LAYER, UNTIL A STABLE STRATIFICATION IS FOUND OR THE BOTTOM LAYER

            DO J=1, K

! FOR MODEL LAYERS GOING DOWN THIS INEQUALITY SHOULD BE REVERSED

               IF(TETA(J) .gt. TETA(J+1)) THEN
                  
! STABLE STRATIFICATION - DO NOT COMBINE LAYERS, ADVANCE INDEX OF
!     COMBINED LAYERS (THE NUMBER OF COMBINED LAYERS CREATED TO THIS
!     POINT), GO TO NEXT MODEL LAYER

                  I=J+1
                  EXIT
               ELSE

! UNSTABLE - COBINE TWO LAYERS [J+1, J] INTO ONE
! LAYER J, REMEMBER ITS INDEX (THE NUMBER OF COMBINED LAYERS
!  

                  PDP(J)=PDP(J+1)+PDP(J)
                  TPP(J)=TPP(J+1)+TPP(J)
                  NML(J)=NML(J+1)+NML(J)
                  TETA(J)=TPP(J)/PDP(J)
                  I=J
               ENDIF
            ENDDO
            K=I
         ENDDO

! RETRIEVE TEMPERATURE FROM POTENTIAL TEMPERATURE OF (COMBINED) LAYERS,
!     SET STARTING MODEL LAYER INDEX

         L=K
         I=L
	 
         DO J=K, 1, -1
! SCAN ALL MODEL LAYER WITHIN EACH NEUTRAL LAYER, RESET STARTING INDEX

            DO L=I,I+NML(J)-1
               T(L)=TETA(J)/R(L)
            ENDDO
            I=L        
         ENDDO
         T2(N, 1:NLEV) =T(1:NLEV)
      ENDDO     !IM
      END SUBROUTINE WAM_DADJ

!***********************************************************************
!***********************************************************************
