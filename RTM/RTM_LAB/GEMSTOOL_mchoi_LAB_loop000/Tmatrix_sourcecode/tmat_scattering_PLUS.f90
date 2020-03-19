 module tmat_scattering_plus

   use tmat_parameters, only : fpk => tmat_fpkind, &
                               cpk => tmat_cpkind, &
                               NPN1, NPN4, NPN6, NPL, NPL1
   use tmat_scattering, only : CCG, CCGIN, FACT, SIGNUM, DIRECT

!  LIST OF SUBROUTINES
!  ===================

!    GSP_PLUS
!    MATR_PLUS

!  Everything PUBLIC here
!  ----------------------

public

contains


SUBROUTINE GSP_plus                               &
     ( NMAX, NLIN, LAM, CSCA, L_CSCA,             & ! Inputs
       TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22,   & ! Inputs
       L_TR11,L_TR12,L_TR21,L_TR22,               & ! Inputs
       L_TI11,L_TI12,L_TI21,L_TI22,               & ! Inputs
       ALF1,ALF2,ALF3,ALF4,BET1,BET2,LMAX,        & ! Outputs
       L_ALF1,L_ALF2,L_ALF3,L_ALF4,L_BET1,L_BET2, & ! Outputs
       fail, message, trace )                       ! Outputs

!********************************************************************
!                                                                   *
!   CALCULATION OF THE EXPANSION COEFFICIENTS FOR (I,Q,U,V) -       *
!   REPRESENTATION.                                                 *
!                                                                   *
!   INPUT PARAMETERS:                                               *
!                                                                   *
!      LAM - WAVELENGTH OF LIGHT                                    *
!      CSCA - SCATTERING CROSS SECTION                              *
!      TR AND TI - ELEMENTS OF THE T-MATRIX. TRANSFERRED THROUGH    *
!                  COMMON /CTM/                                     *
!      NMAX - DIMENSION OF T(M)-MATRICES                            *
!                                                                   *
!   OUTPUT INFORTMATION:                                            *
!                                                                   *
!      ALF1,...,ALF4,BET1,BET2 - EXPANSION COEFFICIENTS             *
!      LMAX - NUMBER OF COEFFICIENTS MINUS 1                        *
!                                                                   *
!********************************************************************
 
!   use tmat_parameters, only : NPN1, NPN4, NPN6, NPL, NPL1

   implicit none

!  input arguments
!  ---------------

!  Numbers

   integer  , intent(in)  :: nmax, nlin
   real(fpk), intent(in)  :: CSCA, LAM
   real(fpk), intent(in)  :: L_CSCA(3)

!  Tmatrix input ( Kind = 4 )

!   real(kind=4), intent(in) :: TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4)
!   real(kind=4), intent(in) :: TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4)
!   real(kind=4), intent(in) :: TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4)
!   real(kind=4), intent(in) :: TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)

!   real(kind=4), intent(in) :: L_TR11(3,NPN6,NPN4,NPN4),L_TR12(3,NPN6,NPN4,NPN4)
!   real(kind=4), intent(in) :: L_TR21(3,NPN6,NPN4,NPN4),L_TR22(3,NPN6,NPN4,NPN4)
!   real(kind=4), intent(in) :: L_TI11(3,NPN6,NPN4,NPN4),L_TI12(3,NPN6,NPN4,NPN4)
!   real(kind=4), intent(in) :: L_TI21(3,NPN6,NPN4,NPN4),L_TI22(3,NPN6,NPN4,NPN4)

   real(fpk), intent(in) :: TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4)
   real(fpk), intent(in) :: TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4)
   real(fpk), intent(in) :: TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4)
   real(fpk), intent(in) :: TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)

   real(fpk), intent(in) :: L_TR11(3,NPN6,NPN4,NPN4),L_TR12(3,NPN6,NPN4,NPN4)
   real(fpk), intent(in) :: L_TR21(3,NPN6,NPN4,NPN4),L_TR22(3,NPN6,NPN4,NPN4)
   real(fpk), intent(in) :: L_TI11(3,NPN6,NPN4,NPN4),L_TI12(3,NPN6,NPN4,NPN4)
   real(fpk), intent(in) :: L_TI21(3,NPN6,NPN4,NPN4),L_TI22(3,NPN6,NPN4,NPN4)

!  Output arguments
!  ----------------

!  Number of coefficients

   integer   , intent(out) :: LMAX

!  Coefficients

   REAL(fpk), intent(out) :: ALF1(NPL),ALF2(NPL),ALF3(NPL)
   REAL(fpk), intent(out) :: ALF4(NPL),BET1(NPL),BET2(NPL)

!  Linearized Coefficients

   REAL(fpk), intent(out) :: L_ALF1(3,NPL),L_ALF2(3,NPL),L_ALF3(3,NPL)
   REAL(fpk), intent(out) :: L_ALF4(3,NPL),L_BET1(3,NPL),L_BET2(3,NPL)

!  Exception handling

   logical       , intent(out) :: fail
   character*(*) , intent(out) :: message
   character*(*) , intent(out) :: trace

!  Local arrays
!  ------------

!  FACT and SIGNUM output

   real(fpk) :: SSIGN(900),FCALC(900)

!  Miscelllaneous numerical factor arrays ( Kind = 8 )

   real(fpk) :: SSI(NPL),SSJ(NPN1)
   real(fpk) :: G1(NPL1,NPN6),G2(NPL1,NPN6)
   real(fpk) :: FR(NPN4,NPN4),FI(NPN4,NPN4),FF(NPN4,NPN4)

!  T arrays. ( Kind = 8 )

   real(fpk) :: TR1(NPL1,NPN4),TR2(NPL1,NPN4)
   real(fpk) :: TI1(NPL1,NPN4),TI2(NPL1,NPN4)
   real(fpk) :: L_TR1(3,NPL1,NPN4),L_TR2(3,NPL1,NPN4)
   real(fpk) :: L_TI1(3,NPL1,NPN4),L_TI2(3,NPL1,NPN4)

!  A arrays. ( Kind = 8 )

   real(fpk) :: AR1(NPN4),AR2(NPN4),AI1(NPN4),AI2(NPN4)
   real(fpk) :: L_AR1(3,NPN4),L_AR2(3,NPN4),L_AI1(3,NPN4),L_AI2(3,NPN4)

!  B arrays. ( Kind = 4 )

!   real(kind=4) :: B1R(NPL1,NPL1,NPN4),B1I(NPL1,NPL1,NPN4)
!   real(kind=4) :: B2R(NPL1,NPL1,NPN4),B2I(NPL1,NPL1,NPN4)
!   real(kind=4) :: L_B1R(3,NPL1,NPL1,NPN4),L_B1I(3,NPL1,NPL1,NPN4)
!   real(kind=4) :: L_B2R(3,NPL1,NPL1,NPN4),L_B2I(3,NPL1,NPL1,NPN4)

   real(fpk) :: B1R(NPL1,NPL1,NPN4),B1I(NPL1,NPL1,NPN4)
   real(fpk) :: B2R(NPL1,NPL1,NPN4),B2I(NPL1,NPL1,NPN4)
   real(fpk) :: L_B1R(3,NPL1,NPL1,NPN4),L_B1I(3,NPL1,NPL1,NPN4)
   real(fpk) :: L_B2R(3,NPL1,NPL1,NPN4),L_B2I(3,NPL1,NPL1,NPN4)

!  D arrays. ( Kind = 4 )

!   real(kind=4) :: D1(NPL1,NPN4,NPN4),D2(NPL1,NPN4,NPN4)
!   real(kind=4) :: D3(NPL1,NPN4,NPN4),D4(NPL1,NPN4,NPN4)
!   real(kind=4) :: D5R(NPL1,NPN4,NPN4),D5I(NPL1,NPN4,NPN4)
!   real(kind=4) :: L_D1(3,NPL1,NPN4,NPN4),L_D2(3,NPL1,NPN4,NPN4)
!   real(kind=4) :: L_D3(3,NPL1,NPN4,NPN4),L_D4(3,NPL1,NPN4,NPN4)
!   real(kind=4) :: L_D5R(3,NPL1,NPN4,NPN4),L_D5I(3,NPL1,NPN4,NPN4)

   real(fpk) :: D1(NPL1,NPN4,NPN4),D2(NPL1,NPN4,NPN4)
   real(fpk) :: D3(NPL1,NPN4,NPN4),D4(NPL1,NPN4,NPN4)
   real(fpk) :: D5R(NPL1,NPN4,NPN4),D5I(NPL1,NPN4,NPN4)
   real(fpk) :: L_D1(3,NPL1,NPN4,NPN4),L_D2(3,NPL1,NPN4,NPN4)
   real(fpk) :: L_D3(3,NPL1,NPN4,NPN4),L_D4(3,NPL1,NPN4,NPN4)
   real(fpk) :: L_D5R(3,NPL1,NPN4,NPN4),L_D5I(3,NPL1,NPN4,NPN4)

!  Complex array. ( Kind = 16 )

   COMPLEX(cpk) :: CIM(NPN1)
 
!  Old commons
!      COMMON /TMAT/ TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22
!      COMMON /CBESS/ B1R,B1I,B2R,B2I    
!      COMMON /SS/ SSIGN
!      EQUIVALENCE ( PLUS1(1),TR11(1,1,1) )
!      EQUIVALENCE (D1(1,1,1),PLUS1(1)),        
!     &            (D2(1,1,1),PLUS1(NPL1*NPN4*NPN4+1)),
!     &            (D3(1,1,1),PLUS1(NPL1*NPN4*NPN4*2+1)),
!     &            (D4(1,1,1),PLUS1(NPL1*NPN4*NPN4*3+1)), 
!     &            (D5R(1,1,1),PLUS1(NPL1*NPN4*NPN4*4+1)) 

!  Other variables
!  ===============

!  Complex Kind = 16

   COMPLEX(cpk) :: CI, CCI, CCJ

!  Double precision

   REAL(fpk) :: SIG, SSS, SXS, HELP, SSS1, SI, SJ, SL, T1, T2, T3, T4
   REAL(fpk) :: TT1, TT2, TT3, TT4, TT5, TT6, TT7, TT8
   REAL(fpk) :: AAR1, AAR2, AAI1, AAI2, RR1, RI1, RR2, RI2
   REAL(fpk) :: BBR1, BBR2, BBI1, BBI2, DK, FFN, XR, XI, XX
   REAL(fpk) :: L_DK(3), L_SL(3)

   REAL(fpk) :: X1, X2, X3, X4, X5, X6, X7, X8
   REAL(fpk) :: DD1, DD2, DD3, DD4, DD5R, DD5I
   REAL(fpk) :: DM1, DM2, DM3, DM4, DM5R, DM5I
   REAL(fpk) :: G1L, G2L, G3L, G4L, G5LR, G5LI

   REAL(fpk) :: L_X1, L_X2, L_X3, L_X4, L_X5, L_X6, L_X7, L_X8
   REAL(fpk) :: L_DD1(3), L_DD2(3), L_DD3(3), L_DD4(3), L_DD5R(3), L_DD5I(3)
   REAL(fpk) :: L_DM1(3), L_DM2(3), L_DM3(3), L_DM4(3), L_DM5R(3), L_DM5I(3)
   REAL(fpk) :: L_G1L(3), L_G2L(3), L_G3L(3), L_G4L(3), L_G5LR(3), L_G5LI(3)

!  Others

   INTEGER   :: N, NN, NNN, N1, NN1, NNMIN, NNMAX, NN1MIN, NN1MAX, NMAX1 
   INTEGER   :: L1MAX, M, M1, MMAX, MMIN, M1MAX, M1MIN, M2, NL, NIFORT
   INTEGER   :: I, I1, J, K1, K2, K3, K4, K5, K6, KN, L, L1, LL
   LOGICAL   :: Calculating

!  Initial calculations
!  ====================

!  Exception handling

      FAIL = .false.
      message = ' '
      TRACE = ' '

!  FACT and SIGNUM

      CALL FACT   (FCALC)
      CALL SIGNUM (SSIGN)

!  Initial setups

      LMAX  = 2*NMAX
      L1MAX = LMAX+1

      CI=(0D0,1D0)
      CIM(1)=CI
      DO 2 I=2,NMAX
         CIM(I)=CIM(I-1)*CI
    2 CONTINUE

      SSI(1)=1D0
      DO 3 I=1,LMAX
         I1=I+1
         SI=DBLE(2*I+1)
         SSI(I1)=SI
         IF(I.LE.NMAX) SSJ(I)=DSQRT(SI)
    3 CONTINUE

      CI=-CI
      DO 5 I=1,NMAX
         SI=SSJ(I)
         CCI=CIM(I)
         DO 4 J=1,NMAX
            SJ=1D0/SSJ(J)
            CCJ=CIM(J)*SJ/CCI
            FR(J,I)=CCJ
            FI(J,I)=CCJ*CI
            FF(J,I)=SI*SJ
    4    CONTINUE
    5 CONTINUE
      NMAX1=NMAX+1
 
!  CALCULATION OF THE ARRAYS B1 AND B2
!  ===================================
 
      K1=1
      K2=0
      K3=0
      K4=1
      K5=1
      K6=2

!  Start 100 loop
 
      DO 100 N=1,NMAX

!  Dummy argument (to satisfy IFORT compiler)

         NIFORT = N

!  CALCULATION OF THE ARRAYS T1 AND T2 
 
         DO 10 NN=1,NMAX
            M1MAX=MIN0(N,NN)+1
            DO 6 M1=1,M1MAX
               M=M1-1

!  Regular

               L1=NPN6+M
               TT1=TR11(M1,N,NN)  ;  TT2=TR12(M1,N,NN)
               TT3=TR21(M1,N,NN)  ;  TT4=TR22(M1,N,NN)
               TT5=TI11(M1,N,NN)  ;  TT6=TI12(M1,N,NN)
               TT7=TI21(M1,N,NN)  ;  TT8=TI22(M1,N,NN)
               T1=TT1+TT2  ;  T2=TT3+TT4
               T3=TT5+TT6  ;  T4=TT7+TT8
               TR1(L1,NN)=T1+T2  ;  TR2(L1,NN)=T1-T2
               TI1(L1,NN)=T3+T4  ;  TI2(L1,NN)=T3-T4
               IF (M.NE.0) THEN
                  L1=NPN6-M
                  T1=TT1-TT2  ;  T2=TT3-TT4
                  T3=TT5-TT6  ;  T4=TT7-TT8
                  TR1(L1,NN)=T1-T2  ;  TR2(L1,NN)=T1+T2
                  TI1(L1,NN)=T3-T4  ;  TI2(L1,NN)=T3+T4
               ENDIF

!  linearized. Restore L1 index (Very important bug, 3/29/11)

               DO LL = 1, NLIN
                  L1=NPN6+M
                  TT1=L_TR11(LL,M1,N,NN)  ;  TT2=L_TR12(LL,M1,N,NN)
                  TT3=L_TR21(LL,M1,N,NN)  ;  TT4=L_TR22(LL,M1,N,NN)
                  TT5=L_TI11(LL,M1,N,NN)  ;  TT6=L_TI12(LL,M1,N,NN)
                  TT7=L_TI21(LL,M1,N,NN)  ;  TT8=L_TI22(LL,M1,N,NN)
                  T1=TT1+TT2  ;  T2=TT3+TT4
                  T3=TT5+TT6  ;  T4=TT7+TT8
                  L_TR1(LL,L1,NN)=T1+T2  ;  L_TR2(LL,L1,NN)=T1-T2
                  L_TI1(LL,L1,NN)=T3+T4  ;  L_TI2(LL,L1,NN)=T3-T4
                  IF (M.NE.0) THEN
                     L1=NPN6-M
                     T1=TT1-TT2  ;  T2=TT3-TT4
                     T3=TT5-TT6  ;  T4=TT7-TT8
                     L_TR1(LL,L1,NN)=T1-T2  ;  L_TR2(LL,L1,NN)=T1+T2
                     L_TI1(LL,L1,NN)=T3-T4  ;  L_TI2(LL,L1,NN)=T3+T4
                  ENDIF
               ENDDO

!  End 6, 10 loops

    6       CONTINUE
   10    CONTINUE

!  debug 3/30/11
!         do nn1 = 1, 19
!            do L1 = 79,83
!               write(68,*)nn1,nn1,l1,n,TR1(L1,NN1)
!               write(69,*)nn1,nn1,l1,n,L_TR1(3,L1,NN1)
!            enddo
!         enddo

!  Main loop (40) for A and B

         NN1MAX=NMAX1+N
         DO 40 NN1=1,NN1MAX
            N1=NN1-1 

!  CALCULATION OF THE ARRAYS A1 AND A2
!  -----------------------------------

!  CLEBSCH-GORDAN call

            CALL CCG ( NIFORT, N1, NMAX, K1, K2, FCALC, SSIGN, & ! Input
                       G1, fail, message )                       ! Output 

            if ( fail ) then
               trace = 'First Call in GSP_plus to CCG, Loop 40 A calculation'
               return
            endif

!  Resume

            NNMAX=MIN0(NMAX,N1+N)
            NNMIN=MAX0(1,IABS(N-N1))
            KN=N+NN1
            DO 15 NN=NNMIN,NNMAX
               NNN=NN+1
               SIG=SSIGN(KN+NN)
               M1MAX=MIN0(N,NN)+NPN6

!  Regular

               AAR1=0D0  ;  AAR2=0D0  ;  AAI1=0D0  ;  AAI2=0D0 
               DO 13 M1=NPN6,M1MAX
                  M=M1-NPN6
                  SSS=G1(M1,NNN)
                  RR1=TR1(M1,NN)  ;  RI1=TI1(M1,NN)
                  RR2=TR2(M1,NN)  ;  RI2=TI2(M1,NN)
                  IF(M.NE.0) THEN
                     M2=NPN6-M
                     RR1=RR1+TR1(M2,NN)*SIG  ;  RI1=RI1+TI1(M2,NN)*SIG
                     RR2=RR2+TR2(M2,NN)*SIG  ;  RI2=RI2+TI2(M2,NN)*SIG
                  ENDIF
                  AAR1=AAR1+SSS*RR1  ;  AAI1=AAI1+SSS*RI1
                  AAR2=AAR2+SSS*RR2  ;  AAI2=AAI2+SSS*RI2
   13          CONTINUE
               XR=FR(NN,N)
               XI=FI(NN,N)
               AR1(NN)=AAR1*XR-AAI1*XI  ;  AI1(NN)=AAR1*XI+AAI1*XR
               AR2(NN)=AAR2*XR-AAI2*XI  ;  AI2(NN)=AAR2*XI+AAI2*XR

!               write(68,*)nn,nn1,m1,n,AI2(NN)

!  Linearized

               DO LL = 1, NLIN
                  AAR1=0D0  ;  AAR2=0D0  ;  AAI1=0D0  ;  AAI2=0D0
                  DO 133 M1=NPN6,M1MAX
                     M=M1-NPN6
                     SSS=G1(M1,NNN)
                     RR1=L_TR1(LL,M1,NN)  ;  RI1=L_TI1(LL,M1,NN)
                     RR2=L_TR2(LL,M1,NN)  ;  RI2=L_TI2(LL,M1,NN)
                     IF(M.NE.0) THEN
                        M2=NPN6-M
                        RR1=RR1+L_TR1(LL,M2,NN)*SIG  ;  RI1=RI1+L_TI1(LL,M2,NN)*SIG
                        RR2=RR2+L_TR2(LL,M2,NN)*SIG  ;  RI2=RI2+L_TI2(LL,M2,NN)*SIG
                     ENDIF
                     AAR1=AAR1+SSS*RR1  ;  AAI1=AAI1+SSS*RI1
                     AAR2=AAR2+SSS*RR2  ;  AAI2=AAI2+SSS*RI2
   133            CONTINUE
                  L_AR1(LL,NN)=AAR1*XR-AAI1*XI  ;  L_AI1(LL,NN)=AAR1*XI+AAI1*XR
                  L_AR2(LL,NN)=AAR2*XR-AAI2*XI  ;  L_AI2(LL,NN)=AAR2*XI+AAI2*XR

!            if ( ll.eq.3)write(69,*)nn,nn1,m1,n,L_AI2(LL,NN)

               ENDDO

   15       CONTINUE
 
!  Finished A1/A2, Now complete B1/B2 computation
!  ----------------------------------------------

!  CLEBSCH-GORDAN call

            CALL CCG ( NIFORT, N1, NMAX, K3, K4, FCALC, SSIGN, & ! Input
                       G2, fail, message )                       ! Output 
            if ( fail ) then
               trace = 'Second Call in GSP to CCG, Loop 40 B calculation'
              return
            endif

!  Resume

            M1=MAX0(-N1+1,-N)
            M2=MIN0(N1+1,N)
            M1MAX=M2+NPN6
            M1MIN=M1+NPN6

!  regular

            DO 30 M1=M1MIN,M1MAX
               BBR1=0D0   ;  BBI1=0D0   ;  BBR2=0D0   ;  BBI2=0D0
               DO 25 NN=NNMIN,NNMAX
                  NNN=NN+1
                  SSS=G2(M1,NNN)
                  BBR1=BBR1+SSS*AR1(NN)   ;  BBI1=BBI1+SSS*AI1(NN)
                  BBR2=BBR2+SSS*AR2(NN)   ;  BBI2=BBI2+SSS*AI2(NN)
   25          CONTINUE
               B1R(NN1,M1,N)=BBR1   ;  B1I(NN1,M1,N)=BBI1
               B2R(NN1,M1,N)=BBR2   ;  B2I(NN1,M1,N)=BBI2

!               write(68,*)m1,nn1,n,n,B1R(NN1,M1,N),&
!               B2R(NN1,M1,N),B1I(NN1,M1,N),B2I(NN1,M1,N)

   30       CONTINUE

!  Linearized

            DO LL = 1, NLIN
               DO 330 M1=M1MIN,M1MAX
                  BBR1=0D0   ;  BBI1=0D0   ;  BBR2=0D0   ;  BBI2=0D0
                  DO 235 NN=NNMIN,NNMAX
                     NNN=NN+1
                     SSS=G2(M1,NNN)
                     BBR1=BBR1+SSS*L_AR1(LL,NN)   ;  BBI1=BBI1+SSS*L_AI1(LL,NN)
                     BBR2=BBR2+SSS*L_AR2(LL,NN)   ;  BBI2=BBI2+SSS*L_AI2(LL,NN)
   235             CONTINUE
                  L_B1R(LL,NN1,M1,N)=BBR1   ;  L_B1I(LL,NN1,M1,N)=BBI1
                  L_B2R(LL,NN1,M1,N)=BBR2   ;  L_B2I(LL,NN1,M1,N)=BBI2

!                  if ( ll.eq.3)write(69,*)m1,nn1,n,n,L_B1R(3,NN1,M1,N),&
!                  L_B2R(3,NN1,M1,N),L_B1I(3,NN1,M1,N),L_B2I(3,NN1,M1,N)

   330         CONTINUE
            ENDDO

!  Finish 40 loop

   40    CONTINUE

!  Finish 100 loop

  100 CONTINUE
  
!      pause'B-arrays Linearized 100'

!  CALCULATION OF THE ARRAYS D1,D2,D3,D4, AND D5
!  =============================================

!      jcount = 0
 
      DO 200 N=1,NMAX
         DO 190 NN=1,NMAX
            M1=MIN0(N,NN)
            M1MAX=NPN6+M1
            M1MIN=NPN6-M1
            NN1MAX=NMAX1+MIN0(N,NN)

            DO 180 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD1   = 0D0   ;   DD2 = 0D0
               L_DD1 = 0.0d0 ; L_DD2 = 0.0d0
               DO 150 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)    ;  X2=B1I(NN1,M1,N)
                  X3=B1R(NN1,M1,NN)   ;  X4=B1I(NN1,M1,NN)
                  X5=B2R(NN1,M1,N)    ;  X6=B2I(NN1,M1,N)
                  X7=B2R(NN1,M1,NN)   ;  X8=B2I(NN1,M1,NN)
                  DD1=DD1+XX*(X1*X3+X2*X4)
                  DD2=DD2+XX*(X5*X7+X6*X8)
                  DO LL = 1, NLIN
                     L_X1=L_B1R(LL,NN1,M1,N)    ;  L_X2=L_B1I(LL,NN1,M1,N)
                     L_X3=L_B1R(LL,NN1,M1,NN)   ;  L_X4=L_B1I(LL,NN1,M1,NN)
                     L_X5=L_B2R(LL,NN1,M1,N)    ;  L_X6=L_B2I(LL,NN1,M1,N)
                     L_X7=L_B2R(LL,NN1,M1,NN)   ;  L_X8=L_B2I(LL,NN1,M1,NN)
                     L_DD1(LL)=L_DD1(LL)+XX*(L_X1*X3+L_X2*X4+X1*L_X3+X2*L_X4)
                     L_DD2(LL)=L_DD2(LL)+XX*(L_X5*X7+L_X6*X8+X5*L_X7+X6*L_X8)
                  ENDDO
  150          CONTINUE
               D1(M1,NN,N)=DD1  ;  D2(M1,NN,N)=DD2
               DO LL = 1, NLIN
                  L_D1(LL,M1,NN,N)=L_DD1(LL)  ;  L_D2(LL,M1,NN,N)=L_DD2(LL)
               ENDDO
  180       CONTINUE

            MMAX=MIN0(N,NN+2)
            MMIN=MAX0(-N,-NN+2)
            M1MAX=NPN6+MMAX
            M1MIN=NPN6+MMIN
            DO 186 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD3=0D0   ;  DD4=0D0     ;  DD5R=0D0     ;  DD5I=0D0
               L_DD3=0D0 ;  L_DD4=0D0   ;  L_DD5R=0D0   ;  L_DD5I=0D0
               M2=-M+2+NPN6
               DO 183 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)    ;  X2=B1I(NN1,M1,N)
                  X3=B2R(NN1,M1,N)    ;  X4=B2I(NN1,M1,N)
                  X5=B1R(NN1,M2,NN)   ;  X6=B1I(NN1,M2,NN)
                  X7=B2R(NN1,M2,NN)   ;  X8=B2I(NN1,M2,NN)
                  DD3=DD3+XX*(X1*X5+X2*X6)
                  DD4=DD4+XX*(X3*X7+X4*X8)
                  DD5R=DD5R+XX*(X3*X5+X4*X6)
                  DD5I=DD5I+XX*(X4*X5-X3*X6)
                  DO LL = 1, NLIN
                     L_X1=L_B1R(LL,NN1,M1,N)    ;  L_X2=L_B1I(LL,NN1,M1,N)
                     L_X3=L_B2R(LL,NN1,M1,N)    ;  L_X4=L_B2I(LL,NN1,M1,N)
                     L_X5=L_B1R(LL,NN1,M2,NN)   ;  L_X6=L_B1I(LL,NN1,M2,NN)
                     L_X7=L_B2R(LL,NN1,M2,NN)   ;  L_X8=L_B2I(LL,NN1,M2,NN)
                     L_DD3(LL) =L_DD3(LL) +XX*(L_X1*X5+L_X2*X6+X1*L_X5+X2*L_X6)
                     L_DD4(LL) =L_DD4(LL) +XX*(L_X3*X7+L_X4*X8+X3*L_X7+X4*L_X8)
                     L_DD5R(LL)=L_DD5R(LL)+XX*(L_X3*X5+L_X4*X6+X3*L_X5+X4*L_X6)
                     L_DD5I(LL)=L_DD5I(LL)+XX*(L_X4*X5-L_X3*X6+X4*L_X5-X3*L_X6)
                  ENDDO
  183          CONTINUE
               D3(M1,NN,N)=DD3     ;  D4(M1,NN,N)=DD4
               D5R(M1,NN,N)=DD5R   ;  D5I(M1,NN,N)=DD5I
               DO LL = 1, NLIN
                  L_D3(LL,M1,NN,N) =L_DD3(LL)   ;  L_D4(LL,M1,NN,N) =L_DD4(LL)
                  L_D5R(LL,M1,NN,N)=L_DD5R(LL)  ;  L_D5I(LL,M1,NN,N)=L_DD5I(LL)
               ENDDO

!               jcount = jcount + 1
!               write(68,*)m1,nn,n,n,D1(M1,NN,N),D2(M1,NN,N),&
!               D3(M1,NN,N),D4(M1,NN,N),D5R(M1,NN,N),D5I(M1,NN,N)
!               write(69,*)m1,nn,n,n,L_D1(3,M1,NN,N),L_D2(3,M1,NN,N),&
!               L_D3(3,M1,NN,N),L_D4(3,M1,NN,N),L_D5R(3,M1,NN,N),L_D5I(3,M1,NN,N)

  186       CONTINUE
  190    CONTINUE
  200 CONTINUE

!      write(*,*)jcount
!      pause'D-arrays Linearized 200'

!  CALCULATION OF THE EXPANSION COEFFICIENTS
!  ========================================= 

!  Normalization
 
      DK = LAM*LAM/(4D0*CSCA*DACOS(-1D0))

!      HELP =  - 2.0d0 * DK / CSCA

      HELP =  - DK / CSCA
      DO LL = 1, NLIN
          L_DK(LL) = HELP * L_CSCA(LL)
      ENDDO

!  Initialize

      Calculating = .true.
      L1 = 0

!  Main loop (Formerly the 300 loop)

!      DO 300 L1=1,L1MAX      (Former Code)

      DO WHILE (Calculating .and. L1.LT.L1MAX)

         L1 = L1 + 1

!  Start

         G1L=0D0     ;  G2L=0D0      ;  G3L=0D0 
         G4L=0D0     ;  G5LR=0D0     ;  G5LI=0D0
         L_G1L=0D0   ;  L_G2L=0D0    ;  L_G3L=0D0 
         L_G4L=0D0   ;  L_G5LR=0D0   ;  L_G5LI=0D0

         L=L1-1
         SL=SSI(L1)*DK
         DO LL = 1, NLIN
            L_SL(LL) = SSI(L1) * L_DK(LL)
         ENDDO

!  290 loop

         DO 290 N=1,NMAX

            NNMIN=MAX0(1,IABS(N-L))
            NNMAX=MIN0(NMAX,N+L)

!  Dummy argument (to satisfy IFORT compiler)

            NIFORT = N

!  All done if this condition holds

            IF(NNMAX.LT.NNMIN) GO TO 290

!  CLEBSCH-GORDAN call # 3

            CALL CCG ( NIFORT, L, NMAX, K1, K2, FCALC, SSIGN, & ! Input
                       G1, fail, message )                      ! Output 
            if ( fail ) then
               trace = 'Third Call in GSP to CCG, Loop 290 calculation'
               return
            endif

!  CLEBSCH-GORDAN call # 4 ( Only for L > 1

            IF ( L.GE.2 ) THEN
               CALL CCG ( NIFORT, L, NMAX, K5, K6, FCALC, SSIGN, & ! Input
                          G2, fail, message )                      ! Output 
               if ( fail ) then
                  trace = 'Fourth Call in GSP to CCG, Loop 290 calculation'
                  return
               endif
            ENDIF

!  Resume

            NL=N+L

!  280 loop

            DO 280  NN=NNMIN,NNMAX

               NNN=NN+1
               MMAX=MIN0(N,NN)
               M1MIN=NPN6-MMAX
               M1MAX=NPN6+MMAX
               SI=SSIGN(NL+NNN)

               DM1  =0D0  ; DM2  =0D0
               L_DM1=0D0  ; L_DM2=0D0
               DO 270 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  IF(M.GE.0) SSS1=G1(M1,NNN)
                  IF(M.LT.0) SSS1=G1(NPN6-M,NNN)*SI
                  DM1=DM1+SSS1*D1(M1,NN,N)
                  DM2=DM2+SSS1*D2(M1,NN,N)
                  DO LL = 1, NLIN
                     L_DM1(LL)=L_DM1(LL)+SSS1*L_D1(LL,M1,NN,N)
                     L_DM2(LL)=L_DM2(LL)+SSS1*L_D2(LL,M1,NN,N)
                  ENDDO
  270          CONTINUE

               FFN=FF(NN,N)
               SSS=G1(NPN6+1,NNN)*FFN
               SXS=SSS*SI
               G1L=G1L+SSS*DM1
               G2L=G2L+SXS*DM2
               DO LL = 1, NLIN
                  L_G1L(LL)=L_G1L(LL)+SSS*L_DM1(LL)
                  L_G2L(LL)=L_G2L(LL)+SXS*L_DM2(LL)
               ENDDO

               IF(L.GE.2) THEN
                  DM3=0D0   ;  DM4=0D0   ;  DM5R=0D0   ;  DM5I=0D0
                  L_DM3=0D0 ;  L_DM4=0D0 ;  L_DM5R=0D0 ;  L_DM5I=0D0
                  MMAX=MIN0(N,NN+2)
                  MMIN=MAX0(-N,-NN+2)
                  M1MAX=NPN6+MMAX
                  M1MIN=NPN6+MMIN

                  DO 275 M1=M1MIN,M1MAX
                     M=M1-NPN6
                     SSS1=G2(NPN6-M,NNN)
                     DM3=DM3+SSS1*D3(M1,NN,N)
                     DM4=DM4+SSS1*D4(M1,NN,N)
                     DM5R=DM5R+SSS1*D5R(M1,NN,N)
                     DM5I=DM5I+SSS1*D5I(M1,NN,N)
                     DO LL = 1, NLIN
                        L_DM3(LL) =L_DM3(LL) +SSS1*L_D3(LL,M1,NN,N)
                        L_DM4(LL) =L_DM4(LL) +SSS1*L_D4(LL,M1,NN,N)
                        L_DM5R(LL)=L_DM5R(LL)+SSS1*L_D5R(LL,M1,NN,N)
                        L_DM5I(LL)=L_DM5I(LL)+SSS1*L_D5I(LL,M1,NN,N)
                     ENDDO
  275             CONTINUE

                  G5LR=G5LR-SSS*DM5R
                  G5LI=G5LI-SSS*DM5I
                  DO LL = 1, NLIN
                     L_G5LR(LL)=L_G5LR(LL)-SSS*L_DM5R(LL)
                     L_G5LI(LL)=L_G5LI(LL)-SSS*L_DM5I(LL)
                  ENDDO

                  SSS=G2(NPN4,NNN)*FFN
                  SXS=SSS*SI
                  G3L=G3L+SSS*DM3
                  G4L=G4L+SXS*DM4
                  DO LL = 1, NLIN
                     L_G3L(LL) =L_G3L(LL) +SSS*L_DM3(LL)
                     L_G4L(LL) =L_G4L(LL) +SXS*L_DM4(LL)
                  ENDDO

               ENDIF

  280       CONTINUE

  290    CONTINUE

!  Final

         DO LL = 1, NLIN
            L_G1L(LL)  = L_G1L(LL)  * SL + G1L  * L_SL(LL)
            L_G2L(LL)  = L_G2L(LL)  * SL + G2L  * L_SL(LL)
            L_G3L(LL)  = L_G3L(LL)  * SL + G3L  * L_SL(LL)
            L_G4L(LL)  = L_G4L(LL)  * SL + G4L  * L_SL(LL)
            L_G5LR(LL) = L_G5LR(LL) * SL + G5LR * L_SL(LL)
            L_G5LI(LL) = L_G5LI(LL) * SL + G5LI * L_SL(LL)
         ENDDO
         G1L=G1L*SL   ;  G2L=G2L*SL    ; G3L=G3L*SL
         G4L=G4L*SL   ; G5LR=G5LR*SL   ; G5LI=G5LI*SL

!        write(68,*)L,L,L1,L1,G1L,G2L,G3L,G4L,G5LR,G5LI
!        write(69,*)m1,nn,n,n,L_G1L(3),L_G2L(3),L_G3L(3),L_G4L(3),L_G5LR(3),L_G5LI(3)

         ALF1(L1)=G1L+G2L
         ALF2(L1)=G3L+G4L
         ALF3(L1)=G3L-G4L
         ALF4(L1)=G1L-G2L
         BET1(L1)=G5LR*2D0
         BET2(L1)=G5LI*2D0
         DO LL = 1, NLIN
            L_ALF1(LL,L1)=L_G1L(LL)+L_G2L(LL)
            L_ALF2(LL,L1)=L_G3L(LL)+L_G4L(LL)
            L_ALF3(LL,L1)=L_G3L(LL)-L_G4L(LL)
            L_ALF4(LL,L1)=L_G1L(LL)-L_G2L(LL)
            L_BET1(LL,L1)=L_G5LR(LL)*2D0
            L_BET2(LL,L1)=L_G5LI(LL)*2D0
         ENDDO

         LMAX = L

!  Examine convergence

         IF ( DABS(G1L).LT.1D-6 ) Calculating = .false. 

      ENDDO        !  Formerly : 300 CONTINUE

!  Finish

      RETURN
END SUBROUTINE GSP_plus


SUBROUTINE MATR_plus                         &
     ( MAXNPA, NPNA, LMAX, NLIN,             & ! Inputs
       A1, A2, A3, A4, B1, B2,               & ! Inputs
       L_A1, L_A2, L_A3, L_A4, L_B1, L_B2,   & ! Inputs
       FMATRIX, L_FMATRIX )                    ! Outputs

!    CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION COEFFICIENTS
 
!    A1,...,B2 - EXPANSION COEFFICIENTS
!    LMAX - NUMBER OF COEFFICIENTS MINUS 1
!    N - NUMBER OF SCATTERING ANGLES
!        THE CORRESPONDING SCATTERING ANGLES ARE GIVEN BY
!        180*(I-1)/(N-1) (DEGREES), WHERE I NUMBERS THE ANGLES

!   use tmat_parameters, only : NPL

   implicit none

!  input arguments
!  ---------------

!  Number of angles

   integer  , intent(in)  :: MAXNPA, NPNA

!  Number of coefficients and derivatives

   integer  , intent(in)  :: LMAX, NLIN

!  Coefficients

   REAL(fpk), intent(in)  :: A1(NPL),A2(NPL),A3(NPL)
   REAL(fpk), intent(in)  :: A4(NPL),B1(NPL),B2(NPL)

!  Linearized Coefficients

   REAL(fpk), intent(in)  :: L_A1(NPL,3),L_A2(NPL,3),L_A3(NPL,3)
   REAL(fpk), intent(in)  :: L_A4(NPL,3),L_B1(NPL,3),L_B2(NPL,3)

!  Outputs
!  -------

!  F-matrix elements
!    Indexing : 1 = 11, 2 = 12, 3 = 22, 4 = 33, 5 = 34, 6 = 44

   REAL(fpk), intent(out)  :: FMATRIX   (MAXNPA,6)
   REAL(fpk), intent(out)  :: L_FMATRIX (MAXNPA,6,3)

!  Local variables
!  ---------------

   integer    :: N, L, L1, L1MAX, I1, qq
   REAL(fpk)  :: DN, DA, DB, DL, DL1, D6, TB, TAA, U
   REAL(fpk)  :: F2, F3, F11, F12, F22, F33, F34, F44
   REAL(fpk)  :: P, P1, P2, P3, P4, PL1, PL2, PL3, PL4, PP1, PP2, PP3, PP4
   REAL(fpk)  :: L_F2(3),  L_F3(3),  L_F11(3), L_F12(3)
   REAL(fpk)  :: L_F22(3), L_F33(3), L_F34(3), L_F44(3)

!  Format sections moved elsewhere, No longer required here

!      PRINT 1000
!      PRINT 1001
!      DO 10 L1=1,L1MAX
!         L=L1-1
!         PRINT 1002,L,A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
!   10 CONTINUE

!      PRINT 1000
!      PRINT 1003
!      DO 500 I1=1,N
! ...............................
!         PRINT 1004,TB,F11,F22,F33,F44,F12,F34
!  500 CONTINUE

! 1000 FORMAT(' ')
! 1001 FORMAT(' ',2X,'S',6X,'ALPHA1',6X,'ALPHA2',6X,'ALPHA3',6X,'ALPHA4',7X,'BETA1',7X,'BETA2')
! 1002 FORMAT(' ',I3,6F12.5)
! 1003 FORMAT(' ',5X,'<',8X,'F11',8X,'F22',8X,'F33',8X,'F44',8X,'F12',8X,'F34')
! 1004 FORMAT(' ',F6.2,6F11.4)

!  Initial

      N=NPNA
      DN=1D0/DBLE(N-1)
      DA=DACOS(-1D0)*DN
      DB=180D0*DN
      L1MAX=LMAX+1

      TB=-DB
      TAA=-DA
      D6=DSQRT(6D0)*0.25D0

!  Start

      DO 500 I1=1,N
         TAA=TAA+DA
         TB=TB+DB
         U=DCOS(TAA)

         F11=0D0   ; F2=0D0    ; F3=0D0    ; F44=0D0    ; F12=0D0    ; F34=0D0
         L_F11=0D0 ; L_F2=0D0  ; L_F3=0D0  ; L_F44=0D0  ; L_F12=0D0  ; L_F34=0D0

         P1=0D0   ; P2=0D0  ; P3=0D0  ; P4=0D0 
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)

         DO 400 L1=1,L1MAX
            L=L1-1
            DL=DBLE(L)
            DL1=DBLE(L1)

            F11=F11+A1(L1)*PP1 ; F44=F44+A4(L1)*PP1
            do qq = 1, nlin
              L_F11(qq) = L_F11(qq) + L_A1(L1,qq)*PP1
              L_F44(qq) = L_F44(qq) + L_A4(L1,qq)*PP1
            enddo

            IF(L.EQ.LMAX) GO TO 350

            PL1=DBLE(2*L+1)
            P=(PL1*U*PP1-DL*P1)/DL1
            P1=PP1
            PP1=P

  350       IF(L.LT.2) GO TO 400

            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            F12=F12+B1(L1)*PP4
            F34=F34+B2(L1)*PP4
            do qq = 1, nlin
              L_F2(qq)  = L_F2(qq)  + (L_A2(L1,qq)+L_A3(L1,qq))*PP2
              L_F3(qq)  = L_F3(qq)  + (L_A2(L1,qq)-L_A3(L1,qq))*PP3
              L_F12(qq) = L_F12(qq) + L_B1(L1,qq)*PP4
              L_F34(qq) = L_F34(qq) + L_B2(L1,qq)*PP4
            enddo

            IF(L.EQ.LMAX) GO TO 400

            PL2=DBLE(L*L1)*U
            PL3=DBLE(L1*(L*L-4))
            PL4=1D0/DBLE(L*(L1*L1-4))
            P=(PL1*(PL2-4D0)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4D0)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DBLE(L*L-4))*P4)/DSQRT(DBLE(L1*L1-4))
            P4=PP4
            PP4=P

  400    CONTINUE

         F22=(F2+F3)*0.5D0
         F33=(F2-F3)*0.5D0
         do qq = 1, nlin
            L_F22(qq)  = ( L_F2(qq) + L_F3(qq) ) * 0.5d0
            L_F33(qq)  = ( L_F2(qq) - L_F3(qq) ) * 0.5d0
         enddo

!  Assign; keep same ordering as in F77 GISS code

         FMATRIX(I1,1) = F11
         FMATRIX(I1,2) = F22
         FMATRIX(I1,3) = F33
         FMATRIX(I1,4) = F44
         FMATRIX(I1,5) = F12
         FMATRIX(I1,6) = F34

         do qq = 1, nlin
            L_FMATRIX(I1,1,qq) = L_F11(qq)
            L_FMATRIX(I1,2,qq) = L_F22(qq)
            L_FMATRIX(I1,3,qq) = L_F33(qq)
            L_FMATRIX(I1,4,qq) = L_F44(qq)
            L_FMATRIX(I1,5,qq) = L_F12(qq)
            L_FMATRIX(I1,6,qq) = L_F34(qq)
         enddo

!  Finish angle loop

  500 CONTINUE

!  Finish

      RETURN
END SUBROUTINE MATR_plus

!  Finish module

end module tmat_scattering_plus
