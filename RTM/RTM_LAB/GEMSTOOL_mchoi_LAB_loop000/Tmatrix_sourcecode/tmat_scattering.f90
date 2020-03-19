 module tmat_scattering

   use tmat_parameters, only : fpk => tmat_fpkind, &
                               cpk => tmat_cpkind, &
                               NPN1, NPN4, NPN6, NPL, NPL1

!  LIST OF SUBROUTINES
!  ===================

!    GSP
!    CCG, CCGIN, FACT, SIGNUM, DIRECT

!    HOVENR

!    MATR

!  Everything PUBLIC here
!  ----------------------

public

contains


SUBROUTINE GSP ( NMAX, LAM, CSCA,                  & ! Inputs
       TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22,    & ! Inputs
       ALF1,ALF2,ALF3,ALF4,BET1,BET2,LMAX,         & ! Outputs
       fail, message, trace )                        ! Outputs

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
 
   use tmat_parameters, only : NPN1, NPN4, NPN6, NPL, NPL1

   implicit none

!  input arguments
!  ---------------

!  Numbers

   integer     , intent(in)  :: nmax
   real(fpk), intent(in)  :: CSCA, LAM

!  Tmatrix input ( Kind = 4 )

!   real(kind=4), intent(in) :: TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4)
!   real(kind=4), intent(in) :: TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4)
!   real(kind=4), intent(in) :: TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4)
!   real(kind=4), intent(in) :: TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)

   real(fpk), intent(in) :: TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4)
   real(fpk), intent(in) :: TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4)
   real(fpk), intent(in) :: TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4)
   real(fpk), intent(in) :: TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)

!  Output arguments
!  ----------------

!  Number of coefficients

   integer     , intent(out) :: LMAX

!  Coefficients

   REAL(fpk), intent(out) :: ALF1(NPL),ALF2(NPL),ALF3(NPL)
   REAL(fpk), intent(out) :: ALF4(NPL),BET1(NPL),BET2(NPL)

!  Exception handling

   logical       , intent(out) :: fail
   character*(*) , intent(out) :: message
   character*(*) , intent(out) :: trace

!  Local arrays
!  ------------

!  FACT and SIGNUM output

   real(fpk) :: SSIGN(900),FCALC(900)

!  Miscellaneous. ( Kind = 8 )

   real(fpk) :: SSI(NPL),SSJ(NPN1)
   real(fpk) :: TR1(NPL1,NPN4),TR2(NPL1,NPN4)
   real(fpk) :: TI1(NPL1,NPN4),TI2(NPL1,NPN4)
   real(fpk) :: G1(NPL1,NPN6),G2(NPL1,NPN6)
   real(fpk) :: FR(NPN4,NPN4),FI(NPN4,NPN4),FF(NPN4,NPN4)

!  A arrays. ( Kind = 8 )

   real(fpk) :: AR1(NPN4),AR2(NPN4),AI1(NPN4),AI2(NPN4)

!  B arrays. ( Kind = 4 )

!   real(kind=4) :: B1R(NPL1,NPL1,NPN4),B1I(NPL1,NPL1,NPN4)
!   real(kind=4) :: B2R(NPL1,NPL1,NPN4),B2I(NPL1,NPL1,NPN4)
   real(fpk) :: B1R(NPL1,NPL1,NPN4),B1I(NPL1,NPL1,NPN4)
   real(fpk) :: B2R(NPL1,NPL1,NPN4),B2I(NPL1,NPL1,NPN4)

!  D arrays. ( Kind = 4 )

!   real(kind=4) :: D1(NPL1,NPN4,NPN4),D2(NPL1,NPN4,NPN4)
!   real(kind=4) :: D3(NPL1,NPN4,NPN4),D4(NPL1,NPN4,NPN4)
!   real(kind=4) :: D5R(NPL1,NPN4,NPN4),D5I(NPL1,NPN4,NPN4)

   real(fpk) :: D1(NPL1,NPN4,NPN4),D2(NPL1,NPN4,NPN4)
   real(fpk) :: D3(NPL1,NPN4,NPN4),D4(NPL1,NPN4,NPN4)
   real(fpk) :: D5R(NPL1,NPN4,NPN4),D5I(NPL1,NPN4,NPN4)

!  Complex array. ( Kind = 16 )

   COMPLEX(cpk) :: CIM(NPN1)
!   COMPLEX(kind=16) :: CIM(NPN1)

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

   REAL(fpk) :: SIG, SSS, SSS1, SI, SJ, SL, T1, T2, T3, T4
   REAL(fpk) :: TT1, TT2, TT3, TT4, TT5, TT6, TT7, TT8
   REAL(fpk) :: AAR1, AAR2, AAI1, AAI2, RR1, RI1, RR2, RI2
   REAL(fpk) :: BBR1, BBR2, BBI1, BBI2, DK, FFN
   REAL(fpk) :: X1, X2, X3, X4, X5, X6, X7, X8, XR, XI, XX
   REAL(fpk) :: DD1, DD2, DD3, DD4, DD5R, DD5I
   REAL(fpk) :: DM1, DM2, DM3, DM4, DM5R, DM5I
   REAL(fpk) :: G1L, G2L, G3L, G4L, G5LR, G5LI

!  Others

   INTEGER      :: N, NN, NNN, N1, NN1, NNMIN, NNMAX, NN1MIN, NN1MAX, NMAX1 
   INTEGER      :: L1MAX, M, M1, MMAX, MMIN, M1MAX, M1MIN, M2, NL, NIFORT
   INTEGER      :: I, I1, J, K1, K2, K3, K4, K5, K6, KN, L, L1
   LOGICAL      :: Calculating

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
               L1=NPN6+M
               TT1=TR11(M1,N,NN)
               TT2=TR12(M1,N,NN)
               TT3=TR21(M1,N,NN)
               TT4=TR22(M1,N,NN)
               TT5=TI11(M1,N,NN)
               TT6=TI12(M1,N,NN)
               TT7=TI21(M1,N,NN)
               TT8=TI22(M1,N,NN)
               T1=TT1+TT2
               T2=TT3+TT4
               T3=TT5+TT6
               T4=TT7+TT8
               TR1(L1,NN)=T1+T2
               TR2(L1,NN)=T1-T2
               TI1(L1,NN)=T3+T4
               TI2(L1,NN)=T3-T4
               IF (M.NE.0) THEN
                  L1=NPN6-M
                  T1=TT1-TT2
                  T2=TT3-TT4
                  T3=TT5-TT6
                  T4=TT7-TT8
                  TR1(L1,NN)=T1-T2
                  TR2(L1,NN)=T1+T2
                  TI1(L1,NN)=T3-T4
                  TI2(L1,NN)=T3+T4
               ENDIF

    6       CONTINUE
   10    CONTINUE

!  debug 3/30/11
!         do nn1 = 1, 19
!            do L1 = 79,83
!               write(67,*)nn1,nn1,l1,n,TR1(L1,NN1)
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
               trace = 'First Call in GSP to CCG, Loop 40 A calculation'
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
               AAR1=0D0
               AAR2=0D0
               AAI1=0D0
               AAI2=0D0
               DO 13 M1=NPN6,M1MAX
                  M=M1-NPN6
                  SSS=G1(M1,NNN)
                  RR1=TR1(M1,NN)
                  RI1=TI1(M1,NN)
                  RR2=TR2(M1,NN)
                  RI2=TI2(M1,NN)
                  IF(M.NE.0) THEN
                     M2=NPN6-M
                     RR1=RR1+TR1(M2,NN)*SIG
                     RI1=RI1+TI1(M2,NN)*SIG
                     RR2=RR2+TR2(M2,NN)*SIG
                     RI2=RI2+TI2(M2,NN)*SIG
                  ENDIF
                  AAR1=AAR1+SSS*RR1
                  AAI1=AAI1+SSS*RI1
                  AAR2=AAR2+SSS*RR2
                  AAI2=AAI2+SSS*RI2
   13          CONTINUE
               XR=FR(NN,N)
               XI=FI(NN,N)
               AR1(NN)=AAR1*XR-AAI1*XI
               AI1(NN)=AAR1*XI+AAI1*XR
               AR2(NN)=AAR2*XR-AAI2*XI
               AI2(NN)=AAR2*XI+AAI2*XR

!               write(67,*)nn,nn1,m1,n,AI2(NN)

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
            DO 30 M1=M1MIN,M1MAX
               BBR1=0D0
               BBI1=0D0
               BBR2=0D0
               BBI2=0D0
               DO 25 NN=NNMIN,NNMAX
                  NNN=NN+1
                  SSS=G2(M1,NNN)
                  BBR1=BBR1+SSS*AR1(NN)
                  BBI1=BBI1+SSS*AI1(NN)
                  BBR2=BBR2+SSS*AR2(NN)
                  BBI2=BBI2+SSS*AI2(NN)
   25          CONTINUE
               B1R(NN1,M1,N)=BBR1
               B1I(NN1,M1,N)=BBI1
               B2R(NN1,M1,N)=BBR2
               B2I(NN1,M1,N)=BBI2

!               write(67,*)m1,nn1,n,n,B1R(NN1,M1,N),&
!               B2R(NN1,M1,N),B1I(NN1,M1,N),B2I(NN1,M1,N)

   30       CONTINUE

!  Finish 40 loop

   40    CONTINUE

!  Finish 100 loop

  100 CONTINUE

!      pause'B-arrays Regular 100'

!  CALCULATION OF THE ARRAYS D1,D2,D3,D4, AND D5
!  =============================================
 
      DO 200 N=1,NMAX
         DO 190 NN=1,NMAX
            M1=MIN0(N,NN)
            M1MAX=NPN6+M1
            M1MIN=NPN6-M1
            NN1MAX=NMAX1+MIN0(N,NN)
            DO 180 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD1=0D0
               DD2=0D0
               DO 150 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B1R(NN1,M1,NN)
                  X4=B1I(NN1,M1,NN)
                  X5=B2R(NN1,M1,N)
                  X6=B2I(NN1,M1,N)
                  X7=B2R(NN1,M1,NN)
                  X8=B2I(NN1,M1,NN)
                  DD1=DD1+XX*(X1*X3+X2*X4)
                  DD2=DD2+XX*(X5*X7+X6*X8)
  150          CONTINUE
               D1(M1,NN,N)=DD1
               D2(M1,NN,N)=DD2
  180       CONTINUE
            MMAX=MIN0(N,NN+2)
            MMIN=MAX0(-N,-NN+2)
            M1MAX=NPN6+MMAX
            M1MIN=NPN6+MMIN
            DO 186 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD3=0D0
               DD4=0D0
               DD5R=0D0
               DD5I=0D0
               M2=-M+2+NPN6
               DO 183 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B2R(NN1,M1,N)
                  X4=B2I(NN1,M1,N)
                  X5=B1R(NN1,M2,NN)
                  X6=B1I(NN1,M2,NN)
                  X7=B2R(NN1,M2,NN)
                  X8=B2I(NN1,M2,NN)
                  DD3=DD3+XX*(X1*X5+X2*X6)
                  DD4=DD4+XX*(X3*X7+X4*X8)
                  DD5R=DD5R+XX*(X3*X5+X4*X6)
                  DD5I=DD5I+XX*(X4*X5-X3*X6)
  183          CONTINUE
               D3(M1,NN,N)=DD3
               D4(M1,NN,N)=DD4
               D5R(M1,NN,N)=DD5R
               D5I(M1,NN,N)=DD5I

!               write(67,*)m1,nn,n,n,D1(M1,NN,N),D2(M1,NN,N),&
!               D3(M1,NN,N),D4(M1,NN,N),D5R(M1,NN,N),D5I(M1,NN,N)

  186       CONTINUE
  190    CONTINUE
  200 CONTINUE

!      pause'D-arrays Regular 200'

!  CALCULATION OF THE EXPANSION COEFFICIENTS
!  ========================================= 

!  Normalization
 
      DK = LAM*LAM/(4D0*CSCA*DACOS(-1D0))

!  Initialize

      Calculating = .true.
      L1 = 0

!  Main loop (Formerly the 300 loop)

!      DO 300 L1=1,L1MAX      (Former Code)
      DO WHILE (Calculating .and. L1.LT.L1MAX)
         L1 = L1 + 1

!  Start

         G1L=0D0
         G2L=0D0
         G3L=0D0
         G4L=0D0
         G5LR=0D0
         G5LI=0D0
         L=L1-1
         SL=SSI(L1)*DK

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
               DM1=0D0
               DM2=0D0

               DO 270 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  IF(M.GE.0) SSS1=G1(M1,NNN)
                  IF(M.LT.0) SSS1=G1(NPN6-M,NNN)*SI
                  DM1=DM1+SSS1*D1(M1,NN,N)
                  DM2=DM2+SSS1*D2(M1,NN,N)
  270          CONTINUE

               FFN=FF(NN,N)
               SSS=G1(NPN6+1,NNN)*FFN
               G1L=G1L+SSS*DM1
               G2L=G2L+SSS*DM2*SI

               IF(L.GE.2) THEN
                  DM3=0D0
                  DM4=0D0
                  DM5R=0D0
                  DM5I=0D0
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
  275             CONTINUE
                  G5LR=G5LR-SSS*DM5R
                  G5LI=G5LI-SSS*DM5I
                  SSS=G2(NPN4,NNN)*FFN
                  G3L=G3L+SSS*DM3
                  G4L=G4L+SSS*DM4*SI
               ENDIF

  280       CONTINUE

  290    CONTINUE

!  Final

         G1L=G1L*SL
         G2L=G2L*SL
         G3L=G3L*SL
         G4L=G4L*SL
         G5LR=G5LR*SL
         G5LI=G5LI*SL

!        write(67,*)L,L,L1,L1,G1L,G2L,G3L,G4L,G5LR,G5LI

         ALF1(L1)=G1L+G2L
         ALF2(L1)=G3L+G4L
         ALF3(L1)=G3L-G4L
         ALF4(L1)=G1L-G2L
         BET1(L1)=G5LR*2D0
         BET2(L1)=G5LI*2D0
         LMAX = L

!  debug for Validation agaianst F77 code
!         write(57,'(2i4,1p6e15.5)')&
!         L1,LMAX,ALF1(L1),ALF2(L1),ALF3(L1),ALF4(L1),BET1(L1),BET2(L1)

!  Examine convergence

         IF ( DABS(G1L).LT.1D-6 ) Calculating = .false. 

      ENDDO        !  Formerly : 300 CONTINUE

!  Finish

      RETURN
END SUBROUTINE GSP


SUBROUTINE FACT (FCALC)

!   CALCULATION OF THE QUANTITIES F(N+1)=0.5*LN(N!)
!   0.LE.N.LE.899
 
   implicit none

!  Argument

   REAL(fpk), intent(out) :: FCALC(900)

!  Local

   Integer :: I, I1

!  Code

   FCALC(1)=0D0
   FCALC(2)=0D0
   DO 2 I=3,900
      I1=I-1
      !FCALC(I)=FCALC(I1)+0.5D0*DLOG(DFLOAT(I1))
      FCALC(I)=FCALC(I1)+0.5D0*DLOG(REAL(I1,fpk))
 2 CONTINUE

!  Finish

   RETURN
END SUBROUTINE FACT 
 
SUBROUTINE SIGNUM(SSIGN)

!   CALCULATION OF THE ARRAY SSIGN(N+1)=SIGN(N)
!   0.LE.N.LE.899
 
   implicit none

!  Argument

   REAL(fpk), intent(out) :: SSIGN(900)

!  Local

   Integer :: N

!  Code

   SSIGN(1)=1D0
   DO 2 N=2,899 
      SSIGN(N)=-SSIGN(N-1)
 2 CONTINUE
   SSIGN(900)=0D0

!  Finish

   RETURN
END SUBROUTINE SIGNUM
 
 
SUBROUTINE CCG ( N, N1, NMAX, K1, K2, FCALC, SSIGN, &
                 GG, fail, message)

!   CALCULATION OF CLEBSCH-GORDAN COEFFICIENTS
!   (N,M:N1,M1/NN,MM)
!   FOR GIVEN N AND N1. M1=MM-M, INDEX MM IS FOUND FROM M AS
!   MM=M*K1+K2
!
!   INPUT PARAMETERS :  N,N1,NMAX,K1,K2
!                               N.LE.NMAX
!                               N.GE.1
!                               N1.GE.0
!                               N1.LE.N+NMAX
!   OUTPUT PARAMETERS : GG(M+NPN6,NN+1) - ARRAY OF THE CORRESPONDING
!                                       COEFFICIENTS
!                               /M/.LE.N
!                               /M1/=/M*(K1-1)+K2/.LE.N1
!                               NN.LE.MIN(N+N1,NMAX)
!                               NN.GE.MAX(/MM/,/N-N1/)
!   IF K1=1 AND K2=0, THEN 0.LE.M.LE.N
 
   use tmat_parameters, only : NPL1, NPN5, NPN6

   implicit none

!  input arguments
!  ---------------

   Integer   , intent(in)    :: NMAX, K1, K2
   Integer   , intent(inout) :: N, N1
   REAL (fpk), intent(in)    :: FCALC(900),SSIGN(900)

!  Outputs
!  -------

   REAL (fpk)   , intent(out)  :: GG(NPL1,NPN6)
   LOGICAL      , intent(out)  :: fail
   CHARACTER*(*), intent(out)  :: message

!  Local variables
!  ---------------

   REAL (fpk) :: CD(0:NPN5),CU(0:NPN5)
   REAL (fpk) :: A, B, C, D, C1, C2
   INTEGER    :: M, MM, MIN, MF, M1, MIND
   INTEGER    :: NN, NNF, NNL, NNM, NNU 
   Logical    :: CHECK

!  Initialize

   GG = 0.0d0
   fail    = .false.
   message = ' '

!  Check

   CHECK =  (NMAX.LE.NPN4.AND.0.LE.N1      &
            .AND.N1.LE.NMAX+N.AND.N.GE.1   &
            .AND.N.LE.NMAX)
   IF ( .not. CHECK ) THEN
      fail = .true.
      message = ' CHECK ERROR AT START of SUBROUTINE CCG'
      return
   endif

!  Code. For now, has been left untouched, ugly Goto constructs and all

      NNF=MIN0(N+N1,NMAX)
      MIN=NPN6-N
      MF=NPN6+N
      IF(K1.EQ.1.AND.K2.EQ.0) MIN=NPN6

      DO 100 MIND=MIN,MF
         M=MIND-NPN6
         MM=M*K1+K2
         M1=MM-M

!          CU = 0D0 ; CD = 0D0

         IF(IABS(M1).GT.N1) GO TO 90
         NNL=MAX0(IABS(MM),IABS(N-N1))
         IF(NNL.GT.NNF) GO TO 90
         NNU=N+N1
         NNM=(NNU+NNL)*0.5D0
         IF (NNU.EQ.NNL) NNM=NNL
         CALL CCGIN(N,N1,M,MM,FCALC,SSIGN,C,fail,message)
         if ( fail ) return
         CU(NNL)=C  
         IF (NNL.EQ.NNF) GO TO 50
         C2=0D0
         C1=C
         DO 7 NN=NNL+1,MIN0(NNM,NNF)
            A=DBLE((NN+MM)*(NN-MM)*(N1-N+NN))
            A=A*DBLE((N-N1+NN)*(N+N1-NN+1)*(N+N1+NN+1))
            A=DBLE(4*NN*NN)/A
            A=A*DBLE((2*NN+1)*(2*NN-1))
            A=DSQRT(A)
            B=0.5D0*DBLE(M-M1)
            D=0D0
            IF(NN.EQ.1) GO TO 5
            B=DBLE(2*NN*(NN-1))
            B=DBLE((2*M-MM)*NN*(NN-1)-MM*N*(N+1)+MM*N1*(N1+1))/B
            D=DBLE(4*(NN-1)*(NN-1))
            D=D*DBLE((2*NN-3)*(2*NN-1))
            D=DBLE((NN-MM-1)*(NN+MM-1)*(N1-N+NN-1))/D
            D=D*DBLE((N-N1+NN-1)*(N+N1-NN+2)*(N+N1+NN))
            D=DSQRT(D)
    5       C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CU(NN)=C
    7    CONTINUE
         IF (NNF.LE.NNM) GO TO 50
         CALL DIRECT(N,M,N1,M1,FCALC,C)
         CD(NNU)=C
         IF (NNU.EQ.NNM+1) GO TO 50
         C2=0D0
         C1=C
         DO 12 NN=NNU-1,NNM+1,-1
            A=DBLE((NN-MM+1)*(NN+MM+1)*(N1-N+NN+1))
            A=A*DBLE((N-N1+NN+1)*(N+N1-NN)*(N+N1+NN+2))
            A=DBLE(4*(NN+1)*(NN+1))/A
            A=A*DBLE((2*NN+1)*(2*NN+3))
            A=DSQRT(A)
            B=DBLE(2*(NN+2)*(NN+1))
            B=DBLE((2*M-MM)*(NN+2)*(NN+1)-MM*N*(N+1) +MM*N1*(N1+1))/B
            D=DBLE(4*(NN+2)*(NN+2))
            D=D*DBLE((2*NN+5)*(2*NN+3))
            D=DBLE((NN+MM+2)*(NN-MM+2)*(N1-N+NN+2))/D
            D=D*DBLE((N-N1+NN+2)*(N+N1-NN-1)*(N+N1+NN+3))
            D=DSQRT(D)
            C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CD(NN)=C
   12    CONTINUE

   50    DO 9 NN=NNL,NNF
            IF (NN.LE.NNM) GG(MIND,NN+1)=CU(NN)
            IF (NN.GT.NNM) GG(MIND,NN+1)=CD(NN)
    9    CONTINUE

   90    CONTINUE
  100 CONTINUE

   RETURN
END SUBROUTINE CCG
 
SUBROUTINE DIRECT (N,M,N1,M1,F,C)

   implicit none

!  inputs

   Integer      , intent(in)   :: N,M,N1,M1
   REAL (fpk), intent(in)   :: F(900)
   REAL (fpk), intent(out)  :: C

!  Code

   C=F(2*N+1)+F(2*N1+1)+F(N+N1+M+M1+1)+F(N+N1-M-M1+1)    
   C=C-F(2*(N+N1)+1)-F(N+M+1)-F(N-M+1)-F(N1+M1+1)-F(N1-M1+1)
   C=DEXP(C)

!  Finish

   RETURN
END SUBROUTINE DIRECT
 
SUBROUTINE CCGIN ( N, N1, M, MM, F, SSIGN, & ! Inputs
                   G, fail, message )        ! Outputs

!   CALCULATION OF THE CLEBCSH-GORDAN COEFFICIENTS
!   G=(N,M:N1,MM-M/NN,MM)
!   FOR GIVEN N,N1,M,MM, WHERE NN=MAX(/MM/,/N-N1/)
!                               /M/.LE.N
!                               /MM-M/.LE.N1
!                               /MM/.LE.N+N1
 
   implicit none

!  inputs

   Integer   , intent(inout)   :: N, N1, M, MM
   REAL (fpk), intent(in)      :: F(900),SSIGN(900)

!  Outputs

   REAL(fpk)    , intent(out)  :: G
   LOGICAL      , intent(out)  :: fail
   CHARACTER*(*), intent(out)  :: message

!  Local

   Integer    :: M1, L1, L2, L3, N2, M2, K, N12, M12
   REAL (fpk) :: A
   Logical    :: CHECK

!  Initialize

   G = 0.0d0
   fail    = .false.
   message = ' '

!  Check

   M1 = MM - M
   CHECK = (N.GE.IABS(M).AND.N1.GE.IABS(M1).AND.IABS(MM).LE.(N+N1))
   IF ( .not. CHECK ) THEN
      fail = .true.
      message = ' ERROR IN SUBROUTINE CCGIN'
      return
   endif

!  Code

   IF (IABS(MM).GT.IABS(N-N1)) then
      A=1D0
      L1=M
      L2=MM
      IF (.not.(MM.GE.0)) then
         MM=-MM
         M=-M
         M1=-M1
         A=SSIGN(MM+N+N1+1)
      ENDIF
      G=A*SSIGN(N+M+1)                                          &
         *DEXP(F(2*MM+2)+F(N+N1-MM+1)+F(N+M+1)+F(N1+M1+1)       &
              -F(N+N1+MM+2)-F(N-N1+MM+1)-F(-N+N1+MM+1)-F(N-M+1) &
              -F(N1-M1+1))
      M=L1
      MM=L2
   ELSE
      L1=N
      L2=N1
      L3=M
      IF (.not.(N1.LE.N)) then
         K=N
         N=N1
         N1=K
         K=M
         M=M1
         M1=K
      ENDIF
      N2=N*2
      M2=M*2
      N12=N1*2
      M12=M1*2
      G=SSIGN(N1+M1+1)                                        &
      *DEXP(F(N+M+1)+F(N-M+1)+F(N12+1)+F(N2-N12+2)-F(N2+2)    &
            -F(N1+M1+1)-F(N1-M1+1)-F(N-N1+MM+1)-F(N-N1-MM+1))
      N=L1
      N1=L2
      M=L3
   ENDIF

!  Finish

   RETURN
END SUBROUTINE CCGIN

  
SUBROUTINE HOVENR ( L1, A1, A2, A3, A4, B1, B2,             & ! Inputs
                    fail_1, fail_2, message_1, message_2 )    ! Outputs

   use tmat_parameters, only : NPL

   implicit none

!  input arguments
!  ---------------

!  Number of coefficients

   integer  , intent(in)  :: L1

!  Coefficients

   REAL(fpk), intent(in)  :: A1(NPL),A2(NPL),A3(NPL)
   REAL(fpk), intent(in)  :: A4(NPL),B1(NPL),B2(NPL)

!  Output arguments
!  ----------------

!  Exception handling

   logical       , intent(out) :: fail_1
   logical       , intent(out) :: fail_2
   character*(*) , intent(out) :: message_1
   character*(*) , intent(out) :: message_2

!  Local variables
!  ---------------

   integer        :: KONTR, L, LL, I
   REAL(fpk)   :: DL, DDL, AA1, AA2, AA3, AA4, BB1, BB2
   REAL(fpk)   :: C, CC, C1, C2, C3
   character*4    :: C4
   character*9    :: C9

!  Initialize output

   fail_1 = .false.
   fail_2 = .false.
   message_1 = ' '
   message_2 = ' '

!  Loop over moments

      DO 100 L=1,L1

         KONTR=1
         LL=L-1
         DL = DBLE(LL)*2D0+1D0
         DDL=DL*0.48D0
         AA1=A1(L)
         AA2=A2(L)
         AA3=A3(L)
         AA4=A4(L)
         BB1=B1(L)
         BB2=B2(L)

!  First test

         IF(LL.GE.1.AND.DABS(AA1).GE.DL) KONTR=2
         IF(DABS(AA2).GE.DL) KONTR=2
         IF(DABS(AA3).GE.DL) KONTR=2
         IF(DABS(AA4).GE.DL) KONTR=2
         IF(DABS(BB1).GE.DDL) KONTR=2
         IF(DABS(BB2).GE.DDL) KONTR=2
         IF(KONTR.EQ.2) THEN
            write(c4,'(I4)')LL
            Message_1 = 'First test of VdM/Hovenier NOT Satisfied, L ='//C4
            Fail_1 = .true.   ; RETURN
         endif

!  Second test

         C=-0.1D0
         DO 50 I=1,11
            C=C+0.1D0
            CC=C*C
            C1=CC*BB2*BB2
            C2=C*AA4
            C3=C*AA3
            IF((DL-C*AA1)*(DL-C*AA2)-CC*BB1*BB1.LE.-1D-4) KONTR=2
            IF((DL-C2)*(DL-C3)+C1.LE.-1D-4) KONTR=2
            IF((DL+C2)*(DL-C3)-C1.LE.-1D-4) KONTR=2
            IF((DL-C2)*(DL+C3)-C1.LE.-1D-4) KONTR=2
            IF(KONTR.EQ.2) THEN
               write(c4,'(I4)')LL
               write(c9,'(d9.2)')C
               Message_2 = 'Second test of VdM/Hovenier NOT Satisfied, L ='//C4//', A = '//C9
               Fail_2 = .true.   ; RETURN
            endif
   50    CONTINUE

  100 CONTINUE

!  Finish

      RETURN
END SUBROUTINE HOVENR


SUBROUTINE MATR                                    &
     ( MAXNPA, NPNA, LMAX, A1, A2, A3, A4, B1, B2, & ! Inputs
       FMATRIX )                                     ! Outputs

!    CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION COEFFICIENTS
 
!    A1,...,B2 - EXPANSION COEFFICIENTS
!    LMAX - NUMBER OF COEFFICIENTS MINUS 1
!    N - NUMBER OF SCATTERING ANGLES
!        THE CORRESPONDING SCATTERING ANGLES ARE GIVEN BY
!        180*(I-1)/(N-1) (DEGREES), WHERE I NUMBERS THE ANGLES

   use tmat_parameters, only : NPL

   implicit none

!  input arguments
!  ---------------

!  Number of angles

   integer  , intent(in)  :: MAXNPA, NPNA

!  Number of coefficients

   integer  , intent(in)  :: LMAX

!  Coefficients

   REAL(fpk), intent(in)  :: A1(NPL),A2(NPL),A3(NPL)
   REAL(fpk), intent(in)  :: A4(NPL),B1(NPL),B2(NPL)

!  Outputs
!  -------

!  F-matrix elements
!    Indexing : 1 = 11, 2 = 12, 3 = 22, 4 = 33, 5 = 34, 6 = 44

   REAL(fpk), intent(out)  :: FMATRIX (MAXNPA,6)

!  Local variables
!  ---------------

   integer    :: N, L, L1, L1MAX, I1
   REAL(fpk)  :: DN, DA, DB, DL, DL1, D6, TB, TAA, U
   REAL(fpk)  :: F2, F3, F11, F12, F22, F33, F34, F44, P
   REAL(fpk)  :: P1, P2, P3, P4, PL1, PL2, PL3, PL4, PP1, PP2, PP3, PP4
 
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
         F11=0D0
         F2=0D0
         F3=0D0
         F44=0D0
         F12=0D0
         F34=0D0
         P1=0D0
         P2=0D0
         P3=0D0
         P4=0D0
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)
         DO 400 L1=1,L1MAX
            L=L1-1
            DL=DBLE(L)
            DL1=DBLE(L1)
            F11=F11+A1(L1)*PP1
            F44=F44+A4(L1)*PP1
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

!  Assign; keep same ordering as in F77 GISS code

         FMATRIX(I1,1) = F11
         FMATRIX(I1,2) = F22
         FMATRIX(I1,3) = F33
         FMATRIX(I1,4) = F44
         FMATRIX(I1,5) = F12
         FMATRIX(I1,6) = F34

!  Finish angle loop

  500 CONTINUE

!  Finish

      RETURN
END SUBROUTINE MATR

!  Finish module

end module tmat_scattering
