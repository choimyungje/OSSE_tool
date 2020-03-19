module tmat_makers_plus

   use tmat_parameters, only : fpk => tmat_fpkind, &
                                NPN1, NPN2, NPNG2

   use lapack_tools
   use tmat_makers, only : VIG

!  LIST OF SUBROUTINES
!  ===================

!    tmatrix_R_plus
!    tmatrix_R0_plus

!       RGQ_Qinv_plus
!       VIG_plus

private  RGQ_Qinv_plus, VIG_plus
public   tmatrix_R_plus, tmatrix_R0_plus

!  Everything PUBLIC here
!  ----------------------

public

contains

subroutine tmatrix_R0_plus                               &
     ( Do_LinearEps, NGAUSS, NMAX, NCHECK, NLIN,         & ! Inputs
       X, W, Le_X, Le_W, AN, ANN, PPI, PIR, PII,         & ! Inputs
       R, DR, DDR, DRR, DRI,                             & ! Inputs
        j_bess,  y_bess,  jr_bess,  ji_bess,             & ! Inputs
       dj_bess, dy_bess, djr_bess, dji_bess,             & ! Inputs
       Le_R, Le_DR, L_DDR, L_DRR, L_DRI,                 & ! Inputs
       L_j_bess,  L_y_bess,  L_jr_bess,  L_ji_bess,      & ! Inputs
       L_dj_bess, L_dy_bess, L_djr_bess, L_dji_bess,     & ! Inputs
       R12, R21, I12, I21, RG12, RG21, IG12, IG21,       & ! Outputs
       L_R12, L_R21, L_I12, L_I21,                       & ! Outputs
       L_RG12, L_RG21, L_IG12, L_IG21,                   & ! Outputs
       TR1, TI1, L_TR1, L_TI1, fail, message, trace )      ! Outputs

   implicit none

!  input arguments
!  ---------------

!  Flag for linear eps

   logical  , intent(in)  :: DO_LinearEps

!  Control

   integer  , intent(in)  :: nmax, ncheck, ngauss, nlin

!  X and W

   REAL(fpk), intent(in)  :: X(NPNG2), W(NPNG2)
   REAL(fpk), intent(in)  :: Le_X(NPNG2), Le_W(NPNG2)

!  Constants

   REAL(fpk), intent(in)  :: AN(NPN1),ANN(NPN1,NPN1)
   real(fpk), intent(in)  :: PPI, PIR, PII

!  R and DR functions and linearizations

   REAL(fpk), intent(in)  :: R(NPNG2), DR(NPNG2)
   real(fpk), intent(in)  :: DDR(NPNG2), DRR(NPNG2), DRI(NPNG2)

   real(fpk), intent(in)  :: Le_R(NPNG2), Le_DR(NPNG2)
   real(fpk), intent(in)  :: L_DDR(3,NPNG2), L_DRR(3,NPNG2), L_DRI(3,NPNG2)

!  Bessel Master output

   real(fpk), intent(in)  :: J_BESS  (NPNG2,NPN1)
   real(fpk), intent(in)  :: Y_BESS  (NPNG2,NPN1)
   real(fpk), intent(in)  :: JR_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: JI_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DJ_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DY_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DJR_BESS(NPNG2,NPN1)
   real(fpk), intent(in)  :: DJI_BESS(NPNG2,NPN1)

!  Linearized Bessel output

   real(fpk), intent(in)  :: L_J_BESS  (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_Y_BESS  (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_JR_BESS (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_JI_BESS (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_DJ_BESS (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_DY_BESS (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_DJR_BESS(3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_DJI_BESS(3,NPNG2,NPN1)

!  Output arguments
!  ----------------

!  Miscellaneous output (Tmatrix)

   REAL(fpk), intent(out) :: R12 (NPN1,NPN1), R21 (NPN1,NPN1)
   REAL(fpk), intent(out) :: I12 (NPN1,NPN1), I21 (NPN1,NPN1)
   REAL(fpk), intent(out) :: RG12(NPN1,NPN1), RG21(NPN1,NPN1)
   REAL(fpk), intent(out) :: IG12(NPN1,NPN1), IG21(NPN1,NPN1)

   REAL(fpk), intent(out) :: L_R12 (3,NPN1,NPN1), L_R21 (3,NPN1,NPN1)
   REAL(fpk), intent(out) :: L_I12 (3,NPN1,NPN1), L_I21 (3,NPN1,NPN1)
   REAL(fpk), intent(out) :: L_RG12(3,NPN1,NPN1), L_RG21(3,NPN1,NPN1)
   REAL(fpk), intent(out) :: L_IG12(3,NPN1,NPN1), L_IG21(3,NPN1,NPN1)

!  T-matrix output

   REAL(fpk), intent(out) :: TR1(NPN2,NPN2), TI1(NPN2,NPN2)
   REAL(fpk), intent(out) :: L_TR1(3,NPN2,NPN2), L_TI1(3,NPN2,NPN2)

!  Exception handling

   logical       , intent(out) :: fail
   character*(*) , intent(out) :: message
   character*(*) , intent(out) :: trace

!  Local variables
!  ---------------

!  Miscellaneous Local arrays

   REAL(fpk)  :: RR(NPNG2), Le_RR(NPNG2), SIG(NPN2)

!  DV1/DV2 and D1/D2 matrices and linearizations (Wigner functions)

   REAL(fpk)  :: DV1(NPN1), DV2(NPN1)
   REAL(fpk)  :: Le_DV1(NPN1),Le_DV2(NPN1)
   REAL(fpk)  :: D1(NPNG2,NPN1), D2(NPNG2,NPN1)
   REAL(fpk)  :: Le_D1(NPNG2,NPN1), Le_D2(NPNG2,NPN1)

!  RGQ and Q matrices

   REAL(fpk)  ::   QR(NPN2,NPN2),  QI(NPN2,NPN2)
   REAL(fpk)  :: RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)

!  Linearized RGQ and Q matrices

   REAL(fpk)  ::   L_QR(3,NPN2,NPN2),  L_QI(3,NPN2,NPN2)
   REAL(fpk)  ::   L_RGQR(3,NPN2,NPN2),L_RGQI(3,NPN2,NPN2)

!  Other forms of the same thing

   REAL(fpk)  ::   TQR(NPN2,NPN2),  TQI(NPN2,NPN2)
   REAL(fpk)  :: TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
   REAL(fpk)  ::   L_TQR(3,NPN2,NPN2),  L_TQI(3,NPN2,NPN2)
   REAL(fpk)  :: L_TRGQR(3,NPN2,NPN2),L_TRGQI(3,NPN2,NPN2)

!  Help variables

   LOGICAL   :: DO_200
   INTEGER   :: MM1, NNMAX, NG, NGSS, NM
   INTEGER   :: I, I1, I2, N, N1, N2, K1, KK1, K2, KK2, P

   REAL(fpk) :: FACTOR, SI, DD1, DD2, DDRI, DRRI, DRII
   REAL(fpk) :: AN1, AN2, AN12, A12, A21, A22, AA1
   REAL(fpk) :: URI, RRI, F1, F21, F22, TPIR, TPII, TPPI

   REAL(fpk) :: AR12, AR21, AI12, AI21
   REAL(fpk) :: GR12, GR21, GI12, GI21

   REAL(fpk) :: TAR12, TAR21, TAI12, TAI21
   REAL(fpk) :: TGR12, TGR21, TGI12, TGI21

   REAL(fpk) :: D1N1, D2N1, D1N2, D2N2
   REAL(fpk) :: QJ1, QY1, QJR2, QJI2
   REAL(fpk) :: QDJR2, QDJI2, QDJ1, QDY1

   REAL(fpk) :: B1R, B1I, B2R, B2I, B3R, B3I, B4R, B4I, B5R, B5I
   REAL(fpk) :: C1R, C1I, C2R, C2I, C3R, C3I, C4R, C4I, C5R, C5I

!  Linearized help variables

   REAL(fpk) :: L_DDRI, L_DRRI, L_DRII
   REAL(fpk) :: Le_A12, Le_A21, Le_A22, Le_AA1
   REAL(fpk) :: Le_URI, Le_RRI, Le_F1, Le_F21, Le_F22, L_TPIR(3), L_TPII(3)

   REAL(fpk) :: L_AR12(3), L_AR21(3), L_AI12(3), L_AI21(3)
   REAL(fpk) :: L_GR12(3), L_GR21(3), L_GI12(3), L_GI21(3)

   REAL(fpk) :: L_TAR12, L_TAR21, L_TAI12, L_TAI21
   REAL(fpk) :: L_TGR12, L_TGR21, L_TGI12, L_TGI21

   REAL(fpk) :: Le_D1N1, Le_D2N1, Le_D1N2, Le_D2N2
   REAL(fpk) :: L_QJ1(3), L_QY1(3), L_QJR2(3), L_QJI2(3)
   REAL(fpk) :: L_QDJR2(3), L_QDJI2(3), L_QDJ1(3), L_QDY1(3)

   REAL(fpk) :: L_B1R, L_B1I, L_B2R, L_B2I, L_B3R, L_B3I, L_B4R, L_B4I, L_B5R, L_B5I
   REAL(fpk) :: L_C1R, L_C1I, L_C2R, L_C2I, L_C3R, L_C3I, L_C4R, L_C4I, L_C5R, L_C5I

!  Initialize
!  ----------

!  Exception handling

   TRACE = ' '

!  Numbers

   MM1   = 1
   NNMAX = NMAX+NMAX
   NG    = 2*NGAUSS
   NGSS  = NG

!  Factor

   FACTOR = 1.0D0
   IF (NCHECK.EQ.1) THEN
      NGSS   = NGAUSS
      FACTOR = 2.0D0
   ENDIF

!  Signs

   SI=1.0D0
   DO N=1,NNMAX
      SI     = -SI
      SIG(N) = SI
   ENDDO

!  R times W

   DO I = 1, NGSS
      RR(I) = W(I) * R(I)
   ENDDO

   if ( Do_LinearEps ) then
      DO I = 1, NGSS
         Le_RR(I) = Le_W(I) * R(I) + W(I) * Le_R(I)
      ENDDO
   endif

!  VIG_plus call
!  -------------

   if ( Do_LinearEps ) then
      DO I=1,NGAUSS
         I1=NGAUSS+I       ;  I2=NGAUSS-I+1
         call VIG_PLUS ( NMAX, X(I1), Le_X(I1), 0, DV1, DV2, Le_DV1, Le_DV2 )
         DO N=1,NMAX
            SI=SIG(N)   
            DD1=DV1(N)     ;  D1(I1,N)=DD1     ;  D1(I2,N)= DD1*SI 
            DD2=DV2(N)     ;  D2(I1,N)=DD2     ;  D2(I2,N)=-DD2*SI
            DD1=Le_DV1(N)  ;  Le_D1(I1,N)=DD1  ;  Le_D1(I2,N)= DD1*SI 
            DD2=Le_DV2(N)  ;  Le_D2(I1,N)=DD2  ;  Le_D2(I2,N)=-DD2*SI
        ENDDO
      ENDDO
   else
      Le_D1 = 0.0d0     ; Le_D2 = 0.0d0
      DO I=1,NGAUSS
         I1=NGAUSS+I       ;  I2=NGAUSS-I+1
         call VIG ( NMAX, X(I1), 0, DV1, DV2 )
         DO  N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)   ;   D1(I1,N)=DD1   ;   D1(I2,N)= DD1*SI 
            DD2=DV2(N)   ;   D2(I1,N)=DD2   ;   D2(I2,N)=-DD2*SI
        ENDDO
      ENDDO
   endif

!  Main calculation loop
!  ---------------------

   DO N1=MM1,NMAX
      AN1=AN(N1)
      DO N2=MM1,NMAX
         AN2=AN(N2)

!  Initialize

         AR12=0.0D0     ; AR21=0.0D0    ; AI12=0.0D0    ; AI21=0.0D0
         GR12=0.0D0     ; GR21=0.0D0    ; GI12=0.0D0    ; GI21=0.0D0
         L_AR12=0.0D0   ; L_AR21=0.0D0  ; L_AI12=0.0D0  ; L_AI21=0.0D0
         L_GR12=0.0D0   ; L_GR21=0.0D0  ; L_GI12=0.0D0  ; L_GI21=0.0D0
           
!  Avoidance condition

         DO_200 = .true.
         IF ( NCHECK.EQ.1.AND.SIG(N1+N2).LT.0.0D0) DO_200 = .false.

!  Do this loop (200)

         IF ( DO_200 ) THEN

            DO I=1,NGSS

!  Initial assign 1

               D1N1 = D1(I,N1)   ;  D2N1 = D2(I,N1)
               D1N2 = D1(I,N2)   ;  D2N2 = D2(I,N2)
               A12  = D1N1*D2N2  ;  A21  = D2N1*D1N2
               A22  = D2N1*D2N2  ;  AA1 = A12+A21

               if ( Do_LinearEps ) then
                  Le_D1N1 = Le_D1(I,N1)  ;  Le_D2N1 = Le_D2(I,N1)
                  Le_D1N2 = Le_D1(I,N2)  ;  Le_D2N2 = Le_D2(I,N2)
                  Le_A12  = Le_D1N1*D2N2 + D1N1*Le_D2N2
                  Le_A21  = Le_D2N1*D1N2 + D2N1*Le_D1N2
                  Le_A22  = Le_D2N1*D2N2 + D2N1*Le_D2N2
                  Le_AA1  = Le_A12 + Le_A21
               endif
 
!  Initial assign 2

               QJ1   = J_BESS(I,N1)    ;  QY1   = Y_BESS(I,N1)
               QJR2  = JR_BESS(I,N2)   ;  QJI2  = JI_BESS(I,N2)
               QDJR2 = DJR_BESS(I,N2)  ;  QDJI2 = DJI_BESS(I,N2)
               QDJ1  = DJ_BESS(I,N1)   ;  QDY1  = DY_BESS(I,N1)

               DO P = 1, nlin
                  L_QJ1(P)   = L_J_BESS(P,I,N1)    ;  L_QY1(P)   = L_Y_BESS(P,I,N1)
                  L_QJR2(P)  = L_JR_BESS(P,I,N2)   ;  L_QJI2(P)  = L_JI_BESS(P,I,N2)
                  L_QDJR2(P) = L_DJR_BESS(P,I,N2)  ;  L_QDJI2(P) = L_DJI_BESS(P,I,N2)
                  L_QDJ1(P)  = L_DJ_BESS(P,I,N1)   ;  L_QDY1(P)  = L_DY_BESS(P,I,N1)
               enddo         

!  Start Regular code

               C1R=QJR2*QJ1
               C1I=QJI2*QJ1
               B1R=C1R-QJI2*QY1
               B1I=C1I+QJR2*QY1

               C2R=QJR2*QDJ1
               C2I=QJI2*QDJ1
               B2R=C2R-QJI2*QDY1
               B2I=C2I+QJR2*QDY1

               DDRI=DDR(I)
               C3R=DDRI*C1R
               C3I=DDRI*C1I
               B3R=DDRI*B1R
               B3I=DDRI*B1I

               C4R=QDJR2*QJ1
               C4I=QDJI2*QJ1
               B4R=C4R-QDJI2*QY1
               B4I=C4I+QDJR2*QY1

               DRRI=DRR(I)
               DRII=DRI(I)
               C5R=C1R*DRRI-C1I*DRII
               C5I=C1I*DRRI+C1R*DRII
               B5R=B1R*DRRI-B1I*DRII
               B5I=B1I*DRRI+B1R*DRII

!              if ( i.lt.5.and.N1.eq.2.and.n2.eq.2)write(*,*)i,c5r,b5r
               URI=DR(I)
               RRI=RR(I)

               F1=RRI*A22
               F21=RRI*URI*AN1*A12
               AR12=AR12+F1*B2R+F21*B3R
               AI12=AI12+F1*B2I+F21*B3I
               GR12=GR12+F1*C2R+F21*C3R
               GI12=GI12+F1*C2I+F21*C3I

               F22=RRI*URI*AN2*A21
               AR21=AR21+F1*B4R+F22*B5R
               AI21=AI21+F1*B4I+F22*B5I
               GR21=GR21+F1*C4R+F22*C5R
               GI21=GI21+F1*C4I+F22*C5I

!  Linearized section
!    There are extra contributions for P = 3 (the EPS linearization)

               do P = 1, nlin

                  L_C1R  = L_QJR2(P)*QJ1 + QJR2*L_QJ1(P)  
                  L_C1I  = L_QJI2(P)*QJ1 + QJI2*L_QJ1(P)
                  L_B1R  = L_C1R - L_QJI2(P)*QY1 - QJI2*L_QY1(P)
                  L_B1I  = L_C1I + L_QJR2(P)*QY1 + QJR2*L_QY1(P)

                  L_C2R  = L_QJR2(P)*QDJ1 + QJR2*L_QDJ1(P)  
                  L_C2I  = L_QJI2(P)*QDJ1 + QJI2*L_QDJ1(P)
                  L_B2R  = L_C2R - L_QJI2(P)*QDY1 - QJI2*L_QDY1(P)
                  L_B2I  = L_C2I + L_QJR2(P)*QDY1 + QJR2*L_QDY1(P)

                  L_DDRI = L_DDR(P,I)
                  L_C3R  = L_DDRI*C1R + DDRI*L_C1R
                  L_C3I  = L_DDRI*C1I + DDRI*L_C1I
                  L_B3R  = L_DDRI*B1R + DDRI*L_B1R
                  L_B3I  = L_DDRI*B1I + DDRI*L_B1I

                  L_C4R  = L_QDJR2(P)*QJ1 + QDJR2*L_QJ1(P)  
                  L_C4I  = L_QDJI2(P)*QJ1 + QDJI2*L_QJ1(P)
                  L_B4R  = L_C4R - L_QDJI2(P)*QY1 - QDJI2*L_QY1(P)
                  L_B4I  = L_C4I + L_QDJR2(P)*QY1 + QDJR2*L_QY1(P)

                  L_DRRI = L_DRR(P,I)
                  L_DRII = L_DRI(P,I)
                  L_C5R  = L_C1R*DRRI-L_C1I*DRII + C1R*L_DRRI-C1I*L_DRII
                  L_C5I  = L_C1I*DRRI+L_C1R*DRII + C1I*L_DRRI+C1R*L_DRII
                  L_B5R  = L_B1R*DRRI-L_B1I*DRII + B1R*L_DRRI-B1I*L_DRII
                  L_B5I  = L_B1I*DRRI+L_B1R*DRII + B1I*L_DRRI+B1R*L_DRII

                  if ( P .eq. 3 ) then
                     Le_URI = Le_DR(I)
                     Le_RRI = Le_RR(I)
                     Le_F1  = Le_RRI*A22 + RRI*Le_A22
                     Le_F21 = AN1 * ( Le_RRI*URI*A12 + RRI*Le_URI*A12 + RRI*URI*Le_A12 )
                     Le_F22 = AN2 * ( Le_RRI*URI*A21 + RRI*Le_URI*A21 + RRI*URI*Le_A21 )
                  endif

                  L_AR12(P) = L_AR12(P) + F1*L_B2R + F21*L_B3R
                  L_AI12(P) = L_AI12(P) + F1*L_B2I + F21*L_B3I
                  L_GR12(P) = L_GR12(P) + F1*L_C2R + F21*L_C3R
                  L_GI12(P) = L_GI12(P) + F1*L_C2I + F21*L_C3I
                  L_AR21(P) = L_AR21(P) + F1*L_B4R + F22*L_B5R
                  L_AI21(P) = L_AI21(P) + F1*L_B4I + F22*L_B5I
                  L_GR21(P) = L_GR21(P) + F1*L_C4R + F22*L_C5R
                  L_GI21(P) = L_GI21(P) + F1*L_C4I + F22*L_C5I

                  IF ( P .eq. 3 )  then
                     L_AR12(P) = L_AR12(P) + Le_F1*B2R + Le_F21*B3R
                     L_AI12(P) = L_AI12(P) + Le_F1*B2I + Le_F21*B3I
                     L_GR12(P) = L_GR12(P) + Le_F1*C2R + Le_F21*C3R
                     L_GI12(P) = L_GI12(P) + Le_F1*C2I + Le_F21*C3I
                     L_AR21(P) = L_AR21(P) + Le_F1*B4R + Le_F22*B5R
                     L_AI21(P) = L_AI21(P) + Le_F1*B4I + Le_F22*B5I
                     L_GR21(P) = L_GR21(P) + Le_F1*C4R + Le_F22*C5R
                     L_GI21(P) = L_GI21(P) + Le_F1*C4I + Le_F22*C5I
                  ENDIF

               ENDDO

            ENDDO
         ENDIF

!  Form the R/I and RG/IG matrices, and their linearizations

         AN12=ANN(N1,N2)*FACTOR
         R12(N1,N2)  = AR12*AN12   ;   R21(N1,N2)  = AR21*AN12
         I12(N1,N2)  = AI12*AN12   ;   I21(N1,N2)  = AI21*AN12
         RG12(N1,N2) = GR12*AN12   ;   RG21(N1,N2) = GR21*AN12
         IG12(N1,N2) = GI12*AN12   ;   IG21(N1,N2) = GI21*AN12
         DO P = 1, NLIN
            L_R12(P,N1,N2)  = L_AR12(P)*AN12   ;   L_R21(P,N1,N2)  = L_AR21(P)*AN12
            L_I12(P,N1,N2)  = L_AI12(P)*AN12   ;   L_I21(P,N1,N2)  = L_AI21(P)*AN12
            L_RG12(P,N1,N2) = L_GR12(P)*AN12   ;   L_RG21(P,N1,N2) = L_GR21(P)*AN12
            L_IG12(P,N1,N2) = L_GI12(P)*AN12   ;   L_IG21(P,N1,N2) = L_GI21(P)*AN12
         ENDDO

!         P = 1
!         write(57,'(2i4,1p8e20.10)')N1,N2,L_R12(P,N1,N2),R12(N1,N2), & 
!                                          L_R21(P,N1,N2),R21(N1,N2), &
!                                          L_I12(P,N1,N2),I12(N1,N2), &
!                                          L_I21(P,N1,N2),I21(N1,N2)
!         write(57,'(2i4,1p8e20.10)')N1,N2,L_RG12(P,N1,N2),RG12(N1,N2), & 
!                                          L_RG21(P,N1,N2),RG21(N1,N2), &
!                                          L_IG12(P,N1,N2),IG12(N1,N2), &
!                                          L_IG21(P,N1,N2),IG21(N1,N2)

      ENDDO
   ENDDO

!   pause'develop'

!  Develop RGQ and Q matrices
!  --------------------------

!  Extended precision ??? (Why the T prefix in this section)

!  constants

   TPIR = PIR 
   TPII = PII
   TPPI = PPI   ! has no linearizations !!!!!!!!!!!!
   L_TPIR(1) = TPIR   ; L_TPIR(2) = 0.0d0  ; L_TPIR(3) = 0.0d0 
   L_TPII(1) = 0.0d0  ; L_TPII(2) = TPII   ; L_TPII(3) = 0.0d0

!  Loops

   NM=NMAX
   DO N1=MM1,NMAX
      K1=N1-MM1+1      ;  KK1=K1+NM
      DO N2=MM1,NMAX
         K2=N2-MM1+1   ;  KK2=K2+NM
 
!  Regular code

         TAR12 =  I12(N1,N2)   ;  TAI12 = -R12(N1,N2)
         TGR12 =  IG12(N1,N2)  ;  TGI12 = -RG12(N1,N2)
         TAR21 = -I21(N1,N2)   ;  TAI21 = R21(N1,N2)
         TGR21 = -IG21(N1,N2)  ;  TGI21 = RG21(N1,N2)

         TQR(K1,K2)   = TPIR*TAR21-TPII*TAI21+TPPI*TAR12
         TQI(K1,K2)   = TPIR*TAI21+TPII*TAR21+TPPI*TAI12
         TRGQR(K1,K2) = TPIR*TGR21-TPII*TGI21+TPPI*TGR12
         TRGQI(K1,K2) = TPIR*TGI21+TPII*TGR21+TPPI*TGI12

         TQR(K1,KK2)   = 0.0D0  ;  TQI(K1,KK2)   = 0.0D0
         TRGQR(K1,KK2) = 0.0D0  ;  TRGQI(K1,KK2) = 0.0D0
         TQR(KK1,K2)   = 0.0D0  ;  TQI(KK1,K2)   = 0.0D0
         TRGQR(KK1,K2) = 0.0D0  ;  TRGQI(KK1,K2) = 0.0D0

         TQR(KK1,KK2)   = TPIR*TAR12-TPII*TAI12+TPPI*TAR21
         TQI(KK1,KK2)   = TPIR*TAI12+TPII*TAR12+TPPI*TAI21
         TRGQR(KK1,KK2) = TPIR*TGR12-TPII*TGI12+TPPI*TGR21
         TRGQI(KK1,KK2) = TPIR*TGI12+TPII*TGR12+TPPI*TGI21

!  Linearized code

         DO P = 1, NLIN

            L_TAR12 =  L_I12(P,N1,N2)   ;  L_TAI12 = -L_R12(P,N1,N2)
            L_TGR12 =  L_IG12(P,N1,N2)  ;  L_TGI12 = -L_RG12(P,N1,N2)
            L_TAR21 = -L_I21(P,N1,N2)   ;  L_TAI21 = L_R21(P,N1,N2)
            L_TGR21 = -L_IG21(P,N1,N2)  ;  L_TGI21 = L_RG21(P,N1,N2)

            L_TQR(P,K1,K2)   = L_TPIR(P)*TAR21 - L_TPII(P)*TAI21 + TPPI*L_TAR12  &
                               + TPIR*L_TAR21  - TPII* L_TAI21
            L_TQI(P,K1,K2)   = L_TPIR(P)*TAI21 + L_TPII(P)*TAR21 + TPPI*L_TAI12  &
                               + TPIR*L_TAI21  + TPII*L_TAR21
            L_TRGQR(P,K1,K2) = L_TPIR(P)*TGR21 - L_TPII(P)*TGI21 + TPPI*L_TGR12  &
                               + TPIR*L_TGR21  - TPII*L_TGI21
            L_TRGQI(P,K1,K2) = L_TPIR(P)*TGI21 + L_TPII(P)*TGR21 + TPPI*L_TGI12  &
                               + TPIR*L_TGI21  + TPII*L_TGR21

            L_TQR(P,K1,KK2)   = 0.0D0  ;  L_TQI(P,K1,KK2)   = 0.0D0
            L_TRGQR(P,K1,KK2) = 0.0D0  ;  L_TRGQI(P,K1,KK2) = 0.0D0
            L_TQR(P,KK1,K2)   = 0.0D0  ;  L_TQI(P,KK1,K2)   = 0.0D0
            L_TRGQR(P,KK1,K2) = 0.0D0  ;  L_TRGQI(P,KK1,K2) = 0.0D0

            L_TQR(P,KK1,KK2)   = L_TPIR(P)*TAR12 - L_TPII(P)*TAI12 + TPPI*L_TAR21  &
                                 + TPIR*L_TAR12  - TPII* L_TAI12
            L_TQI(P,KK1,KK2)   = L_TPIR(P)*TAI12 + L_TPII(P)*TAR12 + TPPI*L_TAI21  &
                                 + TPIR*L_TAI12  + TPII*L_TAR12
            L_TRGQR(P,KK1,KK2) = L_TPIR(P)*TGR12 - L_TPII(P)*TGI12 + TPPI*L_TGR21  &
                                 + TPIR*L_TGR12  - TPII*L_TGI12
            L_TRGQI(P,KK1,KK2) = L_TPIR(P)*TGI12 + L_TPII(P)*TGR12 + TPPI*L_TGI21  &
                                 + TPIR*L_TGI12  + TPII*L_TGR12

         ENDDO

!  End main loops

      ENDDO
   ENDDO
 
!  Copy to Normal RGQ and Q arrays (Why ????)

   NNMAX = 2 * NM
   DO N1 = 1,NNMAX
      DO N2=1,NNMAX
         QR(N1,N2)=TQR(N1,N2)
         QI(N1,N2)=TQI(N1,N2)
         RGQR(N1,N2)=TRGQR(N1,N2)
         RGQI(N1,N2)=TRGQI(N1,N2)
         DO P = 1, NLIN
            L_QR(P,N1,N2)=L_TQR(P,N1,N2)
            L_QI(P,N1,N2)=L_TQI(P,N1,N2)
            L_RGQR(P,N1,N2)=L_TRGQR(P,N1,N2)
            L_RGQI(P,N1,N2)=L_TRGQI(P,N1,N2)
         ENDDO
!         P = 1
!         write(57,'(2i4,1p4e20.10)')n1,n2,L_QR(P,N1,N2),L_QI(P,N1,N2),&
!                                          QR(N1,N2),QI(N1,N2)
!         write(57,'(2i4,1p4e20.10)')n1,n2,L_RGQR(P,N1,N2),L_RGQI(P,N1,N2),&
!                                          RGQR(N1,N2),RGQI(N1,N2)
     ENDDO
   ENDDO

!  Call T-matrix calculator (uses LAPACK inversion)

   call RGQ_Qinv_plus                               &
          ( NMAX, NLIN, QR, QI, RGQR, RGQI,         & ! Inputs
            L_QR, L_QI, L_RGQR, L_RGQI,             & ! Inputs
            TR1, TI1, L_TR1, L_TI1, fail, message )   ! Outputs

!  debug

!        P = 1
!       DO N1=1,NNMAX
!           DO N2=1,NNMAX
!         write(57,'(2i4,1p4e20.10)')n1,n2,L_TR1(P,N1,N2),L_TI1(P,N1,N2),&
!                                          TR1(N1,N2),TI1(N1,N2)
!         write(57,'(2i4,1p2e20.9)')n1,n2,TR1(N1,N2),TI1(N1,N2)
!       enddo
!       enddo
!      pause'57'

!  Exception handling

   if ( fail ) then
      TRACE = ' RGQ_Qinv_plus failed in Tmatrix_R0_plus'
   endif

!  Finish

   return
end subroutine tmatrix_R0_plus


subroutine tmatrix_R_plus                                       &
     ( Do_LinearEps, M, NGAUSS, NMAX, NCHECK, NLIN,             & ! Inputs
       X, W, Le_X, Le_W, AN, ANN, S, SS, Le_S, Le_SS,           & ! Inputs
       PPI, PIR, PII, R, DR, DDR, DRR, DRI,                     & ! Inputs
        j_bess,  y_bess,  jr_bess,  ji_bess,                    & ! Inputs
       dj_bess, dy_bess, djr_bess, dji_bess,                    & ! Inputs
       Le_R, Le_DR, L_DDR, L_DRR, L_DRI,                        & ! Inputs
       L_j_bess,  L_y_bess,  L_jr_bess,  L_ji_bess,             & ! Inputs
       L_dj_bess, L_dy_bess, L_djr_bess, L_dji_bess,            & ! Inputs
       R11, R12, R21, R22, I11, I12, I21, I22,                  & ! Outputs
       RG11,RG12,RG21,RG22,IG11,IG12,IG21,IG22,                 & ! Outputs
       L_R11, L_R12, L_R21, L_R22, L_I11, L_I12, L_I21, L_I22,  & ! Outputs
       L_RG11,L_RG12,L_RG21,L_RG22,L_IG11,L_IG12,L_IG21,L_IG22, & ! Outputs
       TR1, TI1, L_TR1, L_TI1, fail, message, trace )             ! Outputs

   implicit none

!  input arguments
!  ---------------

!  Flag for linear eps

   logical  , intent(in)  :: DO_LinearEps

!  Control

   integer  , intent(in)  :: nmax, ncheck, ngauss, m, nlin

!  X and W and S/SS, plus derivatives

   REAL(fpk), intent(in)  :: X(NPNG2), W(NPNG2)
   REAL(fpk), intent(in)  :: S(NPNG2),SS(NPNG2)
   REAL(fpk), intent(in)  :: Le_X(NPNG2), Le_W(NPNG2)
   REAL(fpk), intent(in)  :: Le_S(NPNG2), Le_SS(NPNG2)

!  Constants

   REAL(fpk), intent(in)  :: AN(NPN1),ANN(NPN1,NPN1)
   real(fpk), intent(in)  :: PPI, PIR, PII

!  R and DR functions and linearizations

   REAL(fpk), intent(in)  :: R(NPNG2), DR(NPNG2)
   real(fpk), intent(in)  :: DDR(NPNG2), DRR(NPNG2), DRI(NPNG2)

   real(fpk), intent(in)  :: Le_R(NPNG2), Le_DR(NPNG2)
   real(fpk), intent(in)  :: L_DDR(3,NPNG2), L_DRR(3,NPNG2), L_DRI(3,NPNG2)

!  Bessel Master output

   real(fpk), intent(in)  :: J_BESS  (NPNG2,NPN1)
   real(fpk), intent(in)  :: Y_BESS  (NPNG2,NPN1)
   real(fpk), intent(in)  :: JR_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: JI_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DJ_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DY_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DJR_BESS(NPNG2,NPN1)
   real(fpk), intent(in)  :: DJI_BESS(NPNG2,NPN1)

!  Linearized Bessel output

   real(fpk), intent(in)  :: L_J_BESS  (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_Y_BESS  (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_JR_BESS (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_JI_BESS (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_DJ_BESS (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_DY_BESS (3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_DJR_BESS(3,NPNG2,NPN1)
   real(fpk), intent(in)  :: L_DJI_BESS(3,NPNG2,NPN1)

!  Output arguments
!  ----------------

!  Miscellaneous output

   REAL(fpk), intent(out) :: R11(NPN1,NPN1),R12(NPN1,NPN1)
   REAL(fpk), intent(out) :: R21(NPN1,NPN1),R22(NPN1,NPN1)
   REAL(fpk), intent(out) :: I11(NPN1,NPN1),I12(NPN1,NPN1)
   REAL(fpk), intent(out) :: I21(NPN1,NPN1),I22(NPN1,NPN1)

   REAL(fpk), intent(out) :: RG11(NPN1,NPN1),RG12(NPN1,NPN1)
   REAL(fpk), intent(out) :: RG21(NPN1,NPN1),RG22(NPN1,NPN1)
   REAL(fpk), intent(out) :: IG11(NPN1,NPN1),IG12(NPN1,NPN1)
   REAL(fpk), intent(out) :: IG21(NPN1,NPN1),IG22(NPN1,NPN1)

   REAL(fpk), intent(out) :: L_R11 (3,NPN1,NPN1), L_R12 (3,NPN1,NPN1)
   REAL(fpk), intent(out) :: L_R21 (3,NPN1,NPN1), L_R22 (3,NPN1,NPN1)
   REAL(fpk), intent(out) :: L_I11 (3,NPN1,NPN1), L_I12 (3,NPN1,NPN1)
   REAL(fpk), intent(out) :: L_I21 (3,NPN1,NPN1), L_I22 (3,NPN1,NPN1)

   REAL(fpk), intent(out) :: L_RG11(3,NPN1,NPN1), L_RG12(3,NPN1,NPN1)
   REAL(fpk), intent(out) :: L_RG21(3,NPN1,NPN1), L_RG22(3,NPN1,NPN1)
   REAL(fpk), intent(out) :: L_IG11(3,NPN1,NPN1), L_IG12(3,NPN1,NPN1)
   REAL(fpk), intent(out) :: L_IG21(3,NPN1,NPN1), L_IG22(3,NPN1,NPN1)

!  T-matrix output

   REAL(fpk), intent(out) :: TR1(NPN2,NPN2), TI1(NPN2,NPN2)
   REAL(fpk), intent(out) :: L_TR1(3,NPN2,NPN2), L_TI1(3,NPN2,NPN2)

!  Exception handling

   logical       , intent(out) :: fail
   character*(*) , intent(out) :: message
   character*(*) , intent(out) :: trace

!  Local variables
!  ---------------

!  Miscellaneous Local arrays

   REAL(fpk)  :: RR(NPNG2), Le_RR(NPNG2), SIG(NPN2)
   REAL(fpk)  :: DS(NPNG2), DSS(NPNG2)
   REAL(fpk)  :: Le_DS(NPNG2)
   REAL(fpk)  :: Le_DSS(NPNG2)

!  DV1/DV2 and D1/D2 matrices and linearizations (Wigner functions)

   REAL(fpk)  :: DV1(NPN1), DV2(NPN1)
   REAL(fpk)  :: Le_DV1(NPN1),Le_DV2(NPN1)
   REAL(fpk)  :: D1(NPNG2,NPN1), D2(NPNG2,NPN1)
   REAL(fpk)  :: Le_D1(NPNG2,NPN1), Le_D2(NPNG2,NPN1)

!  RGQ and Q matrices

   REAL(fpk)  ::   QR(NPN2,NPN2),  QI(NPN2,NPN2)
   REAL(fpk)  :: RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)

!  Linearized RGQ and Q matrices

   REAL(fpk)  ::   L_QR(3,NPN2,NPN2),  L_QI(3,NPN2,NPN2)
   REAL(fpk)  ::   L_RGQR(3,NPN2,NPN2),L_RGQI(3,NPN2,NPN2)

!  Other forms of the same thing

   REAL(fpk)  ::   TQR(NPN2,NPN2),  TQI(NPN2,NPN2)
   REAL(fpk)  :: TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
   REAL(fpk)  ::   L_TQR(3,NPN2,NPN2),  L_TQI(3,NPN2,NPN2)
   REAL(fpk)  :: L_TRGQR(3,NPN2,NPN2),L_TRGQI(3,NPN2,NPN2)

!  Help variables

   INTEGER   :: MM1, NNMAX, NG, NGSS, NM
   INTEGER   :: I, I1, I2, N, N1, N2, K1, KK1, K2, KK2, P
   REAL(fpk) :: FACTOR, SI, DD1, DD2, DDRI, DRRI, DRII, QM, QMM
   REAL(fpk) :: URI, RRI, DSI, DSSI, WR
   REAL(fpk) :: AN1, AN2, AN12, A11, A12, A21, A22, AA1, AA2
   REAL(fpk) :: E1, E2, E3, EH, F1, F21, F22, TPIR, TPII, TPPI

   REAL(fpk) :: AR11, AR12, AR21, AR22
   REAL(fpk) :: AI11, AI12, AI21, AI22
   REAL(fpk) :: GR11, GR12, GR21, GR22
   REAL(fpk) :: GI11, GI12, GI21, GI22

   REAL(fpk) :: TAR11, TAR12, TAR21, TAR22
   REAL(fpk) :: TAI11, TAI12, TAI21, TAI22
   REAL(fpk) :: TGR11, TGR12, TGR21, TGR22
   REAL(fpk) :: TGI11, TGI12, TGI21, TGI22

   REAL(fpk) :: D1N1, D2N1, D1N2, D2N2
   REAL(fpk) :: QJ1, QY1, QJR2, QJI2
   REAL(fpk) :: QDJR2, QDJI2, QDJ1, QDY1

   REAL(fpk) :: B1R, B1I, B2R, B2I, B3R, B3I, B4R, B4I
   REAL(fpk) :: C1R, C1I, C2R, C2I, C3R, C3I, C4R, C4I
   REAL(fpk) :: B5R, B5I, B6R, B6I, B7R, B7I, B8R, B8I
   REAL(fpk) :: C5R, C5I, C6R, C6I, C7R, C7I, C8R, C8I

!  Linearized help variables

   REAL(fpk) :: L_DDRI, L_DRRI, L_DRII
   REAL(fpk) :: Le_A11, Le_A12, Le_A21, Le_A22, Le_AA1, Le_AA2
   REAL(fpk) :: Le_URI, Le_RRI, Le_DSI, Le_DSSI, L_TPIR(3), L_TPII(3)
   REAL(fpk) :: Le_F1, Le_F21, Le_F22, Le_E1, Le_E2, Le_E3, Le_EH, Le_WR

   REAL(fpk) :: L_AR11(3), L_AR12(3), L_AR21(3), L_AR22(3)
   REAL(fpk) :: L_AI11(3), L_AI12(3), L_AI21(3), L_AI22(3)
   REAL(fpk) :: L_GR11(3), L_GR12(3), L_GR21(3), L_GR22(3)
   REAL(fpk) :: L_GI11(3), L_GI12(3), L_GI21(3), L_GI22(3)

   REAL(fpk) :: L_TAR11, L_TAR12, L_TAR21, L_TAR22
   REAL(fpk) :: L_TAI11, L_TAI12, L_TAI21, L_TAI22
   REAL(fpk) :: L_TGR11, L_TGR12, L_TGR21, L_TGR22
   REAL(fpk) :: L_TGI11, L_TGI12, L_TGI21, L_TGI22

   REAL(fpk) :: Le_D1N1, Le_D2N1, Le_D1N2, Le_D2N2
   REAL(fpk) :: L_QJ1(3), L_QY1(3), L_QJR2(3), L_QJI2(3)
   REAL(fpk) :: L_QDJR2(3), L_QDJI2(3), L_QDJ1(3), L_QDY1(3)

   REAL(fpk) :: L_B1R, L_B1I, L_B2R, L_B2I, L_B3R, L_B3I, L_B4R, L_B4I
   REAL(fpk) :: L_C1R, L_C1I, L_C2R, L_C2I, L_C3R, L_C3I, L_C4R, L_C4I
   REAL(fpk) :: L_B5R, L_B5I, L_B6R, L_B6I, L_B7R, L_B7I, L_B8R, L_B8I
   REAL(fpk) :: L_C5R, L_C5I, L_C6R, L_C6I, L_C7R, L_C7I, L_C8R, L_C8I

!  Initialize
!  ----------

!  Exception handling

   TRACE = ' '

!  Numbers

   MM1   = M
   QM    = DBLE(M)
   QMM   = QM*QM

   NM    = NMAX+NMAX
   NG    = 2*NGAUSS
   NGSS  = NG

!  Factor

   FACTOR = 1.0D0
   IF (NCHECK.EQ.1) THEN
      NGSS   = NGAUSS
      FACTOR = 2.0D0
   ENDIF

!  Signs

   SI=1.0D0
   DO N=1,NM
      SI     = -SI
      SIG(N) = SI
   ENDDO

!  R times W, DS and DSS

   DO I = 1, NGSS
      WR     = W(I)  * R(I)
      DS(I)  = S(I)  * QM * WR
      DSS(I) = SS(I) * QMM
      RR(I)  = WR
   ENDDO

   if ( Do_LinearEps ) then
      DO I = 1, NGSS
         WR    = W(I) * R(I)                                ! New
         Le_WR = Le_W(I) * R(I) + W(I) * Le_R(I)
         Le_RR(I) = Le_WR
!         Le_DS(I) = S(I)  * QM * Le_WR                     ! Former code 
         Le_DS(I)  = QM  * ( S(I) * Le_WR + Le_S(I) * WR )  ! New
         Le_DSS(I) = QMM * Le_SS(I)                         ! New
      ENDDO
   endif

!  VIG_plus call
!  -------------

   if ( Do_LinearEps ) then
      DO I=1,NGAUSS
         I1=NGAUSS+I       ;  I2=NGAUSS-I+1
         call VIG_PLUS ( NMAX, X(I1), Le_X(I1), M, DV1, DV2, Le_DV1, Le_DV2 )
         DO N=1,NMAX
            SI=SIG(N)   
            DD1=DV1(N)     ;  D1(I1,N)=DD1     ;  D1(I2,N)= DD1*SI 
            DD2=DV2(N)     ;  D2(I1,N)=DD2     ;  D2(I2,N)=-DD2*SI
            DD1=Le_DV1(N)  ;  Le_D1(I1,N)=DD1  ;  Le_D1(I2,N)= DD1*SI 
            DD2=Le_DV2(N)  ;  Le_D2(I1,N)=DD2  ;  Le_D2(I2,N)=-DD2*SI
        ENDDO
      ENDDO
   else
      Le_D1 = 0.0d0     ; Le_D2 = 0.0d0
      DO I=1,NGAUSS
         I1=NGAUSS+I       ;  I2=NGAUSS-I+1
         call VIG ( NMAX, X(I1), M, DV1, DV2 )
         DO  N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)   ;   D1(I1,N)=DD1   ;   D1(I2,N)= DD1*SI 
            DD2=DV2(N)   ;   D2(I1,N)=DD2   ;   D2(I2,N)=-DD2*SI
        ENDDO
      ENDDO
   endif

!  Main calculation loop

   DO N1=MM1,NMAX
      AN1=AN(N1)
      DO N2=MM1,NMAX

         AN2=AN(N2)
         SI=SIG(N1+N2)

!  Initialize

         AR11=0.0D0     ; AR12=0.0D0    ; AR21=0.0D0    ; AR22=0.0D0
         AI11=0.0D0     ; AI12=0.0D0    ; AI21=0.0D0    ; AI22=0.0D0
         GR11=0.0D0     ; GR12=0.0D0    ; GR21=0.0D0    ; GR22=0.0D0
         GI11=0.0D0     ; GI12=0.0D0    ; GI21=0.0D0    ; GI22=0.0D0

         L_AR11=0.0D0   ; L_AR12=0.0D0  ; L_AR21=0.0D0  ; L_AR22=0.0D0
         L_AI11=0.0D0   ; L_AI12=0.0D0  ; L_AI21=0.0D0  ; L_AI22=0.0D0
         L_GR11=0.0D0   ; L_GR12=0.0D0  ; L_GR21=0.0D0  ; L_GR22=0.0D0
         L_GI11=0.0D0   ; L_GI12=0.0D0  ; L_GI21=0.0D0  ; L_GI22=0.0D0

!  200 loop

         DO I=1,NGSS

!  Initial assign 1

            D1N1 = D1(I,N1)   ;  D2N1 = D2(I,N1)
            D1N2 = D1(I,N2)   ;  D2N2 = D2(I,N2)
            A11  = D1N1*D1N2  ;  A12  = D1N1*D2N2
            A21  = D2N1*D1N2  ;  A22  = D2N1*D2N2
            AA1  = A12+A21    ;  AA2  = A11*DSS(I)+A22

            if ( Do_LinearEps ) then
               Le_D1N1 = Le_D1(I,N1)  ;  Le_D2N1 = Le_D2(I,N1)
               Le_D1N2 = Le_D1(I,N2)  ;  Le_D2N2 = Le_D2(I,N2)
               Le_A11  = Le_D1N1*D1N2 + D1N1*Le_D1N2
               Le_A12  = Le_D1N1*D2N2 + D1N1*Le_D2N2
               Le_A21  = Le_D2N1*D1N2 + D2N1*Le_D1N2
               Le_A22  = Le_D2N1*D2N2 + D2N1*Le_D2N2
               Le_AA1  = Le_A12 + Le_A21
!               Le_AA2  = Le_A11*DSS(I) + Le_A22                 ! Former
               Le_AA2  = Le_A11*DSS(I) + A11*Le_DSS(I) + Le_A22  ! New
            endif
 
!  Initial assign 2

            QJ1   = J_BESS(I,N1)    ;  QY1   = Y_BESS(I,N1)
            QJR2  = JR_BESS(I,N2)   ;  QJI2  = JI_BESS(I,N2)
            QDJR2 = DJR_BESS(I,N2)  ;  QDJI2 = DJI_BESS(I,N2)
            QDJ1  = DJ_BESS(I,N1)   ;  QDY1  = DY_BESS(I,N1)

            DO P = 1, nlin
               L_QJ1(P)   = L_J_BESS(P,I,N1)    ;  L_QY1(P)   = L_Y_BESS(P,I,N1)
               L_QJR2(P)  = L_JR_BESS(P,I,N2)   ;  L_QJI2(P)  = L_JI_BESS(P,I,N2)
               L_QDJR2(P) = L_DJR_BESS(P,I,N2)  ;  L_QDJI2(P) = L_DJI_BESS(P,I,N2)
               L_QDJ1(P)  = L_DJ_BESS(P,I,N1)   ;  L_QDY1(P)  = L_DY_BESS(P,I,N1)
            enddo         
 
!  Regular section

            C1R=QJR2*QJ1
            C1I=QJI2*QJ1
            B1R=C1R-QJI2*QY1
            B1I=C1I+QJR2*QY1

            C2R=QJR2*QDJ1
            C2I=QJI2*QDJ1
            B2R=C2R-QJI2*QDY1
            B2I=C2I+QJR2*QDY1

            DDRI=DDR(I)
            C3R=DDRI*C1R
            C3I=DDRI*C1I
            B3R=DDRI*B1R
            B3I=DDRI*B1I

            C4R=QDJR2*QJ1
            C4I=QDJI2*QJ1
            B4R=C4R-QDJI2*QY1
            B4I=C4I+QDJR2*QY1

            DRRI=DRR(I)
            DRII=DRI(I)
            C5R=C1R*DRRI-C1I*DRII
            C5I=C1I*DRRI+C1R*DRII
            B5R=B1R*DRRI-B1I*DRII
            B5I=B1I*DRRI+B1R*DRII

            C6R=QDJR2*QDJ1
            C6I=QDJI2*QDJ1
            B6R=C6R-QDJI2*QDY1
            B6I=C6I+QDJR2*QDY1

            C7R=C4R*DDRI
            C7I=C4I*DDRI
            B7R=B4R*DDRI
            B7I=B4I*DDRI

            C8R=C2R*DRRI-C2I*DRII
            C8I=C2I*DRRI+C2R*DRII
            B8R=B2R*DRRI-B2I*DRII
            B8I=B2I*DRRI+B2R*DRII

            URI=DR(I)
            DSI=DS(I)
            DSSI=DSS(I)
            RRI=RR(I)

!  Check condition

            IF (NCHECK.EQ.1 ) THEN
               IF ( SI.GT.0.0D0 ) THEN
                  F1=RRI*AA2
                  F21=RRI*URI*AN1*A12
                  AR12=AR12+F1*B2R+F21*B3R
                  AI12=AI12+F1*B2I+F21*B3I
                  GR12=GR12+F1*C2R+F21*C3R
                  GI12=GI12+F1*C2I+F21*C3I
                  F22=RRI*URI*AN2*A21
                  AR21=AR21+F1*B4R+F22*B5R
                  AI21=AI21+F1*B4I+F22*B5I
                  GR21=GR21+F1*C4R+F22*C5R
                  GI21=GI21+F1*C4I+F22*C5I
               ELSE
                  E1=DSI*AA1
                  AR11=AR11+E1*B1R
                  AI11=AI11+E1*B1I
                  GR11=GR11+E1*C1R
                  GI11=GI11+E1*C1I
                  EH=DSI*URI*A11
                  E3=EH*AN2
                  E2=EH*AN1
                  AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                  AI22=AI22+E1*B6I+E2*B7I+E3*B8I
                  GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                  GI22=GI22+E1*C6I+E2*C7I+E3*C8I
               ENDIF
            ENDIF

!  Linearized section
!    There are extra contributions for P = 3 (the EPS linearization)

            do P = 1, nlin

               L_C1R  = L_QJR2(P)*QJ1 + QJR2*L_QJ1(P)  
               L_C1I  = L_QJI2(P)*QJ1 + QJI2*L_QJ1(P)
               L_B1R  = L_C1R - L_QJI2(P)*QY1 - QJI2*L_QY1(P)
               L_B1I  = L_C1I + L_QJR2(P)*QY1 + QJR2*L_QY1(P)

               L_C2R  = L_QJR2(P)*QDJ1 + QJR2*L_QDJ1(P)  
               L_C2I  = L_QJI2(P)*QDJ1 + QJI2*L_QDJ1(P)
               L_B2R  = L_C2R - L_QJI2(P)*QDY1 - QJI2*L_QDY1(P)
               L_B2I  = L_C2I + L_QJR2(P)*QDY1 + QJR2*L_QDY1(P)

               L_DDRI = L_DDR(P,I)
               L_C3R  = L_DDRI*C1R + DDRI*L_C1R
               L_C3I  = L_DDRI*C1I + DDRI*L_C1I
               L_B3R  = L_DDRI*B1R + DDRI*L_B1R
               L_B3I  = L_DDRI*B1I + DDRI*L_B1I

               L_C4R  = L_QDJR2(P)*QJ1 + QDJR2*L_QJ1(P)  
               L_C4I  = L_QDJI2(P)*QJ1 + QDJI2*L_QJ1(P)
               L_B4R  = L_C4R - L_QDJI2(P)*QY1 - QDJI2*L_QY1(P)
               L_B4I  = L_C4I + L_QDJR2(P)*QY1 + QDJR2*L_QY1(P)

               L_DRRI = L_DRR(P,I)
               L_DRII = L_DRI(P,I)
               L_C5R  = L_C1R*DRRI-L_C1I*DRII + C1R*L_DRRI-C1I*L_DRII
               L_C5I  = L_C1I*DRRI+L_C1R*DRII + C1I*L_DRRI+C1R*L_DRII
               L_B5R  = L_B1R*DRRI-L_B1I*DRII + B1R*L_DRRI-B1I*L_DRII
               L_B5I  = L_B1I*DRRI+L_B1R*DRII + B1I*L_DRRI+B1R*L_DRII

               L_C6R  = L_QDJR2(P)*QDJ1 + QDJR2*L_QDJ1(P)  
               L_C6I  = L_QDJI2(P)*QDJ1 + QDJI2*L_QDJ1(P)
               L_B6R  = L_C6R - L_QDJI2(P)*QDY1 - QDJI2*L_QDY1(P)
               L_B6I  = L_C6I + L_QDJR2(P)*QDY1 + QDJR2*L_QDY1(P)

               L_C7R  = L_DDRI*C4R + DDRI*L_C4R
               L_C7I  = L_DDRI*C4I + DDRI*L_C4I
               L_B7R  = L_DDRI*B4R + DDRI*L_B4R
               L_B7I  = L_DDRI*B4I + DDRI*L_B4I

               L_C8R  = L_C2R*DRRI-L_C2I*DRII + C2R*L_DRRI-C2I*L_DRII
               L_C8I  = L_C2I*DRRI+L_C2R*DRII + C2I*L_DRRI+C2R*L_DRII
               L_B8R  = L_B2R*DRRI-L_B2I*DRII + B2R*L_DRRI-B2I*L_DRII
               L_B8I  = L_B2I*DRRI+L_B2R*DRII + B2I*L_DRRI+B2R*L_DRII

               if ( P .eq. 3 ) then
                  Le_URI  = Le_DR(I)
                  Le_RRI  = Le_RR(I)
                  Le_DSI  = Le_DS(I)
                  Le_DSSI = Le_DSS(I)           ! New, not needed though
               endif

!  Check condition

               IF (NCHECK.EQ.1 ) THEN

                  IF ( SI.GT.0.0D0 ) THEN

                     L_AR12(P) = L_AR12(P) + F1*L_B2R + F21*L_B3R
                     L_AI12(P) = L_AI12(P) + F1*L_B2I + F21*L_B3I
                     L_GR12(P) = L_GR12(P) + F1*L_C2R + F21*L_C3R
                     L_GI12(P) = L_GI12(P) + F1*L_C2I + F21*L_C3I
                     L_AR21(P) = L_AR21(P) + F1*L_B4R + F22*L_B5R
                     L_AI21(P) = L_AI21(P) + F1*L_B4I + F22*L_B5I
                     L_GR21(P) = L_GR21(P) + F1*L_C4R + F22*L_C5R
                     L_GI21(P) = L_GI21(P) + F1*L_C4I + F22*L_C5I

                     if ( P .eq. 3 ) then
                        Le_F1 = Le_RRI*AA2 + RRI*Le_AA2
                        Le_F21 = AN1 * ( Le_RRI*URI*A12 + RRI*Le_URI*A12 + RRI*URI*Le_A12 )
                        L_AR12(P) = L_AR12(P) + Le_F1*B2R + Le_F21*B3R
                        L_AI12(P) = L_AI12(P) + Le_F1*B2I + Le_F21*B3I
                        L_GR12(P) = L_GR12(P) + Le_F1*C2R + Le_F21*C3R
                        L_GI12(P) = L_GI12(P) + Le_F1*C2I + Le_F21*C3I
                        Le_F22 = AN2 * ( Le_RRI*URI*A21 + RRI*Le_URI*A21 + RRI*URI*Le_A21 )
                        L_AR21(P) = L_AR21(P) + Le_F1*B4R + Le_F22*B5R
                        L_AI21(P) = L_AI21(P) + Le_F1*B4I + Le_F22*B5I
                        L_GR21(P) = L_GR21(P) + Le_F1*C4R + Le_F22*C5R
                        L_GI21(P) = L_GI21(P) + Le_F1*C4I + Le_F22*C5I
                     ENDIF

                  ELSE

                     L_AR11(P) = L_AR11(P) + E1*L_B1R
                     L_AI11(P) = L_AI11(P) + E1*L_B1I
                     L_GR11(P) = L_GR11(P) + E1*L_C1R
                     L_GI11(P) = L_GI11(P) + E1*L_C1I
                     L_AR22(P) = L_AR22(P) + E1*L_B6R + E2*L_B7R + E3*L_B8R
                     L_AI22(P) = L_AI22(P) + E1*L_B6I + E2*L_B7I + E3*L_B8I
                     L_GR22(P) = L_GR22(P) + E1*L_C6R + E2*L_C7R + E3*L_C8R
                     L_GI22(P) = L_GI22(P) + E1*L_C6I + E2*L_C7I + E3*L_C8I

                     if ( P .eq. 3 ) then
                        Le_E1 = Le_DSI*AA1 + DSI*Le_AA1
                        L_AR11(P) = L_AR11(P) + Le_E1*B1R
                        L_AI11(P) = L_AI11(P) + Le_E1*B1I
                        L_GR11(P) = L_GR11(P) + Le_E1*C1R
                        L_GI11(P) = L_GI11(P) + Le_E1*C1I
                        Le_EH = Le_DSI*URI*A11 + DSI*Le_URI*A11 + DSI*URI*Le_A11
                        Le_E3 = Le_EH * AN2
                        Le_E2 = Le_EH * AN1
                        L_AR22(P) = L_AR22(P) + Le_E1*B6R + Le_E2*B7R + Le_E3*B8R
                        L_AI22(P) = L_AI22(P) + Le_E1*B6I + Le_E2*B7I + Le_E3*B8I
                        L_GR22(P) = L_GR22(P) + Le_E1*C6R + Le_E2*C7R + Le_E3*C8R
                        L_GI22(P) = L_GI22(P) + Le_E1*C6I + Le_E2*C7I + Le_E3*C8I
                     ENDIF

!  End check condition

                  ENDIF
               ENDIF

!  End linear and 200 loops

            ENDDO
         ENDDO

!  Form the R/I and RG/IG matrices, and their linearizations

         AN12=ANN(N1,N2)*FACTOR

         R11(N1,N2)   = AR11*AN12  ;   R12(N1,N2)  = AR12*AN12
         R21(N1,N2)   = AR21*AN12  ;   R22(N1,N2)  = AR22*AN12
         I11(N1,N2)   = AI11*AN12  ;   I12(N1,N2)  = AI12*AN12
         I21(N1,N2)   = AI21*AN12  ;   I22(N1,N2)  = AI22*AN12
         RG11(N1,N2)  = GR11*AN12  ;   RG12(N1,N2) = GR12*AN12
         RG21(N1,N2)  = GR21*AN12  ;   RG22(N1,N2) = GR22*AN12
         IG11(N1,N2)  = GI11*AN12  ;   IG12(N1,N2) = GI12*AN12
         IG21(N1,N2)  = GI21*AN12  ;   IG22(N1,N2) = GI22*AN12

         DO P = 1, NLIN
            L_R11(P,N1,N2)  = L_AR11(P)*AN12   ;   L_R12(P,N1,N2)  = L_AR12(P)*AN12
            L_R21(P,N1,N2)  = L_AR21(P)*AN12   ;   L_R22(P,N1,N2)  = L_AR22(P)*AN12
            L_I11(P,N1,N2)  = L_AI11(P)*AN12   ;   L_I12(P,N1,N2)  = L_AI12(P)*AN12
            L_I21(P,N1,N2)  = L_AI21(P)*AN12   ;   L_I22(P,N1,N2)  = L_AI22(P)*AN12
            L_RG11(P,N1,N2) = L_GR11(P)*AN12   ;   L_RG12(P,N1,N2) = L_GR12(P)*AN12
            L_RG21(P,N1,N2) = L_GR21(P)*AN12   ;   L_RG22(P,N1,N2) = L_GR22(P)*AN12
            L_IG11(P,N1,N2) = L_GI11(P)*AN12   ;   L_IG12(P,N1,N2) = L_GI12(P)*AN12
            L_IG21(P,N1,N2) = L_GI21(P)*AN12   ;   L_IG22(P,N1,N2) = L_GI22(P)*AN12
         ENDDO

!  End N1 and N2 loops

      ENDDO
   ENDDO 

!  Develop RGQ and Q matrices
!  --------------------------

!  Extended precision ??? (Why the T prefix in this section)

!  constants

   TPIR = PIR 
   TPII = PII
   TPPI = PPI   ! has no linearizations !!!!!!!!!!!!
   L_TPIR(1) = TPIR   ; L_TPIR(2) = 0.0d0  ; L_TPIR(3) = 0.0d0 
   L_TPII(1) = 0.0d0  ; L_TPII(2) = TPII   ; L_TPII(3) = 0.0d0

!  Loops

   NM=NMAX-MM1+1
   DO N1=MM1,NMAX
      K1=N1-MM1+1      ;  KK1=K1+NM
      DO N2=MM1,NMAX
         K2=N2-MM1+1   ;  KK2=K2+NM
 
!  Regular code

         TAR11 = -R11(N1,N2)    ;  TAI11 = -I11(N1,N2)
         TGR11 = -RG11(N1,N2)   ;  TGI11 = -IG11(N1,N2)
         TAR12 =  I12(N1,N2)    ;  TAI12 = -R12(N1,N2)
         TGR12 =  IG12(N1,N2)   ;  TGI12 = -RG12(N1,N2)

         TAR21 = -I21(N1,N2)    ;  TAI21 =  R21(N1,N2)
         TGR21 = -IG21(N1,N2)   ;  TGI21 =  RG21(N1,N2)
         TAR22 = -R22(N1,N2)    ;  TAI22 = -I22(N1,N2)
         TGR22 = -RG22(N1,N2)   ;  TGI22 = -IG22(N1,N2)

         TQR(K1,K2)     = TPIR*TAR21-TPII*TAI21+TPPI*TAR12
         TQI(K1,K2)     = TPIR*TAI21+TPII*TAR21+TPPI*TAI12
         TRGQR(K1,K2)   = TPIR*TGR21-TPII*TGI21+TPPI*TGR12
         TRGQI(K1,K2)   = TPIR*TGI21+TPII*TGR21+TPPI*TGI12

         TQR(K1,KK2)    = TPIR*TAR11-TPII*TAI11+TPPI*TAR22
         TQI(K1,KK2)    = TPIR*TAI11+TPII*TAR11+TPPI*TAI22
         TRGQR(K1,KK2)  = TPIR*TGR11-TPII*TGI11+TPPI*TGR22
         TRGQI(K1,KK2)  = TPIR*TGI11+TPII*TGR11+TPPI*TGI22

         TQR(KK1,K2)    = TPIR*TAR22-TPII*TAI22+TPPI*TAR11
         TQI(KK1,K2)    = TPIR*TAI22+TPII*TAR22+TPPI*TAI11
         TRGQR(KK1,K2)  = TPIR*TGR22-TPII*TGI22+TPPI*TGR11
         TRGQI(KK1,K2)  = TPIR*TGI22+TPII*TGR22+TPPI*TGI11

         TQR(KK1,KK2)   = TPIR*TAR12-TPII*TAI12+TPPI*TAR21
         TQI(KK1,KK2)   = TPIR*TAI12+TPII*TAR12+TPPI*TAI21
         TRGQR(KK1,KK2) = TPIR*TGR12-TPII*TGI12+TPPI*TGR21
         TRGQI(KK1,KK2) = TPIR*TGI12+TPII*TGR12+TPPI*TGI21

!  Linearized code

         DO P = 1, NLIN

            L_TAR11 = -L_R11(P,N1,N2)    ;  L_TAI11 = -L_I11(P,N1,N2)
            L_TGR11 = -L_RG11(P,N1,N2)   ;  L_TGI11 = -L_IG11(P,N1,N2)
            L_TAR12 =  L_I12(P,N1,N2)    ;  L_TAI12 = -L_R12(P,N1,N2)
            L_TGR12 =  L_IG12(P,N1,N2)   ;  L_TGI12 = -L_RG12(P,N1,N2)

            L_TAR21 = -L_I21(P,N1,N2)    ;  L_TAI21 =  L_R21(P,N1,N2)
            L_TGR21 = -L_IG21(P,N1,N2)   ;  L_TGI21 =  L_RG21(P,N1,N2)
            L_TAR22 = -L_R22(P,N1,N2)    ;  L_TAI22 = -L_I22(P,N1,N2)
            L_TGR22 = -L_RG22(P,N1,N2)   ;  L_TGI22 = -L_IG22(P,N1,N2)

            L_TQR(P,K1,K2)     = L_TPIR(P)*TAR21 - L_TPII(P)*TAI21 + TPPI*L_TAR12  &
                                 + TPIR*L_TAR21  -   TPII*L_TAI21
            L_TQI(P,K1,K2)     = L_TPIR(P)*TAI21 + L_TPII(P)*TAR21 + TPPI*L_TAI12  &
                                 + TPIR*L_TAI21  +   TPII*L_TAR21
            L_TRGQR(P,K1,K2)   = L_TPIR(P)*TGR21 - L_TPII(P)*TGI21 + TPPI*L_TGR12  &
                                 + TPIR*L_TGR21  -   TPII*L_TGI21
            L_TRGQI(P,K1,K2)   = L_TPIR(P)*TGI21 + L_TPII(P)*TGR21 + TPPI*L_TGI12  &
                                 + TPIR*L_TGI21  +   TPII*L_TGR21

            L_TQR(P,K1,KK2)    = L_TPIR(P)*TAR11 - L_TPII(P)*TAI11 + TPPI*L_TAR22  &
                                 + TPIR*L_TAR11  -   TPII*L_TAI11
            L_TQI(P,K1,KK2)    = L_TPIR(P)*TAI11 + L_TPII(P)*TAR11 + TPPI*L_TAI22  &
                                 + TPIR*L_TAI11  +   TPII*L_TAR11
            L_TRGQR(P,K1,KK2)  = L_TPIR(P)*TGR11 - L_TPII(P)*TGI11 + TPPI*L_TGR22  &
                                 + TPIR*L_TGR11  -   TPII*L_TGI11
            L_TRGQI(P,K1,KK2)  = L_TPIR(P)*TGI11 + L_TPII(P)*TGR11 + TPPI*L_TGI22  &
                                 + TPIR*L_TGI11  +   TPII*L_TGR11

            L_TQR(P,KK1,K2)    = L_TPIR(P)*TAR22 - L_TPII(P)*TAI22 + TPPI*L_TAR11  &
                                 + TPIR*L_TAR22  -   TPII*L_TAI22
            L_TQI(P,KK1,K2)    = L_TPIR(P)*TAI22 + L_TPII(P)*TAR22 + TPPI*L_TAI11  &
                                 + TPIR*L_TAI22  +   TPII*L_TAR22
            L_TRGQR(P,KK1,K2)  = L_TPIR(P)*TGR22 - L_TPII(P)*TGI22 + TPPI*L_TGR11  &
                                 + TPIR*L_TGR22  -   TPII*L_TGI22
            L_TRGQI(P,KK1,K2)  = L_TPIR(P)*TGI22 + L_TPII(P)*TGR22 + TPPI*L_TGI11  &
                                 + TPIR*L_TGI22  +   TPII*L_TGR22

            L_TQR(P,KK1,KK2)   = L_TPIR(P)*TAR12 - L_TPII(P)*TAI12 + TPPI*L_TAR21  &
                                 + TPIR*L_TAR12  -   TPII*L_TAI12
            L_TQI(P,KK1,KK2)   = L_TPIR(P)*TAI12 + L_TPII(P)*TAR12 + TPPI*L_TAI21  &
                                 + TPIR*L_TAI12  +   TPII*L_TAR12
            L_TRGQR(P,KK1,KK2) = L_TPIR(P)*TGR12 - L_TPII(P)*TGI12 + TPPI*L_TGR21  &
                                 + TPIR*L_TGR12  -   TPII*L_TGI12
            L_TRGQI(P,KK1,KK2) = L_TPIR(P)*TGI12 + L_TPII(P)*TGR12 + TPPI*L_TGI21  &
                                 + TPIR*L_TGI12  +   TPII*L_TGR12

         ENDDO

!  End main loops

      ENDDO
   ENDDO

!  Copy to Normal RGQ and Q arrays (Why ????)

   NNMAX = 2 * NM
   DO N1 = 1,NNMAX
      DO N2=1,NNMAX
         QR(N1,N2)=TQR(N1,N2)
         QI(N1,N2)=TQI(N1,N2)
         RGQR(N1,N2)=TRGQR(N1,N2)
         RGQI(N1,N2)=TRGQI(N1,N2)
         DO P = 1, NLIN
            L_QR(P,N1,N2)=L_TQR(P,N1,N2)
            L_QI(P,N1,N2)=L_TQI(P,N1,N2)
            L_RGQR(P,N1,N2)=L_TRGQR(P,N1,N2)
            L_RGQI(P,N1,N2)=L_TRGQI(P,N1,N2)
         ENDDO
!         P = 1
!         write(57,'(2i4,1p2e20.9)')n1,n2,QR(N1,N2),QI(N1,N2)
!         write(57,'(2i4,1p4e20.10)')n1,n2,L_QR(P,N1,N2),L_QI(P,N1,N2),&
!                                          QR(N1,N2),QI(N1,N2)
!         write(57,'(2i4,1p4e20.10)')n1,n2,L_RGQR(P,N1,N2),L_RGQI(P,N1,N2),&
!                                          RGQR(N1,N2),RGQI(N1,N2)
      ENDDO
   ENDDO

!  Call T-matrix calculator (uses LAPACK inversion)

   call RGQ_Qinv_plus                               &
          ( NM, NLIN, QR, QI, RGQR, RGQI,           & ! Inputs
            L_QR, L_QI, L_RGQR, L_RGQI,             & ! Inputs
            TR1, TI1, L_TR1, L_TI1, fail, message )   ! Outputs

!  debug

!        P = 1
!       DO N1=1,NNMAX
!           DO N2=1,NNMAX
!         write(57,'(2i4,1p4e20.10)')n1,n2,L_TR1(P,N1,N2),L_TI1(P,N1,N2),&
!                                          TR1(N1,N2),TI1(N1,N2)
!         write(57,'(2i4,1p2e20.9)')n1,n2,TR1(N1,N2),TI1(N1,N2)
!       enddo
!       enddo
!      pause'57_spadgett'

!  Exception handling

   if ( fail ) then
      TRACE = ' RGQ_Qinv_plus failed in Tmatrix_R_plus'
   endif

!  Finish

   return
end subroutine tmatrix_R_plus

subroutine VIG_PLUS                  &
      ( NMAX, X, L_X, M,             & ! Inputs
        DV1, DV2, L_DV1, L_DV2 )       ! Outputs

!  Polynomial expansions

   implicit none

!  input arguments
!  ---------------

   integer  , intent(in)  :: nmax, m
   REAL(fpk), intent(in)  :: x, L_x

!  Output arguments
!  ----------------

   REAL(fpk), intent(out) :: DV1(NPN1),DV2(NPN1)
   REAL(fpk), intent(out) :: L_DV1(NPN1),L_DV2(NPN1)

!  Local variables
!  ---------------

   INTEGER   :: N, I, I2
   REAL(fpk) :: A, D1, D2, D3, DER, FAC
   REAL(fpk) :: L_A, L_D1, L_D2, L_D3, L_DER, L_QS, L_QS1
   REAL(fpk) :: QS, QS1, QMM, QN, QN1, QN2, QNM, QNM1

!  initialize

   A   = 1.0D0
   QS  = DSQRT(1.0D0-X*X)
   QS1 = 1.0D0/QS

   L_QS  = - X * L_X / QS
   L_QS1 = - QS1 * QS1 * L_QS
   L_A = 0.0d0

   DO N = 1, NMAX
      DV1(N) = 0.0D0
      DV2(N) = 0.0D0
      L_DV1(N) = 0.0D0
      L_DV2(N) = 0.0D0
   ENDDO   

!  Two cases, M = 0, and M =/ 0

   IF ( M.NE.0) then
      QMM = DBLE(M*M)
      DO I = 1,M
         I2 = I*2
         FAC = DSQRT(DBLE(I2-1)/DBLE(I2))
         L_A = ( L_A * QS + A * L_QS ) * FAC
         A   = A * FAC * QS
      ENDDO   

      D1 = 0.0D0
      D2 = A
      L_D1 = 0.0d0
      L_D2 = L_A
 
      DO N = M, NMAX

         QN=DBLE(N)
         QN2=DBLE(2*N+1)
         QN1=DBLE(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)

         D3   = (QN2*X*D2-QNM*D1)/QNM1
         L_D3 = (QN2*(L_X*D2+X*L_D2)-QNM*L_D1)/QNM1        

         DER   = QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2
         L_DER = ( QS1*(-QN1*QNM*L_D1+QN*QNM1*L_D3)+ L_QS1*(-QN1*QNM*D1+QN*QNM1*D3) ) /QN2

         DV1(N)=D2
         DV2(N)=DER
         L_DV1(N)=L_D2
         L_DV2(N)=L_DER
         D1=D2
         D2=D3
         L_D1=L_D2
         L_D2=L_D3

      ENDDO   
   else
      D1 = 1.0D0
      D2 = X  
      L_D1 = 0.0d0
      L_D2 = L_X
      DO N = 1, NMAX  

         QN =DBLE(N)
         QN1=DBLE(N+1)
         QN2=DBLE(2*N+1)
         FAC = (QN1*QN/QN2)

         D3   = (QN2*X*D2-QN*D1)/QN1 
         L_D3 = (QN2*(L_X*D2+X*L_D2)-QN*L_D1)/QN1        

         DER   = FAC * QS1 * (-D1+D3)
         L_DER = FAC * ( QS1 * (-L_D1+L_D3) + L_QS1 * (-D1+D3) )

         DV1(N)=D2
         DV2(N)=DER
         L_DV1(N)=L_D2
         L_DV2(N)=L_DER
         D1=D2
         D2=D3
         L_D1=L_D2
         L_D2=L_D3

      ENDDO   
   endif

!  Finish

   return
end subroutine VIG_PLUS
 
subroutine RGQ_Qinv_plus                            &
          ( NMAX, NLIN, QR, QI, RGQR, RGQI,         & ! Inputs
            L_QR, L_QI, L_RGQR, L_RGQI,             & ! Inputs
            TR1, TI1, L_TR1, L_TI1, fail, message )   ! Outputs

!             T-matrix     TMAT = - RG(Q) * (Q**(-1))

!  Linearized T-matrix  L(TMAT) = - RGQ * L(QINV) - L(RGQ) * QINV
!          where        L(QINV) = - QINV * L(Q) * QINV

!     ===> L(TMAT)    = + RGQ * QINV * L(Q) * QINV - L(RGQ) * QINV
!     ===> L(TMAT)    = - ( TMAT * L(Q) + L(RGQ) ) * QINV

!  Thus, the Matrix operations are

!     ===>   TMAT     = - RGQ  * QINV
!     ===> L(TMAT)    = - WLQ  * QINV
!          where WLQ  =  TMAT * L(Q) + L(RGQ)

   implicit none

!  input arguments
!  ---------------

!  Dimensions

   integer  , intent(in)  :: nmax, nlin

!  RGQ and Q matrices

   REAL(fpk), intent(in)  ::   QR(NPN2,NPN2),  QI(NPN2,NPN2)
   REAL(fpk), intent(in)  :: RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)

!  Linearized RGQ and Q matrices

   REAL(fpk), intent(in)  ::   L_QR(3,NPN2,NPN2),  L_QI(3,NPN2,NPN2)
   REAL(fpk), intent(in)  ::   L_RGQR(3,NPN2,NPN2),L_RGQI(3,NPN2,NPN2)

!  Output arguments
!  ----------------

!  T-matrix output

   REAL(fpk), intent(out) :: TR1(NPN2,NPN2), TI1(NPN2,NPN2)
   REAL(fpk), intent(out) :: L_TR1(3,NPN2,NPN2), L_TI1(3,NPN2,NPN2)

!  Exception handling

   logical       , intent(out) :: fail
   character*(*) , intent(out) :: message

!  Local variables
!  ---------------

!  Local arrays for the LAPACK calculation and final multiplication
!   QINV is stored in ZQ aftern the LAPACK operations

   COMPLEX(fpk) :: ZQ(NPN2,NPN2),ZW(NPN2)
   INTEGER         :: IPIV(NPN2)

!  Auxiliary Matrices

   REAL(fpk)    :: WLQR(NPN2,NPN2), WLQI(NPN2,NPN2)

!  Other variables

   CHARACTER*2  :: C2
   INTEGER      :: I, J, K, M, INFO, NNMAX
   REAL(fpk)    :: AR, AI, ARR, ARI, TR, TI

!  Initialize
!  ----------

!  Exception handling

   FAIL    = .false.
   MESSAGE = ' '

!  actual dimensions

   NNMAX = 2*NMAX
 
!  Develop complex matrix

!mick fix 6/19/2014 - switched CMPLX to f90/f95 form in the two locations below
   !ZQ =   DCMPLX(0.0d0,0.0d0)
   ZQ =   CMPLX(0.0d0,0.0d0,fpk)
   ZW =   0.0d0
   IPIV = 0

   DO I=1,NNMAX
      DO J=1,NNMAX
         !ZQ(I,J) = DCMPLX(QR(I,J),QI(I,J))
         ZQ(I,J) = CMPLX(QR(I,J),QI(I,J),fpk)
      ENDDO
   ENDDO

!  Matrix inversion from LAPACK
!  ----------------------------

!  LU decomposition, with exception handling

   INFO = 0
   CALL ZGETRF(NNMAX,NNMAX,ZQ,NPN2,IPIV,INFO)
   IF (INFO.NE.0) THEN
      write(C2,'(I2)')INFO
      FAIL    = .true.
      MESSAGE = ' ZGETRF call: Degenerative crapumster, Info =  '//C2
      RETURN
   ENDIF

!  Back substitution, with exception handling

   CALL ZGETRI(NNMAX,ZQ,NPN2,IPIV,ZW,NPN2,INFO)
   IF (INFO.NE.0) THEN
      write(C2,'(I2)')INFO
      FAIL    = .true.
      MESSAGE = ' ZGETRI call: Degenerative crapumster, Info =  '//C2
      RETURN
   ENDIF

!  T-matrix completion (Multiplication)

   DO I=1,NNMAX
      DO J=1,NNMAX
         TR=0.0D0
         TI=0.0D0
         DO K=1,NNMAX
            ARR = RGQR(I,K)
            ARI = RGQI(I,K)
            AR  = ZQ(K,J)
!mick fix 6/19/2014 - switched DIMAG to f90/f95 form in the two locations below
            !AI  = DIMAG(ZQ(K,J))
            AI  = AIMAG(ZQ(K,J))
            TR  = TR-ARR*AR+ARI*AI
            TI  = TI-ARR*AI-ARI*AR
         ENDDO
         TR1(I,J)=TR
         TI1(I,J)=TI
      ENDDO
   ENDDO

!  Linearizations

   DO M = 1, NLIN

!  Develop WLQ (Multiplication)

      DO I=1,NNMAX
         DO J=1,NNMAX
            TR=0.0D0
            TI=0.0D0
            DO K=1,NNMAX
               ARR = TR1(I,K)
               ARI = TI1(I,K)
               AR  = L_QR(M,K,J)
               AI  = L_QI(M,K,J)
               TR  = TR + ARR*AR - ARI*AI
               TI  = TI + ARR*AI + ARI*AR
            ENDDO
            WLQR(I,J) = TR + L_RGQR(M,I,J)
            WLQI(I,J) = TI + L_RGQI(M,I,J)
         ENDDO
      ENDDO

!  Linearized T-matrix completion (Multiplication)

      DO I=1,NNMAX
         DO J=1,NNMAX
            TR=0.0D0
            TI=0.0D0
            DO K=1,NNMAX
               ARR = WLQR(I,K)
               ARI = WLQI(I,K)
               AR  = ZQ(K,J)
               !AI  = DIMAG(ZQ(K,J))
               AI  = AIMAG(ZQ(K,J))
               TR  = TR-ARR*AR+ARI*AI
               TI  = TI-ARR*AI-ARI*AR
            ENDDO
            L_TR1(M,I,J)=TR
            L_TI1(M,I,J)=TI
         ENDDO
      ENDDO

!  End parameter loop

   ENDDO

!  Finish

   return
end subroutine RGQ_Qinv_plus

!  End module

end module tmat_makers_plus
