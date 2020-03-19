 module tmat_makers

   use tmat_parameters, only : fpk => tmat_fpkind, &
                                NPN1, NPN2, NPNG2
   use lapack_tools

!  LIST OF SUBROUTINES
!  ===================

!    tmatrix_R 
!    tmatrix_R0

!       RGQ_Qinv
!       VIG

private  RGQ_Qinv
public   tmatrix_R, tmatrix_R0, VIG

contains

subroutine tmatrix_R0                              &
     ( NGAUSS, NMAX, NCHECK, X, W, AN, ANN,        & ! Inputs
       PPI, PIR, PII, R, DR, DDR, DRR, DRI,        & ! Inputs
        j_bess,  y_bess,  jr_bess,  ji_bess,       & ! Inputs
       dj_bess, dy_bess, djr_bess, dji_bess,       & ! Inputs
       R12, R21, I12, I21, RG12, RG21, IG12, IG21, & ! Outputs
       TR1, TI1, fail, message, trace )              ! Outputs

   implicit none

!  input arguments
!  ---------------

   integer  , intent(in)  :: nmax, ncheck, ngauss

   REAL(fpk), intent(in)  :: X(NPNG2), W(NPNG2)
   REAL(fpk), intent(in)  :: AN(NPN1),ANN(NPN1,NPN1)

   real(fpk), intent(in)  :: PPI, PIR, PII

   REAL(fpk), intent(in)  :: R(NPNG2), DR(NPNG2)
   real(fpk), intent(in)  :: DDR(NPNG2), DRR(NPNG2), DRI(NPNG2)

!  Bessel Master output

   real(fpk), intent(in)  :: J_BESS  (NPNG2,NPN1)
   real(fpk), intent(in)  :: Y_BESS  (NPNG2,NPN1)
   real(fpk), intent(in)  :: JR_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: JI_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DJ_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DY_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DJR_BESS(NPNG2,NPN1)
   real(fpk), intent(in)  :: DJI_BESS(NPNG2,NPN1)

!  Output arguments
!  ----------------

!  Miscellaneous output (Tmatrix)

   REAL(fpk), intent(out) :: R12 (NPN1,NPN1), R21 (NPN1,NPN1)
   REAL(fpk), intent(out) :: I12 (NPN1,NPN1), I21 (NPN1,NPN1)
   REAL(fpk), intent(out) :: RG12(NPN1,NPN1), RG21(NPN1,NPN1)
   REAL(fpk), intent(out) :: IG12(NPN1,NPN1), IG21(NPN1,NPN1)

!  T-matrix output

   REAL(fpk), intent(out) :: TR1(NPN2,NPN2), TI1(NPN2,NPN2)

!  Exception handling

   logical       , intent(out) :: fail
   character*(*) , intent(out) :: message
   character*(*) , intent(out) :: trace

!  Local variables
!  ---------------

!  Local arrays

   REAL(fpk)  :: D1(NPNG2,NPN1), D2(NPNG2,NPN1)
   REAL(fpk)  :: RR(NPNG2), SIG(NPN2)
   REAL(fpk)  :: DV1(NPN1), DV2(NPN1)

!  RGQ and Q matrices

   REAL(fpk)  ::   QR(NPN2,NPN2),  QI(NPN2,NPN2)
   REAL(fpk)  :: RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)

   REAL(fpk)  ::   TQR(NPN2,NPN2),  TQI(NPN2,NPN2)
   REAL(fpk)  :: TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)

!  Other variables

   LOGICAL   :: DO_200
   INTEGER   :: MM1, NNMAX, NG, NGSS, NM
   INTEGER   :: I, I1, I2, N, N1, N2, K1, KK1, K2, KK2
   REAL(fpk) :: FACTOR, SI, DD1, DD2, DDRI, DRRI, DRII
   REAL(fpk) :: AN1, AN2, AN12, A12, A21, A22, AA1
   REAL(fpk) :: URI, RRI, F1, F2, TPIR, TPII, TPPI

   REAL(fpk) :: AR12, AR21, AI12, AI21
   REAL(fpk) :: GR12, GR21, GI12, GI21

   REAL(fpk) :: TAR12, TAR21, TAI12, TAI21
   REAL(fpk) :: TGR12, TGR21, TGI12, TGI21

   REAL(fpk) :: D1N1, D2N1, D1N2, D2N2
   REAL(fpk) :: QJ1, QY1, QJR2, QJI2
   REAL(fpk) :: QDJR2, QDJI2, QDJ1, QDY1

   REAL(fpk) :: B1R, B1I, B2R, B2I, B3R, B3I, B4R, B4I, B5R, B5I
   REAL(fpk) :: C1R, C1I, C2R, C2I, C3R, C3I, C4R, C4I, C5R, C5I

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

!  VIG call
!  --------

   DO I=1,NGAUSS
      I1=NGAUSS+I
      I2=NGAUSS-I+1
      CALL VIG ( NMAX, X(I1), 0, DV1, DV2)
      DO  N=1,NMAX
         SI=SIG(N)
         DD1=DV1(N)
         DD2=DV2(N)
         D1(I1,N)=DD1
         D2(I1,N)=DD2
         D1(I2,N)=DD1*SI
         D2(I2,N)=-DD2*SI
     ENDDO
   ENDDO

!  Main calculation loop

   DO N1=MM1,NMAX
      AN1=AN(N1)
      DO N2=MM1,NMAX
         AN2=AN(N2)
         AR12=0.0D0
         AR21=0.0D0
         AI12=0.0D0
         AI21=0.0D0
         GR12=0.0D0
         GR21=0.0D0
         GI12=0.0D0
         GI21=0.0D0

!  Avoidance condition

         DO_200 = .true.
         IF ( NCHECK.EQ.1.AND.SIG(N1+N2).LT.0.0D0) DO_200 = .false.

!  Do this loop (200)

         IF ( DO_200 ) THEN
            DO I=1,NGSS
               D1N1=D1(I,N1)
               D2N1=D2(I,N1)
               D1N2=D1(I,N2)
               D2N2=D2(I,N2)
               A12=D1N1*D2N2
               A21=D2N1*D1N2
               A22=D2N1*D2N2
               AA1=A12+A21
 
               QJ1=J_BESS(I,N1)
               QY1=Y_BESS(I,N1)
               QJR2=JR_BESS(I,N2)
               QJI2=JI_BESS(I,N2)
               QDJR2=DJR_BESS(I,N2)
               QDJI2=DJI_BESS(I,N2)
               QDJ1=DJ_BESS(I,N1)
               QDY1=DY_BESS(I,N1)

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

               URI=DR(I)
               RRI=RR(I)

               F1=RRI*A22
               F2=RRI*URI*AN1*A12
               AR12=AR12+F1*B2R+F2*B3R
               AI12=AI12+F1*B2I+F2*B3I
               GR12=GR12+F1*C2R+F2*C3R
               GI12=GI12+F1*C2I+F2*C3I

               F2=RRI*URI*AN2*A21
               AR21=AR21+F1*B4R+F2*B5R
               AI21=AI21+F1*B4I+F2*B5I
               GR21=GR21+F1*C4R+F2*C5R
               GI21=GI21+F1*C4I+F2*C5I
            ENDDO
         ENDIF

!  Resume

         AN12=ANN(N1,N2)*FACTOR
         R12(N1,N2)=AR12*AN12
         R21(N1,N2)=AR21*AN12
         I12(N1,N2)=AI12*AN12
         I21(N1,N2)=AI21*AN12
         RG12(N1,N2)=GR12*AN12
         RG21(N1,N2)=GR21*AN12
         IG12(N1,N2)=GI12*AN12
         IG21(N1,N2)=GI21*AN12

!         write(56,'(2i4,1p4e20.10)')N1,N2,R12(N1,N2),R21(N1,N2),&
!                                          I12(N1,N2),I21(N1,N2)
!         write(56,'(2i4,1p4e20.10)')N1,N2,RG12(N1,N2),RG21(N1,N2),&
!                                          IG12(N1,N2),IG21(N1,N2)

      ENDDO
   ENDDO

!   pause'develop'

!  Develop RGQ and Q matrices
!  --------------------------

!  Extended precision ??? (Why the T prefix in this section)

   TPIR=PIR
   TPII=PII
   TPPI=PPI
 
   NM=NMAX
   DO N1=MM1,NMAX
      K1=N1-MM1+1
      KK1=K1+NM
      DO N2=MM1,NMAX
         K2=N2-MM1+1
         KK2=K2+NM
 
         TAR12= I12(N1,N2)
         TAI12=-R12(N1,N2)
         TGR12= IG12(N1,N2)
         TGI12=-RG12(N1,N2)

         TAR21=-I21(N1,N2)
         TAI21= R21(N1,N2)
         TGR21=-IG21(N1,N2)
         TGI21= RG21(N1,N2)

         TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
         TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
         TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
         TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12

         TQR(K1,KK2)=0.0D0
         TQI(K1,KK2)=0.0D0
         TRGQR(K1,KK2)=0.0D0
         TRGQI(K1,KK2)=0.0D0

         TQR(KK1,K2)=0.0D0
         TQI(KK1,K2)=0.0D0
         TRGQR(KK1,K2)=0.0D0
         TRGQI(KK1,K2)=0.0D0

         TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
         TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
         TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
         TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
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
!         write(56,'(2i4,1p2e20.10)')n1,n2,QR(N1,N2),QI(N1,N2)
!         write(56,'(2i4,1p2e20.10)')n1,n2,RGQR(N1,N2),RGQI(N1,N2)
      ENDDO
   ENDDO

!  Call T-matrix calculator (uses LAPACK inversion)

   call RGQ_Qinv                          &
     ( NMAX, QR, QI, RGQR, RGQI,          & ! Inputs
       TR1, TI1, fail, message )            ! Outputs

!  debug

!       DO N1=1,NNMAX
!           DO N2=1,NNMAX
!         write(56,'(2i4,1p2e20.10)')n1,n2,TR1(N1,N2),TI1(N1,N2)
!       enddo
!       enddo
!      pause'56'

!  Exception handling

   if ( fail ) then
      TRACE = ' RGQ_Qinv failed in Tmatrix_R0'
   endif

!  Finish

   return
end subroutine tmatrix_R0

subroutine tmatrix_R                                  &
     ( M, NGAUSS, NMAX, NCHECK, X, W, AN, ANN, S, SS, & ! Inputs
       PPI, PIR, PII, R, DR, DDR, DRR, DRI,           & ! Inputs
        j_bess,  y_bess,  jr_bess,  ji_bess,          & ! Inputs
       dj_bess, dy_bess, djr_bess, dji_bess,          & ! Inputs
       R11, R12, R21, R22, I11, I12, I21, I22,        & ! Outputs
       RG11,RG12,RG21,RG22,IG11,IG12,IG21,IG22,       & ! Outputs
       TR1, TI1, fail, message, trace )                 ! Outputs

   implicit none

!  input arguments
!  ---------------

!  Control

   integer  , intent(in)  :: nmax, ncheck, ngauss, m

!  X and W and S/SS

   REAL(fpk), intent(in)  :: X(NPNG2), W(NPNG2)
   REAL(fpk), intent(in)  :: S(NPNG2),SS(NPNG2)

!  Constants

   REAL(fpk), intent(in)  :: AN(NPN1),ANN(NPN1,NPN1)
   real(fpk), intent(in)  :: PPI, PIR, PII

!  R and DR functions

   REAL(fpk), intent(in)  :: R(NPNG2), DR(NPNG2)
   real(fpk), intent(in)  :: DDR(NPNG2), DRR(NPNG2), DRI(NPNG2)

!  Bessel Master output

   real(fpk), intent(in)  :: J_BESS  (NPNG2,NPN1)
   real(fpk), intent(in)  :: Y_BESS  (NPNG2,NPN1)
   real(fpk), intent(in)  :: JR_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: JI_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DJ_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DY_BESS (NPNG2,NPN1)
   real(fpk), intent(in)  :: DJR_BESS(NPNG2,NPN1)
   real(fpk), intent(in)  :: DJI_BESS(NPNG2,NPN1)

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

!  T-matrix output

   REAL(fpk), intent(out) :: TR1(NPN2,NPN2), TI1(NPN2,NPN2)

!  Exception handling

   logical       , intent(out) :: fail
   character*(*) , intent(out) :: message
   character*(*) , intent(out) :: trace

!  Local variables
!  ---------------

!  Local arrays

   REAL(fpk)  :: D1(NPNG2,NPN1), D2(NPNG2,NPN1)
   REAL(fpk)  :: DS(NPNG2), DSS(NPNG2)
   REAL(fpk)  :: RR(NPNG2), SIG(NPN2)
   REAL(fpk)  :: DV1(NPN1), DV2(NPN1)

!  RGQ and Q matrices

   REAL(fpk)  ::   QR(NPN2,NPN2),  QI(NPN2,NPN2)
   REAL(fpk)  :: RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)

   REAL(fpk)  ::   TQR(NPN2,NPN2),  TQI(NPN2,NPN2)
   REAL(fpk)  :: TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)

!  Other variables

   INTEGER   :: MM1, NNMAX, NG, NGSS, NM
   INTEGER   :: I, I1, I2, N, N1, N2, K1, KK1, K2, KK2
   REAL(fpk) :: FACTOR, SI, DD1, DD2, DDRI, DRRI, DRII, QM, QMM
   REAL(fpk) :: URI, RRI, DSI, DSSI, WR
   REAL(fpk) :: AN1, AN2, AN12, A11, A12, A21, A22, AA1, AA2
   REAL(fpk) :: E1, E2, E3, F1, F2, TPIR, TPII, TPPI

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
      DS(I)  = S(I)  * QM*WR
      DSS(I) = SS(I) * QMM
      RR(I)  = W(I)  * R(I)
   ENDDO

!  VIG call
!  --------

   DO I=1,NGAUSS
      I1=NGAUSS+I
      I2=NGAUSS-I+1
      CALL VIG ( NMAX, X(I1), M, DV1, DV2)
      DO  N=1,NMAX
         SI=SIG(N)
         DD1=DV1(N)
         DD2=DV2(N)
         D1(I1,N)=DD1
         D2(I1,N)=DD2
         D1(I2,N)=DD1*SI
         D2(I2,N)=-DD2*SI
     ENDDO
   ENDDO

!  Main calculation loop

   DO N1=MM1,NMAX
      AN1=AN(N1)
      DO N2=MM1,NMAX

         AN2=AN(N2)
         SI=SIG(N1+N2)

!  Zero

         AR11=0.0D0
         AR12=0.0D0
         AR21=0.0D0
         AR22=0.0D0
         AI11=0.0D0
         AI12=0.0D0
         AI21=0.0D0
         AI22=0.0D0
         GR11=0.0D0
         GR12=0.0D0
         GR21=0.0D0
         GR22=0.0D0
         GI11=0.0D0
         GI12=0.0D0
         GI21=0.0D0
         GI22=0.0D0
 
!  200 loop

         DO I=1,NGSS
            D1N1=D1(I,N1)
            D2N1=D2(I,N1)
            D1N2=D1(I,N2)
            D2N2=D2(I,N2)
            A11=D1N1*D1N2
            A12=D1N1*D2N2
            A21=D2N1*D1N2
            A22=D2N1*D2N2
            AA1=A12+A21
            AA2=A11*DSS(I)+A22
            QJ1=J_BESS(I,N1)
            QY1=Y_BESS(I,N1)
            QJR2=JR_BESS(I,N2)
            QJI2=JI_BESS(I,N2)
            QDJR2=DJR_BESS(I,N2)
            QDJI2=DJI_BESS(I,N2)
            QDJ1=DJ_BESS(I,N1)
            QDY1=DY_BESS(I,N1)
 
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
                  F2=RRI*URI*AN1*A12
                  AR12=AR12+F1*B2R+F2*B3R
                  AI12=AI12+F1*B2I+F2*B3I
                  GR12=GR12+F1*C2R+F2*C3R
                  GI12=GI12+F1*C2I+F2*C3I
                  F2=RRI*URI*AN2*A21
                  AR21=AR21+F1*B4R+F2*B5R
                  AI21=AI21+F1*B4I+F2*B5I
                  GR21=GR21+F1*C4R+F2*C5R
                  GI21=GI21+F1*C4I+F2*C5I
               ELSE
                  E1=DSI*AA1
                  AR11=AR11+E1*B1R
                  AI11=AI11+E1*B1I
                  GR11=GR11+E1*C1R
                  GI11=GI11+E1*C1I
                  E2=DSI*URI*A11
                  E3=E2*AN2
                  E2=E2*AN1
                  AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                  AI22=AI22+E1*B6I+E2*B7I+E3*B8I
                  GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                  GI22=GI22+E1*C6I+E2*C7I+E3*C8I
               ENDIF
            ENDIF

         ENDDO

!  Resume after 200 loop

         AN12=ANN(N1,N2)*FACTOR
         R11(N1,N2)=AR11*AN12
         R12(N1,N2)=AR12*AN12
         R21(N1,N2)=AR21*AN12
         R22(N1,N2)=AR22*AN12
         I11(N1,N2)=AI11*AN12
         I12(N1,N2)=AI12*AN12
         I21(N1,N2)=AI21*AN12
         I22(N1,N2)=AI22*AN12
         RG11(N1,N2)=GR11*AN12
         RG12(N1,N2)=GR12*AN12
         RG21(N1,N2)=GR21*AN12
         RG22(N1,N2)=GR22*AN12
         IG11(N1,N2)=GI11*AN12
         IG12(N1,N2)=GI12*AN12
         IG21(N1,N2)=GI21*AN12
         IG22(N1,N2)=GI22*AN12

      ENDDO
   ENDDO 

!  Develop RGQ and Q matrices
!  --------------------------

!  Extended precision ??? (Why the T prefix in this section)

   TPIR=PIR
   TPII=PII
   TPPI=PPI

   NM=NMAX-MM1+1
   DO N1=MM1,NMAX
      K1=N1-MM1+1
      KK1=K1+NM
      DO N2=MM1,NMAX
         K2=N2-MM1+1
         KK2=K2+NM

         TAR11=-R11(N1,N2)
         TAI11=-I11(N1,N2)
         TGR11=-RG11(N1,N2)
         TGI11=-IG11(N1,N2)

         TAR12= I12(N1,N2)
         TAI12=-R12(N1,N2)
         TGR12= IG12(N1,N2)
         TGI12=-RG12(N1,N2)

         TAR21=-I21(N1,N2)
         TAI21= R21(N1,N2)
         TGR21=-IG21(N1,N2)
         TGI21= RG21(N1,N2)

         TAR22=-R22(N1,N2)
         TAI22=-I22(N1,N2)
         TGR22=-RG22(N1,N2)
         TGI22=-IG22(N1,N2)

         TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
         TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
         TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
         TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12

         TQR(K1,KK2)=TPIR*TAR11-TPII*TAI11+TPPI*TAR22
         TQI(K1,KK2)=TPIR*TAI11+TPII*TAR11+TPPI*TAI22
         TRGQR(K1,KK2)=TPIR*TGR11-TPII*TGI11+TPPI*TGR22
         TRGQI(K1,KK2)=TPIR*TGI11+TPII*TGR11+TPPI*TGI22

         TQR(KK1,K2)=TPIR*TAR22-TPII*TAI22+TPPI*TAR11
         TQI(KK1,K2)=TPIR*TAI22+TPII*TAR22+TPPI*TAI11
         TRGQR(KK1,K2)=TPIR*TGR22-TPII*TGI22+TPPI*TGR11
         TRGQI(KK1,K2)=TPIR*TGI22+TPII*TGR22+TPPI*TGI11

         TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
         TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
         TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
         TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
      ENDDO
   ENDDO

!  Copy to Normal RGQ and Q arrays (Why ????)

!   QR = 0.0d0;  QI = 0.0d0;  RGQR = 0.0d0;  RGQI = 0.0d0
   NNMAX = 2 * NM
   DO N1 = 1,NNMAX
      DO N2=1,NNMAX
         QR(N1,N2)=TQR(N1,N2)
         QI(N1,N2)=TQI(N1,N2)
         RGQR(N1,N2)=TRGQR(N1,N2)
         RGQI(N1,N2)=TRGQI(N1,N2)
!         write(56,'(2i4,1p2e20.9)')n1,n2,QR(N1,N2),QI(N1,N2)
!         write(56,'(2i4,1p2e20.9)')n1,n2,RGQR(N1,N2),RGQI(N1,N2)
      ENDDO
   ENDDO

!  Call T-matrix calculator (uses LAPACK inversion)

   call RGQ_Qinv                          &
     ( NM, QR, QI, RGQR, RGQI,            & ! Inputs
       TR1, TI1, fail, message )            ! Outputs

!  debug

!       DO N1=1,NNMAX
!           DO N2=1,NNMAX
!         write(56,'(2i4,1p2e20.10)')n1,n2,TR1(N1,N2),TI1(N1,N2)
!       enddo
!       enddo
!      pause'56 spadgett'

!  Exception handling

   if ( fail ) then
      TRACE = ' RGQ_Qinv failed in Tmatrix_R'
   endif

!  Finish

   return
end subroutine tmatrix_R

subroutine VIG            &
      ( NMAX, X, M,       & ! Inputs
        DV1, DV2 )          ! Outputs

!  Polynomial expansions

   implicit none

!  input arguments
!  ---------------

   integer  , intent(in)  :: nmax, m
   REAL(fpk), intent(in)  :: x

!  Output arguments
!  ----------------

   REAL(fpk), intent(out) :: DV1(NPN1),DV2(NPN1)

!  Local variables
!  ---------------

   INTEGER   :: N, I, I2
   REAL(fpk) :: A, D1, D2, D3, DER
   REAL(fpk) :: QS, QS1, QMM, QN, QN1, QN2, QNM, QNM1

!  initialize

   A   = 1.0D0
   QS  = DSQRT(1.0D0-X*X)
   QS1 = 1.0D0/QS
   DO N = 1, NMAX
      DV1(N) = 0.0D0
      DV2(N) = 0.0D0
   ENDDO   

!  Two cases, M = 0, and M =/ 0

   IF ( M.NE.0) then
      QMM = DBLE(M*M)
      DO I = 1,M
         I2 = I*2
         A = A * DSQRT(DBLE(I2-1)/DBLE(I2))*QS
      ENDDO   
      D1 = 0.0D0
      D2 = A 
      DO N = M, NMAX
         QN=DBLE(N)
         QN2=DBLE(2*N+1)
         QN1=DBLE(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO   
   else
      D1 = 1.0D0
      D2 = X  
      DO N = 1, NMAX  
         QN =DBLE(N)
         QN1=DBLE(N+1)
         QN2=DBLE(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1 
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO   
   endif

!  Finish

   return
end subroutine VIG

subroutine RGQ_Qinv                            &
          ( NMAX, QR, QI, RGQR, RGQI,          & ! Inputs
            TR1, TI1, fail, message )            ! Outputs

!   CALCULATION OF THE MATRIX    T = - RG(Q) * (Q**(-1))

   implicit none

!  input arguments
!  ---------------

!  Dimensions

   integer  , intent(in)  :: nmax

!  RGQ and Q matrices

   REAL(fpk), intent(in)  ::   QR(NPN2,NPN2),  QI(NPN2,NPN2)
   REAL(fpk), intent(in)  :: RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)

!  Output arguments
!  ----------------

!  T-matrix output

   REAL(fpk), intent(out) :: TR1(NPN2,NPN2), TI1(NPN2,NPN2)

!  Exception handling

   logical       , intent(out) :: fail
   character*(*) , intent(out) :: message

!  Local variables
!  ---------------

!  Local arrays for the LAPACK calculation and final multiplication

   COMPLEX(fpk) :: ZQ(NPN2,NPN2),ZW(NPN2)
   INTEGER      :: IPIV(NPN2)

!  Other variables

   CHARACTER*2  :: C2
   INTEGER      :: I, J, K, INFO, NNMAX
   REAL(fpk)    :: AR, AI, ARR, ARI, TR, TI

!  Initialize
!  ----------

!  Exception handling

   FAIL    = .false.
   MESSAGE = ' '

!  actual dimensions

   NNMAX = 2*NMAX
 
!  Develop complex matrix

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
            !AI  = DIMAG(ZQ(K,J))
            AI  = AIMAG(ZQ(K,J))
            TR  = TR-ARR*AR+ARI*AI
            TI  = TI-ARR*AI-ARI*AR
         ENDDO
         TR1(I,J)=TR
         TI1(I,J)=TI
      ENDDO
   ENDDO

!  Finish

   return
end subroutine RGQ_Qinv

!  End module

end module tmat_makers
