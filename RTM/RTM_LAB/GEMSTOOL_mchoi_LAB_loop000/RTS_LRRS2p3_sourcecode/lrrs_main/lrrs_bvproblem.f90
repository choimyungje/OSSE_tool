! ###############################################################
! #                                                             #
! #                    THE LIDORT_RRS MODEL                     #
! #                                                             #
! #      (LInearized Discrete Ordinate Radiative Transfer)      #
! #       --         -        -        -         -              #
! #                 (Rotational Raman Scatter)                  #
! #                  -          -     -                         #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert J. D. Spurr                           #
! #                                                             #
! #  Address :     RT SOLUTIONS Inc.                            #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Version      :  2.3                                        #
! #  Release Date :  March 2011                                 #
! #                                                             #
! ###############################################################

!    #########################################################
!    #                                                       #
!    #   This Version of LIDORT_RRS comes with a GNU-style   #
!    #   license. Please read the license carefully.         #
!    #                                                       #
!    #########################################################

! ###############################################################
! #                                                             #
! #   SUBROUTINES :                                             #
! #                                                             #
! #             SURFACE_HOMOG_SETUP  (*)                        #
! #             SURFACE_BEAM_SETUP   (*)                        #
! #             BVPMATRIX_INIT                                  #
! #             BVPMATRIX_SETUP                                 #
! #             BVPCOLUMN_SETUP                                 #
! #                                                             #
! #   Note (*): Version 2.1 with BRDF setup, November 2007.     #
! #                                                             #
! ###############################################################


      MODULE lrrs_bvproblem

      PRIVATE
      PUBLIC :: SURFACE_HOMOG_SETUP,&
                SURFACE_BEAM_SETUP,&
                BVPMATRIX_INIT,&
                BVPMATRIX_SETUP,&
                BVPCOLUMN_SETUP

      CONTAINS

      SUBROUTINE SURFACE_HOMOG_SETUP &
         ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, &
           R2, FL1, FL2, BIREFLEC, &
           FOURIER, NLAYERS, NSTREAMS, &
           QUAD_STRMWGHT, XPOS, XNEG, &
           R2_HOMP, R2_HOMM )

!  Additional sums for the final surface-reflecting layer

!   Version 2.0 Lambertian only for the RRS model
!   Version 2.1 BRDF treatment introduced. 3 November 2007 by R. Spurr
!   Version 2.2 Treatment standardized. 21 November 2008 by R. Spurr

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  inclusion of surface

      LOGICAL, INTENT(IN) ::          DO_INCLUDE_SURFACE

!  Lambertian flag

      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_SURFACE

!  Fourier number, number of layers and streams

      INTEGER, INTENT(IN) ::          FOURIER, NLAYERS, NSTREAMS

!  2-deltam0 factor, Lambertian fraction factors

      REAL(FPK), INTENT(IN) :: R2, FL1, FL2

!  BRDF function

      REAL(FPK), INTENT(IN) :: BIREFLEC(MAX_STREAMS,MAX_STREAMS)

!  quadrature inputs

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT(MAX_STREAMS)

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS(MAX_2_STREAMS,MAX_STREAMS,MAX_LAYERS)
      REAL(FPK), INTENT(IN) :: XNEG(MAX_2_STREAMS,MAX_STREAMS,MAX_LAYERS)

!  Output arguments
!  ----------------

!  Reflected homogeneous solutions at ground

      REAL(FPK), INTENT(OUT) :: R2_HOMP(MAX_STREAMS,MAX_STREAMS)
      REAL(FPK), INTENT(OUT) :: R2_HOMM(MAX_STREAMS,MAX_STREAMS)

!  local variables
!  ---------------

      REAL(FPK) :: H_1, H_2, FRJ, FRAC1, FRAC2
      INTEGER ::          I, AA, J, LAY
      REAL(FPK) :: SAVP(MAX_STREAMS),SAVM(MAX_STREAMS)

!  Integrate Downward streams of Homogeneous solutions
!     For Lambertian reflectance, all streams are the same
!     Zero output values if surface flag not set

      IF ( DO_INCLUDE_SURFACE ) THEN
        LAY = NLAYERS
        IF ( FOURIER.EQ.0.AND.DO_LAMBERTIAN_SURFACE ) THEN
          FRAC1 = FL1 * R2
          DO AA = 1, NSTREAMS
            DO I = 1, NSTREAMS
              H_1 = ZERO
              H_2 = ZERO
              DO J = 1, NSTREAMS
                H_1 = H_1 + QUAD_STRMWGHT(J)*XPOS(J,AA,LAY)
                H_2 = H_2 + QUAD_STRMWGHT(J)*XNEG(J,AA,LAY)
              ENDDO
              R2_HOMP(I,AA) = H_1 * FRAC1
              R2_HOMM(I,AA) = H_2 * FRAC1
            ENDDO
          ENDDO
        ELSE
          FRAC1 = FL1 * R2
          FRAC2 = FL2 * R2
          DO AA = 1, NSTREAMS
            DO J = 1, NSTREAMS
              SAVP(J) =  QUAD_STRMWGHT(J) * XPOS(J,AA,LAY)
              SAVM(J) =  QUAD_STRMWGHT(J) * XNEG(J,AA,LAY)
            ENDDO
            DO I = 1, NSTREAMS
              H_1 = ZERO
              H_2 = ZERO
              DO J = 1, NSTREAMS
                FRJ = FRAC1 + FRAC2 * BIREFLEC(J,I)
                H_1 = H_1 + SAVP(J) * FRJ
                H_2 = H_2 + SAVM(J) * FRJ
              ENDDO
              R2_HOMP(I,AA) = H_1
              R2_HOMM(I,AA) = H_2
            ENDDO
          ENDDO
        ENDIF
      ELSE
        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            R2_HOMP(I,J) = ZERO
            R2_HOMM(I,J) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE SURFACE_HOMOG_SETUP

!

      SUBROUTINE SURFACE_BEAM_SETUP &
         ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, &
           R2, FL1, FL2, BIREFLEC, &
           FOURIER, NLAYERS, NSTREAMS, &
           QUAD_STRMWGHT, WLOWER, &
           R2_BEAM )

!  Additional sums for the final surface-reflecting layer

!   Version 2.0 Lambertian only for the RRS model
!   Version 2.1 BRDF treatment introduced. 3 November 2007 by R. Spurr
!   Version 2.2 Treatment standardized. 21 November 2008 by R. Spurr

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  inclusion of surface

      LOGICAL, INTENT(IN) ::          DO_INCLUDE_SURFACE

!  Lambertian flag

      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_SURFACE

!  Fourier number, number of layers and streams

      INTEGER, INTENT(IN) ::          FOURIER, NLAYERS, NSTREAMS

!  2-deltam0 factor, Lambertian fraction factors

      REAL(FPK), INTENT(IN) :: R2, FL1, FL2

!  BRDF function

      REAL(FPK), INTENT(IN) :: BIREFLEC(MAX_STREAMS,MAX_STREAMS)

!  quadrature inputs

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT(MAX_STREAMS)

!  particular solution at lower boundary

      REAL(FPK), INTENT(IN) :: WLOWER(MAX_2_STREAMS,MAX_LAYERS)

!  Output arguments
!  ----------------

!  Reflected beam solution at ground

      REAL(FPK), INTENT(OUT) :: R2_BEAM(MAX_STREAMS)

!  local variables
!  ---------------

      INTEGER ::          I, J, LAY
      REAL(FPK) :: SUM, FRJ, FRAC1, FRAC2

!  Integrate Downward streams of particular solutions
!    For Lambertian reflectance, all streams are the same
!    No values if surface flag not set

      IF ( DO_INCLUDE_SURFACE ) THEN
        LAY = NLAYERS
        IF ( FOURIER.EQ.0.AND.DO_LAMBERTIAN_SURFACE ) THEN
          FRAC1 = FL1 * R2
          SUM = ZERO
          DO J = 1, NSTREAMS
            SUM = SUM + QUAD_STRMWGHT(J)*WLOWER(J,LAY)
          ENDDO
          SUM = SUM * FRAC1
          DO I = 1, NSTREAMS
            R2_BEAM(I) = SUM
          ENDDO
        ELSE
          FRAC1 = FL1 * R2
          FRAC2 = FL2 * R2
          DO I = 1, NSTREAMS
            SUM = ZERO
            DO J = 1, NSTREAMS
              FRJ = FRAC1 + FRAC2 * BIREFLEC(J,I)
              SUM = SUM + FRJ * QUAD_STRMWGHT(J) * WLOWER(J,LAY)
            ENDDO
            R2_BEAM(I) = SUM
          ENDDO
        ENDIF
      ELSE
        DO I = 1, NSTREAMS
          R2_BEAM(I) = ZERO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE SURFACE_BEAM_SETUP

!

      SUBROUTINE BVPMATRIX_INIT &
           ( NSTREAMS, NLAYERS, &
             NTOTAL, NSTR2, N_SUBDIAG, N_SUPDIAG, &
             BMAT_ROWMASK, BANDMAT2 )

!  Initialise the compressed band matrix prior to solving BVP

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Computational variables

      INTEGER, INTENT(IN) ::          NSTREAMS, NLAYERS

!  Total size of BVP matrix, number of super and sub-diagonals
!   NSTR2 = 2 * NSTREAMS

      INTEGER, INTENT(IN) ::          NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTR2

!  output arguments
!  ----------------

!  Initialised band matrix

      REAL(FPK), INTENT(OUT) ::        BANDMAT2(MAX_BANDTOTAL,MAX_TOTAL)

!  Index for assignation of band matrix

      INTEGER, INTENT(OUT) ::          BMAT_ROWMASK(MAX_TOTAL,MAX_TOTAL)

!  Local variables
!  ---------------

      INTEGER ::          NMIN(MAX_TOTAL), NMAX(MAX_TOTAL)
      INTEGER ::          I, J, N3, JS, JF, IS, LF, L, KALL, I1

!mick fix 3/31/11 - initialise all components

      BANDMAT2 = ZERO
      BMAT_ROWMASK = ZERO

!  special case of 1 layer, 1 stream

      IF ( NLAYERS .EQ. 1 .AND. NSTREAMS .EQ. 1 ) THEN
        DO I = 1, NTOTAL
          DO J = 1, NTOTAL
            BANDMAT2(I,J)     = ZERO
            BMAT_ROWMASK(I,J) = I
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Crude initialisation

!        DO I = 1, NTOTAL
!          DO J = 1, NTOTAL
!            BANDMAT2(I,J)     = ZERO
!          ENDDO
!        ENDDO

!  compression row indices

      DO J = 1, N_SUPDIAG + 1
        NMIN(J) = 1
      ENDDO
      DO J = N_SUPDIAG + 2, NTOTAL
        NMIN(J) = J - N_SUPDIAG
      ENDDO
      DO J = 1, NTOTAL - N_SUBDIAG
        NMAX(J) = J + N_SUBDIAG
      ENDDO
      DO J = NTOTAL - N_SUBDIAG + 1, NTOTAL
        NMAX(J) = NTOTAL
      ENDDO

!  Set the indexing

      KALL = N_SUBDIAG + N_SUPDIAG + 1
      DO I = 1, NTOTAL
        DO J = 1, NTOTAL
          IF ( (I.GE.NMIN(J)) .AND. (I.LE.NMAX(J)) ) THEN
            BMAT_ROWMASK(I,J) = KALL + I - J
          ENDIF
        ENDDO
      ENDDO

!mick fix 3/31/11 - return now

      RETURN

!  compression matrix zeroing

      N3 = NSTR2 + NSTREAMS
      LF = NLAYERS - 2

!  upper band top

      JS = NSTR2 + 1
      JF = N3 - 1
      DO I = 1, NSTREAMS
        DO J = JS, JF + I
          BANDMAT2(BMAT_ROWMASK(I,J),J) = ZERO
        ENDDO
      ENDDO

!  upper band

      DO L = 1, LF
        IS = L*NSTR2 - NSTREAMS + 1
        JS = IS + N3
        JF = JS - 1
        DO I = 1, NSTR2-1
          I1 = I + IS
          DO J = JS, JF + I
            BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  lower band

      DO L = 1, LF
        IS = L*NSTR2 + NSTREAMS
        JS = IS - N3 + 1
        JF = IS - NSTREAMS
        DO I = 1, NSTR2-1
          I1 = I + IS
          DO J = JS + I, JF
            BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  lower band bottom

      JS = LF * NSTR2 + 1
      IS = JS + N3 - 1
      JF = IS - NSTREAMS
      DO I = 1, NSTREAMS
        I1 = I + IS
        DO J = JS + I, JF
          BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
        ENDDO
      ENDDO

!  finish

      RETURN
      END SUBROUTINE BVPMATRIX_INIT

!

      SUBROUTINE BVPMATRIX_SETUP &
           (  DO_INCLUDE_SURFACE, &
              NSTREAMS, NLAYERS, NSTR2, BMAT_ROWMASK, &
              XPOS, XNEG, R2_HOMP, R2_HOMM, T_DELT_EIGEN, &
              BANDMAT2 )

!  Fills up the compressed BVP band matrix (the "A" as in AX=B)

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Inclusion of a surface term

      LOGICAL, INTENT(IN) ::   DO_INCLUDE_SURFACE

!  Computational variables

      INTEGER, INTENT(IN) ::   NSTREAMS, NLAYERS, NSTR2

!  indexing from the BVP initialization routine

      INTEGER, INTENT(IN) ::   BMAT_ROWMASK(MAX_TOTAL,MAX_TOTAL)

!  eigenvector homogeneous solutions

      REAL(FPK), INTENT(IN) :: XPOS(MAX_2_STREAMS,MAX_STREAMS,MAX_LAYERS)
      REAL(FPK), INTENT(IN) :: XNEG(MAX_2_STREAMS,MAX_STREAMS,MAX_LAYERS)

!  surface-reflected eigensolutions

      REAL(FPK), INTENT(IN) :: R2_HOMP(MAX_STREAMS,MAX_STREAMS)
      REAL(FPK), INTENT(IN) :: R2_HOMM(MAX_STREAMS,MAX_STREAMS)

!  Whole-layer Eigenstream transmittances

      REAL(FPK), INTENT(IN) :: T_DELT_EIGEN(MAX_STREAMS,MAX_LAYERS)

!  output arguments
!  ----------------

!  Filled band matrix

      REAL(FPK), INTENT(INOUT) :: BANDMAT2(MAX_BANDTOTAL,MAX_TOTAL)

!  Local variables
!  ---------------

      INTEGER ::          I,EP,EM,N,N1,I1,AA
      INTEGER ::          C0,CE_OFFSET,CEM,CEP,CEM1,CEP1,CM,CP
      REAL(FPK) ::        XPNET, XMNET

!  Layer 1: top boundary condition: no downward diffuse radiation

      N = 1
      DO I = 1, NSTREAMS
        DO EP = 1, NSTREAMS
          EM = EP + NSTREAMS
          BANDMAT2(BMAT_ROWMASK(I,EP),EP)  = &
                         XPOS(I,EP,N)
          BANDMAT2(BMAT_ROWMASK(I,EM),EM)  = &
                         XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
        ENDDO
      ENDDO

!  intermediate layer boundaries (will not be done if NLAYERS = 1 )

      C0 = - NSTREAMS
      DO N = 2, NLAYERS
        N1 = N - 1
        C0        = C0 + NSTR2
        CE_OFFSET = C0 - NSTREAMS
        DO I = 1, NSTR2
          CM = C0 + I
          DO EP = 1, NSTREAMS
            CEP = CE_OFFSET + EP
            CEM = CEP + NSTREAMS
            CEP1 = CEP + NSTR2
            CEM1 = CEM + NSTR2
            BANDMAT2(BMAT_ROWMASK(CM,CEP),CEP)   = &
                         T_DELT_EIGEN(EP,N1)*XPOS(I,EP,N1)
            BANDMAT2(BMAT_ROWMASK(CM,CEM),CEM)   = &
                         XNEG(I,EP,N1)
            BANDMAT2(BMAT_ROWMASK(CM,CEP1),CEP1) = &
                         -XPOS(I,EP,N)
            BANDMAT2(BMAT_ROWMASK(CM,CEM1),CEM1) = &
                         -T_DELT_EIGEN(EP,N)*XNEG(I,EP,N)
          ENDDO
        ENDDO
      ENDDO

!  bottom BC (with surface reflectivity additions if flagged)

      N  = NLAYERS
      C0 = C0 + NSTR2
      CE_OFFSET = C0 - NSTREAMS

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          CP = C0 + I
          I1 = I + NSTREAMS
          DO AA = 1, NSTREAMS
            CEP = CE_OFFSET + AA
            CEM = CEP + NSTREAMS
            XPNET = XPOS(I1,AA,N) - R2_HOMP(I,AA)
            XMNET = XNEG(I1,AA,N) - R2_HOMM(I,AA)
            BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) = &
                         T_DELT_EIGEN(AA,N) * XPNET
            BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) = &
                         XMNET
          ENDDO
          
        ENDDO
      ELSE
        DO I = 1, NSTREAMS
          CP = C0 + I
          I1 = I + NSTREAMS
          DO AA = 1, NSTREAMS
            CEP = CE_OFFSET + AA
            CEM = CEP + NSTREAMS
            BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) = &
                        T_DELT_EIGEN(AA,N) * XPOS(I1,AA,N)
            BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) = &
                        XNEG(I1,AA,N)
          ENDDO
        ENDDO
      ENDIF

!  finish

      RETURN
      END SUBROUTINE BVPMATRIX_SETUP

!

      SUBROUTINE BVPCOLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, &
             NSTREAMS, NLAYERS, NTOTAL, NSTR2, &
             DIRECT_BEAM, WUPPER, WLOWER, R2_BEAM, &
             COL2 )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Inclusion of a surface term

      LOGICAL, INTENT(IN) ::          DO_INCLUDE_SURFACE

!  Computational variables

      INTEGER, INTENT(IN) ::          NSTREAMS, NLAYERS
      INTEGER, INTENT(IN) ::          NTOTAL, NSTR2

!  Direct beam

      REAL(FPK), INTENT(IN) :: DIRECT_BEAM(MAX_STREAMS)

!  Beam solutions at upper and lower layer boundaries

      REAL(FPK), INTENT(IN) :: WUPPER(MAX_2_STREAMS,MAX_LAYERS)
      REAL(FPK), INTENT(IN) :: WLOWER(MAX_2_STREAMS,MAX_LAYERS)

!  Surface-reflected beam solution

      REAL(FPK), INTENT(IN) :: R2_BEAM(MAX_STREAMS)

!  output arguments
!  ----------------

!  Column vector

      REAL(FPK), INTENT(OUT) :: COL2(MAX_TOTAL,1)

!  Local variables
!  ---------------

      INTEGER ::          I,LAY,LAY1,I1,C0,CM
      REAL(FPK) :: COL2_D(MAX_TOTAL)

!  zero column vector

      DO I = 1, NTOTAL
        COL2(I,1) = ZERO
        COL2_D(I) = ZERO
      ENDDO

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------------

      LAY = 1
      DO I = 1, NSTREAMS
        COL2_D(I)   = - WUPPER(I,LAY)
      ENDDO

!  intermediate layer boundaries (will not be done if NLAYERS = 1 )
!  -----------------------------

      DO LAY = 2, NLAYERS

        LAY1 = LAY - 1
        C0 = LAY1*NSTR2 - NSTREAMS
        DO I = 1, NSTR2
          CM = C0 + I
          COL2_D(CM) = WUPPER(I,LAY) - WLOWER(I,LAY1)
        ENDDO

      ENDDO

!  lowest (surface) reflecting boundary (diffuse radiation terms only)
!  ------------------------------------

      LAY = NLAYERS
      C0 = (LAY-1)*NSTR2 + NSTREAMS

!  with non-zero surface reflectance, need integrated downward reflectan

      IF ( DO_INCLUDE_SURFACE ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          CM = C0 + I
          COL2_D(CM) = - WLOWER(I1,LAY) + R2_BEAM(I)
        ENDDO

!  no albedo, similar code excluding integrated reflectance

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          CM = C0 + I
          COL2_D(CM) = - WLOWER(I1,LAY)
        ENDDO

      ENDIF

!  Add direct beam solution (only to final level)
!  ----------------------------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          COL2_D(CM) = COL2_D(CM) + DIRECT_BEAM(I)
        ENDDO
      ENDIF

!  Add to existing solution

      DO I = 1, NTOTAL
        COL2(I,1) = COL2(I,1) + COL2_D(I)
      ENDDO

!  finish

      RETURN
      END SUBROUTINE BVPCOLUMN_SETUP

      END MODULE lrrs_bvproblem

