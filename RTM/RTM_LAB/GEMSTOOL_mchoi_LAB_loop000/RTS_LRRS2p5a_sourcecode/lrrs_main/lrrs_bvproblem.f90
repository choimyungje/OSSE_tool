
! ###############################################################
! #                                                             #
! #                    THE LIDORT_RRS MODEL                     #
! #                                                             #
! #      (LInearized Discrete Ordinate Radiative Transfer)      #
! #       --         -        -        -         -              #
! #                 (Rotational Raman Scattering)               #
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
! #  Email   :     rtsolutions@verizon.net                      #
! #  Website :     www.rtslidort.com                            #
! #                                                             #
! #  Version  #   :  2.5                                        #
! #  Release Date :  March 2017                                 #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  --- History of the model ------------                      #
! #                                                             #
! #  Version 1.0 : 2005, Fortran 77                             #
! #  Version 1.1 : 2007, F77                                    #
! #                                                             #
! #  Version 2.1 : 2009, F77                                    #
! #       * Linearization for Atmospheric-property Jacobians    #
! #       * Single scatter corrections added                    #
! #                                                             #
! #  Version 2.3 : March 2011, Fortran 90                       #
! #       * Simplified Raman-setup procedure                    #
! #       * F90 Version with Type-structure I/O                 #
! #       * Test package developed for installation             #
! #                                                             #
! #  Version 2.5 : March 2017, F90                              #
! #       * Formal BRDF/SLEAVE supplements developed            #
! #       * New test-bed software for testing supplements       #
! #       * Thread-safe Code for OpenMP applications            #
! #       * Complete revision of Taylor-series modules          #
! #       * New User Guide and Review paper                     #
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

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Bookkeeping improvements (use of "Only", clearer I/O specifications)
!    (2) Specific surface input replaced by general supplement-derived BRDF stuff

      MODULE lrrs_bvproblem_m

!      USE LRRS_PARS_m, Only : SDU

      PRIVATE
      PUBLIC :: SURFACE_HOMOG_SETUP, &
                SURFACE_BEAM_SETUP,  &
                BVPMATRIX_INIT,      &
                BVPMATRIX_SETUP,     &
                BVPCOLUMN_SETUP

      CONTAINS

      SUBROUTINE SURFACE_HOMOG_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR, & ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX,     & ! Inputs
           QUAD_STRMWGHT, ALBEDOS_RANKED, BRDF_F, XPOS, XNEG,   & ! Inputs
           R2_HOMP, R2_HOMM )                                     ! Outputs

!  Additional sums for the final surface-reflecting layer

!   Version 2.0 Lambertian only for the RRS model
!   Version 2.1 BRDF treatment introduced. 3 November 2007 by R. Spurr
!   Version 2.2 Treatment standardized. 21 November 2008 by R. Spurr
!   Version 2.5 Specific surface input replaced by general supplement-derived BRDFs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_MOMENTS, MAX_STREAMS, MAX_2_STREAMS, &
                              MAX_LAYERS, MAX_POINTS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  inclusion of surface

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  BRDF flag (replaces Lamabertian flag, Version 2.5, 9/11/15)

      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE

!  Fourier number, number of layers and streams

      INTEGER  , INTENT(IN) :: FOURIER, NLAYERS, NSTREAMS

!  Indices, Version 2.5, 9/11/15

      INTEGER  , INTENT(IN) :: Raman_IDX, Brdf_IDX

!  2-deltam0 factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Albedos

      REAL(fpk), intent(in)  :: ALBEDOS_RANKED ( MAX_POINTS)

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )

!  quadrature inputs

      REAL(FPK), INTENT(IN)  :: QUAD_STRMWGHT(MAX_STREAMS)

!  Eigenvector solutions

      REAL(FPK), INTENT(IN)  :: XPOS(MAX_2_STREAMS,MAX_STREAMS,MAX_LAYERS)
      REAL(FPK), INTENT(IN)  :: XNEG(MAX_2_STREAMS,MAX_STREAMS,MAX_LAYERS)

!  Output arguments
!  ----------------

!  Reflected homogeneous solutions at ground

      REAL(FPK), INTENT(OUT) :: R2_HOMP(MAX_STREAMS,MAX_STREAMS)
      REAL(FPK), INTENT(OUT) :: R2_HOMM(MAX_STREAMS,MAX_STREAMS)

!  local variables
!  ---------------

      REAL(FPK) :: FACTOR, REFL_P, REFL_M
      INTEGER   :: I, AA, J, LAY

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        DO AA = 1, NSTREAMS
          R2_HOMP(I,AA) = ZERO
          R2_HOMM(I,AA) = ZERO
        ENDDO
      ENDDO

!  Return with Zeroed values if surface flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Layer

      LAY = NLAYERS

!  For Lambertian reflectance, all streams are the same

      IF ( .not. DO_BRDF_SURFACE ) THEN
        FACTOR = SURFACE_FACTOR * ALBEDOS_RANKED(Raman_IDX)
        IF ( FOURIER .EQ. 0 ) THEN
          DO AA = 1, NSTREAMS
            REFL_P = ZERO
            REFL_M = ZERO
            DO J = 1, NSTREAMS
              REFL_P = REFL_P + QUAD_STRMWGHT(J) * XPOS(J,AA,LAY)
              REFL_M = REFL_M + QUAD_STRMWGHT(J) * XNEG(J,AA,LAY)
            ENDDO
            REFL_P = REFL_P * FACTOR
            REFL_M = REFL_M * FACTOR
            DO I = 1, NSTREAMS
              R2_HOMP(I,AA) = REFL_P
              R2_HOMM(I,AA) = REFL_M
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  For BRDF surface

      IF ( DO_BRDF_SURFACE ) THEN
        DO AA = 1, NSTREAMS
          DO I = 1, NSTREAMS
            REFL_P = ZERO
            REFL_M = ZERO
            DO J = 1, NSTREAMS
              REFL_P = REFL_P + QUAD_STRMWGHT(J) * XPOS(J,AA,LAY) * BRDF_F(FOURIER,I,J,Brdf_IDX)
              REFL_M = REFL_M + QUAD_STRMWGHT(J) * XNEG(J,AA,LAY) * BRDF_F(FOURIER,I,J,Brdf_IDX)
            ENDDO
            R2_HOMP(I,AA) = REFL_P * SURFACE_FACTOR
            R2_HOMM(I,AA) = REFL_M * SURFACE_FACTOR
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE SURFACE_HOMOG_SETUP

!

      SUBROUTINE SURFACE_BEAM_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR, & ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX,     & ! Inputs
           QUAD_STRMWGHT, ALBEDOS_RANKED, BRDF_F, WLOWER,       & ! Inputs
           R2_BEAM )                                              ! Output

!  Additional sums for the final surface-reflecting layer

!   Version 2.0 Lambertian only for the RRS model
!   Version 2.1 BRDF treatment introduced. 3 November 2007 by R. Spurr
!   Version 2.2 Treatment standardized. 21 November 2008 by R. Spurr
!   Version 2.5 Specific surface input replaced by general supplement-derived BRDFs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_MOMENTS, MAX_STREAMS, MAX_2_STREAMS, &
                              MAX_LAYERS, MAX_POINTS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  inclusion of surface

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  BRDF flag (replaces Lamabertian flag, Version 2.5, 9/11/15)

      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE

!  Fourier number, number of layers and streams

      INTEGER  , INTENT(IN) :: FOURIER, NLAYERS, NSTREAMS

!  Indices, Version 2.5, 9/11/15

      INTEGER  , INTENT(IN) :: Raman_IDX, Brdf_IDX

!  2-deltam0 factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Albedos

      REAL(fpk), intent(in)  :: ALBEDOS_RANKED ( MAX_POINTS)

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )

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

      REAL(fpk)   :: REFL_B!, FACTOR
      INTEGER     :: I, J, LAY

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        R2_BEAM(I) = ZERO
      ENDDO

!  Return with Zeroed values if surface flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Surface layer

      LAY = NLAYERS

!  For Lambertian reflectance, all streams are the same
!  ----------------------------------------------------

      IF ( .not. DO_BRDF_SURFACE ) THEN
        !FACTOR = SURFACE_FACTOR * ALBEDOS_RANKED(Raman_IDX)
        IF ( FOURIER .EQ. 0 ) THEN
          !mick note - terms of DOT PROD here = elements of LIDORT's H_WLOWER
          !REFL_B = FACTOR * DOT_PRODUCT(QUAD_STRMWGHT(1:NSTREAMS),WLOWER(1:NSTREAMS,LAY))
          !mick note - R2_BEAM here  = LIDORT's R2_PARTIC
          !R2_BEAM(1:NSTREAMS) = REFL_B
          REFL_B = ALBEDOS_RANKED(Raman_IDX) &
                 * DOT_PRODUCT(QUAD_STRMWGHT(1:NSTREAMS),WLOWER(1:NSTREAMS,LAY))
          R2_BEAM(1:NSTREAMS) = REFL_B * SURFACE_FACTOR
        ENDIF
      ENDIF

!  For bidirectional reflecting surface
!  ------------------------------------

      IF ( DO_BRDF_SURFACE ) THEN
        DO I = 1, NSTREAMS
          REFL_B = ZERO
          DO J = 1, NSTREAMS
            REFL_B = REFL_B + BRDF_F(FOURIER,I,J,Brdf_IDX) * QUAD_STRMWGHT(J) * WLOWER(J,LAY)
          ENDDO
          R2_BEAM(I) = REFL_B * SURFACE_FACTOR
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE SURFACE_BEAM_SETUP

!

      SUBROUTINE BVPMATRIX_INIT &
           ( NSTREAMS, NLAYERS,                   & ! Inputs
             NTOTAL, NSTR2, N_SUBDIAG, N_SUPDIAG, & ! Inputs
             BMAT_ROWMASK, BANDMAT2 )

!  Initialise the compressed band matrix prior to solving BVP

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_BANDTOTAL,MAX_TOTAL

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Computational variables

      INTEGER  , INTENT(IN) :: NSTREAMS, NLAYERS

!  Total size of BVP matrix, number of super and sub-diagonals
!   NSTR2 = 2 * NSTREAMS

      INTEGER  , INTENT(IN) :: NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTR2

!  output arguments
!  ----------------

!  Initialised band matrix

      REAL(FPK), INTENT(OUT) :: BANDMAT2(MAX_BANDTOTAL,MAX_TOTAL)

!  Index for assignation of band matrix

      INTEGER  , INTENT(OUT) :: BMAT_ROWMASK(MAX_TOTAL,MAX_TOTAL)

!  Local variables
!  ---------------

      INTEGER :: NMIN(MAX_TOTAL), NMAX(MAX_TOTAL)
      INTEGER :: I, J, N3, JS, JF, IS, LF, L, KALL, I1

!mick fix 3/31/11 - initialise all components

      BANDMAT2     = ZERO
      BMAT_ROWMASK = 0

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
           (  DO_INCLUDE_SURFACE,                         & ! Inputs
              NSTREAMS, NLAYERS, NSTR2, BMAT_ROWMASK,     & ! Inputs
              XPOS, XNEG, R2_HOMP, R2_HOMM, T_DELT_EIGEN, & ! Inputs
              BANDMAT2, SMAT2 )                             ! Output

!  Fills up the compressed BVP band matrix (the "A" as in AX=B)
!  LRRS Version 2.5, Added SMAT2 argument for the single layer case

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_BANDTOTAL, MAX_TOTAL

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Inclusion of a surface term

      LOGICAL  , INTENT(IN) ::   DO_INCLUDE_SURFACE

!  Computational variables

      INTEGER  , INTENT(IN) ::   NSTREAMS, NLAYERS, NSTR2

!  indexing from the BVP initialization routine

      INTEGER  , INTENT(IN) ::   BMAT_ROWMASK(MAX_TOTAL,MAX_TOTAL)

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

!  square matrix for the single layer case

      REAL(fpk), intent(inout) :: SMAT2   (MAX_2_STREAMS,MAX_2_STREAMS)

!  Local variables
!  ---------------

      INTEGER   :: I,EP,EM,N,N1,I1,AA
      INTEGER   :: C0,CE_OFFSET,CEM,CEP,CEM1,CEP1,CM,CP
      REAL(FPK) :: XPNET, XMNET

!   General case

      IF ( NLAYERS .gt. 1 ) THEN

!  Layer 1: top boundary condition: no downward diffuse radiation

         N = 1
         DO I = 1, NSTREAMS
          DO EP = 1, NSTREAMS
            EM = EP + NSTREAMS
            BANDMAT2(BMAT_ROWMASK(I,EP),EP)  = XPOS(I,EP,N)
            BANDMAT2(BMAT_ROWMASK(I,EM),EM)  = XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
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
              BANDMAT2(BMAT_ROWMASK(CM,CEP),CEP)   = T_DELT_EIGEN(EP,N1)*XPOS(I,EP,N1)
              BANDMAT2(BMAT_ROWMASK(CM,CEM),CEM)   = XNEG(I,EP,N1)
              BANDMAT2(BMAT_ROWMASK(CM,CEP1),CEP1) = -XPOS(I,EP,N)
              BANDMAT2(BMAT_ROWMASK(CM,CEM1),CEM1) = -T_DELT_EIGEN(EP,N)*XNEG(I,EP,N)
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
              BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) = T_DELT_EIGEN(AA,N) * XPNET
              BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) = XMNET
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            CP = C0 + I
            I1 = I + NSTREAMS
            DO AA = 1, NSTREAMS
              CEP = CE_OFFSET + AA
              CEM = CEP + NSTREAMS
              BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) = T_DELT_EIGEN(AA,N) * XPOS(I1,AA,N)
              BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) = XNEG(I1,AA,N)
            ENDDO
          ENDDO
        ENDIF

!  Single-layer case
!  =================

      ELSE

!  top BC for layer 1: no downward diffuse radiation

        N = 1
        DO I = 1, NSTREAMS
          DO EP = 1, NSTREAMS
            EM = EP + NSTREAMS
            SMAT2(I,EP) = XPOS(I,EP,N)
            SMAT2(I,EM) = XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
          ENDDO
        ENDDO

!  bottom B! (with albedo additions)

        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            CP = NSTREAMS + I
            DO EP = 1, NSTREAMS
              CEP = EP
              CEM = CEP + NSTREAMS
              XPNET = XPOS(I1,EP,N) - R2_HOMP(I,EP)
              XMNET = XNEG(I1,EP,N) - R2_HOMM(I,EP)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            CP = NSTREAMS + I
            DO EP = 1, NSTREAMS
              CEP = EP
              CEM = CEP + NSTREAMS
              XPNET = XPOS(I1,EP,N)
              XMNET = XNEG(I1,EP,N)
              SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
              SMAT2(CP,CEM) = XMNET
            ENDDO
          ENDDO
        ENDIF

!  End clause

      ENDIF

!  finish

      RETURN
      END SUBROUTINE BVPMATRIX_SETUP

!

      SUBROUTINE BVPCOLUMN_SETUP &
           ( ELASTIC_CALL, DO_INCLUDE_SURFACE,     & ! inputs
             NSTREAMS, NLAYERS, NTOTAL, NSTR2,     & ! inputs
             DIRECT_BEAM, WUPPER, WLOWER, R2_BEAM, & ! inputs
             COL2, SCOL2 )                           ! Outputs

!  Version 2.5, 9/11/15. Added SCOL2. Output declared "InOut" for OpenMP
!               2/6/17.  removed GOTO statement

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_TOTAL

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Inclusion of a surface term

      LOGICAL  , INTENT(IN) :: ELASTIC_CALL

!  Inclusion of a surface term

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  Computational variables

      INTEGER  , INTENT(IN) :: NSTREAMS, NLAYERS
      INTEGER  , INTENT(IN) :: NTOTAL, NSTR2

!  Direct beam

      REAL(FPK), INTENT(IN) :: DIRECT_BEAM(MAX_STREAMS)

!  Beam solutions at upper and lower layer boundaries

      REAL(FPK), INTENT(IN) :: WUPPER(MAX_2_STREAMS,MAX_LAYERS)
      REAL(FPK), INTENT(IN) :: WLOWER(MAX_2_STREAMS,MAX_LAYERS)

!  Surface-reflected beam solution

      REAL(FPK), INTENT(IN) :: R2_BEAM(MAX_STREAMS)

!  output arguments
!  ----------------

!  Column vectors

      REAL(FPK), INTENT(inout) :: COL2 (MAX_TOTAL,1)
      REAL(fpk), intent(inout) :: SCOL2(MAX_2_STREAMS,1)

!  Local variables
!  ---------------

      INTEGER   :: I,LAY,LAY1,I1,C0,CM
      REAL(FPK) :: COL2_D(MAX_TOTAL), SCOL2_D(MAX_2_STREAMS)

!  If Nlayers = 1, special case (below). GOTO 788 removed 2/6/17.
!      IF ( NLAYERS .EQ. 1 ) GOTO 788

!  General case
!  ============

      IF ( NLAYERS .GT. 1 ) THEN

!  Zero column vector

        DO I = 1, NTOTAL
          COL2(I,1) = ZERO
          COL2_D(I) = ZERO
        ENDDO

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------------

        LAY = 1
        DO I = 1, NSTREAMS
          COL2_D(I) = - WUPPER(I,LAY)
        ENDDO

!  Intermediate layer boundaries (will not be done if NLAYERS = 1 )
!  ----------------------------------------------------------------

        DO LAY = 2, NLAYERS
          LAY1 = LAY - 1
          C0 = LAY1*NSTR2 - NSTREAMS
          DO I = 1, NSTR2
            CM = C0 + I
            COL2_D(CM) = WUPPER(I,LAY) - WLOWER(I,LAY1)
          ENDDO
        ENDDO

!  Lowest (surface) reflecting boundary (diffuse radiation terms only)
!  -------------------------------------------------------------------

        LAY = NLAYERS
        C0 = (LAY-1)*NSTR2 + NSTREAMS

!  With non-zero surface reflectance, need integrated downward reflectan

        IF ( DO_INCLUDE_SURFACE ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            CM = C0 + I
            COL2_D(CM) = -WLOWER(I1,LAY) + R2_BEAM(I)
          ENDDO

!  With no surface, similar code excluding integrated reflectance

        ELSE

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            CM = C0 + I
            COL2_D(CM) = -WLOWER(I1,LAY)
          ENDDO

        ENDIF

!  Add direct beam solution (only to final level)
!  ----------------------------------------

        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            CM = C0 + I
            COL2_D(CM) = COL2_D(CM) + DIRECT_BEAM(I)
          ENDDO
        ENDIF

!  Add to existing solution.
!     Note 9/11/15. Why do we need this ???????

        DO I = 1, NTOTAL
          COL2(I,1) = COL2(I,1) + COL2_D(I)
        ENDDO

!  Continuation point for single layer case. Removed, 2/6/17.
!      RETURN
!788   continue

!  Single-layer case
!  =================

      ELSE

!  Zero column vector

        DO I = 1, NTOTAL
          SCOL2_D(I) = ZERO
        ENDDO

!  Upper boundary for layer 1: no downward diffuse radiation

        LAY = 1
        DO I = 1, NSTREAMS
          SCOL2_D(I)   = - WUPPER(I,LAY)
        ENDDO

!  Lowest (surface) boundary with albedo (diffuse radiation terms only)
!     With non-zero albedo, include integrated downward reflectances
!     With no albedo, similar code excluding integrated reflectance

        C0 = NSTREAMS
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            CM = C0 + I
            SCOL2_D(CM) = - WLOWER(I1,LAY) + R2_BEAM(I)
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            CM = C0 + I
            SCOL2_D(CM) = - WLOWER(I1,LAY)
          ENDDO
        ENDIF

!  Add direct beam solution (only to final level)

        IF ( DO_INCLUDE_SURFACE ) THEN
          DO I = 1, NSTREAMS
            CM = C0 + I
            SCOL2_D(CM) = SCOL2_D(CM) + DIRECT_BEAM(I)
          ENDDO
        ENDIF

!  Add to existing solution.
!     Note 9/11/15. Why do we need this ???????

        DO I = 1, NSTR2
          SCOL2(I,1) = SCOL2(I,1) + SCOL2_D(I)
        ENDDO

!  End layer clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BVPCOLUMN_SETUP

!  End module

      END MODULE lrrs_bvproblem_m

