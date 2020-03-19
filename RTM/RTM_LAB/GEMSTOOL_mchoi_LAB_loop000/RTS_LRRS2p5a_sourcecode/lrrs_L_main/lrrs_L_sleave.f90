
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

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #     Top level PUBLIC routines--------------            #
! #            LRRS_LSSL_DBSETUPS  (diffuse)               #
! #            LRRS_LSSL_UDBSETUPS (diffuse, user)         #
! #            LRRS_LSSL_WFS       (diffuse, master)       #
! #                                                        #
! #     High level Jacobian routines  ---------            #
! #            LSSL_UPUSER_SURFACEWF                       #
! #            LSSL_DNUSER_SURFACEWF                       #
! #            LSSL_INTEGRATED_OUTPUT (master)             #
! #                                                        #
! #     Post-processing at user angles --------            #
! #            LSSL_WHOLELAYER_STERM_UP                    #
! #            LSSL_WHOLELAYER_STERM_DN                    #
! #                                                        #
! #     Post-processing at quad angles --------            #
! #            LSSL_QUADWF_LEVEL_UP                        #
! #            LSSL_QUADWF_LEVEL_DN                        #
! #                                                        #
! ##########################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Use of Supplement-derived BRDF/SLEAVE inputs and control
!    (2) Use of Taylor-series control, TAYLOR_ORDER index, passed to PostProcessing
!    (3) Introduction of variable number of surface Jacobians, dimensioning to match
!    (4) Bookkeeping improvements (use of "Only", clearer I/O specifications)

      MODULE lrrs_L_sleave_m

!      USE LRRS_PARS_m, Only : LDU

!  Only 3 routines available publicly to the rest of LRRS...

      PRIVATE
      PUBLIC :: LRRS_LSSL_DBSETUPS, &
                LRRS_LSSL_UDBSETUPS, &
                LRRS_LSSL_WFS

      CONTAINS

      SUBROUTINE LRRS_LSSL_DBSETUPS ( &
        DO_SL_ISOTROPIC, DO_REFLECTED_DIRECTBEAM, FOURIER, & ! Input
        Sleave_IDX, NSTREAMS, N_SLEAVE_WFS,                & ! Input
        FLUX_FACTOR, DELTA_FACTOR,                         & ! Input
        LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,            & ! Input
        LSSL_DIRECT_BEAM )                                   ! Output

      USE LRRS_PARS_m, only : fpk, MAX_SLEAVEWFS, MAX_STREAMS, MAX_USER_STREAMS, &
                              MAX_MOMENTS, MAX_POINTS,  ZERO

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_SL_ISOTROPIC
      LOGICAL, INTENT (IN) ::           DO_REFLECTED_DIRECTBEAM

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           Sleave_IDX
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

      REAL(fpk), INTENT (IN) ::         FLUX_FACTOR, DELTA_FACTOR
      REAL(fpk), INTENT (IN) ::         LSSL_SLTERM_ISOTROPIC &
          ( MAX_SLEAVEWFS, MAX_POINTS )
      REAL(fpk), INTENT (IN) ::         LSSL_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )

!  Outputs

      REAL(fpk), INTENT (OUT) ::        LSSL_DIRECT_BEAM ( MAX_SLEAVEWFS, MAX_STREAMS)

!  Local variables
!  ---------------

      REAL(fpk) :: SL, HELP
      INTEGER   :: I, UI, O1, M, Q

!  Initialize
!  ----------

!  Safety first!

      DO I = 1, NSTREAMS
        LSSL_DIRECT_BEAM(1:N_SLEAVE_WFS,I) = ZERO
      ENDDO

!  Fourier component

      M = FOURIER

!  Corrected implementation, 30 July 2012
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases

      IF ( DO_REFLECTED_DIRECTBEAM ) THEN
!mick fix 7/20/2016 - normalize by PI4 as in DBEAM_SETUP (on hold)
        HELP = FLUX_FACTOR / DELTA_FACTOR
        !HELP = PI4 * FLUX_FACTOR / DELTA_FACTOR
        IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
          DO Q = 1, N_SLEAVE_WFS
            SL = LSSL_SLTERM_ISOTROPIC(Q,Sleave_IDX) * HELP
            LSSL_DIRECT_BEAM(Q,1:NSTREAMS) = SL
          ENDDO
        ELSE
          DO Q = 1, N_SLEAVE_WFS
            DO I = 1, NSTREAMS
              SL = LSSL_SLTERM_F_0(Q,M,I,Sleave_IDX) * HELP
              LSSL_DIRECT_BEAM(Q,I) = SL
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  finish

      RETURN
      END SUBROUTINE LRRS_LSSL_DBSETUPS

!

      SUBROUTINE LRRS_LSSL_UDBSETUPS ( &
        DO_SL_ISOTROPIC, DO_REFLECTED_DIRECTBEAM, FOURIER,  &
        Sleave_IDX, N_USER_STREAMS, N_SLEAVE_WFS,           &
        FLUX_FACTOR, DELTA_FACTOR,                          &
        LSSL_SLTERM_ISOTROPIC, LSSL_USER_SLTERM_F_0,        &
        LSSL_USER_DIRECT_BEAM )

      USE LRRS_PARS_m, only : fpk, MAX_SLEAVEWFS, MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS, ZERO

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_SL_ISOTROPIC
      LOGICAL, INTENT (IN) ::           DO_REFLECTED_DIRECTBEAM

      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           Sleave_IDX
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

      REAL(fpk), INTENT (IN) ::         FLUX_FACTOR, DELTA_FACTOR
      REAL(fpk), INTENT (IN) ::         LSSL_SLTERM_ISOTROPIC &
          ( MAX_SLEAVEWFS, MAX_POINTS )
      REAL(fpk), INTENT (IN) ::         LSSL_USER_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )

!  Outputs

      REAL(fpk), INTENT (OUT) ::        LSSL_USER_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  Local variables
!  ---------------

      REAL(fpk) :: SL, HELP
      INTEGER   :: I, UI, O1, M, Q

!  Initialize
!  ----------

!  Safety first!

      DO UI = 1, N_USER_STREAMS
        LSSL_USER_DIRECT_BEAM(1:N_SLEAVE_WFS,UI) = ZERO
      ENDDO

!  Fourier component

      M = FOURIER

!  Corrected implementation, 30 July 2012
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases

      IF ( DO_REFLECTED_DIRECTBEAM ) THEN
!mick fix 7/20/2016 - normalize by PI4 as in UDBEAM_SETUP (on hold)
        HELP = FLUX_FACTOR / DELTA_FACTOR
        !HELP = PI4 * FLUX_FACTOR / DELTA_FACTOR
        IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
          DO Q = 1, N_SLEAVE_WFS
            SL = LSSL_SLTERM_ISOTROPIC(Q,Sleave_IDX) * HELP
            LSSL_USER_DIRECT_BEAM(Q,1:N_USER_STREAMS) = SL
          ENDDO
        ELSE
          DO Q = 1, N_SLEAVE_WFS
            DO UI = 1, N_USER_STREAMS
              SL = LSSL_USER_SLTERM_F_0(Q,M,UI,Sleave_IDX) * HELP
              LSSL_USER_DIRECT_BEAM(Q,UI) = SL
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  finish

      RETURN
      END SUBROUTINE LRRS_LSSL_UDBSETUPS

!

      SUBROUTINE LRRS_LSSL_WFS ( &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_MVOUT,                              & ! Inputs
        DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_SL_ISOTROPIC,         & ! Inputs
        FOURIER, Raman_IDX, Brdf_IDX, Sleave_IDX,                             & ! Inputs
        NLAYERS, NSTREAMS, N_USER_STREAMS, N_LOUTPUT,                         & ! Inputs
        NSTR2, N_SUBDIAG, N_SUPDIAG, NTOTAL, N_SURFACE_WFS, N_SLEAVE_WFS,     & ! Inputs
        LEVEL_MASK_UP, LEVEL_MASK_DN, ALBEDOS_RANKED, USER_BRDF_F,            & ! Inputs
        SURFACE_FACTOR, FLUX_MULTIPLIER, QUAD_WEIGHTS, QUAD_STRMWGT,          & ! Inputs
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT,                                     & ! Inputs
        T_DELT_EIGEN, SOLA_XPOS, SOLB_XNEG, HMULT_1, HMULT_2,                 & ! Inputs
        T_DELT_USERM, U_XPOS, U_XNEG,                                         & ! Inputs
        LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM,                              & ! Inputs
        LS_ELASTIC_F_UP,    LS_ELASTIC_F_DN,                                  & ! Outputs
        LS_MEAN_ELASTIC_UP, LS_MEAN_ELASTIC_DN,                               & ! Outputs
        LS_FLUX_ELASTIC_UP, LS_FLUX_ELASTIC_DN,                               & ! Outputs
        STATUS, MESSAGE, TRACE )                                                ! Outputs

      USE LRRS_PARS_m, ONLY : fpk, MAX_LOUTPUT, MAX_MOMENTS, MAX_SLEAVEWFS, MAX_SURFACEWFS, MAX_POINTS, &
                              MAX_STREAMS, MAX_2_STREAMS, MAX_USER_STREAMS, MAX_LAYERS, MAX_DIRECTIONS, &
                              MAX_BANDTOTAL, MAX_TOTAL, LRRS_SUCCESS, LRRS_SERIOUS, ZERO, ONE, UPIDX, DNIDX

      USE LRRS_AUX2_m, ONLY : DGBTRS, DGETRS

      IMPLICIT NONE

!  Control

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUT
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_BRDF_SURFACE
      LOGICAL, INTENT (IN) ::          DO_SL_ISOTROPIC

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          Raman_IDX, Brdf_IDX, Sleave_IDX

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          N_LOUTPUT

      INTEGER, INTENT (IN) ::          NSTR2
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SURFACE_WFS, N_SLEAVE_WFS 

      INTEGER, INTENT (IN) ::          LEVEL_MASK_UP  ( MAX_LOUTPUT )
      INTEGER, INTENT (IN) ::          LEVEL_MASK_DN  ( MAX_LOUTPUT )

      REAL(fpk), INTENT(IN)  ::        ALBEDOS_RANKED ( MAX_POINTS )
      REAL(fpk), INTENT (IN) ::        USER_BRDF_F &
          ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Solutions

      REAL(fpk), INTENT (IN) ::        SURFACE_FACTOR
      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER

      REAL(fpk), INTENT (IN) ::        QUAD_WEIGHTS ( MAX_STREAMS )
      REAL(fpk), INTENT (IN) ::        QUAD_STRMWGT ( MAX_STREAMS )

      REAL(fpk), INTENT (IN) ::        BANDMAT2 ( MAX_BANDTOTAL, MAX_TOTAL )
      INTEGER, INTENT (IN) ::          IPIVOT   ( MAX_TOTAL )
      REAL(fpk), INTENT (IN) ::        SMAT2    ( MAX_2_STREAMS, MAX_2_STREAMS )
      INTEGER, INTENT (IN) ::          SIPIVOT  ( MAX_2_STREAMS )

      REAL(fpk), INTENT (IN) ::        T_DELT_EIGEN ( MAX_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        SOLA_XPOS &
          ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        SOLB_XNEG &
          ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        HMULT_1 &
          ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_2 &
          ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        T_DELT_USERM &
          ( MAX_LAYERS, MAX_USER_STREAMS )

      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!mick fix 9/7/2012 - moved MAX_SLEAVEWFS to 1st dimension
      REAL(fpk), INTENT (IN) ::        LSSL_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAX_STREAMS )
      REAL(fpk), INTENT (IN) ::        LSSL_USER_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  Linearized surface output

      REAL(FPK), INTENT(INOUT) :: LS_ELASTIC_F_UP &
          ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS )
      REAL(FPK), INTENT(INOUT) :: LS_ELASTIC_F_DN &
          ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS )

      REAL(FPK), INTENT(INOUT) :: LS_MEAN_ELASTIC_UP &
          ( MAX_SURFACEWFS, MAX_LOUTPUT )
      REAL(FPK), INTENT(INOUT) :: LS_MEAN_ELASTIC_DN &
          ( MAX_SURFACEWFS, MAX_LOUTPUT )

      REAL(FPK), INTENT(INOUT) :: LS_FLUX_ELASTIC_UP &
          ( MAX_SURFACEWFS, MAX_LOUTPUT )
      REAL(FPK), INTENT(INOUT) :: LS_FLUX_ELASTIC_DN &
          ( MAX_SURFACEWFS, MAX_LOUTPUT )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  Linearized BOA terms

      REAL(fpk) ::        LSSL_BOA_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  Linearized BVP solution

      REAL(fpk) ::        COL2_WFSLEAVE  ( MAX_TOTAL, MAX_SLEAVEWFS )
      REAL(fpk) ::        SCOL2_WFSLEAVE ( MAX_2_STREAMS, MAX_SLEAVEWFS )

      REAL(fpk) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )

!  Linearized surface output

      REAL(fpk) ::  SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_DIRECTIONS )
      REAL(fpk) ::  MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_DIRECTIONS )
      REAL(fpk)  :: FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_DIRECTIONS )

!  Other local variables

      INTEGER ::           INFO, M, N, Q, O1
      INTEGER ::           K, K0, K1, K2, KO1, C0, CM, I, UM, LB, UB
      INTEGER ::           IR, IROW, IROW1, IROW_S, IROW1_S
      REAL(fpk) ::         KS
      CHARACTER (LEN=3) :: CI

!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ START OF BLOCK  @@@@@@@@@@@
      INTEGER   ::         O2, J, OM
      REAL(fpk) ::         INTEGRAND ( MAX_SLEAVEWFS, MAX_STREAMS )
      REAL(fpk) ::         SUM_R, SUM_CR, REFLEC, S_REFLEC, REFL_ATTN
      REAL(fpk) ::         H1, H2, NXR, PXR, NXR1, NXR2, PXR1
!  @@@@@@@@@@@@@@@@ END OF BLOCK  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@

!  Initialise status

      STATUS  = LRRS_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Short-hand

      M   = FOURIER

      LB = N_SURFACE_WFS + 1
      UB = N_SURFACE_WFS + N_SLEAVE_WFS

!  Nothing to do if M > 0 and Isotropic

!mick fix 10/19/2015 - this line turned off
      !IF ( .not. DO_INCLUDE_DIRECTBEAM  ) RETURN

!mick fix 7/20/2016 - define linearized output along with the return
      !IF ( M.gt.0 .and. DO_SL_ISOTROPIC ) RETURN
      IF ( M.gt.0 .and. DO_SL_ISOTROPIC ) THEN
        LS_ELASTIC_F_UP ( LB:UB, 1:N_LOUTPUT, 1:N_USER_STREAMS ) = ZERO
        LS_ELASTIC_F_DN ( LB:UB, 1:N_LOUTPUT, 1:N_USER_STREAMS ) = ZERO

        LS_MEAN_ELASTIC_UP ( LB:UB, 1:N_LOUTPUT ) = ZERO
        LS_MEAN_ELASTIC_DN ( LB:UB, 1:N_LOUTPUT ) = ZERO

        LS_FLUX_ELASTIC_UP ( LB:UB, 1:N_LOUTPUT ) = ZERO
        LS_FLUX_ELASTIC_DN ( LB:UB, 1:N_LOUTPUT ) = ZERO

        RETURN
      ENDIF

!  BVP solution for perturbed integration constants
!  ------------------------------------------------

!  Compute the main column B' where AX = B'
!     Regular BVP Solution --->  NO TELESCOPING HERE

!  initialise. Vitally necessary

      DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NTOTAL
          COL2_WFSLEAVE(I,Q) = ZERO
        ENDDO
      ENDDO

!  Last layer : Add direct beam variation due to surface leaving

      N  = NLAYERS
      C0 = N*NSTR2 - NSTREAMS
      DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          IR = I-1
          CM = C0 + IR + 1
!mick fix 9/7/2012 - moved Q to 1st dimension
! Rob Fix @@@ 11 Sep 12, remove SURFACE_FACTOR
          !COL2_WFSLEAVE(CM,Q) = SURFACE_FACTOR * LSSL_DIRECT_BEAM(I,Q)
          !COL2_WFSLEAVE(CM,Q) = SURFACE_FACTOR * LSSL_DIRECT_BEAM(Q,I)
          COL2_WFSLEAVE(CM,Q) =  LSSL_DIRECT_BEAM(Q,I)
        ENDDO
      ENDDO

!  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO Q = 1, N_SLEAVE_WFS
          DO N = 1, NTOTAL
            SCOL2_WFSLEAVE(N,Q) = COL2_WFSLEAVE(N,Q)
          ENDDO
        ENDDO
      ENDIF

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SLEAVE_WFS, &
              BANDMAT2, MAX_BANDTOTAL, IPIVOT, &
              COL2_WFSLEAVE, MAX_TOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (multilayer) in LRRS_SLEAVE WFS'
          STATUS  = LRRS_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_SLEAVE and PCON_SLEAVE, all layer

        DO Q = 1, N_SLEAVE_WFS
          DO N = 1, NLAYERS
            C0 = (N-1)*NSTR2
            DO K = 1, NSTREAMS
              IROW = K
              IROW1 = IROW + NSTREAMS
              NCON_SLEAVE(Q,K,N) = COL2_WFSLEAVE(C0+IROW,Q)
              PCON_SLEAVE(Q,K,N) = COL2_WFSLEAVE(C0+IROW1,Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFSLEAVE

        CALL DGETRS &
           ( 'N', NTOTAL, N_SLEAVE_WFS, SMAT2, MAX_2_STREAMS, &
              SIPIVOT, SCOL2_WFSLEAVE, MAX_2_STREAMS, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (Reg. 1 layer) in LRRS_SLEAVE_WFS'
          STATUS  = LRRS_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_SLEAVE and PCON_SLEAVE, 1 layer

        DO Q = 1, N_SLEAVE_WFS
          N = 1
          DO K = 1, NSTREAMS
            IROW = K
            IROW1 = IROW + NSTREAMS
            NCON_SLEAVE(Q,K,N) = SCOL2_WFSLEAVE(IROW,Q)
            PCON_SLEAVE(Q,K,N) = SCOL2_WFSLEAVE(IROW1,Q)
          ENDDO
        ENDDO

!  end clause

      ENDIF

!  debug------------------------------------------
!        if ( do_debug_write.and.fourier_component.eq.0 ) then
!         DO N = 1, NLAYERS
!          DO K = 1, NSTREAMS
!           write(86,'(3i2,1p6e13.5)')FOURIER,N,K,
!     &                LCON(K,N), MCON(K,N),
!     &                NCON_SLEAVE(1,K,N),PCON_SLEAVE(1,K,N)
!          ENDDO
!         ENDDO
!        ENDIF

!  Get the Post-processed weighting functions
!  ==========================================

!  Upwelling weighting functions
!  -----------------------------

      IF ( DO_UPWELLING ) THEN

!  Derivative of BOA source function

!        KS = SURFACE_FACTOR. Bug rob fix 9/9/14. Should be 1.0 always

        KS = one

!mick fix 9/11/2012 - modified IF condition and added ELSE
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!          IF ( (DO_SL_ISOTROPIC .AND. M.EQ.0) .OR. .NOT.DO_SL_ISOTROPIC ) THEN
        IF ( DO_INCLUDE_DIRECTBEAM .AND. &
             ((DO_SL_ISOTROPIC .AND. M.EQ.0) .OR. .NOT.DO_SL_ISOTROPIC) ) THEN
!mick fix 9/7/2012 - moved Q to 1st dimension
          DO Q = 1, N_SLEAVE_WFS
            DO UM = 1, N_USER_STREAMS
              !LSSL_BOA_SOURCE(Q,UM) = KS * LSSL_USER_DIRECT_BEAM(UM,Q)
              LSSL_BOA_SOURCE(Q,UM) = KS * LSSL_USER_DIRECT_BEAM(Q,UM)
            ENDDO
          ENDDO
        ELSE
          LSSL_BOA_SOURCE = ZERO
        ENDIF


!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ START OF BLOCK  @@@@@@@@@@@

!  Diffuse Term: Contribution due to derivatives of BVP constants
!  -------------------------------------------------------------

!  First compute derivative of downward intensity Integrand at stream an
!        .. reflectance integrand  = a(j).x(j).dI_DOWN(-j)/dS

!  start loops

        N = NLAYERS
        DO Q = 1, N_SLEAVE_WFS
         DO I = 1, NSTREAMS

!  Real homogeneous solutions

           SUM_R = ZERO
           DO K = 1, NSTREAMS
            NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I,K,N)
            PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I,K,N)
            SUM_R = SUM_R + NXR*T_DELT_EIGEN(K,N) + PXR
           ENDDO

!  Final result

           INTEGRAND(Q,I) = QUAD_STRMWGT(I) * SUM_R

!  end loops

         ENDDO
        ENDDO

!  integrated reflectance term
!  ---------------------------

!  Lambertian case, same for all user-streams

        IF ( .NOT. DO_BRDF_SURFACE ) THEN
          IF ( FOURIER.EQ.0 ) THEN
            DO Q = 1, N_SLEAVE_WFS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + INTEGRAND(Q,J)
              ENDDO
              REFLEC = SURFACE_FACTOR * REFLEC * ALBEDOS_RANKED(Raman_Idx)
              DO UM = 1, N_USER_STREAMS
                LSSL_BOA_SOURCE(Q,UM) = LSSL_BOA_SOURCE(Q,UM) + REFLEC
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  BRDF case

        IF ( DO_BRDF_SURFACE ) THEN
          DO Q = 1, N_SLEAVE_WFS
            DO UM = 1, N_USER_STREAMS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                S_REFLEC = INTEGRAND(Q,J) * USER_BRDF_F(M,UM,J,Brdf_IDX)
                REFLEC   = REFLEC + S_REFLEC
              ENDDO
              LSSL_BOA_SOURCE(Q,UM) = LSSL_BOA_SOURCE(Q,UM) + &
                                      REFLEC * SURFACE_FACTOR
            ENDDO
          ENDDO
        ENDIF

!  @@@@@@@@@@@@@@@@ END OF BLOCK    @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@

!  Upwelling Surface WF contribution for SLEAVE

        CALL LSSL_UPUSER_SURFACEWF ( &
          N_SLEAVE_WFS, N_SURFACE_WFS,                     &
          NLAYERS, NSTREAMS, N_LOUTPUT, N_USER_STREAMS,    &
          LEVEL_MASK_UP, FLUX_MULTIPLIER, LSSL_BOA_SOURCE, &
          T_DELT_USERM, U_XPOS, U_XNEG,                    &
          HMULT_1, HMULT_2, NCON_SLEAVE, PCON_SLEAVE,      &
          SURFACEWF_F )

        LS_ELASTIC_F_UP ( LB:UB, 1:N_LOUTPUT, 1:N_USER_STREAMS ) = &
            SURFACEWF_F ( LB:UB, 1:N_LOUTPUT, 1:N_USER_STREAMS, UPIDX )

!mick fix 9/7/2012 - ELSE added
      ELSE
        LS_ELASTIC_F_UP ( LB:UB, 1:N_LOUTPUT, 1:N_USER_STREAMS ) = ZERO
      ENDIF

!  Downwelling Albedo weighting functions
!  --------------------------------------

      IF ( DO_DNWELLING ) THEN
        CALL LSSL_DNUSER_SURFACEWF ( &
          N_SLEAVE_WFS, N_SURFACE_WFS,                & 
          NSTREAMS, N_LOUTPUT, N_USER_STREAMS,        &
          LEVEL_MASK_DN, FLUX_MULTIPLIER,             &
          T_DELT_USERM, U_XNEG, U_XPOS,               &
          HMULT_1, HMULT_2, NCON_SLEAVE, PCON_SLEAVE, &
          SURFACEWF_F )

        LS_ELASTIC_F_DN ( LB:UB, 1:N_LOUTPUT, 1:N_USER_STREAMS) = &
            SURFACEWF_F ( LB:UB, 1:N_LOUTPUT, 1:N_USER_STREAMS, DNIDX)

!mick fix 9/7/2012 - ELSE added
      ELSE
        LS_ELASTIC_F_DN ( LB:UB, 1:N_LOUTPUT, 1:N_USER_STREAMS) = ZERO
      ENDIF

!  Mean value output
!  -----------------

      IF ( DO_INCLUDE_MVOUT ) THEN
        CALL LSSL_INTEGRATED_OUTPUT ( &
          DO_UPWELLING, DO_DNWELLING,   &
          N_SLEAVE_WFS, N_SURFACE_WFS,  &
          NSTREAMS, NLAYERS, N_LOUTPUT, &
          LEVEL_MASK_UP, LEVEL_MASK_DN, FLUX_MULTIPLIER,  &
          QUAD_WEIGHTS, QUAD_STRMWGT, T_DELT_EIGEN,       &
          SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
          MINT_SURFACEWF, FLUX_SURFACEWF )

!mick fix 7/20/2016 - Upwelling & downwelling IF blocks added
        IF ( DO_UPWELLING ) THEN
          LS_MEAN_ELASTIC_UP (LB:UB,1:N_LOUTPUT) = MINT_SURFACEWF (LB:UB,1:N_LOUTPUT,UPIDX)
          LS_FLUX_ELASTIC_UP (LB:UB,1:N_LOUTPUT) = FLUX_SURFACEWF (LB:UB,1:N_LOUTPUT,UPIDX)
        ELSE
          LS_MEAN_ELASTIC_UP (LB:UB,1:N_LOUTPUT) = ZERO
          LS_FLUX_ELASTIC_UP (LB:UB,1:N_LOUTPUT) = ZERO
        ENDIF

        IF ( DO_DNWELLING ) THEN
          LS_MEAN_ELASTIC_DN (LB:UB,1:N_LOUTPUT) = MINT_SURFACEWF (LB:UB,1:N_LOUTPUT,DNIDX)
          LS_FLUX_ELASTIC_DN (LB:UB,1:N_LOUTPUT) = FLUX_SURFACEWF (LB:UB,1:N_LOUTPUT,DNIDX)
        ELSE
          LS_MEAN_ELASTIC_DN (LB:UB,1:N_LOUTPUT) = ZERO
          LS_FLUX_ELASTIC_DN (LB:UB,1:N_LOUTPUT) = ZERO
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LRRS_LSSL_WFS

!

      SUBROUTINE LSSL_UPUSER_SURFACEWF ( &
        N_SLEAVE_WFS, N_SURFACE_WFS,                     &
        NLAYERS, NSTREAMS, N_LOUTPUT, N_USER_STREAMS,    &
        LEVEL_MASK_UP, FLUX_MULTIPLIER, LSSL_BOA_SOURCE, &
        T_DELT_USERM, U_XPOS, U_XNEG,                    &
        HMULT_1, HMULT_2, NCON_SLEAVE, PCON_SLEAVE,      &
        SURFACEWF_F )

      USE LRRS_PARS_m, ONLY : fpk, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_LOUTPUT, &
                              MAX_SLEAVEWFS, MAX_SURFACEWFS, MAX_DIRECTIONS, ZERO, PI4, UPIDX

      IMPLICIT NONE

!  Input control

      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTREAMS

      INTEGER, INTENT (IN) ::          N_LOUTPUT
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

      INTEGER, INTENT (IN) ::          LEVEL_MASK_UP  ( MAX_LOUTPUT )

!  Solution inputs

      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER
      REAL(fpk), INTENT (IN) ::        LSSL_BOA_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

      REAL(fpk), INTENT (IN) ::        T_DELT_USERM &
          ( MAX_LAYERS, MAX_USER_STREAMS )

      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        HMULT_1 &
          ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_2 &
          ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )

!  output

      REAL(fpk), INTENT (INOUT) ::        SURFACEWF_F &
         ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, Q1, UT

      REAL(fpk) ::        HELP

      REAL(fpk) ::        LSSL_CUMUL_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) ::        LSSL_LAYER_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) ::        LSSL_FINAL_SOURCE


!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      DO UTA = 1, N_LOUTPUT
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_SURFACE_WFS
            SURFACEWF_F(Q1,UTA,UM,UPIDX) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to the BOA sum

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, N_SLEAVE_WFS
          LSSL_CUMUL_SOURCE(Q,UM) = LSSL_BOA_SOURCE(Q,UM)
        ENDDO
      ENDDO

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_LOUTPUT, 1, -1

!  Layer index for given optical depth

        NLEVEL = LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        NUT = NLEVEL + 1
        DO N = NSTART, NUT, -1
          CALL LSSL_WHOLELAYER_STERM_UP ( &
              N_SLEAVE_WFS, NSTREAMS, N, N_USER_STREAMS, &
              NCON_SLEAVE, PCON_SLEAVE,         &
              U_XPOS, U_XNEG, HMULT_1, HMULT_2, &
              LSSL_LAYER_SOURCE )

          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SLEAVE_WFS
              LSSL_CUMUL_SOURCE(Q,UM) = LSSL_LAYER_SOURCE(Q,UM) &
                 + T_DELT_USERM(N,UM) * LSSL_CUMUL_SOURCE(Q,UM)
            ENDDO
          ENDDO
        ENDDO

!  Ongrid output
!  -------------

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_SURFACE_WFS
!mick fix 7/20/2016 - normalize by PI4
            !SURFACEWF_F(Q1,UTA,UM,UPIDX) = &
            !           FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,UM)
            HELP = FLUX_MULTIPLIER * PI4
            SURFACEWF_F(Q1,UTA,UM,UPIDX) = &
                       HELP * LSSL_CUMUL_SOURCE(Q,UM)
          ENDDO
        ENDDO

!  Check for updating the recursion

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_UPUSER_SURFACEWF

!

      SUBROUTINE LSSL_DNUSER_SURFACEWF ( &
        N_SLEAVE_WFS, N_SURFACE_WFS,                &
        NSTREAMS, N_LOUTPUT, N_USER_STREAMS,        &
        LEVEL_MASK_DN, FLUX_MULTIPLIER,             &
        T_DELT_USERM, U_XNEG, U_XPOS,               &
        HMULT_1, HMULT_2, NCON_SLEAVE, PCON_SLEAVE, &
        SURFACEWF_F )

      USE LRRS_PARS_m, ONLY : fpk, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_LOUTPUT, &
                              MAX_SLEAVEWFS, MAX_SURFACEWFS, MAX_DIRECTIONS, ZERO, PI4, DNIDX

      IMPLICIT NONE

!  Input control

      INTEGER, INTENT (IN) ::          N_SURFACE_WFS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTREAMS

      INTEGER, INTENT (IN) ::          N_LOUTPUT
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

      INTEGER, INTENT (IN) ::          LEVEL_MASK_DN  ( MAX_LOUTPUT )

!  Solution input variables

      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER

      REAL(fpk), INTENT (IN) ::        T_DELT_USERM &
          ( MAX_LAYERS, MAX_USER_STREAMS )

      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        HMULT_1 &
          ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_2 &
          ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )

!  output

      REAL(fpk), INTENT (INOUT) ::        SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, Q1, UT

      REAL(fpk) ::        HELP

      REAL(fpk) ::        LSSL_CUMUL_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) ::        LSSL_LAYER_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) ::        LSSL_TOA_SOURCE   ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) ::        LSSL_FINAL_SOURCE

!  Zero all Fourier component output

      DO UTA = 1, N_LOUTPUT
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_SURFACE_WFS
            SURFACEWF_F(Q1,UTA,UM,DNIDX) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Initialize post-processing recursion
!  ====================================

!  Get the linearized TOA source terms

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, N_SLEAVE_WFS
          LSSL_TOA_SOURCE(Q,UM)   = ZERO
          LSSL_CUMUL_SOURCE(Q,UM) = LSSL_TOA_SOURCE(Q,UM)
        ENDDO
      ENDDO

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_LOUTPUT

!  Layer index for given optical depth

        NLEVEL = LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        NUT = NLEVEL
        DO N = NSTART, NUT

          CALL LSSL_WHOLELAYER_STERM_DN ( &
              N_SLEAVE_WFS, NSTREAMS, N, N_USER_STREAMS, &
              NCON_SLEAVE, PCON_SLEAVE,                  &
              U_XNEG, U_XPOS, HMULT_1, HMULT_2,          &
              LSSL_LAYER_SOURCE )

          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SLEAVE_WFS
              LSSL_CUMUL_SOURCE(Q,UM) = LSSL_LAYER_SOURCE(Q,UM) &
                + T_DELT_USERM(N,UM) * LSSL_CUMUL_SOURCE(Q,UM)
            ENDDO
          ENDDO

        ENDDO

!  Ongrid output
!  -------------

!  User-defined stream output, just set to the cumulative source term

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_SURFACE_WFS
!mick fix 7/20/2016 - normalize by PI4
            !SURFACEWF_F(Q1,UTA,UM,DNIDX) = &
            !  FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,UM)
            HELP = FLUX_MULTIPLIER * PI4
            SURFACEWF_F(Q1,UTA,UM,DNIDX) = &
              HELP * LSSL_CUMUL_SOURCE(Q,UM)
           ENDDO
        ENDDO


!  Check for updating the recursion

        IF ( NUT.NE. NUT_PREV ) NSTART = NUT + 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_DNUSER_SURFACEWF

!

      SUBROUTINE LSSL_INTEGRATED_OUTPUT ( &
        DO_UPWELLING, DO_DNWELLING,                     &
        N_SLEAVE_WFS, N_SURFACE_WFS,                    &
        NSTREAMS, NLAYERS, N_LOUTPUT,                   &
        LEVEL_MASK_UP, LEVEL_MASK_DN, FLUX_MULTIPLIER,  & 
        QUAD_WEIGHTS, QUAD_STRMWGT, T_DELT_EIGEN,       &
        SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
        MINT_SURFACEWF, FLUX_SURFACEWF )

!  Quadrature output at ongrid optical depths
!  ( Required if mean-value calculations are to be done)

      USE LRRS_PARS_m, ONLY : fpk, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_LOUTPUT, MAX_SLEAVEWFS, &
                              MAX_SURFACEWFS, MAX_DIRECTIONS, ZERO, PI2, HALF, UPIDX, DNIDX

      IMPLICIT NONE

!  INPUT
!  =====

!  Control

      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING

!  Basic numbers

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS, N_SURFACE_WFS
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  Directional and Level output controls

      INTEGER, INTENT (IN) ::          N_LOUTPUT
      INTEGER, INTENT (IN) ::          LEVEL_MASK_UP  ( MAX_LOUTPUT )
      INTEGER, INTENT (IN) ::          LEVEL_MASK_DN  ( MAX_LOUTPUT )

!  stream and multiplier input

      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER
      REAL(fpk), INTENT (IN) ::        QUAD_WEIGHTS ( MAX_STREAMS )
      REAL(fpk), INTENT (IN) ::        QUAD_STRMWGT ( MAX_STREAMS )

!  Solution variables

      REAL(fpk), INTENT (IN) ::        T_DELT_EIGEN ( MAX_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        SOLA_XPOS &
          ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        SOLB_XNEG &
          ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )

!  output
!  ======

      REAL(fpk), INTENT (INOUT) ::        MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_DIRECTIONS )
      REAL(fpk), INTENT (INOUT) ::        FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_DIRECTIONS )

!  Local variables
!  ===============

!  Directional controls

      INTEGER ::          N_DIRECTIONS
      INTEGER ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Local quadrature output

      REAL(fpk) ::        QSLEAVEWF_F &
          ( MAX_SLEAVEWFS, MAX_LOUTPUT, MAX_STREAMS )

!  Help variables

      INTEGER ::          I, IDIR, WDIR, UTA, UT, N, NLEVEL, O1, Q, Q1
      REAL(fpk) ::        SMI, SFX

!  Define directional controls

      N_DIRECTIONS = 0
      IF (DO_UPWELLING) THEN
        N_DIRECTIONS = N_DIRECTIONS + 1
        WHICH_DIRECTIONS(N_DIRECTIONS) = 1
      ENDIF
      IF (DO_DNWELLING) THEN
        N_DIRECTIONS = N_DIRECTIONS + 1
        WHICH_DIRECTIONS(N_DIRECTIONS) = 2
      ENDIF

!  Direction loop

      DO IDIR = 1, N_DIRECTIONS
        WDIR = WHICH_DIRECTIONS(IDIR)

!  Upwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. UPIDX ) THEN
          DO UTA = 1, N_LOUTPUT
            NLEVEL = LEVEL_MASK_UP(UTA)
            CALL LSSL_QUADWF_LEVEL_UP ( &
               UTA, NLEVEL, N_SLEAVE_WFS, NSTREAMS, NLAYERS,   &
               FLUX_MULTIPLIER, T_DELT_EIGEN,                  &
               SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
               QSLEAVEWF_F )
          ENDDO
        ENDIF

!  Downwelling Jacobian output at Quadrature angles

        IF (WDIR .EQ. DNIDX) THEN
          DO UTA = 1, N_LOUTPUT
            NLEVEL = LEVEL_MASK_DN(UTA)
            CALL LSSL_QUADWF_LEVEL_DN  (  &
               UTA, NLEVEL, N_SLEAVE_WFS, NSTREAMS,            &
               FLUX_MULTIPLIER, T_DELT_EIGEN,                  &
               SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
               QSLEAVEWF_F )
          ENDDO
        ENDIF

!  Mean Intensity and Flux output (diffuse term only). If flagged

        DO UTA = 1, N_LOUTPUT
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_SURFACE_WFS
            SMI = ZERO
            SFX = ZERO
            DO I = 1, NSTREAMS
              SMI = SMI + QUAD_WEIGHTS(I) * QSLEAVEWF_F(Q,UTA,I)
              SFX = SFX + QUAD_STRMWGT(I) * QSLEAVEWF_F(Q,UTA,I)
            ENDDO
            MINT_SURFACEWF(Q1,UTA,WDIR) = SMI * HALF
            FLUX_SURFACEWF(Q1,UTA,WDIR) = SFX * PI2
          ENDDO
        ENDDO

!  End directions loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_INTEGRATED_OUTPUT

!

      SUBROUTINE LSSL_WHOLELAYER_STERM_UP ( &
        N_SLEAVE_WFS, NSTREAMS, GIVEN_LAYER, N_USER_STREAMS, &
        NCON_SLEAVE, PCON_SLEAVE,                            &
        U_XPOS, U_XNEG, HMULT_1, HMULT_2,                    &
        LSSL_LAYERSOURCE )

      USE LRRS_PARS_m, ONLY : fpk, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_SLEAVEWFS, ZERO

      IMPLICIT NONE

!  Control inputs

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          GIVEN_LAYER
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  Solution inputs

      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_1 &
          ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_2 &
          ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )

!  outputs

      REAL(fpk), INTENT (OUT) ::        LSSL_LAYERSOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2


!  local indices

      N   = GIVEN_LAYER

!  Homogeneous solutions
!  =====================

!  Loops over user angles, weighting functions and streams

      DO UM = 1, N_USER_STREAMS
       DO Q = 1, N_SLEAVE_WFS

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, NSTREAMS
            NUXR = NCON_SLEAVE(Q,K,N)*U_XPOS(UM,K,N)
            PUXR = PCON_SLEAVE(Q,K,N)*U_XNEG(UM,K,N)
            H1 =  NUXR * HMULT_2(K,UM,N)
            H2 =  PUXR * HMULT_1(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  homogeneous contribution

          LSSL_LAYERSOURCE(Q,UM) = SHOM_R

!  End loops over Q and UM

       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_WHOLELAYER_STERM_UP

!

      SUBROUTINE LSSL_WHOLELAYER_STERM_DN ( &
        N_SLEAVE_WFS, NSTREAMS, GIVEN_LAYER, N_USER_STREAMS, &
        NCON_SLEAVE, PCON_SLEAVE,                            &
        U_XNEG, U_XPOS, HMULT_1, HMULT_2,                    &
        LSSL_LAYERSOURCE )

      USE LRRS_PARS_m, ONLY : fpk, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_SLEAVEWFS, ZERO

      IMPLICIT NONE

!  Control inputs

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          GIVEN_LAYER
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  Solution inputs

      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_1 &
          ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_2 &
          ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )

!  output

      REAL(fpk), INTENT (OUT) ::       LSSL_LAYERSOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  local indices

      N   = GIVEN_LAYER

!  Homogeneous solutions
!  =====================

!  Loop over user angles, weighting functions and STokes

      DO UM = 1, N_USER_STREAMS
       DO Q = 1, N_SLEAVE_WFS

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, NSTREAMS
            NUXR = NCON_SLEAVE(Q,K,N)*U_XNEG(UM,K,N)
            PUXR = PCON_SLEAVE(Q,K,N)*U_XPOS(UM,K,N)
            H1 =  NUXR * HMULT_1(K,UM,N)
            H2 =  PUXR * HMULT_2(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  homogeneous contribution

          LSSL_LAYERSOURCE(Q,UM) = SHOM_R

!  End loops over Q and UM

       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_WHOLELAYER_STERM_DN

!

      SUBROUTINE LSSL_QUADWF_LEVEL_UP ( &
        UTA, NL, N_SLEAVE_WFS, NSTREAMS, NLAYERS,       &
        FLUX_MULTIPLIER, T_DELT_EIGEN,                  &
        SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
        QSURFACEWF_F )

!  Upwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

      USE LRRS_PARS_m, ONLY : fpk, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_SLEAVEWFS, MAX_LOUTPUT, ZERO

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, NL, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTREAMS, NLAYERS
      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER

      REAL(fpk), INTENT (IN) ::        T_DELT_EIGEN ( MAX_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        SOLA_XPOS &
          ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        SOLB_XNEG &
          ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )

!  output

      REAL(fpk), INTENT (INOUT) ::     QSURFACEWF_F &
          ( MAX_SLEAVEWFS, MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, I, I1, O1, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

      N = NL + 1

!  For the lowest level
!  ====================

      IF ( NL .EQ. NLAYERS ) THEN

!  parameter loop, Stokes and streams loops

       DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, NSTREAMS
              NXR = NCON_SLEAVE(Q,K,NL) * SOLA_XPOS(I1,K,NL)
              PXR = PCON_SLEAVE(Q,K,NL) * SOLB_XNEG(I1,K,NL)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,NL) + PXR
            ENDDO

!  collect solution

            QSURFACEWF_F(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R

!  Finish loops

        ENDDO
       ENDDO

      ENDIF

!  For other levels in the atmosphere
!  ==================================

      IF ( NL .NE. NLAYERS ) THEN

!  parameter loop, Stokes and streams loops

       DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, NSTREAMS
              NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I1,K,N)
              PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I1,K,N)
              SHOM_R = SHOM_R + PXR*T_DELT_EIGEN(K,N) + NXR
            ENDDO

!  collect solution

           QSURFACEWF_F(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R

!  Finish loops

        ENDDO
       ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LSSL_QUADWF_LEVEL_UP

!

      SUBROUTINE LSSL_QUADWF_LEVEL_DN ( &
        UTA, NL, N_SLEAVE_WFS, NSTREAMS,                &
        FLUX_MULTIPLIER, T_DELT_EIGEN,                  &
        SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
        QSURFACEWF_F )

      USE LRRS_PARS_m, ONLY : fpk, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_SLEAVEWFS, MAX_LOUTPUT, ZERO

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, NL, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTREAMS
      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER

      REAL(fpk), INTENT (IN) ::        T_DELT_EIGEN ( MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        SOLA_XPOS &
          ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        SOLB_XNEG &
          ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAX_STREAMS, MAX_LAYERS )

!  output

      REAL(fpk), INTENT (INOUT) ::     QSURFACEWF_F &
          ( MAX_SLEAVEWFS, MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, I, O1, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NXR, NXR1, NXR2, PXR, PXR1

!  Downwelling weighting functions at TOA ( or N = 0 ) are zero
!  ============================================================

      IF ( NL .EQ. 0 ) THEN
        DO Q = 1, N_SLEAVE_WFS
          DO I = 1, NSTREAMS
            QSURFACEWF_F(Q,UTA,I) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Downwelling weighting functions at other levels
!  ===============================================

!  Offset

      N = NL

!  parameter loop,  Stokes and streams loops

      DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, NSTREAMS
              NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I,K,N)
              PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I,K,N)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,N) + PXR
            ENDDO

!  collect solution

            QSURFACEWF_F(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R

!  Finish loops

        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_QUADWF_LEVEL_DN

!  Finish module

      END MODULE lrrs_L_sleave_m

