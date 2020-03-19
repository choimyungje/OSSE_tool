
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
! #   MASTER_MODULE :    SOURCES_MASTER_1                       #
! #                        (1 = Whole-layer)                    #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Use of Taylor-series control, TAYLOR_ORDER index, passed to PostProcessing
!    (2) Bookkeeping improvements (use of "Only", clearer I/O specifications)

      MODULE lrrs_sources_master_1_m

!      USE LRRS_PARS_m, Only : SDU

      PRIVATE
      PUBLIC :: SOURCES_MASTER_1

      CONTAINS

      SUBROUTINE SOURCES_MASTER_1 &
       ( FOURIER, ACTUAL_POINT, DO_ENERGY_BALANCING,                  & ! Inputs
         DO_RRS_OVERALL, DO_DIRECTRRS_ONLY, DO_MSMODE_LRRS,           & ! Inputs
         DO_BIN_REALIZATION, DO_MONO_REALIZATION, DO_SSCORR_OUTGOING, & ! Inputs
         DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                 & ! Inputs
         N_RRSBINS, BINMAP, NPOINTS_MONO, OFFSET_INNER, W_EXCIT,      & ! Inputs
         NLAYERS, NSTREAMS, NMOMENTS, N_USER_STREAMS,                 & ! Inputs
         QUAD_WEIGHTS, USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER,      & ! Inputs
         STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, FLUXES_RANKED,       & ! Inputs
         DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, OMEGAMOMS_CABANNES,    & ! Inputs
         OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN, OMEGAMOMS_RRSGAIN,      & ! Inputs
         BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT,                   & ! Inputs
         BEAM_DTRANS, SAVE_TRANS_USERM, PLM_PXI, PLM_MXI,             & ! Inputs
         PLM_PXUI,    PLM_MXUI,    PLM_00_PXI,  PLM_00_MXI,           & ! Inputs
         PLM_00_PXUI, PLM_00_MXUI, PLM_WT_PXI,  PLM_WT_MXI,           & ! Inputs
         L0_KEIGEN,     L0_KTRANS,     L0_WPARTIC,                    & ! Inputs
         L0_WHOM_XPOS,  L0_WHOM_LCON,  L0_WHOM_MCON,                  & ! Inputs
         XPOS, KEIGEN, KTRANS, U_XPOS, U_XNEG,                        & ! Inputs
         EIGENNORM_SAVE, HMULT_1, HMULT_2,                            & ! Inputs
         RAMAN_WUPPER, RAMAN_WLOWER, WLAYER_PIST_UP, WLAYER_PIST_DN,  & ! Output
         FAIL, MESSAGE )                                                ! Output

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, MAX_2_STREAMS, MAX_STREAMS, MAX_USER_STREAMS, &
                              MAX_MOMENTS, MAX_LAYERS, MAX_BINS, MAX_POINTS

      USE LRRS_RTSOLUTIONS_1_m   , Only : GREENFUNC_SOLUTION_1
      USE LRRS_POSTPROCESSING_1_m, Only : POSTPROCESSING_1

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Fourier point, wavelength point

      INTEGER  , INTENT(IN) :: FOURIER
      INTEGER  , INTENT(IN) :: ACTUAL_POINT

!  Flags
!  -----

!  Solution method (options are mutually exclusive)

      LOGICAL  , INTENT(IN) :: DO_ENERGY_BALANCING

!  overall filling flags (directrrs is the solar term only)

      LOGICAL  , INTENT(IN) :: DO_RRS_OVERALL
      LOGICAL  , INTENT(IN) :: DO_DIRECTRRS_ONLY
      LOGICAL  , INTENT(IN) :: DO_MSMODE_LRRS

!  Binning or monochromatic realization
!    (Options are mutually exclusive)

      LOGICAL  , INTENT(IN) :: DO_BIN_REALIZATION
      LOGICAL  , INTENT(IN) :: DO_MONO_REALIZATION

!  SS correction flag; Now incorporates the DB term

      LOGICAL  , INTENT(IN) :: DO_SSCORR_OUTGOING

!  Upwelling and downwelling flags

      LOGICAL  , INTENT(IN) :: DO_UPWELLING
      LOGICAL  , INTENT(IN) :: DO_DNWELLING

!  output for results at user angles (post processing)

      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS

!  Points
!  ------

!  offset for the inner range
!    Depends on the choice of solar spectrum

      INTEGER  , INTENT(IN) :: OFFSET_INNER

!  Bin Mapping quantities

      INTEGER  , INTENT(IN) :: N_RRSBINS  ( MAX_POINTS )
      INTEGER  , INTENT(IN) :: BINMAP ( MAX_BINS, MAX_POINTS )

!  Number of monochromatic points to be calculated = 234

      INTEGER  , INTENT(IN) :: NPOINTS_MONO

!  W_EXCIT is the index of the excitation wavelength in the Rank.
!     ( Monochromatic only )

      INTEGER  , INTENT(IN) :: W_EXCIT

!  Numbers and bookkeeping
!  -----------------------

!  Number of layers

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of discrete ordinate streams

      INTEGER  , INTENT(IN) :: NSTREAMS

!  Number of moments (elastic scattering)

      INTEGER  , INTENT(IN) :: NMOMENTS

!  User stream variables

      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  Taylor-series control. Updated, Version 2.5, 9/10/15

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  Layer masks for doing integrated source terms

      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP ( MAX_LAYERS )
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN ( MAX_LAYERS )

!  Quadrature weights

      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS  ( MAX_STREAMS )

!  Optical properties
!  ------------------

!  Depends on the choice of solar spectrum

      REAL(FPK), INTENT(IN) :: FLUXES_RANKED   ( MAX_POINTS )

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT   ( MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived Cabannes scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_CABANNES &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived rotational raman scattering (Internal variable). Loss term
!    Single scatter albedo * phase function moments
!    Only required for the Energy-balance approximation.

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSLOSS &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived rotational Raman scattering (binning)

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSBIN &
               ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Derived rotational Raman scattering (Monochromatic). Gain terms
!    Single scatter albedo * phase function moments

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSGAIN &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Local setup variables
!  ---------------------

!  Solar beam transmittances, average secant factors

      INTEGER  , INTENT(IN) :: BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(IN) :: SAVE_TRANS_USERM &
            (  MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Legendre polynomial variables
!  -----------------------------

!  Legendre functions at discrete ordinates

      REAL(FPK), INTENT(IN) :: PLM_PXI    ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_MXI    ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Product functions, discrete ordinates and solar stream

      REAL(FPK), INTENT(IN) :: PLM_00_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_00_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Product of Legendre functions with quadrature weights

      REAL(FPK), INTENT(IN) :: PLM_WT_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_WT_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Legendre functions at user defined polar streams

      REAL(FPK), INTENT(IN) :: PLM_PXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_MXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  Product functions, User streams and solar stream

      REAL(FPK), INTENT(IN) :: PLM_00_PXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_00_MXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  local homogeneous solutions
!  ---------------------------

!  (Positive) Eigenvalues

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!   Whole layer transmittance factors for +/- eigenvalues

      REAL(FPK), INTENT(IN) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  user-stream homogeneous solution vectors

      REAL(FPK), INTENT(IN) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Green's function Eigennorms

      REAL(FPK), INTENT(IN) :: EIGENNORM_SAVE ( MAX_STREAMS, MAX_LAYERS )

!  Whole layer multipliers (homogeneous)

      REAL(FPK), INTENT(IN) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Elastic discrete ordinate solutions at shifted wavelengths
!  ----------------------------------------------------------

!  Eigenvalues

      REAL(FPK), INTENT(IN) :: L0_KEIGEN &
              ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Eigentransmittances

      REAL(FPK), INTENT(IN) :: L0_KTRANS &
              ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Particular solution vectors

      REAL(FPK), INTENT(IN) :: L0_WPARTIC &
            ( MAX_2_STREAMS, MAX_LAYERS, MAX_POINTS )

!  homogeneous solution vectors and integration constants

      REAL(FPK), INTENT(IN) :: L0_WHOM_XPOS &
            ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L0_WHOM_LCON &
            ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L0_WHOM_MCON &
            ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  OUTPUT ARGUMENTS
!  ================

!  output arguments from Beam solution modules
!  -------------------------------------------

!  Particular solutions at upper and lower boundaries
!   (as evaluated for the inelastic field by Green's function method)

      REAL(FPK), INTENT(OUT) :: RAMAN_WUPPER  ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(OUT) :: RAMAN_WLOWER  ( MAX_2_STREAMS, MAX_LAYERS )

!  post processed up and down Raman solutions

      REAL(FPK), INTENT(OUT) :: WLAYER_PIST_UP ( MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: WLAYER_PIST_DN ( MAX_LAYERS, MAX_USER_STREAMS )

!  module exception handling
!  -------------------------

      LOGICAL            , INTENT(OUT) :: FAIL
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE

!  Local variables
!  ===============

!  Source function control

      LOGICAL :: CHANGE_EXPONENT
      LOGICAL :: FLIPPER
      INTEGER :: PICUTOFF
      LOGICAL :: DO_MULTIPLIER

!  Source function transmittances, average secant factors

      REAL(FPK) :: ASOURCE
      REAL(FPK) :: INITRANS
      REAL(FPK) :: DELTRANS

!  Source function Q-vectors

      REAL(FPK) :: QSOURCE_QPOS  ( MAX_STREAMS )
      REAL(FPK) :: QSOURCE_QNEG  ( MAX_STREAMS )
      REAL(FPK) :: QSOURCE_UPOS1 ( MAX_USER_STREAMS )
      REAL(FPK) :: QSOURCE_UNEG1 ( MAX_USER_STREAMS )

!  Saved local multipliers

      REAL(FPK) :: WSU_SUSAV ( MAX_STREAMS, MAX_USER_STREAMS )
      REAL(FPK) :: WSU_SDSAV ( MAX_STREAMS, MAX_USER_STREAMS )
      REAL(FPK) :: WSD_SUSAV ( MAX_STREAMS, MAX_USER_STREAMS )
      REAL(FPK) :: WSD_SDSAV ( MAX_STREAMS, MAX_USER_STREAMS )

!  Green's function multipliers

      REAL(FPK) :: MULT_M  ( MAX_STREAMS )
      REAL(FPK) :: MULT_P  ( MAX_STREAMS )

!  Source function Beam multipliers, whole layer

      REAL(FPK) :: EMULT_UP ( MAX_USER_STREAMS )
      REAL(FPK) :: EMULT_DN ( MAX_USER_STREAMS )

!  Green's function arrays

      REAL(FPK) :: ATERM_SAVE ( MAX_STREAMS )
      REAL(FPK) :: BTERM_SAVE ( MAX_STREAMS )

!  Local optical depths

      REAL(FPK) :: DELTAUS

!  Trans UserM

      REAL(FPK) :: TRANS_USERM   ( MAX_USER_STREAMS )

!  Local arrays of omega-moments (IOPs)

      REAL(FPK) :: OMEGA_PHASMOMS     ( MAX_LAYERS, 0: MAX_MOMENTS )
      REAL(FPK) :: OMEGAMOMS_RRSLOCAL ( 0: MAX_MOMENTS )

!  Local arrays of optical depth scaling factors

      REAL(FPK) :: DELTA_SCALING

!  Elastic only flag

      LOGICAL   :: DO_ELASTIC_ONLY

!  Flux comparisons

      REAL(FPK) :: FM_SHIFTED
      REAL(FPK) :: FM_SHIFTED_DIRECT
      REAL(FPK) :: FM_SHIFTED_DIFFUSE
      REAL(FPK) :: FM_NONSHIFT
      REAL(FPK) :: FM_NONSHIFT_DIRECT
      REAL(FPK) :: FM_NONSHIFT_DIFFUSE

!  Dummy initial transmittances

      REAL(FPK) :: DUMTRANS

!  other local variables

      LOGICAL   :: DO_SSCORR_OVERALL
      INTEGER   :: S, P, Z, N, I, ID, ID1, SAVCUT, WW
      INTEGER   :: M, J, UI, L, APT, BPT, CPT
      INTEGER   :: NRRS_MOMENTS, NPOINTS_LOCAL

      REAL(FPK) :: HOLD_RRS(0:2), HOLD_MP(0:2)
      REAL(FPK) :: MC, LC, LP, LM, PHIFUNC, SP, SM

!  debug variables

      LOGICAL  , PARAMETER ::  DO_LOG_OUTPUT    = .FALSE.
!     LOGICAL  , PARAMETER ::  DO_DEBUG_SOURCES = .TRUE.
      LOGICAL  , PARAMETER ::  DO_DEBUG_SOURCES = .FALSE.

      INTEGER  , PARAMETER :: LOG_UNIT = 67

!  Factors (from an earlier version, previous to 2.3)

      INTEGER  , PARAMETER :: Y = -1
      REAL(FPK), PARAMETER :: DECFAC = 1.0D0

!  Start
!  =====

!  initialise output

      FAIL = .FALSE.
      MESSAGE = ' '

!  Overall single scatter correction flag
!    This must now be set for each PI solution
!     Bug detected and fixed, 23 July 2008
!     For SSCORR_OUTGOING, All SS contributions must be turned off
!      DO_SSCORR_OVERALL = DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING

!  Elastic-only flag:
!    RRS Terms not present for Fourier > 2

      DO_ELASTIC_ONLY = .NOT. DO_RRS_OVERALL .OR. &
                           ( DO_RRS_OVERALL .AND. FOURIER .GT. 2 )

!  flux multiplier;  which is correct factor?
!        Turns out to be DECFAC = 1.0

      IF ( DO_ENERGY_BALANCING ) THEN
        FM_NONSHIFT         = ONE
        FM_NONSHIFT_DIRECT  = FM_NONSHIFT
        FM_NONSHIFT_DIFFUSE = FM_NONSHIFT * DECFAC
      ENDIF

!  Fourier proxy

      M = FOURIER

!  Number of Legendre moments = 2

      IF ( NSTREAMS .EQ. 1 .AND. NLAYERS .EQ. 1 ) THEN
        NRRS_MOMENTS = 0
      ELSE
        NRRS_MOMENTS = 2
      ENDIF

!  calculation point
!  For the binning, we need the inner and outer loops

      WW = W_EXCIT
      IF ( DO_BIN_REALIZATION ) THEN
        APT = ACTUAL_POINT
        CPT = APT + OFFSET_INNER
      ELSE
        APT = ACTUAL_POINT
        CPT = W_EXCIT
      ENDIF

!  Number of points or bins, Raman gain terms

      IF ( DO_BIN_REALIZATION ) THEN
        NPOINTS_LOCAL = N_RRSBINS ( APT )
      ELSE
        NPOINTS_LOCAL = NPOINTS_MONO
      ENDIF

!  Dummy values of INITRANS

      DUMTRANS = ONE

!  local IOPS for Elastic scattering only

      IF ( DO_ELASTIC_ONLY ) THEN
        DO N = 1, NLAYERS
          DO L = 0, NMOMENTS
            OMEGA_PHASMOMS(N,L) = OMEGAMOMS_ELASTIC(N,L,CPT)
          ENDDO
        ENDDO
      ENDIF

!  local IOPS for INELASTIC scattering (Energy balancing)

      IF ( .NOT.DO_ELASTIC_ONLY ) THEN
        IF ( DO_ENERGY_BALANCING ) THEN
         DO N = 1, NLAYERS
           DO L = 0, NMOMENTS
            OMEGA_PHASMOMS(N,L) = OMEGAMOMS_ELASTIC(N,L,CPT)
           ENDDO
         ENDDO
        ENDIF
       ENDIF

!  local IOPS for INELASTIC scattering (Cabannes-Raman)

      IF ( .NOT.DO_ELASTIC_ONLY ) THEN
        IF ( .NOT.DO_ENERGY_BALANCING ) THEN
         DO N = 1, NLAYERS
           DO L = 0, NMOMENTS
            OMEGA_PHASMOMS(N,L) = OMEGAMOMS_CABANNES(N,L,APT)
           ENDDO
         ENDDO
        ENDIF
      ENDIF

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!         START LAYER LOOP
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      DO N = 1, NLAYERS

!  Local user-stream transmittances

        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            TRANS_USERM(UI) = SAVE_TRANS_USERM(UI,N,CPT)
          ENDDO
        ENDIF

!  Local deltas at CPT

        DELTAUS = DELTAU_VERT_INPUT(N,CPT)

!  Initialize solution and exponent numbers
!    P  for the exponent number
!    Z  for the regular solution

        P  = 0
        Z  = 0

!   #####################################
!   #####################################
!   #    E L A S T I C   T E R M        #
!   #####################################
!   #####################################

!  Source vector for the Elastic scattering term
!  =============================================

!  increase solution and exponent numbers by 1

        Z = Z + 1
        P = P + 1
        CHANGE_EXPONENT = .TRUE.
        FLIPPER         = .FALSE.
        PICUTOFF        = BEAM_PICUTOFF(CPT)
        DO_MULTIPLIER   = .TRUE.

!  Set single scatter overall flag, set freshly for each

        DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING
!        DO_SSCORR_OVERALL =
!     &     DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING

!  Step 1A. set transmittances and secants for this term
!           Need to re-set this each time.

        ASOURCE  = BEAM_AVSECANT(N,CPT)
        INITRANS = BEAM_ITRANS(N,CPT)
        DELTRANS = BEAM_DTRANS(N,CPT)

!  No step 2 here

!  Step 3A. Regular solution source vector, for quadrature streams

        DO I = 1, NSTREAMS
          SP = ZERO
          SM = ZERO
          DO L = FOURIER, NMOMENTS
            SP = SP + OMEGA_PHASMOMS(N,L) * PLM_00_PXI(I,L)
            SM = SM + OMEGA_PHASMOMS(N,L) * PLM_00_MXI(I,L)
          ENDDO
          QSOURCE_QPOS(I) = SP
          QSOURCE_QNEG(I) = SM
        ENDDO

!  Step 4A. Regular solution source vector for user streams

        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            SP = ZERO
            SM = ZERO
            DO L = FOURIER, NMOMENTS
              SP = SP + OMEGA_PHASMOMS(N,L) * PLM_00_PXUI(UI,L)
              SM = SM + OMEGA_PHASMOMS(N,L) * PLM_00_MXUI(UI,L)
            ENDDO
            QSOURCE_UPOS1(UI) = SP
            QSOURCE_UNEG1(UI) = SM
          ENDDO
        ENDIF

!  Whole layer Greens function Discrete ordinate solutions

        CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE, BTERM_SAVE, & ! Output
              MULT_M, MULT_P )                                      ! Output

!  Call to the post processing routine

        CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Finish layer if elastic only

        IF ( DO_ELASTIC_ONLY ) GO TO 5222

!  For the Cabannes-Raman method, Skip the next section

        IF ( .NOT. DO_ENERGY_BALANCING ) GO TO 7676

!   #####################################
!   #####################################
!   #    R R S   L O S S   T E R M S    #
!   #####################################
!   #####################################

!  Source vectors and PI solutions for non-shifted Solar beam term
!  ===============================================================

!  General Procedure in 3 Steps :
!       (1) set INITRANS, ENDTRANS, ASOURCE (transmittances and secants)
!       (2) source vectors for quadrature streams
!       (3) source vectors for user-angle streams

!  increase solution number by 1, P not changed
!    --- Do the exponent again in this linearized version (November 2008

        Z = Z + 1
!      CHANGE_EXPONENT = .FALSE.               ! commented out 11/11/08
        CHANGE_EXPONENT = .TRUE.
        DO_MULTIPLIER   = .FALSE.
        DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING

!  Step 1... set transmittances and secants (once only for FOURIER = 0 )
!     Already done,

!  No step 2

!  Step 3A. set source vector for quadrature streams

        DO I = 1, NSTREAMS
          SP = ZERO
          SM = ZERO
          DO L = FOURIER, NRRS_MOMENTS
            SP = SP + OMEGAMOMS_RRSLOSS(N,L,APT) * PLM_00_PXI(I,L)
            SM = SM + OMEGAMOMS_RRSLOSS(N,L,APT) * PLM_00_MXI(I,L)
          ENDDO
          QSOURCE_QPOS(I) = - FM_NONSHIFT_DIRECT * SP
          QSOURCE_QNEG(I) = - FM_NONSHIFT_DIRECT * SM
        ENDDO

!  Step 4A. set Regular Solution source vector for user streams

        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            SP = ZERO
            SM = ZERO
            DO L = FOURIER, NRRS_MOMENTS
              SP = SP + OMEGAMOMS_RRSLOSS(N,L,APT)*PLM_00_PXUI(UI,L)
              SM = SM + OMEGAMOMS_RRSLOSS(N,L,APT)*PLM_00_MXUI(UI,L)
            ENDDO
            QSOURCE_UPOS1(UI) = - FM_NONSHIFT_DIRECT * SP
            QSOURCE_UNEG1(UI) = - FM_NONSHIFT_DIRECT * SM
          ENDDO
        ENDIF

!  Whole layer Greens function Discrete ordinate solutions

        CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P )   ! OUTPUT

!  Call to the post processing routine

        CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Sources and PIs for the non-shifted Diffuse-field loss-terms
!  ============================================================

        IF ( .NOT. DO_DIRECTRRS_ONLY ) THEN

!  (b1) Due to diffuse-term particular integral of elastic solution
!       -----------------------------------------------------------

!  General Procedure in 4 Steps :
!       (1) set INITRANS, ENDTRANS, ASOURCE (transmittances and secants)
!       (2) integrate solution for each moment = HOLD(L)
!       (3) source vectors for quadrature streams
!       (4) source vectors for user-angle streams

!  Increase number of solutions by 1, P not changed (same exponent)
!       MUST CHANGE EXPONENT for LINEARIZED

          Z = Z + 1
!          CHANGE_EXPONENT = .FALSE.
          CHANGE_EXPONENT = .TRUE.
!          DO_MULTIPLIER   = .FALSE.
          DO_MULTIPLIER   = .TRUE.
!          DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING   ! Wrong
          DO_SSCORR_OVERALL = .FALSE.

!  Step 1. (set  transmittances and secants for this term)
!    This is the same as previous situation
!  Step 2.... (integrated solution contributions, array HOLD(L) )
!  Step 3.... (set source vector for quadrature streams)
!  Step 4.... (set source vector for user streams)

!  Step 2A. Holding terms, Regular solution

          DO L = FOURIER, NRRS_MOMENTS
            SP = ZERO
            SM = ZERO
            DO ID = 1, NSTREAMS
              ID1 = ID + NSTREAMS
              SP = SP + PLM_WT_PXI(ID,L) * L0_WPARTIC(ID1,N,WW)
              SM = SM + PLM_WT_MXI(ID,L) * L0_WPARTIC(ID,N,WW)
            ENDDO
            HOLD_MP(L)  = SP + SM
            HOLD_RRS(L) = OMEGAMOMS_RRSLOSS(N,L,APT) * HOLD_MP(L)
          ENDDO

!  Step 3A. Regular solution Source vector for quadrature streams

          DO I = 1, NSTREAMS
            SP = ZERO
            SM = ZERO
            DO L = FOURIER, NRRS_MOMENTS
              SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
              SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
            ENDDO
            QSOURCE_QPOS(I) = - FM_NONSHIFT_DIFFUSE * SP
            QSOURCE_QNEG(I) = - FM_NONSHIFT_DIFFUSE * SM
          ENDDO

!  Step 4A. Regular Solution Source vector for user streams

          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
                SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
                SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
              ENDDO
              QSOURCE_UPOS1(UI)  = - FM_NONSHIFT_DIFFUSE * SP
              QSOURCE_UNEG1(UI)  = - FM_NONSHIFT_DIFFUSE * SM
            ENDDO
          ENDIF

!  Whole layer Greens function Discrete ordinate solutions

          CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P )   ! OUTPUT

!  Call to the post processing routine

          CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  (b2) Due to homogeneous functions of non-shifted elastic solution
!       ------------------------------------------------------------

!  start stream loop

          DO J = 1, NSTREAMS

!  homogeneous DN solutions, increase Z and P by 1
!  ===============================================

           Z  = Z + 1
           P  = P + 1
           CHANGE_EXPONENT = .TRUE.
           DO_MULTIPLIER   = .TRUE.
           FLIPPER         = .FALSE.
!           DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING     ! Wrong
           DO_SSCORR_OVERALL = .FALSE.

!  Step 1A. Tranmsittances

           PICUTOFF = NLAYERS
           INITRANS = ONE
           DELTRANS = L0_KTRANS(J,N,W_EXCIT)
           ASOURCE  = L0_KEIGEN(J,N,W_EXCIT)

!  Step 2A. Regular Solution Holding terms

           LC = L0_WHOM_LCON(J,N,WW)
           DO L = FOURIER, NRRS_MOMENTS
            SP = ZERO
            SM = ZERO
            DO ID = 1, NSTREAMS
              ID1 = ID + NSTREAMS
              LP = PLM_WT_PXI(ID,L)
              LM = PLM_WT_MXI(ID,L)
              SP = SP + LC * LP * L0_WHOM_XPOS(ID1,J,N,WW)
              SM = SM + LC * LM * L0_WHOM_XPOS(ID,J,N,WW)
            ENDDO
            HOLD_MP(L)  = SP + SM
            HOLD_RRS(L) = OMEGAMOMS_RRSLOSS(N,L,APT) * HOLD_MP(L)
           ENDDO

!  Step 3A. Regular Solution Source vector for quadrature streams

           DO I = 1, NSTREAMS
            SP = ZERO
            SM = ZERO
            DO L = FOURIER, NRRS_MOMENTS
             SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
             SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
            ENDDO
            QSOURCE_QPOS(I) = - FM_NONSHIFT_DIFFUSE * SP
            QSOURCE_QNEG(I) = - FM_NONSHIFT_DIFFUSE * SM
           ENDDO

!  Step 4A. Regular Solution Source vector for user streams

           IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
             SP = ZERO
             SM = ZERO
             DO L = FOURIER, NRRS_MOMENTS
              SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
              SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
             ENDDO
             QSOURCE_UPOS1(UI)  = - FM_NONSHIFT_DIFFUSE * SP
             QSOURCE_UNEG1(UI)  = - FM_NONSHIFT_DIFFUSE * SM
            ENDDO
           ENDIF

!  Whole layer Greens function Discrete ordinate solutions

           CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P )   ! OUTPUT

!  Call to the post processing routine

           CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  homogeneous UP solutions, increase Z and P by 1
!  ===============================================

           Z = Z + 1
           P = P + 1
           CHANGE_EXPONENT = .TRUE.
           DO_MULTIPLIER   = .TRUE.
           FLIPPER         = .TRUE.
!           DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING         ! Wrong
           DO_SSCORR_OVERALL = .FALSE.

!  Step 1A. Tranmsittances

           PICUTOFF = NLAYERS
           INITRANS = ONE
           DELTRANS = L0_KTRANS(J,N,W_EXCIT)
           ASOURCE  = L0_KEIGEN(J,N,W_EXCIT)

!  Step 2A. Regular Solution Holding terms

           MC = L0_WHOM_MCON(J,N,WW)
           DO L = FOURIER, NRRS_MOMENTS
            SP = ZERO
            SM = ZERO
            DO ID = 1, NSTREAMS
              ID1 = ID + NSTREAMS
              LP = PLM_WT_PXI(ID,L)
              LM = PLM_WT_MXI(ID,L)
              SP = SP + MC * LP * L0_WHOM_XPOS(ID,J,N,WW)
              SM = SM + MC * LM * L0_WHOM_XPOS(ID1,J,N,WW)
            ENDDO
            HOLD_MP(L)  = SP + SM
            HOLD_RRS(L) = OMEGAMOMS_RRSLOSS(N,L,APT) * HOLD_MP(L)
           ENDDO

!  Step 3A. Regular Solution Source vector for quadrature streams

           DO I = 1, NSTREAMS
            SP = ZERO
            SM = ZERO
            DO L = FOURIER, NRRS_MOMENTS
             SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
             SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
            ENDDO
            QSOURCE_QPOS(I) = - FM_NONSHIFT_DIFFUSE * SP
            QSOURCE_QNEG(I) = - FM_NONSHIFT_DIFFUSE * SM
           ENDDO

!  Step 4A. Regular Solution Source vector for user streams

           IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
             SP = ZERO
             SM = ZERO
             DO L = FOURIER, NRRS_MOMENTS
              SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
              SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
             ENDDO
             QSOURCE_UPOS1(UI)  = - FM_NONSHIFT_DIFFUSE * SP
             QSOURCE_UNEG1(UI)  = - FM_NONSHIFT_DIFFUSE * SM
            ENDDO
           ENDIF

!  Whole layer Greens function Discrete ordinate solutions

           CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P )   ! OUTPUT

!  Call to the post processing routine

           CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  End stream loop

          ENDDO

!  End RAMAN LOSS TERM clause

        ENDIF

!  Continuation point for avoiding loss terms.

 7676   CONTINUE

!      pause 'end loss'

!   #####################################
!   #####################################
!   #    R R S  R A M A N   T E R M S   #
!   #####################################
!   #####################################

!  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  For each transition..   The Raman-shifted contributions
!  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!  Start the transitions loop

        DO S = 1, NPOINTS_LOCAL

!  BPT is the Bin pointer       in the  binning      realization
!  BPT is the wavelength itself in the monochromatic realization

         IF ( DO_BIN_REALIZATION ) THEN
          BPT = BINMAP ( S, APT )
         ELSE
          BPT = S
         ENDIF

!  In the monochromatic case, avoid for the actual point

         IF ( DO_MONO_REALIZATION .AND. S.EQ.CPT ) GO TO 5656

!  copy the OMEGAMOMS array (binned or mono)

         IF ( DO_BIN_REALIZATION ) THEN
           DO L = 0, 2
             OMEGAMOMS_RRSLOCAL(L)  = OMEGAMOMS_RRSBIN(N,L,S,APT)
           ENDDO
         ELSE
           DO L = 0, 2
             OMEGAMOMS_RRSLOCAL(L) = OMEGAMOMS_RRSGAIN(N,L,S)
           ENDDO
         ENDIF

!  flux multiplier;  which is correct factor?
!        Turns out to be DECFAC = 1.0

         FM_SHIFTED         = &
              FLUXES_RANKED ( BPT ) / FLUXES_RANKED ( CPT )
         FM_SHIFTED_DIRECT  = FM_SHIFTED
         FM_SHIFTED_DIFFUSE = FM_SHIFTED * DECFAC

!  optical thickness scaling for each shift

         DELTA_SCALING = DELTAU_VERT_INPUT(N,BPT) / &
                              DELTAU_VERT_INPUT(N,CPT)

!  Source vectors for the shifted Solar beam term
!  ==============================================

!  General Procedure in 3 Steps :
!       (1) set INITRANS, ENDTRANS, ASOURCE (transmittances and secants)
!       (2) source vectors for quadrature streams
!       (3) source vectors for user-angle streams

!  increase solution and exponent numbers by 1

         Z = Z + 1
         P = P + 1
         CHANGE_EXPONENT = .TRUE.
         DO_MULTIPLIER   = .TRUE.
         FLIPPER         = .FALSE.
         DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING

!  Step 1A. set transmittances and secants
!     Need to re-set this for each Fourier component

         PICUTOFF = BEAM_PICUTOFF(BPT)
         SAVCUT   = PICUTOFF
         ASOURCE  = BEAM_AVSECANT(N,BPT) * DELTA_SCALING
         INITRANS = BEAM_ITRANS(N,BPT)
         DELTRANS = BEAM_DTRANS(N,BPT)

!  No Step 2.

!  Step 3A, Regular solution source vector for quadrature streams

         DO I = 1, NSTREAMS
           SP = ZERO
           SM = ZERO
           DO L = FOURIER, NRRS_MOMENTS
             PHIFUNC = OMEGAMOMS_RRSLOCAL(L)
             SP = SP + PHIFUNC * PLM_00_PXI(I,L)
             SM = SM + PHIFUNC * PLM_00_MXI(I,L)
           ENDDO
           QSOURCE_QPOS(I) = + FM_SHIFTED_DIRECT * SP
           QSOURCE_QNEG(I) = + FM_SHIFTED_DIRECT * SM
         ENDDO

!  Step 4A, Regular solution source vector for user streams

         IF ( DO_USER_STREAMS ) THEN
           DO UI = 1, N_USER_STREAMS
             SP = ZERO
             SM = ZERO
             DO L = FOURIER, NRRS_MOMENTS
               PHIFUNC = OMEGAMOMS_RRSLOCAL(L)
               SP = SP + PHIFUNC * PLM_00_PXUI(UI,L)
               SM = SM + PHIFUNC * PLM_00_MXUI(UI,L)
             ENDDO
             QSOURCE_UPOS1(UI) = + FM_SHIFTED_DIRECT * SP
             QSOURCE_UNEG1(UI) = + FM_SHIFTED_DIRECT * SM
           ENDDO
         ENDIF

!  Whole layer Greens function Discrete ordinate solutions

         CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P )   ! OUTPUT

!  Call to the post processing routine

         CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Source vectors for the shifted Diffuse-field gain-terms
!  =======================================================

         IF ( .NOT. DO_DIRECTRRS_ONLY ) THEN

!  (b1) Due to particular integral of shifted elastic solution
!       ------------------------------------------------------

!  General Procedure in 4 Steps :
!       (1) set INITRANS, ENDTRANS, ASOURCE (transmittances and secants)
!       (2) integrate solution for each moment = HOLD(L)
!       (3) source vectors for quadrature streams
!       (4) source vectors for user-angle streams

!  increase solution number by 1, do not increase P

           Z = Z + 1
!           CHANGE_EXPONENT = .FALSE.                         ! Safety
           CHANGE_EXPONENT = .TRUE.
!           DO_MULTIPLIER   = .FALSE.
           DO_MULTIPLIER   = .TRUE.
!           DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING            ! Wrong
           DO_SSCORR_OVERALL = .FALSE.

!  Step 1.... (set transmittances and secants for this term)
!    This is the same as previous, only need to do this for Fourier 0

!  Step 2A, regular solution holding array

           DO L = FOURIER, NRRS_MOMENTS
             SP = ZERO
             SM = ZERO
             DO ID = 1, NSTREAMS
               ID1 = ID + NSTREAMS
               SP = SP + PLM_WT_PXI(ID,L) * L0_WPARTIC(ID1,N,BPT)
               SM = SM + PLM_WT_MXI(ID,L) * L0_WPARTIC(ID,N,BPT)
             ENDDO
             PHIFUNC     = OMEGAMOMS_RRSLOCAL(L)
             HOLD_MP(L)  = SP + SM
             HOLD_RRS(L) = PHIFUNC * HOLD_MP(L)
           ENDDO

!  Step 3A. Regular solution source vector, quadrature streams

           DO I = 1, NSTREAMS
             SP = ZERO
             SM = ZERO
             DO L = FOURIER, NRRS_MOMENTS
               SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
               SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
             ENDDO
             QSOURCE_QPOS(I) = + FM_SHIFTED_DIFFUSE * SP
             QSOURCE_QNEG(I) = + FM_SHIFTED_DIFFUSE * SM
           ENDDO

!  Step 4A. Regular solution source vector, user streams

           IF ( DO_USER_STREAMS ) THEN
             DO UI = 1, N_USER_STREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                 SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
                 SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
               ENDDO
               QSOURCE_UPOS1(UI)  = + FM_SHIFTED_DIFFUSE * SP
               QSOURCE_UNEG1(UI)  = + FM_SHIFTED_DIFFUSE * SM
             ENDDO
           ENDIF

!  Whole layer Greens function Discrete ordinate solutions

           CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P )   ! OUTPUT

!  Call to the post processing routine

           CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  (b2) Due to homogeneous functions of shifted elastic solution
!       --------------------------------------------------------

!  start stream loop

           DO J = 1, NSTREAMS

!  homogeneous DN solutions, increase Z and P by 1
!  ===============================================

            Z = Z + 1
            P = P + 1
            CHANGE_EXPONENT = .TRUE.
            DO_MULTIPLIER   = .TRUE.
            FLIPPER         = .FALSE.
!            DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING             !  Wrong
            DO_SSCORR_OVERALL = .FALSE.

!  Step 1A. Transmittances and exponent
!            --------Note scaling on the optical thickness !!

!            PICUTOFF = NLAYERS
            PICUTOFF = savcut
            INITRANS = ONE
            DELTRANS = L0_KTRANS(J,N,BPT)
            ASOURCE  = L0_KEIGEN(J,N,BPT) * DELTA_SCALING

!  Step 2A. Regular solution holding terms

            LC = L0_WHOM_LCON(J,N,BPT)
            DO L = FOURIER, NRRS_MOMENTS
              SP = ZERO
              SM = ZERO
              DO ID = 1, NSTREAMS
                ID1 = ID + NSTREAMS
                LP = PLM_WT_PXI(ID,L)
                LM = PLM_WT_MXI(ID,L)
                SP = SP + LC * LP * L0_WHOM_XPOS(ID1,J,N,BPT)
                SM = SM + LC * LM * L0_WHOM_XPOS(ID,J,N,BPT)
              ENDDO
              PHIFUNC     = OMEGAMOMS_RRSLOCAL(L)
              HOLD_MP(L)  = SP + SM
              HOLD_RRS(L) = PHIFUNC * HOLD_MP(L)
            ENDDO

!  Step 3A. Regular solution source term quadrature streams

            DO I = 1, NSTREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
                SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
                SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
              ENDDO
              QSOURCE_QPOS(I) = + FM_SHIFTED_DIFFUSE * SP
              QSOURCE_QNEG(I) = + FM_SHIFTED_DIFFUSE * SM
            ENDDO

!  Step 4A. Regular solution source vector user streams

            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                SP = ZERO
                SM = ZERO
                DO L = FOURIER, NRRS_MOMENTS
                  SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
                  SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
                ENDDO
                QSOURCE_UPOS1(UI) = + FM_SHIFTED_DIFFUSE * SP
                QSOURCE_UNEG1(UI) = + FM_SHIFTED_DIFFUSE * SM
              ENDDO
            ENDIF

!  Whole layer Greens function Discrete ordinate solutions

            CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P )   ! OUTPUT

!  Call to the post processing routine

            CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  homogeneous UP solutions, increase Z and P by 1
!  ===============================================

            Z = Z + 1
            P = P + 1
            CHANGE_EXPONENT = .TRUE.
            DO_MULTIPLIER   = .TRUE.
            FLIPPER         = .TRUE.
!            DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING            ! Wrong
            DO_SSCORR_OVERALL = .FALSE.

!  Step 1A. Transmittances and exponents

!            PICUTOFF = NLAYERS
            PICUTOFF = savcut
            INITRANS = ONE
            DELTRANS = L0_KTRANS(J,N,BPT)
            ASOURCE  = L0_KEIGEN(J,N,BPT) * DELTA_SCALING

!  Step 2A. Regular solution holding terms

            MC = L0_WHOM_MCON(J,N,BPT)
            DO L = FOURIER, NRRS_MOMENTS
              SP = ZERO
              SM = ZERO
              DO ID = 1, NSTREAMS
                ID1 = ID + NSTREAMS
                LP = PLM_WT_PXI(ID,L)
                LM = PLM_WT_MXI(ID,L)
                SP = SP + MC * LP * L0_WHOM_XPOS(ID,J,N,BPT)
                SM = SM + MC * LM * L0_WHOM_XPOS(ID1,J,N,BPT)
              ENDDO
              PHIFUNC     = OMEGAMOMS_RRSLOCAL(L)
              HOLD_MP(L)  = SP + SM
              HOLD_RRS(L) = PHIFUNC * HOLD_MP(L)
            ENDDO

!  Step 3A. Regular solution source term quadrature streams

            DO I = 1, NSTREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
                SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
                SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
              ENDDO
              QSOURCE_QPOS(I) = + FM_SHIFTED_DIFFUSE * SP
              QSOURCE_QNEG(I) = + FM_SHIFTED_DIFFUSE * SM
            ENDDO

!  Step 4A. Regular solution source vector user streams

            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                SP = ZERO
                SM = ZERO
                DO L = FOURIER, NRRS_MOMENTS
                  SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
                  SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
                ENDDO
                QSOURCE_UPOS1(UI) = + FM_SHIFTED_DIFFUSE * SP
                QSOURCE_UNEG1(UI) = + FM_SHIFTED_DIFFUSE * SM
              ENDDO
            ENDIF

!  Whole layer Greens function Discrete ordinate solutions

            CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P )   ! OUTPUT

!  Call to the post processing routine

            CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  End stream loop

           ENDDO

!  End RRS_CONTROL clause

          ENDIF

!  Continuation point for avoiding calculation (monochromatic only)

 5656     CONTINUE

!  End loop over RRS transitions (binned or mono)

        ENDDO

!  Elastic-only continuation point

 5222   continue

!  End layer loop

      ENDDO

!  Debug 17 November 2008

      if ( do_debug_sources ) then
 5    format(2i4,1p2e25.16)
      do n = 1, 18
        write(52,5)1,n,WLAYER_PIST_UP(N,1),WLAYER_PIST_DN(N,1)
        do i = 1, 16
          write(42,5)i,n,RAMAN_WUPPER(i,n),RAMAN_WLOWER(i,n)
        enddo
      enddo
      endif

!  End of subroutine

      RETURN
      END SUBROUTINE SOURCES_MASTER_1

!  End of module

      END MODULE lrrs_sources_master_1_m

