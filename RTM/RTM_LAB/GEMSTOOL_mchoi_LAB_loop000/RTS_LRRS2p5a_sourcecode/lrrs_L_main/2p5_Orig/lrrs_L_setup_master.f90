
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
! #             L_SETUP_MASTER_BIN  (master)                    #
! #             L_SETUP_MASTER_MONO (master)                    #
! #                                                             #
! #  The BIN module generates complete LRRS set-up, 3 tasks     #
! #     1. Delta-M scaling of Elastic OPs                       #
! #     2. Raman Cross-sections (Spectroscopy)                  #
! #     3. Raman Optical Properties, BIN style                  #
! #                                                             #
! #  The MONO module finishes complete LRRS set-up, 3 tasks     #
! #     1. Delta-M scaling of Elastic OPs                       #
! #     2. Raman Cross-sections (Spectroscopy)                  #
! #     3. Raman Optical Properties, MONO style                 #
! #                                                             #
! #    Monochromatic routine added 8 february 2011              #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Bookkeeping improvements (use of "Only", clearer I/O specifications)

!   -- Rob mod 5/12/17 for 2p5a, remove NMOMENTS_INPUT, add Geometry information (flags, cosines)
!   -- Rob mod 5/12/17 for 2p5a, add PHASE FUNCTION Product input (Elastic) and output (Cabannes)

      MODULE lrrs_L_setup_master_m

!      USE LRRS_PARS_m, Only : LDU

      PRIVATE
      PUBLIC :: L_SETUP_MASTER_BIN, L_SETUP_MASTER_MONO

      CONTAINS

      SUBROUTINE L_SETUP_MASTER_BIN &
        ( DO_DELTAM_SCALING, DO_ELASTIC_ONLY,                        & ! Inputs
          DO_UPWELLING, DO_DNWELLING,                                & ! Inputs
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,                    & ! Inputs
          DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION,         & ! Inputs
          DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF,    & ! Inputs
          DO_NORMALIZED_WFS, NLAYERS, NSTREAMS, NGEOMETRIES,         & ! Inputs
          NPOINTS_INNER, OFFSET_INNER, NPOINTS_OUTER,                & ! Inputs
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                        & ! Inputs
          LAMBDAS_RANKED, FLUXES_RANKED, BINLOWER, BINUPPER,         & ! Inputs
          COSSCAT_UP, COSSCAT_DN, LAYER_TEMPERATURES,                & ! Inputs
          LAYER_AIRCOLUMNS, RAYLEIGH_XSEC, RAYLEIGH_DEPOL,           & ! Inputs
          TEMPERATURES_UNSHIFTED, LAYER_AIRCOLUMNS_dT,               & ! Inputs
          DELTAU_INPUT_UNSCALED,      L_DELTAU_INPUT_UNSCALED,       & ! Inputs
          OMEGAMOMS_ELASTIC_UNSCALED, L_OMEGAMOMS_ELASTIC_UNSCALED,  & ! Inputs
          OMEGAPHASFUNC_ELASTIC_UP,   L_OMEGAPHASFUNC_ELASTIC_UP,    & ! Inputs
          OMEGAPHASFUNC_ELASTIC_DN,   L_OMEGAPHASFUNC_ELASTIC_DN,    & ! Inputs
          N_RRSBINS, BINMAP, RINGSPEC_1,                             & ! Outputs
          NMOMENTS, TRUNC_FACTORS,    L_TRUNC_FACTORS,               & ! Outputs
          DELTAU_VERT_INPUT,          L_DELTAU_VERT_INPUT,           & ! Outputs
          OMEGAMOMS_ELASTIC,          L_OMEGAMOMS_ELASTIC,           & ! Outputs
          OMEGAMOMS_CABANNES,         L_OMEGAMOMS_CABANNES,          & ! Outputs
          OMEGAPHASFUNC_CABANNES_UP,  L_OMEGAPHASFUNC_CABANNES_UP,   & ! Inputs
          OMEGAPHASFUNC_CABANNES_DN,  L_OMEGAPHASFUNC_CABANNES_DN,   & ! Inputs
          OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSLOSS_UNSCALED,  & ! Outputs
          OMEGAMOMS_RRSBIN_UNSCALED,  L_OMEGAMOMS_RRSBIN_UNSCALED,   & ! Outputs
          OMEGAMOMS_RRSLOSS,          L_OMEGAMOMS_RRSLOSS,           & ! Outputs
          OMEGAMOMS_RRSBIN,           L_OMEGAMOMS_RRSBIN,            & ! Outputs
          FAIL, MESSAGE, ACTION )                                      ! Outputs

!  LRRS include files
!  ==================

!   -- Rob mod 5/12/17 for 2p5a, two subroutines from RamanOps and LinRamanOps,
!                                add MAX_GEOMETRIES, remove MAX_MOMENTS_INPUT

      USE LRRS_PARS_m                 ,Only : FPK, MAX_LAYERS, MAX_GEOMETRIES, MAX_MOMENTS, &
                                              MAX_POINTS, MAX_BINS, MAX_ATMOSWFS
      USE LRRS_DELTAMSCALING_m        ,Only : LRRS_DELTAM_SCALING_2P2, LRRS_DELTAM_SCALING_PLUS_2P2
      USE LRRS_RAMAN_SPECTROSCOPY_m   ,Only : RRS_BASIC_SPECTROSCOPY, LRRS_RAMAN_XSECS_BIN
      USE LRRS_L_RAMAN_SPECTROSCOPY_m ,Only : LINLRRS_RAMAN_XSECS_BIN
      USE LRRS_GENERATE_RAMANOPS_m    ,Only : LRRS_RAMAN_OPS_BIN_1, LRRS_RAMAN_OPS_BIN_2
      USE LRRS_L_GENERATE_RAMANOPS_m  ,Only : LRRS_RAMAN_OPS_BIN_PLUS_1, LRRS_RAMAN_OPS_BIN_PLUS_2

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Flags
!  -----

!  Delta-m scaling flag (introduced July 24th, 2006)

      LOGICAL  , INTENT(IN) :: DO_DELTAM_SCALING

!mick fix 7/20/2016 - added flag
!  Elastic-only flag

      LOGICAL  , INTENT(IN) :: DO_ELASTIC_ONLY

!  Directional flags
!   -- Rob mod 5/12/17 for 2p5a, added

      LOGICAL  , INTENT(IN) :: DO_UPWELLING, DO_DNWELLING

!  Solution method (options are mutually exclusive)

      LOGICAL  , INTENT(IN) :: DO_CABANNES_RAMAN
      LOGICAL  , INTENT(IN) :: DO_ENERGY_BALANCING

!  Profile/column linearization control inputs

      LOGICAL  , INTENT(IN) :: DO_PROFILE_LINEARIZATION
      LOGICAL  , INTENT(IN) :: DO_COLUMN_LINEARIZATION

!  Temperature and Air Profile Linearization flags

      LOGICAL  , INTENT(IN) :: DO_TEMPPROFILE_WFS
      LOGICAL  , INTENT(IN) :: DO_AIRPROFILE_WFS

!  Temperature shift Linearization flag (Column WF only)

      LOGICAL  , INTENT(IN) :: DO_TEMPSHIFT_WF

!  Use of normalized weighting functions (Profile WFs only)
!   T-shift Jacobian is automatically normalized

      LOGICAL  , INTENT(IN) :: DO_NORMALIZED_WFS

!  Scattering angles
!  -----------------

!   -- Rob mod 5/12/17 for 2p5a, added

      REAL(FPK), INTENT(IN) :: COSSCAT_UP ( MAX_GEOMETRIES )
      REAL(FPK), INTENT(IN) :: COSSCAT_DN ( MAX_GEOMETRIES )

!  Scattering discretization
!  -------------------------

!  Number of streams, geometries
!   -- Rob mod 5/12/17 for 2p5a, added NGEOMETRIES

      INTEGER  , INTENT(IN) :: NSTREAMS, NGEOMETRIES

!  Atmospheric quantities
!  ----------------------

!  Number of layers

      INTEGER  , INTENT(IN) :: NLAYERS

!  Input layer temperatures, must be in deg K

      REAL(FPK), INTENT(IN) :: LAYER_TEMPERATURES     ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: TEMPERATURES_UNSHIFTED ( MAX_LAYERS )

!  Input layer Air columns, should be in mol/cm^2 or [DU]
!    Derivatives are always absolute values

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS_dT ( MAX_LAYERS )

!  Linearization control

      LOGICAL  , INTENT(IN) :: LAYER_VARY_FLAG   (MAX_LAYERS)
      INTEGER  , INTENT(IN) :: LAYER_VARY_NUMBER (MAX_LAYERS)

!  Wavelengths/Fluxes are defined on outer grid for the binning realizat

      REAL(FPK), INTENT(IN) :: LAMBDAS_RANKED ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: FLUXES_RANKED  ( MAX_POINTS )

!  Basic Rayleigh data
!     Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  For the binning realization
!  ---------------------------

!  outer/inner wavelength range, and offset for the inner range
!    Depends on the choice of solar spectrum

      INTEGER  , INTENT(IN) :: OFFSET_INNER
      INTEGER  , INTENT(IN) :: NPOINTS_INNER
      INTEGER  , INTENT(IN) :: NPOINTS_OUTER

!  Bin upper and lower limits

      REAL(FPK), INTENT(IN) :: BINLOWER ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BINUPPER ( MAX_POINTS )

!  Unscaled elastic scattering properties
!  ---------------------------------------

!  These must be defined on the outer wavelength grid (binning), complete grid (Mono)

!  Number of input phase function Legendre moments. THIS IS A BASIC INPUT THAT USER MUST PROVIDE
!   -- Rob mod 5/12/17 for 2p5a, No Longer required....!!!!!!!!
!      INTEGER  , INTENT(IN) :: NMOMENTS_INPUT

!  Unscaled quantities, Elastic input
!   -- Rob mod 5/12/17 for 2p5a, Changed  Phasmoms dimension to MAX_MOMENTS

      REAL(FPK), INTENT(IN) :: DELTAU_INPUT_UNSCALED      ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC_UNSCALED ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Input elastic scattering phase function products
!   -- Rob mod 5/12/17 for 2p5a, Added

      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Linearized Unscaled Optical properties, Elastic input

      REAL(FPK), INTENT(IN) :: L_DELTAU_INPUT_UNSCALED      ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_ELASTIC_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Output arguments
!  ================

!  Scaled elastic-scattering optical properties
!  --------------------------------------------

!    These must be defined on the outer wavelength grid (binning)

!  Number of moments ( = 2 * NSTREAMS )

      INTEGER  , INTENT(OUT) :: NMOMENTS

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(OUT) :: DELTAU_VERT_INPUT ( MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_ELASTIC ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Truncation factors (zero if no deltam scaling)

      REAL(FPK), INTENT(OUT) :: TRUNC_FACTORS ( MAX_LAYERS, MAX_POINTS )

!  Truncation factors

      REAL(FPK), INTENT(OUT) :: L_TRUNC_FACTORS ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

!  Vertical optical thickness values

      REAL(FPK), INTENT(OUT) :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_ELASTIC ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  RRS spectroscopic quantities
!  ----------------------------

!  Bin Mapping quantities (internal derivation)

      INTEGER  , INTENT(OUT) :: N_RRSBINS  ( MAX_POINTS )
      INTEGER  , INTENT(OUT) :: BINMAP ( MAX_BINS, MAX_POINTS )

!  Chance/Spurr Ring Spectrum. DIAGNOSTIC ONLY.

      REAL(FPK), INTENT(OUT) :: RINGSPEC_1 ( MAX_POINTS )

!  Unscaled Raman quantities
!  -------------------------

!  Output Cabannes-adjusted scattering phase function products
!   -- Rob mod 5/12/17 for 2p5a, Replaces the phasmoms output (OMEGAMOMS_CABANNES_UNSCALED)
!                                Only calculated for the Cabannes-Raman option

      REAL(FPK), INTENT(OUT) :: OMEGAPHASFUNC_CABANNES_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAPHASFUNC_CABANNES_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS_UNSCALED ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSBIN_UNSCALED  ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Linearized

      REAL(FPK), INTENT(OUT) :: L_OMEGAPHASFUNC_CABANNES_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAPHASFUNC_CABANNES_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSLOSS_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSBIN_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Scaled Raman-scattering optical properties
!  ------------------------------------------

!  Derived Cabannes scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived rotational raman scattering (Internal variable). Loss term
!    Single scatter albedo * phase function moments
!    Only required for the Energy-balance approximation.

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived rotational Raman scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSBIN  ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Linearizations

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_CABANNES ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSLOSS  ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSBIN   ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@
!  Exception handling for MAX_BINS Dimension

      LOGICAL          , INTENT(OUT) :: FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE, ACTION

!  Local Arrays
!  ============

!  Logical flag for linearization of RRS spectroscopy

      LOGICAL :: DO_LINRRS

!  shift dimensions

      INTEGER  , PARAMETER :: MAXN2 = 48, MAXO2 = 185

!  Positions

      REAL(FPK) :: O2POS (MAXO2)
      REAL(FPK) :: N2POS (MAXN2)

!  Gammas

      REAL(FPK) :: GAMMA_O2(3)
      REAL(FPK) :: GAMMA_N2(3)

!  Basic coefficients for O2 and N2

      REAL(FPK) :: RRS_N2REF ( MAX_LAYERS, MAXN2 )
      REAL(FPK) :: RRS_O2REF ( MAX_LAYERS, MAXO2 )

!  T-derivatives of Basic coefficients for O2 and N2

      REAL(FPK) :: LINRRS_N2REF ( MAX_LAYERS, MAXN2 )
      REAL(FPK) :: LINRRS_O2REF ( MAX_LAYERS, MAXO2 )

!  Loss term cross-section and number of shifts = N_Stokes + N_Antistoke

      REAL(FPK) :: RRSXSEC_OUT ( MAX_LAYERS, MAX_POINTS )

!  Binned Gain term cross-sections

      REAL(FPK) :: BINXSEC ( MAX_LAYERS, MAX_POINTS, MAX_BINS )

!  T-derivative of Binned Gain term cross-sections

      REAL(FPK) :: BINXSEC_dT ( MAX_LAYERS, MAX_POINTS, MAX_BINS )

!  T-derivative of Loss term cross sections

      REAL(FPK) :: RRSXSEC_OUT_dT ( MAX_LAYERS, MAX_POINTS )

!  debug

!      INTEGER :: n, q, cpt

!  DELTA-M Scaling
!  ===============

!  REMARK: In Version 2.2, this is automatic, controlled only
!          by the DO_DELTAM_SCALING flag.

!  -- Rob mod 5/12/17 for 2p5a, remove MAXMOMENTS_INPUT, NMOMENTS_INPUT, set NMOMENTS first

      NMOMENTS = 2*NSTREAMS - 1

      IF(do_profile_linearization.or.do_column_linearization)THEN
        CALL LRRS_DELTAM_SCALING_PLUS_2P2 &
         ( MAX_LAYERS, MAX_MOMENTS, MAX_POINTS, MAX_ATMOSWFS,     & ! Inputs
           DO_DELTAM_SCALING, NPOINTS_OUTER, NLAYERS, NSTREAMS,   & ! Inputs
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                    & ! Inputs
           DELTAU_INPUT_UNSCALED,   OMEGAMOMS_ELASTIC_UNSCALED,   & ! Inputs
           L_DELTAU_INPUT_UNSCALED, L_OMEGAMOMS_ELASTIC_UNSCALED, & ! Inputs
           TRUNC_FACTORS,         L_TRUNC_FACTORS,                & ! Outputs
           DELTAU_VERT_INPUT,     OMEGAMOMS_ELASTIC,              & ! Outputs
           L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC )               ! Outputs
      ELSE
        CALL LRRS_DELTAM_SCALING_2P2 &
         ( MAX_LAYERS, MAX_MOMENTS, MAX_POINTS,                 & ! Inputs
           DO_DELTAM_SCALING, NPOINTS_OUTER, NLAYERS, NSTREAMS, & ! Inputs
           DELTAU_INPUT_UNSCALED, OMEGAMOMS_ELASTIC_UNSCALED,   & ! Inputs
           TRUNC_FACTORS, DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC )  ! Outputs
      ENDIF

!mick fix 7/20/2016 - added if condition for elastic-only case

      IF ( DO_ELASTIC_ONLY ) RETURN

!  Generation of Raman-scattering Optical Properties (ROPS)
!  ========================================================

!  REMARK: THIS is an Automatic procedure in VERSION 2.2

!  THERE ARE 2 MAIN STEPS
!     (1) SPECTROSCOPY - generate the Raman cross-sections, optional T-derivatives
!     (2) MAKE ROPS    - calculate the additional OPs

!  STEP 1. Raman Cross-sections generator
!  --------------------------------------

!  Get the Raman spectroscopy

      DO_LINRRS = DO_TEMPPROFILE_WFS .or. DO_TEMPSHIFT_WF
      CALL RRS_BASIC_SPECTROSCOPY &
         ( MAX_LAYERS, NLAYERS, DO_LINRRS,                    & ! Inputs
           LAYER_TEMPERATURES, MAXN2, MAXO2,                  & ! Inputs
           N2POS, O2POS, GAMMA_N2, GAMMA_O2,                  & ! Outputs
           RRS_N2REF, RRS_O2REF, LINRRS_N2REF, LINRRS_O2REF )   ! Outputs

!  Given outer and inner windows and wavelengths + bins and solar fluxes
!  this routine generates the following Raman cross-sections:

!   Each wavelength w,                 the number of bins (N_RRSBINS)
!   Each wavelength w, bin b,          the binning map BINMAP(w,b)
!   Each wavelength w, bin b, layer n, Raman gain X-section BINXSEC(w,b,
!   Each wavelength w, layer n,        Raman loss X-section RRSXSEC_OUT(

!  RINGSPEC_1 is the Fraunhofer "Ring spectrum" for this input field.
!    --This is the standard SAO calculation, output as a diagnostic here

!  Two choices, with and without temperature derivatives
!    Chain rule: Temperatures will depend in turn on heights

      IF ( DO_LINRRS ) THEN
        CALL LINLRRS_RAMAN_XSECS_BIN &
        ( MAX_LAYERS, MAX_POINTS, MAX_BINS,                    & ! Inputs
          NLAYERS, NPOINTS_INNER, NPOINTS_OUTER, OFFSET_INNER, & ! Inputs
          BINLOWER, BINUPPER, LAMBDAS_RANKED, FLUXES_RANKED,   & ! Inputs
          MAXN2, MAXO2, RRS_N2REF, RRS_O2REF, LINRRS_N2REF,    & ! Inputs
          LINRRS_O2REF, N2POS, O2POS, GAMMA_N2, GAMMA_O2,      & ! Inputs
          N_RRSBINS, BINMAP, BINXSEC, RRSXSEC_OUT, RINGSPEC_1, & ! Output
          BINXSEC_dT, RRSXSEC_OUT_dT,                          & ! Output
          FAIL, MESSAGE, ACTION )                                ! Output
      ELSE
        CALL LRRS_RAMAN_XSECS_BIN &
        ( MAX_LAYERS, MAX_POINTS, MAX_BINS,                    & ! Inputs
          NLAYERS, NPOINTS_INNER, NPOINTS_OUTER, OFFSET_INNER, & ! Inputs
          BINLOWER, BINUPPER, LAMBDAS_RANKED, FLUXES_RANKED,   & ! Inputs
          MAXN2, MAXO2, RRS_N2REF, RRS_O2REF,                  & ! Inputs
          N2POS, O2POS, GAMMA_N2, GAMMA_O2,                    & ! Inputs
          N_RRSBINS, BINMAP, BINXSEC, RRSXSEC_OUT, RINGSPEC_1, & ! Outputs
          FAIL, MESSAGE, ACTION )                                ! Outputs
      ENDIF

!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@
!  Exception handling for MAX_BINS Dimension. Exit immediately

          IF ( FAIL) RETURN

!  STEP 2. Make the Raman Optical Depths
!  -------------------------------------

!  Rob mod 5/12/17 for 2p5a, Separate routines for Scaled and Unscaled output
!     Unscaled inputs, use LRRS_RAMAN_OPS_BIN_1 or LRRS_RAMAN_OPS_BIN_PLUS_1

!  Formerly, this was the responsibility of the User
!  Now taken care of automatically with these two routines:

!          lrrs_raman_ops_bin_plus (with linearizations);
!          lrrs_raman_ops_bin      (no linearization)

!  These routines are completely stand-alone.  In order to make them work,
!  the user must supply the following:

!      Dimensioning inputs
!      LRRS method (Energy-balancing vs. Cabannes). First is usual.
!      Number of layers, Inner window points and offset, # of bins
!      Number of phase function expansion coefficients (Elastic)
!      Number of product phase functions and scattering geometries (Elastic). Routine 1
!      Number of product phase function expansion coefficients     (Elastic). Routine 2
!      Rayleigh X-secs(**) + depolarizations, Layer Air columns(**)
!      Raman Loss (RRSXSEC_OUT) and binned gain (BINXSEC) cross-sections

!  For the linearized optical properties, need EXTRA user-supplied info:

!      linearization control
!      Linearized Elastic optical properties as defined above

!  (**) Note - cross sections must be [cm^2/mol], Air columns [mol/cm^2]

      if (do_profile_linearization.or.do_column_linearization) then
        CALL LRRS_RAMAN_OPS_BIN_PLUS_1 &
          ( MAX_LAYERS, MAX_POINTS, MAX_GEOMETRIES, MAX_BINS, MAX_ATMOSWFS,            & ! Input dimensions
            DO_UPWELLING, DO_DNWELLING, DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,        & ! Input flags
            DO_NORMALIZED_WFS, DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF, & ! Input flags
            NLAYERS, NGEOMETRIES, NPOINTS_INNER, OFFSET_INNER, LAYER_VARY_NUMBER,      & ! Input integers
            COSSCAT_UP, COSSCAT_DN, LAYER_TEMPERATURES, TEMPERATURES_UNSHIFTED,        & ! Inputs
            LAYER_AIRCOLUMNS, LAYER_AIRCOLUMNS_DT, RAYLEIGH_XSEC, RAYLEIGH_DEPOL,      & ! Inputs
            DELTAU_VERT_INPUT,         L_DELTAU_VERT_INPUT,                            & ! Inputs Optical
            OMEGAPHASFUNC_ELASTIC_UP,  L_OMEGAPHASFUNC_ELASTIC_UP,                     & ! Inputs Optical
            OMEGAPHASFUNC_ELASTIC_DN,  L_OMEGAPHASFUNC_ELASTIC_DN,                     & ! Inputs Optical
            N_RRSBINS, BINXSEC, RRSXSEC_OUT, BINXSEC_DT, RRSXSEC_OUT_DT,               & ! Inputs Raman Spec
            OMEGAPHASFUNC_CABANNES_UP,  L_OMEGAPHASFUNC_CABANNES_UP,                   & ! Outputs
            OMEGAPHASFUNC_CABANNES_DN,  L_OMEGAPHASFUNC_CABANNES_DN,                   & ! Outputs
            OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSLOSS_UNSCALED,                  & ! Outputs
            OMEGAMOMS_RRSBIN_UNSCALED,  L_OMEGAMOMS_RRSBIN_UNSCALED )                    ! Outputs
      else
        call lrrs_raman_ops_bin_1 &
          ( MAX_LAYERS, MAX_POINTS, MAX_GEOMETRIES, MAX_BINS,                   & ! Input dimensions
            DO_UPWELLING, DO_DNWELLING, DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, & ! Input flags
            NLAYERS, NGEOMETRIES, NPOINTS_INNER, OFFSET_INNER,                  & ! Input integers
            COSSCAT_UP, COSSCAT_DN, LAYER_AIRCOLUMNS,              & ! Input 
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL, DELTAU_INPUT_UNSCALED,  & ! Inputs Optical
            OMEGAPHASFUNC_ELASTIC_UP, OMEGAPHASFUNC_ELASTIC_DN,    & ! Inputs Optical
            N_RRSBINS, BINXSEC, RRSXSEC_OUT,                       & ! Inputs Raman Spec.
            OMEGAPHASFUNC_CABANNES_UP, OMEGAPHASFUNC_CABANNES_DN,  & ! Outputs Adjusted Cabannes
            OMEGAMOMS_RRSLOSS_UNSCALED, OMEGAMOMS_RRSBIN_UNSCALED )  ! Outputs Raman Loss//Gain
      endif

!  Call the Raman IOP+ generator, scaled

      if(do_profile_linearization.or.do_column_linearization)then
       CALL LRRS_RAMAN_OPS_BIN_PLUS_2 &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, MAX_BINS, MAX_ATMOSWFS, & ! Inputs
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, DO_NORMALIZED_WFS,   & ! Inputs
          DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF,      & ! Inputs
          NLAYERS, NMOMENTS, NPOINTS_INNER, OFFSET_INNER,              & ! Inputs
          LAYER_VARY_NUMBER, RAYLEIGH_XSEC, RAYLEIGH_DEPOL,            & ! Inputs
          LAYER_TEMPERATURES, TEMPERATURES_UNSHIFTED,                  & ! Inputs
          LAYER_AIRCOLUMNS, LAYER_AIRCOLUMNS_dT,                       & ! Inputs
          DELTAU_VERT_INPUT,     OMEGAMOMS_ELASTIC,                    & ! Inputs
          L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC,                    & ! Inputs
          N_RRSBINS, BINXSEC, RRSXSEC_OUT, BINXSEC_dT, RRSXSEC_OUT_dT, & ! Inputs
          OMEGAMOMS_CABANNES, L_OMEGAMOMS_CABANNES,                    & ! Outputs
          OMEGAMOMS_RRSLOSS,  L_OMEGAMOMS_RRSLOSS,                     & ! Outputs
          OMEGAMOMS_RRSBIN,   L_OMEGAMOMS_RRSBIN )                       ! Outputs
      else
       CALL LRRS_RAMAN_OPS_BIN_2 &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, MAX_BINS,   & ! Inputs
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,          & ! Inputs
          NLAYERS, NMOMENTS, NPOINTS_INNER, OFFSET_INNER,  & ! Inputs
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, LAYER_AIRCOLUMNS, & ! Inputs
          DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC,            & ! Inputs
          N_RRSBINS, BINXSEC, RRSXSEC_OUT,                 & ! Inputs
          OMEGAMOMS_CABANNES, OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN ) ! Outputs
      endif

!  Finish

      END SUBROUTINE L_SETUP_MASTER_BIN

!

      SUBROUTINE L_SETUP_MASTER_MONO &
        ( DO_DELTAM_SCALING, DO_ELASTIC_ONLY,                        & ! Inputs
          DO_UPWELLING, DO_DNWELLING,                                & ! Inputs
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,                    & ! Inputs
          DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION,         & ! Inputs
          DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF,    & ! Inputs
          DO_NORMALIZED_WFS, NLAYERS, NSTREAMS, NGEOMETRIES,         & ! Inputs
          LAMBDA_EXCIT, NPOINTS_MONO, W_EXCIT,                       & ! Inputs
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                        & ! Inputs
          COSSCAT_UP, COSSCAT_DN, LAYER_TEMPERATURES,                & ! Inputs
          TEMPERATURES_UNSHIFTED, LAYER_AIRCOLUMNS,                  & ! Inputs
          LAYER_AIRCOLUMNS_dT, RAYLEIGH_XSEC, RAYLEIGH_DEPOL,        & ! Inputs
          DELTAU_INPUT_UNSCALED,      L_DELTAU_INPUT_UNSCALED,       & ! Inputs
          OMEGAMOMS_ELASTIC_UNSCALED, L_OMEGAMOMS_ELASTIC_UNSCALED,  & ! Inputs
          OMEGAPHASFUNC_ELASTIC_UP,   L_OMEGAPHASFUNC_ELASTIC_UP,    & ! Inputs
          OMEGAPHASFUNC_ELASTIC_DN,   L_OMEGAPHASFUNC_ELASTIC_DN,    & ! Inputs
          NMOMENTS, TRUNC_FACTORS,    L_TRUNC_FACTORS,               & ! Outputs
          DELTAU_VERT_INPUT,          L_DELTAU_VERT_INPUT,           & ! Outputs
          OMEGAMOMS_ELASTIC,          L_OMEGAMOMS_ELASTIC,           & ! Outputs
          OMEGAMOMS_CABANNES,         L_OMEGAMOMS_CABANNES,          & ! Outputs
          OMEGAPHASFUNC_CABANNES_UP,  L_OMEGAPHASFUNC_CABANNES_UP,   & ! Inputs
          OMEGAPHASFUNC_CABANNES_DN,  L_OMEGAPHASFUNC_CABANNES_DN,   & ! Inputs
          OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSLOSS_UNSCALED,  & ! Outputs
          OMEGAMOMS_RRSGAIN_UNSCALED, L_OMEGAMOMS_RRSGAIN_UNSCALED,  & ! Outputs
          OMEGAMOMS_RRSLOSS,          L_OMEGAMOMS_RRSLOSS,           & ! Outputs
          OMEGAMOMS_RRSGAIN,          L_OMEGAMOMS_RRSGAIN )            ! Outputs

!  LRRS include files
!  ==================

!   -- Rob mod 5/12/17 for 2p5a, two subroutines from RamanOps and LinRamanOps,
!                                add MAX_GEOMETRIES, remove MAX_MOMENTS_INPUT

      USE LRRS_PARS_m                 ,Only : FPK, MAX_LAYERS, MAX_GEOMETRIES, MAX_MOMENTS, MAX_POINTS, MAX_ATMOSWFS
      USE LRRS_DELTAMSCALING_m        ,Only : LRRS_DELTAM_SCALING_2P2, LRRS_DELTAM_SCALING_PLUS_2P2
      USE LRRS_RAMAN_SPECTROSCOPY_m   ,Only : RRS_BASIC_SPECTROSCOPY, LRRS_RAMAN_XSECS_MONO
      USE LRRS_GENERATE_RAMANOPS_m    ,Only : LRRS_RAMAN_OPS_MONO_1, LRRS_RAMAN_OPS_MONO_2

      USE LRRS_L_RAMAN_SPECTROSCOPY_m ,Only : LINLRRS_RAMAN_XSECS_MONO
      USE LRRS_L_GENERATE_RAMANOPS_m  ,Only : LRRS_RAMAN_OPS_MONO_PLUS_1, LRRS_RAMAN_OPS_MONO_PLUS_2

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Flags
!  -----

!  Delta-m scaling flag (introduced July 24th, 2006)

      LOGICAL  , INTENT(IN) :: DO_DELTAM_SCALING

!mick fix 7/20/2016 - added flag
!  Elastic-only flag

      LOGICAL  , INTENT(IN) :: DO_ELASTIC_ONLY

!  Directional flags
!   -- Rob mod 5/12/17 for 2p5a, added

      LOGICAL  , INTENT(IN) :: DO_UPWELLING, DO_DNWELLING

!  Solution method (options are mutually exclusive)

      LOGICAL  , INTENT(IN) :: DO_CABANNES_RAMAN
      LOGICAL  , INTENT(IN) :: DO_ENERGY_BALANCING

!  Profile/column linearization control inputs

      LOGICAL  , INTENT(IN) :: do_profile_linearization
      LOGICAL  , INTENT(IN) :: do_column_linearization

!  Temperature and Air Profile Linearization flags

      LOGICAL  , INTENT(IN) :: DO_TEMPPROFILE_WFS
      LOGICAL  , INTENT(IN) :: DO_AIRPROFILE_WFS

!  Temperature shift Linearization flag (Column WF only)

      LOGICAL  , INTENT(IN) :: DO_TEMPSHIFT_WF

!  Use of normalized weighting functions (Profile WFs only)
!   T-shift Jacobian is automatically normalized

      LOGICAL  , INTENT(IN) :: DO_NORMALIZED_WFS

!  Scattering discretization
!  -------------------------

!   -- Rob mod 5/12/17 for 2p5a, added NGEOMETRIES

      INTEGER  , INTENT(IN) :: NSTREAMS, NGEOMETRIES

!  Scattering angles
!  -----------------

!   -- Rob mod 5/12/17 for 2p5a, added

      REAL(FPK), INTENT(IN) :: COSSCAT_UP ( MAX_GEOMETRIES )
      REAL(FPK), INTENT(IN) :: COSSCAT_DN ( MAX_GEOMETRIES )

!  Atmospheric quantities
!  ----------------------

!  Number of layers

      INTEGER  , INTENT(IN) :: NLAYERS

!  Input layer temperatures, must be in deg K

      REAL(FPK), INTENT(IN) :: LAYER_TEMPERATURES     ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: TEMPERATURES_UNSHIFTED ( MAX_LAYERS )

!  Input layer Air columns, should be in mol/cm^2 or [DU]

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS_dT ( MAX_LAYERS )

!  Linearization control

      LOGICAL  , INTENT(IN) :: LAYER_VARY_FLAG   (MAX_LAYERS)
      INTEGER  , INTENT(IN) :: LAYER_VARY_NUMBER (MAX_LAYERS)

!  MONOCHROMATIC Wavelengths
!  -------------------------

!  Excitation wavelength

      REAL(FPK), INTENT(IN) :: LAMBDA_EXCIT

!  Ranked wavelengths and NPOINTS_MONO = 234 points
!  W_EXCIT = position of excitation wavelength in the ranked set
!    NOTE: These are inputs, and should be calculated in the Iopsetup
!          routine RRS_LAMDASRANKED_MONO.

      INTEGER  , INTENT(IN) :: NPOINTS_MONO, W_EXCIT
!      REAL(FPK), INTENT(IN) :: LAMBDAS_RANKED  ( MAX_POINTS )

!  Basic Rayleigh data
!     Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  Unscaled elastic scattering properties
!  ---------------------------------------

!  These must be defined on the outer wavelength grid (binning)

!  Number of input phase function Legendre moments. THIS IS A BASIC INPUT THAT USER MUST PROVIDE
!   -- Rob mod 5/12/17 for 2p5a, No Longer required....!!!!!!!!
!      INTEGER  , INTENT(IN) :: NMOMENTS_INPUT

!  Unscaled quantities, Elastic input
!   -- Rob mod 5/12/17 for 2p5a, Changed  Phasmoms dimension to MAX_MOMENTS

      REAL(FPK), INTENT(IN) :: DELTAU_INPUT_UNSCALED ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC_UNSCALED ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Input elastic scattering phase function products
!   -- Rob mod 5/12/17 for 2p5a, Added

      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Linearized Unscaled Optical properties,Elastic input
!   -- Rob mod 5/12/17 for 2p5a, Changed  Phasmoms dimension to MAX_MOMENTS, Added Phasfuncs

      REAL(FPK), INTENT(IN) :: L_DELTAU_INPUT_UNSCALED      ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_ELASTIC_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Output arguments
!  ================

!  Scaled elastic-scattering optical properties
!  --------------------------------------------

!    These must be defined on the outer wavelength grid (binning), complete grid (mono)

!  Number of moments ( = 2 * NSTREAMS )

      INTEGER  , INTENT(OUT) :: NMOMENTS

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(OUT) :: DELTAU_VERT_INPUT   ( MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_ELASTIC ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Truncation factors (zero if no deltam scaling)

      REAL(FPK), INTENT(OUT) :: TRUNC_FACTORS         ( MAX_LAYERS, MAX_POINTS )

!  Truncation factors

      REAL(FPK), INTENT(OUT) :: L_TRUNC_FACTORS       ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

!  Vertical optical thickness values

      REAL(FPK), INTENT(OUT) :: L_DELTAU_VERT_INPUT   ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_ELASTIC ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Unscaled Raman quantities
!  -------------------------

!  Output Cabannes-adjusted scattering phase function products
!   -- Rob mod 5/12/17 for 2p5a, Replaces the phasmoms output (OMEGAMOMS_CABANNES_UNSCALED)
!                                Only calculated for the Cabannes-Raman option

      REAL(FPK), INTENT(OUT) :: OMEGAPHASFUNC_CABANNES_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAPHASFUNC_CABANNES_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS_UNSCALED ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSGAIN_UNSCALED ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Linearized

      REAL(FPK), INTENT(OUT) :: L_OMEGAPHASFUNC_CABANNES_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAPHASFUNC_CABANNES_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSLOSS_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSGAIN_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Scaled Raman-scattering optical properties
!  ------------------------------------------

!  Derived Cabannes scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived rotational raman scattering (Internal variable). Loss term
!    Single scatter albedo * phase function moments
!    Only required for the Energy-balance approximation.

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived rotational Raman scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSGAIN ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Linearizations

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_CABANNES ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSLOSS  ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSGAIN  ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Local Arrays
!  ============

!  Logical flag for linearization of RRS spectroscopy

      LOGICAL :: DO_LINRRS

!  shift dimensions

      INTEGER  , PARAMETER :: MAXN2 = 48, MAXO2 = 185

!  Positions

      REAL(FPK) :: O2POS (MAXO2)
      REAL(FPK) :: N2POS (MAXN2)

!  Gammas

      REAL(FPK) :: GAMMA_O2(3)
      REAL(FPK) :: GAMMA_N2(3)

!  Basic coefficients for O2 and N2

      REAL(FPK) :: RRS_N2REF ( MAX_LAYERS, MAXN2 )
      REAL(FPK) :: RRS_O2REF ( MAX_LAYERS, MAXO2 )

!  T-derivatives of Basic coefficients for O2 and N2

      REAL(FPK) :: LINRRS_N2REF ( MAX_LAYERS, MAXN2 )
      REAL(FPK) :: LINRRS_O2REF ( MAX_LAYERS, MAXO2 )

!  Gain term cross-sections

      REAL(FPK) :: RRSXSEC_RANKED ( MAX_LAYERS, MAX_POINTS )

!  Loss term cross sections

      REAL(FPK) :: RRSXSEC_OUT_MONO    ( MAX_LAYERS )

!  Temperature linearizations of cross-sections

      REAL(FPK) :: RRSXSEC_RANKED_dT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: RRSXSEC_OUT_MONO_dT    ( MAX_LAYERS )

!  debug
!      INTEGER :: n, w

!  DELTA-M Scaling
!  ===============

!  REMARK: In Version 2.2, this is automatic, controlled only
!          by the DO_DELTAM_SCALING flag.

!  -- Rob mod 5/12/17 for 2p5a, remove MAXMOMENTS_INPUT, NMOMENTS_INPUT, set NMOMENTS first

      NMOMENTS = 2*NSTREAMS - 1

      IF ( do_profile_linearization .or. do_column_linearization ) THEN
        CALL LRRS_DELTAM_SCALING_PLUS_2P2 &
         ( MAX_LAYERS, MAX_MOMENTS, MAX_POINTS, MAX_ATMOSWFS,     & ! Inputs
           DO_DELTAM_SCALING, NPOINTS_MONO, NLAYERS, NSTREAMS,    & ! Inputs
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                    & ! Inputs
           DELTAU_INPUT_UNSCALED,   OMEGAMOMS_ELASTIC_UNSCALED,   & ! Inputs
           L_DELTAU_INPUT_UNSCALED, L_OMEGAMOMS_ELASTIC_UNSCALED, & ! Inputs
           TRUNC_FACTORS,         L_TRUNC_FACTORS,                & ! Outputs
           DELTAU_VERT_INPUT,     OMEGAMOMS_ELASTIC,              & ! Outputs
           L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC )               ! Outputs
      ELSE
        CALL LRRS_DELTAM_SCALING_2P2 &
         ( MAX_LAYERS, MAX_MOMENTS, MAX_POINTS,                  & ! Inputs
           DO_DELTAM_SCALING, NPOINTS_MONO, NLAYERS, NSTREAMS,   & ! Inputs
           DELTAU_INPUT_UNSCALED, OMEGAMOMS_ELASTIC_UNSCALED,    & ! Inputs
           TRUNC_FACTORS, DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC )   ! Outputs
      ENDIF

!mick fix 7/20/2016 - added if condition for elastic-only case

      IF ( DO_ELASTIC_ONLY ) RETURN

!  Generation of Raman-scattering Optical Properties (ROPS)
!  ========================================================

!  STEP 1. Raman Cross-sections generator
!  --------------------------------------

!  REMARK: THIS is an Automatic procedure in VERSION 2.2

!  Get the Raman spectroscopy

      DO_LINRRS = DO_TEMPPROFILE_WFS .or. DO_TEMPSHIFT_WF
      CALL RRS_BASIC_SPECTROSCOPY &
         ( MAX_LAYERS, NLAYERS, DO_LINRRS,                    & ! Inputs
           LAYER_TEMPERATURES, MAXN2, MAXO2,                  & ! Inputs
           N2POS, O2POS, GAMMA_N2, GAMMA_O2,                  & ! Outputs
           RRS_N2REF, RRS_O2REF, LINRRS_N2REF, LINRRS_O2REF )   ! Outputs

!  Given outer and inner windows and wavelengths + bins and solar fluxes
!  this routine generates the following Raman cross-sections:

!   Wavelength LAMBDA_EXCIT, --> Ranked wavelengths, and excitation poin
!                           (Ranked) Raman gain X-sections RRSXSEC_RANKE
!                                    Raman loss X-section RRSXSEC_OUT_MO

!  Two choices, with and without temperature derivatives
!    Chain rule: Temperatures will depend in turn on heights

      IF ( DO_LINRRS ) THEN
        call LINLRRS_RAMAN_XSECS_MONO &
        ( MAX_LAYERS, MAX_POINTS, NLAYERS, LAMBDA_EXCIT,  & ! Inputs
          MAXN2, MAXO2, RRS_N2REF, RRS_O2REF,             & ! Inputs
          LINRRS_N2REF, LINRRS_O2REF,                     & ! Inputs
          N2POS, O2POS, GAMMA_N2, GAMMA_O2,               & ! Inputs
          RRSXSEC_OUT_MONO,    RRSXSEC_RANKED,            & ! Outputs
          RRSXSEC_OUT_MONO_dT, RRSXSEC_RANKED_dT )          ! Outputs
      ELSE
        call LRRS_RAMAN_XSECS_MONO &
        ( MAX_LAYERS, MAX_POINTS, NLAYERS, LAMBDA_EXCIT, & ! Inputs
          MAXN2, MAXO2, RRS_N2REF, RRS_O2REF,            & ! Inputs
          N2POS, O2POS, GAMMA_N2, GAMMA_O2,              & ! Inputs
          RRSXSEC_OUT_MONO, RRSXSEC_RANKED )               ! Outputs
      ENDIF

!  STEP 2. Make the Raman Optical Depths
!  -------------------------------------

!  Formerly, this was the responsibility of the User
!  Now taken care of automatically with these two routines:

!          lrrs_raman_ops_mono_plus (with linearizations);
!          lrrs_raman_ops_mono      (no linearization)

!  These routines are completely stand-alone.  In order to make them work,
!  the user must supply the following:

!      LRRS method (Energy-balancing vs. Cabannes). First is usual.
!      Number of layers, Inner window points and excitation position
!      Number of Input vertical optical depths
!      Number of phase function expansion coefficients (Elastic)
!      Number of product phase functions and scattering geometries (Elastic). Routine 1
!      Number of product phase function expansion coefficients     (Elastic). Routine 2
!      Rayleigh X-secs(**) + depolarizations, Layer Air columns(**)
!      Loss (RRSXSEC_OUT_MONO) and gain (RRSXSEC_RANKED) cross-sections(

!  For the linearized optical properties, need EXTRA user-supplied info:

!      linearization control
!      Linearized Elastic optical properties as defined above

!  (**) Note - cross sections must be [cm^2/mol], Air columns [mol/cm^2]

!  Call the Raman IOP+ generator (delta-m scaled)

      if ( do_profile_linearization .or. do_column_linearization ) then
       CALL LRRS_RAMAN_OPS_MONO_PLUS_1 &
        ( MAX_LAYERS, MAX_POINTS, MAX_GEOMETRIES, MAX_ATMOSWFS,                      & ! Input dimensions
          DO_UPWELLING, DO_DNWELLING, DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,        & ! Input flags
          DO_NORMALIZED_WFS, DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF, & ! Input flags
          NLAYERS, NGEOMETRIES, NPOINTS_MONO, W_EXCIT, LAYER_VARY_NUMBER,            & ! Inputs
          COSSCAT_UP, COSSCAT_DN, LAYER_TEMPERATURES, TEMPERATURES_UNSHIFTED,        & ! Inputs
          LAYER_AIRCOLUMNS, LAYER_AIRCOLUMNS_DT, RAYLEIGH_XSEC, RAYLEIGH_DEPOL,      & ! Inputs
          DELTAU_VERT_INPUT,         L_DELTAU_VERT_INPUT,                             & ! Inputs Optical
          OMEGAPHASFUNC_ELASTIC_UP,  L_OMEGAPHASFUNC_ELASTIC_UP,                      & ! Inputs Optical
          OMEGAPHASFUNC_ELASTIC_DN,  L_OMEGAPHASFUNC_ELASTIC_DN,                      & ! Inputs Optical
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO, RRSXSEC_RANKED_DT, RRSXSEC_OUT_MONO_DT,   & ! Inputs
          OMEGAPHASFUNC_CABANNES_UP,  L_OMEGAPHASFUNC_CABANNES_UP,                    & ! Outputs
          OMEGAPHASFUNC_CABANNES_DN,  L_OMEGAPHASFUNC_CABANNES_DN,                    & ! Outputs
          OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSLOSS_UNSCALED,                   & ! Output
          OMEGAMOMS_RRSGAIN_UNSCALED, L_OMEGAMOMS_RRSGAIN_UNSCALED )                    ! Outputs
      else
       CALL LRRS_RAMAN_OPS_MONO_1 &
        ( MAX_LAYERS, MAX_POINTS, MAX_GEOMETRIES,                             & ! Input dimensions
          DO_UPWELLING, DO_DNWELLING, DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, & ! Input flags
          NLAYERS, NGEOMETRIES, NPOINTS_MONO, W_EXCIT,            & ! Input integers
          COSSCAT_UP, COSSCAT_DN, LAYER_AIRCOLUMNS,               & ! Input 
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, DELTAU_INPUT_UNSCALED,   & ! Inputs Optical
          OMEGAPHASFUNC_ELASTIC_UP, OMEGAPHASFUNC_ELASTIC_DN,     & ! Inputs Optical
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO,                       & ! Inputs Raman Spec.
          OMEGAPHASFUNC_CABANNES_UP, OMEGAPHASFUNC_CABANNES_DN,   & ! Outputs Adjusted Cabannes
          OMEGAMOMS_RRSLOSS_UNSCALED, OMEGAMOMS_RRSGAIN_UNSCALED )  ! Outputs Raman Loss//Gain
      endif

!  Call the Raman IOP generator, scaled
!   Rob mod 5/12/17 for 2p5a, Only the name changed

      if ( do_profile_linearization .or. do_column_linearization ) then
       CALL LRRS_RAMAN_OPS_MONO_PLUS_2 &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, MAX_ATMOSWFS,         & ! Inputs
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, DO_NORMALIZED_WFS, & ! Inputs
          DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF,    & ! Inputs
          NLAYERS, NMOMENTS, NPOINTS_MONO, W_EXCIT,                  & ! Inputs
          LAYER_VARY_NUMBER, RAYLEIGH_XSEC, RAYLEIGH_DEPOL,          & ! Inputs
          LAYER_TEMPERATURES, TEMPERATURES_UNSHIFTED,                & ! Inputs
          LAYER_AIRCOLUMNS, LAYER_AIRCOLUMNS_dT,                     & ! Inputs
          DELTAU_VERT_INPUT,     OMEGAMOMS_ELASTIC,                  & ! Inputs
          L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC,                  & ! Inputs
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO,                          & ! Inputs
          RRSXSEC_RANKED_dT, RRSXSEC_OUT_MONO_dT,                    & ! Inputs
          OMEGAMOMS_CABANNES,L_OMEGAMOMS_CABANNES,                   & ! Outputs
          OMEGAMOMS_RRSLOSS, L_OMEGAMOMS_RRSLOSS,                    & ! Outputs
          OMEGAMOMS_RRSGAIN, L_OMEGAMOMS_RRSGAIN )                     ! Outputs
      else
       CALL LRRS_RAMAN_OPS_MONO_2 &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS,             & ! Inputs
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,          & ! Inputs
          NLAYERS, NMOMENTS, NPOINTS_MONO, W_EXCIT,        & ! Inputs
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, LAYER_AIRCOLUMNS, & ! Inputs
          DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC,            & ! Inputs
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO,                & ! Inputs
          OMEGAMOMS_CABANNES, & ! Outputs
          OMEGAMOMS_RRSLOSS,  & ! Outputs
          OMEGAMOMS_RRSGAIN )   ! Outputs
      endif

!  Finish

      END SUBROUTINE L_SETUP_MASTER_MONO

!  End Module

      END MODULE lrrs_L_setup_master_m

