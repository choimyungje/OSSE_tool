
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
! #   INPUT SUBROUTINES :                                       #
! #                                                             #
! #     Initializing inputs:                                    #
! #                                                             #
! #          1. LRRS_INIT_INPUTS                                #
! #                                                             #
! #          2. LRRS_Sup_Init                                   #
! #          3. LRRS_BRDF_Sup_Init                              #
! #          4. LRRS_SLEAVE_Sup_Init                            #
! #          5. LRRS_SS_Sup_Init                                #
! #                                                             #
! #     Reading inputs from file:                               #
! #                                                             #
! #          6. LRRS_READ_INPUTS                                #
! #                                                             #
! #     These routines are called at start of the Main          #
! #       LRRS modules                                          #
! #                                                             #
! #          7. LRRS_CHECK_INPUT_DIMS                           #
! #          8. LRRS_CHECK_INPUTS                               #
! #          9. LRRS_DERIVE_INPUTS                              #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are:
!    (1) Bookkeeping improvements (use of "Only", clearer I/O specifications)
!    (2) Input for Taylor-series control has been expanded
!    (3) Specific BRDF inputs removed, replaced by general BRDF_SURFACE flag
!    (4) General surface-leaving flags introduced (replaces WATERLEAVING)
!    (5) Addition of LRRS_INIT_INPUTS subroutine to separately handle the
!        initializing of LRRS standard input variables (mick mod 7/20/2016)
!    (6) Addition of the four subroutines to initialize the LRRS supplement
!        standard input variables (mick mod 7/20/2016)

!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT, MAX_MOMENTS_INPUT checks
!   -- Rob mod 5/12/17 for 2p5a, Added PHASFUNC inputs

      MODULE lrrs_io_check_m

!      USE LRRS_PARS_m, Only : SDU

      PRIVATE
      PUBLIC :: LRRS_INIT_INPUTS,      &
                LRRS_Sup_Init,         &
                LRRS_BRDF_Sup_Init,    &
                LRRS_SLEAVE_Sup_Init,  &
                LRRS_SS_Sup_Init,      &     
                LRRS_READ_INPUTS,      &
                LRRS_CHECK_INPUT_DIMS, &
                LRRS_CHECK_INPUTS,     &
                LRRS_DERIVE_INPUTS

      CONTAINS

      SUBROUTINE LRRS_INIT_INPUTS &
        ( LRRS_FixIn, LRRS_ModIn )

!  Module of dimensions and numbers

      USE LRRS_PARS_m
      USE LRRS_IO_DEFS_m

!  Implicit none

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT-RRS standard input structures

      TYPE(LRRS_Fixed_Inputs), INTENT(OUT)    :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT(OUT) :: LRRS_ModIn

!  Initialize inputs
!  =================

!  LRRS Fixed Boolean
!  ------------------

!   -- Rob mod 5/12/17 for 2p5a, SSFULL renamed, SSCORR_NADIR added

      LRRS_FixIn%Bool%DO_RRS_OVERALL        = .FALSE.
      LRRS_FixIn%Bool%DO_DIRECTRRS_ONLY     = .FALSE.
      LRRS_FixIn%Bool%DO_DIRECT_BEAM        = .FALSE.

      LRRS_FixIn%Bool%DO_SSCORR_NADIR       = .FALSE.
      LRRS_FixIn%Bool%DO_SSCORR_OUTGOING    = .FALSE.
      LRRS_FixIn%Bool%DO_SSCORR_ALONE       = .FALSE.  ! renamed

      LRRS_FixIn%Bool%DO_MOLECSCAT_ONLY     = .FALSE.
      LRRS_FixIn%Bool%DO_PLANE_PARALLEL     = .FALSE.
      LRRS_FixIn%Bool%DO_DELTAM_SCALING     = .FALSE.
      LRRS_FixIn%Bool%DO_DOUBLE_CONVTEST    = .FALSE.
      LRRS_FixIn%Bool%DO_UPWELLING          = .FALSE.
      LRRS_FixIn%Bool%DO_DNWELLING          = .FALSE.

!mick chg 7/20/2016 - no partial layer output: hardwired for now 
      !LRRS_FixIn%Bool%DO_LBOUNDARIES        = .FALSE.
      LRRS_FixIn%Bool%DO_LBOUNDARIES        = .TRUE.

      LRRS_FixIn%Bool%DO_MVOUT_ONLY         = .FALSE.
      LRRS_FixIn%Bool%DO_ADDITIONAL_MVOUT   = .FALSE.

!  Surface flags. removed for LRRS Version 2.5
!    Control is now in the new BRDF/SLEAVE supplements
!      LRRS_FixIn%Bool%DO_LAMBERTIAN_SURFACE  = .FALSE.
!      LRRS_FixIn%Bool%DO_WATERLEAVING        = .FALSE.
!      LRRS_FixIn%Bool%DO_SHADOW_EFFECT       = .FALSE.
!      LRRS_FixIn%Bool%DO_GLITTER_DBMS        = .FALSE.
!      LRRS_FixIn%Bool%DO_LAMBERTIAN_FRACTION = .FALSE.

!  Surface control. Added, 9/9/15 for LRRS Version 2.5
!   - Modeled after the LIDORT 3.7 code
      LRRS_FixIn%Bool%DO_BRDF_SURFACE       = .FALSE.

!  Surface leaving Control. Added, 9/9/15 for LRRS Version 2.5
!   - Modeled after the LIDORT 3.7 code
      LRRS_FixIn%Bool%DO_SURFACE_LEAVING    = .FALSE.
      LRRS_FixIn%Bool%DO_SL_ISOTROPIC       = .FALSE.

      LRRS_FixIn%Bool%DO_LRRS_WRITEINPUT    = .FALSE.
      LRRS_FixIn%Bool%DO_LRRS_WRITERESULTS  = .FALSE.
      LRRS_FixIn%Bool%DO_LRRS_WRITESCENARIO = .FALSE.
      LRRS_FixIn%Bool%DO_LRRS_WRITEFOURIER  = .FALSE.

!  LRRS Modified Boolean
!  ---------------------

      LRRS_ModIn%MBool%DO_ELASTIC_ONLY       = .FALSE.
      LRRS_ModIn%MBool%DO_CABANNES_RAMAN     = .FALSE.
      LRRS_ModIn%MBool%DO_ENERGY_BALANCING   = .FALSE.
      LRRS_ModIn%MBool%DO_MSMODE_LRRS        = .FALSE.
      LRRS_ModIn%MBool%DO_BIN_REALIZATION    = .FALSE.
      LRRS_ModIn%MBool%DO_MONO_REALIZATION   = .FALSE.
      LRRS_ModIn%MBool%DO_USER_STREAMS       = .FALSE.
      LRRS_ModIn%MBool%DO_NO_AZIMUTH         = .FALSE.

!  LRRS Fixed Control
!  ------------------

      LRRS_FixIn%Cont%PROGRESS = 0

      LRRS_FixIn%Cont%NLAYERS        = 0
      LRRS_FixIn%Cont%NSTREAMS       = 0
!      LRRS_FixIn%Cont%NMOMENTS_INPUT = 0 !   -- Rob mod 5/12/17 for 2p5a, removed
      LRRS_FixIn%Cont%NFINELAYERS    = 0

      LRRS_FixIn%Cont%FLUX_FACTOR    = ONE
      LRRS_FixIn%Cont%SOLAR_ANGLE    = ZERO
      LRRS_FixIn%Cont%EARTH_RADIUS   = ZERO

      LRRS_FixIn%Cont%LRRS_ACCURACY     = ZERO
      LRRS_FixIn%Cont%LRRS_TAYLOR_ORDER = 0
      LRRS_FixIn%Cont%LRRS_TAYLOR_SMALL = ZERO

      LRRS_FixIn%Cont%DEBUG_FILENAMES = ' '

!  LRRS Fixed UserValues
!  ---------------------

!  Output options - viewer angles

      LRRS_FixIn%UserVal%N_USER_STREAMS = 0
      LRRS_FixIn%UserVal%USER_ANGLES    = ZERO

      LRRS_FixIn%UserVal%N_USER_RELAZMS = 0
      LRRS_FixIn%UserVal%USER_RELAZMS   = ZERO

!  Output options - layer output choices

      LRRS_FixIn%UserVal%N_LOUTPUT          = 0
      LRRS_FixIn%UserVal%LPARTIALS_OUTPUT   = ZERO
      LRRS_FixIn%UserVal%LBOUNDARIES_OUTPUT = 0

!  LRRS Fixed Spectral
!  -------------------

      LRRS_FixIn%Spect%LAMBDAS_RANKED = ZERO
      LRRS_FixIn%Spect%FLUXES_RANKED  = ONE

!  New Isotropic SLTERM. @@@ RobFix 06 sep 12
!     LRRS_FixIn%Spect%LRRS_SLTERM = ZERO

!  ### For monochromatic calculations ###

      LRRS_FixIn%Spect%NPOINTS_MONO = 0
      LRRS_FixIn%Spect%LAMBDA_EXCIT = ZERO
      LRRS_FixIn%Spect%W_EXCIT      = 0

!  ### For binning calculations ###

      LRRS_FixIn%Spect%NPOINTS_INNER = 0
      LRRS_FixIn%Spect%OFFSET_INNER  = 0
      LRRS_FixIn%Spect%NPOINTS_OUTER = 0

      LRRS_FixIn%Spect%BINLOWER = ZERO
      LRRS_FixIn%Spect%BINUPPER = ZERO

!  LRRS Fixed Atmosphere
!  ---------------------

      LRRS_FixIn%Atmos%HEIGHT_GRID = ZERO

      LRRS_FixIn%Atmos%LAYER_TEMPERATURES = ZERO
      LRRS_FixIn%Atmos%LAYER_AIRCOLUMNS   = ZERO

      LRRS_FixIn%Atmos%RAYLEIGH_XSEC  = ZERO
      LRRS_FixIn%Atmos%RAYLEIGH_DEPOL = ZERO

      LRRS_FixIn%Atmos%DELTAU_INPUT_UNSCALED      = ZERO
      LRRS_FixIn%Atmos%OMEGAMOMS_ELASTIC_UNSCALED = ZERO

!   -- Rob mod 5/12/17 for 2p5a, Added PHASFUNC inputs

      LRRS_FixIn%Atmos%OMEGAPHASFUNC_ELASTIC_UP = ZERO
      LRRS_FixIn%Atmos%OMEGAPHASFUNC_ELASTIC_DN = ZERO

!  LRRS Modified Atmosphere
!  ------------------------

      LRRS_ModIn%MAtmos%GEOMETRY_SPECHEIGHT = ZERO

!  LRRS Fixed Surface
!  ------------------

!  All the BRDF stuff has been removed for LRRS Version 2.5. 9/9/15
!    All BRDF control now in the BRDF supplement....
!      LRRS_FixIn%Surf%LAMBERTIAN_FRACTION = ZERO
!      LRRS_FixIn%Surf%NSTREAMS_BRDF = 0
!      LRRS_FixIn%Surf%WHICH_BRDF    = 0
!      LRRS_FixIn%Surf%BRDF_FACTOR   = ZERO
!      LRRS_FixIn%Surf%BRDF_NPARS    = 0
!      LRRS_FixIn%Surf%BRDF_PARS     = ZERO
!      LRRS_FixIn%Surf%BRDF_NAMES    = ' '

!  Ranked surface albedos

      LRRS_FixIn%Surf%ALBEDOS_RANKED = ZERO

!  Finish

      END SUBROUTINE LRRS_INIT_INPUTS

!

      SUBROUTINE LRRS_Sup_Init ( LRRS_Sup )

      USE LRRS_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LRRS supplement input structure

      TYPE(LRRS_Sup_InOut), INTENT(INOUT) :: LRRS_Sup

!  Initialize LRRS supplement inputs
!  =================================

      CALL LRRS_BRDF_Sup_Init   ( LRRS_Sup )
      CALL LRRS_SLEAVE_Sup_Init ( LRRS_Sup )
      CALL LRRS_SS_Sup_Init     ( LRRS_Sup )

!  Finish

      END SUBROUTINE LRRS_Sup_Init

!

      SUBROUTINE LRRS_BRDF_Sup_Init ( LRRS_Sup )

      USE LRRS_PARS_m, Only : ZERO
      USE LRRS_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LRRS supplement input structure

      TYPE(LRRS_Sup_InOut), INTENT(INOUT) :: LRRS_Sup

!  Initialize LRRS brdf supplement inputs
!  ======================================

      LRRS_Sup%BRDF%DO_BRDF_Wav1    = .TRUE.
      LRRS_Sup%BRDF%EXACTDB_BRDFUNC = ZERO
      LRRS_Sup%BRDF%BRDF_F_0        = ZERO
      LRRS_Sup%BRDF%BRDF_F          = ZERO
      LRRS_Sup%BRDF%USER_BRDF_F_0   = ZERO
      LRRS_Sup%BRDF%USER_BRDF_F     = ZERO

!  Finish

      END SUBROUTINE LRRS_BRDF_Sup_Init

!

      SUBROUTINE LRRS_SLEAVE_Sup_Init ( LRRS_Sup )

      USE LRRS_PARS_m, Only : ZERO
      USE LRRS_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LRRS supplement input structure

      TYPE(LRRS_Sup_InOut), INTENT(INOUT) :: LRRS_Sup

!  Initialize LRRS sleave supplement inputs
!  ========================================

      LRRS_Sup%SLEAVE%DO_SLEAVE_Wav1    = .TRUE.
      LRRS_Sup%SLEAVE%SLTERM_ISOTROPIC  = ZERO
      LRRS_Sup%SLEAVE%SLTERM_USERANGLES = ZERO
      LRRS_Sup%SLEAVE%SLTERM_F_0        = ZERO
      LRRS_Sup%SLEAVE%USER_SLTERM_F_0   = ZERO

!  Finish

      END SUBROUTINE LRRS_SLEAVE_Sup_Init

!

      SUBROUTINE LRRS_SS_Sup_Init ( LRRS_Sup )

      USE LRRS_PARS_m, Only : ZERO
      USE LRRS_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LRRS supplement input structure

      TYPE(LRRS_Sup_InOut), INTENT(INOUT) :: LRRS_Sup

!  Initialize LRRS single-scatter supplement inputs
!  ================================================

      LRRS_Sup%SS%ELASTIC_SS_UP = ZERO
      LRRS_Sup%SS%ELASTIC_SS_DN = ZERO
      LRRS_Sup%SS%ELASTIC_DB_UP = ZERO
      LRRS_Sup%SS%RAMAN_SS_UP   = ZERO
      LRRS_Sup%SS%RAMAN_SS_DN   = ZERO
      LRRS_Sup%SS%RAMAN_DB_UP   = ZERO

!  Finish

      END SUBROUTINE LRRS_SS_Sup_Init

!

      SUBROUTINE LRRS_READ_INPUTS &
        ( FILNAM,                 &
          LRRS_FixIn, LRRS_ModIn, &
          STATUS, NMESS, MESS )

!mick mod 7/20/2016 - modified subroutine header description and added list

!  This subroutine initializes all LRRS inputs and reads in most control inputs.
!  It does not read in certain control inputs or optical properties.  Specifically,
!  the following lists the LRRS inputs NOT defined by this subroutine (or only
!  CONDITIONALLY defined - as denoted by a "C").  It is put here for handy reference.

!   -- Rob mod 5/12/17 for 2p5a, List modified --> removed NMOMENTS_INPUT, added PHASFUNC inputs

!  LRRS Fixed Boolean Inputs:
!     LRRS_FixIn%Bool%DO_SSFULL        C
!     LRRS_FixIn%Bool%DO_SL_ISOTROPIC  C
!  LRRS Fixed Control Inputs:
!     LRRS_FixIn%Cont%PROGRESS
!     LRRS_FixIn%Cont%NLAYERS
!     LRRS_FixIn%Cont%NFINELAYERS      C
!  LRRS Fixed Spectral:
!     LRRS_FixIn%Spect%LAMBDAS_RANKED
!     LRRS_FixIn%Spect%FLUXES_RANKED
!     **Mono inputs**
!     LRRS_FixIn%Spect%NPOINTS_MONO
!     LRRS_FixIn%Spect%LAMBDA_EXCIT
!     LRRS_FixIn%Spect%W_EXCIT
!     **Bin inputs**
!     LRRS_FixIn%Spect%NPOINTS_INNER
!     LRRS_FixIn%Spect%OFFSET_INNER
!     LRRS_FixIn%Spect%NPOINTS_OUTER
!     LRRS_FixIn%Spect%BINLOWER
!     LRRS_FixIn%Spect%BINUPPER                   
!  LRRS Fixed Atmosphere:
!     LRRS_FixIn%Atmos%HEIGHT_GRID
!     LRRS_FixIn%Atmos%LAYER_TEMPERATURES
!     LRRS_FixIn%Atmos%LAYER_AIRCOLUMNS
!     LRRS_FixIn%Atmos%RAYLEIGH_XSEC
!     LRRS_FixIn%Atmos%RAYLEIGH_DEPOL
!     LRRS_FixIn%Atmos%DELTAU_INPUT_UNSCALED
!     LRRS_FixIn%Atmos%OMEGAMOMS_ELASTIC_UNSCALED
!      LRRS_FixIn%Atmos%OMEGAPHASFUNC_ELASTIC_UP
!      LRRS_FixIn%Atmos%OMEGAPHASFUNC_ELASTIC_DN
!  LRRS Modified Atmosphere:
!     LRRS_ModIn%MAtmos%GEOMETRY_SPECHEIGHT
!  LRRS Fixed Surface:
!     LRRS_FixIn%Surf%ALBEDOS_RANKED

!  Used modules

      USE LRRS_PARS_m
      USE LRRS_INPUTS_DEF_m
      USE LRRS_AUX2_m , Only : GFINDPAR, LRRS_ERROR_TRACE, LEN_STRING

      IMPLICIT NONE

!  Module input arguments (input filename, error filename)
!  -------------------------------------------------------

      CHARACTER (LEN=*), INTENT(IN) ::  FILNAM

!mick mod 7/20/2016 - removed input of error file to promote more consistency with LIDORT 
      !CHARACTER (LEN=*), INTENT(IN) ::  EFILNAM

!  Module output arguments (read-in from file)
!  -------------------------------------------

!  LIDORT-RRS input structures

      TYPE(LRRS_Fixed_Inputs), INTENT(OUT)    :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT(OUT) :: LRRS_ModIn

!  Overall status

      INTEGER, INTENT(OUT) ::             STATUS
      INTEGER, INTENT(OUT) ::             NMESS
      CHARACTER (LEN=120), INTENT(OUT) :: MESS ( MAX_MESSAGES )

!  Local variables
!  ===============

!mick fix 7/20/2016 - re-organized declaration & initialization of inputs
!                   - insured all inputs now at least initialized 

!  LRRS Fixed Boolean
!  ------------------

!  Top level control

      LOGICAL :: DO_RRS_OVERALL
      LOGICAL :: DO_DIRECTRRS_ONLY

!  Direct beam flag

      LOGICAL :: DO_DIRECT_BEAM

!  Single scattering flags
!   -- Rob mod 5/12/17 for 2p5a, SSFULL renamed, SSCORR_NADIR added

      LOGICAL :: DO_SSCORR_NADIR
      LOGICAL :: DO_SSCORR_OUTGOING
      LOGICAL :: DO_SSCORR_ALONE     ! renamed

!  Molecular scattering only (no aerosols/clouds)

      LOGICAL :: DO_MOLECSCAT_ONLY

!  Plane parallel flag

      LOGICAL :: DO_PLANE_PARALLEL

!  Delta-M scaling (only with aerosols)

      LOGICAL :: DO_DELTAM_SCALING

!  Double convergence test (only with aerosols)

      LOGICAL :: DO_DOUBLE_CONVTEST

!  Upwelling/downwelling flags

      LOGICAL :: DO_UPWELLING
      LOGICAL :: DO_DNWELLING

!  Level boundaries flag

      LOGICAL :: DO_LBOUNDARIES

!  Mean-value output options

      LOGICAL :: DO_MVOUT_ONLY
      LOGICAL :: DO_ADDITIONAL_MVOUT

!  Surface flags. removed for LRRS Version 2.5
!    Control is now in the new BRDF/SLEAVE supplements
!      LOGICAL :: DO_LAMBERTIAN_SURFACE
!      LOGICAL :: DO_WATERLEAVING
!      LOGICAL :: DO_SHADOW_EFFECT
!      LOGICAL :: DO_GLITTER_DBMS
!      LOGICAL :: DO_LAMBERTIAN_FRACTION

!  Surface control. Added, 9/9/15 for LRRS Version 2.5
!   - Modeled after the LIDORT 3.7 code

      LOGICAL :: DO_BRDF_SURFACE

!  Surface leaving Control. Added, 9/9/15 for LRRS Version 2.5
!   - Modeled after the LIDORT 3.7 code

      LOGICAL :: DO_SURFACE_LEAVING
      LOGICAL :: DO_SL_ISOTROPIC

!  Debug write flags

      LOGICAL :: DO_LRRS_WRITEINPUT
      LOGICAL :: DO_LRRS_WRITESCENARIO
      LOGICAL :: DO_LRRS_WRITEFOURIER
      LOGICAL :: DO_LRRS_WRITERESULTS

!  LRRS Modified Boolean
!  ---------------------

!  Top level control

      LOGICAL :: DO_ELASTIC_ONLY

!  Solution method (options are mutually exclusive)

      LOGICAL :: DO_CABANNES_RAMAN
      LOGICAL :: DO_ENERGY_BALANCING

!  MS mode flag

      LOGICAL :: DO_MSMODE_LRRS

!  Binning or monochromatic realization
!    (options are mutually exclusive)

      LOGICAL :: DO_BIN_REALIZATION
      LOGICAL :: DO_MONO_REALIZATION

!  User angles flag

      LOGICAL :: DO_USER_STREAMS

!  No azimuth flag

      LOGICAL :: DO_NO_AZIMUTH

!  LRRS Fixed Control
!  ------------------

!  Progress control (If set to 0, no screen output)

      INTEGER :: PROGRESS

!  Number of layers

      INTEGER :: NLAYERS

!  Number of discrete ordinate streams (in the half-space)

      INTEGER :: NSTREAMS

!  number of moments in the single scatter calculation
!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT
!      INTEGER :: NMOMENTS_INPUT

!  Number of fine layers for sscorr outgoing

      INTEGER :: NFINELAYERS

!  Flux factor

      REAL(FPK) :: FLUX_FACTOR

!  Solar beam input geometry

      REAL(FPK) :: SOLAR_ANGLE

!  Earth Radius

      REAL(FPK) :: EARTH_RADIUS

!  Overall accuracy for convergence criterion

      REAL(FPK) :: LRRS_ACCURACY

!  Small Number control.
!   Taylor ordering parameter, added 9/9/15. Should be set to 2 or 3

      INTEGER   :: LRRS_TAYLOR_ORDER
      REAL(FPK) :: LRRS_TAYLOR_SMALL

!  Filenames for output

      CHARACTER (LEN=60), DIMENSION(4) :: DEBUG_FILENAMES

!  LRRS Fixed UserValues
!  ---------------------

!  User stream variables

      INTEGER :: N_USER_STREAMS
      REAL(FPK), DIMENSION (MAX_USER_STREAMS) :: USER_ANGLES

!  User azimuth variables

      INTEGER :: N_USER_RELAZMS
      REAL(FPK), DIMENSION (MAX_USER_RELAZMS) :: USER_RELAZMS

!  Number of level output choices (all)

      INTEGER :: N_LOUTPUT

!  Off-boundary level choices

      REAL(FPK), DIMENSION (MAX_LOUTPUT) :: LPARTIALS_OUTPUT

!  Layer boundary choices

      INTEGER, DIMENSION (MAX_LOUTPUT) :: LBOUNDARIES_OUTPUT

!  LRRS Fixed Spectral
!  -------------------

!  Ranked wavelengths

      REAL(FPK), DIMENSION (MAX_POINTS) :: LAMBDAS_RANKED

!  Ranked fluxes

      REAL(FPK), DIMENSION (MAX_POINTS) :: FLUXES_RANKED

!  New Isotropic SLTERM. @@@ RobFix 06 sep 12
!      REAL(FPK), DIMENSION (MAX_POINTS) :: LRRS_SLTERM


!  ### For monochromatic calculations ###

!  Number of monochromatic points

      INTEGER :: NPOINTS_MONO

!  Excitation wavelength

      REAL(FPK) :: LAMBDA_EXCIT

!  Excitation index

      INTEGER :: W_EXCIT


!  ### For binning calculations ###

!  Number of inner-window points

      INTEGER :: NPOINTS_INNER

!  Offset for inner-window range

      INTEGER :: OFFSET_INNER

!  Number of outer-window points

      INTEGER :: NPOINTS_OUTER

!  Bin upper and lower limits

      REAL(FPK), DIMENSION (MAX_POINTS) :: BINLOWER
      REAL(FPK), DIMENSION (MAX_POINTS) :: BINUPPER

!  LRRS Fixed Atmosphere
!  ---------------------

!  Height grid

      REAL(FPK), DIMENSION (0:MAX_LAYERS) :: HEIGHT_GRID

!  Input layer temperatures, must be in deg K

      REAL(FPK), DIMENSION (MAX_LAYERS) :: LAYER_TEMPERATURES

!  Input layer Air columns, should be in mol/cm^2 or [DU]

      REAL(FPK), DIMENSION (MAX_LAYERS) :: LAYER_AIRCOLUMNS

!  Rayleigh cross-sections

      REAL(FPK), DIMENSION (MAX_POINTS) :: RAYLEIGH_XSEC

!  Rayleigh depolarization ratios

      REAL(FPK), DIMENSION (MAX_POINTS) :: RAYLEIGH_DEPOL

!  Unscaled optical depth (elastic input)

      REAL(FPK), DIMENSION (MAX_LAYERS, MAX_POINTS) :: DELTAU_INPUT_UNSCALED

!  Product of unscaled single-scatter albedo and phase function moments (elastic input)
!   -- Rob mod 5/12/17 for 2p5a, Changed MAX_MOMENTS dimensioning

      REAL(FPK), DIMENSION (MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS) :: OMEGAMOMS_ELASTIC_UNSCALED

!   -- Rob mod 5/12/17 for 2p5a, added PHASFUNC inputs

      REAL(FPK), DIMENSION (MAX_LAYERS, 0:MAX_GEOMETRIES, MAX_POINTS) :: OMEGAPHASFUNC_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_LAYERS, 0:MAX_GEOMETRIES, MAX_POINTS) :: OMEGAPHASFUNC_ELASTIC_DN

!  LRRS Modified Atmosphere
!  ------------------------

!  A height used when using the outgoing sphericity correction

      REAL(FPK) :: GEOMETRY_SPECHEIGHT

!  LRRS Fixed Surface
!  ------------------

!  All the BRDF stuff has been removed for LRRS Version 2.5. 9/9/15
!    All BRDF control now in the BRDF supplement....
!  Fraction of surface treated as Lambertian
!      REAL(FPK) :: LAMBERTIAN_FRACTION
!  Number of BRDF azimuth streams
!      INTEGER :: NSTREAMS_BRDF
!  BRDF name
!      CHARACTER (LEN=10) :: BRDF_NAMES
!  LRRS index of BRDF
!      INTEGER :: WHICH_BRDF
!  BRDF weighting factor
!      REAL(FPK) :: BRDF_FACTOR
!  Number of parameters associated with BRDF
!      INTEGER :: BRDF_NPARS
!  BRDF parameters
!      REAL(FPK), DIMENSION (MAX_BRDF_PARAMETERS) :: BRDF_PARS

!  Ranked surface albedos

      REAL(FPK), DIMENSION (MAX_POINTS) :: ALBEDOS_RANKED

!  Helper variables
!  ================

      CHARACTER (LEN=6), PARAMETER :: PREFIX = 'LRRS -'

      LOGICAL ::   ERROR
      CHARACTER (LEN=80) :: PAR_STR
      !LOGICAL ::   GFINDPAR
      INTEGER ::   I, FILUNIT
      !INTEGER ::   LEN_STRING
      CHARACTER (LEN=120) :: MAIL! , ACTION
      INTEGER ::   N

!  Initialize status

      STATUS = LRRS_SUCCESS

!  Initialize error message count

      NMESS = 0
      DO N = 1, MAX_MESSAGES
        MESS(N) = ' '
      ENDDO

!  Initialize inputs
!  =================

      CALL LRRS_INIT_INPUTS ( LRRS_FixIn, LRRS_ModIn )

!  Read section
!  ============

!  Open file

      FILUNIT = LRRS_INUNIT
      OPEN(LRRS_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  ALGORITHM CONTROL VARIABLES
!  ===========================

!  Overall filling flag

      PAR_STR = 'Do Filling calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_RRS_OVERALL
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Filling calculation for single scatter only

      PAR_STR = 'Do RRS only with primary scatter?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_DIRECTRRS_ONLY
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Elastic-only flag

      DO_ELASTIC_ONLY = .NOT.DO_RRS_OVERALL .AND. .NOT.DO_DIRECTRRS_ONLY

!  Realization control.  Binning or Monochromatic.
!  options are mutually exclusive
!mick fix 7/20/2016 - removed if condition (need to always read these)

      !IF ( .NOT. DO_ELASTIC_ONLY ) THEN
        PAR_STR = 'Do Inelastic scattering with binning realization?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_BIN_REALIZATION
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

        PAR_STR = &
            'Do Inelastic scattering with monochromatic realization?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_MONO_REALIZATION
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      !ENDIF

!  Solution method. Use Cabannes-Raman or Energy-balancing
!  options are mutually exclusive (checked here). One must be true!

      PAR_STR = 'Do Energy-balancing treatment?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_ENERGY_BALANCING
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      DO_CABANNES_RAMAN = .not. DO_ENERGY_BALANCING

!  Pseudo-spherical control

      PAR_STR = 'Plane-parallel treatment of direct beam?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_PLANE_PARALLEL
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Scatterers and phase function control

      PAR_STR='Molecular scattering atmosphere only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_MOLECSCAT_ONLY
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Multiple scatter source function output control
!   Not enabled in this version
!mick chg 7/20/2016 - initialized with other modified Booleans now

      !DO_MSMODE_LRRS = .FALSE.  

!  Solar beam control

      PAR_STR = 'Include direct beam?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_DIRECT_BEAM
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Outgoing Single scatter correction of whole field
!    ENABLED 29 May 2008.

      PAR_STR = 'Single scatter outgoing correction of all fields?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SSCORR_OUTGOING
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Nadir Single scatter correction of whole field
!   -- Rob mod 5/12/17 for 2p5a,  SSCORR_NADIR added

      PAR_STR = 'Single scatter nadir correction of all fields?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SSCORR_NADIR
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  SSFULL renamed, optional for NADIR as well as OUTGOING
!   -- Rob mod 5/12/17 for 2p5a

      IF ( DO_SSCORR_OUTGOING .or. DO_SSCORR_NADIR ) THEN
        PAR_STR = 'Single scatter correction alone (no MS)?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SSCORR_ALONE
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  Double convergence test

      PAR_STR = 'Perform double convergence test?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_DOUBLE_CONVTEST
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  No azimuth dependence (TRUE means Fourier m = 0 only )

      PAR_STR = 'No azimuth dependence in the solution?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_NO_AZIMUTH
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Delta-M scaling

      PAR_STR = 'Include Delta-M scaling?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_DELTAM_SCALING
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Accuracy input for Fourier series convergence

      PAR_STR = 'Fourier series convergence'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) LRRS_ACCURACY
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Small Number control.
!   Taylor ordering parameter, added 9/9/15. Should be set to 2 or 3
!   Taylor smallness parameter, Suggested value = 1.0D-03

      PAR_STR = 'Small number parameter for Taylor series expansions'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) LRRS_TAYLOR_SMALL
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'Ordering parameter for Taylor series expansions'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) LRRS_TAYLOR_ORDER
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Earth radius

      PAR_STR = 'Earth radius'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) EARTH_RADIUS
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Flux factor

      PAR_STR = 'Flux factor'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FLUX_FACTOR
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Number of fine layers for single scatter outgoing correction

      IF ( DO_SSCORR_OUTGOING ) THEN
        PAR_STR = 'Number of fine layers for single scatter outgoing'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NFINELAYERS
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  2, GEOMETRY AND OUTPUT VARIABLES
!  ================================

!  Directional output control

      PAR_STR = 'Upwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_UPWELLING
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'Downwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_DNWELLING
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  User-defined Stream angle output control

      PAR_STR = 'User-defined stream angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Flag for output for level boundaries
!    - If not set, then output at off-boundary optical depths

      PAR_STR = 'Output at layer boundaries?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) DO_LBOUNDARIES
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      IF ( .NOT. DO_LBOUNDARIES ) THEN
        MAIL   = 'Partial layer output not enabled in LRRS yet'
        NMESS  = NMESS + 1
        MESS(NMESS) = MAIL(1:LEN_STRING(MAIL))
        GO TO 455
      ENDIF

!  Mean-value output control

      PAR_STR = 'Generate only mean value output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_MVOUT_ONLY
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'Generate mean value output additionally?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_ADDITIONAL_MVOUT
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Number of  discrete ordinate streams
!    -- Number of streams is checked now, to see if exceeds dimensioning

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      IF ( NSTREAMS .GT. MAX_STREAMS ) THEN
        MAIL   = 'Number of half-space streams > maximum dimension'
        CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
        GO TO 455
      ENDIF

!  BOA solar zenith angle input
!  9/11/15. LRRS Version 2.5. Changed name from TOA to BOA., otherwise the same

      PAR_STR = 'BOA solar zenith angle (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) SOLAR_ANGLE
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  User defined stream angles (should be positive from 0 to 90)
!    -- Number of angles is checked now, to see if exceeds dimensioning

      IF ( DO_USER_STREAMS ) THEN
        PAR_STR = 'Number of user-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          MAIL   = 'Number of user streams > maximum dimension'
          CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
          GO TO 455
        ENDIF
        PAR_STR = 'User-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
          ENDDO
        ENDIF
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  Relative azimuth input angles if Azimuth flag set
!    -- Number of angles is checked now, to see if exceeds dimensioning
!    -- Set default to 1 value of 0.0 if the flag is not set

      IF ( .NOT. DO_NO_AZIMUTH ) THEN
        PAR_STR = 'Number of user-defined relative azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
          MAIL   = 'Number of relative azimuths > maximum dimension'
          CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
          GO TO 455
        ENDIF
        PAR_STR = 'User-defined relative azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_RELAZMS
            READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
          ENDDO
        ENDIF
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ELSE
        N_USER_RELAZMS = 1
        USER_RELAZMS(1) = ZERO
      ENDIF

!  User defined off-boundary partial layer output choices
!    -- These are specified as follows
!         LPARTIALS = 17.49 means that output will be in Layer 18
!         and .49 of the way down from the top (in height)
!    -- Number is checked now, to see if exceeds dimensioning

      IF ( .NOT.DO_LBOUNDARIES ) THEN
        PAR_STR = 'Number of partial layer output choices'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) N_LOUTPUT
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        IF ( N_LOUTPUT .GT. MAX_LOUTPUT ) THEN
          MAIL   ='Number of output choices > maximum dimension'
          CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
          GO TO 455
        ENDIF
        PAR_STR = 'Partial layer output choices'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_LOUTPUT
            READ (FILUNIT,*,ERR=998) LPARTIALS_OUTPUT(I)
          ENDDO
        ENDIF
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  User defined LAYER boundary optical depths (INTEGER :: input)
!    -- Number is checked now, to see if exceeds dimensioning

      IF ( DO_LBOUNDARIES ) THEN
        PAR_STR = 'Number of layer boundary output choices'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) N_LOUTPUT
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        IF ( N_LOUTPUT .GT. MAX_LOUTPUT ) THEN
          MAIL   ='Number of output choices > maximum dimension'
          CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
          GO TO 455
        ENDIF
        PAR_STR = 'Layer boundaries for output (0 = TOA)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_LOUTPUT
            READ (FILUNIT,*,ERR=998) LBOUNDARIES_OUTPUT(I)
          ENDDO
        ENDIF
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  4. SURFACE VARIABLES
!  ====================

!  BRDF surface. 9/11/15. replaces Lambertian flag, Version 2.5

      PAR_STR = 'BRDF surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_BRDF_SURFACE
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Surface leaving. 9/11/15. New, Version 2.5

      PAR_STR = 'Surface leaving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SURFACE_LEAVING
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      IF ( DO_SURFACE_LEAVING ) THEN
        PAR_STR = 'Isotropic Surface leaving?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SL_ISOTROPIC
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  BRDF input. 9/11/15, removed VErsion 2.5, now part of supplement.

!      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
!  Lambertian fraction to be included ?
!        PAR_STR = 'Do Lambertian component as part of BRDF surface?'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!             READ (FILUNIT,*,ERR=998) DO_LAMBERTIAN_FRACTION
!        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
!  Set waterleaving flag
!        IF ( DO_LAMBERTIAN_FRACTION ) THEN
!          PAR_STR = 'Do water-leaving Lambertian contribution?'
!          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!               READ (FILUNIT,*,ERR=998) DO_WATERLEAVING
!          CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
!        ENDIF
!  Set lambertian fraction
!        IF ( DO_LAMBERTIAN_FRACTION.AND..NOT.DO_WATERLEAVING ) THEN
!          PAR_STR = 'Lambertian fraction as part of BRDF surface'
!          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!               READ (FILUNIT,*,ERR=998) LAMBERTIAN_FRACTION
!          CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
!        ENDIF
!!  Direct beam correction. Now automatic, subsumed in the SS calculation
!!        PAR_STR = 'Do Direct beam correction for BRDF surface?'
!!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
!!     &       READ (FILUNIT,*,ERR=998) DO_DB_CORRECTION
!!        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
!  Number of BRDF azimuth streams, check this value
!        PAR_STR = 'Number of bidirectional reflectance streams'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!             READ (FILUNIT,*,ERR=998) NSTREAMS_BRDF
!        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
!        IF ( NSTREAMS_BRDF .GT. MAX_STREAMS_BRDF ) THEN
!          MAIL = 'Number of BRDF streams > maximum dimension.'
!          ACTION = ' Re-set input value.'
!          STATUS = LRRS_SERIOUS
!! mick fix 3/22/11 - replaced call statement with following two lines
!          !CALL LRRS_ERROR_TRACE ( MAIL, ACTION, STATUS )
!          NMESS = NMESS + 1
!          MESS(NMESS) = &
!            MAIL(1:LEN_STRING(MAIL))//ACTION(1:LEN_STRING(ACTION))
!          RETURN
!        ENDIF
!  Main kernel input
!        PAR_STR = &'Kernel name, index, amplitude, # parameters, parameters'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!          READ (FILUNIT,56,ERR=998) &
!               BRDF_NAMES, WHICH_BRDF, BRDF_FACTOR, BRDF_NPARS,(BRDF_PARS(K),K=1,3)
!        ENDIF
!        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
! 56     FORMAT( A10, I2, F6.2, I2, 3F12.6 )
!  Shadowing input (for Cox-Munk type)
!        IF ( BRDF_NAMES .EQ. 'Cox-Munk  ' ) THEN
!          PAR_STR = 'Do shadow effect for glitter kernels?'
!          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!            READ (FILUNIT,*,ERR=998)DO_SHADOW_EFFECT
!          ENDIF
!          CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
!        ENDIF
!  Multiple reflectance DB correction (for Cox-Munk type)
!        IF ( BRDF_NAMES .EQ. 'Cox-Munk  ' ) THEN
!          PAR_STR = 'Do multiple reflectance for glitter kernels?'
!          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!            READ (FILUNIT,*,ERR=998)DO_GLITTER_DBMS
!          ENDIF
!          CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
!        ENDIF
!      ENDIF

!  5. DEBUG FLAGS and FILE NAMES
!  =============================

!  Output write flags

      PAR_STR = 'Input control write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_LRRS_WRITEINPUT
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'Input scenario write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_LRRS_WRITESCENARIO
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'Fourier component output write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_LRRS_WRITEFOURIER
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'Results write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_LRRS_WRITERESULTS
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Output filenames

      PAR_STR = 'filename for input write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,'(a)',ERR=998) DEBUG_FILENAMES(1)
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'filename for scenario write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,'(a)',ERR=998) DEBUG_FILENAMES(2)
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'filename for Fourier output'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,'(a)',ERR=998) DEBUG_FILENAMES(3)
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'filename for main output'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,'(a)',ERR=998) DEBUG_FILENAMES(4)
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      CLOSE(FILUNIT)

!  Check number of messages at this point

      if ( nmess .gt. 0 ) GO TO 455

!  Define LIDORT-RRS inputs
!  ------------------------

!mick mod 7/20/2016 - initializing of all LRRS type structure input variables now done in
!                     subroutine LRRS_INIT_INPUTS; only those read in from the LRRS config
!                     file possibly modified here.  If conditions applied where appropriate.

!  Fixed Boolean Inputs

      LRRS_FixIn%Bool%DO_RRS_OVERALL         = DO_RRS_OVERALL
      LRRS_FixIn%Bool%DO_DIRECTRRS_ONLY      = DO_DIRECTRRS_ONLY
      LRRS_FixIn%Bool%DO_DIRECT_BEAM         = DO_DIRECT_BEAM

!   -- Rob mod 5/12/17 for 2p5a,  SSCORR_NADIR added, SSFULL renamed

      LRRS_FixIn%Bool%DO_SSCORR_NADIR        = DO_SSCORR_NADIR
      LRRS_FixIn%Bool%DO_SSCORR_OUTGOING     = DO_SSCORR_OUTGOING
      IF ( DO_SSCORR_OUTGOING .or. DO_SSCORR_NADIR ) then
          LRRS_FixIn%Bool%DO_SSCORR_ALONE = DO_SSCORR_ALONE
      endif

      LRRS_FixIn%Bool%DO_MOLECSCAT_ONLY      = DO_MOLECSCAT_ONLY
      LRRS_FixIn%Bool%DO_PLANE_PARALLEL      = DO_PLANE_PARALLEL
      LRRS_FixIn%Bool%DO_DELTAM_SCALING      = DO_DELTAM_SCALING
      LRRS_FixIn%Bool%DO_DOUBLE_CONVTEST     = DO_DOUBLE_CONVTEST
      LRRS_FixIn%Bool%DO_UPWELLING           = DO_UPWELLING
      LRRS_FixIn%Bool%DO_DNWELLING           = DO_DNWELLING
      LRRS_FixIn%Bool%DO_LBOUNDARIES         = DO_LBOUNDARIES
      LRRS_FixIn%Bool%DO_MVOUT_ONLY          = DO_MVOUT_ONLY
      LRRS_FixIn%Bool%DO_ADDITIONAL_MVOUT    = DO_ADDITIONAL_MVOUT

! 9/11/15. Version 2.5, replacement by general inputs
!      LRRS_FixIn%Bool%DO_LAMBERTIAN_SURFACE  = DO_LAMBERTIAN_SURFACE
!      LRRS_FixIn%Bool%DO_WATERLEAVING        = DO_WATERLEAVING
!      LRRS_FixIn%Bool%DO_SHADOW_EFFECT       = DO_SHADOW_EFFECT
!      LRRS_FixIn%Bool%DO_GLITTER_DBMS        = DO_GLITTER_DBMS
!      LRRS_FixIn%Bool%DO_LAMBERTIAN_FRACTION = DO_LAMBERTIAN_FRACTION
      LRRS_FixIn%Bool%DO_BRDF_SURFACE    = DO_BRDF_SURFACE
      LRRS_FixIn%Bool%DO_SURFACE_LEAVING = DO_SURFACE_LEAVING
      IF ( DO_SURFACE_LEAVING ) LRRS_FixIn%Bool%DO_SL_ISOTROPIC = DO_SL_ISOTROPIC

      LRRS_FixIn%Bool%DO_LRRS_WRITEINPUT     = DO_LRRS_WRITEINPUT
      LRRS_FixIn%Bool%DO_LRRS_WRITESCENARIO  = DO_LRRS_WRITESCENARIO
      LRRS_FixIn%Bool%DO_LRRS_WRITEFOURIER   = DO_LRRS_WRITEFOURIER
      LRRS_FixIn%Bool%DO_LRRS_WRITERESULTS   = DO_LRRS_WRITERESULTS

!  Modified Boolean Inputs

      LRRS_ModIn%MBool%DO_ELASTIC_ONLY     = DO_ELASTIC_ONLY
      LRRS_ModIn%MBool%DO_CABANNES_RAMAN   = DO_CABANNES_RAMAN
      LRRS_ModIn%MBool%DO_ENERGY_BALANCING = DO_ENERGY_BALANCING
!mick chg 7/20/2016 - Not enabled in this version / initialized with other modified Booleans
      !LRRS_ModIn%MBool%DO_MSMODE_LRRS      = DO_MSMODE_LRRS
      LRRS_ModIn%MBool%DO_BIN_REALIZATION  = DO_BIN_REALIZATION
      LRRS_ModIn%MBool%DO_MONO_REALIZATION = DO_MONO_REALIZATION
      LRRS_ModIn%MBool%DO_USER_STREAMS     = DO_USER_STREAMS
      LRRS_ModIn%MBool%DO_NO_AZIMUTH       = DO_NO_AZIMUTH

!  Fixed Control Inputs. 9/11/15. Version 2.5, expand Taylor variables

      LRRS_FixIn%Cont%NSTREAMS              = NSTREAMS
      IF ( DO_SSCORR_OUTGOING ) LRRS_FixIn%Cont%NFINELAYERS = NFINELAYERS
      LRRS_FixIn%Cont%FLUX_FACTOR           = FLUX_FACTOR
      LRRS_FixIn%Cont%SOLAR_ANGLE           = SOLAR_ANGLE
      LRRS_FixIn%Cont%EARTH_RADIUS          = EARTH_RADIUS
      LRRS_FixIn%Cont%LRRS_ACCURACY         = LRRS_ACCURACY
      LRRS_FixIn%Cont%LRRS_TAYLOR_ORDER     = LRRS_TAYLOR_ORDER
      LRRS_FixIn%Cont%LRRS_TAYLOR_SMALL     = LRRS_TAYLOR_SMALL
      LRRS_FixIn%Cont%DEBUG_FILENAMES       = DEBUG_FILENAMES

!  Fixed User-Value Inputs

      IF ( DO_USER_STREAMS ) THEN
         LRRS_FixIn%UserVal%N_USER_STREAMS = N_USER_STREAMS
         LRRS_FixIn%UserVal%USER_ANGLES(1:N_USER_STREAMS)   = USER_ANGLES(1:N_USER_STREAMS)
      ENDIF
      LRRS_FixIn%UserVal%N_USER_RELAZMS    = N_USER_RELAZMS
      LRRS_FixIn%UserVal%USER_RELAZMS(1:N_USER_RELAZMS)     = USER_RELAZMS(1:N_USER_RELAZMS)
      LRRS_FixIn%UserVal%N_LOUTPUT         = N_LOUTPUT
      IF ( .NOT.DO_LBOUNDARIES ) THEN
         LRRS_FixIn%UserVal%LPARTIALS_OUTPUT(1:N_LOUTPUT)   = LPARTIALS_OUTPUT(1:N_LOUTPUT)
      ELSE
         LRRS_FixIn%UserVal%LBOUNDARIES_OUTPUT(1:N_LOUTPUT) = LBOUNDARIES_OUTPUT(1:N_LOUTPUT)
      ENDIF

!  Normal return

      RETURN

!mick mod 7/20/2016 - modified error handling to promote more consistency with LIDORT

!  Open file error - abort

300   CONTINUE
      MAIL = 'openfile failure for '//FILNAM(1:LEN_STRING(FILNAM))
      CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
      !CLOSE(FILUNIT)
      GO TO 455

!  Line read error - abort immediately

998   CONTINUE
      MAIL = 'read failure for '//PAR_STR(1:LEN_STRING(PAR_STR))
      CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
      CLOSE(FILUNIT)
      GO TO 455

!  Write collected error messages

455   CONTINUE
      STATUS = LRRS_SERIOUS

      !OPEN ( 44, FILE = EFILNAM, STATUS = 'UNKNOWN' )
      !write(44,'(a,t40,i5)')'Total number of error messages = ',NMESS
      !DO n = 1, NMESS
      !  write(44,'(a4,i5,a)')'  # ',N,MESS(N)
      !ENDDO
      !CLOSE(44)

!  Finish

      RETURN
      END SUBROUTINE LRRS_READ_INPUTS

!

      SUBROUTINE LRRS_CHECK_INPUT_DIMS &
        ( LRRS_FixIn, LRRS_ModIn, &
          STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Purpose: Check input dimensions

!  Use module, dimensions and numbers
!   -- Rob mod 5/12/17 for 2p5a, removed MAX_MOMENTS_INPUT

      USE LRRS_pars_m, only : MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_LOUTPUT, &
                              MAX_STREAMS, MAX_LAYERS, MAX_FINE_LAYERS,        &
                              MAX_MESSAGES, LRRS_SUCCESS, LRRS_SERIOUS

      USE LRRS_Inputs_def_m

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  LRRS input structures

      TYPE(LRRS_Fixed_Inputs)   , INTENT (IN) :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT (IN) :: LRRS_ModIn

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = LRRS_SUCCESS
      NM = NMESSAGES

!  Check LRRS input dimensions against maximum dimensions
!  ========================================================

!  1a. Basic dimensions - always checked

      IF ( LRRS_FixIn%Cont%NSTREAMS .GT. MAX_STREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams NSTREAMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_STREAMS dimension'
        STATUS       = LRRS_SERIOUS
      ENDIF

      IF ( LRRS_FixIn%Cont%NLAYERS .GT. MAX_LAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of layers NLAYERS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_LAYERS dimension'
        STATUS       = LRRS_SERIOUS
      ENDIF

      IF ( LRRS_FixIn%UserVal%N_LOUTPUT .GT. MAX_LOUTPUT ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of user vertical output levels N_LOUTPUT > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_LOUTPUT dimension'
        STATUS       = LRRS_SERIOUS
      ENDIF

!   -- Rob mod 5/12/17 for 2p5a, removed this check on NMOMENTS_INPUT
!      IF ( LRRS_FixIn%Cont%NMOMENTS_INPUT .GT. MAX_MOMENTS_INPUT ) THEN
!        NM = NM + 1
!        MESSAGES(NM) = 'Number of Legendre moments NMOMENTS_INPUT > maximum dimension'
!        ACTIONS(NM)  = 'Re-set input value or increase MAX_MOMENTS_INPUT dimension'
!        STATUS       = LRRS_SERIOUS
!      ENDIF

!  1b. Basic dimensions - conditionally checked
!   -- Rob mod 5/12/17 for 2p5a, removed SSFULL from this check.

      IF ( LRRS_FixIn%Bool%DO_SSCORR_OUTGOING ) THEN
        IF ( LRRS_FixIn%Cont%NFINELAYERS .GT. MAX_FINE_LAYERS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of fine layers NFINELAYERS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_FINE_LAYERS dimension'
          STATUS       = LRRS_SERIOUS
        ENDIF
      ENDIF

!  2a. Geometry dimensions - always checked

      IF ( LRRS_FixIn%UserVal%N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of relative azimuths N_USER_RELAZMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = LRRS_SERIOUS
      ENDIF

!  2b. Geometry dimensions - conditionally checked

      IF ( LRRS_ModIn%MBool%DO_USER_STREAMS ) THEN
        IF ( LRRS_FixIn%UserVal%N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of user streams N_USER_STREAMS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS       = LRRS_SERIOUS
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

      END SUBROUTINE LRRS_CHECK_INPUT_DIMS

!

      SUBROUTINE LRRS_CHECK_INPUTS &
       ( LRRS_FixIn, LRRS_ModIn, N_USER_STREAMS, N_USER_RELAZMS, N_LOUTPUT, & ! Inputs
         SOLAR_ANGLE, USER_ANGLES, USER_RELAZMS, NLAYERS, HEIGHT_GRID,      & ! Inputs
         LBOUNDARIES_OUTPUT, LPARTIALS_OUTPUT, EARTH_RADIUS,                & ! Inputs
         STATUS, N_MESSAGES, MESSAGES )                                       ! Output

!  Version 2.5. removal of specific BRDF information

!  Checking of control inputs

!  Used Modules

      USE LRRS_PARS_m
      USE LRRS_INPUTS_DEF_m

      IMPLICIT NONE

!  Input arguments
!  ---------------

      TYPE(LRRS_Fixed_Inputs), INTENT(IN)    :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT(IN) :: LRRS_ModIn

!  Number of layers

      INTEGER  , INTENT(IN) ::   NLAYERS

!  Solar beam input geometry

      REAL(FPK), INTENT(IN) ::   SOLAR_ANGLE

!  User stream variables

      INTEGER  , INTENT(IN) ::   N_USER_STREAMS
      REAL(FPK), INTENT(IN) ::   USER_ANGLES  ( MAX_USER_STREAMS )

!  User azimuth variables

      INTEGER  , INTENT(IN) ::   N_USER_RELAZMS
      REAL(FPK), INTENT(IN) ::   USER_RELAZMS ( MAX_USER_RELAZMS )

!  Height grid

      REAL(FPK), INTENT(IN) ::   HEIGHT_GRID ( 0: MAX_LAYERS )

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)
!    LBOUNDARIES_OUTPUT = layer boundary choices
!    LPARTIALS_OUTPUT = off-boundary level choices

      INTEGER  , INTENT(IN) ::   N_LOUTPUT
      INTEGER  , INTENT(IN) ::   LBOUNDARIES_OUTPUT ( MAX_LOUTPUT )
      REAL(FPK), INTENT(IN) ::   LPARTIALS_OUTPUT   ( MAX_LOUTPUT )

!  Earth Radius

      REAL(FPK), INTENT(IN) ::   EARTH_RADIUS

!  Surface input, all removed for Versioni 2.5, 9/11/15

!      REAL(FPK), INTENT(IN) :: LAMBERTIAN_FRACTION
!      INTEGER  , INTENT(IN) ::      WHICH_BRDF
!      INTEGER  , INTENT(INOUT) ::   BRDF_NPARS
!      REAL(FPK), INTENT(INOUT) :: BRDF_PARS ( MAX_BRDF_PARAMETERS)
!      CHARACTER (LEN=10), INTENT(IN) ::  BRDF_NAMES

!  Output arguments
!  ----------------

!  Error handling

      INTEGER  , INTENT(OUT) ::  STATUS
      INTEGER  , INTENT(OUT) ::  N_MESSAGES
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGES ( MAX_MESSAGES )

!  Local variables
!  ---------------

!  BRDF flags removed, general surface flags added for Version 2.5
!   -- Rob mod 5/12/17 for 2p5a, added SSCORR_NADIR, renamed SSFULL --> SSCORR_ALONE.

      LOGICAL :: &
         DO_RRS_OVERALL,      DO_DIRECTRRS_ONLY,    DO_ELASTIC_ONLY, &
         DO_CABANNES_RAMAN,   DO_ENERGY_BALANCING,  DO_MSMODE_LRRS, &
         DO_BIN_REALIZATION,  DO_MONO_REALIZATION,  DO_DIRECT_BEAM, &
         DO_SSCORR_NADIR,     DO_SSCORR_OUTGOING,   DO_SSCORR_ALONE, DO_MOLECSCAT_ONLY, &
         DO_PLANE_PARALLEL,   DO_DELTAM_SCALING,    DO_DOUBLE_CONVTEST, &
         DO_UPWELLING,        DO_DNWELLING,         DO_USER_STREAMS, &
         DO_NO_AZIMUTH,       DO_LBOUNDARIES,       DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, &
         DO_LRRS_WRITESCENARIO, DO_LRRS_WRITERESULTS, &
         DO_LRRS_WRITEINPUT,    DO_LRRS_WRITEFOURIER

      INTEGER ::   UTA, N, N1, I, NSTART
      REAL(FPK) :: D1
      CHARACTER (LEN=2) :: C2
      LOGICAL ::   LOOP

!  BRDF Kernel names (check). Removed Version 2.5, 9/11/15
 !     CHARACTER (LEN=10) :: BRDF_CHECK_NAMES ( MAXBRDF_IDX )

!  Define Kernel names. Removed Version 2.5, 9/11/15
!      BRDF_CHECK_NAMES(1) = 'Lambertian'
!      BRDF_CHECK_NAMES(2) = 'Hapke     '
!      BRDF_CHECK_NAMES(3) = 'Rahman    '
!      BRDF_CHECK_NAMES(4) = 'Cox-Munk  '

!  Initialize output status
!  ------------------------

      STATUS     = LRRS_SUCCESS
! mick fix 3/31/11 - initialise these
      N_MESSAGES = 0
      MESSAGES   = ' '

!  Expand type structure variables
!  -------------------------------

!  Fixed Boolean Inputs
!   -- Rob mod 5/12/17 for 2p5a, added SSCORR_NADIR, renamed SSFULL --> SSCORR_ALONE.

      DO_RRS_OVERALL         = LRRS_FixIn%Bool%DO_RRS_OVERALL
      DO_DIRECTRRS_ONLY      = LRRS_FixIn%Bool%DO_DIRECTRRS_ONLY
      DO_DIRECT_BEAM         = LRRS_FixIn%Bool%DO_DIRECT_BEAM

      DO_SSCORR_NADIR        = LRRS_FixIn%Bool%DO_SSCORR_NADIR
      DO_SSCORR_OUTGOING     = LRRS_FixIn%Bool%DO_SSCORR_OUTGOING
      DO_SSCORR_ALONE        = LRRS_FixIn%Bool%DO_SSCORR_ALONE

      DO_MOLECSCAT_ONLY      = LRRS_FixIn%Bool%DO_MOLECSCAT_ONLY
      DO_PLANE_PARALLEL      = LRRS_FixIn%Bool%DO_PLANE_PARALLEL
      DO_DELTAM_SCALING      = LRRS_FixIn%Bool%DO_DELTAM_SCALING
      DO_DOUBLE_CONVTEST     = LRRS_FixIn%Bool%DO_DOUBLE_CONVTEST
      DO_UPWELLING           = LRRS_FixIn%Bool%DO_UPWELLING
      DO_DNWELLING           = LRRS_FixIn%Bool%DO_DNWELLING
      DO_LBOUNDARIES         = LRRS_FixIn%Bool%DO_LBOUNDARIES
      DO_MVOUT_ONLY          = LRRS_FixIn%Bool%DO_MVOUT_ONLY
      DO_ADDITIONAL_MVOUT    = LRRS_FixIn%Bool%DO_ADDITIONAL_MVOUT
      DO_LRRS_WRITESCENARIO  = LRRS_FixIn%Bool%DO_LRRS_WRITESCENARIO
      DO_LRRS_WRITERESULTS   = LRRS_FixIn%Bool%DO_LRRS_WRITERESULTS
      DO_LRRS_WRITEINPUT     = LRRS_FixIn%Bool%DO_LRRS_WRITEINPUT
      DO_LRRS_WRITEFOURIER   = LRRS_FixIn%Bool%DO_LRRS_WRITEFOURIER

!  Surface variables removed. Version 2.5 9/11/15
!      DO_LAMBERTIAN_SURFACE  = LRRS_FixIn%Bool%DO_LAMBERTIAN_SURFACE
!      DO_WATERLEAVING        = LRRS_FixIn%Bool%DO_WATERLEAVING
!      DO_SHADOW_EFFECT       = LRRS_FixIn%Bool%DO_SHADOW_EFFECT
!      DO_GLITTER_DBMS        = LRRS_FixIn%Bool%DO_GLITTER_DBMS
!      DO_LAMBERTIAN_FRACTION = LRRS_FixIn%Bool%DO_LAMBERTIAN_FRACTION

!  Modified Boolean Inputs

      DO_ELASTIC_ONLY     = LRRS_ModIn%MBool%DO_ELASTIC_ONLY
      DO_CABANNES_RAMAN   = LRRS_ModIn%MBool%DO_CABANNES_RAMAN
      DO_ENERGY_BALANCING = LRRS_ModIn%MBool%DO_ENERGY_BALANCING
      DO_MSMODE_LRRS      = LRRS_ModIn%MBool%DO_MSMODE_LRRS
      DO_BIN_REALIZATION  = LRRS_ModIn%MBool%DO_BIN_REALIZATION
      DO_MONO_REALIZATION = LRRS_ModIn%MBool%DO_MONO_REALIZATION
      DO_USER_STREAMS     = LRRS_ModIn%MBool%DO_USER_STREAMS
      DO_NO_AZIMUTH       = LRRS_ModIn%MBool%DO_NO_AZIMUTH

!  Check inputs (both file-read and derived)
!  -----------------------------------------

!mick fix 7/20/2016 - added the following input control checks:
!                     * Raman vs Elastic
!                     * Raman Overall vs Raman Direct-Only
!                     * EB vs CR
!                     * Mono vs Bin

!  Check consistency of Raman vs Elastic input controls

      !mick - 1st form
      !IF ( DO_ELASTIC_ONLY ) THEN
      !  IF ( DO_RRS_OVERALL .OR. DO_DIRECTRRS_ONLY ) THEN
      !    N_MESSAGES = N_MESSAGES + 1
      !    MESSAGES(N_MESSAGES) = &
      !        'Bad input: Cannot have both the Elastic flag and Raman control flag(s) set'
      !  ENDIF
      !ELSE
      !  IF ( .NOT.DO_RRS_OVERALL .AND. .NOT.DO_DIRECTRRS_ONLY ) THEN
      !    N_MESSAGES = N_MESSAGES + 1
      !    MESSAGES(N_MESSAGES) = &
      !        'Bad input: At least one of the Elastic or Raman control flag(s) must be set'
      !  ENDIF
      !END IF

      IF ( ( DO_ELASTIC_ONLY .AND. DO_RRS_OVERALL ) .OR. &
           ( .NOT.DO_ELASTIC_ONLY .AND. .NOT.DO_RRS_OVERALL ) ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = &
            'Bad input: Either the Elastic or Raman Overall control flag (one and only one) must be set'
      ENDIF

!  Check consistency of Raman Overall vs Raman Direct-Only input controls

      IF ( .NOT.DO_RRS_OVERALL .AND. DO_DIRECTRRS_ONLY ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = &
            'Bad input: Must have Raman Overall control flag set if Raman Direct-Only control flag is set'
      ENDIF

!  Check consistency of Energy-Balancing vs Cabannes-Raman calculation mode input controls
!     (if performing Raman calculation)

      IF ( DO_RRS_OVERALL ) THEN
        IF ( ( DO_ENERGY_BALANCING .AND. DO_CABANNES_RAMAN ) .OR. &
             ( .NOT.DO_ENERGY_BALANCING .AND. .NOT.DO_CABANNES_RAMAN ) ) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
              'Bad input: Either the EB or CR control flag (one and only one) must be set when performing a Raman calculation'
        ENDIF
      END IF

!  Check consistency of Mono vs Bin realization input controls

      IF ( ( DO_MONO_REALIZATION .AND. DO_BIN_REALIZATION ) .OR. &
           ( .NOT.DO_MONO_REALIZATION .AND. .NOT.DO_BIN_REALIZATION ) ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = &
            'Bad input: Either the Mono or Bin realization control flag (one and only one) must be set'
      END IF

!  Check consistency of mean value input control

      IF ( DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
              'Bad input: Cannot have both mean-value flags set'
        ENDIF
      ENDIF

      IF ( .NOT.DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          IF ( .NOT.DO_NO_AZIMUTH ) THEN
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
              'Bad input: Mean-value option requires NO azimuth'
          ENDIF
          IF ( DO_USER_STREAMS ) THEN
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
             'Bad input: Mean-value option needs quadratures only'
          ENDIF
        ENDIF
      ENDIF

!  Check directional input

      IF ( .NOT.DO_UPWELLING .AND. .NOT. DO_DNWELLING ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = 'Bad input: no directional input is set'
      ENDIF

!  Check Rayleigh-only input against single scatter correction
!    This is disabled now. Wont make any difference with the Nadir SS
!      IF ( DO_RAYLEIGH_ONLY .AND. DO_SSCORRECTION ) THEN
!        N_MESSAGES = N_MESSAGES + 1
!        MESSAGES(N_MESSAGES) =
!     & 'Bad input: SS correction must be turned off, Rayleigh only'
!      ENDIF

!   -- Rob mod 5/12/17 for 2p5a, added SSCORR_NADIR versus  SSCORR_OUTGOING CHECK

      IF ( DO_SSCORR_OUTGOING .and. DO_SSCORR_NADIR) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = &
             'Bad input: Both SS correction flags have been set - not allowed!'
      ENDIF

!  Check single scattering correction and Do no Azimuth
!    ---WARNING. Do-no-Azimuth Flag turned off
!   -- Rob mod 5/12/17 for 2p5a, added SSCORR_NADIR to this check

      IF ( DO_SSCORR_OUTGOING .or. DO_SSCORR_NADIR) THEN
        IF ( DO_NO_AZIMUTH ) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
             'Bad input: need azimuth dependence for SS corrections'
        ENDIF
      ENDIF

!  Check azimuth-only conditions
!  =============================

!  Check no-Azimuth flag

      IF ( .NOT.DO_NO_AZIMUTH ) THEN
        IF ( DO_USER_STREAMS .AND. N_USER_STREAMS.EQ.1 ) THEN
          IF ( USER_ANGLES(1) .EQ. ZERO ) THEN
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
          'Bad input: single zenith-sky output requires no azimuth'
          ENDIF
        ENDIF
      ENDIF

!  Check viewing geometry input
!  ============================

!  Check earth radius (Chapman function only)

      IF ( EARTH_RADIUS.LT.6320.0D0 .OR. &
             EARTH_RADIUS.GT.6420.0D0 ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = &
          'Bad input: Earth radius not in range [6320,6420]'
       ENDIF

!  Check solar zenith angle input

      IF ( SOLAR_ANGLE.LT.ZERO .OR.SOLAR_ANGLE.GE.90.0 ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = &
           ' Bad input: Sun angle not in range [0,90)'
       ENDIF

!  Check relative azimuths

      LOOP = .TRUE.
      I = 0
       DO WHILE (LOOP .AND. I.LT.N_USER_RELAZMS)
        I = I + 1
        IF ( USER_RELAZMS(I) .GT. 360.0   .OR. &
               USER_RELAZMS(I) .LT. ZERO ) THEN
          WRITE(C2,'(I2)')I
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
           'Bad input: out-of-range azimuth angle, no. '//C2
          LOOP = .FALSE.
        ENDIF
      ENDDO

!  Check user-defined stream angles (should always be [0,90])

      IF ( DO_USER_STREAMS ) THEN
        LOOP = .TRUE.
        I = 0
         DO WHILE (LOOP .AND. I.LT.N_USER_STREAMS)
          I = I + 1
          IF ( USER_ANGLES(I) .GT. 90.0   .OR. &
                 USER_ANGLES(I) .LT. ZERO ) THEN
            WRITE(C2,'(I2)')I
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
             'Bad input: out-of-range user stream, no. '//C2
            LOOP = .FALSE.
          ENDIF
        ENDDO
      ENDIF

!  Check height grid input (Chapman function only)

      LOOP = .TRUE.
      I = 0
       DO WHILE (LOOP .AND. I.LT.NLAYERS)
        I = I + 1
        IF ( HEIGHT_GRID(I-1).LE.HEIGHT_GRID(I) ) THEN
          WRITE(C2,'(I2)')I
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
              'Bad input: Height-grid not monotonically decreasing; Layer '//C2
          LOOP = .FALSE.
        ENDIF
      ENDDO

!  Check layer output choices
!  ==========================

!  Check layer boundary choices (should always be within atmosphere!)

      IF ( DO_LBOUNDARIES ) THEN

!  Number of choices should not exceed number of levels in atmosphere
!   0 = TOA, 1 = bottom of first layer, NLAYERS = surface)

        IF ( N_LOUTPUT .GT. NLAYERS + 1 ) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
              'Bad input: Too many Layer boundary outputs'
        ENDIF

!  Choices should remain in atmosphere, i.e. no values > NLAYERS

        DO I = 1, N_LOUTPUT
          IF ( LBOUNDARIES_OUTPUT(I) .GT. NLAYERS ) THEN
            WRITE(C2,'(I2)')LBOUNDARIES_OUTPUT(I)
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
              'Bad input: Layer boundary level > NLAYERS: '//C2
            LOOP = .FALSE.
          ENDIF
        ENDDO
      ENDIF

!  Check repetition of layer boundary output choices

      IF ( DO_LBOUNDARIES ) THEN
        UTA = 0
        LOOP = .TRUE.
        DO WHILE (LOOP .AND. UTA .LT. N_LOUTPUT )
          UTA = UTA + 1
          N1 = LBOUNDARIES_OUTPUT(UTA)
          NSTART = 0
          DO N = 1, N_LOUTPUT
            IF ( LBOUNDARIES_OUTPUT(N) .EQ. N1 ) NSTART = NSTART + 1
          ENDDO
          IF ( NSTART .GT. 1 ) THEN
            LOOP = .FALSE.
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
              'Bad input: repetition of optical depth input value'
          ENDIF
        ENDDO
      ENDIF

!  Check off-boundary choices of output levels
!     Cannot be greater than NLAYERS, cannot be negative

      IF ( .NOT.DO_LBOUNDARIES ) THEN
        DO I = 1, N_LOUTPUT
          IF ( LPARTIALS_OUTPUT(I) .GT. DBLE(NLAYERS) ) THEN
            WRITE(C2,'(I2)')I
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
              'Bad input: off-boundary level > DBLE(NLAYERS), entry #'//C2
            LOOP = .FALSE.
          ELSE IF ( LPARTIALS_OUTPUT(I) .LT. ZERO ) THEN
            WRITE(C2,'(I2)')I
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
              'Bad input: off-boundary level < 0.0, entry #'//C2
            LOOP = .FALSE.
          ENDIF
        ENDDO
      ENDIF

!  Check repetition of user-defined optical depth input

      IF ( .NOT.DO_LBOUNDARIES ) THEN
        UTA = 0
        LOOP = .TRUE.
        DO WHILE (LOOP .AND. UTA .LT. N_LOUTPUT )
          UTA = UTA + 1
          D1 = LPARTIALS_OUTPUT(UTA)
          NSTART = 0
          DO N = 1, N_LOUTPUT
            IF ( LPARTIALS_OUTPUT(N) .EQ. D1 ) NSTART = NSTART + 1
          ENDDO
          IF ( NSTART .GT. 1 ) THEN
            LOOP = .FALSE.
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
            'Bad input: repetition of off-boundary output choice'
          ENDIF
        ENDDO
      ENDIF

!  Check surface inputs
!  ====================

!  Check DB correction only for non-Lambertian surface
!    ---WARNING. Flag turned off. NOw disabled. DB term is automatic
!      IF ( DO_DB_CORRECTION ) THEN
!        IF ( DO_LAMBERTIAN_SURFACE ) THEN
!          N_MESSAGES = N_MESSAGES + 1
!          MESSAGES(N_MESSAGES) =
!     &          'Bad input: No DB correction for Lambertian only'
!          DO_DB_CORRECTION = .FALSE.
!        ENDIF
!      ENDIF

!  Remaining checks from Version 2.3 Removed for Version 2.5, 9/11/15
!  Set status index

      IF ( N_MESSAGES .NE. 0 ) STATUS = LRRS_SERIOUS

!  Finish

      RETURN
      END SUBROUTINE LRRS_CHECK_INPUTS

!

      SUBROUTINE LRRS_DERIVE_INPUTS &
        ( DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                   & ! Inputs
          DO_BIN_REALIZATION, DO_MOLECSCAT_ONLY, DO_LBOUNDARIES,         & ! Inputs
          NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS, N_LOUTPUT,           & ! Inputs
          NLAYERS, NPOINTS_INNER, SOLAR_ANGLE, USER_ANGLES,              & ! Inputs
          NMOMENTS, N_OUT_STREAMS, LBOUNDARIES_OUTPUT, LPARTIALS_OUTPUT, & ! Inputs
          FLUX_FACTOR, COS_SZA, LEVELMASK_UP, LEVELMASK_DN,              & ! Outputs
          STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, N_CONVTESTS,           & ! Outputs
          QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS )         ! Outputs


!  Derivation of bookkeeping inputs

!  include file of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, THREE, DEG_TO_RAD, &
                              MAX_USER_STREAMS, MAX_STREAMS, MAX_LOUTPUT, MAX_LAYERS

!  Ausiliary subroutines

      USE LRRS_AUX2_m, Only : GETQUAD2, RSSORT, ISSORT

      IMPLICIT NONE

!  Input arguments
!  ===============

!  flags

      LOGICAL  , INTENT(IN) :: &
         DO_BIN_REALIZATION,  DO_MOLECSCAT_ONLY,    DO_LBOUNDARIES, &
         DO_UPWELLING,        DO_DNWELLING,         DO_USER_STREAMS

!  Number of layers

      INTEGER  , INTENT(IN) ::   NLAYERS

!  Number of discrete ordinate streams

      INTEGER  , INTENT(IN) ::   NSTREAMS

!  Solar beam input geometry

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLE

!  User stream variables

      INTEGER  , INTENT(IN) ::   N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_ANGLES  ( MAX_USER_STREAMS )

!  Number of azimuths

      INTEGER  , INTENT(IN) ::   N_USER_RELAZMS

!    N_LOUTPUT = number of level output choices (all)

      INTEGER  , INTENT(IN) ::   N_LOUTPUT

!  Number of inner-window points

      INTEGER  , INTENT(IN) ::   NPOINTS_INNER

!  Output arguments
!  ================

!  number of moments in the multiple scatter calculation

      INTEGER  , INTENT(OUT) ::   NMOMENTS

!  number of output streams (quadrature or user)

      INTEGER  , INTENT(OUT) ::   N_OUT_STREAMS

!  number of convergence tests
!    Added by R. Spurr, 18 November 2005

      INTEGER  , INTENT(OUT) ::   N_CONVTESTS

!  Cosein of the solar zenith angle

      REAL(FPK), INTENT(OUT) :: COS_SZA

!  User polar cosines

      REAL(FPK), INTENT(OUT) :: USER_STREAMS ( MAX_USER_STREAMS )

!  Quadrature streams and weights

      REAL(FPK), INTENT(OUT) :: QUAD_STREAMS  ( MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: QUAD_WEIGHTS  ( MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: QUAD_STRMWGT  ( MAX_STREAMS )

!    LBOUNDARIES_OUTPUT = layer boundary choices

      INTEGER  , INTENT(INOUT) :: LBOUNDARIES_OUTPUT ( MAX_LOUTPUT )
      REAL(FPK), INTENT(INOUT) :: LPARTIALS_OUTPUT   ( MAX_LOUTPUT )

!  output optical depth masks and indices

      INTEGER  , INTENT(OUT) :: LEVELMASK_UP ( MAX_LOUTPUT )
      INTEGER  , INTENT(OUT) :: LEVELMASK_DN ( MAX_LOUTPUT )

!  Layer masks for doing integrated source terms

      LOGICAL  , INTENT(OUT) :: STERM_LAYERMASK_UP ( MAX_LAYERS )
      LOGICAL  , INTENT(OUT) :: STERM_LAYERMASK_DN ( MAX_LAYERS )

!  Flux factor

      REAL(FPK), INTENT(INOUT) :: FLUX_FACTOR

!  Bookkeeping variables not passed
!  --------------------------------

!  output optical depth masks and indices

      LOGICAL :: PARTIALS_OUTFLAG  (MAX_LOUTPUT)
      INTEGER :: PARTIALS_OUTINDEX (MAX_LOUTPUT)
      INTEGER :: PARTIALS_LEVELMASK_UP (MAX_LOUTPUT)
      INTEGER :: PARTIALS_LEVELMASK_DN (MAX_LOUTPUT)

!  off-grid optical depths (values, masks, indices)

      INTEGER   :: N_PARTIALS
      INTEGER   :: PARTIALS_LAYERIDX (MAX_LOUTPUT)
      REAL(FPK) :: PARTIALS_VALUES   (MAX_LOUTPUT)

!  number of whole layer source terms required

!      INTEGER :: N_LAYERSOURCE_UP
!      INTEGER :: N_LAYERSOURCE_DN
      INTEGER :: N_ALLLAYERS_UP
      INTEGER :: N_ALLLAYERS_DN

!  local variables
!  ---------------

      INTEGER   :: UTA, N, N1, I, NT, UT, UT1, UM
      REAL(FPK) :: REM, DT
      INTEGER   :: N_DIRECTIONS

!  Get derived inputs (for Bookkeeping)
!  ====================================

!  Flux factor is always one

      FLUX_FACTOR    = ONE

!  solar zenith cosine

      COS_SZA = COS(SOLAR_ANGLE*DEG_TO_RAD)

!  Number of moments

      NMOMENTS = 2 * NSTREAMS - 1
      IF ( DO_MOLECSCAT_ONLY )  NMOMENTS = 2

!  quadrature

      IF ( NSTREAMS .EQ. 1 ) THEN
        QUAD_STREAMS(1) = DSQRT ( ONE / THREE )
        QUAD_WEIGHTS(1) = ONE
      ELSE
        CALL GETQUAD2( ZERO, ONE, NSTREAMS, QUAD_STREAMS, QUAD_WEIGHTS )
      ENDIF

      DO I = 1, NSTREAMS
        QUAD_STRMWGT(I) = QUAD_WEIGHTS(I) * QUAD_STREAMS(I)
      ENDDO

!  number of output streams

      IF ( DO_USER_STREAMS ) THEN
        N_OUT_STREAMS = N_USER_STREAMS
      ELSE
        N_OUT_STREAMS = NSTREAMS
      ENDIF

!  Set up the streams and sines

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          USER_STREAMS(UM) = COS(USER_ANGLES(UM)*DEG_TO_RAD)
        ENDDO
      ENDIF

!  number of tests to be applied for convergence

      N_DIRECTIONS = 0
      IF ( DO_UPWELLING ) N_DIRECTIONS = N_DIRECTIONS + 1
      IF ( DO_DNWELLING ) N_DIRECTIONS = N_DIRECTIONS + 1
      N_CONVTESTS = N_USER_RELAZMS * N_USER_STREAMS * N_DIRECTIONS
      N_CONVTESTS = N_CONVTESTS * N_LOUTPUT

      IF ( DO_BIN_REALIZATION ) THEN
        N_CONVTESTS = N_CONVTESTS * NPOINTS_INNER
      ENDIF

!        write(*,*)nPOINTS_INNER,N_DIRECTIONS,
!     &    N_USER_RELAZMS,N_USER_STREAMS, N_LOUTPUT

!  Optical depth control (layer boundaries only)
!  =============================================

      IF ( DO_LBOUNDARIES ) THEN

!  Sort out User optical depths
!  ----------------------------

!  assign

        N_PARTIALS = 0
        DO UTA = 1, N_LOUTPUT
          PARTIALS_OUTFLAG(UTA) = .FALSE.
          LPARTIALS_OUTPUT(UTA) = -999.0D0
        ENDDO

!  Sort in ascending order

        IF ( N_LOUTPUT .GT. 1 ) THEN
          CALL ISSORT(N_LOUTPUT,LBOUNDARIES_OUTPUT)
!mick fix 4/4/11 - added if around do loop
          IF (LBOUNDARIES_OUTPUT(1) > LBOUNDARIES_OUTPUT(N_LOUTPUT)) THEN
            DO UTA = 1, N_LOUTPUT/2
              UT =  LBOUNDARIES_OUTPUT(UTA)
              UT1 = LBOUNDARIES_OUTPUT(N_LOUTPUT + 1 - UTA)
              LBOUNDARIES_OUTPUT(UTA) = UT1
              LBOUNDARIES_OUTPUT(N_LOUTPUT + 1 - UTA) = UT
            ENDDO
          END IF
        ENDIF

!  set optical depth mask (for values at layer boundaries)

        DO UTA = 1, N_LOUTPUT
          N1 = LBOUNDARIES_OUTPUT(UTA)
          PARTIALS_LEVELMASK_UP(UTA) = N1
          PARTIALS_LEVELMASK_DN(UTA) = N1
        ENDDO

!  set optical depth mask (for values at layer boundaries)

        DO UTA = 1, N_LOUTPUT
          N1 = LBOUNDARIES_OUTPUT(UTA)
          LEVELMASK_UP(UTA) = N1
          LEVELMASK_DN(UTA) = N1
        ENDDO

!  Partial layer distance-related input

      ELSE

!  Sort in ascending order

        IF ( N_LOUTPUT .GT. 1 ) THEN
          CALL RSSORT(N_LOUTPUT,LPARTIALS_OUTPUT)
        ENDIF

!  mark all optical depths not equal to layer boundary values

        UT = 0
        DO UTA = 1, N_LOUTPUT
          DT = LPARTIALS_OUTPUT(UTA)
          NT = INT(DT)
          REM = DT - DBLE(NT)
          IF ( REM .NE. ZERO ) THEN
            UT = UT + 1
            PARTIALS_OUTFLAG(UTA)      = .TRUE.
            PARTIALS_OUTINDEX(UTA)     = UT
            PARTIALS_LAYERIDX(UT)      = NT + 1
            PARTIALS_LEVELMASK_UP(UTA) = NT + 1
            PARTIALS_LEVELMASK_DN(UTA) = NT
            PARTIALS_VALUES(UT)        = REM
          ELSE
            PARTIALS_OUTFLAG(UTA)      = .FALSE.
            PARTIALS_OUTINDEX(UTA)     = 0
            PARTIALS_LEVELMASK_UP(UTA) = NT
            PARTIALS_LEVELMASK_DN(UTA) = NT
          ENDIF
        ENDDO
        N_PARTIALS = UT

      ENDIF

!      write(*,*)n_offgrid_usertaus
!      do ut = 1, n_offgrid_usertaus
!       write(*,*)ut, PARTIALS_LAYERIDX(UT),
!     &        PARTIALS_VALUES(UT)
!      enddo
!      pause

!  Set masking and number of layer source terms
!  --------------------------------------------

!   .. for upwelling and downwelling

!mick fix 4/6/12 - initialize STERM_LAYERMASK_UP outside if block
      DO N = 1, NLAYERS
        STERM_LAYERMASK_UP(N) = .FALSE.
      ENDDO
      IF ( DO_UPWELLING ) THEN
        !DO N = 1, NLAYERS
        !  STERM_LAYERMASK_UP(N) = .FALSE.
        !ENDDO

!mick fix 4/5/11 - substituted two lines with if block
        !N_LAYERSOURCE_UP = LEVELMASK_UP(1) + 1
        !N_ALLLAYERS_UP   = N_LAYERSOURCE_UP
        IF ( .NOT. PARTIALS_OUTFLAG(1) ) THEN
          N_ALLLAYERS_UP = LEVELMASK_UP(1) + 1
        ELSE
          N_ALLLAYERS_UP = PARTIALS_LAYERIDX(1)
        ENDIF

        DO N = NLAYERS, N_ALLLAYERS_UP, -1
          STERM_LAYERMASK_UP(N) = .TRUE.
        ENDDO
      ENDIF

!mick fix 4/6/12 - initialize STERM_LAYERMASK_DN outside if block
      DO N = 1, NLAYERS
        STERM_LAYERMASK_DN(N) = .FALSE.
      ENDDO
      IF ( DO_DNWELLING ) THEN
        !DO N = 1, NLAYERS
        !  STERM_LAYERMASK_DN(N) = .FALSE.
        !ENDDO

!mick fix 4/5/11 - substituted two lines with if block
        !N_LAYERSOURCE_DN = LEVELMASK_DN(N_LOUTPUT)
        !N_ALLLAYERS_DN   = N_LAYERSOURCE_DN
        IF ( .NOT. PARTIALS_OUTFLAG(N_LOUTPUT) ) THEN
          N_ALLLAYERS_DN = LEVELMASK_DN(N_LOUTPUT)
        ELSE
          N_ALLLAYERS_DN = PARTIALS_LAYERIDX(N_PARTIALS)
        ENDIF

        DO N = 1, N_ALLLAYERS_DN
          STERM_LAYERMASK_DN(N) = .TRUE.
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LRRS_DERIVE_INPUTS

!  End module

      END MODULE lrrs_io_check_m

