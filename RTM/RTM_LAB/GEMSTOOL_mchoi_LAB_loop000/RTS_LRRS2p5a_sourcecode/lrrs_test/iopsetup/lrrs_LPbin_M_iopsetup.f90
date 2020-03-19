
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

!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT, added PHASFUNC inputs
!   -- Rob mod 5/12/17 for 2p5a, Included Accessory Module for phase function calculations

      MODULE LPbin_M_iopsetup_m

      USE LRRS_PARS_m

!  Auxiliary and workhorse modules

      USE LRRS_Iopsetup_Aux_m         , Only : Change_Resolution      
      USE LRRS_Iopsetup_WorkRoutines_m, Only : GET_TV8O3_GM1_All, RAYLEIGH_FUNCTION,&
                                               BREON_O3XSEC_All

      PRIVATE
      PUBLIC :: LRRS_IOPSETUP_LPBIN_M

      CONTAINS

      SUBROUTINE LRRS_IOPSETUP_LPBIN_M ( &
        SOLARFILE, DO_JAC, JAC, DO_FD, ND, EPSFAC,    & !Inputs
        LAMBDA_START, LAMBDA_FINISH,                  & !Inputs
        TSHIFT, O3COLUMN,                             & !Inputs
        DO_AEROSOL, AEROSOL_TOPLEVEL,                 & !Inputs
        AEROSOL_SSALB, AEROSOL_ASYMM, AEROSOL_OPDEP,  & !Inputs
        NPOINTS_SURF, ALBEDO, BRDF_Sup_In,            & !Inputs
        LRRS_FixIn, LRRS_ModIn, LRRS_Sup, LRRS_LinIn, & !InOut
        GASVOD, FAIL, N_MESSAGES, MESSAGES )            !Outputs

!  Used modules
!   -- Rob mod 5/12/17 for 2p5a, removed MAX_MOMENTS_INPUT

      USE LRRS_PARS_m, ONLY: MAX_POINTS, MAX_LAYERS !, MAX_MOMENTS_INPUT
      USE LRRS_IO_DEFS_m
      USE LRRS_LIN_IO_DEFS_m

      USE LRRS_RAMAN_SPECTROSCOPY_m, Only : LRRS_BUFFERING_BIN

      USE BRDF_SUP_MOD_m

!   -- Rob mod 5/12/17 for 2p5a, Included Accessory Modules for phase function calculations

      USE LRRS_Phasfunc_Sup_Accessory_m
      !USE LRRS_Phasfunc_LinSup_Accessory_m  (not currently used)

      IMPLICIT NONE

!  Input arguments
!  ---------------

      CHARACTER (LEN=*), INTENT(IN) :: SOLARFILE

!  Jacobian control
!mick mod 7/20/2016 - added DO_JAC & JAC for absolute Jacobian index control

      LOGICAL, INTENT(IN) ::   DO_JAC ( MAX_ATMOSWFS )
      INTEGER, INTENT(IN) ::   JAC

!  FD input

      LOGICAL, INTENT(IN) ::   DO_FD
      INTEGER, INTENT(IN) ::   ND
      REAL(FPK), INTENT(IN) :: EPSFAC

!  Start and finish wavelengths

      REAL(FPK), INTENT(IN) :: LAMBDA_START, LAMBDA_FINISH

!  Choice of Temperature-shift [K], albedo and total ozone [DU]

      REAL(FPK), INTENT(IN) :: TSHIFT, O3COLUMN

!  Aerosol information - same for all wavelengths.

      LOGICAL, INTENT(IN) ::   DO_AEROSOL
      INTEGER, INTENT(IN) ::   AEROSOL_TOPLEVEL
      REAL(FPK), INTENT(IN) :: AEROSOL_SSALB, AEROSOL_ASYMM, AEROSOL_OPDEP

!  Surface

      INTEGER, INTENT(IN)   :: NPOINTS_SURF
      REAL(FPK), INTENT(IN) :: ALBEDO

!  LIDORT BRDF supplement input structure

      TYPE(BRDF_Sup_Inputs), INTENT(IN) :: BRDF_Sup_In

!  Input/Output arguments
!  ----------------------

!  LIDORT-RRS input structures

      TYPE(LRRS_Fixed_Inputs), INTENT(INOUT)    :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT(INOUT) :: LRRS_ModIn
      TYPE(LRRS_Sup_InOut), INTENT(INOUT)       :: LRRS_Sup
      TYPE(LRRS_LinInputs), INTENT(INOUT)       :: LRRS_LinIn

!  Output arguments
!  ----------------

!  Other output

      REAL(FPK), INTENT(OUT) :: GASVOD ( MAX_POINTS )

!  Error Output

      LOGICAL, INTENT(OUT) ::           FAIL
      INTEGER, INTENT(OUT) ::           N_MESSAGES
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGES ( MAX_MESSAGES )

!  Local variables
!  ===============

!  Data variables
!  --------------

!  Local layer inputs

      INTEGER, PARAMETER :: MAX_LOCAL_LAYERS = 200

      INTEGER ::   DATA_NHEIGHTS
      REAL(FPK) :: DATA_HEIGHTS   ( MAX_LOCAL_LAYERS )
      REAL(FPK) :: DATA_PRESSURES ( MAX_LOCAL_LAYERS )
      REAL(FPK) :: DATA_LOGPRESS  ( MAX_LOCAL_LAYERS )
      REAL(FPK) :: DATA_TEMPS     ( MAX_LOCAL_LAYERS )
      REAL(FPK) :: PSURF

!  LRRS input variables, filled here
!  ---------------------------------

!  Control flags

      LOGICAL :: DO_MOLECSCAT_ONLY, DO_LAMBERTIAN_SURFACE, &
                 DO_DELTAM_SCALING, DO_DOUBLE_CONVTEST

!  Outer/inner wavelength range, and offset for the inner range
!    Depends on the choice of solar spectrum

      INTEGER :: OFFSET_INNER
      INTEGER :: NPOINTS_INNER
      INTEGER :: NPOINTS_OUTER

!  For the binning realization
!  ---------------------------

!  Control flag

      LOGICAL :: DO_BIN_REALIZATION

!  Bin upper and lower limits

      REAL(FPK) :: BINLOWER ( MAX_POINTS )
      REAL(FPK) :: BINUPPER ( MAX_POINTS )

!  Wavelengths and Fluxes
!  ----------------------

!  Local inputs for the solar spectrum

      INTEGER, PARAMETER :: MAX_INPUT_POINTS = 1007

      INTEGER ::   N_INPUT_POINTS
      REAL(FPK) :: INPUT_FLUXES  ( MAX_INPUT_POINTS )
      REAL(FPK) :: INPUT_LAMBDAS ( MAX_INPUT_POINTS )

!  These are defined on the outer grid for the binning realization
!  These are defined on 234 points for the monochromatic

!  Depends on the choice of solar spectrum

      REAL(FPK) :: LAMBDAS_RANKED ( MAX_POINTS )
      REAL(FPK) :: FLUXES_RANKED  ( MAX_POINTS )

!  Atmospheric physical quantities
!  -------------------------------

!  Number of layers

      INTEGER :: NLAYERS

!  Height grid

      REAL(FPK) :: HEIGHT_GRID ( 0:MAX_LAYERS )

!  Input layer temperatures, must be in deg K

      REAL(FPK) :: LAYER_TEMPERATURES     ( MAX_LAYERS )
      REAL(FPK) :: TEMPERATURES_UNSHIFTED ( MAX_LAYERS )

!  Input layer Air columns, should be in mol/cm^2 or [DU]

      REAL(FPK) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )
      REAL(FPK) :: LAYER_AIRCOLUMNS_DT ( MAX_LAYERS )

!  Atmospheric optical properties
!  ------------------------------

!  Basic Rayleigh data
!     Rayleigh Cross-sections and depolarization ratios

      REAL(FPK) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  O3 Trace gas Cross-sections. We are using Breon-Daumont data.

      REAL(FPK) :: O3XSECS ( MAX_LAYERS, MAX_POINTS )

!  Aerosol and T-shift Shape function

      LOGICAL ::   AEROSOL_PRESENT ( MAX_LOCAL_LAYERS )
      REAL(FPK) :: AEROSOL_LOADING ( MAX_LOCAL_LAYERS )
      REAL(FPK) :: SHAPEGRID       ( MAX_LOCAL_LAYERS)

!  Aerosol moments in this example
!   -- Rob mod 5/12/17 for 2p5a

      INTEGER, PARAMETER :: MAX_LOCAL_MOMENTS = 350
      INTEGER ::   NMOMC
      REAL(FPK) :: BETAC(0:MAX_LOCAL_MOMENTS)

!  Number of input phase function Legendre moments
!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT
!      INTEGER ::   NMOMENTS_INPUT

!  Local elastic optical properties
!   -- Rob mod 5/12/17 for 2p5a, Changed dimension in PHASMOMS
!mick mod 11/28/2018 - changed MAX_MOMENTS to MAX_LOCAL_MOMENTS in PHASMOMS_TEMPO and
!                      OMEGAMOMS_ELASTIC_UNSCALED for comparisons with LRRS 2p5

      REAL(FPK) :: DELTAU_TEMPO   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: OMEGA_TEMPO    ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: PHASMOMS_TEMPO ( 0:MAX_LOCAL_MOMENTS, MAX_LAYERS,  MAX_POINTS )

!  Unscaled quantities, Elastic input
!    These must be defined on the outer wavelength grid (binning)
!   -- Rob mod 5/12/17 for 2p5a, Changed dimension in OMEGAMOMS_ELASTIC_UNSCALED

      REAL(FPK) :: DELTAU_INPUT_UNSCALED      ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: OMEGAMOMS_ELASTIC_UNSCALED ( MAX_LAYERS, 0:MAX_LOCAL_MOMENTS, MAX_POINTS )

!   -- Rob mod 5/12/17 for 2p5a, added PHASFUNC inputs

      REAL(FPK) :: OMEGAPHASFUNC_ELASTIC_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: OMEGAPHASFUNC_ELASTIC_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Proxy and Help Variables for the PHASFUNC inputs
!  ------------------------------------------------

!   -- Rob mod 5/12/17 for 2p5a

   logical   :: DO_UPWELLING, DO_DNWELLING
   integer   :: N_USER_ANGLES, N_USER_RELAZMS, N_GEOMETRIES
   real(fpk) :: THETA_BOA, USER_ANGLES(MAX_USER_STREAMS), USER_RELAZMS(MAX_USER_RELAZMS)
   logical   :: DO_SCATANG, DO_LEGCALC
   real(fpk) :: COSSCAT_UP (MAX_GEOMETRIES)
   real(fpk) :: COSSCAT_DN (MAX_GEOMETRIES)

!  Legendre calculations

   real(fpk) :: LEGCALC_UP (MAX_GEOMETRIES,0:MAX_LOCAL_MOMENTS)
   real(fpk) :: LEGCALC_DN (MAX_GEOMETRIES,0:MAX_LOCAL_MOMENTS)

!  True output phase functions

   real(fpk) :: PHASFUNC_UP(MAX_GEOMETRIES)
   real(fpk) :: PHASFUNC_DN(MAX_GEOMETRIES)
   real(fpk) :: RAYFUNC_UP (MAX_GEOMETRIES,MAX_POINTS)
   real(fpk) :: RAYFUNC_DN (MAX_GEOMETRIES,MAX_POINTS)
   real(fpk) :: PHASFUNC_ALL(MAX_GEOMETRIES)
   real(fpk) :: L_PHASFUNC_ALL(MAX_GEOMETRIES)

!  Weighting and Rayleigh moment

   REAL(FPK) :: AERWT(MAX_LAYERS,MAX_POINTS), RAYWT(MAX_LAYERS,MAX_POINTS), RAYMOM(MAX_POINTS)
   REAL(FPK) :: L_AERWT(MAX_ATMOSWFS,MAX_LAYERS,MAX_POINTS)
   REAL(FPK) :: L_RAYWT(MAX_ATMOSWFS,MAX_LAYERS,MAX_POINTS)

!  Surface properties
!  ------------------

!  Lambertian albedo

      !REAL(FPK) :: ALBEDOS_RANKED  ( MAX_POINTS )

!  LRRS BRDF supplement

      TYPE(BRDF_Sup_Outputs)                :: BRDF_Sup_Out
      TYPE(BRDF_Output_Exception_Handling)  :: BRDF_Sup_OutputStatus

      LOGICAL   :: DO_DEBUG_RESTORATION
      INTEGER   :: NUSERS, NAZIMS, NDISOS, NMOMS
      REAL(FPK) :: LAMBDA_SURF ( MAX_POINTS )

!  Linearized inputs
!  -----------------

!  Linearization control

      LOGICAL :: DO_PROFILE_LINEARIZATION, &
                 DO_NORMALIZED_WFS !, DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS

      INTEGER :: N_TOTALCOLUMN_WFS
      LOGICAL :: LAYER_VARY_FLAG   (MAX_LAYERS)
      INTEGER :: LAYER_VARY_NUMBER (MAX_LAYERS)

!  Local linearized elastic optical properties
!   -- Rob mod 5/12/17 for 2p5a, Changed MAX_MOMENTS_INPUT dimension
!mick mod 7/20/2016 - changed 1st dim from 4 to MAX_ATMOSWFS
!mick mod 11/28/2018 - changed MAX_MOMENTS to MAX_LOCAL_MOMENTS in L_PHASMOMS_TEMPO and
!                      L_OMEGAMOMS_ELASTIC_UNSCALED for comparisons with LRRS 2p5

      REAL(FPK) :: L_DELTAU_TEMPO   ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L_OMEGA_TEMPO    ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L_PHASMOMS_TEMPO ( MAX_ATMOSWFS, 0:MAX_LOCAL_MOMENTS, MAX_LAYERS, MAX_POINTS )

!  Linearized unscaled optical properties, elastic input
!   -- Rob mod 5/12/17 for 2p5a, Changed MAX_MOMENTS_INPUT dimension

      REAL(FPK) :: L_DELTAU_INPUT_UNSCALED      ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L_OMEGAMOMS_ELASTIC_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_LOCAL_MOMENTS, MAX_POINTS )

!  Linearized unscaled optical properties, elastic input
!   -- Rob mod 5/12/17 for 2p5a, added PHASFUNC inputs

      REAL(FPK) :: L_OMEGAPHASFUNC_ELASTIC_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: L_OMEGAPHASFUNC_ELASTIC_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Input/output for the TOMS V8 O3 profile data set
!  ------------------------------------------------

!  Latitude and Time. User must specify.

      REAL(FPK) :: LATITUDE
      INTEGER   :: YEAR, MONTH, DAY_OF_MONTH

!  Ozone profile amounts in each layer [mol/cm/cm]
!    [From TOMS V8 climatology, converted from [DU]

      REAL(FPK) :: LAYER_O3COLUMNS    ( MAX_LAYERS )

!  Column derivatives. Additional output, not used further.

      REAL(FPK) :: LAYER_O3COLUMNS_dC ( MAX_LAYERS )

!  O3 cross-section Bass-Paur parameterization

      REAL(FPK) :: O3C0 ( MAX_POINTS )
      REAL(FPK) :: O3C1 ( MAX_POINTS )
      REAL(FPK) :: O3C2 ( MAX_POINTS )

!  Other local variables
!  ---------------------

!   -- Rob mod 5/12/17 for 2p5a, Local number of moments

      INTEGER   :: LOCAL_NMOMENTS

!  Help variables

      INTEGER ::   I, J, L, N, W, Q, NJACS, K, V
      REAL(FPK) :: ABS_DP, SCM_DA, SCM_DT, SCT, ABS_DT
      REAL(FPK) :: TAU, SCM, ABS, SCA, AWT, RWT, CRIT_2, RHO
      REAL(FPK) :: RAY2, OMEGA, OMEGA1, AVT, PROD, AER, HELP, HELP2
      REAL(FPK) :: POT_1, POT_2, RHO_A, TEMP, TEMPSQ, DIFF
      REAL(FPK) :: AIRPROFILE, O3_XSECS_DT, DTAU
      REAL(FPK) :: MOLVOD ( MAX_POINTS ), D1, D0, GASCON

!  Change resolution

      LOGICAL ::   DO_CHANGE_RESOLUTION
      INTEGER ::   NRES
      LOGICAL ::   FINER_RES, HALF_RES, QTR_RES

!  Actual values

      REAL(FPK) :: ACTUAL_OZONE, ACTUAL_SHIFT
      REAL(FPK) :: ACTUAL_OPDEP

!  Local error handling

      CHARACTER (LEN=120) :: LOCAL_MESSAGE

!  Parameters
!  ----------

!  Loschmidt's number (particles/cm3).

      REAL(FPK), PARAMETER :: RHO_STANDARD = 2.68675D+19

!  STP values.

      REAL(FPK), PARAMETER :: PZERO = 1013.25D0, TZERO = 273.15D0
      REAL(FPK), PARAMETER :: PTZERO = PZERO / TZERO

      REAL(FPK), PARAMETER :: RHO_ZERO = RHO_STANDARD*TZERO/PZERO
      REAL(FPK), PARAMETER :: CONSTANT = 1.0D+05*RHO_ZERO

!  const in: Col[cm2] =  2.69e16*col[DU]:

      REAL(FPK), PARAMETER :: DU_TO_CM2 = 2.68668D16

!  CO2 PPMV mixing ratio (for the Rayleigh stuff)

      REAL(FPK), PARAMETER :: CO2_PPMV_MIXRATIO = 390.0d0

!  For debug

      LOGICAL :: DEBUG = .FALSE.
!      LOGICAL :: DEBUG = .TRUE.

!  Start code
!  ----------

!  Set local input variables

      DO_BIN_REALIZATION       = LRRS_ModIn%MBool%DO_BIN_REALIZATION
      DO_LAMBERTIAN_SURFACE    = .not.LRRS_FixIn%Bool%DO_BRDF_SURFACE

      DO_PROFILE_LINEARIZATION = LRRS_LinIn%Cont%DO_PROFILE_WFS

      DO_NORMALIZED_WFS        = LRRS_LinIn%Cont%DO_NORMALIZED_WFS
      !DO_AIRPROFILE_WFS        = LRRS_LinIn%Cont%DO_AIRPROFILE_WFS
      !DO_TEMPPROFILE_WFS       = LRRS_LinIn%Cont%DO_TEMPPROFILE_WFS

!  Initialise subroutine status

      FAIL       = .FALSE.
      N_MESSAGES = 0
      MESSAGES   = ' '

!  Initialise other main output

      NPOINTS_INNER     = 0
      OFFSET_INNER      = 0
      NPOINTS_OUTER     = 0

      BINLOWER = ZERO
      BINUPPER = ZERO

      LAMBDAS_RANKED = ZERO
      FLUXES_RANKED  = ZERO

      NLAYERS            = 0
      HEIGHT_GRID        = ZERO
      LAYER_TEMPERATURES = ZERO
      LAYER_AIRCOLUMNS   = ZERO

      RAYLEIGH_XSEC      = ZERO
      RAYLEIGH_DEPOL     = ZERO

!   -- Rob mod 5/12/17 for 2p5a, removed
!      NMOMENTS_INPUT     = 0

      DELTAU_INPUT_UNSCALED      = ZERO
      OMEGAMOMS_ELASTIC_UNSCALED = ZERO

!   -- Rob mod 5/12/17 for 2p5a, Introducing the phase function product input

      OMEGAPHASFUNC_ELASTIC_UP = ZERO
      OMEGAPHASFUNC_ELASTIC_DN = ZERO

!  Linearization initializing
!   -- Rob mod 5/12/17 for 2p5a, Introducing the phase function product input

      N_TOTALCOLUMN_WFS = 0

      LAYER_VARY_FLAG   = .FALSE.
      LAYER_VARY_NUMBER = 0

      LAYER_AIRCOLUMNS_DT    = ZERO
      TEMPERATURES_UNSHIFTED = ZERO

      L_DELTAU_INPUT_UNSCALED      = ZERO
      L_OMEGAMOMS_ELASTIC_UNSCALED = ZERO

      L_OMEGAPHASFUNC_ELASTIC_UP = ZERO
      L_OMEGAPHASFUNC_ELASTIC_DN = ZERO

!  Change resolution settings

      DO_CHANGE_RESOLUTION = .false.
      FINER_RES = .false. ; NRES    = 0
      HALF_RES  = .false. ; QTR_RES = .false.

!  Define local optical property variables

      Actual_ozone  = O3COLUMN
      Actual_opdep  = AEROSOL_OPDEP
      Actual_shift  = TSHIFT

!  Input data of solar wavelengths and fluxes (USER FILE-READ)
!  ------------------------------------------

      OPEN(1,FILE=SOLARFILE,STATUS='OLD')
      READ(1,*)N_INPUT_POINTS
      DO I = 1, N_INPUT_POINTS
        READ(1,*)INPUT_LAMBDAS(I), INPUT_FLUXES(I)
      ENDDO
      CLOSE(1)

!  Next, Perform buffering in the BINNING realization
! --- -----------------------------------------------

!  This is not required for the monochromatic mode

!   given the starting and ending wavelengths:
!      - Set the inner and outer  windows and the inner window offset
!      - Set the wavelengths and fluxes for the outer window
!      - Get the binning for the outer window

!   The third argument should be 0.0d0 always in this call
!    Fourth argument should be TRUE only for pure elastic output.

      IF ( DO_BIN_REALIZATION ) THEN
        CRIT_2 = 0.0D0
        CALL LRRS_BUFFERING_BIN &
           ( MAX_INPUT_POINTS, MAX_POINTS, CRIT_2, .false., &
             N_INPUT_POINTS, INPUT_LAMBDAS, INPUT_FLUXES, &
             LAMBDA_START, LAMBDA_FINISH, &
             NPOINTS_INNER, NPOINTS_OUTER, OFFSET_INNER, &
             LAMBDAS_RANKED, FLUXES_RANKED, &
             BINLOWER, BINUPPER, FAIL, N_MESSAGES, MESSAGES )
        IF (FAIL )THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = '  These messages from LPbin_M_Iopsetup'
          RETURN
        ENDIF
      ENDIF

!  Next, Read a file of heights and temperatures
!  ---------------------------------------------

!  Get the temperatures and heights from FILE.

!     OPEN(1,FILE='lrrs_test/physics_data/ATMOS/temp_psurf.prf_35',STATUS='OLD')
      OPEN(1,FILE='lrrs_test/physics_data/ATMOS/temp_psurf.prf_17',STATUS='OLD')
      READ(1,*)DATA_NHEIGHTS, psurf
      DO I = 1, DATA_NHEIGHTS
          READ(1,*)DATA_HEIGHTS(DATA_NHEIGHTS-I+1), &
                   DATA_TEMPS(DATA_NHEIGHTS-I+1)
      ENDDO
      CLOSE(1)

!  Change the resolution (if flagged)

      IF ( DO_CHANGE_RESOLUTION ) THEN
         call Change_resolution &
            ( MAX_LOCAL_LAYERS, NRES, FINER_RES, HALF_RES, QTR_RES, &
              DATA_NHEIGHTS, DATA_HEIGHTS, DATA_TEMPS )
      ENDIF

!  Set the pressures using the Hydrostatic Equation
!  ------------------------------------------------

!  Gas constant

      GASCON = 9.81*28.9/8314.0d0

!  Start with the ground surface pressure (psurf)

      DATA_PRESSURES(DATA_NHEIGHTS) = psurf
      DATA_LOGPRESS(DATA_NHEIGHTS)  = dlog(psurf)

!  Work upwards from the surface

      DO I = DATA_NHEIGHTS-1, 1 ,-1
        DIFF = 1000.0d0 * ( DATA_HEIGHTS(I) - DATA_HEIGHTS(I+1) )
        AVT  = 0.5d0 * ((1.0d0/DATA_TEMPS(I))+(1.0d0/DATA_TEMPS(I+1)))
        PROD = GASCON * AVT * DIFF
        DATA_PRESSURES(I) = DATA_PRESSURES(I+1)*dexp(-PROD)
        DATA_LOGPRESS(I)  = DLOG(DATA_PRESSURES(I))
      ENDDO

!  debug
!      do n = 1, data_nheights
!        write(*,'(I4,1P3E16.6)') n,data_heights(n),data_temps(n),data_pressures(n)
!      enddo

!  Next, Set NLAYERS, Height Grid, Layer temperatures, Air columns
!  ---------------------------------------------------------------

!     ****** THESE ARE LRRS INPUT VARIABLES

!  Number of layers, original grid

      NLAYERS = DATA_NHEIGHTS - 1

!  Set TOA height

      HEIGHT_GRID(0) = DATA_HEIGHTS(1)

!  Start layer loop

      DO I = 1, NLAYERS

!  Uniform T-shift

        SHAPEGRID(I) = 1.0d0

!  RHO_A is the average layer air density [mol/m/cm/cm]
!  DIFF  is the layer height thickness in [m]

        POT_1 = DATA_PRESSURES(I)   / DATA_TEMPS(I)
        POT_2 = DATA_PRESSURES(I+1) / DATA_TEMPS(I+1)
        TEMP  = 0.5D0*(DATA_TEMPS(I)+DATA_TEMPS(I+1))
        RHO_A = 0.5D0 * CONSTANT * ( POT_1 + POT_2 )
        DIFF  = DATA_HEIGHTS(I) - DATA_HEIGHTS(I+1)
        AIRPROFILE = DIFF * RHO_A * TEMP

!  These next 3 lines are variables required by LRRS

        TEMPERATURES_UNSHIFTED(I) = TEMP
        LAYER_TEMPERATURES(I)   = TEMP + SHAPEGRID(I) * Actual_shift
        LAYER_AIRCOLUMNS(I)     = AIRPROFILE / LAYER_TEMPERATURES(I)
        HEIGHT_GRID(I)          = DATA_HEIGHTS(I+1)

!  End layer loop

      ENDDO

!  Analytic Jacobian control
!   (1 = Total profile, 2 = aerosol opdep profile, 3 = Air profile, 4 = Temp profile)
!   Require T-derivatives of Air. Always un-normalized

      if ( .not.DO_FD .and. JAC.eq.4 ) &
        LAYER_AIRCOLUMNS_dT(1:NLAYERS) = -LAYER_AIRCOLUMNS(1:NLAYERS) / LAYER_TEMPERATURES(1:NLAYERS)

!  Debug
!      do n = 1, NLAYERS
!        write(*,'(I4,1P3E16.6)') n,height_GRID(n),LAYER_tempERATURES(n),layer_aircolumns(n)
!      enddo

!  Next, get the Ozone profile from a basic data set
!  -------------------------------------------------

!   This is the TOMS Version 8 data set.
!     Interpolation to User-defined PTH grid using LOG(pressures)
!     This will deliver the profile and its derivative w.r.t. total column

!    Latitude 45 N, 15 July 2007

      LATITUDE        = 45.0d0
      YEAR            = 2007
      MONTH           = 7
      DAY_OF_MONTH    = 15

!  Call to extract and interpolate data

      CALL GET_TV8O3_GM1_All &
        ( Actual_ozone, LATITUDE, YEAR, MONTH, DAY_OF_MONTH, &
          MAX_LAYERS, MAX_LOCAL_LAYERS, NLAYERS, data_logpress, psurf, &
          LAYER_O3COLUMNS, LAYER_O3COLUMNS_dC, &
          FAIL, LOCAL_MESSAGE )

!  Debug
!      do n = 1, NLAYERS
!        write(*,'(I4,1P3E16.6)') n,height_GRID(n),layer_o3columns(n)
!      enddo
!      write(*,*)SUM(layer_o3columns(1:nlayers))

!  Exception handling

      IF ( FAIL ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = LOCAL_MESSAGE
        RETURN
      END IF

!  FD Jacobian control - perturb ozone
!mick mod 7/20/2016 - switched from using variable kd to jac
!mick fix 7/20/2016 - split up previous fd jac control if block

      if ( do_fd ) then
         do n = 1, nlayers
            if ( nd.eq.n ) then
               if (jac.eq.1) LAYER_O3COLUMNS(n) = epsfac*LAYER_O3COLUMNS(n)
            endif
         enddo
      endif

!  Convert from [DU] to physical units [mol/cm/cm]

      do n = 1, nlayers
        layer_o3columns(n)    = layer_o3columns(n)    * du_to_cm2
        layer_o3columns_dC(n) = layer_o3columns_dC(n) * du_to_cm2
      enddo

!  Next, get the Rayleigh details (from Bodhaine et al. (1999))
!  ------------------------------------------------------------

!  Call returns output for all wavelengths in outer window

      CALL RAYLEIGH_FUNCTION &
           ( MAX_POINTS, CO2_PPMV_MIXRATIO, &
             NPOINTS_OUTER,  LAMBDAS_RANKED, &
             RAYLEIGH_XSEC, RAYLEIGH_DEPOL )

!  Next, get the ozone cross-sections, using the Breon data
!  --------------------------------------------------------

!    This is parameterized according to the quadratic Temperature law

      CALL BREON_O3XSEC_All &
           ( MAX_POINTS, NPOINTS_OUTER, LAMBDAS_RANKED, &
             O3C0, O3C1, O3C2, &
             FAIL, LOCAL_MESSAGE )

!   Debug
!      do n = 1, npoints_outer
!        write(*,*)n,lambdas_ranked(n),o3c0(n),RAYLEIGH_DEPOL(n)
!      enddo

!  Exception handling

      IF ( FAIL ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = LOCAL_MESSAGE
        RETURN
      END IF

!  Next, aerosol properties
!  ------------------------

!  Initialize

      do i = 1, nlayers
         aerosol_loading(i) = 0.0d0
         aerosol_present(i) = .false.
      enddo

!  Loading: Aerosols between ground level and AEROSOL_TOPLEVEL

      if ( do_aerosol ) then
        dtau = Actual_opdep / &
           ( height_grid(aerosol_toplevel) - height_grid(nlayers) )
        diff = 0.0d0
        do i = 1, nlayers
          if ( i.gt.aerosol_toplevel ) then
            aerosol_present(i) = .true.
            aerosol_loading(i) = dtau*(height_grid(i-1)-height_grid(i))

!  FD Jacobian control - perturb aerosol loading
!mick mod 7/20/2016 - switched from using variable kd to jac
!mick fix 7/20/2016 - split up previous fd jac control if block
            if ( do_fd ) then
               if ( nd.eq.i ) then
                  if (jac.eq.2) aerosol_loading(i) = epsfac*aerosol_loading(i)
               endif
            endif

          endif
          diff = diff + aerosol_loading(i)
        enddo
      endif

!  set aerosol coefficients

      if ( do_aerosol ) then
         nmomc = 150
         betac(0) = 1.0d0
         d0 = 1.0d0
         do k = 1, nmomc
            d1 = dble(2*k+1)
            betac(k) = (d1/d0) * aerosol_asymm * betac(k-1)
            d0 = d1
         enddo
      endif

!  Define number of moments to compute
!mick note 11/28/2018 - note that here inside "iopsetup" we are computing
!                       NMOMC moments and not just 2*NSTREAMS moments
!                       in order to crosscheck with LRRS 2p5 which uses
!                       0:NMOMC moments in its SS solution; however,
!                       for actually passing to LRRS 2p5a, we only pass
!                       2*NSTREAMS moments since it only requires that
!                       for its MS solution

      LOCAL_NMOMENTS = 2 ! Rayleigh_only value
      if ( do_aerosol ) LOCAL_NMOMENTS = NMOMC

!  Set some LRRS flags related to aerosol

      if ( do_aerosol ) then
         DO_MOLECSCAT_ONLY  = .false.
         DO_DELTAM_SCALING  = .true.
         DO_DOUBLE_CONVTEST = .true.
      else
         DO_MOLECSCAT_ONLY  = .true.
         DO_DELTAM_SCALING  = .false.
         DO_DOUBLE_CONVTEST = .false.
      endif

!  FD Jacobian control - perturb the required optical property
!mick mod 7/20/2016 - switched from using variable kd to jac
!mick fix 7/20/2016 - moved perturbations of LAYER_O3COLUMNS & AEROSOL_LOADING

      if ( do_fd ) then
         do n = 1, nlayers
            if ( nd.eq.n ) then
               !if (jac.eq.1) LAYER_O3COLUMNS(n)    = epsfac*LAYER_O3COLUMNS(n)
               !if (jac.eq.2) AEROSOL_LOADING(n)    = epsfac*AEROSOL_LOADING(n)
               if (jac.eq.3) LAYER_AIRCOLUMNS(n)   = epsfac*LAYER_AIRCOLUMNS(n)
               if (jac.eq.4) LAYER_TEMPERATURES(n) = epsfac*LAYER_TEMPERATURES(n)
            endif
         enddo
      endif

!  Analytic Jacobian control
!mick mod 7/20/2016 - define njacs & use it to define
!                     layer_vary_number & layer_vary_flag

      if ( .not. do_fd ) then
         !do N = 1, nlayers
         !   layer_vary_number(n) = 4
         !   layer_vary_flag(n)   = .true.
         !enddo
         njacs = COUNT(DO_JAC)
         if (njacs > 0) then
            layer_vary_flag(1:nlayers)   = .true.
            layer_vary_number(1:nlayers) = njacs
         endif
      endif

!  Make ELASTIC optical properties
!  ===============================

!  Initialize

      DELTAU_TEMPO     = 0.0D0
      OMEGA_TEMPO      = 0.0D0
      PHASMOMS_TEMPO   = 0.0D0
      L_DELTAU_TEMPO   = 0.0D0
      L_OMEGA_TEMPO    = 0.0D0
      L_PHASMOMS_TEMPO = 0.0D0

!  Start outer wavelength loop
!  ===========================

!mick fix 3/31/11 - initialise all
      gasvod = zero

!   -- Rob mod 5/12/17 for 2p5a, Need to initialize RAYWT, AERWT and derivatives

      RAYWT = zero ; L_RAYWT = zero
      AERWT = zero ; L_AERWT = zero

!  Start wavelength loop

      do w = 1, npoints_outer

!  depolarization ratio and  Rayleigh second moment

        rho = RAYLEIGH_DEPOL(W)
        ray2 = (1.0d0 - rho ) / ( 2.0d0 + rho ) ; raymom(W) = ray2

!  debug quantities

        gasvod(w) = 0.0d0
        molvod(w) = 0.0d0

!  Start layer loop

        do n = 1, nlayers

!  Bass-Paur quadratic Temperature dependence O3 cross-sections

          temp = layer_temperatures(n) - 273.15d0
          tempsq = temp * temp
          o3xsecs(n,w) = O3c0(w)+o3c1(w)*temp+o3c2(w)*tempsq

!  Normalized derivative w.r.t Temperature

          o3_xsecs_dT = o3c1(w) + 2.0d0*o3c2(w)*temp
          abs_dT = LAYER_O3COLUMNS(n)     * o3_xsecs_dT
          if ( do_normalized_wfs ) abs_dT = abs_dT*LAYER_TEMPERATURES(N)

!  Molecular absorption and scattering

          abs    = LAYER_O3COLUMNS(n)     * O3xsecs(n,w)
          abs_dP = O3xsecs(n,w)
          if ( do_normalized_wfs ) abs_dP = abs

          scm    = LAYER_AIRCOLUMNS(n)    * RAYLEIGH_XSEC(W)
          scm_dA = RAYLEIGH_XSEC(W)
          if ( do_normalized_wfs ) scm_dA = scm

          scm    = LAYER_AIRCOLUMNS(n)    * RAYLEIGH_XSEC(W)
          scm_dT = LAYER_AIRCOLUMNS_dT(n) * RAYLEIGH_XSEC(W)
          if ( do_normalized_wfs ) scm_dT = scm_dT*LAYER_TEMPERATURES(N)

          molvod(w) = molvod(w) + scm
          gasvod(w) = gasvod(w) + abs

!  Aerosols not present in a given layer
!  -------------------------------------

          if ( .not. aerosol_present(n) ) then

!  tau and omega

             tau = abs + scm
             omega  = scm / tau
             if ( omega .gt. 0.999999d0 ) omega = 0.999999d0
             omega1 = 1.0d0 - omega

!  Basic optical properties

             DELTAU_TEMPO(N,W) = tau
             OMEGA_TEMPO(N,W)  = omega
             PHASMOMS_TEMPO(0,N,W) = 1.0D0
             PHASMOMS_TEMPO(1,N,W) = 0.0D0
             PHASMOMS_TEMPO(2,N,W) = RAY2

!   -- Rob mod 5/12/17 for 2p5a, No need to calculate RAYWT (already initialized)

!mick mod 7/20/2016 - added ".not. do_fd" if condition
!                   - switched Q's role from absolute to relative index for flexibility
!                     in doing different Jacobians during a given test
!                   - added "DO_JAC" if conditions

             if ( .not. do_fd ) then
                Q = 0

!  Linearization w.r.t O3 profile

                !Q = 1
                if ( DO_JAC(1) ) then
                   Q = Q + 1
                   L_DELTAU_TEMPO(Q,N,W) =  abs_dP
                   L_OMEGA_TEMPO (Q,N,W) = -abs_dP * omega / tau
                endif

!  Linearization w.r.t aerosol profile; just increment Q here

                !Q = 2
                if ( DO_JAC(2) ) then
                   Q = Q + 1
                endif

!  Linearization w.r.t Air profile

                !Q = 3
                if ( DO_JAC(3) ) then
                   Q = Q + 1
                   L_DELTAU_TEMPO(Q,N,W) = scm_dA
                   L_OMEGA_TEMPO (Q,N,W) = scm_dA * omega1 / tau
                endif

!  Linearization w.r.t Temp profile

                !Q = 4
                if ( DO_JAC(4) ) then
                   Q = Q + 1
                   L_DELTAU_TEMPO(Q,N,W) = scm_dT + abs_dT
                   L_OMEGA_TEMPO (Q,N,W) = (scm_dT*omega1-abs_dT*omega) / tau
                endif
             endif

          endif

!  Aerosol present in the layer
!  ----------------------------

          if ( aerosol_present(n) ) then

             aer = aerosol_loading(n)
             tau = aer + abs + scm
             sca =  aerosol_ssalb * aer
             sct = sca + scm
             rwt = scm / sct
             awt = sca / sct

!   -- Rob mod 5/12/17 for 2p5a, Need to save rwt and awt

             RAYWT(N,W) = rwt  
             AERWT(N,W) = awt  

!  Toggle

             omega = ( sca + scm ) / tau
             if ( omega .gt. 0.999999d0 ) omega = 0.999999d0
             omega1 = 1.0d0 - omega

!  Basic optical properties

             DELTAU_TEMPO(N,W) = tau
             OMEGA_TEMPO(N,W)  = omega

             PHASMOMS_TEMPO(0,N,W) = 1.0d0
             PHASMOMS_TEMPO(1,N,W) = betac(1) * awt
             PHASMOMS_TEMPO(2,N,W) = betac(2) * awt + ray2 * rwt
             do L = 3, LOCAL_NMOMENTS
               PHASMOMS_TEMPO(L,N,W) = betac(L) * awt
             enddo

!mick mod 7/20/2016 - added ".not. do_fd" if condition
!                   - switched Q's role from absolute to relative index for flexibility
!                     in doing different Jacobians during a given test
!                   - added "DO_JAC" if conditions

             if ( .not. do_fd ) then
                Q = 0

!  Linearization w.r.t O3 profile

                !Q = 1
                if ( DO_JAC(1) ) then
                   Q = Q + 1
                   L_DELTAU_TEMPO(Q,N,W) =  abs_dP
                   L_OMEGA_TEMPO (Q,N,W) = -abs_dP * omega / tau
                endif

!  Linearization w.r.t aerosol profile
!   -- Rob mod 5/12/17 for 2p5a, Need to save rwt and awt derivatives
!mick fix 11/28/2018 - replaced "L_AERWT / L_RAYWT / L_PHASMOMS_TEMPO" portion

                !Q = 2
                if ( DO_JAC(2) ) then
                   Q = Q + 1
                   Help = aer / tau
                   L_DELTAU_TEMPO(Q,N,W) = aer
                   L_OMEGA_TEMPO (Q,N,W) =  ( awt - Help ) * omega

                   !Help = rwt * awt
                   !L_RAYWT(Q,N,W) = - help
                   !L_AERWT(Q,N,W) = help     ! Bug 14 June 2017, you had awt * awt
                   !L_PHASMOMS_TEMPO(Q,1,N,W) = Help * betac(1)
                   !L_PHASMOMS_TEMPO(Q,2,N,W) = Help * ( betac(2) - ray2 )
                   !do L = 3, LOCAL_nmoments
                   !   L_PHASMOMS_TEMPO(Q,L,N,W) = Help * betac(L)
                   !enddo

                   Help  = awt * ( one - awt)
                   Help2 = - awt * rwt
                   L_AERWT(Q,N,W) = Help
                   L_RAYWT(Q,N,W) = Help2
                   L_PHASMOMS_TEMPO(Q,1,N,W) = Help * betac(1)
                   L_PHASMOMS_TEMPO(Q,2,N,W) = Help * betac(2) + Help2 * ray2
                   do L = 3, LOCAL_nmoments
                      L_PHASMOMS_TEMPO(Q,L,N,W) = Help * betac(L)
                   enddo

                endif

!  Linearization w.r.t Air Profile
!   -- Rob mod 5/12/17 for 2p5a, Need to save rwt and awt derivatives

                !Q = 3
                if ( DO_JAC(3) ) then
                   Q = Q + 1
                   L_DELTAU_TEMPO(Q,N,W) = scm_dA
                   L_OMEGA_TEMPO (Q,N,W) = scm_dA*omega1 / tau
                   Help = - scm_dA * awt / sct
                   L_RAYWT(Q,N,W) = - help
                   L_AERWT(Q,N,W) = help
                   L_PHASMOMS_TEMPO(Q,1,N,W) = Help * betac(1)
                   L_PHASMOMS_TEMPO(Q,2,N,W) = Help * ( betac(2) - ray2 )
                   do L = 3, LOCAL_nmoments
                      L_PHASMOMS_TEMPO(Q,L,N,W) = Help * betac(L)
                   enddo
                endif

!  Linearization w.r.t Temp Profile
!   -- Rob mod 5/12/17 for 2p5a, Need to save rwt and awt derivatives

                !Q = 4
                if ( DO_JAC(4) ) then
                   Q = Q + 1
                   L_DELTAU_TEMPO(Q,N,W) = scm_dT + abs_dT
                   L_OMEGA_TEMPO (Q,N,W) = (scm_dT*omega1-abs_dT*omega) / tau
                   Help = - scm_dT * awt / sct
                   L_RAYWT(Q,N,W) = - help
                   L_AERWT(Q,N,W) = help
                   L_PHASMOMS_TEMPO(Q,1,N,W) = Help * betac(1)
                   L_PHASMOMS_TEMPO(Q,2,N,W) = Help * ( betac(2) - ray2 )
                   do L = 3, LOCAL_nmoments
                      L_PHASMOMS_TEMPO(Q,L,N,W) = Help * betac(L)
                   enddo
                endif

!  End linearization with aerosol

             endif

!  End aerosol clause

          endif

!  End layer loop

        enddo

!  End wavelength loop

      enddo

!  Develop LRRS input for Elastic IOPs
!  ===================================

!  Initialize

      DELTAU_INPUT_UNSCALED        = 0.0d0
      OMEGAMOMS_ELASTIC_UNSCALED   = 0.0d0
      L_DELTAU_INPUT_UNSCALED      = 0.0d0
      L_OMEGAMOMS_ELASTIC_UNSCALED = 0.0d0

!  setup the elastic scattering
!   -- Rob mod 5/12/17 for 2p5a, NMOMENTS_INPUT is now disabled, replaced by LOCAL_NMOMENTS

      DO W = 1, NPOINTS_OUTER
        DO N = 1, NLAYERS
          DELTAU_INPUT_UNSCALED(N,W)=DELTAU_TEMPO(N,W)
          DO L = 0, LOCAL_NMOMENTS
            OMEGAMOMS_ELASTIC_UNSCALED (N,L,W) = OMEGA_TEMPO(N,W)*PHASMOMS_TEMPO(L,N,W)
          ENDDO
        ENDDO
      ENDDO

      if ( DO_PROFILE_LINEARIZATION ) then
        do w = 1, npoints_outer
          do n = 1, nlayers
            do q = 1, layer_vary_number(n)
              L_DELTAU_INPUT_UNSCALED(q,n,w) = L_DELTAU_TEMPO(q,n,w)
              DO L = 0, LOCAL_NMOMENTS
                L_OMEGAMOMS_ELASTIC_UNSCALED (Q,N,L,W) = &
                  L_OMEGA_TEMPO(Q,N,W) *   PHASMOMS_TEMPO(L,N,W)   + &
                    OMEGA_TEMPO(N,W)   * L_PHASMOMS_TEMPO(Q,L,N,W)
              ENDDO
            enddo
          enddo
        enddo
      endif

!  Debug check #1 (std)

      IF ( DEBUG ) THEN
        DO W = 1, NPOINTS_OUTER
          write(444,'(3i4,1p3e20.10)')W,NLAYERS,LOCAL_NMOMENTS,&
            Lambdas_ranked(w),fluxes_ranked(w)!,albedos_ranked(w)
          DO N = 1, NLAYERS
            write(444,'(2i4,1pe20.10)')w,n,DELTAU_INPUT_UNSCALED(N,W)
            DO L = 0, LOCAL_NMOMENTS
              write(444,'(i4,1pe20.10)') L,OMEGAMOMS_ELASTIC_UNSCALED(N,L,W)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Debug check #1 (lin)

      IF ( DEBUG .AND. DO_PROFILE_LINEARIZATION ) THEN
        DO W = 1, NPOINTS_OUTER
          write(445,'(3i4,1p3e20.10)')W,NLAYERS,LOCAL_NMOMENTS,&
            Lambdas_ranked(w),fluxes_ranked(w)!,albedos_ranked(w)
          DO N = 1, NLAYERS
            DO Q = 1, LAYER_VARY_NUMBER(N)
              write(445,'(3i4,1pe20.10)')w,n,q,L_DELTAU_INPUT_UNSCALED(Q,N,W)
              DO L = 0, LOCAL_NMOMENTS
                write(445,'(i4,1pe20.10)') L,L_OMEGAMOMS_ELASTIC_UNSCALED (Q,N,L,W)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!stop 'at debug check stop 1'

!  -- Rob mod 5/12/17 for 2p5a, New routine for generating PHASFUNC inputs
!  =======================================================================

!  Set the routine up

      DO_SCATANG     = .false.
      DO_LEGCALC     = .false.
      DO_UPWELLING   = LRRS_FixIn%Bool%DO_UPWELLING
      DO_DNWELLING   = LRRS_FixIn%Bool%DO_DNWELLING
      N_USER_ANGLES  = LRRS_FixIn%UserVal%N_USER_STREAMS
      N_USER_RELAZMS = LRRS_FixIn%UserVal%N_USER_RELAZMS
      THETA_BOA      = LRRS_FixIn%Cont%SOLAR_ANGLE
      USER_ANGLES    = LRRS_FixIn%UserVal%USER_ANGLES
      USER_RELAZMS   = LRRS_FixIn%UserVal%USER_RELAZMS
      N_GEOMETRIES   = N_USER_ANGLES * N_USER_RELAZMS

!  Rayleigh first, all points

      CALL Rayfunc_Sup_Accessory &
        ( DO_UPWELLING, DO_DNWELLING,                   & ! Input
          N_USER_ANGLES, N_USER_RELAZMS, NPOINTS_OUTER, & ! Input
          THETA_BOA, USER_ANGLES, USER_RELAZMS, RAYMOM, & ! Input
          DO_SCATANG, COSSCAT_UP, COSSCAT_DN,           & ! Modified input/output
          RAYFUNC_UP, RAYFUNC_DN )                        ! Modified I/O and True output

!  Single set of aerosol coefficients

      if ( do_aerosol ) then
         CALL Phasfunc_Sup_Accessory &
            ( MAX_LOCAL_MOMENTS, DO_UPWELLING, DO_DNWELLING,    & ! Input
             NMOMC, N_USER_ANGLES, N_USER_RELAZMS,              & ! Input
             THETA_BOA, USER_ANGLES, USER_RELAZMS, BETAC,       & ! Input
             DO_LEGCALC, DO_SCATANG, COSSCAT_UP, COSSCAT_DN,    & ! Modified input/output
             LEGCALC_UP, LEGCALC_DN, PHASFUNC_UP, PHASFUNC_DN )   ! Modified I/O and True output
      endif

!  Loop over layers 

      DO N = 1, NLAYERS

!  Mask
!mick fix 11/28/2018 - added "DO_PROFILE_LINEARIZATION" IF blocks

         IF ( aerosol_present(n) ) THEN
            IF ( DO_UPWELLING ) THEN
               DO W = 1, NPOINTS_OUTER
                 DO V = 1, N_GEOMETRIES
                   PHASFUNC_ALL(V) = RAYWT(N,W) * RAYFUNC_UP(V,W) + AERWT(N,W) * PHASFUNC_UP(V)
                   OMEGAPHASFUNC_ELASTIC_UP(N,V,W) = OMEGA_TEMPO(N,W) * PHASFUNC_ALL(V)
                   if ( DO_PROFILE_LINEARIZATION ) then
                     do q = 1, layer_vary_number(n)
                       L_PHASFUNC_ALL(V) = L_RAYWT(Q,N,W) * RAYFUNC_UP(V,W) + L_AERWT(Q,N,W) * PHASFUNC_UP(V)
                       L_OMEGAPHASFUNC_ELASTIC_UP(Q,N,V,W) = &
                              L_OMEGA_TEMPO(Q,N,W) * PHASFUNC_ALL(V) + OMEGA_TEMPO(N,W) * L_PHASFUNC_ALL(V)
                     enddo
                   endif
                 ENDDO
               ENDDO
            ENDIF
            IF ( DO_DNWELLING ) THEN
               DO W = 1, NPOINTS_OUTER
                 DO V = 1, N_GEOMETRIES
                   PHASFUNC_ALL(V) = RAYWT(N,W) * RAYFUNC_DN(V,W) + AERWT(N,W) * PHASFUNC_DN(V)
                   OMEGAPHASFUNC_ELASTIC_DN(N,V,W) = OMEGA_TEMPO(N,W) * PHASFUNC_ALL(V)
                   if ( DO_PROFILE_LINEARIZATION ) then
                     do q = 1, layer_vary_number(n)
                       L_PHASFUNC_ALL(V) = L_RAYWT(Q,N,W) * RAYFUNC_DN(V,W) + L_AERWT(Q,N,W) * PHASFUNC_DN(V)
                       L_OMEGAPHASFUNC_ELASTIC_DN(Q,N,V,W) = &
                              L_OMEGA_TEMPO(Q,N,W) * PHASFUNC_ALL(V) + OMEGA_TEMPO(N,W) * L_PHASFUNC_ALL(V)
                     enddo
                   endif
                 ENDDO
               ENDDO
            ENDIF
         ELSE
            IF ( DO_UPWELLING ) THEN
               DO W = 1, NPOINTS_OUTER
                 DO V = 1, N_GEOMETRIES
                   OMEGAPHASFUNC_ELASTIC_UP(N,V,W) = OMEGA_TEMPO(N,W) * RAYFUNC_UP(V,W)
                   if ( DO_PROFILE_LINEARIZATION ) then
                     do q = 1, layer_vary_number(n)
                        L_OMEGAPHASFUNC_ELASTIC_UP(Q,N,V,W) = L_OMEGA_TEMPO(Q,N,W) * RAYFUNC_UP(V,W)
                     enddo
                   endif
                 ENDDO
               ENDDO
            ENDIF
            IF ( DO_DNWELLING ) THEN
               DO W = 1, NPOINTS_OUTER
                 DO V = 1, N_GEOMETRIES
                   OMEGAPHASFUNC_ELASTIC_DN(N,V,W) = OMEGA_TEMPO(N,W) * RAYFUNC_DN(V,W)
                   if ( DO_PROFILE_LINEARIZATION ) then
                     do q = 1, layer_vary_number(n)
                        L_OMEGAPHASFUNC_ELASTIC_DN(Q,N,V,W) = L_OMEGA_TEMPO(Q,N,W) * RAYFUNC_DN(V,W)
                     enddo
                   endif
                ENDDO
              ENDDO
            ENDIF
         ENDIF

      ENDDO

!  Debug check #2 (std)

      IF ( DEBUG ) THEN
        DO W = 1, NPOINTS_OUTER
          DO N = 1, NLAYERS
            DO V = 1, N_GEOMETRIES
              write(446,'(3i4,1p2e20.10)')W,N,V,&
                OMEGAPHASFUNC_ELASTIC_UP(N,V,W),OMEGAPHASFUNC_ELASTIC_DN(N,V,W)
            ENDDO
          ENDDO
        ENDDO
      END IF

!  Debug check #2 (lin)

      IF ( DEBUG .AND. DO_PROFILE_LINEARIZATION ) THEN
        DO W = 1, NPOINTS_OUTER
          DO N = 1, NLAYERS
            DO V = 1, N_GEOMETRIES
              DO Q = 1, LAYER_VARY_NUMBER(N)
                write(447,'(4i4,1p2e20.10)')W,N,V,Q,&
                  L_OMEGAPHASFUNC_ELASTIC_UP(Q,N,V,W),L_OMEGAPHASFUNC_ELASTIC_DN(Q,N,V,W)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!stop 'at debug check stop 2'

!  Finally, get the surface properties
!    Might want to later call the Water-leaving module.......

      IF ( (NPOINTS_SURF .NE. 1) .AND. (NPOINTS_SURF .NE. NPOINTS_OUTER) ) THEN
        FAIL = .TRUE.
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = 'Error: NPOINTS_SURF should be = 1 or NPOINTS_OUTER'
        RETURN
      ENDIF

      IF ( .NOT.DO_LAMBERTIAN_SURFACE ) THEN
        !Define ranked BRDF(s)
        NUSERS = BRDF_SUP_IN%BS_N_USER_STREAMS
        NAZIMS = BRDF_SUP_IN%BS_N_USER_RELAZMS
        NDISOS = BRDF_SUP_IN%BS_NSTREAMS
        NMOMS  = 2*NDISOS - 1

        DO_DEBUG_RESTORATION = .FALSE.
        IF (NPOINTS_SURF .EQ. 1) THEN
          !Define values at elastic wvls using value @ avg Raman wvl
          LAMBDA_SURF(1) = (LAMBDA_FINISH - LAMBDA_START)/2.0D0
          CALL BRDF_MAINMASTER ( &
            DO_DEBUG_RESTORATION, NMOMS,  & !Inputs
            NPOINTS_SURF, LAMBDA_SURF(1), & !Input wavelengths [nm] for LRRS
            BRDF_Sup_In,                  & !Inputs
            BRDF_Sup_Out,                 & !Outputs
            BRDF_Sup_OutputStatus )         !Output exception handling
        ELSE IF (NPOINTS_SURF .EQ. NPOINTS_OUTER) THEN
          !Define values at elastic wvls using values @ their own wvls (LAMBDAS_RANKED)
          CALL BRDF_MAINMASTER ( &
            DO_DEBUG_RESTORATION, NMOMS,  & !Inputs
            NPOINTS_SURF, LAMBDAS_RANKED, & !Input wavelengths [nm] for LRRS
            BRDF_Sup_In,                  & !Inputs
            BRDF_Sup_Out,                 & !Outputs
            BRDF_Sup_OutputStatus )         !Output exception handling
        ENDIF
      ENDIF

!  Define these LIDORT-RRS inputs
!  ------------------------------

!  Fixed Boolean Inputs

      LRRS_FixIn%Bool%DO_MOLECSCAT_ONLY  = DO_MOLECSCAT_ONLY
      LRRS_FixIn%Bool%DO_DELTAM_SCALING  = DO_DELTAM_SCALING
      LRRS_FixIn%Bool%DO_DOUBLE_CONVTEST = DO_DOUBLE_CONVTEST

!  Fixed Control Inputs
!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT
!      LRRS_FixIn%Cont%NMOMENTS_INPUT  = NMOMENTS_INPUT

      LRRS_FixIn%Cont%NLAYERS        = NLAYERS

!  Fixed Spectral Inputs

      LRRS_FixIn%Spect%LAMBDAS_RANKED = LAMBDAS_RANKED
      LRRS_FixIn%Spect%FLUXES_RANKED  = FLUXES_RANKED

!  For binning calculations

      LRRS_FixIn%Spect%NPOINTS_INNER = NPOINTS_INNER
      LRRS_FixIn%Spect%OFFSET_INNER  = OFFSET_INNER
      LRRS_FixIn%Spect%NPOINTS_OUTER = NPOINTS_OUTER
      LRRS_FixIn%Spect%BINLOWER      = BINLOWER
      LRRS_FixIn%Spect%BINUPPER      = BINUPPER

!  Fixed Atmosphere Inputs
!mick mod 11/28/2018 - limited dims on OMEGAMOMS_ELASTIC_UNSCALED

      LRRS_FixIn%Atmos%HEIGHT_GRID        = HEIGHT_GRID
      LRRS_FixIn%Atmos%LAYER_TEMPERATURES = LAYER_TEMPERATURES
      LRRS_FixIn%Atmos%LAYER_AIRCOLUMNS   = LAYER_AIRCOLUMNS

      LRRS_FixIn%Atmos%RAYLEIGH_XSEC  = RAYLEIGH_XSEC
      LRRS_FixIn%Atmos%RAYLEIGH_DEPOL = RAYLEIGH_DEPOL

      LRRS_FixIn%Atmos%DELTAU_INPUT_UNSCALED      = DELTAU_INPUT_UNSCALED
      !LRRS_FixIn%Atmos%OMEGAMOMS_ELASTIC_UNSCALED = OMEGAMOMS_ELASTIC_UNSCALED
      LRRS_FixIn%Atmos%OMEGAMOMS_ELASTIC_UNSCALED(1:NLAYERS,0:MAX_MOMENTS,1:NPOINTS_OUTER) = &
                       OMEGAMOMS_ELASTIC_UNSCALED(1:NLAYERS,0:MAX_MOMENTS,1:NPOINTS_OUTER)

!   -- Rob mod 5/12/17 for 2p5a, Introducing the phase function product input

      LRRS_FixIn%Atmos%OMEGAPHASFUNC_ELASTIC_UP = OMEGAPHASFUNC_ELASTIC_UP
      LRRS_FixIn%Atmos%OMEGAPHASFUNC_ELASTIC_DN = OMEGAPHASFUNC_ELASTIC_DN

!  Modified Atmosphere Inputs
!mick fix 7/20/2016 - define GEOMETRY_SPECHEIGHT

      LRRS_ModIn%MAtmos%GEOMETRY_SPECHEIGHT = HEIGHT_GRID(NLAYERS)

!  Fixed Surface Inputs

!  Note: LRRS surface arrays initialized elsewhere;
!        appropriate portions just re-defined here

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        LRRS_FixIn%Surf%ALBEDOS_RANKED(1:NPOINTS_OUTER) = ALBEDO
      ELSE
        DO I=1,NPOINTS_OUTER
          IF (NPOINTS_SURF .EQ. 1) THEN
            J = 1
          ELSE
            J = I
          ENDIF
          LRRS_Sup%BRDF%EXACTDB_BRDFUNC     (1:NUSERS,1:NAZIMS,I)         = &
            BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC (1:NUSERS,1:NAZIMS,J)
          LRRS_Sup%BRDF%BRDF_F_0            (0:NMOMS,1:NDISOS,I)          = &
            BRDF_Sup_Out%BS_BRDF_F_0        (0:NMOMS,1:NDISOS,J)
          LRRS_Sup%BRDF%BRDF_F              (0:NMOMS,1:NDISOS,1:NDISOS,I) = &
            BRDF_Sup_Out%BS_BRDF_F          (0:NMOMS,1:NDISOS,1:NDISOS,J)
          LRRS_Sup%BRDF%USER_BRDF_F_0       (0:NMOMS,1:NUSERS,I)          = &
            BRDF_Sup_Out%BS_USER_BRDF_F_0   (0:NMOMS,1:NUSERS,J)
          LRRS_Sup%BRDF%USER_BRDF_F         (0:NMOMS,1:NUSERS,1:NDISOS,I) = &
            BRDF_Sup_Out%BS_USER_BRDF_F     (0:NMOMS,1:NUSERS,1:NDISOS,J)
        ENDDO
      ENDIF

!  Linearized Atmosphere Inputs
!mick mod 11/28/2018 - limited dims on L_OMEGAMOMS_ELASTIC_UNSCALED

      LRRS_LinIn%Cont%N_TOTALCOLUMN_WFS = N_TOTALCOLUMN_WFS

      LRRS_LinIn%Cont%LAYER_VARY_FLAG   = LAYER_VARY_FLAG
      LRRS_LinIn%Cont%LAYER_VARY_NUMBER = LAYER_VARY_NUMBER

      LRRS_LinIn%Optical%LAYER_AIRCOLUMNS_DT    = LAYER_AIRCOLUMNS_DT
      LRRS_LinIn%Optical%TEMPERATURES_UNSHIFTED = TEMPERATURES_UNSHIFTED

      LRRS_LinIn%Optical%L_DELTAU_INPUT_UNSCALED      = L_DELTAU_INPUT_UNSCALED
      !LRRS_LinIn%Optical%L_OMEGAMOMS_ELASTIC_UNSCALED = L_OMEGAMOMS_ELASTIC_UNSCALED
      LRRS_LinIn%Optical%L_OMEGAMOMS_ELASTIC_UNSCALED(:,1:NLAYERS,0:MAX_MOMENTS,1:NPOINTS_OUTER) = &
                         L_OMEGAMOMS_ELASTIC_UNSCALED(:,1:NLAYERS,0:MAX_MOMENTS,1:NPOINTS_OUTER)

!   -- Rob mod 5/12/17 for 2p5a, Introducing the phase function product input

      LRRS_LinIn%Optical%L_OMEGAPHASFUNC_ELASTIC_UP = L_OMEGAPHASFUNC_ELASTIC_UP
      LRRS_LinIn%Optical%L_OMEGAPHASFUNC_ELASTIC_DN = L_OMEGAPHASFUNC_ELASTIC_DN

!  End

      RETURN
      END SUBROUTINE LRRS_IOPSETUP_LPBIN_M

!  Finish module

      END MODULE LPbin_M_iopsetup_m
