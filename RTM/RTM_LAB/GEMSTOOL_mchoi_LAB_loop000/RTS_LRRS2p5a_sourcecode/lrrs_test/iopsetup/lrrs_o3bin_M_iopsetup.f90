
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

      MODULE o3bin_M_iopsetup_NoLin_m

!  Auxiliary and workhorse modules

      USE LRRS_PARS_m
      USE LRRS_Iopsetup_WorkRoutines_m, Only : GET_HEIGHTS_A, GET_PTAO3, &
                                               RAYLEIGH_FUNCTION, GET_O3XSEC_GOME_All

      PRIVATE
      PUBLIC :: LRRS_IOPSETUP_O3BIN_M

      CONTAINS

      SUBROUTINE LRRS_IOPSETUP_O3BIN_M ( &
        SOLARFILE, NOO3,                   & !Inputs
        LAMBDA_START, LAMBDA_FINISH,       & !Inputs
        DO_SLITFUNC_WEIGHTING,             & !Inputs
        SLITFUNC_FWHM, SLITFUNC_CUTOFF,    & !Inputs
        NPOINTS_SURF, ALBEDO, BRDF_Sup_In, & !Inputs
        LRRS_FixIn, LRRS_ModIn, LRRS_Sup,  & !InOut
        GASVOD, FAIL, N_MESSAGES, MESSAGES ) !Outputs

!  This subroutine sets up the optical property inputs required for LRRS

!  Version 2.5 setup routine, without Linearization
!    Lambertian albedo only

!  Used modules
!   -- Rob mod 5/12/17 for 2p5a, removed MAX_MOMENTS_INPUT

      USE LRRS_PARS_m, ONLY: MAX_POINTS, MAX_LAYERS !, MAX_MOMENTS_INPUT
      USE LRRS_IO_DEFS_m
      USE LRRS_RAMAN_SPECTROSCOPY_m, Only : LRRS_BUFFERING_BIN

      USE BRDF_SUP_MOD_m

!   -- Rob mod 5/12/17 for 2p5a, Included Accessory Module for phase function calculations
      USE LRRS_Phasfunc_Sup_Accessory_m

      IMPLICIT NONE

!  Input arguments
!  ---------------

      CHARACTER (LEN=*), INTENT(IN) :: SOLARFILE

!  Ozone control

      LOGICAL, INTENT(IN) ::   NOO3

!  Start and finish wavelengths

      REAL(FPK), INTENT(IN) :: LAMBDA_START, LAMBDA_FINISH

!  Slit function options

      LOGICAL, INTENT(IN) ::   DO_SLITFUNC_WEIGHTING
      REAL(FPK), INTENT(IN) :: SLITFUNC_FWHM, SLITFUNC_CUTOFF

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

      INTEGER, PARAMETER  :: MAX_LOCAL_LAYERS = 50

      INTEGER ::   DATA_NHEIGHTS
      INTEGER ::   CLOUDMASK(MAX_LOCAL_LAYERS)
      REAL(FPK) :: DATA_HEIGHTS(MAX_LOCAL_LAYERS)
      REAL(FPK) :: DATA_PRESSURES(MAX_LOCAL_LAYERS)
      REAL(FPK) :: DATA_TEMPERATURES(MAX_LOCAL_LAYERS)
      REAL(FPK) :: HEIGHT(MAX_LOCAL_LAYERS),TEMPR(MAX_LOCAL_LAYERS)
      REAL(FPK) :: PSURF

!  LRRS input variables, filled here
!  ---------------------------------

!  Outer/inner wavelength range, and offset for the inner range
!    Depends on the choice of solar spectrum

      INTEGER ::   OFFSET_INNER
      INTEGER ::   NPOINTS_INNER
      INTEGER ::   NPOINTS_OUTER

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

      REAL(FPK) :: LAMBDAS_RANKED  ( MAX_POINTS )
      REAL(FPK) :: FLUXES_RANKED   ( MAX_POINTS )

!  Atmospheric physical quantities
!  -------------------------------

!  Number of layers

      INTEGER ::   NLAYERS

!  Height grid

      REAL(FPK) :: HEIGHT_GRID ( 0:MAX_LAYERS )

!  Input layer temperatures, must be in deg K

      REAL(FPK) :: LAYER_TEMPERATURES ( MAX_LAYERS )

!  Input layer Air columns, should be in mol/cm^2 or [DU]

      REAL(FPK) :: LAYER_AIRCOLUMNS ( MAX_LAYERS )

!  Total ozone amount

      REAL(FPK) :: ACTUAL_O3COLUMN

!  Profiles

      REAL(FPK) :: LAYER_O3UMKEHRS ( MAX_LAYERS )

!  Atmospheric optical properties
!  ------------------------------

!  Basic Rayleigh data
!     Rayleigh Cross-sections and depolarization ratios

      REAL(FPK) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  Ozone cross-sections

      REAL(FPK) :: O3XSECS ( MAX_POINTS, MAX_LAYERS )

!  Local cloud inputs

      INTEGER, PARAMETER :: MAX_LOCAL_MOMENTS = 350

      INTEGER ::   NMOMS_1, NMOMS_2
      REAL(FPK) :: CLOUD_MOM_1(0:MAX_LOCAL_MOMENTS)
      REAL(FPK) :: CLOUD_MOM_2(0:MAX_LOCAL_MOMENTS)

!  Number of input phase function Legendre moments
!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT
!      INTEGER ::   NMOMENTS_INPUT

!  Local elastic optical properties
!   -- Rob mod 5/12/17 for 2p5a, used MAX_MOMENTS dimension
!mick mod 11/28/2018 - changed MAX_MOMENTS to MAX_LOCAL_MOMENTS in PHASMOMS_TEMPO and
!                      OMEGAMOMS_ELASTIC_UNSCALED for comparisons with LRRS 2p5

      REAL(FPK) :: DELTAU_TEMPO   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: OMEGA_TEMPO    ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: PHASMOMS_TEMPO ( MAX_LAYERS, 0:MAX_LOCAL_MOMENTS, MAX_POINTS )

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

!  Weighting and Rayleigh moment

   REAL(FPK) :: CLD_WT(MAX_POINTS,2), RAY_WT(MAX_POINTS,2), RAYMOM(MAX_POINTS)

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

!  Other local variables
!  ---------------------

!  File name control

      CHARACTER (LEN=256) :: FILENAME

!  Local error handling

      CHARACTER (LEN=120) :: LOCAL_MESSAGE

!  Help variables

      LOGICAL   :: DO_LAMBERTIAN_SURFACE
      LOGICAL   :: USE_STD_INFILES
      INTEGER   :: I, J, L, N, N1, N2, NL, W, IDUM, NC, V

!   -- Rob mod 5/12/17 for 2p5a, Local number of moments and overall cloud flag

      INTEGER   :: LOCAL_NMOMENTS
      LOGICAL   :: OVERALL_CLOUD

      REAL(FPK) :: RAY_OPDEP, GAS_OPDEP, MOL_OPDEP
      REAL(FPK) :: CLD_OPDEP, CLD_SCDEP
      REAL(FPK) :: TOT_OPDEP, TOT_SCDEP, DN1, DN2
      REAL(FPK) :: RHO_1, RHO_2, RHO_A, TEMP, DIFF
      REAL(FPK) :: CLOUD_TAU_1, CLOUD_TAU_2
      REAL(FPK) :: CLOUD_SSA_1, CLOUD_SSA_2
      REAL(FPK) :: ALBEDO_CHOICE, G

      REAL(FPK) :: SIGMA, CRIT_1, CRIT_2

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

      REAL(FPK), PARAMETER :: CO2_PPMV_MIXRATIO = 380.0d0

!  Start code
!  ----------

!  Set local input variables

      DO_BIN_REALIZATION       = LRRS_ModIn%MBool%DO_BIN_REALIZATION
      DO_LAMBERTIAN_SURFACE    = .not.LRRS_FixIn%Bool%DO_BRDF_SURFACE

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

      OMEGAPHASFUNC_ELASTIC_UP = zero
      OMEGAPHASFUNC_ELASTIC_DN = zero

!  Slit function stuff

      IF ( DO_SLITFUNC_WEIGHTING ) THEN
        SIGMA = SLITFUNC_FWHM / DSQRT(DLOG(2.0D0))
        CRIT_1 = DSQRT(- DLOG (SLITFUNC_CUTOFF))
        CRIT_2 = CRIT_1 * SIGMA * 1.0001
      ELSE
        CRIT_2 = 0.0D0
      ENDIF

!  Input data of solar wavelengths and fluxes (USER FILE-READ)
!  ------------------------------------------

      OPEN(1,FILE=Adjustl(Trim(SOLARFILE)),STATUS='OLD')
      READ(1,*)N_INPUT_POINTS
      DO I = 1, N_INPUT_POINTS
        READ(1,*)INPUT_LAMBDAS(I), INPUT_FLUXES(I)
      ENDDO
      CLOSE(1)

!  Bin Realization: Perform buffering (AUTOMATIC)
!  ----------------------------------------------

!  This is not required for the monochromatic mode

!   given the starting and ending wavelengths:
!      - set the inner and outer   windows and the inner window offset
!      - Set the wavelengths and fluxes on the outer window
!      - Get the binning limits for the outer window

      IF ( DO_BIN_REALIZATION ) THEN
        CALL LRRS_BUFFERING_BIN &
           ( MAX_INPUT_POINTS, MAX_POINTS, CRIT_2, .FALSE.,     & !Inputs
             N_INPUT_POINTS, INPUT_LAMBDAS, INPUT_FLUXES,       & !Inputs
             LAMBDA_START, LAMBDA_FINISH,                       & !Inputs
             NPOINTS_INNER, NPOINTS_OUTER, OFFSET_INNER,        & !Outputs
             LAMBDAS_RANKED, FLUXES_RANKED, BINLOWER, BINUPPER, & !Outputs
             FAIL, N_MESSAGES, MESSAGES )                         !Outputs
        IF (FAIL) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = '  These messages from o3bin_M_Iopsetup'
          RETURN
        END IF
      ENDIF

!  Task 5. Get the physics of the problem (USER-DEFINED INPUTS)
!  --------------------------------------

      !USE_STD_INFILES = .TRUE. !use a standard "temp_psurf.prf*" file
      USE_STD_INFILES = .FALSE. !use alternate "input_heights.dat" file
      IF ( USE_STD_INFILES ) THEN

!  Use standard infiles
!  ====================

!  Get the pressures, temperatures and heights from FILE.

       OPEN(1,FILE='lrrs_test/physics_data/ATMOS/temp_psurf.prf_35',STATUS='OLD')
!        OPEN(1,FILE='lrrs_test/physics_data/ATMOS/temp_psurf.prf_17',STATUS='OLD')
        READ(1,*)DATA_NHEIGHTS,psurf
        DO I = 1, DATA_NHEIGHTS
            READ(1,*)HEIGHT(I), TEMPR(I)
        ENDDO
        CLOSE(1)

!  Set the pressures

        DATA_PRESSURES(DATA_NHEIGHTS)=psurf
        DO I = DATA_NHEIGHTS-1, 1 ,-1
          DATA_PRESSURES(I) = DATA_PRESSURES( &
          I+1)*exp(-9.81*28.9/8314.0*(height(DATA_NHEIGHTS-I+1)- &
          height(DATA_NHEIGHTS-I))*1.0E3/2.0*(1.0/tempr(DATA_NHEIGHTS-I+1) &
          +1.0/tempr(DATA_NHEIGHTS-I)))
        ENDDO

!  Set the data heights

        DO I = 1, DATA_NHEIGHTS
          DATA_HEIGHTS(I)=HEIGHT(DATA_NHEIGHTS-I+1)
          DATA_TEMPERATURES(I)=TEMPR(DATA_NHEIGHTS-I+1)
        ENDDO

!  Set the height grid, number of layers, Layer temperatures and Aircolum

        NLAYERS = DATA_NHEIGHTS - 1
        HEIGHT_GRID(0) = DATA_HEIGHTS(1)
        DO I = 1, NLAYERS
          RHO_1 = DATA_PRESSURES(I) / DATA_TEMPERATURES(I)
          RHO_2 = DATA_PRESSURES(I+1) / DATA_TEMPERATURES(I+1)
          TEMP  = 0.5D0*(DATA_TEMPERATURES(I)+DATA_TEMPERATURES(I+1))
          RHO_A = 0.5D0 * CONSTANT * ( RHO_1 + RHO_2 )
          DIFF  = DATA_HEIGHTS(I) - DATA_HEIGHTS(I+1)
          HEIGHT_GRID(I) = DATA_HEIGHTS(I+1)
          LAYER_TEMPERATURES(I) = TEMP
          LAYER_AIRCOLUMNS(I)   = DIFF * RHO_A
        ENDDO

      ELSE

!  Use alternative input files
!  ===========================

!  Get heights

        FILENAME = 'lrrs_test/physics_data/ATMOS/input_heights.dat'
        CALL GET_HEIGHTS_A &
           ( FILENAME, MAX_LAYERS,        & !Inputs
             HEIGHT_GRID,                 & !Inout
             NLAYERS, FAIL, LOCAL_MESSAGE ) !Outputs
        IF ( FAIL ) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = LOCAL_MESSAGE
          RETURN
        ENDIF

!  Get PTH ozone

        ACTUAL_O3COLUMN = 350.0d0
        FILENAME = 'lrrs_test/physics_data/ATMOS/iup_doc_nhmid_sf_o3vmr.dat'
        CALL GET_PTAO3 &
           ( FILENAME, ACTUAL_O3COLUMN,                             & !Inputs
             MAX_LAYERS, NLAYERS, HEIGHT_GRID,                      & !Inputs
             LAYER_TEMPERATURES, LAYER_AIRCOLUMNS, LAYER_O3UMKEHRS, & !Inputs
             FAIL, LOCAL_MESSAGE )                                    !Outputs
        IF ( FAIL ) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = LOCAL_MESSAGE
          RETURN
        ENDIF

      ENDIF

!  set the Rayleigh cross-sections and depolarization ratios
!    (CO2 mixing ratio is set to 380)

      CALL RAYLEIGH_FUNCTION &
           ( MAX_POINTS, CO2_PPMV_MIXRATIO,  & !Inputs
             NPOINTS_OUTER,  LAMBDAS_RANKED, & !Inputs
             RAYLEIGH_XSEC, RAYLEIGH_DEPOL   ) !Outputs

!  Get the ozone cross-sections. All.

      FILENAME = 'lrrs_test/physics_data/XSECS/gome_o3xsecfm_305365.pre'
      CALL GET_O3XSEC_GOME_All &
          ( FILENAME,                                       & !Inputs
            MAX_POINTS, MAX_LAYERS, NPOINTS_OUTER, NLAYERS, & !Inputs
            LAMBDAS_RANKED, LAYER_TEMPERATURES,             & !Inputs
            O3XSECS, FAIL, LOCAL_MESSAGE )                    !Outputs
      IF ( FAIL ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = LOCAL_MESSAGE
        RETURN
      ENDIF

!  Get cloud
!  ---------

!  Set Cloud Mask and count number of cloud-filled layers
!    0 = clear, 1 = Low Cloud, 2 = high Cloud

      N1 = 0
      N2 = 0
      OVERALL_CLOUD = .false.
      DO N = 1, NLAYERS
        NL = NLAYERS  + 1 - N
        CLOUDMASK(N) = 0
!        if(n .eq. nlayers-1) CLOUDMASK(N) = 1
!        if(n .eq. nlayers-4) CLOUDMASK(N) = 1
!        if(nl .eq. 9) CLOUDMASK(N) = 2
!      print *,'hg',HEIGHT_GRID(N-1),HEIGHT_GRID(N)
!        IF ( NL.GT.2.and.NL.LT.6   ) CLOUDMASK(N) = 1
!        IF ( NL.GT.10.and.NL.LT.12 ) CLOUDMASK(N) = 2
        IF ( CLOUDMASK(N) .EQ. 1 ) N1 = N1 + 1
        IF ( CLOUDMASK(N) .EQ. 2 ) N2 = N2 + 1
        OVERALL_CLOUD = ((CLOUDMASK(N).ne.0).or.OVERALL_CLOUD)
      ENDDO
      DN1 = DBLE(N1)
      DN2 = DBLE(N2)

!  Set cloud total optical thickness and single scatter albedo

!      CLOUD_TAU_1 = 40.0D0
!      CLOUD_TAU_1 = 10.0D0
!      CLOUD_TAU_1 = 5.0D0
!      CLOUD_TAU_1 = 2.0D0
!      CLOUD_TAU_1 = 1.0D0
!      CLOUD_TAU_1 = 0.5D0

      CLOUD_TAU_1 = 0.0D0
      CLOUD_TAU_2 = 0.0D0

      CLOUD_SSA_1 = 0.9999999D0
      CLOUD_SSA_2 = 0.99999D0

!  Turn off the mask if no clouds

      IF ( CLOUD_TAU_1.EQ.0.0D0 ) THEN
        DO N = 1, NLAYERS
          CLOUDMASK(N) = 0
        ENDDO
      ENDIF

!  Get Cloud moments (Open file and read). Low cloud

      OPEN(1,FILE='lrrs_test/physics_data/XSECS/C1_pmom_lesha',STATUS='OLD')
      READ(1,*)NMOMS_1
      DO L = 0, NMOMS_1 - 1
        READ(1,*)IDUM, CLOUD_MOM_1(L)
      ENDDO
      CLOSE(1)

!  Get High Cloud moments (Open file and read). Copy low cloud values

!      OPEN(1,FILE='lrrs_test/physics_data/XSECS/cirrus_m550_legmoms.dat',STATUS='OLD')
!      READ(1,*)NMOMS_2
!      DO L = 0, NMOMS_2 - 1
!        READ(1,*)IDUM, CLOUD_MOM_2(L), DUM
!      ENDDO
!      CLOSE(1)
!
      NMOMS_2=NMOMS_1
      DO L = 0, NMOMS_2 - 1
        CLOUD_MOM_2(L)=CLOUD_MOM_1(L)
      ENDDO

!  Temporary fix, Henyey-Greenstein clouds
!      G = 0.85d0
!      CLOUD_MOM_1(0) = 1.0d0
!      DO L = 1, NMOMS_1 - 1
!        CLOUD_MOM_1(L) = DBLE(2*L+1)*G**DBLE(L)
!      ENDDO

!  Choose the number of moments you want
!   -- Rob mod 5/12/17 for 2p5a, NMOMENTS_INPUT is now disabled, replaced by LOCAL_NMOMENTS

!      NMOMENTS_INPUT = 40
!      NMOMENTS_INPUT = 2 * NSTREAMS
!      NMOMENTS_INPUT = 2

!  Temporary fix, Henyey-Greenstein clouds
!  --------------------------------------

      G = 0.85d0 ; NMOMS_1 = 150
      CLOUD_MOM_1(0) = 1.0d0
      DO L = 1, NMOMS_1
        CLOUD_MOM_1(L) = DBLE(2*L+1)*G**DBLE(L)
      ENDDO
      CLOUD_TAU_1 = 0.5D0
      CLOUD_SSA_1 = 0.999D0
      CLOUDMASK(NLAYERS-2) = 1  ! --> cloud in layer 21
      N1 = 1 ; DN1 = DBLE(N1)

!  Define number of moments to compute
!mick note 11/28/2018 - note that here inside "iopsetup" we are computing
!                       NMOMS_1 moments and not just 2*NSTREAMS moments
!                       in order to crosscheck with LRRS 2p5 which uses
!                       0:NMOMS_1 moments in its SS solution; however,
!                       for actually passing to LRRS 2p5a, we only pass
!                       2*NSTREAMS moments since it only requires that
!                       for its MS solution

      LOCAL_NMOMENTS = 2  ! Rayleigh_only value
      OVERALL_CLOUD = .TRUE.
      IF ( OVERALL_CLOUD ) LOCAL_NMOMENTS = NMOMS_1

!  Get the optical properties
!  ==========================

!mick fix 3/31/2011 - initialize all elements of "GASVOD"
!mick fix 7/20/2016 - modified "OMEGA_TEMPO" IF condition

      GASVOD = ZERO ; RAY_WT = zero ; CLD_WT = zero

      DO W = 1, NPOINTS_OUTER

        RAYMOM(W) = ( 1.0d0 - RAYLEIGH_DEPOL(W) ) &
                  / ( 2.0D0 - RAYLEIGH_DEPOL(W) )

        DO N = 1, NLAYERS
          RAY_OPDEP = RAYLEIGH_XSEC(W) * LAYER_AIRCOLUMNS(N)
          GAS_OPDEP = 0.0D0
          IF (.NOT.NOO3) GAS_OPDEP = LAYER_O3UMKEHRS(N) * O3XSECS(W,N)
          GASVOD(W) = GASVOD(W) + GAS_OPDEP
          MOL_OPDEP =  RAY_OPDEP + GAS_OPDEP
          IF ( CLOUDMASK(N).EQ.0 ) THEN
            !No cloud / Rayleigh only
            DELTAU_TEMPO(N,W) = MOL_OPDEP
            OMEGA_TEMPO(N,W)  = RAY_OPDEP / MOL_OPDEP
            !IF( NOO3 .or. LAYER_o3umkehrs(n) .eq. 0.0 ) &
            !    OMEGA_TEMPO(N,W) = 0.99999999D0
            IF ( NOO3 ) THEN
              OMEGA_TEMPO(N,W) = 0.99999999D0
            ELSE IF ( LAYER_O3UMKEHRS(N) .eq. 0.0 ) THEN
              OMEGA_TEMPO(N,W) = 0.99999999D0
            ENDIF
            PHASMOMS_TEMPO(N,0,W) = 1.0D0
            PHASMOMS_TEMPO(N,1,W) = 0.0D0
            PHASMOMS_TEMPO(N,2,W) = RAYMOM(W)
            DO L = 3, LOCAL_NMOMENTS
              PHASMOMS_TEMPO(N,L,W) = 0.0D0
            ENDDO
          ELSE IF ( CLOUDMASK(N).EQ.1 ) THEN
            !Cloud - type #1
            CLD_OPDEP = CLOUD_TAU_1 / DN1
            CLD_SCDEP = CLD_OPDEP * CLOUD_SSA_1
            TOT_OPDEP = CLD_OPDEP + MOL_OPDEP
            TOT_SCDEP = CLD_SCDEP + RAY_OPDEP
            DELTAU_TEMPO(N,W) = TOT_OPDEP
            OMEGA_TEMPO(N,W)  = TOT_SCDEP / TOT_OPDEP
            IF( OMEGA_TEMPO(N,W) .GT. 0.9999999D0 ) &
                OMEGA_TEMPO(N,W)=0.9999999D0
            RAY_WT(W,1) = RAY_OPDEP /  TOT_SCDEP
            CLD_WT(W,1) = CLD_SCDEP /  TOT_SCDEP
            PHASMOMS_TEMPO(N,0,W) = 1.0D0
            PHASMOMS_TEMPO(N,1,W) = CLD_WT(W,1) * CLOUD_MOM_1(1)
            PHASMOMS_TEMPO(N,2,W) = RAY_WT(W,1) * RAYMOM(W) + &
                                    CLD_WT(W,1) * CLOUD_MOM_1(2)
            DO L = 3, LOCAL_NMOMENTS
              PHASMOMS_TEMPO(N,L,W) = CLD_WT(W,1) * CLOUD_MOM_1(L)
            ENDDO
          ELSE IF (CLOUDMASK(N).EQ.2) THEN
            !Cloud - type #2
            CLD_OPDEP = CLOUD_TAU_2 / DN2
            CLD_SCDEP = CLD_OPDEP * CLOUD_SSA_2
            TOT_OPDEP = CLD_OPDEP + MOL_OPDEP
            TOT_SCDEP = CLD_SCDEP + RAY_OPDEP
            DELTAU_TEMPO(N,W) = TOT_OPDEP
            OMEGA_TEMPO(N,W)  = TOT_SCDEP / TOT_OPDEP
            IF( OMEGA_TEMPO(N,W) .GT. 0.99999999D0 ) &
                OMEGA_TEMPO(N,W)=0.99999999D0
            RAY_WT(W,2) = RAY_OPDEP /  TOT_SCDEP
            CLD_WT(W,2) = CLD_SCDEP /  TOT_SCDEP
            PHASMOMS_TEMPO(N,0,W) = 1.0D0
            PHASMOMS_TEMPO(N,1,W) = CLD_WT(W,2) * CLOUD_MOM_2(1)
            PHASMOMS_TEMPO(N,2,W) = RAY_WT(W,2) * RAYMOM(W) + &
                                    CLD_WT(W,2) * CLOUD_MOM_2(2)
            DO L = 3, LOCAL_NMOMENTS
              PHASMOMS_TEMPO(N,L,W) = CLD_WT(W,2) * CLOUD_MOM_2(L)
            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  Setup the elastic scattering

      DO W = 1, NPOINTS_OUTER
        DO N = 1, NLAYERS
          DELTAU_INPUT_UNSCALED(N,W)=DELTAU_TEMPO(N,W)
          DO L = 0, LOCAL_NMOMENTS
            OMEGAMOMS_ELASTIC_UNSCALED (N,L,W) = OMEGA_TEMPO(N,W)*PHASMOMS_TEMPO(N,L,W)
          ENDDO
        ENDDO
      ENDDO

!  Debug check #1

!      DO W = 1, NPOINTS_OUTER
!        !write(444,'(3i4,1p3e20.10)')W,NLAYERS,NMOMENTS_INPUT,&
!        !  Lambdas_ranked(w),fluxes_ranked(w)!,albedos_ranked(w)
!        !DO N = 1, NLAYERS
!        !  write(444,'(2i4,1p4e20.10)')w,n,DELTAU_INPUT_UNSCALED(N,W),&
!        !    (OMEGAMOMS_ELASTIC_UNSCALED (N,L,W),L=0,2)
!        !ENDDO
!
!        !Using "LOCAL_NMOMENTS" here
!        write(444,'(3i4,1p3e20.10)')W,NLAYERS,LOCAL_NMOMENTS,&
!          Lambdas_ranked(w),fluxes_ranked(w)!,albedos_ranked(w)
!        DO N = 1, NLAYERS
!          write(444,'(2i4,1pe20.10)')w,n,DELTAU_INPUT_UNSCALED(N,W)
!          DO L = 0, LOCAL_NMOMENTS
!            write(444,'(i4,1pe20.10)') L,OMEGAMOMS_ELASTIC_UNSCALED(N,L,W)
!          ENDDO
!        ENDDO
!      ENDDO

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

!  Loop over layers 

      DO N = 1, NLAYERS

!  Mask

         NC = CLOUDMASK(N)

!  Construct the Cloud phase functions, if present

         IF ( CLOUDMASK(N).EQ.1 ) THEN
            CALL Phasfunc_Sup_Accessory &
                ( MAX_LOCAL_MOMENTS, DO_UPWELLING, DO_DNWELLING,     & ! Input
                  NMOMS_1, N_USER_ANGLES, N_USER_RELAZMS,            & ! Input
                  THETA_BOA, USER_ANGLES, USER_RELAZMS, CLOUD_MOM_1, & ! Input
                  DO_LEGCALC, DO_SCATANG, COSSCAT_UP, COSSCAT_DN,    & ! Modified input/output
                  LEGCALC_UP, LEGCALC_DN, PHASFUNC_UP, PHASFUNC_DN )   ! Modified I/O and True output
         ELSE IF ( CLOUDMASK(N).EQ.2 ) THEN
            CALL Phasfunc_Sup_Accessory &
                ( MAX_LOCAL_MOMENTS, DO_UPWELLING, DO_DNWELLING,     & ! Input
                  NMOMS_2, N_USER_ANGLES, N_USER_RELAZMS,            & ! Input
                  THETA_BOA, USER_ANGLES, USER_RELAZMS, CLOUD_MOM_2, & ! Input
                  DO_LEGCALC, DO_SCATANG, COSSCAT_UP, COSSCAT_DN,    & ! Modified input/output
                  LEGCALC_UP, LEGCALC_DN, PHASFUNC_UP, PHASFUNC_DN )   ! Modified I/O and True output
         ENDIF

!  Set the overall input, including Rayleigh

         IF ( NC.EQ.0 ) THEN
            !Rayleigh only
            IF ( DO_UPWELLING ) THEN
               DO W = 1, NPOINTS_OUTER
                 DO V = 1, N_GEOMETRIES
                   OMEGAPHASFUNC_ELASTIC_UP(N,V,W) = OMEGA_TEMPO(N,W) * RAYFUNC_UP(V,W)
                 ENDDO
               ENDDO
            ENDIF
            IF ( DO_DNWELLING ) THEN
               DO W = 1, NPOINTS_OUTER
                 DO V = 1, N_GEOMETRIES
                   OMEGAPHASFUNC_ELASTIC_DN(N,V,W) = OMEGA_TEMPO(N,W) * RAYFUNC_DN(V,W)
                 ENDDO
               ENDDO
            ENDIF
         ELSE
            !Rayleigh + Cloud
            IF ( DO_UPWELLING ) THEN
               DO W = 1, NPOINTS_OUTER
                 DO V = 1, N_GEOMETRIES
                   OMEGAPHASFUNC_ELASTIC_UP(N,V,W) = OMEGA_TEMPO(N,W) * &
                      (RAY_WT(W,NC) * RAYFUNC_UP(V,W) &
                     + CLD_WT(W,NC) * PHASFUNC_UP(V))
                 ENDDO
               ENDDO
            ENDIF
            IF ( DO_DNWELLING ) THEN
               DO W = 1, NPOINTS_OUTER
                 DO V = 1, N_GEOMETRIES
                   OMEGAPHASFUNC_ELASTIC_DN(N,V,W) = OMEGA_TEMPO(N,W) * &
                      (RAY_WT(W,NC) * RAYFUNC_DN(V,W) &
                     + CLD_WT(W,NC) * PHASFUNC_DN(V))
                 ENDDO
               ENDDO
            ENDIF
         ENDIF

      ENDDO

!  Debug check #2

!      DO W = 1, NPOINTS_OUTER
!        DO N = 1, NLAYERS
!          DO V = 1, N_GEOMETRIES
!            write(445,'(3i4,1p2e20.10)')W,N,V,&
!              OMEGAPHASFUNC_ELASTIC_UP(N,V,W),OMEGAPHASFUNC_ELASTIC_DN(N,V,W)
!          ENDDO
!        ENDDO
!      ENDDO

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

!  Fixed Control Inputs
!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT
!      LRRS_FixIn%Cont%NMOMENTS_INPUT  = NMOMENTS_INPUT

      LRRS_FixIn%Cont%NLAYERS         = NLAYERS

!  Fixed Spectral Inputs

      LRRS_FixIn%Spect%LAMBDAS_RANKED = LAMBDAS_RANKED
      LRRS_FixIn%Spect%FLUXES_RANKED  = FLUXES_RANKED

!  For binning calculations

      LRRS_FixIn%Spect%NPOINTS_INNER  = NPOINTS_INNER
      LRRS_FixIn%Spect%OFFSET_INNER   = OFFSET_INNER
      LRRS_FixIn%Spect%NPOINTS_OUTER  = NPOINTS_OUTER
      LRRS_FixIn%Spect%BINLOWER       = BINLOWER
      LRRS_FixIn%Spect%BINUPPER       = BINUPPER

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

!  Surface Inputs

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

!  End

      RETURN
      END SUBROUTINE LRRS_IOPSETUP_O3BIN_M

!  End module

      END MODULE o3bin_M_iopsetup_NoLin_m
