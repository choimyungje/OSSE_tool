module GEMSTOOL_RTInitialize_m

!  Version 3. Rob Spurr 10/17/16 - 10/28/16
!   - Upgrade to Version 2.7 of VLIDORT. Minor changes throughout.
!   - Introduced SIF and BRDF flags, set from  GEMSTOOL_Inputs

use GEMSTOOL_Input_Types_m

USE VLIDORT_PARS
USE VLIDORT_IO_DEFS

contains

subroutine GEMSTOOL_VLIDORT_Initialize &
   ( GEMSTOOL_INPUTS, VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, fail, message )

   implicit none

! GEMSTOOL inputs, structure. Intent(in) here

   TYPE(GEMSTOOL_Config_Inputs), INTENT(IN) :: GEMSTOOL_INPUTS

!  VLIDORT settings structure, Intent(out) here

   TYPE(VLIDORT_Fixed_Inputs), INTENT(INOUT)       :: VLIDORT_FixIn
   TYPE(VLIDORT_Modified_Inputs), INTENT(INOUT)    :: VLIDORT_ModIn

!  VLIDORT supplements i/o structure

   TYPE(VLIDORT_Sup_InOut), INTENT(INOUT)          :: VLIDORT_Sup

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local

   integer :: g

!  ##################################################################
!  ##################################################################
!      V L I D O R T    I N I T I A L I Z A T I O N   S E C T I O N
!  ##################################################################
!  ##################################################################

!  Intialize exception handling

   fail = .false.
   message = ' '

!  SECTION 1 (One-time inputs)
!  ===========================

!  Structure 1: Fixed Boolean inputs
!  ----------------------------------

!    -- No full_up single scattering or SS Truncation
!    -- Pseudo-spherical, no plane parallel
!    -- Upwelling OR downwelling, but not both !!!!!!!!!!!!!
!    -- No thermal or surface emission
!    -- Lambertian Surface (Actual Albedo comes from Closure Settings)
!    -- No specialist options
!    -- No surface leaving

!  Set MS-mode operations in VLIDORT, depenidng on first-order choices
!    if using FO code, then VLIDORT operates only in Multiple-scatter mode.....

   if ( GEMSTOOL_INPUTS%RTcontrol%do_firstorder_option ) then
      VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE     = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL      = .TRUE.
   else
      VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE     = .TRUE.
      VLIDORT_FixIn%Bool%TS_DO_SS_EXTERNAL      = .FALSE.
   endif

   VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION   = .FALSE.
   VLIDORT_FixIn%Bool%TS_DO_SSFULL              = .FALSE.

   VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION    = .FALSE.
   VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION    = .FALSE.

   VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL      = .FALSE. !mchoi; default=False

   VLIDORT_FixIn%Bool%TS_DO_UPWELLING  = GEMSTOOL_INPUTS%RTControl%do_Upwelling
   VLIDORT_FixIn%Bool%TS_DO_DNWELLING  = .not.VLIDORT_FixIn%Bool%TS_DO_UPWELLING

!  Need to turn on the Quadrature output for Spherical Albedo calculation

!   if ( GEMSTOOL_INPUTS%RTcontrol%do_SphericalAlbedo ) then
!      VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT         = .TRUE.
!   else
      VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT         = .FALSE.
!   endif

!  Rob Fix 10/25/16. Now controlled by GEMSTOOL BRDF flag
!    LAMBERTIAN_ALBEDO is set inside wavelength loop; zero here.
!    12/07/16. Note the BRDF option will be turned on for water-leaving

   VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO      = 0.0d0
   if ( GEMSTOOL_INPUTS%BRDF%do_Gemstool_BRDF ) then
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE  = .FALSE.
   else
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE  = .TRUE.
   ENDIF

!  These are specialist operations, do not need them

   VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS        = .FALSE.
   VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_1 = .FALSE.
   VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_2 = .FALSE.
   VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_3 = .FALSE.

!  Fluorescence introduced 10/18/16

   if ( GEMSTOOL_INPUTS%Sleave%do_SIF ) then
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING     = .TRUE.
      VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC        = .TRUE.
   else
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING     = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC        = .FALSE.
   endif

!  12/07/16. For water-leaving control, we need
!        BRDF surface (NewCM), so Lambertian flag is OFF
!        Surface leaving settings, non isotropic

   if ( GEMSTOOL_INPUTS%Sleave%DO_WaterLeaving ) then
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE  = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING     = .TRUE.
      VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC        = .FALSE.
   endif

!  Structure 2: Modified Boolean inputs
!  ------------------------------------

!    -- Outgoing (not "nadir") single scattering option (Might be external)
!    -- Solar sources
!    -- Chapman function with no refractive geometry
!    -- No Mean-value output
!    -- No thermal Transonly output
!    -- User-stream value output

!  VLIDORT Single-scattering flags depend on the GEMS RT-control configurations

   if ( GEMSTOOL_INPUTS%RTcontrol%do_firstorder_option ) then
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR    = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING = .FALSE.
   else
      if ( GEMSTOOL_INPUTS%RTcontrol%FO_do_regular_ps ) then
         VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING    = .FALSE.
         VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR       = .TRUE.
      else
         VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR       = .FALSE.
         VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING    = .TRUE.
      endif
   endif

!######## Rob Upgrade 10/17/16 ######################################

!  Flag for using FO code to compute single-scatter solution.
!   Consider using this for GEMSTOOL upgrades

   VLIDORT_ModIn%MBool%TS_DO_FO_CALC = .false.
!   VLIDORT_ModIn%MBool%TS_DO_FO_CALC = .true.

!  WARNING THOUGH, FO Calculation has no surface-leaving !!!!

!####################################################################

   VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = .TRUE.

   VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = .FALSE.
   VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = .TRUE.

!  Rayleigh-only flag depends on Configuration-file input for aerosols and clouds
!  This in turn determines the convtest, scaling and performance flags
!    Rayleigh only if (1) Clear sky or (2) no aerosols and/or clouds as reflectors 

   if ( GEMSTOOL_INPUTS%Atmosph%do_aerosols .or. &
         GEMSTOOL_INPUTS%Atmosph%do_clouds .and. GEMSTOOL_INPUTS%Clouds%do_Scattering_Clouds ) then
      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY = .false.
   else
      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY = .true.
   endif

!  Rayleigh only, no delta-M scaling, no performance enhancements

   IF ( VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY ) THEN
      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST = .false.
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING  = .false.
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = .false.
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = .false.
   ENDIF

!  Rayleigh+Aerosols, delta-M scaling, no performance enhancements
!      Performance enhancements need to used with care......... default, turn them off

   IF ( .not.VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY ) THEN
      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST = .true.
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING  = .true.
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = .false.
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = .false.
   ENDIF

!  Alwasy need this flag

   VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES       = .TRUE.

!  Need to turn on Mean-value (flux) output for Spherical Albedo calculation

   if ( GEMSTOOL_INPUTS%RTcontrol%do_SphericalAlbedo ) then
      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = .TRUE.
   else
      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT    = .FALSE.
   endif

!  Not doing Flux-alone calculation

   VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY          = .FALSE.

!  No thermal

   VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY   = .FALSE.

!  New flag (Observational Geometry). Default - Use Observational geometry.

   VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = .TRUE.

!  Structure 3: Fixed control inputs
!  ---------------------------------

!  Basic control integers
!mick fix 1/24/2017 - define order of Taylor-series expansions

!    -- Taylor ordering parameter (hardwired)
!    -- Number of Stokes parameters  from GEMSTOOL RT control configuration input
!    -- Number of discrete ordinates from GEMSTOOL RT control configuration input

   VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER     = 3
   VLIDORT_FixIn%Cont%TS_NSTOKES          = GEMSTOOL_INPUTS%RTcontrol%NVlidort_nstokes
   VLIDORT_FixIn%Cont%TS_NSTREAMS         = GEMSTOOL_INPUTS%RTcontrol%NVlidort_nstreams

!  Check this input

   if ( VLIDORT_FixIn%Cont%TS_NSTREAMS .gt. MAXSTREAMS ) then
       message = 'Number of discrete ordinates too large, exceeds VLIDORT dimension'
       fail = .true. ; return
   endif

!    -- Number of layers determined from PTH profile.

   VLIDORT_FixIn%Cont%TS_NLAYERS          = 0

!    -- Number of fine-layer subdivisions in single-scattering = 2 (Not FO). DEFAULT VALUE
! Set to default value for internal FO calculation (Version 3, 10/17/16)

!   VLIDORT_FixIn%Cont%TS_NFINELAYERS      = 2    ! Not required for FO calculation
   VLIDORT_FixIn%Cont%TS_NFINELAYERS      = 6 

!    -- No thermal coefficients

   VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = 0

!    -- Accuracy. DEFAULT VALUE - Might want to introduce some configuration input here 

   VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY = 0.0001d0

!    -- not used

   VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS     = 0
   VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF   = 0

!  Structure 4: Modified control inputs
!  ------------------------------------

!    -- Number of expansion coefficients = 2     for Rayleigh-only
!    -- Number of expansion coefficients = Input for Rayleigh+Aerosol (Not known yet, so set to zero)

   VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 0
   IF ( VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY ) THEN
      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 2
   ELSE
      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 0   !  set to zero here, will be determined later on. 
   ENDIF

!  Structure 5: Beam inputs
!  ------------------------

!  Flux Factor
!   -- Set to 1.0 if you want sun-normalized only. Set later with solar spectrum (if used).

   VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR = 1.0d0

!  -- Beam input; set from GEMSTOOL Configuration file - simple copy.

   VLIDORT_ModIn%MSunrays%TS_N_SZANGLES  = GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries
   do g = 1, GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(g)   = GEMSTOOL_INPUTS%Geometry%GEMS_szas(g)
   enddo

!  Structure 6: User Value inputs
!  ------------------------------

!  -- User VZA and AZM input; set from GEMSTOOL Configuration file - simple copies.

   VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS  = GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries
   do g = 1, GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(g)   = GEMSTOOL_INPUTS%Geometry%GEMS_azms(g)
   enddo

   VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES  = GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries
   do g = 1, GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(g)   = GEMSTOOL_INPUTS%Geometry%GEMS_vzas(g)
   enddo

!    -- One Level only, Level will be set later.

   VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS = 1           ! One level
   VLIDORT_ModIn%MUserVal%TS_USER_LEVELS  = zero        ! Zero for now

!    -- Geometry specification height = Bottom of height grid (set later, zeroed here)

   VLIDORT_ModIn%MUserval%TS_GEOMETRY_SPECHEIGHT = zero

!  Observational geometry input. Default = set from GEMSTOOL Configuration file - simple copies.

   VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS  = GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries
   do g = 1, GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(g,1)   = GEMSTOOL_INPUTS%Geometry%GEMS_szas(g)
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(g,2)   = GEMSTOOL_INPUTS%Geometry%GEMS_vzas(g)
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(g,3)   = GEMSTOOL_INPUTS%Geometry%GEMS_azms(g)
   enddo

!  Structure 7: Fixed Chapman function inputs
!  ------------------------------------------

!    -- Set later after PTH profiles have been formed, Zeroed here

   VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID     = zero
   
!    -- These are not required
   
   VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID    = zero
   VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID = zero
   VLIDORT_FixIn%Chapman%TS_FINEGRID         = 0

!    -- earth radius. Default - use the GEMSTOOL Configuration input value
!    -- Refractive geometry parameter (not required)
   
   VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS     = GEMSTOOL_INPUTS%TimePos%earthradius
   VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = zero

!  Structure 8: Fixed optical inputs
!  ---------------------------------

!    -- These are set later, zeroed here

   VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT     = zero
   VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT    = zero
   VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT  = zero

!    -- No thermal values here

   VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT      = zero
   VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT      = zero

!    -- This is set later, zeroed here

   VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO     = zero

!    -- Atmospheric wavelength, set later on, zeroed here. New Version 2.7, 10/17/16
!          --> Wavelength (Microns) as a Diagnostic
!          --> This should always be set, even if not used
!          -->  This MUST BE SET when using VBRDF and/or VSLEAVE supplements 

   VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH     = zero

!    -- not used. No longer present in Version 2.7. Removed, 10/17/16
!   VLIDORT_FixIn%Optical%TS_LTE_DELTAU_VERT_INPUT = zero
!   VLIDORT_FixIn%Optical%TS_LTE_THERMAL_BB_INPUT  = zero

!  Structure 9: Fixed write inputs
!  -------------------------------

!  No debug or write-to-file options

   VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE          = .FALSE.
   VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT          = .FALSE.
   VLIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME    = ' '

   VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO       = .FALSE.
   VLIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME = ' '

   VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER        = .FALSE.
   VLIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME  = ' '

   VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS        = .FALSE.
   VLIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME  = ' '

!  Fixed linearized control inputs
!  -------------------------------

!  Depends upon the Configuration file input (Jacobians for GEMSTOOL)

!  These quantities are not needed

   !VLIDORT_LinFixIn%Cont%TS_DO_SURFBB_LINEARIZATION = .false.
   !VLIDORT_LinFixIn%Cont%TS_DO_LTE_LINEARIZATION    = .false.

   !VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS  = 0
   !VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES     = ' '

   !VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS     = .false.
   !VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS       = 0

!  These quantities are just initialized here

   !VLIDORT_LinFixIn%Cont%TS_DO_SIMULATION_ONLY = .false.
   !VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG    = .false.
   !VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER  = 0
   !VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = 0
   !VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS      = 0
   !VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES    = ' '

!  Modified linearized control inputs
!  ----------------------------------

!  These things are set here
!    -- profile, column and Surface Jacobians, flags set by GEMSTOOL Configuration inputs

!   VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = .false.
!   VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = .false.
!   VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .false.
!   VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .false.
!   VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = .false.

! HERE ARE SOME HINTS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Set the linearization control for profile, column and LTE Jacobians

!   if ( do_JACOBIANS ) then
!      VLIDORT_LinModIn%MCont%TS_do_profile_linearization = .true.
!      VLIDORT_LinModIn%MCont%TS_do_atmos_linearization   = .true.
!      VLIDORT_LinModIn%MCont%TS_do_linearization         = .true.
!      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = 1
!      VLIDORT_LinFixIn%Cont%TS_profilewf_names(1)        = '-Trace Gas Volume Mixing Ratio-'
!      if ( do_T_Jacobians ) then
!         VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS    = 2
!         VLIDORT_LinFixIn%Cont%TS_profilewf_names(2)    = '------Layer Temperatures-------'
!      endif
!      do n = 1, GC_nlayers
!        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(n)   = .true.
!        VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(n) = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
!      enddo
!   else
!      VLIDORT_LinModIn%MCont%TS_do_profile_linearization = .false.
!      VLIDORT_LinModIn%MCont%TS_do_atmos_linearization   = .false.
!      VLIDORT_LinModIn%MCont%TS_do_linearization         = .false.
!      VLIDORT_LinFixIn%Cont%TS_do_simulation_only        = .true.
!   endif

!  Set the linearization control for surface (albedo) linearization

!   if ( do_JACOBIANS ) then
!      VLIDORT_LinModIn%MCont%TS_do_surface_linearization = .true.
!      VLIDORT_LinFixIn%Cont%TS_n_surface_wfs             = 1
!   endif

!  Zero the supplemental variables
!  -------------------------------
!mick fix 1/24/2017 - added VLIDORT linearized SLEAVE input initializations here for
!                     completeness; however, note that they are commented out along with the
!                     others as the initializing of VLIDORT linearized inputs is currently
!                     done in the applicable GEMSTOOL linearized driver

!  Supplemental BRDF and SLEAVE variables

   VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC  = ZERO
   VLIDORT_Sup%BRDF%TS_BRDF_F_0         = ZERO
   VLIDORT_Sup%BRDF%TS_BRDF_F           = ZERO
   VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0    = ZERO
   VLIDORT_Sup%BRDF%TS_USER_BRDF_F      = ZERO

   !VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC = ZERO
   !VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0        = ZERO
   !VLIDORT_LinSup%BRDF%TS_LS_BRDF_F          = ZERO
   !VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0   = ZERO
   !VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F     = ZERO

   VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC   = ZERO
   VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES  = ZERO
   VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0         = ZERO
   VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0    = ZERO

   !VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC  = ZERO
   !VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES = ZERO
   !VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0        = ZERO
   !VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0   = ZERO

!  Emissivity is zeroed here

   VLIDORT_Sup%BRDF%TS_EMISSIVITY       = ZERO
   VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY  = ZERO

   !VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY      = ZERO
   !VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY = ZERO

!  Done

   return
end subroutine GEMSTOOL_VLIDORT_Initialize

!  End module

end module GEMSTOOL_RTInitialize_m

