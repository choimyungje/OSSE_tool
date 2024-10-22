
module GEMSTOOL_Input_Types_m

!  10/18/16. Added SIF control (new sub-structure GEMSTOOL_SurfaceLeaving)
!            Expanded to include Water-leaving inputs (10/25/16 - see below)

!  10/24/16. Added new substructure for using BRDF options - this is the standard
!            BRDF multi-kernel VLIDORT system, should not be used for Ocean Color.

!  10/25/16. Added Control for Water-leaving inputs in sub-structure GEMSTOOL_SurfaceLeaving
!            This also includes use of an optional BRDF Cox-Munk glitter function.

!  10/26/16. Added new substructure for using Instrument band data. Options
!             limited to Himawari-8 at the moment, and only for UVN tool

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  10/19/16.
!     Number of aerosol modes (Local parameter setting here - could be changed later)
!     2 = Bimodal. Should be the same as in the AerProperties Modules.

   integer, parameter :: maxAerModes = 2

!  10/26/16. Satellite data
!   Inst_maxlambdas = maximum channel wavelengths you think you will need.

   integer, parameter :: Inst_maxlambdas = 5001

   public
   private :: fpk, maxAerModes, Inst_maxlambdas

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_Atmosphere

!  Constituent flags; Configuration Task 1

!  Overall flags

   LOGICAL                      :: DO_AEROSOLS
   LOGICAL                      :: DO_CLOUDS
   LOGICAL                      :: DO_TRACEGASES

!  Aerosol Mie/Tmat flags. User flag added 10/19/16.

   LOGICAL            :: do_Mie_aerosols
   LOGICAL            :: do_Tmat_aerosols
   LOGICAL            :: do_User_aerosols

   END TYPE GEMSTOOL_Atmosphere

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_Lambdas

!  wavelength inputs, Configuration Task 2a

!  Monochromatic flag. 10/26/16.
!   ---- If not set, will do an instrumental band calculation

   LOGICAL             :: do_Monochromatic

!  wavelength start, finish, resolution

   real(kind=fpk)      :: lambda_start, lambda_finish, lambda_wavres

   END TYPE GEMSTOOL_Lambdas

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_Wavenums

!  wavenumber inputs, Configuration Task 2b

   real(kind=fpk)      :: wavnum_start, wavnum_finish, wavnum_res

   END TYPE GEMSTOOL_Wavenums

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_Instruments

!  Instrument/Band inputs, Configuration Task 2c

!  Convolution and Pre-convolution options
!  --- If Pre-convolution set, then Cross-sections are pre-convolved.

   Logical      :: do_Convolution
   Logical      :: do_Preconvolution

!  Instrument and Band choices
!    1 = HIMAWARI-8. Only one implemented so far
!    2 = GOCI-1      TBD
!    3 = GOCI-2      TBD

   integer      :: Instrument_choice
   character*25 :: Instrument_name
   integer      :: Instrument_band

!  Band information.
!   This will be taken from pre-pared file, e.g. HIMAWARI-8

   real(kind=fpk)  :: band_start
   real(kind=fpk)  :: band_finish
   real(kind=fpk)  :: band_median
   real(kind=fpk)  :: band_fineres
   real(kind=fpk)  :: band_extension
   real(kind=fpk)  :: band_slitpars(2)
  
! RSR values.
!   Inst_maxlambdas = maximum channel wavelengths you think you will need.

   integer        :: rsr_nlambda
   real(kind=fpk) :: rsr_ctrwl
   real(kind=fpk) :: rsr_lambda (inst_maxlambdas)
   real(kind=fpk) :: rsr_vals   (inst_maxlambdas)

   END TYPE GEMSTOOL_Instruments

! #####################################################################
! #####################################################################

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_RTControl

!  RT Control, Configuration Task 3

   integer :: NVlidort_nstreams     ! Number of Discrete Ordinates in VLIDORT
   integer :: NVlidort_nstokes      ! number of stokes parameters in VLIDORT  
   logical :: do_Upwelling          ! turn on flag for Upwelling
   ! integer :: observation_height    ! select user level
   real(kind=fpk) :: observation_height    ! select user level
   logical :: do_PathRadiance       ! turn on flag for "Path Radiance"      (TBD)
   logical :: do_SphericalAlbedo    ! turn on flag for "Spherical Albedo"   (TBD)
   logical :: do_2wayTransmittance  ! turn on flag for "2-wayTransmittance" (TBD)
   logical :: do_firstorder_option  ! turn on, using new FO code
   logical :: FO_do_regular_ps      ! turn on, using FO code regular PS mode

   END TYPE GEMSTOOL_RTControl

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_TimePos

!  Geographical coordinates and time, Configuration Task 4

   real(kind=fpk) :: latitude
   real(kind=fpk) :: longitude
   real(kind=fpk) :: earthradius

   integer :: year
   integer :: month
   integer :: day_of_month

!  These variables added for the specification of the SIF Epoch
!    10/18/16. Rob Spurr

   integer :: Hour
   integer :: Minute
   integer :: Second

   END TYPE GEMSTOOL_TimePos

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_TraceGas

!  Up to 8 trace gases assumed here

   integer :: ngases
   character(Len=4)   :: which_gases (8)
   Logical            :: do_gases    (8)
   Logical            :: do_gas_wfs  (8)
   character(Len=256) :: Gas_Profile_names (8)

!  @@ Rob fix 7/23/14, Updated to include H2O scaling flag and value

   Logical            :: do_H2OScaling
   real(kind=fpk)     :: H2OScaling

!  @@ Y.Jung fix 2/1/15, Updated to include CH4 scaling flag and value

   Logical            :: do_CH4Scaling
   real(kind=fpk)     :: CH4Scaling

   END TYPE GEMSTOOL_TraceGas

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_AerLoad

   Logical        :: Loading_DoSurfboundary              ! Flag for using surface as lower boundary
   real(kind=fpk) :: Loading_upperboundary               ! Upper height limit of aerosol loading [km]
   real(kind=fpk) :: Loading_lowerboundary               ! Lower height limit of aerosol loading [km]
   Logical        :: Loading_DoFinegridding              ! Flag for using Fine gridding
   real(kind=fpk) :: Loading_Finegridding                ! Approximate Height resolution [km] for fine-gridding

   Integer        :: Loading_case                        ! Type of Loading

!  AOD and eference wavelength

   real(kind=fpk) :: aertau_input_w0                     ! Total TAU loading at reference wavelength w0
   real(kind=fpk) :: reference_w0                        ! reference wavelength w0 [nm]

!  Parameters for the Exponential and GDF profiles

   real(kind=fpk) :: exploading_relaxation               ! Relaxation parameter for exponential loading (case 2)
   real(kind=fpk) :: gdfloading_peakheight               ! Peak height [km] parameter for GDF loading (case 3)
   real(kind=fpk) :: gdfloading_halfwidth                ! Half width  [km] parameter for GDF loading (case 3)

   END TYPE GEMSTOOL_AerLoad

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_MieTmatUser

!  User-defined aerosol, filenames. Added 10/19/16

   character*256  :: User_Aerofile_names(maxAerModes)

!  Bimodal control (Mie or Tmat)

   Logical        :: do_bimodal                          ! Flag for Bimodal aerosol
   real(kind=fpk) :: bimodal_fraction                    ! Bimodal fraction

!  PSD parameters (Mie or Tmat)

   Integer        :: PSDindex(maxAerModes)   ! Mie/Tmat Particle Size Distribution (PSD) indices
   real(kind=fpk) :: PSDpars(3,maxAerModes)  ! Mie/Tmat PSD parameter values

!  PSD Integration control (Mie or Tmat)

   Logical        :: FixR1R2(maxAerModes)           ! Mie/Tmat Use particle radius cutoff to fix Rmin/Rmax
   real(kind=fpk) :: R1R2_cutoff(maxAerModes)       ! Mie/Tmat Particle radius cutoff 
   real(kind=fpk) :: R1(maxAerModes)                ! Mie/Tmat Minimum particle radii (Microns)
   real(kind=fpk) :: R2(maxAerModes)                ! Mie/Tmat Maximum particle radii (Microns)

!  refractive indices (Mie or Tmat)

   real(kind=fpk) :: nreal(maxAerModes)      ! Mie/Tmat Real parts of refractive index
   real(kind=fpk) :: nimag(maxAerModes)      ! Mie/Tmat Imaginary part of refractive index

!  Orientation control ( TMAT ONLY )

   logical :: Do_EqSaSphere                          ! Tmat flag for equivalent-surface-sphere representation
   integer :: Tmat_Sphtype                           ! Tmat Type of Spheroidal particle
   integer :: Tmat_ndgs(maxAerModes)                           ! Tmat Gaussian surface numbers
   real(kind=fpk) :: Tmat_accuracy                     ! Tmat Accuracy for the computation

!  PSD Integration control ( Mie ONLY )

   Integer      :: nblocks(maxAerModes)           ! Mie Number of blocks covering the PSDs
   Integer      :: nweights(maxAerModes)          ! Mie Number of Quadrature weights in each block
   real(kind=fpk) :: xparticle_limit      ! Mie Particle size limit

!  PSD Integration control ( TMAT ONLY )

   integer :: Tmat_nkmax(maxAerModes)          ! Tmat Gaussian PSD numbers

!  shape factor ( TMAT ONLY )

   real(kind=fpk) :: Tmat_eps(maxAerModes)             ! Tmat shape factors

   END TYPE GEMSTOOL_MieTmatUser

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_Clouds

!  Cloud physical control
!  ----------------------

!  flag for choosing Clouds-as-layers (CAL) = Method 2

   logical        :: do_Scattering_Clouds

!  Effective cloud fraction

   real(kind=fpk) :: cloud_fraction                  ! Cloud fraction for IPA calculation (Methods 1 and 2)

!  Control for using level boundaires to define clouds. SIMPLEST CHOICE
!    If we choose this way, then the Cloud pressures/heights come from the atmopsheric Pressure/height layering

   logical        :: Cloud_boundaries                ! Flag for using level boundaries for cloud input
   integer        :: ctop_level                      ! cloud-top    level boundary            Methods 1 and 2
   integer        :: cbot_level                      ! cloud-bottom level boundary            Only for method 2

!  Control for defining Cloud-top and bottom pressures or heights. MORE COMPLICATED CHOICE

   logical        :: cloud_HeightOrPressure           ! Flag for choosing height or pressure input (T = height)
   real(kind=fpk) :: ctop_height                      ! cloud-top    height   [km]             Methods 1 and 2
   real(kind=fpk) :: cbot_height                      ! cloud-bottom height   [km]             Only for method 2
   real(kind=fpk) :: ctop_pressure                    ! cloud-top pressure    [mb]             Methods 1 and 2
   real(kind=fpk) :: cbot_pressure                    ! cloud-bottom pressure [mb]             Only for method 2

!  Optical depth of cloud at reference wavelength w1

   real(kind=fpk) :: cldtau_input_w1                  ! cloud optical depth at wavelength w1   Only for method 2
   real(kind=fpk) :: reference_w1                     ! reference wavelength w1 [nm]           Only for method 2

!  Cloud albedo

   real(kind=fpk) :: cloud_albedo                     ! Method  (Reflecting clouds) specify cloud albedo (~0.8)

!  Henyey-Greenstein flag

   logical        :: do_Henyey_Greenstein             ! Method 2 only

!  Cloud Mie calculation inputs
!  ----------------------------

!  PSD parameters

   Integer        :: CLD_PSDindex   ! Particle Size Distribution (PSD) indices
   real(kind=fpk) :: CLD_PSDpars(3) ! PSD parameter values

!  PSD Integration control

   Logical        :: CLD_FixR1R2           ! Mie/Tmat Use particle radius cutoff to fix Rmin/Rmax
   real(kind=fpk) :: CLD_R1R2_cutoff       ! Mie/Tmat Particle radius cutoff 
   real(kind=fpk) :: CLD_R1                ! Mie/Tmat Minimum particle radii (Microns)
   real(kind=fpk) :: CLD_R2                ! Mie/Tmat Maximum particle radii (Microns)

!  refractive index

   real(kind=fpk) :: CLD_nreal      ! Real part of refractive index
   real(kind=fpk) :: CLD_nimag      ! Imaginary part of refractive index

!  PSD Integration control

   Integer        :: CLD_nblocks              ! Mie Number of blocks covering the PSDs
   Integer        :: CLD_nweights             ! Mie Number of Quadrature weights in each block
   real(kind=fpk) :: CLD_xparticle_limit      ! Mie Particle size limit

   END TYPE GEMSTOOL_Clouds

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_BRDF

!  top level flag

   logical        :: do_Gemstool_BRDF

!  Multi-Kernel model choices (up to 4 kernels)

   integer        :: n_kernels
   integer        :: Kernel_indices(4)
   integer        :: Kernel_n_parameters(4)
   real(kind=fpk) :: Kernel_factors(4)
   real(kind=fpk) :: Kernel_parameters(4,3)

!  For reference, here is the list of Kernel indices, with names, and scalar/vector status
!   THESE ARE THE VLIDORT KERNELS

!  Index     Name      Scalar/Vector   # parameters     Remarks

!    1    Lambertian       Scalar          0            MODIS system
!    2    Ross-thin        Scalar          0            MODIS system
!    3    Ross-thick       Scalar          0            MODIS system
!    4    Li-sparse        Scalar          2            MODIS system 
!    5    Li_dense         Scalar          2            MODIS system 
!    6    Hapke            Scalar          3            not recommended
!    7    Roujean          Scalar          0            not recommended
!    8    Rahman           Scalar          3            Also known as RPV
!    9    Cox-Munk         Scalar          3            Ocean Glitter
!    10   GissCoxMunk      Vector          3            Ocean Glitter, polarized
!    11   GCMComplex       Vector          3            Ocean Glitter, polarized with complex RI
!    12   BPDF-Vegn        Vector          1            Land BRDF with polarization
!    13   BPDF-Soil        Vector          1            Land BRDF with polarization
!    14   BPDF-NDVI        Vector          3            Land BRDF with polarization (recommended)
!    15   NewCMGlint       Scalar          3            Ocean Glitter with water-leaving
!    16   NewGCMGlint      Vector          3            Ocean Glitter with water-leaving, polarized

   END TYPE GEMSTOOL_BRDF

! #####################################################################
! #####################################################################

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_Closure

!  Up to 20 Closure bands assumed here

   Integer        :: n_closure_bands
   real(kind=fpk) :: closure_start  (20)
   real(kind=fpk) :: closure_finish (20)
   integer        :: closure_ncoeffs(20)
   real(kind=fpk) :: closure_coeffs (20,5)

   END TYPE GEMSTOOL_Closure

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_Geometry

!  Up to 40 Observational geometry triplets assumed here

   Integer      :: n_GEMS_geometries
   real(kind=fpk) :: GEMS_szas (40)
   real(kind=fpk) :: GEMS_vzas (40)
   real(kind=fpk) :: GEMS_azms (40)

   END TYPE GEMSTOOL_Geometry

! #####################################################################
! #####################################################################

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_SurfaceLeaving

!  Water-leaving control
!  ---------------------

!  water-leaving flag. Mutually exclusive from the SIF flag (this is checked!)

   LOGICAL                      :: DO_WaterLeaving

!  Water-leaving and glitter - if set this will trigger a glint calculation
!                              using the NewCM options in the BRDF supplement

   LOGICAL                      :: DO_OceanColor

!  Chlorophyll concentration in [mg/M] and Salinity in [ppt]

   REAL(fpk)                    :: ChlorConc
   REAL(fpk)                    :: Salinity

!  Windspeed and direction (relative to 1 solar beam)

   REAL(fpk)                    :: WindSpeed
   REAL(fpk)                    :: WindDir

!  Flags for glint shadowing, Foam Correction, facet Isotropy

   LOGICAL                      :: DO_GlintShadow
   LOGICAL                      :: DO_FoamOption
   LOGICAL                      :: DO_FacetIsotropy

!  NOTES - Wavelength will be derived from the UVN Driver
!        - WaterLeaving will only be used with the UVN tool

!  Fluorescence control
!  --------------------

!  SIF flag. (will be assumed isotropic, in the treatment)

   LOGICAL                      :: DO_SIF

!  Approximate method (using exact-only Surface-leaving)

   LOGICAL                      :: DO_SIF_ExactOnly

!  Fluorescence Amplitude-Scaling on the Database value at 755 nm
!    -- This is the basic control variable for VLIDORT
!    -- ( 11/30/16) Only to be used with the Gaussian parameterization

   REAL(fpk)                    :: SIF755_Amplitude

!  Linear parameterization introduced 11/30/16

   LOGICAL                      :: DO_SIF_LinearParam
   REAL(fpk)                    :: SIF755_Scaling_Constant
   REAL(fpk)                    :: SIF755_Scaling_Gradient

!  NOTES - Wavelength will be derived from the UVN Driver
!        - SIF will only be used with the NSW tool

   END TYPE GEMSTOOL_SurfaceLeaving

! #####################################################################
! #####################################################################

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_LinControl

!  Flags for profile Jacobians

   Logical        :: do_GasProfile_Jacobians
   Logical        :: do_AerOpdepProfile_Jacobians

!  Surface Jacobian

   Logical        :: do_Surface_Jacobians

!  Tshift and surface pressure

   Logical        :: do_Tshift_Jacobian
   Logical        :: do_SurfPress_Jacobian

!  Added by Myungje Choi (mchoi) for read Tshift
   Real(fpk)      :: Tshift_read

!  @@ Rob fix 7/23/14, add H2O scaling Jacobian flag

   Logical        :: do_H2OScaling_Jacobian

!  @@ Y.Jung fix 2/1/15, add CH4 scaling Jacobian flag

   Logical        :: do_CH4Scaling_Jacobian

!  ## Rob Fix 10/18/16. Flag for SIF  Jacobian
!             11/30/16. Renamed, now with linear parameterization in SIF

!   Logical        :: do_SIFScaling_Jacobian
   Logical        :: do_SIF_Jacobians

!  Aerosol Bulk-property Jacobians. New, June 2014.
!  -----------------------------------------------

!  Overall flag

   Logical        :: do_AerBulk_Jacobians

!     do_AerBulk_LoadPars_Jacobians --> w.r.t. Profile loading parameters
!     do_AerBulk_RefIndex_Jacobians --> w.r.t. Refractive Index components
!     do_AerBulk_ShapeFac_Jacobians --> w.r.t. Shape Factor (T-matrix only)
!     do_AerBulk_SizeDist_Jacobians --> w.r.t. Size Distribution Parameters
!     do_AerBulk_BmodFrac_Jacobian  --> w.r.t. Bimodal Fraction

   Logical        :: do_AerBulk_LoadPars_Jacobians
   Logical        :: do_AerBulk_RefIndex_Jacobians
   Logical        :: do_AerBulk_ShapeFac_Jacobians
   Logical        :: do_AerBulk_SizeDist_Jacobians
   Logical        :: do_AerBulk_BmodFrac_Jacobian

!  Flag for generating normalized weighting function output

   logical        :: do_normalized_wfoutput

! Flag for using HITRAN

   logical        :: do_HITRAN

   END TYPE GEMSTOOL_LinControl

! #####################################################################
! #####################################################################
!  And last, the nested structure
!   10/18/16. Sleave structure is new
!   10/19/16. MieTmat structure renamed
!   10/23/16. BRDF structure is new
!   10/26/16. Instruments structure is new

  TYPE GEMSTOOL_Config_Inputs

      TYPE(GEMSTOOL_Atmosphere)     :: Atmosph
      TYPE(GEMSTOOL_Lambdas)        :: Lambdas
      TYPE(GEMSTOOL_Wavenums)       :: Wavenums
      TYPE(GEMSTOOL_Instruments)    :: Instruments
      TYPE(GEMSTOOL_RTControl)      :: RTControl
      TYPE(GEMSTOOL_TimePos)        :: TimePos
      TYPE(GEMSTOOL_TraceGas)       :: TraceGas
      TYPE(GEMSTOOL_AerLoad)        :: AerLoad
      TYPE(GEMSTOOL_MieTmatUser)    :: MieTmatUser
      TYPE(GEMSTOOL_Clouds)         :: Clouds
      TYPE(GEMSTOOL_BRDF)           :: BRDF
      TYPE(GEMSTOOL_Closure)        :: Closure
      TYPE(GEMSTOOL_SurfaceLeaving) :: Sleave
      TYPE(GEMSTOOL_Geometry)       :: Geometry
      TYPE(GEMSTOOL_LinControl)     :: LinControl

  END TYPE GEMSTOOL_Config_Inputs

!  End module

end module GEMSTOOL_Input_Types_m
