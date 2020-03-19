
module GEMSTOOL_Input_Types_m

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

   public
   private :: fpk

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_Atmosphere

!  Constituent flags; Configuration Task 1

!  Overall flags

   LOGICAL                      :: DO_AEROSOLS
   LOGICAL                      :: DO_CLOUDS
   LOGICAL                      :: DO_TRACEGASES

!  Aerosol Mie/Tmat flags

   LOGICAL            :: do_Mie_aerosols
   LOGICAL            :: do_Tmat_aerosols

   END TYPE GEMSTOOL_Atmosphere

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_Lambdas

!  wavelength inputs, Configuration Task 2

   real(kind=fpk)      :: lambda_start, lambda_finish, lambda_wavres

   END TYPE GEMSTOOL_Lambdas

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_Wavenums

!  wavenumber inputs, Configuration Task 2

   real(kind=fpk)      :: wavnum_start, wavnum_finish, wavnum_res

   END TYPE GEMSTOOL_Wavenums

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_RTControl

!  RT Control, Configuration Task 3

   integer :: NVlidort_nstreams     ! Number of Discrete Ordinates in VLIDORT
   integer :: NVlidort_nstokes      ! number of stokes parameters in VLIDORT  
   logical :: do_Upwelling          ! turn on flag for Upwelling
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

   END TYPE GEMSTOOL_TraceGas

! #####################################################################
! #####################################################################

   TYPE GEMSTOOL_AerLoad

   Logical      :: Loading_DoSurfboundary              ! Flag for using surface as lower boundary
   real(kind=fpk) :: Loading_upperboundary               ! Upper height limit of aerosol loading [km]
   real(kind=fpk) :: Loading_lowerboundary               ! Lower height limit of aerosol loading [km]
   Logical      :: Loading_DoFinegridding              ! Flag for using Fine gridding
   real(kind=fpk) :: Loading_Finegridding                ! Approximate Height resolution [km] for fine-gridding

   Integer      :: Loading_case                        ! Type of Loading

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

   TYPE GEMSTOOL_MieTmat

!  Bimodal control (Mie or Tmat)

   Logical      :: do_bimodal                          ! Flag for Bimodal aerosol
   real(kind=fpk) :: bimodal_fraction                    ! Bimodal fraction

!  PSD parameters (Mie or Tmat)

   Integer      :: PSDindex(2)   ! Mie/Tmat Particle Size Distribution (PSD) indices
   real(kind=fpk) :: PSDpars(3,2)  ! Mie/Tmat PSD parameter values

!  PSD Integration control (Mie or Tmat)

   Logical      :: FixR1R2(2)           ! Mie/Tmat Use particle radius cutoff to fix Rmin/Rmax
   real(kind=fpk) :: R1R2_cutoff(2)       ! Mie/Tmat Particle radius cutoff 
   real(kind=fpk) :: R1(2)                ! Mie/Tmat Minimum particle radii (Microns)
   real(kind=fpk) :: R2(2)                ! Mie/Tmat Maximum particle radii (Microns)

!  refractive indices (Mie or Tmat)

   real(kind=fpk) :: nreal(2)      ! Mie/Tmat Real parts of refractive index
   real(kind=fpk) :: nimag(2)      ! Mie/Tmat Imaginary part of refractive index

!  Orientation control ( TMAT ONLY )

   logical :: Do_EqSaSphere                          ! Tmat flag for equivalent-surface-sphere representation
   integer :: Tmat_Sphtype                           ! Tmat Type of Spheroidal particle
   integer :: Tmat_ndgs(2)                           ! Tmat Gaussian surface numbers
   real(kind=fpk) :: Tmat_accuracy                     ! Tmat Accuracy for the computation

!  PSD Integration control ( Mie ONLY )

   Integer      :: nblocks(2)           ! Mie Number of blocks covering the PSDs
   Integer      :: nweights(2)          ! Mie Number of Quadrature weights in each block
   real(kind=fpk) :: xparticle_limit      ! Mie Particle size limit

!  PSD Integration control ( TMAT ONLY )

   integer :: Tmat_nkmax(2)          ! Tmat Gaussian PSD numbers

!  shape factor ( TMAT ONLY )

   real(kind=fpk) :: Tmat_eps(2)             ! Tmat shape factors

   END TYPE GEMSTOOL_MieTmat

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

   TYPE GEMSTOOL_LinControl

!  Flags for profile Jacobians

   Logical        :: do_GasProfile_Jacobians
   Logical        :: do_AerOpdepProfile_Jacobians

!  Surface Jacobian

   Logical        :: do_Surface_Jacobians

!  Tshift and surface pressure

   Logical        :: do_Tshift_Jacobian
   Logical        :: do_SurfPress_Jacobian

!  Flag for generating normalized weighting function output

   logical        :: do_normalized_wfoutput

   END TYPE GEMSTOOL_LinControl

! #####################################################################
! #####################################################################
!  And last, the nested structure

  TYPE GEMSTOOL_Config_Inputs

      TYPE(GEMSTOOL_Atmosphere)     :: Atmosph
      TYPE(GEMSTOOL_Lambdas)        :: Lambdas
      TYPE(GEMSTOOL_Wavenums)       :: Wavenums
      TYPE(GEMSTOOL_RTControl)      :: RTControl
      TYPE(GEMSTOOL_TimePos)        :: TimePos
      TYPE(GEMSTOOL_TraceGas)       :: TraceGas
      TYPE(GEMSTOOL_AerLoad)        :: AerLoad
      TYPE(GEMSTOOL_MieTmat)        :: MieTmat
      TYPE(GEMSTOOL_Clouds)         :: Clouds
      TYPE(GEMSTOOL_Closure)        :: Closure
      TYPE(GEMSTOOL_Geometry)       :: Geometry
      TYPE(GEMSTOOL_LinControl)     :: LinControl

  END TYPE GEMSTOOL_Config_Inputs

!  End module

end module GEMSTOOL_Input_Types_m
