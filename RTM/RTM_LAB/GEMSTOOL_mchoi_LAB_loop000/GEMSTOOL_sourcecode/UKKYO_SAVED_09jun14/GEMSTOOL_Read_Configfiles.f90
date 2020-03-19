module GEMSTOOL_Read_Configfiles_m

use GEMSTOOL_Input_Types_m

private
public  GEMSTOOL_Read_Configfiles

contains

subroutine GEMSTOOL_Read_Configfiles &
    ( do_wavnums, ConfigPath, GEMSTOOL_INPUTS, fail, message   )

   implicit none

!  Control for using wavenumbers

   logical, intent(in)       :: do_wavnums

!  path for configuration files

   CHARACTER*(*), intent(in) :: ConfigPath

! GEOMTOOL inputs, structure. Intent(out) here

   TYPE(GEMSTOOL_Config_Inputs), INTENT(INOUT) :: GEMSTOOL_INPUTS

!  Excpetion handling
   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variable

   integer       :: g, k, kdum
   character*100 :: filnam, filnam_save
   data filnam_save /'                                                                                                    '/

!  Tasks as follows

!  NOTES
!  =====

!  1. Aerosol Mie or T-matrix calculations will be done at every 10 nm
!     Range of Mie/Tmat points is determined automatically from the wavelength input

!  Initialize all variables
!  ========================

   CALL GEMSTOOL_Initialize_ConfigInputs ( GEMSTOOL_INPUTS )
   fail = .false. ; message = ' '

!  FILE READS
!  ==========

!  1. Atmospheric constituents
!  ---------------------------

   filnam = filnam_save
   filnam = adjustl(trim(ConfigPath))//'GEMSTOOL_Atmosphere.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_Tracegases         ! turn on, or off do trace gases
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_aerosols           ! turn on, or off do aerosols calculation
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols       ! turn on, or off do aerosols Mie calculation
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols      ! turn on, or off do aerosols T-matrix calculation
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_clouds             ! turn on, or off General Cloud flag
!   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_fileread_wavenums             ! turn on special spectrum
   CLOSE(1)

!  Checking........

   if ( GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols.and.GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols) then
      message = 'Cannot have Mie and Tmatrix aerosols, error in GEMSTOOL_Atmosphere.cfg'
      fail = .true. ;  return
   endif

   if ( GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols.and. &
              (.not.GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols.and..not.GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols)) then
      message = 'Neither Mie and Tmatrix aerosolare specified, even though you have aerosols. error in GEMSTOOL_Atmosphere.cfg'
      fail = .true. ;  return
   endif

!  2. Wavelengths or wavenumbers
!  -----------------------------

   filnam = filnam_save
   if ( do_wavnums ) then
!      if ( GEMSTOOL_INPUTS%Atmosph%do_fileread_wavenums ) then
!........... READ YOUR OWN WAVENUMBERS from ASCII file
!      else
         filnam = adjustl(trim(ConfigPath))//'GEMSTOOL_Wavenums.cfg'
         OPEN(1,file=trim(filnam),err=90,status='old')
         read(1,*)GEMSTOOL_INPUTS%Wavenums%wavnum_start        ! Initial wavenumber
         read(1,*)GEMSTOOL_INPUTS%Wavenums%wavnum_finish       ! Final wavenumber
         read(1,*)GEMSTOOL_INPUTS%Wavenums%wavnum_res          ! Resolution
         CLOSE(1)
!      endif
   else
      filnam = adjustl(trim(ConfigPath))//'GEMSTOOL_Lambdas.cfg'
      OPEN(1,file=trim(filnam),err=90,status='old')
      read(1,*)GEMSTOOL_INPUTS%Lambdas%lambda_start        ! Initial wavelength
      read(1,*)GEMSTOOL_INPUTS%Lambdas%lambda_finish       ! Final wavelength
      read(1,*)GEMSTOOL_INPUTS%Lambdas%lambda_wavres       ! Resolution
      CLOSE(1)
   endif

!  3. Radiative transfer control
!  -----------------------------

   filnam = filnam_save
   filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_RTControl.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')
   read(1,*)GEMSTOOL_INPUTS%RTControl%NVlidort_nstreams     ! Number of Discrete Ordinates in VLIDORT
   read(1,*)GEMSTOOL_INPUTS%RTControl%NVlidort_nstokes      ! number of stokes parameters in VLIDORT  
   read(1,*)GEMSTOOL_INPUTS%RTControl%do_Upwelling          ! turn on flag for Upwelling
   read(1,*)GEMSTOOL_INPUTS%RTControl%do_PathRadiance       ! turn on flag for "Path Radiance"      (TBD)
   read(1,*)GEMSTOOL_INPUTS%RTControl%do_SphericalAlbedo    ! turn on flag for "Spherical Albedo"   (TBD)
   read(1,*)GEMSTOOL_INPUTS%RTControl%do_2wayTransmittance  ! turn on flag for "2-wayTransmittance" (TBD)
   read(1,*)GEMSTOOL_INPUTS%RTControl%do_firstorder_option  ! turn on, using new FO code
   read(1,*)GEMSTOOL_INPUTS%RTControl%FO_do_regular_ps      ! turn on, using FO code regular PS mode
   close(1)

!  4. Time/Location
!  ----------------

   filnam = filnam_save
   filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_TimePosition.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')
   READ(1,*)GEMSTOOL_INPUTS%TimePos%latitude
   READ(1,*)GEMSTOOL_INPUTS%TimePos%longitude
   READ(1,*)GEMSTOOL_INPUTS%TimePos%earthradius
   READ(1,*)GEMSTOOL_INPUTS%TimePos%year
   READ(1,*)GEMSTOOL_INPUTS%TimePos%month
   READ(1,*)GEMSTOOL_INPUTS%TimePos%day_of_month
   CLOSE(1)

!  5. Trace gas information
!  ------------------------

!  example--------
!        2
!        O3  T T uu.prf                       Brion_malicet_daumont
!        NO2 T T no2.prf                      sao_no2_xec

   if ( GEMSTOOL_INPUTS%Atmosph%do_Tracegases ) then
      filnam = filnam_save
      filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_TraceGases.cfg'
      OPEN(1,file=trim(filnam),err=90,status='old')
      read(1,*)GEMSTOOL_INPUTS%Tracegas%ngases
      do g = 1, GEMSTOOL_INPUTS%Tracegas%ngases
         read(1,'(A4,2L2,1X,A)')GEMSTOOL_INPUTS%Tracegas%which_gases(g), &
                                GEMSTOOL_INPUTS%Tracegas%do_gases(g),    &
                                GEMSTOOL_INPUTS%Tracegas%do_gas_wfs(g),  &
                                GEMSTOOL_INPUTS%Tracegas%Gas_Profile_names(g)
      enddo
      CLOSE(1)
   endif

!  6. Aerosol Loading
!  ------------------

   if (  GEMSTOOL_INPUTS%Atmosph%do_aerosols ) then
      if ( GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols .or. GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols) then

         filnam = filnam_save
         filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_AerosolLoading.cfg'
         OPEN(1,file=trim(filnam),err=90,status='old')

!  Loading control

         read(1,*)GEMSTOOL_INPUTS%AerLoad%Loading_DoSurfboundary              ! Flag for using surface as lower boundary
         read(1,*)GEMSTOOL_INPUTS%AerLoad%Loading_upperboundary               ! Upper height limit of aerosol loading [km]
         read(1,*)GEMSTOOL_INPUTS%AerLoad%Loading_lowerboundary               ! Lower height limit of aerosol loading [km]
         read(1,*)GEMSTOOL_INPUTS%AerLoad%Loading_DoFinegridding              ! Flag for using Fine gridding
         read(1,*)GEMSTOOL_INPUTS%AerLoad%Loading_Finegridding                ! Approximate Height resolution [km] for fine-gridding

! Loading case ( 1 = uniform, 2 = Exponential, 3 = GDF )

         read(1,*)GEMSTOOL_INPUTS%AerLoad%Loading_case                        ! Type of Loading

!  AOD and eference wavelength

         read(1,*)GEMSTOOL_INPUTS%AerLoad%aertau_input_w0                     ! Total TAU loading at reference wavelength w0
         read(1,*)GEMSTOOL_INPUTS%AerLoad%reference_w0                        ! reference wavelength w0 [nm]

!  Parameters for the Exponential and GDF profiles

         read(1,*)GEMSTOOL_INPUTS%AerLoad%exploading_relaxation               ! Relaxation parameter for exponential loading (case 2)
         read(1,*)GEMSTOOL_INPUTS%AerLoad%gdfloading_peakheight               ! Peak height [km] parameter for GDF loading (case 3)
         read(1,*)GEMSTOOL_INPUTS%AerLoad%gdfloading_halfwidth                ! Half width  [km] parameter for GDF loading (case 3)

!  End aerosol loading file-read

         close(1)
      endif
   endif

!  7. Mie aerosol control
!  ----------------------

   if ( GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols ) then
      filnam = filnam_save
      filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_MieAerosols.cfg'
      OPEN(1,file=trim(filnam),err=90,status='old')

 !  Bimodal control

      read(1,*)GEMSTOOL_INPUTS%MieTmat%do_bimodal                          ! Flag for Bimodal aerosol
      read(1,*)GEMSTOOL_INPUTS%MieTmat%bimodal_fraction                    ! Bimodal fraction

!  PSD parameters

      read(1,*)GEMSTOOL_INPUTS%MieTmat%PSDindex(1),  GEMSTOOL_INPUTS%MieTmat%PSDindex(2)  ! Mie Particle Size Distribution (PSD) indices
      read(1,*)GEMSTOOL_INPUTS%MieTmat%PSDpars(1,1), GEMSTOOL_INPUTS%MieTmat%PSDpars(1,2) ! Mie PSD First  parameter values
      read(1,*)GEMSTOOL_INPUTS%MieTmat%PSDpars(2,1), GEMSTOOL_INPUTS%MieTmat%PSDpars(2,2) ! Mie PSD Second parameter values
      read(1,*)GEMSTOOL_INPUTS%MieTmat%PSDpars(3,1), GEMSTOOL_INPUTS%MieTmat%PSDpars(3,2) ! Mie PSD Third  parameter values

!  PSD Integration control

      read(1,*)GEMSTOOL_INPUTS%MieTmat%nblocks(1),    GEMSTOOL_INPUTS%MieTmat%nblocks(2)           ! Mie Number of blocks covering the PSDs
      read(1,*)GEMSTOOL_INPUTS%MieTmat%nweights(1),   GEMSTOOL_INPUTS%MieTmat%nweights(2)          ! Mie Number of Quadrature weights in each block
      read(1,*)GEMSTOOL_INPUTS%MieTmat%xparticle_limit                              ! Mie Particle size limit
      read(1,*)GEMSTOOL_INPUTS%MieTmat%FixR1R2(1),    GEMSTOOL_INPUTS%MieTmat%FixR1R2(2)           ! Mie Use particle radius cutoff to fix Rmin/Rmax
      read(1,*)GEMSTOOL_INPUTS%MieTmat%R1R2_cutoff(1),GEMSTOOL_INPUTS%MieTmat%R1R2_cutoff(2)       ! Mie Particle radius cutoff 
      read(1,*)GEMSTOOL_INPUTS%MieTmat%R1(1), GEMSTOOL_INPUTS%MieTmat%R1(2)                        ! Mie Minimum particle radii (Microns)
      read(1,*)GEMSTOOL_INPUTS%MieTmat%R2(1), GEMSTOOL_INPUTS%MieTmat%R2(2)                        ! Mie/Maximum particle radii (Microns)

!  refractive indices

      read(1,*)GEMSTOOL_INPUTS%MieTmat%nreal(1), GEMSTOOL_INPUTS%MieTmat%nreal(2)      ! Mie Real parts of refractive index
      read(1,*)GEMSTOOL_INPUTS%MieTmat%nimag(1), GEMSTOOL_INPUTS%MieTmat%nimag(2)      ! Mie Imaginary part of refractive index

! End Mie aerosols

      close(1)
   endif

!  8. T-Matrix aerosol control
!  ---------------------------

   if ( GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols ) then
      filnam = filnam_save
      filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_TmatAerosols.cfg'
      OPEN(1,file=trim(filnam),err=90,status='old')

 !  Bimodal control

      read(1,*)GEMSTOOL_INPUTS%MieTmat%do_bimodal                          ! Flag for Bimodal aerosol
      read(1,*)GEMSTOOL_INPUTS%MieTmat%bimodal_fraction                    ! Bimodal fraction

!  Orientation control

      read(1,*)GEMSTOOL_INPUTS%MieTmat%Do_EqSaSphere                          ! Tmat flag for equivalent-surface-sphere representation
      read(1,*)GEMSTOOL_INPUTS%MieTmat%Tmat_Sphtype                           ! Tmat Type of Spheroidal particle
      read(1,*)GEMSTOOL_INPUTS%MieTmat%Tmat_ndgs(1),GEMSTOOL_INPUTS%MieTmat%Tmat_ndgs(2)              ! Tmat Gaussian surface numbers
      read(1,*)GEMSTOOL_INPUTS%MieTmat%Tmat_accuracy                          ! Tmat Accuracy for the computation

!  PSD parameters

      read(1,*)GEMSTOOL_INPUTS%MieTmat%PSDindex(1),  GEMSTOOL_INPUTS%MieTmat%PSDindex(2)  ! Tmat Particle Size Distribution (PSD) indices
      read(1,*)GEMSTOOL_INPUTS%MieTmat%PSDpars(1,1), GEMSTOOL_INPUTS%MieTmat%PSDpars(1,2) ! Tmat PSD First  parameter values
      read(1,*)GEMSTOOL_INPUTS%MieTmat%PSDpars(2,1), GEMSTOOL_INPUTS%MieTmat%PSDpars(2,2) ! Tmat PSD Second parameter values
      read(1,*)GEMSTOOL_INPUTS%MieTmat%PSDpars(3,1), GEMSTOOL_INPUTS%MieTmat%PSDpars(3,2) ! Tmat PSD Third  parameter values

!  PSD Integration control

      read(1,*)GEMSTOOL_INPUTS%MieTmat%Tmat_nkmax(1),GEMSTOOL_INPUTS%MieTmat%Tmat_nkmax(2)          ! Tmat Gaussian PSD numbers
      read(1,*)GEMSTOOL_INPUTS%MieTmat%FixR1R2(1),    GEMSTOOL_INPUTS%MieTmat%FixR1R2(2)           ! Tmat Use particle radius cutoff to fix Rmin/Rmax
      read(1,*)GEMSTOOL_INPUTS%MieTmat%R1R2_cutoff(1),GEMSTOOL_INPUTS%MieTmat%R1R2_cutoff(2)       ! Tmat Particle radius cutoff 
      read(1,*)GEMSTOOL_INPUTS%MieTmat%R1(1), GEMSTOOL_INPUTS%MieTmat%R1(2)                        ! Tmat Minimum particle radii (Microns)
      read(1,*)GEMSTOOL_INPUTS%MieTmat%R2(1), GEMSTOOL_INPUTS%MieTmat%R2(2)                        ! Tmat/Maximum particle radii (Microns)

!  refractive indices, shape factor

      read(1,*)GEMSTOOL_INPUTS%MieTmat%nreal(1), GEMSTOOL_INPUTS%MieTmat%nreal(2)      ! Tmat Real parts of refractive index
      read(1,*)GEMSTOOL_INPUTS%MieTmat%nimag(1), GEMSTOOL_INPUTS%MieTmat%nimag(2)      ! Tmat Imaginary part of refractive index
      read(1,*)GEMSTOOL_INPUTS%MieTmat%Tmat_eps(1),GEMSTOOL_INPUTS%MieTmat%Tmat_eps(2)             ! Tmat shape factors

! End T-matrix aerosols

      close(1)
   endif

!  9. Cloud control
!  ----------------

!  Reflecting Cloud (do_Scattering_Clouds = .FALSE. )
!     Cloud_method = 1. Cloudy-sky calculation, Lambertian cloud-top. Reduced number of layers to CTOP
!                       Clear-sky  calculation  with 1 extra layer. Uses CTOP_PRESSURE

!  Scattering Cloud (do_Scattering_Clouds = .TRUE. )
!     Cloud_method = 2. Cloudy-sky calculation, Mie scattering layer. Extra number of layers to CTOP/CBOT
!                       Clear-sky  calculation  with 2 extra layers defined by CTOP/CBOT

   if ( GEMSTOOL_INPUTS%Atmosph%do_clouds ) then
      filnam = filnam_save
      filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_Clouds.cfg'
      OPEN(1,file=trim(filnam),err=90,status='old')

!  Part 1. Cloud Physical Control
!  ------------------------------

      read(1,*)  ; read(1,*) ; read(1,*) ! header lines

      read(1,*)GEMSTOOL_INPUTS%Clouds%do_Scattering_Clouds             ! See above

      read(1,*)GEMSTOOL_INPUTS%Clouds%cloud_fraction                   ! Cloud fraction for IPA calculation (Methods 1 and 2)

      read(1,*)GEMSTOOL_INPUTS%Clouds%Cloud_boundaries                 ! Flag for using level boundaries for cloud input
      read(1,*)GEMSTOOL_INPUTS%Clouds%ctop_level                       ! cloud-top    level boundary            Methods 1 and 2
      read(1,*)GEMSTOOL_INPUTS%Clouds%cbot_level                       ! cloud-bottom level boundary            Only for method 2

      read(1,*)GEMSTOOL_INPUTS%Clouds%cloud_HeightOrPressure           ! Flag for choosing height or pressure input (T = height)
      read(1,*)GEMSTOOL_INPUTS%Clouds%ctop_height                      ! cloud-top    height   [km]             Methods 1 and 2
      read(1,*)GEMSTOOL_INPUTS%Clouds%cbot_height                      ! cloud-bottom height   [km]             Only for method 2
      read(1,*)GEMSTOOL_INPUTS%Clouds%ctop_pressure                    ! cloud-top pressure    [mb]             Methods 1 and 2
      read(1,*)GEMSTOOL_INPUTS%Clouds%cbot_pressure                    ! cloud-bottom pressure [mb]             Only for method 2

      read(1,*)GEMSTOOL_INPUTS%Clouds%cldtau_input_w1                  ! cloud optical depth at wavelength w1   Only for method 2
      read(1,*)GEMSTOOL_INPUTS%Clouds%reference_w1                     ! reference wavelength w1 [nm]           Only for method 2

      read(1,*)GEMSTOOL_INPUTS%Clouds%cloud_albedo                     ! Method  (Reflecting clouds) specify cloud albedo (~0.8)

      read(1,*)GEMSTOOL_INPUTS%Clouds%do_Henyey_Greenstein             ! Option for use of H-G phase function (method 2 only)

!  Part 2. Cloud Mie program control
!  ---------------------------------

      read(1,*)  ; read(1,*) ; read(1,*) ! header lines

      if ( GEMSTOOL_INPUTS%Clouds%do_Scattering_CLouds ) then
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_PSDIndex           ! Method 2 (Scattering Clouds), Distribution index
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_PSDpars(1)         ! Method 2 (Scattering Clouds), First  parameter value
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_PSDpars(2)         ! Method 2 (Scattering Clouds), Second parameter value
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_PSDpars(3)         ! Method 2 (Scattering Clouds), Third  parameter value
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_nblocks            ! Method 2 (Scattering Clouds), Number of blocks covering the PSDs
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_nweights           ! Method 2 (Scattering Clouds), Number of Quadrature weights in each block
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_xparticle_limit    ! Method 2 (Scattering Clouds), Particle size limit
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_FixR1R2            ! Method 2 (Scattering Clouds), Use particle radius cutoff to fix Rmin/Rmax
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_R1R2_cutoff        ! Method 2 (Scattering Clouds), Particle radius cutoff 
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_r1, &
                  GEMSTOOL_INPUTS%Clouds%cld_r2                 ! Method 2 (Scattering Clouds), Minimum/Maximum particle radii (Microns)
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_nreal              ! Method 2 (Scattering Clouds), Real part of refractive index
         read(1,*)GEMSTOOL_INPUTS%Clouds%cld_nimag              ! Method 2 (Scattering Clouds), Imaginary part of refractive index
      endif

!  Close file

      CLOSE(1)
   else
      GEMSTOOL_INPUTS%Clouds%do_Scattering_Clouds = .false.
   endif

!  error handling

   if ( GEMSTOOL_INPUTS%Atmosph%do_clouds ) then
      if ( GEMSTOOL_INPUTS%Clouds%cloud_fraction.lt.0.0 .or. GEMSTOOL_INPUTS%Clouds%cloud_fraction.gt.1.0 ) then
         message = 'Cloud fraction out of range, must be between 0 and 1'
         fail = .true. ;  return
      endif
      if ( GEMSTOOL_INPUTS%Clouds%Cloud_boundaries ) then
         if ( GEMSTOOL_INPUTS%Clouds%ctop_level.ge. GEMSTOOL_INPUTS%Clouds%cbot_level ) then
            message = 'Cloud Top must be higher than cloud bottom'
            fail = .true. ;  return
         endif
      endif
   endif

!  10. Albedo CLosure control
!  --------------------------

   filnam = filnam_save
   filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_AlbedoClosure.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')
   read(1,*)GEMSTOOL_INPUTS%Closure%n_closure_bands
   do g = 1, GEMSTOOL_INPUTS%Closure%n_closure_bands
      read(1,*)GEMSTOOL_INPUTS%Closure%closure_start(g),  &
                GEMSTOOL_INPUTS%Closure%closure_finish(g), &
                GEMSTOOL_INPUTS%Closure%closure_ncoeffs(g)
      do k = 1, GEMSTOOL_INPUTS%Closure%closure_ncoeffs(g)
         read(1,*)kdum,GEMSTOOL_INPUTS%Closure%closure_coeffs(g,k)
      enddo
   enddo
   close(1)

!  11. Geometrical inputs
!  ----------------------

!  Default (22 July 2013). Use Observational Geometry option

   filnam = filnam_save
   filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_Geometries.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')
   read(1,*)GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries
   do g = 1, GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries
      read(1,*)GEMSTOOL_INPUTS%Geometry%GEMS_szas(g), &
                GEMSTOOL_INPUTS%Geometry%GEMS_vzas(g), &
                GEMSTOOL_INPUTS%Geometry%GEMS_azms(g)
   enddo
   CLOSE(1)

!  12. Linearization Control
!  -------------------------

!  Introduced, 24 October 2013

   filnam = filnam_save
   filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_LinControl.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_GasProfile_Jacobians
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_AerOpdepProfile_Jacobians
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_Surface_Jacobians
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_Tshift_Jacobian
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_SurfPress_Jacobian
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_normalized_wfoutput
   CLOSE(1)

!  Exception handling

   if ( GEMSTOOL_INPUTS%LinControl%do_GasProfile_Jacobians.or.GEMSTOOL_INPUTS%LinControl%do_AerOpdepProfile_Jacobians ) then
      if (GEMSTOOL_INPUTS%LinControl%do_Tshift_Jacobian.or.GEMSTOOL_INPUTS%LinControl%do_SurfPress_Jacobian ) then
         message = 'Cannot have Column and Profile Jacobians together - change the flags'
         fail = .true. ; return
      endif
   endif

!  Normal finish

   return

!  Error finish

90 continue
   fail    = .true.
   message = trim(filnam)//': File not found'

!  Finish

   return
end subroutine GEMSTOOL_Read_Configfiles

subroutine GEMSTOOL_Initialize_ConfigInputs ( GEMSTOOL_INPUTS )

   implicit none

! GEOMTOOL inputs, structure. Intent(out) here

   TYPE(GEMSTOOL_Config_Inputs), INTENT(INOUT) :: GEMSTOOL_INPUTS

!  1. Atmospheric constituents
!  ---------------------------

   GEMSTOOL_INPUTS%Atmosph%do_Tracegases    = .false.
   GEMSTOOL_INPUTS%Atmosph%do_aerosols      = .false.
   GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols  = .false.
   GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols = .false.
   GEMSTOOL_INPUTS%Atmosph%do_clouds        = .false.

!  2. Wavelengths and Wavenumbers
!  ------------------------------

   GEMSTOOL_INPUTS%Lambdas%lambda_start  = 0.0d0
   GEMSTOOL_INPUTS%Lambdas%lambda_finish = 0.0d0
   GEMSTOOL_INPUTS%Lambdas%lambda_wavres = 0.0d0

   GEMSTOOL_INPUTS%Wavenums%wavnum_start  = 0.0d0
   GEMSTOOL_INPUTS%Wavenums%wavnum_finish = 0.0d0
   GEMSTOOL_INPUTS%Wavenums%wavnum_res    = 0.0d0

!  3. Radiative transfer control
!  -----------------------------

   GEMSTOOL_INPUTS%RTControl%NVlidort_nstreams     = 0
   GEMSTOOL_INPUTS%RTControl%NVlidort_nstokes      = 0
   GEMSTOOL_INPUTS%RTControl%do_Upwelling          = .false.
   GEMSTOOL_INPUTS%RTControl%do_PathRadiance       = .false.
   GEMSTOOL_INPUTS%RTControl%do_SphericalAlbedo    = .false.
   GEMSTOOL_INPUTS%RTControl%do_2wayTransmittance  = .false.
   GEMSTOOL_INPUTS%RTControl%do_firstorder_option  = .false.
   GEMSTOOL_INPUTS%RTControl%FO_do_regular_ps      = .false. 

!  4. Time/Location
!  ----------------

   GEMSTOOL_INPUTS%TimePos%latitude     = 0.0d0
   GEMSTOOL_INPUTS%TimePos%longitude    = 0.0d0
   GEMSTOOL_INPUTS%TimePos%earthradius  = 0.0d0
   GEMSTOOL_INPUTS%TimePos%year         = 0
   GEMSTOOL_INPUTS%TimePos%month        = 0
   GEMSTOOL_INPUTS%TimePos%day_of_month = 0

!  5. Trace gas information
!  ------------------------

   GEMSTOOL_INPUTS%Tracegas%ngases      = 0
   GEMSTOOL_INPUTS%Tracegas%which_gases = '    '
   GEMSTOOL_INPUTS%Tracegas%do_gases    = .false.
   GEMSTOOL_INPUTS%Tracegas%do_gas_wfs  = .false.
   GEMSTOOL_INPUTS%Tracegas%Gas_Profile_names = ' '

!  AND SO ON

   return
end subroutine GEMSTOOL_Initialize_ConfigInputs

!  End module

end module GEMSTOOL_Read_Configfiles_m
