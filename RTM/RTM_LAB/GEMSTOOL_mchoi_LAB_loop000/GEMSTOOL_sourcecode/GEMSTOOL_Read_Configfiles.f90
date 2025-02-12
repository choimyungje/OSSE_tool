module GEMSTOOL_Read_Configfiles_m

!  Rob Fix, 10/17/16, 10/19/16, 10/20/16, 10/24/16, 10/25/16
!    -- New Type structure with variables for SurfaceLeaving and BRDFs
!    -- "MieTmat" structure renamed "MieTmatUser"
!    -- "MieTmatUser" contains new "User-defined" aerosol inputs

!  Aerosol I/O routines--
!    -- Subroutines from CODEC, written by M. Christi, adapted for GEMSTOOL
!        ** Read_aerosol_oprop_file            - reads  user-defined aerosols
!        ** Write_aerosol_oprop_file           - writes user-defined aerosols
!        ** Read_aerosol_oprop_file_plus       - read  for linearized user aerosols
!        ** Read_aerosol_oprop_file_plus_short - read  for linearized user aerosols (just headers)
!        ** Write_aerosol_oprop_file_plus      - write for linearized user aerosols
!        ** Get_fileUnit (Auxiliary routine)

use GEMSTOOL_Input_Types_m

private :: GEMSTOOL_Initialize_ConfigInputs
public  :: GEMSTOOL_Read_Configfiles,          &
           Get_fileUnit,                       &
           Read_aerosol_oprop_file,            &
           Read_aerosol_oprop_file_plus,       &
           Read_aerosol_oprop_file_plus_short, &
           Write_aerosol_oprop_file,           &
           Write_aerosol_oprop_file_plus

contains

subroutine GEMSTOOL_Read_Configfiles &
    ( do_wavnums, ConfigPath, GEMSTOOL_INPUTS, fail, message   )

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Control for using wavenumbers

   logical, intent(in)        :: do_wavnums

!  path for configuration files

   CHARACTER*(*), intent(in)  :: ConfigPath

! GEOMTOOL inputs, structure. Intent(out) here

   TYPE(GEMSTOOL_Config_Inputs), intent(inout) :: GEMSTOOL_INPUTS

!  Excpetion handling
   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables

   logical       :: no_H2O, no_CH4
   logical       :: do_Instrument, trawl
   integer       :: g, k, kdum, nc, nd, nbands, ifail
   real(fpk)     :: banddata(7)
   character*7   :: bandid
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

!  Add User defined aerosol flag, 10/20/16

   filnam = filnam_save
   filnam = adjustl(trim(ConfigPath))//'GEMSTOOL_Atmosphere.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_Tracegases         ! turn on, or off do trace gases
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_aerosols           ! turn on, or off do aerosols calculation
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols       ! turn on, or off do aerosols Mie calculation
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols      ! turn on, or off do aerosols T-matrix calculation
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_User_aerosols      ! turn on, or off User-defined aerosols
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_clouds             ! turn on, or off General Cloud flag
!   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_fileread_wavenums             ! turn on special spectrum
   CLOSE(1)

!  Aerosol flags, Checking expanded to include User option. 10/20/16

   if ( GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols.and.GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols) then
      message = 'Cannot have Mie and Tmatrix aerosols, error in GEMSTOOL_Atmosphere.cfg'
      fail = .true. ;  return
   endif
   if ( GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols.and.GEMSTOOL_INPUTS%Atmosph%do_User_aerosols) then
      message = 'Cannot have Mie and User-defined aerosols, error in GEMSTOOL_Atmosphere.cfg'
      fail = .true. ;  return
   endif
   if ( GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols.and.GEMSTOOL_INPUTS%Atmosph%do_User_aerosols) then
      message = 'Cannot have Tmat and User-defined aerosols, error in GEMSTOOL_Atmosphere.cfg'
      fail = .true. ;  return
   endif

!  Error in previous code, expansion to include user-defined aerosols

   if ( GEMSTOOL_INPUTS%Atmosph%do_aerosols.and. &
              ( .not.GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols  .and. &
                .not.GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols .and. &
                .not.GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) ) then
      message = 'One of Mie/Tmatrix/User aerosol must be set when general aerosol flag set '//&
                'in GEMSTOOL_Atmosphere.cfg'
      fail = .true. ;  return
   endif

!  2a/2b. Wavelengths or wavenumbers
!  ---------------------------------

!  Monochromatic flag introduced 10/26/16.
!    - Only for wavelengths (UVN); if not set --> will do an Instrument simulation

   filnam = filnam_save
   if ( do_wavnums ) then
!mick fix 1/24/2017 - initialize "do_Monochromatic" for this case
       GEMSTOOL_INPUTS%Lambdas%do_Monochromatic = .false.
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
      read(1,*)GEMSTOOL_INPUTS%Lambdas%do_Monochromatic        ! Flag for standard calculation 10/26/16
      read(1,*)GEMSTOOL_INPUTS%Lambdas%lambda_start            ! Initial wavelength
      read(1,*)GEMSTOOL_INPUTS%Lambdas%lambda_finish           ! Final wavelength
      read(1,*)GEMSTOOL_INPUTS%Lambdas%lambda_wavres           ! Resolution
      CLOSE(1)
   endif

!  2c. Satellite simulation
!  ------------------------

!  Introduced 10/26/16. Only for UVN (wavelengths) and non-monochromatic
!    choices are HIMAWARI-8 (1), GOCI-1 (2) etc....

   do_Instrument = .false.
   filnam = filnam_save
   if ( .not.do_wavnums .and. .not.GEMSTOOL_INPUTS%Lambdas%do_Monochromatic ) then
      do_Instrument = .true.
      filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_Instrument.cfg'
      OPEN(1,file=trim(filnam),err=90,status='old')
      read(1,*)GEMSTOOL_INPUTS%Instruments%Instrument_choice     ! HIMAWARI-8 etc...
      read(1,*)GEMSTOOL_INPUTS%Instruments%Instrument_name
      read(1,*)GEMSTOOL_INPUTS%Instruments%Instrument_band
      read(1,*)GEMSTOOL_INPUTS%Instruments%do_Convolution
      read(1,*)GEMSTOOL_INPUTS%Instruments%do_Preconvolution
   endif

!  Check inputs 10/26/16
!    Currently only HIMIWARI-8, bands 1-2 (visible)

   if ( do_Instrument ) then
      if ( GEMSTOOL_INPUTS%Instruments%Instrument_choice .ne. 1 ) then
         message = 'Instrument  is not HIMAWARI-8, currently the only allowed choice (10/26/16)!'
         fail = .true. ;  return
      else 
         if ( GEMSTOOL_INPUTS%Instruments%Instrument_band .gt. 2 ) then
            message = 'HIMAWARI-8 band number must be 1 or 2 for Visible range(10/26/16)!'
            fail = .true. ;  return
        endif
      endif
   endif

!  Now read the band data
!  --> HIMAWARI is always in wavenumber, will be converted later. Use 1% RSR threshold.

   if ( do_Instrument ) then
     if ( GEMSTOOL_INPUTS%Instruments%Instrument_choice .eq. 1 ) then
       filnam = 'InstrumentFiles/HIMAWARI-8'
       open(1,file=Trim(filnam), iostat = IFAIL, status='old')
       if (IFAIL .ne. 0 ) then
         message = ' --> '//Trim(filnam)//' <-- NOT FOUND, abandon'
         close(1) ; fail = .true. ; return
       else  
         read(1,*) ; read(1,*)nbands ; read(1,*) ; trawl = .true. ; nc = 0
         do while (trawl .and.nc.lt.nbands)
           nc = nc + 1 ; read(1,88)nd, bandid, banddata(1:7)
           if ( nd.eq.GEMSTOOL_INPUTS%Instruments%Instrument_band ) trawl = .false. 
         enddo
         GEMSTOOL_INPUTS%Instruments%band_start     = banddata(2)
         GEMSTOOL_INPUTS%Instruments%band_finish    = banddata(1)
         GEMSTOOL_INPUTS%Instruments%band_median    = 0.5d0 * ( banddata(2) + banddata(1) )
         GEMSTOOL_INPUTS%Instruments%band_fineres   = banddata(5)
         GEMSTOOL_INPUTS%Instruments%band_extension = banddata(6)
         close(1)
         GEMSTOOL_INPUTS%Instruments%band_slitpars(1) = banddata(5)   ! Gaussian slit HWHM = Resolution
         GEMSTOOL_INPUTS%Instruments%band_slitpars(2) = 0.0_fpk
       endif
     else if ( GEMSTOOL_INPUTS%Instruments%Instrument_choice .eq. 2 ) then
       filnam = 'InstrumentFiles/GOCI-1' ; ! PLACEHOLDER
     else if ( GEMSTOOL_INPUTS%Instruments%Instrument_choice .eq. 3 ) then
       filnam = 'InstrumentFiles/GOCI-2' ; ! PLACEHOLDER
     endif
   endif      
88 FORMAT(I2,1x,A6,7F9.1)

!  3. Radiative transfer control
!  -----------------------------

   filnam = filnam_save
   filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_RTControl.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')
   read(1,*)GEMSTOOL_INPUTS%RTControl%NVlidort_nstreams     ! Number of Discrete Ordinates in VLIDORT
   read(1,*)GEMSTOOL_INPUTS%RTControl%NVlidort_nstokes      ! number of stokes parameters in VLIDORT  
   read(1,*)GEMSTOOL_INPUTS%RTControl%do_Upwelling          ! turn on flag for Upwelling
   read(1,*)GEMSTOOL_INPUTS%RTControl%observation_height    ! New part; mchoi
   read(1,*)GEMSTOOL_INPUTS%RTControl%do_PathRadiance       ! turn on flag for "Path Radiance"      (TBD)
   read(1,*)GEMSTOOL_INPUTS%RTControl%do_SphericalAlbedo    ! turn on flag for "Spherical Albedo"   (TBD)
   read(1,*)GEMSTOOL_INPUTS%RTControl%do_2wayTransmittance  ! turn on flag for "2-wayTransmittance" (TBD)
   read(1,*)GEMSTOOL_INPUTS%RTControl%do_firstorder_option  ! turn on, using new FO code
   read(1,*)GEMSTOOL_INPUTS%RTControl%FO_do_regular_ps      ! turn on, using FO code regular PS mode
   close(1)

!  4. Time/Location
!  ----------------

!  Hour/Minute/second variables added for specification of SIF Epoch
!    10/18/16. Rob Spurr

   filnam = filnam_save
   filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_TimePosition.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')
   READ(1,*)GEMSTOOL_INPUTS%TimePos%latitude
   READ(1,*)GEMSTOOL_INPUTS%TimePos%longitude
   READ(1,*)GEMSTOOL_INPUTS%TimePos%earthradius
   READ(1,*)GEMSTOOL_INPUTS%TimePos%year
   READ(1,*)GEMSTOOL_INPUTS%TimePos%month
   READ(1,*)GEMSTOOL_INPUTS%TimePos%day_of_month
   READ(1,*)GEMSTOOL_INPUTS%TimePos%Hour
   READ(1,*)GEMSTOOL_INPUTS%TimePos%Minute
   READ(1,*)GEMSTOOL_INPUTS%TimePos%Second
   CLOSE(1)

!  5. Trace gas information
!  ------------------------

!  @@ Rob fix 7/23/14, add H2O scaling flag and value
!  @@ Y.Jung fix 2/1/15, add CH4 scaling flag and value

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
      read(1,'(L2,1X,f10.5)')GEMSTOOL_INPUTS%Tracegas%do_H2OScaling, &
                         GEMSTOOL_INPUTS%Tracegas%H2OScaling
      read(1,'(L2,1X,f10.5)')GEMSTOOL_INPUTS%Tracegas%do_CH4Scaling, &
                         GEMSTOOL_INPUTS%Tracegas%CH4Scaling

      CLOSE(1)
   endif

!  @@ Rob fix 7/23/14, Exception handling for the H2O scaling

   if ( GEMSTOOL_INPUTS%Tracegas%do_H2OScaling ) then
      no_H2O = .true.
      do g = 1, GEMSTOOL_INPUTS%Tracegas%ngases
         if (GEMSTOOL_INPUTS%Tracegas%which_gases(g).eq.'H2O ' ) no_H2O = .false. !! CHECK
      enddo
      if (no_H2O) then
         message = 'Cannot have H2O scaling if H2O is not included in the gas list'
         fail = .true. ;  return
      endif
   endif

!  @@ Y.Jung fix 2/1/15,  Exception handling for the CH4 scaling

   if ( GEMSTOOL_INPUTS%Tracegas%do_CH4Scaling ) then
      no_CH4 = .true.
      do g = 1, GEMSTOOL_INPUTS%Tracegas%ngases
         if (GEMSTOOL_INPUTS%Tracegas%which_gases(g).eq.'CH4 ' ) no_CH4 = .false. !! CHECK
      enddo
      if (no_CH4) then
         message = 'Cannot have CH4 scaling if CH4 is not included in the gas list'
         fail = .true. ;  return
      endif
   endif

!  6. Aerosol Loading
!  ------------------

!  Exapnded to include User aerosol too.

   if (  GEMSTOOL_INPUTS%Atmosph%do_aerosols ) then
      if ( GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols  .or. &
           GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols .or. &
           GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then

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

      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%do_bimodal                                                 ! Flag for Bimodal aerosol
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%bimodal_fraction                                           ! Bimodal fraction

!  PSD parameters

      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%PSDindex(1),  GEMSTOOL_INPUTS%MieTmatUser%PSDindex(2)      ! Mie Particle Size Distribution (PSD) indices
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1,1), GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1,2)     ! Mie PSD First  parameter values
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%PSDpars(2,1), GEMSTOOL_INPUTS%MieTmatUser%PSDpars(2,2)     ! Mie PSD Second parameter values
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%PSDpars(3,1), GEMSTOOL_INPUTS%MieTmatUser%PSDpars(3,2)     ! Mie PSD Third  parameter values

!  PSD Integration control

      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%nblocks(1),    GEMSTOOL_INPUTS%MieTmatUser%nblocks(2)      ! Mie Number of blocks covering the PSDs
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%nweights(1),   GEMSTOOL_INPUTS%MieTmatUser%nweights(2)     ! Mie Number of Quadrature weights in each block
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%xparticle_limit                                            ! Mie Particle size limit
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%FixR1R2(1),    GEMSTOOL_INPUTS%MieTmatUser%FixR1R2(2)      ! Mie Use particle radius cutoff to fix Rmin/Rmax
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%R1R2_cutoff(1),GEMSTOOL_INPUTS%MieTmatUser%R1R2_cutoff(2)  ! Mie Particle radius cutoff 
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%R1(1), GEMSTOOL_INPUTS%MieTmatUser%R1(2)                   ! Mie Minimum particle radii (Microns)
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%R2(1), GEMSTOOL_INPUTS%MieTmatUser%R2(2)                   ! Mie/Maximum particle radii (Microns)

!  refractive indices

      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%nreal(1), GEMSTOOL_INPUTS%MieTmatUser%nreal(2)             ! Mie Real parts of refractive index
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%nimag(1), GEMSTOOL_INPUTS%MieTmatUser%nimag(2)             ! Mie Imaginary part of refractive index

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

      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%do_bimodal                                                 ! Flag for Bimodal aerosol
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%bimodal_fraction                                           ! Bimodal fraction

!  Orientation control

      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%Do_EqSaSphere                                              ! Tmat flag for equivalent-surface-sphere representation
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%Tmat_Sphtype                                               ! Tmat Type of Spheroidal particle
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%Tmat_ndgs(1),GEMSTOOL_INPUTS%MieTmatUser%Tmat_ndgs(2)      ! Tmat Gaussian surface numbers
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%Tmat_accuracy                                              ! Tmat Accuracy for the computation

!  PSD parameters

      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%PSDindex(1),  GEMSTOOL_INPUTS%MieTmatUser%PSDindex(2)      ! Tmat Particle Size Distribution (PSD) indices
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1,1), GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1,2)     ! Tmat PSD First  parameter values
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%PSDpars(2,1), GEMSTOOL_INPUTS%MieTmatUser%PSDpars(2,2)     ! Tmat PSD Second parameter values
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%PSDpars(3,1), GEMSTOOL_INPUTS%MieTmatUser%PSDpars(3,2)     ! Tmat PSD Third  parameter values

!  PSD Integration control

      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%Tmat_nkmax(1), GEMSTOOL_INPUTS%MieTmatUser%Tmat_nkmax(2)   ! Tmat Gaussian PSD numbers
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%FixR1R2(1),    GEMSTOOL_INPUTS%MieTmatUser%FixR1R2(2)      ! Tmat Use particle radius cutoff to fix Rmin/Rmax
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%R1R2_cutoff(1),GEMSTOOL_INPUTS%MieTmatUser%R1R2_cutoff(2)  ! Tmat Particle radius cutoff 
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%R1(1), GEMSTOOL_INPUTS%MieTmatUser%R1(2)                   ! Tmat Minimum particle radii (Microns)
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%R2(1), GEMSTOOL_INPUTS%MieTmatUser%R2(2)                   ! Tmat/Maximum particle radii (Microns)

!  refractive indices, shape factor

      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%nreal(1), GEMSTOOL_INPUTS%MieTmatUser%nreal(2)             ! Tmat Real parts of refractive index
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%nimag(1), GEMSTOOL_INPUTS%MieTmatUser%nimag(2)             ! Tmat Imaginary part of refractive index
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%Tmat_eps(1),GEMSTOOL_INPUTS%MieTmatUser%Tmat_eps(2)        ! Tmat shape factors

! End T-matrix aerosols

      close(1)
   endif

!  9. User-defined aerosol control
!  -------------------------------

   if ( GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then
      filnam = filnam_save
      filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_UserAerosols.cfg'
      OPEN(1,file=trim(filnam),err=90,status='old')

 !  Bimodal control

      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%do_bimodal                          ! Flag for Bimodal aerosol
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%bimodal_fraction                    ! Bimodal fraction

!  file names
      
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%User_Aerofile_names(1)
      read(1,*)GEMSTOOL_INPUTS%MieTmatUser%User_Aerofile_names(2)

! close file

      close(1)
   endif

!  10. Cloud control
!  -----------------

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

!  11. Albedo Closure control
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

!  12. BRDF Surface Control
!  ------------------------

!  Added, 10/24/16. Rob Spurr
!    1. General flag for using BRDF surface instead of Lambertian
!    2. Kernel model choices, as with VLIDORT's BRDF supplement
!    3. Kernel information not used, if no BRDF --> use Lambertian surface

   filnam = filnam_save
   filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_BRDF.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')
   read(1,*)GEMSTOOL_INPUTS%BRDF%do_Gemstool_BRDF
   if  ( GEMSTOOL_INPUTS%BRDF%do_Gemstool_BRDF ) then
      read(1,*)GEMSTOOL_INPUTS%BRDF%n_Kernels
      do g = 1, GEMSTOOL_INPUTS%BRDF%n_Kernels
         read(1,*)GEMSTOOL_INPUTS%BRDF%Kernel_indices(g),     &
                  GEMSTOOL_INPUTS%BRDF%Kernel_n_parameters(g),&
                  GEMSTOOL_INPUTS%BRDF%Kernel_factors(g),     &
                  GEMSTOOL_INPUTS%BRDF%Kernel_parameters(g,1:3)
      enddo
   endif
   close(1)

   ! write(*,*) 'Read_Configfiles'
   ! write(*,*) GEMSTOOL_INPUTS%BRDF%Kernel_indices
   ! write(*,*) GEMSTOOL_INPUTS%BRDF%Kernel_n_parameters
   ! write(*,*) GEMSTOOL_INPUTS%BRDF%Kernel_factors
   ! write(*,*) GEMSTOOL_INPUTS%BRDF%Kernel_parameters
   ! pause

!  13. Surface-Leaving control
!  ---------------------------

!  SIF Added, 10/18/16, updated 11/30/16 for linear parameterization. Rob Spurr
!    1. General flag for including SIF
!    2. Amplitude-Scaling of database 755 nm Fluorescence
!     - SIF will be assumed Isotropic
!     - EITHER Gaussian shapes from the database, OR Linearily parameterized
!     - GAUSSIAN Scaling-amplitude  default value should be 1.0
!     - LINEAR   Scaling Constant,  default value should be 1.0
!     - LINEAR   Scaling Gradient,  default value should be 27.0 (fit-to-Gaussians, 755-775 nm)

   filnam = filnam_save
   filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_SurfaceLeaving.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')

!  top level control : SIF or Waterleaving. Cannot have both

   read(1,*)  ! Header
   read(1,*)GEMSTOOL_INPUTS%Sleave%DO_WaterLeaving
   read(1,*)GEMSTOOL_INPUTS%Sleave%DO_OceanColor
   read(1,*)GEMSTOOL_INPUTS%Sleave%DO_SIF

   if ( GEMSTOOL_INPUTS%Sleave%DO_WaterLeaving .and.GEMSTOOL_INPUTS%Sleave%DO_SIF ) then
      message = 'Cannot have goth Fluorescence and Water Leavgin- Turn off one of them'
      fail = .true. ; return
   endif

!  Water-leaving - read header and next 7 lines then close
!  SIF           - Skip 8 lines, read header and inputs, then close
!                - (11/30/16). Extra control for linear parameterization, Exact-Only solution

   IF ( GEMSTOOL_INPUTS%Sleave%DO_WaterLeaving ) THEN
      read(1,*)  ! Header
      read(1,*)GEMSTOOL_INPUTS%Sleave%Chlorconc          ! units [mg/M]
      read(1,*)GEMSTOOL_INPUTS%Sleave%Salinity           ! units [ppt]
      read(1,*)GEMSTOOL_INPUTS%Sleave%WindSpeed          ! units Metres/second
      read(1,*)GEMSTOOL_INPUTS%Sleave%WindDir            ! azimuth relative to solar
      read(1,*)GEMSTOOL_INPUTS%Sleave%DO_GlintShadow     ! recommended T
      read(1,*)GEMSTOOL_INPUTS%Sleave%DO_FoamOption      ! recommended T
      read(1,*)GEMSTOOL_INPUTS%Sleave%DO_FacetIsotropy   ! recommended T, if set ignores WindDir
   ELSE
      read(1,*) ; read(1,*) ; read(1,*) ; read(1,*) ; read(1,*) ; read(1,*) ; read(1,*) ; read(1,*)
      read(1,*)  ! Header
      read(1,*)GEMSTOOL_INPUTS%Sleave%DO_SIF_ExactOnly
      read(1,*)GEMSTOOL_INPUTS%Sleave%DO_SIF_LinearParam
      read(1,*)GEMSTOOL_INPUTS%Sleave%SIF755_Scaling_Constant
      read(1,*)GEMSTOOL_INPUTS%Sleave%SIF755_Scaling_Gradient
      read(1,*)GEMSTOOL_INPUTS%Sleave%SIF755_Amplitude
   endif

   close(1)

!  14. Geometrical inputs
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

!  15. Linearization Control
!  -------------------------

!  NOTE - these quantities always read-in, but not used in Radiance-only wrapper.

!  Introduced, 24 October 2013
!  @@ Rob   fix 7/23/14 , add H2O scaling Jacobian flag
!  @@ YJung fix 2/1/15  , add CH4 scaling Jacobian flag
!  @@ Rob   fix 10/18/16, add SIF Jacobian flag, updated 11/30/16

   filnam = filnam_save
   filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_LinControl.cfg'
   OPEN(1,file=trim(filnam),err=90,status='old')
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_GasProfile_Jacobians          ! 1
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_AerOpdepProfile_Jacobians     ! 2
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_AerBulk_Jacobians             ! New, June 2014
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_Surface_Jacobians             ! 4
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_Tshift_Jacobian               ! 5
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_SurfPress_Jacobian            ! 6
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_H2OScaling_Jacobian           ! New July 2014
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_CH4Scaling_Jacobian           ! New Feb 2015
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_SIF_Jacobians                 ! New Oct/Nov 2016
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_normalized_wfoutput           ! 9
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_hitran                        ! New Feb 2015

   print*,'do_HITRAN=', GEMSTOOL_INPUTS%LinControl%do_hitran
        
   CLOSE(1)

!  @@ Rob fix 7/23/14, check compatibility of H2O scaling flags

   if ( GEMSTOOL_INPUTS%LinControl%do_H2OScaling_Jacobian .and..not.GEMSTOOL_INPUTS%Tracegas%do_H2OScaling ) then
      message = 'H2O scaling: Jacobian flag turned on, but scaling flag is false'
      fail = .true. ; return
   endif

!  @@ Y.Jung fix 2/1/15, check compatibility of CH4 scaling flags

   if ( GEMSTOOL_INPUTS%LinControl%do_CH4Scaling_Jacobian .and..not.GEMSTOOL_INPUTS%Tracegas%do_CH4Scaling ) then
      message = 'CH4 scaling: Jacobian flag turned on, but scaling flag is false'
      fail = .true. ; return
   endif

!  New file, Introduced June 2014 for the Aerosol Bulk properties

   if ( GEMSTOOL_INPUTS%LinControl%do_AerBulk_Jacobians ) then
      filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_LinControl_AerBulk.cfg'
      OPEN(1,file=trim(filnam),err=90,status='old')
      read(1,*)GEMSTOOL_INPUTS%LinControl%do_AerBulk_LoadPars_Jacobians
      read(1,*)GEMSTOOL_INPUTS%LinControl%do_AerBulk_RefIndex_Jacobians
      read(1,*)GEMSTOOL_INPUTS%LinControl%do_AerBulk_ShapeFac_Jacobians
      read(1,*)GEMSTOOL_INPUTS%LinControl%do_AerBulk_SizeDist_Jacobians
      read(1,*)GEMSTOOL_INPUTS%LinControl%do_AerBulk_BmodFrac_Jacobian
      CLOSE(1)
   endif

!  New file, Introduced by Myungje Choi for Tshift input read

   if ( GEMSTOOL_INPUTS%LinControl%do_Tshift_Jacobian ) then
      filnam=adjustl(trim(ConfigPath))//'GEMSTOOL_Tshift.cfg'
      OPEN(1,file=trim(filnam),err=90,status='old')
      read(1,*)GEMSTOOL_INPUTS%LinControl%Tshift_read
      CLOSE(1)
   endif



!  Exception handling
!  @@ Rob fix 7/23/14, Updated to include H2O scaling Jacobian flag
!  @@ Y.Jung fix 2/1/15, Updated to include CH4 scaling Jacobian flag

   if ( GEMSTOOL_INPUTS%LinControl%do_GasProfile_Jacobians.or.GEMSTOOL_INPUTS%LinControl%do_AerOpdepProfile_Jacobians ) then
      if (  GEMSTOOL_INPUTS%LinControl%do_Tshift_Jacobian     .or. &
            GEMSTOOL_INPUTS%LinControl%do_SurfPress_Jacobian  .or. &
            GEMSTOOL_INPUTS%LinControl%do_H2OScaling_Jacobian .or. &
            GEMSTOOL_INPUTS%LinControl%do_CH4Scaling_Jacobian .or. &
            GEMSTOOL_INPUTS%LinControl%do_AerBulk_Jacobians )  then
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

!  routine has been upgraded to include new type structures
!  Also, everything is now properly initialized

   implicit none

! GEOMTOOL inputs, structure. Intent(out) here

   TYPE(GEMSTOOL_Config_Inputs), INTENT(INOUT) :: GEMSTOOL_INPUTS

!  1. Atmospheric constituents
!  ---------------------------

   GEMSTOOL_INPUTS%Atmosph%do_Tracegases    = .false.
   GEMSTOOL_INPUTS%Atmosph%do_aerosols      = .false.
   GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols  = .false.
   GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols = .false.
   GEMSTOOL_INPUTS%Atmosph%do_User_aerosols = .false.
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
   GEMSTOOL_INPUTS%RTControl%observation_height    = 0.0d0
   GEMSTOOL_INPUTS%RTControl%do_PathRadiance       = .false.
   GEMSTOOL_INPUTS%RTControl%do_SphericalAlbedo    = .false.
   GEMSTOOL_INPUTS%RTControl%do_2wayTransmittance  = .false.
   GEMSTOOL_INPUTS%RTControl%do_firstorder_option  = .false.
   GEMSTOOL_INPUTS%RTControl%FO_do_regular_ps      = .false. 

!  4. Time/Location
!  ----------------

!  Hour/Minute/Second  variables added for  specification of SIF Epoch
!    10/18/16. Rob Spurr

   GEMSTOOL_INPUTS%TimePos%latitude     = 0.0d0
   GEMSTOOL_INPUTS%TimePos%longitude    = 0.0d0
   GEMSTOOL_INPUTS%TimePos%earthradius  = 0.0d0
   GEMSTOOL_INPUTS%TimePos%year         = 0
   GEMSTOOL_INPUTS%TimePos%month        = 0
   GEMSTOOL_INPUTS%TimePos%day_of_month = 0
   GEMSTOOL_INPUTS%TimePos%Hour         = 0
   GEMSTOOL_INPUTS%TimePos%Minute       = 0
   GEMSTOOL_INPUTS%TimePos%second       = 0

!  5. Trace gas information
!  ------------------------

!  @@ Rob fix 7/23/14, Updated to include H2O scaling flag and value
!  @@ Y.Jung fix 2/1/15, Updated to include CH4 scaling flag and value

   GEMSTOOL_INPUTS%Tracegas%ngases      = 0
   GEMSTOOL_INPUTS%Tracegas%which_gases = '    '
   GEMSTOOL_INPUTS%Tracegas%do_gases    = .false.
   GEMSTOOL_INPUTS%Tracegas%do_gas_wfs  = .false.
   GEMSTOOL_INPUTS%Tracegas%Gas_Profile_names = ' '

   GEMSTOOL_INPUTS%Tracegas%do_H2OScaling = .false.
   GEMSTOOL_INPUTS%Tracegas%H2OScaling    = 0.0d0

   GEMSTOOL_INPUTS%Tracegas%do_CH4Scaling = .false.
   GEMSTOOL_INPUTS%Tracegas%CH4Scaling    = 0.0d0

!  12. SIF inputs
!  --------------

! @ Rob added 10/18/16
! @ Rob added 11/30/16. Linear parameterization control, Exact-only control

   GEMSTOOL_INPUTS%Sleave%DO_SIF = .false.
   GEMSTOOL_INPUTS%Sleave%DO_SIF_ExactOnly   = .false.
   GEMSTOOL_INPUTS%Sleave%DO_SIF_LinearParam = .false.
   GEMSTOOL_INPUTS%Sleave%SIF755_Scaling_Constant = 0.0d0
   GEMSTOOL_INPUTS%Sleave%SIF755_Scaling_Gradient = 0.0d0
   GEMSTOOL_INPUTS%Sleave%SIF755_AMPLITUDE        = 0.0d0

!  AND SO ON

   return
end subroutine GEMSTOOL_Initialize_ConfigInputs

!

subroutine Get_FileUnit ( UnitNum )

   implicit none

!  Arguments
!  =========

   integer, intent(out) :: UnitNum

!  Local
!  =====

   logical :: InUse

!  Obtain a file unit number not already in use

   UnitNum=1
   inquire (unit=UnitNum,opened=InUse)
   do while (InUse)
      UnitNum=UnitNum+1
      inquire(unit=UnitNum,opened=InUse)
   enddo

!  Finish

   return
end subroutine Get_FileUnit

!

subroutine Read_aerosol_oprop_file ( &
   max_user_angles, InFileUnit, Micron_wavelength_in, RefWvl, & ! Inputs
   maxAerModes, nAerMode, nPLines, DoPreamble,                & ! I/O
   MTU_dist, MTU_bulk, MTU_asymm, MTU_ncoeffs, MTU_expcoeffs, & ! Outputs
   fail, MTU_messages )                                         ! Exception handling

!  Introduced 10/19/16. Adapted from CODECTOOL routine of the same name.
!    Major new thing: aerosol distribution characteristics for Ref Wavl.

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   integer, intent(in)    :: max_user_angles
   integer, intent(in)    :: InFileUnit
   real(fpk), intent(in)  :: Micron_wavelength_in
   logical, intent(in)    :: RefWvl

!mick mod 1/24/2017 - added aerosol file preamble capacity
   integer, intent(in)           :: maxAerModes, nAerMode, nPLines
   logical, intent(inout)        :: DoPreamble(maxAerModes)

!  Distribution parameters

   real(fpk), intent(out) :: MTU_dist(5)

!  Bulk optical  parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk), intent(out) :: MTU_bulk(3)

!  Asymmetry parameter and number of expansion coefficients

   real(fpk), intent(out) :: MTU_asymm

!  Number and value of expansion coefficients

   integer, intent(out)    :: MTU_ncoeffs
   real(fpk), intent(out)  :: MTU_expcoeffs(6,0:max_user_angles)

!  Error handling

   logical, intent(out)      :: fail
   character*90, intent(out) :: MTU_messages(3)

!  Local variables
!  ---------------

   integer    :: i,l,ldud,w
   real(fpk)  :: DiffThreshold,Micron_wavelength
   character(len=4) :: RefID

   fail         = .false.
   MTU_messages = ' '

!  Read input from file

!mick mod 1/24/2017 - added aerosol file preamble capacity
   if ( DoPreamble(nAerMode) ) then
      do i=1,nPLines
         read(InFileUnit,*) !Preamble(i)
      enddo
      DoPreamble(nAerMode) = .false.
   endif

   read(InFileUnit,*) ! Divider: '--------------------------------------------'

   read(InFileUnit,*) ! Header : 'Wvl idx','Wvl (um)'
   if ( .not.RefWvl ) then
      read(InFileUnit,*) w,Micron_wavelength
   else
      read(InFileUnit,*) RefID,Micron_wavelength
   endif

   DiffThreshold = 1.0d-10
   if ( dabs(Micron_wavelength - Micron_wavelength_in) > DiffThreshold ) then
!mick mod 1/24/2017 - modified error handling
      !!Read unexpected wavelength --> abort!
      !MTU_messages(1) = 'Failure: Aerosol optical property input file, short version'
      !!MTU_messages(2) --> file name - defined outside aubroutine
      !MTU_messages(3) = 'Problem: Read unexpected wavelength'
      write(*,*)
      write(*,'(a)')      'Error: mismatch in wavelengths'
      write(*,'(a,f6.1)') 'Wavelength input to subroutine  (in um) = ',Micron_wavelength_in
      write(*,'(a,f6.1)') 'Wavelength input from aero file (in um) = ',Micron_wavelength
      !write(*,*) 'Difference Threshold                       = ',DiffThreshold
      !write(*,*) 'Magnitude of actual difference             = ',dabs(Micron_wavelength - Micron_wavelength_in)

      !Read unexpected wavelength --> abort!
      MTU_messages(1) = 'Failure: Aero optical property input file'
      !MTU_messages(2) --> file name - defined outside subroutine
      MTU_messages(3) = 'Problem: Mismatch between spectral cfg file & aero interpolate flag with aero file wvls'
      fail = .true. ; return
   end if

   if ( .not.RefWvl ) then
      read(InFileUnit,*) ! Header : 'Ext Coeff','Scat Coeff','SSA','Asym Par','Num of Expan Coeffs'
      read(InFileUnit,*) &
         MTU_bulk(1),&  ! Extinction coefficient
         MTU_bulk(2),&  ! Scattering coefficient
         MTU_bulk(3),&  ! Single scattering albedo
         MTU_asymm  ,&  ! Asymmetry parameter
         MTU_ncoeffs    ! Number of Expansion coefficients
      if ( MTU_ncoeffs > max_user_angles ) then
         !# of coefficients exceeds max allowed --> abort!
         MTU_messages(1) = 'Failure: Aerosol optical property input file'
         !MTU_messages(2) --> file name - defined outside aubroutine
         MTU_messages(3) = 'Problem: # of coefficients exceeds max allowed for current CODEC settings'
         fail = .true. ; return
      end if
      read(InFileUnit,*) ! Header : 'Moment','Expan Coeffs (a1,a2,a3,a4,b1,b2)'
      do l=0,MTU_ncoeffs
         read(InFileUnit,*) &
            ldud,(MTU_expcoeffs(i,l),i=1,6) !Index, Expansion coefficient input (a1,a2,a3,a4,b1,b2)
      enddo
   else
      read(InFileUnit,*) ! Header : 'Ext Coeff','Scat Coeff','SSA'
      read(InFileUnit,*) &
         MTU_bulk(1),&  ! Extinction coefficient
         MTU_bulk(2),&  ! Scattering coefficient
         MTU_bulk(3)    ! Single scattering albedo
      read(InFileUnit,*) ! Header : 'Particle Size Distribution Characteristics'
      read(InFileUnit,*) ! Header : 'Normalization','Geom Xsection','Geom. Volume', 'REFF', 'VEFF'
      read(InFileUnit,*) MTU_dist(1:5)

   endif

end subroutine Read_aerosol_oprop_file

!

subroutine Write_aerosol_oprop_file ( &
   max_user_angles, OutFileUnit, w, Micron_wavelength, RefWvl,  & ! Inputs
   maxAerModes, nAerMode, nPLines, DoPreamble, Preamble,        & ! I/O
   MTU_dist, MTU_bulk, MTU_asymm, MTU_ncoeffs, MTU_expcoeffs )    ! Inputs

!  Introduced 10/19/16. Adapted from CODECTOOL routine of the same name.
!    Major new thing: aerosol distribution characteristics for Ref Wavl.

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   integer, intent(in)           :: max_user_angles
   integer, intent(in)           :: OutFileUnit, w
   real(fpk), intent(in)         :: Micron_wavelength
   logical, intent(in)           :: RefWvl

!mick mod 1/24/2017 - added aerosol file preamble capacity
!  Rob Fix 3/7/17. Preamble is now mode-dependent
   integer, intent(in)           :: maxAerModes, nAerMode, nPLines
   logical, intent(inout)        :: DoPreamble(maxAerModes)
   character(len=*), intent(in)  :: Preamble(nPLines,maxAerModes)

!  Distribution parameters

   real(fpk), intent(in)  :: MTU_dist(5)

!  Bulk optical parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk), intent(in)  :: MTU_bulk(3)

!  Asymmetry parameter and number of expansion coefficients

   real(fpk), intent(in)  :: MTU_asymm

!  Number and value of expansion coefficients

   integer, intent(in)    :: MTU_ncoeffs
   real(fpk), intent(in)  :: MTU_expcoeffs(6,0:max_user_angles)

!  Local variables
!  ---------------

   integer :: i,l, Preamble_unit
   character*2 :: c2

!  Write output to file

!mick mod 1/24/2017 - added aerosol file preamble capacity
!Rob mod  3/3/2017  -  always print out preamble file to UserAeroFiles/out/
!Rob mod  3/7/2017  -  Preamble lines are mode-dependent

   if ( DoPreamble(nAerMode) ) then
      Preamble_unit = 1111 ; c2 = '00' ; write(c2(2:2),'(I1)')nAerMode
      open(Preamble_unit,file='UserAeroFiles/out/AerFile_Preamble_Mode_'//c2//'.dat',status='replace')
      do i=1,nPLines
         write(Preamble_unit,'(a)') trim(Preamble(i,nAerMode))
         write(OutFileUnit,'(a)') trim(Preamble(i,nAerMode))
      enddo
      DoPreamble(nAerMode) = .false.
      close(Preamble_unit)
   endif

   write(OutFileUnit,'(3x,a)') '--------------------------------------------' // &
                               '--------------------------------------------' // &
                               '--------------------------------------------'
   write(OutFileUnit,'(6x,a,7x,a)') 'Wvl idx','Wvl (um)'

   if ( .not.RefWvl ) then

!  Header
!  ------

      write(OutFileUnit,'(8x,i4.4,e20.11)') w,Micron_wavelength

!  Bulk properties
!  ---------------

      !Sample:
      !      Ext Coeff           Scat Coeff              SSA               Asym Par       Num of Expan Coeffs
      !   0.30186852870E-01   0.29730300500E-01   0.98487578775E+00   0.59090185717E+00           0
      write(OutFileUnit,'(9x,a,11x,a,14x,a,15x,a,7x,a)') 'Ext Coeff','Scat Coeff','SSA','Asym Par','Num of Expan Coeffs'
      write(OutFileUnit,'(3x,4e20.11,9x,i4)') &  ! AMS 09 JAN 2015 Updated expansion coeffs from i3 to i4 as sometimes we have >1000.
         MTU_bulk(1),& ! Extinction coefficient
         MTU_bulk(2),& ! Scattering coefficient
         MTU_bulk(3),& ! Single scattering albedo
         MTU_asymm  ,& ! Asymmetry parameter
         MTU_ncoeffs   ! Number of Expansion coefficients

!  Coefficients
!  ------------

      !Sample:
      !      Moment                                             Expan Coeffs (a1,a2,a3,a4,b1,b2)
      !       0000   0.10000000000E+01   0.00000000000E+00   0.00000000000E+00   0.85192295256E+00   0.00000000000E+00  -0.00000000000E+00
      write(OutFileUnit,'(6x,a,45x,a)') 'Moment','Expan Coeffs (a1,a2,a3,a4,b1,b2)'
      do l=0,MTU_ncoeffs
         write(OutFileUnit,'(7x,i4.4,6e20.11)') &
            l,(MTU_expcoeffs(i,l),i=1,6) ! Index, Expansion coefficient output (a1,a2,a3,a4,b1,b2)
      enddo

!  REFERENCE WAVELENGTH ONLY
!  =========================

   else

!  Header
!  ------

      write(OutFileUnit,'(8x,a,e20.11)') 'Ref#',Micron_wavelength

!  Bulk properties only
!  --------------------

      !Sample:
      !      Ext Coeff           Scat Coeff              SSA
      !   0.30186852870E-01   0.29730300500E-01   0.98487578775E+00
      write(OutFileUnit,'(9x,a,11x,a,14x,a)') 'Ext Coeff','Scat Coeff','SSA'
      write(OutFileUnit,'(3x,4e20.11,9x,i3)') &
         MTU_bulk(1),& ! Extinction coefficient
         MTU_bulk(2),& ! Scattering coefficient
         MTU_bulk(3)   ! Single scattering albedo

!  Size Distribution parameters (for reference only)
!  ----------------------------

      !Sample:
      !           ----------------Particle Size Distribution Characteristics----------------
      !     Normalization       Geom Xsection        Geom Volume      Effective radius    Effective variance
      !   0.30186852870E-01   0.29730300500E-01   0.98487578775E+00   0.98487578775E+00   0.98487578775E+00
      write(OutFileUnit,'(17x,a)') '----------------Particle Size Distribution Characteristics----------------'
      write(OutFileUnit,'(8x,a,7x,a,8x,a,5x,a,4x,a)') 'Normalization','Geom Xsection','Geom. Volume',&
                                                      'Effective radius','Effective variance'
      write(OutFileUnit,'(3x,5e20.11)') MTU_dist(1:5)

   endif

end subroutine Write_aerosol_oprop_file

!

subroutine Read_aerosol_oprop_file_plus ( &
   max_user_angles, InFileUnit, Micron_wavelength_in, RefWvl,        & ! Inputs
   maxAerModes, nAerMode, nPLines, DoPreamble,                       & ! I/O
   do_LinearRef, do_LinearPSD, NLIN, NPSD, PARSLIST,                 & ! Outputs
   MTU_dist, MTU_bulk, MTU_asymm, MTU_ncoeffs, MTU_expcoeffs,        & ! Outputs
   LRFE_MTU_bulk, LRFE_MTU_asymm, LRFE_MTU_expcoeffs,                & ! Outputs
   LPSD_MTU_bulk, LPSD_MTU_asymm, LPSD_MTU_expcoeffs, LPSD_MTU_dist, & ! Outputs
   fail, MTU_messages )                                                ! Exception handling

!  Introduced 10/19/16. Adapted from CODECTOOL routine of the same name.
!    Major new thing: aerosol distribution characteristics for Ref Wavl.

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   integer, intent(in)    :: max_user_angles
   integer, intent(in)    :: InFileUnit
   real(fpk), intent(in)  :: Micron_wavelength_in
   logical, intent(in)    :: RefWvl

!mick mod 1/24/2017 - added aerosol file preamble capacity
   integer, intent(in)           :: maxAerModes, nAerMode, nPLines
   logical, intent(inout)        :: DoPreamble(maxAerModes)

!  Linearization control

   logical  , intent(out)    :: do_LinearPSD, do_LinearRef
   integer  , intent(out)    :: NPSD, NLIN
   real(fpk), intent(out)    :: PARSLIST(10)

!  Distribution parameters

   real(fpk), intent(out) :: MTU_dist(5)

!  Bulk optical  parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk), intent(out) :: MTU_bulk(3)

!  Asymmetry parameter and number of expansion coefficients

   real(fpk), intent(out) :: MTU_asymm

!  Number and value of expansion coefficients

   integer, intent(out)    :: MTU_ncoeffs
   real(fpk), intent(out)  :: MTU_expcoeffs(6,0:max_user_angles)

!  Linearized quantities

   real(fpk), intent(out)   :: LRFE_MTU_bulk (3,3)     ! 3 in second dimension, could be used by Tmat
   real(fpk), intent(out)   :: LPSD_MTU_bulk (3,3)
  
   real(fpk), intent(out)   :: LRFE_MTU_asymm(3)       ! 3 in second dimension, could be used by Tmat
   real(fpk), intent(out)   :: LPSD_MTU_asymm(3)

   real(fpk), intent(out)   :: LPSD_MTU_expcoeffs (6,0:max_user_angles,3)
   real(fpk), intent(out)   :: LRFE_MTU_expcoeffs (6,0:max_user_angles,3)       ! 3 in second dimension, could be used by Tmat

!  linearizations w.r.t. PSD parameters

   real(fpk), intent(out)   :: LPSD_MTU_dist (5,3)

!  Error handling

   logical, intent(out)       :: fail
   character*(*), intent(out) :: MTU_messages(3)

!  Local variables
!  ---------------

   integer    :: i,l,ldud,w,q,qdum,npars
   real(fpk)  :: DiffThreshold,Micron_wavelength
   character(len=4) :: RefID

   fail         = .false.
   MTU_messages = ' '

!  Read input from file

!mick mod 1/24/2017 - added aerosol file preamble capacity
   if ( DoPreamble(nAerMode) ) then
      do i=1,nPLines
         read(InFileUnit,*) !Preamble(i)
      enddo
      DoPreamble(nAerMode) = .false.
   endif

   read(InFileUnit,*) ! Divider: '--------------------------------------------'

   read(InFileUnit,*) ! Header : 'Wvl idx','Wvl (um)'
   if ( .not.RefWvl ) then
      read(InFileUnit,*) w,Micron_wavelength
   else
      read(InFileUnit,*) RefID,Micron_wavelength
   endif

   DiffThreshold = 1.0d-10
   if ( dabs(Micron_wavelength - Micron_wavelength_in) > DiffThreshold ) then
!mick mod 1/24/2017 - modified error handling
      !!Read unexpected wavelength --> abort!
      !MTU_messages(1) = 'Failure: Aerosol optical property input file, short version'
      !!MTU_messages(2) --> file name - defined outside aubroutine
      !MTU_messages(3) = 'Problem: Read unexpected wavelength'
      write(*,*)
      write(*,'(a)')      'Error: mismatch in wavelengths'
      write(*,'(a,f6.1)') 'Wavelength input to subroutine  (in um) = ',Micron_wavelength_in
      write(*,'(a,f6.1)') 'Wavelength input from aero file (in um) = ',Micron_wavelength
      !write(*,*) 'Difference Threshold                       = ',DiffThreshold
      !write(*,*) 'Magnitude of actual difference             = ',dabs(Micron_wavelength - Micron_wavelength_in)

      !Read unexpected wavelength --> abort!
      MTU_messages(1) = 'Failure: Aero optical property input file'
      !MTU_messages(2) --> file name - defined outside subroutine
      MTU_messages(3) = 'Problem: Mismatch between spectral cfg file & aero interpolate flag with aero file wvls'
      fail = .true. ; return
   end if

!  Linearization information headers

   Read(InFileUnit,*) ! Header 'Linearization flags/numbers/parameters'
   Read(InFileUnit,*)do_LinearREF, do_LinearPSD, NLIN, NPSD ; npars = NLIN + NPSD
   Read(InFileUnit,*)PARSLIST(1:npars)

!  general case
!  ============

   if ( .not.RefWvl ) then

!  Bulk quantities and linearizations
!  ----------------------------------

      read(InFileUnit,*) ! Header : 'Ext Coeff','Scat Coeff','SSA','Asym Par','Num of Expan Coeffs'
      read(InFileUnit,*) &
         MTU_bulk(1),&  ! Extinction coefficient
         MTU_bulk(2),&  ! Scattering coefficient
         MTU_bulk(3),&  ! Single scattering albedo
         MTU_asymm  ,&  ! Asymmetry parameter
         MTU_ncoeffs    ! Number of Expansion coefficients
      if ( MTU_ncoeffs > max_user_angles ) then
         !# of coefficients exceeds max allowed --> abort!
         MTU_messages(1) = 'Failure: Aerosol optical property input file'
         !MTU_messages(2) --> file name - defined outside subroutine
         MTU_messages(3) = 'Problem: # of coefficients exceeds max allowed for current CODEC settings'
         fail = .true. ; return
      end if

      if ( do_LinearRef ) then
         Read(InFileUnit,*)  ! Header : 'LRFE ExtCoeff','LRFE ScaCoeff','LRFE SSA','LRFE AsymPar'
         do q = 1, NLIN
            Read(InFileUnit,*)QDUM, LRFE_MTU_bulk(1:3,q), LRFE_MTU_asymm(q)
         enddo
      endif
      if ( do_LinearPSD ) then
         Read(InFileUnit,*)  ! Header : 'LPSD ExtCoeff','LPSD ScaCoeff','LPSD SSA','LPSD AsymPar'
         do q = 1, NPSD
            Read(InFileUnit,*)QDUM, LPSD_MTU_bulk(1:3,q), LPSD_MTU_asymm(q)
         enddo
      endif

!  Coefficients and linearizations
!  -------------------------------

      read(InFileUnit,*) ! Header : 'Moment','Expan Coeffs (a1,a2,a3,a4,b1,b2)'
      do l=0,MTU_ncoeffs
         read(InFileUnit,*) &
            ldud,(MTU_expcoeffs(i,l),i=1,6) !Index, Expansion coefficient input (a1,a2,a3,a4,b1,b2)
      enddo
      if ( do_LinearRef ) then
         read(InFileUnit,*) ! Header : 'WF#','Moment','Linearized RFE Expan Coeffs (a1,a2,a3,a4,b1,b2)'
         do l=0,MTU_ncoeffs
            do q = 1, NLIN
               read(InFileUnit,*)QDUM, ldud,(LRFE_MTU_expcoeffs(i,l,q),i=1,6)
            enddo
         enddo
      endif
      if ( do_LinearPSD ) then
         read(InFileUnit,*) ! Header : 'WF#','Moment','Linearized PSD Expan Coeffs (a1,a2,a3,a4,b1,b2)'
         do l=0,MTU_ncoeffs
            do q = 1, NPSD
               read(InFileUnit,*) QDUM, ldud,(LPSD_MTU_expcoeffs(i,l,q),i=1,6)
            enddo
         enddo
      endif

!  REFERENCE WAVELENGTH
!  ====================

   else

!  Bulk properties only

      read(InFileUnit,*) ! Header : 'Ext Coeff','Scat Coeff','SSA'
      read(InFileUnit,*) &
         MTU_bulk(1),&  ! Extinction coefficient
         MTU_bulk(2),&  ! Scattering coefficient
         MTU_bulk(3)    ! Single scattering albedo
      if ( do_LinearRef ) then
         Read(InFileUnit,*)  ! Header : 'LRFE ExtCoeff','LRFE ScaCoeff','LRFE SSA'
         do q = 1, NLIN
            Read(InFileUnit,*)QDUM, LRFE_MTU_bulk(1:3,q)
         enddo
      endif
      if ( do_LinearPSD ) then
         Read(InFileUnit,*)  ! Header : 'LPSD ExtCoeff','LPSD ScaCoeff','LPSD SSA'
         do q = 1, NPSD
            Read(InFileUnit,*)QDUM, LPSD_MTU_bulk(1:3,q)
         enddo
      endif

!  Particle size distribution

      read(InFileUnit,*) ! Header : 'Particle Size Distribution Characteristics'
      read(InFileUnit,*) ! Header : 'Normalization','Geom Xsection','Geom. Volume', 'REFF', 'VEFF'
      read(InFileUnit,*) MTU_dist(1:5)
      if ( do_LinearPSD ) then
         do q = 1, NPSD
            Read(InFileUnit,*) QDUM, LPSD_MTU_dist(1:5,q)
         enddo
      endif

   endif

end subroutine Read_aerosol_oprop_file_plus

!

subroutine Read_aerosol_oprop_file_plus_short (       &
    InFileUnit, Micron_wavelength_in, RefWvl,         & ! Input
    maxAerModes, nAerMode, nPLines, DoPreamble,       & ! I/O
    do_LinearRef, do_LinearPSD, NLIN, NPSD, PARSLIST, & ! Ouput
    fail, MTU_messages )                                ! Ouput

!  Introduced 10/19/16. Adapted from CODECTOOL routine of the same name.
!    Major new thing: aerosol distribution characteristics for Ref Wavl.

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   integer, intent(in)    :: InFileUnit
   real(fpk), intent(in)  :: Micron_wavelength_in
   logical, intent(in)    :: RefWvl

!mick mod 1/24/2017 - added aerosol file preamble capacity
   integer, intent(in)    :: maxAerModes, nAerMode, nPLines
   logical, intent(inout) :: DoPreamble(maxAerModes)

!  Linearization control

   logical  , intent(out)    :: do_LinearPSD, do_LinearRef
   integer  , intent(out)    :: NPSD, NLIN
   real(fpk), intent(out)    :: PARSLIST(10)

!  Error handling

   logical, intent(out)      :: fail
   character*(*), intent(out) :: MTU_messages(3)

!  Local variables
!  ---------------

   integer     :: i,w,npars
   real(fpk)   :: DiffThreshold,Micron_wavelength
   character*4 :: RefID

   fail         = .false.
   MTU_messages = ' '

!  Read input from file

!mick mod 1/24/2017 - added aerosol file preamble capacity
   if ( DoPreamble(nAerMode) ) then
      do i=1,nPLines
         read(InFileUnit,*) !Preamble(i)
      enddo
      DoPreamble(nAerMode) = .false.
   endif

   read(InFileUnit,*) ! Divider: '--------------------------------------------'

   read(InFileUnit,*) ! Header : 'Wvl idx','Wvl (um)'
   if ( .not.RefWvl ) then
      read(InFileUnit,*) w,Micron_wavelength
   else
      read(InFileUnit,*) RefID,Micron_wavelength
   endif
   DiffThreshold = 1.0d-10
   if ( dabs(Micron_wavelength - Micron_wavelength_in) > DiffThreshold ) then
!mick mod 1/24/2017 - modified error handling
      !!Read unexpected wavelength --> abort!
      !MTU_messages(1) = 'Failure: Aerosol optical property input file, short version'
      !!MTU_messages(2) --> file name - defined outside aubroutine
      !MTU_messages(3) = 'Problem: Read unexpected wavelength'
      write(*,*)
      write(*,'(a)')      'Error: mismatch in wavelengths'
      write(*,'(a,f6.1)') 'Wavelength input to subroutine  (in um) = ',Micron_wavelength_in
      write(*,'(a,f6.1)') 'Wavelength input from aero file (in um) = ',Micron_wavelength
      !write(*,*) 'Difference Threshold                       = ',DiffThreshold
      !write(*,*) 'Magnitude of actual difference             = ',dabs(Micron_wavelength - Micron_wavelength_in)

      !Read unexpected wavelength --> abort!
      MTU_messages(1) = 'Failure: Aero optical property input file'
      !MTU_messages(2) --> file name - defined outside subroutine
      MTU_messages(3) = 'Problem: Mismatch between spectral cfg file & aero interpolate flag with aero file wvls'
      fail = .true. ; return
   endif

!  Linearization information headers

   Read(InFileUnit,*) ! Header 'Linearization flags/numbers/parameters'
   Read(InFileUnit,*)do_LinearREF, do_LinearPSD, NLIN, NPSD ; npars = NLIN + NPSD
   Read(InFileUnit,*)PARSLIST(1:npars)

!  Done

   return
end subroutine Read_aerosol_oprop_file_plus_short

!

subroutine Write_aerosol_oprop_file_plus ( &
   max_user_angles, OutFileUnit, w, Micron_wavelength, RefWvl,       & ! Inputs
   maxAerModes, nAerMode, nPLines, DoPreamble, Preamble,             & ! I/O
   do_LinearRef, do_LinearPSD, NLIN, NPSD, PARSLIST,                 & ! Inputs
   MTU_dist, MTU_bulk, MTU_asymm, MTU_ncoeffs, MTU_expcoeffs,        & ! Inputs
   LRFE_MTU_bulk, LRFE_MTU_asymm, LRFE_MTU_expcoeffs,                & ! Inputs
   LPSD_MTU_bulk, LPSD_MTU_asymm, LPSD_MTU_expcoeffs, LPSD_MTU_dist  ) ! Inputs

!  Introduced 10/19/16. Adapted from CODECTOOL routine of the same name.
!    Major new thing: aerosol distribution characteristics for Ref Wavl.

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  basic

   integer, intent(in)    :: max_user_angles
   integer, intent(in)    :: OutFileUnit, w
   real(fpk), intent(in)  :: Micron_wavelength
   logical, intent(in)    :: RefWvl

!mick mod 1/24/2017 - added aerosol file preamble capacity
!  Rob Fix 3/7/17. Preamble is now mode-dependent
   integer, intent(in)           :: maxAerModes, nAerMode, nPLines
   logical, intent(inout)        :: DoPreamble(maxAerModes)
   character(len=*), intent(in)  :: Preamble(nPLines,maxAerModes)

!  Linearization control

   logical  , intent(in)    :: do_LinearPSD, do_LinearRef
   integer  , intent(in)    :: NPSD, NLIN
   real(fpk), intent(in)    :: PARSLIST(10)

!  Distribution parameters

   real(fpk), intent(in)  :: MTU_dist(5)

!  Bulk optical parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk), intent(in)  :: MTU_bulk(3)

!  Asymmetry parameter and number of expansion coefficients

   real(fpk), intent(in)  :: MTU_asymm

!  Number and value of expansion coefficients

   integer, intent(in)    :: MTU_ncoeffs
   real(fpk), intent(in)  :: MTU_expcoeffs(6,0:max_user_angles)

!  Linearized inputs

   real(fpk), intent(in)   :: LRFE_MTU_bulk (3,3)     ! 3 in second dimension, could be used by Tmat
   real(fpk), intent(in)   :: LPSD_MTU_bulk (3,3)
  
   real(fpk), intent(in)   :: LRFE_MTU_asymm(3)       ! 3 in second dimension, could be used by Tmat
   real(fpk), intent(in)   :: LPSD_MTU_asymm(3)

   real(fpk), intent(in)   :: LPSD_MTU_expcoeffs (6,0:max_user_angles,3)
   real(fpk), intent(in)   :: LRFE_MTU_expcoeffs (6,0:max_user_angles,3)       ! 3 in second dimension, could be used by Tmat

!  linearizations w.r.t. PSD parameters

   real(fpk), intent(in)   :: LPSD_MTU_dist (5,3)

!  Local variables
!  ---------------

   integer :: i,l,q,npars, Preamble_unit
   character*2 :: c2

!  Write output to file

!mick mod 1/24/2017 - added aerosol file preamble capacity
!Rob mod  3/3/2017  -  always print out preamble file to UserAeroFiles/out/
!  Rob Fix 3/7/17   -  Preamble is now mode-dependent

   if ( DoPreamble(nAerMode) ) then
      Preamble_unit = 1111 ; c2 = '00' ; write(c2(2:2),'(I1)')nAerMode
      open(Preamble_unit,file='UserAeroFiles/out/AerFile_Preamble_Mode_'//c2//'.dat',status='replace')
      do i=1,nPLines
         write(Preamble_unit,'(a)') trim(Preamble(i,nAerMode))
         write(OutFileUnit,'(a)') trim(Preamble(i,nAerMode))
      enddo
      DoPreamble(nAerMode) = .false.
      close(Preamble_unit)
   endif

   write(OutFileUnit,'(3x,a)') '--------------------------------------------' // &
                               '--------------------------------------------' // &
                               '--------------------------------------------'
   write(OutFileUnit,'(6x,a,7x,a)') 'Wvl idx','Wvl (um)'

   if ( .not.RefWvl ) then

!  Headers
!  -------

      write(OutFileUnit,'(8x,i4.4,e20.11)') w,Micron_wavelength
      write(OutFileUnit,'(6x,a)')          'Linearization flags/Numbers/Parameters --->   ' ; NPARS = NLIN + NPSD
      write(OutFileUnit,'(5x,2L2,4x,2I2)')  do_LinearREF, do_LinearPSD, NLIN, NPSD
      write(OutFileUnit,'(3x,6e20.11)')     PARSLIST(1:NPARS)

!  Bulk properties
!  ---------------

      !Sample:
      !      Ext Coeff           Scat Coeff              SSA               Asym Par       Num of Expan Coeffs
      !   0.30186852870E-01   0.29730300500E-01   0.98487578775E+00   0.59090185717E+00           0
      write(OutFileUnit,'(9x,a,11x,a,14x,a,15x,a,7x,a)') 'Ext Coeff','Scat Coeff','SSA','Asym Par','Num of Expan Coeffs'
      write(OutFileUnit,'(3x,4e20.11,9x,i4)') &  ! AMS 09 JAN 2015 Updated expansion coeffs from i3 to i4 as sometimes we have >1000.
         MTU_bulk(1),& ! Extinction coefficient
         MTU_bulk(2),& ! Scattering coefficient
         MTU_bulk(3),& ! Single scattering albedo
         MTU_asymm  ,& ! Asymmetry parameter
         MTU_ncoeffs   ! Number of Expansion coefficients

      if ( do_LinearRef ) then
         !Sample:
         ! WF#      LRFE ExtCoeff       LRFE ScaCoeff          LRFE SSA           LRFE AsymPar
         !  1     0.30186852870E-01   0.29730300500E-01   0.98487578775E+00    0.59090185717E+00
         !  2     0.30186852870E-01   0.29730300500E-01   0.98487578775E+00    0.59090185717E+00 ...etc.
         write(OutFileUnit,'(1x,a,4x,a,6x,a,10x,a,11x,a)') 'WF#','LRFE ExtCoeff','LRFE ScaCoeff','LRFE SSA','LRFE AsymPar'
         do q = 1, NLIN
            write(OutFileUnit,'(1x,i2,4e20.11)') Q, LRFE_MTU_bulk(1:3,q), LRFE_MTU_asymm(q)
         enddo
      endif
      if ( do_LinearPSD ) then
         !Sample:
         ! WF#      LPSD ExtCoeff       LPSD ScaCoeff          LPSD SSA           LPSD AsymPar
         !  1     0.30186852870E-01   0.29730300500E-01   0.98487578775E+00    0.59090185717E+00
         !  2     0.30186852870E-01   0.29730300500E-01   0.98487578775E+00    0.59090185717E+00 ...etc.
         write(OutFileUnit,'(1x,a,4x,a,6x,a,10x,a,11x,a)') 'WF#', 'LPSD ExtCoeff','LPSD ScaCoeff','LPSD SSA','LPSD AsymPar'
         do q = 1, NPSD
            write(OutFileUnit,'(1x,i2,4e20.11)') Q, LPSD_MTU_bulk(1:3,q), LPSD_MTU_asymm(q)
         enddo
      endif

!  Coefficients
!  ------------

      !Sample:
      !      Moment                                             Expan Coeffs (a1,a2,a3,a4,b1,b2)
      !       0000   0.10000000000E+01   0.00000000000E+00   0.00000000000E+00   0.85192295256E+00   0.00000000000E+00  -0.00000000000E+00
      write(OutFileUnit,'(6x,a,45x,a)') 'Moment','Expan Coeffs (a1,a2,a3,a4,b1,b2)'
      do l=0,MTU_ncoeffs
         write(OutFileUnit,'(7x,i4.4,6e20.11)') &
            l,(MTU_expcoeffs(i,l),i=1,6) ! Index, Expansion coefficient output (a1,a2,a3,a4,b1,b2)
      enddo

      if ( do_LinearRef ) then
         !Sample:
         ! WF#  Moment                                     Linearized RFE Expan Coeffs (a1,a2,a3,a4,b1,b2)
         !  1    0000   0.10000000000E+01   0.00000000000E+00   0.00000000000E+00   0.85192295256E+00   0.00000000000E+00  -0.00000000000E+00
         !  2    0000   0.10000000000E+01   0.00000000000E+00   0.00000000000E+00   0.85192295256E+00   0.00000000000E+00  -0.00000000000E+00
         write(OutFileUnit,'(1x,a,2x,a,37x,a)') 'WF#','Moment','Linearized RFE Expan Coeffs (a1,a2,a3,a4,b1,b2)'
         do l=0,MTU_ncoeffs
            do q = 1, NLIN
               write(OutFileUnit,'(1x,i2,4x,i4.4, 6e20.11)') Q, l, (LRFE_MTU_expcoeffs(i,l,q),i=1,6)
            enddo
         enddo
      endif

      if ( do_LinearPSD ) then
         !Sample:
         ! WF#  Moment                                     Linearized PSD Expan Coeffs (a1,a2,a3,a4,b1,b2)
         !  1    0000   0.10000000000E+01   0.00000000000E+00   0.00000000000E+00   0.85192295256E+00   0.00000000000E+00  -0.00000000000E+00
         !  2    0000   0.10000000000E+01   0.00000000000E+00   0.00000000000E+00   0.85192295256E+00   0.00000000000E+00  -0.00000000000E+00
         write(OutFileUnit,'(1x,a,2x,a,37x,a)') 'WF#','Moment','Linearized PSD Expan Coeffs (a1,a2,a3,a4,b1,b2)'
         do l=0,MTU_ncoeffs
            do q = 1, NPSD
               write(OutFileUnit,'(1x,i2,4x,i4.4,6e20.11)') Q, l,(LPSD_MTU_expcoeffs(i,l,q),i=1,6)
            enddo
         enddo
      endif

!  REFERENCE WAVELENGTH ONLY
!  =========================

   else

!  Headers
!  -------

      write(OutFileUnit,'(8x,a,e20.11)')  'Ref#',Micron_wavelength
      write(OutFileUnit,'(6x,a)')         'Linearization flags/Numbers/Parameters --->   ' ; NPARS = NLIN + NPSD
      write(OutFileUnit,'(5x,2L2,4x,2I2)') do_LinearREF, do_LinearPSD, NLIN, NPSD
      write(OutFileUnit,'(3x,6e20.11)')    PARSLIST(1:NPARS)

!  Bulk properties only
!  --------------------

      !Sample:
      !      Ext Coeff           Scat Coeff              SSA
      !   0.30186852870E-01   0.29730300500E-01   0.98487578775E+00
      write(OutFileUnit,'(9x,a,11x,a,14x,a)') 'Ext Coeff','Scat Coeff','SSA'
      write(OutFileUnit,'(3x,4e20.11,9x,i3)') &
         MTU_bulk(1),& ! Extinction coefficient
         MTU_bulk(2),& ! Scattering coefficient
         MTU_bulk(3)   ! Single scattering albedo

      if ( do_LinearRef ) then
         !Sample:
         ! WF#      LRFE ExtCoeff       LRFE ScaCoeff          LRFE SSA
         !  1     0.30186852870E-01   0.29730300500E-01   0.98487578775E+00
         !  2     0.30186852870E-01   0.29730300500E-01   0.98487578775E+00 ...etc.
         write(OutFileUnit,'(1x,a,4x,a,6x,a,10x,a)') 'WF#','LRFE ExtCoeff','LRFE ScaCoeff','LRFE SSA'
         do q = 1, NLIN
            write(OutFileUnit,'(1x,i2,3e20.11)') Q, LRFE_MTU_bulk(1:3,q)
         enddo
      endif

      if ( do_LinearPSD ) then
         !Sample:
         ! WF#      LPSD ExtCoeff       LPSD ScaCoeff          LPSD SSA
         !  1     0.30186852870E-01   0.29730300500E-01   0.98487578775E+00
         !  2     0.30186852870E-01   0.29730300500E-01   0.98487578775E+00 ...etc.
         write(OutFileUnit,'(1x,a,4x,a,6x,a,10x,a)') 'WF#','LPSD ExtCoeff','LPSD ScaCoeff','LPSD SSA'
         do q = 1, NPSD
            write(OutFileUnit,'(1x,i2,3e20.11)') Q, LPSD_MTU_bulk(1:3,q)
         enddo
      endif

!  Size Distribution parameters and linearizations (for reference only)
!  -----------------------------------------------

      !Sample:
      !                ----------------Particle Size Distribution Characteristics----------------
      ! WF#      Normalization       Geom Xsection        Geom Volume      Effective radius    Effective variance
      !        0.30186852870E-01   0.29730300500E-01   0.98487578775E+00   0.98487578775E+00   0.98487578775E+00
      !  1     0.30186852870E-01   0.29730300500E-01   0.98487578775E+00   0.98487578775E+00   0.98487578775E+00
      !  2     0.30186852870E-01   0.29730300500E-01   0.98487578775E+00   0.98487578775E+00   0.98487578775E+00....etc
      write(OutFileUnit,'(17x,a)') '----------------Particle Size Distribution Characteristics----------------'
      write(OutFileUnit,'(1x,a,4x,a,7x,a,8x,a,5x,a,4x,a)') 'WF#','Normalization','Geom Xsection','Geom. Volume',&
                                                           'Effective radius','Effective variance'
      write(OutFileUnit,'(3x,5e20.11)') MTU_dist(1:5)

      if ( do_LinearPSD ) then
         do q = 1, NPSD
            write(OutFileUnit,'(1x,i2,5e20.11)') Q, LPSD_MTU_dist(1:5,q)
         enddo
      endif

   endif

end subroutine Write_aerosol_oprop_file_plus

!  End module

end module GEMSTOOL_Read_Configfiles_m
