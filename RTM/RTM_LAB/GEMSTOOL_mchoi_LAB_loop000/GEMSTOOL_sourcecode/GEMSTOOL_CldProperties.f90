Module GEMSTOOL_CldProperties_m

!  Use Type structures (GEMSTOOL)

   use GEMSTOOL_Input_Types_m


!  Use Mie code
!  ------------

!  parameters

   use RTSMie_parameters_m

!  Single call

   use RTSMie_sourcecode_m

!  All routines are public

private
public :: GEMSTOOL_cld_PROPERTIES

contains
!        PSD_Index, PSD_pars, n_real, n_imag, FixR1R2, R1, R2, epsnum, epsfac,   & ! INPUT, Mie/Tmat Parameters

subroutine GEMSTOOL_CLD_PROPERTIES &
      ( maxlayers, maxwav, maxcldmoms, interpolate_clouds, do_wavnums,                 & ! Dimensions and Flags 
        nlayers, nmuller, nwav, wav, height_grid, GEMSTOOL_INPUTS, momsize_cutoff,     & ! Inputs
        cldlayerflags, Loading, n_scatmoms, cldtau_unscaled,                           & ! OUTPUT, cloud control
        cod_scaling, cloud_deltau, cloud_ssalbs, cloud_scatmoms,                       & ! OUTPUT, cloud optical
        fail1, Messages_Optical )

!  =============================================================================
!                            CLOUD REGULAR CREATION
!  =============================================================================

!  CLoud Loading = Uniform cloud
!  Generation of Cloud optical properties
!     1. Call to the Mie program
!     2. Convert Mie output (Microsopic) to IOP output (macroscopic) 

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)
  
!  Input variables
!  ---------------

!  External Dimensioning

   integer, INTENT (IN) :: maxlayers, maxwav
   integer, INTENT (IN) :: maxcldmoms

!  Flag for interpolation of clouds

   logical, INTENT(in)  :: interpolate_clouds

!  Flag for using wavnumber output

   logical, INTENT(in)  :: do_wavnums

!  Numbers

   integer, INTENT (IN) ::  nlayers
   integer, INTENT (IN) ::  nwav

!  Heights and wavelengths

   REAL    (fpk),    INTENT (IN)   :: wav(maxwav)
   REAL    (fpk),    INTENT (IN)   :: height_grid(0:maxlayers)

!  nmuller = 1 (scalar code), = 6 (vector code)

   integer, INTENT (IN) ::  nmuller

!  cloud moment size cutoff. DEFAULT = 0.001

   REAL    (fpk),    INTENT (IN)   :: momsize_cutoff

!  Type Structure inputs

   TYPE(GEMSTOOL_Config_Inputs) :: GEMSTOOL_INPUTS

!  Mie PSD inputs
!  ==============

!  PSD inputs (distribution index, PSD parameters)

!      psd_Index      - Index for particle size distribution of spheres
!      psd_pars       - Parameters characterizing PSD (up to 3 allowed)

!    PSD_index = 1 : TWO-PARAMETER GAMMA with alpha and b given
!    PSD_index = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given
!    PSD_index = 3 : BIMODAL GAMMA with equal mode weights
!    PSD_index = 4 : LOG-NORMAL with rg and sigma given
!    PSD_index = 5 : LOG-NORMAL with reff and veff given
!    PSD_index = 6 : POWER LAW
!    PSD_index = 7 : MODIFIED GAMMA with alpha, rc and gamma given
!    PSD_index = 8 : MODIFIED GAMMA with alpha, b and gamma given

!  FixR1R2 : If  set, Use Internal routine to calculate R1 and R2 (outputs)
!            If Not set, Use Input R1 and R2 for PSD limits.
!  R1, R2         - Minimum and Maximum radii (Microns)
!  N_REAL, N_IMAG - real and imaginary parts, refractive index (N-i.GE.0)

!  Mie-specific inputs
!  ===================

!  Limiting particle size value. Set to 10000.0 default.
!   If you exceed this, program will tell you to increase dimensioning.

!  R1R2_cutoff particle size for setting R1 and R2 internally

!  PSD quadrature control
!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.

!  output
!  ======

!  Loading output

    real(fpk),    dimension ( maxlayers ),  INTENT (OUT)      :: Loading

!  cloud layer flags

    LOGICAL,      DIMENSION ( maxlayers ), intent(out)  :: cldLAYERFLAGS 

!  AOP output: Number of exapnsion coefficients to be used

   integer, INTENT (OUT)   ::  n_scatmoms

!  Unscaled profiles of optical depth

    real(fpk),    INTENT (OUT)  :: cldtau_unscaled(maxlayers)

!  cod scaling

    real(fpk),    INTENT (OUT)  :: cod_scaling(maxwav)

!  cod output, final

    REAL(fpk),    DIMENSION( maxlayers, maxwav )      , INTENT (OUT)  :: cloud_DELTAU

!  AOPS

    REAL(fpk),    DIMENSION( maxwav )                 , INTENT (OUT)  :: cloud_SSALBS 
    REAL(fpk),    DIMENSION( 6, 0:maxcldmoms, maxwav ), INTENT (OUT)  :: cloud_SCATMOMS

!  Exception handling

   logical,        INTENT (OUT)           :: fail1
   character*(*),  INTENT (OUT)           :: Messages_Optical(4)

!  LOCAL VARIABLES
!  @@@@@@@@@@@@@@@

!  Mie/Tmatrix LOCAL INPUT variables (Not part of module input)
!  ============================================================

!      Do_Expcoeffs      - Boolean flag for computing Expansion Coefficients
!      Do_Fmatrix        - Boolean flag for computing F-matrix at equal-angles

   logical    :: Do_Expcoeffs
   logical    :: Do_Fmatrix

!      Do_Monodisperse   - Boolean flag for Doing a Monodisperse calculation
!                          If set, the PSD stuff will be turned off internally

   LOGICAL    :: do_Monodisperse

!  F-matrix Angular control input (NOT REQUIRED HERE)

!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles. (NPNA)
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

!  NPNA - number of equidistant scattering angles (from 0
!             to 180 deg) for which the scattering matrix is
!             calculated.

!    NOT REQUIRED HERE

   INTEGER           :: n_Fmatrix_angles
   REAL    (fpk)     :: Fmatrix_angles(max_Mie_angles)

!  Monoradius     - Monodisperse radius size (Microns)

   real    (fpk)     :: Monoradius

!  Mie code OUTPUT variables
!  =========================

!  Bulk distribution parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk)    :: Smie_bulk (3)

!  Expansion coefficients and Asymmetry parameter

   integer      :: Smie_ncoeffs
   real(fpk)    :: Smie_expcoeffs (6,0:max_Mie_angles)
   real(fpk)    :: Smie_asymm

!  F-matrix,  optional output

   real(fpk)    :: Smie_Fmatrix(4,max_Mie_angles)

!  Distribution parameters
!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(fpk)     :: Smie_dist (5)

!  Exception handling

   character*90  :: Mie_Bmessages(3)

!  Local variables
!  ===============

!  Local cloud properties ! New, @@@ RTS 09/11/12

!   integer   :: PSD_index(2)              ! New, @@@ RTS 09/11/12
!   real(fpk) :: PSD_pars(3,2)             ! New, @@@ RTS 09/11/12
!   real(fpk) :: n_real(2), n_imag(2)      ! New, @@@ RTS 09/11/12

!  COD output: Reference values
!    * These are the values at reference wavelength w0
!    * They are set for wavelength w0, then used again

    real(fpk)     :: extinction_ref

!  Help Variables

   character*3  :: cwav
   integer      :: n, k, l, m, istatus, point_index, w, wc, wavmask(maxwav)
   integer      :: n_scatmoms_w, n_cldwavs, local_nwav, wa, wa1, wa2, wastart
   real(fpk)    :: wavelength, Mie_wavelength, extinction, fa1, fa2, lam, lamstart, lamfinish, cloud_ext, diff
   logical      :: trawl

!  Interpolation arrays
!     UV-Vis : 56  is enough wavelengths for 10 nm intervals over a 250-800  nm range
!     SWIr   : 136 is enough wavelengths for 10 nm intervals over a 750-2100 nm range
!    NOTE    : CAN MAKE THESE ARRAYS ALLOCATABLE ----------THINK ABOUT IT !!!!!!!!!

   real(fpk)    :: cldwav(136),local_codscaling(136),local_cldssalbs(136), local_scatmoms(6,0:maxcldmoms,136)

!  Local control

   logical, parameter :: do_iopchecker    = .false.
!   logical, parameter :: do_iopchecker    = .true.
   logical, parameter :: do_cld_Jacobians = .false.

!  Initialize output
!  =================

!  Initialize exception handling

   fail1 = .false.
   Messages_Optical  = ' '

!  initialize cloud Loading

   Loading        = 0.0d0
   cldlayerflags = .false.

!  Set local Mie/Tmatrix variables

   do_monodisperse = .false.             ! Always false
   Do_Fmatrix      = .false.             ! Always false
   DO_Expcoeffs    = .false.             ! this will be set later on
   monoradius      = 0.0d0

!  Initialize optical properties

   extinction_ref    = 0.0d0
   cod_scaling       = 0.0d0
   cloud_deltau    = 0.0d0
   cloud_ssalbs    = 0.0d0
   cloud_scatmoms  = 0.0d0
   cldtau_unscaled   = 0.0d0

!  Initialize local cloud properties (Just a precaution)! New, @@@ RTS 09/11/12

!   n_real = 0.0_fpk ; n_imag   = 0.0_fpk
!   PSDIndex = 0    ; PSD_pars = 0.0_fpk

!  Now form the cloud Loading (uniform)
!  ------------------------------------

!  Extinction

   Cloud_ext = GEMSTOOL_INPUTS%Clouds%cldtau_input_w1 / &
                   ( GEMSTOOL_INPUTS%Clouds%ctop_height - GEMSTOOL_INPUTS%Clouds%cbot_height )

!  Fill up cloud layers

   do n = 1, nlayers
      if (n.gt.GEMSTOOL_INPUTS%Clouds%ctop_level .and. n.le.GEMSTOOL_INPUTS%Clouds%cbot_level) then
         cldlayerflags(n) = .true.
         diff = height_grid(n-1) - height_grid(n)
         loading(n) = diff * Cloud_ext
         cldtau_unscaled(n) = loading(n)
      endif
   enddo

!  Interpolation setup
!  ===================

   if ( interpolate_clouds ) then
      trawl = .true. ; wa = 0
      if ( .not. do_wavnums ) then
         lamstart = 200.0d0 ; lam = lamstart        !  Smallest wavelength is 200 nm (UVN application)
         do while (trawl)
            lam = lam + 20.0d0
            if ( lam.ge.wav(1) ) then
               wa = wa + 1 ; cldwav(wa) = lam - 20.0d0
               if ( lam .gt. wav(nwav) ) then
                  wa = wa + 1 ; cldwav(wa) = lam ; trawl = .false.
               endif
            endif
         enddo
          n_cldwavs = wa ; local_nwav = n_cldwavs
      else
         lamfinish = 2100.0d0 ; lam = lamfinish     !  Largest wavelength is 2100 nm (NSW application)
         do while (trawl)
            lam = lam - 20.0d0
            if ( lam.le.1.0d+07/wav(1) ) then
               wa = wa + 1 ; cldwav(wa) = lam + 20.0d0
               if ( lam .lt. 1.0d+07/wav(nwav) ) then
                  wa = wa + 1 ; cldwav(wa) = lam ; trawl = .false.
               endif
            endif
         enddo
          n_cldwavs = wa ; local_nwav = n_cldwavs
      endif
   else
      local_nwav = nwav
   endif

!  Check to see if reference wavelength is one of the set. Initialize Mask.

   point_index = 0
   if ( .not. do_wavnums ) then
      if ( interpolate_clouds ) then
         do w = 1, n_cldwavs
            wavmask(w) = w ; if ( cldwav(w) .eq. GEMSTOOL_INPUTS%Clouds%reference_w1 ) point_index = w
         enddo
      else
         do w = 1, nwav
            wavmask(w) = w ; if ( wav(w)    .eq. GEMSTOOL_INPUTS%Clouds%reference_w1 ) point_index = w
         enddo
      endif
   else
      if ( interpolate_clouds ) then
         do w = 1, n_cldwavs
            wavmask(w) = w ; if ( cldwav(w) .eq. GEMSTOOL_INPUTS%Clouds%reference_w1 ) point_index = w
         enddo
      else
         do w = 1, nwav
            wavmask(w) = w
         enddo
      endif
   endif

!   write(*,*)local_nwav, n_cldwavs, lambda_index

!  Mask to use if reference wavelength is one of list of wavelengths. [UVN only]

   if ( point_index .ne. 0 ) then
     wavmask(1) = point_index
     wc = 1
     do w = 1, point_index - 1
       wc = wc + 1 ; wavmask(wc) = w
     enddo
     do w = point_index + 1, local_nwav
       wc = wc + 1 ; wavmask(wc) = w
     enddo
   endif

!  debug
!   do w = 1, local_nwav
!      write(*,*)w,wavmask(w),cldwav(w)
!   enddo
!   pause'interp'

!  Initialize

   n_scatmoms = 50000

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      M I E   C a l c u l a t i o n 
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Prepare the reference-wavelength Mie Inputs
!  ===========================================

   if ( point_index .eq. 0 ) then

!  Only require extinction coefficient if flagged
!  Set the local Mie program inputs (bulk properties only)

      Do_Expcoeffs     = .FALSE.

!  reference wavelength

      wavelength     = GEMSTOOL_INPUTS%Clouds%reference_w1
      Mie_wavelength = wavelength/1000.0d0

!  progress

      write(*,'(4x,a,5x,f12.5)') 'Regular Mie cloud calculation - Doing reference wavelength   : ',wavelength

!  SINGLE CALL

      call RTSMie_main  & !---MIE CALL
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,                                                                   & ! I
                GEMSTOOL_INPUTS%Clouds%CLD_PSDIndex,     GEMSTOOL_INPUTS%Clouds%CLD_PSDpars, MonoRadius,                     & ! I
                GEMSTOOL_INPUTS%Clouds%CLD_R1,           GEMSTOOL_INPUTS%Clouds%CLD_R2,  GEMSTOOL_INPUTS%Clouds%CLD_FixR1R2, & ! I
                GEMSTOOL_INPUTS%Clouds%CLD_nblocks,      GEMSTOOL_INPUTS%Clouds%CLD_nweights,                                & ! I
                GEMSTOOL_INPUTS%Clouds%CLD_xparticle_limit,  GEMSTOOL_INPUTS%Clouds%CLD_R1R2_cutoff,                         & ! I
                n_Fmatrix_angles, Fmatrix_angles, Mie_wavelength,                                                            & ! I
                GEMSTOOL_INPUTS%Clouds%CLD_nreal,       GEMSTOOL_INPUTS%Clouds%CLD_nimag,                                    & ! I
                Smie_bulk, Smie_asymm, Smie_ncoeffs, Smie_expcoeffs, Smie_Fmatrix, Smie_dist,                    & ! O
                fail1, istatus, Mie_Bmessages(1), Mie_Bmessages(2), Mie_Bmessages(3) )                             ! O

!  Exception handling on everything

      if ( Fail1 ) then  
         do m = 1, 3   
            Messages_Optical(m) = adjustl(trim(Mie_Bmessages(m)))
         enddo
         Messages_Optical(4) = 'Single Mie Call'
         return
      endif

!  Set the reference quantity

      extinction_ref    = Smie_bulk(1)

!  End Mie reference-wavelength calculation

   endif

!  Prepare General (all-wavelength, all-wavenumber) Mie inputs
!  ===========================================================

!  wavelength loop. First wavelength will be the reference, if in  list.
!  Wavnumber  loop. 

   do wc = 1, local_nwav

!  Wavelengths [nm]

      w = wavmask(wc)
      if ( interpolate_clouds ) then
         wavelength = cldwav(w)
      else
        if ( do_wavnums      ) wavelength = 1.0d+07/wav(w)
        if ( .not.do_wavnums ) wavelength = wav(w)
      endif

!  progress

      write(*,'(4x,a,1x,i3,1x,f12.5)') 'Regular Mie cloud calculation - Doing point, wavelength      : ', wc, wavelength

!  wavelength for the Mie code (micron unit)  

      Mie_wavelength = wavelength/1000.0d0

!  Set the local Mie program inputs (general)

      Do_Expcoeffs     = .TRUE.

!  BIMODAL vs. SINGLE CALL


      call RTSMie_main  & !---MIE CALL
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,                                                                   & ! I
                GEMSTOOL_INPUTS%Clouds%CLD_PSDIndex,     GEMSTOOL_INPUTS%Clouds%CLD_PSDpars, MonoRadius,                     & ! I
                GEMSTOOL_INPUTS%Clouds%CLD_R1,           GEMSTOOL_INPUTS%Clouds%CLD_R2,  GEMSTOOL_INPUTS%Clouds%CLD_FixR1R2, & ! I
                GEMSTOOL_INPUTS%Clouds%CLD_nblocks,      GEMSTOOL_INPUTS%Clouds%CLD_nweights,                                & ! I
                GEMSTOOL_INPUTS%Clouds%CLD_xparticle_limit,  GEMSTOOL_INPUTS%Clouds%CLD_R1R2_cutoff,                         & ! I
                n_Fmatrix_angles, Fmatrix_angles, Mie_wavelength,                                                            & ! I
                GEMSTOOL_INPUTS%Clouds%CLD_nreal,       GEMSTOOL_INPUTS%Clouds%CLD_nimag,                                    & ! I
                Smie_bulk, Smie_asymm, Smie_ncoeffs, Smie_expcoeffs, Smie_Fmatrix, Smie_dist,                    & ! O
                fail1, istatus, Mie_Bmessages(1), Mie_Bmessages(2), Mie_Bmessages(3) )                             ! O

!  Exception handling on everything

      if ( Fail1 ) then
         write(cwav,'(I3)')wc
         do m = 1, 3
            Messages_Optical(m) = adjustl(trim(Mie_Bmessages(m)))
         enddo
         Messages_Optical(4) = 'Single Mie Call'
         return
      endif
      
!  Set the reference quantities, if reference wavelength is in the list.
!   Values for the first (masked) wavelength

      if ( point_index .ne. 0 .and. wc.eq.1 ) then
         extinction_ref = Smie_bulk(1)
      endif

!  For the interpolation case, save output to local arrays

      if ( interpolate_clouds ) then
         local_codscaling(w)    = Smie_bulk(1) / extinction_ref
         local_cldssalbs(w)     = Smie_bulk(3)
         l = 0 ; local_scatmoms(1,0,w) = 1.0d0

         if ( GEMSTOOL_INPUTS%Clouds%do_Henyey_Greenstein ) then
            do while (local_scatmoms(1,l,w).gt.momsize_cutoff.and.l.lt.maxcldmoms )
               l = l + 1 ; local_scatmoms(1,l,w) = real(2*L+1,fpk) * Smie_asymm ** L     ! Henyey-Greenstein Clouds
            enddo
            n_scatmoms_w = l
            do l = 0, n_scatmoms_w
                local_scatmoms(2:nmuller,l,w) = 0.0_fpk                                  ! Henyey-Greenstein Clouds
            enddo
         else
            do while (local_scatmoms(1,l,w).gt.momsize_cutoff.and.l.lt.maxcldmoms )
               l = l + 1 ; local_scatmoms(1,l,w) = Smie_expcoeffs(1,l)
            enddo
            n_scatmoms_w = l
            do l = 0, n_scatmoms_w
               local_scatmoms(2:nmuller,l,w) = Smie_expcoeffs(2:nmuller,l) 
            enddo
         endif

         n_scatmoms = min(n_scatmoms, n_scatmoms_w)
         go to 677   
      endif

!  NOW set the optical property output
!  ===================================

!  Extinction and its scaling factor. All wavelengths.

      extinction = Smie_bulk(1)
      cod_scaling(w) = extinction / extinction_ref

!  Assign SSAs and expansion coefficients, single/bimodal cloud type

      cloud_ssalbs(w) = Smie_bulk(3)
      l = 0
      cloud_scatmoms(1,0,w) = 1.0d0

      if ( GEMSTOOL_INPUTS%Clouds%do_Henyey_Greenstein ) then
         do while (local_scatmoms(1,l,w).gt.momsize_cutoff.and.l.lt.maxcldmoms )
            l = l + 1 ; cloud_scatmoms(1,l,w) = real(2*L+1,fpk) * Smie_asymm ** L     ! Henyey-Greenstein Clouds
         enddo
         n_scatmoms_w = l
         do l = 0, n_scatmoms_w
            cloud_scatmoms(2:nmuller,l,w) = 0.0_fpk                                  ! Henyey-Greenstein Clouds
         enddo
      else
         do while (local_scatmoms(1,l,w).gt.momsize_cutoff.and.l.lt.maxcldmoms )
            l = l + 1 ; cloud_scatmoms(1,l,w) = Smie_expcoeffs(1,l)
         enddo
         n_scatmoms_w = l
         do l = 0, n_scatmoms_w
            cloud_scatmoms(2:nmuller,l,w) = Smie_expcoeffs(2:nmuller,l) 
         enddo
      endif

!  Update n_scatmoms

      n_scatmoms = min(n_scatmoms, n_scatmoms_w)

!  Apply scalings to loadings

      do n = 1, nlayers
         cloud_deltau(n,w) = cldtau_unscaled(n) * cod_scaling(w)
      enddo

!  debug cloud optical properties. VERSION TWO only

      if ( do_iopchecker ) then
         do n = 1, nlayers
           if (cldlayerflags(N).and.n.eq.107 ) then
              write(999,'(i4,1p6e20.10)')n,cloud_deltau(n,w),cloud_ssalbs(w)
              do l = 0, n_scatmoms_w
                write(999,'(2i5,1p6e20.10)')n,l,(cloud_scatmoms(k,l,w),k=1,1)
              enddo
!            else
!              write(999,'(i4,1p6e20.10)')n,cloud_deltau(n,w)
           endif
         enddo
!           pause'Reg 999'
      endif

!  continuation point for avoiding the exact monochromatic solution

677   continue

!  End wavelength/wavenumber loop

   enddo

!  Monochromatic solution, return

   if ( .not. interpolate_clouds ) return

!  Interpolation, wavelength regime

   if ( .not. do_wavnums ) then
     wastart = 1
     do w = 1, nwav
       wa = wastart ; trawl = .true.
!       write(*,*)w,wav(w),wa,cldwav(wa),cldwav(wa+1)
       do while (trawl)
         if ( wav(w) .ge. cldwav(wa) .and. wav(w) .le.cldwav(wa+1) ) trawl = .false.
       enddo
       wa1 = wa ; wa2 = wa + 1 ; fa1 = ( cldwav(wa2) - wav(w) ) / ( cldwav(wa2) -  cldwav(wa1) ) ; fa2 = 1.0d0 - fa1
       wastart = wa1
       if ( w.lt.nwav) then
          if(wav(w+1).ge.cldwav(wa+1))wastart = wa2
       endif
       cod_scaling(w) = fa1 * local_codscaling(wa1) + fa2 * local_codscaling(wa2)
       do n = 1, nlayers
         cloud_deltau(n,w) = cldtau_unscaled(n) * cod_scaling(w)
       enddo
       cloud_ssalbs(w) = fa1 * local_cldssalbs(wa1) + fa2 * local_cldssalbs(wa2)
       do l = 0, n_scatmoms
         do k = 1, nmuller
            cloud_scatmoms(1:nmuller,l,w) = fa1 * local_scatmoms(1:nmuller,l,wa1) + fa2 * local_scatmoms(1:nmuller,l,wa2)
         enddo
!        if ( w.eq.1.and.l.lt.50)write(*,*)l,cloud_scatmoms(1:nmuller,l,w) 
       enddo
       cloud_scatmoms(1,0,w) = 1.0d0
!       write(*,*)w,wav(w),fa1,fa2,cloud_ssalbs(w),n_scatmoms,sum(cloud_deltau(1:nlayers,w))
     enddo
   endif

!  Interpolation, wavenumber regime
!   12 August 2013 --> KLUTZY CODE here,,,,,,,,,,Improve it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if ( do_wavnums ) then
     wastart = 1
     do w = 1, nwav
       wa = wastart ; trawl = .true. ; lam =  1.0d+07/wav(w)
54     continue
       do while (trawl)
         if ( lam .lt. cldwav(wa) ) then
            if ( lam .ge. cldwav(wa+1) ) then
               trawl = .false.
            else if (lam .lt. cldwav(wa+1) ) then
               wa = wa + 1 ; go to 54
            endif
         endif
       enddo
       wa1 = wa ; wa2 = wa + 1 ; fa1 = ( cldwav(wa2) - lam ) / ( cldwav(wa2) -  cldwav(wa1) ) ; fa2 = 1.0d0 - fa1
       wastart = wa1
       cod_scaling(w) = fa1 * local_codscaling(wa1) + fa2 * local_codscaling(wa2)
       do n = 1, nlayers
         cloud_deltau(n,w) = cldtau_unscaled(n) * cod_scaling(w)
       enddo
       cloud_ssalbs(w) = fa1 * local_cldssalbs(wa1) + fa2 * local_cldssalbs(wa2)
       do l = 0, n_scatmoms
         do k = 1, nmuller
            cloud_scatmoms(1:nmuller,l,w) = fa1 * local_scatmoms(1:nmuller,l,wa1) + fa2 * local_scatmoms(1:nmuller,l,wa2)
         enddo
!        if ( w.eq.1.and.l.lt.50)write(*,*)l,cloud_scatmoms(1:nmuller,l,w) 
       enddo
       cloud_scatmoms(1,0,w) = 1.0d0
!       write(*,*)w,wav(w),fa1,fa2,cloud_ssalbs(w),n_scatmoms,sum(cloud_deltau(1:nlayers,w))
     enddo
   endif

!  Finish 

   return

end subroutine GEMSTOOL_CLD_PROPERTIES

!  End module

end Module GEMSTOOL_CldProperties_m
