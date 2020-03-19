Module GEMSTOOL_AerProperties_m

!  Rob Fix, 10/19/16.
!  New facility for reading and Writing User-defined optical properties
    
!  Use Type structures (GEMSTOOL)

   use GEMSTOOL_Input_Types_m

!  Auxiliary files from the GEMSTOOL_Read_Configfiles_m

   use GEMSTOOL_Read_Configfiles_m, only : Get_FileUnit, &
                                           Read_aerosol_oprop_file, &
                                           Write_aerosol_oprop_file

!  Use Loading routines
!  --------------------

   use GEMSTOOL_aerload_routines_m

!  Use Mie code
!  ------------

!  parameters, sourcecode (single-mode master)

   use RTSMie_parameters_m
   use RTSMie_sourcecode_m

!  Use Tmatrix code
!  ----------------

!  parameters, single-mode master

   use tmat_parameters
   use tmat_master_m

!  10/19/16 Main routine is public, Read/Write User aerosol routines are private

public  :: GEMSTOOL_AER_PROPERTIES
private :: fpk

contains

!        PSD_Index, PSD_pars, n_real, n_imag, FixR1R2, R1, R2, epsnum, epsfac,   & ! INPUT, Mie/Tmat Parameters

subroutine GEMSTOOL_AER_PROPERTIES &
      ( maxlayers, maxwav, maxaermoms,                                             & ! Dimensions
        interpolate_aerosols, do_wavnums, FDAer, FDLay, FDBul, FDeps,              & ! Flags and FD control
        nlayers, nmuller, nwav, wav, height_grid, GEMSTOOL_INPUTS, momsize_cutoff, & ! Inputs
        aerlayerflags, Loading, n_scatmoms, aertau_unscaled, aod_scaling,          & ! OUTPUT, aerosol optical properties
        aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms, aerosol_distchars,       & ! OUTPUT, aerosol optical properties
        fail1, fail2, Message_Loading, Messages_Optical )                            ! Exception handling

!  =============================================================================
!                            AEROSOL REGULAR CREATION
!  =============================================================================

!  aerosol Loading:
!    Loading = optical depth profile
!    Loading Jacobians (Not used here)
!       Cases 1-3 : dloading_Dtau      = derivative of profile w.r.t the total aerosol optical depth at wavelength w0
!       Case 2    : dloading_Dpars(1)  = derivative of profile w.r.t the relaxation parameter (Exponential)
!       Case 3    : dloading_Dpars(1)  = derivative of profile w.r.t the GDF Peak height
!       Case 3    : dloading_Dpars(2)  = derivative of profile w.r.t the GDF Half Width 

!  Generation of Aerosol optical properties
!     1. Call to the Mie/Tmatrix program
!     2. Convert Mie/Tmatrix output (Microsopic) to IOP output (macroscopic) 

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)
  
!  Input variables
!  ---------------

!  External Dimensioning

   integer, INTENT (IN) :: maxlayers, maxwav
   integer, INTENT (IN) :: maxaermoms

!  Flag for interpolation of aerosols

   logical, INTENT(in)  :: interpolate_aerosols

!  Flag for using wavnumber output

   logical, INTENT(in)  :: do_wavnums

!  FD perturbation control

   logical  , INTENT(in)  :: FDAer
   integer  , INTENT(in)  :: FDLay
   integer  , INTENT(in)  :: FDBul
   REAL(fpk), INTENT (IN) :: FDeps

!  Numbers

   integer, INTENT (IN) ::  nlayers
   integer, INTENT (IN) ::  nwav

!  Heights and wavelengths

   REAL    (fpk),    INTENT (IN)   :: wav(maxwav)
   REAL    (fpk),    INTENT (IN)   :: height_grid(0:maxlayers)

!  nmuller = 1 (scalar code), = 6 (vector code)

   integer, INTENT (IN) ::  nmuller

!  Aerosol moment size cutoff. DEFAULT = 0.001

   REAL    (fpk),    INTENT (IN)   :: momsize_cutoff

!  Type Structure inputs

   TYPE(GEMSTOOL_Config_Inputs) :: GEMSTOOL_INPUTS

!  Mie/Tmatrix PSD inputs
!  ======================

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

!  Tmatrix-specific inputs
!  =======================

!  Logical flag for using equal-surface-area sepcification

!      NKMAX.LE.988 is such that NKMAX+2 is the                        
!           number of Gaussian quadrature points used in               
!           integrating over the size distribution for particles
!           MKMAX should be set to -1 for Monodisperse

!      NDGS - parameter controlling the number of division points      
!             in computing integrals over the particle surface.        
!             For compact particles, the recommended value is 2.       
!             For highly aspherical particles larger values (3, 4,...) 
!             may be necessary to obtain convergence.                  
!             The code does not check convergence over this parameter. 
!             Therefore, control comparisons of results obtained with  
!             different NDGS-values are recommended.

!      EPS (Shape_factor) and NP (Spheroid type) - specify the shape of the particles.                
!             For spheroids NP=-1 and EPS is the ratio of the          
!                 horizontal to rotational axes.  EPS is larger than   
!                 1 for oblate spheroids and smaller than 1 for       
!                 prolate spheroids.                                   
!             For cylinders NP=-2 and EPS is the ratio of the          
!                 diameter to the length.                              
!             For Chebyshev particles NP must be positive and 
!                 is the degree of the Chebyshev polynomial, while     
!                 EPS is the deformation parameter                     

!      Accuracy       - accuracy of the computations

!  output
!  ======

!  Loading output

    real(fpk),    dimension ( maxlayers ),  INTENT (OUT)      :: Loading

!  aerosol layer flags

    LOGICAL,      DIMENSION ( maxlayers ), intent(out)  :: AERLAYERFLAGS 

!  AOP output: Number of exapnsion coefficients to be used
!    -- Now wavelength dependent, 10/19/16

    integer, INTENT (OUT)   ::  n_scatmoms(maxwav)

!  Unscaled profiles of optical depth

    real(fpk),    INTENT (OUT)  :: aertau_unscaled(maxlayers)

!  AOD scaling

    real(fpk),    INTENT (OUT)  :: aod_scaling(maxwav)

!  Aod output, final

    REAL(fpk),    DIMENSION( maxlayers, maxwav )      , INTENT (OUT)  :: AEROSOL_DELTAU

!  AOPS

    REAL(fpk),    DIMENSION( maxwav )                 , INTENT (OUT)  :: AEROSOL_SSALBS 
    REAL(fpk),    DIMENSION( 6, 0:maxaermoms, maxwav ), INTENT (OUT)  :: AEROSOL_SCATMOMS

!  Aerosol distribution characteristics, ONLY FOR THE REFERENCE WAVELENGTH

   real(fpk),     DIMENSION( 5, 2 )           , INTENT (OUT)  :: AEROSOL_DISTCHARS

!  Exception handling

   logical,        INTENT (OUT)           :: fail1, fail2
   character*(*),  INTENT (OUT)           :: message_Loading(3)
   character*(*),  INTENT (OUT)           :: Messages_Optical(5)

!  LOCAL VARIABLES
!  @@@@@@@@@@@@@@@

!  Loading derivatives (NOT REQUIRED for OUTPUT, here)

    real(fpk),    dimension ( maxlayers )      :: DLoading_Dtau
    real(fpk),    dimension ( maxlayers, 2 )   :: Dloading_Dpars

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

   INTEGER            :: n_Fangles
   REAL    (fpk)      :: Fangles(max_Mie_angles)

!  Monoradius     - Monodisperse radius size (Microns)

   real    (fpk)     :: Monoradius

!  NPNA - number of equidistant scattering angles (from 0
!             to 180 deg) for which the scattering matrix is
!             calculated.

!  (Tmatrix only). Style flag.
!    * This is checked and re-set (if required) for Monodisperse case

   logical           :: Do_psd_OldStyle

!  Aerosol codes OUTPUT variables
!  ==============================

!  10/19/16. Maximum 2 modes for GEMSTOOL.
!  same settings as in the Type-structures file for "GEMSTOOL_MieTmatUser"

   integer, parameter :: maxAerModes = 2

!  10/19/16.Local dimension for the User aerosols.

   integer, parameter :: max_user_angles = max(max_Mie_angles,NPL1)

!  Bulk optical parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk)    :: MTU_bulk(3,maxAerModes)

!jlee added, extinction coefficients at reference wavelength for each mode and regime

   real(fpk)    :: MTU_bulk_ext_ref(maxAerModes) 

!  Number of Expansion coefficients and Asymmetry parameter

   integer      :: MTU_ncoeffs(maxAerModes)
   real(fpk)    :: MTU_asymm(maxAerModes)

!  Distribution parameters
!    1 = Normalization.  2 = Cross-section. 3 = Volume.  4 = REFF.  5 = VEFF

   real(fpk)    :: MTU_dist(5,maxAerModes)

!  Exception handling

   character*90  :: MTU_messages(3)

!  Aerosol data from Mie subroutine: Expansion coefficients and F-matrix (optional)

   real(fpk)    :: Mie_expcoeffs(6,0:max_Mie_angles,maxAerModes)
   real(fpk)    :: Mie_Fmatrix(4,max_Mie_angles,maxAerModes)

!  Aerosol data from Tmatrix subroutine: Expansion coefficients and F-matrix (optional)

   real(fpk)    :: Tmat_expcoeffs(NPL1,6,maxAerModes)
   real(fpk)    :: Tmat_Fmatrix(MAXNPA,6,maxAerModes)

!  Other aerosol data-related variables
!  ====================================

!  Aerosol data from user input file: Expansion coefficients

   real(fpk)    :: User_expcoeffs(6,0:max_user_angles,maxAerModes)

!  Universal aerosol file output

   real(fpk)    :: MTU_expcoeffs(6,0:max_user_angles,maxAerModes)

!  Aerosol I/O data file names

   character*90  :: AeroInFile,AeroOutFile

!  Local variables
!  ===============

!  fixed Proxies

   integer      :: Tmat_Sphtype, Tmat_ndgs, Tmat_nkmax, nblocks, nweights
   real(fpk)    :: Tmat_accuracy, R1, R2, R1R2_cut, xlimit
   logical      :: FixR1R2, Do_EqSaSphere

!  regime-dependent Local aerosol properties

   integer   :: PSDindex
   real(fpk) :: PSDpars(3)
   real(fpk) :: nreal, nimag
   real(fpk) :: Tmateps

!  AOP output: Reference values
!    * These are the values at reference wavelength w0
!    * They are set for wavelength w0, then used again

   real(fpk)    :: extinction_ref, extinction

!  Help Variables
!  --------------

!  revised 10/19/16 to include many new variables

   character*2  :: cnm
   character*4  :: cwav
   character*9  :: AerChar, CharSrc(3)

   logical      :: RefWvl
   logical      :: MakeAeroFile(maxAerModes)
!   logical      :: ReadAeroFile(maxAerModes)

   integer      :: k, n, nm, nmodes, l, istatus, point_index, w, wc, m, wavmask(maxwav)
   integer      :: OpropSrc(maxAerModes), InFileUnit(maxAerModes), OutFileUnit(maxAerModes)
   real(fpk)    :: Csca_total
   real(fpk)    :: frac(maxAerModes),Csca_frac(maxAerModes),sca_wgt(maxAerModes)
   real(fpk)    :: wavelength, Micron_wavelength

   integer      :: n_scatmoms_w, n_aerwavs, local_nwav, wa, wa1, wa2, wastart, wdum
   real(fpk)    :: fa1, fa2, lam, lamstart, lamfinish, epsfac
   real(fpk)    :: dum1, dum2
   logical      :: trawl, do_bimodal

!  Interpolation arrays
!     UV-Vis : 56  is enough wavelengths for 10 nm intervals over a 250-800  nm range
!     SWIr   : 136 is enough wavelengths for 10 nm intervals over a 750-2100 nm range
!    NOTE    : CAN MAKE THESE ARRAYS ALLOCATABLE ----------THINK ABOUT IT !!!!!!!!!

   real(fpk)    :: aerwav(136),local_aodscaling(136),local_aerssalbs(136), local_scatmoms(6,0:maxaermoms,136)
   integer      :: local_nscatmoms(136)

!mick mod 1/24/2017 - added aerosol file preamble capacity
!  Rob mod 3/7/17 - expanded preamble to include ref wavelength and PSD information (11 lines in all)
!  Rob fix, 3/7/17 - AeroFilePreamble needs to have a mode index
!  Aerosol output file description header
!   integer, parameter :: nPLines = 6
   integer, parameter :: nPLines = 11
   real(fpk)          :: aerwavlo,aerwavhi,wavlo,wav2,wavhi,wavres
   logical            :: DoAeroFilePreamble(maxAerModes)
!   character(len=80)  :: AeroFilePreamble(nPLines)
   character(len=80)  :: AeroFilePreamble(nPLines,maxAerModes)

!  Local control parameters
!  ------------------------

   integer :: dunit
   logical, parameter :: do_iopchecker    = .false.
!   logical, parameter :: do_iopchecker    = .true.

   logical, parameter :: do_aer_Jacobians = .false.

!  New Local control parameters 10/19/16
!  -------------------------------------

   logical, parameter :: jlee_method    = .false.
!   logical, parameter :: jlee_method    = .true.

!  verbose and debug parameters

   logical :: Verbose  = .true.
   logical :: Debug    = .true.

! Tmat_Verbose flag replaces Tmat_Progress

   logical :: Tmat_Verbose  = .true.

!  Initialize output
!  =================

!  12 FD tests, In order of appearance
!  -----------------------------------

   epsfac = 1.0d0 + FDeps
   if (FDBul.ne.0 ) then

!  Mie parameters

     If (FDBul.eq.1) GEMSTOOL_INPUTS%MieTmatUser%nreal(1)          = epsfac * GEMSTOOL_INPUTS%MieTmatUser%nreal(1) 

     If (FDBul.eq.2) GEMSTOOL_INPUTS%MieTmatUser%nimag(1)          = epsfac * GEMSTOOL_INPUTS%MieTmatUser%nimag(1)
     If (FDBul.eq.3) GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1,1)      = epsfac * GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1,1)

     If (FDBul.eq.4) GEMSTOOL_INPUTS%MieTmatUser%PSDpars(2,1)      = epsfac * GEMSTOOL_INPUTS%MieTmatUser%PSDpars(2,1)
     If (FDBul.eq.5) GEMSTOOL_INPUTS%MieTmatUser%nreal(2)          = epsfac * GEMSTOOL_INPUTS%MieTmatUser%nreal(2)
     If (FDBul.eq.6) GEMSTOOL_INPUTS%MieTmatUser%nimag(2)          = epsfac * GEMSTOOL_INPUTS%MieTmatUser%nimag(2)

     If (FDBul.eq.7) GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1,2)      = epsfac * GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1,2)
     If (FDBul.eq.8) GEMSTOOL_INPUTS%MieTmatUser%PSDpars(2,2)      = epsfac * GEMSTOOL_INPUTS%MieTmatUser%PSDpars(2,2)
     If (FDBul.eq.9) GEMSTOOL_INPUTS%MieTmatUser%bimodal_fraction  = epsfac * GEMSTOOL_INPUTS%MieTmatUser%bimodal_fraction

!  Done, aerosol_deltau only

     If (FDBul.eq.10) GEMSTOOL_INPUTS%AerLoad%aertau_input_w0       = epsfac * GEMSTOOL_INPUTS%AerLoad%aertau_input_w0
     If (FDBul.eq.11) GEMSTOOL_INPUTS%AerLoad%gdfloading_peakheight = epsfac * GEMSTOOL_INPUTS%AerLoad%gdfloading_peakheight
     If (FDBul.eq.12) GEMSTOOL_INPUTS%AerLoad%gdfloading_halfwidth  = epsfac * GEMSTOOL_INPUTS%AerLoad%gdfloading_halfwidth

   endif

!  Initialize exception handling

   fail1 = .false.
   fail2 = .false.
   Message_Loading   = ' '
   Messages_Optical  = ' '

!  Initialize Aerosol Loading

   Loading        = 0.0d0
   Dloading_Dtau  = 0.0d0     ! Not required
   Dloading_Dpars = 0.0d0
   aerlayerflags = .false.

!  Initialize optical properties

   extinction_ref    = 0.0d0
   aod_scaling       = 0.0d0
   aerosol_deltau    = 0.0d0
   aerosol_ssalbs    = 0.0d0
   aerosol_scatmoms  = 0.0d0
   aerosol_distchars = 0.0d0
   aertau_unscaled   = 0.0d0

!  Initialize nscatmoms

   n_scatmoms = 50000

!  Set local Mie/Tmatrix variables

   do_monodisperse = .false.             ! Always false
   Do_Fmatrix      = .false.             ! Always false
   DO_Expcoeffs    = .false.             ! this will be set later on
   monoradius      =  0.0d0
   do_psd_oldstyle = .false.

!  Preliminary bookkeeping, necessary before we start
!  --------------------------------------------------

!  Number of modes (1 or 2). Unimoal or Bimodal only.

   nmodes = 1 ; frac(1) = 1.0d0
   do_bimodal = GEMSTOOL_INPUTS%MieTmatUser%do_bimodal
   if ( do_bimodal ) then
      nmodes = 2
      frac(1) =  GEMSTOOL_INPUTS%MieTmatUser%bimodal_fraction
      frac(2) =  1.0d0 - frac(1)
   endif

!  Initialize local aerosol properties (Just a precaution)! New, @@@ RTS 09/11/12

!   n_real = 0.0_fpk ; n_imag   = 0.0_fpk
!   PSDIndex = 0    ; PSD_pars = 0.0_fpk

!  ============================
!  Now form the aerosol Loading
!  ============================

!    Loading = optical depth profile

!    Derivatives : 
!       Cases 1-3 : dloading_Dtau      = derivative of profile w.r.t the total aerosol optical depth at wavelength w0
!       Case 2    : dloading_Dpars(1)  = derivative of profile w.r.t the relaxation parameter (Exponential)
!       Case 3    : dloading_Dpars(1)  = derivative of profile w.r.t the GDF Peak height
!       Case 3    : dloading_Dpars(2)  = derivative of profile w.r.t the GDF Half Width 

!  Case 1: Uniform layer of aerosol
!  ********************************

!  write(*,*)GEMSTOOL_INPUTS%AerLoad%loading_case ; pause 'GRONK'

   if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 1 ) then

      CALL profiles_uniform  &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,  & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_upperboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_lowerboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%aertau_input_w0,            & ! Inputs
            Loading, Dloading_Dtau,                             & ! output
            fail1, message_Loading(1), message_Loading(2) )       ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'Uniform aerosol Loading failed'

!  Case 2: Exponential decay profile
!  *********************************

   else if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 2 ) then

      CALL profiles_expone &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,  & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_upperboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_lowerboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%exploading_relaxation,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%aertau_input_w0,            & ! Inputs
            Loading, Dloading_Dtau, Dloading_Dpars(:,1),        & ! output
            fail1, message_Loading(1), message_Loading(2) )       ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'Exponential aerosol Loading failed'

!  Case 3: GDF (quasi-Gaussian) profile
!  ************************************

   else if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 3 ) then

      CALL profiles_gdfone &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,  & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_upperboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%gdfloading_peakheight,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_lowerboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%gdfloading_halfwidth,       & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%aertau_input_w0,            & ! Inputs
            Loading, Dloading_Dtau, Dloading_Dpars(:,1), Dloading_Dpars(:,2),     & ! output
            fail1, message_Loading(1), message_Loading(2) )                         ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'GDF aerosol Loading failed'

!  Case 4: User-defined profile
!  ****************************

    else if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 4 ) then
        
        print *, 'User-defined AEROSOL PROFILES !!! '
        print *, nlayers
        
        open(11, file='aod_loading_user.input', status = 'old')

        do n = 1, nlayers
            read(11,*) dum1, dum2
            loading(n) = dum2
            ! print *, n, loading(n)
        enddo
        close(11)
   endif


! mchoi. Overwrite AOD loading profiles to aod_loading_user.input file
   open(11, file='aod_loading_user.input', status = 'replace')
   do n = 1, nlayers
      write(11,900) n-1, loading(n)
   enddo
   close(11)
   ! pause
900 format(i7,1p1e20.10e3)


!  Return if failure at this stage

   if ( fail1 ) return



!  Assign Unscaled loadings

   do n = 1, nlayers
      aerlayerflags(n) =  ( Loading(n) .ne. 0.0d0  )
      !print*, n, loading(n), aerlayerflags(n)
      aertau_unscaled(n) = loading(n)
      if ( aerlayerflags(n) .and. FDAer .and. FDLay.eq.n ) aertau_unscaled(n) = aertau_unscaled(n) * (1.0d0 + FDeps )
!      write(*,*)n,(height_grid(n-1)+heightgrid(n))*0.5d0,aerlayerflags(n),loading(n)
   enddo
!  write(*,*)'TOTAL',sum(loading(1:nlayers)) ,GEMSTOOL_INPUTS%AerLoad%aertau_input_w0

!   do n = 1, nlayers
!      if (aerlayerflags(n) ) write(*,*)n,aertau_unscaled(n)
!   enddo
!   pause' aerosol loading, FD check'

!  Interpolation setup
!  ===================

   if ( interpolate_aerosols ) then
!mick mod 1/24/2017 - provide user feedback in interpolated aerosol case
!                     (here and below)
      !write(*,*)
      write(*,'(a)') '    Interpolating aerosols ...'
      trawl = .true. ; wa = 0
      if ( .not. do_wavnums ) then
         !Wavelength set
         lamstart = 200.0d0 ; lam = lamstart        !  Smallest wavelength is 200 nm (UVN application)
         do while (trawl)
            lam = lam + 20.0d0
            if ( lam .ge. wav(1) ) then
               wa = wa + 1 ; aerwav(wa) = lam - 20.0d0
               if ( lam .gt. wav(nwav) ) then
                  wa = wa + 1 ; aerwav(wa) = lam ; trawl = .false.
               endif
            endif
         enddo
          n_aerwavs = wa ; local_nwav = n_aerwavs
      else
         !Wavenumber set
         lamfinish = 2100.0d0 ; lam = lamfinish     !  Largest wavelength is 2100 nm (NSW application)
         do while (trawl)
            lam = lam - 20.0d0
            if ( lam .le. 1.0d+07/wav(1) ) then
               wa = wa + 1 ; aerwav(wa) = lam + 20.0d0
               if ( lam .lt. 1.0d+07/wav(nwav) ) then
                  wa = wa + 1 ; aerwav(wa) = lam ; trawl = .false.
               endif
            endif
         enddo
          n_aerwavs = wa ; local_nwav = n_aerwavs
      endif
      write(*,'(a,i2)')   '    Number of aerosol optical property wavelengths: ',n_aerwavs
      if ( .not. do_wavnums ) then
         write(*,'(a,f8.3)') '    Lo wavelength (nm): ',aerwav(1)
         write(*,'(a,f8.3)') '    Hi wavelength (nm): ',aerwav(n_aerwavs)
      else
         write(*,'(a,f8.3)') '    Lo wavelength (nm): ',aerwav(n_aerwavs)
         write(*,'(a,f8.3)') '    Hi wavelength (nm): ',aerwav(1)
      endif
      write(*,'(a,f8.3)') '    Resolution    (nm): ',20.0d0
   else
      local_nwav = nwav
   endif

!  Check to see if reference wavelength is one of the set. Initialize Mask.

   point_index = 0
   if ( .not. do_wavnums ) then
      if ( interpolate_aerosols ) then
         do w = 1, n_aerwavs
            wavmask(w) = w ; if ( aerwav(w) .eq. GEMSTOOL_INPUTS%AerLoad%reference_w0 ) point_index = w
         enddo
      else
         do w = 1, nwav
            wavmask(w) = w ; if ( wav(w)    .eq. GEMSTOOL_INPUTS%AerLoad%reference_w0 ) point_index = w
         enddo
      endif
   else
      if ( interpolate_aerosols ) then
         do w = 1, n_aerwavs
            wavmask(w) = w ; if ( aerwav(w) .eq. GEMSTOOL_INPUTS%AerLoad%reference_w0 ) point_index = w
         enddo
      else
         do w = 1, nwav
            wavmask(w) = w
         enddo
      endif
   endif

!  Mask to use IF reference wavelength IS one of list of wavelengths. [UVN only]

!   if ( .not.do_Wavnums ) then
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
!   endif

!  ==================
!  OPTICAL PROPERTIES
!  ==================

!  Make aerosol file flags. These will be set if you are doing Mie or Tmatrix

    MakeAeroFile = .true.
    if ( GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then
       MakeAeroFile(1:nmodes) = .false.
    endif
!mick mod 1/24/2017 - added aerosol file preamble capacity
    DoAeroFilePreamble = .false.

!  Optical sources

   OpropSrc = 0 ; CharSrc = ' '
   if ( GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols ) then
      OpropSrc(1) = 1 ; if ( GEMSTOOL_INPUTS%MieTmatUser%do_bimodal) OpropSrc(2) = 1 ; CharSrc(1) = '(MieScat)'
   else if ( GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols ) then
      OpropSrc(1) = 2 ; if ( GEMSTOOL_INPUTS%MieTmatUser%do_bimodal) OpropSrc(2) = 2 ; CharSrc(2) = '(Tmatrix)'
   else
      OpropSrc(1) = 3 ; if ( GEMSTOOL_INPUTS%MieTmatUser%do_bimodal) OpropSrc(2) = 3 ; CharSrc(3) = '(UserDef)'
   endif
 
!  First, open files for READING user-defined aerosol optical property input if requested
!  (note: these files closed @ ~line 1340)
!mick mod 1/24/2017 - added filename "Debug" output

   if ( GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then
      do nm = 1, nmodes
         call Get_FileUnit(InFileUnit(nm))
         AeroInFile='UserAeroFiles/in/' // trim(GEMSTOOL_INPUTS%MieTmatUser%User_Aerofile_names(nm))
         open(unit=InFileUnit(nm), file=trim(AeroInFile),iostat=istatus,status='old')
         if ( Debug ) then
           if (nm .eq. 1) write(*,*)
           write(*,*) '   Using ' // trim(AeroInFile)
         endif
         if ( istatus > 0 ) then
            Messages_Optical(1) = 'Failure: Aerosol optical property input file'
            Messages_Optical(2) = 'File   : ' // trim(AeroInFile)
            Messages_Optical(3) = 'Problem: Error opening file'
            fail2 = .true. ; return
         endif

!mick mod 1/24/2017 - added aerosol file preamble capacity
         DoAeroFilePreamble(nm) = .true.
      enddo
   endif

!  Second, open files for WRITING Mie or T-matrix aerosol optical property output
!   - This is independent of the above read-user code

   if ( .not. GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then
     do nm = 1, nmodes
       if ( MakeAeroFile(nm) ) then
         call Get_FileUnit(OutFileUnit(nm)) ; write(cnm,'(i2.2)') nm
         if ( OpropSrc(nm) .eq. 1 ) then
            AeroOutFile='UserAeroFiles/out/' // 'TestProps_mode' // cnm // '_Mie_prop.dat'
         elseif ( OpropSrc(nm) .eq. 2 ) then
            AeroOutFile='UserAeroFiles/out/' // 'TestProps_mode' // cnm // '_Tmat_prop.dat'
         else
            AeroOutFile='UserAeroFiles/out/' // 'TestProps_mode' // cnm // '_User_prop.dat'
         endif
         open(unit=OutFileUnit(nm), file=trim(AeroOutFile), status='replace')

!mick mod 1/24/2017 - added aerosol file preamble capacity
!  Define aerosol optical property file preamble

         if ( .not. do_wavnums ) then
            !already in nm
            wavlo = wav(1)
            wav2  = wav(2)
            wavhi = wav(nwav)

            aerwavlo = aerwav(1)
            aerwavhi = aerwav(n_aerwavs)
         else
            !convert from wavenumbers in cm^-1 to wavelengths in nm
            wavlo = 1.0d+07/wav(nwav)   
            wav2  = 1.0d+07/wav(nwav-1)
            wavhi = 1.0d+07/wav(1)

            aerwavlo = aerwav(n_aerwavs)
            aerwavhi = aerwav(1)
         endif
         wavres = dabs(wav2-wavlo)

         DoAeroFilePreamble(nm) = .true.

!  Original 6-lines Preamble.
!         write(AeroFilePreamble(1),'(6x,a)')    'Inputs used to create aerosol file (all wvls in nm):'
!         write(AeroFilePreamble(2),'(6x,a,l1)') 'Aerosols interpolated: ',interpolate_aerosols
!         if ( interpolate_aerosols ) then
!            write(AeroFilePreamble(3),'(6x,2(a,f8.3))') 'Band Lo Wvl  : ',wavlo ,'   Aero Lo Wvl  : ',aerwavlo
!            write(AeroFilePreamble(4),'(6x,2(a,f8.3))') 'Band Hi Wvl  : ',wavhi ,'   Aero Hi Wvl  : ',aerwavhi
!            write(AeroFilePreamble(5),'(6x,2(a,f8.3))') 'Band Res     : ',wavres,'   Aero Res     : ',20.0d0
!            write(AeroFilePreamble(6),'(6x,2(a,i8  ))') 'Band Num Wvls: ',nwav  ,'   Aero Num Wvls: ',n_aerwavs
!         else 
!            write(AeroFilePreamble(3),'(6x,a,f8.3)')    'Band Lo Wvl  : ',wavlo
!            write(AeroFilePreamble(4),'(6x,a,f8.3)')    'Band Hi Wvl  : ',wavhi
!            write(AeroFilePreamble(5),'(6x,a,f8.3)')    'Band Res     : ',wavres
!            write(AeroFilePreamble(6),'(6x,a,i8  )')    'Band Num Wvls: ',nwav
!         endif

!  New 11-line preamble

         write(AeroFilePreamble(1,nm),'(6x,a)')    'Inputs used to create aerosol file (all wvls in nm):'
         write(AeroFilePreamble(2,nm),'(6x,a,l1)') 'Aerosols interpolated: ',interpolate_aerosols
         write(AeroFilePreamble(3,nm),'(6x,a,f8.3)') 'Aerosol RefWavelength: ',GEMSTOOL_INPUTS%AerLoad%reference_w0
         if ( interpolate_aerosols ) then
            write(AeroFilePreamble(4,nm),'(6x,2(a,f8.3))') 'Band Lo Wvl  : ',wavlo ,'   Aero Lo Wvl  : ',aerwavlo
            write(AeroFilePreamble(5,nm),'(6x,2(a,f8.3))') 'Band Hi Wvl  : ',wavhi ,'   Aero Hi Wvl  : ',aerwavhi
            write(AeroFilePreamble(6,nm),'(6x,2(a,f8.3))') 'Band Res     : ',wavres,'   Aero Res     : ',20.0d0
            write(AeroFilePreamble(7,nm),'(6x,2(a,i8  ))') 'Band Num Wvls: ',nwav  ,'   Aero Num Wvls: ',n_aerwavs
         else 
            write(AeroFilePreamble(4,nm),'(6x,a,f8.3)')    'Band Lo Wvl  : ',wavlo
            write(AeroFilePreamble(5,nm),'(6x,a,f8.3)')    'Band Hi Wvl  : ',wavhi
            write(AeroFilePreamble(6,nm),'(6x,a,f8.3)')    'Band Res     : ',wavres
            write(AeroFilePreamble(7,nm),'(6x,a,i8  )')    'Band Num Wvls: ',nwav
         endif
         PSDIndex = GEMSTOOL_INPUTS%MieTmatUser%PSDIndex(nm)
         PSDpars  = GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1:3,nm)
         write(AeroFilePreamble(8,nm), '(6x,a,i9)')      'PSD Index      : ',PSDIndex
         write(AeroFilePreamble(9,nm), '(6x,a,f9.5)')    'PSD parameter 1: ',PSDpars(1)
         write(AeroFilePreamble(10,nm),'(6x,a,f9.5)')    'PSD parameter 2: ',PSDpars(2)
         write(AeroFilePreamble(11,nm),'(6x,a,f9.5)')    'PSD parameter 3: ',PSDpars(3)

       endif
     enddo
   endif

!  Prepare the reference-wavelength Mie/Tmatrix/User Inputs
!  ========================================================

   if ( point_index .eq. 0 ) then

!  Only require extinction coefficient if flagged
!  Set the local Mie program inputs (bulk properties only)

      Do_Expcoeffs = .false.
      RefWvl       = .true.

!  Dummy index for writing

      wdum = 0

!  Reference wavelength

      wavelength        = GEMSTOOL_INPUTS%AerLoad%reference_w0
      Micron_wavelength = wavelength/1000.0d0

!  Progress. 2/2/15 Format changed from f8.3 to f12.5

      if ( Verbose .or. Debug ) then
         if ( Debug ) write(*,*)
         write(*,'(4x,a,5x,f12.5)') 'Aerosol calculation - Doing reference wavelength   : ',wavelength
      endif

!  Reference extinction

      extinction_ref = 0.0_fpk

!  Mode loop

      do nm = 1, nmodes

!  Aerosol type

         AerChar = CharSrc(OPropSrc(nm))

!  Mie/Tmat Proxies, mode dependent

         if ( OPropSrc(nm) .eq. 1 .or.  OPropSrc(nm) .eq. 2 ) then
            PSDIndex = GEMSTOOL_INPUTS%MieTmatUser%PSDIndex(nm)
            PSDpars  = GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1:3,nm)
            R1       = GEMSTOOL_INPUTS%MieTmatUser%R1(nm)
            R2       = GEMSTOOL_INPUTS%MieTmatUser%R2(nm)
            FixR1R2  = GEMSTOOL_INPUTS%MieTmatUser%FixR1R2(nm)
            nreal    = GEMSTOOL_INPUTS%MieTmatUser%nreal(nm)
            nimag    = GEMSTOOL_INPUTS%MieTmatUser%nimag(nm)
         endif

         if ( OPropSrc(nm) .eq. 1 ) then
            nblocks   = GEMSTOOL_INPUTS%MieTmatUser%nblocks(nm)
            nweights  = GEMSTOOL_INPUTS%MieTmatUser%nweights(nm)
            xlimit    = GEMSTOOL_INPUTS%MieTmatUser%xparticle_limit
            R1R2_cut  = GEMSTOOL_INPUTS%MieTmatUser%R1R2_cutoff(nm)
         else if ( OPropSrc(nm) .eq. 2 ) then
            Do_EqSaSphere = GEMSTOOL_INPUTS%MieTmatUser%Do_EqSaSphere
            Tmat_Sphtype  = GEMSTOOL_INPUTS%MieTmatUser%Tmat_Sphtype
            Tmateps       = GEMSTOOL_INPUTS%MieTmatUser%Tmat_eps(nm)
            Tmat_nkmax    = GEMSTOOL_INPUTS%MieTmatUser%Tmat_nkmax(nm)
            Tmat_ndgs     = GEMSTOOL_INPUTS%MieTmatUser%Tmat_ndgs(nm)
            Tmat_accuracy = GEMSTOOL_INPUTS%MieTmatUser%Tmat_accuracy  
         endif
         
!  Progress

         if ( Verbose .or. Debug ) then
            write(*,'(5x,a,1x,i3,1x,A9)') '-- Doing mode # ',nm, charSrc(OpropSrc(nm))
         endif

!  SINGLE CALL. [Tmat_Verbose flag introduced 6/6/14]

         if ( OpropSrc(nm) .eq. 1 ) then  !Mie
            call RTSMie_main  &
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,                                     & ! Input top level flags
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, nblocks, nweights,             & ! Inputs (PSD stuff)
                xlimit, R1R2_cut, n_Fangles, Fangles, Micron_wavelength, nreal, nimag,         & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm), Mie_expcoeffs(:,:,nm),         & ! Outputs
                Mie_Fmatrix(:,:,nm), MTU_dist(:,nm),                                           & ! Outputs
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exeception handling
         elseif ( OpropSrc(nm) .eq. 2 ) then !Tmat
            call tmat_master ( Tmat_Verbose, &
                Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_EqSaSphere, Do_psd_OldStyle,     & ! Input top level flags
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, Tmat_Sphtype, Tmat_nkmax,      & ! Inputs (PSD stuff)
                n_Fangles, Tmat_ndgs, Tmateps, Tmat_accuracy, Micron_wavelength, nreal, nimag, & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm), Tmat_expcoeffs(:,:,nm),        & ! Outputs
                Tmat_Fmatrix(:,:,nm), MTU_dist(:,nm),                                          & ! Outputs
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
         else  !User
           call Read_aerosol_oprop_file  &
              ( max_user_angles, InFileUnit(nm), Micron_wavelength, RefWvl,                             & ! Inputs
                maxAerModes, nm, nPLines, DoAeroFilePreamble,                                           & ! I/O
                MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), User_expcoeffs(1,0,nm), & ! Outputs
                fail2, MTU_messages )                                                                     ! Exception handling
            if ( fail2 ) MTU_messages(2) = 'File   : ' // trim(GEMSTOOL_INPUTS%MieTmatUser%User_Aerofile_names(nm))
         endif

!  Exception handling on everything

         if ( Fail2 ) then
            do m = 1, 3
               Messages_Optical(m) = adjustl(trim(MTU_messages(m)))
            enddo
            if ( OpropSrc(nm) .eq. 1 .or. OpropSrc(nm) .eq. 2 ) then
               Messages_Optical(4) = 'First call to regular '//AerChar//' program, reference wavelength'
            else
               Messages_Optical(4) = 'First read from regular '//AerChar//' file, reference wavelength'
            endif
            Messages_Optical(5) = 'Failure for reference wavelength calculation'
            return
         endif

!  Coefficients Must be copied. Note special care over Tmat output.
!   This is the reference wavelength calculation, do not need coefficients

         MTU_expcoeffs(1,0,nm) = 0.0_fpk 
!         if ( OpropSrc(nm) .eq. 1 ) then
!            do k = 1, 6
!               MTU_expcoeffs(k,0:MTU_ncoeffs(nm),nm) = Mie_expcoeffs(k,0:MTU_ncoeffs(nm),nm)
!            enddo
!         else if ( OpropSrc(nm) .eq. 2 ) then
!            do k = 1, 6
!               MTU_expcoeffs(k,0:MTU_ncoeffs(nm)-1,nm) = Tmat_expcoeffs(1:MTU_ncoeffs(nm),k,nm)
!            enddo
!            MTU_ncoeffs(nm) = MTU_ncoeffs(nm)-1
!         else
!            do k = 1, 6
!               MTU_expcoeffs(k,0:MTU_ncoeffs(nm),nm) = User_expcoeffs(k,0:MTU_ncoeffs(nm),nm)
!            enddo
!         endif

!  Write aerosol optical properties to file (optional)
!     - Coefficients Do not need to be copied here

         if ( MakeAeroFile(nm) ) then
            call Write_aerosol_oprop_file  &
                 ( max_user_angles, OutFileUnit(nm), w, Micron_wavelength, RefWvl,                       & ! Inputs
                   maxAerModes, nm, nPLines, DoAeroFilePreamble, AeroFilePreamble,                       & ! I/O
                   MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), MTU_expcoeffs(1,0,nm) ) ! Inputs
         endif

!  Set the reference quantity

         !Composite extinction @ REF wvl produced by mode extinctions from
         !  the Tmat or Mie code weighted by input mode fractions.
         !  --> Case: ref wvl is SEPARATE FROM the current spectral band 

         MTU_bulk_ext_ref(nm) = MTU_bulk(1,nm) ! jlee added, extinction coefficients at reference wavelength
         extinction_ref = extinction_ref + frac(nm)*MTU_bulk(1,nm) 

!  Set the Distribution

         AEROSOL_DISTCHARS(1:5,nm) = MTU_dist(1:5,nm)

!  End mode loop

      enddo

!  Progress

      if ( Verbose .or. Debug ) then
         write(*,'(4x,a,5x,f12.5)') &
            'Aerosol calculation - Finished reference wavelength: ',wavelength
      endif

!  End reference-wavelength calculation

   endif

!  Prepare General (all-wavelength, all-wavenumber) Mie/Tmat/User inputs
!  =====================================================================

!  Set the local inputs (general)

   Do_Expcoeffs = .true.
   RefWvl       = .false.

!  wavelength loop. First wavelength will be the reference, if in  list.
!  Wavnumber  loop. 

   do wc = 1, local_nwav

!  Wavelengths [nm]

      w = wavmask(wc)

!mick note:
!      aerwav(w) - comes from interpolation setup code above
!      wav(w)    - comes from spectral config file (either "GEMSTOOL_Lambdas.cfg" or "GEMSTOOL_Wavenums.cfg")

      if ( interpolate_aerosols ) then
         wavelength = aerwav(w)
      else
        if ( do_wavnums      ) wavelength = 1.0d+07/wav(w) !convert from cm^-1 to nm
        if ( .not.do_wavnums ) wavelength = wav(w)         !in nm
      endif
      Micron_wavelength = wavelength/1000.0d0 !in um

!  Progress.  2/2/15 Format changed from f8.3 to f12.5

      if ( Verbose.or.Debug ) then
         if ( Debug ) write(*,*)
         write(*,'(4x,a,1x,i3,1x,f12.5)') 'Aerosol calculation - Doing point, wavelength      : ', wc, wavelength
      endif

!  Mode loop

      do nm = 1, nmodes

!  Aerosol type

         AerChar = CharSrc(OPropSrc(nm)) ; write(cwav,'(I4)')wc ; write(cnm,'(I2)') nm

!  Mie/Tmat Proxies, mode dependent

         if ( OPropSrc(nm) .eq. 1 .or.  OPropSrc(nm) .eq. 2 ) then
            PSDIndex = GEMSTOOL_INPUTS%MieTmatUser%PSDIndex(nm)
            PSDpars  = GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1:3,nm)
            R1       = GEMSTOOL_INPUTS%MieTmatUser%R1(nm)
            R2       = GEMSTOOL_INPUTS%MieTmatUser%R2(nm)
            FixR1R2  = GEMSTOOL_INPUTS%MieTmatUser%FixR1R2(nm)
            nreal    = GEMSTOOL_INPUTS%MieTmatUser%nreal(nm)
            nimag    = GEMSTOOL_INPUTS%MieTmatUser%nimag(nm)
         endif
         if ( OPropSrc(nm) .eq. 1 ) then
            nblocks   = GEMSTOOL_INPUTS%MieTmatUser%nblocks(nm)
            nweights  = GEMSTOOL_INPUTS%MieTmatUser%nweights(nm)
            xlimit    = GEMSTOOL_INPUTS%MieTmatUser%xparticle_limit
            R1R2_cut  = GEMSTOOL_INPUTS%MieTmatUser%R1R2_cutoff(nm)
         else if ( OPropSrc(nm) .eq. 2 ) then
!mick fix 1/24/2017 - define "Do_EqSaSphere"
            Do_EqSaSphere = GEMSTOOL_INPUTS%MieTmatUser%Do_EqSaSphere
            Tmat_Sphtype  = GEMSTOOL_INPUTS%MieTmatUser%Tmat_Sphtype
            Tmateps       = GEMSTOOL_INPUTS%MieTmatUser%Tmat_eps(nm)
            Tmat_nkmax    = GEMSTOOL_INPUTS%MieTmatUser%Tmat_nkmax(nm)
            Tmat_ndgs     = GEMSTOOL_INPUTS%MieTmatUser%Tmat_ndgs(nm)
            Tmat_accuracy = GEMSTOOL_INPUTS%MieTmatUser%Tmat_accuracy  
         endif

!  Progress

         if ( Verbose .or. Debug ) then
            write(*,'(5x,a,1x,i3,1x,A9)') '-- Doing mode # ',nm, charSrc(OpropSrc(nm))
         endif

!  SINGLE CALL. [Tmat_Verbose flag introduced 6/6/14]

         if ( OpropSrc(nm) .eq. 1 ) then  !Mie
            call RTSMie_main  &
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,                                     & ! Input top level flags
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, nblocks, nweights,             & ! Inputs (PSD stuff)
                xlimit, R1R2_cut, n_Fangles, Fangles, Micron_wavelength, nreal, nimag,         & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm), Mie_expcoeffs(:,:,nm),         & ! Outputs
                Mie_Fmatrix(:,:,nm), MTU_dist(:,nm),                                           & ! Outputs
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exeception handling
         elseif ( OpropSrc(nm) .eq. 2 ) then !Tmat
            call tmat_master ( Tmat_Verbose, &
                Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_EqSaSphere, Do_psd_OldStyle,     & ! Input top level flags
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, Tmat_Sphtype, Tmat_nkmax,      & ! Inputs (PSD stuff)
                n_Fangles, Tmat_ndgs, Tmateps, Tmat_accuracy, Micron_wavelength, nreal, nimag, & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm), Tmat_expcoeffs(:,:,nm),        & ! Outputs
                Tmat_Fmatrix(:,:,nm), MTU_dist(:,nm),                                          & ! Outputs
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
         else  !User
           call Read_aerosol_oprop_file  &
              ( max_user_angles, InFileUnit(nm), Micron_wavelength, RefWvl,                             & ! Inputs
                maxAerModes, nm, nPLines, DoAeroFilePreamble,                                           & ! I/O
                MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), User_expcoeffs(1,0,nm), & ! Outputs
                fail2, MTU_messages )                                                                     ! Exception handling
            if ( fail2 ) MTU_messages(2) = 'File   : ' // trim(GEMSTOOL_INPUTS%MieTmatUser%User_Aerofile_names(nm))
         endif

!  Exception handling on everything

         if ( Fail2 ) then
            do m = 1, 3
               Messages_Optical(m) = adjustl(trim(MTU_messages(m)))
            enddo
            if ( OpropSrc(nm) .eq. 1 .or. OpropSrc(nm) .eq. 2 ) then
               Messages_Optical(4) = 'Second call to regular '//AerChar//' program, mode # '//cnm//', wavelength # '//cwav
            else
               Messages_Optical(4) = 'Second read from regular '//AerChar//' file, mode # '//cnm//', wavelength # '//cwav
            endif
            Messages_Optical(5) = 'Failure for  mode # '//cnm
            return
         endif

!  Coefficients Must be copied. Note special care over Tmat output.

         if ( OpropSrc(nm) .eq. 1 ) then
            do k = 1, 6
               MTU_expcoeffs(k,0:MTU_ncoeffs(nm),nm) = Mie_expcoeffs(k,0:MTU_ncoeffs(nm),nm)
            enddo
         else if ( OpropSrc(nm) .eq. 2 ) then
            do k = 1, 6
               MTU_expcoeffs(k,0:MTU_ncoeffs(nm)-1,nm) = Tmat_expcoeffs(1:MTU_ncoeffs(nm),k,nm)
            enddo
            MTU_ncoeffs(nm) = MTU_ncoeffs(nm)-1
         else
            do k = 1, 6
               MTU_expcoeffs(k,0:MTU_ncoeffs(nm),nm) = User_expcoeffs(k,0:MTU_ncoeffs(nm),nm)
            enddo
         endif

!  Write aerosol optical properties to file (optional)
 
         if ( MakeAeroFile(nm) ) then
            call Write_aerosol_oprop_file  &
                 ( max_user_angles, OutFileUnit(nm), w, Micron_wavelength, RefWvl,                       & ! Inputs
                   maxAerModes, nm, nPLines, DoAeroFilePreamble, AeroFilePreamble,                       & ! I/O
                   MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), MTU_expcoeffs(1,0,nm) ) ! Inputs
         endif

!  Scattering fraction

         if ( .not. jlee_method ) then
            Csca_frac(nm) = frac(nm)*MTU_bulk(2,nm)
         endif

!  Debug
!         write(445,*)wc,w,nm,Csca_frac(nm),MTU_bulk(1:3,nm)

!  End mode loop

      enddo

!  Scattering weights

      if ( .not. jlee_method ) then
         Csca_total  = sum( Csca_frac(1:nmodes) )           
         sca_wgt(1:nmodes) = Csca_frac(1:nmodes) / Csca_total 
      endif

!  Set the reference quantities, if reference wavelength is in the list.
!    Values for the first (masked) wavelength
!  Composite extinction @ REF wvl produced by mode extinctions from
!    the Tmat or Mie code weighted by input mode fractions.
!    --> Case: ref wvl IS A MEMBER of the current spectral band

      if ( point_index.ne.0 .and. wc.eq.1 ) then
         extinction_ref = 0.0_fpk
         do nm = 1, nmodes
            MTU_bulk_ext_ref(nm) = MTU_bulk(1,nm)
            extinction_ref = extinction_ref + frac(nm)*MTU_bulk(1,nm)
         enddo
      endif

! jlee added, RSpurr option returned
      if ( jlee_method ) then
         do nm = 1, nmodes
            Csca_frac(nm) = frac(nm)*MTU_bulk(2,nm)/MTU_bulk_ext_ref(nm)
         enddo
         Csca_total        = sum( Csca_frac(1:nmodes) )            ! sum of normalized scattering AODs at calc. wavelen.
         sca_wgt(1:nmodes) = Csca_frac(1:nmodes) / Csca_total ! scattering AOD fraction of each mode at calc. wavelen
      endif
! end jlee added

!  NOW set the optical property output
!  ===================================

!  Composite extinction @ CURRENT wvl produced by mode extinctions from
!     the Tmat or Mie code weighted by input mode fractions
!     * jlee Method  --> sum of normalized extinction AODs at calculation wavelength

!  Scaling factor. All wavelengths, Using the composite Tmat/Mie extinctions @ CURRENT wvl and REF wvl,
!       define a scale factor to make the resulting aero tau @ CURRENT wvl  compatible with the input aero tau @ REF wvl

!  1a. For the interpolation case, Create output on local arrays
!  -------------------------------------------------------------
 
      if ( interpolate_aerosols ) then

!  Extinction and scaling

         extinction = 0.0_fpk
         if ( point_index .ne. 0 .and. wc.eq.1 ) then
            extinction = Extinction_ref
            local_aodscaling(w) = 1.0_fpk
            if ( jlee_method ) local_aodscaling(w) =  extinction
         else
            if ( jlee_method ) then
               do nm = 1, nmodes
                  extinction = extinction + frac(nm)*MTU_bulk(1,nm)/MTU_bulk_ext_ref(nm) 
               enddo
               local_aodscaling(w) = extinction
            else
               extinction = dot_product(frac(1:nmodes),MTU_bulk(1,1:nmodes))
               local_aodscaling(w) = extinction / extinction_ref
            endif
          endif
 
!  SSAs and Expansion coefficients

         local_aerssalbs(w) = Csca_total / extinction
         l = 0 ; local_scatmoms(1,0,w) = 1.0d0
         do while ( local_scatmoms(1,l,w).gt.momsize_cutoff .and. l.lt.maxaermoms )
            l = l + 1 ; local_scatmoms(1,l,w) = dot_product(sca_wgt(1:nmodes),MTU_expcoeffs(1,l,1:nmodes))
         enddo
         n_scatmoms_w = l
         do l = 0, n_scatmoms_w
            do k = 2, nmuller
               local_scatmoms(k,l,w) = dot_product(sca_wgt(1:nmodes),MTU_expcoeffs(k,l,1:nmodes))
            enddo
         enddo
         local_nscatmoms(w) = n_scatmoms_w

!  Debug aerosol optical properties. VERSION TWO only

!         if ( do_iopchecker ) then
!            dunit = 999 ; if ( GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) dunit = 998
!            write(dunit,'(2i4,f10.4,1p2e20.10,i5)')wc,w,wavelength,&
!                                     local_aodscaling(w),local_aerssalbs(w),local_nscatmoms(w) 
!            do l = 0, local_nscatmoms(w) 
!               write(dunit,'(2i5,1p6e20.10)')w,l,(local_scatmoms(k,l,w),k=1,6)
!            enddo
!         endif

!  End "interpolate aerosols" if block

      endif

!  1b. For the Monochromatic case, Create output on final arrays
!  -------------------------------------------------------------

      if ( .not. interpolate_aerosols ) then

!  Extinction and scaling

         extinction = 0.0_fpk
         if ( jlee_method ) then
            do nm = 1, nmodes
               extinction = extinction + frac(nm)*MTU_bulk(1,nm)/MTU_bulk_ext_ref(nm) 
            enddo
            aod_scaling(w) = extinction
         else
            do nm = 1, nmodes
               extinction = extinction + frac(nm)*MTU_bulk(1,nm)
            enddo
            aod_scaling(w) = extinction / extinction_ref
         endif

!  SSAs and assymetry parameters (Not doing the latter)
!   Rob Fix 8/25/14. Total SSA is not scatter-weighted sum of individual mode SSAs

         aerosol_ssalbs(w) = Csca_total / extinction
!         aerosol_Asymms(w) = 0.0_fpk
!         do nm = 1, nmodes
!            aerosol_Asymms(w)   = aerosol_Asymms(w)  + sca_wgt(nm)*MTU_Asymm(nm)
!         enddo

         l = 0 ;  aerosol_scatmoms(1,0,w) = 1.0_fpk
         do while ( aerosol_scatmoms(1,l,w).gt.momsize_cutoff .and. l.lt.maxaermoms )
            l = l + 1
            do nm = 1, nmodes
               aerosol_scatmoms(1,l,w) = aerosol_scatmoms(1,l,w) + sca_wgt(nm)*MTU_expcoeffs(1,l,nm)
            enddo
         enddo
         n_scatmoms_w = l
         do l = 0, n_scatmoms_w
            do k = 2, nmuller
               do nm = 1, nmodes
                  aerosol_scatmoms(k,l,w) = aerosol_scatmoms(k,l,w) + sca_wgt(nm)*MTU_expcoeffs(k,l,nm)
               enddo
            enddo
         enddo
         n_scatmoms(w) = n_scatmoms_w

!  Apply scalings to loadings (i.e. define the layer aero tau @ CURRENT wvl so it's
!                              compatible with the layer aero tau @ REF wvl using the
!                              previously defined scale factors)

         do n = 1, nlayers
            aerosol_deltau(n,w) = aertau_unscaled(n) * aod_scaling(w)
         enddo

!  Debug aerosol optical properties. VERSION TWO only

         if ( do_iopchecker ) then
            do n = 1, nlayers
              if (aerlayerflags(N).and.n.eq.nlayers ) then
                 write(999,'(i4,1p6e20.10)')n,aerosol_deltau(n,w),aerosol_ssalbs(w)
                 do l = 0, n_scatmoms_w
                    write(999,'(2i5,1p6e20.10)')n,l,(aerosol_scatmoms(k,l,w),k=1,1)
                 enddo
              endif
            enddo
         endif

!  End monochromatic ("not interpolate aerosols") if block
!    (replaces continuation point 677)

      endif

!  continuation point for avoiding the exact monochromatic solution
!677   continue

!  Progress

      if ( Verbose.or.Debug ) then
         write(*,'(4x,a,1x,i3,1x,f12.5,2(1x,i3))') &
            'Aerosol calculation - Finished point, wavelength   : ',wc, wavelength
      endif

!  End wavelength/wavenumber loop

   enddo

!  1c. carry out interpolation to get final results
!  ================================================

   if ( interpolate_aerosols ) then

!  Interpolation, wavelength regime

     if ( .not. do_wavnums ) then

       wastart = 1
       do w = 1, nwav
         wa = wastart ; trawl = .true.
         do while (trawl)
           if ( wav(w) .ge. aerwav(wa) .and. wav(w) .le.aerwav(wa+1) ) trawl = .false.
         enddo
         wa1 = wa ; wa2 = wa + 1
         fa1 = ( aerwav(wa2) - wav(w) ) / ( aerwav(wa2) -  aerwav(wa1) ) ; fa2 = 1.0d0 - fa1
         wastart = wa1
         if ( w.lt.nwav) then
           if(wav(w+1).ge.aerwav(wa+1))wastart = wa2
         endif
         aod_scaling(w) = fa1 * local_aodscaling(wa1) + fa2 * local_aodscaling(wa2)
         do n = 1, nlayers
           aerosol_deltau(n,w) = aertau_unscaled(n) * aod_scaling(w) 
         enddo
         aerosol_ssalbs(w) = fa1 * local_aerssalbs(wa1) + fa2 * local_aerssalbs(wa2)
         n_scatmoms(w) = min(local_nscatmoms(wa1),local_nscatmoms(wa2))
         do l = 0, n_scatmoms(w)
!mick mod 1/24/2017 - turn off do loop since f90 vector ops already in use here
           !do k = 1, nmuller
              aerosol_scatmoms(1:nmuller,l,w) = fa1 * local_scatmoms(1:nmuller,l,wa1) &
                                              + fa2 * local_scatmoms(1:nmuller,l,wa2)
           !enddo
         enddo

!         aerosol_asymms(w) = aerosol_scatmoms(1,1,w) / 3.0d0
         aerosol_scatmoms(1,0,w) = 1.0d0

!  Debug aerosol optical properties. VERSION TWO only
         if ( do_iopchecker ) then
            do n = 1, nlayers
              if (aerlayerflags(N).and.n.eq.nlayers ) then
                write(999,'(2i5,1p2e20.10,i5)')w,n,aerosol_deltau(n,w),aerosol_ssalbs(w), n_scatmoms(w)
                 do l = 0, n_scatmoms(w)
                    write(999,'(2i5,1p6e20.10)')n,l,(aerosol_scatmoms(k,l,w),k=1,1)
                 enddo
              endif
            enddo
         endif
!  End debug aerosol optical properties

       enddo
     endif

!  Interpolation, wavenumber regime
!   12 August 2013 --> KLUTZY CODE here,,,,,,,,,,Improve it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if ( do_wavnums ) then
       wastart = 1

!  Start wavenumber loop

       do w = 1, nwav

!  Define weights to interpolate from aerosol grid to RT grid

         wa = wastart ; trawl = .true. ; lam =  1.0d+07/wav(w)
54       continue
         do while (trawl)
           if ( lam .lt. aerwav(wa) ) then
             if ( lam .ge. aerwav(wa+1) ) then
               trawl = .false.
             else if (lam .lt. aerwav(wa+1) ) then
               wa = wa + 1 ; go to 54
             endif
           endif
         enddo
         wa1 = wa ; wa2 = wa + 1
         fa1 = ( aerwav(wa2) - lam ) / ( aerwav(wa2) -  aerwav(wa1) ) ; fa2 = 1.0d0 - fa1
         wastart = wa1

!  Interpolate regular quantities

         aod_scaling(w) = fa1 * local_aodscaling(wa1) + fa2 * local_aodscaling(wa2)
         do n = 1, nlayers
           aerosol_deltau(n,w) = aertau_unscaled(n) * aod_scaling(w)
         enddo
         aerosol_ssalbs(w) = fa1 * local_aerssalbs(wa1) + fa2 * local_aerssalbs(wa2)
         n_scatmoms(w) = min(local_nscatmoms(wa1),local_nscatmoms(wa2))
         do l = 0, n_scatmoms(w)
           do k = 1, nmuller
             aerosol_scatmoms(1:nmuller,l,w) = fa1 * local_scatmoms(1:nmuller,l,wa1) + fa2 * local_scatmoms(1:nmuller,l,wa2)
           enddo
         enddo
!         aerosol_asymms(w) = aerosol_scatmoms(1,1,w) / 3.0d0
         aerosol_scatmoms(1,0,w) = 1.0d0

!  End wavenumber loop

       enddo

!  End wavenumber regime

     endif

!  Finish aerosol interpolation

   endif

!  Safety - File Closures, before leaving routine
!  (note: these files opened @ ~line 700)

   if ( GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then
      do nm = 1, nmodes
         close(InFileUnit(nm))
      enddo
   else
      do nm = 1, nmodes
         if ( MakeAeroFile(nm) ) close(OutFileUnit(nm))
      enddo
   endif

!  End subroutine

   return

end subroutine GEMSTOOL_AER_PROPERTIES

!  End module

end Module GEMSTOOL_AerProperties_m
