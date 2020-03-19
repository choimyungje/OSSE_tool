Module aerosol_regular_creation_m

!  Notes: 18 February 2013
!  -----------------------

!     * profiles_lidar routine added
!        - Read LIDAR Extinction [km-1] and height [km] profiles from FILE
!        - Parcel the entire LIDAR profile into the output array
!        - Ignores z_upperlimit, z_lowerlimit
!        - Only one offline test so far (19 february, 2013)

!     * Add AEROSOL_DISTCHARS to output

!  Loading routines

   use aerosol_loading_routines_m, only: profiles_uniform, &
                                         profiles_expone,  &
                                         profiles_gdfone,  &
                                         profiles_lidar

!  Use Mie code
!  ------------

!  parameters

   use RTSMie_parameters_m

!  Single call

   use RTSMie_sourcecode_m

!  Bimodal call

   use RTSMie_master_bimodal_m

!  Use Tmatrix code
!  ----------------

!  parameters

   use tmat_parameters

!  Single call

   use tmat_master_m

!  Bimodal call

   use tmat_master_bimodal_m

!  All routines are public

private
public :: aerosol_regular_creation, get_variable_data, get_variable_props

contains
!        PSD_Index, PSD_pars, n_real, n_imag, FixR1R2, R1, R2, epsnum, epsfac,   & ! INPUT, Mie/Tmat Parameters

subroutine aerosol_regular_creation &
      ( maxlayers, maxwav, maxaermoms, do_bimodal, Loading_case, TmatMie_case,            & ! INPUT, dimensions, flags
        loading_upperboundary, loading_lowerboundary, aertau_input_w0, reference_w0,      & ! INPUT, Loading control
        exploading_relaxation, gdfloading_peakheight, gdfloading_halfwidth, LIDAR_file,   & ! INPUT, Loading control
        nlayers, nmuller, nwav, wav, height_grid, bimodal_fraction, momsize_cutoff,       & ! INPUT, Numbers, misc.
        do_Varprops, Max_Varpts, N_Varpts, Aerwavs_Var, FixR1R2, R1, R2,  & ! INPUT, Mie/Tmat Parameters            @@@ RTS 09/11/12
        PSDindex_fixed, PSDpars_fixed, nreal_fixed, nimag_fixed,          & ! INPUT, Mie/Tmat Parameters (fixed)    @@@ RTS 09/11/12
        PSDindex_var,   PSDpars_var,   nreal_var,   nimag_var,            & ! INPUT, Mie/Tmat Parameters (Variable) @@@ RTS 09/11/12
        R1R2_cutoff, nblocks, nweights, xparticle_limit,                                  & ! INPUT, Mie     PSD control
        do_EqSaSphere, spheroid_type, nkmax, ndgs, shape_factor, Tmat_accuracy,           & ! INPUT, Tmatrix PSD control 
        aerlayerflags, Loading, n_scatmoms, aertau_unscaled, aod_scaling,                 & ! OUTPUT, aerosol optical properties
        aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms, aerosol_distchars,              & ! OUTPUT, aerosol optical properties
        fail1, fail2, Message_Loading, Messages_Optical )                                   ! Exception handling

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
!       Case 4    : dloading_Dpars     = 0 (No Derivatives allowed, LIDAR profile). Added 18 February 2013

!  Generation of Aerosol optical properties
!     1. Call to the Mie/Tmatrix program
!     2. Convert Mie/Tmatrix output (Microsopic) to IOP output (macroscopic) 

   implicit none

!  Input variables
!  ---------------

!  External Dimensioning

   integer, INTENT (IN) :: maxlayers, maxwav
   integer, INTENT (IN) :: maxaermoms

!  top level flags

   logical, INTENT (IN) :: do_bimodal

!  top level index (1 = Tmat, 2 = Mie)

   integer, INTENT(IN) :: TmatMie_case

!  Aerosol loading control
!  TOTAL AOT = aerosol optical thickness (value at reference wavelength )
!    Case 1 = uniform aerosol 
!    Case 2 = Exponential Loading
!    Case 3 = GDF loading

   integer, INTENT(IN)       :: loading_case
   real(fpk),    INTENT (IN) :: aertau_input_w0
   real(fpk),    INTENT (IN) :: reference_w0

   real(fpk),    INTENT (IN) :: exploading_relaxation
   real(fpk),    INTENT (IN) :: gdfloading_peakheight
   real(fpk),    INTENT (IN) :: gdfloading_halfwidth
   real(fpk),    INTENT (IN) :: loading_upperboundary
   real(fpk),    INTENT (IN) :: loading_lowerboundary

!  @@@ New inputs (18 February 2013)
!      Character String for file of LIDAR height/extinction profile data (Loading case 4)

   character*(*), INTENT(IN) :: LIDAR_file

!  Numbers

   integer, INTENT (IN) ::  nlayers
   integer, INTENT (IN) ::  nwav

!  Heights and wavelengths

   REAL    (fpk),    INTENT (IN)   :: wav(maxwav)
   REAL    (fpk),    INTENT (IN)   :: height_grid(0:maxlayers)

!  nmuller = 1 (scalar code), = 6 (vector code)

   integer, INTENT (IN) ::  nmuller

!  Bimodal fraction F (that is, total = F.Mie_1 + (1-F).Mie_2 )

   real(fpk),    INTENT (IN) :: bimodal_fraction

!  Cutoff size for smallest phase function moment

   real(fpk),    INTENT (IN) :: momsize_cutoff

!  Fd testing

!   integer     , intent(in)  :: epsnum
!   real(kind=8), INTENT (IN) :: epsfac

!  Mie/Tmatrix PSD inputs
!  ======================

!  Number of points for interpolation of Variable aerosol properties
!    ( @@@ RTS 09/11/12, Use of variable aerosol input )

   logical, INTENT (IN)  :: do_Varprops             ! New @@@ RTS 09/11/12
   integer  , intent(in) :: Max_Varpts, N_Varpts    ! New @@@ RTS 09/11/12
   real(fpk), intent(in) :: Aerwavs_Var(Max_Varpts) ! New @@@ RTS 09/11/12

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

   integer         , intent(in)  :: PSDIndex_fixed(2)               ! Renamed, @@@ RTS 09/11/12
   real    (fpk),    intent(in)  :: PSDpars_fixed (3,2)             ! Renamed, @@@ RTS 09/11/12
   integer         , intent(in)  :: PSDIndex_Var(Max_Varpts,2)      ! New, @@@ RTS 09/11/12
   real    (fpk),    intent(in)  :: PSDpars_Var (Max_Varpts,3,2)    ! New, @@@ RTS 09/11/12

!  FixR1R2 : If  set, Use Internal routine to calculate R1 and R2 (outputs)
!            If Not set, Use Input R1 and R2 for PSD limits.

   logical         , intent(inout)  :: FixR1R2(2)

!  R1, R2         - Minimum and Maximum radii (Microns)

   real    (fpk),    intent(inout)  :: R1(2), R2(2)

!  N_REAL, N_IMAG - real and imaginary parts, refractive index (N-i.GE.0)

   real    (fpk),    intent(in)  :: nreal_fixed(2), nimag_fixed(2)                   ! Renamed, @@@ RTS 09/11/12
   real    (fpk),    intent(in)  :: nreal_Var(Max_Varpts,2), nimag_Var(Max_Varpts,2) ! New, @@@ RTS 09/11/12

!  Mie-specific inputs
!  ===================

!  Limiting particle size value. Set to 10000.0 default.
!   If you exceed this, program will tell you to increase dimensioning.

   REAL    (fpk),    intent(inout)   :: xparticle_limit

!  R1R2_cutoff particle size for setting R1 and R2 internally

   REAL    (fpk),    intent(in)  :: R1R2_cutoff(2)

!  PSD quadrature control
!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.

   INTEGER         , intent(in)  :: nblocks(2)
   INTEGER         , intent(in)  :: nweights(2)

!  Tmatrix-specific inputs
!  =======================

!  Logical flag for using equal-surface-area sepcification

   logical  , intent(in)  :: Do_EqSaSphere

!      NKMAX.LE.988 is such that NKMAX+2 is the                        
!           number of Gaussian quadrature points used in               
!           integrating over the size distribution for particles
!           MKMAX should be set to -1 for Monodisperse

   integer  , intent(inout)  :: nkmax(2)

!      NDGS - parameter controlling the number of division points      
!             in computing integrals over the particle surface.        
!             For compact particles, the recommended value is 2.       
!             For highly aspherical particles larger values (3, 4,...) 
!             may be necessary to obtain convergence.                  
!             The code does not check convergence over this parameter. 
!             Therefore, control comparisons of results obtained with  
!             different NDGS-values are recommended.

   integer  , intent(in)     :: ndgs(2)

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

   integer  , intent(in)  :: spheroid_type
   real(fpk), intent(in)  :: shape_factor(2)

!      Accuracy       - accuracy of the computations

   real(fpk), intent(in)  :: Tmat_accuracy


!  output
!  ======

!  Loading output

    real(fpk),    dimension ( maxlayers ),  INTENT (OUT)      :: Loading

!  aerosol layer flags

    LOGICAL,      DIMENSION ( maxlayers ), intent(out)  :: AERLAYERFLAGS 

!  AOP output: Number of exapnsion coefficients to be used

   integer, INTENT (OUT)   ::  n_scatmoms

!  Unscaled profiles of optical depth

    real(fpk),    INTENT (OUT)  :: aertau_unscaled(maxlayers)

!  AOD scaling

    real(fpk),    INTENT (OUT)  :: aod_scaling(maxwav)

!  Aod output, final

    REAL(fpk),    DIMENSION( maxlayers, maxwav )      , INTENT (OUT)  :: AEROSOL_DELTAU

!  AOPS

    REAL(fpk),    DIMENSION( maxwav )                 , INTENT (OUT)  :: AEROSOL_SSALBS 
    REAL(fpk),    DIMENSION( 6, 0:maxaermoms, maxwav ), INTENT (OUT)  :: AEROSOL_SCATMOMS

!  >>>>>>>>>>>> Added, 18 February 2013
!  Aerosol distribution characterisstics.
!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF
   real(fpk),     DIMENSION( 5, 2, maxwav)            , INTENT (OUT)  :: AEROSOL_DISTCHARS
!  >>>>>>>>>>>> Added, 18 February 2013

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

!  NPNA - number of equidistant scattering angles (from 0
!             to 180 deg) for which the scattering matrix is
!             calculated.

!    NOT REQUIRED HERE

   INTEGER           :: n_Fmatrix_angles
   REAL    (fpk)     :: Fmatrix_angles(max_Mie_angles)

!  Monoradius     - Monodisperse radius size (Microns)

   real    (fpk)     :: Monoradius

!  (Tmatrix only). Style flag.
!    * This is checked and re-set (if required) for Monodisperse case

   logical           :: Do_psd_OldStyle

!  Mie code OUTPUT variables
!  =========================

!  Bulk distribution parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk)    :: BMie_bulk (3)

!  Expansion coefficients and Asymmetry parameter

   integer      :: BMie_ncoeffs
   real(fpk)    :: BMie_expcoeffs (6,0:max_Mie_angles)
   real(fpk)    :: BMie_asymm

!  F-matrix,  optional output

   real(fpk)    :: BMie_Fmatrix(4,max_Mie_angles)

!  Distribution parameters
!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(fpk)     :: BMie_dist (5,2)

!  Exception handling

   character*90  :: Mie_Bmessages(3)
   character*90  :: Mie_trace_3

!  Tmatrix code OUTPUT variables (Notation similar)
!  =============================

!  Bulk distribution parameters

   real(fpk)    :: BTmat_bulk (3)

!  Expansion coefficients and Asymmetry parameter

   integer      :: BTmat_ncoeffs
   real(fpk)    :: BTmat_expcoeffs (NPL1,6)
   real(fpk)    :: BTmat_asymm

!  F-matrix,  optional output

   real(fpk)    :: BTmat_Fmatrix(MAXNPA,6)

!  Distribution parameters

   real(fpk)     :: BTmat_dist (5,2)

!  Exception handling

   character*90  :: Tmat_message
   character*90  :: Tmat_trace
   character*90  :: Tmat_trace_2
   character*90  :: Tmat_trace_3

!  Local variables
!  ===============

!  Local aerosol properties ! New, @@@ RTS 09/11/12

   integer   :: PSD_index(2)              ! New, @@@ RTS 09/11/12
   real(fpk) :: PSD_pars(3,2)             ! New, @@@ RTS 09/11/12
   real(fpk) :: n_real(2), n_imag(2)      ! New, @@@ RTS 09/11/12

!  AOP output: Reference values
!    * These are the values at reference wavelength w0
!    * They are set for wavelength w0, then used again

    real(fpk)     :: extinction_ref

!  Help Variables

   character*3  :: cwav
   integer      :: n, k, l, m, istatus, lambda_index, w, wc, wavmask(maxwav), n_scatmoms_w
   real(fpk)    :: wavelength, Mie_wavelength, Tmat_wavelength, extinction

!  Local control

   logical, parameter :: do_iopchecker    = .false.
!   logical, parameter :: do_iopchecker    = .true.
   logical, parameter :: do_aer_Jacobians = .false.

!  Initialize output
!  =================

!  Initialize exception handling

   fail1 = .false.
   fail2 = .false.
   Message_Loading   = ' '
   Messages_Optical  = ' '

!  initialize Aerosol Loading

   Loading        = 0.0d0
   Dloading_Dtau  = 0.0d0     ! Not required
   Dloading_Dpars = 0.0d0
   aerlayerflags = .false.

!  Set local Mie/Tmatrix variables

   do_monodisperse = .false.             ! Always false
   Do_Fmatrix      = .false.             ! Always false
   DO_Expcoeffs    = .false.             ! this will be set later on
   monoradius      = 0.0d0
   do_psd_oldstyle = .false.

!  Initialize optical properties

   extinction_ref    = 0.0d0
   aod_scaling       = 0.0d0
   aerosol_deltau    = 0.0d0
   aerosol_ssalbs    = 0.0d0
   aerosol_scatmoms  = 0.0d0
   aerosol_distchars = 0.0d0

!  Initialize local aerosol properties (Just a precaution)! New, @@@ RTS 09/11/12

   n_real = 0.0_fpk ; n_imag   = 0.0_fpk
   PSD_Index = 0    ; PSD_pars = 0.0_fpk

!  Now form the aerosol Loading
!  ----------------------------

!  @@@ Notes: 18 February 2013
!        profiles_lidar routine added (loading case 4)
!         - Read LIDAR Extinction [km-1] and height [km] profiles from FILE
!         - Parcel the entire LIDAR profile into the output array
!         - Ignores z_upperlimit, z_lowerlimit
!         - Only one offline test so far (19 february, 2013)

!    Loading = optical depth profile

!    Derivatives : 
!       Cases 1-3 : dloading_Dtau      = derivative of profile w.r.t the total aerosol optical depth at wavelength w0
!       Case 2    : dloading_Dpars(1)  = derivative of profile w.r.t the relaxation parameter (Exponential)
!       Case 3    : dloading_Dpars(1)  = derivative of profile w.r.t the GDF Peak height
!       Case 3    : dloading_Dpars(2)  = derivative of profile w.r.t the GDF Half Width 
!       Case 4    : dloading_Dpars     = 0 (No Derivatives allowed, LIDAR profile)

!  Case 1: Uniform layer of aerosol
!  ********************************

   if ( loading_case .eq. 1 ) then

      CALL profiles_uniform                                                 &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,              & ! Inputs
            loading_upperboundary, loading_lowerboundary, aertau_input_w0,  & ! Inputs
            Loading, Dloading_Dtau,                                         & ! output
            fail1, message_Loading(1), message_Loading(2) )                   ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'Uniform aerosol Loading failed'

!  Case 2: Exponential decay profile
!  *********************************

   else if ( loading_case .eq. 2 ) then

      CALL profiles_expone                                                  &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,              & ! Inputs
            loading_upperboundary, loading_lowerboundary,                   & ! Inputs
            exploading_relaxation, aertau_input_w0,                         & ! Inputs
            Loading, Dloading_Dtau, Dloading_Dpars(:,1),                    & ! output
            fail1, message_Loading(1), message_Loading(2) )                   ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'Exponential aerosol Loading failed'

!  Case 3: GDF (quasi-Gaussian) profile
!  ************************************

   else if ( loading_case .eq. 3 ) then

      CALL profiles_gdfone                                                        &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,                    & ! Inputs
            loading_upperboundary, gdfloading_peakheight, loading_lowerboundary,  & ! Inputs
            gdfloading_halfwidth, aertau_input_w0,                                & ! Inputs
            Loading, Dloading_Dtau, Dloading_Dpars(:,1), Dloading_Dpars(:,2),     & ! output
            fail1, message_Loading(1), message_Loading(2) )                         ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'GDF aerosol Loading failed'

!  Case 4: LIDAR profile from file
!  *******************************

   else if ( loading_case .eq. 4 ) then

      CALL profiles_lidar &
      ( maxlayers, nlayers, LIDAR_file, height_grid,    & ! Inputs
        Loading, Dloading_Dtau, Dloading_Dpars,         & ! Outputs
        fail1, message_Loading(1), message_Loading(2) )   ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'LIDAR aerosol profile (Case 4) failed'

   endif

!  Return if failure at this stage

   if ( fail1 ) return

!  pause' aerosol loading'

!  Debug loading derivatives
!     do n = 91, 113,10
!     write(*,466)n,Loading(n)
!     enddo
!466 format(i4,1pe20.10)

!  Check to see if reference wavelength is one of the set. Initialize Mask.

   lambda_index = 0
   do w = 1, nwav
     wavmask(w) = w
     if ( wav(w) .eq. reference_w0 ) lambda_index = w
   enddo
!   write(*,*)lambda_index

!  Mask to use if reference wavelength is one of list of wavelengths

   if ( lambda_index .ne. 0 ) then
     wavmask(1) = lambda_index
     wc = 1
     do w = 1, lambda_index - 1
       wc = wc + 1
       wavmask(wc) = w
     enddo
     do w = lambda_index + 1, nwav
       wc = wc + 1
       wavmask(wc) = w
     enddo
   endif

!  Initialize

   n_scatmoms = 50000

!  Skip Mie section if Doing Tmatrix

   if ( TmatMie_case .eq. 1 ) go to 67

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      M I E   C a l c u l a t i o n 
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  FD perturbation for the Mie aerosol inputs.

!   if ( epsnum .ne. 0 ) then
!      if ( epsnum .eq. 1 ) then
!        n_real(1) = n_real(1) * epsfac
!      else if ( epsnum .eq. 2 ) then
!        n_imag(1) = n_imag(1) * epsfac
!      else if ( epsnum .eq. 3 ) then
!        psd_pars(1,1) = psd_pars(1,1) * epsfac
!      else if ( epsnum .eq. 4 ) then
!        psd_pars(2,1) = psd_pars(2,1) * epsfac
!      else if ( epsnum .eq. 5  .and. do_bimodal ) then
!        n_real(2) = n_real(2) * epsfac
!      else if ( epsnum .eq. 6  .and. do_bimodal ) then
!        n_imag(2) = n_imag(2) * epsfac
!      else if ( epsnum .eq. 7  .and. do_bimodal ) then
!        psd_pars(1,2) = psd_pars(1,2) * epsfac
!      else if ( epsnum .eq. 8  .and. do_bimodal ) then
!        psd_pars(2,2) = psd_pars(2,2) * epsfac
!      else if ( epsnum .eq. 9  .and. do_bimodal ) then
!        bimodal_fraction = bimodal_fraction * epsfac
!      endif
!   endif

!  Prepare the reference-wavelength Mie Inputs
!  ===========================================

   if ( TmatMie_case .eq. 2 .and. lambda_index .eq. 0 ) then

!  Only require extinction coefficient if flagged
!  Set the local Mie program inputs (bulk properties only)

      Do_Expcoeffs     = .FALSE.

!  reference wavelength

      wavelength     = reference_w0
      Mie_wavelength = wavelength/1000.0d0

!  Call to assign local aerosol properties ! New, @@@ RTS 09/11/12

      call get_variable_props &
        ( do_Varprops, do_bimodal, Max_Varpts, N_Varpts, wavelength,    &
          nreal_fixed, nimag_fixed, psdindex_fixed, psdpars_fixed,      &
          Aerwavs_Var, nreal_var, nimag_var, psdindex_var, psdpars_var, &
          PSD_Index, PSD_pars, n_real, n_imag )

!  progress

!      write(*,*)'Regular Mie Aerosol calculation, reference wavelength # '

!  BIMODAL vs. SINGLE CALL

      if ( do_bimodal) then
         call RTSMie_master_bimodal                               &
            ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
              PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,   & ! I
              nblocks, nweights, xparticle_limit, R1R2_cutoff,    & ! I
              n_Fmatrix_angles, Fmatrix_angles,                   & ! I
              Mie_wavelength, n_real, n_imag, bimodal_fraction,   & ! I
              BMie_bulk, BMie_asymm, BMie_ncoeffs,                & ! O
              BMie_expcoeffs, BMie_Fmatrix, BMie_dist,            & ! O
              fail2, istatus, Mie_Bmessages, Mie_trace_3 )          ! O
      else
         BMie_dist(:,2) = 0.0d0
         call RTSMie_main                                                                          & !---MIE CALL
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,                                         & ! I
                PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),                 & ! I
                nblocks(1), nweights(1), xparticle_limit, R1R2_cutoff(1),                          & ! I
                n_Fmatrix_angles, Fmatrix_angles, Mie_wavelength, n_real(1), n_imag(1),            & ! I
                BMie_bulk, BMie_asymm, BMie_ncoeffs, BMie_expcoeffs, BMie_Fmatrix, BMie_dist(:,1), & ! O
                fail2, istatus, Mie_Bmessages(1), Mie_Bmessages(2), Mie_Bmessages(3) )               ! O
      endif

!  Exception handling on everything

      if ( Fail2 ) then  
         do m = 1, 3   
            Messages_Optical(m) = adjustl(trim(Mie_Bmessages(m)))
         enddo
         Messages_Optical(4) = 'Single Mie Call'; if ( do_bimodal ) Messages_Optical(4) = adjustl(trim(Mie_trace_3))
         Messages_Optical(5) = 'First call to the Regular Mie program in AEROSOL CREATION, reference wavelength'
         return
      endif

!  Set the reference quantity

      extinction_ref = BMie_bulk(1)

!      write(*,*)do_Varprops, do_bimodal, Max_Varpts, N_Varpts
!      write(*,*)n_real(1),n_real(2),n_imag(1),n_imag(2)
!      write(*,*)wavelength,PSD_Index, PSD_pars(1,1),PSD_pars(2,1)

!     write(*,*)'Ext Ref 0',extinction_ref; pause

!  End Mie reference-wavelength calculation

   endif

!  Prepare General (all-wavelength) Mie inputs
!  ===========================================

   if ( TmatMie_case .eq. 2 ) then

!  wavelength loop. First wavelength will be the reference, if in  list.

      do wc = 1, nwav

!  Wavelengths [nm]

         w = wavmask(wc)
         wavelength = wav(w)

!  Call to assign local aerosol properties ! New, @@@ RTS 09/11/12

         call get_variable_props &
        ( do_Varprops, do_bimodal, Max_Varpts, N_Varpts, wavelength,    &
          nreal_fixed, nimag_fixed, psdindex_fixed, psdpars_fixed,      &
          Aerwavs_Var, nreal_var, nimag_var, psdindex_var, psdpars_var, &
          PSD_Index, PSD_pars, n_real, n_imag )

!  progress

         write(*,*)'Regular Mie Aerosol calculation, doing wavelength # ',wc,wav(w)

!  wavelength for the Mie code (micron unit)  

         Mie_wavelength = wavelength/1000.0d0

!  Set the local Mie program inputs (general)

         Do_Expcoeffs     = .TRUE.

!  BIMODAL vs. SINGLE CALL

         if ( do_bimodal) then
            call RTSMie_master_bimodal                            &
            ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
              PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,   & ! I
              nblocks, nweights, xparticle_limit, R1R2_cutoff,    & ! I
              n_Fmatrix_angles, Fmatrix_angles,                   & ! I
              Mie_wavelength, n_real, n_imag, bimodal_fraction,   & ! I
              BMie_bulk, BMie_asymm, BMie_ncoeffs,                & ! O
              BMie_expcoeffs, BMie_Fmatrix, BMie_dist,            & ! O
              fail2, istatus, Mie_Bmessages, Mie_trace_3 )          ! O
         else
           BMie_dist(:,2) = 0.0d0
           call RTSMie_main                                                                       & !---MIE CALL
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,                                         & ! I
                PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),                 & ! I
                nblocks(1), nweights(1), xparticle_limit, R1R2_cutoff(1),                          & ! I
                n_Fmatrix_angles, Fmatrix_angles, Mie_wavelength, n_real(1), n_imag(1),            & ! I
                BMie_bulk, BMie_asymm, BMie_ncoeffs, BMie_expcoeffs, BMie_Fmatrix, BMie_dist(:,1), & ! O
                fail2, istatus, Mie_Bmessages(1), Mie_Bmessages(2), Mie_Bmessages(3) )               ! O
        endif

!  Exception handling on everything

        if ( Fail2 ) then
          write(cwav,'(I3)')wc
          do m = 1, 3
            Messages_Optical(m) = adjustl(trim(Mie_Bmessages(m)))
          enddo
          Messages_Optical(4) = 'Single Mie Call'; if ( do_bimodal ) Messages_Optical(4) = adjustl(trim(Mie_trace_3))
          Messages_Optical(5) = 'First call to the Regular Mie program in AEROSOL CREATION, wavelength # '//cwav
          return
        endif

!  NOW set the optical property output
!  ===================================

!  >>>>>>>>>>>>>> New, 18 February 2013
!  set the distribution characteristics

         aerosol_distchars(:,:,w) = BMie_dist

!  Set the reference quantities, if reference wavelength is in the list.
!   Values for the first (masked) wavelength

         if ( lambda_index .ne. 0 .and. wc.eq.1 ) then
            extinction_ref = BMie_bulk(1)
         endif

!  Extinction and its scaling factor. All wavelengths.

         extinction = BMie_bulk(1)
         aod_scaling(w) = extinction / extinction_ref

!  Assign SSAs and expansion coefficients, single/bimodal aerosol type

         aerosol_ssalbs(w) = BMie_bulk(3)
         l = 0
         aerosol_scatmoms(1,0,w) = 1.0d0
         do while (aerosol_scatmoms(1,l,w).gt.momsize_cutoff.and.l.lt.maxaermoms )
            l = l + 1
            aerosol_scatmoms(1,l,w) = Bmie_expcoeffs(1,l)
         enddo
         n_scatmoms_w = l
         do l = 0, n_scatmoms_w
            do k = 2, nmuller
               aerosol_scatmoms(k,l,w) = BMie_expcoeffs(k,l) 
            enddo
         enddo

!  Update n_scatmoms

         n_scatmoms = min(n_scatmoms, n_scatmoms_w)

!  Assign Unscaled loadings

         aertau_unscaled = 0.0d0
         do n = 1, nlayers
            aerlayerflags(n) =  ( Loading(n) .ne. 0.0d0  )
            aertau_unscaled(n) = loading(n)
         enddo

!  Apply scalings to loadings

         do n = 1, nlayers
            aerosol_deltau(n,w) = aertau_unscaled(n) * aod_scaling(w)
         enddo

!  debug aerosol optical properties. VERSION TWO only

         if ( do_iopchecker ) then
           do n = 1, nlayers
            if (aerlayerflags(N).and.n.eq.107 ) then
              write(999,'(i4,1p6e20.10)')n,aerosol_deltau(n,w),aerosol_ssalbs(w)
              do l = 0, n_scatmoms_w
                write(999,'(2i5,1p6e20.10)')n,l,(aerosol_scatmoms(k,l,w),k=1,1)
              enddo
!            else
!              write(999,'(i4,1p6e20.10)')n,aerosol_deltau(n,w)
            endif
           enddo
!           pause'Reg 999'
         endif

!  End wavelength loop

      enddo

!  End Mie aerosol clause

   endif

!  Mie return

   return

!  Continuation point for Tmatrix calculation

67 continue

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      T M A T R I X    C a l c u l a t i o n 
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  FD perturbation for the Tmat aerosol inputs.

!   if ( epsnum .ne. 0 ) then
!      if ( epsnum .eq. 1 ) then
!        n_real(1) = n_real(1) * epsfac
!      else if ( epsnum .eq. 2 ) then
!        n_imag(1) = n_imag(1) * epsfac
!      else if ( epsnum .eq. 3 ) then
!        Tmat_eps(1) = Tmat_eps(1) * epsfac
!      else if ( epsnum .eq. 4 ) then
!        psd_pars(1,1) = psd_pars(1,1) * epsfac
!      else if ( epsnum .eq. 5 ) then
!        psd_pars(2,1) = psd_pars(2,1) * epsfac
!      else if ( epsnum .eq. 6 .and. do_bimodal ) then
!        n_real(2) = n_real(2) * epsfac
!      else if ( epsnum .eq. 7 .and. do_bimodal ) then
!        n_imag(2) = n_imag(2) * epsfac
!      else if ( epsnum .eq. 8 .and. do_bimodal ) then
!        Tmat_eps(2) = Tmat_eps(2) * epsfac
!      else if ( epsnum .eq. 9 .and. do_bimodal ) then
!        psd_pars(1,2) = psd_pars(1,2) * epsfac
!      else if ( epsnum .eq. 10 .and. do_bimodal ) then
!        psd_pars(2,2) = psd_pars(2,2) * epsfac
!      else if ( epsnum .eq. 11 .and. do_bimodal ) then
!        bimodal_fraction = bimodal_fraction * epsfac
!      endif
!   endif

!  Prepare the reference-wavelength Tmatrix Inputs
!  ===============================================

   if ( TmatMie_case .eq. 1 .and. lambda_index .eq. 0 ) then

!  Only require extinction coefficient if flagged
!  Set the local Tmatrix program inputs (bulk properties only)

      Do_Expcoeffs     = .FALSE.

!  reference wavelength

      wavelength      = reference_w0
      Tmat_wavelength = wavelength/1000.0d0

!  Call to assign local aerosol properties ! New, @@@ RTS 09/11/12

      call get_variable_props &
        ( do_Varprops, do_bimodal, Max_Varpts, N_Varpts, wavelength,    &
          nreal_fixed, nimag_fixed, psdindex_fixed, psdpars_fixed,      &
          Aerwavs_Var, nreal_var, nimag_var, psdindex_var, psdpars_var, &
          PSD_Index, PSD_pars, n_real, n_imag )

!  progress

      write(*,*)'   Regular Tmatrix Aerosol calculation, reference wavelength # '

!  BIMODAL vs. SINGLE CALL

      if ( do_bimodal) then
         call tmat_master_bimodal                 &
               ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere,              & ! Inputs (Flags)
                 Do_psd_OldStyle, psd_Index, psd_pars, MonoRadius, R1, R2, FixR1R2,     & ! Inputs (PSD)
                 bimodal_fraction, Spheroid_type, nkmax, n_Fmatrix_angles, ndgs,        & ! Inputs (General)
                 shape_factor, Tmat_accuracy, Tmat_wavelength, n_real, n_imag,          & ! Inputs (Optical)
                 Btmat_bulk, Btmat_asymm, Btmat_ncoeffs,                                & ! Outputs (Tmat)
                 Btmat_expcoeffs, Btmat_Fmatrix, Btmat_dist,                            & ! Outputs (PSD)
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2, Tmat_trace_3 )   ! Outputs (status)
      else
         Btmat_dist(:,2) = 0.0d0
         call tmat_master                         &
               ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere, Do_psd_OldStyle, & ! Inputs (Flags)
                 PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),         & ! Inputs (PSD)
                 Spheroid_type, nkmax(1), n_Fmatrix_angles, ndgs(1),                        & ! Inputs (General)
                 shape_factor(1), Tmat_accuracy, Tmat_wavelength, n_real(1), n_imag(1),     & ! Inputs (Optical)
                 Btmat_bulk, Btmat_asymm, Btmat_ncoeffs,                                    & ! Outputs (Tmat)
                 Btmat_expcoeffs, Btmat_Fmatrix, Btmat_dist(:,1),                           & ! Outputs (PSD)
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2 )                     ! Outputs 
      endif

!  Exception handling on everything

      if ( Fail2 ) then
         Messages_Optical(1) = adjustl(trim(Tmat_message))
         Messages_Optical(2) = adjustl(trim(Tmat_trace))
         Messages_Optical(3) = adjustl(trim(Tmat_trace_2))
         Messages_Optical(4) = 'Single Tmatrix Call'; if ( do_bimodal ) Messages_Optical(4) = adjustl(trim(Tmat_trace_3))
         Messages_Optical(5) = 'First call to the Regular Tmatrix program in AEROSOL CREATION, reference wavelength'
         return
      endif

!  Set the reference quantity

      extinction_ref = BTmat_bulk(1)

!  End Tmat reference-wavelength calculation

   endif

!  Prepare General (all-wavelength) Tmat inputs
!  ===========================================

   if ( TmatMie_case .eq. 2 ) then

!  wavelength loop. First wavelength will be the reference, if in  list.

      do wc = 1, nwav

!  Wavelengths [nm]

         w = wavmask(wc)
         wavelength = wav(w)

!  progress

         write(*,*)'   Regular Tmat Aerosol calculation, doing wavelength # ',wc

!  wavelength for the Tmat code (micron unit)  

         Tmat_wavelength = wavelength/1000.0d0

!  Call to assign local aerosol properties ! New, @@@ RTS 09/11/12

         call get_variable_props &
        ( do_Varprops, do_bimodal, Max_Varpts, N_Varpts, wavelength,    &
          nreal_fixed, nimag_fixed, psdindex_fixed, psdpars_fixed,      &
          Aerwavs_Var, nreal_var, nimag_var, psdindex_var, psdpars_var, &
          PSD_Index, PSD_pars, n_real, n_imag )

!  Set the local Tmatrix program inputs (general)

         Do_Expcoeffs     = .TRUE.

!  BIMODAL vs. SINGLE CALL

         if ( do_bimodal) then
            call tmat_master_bimodal                 &
               ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere,              & ! Inputs (Flags)
                 Do_psd_OldStyle, psd_Index, psd_pars, MonoRadius, R1, R2, FixR1R2,     & ! Inputs (PSD)
                 bimodal_fraction, Spheroid_type, nkmax, n_Fmatrix_angles, ndgs,        & ! Inputs (General)
                 shape_factor, Tmat_accuracy, Tmat_wavelength, n_real, n_imag,          & ! Inputs (Optical)
                 Btmat_bulk, Btmat_asymm, Btmat_ncoeffs,                                & ! Outputs (Tmat)
                 Btmat_expcoeffs, Btmat_Fmatrix, Btmat_dist,                            & ! Outputs (PSD)
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2, Tmat_trace_3 )   ! Outputs (status)
         else
            Btmat_dist(:,2) = 0.0d0
            call tmat_master                         &
               ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere, Do_psd_OldStyle, & ! Inputs (Flags)
                 PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),         & ! Inputs (PSD)
                 Spheroid_type, nkmax(1), n_Fmatrix_angles, ndgs(1),                        & ! Inputs (General)
                 shape_factor(1), Tmat_accuracy, Tmat_wavelength, n_real(1), n_imag(1),     & ! Inputs (Optical)
                 Btmat_bulk, Btmat_asymm, Btmat_ncoeffs,                                    & ! Outputs (Tmat)
                 Btmat_expcoeffs, Btmat_Fmatrix, Btmat_dist(:,1),                           & ! Outputs (PSD)
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2 )                     ! Outputs 
         endif

!  Exception handling on everything

         if ( Fail2 ) then     
            write(cwav,'(I3)')wc
            Messages_Optical(1) = adjustl(trim(Tmat_message))
            Messages_Optical(2) = adjustl(trim(Tmat_trace))
            Messages_Optical(3) = adjustl(trim(Tmat_trace_2))
            Messages_Optical(4) = 'Single Tmatrix Call'; if ( do_bimodal ) Messages_Optical(4) = adjustl(trim(Tmat_trace_3))
            Messages_Optical(5) = 'First call to the Regular Tmatrix program in AEROSOL CREATION, reference wavelength'
            return
         endif

!  NOW set the optical property output
!  ===================================

!  >>>>>>>>>>>>>> New, 18 February 2013
!  set the distribution characteristics

         aerosol_distchars(:,:,w) = Btmat_dist

!  Set the reference quantities, if reference wavelength is in the list.
!   Values for the first (masked) wavelength

         if ( lambda_index .ne. 0 .and. wc.eq.1 ) then
            extinction_ref = BTmat_bulk(1)
         endif

!  Extinction and its scaling factor. All wavelengths.

         extinction = BTmat_bulk(1)
         aod_scaling(w) = extinction / extinction_ref

!  Assign SSAs and expansion coefficients, single/bimodal aerosol type

         aerosol_ssalbs(w) = BTmat_bulk(3)
         l = 0
         aerosol_scatmoms(1,0,w) = 1.0d0
         do while (aerosol_scatmoms(1,l,w).gt.momsize_cutoff.and.l.lt.maxaermoms )
            l = l + 1
            aerosol_scatmoms(1,l,w) = BTmat_expcoeffs(1,l)
         enddo
         n_scatmoms_w = l
         do l = 0, n_scatmoms_w
            do k = 2, nmuller
               aerosol_scatmoms(k,l,w) = BTmat_expcoeffs(k,l) 
            enddo
         enddo

!  Update n_scatmoms

         n_scatmoms = min(n_scatmoms, n_scatmoms_w)

!  Assign Unscaled loadings

         aertau_unscaled = 0.0d0
         do n = 1, nlayers
            aerlayerflags(n) =  ( Loading(n) .ne. 0.0d0  )
            aertau_unscaled(n) = loading(n)
         enddo

!  Apply scalings to loadings

         do n = 1, nlayers
            aerosol_deltau(n,w) = aertau_unscaled(n) * aod_scaling(w)
         enddo

!  debug aerosol optical properties. VERSION TWO only

         if ( do_iopchecker ) then
           do n = 1, nlayers
            if (aerlayerflags(N) ) then
              write(888,'(i4,1p6e20.10)')n,aerosol_deltau(n,w),aerosol_ssalbs(w)
              do l = 0, n_scatmoms_w
                write(888,'(2i5,1p6e20.10)')n,l,(aerosol_scatmoms(k,l,w),k=1,1)
              enddo
            else
              write(888,'(i4,1p6e20.10)')n,aerosol_deltau(n,w)
            endif
           enddo
           pause'Reg 888'
         endif

!  End wavelength loop

      enddo

!  End Tmat aerosol clause

   endif

!  Finish

   return
end subroutine aerosol_regular_creation

subroutine get_variable_data &
   ( do_Varprops, do_bimodal, Max_vpts, Files,                           &
     N_vpts, wavs_var, nreal_var, nimag_var, psdindex_var, psdpars_var,  &
     fail, message )

!  Subroutine for extracting variable-wavelength aerosol properties from data

   implicit none

!  inputs
!  ------

!  Control

   character*(*), intent(in) :: files(2)
   logical  , intent(in) :: do_Varprops, do_bimodal
   integer  , intent(in) :: Max_vpts

!  Outputs
!  -------

!  Variable-property

   integer  , intent(out) :: N_vpts
   integer  , intent(out) :: psdindex_var(Max_vpts,2)
   real(fpk), intent(out) :: psdpars_var(Max_vpts,3,2)
   real(fpk), intent(out) :: nreal_var(Max_vpts,2), nimag_var(Max_vpts,2)
   real(fpk), intent(out) :: wavs_var(Max_vpts)

!  Exception handling

   logical, intent(inout)     :: fail
   character*(*), intent(out) :: message 
 
!  Local
!  -----

   character*1 :: b
   integer     :: q, w, npars, N_vpts_2
   real*8      :: dum

!  Zero output

   nreal_var= 0.0_fpk ; nimag_var  = 0.0_fpk
   psdindex_var = 0   ; psdpars_var= 0.0_fpk
   fail = .false. ; message = ' '
   if ( .not. do_Varprops ) return

!  First file

   b = '1' ; npars = 0
   open(1,file='physics_data/'//adjustl(trim(files(1))),err=98,status='old')
   read(1,*,err=97)N_vpts
   if ( N_vpts .lt. 1 .or. N_vpts .gt. Max_vpts ) then
      message = 'Number of variable points is wrong (too large or zero), first file'
      fail = .true. ; close(1) ; return
   endif
   do w = 1, N_vpts
      read(1,*,err=97)wavs_var(w),nreal_var(w,1),nimag_var(w,1),psdindex_var(w,1),npars,(psdpars_var(w,q,1),q=1,npars)
   enddo
   close(1)

!  second file - Should have the same wavelength points as the first.

   if ( do_bimodal ) then
      b = '2'
      open(1,file='physics_data/'//adjustl(trim(files(2))),err=98,status='old')
      read(1,*,err=97)N_vpts_2
      if ( N_vpts_2 .lt. 1 .or. N_vpts_2 .gt. Max_vpts ) then
         message = 'Number of variable points is wrong (too large or zero), second file'
         fail = .true. ; close(1) ; return
      endif
      if ( N_vpts_2 .ne. N_vpts ) then
         message = 'Number of variable points is not equal to first value, second file'
         fail = .true. ; close(1) ; return
      endif
      do w = 1, N_vpts
         read(1,*,err=97)dum,nreal_var(w,2),nimag_var(w,2),psdindex_var(w,2),npars,(psdpars_var(w,q,2),q=1,npars)
      enddo
      close(1)
   endif

!  Finish

   return

!  Error finishes

98 continue
   message = 'Aerosol NRNI/PSD properties: Open file error for source number '//b
   fail = .true. ; return
97 continue
   message = 'Aerosol NRNI/PSD properties: Read file error for source number '//b
   fail = .true. ; return
 
end subroutine get_variable_data

subroutine get_variable_props &
   ( do_Varprops, do_bimodal, Max_vpts, N_vpts, lambda,        &
     nreal_fixed, nimag_fixed, psdindex_fixed, psdpars_fixed,  &
     wavs_var, nreal_var, nimag_var, psdindex_var, psdpars_var,  &
     PSD_Index, PSD_pars, n_real, n_imag )

!  Subroutine for using Fixed or variable-wavelength aerosol properties

   implicit none

!  inputs
!  ------

!  Control

   logical  , intent(in) :: do_Varprops, do_bimodal
   integer  , intent(in) :: Max_vpts, N_vpts
   real(fpk), intent(in) :: lambda

!  Fixed inputs

   integer  , intent(in) :: psdindex_fixed(2)
   real(fpk), intent(in) :: psdpars_fixed(3,2)
   real(fpk), intent(in) :: nreal_fixed(2), nimag_fixed(2)

!  Variable-property inputs

   integer  , intent(in) :: psdindex_var(Max_vpts,2)
   real(fpk), intent(in) :: psdpars_var(Max_vpts,3,2)
   real(fpk), intent(in) :: nreal_var(Max_vpts,2), nimag_var(Max_vpts,2)
   real(fpk), intent(in) :: wavs_var(Max_vpts)

!  Outputs
!  -------

   integer  , intent(out) :: PSD_index(2)
   real(fpk), intent(out) :: PSD_pars(3,2)
   real(fpk), intent(out) :: n_real(2), n_imag(2)

!  Local
!  -----

   logical   :: trawl
   integer   :: m, b, w, w1, w2, q
   real(fpk) :: f1, f2

!  Zero output

   n_real = 0.0_fpk ; n_imag   = 0.0_fpk
   PSD_Index = 0    ; PSD_pars = 0.0_fpk

   b = 1 ; if ( do_bimodal ) b = 2

!  If variable properties, trawl to get right wavelength
!    then interpolate (or extrapolate) linearly.

   if ( do_Varprops ) then
      w = 0 ; trawl = .true.
      do while ( trawl .and. w.lt. N_vpts )
         w = w + 1 ; trawl = ( lambda .gt. wavs_var(w))
      enddo
      w2 = w ; if ( w.eq.1 ) w2 = 2 ; w1 = w2 - 1
      f1 = ( wavs_var(w2) - lambda ) / ( wavs_var(w2) - wavs_var(w1) )
      f2 = 1.0_fpk - f1
      do m = 1, b
         PSD_Index(m)     = f1*psdindex_var(w1,m) + f2*psdindex_var(w2,m)
         do q = 1, 3
            PSD_pars(q,m) = f1*psdpars_var(w1,q,m) + f2*psdpars_var(w2,q,m)
         enddo
         n_real(m) = f1*nreal_var(w1,m) + f2*nreal_var(w2,m)
         n_imag(m) = f1*nimag_var(w1,m) + f2*nimag_var(w2,m)
      enddo
   endif

!  If fixed properties, copy input values directly (original default)

   if ( .not. do_Varprops ) then
      do m = 1, b
         PSD_Index(m)    = psdindex_fixed(m)
         PSD_pars(1:3,m) = psdpars_fixed(1:3,m)
         n_real(m) = nreal_fixed(m)
         n_imag(m) = nimag_fixed(m)
      enddo
   endif

!  Finish

   return
end subroutine get_variable_props

!  End module
 
End Module aerosol_regular_creation_m

