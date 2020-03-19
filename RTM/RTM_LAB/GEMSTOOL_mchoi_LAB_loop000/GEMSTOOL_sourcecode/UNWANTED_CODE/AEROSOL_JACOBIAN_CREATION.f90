Module aerosol_jacobian_creation_m

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
   use RTSMie_sourcecode_plus_m

!  Bimodal call

   use RTSMie_master_bimodal_m
   use RTSMie_master_bimodal_PLUS_m

!  Use Tmatrix code
!  ----------------

!  parameters

   use tmat_parameters

!  Single call

   use tmat_master_m
   use tmat_master_PLUS_m

!  Bimodal call

   use tmat_master_bimodal_m
   use tmat_master_bimodal_PLUS_m

!  Property call (@@ New, 11 Feb 2013 update)

   use aerosol_regular_creation_m, only : get_variable_props

!  Main routine is public

private
public :: aerosol_jacobian_creation

contains
!        PSD_Index, PSD_pars, n_real, n_imag, FixR1R2, R1, R2, epsnum, epsfac,   & ! INPUT, Mie/Tmat Parameters

subroutine aerosol_jacobian_creation &
      ( maxlayers, maxwav, maxaermoms, maxaerwfs, do_aer_jacobians_ref,                  & ! dimensioning
        do_aer_jacobians_psd, do_bimodal, Loading_case, TmatMie_case,                    & ! flags
        loading_upperboundary, loading_lowerboundary, aertau_input_w0, reference_w0,     & ! INPUT, Loading control
        exploading_relaxation, gdfloading_peakheight, gdfloading_halfwidth, LIDAR_file,  & ! INPUT, Loading control
        nlayers, nmuller, nwav, wav, height_grid, bimodal_fraction, momsize_cutoff,      & ! INPUT, Numbers, misc.
        do_Varprops, Max_Varpts, N_Varpts, Aerwavs_Var, FixR1R2, R1, R2,  & ! INPUT, Mie/Tmat Parameters            @@@ RTS 09/11/12
        PSDindex_fixed, PSDpars_fixed, nreal_fixed, nimag_fixed,          & ! INPUT, Mie/Tmat Parameters (fixed)    @@@ RTS 09/11/12
        PSDindex_var,   PSDpars_var,   nreal_var,   nimag_var,            & ! INPUT, Mie/Tmat Parameters (Variable) @@@ RTS 09/11/12
        R1R2_cutoff, nblocks, nweights, xparticle_limit,                                 & ! INPUT, Mie     PSD control
        do_EqSaSphere, spheroid_type, nkmax, ndgs, shape_factor, Tmat_accuracy,          & ! INPUT, Tmatrix PSD control 
        aerlayerflags, Loading, Dloading_Dtau, Dloading_Dpars,                           & ! OUTPUT, aerosol Loading
        n_aerosol_wfs, n_Aer_wfs_1, n_Aer_wfs_2, aerosol_pars, n_scatmoms,               & ! OUTPUT, aerosol WF control
        aertau_unscaled, aod_scaling, aerosol_deltau,                                    & ! OUTPUT, aerosol optical properties
        aerosol_ssalbs, aerosol_scatmoms, aerosol_distchars,                             & ! OUTPUT, aerosol optical properties
        L_aertau_unscaled, L_aod_scaling, L_aerosol_deltau,                              & ! OUTPUT, aerosol Linearized OPs
        L_aerosol_ssalbs, L_aerosol_scatmoms,                                            & ! OUTPUT, aerosol Linearized OPs
        fail1, fail2, Message_Loading, Messages_Optical )                                  ! Exception handling

!  =============================================================================
!                        AEROSOL JACOBIAN CREATION
!  =============================================================================

!  aerosol Loading:
!    Loading = optical depth profile
!    Loading Jacobians
!       Cases 1-3 : dloading_Dtau      = derivative of profile w.r.t the total aerosol optical depth at wavelength w0
!       Case 2    : dloading_Dpars(1)  = derivative of profile w.r.t the relaxation parameter (Exponential)
!       Case 3    : dloading_Dpars(1)  = derivative of profile w.r.t the GDF Peak height
!       Case 3    : dloading_Dpars(2)  = derivative of profile w.r.t the GDF Half Width 
!       Case 4    : dloading_Dpars     = 0 (No Derivatives allowed, LIDAR profile). Added 18 February 2013

!  Generation of Aerosol optical properties
!     1. Call to the Linearized Mie/Tmatrix program
!     2. Convert Mie/Tmatrix output (Microsopic) to IOP output (macroscopic) 

!  Here is a breakdown of the aerosol weighting functions
!    Bimodal, both Log-normal (2-parameter)
!   #               WEIGHTING FUNCTION

!   1         Real part of the refractive index for    Mode-1
!   2    Imaginary part of the refractive index for    Mode-1
!   3    Mode radius of Particle Size Distribution for Mode-1
!   4    Std. dev.   of Particle Size Distribution for Mode-1
!   5         Real part of the refractive index for    Mode-2
!   6    Imaginary part of the refractive index for    Mode-2
!   7    Mode radius of Particle Size Distribution for Mode-2
!   8    Std. dev.   of Particle Size Distribution for Mode-2
!   9    Bimodal fraction F, total = (1-f).Mode-1 + F.Mode-2
!   10   Total optical depth of aerosol (aerau_input_w0) at reference wavelength w0
!   11   EITHER Relaxation parameter in exponential loading; OR peakheight parameter in GDF loading
!   12   halfwidth parameter in GDF Loading

!  Here is a breakdown of the aerosol weighting functions
!    Bimodal, both MODIFIED GAMMA (3-parameter)

!   #               WEIGHTING FUNCTION

!   1         Real part of the refractive index for    Mode-1
!   2    Imaginary part of the refractive index for    Mode-1
!   3    Parameter 1 of Particle Size Distribution for Mode-1
!   4    Parameter 2 of Particle Size Distribution for Mode-1
!   5    Parameter 3 of Particle Size Distribution for Mode-1
!   6         Real part of the refractive index for    Mode-2
!   7    Imaginary part of the refractive index for    Mode-2
!   8    Parameter 1 of Particle Size Distribution for Mode-2
!   9    Parameter 2 of Particle Size Distribution for Mode-2
!   10   Parameter 3 of Particle Size Distribution for Mode-2
!   11   Bimodal fraction F, total = (1-f).Mode-1 + F.Mode-2
!   12   Total optical depth of aerosol (aerau_input_w0) at reference wavelength w0
!   13   EITHER Relaxation parameter in exponential loading; OR peakheight parameter in GDF loading
!   14   halfwidth parameter in GDF Loading

   implicit none

!  Input variables
!  ---------------

!  External Dimensioning

   integer, INTENT (IN) :: maxlayers, maxwav
   integer, INTENT (IN) :: maxaermoms, maxaerwfs

!  top level flags

   logical, INTENT (IN) :: do_bimodal
   logical, INTENT (IN) :: do_aer_Jacobians_ref
   logical, INTENT (IN) :: do_aer_Jacobians_psd

!  top level index (1 = Tmat, 2 = Mie)

   integer, INTENT(IN) :: TmatMie_case

!  Aerosol loading control
!  TOTAL AOT = aerosol optical thickness (value at reference wavelength )
!    Case 1 = uniform aerosol
!    Case 2 = Exponential Loading
!    Case 3 = GDF loading

   integer, INTENT(IN)       :: loading_case
   real(fpk)   , INTENT (IN) :: aertau_input_w0
   real(fpk)   , INTENT (IN) :: reference_w0

   real(fpk)   , INTENT (IN) :: exploading_relaxation
   real(fpk)   , INTENT (IN) :: gdfloading_peakheight
   real(fpk)   , INTENT (IN) :: gdfloading_halfwidth
   real(fpk)   , INTENT (IN) :: loading_upperboundary
   real(fpk)   , INTENT (IN) :: loading_lowerboundary

!  @@@ New inputs (18 February 2013)
!      Character String for file of LIDAR height/extinction profile data (Loading case 4)

   character*(*), INTENT(IN) :: LIDAR_file

!  Numbers

   integer, INTENT (IN) ::  nlayers
   integer, INTENT (IN) ::  nwav

!  Heights and wavelengths

   REAL    (fpk)   , INTENT (IN)   :: wav(maxwav)
   REAL    (fpk)   , INTENT (IN)   :: height_grid(0:maxlayers)

!    nmuller = 1 (scalar code), = 6 (vector code)

   integer, INTENT (IN) ::  nmuller


!  Bimodal fraction F (that is, total = F.Mie_1 + (1-F).Mie_2 )

   real(fpk)   , INTENT (IN) :: bimodal_fraction

!  Cutoff size for smallest phase function moment

   real(fpk)   , INTENT (IN) :: momsize_cutoff

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

   real    (fpk)   , intent(inout)  :: R1(2), R2(2)

!  N_REAL, N_IMAG - real and imaginary parts, refractive index (N-i.GE.0)

   real    (fpk),    intent(in)  :: nreal_fixed(2), nimag_fixed(2)                   ! Renamed, @@@ RTS 09/11/12
   real    (fpk),    intent(in)  :: nreal_Var(Max_Varpts,2), nimag_Var(Max_Varpts,2) ! New, @@@ RTS 09/11/12

!  Mie-specific inputs
!  ===================

!  Limiting particle size value. Set to 10000.0 default.
!   If you exceed this, program will tell you to increase dimensioning.

   REAL    (fpk)   , intent(inout)   :: xparticle_limit

!  R1R2_cutoff particle size for setting R1 and R2 internally

   REAL    (fpk)   , intent(in)  :: R1R2_cutoff(2)

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

    real(fpk)   , dimension ( maxlayers ),  INTENT (OUT)      :: Loading
    real(fpk)   , dimension ( maxlayers ),  INTENT (OUT)      :: DLoading_Dtau
    real(fpk)   , dimension ( maxlayers, 2 ),  INTENT (OUT)   :: Dloading_Dpars

!  aerosol layer flags

    LOGICAL,      DIMENSION ( maxlayers ), intent(out)  :: AERLAYERFLAGS 

!  General output
!    Total number of aerosol weighting functions, aerosol parameters
!    n_Aer_wfs = up to 5/6, from the linearized Mie/Tmatrix programs

   INTEGER         , INTENT (OUT)   :: n_aerosol_wfs, n_Aer_wfs_1, n_Aer_wfs_2
   REAL    (fpk)   , INTENT (OUT)   :: aerosol_pars(maxaerwfs)

!  AOP output: Number of exapnsion coefficients to be used

   integer, INTENT (OUT)   ::  n_scatmoms

!  Unscaled profiles of optical depth

    real(fpk)   , INTENT (OUT)  :: aertau_unscaled(maxlayers)
    real(fpk)   , INTENT (OUT)  :: L_aertau_unscaled(maxaerwfs,maxlayers)

!  AOD scaling

    real(fpk)   , INTENT (OUT) :: aod_scaling(maxwav), L_aod_scaling(maxwav,maxaerwfs)

!  Aod output, final

    REAL(fpk)   , DIMENSION( maxlayers, maxwav )           , INTENT (OUT)  :: AEROSOL_DELTAU
    REAL(fpk)   , DIMENSION( maxaerwfs, maxlayers, maxwav ), INTENT (OUT)  :: L_AEROSOL_DELTAU

!  AOPS

    REAL(fpk)   , DIMENSION( maxwav )                 , INTENT (OUT)  :: AEROSOL_SSALBS 
    REAL(fpk)   , DIMENSION( 6, 0:maxaermoms, maxwav ), INTENT (OUT)  :: AEROSOL_SCATMOMS

!  >>>>>>>>>>>> Added, 18 February 2013
!  Aerosol distribution characterisstics.
!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF
   real(fpk),     DIMENSION( 5, 2, maxwav)            , INTENT (OUT)  :: AEROSOL_DISTCHARS
!  >>>>>>>>>>>> Added, 18 February 2013

!  Linearized AOPs

    REAL(fpk)   , DIMENSION( maxaerwfs, maxwav )                  , INTENT (OUT) :: L_AEROSOL_SSALBS 
    REAL(fpk)   , DIMENSION( maxaerwfs, 6,  0:maxaermoms, maxwav ), INTENT (OUT) :: L_AEROSOL_SCATMOMS

!  Exception handling

   logical,        INTENT (OUT)           :: fail1, fail2
   character*(*),  INTENT (OUT)           :: message_Loading(3)
   character*(*),  INTENT (OUT)           :: Messages_Optical(5)

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

!  (Tmatrix only). Style flag.
!    * This is checked and re-set (if required) for Monodisperse case

   logical           :: Do_psd_OldStyle

!  linearization w.r.t refractive index components

   logical    :: Do_LinearRef

!  linearization w.r.t shape factor (T-matrix only)

   logical    :: Do_LinearEps

!  Linearization for PSD parameters

   logical    :: Do_LinearPSD

!  Mie code OUTPUT variables
!  =========================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk)    :: BMie_bulk (3)

!  linearizations w.r.t. PSD parameters

   real(fpk)    :: LPSD_BMie_bulk (3,3,2)

!  linearizations w.r.t. Refractive index parameters

   real(fpk)    :: LRFE_BMie_bulk (3,2,2)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

!  Regular quantities

   integer      :: BMie_ncoeffs
   real(fpk)    :: BMie_expcoeffs (6,0:max_Mie_angles)
   real(fpk)    :: BMie_asymm

!  linearizations w.r.t. PSD parameters

   real(fpk)     :: LPSD_BMie_expcoeffs (6,0:max_Mie_angles,3,2)
   real(fpk)     :: LPSD_BMie_asymm(3,2)

!  linearizations w.r.t. Refractive index parameters

   real(fpk)    :: LRFE_BMie_expcoeffs (6,0:max_Mie_angles,2,2)
   real(fpk)    :: LRFE_BMie_asymm(2,2)

!  F-matrix,  optional output
!  --------------------------

!  F-matrix

   real(fpk)    :: BMie_Fmatrix(4,max_Mie_angles)

!  Linearizations of F-matrix

   real(fpk)    :: LPSD_BMie_Fmatrix(4,max_Mie_angles,3,2)
   real(fpk)    :: LRFE_BMie_Fmatrix(4,max_Mie_angles,2,2)

!  Fraction Jacobian
!  -----------------

   real(fpk)     :: LFRC_BMie_bulk (3)
   real(fpk)     :: LFRC_BMie_expcoeffs (6,0:max_Mie_angles)
   real(fpk)     :: LFRC_BMie_asymm
   real(fpk)     :: LFRC_BMie_Fmatrix (4,max_Mie_angles)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(fpk)     :: BMie_dist (5,2)
   real(fpk)     :: LPSD_BMie_dist (5,3,2)

!  Exception handling
!  ------------------

   character*90  :: Mie_Bmessages(3)
   character*90  :: Mie_trace_3

!  Tmatrix code OUTPUT variables (Notation similar)
!  =============================

!  Bulk distribution parameters

   real(fpk)    :: BTmat_bulk (3)

!  linearizations w.r.t. PSD parameters

   real(fpk)    :: LPSD_BTmat_bulk (3,3,2)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(fpk)    :: LRFE_BTmat_bulk (3,3,2)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

!  Regular quantities

   integer       :: BTmat_ncoeffs
   real(fpk)     :: BTmat_expcoeffs (NPL1,6)
   real(fpk)     :: BTmat_asymm

!  linearizations w.r.t. PSD parameters

   real(fpk)     :: LPSD_BTmat_expcoeffs (NPL1,6,3,2)
   real(fpk)     :: LPSD_BTmat_asymm(3,2)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(fpk)     :: LRFE_BTmat_expcoeffs (NPL1,6,3,2)
   real(fpk)     :: LRFE_BTmat_asymm(3,2)

!  F-matrix,  optional output
!  --------------------------

!  F-matrix

   real(fpk)     :: BTmat_Fmatrix (MAXNPA,6)

!  Linearizations of F-matrix

   real(fpk)     :: LPSD_BTmat_Fmatrix (MAXNPA,6,3,2)
   real(fpk)     :: LRFE_BTmat_Fmatrix (MAXNPA,6,3,2)

!  Fraction Jacobian
!  -----------------

   real(fpk)     :: LFRC_BTmat_bulk (3)
   real(fpk)     :: LFRC_BTmat_expcoeffs (NPL1,6)
   real(fpk)     :: LFRC_BTmat_asymm
   real(fpk)     :: LFRC_BTmat_Fmatrix (MAXNPA,6)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(fpk)     :: BTmat_dist (5,2)
   real(fpk)     :: LPSD_BTmat_dist (5,3,2)

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

    real(fpk)     :: extinction_ref, L_extinction_ref(maxaerwfs)

!  Local mask

   integer :: distpars_mask(9), npsd(2)
   data distpars_mask / 2,2,3,2,3,3,3,3,0 /

!  Help Variables

   character*10 :: ctype
   character*3  :: cwav
   integer      :: istatus, NLIN
   integer      :: n, k, l, m, kd, q, qm, qm1, qm2, qmt, lambda_index, w, wc, wavmask(maxwav), n_scatmoms_w
   real(fpk)    :: wavelength, Mie_wavelength, Tmat_wavelength, extinction, L_extinction(maxaerwfs)
   logical      :: do_aer_Jacobians, Local_do_Varprops             ! New @@@ RTS 09/11/12

   real(fpk)    :: normfac(maxaerwfs)

!  Local control

!   logical, parameter :: do_iopchecker   = .true.
   logical, parameter :: do_iopchecker   = .false.

!  Initialize output
!  =================

!  Initialize exception handling

   fail1 = .false.
   fail2 = .false.
   Message_Loading  = ' '
   Messages_optical = ' '

!  initialize Aerosol Loading

   Loading        = 0.0d0
   Dloading_Dtau  = 0.0d0
   Dloading_Dpars = 0.0d0
   aerlayerflags = .false.

!  Set local Mie/Tmatrix variables

   do_monodisperse = .false.             ! Always false
   Do_Fmatrix      = .false.             ! Always false
   DO_Expcoeffs    = .false.             ! this will be set later on
   monoradius      = 0.0d0
   do_psd_oldstyle = .false.

!  Set local Mie/Tmatrix Linearization variables

   do_LinearRef = do_aer_Jacobians_Ref
   do_LinearPSD = do_aer_Jacobians_PSD
   do_LinearEps = (do_aer_Jacobians_Ref .and. TmatMie_case .eq. 1)

!  overall Jacobian flag

   do_aer_Jacobians = do_aer_Jacobians_Ref .or.do_aer_Jacobians_PSD

!  Initialize optical properties

   extinction_ref    = 0.0d0
   aod_scaling       = 0.0d0
   aerosol_deltau    = 0.0d0
   aerosol_ssalbs    = 0.0d0
   aerosol_scatmoms  = 0.0d0
   aerosol_distchars = 0.0d0

   L_aod_scaling      = 0.0d0
   L_aerosol_deltau   = 0.0d0
   L_aerosol_ssalbs   = 0.0d0
   L_aerosol_scatmoms = 0.0d0

!  Initialize local aerosol properties (Just a precaution)! New, @@@ RTS 09/11/12

   n_real = 0.0_fpk ; n_imag   = 0.0_fpk
   PSD_Index = 0    ; PSD_pars = 0.0_fpk

!  Important Check
!  ===============

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

!  Linearization with wavelength-variable microscopic properties is NOT ALLOWED
!    Otherwise bookkeeping gets reallyl complicated, requiring many changes (11 Feb 2013)

   if ( do_aer_Jacobians .and. Do_varprops ) then
      fail2 = .true.
      Messages_Optical(1) = 'Jacobian Calculation not allowed for the following condition'
      Messages_Optical(2) = '  * Do_varprops is true '
      Messages_Optical(3) = '  * Do_varprops is true '
      Messages_Optical(4) = 'Action: Turn off var-props Flag in aerosol inputs'
      Messages_optical(5) = 'Trace : Initial check inin AEROSOL_JACOBIAN_CREATION'
      return
   endif

!  LATER - We WANT TO INTRODUCE ANOTHER CONDITION..........(only when we are sure)
 !  Linearization with wavelength-variable microscopic properties is only possible
!   if the Reference Wavelength is one of the list of established wavelengths.
!    That is, Do_varprops and lambda_index = 0 cannot both be true.

!   if ( do_aer_Jacobians .and. (Do_varprops. and. lambda_index.eq.0) ) then
!      fail2 = .true.
!      Messages_Optical(1) = 'Jacobian Calculation not allowed for the following condition'
!      Messages_Optical(2) = '  * Do_varprops. and. lambda_index.eq.0 '
!      Messages_Optical(3) = '  * That is, Reference wavelength must be one of the output wavelengths'
!      Messages_Optical(4) = 'Action: Set reference wavelength one of the outputs'
!      Messages_optical(5) = 'Trace : Initial check inin AEROSOL_JACOBIAN_CREATION'
!      return
!   endif

!  Now form the aerosol Loading
!  ----------------------------

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

!  Bookkeeping
!  @@@@@@@@@@@

!  Initialize number of weighting functions

   n_Aer_wfs_1 = 0
   n_Aer_wfs_2 = 0

!  Methodology, 11 february 2013. No variable properties for Jacobians.
!     This way, we can do the bookkeeping once and for all, here.

!  Get properties (copied from "fixed" quantities)
!  Call to assign local aerosol properties ! New, @@@ RTS 09/11/12, 02/11/13

   Local_do_Varprops = .false.
   call get_variable_props &
        ( local_do_Varprops, do_bimodal, Max_Varpts, N_Varpts, wav(1),  &
          nreal_fixed, nimag_fixed, psdindex_fixed, psdpars_fixed,      &
          Aerwavs_Var, nreal_var, nimag_var, psdindex_var, psdpars_var, &
          PSD_Index, PSD_pars, n_real, n_imag )

!  No bookkeeping if no WFs

   if ( .not. do_aer_Jacobians ) go to 433

!  Set up aerosol parameters and total number of aerosol WFs, First aerosol mode

   npsd(1) = distpars_mask(PSD_index(1))
   m = 0
   if ( do_LinearRef ) then
      aerosol_pars(1) = n_real(1) ; aerosol_pars(2) = n_imag(1)
      m = m + 2
   endif
   if ( do_LinearEps ) then
       aerosol_pars(3) = Shape_factor(1) ; m = m + 1
   endif
   if ( do_LinearPSD ) then
      do kd = 1, npsd(1)
         m = m + 1 ; aerosol_pars(m) = psd_pars(kd,1)
      end do
   endif
   n_Aer_wfs_1 = m

!  Set up aerosol parameters and total number of aerosol WFs, Second aerosol mode

   if ( do_bimodal ) then
      npsd(2) = distpars_mask(PSD_index(2))
      m = n_Aer_wfs_1
      if ( do_LinearRef ) then
         aerosol_pars(m+1) = n_real(2) ; aerosol_pars(m+2) = n_imag(2)
         m = m + 2
      endif
      if ( do_LinearEps ) then
          aerosol_pars(m+1) = Shape_factor(2) ; m = m + 1
      endif
      if ( do_LinearPSD ) then
         do kd = 1, npsd(2)
            m = m + 1 ; aerosol_pars(m) = psd_pars(kd,2)
         end do
      endif
      n_Aer_wfs_2 = m - n_Aer_wfs_1
      aerosol_pars(m+1) = bimodal_fraction
   endif

!  Total number of WFs from the Mie/Tmatrix parts

    qm1 = n_Aer_wfs_1
    if ( do_bimodal )  qm1 = qm1 + n_Aer_wfs_2 + 1

!  Add WFs for the loading
!    Total loading + relaxation for the Exponential profile
!    Total loading + peak height + half width for the GDF profile
!    Total loading only for the Uniform profile.

    if ( loading_case .eq. 2) then
      qmt = qm1 + 2
    else if ( loading_case .eq. 3 ) then
      qmt = qm1 + 3
    else
      qmt = qm1 + 1
    endif

!  total number of weighing functions

   n_aerosol_wfs = qmt

!  Add 1 or 2-3 more parameters for loading 

   qm2 = qm1 + 1
   aerosol_pars(qm2) = aertau_input_w0
   if ( loading_case .eq. 2 ) then
     aerosol_pars(qmt) = exploading_relaxation
   else if ( loading_case .eq. 3 ) then
     aerosol_pars(qm2+1) = gdfloading_peakheight
     aerosol_pars(qmt)   = gdfloading_halfwidth
   endif

!  debug Loading derivatives. OK - 26 June 2012
!     do n = 91, 113, 10
!     write(*,466)n,Loading(n),Dloading_Dtau(n)*aertau_input_w0,&
!                 Dloading_Dpars(n,1)*gdfloading_peakheight,&
!                 Dloading_Dpars(n,2)*gdfloading_halfwidth
!     enddo
!466 format(i4,1p4e20.10)

!  Continuation point for skipping Jacobian bookkeeping

433 continue

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

!  Only require extinction coefficient and its linearization if flagged
!  Set the local Mie program inputs (bulk properties only)

      Do_Expcoeffs     = .FALSE.

!  reference wavelength

      wavelength     = reference_w0
      Mie_wavelength = wavelength/1000.0d0

!  Local properties. Limited to the "Fixed option" here. 02/11/13

      Local_do_Varprops = .false.
      call get_variable_props &
        ( local_do_Varprops, do_bimodal, Max_Varpts, N_Varpts, wavelength,    &
          nreal_fixed, nimag_fixed, psdindex_fixed, psdpars_fixed,      &
          Aerwavs_Var, nreal_var, nimag_var, psdindex_var, psdpars_var, &
          PSD_Index, PSD_pars, n_real, n_imag )

!  progress

      if ( do_aer_jacobians ) then
         write(*,*)'Linearized Mie Aerosol calculation, reference wavelength # '
      else 
         write(*,*)'Regular Mie Aerosol calculation, reference wavelength # '
      endif

!  BIMODAL CALL

      if ( do_bimodal) then
         if ( do_aer_Jacobians ) then
            call RTSMie_master_bimodal_plus                       &
            ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
              Do_LinearRef, Do_LinearPSD,                         & ! I
              PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,   & ! I
              nblocks, nweights, xparticle_limit, R1R2_cutoff,    & ! I
              n_Fmatrix_angles, Fmatrix_angles,                   & ! I
              Mie_wavelength, n_real, n_imag, bimodal_fraction,   & ! I
              BMie_bulk, BMie_asymm, BMie_ncoeffs,      & ! O
              BMie_expcoeffs, BMie_Fmatrix,             & ! O
              LPSD_BMie_bulk, LPSD_BMie_asymm,          & ! O
              LPSD_BMie_expcoeffs, LPSD_BMie_Fmatrix,   & ! O
              LRFE_BMie_bulk, LRFE_BMie_asymm,          & ! O
              LRFE_BMie_expcoeffs, LRFE_BMie_Fmatrix,   & ! O
              LFRC_BMie_bulk, LFRC_BMie_asymm,          & ! O
              LFRC_BMie_expcoeffs, LFRC_BMie_Fmatrix,   & ! O
              BMie_dist, LPSD_BMie_dist,                & ! O
              fail2, istatus, Mie_Bmessages, Mie_trace_3 )  ! O
         else
            call RTSMie_master_bimodal                            &
            ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
              PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,   & ! I
              nblocks, nweights, xparticle_limit, R1R2_cutoff,    & ! I
              n_Fmatrix_angles, Fmatrix_angles,                   & ! I
              Mie_wavelength, n_real, n_imag, bimodal_fraction,   & ! I
              BMie_bulk, BMie_asymm, BMie_ncoeffs,                & ! O
              BMie_expcoeffs, BMie_Fmatrix, BMie_dist,            & ! O
              fail2, istatus, Mie_Bmessages, Mie_trace_3 )          ! O
         endif
      endif

!  SINGLE CALL

      if ( .not. do_bimodal ) then
         BMie_dist(:,2) = 0.0d0
         if ( do_aer_Jacobians ) then
            call RTSMie_main_plus                                                               &  !---MIE CALL
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,  Do_LinearRef, Do_LinearPSD,         & ! I
                PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),              & ! I
                nblocks(1), nweights(1), xparticle_limit, R1R2_cutoff(1),                       & ! I
                n_Fmatrix_angles, Fmatrix_angles, Mie_wavelength, n_real(1), n_imag(1),         & ! I
                BMie_bulk, BMie_asymm, BMie_ncoeffs, BMie_expcoeffs, BMie_Fmatrix, BMie_dist(:,1),   & ! O
                LPSD_BMie_bulk(:,:,1), LPSD_BMie_asymm(:,1),  LPSD_BMie_expcoeffs(:,:,:,1),     & ! O
                LPSD_BMie_Fmatrix(:,:,:,1), LPSD_BMie_dist(:,:,1), LRFE_BMie_bulk(:,:,1),       & ! O
                LRFE_BMie_asymm(:,1), LRFE_BMie_expcoeffs(:,:,:,1), LRFE_BMie_Fmatrix(:,:,:,1), & ! O
                fail2, istatus, Mie_Bmessages(1), Mie_Bmessages(2), Mie_Bmessages(3) )            ! O
         else
            call RTSMie_main                                                                   & !---MIE CALL
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,                                     & ! I
                PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),             & ! I
                nblocks(1), nweights(1), xparticle_limit, R1R2_cutoff(1),                      & ! I
                n_Fmatrix_angles, Fmatrix_angles, Mie_wavelength, n_real(1), n_imag(1),        & ! I
                BMie_bulk, BMie_asymm, BMie_ncoeffs, BMie_expcoeffs, BMie_Fmatrix, BMie_dist(:,1),  & ! O
                fail2, istatus, Mie_Bmessages(1), Mie_Bmessages(2), Mie_Bmessages(3) )           ! O
         endif
      endif

!  Exception handling on everything

      if ( Fail2 ) then  
         do m = 1, 3   
            Messages_Optical(m) = adjustl(trim(Mie_Bmessages(m)))
         enddo
         Messages_Optical(4) = 'Single Mie Call'; if ( do_bimodal ) Messages_Optical(4) = adjustl(trim(Mie_trace_3))
         ctype = 'Regular   ' ; if ( do_aer_Jacobians ) ctype = 'Linearized'
         Messages_optical(5) = 'First call to the '//ctype//' Mie program in AEROSOL CREATION, reference wavelength'
         return
      endif

!  Set the reference quantities. NOTE THE Bookkeeping.

      extinction_ref = BMie_bulk(1) ; qm  = 0
      if ( do_aer_Jacobians ) then
         if ( Do_LinearRef) then
            L_extinction_ref(qm+1:qm+2) = LRFE_BMie_bulk(1,1:2,1) ; qm = qm + 2
         endif
         if ( Do_LinearPSD) then
            L_extinction_ref(qm+1:qm+npsd(1)) = LPSD_BMie_bulk(1,1:npsd(1),1) ; qm = qm + npsd(1)
         endif
         if ( do_bimodal ) then
            if ( Do_LinearRef) then
               L_extinction_ref(qm+1:qm+2) = LRFE_BMie_bulk(1,1:2,2) ; qm = qm + 2
            endif
            if ( Do_LinearPSD) then
               L_extinction_ref(qm+1:qm+npsd(2)) = LPSD_BMie_bulk(1,1:npsd(2),2) ; qm = qm + npsd(2)
            endif
            L_extinction_ref(qm+1) = LFRC_BMie_bulk(1)
         endif
      endif

!      do q = 1, 9
!        write(*,*)'Ext Ref',q,extinction_ref, L_extinction_ref(q)*aerosol_pars(q)
!      enddo

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

!  Local properties. Limited to the "Fixed option" here. 02/11/13

         Local_do_Varprops = .false.
         call get_variable_props &
           ( local_do_Varprops, do_bimodal, Max_Varpts, N_Varpts, wavelength,    &
             nreal_fixed, nimag_fixed, psdindex_fixed, psdpars_fixed,      &
             Aerwavs_Var, nreal_var, nimag_var, psdindex_var, psdpars_var, &
             PSD_Index, PSD_pars, n_real, n_imag )

!  progress

         if ( do_aer_jacobians ) then
            write(*,*)'Linearized Mie Aerosol calculation, doing wavelength # ',wc
         else 
            write(*,*)'Regular Mie Aerosol calculation, doing wavelength # ',wc
         endif

!  wavelength for the Mie code (micron unit)  

         Mie_wavelength = wavelength/1000.0d0

!  First Call to develop aerosol Mie properties
!  --------------------------------------------

!  Set the local Mie program inputs (general)

         Do_Expcoeffs     = .TRUE.

!  BIMODAL CALL

         write(*,*)'do_bimodal',do_bimodal
         if ( do_bimodal) then
          if ( do_aer_Jacobians ) then
            call RTSMie_master_bimodal_plus                       &
            ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
              Do_LinearRef, Do_LinearPSD,                         & ! I
              PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,   & ! I
              nblocks, nweights, xparticle_limit, R1R2_cutoff,    & ! I
              n_Fmatrix_angles, Fmatrix_angles,                   & ! I
              Mie_wavelength, n_real, n_imag, bimodal_fraction,   & ! I
              BMie_bulk, BMie_asymm, BMie_ncoeffs,      & ! O
              BMie_expcoeffs, BMie_Fmatrix,             & ! O
              LPSD_BMie_bulk, LPSD_BMie_asymm,          & ! O
              LPSD_BMie_expcoeffs, LPSD_BMie_Fmatrix,   & ! O
              LRFE_BMie_bulk, LRFE_BMie_asymm,          & ! O
              LRFE_BMie_expcoeffs, LRFE_BMie_Fmatrix,   & ! O
              LFRC_BMie_bulk, LFRC_BMie_asymm,          & ! O
              LFRC_BMie_expcoeffs, LFRC_BMie_Fmatrix,   & ! O
              BMie_dist, LPSD_BMie_dist,                & ! O
              fail2, istatus, Mie_Bmessages, Mie_trace_3 )  ! O
          else
            call RTSMie_master_bimodal                            &
            ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
              PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,   & ! I
              nblocks, nweights, xparticle_limit, R1R2_cutoff,    & ! I
              n_Fmatrix_angles, Fmatrix_angles,                   & ! I
              Mie_wavelength, n_real, n_imag, bimodal_fraction,   & ! I
              BMie_bulk, BMie_asymm, BMie_ncoeffs,                & ! O
              BMie_expcoeffs, BMie_Fmatrix, BMie_dist,            & ! O
              fail2, istatus, Mie_Bmessages, Mie_trace_3 )          ! O
          endif
        endif

!  SINGLE CALL

        if ( .not. do_bimodal ) then
          BMie_dist(:,2) = 0.0d0
          if ( do_aer_Jacobians ) then
            call RTSMie_main_plus                                                               &  !---MIE CALL
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,  Do_LinearRef, Do_LinearPSD,         & ! I
                PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),              & ! I
                nblocks(1), nweights(1), xparticle_limit, R1R2_cutoff(1),                       & ! I
                n_Fmatrix_angles, Fmatrix_angles, Mie_wavelength, n_real(1), n_imag(1),         & ! I
                BMie_bulk, BMie_asymm, BMie_ncoeffs, BMie_expcoeffs, BMie_Fmatrix, BMie_dist(:,1),   & ! O
                LPSD_BMie_bulk(:,:,1), LPSD_BMie_asymm(:,1),  LPSD_BMie_expcoeffs(:,:,:,1),     & ! O
                LPSD_BMie_Fmatrix(:,:,:,1), LPSD_BMie_dist(:,:,1), LRFE_BMie_bulk(:,:,1),        & ! O
                LRFE_BMie_asymm(:,1), LRFE_BMie_expcoeffs(:,:,:,1), LRFE_BMie_Fmatrix(:,:,:,1),  & ! O
                fail2, istatus, Mie_Bmessages(1), Mie_Bmessages(2), Mie_Bmessages(3) )           ! O
          else
            call RTSMie_main                                                                   & !---MIE CALL
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,                                     & ! I
                PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),             & ! I
                nblocks(1), nweights(1), xparticle_limit, R1R2_cutoff(1),                      & ! I
                n_Fmatrix_angles, Fmatrix_angles, Mie_wavelength, n_real(1), n_imag(1),        & ! I
                BMie_bulk, BMie_asymm, BMie_ncoeffs, BMie_expcoeffs, BMie_Fmatrix, BMie_dist(:,1),  & ! O
                fail2, istatus, Mie_Bmessages(1), Mie_Bmessages(2), Mie_Bmessages(3) )           ! O
          endif
        endif

!  Exception handling on everything

        if ( Fail2 ) then
          write(cwav,'(I3)')wc
          do m = 1, 3
            Messages_Optical(m) = adjustl(trim(Mie_Bmessages(m)))
          enddo
          Messages_Optical(4) = 'Single Mie Call'; if ( do_bimodal ) Messages_Optical(4) = adjustl(trim(Mie_trace_3))
          ctype = 'Regular   ' ; if ( do_aer_Jacobians ) ctype = 'Linearized'
          Messages_Optical(5) = 'First call to the '//ctype//' Mie program in AEROSOL CREATION, wavelength # '//cwav
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
            extinction_ref = BMie_bulk(1) ; qm  = 0
            if ( do_aer_Jacobians ) then
               if ( Do_LinearRef) then
                  L_extinction_ref(qm+1:qm+2) = LRFE_BMie_bulk(1,1:2,1) ; qm = qm + 2
               endif
               if ( Do_LinearPSD) then
                  L_extinction_ref(qm+1:qm+npsd(1)) = LPSD_BMie_bulk(1,1:npsd(1),1) ; qm = qm + npsd(1)
               endif
               if ( do_bimodal ) then
                  if ( Do_LinearRef) then
                     L_extinction_ref(qm+1:qm+2) = LRFE_BMie_bulk(1,1:2,2) ; qm = qm + 2
                  endif
                  if ( Do_LinearPSD) then
                     L_extinction_ref(qm+1:qm+npsd(2)) = LPSD_BMie_bulk(1,1:npsd(2),2) ; qm = qm + npsd(2)
                  endif
                  L_extinction_ref(qm+1) = LFRC_BMie_bulk(1)
               endif
            endif
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

!  Convert Mie derivatives to IOP derivatives
!  ==========================================

         if ( do_aer_Jacobians ) then

!  Linearized extinction

            qm = 0
            if ( Do_LinearRef) then
               L_extinction(qm+1:qm+2) = LRFE_BMie_bulk(1,1:2,1) ; qm = qm + 2
            endif
            if ( Do_LinearPSD) then
               L_extinction(qm+1:qm+npsd(1)) = LPSD_BMie_bulk(1,1:npsd(1),1) ; qm = qm + npsd(1)
            endif
            if ( do_bimodal ) then
               if ( Do_LinearRef) then
                  L_extinction(qm+1:qm+2) = LRFE_BMie_bulk(1,1:2,2) ; qm = qm + 2
               endif
               if ( Do_LinearPSD) then
                  L_extinction(qm+1:qm+npsd(2)) = LPSD_BMie_bulk(1,1:npsd(2),2) ; qm = qm + npsd(2)
               endif
               L_extinction(qm+1) = LFRC_BMie_bulk(1)
            endif

!  Linearization of the scaling factor

            do qm = 1, qm1
              L_aod_scaling(w,qm) = aod_scaling(w) * ( (L_extinction(qm)/extinction)        &
                                           - (L_extinction_ref(qm)/extinction_ref) )
            enddo
            qm2 = qm1 + 1
            L_aod_scaling(w,qm2) = aod_scaling(w)
            if ( loading_case .eq. 2 ) then
               L_aod_scaling(w,qmt) = aod_scaling(w)
            else if ( loading_case .eq. 3 ) then
               L_aod_scaling(w,qm2+1) = aod_scaling(w)
               L_aod_scaling(w,qmt)   = aod_scaling(w)
            endif

            qm = 0
            if ( Do_LinearRef) then
               do l = 0, n_scatmoms_w
                  do k = 1, nmuller
                     L_aerosol_scatmoms(qm+1:qm+2,k,l,w) = LRFE_BMie_expcoeffs(k,l,1:2,1)
                  enddo
               enddo
               L_aerosol_ssalbs(qm+1:qm+2,w) = LRFE_BMie_bulk(3,1:2,1) ; qm = qm + 2
            endif
            if ( Do_LinearPSD) then
               do l = 0, n_scatmoms_w
                  do k = 1, nmuller
                     L_aerosol_scatmoms(qm+1:qm+npsd(1),k,l,w) = LPSD_BMie_expcoeffs(k,l,1:npsd(1),1)
                  enddo
               enddo
               L_aerosol_ssalbs(qm+1:qm+npsd(1),w) = LPSD_BMie_bulk(3,1:npsd(1),1) ; qm = qm + npsd(1)
            endif

            if ( do_bimodal ) then
               if ( Do_LinearRef) then
                  do l = 0, n_scatmoms_w
                     do k = 1, nmuller
                        L_aerosol_scatmoms(qm+1:qm+2,k,l,w) = LRFE_BMie_expcoeffs(k,l,1:2,2)
                     enddo
                  enddo
                  L_aerosol_ssalbs(qm+1:qm+2,w) = LRFE_BMie_bulk(3,1:2,2) ; qm = qm + 2
               endif
               if ( Do_LinearPSD) then
                  do l = 0, n_scatmoms_w
                     do k = 1, nmuller
                        L_aerosol_scatmoms(qm+1:qm+npsd(2),k,l,w) = LPSD_BMie_expcoeffs(k,l,1:npsd(2),2)
                     enddo
                  enddo
                  L_aerosol_ssalbs(qm+1:qm+npsd(2),w) = LPSD_BMie_bulk(3,1:npsd(2),2) ; qm = qm + npsd(2)
               endif
               do l = 0, n_scatmoms_w
                  do k = 1, nmuller
                     L_aerosol_scatmoms(qm1,k,l,w) = LFRC_BMie_expcoeffs(k,l)
                  enddo
               enddo
               L_aerosol_ssalbs(qm1,w) = LFRC_BMie_bulk(3)
            endif

!  Phase function first moment

            do qm = 1, qm1
               L_aerosol_scatmoms(qm,1,0,w) = 0.0d0
            enddo

!  End Jacobians clause

         endif

!  Update n_scatmoms

         n_scatmoms = min(n_scatmoms, n_scatmoms_w)

!  Assign Unscaled loadings + linearizations. VERSION 2

         aertau_unscaled = 0.0d0 ; L_aertau_unscaled = 0.0d0
         do n = 1, nlayers
            aerlayerflags(n) =  ( Loading(n) .ne. 0.0d0  )
            aertau_unscaled(n) = loading(n)
         enddo
         if ( do_aer_Jacobians ) then
            do n = 1, nlayers
               aertau_unscaled(n) = loading(n)
               do q = 1, qm1
                  L_aertau_unscaled(q,n) = loading(n)
               enddo
               L_aertau_unscaled(qm1+1,n) = Dloading_Dtau(n)
               if ( loading_case .eq. 2 ) then
                  L_aertau_unscaled(qm1+2,n) = Dloading_Dpars(n,1)
               else if ( loading_case .eq. 3 ) then
                  L_aertau_unscaled(qm1+2,n) = Dloading_Dpars(n,1)
                  L_aertau_unscaled(qm1+3,n) = Dloading_Dpars(n,2)
               endif
            enddo
         endif

!  Apply scalings to loadings and linearizations

         do n = 1, nlayers
            aerosol_deltau(n,w) = aertau_unscaled(n) * aod_scaling(w)
         enddo
         if ( do_aer_jacobians ) then
            do n = 1, nlayers
               do q = 1, n_aerosol_wfs
                  L_aerosol_deltau(q,n,w) = L_aertau_unscaled(q,n) * L_aod_scaling(w,q)
               enddo
            enddo
         endif

!  debug aerosol optical properties. VERSION TWO only

         if ( do_iopchecker ) then
          if ( do_aer_jacobians ) then
           do n = 1, nlayers
            if (aerlayerflags(n).and.n.eq.107 ) then
              write(777,'(i4,1p6e20.10)')n,aerosol_deltau(n,w),aerosol_ssalbs(w)
              write(777,'(i4,1p11e20.10)')n,(L_aerosol_deltau(q,n,w),q=1,11)
              write(777,'(i4,1p11e20.10)')n,(L_aerosol_ssalbs(q,w),q=1,11)
              do l = 0, n_scatmoms_w
                write(777,'(2i5,1pe20.10,1p11e15.6)')n,l,aerosol_scatmoms(1,l,w),(l_aerosol_scatmoms(q,1,l,w),q=1,11)
              enddo
            endif
           enddo
!           pause'Lin 777'
          else
           do n = 1, nlayers
            if (aerlayerflags(N).and.n.eq.107 ) then
              write(888,'(i4,1p6e20.10)')n,aerosol_deltau(n,w),aerosol_ssalbs(w)
              do l = 0, n_scatmoms_w
                write(888,'(2i5,1p6e20.10)')n,l,(aerosol_scatmoms(k,l,w),k=1,1)
              enddo
!            else
!              write(888,'(i4,1p6e20.10)')n,aerosol_deltau(n,w)
            endif
           enddo
!           pause'Reg 888'
          endif
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

!  Prepare the reference-wavelength Tmatrix Inputs
!  ===============================================

   if ( TmatMie_case .eq. 2 .and. lambda_index .eq. 0 ) then

!  Only require extinction coefficient and its linearization if flagged
!  Set the local Tmat program inputs (bulk properties only)

      Do_Expcoeffs     = .FALSE.

!  reference wavelength

      wavelength      = reference_w0
      Tmat_wavelength = wavelength/1000.0d0

!  Local properties. Limited to the "Fixed option" here. 02/11/13

      Local_do_Varprops = .false.
      call get_variable_props &
           ( local_do_Varprops, do_bimodal, Max_Varpts, N_Varpts, wavelength,    &
             nreal_fixed, nimag_fixed, psdindex_fixed, psdpars_fixed,      &
             Aerwavs_Var, nreal_var, nimag_var, psdindex_var, psdpars_var, &
             PSD_Index, PSD_pars, n_real, n_imag )

!  progress

      if ( do_aer_jacobians ) then
         write(*,*)'Linearized Tmat Aerosol calculation, reference wavelength # '
      else 
         write(*,*)'   Regular Tmat Aerosol calculation, reference wavelength # '
      endif

!  BIMODAL CALL

      if ( do_bimodal) then
         if ( do_aer_Jacobians ) then
            call tmat_master_bimodal_PLUS                       &
               ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere,              & ! Inputs (Flags)
                 Do_LinearRef, Do_LinearEps, Do_LinearPSD,                              & ! Inputs (Flags)
                 Do_psd_OldStyle, psd_Index, psd_pars, MonoRadius, R1, R2, FixR1R2,     & ! Inputs (PSD)
                 bimodal_fraction, Spheroid_type, nkmax, n_Fmatrix_angles, ndgs,        & ! Inputs (General)
                 shape_factor, Tmat_accuracy, Tmat_wavelength, n_real, n_imag,          & ! Inputs (Optical)
                 BTmat_bulk, BTmat_asymm, BTmat_ncoeffs,     & ! O
                 BTmat_expcoeffs, BTmat_Fmatrix,             & ! O
                 LPSD_BTmat_bulk, LPSD_BTmat_asymm,          & ! O
                 LPSD_BTmat_expcoeffs, LPSD_BTmat_Fmatrix,   & ! O
                 LRFE_BTmat_bulk, LRFE_BTmat_asymm,          & ! O
                 LRFE_BTmat_expcoeffs, LRFE_BTmat_Fmatrix,   & ! O
                 LFRC_BTmat_bulk, LFRC_BTmat_asymm,          & ! O
                 LFRC_BTmat_expcoeffs, LFRC_BTmat_Fmatrix,   & ! O
                 NLIN, BTmat_dist, LPSD_BTmat_dist,          & ! O
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2, Tmat_trace_3 )   ! Outputs (status)
         else
            call tmat_master_bimodal                 &
               ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere,              & ! Inputs (Flags)
                 Do_psd_OldStyle, psd_Index, psd_pars, MonoRadius, R1, R2, FixR1R2,     & ! Inputs (PSD)
                 bimodal_fraction, Spheroid_type, nkmax, n_Fmatrix_angles, ndgs,        & ! Inputs (General)
                 shape_factor, Tmat_accuracy, Tmat_wavelength, n_real, n_imag,          & ! Inputs (Optical)
                 Btmat_bulk, Btmat_asymm, Btmat_ncoeffs,                                & ! Outputs (Tmat)
                 Btmat_expcoeffs, Btmat_Fmatrix, Btmat_dist,                            & ! Outputs (PSD)
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2, Tmat_trace_3 )   ! Outputs (status)
         endif
      endif

!  SINGLE CALL

      if ( .not. do_bimodal ) then
         Btmat_dist(:,2) = 0.0d0
         if ( do_aer_Jacobians ) then
            call tmat_master_PLUS                    &
              ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere,                  & ! Inputs (Flags) 
                Do_LinearRef, Do_LinearEps, Do_LinearPSD, Do_psd_OldStyle,                 & ! Inputs (Flags) 
                 PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),        & ! Inputs (PSD)
                 Spheroid_type, nkmax(1), n_Fmatrix_angles, ndgs(1),                       & ! Inputs (General)
                 shape_factor(1), Tmat_accuracy, Tmat_wavelength, n_real(1), n_imag(1),    & ! Inputs (Optical)
                 BTmat_bulk, BTmat_asymm, BTmat_ncoeffs, BTmat_expcoeffs, BTmat_Fmatrix,   & ! O
                 LPSD_BTmat_bulk(:,:,1), LPSD_BTmat_asymm(:,1),                & ! O
                 LPSD_BTmat_expcoeffs(:,:,:,1), LPSD_BTmat_Fmatrix(:,:,:,1),   & ! O
                 LRFE_BTmat_bulk(:,:,1), LRFE_BTmat_asymm(:,1),                & ! O
                 LRFE_BTmat_expcoeffs(:,:,:,1), LRFE_BTmat_Fmatrix(:,:,:,1),   & ! O
                 NLIN, BTmat_dist(:,1), LPSD_BTmat_dist(:,:,1),                & ! O
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2 )                     ! Outputs 
         else
            call tmat_master                         &
               ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere, Do_psd_OldStyle, & ! Inputs (Flags)
                 PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),         & ! Inputs (PSD)
                 Spheroid_type, nkmax(1), n_Fmatrix_angles, ndgs(1),                        & ! Inputs (General)
                 shape_factor(1), Tmat_accuracy, Tmat_wavelength, n_real(1), n_imag(1),     & ! Inputs (Optical)
                 Btmat_bulk, Btmat_asymm, Btmat_ncoeffs,                                    & ! Outputs (Tmat)
                 Btmat_expcoeffs, Btmat_Fmatrix, Btmat_dist(:,1),                           & ! Outputs (PSD)
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2 )                     ! Outputs 
         endif
      endif

!  Exception handling on everything

      if ( Fail2 ) then  
         Messages_Optical(1) = adjustl(trim(Tmat_message))
         Messages_Optical(2) = adjustl(trim(Tmat_trace))
         Messages_Optical(3) = adjustl(trim(Tmat_trace_2))
         Messages_Optical(4) = 'Single Tmatrix Call'; if ( do_bimodal ) Messages_Optical(4) = adjustl(trim(Tmat_trace_3))
         ctype = 'Regular   ' ; if ( do_aer_Jacobians ) ctype = 'Linearized'
         Messages_optical(5) = 'First call to the '//ctype//' Tmat program in AEROSOL CREATION, reference wavelength'
         return
      endif

!  Set the reference quantities. NOTE THE Bookkeeping.

      extinction_ref = BTmat_bulk(1) ; qm  = 0
      if ( do_aer_Jacobians ) then
         if ( Do_LinearRef) then
            L_extinction_ref(qm+1:qm+2) = LRFE_BTmat_bulk(1,1:2,1) ; qm = qm + 2
         endif
         if ( Do_LinearEps) then
            L_extinction_ref(qm+1:qm+1) = LRFE_BTmat_bulk(1,3,1) ; qm = qm + 1
         endif
         if ( Do_LinearPSD) then
            L_extinction_ref(qm+1:qm+npsd(1)) = LPSD_BTmat_bulk(1,1:npsd(1),1) ; qm = qm + npsd(1)
         endif
         if ( do_bimodal ) then
            if ( Do_LinearRef) then
               L_extinction_ref(qm+1:qm+2) = LRFE_BTmat_bulk(1,1:2,2) ; qm = qm + 2
            endif
            if ( Do_LinearEps) then
               L_extinction_ref(qm+1:qm+1) = LRFE_BTmat_bulk(1,3,2) ; qm = qm + 1
            endif
            if ( Do_LinearPSD) then
               L_extinction_ref(qm+1:qm+npsd(2)) = LPSD_BTmat_bulk(1,1:npsd(2),2) ; qm = qm + npsd(2)
            endif
            L_extinction_ref(qm+1) = LFRC_BTmat_bulk(1)
         endif
      endif

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

!  Local properties. Limited to the "Fixed option" here. 02/11/13

         Local_do_Varprops = .false.
         call get_variable_props &
           ( local_do_Varprops, do_bimodal, Max_Varpts, N_Varpts, wavelength,    &
             nreal_fixed, nimag_fixed, psdindex_fixed, psdpars_fixed,      &
             Aerwavs_Var, nreal_var, nimag_var, psdindex_var, psdpars_var, &
             PSD_Index, PSD_pars, n_real, n_imag )

!  progress

         if ( do_aer_jacobians ) then
            write(*,*)'Linearized Tmat Aerosol calculation, doing wavelength # ',wc
         else 
            write(*,*)'   Regular Tmat Aerosol calculation, doing wavelength # ',wc
         endif

!  wavelength for the Tmat code (micron unit)  

         Tmat_wavelength = wavelength/1000.0d0

!  First Call to develop aerosol Tmat properties
!  --------------------------------------------

!  Set the local Tmat program inputs (general)

         Do_Expcoeffs     = .TRUE.

!  BIMODAL CALL

         if ( do_bimodal) then
          if ( do_aer_Jacobians ) then
            call tmat_master_bimodal_PLUS                       &
               ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere,              & ! Inputs (Flags)
                 Do_LinearRef, Do_LinearEps, Do_LinearPSD,                              & ! Inputs (Flags)
                 Do_psd_OldStyle, psd_Index, psd_pars, MonoRadius, R1, R2, FixR1R2,     & ! Inputs (PSD)
                 bimodal_fraction, Spheroid_type, nkmax, n_Fmatrix_angles, ndgs,        & ! Inputs (General)
                 shape_factor, Tmat_accuracy, Tmat_wavelength, n_real, n_imag,          & ! Inputs (Optical)
                 BTmat_bulk, BTmat_asymm, BTmat_ncoeffs,     & ! O
                 BTmat_expcoeffs, BTmat_Fmatrix,             & ! O
                 LPSD_BTmat_bulk, LPSD_BTmat_asymm,          & ! O
                 LPSD_BTmat_expcoeffs, LPSD_BTmat_Fmatrix,   & ! O
                 LRFE_BTmat_bulk, LRFE_BTmat_asymm,          & ! O
                 LRFE_BTmat_expcoeffs, LRFE_BTmat_Fmatrix,   & ! O
                 LFRC_BTmat_bulk, LFRC_BTmat_asymm,          & ! O
                 LFRC_BTmat_expcoeffs, LFRC_BTmat_Fmatrix,   & ! O
                 NLIN, BTmat_dist, LPSD_BTmat_dist,          & ! O
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2, Tmat_trace_3 )   ! Outputs (status)
          else
            call tmat_master_bimodal                 &
               ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere,              & ! Inputs (Flags)
                 Do_psd_OldStyle, psd_Index, psd_pars, MonoRadius, R1, R2, FixR1R2,     & ! Inputs (PSD)
                 bimodal_fraction, Spheroid_type, nkmax, n_Fmatrix_angles, ndgs,        & ! Inputs (General)
                 shape_factor, Tmat_accuracy, Tmat_wavelength, n_real, n_imag,          & ! Inputs (Optical)
                 Btmat_bulk, Btmat_asymm, Btmat_ncoeffs,                                & ! Outputs (Tmat)
                 Btmat_expcoeffs, Btmat_Fmatrix, Btmat_dist,                            & ! Outputs (PSD)
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2, Tmat_trace_3 )   ! Outputs (status)
          endif
        endif

!  SINGLE CALL

        if ( .not. do_bimodal ) then
          Btmat_dist(:,2) = 0.0d0
          if ( do_aer_Jacobians ) then
            call tmat_master_PLUS                    &
              ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere,                  & ! Inputs (Flags) 
                Do_LinearRef, Do_LinearEps, Do_LinearPSD, Do_psd_OldStyle,                 & ! Inputs (Flags) 
                 PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),        & ! Inputs (PSD)
                 Spheroid_type, nkmax(1), n_Fmatrix_angles, ndgs(1),                       & ! Inputs (General)
                 shape_factor(1), Tmat_accuracy, Tmat_wavelength, n_real(1), n_imag(1),    & ! Inputs (Optical)
                 BTmat_bulk, BTmat_asymm, BTmat_ncoeffs, BTmat_expcoeffs, BTmat_Fmatrix,   & ! O
                 LPSD_BTmat_bulk(:,:,1), LPSD_BTmat_asymm(:,1),                & ! O
                 LPSD_BTmat_expcoeffs(:,:,:,1), LPSD_BTmat_Fmatrix(:,:,:,1),   & ! O
                 LRFE_BTmat_bulk(:,:,1), LRFE_BTmat_asymm(:,1),                & ! O
                 LRFE_BTmat_expcoeffs(:,:,:,1), LRFE_BTmat_Fmatrix(:,:,:,1),   & ! O
                 NLIN, BTmat_dist(:,1), LPSD_BTmat_dist(:,:,1),                & ! O
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2 )                     ! Outputs 
          else
            call tmat_master                         &
               ( Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse, Do_EqSaSphere, Do_psd_OldStyle, & ! Inputs (Flags)
                 PSD_Index(1), PSD_pars(:,1), MonoRadius, R1(1), R2(1), FixR1R2(1),         & ! Inputs (PSD)
                 Spheroid_type, nkmax(1), n_Fmatrix_angles, ndgs(1),                        & ! Inputs (General)
                 shape_factor(1), Tmat_accuracy, Tmat_wavelength, n_real(1), n_imag(1),     & ! Inputs (Optical)
                 Btmat_bulk, Btmat_asymm, Btmat_ncoeffs,                                    & ! Outputs (Tmat)
                 Btmat_expcoeffs, Btmat_Fmatrix, Btmat_dist(:,1),                           & ! Outputs (PSD)
                 fail2, istatus, Tmat_message, Tmat_trace, Tmat_trace_2 )                     ! Outputs 
          endif
        endif

!  Exception handling on everything

        if ( Fail2 ) then
          write(cwav,'(I3)')wc
          Messages_Optical(1) = adjustl(trim(Tmat_message))
          Messages_Optical(2) = adjustl(trim(Tmat_trace))
          Messages_Optical(3) = adjustl(trim(Tmat_trace_2))
          Messages_Optical(4) = 'Single Tmatrix Call'; if ( do_bimodal ) Messages_Optical(4) = adjustl(trim(Tmat_trace_3))
          ctype = 'Regular   ' ; if ( do_aer_Jacobians ) ctype = 'Linearized'
          Messages_Optical(5) = 'First call to the '//ctype//' Tmat program in AEROSOL CREATION, wavelength # '//cwav
          return
        endif

!  NOW set the optical property output
!  ===================================

!  >>>>>>>>>>>>>> New, 18 February 2013
!  set the distribution characteristics

         aerosol_distchars(:,:,w) = Btmat_dist

!  Set the reference quantities, if reference wavelength is in the list.
!   Values for the first (masked) wavelength.

         if ( lambda_index .ne. 0 .and. wc.eq.1 ) then
            extinction_ref = BTmat_bulk(1) ; qm  = 0
            if ( do_aer_Jacobians ) then
               if ( Do_LinearRef) then
                  L_extinction_ref(qm+1:qm+2) = LRFE_BTmat_bulk(1,1:2,1) ; qm = qm + 2
               endif
               if ( Do_LinearEps) then
                  L_extinction_ref(qm+1:qm+1) = LRFE_BTmat_bulk(1,3,1) ; qm = qm + 1
               endif
               if ( Do_LinearPSD) then
                  L_extinction_ref(qm+1:qm+npsd(1)) = LPSD_BTmat_bulk(1,1:npsd(1),1) ; qm = qm + npsd(1)
               endif
               if ( do_bimodal ) then
                  if ( Do_LinearRef) then
                     L_extinction_ref(qm+1:qm+2) = LRFE_BTmat_bulk(1,1:2,2) ; qm = qm + 2
                  endif
                  if ( Do_LinearEps) then
                     L_extinction_ref(qm+1:qm+1) = LRFE_BTmat_bulk(1,3,2) ; qm = qm + 1
                  endif
                  if ( Do_LinearPSD) then
                     L_extinction_ref(qm+1:qm+npsd(2)) = LPSD_BTmat_bulk(1,1:npsd(2),2) ; qm = qm + npsd(2)
                  endif
                  L_extinction_ref(qm+1) = LFRC_BTmat_bulk(1)
               endif
            endif
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

!  Convert Tmat derivatives to IOP derivatives
!  ==========================================

         if ( do_aer_Jacobians ) then

!  Division factor (necessary because T-matrix output is single-normalized, except for fraction Jacobian)
!      Have to unravel the single normalization here, in order to re-apply it later !!!!

            do qm = 1, qm1 - 1
               normfac(qm) = 1.0d0 / aerosol_pars(qm)
            enddo
            normfac(qm1) = 1.0d0

!  Linearized extinction

            qm = 0
            if ( Do_LinearRef) then
               L_extinction(qm+1:qm+2) = LRFE_BTmat_bulk(1,1:2,1) ; qm = qm + 2
            endif
            if ( Do_LinearEps) then
               L_extinction(qm+1:qm+1) = LRFE_BTmat_bulk(1,3,1) ; qm = qm + 1
            endif
            if ( Do_LinearPSD) then
               L_extinction(qm+1:qm+npsd(1)) = LPSD_BTmat_bulk(1,1:npsd(1),1) ; qm = qm + npsd(1)
            endif
            if ( do_bimodal ) then
               if ( Do_LinearRef) then
                  L_extinction(qm+1:qm+2) = LRFE_BTmat_bulk(1,1:2,2) ; qm = qm + 2
               endif
               if ( Do_LinearEps) then
                  L_extinction(qm+1:qm+1) = LRFE_BTmat_bulk(1,3,2) ; qm = qm + 1
               endif
               if ( Do_LinearPSD) then
                  L_extinction(qm+1:qm+npsd(2)) = LPSD_BTmat_bulk(1,1:npsd(2),2) ; qm = qm + npsd(2)
               endif
               L_extinction(qm+1) = LFRC_BTmat_bulk(1)
            endif

!  Linearization of the scaling factor

            do qm = 1, qm1
              L_aod_scaling(w,qm) = aod_scaling(w) * ( (L_extinction(qm)/extinction)        &
                                           - (L_extinction_ref(qm)/extinction_ref) )
            enddo
            qm2 = qm1 + 1
            L_aod_scaling(w,qm2) = aod_scaling(w)
            if ( loading_case .eq. 2 ) then
               L_aod_scaling(w,qmt) = aod_scaling(w)
            else if ( loading_case .eq. 3 ) then
               L_aod_scaling(w,qm2+1) = aod_scaling(w)
               L_aod_scaling(w,qmt)   = aod_scaling(w)
            endif

!  Assign Linearized SSAs and expansion coefficients

            qm = 0
            if ( Do_LinearRef) then
               do l = 0, n_scatmoms_w
                  do k = 1, nmuller
                     L_aerosol_scatmoms(qm+1:qm+2,k,l,w) = LRFE_BTmat_expcoeffs(k,l,1:2,1)
                  enddo
               enddo
               L_aerosol_ssalbs(qm+1:qm+2,w) = LRFE_BTmat_bulk(3,1:2,1) ; qm = qm + 2
            endif
            if ( Do_LinearEps) then
               do l = 0, n_scatmoms_w
                  do k = 1, nmuller
                     L_aerosol_scatmoms(qm+1:qm+1,k,l,w) = LRFE_BTmat_expcoeffs(k,l,3,1)
                  enddo
               enddo
               L_aerosol_ssalbs(qm+1:qm+1,w) = LRFE_BTmat_bulk(3,3,1) ; qm = qm + 1
            endif
            if ( Do_LinearPSD) then
               do l = 0, n_scatmoms_w
                  do k = 1, nmuller
                     L_aerosol_scatmoms(qm+1:qm+npsd(1),k,l,w) = LPSD_BTmat_expcoeffs(k,l,1:npsd(1),1)
                  enddo
               enddo
               L_aerosol_ssalbs(qm+1:qm+npsd(1),w) = LPSD_BTmat_bulk(3,1:npsd(1),1) ; qm = qm + npsd(1)
            endif

            if ( do_bimodal ) then
               if ( Do_LinearRef) then
                  do l = 0, n_scatmoms_w
                     do k = 1, nmuller
                        L_aerosol_scatmoms(qm+1:qm+2,k,l,w) = LRFE_BTmat_expcoeffs(k,l,1:2,2)
                     enddo
                  enddo
                  L_aerosol_ssalbs(qm+1:qm+2,w) = LRFE_BTmat_bulk(3,1:2,2) ; qm = qm + 2
               endif
               if ( Do_LinearEps) then
                  do l = 0, n_scatmoms_w
                     do k = 1, nmuller
                        L_aerosol_scatmoms(qm+1:qm+1,k,l,w) = LRFE_BTmat_expcoeffs(k,l,3,2)
                     enddo
                  enddo
                  L_aerosol_ssalbs(qm+1:qm+1,w) = LRFE_BTmat_bulk(3,3,2) ; qm = qm + 1
               endif
               if ( Do_LinearPSD) then
                  do l = 0, n_scatmoms_w
                     do k = 1, nmuller
                        L_aerosol_scatmoms(qm+1:qm+npsd(2),k,l,w) = LPSD_BTmat_expcoeffs(k,l,1:npsd(2),2)
                     enddo
                  enddo
                  L_aerosol_ssalbs(qm+1:qm+npsd(2),w) = LPSD_BTmat_bulk(3,1:npsd(2),2) ; qm = qm + npsd(2)
               endif
               do l = 0, n_scatmoms_w
                  do k = 1, nmuller
                     L_aerosol_scatmoms(qm1,k,l,w) = LFRC_BTmat_expcoeffs(k,l)
                  enddo
               enddo
               L_aerosol_ssalbs(qm1,w) = LFRC_BTmat_bulk(3)
            endif

!  Phase function first moment

            do qm = 1, qm1
               L_aerosol_scatmoms(qm,1,0,w) = 0.0d0
            enddo

!  End Jacobians clause

         endif

!  Update n_scatmoms

         n_scatmoms = min(n_scatmoms, n_scatmoms_w)

!  Un-normalize optical properties

         if ( do_aer_Jacobians ) then
            do qm = 1, qm1 - 1
               L_aod_scaling(w,qm) = normfac(qm) *  L_aod_scaling(w,qm)
               L_aerosol_ssalbs(qm,w) = L_aerosol_ssalbs(qm,w) * normfac(qm)
               do l = 0, n_scatmoms
                  do k = 1, nmuller
                     L_aerosol_scatmoms(qm,k,l,w) = L_aerosol_scatmoms(qm,k,l,w) * normfac(qm)
                  enddo
               enddo
            enddo
         endif

!  Assign Unscaled loadings + linearizations. VERSION 2

         aertau_unscaled = 0.0d0 ; L_aertau_unscaled = 0.0d0
         do n = 1, nlayers
            aerlayerflags(n) =  ( Loading(n) .ne. 0.0d0  )
            aertau_unscaled(n) = loading(n)
         enddo

         if ( do_aer_Jacobians ) then
            do n = 1, nlayers
               aertau_unscaled(n) = loading(n)
               do q = 1, qm1
                  L_aertau_unscaled(q,n) = loading(n)
               enddo
               L_aertau_unscaled(qm1+1,n) = Dloading_Dtau(n)
               if ( loading_case .eq. 2 ) then
                  L_aertau_unscaled(qm1+2,n) = Dloading_Dpars(n,1)
               else if ( loading_case .eq. 3 ) then
                  L_aertau_unscaled(qm1+2,n) = Dloading_Dpars(n,1)
                  L_aertau_unscaled(qm1+3,n) = Dloading_Dpars(n,2)
               endif
            enddo
         endif

!  Apply scalings to loadings and linearizations

         do n = 1, nlayers
            aerosol_deltau(n,w) = aertau_unscaled(n) * aod_scaling(w)
         enddo
         if ( do_aer_jacobians ) then
            do n = 1, nlayers
               do q = 1, n_aerosol_wfs
                  L_aerosol_deltau(q,n,w) = L_aertau_unscaled(q,n) * L_aod_scaling(q,w)
               enddo
            enddo
         endif

!  debug aerosol optical properties. VERSION TWO only

         if ( do_iopchecker ) then
          if ( do_aer_jacobians ) then
           do n = 1, nlayers
            if (aerlayerflags(n).and.n.eq.29 ) then
              write(777,'(i4,1p6e20.10)')n,aerosol_deltau(n,w),aerosol_ssalbs(w)
              write(777,'(i4,1p11e20.10)')n,(L_aerosol_deltau(q,n,w),q=1,11)
              write(777,'(i4,1p11e20.10)')n,(L_aerosol_ssalbs(q,w),q=1,11)
              do l = 0, n_scatmoms_w
                write(777,'(2i5,1pe20.10,1p11e15.6)')n,l,aerosol_scatmoms(1,l,w),(l_aerosol_scatmoms(q,1,l,w),q=1,11)
              enddo
            endif
           enddo
!           pause'Lin 777'
          else
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
!           pause'Reg 888'
          endif
         endif

!  End wavelength loop

      enddo

!  End Tmatrix aerosol clause

   endif


!  Finish

   return
end subroutine aerosol_jacobian_creation

!  End module

End Module aerosol_jacobian_creation_m

