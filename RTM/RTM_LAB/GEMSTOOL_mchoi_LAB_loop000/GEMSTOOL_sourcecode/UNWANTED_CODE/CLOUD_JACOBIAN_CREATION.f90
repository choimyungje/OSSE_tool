Module cloud_jacobian_creation_m

!  Type definitions

   use vlidort_pars, only : fpk

!  Use Mie code
!  ------------

!  parameters

   use RTSMie_parameters_m

!  Single call

   use RTSMie_sourcecode_m
   use RTSMie_sourcecode_plus_m

!  All routines are Public

private
public :: cloud_jacobian_creation

contains

subroutine cloud_jacobian_creation                                         &
          ( maxlayers, maxwav, maxmoments, maxcldwfs, do_cld_Jacobians,    & ! Dimensioning, flags and numbers
            cldtau_input_w1, reference_w1, cldlayerflags,                  & ! cldtau and flagging
            nlayers, nmuller, nwav, wav, height_grid,                      & ! INPUT, Numbers, misc.
            PSD_Index, PSD_pars, n_real, n_imag, momsize_cutoff,           & ! Mie Parameter input
            nblocks, nweights, FixR1R2, R1R2_cutoff, R1, R2, xlimit,       & ! Mie PSD control 
            Loading, Dloading_Dtau, n_cloud_wfs, n_Cld_wfs_1,  cloud_pars, & ! OUTPUT
            n_scatmoms, cloud_deltau,     cloud_ssalbs,   cloud_scatmoms,  & ! OUTPUT
            L_cloud_deltau, L_cloud_ssalbs, L_cloud_scatmoms,              & ! OUTPUT
            fail2, Messages_Mie )                                            ! Exception handling

!  =============================================================================
!                        CLOUD CREATION: METHOD 2 only
!  =============================================================================

!  Here is a breakdown of the Cloud weighting functions
!    Log-normal (2-parameter). Uniform loading

!   #               WEIGHTING FUNCTION

!   1         Real part of the refractive index for    Mode-1
!   2    Imaginary part of the refractive index for    Mode-1
!   3    Mode radius of Particle Size Distribution for Mode-1
!   4    Std. dev.   of Particle Size Distribution for Mode-1
!   5    Total optical depth of cloud at first wavelength

!  Here is a breakdown of the cloud weighting functions
!    MODIFIED GAMMA (3-parameter), Uniform loading

!   #               WEIGHTING FUNCTION

!   1         Real part of the refractive index for    Mode-1
!   2    Imaginary part of the refractive index for    Mode-1
!   3    Parameter 1 of Particle Size Distribution for Mode-1
!   4    Parameter 2 of Particle Size Distribution for Mode-1
!   5    Parameter 3 of Particle Size Distribution for Mode-1
!   6    Total optical depth of cloud at first wavelength

!  Generation of Cloud optical properties
!     1. Call to the Linearized Mie program
!     2. Convert Mie output (Microsopic) to IOP output (macroscopic) 

  implicit none

!  Input variables
!  ---------------

!  External Dimensioning

   integer, INTENT (IN) ::  maxlayers, maxwav, maxmoments, maxcldwfs

!  top level flags

   logical, INTENT (IN) :: do_cld_Jacobians

!  Numbers

   integer, INTENT (IN) ::  nlayers
   integer, INTENT (IN) ::  nwav

!  Heights and wavelengths

   REAL    (fpk),    INTENT (IN)   :: wav(maxwav)
   REAL    (fpk),    INTENT (IN)   :: height_grid(0:maxlayers)

!    nmuller = 1 (scalar code), = 6 (vector code)

   integer, INTENT (IN) ::  nmuller

!  Cloud loading control
!  TOTAL AOT = Reference cloud optical thickness (value at reference wavelength w1 )

   real(fpk),    INTENT (IN) :: cldtau_input_w1
   real(fpk),    INTENT (IN) :: reference_w1

!  Cloud loading flags

   logical,      dimension (maxlayers), INTENT (IN) :: cldlayerflags

!  Cutoff size for smallest phase function moment

   real(fpk),    INTENT (IN) :: momsize_cutoff, xlimit

!  Mie PSD inputs
!  --------------

!  Mie inputs (distribution index, PSD parameters)

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

   integer         , intent(in)  :: PSD_Index
   real    (fpk),    intent(in)  :: psd_pars (3)

!  FixR1R2 : If  set, Use "R1R2_cutoff" for smallest particle size
!                     then Internal routine to calculate R1 and R2 (outputs)
!            If Not set, Use Input R1 and R2 for PSD limits.

   logical         , intent(inout)  :: FixR1R2

!  R1, R2         - Minimum and Maximum radii (Microns)

   real    (fpk),    intent(inout)  :: R1, R2

!  R1R2_cutoff particle size for setting R1 and R2 internally

   REAL    (fpk),    intent(in)  :: R1R2_cutoff

!  PSD quadrature control
!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.

   INTEGER         , intent(in)  :: nblocks
   INTEGER         , intent(in)  :: nweights

!  N_REAL, N_IMAG - real and imaginary parts, refractive index (N-i.GE.0)

   real    (fpk),    intent(in)  :: n_real, n_imag

!  output
!  ------

!  Loading output

   real(fpk),    dimension (maxlayers), INTENT (OUT) :: loading
   real(fpk),    dimension (maxlayers), INTENT (OUT) :: Dloading_Dtau

!  General output
!    Total number of cloud weighting functions, cloud parameters. Up to 6 ( 5 + 1)
!    n_Cld_wfs = up to 5, from the linearized Mie program

   INTEGER         , INTENT (OUT)   :: n_cloud_wfs, n_Cld_wfs_1
   REAL    (fpk),    INTENT (OUT)   :: cloud_pars(maxcldwfs)

!  AOP output: Number of exapnsion coefficients to be used

   integer, INTENT (OUT)   ::  n_scatmoms

!  AOPS

    REAL(fpk),    DIMENSION( maxlayers, maxwav )      , INTENT (OUT)  :: CLOUD_DELTAU
    REAL(fpk),    DIMENSION( maxwav )                 , INTENT (OUT)  :: CLOUD_SSALBS 
    REAL(fpk),    DIMENSION( 6, 0:maxmoments, maxwav ), INTENT (OUT)  :: CLOUD_SCATMOMS

!  Linearized AOPs

    REAL(fpk),    DIMENSION( maxcldwfs, maxlayers, maxwav )       , INTENT (OUT) :: L_CLOUD_DELTAU
    REAL(fpk),    DIMENSION( maxcldwfs, maxwav )                  , INTENT (OUT) :: L_CLOUD_SSALBS 
    REAL(fpk),    DIMENSION( maxcldwfs, 6,  0:maxmoments, maxwav ), INTENT (OUT) :: L_CLOUD_SCATMOMS

!  Exception handling

   logical,        INTENT (OUT)           :: fail2
   character*(*),  INTENT (OUT)           :: Messages_Mie(4)

!  Mie code LOCAL INPUT variables (Not part of module input)
!  ==============================

!  Flag inputs
!  -----------

!      Do_Expcoeffs      - Boolean flag for computing Expansion Coefficients
!      Do_Fmatrix        - Boolean flag for computing F-matrix at equal-angles

   logical    :: Do_Expcoeffs
   logical    :: Do_Fmatrix

!      Do_Monodisperse   - Boolean flag for Doing a Monodisperse calculation
!                          If set, the PSD stuff will be turned off internally

   LOGICAL    :: do_Monodisperse

!  linearization w.r.t refractive index components

   logical    :: Do_LinearRef

!  Linearization for PSD parameters

   logical    :: Do_LinearPSD

!  F-matrix Angular control input
!  ------------------------------

!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]
!    NOT REQUIRED HERE

   INTEGER           :: n_Fmatrix_angles
   REAL    (fpk)     :: Fmatrix_angles(max_Mie_angles)

!  general variable applicable to both calls
!  -----------------------------------------

!  Limiting particle size value. Set to 10000.0 default.
!   If you exceed this, program will tell you to increase dimensioning.

   REAL    (fpk)      :: xparticle_limit

!  Monoradius     - Monodisperse radius size (Microns)

   real    (fpk)     :: Monoradius

!  Mie code OUTPUT variables
!  =========================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk)    :: BMie_bulk (3)

!  linearizations w.r.t. PSD parameters

   real(fpk)    :: LPSD_BMie_bulk (3,3)

!  linearizations w.r.t. Refractive index parameters

   real(fpk)    :: LRFE_BMie_bulk (3,2)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

!  Regular quantities

   integer      :: BMie_ncoeffs
   real(fpk)    :: BMie_expcoeffs (6,0:max_Mie_angles)
   real(fpk)    :: BMie_asymm

!  linearizations w.r.t. PSD parameters

   real(fpk)    :: LPSD_BMie_expcoeffs (6,0:max_Mie_angles,3)
   real(fpk)    :: LPSD_BMie_asymm(3)

!  linearizations w.r.t. Refractive index parameters

   real(fpk)    :: LRFE_BMie_expcoeffs (6,0:max_Mie_angles,2)
   real(fpk)    :: LRFE_BMie_asymm(2)

!  F-matrix,  optional output
!  --------------------------

!  F-matrix

   real(fpk)    :: BMie_Fmatrix(4,max_Mie_angles)

!  Linearizations of F-matrix

   real(fpk)    :: LPSD_BMie_Fmatrix(4,max_Mie_angles,3)
   real(fpk)    :: LRFE_BMie_Fmatrix(4,max_Mie_angles,2)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(fpk)     :: BMie_dist (5)
   real(fpk)     :: LPSD_BMie_dist (5,3)

!  Exception handling
!  ------------------

   integer       :: istatus
   character*90  :: Bmessages(3)


!  Local mask

   integer, dimension(9) :: distpars_mask
   data distpars_mask / 2,2,3,2,3,3,3,3,0 /

!  Local Variables
!  ===============

!  AOP output: Reference values
!    * These are the values at reference wavelength w1
!    * They are set for wavelength w1, then used again

    real(fpk)    :: extinction_ref, L_extinction_ref(maxcldwfs)

!  Help variables

    character*10 :: ctype
    character*3  :: cwav
    integer      :: n, k, l, m, kd, q, qm, qm1, qm2, qmt, lambda_index, w, wc, wavmask(maxwav),n_scatmoms_w,npsd
    real(fpk)    :: extent, density, scaling, wavelength, Mie_wavelength

!    logical, parameter :: do_iopchecker = .true.
    logical, parameter :: do_iopchecker = .false.

!  scaling

    real(fpk)    :: extinction,  L_extinction(maxcldwfs)
    real(fpk)    :: cod_scaling(Maxwav), L_cod_scaling(maxcldwfs,maxwav)
    real(fpk)    :: cldtau_unscaled(maxlayers), L_cldtau_unscaled(maxcldwfs,maxlayers)

!  Initialize output
!  =================

!  Initialize exception handling

   fail2 = .false.
   Messages_Mie    = ' '

!  initialize Loading

!    Loading = optical depth profile
!    Dloading_Dtau = derivative of profile w.r.t the total loading CLOUD_OPDEP

   loading  = 0.0d0
   Dloading_Dtau = 0.0d0

   extinction_ref = 0.0d0

   cloud_deltau = 0.0d0
   cloud_ssalbs = 0.0d0
   cloud_scatmoms = 0.0d0

   if ( do_cld_Jacobians ) then
     L_cloud_deltau = 0.0d0
     L_cloud_ssalbs = 0.0d0
     L_cloud_scatmoms = 0.0d0
   endif

!  Set local Mie variables

   do_monodisperse = .false.

   do_LinearRef = do_cld_Jacobians
   do_LinearPSD = do_cld_Jacobians

   Do_Fmatrix   = .false.             ! Always false
   DO_Expcoeffs = .false.             ! this will be set later on

   monoradius = 0.0d0
   xparticle_limit = xlimit

!  Do the loading
!  ==============

   extent = 0.0d0
   do n = 1, nlayers
      if ( cldlayerflags(n) ) then
         extent = extent + height_grid(n-1) - height_grid(n)
      endif
   enddo
   density = cldtau_input_w1 / extent
   do n = 1, nlayers
      if ( cldlayerflags(n) ) then
         scaling = (height_grid(n-1) - height_grid(n) ) / extent
         Loading(n) = density *  ( height_grid(n-1) - height_grid(n) )
         DLoading_Dtau(n) = scaling
      endif
   enddo

!  Bookkeeping
!  @@@@@@@@@@@

   n_Cld_wfs_1 = 0

!  Set up cloud parameters and total number of cloud WFs

   npsd = distpars_mask(PSD_index)
   m = 0
   if ( do_LinearRef ) then
      cloud_pars(1) = n_real ; cloud_pars(2) = n_imag
      m = m + 2
   endif
   if ( do_LinearPSD ) then
      do kd = 1, npsd
         m = m + 1 ; cloud_pars(m) = psd_pars(kd)
      end do
   endif
   n_Cld_wfs_1 = m ; qm1 = m

!  Add WFs for the loading
!    Total loading only for the Uniform profile.

   qmt = qm1 + 1

!  total number of weighing functions

   n_cloud_wfs = qmt

!  Add 1 parameter for loading 

   qm2 = qm1 + 1
   cloud_pars(qm2) = cldtau_input_w1

!  n_external_angles = 91
!  DO j = 1, n_external_angles
!     ang = d_two * dble(j-1)
!     external_angle_cosines(j) = COS(dacos(-1.0d0)*ang/180.0d0)
!  END DO

!  Check to see if reference wavelength is one of the set. Initialize Mask.

   lambda_index = 0
   do w = 1, nwav
     wavmask(w) = w
     if ( wav(w) .eq. reference_w1 ) lambda_index = w
   enddo

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

!  initialize

   n_scatmoms = 50000

!  Prepare the reference-wavelength Mie Inputs
!  ===========================================

   if ( lambda_index .eq. 0 ) then

!  Only require extinction coefficient and its linearization if flagged
!  Set the local Mie program inputs (bulk properties only)

      Do_Expcoeffs     = .FALSE.

!  reference wavelength

      wavelength           = reference_w1

!  progress

      if ( do_cld_jacobians ) then
         write(*,*)'Linearized Mie Cloud calculation, at reference wavelength'
      else 
         write(*,*)'Regular Mie Cloud calculation, at reference wavelength'
      endif

!  wavelength for the Mie code (micron unit)  

      Mie_wavelength = wavelength/1000.0d0

!  Call to the linearized or Regular Mie package

      if ( do_cld_Jacobians ) then
         call RTSMie_main_plus                                                     &  !---MIE CALL
         ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_LinearRef,                & ! I
           Do_LinearPSD, PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,         & ! I
           nblocks, nweights, xparticle_limit, R1R2_cutoff,  n_Fmatrix_angles,     & ! I
           Fmatrix_angles, Mie_wavelength, n_real, n_imag,                         & ! I
           BMie_bulk, BMie_asymm, BMie_ncoeffs,                            & ! O
           BMie_expcoeffs, BMie_Fmatrix, BMie_dist,                        & ! O
           LPSD_BMie_bulk, LPSD_BMie_asymm,  LPSD_BMie_expcoeffs,          & ! O
           LPSD_BMie_Fmatrix, LPSD_BMie_dist, LRFE_BMie_bulk,              & ! O
           LRFE_BMie_asymm, LRFE_BMie_expcoeffs, LRFE_BMie_Fmatrix,        & ! O
           fail2, istatus, Bmessages(1), Bmessages(2), Bmessages(3) )        ! O
      else
         call RTSMie_main                                                                  & !---MIE CALL
          ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, PSD_Index, PSD_pars,                & ! I
            MonoRadius, R1, R2, FixR1R2, nblocks, nweights, xparticle_limit, R1R2_cutoff,  & ! I
            n_Fmatrix_angles, Fmatrix_angles, Mie_wavelength, n_real, n_imag,              & ! I
            BMie_bulk, BMie_asymm, BMie_ncoeffs, BMie_expcoeffs, BMie_Fmatrix, BMie_dist,  & ! O
            fail2, istatus, Bmessages(1), Bmessages(2), Bmessages(3) )                       ! O
      endif

!  Exception handling

      if ( Fail2 ) then  
         do m = 1, 3   
            Messages_Mie(m) = adjustl(trim(Bmessages(m)))
         enddo
         ctype = 'Regular   ' ; if ( do_cld_Jacobians ) ctype = 'Linearized'
         Messages_Mie(4) = 'First call to the '//ctype//' Mie program in CLOUD CREATION, reference wavelength'
         return
      endif

!  Set the reference quantities

      extinction_ref = BMie_bulk(1) ; qm  = 0
      if ( do_cld_Jacobians ) then
         if ( Do_LinearRef) then
            L_extinction_ref(qm+1:qm+2) = LRFE_BMie_bulk(1,1:2) ; qm = qm + 2
         endif
         if ( Do_LinearPSD) then
            L_extinction_ref(qm+1:qm+npsd) = LPSD_BMie_bulk(1,1:npsd) ; qm = qm + npsd
         endif
      endif

!  End Mie calculation at reference wavelength

   endif

!  Prepare the multi-wavelength Mie Inputs
!  =======================================

!  wavelength loop. First wavelength will be the reference, if in  list.

   do wc = 1, nwav

!  Wavelengths [nm]

      w = wavmask(wc)
      wavelength = wav(w)

!  progress

      if ( do_cld_jacobians ) then
         write(*,*)'Linearized Mie Cloud calculation, doing wavelength # ',wc
      else 
         write(*,*)'Regular Mie Cloud calculation, doing wavelength # ',wc
      endif

!  wavelength for the Mie code (micron unit)  

      Mie_wavelength = wavelength/1000.0d0

!  Set the local Mie program inputs (general)

      Do_Expcoeffs     = .TRUE.

!  Call to the linearized or Regular Mie package

      if ( do_cld_Jacobians ) then
         call RTSMie_main_plus                                                     &  !---MIE CALL
         ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_LinearRef,                & ! I
           Do_LinearPSD, PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,         & ! I
           nblocks, nweights, xparticle_limit, R1R2_cutoff,  n_Fmatrix_angles,     & ! I
           Fmatrix_angles, Mie_wavelength, n_real, n_imag,                         & ! I
           BMie_bulk, BMie_asymm, BMie_ncoeffs,                            & ! O
           BMie_expcoeffs, BMie_Fmatrix, BMie_dist,                        & ! O
           LPSD_BMie_bulk, LPSD_BMie_asymm,  LPSD_BMie_expcoeffs,          & ! O
           LPSD_BMie_Fmatrix, LPSD_BMie_dist, LRFE_BMie_bulk,              & ! O
           LRFE_BMie_asymm, LRFE_BMie_expcoeffs, LRFE_BMie_Fmatrix,        & ! O
           fail2, istatus, Bmessages(1), Bmessages(2), Bmessages(3) )        ! O
      else
         call RTSMie_main                                                                  & !---MIE CALL
          ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, PSD_Index, PSD_pars,                & ! I
            MonoRadius, R1, R2, FixR1R2, nblocks, nweights, xparticle_limit, R1R2_cutoff,  & ! I
            n_Fmatrix_angles, Fmatrix_angles, Mie_wavelength, n_real, n_imag,              & ! I
            BMie_bulk, BMie_asymm, BMie_ncoeffs, BMie_expcoeffs, BMie_Fmatrix, BMie_dist,  & ! O
            fail2, istatus, Bmessages(1), Bmessages(2), Bmessages(3) )                       ! O
      endif

!  Exception handling

      if ( Fail2 ) then     
         write(cwav,'(I3)')wc
         do m = 1, 3   
            Messages_Mie(m) = adjustl(trim(Bmessages(m)))
         enddo
         ctype = 'Regular   ' ; if ( do_cld_Jacobians ) ctype = 'Linearized'
         Messages_Mie(4) = 'First call to the '//ctype//' Mie program in CLOUD CREATION, wavelength # '//cwav
         return
      endif

!  Part 2: CONVERT MIE output to COPs
!  ==================================

!  Set the reference quantities, if reference wavelength is in the list.
!   Values for the first (masked) wavelength

      if ( lambda_index .ne. 0 .and. wc.eq.1 ) then
         extinction_ref = BMie_bulk(1) ; qm  = 0
         if ( do_cld_Jacobians ) then
            if ( Do_LinearRef) then
               L_extinction_ref(qm+1:qm+2) = LRFE_BMie_bulk(1,1:2) ; qm = qm + 2
            endif
            if ( Do_LinearPSD) then
               L_extinction_ref(qm+1:qm+npsd) = LPSD_BMie_bulk(1,1:npsd) ; qm = qm + npsd
            endif
         endif
      endif

!  Extinction and its scaling factor

      extinction = BMie_bulk(1)
      cod_scaling(w) = extinction / extinction_ref

!  Assign SSAs and expansion coefficients, single cloud type

      cloud_ssalbs(w) = BMie_bulk(3)
      l = 0
      cloud_scatmoms(1,l,w) = 1.0d0
      do while (cloud_scatmoms(1,l,w).gt.momsize_cutoff.and.l.lt.maxmoments )
         l = l + 1
         cloud_scatmoms(1,l,w) = BMie_expcoeffs(1,l)
      enddo
      n_scatmoms_w = l - 1
      do l = 0, n_scatmoms_w
         do k = 2, nmuller
            cloud_scatmoms(k,l,w) = BMie_expcoeffs(k,l)
         enddo
      enddo

!  Update number of scattering moments

      n_scatmoms = min(n_scatmoms, n_scatmoms_w)

!  Convert Mie derivatives to IOP derivatives
!  ==========================================

      if ( do_cld_Jacobians ) then

!  Linearized extinction

         qm = 0
         if ( Do_LinearRef) then
            L_extinction(qm+1:qm+2) = LRFE_BMie_bulk(1,1:2) ; qm = qm + 2
         endif
         if ( Do_LinearPSD) then
            L_extinction(qm+1:qm+npsd) = LPSD_BMie_bulk(1,1:npsd) ; qm = qm + npsd
         endif

!  Linearization of the scaling factor

         do qm = 1, qm1
           L_cod_scaling(qm,w) = cod_scaling(w) * ( (L_extinction(qm)/extinction)        &
                                        - (L_extinction_ref(qm)/extinction_ref) )
         enddo
         qm2 = qm1 + 1
         L_cod_scaling(qm2,w) = cod_scaling(w)

!  Assign Linearized SSAs and expansion coefficients

         qm = 0
         if ( Do_LinearRef) then
            do l = 0, n_scatmoms_w
               do k = 1, nmuller
                  L_cloud_scatmoms(qm+1:qm+2,k,l,w) = LRFE_BMie_expcoeffs(k,l,1:2)
               enddo
            enddo
            L_cloud_ssalbs(qm+1:qm+2,w) = LRFE_BMie_bulk(3,1:2) ; qm = qm + 2
         endif
         if ( Do_LinearPSD) then
            do l = 0, n_scatmoms_w
               do k = 1, nmuller
                  L_cloud_scatmoms(qm+1:qm+npsd,k,l,w) = LPSD_BMie_expcoeffs(k,l,1:npsd)
               enddo
            enddo
            L_cloud_ssalbs(qm+1:qm+npsd,w) = LPSD_BMie_bulk(3,1:npsd) ; qm = qm + npsd
         endif

!  Phase function first moment

         do qm = 1, qm1
            L_cloud_scatmoms(qm,1,0,w) = d_zero
         enddo

!  End Jacobians clause

      endif

!  Assign Unscaled loadings + linearizations. VERSION 2

      cldtau_unscaled = 0.0d0 ; L_cldtau_unscaled = 0.0d0
      do n = 1, nlayers
         cldtau_unscaled(n) = loading(n)
      enddo
      if ( do_cld_Jacobians ) then
         qm1 = n_Cld_wfs_1
         do n = 1, nlayers
            cldtau_unscaled(n) = loading(n)
            do q = 1, qm1
               L_cldtau_unscaled(q,n) = loading(n)
            enddo
            L_cldtau_unscaled(qm1+1,n) = Dloading_Dtau(n)
         enddo
      endif
         
!  Apply scalings to loadings and linearizations

      do n = 1, nlayers
         cloud_deltau(n,w) = cldtau_unscaled(n) * cod_scaling(w)
      enddo
      if ( do_cld_jacobians ) then
         do n = 1, nlayers
            do q = 1, n_cloud_wfs
               L_cloud_deltau(q,n,w) = L_cldtau_unscaled(q,n) * L_cod_scaling(q,w)
            enddo
         enddo
      endif

!  debug cloud optical properties

      if ( do_iopchecker ) then
         if ( do_cld_jacobians ) then
            do n = 1, nlayers
               if (cldlayerflags(n).and.n.eq.29 ) then
                  write(67,'(i4,1p6e20.10)')n,cloud_deltau(n,w),cloud_ssalbs(w)
                  write(67,'(i4,1p11e20.10)')n,(L_cloud_deltau(q,n,w),q=1,11)
                  write(67,'(i4,1p11e20.10)')n,(L_cloud_ssalbs(q,w),q=1,11)
                  do l = 0, n_scatmoms_w
                     write(67,'(2i5,1pe20.10,1p11e15.6)')n,l,cloud_scatmoms(1,l,w),(l_cloud_scatmoms(q,1,l,w),q=1,11)
                  enddo
               endif
            enddo
!           pause'Lin 67'
         else
            do n = 1, nlayers
               if (cldlayerflags(N) ) then
                  write(68,'(i4,1p6e20.10)')n,cloud_deltau(n,w),cloud_ssalbs(w)
                  do l = 0, n_scatmoms_w
                     write(68,'(2i5,1p6e20.10)')n,l,(cloud_scatmoms(k,l,w),k=1,1)
                  enddo
               else
                  write(68,'(i4,1p6e20.10)')n,cloud_deltau(n,w)
               endif
            enddo
!           pause'Reg 68'
         endif
      endif

!  End wavelength loop

   enddo

!  Finish

   return
end subroutine cloud_jacobian_creation

!  End module

End module cloud_jacobian_creation_m
