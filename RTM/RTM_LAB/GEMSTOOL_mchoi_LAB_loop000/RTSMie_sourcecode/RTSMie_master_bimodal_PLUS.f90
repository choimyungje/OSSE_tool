module RTSMie_master_bimodal_plus_m

!  This is the Bimodal master #1 for RTS MIE code
!    ** RT Solutions, Version 1.4, 30 June     2011 (Bimodal control Tmatrix)
!    ** RT Solutions, Version 1.5, 17 August   2011 (Bimodal control Mie)

  use RTSMie_parameters_m
  use RTSMie_sourcecode_plus_m

!  Everything PUBLIC here
!  ----------------------

public

contains

subroutine RTSMie_master_bimodal_plus                        &
       ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
         Do_LinearRef, Do_LinearPSD,                         & ! I
         PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,   & ! I
         nblocks, nweights, xparticle_limit, R1R2_cutoff,    & ! I
         n_Fmatrix_angles, Fmatrix_angles,                   & ! I
         lambda, n_real, n_imag, fraction,                   & ! I
         BMie_bulk, BMie_asymm, BMie_ncoeffs,      & ! O
         BMie_expcoeffs, BMie_Fmatrix,             & ! O
         LPSD_BMie_bulk, LPSD_BMie_asymm,          & ! O
         LPSD_BMie_expcoeffs, LPSD_BMie_Fmatrix,   & ! O
         LRFE_BMie_bulk, LRFE_BMie_asymm,          & ! O
         LRFE_BMie_expcoeffs, LRFE_BMie_Fmatrix,   & ! O
         LFRC_BMie_bulk, LFRC_BMie_asymm,          & ! O
         LFRC_BMie_expcoeffs, LFRC_BMie_Fmatrix,   & ! O
         BMie_dist, LPSD_BMie_dist,                & ! O
         fail, istatus, Bmessages, trace_3 )         ! O

   USE RTSMie_parameters_m

!  implicit none statement

   IMPLICIT NONE

!  List of Inputs
!  ==============

!  Flag inputs
!  -----------

!      Do_Expcoeffs      - Boolean flag for computing Expansion Coefficients
!      Do_Fmatrix        - Boolean flag for computing F-matrix at equal-angles

   logical  , intent(in)  :: Do_Expcoeffs
   logical  , intent(in)  :: Do_Fmatrix

!      Do_Monodisperse   - Boolean flag for Doing a Monodisperse calculation
!                          If set, the PSD stuff will be turned off internally

   LOGICAL  , INTENT (IN) :: do_Monodisperse

!  Linearization control
!      Do_LinearRef      - Boolean Flag for doing Refractive Index linearization
!      Do_LinearPSD      - Boolean Flag for doing PSD linearization
!                          This is checked and turned off for Monodisperse

   logical  , intent(in)     :: Do_LinearRef
   logical  , intent(inout)  :: Do_LinearPSD

!  PSD inputs
!  ----------

!  FixR1R2 : If  set, Use "R1R2_cutoff" for smallest particle size
!                     then Internal routine to calculate R1 and R2 (outputs)
!            If Not set, Use Input R1 and R2 for PSD limits.

   logical, intent(inout)  :: FixR1R2(2)

!  R1, R2         - Minimum and Maximum radii (Microns)

   real    (KIND=dp), intent(inout)  :: R1(2), R2(2)

!  Limiting particle size value. Set to 10000.0 default.
!   If you exceed this, program will tell you to increase dimensioning.

   REAL    (KIND=dp), INTENT (IN) :: xparticle_limit

!      Monoradius     - Monodisperse radius size (Microns)

   real    (KIND=dp), intent(in)  :: Monoradius

!      psd_Index      - Index for particle size distribution of spheres
!      psd_pars       - Parameters characterizing PSD (up to 3 allowed)

!  Mie inputs (distribution index, PSD parameters)
!    PSD_index = 1 : TWO-PARAMETER GAMMA with alpha and b given
!    PSD_index = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given
!    PSD_index = 3 : BIMODAL GAMMA with equal mode weights
!    PSD_index = 4 : LOG-NORMAL with rg and sigma given
!    PSD_index = 5 : LOG-NORMAL with reff and veff given
!    PSD_index = 6 : POWER LAW
!    PSD_index = 7 : MODIFIED GAMMA with alpha, rc and gamma given
!    PSD_index = 8 : MODIFIED GAMMA with alpha, b and gamma given

   integer          , intent(in)  :: PSD_Index(2)
   real    (KIND=dp), intent(in)  :: PSD_pars (3,2)

!  PSD quadrature control
!  ----------------------

!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.
!    R1R2_cutoff particle size for setting R1 and R2 internally

   INTEGER          , INTENT (IN) :: nblocks(2)
   INTEGER          , INTENT (IN) :: nweights(2)
   REAL    (KIND=dp), INTENT (IN) :: R1R2_cutoff(2)

!  Optical: Wavelength, refractive index
!  -------------------------------------

!      LAMBDA         - wavelength of light (microns)
!      N_REAL, N_IMAG - real and imaginary parts, refractive index (N-i.GE.0)

   real    (KIND=dp), intent(in)  :: lambda, n_real(2), n_imag(2)

!  Fraction (by volume weight, First PSD mode)

   real    (KIND=dp), intent(in)  :: fraction

!  F-matrix Angular control input
!  ------------------------------

!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER          , INTENT (IN) :: n_Fmatrix_angles
   REAL    (KIND=dp), INTENT (IN) :: Fmatrix_angles(max_Mie_angles)

!  Output arguments
!  ================

!  Bulk distribution parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(KIND=dp), intent(out) :: BMie_bulk (3)

!  linearizations w.r.t. PSD parameters

   real(KIND=dp), intent(out) :: LPSD_BMie_bulk (3,3,2)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(KIND=dp), intent(out) ::  LRFE_BMie_bulk (3,2,2)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

!  Regular quantities

   integer      , intent(out) :: BMie_ncoeffs
   real(KIND=dp), intent(out) :: BMie_expcoeffs (6,0:max_Mie_angles)
   real(KIND=dp), intent(out) :: BMie_asymm

!  linearizations w.r.t. PSD parameters

   real(KIND=dp), intent(out) :: LPSD_BMie_expcoeffs (6,0:max_Mie_angles,3,2)
   real(KIND=dp), intent(out) :: LPSD_BMie_asymm(3,2)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(KIND=dp), intent(out) :: LRFE_BMie_expcoeffs (6,0:max_Mie_angles,2,2)
   real(KIND=dp), intent(out) :: LRFE_BMie_asymm(2,2)

!  F-matrix,  optional output
!  --------------------------

!  F-matrix

   real(KIND=dp), intent(out) :: BMie_Fmatrix(4,max_Mie_angles)

!  Linearizations of F-matrix

   real(KIND=dp), intent(out) :: LPSD_BMie_Fmatrix(4,max_Mie_angles,3,2)
   real(KIND=dp), intent(out) :: LRFE_BMie_Fmatrix(4,max_Mie_angles,2,2)

!  Fraction Jacobian
!  -----------------

   real(kind=dp), intent(out)  :: LFRC_BMie_bulk (3)
   real(kind=dp), intent(out)  :: LFRC_BMie_expcoeffs (6,0:max_Mie_angles)
   real(kind=dp), intent(out)  :: LFRC_BMie_asymm
   real(kind=dp), intent(out)  :: LFRC_BMie_Fmatrix (4,max_Mie_angles)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(KIND=dp), intent(out) :: BMie_dist (5,2)
   real(KIND=dp), intent(out) :: LPSD_Bmie_dist (5,3,2)

!  Exception handling
!  ------------------

   LOGICAL          , INTENT (OUT)   :: fail
   INTEGER          , INTENT (OUT)   :: istatus
   CHARACTER*(*)    , INTENT (OUT)   :: Bmessages(3), trace_3

!  Local Arrays
!  ============

!  Bulk distribution parameters
!  Expansion coefficients and Asymmetry parameter
!  F-matrix,  optional output

   real(KIND=dp) :: Mie1_bulk (3)
   integer       :: Mie1_ncoeffs
   real(KIND=dp) :: Mie1_expcoeffs (6,0:max_Mie_angles)
   real(KIND=dp) :: Mie1_asymm
   real(KIND=dp) :: Mie1_Fmatrix (4,max_Mie_angles)

   real(KIND=dp) :: Mie2_bulk (3)
   integer       :: Mie2_ncoeffs
   real(KIND=dp) :: Mie2_expcoeffs (6,0:max_Mie_angles)
   real(KIND=dp) :: Mie2_asymm
   real(KIND=dp) :: Mie2_Fmatrix (4,max_Mie_angles)

   real(KIND=dp) :: LRFE_Mie1_bulk (3,2)
   real(KIND=dp) :: LRFE_Mie1_expcoeffs (6,0:max_Mie_angles,2)
   real(KIND=dp) :: LRFE_Mie1_asymm(2)
   real(KIND=dp) :: LRFE_Mie1_Fmatrix (4,max_Mie_angles,2)

   real(KIND=dp) :: LRFE_Mie2_bulk (3,2)
   real(KIND=dp) :: LRFE_Mie2_expcoeffs (6,0:max_Mie_angles,2)
   real(KIND=dp) :: LRFE_Mie2_asymm(2)
   real(KIND=dp) :: LRFE_Mie2_Fmatrix (4,max_Mie_angles,2)

   real(KIND=dp) :: LPSD_Mie1_bulk (3,3)
   real(KIND=dp) :: LPSD_Mie1_expcoeffs (6,0:max_Mie_angles,3)
   real(KIND=dp) :: LPSD_Mie1_asymm(3)
   real(KIND=dp) :: LPSD_Mie1_Fmatrix (4,max_Mie_angles,3)

   real(KIND=dp) :: LPSD_Mie2_bulk (3,3)
   real(KIND=dp) :: LPSD_Mie2_expcoeffs (6,0:max_Mie_angles,3)
   real(KIND=dp) :: LPSD_Mie2_asymm(3)
   real(KIND=dp) :: LPSD_Mie2_Fmatrix (4,max_Mie_angles,3)

!  Other local variables
!  ---------------------

   integer       :: k, L, q, nlin
   real(KIND=dp) :: FF1, FF2, WW1, WW2, TERM1, TERM2
   real(KIND=dp) :: Csca1_FF1, Csca2_FF2, Csca_total, D1_Csca1(3), D2_Csca2(3)
   real(KIND=dp) :: D1_WW1(3), D2_WW1(3), D1_WW2(3), D2_WW2(3), DF_WW1

!  Zero the output
!  ---------------

   BMie_bulk      = d_zero
   BMie_Fmatrix   = d_zero
   BMie_expcoeffs = d_zero
   BMie_asymm     = d_zero
   BMie_ncoeffs   = 0

   LPSD_BMie_bulk      = d_zero
   LPSD_BMie_Fmatrix   = d_zero
   LPSD_BMie_expcoeffs = d_zero
   LPSD_BMie_asymm     = d_zero

   LRFE_BMie_bulk      = d_zero
   LRFE_BMie_Fmatrix   = d_zero
   LRFE_BMie_expcoeffs = d_zero
   LRFE_BMie_asymm     = d_zero

   LFRC_BMie_bulk      = d_zero
   LFRC_BMie_Fmatrix   = d_zero
   LFRC_BMie_expcoeffs = d_zero
   LFRC_BMie_asymm     = d_zero

   BMie_dist      = d_zero
   LPSD_BMie_dist = d_zero

   FF1 = fraction
   FF2 = d_one - FF1

   fail           = .false.
   istatus        = 0
   Bmessages(1:3) = ' '
   trace_3        = ' '

!  Check: No Monodisperse here !
!  -----------------------------

   if ( do_monodisperse ) then
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal_PLUS module: Input error: MONODISPERSE FLAG must be Off!'
      return
   endif

!  First Call
!  ---------

   k = 1 ; write(*,*)' ** Doing Mie Linearized for PSD # 1 ----------------------'

   CALL RTSMie_main_plus                                       & 
       ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,            & ! I
         Do_LinearRef, Do_LinearPSD,                           & ! I
         PSD_Index(k), PSD_pars(:,k), MonoRadius,              & ! I
         R1(k), R2(k), FixR1R2(k), nblocks(k), nweights(k),    & ! I
         xparticle_limit, R1R2_cutoff(k),                      & ! I
         n_Fmatrix_angles, Fmatrix_angles,                     & ! I
         lambda, n_real(k), n_imag(k),                         & ! I
         Mie1_bulk, Mie1_asymm, Mie1_ncoeffs,                    & ! O
         Mie1_expcoeffs, Mie1_Fmatrix, BMie_dist(:,k),           & ! O
         LPSD_Mie1_bulk, LPSD_Mie1_asymm,                        & ! O
         LPSD_Mie1_expcoeffs, LPSD_Mie1_Fmatrix, LPSD_BMie_dist(:,:,k),  & ! O
         LRFE_Mie1_bulk, LRFE_Mie1_asymm,                        & ! O
         LRFE_Mie1_expcoeffs, LRFE_Mie1_Fmatrix,                 & ! O
         fail, istatus, Bmessages(1), Bmessages(2), Bmessages(3) ) ! O

!  Exception handling

   if ( fail ) then
      trace_3 = 'RTSMie_master_bimodal_PLUS module: First PSD call, Error'
      if ( Istatus .eq. 1 ) return
   endif

!  Second call
!  -----------

   k = 2 ; write(*,*)' ** Doing Mie Linearized for PSD # 2 ----------------------'

   CALL RTSMie_main_plus                                       & 
       ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,            & ! I
         Do_LinearRef, Do_LinearPSD,                           & ! I
         PSD_Index(k), PSD_pars(:,k), MonoRadius,              & ! I
         R1(k), R2(k), FixR1R2(k), nblocks(k), nweights(k),    & ! I
         xparticle_limit, R1R2_cutoff(k),                      & ! I
         n_Fmatrix_angles, Fmatrix_angles,                     & ! I
         lambda, n_real(k), n_imag(k),                         & ! I
         Mie2_bulk, Mie2_asymm, Mie2_ncoeffs,                    & ! O
         Mie2_expcoeffs, Mie2_Fmatrix, BMie_dist(:,k),           & ! O
         LPSD_Mie2_bulk, LPSD_Mie2_asymm,                        & ! O
         LPSD_Mie2_expcoeffs, LPSD_Mie2_Fmatrix, LPSD_BMie_dist(:,:,k),  & ! O
         LRFE_Mie2_bulk, LRFE_Mie2_asymm,                        & ! O
         LRFE_Mie2_expcoeffs, LRFE_Mie2_Fmatrix,                 & ! O
         fail, istatus, Bmessages(1), Bmessages(2), Bmessages(3) ) ! O

!  Exception handling

   if ( fail ) then
      trace_3 = 'RTSMie_master_bimodal_PLUS module: second PSD call, Error'
      if ( Istatus .eq. 1 ) return
   endif

!  Bimodal determination
!  =====================

!  Revision 20 September 2011
!    Correct definition for Expcoeffs/Fmatrix: WW1/WW2 in place of FF1/FF2

   Csca1_FF1  = FF1 * Mie1_bulk(2)   
   Csca2_FF2  = FF2 * Mie2_bulk(2)   
   Csca_total =  Csca1_FF1 + Csca2_FF2
   WW1   = Csca1_FF1 / Csca_total
   WW2   = Csca2_FF2 / Csca_total

!  Bimodal quantities

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong
   BMie_bulk(1:2) = FF1 * Mie1_bulk(1:2) + FF2 * Mie2_bulk(1:2)
   BMie_bulk(3)   = BMie_bulk(2) / BMie_bulk(1)
!  original code
!   BMie_bulk(1:3) = FF1 * Mie1_bulk(1:3) + FF2 * Mie2_bulk(1:3)

   if ( Do_Expcoeffs ) then
      BMie_asymm = WW1 * Mie1_asymm + WW2 * Mie2_asymm
      BMie_ncoeffs = max(Mie1_ncoeffs,Mie2_ncoeffs)      
      do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
         BMie_expcoeffs(1:6,L) = WW1 * Mie1_expcoeffs(1:6,L) + & 
                                 WW2 * Mie2_expcoeffs(1:6,L)
      enddo
      if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
          do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
              BMie_expcoeffs(1:6,L) = WW2 * Mie2_expcoeffs(1:6,L)
          enddo
      else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
          do L = Mie2_ncoeffs + 1,Mie1_ncoeffs
              BMie_expcoeffs(1:6,L) = WW1 * Mie1_expcoeffs(1:6,L)
          enddo
      endif
   endif
   if ( Do_Fmatrix ) then
      do L = 1, n_Fmatrix_angles
         BMie_Fmatrix(1:4,L) = WW1 * Mie1_Fmatrix(1:4,L) + & 
                               WW2 * Mie2_Fmatrix(1:4,L)
      enddo
   endif 
  
!  LRFE Linearizations
!  ===================

   if ( do_LinearRef ) then
      nlin = 2

!  help variables
!  --------------

      D1_Csca1(1:nlin) = LRFE_Mie1_bulk(2,1:nlin)
      D2_Csca2(1:nlin) = LRFE_Mie2_bulk(2,1:nlin)
      D1_WW1(1:nlin) = FF1 * D1_Csca1(1:nlin) * WW2 / Csca_total
      D2_WW2(1:nlin) = FF2 * D2_Csca2(1:nlin) * WW1 / Csca_total
      D1_WW2(1:nlin) = - D1_WW1(1:nlin)
      D2_WW1(1:nlin) = - D2_WW2(1:nlin)

!  Bulk variables
!  --------------

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong
      LRFE_BMie_bulk(1:2,1:nlin,1) = FF1 * LRFE_Mie1_bulk(1:2,1:nlin)
      LRFE_BMie_bulk(1:2,1:nlin,2) = FF2 * LRFE_Mie2_bulk(1:2,1:nlin)
      LRFE_BMie_bulk(3,1:nlin,1) = ( LRFE_BMie_bulk(2,1:nlin,1) &
                 - BMie_bulk(3) * LRFE_BMie_bulk(1,1:nlin,1) ) / BMie_bulk(1)
      LRFE_BMie_bulk(3,1:nlin,2) = ( LRFE_BMie_bulk(2,1:nlin,2)  &
                 - BMie_bulk(3) * LRFE_BMie_bulk(1,1:nlin,2) ) / BMie_bulk(1)
!  original code
!      LRFE_BMie_bulk(1:3,1:nlin,1) = FF1 * LRFE_Mie1_bulk(1:3,1:nlin)
!      LRFE_BMie_bulk(1:3,1:nlin,2) = FF2 * LRFE_Mie2_bulk(1:3,1:nlin)

!      write(16,68) 'Bulk-B-original',(BMie_bulk(k),k=1,3)
!      write(16,68) 'dBulk_Dist1_dN1',(n_real(1)*LRFE_BMie_bulk(k,1,1),k=1,3)
!      write(16,68) 'dBulk_Dist1_dN2',(n_imag(1)*LRFE_BMie_bulk(k,2,1),k=1,3)
!      write(16,68) 'dBulk_Dist2_dN1',(n_real(2)*LRFE_BMie_bulk(k,1,2),k=1,3)
!      write(16,68) 'dBulk_Dist2_dN2',(n_imag(2)*LRFE_BMie_bulk(k,2,2),k=1,3)

!  Asymmetry parameter
!  -------------------

      if ( Do_Expcoeffs ) then
         do q = 1, nlin
            TERM1 = WW1 * LRFE_Mie1_asymm(q)
            TERM2 = D1_WW1(q) * Mie1_asymm + D1_WW2(q) * Mie2_asymm
            LRFE_BMie_asymm(q,1) = TERM1 + TERM2
            TERM1 = WW2 * LRFE_Mie2_asymm(q)
            TERM2 = D2_WW1(q) * Mie1_asymm + D2_WW2(q) * Mie2_asymm
            LRFE_BMie_asymm(q,2) = TERM1 + TERM2
         enddo
      endif

!  Expansion coefficients
!  ----------------------

      if ( Do_Expcoeffs ) then

!  w.r.t  microphysical parameters for Particles in Mode 1

         do q = 1, nlin
            do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW1 * LRFE_Mie1_expcoeffs(k,L,q)
                  TERM2 = D1_WW1(q) * Mie1_expcoeffs(k,L) + D1_WW2(q) * Mie2_expcoeffs(k,L)
                  LRFE_BMie_expcoeffs(k,L,q,1) = TERM1 + TERM2
               enddo
            enddo
            if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
               do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
                  do k = 1, 6
                     LRFE_BMie_expcoeffs(k,L,q,1) = D1_WW2(q) * Mie2_expcoeffs(k,L)
                  enddo
               enddo  
            else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
               do L = Mie2_ncoeffs + 1,Mie1_ncoeffs
                  do k = 1, 6
                     TERM1 = WW1 * LRFE_Mie1_expcoeffs(k,L,q)
                     LRFE_BMie_expcoeffs(k,L,q,1) = TERM1 + D1_WW1(q) * Mie1_expcoeffs(k,L)
                  enddo
               enddo  
            endif
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW2 * LRFE_Mie2_expcoeffs(k,L,q)
                  TERM2 = D2_WW1(q) * Mie1_expcoeffs(k,L) + D2_WW2(q) * Mie2_expcoeffs(k,L)
                  LRFE_BMie_expcoeffs(k,L,q,2) = TERM1 + TERM2
               enddo
            enddo
            if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
               do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
                  do k = 1, 6
                     TERM1 = WW2 * LRFE_Mie2_expcoeffs(k,L,q)
                     LRFE_BMie_expcoeffs(k,L,q,2) = TERM1 + D2_WW2(q) * Mie2_expcoeffs(k,L)
                  enddo
               enddo  
            else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
               do L = Mie2_ncoeffs + 1,Mie1_ncoeffs
                  do k = 1, 6
                     LRFE_BMie_expcoeffs(k,L,q,2) = D2_WW1(q) * Mie1_expcoeffs(k,L)
                  enddo
               enddo  
            endif
         enddo

      endif

!  Fmatrix
!  -------

      if ( Do_Fmatrix ) then

!  w.r.t  microphysical parameters for Particles in Mode 1

         do q = 1, nlin
            do L = 1, n_Fmatrix_angles
               do k = 1, 4
                  TERM1 = WW1 * LRFE_Mie1_Fmatrix(k,L,q)
                  TERM2 = D1_WW1(q) * Mie1_Fmatrix(k,L) + D1_WW2(q) * Mie2_Fmatrix(k,L)
                  LRFE_BMie_Fmatrix(k,L,q,1) = TERM1 + TERM2
               enddo
            enddo
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 1, n_Fmatrix_angles
               do k = 1, 4
                  TERM1 = WW2 * LRFE_Mie2_Fmatrix(k,L,q)
                  TERM2 = D2_WW1(q) * Mie1_Fmatrix(k,L) + D2_WW2(q) * Mie2_Fmatrix(k,L)
                  LRFE_BMie_Fmatrix(k,L,q,2) = TERM1 + TERM2
               enddo
            enddo
         enddo

!  End clauses

      endif   
   endif

!  LPSD Linearizations
!  ===================

   if ( do_LinearPSD ) then

!  help variables
!  --------------

!  Bug discovered here by R. Spurr, June 2014

      nlin = 3
      D1_Csca1(1:nlin) = LPSD_Mie1_bulk(2,1:nlin)
      D2_Csca2(1:nlin) = LPSD_Mie2_bulk(2,1:nlin)
      D1_WW1(1:nlin) = FF1 * D1_Csca1(1:nlin) * WW2 / Csca_total
      D2_WW2(1:nlin) = FF2 * D2_Csca2(1:nlin) * WW1 / Csca_total
      D1_WW2(1:nlin) = - D1_WW1(1:nlin)
      D2_WW1(1:nlin) = - D2_WW2(1:nlin)       ! This line was improperly ordered

!  Bulk variables
!  --------------

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong

      do q = 1, nlin
         LPSD_BMie_bulk(1:2,q,1) = FF1 * LPSD_Mie1_bulk(1:2,q)
         LPSD_BMie_bulk(1:2,q,2) = FF2 * LPSD_Mie2_bulk(1:2,q)
         LPSD_BMie_bulk(3,q,1) = ( LPSD_BMie_bulk(2,q,1) &
                 - BMie_bulk(3) * LPSD_BMie_bulk(1,q,1) ) / BMie_bulk(1)
         LPSD_BMie_bulk(3,q,2) = ( LPSD_BMie_bulk(2,q,2) &
                 - BMie_bulk(3) * LPSD_BMie_bulk(1,q,2) ) / BMie_bulk(1)
      enddo

!      write(17,68) 'Bulk-B-original',(BMie_bulk(k),k=1,3)
!      write(17,68) 'dBulk_Dist1_dP1',(PSD_pars(1,1)*LPSD_BMie_bulk(k,1,1),k=1,3)
!      write(17,68) 'dBulk_Dist1_dP2',(PSD_pars(2,1)*LPSD_BMie_bulk(k,2,1),k=1,3)
!      write(17,68) 'dBulk_Dist2_dP1',(PSD_pars(1,2)*LPSD_BMie_bulk(k,1,2),k=1,3)
!      write(17,68) 'dBulk_Dist2_dP2',(PSD_pars(2,2)*LPSD_BMie_bulk(k,2,2),k=1,3)
!68 format(a15,1p3e20.10)

!  original code
!      LPSD_BMie_bulk(1:3,1:nlin,1) = FF1 * LPSD_Mie1_bulk(1:3,1:nlin)
!      LPSD_BMie_bulk(1:3,1:nlin,2) = FF2 * LPSD_Mie2_bulk(1:3,1:nlin)

!  Asymmetry parameter
!  -------------------

      if ( Do_Expcoeffs ) then
         do q = 1, nlin
            TERM1 = WW1 * LPSD_Mie1_asymm(q)
            TERM2 = D1_WW1(q) * Mie1_asymm + D1_WW2(q) * Mie2_asymm
            LPSD_BMie_asymm(q,1) = TERM1 + TERM2
            TERM1 = WW2 * LPSD_Mie2_asymm(q)
            TERM2 = D2_WW1(q) * Mie1_asymm + D2_WW2(q) * Mie2_asymm
            LPSD_BMie_asymm(q,2) = TERM1 + TERM2
         enddo
      endif

!  Expansion coefficients
!  ----------------------

      if ( Do_Expcoeffs ) then

!  w.r.t  microphysical parameters for Particles in Mode 1

         do q = 1, nlin
            do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW1 * LPSD_Mie1_expcoeffs(k,L,q)
                  TERM2 = D1_WW1(q) * Mie1_expcoeffs(k,L) + D1_WW2(q) * Mie2_expcoeffs(k,L)
                  LPSD_BMie_expcoeffs(k,L,q,1) = TERM1 + TERM2
               enddo
            enddo
            if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
               do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
                  do k = 1, 6
                     LPSD_BMie_expcoeffs(k,L,q,1) = D1_WW2(q) * Mie2_expcoeffs(k,L)
                  enddo
               enddo  
            else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
               do L = Mie2_ncoeffs + 1,Mie1_ncoeffs
                  do k = 1, 6
                     TERM1 = WW1 * LPSD_Mie1_expcoeffs(k,L,q)
                     LPSD_BMie_expcoeffs(k,L,q,1) = TERM1 + D1_WW1(q) * Mie1_expcoeffs(k,L)
                  enddo
               enddo  
            endif
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW2 * LPSD_Mie2_expcoeffs(k,L,q)
                  TERM2 = D2_WW1(q) * Mie1_expcoeffs(k,L) + D2_WW2(q) * Mie2_expcoeffs(k,L)
                  LPSD_BMie_expcoeffs(k,L,q,2) = TERM1 + TERM2
               enddo
            enddo
            if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
               do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
                  do k = 1, 6
                     TERM1 = WW2 * LPSD_Mie2_expcoeffs(k,L,q)
                     LPSD_BMie_expcoeffs(k,L,q,2) = TERM1 + D2_WW2(q) * Mie2_expcoeffs(k,L)
                  enddo
               enddo  
            else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
               do L = Mie2_ncoeffs + 1,Mie1_ncoeffs
                  do k = 1, 6
                     LPSD_BMie_expcoeffs(k,L,q,2) = D2_WW1(q) * Mie1_expcoeffs(k,L)
                  enddo
               enddo  
            endif
         enddo

      endif

!  Fmatrix
!  -------

      if ( Do_Fmatrix ) then

!  w.r.t  microphysical parameters for Particles in Mode 1

         do q = 1, nlin
            do L = 1, n_Fmatrix_angles
               do k = 1, 4
                  TERM1 = WW1 * LPSD_Mie1_Fmatrix(k,L,q)
                  TERM2 = D1_WW1(q) * Mie1_Fmatrix(k,L) + D1_WW2(q) * Mie2_Fmatrix(k,L)
                  LPSD_BMie_Fmatrix(k,L,q,1) = TERM1 + TERM2
               enddo
            enddo
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 1, n_Fmatrix_angles
               do k = 1, 4
                  TERM1 = WW2 * LPSD_Mie2_Fmatrix(k,L,q)
                  TERM2 = D2_WW1(q) * Mie1_Fmatrix(k,L) + D2_WW2(q) * Mie2_Fmatrix(k,L)
                  LPSD_BMie_Fmatrix(k,L,q,2) = TERM1 + TERM2
               enddo
            enddo
         enddo

!  End clauses

      endif   
   endif

!  Fractional linearization. NOT NORMALIZED
!  ========================================

!  Help variables

   DF_WW1 = ( WW1 * Mie2_bulk(2) + WW2 * Mie1_bulk(2) ) / Csca_total 

!  Bulk

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong
   LFRC_BMie_bulk(1:2) = Mie1_bulk(1:2) - Mie2_bulk(1:2)
   LFRC_BMie_bulk(3) = ( LFRC_BMie_bulk(2) &
                 - BMie_bulk(3) * LFRC_BMie_bulk(1) ) / BMie_bulk(1)
!  original code
!   LFRC_BMie_bulk(1:3) = Mie1_bulk(1:3) - Mie2_bulk(1:3)

!  Coefficients, asymmetry parameter

   if ( Do_Expcoeffs ) then
      LFRC_BMie_asymm = DF_WW1 * ( Mie1_asymm - Mie2_asymm )
      do L = 0, BMie_ncoeffs
         LFRC_BMie_expcoeffs(1:6,L) = DF_WW1 * ( Mie1_expcoeffs(1:6,L) - Mie2_expcoeffs(1:6,L) )
      enddo
   endif

!  Fmatrix

   if ( Do_Fmatrix ) then
      do L = 1, n_Fmatrix_angles
         LFRC_BMie_Fmatrix(1:4,L) = DF_WW1 * ( Mie1_Fmatrix(1:4,L) - Mie2_Fmatrix(1:4,L) )
      enddo
   endif   

!  Finish

   return
end subroutine RTSMie_master_bimodal_plus

!  End module

end module RTSMie_master_bimodal_plus_m



