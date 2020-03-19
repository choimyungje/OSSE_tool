module tmat_master_bimodal_plus_m

!  This is the LINEARIZED Bimodal master #2 for Tmatrix code
!    ** RT Solutions, Version 1.0, 21 December 2010
!    ** RT Solutions, Version 1.1, 07 January  2011
!    ** RT Solutions, Version 1.2, 29 March    2011
!    ** RT Solutions, Version 1.3, 24 June     2011 (mono    control)
!    ** RT Solutions, Version 1.4, 30 June     2011 (Bimodal control)
!    ** RT Solutions, Version 1.5, 25 August   2011 (Bimodal control, more)

use tmat_parameters, only : fpk => tmat_fpkind, d_one, d_zero, NPL1, MAXNPA

use tmat_master_plus_m

!  Everything PUBLIC here
!  ----------------------

public

contains

subroutine tmat_master_bimodal_PLUS ( Tmat_Verbose, &
     Do_Expcoeffs, Do_Fmatrix,                   & ! Gp 1 Inputs (Flags)
     Do_Monodisperse, Do_EqSaSphere,             & ! Gp 1   Inputs (Flags)
     Do_LinearRef, Do_LinearEps, Do_LinearPSD,   & ! Gp 1 Inputs (Flags)
     Do_psd_OldStyle, psd_Index, psd_pars,       & ! Gp 1/2 Inputs (PSD)
     MonoRadius, R1, R2, FixR1R2, fraction,      & ! Gp 2   Inputs (PSD)
     np, nkmax, npna, ndgs, eps, accuracy,       & ! Gp 3 Inputs (General)
     lambda, n_real, n_imag,                     & ! Gp 4 Inputs (Optical)
     Btmat_bulk, Btmat_asymm, Btmat_ncoeffs,     & ! Outputs (Tmat)
     Btmat_expcoeffs, Btmat_Fmatrix,             & ! Outputs (Tmat)
     LPSD_Btmat_bulk, LPSD_Btmat_asymm,          & ! Outputs (LinTmat)
     LPSD_Btmat_expcoeffs, LPSD_Btmat_Fmatrix,   & ! Outputs (LinTmat)
     LRFE_Btmat_bulk, LRFE_Btmat_asymm,          & ! Outputs (LinTmat)
     LRFE_Btmat_expcoeffs, LRFE_Btmat_Fmatrix,   & ! Outputs (LinTmat)
     LFRC_Btmat_bulk, LFRC_Btmat_asymm,          & ! Outputs (Bimodal frac)
     LFRC_Btmat_expcoeffs, LFRC_Btmat_Fmatrix,   & ! Outputs (Bimodal frac)
     NLIN, Tmat_dist, LPSD_Tmat_dist,            & ! Outputs (PSD)
     fail, istatus, message, trace, trace_2, trace_3 )     ! Outputs (status)

!  List of Inputs
!  ==============

!  Flag inputs
!  -----------

!      Do_Expcoeffs      - Boolean flag for computing Expansion Coefficients
!      Do_Fmatrix        - Boolean flag for computing F-matrix at equal-angles

!      Do_Monodisperse   - Boolean flag for Doing a Monodisperse calculation
!                          If set, the PSD stuff will be turned off internally

!      Do_EqSaSphere     - Boolean flag for specifying particle size in terms
!                          of the  equal-surface-area-sphere radius

!      Do_psd_OldStyle   - Boolean flag for using original PSD specifications

!      Do_LinearRef      - Boolean Flag for doing Refractive Index linearization
!      Do_LinearEps      - Boolean Flag for doing (Aspect ratio)   linearization

!      Do_LinearPSD      - Boolean Flag for doing PSD linearization

!  General inputs
!  --------------

!      NKMAX.LE.988 is such that NKMAX+2 is the                        
!           number of Gaussian quadrature points used in               
!           integrating over the size distribution for particles

!      NDGS - parameter controlling the number of division points      
!             in computing integrals over the particle surface.        
!             For compact particles, the recommended value is 2.       
!             For highly aspherical particles larger values (3, 4,...) 
!             may be necessary to obtain convergence.                  
!             The code does not check convergence over this parameter. 
!             Therefore, control comparisons of results obtained with  
!             different NDGS-values are recommended.

!      NPNA - number of equidistant scattering angles (from 0      
!             to 180 deg) for which the scattering matrix is           
!             calculated.                                              
            
!      EPS and NP - specify the shape of the particles.                
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

!  optical inputs
!  --------------

!      LAMBDA         - wavelength of light (microns)
!      N_REAL, N_IMAG - real and imaginary parts, refractive index (N-i.GE.0)   

!  PSD inputs
!  ----------

!      psd_Index      - Index for particle size distribution of spheres
!      psd_pars       - Parameters characterizing PSD (up to 3 allowed)

!      Monoradius     - Monodisperse radius size (Microns)

!      R1, R2         - Minimum and Maximum radii (Microns)
!      FixR1R2        - Boolean flag for allowing internal calculation of R1/R2

   implicit none

!  Boolean Input arguments
!  -----------------------

!  Verbose flag now passed, 10/19/16

   logical  , intent(in)  :: Tmat_Verbose

!  Flags for Expansion Coefficient and Fmatrix calculations

   logical, intent(in)  :: Do_Expcoeffs
   logical, intent(in)  :: Do_Fmatrix

!  Logical flag for using equal-surface-area sepcification

   logical, intent(in)  :: Do_EqSaSphere

!  Logical flag for Monodisperse calculation

   logical, intent(in)  :: Do_monodisperse

!  linearization of other quantitites

   logical, intent(in)  :: Do_LinearRef
   logical, intent(in)  :: Do_LinearEps

!  PSD Linearization
!    * This is checked and re-set (if required) for Monodisperse case

   logical  , intent(inout)  :: Do_LinearPSD

!  Style flag
!    * This is checked and re-set (if required) for Monodisperse case

   logical  , intent(inout)  :: Do_psd_OldStyle

!  General Input arguments
!  -----------------------

!  integers (nkmax may be re-set for Monodisperse case)

   integer  , intent(in)     :: np, ndgs(2), npna
   integer  , intent(inout)  :: nkmax(2)

!  Accuracy and aspect ratio

   real(fpk), intent(in)  :: accuracy, eps(2)

!  Optical: Wavelength, refractive index
!  -------------------------------------

   real(fpk), intent(in)  :: lambda, n_real(2), n_imag(2)

!  PSD inputs
!  ----------

!  Flag for making an internal Fix of R1 and R2
!    ( Not relevant for the Old distribution

   logical, intent(inout)  :: FixR1R2(2)

!  Monodisperse radius (input)

   real(fpk), intent(in)   :: Monoradius

!  R1 and R2 (intent(inout))

   real(fpk), intent(inout)  :: R1(2), R2(2)

!  PSD index and parameters

   integer  , intent(in)   :: psd_Index(2)
   real(fpk), intent(in)   :: psd_pars (3,2)

!  Fraction

   real(fpk), intent(in)   :: fraction

!  Output arguments --------> BIMODAL
!  ================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk), intent(out) :: BTmat_bulk (3)

!  linearizations w.r.t. PSD parameters

   real(fpk), intent(out) :: LPSD_BTmat_bulk (3,3,2)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(fpk), intent(out) :: LRFE_BTmat_bulk (3,3,2)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

!  Regular quantities

   integer  , intent(out) :: Btmat_ncoeffs
   real(fpk), intent(out) :: BTmat_expcoeffs (NPL1,6)
   real(fpk), intent(out) :: BTmat_asymm

!  linearizations w.r.t. PSD parameters

   real(fpk), intent(out) :: LPSD_BTmat_expcoeffs (NPL1,6,3,2)
   real(fpk), intent(out) :: LPSD_BTmat_asymm(3,2)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(fpk), intent(out) :: LRFE_BTmat_expcoeffs (NPL1,6,3,2)
   real(fpk), intent(out) :: LRFE_BTmat_asymm(3,2)

!  F-matrix,  optional output
!  --------------------------

!  F-matrix

   real(fpk), intent(out) :: BTmat_Fmatrix (MAXNPA,6)

!  Linearizations of F-matrix

   real(fpk), intent(out) :: LPSD_BTmat_Fmatrix (MAXNPA,6,3,2)
   real(fpk), intent(out) :: LRFE_BTmat_Fmatrix (MAXNPA,6,3,2)

!  Fraction Jacobian
!  -----------------

   real(fpk), intent(out)  :: LFRC_BTmat_bulk (3)
   real(fpk), intent(out)  :: LFRC_BTmat_expcoeffs (NPL1,6)
   real(fpk), intent(out)  :: LFRC_BTmat_asymm
   real(fpk), intent(out)  :: LFRC_BTmat_Fmatrix (MAXNPA,6)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(fpk), intent(out) :: Tmat_dist (5,2)
   real(fpk), intent(out) :: LPSD_Tmat_dist (5,3,2)

!  Number of RFE parameters

   integer       , intent(out) :: nlin

!  Exception handling
!  ------------------

   logical       , intent(out) :: fail
   integer       , intent(out) :: istatus
   character*(*) , intent(out) :: message
   character*(*) , intent(out) :: trace
   character*(*) , intent(out) :: trace_2
   character*(*) , intent(out) :: trace_3

!  Local Arrays
!  ============

!  Bulk distribution parameters
!  Expansion coefficients and Asymmetry parameter
!  F-matrix,  optional output

   real(fpk) :: Tmat1_bulk (3)
   integer   :: Tmat1_ncoeffs
   real(fpk) :: Tmat1_expcoeffs (NPL1,6)
   real(fpk) :: Tmat1_asymm
   real(fpk) :: Tmat1_Fmatrix (MAXNPA,6)

   real(fpk) :: Tmat2_bulk (3)
   integer   :: Tmat2_ncoeffs
   real(fpk) :: Tmat2_expcoeffs (NPL1,6)
   real(fpk) :: Tmat2_asymm
   real(fpk) :: Tmat2_Fmatrix (MAXNPA,6)


   real(fpk) :: LPSD_Tmat1_bulk (3,3)
   real(fpk) :: LPSD_Tmat1_Fmatrix (MAXNPA,6,3)
   real(fpk) :: LPSD_Tmat1_expcoeffs (NPL1,6,3)
   real(fpk) :: LPSD_Tmat1_asymm(3)

   real(fpk) :: LRFE_Tmat1_bulk (3,3)
   real(fpk) :: LRFE_Tmat1_Fmatrix (MAXNPA,6,3)
   real(fpk) :: LRFE_Tmat1_expcoeffs (NPL1,6,3)
   real(fpk) :: LRFE_Tmat1_asymm(3)

   real(fpk) :: LPSD_Tmat2_bulk (3,3)
   real(fpk) :: LPSD_Tmat2_Fmatrix (MAXNPA,6,3)
   real(fpk) :: LPSD_Tmat2_expcoeffs (NPL1,6,3)
   real(fpk) :: LPSD_Tmat2_asymm(3)

   real(fpk) :: LRFE_Tmat2_bulk (3,3)
   real(fpk) :: LRFE_Tmat2_Fmatrix (MAXNPA,6,3)
   real(fpk) :: LRFE_Tmat2_expcoeffs (NPL1,6,3)
   real(fpk) :: LRFE_Tmat2_asymm(3)

!  Other local variables
!  ---------------------

   integer   :: k, L, q
   real(fpk) :: FF1, FF2, WW1, WW2, TERM1, TERM2
   real(fpk) :: Csca1_FF1, Csca2_FF2, Csca_total, D1_Csca1(3), D2_Csca2(3)
   real(fpk) :: D1_WW1(3), D2_WW1(3), D1_WW2(3), D2_WW2(3), DF_WW1

!  Zero the output
!  ---------------

   BTmat_bulk      = d_zero
   BTmat_Fmatrix   = d_zero
   BTmat_expcoeffs = d_zero
   BTmat_asymm     = d_zero
   BTmat_ncoeffs   = 0

   LRFE_BTmat_bulk      = d_zero
   LRFE_BTmat_Fmatrix   = d_zero
   LRFE_BTmat_expcoeffs = d_zero
   LRFE_BTmat_asymm     = d_zero

   LPSD_BTmat_bulk      = d_zero
   LPSD_BTmat_Fmatrix   = d_zero
   LPSD_BTmat_expcoeffs = d_zero
   LPSD_BTmat_asymm     = d_zero

   LFRC_BTmat_bulk      = d_zero
   LFRC_BTmat_Fmatrix   = d_zero
   LFRC_BTmat_expcoeffs = d_zero
   LFRC_BTmat_asymm     = d_zero

   Tmat_dist       = d_zero
   LPSD_Tmat_dist  = d_zero

   FF1 = fraction
   FF2 = d_one - FF1

   trace_3 = ' '

!  Check: No Monodisperse here !
!  -----------------------------

   if ( do_monodisperse ) then
      fail = .true.; istatus = 2
      trace_3 = 'tmat_master_bimodal_PLUS module: Input error: MONODISPERSE FLAG must be Off!'
      return
   endif

!  First Call
!  ---------

   k = 1 ; write(*,*)' ** Doing Tmatrix PLUS for PSD # 1 ----------------------'
   call tmat_master_PLUS ( Tmat_Verbose,             &
     Do_Expcoeffs, Do_Fmatrix,                       & ! Gp 1 Inputs (Flags)
     Do_Monodisperse, Do_EqSaSphere,                 & ! Gp 1   Inputs (Flags)
     Do_LinearRef, Do_LinearEps, Do_LinearPSD,       & ! Gp 1 Inputs (Flags)
     Do_psd_OldStyle, psd_Index(k), psd_pars(:,k),   & ! Gp 1/2 Inputs (PSD)
     MonoRadius, R1(k), R2(k), FixR1R2(k),           & ! Gp 2   Inputs (PSD)
     np, nkmax(k), npna, ndgs(k), eps(k), accuracy,  & ! Gp 3   Inputs (General)
     lambda, n_real(k), n_imag(k),                   & ! Gp 4   Inputs (Optical)
     tmat1_bulk, tmat1_asymm, tmat1_ncoeffs,         & ! Outputs (Tmat)
     tmat1_expcoeffs, tmat1_Fmatrix,                 & ! Outputs (Tmat)
     LPSD_tmat1_bulk, LPSD_tmat1_asymm,              & ! Outputs (LinTmat)
     LPSD_tmat1_expcoeffs, LPSD_tmat1_Fmatrix,       & ! Outputs (LinTmat)
     LRFE_tmat1_bulk, LRFE_tmat1_asymm,              & ! Outputs (LinTmat)
     LRFE_tmat1_expcoeffs, LRFE_tmat1_Fmatrix,       & ! Outputs (LinTmat)
     NLIN, Tmat_dist(:,k), LPSD_Tmat_dist(:,:,k),    & ! Outputs (PSD)
     fail, istatus, message, trace, trace_2 )          ! Outputs (status)

!  Exception handling

   if ( fail ) then
      trace_3 = 'tmat_master_bimodal_PLUS module: First PSD call, Warning or Error'
      if ( Istatus .eq. 2 ) return
   endif

!  Second call
!  -----------

   k = 2 ; write(*,*)' ** Doing Tmatrix PLUS for PSD # 2 ----------------------'
   call tmat_master_PLUS ( Tmat_Verbose,             &
     Do_Expcoeffs, Do_Fmatrix,                       & ! Gp 1 Inputs (Flags)
     Do_Monodisperse, Do_EqSaSphere,                 & ! Gp 1   Inputs (Flags)
     Do_LinearRef, Do_LinearEps, Do_LinearPSD,       & ! Gp 1 Inputs (Flags)
     Do_psd_OldStyle, psd_Index(k), psd_pars(:,k),   & ! Gp 1/2 Inputs (PSD)
     MonoRadius, R1(k), R2(k), FixR1R2(k),           & ! Gp 2   Inputs (PSD)
     np, nkmax(k), npna, ndgs(k), eps(k), accuracy,  & ! Gp 3   Inputs (General)
     lambda, n_real(k), n_imag(k),                   & ! Gp 4   Inputs (Optical)
     tmat2_bulk, tmat2_asymm, tmat2_ncoeffs,         & ! Outputs (Tmat)
     tmat2_expcoeffs, tmat2_Fmatrix,                 & ! Outputs (Tmat)
     LPSD_tmat2_bulk, LPSD_tmat2_asymm,              & ! Outputs (LinTmat)
     LPSD_tmat2_expcoeffs, LPSD_tmat2_Fmatrix,       & ! Outputs (LinTmat)
     LRFE_tmat2_bulk, LRFE_tmat2_asymm,              & ! Outputs (LinTmat)
     LRFE_tmat2_expcoeffs, LRFE_tmat2_Fmatrix,       & ! Outputs (LinTmat)
     NLIN, Tmat_dist(:,k), LPSD_Tmat_dist(:,:,k),    & ! Outputs (PSD)
     fail, istatus, message, trace, trace_2 )          ! Outputs (status)

!  Exception handling

   if ( fail ) then
      trace_3 = 'tmat_master_bimodal_PLUS module: Second PSD call, Warning or Error'
      if ( Istatus .eq. 2 ) return
   endif

!  Bimodal determination
!  =====================

!  Revision 20 September 2011
!    Correct definition for Expcoeffs/Fmatrix: WW1/WW2 in place of FF1/FF2

   Csca1_FF1  = FF1 * Tmat1_bulk(2)   
   Csca2_FF2  = FF2 * Tmat2_bulk(2)   
   Csca_total =  Csca1_FF1 + Csca2_FF2
   WW1   = Csca1_FF1 / Csca_total
   WW2   = Csca2_FF2 / Csca_total

!  Bimodal quantities

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong
   BTmat_bulk(1:2) = FF1 * Tmat1_bulk(1:2) + FF2 * Tmat2_bulk(1:2)
   BTmat_bulk(3)   = BTmat_bulk(2) / BTmat_bulk(1)
!  original code
!   BTmat_bulk(1:3) = FF1 * Tmat1_bulk(1:3) + FF2 * Tmat2_bulk(1:3)

   if ( Do_Expcoeffs ) then
      BTmat_asymm = WW1 * Tmat1_asymm + WW2 * Tmat2_asymm
      BTmat_ncoeffs = max(tmat1_ncoeffs,tmat2_ncoeffs)
      do L = 1, min(tmat1_ncoeffs,tmat2_ncoeffs)
         Btmat_expcoeffs(L,1:6) = WW1 * tmat1_expcoeffs(L,1:6) + & 
                                  WW2 * tmat2_expcoeffs(L,1:6)
      enddo
      if ( tmat1_ncoeffs .lt. tmat2_ncoeffs ) then
          do L = tmat1_ncoeffs + 1,tmat2_ncoeffs
              Btmat_expcoeffs(L,1:6) = WW2 * tmat2_expcoeffs(L,1:6)
          enddo
      else if ( tmat1_ncoeffs .gt. tmat2_ncoeffs ) then
          do L = tmat2_ncoeffs + 1,tmat1_ncoeffs
              Btmat_expcoeffs(L,1:6) = WW1 * tmat1_expcoeffs(L,1:6)
          enddo
      endif
   endif
   if ( Do_Fmatrix ) then
      do L = 1, npna
         Btmat_Fmatrix(L,1:6) = WW1 * tmat1_Fmatrix(L,1:6) + & 
                                WW2 * tmat2_Fmatrix(L,1:6)
      enddo
   endif   

!  Revision 20 September 2011
!    Correct definition for Expcoeffs/Fmatrix: WW1/WW2 in place of FF1/FF2

   Csca1_FF1  = FF1 * Tmat1_bulk(2)   
   Csca2_FF2  = FF2 * Tmat2_bulk(2)   
   Csca_total =  Csca1_FF1 + Csca2_FF2
   WW1   = Csca1_FF1 / Csca_total
   WW2   = Csca2_FF2 / Csca_total

!  LRFE Linearizations
!  ===================

   if ( nlin .gt. 0 ) then

!  help variables
!  --------------

      D1_Csca1(1:nlin) = LRFE_Tmat1_bulk(2,1:nlin)
      D2_Csca2(1:nlin) = LRFE_Tmat2_bulk(2,1:nlin)
      D1_WW1(1:nlin) = FF1 * D1_Csca1(1:nlin) * WW2 / Csca_total
      D2_WW2(1:nlin) = FF2 * D2_Csca2(1:nlin) * WW1 / Csca_total
      D1_WW2(1:nlin) = - D1_WW1(1:nlin)
      D2_WW1(1:nlin) = - D2_WW2(1:nlin)

!  Bulk variables
!  --------------

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong
      LRFE_BTmat_bulk(1:2,1:nlin,1) = FF1 * LRFE_Tmat1_bulk(1:2,1:nlin)
      LRFE_BTmat_bulk(1:2,1:nlin,2) = FF2 * LRFE_Tmat2_bulk(1:2,1:nlin)
      LRFE_BTmat_bulk(3,1:nlin,1) = ( LRFE_BTmat_bulk(2,1:nlin,1) &
                 - BTmat_bulk(3) * LRFE_BTmat_bulk(1,1:nlin,1) ) / BTmat_bulk(1)
      LRFE_BTmat_bulk(3,1:nlin,2) = ( LRFE_BTmat_bulk(2,1:nlin,2) &
                 - BTmat_bulk(3) * LRFE_BTmat_bulk(1,1:nlin,2) ) / BTmat_bulk(1)
!  original code
!      LRFE_BTmat_bulk(1:3,1:nlin,1) = FF1 * LRFE_Tmat1_bulk(1:3,1:nlin)
!      LRFE_BTmat_bulk(1:3,1:nlin,2) = FF2 * LRFE_Tmat2_bulk(1:3,1:nlin)


!  Asymmetry parameter
!  -------------------

      if ( Do_Expcoeffs ) then
         do q = 1, nlin
            TERM1 = WW1 * LRFE_Tmat1_asymm(q)
            TERM2 = D1_WW1(q) * Tmat1_asymm + D1_WW2(q) * Tmat2_asymm
            LRFE_BTmat_asymm(q,1) = TERM1 + TERM2
            TERM1 = WW2 * LRFE_Tmat2_asymm(q)
            TERM2 = D2_WW1(q) * Tmat1_asymm + D2_WW2(q) * Tmat2_asymm
            LRFE_BTmat_asymm(q,2) = TERM1 + TERM2
         enddo
      endif

!  Expansion coefficients
!  ----------------------

      if ( Do_Expcoeffs ) then

!  w.r.t  microphysical parameters for Particles in Mode 1

         do q = 1, nlin
            do L = 1, min(tmat1_ncoeffs,tmat2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW1 * LRFE_tmat1_expcoeffs(L,k,q)
                  TERM2 = D1_WW1(q) * Tmat1_expcoeffs(L,k) + D1_WW2(q) * Tmat2_expcoeffs(L,k)
                  LRFE_Btmat_expcoeffs(L,k,q,1) = TERM1 + TERM2
               enddo
            enddo
            if ( tmat1_ncoeffs .lt. tmat2_ncoeffs ) then
               do L = tmat1_ncoeffs + 1,tmat2_ncoeffs
                  do k = 1, 6
                     LRFE_Btmat_expcoeffs(L,k,q,1) = D1_WW2(q) * tmat2_expcoeffs(L,k)
                  enddo
               enddo  
            else if ( tmat1_ncoeffs .gt. tmat2_ncoeffs ) then
               do L = tmat2_ncoeffs + 1,tmat1_ncoeffs
                  do k = 1, 6
                     TERM1 = WW1 * LRFE_tmat1_expcoeffs(L,k,q)
                     LRFE_Btmat_expcoeffs(L,k,q,1) = TERM1 + D1_WW1(q) * tmat1_expcoeffs(L,k)
                  enddo
               enddo  
            endif
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 1, min(tmat1_ncoeffs,tmat2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW2 * LRFE_tmat2_expcoeffs(L,k,q)
                  TERM2 = D2_WW1(q) * Tmat1_expcoeffs(L,k) + D2_WW2(q) * Tmat2_expcoeffs(L,k)
                  LRFE_Btmat_expcoeffs(L,k,q,2) = TERM1 + TERM2
               enddo
            enddo
            if ( tmat1_ncoeffs .lt. tmat2_ncoeffs ) then
               do L = tmat1_ncoeffs + 1,tmat2_ncoeffs
                  do k = 1, 6
                     TERM1 = WW2 * LRFE_tmat2_expcoeffs(L,k,q)
                     LRFE_Btmat_expcoeffs(L,k,q,2) = TERM1 + D2_WW2(q) * tmat2_expcoeffs(L,k)
                  enddo
               enddo  
            else if ( tmat1_ncoeffs .gt. tmat2_ncoeffs ) then
               do L = tmat2_ncoeffs + 1,tmat1_ncoeffs
                  do k = 1, 6
                     LRFE_Btmat_expcoeffs(L,k,q,2) = D2_WW1(q) * tmat1_expcoeffs(L,k)
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
            do L = 1, npna
               do k = 1, 6
                  TERM1 = WW1 * LRFE_tmat1_Fmatrix(L,k,q)
                  TERM2 = D1_WW1(q) * Tmat1_Fmatrix(L,k) + D1_WW2(q) * Tmat2_Fmatrix(L,k)
                  LRFE_Btmat_Fmatrix(L,k,q,1) = TERM1 + TERM2
               enddo
            enddo
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 1, npna
               do k = 1, 6
                  TERM1 = WW2 * LRFE_tmat2_Fmatrix(L,k,q)
                  TERM2 = D2_WW1(q) * Tmat1_Fmatrix(L,k) + D2_WW2(q) * Tmat2_Fmatrix(L,k)
                  LRFE_Btmat_Fmatrix(L,k,q,2) = TERM1 + TERM2
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

!  @@@ Rob Fix 05June 14, Ordering was wrong for D2_WW1

      nlin = 3
      D1_Csca1(1:nlin) = LPSD_Tmat1_bulk(2,1:nlin)
      D2_Csca2(1:nlin) = LPSD_Tmat2_bulk(2,1:nlin)
      D1_WW1(1:nlin) = FF1 * D1_Csca1(1:nlin) * WW2 / Csca_total
      D2_WW2(1:nlin) = FF2 * D2_Csca2(1:nlin) * WW1 / Csca_total
      D1_WW2(1:nlin) = - D1_WW1(1:nlin)
      D2_WW1(1:nlin) = - D2_WW2(1:nlin)

!  Bulk variables
!  --------------

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong
      LPSD_BTmat_bulk(1:2,1:nlin,1) = FF1 * LPSD_Tmat1_bulk(1:2,1:nlin)
      LPSD_BTmat_bulk(1:2,1:nlin,2) = FF2 * LPSD_Tmat2_bulk(1:2,1:nlin)
      LPSD_BTmat_bulk(3,1:nlin,1) = ( LPSD_BTmat_bulk(2,1:nlin,1) &
                 - BTmat_bulk(3) * LPSD_BTmat_bulk(1,1:nlin,1) ) / BTmat_bulk(1)
      LPSD_BTmat_bulk(3,1:nlin,2) = ( LPSD_BTmat_bulk(2,1:nlin,2) &
                 - BTmat_bulk(3) * LPSD_BTmat_bulk(1,1:nlin,2) ) / BTmat_bulk(1)
!  original code
!      LPSD_BTmat_bulk(1:3,1:nlin,1) = FF1 * LPSD_Tmat1_bulk(1:3,1:nlin)
!      LPSD_BTmat_bulk(1:3,1:nlin,2) = FF2 * LPSD_Tmat2_bulk(1:3,1:nlin)

!  Asymmetry parameter
!  -------------------

      if ( Do_Expcoeffs ) then
         do q = 1, nlin
            TERM1 = WW1 * LPSD_Tmat1_asymm(q)
            TERM2 = D1_WW1(q) * Tmat1_asymm + D1_WW2(q) * Tmat2_asymm
            LPSD_BTmat_asymm(q,1) = TERM1 + TERM2
            TERM1 = WW2 * LPSD_Tmat2_asymm(q)
            TERM2 = D2_WW1(q) * Tmat1_asymm + D2_WW2(q) * Tmat2_asymm
            LPSD_BTmat_asymm(q,2) = TERM1 + TERM2
         enddo
      endif

!  Expansion coefficients
!  ----------------------

      if ( Do_Expcoeffs ) then

!  w.r.t  microphysical parameters for Particles in Mode 1

         do q = 1, nlin
            do L = 1, min(tmat1_ncoeffs,tmat2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW1 * LPSD_tmat1_expcoeffs(L,k,q)
                  TERM2 = D1_WW1(q) * Tmat1_expcoeffs(L,k) + D1_WW2(q) * Tmat2_expcoeffs(L,k)
                  LPSD_Btmat_expcoeffs(L,k,q,1) = TERM1 + TERM2
               enddo
            enddo
            if ( tmat1_ncoeffs .lt. tmat2_ncoeffs ) then
               do L = tmat1_ncoeffs + 1,tmat2_ncoeffs
                  do k = 1, 6
                     LPSD_Btmat_expcoeffs(L,k,q,1) = D1_WW2(q) * tmat2_expcoeffs(L,k)
                  enddo
               enddo  
            else if ( tmat1_ncoeffs .gt. tmat2_ncoeffs ) then
               do L = tmat2_ncoeffs + 1,tmat1_ncoeffs
                  do k = 1, 6
                     TERM1 = WW1 * LPSD_tmat1_expcoeffs(L,k,q)
                     LPSD_Btmat_expcoeffs(L,k,q,1) = TERM1 + D1_WW1(q) * tmat1_expcoeffs(L,k)
                  enddo
               enddo  
            endif
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 1, min(tmat1_ncoeffs,tmat2_ncoeffs)
               do k = 1, 6
                  TERM1 = WW2 * LPSD_tmat2_expcoeffs(L,k,q)
                  TERM2 = D2_WW1(q) * Tmat1_expcoeffs(L,k) + D2_WW2(q) * Tmat2_expcoeffs(L,k)
                  LPSD_Btmat_expcoeffs(L,k,q,2) = TERM1 + TERM2
               enddo
            enddo
            if ( tmat1_ncoeffs .lt. tmat2_ncoeffs ) then
               do L = tmat1_ncoeffs + 1,tmat2_ncoeffs
                  do k = 1, 6
                     TERM1 = WW2 * LPSD_tmat2_expcoeffs(L,k,q)
                     LPSD_Btmat_expcoeffs(L,k,q,2) = TERM1 + D2_WW2(q) * tmat2_expcoeffs(L,k)
                  enddo
               enddo  
            else if ( tmat1_ncoeffs .gt. tmat2_ncoeffs ) then
               do L = tmat2_ncoeffs + 1,tmat1_ncoeffs
                  do k = 1, 6
                     LPSD_Btmat_expcoeffs(L,k,q,2) = D2_WW1(q) * tmat1_expcoeffs(L,k)
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
            do L = 1, npna
               do k = 1, 6
                  TERM1 = WW1 * LPSD_tmat1_Fmatrix(L,k,q)
                  TERM2 = D1_WW1(q) * Tmat1_Fmatrix(L,k) + D1_WW2(q) * Tmat2_Fmatrix(L,k)
                  LPSD_Btmat_Fmatrix(L,k,q,1) = TERM1 + TERM2
               enddo
            enddo
         enddo

!  w.r.t  microphysical parameters for Particles in Mode 2

         do q = 1, nlin
            do L = 1, npna
               do k = 1, 6
                  TERM1 = WW2 * LPSD_tmat2_Fmatrix(L,k,q)
                  TERM2 = D2_WW1(q) * Tmat1_Fmatrix(L,k) + D2_WW2(q) * Tmat2_Fmatrix(L,k)
                  LPSD_Btmat_Fmatrix(L,k,q,2) = TERM1 + TERM2
               enddo
            enddo
         enddo

!  End clauses

      endif   
   endif

!  Fractional linearization. NOT NORMALIZED
!  ========================================

!  Help variables

   DF_WW1 = ( WW1 * Tmat2_bulk(2) + WW2 * Tmat1_bulk(2) ) / Csca_total 

!  Bulk

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong
   LFRC_BTmat_bulk(1:2) = Tmat1_bulk(1:2) - Tmat2_bulk(1:2)
   LFRC_BTmat_bulk(3) = ( LFRC_BTmat_bulk(2) &
                 - BTmat_bulk(3) * LFRC_BTmat_bulk(1) ) / BTmat_bulk(1)
!  original code
!   LFRC_BTmat_bulk(1:3) = Tmat1_bulk(1:3) - Tmat2_bulk(1:3)

!  Coefficients, asymmetry parameter

   if ( Do_Expcoeffs ) then
      LFRC_BTmat_asymm = DF_WW1 * ( Tmat1_asymm - Tmat2_asymm )
      do L = 1, BTmat_ncoeffs
         LFRC_Btmat_expcoeffs(L,1:6) = DF_WW1 * ( tmat1_expcoeffs(L,1:6) - tmat2_expcoeffs(L,1:6) )
      enddo
   endif

!  Fmatrix

   if ( Do_Fmatrix ) then
      do L = 1, npna
         LFRC_Btmat_Fmatrix(L,1:6) = DF_WW1 * ( tmat1_Fmatrix(L,1:6) - tmat2_Fmatrix(L,1:6) )
      enddo
   endif   

!  Finish

   return
end subroutine tmat_master_bimodal_PLUS

!  End module

end module tmat_master_bimodal_plus_m
