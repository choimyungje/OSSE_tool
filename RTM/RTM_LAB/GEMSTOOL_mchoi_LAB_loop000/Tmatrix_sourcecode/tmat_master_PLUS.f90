module tmat_master_plus_m

!  This is the LINEARIZED master #2 for Tmatrix code
!    ** RT Solutions, Version 1.0, 21 December 2010
!    ** RT Solutions, Version 1.1, 07 January  2011
!    ** RT Solutions, Version 1.2, 29 March    2011
!    ** RT Solutions, Version 1.3, 24 June     2011 (mono control)

use tmat_parameters, only : NPNG1, NPNG2, NPN2, NPN1, NPN4, NPN6, MAXNPA, &
                            NPL, NPL1, fpk => tmat_fpkind, d_one

use tmat_distributions
use tmat_functions
use tmat_functions_plus
use tmat_makers
use tmat_makers_plus
use tmat_scattering
use tmat_scattering_plus

!  Everything PUBLIC here
!  ----------------------

public

contains

subroutine tmat_master_PLUS ( Tmat_verbose,    &
     Do_Expcoeffs, Do_Fmatrix,                 & ! Gp 1 Inputs (Flags)
     Do_Monodisperse, Do_EqSaSphere,           & ! Gp 1   Inputs (Flags)
     Do_LinearRef, Do_LinearEps, Do_LinearPSD, & ! Gp 1 Inputs (Flags)
     Do_psd_OldStyle, psd_Index, psd_pars,     & ! Gp 1/2 Inputs (PSD)
     MonoRadius, R1, R2, FixR1R2,              & ! Gp 2   Inputs (PSD)
     np, nkmax, npna, ndgs, eps, accuracy,     & ! Gp 3 Inputs (General)
     lambda, n_real, n_imag,                   & ! Gp 4 Inputs (Optical)
     tmat_bulk, tmat_asymm, tmat_ncoeffs,      & ! Outputs (Tmat)
     tmat_expcoeffs, tmat_Fmatrix,             & ! Outputs (Tmat)
     LPSD_tmat_bulk, LPSD_tmat_asymm,          & ! Outputs (LinTmat)
     LPSD_tmat_expcoeffs, LPSD_tmat_Fmatrix,   & ! Outputs (LinTmat)
     LRFE_tmat_bulk, LRFE_tmat_asymm,          & ! Outputs (LinTmat)
     LRFE_tmat_expcoeffs, LRFE_tmat_Fmatrix,   & ! Outputs (LinTmat)
     NLIN, Tmat_dist, LPSD_Tmat_dist,          & ! Outputs (PSD)
     fail, istatus, message, trace, trace_2 )    ! Outputs (status)

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

!  Verbose output flag

   logical  , intent(in)  :: Tmat_verbose

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

   integer  , intent(in)     :: np, ndgs, npna
   integer  , intent(inout)  :: nkmax

!  Accuracy and aspect ratio

   real(fpk), intent(in)  :: accuracy, eps

!  Optical: Wavelength, refractive index
!  -------------------------------------

   real(fpk), intent(in)  :: lambda, n_real, n_imag

!  PSD inputs
!  ----------

!  Flag for making an internal Fix of R1 and R2
!    ( Not relevant for the Old distribution

   logical, intent(inout)  :: FixR1R2

!  Monodisperse radius (input)

   real(fpk), intent(in)   :: Monoradius

!  R1 and R2 (intent(inout))

   real(fpk), intent(inout)  :: R1, R2

!  PSD index and parameters

   integer  , intent(in)   :: psd_Index
   real(fpk), intent(in)   :: psd_pars (3)

!  Output arguments
!  ================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk), intent(out) :: Tmat_bulk (3)

!  linearizations w.r.t. PSD parameters

   real(fpk), intent(out) :: LPSD_Tmat_bulk (3,3)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(fpk), intent(out) :: LRFE_Tmat_bulk (3,3)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

!  Regular quantities

   integer  , intent(out) :: tmat_ncoeffs
   real(fpk), intent(out) :: Tmat_expcoeffs (NPL1,6)
   real(fpk), intent(out) :: Tmat_asymm

!  linearizations w.r.t. PSD parameters

   real(fpk), intent(out) :: LPSD_Tmat_expcoeffs (NPL1,6,3)
   real(fpk), intent(out) :: LPSD_Tmat_asymm(3)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(fpk), intent(out) :: LRFE_Tmat_expcoeffs (NPL1,6,3)
   real(fpk), intent(out) :: LRFE_Tmat_asymm(3)

!  F-matrix,  optional output
!  --------------------------

!  F-matrix

   real(fpk), intent(out) :: Tmat_Fmatrix (MAXNPA,6)

!  Linearizations of F-matrix

   real(fpk), intent(out) :: LPSD_Tmat_Fmatrix (MAXNPA,6,3)
   real(fpk), intent(out) :: LRFE_Tmat_Fmatrix (MAXNPA,6,3)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(fpk), intent(out) :: Tmat_dist (5)
   real(fpk), intent(out) :: LPSD_Tmat_dist (5,3)

!  Number of RFE parameters

   integer       , intent(out) :: nlin

!  Exception handling
!  ------------------

   logical       , intent(out) :: fail
   integer       , intent(out) :: istatus
   character*(*) , intent(out) :: message
   character*(*) , intent(out) :: trace
   character*(*) , intent(out) :: trace_2

!  Local Arrays
!  ============

!  Constants 

   real(fpk)  :: AN(NPN1),ANN(NPN1,NPN1)

! Various quadratures

   real(fpk)  :: X(NPNG2),W(NPNG2)
   REAL(fpk)  :: Le_X(NPNG2), Le_W(NPNG2)
   real(fpk)  :: S(NPNG2),SS(NPNG2)
   real(fpk)  :: Le_S(NPNG2),Le_SS(NPNG2)

   real(fpk)  :: XG(1000) ,WG(1000)
   real(fpk)  :: XG1(2000),WG1(2000)
   real(fpk)  :: WG1_derivs(2000,3)

!  R and DR functions and linearizations

   REAL(fpk)  :: R(NPNG2), DR(NPNG2)
   real(fpk)  :: DDR(NPNG2), DRR(NPNG2), DRI(NPNG2)

   real(fpk)  :: Le_R(NPNG2), Le_DR(NPNG2)
   real(fpk)  :: L_DDR(3,NPNG2), L_DRR(3,NPNG2), L_DRI(3,NPNG2)

!  Bessel Master output

   real(fpk)  :: J_BESS  (NPNG2,NPN1)
   real(fpk)  :: Y_BESS  (NPNG2,NPN1)
   real(fpk)  :: JR_BESS (NPNG2,NPN1)
   real(fpk)  :: JI_BESS (NPNG2,NPN1)
   real(fpk)  :: DJ_BESS (NPNG2,NPN1)
   real(fpk)  :: DY_BESS (NPNG2,NPN1)
   real(fpk)  :: DJR_BESS(NPNG2,NPN1)
   real(fpk)  :: DJI_BESS(NPNG2,NPN1)

!  Linearized Bessel output

   real(fpk)  :: L_J_BESS  (3,NPNG2,NPN1)
   real(fpk)  :: L_Y_BESS  (3,NPNG2,NPN1)
   real(fpk)  :: L_JR_BESS (3,NPNG2,NPN1)
   real(fpk)  :: L_JI_BESS (3,NPNG2,NPN1)
   real(fpk)  :: L_DJ_BESS (3,NPNG2,NPN1)
   real(fpk)  :: L_DY_BESS (3,NPNG2,NPN1)
   real(fpk)  :: L_DJR_BESS(3,NPNG2,NPN1)
   real(fpk)  :: L_DJI_BESS(3,NPNG2,NPN1)

!  Miscellaneous output

   REAL(fpk)  :: R11(NPN1,NPN1),R12(NPN1,NPN1)
   REAL(fpk)  :: R21(NPN1,NPN1),R22(NPN1,NPN1)
   REAL(fpk)  :: I11(NPN1,NPN1),I12(NPN1,NPN1)
   REAL(fpk)  :: I21(NPN1,NPN1),I22(NPN1,NPN1)

   REAL(fpk)  :: RG11(NPN1,NPN1),RG12(NPN1,NPN1)
   REAL(fpk)  :: RG21(NPN1,NPN1),RG22(NPN1,NPN1)
   REAL(fpk)  :: IG11(NPN1,NPN1),IG12(NPN1,NPN1)
   REAL(fpk)  :: IG21(NPN1,NPN1),IG22(NPN1,NPN1)

   REAL(fpk)  :: L_R11 (3,NPN1,NPN1), L_R12 (3,NPN1,NPN1)
   REAL(fpk)  :: L_R21 (3,NPN1,NPN1), L_R22 (3,NPN1,NPN1)
   REAL(fpk)  :: L_I11 (3,NPN1,NPN1), L_I12 (3,NPN1,NPN1)
   REAL(fpk)  :: L_I21 (3,NPN1,NPN1), L_I22 (3,NPN1,NPN1)

   REAL(fpk)  :: L_RG11(3,NPN1,NPN1), L_RG12(3,NPN1,NPN1)
   REAL(fpk)  :: L_RG21(3,NPN1,NPN1), L_RG22(3,NPN1,NPN1)
   REAL(fpk)  :: L_IG11(3,NPN1,NPN1), L_IG12(3,NPN1,NPN1)
   REAL(fpk)  :: L_IG21(3,NPN1,NPN1), L_IG22(3,NPN1,NPN1)

!  Tmatrix input ( Kind = 4 )

!   real(kind=4)  :: TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4)
!   real(kind=4)  :: TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4)
!   real(kind=4)  :: TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4)
!   real(kind=4)  :: TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)

!   real(kind=4)  :: L_TR11(3,NPN6,NPN4,NPN4),L_TR12(3,NPN6,NPN4,NPN4)
!   real(kind=4)  :: L_TR21(3,NPN6,NPN4,NPN4),L_TR22(3,NPN6,NPN4,NPN4)
!   real(kind=4)  :: L_TI11(3,NPN6,NPN4,NPN4),L_TI12(3,NPN6,NPN4,NPN4)
!   real(kind=4)  :: L_TI21(3,NPN6,NPN4,NPN4),L_TI22(3,NPN6,NPN4,NPN4)

   real(fpk)  :: TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4)
   real(fpk)  :: TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4)
   real(fpk)  :: TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4)
   real(fpk)  :: TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)

   real(fpk)  :: L_TR11(3,NPN6,NPN4,NPN4),L_TR12(3,NPN6,NPN4,NPN4)
   real(fpk)  :: L_TR21(3,NPN6,NPN4,NPN4),L_TR22(3,NPN6,NPN4,NPN4)
   real(fpk)  :: L_TI11(3,NPN6,NPN4,NPN4),L_TI12(3,NPN6,NPN4,NPN4)
   real(fpk)  :: L_TI21(3,NPN6,NPN4,NPN4),L_TI22(3,NPN6,NPN4,NPN4)

!  Coefficients

   REAL(fpk)  :: ALPH1(NPL),ALPH2(NPL),ALPH3(NPL)
   REAL(fpk)  :: ALPH4(NPL),BET1(NPL), BET2(NPL)

   REAL(fpk)  :: AL1(NPL),AL2(NPL),AL3(NPL)
   REAL(fpk)  :: AL4(NPL),BE1(NPL),BE2(NPL)
   REAL(fpk)  :: L_AL1(3,NPL),L_AL2(3,NPL),L_AL3(3,NPL)
   REAL(fpk)  :: L_AL4(3,NPL),L_BE1(3,NPL),L_BE2(3,NPL)

!  T-matrix results

   real(fpk)  :: TR1(NPN2,NPN2),TI1(NPN2,NPN2)
   REAL(fpk)  :: L_TR1(3,NPN2,NPN2), L_TI1(3,NPN2,NPN2)

!  Linearized coefficients

   REAL(fpk)  :: LPSD_ALPH1(NPL,3),LPSD_ALPH2(NPL,3),LPSD_ALPH3(NPL,3)
   REAL(fpk)  :: LPSD_ALPH4(NPL,3),LPSD_BET1(NPL,3), LPSD_BET2(NPL,3)

   REAL(fpk)  :: LRFE_ALPH1(NPL,3),LRFE_ALPH2(NPL,3),LRFE_ALPH3(NPL,3)
   REAL(fpk)  :: LRFE_ALPH4(NPL,3),LRFE_BET1(NPL,3), LRFE_BET2(NPL,3)

!  Local variables
!  ---------------

!  Debug flag and unit

   logical    :: do_debug_output
   integer    :: debug_unit, du

!  PSD functions

   REAL(fpk)  :: PSDF(2000), PSDF_DERIVS(2000,3)

!  Local integers, etc.

   integer       :: LMAX, L1MAX, L1M, L1, kpsd, LL
   integer       :: I, NK, INK, NCHECK, IXXX, II, INM1, M1, M
   integer       :: N, N1, N2, NN1, NN2, N11, N22, NM1, NNM, NM
   integer       :: NMA, NMAX, NMAX1, NMIN, MMAX
   integer       :: NGAUSS, NNNGGG, NGAUS, NGGG
   logical       :: first_loop, second_loop, third_loop
   logical       :: fail_1, fail_2, faildist, do_rfe_linearize
   character*90  :: message_1, message_2
   character*4   :: c4
   character*1   :: c1

!  Local FP Variables

   real(fpk)  :: Greek_pie, DDELT, RAT, Z1, Z2, Z3, WALB
   real(fpk)  :: PI, PPI, PIR, PII, radius, XEV, REFF, VEFF
   real(fpk)  :: NDENS, GXSEC, VOLUME
   real(fpk)  :: CSCAT, CEXT, COEFF1, CSCA, CEXTIN
   real(fpk)  :: DQSCA, DEXT, DSCA, WGSC, WGXT
   real(fpk)  :: QXT, QSC, QSCA1, QEXT1, QSCA, QEXT, DQEXT

   real(fpk)  :: LPSD_CSCAT(3), LPSD_CEXTIN(3)
   real(fpk)  :: LRFE_CSCAT(3), LRFE_CEXTIN(3)
   real(fpk)  :: L_CSCAT, L_Walb
   real(fpk)  :: Le_RAT, Le_radius, L_QSC(3), L_QXT(3)
   real(fpk)  :: L_QSCA(3), L_QEXT(3), L_QSCA1(3), L_QEXT1(3)
   real(fpk)  :: L_CSCA(3), L_CEXT(3)

   real(fpk)  :: TR1NN, TI1NN, TR1NN1, TI1NN1, DN1, WGI, WGII, WGIID
   real(fpk)  :: ZZ1, ZZ2, ZZ3, ZZ4, ZZ5, ZZ6, ZZ7, ZZ8
   real(fpk)  :: L_ZZ1, L_ZZ2, L_ZZ3, L_ZZ4, L_ZZ5, L_ZZ6, L_ZZ7, L_ZZ8
   real(fpk)  :: L_TR1NN, L_TI1NN, L_TR1NN1, L_TI1NN1

!  Report PSD calculation progress. Now linked to Tmat_Verbose flag

   logical    :: report_progress
   !logical    :: report_progress=.true.

!  Start Code
!  ----------

!  Report progress flag

   report_progress = Tmat_verbose

!  debug output flag

   do_debug_output = .false. ; du = 23
   debug_unit      = 23 ; du = debug_unit
   if ( do_debug_output ) then
      open(du,file='tmat_debug_output.log', status = 'unknown' )
   endif

!  Pie

   Greek_pie = acos(-d_one)

!  Initialize exception handling

   fail    = .false.
   istatus = 0
   message = ' '
   trace   = ' '
   trace_2 = ' '

!  Initialize Basic output

   tmat_dist      = 0.0d0
   LPSD_tmat_dist = 0.0d0

   tmat_bulk      = 0.0d0
   LPSD_tmat_bulk = 0.0d0

!  Initialize Coefficient output

   if ( Do_Expcoeffs ) then
      tmat_ncoeffs = 0
      tmat_expcoeffs = 0.0d0       ; tmat_asymm      = 0.0d0
      LPSD_tmat_expcoeffs = 0.0d0  ; LPSD_tmat_asymm = 0.0d0
   endif

!  Initialize F-matrix

   if ( Do_Expcoeffs.and. Do_Fmatrix ) then
      tmat_Fmatrix      = 0.0d0
      LPSD_tmat_Fmatrix = 0.0d0
   endif

!  Some elementary checks
!  ----------------------

!  Cannot PSD linearize with old-style PSDs !

   if ( Do_psd_OldStyle .and. Do_LinearPSD ) then
      message = 'Cannot perform PSD linearization with Old PSD method'
      trace   = 'tmat_matrix_PLUS module: Elementary check number 1 - abort'
      trace_2 = 'Action: change input flags'
      istatus = 2
      return
   endif

!  Check integer
!  -------------

!    NCHECK = 1, Except for Chebyshev Odd-power and spherical

   NCHECK=0
   IF (NP.EQ.-1.OR.NP.EQ.-2)      NCHECK=1   !  Spheroids, Cylinders
   IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK=1   !  Even-powered Chebyshevs

!  Equivalent-sphere Radius factor
!  -------------------------------

!    RAT = 1 for Volumes, Calculated RAT for Surfaces (3 choices)

!  Initial values (Good for the Volume representation)

   RAT = 1.0d0 ; Le_RAT = 0.0d0

!  For the ESAS representation, alculated RAT for Surfaces (3 choices)

   if ( Do_EqSaSphere ) then
      if ( NP.EQ.-1 ) then
         if ( Do_LinearEps ) then
            call  Eqv_radius_spheroids_plus ( eps, rat, Le_rat )
         else
            call Eqv_radius_spheroids ( eps, rat )
         endif
      else if ( NP.EQ.-2 ) then
         if ( Do_LinearEps ) then
            call  Eqv_radius_cylinder_plus ( eps, rat, Le_rat )
         else
            call Eqv_radius_cylinder ( eps, rat )
         endif
      else if ( NP.GT.0) then
         if ( Do_LinearEps ) then
            call  Eqv_radius_chebyshev_plus ( np, eps, rat, Le_rat )
         else
            call Eqv_radius_chebyshev ( np, eps, rat )
         endif
      else
         RAT    = 1.0d0
      endif
   endif

!  debug output

   if ( do_debug_output ) then
      write(du,8000)RAT 
      IF(NP.EQ.-1.AND.EPS.GE.1D0) write(du,7000)EPS
      IF(NP.EQ.-1.AND.EPS.LT.1D0) write(du,7001)EPS
      IF(NP.GE.0) write(du,7100)NP,EPS
      IF(NP.EQ.-2.AND.EPS.GE.1D0) write(du,7150)EPS
      IF(NP.EQ.-2.AND.EPS.LT.1D0) write(du,7151)EPS
      write(du,7400) LAMBDA,N_REAL,N_IMAG
      write(du,7200) ACCURACY
 8000 FORMAT ('RAT=',F8.6)
 7000 FORMAT('RANDOMLY ORIENTED OBLATE SPHEROIDS, A/B=',F11.7)
 7001 FORMAT('RANDOMLY ORIENTED PROLATE SPHEROIDS, A/B=',F11.7)
 7100 FORMAT('RANDOMLY ORIENTED CHEBYSHEV PARTICLES, T',I1,'(',F5.2,')')
 7150 FORMAT('RANDOMLY ORIENTED OBLATE CYLINDERS, D/L=',F11.7)
 7151 FORMAT('RANDOMLY ORIENTED PROLATE CYLINDERS, D/L=',F11.7)
 7200 FORMAT ('ACCURACY OF COMPUTATIONS DDELT = ',D8.2)
 7400 FORMAT('LAMBDA=',F10.6,3X,'N_REAL=',D10.4,3X,'N_IMAG=',D10.4)
   endif

!  Accuracy (Why 0.1d0 ???)

   DDELT = 0.1D0*ACCURACY

!  Various quantities

   PI  = 2.0d0 * Greek_pie / lambda
   PPI = PI * PI
   PIR = PPI * n_real
   PII = PPI * n_imag

!  Monodisperse
!  ------------

!    Multiply RAT and Le_RAT (maybe 0) by monodisperse radius.
!      Set trivial quadratures

   if ( Do_monodisperse ) then
      Do_psd_Oldstyle = .false.
      Do_LinearPSD    = .false.
      RAT    =    RAT * Monoradius
      Le_RAT = Le_RAT * Monoradius
      NKMAX = -1 ; NK = 1
      XG1(1) = 1.0d0 ; WG1(1) = 1.0d0
      WG1_derivs = 0.0d0
   endif

!  Skip PSD stuff if Monodisperse

   if ( Do_monodisperse ) GO TO 77

!  Particle Radii and PSD values
!  -----------------------------

!  Number of distribution points

   NK = NKMAX + 2

!  Power law, Fix the R1 and R2

   IF ( Do_psd_Oldstyle ) then
      IF (psd_Index.EQ.3 .and.FixR1R2 ) then
         CALL POWER (psd_pars(1),psd_pars(2),R1,R2)
      ENDIF
   ENDIF

!  Distribution quadratures

!   CALL GAULEG_wrong(-d_one,d_one,XG,WG,NK)
   CALL GAULEG_right(NK,0,0,XG,WG)

!  Normalize Quadratures to the given range [R1,R2]

   Z1=(R2-R1)*0.5D0
   Z2=(R1+R2)*0.5D0
   Z3=R1*0.5D0

!  For modified Power Law (NDISTR = 5 ) double the number of Quadratures

   IF ( Do_psd_Oldstyle ) then
      IF (psd_Index.EQ.5) THEN
         DO I=1,NK
            XG1(I)=Z3*XG(I)+Z3
            WG1(I)=WG(I)*Z3
         ENDDO
         DO I=NK+1,2*NK
            II=I-NK
            XG1(I)=Z1*XG(II)+Z2
            WG1(I)=WG(II)*Z1
         ENDDO
         NK=NK*2
      ELSE
         DO I=1,NK
            XG1(I)=Z1*XG(I)+Z2
            WG1(I)=WG(I)*Z1
         ENDDO
      ENDIF
   ENDIF

!  New Style, Set physical quantities

   IF ( .not. Do_psd_Oldstyle ) then
      DO I=1,NK
         XG1(I)=Z1*XG(I)+Z2
         WG1(I)=WG(I)*Z1
      ENDDO
   ENDIF

!  Debug
!   do i=1,nk
!      write(*,*)xg(i),wg(i), XG1(I), WG1(I)
!   enddo
!   pause

!  Distribution values. Old Scheme
!  ===============================

   IF ( Do_psd_OldStyle ) then

!  Psd old-style

      CALL DISTRB                                  &
       ( NK, XG1, psd_Index, R1, DO_DEBUG_OUTPUT,  & ! Inputs
         psd_pars(1), psd_pars(2), psd_pars(3),    & ! Inputs
         WG1, NDENS, GXSEC, VOLUME, REFF, VEFF)      ! Outputs

!  debug
!   do i = 1, nk
!      write(*,*) XG1(I), WG1(i)
!   enddo
!   write(*,*)reff,veff
!   pause

!  Distribution outputs

      tmat_dist(1) = NDENS
      tmat_dist(2) = GXSEC
      tmat_dist(3) = VOLUME
      tmat_dist(4) = REFF
      tmat_dist(5) = VEFF

   ENDIF

!  Distribution values, New scheme
!  ===============================

!    2000 is the dimension given for XG1 and WG1

   IF ( .not. Do_psd_OldStyle ) then

!  Linearized call
!  ---------------

      IF ( Do_LinearPSD ) then

!  Get the PSD and its derivatives

         CALL DISTRB_new_plus                    & 
       ( 2000, psd_Index, psd_pars, XG1, NK,     & ! Inputs
         PSDF, PSDF_derivs, message, faildist )    ! Outputs

!  Exception handling

         if ( faildist ) then
            write(c1,'(i1)')psd_index
            trace   = 'call from DISTRB_new_plus, PSD type '//c1
            trace_2 = 'tmat_matrix_PLUS module; distribution_Plus failure'
            fail = .true. ; istatus = 2 ; return
         endif

!  Bulk properties

         call DISTRB_bulk_plus                            &
          ( 2000, XG1, NK, psd_pars, PSDF, PSDF_derivs,   & ! Inputs
            WG1, WG1_derivs, tmat_dist, LPSD_tmat_dist )    ! Outputs

      endif

!
!      write(*,*)tmat_dist(1),LPSD_tmat_dist(1,1),LPSD_tmat_dist(1,2)
!      write(*,*)tmat_dist(2),LPSD_tmat_dist(2,1),LPSD_tmat_dist(2,2)
!      write(*,*)tmat_dist(3),LPSD_tmat_dist(3,1),LPSD_tmat_dist(3,2)
!      write(*,*)tmat_dist(4),LPSD_tmat_dist(4,1),LPSD_tmat_dist(4,2)
!      write(*,*)tmat_dist(5),LPSD_tmat_dist(5,1),LPSD_tmat_dist(5,2)
!      pause'kkkkkkkkkkkkkkkkkkkkk'

!  Call without linearization
!  --------------------------

      IF ( .not. Do_LinearPSD ) then

!  Get the PSD

         CALL DISTRB_new                     & 
       ( 2000, psd_Index, psd_pars, XG1, NK, & ! Inputs
         PSDF, message, faildist )             ! Outputs

!  Exception handling

         if ( faildist ) then
            write(c1,'(i1)')psd_index
            trace   = 'call from DISTRB_new, PSD type '//c1
            trace_2 = 'tmat_matrix_PLUS module; distribution failure'
            fail = .true. ; istatus = 2 ; return
         endif

!  Bulk properties

         call DISTRB_bulk        &
          ( 2000, XG1, NK, PSDF, & ! Inputs
            WG1, tmat_dist )       ! Outputs

      endif

!  Finish New distribution clause

   endif

!  debug. FD testing successful for Log-normal, 6 January 2011
!   kpsd = 1
!   write(*,*)tmat_dist(1),LPSD_tmat_dist(1,kpsd)
!   write(*,*)tmat_dist(2),LPSD_tmat_dist(2,kpsd)
!   write(*,*)tmat_dist(3),LPSD_tmat_dist(3,kpsd)
!   write(*,*)tmat_dist(4),LPSD_tmat_dist(4,kpsd)
!   write(*,*)tmat_dist(5),LPSD_tmat_dist(5,kpsd)
!   pause'hello linear'

!  Debug (adapted from the original)

   if ( do_debug_output ) then
      write(du,8002)R1,R2
 8002 FORMAT('R1=',F10.6,'   R2=',F10.6)
      IF (DABS(RAT-1D0).LE.1D-6) write(du,8003) tmat_dist(4) ,tmat_dist(5) 
      IF (DABS(RAT-1D0).GT.1D-6) write(du,8004) tmat_dist(4) ,tmat_dist(5) 
 8003 FORMAT('EQUAL-VOLUME-SPHERE REFF=',F8.4,'   VEFF=',F7.4)
 8004 FORMAT('EQUAL-SURFACE-AREA-SPHERE REFF=',F8.4, '   VEFF=',F7.4)
      write(du,7250)NK
 7250 FORMAT('NUMBER OF GAUSSIAN QUADRATURE POINTS ', 'IN SIZE AVERAGING =',I4)
   endif

!  ***********************
!       Main section
!  ***********************

!  Continuation point for Avoiding PSD

77 continue

!  Initialize
!  ----------

!  Initialize Bulk quantities

   CSCAT  = 0D0
   CEXTIN = 0D0

!  Initialize PSD Linearized Bulk quantities

   if ( .not. Do_monodisperse ) then
      if ( .not. Do_psd_Oldstyle .and. Do_LinearPSD ) then
         do kpsd = 1, 3
            LPSD_cscat(kpsd)  = 0.0d0
            LPSD_cextin(kpsd) = 0.0d0
         enddo
      endif
   endif

!  Initialize RFE linearization

   Do_RFE_linearize = Do_LinearRef .or. Do_LinearEps

   nlin = 0
   if ( Do_LinearRef ) nlin = 2
   if ( Do_LinearEps ) nlin = nlin + 1

   if ( Do_RFE_linearize ) then
      do LL = 1, nlin
         LRFE_cscat(LL)  = 0.0d0
         LRFE_cextin(LL) = 0.0d0
      enddo
   endif

!  Initialize Expansion Coefficients

   L1MAX  = 0
   if ( Do_Expcoeffs ) then
      DO I=1,NPL
         ALPH1(I)=0D0
         ALPH2(I)=0D0
         ALPH3(I)=0D0
         ALPH4(I)=0D0
         BET1(I)=0D0
         BET2(I)=0D0
      ENDDO      
   endif

!  Initialize PSD-Linearized Expansion coefficients

   if ( Do_Expcoeffs ) then
      if ( .not. Do_monodisperse ) then
         if ( .not. Do_psd_Oldstyle .and. Do_LinearPSD ) then
            do kpsd = 1, 3
               DO I=1,NPL
                  LPSD_ALPH1(I,kpsd)=0D0
                  LPSD_ALPH2(I,kpsd)=0D0
                  LPSD_ALPH3(I,kpsd)=0D0
                  LPSD_ALPH4(I,kpsd)=0D0
                  LPSD_BET1(I,kpsd)=0D0
                  LPSD_BET2(I,kpsd)=0D0
               ENDDO      
            enddo
         endif
      endif
   endif

!  Initialize REF-Linearized Expansion coefficients
!    PLACEHOLDER

   if ( Do_Expcoeffs ) then
      if ( Do_RFE_linearize ) then
         do LL = 1, 3
            DO I=1,NPL
               LRFE_ALPH1(I,LL)=0D0
               LRFE_ALPH2(I,LL)=0D0
               LRFE_ALPH3(I,LL)=0D0
               LRFE_ALPH4(I,LL)=0D0
               LRFE_BET1(I,LL)=0D0
               LRFE_BET2(I,LL)=0D0
            ENDDO      
         enddo
      endif
   endif

!  Start distribution point loop
!  -----------------------------

!   DO 56 INK=1,1         ! debug only

   DO 56 INK=1,NK

      I=NK-INK+1

!  Radius

      radius = RAT*XG1(I)
      if ( Do_LinearEps ) Le_radius = Le_RAT * XG1(I)

!  Check maximum particle size

      XEV  = 2D0 * Greek_pie * radius / lambda
      IXXX = XEV+4.05D0*XEV**0.333333D0
      INM1 = MAX0(4,IXXX)

      IF (INM1.GE.NPN1) then
         istatus = 2
         write(c4,'(I4)')inm1+2
         message ='INM1 >/= NPN1. Execution Terminated. Action: increase NPN1 to at least '//c4
         trace   = 'check maximum particle size'
         trace_2 = 'tmat_master_PLUS module; start of POint loop'
         fail = .true. ;return
      endif

!  initialize coefficients

      QEXT1=0D0
      QSCA1=0D0

      if ( Do_RFE_linearize ) then
         do LL = 1, nlin
            L_qext1(LL) = 0.0d0
            L_qsca1(LL) = 0.0d0
         enddo
      endif

!  Round ONE
!  =========

!  First Coefficient loop

      NMA = INM1 - 1
      FIRST_LOOP = .true.

      DO while ( first_loop .and. nma .lt. NPN1 )

         NMA = NMA + 1
         NMAX=NMA
         MMAX=1
         NGAUSS=NMAX*NDGS

!  Check NGAUSS value

         IF (NGAUSS.GE.NPNG1) then
            istatus = 2
            write(c4,'(I4)')ngauss+1
            message ='No convergence possible, local NPNG1 = '//c4//' exceeded'
            trace   = 'NGAUSS value too big: just before call to tmatrix_constants'
            trace_2 = 'tmat_master_PLUS module'
            fail    = .true. ; return
         endif

!  Set up constants

         if ( Do_LinearEps ) then
            call tmatrix_constants_plus                         &
              ( ngauss, nmax, np, eps,                          & ! inputs
                x, w, Le_x, Le_w, an, ann, s, ss, Le_s, Le_ss )   ! outputs
         else
            call tmatrix_constants             &
             ( ngauss, nmax, np, eps,          & ! inputs
               x, w, an, ann, s, ss )            ! outputs
         endif

!               if ( nma .eq. 20 ) then
!                write(*,*)nma,Radius, Le_Radius
!               do n = 1, ngauss*2
!                  write(51,*)n,s(n),ss(n),Le_s(n),Le_ss(n)
!               enddo
!               endif
!               pause'2'

!  Set up Bessel functions

         if ( Do_rfe_linearize ) then
            call tmatrix_vary_plus                               &
             ( np, ngauss, nmax, n_real, n_imag, eps,            & ! inputs
               Do_LinearEps, PI, radius, x, Le_radius, Le_x,     & ! inputs
               R, DR, DDR, DRR, DRI,                             & ! outputs
               j_bess,  y_bess,  jr_bess,  ji_bess,              & ! Outputs
               dj_bess, dy_bess, djr_bess, dji_bess,             &  ! Outputs
               Le_R, Le_DR, L_DDR, L_DRR, L_DRI, NLIN,           & ! outputs
               L_j_bess,  L_y_bess,  L_jr_bess,  L_ji_bess,      & ! Outputs
               L_dj_bess, L_dy_bess, L_djr_bess, L_dji_bess )      ! Outputs
         else
            call tmatrix_vary                          &
           ( n_real, n_imag, radius, PI, x,            & ! inputs
             eps, np, ngauss, nmax,                    & ! inputs
             R, DR, DDR, DRR, DRI,                     & ! outputs
             j_bess,  y_bess,  jr_bess,  ji_bess,      & ! Outputs
             dj_bess, dy_bess, djr_bess, dji_bess )      ! Outputs
         endif

!  Debug. FD tested, 7 January 2011
!      LL = 3
!      write(0,*)'First loop test',nmax,2*ngauss
!      do i = 1, 2*ngauss
!        write(77,116)i,ddr(i),drr(i),dri(i)
!        do n = 1, nmax
!           write(77,117)i,n,j_bess(i,n),jr_bess(i,n),ji_bess(i,n),y_bess(i,n)
!           write(77,117)i,n,dj_bess(i,n),djr_bess(i,n),dji_bess(i,n),dy_bess(i,n)
!         enddo
!      enddo
!      do i = 1, 2*ngauss
!        write(78,116)i,L_ddr(LL,i),L_drr(LL,i),L_dri(LL,i)
!        do n = 1, nmax
!           write(78,117)i,n,L_j_bess(LL,i,n),L_jr_bess(LL,i,n),L_ji_bess(LL,i,n),L_y_bess(LL,i,n)
!           write(78,117)i,n,L_dj_bess(LL,i,n),L_djr_bess(LL,i,n),L_dji_bess(LL,i,n),L_dy_bess(LL,i,n)
!         enddo
!      enddo
!      pause'First loop test, Write fort 77 and 78'
!116   format(i5 ,1p3e20.10)
!117   format(2i5,1p4e20.10)

!  First Call to TMATRIX-R0
!  ------------------------

         if ( do_rfe_linearize ) then

            call tmatrix_R0_plus                                 &
             ( Do_LinearEps, NGAUSS, NMAX, NCHECK, NLIN,         & ! Inputs
               X, W, Le_X, Le_W, AN, ANN, PPI, PIR, PII,         & ! Inputs
               R, DR, DDR, DRR, DRI,                             & ! Inputs
                j_bess,  y_bess,  jr_bess,  ji_bess,             & ! Inputs
               dj_bess, dy_bess, djr_bess, dji_bess,             & ! Inputs
               Le_R, Le_DR, L_DDR, L_DRR, L_DRI,                 & ! Inputs
               L_j_bess,  L_y_bess,  L_jr_bess,  L_ji_bess,      & ! Inputs
               L_dj_bess, L_dy_bess, L_djr_bess, L_dji_bess,     & ! Inputs
               R12, R21, I12, I21, RG12, RG21, IG12, IG21,       & ! Outputs
               L_R12, L_R21, L_I12, L_I21,                       & ! Outputs
               L_RG12, L_RG21, L_IG12, L_IG21,                   & ! Outputs
               TR1, TI1, L_TR1, L_TI1, fail, message, trace )      ! Outputs

            if ( fail ) then
               istatus = 2
               write(c4,'(i4)')nma
               trace_2 = 'tmat_master_PLUS module: First Call to Tmatrix_r0_plus, NMA = '//c4
               return
            endif

         else

            call tmatrix_R0                                &
             ( NGAUSS, NMAX, NCHECK, X, W, AN, ANN,        & ! Inputs
               PPI, PIR, PII, R, DR, DDR, DRR, DRI,        & ! Inputs
                j_bess,  y_bess,  jr_bess,  ji_bess,       & ! Inputs
               dj_bess, dy_bess, djr_bess, dji_bess,       & ! Inputs
               R12, R21, I12, I21, RG12, RG21, IG12, IG21, & ! Outputs
               TR1, TI1, fail, message, trace )              ! Outputs

            if ( fail ) then
               istatus = 2
               write(c4,'(i4)')nma
               trace_2 = 'tmat_master_PLUS module: First Call to Tmatrix_r0, NMA = '//c4
               return
            endif

         endif

!  Local computations, and absolute differences
!  --------------------------------------------

!  Initialize

         QEXT=0D0
         QSCA=0D0
         if ( Do_RFE_linearize ) then
            do LL = 1, nlin
               L_qext(LL) = 0.0d0
               L_qsca(LL) = 0.0d0
            enddo
         endif

!  Summation loop

         DO N=1,NMAX

            N1=N+NMAX
            !DN1=DFLOAT(2*N+1)
            DN1=REAL(2*N+1,fpk)

            TR1NN  = TR1(N,N)    ;  TI1NN  = TI1(N,N)
            TR1NN1 = TR1(N1,N1)  ;  TI1NN1 = TI1(N1,N1)
            QEXT = QEXT + DN1 * ( TR1NN+TR1NN1 )
            QSCA = QSCA + DN1 * ( TR1NN*TR1NN   + TI1NN*TI1NN     &
                                + TR1NN1*TR1NN1 + TI1NN1*TI1NN1 )

            if ( do_rfe_linearize ) then
               do LL = 1, nlin
                  L_TR1NN  = L_TR1(LL,N,N)    ;  L_TI1NN  = L_TI1(LL,N,N)
                  L_TR1NN1 = L_TR1(LL,N1,N1)  ;  L_TI1NN1 = L_TI1(LL,N1,N1)
                  L_QEXT(LL) = L_QEXT(LL) + DN1 * ( L_TR1NN+L_TR1NN1 )
                  L_QSCA(LL) = L_QSCA(LL) + DN1 * &
                         ( L_TR1NN*TR1NN   + TR1NN*L_TR1NN     &
                         + L_TI1NN*TI1NN   + TI1NN*L_TI1NN     &
                         + L_TR1NN1*TR1NN1 + TR1NN1*L_TR1NN1   &
                         + L_TI1NN1*TI1NN1 + TI1NN1*L_TI1NN1 )
               enddo
            endif

         ENDDO

!  Difference

         DSCA=DABS((QSCA1-QSCA)/QSCA)
         DEXT=DABS((QEXT1-QEXT)/QEXT)

!  Upgrade Local values

         QEXT1=QEXT
         QSCA1=QSCA
         if ( Do_RFE_linearize ) then
            do LL = 1, nlin
               L_qext1(LL) = L_qext(LL)
               L_qsca1(LL) = L_qsca(LL)
            enddo
         endif

!  Second loop convergence, to establish how many terms

         NMIN=DBLE(NMAX)/2D0+1D0
         second_loop = .true.
         N = Nmin - 1

         DO while ( second_loop .and. N.lt.nmax )

            N = N + 1
            N1=N+NMAX
            !DN1=DFLOAT(2*N+1)
            DN1=REAL(2*N+1,fpk)

            TR1NN  = TR1(N,N)    ;  TI1NN  = TI1(N,N)
            TR1NN1 = TR1(N1,N1)  ;  TI1NN1 = TI1(N1,N1)
            DQEXT = DN1 * ( TR1NN+TR1NN1 )
            DQSCA = DN1 * ( TR1NN*TR1NN   + TI1NN*TI1NN     &
                          + TR1NN1*TR1NN1 + TI1NN1*TI1NN1 )

            DQSCA=DABS(DQSCA/QSCA)
            DQEXT=DABS(DQEXT/QEXT)
            NMAX1=N
            IF (DQSCA.LE.DDELT.AND.DQEXT.LE.DDELT) second_loop = .false.

         ENDDO

!         write(0,*)'After 12',DQSCA,DQEXT
!         pause

!  Exception handling

         IF (NMA.EQ.NPN1) THEN
            fail = .true.
            message = 'No convergence, NPN1=   . EXECUTION TERMINATED'
            trace   = 'First use of Tmatrix for DQSCA and DQEXT' 
            trace_2 = 'tmat_master_PLUS module'
            istatus = 2
            return
         endif

!  First loop covergence

         IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) first_loop = .false.

      ENDDO

!  Warning flag

      IF (NGAUSS.EQ.NPNG1) THEN
         message = 'WARNING: NGAUSS=NPNG1'
         trace = 'After First use of Tmatrix for DQSCA and DQEXT' 
         trace_2 = 'tmat_master_PLUS module'
         istatus = 1
      endif

!  Round Two
!  =========

!  third Coefficient loop

      NNNGGG=NGAUSS+1
      NGAUS = NNNGGG - 1
      THIRD_LOOP = .true.
      MMAX=NMAX1

      DO while ( third_loop .and. NGAUS .lt. NPNG1 )

         NGAUS = NGAUS + 1
         NGAUSS=NGAUS
         NGGG=2*NGAUSS

!  Set up constants

         if ( Do_LinearEps ) then
            call tmatrix_constants_plus             &
              ( ngauss, nmax, np, eps,                          & ! inputs
                x, w, Le_x, Le_w, an, ann, s, ss, Le_s, Le_ss )   ! outputs
         else
            call tmatrix_constants             &
             ( ngauss, nmax, np, eps,          & ! inputs
               x, w, an, ann, s, ss )            ! outputs
         endif

!  Set up Bessel functions

         if ( Do_rfe_linearize ) then
            call tmatrix_vary_plus                               &
             ( np, ngauss, nmax, n_real, n_imag, eps,            & ! inputs
               Do_LinearEps, PI, radius, x, Le_radius, Le_x,     & ! inputs
               R, DR, DDR, DRR, DRI,                             & ! outputs
               j_bess,  y_bess,  jr_bess,  ji_bess,              & ! Outputs
               dj_bess, dy_bess, djr_bess, dji_bess,             &  ! Outputs
               Le_R, Le_DR, L_DDR, L_DRR, L_DRI, NLIN,           & ! outputs
               L_j_bess,  L_y_bess,  L_jr_bess,  L_ji_bess,      & ! Outputs
               L_dj_bess, L_dy_bess, L_djr_bess, L_dji_bess )      ! Outputs
         else
            call tmatrix_vary                          &
           ( n_real, n_imag, radius, PI, x,            & ! inputs
             eps, np, ngauss, nmax,                    & ! inputs
             R, DR, DDR, DRR, DRI,                     & ! outputs
             j_bess,  y_bess,  jr_bess,  ji_bess,      & ! Outputs
             dj_bess, dy_bess, djr_bess, dji_bess )      ! Outputs
         endif

!  Second Call to TMATRIX-R0
!  -------------------------

         if ( do_rfe_linearize ) then

            call tmatrix_R0_plus                                 &
             ( Do_LinearEps, NGAUSS, NMAX, NCHECK, NLIN,         & ! Inputs
               X, W, Le_X, Le_W, AN, ANN, PPI, PIR, PII,         & ! Inputs
               R, DR, DDR, DRR, DRI,                             & ! Inputs
                j_bess,  y_bess,  jr_bess,  ji_bess,             & ! Inputs
               dj_bess, dy_bess, djr_bess, dji_bess,             & ! Inputs
               Le_R, Le_DR, L_DDR, L_DRR, L_DRI,                 & ! Inputs
               L_j_bess,  L_y_bess,  L_jr_bess,  L_ji_bess,      & ! Inputs
               L_dj_bess, L_dy_bess, L_djr_bess, L_dji_bess,     & ! Inputs
               R12, R21, I12, I21, RG12, RG21, IG12, IG21,       & ! Outputs
               L_R12, L_R21, L_I12, L_I21,                       & ! Outputs
               L_RG12, L_RG21, L_IG12, L_IG21,                   & ! Outputs
               TR1, TI1, L_TR1, L_TI1, fail, message, trace )      ! Outputs

            if ( fail ) then
               istatus = 2
               write(c4,'(i4)')nma
               trace_2 = 'tmat_master_PLUS module: Second Call to Tmatrix_r0_plus, NMA = '//c4
               return
            endif

         else

            call tmatrix_R0                                &
             ( NGAUSS, NMAX, NCHECK, X, W, AN, ANN,        & ! Inputs
               PPI, PIR, PII, R, DR, DDR, DRR, DRI,        & ! Inputs
                j_bess,  y_bess,  jr_bess,  ji_bess,       & ! Inputs
               dj_bess, dy_bess, djr_bess, dji_bess,       & ! Inputs
               R12, R21, I12, I21, RG12, RG21, IG12, IG21, & ! Outputs
               TR1, TI1, fail, message, trace )              ! Outputs

            if ( fail ) then
               istatus = 2
               write(c4,'(i4)')nma
               trace_2 = 'tmat_master_PLUS module : Second Call to Tmatrix_r0, NMA = '//c4
               return
            endif

         endif

!  calculate again
!  ---------------

!  Initialize

         QEXT=0D0
         QSCA=0D0
         if ( Do_RFE_linearize ) then
            do LL = 1, nlin
               L_qext(LL) = 0.0d0
               L_qsca(LL) = 0.0d0
            enddo
         endif

!  Summation loop

         DO N=1,NMAX

            N1=N+NMAX
            !DN1=DFLOAT(2*N+1)
            DN1=REAL(2*N+1,fpk)

            TR1NN  = TR1(N,N)    ;  TI1NN  = TI1(N,N)
            TR1NN1 = TR1(N1,N1)  ;  TI1NN1 = TI1(N1,N1)
            QEXT = QEXT + DN1 * ( TR1NN+TR1NN1 )
            QSCA = QSCA + DN1 * ( TR1NN*TR1NN   + TI1NN*TI1NN     &
                                + TR1NN1*TR1NN1 + TI1NN1*TI1NN1 )

            if ( do_rfe_linearize ) then
               do LL = 1, nlin
                  L_TR1NN  = L_TR1(LL,N,N)    ;  L_TI1NN  = L_TI1(LL,N,N)
                  L_TR1NN1 = L_TR1(LL,N1,N1)  ;  L_TI1NN1 = L_TI1(LL,N1,N1)
                  L_QEXT(LL) = L_QEXT(LL) + DN1 * ( L_TR1NN+L_TR1NN1 )
                  L_QSCA(LL) = L_QSCA(LL) + DN1 * &
                         ( L_TR1NN*TR1NN   + TR1NN*L_TR1NN     &
                         + L_TI1NN*TI1NN   + TI1NN*L_TI1NN     &
                         + L_TR1NN1*TR1NN1 + TR1NN1*L_TR1NN1   &
                         + L_TI1NN1*TI1NN1 + TI1NN1*L_TI1NN1 )
               enddo
            endif

         ENDDO

!  Difference

         DSCA=DABS((QSCA1-QSCA)/QSCA)
         DEXT=DABS((QEXT1-QEXT)/QEXT)

!  Upgrade Local values

         QEXT1=QEXT
         QSCA1=QSCA
         if ( Do_RFE_linearize ) then
            do LL = 1, nlin
               L_qext1(LL) = L_qext(LL)
               L_qsca1(LL) = L_qsca(LL)
            enddo
         endif

!  Convergence on third loop

         IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) third_loop = .false.

!  Warning flag
  
         IF (NGAUS.EQ.NPNG1) THEN
            message = 'WARNING: NGAUS=NPNG1'
            trace   = 'Occuring at the end of third Tmatrix loop'
            trace_2 = 'tmat_master_PLUS module'
            istatus = 1
         endif

!  End 3rd loop

      enddo

!  Exception handling

      IF (NMAX1.GT.NPN4) THEN
         write(c4,'(i4)')nmax1+1
         message ='NMAX1 > NPN4. No convergence, execution terminated. Action: increase NPN4 to at least '//c4
         trace   = 'After second application of TMAT'
         trace_2 = 'tmat_master_PLUS module'
         fail = .true.
         istatus = 2
         return
      endif

!  CALCULATION PROPER (Results for PHI-angle M = 0 Fourier Component )
!  ==================

!  Initialize

      QEXT=0D0
      QSCA=0D0
      if ( Do_RFE_linearize ) then
         do LL = 1, nlin
            L_qext(LL) = 0.0d0
            L_qsca(LL) = 0.0d0
         enddo
      endif

!  Again for Qext and its linearization

      NNM=NMAX*2
      DO N=1,NNM
         QEXT = QEXT + TR1(N,N)
      ENDDO

      if ( Do_RFE_linearize ) then
         do LL = 1, nlin
            DO N=1,NNM
               L_qext(LL) = L_qext(LL) + L_TR1(LL,N,N)
            ENDDO
         enddo
      endif

!  Again for QSca and its linearization

      DO N2=1,NMAX1
         NN2=N2+NMAX
         DO N1=1,NMAX1
            NN1=N1+NMAX

            ZZ1=TR1(N1,N2)      ;  TR11(1,N1,N2)=ZZ1
            ZZ2=TI1(N1,N2)      ;  TI11(1,N1,N2)=ZZ2
            ZZ3=TR1(N1,NN2)     ;  TR12(1,N1,N2)=ZZ3
            ZZ4=TI1(N1,NN2)     ;  TI12(1,N1,N2)=ZZ4
            ZZ5=TR1(NN1,N2)     ;  TR21(1,N1,N2)=ZZ5
            ZZ6=TI1(NN1,N2)     ;  TI21(1,N1,N2)=ZZ6
            ZZ7=TR1(NN1,NN2)    ;  TR22(1,N1,N2)=ZZ7
            ZZ8=TI1(NN1,NN2)    ;  TI22(1,N1,N2)=ZZ8
            QSCA=QSCA+ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4 + &
                      ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8

            if ( Do_RFE_linearize ) then
               DO LL = 1, nlin
                  L_ZZ1=L_TR1(LL,N1,N2)      ;  L_TR11(LL,1,N1,N2)=L_ZZ1
                  L_ZZ2=L_TI1(LL,N1,N2)      ;  L_TI11(LL,1,N1,N2)=L_ZZ2
                  L_ZZ3=L_TR1(LL,N1,NN2)     ;  L_TR12(LL,1,N1,N2)=L_ZZ3
                  L_ZZ4=L_TI1(LL,N1,NN2)     ;  L_TI12(LL,1,N1,N2)=L_ZZ4
                  L_ZZ5=L_TR1(LL,NN1,N2)     ;  L_TR21(LL,1,N1,N2)=L_ZZ5
                  L_ZZ6=L_TI1(LL,NN1,N2)     ;  L_TI21(LL,1,N1,N2)=L_ZZ6
                  L_ZZ7=L_TR1(LL,NN1,NN2)    ;  L_TR22(LL,1,N1,N2)=L_ZZ7
                  L_ZZ8=L_TI1(LL,NN1,NN2)    ;  L_TI22(LL,1,N1,N2)=L_ZZ8
                  L_QSCA(LL) = L_QSCA(LL) + 2.0d0 * &
                     ( L_ZZ1*ZZ1+L_ZZ2*ZZ2+L_ZZ3*ZZ3+L_ZZ4*ZZ4 + &
                       L_ZZ5*ZZ5+L_ZZ6*ZZ6+L_ZZ7*ZZ7+L_ZZ8*ZZ8 )
               ENDDO
            ENDIF

         ENDDO
      ENDDO

!  ROUND THREE (PHI-angle Fourier components M > 0 )
!  ===========

!  Start Fourier loop

      DO 220 M=1,NMAX1

!  TMAT_R call

         if ( do_rfe_linearize ) then

            call tmatrix_R_plus                                        &
            ( Do_LinearEps, M, NGAUSS, NMAX, NCHECK, NLIN,             & ! Inputs
              X, W, Le_X, Le_W, AN, ANN, S, SS, Le_S, Le_SS,           & ! Inputs
              PPI, PIR, PII, R, DR, DDR, DRR, DRI,                     & ! Inputs
               j_bess,  y_bess,  jr_bess,  ji_bess,                    & ! Inputs
              dj_bess, dy_bess, djr_bess, dji_bess,                    & ! Inputs
              Le_R, Le_DR, L_DDR, L_DRR, L_DRI,                        & ! Inputs
              L_j_bess,  L_y_bess,  L_jr_bess,  L_ji_bess,             & ! Inputs
              L_dj_bess, L_dy_bess, L_djr_bess, L_dji_bess,            & ! Inputs
              R11, R12, R21, R22, I11, I12, I21, I22,                  & ! Outputs
              RG11,RG12,RG21,RG22,IG11,IG12,IG21,IG22,                 & ! Outputs
              L_R11, L_R12, L_R21, L_R22, L_I11, L_I12, L_I21, L_I22,  & ! Outputs
              L_RG11,L_RG12,L_RG21,L_RG22,L_IG11,L_IG12,L_IG21,L_IG22, & ! Outputs
              TR1, TI1, L_TR1, L_TI1, fail, message, trace )             ! Outputs

            if ( fail ) then
               write(c4,'(i4)')m
               trace_2 = 'tmat_master_PLUS module: Third Call to Tmatrix_R_plus, M = '//c4
               istatus = 2
               return
            endif

         else
          
            call tmatrix_R                                   &
            ( M, NGAUSS, NMAX, NCHECK, X, W, AN, ANN, S, SS, & ! Inputs
              PPI, PIR, PII, R, DR, DDR, DRR, DRI,           & ! Inputs
               j_bess,  y_bess,  jr_bess,  ji_bess,          & ! Inputs
              dj_bess, dy_bess, djr_bess, dji_bess,          & ! Inputs
              R11, R12, R21, R22, I11, I12, I21, I22,        & ! Outputs
              RG11,RG12,RG21,RG22,IG11,IG12,IG21,IG22,       & ! Outputs
              TR1, TI1, fail, message, trace )                 ! Outputs

            if ( fail ) then
               write(c4,'(i4)')m
               trace_2 = 'tmat_master_PLUS module: Third Call to Tmatrix_R, M = '//c4
               istatus = 2
               return
            endif

         endif

!  Count

         NM=NMAX-M+1
         NM1=NMAX1-M+1
         M1=M+1

!  Initialize

         QXT=0D0
         QSC=0D0
         if ( Do_RFE_linearize ) then
            do LL = 1, nlin
               L_qxt(LL) = 0.0d0
               L_qsc(LL) = 0.0d0
            enddo
         endif

!  QExt calculation for this component

         NNM=2*NM
         DO N=1,NNM
            QXT = QXT + TR1(N,N) * 2.0d0
         ENDDO

         if ( Do_RFE_linearize ) then
            do LL = 1, nlin
               DO N=1,NNM
                  L_qxt(LL) = L_qxt(LL) + L_TR1(LL,N,N) * 2.0d0
               ENDDO
            enddo
         endif

!  QSCa calculation for this component

         DO N2=1,NM1
            NN2=N2+M-1
            N22=N2+NM
            DO N1=1,NM1
               NN1=N1+M-1
               N11=N1+NM

               ZZ1=TR1(N1,N2)    ;   TR11(M1,NN1,NN2)=ZZ1
               ZZ2=TI1(N1,N2)    ;   TI11(M1,NN1,NN2)=ZZ2
               ZZ3=TR1(N1,N22)   ;   TR12(M1,NN1,NN2)=ZZ3
               ZZ4=TI1(N1,N22)   ;   TI12(M1,NN1,NN2)=ZZ4
               ZZ5=TR1(N11,N2)   ;   TR21(M1,NN1,NN2)=ZZ5
               ZZ6=TI1(N11,N2)   ;   TI21(M1,NN1,NN2)=ZZ6
               ZZ7=TR1(N11,N22)  ;   TR22(M1,NN1,NN2)=ZZ7
               ZZ8=TI1(N11,N22)  ;   TI22(M1,NN1,NN2)=ZZ8
               QSC=QSC + ( ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4 &
                       +   ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8 ) * 2D0

               if ( Do_RFE_linearize ) then
                  DO LL = 1, nlin
                     L_ZZ1=L_TR1(LL,N1,N2)      ;  L_TR11(LL,M1,NN1,NN2)=L_ZZ1
                     L_ZZ2=L_TI1(LL,N1,N2)      ;  L_TI11(LL,M1,NN1,NN2)=L_ZZ2
                     L_ZZ3=L_TR1(LL,N1,N22)     ;  L_TR12(LL,M1,NN1,NN2)=L_ZZ3
                     L_ZZ4=L_TI1(LL,N1,N22)     ;  L_TI12(LL,M1,NN1,NN2)=L_ZZ4
                     L_ZZ5=L_TR1(LL,N11,N2)     ;  L_TR21(LL,M1,NN1,NN2)=L_ZZ5
                     L_ZZ6=L_TI1(LL,N11,N2)     ;  L_TI21(LL,M1,NN1,NN2)=L_ZZ6
                     L_ZZ7=L_TR1(LL,N11,N22)    ;  L_TR22(LL,M1,NN1,NN2)=L_ZZ7
                     L_ZZ8=L_TI1(LL,N11,N22)    ;  L_TI22(LL,M1,NN1,NN2)=L_ZZ8
                     L_QSC(LL) = L_QSC(LL) + 4.0d0 * &
                           ( L_ZZ1*ZZ1+L_ZZ2*ZZ2+L_ZZ3*ZZ3+L_ZZ4*ZZ4 + &
                             L_ZZ5*ZZ5+L_ZZ6*ZZ6+L_ZZ7*ZZ7+L_ZZ8*ZZ8 )
                  ENDDO
               ENDIF

            ENDDO
         ENDDO

!  Add Fourier component to total

         QEXT = QEXT + QXT
         QSCA = QSCA + QSC
         if ( Do_RFE_linearize ) then
            do LL = 1, nlin
               L_qext(LL) = L_qext(LL) + L_qxt(LL)
               L_qsca(LL) = L_qsca(LL) + L_qsc(LL)
            enddo
         endif

!  End fourier loop

  220 CONTINUE

!  Final section
!  =============

!  multiply coefficients by Lam^2/2pi

      COEFF1 = lambda * lambda * 0.5D0 / Greek_pie

!  Local Bulk values

      CSCA =  QSCA*COEFF1 
      CEXT = -QEXT*COEFF1
      if ( Do_RFE_linearize ) then
         do LL = 1, nlin
            L_CSCA(LL) = + L_QSCA(LL) *COEFF1
            L_CEXT(LL) = - L_QEXT(LL) *COEFF1
         enddo
      endif

!  Local Coefficients
!  ------------------

      if ( Do_Expcoeffs ) then

!  Main routine for linearized + regular coefficients

         if ( do_rfe_linearize ) then

            call GSP_plus                                 &
             ( NMAX1, NLIN, LAMBDA, CSCA, L_CSCA,         & ! Inputs
               TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22,   & ! Inputs
               L_TR11,L_TR12,L_TR21,L_TR22,               & ! Inputs
               L_TI11,L_TI12,L_TI21,L_TI22,               & ! Inputs
               AL1,AL2,AL3,AL4,BE1,BE2,LMAX,              & ! Outputs
               L_AL1,L_AL2,L_AL3,L_AL4,L_BE1,L_BE2,       & ! Outputs
               fail, message, trace )                       ! Outputs

            if ( fail ) then
               trace_2 = 'tmat_master_PLUS module: Call to GSP_plus Coefficients routine'
               istatus = 2
               return
            endif

         else

            CALL GSP ( NMAX1, LAMBDA, CSCA,              & ! Inputs
              TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22,   & ! Inputs
              AL1,AL2,AL3,AL4,BE1,BE2,LMAX,              & ! Outputs
              fail, message, trace )                       ! Outputs

            if ( fail ) then
               trace_2 = 'tmat_master_PLUS module: Call to GSP Coefficients routine'
               istatus = 2
               return
            endif

         endif

!         pause'1'

!  Local maximum count

         L1M=LMAX+1
         L1MAX=MAX(L1MAX,L1M)

      endif

!  Polydisperse Summed values (Trivial, if Monodisperse)
!  ==========================

!  Bulk quantities.  WGII = 1.0 for MONO

      WGII = WG1(I)
      WGXT = WGII * CEXT
      WGSC = WGII * CSCA
      CSCAT  = CSCAT  + WGSC
      CEXTIN = CEXTIN + WGXT

!  PSD-linearized Bulk values

      if ( .not. do_monodisperse ) then
         if ( .not. Do_psd_OldStyle .and. Do_LinearPSD ) then
            do kpsd = 1, 3
               WGIID = WG1_derivs(I,kpsd)
               LPSD_cscat(kpsd)  = LPSD_cscat(kpsd)  + CSCA * WGIID
               LPSD_cextin(kpsd) = LPSD_cextin(kpsd) + CEXT * WGIID
            enddo
         endif
      endif

!  RFE-linearized Bulk values

      if ( Do_rfe_Linearize ) then
         do LL = 1, nlin
            LRFE_cscat(LL)  = LRFE_cscat (LL) + L_CSCA(LL) * WGII
            LRFE_cextin(LL) = LRFE_cextin(LL) + L_CEXT(LL) * WGII
         enddo
      endif

!  debug

!   write(*,*)cscat, cextin
!   write(*,*)LRFE_cscat(1), LRFE_cextin(1)
!   pause'res1'

!  Expansion coefficients
!  ----------------------

!  Regular

      if ( Do_Expcoeffs ) then
         DO 250 L1=1,L1M
            ALPH1(L1)=ALPH1(L1)+AL1(L1)*WGSC
            ALPH2(L1)=ALPH2(L1)+AL2(L1)*WGSC
            ALPH3(L1)=ALPH3(L1)+AL3(L1)*WGSC
            ALPH4(L1)=ALPH4(L1)+AL4(L1)*WGSC
            BET1(L1)=BET1(L1)+BE1(L1)*WGSC
            BET2(L1)=BET2(L1)+BE2(L1)*WGSC
  250    CONTINUE
      endif

!  PSD-linearized expansion coefficients

      if ( Do_Expcoeffs ) then
         if ( .not. do_monodisperse ) then
            if ( .not. Do_psd_OldStyle .and. Do_LinearPSD ) then
               DO 251 L1=1,L1M
                  do kpsd = 1, 3
                     WGSC = CSCA * WG1_derivs(I,kpsd)
                     LPSD_ALPH1(L1,kpsd)= LPSD_ALPH1(L1,kpsd)+AL1(L1)*WGSC
                     LPSD_ALPH2(L1,kpsd)= LPSD_ALPH2(L1,kpsd)+AL2(L1)*WGSC
                     LPSD_ALPH3(L1,kpsd)= LPSD_ALPH3(L1,kpsd)+AL3(L1)*WGSC
                     LPSD_ALPH4(L1,kpsd)= LPSD_ALPH4(L1,kpsd)+AL4(L1)*WGSC
                     LPSD_BET1(L1,kpsd) = LPSD_BET1(L1,kpsd) +BE1(L1)*WGSC
                     LPSD_BET2(L1,kpsd) = LPSD_BET2(L1,kpsd) +BE2(L1)*WGSC
                  enddo
  251          CONTINUE
            endif
         endif
      endif

!  RFE-linearized expansion coefficients

      if ( Do_Expcoeffs ) then
         if ( Do_rfe_Linearize ) then
            DO 259 L1 = 1, L1M
               do  LL = 1, nlin
                  WGI = WG1(I)
                  LRFE_ALPH1(L1,LL)= LRFE_ALPH1(L1,LL) + (AL1(L1) * L_CSCA(LL) + L_AL1(LL,L1) * CSCA ) * WGI
                  LRFE_ALPH2(L1,LL)= LRFE_ALPH2(L1,LL) + (AL2(L1) * L_CSCA(LL) + L_AL2(LL,L1) * CSCA ) * WGI
                  LRFE_ALPH3(L1,LL)= LRFE_ALPH3(L1,LL) + (AL3(L1) * L_CSCA(LL) + L_AL3(LL,L1) * CSCA ) * WGI
                  LRFE_ALPH4(L1,LL)= LRFE_ALPH4(L1,LL) + (AL4(L1) * L_CSCA(LL) + L_AL4(LL,L1) * CSCA ) * WGI
                  LRFE_BET1(L1,LL) = LRFE_BET1(L1,LL)  + (BE1(L1) * L_CSCA(LL) + L_BE1(LL,L1) * CSCA ) * WGI
                  LRFE_BET2(L1,LL) = LRFE_BET2(L1,LL)  + (BE2(L1) * L_CSCA(LL) + L_BE2(LL,L1) * CSCA ) * WGI
               enddo
  259       CONTINUE
         endif
      endif

!  progress

      if (report_progress) then
         if ( do_monodisperse ) then
            write(*,'(a,i5)')'     -- Done Monodisperse, NMAX = ',NMAX
         else
            write(*,'(a,i4,a,i5)')'     -- Done PSD point # ', INK,', NMAX = ',NMAX
         endif
      endif

!  debug
!      if (INK.eq.1)pause

!  End PSD loop

56 CONTINUE

!  Output generation for Bulk values
!  =================================

!  Regular

   WALB=CSCAT/CEXTIN
   tmat_bulk(1) = CEXTIN
   tmat_bulk(2) = CSCAT
   tmat_bulk(3) = WALB

!  Warning on WALB

   if ( WALB  > 1.0d0 ) then
      fail = .true.
      istatus = 1
      message = 'WARNING: W IS GREATER THAN 1'
      trace   = 'Output section (bulk)'
      trace_2 = 'tmat_master_PLUS module' 
   endif

!  PSD-linearized outputs

   if ( .not. do_monodisperse ) then
      if ( .not. Do_psd_OldStyle .and. Do_LinearPSD ) then
         do kpsd = 1, 3
            L_Walb  = Walb * &
             ( ( LPSD_cscat(kpsd) / CSCAT ) - ( LPSD_cextin(kpsd) / CEXTIN ) )
            LPSD_tmat_bulk(1,kpsd)  = LPSD_cextin(kpsd)
            LPSD_tmat_bulk(2,kpsd)  = LPSD_cscat(kpsd)
            LPSD_tmat_bulk(3,kpsd)  = L_Walb
         enddo
      endif
   endif

!  RFE-linearized outputs

   if ( Do_rfe_Linearize ) then
      do LL = 1, nlin
         L_Walb  = Walb * &
          ( ( LRFE_cscat(LL) / CSCAT ) - ( LRFE_cextin(LL) / CEXTIN ) )
         LRFE_tmat_bulk(1,LL)  = LRFE_cextin(LL)
         LRFE_tmat_bulk(2,LL)  = LRFE_cscat(LL)
         LRFE_tmat_bulk(3,LL)  = L_Walb
      enddo
   endif

!  Output generation for Expansion coefficients
!  --------------------------------------------

!  Only if flag set  for Expansion coefficients output

   if ( Do_Expcoeffs ) then

!  Normalize output

      DO 510 L1=1,L1MAX
         ALPH1(L1)=ALPH1(L1)/CSCAT
         ALPH2(L1)=ALPH2(L1)/CSCAT
         ALPH3(L1)=ALPH3(L1)/CSCAT
         ALPH4(L1)=ALPH4(L1)/CSCAT
         BET1(L1)=BET1(L1)/CSCAT
         BET2(L1)=BET2(L1)/CSCAT
  510 CONTINUE

!  Normalize PSD-Linearized output

      if ( .not. do_monodisperse ) then
         if ( .not. Do_psd_OldStyle .and. Do_LinearPSD ) then
            do kpsd = 1, 3
               L_CSCAT = LPSD_cscat(kpsd)
               DO 511 L1=1,L1MAX
                  LPSD_ALPH1(L1,kpsd) = &
                ( LPSD_ALPH1(L1,kpsd) - L_CSCAT * ALPH1(L1) ) / CSCAT
                  LPSD_ALPH2(L1,kpsd) = &
                ( LPSD_ALPH2(L1,kpsd) - L_CSCAT * ALPH2(L1) ) / CSCAT
                  LPSD_ALPH3(L1,kpsd) = &
                ( LPSD_ALPH3(L1,kpsd) - L_CSCAT * ALPH3(L1) ) / CSCAT
                  LPSD_ALPH4(L1,kpsd) = &
                ( LPSD_ALPH4(L1,kpsd) - L_CSCAT * ALPH4(L1) ) / CSCAT
                  LPSD_BET1(L1,kpsd) = &
                ( LPSD_BET1(L1,kpsd) - L_CSCAT * BET1(L1) ) / CSCAT
                  LPSD_BET2(L1,kpsd) = &
                ( LPSD_BET2(L1,kpsd) - L_CSCAT * BET2(L1) ) / CSCAT
  511          CONTINUE
            enddo
         endif
      endif

!  Normalize RFE-Linearized output

      if ( Do_rfe_Linearize ) then
         do LL = 1, nlin
            L_CSCAT = LRFE_cscat(LL)
            DO 531 L1=1,L1MAX
               LRFE_ALPH1(L1,LL) = &
             ( LRFE_ALPH1(L1,LL) - L_CSCAT * ALPH1(L1) ) / CSCAT
               LRFE_ALPH2(L1,LL) = &
             ( LRFE_ALPH2(L1,LL) - L_CSCAT * ALPH2(L1) ) / CSCAT
               LRFE_ALPH3(L1,LL) = &
             ( LRFE_ALPH3(L1,LL) - L_CSCAT * ALPH3(L1) ) / CSCAT
               LRFE_ALPH4(L1,LL) = &
             ( LRFE_ALPH4(L1,LL) - L_CSCAT * ALPH4(L1) ) / CSCAT
               LRFE_BET1(L1,LL) = &
             ( LRFE_BET1(L1,LL) - L_CSCAT * BET1(L1) ) / CSCAT
               LRFE_BET2(L1,LL) = &
             ( LRFE_BET2(L1,LL) - L_CSCAT * BET2(L1) ) / CSCAT
  531       CONTINUE
         enddo
      endif

!  First, do the Hovenier and van der Mee check

      CALL HOVENR ( L1MAX,ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2, & ! Inputs
                    fail_1, fail_2, message_1, message_2 )     ! Outputs

!  Exception handling on the check

      if ( fail_1 .or. fail_2 ) then
         fail = .true.
         if ( fail_1 ) message = TRIM(message_1)
         if ( fail_2 ) message = TRIM(message_2)
         trace   = 'VanderMee/Hovenier check failed'
         trace_2 = 'tmat_master_PLUS module' 
         istatus = 2
         return
      endif

!  Asymmetry parameter Assignations.

      tmat_asymm = alph1(2) / 3.0d0

      if ( .not. do_monodisperse ) then
         if ( .not. Do_psd_OldStyle .and. Do_LinearPSD ) then
            do kpsd = 1, 3
               LPSD_tmat_asymm(kpsd) = LPSD_alph1(2,kpsd) / 3.0d0
            enddo
         endif
      endif

      if ( Do_rfe_Linearize ) then
         do LL = 1, 3
            LRFE_tmat_asymm(LL) = LRFE_alph1(2,LL) / 3.0d0
         enddo
      endif

!  Expansion coefficients assignations

      tmat_ncoeffs = L1MAX
      DO 512 L1=1,L1MAX
         tmat_expcoeffs(L1,1) = ALPH1(L1)
         tmat_expcoeffs(L1,2) = ALPH2(L1)
         tmat_expcoeffs(L1,3) = ALPH3(L1)
         tmat_expcoeffs(L1,4) = ALPH4(L1)
         tmat_expcoeffs(L1,5) = BET1(L1)
         tmat_expcoeffs(L1,6) = BET2(L1)
  512 CONTINUE

!  PSD-Linearized

      if ( .not. do_monodisperse ) then
         if ( .not. Do_psd_OldStyle .and. Do_LinearPSD ) then
            do kpsd = 1, 3
               DO 515 L1=1,L1MAX
                  LPSD_tmat_expcoeffs(L1,1,kpsd) = LPSD_ALPH1(L1,kpsd)
                  LPSD_tmat_expcoeffs(L1,2,kpsd) = LPSD_ALPH2(L1,kpsd)
                  LPSD_tmat_expcoeffs(L1,3,kpsd) = LPSD_ALPH3(L1,kpsd)
                  LPSD_tmat_expcoeffs(L1,4,kpsd) = LPSD_ALPH4(L1,kpsd)
                  LPSD_tmat_expcoeffs(L1,5,kpsd) = LPSD_BET1(L1,kpsd)
                  LPSD_tmat_expcoeffs(L1,6,kpsd) = LPSD_BET2(L1,kpsd)
  515          CONTINUE
               LPSD_tmat_expcoeffs(1,1,kpsd) = 0.0d0
            enddo
         endif
      endif

!  RFE-Linearized

      if ( Do_rfe_Linearize ) then
         do LL = 1, 3
            DO 535 L1=1,L1MAX
               LRFE_tmat_expcoeffs(L1,1,LL) = LRFE_ALPH1(L1,LL)
               LRFE_tmat_expcoeffs(L1,2,LL) = LRFE_ALPH2(L1,LL)
               LRFE_tmat_expcoeffs(L1,3,LL) = LRFE_ALPH3(L1,LL)
               LRFE_tmat_expcoeffs(L1,4,LL) = LRFE_ALPH4(L1,LL)
               LRFE_tmat_expcoeffs(L1,5,LL) = LRFE_BET1(L1,LL)
               LRFE_tmat_expcoeffs(L1,6,LL) = LRFE_BET2(L1,LL)
  535       CONTINUE
            LRFE_tmat_expcoeffs(1,1,LL) = 0.0d0
         enddo
      endif

!  F-matrix calculation
!    Some redundancy here.          January 6th, 2011

      if ( Do_Fmatrix ) then

         LMAX=L1MAX-1

!  with no linearization

         if ( .not. Do_LinearPSD .and. .not. Do_rfe_Linearize ) then

            CALL MATR ( MAXNPA, NPNA, LMAX,                & ! Inputs
                     ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2,    & ! Inputs
                     Tmat_FMATRIX )                          ! Outputs
         endif

!  PSD-Linearized F matrix included

         if ( .not. do_monodisperse ) then
            if ( .not. Do_psd_OldStyle .and. Do_LinearPSD ) then
               CALL MATR_plus ( MAXNPA, NPNA, LMAX, 3,     & ! Inputs
                     ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2,    & ! Inputs
                     LPSD_ALPH1,LPSD_ALPH2,LPSD_ALPH3,     & ! Inputs
                     LPSD_ALPH4,LPSD_BET1,LPSD_BET2,       & ! Inputs
                     Tmat_FMATRIX, LPSD_Tmat_FMATRIX )       ! Outputs
            endif
         endif

!  RFE-Linearized F matrix included

         if ( Do_rfe_Linearize ) then
            CALL MATR_plus ( MAXNPA, NPNA, LMAX, NLIN,     & ! Inputs
                     ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2,    & ! Inputs
                     LRFE_ALPH1,LRFE_ALPH2,LRFE_ALPH3,     & ! Inputs
                     LRFE_ALPH4,LRFE_BET1,LRFE_BET2,       & ! Inputs
                     Tmat_FMATRIX, LRFE_Tmat_FMATRIX )       ! Outputs
         endif

!  End F-matrix clause

      endif

!  End coefficients clause

   endif

!  Finish

   return
end subroutine tmat_master_PLUS

!  End module

end module tmat_master_plus_m
