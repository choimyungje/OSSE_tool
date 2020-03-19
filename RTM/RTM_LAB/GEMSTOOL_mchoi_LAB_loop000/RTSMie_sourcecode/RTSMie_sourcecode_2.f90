MODULE RTSMie_sourcecode_m

use RTSMie_parameters_m
use RTSMie_distributions_m, ONLY : Mie_gauleg, rminmax, sizedis

private
public RTSMie_main

!  Following routines are here
!    RTSMie_main         - Main subroutine
!    mie_AnBncoeffs      - An + Bn + SP/Sm matrices
!    develop             - Develop coefficients
!    expand              - Make F-matrix from coefficient

contains

SUBROUTINE RTSMie_main                                       & 
       ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
         PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,   & ! I
         nblocks, nweights, xparticle_limit, R1R2_cutoff,    & ! I
         n_Fmatrix_angles, Fmatrix_angles,                   & ! I
         lambda, n_real, n_imag,                             & ! I
         Mie_bulk, Mie_asymm, Mie_ncoeffs,                   & ! O
         Mie_expcoeffs, Mie_Fmatrix, Mie_dist,               & ! O
         fail, istatus, message, trace, action )               ! O

!  Parameter module

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

!  PSD inputs
!  ----------

!  FixR1R2 : If  set, Use "R1R2_cutoff" for smallest particle size
!                     then Internal routine to calculate R1 and R2 (outputs)
!            If Not set, Use Input R1 and R2 for PSD limits.

   logical, intent(inout)  :: FixR1R2

!  R1, R2         - Minimum and Maximum radii (Microns)

   real    (KIND=dp), intent(inout)  :: R1, R2

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

   integer          , intent(in)  :: PSD_Index
   real    (KIND=dp), intent(in)  :: PSD_pars (3)

!  PSD quadrature control
!  ----------------------

!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.
!    R1R2_cutoff particle size for setting R1 and R2 internally

   INTEGER          , INTENT (IN) :: nblocks
   INTEGER          , INTENT (IN) :: nweights
   REAL    (KIND=dp), INTENT (IN) :: R1R2_cutoff

!  Optical: Wavelength, refractive index
!  -------------------------------------

!      LAMBDA         - wavelength of light (microns)
!      N_REAL, N_IMAG - real and imaginary parts, refractive index (N-i.GE.0)

   real    (KIND=dp), intent(in)  :: lambda, n_real, n_imag

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

   real(KIND=dp), intent(out) :: Mie_bulk (3)

!  Expansion coefficients and Asymmetry parameter, optional output

   integer      , intent(out) :: Mie_ncoeffs
   real(KIND=dp), intent(out) :: Mie_expcoeffs (6,0:max_Mie_angles)
   real(KIND=dp), intent(out) :: Mie_asymm

!  F-matrix, optional output

   real(KIND=dp), intent(out) :: Mie_Fmatrix(4,max_Mie_angles)

!  Distribution parameters
!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(KIND=dp), intent(out) :: Mie_dist (5)

!  Exception handling

   LOGICAL          , INTENT (OUT)   :: fail
   INTEGER          , INTENT (OUT)   :: istatus
   CHARACTER*(*)    , INTENT (OUT)   :: message, trace, action

!  Derived and Local Quantities
!  ============================
 
!  Complex refractive index

   COMPLEX (KIND=dp) :: m_complex

!  Qudrature for evaluating Coefficients

   INTEGER           :: n_coeffct_angles
   REAL    (KIND=dp) :: coeff_cosines(max_Mie_angles)
   REAL    (KIND=dp) :: coeff_weights(max_Mie_angles)

!  Local output from Mie_coeffs routine

   REAL    (KIND=dp), DIMENSION (max_Mie_sizes)                  :: q_ext, q_sca, asym
   COMPLEX (KIND=dp), DIMENSION (max_Mie_angles, max_Mie_sizes)  :: splus,sminus

!  Local Fmatrix for coefficients only

   REAL (KIND=dp), DIMENSION (4,max_Mie_angles)  :: Local_Fmatrix

!  Help variables

   CHARACTER*5 :: char5
   LOGICAL     :: do_angular_variation, do_polydisperse
   LOGICAL     :: failmm, faild
   INTEGER     :: i, angle, n_angles
   INTEGER     :: iblock, n_sizes, kf
 
   REAL    (KIND=dp)  :: factor_0, factor_1, d_pi, dtr
   REAL    (KIND=dp)  :: rstart, rfinis, help, rblock
   REAL    (KIND=dp)  :: quad, quadr2, quadr3, quadr4, mom2, mom3, mom4
   REAL    (KIND=dp)  :: ndens, gxsec, reff, volume, veff
   REAL    (KIND=dp)  :: Qext, Qsca, Qsca2, Qasy, ssalbedo
   REAL    (KIND=dp)  :: f(4)
   REAL    (KIND=dp)  :: angle_cosines(max_Mie_angles)
   REAL    (KIND=dp)  :: Fmatrix_cosines(max_Mie_angles)
   REAL    (KIND=dp), DIMENSION (max_Mie_sizes)   :: particle_sizes
   REAL    (KIND=dp), DIMENSION (max_Mie_sizes)   :: rquad, weights, nr

   COMPLEX (KIND=dp)  :: sp, sm, csp, csm
   COMPLEX (KIND=dp)  :: c_i, c_mi

   REAL    (KIND=dp)  :: xpart_root3, xparticle
   COMPLEX (KIND=dp)  :: y_argument

!  Maximum number of coefficients
!   (max_Mie_points is now a variable used to perform ALLOCATION)

   INTEGER            :: limmax, limstop, limsize
   INTEGER            :: max_Mie_points

!  Initial section
!  ===============

!  Output Zeroing
!  --------------

!  Main output

   Mie_dist   = d_zero
   Mie_bulk   = d_zero
   Mie_asymm  = d_zero

   Mie_Fmatrix    = d_zero
   Mie_ncoeffs    = 0
   Mie_expcoeffs  = d_zero

!  Exception handling

   trace   = ' '
   message = ' '
   action  = ' '
   fail    = .FALSE.
   istatus = 0

!  Constants
!  ---------

   c_i = ( 0.0_dp, 1.0_dp )
   c_mi = - c_i
   n_sizes = nweights

   d_pi = 4.0_dp * ATAN(d_one)
   dtr  = d_pi / 180.0_dp
   factor_0 = d_two * d_pi / lambda
   factor_1 = lambda / factor_0

   m_complex  = n_real * ( 1.0_dp, 0.0_dp ) + n_imag * c_i

!  Local values
!  ------------

   failmm = .false.

!  Local Bulk properties, initialize

   Qext = d_zero
   Qsca = d_zero
   Qasy = d_zero

!  Local Fmatrix, initialize

   Local_Fmatrix = d_zero

!  Local Distribution, initialize

   ndens  = d_zero
   gxsec  = d_zero ; mom2 = d_zero
   reff   = d_zero ; mom3 = d_zero
   veff   = d_zero ; mom4 = d_zero
   volume = d_zero

!  Local polydisperse flag

   do_polydisperse = .not. do_monodisperse

!  limiting radii calculation (Polydisperse)
!  -----------------------------------------

!  calculate R1 and R2, if the flag is set

   if ( FixR1R2 .and. do_polydisperse ) then
     CALL rminmax ( PSD_index, PSD_pars, R1R2_cutoff, R1, R2, message, failmm )
     IF ( failmm ) THEN
         fail    = .TRUE.;  istatus = 1
         trace   = 'Trace   : First Check in Mie Main. Failed to find radii extremes R1 and R2'
         action  = 'Action  : Consult with R. Spurr' 
         RETURN
      END IF
   endif

!  Check limiting radii if set externally (Polydisperse)

   if ( .not. FixR1R2 .and. do_polydisperse ) then
      if ( R1.lt.0.0d0 ) then
         failmm = .true.
         message = 'External R1 value < 0, out of bounds'
      else if ( R2 .le. 0.0d0 ) then
         failmm = .true.
         message = 'External R2 value =< 0, out of bounds'
      else if ( R1 .ge. R2 ) then
         failmm = .true.
         message = 'External R1 >= R2, Cannot be possible!'
      endif
      if ( failmm ) then
         fail    = .TRUE.;  istatus = 1
         trace   = 'Trace   : First Check in Mie Main. User R1/R2 wrong'
         action  = 'Action  : Change input values of R1 and R2' 
         RETURN
      END IF
   endif    

!  Limiting number of Mie Coefficients
!  -----------------------------------

!  number of blocks (Polydisperse). R2 = MonoRadius (Monodisperse)

   if (  do_polydisperse ) then
      rblock = ( R2 - R1 ) / DBLE(nblocks)
   else
      R2 = MonoRadius
   endif

!  limiting number of terms for coefficient computation
!    Applies to Mono or Polydisperse

   xparticle  = factor_0  * R2
   y_argument = xparticle * m_complex
   limstop    = 2
   IF ( xparticle > 0.02) THEN
      xpart_root3 = xparticle ** ( d_one / d_three )
      IF ( xparticle <= 8.0_dp ) THEN
         limstop = xparticle + 4.0_dp * xpart_root3 + d_two
      ELSE IF ( xparticle < 4200.0_dp ) THEN
         limstop = xparticle + 4.05_dp * xpart_root3 + d_two
      ELSE
         limstop = xparticle + 4.0_dp * xpart_root3 + d_two
      END IF
   END IF
   limmax = nint(max(DBLE(limstop),ABS(y_argument)) + 15.0_dp)

!  Set max_Mie_points for Memory allocation

   max_Mie_points = max(limstop,limmax)

!  Dimensioning and exception handling checks
!  ------------------------------------------

!  return if size limit exceeded

   IF ( xparticle > xparticle_limit ) THEN
      fail    = .TRUE.;  istatus = 1
      limsize = int(xparticle) + 1
      write(char5,'(i5)')limsize
      message = 'Message : error size parameter overflow'
      trace   = 'Trace   : Second check in Mie_main'
      action  = 'Action  : In configuration file, increase cutoff or '// &
                           'increase xparticle_limit to at least '//char5
      RETURN
   END IF

!  return if maximum number of terms too great 
!         (Not required now, Using Memory Allocation)

!  IF ( limstop > max_Mie_points )  THEN
!     failmie = .TRUE.
!     write(char5,'(i5)')limstop
!     message = 'Message : Insufficient dimensioning for maximum number of terms!'
!     trace   = 'Trace   : Third check in Mie_main'
!     action  = 'Action  : Increase max_Mie_points in calling program to at least '//char5
!     RETURN
!  END IF

!  And again, Dave recurrence
!         (Not required now, Using Memory Allocation)

!  IF ( limmax > max_Mie_points )  THEN
!     failmie = .TRUE.
!     write(char5,'(i5)')limmax
!     message = 'Message : Insufficient dimensioning for maximum number of terms (Dave recurrence)'
!     trace   = 'Trace   : Fourth check in Mie_main'
!     action  = 'Action  : Increase max_Mie_points in calling program to at least '//char5
!     RETURN
!  END IF

!  Compute the number of angles required for coefficient computation

   IF ( Do_Expcoeffs ) THEN
      n_coeffct_angles = 2*limmax + 2
      if ( n_coeffct_angles > max_Mie_angles ) then
         fail    = .TRUE.;  istatus = 1
         write(char5,'(i5)')n_coeffct_angles
         message = 'Message : Dimensioning error for number of terms for coefficient computation'
         trace   = 'Trace   : Third check in Mie Main'
         action  = 'Action  : Increase value of max_Mie_angles in calling program to at least '//char5
         return
      endif
   ENDIF

!  Compute the angles required for coefficient computation
!  -------------------------------------------------------

!  Quadrature (only need to do it once)

   IF ( Do_Expcoeffs ) THEN
      n_coeffct_angles = 2*limmax + 2
      CALL mie_gauleg ( max_Mie_angles, n_coeffct_angles, -1.0_dp, 1.0_dp, & ! Input
                        coeff_cosines, coeff_weights )                       ! Output
      Mie_ncoeffs = n_coeffct_angles
   END IF

! Overall cosines

   n_angles = 0
   do_angular_variation = .FALSE.

   IF ( Do_Expcoeffs ) THEN
      n_angles = n_coeffct_angles
      DO angle = 1, n_angles
         angle_cosines(angle) = coeff_cosines(angle)
      END DO
      do_angular_variation = .TRUE.
      IF ( Do_Fmatrix ) THEN
         DO angle = 1, n_Fmatrix_angles
            Fmatrix_cosines(angle) = dcos(Fmatrix_angles(angle)*dtr)
         END DO
      ENDIF
   ELSE
      IF ( Do_Fmatrix ) THEN
         n_angles = n_Fmatrix_angles
         DO angle = 1, n_angles
            angle_cosines(angle) = dcos(Fmatrix_angles(angle)*dtr)
         END DO
         do_angular_variation = .TRUE.
      END IF
   ENDIF

!  Monodisperse Calculation
!  ========================

   if ( do_Monodisperse ) then

!  particle size

      n_sizes           = 1
      particle_sizes(1) = factor_0 * MonoRadius

!  geometric cross-section, volume, effective radius/variance, density

      ndens  = d_one
      reff   = MonoRadius
      gxsec  = d_pi * reff * reff
      volume = (4.0_dp/3.0_dp) * gxsec * reff
 
!  Call to compute A_n and B_n coefficients and optical output
!    WARNING - easy to get segmentation fault before this call
!              Check dimensioning first, Use lots of memory    !!!!!!!

  
      CALL mie_AnBncoeffs                                        &
       ( max_Mie_angles, max_Mie_sizes, max_Mie_points,          & ! Dimensioning
         do_angular_variation,                                   & ! Input
         n_angles, n_sizes, m_complex,                           & ! Input
         particle_sizes, angle_cosines,                          & ! Input
         q_ext, q_sca, asym, splus, sminus )                       ! Output

!  Basic coefficients

      Qext  = Qext  + q_ext(1)
      Qsca  = Qsca  + q_sca(1)
      Qasy  = Qasy  + asym(1)

!  Local F-matrix

      IF ( do_angular_variation ) THEN
         Qsca2 =  Qsca / d_half 
         DO angle = 1, n_angles
            sp  = splus(angle,1)
            sm  = sminus(angle,1)
            csp = CONJG(sp)
            csm = CONJG(sm)
            f(1) =   REAL   ( sp * csp + sm * csm )
            f(2) = - REAL   ( sm * csp + sp * csm )
            f(3) =   REAL   ( sp * csp - sm * csm )
            f(4) =   REAL ( ( sm * csp - sp * csm ) * c_mi )
            DO kf = 1, 4
               Local_Fmatrix(kf,angle) = f(kf) / Qsca2
            ENDDO
         ENDDO
      END IF

!  asymmetry parameter and single scattering albedo

      Qasy     = d_two * Qasy / Qsca
      ssalbedo = Qsca/Qext

!  basic coefficients

      Qsca  = Qsca * factor_1
      Qext  = Qext * factor_1

!  Bulk properties

      MIE_BULK(1) = Qext
      MIE_BULK(2) = Qsca
      MIE_BULK(3) = ssalbedo
      MIE_ASYMM   = Qasy

!  Distribution properties

      MIE_DIST(1) = ndens
      MIE_DIST(2) = gxsec
      MIE_DIST(3) = volume
      MIE_DIST(4) = reff
      MIE_DIST(5) = veff

!  End Monodisperse

   ENDIF

!  Polydisperse
!  ============

   if ( do_Polydisperse ) then

!  start integration
!  -----------------

      DO iblock = 1, nblocks

!  Limiting radius for this block

         rstart = R1 + ( iblock-1) * rblock
         rfinis = rstart + rblock

!  Quadrature for this block

         CALL mie_gauleg ( max_Mie_sizes, n_sizes, rstart, rfinis, rquad, weights )

!  prepare particle sizes
    
         DO i = 1, n_sizes
            particle_sizes(i) = factor_0 * rquad(i)
         ENDDO

!  Call to compute A_n and B_n coefficients and optical output
!    WARNING - easy o get segmentation fault before this call
!              Check dimensioning first, Use lots of memory    !!!!!!!

         CALL mie_AnBncoeffs                                     &
       ( max_Mie_angles, max_Mie_sizes, max_Mie_points,          & ! Dimensioning
         do_angular_variation,                                   & ! Input
         n_angles, n_sizes, m_complex,                           & ! Input
         particle_sizes, angle_cosines,                          & ! Input
         q_ext, q_sca, asym, splus, sminus )                       ! Output
 
!  size distribution values

         CALL sizedis &
       ( max_Mie_sizes, PSD_index, PSD_pars, rquad, n_sizes, &
               nr,  message, faild )

         IF ( faild ) THEN
            fail = faild; istatus = 1
            write(char5,'(i5)')iblock
            trace  = 'Trace   : Fourth check in Mie_main. Subroutine sizedis failed for block number '//char5
            action = 'Action  : Consult with R. Spurr' 
            RETURN
         END IF


!  Integration over particle sizes within block
!  --------------------------------------------

         DO i = 1, n_sizes

!  Number density, geometric cross-section, 3rd and 4th powers

            quad   = nr(i) * weights(i)
            quadr2 = quad   * rquad(i) * rquad(i)
            quadr3 = quadr2 * rquad(i)
            quadr4 = quadr3 * rquad(i)
            ndens  = ndens + quad
            mom2   = mom2  + quadr2
            mom3   = mom3  + quadr3
            mom4   = mom4  + quadr4

!  Basic coefficients

            Qext  = Qext  + quad * q_ext(i)
            Qsca  = Qsca  + quad * q_sca(i)
            Qasy  = Qasy  + quad * asym(i)

!  angular variation loop

            IF ( do_angular_variation ) THEN
               DO angle = 1, n_angles
                  sp  = splus(angle,i)
                  sm  = sminus(angle,i)
                  csp = CONJG(sp)
                  csm = CONJG(sm)
                  f(1) =   REAL   ( sp * csp + sm * csm )
                  f(2) = - REAL   ( sm * csp + sp * csm )
                  f(3) =   REAL   ( sp * csp - sm * csm )
                  f(4) =   REAL ( ( sm * csp - sp * csm ) * c_mi )
                  DO kf = 1, 4
                     Local_Fmatrix(kf,angle) = Local_Fmatrix(kf,angle) + quad*f(kf) 
                  ENDDO
               END DO
            END IF
   
!  Finish PSD integration loops

         END DO
      END DO

!  Final Assignations
!  ------------------

!  Normalize F-matrix with scattering coefficient

      IF ( do_angular_variation ) THEN
         Qsca2 =  Qsca / d_half 
         DO angle = 1, n_angles
            DO kf = 1, 4
               Local_Fmatrix(kf,angle) = Local_Fmatrix(kf,angle) / Qsca2
            END DO
         END DO
      END IF

!  geometric cross-section

      gxsec   = d_pi * mom2 / ndens

!  asymmetry parameter and single scattering albedo

      Qasy     = d_two * Qasy / Qsca
      ssalbedo = Qsca/Qext

!  basic coefficients

      help  = factor_1 / ndens
      Qsca  = Qsca * help
      Qext  = Qext * help

!  geometrical quantities

      volume = (4.0_dp/3.0_dp) * d_pi * mom3 / ndens
      reff   = mom3 / mom2

!  Variance output

      help  = mom2 / mom3 / mom3 
      veff  = help * mom4
      veff  = veff - d_one

!  Bulk properties

      MIE_BULK(1) = Qext
      MIE_BULK(2) = Qsca
      MIE_BULK(3) = ssalbedo
      MIE_ASYMM   = Qasy

!  Distribution properties

      MIE_DIST(1) = ndens
      MIE_DIST(2) = gxsec
      MIE_DIST(3) = volume
      MIE_DIST(4) = reff
      MIE_DIST(5) = veff

!  End polydisperse

   ENDIF

!  Coefficient development and expansion
!  =====================================

!  Call to develop expansion coefficients

   if ( do_Expcoeffs ) then
      CALL develop                                                      & !----DEVELOP CALL
         ( max_Mie_angles, n_coeffct_angles, n_coeffct_angles,          & ! I
           coeff_cosines, coeff_weights, Local_Fmatrix, Mie_Expcoeffs )   ! I/O
   endif

!  F-matrix output.
!    If Coefficients set, use EXPAND routine to generate F-matrix.
!    Otherwise, just copy the already-calculated F-matrix.

   if ( do_Fmatrix ) then
      if ( do_Expcoeffs ) then
         CALL expand                                                    & !----EXPAND CALL
           ( max_Mie_angles, n_coeffct_angles, n_Fmatrix_angles,        & ! I
             Fmatrix_cosines, Mie_Expcoeffs, Mie_Fmatrix )                ! I/O 
      else
         Mie_Fmatrix(1:4,1:n_Fmatrix_angles) = Local_Fmatrix(1:4,1:n_Fmatrix_angles)
      endif
   endif

!  Finish
!  ======

   RETURN
END SUBROUTINE RTSMie_main

!

SUBROUTINE mie_AnBncoeffs                                        & !
       ( max_Mie_angles, max_Mie_sizes, max_Mie_points,          & ! Dimensioning
         do_angular_variation, n_angles, n_sizes, m_complex,     & ! Input
         particle_sizes, angle_cosines,                          & ! Input
         q_ext,  q_sca,  asym,  splus,  sminus )                   ! Output

! name:
!       mie_AnBncoeffs

! purpose:
!       calculates the scattering parameters of a series of particles
!       using the mie scattering theory. FOR USE WITH POLYDISPERSE COD

! inputs:
!       particle_sizes:       array of particle size parameters
!       angle_cosines:        array of angle cosines
!       m_complex:            the complex refractive index of the particles
!       n_angles, n_sizes:    number of scattering angles, number of particle sizes
!       do_angular_variation: flag for S+/S- output


! outputs (1):
!       q_ext:       the extinction efficiency
!       q_sca:       the scattering efficiency
!       asym:        the asymmetry parameter
!       splus:       the first amplitude function
!       sminus:      the second amplitude function

! modification history
!       g. thomas IDL Mie code       (February  2004). Basic Monodisperse derivatives.
!       r. spurr  F90 Mie code       ( October  2004). Extension all Polydisperse derivatives.
!       r. spurr  Exception handling (September 2008). Exception handling removed.

!  modules

  USE RTSMie_parameters_m, ONLY : dp, d_zero, d_half, d_one, d_two, d_three

!  implicit none statement

  IMPLICIT NONE

!  Dimensioning input

  INTEGER, INTENT (IN) :: max_Mie_angles, max_Mie_sizes      !  Fixed
  INTEGER, INTENT (IN) :: max_Mie_points                     !  Allocatable

!  input

  LOGICAL          , INTENT (IN) :: do_angular_variation
  INTEGER          , INTENT (IN) :: n_angles, n_sizes 
  COMPLEX (KIND=dp), INTENT (IN) :: m_complex

  REAL    (KIND=dp), DIMENSION (max_Mie_sizes),  INTENT (IN) :: particle_sizes
  REAL    (KIND=dp), DIMENSION (max_Mie_angles), INTENT (IN) :: angle_cosines

!  output (1)

  REAL    (KIND=dp), DIMENSION (max_Mie_sizes),                 INTENT (OUT) :: q_ext, q_sca, asym
  COMPLEX (KIND=dp), DIMENSION (max_Mie_angles, max_Mie_sizes), INTENT (OUT) :: splus, sminus

!  local variables for Mie code

  INTEGER            :: size, angle, n, nm1, nstop(max_Mie_sizes), nmax, maxstop
  REAL    (KIND=dp)  :: xparticle, xpart_root3
  REAL    (KIND=dp)  :: xinv, xinvsq, two_d_xsq, xinv_dx
  REAL    (KIND=dp)  :: dn, dnp1, dnm1, dnsq, dnnp1, tnp1, tnm1, hnp1, hnm1
  REAL    (KIND=dp)  :: cos_x, sin_x, psi0, psi1, chi0, chi1, psi, chi
  REAL    (KIND=dp)  :: s, t, tau_n, factor, forward, bckward

  COMPLEX (KIND=dp)  :: inverse_m, y_argument, yinv, yinvsq, a1, zeta, zeta1
  COMPLEX (KIND=dp)  :: an, bn, an_star, bn_star, anm1, bnm1, bnm1_star
  COMPLEX (KIND=dp)  :: biga_divs_m, biga_mult_m, noverx, aterm, bterm
  COMPLEX (KIND=dp)  :: facplus, facminus, c_zero, c_one, c_i, c_mi
  COMPLEX (KIND=dp)  :: an_denom, bn_denom

!  redundant variables
!  REAL    (KIND=dp), INTENT (IN) :: xparticle_limit           ! subroutine argument
!  REAL    (KIND=dp)  :: four_d_xsq
!  CHARACTER*4        :: char4
!  COMPLEX (KIND=dp)  :: common, an_denom_dm, s1, s2
!  INTEGER            :: nmax_end

!  local arrays. Note Use of Allocatables, 25 March 2011

  REAL    (KIND=dp), DIMENSION (max_Mie_angles)     :: pi_n, pi_nm1
  COMPLEX (KIND=dp), DIMENSION (:),   ALLOCATABLE   :: biga
  REAL    (KIND=dp), DIMENSION (:,:), ALLOCATABLE   :: polyplus, polyminus

!  COMPLEX (KIND=dp), DIMENSION (max_Mie_points)              :: biga
!  REAL  (KIND=dp), DIMENSION (max_Mie_angles,max_Mie_points) :: polyplus, polyminus

!  Initial section
!  ---------------

!  Allocate memory

  ALLOCATE (polyplus (max_Mie_angles,max_Mie_points))
  ALLOCATE (polyminus(max_Mie_angles,max_Mie_points))
  ALLOCATE (biga(max_Mie_points))

!  Complex Constants

  c_zero = ( 0.0_dp, 0.0_dp )
  c_one  = ( 1.0_dp, 0.0_dp )
  c_i    = ( 0.0_dp, 1.0_dp )
  c_mi   = -c_i

!  assign number of terms and maximum 

  maxstop = 0
  DO size = 1, n_sizes
    xparticle  = particle_sizes (size)
    IF ( xparticle < 0.02) THEN
      nstop(size) = 2
    ELSE
      xpart_root3 = xparticle ** ( d_one / d_three )
      IF ( xparticle <= 8.0_dp ) THEN
        nstop(size) = xparticle + 4.0_dp * xpart_root3 + d_two
      ELSE IF ( xparticle < 4200.0_dp ) THEN
        nstop(size) = xparticle + 4.05_dp * xpart_root3 + d_two
      ELSE
      nstop(size) = xparticle + 4.0_dp * xpart_root3 + d_two
      END IF
    END IF
    maxstop = max(nstop(size),maxstop)
  END DO

!  phase function expansion polynomials
!    ---> initialise phase function Legendre polynomials
!    ---> Recurrence phase function Legendre polynomials

  IF ( do_angular_variation ) THEN

    DO angle = 1, n_angles
      pi_nm1(angle) = d_zero
      pi_n(angle)   = d_one
    END DO

    DO n = 1, maxstop
      nm1 = n - 1
      dn    = dble(n)
      dnp1  = dn   + d_one
      forward = dnp1 / dn
      DO angle = 1, n_angles
        s = angle_cosines(angle) * pi_n(angle)
        t = s - pi_nm1(angle)
        tau_n = dn*t - pi_nm1(angle)
        polyplus(angle,n)  = pi_n(angle) + tau_n
        polyminus(angle,n) = pi_n(angle) - tau_n
        pi_nm1(angle) = pi_n(angle)
        pi_n(angle)   = s + t*forward
      END DO
    END DO

  END IF

!  start loop over particle sizes
!  ------------------------------

  DO size = 1, n_sizes

!  initialize output

    asym(size)  = d_zero
    q_ext(size) = d_zero
    q_sca(size) = d_zero

!  some auxiliary quantities

    xparticle  = particle_sizes (size)
    xinv       = d_one / xparticle
    xinvsq     = xinv * xinv
    two_d_xsq  = d_two * xinvsq
    xinv_dx    = - d_two * xinv
 
    inverse_m  =     c_one / m_complex
    y_argument = xparticle * m_complex
    yinv       = d_one / y_argument 
    yinvsq     = yinv * yinv

!  Biga = ratio derivative, recurrence due to J. Dave

    nmax = nint(max(dble(nstop(size)),abs(y_argument)) + 15.0_dp)
    biga(nmax) = c_zero
    DO n = nmax-1, 1,-1
       a1      = dble(n+1) / y_argument
       biga(n) = a1 - c_one / (a1+biga(n+1))
    END DO

! initialize Riccati-Bessel functions

    tnp1 = d_one
    cos_x = COS(xparticle)
    sin_x = SIN(xparticle)
    psi0 = cos_x
    psi1 = sin_x
    chi1 =-cos_x
    chi0 = sin_x
    zeta1 = CMPLX(psi1,chi1,kind=dp)

!  initialise sp and sm

    IF ( do_angular_variation ) THEN
      DO angle = 1, n_angles
         splus(angle,size)  = c_zero
         sminus(angle,size) = c_zero
      END DO
    END IF

!  main loop

    DO n = 1, nstop(size)

!  various factors

      dn    = dble(n)
      dnp1  = dn   + d_one
      dnm1  = dn   - d_one
      tnp1  = tnp1 + d_two
      tnm1  = tnp1 - d_two

      dnsq  = dn * dn
      dnnp1 = dnsq + dn
      factor  = tnp1 / dnnp1
      bckward = dnm1 / dn

!  Ricatti - Bessel recurrence

      psi = tnm1 * psi1/xparticle - psi0
      chi = tnm1 * chi1/xparticle - chi0
      zeta = CMPLX(psi,chi,kind=dp)

!  a(n) and b(n)

      biga_divs_m = biga(n) * inverse_m
      biga_mult_m = biga(n) * m_complex
      noverx = CMPLX(dn/xparticle,d_zero,kind=dp)
      aterm = biga_divs_m + noverx
      bterm = biga_mult_m + noverx

      an_denom = (aterm * zeta - zeta1)
      bn_denom = (bterm * zeta - zeta1)
      an =   ( aterm*psi-psi1 ) / an_denom
      bn =   ( bterm*psi-psi1 ) / bn_denom
      an_star = CONJG(an)
      bn_star = CONJG(bn)

!  basic coefficients
!  ------------------

!  Q coefficients

      q_ext(size) = q_ext(size) + tnp1 * REAL ( an + bn )   
      q_sca(size) = q_sca(size) + tnp1 * REAL ( an*CONJG(an) + bn*CONJG(bn) )

!  asymmetry parameter

      IF ( n > 1 ) THEN
        hnp1 = bckward * dnp1
        hnm1 = tnm1  / (dnsq - dn)
        asym(size) = asym(size) &
           +  hnp1 * REAL ( anm1*an_star + bnm1*bn_star) &
           +  hnm1 * REAL ( anm1*bnm1_star) 
      END IF

!  Upgrades
!  --------

!  upgrade an/bn recurrences (only for asymmetry parameter)

      anm1      = an
      bnm1      = bn
      bnm1_star = bn_star

!  upgrade Ricatti-Bessel recurrences

      psi0 = psi1
      psi1 = psi
      chi0 = chi1
      chi1 = chi
      zeta1 = CMPLX(psi1,chi1,kind=dp)

!  S+/S- function stuff
!  --------------------

      IF ( do_angular_variation ) THEN
        facplus  = factor * ( an + bn )
        facminus = factor * ( an - bn )
        DO angle = 1, n_angles
          splus(angle,size)  = splus(angle,size)  + facplus  * polyplus(angle,n)
          sminus(angle,size) = sminus(angle,size) + facminus * polyminus(angle,n)
        END DO
      END IF

!  end sum loop

    END DO

!  End loop and finish
!  -------------------

!  end loop over particle sizes

  END DO

!  de-allocate the large arrays

  DEALLOCATE ( polyplus, polyminus, biga )

!  finish

  RETURN
END SUBROUTINE mie_AnBncoeffs


SUBROUTINE develop ( max_Mie_angles, ncoeffs, nangles, &
                     cosines, weights, FMAT, expcoeffs )

!  Based on the Meerhoff Mie code

!************************************************************************
!*  Calculate the expansion coefficients of the scattering matrix in    *
!*  generalized spherical functions by numerical integration over the   *
!*  scattering angle.                                   *
!************************************************************************

!  modules

  USE RTSMie_parameters_m, ONLY : dp, d_zero, d_half, d_one, d_two, d_three, d_four

!  implicit none statement

  IMPLICIT NONE

!  input

  INTEGER          , INTENT (IN) :: max_Mie_angles
  INTEGER          , INTENT (IN) :: ncoeffs, nangles
 
  REAL    (KIND=dp), INTENT (IN) :: cosines(max_Mie_angles)
  REAL    (KIND=dp), INTENT (IN) :: weights(max_Mie_angles)
  REAL    (KIND=dp), INTENT (IN) :: FMAT(4,max_Mie_angles)

!  output. Already initialized

  REAL    (KIND=dp), INTENT (INOUT) :: expcoeffs(6,0:max_Mie_angles)

!  local variables

  REAL    (KIND=dp) :: P00(max_Mie_angles,2)
  REAL    (KIND=dp) :: P02(max_Mie_angles,2)
  REAL    (KIND=dp) :: P22(max_Mie_angles,2)
  REAL    (KIND=dp) :: P2m2(max_Mie_angles,2)
  REAL    (KIND=dp) :: fmatw(4,max_Mie_angles)

  INTEGER           :: i, j, l, lnew, lold, itmp
  INTEGER           :: index_11, index_12, index_22, index_33, index_34, index_44 
  REAL    (KIND=dp) :: dl, dl2, qroot6, fac1, fac2, fac3, fl,&
                       sql4, sql41, twol1, tmp1, tmp2, denom, &
                       alfap, alfam, ps

!  Initialization

  qroot6 = -0.25_dp*SQRT(6.0_dp)

!  Old indexing

!  ps = 1.0d0
!  index_11 = 1
!  index_12 = 2
!  index_22 = 3
!  index_33 = 4
!  index_34 = 5
!  index_44 = 6

!  New indexing consistent with Tmatrix output

  ps = -1.0d0
  index_11 = 1
  index_22 = 2
  index_33 = 3
  index_44 = 4
  index_12 = 5
  index_34 = 6

!  Multiply the scattering matrix F with the weights w for all angles  *
!  We do this here because otherwise it should be done for each l      *

  DO i = 1, 4
    DO j = 1, nangles
      fmatw(i,j) = weights(j)*FMAT(i,j)
    END DO
  END DO

!  Start loop over the coefficient index l                             *
!  first update generalized spherical functions, then calculate coefs. *
!  lold and lnew are pointer-like indices used in recurrence           *

  lnew = 1
  lold = 2

  DO l = 0, ncoeffs

    IF (l == 0) THEN

      dl   = d_zero
      DO  i=1, nangles
        P00(i,lold) = d_one
        P00(i,lnew) = d_zero
        P02(i,lold) = d_zero
        P22(i,lold) = d_zero
        P2m2(i,lold)= d_zero
        P02(i,lnew) = d_zero
        P22(i,lnew) = d_zero
        P2m2(i,lnew)= d_zero
      END DO

    ELSE

      dl   = DBLE(l)
      dl2  = dl * dl
      fac1 = (d_two*dl-d_one)/dl
      fac2 = (dl-d_one)/dl
      DO  i=1, nangles
        P00(i,lold) = fac1*cosines(i)*P00(i,lnew) - fac2*P00(i,lold)
      END DO

    ENDIF

    IF (l == 2) THEN

      DO  i=1, nangles
        P02(i,lold) = qroot6*(d_one-cosines(i)*cosines(i))
        P22(i,lold) = 0.25_dp*(d_one+cosines(i))*(d_one+cosines(i))
        P2m2(i,lold)= 0.25_dp*(d_one-cosines(i))*(d_one-cosines(i))
        P02(i,lnew) = d_zero
        P22(i,lnew) = d_zero
        P2m2(i,lnew)= d_zero
      END DO
      sql41 = d_zero

    ELSE IF (l > 2) THEN

      sql4  = sql41
      sql41 = dsqrt(dl2-d_four)
      twol1 = 2.D0*dl - d_one
      tmp1  = twol1/sql41
      tmp2  = sql4/sql41
      denom = (dl-d_one)*(dl2-d_four)
      fac1  = twol1*(dl-d_one)*dble(l)/denom
      fac2  = 4.D0*twol1/denom
      fac3  = dl*((dl-d_one)*(dl-d_one)-d_four)/denom
      DO i=1, nangles
        P02(i,lold) = tmp1*cosines(i)*P02(i,lnew)         - tmp2*P02(i,lold)
        P22(i,lold) = (fac1*cosines(i)-fac2)*P22(i,lnew)  - fac3*P22(i,lold)
        P2m2(i,lold)= (fac1*cosines(i)+fac2)*P2m2(i,lnew) - fac3*P2m2(i,lold)
      END DO

    END IF

    itmp = lnew
    lnew = lold
    lold = itmp
    alfap = d_zero
    alfam = d_zero

    fl = dl+d_half
    do i=1, nangles
      expcoeffs(index_11,l) = expcoeffs(index_11,l) + P00(i,lnew)*fmatw(1,i)
      alfap = alfap + P22(i,lnew)  * (fmatw(1,i)+fmatw(3,i))
      alfam = alfam + P2m2(i,lnew) * (fmatw(1,i)-fmatw(3,i))
      expcoeffs(index_44,l) = expcoeffs(index_44,l) + P00(i,lnew)*fmatw(3,i)
      expcoeffs(index_12,l) = expcoeffs(index_12,l) + P02(i,lnew)*fmatw(2,i)
      expcoeffs(index_34,l) = expcoeffs(index_34,l) + P02(i,lnew)*fmatw(4,i)
    END DO
    expcoeffs(index_11,l) =  fl*expcoeffs(index_11,l)
    expcoeffs(index_22,l) =  fl*d_half*(alfap+alfam)
    expcoeffs(index_33,l) =  fl*d_half*(alfap-alfam)
    expcoeffs(index_44,l) =  fl*expcoeffs(index_44,l)
    expcoeffs(index_12,l) =  fl*expcoeffs(index_12,l)
    expcoeffs(index_34,l) =  ps * fl*expcoeffs(index_34,l)
  END DO

!  Phase function normalization

  expcoeffs(index_11,0)        = d_one

!  Finish

  RETURN
END SUBROUTINE develop



SUBROUTINE expand ( max_Mie_angles, ncoeffs, nangles, cosines, expcoeffs, FMAT )

!  Based on the Meerhoff Mie code

!  Use the expansion coefficients of the scattering matrix in 
!  generalized spherical functions to expand F matrix

!  modules

  USE RTSMie_parameters_m, ONLY : dp, d_zero, d_one, d_two, d_four

!  implicit none statement

  IMPLICIT NONE

!  input

  INTEGER          , INTENT (IN) :: max_Mie_angles
  INTEGER          , INTENT (IN) :: ncoeffs, nangles
  REAL    (KIND=dp), INTENT (IN) :: cosines(max_Mie_angles)
  REAL    (KIND=dp), INTENT (IN) :: expcoeffs(6,0:max_Mie_angles)

!  output, already initialized

  REAL    (KIND=dp), INTENT (INOUT) :: FMAT(4,max_Mie_angles)

!  local variables

  REAL    (KIND=dp) :: P00(max_Mie_angles,2)
  REAL    (KIND=dp) :: P02(max_Mie_angles,2)

  INTEGER           :: i, l, lnew, lold, itmp, findex(4)
  INTEGER           :: index_11, index_12, index_34, index_44 
  REAL    (KIND=dp) :: dl, qroot6, fac1, fac2, sql4, sql41, tmp1, tmp2

!  Initialization

  qroot6 = -0.25_dp*SQRT(6.0_dp)

!  Old indexing

!  index_11 = 1 ; findex(1) = 1
!  index_12 = 2 ; findex(2) = 2
!  index_44 = 6 ; findex(3) = 3
!  index_34 = 5 ; findex(4) = 4

!  New indexing

  index_11 = 1 ; findex(1) = 1
  index_12 = 5 ; findex(2) = 3
  index_44 = 4 ; findex(3) = 2
  index_34 = 6 ; findex(4) = 4

!  Start loop over the coefficient index l
!  first update generalized spherical functions, then calculate coefs.
!  lold and lnew are pointer-like indices used in recurrence 

  lnew = 1
  lold = 2

  DO l = 0, ncoeffs

     IF ( l == 0) THEN

!  Adding paper Eqs. (76) and (77) with m=0

        DO i=1, nangles
           P00(i,lold) = d_one
           P00(i,lnew) = d_zero
           P02(i,lold) = d_zero
           P02(i,lnew) = d_zero
        END DO

     ELSE

        dl   = DBLE(l)
        fac1 = (d_two*dl-d_one)/dl
        fac2 = (dl-d_one)/dl

! Adding paper Eq. (81) with m=0

        DO i=1, nangles
           P00(i,lold) = fac1*cosines(i)*P00(i,lnew) - fac2*P00(i,lold)
        END DO

     END IF

     IF ( l == 2) THEN

! Adding paper Eq. (78)  
! sql4 contains the factor dsqrt((l+1)*(l+1)-4) needed in
! the recurrence Eqs. (81) and (82)

        DO i=1, nangles
           P02(i,lold) = qroot6*(d_one-cosines(i)*cosines(i))
           P02(i,lnew) = d_zero
        END DO
        sql41 = d_zero

     ELSE IF ( l > 2) THEN

! Adding paper Eq. (82) with m=0

        sql4  = sql41
        sql41 = dsqrt(dl*dl-d_four)
        tmp1  = (d_two*dl-d_one)/sql41
        tmp2  = sql4/sql41

        DO i=1, nangles
           P02(i,lold) = tmp1*cosines(i)*P02(i,lnew) - tmp2*P02(i,lold)
        END DO

     END IF

! Switch indices so that lnew indicates the function with
! the present index value l, this mechanism prevents swapping
! of entire arrays.

     itmp = lnew
     lnew = lold
     lold = itmp

! Now add the l-th term to the scattering matrix.
! See de Haan et al. (1987) Eqs. (68)-(73).
! Remember for Mie scattering : F11 = F22 and F33 = F44

     DO i=1, nangles
        FMAT(findex(1),i) = FMAT(findex(1),i) + expcoeffs(index_11,l)*P00(i,lnew)
        FMAT(findex(2),i) = FMAT(findex(2),i) + expcoeffs(index_12,l)*P02(i,lnew)
        FMAT(findex(3),i) = FMAT(findex(3),i) + expcoeffs(index_44,l)*P00(i,lnew)
        FMAT(findex(4),i) = FMAT(findex(4),i) + expcoeffs(index_34,l)*P02(i,lnew)
     END DO

  END DO

  RETURN
END SUBROUTINE expand

!  End Module

end MODULE RTSMie_sourcecode_m

