
module FO_Vector_RTCalcs_m

!  For given wavelength, routine calculates First-Order upwelling+downwelling Stokes vector
!     (1) For the Atmospheric Single-scatter and Surface Direct-Beam solar sources,
!     (2) For the Atmospheric and Surface Thermal emission Direct solutions

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS (LOS-path sphericity) calculations
!  This will perform Regular-PS  (average-secant pseudo-spherical) calculations

!  This is Version 1.2, without Partials. Code is stand alone with no dependencies.
!    01 December 2011, R. Spurr, RT Solutions Inc.
!    13 February 2012, R. Spurr, RT Solutions Inc.
!    21 June     2012, R. Spurr, RT Solutions Inc.

!  For Solar sources, the subroutines are
!       SSV_Integral_1_UP   (Upwelling only)
!       SSV_Integral_1_DN   (Downwelling only)
!       SSV_Integral_1_UPDN (Upwelling and Downwelling)

!  For Thermal Emission sources, the subroutines are
!       DTEV_Integral_1_UP   (Upwelling only)
!       DTEV_Integral_1_DN   (Downwelling only)
!       DTEV_Integral_1_UPDN (Upwelling and Downwelling)

!  All subroutines public

public

contains


subroutine SSV_Integral_1_UP &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels, do_sunlight,   & ! Inputs (dimension)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, do_lambertian, doNadir,   & ! Inputs (Flags)
     nstokes, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,    & ! Inputs (control)
     reflec, extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,      & ! Inputs (Optical)
     Mu0, Mu1, GenSpher, Rotations, NCrit, xfine, wfine, csqfine, cotfine,       & ! Inputs (Geometry)
     Raycon, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,           & ! Inputs (Geometry)
     stokes_up, stokes_db, cumsource_up )                                          ! Outputs

!  Stand alone routine for SS field with Solar sources alone
!    No partials

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

      INTEGER, Intent(in) :: maxlayers
      INTEGER, Intent(in) :: maxfinelayers
      INTEGER, Intent(in) :: maxmoments_input
      INTEGER, Intent(in) :: max_user_levels

!  flags

      LOGICAL, Intent(in) ::  DO_SUNLIGHT
      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DO_LAMBERTIAN
      LOGICAL, Intent(in) ::  DONADIR

!  Numbers

      INTEGER, Intent(in) ::  NSTOKES
      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  NMOMENTS_INPUT

      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

      REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
      REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
      REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
      REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
      REAL(fpk), Intent(in) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4 )

!  Solar Flux and Surface reflectivity (Could be the albedo)

      REAL(fpk), Intent(in) :: REFLEC(4,4), FLUX, FLUXVEC(4)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      integer  , Intent(in)  :: NCrit
      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers)
      real(fpk), Intent(in)  :: Mu0, Mu1

!  solar paths 

      integer  , Intent(in)  :: ntraverse  (0:maxlayers)
      real(fpk), Intent(in)  :: sunpaths   (0:maxlayers,maxlayers)
      integer  , Intent(in)  :: ntraverse_fine(maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)

!  Generalized spherical functions.
!    Rotations(1-4)    = C1, S1, C2, S2

      REAL(fpk), Intent(in)  :: GenSpher(0:maxmoments_input, 4)
      REAL(fpk), Intent(in)  :: Rotations(4)

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: stokes_up     ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: stokes_db     ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: cumsource_up  ( 0:maxlayers,     4 )

!  LOCAL
!  -----

!  Attenuations

      real(fpk)  :: suntau            (0:maxlayers)
      real(fpk)  :: attenuations      (0:maxlayers)
      real(fpk)  :: attenuations_fine (maxlayers,maxfinelayers)

!  Scattering

      real(fpk)  :: tms (maxlayers)
      real(fpk)  :: exactscat_up (maxlayers,4,4)

!  Source function integration results

      real(fpk)  :: sources_up       ( maxlayers, 4 )
      real(fpk)  :: lostrans_up      ( maxlayers )
      real(fpk)  :: multiplier_up    ( maxlayers )

!  Help

      integer    :: n, ns, uta, nstart, nc, nut, nut_prev, j, k, L, o1, o2
      logical    :: layermask_up(maxlayers)

      real(fpk)  :: sumd, help, sum, tran_1, tran, func, kn, ke, xjkn
      real(fpk)  :: cot_1, cot_2, lostau, factor1, factor2
      real(fpk)  :: cumsource_db(4), fmat(6), sum23, dif23
      real(fpk)  :: help3c1, help3s1, help4c1, help4s1, m4

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      CUMSOURCE_UP  = zero ; STOKES_UP  = zero ; STOKES_DB    = zero
      lostrans_up   = zero ; sources_up = zero ; exactscat_up = zero ; multiplier_up = zero

!  Bookkeeping

      ns = nstokes ; fmat = zero
      NUT = USER_LEVELS(1) + 1
      LAYERMASK_UP = .false.
      LAYERMASK_UP(NUT:NLAYERS) = .true.

!  TMS factors

      do n = 1, nlayers
         if ( do_deltam_scaling ) then
            help = one - truncfac(n) * omega(n)
            tms(n) = omega(n) / help
         else
            tms(n) = omega(n)
         endif
      enddo

!  Scattering functions
!  --------------------

!  Scalar only

      if ( nstokes .eq. 1 ) then
         do n = 1, nlayers
            if ( layermask_up(n) ) then
              sum = zero
               do L = 0, nmoments_input
                  sum = sum + GenSpher(L,1) * Greekmat(n,L,1,1)
               enddo
               exactscat_up(n,1,1) = sum * tms(n)
            endif
         enddo
      endif

!  Vector with Sunlight

      if ( nstokes .gt. 1 .and. do_sunlight ) then
         do n = 1, nlayers
            if ( layermask_up(n) ) then
               fmat(1:2) = zero
               do L = 0, nmoments_input
                  fmat(1) = fmat(1) + Genspher(L,1) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2) * greekmat(n,L,1,2)
               enddo
               exactscat_up(n,1,1) = + fmat(1)
               exactscat_up(n,2,1) = - fmat(2) * Rotations(3)
               exactscat_up(n,3,1) = + fmat(2) * Rotations(4)
               exactscat_up(n,1:ns,1) = tms(n) * exactscat_up(n,1:ns,1) 
            endif
         enddo
      endif

!  Vector General case, USE FULL 4X4 MATRIX
!    CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

      if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
         do n = 1, nlayers
            if ( layermask_up(n) ) then
               fmat = zero
               do L = 0, nmoments_input
                  fmat(1) = fmat(1) + Genspher(L,1) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2) * greekmat(n,L,1,2)
                  sum23 = greekmat(n,L,2,2) + greekmat(n,L,3,3)
                  dif23 = greekmat(n,L,2,2) - greekmat(n,L,3,3)
                  fmat(3) = fmat(3) + Genspher(L,3) * sum23
                  fmat(4) = fmat(4) + Genspher(L,4) * dif23
               enddo
               fmat(3) = ( fmat(3) + fmat(4) ) * 0.5_fpk
               fmat(4) = ( fmat(3) - fmat(4) )
               if ( nstokes.eq.4) then
                  do L = 0, nmoments_input
                     fmat(5) = fmat(5) + Genspher(L,2) * greekmat(n,L,3,4)
                     fmat(6) = fmat(6) + Genspher(L,1) * greekmat(n,L,4,4)
                  enddo
               endif
               help3c1 = fmat(3) * Rotations(1)
               help3s1 = fmat(3) * Rotations(2)
               help4c1 = fmat(4) * Rotations(1)
               help4s1 = fmat(4) * Rotations(2)
               exactscat_up(n,1,1) = + fmat(1)
               exactscat_up(n,2,1) = - fmat(2) * Rotations(3)
               exactscat_up(n,3,1) = + fmat(2) * Rotations(4)
               exactscat_up(n,1,2) = + fmat(2) * Rotations(1)
               exactscat_up(n,1,3) = - fmat(2) * Rotations(2)
               exactscat_up(n,2,2) = + help3c1 * Rotations(3) - help4s1 * Rotations(4)
               exactscat_up(n,2,3) = - help3s1 * Rotations(3) - help4c1 * Rotations(4)
               exactscat_up(n,3,2) = + help3c1 * Rotations(4) + help4s1 * Rotations(3)
               exactscat_up(n,3,3) = - help3s1 * Rotations(4) + help4c1 * Rotations(3)
               if ( nstokes .eq. 4 ) then
                  exactscat_up(n,2,4) = - fmat(5) * Rotations(4) 
                  exactscat_up(n,4,2) = - fmat(5) * Rotations(2) 
                  exactscat_up(n,3,4) = + fmat(5) * Rotations(3) 
                  exactscat_up(n,4,3) = - fmat(5) * Rotations(1) 
                  exactscat_up(n,4,4) = + fmat(6)
               endif
               exactscat_up(n,1:ns,1:ns) = tms(n)*exactscat_up(n,1:ns,1:ns)
            endif
         enddo
      endif

!  Attenuations
!  ============

!  Initialize, only to layer Ncrit if applicable

      Attenuations = zero ; Attenuations_fine = zero ; Suntau = zero 
      nstart = nlayers ; if (Ncrit.ne.0) nstart = nCrit


!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         sumd = zero
         do k = 1, ntraverse(n)
            sumd = sumd + extinction(k) * sunpaths(n,k)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
!         if ( n.gt.87 ) write(*,*)n,sumd,Attenuations(n)
      enddo

!  Adjust nstart

      do n = 1, nlayers
         if ( layermask_up(n) .and. attenuations(n-1).ne.zero )  nstart = n
      enddo

!  Enhanced-spherical, fine-layer attenuations

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_up(n) ) then
               do j = 1, nfinedivs(n)
                  sumd = ZERO
                  do k = 1, ntraverse_fine(n,j)
                     sumd = sumd + extinction(k) * sunpaths_fine(n,k,j)
                  enddo
                  if (sumd .lt. cutoff ) Attenuations_fine(n,j) = exp( - sumd )
               enddo
            endif
         enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Regular PS multiplier

      if ( do_regular_ps ) then
         if ( Mu1 .eq. zero ) then
            do n = nlayers, 1, -1
               if ( layermask_up(n) .and. n.le.nstart ) then
                  factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                  multiplier_up(n) = factor1 
               endif
            enddo 
         else
            do n = nlayers, 1, -1
               lostau = deltaus(n) / Mu1
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( layermask_up(n) .and. n.le.nstart ) then
                  factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                  factor2 = (suntau(n) - suntau(n-1))/lostau
                  multiplier_up(n) = factor1 / (factor2 + one)
               endif
!               write(76,*)n,lostau, lostrans_up(n),multiplier_up(n)
            enddo
         endif
      endif

!  Enhanced PS multiplier: special case (nadir viewing)

      if ( do_enhanced_ps .and. doNadir ) then
         do n = nlayers, 1, -1
            kn = extinction(n)
            lostrans_up(n)  = exp ( - deltaus(n))
            if ( layermask_up(n) .and. n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  xjkn = xfine(n,j) * kn
                  func = attenuations_fine(n,j) * exp ( - xjkn )
                  sum = sum + func * wfine(n,j)
               enddo
               multiplier_up(n) = sum * kn
            endif
         enddo
      endif

!  Enhanced PS multiplier: General case

      if ( do_enhanced_ps .and. .not. doNadir ) then
         do n = nlayers, 1, -1
            cot_2 = cota(n-1) ; cot_1 = cota(n)
            kn = extinction(n) ;  ke = raycon * kn
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_up(n) = tran_1
            if ( n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  tran = exp ( - ke * ( cot_2 - cotfine(n,j) ) )
                  func = attenuations_fine(n,j) * csqfine(n,j) * tran
                  sum  = sum + func * wfine(n,j)
               enddo
               multiplier_up(n) = sum * ke 
            endif        
         enddo
      endif

!  Layer sources

      do n = nlayers, 1, -1
         if ( layermask_up(n) .and. n.le.nstart  ) then
            if ( do_sunlight ) then
               do o1 = 1, nstokes
                  sources_up(n,o1) = exactscat_up(n,o1,1) * multiplier_up(n) * fluxvec(1)
               enddo
            else
               do o1 = 1, nstokes
                  sources_up(n,o1) = dot_product(exactscat_up(n,o1,1:ns),fluxvec(1:ns)) * multiplier_up(n)
               enddo
            endif
         endif
      enddo

!  Source function integration
!  ===========================

!  initialize recursion ( FOr Direct Beam, use PI.mu0.R.Atten )

      NC =  0
      CUMSOURCE_UP(NC,:) = zero
      CUMSOURCE_DB       = zero
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Surface term

      M4 = 4.0d0 * MU0

      if ( DO_LAMBERTIAN ) then
         CUMSOURCE_DB(1) = M4 * REFLEC(1,1) * attenuations(nlayers) * fluxvec(1)
      else
         do o1 = 1, nstokes
            sum = zero
            do o2 = 1, nstokes
               sum = sum + reflec(o1,o2) * fluxvec(o2)
            enddo
            CUMSOURCE_DB(o1) = M4 * sum
         enddo
      endif

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            do o1 = 1, nstokes
               CUMSOURCE_DB(O1)     = LOSTRANS_UP(N) * CUMSOURCE_DB(O1)
               CUMSOURCE_UP(NC,O1)  = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1,O1) + SOURCES_UP(N,O1)
            enddo
         ENDDO
         do o1 = 1, nstokes
            STOKES_UP(UTA,O1) = FLUX * CUMSOURCE_UP(NC,O1)
            STOKES_DB(UTA,O1) = FLUX * CUMSOURCE_DB(O1)
         enddo
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Finish

      return
end subroutine SSV_Integral_1_UP

!

subroutine SSV_Integral_1_DN &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels, do_sunlight,      & ! Inputs (dimension)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,                     & ! Inputs (Flags)
     nstokes, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,       & ! Inputs (control)
     extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,                 & ! Inputs (Optical)
     Mu1, GenSpher, Rotations, NCrit, RadCrit, CotCrit,                             & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                           & ! Inputs (Geometry)
     sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                            & ! Inputs (Geometry)
     stokes_dn, cumsource_dn )                                                        ! Outputs

!  Stand alone routine for SS field with Solar sources alone
!    No partials

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

      INTEGER, Intent(in) :: maxlayers
      INTEGER, Intent(in) :: maxfinelayers
      INTEGER, Intent(in) :: maxmoments_input
      INTEGER, Intent(in) :: max_user_levels

!  flags

      LOGICAL, Intent(in) ::  DO_SUNLIGHT
      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Numbers

      INTEGER, Intent(in) ::  NSTOKES
      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  NMOMENTS_INPUT

      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

      REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
      REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
      REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
      REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
      REAL(fpk), Intent(in) :: GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4 )

!  Solar Flux 

      REAL(fpk), Intent(in) :: FLUX, FLUXVEC(4)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      integer  , Intent(in)  :: NCrit
      real(fpk), Intent(in)  :: RadCrit, CotCrit
      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers), radii(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1

!  solar paths 

      integer  , Intent(in)  :: ntraverse  (0:maxlayers)
      real(fpk), Intent(in)  :: sunpaths   (0:maxlayers,maxlayers)
      integer  , Intent(in)  :: ntraverse_fine(maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)

!  Generalized spherical functions.
!    Rotations(1-4)    = C1, S1, C2, S2

      REAL(fpk), Intent(in)  :: GenSpher(0:maxmoments_input, 4)
      REAL(fpk), Intent(in)  :: Rotations(4)

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: stokes_dn     ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: cumsource_dn  ( 0:maxlayers,     4 )

!  LOCAL
!  -----

!  Attenuations

      real(fpk)  :: suntau            (0:maxlayers)
      real(fpk)  :: attenuations      (0:maxlayers)
      real(fpk)  :: attenuations_fine (maxlayers,maxfinelayers)

!  Scattering

      real(fpk)  :: tms (maxlayers)
      real(fpk)  :: exactscat_dn (maxlayers,4,4)

!  Source function integration results

      real(fpk)  :: sources_dn       ( maxlayers, 4 )
      real(fpk)  :: lostrans_dn      ( maxlayers )
      real(fpk)  :: multiplier_dn    ( maxlayers )

!  Help

      integer    :: n, ns, uta, nstart, nc, nut, nut_prev, j, k, L, O1, O2
      logical    :: layermask_dn(maxlayers)
      real(fpk)  :: sumd, help, sum, tran_1, tran, func, kn, ke, xjkn, trand
      real(fpk)  :: cot_1, cot_2, factor1, factor2, factor3, lostau, rdiff, cot_c
      real(fpk)  :: fmat(6), sum23, dif23, help3c1, help3s1, help4c1, help4s1

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      CUMSOURCE_DN = zero ; STOKES_DN  = zero
      lostrans_dn  = zero ; sources_dn = zero ; exactscat_dn = zero ; multiplier_dn = zero

!  Bookkeeping

      ns = nstokes ; fmat = zero
      NUT = USER_LEVELS(N_USER_LEVELS) + 1
      IF ( NUT > NLAYERS ) NUT = NLAYERS
      LAYERMASK_DN = .false.
      LAYERMASK_DN(1:NUT) = .true.

!  TMS factors

      do n = 1, nlayers
         if ( do_deltam_scaling ) then
            help = one - truncfac(n) * omega(n)
            tms(n) = omega(n) / help
         else
            tms(n) = omega(n)
         endif
      enddo

!  Scattering functions
!  --------------------

!  Scalar only

      if ( nstokes .eq. 1 ) then
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
              sum = zero
               do L = 0, nmoments_input
                  sum = sum + GenSpher(L,1) * Greekmat(n,L,1,1)
               enddo
               exactscat_dn(n,1,1) = sum * tms(n)
            endif
         enddo
      endif

!  Vector with Sunlight

      if ( nstokes .gt. 1 .and. do_sunlight ) then
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
               fmat(1:2) = zero
               do L = 0, nmoments_input
                  fmat(1) = fmat(1) + Genspher(L,1) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2) * greekmat(n,L,1,2)
               enddo
               exactscat_dn(n,1,1) = + fmat(1)
               exactscat_dn(n,2,1) = - fmat(2) * Rotations(3)
               exactscat_dn(n,3,1) = + fmat(2) * Rotations(4)
               exactscat_dn(n,1:ns,1) = tms(n) * exactscat_dn(n,1:ns,1) 
            endif
         enddo
      endif

!  Vector General case, USE FULL 4X4 MATRIX
!    CODE INTRODUCED BUT NOT TESTED, 05 OCTOBER 2010

      if ( nstokes .gt. 1 .and. .not. do_sunlight ) then
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
               fmat = zero
               do L = 0, nmoments_input
                  fmat(1) = fmat(1) + Genspher(L,1) * greekmat(n,L,1,1)
                  fmat(2) = fmat(2) + Genspher(L,2) * greekmat(n,L,1,2)
                  sum23 = greekmat(n,L,2,2) + greekmat(n,L,3,3)
                  dif23 = greekmat(n,L,2,2) - greekmat(n,L,3,3)
                  fmat(3) = fmat(3) + Genspher(L,3) * sum23
                  fmat(4) = fmat(4) + Genspher(L,4) * dif23
               enddo
               fmat(3) = ( fmat(3) + fmat(4) ) * 0.5_fpk
               fmat(4) = ( fmat(3) - fmat(4) )
               if ( nstokes.eq.4) then
                  do L = 0, nmoments_input
                     fmat(5) = fmat(5) + Genspher(L,2) * greekmat(n,L,3,4)
                     fmat(6) = fmat(6) + Genspher(L,1) * greekmat(n,L,4,4)
                  enddo
               endif
               help3c1 = fmat(3) * Rotations(1)
               help3s1 = fmat(3) * Rotations(2)
               help4c1 = fmat(4) * Rotations(1)
               help4s1 = fmat(4) * Rotations(2)
               exactscat_dn(n,1,1) = + fmat(1)
               exactscat_dn(n,2,1) = - fmat(2) * Rotations(3)
               exactscat_dn(n,3,1) = + fmat(2) * Rotations(4)
               exactscat_dn(n,1,2) = + fmat(2) * Rotations(1)
               exactscat_dn(n,1,3) = - fmat(2) * Rotations(2)
               exactscat_dn(n,2,2) = + help3c1 * Rotations(3) - help4s1 * Rotations(4)
               exactscat_dn(n,2,3) = - help3s1 * Rotations(3) - help4c1 * Rotations(4)
               exactscat_dn(n,3,2) = + help3c1 * Rotations(4) + help4s1 * Rotations(3)
               exactscat_dn(n,3,3) = - help3s1 * Rotations(4) + help4c1 * Rotations(3)
               if ( nstokes .eq. 4 ) then
                  exactscat_dn(n,2,4) = - fmat(5) * Rotations(4) 
                  exactscat_dn(n,4,2) = - fmat(5) * Rotations(2) 
                  exactscat_dn(n,3,4) = + fmat(5) * Rotations(3) 
                  exactscat_dn(n,4,3) = - fmat(5) * Rotations(1) 
                  exactscat_dn(n,4,4) = + fmat(6)
               endif
               exactscat_dn(n,1:ns,1:ns) = tms(n)*exactscat_dn(n,1:ns,1:ns)
            endif
         enddo
      endif

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations = ZERO ; Attenuations_fine = ZERO
      nstart = nlayers ; if (Ncrit.ne.0) nstart = nCrit

!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         sumd = zero
         do k = 1, ntraverse(n)
            sumd = sumd + extinction(k) * sunpaths(n,k)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
      enddo

!  Adjust nstart

      do n = 1, nlayers
         if ( layermask_dn(n) .and. attenuations(n-1).ne.zero )  nstart = n
      enddo

!  Enhanced-spherical, fine-layer attenuations

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_dn(n) ) then
               do j = 1, nfinedivs(n)
                  sumd = zero
                  do k = 1, ntraverse_fine(n,j)
                     sumd = sumd + extinction(k) * sunpaths_fine(n,k,j)
                  enddo
                  if (sumd .lt. cutoff ) Attenuations_fine(n,j) = exp( - sumd )
               enddo
            endif
         enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Regular PS multipliers and LOSTRANS

      if ( do_regular_ps ) then
         if( Mu1 .eq. zero ) then
            do n = 1, nlayers
               if ( layermask_dn(n).and. n.le.nstart ) then
                  if ( attenuations(n-1).ne.zero ) then
                     factor1 = Attenuations(n)
                     multiplier_dn(n) = factor1
                  endif
               endif
            enddo
         else
            do n = 1, nlayers
               lostau = deltaus(n) / Mu1
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               if ( layermask_dn(n).and. n.le.nstart ) then
                  factor1 = lostrans_dn(n) * Attenuations(n-1) - Attenuations(n)
                  factor2 = (suntau(n) - suntau(n-1))/lostau
                  factor3 = factor2 - one
                  multiplier_dn(n) = factor1 / factor3
               endif
            enddo
         endif
      endif
!  Enhanced PS: special case (nadir viewing)

      if ( do_enhanced_ps .and. doNadir ) then
         do n = nlayers, 1, -1
            kn = extinction(n)
            lostrans_dn(n)  = exp ( - deltaus(n))
            rdiff = radii(n-1) - radii(n) ; if ( n.eq.NCrit) rdiff = radii(n-1) - RadCrit
            trand = one ; if ( n.eq.NCrit) trand = exp ( -kn * (RadCrit -radii(n) ) )
            if ( layermask_dn(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  xjkn = ( rdiff - xfine(n,j) ) * kn
                  func = attenuations_fine(n,j) * exp ( - xjkn )
                  sum = sum + func * wfine(n,j)
               enddo
               multiplier_dn(n) = sum * kn
            endif
         enddo
      endif

!  Enhanced PS: General case

      if ( do_enhanced_ps .and. .not. doNadir ) then
         do n = nlayers, 1, -1
            cot_2 = cota(n-1) ; cot_1 = cota(n)
            cot_c = cot_1  ; if ( n.eq.NCrit ) cot_c = CotCrit
            kn = extinction(n) ;  ke = raycon * kn
            trand = one  ; if ( n.eq.NCrit ) trand = exp ( - ke * ( CotCrit - cot_1 ) )
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_dn(n) = tran_1
            if ( layermask_dn(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  tran = exp ( - ke * ( cotfine(n,j) - cot_c ) )   !  Down
!                 tran = exp ( - ke * ( cot_2 - cotfine(n,j) ) )   !  Up
                  func = Attenuations_fine(n,j) * csqfine(n,j) * tran
                  sum  = sum + func * wfine(n,j)
               enddo
               multiplier_dn(n) = sum * ke * trand
            endif        
         enddo
      endif

!  Layer sources

      do n = nlayers, 1, -1
         if ( layermask_dn(n) .and. n.le.nstart ) then
            do o1 = 1, nstokes
               sum = zero
               do o2 = 1, nstokes
                  sum = sum + exactscat_dn(n,o1,o2) * fluxvec(o2)
               enddo
               sources_dn(n,o1) = sum * multiplier_dn(n)
            enddo
         endif
      enddo

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC,:) = zero
      NSTART = 1 ; NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working Downn from NSTART to NUT
!     Check for dndating the recursion

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            do o1 = 1, nstokes
               CUMSOURCE_DN(NC,O1)  = LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1,O1) + SOURCES_DN(N,O1)
            enddo
         ENDDO
         do o1 = 1, nstokes
            STOKES_DN(UTA,O1) = FLUX * CUMSOURCE_DN(NC,O1)
         enddo
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  Finish

      return
end subroutine SSV_Integral_1_DN

!

subroutine SSV_Integral_1_UPDN   &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels,           & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_deltam_scaling, do_sunlight,            & ! Inputs (Flags)
     do_lambertian, do_regular_ps, do_enhanced_ps, doNadir, nstokes,        & ! Inputs (Flags)
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,        & ! Inputs (control output)
     reflec, extinction, deltaus, omega, truncfac, Greekmat, flux, fluxvec, & ! Inputs (Optical)
     Mu0, Mu1, GenSpher_up, GenSpher_dn, Rotations_up, Rotations_dn, NCrit, & ! Inputs (Geometry)
     RadCrit, CotCrit, xfine, wfine, csqfine, cotfine, Raycon, radii, cota, & ! Inputs (Geometry)
     sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up,        & ! Inputs (Geometry)
     sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,        & ! Inputs (Geometry)
     stokes_up, stokes_db, cumsource_up, stokes_dn, cumsource_dn )            ! Outputs

!  Stand alone routine for SS field with Solar sources alone
!    No partials

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

      INTEGER, Intent(in) :: maxlayers
      INTEGER, Intent(in) :: maxfinelayers
      INTEGER, Intent(in) :: maxmoments_input
      INTEGER, Intent(in) :: max_user_levels

!  flags

      LOGICAL, Intent(in) ::  DO_SUNLIGHT
      LOGICAL, Intent(in) ::  DO_UPWELLING
      LOGICAL, Intent(in) ::  DO_DNWELLING
      LOGICAL, Intent(in) ::  DO_LAMBERTIAN
      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Numbers

      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  NMOMENTS_INPUT, NSTOKES

      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

      REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
      REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
      REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
      REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
      REAL(fpk), Intent(in) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4 )

!  Solar Flux and Surface reflectivity (Could be the albedo)

      REAL(fpk), Intent(in) :: REFLEC(4,4), FLUX, FLUXVEC(4)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      integer  , Intent(in)  :: NCrit
      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers), radii(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1, Mu0, RadCrit, CotCrit

!  solar paths 

      integer  , Intent(in)  :: ntraverse_up  (0:maxlayers)
      real(fpk), Intent(in)  :: sunpaths_up   (0:maxlayers,maxlayers)
      integer  , Intent(in)  :: ntraverse_fine_up(maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: sunpaths_fine_up (maxlayers,maxlayers,maxfinelayers)

      integer  , Intent(in)  :: ntraverse_dn  (0:maxlayers)
      real(fpk), Intent(in)  :: sunpaths_dn   (0:maxlayers,maxlayers)
      integer  , Intent(in)  :: ntraverse_fine_dn(maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: sunpaths_fine_dn (maxlayers,maxlayers,maxfinelayers)

!  Generalized spherical functions.
!    Rotations(1-4)    = C1, S1, C2, S2

      REAL(fpk), Intent(in)  :: GenSpher_up(0:maxmoments_input, 4)
      REAL(fpk), Intent(in)  :: GenSpher_dn(0:maxmoments_input, 4)
      REAL(fpk), Intent(in)  :: Rotations_up(4)
      REAL(fpk), Intent(in)  :: Rotations_dn(4)

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: stokes_up     ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: stokes_db     ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: cumsource_up  ( 0:maxlayers, 4 )
      real(fpk), Intent(Out)  :: stokes_dn     ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: cumsource_dn  ( 0:maxlayers, 4 )

!  Upwelling
!  ---------

   if ( do_upwelling ) then
       call SSV_Integral_1_UP &
    ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels, do_sunlight,     & ! Inputs (dimension)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, do_lambertian, doNadir,      & ! Inputs (Flags)
     nstokes, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,       & ! Inputs (control)
     reflec, extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,         & ! Inputs (Optical)
     Mu0, Mu1, GenSpher_up, Rotations_up, NCrit, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
     Raycon, cota, sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up,  & ! Inputs (Geometry)
     stokes_up, stokes_db, cumsource_up )                                             ! Outputs
   endif

   if ( do_dnwelling ) then
       call SSV_Integral_1_DN &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels, do_sunlight,      & ! Inputs (dimension)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,                     & ! Inputs (Flags)
     nstokes, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,       & ! Inputs (control)
     extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,                 & ! Inputs (Optical)
     Mu1, GenSpher_dn, Rotations_dn, NCrit, RadCrit, CotCrit,                       & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                           & ! Inputs (Geometry)
     sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,                & ! Inputs (Geometry)
     stokes_dn, cumsource_dn )                                                        ! Outputs
   endif

!  Finish

   return
end subroutine SSV_Integral_1_UPDN



subroutine DTEV_Integral_1_UP &
   ( maxlayers, maxfinelayers, max_user_levels, Do_Thermset, & ! Inputs (dimensioning)
     do_dmscaling, do_regular_ps, do_enhanced_ps, doNadir,   & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,         & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                      & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                   & ! Inputs (Optical)
     Mu1, Raycon, cota, xfine, wfine, csqfine, cotfine,      & ! Inputs (Geometry)
     intensity_date_up, intensity_dste, cumsource_up, tcom1 )  ! Outputs

!  Stand alone routine for DTE field with Thermal Emission sources
!    No partials

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

      INTEGER, Intent(in) :: maxlayers
      INTEGER, Intent(in) :: maxfinelayers
      INTEGER, Intent(in) :: max_user_levels

!  Thermal setup flag (for TCOM1)

      LOGICAL, Intent(inout) ::  Do_Thermset

!  flags

      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DMSCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Numbers

      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere extinction and deltaus

      REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
      REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
      REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
      REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric BB functions and Surface BB and emissivity

      REAL(fpk), Intent(in) :: SURFBB, USER_EMISSIVITY
      REAL(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!   Ray constant, Secant of LOS angle (for Regular PS only)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      real(fpk), Intent(in)  :: Mu1
      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers)

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_date_up ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: intensity_dste    ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: cumsource_up      ( 0:maxlayers )

!  Thermal setup

      real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)

!  LOCAL
!  -----

!  Source function integration results

      real(fpk)  :: sources_up       ( maxlayers )
      real(fpk)  :: lostrans_up      ( maxlayers )

!  Help

      integer    :: n, uta, nstart, nc, nut, nut_prev, j
      logical    :: layermask_up(maxlayers)
      real(fpk)  :: help, sum, tran_1, tran, func, kn, ke, xjkn, cot_1, cot_2
      real(fpk)  :: cumsource_dste, t_mult_up(0:2), thermcoeffs(2), tms

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      CUMSOURCE_UP = zero ; INTENSITY_DATE_up = zero ; INTENSITY_DSTE = zero
      lostrans_up = zero  ; sources_up = zero

!  Bookkeeping

      NUT = USER_LEVELS(1) + 1
      LAYERMASK_UP = .false.
      LAYERMASK_UP(NUT:NLAYERS) = .true.

!  Thermal setup factors
!     TMS, Initial set of thermal coefficients and TCOM1 variable

      if ( do_Thermset ) then
         tcom1 = zero
         do n = 1, nlayers
            tms = one - omega(n) 
            if ( do_dmscaling ) then
               help = one - truncfac(n) * omega(n)
               tms = tms / help
            endif
            thermcoeffs(1)  = bb_input(n-1)
            thermcoeffs(2)  = (bb_input(n)-bb_input(n-1)) / deltaus(n)
            tcom1(n,1) = thermcoeffs(1) * tms
            tcom1(n,2) = thermcoeffs(2) * tms
         ENDDO
         do_Thermset = .false.
      endif

!  Layer integrated Solar sources
!  ==============================

!  Regular PS

      if ( do_regular_ps ) then
        if ( Mu1 .eq. zero ) then
          DO n = 1, nlayers
            t_mult_up(2) = tcom1(n,2)
            t_mult_up(1) = tcom1(n,1)
            sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
            t_mult_up(0) = - sum
            sources_up(n) = t_mult_up(1)
          enddo
        else
          DO n = 1, nlayers
            lostrans_up(n) = exp ( - deltaus(n) / Mu1 )
            t_mult_up(2) = tcom1(n,2)
            t_mult_up(1) = tcom1(n,1) + t_mult_up(2) * Mu1
            sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
            t_mult_up(0) = - sum
            sources_up(n) = t_mult_up(0) * lostrans_up(n)  + t_mult_up(1)
 !               if (n.gt.67)write(*,*)n,sources_up(n)
          enddo
        endif
      endif

!  Enhanced PS: special case (nadir viewing)

      if ( do_enhanced_ps .and. doNadir ) then
         do n = nlayers, 1, -1
            kn = extinction(n)
            lostrans_up(n)  = exp ( - deltaus(n))
            sum = zero
            do j = 1, nfinedivs(n)
               xjkn = xfine(n,j) * kn
               func = ( tcom1(n,1) + xjkn * tcom1(n,2) ) * exp ( - xjkn )
               sum = sum + func * wfine(n,j)
            enddo
            sources_up(n) = sum * kn
         enddo
      endif

!  Enhanced PS: General case

      if ( do_enhanced_ps .and. .not. doNadir ) then
         do n = nlayers, 1, -1
            cot_2 = cota(n-1) ; cot_1 = cota(n)
            kn = extinction(n) ;  ke = raycon * kn
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_up(n) = tran_1
            sum = zero
            do j = 1, nfinedivs(n)
               xjkn = xfine(n,j) * kn
               tran = exp ( - ke * ( cot_2 - cotfine(n,j) ) )
               func = ( tcom1(n,1) + xjkn * tcom1(n,2) ) * csqfine(n,j) * tran
               sum  = sum + func * wfine(n,j)
            enddo
            sources_up(n) = sum * ke 
         enddo
      endif

!  Source function integration
!  ===========================

!  start recursion ( For DSTE term, Use surface emissivity )

      NC =  0
      CUMSOURCE_UP(NC) = zero
      CUMSOURCE_DSTE   = SURFBB * USER_EMISSIVITY
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_DSTE   = LOSTRANS_UP(N) * CUMSOURCE_DSTE
            CUMSOURCE_UP(NC) = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1) + SOURCES_UP(N)
!            write(45,*)n,sources_up(n)
        ENDDO
         INTENSITY_DATE_UP(UTA,1) = CUMSOURCE_UP(NC)
         INTENSITY_DSTE(UTA,1)    = CUMSOURCE_DSTE
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Finish

      return
end subroutine DTEV_Integral_1_UP

!

subroutine DTEV_Integral_1_DN &
   ( maxlayers, maxfinelayers, max_user_levels, do_Thermset,      & ! Inputs (dimensioning)
     do_dmscaling, do_regular_ps, do_enhanced_ps, doNadir,        & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,              & ! Inputs (control output)
     BB_input, extinction, deltaus, omega, truncfac,              & ! Inputs (Optical)
     Mu1, Raycon, radii, cota, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
     intensity_date_dn, cumsource_dn, tcom1 )                       ! Outputs

!  Stand alone routine for DTE field with Thermal Emission sources
!    No partials

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

      INTEGER, Intent(in) :: maxlayers
      INTEGER, Intent(in) :: maxfinelayers
      INTEGER, Intent(in) :: max_user_levels

!  Thermal setup flag (for TCOM1)

      LOGICAL, Intent(inout) ::  Do_Thermset

!  flags

      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DMSCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Numbers

      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

      REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
      REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
      REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
      REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric thermal BB functions

      REAL(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers), radii(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_date_dn( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: cumsource_dn     ( 0:maxlayers )

!  Thermal setup

      real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)

!  LOCAL
!  -----

!  Source function integration results

      real(fpk)  :: sources_dn       ( maxlayers )
      real(fpk)  :: lostrans_dn      ( maxlayers )

!  Help

      integer    :: n, uta, nstart, nc, nut, nut_prev, j
      logical    :: layermask_dn(maxlayers)
      real(fpk)  :: help, sum, tran_1, tran, func, kn, ke, xjkn
      real(fpk)  :: cot_1, cot_2, rdiff
      real(fpk)  :: t_mult_dn(0:2), thermcoeffs(2), tms

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      CUMSOURCE_DN = zero ; INTENSITY_DATE_DN = zero
      lostrans_dn  = zero ; sources_dn   = zero

!  Bookkeeping

      NUT = USER_LEVELS(N_USER_LEVELS) + 1
      LAYERMASK_DN = .false.
      LAYERMASK_DN(1:NUT) = .true.

!  Thermal setup factors
!     TMS, Initial set of thermal coefficients and TCOM1 variable

      if ( do_Thermset ) then
         tcom1 = zero
         do n = 1, nlayers
            tms = one - omega(n) 
            if ( do_dmscaling ) then
               help = one - truncfac(n) * omega(n)
               tms = tms / help
            endif
            thermcoeffs(1)  = bb_input(n-1)
            thermcoeffs(2)  = (bb_input(n)-bb_input(n-1)) / deltaus(n)
            tcom1(n,1) = thermcoeffs(1) * tms
            tcom1(n,2) = thermcoeffs(2) * tms
         ENDDO
         do_Thermset = .false.
      endif

!  Layer integrated Thermal sources
!  ================================

!  Regular PS

      if ( do_regular_ps ) then
         if ( Mu1 .eq. zero ) then
            DO n = 1, nlayers
               t_mult_dn(2)   = tcom1(n,2)
               t_mult_dn(1)   = tcom1(n,1)
               t_mult_dn(0)   = - t_mult_dn(1)
               sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
               sources_dn(n)  = sum
            END DO
         else
            DO n = 1, nlayers
               lostrans_dn(n) = exp ( - deltaus(n) / Mu1 )
               t_mult_dn(2)   = tcom1(n,2)
               t_mult_dn(1)   = tcom1(n,1) - t_mult_dn(2) * Mu1
               t_mult_dn(0)   = - t_mult_dn(1)
               sources_dn(n)  = t_mult_dn(0) * lostrans_dn(n)
               sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
               sources_dn(n)  = sources_dn(n) + sum
            END DO
         endif
      endif

!  Enhanced PS: special case (nadir viewing)

      if ( do_enhanced_ps .and. doNadir ) then
         do n = nlayers, 1, -1
            kn = extinction(n)
            lostrans_dn(n)  = exp ( - deltaus(n))
            rdiff = radii(n-1) - radii(n)
            do j = 1, nfinedivs(n)
               xjkn = ( rdiff - xfine(n,j) ) * kn
               func = ( tcom1(n,1) + xjkn * tcom1(n,2) ) * exp ( - xjkn )
               sum = sum + func * wfine(n,j)
            enddo
            sources_dn(n) = sum * kn
         enddo
      endif

!  Enhanced PS: General case

      if ( do_enhanced_ps .and. .not. doNadir ) then
         do n = nlayers, 1, -1
            cot_2 = cota(n-1) ; cot_1 = cota(n)
            kn = extinction(n) ;  ke = raycon * kn
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_dn(n) = tran_1
            sum = zero
            do j = 1, nfinedivs(n)
               tran = exp ( - ke * ( cotfine(n,j) - cot_1 ) )   !  Down
               func = ( tcom1(n,1) + xfine(n,j) * tcom1(n,2) ) * csqfine(n,j) * tran
               sum  = sum + func * wfine(n,j)
            enddo
            sources_dn(n) = sum * ke
         enddo
      endif

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working Downn from NSTART to NUT
!     Check for dndating the recursion

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            CUMSOURCE_DN(NC) = SOURCES_DN(N) + LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1)
         ENDDO
         INTENSITY_DATE_DN(UTA,1) = CUMSOURCE_DN(NC)
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  Finish

      return
end subroutine DTEV_Integral_1_DN



subroutine DTEV_Integral_1_UPDN   &
   ( maxlayers, maxfinelayers, max_user_levels, do_Thermset,  & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_deltam_scaling,           & ! Inputs (Flags)
     do_regular_ps, do_enhanced_ps, doNadir,                  & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,          & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                       & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                    & ! Inputs (Optical)
     Mu1, Raycon, radii, cota, xfine, wfine, csqfine, cotfine,& ! Inputs (Geometry)
     intensity_date_up, intensity_dste, cumsource_up,         & ! Outputs
     intensity_date_dn, cumsource_dn, tcom1 )                   ! Outputs

!  Stand alone routine for SS field with Solar sources alone
!    No partials

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

      INTEGER, Intent(in) :: maxlayers
      INTEGER, Intent(in) :: maxfinelayers
      INTEGER, Intent(in) :: max_user_levels

!  Thermal setup flag (for TCOM1)

      LOGICAL, Intent(inout) ::  Do_Thermset

!  flags

      LOGICAL, Intent(in) ::  DO_UPWELLING
      LOGICAL, Intent(in) ::  DO_DNWELLING
      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Numbers

      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

      REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
      REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
      REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
      REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric BB functions and Surface BB and emissivity

      REAL(fpk), Intent(in) :: SURFBB, USER_EMISSIVITY
      REAL(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers), radii(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_date_up     ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: intensity_dste        ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: cumsource_up          ( 0:maxlayers )
      real(fpk), Intent(Out)  :: intensity_date_dn     ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: cumsource_dn          ( 0:maxlayers )

!  Thermal setup

      real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)

!  Upwelling
!  ---------

   if ( do_upwelling ) then
      call DTEV_Integral_1_UP &
   ( maxlayers, maxfinelayers, max_user_levels, Do_Thermset,   & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,& ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,           & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                        & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                     & ! Inputs (Optical)
     Mu1, Raycon, cota, xfine, wfine, csqfine, cotfine,        & ! Inputs (Geometry)
     intensity_date_up, intensity_dste, cumsource_up, tcom1 )    ! Outputs
   endif

   if ( do_dnwelling ) then
       call DTEV_Integral_1_DN &
   ( maxlayers, maxfinelayers, max_user_levels, do_Thermset,      & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,   & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,              & ! Inputs (control output)
     BB_input, extinction, deltaus, omega, truncfac,              & ! Inputs (Optical)
     Mu1, Raycon, radii, cota, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
     intensity_date_dn, cumsource_dn, tcom1 )                       ! Outputs
   endif

!  Finish

   return
end subroutine DTEV_Integral_1_UPDN

!  End module

end module FO_Vector_RTCalcs_m
