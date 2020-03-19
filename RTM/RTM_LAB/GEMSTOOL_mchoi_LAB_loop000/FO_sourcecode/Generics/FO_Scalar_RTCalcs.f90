! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! #  This Version :   3.6 F90                               #
! #  Release Date :   March 2012                            #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       New Direct-Thermal & Single-scatter codes (3.6)   #
! #       New Fast Geometry + Precalculation choice (3.6)   #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

module FO_Scalar_RTCalcs_m

!  For a given wavelength, this routine will calculate Corrected Intensities:

!     (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.
!     (2) For the Atmospheric and Surface Direct Thermal Emission (DTE) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel or incoming solar pseudo-spherical)

!  This is Versions 1-2, without Partials. Code is stand alone with no dependencies.
!    Version 1, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 2, 01 June     2012, R. Spurr, RT Solutions Inc.

!  Option (1). For Solar sources, the subroutines are
!       SS_Integral_1_UP   (Upwelling only)
!       SS_Integral_1_DN   (Downwelling only)
!       SS_Integral_1_UPDN (Upwelling and Downwelling)

!  Option (2). For Thermal Emission sources, the subroutines are
!       DTE_Integral_1_UP   (Upwelling only)
!       DTE_Integral_1_DN   (Downwelling only)
!       DTE_Integral_1_UPDN (Upwelling and Downwelling)

!  All subroutines public

public

contains


subroutine SS_Integral_1_UP &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels,           & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,             & ! Inputs (Flags)
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,        & ! Inputs (control output)
     reflec, extinction, deltaus, omega, truncfac, phasmoms, flux,          & ! Inputs (Optical)
     Mu0, Mu1, LegPoly_up, NCrit, xfine, wfine, csqfine, cotfine,           & ! Inputs (Geometry)
     Raycon, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,      & ! Inputs (Geometry)
     intensity_up, intensity_db, cumsource_up )                               ! Outputs

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of radiance. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012

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

      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Numbers

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
      REAL(fpk), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )

!  Solar Flux and Surface reflectivity (Could be the albedo)

      REAL(fpk), Intent(in) :: REFLEC, FLUX

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

!  Legendres

      REAL(fpk), Intent(in)  :: LegPoly_up(0:maxmoments_input)

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels )
      real(fpk), Intent(Out)  :: intensity_db     ( max_user_levels )
      real(fpk), Intent(Out)  :: cumsource_up     ( 0:maxlayers )

!  LOCAL
!  -----

!  Attenuations

      real(fpk)  :: attenuations      (0:maxlayers)
      real(fpk)  :: attenuations_fine (maxlayers,maxfinelayers)

!  Solutions

      real(fpk)  :: Solutions_fine (maxlayers,maxfinelayers)
      real(fpk)  :: Solutions (0:maxlayers)

!  Scattering

      real(fpk)  :: tms (maxlayers)
      real(fpk)  :: exactscat_up (maxlayers)

!  Source function integration results

      real(fpk)  :: sources_up       ( maxlayers )
      real(fpk)  :: lostrans_up      ( maxlayers )

!  Average secant type arrays (REGULAR-PS only)

      real(fpk)  :: factor1       ( maxlayers )
      real(fpk)  :: factor2       ( maxlayers )

!  Help

      integer    :: n, uta, nstart, nc, nut, nut_prev, j, k, L
      logical    :: layermask_up(maxlayers)
      real(fpk)  :: sumd, help, sum, tran_1, tran, func, kn, ke, xjkn
      real(fpk)  :: cot_1, cot_2, multiplier, suntau(0:maxlayers), lostau
      real(fpk)  :: cumsource_db

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      CUMSOURCE_UP = zero ; INTENSITY_UP = zero ; INTENSITY_DB = zero
      lostrans_up = zero  ; sources_up   = zero

!  Bookkeeping

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

!  Scattering function

      do n = 1, nlayers
         if ( layermask_up(n) ) then
            sum = zero
            do L = 0, nmoments_input
              sum = sum + LegPoly_Up(L) * phasmoms(n,L)
            enddo
            exactscat_up(n) = sum * tms(n)
         endif
      enddo

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations = ZERO ; Attenuations_fine = ZERO
      Solutions    = zero ; Solutions_fine    = zero
      nstart = nlayers ; if (Ncrit.ne.0) nstart = nCrit

!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         sumd = ZERO
         do k = 1, ntraverse(n)
            sumd = sumd + extinction(k) * sunpaths(n,k)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
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
                  Solutions_fine(n,j) = exactscat_up(n) * Attenuations_fine(n,j)
!            if(n.eq.Ncrit.and.j.eq.2)write(*,*)'FD',n, Solutions_fine(n,j)
               enddo
            endif
         enddo
      endif

!  Regular-PS (Average secant formulation)
!    Special treatment for the horizonal case --> Factor2 = 0, lostrans = 0

      if ( do_regular_ps ) then
         factor1 = zero ; factor2 = zero
         if ( Mu1 .eq. zero ) then
            do n = 1, nlayers
               Solutions(n) = exactscat_up(n) * Attenuations(n-1)
               if ( layermask_up(n) ) then
                  if ( attenuations(n-1).ne.zero ) then
                     factor1(n) = Attenuations(n)/Attenuations(n-1)
                     nstart = n
                  endif
               endif
            enddo
         else
            do n = 1, nlayers
               lostau = deltaus(n) / Mu1
               Solutions(n) = exactscat_up(n) * Attenuations(n-1)
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( layermask_up(n) ) then
                  if ( attenuations(n-1).ne.zero ) then
                     factor1(n) = Attenuations(n)/Attenuations(n-1)
                     factor2(n) = (suntau(n) - suntau(n-1))/lostau
                     nstart = n
                  endif
               endif
            enddo
         endif
      endif

!  Layer integrated Solar sources
!  ==============================

!  Regular PS, all cases (Multiplier = 1.0 for Horizontal case)

      if ( do_regular_ps ) then
         do n = nlayers, 1, -1
            if ( n.le.nstart ) then
               multiplier = ( one - Factor1(n)*lostrans_up(n) ) / (factor2(n) + one)
            else
               multiplier = zero
            endif  
            sources_up(n) = solutions(n) * multiplier
         enddo
      endif

!  Enhanced PS: special case (nadir viewing)

      if ( do_enhanced_ps .and. doNadir ) then
         do n = nlayers, 1, -1
            kn = extinction(n)
            lostrans_up(n)  = exp ( - deltaus(n))
            if ( n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  xjkn = xfine(n,j) * kn
                  func = solutions_fine(n,j) * exp ( - xjkn )
                  sum = sum + func * wfine(n,j)
               enddo
               sources_up(n) = sum * kn
            endif
         enddo
      endif

!  Enhanced PS: General case

      if ( do_enhanced_ps .and. .not. doNadir ) then
         do n = nlayers, 1, -1
            cot_2 = cota(n-1) ; cot_1 = cota(n)
            kn = extinction(n) ;  ke = raycon * kn
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_up(n) = tran_1
!            if(n.eq.NCrit)write(*,*)n,lostrans_up(n)
            if ( n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  tran = exp ( - ke * ( cot_2 - cotfine(n,j) ) )
                  func = solutions_fine(n,j) * csqfine(n,j) * tran
                  sum  = sum + func * wfine(n,j)
               enddo
               sources_up(n) = sum * ke 
            endif        
         enddo
      endif

!  Source function integration
!  ===========================

!  start recursion ( FOr Direct Beam, use PI.mu0.R.Atten )

      NC =  0
      CUMSOURCE_UP(NC) = zero
      CUMSOURCE_DB     = 4.0_fpk * Mu0 * REFLEC * attenuations(nlayers)
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
            CUMSOURCE_DB     = LOSTRANS_UP(N) * CUMSOURCE_DB
            CUMSOURCE_UP(NC) = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1) + SOURCES_UP(N)
         ENDDO
         INTENSITY_UP(UTA) = FLUX * CUMSOURCE_UP(NC)
         INTENSITY_DB(UTA) = FLUX * CUMSOURCE_DB
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Finish

      return
end subroutine SS_Integral_1_UP

!

subroutine SS_Integral_1_DN &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels,              & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,                & ! Inputs (Flags)
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,           & ! Inputs (control output)
     extinction, deltaus, omega, truncfac, phasmoms, flux,                     & ! Inputs (Optical)
     Mu1, LegPoly_dn, NCrit, RadCrit, CotCrit, xfine, wfine, csqfine, cotfine, & ! Inputs (Geometry)
     Raycon, radii, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,  & ! Inputs (Geometry)
     intensity_dn, cumsource_dn )                                                ! Outputs

!  Stand-alone routine for Downwelling Solar-beam Single-scatter (SS)
!    computation of radiance. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012

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

      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  doNadir

!  Numbers

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
      REAL(fpk), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )

!  Solar Flux

      REAL(fpk), Intent(in) ::  FLUX

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

!  Legendres

      REAL(fpk), Intent(in)  :: LegPoly_dn(0:maxmoments_input)

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_dn     ( max_user_levels )
      real(fpk), Intent(Out)  :: cumsource_dn     ( 0:maxlayers )

!  LOCAL
!  -----

!  Attenuations

      real(fpk)  :: attenuations      (0:maxlayers)
      real(fpk)  :: attenuations_fine (maxlayers,maxfinelayers)

!  Solutions

      real(fpk)  :: Solutions_fine (maxlayers,maxfinelayers)
      real(fpk)  :: Solutions (0:maxlayers)

!  Scattering

      real(fpk)  :: tms (maxlayers)
      real(fpk)  :: exactscat_dn (maxlayers)

!  Source function integration results

      real(fpk)  :: sources_dn       ( maxlayers )
      real(fpk)  :: lostrans_dn      ( maxlayers )

!  Average secant type arrays (REGULAR-PS only)

      real(fpk)  :: factor1       ( maxlayers )
      real(fpk)  :: factor2       ( maxlayers )

!  Help

      integer    :: n, uta, nstart, nc, nut, nut_prev, j, k, L
      logical    :: layermask_dn(maxlayers)
      real(fpk)  :: sumd, help, sum, tran_1, tran, func, kn, ke, xjkn, trand
      real(fpk)  :: cot_1, cot_2, multiplier, suntau(0:maxlayers), lostau, rdiff, cot_c

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      CUMSOURCE_DN = zero ; INTENSITY_DN = zero
      lostrans_dn  = zero ; sources_dn   = zero

!  Bookkeeping

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

!  Scattering function

      do n = 1, nlayers
         if ( layermask_dn(n) ) then
            sum = zero
            do L = 0, nmoments_input
              sum = sum + LegPoly_Dn(L) * phasmoms(n,L)
            enddo
            exactscat_dn(n) = sum * tms(n)
         endif
      enddo

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations = ZERO ; Attenuations_fine = ZERO
      Solutions    = zero ; Solutions_fine    = zero
      nstart = nlayers ; if (Ncrit.ne.0) nstart = nCrit

!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         sumd = ZERO
         do k = 1, ntraverse(n)
            sumd = sumd + extinction(k) * sunpaths(n,k)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
      enddo

!  Enhanced-spherical, fine-layer attenuations

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_dn(n) ) then
               do j = 1, nfinedivs(n)
                  sumd = ZERO
                  do k = 1, ntraverse_fine(n,j)
                     sumd = sumd + extinction(k) * sunpaths_fine(n,k,j)
                  enddo
                  if (sumd .lt. cutoff ) Attenuations_fine(n,j) = exp( - sumd )
                  Solutions_fine(n,j) = exactscat_dn(n) * Attenuations_fine(n,j)
               enddo
            endif
         enddo
      endif

!  Regular-PS (Average secant formulation)
!    Special treatment for the horizonal case --> Factor2 = 0, lostrans = 0

!      if ( do_regular_ps ) then
!         factor1 = zero ; factor2 = zero ; Solutions = zero
!         if( Mu1 .eq. zero ) then
!            do n = 1, nlayers
!               Solutions(n) = exactscat_dn(n) * Attenuations(n-1)
!               if ( layermask_dn(n) ) then
!                  if ( attenuations(n-1).ne.zero ) then
!                     factor1(n) = Attenuations(n)/Attenuations(n-1)
!                     nstart = n
!                  endif
!               endif
!            enddo
!        else
!           do n = 1, nlayers
!               lostau = deltaus(n) / Mu1
!               Solutions(n) = exactscat_dn(n) * Attenuations(n-1)
!               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
!               if ( layermask_dn(n) ) then
!                  if ( attenuations(n-1).ne.zero ) then
!                     factor1(n) = Attenuations(n)/Attenuations(n-1)
!                     factor2(n) = (suntau(n) - suntau(n-1))/lostau
!                     nstart = n
!                  endif
!               endif
!            enddo
!         endif
!      endif

!  Layer integrated Solar sources
!  ==============================

!  Regular PS

      if ( do_regular_ps ) then
        factor1 = zero ; factor2 = zero ; Solutions = zero
        do n = nlayers, 1, -1 
          if ( Mu1 .gt. zero ) then
            lostau = deltaus(n) / Mu1
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
          endif
          if ( layermask_dn(n) .and. n.le.nstart  ) then
            if ( Mu1 .gt. zero ) then
               factor1(n) = Attenuations(n-1)*lostrans_dn(n) - Attenuations(n)
               factor2(n) = ((suntau(n) - suntau(n-1))/lostau) - one
               multiplier = factor1(n) / factor2(n)
               sources_dn(n) = exactscat_dn(n) * multiplier
            else if ( Mu1 .eq. zero ) then
               factor1(n) = Attenuations(n)
               sources_dn(n) = exactscat_dn(n) * factor1(n)
            endif
!            write(*,*)'FD1',n,exactscat_dn(n),multiplier
          endif
        enddo
      endif

!      if ( do_regular_ps ) then
!         do n = nlayers, 1, -1
!            if ( n.le.nstart ) then
!               multiplier = ( lostrans_dn(n) - Factor1(n) ) / (factor2(n) - one!)
!            else
!               multiplier = zero
!            endif  
!            sources_dn(n) = solutions(n) * multiplier
!            write(*,*)'FD1',n,solutions(n),multiplier
!         enddo
!      endif

!  Enhanced PS: special case (nadir viewing)

      if ( do_enhanced_ps .and. doNadir ) then
         do n = nlayers, 1, -1
            kn = extinction(n)
            lostrans_dn(n)  = exp ( - deltaus(n))
            rdiff = radii(n-1) - radii(n) ; if ( n.eq.NCrit) rdiff = radii(n-1) - RadCrit
            trand = one ; if ( n.eq.NCrit) trand = exp ( -kn * (RadCrit -radii(n) ) )
            if ( n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  xjkn = ( rdiff - xfine(n,j) ) * kn
                  func = solutions_fine(n,j) * exp ( - xjkn )
                  sum = sum + func * wfine(n,j)
               enddo
               sources_dn(n) = sum * kn
!  @@ Robfix, add following line
               if ( n.eq.NCrit ) sources_dn(n) = sources_dn(n) * trand
!  @@ End Robfix add line
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
            if ( n.le.nstart ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  tran = exp ( - ke * ( cotfine(n,j) - cot_c ) )   !  Down
!                 tran = exp ( - ke * ( cot_2 - cotfine(n,j) ) )   !  Up
                  func = solutions_fine(n,j) * csqfine(n,j) * tran
                  sum  = sum + func * wfine(n,j)
               enddo
               sources_dn(n) = sum * ke * trand
            endif        
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
         INTENSITY_DN(UTA) = FLUX * CUMSOURCE_DN(NC)
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  Finish

      return
end subroutine SS_Integral_1_DN

!

subroutine SS_Integral_1_UPDN   &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels,           & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_deltam_scaling,                         & ! Inputs (Flags)
     do_regular_ps, do_enhanced_ps, doNadir,                                & ! Inputs (Flags)
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,        & ! Inputs (control output)
     reflec, extinction, deltaus, omega, truncfac, phasmoms, flux,          & ! Inputs (Optical)
     Mu0, Mu1, LegPoly_up, LegPoly_dn, NCrit, RadCrit, CotCrit,             & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                   & ! Inputs (Geometry)
     sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up,        & ! Inputs (Geometry)
     sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,        & ! Inputs (Geometry)
     intensity_up, intensity_db, cumsource_up, intensity_dn, cumsource_dn )   ! Outputs

!  Stand-alone routine for Upwelling and Downwelling Solar-beam Single-scatter (SS)
!    computation of radiance. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012

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

      LOGICAL, Intent(in) ::  DO_UPWELLING
      LOGICAL, Intent(in) ::  DO_DNWELLING
      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Numbers

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
      REAL(fpk), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )

!  Solar Flux and Surface reflectivity (Could be the albedo)

      REAL(fpk), Intent(in) :: REFLEC, FLUX

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu0 = cos(theta_boa), required for the surface term (both regular and enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      integer  , Intent(in)  :: NCrit
      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers), radii(0:maxlayers)
      real(fpk), Intent(in)  :: Mu0, Mu1, RadCrit, CotCrit

!  solar paths 

      integer  , Intent(in)  :: ntraverse_up  (0:maxlayers)
      real(fpk), Intent(in)  :: sunpaths_up   (0:maxlayers,maxlayers)
      integer  , Intent(in)  :: ntraverse_fine_up(maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: sunpaths_fine_up (maxlayers,maxlayers,maxfinelayers)

      integer  , Intent(in)  :: ntraverse_dn  (0:maxlayers)
      real(fpk), Intent(in)  :: sunpaths_dn   (0:maxlayers,maxlayers)
      integer  , Intent(in)  :: ntraverse_fine_dn(maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: sunpaths_fine_dn (maxlayers,maxlayers,maxfinelayers)

!  Legendres

      REAL(fpk), Intent(in)  :: LegPoly_up(0:maxmoments_input)
      REAL(fpk), Intent(in)  :: LegPoly_dn(0:maxmoments_input)

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels )
      real(fpk), Intent(Out)  :: intensity_db     ( max_user_levels )
      real(fpk), Intent(Out)  :: cumsource_up     ( 0:maxlayers )
      real(fpk), Intent(Out)  :: intensity_dn     ( max_user_levels )
      real(fpk), Intent(Out)  :: cumsource_dn     ( 0:maxlayers )

!  Upwelling
!  ---------

   if ( do_upwelling ) then
       call SS_Integral_1_UP &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels,               & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,                 & ! Inputs (Flags)
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,            & ! Inputs (control output)
     reflec, extinction, deltaus, omega, truncfac, phasmoms, flux,              & ! Inputs (Optical)
     Mu0, Mu1, LegPoly_up, NCrit, xfine, wfine, csqfine, cotfine, Raycon, cota, & ! Inputs (Geometry)
     sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up,            & ! Inputs (Geometry)
     intensity_up, intensity_db, cumsource_up )                                   ! Outputs
   endif

   if ( do_dnwelling ) then
       call SS_Integral_1_DN &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels,              & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,                & ! Inputs (Flags)
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,           & ! Inputs (control output)
     extinction, deltaus, omega, truncfac, phasmoms, flux,                     & ! Inputs (Optical)
     Mu1, LegPoly_dn, NCrit, RadCrit, CotCrit,                                 & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                      & ! Inputs (Geometry)
     sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,           & ! Inputs (Geometry)
     intensity_dn, cumsource_dn )                                                ! Outputs
   endif

!  Finish

   return
end subroutine SS_Integral_1_UPDN

!

subroutine DTE_Integral_1_UP &
   ( maxlayers, maxfinelayers, max_user_levels, Do_Thermset,   & ! Inputs (dimensioning)
     do_dmscaling, do_regular_ps, do_enhanced_ps, doNadir,     & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,           & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                        & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                     & ! Inputs (Optical)
     Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine, & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, cumsource_up, tcom1 )      ! Outputs

!  Stand alone routine for Upwelling Direct-thermal-emission (DTE)
!    computation of radiance. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012

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

!  Ray constant, Cotangents, Critical layer
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      integer  , Intent(in)  :: NCrit
      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_dta_up ( max_user_levels )
      real(fpk), Intent(Out)  :: intensity_dts    ( max_user_levels )
      real(fpk), Intent(Out)  :: cumsource_up      ( 0:maxlayers )

!  Thermal setup

      real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)

!  LOCAL
!  -----

!  Local solutions (enhanced_ps case)

      real(fpk)  :: Solutions_fine (maxlayers,maxfinelayers)
      real(fpk)  :: Wtrans_fine    (maxlayers,maxfinelayers)

!  Source function integration results

      real(fpk)  :: sources_up       ( maxlayers )
      real(fpk)  :: lostrans_up      ( maxlayers )

!  Help

      integer    :: n, uta, nstart, nc, nut, nut_prev, j, nj
      logical    :: layermask_up(maxlayers)
      real(fpk)  :: help, sum, tran, kn, ke, xjkn, cot_1, cot_2
      real(fpk)  :: cumsource_dste, t_mult_up(0:2), thermcoeffs(2), tms, lostau

      real(fpk), parameter  :: cutoff = 88.0d0
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      CUMSOURCE_UP = zero ; INTENSITY_dta_up = zero ; INTENSITY_dts = zero
      lostrans_up = zero  ; sources_up = zero

!  Bookkeeping

      NUT = USER_LEVELS(1) + 1
      LAYERMASK_UP = .false.
      LAYERMASK_UP(NUT:NLAYERS) = .true.

      nstart = nlayers ; if (Ncrit.ne.0) nstart = nCrit
      if (nstart.lt.nlayers) LAYERMASK_UP(nstart+1:nlayers) = .false.

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

!  Plane/Parallel Layer integrated source terms
!  ============================================

      if ( do_regular_ps ) then
        if ( doNadir ) then
          DO n = 1, nlayers
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
            if ( layermask_up(n) ) then
              t_mult_up(1:2) = tcom1(n,1:2)
              sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
              t_mult_up(0) = - sum
              sources_up(n) = t_mult_up(1)
            endif
          enddo
        else
          DO n = 1, nlayers
            lostau = deltaus(n) / Mu1
            if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
            if ( layermask_up(n) ) then
              t_mult_up(2) = tcom1(n,2)
              t_mult_up(1) = tcom1(n,1) + t_mult_up(2) * Mu1
              sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
              t_mult_up(0) = - sum
              sources_up(n) = t_mult_up(0) * lostrans_up(n)  + t_mult_up(1)
            endif
          enddo
        endif
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

      if ( do_enhanced_ps ) then
         do n = nlayers, 1, -1
            kn = extinction(n) ; nj = nfinedivs(n)
            if (  doNadir  .and. layermask_up(n) ) then
               lostau = deltaus(n)
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               do j = 1, nj
                  xjkn = xfine(n,j) * kn
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = kn * exp ( - xjkn )* wfine(n,j)
              enddo
            else if ( .not. doNadir .and. layermask_up(n) ) then
               cot_2 = cota(n-1) ; cot_1 = cota(n)
               ke = raycon * kn  ; lostau = ke * ( cot_2 - cot_1 )
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               do j = 1, nj
                  xjkn = xfine(n,j) * kn
                  tran = exp ( - ke * ( cot_2 - cotfine(n,j) ) )
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = ke * tran * csqfine(n,j) * wfine(n,j)
               enddo
            endif
            sources_up(n) = dot_product(solutions_fine(n,1:nj),wtrans_fine(n,1:nj))
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
        ENDDO
         INTENSITY_DTA_UP(UTA) = CUMSOURCE_UP(NC)
         INTENSITY_DTS(UTA)    = CUMSOURCE_DSTE
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Finish

      return
end subroutine DTE_Integral_1_UP

!

subroutine DTE_Integral_1_DN &
   ( maxlayers, maxfinelayers, max_user_levels, do_Thermset,      & ! Inputs (dimensioning)
     do_dmscaling, do_regular_ps, do_enhanced_ps, doNadir,        & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,              & ! Inputs (control output)
     BB_input, extinction, deltaus, omega, truncfac,              & ! Inputs (Optical)
     Mu1, NCrit, RadCrit, CotCrit, Raycon, radii,                 & ! Inputs (Geometry)
     cota, xfine, wfine, csqfine, cotfine,                        & ! Inputs (Geometry)
     intensity_dta_dn, cumsource_dn, tcom1 )                       ! Outputs

!  Stand alone routine for Downwelling Direct-thermal-emission (DTE)
!    computation of radiance. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012

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

!  Ray constant, Cotangents, Critical layer
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      integer  , Intent(in)  :: NCrit
      real(fpk), Intent(in)  :: RadCrit, CotCrit
      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers), radii(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_dta_dn     ( max_user_levels )
      real(fpk), Intent(Out)  :: cumsource_dn     ( 0:maxlayers )

!  Thermal setup

      real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)

!  LOCAL
!  -----

!  Local solutions (enhanced_ps case)

      real(fpk)  :: Solutions_fine (maxlayers,maxfinelayers)
      real(fpk)  :: Wtrans_fine    (maxlayers,maxfinelayers)

!  Source function integration results

      real(fpk)  :: sources_dn       ( maxlayers )
      real(fpk)  :: lostrans_dn      ( maxlayers )

!  Help

      integer    :: n, uta, nstart, nc, nut, nut_prev, j, nj
      logical    :: layermask_dn(maxlayers)
      real(fpk)  :: help, sum, tran, kn, ke, xjkn, trand, lostau
      real(fpk)  :: cot_1, cot_2, rdiff, cot_c
      real(fpk)  :: t_mult_dn(0:2), thermcoeffs(2), tms

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      CUMSOURCE_DN = zero ; INTENSITY_DTA_DN = zero
      lostrans_dn  = zero ; sources_dn   = zero

!  Bookkeeping

      NUT = USER_LEVELS(N_USER_LEVELS) + 1
      IF ( NUT > NLAYERS ) NUT = NLAYERS
      LAYERMASK_DN = .false.
      LAYERMASK_DN(1:NUT) = .true.

      nstart = nlayers ; if (Ncrit.ne.0) nstart = nCrit
      if (nstart.lt.nlayers) LAYERMASK_DN(nstart+1:nlayers) = .false.

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

!  Plane/Parallel Layer integrated source terms
!  ============================================

      if ( do_regular_ps ) then
        if ( doNadir ) then
          DO n = 1, nlayers
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
            if ( layermask_dn(n) ) then
              t_mult_dn(1:2) = tcom1(n,1:2)
              t_mult_dn(0)   = - t_mult_dn(1)
              sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
              sources_dn(n)  = sum
            endif
          enddo
        else
          DO n = 1, nlayers
            lostau = deltaus(n) / Mu1
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
            if ( layermask_dn(n) ) then
              t_mult_dn(2)   = tcom1(n,2)
              t_mult_dn(1)   = tcom1(n,1) - t_mult_dn(2) * Mu1
              t_mult_dn(0)   = - t_mult_dn(1)
              sources_dn(n)  = t_mult_dn(0) * lostrans_dn(n)
              sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
              sources_dn(n)  = sources_dn(n) + sum
            endif
          enddo
        endif
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

      if ( do_enhanced_ps ) then
         do n = nlayers, 1, -1
            kn = extinction(n) ; nj = nfinedivs(n)
            if (  doNadir .and. layermask_dn(n) ) then
               rdiff = radii(n-1) - radii(n) ; if ( n.eq.NCrit) rdiff = radii(n-1) - RadCrit
               trand = one ; if ( n.eq.NCrit) trand = exp ( -kn * (RadCrit -radii(n) ) )
               lostau = deltaus(n)
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               do j = 1, nj
                  xjkn = ( rdiff - xfine(n,j) ) * kn
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = kn * exp ( - xjkn )* wfine(n,j)
               enddo
            else if ( .not. doNadir .and. layermask_dn(n)  ) then
               cot_2 = cota(n-1) ; cot_1 = cota(n)
               cot_c = cot_1     ; if ( n.eq.NCrit ) cot_c = CotCrit
               ke = raycon * kn  ; lostau = ke * ( cot_2 - cot_1 )
               trand = one  ; if ( n.eq.NCrit ) trand = exp ( - ke * ( CotCrit - cot_1 ) )
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               do j = 1, nj
                  xjkn = xfine(n,j) * kn
                  tran = exp ( - ke * ( cotfine(n,j) - cot_1 ) )   !  Down
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = ke * tran * csqfine(n,j) * wfine(n,j)
               enddo
            endif
            sources_dn(n) = dot_product(solutions_fine(n,1:nj),wtrans_fine(n,1:nj))
            if ( n.eq.NCrit ) sources_dn(n) = sources_dn(n) * trand         !@@ Robfix
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
         INTENSITY_DTA_DN(UTA) = CUMSOURCE_DN(NC)
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  Finish

      return
end subroutine DTE_Integral_1_DN



subroutine DTE_Integral_1_UPDN   &
   ( maxlayers, maxfinelayers, max_user_levels, do_Thermset,   & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_deltam_scaling,            & ! Inputs (Flags)
     do_regular_ps, do_enhanced_ps, doNadir,                   & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,           & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                        & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                     & ! Inputs (Optical)
     Mu1, NCrit, RadCrit, CotCrit, Raycon, radii,              & ! Inputs (Geometry)
     cota, xfine, wfine, csqfine, cotfine,                     & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, cumsource_up,            & ! Outputs
     intensity_dta_dn, cumsource_dn, tcom1 )                     ! Outputs

!  Stand alone routine for Upwelling and Downwelling Direct-thermal-emission (DTE)
!    computation of radiance. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012

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

!  Ray constant, Cotangents, Critical layer
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      integer  , Intent(in)  :: NCrit
      real(fpk), Intent(in)  :: RadCrit, CotCrit
      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers), radii(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_dta_up     ( max_user_levels )
      real(fpk), Intent(Out)  :: intensity_dts        ( max_user_levels )
      real(fpk), Intent(Out)  :: cumsource_up         ( 0:maxlayers )
      real(fpk), Intent(Out)  :: intensity_dta_dn     ( max_user_levels )
      real(fpk), Intent(Out)  :: cumsource_dn         ( 0:maxlayers )

!  Thermal setup

      real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)

!  Upwelling
!  ---------

   if ( do_upwelling ) then
      call DTE_Integral_1_UP &
   ( maxlayers, maxfinelayers, max_user_levels, Do_Thermset,   & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,& ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,           & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                        & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                     & ! Inputs (Optical)
     Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine, & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, cumsource_up, tcom1 )      ! Outputs
   endif

   if ( do_dnwelling ) then
       call DTE_Integral_1_DN &
   ( maxlayers, maxfinelayers, max_user_levels, do_Thermset,      & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,   & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,              & ! Inputs (control output)
     BB_input, extinction, deltaus, omega, truncfac,              & ! Inputs (Optical)
     Mu1, NCrit, RadCrit, CotCrit, Raycon, radii,                 & ! Inputs (Geometry)
     cota, xfine, wfine, csqfine, cotfine,                        & ! Inputs (Geometry)
     intensity_dta_dn, cumsource_dn, tcom1 )                       ! Outputs
   endif

!  Finish

   return
end subroutine DTE_Integral_1_UPDN

!  End module

end module FO_Scalar_RTCalcs_m

