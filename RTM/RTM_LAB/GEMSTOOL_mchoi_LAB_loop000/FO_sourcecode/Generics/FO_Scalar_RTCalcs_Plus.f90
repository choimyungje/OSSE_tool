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


module FO_Scalar_RTCalcs_Plus_m

!  For a given wavelength, this routine will calculate Corrected upwelling and downwelling
!  First Order Intensities, and any number of LCLPLS Jacobians (column/profile/surface)

!     (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.
!     (2) For the Atmospheric and Surface Direct Thermal Emission (DTE) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel or incoming solar pseudo-spherical)

!  This is Versions 1-2, without Partials. Code is stand alone with no dependencies.
!    Version 1, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1, 02 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2, 01 June     2012, R. Spurr, RT Solutions Inc.

!  Option (1). For Solar sources, the subroutines are
!       SS_Integral_LPLCLS_Plus_1_UP   (Upwelling only)
!       SS_Integral_LPLCLS_Plus_1_DN   (Downwelling only)
!       SS_Integral_LPLCLS_Plus_1_UPDN (Upwelling and Downwelling)

!  Option (2). For Thermal Emission sources, the subroutines are
!       DTE_Integral_LPLCLS_Plus_1_UP   (Upwelling only)
!       DTE_Integral_LPLCLS_Plus_1_DN   (Downwelling only)
!       DTE_Integral_LPLCLS_Plus_1_UPDN (Upwelling and Downwelling)

!  All subroutines public

public

contains


subroutine SS_Integral_LPLCLS_Plus_1_UP &
   ( maxlayers, maxfinelayers, maxmoments_input,                                & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs,                             & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,                 & ! Inputs (Flags - General  )
     do_columnwfs, do_profilewfs, do_surfacewfs,                                & ! Inputs (Flags - Jacobian )
     nlayers, nfinedivs, NCrit, nmoments_input, n_user_levels, user_levels,     & ! Inputs (Layers/Levels control)
     n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,               & ! Inputs (Jacobian control)
     Mu0, Mu1, LegPoly_up, Raycon, cota, sunpaths, ntraverse,                   & ! Inputs (Whole layer Geometry)
     xfine, wfine, csqfine, cotfine, sunpaths_fine, ntraverse_fine,             & ! Inputs (Fine  layer Geometry)
     Flux, reflec, extinction, deltaus, omega, truncfac, phasmoms,              & ! Inputs (Optical - Regular)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,       & ! Inputs (Optical - Linearized)
     intensity_up, intensity_db_up, LC_Jacobians_up, LP_Jacobians_up,           & ! Output
     LC_Jacobians_db_up, LP_Jacobians_db_up, LS_Jacobians_db_up )                 ! Output

!  Stand alone routine for SS field (Radiances and Jacobians) with Solar sources alone
!    No partials

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Lines 1-2. Dimensions

      INTEGER, Intent(in) :: maxlayers
      INTEGER, Intent(in) :: maxfinelayers
      INTEGER, Intent(in) :: maxmoments_input
      INTEGER, Intent(in) :: max_user_levels
      INTEGER, Intent(in) :: max_atmoswfs
      INTEGER, Intent(in) :: max_surfacewfs

!  Line 3. General flags

      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Line 4. Jacobian Flags

      LOGICAL, Intent(in) ::  do_surfacewfs
      LOGICAL, Intent(in) ::  do_columnwfs
      LOGICAL, Intent(in) ::  do_profilewfs

!  Line 5. Layer and Level Control Numbers, Number of Moments

      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS), NCRIT
      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )
      INTEGER, Intent(in) ::  NMOMENTS_INPUT

!  Line 6. Jacobian control

      INTEGER, Intent(in) ::  n_columnwfs
      INTEGER, Intent(in) ::  n_surfacewfs
      LOGICAL, Intent(in) ::  Lvaryflags(maxlayers)
      INTEGER, Intent(in) ::  Lvarynums (maxlayers)
      LOGICAL, Intent(in) ::  Lvarymoms (maxlayers,max_atmoswfs)

!  Line 7. Miscellaneous Geometrical inputs
!       Ray constant, Cotangents
!       Mu0 = cos(theta_boa), required for the surface term (both regular and enhanced)
!       Mu1 = cos(alpha_boa), required for the Regular PS only
!       solar paths, Legendre Polynomials

      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1, Mu0
      integer  , Intent(in)  :: ntraverse  (0:maxlayers)
      real(fpk), Intent(in)  :: sunpaths   (0:maxlayers,maxlayers)
      REAL(fpk), Intent(in)  :: LegPoly_up(0:maxmoments_input)

!  Line 8. Fine-layering for Enhanced PS Calculation

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)
      integer  , Intent(in)  :: ntraverse_fine(maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)

!  Line 9. Optical inputs for Atmosphere, Surface reflectivity

      REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
      REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
      REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
      REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
      REAL(fpk), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )
      REAL(fpk), Intent(in) :: REFLEC, FLUX

!  Line 10. Linearized optical inputs

      REAL(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT, max_atmoswfs )
      REAL(fpk), Intent(in) :: LS_REFLEC     ( max_surfacewfs)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels )
      real(fpk), Intent(Out)  :: intensity_db_up  ( max_user_levels )
      real(fpk), Intent(Out)  :: LC_Jacobians_up  ( max_user_levels, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_up  ( max_user_levels, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LC_Jacobians_db_up  ( max_user_levels, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_db_up  ( max_user_levels, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LS_Jacobians_db_up  ( max_user_levels, max_surfacewfs )

!  LOCAL
!  -----

!  Attenuations

      real(fpk)  :: attenuations     (0:maxlayers)
      real(fpk)  :: L_attenuations   (0:maxlayers,0:maxlayers,max_atmoswfs)

!  Solutions for Enhaced-PS

      real(fpk)  :: Solutions_fine   (maxlayers,maxfinelayers)
      real(fpk)  :: L_solutions_fine (maxlayers,maxfinelayers,0:maxlayers,max_atmoswfs)

!  Scattering

      real(fpk)  :: tms            (maxlayers)
      real(fpk)  :: exactscat_up   (maxlayers)
      real(fpk)  :: L_tms          (maxlayers,max_atmoswfs)
      real(fpk)  :: L_exactscat_up (maxlayers,max_atmoswfs)

!  Source function integration results

      real(fpk) :: sources_up         ( maxlayers )
      real(fpk) :: L_sources_up       ( maxlayers,0:maxlayers,max_atmoswfs )

      real(fpk) :: lostrans_up        ( maxlayers )
      real(fpk) :: L_lostrans_up      ( maxlayers,max_atmoswfs )

      real(fpk) :: cumsource_db      ( 0:maxlayers )
      real(fpk) :: cumsource_up      ( 0:maxlayers )
      real(fpk) :: L_cumsource       ( max_atmoswfs )
      real(fpk) :: LS_cumsource      ( max_surfacewfs )

!  Help

      integer    :: n, k, kc, j, q, L, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
      logical    :: do_atmoswfs, layermask_up(maxlayers), Qvary(maxlayers)

      real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers), func(maxfinelayers)

      real(fpk)  :: cons, help, sum, tran_1, kn, ke, factor1, factor2, m4, term1
      real(fpk)  :: L_help, L_sum, L_tran, L_func, L_factor1, L_factor2
      real(fpk)  :: cot_1, cot_2, multiplier, suntau(0:maxlayers), lostau
      real(fpk)  :: L_multiplier, L_suntau(0:maxlayers,0:maxlayers,max_atmoswfs), L_lostau(max_atmoswfs)
      real(fpk)  :: attenuations_fine, L_attenuations_fine, sumd, L_sumd, L_exactscat

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      INTENSITY_UP = zero
      LP_JACOBIANS_UP = zero
      LC_JACOBIANS_UP = zero

      INTENSITY_DB_UP = zero
      LP_JACOBIANS_DB_UP = zero
      LC_JACOBIANS_DB_UP = zero
      LS_JACOBIANS_DB_UP = zero

!  Zero local sources

      lostrans_up   = zero  ; sources_up   = zero ; cumsource_up = zero
      L_lostrans_up = zero  ; L_sources_up = zero

!  Bookkeeping

      NUT = USER_LEVELS(1) + 1
      LAYERMASK_UP = .false.
      LAYERMASK_UP(NUT:NLAYERS) = .true.

      do_atmoswfs = do_profilewfs .or. do_columnwfs
      Qvary = .false. ; QNums = 0
      if ( do_profilewfs ) then
         Qvary(1:nlayers) = Lvaryflags(1:nlayers)
         QNums(1:nlayers) = Lvarynums (1:nlayers)
      else if ( do_columnwfs ) then
         Qvary(1:nlayers) = .true.
         QNums(1:nlayers) =  n_columnwfs
      endif
      kc = 0

!  TMS factors and linearizations

      if ( do_deltam_scaling ) then
         do n = 1, nlayers
            help = one - truncfac(n) * omega(n)
            tms(n) = omega(n) / help
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
                  L_tms(n,q) = ( L_omega(n,q) - tms(n)*L_help ) / help
               enddo
            endif
         enddo
      else
         do n = 1, nlayers
            tms(n) = omega(n)
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_tms(n,q) = L_omega(n,q)
               enddo
            endif
         enddo
      endif
 
!  Scattering functions and Linearization

      do n = 1, nlayers
         if ( layermask_up(n) ) then
            sum = zero
            do L = 0, nmoments_input
               sum = sum + LegPoly_Up(L) * phasmoms(n,L)
            enddo
            exactscat_up(n) = sum * tms(n)
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) then
                     L_sum = zero
                     do L = 0, nmoments_input
                        L_sum = L_sum + LegPoly_Up(L) * L_phasmoms(n,L,q)
                     enddo
                     L_exactscat_up(n,q) = L_sum * tms(n) + sum * L_tms(n,q)
                  else
                     L_exactscat_up(n,q) = sum * L_tms(n,q)
                  endif
               enddo
            endif               
         endif
      enddo

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations = ZERO   ; L_Attenuations = ZERO
      Suntau       = ZERO   ; L_suntau       = ZERO
      Solutions_fine = zero ; L_Solutions_fine = zero

      nstart = nlayers ; if (Ncrit.ne.0) nstart = nCrit

!  Attenuations to End points (including TOA). Both PS representations
!  ===================================================================

      do n = 0, nlayers
         sumd = ZERO
         do k = 1, ntraverse(n)
            sumd = sumd + extinction(k) * sunpaths(n,k)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if ( do_profilewfs ) then
            do k = 1, nlayers
               if ( Qvary(k) .and. k.le.ntraverse(n) ) then
                  do q = 1, Qnums(k)
                     L_suntau(n,k,q) = L_extinction(k,q) * sunpaths(n,k)
                     L_Attenuations(n,k,q) = - Attenuations(n) * L_suntau(n,k,q)
                  enddo
               endif
            enddo
         else if ( do_columnwfs ) then
            do q = 1, n_columnwfs
               L_sumd = ZERO
               do k = 1, ntraverse(n)
                  L_sumd = L_sumd + L_extinction(k,q) * sunpaths(n,k)
               enddo
               L_suntau(n,kc,q) = L_sumd
               L_Attenuations(n,kc,q) = - Attenuations(n) * L_sumd
            enddo
         endif
      enddo

!  Enhanced-spherical, fine-layer attenuations
!  ===========================================

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_up(n) ) then
               do j = 1, nfinedivs(n)
                  sumd = ZERO
                  do k = 1, ntraverse_fine(n,j)
                     sumd = sumd + extinction(k) * sunpaths_fine(n,k,j)
                  enddo
                  if (sumd .lt. cutoff ) Attenuations_fine = exp( - sumd )
                  Solutions_fine(n,j) = exactscat_up(n) * Attenuations_fine
                  if ( do_profilewfs ) then
                     do k = 1, nlayers
                        if ( Qvary(k) .and. k.le.ntraverse_fine(n,j) ) then
                           do q = 1, Qnums(k)
                               L_exactscat = zero ; if ( k.eq.n) L_exactscat = L_exactscat_up(n,q)
                               L_sumd = L_extinction(k,q) * sunpaths_fine(n,k,j)
                               L_Attenuations_fine = - Attenuations_fine * L_sumd 
                               L_Solutions_fine(n,j,k,q) = L_exactscat       *   Attenuations_fine   + &
                                                             exactscat_up(n) * L_Attenuations_fine
                           enddo
                        endif
                     enddo
                  else if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_sumd = ZERO
                        do k = 1, ntraverse_fine(n,j)
                           L_sumd = L_sumd + L_extinction(k,q) * sunpaths_fine(n,k,j)
                        enddo
                        L_Attenuations_fine = - Attenuations_fine * L_sumd 
                        L_Solutions_fine(n,j,kc,q) = L_exactscat_up(n,q) *   Attenuations_fine   + &
                                                       exactscat_up(n)   * L_Attenuations_fine
!            if(n.eq.Ncrit.and.j.eq.2)write(*,*)n,q,j,L_Solutions_fine(n,j,kc,q), Solutions_fine(n,j)
                     enddo
                  endif
               enddo
            endif
         enddo
      endif

!  Regular-PS (Average secant formulation) Layer integrated Solar sources
!  ======================================================================

      if ( do_regular_ps ) then
         do n = nlayers, 1, -1

!  LOS transmittance (not for the Horizontal case)

            if ( Mu1 .gt. zero ) then
               lostau = deltaus(n)  / Mu1
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( do_atmoswfs ) then
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_lostau(q)        = L_deltaus(n,q) / Mu1
                        L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                     enddo
                  endif
               endif
            endif

!  Sources, general case

            if ( layermask_up(n) .and. n.le.nstart  ) then
              if ( Mu1 .gt. zero ) then
                factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                factor2 = one + (suntau(n) - suntau(n-1))/lostau
                multiplier = factor1 / factor2
                sources_up(n) = exactscat_up(n) * multiplier
                if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           L_factor1 = L_Attenuations(n-1,k,q) - L_Attenuations(n,k,q)*lostrans_up(n)
                           L_factor2 = (L_suntau(n,k,q) - L_suntau(n-1,k,q))/lostau
                           if ( k.eq.n) then
                              L_factor1 = L_factor1 - Attenuations(n)*L_lostrans_up(n,q)
                              L_factor2 = L_factor2 - (factor2 - one)*L_lostau(q)/lostau 
                              L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                              L_sources_up(n,k,q) = L_exactscat_up(n,q) *   multiplier + &
                                                      exactscat_up(n)   * L_multiplier
                           else
                              L_multiplier = ( L_factor1 - multiplier * L_factor2 ) / factor2
                              L_sources_up(n,k,q) = exactscat_up(n)   * L_multiplier
                           endif
                        enddo
                     endif
                  enddo
                else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_factor1 = L_Attenuations(n-1,kc,q) - L_Attenuations(n,kc,q)*lostrans_up(n)
                     L_factor2 = (L_suntau(n,kc,q) - L_suntau(n-1,kc,q))/lostau
                     L_factor1 = L_factor1 - Attenuations(n)*L_lostrans_up(n,q)
                     L_factor2 = L_factor2 - (factor2 - one)*L_lostau(q)/lostau 
                     L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                     L_sources_up(n,kc,q) = L_exactscat_up(n,q) *   multiplier + &
                                              exactscat_up(n)   * L_multiplier
                  enddo
                endif
              endif
            endif

!  Sources, special case (horizonal view)

            if ( layermask_up(n) .and. n.le.nstart  ) then
              if ( Mu1 .eq. zero ) then
                factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                sources_up(n) = exactscat_up(n) * factor1
                if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           L_factor1 = L_Attenuations(n-1,k,q)
                           if ( k.eq.n) then
                              L_sources_up(n,k,q) = L_exactscat_up(n,q) *   factor1  + &
                                                      exactscat_up(n)   * L_factor1
                           else
                              L_sources_up(n,k,q) = exactscat_up(n)   * L_factor1
                           endif
                        enddo
                     endif
                  enddo
                else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_factor1 = L_Attenuations(n-1,kc,q)
                     L_sources_up(n,kc,q) = L_exactscat_up(n,q) *   factor1 + &
                                              exactscat_up(n)   * L_factor1
                  enddo
                endif
              endif
            endif

!  End layers and regular-PS formulation

         enddo
      endif

!  Enhanced PS: special case (nadir viewing). Layer integrated Solar sources
!  =========================================================================

      if ( do_enhanced_ps .and. doNadir ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            kn     = extinction(n)
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
            if ( do_atmoswfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_deltaus(n,q)
                     L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
                  enddo
               endif
            endif

!  Sources

            if ( layermask_up(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  argum(j) = xfine(n,j)
                  tran(j)  = exp ( - argum(j) * kn )
                  func(j)  = solutions_fine(n,j) * tran(j)
                  sum = sum + func(j) * wfine(n,j)
               enddo
               sources_up(n) = sum * kn
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = L_solutions_fine(n,j,k,q) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_sources_up(n,k,q)  = L_sum * kn + L_extinction(N,q) * sum
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_func = L_solutions_fine(n,j,k,q)  * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_sources_up(n,k,q)  = L_sum * kn
                           endif
                        enddo
                     endif
                  enddo
               else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = L_solutions_fine(n,j,kc,q) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j)
                     enddo
                     L_sources_up(n,kc,q)  = L_sum * kn + L_extinction(N,q) * sum
                  enddo
               endif
            endif

!  End layer loop and Nadir Enhanced PS case

         enddo
      endif

!  Enhanced PS: General case. Layer integrated Solar sources
!  =========================================================

      if ( do_enhanced_ps .and. .not. doNadir ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            cot_2 = cota(n-1) ; cot_1 = cota(n)
            kn = extinction(n) ;  ke = raycon * kn ; cons = raycon * ( cot_2 - cot_1 )
            tran_1 = kn * cons
            if ( tran_1 .lt. cutoff ) lostrans_up(n) = exp ( - tran_1 )
            if ( do_atmoswfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_extinction(n,q) * cons
                     L_lostrans_up(n,q) = - L_lostau(q) * lostrans_up(n)
!                     if(n.eq.NCrit)write(*,*)n,q,L_lostrans_up(n,q),lostrans_up(n)
                  enddo
               endif
            endif

!  Sources

            if ( layermask_up(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  argum(j) = Raycon * ( cot_2 - cotfine(n,j) )
                  tran(j)  = exp ( - kn * argum(j) )
                  func(j)  = solutions_fine(n,j) * csqfine(n,j) * tran(j)
                  sum      = sum + func(j) * wfine(n,j)
               enddo
               sources_up(n) = sum * ke 
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = L_solutions_fine(n,j,k,q) * csqfine(n,j) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_sources_up(n,k,q)  = L_sum * ke + L_extinction(N,q) * Raycon * sum
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_func = L_solutions_fine(n,j,k,q) * csqfine(n,j) * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_sources_up(n,k,q)  = L_sum * ke
                           endif
                        enddo
                     endif
                  enddo
               else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = L_solutions_fine(n,j,kc,q) * csqfine(n,j) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j)
                     enddo
                     L_sources_up(n,kc,q)  = L_sum * ke + L_extinction(N,q) * Raycon * sum
                  enddo
               endif
            endif

!  End layer loop and general Enhanced PS case

         enddo
      endif

!  Source function integration
!  ===========================

!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!  Check for updating the recursion

!  INTENSITY Main loop over all output optical depths
!          Cumulative source term will be saved

      NC = 0 ; M4 = 4.0_fpk * Mu0 
      CUMSOURCE_UP(NC) = zero
      CUMSOURCE_DB(NC) = M4 * REFLEC * attenuations(nlayers)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1
      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_DB(NC) = LOSTRANS_UP(N) * CUMSOURCE_DB(NC-1)
            CUMSOURCE_UP(NC) = SOURCES_UP(N) + LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1)
         ENDDO
         INTENSITY_DB_UP(UTA) = FLUX * CUMSOURCE_DB(NC)
         INTENSITY_UP(UTA)    = FLUX * CUMSOURCE_UP(NC)
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Surface WFs

      if ( do_surfacewfs ) then
         Term1 = M4  * attenuations(nlayers)
         do q = 1, n_surfacewfs
            LS_cumsource(q) = term1 * LS_reflec(q)
         enddo
         NSTART = NLAYERS
         NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_surfacewfs
                  LS_cumsource(q) = LOSTRANS_UP(N) * LS_CUMSOURCE(Q)
               enddo
            ENDDO
            do q = 1, n_surfacewfs
               LS_JACOBIANS_DB_UP(UTA,Q) = FLUX * LS_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Profile Wfs (atmospheric term)

      if ( do_profilewfs ) then
         do k = 1, nlayers
            if ( Qvary(k) ) then
               L_CUMSOURCE = zero
               NSTART = NLAYERS
               NUT_PREV = NSTART + 1
               DO UTA = N_USER_LEVELS, 1, -1
                  NUT    = USER_LEVELS(UTA) + 1
                  DO N = NSTART, NUT, -1
                     NC = NLAYERS + 1 - N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(q) = L_SOURCES_UP(N,K,Q)         + &
                                L_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1) + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = L_SOURCES_UP(N,K,Q)         + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_JACOBIANS_UP(UTA,K,Q) = FLUX * L_CUMSOURCE(Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Profile Wfs (direct beam term)

      if ( do_profilewfs ) then
         term1 = M4 * reflec
         do k = 1, nlayers
            if ( Qvary(k) ) then
               do q = 1, Qnums(k)
                  L_CUMSOURCE(q) = Term1 * L_attenuations(nlayers,k,q)
               enddo
               NSTART = NLAYERS
               NUT_PREV = NSTART + 1
               DO UTA = N_USER_LEVELS, 1, -1
                  NUT    = USER_LEVELS(UTA) + 1
                  DO N = NSTART, NUT, -1
                     NC = NLAYERS + 1 - N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(q) =  L_LOSTRANS_UP(N,Q) * CUMSOURCE_DB(NC-1) + &
                                               LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LOSTRANS_UP(N) * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_JACOBIANS_DB_UP(UTA,K,Q) = FLUX * L_CUMSOURCE(Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Column Wfs (Atmospheric term)

      if ( do_columnwfs ) then
         L_CUMSOURCE = zero
         NSTART = NLAYERS
         NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(q) = L_SOURCES_UP(N,KC,Q)         + &
                                L_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1) + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_JACOBIANS_UP(UTA,Q) = FLUX * L_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Column Wfs (Surface term)

      if ( do_columnwfs ) then
         term1 = M4 * reflec
         do q = 1, n_columnwfs
            L_CUMSOURCE(q) = Term1 * L_attenuations(nlayers,kc,q)
         enddo
         NSTART = NLAYERS
         NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(q) =  L_LOSTRANS_UP(N,Q) * CUMSOURCE_DB(NC-1) + &
                                      LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_JACOBIANS_DB_UP(UTA,Q) = FLUX * L_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Finish

      return
end subroutine SS_Integral_LPLCLS_Plus_1_UP

!

subroutine SS_Integral_LPLCLS_Plus_1_DN &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels, max_atmoswfs,               & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir, do_columnwfs, do_profilewfs,  & ! Inputs (Flags - Jacobian )
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,          & ! Inputs (Layers/Levels control)
     n_columnwfs, Lvaryflags, Lvarynums, Lvarymoms, NCrit, RadCrit, CotCrit,  & ! Inputs (Jacobian control)
     Mu1, LegPoly_dn, Raycon, radii, cota, sunpaths, ntraverse,               & ! Inputs (Whole layer Geometry)
     xfine, wfine, csqfine, cotfine, sunpaths_fine, ntraverse_fine,           & ! Inputs (Fine  layer Geometry)
     Flux, extinction, deltaus, omega, truncfac, phasmoms,                    & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,                & ! Inputs (Optical - Linearized)
     intensity_dn, LC_Jacobians_dn, LP_Jacobians_dn )                           ! Output

!  Stand alone routine for SS field (Radiances and Jacobians) with Solar sources alone
!    No partials

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Line 1. Dimensions

      INTEGER, Intent(in) :: maxlayers
      INTEGER, Intent(in) :: maxfinelayers
      INTEGER, Intent(in) :: maxmoments_input
      INTEGER, Intent(in) :: max_user_levels
      INTEGER, Intent(in) :: max_atmoswfs

!  Line 2. General and Jacobianflags

      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DONADIR
      LOGICAL, Intent(in) ::  do_columnwfs
      LOGICAL, Intent(in) ::  do_profilewfs

!  Line 3. Layer and Level Control Numbers, Number of Moments

      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )
      INTEGER, Intent(in) ::  NMOMENTS_INPUT

!  Line 4. Jacobian control

      INTEGER, Intent(in) ::  n_columnwfs
      LOGICAL, Intent(in) ::  Lvaryflags(maxlayers)
      INTEGER, Intent(in) ::  Lvarynums (maxlayers)
      LOGICAL, Intent(in) ::  Lvarymoms (maxlayers,max_atmoswfs)

!  Line 4. Crit layer control

      INTEGER  , Intent(in) :: NCRIT
      real(fpk), Intent(in) :: Radcrit, CotCrit


!  Line 5. Miscellaneous Geometrical inputs
!       Ray constant, Cotangents, radii
!       Mu1 = cos(alpha_boa), required for the Regular PS only
!       solar paths, Legendre Polynomials

      real(fpk), Intent(in)  :: Raycon, radii(0:maxlayers), cota(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1
      integer  , Intent(in)  :: ntraverse  (0:maxlayers)
      real(fpk), Intent(in)  :: sunpaths   (0:maxlayers,maxlayers)
      REAL(fpk), Intent(in)  :: LegPoly_dn(0:maxmoments_input)

!  Line 6. Fine-layering for Enhanced PS Calculation

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)
      integer  , Intent(in)  :: ntraverse_fine(maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)

!  Line 7. Optical inputs for Atmosphere

      REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
      REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
      REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
      REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
      REAL(fpk), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )
      REAL(fpk), Intent(in) :: FLUX

!  Line 8. Linearized optical inputs

      REAL(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT, max_atmoswfs )

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_dn     ( max_user_levels )
      real(fpk), Intent(Out)  :: LC_Jacobians_dn  ( max_user_levels, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dn  ( max_user_levels, maxlayers, max_atmoswfs )

!  LOCAL
!  -----

!  Attenuations

      real(fpk)  :: attenuations     (0:maxlayers)
      real(fpk)  :: L_attenuations   (0:maxlayers,0:maxlayers,max_atmoswfs)

!  Solutions for Enhaced-PS

      real(fpk)  :: Solutions_fine   (maxlayers,maxfinelayers)
      real(fpk)  :: L_solutions_fine (maxlayers,maxfinelayers,0:maxlayers,max_atmoswfs)

!  Scattering

      real(fpk)  :: tms            (maxlayers)
      real(fpk)  :: exactscat_dn   (maxlayers)
      real(fpk)  :: L_tms          (maxlayers,max_atmoswfs)
      real(fpk)  :: L_exactscat_dn (maxlayers,max_atmoswfs)

!  Source function integration results

      real(fpk) :: sources_dn         ( maxlayers )
      real(fpk) :: L_sources_dn       ( maxlayers,0:maxlayers,max_atmoswfs )

      real(fpk) :: lostrans_dn        ( maxlayers )
      real(fpk) :: L_lostrans_dn      ( maxlayers,max_atmoswfs )

      real(fpk) :: cumsource_dn      ( 0:maxlayers )
      real(fpk) :: L_cumsource       ( max_atmoswfs )

!  Help

      integer    :: n, k, kc, j, q, L, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
      logical    :: do_atmoswfs, layermask_dn(maxlayers), Qvary(maxlayers)

      real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers), func(maxfinelayers)

      real(fpk)  :: cons, help, sum, tran_1, kn, ke, factor1, factor2
      real(fpk)  :: L_help, L_sum, L_tran, L_func, L_factor1, L_factor2
      real(fpk)  :: cot_1, cot_2, multiplier, suntau(0:maxlayers), lostau, rdiff, cot_c
      real(fpk)  :: L_multiplier, L_suntau(0:maxlayers,0:maxlayers,max_atmoswfs), L_lostau(max_atmoswfs)
      real(fpk)  :: attenuations_fine, L_attenuations_fine, sumd, L_sumd, L_exactscat

!  @@ Robfix, add following line (New variables)
      real(fpk)  :: consc, trand
!  @@ End Robfix add line

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      INTENSITY_DN = zero
      LP_JACOBIANS_DN = zero
      LC_JACOBIANS_DN = zero

!  Zero local sources

      lostrans_dn   = zero  ; sources_dn   = zero ; cumsource_dn = zero
      L_lostrans_dn = zero  ; L_sources_dn = zero

!  Bookkeeping

      NUT = USER_LEVELS(N_USER_LEVELS) + 1
      IF ( NUT > NLAYERS ) NUT = NLAYERS
      LAYERMASK_DN = .false.
      LAYERMASK_DN(1:NUT) = .true.

      do_atmoswfs = do_profilewfs .or. do_columnwfs
      Qvary = .false. ; QNums = 0
      if ( do_profilewfs ) then
         Qvary(1:nlayers) = Lvaryflags(1:nlayers)
         QNums(1:nlayers) = Lvarynums (1:nlayers)
      else if ( do_columnwfs ) then
         Qvary(1:nlayers) = .true.
         QNums(1:nlayers) =  n_columnwfs
      endif
      kc = 0

!  TMS factors and linearizations

      if ( do_deltam_scaling ) then
         do n = 1, nlayers
            help = one - truncfac(n) * omega(n)
            tms(n) = omega(n) / help
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
                  L_tms(n,q) = ( L_omega(n,q) - tms(n)*L_help ) / help
               enddo
            endif
         enddo
      else
         do n = 1, nlayers
            tms(n) = omega(n)
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_tms(n,q) = L_omega(n,q)
               enddo
            endif
         enddo
      endif
 
!  Scattering functions and Linearization

      do n = 1, nlayers
         if ( layermask_dn(n) ) then
            sum = zero
            do L = 0, nmoments_input
               sum = sum + LegPoly_dn(L) * phasmoms(n,L)
            enddo
            exactscat_dn(n) = sum * tms(n)
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  if ( Lvarymoms(n,q) ) then
                     L_sum = zero
                     do L = 0, nmoments_input
                        L_sum = L_sum + LegPoly_dn(L) * L_phasmoms(n,L,q)
                     enddo
                     L_exactscat_dn(n,q) = L_sum * tms(n) + sum * L_tms(n,q)
                  else
                     L_exactscat_dn(n,q) = sum * L_tms(n,q)
                  endif
               enddo
            endif               
         endif
      enddo

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      Attenuations = ZERO   ; L_Attenuations = ZERO
      Suntau       = ZERO   ; L_suntau       = ZERO
      Solutions_fine = zero ; L_Solutions_fine = zero

      nstart = nlayers ; if (Ncrit.ne.0) nstart = nCrit

!  @@ Robfix, add following line (for initialization; do not leave unassigned)
      consc = zero ; trand = one
!  @@ End Robfix add line

!  Attenuations to End points (including TOA). Both PS representations
!  ===================================================================

      do n = 0, nlayers
         sumd = ZERO
         do k = 1, ntraverse(n)
            sumd = sumd + extinction(k) * sunpaths(n,k)
         enddo
         suntau(n) = sumd
         If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if ( do_profilewfs ) then
            do k = 1, nlayers
               if ( Qvary(k) .and. k.le.ntraverse(n) ) then
                  do q = 1, Qnums(k)
                     L_suntau(n,k,q) = L_extinction(k,q) * sunpaths(n,k)
                     L_Attenuations(n,k,q) = - Attenuations(n) * L_suntau(n,k,q)
                  enddo
               endif
            enddo
         else if ( do_columnwfs ) then
            do q = 1, n_columnwfs
               L_sumd = ZERO
               do k = 1, ntraverse(n)
                  L_sumd = L_sumd + L_extinction(k,q) * sunpaths(n,k)
               enddo
               L_suntau(n,kc,q) = L_sumd
               L_Attenuations(n,kc,q) = - Attenuations(n) * L_sumd
            enddo
         endif
      enddo

!  Enhanced-spherical, fine-layer attenuations
!  ===========================================

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_dn(n) ) then
               do j = 1, nfinedivs(n)
                  sumd = ZERO
                  do k = 1, ntraverse_fine(n,j)
                     sumd = sumd + extinction(k) * sunpaths_fine(n,k,j)
                  enddo
                  if (sumd .lt. cutoff ) Attenuations_fine = exp( - sumd )
                  Solutions_fine(n,j) = exactscat_dn(n) * Attenuations_fine
                  if ( do_profilewfs ) then
                     do k = 1, nlayers
                        if ( Qvary(k) .and. k.le.ntraverse_fine(n,j) ) then
                           do q = 1, Qnums(k)
                               L_exactscat = zero ; if ( k.eq.n) L_exactscat = L_exactscat_dn(n,q)
                               L_sumd = L_extinction(k,q) * sunpaths_fine(n,k,j)
                               L_Attenuations_fine = - Attenuations_fine * L_sumd 
                               L_Solutions_fine(n,j,k,q) = L_exactscat       *   Attenuations_fine   + &
                                                             exactscat_dn(n) * L_Attenuations_fine
                           enddo
                        endif
                     enddo
                  else if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_sumd = ZERO
                        do k = 1, ntraverse_fine(n,j)
                           L_sumd = L_sumd + L_extinction(k,q) * sunpaths_fine(n,k,j)
                        enddo
                        L_Attenuations_fine = - Attenuations_fine * L_sumd 
                        L_Solutions_fine(n,j,kc,q) = L_exactscat_dn(n,q) *   Attenuations_fine   + &
                                                       exactscat_dn(n)   * L_Attenuations_fine
                     enddo
                  endif
               enddo
            endif
         enddo
      endif

!  Regular-PS (Average secant formulation) Layer integrated Solar sources
!  ======================================================================

      if ( do_regular_ps ) then
         do n = nlayers, 1, -1

!  LOS transmittance (Not for the Horizontal View)

            if ( Mu1 .gt. zero ) then
               lostau = deltaus(n) / Mu1
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               if ( do_atmoswfs ) then
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_lostau(q)        = L_deltaus(n,q)  / Mu1
                        L_lostrans_dn(n,q) = - L_lostau(q) * lostrans_dn(n)
                     enddo
                  endif
               endif
            endif

!  Sources, general case

            if ( layermask_dn(n) .and. n.le.nstart  ) then
              if ( Mu1 .gt. zero ) then
                factor1 = Attenuations(n-1)*lostrans_dn(n) - Attenuations(n)
                factor2 = ((suntau(n) - suntau(n-1))/lostau) - one
                multiplier = factor1 / factor2
                sources_dn(n) = exactscat_dn(n) * multiplier
!                write(*,*)'Line',n,exactscat_dn(n),multiplier
                if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           L_factor1 = L_Attenuations(n-1,k,q)*lostrans_dn(n) - L_Attenuations(n,k,q)
                           L_factor2 = (L_suntau(n,k,q) - L_suntau(n-1,k,q))/lostau
                           if ( k.eq.n) then
                              L_factor1 = L_factor1 + Attenuations(n-1)*L_lostrans_dn(n,q)
                              L_factor2 = L_factor2 - (factor2 + one)*L_lostau(q)/lostau 
                              L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                              L_sources_dn(n,k,q) = L_exactscat_dn(n,q) *   multiplier + &
                                                      exactscat_dn(n)   * L_multiplier
                           else
                              L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                              L_sources_dn(n,k,q) = exactscat_dn(n)   * L_multiplier
                           endif
                        enddo
                     endif
                  enddo
                else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_factor1 = L_Attenuations(n-1,kc,q)*lostrans_dn(n) - L_Attenuations(n,kc,q)
                     L_factor2 = (L_suntau(n,kc,q) - L_suntau(n-1,kc,q))/lostau
                     L_factor1 = L_factor1 + Attenuations(n-1)*L_lostrans_dn(n,q)
                     L_factor2 = L_factor2 - (factor2 + one)*L_lostau(q)/lostau 
                     L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                     L_sources_dn(n,kc,q) = L_exactscat_dn(n,q) *   multiplier + &
                                              exactscat_dn(n)   * L_multiplier
                  enddo
                endif
              endif
            endif

!  Sources, Special case

            if ( layermask_dn(n) .and. n.le.nstart  ) then
              if ( Mu1 .eq. zero ) then
                factor1 = Attenuations(n)
                sources_dn(n) = exactscat_dn(n) * factor1
                if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           L_factor1 = L_Attenuations(n,k,q)
                           if ( k.eq.n) then
                              L_sources_dn(n,k,q) = L_exactscat_dn(n,q) *   factor1 + &
                                                      exactscat_dn(n)   * L_factor1
                           else
                              L_sources_dn(n,k,q) = exactscat_dn(n)   * L_factor1
                           endif
                        enddo
                     endif
                  enddo
                else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_factor1 = L_Attenuations(n,kc,q)
                     L_sources_dn(n,kc,q) = L_exactscat_dn(n,q) *   factor1 + &
                                              exactscat_dn(n)   * L_factor1
                  enddo
                endif
              endif
            endif

!  End layers and regular-PS formulation

         enddo
      endif

!  Enhanced PS: special case (nadir viewing). Layer integrated Solar sources
!  =========================================================================

      if ( do_enhanced_ps .and. doNadir ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            kn     = extinction(n)
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
            if ( do_atmoswfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_deltaus(n,q)
                     L_lostrans_dn(n,q) = - L_lostau(q) * lostrans_dn(n)
                  enddo
               endif
            endif

!  Sources

            if ( layermask_dn(n) .and. n.le.nstart  ) then
               rdiff = radii(n-1) - radii(n)
               if ( n.eq.NCrit) rdiff = radii(n-1) - RadCrit
               sum = zero
               do j = 1, nfinedivs(n)
                  argum(j) = rdiff - xfine(n,j)
                  tran(j)  = exp ( - argum(j) * kn )
                  func(j)  = solutions_fine(n,j) * tran(j)
                  sum = sum + func(j) * wfine(n,j)
               enddo
               sources_dn(n) = sum * kn

!  @@ Robfix, add following 4 lines
               if ( n.eq.NCrit ) then
                  trand = exp ( - kn * rdiff )
                  sources_dn(n) = sources_dn(n) * trand
               endif
!  @@ End Robfix add 4 lines

               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = L_solutions_fine(n,j,k,q) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_sources_dn(n,k,q)  = L_sum * kn + L_extinction(N,q) * sum

!  @@ Robfix, add following 3/4 lines
                              if ( n.eq.NCrit ) then
                                 L_sources_dn(n,k,q) =  L_sources_dn(n,k,q) * trand - &
                                           sources_dn(n)* L_extinction(N,q) * rdiff
                              endif
!  @@ End Robfix add 3/4 lines

                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_func = L_solutions_fine(n,j,k,q)  * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_sources_dn(n,k,q)  = L_sum * kn

!  @@ Robfix, add following line
                              if ( n.eq.NCrit ) L_sources_dn(n,k,q) = L_sources_dn(n,k,q) * trand
!  @@ End Robfix add line

                           endif
                        enddo
                     endif
                  enddo
               else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = L_solutions_fine(n,j,kc,q) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j)
                     enddo
                     L_sources_dn(n,kc,q)  = L_sum * kn + L_extinction(N,q) * sum

!  @@ Robfix, add following 3/4 lines
                              if ( n.eq.NCrit ) then
                                 L_sources_dn(n,kc,q) =  L_sources_dn(n,kc,q) * trand - &
                                           sources_dn(n)* L_extinction(N,q) * rdiff
                              endif
!  @@ End Robfix add 3/4 lines

                  enddo
               endif
            endif

!  End layer loop and Nadir Enhanced PS case

         enddo
      endif

!  Enhanced PS: General case. Layer integrated Solar sources
!  =========================================================

      if ( do_enhanced_ps .and. .not. doNadir ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            cot_2 = cota(n-1) ; cot_1 = cota(n)
            cot_c = cot_1  ; if ( n.eq.NCrit ) cot_c = CotCrit
            kn = extinction(n) ;  ke = raycon * kn ; cons = raycon * ( cot_2 - cot_1 )
            tran_1 = kn * cons
            if ( tran_1 .lt. cutoff ) lostrans_dn(n) = exp ( - tran_1 )

!  @@ Robfix, add following 4 lines
            if ( n.eq.NCrit ) then
               consc = raycon * ( CotCrit - cot_1 )
               trand = exp ( - kn * consc )
            endif
!  @@ End Robfix add 4 lines

            if ( do_atmoswfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_extinction(n,q) * cons
                     L_lostrans_dn(n,q) = - L_lostau(q) * lostrans_dn(n)
                  enddo
               endif
            endif

!  Sources

            if ( layermask_dn(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  argum(j) = Raycon * ( cotfine(n,j) - cot_c )
                  tran(j)  = exp ( - kn * argum(j) )
                  func(j)  = solutions_fine(n,j) * csqfine(n,j) * tran(j)
                  sum      = sum + func(j) * wfine(n,j)
               enddo
               sources_dn(n) = sum * ke

!  @@ Robfix, add following line
               if ( n.eq.NCrit ) sources_dn(n) = sources_dn(n) * trand
!  @@ End Robfix add line

               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = L_solutions_fine(n,j,k,q) * csqfine(n,j) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_sources_dn(n,k,q)  = L_sum * ke + L_extinction(N,q) * Raycon * sum

!  @@ Robfix, add following 3/4 lines
                              if ( n.eq.NCrit ) then
                                 L_sources_dn(n,k,q) =  L_sources_dn(n,k,q) * trand - &
                                           sources_dn(n)* L_extinction(N,q) * consc
                              endif
!  @@ End Robfix add 3/4 lines

                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_func = L_solutions_fine(n,j,k,q) * csqfine(n,j) * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_sources_dn(n,k,q)  = L_sum * ke

!  @@ Robfix, add following line
                              if ( n.eq.NCrit ) L_sources_dn(n,k,q) = L_sources_dn(n,k,q) * trand
!  @@ End Robfix add line

                           endif
                        enddo
                     endif
                  enddo
               else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = L_solutions_fine(n,j,kc,q) * csqfine(n,j) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j)
                     enddo
                     L_sources_dn(n,kc,q)  = L_sum * ke + L_extinction(N,q) * Raycon * sum

!  @@ Robfix, add following 3/4 lines
                              if ( n.eq.NCrit ) then
                                 L_sources_dn(n,kc,q) =  L_sources_dn(n,kc,q) * trand - &
                                           sources_dn(n)* L_extinction(N,q) * consc
                              endif
!  @@ End Robfix add 3/4 lines

                  enddo
               endif
            endif

!  End layer loop and general Enhanced PS case

         enddo
      endif

!  Source function integration
!  ===========================

!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!  Check for updating the recursion

!  INTENSITY Main loop over all output optical depths
!          Cumulative source term will be saved

      NC = 0
      CUMSOURCE_DN(NC) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1
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

!  Profile Wfs

      if ( do_profilewfs ) then
         do k = 1, nlayers
            if ( Qvary(k) ) then
               do q = 1, Qnums(k)
                  L_CUMSOURCE(q) = ZERO
               enddo
               NSTART = 1
               NUT_PREV = NSTART - 1
               DO UTA = 1, N_USER_LEVELs
                  NUT    = USER_LEVELS(UTA)
                  DO N = NSTART, NUT
                     NC = N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(q) = L_SOURCES_DN(N,K,Q)         + &
                                L_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1) + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = L_SOURCES_DN(N,K,Q)         + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_JACOBIANS_DN(UTA,K,Q) = FLUX * L_CUMSOURCE(Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
                  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Column Wfs

      if ( do_columnwfs ) then
         do q = 1, n_columnwfs
            L_CUMSOURCE(q) = ZERO
         enddo
         NSTART = 1
         NUT_PREV = NSTART - 1
         DO UTA = 1, N_USER_LEVELS
            NUT    = USER_LEVELS(UTA)
            DO N = NSTART, NUT
               NC = N
               do q = 1, n_columnwfs
                  L_cumsource(q) = L_SOURCES_DN(N,KC,Q)         + &
                                L_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1) + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(Q)
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_JACOBIANS_DN(UTA,Q) = FLUX * L_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Finish

      return
end subroutine SS_Integral_LPLCLS_Plus_1_DN

!

subroutine SS_Integral_LPLCLS_Plus_1_UPDN &
   ( maxlayers, maxfinelayers, maxmoments_input,                                & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs,                             & ! Inputs (dimensioning)
     do_deltam_scaling, do_upwelling, do_dnwelling, do_regular_ps,              & ! Inputs (Flags)
     do_enhanced_ps, doNadir, do_columnwfs, do_profilewfs, do_surfacewfs,       & ! Inputs (Flags)
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,            & ! Inputs (Layers/Levels control)
     n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,               & ! Inputs (Jacobian control)
     NCrit, RadCrit, CotCrit, Mu0, Mu1, LegPoly_up, LegPoly_dn,                 & ! Inputs (Misc. Geometry)
     Raycon, radii, cota, xfine, wfine, csqfine, cotfine,                       & ! Inputs (Whole layer Geometry)
     sunpaths_up, ntraverse_up, sunpaths_dn, ntraverse_dn,                      & ! Inputs (Whole layer Geometry)
     sunpaths_fine_up, ntraverse_fine_up, sunpaths_fine_dn, ntraverse_fine_dn,  & ! Inputs (Fine  layer Geometry)
     Flux, reflec, extinction, deltaus, omega, truncfac, phasmoms,              & ! Inputs (Optical - Regular)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,       & ! Inputs (Optical - Linearized)
     intensity_up, intensity_db_up, LC_Jacobians_up, LP_Jacobians_up,           & ! Output
     LC_Jacobians_db_up, LP_Jacobians_db_up, LS_Jacobians_db_up,                & ! Output
     intensity_dn, LC_Jacobians_dn, LP_Jacobians_dn )                             ! Output

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Dimensions

      INTEGER, Intent(in) :: maxlayers
      INTEGER, Intent(in) :: maxfinelayers
      INTEGER, Intent(in) :: maxmoments_input
      INTEGER, Intent(in) :: max_user_levels
      INTEGER, Intent(in) :: max_atmoswfs
      INTEGER, Intent(in) :: max_surfacewfs

!  General flags

      LOGICAL, Intent(in) ::  DO_UPWELLING
      LOGICAL, Intent(in) ::  DO_DNWELLING
      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Jacobian Flags

      LOGICAL, Intent(in) ::  do_surfacewfs
      LOGICAL, Intent(in) ::  do_columnwfs
      LOGICAL, Intent(in) ::  do_profilewfs

!  Layer and Level Control Numbers, Number of Moments

      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )
      INTEGER, Intent(in) ::  NMOMENTS_INPUT

!  Jacobian control

      INTEGER, Intent(in) ::  n_columnwfs
      INTEGER, Intent(in) ::  n_surfacewfs
      LOGICAL, Intent(in) ::  Lvaryflags(maxlayers)
      INTEGER, Intent(in) ::  Lvarynums (maxlayers)
      LOGICAL, Intent(in) ::  Lvarymoms (maxlayers,max_atmoswfs)

!  Miscellaneous Geometrical inputs
!       Ray constant, Cotangents, Radii
!       Mu0 = cos(theta_boa), required for the surface term (both regular and enhanced)
!       Mu1 = cos(alpha_boa), required for the Regular PS only
!       Legendre Polynomials

      integer  , Intent(in)  :: NCrit
      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers), radii(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1, Mu0, RadCrit, CotCrit
      REAL(fpk), Intent(in)  :: LegPoly_up(0:maxmoments_input)
      REAL(fpk), Intent(in)  :: LegPoly_dn(0:maxmoments_input)

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  solar paths 

      integer  , Intent(in)  :: ntraverse_up  (0:maxlayers)
      real(fpk), Intent(in)  :: sunpaths_up   (0:maxlayers,maxlayers)
      integer  , Intent(in)  :: ntraverse_fine_up(maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: sunpaths_fine_up (maxlayers,maxlayers,maxfinelayers)
      integer  , Intent(in)  :: ntraverse_dn  (0:maxlayers)
      real(fpk), Intent(in)  :: sunpaths_dn   (0:maxlayers,maxlayers)
      integer  , Intent(in)  :: ntraverse_fine_dn(maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: sunpaths_fine_dn (maxlayers,maxlayers,maxfinelayers)

!  Optical inputs for Atmosphere, Surface reflectivity

      REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
      REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
      REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
      REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
      REAL(fpk), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )
      REAL(fpk), Intent(in) :: REFLEC, FLUX

!  Linearized optical inputs

      REAL(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT, max_atmoswfs )
      REAL(fpk), Intent(in) :: LS_REFLEC     ( max_surfacewfs)

!  outputs
!  -------

      real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels )
      real(fpk), Intent(Out)  :: intensity_db_up  ( max_user_levels )
      real(fpk), Intent(Out)  :: LC_Jacobians_up  ( max_user_levels, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_up  ( max_user_levels, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LC_Jacobians_db_up  ( max_user_levels, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_db_up  ( max_user_levels, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LS_Jacobians_db_up  ( max_user_levels, max_surfacewfs )

      real(fpk), Intent(Out)  :: intensity_dn     ( max_user_levels )
      real(fpk), Intent(Out)  :: LC_Jacobians_dn  ( max_user_levels, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dn  ( max_user_levels, maxlayers, max_atmoswfs )

   if ( do_upwelling  ) then
      call SS_Integral_LPLCLS_Plus_1_UP &
   ( maxlayers, maxfinelayers, maxmoments_input,                                & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs,                             & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,                 & ! Inputs (Flags - General  )
     do_columnwfs, do_profilewfs, do_surfacewfs,                                & ! Inputs (Flags - Jacobian )
     nlayers, nfinedivs, NCrit, nmoments_input, n_user_levels, user_levels,     & ! Inputs (Layers/Levels control)
     n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,               & ! Inputs (Jacobian control)
     Mu0, Mu1, LegPoly_up,  Raycon, cota, sunpaths_up, ntraverse_up,            & ! Inputs (Whole layer Geometry)
     xfine, wfine, csqfine, cotfine, sunpaths_fine_up, ntraverse_fine_up,       & ! Inputs (Fine  layer Geometry)
     Flux, reflec, extinction, deltaus, omega, truncfac, phasmoms,              & ! Inputs (Optical - Regular)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,       & ! Inputs (Optical - Linearized)
     intensity_up, intensity_db_up, LC_Jacobians_up, LP_Jacobians_up,           & ! Output
     LC_Jacobians_db_up, LP_Jacobians_db_up, LS_Jacobians_db_up )                 ! Output
   endif

   if ( do_dnwelling  ) then
      call SS_Integral_LPLCLS_Plus_1_DN &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels, max_atmoswfs,               & ! Inputs (dimensioning)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir, do_columnwfs, do_profilewfs,  & ! Inputs (Flags - Jacobian )
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,          & ! Inputs (Layers/Levels control)
     n_columnwfs, Lvaryflags, Lvarynums, Lvarymoms, NCrit, RadCrit, CotCrit,  & ! Inputs (Jacobian control)
     Mu1, LegPoly_dn, Raycon, radii, cota, sunpaths_dn, ntraverse_dn,         & ! Inputs (Whole layer Geometry)
     xfine, wfine, csqfine, cotfine, sunpaths_fine_dn, ntraverse_fine_dn,     & ! Inputs (Fine  layer Geometry)
     Flux, extinction, deltaus, omega, truncfac, phasmoms,                    & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,                & ! Inputs (Optical - Linearized)
     intensity_dn, LC_Jacobians_dn, LP_Jacobians_dn )                           ! Output
   endif

!  Finish

   return
end subroutine SS_Integral_LPLCLS_Plus_1_UPDN

subroutine DTE_Integral_LPLCLS_Plus_1_UP &
   ( maxlayers, maxfinelayers, max_user_levels, max_atmoswfs, max_surfacewfs, & ! Inputs (dimensioning)
     Do_Thermset, do_dmscaling, do_regular_ps, do_enhanced_ps,                & ! Inputs (Flags)
     doNadir, do_columnwfs, do_profilewfs, do_surfacewfs,                     & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,                          & ! Inputs (control output)
     n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums,                        & ! Inputs (Jacobian control)
     bb_input, surfbb, user_emissivity, extinction, deltaus, omega, truncfac, & ! Inputs (Optical)
     LS_user_emissivity, L_extinction, L_deltaus, L_omega, L_truncfac,        & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,                & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, tcom1, L_tcom1,                         & ! Outputs
     LC_Jacobians_dta_up, LP_Jacobians_dta_up, LS_Jacobians_dts )               ! Output

!  Stand alone routine for Upwelling Direct-thermal-emission (DTE)
!    computation of radiances and Jacobians. Can be derived from input Planck functions.

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
      INTEGER, Intent(in) :: max_atmoswfs
      INTEGER, Intent(in) :: max_surfacewfs

!  Thermal setup flag (for TCOM1)

      LOGICAL, Intent(inout) ::  Do_Thermset

!  flags

      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DMSCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Jacobian Flags

      LOGICAL, Intent(in) ::  do_surfacewfs
      LOGICAL, Intent(in) ::  do_columnwfs
      LOGICAL, Intent(in) ::  do_profilewfs

!  Layer and Level Control Numbers, Number of Moments

      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

      INTEGER, Intent(in) ::  n_columnwfs
      INTEGER, Intent(in) ::  n_surfacewfs
      LOGICAL, Intent(in) ::  Lvaryflags(maxlayers)
      INTEGER, Intent(in) ::  Lvarynums (maxlayers)

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

!  Linearized optical inputs

      REAL(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: LS_USER_EMISSIVITY  ( max_surfacewfs)

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

!  Radiances

      real(fpk), Intent(Out)  :: intensity_dta_up ( max_user_levels )
      real(fpk), Intent(Out)  :: intensity_dts    ( max_user_levels )

      real(fpk), Intent(Out)  :: LC_Jacobians_dta_up  ( max_user_levels, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dta_up  ( max_user_levels, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LS_Jacobians_dts     ( max_user_levels, max_surfacewfs )

!  Thermal setup

      real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)
      real(fpk), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)

!  LOCAL
!  -----

!  Local solutions (enhanced_ps case)

      real(fpk)  :: Solutions_fine (maxlayers,maxfinelayers)
      real(fpk)  :: Wtrans_fine    (maxlayers,maxfinelayers)
      real(fpk)  :: L_Solutions_fine (maxlayers,maxfinelayers,max_atmoswfs)
      real(fpk)  :: L_Wtrans_fine    (maxlayers,maxfinelayers,max_atmoswfs)

!  Source function integration results

      real(fpk)  :: sources_up       ( maxlayers )
      real(fpk)  :: lostrans_up      ( maxlayers )
      real(fpk)  :: cumsource_up     ( 0:maxlayers )

      real(fpk)  :: L_lostrans_up    ( maxlayers,max_atmoswfs )
      real(fpk)  :: L_cumsource      ( max_atmoswfs )
      real(fpk)  :: L_sources_up     ( maxlayers,max_atmoswfs )
      real(fpk)  :: LS_cumsource     ( max_surfacewfs )

!  Help

      integer    :: n, k, kc, j, nj, q, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
      logical    :: do_atmoswfs, layermask_up(maxlayers), Qvary(maxlayers)

      real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers)
      real(fpk)  :: help, sum, kn, ke, xjkn, cot_1, cot_2, arg, tcw
      real(fpk)  :: L_help, L_sum, L_tran, L_kn
      real(fpk)  :: cumsource_dste, t_mult_up(0:2), L_t_mult_up(0:2), thermcoeffs(2)
      real(fpk)  :: tms, L_tms(max_atmoswfs), lostau

      real(fpk), parameter  :: cutoff = 88.0d0
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

      INTENSITY_dta_up = zero
      INTENSITY_dts    = zero

      LP_JACOBIANS_dta_UP = zero
      LC_JACOBIANS_dta_UP = zero
      LS_JACOBIANS_dts    = zero

!  Zero local sources

      lostrans_up   = zero  ; sources_up   = zero ; cumsource_up = zero
      L_lostrans_up = zero  ; L_sources_up = zero

!  Bookkeeping

      NUT = USER_LEVELS(1) + 1
      LAYERMASK_UP = .false.
      LAYERMASK_UP(NUT:NLAYERS) = .true.

      nstart = nlayers ; if (Ncrit.ne.0) nstart = nCrit
      if (nstart.lt.nlayers) LAYERMASK_UP(nstart+1:nlayers) = .false.

      do_atmoswfs = do_profilewfs .or. do_columnwfs
      Qvary = .false. ; QNums = 0
      if ( do_profilewfs ) then
         Qvary(1:nlayers) = Lvaryflags(1:nlayers)
         QNums(1:nlayers) = Lvarynums (1:nlayers)
      else if ( do_columnwfs ) then
         Qvary(1:nlayers) = .true.
         QNums(1:nlayers) =  n_columnwfs
      endif
      kc = 0

!  Thermal setup factors and linearizations
!     TMS, Initial set of thermal coefficients and TCOM1 variable

      if ( do_Thermset ) then
         tcom1 = zero
         do n = 1, nlayers
            tms = one - omega(n) 
            if ( do_dmscaling ) then
               help = one - truncfac(n) * omega(n)
               tms = tms / help
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
                     L_tms(q) = - ( L_omega(n,q) - tms * L_help ) / help
                  enddo
               endif
            endif
            thermcoeffs(1)  = bb_input(n-1)
            thermcoeffs(2)  = (bb_input(n)-bb_input(n-1)) / deltaus(n)
            tcom1(n,1) = thermcoeffs(1) * tms
            tcom1(n,2) = thermcoeffs(2) * tms
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_tcom1(n,1,q) = thermcoeffs(1) * L_tms(q)
                  L_tcom1(n,2,q) = thermcoeffs(2) * L_tms(q)
               enddo
            endif
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
            if ( Qvary(n) ) L_lostrans_up(n,1:Qnums(n)) = -lostrans_up(n) * L_deltaus(n,1:Qnums(n))
            if ( layermask_up(n) ) then
              t_mult_up(1:2) = tcom1(n,1:2)
              sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
              t_mult_up(0) = - sum
              sources_up(n) = t_mult_up(1)           ! Check this formula
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                   L_t_mult_up(1:2) = L_tcom1(n,1:2,q)
                   L_sum = L_t_mult_up(1) + t_mult_up(2) * L_deltaus(n,q) + L_t_mult_up(2) * deltaus(n)
                   L_t_mult_up(0) = - sum
                   L_sources_up(n,q) = L_t_mult_up(1)           ! Check this formula
                enddo
              endif
            endif
          enddo
        else
          DO n = 1, nlayers
            lostau = deltaus(n) / Mu1
            if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
            if ( Qvary(n) ) L_lostrans_up(n,1:Qnums(n)) = - lostrans_up(n) * L_deltaus(n,1:Qnums(n)) / Mu1
            if ( layermask_up(n) ) then
              t_mult_up(2) = tcom1(n,2)
              t_mult_up(1) = tcom1(n,1) + t_mult_up(2) * Mu1
              sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
              t_mult_up(0) = - sum
              sources_up(n) = t_mult_up(0) * lostrans_up(n)  + t_mult_up(1)
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                   L_t_mult_up(2) = L_tcom1(n,2,q)
                   L_t_mult_up(1) = L_tcom1(n,1,q) + L_t_mult_up(2) * Mu1
                   L_sum = L_t_mult_up(1) + t_mult_up(2) * L_deltaus(n,q) + L_t_mult_up(2) * deltaus(n)
                   L_t_mult_up(0) = - sum
                   L_sources_up(n,q) = L_t_mult_up(0) *   lostrans_up(n)   + &
                                         t_mult_up(0) * L_lostrans_up(n,q) + L_t_mult_up(1)
                enddo
              endif
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
               if ( Qvary(n) ) L_lostrans_up(n,1:Qnums(n)) = -lostrans_up(n) * L_deltaus(n,1:Qnums(n))
               do j = 1, nj
                  argum(j) = xfine(n,j) ; xjkn = argum(j) * kn
                  tran(j)  = exp ( -xjkn )
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = kn * tran(j) * wfine(n,j)
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_kn = L_extinction(n,q) ; L_tran = argum(j) * L_kn
                        L_solutions_fine(n,j,q) = L_tcom1(n,1,q) + xjkn * L_tcom1(n,2,q) + L_tran * tcom1(n,2)
                        L_wtrans_fine(n,j,q)    = ( L_kn - kn * L_tran ) * tran(j) * wfine(n,j)
                     enddo
                  endif
               enddo
            else if ( .not. doNadir .and. layermask_up(n) ) then
               cot_2 = cota(n-1) ; cot_1 = cota(n)
               ke = raycon * kn  ; arg = raycon * ( cot_2 - cot_1 ) ; lostau = kn * arg
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               if ( Qvary(n) ) L_lostrans_up(n,1:Qnums(n)) = - arg * lostrans_up(n) * L_extinction(n,1:Qnums(n))
               do j = 1, nj
                  argum(j) = Raycon * ( cot_2 - cotfine(n,j) )  ; xjkn = xfine(n,j) * kn
                  tran(j)  = exp ( - kn * argum(j) ) ; tcw = tran(j) * csqfine(n,j) * wfine(n,j)
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = kn * tcw
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_kn = L_extinction(n,q) ; L_tran = argum(j) * L_kn
                        L_solutions_fine(n,j,q) = L_tcom1(n,1,q) + xfine(n,j)*(kn*L_tcom1(n,2,q) + L_kn*tcom1(n,2))
                        L_wtrans_fine(n,j,q)    = ( L_kn - kn * L_tran ) * tcw
                     enddo
                  endif
               enddo
            endif
            sources_up(n) = dot_product(solutions_fine(n,1:nj),wtrans_fine(n,1:nj))
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_sources_up(n,q) = dot_product(L_solutions_fine(n,1:nj,q),  wtrans_fine(n,1:nj)) + &
                                      dot_product(  solutions_fine(n,1:nj),  L_wtrans_fine(n,1:nj,q))
               enddo
            endif
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

!  Surface WFs

      if ( do_surfacewfs ) then
         do q = 1, n_surfacewfs
            LS_cumsource(q) = SURFBB * LS_USER_EMISSIVITY(q)
         enddo
         NSTART = NLAYERS
         NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_surfacewfs
                  LS_cumsource(q) = LOSTRANS_UP(N) * LS_CUMSOURCE(Q)
               enddo
            ENDDO
            do q = 1, n_surfacewfs
               LS_Jacobians_dts(UTA,Q) = LS_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Profile Wfs (atmospheric term)

      if ( do_profilewfs ) then
         do k = 1, nlayers
            if ( Qvary(k) ) then
               L_CUMSOURCE = zero
               NSTART = NLAYERS
               NUT_PREV = NSTART + 1
               DO UTA = N_USER_LEVELS, 1, -1
                  NUT    = USER_LEVELS(UTA) + 1
                  DO N = NSTART, NUT, -1
                     NC = NLAYERS + 1 - N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(q) = L_SOURCES_UP(N,Q)           + &
                                L_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1) + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LOSTRANS_UP(N) * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_Jacobians_dta_up(UTA,K,Q) = L_CUMSOURCE(Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Column Wfs (Atmospheric term)

      if ( do_columnwfs ) then
         L_CUMSOURCE = zero
         NSTART = NLAYERS
         NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(q) = L_SOURCES_UP(N,Q)                     + &
                                 L_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1) + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_Jacobians_dta_up(UTA,Q) = L_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Finish

      return
end subroutine DTE_Integral_LPLCLS_Plus_1_UP

!

subroutine DTE_Integral_LPLCLS_Plus_1_DN &
   ( maxlayers, maxfinelayers, max_user_levels, max_atmoswfs,                 & ! Inputs (dimensioning)
     Do_Thermset, do_dmscaling, do_regular_ps, do_enhanced_ps,                & ! Inputs (Flags)
     doNadir, do_columnwfs, do_profilewfs, nlayers, nfinedivs,                & ! Inputs (Flags/Numbers)
     n_user_levels, user_levels,  n_columnwfs, Lvaryflags, Lvarynums,         & ! Inputs (Numbers,Jacobian control)
     bb_input, extinction, deltaus, omega, truncfac,                          & ! Inputs (Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac,                            & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Radcrit, CotCrit, Raycon, cota,                              & ! Inputs (Geometry)
     radii, xfine, wfine, csqfine, cotfine,                                   & ! Inputs (Geometry)
     intensity_dta_dn, tcom1, L_tcom1, LC_Jacobians_dta_dn, LP_Jacobians_dta_dn )    ! Output

!  Stand alone routine for Downwelling Direct-thermal-emission (DTE)
!    computation of radiances and Jacobians. Can be derived from input Planck functions.

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
      INTEGER, Intent(in) :: max_atmoswfs

!  Thermal setup flag (for TCOM1)

      LOGICAL, Intent(inout) ::  Do_Thermset

!  flags

      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DMSCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Jacobian Flags

      LOGICAL, Intent(in) ::  do_columnwfs
      LOGICAL, Intent(in) ::  do_profilewfs

!  Numbers

      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

      INTEGER, Intent(in) ::  n_columnwfs
      LOGICAL, Intent(in) ::  Lvaryflags(maxlayers)
      INTEGER, Intent(in) ::  Lvarynums (maxlayers)

!  optical inputs
!  --------------

!  Atmosphere extinction and deltaus

      REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
      REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
      REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
      REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric BB functions

      REAL(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Linearized optical inputs

      REAL(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      integer  , Intent(in)  :: NCrit
      real(fpk), Intent(in)  :: Radcrit, CotCrit
      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers), radii(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

!  Radiances

      real(fpk), Intent(Out)  :: intensity_dta_dn ( max_user_levels )

      real(fpk), Intent(Out)  :: LC_Jacobians_dta_dn  ( max_user_levels, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dta_dn  ( max_user_levels, maxlayers, max_atmoswfs )

!  Thermal setup

      real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)
      real(fpk), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)

!  LOCAL
!  -----

!  Local solutions (enhanced_ps case)

      real(fpk)  :: Solutions_fine (maxlayers,maxfinelayers)
      real(fpk)  :: Wtrans_fine    (maxlayers,maxfinelayers)
      real(fpk)  :: L_Solutions_fine (maxlayers,maxfinelayers,max_atmoswfs)
      real(fpk)  :: L_Wtrans_fine    (maxlayers,maxfinelayers,max_atmoswfs)

!  Source function integration results

      real(fpk)  :: sources_dn       ( maxlayers )
      real(fpk)  :: lostrans_dn      ( maxlayers )
      real(fpk)  :: cumsource_dn     ( 0:maxlayers )

      real(fpk)  :: L_lostrans_dn    ( maxlayers,max_atmoswfs )
      real(fpk)  :: L_cumsource      ( max_atmoswfs )
      real(fpk)  :: L_sources_dn     ( maxlayers,max_atmoswfs )

!  Help

      integer    :: n, k, kc, j, nj, q, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
      logical    :: do_atmoswfs, layermask_dn(maxlayers), Qvary(maxlayers)

      real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers)
      real(fpk)  :: help, sum, kn, ke, xjkn, cot_1, cot_2, cot_c, rdiff, arg0, arg, trand, tcw
      real(fpk)  :: L_help, L_sum, L_tran, L_kn, L_trand
      real(fpk)  :: t_mult_dn(0:2), L_t_mult_dn(0:2), thermcoeffs(2)
      real(fpk)  :: tms, L_tms(max_atmoswfs), lostau

      real(fpk), parameter  :: cutoff = 88.0d0
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

      INTENSITY_dta_dn = zero

      LP_JACOBIANS_dta_dn = zero
      LC_JACOBIANS_dta_dn = zero

!  Zero local sources

      lostrans_dn   = zero  ; sources_dn   = zero ; cumsource_dn = zero
      L_lostrans_dn = zero  ; L_sources_dn = zero

      trand = one ; L_trand = zero

!  Bookkeeping

      NUT = USER_LEVELS(1) + 1
      IF ( NUT > NLAYERS ) NUT = NLAYERS
      LAYERMASK_DN = .false.
      LAYERMASK_DN(1:NUT) = .true.

      nstart = nlayers ; if (Ncrit.ne.0) nstart = nCrit
      if (nstart.lt.nlayers) LAYERMASK_DN(nstart+1:nlayers) = .false.

      do_atmoswfs = do_profilewfs .or. do_columnwfs
      Qvary = .false. ; QNums = 0
      if ( do_profilewfs ) then
         Qvary(1:nlayers) = Lvaryflags(1:nlayers)
         QNums(1:nlayers) = Lvarynums (1:nlayers)
      else if ( do_columnwfs ) then
         Qvary(1:nlayers) = .true.
         QNums(1:nlayers) =  n_columnwfs
      endif
      kc = 0

!  Thermal setup factors and linearizations
!     TMS, Initial set of thermal coefficients and TCOM1 variable

      if ( do_Thermset ) then
         tcom1 = zero
         do n = 1, nlayers
            tms = one - omega(n) 
            if ( do_dmscaling ) then
               help = one - truncfac(n) * omega(n)
               tms = tms / help
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
                     L_tms(q) = - ( L_omega(n,q) - tms * L_help ) / help
                  enddo
               endif
            endif
            thermcoeffs(1)  = bb_input(n-1)
            thermcoeffs(2)  = (bb_input(n)-bb_input(n-1)) / deltaus(n)
            tcom1(n,1) = thermcoeffs(1) * tms
            tcom1(n,2) = thermcoeffs(2) * tms
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_tcom1(n,1,q) = thermcoeffs(1) * L_tms(q)
                  L_tcom1(n,2,q) = thermcoeffs(2) * L_tms(q)
               enddo
            endif
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
            if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = -lostrans_dn(n) * L_deltaus(n,1:Qnums(n))
            if ( layermask_dn(n) ) then
              t_mult_dn(1:2) = tcom1(n,1:2)
              t_mult_dn(0)   = - t_mult_dn(1)
              sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
              sources_dn(n)  = sum
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                   L_t_mult_dn(1:2) = L_tcom1(n,1:2,q)
                   L_t_mult_dn(0)   = - L_t_mult_dn(1)
                   L_sum = L_t_mult_dn(1) + t_mult_dn(2) * L_deltaus(n,q) + L_t_mult_dn(2) * deltaus(n)
                   L_sources_dn(n,q) = L_sum          ! Check this formula
                enddo
              endif
            endif
          enddo
        else
          DO n = 1, nlayers
            lostau = deltaus(n) / Mu1
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
            if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = - lostrans_dn(n) * L_deltaus(n,1:Qnums(n)) / Mu1
            if ( layermask_dn(n) ) then
              t_mult_dn(2)   = tcom1(n,2)
              t_mult_dn(1)   = tcom1(n,1) - t_mult_dn(2) * Mu1
              t_mult_dn(0)   = - t_mult_dn(1)
              sources_dn(n)  = t_mult_dn(0) * lostrans_dn(n)
              sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
              sources_dn(n)  = sources_dn(n) + sum
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                   L_t_mult_dn(2) = L_tcom1(n,2,q)
                   L_t_mult_dn(1) = L_tcom1(n,1,q) - L_t_mult_dn(2) * Mu1
                   L_t_mult_dn(0) = - L_t_mult_dn(1)
                   L_sources_dn(n,q)  = L_t_mult_dn(0) * lostrans_dn(n) + t_mult_dn(0) * L_lostrans_dn(n,q)
                   L_sum = L_t_mult_dn(1) + t_mult_dn(2) * L_deltaus(n,q) + L_t_mult_dn(2) * deltaus(n)
                   L_sources_dn(n,q) = L_sources_dn(n,q) + L_sum
                enddo
              endif
            endif
          enddo
        endif
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

      if ( do_enhanced_ps ) then
         do n = nlayers, 1, -1
            kn = extinction(n) ; nj = nfinedivs(n)
            if (  doNadir  .and. layermask_dn(n) ) then
               rdiff = radii(n-1) - radii(n) ; if ( n.eq.NCrit) rdiff = radii(n-1) - RadCrit
               trand = one ; if ( n.eq.NCrit) trand = exp ( -kn * (RadCrit -radii(n) ) )
               lostau = deltaus(n)
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = -lostrans_dn(n) * L_deltaus(n,1:Qnums(n))
               do j = 1, nj
                  argum(j) = rdiff - xfine(n,j) ; xjkn = xfine(n,j) * kn
                  tran(j)  = exp ( -xjkn )
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = kn * tran(j) * wfine(n,j)
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_kn = L_extinction(n,q) ; L_tran = argum(j) * L_kn
                        L_solutions_fine(n,j,q) = L_tcom1(n,1,q) + xjkn * L_tcom1(n,2,q) + L_tran * tcom1(n,2)
                        L_wtrans_fine(n,j,q)    = ( L_kn - kn * L_tran ) * tran(j) * wfine(n,j)
                     enddo
                  endif
               enddo
            else if ( .not. doNadir .and. layermask_dn(n) ) then
               cot_2 = cota(n-1) ; cot_1 = cota(n)
               cot_c = cot_1     ; if ( n.eq.NCrit ) cot_c = CotCrit
               ke = raycon * kn  ; arg = raycon * ( cot_2 - cot_1 ) ; lostau = kn * arg
               if ( n.eq.NCrit ) then
                  arg0 = Raycon * ( CotCrit - cot_1 ) ; trand = exp ( - kn * arg0 )
               endif
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = - arg * lostrans_dn(n) * L_extinction(n,1:Qnums(n))
               do j = 1, nj
                  argum(j) = raycon * ( cotfine(n,j) - cot_1 )  ; xjkn = xfine(n,j) * kn
                  tran(j)  = exp ( - kn * argum(j) ) ; tcw = tran(j) * csqfine(n,j) * wfine(n,j)
                  solutions_fine(n,j) = tcom1(n,1) + xjkn * tcom1(n,2)
                  wtrans_fine(n,j)    = kn * tcw
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_kn = L_extinction(n,q) ; L_tran = argum(j) * L_kn
                        L_solutions_fine(n,j,q) = L_tcom1(n,1,q) + xfine(n,j)*(kn*L_tcom1(n,2,q) + L_kn*tcom1(n,2))
                        L_wtrans_fine(n,j,q)    = ( L_kn - kn * L_tran ) * tcw
                     enddo
                  endif
               enddo
            endif
            sources_dn(n) = dot_product(solutions_fine(n,1:nj),wtrans_fine(n,1:nj))
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_sources_dn(n,q) = dot_product(L_solutions_fine(n,1:nj,q),  wtrans_fine(n,1:nj)) + &
                                      dot_product(  solutions_fine(n,1:nj),  L_wtrans_fine(n,1:nj,q) )
                  if ( n.eq.NCrit ) then
                     L_trand =  - L_extinction(n,q) * trand * arg0
                     L_sources_dn(n,q) = trand  * L_sources_dn(n,q) + L_trand * sources_dn(n)    !@@ Robfix
                  endif
               enddo
            endif
            if ( n.eq.NCrit ) sources_dn(n) = sources_dn(n) * trand         !@@ Robfix
         enddo
      endif

!  Source function integration
!  ===========================

!  start recursion ( For DSTE term, Use surface emissivity )

      NC =  0
      CUMSOURCE_DN(NC) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

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

!  Profile Wfs (atmospheric term)

      if ( do_profilewfs ) then
         do k = 1, nlayers
            if ( Qvary(k) ) then
               L_CUMSOURCE = zero
               NSTART = 1
               NUT_PREV = NSTART - 1
               DO UTA = 1, N_USER_LEVELS
                  NUT    = USER_LEVELS(UTA)
                  DO N = NSTART, NUT
                     NC = N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(q) = L_SOURCES_DN(N,Q)           + &
                                L_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1) + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LOSTRANS_DN(N) * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_Jacobians_dta_DN(UTA,K,Q) = L_CUMSOURCE(Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
                  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Column Wfs (Atmospheric term)

      if ( do_columnwfs ) then
         L_CUMSOURCE = zero
         NSTART = 1
         NUT_PREV = NSTART - 1
         DO UTA = 1, N_USER_LEVELS
            NUT    = USER_LEVELS(UTA)
            DO N = NSTART, NUT
               NC = N
               do q = 1, n_columnwfs
                  L_cumsource(q) = L_SOURCES_DN(N,Q)                     + &
                                 L_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1) + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(Q)
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_Jacobians_dta_dn(UTA,Q) = L_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Finish

      return
end subroutine DTE_Integral_LPLCLS_Plus_1_DN

!

subroutine DTE_Integral_LPLCLS_Plus_1_UPDN &
   ( maxlayers, maxfinelayers, max_user_levels, max_atmoswfs, max_surfacewfs, & ! Inputs (dimensioning)
     Do_Thermset, do_upwelling, do_dnwelling, do_dmscaling, do_regular_ps,    & ! Inputs (Flags)
     do_enhanced_ps, doNadir, do_columnwfs, do_profilewfs, do_surfacewfs,     & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,                          & ! Inputs (control output)
     n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums,                        & ! Inputs (Jacobian control)
     bb_input, surfbb, user_emissivity, extinction, deltaus, omega, truncfac, & ! Inputs (Optical)
     LS_user_emissivity, L_extinction, L_deltaus, L_omega, L_truncfac,        & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Radcrit, CotCrit, Raycon, cota,                              & ! Inputs (Geometry)
     radii, xfine, wfine, csqfine, cotfine,                                   & ! Inputs (Geometry)
     intensity_dta_up, intensity_dta_dn, intensity_dts, tcom1, L_tcom1,       & ! Outputs
     LC_Jacobians_dta_up, LP_Jacobians_dta_up, LS_Jacobians_dts,              & ! Output
     LC_Jacobians_dta_dn, LP_Jacobians_dta_dn )                                 ! Output

!  Stand alone routine for Upwelling and downwelling Direct-thermal-emission (DTE)
!    computation of radiances and Jacobians. Can be derived from input Planck functions.

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
      INTEGER, Intent(in) :: max_atmoswfs
      INTEGER, Intent(in) :: max_surfacewfs

!  Thermal setup flag (for TCOM1)

      LOGICAL, Intent(inout) ::  Do_Thermset

!  flags

      LOGICAL, Intent(in) ::  DO_UPWELLING
      LOGICAL, Intent(in) ::  DO_DNWELLING
      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DMSCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Jacobian Flags

      LOGICAL, Intent(in) ::  do_surfacewfs
      LOGICAL, Intent(in) ::  do_columnwfs
      LOGICAL, Intent(in) ::  do_profilewfs

!  Numbers

      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

      INTEGER, Intent(in) ::  n_columnwfs
      INTEGER, Intent(in) ::  n_surfacewfs
      LOGICAL, Intent(in) ::  Lvaryflags(maxlayers)
      INTEGER, Intent(in) ::  Lvarynums (maxlayers)

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

!  Linearized optical inputs

      REAL(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: LS_USER_EMISSIVITY  ( max_surfacewfs)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!    Mu1 = cos(alpha_boa), required for the Regular PS only

      integer  , Intent(in)  :: NCrit
      real(fpk), Intent(in)  :: Radcrit, CotCrit
      real(fpk), Intent(in)  :: Raycon, cota(0:maxlayers), radii(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1

!  LOS Quadratures for Enhanced PS

      real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers)
      real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers)

!  outputs
!  -------

!  Radiances

      real(fpk), Intent(Out)  :: intensity_dta_dn ( max_user_levels )
      real(fpk), Intent(Out)  :: intensity_dta_up ( max_user_levels )
      real(fpk), Intent(Out)  :: intensity_dts    ( max_user_levels )

      real(fpk), Intent(Out)  :: LC_Jacobians_dta_dn  ( max_user_levels, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dta_dn  ( max_user_levels, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LC_Jacobians_dta_up  ( max_user_levels, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dta_up  ( max_user_levels, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LS_Jacobians_dts     ( max_user_levels, max_surfacewfs )

!  Thermal setup

      real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)
      real(fpk), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)


   if ( do_upwelling ) then
      call DTE_Integral_LPLCLS_Plus_1_UP &
   ( maxlayers, maxfinelayers, max_user_levels, max_atmoswfs, max_surfacewfs, & ! Inputs (dimensioning)
     Do_Thermset, do_dmscaling, do_regular_ps, do_enhanced_ps,                & ! Inputs (Flags)
     doNadir, do_columnwfs, do_profilewfs, do_surfacewfs,                     & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,                          & ! Inputs (control output)
     n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums,                        & ! Inputs (Jacobian control)
     bb_input, surfbb, user_emissivity, extinction, deltaus, omega, truncfac, & ! Inputs (Optical)
     LS_user_emissivity, L_extinction, L_deltaus, L_omega, L_truncfac,        & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,                & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, tcom1, L_tcom1,                         & ! Outputs
     LC_Jacobians_dta_up, LP_Jacobians_dta_up, LS_Jacobians_dts )               ! Output
   endif

   if ( do_dnwelling ) then
      call DTE_Integral_LPLCLS_Plus_1_DN &
   ( maxlayers, maxfinelayers, max_user_levels, max_atmoswfs,                 & ! Inputs (dimensioning)
     Do_Thermset, do_dmscaling, do_regular_ps, do_enhanced_ps,                & ! Inputs (Flags)
     doNadir, do_columnwfs, do_profilewfs, nlayers, nfinedivs,                & ! Inputs (Flags/Numbers)
     n_user_levels, user_levels,  n_columnwfs, Lvaryflags, Lvarynums,         & ! Inputs (Numbers,Jacobian control)
     bb_input, extinction, deltaus, omega, truncfac,                          & ! Inputs (Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac,                            & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Radcrit, CotCrit, Raycon, cota,                              & ! Inputs (Geometry)
     radii, xfine, wfine, csqfine, cotfine,                                   & ! Inputs (Geometry)
     intensity_dta_dn, tcom1, L_tcom1, LC_Jacobians_dta_dn, LP_Jacobians_dta_dn )    ! Output
   endif

!  Finish

   return
end subroutine DTE_Integral_LPLCLS_Plus_1_UPDN

!  End module

end module FO_Scalar_RTCalcs_Plus_m


