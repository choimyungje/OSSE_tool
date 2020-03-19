
module FO_Vector_RTCalcs_Plus_m

!  For given wavelength, routine calculates First Order upwelling and downwelling
!  Stokes vectors and any number of LCLPLS Jacobians (column/profile/surface)
!     (1) For the Atmospheric Single-scatter and Surface Direct-Beam solar sources,
!     (2) For the Atmospheric and Surface Thermal emission Direct solutions

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS (LOS-path sphericity) calculations
!  This will perform Regular-PS  (average-secant pseudo-spherical) calculations

!  This is Version 1.2, without Partials. Code is stand alone with no dependencies.
!    01 December 2011, R. Spurr, RT Solutions Inc.
!    02 February 2012, R. Spurr, RT Solutions Inc.
!    13 February 2012, R. Spurr, RT Solutions Inc.
!    21 June     2012, R. Spurr, RT Solutions Inc.

!  For Solar sources, the subroutines are
!       SSV_Integral_Plus_1_UP   (Upwelling only)
!       SSV_Integral_Plus_1_DN   (Downwelling only)
!       SSV_Integral_Plus_1_UPDN (Upwelling and Downwelling)

!  For Thermal Emission sources, the subroutines are
!       DTEV_Integral_Plus_1_UP   (Upwelling only)
!       DTEV_Integral_Plus_1_DN   (Downwelling only)
!       DTEV_Integral_Plus_1_UPDN (Upwelling and Downwelling)

!  All subroutines public

public

contains


subroutine SSV_Integral_Plus_1_UP &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels, do_sunlight,  & ! Inputs (dimensioning/flag)
     max_atmoswfs, max_surfacewfs, do_columnwfs, do_profilewfs, do_surfacewfs,  & ! Inputs (dimensioning/flags)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, do_lambertian, doNadir,  & ! Inputs (Flags - General  )
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,            & ! Inputs (Layers/Levels control)
     nstokes, n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,      & ! Inputs (Jacobian control)
     reflec, extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,     & ! Inputs (Optical)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat,       & ! Inputs (Optical - Linearized)
     Mu0, Mu1, GenSpher, Rotations, NCrit, xfine, wfine, csqfine, cotfine,      & ! Inputs (Geometry)
     Raycon, cota, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,          & ! Inputs (Geometry)
     Stokes_up, Stokes_db_up, LC_Jacobians_up, LP_Jacobians_up,                 & ! Output
     LC_Jacobians_db_up, LP_Jacobians_db_up, LS_Jacobians_db_up )                 ! Output

!  Stand alone routine for SS field (Radiances and Jacobians) with Solar sources alone
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

      INTEGER, Intent(in) :: max_atmoswfs
      INTEGER, Intent(in) :: max_surfacewfs

!  General flags

      LOGICAL, Intent(in) ::  DO_SUNLIGHT

      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DO_LAMBERTIAN
      LOGICAL, Intent(in) ::  DONADIR

!  Jacobian Flags

      LOGICAL, Intent(in) ::  do_surfacewfs
      LOGICAL, Intent(in) ::  do_columnwfs
      LOGICAL, Intent(in) ::  do_profilewfs

!  Numbers

      INTEGER, Intent(in) ::  NSTOKES
      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  NMOMENTS_INPUT

      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

      INTEGER, Intent(in) ::  n_columnwfs
      INTEGER, Intent(in) ::  n_surfacewfs
      LOGICAL, Intent(in) ::  Lvaryflags(maxlayers)
      INTEGER, Intent(in) ::  Lvarynums (maxlayers)
      LOGICAL, Intent(in) ::  Lvarymoms (maxlayers,max_atmoswfs)

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

!  Linearized optical inputs

      REAL(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4, max_atmoswfs )
      REAL(fpk), Intent(in) :: LS_REFLEC     ( 4, 4, max_surfacewfs )

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

      real(fpk), Intent(Out)  :: Stokes_up     ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: Stokes_db_up  ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: LC_Jacobians_up  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_up  ( max_user_levels, 4, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LC_Jacobians_db_up  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_db_up  ( max_user_levels, 4, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LS_Jacobians_db_up  ( max_user_levels, 4, max_surfacewfs )

!  LOCAL
!  -----

!  Attenuations

      real(fpk)  :: suntau              (0:maxlayers)
      real(fpk)  :: attenuations        (0:maxlayers)
      real(fpk)  :: attenuations_fine   (maxlayers,maxfinelayers)

      real(fpk)  :: L_suntau            (0:maxlayers,0:maxlayers,max_atmoswfs)
      real(fpk)  :: L_attenuations      (0:maxlayers,0:maxlayers,max_atmoswfs)
      real(fpk)  :: L_attenuations_fine (maxlayers,maxfinelayers,0:maxlayers,max_atmoswfs)

!  Scattering

      real(fpk)  :: tms            (maxlayers)
      real(fpk)  :: exactscat_up   (maxlayers,4,4)

      real(fpk)  :: L_tms          (maxlayers,max_atmoswfs)
      real(fpk)  :: L_exactscat_up (maxlayers,4,4, max_atmoswfs)

!  Source function integration results

      real(fpk)  :: sources_up       ( maxlayers, 4 )
      real(fpk)  :: lostrans_up      ( maxlayers )
      real(fpk)  :: multiplier_up    ( maxlayers )

      real(fpk)  :: L_sources_up      ( maxlayers,4,0:maxlayers,max_atmoswfs )
      real(fpk)  :: L_lostrans_up     ( maxlayers, max_atmoswfs )
      real(fpk)  :: L_multiplier_up   ( maxlayers,0:maxlayers,max_atmoswfs )

!  Local cumulative source terms

      real(fpk)  :: cumsource_db      ( 0:maxlayers, 4 )
      real(fpk)  :: cumsource_up      ( 0:maxlayers, 4 )

      real(fpk)  :: L_cumsource       ( 4, max_atmoswfs )
      real(fpk)  :: LS_cumsource      ( 4, max_surfacewfs )

!  Help

      integer    :: n, ns, k, kc, j, q, L, o1, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
      logical    :: do_atmoswfs, layermask_up(maxlayers), Qvary(maxlayers)

      real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers), func(maxfinelayers)

      real(fpk)  :: cons, help, sum, tran_1, kn, ke, factor1, factor2, factor3, m4, m4a, rhelp(4), shelp(4)
      real(fpk)  :: cot_1, cot_2, L_help, L_sum, L_tran, L_func, L_factor1, L_factor2, sumd, L_sumd
      real(fpk)  :: lostau, L_lostau(max_atmoswfs), fmat(6), L_fmat(6), sum23, dif23, L_Shelp
      real(fpk)  :: help3c1, help3s1, help4c1, help4s1

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      STOKES_UP       = zero
      LP_JACOBIANS_UP = zero
      LC_JACOBIANS_UP = zero

      STOKES_DB_UP       = zero
      LP_JACOBIANS_DB_UP = zero
      LC_JACOBIANS_DB_UP = zero
      LS_JACOBIANS_DB_UP = zero

!  Zero local sources

      lostrans_up   = zero  ; sources_up   = zero ; exactscat_up   = zero ; multiplier_up   = zero
      L_lostrans_up = zero  ; L_sources_up = zero ; L_exactscat_up = zero ; L_multiplier_up = zero

!  Bookkeeping

      ns = nstokes ; fmat = zero ; L_fmat = zero

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
!  --------------------------------------

!  Scalar only

      if ( nstokes .eq. 1 ) then
         do n = 1, nlayers
            if ( layermask_up(n) ) then
              sum = zero
               do L = 0, nmoments_input
                  sum = sum + GenSpher(L,1) * Greekmat(n,L,1,1)
               enddo
               exactscat_up(n,1,1) = sum * tms(n)
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_sum = zero
                        do L = 0, nmoments_input
                           L_sum = L_sum + GenSpher(L,1)  * L_Greekmat(n,L,1,1,q)
                        enddo
                        L_exactscat_up(n,1,1,q) = L_sum * tms(n) + sum * L_tms(n,q)
                     else
                        L_exactscat_up(n,1,1,q) = sum * L_tms(n,q)
                     endif
                  enddo
               endif
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
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_fmat(1:2) = zero
                        do L = 0, nmoments_input
                           L_fmat(1) = L_fmat(1) + GenSpher(L,1)  * L_Greekmat(n,L,1,1,q)
                           L_fmat(2) = L_fmat(2) + GenSpher(L,2)  * L_Greekmat(n,L,1,2,q)
                        enddo
                        L_exactscat_up(n,1,1,q) = + L_fmat(1)
                        L_exactscat_up(n,2,1,q) = - L_fmat(2) * Rotations(3)
                        L_exactscat_up(n,3,1,q) = + L_fmat(2) * Rotations(4)
                     endif
                     L_exactscat_up(n,1:ns,1,q) = L_exactscat_up(n,1:ns,1,q) *   tms(n) &
                                                  + exactscat_up(n,1:ns,1)   * L_tms(n,q)
                  enddo
               endif
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
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_fmat = zero
                        do L = 0, nmoments_input
                           L_fmat(1) = L_fmat(1) + GenSpher(L,1)  * L_Greekmat(n,L,1,1,q)
                           L_fmat(2) = L_fmat(2) + GenSpher(L,2)  * L_Greekmat(n,L,1,2,q)
                           sum23 = L_greekmat(n,L,2,2,q) + L_greekmat(n,L,3,3,q)
                           dif23 = L_greekmat(n,L,2,2,q) - L_greekmat(n,L,3,3,q)
                           L_fmat(3) = L_fmat(3) + Genspher(L,3) * sum23
                           L_fmat(4) = L_fmat(4) + Genspher(L,4) * dif23
                       enddo
                       L_fmat(3) = ( L_fmat(3) + L_fmat(4) ) * 0.5_fpk
                       L_fmat(4) = ( L_fmat(3) - L_fmat(4) )
                       if ( nstokes.eq.4) then
                          do L = 0, nmoments_input
                             L_fmat(5) = L_fmat(5) + Genspher(L,2) * L_greekmat(n,L,3,4,q)
                             L_fmat(6) = L_fmat(6) + Genspher(L,1) * L_greekmat(n,L,4,4,q)
                          enddo
                       endif
                       help3c1 = L_fmat(3) * Rotations(1)
                       help3s1 = L_fmat(3) * Rotations(2)
                       help4c1 = L_fmat(4) * Rotations(1)
                       help4s1 = L_fmat(4) * Rotations(2)
                       L_exactscat_up(n,1,1,q) = + L_fmat(1)
                       L_exactscat_up(n,2,1,q) = - L_fmat(2) * Rotations(3)
                       L_exactscat_up(n,3,1,q) = + L_fmat(2) * Rotations(4)
                       L_exactscat_up(n,2,2,q) = + help3c1 * Rotations(3) - help4s1 * Rotations(4)
                       L_exactscat_up(n,2,3,q) = - help3s1 * Rotations(3) - help4c1 * Rotations(4)
                       L_exactscat_up(n,3,2,q) = + help3c1 * Rotations(4) + help4s1 * Rotations(3)
                       L_exactscat_up(n,3,3,q) = - help3s1 * Rotations(4) + help4c1 * Rotations(3)
                       if ( nstokes .eq. 4 ) then
                          L_exactscat_up(n,2,4,q) = - L_fmat(5) * Rotations(4) 
                          L_exactscat_up(n,4,2,q) = - L_fmat(5) * Rotations(2) 
                          L_exactscat_up(n,3,4,q) = + L_fmat(5) * Rotations(3) 
                          L_exactscat_up(n,4,3,q) = - L_fmat(5) * Rotations(1) 
                          L_exactscat_up(n,4,4,q) = + L_fmat(6)
                       endif
                     endif
                     L_exactscat_up(n,1:ns,1:ns,q) = L_exactscat_up(n,1:ns,1:ns,q) *   tms(n) &
                                                     + exactscat_up(n,1:ns,1:ns)   * L_tms(n,q)
                  enddo
               endif
               exactscat_up(n,1:ns,1:ns) = tms(n)*exactscat_up(n,1:ns,1:ns)
            endif
         enddo
      endif

!  Attenuations
!  ============

!  Initialize, only to layer Ncrit if applicable

      Attenuations   = zero ; Attenuations_fine   = zero
      L_Attenuations = zero ; L_Attenuations_fine = zero
      Suntau         = zero ; L_suntau            = zero
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

!  Adjust nstart

      do n = 1, nlayers
         if ( layermask_up(n) .and. attenuations(n-1).ne.zero )  nstart = n
      enddo

!  Enhanced-spherical, fine-layer attenuations

      if ( do_enhanced_ps ) then
         do n = 1, nstart
            if ( layermask_up(n) ) then
               do j = 1, nfinedivs(n)
                  sumd = zero
                  do k = 1, ntraverse_fine(n,j)
                     sumd = sumd + extinction(k) * sunpaths_fine(n,k,j)
                  enddo
                  if (sumd .lt. cutoff ) Attenuations_fine(n,j) = exp( - sumd )
                  if ( do_profilewfs ) then
                     do k = 1, nlayers
                        if ( Qvary(k) .and. k.le.ntraverse_fine(n,j) ) then
                           do q = 1, Qnums(k)
                               L_sumd = L_extinction(k,q) * sunpaths_fine(n,k,j)
                               L_Attenuations_fine(n,j,k,q) = - Attenuations_fine(n,j) * L_sumd 
                           enddo
                        endif
                     enddo
                  else if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_sumd = zero
                        do k = 1, ntraverse_fine(n,j)
                           L_sumd = L_sumd + L_extinction(k,q) * sunpaths_fine(n,k,j)
                        enddo
                        L_Attenuations_fine(n,j,kc,q) = - Attenuations_fine(n,j) * L_sumd 
                     enddo
                  endif
               enddo
            endif
         enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Regular PS multipliers and LOSTRANS
!  -----------------------------------

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

!  Multipliers

            if ( layermask_up(n) .and. n.le.nstart  ) then
              if ( Mu1 .gt. zero ) then
                if ( attenuations(n-1).ne.zero ) then
                  factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(n)
                  factor2 = (suntau(n) - suntau(n-1))/lostau
                  factor3 = one + factor2
                  multiplier_up(n) = factor1 / factor3
                  if ( do_profilewfs ) then
                     do k = 1, nlayers
                        if ( Qvary(k).and.k.le.ntraverse(n) ) then
                           do q = 1, Qnums(k)
                              L_factor1 = L_Attenuations(n-1,k,q) - L_Attenuations(n,k,q)*lostrans_up(n)
                              L_factor2 = ( L_suntau(n,k,q) - L_suntau(n-1,k,q) ) / lostau
                              if ( k.eq.n ) then
                                 L_factor1 = L_factor1 - Attenuations(n) * L_lostrans_up(n,q)
                                 L_factor2 = L_factor2 - factor2 *  L_lostau(q) / lostau
                              endif
                              L_multiplier_up(n,k,q) = ( L_factor1 - multiplier_up(n) * L_factor2 ) / factor3
                           enddo
                        endif
                     enddo
                  else if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_factor1 = L_Attenuations(n-1,kc,q) - L_Attenuations(n,kc,q)*lostrans_up(n)
                        L_factor2 = ( L_suntau(n,kc,q) - L_suntau(n-1,kc,q) ) / lostau
                        L_factor1 = L_factor1 - Attenuations(n) * L_lostrans_up(n,q)
                        L_factor2 = L_factor2 - factor2 *  L_lostau(q) / lostau
                        L_multiplier_up(n,kc,q) = ( L_factor1 - multiplier_up(n) * L_factor2 ) / factor3
                     enddo
                  endif
                endif
              else
                if ( attenuations(n-1).ne.zero ) then
                  multiplier_up(n) = Attenuations(n-1)
                  if ( do_profilewfs ) then
                     do k = 1, nlayers
                        if ( Qvary(k).and.k.le.ntraverse(n) ) then
                           do q = 1, Qnums(k)
                              L_multiplier_up(n,k,q) = L_Attenuations(n-1,k,q)
                           enddo
                        endif
                     enddo
                  else if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_multiplier_up(n,kc,q) = L_Attenuations(n-1,kc,q)
                     enddo
                  endif
                endif
              endif
            endif

!  End layers and regular-PS formulation

         enddo
      endif

!  Enhanced PS multipliers and LOSTRANS: special case (nadir viewing)
!  ------------------------------------------------------------------

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

!  Multipliers

            if ( layermask_up(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  argum(j) = xfine(n,j)
                  tran(j)  = exp ( - argum(j) * kn )
                  func(j)  = attenuations_fine(n,j) * tran(j)
                  sum = sum + func(j) * wfine(n,j)
               enddo
               multiplier_up(n) = sum * kn
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = L_attenuations_fine(n,j,k,q) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_multiplier_up(n,k,q)  = L_sum * kn + sum * L_extinction(n,q)
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_func = L_attenuations_fine(n,j,k,q)  * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_multiplier_up(n,k,q)  = L_sum * kn
                           endif
                        enddo
                     endif
                  enddo
               else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = L_attenuations_fine(n,j,kc,q) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j)
                     enddo
                     L_multiplier_up(n,kc,q)  = L_sum * kn + L_extinction(n,q) * sum
                  enddo
               endif
            endif

!  End layer loop and Nadir Enhanced PS case

         enddo
      endif

!  Enhanced PS multipliers and LOSTRANS: General case
!  --------------------------------------------------

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
                  enddo
               endif
            endif

!  Multipliers

            if ( layermask_up(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  argum(j) = Raycon * ( cot_2 - cotfine(n,j) )
                  tran(j)  = exp ( - kn * argum(j) )
                  func(j)  = attenuations_fine(n,j) * csqfine(n,j) * tran(j)
                  sum      = sum + func(j) * wfine(n,j)
               enddo
               Multiplier_up(n) = sum * ke
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = L_attenuations_fine(n,j,k,q) * csqfine(n,j) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_multiplier_up(n,k,q)  = L_sum * ke + L_extinction(N,q) * Raycon * sum
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_func = L_attenuations_fine(n,j,k,q) * csqfine(n,j) * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_multiplier_up(n,k,q)  = L_sum * ke
                           endif
                        enddo
                     endif
                  enddo
               else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = L_attenuations_fine(n,j,kc,q) * csqfine(n,j) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j)
                     enddo
                     L_multiplier_up(n,kc,q)  = L_sum * ke + L_extinction(N,q) * Raycon * sum
                  enddo
               endif
            endif

!  End layer loop and general Enhanced PS case

         enddo
      endif

!  Layer sources
!  -------------

      do n = nlayers, 1, -1
         if ( layermask_up(n) .and. n.le.nstart  ) then
            do o1 = 1, nstokes
               shelp(o1) = dot_product(exactscat_up(n,o1,1:ns),fluxvec(1:ns))
               sources_up(n,o1) = shelp(o1) * multiplier_up(n)
            enddo
            if ( do_profilewfs ) then
               do k = 1, nlayers
                  if ( Qvary(k) ) then
                     do q = 1, Qnums(k)
                        do o1 = 1, nstokes
                           L_sources_up(n,o1,k,q) = shelp(o1) * L_multiplier_up(n,k,q)
                           if ( k.eq.n ) then
                              L_Shelp = dot_product(L_exactscat_up(n,o1,1:ns,q),fluxvec(1:ns))
                              L_sources_up(n,o1,k,q) =  L_sources_up(n,o1,k,q) + L_Shelp * multiplier_up(n)
                           endif
                        enddo
                     enddo
                  endif
               enddo
            else if ( do_columnwfs ) then
               do q = 1, n_columnwfs
                  do o1 = 1, nstokes
                     L_sources_up(n,o1,kc,q) = shelp(o1) * L_multiplier_up(n,kc,q)
                     L_Shelp = dot_product(L_exactscat_up(n,o1,1:ns,q),fluxvec(1:ns))
                     L_sources_up(n,o1,kc,q) =  L_sources_up(n,o1,kc,q) + L_Shelp * multiplier_up(n)
                  enddo
               enddo
            endif
         endif
      enddo

!  Source function integration
!  ===========================

!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!  Check for updating the recursion

!  INTENSITY Main loop over all output optical depths
!          Cumulative source term will be saved

      NC = 0 
      CUMSOURCE_UP(NC,:) = zero
      CUMSOURCE_DB(NC,:) = zero
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Surface term

      RHELP = zero; M4 = 4.0_fpk * MU0 ; M4A = M4 * attenuations(nlayers)
      if ( DO_LAMBERTIAN ) then
         RHELP(1) = M4 * REFLEC(1,1) * Fluxvec(1)
         CUMSOURCE_DB(NC,1) = RHELP(1) * attenuations(nlayers)
      else
         do o1 = 1, nstokes
            RHELP(O1) = M4 * dot_product(REFLEC(O1,1:ns),fluxvec(1:ns))
            CUMSOURCE_DB(NC,o1) = RHELP(O1) * attenuations(nlayers)
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
               CUMSOURCE_DB(NC,O1) = LOSTRANS_UP(N) * CUMSOURCE_DB(NC-1,O1)
               CUMSOURCE_UP(NC,O1) = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1,O1) + SOURCES_UP(N,O1)
            enddo
         ENDDO
         do o1 = 1, nstokes
            STOKES_UP   (UTA,O1) = FLUX * CUMSOURCE_UP(NC,O1)
            STOKES_DB_UP(UTA,O1) = FLUX * CUMSOURCE_DB(NC,O1)
         enddo
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
      ENDDO

!  Surface WFs

      if ( do_surfacewfs ) then
         LS_CUMSOURCE = zero
         do q = 1, n_surfacewfs
            if ( DO_LAMBERTIAN ) then
               LS_cumsource(1,q) = M4A * LS_REFLEC(1,1,q)
            else
               do o1 = 1, nstokes
                  LS_cumsource(o1,q) = M4A * dot_product(LS_REFLEC(O1,1:ns,q),fluxvec(1:ns))
               enddo
            endif
         enddo
         NSTART = NLAYERS ; NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               do q = 1, n_surfacewfs
                  LS_cumsource(1:ns,q) = LOSTRANS_UP(N) * LS_CUMSOURCE(1:ns,q)
               enddo
            ENDDO
            do q = 1, n_surfacewfs
               LS_JACOBIANS_DB_UP(UTA,1:ns,Q) = FLUX * LS_CUMSOURCE(1:ns,Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
         ENDDO
      endif

!  Profile Wfs (atmospheric term)

      if ( do_profilewfs ) then
         do k = 1, nlayers
            if ( Qvary(k) ) then
               L_CUMSOURCE = zero
               NSTART = NLAYERS ; NUT_PREV = NSTART + 1
               DO UTA = N_USER_LEVELS, 1, -1
                  NUT    = USER_LEVELS(UTA) + 1
                  DO N = NSTART, NUT, -1
                     NC = NLAYERS + 1 - N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(1:ns,q) = L_SOURCES_UP(N,1:ns,K,Q)    + &
                                L_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1,1:ns) + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(1:ns,Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(1:ns,q) = L_SOURCES_UP(N,1:ns,K,Q)  + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(1:ns,Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_JACOBIANS_UP(UTA,1:ns,K,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Profile Wfs (direct beam term)

      if ( do_profilewfs ) then
         do k = 1, nlayers
            if ( Qvary(k) ) then
               do q = 1, Qnums(k)
                  L_CUMSOURCE(1:ns,q) = RHELP(1:ns) * L_attenuations(nlayers,k,q)
               enddo
               NSTART = NLAYERS ; NUT_PREV = NSTART + 1
               DO UTA = N_USER_LEVELS, 1, -1
                  NUT    = USER_LEVELS(UTA) + 1
                  DO N = NSTART, NUT, -1
                     NC = NLAYERS + 1 - N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(1:ns,q) =  L_LOSTRANS_UP(N,Q) * CUMSOURCE_DB(NC-1,1:ns) + &
                                                    LOSTRANS_UP(N)   * L_CUMSOURCE(1:ns,Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(1:ns,q) = LOSTRANS_UP(N) * L_CUMSOURCE(1:ns,Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_JACOBIANS_DB_UP(UTA,1:ns,K,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Column Wfs (Atmospheric term)

      if ( do_columnwfs ) then
         L_CUMSOURCE = zero
         NSTART = NLAYERS  ;  NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) = L_SOURCES_UP(N,1:ns,KC,Q) + &
                                L_LOSTRANS_UP(N,Q) * CUMSOURCE_UP(NC-1,1:ns) + &
                                  LOSTRANS_UP(N)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_JACOBIANS_UP(UTA,1:ns,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ;  NUT_PREV = NUT
         ENDDO
      endif

!  Column Wfs (Surface term)

      if ( do_columnwfs ) then
         do q = 1, n_columnwfs
            L_CUMSOURCE(1:ns,q) = RHELP(1:ns) * L_attenuations(nlayers,kc,q)
         enddo
         NSTART = NLAYERS  ;  NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) =  L_LOSTRANS_UP(N,Q) * CUMSOURCE_DB(NC-1,1:ns) + &
                                      LOSTRANS_UP(N)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_JACOBIANS_DB_UP(UTA,1:ns,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1 ; NUT_PREV = NUT
         ENDDO
      endif

!  Finish

      return
end subroutine SSV_Integral_Plus_1_UP

!

subroutine SSV_Integral_Plus_1_DN &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels, do_sunlight,  & ! Inputs (dimensioning/flag)
     max_atmoswfs, do_columnwfs, do_profilewfs,                                 & ! Inputs (dimensioning/flags)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,                 & ! Inputs (Flags - General  )
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,            & ! Inputs (Layers/Levels control)
     nstokes, n_columnwfs, Lvaryflags, Lvarynums, Lvarymoms,                    & ! Inputs (Jacobian control)
     extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,             & ! Inputs (Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat,                  & ! Inputs (Optical - Linearized)
     Mu1, GenSpher, Rotations, NCrit, RadCrit, CotCrit,                         & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                       & ! Inputs (Geometry)
     sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,                        & ! Inputs (Geometry)
     Stokes_dn, LC_Jacobians_dn, LP_Jacobians_dn  )                               ! Output

!  Stand alone routine for SS field (Radiances and Jacobians) with Solar sources alone
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

      INTEGER, Intent(in) :: max_atmoswfs

!  General flags

      LOGICAL, Intent(in) ::  DO_SUNLIGHT

      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DONADIR

!  Jacobian Flags

      LOGICAL, Intent(in) ::  do_columnwfs
      LOGICAL, Intent(in) ::  do_profilewfs

!  Numbers

      INTEGER, Intent(in) ::  NSTOKES
      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  NMOMENTS_INPUT

      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

      INTEGER, Intent(in) ::  n_columnwfs
      LOGICAL, Intent(in) ::  Lvaryflags(maxlayers)
      INTEGER, Intent(in) ::  Lvarynums (maxlayers)
      LOGICAL, Intent(in) ::  Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Atmosphere

      REAL(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
      REAL(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
      REAL(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
      REAL(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
      REAL(fpk), Intent(in) :: GREEKMAT    ( MAXLAYERS, 0:MAXMOMENTS_INPUT, 4, 4 )

!  Solar Flux

      REAL(fpk), Intent(in) :: FLUX, FLUXVEC(4)

!  Linearized optical inputs

      REAL(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4, max_atmoswfs )

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!       Mu1 = cos(alpha_boa), required for the Regular PS only

      integer  , Intent(in)  :: NCrit
      real(fpk), Intent(in)  :: Raycon, radii(0:maxlayers),cota(0:maxlayers)
      real(fpk), Intent(in)  :: Mu1
      real(fpk), Intent(in)  :: Radcrit, CotCrit

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

      real(fpk), Intent(Out)  :: Stokes_dn        ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: LC_Jacobians_dn  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dn  ( max_user_levels, 4, maxlayers, max_atmoswfs )

!  LOCAL
!  -----

!  Attenuations

      real(fpk)  :: suntau              (0:maxlayers)
      real(fpk)  :: attenuations        (0:maxlayers)
      real(fpk)  :: attenuations_fine   (maxlayers,maxfinelayers)

      real(fpk)  :: L_suntau            (0:maxlayers,0:maxlayers,max_atmoswfs)
      real(fpk)  :: L_attenuations      (0:maxlayers,0:maxlayers,max_atmoswfs)
      real(fpk)  :: L_attenuations_fine (maxlayers,maxfinelayers,0:maxlayers,max_atmoswfs)

!  Scattering

      real(fpk)  :: tms            (maxlayers)
      real(fpk)  :: exactscat_dn   (maxlayers,4,4)

      real(fpk)  :: L_tms          (maxlayers,max_atmoswfs)
      real(fpk)  :: L_exactscat_dn (maxlayers,4,4, max_atmoswfs)

!  Source function integration results

      real(fpk)  :: sources_dn       ( maxlayers, 4 )
      real(fpk)  :: lostrans_dn      ( maxlayers )
      real(fpk)  :: multiplier_dn    ( maxlayers )

      real(fpk)  :: L_sources_dn      ( maxlayers,4,0:maxlayers,max_atmoswfs )
      real(fpk)  :: L_lostrans_dn     ( maxlayers, max_atmoswfs )
      real(fpk)  :: L_multiplier_dn   ( maxlayers,0:maxlayers,max_atmoswfs )

!  Local cumulative source terms

      real(fpk)  :: cumsource_dn      ( 0:maxlayers, 4 )
      real(fpk)  :: L_cumsource       ( 4, max_atmoswfs )

!  Help

      integer    :: n, ns, k, kc, j, q, L, o1, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers)
      logical    :: do_atmoswfs, layermask_dn(maxlayers), Qvary(maxlayers)

      real(fpk)  :: argum(maxfinelayers), tran(maxfinelayers), func(maxfinelayers)

      real(fpk)  :: cons, help, sum, tran_1, kn, ke, factor1, factor2, factor3, shelp(4)
      real(fpk)  :: cot_1, cot_2, L_help, L_sum, L_tran, L_func, L_factor1, L_factor2, sumd, L_sumd
      real(fpk)  :: lostau, L_lostau(max_atmoswfs), fmat(6), L_fmat(6), sum23, dif23, L_Shelp
      real(fpk)  :: help3c1, help3s1, help4c1, help4s1, cot_c, rdiff

      real(fpk), parameter  :: cutoff = 88.0_fpk
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

      STOKES_DN       = zero
      LP_JACOBIANS_DN = zero
      LC_JACOBIANS_DN = zero

!  Zero local sources

      lostrans_dn   = zero  ; sources_dn   = zero ; exactscat_dn   = zero ; multiplier_dn   = zero
      L_lostrans_dn = zero  ; L_sources_dn = zero ; L_exactscat_dn = zero ; L_multiplier_dn = zero

!  Bookkeeping

      ns = nstokes ; fmat = zero ; L_fmat = zero

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
!  --------------------------------------

!  Scalar only

      if ( nstokes .eq. 1 ) then
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
              sum = zero
               do L = 0, nmoments_input
                  sum = sum + GenSpher(L,1) * Greekmat(n,L,1,1)
               enddo
               exactscat_dn(n,1,1) = sum * tms(n)
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_sum = zero
                        do L = 0, nmoments_input
                           L_sum = L_sum + GenSpher(L,1)  * L_Greekmat(n,L,1,1,q)
                        enddo
                        L_exactscat_dn(n,1,1,q) = L_sum * tms(n) + sum * L_tms(n,q)
                     else
                        L_exactscat_dn(n,1,1,q) = sum * L_tms(n,q)
                     endif
                  enddo
               endif               
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
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_fmat(1:2) = zero
                        do L = 0, nmoments_input
                           L_fmat(1) = L_fmat(1) + GenSpher(L,1)  * L_Greekmat(n,L,1,1,q)
                           L_fmat(2) = L_fmat(2) + GenSpher(L,2)  * L_Greekmat(n,L,1,2,q)
                        enddo
                        L_exactscat_dn(n,1,1,q) = + L_fmat(1)
                        L_exactscat_dn(n,2,1,q) = - L_fmat(2) * Rotations(3)
                        L_exactscat_dn(n,3,1,q) = + L_fmat(2) * Rotations(4)
                     endif
                     L_exactscat_dn(n,1:ns,1,q) = L_exactscat_dn(n,1:ns,1,q) *   tms(n) &
                                                  + exactscat_dn(n,1:ns,1)   * L_tms(n,q)
                  enddo
               endif               
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
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     if ( Lvarymoms(n,q) ) then
                        L_fmat = zero
                        do L = 0, nmoments_input
                           L_fmat(1) = L_fmat(1) + GenSpher(L,1)  * L_Greekmat(n,L,1,1,q)
                           L_fmat(2) = L_fmat(2) + GenSpher(L,2)  * L_Greekmat(n,L,1,2,q)
                           sum23 = L_greekmat(n,L,2,2,q) + L_greekmat(n,L,3,3,q)
                           dif23 = L_greekmat(n,L,2,2,q) - L_greekmat(n,L,3,3,q)
                           L_fmat(3) = L_fmat(3) + Genspher(L,3) * sum23
                           L_fmat(4) = L_fmat(4) + Genspher(L,4) * dif23
                       enddo
                       L_fmat(3) = ( L_fmat(3) + L_fmat(4) ) * 0.5_fpk
                       L_fmat(4) = ( L_fmat(3) - L_fmat(4) )
                       if ( nstokes.eq.4) then
                          do L = 0, nmoments_input
                             L_fmat(5) = L_fmat(5) + Genspher(L,2) * L_greekmat(n,L,3,4,q)
                             L_fmat(6) = L_fmat(6) + Genspher(L,1) * L_greekmat(n,L,4,4,q)
                          enddo
                       endif
                       help3c1 = L_fmat(3) * Rotations(1)
                       help3s1 = L_fmat(3) * Rotations(2)
                       help4c1 = L_fmat(4) * Rotations(1)
                       help4s1 = L_fmat(4) * Rotations(2)
                       L_exactscat_dn(n,1,1,q) = + L_fmat(1)
                       L_exactscat_dn(n,2,1,q) = - L_fmat(2) * Rotations(3)
                       L_exactscat_dn(n,3,1,q) = + L_fmat(2) * Rotations(4)
                       L_exactscat_dn(n,2,2,q) = + help3c1 * Rotations(3) - help4s1 * Rotations(4)
                       L_exactscat_dn(n,2,3,q) = - help3s1 * Rotations(3) - help4c1 * Rotations(4)
                       L_exactscat_dn(n,3,2,q) = + help3c1 * Rotations(4) + help4s1 * Rotations(3)
                       L_exactscat_dn(n,3,3,q) = - help3s1 * Rotations(4) + help4c1 * Rotations(3)
                       if ( nstokes .eq. 4 ) then
                          L_exactscat_dn(n,2,4,q) = - L_fmat(5) * Rotations(4) 
                          L_exactscat_dn(n,4,2,q) = - L_fmat(5) * Rotations(2) 
                          L_exactscat_dn(n,3,4,q) = + L_fmat(5) * Rotations(3) 
                          L_exactscat_dn(n,4,3,q) = - L_fmat(5) * Rotations(1) 
                          L_exactscat_dn(n,4,4,q) = + L_fmat(6)
                       endif
                     endif
                     L_exactscat_dn(n,1:ns,1:ns,q) = L_exactscat_dn(n,1:ns,1:ns,q) *   tms(n) &
                                                     + exactscat_dn(n,1:ns,1:ns)   * L_tms(n,q)
                  enddo
               endif               
               exactscat_dn(n,1:ns,1:ns) = tms(n)*exactscat_dn(n,1:ns,1:ns)
            endif
         enddo
      endif

!  Attenuations
!  ============

!  Initialize, only to layer Ncrit if applicable

      Attenuations   = zero ; Attenuations_fine   = zero
      L_Attenuations = zero ; L_Attenuations_fine = zero
      Suntau         = zero ; L_suntau            = zero
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
                  if ( do_profilewfs ) then
                     do k = 1, nlayers
                        if ( Qvary(k) .and. k.le.ntraverse_fine(n,j) ) then
                           do q = 1, Qnums(k)
                               L_sumd = L_extinction(k,q) * sunpaths_fine(n,k,j)
                               L_Attenuations_fine(n,j,k,q) = - Attenuations_fine(n,j) * L_sumd 
                           enddo
                        endif
                     enddo
                  else if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_sumd = zero
                        do k = 1, ntraverse_fine(n,j)
                           L_sumd = L_sumd + L_extinction(k,q) * sunpaths_fine(n,k,j)
                        enddo
                        L_Attenuations_fine(n,j,kc,q) = - Attenuations_fine(n,j) * L_sumd 
                     enddo
                  endif
               enddo
            endif
         enddo
      endif

!  Regular PS multipliers and LOSTRANS
!  -----------------------------------

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

!  Multipliers

            if ( layermask_dn(n) .and. n.le.nstart  ) then
              if ( Mu1 .gt. zero ) then
                if ( attenuations(n-1).ne.zero ) then
                  factor1 = Attenuations(n-1)*lostrans_dn(n) - Attenuations(n)
                  factor2 = (suntau(n) - suntau(n-1))/lostau
                  factor3 = factor2 - one
                  multiplier_dn(n) = factor1 / factor3
                  if ( do_profilewfs ) then
                     do k = 1, nlayers
                        if ( Qvary(k).and.k.le.ntraverse(n) ) then
                           do q = 1, Qnums(k)
                              L_factor1 = L_Attenuations(n-1,k,q)*lostrans_dn(n) - L_Attenuations(n,k,q)
                              L_factor2 = ( L_suntau(n,k,q) - L_suntau(n-1,k,q) ) / lostau
                              if ( k.eq.n ) then
                                 L_factor1 = L_factor1 + Attenuations(n-1) * L_lostrans_dn(n,q)
                                 L_factor2 = L_factor2 - factor2 *  L_lostau(q) / lostau
                              endif
                              L_multiplier_dn(n,k,q) = ( L_factor1 - multiplier_dn(n) * L_factor2 ) / factor3
                           enddo
                        endif
                     enddo
                  else if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_factor1 = L_Attenuations(n-1,kc,q)*lostrans_dn(n) - L_Attenuations(n,kc,q)
                        L_factor2 = ( L_suntau(n,kc,q) - L_suntau(n-1,kc,q) ) / lostau
                        L_factor1 = L_factor1 + Attenuations(n-1) * L_lostrans_dn(n,q)
                        L_factor2 = L_factor2 - factor2 *  L_lostau(q) / lostau
                        L_multiplier_dn(n,kc,q) = ( L_factor1 - multiplier_dn(n) * L_factor2 ) / factor3
                     enddo
                  endif
                endif
              else
                if ( attenuations(n-1).ne.zero ) then
                  multiplier_dn(n) = Attenuations(n)
                  if ( do_profilewfs ) then
                     do k = 1, nlayers
                        if ( Qvary(k).and.k.le.ntraverse(n) ) then
                           do q = 1, Qnums(k)
                              L_multiplier_dn(n,k,q) = L_Attenuations(n-1,k,q)
                           enddo
                        endif
                     enddo
                  else if ( do_columnwfs ) then
                     do q = 1, n_columnwfs
                        L_multiplier_dn(n,kc,q) = L_Attenuations(n-1,kc,q)
                     enddo
                  endif
                endif
              endif
            endif

!  End layers and regular-PS formulation

         enddo
      endif

!  Enhanced PS multipliers and LOSTRANS: special case (nadir viewing)
!  ------------------------------------------------------------------

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

!  Multipliers

            if ( layermask_dn(n) .and. n.le.nstart  ) then
               rdiff = radii(n-1) - radii(n)
               if ( n.eq.NCrit) rdiff = radii(n-1) - RadCrit
               sum = zero
               do j = 1, nfinedivs(n)
                  argum(j) = rdiff - xfine(n,j)
                  tran(j)  = exp ( - argum(j) * kn )
                  func(j)  = attenuations_fine(n,j) * tran(j)
                  sum = sum + func(j) * wfine(n,j)
               enddo
               multiplier_dn(n) = sum * kn
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = L_attenuations_fine(n,j,k,q) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_multiplier_dn(n,k,q)  = L_sum * kn + L_extinction(N,q) * sum
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_func = L_attenuations_fine(n,j,k,q)  * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_multiplier_dn(n,k,q)  = L_sum * kn
                           endif
                        enddo
                     endif
                  enddo
               else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = L_attenuations_fine(n,j,kc,q) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j)
                     enddo
                     L_multiplier_dn(n,kc,q)  = L_sum * kn + L_extinction(N,q) * sum
                  enddo
               endif
            endif

!  End layer loop and Nadir Enhanced PS case

         enddo
      endif

!  Enhanced PS multipliers and LOSTRANS: General case
!  --------------------------------------------------

      if ( do_enhanced_ps .and. .not. doNadir ) then
         do n = nlayers, 1, -1

!  LOS transmittance

            cot_2 = cota(n-1) ; cot_1 = cota(n)
            cot_c = cot_1  ; if ( n.eq.NCrit ) cot_c = CotCrit
            kn = extinction(n) ;  ke = raycon * kn ; cons = raycon * ( cot_2 - cot_1 )
            tran_1 = kn * cons
            if ( tran_1 .lt. cutoff ) lostrans_dn(n) = exp ( - tran_1 )
            if ( do_atmoswfs ) then
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_lostau(q)        = L_extinction(n,q) * cons
                     L_lostrans_dn(n,q) = - L_lostau(q) * lostrans_dn(n)
                  enddo
               endif
            endif

!  multiplier

            if ( layermask_dn(n) .and. n.le.nstart  ) then
               sum = zero
               do j = 1, nfinedivs(n)
                  argum(j) = Raycon * ( cotfine(n,j) - cot_c )
                  tran(j)  = exp ( - kn * argum(j) )
                  func(j)  = attenuations_fine(n,j) * csqfine(n,j) * tran(j)
                  sum      = sum + func(j) * wfine(n,j)
               enddo
               multiplier_dn(n) = sum * ke 
               if ( do_profilewfs ) then
                  do k = 1, nlayers
                     if ( Qvary(k).and.k.le.ntraverse(n) ) then
                        do q = 1, Qnums(k)
                           if ( k.eq.n ) then
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_tran = - argum(j) * L_extinction(n,q)
                                 L_func = L_attenuations_fine(n,j,k,q) * csqfine(n,j) * tran(j) + L_tran * func(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_multiplier_dn(n,k,q)  = L_sum * ke + L_extinction(N,q) * Raycon * sum
                           else
                              L_sum = zero
                              do j = 1, nfinedivs(n)
                                 L_func = L_attenuations_fine(n,j,k,q) * csqfine(n,j) * tran(j)
                                 L_sum  = L_sum + L_func * wfine(n,j)
                              enddo
                              L_multiplier_dn(n,k,q)  = L_sum * ke
                           endif
                        enddo
                     endif
                  enddo
               else if ( do_columnwfs ) then
                  do q = 1, n_columnwfs
                     L_sum = zero
                     do j = 1, nfinedivs(n)
                        L_tran = - argum(j) * L_extinction(n,q)
                        L_func = L_attenuations_fine(n,j,kc,q) * csqfine(n,j) * tran(j) + L_tran * func(j)
                        L_sum  = L_sum + L_func * wfine(n,j)
                     enddo
                     L_multiplier_dn(n,kc,q)  = L_sum * ke + L_extinction(N,q) * Raycon * sum
                  enddo
               endif
            endif

!  End layer loop and general Enhanced PS case

         enddo
      endif

!  Layer sources
!  -------------

      do n = nlayers, 1, -1
         if ( layermask_dn(n) .and. n.le.nstart  ) then
            do o1 = 1, nstokes
               shelp(o1) = dot_product(exactscat_dn(n,o1,1:ns),fluxvec(1:ns))
               sources_dn(n,o1) = shelp(o1) * multiplier_dn(n)
            enddo
            if ( do_profilewfs ) then
               do k = 1, nlayers
                  if ( Qvary(k) ) then
                     do q = 1, Qnums(k)
                        do o1 = 1, nstokes
                           L_sources_dn(n,o1,k,q) = shelp(o1) * L_multiplier_dn(n,k,q)
                           if ( k.eq.n ) then
                              L_Shelp = dot_product(L_exactscat_dn(n,o1,1:ns,q),fluxvec(1:ns))
                              L_sources_dn(n,o1,k,q) =  L_sources_dn(n,o1,k,q) + L_Shelp * multiplier_dn(n)
                           endif
                        enddo
                     enddo
                  endif
               enddo
            else if ( do_columnwfs ) then
               do q = 1, n_columnwfs
                  do o1 = 1, nstokes
                     L_sources_dn(n,o1,kc,q) = shelp(o1) * L_multiplier_dn(n,kc,q)
                     L_Shelp = dot_product(L_exactscat_dn(n,o1,1:ns,q),fluxvec(1:ns))
                     L_sources_dn(n,o1,kc,q) =  L_sources_dn(n,o1,kc,q) + L_Shelp * multiplier_dn(n)
                  enddo
               enddo
            endif
         endif
      enddo


!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC,:) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

      NC = 0
      CUMSOURCE_DN(NC,:) = zero
      NSTART = 1 ; NUT_PREV = NSTART - 1
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
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1 ; NUT_PREV = NUT
      ENDDO

!  Profile Wfs

      if ( do_profilewfs ) then
         do k = 1, nlayers
            if ( Qvary(k) ) then
               L_CUMSOURCE = zero
               NSTART = 1 ; NUT_PREV = NSTART - 1
               DO UTA = 1, N_USER_LEVELs
                  NUT    = USER_LEVELS(UTA)
                  DO N = NSTART, NUT
                     NC = N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(1:ns,q) = L_SOURCES_DN(N,1:ns,K,Q)    + &
                                L_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1,1:ns) + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(1:ns,Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(1:ns,q) = L_SOURCES_DN(N,1:ns,K,Q)  + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(1:ns,Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_JACOBIANS_DN(UTA,1:ns,K,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1  ;  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Column Wfs

      if ( do_columnwfs ) then
         L_CUMSOURCE = zero
         NSTART = 1 ; NUT_PREV = NSTART - 1
         DO UTA = 1, N_USER_LEVELS
            NUT    = USER_LEVELS(UTA)
            DO N = NSTART, NUT
               NC = N
               do q = 1, n_columnwfs
                  L_cumsource(1:ns,q) = L_SOURCES_DN(N,1:ns,KC,Q)         + &
                                L_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1,1:ns) + &
                                  LOSTRANS_DN(N)   * L_CUMSOURCE(1:ns,Q)
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_JACOBIANS_DN(UTA,1:ns,Q) = FLUX * L_CUMSOURCE(1:ns,Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1 ;  NUT_PREV = NUT
         ENDDO
      endif

!  Finish

      return
end subroutine SSV_Integral_Plus_1_DN


!

subroutine SSV_Integral_Plus_1_UPDN &
   ( maxlayers, maxfinelayers, maxmoments_input,                                  & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs, do_sunlight,                  & ! Inputs (dimensioning)
     do_deltam_scaling, do_upwelling, do_dnwelling, do_regular_ps, do_lambertian, & ! Inputs (Flags)
     do_enhanced_ps, doNadir, do_columnwfs, do_profilewfs, do_surfacewfs,         & ! Inputs (Flags)
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,              & ! Inputs (Layers/Levels control)
     nstokes, n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,        & ! Inputs (Jacobian control)
     NCrit, Mu0, Mu1, Genspher_up, Genspher_dn, Rotations_up, Rotations_dn,       & ! Inputs (Misc. Geometry)
     Raycon, radii, cota, xfine, wfine, csqfine, cotfine,                         & ! Inputs (Whole layer Geometry)
     RadCrit, CotCrit, sunpaths_up, ntraverse_up, sunpaths_dn, ntraverse_dn,      & ! Inputs (Whole layer Geometry)
     sunpaths_fine_up, ntraverse_fine_up, sunpaths_fine_dn, ntraverse_fine_dn,    & ! Inputs (Fine  layer Geometry)
     Flux, FLuxvec, reflec, extinction, deltaus, omega, truncfac, greekmat,       & ! Inputs (Optical - Regular)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat,         & ! Inputs (Optical - Linearized)
     Stokes_up, Stokes_db_up, LC_Jacobians_up, LP_Jacobians_up,                   & ! Output
     LC_Jacobians_db_up, LP_Jacobians_db_up, LS_Jacobians_db_up,                  & ! Output
     Stokes_dn, LC_Jacobians_dn, LP_Jacobians_dn )                                  ! Output

   Implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

      INTEGER, Intent(in) :: maxlayers
      INTEGER, Intent(in) :: maxfinelayers
      INTEGER, Intent(in) :: maxmoments_input
      INTEGER, Intent(in) :: max_user_levels

      INTEGER, Intent(in) :: max_atmoswfs
      INTEGER, Intent(in) :: max_surfacewfs

!  General flags

      LOGICAL, Intent(in) ::  DO_SUNLIGHT
      LOGICAL, Intent(in) ::  DO_UPWELLING
      LOGICAL, Intent(in) ::  DO_DNWELLING
      LOGICAL, Intent(in) ::  DO_REGULAR_PS
      LOGICAL, Intent(in) ::  DO_ENHANCED_PS
      LOGICAL, Intent(in) ::  DO_DELTAM_SCALING
      LOGICAL, Intent(in) ::  DO_LAMBERTIAN
      LOGICAL, Intent(in) ::  DONADIR

!  Jacobian Flags

      LOGICAL, Intent(in) ::  do_surfacewfs
      LOGICAL, Intent(in) ::  do_columnwfs
      LOGICAL, Intent(in) ::  do_profilewfs

!  Numbers

      INTEGER, Intent(in) ::  NSTOKES
      INTEGER, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS)
      INTEGER, Intent(in) ::  NMOMENTS_INPUT

      INTEGER, Intent(in) ::  N_USER_LEVELS
      INTEGER, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Jacobian control

      INTEGER, Intent(in) ::  n_columnwfs
      INTEGER, Intent(in) ::  n_surfacewfs
      LOGICAL, Intent(in) ::  Lvaryflags(maxlayers)
      INTEGER, Intent(in) ::  Lvarynums (maxlayers)
      LOGICAL, Intent(in) ::  Lvarymoms (maxlayers,max_atmoswfs)

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

!  Linearized optical inputs

      REAL(fpk), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
      REAL(fpk), Intent(in) :: L_GREEKMAT    ( MAXLAYERS,0:MAXMOMENTS_INPUT, 4, 4, max_atmoswfs )
      REAL(fpk), Intent(in) :: LS_REFLEC     ( 4, 4, max_surfacewfs )

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

      real(fpk), Intent(Out)  :: Stokes_up     ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: Stokes_db_up  ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: LC_Jacobians_up  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_up  ( max_user_levels, 4, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LC_Jacobians_db_up  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_db_up  ( max_user_levels, 4, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LS_Jacobians_db_up  ( max_user_levels, 4, max_surfacewfs )

      real(fpk), Intent(Out)  :: Stokes_dn        ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: LC_Jacobians_dn  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dn  ( max_user_levels, 4, maxlayers, max_atmoswfs )

!  Code

   if ( do_upwelling  ) then
      call SSV_Integral_Plus_1_UP &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels, do_sunlight,      & ! Inputs (dimensioning/flag)
     max_atmoswfs, max_surfacewfs, do_columnwfs, do_profilewfs, do_surfacewfs,      & ! Inputs (dimensioning/flags)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, do_lambertian, doNadir,      & ! Inputs (Flags - General)
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,                & ! Inputs (Layers/Levels)
     nstokes, n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,          & ! Inputs (Jacobian control)
     reflec, extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,         & ! Inputs (Optical)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat,           & ! Inputs (Optical Linearized)
     Mu0, Mu1, GenSpher_up, Rotations_up, NCrit, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
     Raycon, cota, sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up,  & ! Inputs (Geometry)
     Stokes_up, Stokes_db_up, LC_Jacobians_up, LP_Jacobians_up,                     & ! Output
     LC_Jacobians_db_up, LP_Jacobians_db_up, LS_Jacobians_db_up )                     ! Output
   endif

   if ( do_dnwelling  ) then
      call SSV_Integral_Plus_1_DN &
   ( maxlayers, maxfinelayers, maxmoments_input, max_user_levels, do_sunlight,  & ! Inputs (dimensioning/flag)
     max_atmoswfs, do_columnwfs, do_profilewfs,                                 & ! Inputs (dimensioning/flags)
     do_deltam_scaling, do_regular_ps, do_enhanced_ps, doNadir,                 & ! Inputs (Flags - General  )
     nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,            & ! Inputs (Layers/Levels control)
     nstokes, n_columnwfs, Lvaryflags, Lvarynums, Lvarymoms,                    & ! Inputs (Jacobian control)
     extinction, deltaus, omega, truncfac, greekmat, flux, fluxvec,             & ! Inputs (Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_greekmat,                  & ! Inputs (Optical - Linearized)
     Mu1, GenSpher_dn, Rotations_dn, NCrit, RadCrit, CotCrit,                   & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                       & ! Inputs (Geometry)
     sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn,            & ! Inputs (Geometry)
     Stokes_dn, LC_Jacobians_dn, LP_Jacobians_dn  )                               ! Output
   endif

!  Finish

   return
end subroutine SSV_Integral_Plus_1_UPDN


subroutine DTEV_Integral_Plus_1_UP &
   ( maxlayers, maxfinelayers, max_user_levels, max_atmoswfs, max_surfacewfs, & ! Inputs (dimensioning)
     Do_Thermset, do_dmscaling, do_regular_ps, do_enhanced_ps,                & ! Inputs (Flags)
     doNadir, do_columnwfs, do_profilewfs, do_surfacewfs,                     & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,                          & ! Inputs (control output)
     n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums,                        & ! Inputs (Jacobian control)
     bb_input, surfbb, user_emissivity, extinction, deltaus, omega, truncfac, & ! Inputs (Optical)
     LS_user_emissivity, L_extinction, L_deltaus, L_omega, L_truncfac,        & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,                & ! Inputs (Geometry)
     Stokes_dta_up, Stokes_dts, tcom1, L_tcom1,                               & ! Outputs
     LC_Jacobians_dta_up, LP_Jacobians_dta_up,                                & ! Outputs
     LC_Jacobians_dts_up, LP_Jacobians_dts_up, LS_Jacobians_dts )               ! Output

!  Stand alone routine for Upwelling Direct-thermal-emission (DTEV)
!    computation of radiances and Jacobians. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012
!  This version, revised by R. Spurr, 21 June 2012

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

      real(fpk), Intent(Out)  :: Stokes_dta_up ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: Stokes_dts    ( max_user_levels, 4 )

      real(fpk), Intent(Out)  :: LC_Jacobians_dta_up  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dta_up  ( max_user_levels, 4, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LC_Jacobians_dts_up  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dts_up  ( max_user_levels, 4, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LS_Jacobians_dts     ( max_user_levels, 4, max_surfacewfs )

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
      real(fpk)  :: tms, L_tms(max_atmoswfs), lostau, L_thermcoeffs2

      real(fpk), parameter  :: cutoff = 88.0d0
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

      STOKES_dta_up = zero
      STOKES_dts    = zero

      LP_JACOBIANS_dta_UP = zero
      LC_JACOBIANS_dta_UP = zero
      LP_JACOBIANS_dts_UP = zero
      LC_JACOBIANS_dts_UP = zero
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
            else
               if ( Qvary(n) ) then
                  do q = 1, Qnums(n)
                     L_tms(q) = - L_omega(n,q)
                  enddo
               endif
            endif
            thermcoeffs(1)  = bb_input(n-1)
            thermcoeffs(2)  = (bb_input(n)-bb_input(n-1)) / deltaus(n)
            tcom1(n,1) = thermcoeffs(1) * tms
            tcom1(n,2) = thermcoeffs(2) * tms
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_thermcoeffs2 = - thermcoeffs(2) * L_deltaus(n,q) / deltaus(n)
                  L_tcom1(n,1,q) = thermcoeffs(1) * L_tms(q)
                  L_tcom1(n,2,q) = thermcoeffs(2) * L_tms(q) + L_thermcoeffs2 * tms
               enddo
            endif
         ENDDO
         do_Thermset = .false.
      endif

!  Plane/Parallel Layer integrated source terms
!  ============================================

      if ( do_regular_ps ) then
        if ( Mu1 .eq. zero ) then
          DO n = 1, nlayers
            if ( layermask_up(n) ) then
              t_mult_up(1:2) = tcom1(n,1:2)
              sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
              t_mult_up(0) = - sum
              sources_up(n) = t_mult_up(1)           ! Check this formula
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                   L_t_mult_up(1:2) = L_tcom1(n,1:2,q)
                   L_sum = L_t_mult_up(1) + t_mult_up(2) * L_deltaus(n,q) + L_t_mult_up(2) * deltaus(n)
                   L_t_mult_up(0) = - L_sum
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
                   L_t_mult_up(0) = - L_sum
                   L_sources_up(n,q) = L_t_mult_up(0) *   lostrans_up(n)   + &
                                         t_mult_up(0) * L_lostrans_up(n,q) + L_t_mult_up(1)
                enddo
!                if (n.gt.67)write(*,*)n,L_sources_up(n,1), sources_up(n)
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
                  wtrans_fine(n,j)    = ke * tcw
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_kn = L_extinction(n,q) ; L_tran = argum(j) * L_kn
                        L_solutions_fine(n,j,q) = L_tcom1(n,1,q) + xfine(n,j)*(kn*L_tcom1(n,2,q) + L_kn*tcom1(n,2))
                        L_wtrans_fine(n,j,q)    = ( L_kn - kn * L_tran ) * tcw * Raycon
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
         STOKES_DTA_UP(UTA,1) = CUMSOURCE_UP(NC)
         STOKES_DTS(UTA,1)    = CUMSOURCE_DSTE
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
               LS_Jacobians_dts(UTA,1,Q) = LS_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Profile Wfs (atmospheric emission terms)

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
                     LP_Jacobians_dta_up(UTA,1,K,Q) = L_CUMSOURCE(Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Column Wfs (Atmospheric emission terms)

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
               LC_Jacobians_dta_up(UTA,1,Q) = L_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Profile Wfs (surface emission term)

      if ( do_profilewfs ) then
         CUMSOURCE_DSTE   = SURFBB * USER_EMISSIVITY
         do k = 1, nlayers
            if ( Qvary(k) ) then
               L_CUMSOURCE = CUMSOURCE_DSTE
               NSTART = NLAYERS
               NUT_PREV = NSTART + 1
               DO UTA = N_USER_LEVELS, 1, -1
                  NUT    = USER_LEVELS(UTA) + 1
                  DO N = NSTART, NUT, -1
                     NC = NLAYERS + 1 - N
                     if ( k.eq.n ) then
                        do q = 1, Qnums(k)
                           L_cumsource(q) =  L_LOSTRANS_UP(N,Q) * L_CUMSOURCE(Q)
                        enddo
                     else
                        do q = 1, Qnums(k)
                           L_cumsource(q) = LOSTRANS_UP(N) * L_CUMSOURCE(Q)
                        enddo
                     endif
                  ENDDO
                  do q = 1, Qnums(k)
                     LP_Jacobians_dts_up(UTA,1,K,Q) = L_CUMSOURCE(Q)
                  enddo
                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                  NUT_PREV = NUT
               ENDDO
            endif
         enddo
      endif

!  Column Wfs (surface emission term)

      if ( do_columnwfs ) then
         CUMSOURCE_DSTE   = SURFBB * USER_EMISSIVITY
         L_CUMSOURCE     = CUMSOURCE_DSTE
          CUMSOURCE_DSTE   = LOSTRANS_UP(N) * CUMSOURCE_DSTE
         NSTART = NLAYERS
         NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N
               do q = 1, n_columnwfs
                  L_cumsource(q) =  L_LOSTRANS_UP(N,Q) * CUMSOURCE_DSTE + &
                                      LOSTRANS_UP(N)   * L_CUMSOURCE(Q)
                  CUMSOURCE_DSTE   = LOSTRANS_UP(N) * CUMSOURCE_DSTE
               enddo
            ENDDO
            do q = 1, n_columnwfs
               LC_Jacobians_dts_up(UTA,1,Q) = L_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Finish

      return
end subroutine DTEV_Integral_Plus_1_UP

!

subroutine DTEV_Integral_Plus_1_DN &
   ( maxlayers, maxfinelayers, max_user_levels, max_atmoswfs,                 & ! Inputs (dimensioning)
     Do_Thermset, do_dmscaling, do_regular_ps, do_enhanced_ps,                & ! Inputs (Flags)
     doNadir, do_columnwfs, do_profilewfs, nlayers, nfinedivs,                & ! Inputs (Flags/Numbers)
     n_user_levels, user_levels,  n_columnwfs, Lvaryflags, Lvarynums,         & ! Inputs (Numbers,Jacobian control)
     bb_input, extinction, deltaus, omega, truncfac,                          & ! Inputs (Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac,                            & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Radcrit, CotCrit, Raycon, cota,                              & ! Inputs (Geometry)
     radii, xfine, wfine, csqfine, cotfine,                                   & ! Inputs (Geometry)
     Stokes_dta_dn, tcom1, L_tcom1, LC_Jacobians_dta_dn, LP_Jacobians_dta_dn )    ! Output

!  Stand alone routine for Downwelling Direct-thermal-emission (DTE)
!    computation of radiances and Jacobians. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012
!  This version, revised by R. Spurr, 21 June 2012

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

      real(fpk), Intent(Out)  :: Stokes_dta_dn ( max_user_levels, 4 )

      real(fpk), Intent(Out)  :: LC_Jacobians_dta_dn  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dta_dn  ( max_user_levels, 4, maxlayers, max_atmoswfs )

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
      real(fpk)  :: tms, L_tms(max_atmoswfs), lostau, L_thermcoeffs2

      real(fpk), parameter  :: cutoff = 88.0d0
      real(fpk), parameter  :: zero   = 0.0_fpk
      real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

      STOKES_dta_dn = zero

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
                  L_thermcoeffs2 = - thermcoeffs(2) * L_deltaus(n,q) / deltaus(n)
                  L_tcom1(n,1,q) = thermcoeffs(1) * L_tms(q)
                  L_tcom1(n,2,q) = thermcoeffs(2) * L_tms(q) + L_thermcoeffs2 * tms
               enddo
            endif
         ENDDO
         do_Thermset = .false.
      endif

!  Plane/Parallel Layer integrated source terms
!  ============================================

      if ( do_regular_ps ) then
        if ( Mu1 .eq. zero  ) then
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
                  wtrans_fine(n,j)    = ke * tcw
                  if ( Qvary(n) ) then
                     do q = 1, Qnums(n)
                        L_kn = L_extinction(n,q) ; L_tran = argum(j) * L_kn
                        L_solutions_fine(n,j,q) = L_tcom1(n,1,q) + xfine(n,j)*(kn*L_tcom1(n,2,q) + L_kn*tcom1(n,2))
                        L_wtrans_fine(n,j,q)    = ( L_kn - kn * L_tran ) * tcw * Raycon
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
         STOKES_DTA_DN(UTA,1) = CUMSOURCE_DN(NC)
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
                     LP_Jacobians_dta_DN(UTA,1,K,Q) = L_CUMSOURCE(Q)
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
               LC_Jacobians_dta_dn(UTA,1,Q) = L_CUMSOURCE(Q)
            enddo
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Finish

      return
end subroutine DTEV_Integral_Plus_1_DN

!

subroutine DTEV_Integral_Plus_1_UPDN &
   ( maxlayers, maxfinelayers, max_user_levels, max_atmoswfs, max_surfacewfs, & ! Inputs (dimensioning)
     Do_Thermset, do_upwelling, do_dnwelling, do_dmscaling, do_regular_ps,    & ! Inputs (Flags)
     do_enhanced_ps, doNadir, do_columnwfs, do_profilewfs, do_surfacewfs,     & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,                          & ! Inputs (control output)
     n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums,                        & ! Inputs (Jacobian control)
     bb_input, surfbb, user_emissivity, extinction, deltaus, omega, truncfac, & ! Inputs (Optical)
     LS_user_emissivity, L_extinction, L_deltaus, L_omega, L_truncfac,        & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Radcrit, CotCrit, Raycon, cota,                              & ! Inputs (Geometry)
     radii, xfine, wfine, csqfine, cotfine,                                   & ! Inputs (Geometry)
     Stokes_dta_up, Stokes_dta_dn, Stokes_dts, tcom1, L_tcom1,                & ! Outputs
     LC_Jacobians_dta_up, LP_Jacobians_dta_up,                                & ! Output
     LC_Jacobians_dta_dn, LP_Jacobians_dta_dn,                                & ! Outputs
     LC_Jacobians_dts_up, LP_Jacobians_dts_up, LS_Jacobians_dts )               ! Output

!  Stand alone routine for Upwelling and downwelling Direct-thermal-emission (DTE)
!    computation of radiances and Jacobians. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012
!  This version, revised by R. Spurr, 21 June 2012

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

      real(fpk), Intent(Out)  :: Stokes_dta_dn ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: Stokes_dta_up ( max_user_levels, 4 )
      real(fpk), Intent(Out)  :: Stokes_dts    ( max_user_levels, 4 )

      real(fpk), Intent(Out)  :: LC_Jacobians_dta_dn  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dta_dn  ( max_user_levels, 4, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LC_Jacobians_dta_up  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dta_up  ( max_user_levels, 4, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LC_Jacobians_dts_up  ( max_user_levels, 4, max_atmoswfs )
      real(fpk), Intent(Out)  :: LP_Jacobians_dts_up  ( max_user_levels, 4, maxlayers, max_atmoswfs )
      real(fpk), Intent(Out)  :: LS_Jacobians_dts     ( max_user_levels, 4, max_surfacewfs )

!  Thermal setup

      real(fpk), Intent(InOut)   :: tcom1(maxlayers,2)
      real(fpk), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)


   if ( do_upwelling ) then
      call DTEV_Integral_Plus_1_UP &
   ( maxlayers, maxfinelayers, max_user_levels, max_atmoswfs, max_surfacewfs, & ! Inputs (dimensioning)
     Do_Thermset, do_dmscaling, do_regular_ps, do_enhanced_ps,                & ! Inputs (Flags)
     doNadir, do_columnwfs, do_profilewfs, do_surfacewfs,                     & ! Inputs (Flags)
     nlayers, nfinedivs, n_user_levels, user_levels,                          & ! Inputs (control output)
     n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums,                        & ! Inputs (Jacobian control)
     bb_input, surfbb, user_emissivity, extinction, deltaus, omega, truncfac, & ! Inputs (Optical)
     LS_user_emissivity, L_extinction, L_deltaus, L_omega, L_truncfac,        & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,                & ! Inputs (Geometry)
     Stokes_dta_up, Stokes_dts, tcom1, L_tcom1,                               & ! Outputs
     LC_Jacobians_dta_up, LP_Jacobians_dta_up,                                & ! Output
     LC_Jacobians_dts_up, LP_Jacobians_dts_up, LS_Jacobians_dts )                     ! Output
   endif

   if ( do_dnwelling ) then
      call DTEV_Integral_Plus_1_DN &
   ( maxlayers, maxfinelayers, max_user_levels, max_atmoswfs,                 & ! Inputs (dimensioning)
     Do_Thermset, do_dmscaling, do_regular_ps, do_enhanced_ps,                & ! Inputs (Flags)
     doNadir, do_columnwfs, do_profilewfs, nlayers, nfinedivs,                & ! Inputs (Flags/Numbers)
     n_user_levels, user_levels,  n_columnwfs, Lvaryflags, Lvarynums,         & ! Inputs (Numbers,Jacobian control)
     bb_input, extinction, deltaus, omega, truncfac,                          & ! Inputs (Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac,                            & ! Inputs (Optical - Linearized)
     Mu1, NCrit, Radcrit, CotCrit, Raycon, cota,                              & ! Inputs (Geometry)
     radii, xfine, wfine, csqfine, cotfine,                                   & ! Inputs (Geometry)
     Stokes_dta_dn, tcom1, L_tcom1, LC_Jacobians_dta_dn, LP_Jacobians_dta_dn )    ! Output
   endif

!  Finish

   return
end subroutine DTEV_Integral_Plus_1_UPDN

!  End module

end module FO_Vector_RTCalcs_Plus_m


