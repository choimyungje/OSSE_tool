
module FO_Vector_SS_Master_Plus_m

!  Interface based on 10 May 2012 FO code package.
!  Mk1, developed 14 June with TOA-only fixed inputs
!  Mk2, developed 20 June with more flexible inputs
!  Mk3, developed 21 June with update FO code from June4-8, 2012

!  Restricted to Lambertian only

   use FO_geometry_SSonly_m, only : SS_Geometry_1
   use FO_Vector_spherfuncs_m
   use FO_Vector_RTCalcs_Plus_m      , only : SSV_Integral_Plus_1_UP, SSV_Integral_Plus_1_DN

contains

subroutine FO_Vector_SS_Master_Plus &
        ( FO_maxfine, FO_nfineinput, FO_do_regular_ps, FO_ACrit,         & ! FO-specific inputs
          do_deltam_scaling, do_upwelling, do_dnwelling,                 & ! VLIDORT proxy inputs (flags)
          do_columnwfs, do_profilewfs, do_surfacewfs,                    & ! VLIDORT proxy inputs (Lin flags)
          nlayers, nstokes, nmoments, nstreams, nuserlevels, userlevels, & ! VLIDORT proxy inputs (control)
          n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,   & ! VLIDORT proxy inputs (Lin control)
          alpha_boa, theta_boa, phi_boa, eradius, heights,               & ! VLIDORT proxy inputs (geometry)
          Fluxfac, lambertian_albedo, deltaus, omegas, greekmat,         & ! VLIDORT proxy inputs (optical)
          L_deltaus, L_omegas, L_greekmat,                               & ! VLIDORT proxy inputs (Lin optical)
          FO_stokes_ss, FO_stokes_db, FO_LC_Jacobians_ss, FO_LP_Jacobians_ss,    & ! Output
          FO_LC_Jacobians_db, FO_LP_Jacobians_db, FO_LS_Jacobians_db,            & ! Output
          FO_fail, FO_message, FO_trace )                                          ! Output

!  VLIDORT parameter file

   use VLIDORT_PARS

!  Implicit none

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Dimensions - set them to the VLIDORT values (not finelayers)

   integer, parameter :: FO_maxlayers     = maxlayers
   integer, parameter :: FO_maxmoments    = maxmoments_input
   INTEGER, parameter :: FO_maxuserlevels = max_user_levels

   integer, parameter :: FO_max_atmoswfs   = max_atmoswfs
   integer, parameter :: FO_max_surfacewfs = max_surfacewfs

!  Inputs (FO-specific)
!  ====================

!  Fine layer dimensioning (should be independent of VLIDORT quantity)

   integer, intent(in) :: FO_maxfine

!  FO Regular psuedo-spherical flag

   logical, intent(in) :: FO_do_regular_ps

!  FO finelayer input, starting value

   integer, intent(in) :: FO_nfineinput

!  FO Critical attenuation

   REAL(ffp), Intent(in) ::  FO_ACrit

!  Inputs (VLIDORT PROXIES)
!  ========================

!  Control and numbers

   logical  , Intent(in) ::  do_deltam_scaling, do_upwelling, do_dnwelling
   integer  , Intent(in) ::  nstokes, nlayers, nmoments, nstreams, nuserlevels
   REAL(ffp), Intent(in) ::  eradius, heights(0:MAXLAYERS), userlevels(max_user_levels)
   REAL(ffp), Intent(in) ::  alpha_boa, theta_boa, phi_boa

!  Linearization Control and numbers

   logICAL, Intent(in) ::  do_surfacewfs
   logICAL, Intent(in) ::  do_columnwfs
   logICAL, Intent(in) ::  do_profilewfs

   INTEGER, Intent(in) ::  n_columnwfs
   INTEGER, Intent(in) ::  n_surfacewfs
   logICAL, Intent(in) ::  Lvaryflags(maxlayers)
   INTEGER, Intent(in) ::  Lvarynums (maxlayers)
   logICAL, Intent(in) ::  Lvarymoms (maxlayers,max_atmoswfs)

!  Optical inputs

   REAL(ffp), Intent(in) :: Fluxfac
   REAL(ffp), Intent(in) :: lambertian_albedo
   REAL(ffp), Intent(in) :: deltaus    ( MAXLAYERS )
   REAL(ffp), Intent(in) :: omegas     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: greekmat   ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Linearized optical inputs

   REAL(ffp), Intent(in) :: L_DELTAUS     ( max_atmoswfs, MAXLAYERS )
   REAL(ffp), Intent(in) :: L_OMEGAS      ( max_atmoswfs, MAXLAYERS )
   REAL(ffp), Intent(in) :: L_GREEKMAT    ( max_atmoswfs, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Outputs
!  =======

!  Stokes vectors

   real(ffp), intent(out)  :: FO_stokes_ss     ( max_user_levels, 4 )
   real(ffp), intent(out)  :: FO_stokes_db     ( max_user_levels, 4 )

!  Linearized SS Stokes vectors (Up or Down)

   real(ffp), Intent(Out)  :: FO_LC_Jacobians_ss  ( max_user_levels, 4, max_atmoswfs )
   real(ffp), Intent(Out)  :: FO_LP_Jacobians_ss  ( max_user_levels, 4, maxlayers, max_atmoswfs )

!  DB dependencies (Up only)

   real(ffp), Intent(Out)  :: FO_LC_Jacobians_db  ( max_user_levels, 4, max_atmoswfs )
   real(ffp), Intent(Out)  :: FO_LP_Jacobians_db  ( max_user_levels, 4, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: FO_LS_Jacobians_db  ( max_user_levels, 4, max_surfacewfs )

!  Exception handling from Geometry routine

   logical      , intent(out)   :: FO_fail
   character*(*), intent(out)   :: FO_message
   character*(*), intent(out)   :: FO_trace

!  Local
!  =====

!  Layer control. Finelayer divisions may be changed

   integer    :: FO_nlayers
   integer    :: FO_nfinedivs(FO_maxlayers)

!  Number of user levels

   integer    :: FO_n_user_levels
   integer    :: FO_user_levels(FO_maxuserlevels)

!  Radius + heights

   real(ffp)  :: FO_eradius, FO_heights (0:FO_maxlayers)

!  input angles (Degrees), dtr = degrees-to-Radians. VSIGN = +1 (Up); -1(Down)

   real(ffp)  :: FO_dtr, FO_Pie, FO_vsign
   real(ffp)  :: FO_alpha_boa, FO_theta_boa, FO_phi_boa

!  Geometry Flags

   logical    :: FO_do_Chapman
   logical    :: FO_do_LOSpaths

!  Derived flags, optical settings

   logical    :: FO_do_enhanced_ps
   logical    :: FO_do_planpar
   logical    :: FO_do_deltam_scaling
   logical    :: FO_do_lambertian

!  Linearized

   logICAL    :: FO_do_surfacewfs
   logICAL    :: FO_do_columnwfs
   logICAL    :: FO_do_profilewfs

   INTEGER    :: FO_n_columnwfs
   INTEGER    :: FO_n_surfacewfs
   logICAL    :: FO_Lvaryflags(FO_maxlayers)
   INTEGER    :: FO_Lvarynums (FO_maxlayers)
   logICAL    :: FO_Lvarymoms (FO_maxlayers,FO_max_atmoswfs)

!  Critical adjustment for cloud layers

   logical    :: FO_doCrit
   real(ffp)  :: FO_extincs(FO_maxlayers)

!  Flag for the Nadir case. Intent(inout), input if DO_LOSpaths set)

   logical    :: FO_doNadir
 
!  Alphas,  Cotangents, Radii, Ray constant. Intent(inout), input if DO_LOSpaths set)

   real(ffp)  :: FO_radii    (0:FO_maxlayers)
   real(ffp)  :: FO_Raycon
   real(ffp)  :: FO_alpha    (0:FO_maxlayers)
   real(ffp)  :: FO_cota     (0:FO_maxlayers)

!  LOS Quadratures for Enhanced PS, input if DO_LOSpaths set)

   real(ffp)  :: FO_xfine    (FO_maxlayers,FO_maxfine)
   real(ffp)  :: FO_wfine    (FO_maxlayers,FO_maxfine)
   real(ffp)  :: FO_csqfine  (FO_maxlayers,FO_maxfine)
   real(ffp)  :: FO_cotfine  (FO_maxlayers,FO_maxfine)

!  Fine layering output, input if DO_LOSpaths set)

   real(ffp)  :: FO_alphafine (FO_maxlayers,FO_maxfine)
   real(ffp)  :: FO_radiifine (FO_maxlayers,FO_maxfine)

!  Critical layer (Intent out)

   integer   :: FO_Ncrit
   real(ffp) :: FO_RadCrit, FO_CotCrit

!  solar paths, Intent out

   integer   :: FO_ntraverse  (0:FO_maxlayers)
   real(ffp) :: FO_sunpaths   (0:FO_maxlayers,FO_maxlayers)
   integer   :: FO_ntraverse_fine(FO_maxlayers,FO_maxfine)
   real(ffp) :: FO_sunpaths_fine (FO_maxlayers,FO_maxlayers,FO_maxfine)
   real(ffp) :: FO_Chapfacs   (FO_maxlayers,FO_maxlayers)

!  Cosine scattering angle and other cosines

   real(ffp) :: FO_cosscat
   real(ffp) :: FO_Mu1
   real(ffp) :: FO_Mu0

!  SPherical Functions: Inputs. [Starter flag may be re-set]

   logICAL   :: FO_STARTER
   logICAL   :: FO_DO_SUNLIGHT
   intEGER   :: FO_NMOMENTS, FO_NSTOKES

!  SPherical Functions: Outputs

   real(ffp) :: FO_ROTATIONS(4)
   real(ffp) :: FO_GENSPHER(0:FO_MAXMOMENTS,4)
   real(ffp) :: FO_GSHELP(7,FO_MAXMOMENTS)

!  optical inputs Atmosphere

   real(ffp) :: FO_DELTAUS     ( FO_MAXLAYERS )
   real(ffp) :: FO_OMEGAS      ( FO_MAXLAYERS )
   real(ffp) :: FO_TRUNCFAC    ( FO_MAXLAYERS )
   real(ffp) :: FO_GREEKMAT    ( FO_MAXLAYERS, 0:FO_MAXMOMENTS, 4, 4 )

!  linearized optical inputs (Atmosphere)

   REAL(fpk) :: FO_L_DELTAUS  ( FO_MAXLAYERS, FO_max_atmoswfs )
   REAL(fpk) :: FO_L_EXTINCS  ( FO_MAXLAYERS, FO_max_atmoswfs )
   REAL(ffp) :: FO_L_OMEGAS   ( FO_MAXLAYERS, FO_max_atmoswfs )
   REAL(ffp) :: FO_L_TRUNCFAC ( FO_MAXLAYERS, FO_max_atmoswfs )
   REAL(ffp) :: FO_L_GREEKMAT ( FO_MAXLAYERS, 0:FO_MAXMOMENTS, 4, 4, FO_max_atmoswfs )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(ffp) :: FO_REFLEC(4,4), FO_FLUX, FO_FLUXVEC(4)
   REAL(ffp) :: FO_LS_REFLEC(4,4,FO_max_surfacewfs )

!  First order output

!   real(ffp) :: FO_cumsource  ( 0:FO_maxlayers, 4 )

!  Local variables
!  ===============

   logical   :: loop
   integer   :: NA, L, M1, O1, O2, MM, Q, M21, NM1
   real(ffp) :: scale, diff, mdiff, L_scale

!  settings from VLIDORT Proxies
!  =============================

!  numbers

   FO_nlayers  = nlayers
   FO_nmoments = nmoments
   FO_nstokes  = nstokes

!  Heights and eradius

   FO_eradius  = eradius
   FO_heights(0:nlayers) = heights(0:nlayers)

!  flags

   FO_do_deltam_scaling  = DO_DELTAM_SCALING

!  Output levels (integer-valued in the FO code, real-valued in VLIDORT)

   FO_n_user_levels = nuserlevels
   do mm = 1, nuserlevels
      FO_user_levels(mm) = int(userlevels(mm))
   enddo

!  Trig and Flux

!   FO_Pie = acos(-1.0_ffp)   ; FO_dtr = FO_Pie / 180.0_ffp

   FO_Pie = PIE              ; FO_dtr = DEG_TO_RAD
   FO_FLUX = 0.25_ffp/FO_Pie ; FO_do_sunlight = .true.
   FO_FLUXVEC = 0.0_ffp      ; FO_FLUXVEC(1) = Fluxfac

!  Linearization control

   FO_lvaryflags(1:nlayers) = lvaryflags(1:nlayers)
   FO_lvarynums (1:nlayers) = lvarynums (1:nlayers)

   FO_n_columnwfs   = n_columnwfs
   FO_do_columnwfs  = do_columnwfs
   FO_do_profilewfs = do_profilewfs
   FO_do_surfacewfs = do_surfacewfs 
   FO_n_surfacewfs  = n_surfacewfs

!  Set Optical properties from VLIDORT proxies
!  ===========================================

!  Initialize

   FO_REFLEC   = 0.0_ffp ; FO_REFLEC(1,1) = lambertian_albedo
   FO_TRUNCFAC = 0.0_ffp ; FO_GREEKMAT = 0.0_ffp
   if ( FO_do_deltam_scaling ) then
      M21 = 2*nstreams ; mdiff = real(2*M21+1,ffp)
   endif

!  regular properties

   do na = 1, nlayers
      FO_OMEGAS(na)  = omegas(na)
      if ( FO_do_deltam_scaling ) FO_TRUNCFAC(na) = greekmat(M21,na,1) / mdiff
      scale = 1.0_ffp - FO_TRUNCFAC(na)*FO_omegas(na)
      FO_DELTAUS(na)  = scale * deltaus(na)
      diff = FO_heights(na-1) - FO_heights(na)
      FO_EXTINCS(na) = FO_DELTAUS(na) / diff
      do L = 0, nmoments
         DO O1 = 1, nstokes
            M1 = 4 * ( O1 - 1 )
            DO O2 = 1, nstokes
               MM = M1 + O2
               FO_GREEKMAT(NA,L,O1,O2) = greekmat(L,na,MM) 
            enddo
         enddo
      enddo
   enddo

!  Linearized properties
!    --- FO linearized stuff is Single normalized
!    --- FO Lvary moms is computed from scratch.

   FO_L_TRUNCFAC = 0.0_ffp ; FO_lvarymoms  = .false.
   nm1 = 2*nstreams-1 ; if ( FO_do_deltam_scaling ) nm1 = nm1 + 1

!mick fix
   FO_L_EXTINCS = 0.0_ffp

   do na = 1, nlayers
     if ( lvaryflags(na) ) then
       do q = 1, lvarynums(na)
         LOOP = .TRUE. ;  L = 0
         DO WHILE (LOOP.AND.L.LT.NM1)
           L = L + 1 ; loop = (abs(L_greekmat(q,L,na,1)).lt.1.0d-08)
         ENDDO
         FO_lvarymoms(na,q) = .not.loop
         FO_L_OMEGAS(na,q)  = L_omegas(q,na) * omegas(na)
         if ( FO_do_deltam_scaling ) FO_L_TRUNCFAC(na,q) = L_greekmat(q,M1,na,1) * FO_TRUNCFAC(na)
         L_scale = FO_L_TRUNCFAC(na,q)*FO_omegas(na) - FO_TRUNCFAC(na)*FO_L_omegas(na,q)
         FO_L_DELTAUS(na,q)  = ( L_scale + scale * L_deltaus(q,na) ) * deltaus(na)
         diff = FO_heights(na-1) - FO_heights(na)
         FO_L_EXTINCS(na,q) = FO_L_DELTAUS(na,q) / diff
         do L = 0, nmoments
           DO O1 = 1, nstokes
             M1 = 4 * ( O1 - 1 )
             DO O2 = 1, nstokes
               MM = M1 + O2
               FO_L_GREEKMAT(NA,L,O1,O2,q) = L_greekmat(q,L,na,MM) * greekmat(L,na,MM) 
             enddo
           enddo
         enddo
       enddo
     endif
   enddo

!  initial settings
!  ----------------

!  LOSPATHS variable should be false

   FO_do_LOSPATHS = .false.
   FO_starter     = .true.

!  other fixed flags

   FO_do_lambertian = .true.
   FO_do_planpar    = .false.

!  pseudo-spherical flag

   FO_do_enhanced_PS     = .not.FO_do_regular_ps

!  Angles

   FO_alpha_boa   = alpha_boa 
   FO_Mu1         = cos(FO_alpha_boa * FO_dtr)
   FO_theta_boa   = theta_boa
   FO_Mu0         = cos(FO_theta_boa * FO_dtr)
   FO_phi_boa     = phi_boa

!  Scattering angle control

   if ( do_upwelling ) FO_vsign =  1.0_ffp
   if ( do_dnwelling ) FO_vsign = -1.0_ffp

!  Attenuation control

   FO_nfinedivs = FO_nfineinput                ! Now an input
   FO_doCrit    = .not.FO_do_regular_ps
!   FO_Acrit     = 1.0e-10_ffp                 ! Now an input

!  Geometry call
!  -------------

   call SS_Geometry_1 &
              ( FO_maxlayers, FO_maxfine, FO_nlayers, FO_nfinedivs, FO_dtr, FO_Pie, FO_vsign,     & ! Input
                FO_eradius, FO_heights, FO_alpha_boa, FO_theta_boa, FO_phi_boa, FO_do_Chapman,    & ! Input
                FO_do_LOSpaths, FO_do_planpar, FO_do_regular_ps, FO_doCrit, FO_Acrit, FO_extincs, & ! Input
                FO_doNadir, FO_Raycon, FO_radii, FO_alpha, FO_cota,                               & ! Input/Output
                FO_xfine, FO_wfine, FO_csqfine, FO_cotfine, FO_alphafine, FO_radiifine,           & ! Input/Output
                FO_NCrit, FO_RadCrit, FO_CotCrit, FO_cosscat, FO_chapfacs,                        & ! Output
                FO_sunpaths, FO_ntraverse, FO_sunpaths_fine, FO_ntraverse_fine,                   & ! Output
                FO_fail, FO_message, FO_trace )                                                     ! Output

!  Exception handling on geometry

   if ( FO_Fail ) then
      write(*,*)FO_message ; write(*,*)FO_trace ; return
   endif

!  Vector Spherical functions
!  --------------------------

   CALL FO_Vector_spherfuncs &
              ( FO_STARTER, FO_MAXMOMENTS, FO_NMOMENTS, FO_NSTOKES, FO_DTR,                   & ! Inputs
                FO_DO_SUNLIGHT, FO_THETA_BOA, FO_ALPHA_BOA, FO_PHI_BOA, FO_COSSCAT, FO_VSIGN, & ! Inputs
                FO_ROTATIONS, FO_GSHELP, FO_GENSPHER )                                          ! Outputs

!  Do the source function integration. UPWELLING
!  =============================================

   if ( do_upwelling ) then
      call SSV_Integral_Plus_1_UP &
           ( FO_maxlayers, FO_maxfine, FO_maxmoments, FO_maxuserlevels, FO_do_sunlight,                       & ! Inputs (dimension)
             FO_max_atmoswfs, FO_max_surfacewfs, FO_do_columnwfs, FO_do_profilewfs, FO_do_surfacewfs,         & ! Inputs (dimensioning/flags)
             FO_do_deltam_scaling, FO_do_regular_ps, FO_do_enhanced_ps, FO_do_lambertian, FO_doNadir,         & ! Inputs (Flags)
             FO_nlayers, FO_nfinedivs, FO_nmoments, FO_n_user_levels, FO_user_levels,                         & ! Inputs (control)
             FO_nstokes, FO_n_columnwfs, FO_n_surfacewfs, FO_Lvaryflags, FO_Lvarynums, FO_Lvarymoms,          & ! Inputs (Jacobian control)
             FO_reflec, FO_extincs, FO_deltaus, FO_omegas, FO_truncfac, FO_greekmat, FO_flux, FO_fluxvec,     & ! Inputs (Optical)
             FO_LS_reflec, FO_L_extincs, FO_L_deltaus, FO_L_omegas, FO_L_truncfac, FO_L_greekmat,             & ! Inputs (Optical - Linearized)
             FO_Mu0, FO_Mu1, FO_GenSpher, FO_Rotations, FO_NCrit, FO_xfine, FO_wfine, FO_csqfine, FO_cotfine, & ! Inputs (Geometry)
             FO_Raycon, FO_cota, FO_sunpaths, FO_ntraverse, FO_sunpaths_fine, FO_ntraverse_fine,              & ! Inputs (Geometry)
             FO_stokes_ss, FO_stokes_db, FO_LC_Jacobians_ss, FO_LP_Jacobians_ss,                              & ! Output
             FO_LC_Jacobians_db, FO_LP_Jacobians_db, FO_LS_Jacobians_db )                                       ! Output
   endif

!  Do the source function integration. DNWELLING
!  =============================================

   if ( do_dnwelling ) then
      call SSV_Integral_Plus_1_DN &
           ( FO_maxlayers, FO_maxfine, FO_maxmoments, FO_maxuserlevels, FO_do_sunlight,            & ! Inputs (dimension)
             FO_max_atmoswfs, FO_do_columnwfs, FO_do_profilewfs,                                   & ! Inputs (dimensioning/flags)
             FO_do_deltam_scaling, FO_do_regular_ps, FO_do_enhanced_ps, FO_doNadir,                & ! Inputs (Flags)
             FO_nlayers, FO_nfinedivs, FO_nmoments, FO_n_user_levels, FO_user_levels,              & ! Inputs (control)
             FO_nstokes, FO_n_columnwfs, FO_Lvaryflags, FO_Lvarynums, FO_Lvarymoms,                & ! Inputs (Jacobian control)
             FO_extincs, FO_deltaus, FO_omegas, FO_truncfac, FO_greekmat, FO_flux, FO_fluxvec,     & ! Inputs (Optical)
             FO_L_extincs, FO_L_deltaus, FO_L_omegas, FO_L_truncfac, FO_L_greekmat,                & ! Inputs (Optical - Linearized)
             FO_Mu1, FO_GenSpher, FO_Rotations, FO_NCrit, FO_RadCrit, FO_CotCrit,                  & ! Inputs (Geometry)
             FO_xfine, FO_wfine, FO_csqfine, FO_cotfine, FO_Raycon, FO_radii, FO_cota,             & ! Inputs (Geometry)
             FO_sunpaths, FO_ntraverse, FO_sunpaths_fine, FO_ntraverse_fine,                       & ! Inputs (Geometry)
             FO_stokes_ss, FO_LC_Jacobians_ss, FO_LP_Jacobians_ss )                                  ! Output
   endif

!  Finish

   return
end subroutine FO_Vector_SS_Master_Plus

!  Finish module

end module FO_Vector_SS_Master_Plus_m

