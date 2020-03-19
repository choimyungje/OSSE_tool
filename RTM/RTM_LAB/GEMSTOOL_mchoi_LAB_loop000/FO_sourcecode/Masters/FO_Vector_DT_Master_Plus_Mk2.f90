
module FO_Vector_DT_Master_Plus_m

!  Interface based on 10 May 2012 FO code package.
!  Mk1, developed 14 June with TOA-only fixed inputs
!  Mk2, developed 20 June with more flexible inputs
!  Mk3, developed 21 June with update FO code from June4-8, 2012

!  Restricted to Lambertian only

   use FO_geometry_DTonly_m
   use FO_Vector_RTCalcs_Plus_m , only : DTEV_Integral_Plus_1_UP, DTEV_Integral_Plus_1_DN

contains

subroutine FO_Vector_DT_Master_Plus &
        ( FO_maxfine, FO_nfineinput, FO_do_regular_ps, FO_ACrit,          & ! FO-specific inputs
          do_deltam_scaling, do_upwelling, do_dnwelling,                  & ! VLIDORT proxy inputs (flags)
          do_columnwfs, do_profilewfs, do_surfacewfs,                     & ! VLIDORT proxy inputs (Lin flags)
          nlayers, nstreams, nuserlevels, userlevels,                     & ! VLIDORT proxy inputs (control)
          n_columnwfs, n_surfacewfs, Lvaryflags, Lvarynums, alpha_boa,    & ! VLIDORT proxy inputs (Lin control)
          eradius, heights, bb_input, surfbb, user_emissivity, LS_user_emissivity,    & ! VLIDORT proxy inputs (Geom/Emission)
          deltaus, omegas, greekmat, L_deltaus, L_omegas, L_greekmat,                 & ! VLIDORT proxy inputs (optical)
          FO_stokes_dta, FO_stokes_dts, FO_LC_Jacobians_dta, FO_LP_Jacobians_dta,     & ! Output
          FO_LC_Jacobians_dts, FO_LP_Jacobians_dts, FO_LS_Jacobians_dts,              & ! Output
          FO_fail, FO_message, FO_trace )                                               ! Output

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
   integer  , Intent(in) ::  nlayers, nstreams, nuserlevels
   REAL(ffp), Intent(in) ::  eradius, heights(0:MAXLAYERS), userlevels(max_user_levels)
   REAL(ffp), Intent(in) ::  alpha_boa

!  Linearization Control and numbers

   logICAL, Intent(in) ::  do_surfacewfs
   logICAL, Intent(in) ::  do_columnwfs
   logICAL, Intent(in) ::  do_profilewfs

   INTEGER, Intent(in) ::  n_columnwfs
   INTEGER, Intent(in) ::  n_surfacewfs
   logICAL, Intent(in) ::  Lvaryflags(maxlayers)
   INTEGER, Intent(in) ::  Lvarynums (maxlayers)

!  Optical inputs

   REAL(ffp), Intent(in) :: deltaus    ( MAXLAYERS )
   REAL(ffp), Intent(in) :: omegas     ( MAXLAYERS )
   REAL(ffp), Intent(in) :: greekmat   ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Linearized optical inputs

   REAL(ffp), Intent(in) :: L_DELTAUS     ( max_atmoswfs, MAXLAYERS )
   REAL(ffp), Intent(in) :: L_OMEGAS      ( max_atmoswfs, MAXLAYERS )
   REAL(ffp), Intent(in) :: L_GREEKMAT    ( max_atmoswfs, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Emission inputs

   REAL(ffp), Intent(in) :: bb_input    ( 0:MAXLAYERS )
   REAL(ffp), Intent(in) :: surfbb, user_emissivity, LS_user_emissivity (MAX_SURFACEWFS)

!  Outputs
!  =======

!  Stokes vectors

   real(ffp), intent(out)  :: FO_stokes_dta     ( max_user_levels, 4 )
   real(ffp), intent(out)  :: FO_stokes_dts     ( max_user_levels, 4 )

!  Linearized SS Stokes vectors (Up or Down)

   real(ffp), Intent(Out)  :: FO_LC_Jacobians_dta  ( max_user_levels, 4, max_atmoswfs )
   real(ffp), Intent(Out)  :: FO_LP_Jacobians_dta  ( max_user_levels, 4, maxlayers, max_atmoswfs )

!  DB dependencies (Up only)

   real(ffp), Intent(Out)  :: FO_LC_Jacobians_dts  ( max_user_levels, 4, max_atmoswfs )
   real(ffp), Intent(Out)  :: FO_LP_Jacobians_dts  ( max_user_levels, 4, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: FO_LS_Jacobians_dts  ( max_user_levels, 4, max_surfacewfs )

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

!  input angles (Degrees), dtr = degrees-to-Radians. 

   real(ffp)  :: FO_dtr, FO_alpha_boa

!  Geometry Flags

   logical    :: FO_do_LOSpaths

!  Derived flags, optical settings

   logical    :: FO_do_enhanced_ps
   logical    :: FO_do_deltam_scaling
   logical    :: FO_do_Thermset

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

!  other cosines

   real(ffp) :: FO_Mu1

!  optical inputs Atmosphere

   real(ffp) :: FO_DELTAUS     ( FO_MAXLAYERS )
   real(ffp) :: FO_OMEGAS      ( FO_MAXLAYERS )
   real(ffp) :: FO_TRUNCFAC    ( FO_MAXLAYERS )

!  linearized optical inputs (Atmosphere)

   REAL(fpk) :: FO_L_DELTAUS  ( FO_MAXLAYERS, FO_max_atmoswfs )
   REAL(fpk) :: FO_L_EXTINCS  ( FO_MAXLAYERS, FO_max_atmoswfs )
   REAL(ffp) :: FO_L_OMEGAS   ( FO_MAXLAYERS, FO_max_atmoswfs )
   REAL(ffp) :: FO_L_TRUNCFAC ( FO_MAXLAYERS, FO_max_atmoswfs )

!  Emission inputs

   real(ffp) :: FO_BB_INPUT (0:FO_MAXLAYERS)
   real(ffp) :: FO_SURFBB, FO_USER_EMISSIVITY
   real(ffp) :: FO_LS_USER_EMISSIVITY(FO_max_surfacewfs)

!  Thermal setup

   real(ffp) :: FO_tcom1(FO_maxlayers,2)
   real(ffp) :: FO_L_tcom1(FO_maxlayers,2,FO_max_atmoswfs)

!  Local variables
!  ===============

   logical   :: loop
   integer   :: NA, MM, Q, M21, NM1, L
   real(ffp) :: scale, diff, mdiff, L_scale

!  settings from VLIDORT Proxies
!  =============================

!  Thermset flag (true here)

   FO_do_Thermset = .true.

!  numbers

   FO_nlayers  = nlayers

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

   FO_dtr = DEG_TO_RAD

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

   FO_TRUNCFAC = 0.0_ffp
   if ( FO_do_deltam_scaling ) then
      M21 = 2*nstreams ; mdiff = real(2*M21+1,ffp)
   endif

!  Emission

   FO_surfbb   = SURFBB
   FO_bb_input = BB_INPUT
   FO_user_emissivity    = USER_EMISSIVITY
   FO_LS_user_emissivity = LS_USER_EMISSIVITY

!  regular properties

   do na = 1, nlayers
      FO_OMEGAS(na)  = omegas(na)
      if ( FO_do_deltam_scaling ) FO_TRUNCFAC(na) = greekmat(M21,na,1) / mdiff
      scale = 1.0_ffp - FO_TRUNCFAC(na)*FO_omegas(na)
      FO_DELTAUS(na)  = scale * deltaus(na)
      diff = FO_heights(na-1) - FO_heights(na)
      FO_EXTINCS(na) = FO_DELTAUS(na) / diff
   enddo

!  Linearized properties
!    --- FO linearized stuff is Single normalized
!    --- FO Lvary moms is computed from scratch.

   FO_L_TRUNCFAC = 0.0_ffp ; FO_lvarymoms  = .false.
   nm1 = 2*nstreams-1 ; if ( FO_do_deltam_scaling ) nm1 = nm1 + 1

   do na = 1, nlayers
     if ( lvaryflags(na) ) then
       do q = 1, lvarynums(na)
         LOOP = .TRUE. ;  L = 0
         DO WHILE (LOOP.AND.L.LT.NM1)
           L = L + 1 ; loop = (abs(L_greekmat(q,L,na,1)).lt.1.0d-08)
         ENDDO
         FO_lvarymoms(na,q) = .not.loop
         FO_L_OMEGAS(na,q)  = L_omegas(q,na) * omegas(na)
         if ( FO_do_deltam_scaling ) FO_L_TRUNCFAC(na,q) = L_greekmat(q,M21,na,1) * FO_TRUNCFAC(na)
         L_scale = FO_L_TRUNCFAC(na,q)*FO_omegas(na) - FO_TRUNCFAC(na)*FO_L_omegas(na,q)
         FO_L_DELTAUS(na,q)  = ( L_scale + scale * L_deltaus(q,na) ) * deltaus(na)
         diff = FO_heights(na-1) - FO_heights(na)
         FO_L_EXTINCS(na,q) = FO_L_DELTAUS(na,q) / diff
       enddo
     endif
   enddo

!  initial settings
!  ----------------

!  LOSPATHS variable should be false

   FO_do_LOSPATHS = .false.

!  other fixed flags

   FO_do_Thermset = .true.

!  pseudo-spherical flag

   FO_do_enhanced_PS     = .not.FO_do_regular_ps

!  Angles

   FO_alpha_boa   = alpha_boa 
   FO_Mu1         = cos(FO_alpha_boa * FO_dtr)

!  Attenuation control

   FO_nfinedivs = FO_nfineinput                ! Now an input
   FO_doCrit    = .not.FO_do_regular_ps
!   FO_Acrit     = 1.0e-10_ffp                 ! Now an input

!  Geometry call
!  -------------

   call DT_Geometry           &
          ( FO_maxlayers, FO_maxfine, FO_nlayers, FO_nfinedivs,                      & ! Input
            FO_dtr, FO_eradius, FO_heights, FO_alpha_boa,                            & ! Input
            FO_do_LOSpaths, FO_do_regular_ps, FO_doCrit, FO_Acrit, FO_extincs,       & ! Input
            FO_doNadir, FO_Raycon, FO_radii, FO_alpha, FO_cota,                      & ! Input/Output
            FO_xfine, FO_wfine, FO_csqfine, FO_cotfine, FO_alphafine, FO_radiifine,  & ! Input/Output
            FO_NCrit, FO_RadCrit, FO_CotCrit, FO_fail, FO_message, FO_trace )          ! Output

!  Exception handling on geometry

   if ( FO_Fail ) then
      write(*,*)FO_message ; write(*,*)FO_trace ; return
   endif

!  Do the source function integration. UPWELLING
!  =============================================

   if ( do_upwelling ) then
      call DTEV_Integral_Plus_1_UP &
           ( FO_maxlayers, FO_maxfine, FO_maxuserlevels, FO_max_atmoswfs, FO_max_surfacewfs,     & ! Inputs (dimensioning)
             FO_do_Thermset, FO_do_deltam_scaling, FO_do_regular_ps, FO_do_enhanced_ps,          & ! Inputs (Flags)
             FO_doNadir, FO_do_columnwfs, FO_do_profilewfs, FO_do_surfacewfs,                    & ! Inputs (flags)
             FO_nlayers, FO_nfinedivs, FO_n_user_levels, FO_user_levels,                         & ! Inputs (control)
             FO_n_columnwfs, FO_n_surfacewfs, FO_Lvaryflags, FO_Lvarynums,                       & ! Inputs (Jacobian control)
             FO_bb_input, FO_surfbb, FO_user_emissivity, FO_extincs,                             & ! Inputs (Optical)
             FO_deltaus, FO_omegas, FO_truncfac, FO_LS_user_emissivity,                          & ! Inputs (Optical)
             FO_L_extincs, FO_L_deltaus, FO_L_omegas, FO_L_truncfac,                             & ! Inputs (Optical)
             FO_Mu1, FO_NCrit, FO_Raycon, FO_cota, FO_xfine, FO_wfine, FO_csqfine, FO_cotfine,   & ! Inputs (Geometry)
             FO_stokes_dta, FO_stokes_dts, FO_tcom1, FO_L_tcom1,                                 & ! Output
             FO_LC_Jacobians_dta, FO_LP_Jacobians_dta,                                           & ! Output
             FO_LC_Jacobians_dts, FO_LP_Jacobians_dts, FO_LS_Jacobians_dts )                       ! Output
   endif

!  Do the source function integration. DNWELLING
!  =============================================

   if ( do_dnwelling ) then
      call DTEV_Integral_Plus_1_DN &
           ( FO_maxlayers, FO_maxfine, FO_maxuserlevels, FO_max_atmoswfs,                        & ! Inputs (dimensioning)
             FO_do_Thermset, FO_do_deltam_scaling, FO_do_regular_ps, FO_do_enhanced_ps,          & ! Inputs (Flags)
             FO_doNadir, FO_do_columnwfs, FO_do_profilewfs,                                      & ! Inputs (flags)
             FO_nlayers, FO_nfinedivs, FO_n_user_levels, FO_user_levels,                         & ! Inputs (control)
             FO_n_columnwfs,  FO_Lvaryflags, FO_Lvarynums,                                       & ! Inputs (Jacobian control)
             FO_bb_input, FO_extincs, FO_deltaus, FO_omegas, FO_truncfac,                        & ! Inputs (Optical)
             FO_L_extincs, FO_L_deltaus, FO_L_omegas, FO_L_truncfac,                             & ! Inputs (Optical)
             FO_Mu1, FO_NCrit, FO_Radcrit, FO_CotCrit, FO_Raycon, FO_cota,                       & ! Inputs (Geometry)
             FO_Radii, FO_xfine, FO_wfine, FO_csqfine, FO_cotfine,                               & ! Inputs (Geometry)
             FO_stokes_dta, FO_tcom1, FO_L_tcom1, FO_LC_Jacobians_dta, FO_LP_Jacobians_dta )       ! Output
   endif

!  Finish

   return
end subroutine FO_Vector_DT_Master_Plus

!  Finish module

end module FO_Vector_DT_Master_Plus_m

