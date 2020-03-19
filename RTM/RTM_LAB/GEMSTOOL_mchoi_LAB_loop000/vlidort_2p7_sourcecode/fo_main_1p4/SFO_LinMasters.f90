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
! ###########################################################

! ###########################################################
! #                                                         #
! #              FIRST-ORDER VECTOR MODEL                   #
! #     (EXACT SINGLE-SCATTERING and DIRECT-THERMAL)        #
! #                                                         #
! #  This Version :   1.4 F90                               #
! #  Release Date :   August 2013                           #
! #                                                         #
! #   Version 1.1,  13 February 2012, First Code            #
! #   Version 1.2,  01 June     2012, Modularization        #
! #   Version 1.3a, 29 October  2012, Obsgeom Multi-geom.   #
! #   Version 1.3b, 24 January  2013, BRDF/SL Supplements   #
! #   Version 1.4,  31 July     2013, Lattice Multi-geom.   #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

!  This module prepared for VLIDORT ingestion, R.Spurr 3/19/15

module SFO_LinMasters_m

!  All subroutines public

public  :: SFO_LPS_MASTER, &
           SFO_LCS_MASTER

contains

subroutine SFO_LPS_MASTER &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxfine,             & ! Input max dims
         maxmoments_input, max_user_levels,  max_atmoswfs, max_surfacewfs,    & ! Input max dims
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, nmoments_input,         & ! Input dims
         n_user_levels, n_surfacewfs,                                         & ! Input dims
         do_solar_sources, do_thermal_emission, do_surface_emission,          & ! Input flags
         do_planpar, do_regular_ps, do_enhanced_ps, do_deltam_scaling,        & ! Input flags
         do_upwelling, do_dnwelling, do_ObsGeom, do_profilewfs, do_surfacewfs,& ! Input flags
         dtr, Pie, doCrit, Acrit, eradius, heights, user_levels,              & ! Input general
         obsgeom_boa, theta_boa, alpha_boa, phi_boa,                          & ! Input geometry
         flux, extinction, deltaus, omega, phasmoms, truncfac,                & ! Input - std atmos optical
         bb_input, reflec, surfbb, emiss,                                     & ! Input - std surf optical
         Lvaryflags, Lvarynums, Lvarymoms,                                    & ! Input - lin control
         L_extinction, L_deltaus, L_omega, L_phasmoms, L_truncfac,            & ! Input - lin atmos optical
         LS_reflec, LS_emiss,                                                 & ! Input - lin surf optical
         fo_intensity_ss, fo_intensity_db, fo_intensity_dta, fo_intensity_dts,            & ! Output - std
         fo_intensity,                                                           & ! Output - std
         fo_profilewf_ss,  fo_profilewf_db,  fo_surfacewf_db,                 & ! Output - lin
         fo_profilewf_dta, fo_profilewf_dts, fo_surfacewf_dts,                & ! Output - lin
         fo_profilewf, fo_surfacewf,                                          & ! Output - lin
         Master_fail, message, trace_1, trace_2 )                               ! Output (Exception handling)

!  Use modules

   USE FO_SSGeometry_Master_m
   USE FO_DTGeometry_Master_m

   USE FO_ScalarSS_spherfuncs_m
   USE FO_ScalarSS_RTCalcs_ILPS_m

   USE FO_Thermal_RTCalcs_ILPS_m

   implicit none

!  parameter arguments

   integer, parameter :: ffp = selected_real_kind(15),&
                         max_directions = 2, upidx = 1, dnidx = 2

!  Subroutine inputs
!  =================

!  Max dimensions
!  --------------

   integer  :: maxgeoms, maxszas, maxvzas, maxazms
   integer  :: maxlayers, maxfine, maxmoments_input, max_user_levels
   integer  :: max_atmoswfs, max_surfacewfs

!  Dimensions
!  ----------

!  Layer and geometry control. Finelayer divisions may be changed

   integer, intent(in)    :: ngeoms, nszas, nvzas, nazms, nlayers, nfine
   integer, intent(in)    :: nmoments_input
   integer, intent(in)    :: n_user_levels

!  Number of surface weighting functions

   integer, intent(in)    :: n_surfacewfs

!  Configuration inputs
!  --------------------

!  Sources control, including thermal

   logical, intent(in)  :: do_solar_sources
   logical, intent(in)  :: do_thermal_emission
   logical, intent(in)  :: do_surface_emission

!  Flags (sphericity flags should be mutually exclusive)

   logical, intent(in)  :: do_planpar
   logical, intent(in)  :: do_regular_ps
   logical, intent(in)  :: do_enhanced_ps

!  deltam scaling flag

   logical, intent(in)  :: do_deltam_scaling

!  Directional Flags

   logical, intent(in)  :: do_upwelling, do_dnwelling

!  Obsgeom flag

   logical, intent(in)     :: do_Obsgeom

!  Linearization flags

   logical, intent(in)  :: do_surfacewfs, do_profilewfs

!  General inputs
!  --------------

!  DTR = degrees-to-Radians. Pie = 3.14159...

   real(ffp), intent(in) :: dtr, Pie

!  Critical adjustment for cloud layers

   real(ffp), intent(in)  :: Acrit
   logical, intent(inout)  :: doCrit

!  Earth radius + heights

   real(ffp), intent(in) :: eradius, heights (0:maxlayers)

!  Output levels

   integer, intent(in)   :: user_levels ( max_user_levels )

!  Geometry inputs
!  ---------------

!  input angles (Degrees). Enough information for Lattice or Obsgeom.
!   Convention for ObsGeom = same as VLIDORT/LIDORT (1=sza,2=vza,3=azm)
!    In both cases, the Phi angle may be changed.....

   real(ffp), intent(inout)  :: Obsgeom_boa(maxgeoms,3)
   real(ffp), intent(inout)  :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  Optical inputs
!  --------------

!  Solar flux

   real(ffp), intent(in) :: flux

!  Atmosphere

   real(ffp), intent(in) :: extinction  ( maxlayers )
   real(ffp), intent(in) :: deltaus     ( maxlayers )
   real(ffp), intent(in) :: omega       ( maxlayers )
   real(ffp), intent(in) :: phasmoms    ( maxlayers, 0:maxmoments_input )

!  For TMS correction

   real(ffp), intent(in) :: truncfac    ( maxlayers )

!  Thermal inputs

   real(ffp), intent(in) :: bb_input ( 0:maxlayers )

!  Surface properties - reflective (could be the albedo)

   real(ffp), intent(in) :: reflec(maxgeoms)

!  Surface properties - emissive

   real(ffp), intent(in) :: surfbb
   real(ffp), intent(in) :: emiss ( maxgeoms )

!  Linearized inputs
!  -----------------

!  Linearization control

   logical, intent(in)   :: Lvaryflags(maxlayers)
   integer, intent(in)   :: Lvarynums (maxlayers)
   logical, intent(in)   :: Lvarymoms (maxlayers,max_atmoswfs)

!  Linearized optical inputs

   real(ffp), intent(in) :: L_extinction  ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_deltaus     ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_omega       ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_phasmoms    ( maxlayers, 0:maxmoments_input, max_atmoswfs )

!  Linearized TMS correction

   real(ffp), intent(in) :: L_truncfac    ( maxlayers, max_atmoswfs )

!  Surface properties - reflective

   real(ffp), intent(in) :: LS_reflec ( maxgeoms, max_surfacewfs )

!  Surface properties - emissive

   real(ffp), intent(in) :: LS_emiss ( maxgeoms, max_surfacewfs )

!  Subroutine outputs
!  ==================

!  Solar

   real(ffp), intent(out) :: fo_intensity_ss ( max_user_levels,maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_intensity_db ( max_user_levels,maxgeoms )

!  Thermal

   real(ffp), intent(out) :: fo_intensity_dta ( max_user_levels,maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_intensity_dts ( max_user_levels,maxgeoms )

!  Composite

   real(ffp), intent(out) :: fo_intensity       ( max_user_levels,maxgeoms,max_directions )

!  Solar

   real(ffp), intent(out) :: fo_profilewf_ss ( max_atmoswfs,maxlayers,max_user_levels,&
                                               maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_profilewf_db ( max_atmoswfs,maxlayers,max_user_levels,&
                                               maxgeoms )
   real(ffp), intent(out) :: fo_surfacewf_db ( max_surfacewfs,max_user_levels,&
                                               maxgeoms )

!  Thermal

   real(ffp), intent(out) :: fo_profilewf_dta ( max_atmoswfs,maxlayers,max_user_levels,&
                                                maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_profilewf_dts ( max_atmoswfs,maxlayers,max_user_levels,&
                                                maxgeoms )
   real(ffp), intent(out) :: fo_surfacewf_dts ( max_surfacewfs,max_user_levels,&
                                                maxgeoms )

!  Composite

   real(ffp), intent(out) :: fo_profilewf       ( max_atmoswfs,maxlayers,max_user_levels,&
                                                  maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_surfacewf       ( max_surfacewfs,max_user_levels,&
                                                  maxgeoms )

!  Exception handling

   logical, intent(out)           :: Master_fail
   character (len=*), intent(out) :: message
   character (len=*), intent(out) :: trace_1, trace_2

!  Other variables
!  ===============

!  Geometry routine outputs
!  ------------------------

  !  VSIGN = +1 (Up); -1(Down)

   real(ffp)  :: vsign

!  LOSPATHS flag has been removed now.  8/1/13

!  Flag for the Nadir case

   logical    :: doNadir(maxgeoms)
  
!  Alphas,  Cotangents, Radii, Ray constant. 

   real(ffp)  :: radii    (0:maxlayers)
   real(ffp)  :: Raycon   (maxgeoms)
   real(ffp)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp)  :: cota     (0:maxlayers,maxgeoms)

!  Critical layer

   integer    :: Ncrit(maxgeoms)
   real(ffp)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  solar paths 

   integer    :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)

   integer    :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)
   real(ffp)  :: Chapfacs      (maxlayers,maxlayers,maxgeoms)

!  Cosine scattering angle, other cosines

   real(ffp)  :: cosscat (maxgeoms)
   real(ffp)  :: Mu0     (maxgeoms)
   real(ffp)  :: Mu1     (maxgeoms)

!  LOS Quadratures for Enhanced PS

   integer    :: nfinedivs(maxlayers,maxgeoms)
   real(ffp)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: wfine   (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: csqfine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: cotfine (maxlayers,maxfine,maxgeoms)

!  Fine layering output

   real(ffp)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: radiifine (maxlayers,maxfine,maxgeoms)

!  Spherfunc routine outputs
!  -------------------------

!  Help variables

   real(ffp) :: DF1(MAXMOMENTS_INPUT)
   real(ffp) :: DF2(MAXMOMENTS_INPUT)

!  Legendre Polynomials

   real(ffp)  :: LegPoly_up(0:maxmoments_input,maxgeoms)
   real(ffp)  :: LegPoly_dn(0:maxmoments_input,maxgeoms)

!  RT Calculation outputs
!  ----------------------

!  SS routines output

   real(ffp)  :: intensity_up    ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_dn    ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_db    ( max_user_levels,maxgeoms )

   real(ffp)  :: LP_Jacobians_up  ( max_user_levels,  maxgeoms, maxlayers, max_atmoswfs )
   real(ffp)  :: LP_Jacobians_dn  ( max_user_levels,  maxgeoms, maxlayers, max_atmoswfs )
   real(ffp)  :: LP_Jacobians_db  ( max_user_levels,  maxgeoms, maxlayers, max_atmoswfs )

   real(ffp)  :: LS_Jacobians_db  ( max_user_levels,  maxgeoms, max_surfacewfs )

!  Thermal routines output

!   real(ffp)  :: intensity_dta_up ( max_user_levels,maxgeoms )
!   real(ffp)  :: intensity_dta_dn ( max_user_levels,maxgeoms )
!   real(ffp)  :: intensity_dts    ( max_user_levels,maxgeoms )

!   real(ffp)  :: LP_Jacobians_dta_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
!   real(ffp)  :: LP_Jacobians_dta_dn  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
!   real(ffp)  :: LP_Jacobians_dts_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )

!   real(ffp)  :: LS_Jacobians_dts     ( max_user_levels, maxgeoms, max_surfacewfs )

!  Intermediate RT products
!  ------------------------

!  Composite

   real(ffp)  :: fo_intensity_atmos ( max_user_levels,maxgeoms,max_directions )
   real(ffp)  :: fo_intensity_surf  ( max_user_levels,maxgeoms )

   real(ffp)  :: fo_profilewf_atmos ( max_atmoswfs,maxlayers,max_user_levels,&
                                      maxgeoms,max_directions )
   real(ffp)  :: fo_profilewf_surf  ( max_atmoswfs,maxlayers,max_user_levels,&
                                      maxgeoms )
   real(ffp)  :: fo_surfacewf_surf  ( max_surfacewfs,max_user_levels,&
                                      maxgeoms )

!  LOS VARIABLES (THERMAL SOLUTION)
!  --------------------------------

   real(ffp)  :: intensity_dta_up_LOS ( max_user_levels,maxvzas )
   real(ffp)  :: intensity_dta_dn_LOS ( max_user_levels,maxvzas )
   real(ffp)  :: intensity_dts_LOS    ( max_user_levels,maxvzas )

   real(ffp)  :: LP_Jacobians_dta_up_LOS  ( max_user_levels, maxvzas, maxlayers, max_atmoswfs )
   real(ffp)  :: LP_Jacobians_dta_dn_LOS  ( max_user_levels, maxvzas, maxlayers, max_atmoswfs )
   real(ffp)  :: LP_Jacobians_dts_LOS     ( max_user_levels, maxvzas, maxlayers, max_atmoswfs )
   real(ffp)  :: LS_Jacobians_dts_LOS     ( max_user_levels, maxvzas, max_surfacewfs )

   logical    :: doNadir_LOS  (maxvzas)
   real(ffp)  :: Raycon_LOS   (maxvzas)
   real(ffp)  :: alpha_LOS    (0:maxlayers,maxvzas)
   real(ffp)  :: cota_LOS     (0:maxlayers,maxvzas)
   integer    :: Ncrit_LOS(maxvzas)
   real(ffp)  :: RadCrit_LOS(maxvzas), CotCrit_LOS(maxvzas)
   real(ffp)  :: Mu1_LOS(maxvzas)
   integer    :: nfinedivs_LOS(maxlayers,maxvzas)
   real(ffp)  :: xfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: wfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: csqfine_LOS (maxlayers,maxfine,maxvzas)
   real(ffp)  :: cotfine_LOS (maxlayers,maxfine,maxvzas)
   real(ffp)  :: alphafine_LOS (maxlayers,maxfine,maxvzas)
   real(ffp)  :: radiifine_LOS (maxlayers,maxfine,maxvzas)

!  Other products
!  --------------

!  Thermal setup and linearization

   real(ffp)  :: tcom1(maxlayers,2)
   real(ffp)  :: L_tcom1(maxlayers,2,max_atmoswfs)

!  Dummies

!   real(ffp)  :: SScumsource_up     ( 0:maxlayers,maxgeoms )
!   real(ffp)  :: SScumsource_dn     ( 0:maxlayers,maxgeoms )
!   real(ffp)  :: DTcumsource_up     ( 0:maxlayers,maxgeoms )
!   real(ffp)  :: DTcumsource_dn     ( 0:maxlayers,maxgeoms )

!  LOCAL HELP VARIABLES
!  --------------------

!  numbers

   real(ffp), parameter :: zero = 0.0_ffp, one = 1.0_ffp

!  help variables

   integer   :: ns, nv, na, g, n, spar, par, lev, nv_offset(maxszas), na_offset(maxszas,maxvzas)
   logical   :: STARTER, do_Thermset, fail, do_Chapman

!  Initialize output

   fo_intensity_ss  = zero
   fo_intensity_db  = zero
   fo_intensity_dta = zero
   fo_intensity_dts = zero
   fo_intensity     = zero

   fo_profilewf_ss  = zero
   fo_profilewf_db  = zero
   fo_surfacewf_db  = zero
   fo_profilewf_dta = zero
   fo_profilewf_dts = zero
   fo_surfacewf_dts = zero
   fo_profilewf     = zero
   fo_surfacewf     = zero

   Master_fail = .false.
   message = ' '
   trace_1 = ' '
   trace_2 = ' '

!  Flags to be set for each calculation (safety)

   do_Chapman  = .false.
   do_Thermset = .true.
   starter     = .true.

!  Offsets

   na_offset = 0 ; nv_offset = 0
   if ( .not. do_obsgeom ) then
     do ns = 1, nszas
       nv_offset(ns) = nvzas * nazms * (ns - 1) 
       do nv = 1, nvzas
         na_offset(ns,nv) = nv_offset(ns) + nazms * (nv - 1)
       enddo
     enddo
   endif

!  Solar sources run (NO THERMAL)
!  ------------------------------

   if ( do_solar_sources ) then

!  Upwelling

     if ( do_upwelling ) then
       vsign =  1.0_ffp

!  Geometry call

       call FO_SSGeometry_Master &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxfine,            & ! Input dimensions
         do_obsgeom, do_Chapman, do_planpar, do_enhanced_ps,                 & ! Input flags
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, dtr, Pie, vsign,       & ! Input control and constants
         eradius, heights, obsgeom_boa, alpha_boa, theta_boa, phi_boa,       & ! Input
         doNadir, doCrit, Acrit, extinction, Raycon, radii, alpha, cota,     & ! Input/Output(level)
         nfinedivs, xfine, wfine, csqfine, cotfine, alphafine, radiifine,    & ! Output(Fine)
         NCrit, RadCrit, CotCrit, Mu0, Mu1, cosscat, chapfacs,               & ! Output(Crit/scat)
         sunpaths, ntraverse, sunpathsfine, ntraversefine,                   & ! Output(Sunpaths)
         fail, message, trace_1 )                                                ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_SSGeometry_Master, Solar Sources, Upwelling calculation'
         Master_fail = .true. ; return
       endif

!  Spherical functions call

       Call FO_ScalarSS_spherfuncs ( STARTER, MAXMOMENTS_INPUT, MAXGEOMS, &
                 NMOMENTS_INPUT, NGEOMS, DF1, DF2, COSSCAT, LEGPOLY_UP )

!  RT Call Solar only

       call SS_Integral_ILPS_UP &
        ( maxgeoms, maxlayers, maxfine, maxmoments_input,                               & ! Inputs (dimensioning)
          max_user_levels, max_atmoswfs, max_surfacewfs,                                & ! Inputs (dimensioning)
          do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,        & ! Inputs (Flags - General)
          do_profilewfs, do_surfacewfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms, & ! Inputs (control, Jacobian )
          ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,       & ! Inputs (control,  output)
          reflec, extinction, deltaus, omega, truncfac, phasmoms, flux,                 & ! Inputs (Optical)
          LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,          & ! Inputs (Optical - Linearized)
          Mu0, Mu1, LegPoly_up, NCrit, xfine, wfine, csqfine, cotfine,                  & ! Inputs (Geometry)
          Raycon, cota, sunpaths, ntraverse, sunpathsfine, ntraversefine,               & ! Inputs (Geometry)
          Intensity_up, Intensity_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db )       ! Output

!  Save results

       do g = 1, ngeoms
          do lev=1,n_user_levels
             fo_intensity_ss(lev,g,upidx) = intensity_up(lev,g)
             fo_intensity_db(lev,g)       = intensity_db(lev,g)
             fo_intensity(lev,g,upidx)    = fo_intensity_ss(lev,g,upidx) &
                                          + fo_intensity_db(lev,g)
          enddo
       enddo

       if ( do_profilewfs ) then
          do g = 1, ngeoms
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_ss(par,n,lev,g,upidx) = &
                        LP_Jacobians_up(lev,g,n,par)
                     fo_profilewf_db(par,n,lev,g)       = &
                        LP_Jacobians_db(lev,g,n,par)
                     fo_profilewf(par,n,lev,g,upidx)    = &
                        fo_profilewf_ss(par,n,lev,g,upidx) &
                      + fo_profilewf_db(par,n,lev,g)
                   enddo
                 endif
               enddo
             enddo
          enddo
       endif

       if ( do_surfacewfs ) then
          do g = 1, ngeoms
             do lev=1,n_user_levels
               do spar=1,n_surfacewfs
                 fo_surfacewf_db(spar,lev,g) = &
                    LS_Jacobians_db(lev,g,spar)
                 fo_surfacewf(spar,lev,g)    = &
                    LS_Jacobians_db(lev,g,spar)
               enddo
             enddo
          enddo
       endif

! End upwelling

     end if

!  Donwelling

     if ( do_dnwelling ) then
       vsign =  -one

!  Geometry call

       call FO_SSGeometry_Master &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxfine,            & ! Input dimensions
         do_obsgeom, do_Chapman, do_planpar, do_enhanced_ps,                 & ! Input flags
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, dtr, Pie, vsign,       & ! Input control and constants
         eradius, heights, obsgeom_boa, alpha_boa, theta_boa, phi_boa,       & ! Input
         doNadir, doCrit, Acrit, extinction, Raycon, radii, alpha, cota,     & ! Input/Output(level)
         nfinedivs, xfine, wfine, csqfine, cotfine, alphafine, radiifine,    & ! Output(Fine)
         NCrit, RadCrit, CotCrit, Mu0, Mu1, cosscat, chapfacs,               & ! Output(Crit/scat)
         sunpaths, ntraverse, sunpathsfine, ntraversefine,                   & ! Output(Sunpaths)
         fail, message, trace_1 )                                                ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_SSGeometry_Master, Solar Sources, Downwelling calculation'
         Master_fail = .true. ; return
       endif

!  Spherical functions call

       Call FO_ScalarSS_spherfuncs ( STARTER, MAXMOMENTS_INPUT, MAXGEOMS, &
                 NMOMENTS_INPUT, NGEOMS, DF1, DF2, COSSCAT, LEGPOLY_DN )

! RT call - solar only

       call SS_Integral_ILPS_DN &
        ( maxgeoms, maxlayers, maxfine, maxmoments_input,                          & ! Inputs (dimensioning)
          max_user_levels, max_atmoswfs,                                           & ! Inputs (dimensioning)
          do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,   & ! Inputs (Flags)
          do_profilewfs, Lvaryflags, Lvarynums, Lvarymoms,                         & ! Inputs (control, Jacobian )
          ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,  & ! Inputs (control,  output)
          extinction, deltaus, omega, truncfac, phasmoms, flux,                    & ! Inputs (Optical - Regular)
          L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,                & ! Inputs (Optical - Linearized)
          Mu1, LegPoly_dn, radii, NCrit, RadCrit, CotCrit, xfine, wfine, csqfine,  & ! Inputs (Geometry)
          cotfine, Raycon, cota, sunpaths, ntraverse, sunpathsfine, ntraversefine, & ! Inputs (Geometry)
          Intensity_dn, LP_Jacobians_dn )

!  Save results

       do g = 1, ngeoms
          do lev=1,n_user_levels
             fo_intensity_ss(lev,g,dnidx) = intensity_dn(lev,g)
             fo_intensity(lev,g,dnidx)    = fo_intensity_ss(lev,g,dnidx)
          enddo
       enddo

       if ( do_profilewfs ) then
          do g = 1, ngeoms
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_ss(par,n,lev,g,dnidx) = LP_Jacobians_dn(lev,g,n,par)
                     fo_profilewf(par,n,lev,g,dnidx)    = fo_profilewf_ss(par,n,lev,g,dnidx)
                   enddo
                 endif
               enddo
             enddo
          enddo
       endif

!  End downwelling

     endif

!  End solar run

   endif

!  Thermal sources run
!  -------------------

   if ( do_thermal_emission.and.do_surface_emission ) then

!  DT Geometry call

     call FO_DTGeometry_Master  &
          ( maxvzas, maxlayers, maxfine, do_planpar, do_enhanced_ps,  & ! Input dimensions/flags
            nvzas, nlayers, nfine, dtr, eradius, heights, alpha_boa,  & ! Inputs
            doNadir_LOS, doCrit, Acrit, extinction,                   & ! Input/Output(level)
            Raycon_LOS, radii, alpha_LOS, cota_LOS,                   & ! Input/Output(level)
            nfinedivs_LOS, xfine_LOS, wfine_LOS, csqfine_LOS,         & ! Output(Fine)
            cotfine_LOS, alphafine_LOS, radiifine_LOS,                & ! Output(Fine)
            NCrit_LOS, RadCrit_LOS, CotCrit_LOS, Mu1_LOS,             & ! Output(Critical/Cos)
            fail, message, trace_1 )                                    ! Output(Status)

     if ( fail ) then
       trace_2 = 'Failure from FO_DTGeometry_Master'
       Master_fail = .true. ; return
     endif

!  Upwelling
!  ---------

     if ( do_upwelling ) then

!  Direct thermal, calculate

        call DTE_Integral_ILPS_UP &
         ( maxvzas, maxlayers, maxfine, max_user_levels,                      & ! Inputs (dimensioning)
           max_atmoswfs, max_surfacewfs,                                      & ! Inputs (dimensioning)
           Do_Thermset, do_deltam_scaling, do_PlanPar,                        & ! Inputs (Flags)
           do_regular_ps, do_enhanced_ps, doNadir_LOS,                        & ! Inputs (Flags)
           do_profilewfs, do_surfacewfs, Lvaryflags, Lvarynums, n_surfacewfs, & ! Inputs (Control, Jacobians)
           nvzas, nlayers, nfinedivs_LOS, n_user_levels, user_levels,         & ! Inputs (control, output)
           bb_input, surfbb, emiss, LS_emiss,                                 & ! Inputs (Thermal)
           extinction, deltaus, omega, truncfac,                              & ! Inputs (Optical - Regular)
           L_extinction, L_deltaus, L_omega, L_truncfac,                      & ! Inputs (Optical - Linearized)
           Mu1_LOS, NCrit_LOS, Raycon_LOS, cota_LOS,                          & ! Inputs (Geometry)
           xfine_LOS, wfine_LOS, csqfine_LOS, cotfine_LOS,                    & ! Inputs (Geometry)
           intensity_dta_up_LOS, intensity_dts_LOS, LP_Jacobians_dta_up_LOS,  & ! Outputs
           LP_Jacobians_dts_LOS, LS_Jacobians_dts_LOS, tcom1, L_tcom1  )        ! Output

!  Save results

       if ( do_obsgeom ) then
         do g = 1, nvzas
           do lev=1,n_user_levels
             fo_intensity_dta(lev,g,upidx) = intensity_dta_up_LOS(lev,g)
             fo_intensity_dts(lev,g)       = intensity_dts_LOS(lev,g)
             fo_intensity(lev,g,upidx)     = fo_intensity_dta(lev,g,upidx) + fo_intensity_dts(lev,g)
           enddo
         enddo
       else
         do nv = 1, nvzas
           do ns = 1, nszas
             do na = 1, nazms
               g = na_offset(ns,nv) + na
               do lev=1,n_user_levels
                 fo_intensity_dta(lev,g,upidx) = intensity_dta_up_LOS(lev,nv)
                 fo_intensity_dts(lev,g)       = intensity_dts_LOS(lev,nv)
                 fo_intensity(lev,g,upidx)     = fo_intensity_dta(lev,g,upidx) + fo_intensity_dts(lev,g)
               enddo
             enddo
           enddo
         enddo
       endif

       if ( do_profilewfs ) then
         if ( do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_dta(par,n,lev,g,upidx) = LP_Jacobians_dta_up_LOS(lev,g,n,par)
                     fo_profilewf_dts(par,n,lev,g) = LP_Jacobians_dts_LOS(lev,g,n,par)
                     fo_profilewf(par,n,lev,g,upidx)    = &
                        fo_profilewf_dta(par,n,lev,g,upidx) + fo_profilewf_dts(par,n,lev,g)
                   enddo
                 endif
               enddo
             enddo
           enddo
         else
           do nv = 1, nvzas
             do ns = 1, nszas
               do na = 1, nazms
                 g = na_offset(ns,nv) + na
                 do lev=1,n_user_levels
                   do n = 1, nlayers
                     if ( Lvaryflags(n) ) then
                       do par=1,Lvarynums(n)
                         fo_profilewf_dta(par,n,lev,g,upidx) = LP_Jacobians_dta_up_LOS(lev,nv,n,par)
                         fo_profilewf_dts(par,n,lev,g) = LP_Jacobians_dts_LOS(lev,nv,n,par)
                         fo_profilewf(par,n,lev,g,upidx)    = &
                            fo_profilewf_dta(par,n,lev,g,upidx) + fo_profilewf_dts(par,n,lev,g)
                       enddo
                     endif
                   enddo
                 enddo
               enddo
             enddo
           enddo
         endif
       endif

       if ( do_surfacewfs ) then
         if ( do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do spar=1,n_surfacewfs
                 fo_surfacewf_dts(spar,lev,g) = LS_Jacobians_dts_LOS(lev,g,spar)
                 fo_surfacewf(spar,lev,g)     = fo_surfacewf_dts(spar,lev,g)
               enddo
             enddo
           enddo
         else
           do nv = 1, nvzas
             do ns = 1, nszas
               do na = 1, nazms
                 g = na_offset(ns,nv) + na
                 do spar=1,n_surfacewfs
                   fo_surfacewf_dts(spar,lev,g) = LS_Jacobians_dts_LOS(lev,nv,spar)
                   fo_surfacewf(spar,lev,g)     = fo_surfacewf_dts(spar,lev,g)
                 enddo
               enddo
             enddo
           enddo
         endif
       endif

!  End upwelling

     endif

!  Downwelling
!  -----------

     if ( do_dnwelling ) then

!  Direct thermal, calculate

!mick fix 4/3/2015 - swapped radii & cota_LOS in argument list
        call DTE_Integral_ILPS_DN &
         ( maxvzas, maxlayers, maxfine, max_user_levels, max_atmoswfs,         & ! Inputs (dimensioning)
           Do_Thermset, do_deltam_scaling, do_PlanPar,                         & ! Inputs (Flags)
           do_regular_ps, do_enhanced_ps, doNadir_LOS,                         & ! Inputs (Flags)
           do_profilewfs, Lvaryflags, Lvarynums,                               & ! Inputs (Control, Jacobians)
           nvzas, nlayers, nfinedivs_LOS, n_user_levels, user_levels,          & ! Inputs (control, output)
           bb_input, extinction, deltaus, omega, truncfac,                     & ! Inputs (Optical - Regular)
           L_extinction, L_deltaus, L_omega, L_truncfac,                       & ! Inputs (Optical - Linearized)
           Mu1_LOS, NCrit_LOS, RadCrit_LOS, CotCrit_LOS, Raycon_LOS, cota_LOS, & ! Inputs (Geometry)
           radii, xfine_LOS, wfine_LOS, csqfine_LOS, cotfine_LOS,              & ! Inputs (Geometry) 
           intensity_dta_dn_LOS, LP_Jacobians_dta_dn_LOS, tcom1, L_tcom1 )       ! Output

!  Save results

       if ( do_obsgeom ) then
         do g = 1, nvzas
           do lev=1,n_user_levels
             fo_intensity_dta(lev,g,dnidx) = intensity_dta_dn_LOS(lev,g)
             fo_intensity(lev,g,dnidx)     = fo_intensity_dta(lev,g,dnidx)
           enddo
         enddo
       else
         do nv = 1, nvzas
           do ns = 1, nszas
             do na = 1, nazms
               g = na_offset(ns,nv) + na
               do lev=1,n_user_levels
                 fo_intensity_dta(lev,g,dnidx) = intensity_dta_dn_LOS(lev,nv)
                 fo_intensity(lev,g,dnidx)     = fo_intensity_dta(lev,g,dnidx)
               enddo
             enddo
           enddo
         enddo
       endif

       if ( do_profilewfs ) then
          if (do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do n = 1, nlayers
                 if ( Lvaryflags(n) ) then
                   do par=1,Lvarynums(n)
                     fo_profilewf_dta(par,n,lev,g,dnidx) = LP_Jacobians_dta_dn_LOS(lev,g,n,par)
                     fo_profilewf(par,n,lev,g,dnidx)     = fo_profilewf_dta(par,n,lev,g,dnidx)
                   enddo
                 endif
               enddo
             enddo
           enddo
         else
           do nv = 1, nvzas
             do ns = 1, nszas
               do na = 1, nazms
                 g = na_offset(ns,nv) + na
                 do lev=1,n_user_levels
                   do n = 1, nlayers
                     if ( Lvaryflags(n) ) then
                       do par=1,Lvarynums(n)
                         fo_profilewf_dta(par,n,lev,g,dnidx) = LP_Jacobians_dta_dn_LOS(lev,nv,n,par)
                         fo_profilewf(par,n,lev,g,dnidx)     = fo_profilewf_dta(par,n,lev,g,dnidx)
                       enddo
                     endif
                   enddo
                 enddo
               enddo
             enddo
           enddo
         endif
       endif

!  end downwelling

      endif

!  End Thermal run

   endif

!  Solar and Thermal sources run
!  -----------------------------

! Add solar and thermal components

   if ( do_solar_sources.and.(do_thermal_emission.and.do_surface_emission)  ) then

!  Upwelling

     if ( do_upwelling ) then

       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_intensity_atmos(lev,g,upidx) = fo_intensity_ss(lev,g,upidx) + fo_intensity_dta(lev,g,upidx)
           fo_intensity_surf(lev,g)        = fo_intensity_db(lev,g)       + fo_intensity_dts(lev,g)
           fo_intensity(lev,g,upidx)       = fo_intensity_atmos(lev,g,upidx) + fo_intensity_surf(lev,g)
         enddo
       enddo

       if ( do_profilewfs ) then
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do n = 1, nlayers
               if ( Lvaryflags(n) ) then
                 do par=1,Lvarynums(n)
                   fo_profilewf_atmos(par,n,lev,g,upidx) = &
                      fo_profilewf_ss(par,n,lev,g,upidx) + fo_profilewf_dta(par,n,lev,g,upidx)
                   fo_profilewf_surf(par,n,lev,g)   = &
                      fo_profilewf_db(par,n,lev,g) + fo_profilewf_dts(par,n,lev,g)
                   fo_profilewf(par,n,lev,g,upidx)  = &
                      fo_profilewf_atmos(par,n,lev,g,upidx) + fo_profilewf_surf(par,n,lev,g)
                 enddo
               endif
             enddo
           enddo
         enddo
       endif

       if ( do_surfacewfs ) then
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do spar=1,n_surfacewfs
               fo_surfacewf_surf(spar,lev,g) = &
                  fo_surfacewf_db(spar,lev,g) + fo_surfacewf_dts(spar,lev,g)
               fo_surfacewf(spar,lev,g) = fo_surfacewf_surf(spar,lev,g)
             enddo
           enddo
         enddo
       endif

     endif

!   Downwelling

     if ( do_dnwelling ) then

       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_intensity_atmos(lev,g,dnidx) = fo_intensity_ss(lev,g,dnidx) + fo_intensity_dta(lev,g,dnidx)
           fo_intensity(lev,g,dnidx)       = fo_intensity_atmos(lev,g,dnidx)
         enddo
       enddo

       if ( do_profilewfs ) then
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do n = 1, nlayers
               if ( Lvaryflags(n) ) then
                 do par=1,Lvarynums(n)
                   fo_profilewf_atmos(par,n,lev,g,dnidx) = &
                      fo_profilewf_ss(par,n,lev,g,dnidx) + fo_profilewf_dta(par,n,lev,g,dnidx)
                   fo_profilewf(par,n,lev,g,dnidx)    = fo_profilewf_atmos(par,n,lev,g,dnidx)
                 enddo
               endif
             enddo
           enddo
         enddo
       endif

     endif

!  End solar+thermal

   endif

!  Finish

   return
end subroutine SFO_LPS_MASTER

!

subroutine SFO_LCS_MASTER &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxfine,            & ! Input max dims
         maxmoments_input, max_user_levels,  max_atmoswfs, max_surfacewfs,   & ! Input max dims
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, nmoments_input,        & ! Input dims
         n_user_levels, n_columnwfs, n_surfacewfs,                           & ! Input dims
         do_solar_sources, do_thermal_emission, do_surface_emission,         & ! Input flags
         do_planpar, do_regular_ps, do_enhanced_ps, do_deltam_scaling,       & ! Input flags
         do_upwelling, do_dnwelling, do_ObsGeom, do_columnwfs, do_surfacewfs,& ! Input flags
         dtr, Pie, doCrit, Acrit, eradius, heights, user_levels,             & ! Input general
         obsgeom_boa, theta_boa, alpha_boa, phi_boa,                         & ! Input geometry
         flux, extinction, deltaus, omega, phasmoms, truncfac,               & ! Input - std atmos optical
         bb_input, reflec, surfbb, emiss,                                    & ! Input - std surf optical
         Lvarymoms,                                                          & ! Input - lin control
         L_extinction, L_deltaus, L_omega, L_phasmoms, L_truncfac,           & ! Input - lin atmos optical
         LS_reflec, LS_emiss,                                                & ! Input - lin surf optical
         fo_intensity_ss, fo_intensity_db, fo_intensity_dta, fo_intensity_dts,           & ! Output - std
         fo_intensity,                                                          & ! Output - std
         fo_columnwf_ss,  fo_columnwf_db,  fo_surfacewf_db,                  & ! Output - lin
         fo_columnwf_dta, fo_columnwf_dts, fo_surfacewf_dts,                 & ! Output - lin
         fo_columnwf, fo_surfacewf,                                          & ! Output - lin
         Master_fail, message, trace_1, trace_2 )                              ! Output (Exception handling)

!  Use modules

   USE FO_SSGeometry_Master_m
   USE FO_DTGeometry_Master_m

   USE FO_ScalarSS_spherfuncs_m
   USE FO_ScalarSS_RTCalcs_ILCS_m

   USE FO_Thermal_RTCalcs_ILCS_m

   implicit none

!  parameter arguments

   integer, parameter :: ffp = selected_real_kind(15),&
                         max_directions = 2, upidx = 1, dnidx = 2

!  Subroutine inputs
!  =================

!  Max dimensions
!  --------------

   integer  :: maxgeoms, maxszas, maxvzas, maxazms
   integer  :: maxlayers, maxfine, maxmoments_input, max_user_levels
   integer  :: max_atmoswfs, max_surfacewfs

!  Dimensions
!  ----------

!  Layer and geometry control. Finelayer divisions may be changed

   integer, intent(in)    :: ngeoms, nszas, nvzas, nazms, nlayers, nfine
   integer, intent(in)    :: nmoments_input
   integer, intent(in)    :: n_user_levels

!  Number of column and surface weighting functions

   integer, intent(in)    :: n_columnwfs, n_surfacewfs

!  Configuration inputs
!  --------------------

!  Sources control, including thermal

   logical, intent(in)  :: do_solar_sources
   logical, intent(in)  :: do_thermal_emission
   logical, intent(in)  :: do_surface_emission

!  Flags (sphericity flags should be mutually exclusive)

   logical, intent(in)  :: do_planpar
   logical, intent(in)  :: do_regular_ps
   logical, intent(in)  :: do_enhanced_ps

!  deltam scaling flag

   logical, intent(in)  :: do_deltam_scaling

!  Directional Flags

   logical, intent(in)  :: do_upwelling, do_dnwelling

!  Obsgeom flag

   logical, intent(in)  :: do_Obsgeom

!  Linearization flags

   logical, intent(in)  :: do_surfacewfs, do_columnwfs

!  General inputs
!  --------------

!  DTR = degrees-to-Radians. Pie = 3.14159...

   real(ffp), intent(in) :: dtr, Pie

!  Critical adjustment for cloud layers

   real(ffp), intent(in)   :: Acrit
   logical, intent(inout)  :: doCrit

!  Earth radius + heights

   real(ffp), intent(in) :: eradius, heights (0:maxlayers)

!  Output levels

   integer, intent(in)   :: user_levels ( max_user_levels )

!  Geometry inputs
!  ---------------

!  input angles (Degrees). Enough information for Lattice or Obsgeom.
!   Convention for ObsGeom = same as VLIDORT/LIDORT (1=sza,2=vza,3=azm)
!    In both cases, the Phi angle may be changed.....

   real(ffp), intent(inout)  :: Obsgeom_boa(maxgeoms,3)
   real(ffp), intent(inout)  :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  Optical inputs
!  --------------

!  Solar flux

   real(ffp), intent(in) :: flux

!  Atmosphere

   real(ffp), intent(in) :: extinction  ( maxlayers )
   real(ffp), intent(in) :: deltaus     ( maxlayers )
   real(ffp), intent(in) :: omega       ( maxlayers )
   real(ffp), intent(in) :: phasmoms    ( maxlayers, 0:maxmoments_input )

!  For TMS correction

   real(ffp), intent(in) :: truncfac    ( maxlayers )

!  Thermal inputs

   real(ffp), intent(in) :: bb_input ( 0:maxlayers )

!  Surface properties - reflective (could be the albedo)

   real(ffp), intent(in) :: reflec(maxgeoms)

!  Surface properties - emissive

   real(ffp), intent(in) :: surfbb
   real(ffp), intent(in) :: emiss ( maxgeoms )

!  Linearized inputs
!  -----------------

!  Linearization control

   logical, intent(in)   :: Lvarymoms (maxlayers,max_atmoswfs)

!  Linearized optical inputs

   real(ffp), intent(in) :: L_extinction  ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_deltaus     ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_omega       ( maxlayers, max_atmoswfs )
   real(ffp), intent(in) :: L_phasmoms    ( maxlayers, 0:maxmoments_input, max_atmoswfs )

!  Linearized TMS correction

   real(ffp), intent(in) :: L_truncfac    ( maxlayers, max_atmoswfs )

!  Surface properties - reflective

   real(ffp), intent(in) :: LS_reflec (   &
                                        maxgeoms, max_surfacewfs )

!  Surface properties - emissive

   real(ffp), intent(in) :: LS_emiss ( maxgeoms, max_surfacewfs )

!  Subroutine outputs
!  ==================

!  Solar

   real(ffp), intent(out) :: fo_intensity_ss ( max_user_levels,maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_intensity_db ( max_user_levels,maxgeoms )

!  Thermal

   real(ffp), intent(out) :: fo_intensity_dta ( max_user_levels,maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_intensity_dts ( max_user_levels,maxgeoms )

!  Composite

   real(ffp), intent(out) :: fo_intensity       ( max_user_levels,maxgeoms,max_directions )

!  Solar

   real(ffp), intent(out) :: fo_columnwf_ss  ( max_atmoswfs,max_user_levels,&
                                               maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_columnwf_db  ( max_atmoswfs,max_user_levels,&
                                               maxgeoms )
   real(ffp), intent(out) :: fo_surfacewf_db ( max_surfacewfs,max_user_levels,&
                                               maxgeoms )

!  Thermal

   real(ffp), intent(out) :: fo_columnwf_dta  ( max_atmoswfs,max_user_levels,&
                                                maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_columnwf_dts  ( max_atmoswfs,max_user_levels,&
                                                maxgeoms )
   real(ffp), intent(out) :: fo_surfacewf_dts ( max_surfacewfs,max_user_levels,&
                                                maxgeoms )

!  Composite

   real(ffp), intent(out) :: fo_columnwf        ( max_atmoswfs,max_user_levels,&
                                                  maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_surfacewf       ( max_surfacewfs,max_user_levels,&
                                                  maxgeoms )

!  Exception handling

   logical, intent(out)           :: Master_fail
   character (len=*), intent(out) :: message
   character (len=*), intent(out) :: trace_1, trace_2

!  Other variables
!  ===============

!  Geometry routine outputs
!  ------------------------

  !  VSIGN = +1 (Up); -1(Down)

   real(ffp)  :: vsign

!  LOSPATHS flag has been removed now.  8/1/13

!  Flag for the Nadir case

   logical    :: doNadir(maxgeoms)
  
!  Alphas,  Cotangents, Radii, Ray constant. 

   real(ffp)  :: radii    (0:maxlayers)
   real(ffp)  :: Raycon   (maxgeoms)
   real(ffp)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp)  :: cota     (0:maxlayers,maxgeoms)

!  Critical layer

   integer    :: Ncrit(maxgeoms)
   real(ffp)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  solar paths 

   integer    :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)

   integer    :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)
   real(ffp)  :: Chapfacs      (maxlayers,maxlayers,maxgeoms)

!  Cosine scattering angle, other cosines

   real(ffp)  :: cosscat (maxgeoms)
   real(ffp)  :: Mu0     (maxgeoms)
   real(ffp)  :: Mu1     (maxgeoms)

!  LOS Quadratures for Enhanced PS

   integer    :: nfinedivs(maxlayers,maxgeoms)
   real(ffp)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: wfine   (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: csqfine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: cotfine (maxlayers,maxfine,maxgeoms)

!  Fine layering output

   real(ffp)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: radiifine (maxlayers,maxfine,maxgeoms)

!  Spherfunc routine outputs
!  -------------------------

!  Help variables

   real(ffp) :: DF1(MAXMOMENTS_INPUT)
   real(ffp) :: DF2(MAXMOMENTS_INPUT)

!  Legendre Polynomials

   real(ffp)  :: LegPoly_up(0:maxmoments_input,maxgeoms)
   real(ffp)  :: LegPoly_dn(0:maxmoments_input,maxgeoms)

!  RT Calculation outputs
!  ----------------------

!  SS routines output

   real(ffp)  :: intensity_up    ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_dn    ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_db    ( max_user_levels,maxgeoms )

   real(ffp)  :: LC_Jacobians_up  ( max_user_levels,  maxgeoms, max_atmoswfs )
   real(ffp)  :: LC_Jacobians_dn  ( max_user_levels,  maxgeoms, max_atmoswfs )
   real(ffp)  :: LC_Jacobians_db  ( max_user_levels,  maxgeoms, max_atmoswfs )

   real(ffp)  :: LS_Jacobians_db  ( max_user_levels,  maxgeoms, max_surfacewfs )

!  Thermal routines output

!   real(ffp)  :: intensity_dta_up ( max_user_levels,maxgeoms )
!   real(ffp)  :: intensity_dta_dn ( max_user_levels,maxgeoms )
!   real(ffp)  :: intensity_dts    ( max_user_levels,maxgeoms )

!   real(ffp)  :: LC_Jacobians_dta_up  ( max_user_levels, maxgeoms, max_atmoswfs )
!   real(ffp)  :: LC_Jacobians_dta_dn  ( max_user_levels, maxgeoms, max_atmoswfs )
!   real(ffp)  :: LC_Jacobians_dts_up  ( max_user_levels, maxgeoms, max_atmoswfs )

!   real(ffp)  :: LS_Jacobians_dts     ( max_user_levels, maxgeoms, max_surfacewfs )

!  Intermediate RT products
!  ------------------------

!  Composite

   real(ffp)  :: fo_intensity_atmos ( max_user_levels,maxgeoms,max_directions )
   real(ffp)  :: fo_intensity_surf  ( max_user_levels,maxgeoms )

   real(ffp)  :: fo_columnwf_atmos ( max_atmoswfs,max_user_levels,&
                                      maxgeoms,max_directions )
   real(ffp)  :: fo_columnwf_surf  ( max_atmoswfs,max_user_levels,&
                                      maxgeoms )
   real(ffp)  :: fo_surfacewf_surf  ( max_surfacewfs,max_user_levels,&
                                      maxgeoms )

!  LOS VARIABLES (THERMAL SOLUTION)
!  --------------------------------

   real(ffp)  :: intensity_dta_up_LOS ( max_user_levels,maxvzas )
   real(ffp)  :: intensity_dta_dn_LOS ( max_user_levels,maxvzas )
   real(ffp)  :: intensity_dts_LOS    ( max_user_levels,maxvzas )

   real(ffp)  :: LC_Jacobians_dta_up_LOS  ( max_user_levels, maxvzas, max_atmoswfs )
   real(ffp)  :: LC_Jacobians_dta_dn_LOS  ( max_user_levels, maxvzas, max_atmoswfs )
   real(ffp)  :: LC_Jacobians_dts_LOS     ( max_user_levels, maxvzas, max_atmoswfs )
   real(ffp)  :: LS_Jacobians_dts_LOS     ( max_user_levels, maxvzas, max_surfacewfs )

   logical    :: doNadir_LOS  (maxvzas)
   real(ffp)  :: Raycon_LOS   (maxvzas)
   real(ffp)  :: alpha_LOS    (0:maxlayers,maxvzas)
   real(ffp)  :: cota_LOS     (0:maxlayers,maxvzas)
   integer    :: Ncrit_LOS(maxvzas)
   real(ffp)  :: RadCrit_LOS(maxvzas), CotCrit_LOS(maxvzas)
   real(ffp)  :: Mu1_LOS(maxvzas)
   integer    :: nfinedivs_LOS(maxlayers,maxvzas)
   real(ffp)  :: xfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: wfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: csqfine_LOS (maxlayers,maxfine,maxvzas)
   real(ffp)  :: cotfine_LOS (maxlayers,maxfine,maxvzas)
   real(ffp)  :: alphafine_LOS (maxlayers,maxfine,maxvzas)
   real(ffp)  :: radiifine_LOS (maxlayers,maxfine,maxvzas)

!  Other products
!  --------------

!  Thermal setup and linearization

   real(ffp)  :: tcom1(maxlayers,2)
   real(ffp)  :: L_tcom1(maxlayers,2,max_atmoswfs)

!  Dummies

!   real(ffp)  :: SScumsource_up     ( 0:maxlayers,maxgeoms )
!   real(ffp)  :: SScumsource_dn     ( 0:maxlayers,maxgeoms )
!   real(ffp)  :: DTcumsource_up     ( 0:maxlayers,maxgeoms )
!   real(ffp)  :: DTcumsource_dn     ( 0:maxlayers,maxgeoms )

!  LOCAL HELP VARIABLES
!  --------------------

!  numbers

   real(ffp), parameter :: zero = 0.0_ffp, one = 1.0_ffp

!  help variables

   integer   :: ns, nv, na, g, par, spar, lev, nv_offset(maxszas), na_offset(maxszas,maxvzas)
   logical   :: STARTER, do_Thermset, fail, do_Chapman

!  Initialize output

   fo_intensity_ss  = zero
   fo_intensity_db  = zero
   fo_intensity_dta = zero
   fo_intensity_dts = zero
   fo_intensity     = zero

   fo_columnwf_ss   = zero
   fo_columnwf_db   = zero
   fo_surfacewf_db  = zero
   fo_columnwf_dta  = zero
   fo_columnwf_dts  = zero
   fo_surfacewf_dts = zero
   fo_columnwf      = zero
   fo_surfacewf     = zero

   Master_fail = .false.
   message = ' '
   trace_1 = ' '
   trace_2 = ' '

!  Flags to be set for each calculation (safety)

   do_Chapman  = .false.
   do_Thermset = .true.
   starter     = .true.

!  Offsets

   na_offset = 0 ; nv_offset = 0
   if ( .not. do_obsgeom ) then
     do ns = 1, nszas
       nv_offset(ns) = nvzas * nazms * (ns - 1) 
       do nv = 1, nvzas
         na_offset(ns,nv) = nv_offset(ns) + nazms * (nv - 1)
       enddo
     enddo
   endif

!  Solar sources run (NO THERMAL)
!  ------------------------------

   if ( do_solar_sources ) then

!  Upwelling

     if ( do_upwelling ) then
       vsign =  1.0_ffp

!  Geometry call

       call FO_SSGeometry_Master &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxfine,            & ! Input dimensions
         do_obsgeom, do_Chapman, do_planpar, do_enhanced_ps,                 & ! Input flags
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, dtr, Pie, vsign,       & ! Input control and constants
         eradius, heights, obsgeom_boa, alpha_boa, theta_boa, phi_boa,       & ! Input
         doNadir, doCrit, Acrit, extinction, Raycon, radii, alpha, cota,     & ! Input/Output(level)
         nfinedivs, xfine, wfine, csqfine, cotfine, alphafine, radiifine,    & ! Output(Fine)
         NCrit, RadCrit, CotCrit, Mu0, Mu1, cosscat, chapfacs,               & ! Output(Crit/scat)
         sunpaths, ntraverse, sunpathsfine, ntraversefine,                   & ! Output(Sunpaths)
         fail, message, trace_1 )                                                ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_SSGeometry_Master, Solar Sources, Upwelling calculation'
         Master_fail = .true. ; return
       endif

!  Spherical functions call

       Call FO_ScalarSS_spherfuncs ( STARTER, MAXMOMENTS_INPUT, MAXGEOMS, &
                 NMOMENTS_INPUT, NGEOMS, DF1, DF2, COSSCAT, LEGPOLY_UP )

!  RT Call Solar only

       call SS_Integral_ILCS_UP &
   ( maxgeoms, maxlayers, maxfine, maxmoments_input,                            & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs,                             & ! Inputs (dimensioning)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,     & ! Inputs (Flags - General)
     do_columnwfs, do_surfacewfs, n_columnwfs, n_surfacewfs, Lvarymoms,         & ! Inputs (control, Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,    & ! Inputs (control,  output)
     reflec, extinction, deltaus, omega, truncfac, phasmoms, flux,              & ! Inputs (Optical)
     LS_reflec, L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,       & ! Inputs (Optical - Linearized)
     Mu0, Mu1, LegPoly_up, NCrit, xfine, wfine, csqfine, cotfine,               & ! Inputs (Geometry)
     Raycon, cota, sunpaths, ntraverse, sunpathsfine, ntraversefine,            & ! Inputs (Geometry)
     Intensity_up, Intensity_db, LC_Jacobians_up, LC_Jacobians_db, LS_Jacobians_db )    ! Output

!  Save results

       do g = 1, ngeoms
          do lev=1,n_user_levels
             fo_intensity_ss(lev,g,upidx) = intensity_up(lev,g)
             fo_intensity_db(lev,g)       = intensity_db(lev,g)
             fo_intensity(lev,g,upidx)    = fo_intensity_ss(lev,g,upidx) &
                                          + fo_intensity_db(lev,g)
          enddo
       enddo

       if ( do_columnwfs ) then
          do g = 1, ngeoms
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_ss(par,lev,g,upidx) = LC_Jacobians_up(lev,g,par)
                 fo_columnwf_db(par,lev,g)       = LC_Jacobians_db(lev,g,par)
                 fo_columnwf(par,lev,g,upidx)    = &
                    fo_columnwf_ss(par,lev,g,upidx) + fo_columnwf_db(par,lev,g)
               enddo
             enddo
          enddo
       endif

       if ( do_surfacewfs ) then
          do g = 1, ngeoms
             do lev=1,n_user_levels
               do spar=1,n_surfacewfs
                 fo_surfacewf_db(spar,lev,g) = LS_Jacobians_db(lev,g,spar)
                 fo_surfacewf(spar,lev,g)    = LS_Jacobians_db(lev,g,spar)
               enddo
             enddo
          enddo
       endif

! End upwelling

     end if

!  Donwelling

     if ( do_dnwelling ) then
       vsign =  -one

!  Geometry call

       call FO_SSGeometry_Master &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxfine,            & ! Input dimensions
         do_obsgeom, do_Chapman, do_planpar, do_enhanced_ps,                 & ! Input flags
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, dtr, Pie, vsign,       & ! Input control and constants
         eradius, heights, obsgeom_boa, alpha_boa, theta_boa, phi_boa,       & ! Input
         doNadir, doCrit, Acrit, extinction, Raycon, radii, alpha, cota,     & ! Input/Output(level)
         nfinedivs, xfine, wfine, csqfine, cotfine, alphafine, radiifine,    & ! Output(Fine)
         NCrit, RadCrit, CotCrit, Mu0, Mu1, cosscat, chapfacs,               & ! Output(Crit/scat)
         sunpaths, ntraverse, sunpathsfine, ntraversefine,                   & ! Output(Sunpaths)
         fail, message, trace_1 )                                                ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_SSGeometry_Master, Solar Sources, Downwelling calculation'
         Master_fail = .true. ; return
       endif

!  Spherical functions call

       Call FO_ScalarSS_spherfuncs ( STARTER, MAXMOMENTS_INPUT, MAXGEOMS, &
                 NMOMENTS_INPUT, NGEOMS, DF1, DF2, COSSCAT, LEGPOLY_DN )

! RT call - solar only

          call SS_Integral_ILCS_DN &
   ( maxgeoms, maxlayers, maxfine, maxmoments_input,                            & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs,                                             & ! Inputs (dimensioning)
     do_deltam_scaling, do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,     & ! Inputs (Flags - General)
     do_columnwfs, n_columnwfs, Lvarymoms,                                      & ! Inputs (control, Jacobian )
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,    & ! Inputs (control,  output)
     extinction, deltaus, omega, truncfac, phasmoms, flux,                      & ! Inputs (Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms,                  & ! Inputs (Optical - Linearized)
     Mu1, LegPoly_dn, radii, NCrit, RadCrit, CotCrit, xfine, wfine, csqfine,    & ! Inputs (Geometry)
     cotfine, Raycon, cota, sunpaths, ntraverse, sunpathsfine, ntraversefine,   & ! Inputs (Geometry)
     Intensity_dn, LC_Jacobians_dn )                                              ! Output

!  Save results

       do g = 1, ngeoms
          do lev=1,n_user_levels
             fo_intensity_ss(lev,g,dnidx) = intensity_dn(lev,g)
             fo_intensity(lev,g,dnidx)    = fo_intensity_ss(lev,g,dnidx)
          enddo
       enddo

       if ( do_columnwfs ) then
          do g = 1, ngeoms
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_ss(par,lev,g,dnidx) = LC_Jacobians_dn(lev,g,par)
                 fo_columnwf(par,lev,g,dnidx)    = fo_columnwf_ss(par,lev,g,dnidx)
               enddo
             enddo
          enddo
       endif

!  End downwelling

     endif

!  End solar run

   endif

!  Thermal sources run
!  -------------------

   if ( do_thermal_emission.and.do_surface_emission ) then

!  DT Geometry call

     call FO_DTGeometry_Master  &
          ( maxvzas, maxlayers, maxfine, do_planpar, do_enhanced_ps,  & ! Input dimensions/flags
            nvzas, nlayers, nfine, dtr, eradius, heights, alpha_boa,  & ! Inputs
            doNadir_LOS, doCrit, Acrit, extinction,                   & ! Input/Output(level)
            Raycon_LOS, radii, alpha_LOS, cota_LOS,                   & ! Input/Output(level)
            nfinedivs_LOS, xfine_LOS, wfine_LOS, csqfine_LOS,         & ! Output(Fine)
            cotfine_LOS, alphafine_LOS, radiifine_LOS,                & ! Output(Fine)
            NCrit_LOS, RadCrit_LOS, CotCrit_LOS, Mu1_LOS,             & ! Output(Critical/Cos)
            fail, message, trace_1 )                                    ! Output(Status)

     if ( fail ) then
       trace_2 = 'Failure from FO_DTGeometry_Master'
       Master_fail = .true. ; return
     endif

!  Upwelling
!  ---------

     if ( do_upwelling ) then

!  Direct thermal, calculate

        call DTE_Integral_ILCS_UP &
         ( maxvzas, maxlayers, maxfine, max_user_levels,                    & ! Inputs (dimensioning)
           max_atmoswfs, max_surfacewfs,                                    & ! Inputs (dimensioning)
           Do_Thermset, do_deltam_scaling, do_PlanPar,                      & ! Inputs (Flags)
           do_regular_ps, do_enhanced_ps, doNadir_LOS,                      & ! Inputs (Flags)
           do_columnwfs, do_surfacewfs, n_columnwfs, n_surfacewfs,          & ! Inputs (Control, Jacobians)
           nvzas, nlayers, nfinedivs_LOS, n_user_levels, user_levels,       & ! Inputs (control, output)
           bb_input, surfbb, emiss, LS_emiss,                               & ! Inputs (Thermal)           
           extinction, deltaus, omega, truncfac,                            & ! Inputs (Optical - Regular)
           L_extinction, L_deltaus, L_omega, L_truncfac,                    & ! Inputs (Optical - Linearized)
           Mu1_LOS, NCrit_LOS, Raycon_LOS, cota_LOS,                        & ! Inputs (Geometry)
           xfine_LOS, wfine_LOS, csqfine_LOS, cotfine_LOS,                  & ! Inputs (Geometry)
           intensity_dta_up_LOS, intensity_dts_LOS, LC_Jacobians_dta_up_LOS,& ! Outputs
           LC_Jacobians_dts_LOS, LS_Jacobians_dts_LOS, tcom1, L_tcom1  )      ! Output

!  Save results

       if ( do_obsgeom ) then
         do g = 1, nvzas
           do lev=1,n_user_levels
             fo_intensity_dta(lev,g,upidx) = intensity_dta_up_LOS(lev,g)
             fo_intensity_dts(lev,g)       = intensity_dts_LOS(lev,g)
             fo_intensity(lev,g,upidx)     = fo_intensity_dta(lev,g,upidx) + fo_intensity_dts(lev,g)
           enddo
         enddo
       else
         do nv = 1, nvzas
           do ns = 1, nszas
             do na = 1, nazms
               g = na_offset(ns,nv) + na
               do lev=1,n_user_levels
                 fo_intensity_dta(lev,g,upidx) = intensity_dta_up_LOS(lev,nv)
                 fo_intensity_dts(lev,g)       = intensity_dts_LOS(lev,nv)
                 fo_intensity(lev,g,upidx)     = fo_intensity_dta(lev,g,upidx) + fo_intensity_dts(lev,g)
               enddo
             enddo
           enddo
         enddo
       endif

       if ( do_columnwfs ) then
         if ( do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_dta(par,lev,g,upidx) = LC_Jacobians_dta_up_LOS(lev,g,par)
                 fo_columnwf_dts(par,lev,g) = LC_Jacobians_dts_LOS(lev,g,par)
                 fo_columnwf(par,lev,g,upidx)    = &
                    fo_columnwf_dta(par,lev,g,upidx) + fo_columnwf_dts(par,lev,g)
               enddo
             enddo
           enddo
         else
           do nv = 1, nvzas
             do ns = 1, nszas
               do na = 1, nazms
                 g = na_offset(ns,nv) + na
                 do lev=1,n_user_levels
                   do par=1,n_columnwfs
                     fo_columnwf_dta(par,lev,g,upidx) = LC_Jacobians_dta_up_LOS(lev,nv,par)
                     fo_columnwf_dts(par,lev,g)       = LC_Jacobians_dts_LOS(lev,nv,par)
                     fo_columnwf(par,lev,g,upidx)    = &
                        fo_columnwf_dta(par,lev,g,upidx) + fo_columnwf_dts(par,lev,g)
                   enddo
                 enddo
               enddo
             enddo
           enddo
         endif
       endif

       if ( do_surfacewfs ) then
         if ( do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do spar=1,n_surfacewfs
                 fo_surfacewf_dts(spar,lev,g) = LS_Jacobians_dts_LOS(lev,g,spar)
                 fo_surfacewf(spar,lev,g)     = fo_surfacewf_dts(spar,lev,g)
               enddo
             enddo
           enddo
         else
           do nv = 1, nvzas
             do ns = 1, nszas
               do na = 1, nazms
                 g = na_offset(ns,nv) + na
                 do spar=1,n_surfacewfs
                   fo_surfacewf_dts(spar,lev,g) = LS_Jacobians_dts_LOS(lev,nv,spar)
                   fo_surfacewf(spar,lev,g)     = fo_surfacewf_dts(spar,lev,g)
                 enddo
               enddo
             enddo
           enddo
         endif
       endif

!  End upwelling

     endif

!  Downwelling
!  -----------

     if ( do_dnwelling ) then

!  Direct thermal, calculate

!mick fix 4/3/2015 - swapped radii & cota_LOS in argument list
       call DTE_Integral_ILCS_DN &
         ( maxvzas, maxlayers, maxfine, max_user_levels, max_atmoswfs,         & ! Inputs (dimensioning)
           Do_Thermset, do_deltam_scaling, do_PlanPar,                         & ! Inputs (Flags)
           do_regular_ps, do_enhanced_ps, doNadir_LOS,                         & ! Inputs (Flags)
           do_columnwfs, n_columnwfs,                                          & ! Inputs (Control, Jacobians)
           nvzas, nlayers, nfinedivs_LOS, n_user_levels, user_levels,          & ! Inputs (control, output)
           bb_input, extinction, deltaus, omega, truncfac,                     & ! Inputs (Optical - Regular)
           L_extinction, L_deltaus, L_omega, L_truncfac,                       & ! Inputs (Optical - Linearized)
           Mu1_LOS, NCrit_LOS, RadCrit_LOS, CotCrit_LOS, Raycon_LOS, cota_LOS, & ! Inputs (Geometry)
           radii, xfine_LOS, wfine_LOS, csqfine_LOS, cotfine_LOS,              & ! Inputs (Geometry)
           intensity_dta_dn_LOS, LC_Jacobians_dta_dn_LOS, tcom1, L_tcom1 )       ! Output

!  Save results

       if ( do_obsgeom ) then
         do g = 1, nvzas
           do lev=1,n_user_levels
             fo_intensity_dta(lev,g,dnidx) = intensity_dta_dn_LOS(lev,g)
             fo_intensity(lev,g,dnidx)     = fo_intensity_dta(lev,g,dnidx)
           enddo
         enddo
       else
         do nv = 1, nvzas
           do ns = 1, nszas
             do na = 1, nazms
               g = na_offset(ns,nv) + na
               do lev=1,n_user_levels
                 fo_intensity_dta(lev,g,dnidx) = intensity_dta_dn_LOS(lev,nv)
                 fo_intensity(lev,g,dnidx)     = fo_intensity_dta(lev,g,dnidx)
               enddo
             enddo
           enddo
         enddo
       endif

       if ( do_columnwfs ) then
          if (do_ObsGeom ) then
           do g = 1, ngeoms
             do lev=1,n_user_levels
               do par=1,n_columnwfs
                 fo_columnwf_dta(par,lev,g,dnidx) = LC_Jacobians_dta_dn_LOS(lev,g,par)
                 fo_columnwf(par,lev,g,dnidx)     = fo_columnwf_dta(par,lev,g,dnidx)
               enddo
             enddo
           enddo
         else
           do nv = 1, nvzas
             do ns = 1, nszas
               do na = 1, nazms
                 g = na_offset(ns,nv) + na
                 do lev=1,n_user_levels
                   do par=1,n_columnwfs
                     fo_columnwf_dta(par,lev,g,dnidx) = LC_Jacobians_dta_dn_LOS(lev,nv,par)
                     fo_columnwf(par,lev,g,dnidx)     = fo_columnwf_dta(par,lev,g,dnidx)
                   enddo
                 enddo
               enddo
             enddo
           enddo
         endif
       endif

!  end downwelling

      endif

!  End thermal run

   endif

!  Solar and Thermal sources run
!  -----------------------------

! Add solar and thermal components

   if ( do_solar_sources.and.(do_thermal_emission.and.do_surface_emission)  ) then

!  Upwelling

     if ( do_upwelling ) then

       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_intensity_atmos(lev,g,upidx) = fo_intensity_ss(lev,g,upidx) + fo_intensity_dta(lev,g,upidx)
           fo_intensity_surf(lev,g)        = fo_intensity_db(lev,g)       + fo_intensity_dts(lev,g)
           fo_intensity(lev,g,upidx)       = fo_intensity_atmos(lev,g,upidx) + fo_intensity_surf(lev,g)
         enddo
       enddo

       if ( do_columnwfs ) then
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do par=1,n_columnwfs
               fo_columnwf_atmos(par,lev,g,upidx) = &
                  fo_columnwf_ss(par,lev,g,upidx) + fo_columnwf_dta(par,lev,g,upidx)
               fo_columnwf_surf(par,lev,g)   = &
                  fo_columnwf_db(par,lev,g) + fo_columnwf_dts(par,lev,g)
               fo_columnwf(par,lev,g,upidx)  = &
                  fo_columnwf_atmos(par,lev,g,upidx) + fo_columnwf_surf(par,lev,g)
             enddo
           enddo
         enddo
       endif

       if ( do_surfacewfs ) then
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do spar=1,n_surfacewfs
               fo_surfacewf_surf(spar,lev,g) = &
                  fo_surfacewf_db(spar,lev,g) + fo_surfacewf_dts(spar,lev,g)
               fo_surfacewf(spar,lev,g) = fo_surfacewf_surf(spar,lev,g)
             enddo
           enddo
         enddo
       endif

     endif

!   Downwelling

     if ( do_dnwelling ) then

       do g = 1, ngeoms
         do lev=1,n_user_levels
           fo_intensity_atmos(lev,g,dnidx) = fo_intensity_ss(lev,g,dnidx) + fo_intensity_dta(lev,g,dnidx)
           fo_intensity(lev,g,dnidx)       = fo_intensity_atmos(lev,g,dnidx)
         enddo
       enddo

       if ( do_columnwfs ) then
         do g = 1, ngeoms
           do lev=1,n_user_levels
             do par=1,n_columnwfs
               fo_columnwf_atmos(par,lev,g,dnidx) = &
                  fo_columnwf_ss(par,lev,g,dnidx) + fo_columnwf_dta(par,lev,g,dnidx)
               fo_columnwf(par,lev,g,dnidx)    = fo_columnwf_atmos(par,lev,g,dnidx)
             enddo
           enddo
         enddo
       endif

     endif

!  End solar+thermal

   endif

!  Finish

   return
end subroutine SFO_LCS_MASTER

end module SFO_LinMasters_m

