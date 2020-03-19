
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

module FirstOrder_geometry_Mk2_20jan12_m

!  Geometry Calculations
!    Assembled 05 December 2011 by R. Spurr, RT SOLUTIONS Inc.
!    revised   16 January  2012 by R. Spurr, RT SOLUTIONS Inc.
!    revised   20 January  2012 by R. Spurr, RT SOLUTIONS Inc. (Dovetailing DTE and SS)

!  PUBLIC Routines: 

!      DTE_GEOMETRY_1; Does LOS outgoing calculations (Up and Down)

!      SS_GEOMETRY_1 ; Does LOS outgoing and Solar Incoming geometrical calculations (Up or  Down)
!      SS_GEOMETRY_2 ; Does LOS outgoing and Solar Incoming geometrical calculations (Up and Down)
!      LEGENDRE_NG   ; Does (for solar beam scattering) Legndre polynomials (Up or  Down)
!      RegularPS_sphergeom  ; Does Regular Pseudo-spherical calculation

!  Private routines

!      FindSunPaths_D
!      FindSunPaths_T
!      FindSun
!      FindSunPath
!      Gauleg_ng

!      STD_outgoing_sphergeom_Initial
!      STD_outgoing_sphergeom_Qbasic

!      SD_incoming_FindCritLayer
!      STD_outgoing_sphergeom_Qadjusted
!      SD_incoming_sphergeom

!  ==================================================================================
!  ==================================================================================

!  DTE_GEOMETRY_1
!  --------------
  
!  [ Regular Pseudo-spherical geometry ]

!       Code included, no calls

!  [ Enhanced Pseudo-spherical geometry ]
!       there are 2 Main Steps, calling the following sequence

!       1. STD_outgoing_sphergeom_Initial

!       2. STD_outgoing_sphergeom_Qbasic
!                  --- Gauleg_ng

!  SS_GEOMETRY_1
!  -------------
  
!  [ Regular Pseudo-spherical geometry ]

!       RegularPS_sphergeom
!                  --- FindSunPaths_D

!  [ Enhanced Pseudo-spherical geometry ]
!       there are 4 Main Steps, calling the following sequence

!       1. STD_outgoing_sphergeom_Initial

!       2. SD_incoming_FindCritLayer
!                  --- FindSunPath

!       3. STD_outgoing_sphergeom_Qadjusted
!                  --- Gauleg_ng

!       4. SD_incoming_sphergeom
!                  --- FindSun
!                  --- FindSunPaths_D
!                  --- FindSunPaths_T

!  SS_GEOMETRY_2
!  -------------

!  [ Regular Pseudo-spherical geometry ]

!       1. RegularPS_sphergeom (Downwelling and upwelling)

!  [ Enhanced Pseudo-spherical geometry ]

!       1.  STD_outgoing_sphergeom_Initial
!       2.  SD_incoming_FindCritLayer
!       3.  STD_outgoing_sphergeom_Qadjusted
!       4a. SD_incoming_sphergeom (Upwelling   solar scattering)
!       4b. SD_incoming_sphergeom (Downwelling solar scattering)

!  =====================================================================================

private
public SS_Geometry_1, SS_Geometry_2, DTE_Geometry_1, Legendre_ng, RegularPS_sphergeom

contains

subroutine DTE_Geometry_1                                      &
       ( maxlayers, maxfine, nlayers, nfinedivs, dtr,          & ! Input
         eradius, heights, alpha_boa, do_regular_ps,           & ! Input
         doNadir, Raycon, radii, alpha, cota,                  & ! Output
         xfine, wfine, csqfine, cotfine, alphafine, radiifine, & ! Output
         fail, message, trace )                                  ! Output

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

      integer  , intent(in)     :: maxlayers, maxfine

!  Layer control. Finelayer divisions may be changed

      integer  , intent(in)     :: nlayers
      integer  , intent(inout)  :: nfinedivs(maxlayers)

! dtr = degrees-to-Radians

      real(fpk), intent(in)     :: dtr

!  Radius + heights

      real(fpk), intent(in)     :: eradius, heights (0:maxlayers)

!  input angle (Degrees)

      real(fpk), intent(InOut)  :: alpha_boa

!  Flags

      logical  , intent(in)     :: do_regular_ps

!  Output arguments
!  ================

!  Flag for the Nadir case

      logical  , intent(inout)  :: doNadir
  
!  Alphas,  Cotangents, Radii, Ray constant

      real(fpk), intent(inout)  :: radii    (0:maxlayers)
      real(fpk), intent(inout)  :: Raycon
      real(fpk), intent(inout)  :: alpha    (0:maxlayers)
      real(fpk), intent(inout)  :: cota     (0:maxlayers)

!  LOS Quadratures for Enhanced PS

      real(fpk), intent(inout)  :: xfine    (maxlayers,maxfine)
      real(fpk), intent(inout)  :: wfine    (maxlayers,maxfine)
      real(fpk), intent(inout)  :: csqfine  (maxlayers,maxfine)
      real(fpk), intent(inout)  :: cotfine  (maxlayers,maxfine)

!  Fine layering output

      real(fpk), intent(inout)  :: alphafine (maxlayers,maxfine)
      real(fpk), intent(inout)  :: radiifine (maxlayers,maxfine)

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message
      character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

      real(fpk)  :: Lospaths (maxlayers)

!  Help variables

      integer    :: n
      real(fpk)  :: alpha_boa_R
      real(fpk), parameter   :: zero  = 0.0

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '

!  check range of inputs

   if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.zero ) then
      message = 'boa LOS angle outside range [0,90]); Check it!'
      trace = 'Initial Angle Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  Regular PS, trivial
!  -------------------

   if ( do_regular_ps ) then
      alpha_boa_R = alpha_boa * dtr
      Raycon      = dsin(alpha_boa_R) * radii(nlayers)
      do n = 0, nlayers
        radii(n) = eradius + heights(n)
        alpha(n) = alpha_boa_R ; cota(n) = zero
      enddo
      return
   endif

!  Enhanced PS; proceed in 2 Steps
!  -------------------------------

!  Step 1; Initial LOS-path quantities
!    Given heights and BOA LOS angle, compute path angles and radii

   CALL STD_outgoing_sphergeom_Initial                           &
       ( maxlayers, nlayers, heights, eradius, alpha_boa, dtr,   & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, cota )           ! Output

!  Step 2. LOS fine-layer quadratures (Regular, non-adjusted)

   CALL STD_outgoing_sphergeom_Qbasic                    &
       ( maxlayers, maxfine, nlayers, nfinedivs,         & ! Input
         doNadir, radii, alpha, Raycon,                  & ! Input
         radiifine, alphafine, xfine, wfine,             & ! Output
         csqfine, cotfine, fail, message )                 ! Output

!  Exception handling

   if ( Fail ) then
      trace = 'Error from STD_outgoing_sphergeom_Qbasic in DTE_Geometry_1' ; return
   endif

!  Finish

   return
end subroutine DTE_Geometry_1

!

subroutine SS_Geometry_1                                                 &
       ( maxlayers, maxfine, nlayers, nfinedivs, dtr, Pie, vsign,        & ! Input
         eradius, heights, alpha_boa, theta_boa, phi_boa,                & ! Input
         do_Chapman, do_LOSpaths, do_regular_ps, doCrit, Acrit, extinc,  & ! Input
         doNadir, Raycon, radii, alpha, cota,                            & ! Input/Output
         xfine, wfine, csqfine, cotfine, alphafine, radiifine,           & ! Input/Output
         NCrit, RadCrit, CotCrit, cosscat, chapfacs,                     & ! Output
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,             & ! Output
         fail, message, trace )                                            ! Output

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

      integer  , intent(in)     :: maxlayers, maxfine

!  Layer control. Finelayer divisions may be changed

      integer  , intent(in)     :: nlayers
      integer  , intent(inout)  :: nfinedivs(maxlayers)

!  Radius + heights

      real(fpk), intent(in)     :: eradius, heights (0:maxlayers)

!  input angles (Degrees), dtr = degrees-to-Radians. VSIGN = +1 (Up); -1(Down)

      real(fpk), intent(in)     :: dtr, Pie, vsign
      real(fpk), intent(InOut)  :: alpha_boa, theta_boa, phi_boa

!  Flags

      logical  , intent(in)     :: do_Chapman
      logical  , intent(in)     :: do_regular_ps

!    LOS paths flag is new, 20 January 2012

      logical  , intent(in)     :: do_LOSpaths

!  Critical adjustment for cloud layers

      logical  , intent(inout)  :: doCrit
      real(fpk), intent(in)     :: extinc(maxlayers)
      real(fpk), intent(in)     :: Acrit

!  Input/Output Arguments (These may already be set if do_LOSpaths )
!  ======================

!  Flag for the Nadir case. Intent(inout), input if DO_LOSpaths set)

      logical  , intent(inout)  :: doNadir
  
!  Alphas,  Cotangents, Radii, Ray constant. Intent(inout), input if DO_LOSpaths set)

      real(fpk), intent(inout)  :: radii    (0:maxlayers)
      real(fpk), intent(inout)  :: Raycon
      real(fpk), intent(inout)  :: alpha    (0:maxlayers)
      real(fpk), intent(inout)  :: cota     (0:maxlayers)

!  LOS Quadratures for Enhanced PS

      real(fpk), intent(inout)  :: xfine    (maxlayers,maxfine)
      real(fpk), intent(inout)  :: wfine    (maxlayers,maxfine)
      real(fpk), intent(inout)  :: csqfine  (maxlayers,maxfine)
      real(fpk), intent(inout)  :: cotfine  (maxlayers,maxfine)

!  Fine layering output

      real(fpk), intent(inout)  :: alphafine (maxlayers,maxfine)
      real(fpk), intent(inout)  :: radiifine (maxlayers,maxfine)

!  Output arguments
!  ================

!  Critical layer

      integer  , intent(out)  :: Ncrit
      real(fpk), intent(out)  :: RadCrit, CotCrit

!  solar paths 

      integer  , Intent(out)  :: ntraverse  (0:maxlayers)
      real(fpk), Intent(out)  :: sunpaths   (0:maxlayers,maxlayers)
      integer  , Intent(out)  :: ntraverse_fine(maxlayers,maxfine)
      real(fpk), Intent(out)  :: sunpaths_fine (maxlayers,maxlayers,maxfine)
      real(fpk), Intent(out)  :: Chapfacs   (maxlayers,maxlayers)

!  Cosin scattering angle

      REAL(fpk), Intent(out)  :: cosscat

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message
      character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

      real(fpk)  :: Lospaths (maxlayers)

!  Other angles

      real(fpk)  :: theta_all  (0:maxlayers)
      real(fpk)  :: phi_all    (0:maxlayers)

!  Critical values

      integer    :: NCritfine
      real(fpk)  :: AlphaCrit

!  Help variables

      real(fpk), parameter   :: zero  = 0.0d0
      real(fpk)  :: term1, term2

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   NCrit = 0 ; NCritfine = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero

!  check range of inputs

   if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
      message = 'boa LOS angle outside range [0,90]); Check it!'
      trace = 'Initial Angle Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

   if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
   if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa

   if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
      message = 'boa SZA angle outside range [0,90]); Check it!'
      trace = 'Initial Angle Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  Regular PS, One routine only
!  ----------------------------

   if ( do_regular_ps ) then
      CALL RegularPS_sphergeom &
       ( maxlayers, nlayers, heights, eradius, do_Chapman,              & ! Inputs
         alpha_boa, theta_boa, phi_boa, vsign, dtr,                     & ! Inputs
         Raycon, radii, alpha, sunpaths, ntraverse,                     & ! Outputs
         chapfacs, cosscat, term1, term2 )                                ! Outputs
      return
   endif

!  Enhanced PS; proceed in 4 Steps
!  -------------------------------

!  Step 1; Initial LOS-path quantities
!    Given heights and BOA LOS angle, compute path angles and radii
!   Only need to do this if LOSpaths flag is not set.

   if ( .not. do_LOSpaths ) then
      CALL STD_outgoing_sphergeom_Initial                        &
       ( maxlayers, nlayers, heights, eradius, alpha_boa, dtr,   & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, cota )           ! Output
   endif

!  Step 2; Find Critical-layer adjustments (Optional)

   if ( doCrit ) then
      CALL SD_incoming_FindCritLayer                            &
       ( maxlayers, nlayers, doNadir, alpha, radii,             & ! Input
         dtr, extinc, theta_boa, Acrit, doCrit,                 & ! Inputs
         Ncrit, NCritfine, AlphaCrit, RadCrit, CotCrit, fail, message )    ! Outputs
      if ( Fail ) then
         trace = 'Error from SD_incoming_FindCritLayer in SS_Geometry_1' ; return
      endif
   endif

!  Step 3. LOS fine-layer quadratures
!           Critical-layer adjustment of Qudarature done here.
!           Regular quadratures may be avoided if do_LOSpaths is set

   CALL STD_outgoing_sphergeom_Qadjusted                 &
       ( maxlayers, maxfine, nlayers, nfinedivs,         & ! Input
         do_LOSpaths, doNadir, radii, alpha, Raycon,     & ! Input
         doCrit, Ncrit, NCritfine, AlphaCrit, RadCrit,   & ! Input
         radiifine, alphafine, xfine, wfine,             & ! Output
         csqfine, cotfine, fail, message )                 ! Output
   if ( Fail ) then
      trace = 'Error from STD_outgoing_sphergeom_Qadjusted in SS_Geometry_1' ; return
   endif

!  Step 4. Solar path lengths

   CALL SD_incoming_sphergeom                                      &
       ( maxlayers, maxfine, nlayers, nfinedivs, do_Chapman,       & ! Input
         DoCrit, NCrit, theta_boa, phi_boa, radii, alpha, vsign,   & ! Input
         dtr, Pie, RadCrit, AlphaCrit, radiifine, alphafine,       & ! Input
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,       & ! Output
         chapfacs, cosscat, theta_all, phi_all )                     ! Output

!  Finish

   return
end subroutine SS_Geometry_1

!

subroutine SS_Geometry_2                                                 &
       ( maxlayers, maxfine, nlayers, nfinedivs, dtr, Pie, vsign,        & ! Input
         eradius, heights, alpha_boa, theta_boa, phi_boa,                & ! Input
         do_Chapman, do_LOSpaths, do_regular_ps, doCrit, Acrit, extinc,  & ! Input
         doNadir, Raycon, radii, alpha, cota,                            & ! Input/Output
         xfine, wfine, csqfine, cotfine, alphafine, radiifine,           & ! Input/Output
         NCrit, RadCrit, CotCrit,                                        & ! Output
         cosscat_up, chapfacs_up, cosscat_dn, chapfacs_dn,               & ! Output
         sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up, & ! Output
         sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn, & ! Output
         fail, message, trace )                                            ! Output

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

      integer  , intent(in)     :: maxlayers, maxfine

!  Layer control. Finelayer divisions may be changed

      integer  , intent(in)     :: nlayers
      integer  , intent(inout)  :: nfinedivs(maxlayers)

!  Radius + heights

      real(fpk), intent(in)     :: eradius, heights (0:maxlayers)

!  input angles (Degrees), dtr = degrees-to-Radians

      real(fpk), intent(in)     :: dtr, Pie
      real(fpk), intent(InOut)  :: alpha_boa, theta_boa, phi_boa

!  Flags

      logical  , intent(in)     :: do_Chapman
      logical  , intent(in)     :: do_regular_ps

!    LOS paths flag is new, 20 January 2012

      logical  , intent(in)     :: do_LOSpaths

!  Critical adjustment for cloud layers

      logical  , intent(inout)  :: doCrit
      real(fpk), intent(in)     :: extinc(maxlayers)
      real(fpk), intent(in)     :: Acrit

!  Input/Output Arguments (These may already be set if do_LOSpaths )
!  ======================

!  Flag for the Nadir case. Intent(inout), input if DO_LOSpaths set)

      logical  , intent(inout)  :: doNadir
  
!  Alphas,  Cotangents, Radii, Ray constant. Intent(inout), input if DO_LOSpaths set)

      real(fpk), intent(inout)  :: radii    (0:maxlayers)
      real(fpk), intent(inout)  :: Raycon
      real(fpk), intent(inout)  :: alpha    (0:maxlayers)
      real(fpk), intent(inout)  :: cota     (0:maxlayers)

!  LOS Quadratures for Enhanced PS

      real(fpk), intent(inout)  :: xfine    (maxlayers,maxfine)
      real(fpk), intent(inout)  :: wfine    (maxlayers,maxfine)
      real(fpk), intent(inout)  :: csqfine  (maxlayers,maxfine)
      real(fpk), intent(inout)  :: cotfine  (maxlayers,maxfine)

!  Fine layering output

      real(fpk), intent(inout)  :: alphafine (maxlayers,maxfine)
      real(fpk), intent(inout)  :: radiifine (maxlayers,maxfine)

!  Output arguments
!  ================

!  Critical layer

      integer  , intent(out)  :: Ncrit
      real(fpk), intent(out)  :: RadCrit, CotCrit

!  solar paths upwelling

      integer  , Intent(out)  :: ntraverse_up  (0:maxlayers)
      real(fpk), Intent(out)  :: sunpaths_up   (0:maxlayers,maxlayers)
      integer  , Intent(out)  :: ntraverse_fine_up(maxlayers,maxfine)
      real(fpk), Intent(out)  :: sunpaths_fine_up (maxlayers,maxlayers,maxfine)
      real(fpk), Intent(out)  :: Chapfacs_up  (maxlayers,maxlayers)

!  solar paths downwelling

      integer  , Intent(out)  :: ntraverse_dn  (0:maxlayers)
      real(fpk), Intent(out)  :: sunpaths_dn   (0:maxlayers,maxlayers)
      integer  , Intent(out)  :: ntraverse_fine_dn(maxlayers,maxfine)
      real(fpk), Intent(out)  :: sunpaths_fine_dn (maxlayers,maxlayers,maxfine)
      real(fpk), Intent(out)  :: Chapfacs_dn  (maxlayers,maxlayers)

!  Cosine scattering angles

      REAL(fpk), Intent(out)  :: cosscat_up, cosscat_dn

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message
      character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

      real(fpk)  :: Lospaths (maxlayers)

!  Other angles

      real(fpk)  :: theta_all  (0:maxlayers)
      real(fpk)  :: phi_all    (0:maxlayers)

!  Critical values

      integer    :: NCritfine
      real(fpk)  :: AlphaCrit

!  Help variables
!    VSIGN = +1 for Upwelling, -1 for Downwelling

      real(fpk), parameter   :: zero  = 0.0d0
      real(fpk)  :: vsign, term1, term2

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   NCrit = 0 ; NCritfine = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero

!  check range of inputs

   if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
      message = 'boa LOS angle outside range [0,90]); Check it!'
      trace = 'Initial Angle Check in SS_Geometry_2'
      fail    = .true. ;  return
   endif

   if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
   if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa

   if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
      message = 'boa SZA angle outside range [0,90]); Check it!'
      trace = 'Initial Angle Check in SS_Geometry_2'
      fail    = .true. ;  return
   endif

!  Regular PS, One routine only
!  ----------------------------

!  Single call to upwelling (VSIGN = +1)
!  Downwelling only requires new COSSCAT, rest is copied

   if ( do_regular_ps ) then
      vsign = +1.0d0
      CALL RegularPS_sphergeom &
       ( maxlayers, nlayers, heights, eradius, do_Chapman,              & ! Inputs
         alpha_boa, theta_boa, phi_boa, vsign, dtr,                     & ! Inputs
         Raycon, radii, alpha, sunpaths_up, ntraverse_up,               & ! Outputs
         chapfacs_up, cosscat_up, term1, term2 )                          ! Outputs
      vsign = -1.0d0  ;  cosscat_dn = -vsign * term1 + term2
      ntraverse_dn = ntraverse_up
      sunpaths_dn  = sunpaths_up
      chapfacs_dn  = chapfacs_up
      return
   endif

!  Enhanced PS; proceed in 4 Steps
!  -------------------------------

!  Step 1; Basic LOS quantities
!    Given heights and BOA LOS angle, compute path angles and radii
!   Only need to do this if LOSpaths flag is not set.

   if ( .not. do_LOSpaths ) then
      CALL STD_outgoing_sphergeom_Initial                        &
       ( maxlayers, nlayers, heights, eradius, alpha_boa, dtr,   & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, cota )           ! Output
   endif

!  Step 2; Find Critical-layer adjustments (Optional)

   if ( doCrit ) then
      CALL SD_incoming_FindCritLayer                            &
       ( maxlayers, nlayers, doNadir, alpha, radii,             & ! Input
         dtr, extinc, theta_boa, Acrit, doCrit,                 & ! Inputs
         Ncrit, NCritfine, AlphaCrit, RadCrit, CotCrit, fail, message )    ! Outputs
      if ( Fail ) then
         trace = 'Error from SD_incoming_FindCritLayer in SS_Geometry_2' ; return
      endif
   endif

!  Step 3. LOS fine-layer quadratures
!           Regular quadratures may be avoided if do_LOSpaths is set

   CALL STD_outgoing_sphergeom_Qadjusted                 &
       ( maxlayers, maxfine, nlayers, nfinedivs,         & ! Input
         do_LOSpaths, doNadir, radii, alpha, Raycon,     & ! Input
         doCrit, Ncrit, NCritfine, AlphaCrit, RadCrit,   & ! Input
         radiifine, alphafine, xfine, wfine,             & ! Output
         csqfine, cotfine, fail, message )                 ! Output
   if ( Fail ) then
      trace = 'Error from STD_outgoing_sphergeom_Qadjusted in SS_Geometry_2' ; return
   endif

!  Step 4a. Solar path lengths, Upwelling

   vsign = + 1.0d0
   CALL SD_incoming_sphergeom                                            &
       ( maxlayers, maxfine, nlayers, nfinedivs, do_Chapman,             & ! Input
         DoCrit, NCrit, theta_boa, phi_boa, radii, alpha, vsign,         & ! Input
         dtr, Pie, RadCrit, AlphaCrit, radiifine, alphafine,             & ! Input
         sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up, & ! Output
         chapfacs_up, cosscat_up, theta_all, phi_all )                     ! Output

!  Step 4b. Solar path lengths, Downwelling

   vsign = - 1.0d0
   CALL SD_incoming_sphergeom                                            &
       ( maxlayers, maxfine, nlayers, nfinedivs, do_Chapman,             & ! Input
         DoCrit, NCrit, theta_boa, phi_boa, radii, alpha, vsign,         & ! Input
         dtr, Pie, RadCrit, AlphaCrit, radiifine, alphafine,             & ! Input
         sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn, & ! Output
         chapfacs_dn, cosscat_dn, theta_all, phi_all )                     ! Output

!  Finish

   return
end subroutine SS_Geometry_2

!

subroutine STD_outgoing_sphergeom_Initial                        &
       ( maxlayers, nlayers, heights, eradius, alpha_boa, dtr,   & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, cota )           ! Output

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!  This routine: Initial LOS path setup

!    starting inputs are - BOA value of VZA (alpha_boa), in degrees
!                        - height grid, earth radius

      implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

      integer  , intent(in)   :: maxlayers

!  Layer control

      integer  , intent(in)   :: nlayers
      real(fpk), intent(in)   :: eradius, heights (0:maxlayers)

!  input angle

      real(fpk), intent(in)   :: alpha_boa, dtr

!  Flag for the Nadir case

      logical  , intent(out)  :: doNadir
  
!  Alphas, Radii, Ray constant, Lospaths

      real(fpk), intent(out)  :: radii    (0:maxlayers)
      real(fpk), intent(out)  :: Raycon
      real(fpk), intent(out)  :: Lospaths (maxlayers)
      real(fpk), intent(out)  :: alpha    (0:maxlayers)
      real(fpk), intent(out)  :: cota     (0:maxlayers)

!  Local
!  -----

      integer         :: n, n1
      real(fpk)       :: salpha_boa, difh, alpha_boa_R
      real(fpk)       :: salpha, calpha, calpha1
      real(fpk), parameter :: zero = 0.0d0

!  Zero output

      Lospaths = zero ; cota = zero ; alpha = zero ; Raycon = zero

!  Radii

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Special case

      doNadir = .false.
      if ( alpha_boa.eq.0.0d0 ) doNadir = .true.

!  Special case. Direct nadir viewing. Compute everything and Exit.

      if ( doNadir ) then
        do n = nlayers,1,-1
          difh = radii(n-1) - radii(n) ; Lospaths(n) = difh
        enddo
        return
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  start at BOA

      alpha_boa_R    = alpha_boa * DTR
      salpha_boa     = dsin(alpha_boa_R)
      calpha1        = dcos(alpha_boa_R)
      cota(nlayers)  = calpha1 / salpha_boa
      alpha(nlayers) = alpha_boa_R

!  Ray constant

      Raycon = salpha_boa * radii(nlayers)

!  Whole layer values

      do n = nlayers - 1, 0, -1
         n1 = n + 1
         salpha  = Raycon / radii(n) ; alpha(n) = dasin(salpha)
         calpha  = dcos(alpha(n))
         cota(n) = calpha / salpha
         Lospaths(n1) = radii(n)*calpha - radii(n1)*calpha1
         calpha1 = calpha
      enddo

!  Finish

      return
end subroutine STD_outgoing_sphergeom_initial

!

subroutine SD_incoming_FindCritLayer                            &
       ( maxlayers, nlayers, doNadir, alpha, radii,             & ! Input
         dtr, extinctions, theta_boa, Acrit, doCrit,            & ! Inputs
         Ncrit, NCritfine, AlphaCrit, RadCrit, CotCrit, fail, message )    ! Outputs

!  Purpose: Given a list of Maximum extinctions and solar angle at BOA
!              Then find Critical layer (NCrit) and point where attenuation wipe-out (Acrit) is achieved
!              Then find the LOS angle and Radius (AlphaCrit,RadCrit) for this Critical Point

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ------

!  Dimensioning

   integer, intent(in) :: maxlayers

!  Layer control

   integer, intent(in) :: nlayers

!  Special case, Nadir viewing

   logical, intent(in)  :: doNadir

!  View angles and Radii at layer boundaries

   real(fpk), intent(in)  :: alpha(0:maxlayers)
   real(fpk), intent(in)  :: radii(0:maxlayers)

!  Extinctions

   real(fpk), intent(in)  :: extinctions(maxlayers)

!  Solar control and angle

   real(fpk), intent(in)  :: Acrit, theta_boa, dtr

!  Overall control (May be switched off if Critical test is negative)

   logical, intent(inout) :: doCrit

!  outputs
!  -------

!  Critical layer, Number of Fine divisions for this layer

   integer  , intent(out)  :: Ncrit
   integer  , intent(out)  :: Ncritfine

!  Critical angle and radius and cotangent

   real(fpk), intent(out)  :: AlphaCrit
   real(fpk), intent(out)  :: RadCrit, CotCrit

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

!  Solar direction and critical path

   real(fpk)  :: SolarDirection(3)
   real(fpk)  :: CritSunPath

!  Bisection accuracy

   real(fpk), parameter  :: BisectionAccuracy = 1.0d-10
   integer   , parameter :: jmax = 40

!  Other variables

   logical    ::  Finding, NonCrit
   integer    ::  J, N
   real(fpk)  ::  xboa, x1, x2, xmid, rtop, raycon, sx, arg, attn, t0, t1, cons
   real(fpk)  ::  dist2, distmid, theta0, rtbis, dx, f, fmid, s0

!  Initialize

   Ncrit     = 0
   NCritfine = 0
   AlphaCrit = 0.0d0
   RadCrit   = 0.0d0
   CotCrit   = 0.0d0

   fail    = .false.
   message = ' '

!  Solar direction chosen to be maximum in the forward direction

   t0 = theta_boa * dtr
   SolarDirection(1) = - dsin(t0)
   SolarDirection(2) = 0.0d0
   SolarDirection(3) = - dcos(t0)

!  set number of fine layers

!   Ncritfine = 10
!   if ( Acrit .gt. 1.0d-07 ) then
!      Ncritfine = 7
!   else
!      if ( Acrit .gt. 1.0d-08 ) NCritfine = 8
!   endif

   Ncritfine = INT(-DLOG10(Acrit))

!  Special case
!  ------------

   if ( doNadir ) then

!  Trawl

      NonCrit = .true.; n = 0 ; s0 = -SolarDirection(1)
      do while (n.lt.nlayers .and. NonCrit)
         n = n + 1
         Cons  = radii(n-1) / s0 ; t1 = dasin(radii(n)/cons)
         dist2 = dsin(t0-t1)*cons
         Arg = dist2 * extinctions(n)
         Attn = 0.0d0 ;  if (arg .lt. 80.0d0 ) Attn = dexp ( - Arg )
         if ( Attn .lt. Acrit ) NonCrit = .false.
      enddo

!  return if no critical layer

      doCrit = (.not.NonCrit) ; if ( NonCrit ) return

!  Set Critical layer and sunpath distance

      if ( doCrit ) then
         Ncrit = n
         CritSunPath = -dlog(Acrit)/extinctions(n)
         t1 = t0 - dasin(CritSunPath /cons)
         RadCrit = dsin(t1)*cons
      endif

!  Return special case

      return

   endif

!  General case
!  ------------

!  setups

   xboa   = alpha(nlayers)
   raycon = radii(nlayers) * dsin(xboa)

!  Find NCrit

   NonCrit = .true.; n = 0
   do while (n.lt.nlayers .and. NonCrit)
      n = n + 1
      call FindSunPath ( alpha(n), xboa, radii(n-1), raycon, SolarDirection, dist2, theta0 )
      Arg = dist2 * extinctions(n)
      Attn = 0.0d0 ;  if (arg .lt. 80.0d0 ) Attn = dexp ( - Arg )
      if ( Attn .lt. Acrit ) NonCrit = .false.
   enddo

!  return if no critical layer

   doCrit = (.not.NonCrit) ; if ( NonCrit ) return

!  Set Critical layer and sunpath distance

   if ( doCrit ) then
      Ncrit = n
      CritSunPath = -dlog(Acrit)/extinctions(n)
   endif

!  Now find the Critical angle and Radius for this Length in layer Ncrit
!   ++ Use Bisection Function F(x) = SunPath(x) - CritSunPath

   rtop   = radii(NCrit-1)

!  Lowest value of Function (layer top)

   x1    = alpha(NCrit-1)
   F     = - CritSunPath

!  Highest value of function (layer bottom)

   x2   = alpha(NCrit)
   FMID = dist2 - CritSunPath

!  No bisection

   IF(F*FMID.GE.0.) then
      fail = .true. ; message = 'Root must be bracketed for bisection.' ; return
   ENDIF

!  Bisect the function setup

   IF(F.LT.0.)THEN
      RTBIS=X1
      DX=X2-X1
   ELSE
      RTBIS=X2
      DX=X1-X2
   ENDIF

!  Iterate

   Finding = .true. ; J = 0
   DO While (Finding .and. j .lt. JMAX)
      J = J + 1 ; dx = 0.5d0 * dx ; XMID = RTBIS + DX
      call FindSunPath ( xmid, xboa, rtop, raycon, SolarDirection, distmid, theta0 )
!      write(*,*)J,xmid*180.0d0/dacos(-1.0d0),distmid-CritSunPath
      FMID = distmid - CritSunPath
      IF ( FMID.LE.0.0d0 ) RTBIS = XMID
      IF(DABS(DX).LT.BisectionAccuracy .OR. FMID.EQ.0.) Finding = .false.
   ENDDO

!  Exception (too many bisections)

   if ( Finding ) Then
      fail = .true.
      message = 'Too many Bisections (>40); Root not found'
      return
   endif

!  Set final output if successful

   AlphaCrit = RTBIS ;  SX = dsin(AlphaCrit)
   RadCrit   = Raycon / SX
   CotCrit   = dsqrt(1.0d0-SX*SX)/SX

!  Finish

   RETURN
END subroutine SD_incoming_FindCritLayer
  
!

subroutine STD_outgoing_sphergeom_Qadjusted               &
       ( maxlayers, maxfine, nlayers, nfinedivs,         & ! Input
         do_LOSpaths, doNadir, radii, alpha, Raycon,     & ! Input
         doCrit, Ncrit, NCritfine, AlphaCrit, RadCrit,   & ! Input
         radii_fine, alpha_fine, xfine, wfine,           & ! Output/Input
         csqfine, cotfine, fail, message )                 ! Output/Input

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!    starting inputs are - BOA value of VZA (alpha_boa), in degrees
!                        - height grid, earth radius, Layer control
!                        - Critical Layer control

!  Regular Quadrature need not be done if LOSPATHS is set

      implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

      integer, intent(in)          :: maxlayers, maxfine

!  Layer control. Finelayer divisions may be changed

      integer, intent(in)          :: nlayers
      integer, intent(inout)       :: nfinedivs(maxlayers)

!  Flag for the Nadir case

      logical  , intent(in)     :: doNadir
  
!  Flag for pre-existing calculations

      logical  , intent(in)     :: do_LOSpaths
  
!  Alphas, Radii, Ray constant

      real(fpk), intent(in)  :: alpha      (0:maxlayers)
      real(fpk), intent(in)  :: radii      (0:maxlayers)
      real(fpk), intent(in)  :: Raycon

!  Critical stuff

   logical  , intent(in)  :: doCrit
   integer  , intent(in)  :: Ncrit
   integer  , intent(in)  :: Ncritfine
   real(fpk), intent(in)  :: AlphaCrit
   real(fpk), intent(in)  :: RadCrit

!  Outputs
!  =======

!  Fine layering

      real(fpk), intent(out)  :: alpha_fine (maxlayers,maxfine)
      real(fpk), intent(out)  :: radii_fine (maxlayers,maxfine)

!  Quadratures

      real(fpk), intent(out)  :: xfine    (maxlayers,maxfine)
      real(fpk), intent(out)  :: wfine    (maxlayers,maxfine)

!  Local geoemetry arrays

      real(fpk), intent(out)  :: csqfine  (maxlayers,maxfine)
      real(fpk), intent(out)  :: cotfine  (maxlayers,maxfine)

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message

!  Local
!  -----

      integer, parameter :: Localmaxfine = 20
      integer            :: n, n1, j, nfine
      real(fpk)          :: difh, csfine
      real(fpk)          :: tfine(Localmaxfine), afine(Localmaxfine)

      real(fpk), parameter :: zero = 0.0d0

!  Zero output

      if ( .not. do_LOSpaths ) then
         alpha_fine = zero    ; radii_fine = zero
         xfine      = zero    ; wfine      = zero
         cotfine    = zero    ; csqfine    = zero
      endif

      fail       = .false. ; message    = ' '

!  Check

      if ( maxfine.gt.localmaxfine ) then
         fail = .true.
         message = 'Failure! Localmaxfine dimension not large enough'
         return
      endif

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit. Qudratures are height-oriented
!    (This should be the same as the regular pseudo-spherical )

      if ( doNadir ) then

!  For normal atmosphere, Regular quadrature only if flagged

         if ( .not. doCrit .or. NCrit .eq. 0 ) then
            if ( .not. do_LOSpaths ) then
               do n = nlayers,1,-1
                  difh  = radii(n-1) - radii(n)
                  nfine = nfinedivs(n)
                  call gauleg_ng (0.0d0,difh,tfine,afine,nfine,localmaxfine)
                  do j = 1, nfine
                     radii_fine(n,j) = radii(n-1) - tfine(j)
                     xfine(n,j) = tfine(j)
                     wfine(n,j) = afine(j)
                  enddo
               enddo
            endif

!  Otherwise.....

         else

!    -- Adjust quadrature for the Critical layer
            
            n = NCrit ; difh = radii(n-1) - Radcrit ; nfine = NCritfine ; nfinedivs(n) = nfine
            call gauleg_ng (0.0d0,difh,tfine,afine,nfine,localmaxfine)
            do j = 1, nfine
               radii_fine(n,j) = radii(n-1) - tfine(j)
               xfine(n,j) = tfine(j)
               wfine(n,j) = afine(j)
            enddo

!    -- For all layers above Critical layer, Regular quadrature only if flagged

            if ( .not. do_LOSpaths ) then
               do n = NCrit+1,1,-1
                  difh  = radii(n-1) - radii(n) ; nfine = nfinedivs(n)
                  call gauleg_ng (0.0d0,difh,tfine,afine,nfine,localmaxfine)
                  do j = 1, nfine
                     radii_fine(n,j) = radii(n-1) - tfine(j)
                     xfine(n,j) = tfine(j)
                     wfine(n,j) = afine(j)
                  enddo
               enddo
            endif
         endif

!  Done Nadir case so return

         return
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  For normal atmosphere, Regular quadrature only if flagged

      if ( .not. doCrit .or. NCrit .eq. 0 ) then
         if ( .not. do_LOSpaths ) then
            do n = nlayers, 1, -1
               n1 = n - 1
               nfine = nfinedivs(n)
               call gauleg_ng (alpha(n1),alpha(n),tfine,afine,nfine,localmaxfine)
               do j = 1,  nfine
                  csfine = 1.0d0 / dsin(tfine(j))
                  radii_fine(n,j) = raycon * csfine
                  alpha_fine(n,j) = tfine(j)
                  xfine(n,j)   = radii(n1) - radii_fine(n,j)
                  wfine(n,j)   = afine(j)
                  cotfine(n,j) = dcos(tfine(j)) * csfine
                  csqfine(n,j) = csfine * csfine
               enddo
            enddo
         endif

!  Otherwise

      else

!    -- Adjust quadrature for the Critical layer

         n = NCrit ; n1 = n - 1 ; nfine = NCritfine ; nfinedivs(n) = nfine
         call gauleg_ng (alpha(n1),AlphaCrit,tfine,afine,nfine,localmaxfine)
         do j = 1,  nfine
            csfine = 1.0d0 / dsin(tfine(j))
            radii_fine(n,j) = raycon * csfine
            alpha_fine(n,j) = tfine(j)
            xfine(n,j)   = radii(n1) - radii_fine(n,j)
            wfine(n,j)   = afine(j)
            cotfine(n,j) = dcos(tfine(j)) * csfine
            csqfine(n,j) = csfine * csfine
         enddo

!    -- For all layers above Critical layer, Regular quadrature only if flagged

         if ( .not. do_LOSpaths ) then
            do n = NCrit+1,1,-1
               n1 = n - 1
               nfine = nfinedivs(n)
               call gauleg_ng (alpha(n1),alpha(n),tfine,afine,nfine,localmaxfine)
               do j = 1,  nfine
                  csfine = 1.0d0 / dsin(tfine(j))
                  radii_fine(n,j) = raycon * csfine
                  alpha_fine(n,j) = tfine(j)
                  xfine(n,j)   = radii(n1) - radii_fine(n,j)
                  wfine(n,j)   = afine(j)
                  cotfine(n,j) = dcos(tfine(j)) * csfine
                  csqfine(n,j) = csfine * csfine
               enddo
            enddo
         endif

!  Done

      endif

!  Finish

      return
end subroutine STD_outgoing_sphergeom_Qadjusted

!
subroutine STD_outgoing_sphergeom_Qbasic                  &
       ( maxlayers, maxfine, nlayers, nfinedivs,         & ! Input
         doNadir, radii, alpha, Raycon,                  & ! Input
         radii_fine, alpha_fine, xfine, wfine,           & ! Output
         csqfine, cotfine, fail, message )                 ! Output

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!    starting inputs are - BOA value of VZA (alpha_boa), in degrees
!                        - height grid, earth radius, Layer control

      implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

      integer, intent(in)          :: maxlayers, maxfine

!  Layer control. Finelayer divisions is Strictly input

      integer, intent(in)       :: nlayers
      integer, intent(in)       :: nfinedivs(maxlayers)

!  Flag for the Nadir case

      logical  , intent(in)     :: doNadir
  
!  Alphas, Radii, Ray constant

      real(fpk), intent(in)  :: alpha      (0:maxlayers)
      real(fpk), intent(in)  :: radii      (0:maxlayers)
      real(fpk), intent(in)  :: Raycon

!  Outputs
!  =======

!  Fine layering

      real(fpk), intent(out)  :: alpha_fine (maxlayers,maxfine)
      real(fpk), intent(out)  :: radii_fine (maxlayers,maxfine)

!  Quadratures

      real(fpk), intent(out)  :: xfine    (maxlayers,maxfine)
      real(fpk), intent(out)  :: wfine    (maxlayers,maxfine)

!  Local geoemetry arrays

      real(fpk), intent(out)  :: csqfine  (maxlayers,maxfine)
      real(fpk), intent(out)  :: cotfine  (maxlayers,maxfine)

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message

!  Local
!  -----

      integer, parameter :: Localmaxfine = 20
      integer            :: n, n1, j, nfine
      real(fpk)          :: difh, csfine
      real(fpk)          :: tfine(Localmaxfine), afine(Localmaxfine)

      real(fpk), parameter :: zero = 0.0d0

!  Zero output

      alpha_fine = zero    ; radii_fine = zero
      xfine      = zero    ; wfine      = zero
      cotfine    = zero    ; csqfine    = zero
      fail       = .false. ; message    = ' '

!  Check

      if ( maxfine.gt.localmaxfine ) then
         fail = .true.
         message = 'Failure! Localmaxfine dimension not large enough'
         return
      endif

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit. Qudratures are height-oriented
!    (This should be the same as the regular pseudo-spherical )

      if ( doNadir ) then
         do n = nlayers,1,-1
            difh  = radii(n-1) - radii(n)
            nfine = nfinedivs(n)
            call gauleg_ng (0.0d0,difh,tfine,afine,nfine,localmaxfine)
            do j = 1, nfine
                radii_fine(n,j) = radii(n-1) - tfine(j)
               xfine(n,j) = tfine(j)
               wfine(n,j) = afine(j)
            enddo
         enddo
         return
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  Whole layer values

      do n = nlayers, 1, -1
         n1 = n - 1
         nfine = nfinedivs(n)
         call gauleg_ng (alpha(n1),alpha(n),tfine,afine,nfine,localmaxfine)
         do j = 1,  nfine
            csfine = 1.0d0 / dsin(tfine(j))
            radii_fine(n,j) = raycon * csfine
            alpha_fine(n,j) = tfine(j)
            xfine(n,j)   = radii(n1) - radii_fine(n,j)
            wfine(n,j)   = afine(j)
            cotfine(n,j) = dcos(tfine(j)) * csfine
            csqfine(n,j) = csfine * csfine
         enddo
      enddo

!  Finish

      return
end subroutine STD_outgoing_sphergeom_Qbasic

!

subroutine SD_incoming_sphergeom                                 &
       ( maxlayers, maxfine, nlayers, nfinedivs, do_Chapman,       & ! Input
         DoCrit, NCrit, theta_boa, phi_boa, radii, alpha, vsign,   & ! Input
         dtr, Pie, RadCrit, AlphaCrit, radii_fine, alpha_fine,     & ! Input
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,       & ! Output
         chapfacs, cosscat, theta_all, phi_all )                     ! Output

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the incoming Solar Beam
!     This is applicable to Both Upwelling and Downwelling LOS-path geometries
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control
!    need also the complete values of all VZAs along outgoing path

!  This routine has the fine gridding treatment

      implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  inputs

      integer  , intent(In)    :: maxlayers, maxfine
      integer  , intent(In)    :: nlayers, nfinedivs(maxlayers), NCrit
      real(fpk), intent(InOut) :: theta_boa, phi_boa
      real(fpk), intent(In)    :: vsign
      logical  , intent(In)    :: do_Chapman, DoCrit

!  Los geometry

      real(fpk), intent(In)   :: alpha         (0:maxlayers)
      real(fpk), intent(In)   :: alpha_fine    (maxlayers,maxfine)
      real(fpk), intent(In)   :: radii         (0:maxlayers)
      real(fpk), intent(In)   :: radii_fine    (maxlayers,maxfine)
      real(fpk), intent(In)   :: AlphaCrit, RadCrit, dtr, pie

!  main outputs (geometry)

      integer  , intent(Out)  :: ntraverse  (0:maxlayers)
      real(fpk), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers)
      real(fpk), intent(Out)  :: chapfacs   (maxlayers,maxlayers)

!  Fine level output (geometry)

      integer  , intent(Out)  :: ntraverse_fine(maxlayers,maxfine)
      real(fpk), intent(Out)  :: sunpaths_fine (maxlayers,maxlayers,maxfine)

!  scattering angle and associated angles

      real(fpk), intent(Out)  :: cosscat
      real(fpk), intent(Out)  :: theta_all  (0:maxlayers)
      real(fpk), intent(Out)  :: phi_all    (0:maxlayers)

!  Local

      logical       :: DirectSun, Do_OverheadSun, Do_ZeroSunBOA, Do_Normal
      integer       :: n, j, k
      real(fpk)     :: SolarDirection(3), Radstart
      real(fpk)     :: salpha_boa, calpha_boa, phi_boa_R, sphi_boa
      real(fpk)     :: theta_boa_R, stheta_boa, ctheta_boa, cphi_boa
      real(fpk)     :: ctheta, stheta, calpha, salpha, cphi, CumAngle

!  Check
!      real(fpk)  :: sumd, sume, sth1

!  Local arrays associated with fine grid output

      logical            :: DirectSunf(maxfine)
      real(fpk)       :: thetaf(maxfine)
      real(fpk)       :: sthetaf(maxfine)
      real(fpk)       :: cthetaf(maxfine)

!  Initialise output

      ntraverse = 0     ; ntraverse_fine = 0
      sunpaths  = 0.0d0 ; sunpaths_fine  = 0.0d0
      chapfacs  = 0.0d0

      phi_all = 0.0d0   ; theta_all = 0.0d0 ; cosscat = 0.0d0

!  Nominal traverse paths for Full illumination

      ntraverse(0) = 0
      do n = 1, nlayers
         ntraverse(n) = n
         do j = 1, nfinedivs(n)
            ntraverse_fine(n,j) = n
         enddo
      enddo

!  check range of inputs

!      if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
!        message = 'boa LOS angle outside range [0,90]); Check it!'
!        fail    = .true.
!        return
!      endif
!      if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
!      if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa
!      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
!        message = 'boa SZA angle outside range [0,90]); Check it!'
!        fail    = .true.
!        return
!      endif
!      if ( maxfine.gt.maxlocalfine ) then
!         message = 'local finelayer dimensioning insufficient'
!         fail    = .true.
!         return
!      endif 

!  start at BOA

      theta_boa_R = theta_boa * dtr
      phi_boa_R   = phi_boa * dtr
      salpha_boa  = dsin(alpha(nlayers))
      calpha_boa  = dcos(alpha(nlayers))
      stheta_boa  = dsin(theta_boa_R)
      ctheta_boa  = dcos(theta_boa_R)
      cphi_boa    = dcos(phi_boa_R)
      sphi_boa    = dsin(phi_boa_R)

!  Cosine of scattering angle at boa

      cosscat = - vsign * calpha_boa * ctheta_boa + salpha_boa * stheta_boa * cphi_boa 

!  General case: LOS path in spherical geometry
!  ============================================

!  define Unit solar vector at BOA

      SolarDirection(1) = - stheta_boa * cphi_boa * vsign
      SolarDirection(2) = - stheta_boa * sphi_boa
      SolarDirection(3) = - ctheta_boa

!  Special cases

      Do_OverheadSun = stheta_boa.eq.0.0d0 

!  Start loop over all layers

      do n = nlayers, 1, -1

!  Special cases

        DO_ZeroSunBOA  = (Do_OverheadSun .and. n.eq. nlayers )
        DO_Normal      = .not. doCrit .or. ( doCrit .and. n.le. NCrit )

!  Layer boundary Sun position
!     * Local save of angles, cosines, sines and  illumination flags
!     * Use critical ALPHA and RADIUS if N = NCrit
!     * Use Bottom-of-layer values if N < NCrit (BOA values if illuminated)

        if ( Do_Normal ) then
           if ( doCrit .and. n .eq. NCrit ) then
              CumAngle = alpha(nlayers) - AlphaCrit ; Radstart = RadCrit
              call FindSun(Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_all(n),stheta,ctheta,DirectSun)
           else
              Radstart = radii(n)
              if ( n.eq. nlayers ) then
                 theta_all(n) = theta_boa*dtr ; stheta = stheta_boa ; ctheta = ctheta_boa ; DirectSun = .true.
              else
                 CumAngle = alpha(nlayers) - alpha(n)
                 call FindSun(Do_OverheadSun,radii(n),SolarDirection,CumAngle,theta_all(n),stheta,ctheta,DirectSun)
              endif
           endif
        endif

!  Fine-layer sun positions

        if ( Do_Normal ) then
           do j = 1, nfinedivs(n)
              CumAngle = alpha(nlayers) - alpha_fine(n,j)
              call FindSun(Do_OverheadSun,radii_fine(n,j),SolarDirection,CumAngle,thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
           enddo
        endif

!  Sun paths in layer

        if ( Do_Normal ) then
           if ( DirectSun ) then
              call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,theta_all(n),stheta,N,sunpaths(n,:))
           else
              call FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,theta_all(n),stheta,N,sunpaths(n,:),ntraverse(n))
           endif
           do j = 1, nfinedivs(n) 
              if ( DirectSunf(j) ) then
                 call FindSunPaths_D &
                  (Do_ZeroSunBOA,Maxlayers,Radii_fine(n,j),Radii,thetaf(j),sthetaf(j),N,sunpaths_fine(n,:,J))
              else
                 call FindSunPaths_T &
                  (Maxlayers,Pie,Radii_fine(n,j),Radii,thetaf(j),sthetaf(j),N,sunpaths_fine(n,:,J),ntraverse_fine(n,J))
              endif
!             if ( n.eq.14 ) write(*,*)j,n,Radii_fine(n,j)-radii(n)
           enddo
        endif

!  debugging

!        if ( n.eq.14) then
!       sumd = SUM(sunpaths(n,1:ntraverse(n)))
!       sth1 = stheta*RadCrit/radii(0)
!       sume = dsin(theta_all(n) - dasin(sth1))*radii(0)/stheta
!       write(*,*)n,sumd,sume
!       do j = 1, nfinedivs(n)
!         sumd = SUM(sunpaths_fine(n,1:ntraverse_fine(n,j),j))
!         sth1 = sthetaf(j)*radii_fine(n,j)/radii(0)
!         sume = dsin(thetaf(j) - dasin(sth1))*radii(0)/sthetaf(j)
!         write(*,*)j,sumd,sume
!       enddo
!       pause
!      endif

!  Fix phi by using constancy of scatter angle


        if (stheta_boa.eq.0.0d0.or.salpha_boa.eq.0.0d0 ) then
           phi_all(n)     = phi_boa * dtr
        else
           if ( do_Normal ) then
              if ( doCrit .and. n .eq. NCrit ) then
                 salpha = dsin(AlphaCrit)
                 calpha = dcos(AlphaCrit)
              else
                 salpha = dsin(alpha(n))
                 calpha = dcos(alpha(n))
              endif
              cphi = (cosscat+vsign*calpha*ctheta)/stheta/salpha
              if ( cphi.gt.1.0d0)  cphi = 1.0d0
              if ( cphi.lt.-1.0d0) cphi = -1.0d0
              phi_all(n)     = dacos(cphi)
           endif
        endif

!  End layer loop

      enddo

!  TOA Sun angle sunpaths and PHI.
!    (No sunpaths if directly illuminated)

      CumAngle = alpha(nlayers) - alpha(0) ; Radstart = radii(0)
      call FindSun(Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_all(0),stheta,ctheta,DirectSun)
      if (.not.DirectSun ) then
          call FindSunPaths_T (Maxlayers,Pie,Radii(0),Radii,theta_all(0),stheta,1,sunpaths(0,:),ntraverse(0))
      endif
      if (stheta_boa.eq.0.0d0 ) then
         phi_all(0)     = phi_boa * dtr
      else
         cphi = (cosscat+vsign*calpha*ctheta)/stheta/salpha
         if ( cphi.gt.1.0d0)  cphi = 1.0d0 ; if ( cphi.lt.-1.0d0) cphi = -1.0d0
         phi_all(0)     = dacos(cphi)
      endif

!  Chapman factor calculations
!  ---------------------------

      Do_OverheadSun = stheta_boa.eq.0.0d0 
      if ( do_Chapman ) then
         do n = 1, nlayers
            call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii(n),Radii,theta_boa_R,stheta_boa,N,chapfacs(n,:))
            do k = 1, n
               chapfacs(n,k) = chapfacs(n,k)/(radii(k-1)-radii(k))
            enddo
         enddo
      endif

!  Finish

      return
end subroutine SD_incoming_sphergeom
!

subroutine RegularPS_sphergeom                                          &
       ( maxlayers, nlayers, heights, eradius, do_Chapman,              & ! Inputs
         alpha_boa, theta_boa, phi_boa, vsign, dtr,                     & ! Inputs
         Raycon, radii, alpha, sunpaths, ntraverse,                     & ! Outputs
         chapfacs, cosscat, term1, term2 )                                ! Outputs

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Regular PS choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometry
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

      implicit none

!  parameter argument

      integer, parameter :: fpk = selected_real_kind(15)

!  inputs

      integer  , intent(In)    :: maxlayers
      integer  , intent(In)    :: nlayers
      real(fpk), intent(InOut) :: alpha_boa, theta_boa, phi_boa
      logical  , intent(In)    :: do_Chapman
      real(fpk), intent(In)    :: dtr, eradius, heights (0:maxlayers), vsign

!  Los geometry

      real(fpk), intent(Out)  :: alpha         (0:maxlayers)
      real(fpk), intent(Out)  :: radii         (0:maxlayers)
      real(fpk), intent(Out)  :: Raycon

!  main outputs (geometry)

      integer  , intent(Out)  :: ntraverse  (0:maxlayers)
      real(fpk), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers)
      real(fpk), intent(Out)  :: chapfacs   (maxlayers,maxlayers)
      real(fpk), intent(Out)  :: cosscat, term1, term2

!  Local

      logical       :: Do_OverheadSun
      integer       :: n, k
      real(fpk)     :: alpha_boa_R, theta_boa_R
      real(fpk)     :: salpha_boa, calpha_boa
      real(fpk)     :: stheta_boa, ctheta_boa, cphi_boa

!  Initialise output

      radii = 0.0d0 ; alpha = 0.0d0 ; Raycon = 0.0d0
      ntraverse = 0
      sunpaths  = 0.0d0 ; chapfacs  = 0.0d0
      cosscat   = 0.0d0 ; term1 = 0.0d0 ; term2 = 0.0d0

!  BOA angles

      alpha_boa_R    = alpha_boa * DTR
      theta_boa_R    = theta_boa * DTR
      calpha_boa     = dcos(alpha_boa_R)
      salpha_boa     = dsin(alpha_boa_R)
      stheta_boa     = dsin(theta_boa_R)
      ctheta_boa     = dcos(theta_boa_R)
      cphi_boa       = dcos(phi_boa * dtr)

!  Nominal traverse paths for Full illumination

      do n = 1, nlayers
         ntraverse(n) = n
      enddo

!  Radii

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Set Alpha, ray constant, scattering angle

      alpha(1:nlayers) = alpha_boa_R
      Raycon           = salpha_boa * radii(nlayers)
      term1 = salpha_boa * stheta_boa * cphi_boa ; term2 = calpha_boa * ctheta_boa
      cosscat = - vsign * term2 + term1 

!  Sunpath/Chapman factor calculations

      Do_OverheadSun = stheta_boa.eq.0.0d0 
      do n = 1, nlayers
         call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii(n),Radii,theta_boa_R,stheta_boa,N,sunpaths(n,:))
         if ( do_Chapman ) then
            do k = 1, n
               chapfacs(n,k) = sunpaths(n,k)/(radii(k-1)-radii(k))
            enddo
         endif
      enddo

!  Finish

      return
end subroutine RegularPS_sphergeom

!  General Routines for Sun positioning

subroutine FindSun(Do_OverheadSun,Radius,SolarDirection,CumAngle,theta,stheta,ctheta,DirSun)

!  Find the solar anlge along the LOS path, for given radius and cumulative angle from BOA
!    SolarDirection is defined at BOA, with azimuth relative to the LOS direction.

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs

   logical   , Intent(In)    :: Do_OverheadSun
   real(fpk) , Intent(in)    :: Radius,SolarDirection(3),CumAngle

!  Outputs

   real(fpk) , Intent(out)   :: theta,stheta,ctheta
   logical   , Intent(InOut) :: DirSun

!  Local

   real(fpk) :: px(3),b

!  Calculation (overhead sun)

   if ( Do_OverheadSun ) then
      ctheta = dcos(CumAngle)
      DirSun = .true.
      stheta = dsin(CumAngle)
      theta  = CumAngle
      return
   endif

!  Calculation (General)

   px(1) = - Radius * dsin(CumAngle)
   px(2) = 0.0d0
   px(3) =   Radius * dcos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta = -b/Radius
   DirSun = (ctheta.ge.0.d0)
   stheta = dsqrt(1.0d0-ctheta*ctheta)
   theta  = dacos(ctheta)

!  Done

   return
end subroutine FindSun


subroutine FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,thstart,sthstart,N,sunpaths)

!  Sunpaths for the Direct-sun illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N
!  Special case = Overhead sun at BOA

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  inputs

   LOGICAL   , Intent(In)   :: Do_ZeroSunBOA
   INTEGER   , Intent(In)   :: maxlayers, N
   real(fpk) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(fpk) , Intent(In)   :: thstart,sthstart

!  Output

   real(fpk), Intent(InOut) :: Sunpaths(maxlayers)

!  Local

   integer    :: n1, k
   real(fpk)  :: sth0, th0, sth1, th1, ks1

!  Layer boundary upper

   N1 = N - 1

!  SBOA condition

   if ( Do_ZeroSunBOA ) then
      sunpaths(n) = radii(n1) - Radstart
      do k = n1, 1, -1
         sunpaths(n) = radii(k-1) - radii(k)
      enddo
      return
   endif

!  First layer

   sth0 = sthstart
   th0  = thstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = dasin(sth1)
   ks1  = th0-th1
   sunpaths(n) = dsin(ks1)*Radstart/sth1

!  Other layers to TOA

   sth0 = sth1
   th0  = th1
   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = dasin(sth1)
      ks1  = th0-th1
      sunpaths(k) = dsin(ks1)*radii(k)/sth1 
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_D

subroutine FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,thstart,sthstart,N,sunpaths,NT)

!  Sunpaths for the Tangent-height illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  inputs

   INTEGER   , Intent(In)   :: maxlayers, N
   real(fpk) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(fpk) , Intent(In)   :: thstart,sthstart,Pie

!  Output

   INTEGER   , Intent(InOut)   :: NT
   real(fpk), Intent(InOut)    :: Sunpaths(maxlayers)

!  Local

   logical    :: trawl
   integer    :: n1, k
   real(fpk)  :: sth0, th0, sth1, th1, ks1, tanr

!  Layer boundary upper

   N1 = N - 1

!  tangent height, Find which layer NT

   NT = N
   tanr = sthstart * Radstart
   k = n1 ; trawl = .true.
   do while (k.ge.n1.and.trawl)
      trawl = (radii(k).gt.tanr) ; k = k + 1
   enddo
   nt = k-1 !; write(*,*)n,nt

!  Distances for layers N and below to NT

   if ( nt.gt.n ) then
      th0  = pie - thstart ; sth0 = sthstart
      sth1 = sth0*Radstart/radii(n)
      th1  = dasin(sth1) ; ks1  = th0-th1
      sunpaths(n) = 2.0d0 * dsin(ks1)*Radstart/sth1
      sth0 = sth1 ; th0 = th1
      do k = n+1,nt-1
        sth1 = sth0*radii(k-1)/radii(k)
        th1  = dasin(sth1) ; ks1  = th0-th1
        sunpaths(k) = 2.0d0 * dsin(ks1)*radii(k)/sth0
        sth0 = sth1 ; th0 = th1
      enddo
      sth1 = 1.0d0 ; ks1 = 0.5*pie - th0
      sunpaths(nt) = 2.0d0 * dsin(ks1)*radii(nt-1)
   else if ( nt.eq.n ) then
      sunpaths(n) = - 2.0d0 * Radstart * dcos(thstart)
   endif

!  Rest of layer n up to the upper boundary

   th0 = pie - thstart ; sth0 = sthstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = dasin(sth1) ; ks1  = th0-th1
   sunpaths(n) = sunpaths(n) + dsin(ks1)*Radstart/sth1
   sth0 = sth1 ; th0 = th1

!  Trawl up from layers above n, to TOA

   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = dasin(sth1)
      ks1  = th0-th1
      sunpaths(k) = dsin(ks1)*radii(k)/sth1 
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_T

SUBROUTINE GAULEG_NG(X1,X2,X,W,N,NMAX)

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Input/Output

      INTEGER  , intent(in)  :: N,NMAX
      REAL(fpk), intent(in)  :: X1, X2
      REAL(fpk), intent(out) :: X(NMAX),W(NMAX)

      INTEGER     :: I, M, J
      REAL(fpk)   :: XM,XL,P1,P2,P3,PP,Z,Z1
      REAL(fpk), PARAMETER :: EPS = 3.0D-14

      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)

      DO I=1,M
            Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
            Z1 = 0.0d0
            DO WHILE (DABS(Z-Z1).GT.EPS)
                  P1=1.D0
                  P2=0.D0
                  DO J=1,N
                        P3=P2
                        P2=P1
                        P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
                  ENDDO
                  PP=N*(Z*P1-P2)/(Z*Z-1.D0)
                  Z1=Z
                  Z=Z1-P1/PP
            ENDDO
            X(I)=XM-XL*Z
            X(N+1-I)=XM+XL*Z
            W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
            W(N+1-I)=W(I)
      ENDDO
      RETURN
END SUBROUTINE GAULEG_NG


subroutine FindSunPath  ( x, xboa, rtop, raycon, sundir, sundist, theta0 )

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  I/O

   real(fpk), intent(in)  ::  x, xboa, rtop, raycon, sundir(3)
   real(fpk), intent(out) ::  sundist, theta0

!  Local

   real(fpk) :: xicum, sinx, rad, c0, s0, s1, c1

!  Subroutine for the quick calculation of sunpath from point X on View Path to point at Top of layer

   xicum = xboa - x
   sinx = dsin(x)
   rad   = Raycon / sinx
   c0 = sundir(1) * dsin(xicum) - sundir(3) * dcos(xicum)
   theta0 = - dacos(c0)
   s0 = dsqrt(1.0d0-c0*c0)
   s1 = s0 * rad / rtop
   c1 = dsqrt(1.0d0-s1*s1)
   sundist = -rtop * (s1*c0-s0*c1)/s0

!  finish

   return
end subroutine FindSunPath



SUBROUTINE LEGENDRE_NG ( STARTER, MAXMOMS, NMOMS, DF1, DF2, MU, SS_PLEG )

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  I/O

      LOGICAL  , intent(inout) :: STARTER
      INTEGER  , intent(in)    :: MAXMOMS, NMOMS
      REAL(fpk), intent(in)    :: MU
      REAL(fpk), intent(out)   :: SS_PLEG(0:MAXMOMS)
      REAL(fpk), intent(inout) :: DF1(MAXMOMS)
      REAL(fpk), intent(inout) :: DF2(MAXMOMS)

!  Local

      integer  :: L
     
!  Help arrays

      IF ( STARTER ) THEN
         DF1(1) = 0.0d0 ; DF2(1) = 0.0d0
         DO L = 2, NMOMS
            DF1(L) = DBLE(2*L-1) / DBLE(L)
            DF2(L) = DBLE(L-1)   / DBLE(L)
         ENDDO
         STARTER = .false.
      ENDIF

!  Legendre

      SS_PLEG(0) = 1.0d0
      SS_PLEG(1) = MU
      DO L = 2, NMOMS
         SS_PLEG(L) = DF1(L) * SS_PLEG(L-1) * MU - DF2(L) * SS_PLEG(L-2)
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LEGENDRE_NG

!  Finish

end module FirstOrder_geometry_Mk2_20jan12_m

