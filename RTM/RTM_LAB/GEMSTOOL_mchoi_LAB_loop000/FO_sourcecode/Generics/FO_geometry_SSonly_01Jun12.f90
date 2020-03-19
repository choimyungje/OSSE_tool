

module FO_geometry_SSonly_m

   use FO_geometry_Pool_m

private
public SS_Geometry_1, SS_Geometry_2, RegularPS_sphergeom, PlanePargeom

contains

subroutine SS_Geometry_1                                                 &
       ( maxlayers, maxfine, nlayers, nfinedivs, dtr, Pie, vsign,        & ! Input
         eradius, heights, alpha_boa, theta_boa, phi_boa, do_Chapman,    & ! Input
         do_LOSpaths, do_planpar, do_regular_ps, doCrit, Acrit, extinc,  & ! Input
         doNadir, Raycon, radii, alpha, cota,                            & ! Input/Output
         xfine, wfine, csqfine, cotfine, alphafine, radiifine,           & ! Input/Output
         NCrit, RadCrit, CotCrit, cosscat, chapfacs,                     & ! Output
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,             & ! Output
         fail, message, trace )                                            ! Output

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

      integer  , intent(in)     :: maxlayers, maxfine

!  Layer control. Finelayer divisions may be changed

      integer  , intent(in)     :: nlayers
      integer  , intent(inout)  :: nfinedivs(maxlayers)

!  Radius + heights

      real(ffp), intent(in)     :: eradius, heights (0:maxlayers)

!  input angles (Degrees), dtr = degrees-to-Radians. VSIGN = +1 (Up); -1(Down)

      real(ffp), intent(in)     :: dtr, Pie, vsign
      real(ffp), intent(InOut)  :: alpha_boa, theta_boa, phi_boa

!  Flags

      logical  , intent(in)     :: do_Chapman
      logical  , intent(in)     :: do_regular_ps
      logical  , intent(in)     :: do_planpar

!    LOS paths flag is new, 20 January 2012

      logical  , intent(in)     :: do_LOSpaths

!  Critical adjustment for cloud layers

      logical  , intent(inout)  :: doCrit
      real(ffp), intent(in)     :: extinc(maxlayers)
      real(ffp), intent(in)     :: Acrit

!  Input/Output Arguments (These may already be set if do_LOSpaths )
!  ======================

!  Flag for the Nadir case. Intent(inout), input if DO_LOSpaths set)

      logical  , intent(inout)  :: doNadir
  
!  Alphas,  Cotangents, Radii, Ray constant. Intent(inout), input if DO_LOSpaths set)

      real(ffp), intent(inout)  :: radii    (0:maxlayers)
      real(ffp), intent(inout)  :: Raycon
      real(ffp), intent(inout)  :: alpha    (0:maxlayers)
      real(ffp), intent(inout)  :: cota     (0:maxlayers)

!  LOS Quadratures for Enhanced PS

      real(ffp), intent(inout)  :: xfine    (maxlayers,maxfine)
      real(ffp), intent(inout)  :: wfine    (maxlayers,maxfine)
      real(ffp), intent(inout)  :: csqfine  (maxlayers,maxfine)
      real(ffp), intent(inout)  :: cotfine  (maxlayers,maxfine)

!  Fine layering output

      real(ffp), intent(inout)  :: alphafine (maxlayers,maxfine)
      real(ffp), intent(inout)  :: radiifine (maxlayers,maxfine)

!  Output arguments
!  ================

!  Critical layer

      integer  , intent(out)  :: Ncrit
      real(ffp), intent(out)  :: RadCrit, CotCrit

!  solar paths 

      integer  , Intent(out)  :: ntraverse  (0:maxlayers)
      real(ffp), Intent(out)  :: sunpaths   (0:maxlayers,maxlayers)
      integer  , Intent(out)  :: ntraverse_fine(maxlayers,maxfine)
      real(ffp), Intent(out)  :: sunpaths_fine (maxlayers,maxlayers,maxfine)
      real(ffp), Intent(out)  :: Chapfacs   (maxlayers,maxlayers)

!  Cosin scattering angle

      REAL(ffp), Intent(out)  :: cosscat

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message
      character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

      real(ffp)  :: Lospaths (maxlayers)

!  Other angles

      real(ffp)  :: theta_all  (0:maxlayers)
      real(ffp)  :: phi_all    (0:maxlayers)
      real(ffp)  :: cosa       (0:maxlayers)
      real(ffp)  :: sina       (0:maxlayers)

!  Critical values

      real(ffp)  :: AlphaCrit

!  Help variables

!      integer                :: n
      real(ffp), parameter   :: zero    = 0.0d0
      real(ffp)              :: term1, term2, cutoff

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   NCrit = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero
   cutoff = -log(ACrit)

!  check range of inputs
!  ---------------------

!  Cannot have Plane-parallel and Regular PS

   if ( do_planpar .and. do_regular_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Regular PS options'
      trace   = 'Initial Flag Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  VZA can be 0-90 degrees inclusive, but not outside this range

   if ( alpha_boa.gt.90.0d0.or.alpha_boa.lt.0.0d0 ) then
      message = 'Boa LOS angle outside range [0,90]); Check it!'
      trace   = 'Initial Angle Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  PHI is not limited to <= 180 degs. Also, not negative.
!     Old Code :     if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa

   if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
   if ( phi_boa.gt.360.0d0 ) phi_boa = 360.0d0 - phi_boa - 360.0d0

!  SZA can be 0-90 degrees inclusive, but not outside this range
!    For plane-parallel, 90 degrees is not allowed

   if ( do_planpar ) then
      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
         message = 'Plane-parallel: Boa SZA angle outside range [0,90)); Check it!'
         trace   = 'Initial Angle Check in SS_Geometry_1'
         fail    = .true. ;  return
      endif
   else
      if ( theta_boa.gt.90.0d0.or.theta_boa.lt.0.0d0 ) then
         message = 'Pseudo-spherical : Boa SZA angle outside range [0,90]); Check it!'
         trace   = 'Initial Angle Check in SS_Geometry_1'
         fail    = .true. ;  return
      endif
   endif

!  Plane-parallel, One routine only
!  --------------------------------

   if ( do_planpar ) then
      CALL PlanePargeom &
       ( maxlayers, nlayers, heights, do_Chapman,         & ! Inputs
         alpha_boa, theta_boa, phi_boa, vsign, dtr,       & ! Inputs
         alpha, sunpaths, ntraverse,                      & ! Outputs
         chapfacs, cosscat, term1, term2 )                  ! Outputs
      return
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
!  ===============================

!  Step 1; Initial LOS-path quantities
!  -----------------------------------

!    Given heights and BOA LOS angle, compute path angles and radii
!   Only need to do this if LOSpaths flag is not set.

   if ( .not. do_LOSpaths ) then
      CALL STD_outgoing_sphergeom_Initial                            &
       ( maxlayers, nlayers, heights, eradius, alpha_boa, dtr,       & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, sina, cosa, cota )   ! Output
   endif

!  Step 2, Adjust for Criticality
!  ------------------------------
   
!  Step 2a; Outgoing, Find Critical-layer adjustments (Optional)

   if ( doCrit) then
      CALL STD_outgoing_FindCritlayer                           &
       ( maxlayers, nlayers, zero, Acrit, Cutoff, doNadir,      & ! Inputs
         extinc, Lospaths, sina, cosa, radii, nfinedivs,        & ! Input
         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )      ! Outputs
      if ( Fail ) then
         trace = 'Error from STD_Outgoing_FindCritLayer in SS_Geometry_1' ; return
      endif
   endif

!  Step 2b; Incoming, Find Critical-layer adjustments (Optional)

   if ( doCrit) then
      call SD_incoming_FindCritLayer                            &
       ( maxlayers, nlayers, doNadir, zero, dtr, Acrit,         & ! Input
         cutoff, alpha, radii, extinc, Raycon, theta_boa,       & ! Inputs
         doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit, & ! Outputs
         fail, message )                                          ! Outputs
      if ( Fail ) then
         trace = 'Error from SD_incoming_FindCritLayer in SS_Geometry_1' ; return
      endif
   endif

!  Step 3. Set Quadratures
!  -----------------------

!  Step 3a. LOS fine-layer quadratures (Regular, non-adjusted, no Criticality)
!           Regular quadratures may be avoided if do_LOSpaths is set

!   write(*,*)'here 1',doCrit
   if ( .not. doCrit) then
      CALL STD_outgoing_sphergeom_Qbasic                    &
       ( maxlayers, maxfine, nlayers, nfinedivs,            & ! Input
         doNadir, radii, alpha, Raycon,                     & ! Input
         radiifine, alphafine, xfine, wfine,                & ! Output
         csqfine, cotfine, fail, message )                    ! Output
   endif
!   write(*,*)'here 2',doCrit

!  Step 3b. LOS fine-layer quadratures
!           Critical-layer adjustment of Qudarature done here.
!           Regular quadratures may be avoided if do_LOSpaths is set

   if ( doCrit) then
      CALL STD_outgoing_sphergeom_Qadjusted                 &
       ( maxlayers, maxfine, nlayers, nfinedivs,         & ! Input
         do_LOSpaths, doNadir, radii, alpha, Raycon,     & ! Input
         doCrit, Ncrit, AlphaCrit, RadCrit,              & ! Input
         radiifine, alphafine, xfine, wfine,             & ! Output
         csqfine, cotfine, fail, message )                 ! Output
      if ( Fail ) then
         trace = 'Error from STD_outgoing_sphergeom_Qadjusted in SS_Geometry_1' ; return
      endif
   endif
!   write(*,*)'here 3',doCrit,doNadir

!  Step 4. Solar path lengths
!  --------------------------

   CALL SD_incoming_sphergeom                                         &
       ( maxlayers, maxfine, nlayers, nfinedivs, do_Chapman, doNadir, & ! Input
         DoCrit, NCrit, alpha_boa, theta_boa, phi_boa, radii, alpha,  & ! Input
         vsign, dtr, Pie, RadCrit, AlphaCrit, radiifine, alphafine,   & ! Input
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,          & ! Output
         chapfacs, cosscat, theta_all, phi_all )                        ! Output

!  Finish

   return
end subroutine SS_Geometry_1

!

subroutine SS_Geometry_2                                                 &
       ( maxlayers, maxfine, nlayers, nfinedivs, dtr, Pie, vsign,        & ! Input
         eradius, heights, alpha_boa, theta_boa, phi_boa, do_Chapman,    & ! Input
         do_LOSpaths, do_planpar, do_regular_ps, doCrit, Acrit, extinc,  & ! Input
         doNadir, Raycon, radii, alpha, cota,                            & ! Input/Output
         xfine, wfine, csqfine, cotfine, alphafine, radiifine,           & ! Input/Output
         NCrit, RadCrit, CotCrit,                                        & ! Output
         cosscat_up, chapfacs_up, cosscat_dn, chapfacs_dn,               & ! Output
         sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up, & ! Output
         sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn, & ! Output
         fail, message, trace )                                            ! Output

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

      integer  , intent(in)     :: maxlayers, maxfine

!  Layer control. Finelayer divisions may be changed

      integer  , intent(in)     :: nlayers
      integer  , intent(inout)  :: nfinedivs(maxlayers)

!  Radius + heights

      real(ffp), intent(in)     :: eradius, heights (0:maxlayers)

!  input angles (Degrees), dtr = degrees-to-Radians

      real(ffp), intent(in)     :: dtr, Pie
      real(ffp), intent(InOut)  :: alpha_boa, theta_boa, phi_boa

!  Flags

      logical  , intent(in)     :: do_Chapman
      logical  , intent(in)     :: do_regular_ps
      logical  , intent(in)     :: do_planpar

!    LOS paths flag is new, 20 January 2012

      logical  , intent(in)     :: do_LOSpaths

!  Critical adjustment for cloud layers

      logical  , intent(inout)  :: doCrit
      real(ffp), intent(in)     :: extinc(maxlayers)
      real(ffp), intent(in)     :: Acrit

!  Input/Output Arguments (These may already be set if do_LOSpaths )
!  ======================

!  Flag for the Nadir case. Intent(inout), input if DO_LOSpaths set)

      logical  , intent(inout)  :: doNadir
  
!  Alphas,  Cotangents, Radii, Ray constant. Intent(inout), input if DO_LOSpaths set)

      real(ffp), intent(inout)  :: radii    (0:maxlayers)
      real(ffp), intent(inout)  :: Raycon
      real(ffp), intent(inout)  :: alpha    (0:maxlayers)
      real(ffp), intent(inout)  :: cota     (0:maxlayers)

!  LOS Quadratures for Enhanced PS

      real(ffp), intent(inout)  :: xfine    (maxlayers,maxfine)
      real(ffp), intent(inout)  :: wfine    (maxlayers,maxfine)
      real(ffp), intent(inout)  :: csqfine  (maxlayers,maxfine)
      real(ffp), intent(inout)  :: cotfine  (maxlayers,maxfine)

!  Fine layering output

      real(ffp), intent(inout)  :: alphafine (maxlayers,maxfine)
      real(ffp), intent(inout)  :: radiifine (maxlayers,maxfine)

!  Output arguments
!  ================

!  Critical layer

      integer  , intent(out)  :: Ncrit
      real(ffp), intent(out)  :: RadCrit, CotCrit

!  solar paths upwelling

      integer  , Intent(out)  :: ntraverse_up  (0:maxlayers)
      real(ffp), Intent(out)  :: sunpaths_up   (0:maxlayers,maxlayers)
      integer  , Intent(out)  :: ntraverse_fine_up(maxlayers,maxfine)
      real(ffp), Intent(out)  :: sunpaths_fine_up (maxlayers,maxlayers,maxfine)
      real(ffp), Intent(out)  :: Chapfacs_up  (maxlayers,maxlayers)

!  solar paths downwelling

      integer  , Intent(out)  :: ntraverse_dn  (0:maxlayers)
      real(ffp), Intent(out)  :: sunpaths_dn   (0:maxlayers,maxlayers)
      integer  , Intent(out)  :: ntraverse_fine_dn(maxlayers,maxfine)
      real(ffp), Intent(out)  :: sunpaths_fine_dn (maxlayers,maxlayers,maxfine)
      real(ffp), Intent(out)  :: Chapfacs_dn  (maxlayers,maxlayers)

!  Cosine scattering angles

      REAL(ffp), Intent(out)  :: cosscat_up, cosscat_dn

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message
      character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

      real(ffp)  :: Lospaths (maxlayers)

!  Other angles

      real(ffp)  :: theta_all  (0:maxlayers)
      real(ffp)  :: phi_all    (0:maxlayers)
      real(ffp)  :: cosa       (0:maxlayers)
      real(ffp)  :: sina       (0:maxlayers)

!  Critical values

      real(ffp)  :: AlphaCrit

!  Help variables
!    VSIGN = +1 for Upwelling, -1 for Downwelling

!      integer                :: n
      real(ffp), parameter   :: zero    = 0.0d0
      real(ffp)              :: vsign, term1, term2, cutoff

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   NCrit = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero
   cutoff = -log(ACrit)

!  check range of inputs
!  ---------------------

!  Cannot have Plane-parallel and Regular PS

   if ( do_planpar .and. do_regular_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Regular PS options'
      trace   = 'Initial Flag Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  VZA can be 0-90 degrees inclusive, but not outside this range

   if ( alpha_boa.gt.90.0d0.or.alpha_boa.lt.0.0d0 ) then
      message = 'Boa LOS angle outside range [0,90]); Check it!'
      trace   = 'Initial Angle Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  PHI is not limited to <= 180 degs. Also, not negative.
!     Old Code :     if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa

   if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
   if ( phi_boa.gt.360.0d0 ) phi_boa = 360.0d0 - phi_boa - 360.0d0

!  SZA can be 0-90 degrees inclusive, but not outside this range
!    For plane-parallel, 90 degrees is not allowed

   if ( do_planpar ) then
      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
         message = 'Plane-parallel: Boa SZA angle outside range [0,90)); Check it!'
         trace   = 'Initial Angle Check in SS_Geometry_1'
         fail    = .true. ;  return
      endif
   else
      if ( theta_boa.gt.90.0d0.or.theta_boa.lt.0.0d0 ) then
         message = 'Pseudo-spherical : Boa SZA angle outside range [0,90]); Check it!'
         trace   = 'Initial Angle Check in SS_Geometry_1'
         fail    = .true. ;  return
      endif
   endif

!  Plane-parallel, One routine only
!  --------------------------------

!  Single call to upwelling (VSIGN = +1)
!  Downwelling only requires new COSSCAT, rest is copied

   if ( do_planpar ) then
      vsign = +1.0d0
      CALL PlanePargeom &
       ( maxlayers, nlayers, heights, do_Chapman,         & ! Inputs
         alpha_boa, theta_boa, phi_boa, vsign, dtr,       & ! Inputs
         alpha, sunpaths_up, ntraverse_up,                & ! Outputs
         chapfacs_up, cosscat_up, term1, term2 )            ! Outputs
      vsign = -1.0d0  ;  cosscat_dn = -vsign * term1 + term2
      ntraverse_dn = ntraverse_up
      sunpaths_dn  = sunpaths_up
      chapfacs_dn  = chapfacs_up
      return
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
!  ===============================

!  Step 1; Basic LOS quantities
!  ----------------------------

!    Given heights and BOA LOS angle, compute path angles and radii
!   Only need to do this if LOSpaths flag is not set.

   if ( .not. do_LOSpaths ) then
      CALL STD_outgoing_sphergeom_Initial                            &
       ( maxlayers, nlayers, heights, eradius, alpha_boa, dtr,       & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, sina, cosa, cota )   ! Output
   endif

!  Step 2a; Outgoing, Find Critical-layer adjustments (Optional)

   if ( doCrit) then
      CALL STD_outgoing_FindCritlayer                           &
       ( maxlayers, nlayers, zero, Acrit, Cutoff, doNadir,      & ! Inputs
         extinc, Lospaths, sina, cosa, radii, nfinedivs,        & ! Input
         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )      ! Outputs
      if ( Fail ) then
         trace = 'Error from STD_Outgoing_FindCritLayer in SS_Geometry_1' ; return
      endif
   endif

!  Step 2b; Incoming, Find Critical-layer adjustments (Optional)

   if ( doCrit) then
      call SD_incoming_FindCritLayer                            &
       ( maxlayers, nlayers, doNadir, zero, dtr, Acrit,         & ! Input
         cutoff, alpha, radii, extinc, Raycon, theta_boa,       & ! Inputs
         doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit, & ! Outputs
         fail, message )                                          ! Outputs
      if ( Fail ) then
         trace = 'Error from SD_incoming_FindCritLayer in SS_Geometry_2' ; return
      endif
   endif

!  Step 3. Set Quadratures
!  -----------------------

!  Step 3a. LOS fine-layer quadratures (Regular, non-adjusted, no Criticality)
!           Regular quadratures may be avoided if do_LOSpaths is set

   if ( .not. doCrit) then
      CALL STD_outgoing_sphergeom_Qbasic                    &
       ( maxlayers, maxfine, nlayers, nfinedivs,            & ! Input
         doNadir, radii, alpha, Raycon,                     & ! Input
         radiifine, alphafine, xfine, wfine,                & ! Output
         csqfine, cotfine, fail, message )                    ! Output
   endif

!  Step 3b. LOS fine-layer quadratures
!           Critical-layer adjustment of Qudarature done here.
!           Regular quadratures may be avoided if do_LOSpaths is set

   if ( doCrit) then
      CALL STD_outgoing_sphergeom_Qadjusted                 &
       ( maxlayers, maxfine, nlayers, nfinedivs,         & ! Input
         do_LOSpaths, doNadir, radii, alpha, Raycon,     & ! Input
         doCrit, Ncrit, AlphaCrit, RadCrit,              & ! Input
         radiifine, alphafine, xfine, wfine,             & ! Output
         csqfine, cotfine, fail, message )                 ! Output
      if ( Fail ) then
         trace = 'Error from STD_outgoing_sphergeom_Qadjusted in SS_Geometry_2' ; return
      endif
   endif

!  Step 4. Solar path lengths
!  --------------------------

!  Step 4a. Solar path lengths, Upwelling

   vsign = + 1.0d0
   CALL SD_incoming_sphergeom                                            &
       ( maxlayers, maxfine, nlayers, nfinedivs, do_Chapman, doNadir,    & ! Input
         DoCrit, NCrit, alpha_boa, theta_boa, phi_boa, radii, alpha,     & ! Input
         vsign, dtr, Pie, RadCrit, AlphaCrit, radiifine, alphafine,      & ! Input
         sunpaths_up, ntraverse_up, sunpaths_fine_up, ntraverse_fine_up, & ! Output
         chapfacs_up, cosscat_up, theta_all, phi_all )                     ! Output

!  Step 4b. Solar path lengths, Downwelling

   vsign = - 1.0d0
   CALL SD_incoming_sphergeom                                            &
       ( maxlayers, maxfine, nlayers, nfinedivs, do_Chapman, doNadir,    & ! Input
         DoCrit, NCrit, alpha_boa, theta_boa, phi_boa, radii, alpha,     & ! Input
         vsign, dtr, Pie, RadCrit, AlphaCrit, radiifine, alphafine,      & ! Input
         sunpaths_dn, ntraverse_dn, sunpaths_fine_dn, ntraverse_fine_dn, & ! Output
         chapfacs_dn, cosscat_dn, theta_all, phi_all )                     ! Output

!  Finish

   return
end subroutine SS_Geometry_2

!

subroutine Planepargeom                               &
       ( maxlayers, nlayers, heights, do_Chapman,     & ! Inputs
         alpha_boa, theta_boa, phi_boa, vsign, dtr,   & ! Inputs
         alpha, sunpaths, ntraverse,                  & ! Outputs
         chapfacs, cosscat, term1, term2 )              ! Outputs

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Regular PS choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometry
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

      implicit none

!  parameter argument

      integer, parameter :: ffp = selected_real_kind(15)

!  inputs

      integer  , intent(In)    :: maxlayers
      integer  , intent(In)    :: nlayers
      real(ffp), intent(InOut) :: alpha_boa, theta_boa, phi_boa
      logical  , intent(In)    :: do_Chapman
      real(ffp), intent(In)    :: dtr, heights (0:maxlayers), vsign

!  Los geometry

      real(ffp), intent(Out)  :: alpha         (0:maxlayers)

!  main outputs (geometry)

      integer  , intent(Out)  :: ntraverse  (0:maxlayers)
      real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers)
      real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers)
      real(ffp), intent(Out)  :: cosscat, term1, term2

!  Local

      logical       :: Do_OverheadSun
      integer       :: n
      real(ffp)     :: alpha_boa_R, theta_boa_R
      real(ffp)     :: salpha_boa, calpha_boa
      real(ffp)     :: stheta_boa, ctheta_boa, utheta_boa, cphi_boa, diffhts(maxlayers)

!  Initialise output

      alpha = 0.0d0
      ntraverse = 0
      sunpaths  = 0.0d0 ; chapfacs  = 0.0d0
      cosscat   = 0.0d0 ; term1 = 0.0d0 ; term2 = 0.0d0

!  BOA angles

      alpha_boa_R    = alpha_boa * DTR
      if ( alpha_boa.eq.90.0d0 ) then
         calpha_boa     = 0.0d0
         salpha_boa     = 1.0d0
      else
         calpha_boa     = dcos(alpha_boa_R)
         salpha_boa     = dsin(alpha_boa_R)
      endif

      theta_boa_R    = theta_boa * DTR
      stheta_boa     = dsin(theta_boa_R)
      ctheta_boa     = dcos(theta_boa_R)

      cphi_boa       = dcos(phi_boa * dtr)

!  Nominal traverse paths for Full illumination. Difference heights

      do n = 1, nlayers
         diffhts(n) = heights(n-1) - heights(n)
         ntraverse(n) = n
      enddo

!  Overhead Sun

      Do_OverheadSun = (theta_boa.eq.0.0d0) 

!  Set Alpha, scattering angle

      alpha(1:nlayers) = alpha_boa_R
      if ( Do_OverheadSun ) then
         term1 = 0.0d0
         term2 = calpha_boa
         cosscat = - vsign * term2 ; if (term2.eq.0.0d0) cosscat = term2
      else
         term1 = salpha_boa * stheta_boa * cphi_boa
         term2 = calpha_boa * ctheta_boa
         cosscat = - vsign * term2 + term1 
      endif

!  Sunpath/Chapman factor calculations

      utheta_boa     = 1.0d0 / ctheta_boa
      do n = 1, nlayers
         sunpaths(n,1:n) = diffhts(1:n) * utheta_boa
         if ( do_Chapman ) chapfacs(n,1:n) = utheta_boa
      enddo

!  Finish

      return
end subroutine Planepargeom

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

      integer, parameter :: ffp = selected_real_kind(15)

!  inputs

      integer  , intent(In)    :: maxlayers
      integer  , intent(In)    :: nlayers
      real(ffp), intent(InOut) :: alpha_boa, theta_boa, phi_boa
      logical  , intent(In)    :: do_Chapman
      real(ffp), intent(In)    :: dtr, eradius, heights (0:maxlayers), vsign

!  Los geometry

      real(ffp), intent(Out)  :: alpha         (0:maxlayers)
      real(ffp), intent(Out)  :: radii         (0:maxlayers)
      real(ffp), intent(Out)  :: Raycon

!  main outputs (geometry)

      integer  , intent(Out)  :: ntraverse  (0:maxlayers)
      real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers)
      real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers)
      real(ffp), intent(Out)  :: cosscat, term1, term2

!  Local

      logical       :: Do_OverheadSun
      integer       :: n, k
      real(ffp)     :: alpha_boa_R, theta_boa_R
      real(ffp)     :: salpha_boa, calpha_boa
      real(ffp)     :: stheta_boa, ctheta_boa, cphi_boa, sunpaths_local(maxlayers)

!  Initialise output

      radii = 0.0d0 ; alpha = 0.0d0 ; Raycon = 0.0d0
      ntraverse = 0
      sunpaths  = 0.0d0 ; chapfacs  = 0.0d0
      cosscat   = 0.0d0 ; term1 = 0.0d0 ; term2 = 0.0d0

!  BOA angles

      alpha_boa_R    = alpha_boa * DTR
      if ( alpha_boa.eq.90.0d0 ) then
         calpha_boa     = 0.0d0
         salpha_boa     = 1.0d0
      else
         calpha_boa     = dcos(alpha_boa_R)
         salpha_boa     = dsin(alpha_boa_R)
      endif

      theta_boa_R    = theta_boa * DTR
      if ( theta_boa.eq.90.0d0 ) then
         ctheta_boa     = 0.0d0
         stheta_boa     = 1.0d0
      else
         stheta_boa     = dsin(theta_boa_R)
         ctheta_boa     = dcos(theta_boa_R)
      endif

      cphi_boa       = dcos(phi_boa * dtr)

!  Nominal traverse paths for Full illumination

      do n = 1, nlayers
         ntraverse(n) = n
      enddo

!  Radii

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Overhead Sun

      Do_OverheadSun = (theta_boa.eq.0.0d0) 

!  Set Alpha, ray constant, scattering angle

      alpha(1:nlayers) = alpha_boa_R
      Raycon           = salpha_boa * radii(nlayers)
      if ( Do_OverheadSun ) then
         term1 = 0.0d0
         term2 = calpha_boa
         cosscat = - vsign * term2 ; if (term2.eq.0.0d0) cosscat = term2
      else
         term1 = salpha_boa * stheta_boa * cphi_boa
         term2 = calpha_boa * ctheta_boa
         cosscat = - vsign * term2 + term1
      endif

!  Sunpath/Chapman factor calculations

!mick fix 4/12/12 - adjusted dimension of array "sunpaths" assignments to be
!                   compatible with subroutine "FindSunPaths_D" computations
      !sunpaths(0,1:maxlayers) = 0.0d0
      sunpaths(0,1:nlayers) = 0.0d0
      do n = 1, nlayers
         call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii(n),Radii,&
           theta_boa_R,stheta_boa,N,sunpaths_local)
         !sunpaths(n,1:maxlayers) = sunpaths_local(1:maxlayers)
         sunpaths(n,1:n) = sunpaths_local(1:n)
         if ( do_Chapman ) then
            do k = 1, n
               chapfacs(n,k) = sunpaths(n,k)/(radii(k-1)-radii(k))
            enddo
         endif
      enddo

!  Finish

      return
end subroutine RegularPS_sphergeom


!  Finish

end module FO_geometry_SSonly_m

