

module FirstOrder_geometry_Mk3_m

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

!

subroutine STD_outgoing_sphergeom_Initial                           &
       ( maxlayers, nlayers, heights, eradius, alpha_boa, dtr,      & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, sina, cosa, cota )  ! Output

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!  This routine: Initial LOS path setup

!    starting inputs are - BOA value of VZA (alpha_boa), in degrees
!                        - height grid, earth radius

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

      integer  , intent(in)   :: maxlayers

!  Layer control

      integer  , intent(in)   :: nlayers
      real(ffp), intent(in)   :: eradius, heights (0:maxlayers)

!  input angle

      real(ffp), intent(in)   :: alpha_boa, dtr

!  Flag for the Nadir case

      logical  , intent(out)  :: doNadir
  
!  Alphas, Radii, Ray constant, Lospaths

      real(ffp), intent(out)  :: radii    (0:maxlayers)
      real(ffp), intent(out)  :: Raycon
      real(ffp), intent(out)  :: Lospaths (maxlayers)
      real(ffp), intent(out)  :: alpha    (0:maxlayers)
      real(ffp), intent(out)  :: cosa     (0:maxlayers)
      real(ffp), intent(out)  :: sina     (0:maxlayers)
      real(ffp), intent(out)  :: cota     (0:maxlayers)

!  Local
!  -----

      integer         :: n, n1
      real(ffp)       :: salpha_boa, difh, alpha_boa_R
      real(ffp)       :: calpha, calpha1
      real(ffp), parameter :: zero = 0.0_ffp

!  Zero output

      Lospaths = zero ; cota = zero ; cosa = zero ; sina = zero ; alpha = zero ; Raycon = zero

!  Radii

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Special case

      doNadir = .false.
      if ( alpha_boa.eq.zero ) doNadir = .true.

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
      if ( alpha_boa .eq. 90.0_ffp ) then
         salpha_boa     = 1.0_ffp
         calpha1        = 0.0_ffp
      else
         salpha_boa     = dsin(alpha_boa_R)
         calpha1        = dcos(alpha_boa_R)
      endif

      cosa(nlayers)  = calpha1
      sina(nlayers)  = salpha_boa
      cota(nlayers)  = calpha1 / salpha_boa
      alpha(nlayers) = alpha_boa_R

!  Ray constant

      Raycon = salpha_boa * radii(nlayers)

!  Whole layer values

      do n = nlayers - 1, 0, -1
         n1 = n + 1
         sina(n) = Raycon / radii(n) ; alpha(n) = dasin(sina(n))
         calpha  = dcos(alpha(n)) ; cosa(n) = calpha 
         cota(n) = cosa(n) / sina(n)
         Lospaths(n1) = radii(n)*calpha - radii(n1)*calpha1
         calpha1 = calpha
      enddo

!  Finish

      return
end subroutine STD_outgoing_sphergeom_initial

!

SUBROUTINE STD_outgoing_FindCritlayer                           &
       ( maxlayers, nlayers, zero, Acrit, Cutoff, doNadir,      & ! Inputs
         extinc, Lospaths, sina, cosa, radii, nfinedivs,        & ! Input
         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )      ! Outputs
!       ( maxlayers, nlayers, zero, dtr, Acrit, Cutoff, doNadir, & ! Inputs
!         extinc, Lospaths, alpha, sina, cosa, radii, nfinedivs, & ! Input
!         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )      ! Outputs

!  Purpose: Given a list of Maximum extinctions and solar
!     Then find Critical layer (NCrit) and point where LOS attenuation wipe-out (Acrit) is achieved
!     Then find the LOS angle and Radius (AlphaCrit,RadCrit) for this Critical Point

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Dimensioning

   integer, intent(in) :: maxlayers

!  Layer control

   integer, intent(in) :: nlayers

!  Attenuation and other parameters

   real(ffp), intent(in)  :: zero, Acrit, Cutoff !, dtr

!  Special case, Nadir viewing

   logical, intent(in)  :: doNadir

!  View angles and Radii at layer boundaries

!   real(ffp), intent(in)  :: alpha(0:maxlayers)
   real(ffp), intent(in)  :: sina (0:maxlayers)
   real(ffp), intent(in)  :: cosa (0:maxlayers)
   real(ffp), intent(in)  :: radii(0:maxlayers)

!  Extinctions

   real(ffp), intent(in)  :: Lospaths(maxlayers)
   real(ffp), intent(in)  :: extinc(maxlayers)

!  Modified inputs
!  ---------------

!  Number of Fine divisions

   integer, intent(inout) :: nfinedivs(maxlayers)

!  outputs
!  -------

!  Critical layer, Number of Fine divisions for this layer

   integer  , intent(out)  :: Ncrit

!  Critical angle and radius and cotangent

   real(ffp), intent(out)  :: AlphaCrit
   real(ffp), intent(out)  :: RadCrit, CotCrit

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

!  Other variables

   logical    ::  trawl
   integer    ::  N, N1, ntrans
   real(ffp)  ::  opdep, cumtrans_0, cumtrans_1,trans,transcrit,dcrit,TanCrit

!  Initialize

   Ncrit     = 0
   AlphaCrit = zero
   RadCrit   = zero
   CotCrit   = zero

   fail    = .false.
   message = ' '

!  Set trawl
!   Tested March 17th

   trawl = .true. ; n = 0 ; cumtrans_0 = 1.0d0
   do while (trawl .and.n.lt.nlayers) 
      n = n + 1
      opdep = Lospaths(n) * extinc(n) ; trans = zero
      if ( opdep .lt. cutoff ) trans = exp ( - opdep )
      cumtrans_1 = cumtrans_0 * trans
      if ( cumtrans_1 .gt. Acrit ) then
         ntrans = int(-log10(trans) + 1)
         nfinedivs(n) = max(nfinedivs(n),ntrans)
         cumtrans_0 = cumtrans_1
      else
         NCrit = n ; trawl = .false.
         transcrit    = Acrit / cumtrans_0
         ntrans = int(-log10(transcrit) + 1)
         nfinedivs(n) = max(nfinedivs(n),ntrans)
         dcrit        = - log(transcrit) / extinc(n)
         if ( doNadir ) then
            Radcrit    = radii(n-1) - dcrit      
         else
            n1 = n-1 ; TanCrit = radii(n1)*sina(n1)/(radii(n1)*cosa(n1)-dcrit)
            Cotcrit    = 1.0_ffp / Tancrit
            alphacrit  = atan( TanCrit)
            radcrit    = sina(n) * radii(n) / sin(alphacrit)  
         endif
      endif
   enddo

!  Zero the rest

   if ( NCrit .ne. 0 ) nfinedivs(NCrit+1:nlayers) = 0

!  Finish

   return
end subroutine STD_outgoing_FindCritlayer


subroutine SD_incoming_FindCritLayer                            &
       ( maxlayers, nlayers, doNadir, zero, dtr, Acrit,         & ! Input
         cutoff, alpha, radii, extinc, Raycon, theta_boa,       & ! Inputs
         doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit, & ! Outputs
         fail, message )                                          ! Outputs

!  Purpose: Given a list of Maximum extinctions and solar angle at BOA
!           Then find Critical layer (NCrit) and point where TOA attenuation wipe-out (Acrit) is achieved
!           Then find the LOS angle and Radius (AlphaCrit,RadCrit) for this Critical Point
!           Nadir case, Alpha = 0.0, find only the radius (RadCrit)

!  Find the Critical Radius (or angle) in layer Ncrit_i, Use Bisection 
!      based on the Function F(x) = opdep(x) - Crit_opdep

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Dimensioning

   integer, intent(in) :: maxlayers

!  Layer control

   integer, intent(in) :: nlayers

!  Special case, Nadir viewing

   logical, intent(in)  :: doNadir

!  View angles and Radii at layer boundaries, Ray constant

   real(ffp), intent(in)  :: alpha(0:maxlayers)
   real(ffp), intent(in)  :: radii(0:maxlayers)
   real(ffp), intent(in)  :: Raycon

!  Extinctions

   real(ffp), intent(in)  :: extinc(maxlayers)

!  Solar control and other parameters

   real(ffp), intent(in)  :: Acrit, theta_boa, dtr, zero, cutoff

!  Modified inputs (outputs)
!  -------------------------

!  Overall control (May be switched off if Critical test is negative)

   logical  , intent(inout)  :: doCrit

!  Critical layer

   integer  , intent(inout)  :: Ncrit

!  Number of Fine divisions. This is updated according to Criticality

   integer  , intent(inout)  :: nfinedivs(maxlayers)

!  Critical angle and radius and cotangent

   real(ffp), intent(inout)  :: AlphaCrit
   real(ffp), intent(inout)  :: RadCrit, CotCrit

!  Outputs
!  -------

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

!  Bisection accuracy (lower for the Nadir case, as using distances)

   real(ffp), parameter  :: BisectionAccuracy_Nadir   = 1.0d-6
   real(ffp), parameter  :: BisectionAccuracy_General = 1.0d-12
   integer  , parameter  :: jmax = 100

!  Other variables

   logical    ::  Finding, trawl, do_ZeroSunBOA
   integer    ::  j, n, ntrans, ncrit_i, k
   real(ffp)  ::  s0, x1, x2, xmid, rtbis, dx, f, fmid, suncon, radii_x, sx, accuracy
   real(ffp)  ::  dist, theta_0, theta_1, stheta_0, stheta_1, ground_R, sunpaths(maxlayers)
   real(ffp)  ::  opdep, trans, atten_0, atten_1, transcrit, theta_n, stheta_n, theta_boa_R

!  Initialize

   fail    = .false.
   message = ' '

!  Initial setups

   atten_0 = 1.0d0 ; theta_boa_R = theta_boa * dtr
   if ( doNadir ) then
      s0 = sin(theta_boa_R) ; ground_R = zero ; accuracy = BisectionAccuracy_Nadir
   else
      s0 = zero ; ground_R = alpha(nlayers) + theta_boa_R ; accuracy = BisectionAccuracy_General
   endif

!  Trawl through layers until Critical layer is reached. Nfinedivs is updated.
!    Only go down to Initial (LOS-path) Critical layer 
!       Condition on ZerosunBOA changed 27 March 2012

   NCrit_i = 0 ; trawl = .true. ; n = 0
   do while ( trawl .and.( (Ncrit.eq.0.and.n.lt.nlayers).or.(NCrit.ne.0.and.n.lt.NCrit) ) )
      n = n + 1
      do_ZeroSunBOA = (n.eq.nlayers.and.theta_boa.eq.0.0_ffp).or.(doNadir.and.theta_boa.eq.0.0_ffp)
      if ( doNadir ) then
         theta_n = theta_boa_R ; stheta_n = s0
      else
         theta_n = ground_R - alpha(n) ; stheta_n = sin(theta_n)
      endif
      call FindSunPaths_D (do_ZeroSunBOA,Maxlayers,radii(n),Radii,&
        theta_n,stheta_n,n,sunpaths)
      opdep = sum(extinc(1:n)*sunpaths(1:n)) ; trans = zero
      atten_1 = 0.0d0; if ( opdep .lt. cutoff ) atten_1 = exp ( - opdep )
      if ( atten_1 .gt. Acrit ) then
         trans = atten_1 / atten_0
         ntrans = int(-log10(trans) + 1)
         nfinedivs(n) = max(nfinedivs(n),ntrans)
         atten_0 = atten_1
      else
         NCrit_i = n ; trawl = .false.
         transcrit    = Acrit / atten_0
         ntrans = int(-log10(transcrit) + 1)
         nfinedivs(n) = max(nfinedivs(n),ntrans)
      endif
   enddo

!  Nothing to do if No criticality (previous Critical values are unchanged)

   if ( trawl .and. NCrit_i.eq. 0 ) then
      if ( trawl .and. NCrit  .eq. 0 ) DoCrit = .false.
      return
   endif

!  Bisection: set Highest/Lowest value of Function (layer bottom/top). 

   if ( doNadir ) then
      x1 = zero             ; x2 = radii(NCrit_i-1) - radii(NCrit_i)
   else
      x1 = alpha(NCrit_i-1) ; x2 = alpha(NCrit_i) 
   endif
   F     = log(Atten_0) - cutoff
   FMID  = opdep        - cutoff

!  Bisection: Check bracketing, if OK, perform bisection

   IF(F*FMID.GE.0.) then
       fail = .true. ; message = 'Root must be bracketed for bisection.' ; return
   ENDIF
   IF(F.LT.0.)THEN
      RTBIS=X1 ; DX=X2-X1
   ELSE
      RTBIS=X2 ; DX=X1-X2
   ENDIF

!  Bisection: Iterate to find the answer

   Finding = .true. ; J = 0
   DO While (Finding .and. j .lt. JMAX)
      J = J + 1 ; dx = 0.5d0 * dx ; XMID = RTBIS + DX
      if ( doNadir ) then
         theta_0 = theta_boa_R ; stheta_0 = s0
         radii_x = radii(NCrit_i-1) - xmid 
      else
         theta_0 = ground_R - xmid ; stheta_0 = sin(theta_0)
         radii_x = Raycon / sin(xmid)
      endif
      suncon = radii_x * stheta_0
      stheta_1 = suncon / radii(NCrit_i-1) ;  theta_1 = asin(stheta_1)
      dist = radii(NCrit_i-1) * sin(theta_0-theta_1) / stheta_0
      opdep = dist * extinc(NCrit_i)
      theta_0 = theta_1 ; stheta_1 = stheta_0
      do k = n - 1, 1, -1
         stheta_1 = suncon / radii(k-1) ; theta_1 = asin(stheta_1)
         dist = radii(k-1) * sin(theta_0-theta_1) / stheta_0
         opdep = opdep + dist * extinc(k)
         theta_0 = theta_1 ; stheta_0 = stheta_1
      enddo
      fmid = opdep - cutoff
      IF ( FMID.LE.0.0d0 ) RTBIS = XMID
      IF(ABS(DX).LT.Accuracy .OR. FMID.EQ.0.) Finding = .false.
   ENDDO

!  Exception (too many bisections)

   if ( Finding ) Then
      fail = .true.
      message = 'Too many Bisections (540); Root not found'
      return
   endif

!  Set final output if successful

   if ( doNadir ) then
      RadCrit   =  radii(NCrit_i-1) - RTBIS
   else
      AlphaCrit = RTBIS ;  SX = dsin(AlphaCrit)
      RadCrit   = Raycon / SX
      CotCrit   = dsqrt(1.0d0-SX*SX)/SX
   endif
   NCrit     = NCrit_i
   nfinedivs(NCrit+1:nlayers) = 0

!  Finish

   return
end subroutine SD_incoming_FindCritLayer

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

   integer, parameter :: ffp = selected_real_kind(15)

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

      real(ffp), intent(in)  :: alpha      (0:maxlayers)
      real(ffp), intent(in)  :: radii      (0:maxlayers)
      real(ffp), intent(in)  :: Raycon

!  Outputs
!  =======

!  Fine layering

      real(ffp), intent(out)  :: alpha_fine (maxlayers,maxfine)
      real(ffp), intent(out)  :: radii_fine (maxlayers,maxfine)

!  Quadratures

      real(ffp), intent(out)  :: xfine    (maxlayers,maxfine)
      real(ffp), intent(out)  :: wfine    (maxlayers,maxfine)

!  Local geoemetry arrays

      real(ffp), intent(out)  :: csqfine  (maxlayers,maxfine)
      real(ffp), intent(out)  :: cotfine  (maxlayers,maxfine)

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message

!  Local
!  -----

      integer, parameter :: Localmaxfine = 20
      integer            :: n, n1, j, nfine
      real(ffp)          :: difh, csfine
      real(ffp)          :: tfine(Localmaxfine), afine(Localmaxfine)

      real(ffp), parameter :: zero = 0.0d0

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

subroutine STD_outgoing_sphergeom_Qadjusted              &
       ( maxlayers, maxfine, nlayers, nfinedivs,         & ! Input
         do_LOSpaths, doNadir, radii, alpha, Raycon,     & ! Input
         doCrit, Ncrit, AlphaCrit, RadCrit,              & ! Input
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

   integer, parameter :: ffp = selected_real_kind(15)

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

      real(ffp), intent(in)  :: alpha      (0:maxlayers)
      real(ffp), intent(in)  :: radii      (0:maxlayers)
      real(ffp), intent(in)  :: Raycon

!  Critical stuff

   logical  , intent(in)  :: doCrit
   integer  , intent(in)  :: Ncrit
   real(ffp), intent(in)  :: AlphaCrit
   real(ffp), intent(in)  :: RadCrit

!  Outputs
!  =======

!  Fine layering

      real(ffp), intent(out)  :: alpha_fine (maxlayers,maxfine)
      real(ffp), intent(out)  :: radii_fine (maxlayers,maxfine)

!  Quadratures

      real(ffp), intent(out)  :: xfine    (maxlayers,maxfine)
      real(ffp), intent(out)  :: wfine    (maxlayers,maxfine)

!  Local geoemetry arrays

      real(ffp), intent(out)  :: csqfine  (maxlayers,maxfine)
      real(ffp), intent(out)  :: cotfine  (maxlayers,maxfine)

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message

!  Local
!  -----

      integer, parameter :: Localmaxfine = 20
      integer            :: n, n1, j, nfine
      real(ffp)          :: difh, csfine
      real(ffp)          :: tfine(Localmaxfine), afine(Localmaxfine)

      real(ffp), parameter :: zero = 0.0d0

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
            
            n = NCrit ; difh = radii(n-1) - Radcrit ; nfine = nfinedivs(n)
            call gauleg_ng (0.0d0,difh,tfine,afine,nfine,localmaxfine)
            do j = 1, nfine
               radii_fine(n,j) = radii(n-1) - tfine(j)
               xfine(n,j) = tfine(j)
               wfine(n,j) = afine(j)
            enddo

!    -- For all layers above Critical layer, Regular quadrature only if flagged

            if ( .not. do_LOSpaths ) then
               do n = NCrit-1,1,-1
!               do n = NCrit+1,1,-1
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

         n = NCrit ; n1 = n - 1 ; nfine = nfinedivs(n)
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
            do n = NCrit-1,1,-1
!            do n = NCrit+1,1,-1
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

subroutine SD_incoming_sphergeom                                      &
       ( maxlayers, maxfine, nlayers, nfinedivs, do_Chapman, doNadir, & ! Input
         DoCrit, NCrit, alpha_boa, theta_boa, phi_boa, radii, alpha,  & ! Input
         vsign, dtr, Pie, RadCrit, AlphaCrit, radii_fine, alpha_fine, & ! Input
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,          & ! Output
         chapfacs, cosscat, theta_all, phi_all )                        ! Output

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

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

      integer  , intent(In)    :: maxlayers, maxfine
      integer  , intent(In)    :: nlayers, nfinedivs(maxlayers), NCrit
      real(ffp), intent(InOut) :: alpha_boa, theta_boa, phi_boa
      real(ffp), intent(In)    :: vsign
      logical  , intent(In)    :: do_Chapman, DoCrit, doNadir

!  Los geometry

      real(ffp), intent(In)   :: alpha         (0:maxlayers)
      real(ffp), intent(In)   :: alpha_fine    (maxlayers,maxfine)
      real(ffp), intent(In)   :: radii         (0:maxlayers)
      real(ffp), intent(In)   :: radii_fine    (maxlayers,maxfine)
      real(ffp), intent(In)   :: AlphaCrit, RadCrit, dtr, pie

!  main outputs (geometry)

      integer  , intent(Out)  :: ntraverse  (0:maxlayers)
      real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers)
      real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers)

!  Fine level output (geometry)

      integer  , intent(Out)  :: ntraverse_fine(maxlayers,maxfine)
      real(ffp), intent(Out)  :: sunpaths_fine (maxlayers,maxlayers,maxfine)

!  scattering angle and associated angles

      real(ffp), intent(Out)  :: cosscat
      real(ffp), intent(Out)  :: theta_all  (0:maxlayers)
      real(ffp), intent(Out)  :: phi_all    (0:maxlayers)

!  Local

      logical       :: DirectSun, Do_OverheadSun, Do_ZeroSunBOA, Do_Normal
      integer       :: n, j, k
      real(ffp)     :: SolarDirection(3), Radstart, term1, term2
      real(ffp)     :: salpha_boa, calpha_boa, phi_boa_R, sphi_boa
      real(ffp)     :: theta_boa_R, stheta_boa, ctheta_boa, cphi_boa
      real(ffp)     :: ctheta, stheta, calpha, salpha, cphi, CumAngle

!  Check
!      real(ffp)  :: sumd, sume, sth1

!  Local arrays associated with fine grid output

      logical            :: DirectSunf(maxfine)
      real(ffp)       :: thetaf(maxfine)
      real(ffp)       :: sthetaf(maxfine)
      real(ffp)       :: cthetaf(maxfine)

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

!  Special case

      Do_OverheadSun = theta_boa.eq.0.0d0

!  BOA angles

      if ( alpha_boa.eq.90.0d0 ) then
         calpha_boa     = 0.0d0
         salpha_boa     = 1.0d0
      else
         salpha_boa  = dsin(alpha(nlayers))
         calpha_boa  = dcos(alpha(nlayers))
      endif

      theta_boa_R    = theta_boa * DTR
      if ( theta_boa.eq.90.0d0 ) then
         ctheta_boa     = 0.0d0
         stheta_boa     = 1.0d0
      else
         stheta_boa     = dsin(theta_boa_R)
         ctheta_boa     = dcos(theta_boa_R)
      endif

      phi_boa_R   = phi_boa * dtr
      cphi_boa    = dcos(phi_boa_R)
      sphi_boa    = dsin(phi_boa_R)
      cphi_boa       = dcos(phi_boa * dtr)

!  define Unit solar vector at BOA

      if ( Do_OverheadSun ) then
         SolarDirection = 0.0_ffp
      else
         SolarDirection(1) = - stheta_boa * cphi_boa * vsign
         SolarDirection(2) = - stheta_boa * sphi_boa
         SolarDirection(3) = - ctheta_boa
      endif

!  Cosine of scattering angle at boa

      if ( Do_OverheadSun ) then
         term1 = 0.0d0
         term2 = calpha_boa
         cosscat = - vsign * term2 ; if (term2.eq.0.0d0) cosscat = term2
      else
         term1 = salpha_boa * stheta_boa * cphi_boa
         term2 = calpha_boa * ctheta_boa
         cosscat = - vsign * term2 + term1 
      endif

!  General case: LOS path in spherical geometry
!  ============================================

!  Start loop over all layers

      do n = nlayers, 1, -1

!  Special cases

        DO_ZeroSunBOA  = Do_OverheadSun.and.(n.eq.nlayers.or.doNadir)
        DO_Normal      = .not. doCrit .or. ( doCrit .and. n.le. NCrit )

!  Layer boundary Sun position
!     * Local save of angles, cosines, sines and  illumination flags
!     * Use critical ALPHA and RADIUS if N = NCrit
!     * Use Bottom-of-layer values if N < NCrit (BOA values if illuminated)

!        if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
           if ( doCrit .and. n .eq. NCrit ) then
              CumAngle = alpha(nlayers) - AlphaCrit ; Radstart = RadCrit
              call FindSun(DoNadir,Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R,&
                           theta_all(n),stheta,ctheta,DirectSun)
           else
              Radstart = radii(n)
              if ( n.eq. nlayers ) then
                 theta_all(n) = theta_boa*dtr ; stheta = stheta_boa ; ctheta = ctheta_boa ; DirectSun = .true.
              else
                 CumAngle = alpha(nlayers) - alpha(n)
                 call FindSun(DoNadir,Do_OverheadSun,radii(n),SolarDirection,CumAngle,theta_boa_R,&
                              theta_all(n),stheta,ctheta,DirectSun)
              endif
           endif
!        endif                                               !   @@RTSFix 9/5/12 (Comment out line)

!  Fine-layer sun positions

        if ( Do_Normal ) then
           do j = 1, nfinedivs(n)
              CumAngle = alpha(nlayers) - alpha_fine(n,j)
              call FindSun(DoNadir,Do_OverheadSun,radii_fine(n,j),SolarDirection,CumAngle,theta_boa_R,&
                           thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
           enddo
        endif

!  Sun paths in layer

!        if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
           if ( DirectSun ) then
              call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,&
                theta_all(n),stheta,N,sunpaths(n,:))
           else
              call FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,theta_all(n),stheta,N,sunpaths(n,:),ntraverse(n))
           endif
        if ( Do_Normal ) then                                !   @@RTSFix 9/5/12 (Addline)
           do j = 1, nfinedivs(n) 
              if ( DirectSunf(j) ) then
                 call FindSunPaths_D &
                  (Do_ZeroSunBOA,Maxlayers,Radii_fine(n,j),Radii,&
                   thetaf(j),sthetaf(j),N,sunpaths_fine(n,:,J))
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
!     If AZM > 180, Subtract from 360 for consistency. (VLIDORT code, 10 October 2011)

        if (Do_OverheadSun.or.doNadir ) then
           phi_all(n)     = phi_boa * dtr
        else
!           if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
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
              if ( phi_boa.gt.180.0d0) phi_all(n) = 2.0d0 * Pie - phi_all(n)
!           endif                                               !   @@RTSFix 9/5/12 (Comment out line)
        endif

!  End layer loop

      enddo

!  TOA Sun angle sunpaths and PHI.
!    (No sunpaths if directly illuminated)

      DO_ZeroSunBOA  = Do_OverheadSun.and.doNadir
      CumAngle = alpha(nlayers) - alpha(0) ; Radstart = radii(0)
      call FindSun(DoNadir,Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R,&
                   theta_all(0),stheta,ctheta,DirectSun)
      if (.not.DirectSun ) then
          call FindSunPaths_T (Maxlayers,Pie,Radii(0),Radii,theta_all(0),stheta,1,sunpaths(0,:),ntraverse(0))
      endif
      if ( Do_OverheadSun .or. doNadir ) then
         phi_all(0)     = phi_boa * dtr
      else
         cphi = (cosscat+vsign*calpha*ctheta)/stheta/salpha
         if ( cphi.gt.1.0d0)  cphi = 1.0d0 ; if ( cphi.lt.-1.0d0) cphi = -1.0d0
         phi_all(0)     = dacos(cphi)
         if ( phi_boa.gt.180.0d0) phi_all(0) = 2.0d0 * Pie - phi_all(0)
      endif

!  Chapman factor calculations
!  ---------------------------

      if ( do_Chapman ) then
         do n = 1, nlayers
            call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii(n),Radii,&
              theta_boa_R,stheta_boa,N,chapfacs(n,:))
            do k = 1, n
               chapfacs(n,k) = chapfacs(n,k)/(radii(k-1)-radii(k))
            enddo
         enddo
      endif

!  Finish

      return
end subroutine SD_incoming_sphergeom

!

!  General Routines for Sun positioning

subroutine FindSun(DoNadir,Do_OverheadSun,Radius,SolarDirection,CumAngle,theta_boa_R,theta,stheta,ctheta,DirSun)

!  Find the solar anlge along the LOS path, for given radius and cumulative angle from BOA
!    SolarDirection is defined at BOA, with azimuth relative to the LOS direction.

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs

   logical   , Intent(In)    :: DoNadir,Do_OverheadSun
   real(ffp) , Intent(in)    :: Radius,SolarDirection(3),CumAngle,theta_boa_R

!  Outputs

   real(ffp) , Intent(out)   :: theta,stheta,ctheta
   logical   , Intent(InOut) :: DirSun

!  Local

   real(ffp) :: px(3),b

!  Calculation (Nadir view scenario)

   if ( doNadir ) then
      DirSun = .true.
      theta = theta_boa_R
      ctheta = dcos(theta_boa_R)
      stheta = dsin(theta_boa_R)
      return
   endif

!  Calculation (overhead sun)

   if ( Do_OverheadSun ) then
      DirSun = .true.
      ctheta = dcos(CumAngle)
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


subroutine FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,&
                           thstart,sthstart,N,sunpaths)

!  Sunpaths for the Direct-sun illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N
!  Special case = Overhead sun at BOA

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   LOGICAL   , Intent(In)   :: Do_ZeroSunBOA
   INTEGER   , Intent(In)   :: maxlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart

!  Output

   real(ffp), Intent(InOut) :: Sunpaths(maxlayers)

!  Local

   integer    :: n1, k
   real(ffp)  :: sth0, th0, sth1, th1, ks1

!  Layer boundary upper

   N1 = N - 1

!  SBOA condition

   if ( Do_ZeroSunBOA ) then
      sunpaths(n) = radii(n1) - Radstart
      do k = n1, 1, -1
         sunpaths(k) = radii(k-1) - radii(k)
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

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   INTEGER   , Intent(In)   :: maxlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart,Pie

!  Output

   INTEGER   , Intent(InOut)   :: NT
   real(ffp), Intent(InOut)    :: Sunpaths(maxlayers)

!  Local

   logical    :: trawl
   integer    :: n1, k
   real(ffp)  :: sth0, th0, sth1, th1, ks1, tanr

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

   integer, parameter :: ffp = selected_real_kind(15)

!  Input/Output

      INTEGER  , intent(in)  :: N,NMAX
      REAL(ffp), intent(in)  :: X1, X2
      REAL(ffp), intent(out) :: X(NMAX),W(NMAX)

      INTEGER     :: I, M, J
      REAL(ffp)   :: XM,XL,P1,P2,P3,PP,Z,Z1
      REAL(ffp), PARAMETER :: EPS = 3.0D-14

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

   integer, parameter :: ffp = selected_real_kind(15)

!  I/O

   real(ffp), intent(in)  ::  x, xboa, rtop, raycon, sundir(3)
   real(ffp), intent(out) ::  sundist, theta0

!  Local

   real(ffp) :: xicum, sinx, rad, c0, s0, s1, c1

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

!  Finish

end module FirstOrder_geometry_Mk3_m

