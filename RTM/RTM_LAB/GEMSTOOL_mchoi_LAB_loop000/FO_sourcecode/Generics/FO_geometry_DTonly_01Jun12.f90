

module FO_geometry_DTonly_m

   use FO_geometry_Pool_m, only : STD_outgoing_sphergeom_initial, STD_outgoing_FindCritlayer, &
                                  STD_outgoing_sphergeom_Qbasic , STD_outgoing_sphergeom_Qadjusted

!  Routines for Geometry with Thermal only

public DT_Geometry

contains

subroutine DT_Geometry                                                         &
       ( maxlayers, maxfine, nlayers, nfinedivs, dtr, eradius, heights,   & ! Input
         alpha_boa,  do_LOSpaths, do_planpar, doCrit, Acrit, extinc,      & ! Input
         doNadir, Raycon, radii, alpha, cota, xfine, wfine,               & ! Input/Output
         csqfine, cotfine, alphafine, radiifine, NCrit, RadCrit, CotCrit, & ! Output
         fail, message, trace )                                             ! Output

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

!  input angle (Degrees), dtr = degrees-to-Radians

      real(ffp), intent(in)     :: dtr
      real(ffp), intent(InOut)  :: alpha_boa

!  Flags

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

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message
      character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

      real(ffp)  :: Lospaths (maxlayers)

!  Other angles

      real(ffp)  :: cosa       (0:maxlayers)
      real(ffp)  :: sina       (0:maxlayers)

!  Critical values

      real(ffp)  :: AlphaCrit

!  Help variables

      real(ffp), parameter   :: zero    = 0.0d0
      real(ffp)              :: cutoff, alpha_boa_R

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   NCrit = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero
   cutoff = -log(ACrit)

!  check range of inputs
!  ---------------------

!  VZA can be 0-90 degrees inclusive, but not outside this range

   if ( alpha_boa.gt.90.0d0.or.alpha_boa.lt.0.0d0 ) then
      message = 'Boa LOS angle outside range [0,90]); Check it!'
      trace   = 'Initial Angle Check in SS_Geometry_1'
      fail    = .true. ;  return
   endif

!  Plane-parallel, just set alpha and return
!  -----------------------------------------

   if ( do_planpar ) then
      alpha = 0.0d0
      alpha_boa_R      = alpha_boa * DTR
      alpha(1:nlayers) = alpha_boa_R
      return
   endif

!  Enhanced PS; proceed in 3 Steps
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
         trace = 'Error from STD_Outgoing_FindCritLayer in DT_Geometry_1' ; return
      endif
   endif

!  Step 3. Set Quadratures
!  -----------------------

!  Step 3a. LOS fine-layer quadratures (Regular, non-adjusted, no Criticality)
!           Regular quadratures may be avoided if do_LOSpaths is set

   CALL STD_outgoing_sphergeom_Qbasic                       &
       ( maxlayers, maxfine, nlayers, nfinedivs,            & ! Input
         doNadir, radii, alpha, Raycon,                     & ! Input
         radiifine, alphafine, xfine, wfine,                & ! Output
         csqfine, cotfine, fail, message )                    ! Output

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
         trace = 'Error from STD_outgoing_sphergeom_Qadjusted in DT_Geometry_1' ; return
      endif
   endif

!  Finish

   return
end subroutine DT_Geometry

!  Finish

end module FO_geometry_DTonly_m

