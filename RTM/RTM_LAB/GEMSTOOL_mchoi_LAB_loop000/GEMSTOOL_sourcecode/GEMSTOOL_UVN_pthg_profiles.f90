module GEMSTOOL_UVN_pthgas_profiles_m

!  UVN-specific Type Structures. 12 August 2013
!  Should change 1) nlayer at line 140 2)filename at line 179


private
public :: GEMSTOOL_UVN_PTHGAS_PROFILES

contains

subroutine GEMSTOOL_UVN_PTHGAS_PROFILES &
           ( Maxlay, MaxGases, ngases, which_gases, Tshift, FDGas, FDLev, FDsfp, FDeps,  &  ! Input
             PhysicsDataPath, PTH_profile_name, GAS_profile_names,                       &  ! Input
             nlevels, level_heights, level_temps, level_press, level_vmrs,               &  ! Output (LEVELS)
             nlayers, aircolumns, daircolumns_dS, daircolumn_dP, Tshape,                 &  ! Output (LAYERS). Air Density
             gascolumns, dgascolumns_dS, dgascolumn_dP, dgascolumns_dv,                  &  ! Output (LAYERS). Gas densities
             fail, message )                                                                ! output

  implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  INPUTS
!  ======

!  dimensioning

   integer, intent(in)       :: Maxlay, MaxGases

!  Gas control

   integer                                 , intent(in) :: NGASES
   character(LEN=4), dimension ( MaxGases ), intent(in) :: WHICH_GASES

!  temperature shift

   real(fpk), intent(in) :: Tshift

!  Indices for Finite Difference testing, Perturbation FDeps
!    FDGas Must be 1, 2, 3.... up to ngases
!    FDLev Must be 0, 1, 2, 3 ... up to nlayers

   logical  , intent(in) :: FDsfp
   integer  , intent(in) :: FDGas, FDLev
   real(fpk), intent(in) :: FDeps

!  File paths and file names

   character*(*), intent(in) :: PhysicsDataPath, PTH_profile_name, GAS_profile_names ( Maxgases)

!  OUTPUTS
!  =======

    !      Trace gas Layer column amounts and cross-sections
    !      Units should be consistent, e.g:

    !      gas X-sections   are in cm^2/mol  or inverse [DU]
    !      gas columns      are in mol/cm^2  or [DU]

    !      O2O2 X-sections  are in cm^5/mol^2    ! ALWAYS
    !      O2O2 columns     are in mol^2/cm^5    ! ALWAYS

!  Levels output

    integer                                      , intent(out)     :: NLEVELS
    REAL(fpk)   , DIMENSION( 0:maxlay )          , intent(out)     :: LEVEL_TEMPS
    REAL(fpk)   , DIMENSION( 0:maxlay )          , intent(out)     :: LEVEL_HEIGHTS
    REAL(fpk)   , DIMENSION( 0:maxlay )          , intent(out)     :: LEVEL_PRESS
    REAL(fpk)   , DIMENSION( 0:maxlay, maxgases ), intent(out)     :: LEVEL_VMRS

!  Layer outputs

    integer                                       , intent(out)    :: NLAYERS
    REAL(fpk)   , DIMENSION( maxlay )             , intent(out)    :: AIRCOLUMNS
    REAL(fpk)   , DIMENSION( maxlay )             , intent(out)    :: dAIRCOLUMNS_dS
    REAL(fpk)                                     , intent(out)    :: dAIRCOLUMN_dP

    REAL(fpk)   , DIMENSION( maxlay, maxgases, 2 ), intent(out)    :: GASCOLUMNS
    REAL(fpk)   , DIMENSION( maxlay, maxgases, 2 ), intent(out)    :: dGASCOLUMNS_dS
    REAL(fpk)   , DIMENSION( maxlay, maxgases, 2 ), intent(out)    :: dGASCOLUMNS_dv
    REAL(fpk)   , DIMENSION( maxgases )           , intent(out)    :: dGASCOLUMN_dP

!  Shape function (needed for later)

    REAL(fpk)   , DIMENSION ( 0:maxlay )          , intent(out)    :: TSHAPE

!  Excpetion handling

    logical      , intent(out) :: fail
    character*(*), intent(out) :: message

!  Local variables
!  ---------------

!  help

    integer            :: n, n1, ndum, g, g1, k, score
    real(fpk)          :: constant, O2O2cons, diff, hcon, hcon_O2O2, hcon_O2, refpress, reftemp, dum, gasdum(6), zsur, avit, ccon
    real(fpk)          :: rho1, rho2, drho1_dS, drho2_dS,  vmr1(maxgases), vmr2(maxgases)
    character*256      :: file_pth, file_gas

!  Loschmidt's number (particles/cm3), STP values

    real(fpk)   , PARAMETER :: RHO_STANDARD = 2.68675D+19
    real(fpk)   , PARAMETER :: PZERO = 1013.25D0
    real(fpk)   , PARAMETER :: TZERO = 273.15D0

    real(fpk)   , PARAMETER :: xo2 = 0.209476_fpk, kb = 1.381e-23_fpk

!  Zero the output
!  ===============

   fail = .false. ; message = ' '

   nlevels = 0 ; nlayers = 0

   level_heights = 0.0_fpk
   level_press   = 0.0_fpk
   level_temps   = 0.0_fpk
   level_vmrs    = 0.0_fpk

   aircolumns     = 0.0_fpk ; gascolumns     = 0.0_fpk
   daircolumns_dS = 0.0_fpk ; dgascolumns_dS = 0.0_fpk ; dgascolumns_dV = 0.0_fpk
   daircolumn_dP  = 0.0_fpk ; dgascolumn_dP  = 0.0_fpk

!  Uniform T-shift profile shape

   Tshape = 1.0_fpk

!  PTH profiles
!  ============

!  Open PTH file and read the data.
!      Assume DISCOVER AQ 75-level data, then we get VMRs at the same Pth levels/

!  Set the number of layers, check dimensions immediately]

   nlevels = 36 ; nlayers = nlevels - 1 
   !nlevels = 75 ; nlayers = nlevels - 1 
   if ( nlayers .gt. maxlay ) then
      message = 'NLAYERS > Dimension MAXLAYERS ==> Increase VLIDORT parameter MAXLAYERS'
      fail = .true. ; return
   endif

!  Now read the file

   file_pth = adjustl(trim(PhysicsDataPath))//'PTH_PROFILES/'//adjustl(trim(Pth_profile_name))

   open(1,file=trim(file_pth),err=90,status='old')
   do k = 1, 11
      read(1,*)
   enddo
   do n = nlevels, 1, -1
      n1 = n - 1 ; read(1,*)k,refpress, reftemp, dum,dum,dum,dum,dum,dum
      level_press(n1) = refpress
!      level_temps(n1) = reftemp + Tshift *  Tshape(n1)
      level_temps(n1) = reftemp
   enddo

!  heights by Hydrostatic Eqn. (includes TOA)
!    For now, zsur = surface height in [m] ------Topography data set ????

!   zsur = 0.0d0
!   Level_heights(nlayers) = zsur / 1000.0d0
!   ccon = - 9.81d0 * 28.9d0 / 8314.0d0 * 500.0d0
!   do n = nlayers, 1, -1
!      avit = (1.0d0/level_temps(n-1))+(1.0d0/level_temps(n))
!      level_heights(n-1) = level_heights(n) - log(Level_press(n)/Level_press(n-1))/avit/ccon
!   enddo

!  read from file

!   open(1,file='UVN_heightgrid_1.dat',status='replace')
!   open(1,file='UVN_heightgrid_1.dat',status='old')
!    open(1,file='UVN_heightgrid_1_mchoi_save.dat',status='replace')
!    open(1,file='UVN_heightgrid_0km_US76_fin_mchoi.dat',status='old')
   open(1,file='UVN_heightgrid_0km_10km_US76_ver20181108_mchoi.dat',status='old')
   do n = 0, nlayers
!      write(1,*)n1,level_heights(n)
      read(1,*)n1,level_heights(n)
   enddo
   close(1)

!  Perturbation test

   if ( FDsfp ) then
      level_press(nlayers) = level_press(nlayers) * ( 1.0d0 + FDeps )
   endif

!  Open GAS files and read the data
!    @@@@ 7/25/13. Assume for now that GAS profiles use same level scheme as PTH
!                  Assume everything is VMR (not ppmv)
!                  Assume DISCOVER AQ data, then we get VMRs at the same Pth levels/
!                  Assume only 5 gases available ---> NO BRO

   score = 0
   do g = 1, ngases
      if ( which_gases(g).eq.'H2O ')score = score + 1
      if ( which_gases(g).eq.'O3  ')score = score + 1
      if ( which_gases(g).eq.'NO2 ')score = score + 1
      if ( which_gases(g).eq.'HCHO')score = score + 1
      if ( which_gases(g).eq.'SO2 ')score = score + 1
      if ( which_gases(g).eq.'O2O2' .or. which_gases(g).eq.'O2  ' ) score = score + 1
   enddo
   if ( score .ne. ngases ) then
       message = 'One of the gases is not on the LIST'
       fail = .true. ; return
   endif

! Data

   level_vmrs = 0.0d0
   do g = 1, ngases
      if ( which_gases(g).ne.'O2O2' .or. which_gases(g).ne.'O2  ' ) then
         g1 = g
         file_gas = adjustl(trim(PhysicsDataPath))//'GAS_PROFILES/'//adjustl(trim(Gas_profile_names(g)))
         open(1,file=trim(file_gas),err=91,status='old')
         do k = 1, 11
            read(1,*)
         enddo
         do n = nlevels, 1, -1
            n1 = n - 1 ; read(1,*)ndum,dum,dum,(gasdum(k),k=1,6)
            if ( which_gases(g).eq.'H2O ')level_vmrs(n1,g) = gasdum(1)
            if ( which_gases(g).eq.'O3  ')level_vmrs(n1,g) = gasdum(2)
            if ( which_gases(g).eq.'NO2 ')level_vmrs(n1,g) = gasdum(3)
            if ( which_gases(g).eq.'HCHO')level_vmrs(n1,g) = gasdum(4)
            if ( which_gases(g).eq.'SO2 ')level_vmrs(n1,g) = gasdum(5)
         enddo
         close(1)
      endif
      if ( which_gases(g).eq.'O2  '  ) level_vmrs(0:nlayers,g) = xo2             ! Not really used
      if ( which_gases(g).eq.'O2O2'  ) level_vmrs(0:nlayers,g) = xo2 * xo2       ! Dummy, not used
      if ( which_gases(g).eq.'BRO ') then
         message = 'BRO not available for this data set --> Remove from list of gases!'
         fail = .true. ; return
      endif
   enddo

!  No BRO in the
!  Construct the air and gas columns, Tshift, SurfPress and VMR derivatives
!  ========================================================================

!  Number of layers = number of levels MINUS 1

   nlayers = nlevels - 1

!  Gas constant, O2O2 constant

   CONSTANT  = 1.0D+05 * RHO_STANDARD * TZERO / PZERO
   O2O2CONS = xo2 * 1.0D-4 / kb

!  Initial values at TOA ( Level 0 )

   rho1     = level_press(0) / level_temps(0)
   drho1_dS = - rho1 * Tshape(0) / level_temps(0) 
   vmr1(1:ngases) = level_vmrs(0,1:ngases)            ! @@ OK for all gases, since level_vmrs is zeroed

   DO N = 1, NLAYERS
      DIFF = LEVEL_HEIGHTS(n-1) - LEVEL_HEIGHTS(n)
      HCON = 0.5D0 * CONSTANT * DIFF

!  Values at level N
!      rho1     = level_press(n-1) / level_temps(n-1)  !for NGST
      rho2     = level_press(n) / level_temps(n)
      drho2_dS = - rho2 * Tshape(n) / level_temps(n) 
      vmr2(1:ngases) = level_vmrs(n,1:ngases)   ! @@ OK for all gases, since level_vmrs is zeroed

!  Air columns + Tshift derivative

      AIRCOLUMNS(N)      = HCON * ( rho1 + rho2 )
      dAIRCOLUMNS_dS(N)  = HCON * ( dRHO1_dS + dRHO2_dS )

!  gas columns and VMR + TSHIFT Derivatives

      do g = 1, ngases
         if ( which_gases(g).eq.'O2O2' ) then
            HCON_O2O2 =  0.5D0 * O2O2CONS * O2O2CONS * DIFF * 1.0D+5
            GASCOLUMNS(N,G,1)     = HCON_O2O2 * RHO1 * RHO1
            GASCOLUMNS(N,G,2)     = HCON_O2O2 * RHO2 * RHO2
            dGASCOLUMNS_dS(N,G,1) = HCON_O2O2 * RHO1 * dRHO1_dS * 2.0d0
            dGASCOLUMNS_dS(N,G,2) = HCON_O2O2 * RHO2 * dRHO2_dS * 2.0d0
         else if ( which_gases(g).eq.'O2' ) then
            HCON_O2 = XO2 * HCON
            GASCOLUMNS(N,G,1)     = HCON_O2 * RHO1
            GASCOLUMNS(N,G,2)     = HCON_O2 * RHO2
            dGASCOLUMNS_dS(N,G,1) = HCON_O2 * dRHO1_dS
            dGASCOLUMNS_dS(N,G,2) = HCON_O2 * dRHO2_dS
         else
            GASCOLUMNS(N,G,1)     = HCON * RHO1 * VMR1(g)
            GASCOLUMNS(N,G,2)     = HCON * RHO2 * VMR2(g)
            dGASCOLUMNS_dS(N,G,1) = HCON * dRHO1_dS * VMR1(g)
            dGASCOLUMNS_dS(N,G,2) = HCON * dRHO2_dS * VMR2(g)
            dGASCOLUMNS_dV(N,G,1) = HCON * RHO1      ! upper boundary
            dGASCOLUMNS_dV(N,G,2) = HCON * RHO2      ! Lower boundary
         endif
      enddo

!  Surface pressure derivatives

      if ( n.eq.nlayers) then
         DAIRCOLUMN_dP = HCON / Level_Temps(n)
         DGASCOLUMN_dP(1:ngases) = DAIRCOLUMN_dP * VMR2(1:ngases) 
      endif

!  Update the level information

      rho1 = rho2 ; drho1_dS = drho2_dS
      vmr1(1:ngases) = vmr2(1:ngases)          ! @@ OK for all gases, since level_vmrs is zeroed

!  ENd layer loop

   ENDDO

!   write(*,*)sum(Aircolumns(1:nlayers)) 

!  normal Finish

   return

!  error finish

90 continue
   fail    = .true.
   message = 'PTH file not found' 
   return
91 continue
   fail    = .true.
   message = 'GAS file not found : '//which_gases(g1) 
   return

!  End

end subroutine GEMSTOOL_UVN_PTHGAS_PROFILES

!  End module

end module GEMSTOOL_UVN_pthgas_profiles_m


