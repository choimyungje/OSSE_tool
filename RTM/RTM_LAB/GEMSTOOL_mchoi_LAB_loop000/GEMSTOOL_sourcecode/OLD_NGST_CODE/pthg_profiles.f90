module pthg_profiles_m

!  New, 18 february 2013
!  ---------------------

!  TOMS Standard Profile PTH and O3 (2 routines). J. Park

!  Type definitions

   use vlidort_pars, only : fpk,vlidort_spkind, zero, one, half, two

!  Use numerical and dataread_101 Modules

   use dataread_101_m, only : read_101prf
   use numerical_m   , only : lintp2, linear_trawler_ascending

   integer, parameter :: spk = vlidort_spkind

!  Only main routine is public, private routines are
!        pthg_profiles_z101
!        extract_usadata_pth
!        develop_geoschem_pth
!        extract_geoschem_so2
!        develop_geoschem_so2
!        extract_tomsv8_o3
!        develop_tomsv8_o3
!        TOMS_STDPRF_T, TOMS_STDPRF_O3   ( New, 18 February, after J. Park)

private
public :: pthg_profiles, pthg_finer

contains

subroutine pthg_profiles                                               &
  ( maxlay, maxlev, maxgas, do_gas_preserve, do_gas_jacobians,         & ! INPUT
    do_profile_linearization, do_column_linearization, PTH_SOURCES,    & ! INPUT
    WHICH_GASES, NGASES, GAS_SOURCES, GAS_INPUT_COLUMNS,               & ! INPUT
    LATITUDE, LONGITUDE, YEAR, MONTH, DAY_OF_MONTH, surface_pressure,  & ! INPUT
    LEVEL_HEIGHTS, LEVEL_PRESS, LAYER_AIRCOLUMNS, LAYER_TEMPERATURES,  & ! OUTPUT
    GAS_PROFILES, GAS_PROFILES_DERIV,                                  & ! OUTPUT
    GAS_OUTPUT_COLUMNS, NLAYERS,                                       & ! OUTPUT
    FAIL, MESSAGE )

  implicit none

!  Key for Sources

!     for PTH : PTH_SOURCES    = 1 --> TOMS V8 CLIMATOLOGY pressure grid  + interpolated USA heights/temperatures
!     for PTH : PTH_SOURCES    = 2 --> z101, ozone only
!     for PTH : PTH_SOURCES    = 3 --> GEOSCHEM pressure grid + interpolated USA heights/temperatures
!     for PTH : PTH_SOURCES    = 4 --> GEOSCHEM pressure grid + interpolated TOMS temperatures, hydrostatic H
!     for PTH : PTH_SOURCES    = 5 --> TOMS V8 STANDARD PROFILES   @@@@@@ New 2/18/13

!     for O3 : GAS_SOURCES(G)  = 1 --> TOMS V8 CLIMATOLOGY ozone
!     for O3 : GAS_SOURCES(G)  = 2 --> z101, ozone only
!     for O3 : GAS_SOURCES(G)  = 3 --> TOMS V8 STANDARD PROFILES, O3 only  @@@@@@ New 2/18/13

!     for SO2: GAS_SOURCES(G) = 1 --> GEOS_CHEM So2 shape factors + profiles.

!  New, 18 february 2013
!  ---------------------

!  TOMS Standard Profile PTH and O3 (2 routines). J. Park

!  INPUTS
!  ======

!  dimensioning

   integer       :: maxlay, maxlev, maxgas

!  Control flags

   logical       :: do_gas_preserve
   logical       :: do_gas_jacobians
   logical       :: do_profile_linearization
   logical       :: do_column_linearization

!  Gas and PTHG control

   integer                                :: NGASES
   character(LEN=4), dimension ( maxgas ) :: WHICH_GASES
   integer,          dimension ( maxgas ) :: GAS_SOURCES
   real(fpk)   ,     dimension ( maxgas ) :: GAS_INPUT_COLUMNS

   integer                                :: PTH_SOURCES

!  Time and position

   real(fpk)    :: latitude, longitude
   integer      :: year, month, day_of_month

!  surface pressure

   real(fpk)    :: surface_pressure

!  OUTPUTS
!  =======

    !      Trace gas Layer column amounts and cross-sections
    !      Units should be consistent, e.g:
    !      gas X-sections   are in cm^2/mol  or inverse [DU]
    !      gas columns      are in mol/cm^2  or [DU]
    !      O2O2 X-sections  are in cm^5/mol^2
    !      O2O2 columns     are in mol^2/cm^5

    integer                                           :: NLAYERS
    REAL(fpk)   , DIMENSION( maxlay )                 :: LAYER_AIRCOLUMNS
    REAL(fpk)   , DIMENSION( maxlay )                 :: LAYER_TEMPERATURES
    REAL(fpk)   , DIMENSION( 0:maxlay )               :: LEVEL_HEIGHTS
    REAL(fpk)   , DIMENSION( 0:maxlay )               :: LEVEL_PRESS

    REAL(fpk)   , DIMENSION( maxlay, maxgas )         :: GAS_PROFILES
    REAL(fpk)   , DIMENSION( maxlay, maxgas )         :: GAS_PROFILES_DERIV
    REAL(fpk)   , DIMENSION( maxgas )                 :: GAS_OUTPUT_COLUMNS 

!  Actually, this is local

    REAL(fpk)   , DIMENSION( maxlev   )               :: LEVEL_LOGPRESS

!  Excpetion handling

    logical       :: fail
    character*(*) :: message

!  Local variables
!  ===============

!  help

    integer            :: G, n
    real(fpk)          :: total_o3, sum
    CHARACTER*100      :: IFILE

!  SO2 GEOSCHEM extraction variables
!  ---------------------------------

!  Flag for first read of the Geoschem data

      logical          :: do_first_read

!  GEOSCHEM Data outputs

      integer          :: geoschem_month
      real(fpk)        :: so2gc(91,144,30)
      real(fpk)        :: etagc(30),latgc(91),longc(144)

!  Extracted output

      real(fpk)        :: so2gc_shapefactors(30)

!  Vacuum pressure (here, hardwired to 0.01 hPa )

    real(fpk)   , parameter :: vacuum_pressure = 0.01d0

!  USA data
!  --------

!  Do standard flag (here, not set)

    logical, parameter :: do_usa_standard = .false.

!  PTH and number of levels ( = 46 )

    integer                     :: USA_NLEVELS
    real(fpk)   , dimension(46) :: USA_HEIGHTS     (46)
    real(fpk)   , dimension(46) :: USA_LOGPRESSURE (46)
    real(fpk)   , dimension(46) :: USA_TEMPERATURE (46)

!  z101 extraction
!  ---------------

    real(fpk)   , dimension(maxlay) :: UMKEHR_O3PROFILE
    real(fpk)   , dimension(maxlay) :: DC_UMKEHR_O3PROFILE

!  Local O2O2 column
!  -----------------

    REAL(fpk)   , DIMENSION(maxlay) :: LAYER_O2O2COLUMNS

!  PTH Sources
!  ===========

!  PTH source = 5 (TOMS Standard Temperatures)

   IF ( PTH_SOURCES .EQ. 5 ) THEN
      IFILE = 'physics_data/toms_stdprf_t.dat'
      CALL TOMS_STDPRF_T                                   &
             ( maxlay, IFILE, SURFACE_PRESSURE, NLAYERS, LEVEL_PRESS,    &
               LEVEL_LOGPRESS, LEVEL_HEIGHTS, LAYER_TEMPERATURES, LAYER_AIRCOLUMNS, FAIL, MESSAGE )
        if ( fail ) return
   ENDIF

!  PTH source = 2, Use the complete z101 specification

   IF ( PTH_SOURCES .EQ. 2 ) THEN
     call pthg_profiles_z101                                           &
            ( maxlay, maxlev, do_gas_jacobians, surface_pressure,      & 
              do_profile_linearization, do_column_linearization,       & 
              LEVEL_HEIGHTS, LEVEL_PRESS, LEVEL_LOGPRESS,              &
              LAYER_AIRCOLUMNS, LAYER_O2O2COLUMNS, LAYER_TEMPERATURES, & 
              UMKEHR_O3PROFILE, DC_UMKEHR_O3PROFILE,                   &
              total_o3, NLAYERS )
   ENDIF

!  PTH source = 1 or 3, Need to get the USA atmospheres

   IF ( PTH_SOURCES == 1 .or. PTH_SOURCES == 3 ) THEN
     call extract_usadata_pth ( do_usa_standard, latitude, month, &
            USA_NLEVELS,USA_HEIGHTS,USA_LOGPRESSURE, USA_TEMPERATURE)
   ENDIF

!  PTH source = 3, Extract GEOSCHEM pressures (and SO2 shape factors)

   IF ( PTH_SOURCES .EQ. 3 ) THEN
      geoschem_month = month
      geoschem_month = 7             ! Tempo
      do_first_read = .true.
      call extract_geoschem_so2                                &
        ( latitude, longitude, geoschem_month, do_first_read,  & ! input
          etagc, latgc, longc, so2gc,                          & ! Output
          so2gc_shapefactors, fail, message )                    ! Output
        if ( fail ) return
      call develop_geoschem_pth                                           &
            ( maxlay, maxlev, surface_pressure, vacuum_pressure, etagc,   & ! Input
              usa_nlevels, usa_heights, usa_logpressure, usa_temperature, & ! Input
              nlayers, level_heights, level_logpress,                     & ! Output
              layer_aircolumns, layer_o2o2columns, layer_temperatures )     ! Output
   ENDIF

!  debug PTH
!     do g = 1, nlayers
!        write(*,*)g,layer_temperatures(g),dexp(level_logpress(g))
!     enddo
!     write(*,*)PTH_SOURCES, NGASES, NLAYERS
!    pause

!  Gas sources
!  ===========

!  initialize

     GAS_PROFILES       = 0.0d0
     GAS_PROFILES_DERIV = 0.0d0

!  Start loop

     DO G = 1, NGASES

!  OZONE (O3)
!  ----------

      IF ( WHICH_GASES(G) == 'O3  ' ) THEN

!  ozone from TOMS standard profile

         IF ( GAS_SOURCES(G) == 3 ) then
            IFILE = 'physics_data/toms_stdprf_o3.dat'
            CALL TOMS_STDPRF_O3  &
              ( GAS_INPUT_COLUMNS(G), LATITUDE, IFILE, MAXLAY, NLAYERS,    & ! Input
                GAS_PROFILES(:,G), GAS_PROFILES_DERIV(:,G), FAIL, MESSAGE )  ! Output
            if ( fail ) return
         ENDIF

!  ozone from z101, just copy into output

         IF ( GAS_SOURCES(G) == 2 ) then
           GAS_OUTPUT_COLUMNS(G)           = total_o3
           GAS_PROFILES      (1:NLAYERS,G) =    UMKEHR_O3PROFILE(1:NLAYERS)
           GAS_PROFILES_DERIV(1:NLAYERS,G) = DC_UMKEHR_O3PROFILE(1:NLAYERS)
         ENDIF

!  Ozone from TOMS V8. Also Get PTH from TOMS V8 pressures if using that source
 
         IF ( GAS_SOURCES(G) == 1 ) then
            IF ( PTH_SOURCES .EQ. 1 ) THEN
              call develop_tomsv8_o3pth                                     &
              ( GAS_INPUT_COLUMNS(G), LATITUDE, YEAR, MONTH, DAY_OF_MONTH,  & ! Input
                MAXLAY, MAXLEV, surface_pressure, do_gas_preserve,          & ! Input
                usa_nlevels, usa_heights, usa_logpressure, usa_temperature, & ! Input
                nlayers, level_heights, level_press, level_logpress,        & ! Output
                layer_aircolumns, layer_o2o2columns, layer_temperatures,    & ! Output
                GAS_PROFILES(:,G), GAS_PROFILES_DERIV(:,G),                 & ! Output
                FAIL,  MESSAGE )
            ELSE
              call develop_tomsv8_o3                           &
              ( GAS_INPUT_COLUMNS(G), LATITUDE, YEAR, MONTH,   & ! Input
                DAY_OF_MONTH, MAXLAY, MAXLEV, do_gas_preserve, & ! Input
                NLAYERS, level_logpress, surface_pressure,     & ! Input
                GAS_PROFILES(:,G), GAS_PROFILES_DERIV(:,G),    & ! Output
                FAIL,  MESSAGE )
            ENDIF
            if ( fail ) return
         ENDIF

!  profile linearization, do not require the column derivatives

         if ( do_profile_linearization ) then
            do n = 1, nlayers
               gas_profiles_deriv(n,g) = 1.0d0
            enddo
         endif

!  debug O3

!     sum = 0.0d0
!     do n = 1, nlayers
!        sum = sum + gas_profiles(n,g)
!        write(*,*)n,gas_profiles(n,g),gas_profiles_deriv(n,g)
!     enddo
!     write(*,*)sum/2.66868d+16
!     pause'o3'

!  SULFUR DIOXIDE
!  --------------

      ELSE IF ( WHICH_GASES(G) == 'SO2 ' ) THEN

!  Geoschem - develop Umkehr profile from GEOSCHEM shape-factors
!    AMOUNT of SO2 is always preserved (03 September 2012)

         IF ( GAS_SOURCES(G) == 1 ) then
           call develop_geoschem_so2                         &
         ( MAXLAY, do_column_linearization,                  & ! Input
           NLAYERS, GAS_INPUT_COLUMNS(G),                    & ! Input
           layer_aircolumns, so2gc_shapefactors,             & ! Input
           GAS_PROFILES(:,G), GAS_PROFILES_DERIV(:,G) )        ! Output
         ENDIF

!  debug SO2
!     do n = 1, nlayers
!      write(*,*)n,layer_aircolumns(n), gas_profiles(n,g)
!     enddo
!     pause

!  Other SO2 options (plume)

!  O2-O2 absoprtion (introduced 03 September 2012)
!  -----------------------------------------------

!  Just copy the automatically calculated profile

      ELSE IF ( WHICH_GASES(G) == 'O2O2' ) THEN

         do N = 1, nlayers
            gas_PROFILES(n,G) = layer_o2o2columns(n)
         enddo
     
!  OTHER GASES
!  -----------

      ELSE

!                      P L A C E H O L D E R

      ENDIF

!  Output columns (this is just a check)

      SUM = 0.0d0
      do N = 1, nlayers
        sum = sum + gas_PROFILES(n,G)
      enddo
      gas_output_columns(g) = sum

!  END GAS LOOP
!  ------------

   ENDDO

!  Make Output pressures

   do n = 0, nlayers
     level_press(n) = dexp(level_logpress(n+1))
   enddo

!  Normal FInish

   RETURN

END SUBROUTINE pthg_profiles



subroutine pthg_profiles_z101                                &
  ( maxlev, maxlev1, do_gas_jacobians, surface_pressure,     & ! INPUT
    do_profile_linearization, do_column_linearization,       & ! INPUT
    LEVEL_HEIGHTS, LEVEL_PRESS, LEVEL_LOGPRESS,              & ! OUTPUT
    LAYER_AIRCOLUMNS, LAYER_O2O2COLUMNS, LAYER_TEMPERATURES, & ! OUTPUT
    LAYER_GASCOLUMNS, LAYER_GASDERIVS,                       & ! OUTPUT
    TOTAL_COLUMN, NLAYERS )                                    ! OUTPUT

!  Delivers PTH + O3 (from fixed file). One gas only

      implicit none

!  INPUTS
!  ======

!  dimensioning

   integer       :: maxlev, maxlev1

!  Control flags

   logical       :: do_gas_jacobians
   logical       :: do_profile_linearization
   logical       :: do_column_linearization

!  surface pressure

    real(fpk)    :: surface_pressure

!  OUTPUTS
!  =======

    !      Trace gas Layer column amounts and cross-sections
    !      Units should be consistent, e.g:
    !      X-sections  are in cm^2/mol  or [DU]
    !      gas columns are in mol/cm^2  or inverse [DU]

    integer                                   :: NLAYERS
    REAL(fpk)   , DIMENSION( maxlev )         :: LAYER_AIRCOLUMNS
    REAL(fpk)   , DIMENSION( maxlev )         :: LAYER_O2O2COLUMNS
    REAL(fpk)   , DIMENSION( maxlev )         :: LAYER_TEMPERATURES
    REAL(fpk)   , DIMENSION( 0:maxlev )       :: LEVEL_PRESS
    REAL(fpk)   , DIMENSION( 0:maxlev )       :: LEVEL_HEIGHTS
    REAL(fpk)   , DIMENSION( maxlev1 )        :: LEVEL_LOGPRESS

    REAL(fpk)   , DIMENSION( maxlev )         :: LAYER_GASCOLUMNS 
    REAL(fpk)   , DIMENSION( maxlev )         :: LAYER_GASDERIVS 
    REAL(fpk)                                 :: TOTAL_COLUMN

!  LOCAL VARIABLES
!  ===============

!  Allocatable local arrays

!    real(fpk)   , dimension(:), allocatable :: z101, z101_old, p101, t101,  x101, x101c
!    real(fpk)   , dimension(:), allocatable :: ploginv, p101log, x101log, x101clog
!    real(fpk)   , dimension(:), allocatable :: pprf, zprf, tprf, xprf

    real(fpk)   , dimension(101) :: z101, z101_old, p101, t101,  x101, x101c
    real(fpk)   , dimension(101) :: ploginv, p101log, x101log, x101clog
    real(fpk)   , dimension(101) :: pprf, zprf, tprf, xprf

!  parameters

    real(fpk)   , PARAMETER :: rgas=8314.32    ! gas constant [N*m]/[mol*K]
    real(fpk)   , PARAMETER :: mair=28.9664    ! mean molecular mass of air [g/mol]
    real(fpk)   , PARAMETER :: grav=9.80665    ! accel due to gravity at surface [m/s^2]

!  help variables

    real(fpk)    :: mg_over_r, frac, psfc2, z, t, dlogp, dzhalf, dz, r0, aconst
    real(fpk)    :: dens1, dens2, hdiff, toa_height, rho_stp, xo2_kb
    integer      :: i, nplevprf, nplev, idx
    logical  , parameter :: new_treatment = .true.
    real(fpk), parameter :: xo2 = 0.209476_fpk, kb = 1.381e-23_fpk

!  read profiles
!  -------------

!  get the data

  nplev=101
!  allocate(pprf(nplev),zprf(nplev),tprf(nplev),xprf(nplev))
  call read_101prf( 'physics_data/uu.prf', pprf, zprf, tprf, xprf, nplevprf )

!  Convert temperatures to deg K
  
  tprf=tprf+273.15d0

  ! Test for nplev = nplevprf

  ! Define atmosphere's pressure grid straight from supplied input profile.
  ! Assign corresponding temperatures for each pressure level, and ozone
  ! amounts for each layer from supplied input profile as well.

  ! Pressure and temperature at level i corresponds to bottom of layer i.
  ! Ozone amount is specified for layer i and converted to atm-cm.
  ! p(1) is at the top of the atmosphere.
  
!  allocate(p101(nplev),t101(nplev),x101(nplev),z101(nplev),z101_old(nplev))
!  allocate(p101log(nplev), x101log(nplev), x101c(nplev), x101clog(nplev), ploginv(nplev))

  p101 = pprf
  p101log = dlog(p101)
  t101 = tprf
  x101 = xprf/1000.
  x101log = dlog(x101)
  
  ! Compute cumulative ozone at all pressure levels. Units are atm-cm.
  
  !x101c(1) = xprf(1)
  !do i=2, nplev
  !   x101c(i) = x101c(i-1)+xprf(i)
  !end do
  !x101clog = dlog(x101c)

! ####################################################################################

!  Apply Hydrostatic Eq.
!  ---------------------

  ! Compute heights z in km. corresponding to the atmospheric pressure levels p
  !     in inverse order (1 <- 101) using equations :
  !     dz/dlogp=n*r*t/g
  !     g=g0*(r0/(r0+z))**2
  !
  !     constant is g0/(n*r); z(nplev) is at bottom of atmosphere
  !

 ! apply hydrostatic equation
 ! dp/p= -(m*g)/(r*t) dz

   mg_over_r=(mair*grav)/rgas
   aconst=1/mg_over_r
  
 ! first convert g to km
   aconst=aconst/1000.

 ! r0: 6356.76  ! radius of earth [km]
  r0 = 6356.76   

!-------------------------------------
! Consistent z-heights with Dave Haffner.
!-------------------------------------
  
  z101_old(1)=0.d0
  dzhalf = 0.482d0
!
  do i= 2, nplev
     dlogp = p101log(i) - p101log(i-1)
     z  = z101_old(i-1)+dzhalf
     t  = t101(i-1)
     dz = -aconst*((r0+z)/r0)**2*t*dlogp
     z101_old(i)=z101_old(i-1)+dz
     dzhalf=dz/2.0
  end do
! bottom(1) --> top height(101) arrangements for z101_old

!------------------------------------------------------------------------------------
! Redefine P,T,H and gas profiles according to surface pressure (D.Hafffner method).
!------------------------------------------------------------------------------------

  psfc2 = surface_pressure 
  idx = 1
  do i=1, nplev
     if (psfc2 .ge. p101(i)) then
         exit
     else
        idx = i
     endif
  enddo
!  write(*,*) 'idx, p101(i), psfc2, p101(i+1)', idx, p101(idx), psfc2, p101(idx+1)
!  pause

  frac = (DLOG(psfc2)-DLOG(p101(idx+1))) / (DLOG(p101(idx))-DLOG(p101(idx+1)))
  z = (z101_old(idx+1) - z101_old(idx))*0.5d0
  dz = - aconst * t101(idx) * ((r0+z)/r0)**2.0d0*dlog(psfc2/p101(idx))
  z101_old(idx) = z101_old(idx) + dz
  !
  tprf(idx) = tprf(idx)*(1.0d0-frac)+tprf(idx+1)*frac
  t101(idx) = tprf(idx)
  p101(idx) = psfc2

!write(*,*) p101
!write(*,*) frac, z101_old(idx), dz, z101_old(idx-1) 
!write(*,*) tprf

! top(0) --> bottom height(101) arrangements for height_grid

!   Continuation point

   if ( .not. New_treatment ) go to 456

!  NEW TREATMENT, cuts off at SURFACE PRESSURE
!  ===========================================

  p101log(idx) = log(psfc2)
  nlayers = nplev - idx + 1
  toa_height = 2*z101_old(nplev) - z101_old(nplev-1)
  level_heights(0)  = toa_height
  level_logpress(1) = -50.0d0
  do i = 1, nlayers
     level_heights(i)    = z101_old(nplev + 1 - i)
     level_logpress(i+1) = p101log(nplev + 1 - i)
     layer_temperatures(i) = tprf (nplev + 1 - i)
  enddo

!  Air density, molecules/cm^2
!  ---------------------------

  rho_stp  = 2.68675D+24 * 273.15d0 / 1013.25D0
  dens1 = 0.0d0
  dens2 = p101(nplev) *  rho_stp / t101(nplev)
  hdiff = toa_height - z101_old(nplev)
  layer_aircolumns(1) = hdiff * 0.5d0 * ( dens1 + dens2 )
  dens1 = dens2
  do i = nplev - 1, idx, -1
     dens2 = p101(i) *  rho_stp / t101(i)
     hdiff = z101_old(i+1) - z101_old(i)
     layer_aircolumns(nplev + 1 - i) = hdiff * 0.5d0 * ( dens1 + dens2 )
     dens1 = dens2
  enddo

!  O2O2 density, molecules^2/cm^5
!  ------------------------------

  xo2_kb  = xo2 * 1.0D-4 / kb
  dens1 = 0.0d0
  dens2 = p101(nplev) * xo2_kb / t101(nplev)
  hdiff = toa_height - z101_old(nplev)
  layer_o2o2columns(1) = hdiff * 0.5d0 * ( dens1*dens1 + dens2*dens2 ) * 1.0D+5
  dens1 = dens2
  do i = nplev - 1, idx, -1
     dens2 = p101(i) *  xo2_kb / t101(i)
     hdiff = z101_old(i+1) - z101_old(i)
     layer_o2o2columns(nplev + 1 - i) = hdiff * 0.5d0 * ( dens1*dens1 + dens2*dens2 ) * 1.0D+5
     dens1 = dens2
  enddo

!  DEBUG

!  do i = 1, nlayers
!     if (i.gt.90)write(*,*)i,level_heights(i),layer_aircolumns(i),exp(level_logpress(i+1)),layer_temperatures(i)
!  enddo

!  Done - No gas here   !!!!

   return

!   pause'alt 101 profiles'

!  OLD TREATMENT, KEEPS ALL LAYERS TO GROUND
!  =========================================

456 continue

!  height_grid
!  -----------

  nlayers = nplev
  toa_height = 2*z101_old(nplev) - z101_old(nplev-1)
  level_heights(0)  = toa_height
  level_logpress(1) = -50.0d0
  do i = 1, nplev
     level_heights(i)    = z101_old(nplev + 1 - i)
     level_logpress(i+1) = p101log(nplev + 1 - i)
  enddo

!  temperatures
!  ------------

  do i = 1, nplev
    layer_temperatures(i) = tprf (nplev + 1 - i)
  enddo

!  Air density, molecules/cm^2
!  ---------------------------

  rho_stp  = 2.68675D+24 * 273.15d0 / 1013.25D0
  dens1 = 0.0d0
  dens2 = p101(nplev) *  rho_stp / t101(nplev)
  hdiff = toa_height - z101_old(nplev)
  layer_aircolumns(1) = hdiff * 0.5d0 * ( dens1 + dens2 )
  dens1 = dens2
  do i = nplev-1, 1, -1
     dens2 = p101(i) *  rho_stp / t101(i)
     hdiff = z101_old(i+1) - z101_old(i)
     layer_aircolumns(nplev+1 - i) = hdiff * 0.5d0 * ( dens1 + dens2 )
     dens1 = dens2
  enddo

!  O2O2 density, molecules^2/cm^5
!  ------------------------------

  xo2_kb  = xo2 * 1.0D-4 / kb
  dens1 = 0.0d0
  dens2 = p101(nplev) * xo2_kb / t101(nplev)
  hdiff = toa_height - z101_old(nplev)
  layer_o2o2columns(1) = hdiff * 0.5d0 * ( dens1*dens1 + dens2*dens2 ) * 1.0D+5
  dens1 = dens2
  do i = nplev - 1, 1, -1
     dens2 = p101(i) *  xo2_kb / t101(i)
     hdiff = z101_old(i+1) - z101_old(i)
     layer_o2o2columns(nplev + 1 - i) = hdiff * 0.5d0 * ( dens1*dens1 + dens2*dens2 ) * 1.0D+5
     dens1 = dens2
  enddo

  do i = 1, nlayers
     if (i.gt.90)write(*,*)i,level_heights(i),layer_aircolumns(i),exp(level_logpress(i+1)),layer_temperatures(i)
  enddo

  do i = 0, nlayers
     level_press(i) = exp(level_logpress(i+1))
  enddo

!  Layer gas columns for ozone
!  ---------------------------

!   Use units "as is"

  layer_gascolumns(1) = xprf(nplev)
  do i = 2, nplev
     layer_gascolumns(i) = xprf(nplev+1-i)-xprf(nplev+2-i)
  enddo
  total_column = xprf(1)

!   profile derivatives
  if ( do_gas_jacobians .and. do_profile_linearization ) then
    layer_gasderivs(1)  = 1.0d0
    do i = 2, nplev
     layer_gasderivs(i)  = 1.0d0
    enddo
  endif

! column derivatives
  if ( do_gas_jacobians .and. do_column_linearization ) then
    layer_gasderivs(1)  = xprf(nplev) / xprf(1) / 1000.0d0
    do i = 2, nplev
     layer_gasderivs(i)  = layer_gascolumns(i) / xprf(1) / 1000.0d0
    enddo
  endif

!  Debug

!   do i = 1, nplev
!     write(*,*)i,layer_gascolumns(i)
!   enddo

!  Finish
!  ------

return  
end subroutine pthg_profiles_z101

subroutine develop_tomsv8_o3pth                                             &
              ( ACTUAL_O3COLUMN, LATITUDE, YEAR, MONTH, DAY_OF_MONTH,       & ! Input
                MAXLAYERS, MAXLEVELS, surface_press, preserve_gas,          & ! Input
                usa_nlevels, usa_heights, usa_logpressure, usa_temperature, & ! Input
                nlayers, level_heights, level_press, level_logpress,        & ! Output
                layer_aircolumns, layer_o2o2columns, layer_temperatures,    & ! Output
                o3umkehrs, dc_o3umkehrs,                                    & ! Output
                FAIL,  MESSAGE )

      implicit none

!     Input arguments
!     ---------------

!  Actual O3 column

      REAL(fpk)        :: ACTUAL_O3COLUMN

!     Latitude and month

      REAL(fpk)        :: LATITUDE
      INTEGER          :: YEAR, MONTH, DAY_OF_MONTH

!  Dimensioning

      INTEGER          :: MAXLAYERS, MAXLEVELS

! USA profiles

    integer                     :: USA_NLEVELS
    real(fpk)   , dimension(46) :: USA_HEIGHTS     (46)
    real(fpk)   , dimension(46) :: USA_LOGPRESSURE (46)
    real(fpk)   , dimension(46) :: USA_TEMPERATURE (46)

!  surface pressure input

      REAL(fpk)        :: SURFACE_PRESS

!  Gas preservation flag
!    (If set, uses Sigma coordinates to preserve total amount of O3)

     LOGICAL           :: preserve_gas

!  outputs
!  -------

! Trace gas Layer column amounts and cross-sections
!   Units should be consistent, e.g:
!      X-sections  are in cm^2/mol  or inverse [DU]
!      gas columns are in mol/cm^2  or [DU]

!  Atmospheric grid output

    integer                                              :: NLAYERS
    REAL(fpk)   , DIMENSION( MAXLAYERS )                 :: LAYER_AIRCOLUMNS
    REAL(fpk)   , DIMENSION( MAXLAYERS )                 :: LAYER_O2O2COLUMNS
    REAL(fpk)   , DIMENSION( MAXLAYERS )                 :: LAYER_TEMPERATURES
    REAL(fpk)   , DIMENSION( 0:MAXLAYERS )               :: LEVEL_HEIGHTS
    REAL(fpk)   , DIMENSION( 0:MAXLAYERS )               :: LEVEL_PRESS
    REAL(fpk)   , DIMENSION( MAXLEVELS   )               :: LEVEL_LOGPRESS

!  gas column amounts

      REAL(fpk)        :: O3UMKEHRS    ( MAXLAYERS )
      REAL(fpk)        :: DC_O3UMKEHRS ( MAXLAYERS )

!  status

      LOGICAL          :: FAIL
      CHARACTER*(*)    :: MESSAGE

!     Local variables
!     ---------------

!  conversion factor

      real(fpk)   , parameter :: du_to_cm2 = 2.66868D+16

!     output mapping coefficients for clear-sky cases
!     --> profile sets which define the mapping domain
!     reference columns spanning the mapping domain

      REAL(fpk)       :: CUMO3_MAPPING(14,10)
      INTEGER         :: O3MAPPING_N_REFCOLS
      REAL(fpk)       :: O3MAPPING_REFCOLS (10)

!  Interpolation variables

      integer         :: N1, N2, NLEVELS, K1, K, KC
      REAL(fpk)       :: df1, df2, gradn, cgrad, cgradp, diff, xfac, local, sum
      real(fpk)       :: psv, pmid, rho_1, rho_2, rho_a, temp, prof1, prof2, xo2_kb, tp1, tp2
      real(fpk)       :: level_heights1(201),level_temperatures(201), constant, redfac(maxlayers)

!  Time/latitude help variables

      LOGICAL         ::  LOOP, DO_INTERP
      REAL(spk)       ::  SLOPE, DIST, O3(11,10,2)
      INTEGER         ::  N, ZONE, M, REF, TIME, NBDAY(12)
      DATA NBDAY /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      REAL(fpk)       ::  COL, DEL, PLAT1, PSURF, PVAC, MULT
      REAL(fpk)       ::  PROFILE (14), CUMO3_LNPRESS(14),LATZONES(19), SIGMA_TOMS(14)
      DATA LATZONES /     -90.0, -80.0, -70.0, -60.0, -50.0, -40.0, -30.0, -20.0, -10.0, 0.0, &
                           10.0,  20.0,  30.0,  40.0, 50.0, 60.0,  70.0,  80.0,  90.0 /

!  output TOMS version 8 profiles (New versions)

      logical         ::  fail_open, fail_read
      REAL(spk)       ::  TV8O3NEW_DATA_PROFILES   (12,18,11,10)
      REAL(spk)       ::  TV8O3NEW_DATA_COVMATRIX  (11, 11)
      REAL(spk)       ::  TV8O3NEW_DATA_COLUMNS    (12,18,10)
      INTEGER         ::  TV8O3NEW_DATA_COLNUMBERS (12,18)
      REAL(spk)       ::  TV8O3NEW_DATA_PRESSURES  (14)

!  Loschmidt's number (particles/cm3), STP values

    real(fpk)   , PARAMETER :: RHO_STANDARD = 2.68675D+19
    real(fpk)   , PARAMETER :: PZERO = 1013.25D0
    real(fpk)   , PARAMETER :: TZERO = 273.15D0
    real(fpk)   , PARAMETER :: xo2 = 0.209476_fpk, kb = 1.381e-23_fpk


!  First, read the data

      CALL extract_tomsv8_o3                &
               ( TV8O3NEW_DATA_PROFILES,    &
                 TV8O3NEW_DATA_COVMATRIX,   &
                 TV8O3NEW_DATA_COLUMNS,     &
                 TV8O3NEW_DATA_COLNUMBERS,  &
                 TV8O3NEW_DATA_PRESSURES,   &
                 FAIL_OPEN, FAIL_READ )

!  Exception handling

      FAIL = FAIL_OPEN .OR. FAIL_READ
      IF ( FAIL ) THEN
         message = ' extract TOMSV8 failed - abort'
         return
      endif

!  Next, do the time/geography interpolation

      DO_INTERP=.TRUE.
         
         IF (DO_INTERP) THEN
            IF (DAY_OF_MONTH .LE. 15) MONTH=MONTH-1
            IF (MONTH .EQ. 0) MONTH=12
            IF (YEAR .EQ. 1996 .OR. YEAR .EQ. 2000 .OR. YEAR .EQ. 2004 &
                .OR. YEAR .EQ. 2008 .OR. YEAR .EQ. 2012 .OR. &
                YEAR .EQ. 2016) NBDAY(2)=29
            DO TIME=1,2
               IF (TIME .EQ. 2) MONTH=MONTH+1
               IF (MONTH .EQ. 13) MONTH=1

               IF (LATITUDE .LE. -85) THEN
                  DIST=LATITUDE+85
                  O3MAPPING_N_REFCOLS=TV8O3NEW_DATA_COLNUMBERS(MONTH,1)                    
                  DO N = 1, O3MAPPING_N_REFCOLS
                     COL=TV8O3NEW_DATA_COLUMNS(MONTH,1,N)
                     M=1
                     DO WHILE(TV8O3NEW_DATA_COLUMNS(MONTH,2,M)/=COL.AND. &
                         M .LT. TV8O3NEW_DATA_COLNUMBERS(MONTH,2) )
                        M = M+1
                     ENDDO
                     
                     DO K=1,11
                        IF (TV8O3NEW_DATA_COLUMNS(MONTH,2,M)==COL) THEN 
                           SLOPE=(TV8O3NEW_DATA_PROFILES(MONTH,2,K,M) - &
                               TV8O3NEW_DATA_PROFILES(MONTH,1,K,N))/10
                           O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,1,K,N)+ SLOPE*DIST
                        ELSE
                           O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,1,K,N)
                        ENDIF
                     ENDDO              
                  ENDDO
                  
               ELSEIF (LATITUDE .GT. 85) THEN
                  DIST=LATITUDE-85
                  O3MAPPING_N_REFCOLS=TV8O3NEW_DATA_COLNUMBERS(MONTH,18)
                  DO N = 1, O3MAPPING_N_REFCOLS
                     COL=TV8O3NEW_DATA_COLUMNS(MONTH,18,N)
                     M=1
                     DO WHILE(TV8O3NEW_DATA_COLUMNS(MONTH,17,M)/=COL  .AND. &
                         M .LT. TV8O3NEW_DATA_COLNUMBERS(MONTH,17))
                        M = M+1
                     ENDDO
                     
                     DO K=1,11
                        IF (TV8O3NEW_DATA_COLUMNS(MONTH,17,M)==COL) THEN 
                           SLOPE=(TV8O3NEW_DATA_PROFILES(MONTH,18,K,N) - &
                               TV8O3NEW_DATA_PROFILES(MONTH,17,K,M))/10
                           O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,18,K,N)+ SLOPE*DIST
                        ELSE
                           O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,18,K,N)
                        ENDIF
                     ENDDO              
                  ENDDO                 
               ELSE             !IF LATITUDE   
                  ZONE = 0
                  LOOP = .TRUE.
                  DO WHILE (LOOP)
                     ZONE = ZONE + 1
                     IF ( LATITUDE .GT. LATZONES(ZONE)+5 .AND. LATITUDE .LE. LATZONES(ZONE+1)+5 )LOOP = .FALSE.
                  ENDDO
                  DIST=LATITUDE-(LATZONES(ZONE)+5)                 

                  IF( TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE) >= &
                      TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1)) THEN
                     REF=ZONE
                     O3MAPPING_N_REFCOLS =TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE)
                     
                     DO N = 1, O3MAPPING_N_REFCOLS
                        COL=TV8O3NEW_DATA_COLUMNS(MONTH,REF,N)
                        M=1              
                        DO WHILE( TV8O3NEW_DATA_COLUMNS(MONTH,REF+1,M)/=COL  &
                           .AND. M .LT. TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1) )
                           M = M+1
                        ENDDO
                        
                        DO K=1,11
                           IF (TV8O3NEW_DATA_COLUMNS(MONTH,REF+1,M) ==COL) THEN 
                              SLOPE= (TV8O3NEW_DATA_PROFILES(MONTH,REF+1,K,M) - &
                                      TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N))/10
                              O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N)  + SLOPE*DIST
                           ELSE
                              O3(K,N,TIME)=TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N)
                           ENDIF
                        ENDDO
                     ENDDO
                  ELSE
                     REF=ZONE+1
                     O3MAPPING_N_REFCOLS= TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1)
                     
                     DO N = 1, O3MAPPING_N_REFCOLS
                        COL=TV8O3NEW_DATA_COLUMNS(MONTH,REF,N)
                        M=1              
                        DO WHILE(TV8O3NEW_DATA_COLUMNS(MONTH,REF-1,M)/=COL  & 
                            .AND. M .LT. TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE-1))
                           M = M+1                  
                        ENDDO  

                        DO K=1,11
                           IF (TV8O3NEW_DATA_COLUMNS(MONTH,REF-1,M) ==COL) THEN 
                              SLOPE=(TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N) - &
                                    TV8O3NEW_DATA_PROFILES(MONTH,REF-1,K,M))/10
                              O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES (MONTH,REF-1,K,M) + SLOPE*DIST
                           ELSE
                              O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N)
                           ENDIF
                        ENDDO

                     ENDDO
                  ENDIF
                  
               ENDIF            !ENDIF LATITUDE
            ENDDO               ! END LOOP TIME

            ZONE=1
            IF (DAY_OF_MONTH .GT. 15) MONTH=MONTH-1
            IF (MONTH .EQ. 0) MONTH=12
            DO N=1,O3MAPPING_N_REFCOLS
               DO K=1,11
                  IF (DAY_OF_MONTH .LE. 15) THEN
                     IF (MONTH .EQ. 1) THEN
                        TV8O3NEW_DATA_PROFILES(MONTH,1,K,N) = O3(K,N,1)+ &
                            (O3(K,N,2)-O3(K,N,1))/NBDAY(12)*             &
                            (NBDAY(12)-15+DAY_OF_MONTH)
                     ELSE
                        TV8O3NEW_DATA_PROFILES(MONTH,1,K,N) = O3(K,N,1)+ &
                            (O3(K,N,2)-O3(K,N,1))/NBDAY(MONTH-1)*        &
                            (NBDAY(MONTH-1)-15+DAY_OF_MONTH)
                     ENDIF
                  ELSE
                     TV8O3NEW_DATA_PROFILES(MONTH,1,K,N) = O3(K,N,1)+    &
                         (O3(K,N,2)-O3(K,N,1))/                          &
                         NBDAY(MONTH)*(DAY_OF_MONTH-15)
                  ENDIF
               ENDDO
            ENDDO

         ELSE                   ! DO_INTERP

!     get the latitude zone

            IF ( LATITUDE .LE. -80.0D0 ) THEN
               ZONE = 1
            ELSE IF ( LATITUDE .GT. +80.0D0 ) THEN
               ZONE = 18
            ELSE
               ZONE = 1
               LOOP = .TRUE.
               DO WHILE (LOOP)
                  ZONE = ZONE + 1
                  IF ( LATITUDE .GT. LATZONES(ZONE) .AND. &
                       LATITUDE .LE. LATZONES(ZONE+1) ) LOOP = .FALSE.
               END DO
            ENDIF

!     Get the reference columns

            O3MAPPING_N_REFCOLS = TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE)

         ENDIF                  !  DO_INTERP

!  Get the basic profiles
!  ----------------------

!     For each of the column set:

      DO N = 1, O3MAPPING_N_REFCOLS

!     create TOMS V8 profile, parcel out top 3 layers

        DO K = 1, 11
          PLAT1 = DBLE ( TV8O3NEW_DATA_PROFILES(MONTH,ZONE,K,N) )
          PROFILE(K+2) = PLAT1
        ENDDO
        DEL = PROFILE(3)
        PROFILE(1) = 1.0D0 * DEL / 16.0D0
        PROFILE(2) = 3.0d0 * DEL / 16.0D0
        PROFILE(3) = 3.0d0 * DEL / 4.0D0

!     set profile mapping, cumulative ozone

        CUMO3_MAPPING(1,N) = 0.0d0
        DO K = 2, 14
          CUMO3_MAPPING(K,N) = CUMO3_MAPPING(K-1,N)+PROFILE(K-1)
        ENDDO

!     Set the reference columns to clear sky columns (MVR, 20.3.2006)

        O3MAPPING_REFCOLS(N) = CUMO3_MAPPING(14,N)

      ENDDO

!  Preserve ozone
!  ==============

!  Skip if not flagged

      if ( .not. preserve_gas ) go to 633

!  set Log pressure grid

      DO K = 1, 14
        CUMO3_LNPRESS(K)     = DLOG(DBLE(TV8O3NEW_DATA_PRESSURES(K)))
      ENDDO
      CUMO3_LNPRESS(14) = DLOG(SURFACE_PRESS)

!  TOMS Sigma coordinates

      PSURF = EXP(DLOG(DBLE(TV8O3NEW_DATA_PRESSURES(14))))
      PVAC  = EXP(DLOG(DBLE(TV8O3NEW_DATA_PRESSURES(1))))
      MULT = 1.0d0 / ( PSURF - PVAC )
      SIGMA_TOMS(1) = 0.0d0
      DO K = 2, 13
         SIGMA_TOMS(K) = ( EXP(CUMO3_LNPRESS(K)) - PVAC ) * MULT
      ENDDO
      SIGMA_TOMS(14) = 1.0d0

!  Redefine the Logpressure grid

      PSURF = SURFACE_PRESS
      DO K = 2, 14
         CUMO3_LNPRESS(K) = DLOG ( SIGMA_TOMS(K) * PSURF  + PVAC * ( 1.0d0 - SIGMA_TOMS(K) ) )
      ENDDO

!  debug Sigma-scaling

!      DO K = 1, 14
!        write(*,*)'new',K,CUMO3_LNPRESS(K),EXP(CUMO3_LNPRESS(K))
!      ENDDO
!      pause'sigma_toms'

!  set the pressures and number of layers

      nlayers = 13 ; nlevels = 14
      DO K = 1, NLAYERS + 1
         level_Logpress(k) = CUMO3_LNPRESS(K)
         level_press(k-1)  = dexp(CUMO3_LNPRESS(K))
      enddo

!  Set the ozone, by Interpolating TOMSV8 Ozone to our incoming column

      CALL LINEAR_TRAWLER_ASCENDING                 &
            ( 10, O3MAPPING_N_REFCOLS,              &
              O3MAPPING_REFCOLS, ACTUAL_O3COLUMN,   &
              DF1, DF2, N1, N2 )

!  Set profile

      DIFF = O3MAPPING_REFCOLS(N2) - O3MAPPING_REFCOLS(N1)
      XFAC = ACTUAL_O3COLUMN - O3MAPPING_REFCOLS(N1)
      DO K = 1, NLAYERS
         K1 = K + 1
         Prof1 = CUMO3_MAPPING(K1,N1) - CUMO3_MAPPING(K,N1)
         Prof2 = CUMO3_MAPPING(K1,N2) - CUMO3_MAPPING(K,N2)
         CGRAD = ( Prof2 - Prof1 ) / Diff
         o3umkehrs(k)    = Prof1 + XFAC * cgrad
         dc_o3umkehrs(k) = cgrad
         if ( dabs(cgrad).lt.1.0d-12)dc_o3umkehrs(k) =0.0d0
      ENDDO

!  continuation point

      go to 634

!  Non-preserved ozone
!  ===================

!  continuation point

633   continue

!  set Log pressure grid

      KC = 0
      DO K = 1, 14
        TP1 = DBLE(TV8O3NEW_DATA_PRESSURES(K))
        if ( SURFACE_PRESS .ge. TP1 ) then
           KC = KC + 1 ; CUMO3_LNPRESS(KC) = DLOG(TP1)
        endif
      ENDDO
      NLEVELS = KC + 1 ; if (KC.eq.14) NLEVELS = KC
      CUMO3_LNPRESS(NLEVELS) = DLOG(SURFACE_PRESS)
      NLAYERS = NLEVELS - 1

      DO K = 1, NLEVELS
         level_Logpress(k) = CUMO3_LNPRESS(K)
         level_press(k-1)  = dexp(CUMO3_LNPRESS(K))
      enddo

!  Reduction factor for profiles

      REDFAC = 1.0_fpk
      TP1 = DBLE(TV8O3NEW_DATA_PRESSURES(NLAYERS))
      TP2 = DBLE(TV8O3NEW_DATA_PRESSURES(NLEVELS))
      REDFAC(NLAYERS) =  ( SURFACE_PRESS - TP1 ) / ( TP2 - TP1 )

!  Set the ozone, by Interpolating TOMSV8 Ozone to our incoming column

      CALL LINEAR_TRAWLER_ASCENDING                 &
            ( 10, O3MAPPING_N_REFCOLS,              &
              O3MAPPING_REFCOLS, ACTUAL_O3COLUMN,   &
              DF1, DF2, N1, N2 )

!  Set profile

      DIFF = O3MAPPING_REFCOLS(N2) - O3MAPPING_REFCOLS(N1)
      XFAC = ACTUAL_O3COLUMN - O3MAPPING_REFCOLS(N1)
      DO K = 1, NLAYERS
         K1 = K + 1
         Prof1 = CUMO3_MAPPING(K1,N1) - CUMO3_MAPPING(K,N1)
         Prof2 = CUMO3_MAPPING(K1,N2) - CUMO3_MAPPING(K,N2)
         CGRAD = REDFAC(k) * ( Prof2 - Prof1 ) / Diff
         o3umkehrs(k)    = REDFAC(k) * Prof1 + XFAC * cgrad
         dc_o3umkehrs(k) = cgrad
         if ( dabs(cgrad).lt.1.0d-12)dc_o3umkehrs(k) =0.0d0
      ENDDO

!  Convert Profile to mol/cm/cm
!  ============================

!  continuation point

634   continue

!    Want gradient to be un-normalized

      cgrad = 0.0d0
      do n = 1, nlayers
         cgrad = cgrad + o3umkehrs(n)
         o3umkehrs(n)    =  o3umkehrs(n)    * du_to_cm2
         dc_o3umkehrs(n) =  dc_o3umkehrs(n) * du_to_cm2
      enddo
!     write(*,*)ACTUAL_O3COLUMN,cgrad ; pause

!  set the temperatures, Air density and O2O2 density
!  ==================================================

!  Interpolate USA height and temperature to the TOMS grid

      call lintp2 &
          (usa_nlevels,usa_logpressure,usa_temperature,&
           nlevels, level_logpress, level_temperatures )

      call lintp2 &
          (usa_nlevels,usa_logpressure,usa_heights,&
           nlevels, level_logpress, level_heights1 )

!  gas constant

      CONSTANT = 1.0D+05 * RHO_STANDARD * TZERO / PZERO

!  Develop layer temperatures, layer air densities, assign height grid

      LEVEL_HEIGHTS(0) = LEVEL_HEIGHTS1(1)
      DO K = 1, NLAYERS
        RHO_1 = LEVEL_PRESS(K-1)   / LEVEL_TEMPERATURES(K)
        RHO_2 = LEVEL_PRESS(K)     / LEVEL_TEMPERATURES(K+1)
        TEMP  = 0.5D0*(LEVEL_TEMPERATURES(K)+LEVEL_TEMPERATURES(K+1))
        RHO_A = 0.5D0 * CONSTANT * ( RHO_1 + RHO_2 )
        LEVEL_HEIGHTS(K) = LEVEL_HEIGHTS1(K+1)
        DIFF  = LEVEL_HEIGHTS1(K) - LEVEL_HEIGHTS1(K+1)
        LAYER_TEMPERATURES(K)   = TEMP
        LAYER_AIRCOLUMNS(K)     = DIFF * RHO_A
      ENDDO

!  Develop O2O2 densities

      xo2_kb  = xo2 * 1.0D-4 / kb
      DO K = 1, NLAYERS
        RHO_1 = xo2_kb * LEVEL_PRESS(K-1)   / LEVEL_TEMPERATURES(K)
        RHO_2 = xo2_kb * LEVEL_PRESS(K)     / LEVEL_TEMPERATURES(K+1)
        RHO_A = 0.5D0 * ( RHO_1 * RHO_1 + RHO_2 * RHO_2 )
        DIFF  = LEVEL_HEIGHTS1(K) - LEVEL_HEIGHTS1(K+1)
        LAYER_O2O2COLUMNS(K)     = DIFF * RHO_A * 1.0E+5_fpk
      ENDDO

!  Finish

      return
end subroutine develop_tomsv8_o3pth

subroutine develop_tomsv8_o3                                     &
         ( ACTUAL_O3COLUMN, LATITUDE, YEAR, MONTH, DAY_OF_MONTH, & ! Input
           MAXLAYERS, MAXLEVELS, preserve_gas, NLAYERS,          & ! Input
           level_logpressures, surface_press,                    & ! Input
           o3umkehrs, dc_o3umkehrs,                              & ! Output
           FAIL,  MESSAGE )

      implicit none

!     Input arguments
!     ---------------

!  Actual O3 column

      REAL(fpk)        :: ACTUAL_O3COLUMN

!     Latitude and month

      REAL(fpk)        :: LATITUDE
      INTEGER          :: YEAR, MONTH, DAY_OF_MONTH

!  Dimensioning

      INTEGER          :: MAXLAYERS, MAXLEVELS

!  Gas preservation flag
!    (If set, uses Sigma coordinates to preserve total amount of O3)

     LOGICAL           :: preserve_gas

!  pressure grid input

      INTEGER          :: NLAYERS
      REAL(fpk)        :: LEVEL_LOGPRESSURES    ( MAXLEVELS )
      REAL(fpk)        :: SURFACE_PRESS

!  outputs
!  -------

!  gas column amounts

      REAL(fpk)        :: O3UMKEHRS    ( MAXLAYERS )
      REAL(fpk)        :: DC_O3UMKEHRS ( MAXLAYERS )

!  status

      LOGICAL          :: FAIL
      CHARACTER*(*)    :: MESSAGE

!     Local variables
!     ---------------

!  conversion factor

      real(fpk)   , parameter :: du_to_cm2 = 2.66868D+16

!     output mapping coefficients for clear-sky cases
!     --> profile sets which define the mapping domain
!     reference columns spanning the mapping domain

      REAL(fpk)       :: CUMO3_MAPPING(14,10)
      INTEGER         :: O3MAPPING_N_REFCOLS
      REAL(fpk)       :: O3MAPPING_REFCOLS (10)

!  Interpolation variables

      integer         :: N1, N2, NLEVELS
      REAL(fpk)       :: df1, df2, gradn, cgrad, cgradp, diff, xfac, local
      REAL(fpk)       :: local_o31(140),local_o32(140), sum

!  Time/latitude help variables

      LOGICAL         ::  LOOP, DO_INTERP
      REAL(spk)       ::  SLOPE, DIST, O3(11,10,2)
      INTEGER         ::  K, N, ZONE, M, REF, TIME,NBDAY(12)
      DATA NBDAY /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      REAL(fpk)       ::  COL, DEL, PLAT1, PSURF, PVAC, MULT
      REAL(fpk)       ::  PROFILE (14), CUMO3_LNPRESS(14),LATZONES(19), SIGMA_TOMS(14)
      DATA LATZONES /     -90.0, -80.0, -70.0, -60.0, -50.0, -40.0, -30.0, -20.0, -10.0, 0.0, &
                           10.0,  20.0,  30.0,  40.0, 50.0, 60.0,  70.0,  80.0,  90.0 /

!  output TOMS version 8 profiles (New versions)

      logical         ::  fail_open, fail_read
      REAL(spk)       ::  TV8O3NEW_DATA_PROFILES   (12,18,11,10)
      REAL(spk)       ::  TV8O3NEW_DATA_COVMATRIX  (11, 11)
      REAL(spk)       ::  TV8O3NEW_DATA_COLUMNS    (12,18,10)
      INTEGER         ::  TV8O3NEW_DATA_COLNUMBERS (12,18)
      REAL(spk)       ::  TV8O3NEW_DATA_PRESSURES  (14)

!  First, read the data

      CALL extract_tomsv8_o3                &
               ( TV8O3NEW_DATA_PROFILES,    &
                 TV8O3NEW_DATA_COVMATRIX,   &
                 TV8O3NEW_DATA_COLUMNS,     &
                 TV8O3NEW_DATA_COLNUMBERS,  &
                 TV8O3NEW_DATA_PRESSURES,   &
                 FAIL_OPEN, FAIL_READ )

!  Exception handling

      FAIL = FAIL_OPEN .OR. FAIL_READ
      IF ( FAIL ) THEN
         message = ' extract TOMSV8 failed - abort'
         return
      endif

!  Next, do the time/geography interpolation

      DO_INTERP=.TRUE.
         
         IF (DO_INTERP) THEN
            IF (DAY_OF_MONTH .LE. 15) MONTH=MONTH-1
            IF (MONTH .EQ. 0) MONTH=12
            IF (YEAR .EQ. 1996 .OR. YEAR .EQ. 2000 .OR. YEAR .EQ. 2004 &
                .OR. YEAR .EQ. 2008 .OR. YEAR .EQ. 2012 .OR. &
                YEAR .EQ. 2016) NBDAY(2)=29
            DO TIME=1,2
               IF (TIME .EQ. 2) MONTH=MONTH+1
               IF (MONTH .EQ. 13) MONTH=1

               IF (LATITUDE .LE. -85) THEN
                  DIST=LATITUDE+85
                  O3MAPPING_N_REFCOLS=TV8O3NEW_DATA_COLNUMBERS(MONTH,1)                    
                  DO N = 1, O3MAPPING_N_REFCOLS
                     COL=TV8O3NEW_DATA_COLUMNS(MONTH,1,N)
                     M=1
                     DO WHILE(TV8O3NEW_DATA_COLUMNS(MONTH,2,M)/=COL.AND. &
                         M .LT. TV8O3NEW_DATA_COLNUMBERS(MONTH,2) )
                        M = M+1
                     ENDDO
                     
                     DO K=1,11
                        IF (TV8O3NEW_DATA_COLUMNS(MONTH,2,M)==COL) THEN 
                           SLOPE=(TV8O3NEW_DATA_PROFILES(MONTH,2,K,M) - &
                               TV8O3NEW_DATA_PROFILES(MONTH,1,K,N))/10
                           O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,1,K,N)+ SLOPE*DIST
                        ELSE
                           O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,1,K,N)
                        ENDIF
                     ENDDO              
                  ENDDO
                  
               ELSEIF (LATITUDE .GT. 85) THEN
                  DIST=LATITUDE-85
                  O3MAPPING_N_REFCOLS=TV8O3NEW_DATA_COLNUMBERS(MONTH,18)
                  DO N = 1, O3MAPPING_N_REFCOLS
                     COL=TV8O3NEW_DATA_COLUMNS(MONTH,18,N)
                     M=1
                     DO WHILE(TV8O3NEW_DATA_COLUMNS(MONTH,17,M)/=COL  .AND. &
                         M .LT. TV8O3NEW_DATA_COLNUMBERS(MONTH,17))
                        M = M+1
                     ENDDO
                     
                     DO K=1,11
                        IF (TV8O3NEW_DATA_COLUMNS(MONTH,17,M)==COL) THEN 
                           SLOPE=(TV8O3NEW_DATA_PROFILES(MONTH,18,K,N) - &
                               TV8O3NEW_DATA_PROFILES(MONTH,17,K,M))/10
                           O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,18,K,N)+ SLOPE*DIST
                        ELSE
                           O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,18,K,N)
                        ENDIF
                     ENDDO              
                  ENDDO                 
               ELSE             !IF LATITUDE   
                  ZONE = 0
                  LOOP = .TRUE.
                  DO WHILE (LOOP)
                     ZONE = ZONE + 1
                     IF ( LATITUDE .GT. LATZONES(ZONE)+5 .AND. LATITUDE .LE. LATZONES(ZONE+1)+5 )LOOP = .FALSE.
                  ENDDO
                  DIST=LATITUDE-(LATZONES(ZONE)+5)                 

                  IF( TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE) >= &
                      TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1)) THEN
                     REF=ZONE
                     O3MAPPING_N_REFCOLS =TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE)
                     
                     DO N = 1, O3MAPPING_N_REFCOLS
                        COL=TV8O3NEW_DATA_COLUMNS(MONTH,REF,N)
                        M=1              
                        DO WHILE( TV8O3NEW_DATA_COLUMNS(MONTH,REF+1,M)/=COL  &
                           .AND. M .LT. TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1) )
                           M = M+1
                        ENDDO
                        
                        DO K=1,11
                           IF (TV8O3NEW_DATA_COLUMNS(MONTH,REF+1,M) ==COL) THEN 
                              SLOPE= (TV8O3NEW_DATA_PROFILES(MONTH,REF+1,K,M) - &
                                      TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N))/10
                              O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N)  + SLOPE*DIST
                           ELSE
                              O3(K,N,TIME)=TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N)
                           ENDIF
                        ENDDO
                     ENDDO
                  ELSE
                     REF=ZONE+1
                     O3MAPPING_N_REFCOLS= TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1)
                     
                     DO N = 1, O3MAPPING_N_REFCOLS
                        COL=TV8O3NEW_DATA_COLUMNS(MONTH,REF,N)
                        M=1              
                        DO WHILE(TV8O3NEW_DATA_COLUMNS(MONTH,REF-1,M)/=COL  & 
                            .AND. M .LT. TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE-1))
                           M = M+1                  
                        ENDDO  

                        DO K=1,11
                           IF (TV8O3NEW_DATA_COLUMNS(MONTH,REF-1,M) ==COL) THEN 
                              SLOPE=(TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N) - &
                                    TV8O3NEW_DATA_PROFILES(MONTH,REF-1,K,M))/10
                              O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES (MONTH,REF-1,K,M) + SLOPE*DIST
                           ELSE
                              O3(K,N,TIME)= TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N)
                           ENDIF
                        ENDDO

                     ENDDO
                  ENDIF
                  
               ENDIF            !ENDIF LATITUDE
            ENDDO               ! END LOOP TIME

            ZONE=1
            IF (DAY_OF_MONTH .GT. 15) MONTH=MONTH-1
            IF (MONTH .EQ. 0) MONTH=12
            DO N=1,O3MAPPING_N_REFCOLS
               DO K=1,11
                  IF (DAY_OF_MONTH .LE. 15) THEN
                     IF (MONTH .EQ. 1) THEN
                        TV8O3NEW_DATA_PROFILES(MONTH,1,K,N) = O3(K,N,1)+ &
                            (O3(K,N,2)-O3(K,N,1))/NBDAY(12)*             &
                            (NBDAY(12)-15+DAY_OF_MONTH)
                     ELSE
                        TV8O3NEW_DATA_PROFILES(MONTH,1,K,N) = O3(K,N,1)+ &
                            (O3(K,N,2)-O3(K,N,1))/NBDAY(MONTH-1)*        &
                            (NBDAY(MONTH-1)-15+DAY_OF_MONTH)
                     ENDIF
                  ELSE
                     TV8O3NEW_DATA_PROFILES(MONTH,1,K,N) = O3(K,N,1)+    &
                         (O3(K,N,2)-O3(K,N,1))/                          &
                         NBDAY(MONTH)*(DAY_OF_MONTH-15)
                  ENDIF
               ENDDO
            ENDDO

         ELSE                   ! DO_INTERP

!     get the latitude zone

            IF ( LATITUDE .LE. -80.0D0 ) THEN
               ZONE = 1
            ELSE IF ( LATITUDE .GT. +80.0D0 ) THEN
               ZONE = 18
            ELSE
               ZONE = 1
               LOOP = .TRUE.
               DO WHILE (LOOP)
                  ZONE = ZONE + 1
                  IF ( LATITUDE .GT. LATZONES(ZONE) .AND. &
                       LATITUDE .LE. LATZONES(ZONE+1) ) LOOP = .FALSE.
               END DO
            ENDIF

!     Get the reference columns

            O3MAPPING_N_REFCOLS = TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE)

         ENDIF                  !  DO_INTERP

!  Get the basic profiles
!  ----------------------

!     For each of the column set:

      DO N = 1, O3MAPPING_N_REFCOLS

!     create TOMS V8 profile, parcel out top 3 layers

        DO K = 1, 11
          PLAT1 = DBLE ( TV8O3NEW_DATA_PROFILES(MONTH,ZONE,K,N) )
          PROFILE(K+2) = PLAT1
        ENDDO
        DEL = PROFILE(3)
        PROFILE(1) = 1.0D0 * DEL / 16.0D0
        PROFILE(2) = 3.0d0 * DEL / 16.0D0
        PROFILE(3) = 3.0d0 * DEL / 4.0D0

!     set profile mapping, cumulative ozone

        CUMO3_MAPPING(1,N) = 0.0d0
        DO K = 2, 14
          CUMO3_MAPPING(K,N) = CUMO3_MAPPING(K-1,N)+PROFILE(K-1)
        ENDDO

!     Set the reference columns to clear sky columns (MVR, 20.3.2006)

        O3MAPPING_REFCOLS(N) = CUMO3_MAPPING(14,N)

      ENDDO

!  Preserve ozone
!  ==============

!  Skip if not flagged

      if ( .not. preserve_gas ) go to 633

!  set Log pressure grid

      DO K = 1, 14
        CUMO3_LNPRESS(K)     = DLOG(DBLE(TV8O3NEW_DATA_PRESSURES(K)))
      ENDDO
      CUMO3_LNPRESS(14) = DLOG(SURFACE_PRESS)

!  Make sure there's no extrapolation

      if ( CUMO3_LNPRESS(1).gt.LEVEL_LOGPRESSURES(1)) THEN
         CUMO3_LNPRESS(1) = LEVEL_LOGPRESSURES(1) - 0.000001d0
      endif

!  TOMS Sigma coordinates

      PSURF = EXP(DLOG(DBLE(TV8O3NEW_DATA_PRESSURES(14))))
      PVAC  = EXP(LEVEL_LOGPRESSURES(1)) ; MULT = 1.0d0 / ( PSURF - PVAC )
      SIGMA_TOMS(1) = 0.0d0
      DO K = 2, 13
         SIGMA_TOMS(K) = ( EXP(CUMO3_LNPRESS(K)) - PVAC ) * MULT
      ENDDO
      SIGMA_TOMS(14) = 1.0d0

!  Redefine the Logpressure grid

      PSURF = SURFACE_PRESS
      DO K = 2, 14
         CUMO3_LNPRESS(K) = DLOG ( SIGMA_TOMS(K) * PSURF  + PVAC * ( 1.0d0 - SIGMA_TOMS(K) ) )
      ENDDO

!  Set the ozone, by Interpolating TOMSV8 Ozone to our incoming column

      CALL LINEAR_TRAWLER_ASCENDING                 &
            ( 10, O3MAPPING_N_REFCOLS,              &
              O3MAPPING_REFCOLS, ACTUAL_O3COLUMN,   &
              DF1, DF2, N1, N2 )

!     write(*,*)DEXP(LEVEL_LOGPRESSURES(1))
!     do k = 1, 13
!       write(*,*)DEXP(CUMO3_LNPRESS(k+1)),CUMO3_MAPPING(k+1,N1)-CUMO3_MAPPING(k,N1),CUMO3_MAPPING(k+1,N2)-CUMO3_MAPPING(k,N2)
!     enddo

!  interpolate Cumulative O3 to incoming pressure grid

      nlevels = nlayers + 1
      CALL LINTP2 ( 14,CUMO3_LNPRESS,CUMO3_MAPPING(1,N1),    &
                    NLEVELS, LEVEL_LOGPRESSURES, LOCAL_O31 )
      CALL LINTP2 ( 14,CUMO3_LNPRESS,CUMO3_MAPPING(1,N2),   &
                    NLEVELS, LEVEL_LOGPRESSURES, LOCAL_O32 )

      DIFF = O3MAPPING_REFCOLS(N2) - O3MAPPING_REFCOLS(N1)
      XFAC = ACTUAL_O3COLUMN - O3MAPPING_REFCOLS(N1)
      sum = 0.0d0
      do n = 1, nlayers
        if ( n.eq.1 ) then
          cgrad =  ( LOCAL_O32(2) - LOCAL_O31(2) ) / DIFF
          o3umkehrs(n)    = LOCAL_O31(2) + XFAC * cgrad
          dc_o3umkehrs(n) = cgrad
          cgradp = cgrad
        else
          cgrad = ( LOCAL_O32(N+1) - LOCAL_O31(N+1) ) / DIFF
          gradn = cgrad - cgradp
          local = LOCAL_O31(N+1) - LOCAL_O31(N)
          o3umkehrs(n)    =  local + XFAC * GRADN
          dc_o3umkehrs(n) = GRADN
          if ( dabs(gradn).lt.1.0d-12)dc_o3umkehrs(n) =0.0d0
          cgradp = cgrad
        endif
        sum = sum + o3umkehrs(n)
      enddo

!  continuation point

      go to 634

!  Non-preserved ozone
!  ===================

!  continuation point

633   continue

!  Adjust lowest

      nlevels = nlayers + 1
      DO K = 1, 14
        CUMO3_LNPRESS(K)     = DLOG(DBLE(TV8O3NEW_DATA_PRESSURES(K)))
      ENDDO
      CUMO3_LNPRESS(14) = MAX(LEVEL_LOGPRESSURES(NLEVELS),CUMO3_LNPRESS(14))

!  Set the ozone, by Interpolating TOMSV8 Ozone to our incoming column

      CALL LINEAR_TRAWLER_ASCENDING                 &
            ( 10, O3MAPPING_N_REFCOLS,              &
              O3MAPPING_REFCOLS, ACTUAL_O3COLUMN,   &
              DF1, DF2, N1, N2 )

!  interpolate Cumulative O3 to incoming pressure grid

      nlevels = nlayers + 1
      CALL LINTP2 ( 14,CUMO3_LNPRESS,CUMO3_MAPPING(1,N1),    &
                    NLEVELS, LEVEL_LOGPRESSURES, LOCAL_O31 )
      CALL LINTP2 ( 14,CUMO3_LNPRESS,CUMO3_MAPPING(1,N2),   &
                    NLEVELS, LEVEL_LOGPRESSURES, LOCAL_O32 )

      DIFF = O3MAPPING_REFCOLS(N2) - O3MAPPING_REFCOLS(N1)
      XFAC = ACTUAL_O3COLUMN - O3MAPPING_REFCOLS(N1)
      sum = 0.0d0
      do n = 1, nlayers
        if ( n.eq.1 ) then
          cgrad =  ( LOCAL_O32(2) - LOCAL_O31(2) ) / DIFF
          o3umkehrs(n)    = LOCAL_O31(2) + XFAC * cgrad
          dc_o3umkehrs(n) = cgrad
          cgradp = cgrad
        else
          cgrad = ( LOCAL_O32(N+1) - LOCAL_O31(N+1) ) / DIFF
          gradn = cgrad - cgradp
          local = LOCAL_O31(N+1) - LOCAL_O31(N)
          o3umkehrs(n)    =  local + XFAC * GRADN
          dc_o3umkehrs(n) = GRADN
          if ( dabs(gradn).lt.1.0d-12)dc_o3umkehrs(n) =0.0d0
          cgradp = cgrad
        endif
        sum = sum + o3umkehrs(n)
      enddo

!  Convert Profile to mol/cm/cm
!  ============================

!  Continuation point

634   continue
      cgrad = 0.0d0
      do n = 1, nlayers
        cgrad = cgrad + o3umkehrs(n)
        o3umkehrs(n)    =  o3umkehrs(n)    * du_to_cm2
        dc_o3umkehrs(n) =  dc_o3umkehrs(n) * du_to_cm2
      enddo
!     write(*,*)ACTUAL_O3COLUMN,cgrad ; pause

!  Finish

      return
end subroutine develop_tomsv8_o3



SUBROUTINE extract_tomsv8_o3        &
               ( TV8O3NEW_DATA_PROFILES,    &
                 TV8O3NEW_DATA_COVMATRIX,   &
                 TV8O3NEW_DATA_COLUMNS,     &
                 TV8O3NEW_DATA_COLNUMBERS,  &
                 TV8O3NEW_DATA_PRESSURES,   &
                 FAIL_OPEN, FAIL_READ )

      implicit none

!  Read the TOMS Version 8 profile data - New data set.
!  R. Spurr, 15 december 2003
!  R. Spurr, 18 May 2009 (F90 translation)

!  output TOMS version 8 profiles (New versions)

      REAL(spk)    ::  TV8O3NEW_DATA_PROFILES   (12,18,11,10)
      REAL(spk)    ::  TV8O3NEW_DATA_COVMATRIX  (11, 11)
      REAL(spk)    ::  TV8O3NEW_DATA_COLUMNS    (12,18,10)
      INTEGER      ::  TV8O3NEW_DATA_COLNUMBERS (12,18)
      REAL(spk)    ::  TV8O3NEW_DATA_PRESSURES  (14)

!  status information

      LOGICAL      ::  FAIL_OPEN
      LOGICAL      ::  FAIL_READ

!  local

      INTEGER      ::  I, K, NC, C, M, Z, UNUM
      REAL(spk)    ::  COL, PROF(11)

!  initialise file read

      FAIL_OPEN = .FALSE.
      FAIL_READ = .FALSE.
      
!  Open files

      UNUM = 44
      OPEN ( UNIT   = UNUM, FILE   = 'physics_data/tomsv8_ozoneprofiles_new.dat', &
            ERR    = 900, STATUS = 'OLD' )
      READ (UNUM,*)
      READ (UNUM,*)
      READ (UNUM,*,ERR=901)(TV8O3NEW_DATA_PRESSURES(K),K=1,14)
      READ (UNUM,*)
      DO I = 1, 11
         READ (UNUM,*,ERR=901)(TV8O3NEW_DATA_COVMATRIX(I,K),K=1,11)
      ENDDO
      DO M = 1, 12
         DO Z = 1, 18
            READ(UNUM,*)        ! The "Month = M Latitude =  Z" lines
            NC = 0
            DO C = 1, 10
               READ(UNUM,'(F7.1,11F7.2)')COL,(PROF(K),K=11,1,-1)
               IF ( COL .NE. 999.0 ) THEN
                  NC = NC + 1
                  TV8O3NEW_DATA_COLUMNS(M,Z,NC) = COL
                  DO K = 1, 11
                     TV8O3NEW_DATA_PROFILES(M,Z,K,NC) = PROF(K)
                  ENDDO
               ENDIF
               TV8O3NEW_DATA_COLNUMBERS(M,Z) = NC
            ENDDO
         ENDDO
      ENDDO
      CLOSE (UNUM)

!  Normal return

      RETURN

!  Open-file error return

 900  CONTINUE
      FAIL_OPEN = .TRUE.
      CLOSE(UNUM)
      RETURN

!  open-read error return

 901  CONTINUE
      FAIL_READ = .TRUE.
      CLOSE(UNUM)
      RETURN

!  finish

END SUBROUTINE extract_tomsv8_o3

subroutine extract_geoschem_so2                      &
        ( latitude, longitude, month, do_first_read, & ! Input
          eta, latgc, longc, so2gc,                  & ! Output
          so2profile, fail, message )                  ! Output

      implicit none

!  Inputs

      real(fpk)        :: latitude, longitude
      integer          :: month

!  Modified input

      logical          :: do_first_read

!  Data outputs

      real(fpk)        :: so2gc(91,144,30)
      real(fpk)        :: eta(30),latgc(91),longc(144)

!  Extracted output

      real(fpk)        :: so2profile(30)

!  Exception handling

      logical          :: fail
      character*(*)    :: message

!  Local

      integer          :: k, i, j, m1, m2, n1, n2
      real(fpk)        :: qlat(2),qlon(2),lon1,lon2,lat1,lat2
      real(fpk)        :: a11,a12,a21,a22,square,terp

      character(LEN=14)     :: c14
      character(LEN=36)     :: filenames (1)
      data filenames / 'physics_data/shape_factor_200607.dat' /

!  Initialise

      fail = .false.
      message = ' '
      if ( month.ne.7) then
        message = 'Extract Geoschem SO2; must have Month = 7 (only data available)'
        fail = .true.
        return
      endif

!  read from file
      
      if (do_first_read) then
       open(1,file = filenames(1),err=990,status='old')
       read(1,'(a14,30f11.6)')c14,(eta(k),k=1,30)
       do i = 1, 91
         do j = 1, 144
           read(1,'(f6.2,f8.2,30f11.6)')latgc(i),longc(j),(so2gc(i,j,k),k=1,30)
         enddo
       enddo
       close(1)
      endif
      do_first_read = .false.
      go to 677

!  fileread error

 990  continue
      fail = .true.
      message = 'so2 goeschem shape-fectors: FILE NOT FOUND'
      return

!  Continue with extraction (bilinear interpolation)

 677  continue

!  Latitude trawling

      if ( latitude .le. latgc(1) ) then
        n1 = 1
        n2 = 2
      else if ( latitude .ge. latgc(91) ) then
        n1 = 90
        n2 = 91
      else
        n1 = 1
        do while (latitude .ge. latgc(n1) )
          n1 = n1 + 1
        enddo
        n1 = n1 - 1
        n2 = n1 + 1
      endif
      qlat(1) = latgc(n1)
      qlat(2) = latgc(n2)

!  Longitude trawling

      if ( longitude .ge. longc(144) ) then
        m1 = 144
        m2 = 1
        qlon(1) = longc(144)
        qlon(2) = 180.0d0
      else
        m1 = 1
        do while (longitude .ge. longc(m1) )
          m1 = m1 + 1
        enddo
        m1 = m1 - 1
        m2 = m1 + 1
        qlon(1) = longc(m1)
        qlon(2) = longc(m2)
      endif

!  Bilinear interpolation

      square = 1.0d0 / (qlon(2)-qlon(1)) / (qlat(2) - qlat(1))
      lon1 = longitude - qlon(1)
      lon2 = qlon(2) - longitude
      lat1 = latitude - qlat(1)
      lat2 = qlat(2) - latitude
      a11  = lon1 * lat1 * square
      a12  = lon1 * lat2 * square
      a21  = lon2 * lat1 * square
      a22  = lon2 * lat2 * square

      do k = 1, 30
        terp = a11 * so2gc(n1,m1,k) + a12 * so2gc(n1,m2,k) + &
               a21 * so2gc(n2,m1,k) + a22 * so2gc(n2,m2,k)
        so2profile(k) = terp
      enddo

!  finish

      return
end subroutine extract_geoschem_so2

subroutine develop_geoschem_pth                                          &
          ( maxlay, maxlev, surface_pressure, vacuum_pressure, etagc,    & 
            usa_nlevels, usa_heights, usa_logpressure, usa_temperature,  &
            nlayers, level_heights, level_logpress,                      &
            layer_aircolumns, layer_o2o2columns, layer_temperatures )

      implicit none

!  inputs
!  ------

!  dimensioning

   integer :: maxlay, maxlev

!  Bounding pressures

   real(fpk)    :: surface_pressure, vacuum_pressure

!  GEOSCHEM eta-coordinates:   eta = ( P - Pvac) / (Psurf - Pvac)

   real(fpk)    :: etagc(30)

! USA profiles

    integer                     :: USA_NLEVELS
    real(fpk)   , dimension(46) :: USA_HEIGHTS     (46)
    real(fpk)   , dimension(46) :: USA_LOGPRESSURE (46)
    real(fpk)   , dimension(46) :: USA_TEMPERATURE (46)

!  OUTPUTS
!  =======

    !      Trace gas Layer column amounts and cross-sections
    !      Units should be consistent, e.g:
    !      X-sections  are in cm^2/mol  or [DU]
    !      gas columns are in mol/cm^2  or inverse [DU]

    integer                                           :: NLAYERS
    REAL(fpk)   , DIMENSION( maxlay )                 :: LAYER_AIRCOLUMNS
    REAL(fpk)   , DIMENSION( maxlay )                 :: LAYER_O2O2COLUMNS
    REAL(fpk)   , DIMENSION( maxlay )                 :: LAYER_TEMPERATURES
    REAL(fpk)   , DIMENSION( 0:maxlay )               :: LEVEL_HEIGHTS
    REAL(fpk)   , DIMENSION( maxlev   )               :: LEVEL_LOGPRESS

!  Local
!  -----

    integer            :: nlevels, j, i
    real(fpk)          :: psv, pmid, rho_1, rho_2, rho_a, diff, temp, xo2_kb
    real(fpk)          :: level_pressures(31),level_temperatures(31), constant

!  Loschmidt's number (particles/cm3), STP values

    real(fpk)   , PARAMETER :: RHO_STANDARD = 2.68675D+19
    real(fpk)   , PARAMETER :: PZERO = 1013.25D0
    real(fpk)   , PARAMETER :: TZERO = 273.15D0
    real(fpk)   , PARAMETER :: xo2 = 0.209476_fpk, kb = 1.381e-23_fpk

!  Get the pressures from GEOSCHEM
!    nlayers = 30 (hardwired, equal to the GEOSCHEM value)
!    pmid = pressure corresponding to ETA value (layer midpoint)

      psv = surface_pressure - vacuum_pressure
      nlayers = 30
      nlevels = nlayers + 1
      level_pressures(31) = surface_pressure
      do j = 1, 30
        pmid = psv * etagc(j) + vacuum_pressure
        level_pressures(31-j) = 2.0d0*pmid - level_pressures(32-j)
      enddo
      do j = 1, 31
        level_logpress(j) = dlog(level_pressures(j))
      enddo

!  Interpolate USA height and temperature to the GEOSCHEM grid

      call lintp2 &
          (usa_nlevels,usa_logpressure,usa_temperature,&
           nlevels, level_logpress, level_temperatures )

      call lintp2 &
          (usa_nlevels,usa_logpressure,usa_heights,&
           nlevels, level_logpress, level_heights )

!  gas constant

      CONSTANT = 1.0D+05 * RHO_STANDARD * TZERO / PZERO

!  Develop layer temperatures, layer air densities, assign height grid

      DO I = 1, NLAYERS
        RHO_1 = LEVEL_PRESSURES(I)   / LEVEL_TEMPERATURES(I)
        RHO_2 = LEVEL_PRESSURES(I+1) / LEVEL_TEMPERATURES(I+1)
        TEMP  = 0.5D0*(LEVEL_TEMPERATURES(I)+LEVEL_TEMPERATURES(I+1))
        RHO_A = 0.5D0 * CONSTANT * ( RHO_1 + RHO_2 )
        DIFF  = LEVEL_HEIGHTS(I) - LEVEL_HEIGHTS(I+1)
        LAYER_TEMPERATURES(I)   = TEMP
        LAYER_AIRCOLUMNS(I)     = DIFF * RHO_A
      ENDDO

!  Develop O2O2 densities

      xo2_kb  = xo2 * 1.0D-4 / kb
      DO I= 1, NLAYERS
        RHO_1 = xo2_kb * LEVEL_PRESSURES(I)   / LEVEL_TEMPERATURES(I)
        RHO_2 = xo2_kb * LEVEL_PRESSURES(I+1) / LEVEL_TEMPERATURES(I+1)
        RHO_A = 0.5D0 * ( RHO_1 * RHO_1 + RHO_2 * RHO_2 )
        DIFF  = LEVEL_HEIGHTS(I) - LEVEL_HEIGHTS(I+1)
        LAYER_O2O2COLUMNS(I)     = DIFF * RHO_A * 1.0E+5_fpk
      ENDDO

!  Finish

      RETURN
end subroutine develop_geoschem_pth

subroutine develop_geoschem_so2                            & ! Input
         ( MAXLAY, do_column_linearization,                & ! Input
           NLAYERS, actual_so2column,                      & ! Input
           layer_aircolumns, so2gc_shapefactors,           & ! Input
           so2_profile, so2_profile_deriv )                  ! Output

      implicit none

!  Inputs
!  ------

!  dimensioning

   integer :: maxlay

!  Column linearization flag

   logical :: do_column_linearization

!  Number of layers and air columns

    integer                                           :: NLAYERS
    REAL(fpk)   , DIMENSION( maxlay )                 :: LAYER_AIRCOLUMNS

!  Actual SO2 column (in [DU])

   real(fpk)    ::  actual_so2column

!  GEOSCHEM SO2 shap factors

   real(fpk)    ::  so2gc_shapefactors(30)

!  Outputs
!  -------

!  profiles

    REAL(fpk)   , DIMENSION( maxlay )                 :: so2_profile
    REAL(fpk)   , DIMENSION( maxlay )                 :: so2_profile_deriv

!  Local
!  -----

   integer      :: i, i1
   real(fpk)    :: total_air, total_so2, rr, sum
   real(fpk)   , parameter :: du_to_cm2 = 2.66868D+16

!  Total Air

      total_air = 0.0d0
      do i = 1, nlayers
        total_air = total_air + layer_aircolumns(i)
      enddo

!  Assign So2 columns in mol/cm^2.
!    AVERAGE VMR for the layer, uses the shape factor "as is"

      TOTAL_SO2 = actual_so2column * DU_TO_CM2
      RR = TOTAL_SO2 / TOTAL_AIR
      sum = 0.0d0
      DO I = 1, NLAYERS
        I1 = NLAYERS + 1 - I
        so2_profile(I) = RR*LAYER_AIRCOLUMNS(I) * so2gc_shapefactors(i1)
        sum = sum +  so2_profile(I)
      ENDDO

!  rescale profile so that total SO2 = value at the start !!!!
!   Normalized Form of the Column Derivatives

      RR = du_to_cm2 * actual_so2column / sum
      sum = 0.0d0
      DO I = 1, NLAYERS
        so2_profile(I) = so2_profile(I) * RR
        sum = sum + so2_profile(I)
      ENDDO
      
!  Derivative
!    In the column case, this is an unnormalized derivative

      if ( do_column_linearization ) THEN
        do I = 1, Nlayers
          so2_profile_deriv(I) = so2_profile(I)/TOTAL_SO2
        enddo
      else
        do I = 1, Nlayers
          so2_profile_deriv(I) = 1.0d0
        enddo
      endif

!  Finish

      return
end subroutine develop_geoschem_so2

                  
subroutine extract_usadata_pth ( do_standard, latitude, month, &
            USA_NLEVELS,USA_HEIGHTS,USA_LOGPRESSURE, USA_TEMPERATURE)

      implicit none

!  Inputs
!  ------

!  Variable for using Standard Atmosphere

   logical      :: do_standard

!  Time and position

   real(fpk)    :: latitude
   integer      :: month

!  Output
!  ------

!  PTH and number of levels ( = 46 )

    integer                     :: USA_NLEVELS
    real(fpk)   , dimension(46) :: USA_HEIGHTS     (46)
    real(fpk)   , dimension(46) :: USA_LOGPRESSURE (46)
    real(fpk)   , dimension(46) :: USA_TEMPERATURE (46)

!  Local
!  -----

    integer                     :: usadata_nlevels, N, K, KS, band, season
    real(fpk)   , dimension(46) :: USADATA_PRESSURES    (46,6)
    real(fpk)   , dimension(46) :: USADATA_TEMPERATURES (46,6)

!  Get Reference pressures and temperatures from USA

      OPEN ( 1, file = 'physics_data/pth_usa.dat', STATUS = 'OLD' )
      READ (1,*)
      READ (1,*) USADATA_NLEVELS
      DO N = 1, USADATA_NLEVELS
         READ(1,*)USA_HEIGHTS(N),(USADATA_PRESSURES(N,K),K=1,6)
      ENDDO
      DO N = 1, USADATA_NLEVELS
         READ(1,*)USA_HEIGHTS(N),(USADATA_TEMPERATURES(N,K),K=1,6)
      ENDDO
      CLOSE (1)

!  Select atmosphere

      if ( do_standard ) then
        KS = 6
      else
        if ( dabs(latitude) .lt. 30.0d0 ) then
          ks = 1
        else
          band = 1
          if ( dabs(latitude) .ge. 60.0d0 ) band = 2
          season = 1
          if ( month .gt. 3 .and. month .lt. 10 ) season = 2
          ks = 1 + 2*(band-1) + season
        endif
      endif

!  set output

      USA_NLEVELS = USADATA_NLEVELS
      DO N = 1, USADATA_NLEVELS
        USA_LOGPRESSURE(N) = DLOG(USADATA_PRESSURES(N,KS))
        USA_TEMPERATURE(N) = USADATA_TEMPERATURES(N,KS)
      ENDDO

!  finish

      RETURN
end subroutine extract_usadata_pth

! @@@@@@@@@@@@ NEW ROUTINES ADDED, 18 FEBRUARY 2013

!JPARK##########################################################################################
      SUBROUTINE TOMS_STDPRF_T  &
           ( MXLAY, IFILE, SURFACE_PRESS, NLAYS,                & ! Input
             LEVEL_PRESS, LEVEL_LOGPRESS, LEVEL_HEIGHTS,        & ! Output
             LAYER_TEMPERATURES, LAYER_AIRCOLUMNS, FAIL, MESSAGE) ! Output


      IMPLICIT NONE

      INTEGER,                       INTENT (IN) :: MXLAY
      REAL(fpk),                     INTENT (IN) :: SURFACE_PRESS
      CHARACTER*(*),                 intent(in)  :: IFILE 


      INTEGER,                       INTENT(OUT) :: NLAYS
      REAL(fpk),   DIMENSION(MXLAY), INTENT(INOUT) :: LAYER_TEMPERATURES, &
                                                    LAYER_AIRCOLUMNS
      REAL(fpk), DIMENSION(0:MXLAY), INTENT(INOUT) :: LEVEL_PRESS, &
                                                    LEVEL_HEIGHTS
      REAL(fpk), DIMENSION(0:MXLAY), INTENT(INOUT) :: LEVEL_LOGPRESS
      LOGICAL,                       INTENT(OUT) :: FAIL
      CHARACTER*(*),                 INTENT(OUT) :: MESSAGE

      !-- local variables
      REAL(fpk),                DIMENSION(MXLAY) :: TEM0
      REAL(fpk),                       PARAMETER :: RGAS = 8314.32D0 
      REAL(fpk),                       PARAMETER :: MAIR = 28.9664D0
      REAL(fpk),                       PARAMETER :: GRAV = 9.80665D0
      REAL(fpk),                       PARAMETER :: ERAD = 6356.76D0  
      REAL(fpk) :: MG_OVER_R, ACONST, dZHALF, dLOGP, ZZ, dZ, TOA_HEIGHT, DENS1, DENS2, HDIFF, RHO_STP
      INTEGER           :: LUN, IOS, I


!  Initialize

      FAIL = .FALSE. ; MESSAGE = ' '
      NLAYS = 11

      !-- get layer temperature
      LUN=9
      OPEN(UNIT=LUN, FILE=TRIM(IFILE), IOSTAT=IOS, STATUS='OLD', ACTION='READ')
      IF (IOS/=0) THEN 
        FAIL = .true.
        MESSAGE = 'I/O Error for TOMS Standard temperature profile file, where is it?'
      ELSE
        READ(LUN,"(4X,11F7.2)") TEM0(1:NLAYS)  !++ BOA --> TOA
      ENDIF
      CLOSE (LUN)
      if ( fail) return

      DO i=1,NLAYS  !++ TOA --> BOA
        LAYER_TEMPERATURES(i) = TEM0(NLAYS+1-i)
      ENDDO

      !-- get pressure
      DO i=NLAYS,0,-1  !++ TOA --> BOA
        LEVEL_LOGPRESS(NLAYS-i) = LOG(SURFACE_PRESS) -REAL(i,fpk) *LOG(two)
      ENDDO
      LEVEL_PRESS(0:NLAYS) = EXP(LEVEL_LOGPRESS(0:NLAYS)) 

      !-- calc heights using hydrostatic eq.
      MG_OVER_R = (MAIR *GRAV) /RGAS
      ACONST = (one /MG_OVER_R) /1000.0_fpk

      LEVEL_HEIGHTS(NLAYS) = zero  !++ @TOA
      dZHALF = 0.482_fpk
      DO i=NLAYS,1,-1
        dLOGP = LEVEL_LOGPRESS(i-1) -LEVEL_LOGPRESS(i)
        ZZ = LEVEL_HEIGHTS(i) +dZHALF
        dZ = -ACONST *((ERAD +ZZ) /ERAD)**2 *LAYER_TEMPERATURES(i) *dLOGP
        LEVEL_HEIGHTS(i-1) = LEVEL_HEIGHTS(i) +dZ
        dZHalf = dZ /two
      ENDDO

      !-- find air density, mole/cm^2
      TOA_HEIGHT = two *LEVEL_HEIGHTS(0) -LEVEL_HEIGHTS(1)
!!!!  LEVEL_HEIGHTS(NLAYS) = TOA_HEIGHT

      RHO_STP  = 2.68675E+24_fpk *273.15_fpk /1013.25_fpk
      DENS1 = 0.0d0
      DENS2 = LEVEL_PRESS(1) *RHO_STP /LAYER_TEMPERATURES(1)
      HDIFF = TOA_HEIGHT -LEVEL_HEIGHTS(1)
      LAYER_AIRCOLUMNS(1) = HDIFF * half *(DENS1 +DENS2)
      DENS1 = DENS2
      DO i=2,NLAYS
        DENS2 = LEVEL_PRESS(i) *RHO_STP /LAYER_TEMPERATURES(i)
        HDIFF = LEVEL_HEIGHTS(i-1) - LEVEL_HEIGHTS(i)
        LAYER_AIRCOLUMNS(i) = HDIFF * half *(DENS1 +DENS2)
        DENS1 = DENS2
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE TOMS_STDPRF_T

!JPARK#########################################################################

      SUBROUTINE TOMS_STDPRF_O3                                 &
              (INPUT_COLUMN_O3, LATITUDE, IFILE, MXLAY, NLAYS,  & ! Input
               TOMS_STD_O3, TOMS_STD_O3_DERIV, FAIL, MESSAGE)     ! Output

      IMPLICIT NONE
 
      !-- input variables
      INTEGER,                       INTENT (IN) :: MXLAY, &
                                                    NLAYS
      REAL(fpk),                     INTENT (IN) :: INPUT_COLUMN_O3, &
                                                    LATITUDE

      CHARACTER*(*),                 intent(in)  :: IFILE 

      !-- output variables
      REAL(fpk),   DIMENSION(MXLAY), INTENT(INOUT) :: TOMS_STD_O3, &
                                                      TOMS_STD_O3_DERIV
      LOGICAL,                       INTENT(OUT) :: FAIL
      CHARACTER*(*),                 INTENT(OUT) :: MESSAGE
      !-- local variables
      INTEGER,                         PARAMETER :: NPRFS=21
      REAL(fpk),                       PARAMETER :: DU2CM2=2.66868D+16
      REAL(fpk),          DIMENSION(NPRFS,NLAYS) :: STD_O3_PRF
      REAL(fpk),     DIMENSION(:,:), ALLOCATABLE :: ROI_O3_PRF
      REAL(fpk),     DIMENSION(10)               :: COLUMN_O3
      INTEGER,                      DIMENSION(4) :: NGP
      INTEGER                                    :: LUN, IOS, I, J, GP, CID
      DATA  NGP   /1, 3, 8, 11/, &
           COLUMN_O3 /125, 175, 225, 275, 325, 375, 425, 475, 525, 575/

!  Few extras

      LOGICAL   :: LOOP
      INTEGER   :: N, NGRID, K1, K2
      real(fpk) :: X, GRID(10), DIFF, XFAC, CGRAD, PROF1, PROF2, F1, F2

!  Initialize

      FAIL = .FALSE. ; MESSAGE = ' '

!-- check geophysical location
      GP  = 3  !++ default mid-latitude
      CID = 3  !++        "
      IF (ABS(LATITUDE).LT.30.) THEN       !++ low latitude
        GP = 2
      ELSE IF (ABS(LATITUDE).GE.60.) THEN  !++ high latitude
        GP = 4
        CID = 1
      ENDIF

      !-- read TOMS Standard Ozone Profiles 
      LUN = 9
      OPEN(UNIT=LUN, FILE=TRIM(IFILE), IOSTAT=IOS, STATUS='OLD', ACTION='READ')
      IF (IOS/=0) THEN
        FAIL = .true.
        MESSAGE = 'I/O Error for TOMS Standard ozone profile file, where is it?'
      ELSE
        READ(LUN,"(4X,11F7.2)") (STD_O3_PRF(I,1:NLAYS), I=1,NPRFS)
      ENDIF
      CLOSE (LUN)
      if ( FAIL ) RETURN

!  Allocate

      ALLOCATE( ROI_O3_PRF(NGP(GP),NLAYS) )
      DO i=1,NGP(GP)
        DO j=1,NLAYS
          ROI_O3_PRF(i,j) = STD_O3_PRF(NGP(GP-1)+i,j)
        ENDDO
      ENDDO

      !-- just find same profles, not interpolate
 !     DO i=1,NGP(GP)
 !       IF (INPUT_COLUMN_O3 .EQ. COLUMN_O3(CID-1+i)) THEN
 !         DO j=1,NLAYS
 !           TOMS_STD_O3(j) = ROI_O3_PRF(i,NLAYS+1-j) *DU2CM2
 !         ENDDO
 !         FAIL = .FALSE.
 !         EXIT
 !       ENDIF
 !     ENDDO

!  Save all profile for this latitude zone

      NGRID = NGP(GP)
      DO i=1,NGP(GP)
         GRID(I) = COLUMN_O3(CID-1+i)
      ENDDO
      X = INPUT_COLUMN_O3

!  Set the ozone, by Interpolating TOMSV8 Ozone to our incoming column

      IF ( X .LT. GRID(2) ) THEN
        K1 = 1
      ELSE IF ( X .GE. GRID(NGRID-1) ) THEN
        K1 = NGRID - 1
      ELSE
        LOOP = .TRUE.
        K1 = 0
        DO WHILE (LOOP)
          K1 = K1 + 1
          IF ( X .GE. GRID(K1) .AND. X .LT. GRID(K1+1) ) LOOP=.FALSE.
        ENDDO
      ENDIF

      K2 = K1 + 1
      F1 = ( X - GRID(K1) ) / ( GRID(K2) - GRID(K1) )
      F2 = 1.0D0 - F1

!  Set profile and derivative, using interpolation

      DIFF = GRID(K2) - GRID(K1)
      XFAC = INPUT_COLUMN_O3 - GRID(K1)
      DO N = 1, NLAYS
         Prof1 = ROI_O3_PRF(K1,NLAYS+1-N)
         Prof2 = ROI_O3_PRF(K2,NLAYS+1-N)
         CGRAD = ( Prof2 - Prof1 ) / Diff
         TOMS_STD_O3(n)    = Prof1 + XFAC * cgrad
         TOMS_STD_O3_DERIV(n) = cgrad
         if ( dabs(cgrad).lt.1.0d-12)TOMS_STD_O3_DERIV(n) =0.0d0
         f1 = f1 + prof1 ; f2 = f2 + prof2
      ENDDO

!  Convert Profiles to Mol/cm/cm, Make normalized derivative

      DO N = 1, nlays
         TOMS_STD_O3(n) = TOMS_STD_O3(n) * du2cm2
         TOMS_STD_O3_DERIV(n) = TOMS_STD_O3_DERIV(n) * INPUT_COLUMN_O3 * du2cm2
      ENDDO

!  USEFUL Check

!      write(*,*)sum(toms_std_o3(1:nlays))/du2cm2, input_column_o3

!  Former code
!      TOMS_STD_O3_DERIV = TOMS_STD_O3 /INPUT_COLUMN_O3 /1000.D0  !++ not important to cal. radiance

      IF (ALLOCATED( ROI_O3_PRF )) DEALLOCATE( ROI_O3_PRF )

!  Normal return

      RETURN
      END SUBROUTINE TOMS_STDPRF_O3
!####################################################################################################


subroutine pthg_finer                                           &
  ( maxlay, maxgas, do_gas_jacobians, NGASES,                   & ! INPUT
    Lower_boundary, Upper_Boundary, modulus,                    & ! INPUT
    NLAYERS, LEVEL_HEIGHTS, LEVEL_PRESS, LAYER_TEMPERATURES,    & ! OUTPUT
    LAYER_AIRCOLUMNS, GAS_PROFILES, GAS_PROFILES_DERIV,         & ! OUTPUT
    FAIL, MESSAGE )                                               ! OUTPUT

  implicit none

!  Purpose.

!   Gien a PTHG atmosphere, we want to make finer layer subdivisions
!   for that part of the atmosphere lying between the Boundary limits.

!  INPUTS
!  ======

!  dimensioning

   integer, intent(in)       :: maxlay, maxgas

!  Gas control

   logical, intent(in)       :: do_gas_jacobians
   integer, intent(in)       :: NGASES

!  Finer subdivision control

   real(fpk)   , intent(in)  :: Lower_boundary, Upper_Boundary
   integer, intent(in)       :: modulus

!  OUTPUTS (modified PTHG arrays)
!  =======

    !      Trace gas Layer column amounts and cross-sections
    !      Units should be consistent, e.g:
    !      X-sections  are in cm^2/mol  or [DU]
    !      gas columns are in mol/cm^2  or inverse [DU]

    integer                                  , intent(inout)  :: NLAYERS
    REAL(fpk)   , DIMENSION( maxlay )        , intent(inout)  :: LAYER_AIRCOLUMNS
    REAL(fpk)   , DIMENSION( maxlay )        , intent(inout)  :: LAYER_TEMPERATURES
    REAL(fpk)   , DIMENSION( 0:maxlay )      , intent(inout)  :: LEVEL_HEIGHTS
    REAL(fpk)   , DIMENSION( 0:maxlay )      , intent(inout)  :: LEVEL_PRESS
    REAL(fpk)   , DIMENSION( maxlay, maxgas ), intent(inout)  :: GAS_PROFILES
    REAL(fpk)   , DIMENSION( maxlay, maxgas ), intent(inout)  :: GAS_PROFILES_DERIV

!  Exception handling

    logical      , intent(out) :: fail 
    character*60 , intent(out) :: message

!  Local parameter, must be at least 4 + 2*maxgases

   integer, parameter :: maxhelp = 8

!  Local

   integer :: n, ntop, nbot, nc, n41, nlayers_new, m, g
   logical :: trawl
   real*8  :: h1, h2, p1, p2, q1, q2, pfac, hdiff, pdiff, qdiff, qgrad
   real*8  :: help(maxlay,maxhelp), diff
   character*3 :: c3

!  intiaizlie

   fail    = .false.
   message = ' '

   if (  4 + 2*ngases .gt. maxhelp ) then
      fail = .true.
      message = 'Local GAS dimension maxhelp, not big enough'
      return
   endif

!  Holding arrays

   n41 = 4 + ngases
   do n = 1, nlayers
      Help(n,1) = LAYER_AIRCOLUMNS(n)
      Help(n,2) = LAYER_TEMPERATURES(n)
      Help(n,3) = LEVEL_HEIGHTS(n)
      Help(n,4) = LEVEL_PRESS(n)
      do g = 1, ngases
        Help(n,4+g) = GAS_PROFILES(n,g)
      enddo
   enddo
   if ( do_gas_jacobians ) then
      do n = 1, nlayers
         do g = 1, ngases
           Help(n,n41+g) = GAS_PROFILES_DERIV(n,g)
         enddo
      enddo
   endif

!  Find the grid levels adjacent to the limit boundaries

   trawl = .true. ; n = 0
   do while (trawl .and.n.lt.nlayers)
     n = n + 1
     if ( level_heights(n) .le. Upper_boundary )  trawl = .false.
   enddo
   ntop = n - 1; if ( level_heights(n) .eq. Upper_boundary ) ntop = n

   trawl = .true. ; n = ntop - 1
   do while (trawl .and.n.lt.nlayers)
     n = n + 1
     if ( level_heights(n) .le. Lower_boundary )  trawl = .false.
   enddo
   nbot = n

!  Exception handling

   if ( trawl ) then
      fail    = .true.
      message = 'Boundary limits not found within present height grid'
      return
   endif

!  Check new total number of layers is not exceeding dimensions

   nlayers_new = nlayers + (modulus-1) * ( nbot-ntop )
   if ( nlayers_new .gt. maxlay ) then
      write(c3,'(I3)')nlayers_new+2
      fail    = .true.
      message = 'Too many layers; increase dimension maxlay to at least '//c3
      return
   endif

!  From ntop to nbot, use finer subdivisions

   nc = ntop
   h1 = help(ntop,3) ; p1 = help(ntop,4) ; q1 = dlog(p1)
   do n = ntop + 1, nbot
      hdiff = help(n-1,3) - help(n,3)
      qdiff = dlog(help(n-1,4)/help(n,4))
      pdiff = help(n,4) - help(n-1,4)
      qgrad = - qdiff / hdiff
      diff  = (help(n-1,3) - help(n,3))/modulus 
      do m = 1, modulus
         nc = nc + 1
         h2 = h1 - diff
         q2 = q1 + qgrad * diff
         p2 = dexp(q2)
         pfac = ( p2 - p1 ) / pdiff
         level_heights(nc) = h2
         level_press(nc)   = p2
         layer_temperatures(nc) = help(n,2)
         layer_aircolumns(nc)   = pfac * help(n,1)
         do g = 1, ngases
            GAS_PROFILES(nc,g) = pfac * Help(n,4+g)
         enddo
         if ( do_gas_jacobians ) then
            do g = 1, ngases
               GAS_PROFILES_DERIV(nc,g) = pfac * Help(n,n41+g)
            enddo
         endif
         p1 = p2
         h1 = h2
         q1 = q2
      enddo
   enddo

!  Copy down to bottom

   do n = nbot + 1, nlayers
      nc = nc + 1
      level_heights(nc) = help(n,3)
      level_press(nc)   = help(n,4)
      layer_temperatures(nc) = help(n,2)
      layer_aircolumns(nc)   = help(n,1)
      do g = 1, ngases
         GAS_PROFILES(nc,g) = Help(n,4+g)
      enddo
      if ( do_gas_jacobians ) then
         do g = 1, ngases
            GAS_PROFILES_DERIV(nc,g) = Help(n,n41+g)
         enddo
      endif
   enddo

!  Final nlayers

    if (nbot.eq.nlayers)level_heights(nc) = help(nbot,3)
    nlayers = nc

!  Finish

   return
end subroutine pthg_finer

!  End module

end module pthg_profiles_m

