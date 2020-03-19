subroutine geocape_profile_reader_1 &
           ( filename,  &     ! input
             profile_data, footprint_data, nlayers, &
             zsur, time_data, utcsec, &
             message, fail ) ! output

!  The reader to read in profiles in a .asc file

!  input: filename
!  output: profile_data(75,8), footprint_data(3)

!  profile_data(75,8) has following 8-column attributes
!  at 75 pressure levels. The fill value is -999.

!    1.   LEVELS Pressure (hPa),
!    2.   LEVELS Atmospheric temperature (K),
!    3.   H2O (vmr),
!    4.   O3 (vmr),
!    5.   NO2 (vmr),
!    6.   HCHO (vmr),
!    7.   SO2 (vmr),
!    8.   CO (vmr)

!  footprint_data(3) has following data
!    1.   Latitude,
!    2.   Longitude,
!    3.   Surface temperature (K)

   implicit none

!  input filename

   character(len=*),                intent(in)  :: filename

!  Main output

   integer,                        intent(out) :: nlayers
   real(kind=8), dimension(75,8),  intent(out) :: profile_data
   real(kind=8), dimension(3),     intent(out) :: footprint_data
   real(kind=8),                   intent(out) :: zsur
   integer(kind=4), dimension(3),  intent(out) :: time_data
   real(kind=8),                   intent(out) :: utcsec

!  Exception handling

   logical,       intent(INOUT) :: fail
   character*(*), intent(INOUT) :: message

!  Local variables

   integer      :: i, j, nr, nc, np
   real(kind=8), dimension(75,9) :: tmp_data
   real(kind=8) :: surfTemp, lat, lon
   character*80 :: dummy, time_string

!  initialize

   fail    = .false.
   message = ' '
   do i = 1, 75
      do j = 1,8
         profile_data(i,j) = -999.
      enddo
   enddo

   np = 1
   do while (filename(np:np).ne. ' ')
     np = np + 1 
   enddo
   np = np - 1

!  Open and read file

   open(1,file=filename(1:np),err=90, status='old')
   read(1,*)
   read(1,*) dummy,dummy,nr,dummy,nc
   read(1,*) dummy,dummy,surfTemp
   read(1,*) dummy,dummy,zsur
   read(1,*) dummy,dummy,lat
   read(1,*) dummy,dummy,lon
   read(1,*) dummy,dummy,time_string
   read(1,*) dummy,dummy, utcsec
   do i = 1,3
      read(1,*)
   enddo

   footprint_data(1) = lat
   footprint_data(2) = lon
   footprint_data(3) = surfTemp

!  Extract year/month/date from the "time_string" variable

   read(time_string(1:4),'(I4)') time_data(1)
   read(time_string(5:6),'(I2)') time_data(2)
   read(time_string(7:8),'(I2)') time_data(3)

   do i = 1, nr
     read(1,*) (tmp_data(i,j),j=1,nc)
     do j = 1,2
        profile_data(i,j) = tmp_data(i,j+1)
     enddo
     do j = 3,nc-1
! turn species VMR to ppm
        profile_data(i,j) = tmp_data(i,j+1)*1.e6
     enddo
   enddo

   nlayers = nr-1
   
   close(1)

   return

! error return

90 continue
   fail = .true.
   message = 'Open failure for profile file = '//filename(1:LEN(filename))
   return

end subroutine geocape_profile_reader_1

subroutine geocape_profile_setter_1                  &
    ( MaxLayers, MaxGases, nlayers, ngases,          &  ! Input
      which_gases, zsur, profile_data,               &  ! Input
      heights, temperatures, pressures,              &  ! Output
      aircolumns, daircolumns_dT,                    &  ! Output
      gas_partialcolumns, gas_totalcolumns,          &  ! Output
      fail, message )                                  ! Output

   implicit none

!  inputs
!  ------

!  Dimensioning

!  Dimensioning

   integer, intent(in) :: MaxLayers, MaxGases
   integer, intent(in) :: nlayers, ngases

!  Data from file-read

   real(kind=8), dimension(75,8),  intent(in) :: profile_data

!  Trace gas control

   character(Len=4), dimension ( MaxGases ), intent(in) :: which_gases

!  Surface elevation

   real(kind=8),                           intent(in) :: zsur

!  Output
!  ------

!  Atmospheric quantities (PTH)

   real(kind=8), dimension ( 0:MaxLayers ), intent(out)  :: heights
   real(kind=8), dimension (   MaxLayers ), intent(out)  :: temperatures
   real(kind=8), dimension ( 0:MaxLayers ), intent(out)  :: pressures

!  Air density Partial columns and T-derivative

   real(kind=8), dimension (MaxLayers), intent(out)          :: aircolumns
   real(kind=8), dimension (MaxLayers), intent(out)          :: daircolumns_dT

!  Trace gas partial columns (profile) and total columns

   real(kind=8), dimension (MaxLayers,MaxGases), intent(out) :: gas_partialcolumns
   real(kind=8), dimension (MaxGases), intent(out)           :: gas_totalcolumns

!  Exception handling

   logical, intent(INOUT)       :: fail
   character*(*), intent(INOUT) :: message

!  Local variables
!  ---------------

!  Array of level temperatures

   real(kind=8), dimension(0:MaxLayers) :: leveltemp

!  Array of derived gas constants

   real(kind=8), dimension(MaxLayers)   :: gasconstants

!  help variables

   integer       :: n, n1, g, ngas_check
   real(kind=8)  :: ccon, rho1, rho2, col, pp, delp, temp, airc, avit

!  Parameters: Loschmidt's number (particles/cm3), STP parameters

   real(kind=8), parameter ::  RHO_STAND = 2.68675D+19
   real(kind=8), parameter ::  PZERO     = 1013.25D0
   real(kind=8), parameter ::  TZERO     = 273.15D0
   real(kind=8), parameter ::  RHO_ZERO  = RHO_STAND * TZERO / PZERO
   real(kind=8), parameter ::  CONST     = 1.0D+05 * RHO_ZERO
   real(kind=8), parameter ::  DU_TO_CM2 = 2.68668D16
   real(kind=8), parameter ::  O2RATIO   = 0.2095D0

!  initialize

   fail    = .false.
   message = ' '

!  Set level data for P, T and H
!    --- Convert to [km] units (height), hPa (pressures)

   do n = 0, nlayers
      n1 = nlayers - n + 1
      leveltemp(n) = profile_data(n1,2)
      pressures(n) = profile_data(n1,1)
   enddo

!  heights by Hydrostatic Eqn. (includes TOA)
     
   heights(nlayers) = zsur/1000.d0
   ccon = - 9.81d0 * 28.9d0 / 8314.0d0 * 500.0d0
   do n = nlayers, 1, -1
      avit = (1.0d0/leveltemp(n-1))+(1.0d0/leveltemp(n))
      heights(n-1) = heights(n) - dlog(pressures(n)/pressures(n-1))/avit/ccon
   enddo

!  develop air density
!  -------------------

!    Fiddle "pressure-difference method)" 
!         for derivative of Air density w.r.t Temperature

   do n = 1, nlayers
      n1 = n - 1
      rho1 = pressures(n1)/ leveltemp(n1)
      rho2 = pressures(n)/ leveltemp(n)
      temp = 0.5d0 * (leveltemp(n1)+leveltemp(n))
      airc = 0.5d0 * const * ( rho1 + rho2 ) * (heights(n1)-heights(n))
      delp = pressures(n) - pressures(n1)
      gasconstants(n)   = airc * temp / delp
      temperatures(n)   = temp
      aircolumns(n)     = gasconstants(n) * delp / temperatures(n)
      daircolumns_dT(n) = - aircolumns(n) / temperatures(n)
   enddo

!  Develop gas partial columns
!  ---------------------------

!    First default, 3 June 2009. 4 UV gases (O3, NO2, HCHO, SO2)
!    T-derivatives not required explicitly, handled by air column T-deriv above.

   pp = 1.0d-06
   ngas_check = 0
   do g = 1, ngases
      if ( which_gases(g) .eq. 'O3  ' ) then
         ngas_check = ngas_check + 1
         do n = 1, nlayers
            n1 = nlayers + 1 - n
            gas_partialcolumns(n,g) = pp * aircolumns(n) * profile_data(n1,4)  ! O3
         enddo
      else if ( which_gases(g) .eq. 'NO2 ' ) then
         ngas_check = ngas_check + 1
         do n = 1, nlayers
            n1 = nlayers + 1 - n
            gas_partialcolumns(n,g) = pp * aircolumns(n) * profile_data(n1,5) ! NO2
         enddo
      else if ( which_gases(g) .eq. 'HCHO' ) then
         ngas_check = ngas_check + 1
         do n = 1, nlayers
            n1 = nlayers + 1 - n
            gas_partialcolumns(n,g) = pp * aircolumns(n) * profile_data(n1,6) ! HCHO
         enddo
      else if ( which_gases(g) .eq. 'SO2 ' ) then
         ngas_check = ngas_check + 1
         do n = 1, nlayers
            n1 = nlayers + 1 - n
            gas_partialcolumns(n,g) = pp * aircolumns(n) * profile_data(n1,7) ! SO2
         enddo
      else if ( which_gases(g) .eq. 'H2O ' ) then
         ngas_check = ngas_check + 1
         do n = 1, nlayers
            n1 = nlayers + 1 - n
            gas_partialcolumns(n,g) = pp * aircolumns(n) * profile_data(n1,3)  ! H2O
         enddo
      else if ( which_gases(g) .eq. 'O4  ' ) then
         ngas_check = ngas_check + 1
         do n = 1, nlayers
            n1 = nlayers + 1 - n
            gas_partialcolumns(n,g) = (aircolumns(n) * O2RATIO) ** 2.0 / (heights(n-1)-heights(n)) / 1.0D5 ! O4
         enddo
      endif
   enddo

!  Check that All input gases have been found

   if ( ngas_check .ne. ngases ) then
      message = 'Not all desired trace gases are present in data set: Reduce choice!'
      fail = .true.
      return
   endif

!  Set non-physical entries to zero.

   do g = 1, ngases
      do n = 1, nlayers
         if (gas_partialcolumns(n,g).lt.0.0d0)gas_partialcolumns(n,g)=0.0d0
      enddo
   enddo

!  Develop total columns in [DU]

   do g = 1, ngases
      col = 0.0d0
      do n = 1, nlayers
         col = col + gas_partialcolumns(n,g)
      enddo
      gas_totalcolumns(g) = col / du_to_cm2
   enddo

!  Finish

   return
end subroutine geocape_profile_setter_1

