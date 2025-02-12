module GEMSTOOL_NSW_cross_sections_m

!  Original construction for GEMS, 23 July 2013
!  UVN-specific Version for GEMS, 12 August 2013.
!  NSW-specific Version for GEMS, 12 August 2013. (This module)

!  Older notes
!  -----------

!  Revision. 03 September 2012
!   - GET_O2O2XSEC_1, introduced for O2-O2 absorption (Greenblatt 1990 data)

!  Revision. 18 February 2013
!   - Added several more Cross-sections, from GEOCAPE (SAO) compilation.
!   - separated Cross-section data files in Physics_data

!  Revision. 21 October 2013
!     Rayleigh function separated

!  Use modules
!  -----------

!  Use numerical Subroutines

   use GEMSTOOL_numerical_m   , only : spline, splint

!  Hitran modules

   use get_hitran_crs_m

!  Rayleigh function

   use Rayleigh_function_m

implicit none

!  Precision

integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Only main routine is public

public  :: GEMSTOOL_NSW_CROSS_SECTIONS
private :: reverse_1, fpk

contains

SUBROUTINE GEMSTOOL_NSW_CROSS_SECTIONS  &
          ( MAXWAVNUMS, MAXLAYERS, MAXGASES, CO2_PPMV_MIXRATIO,   & ! INPUT
            NLAYERS, NWAVNUMS, LEVEL_TEMPS, LEVEL_PRESS, WAVNUMS, & ! INPUT
            NGASES, WHICH_GASES, HITRAN_PATH,                     & ! INPUT   
            DO_GASES,                                             & ! IN/OUT    
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL, GAS_XSECS,             & ! OUTPUT
            FAIL, MESSAGE )                                         ! OUTPUT

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input arguments
!  ---------------

!  Dimensioning

      INTEGER     :: MAXWAVNUMS, MAXLAYERS, MAXGASES

!  wavelengths
 
      INTEGER     :: NWAVNUMS
      real(fpk),    dimension ( MAXWAVNUMS ) :: WAVNUMS 

!  CO2 mixing ratio

      real(fpk)    :: CO2_PPMV_MIXRATIO

!  Layer control

      INTEGER     :: NLAYERS

!  Level P and T

      real(fpk),    dimension ( 0:MAXLAYERS ) :: LEVEL_TEMPS 
      real(fpk),    dimension ( 0:MAXLAYERS ) :: LEVEL_PRESS 

!  Gas control

      integer                                  :: NGASES
      character(LEN=4), dimension ( maxgases ) :: WHICH_GASES

!  Hitran path name

      character*(*), intent(in) :: HITRAN_PATH

!  Input/Output arguments
!  ----------------------

      logical, dimension ( maxgases )          :: DO_GASES

!  Output arguments
!  ----------------

!  Rayleigh cross-sections and depolarization output

      real(fpk),    dimension ( MAXWAVNUMS ) :: RAYLEIGH_XSEC 
      real(fpk),    dimension ( MAXWAVNUMS ) :: RAYLEIGH_DEPOL

!  Gas cross-sections, all levels, all gases

      REAL(fpk),    DIMENSION ( MAXWAVNUMS, 0:maxlayers, maxgases ) :: GAS_XSECS  

!  status

      LOGICAL       ::        FAIL
      CHARACTER*(*) ::        MESSAGE

!  Local variables
!  ---------------

      integer            :: G, NLEVELS, W, W1

!  Hitran inputs (local); dummy output

      CHARACTER (LEN=6)  :: the_molecule
      LOGICAL            :: is_wavenum    ! *** T: wavenumber, F: nm ***
      REAL(KIND=FPK)     :: fwhm          ! *** in cm^-1 or nm ***
      real(KIND=FPK)     :: ATM_PRESS ( MAXLAYERS + 1 )
      INTEGER            :: errstat  ! Error status
      CHARACTER(Len=100) :: message_hitran  !
      REAL(KIND=FPK)     :: dummy (MAXWAVNUMS,1:maxlayers+1)

!  Local wavelength array

      REAL(KIND=FPK)     :: lambdas (MAXWAVNUMS)

!  Initialize exception

   fail = .false. ; message = ' '

!  set the pressures in bars

   nlevels = nlayers + 1
   atm_press(1:nlevels) = level_press(0:nlayers)  / 1013.25d0

!  Use wavenumbers !!!!

   fwhm = 0.0d0 ; is_wavenum = .true.

!  Gases

   DO G = 1, NGASES

      if ( which_gases(g) .eq. 'H2O ' ) the_molecule = 'H2O   '
      if ( which_gases(g) .eq. 'O2  ' ) the_molecule = 'O2    '
      if ( which_gases(g) .eq. 'CO2 ' ) the_molecule = 'CO2   '
      if ( which_gases(g) .eq. 'CH4 ' ) the_molecule = 'CH4   '

      call get_hitran_crs &
              ( Hitran_path, the_molecule, nwavnums, wavnums(1:nwavnums), &
                is_wavenum, nlevels, atm_press(1:nlevels), level_temps(0:nlayers), fwhm,  &
                GAS_XSECS(1:nwavnums,0:nlayers,g), errstat, message_hitran, dummy(1:nwavnums,1:nlevels) )
!        Fatal error --> Exit. Warning error --> write to screen and move on !!!!

      !   write(*,*)'Done Hitran calculation for gas = '//which_gases(g)

      if (errstat.eq.1 ) then
         fail = .true.
         message = adjustl(trim(message_hitran))
         return
      else if ( errstat.eq.2 ) then
         write(*,*)' * '//adjustl(trim(message_hitran))//', for gas = '//which_gases(g)
!mick fix - added logical statement
         do_gases(g) = .false.
      endif

   enddo

!  debug 24 October 2013, FD testing
!   do w = 1, nwavnums
!      write(87,*)w,GAS_XSECS(w,nlayers,1)
!   enddo

!  Get the Rayleigh. Input wavelengths in [nm], NOT ANGSTROMS

   do w = 1, nwavnums
      w1 = nwavnums + 1 - w ; lambdas(w) = 1.0d+07 / wavnums(w1)
   enddo
   CALL RAYLEIGH_FUNCTION                                &
     ( MAXWAVNUMS, CO2_PPMV_MIXRATIO, NWAVNUMS, lambdas, &
       RAYLEIGH_XSEC, RAYLEIGH_DEPOL )
   call reverse_1 ( RAYLEIGH_XSEC (1:nwavnums), nwavnums )
   call reverse_1 ( RAYLEIGH_DEPOL(1:nwavnums), nwavnums )

!  debug Rayleigh
!   do w = 1, nwavnums
!      write(*,*)w,wavnums(w), lambdas(nwavnums+1-w), RAYLEIGH_XSEC(w), RAYLEIGH_DEPOL(w)
!   enddo
!   pause'Ray'

!  FInish

  RETURN
END SUBROUTINE GEMSTOOL_NSW_CROSS_SECTIONS

!

SUBROUTINE reverse_1 ( inarr, num )
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  INTEGER, INTENT(IN) :: num
  INTEGER             :: i
  REAL (KIND=dp), DIMENSION(1: num), INTENT(INOUT) :: inarr
  REAL (KIND=dp), DIMENSION(1: num)                :: temp

  DO i = 1, num
     temp(i) = inarr(num - i + 1)
  ENDDO
  inarr = temp

  RETURN
END SUBROUTINE reverse_1

!  End module

End module GEMSTOOL_NSW_cross_sections_m

