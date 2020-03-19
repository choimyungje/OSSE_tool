module GEMSTOOL_Solarspec_m

!  Use numerical Subroutines

   use GEMSTOOL_Numerical_m, only : spline, splint

!  Original: GEMSTOOL UVN-specific Routine. 12 August 2013
!  Adapted for CODECTOOL, 29 September 2013
!  Adapted back GEMSTOOL UVN-specific Routine, 10/26/16

!  Precision

integer, parameter :: fpk = SELECTED_REAL_KIND(15)

private :: fpk
public  :: GEMSTOOL_Solarspec

contains

subroutine GEMSTOOL_Solarspec ( &
   MaxLambdas, Physics_DataPath, SunfileName, nlambdas, Lambdas, & ! Input
   SolarSpec, fail, message )                                      ! Output

   implicit none

!  INPUTS
!  ======

   integer  , intent(in) :: MaxLambdas

!  path and filename

   character*(*), intent(in) :: Physics_DataPath
   character*(*), intent(in) :: SunfileName

!  Wavelengths

   integer  , intent(in) :: nlambdas
   real(fpk), intent(in) :: Lambdas(MaxLambdas)

!  OUTPUTS
!  =======

!  Solar Spectrum

   Real(fpk), intent(out) :: SolarSpec(MaxLambdas)

!  Exception handling

    logical      , intent(out) :: fail
    character*(*), intent(out) :: message

!  LOCAL
!  =====

!mick fix 12/10/2013 - added flexible dimensioning for SUNWAV & SUNFLUX

   integer         :: I, NI, n, w, nwav
   real(kind=fpk)  :: X, Y
   real(kind=fpk)  :: INPUT_WAVNUM, INPUT_FLUX, wavnum_start, wavnum_finis
   real(kind=fpk)  :: INPUT_WAVNUM_PREV, INPUT_FLUX_PREV, wavstart, wavfinis
   real(kind=fpk), parameter :: small = 1.0d-30

   real(kind=fpk), allocatable :: SUNWAV(:),SUNFLUX(:)
   real(kind=fpk), allocatable :: SUNX(:),SUNY(:),SUN2D(:)

!  Initialize output

   fail    = .false.
   message = ' '
   SolarSpec = 0.0_fpk

!  Set limits (converting wvl in nm to wvn in cm^-1)

   wavstart = Lambdas(1) ; wavfinis = Lambdas(nLambdas)
   wavnum_start = 1.0d+07/wavfinis - 2.0
   wavnum_finis = 1.0d+07/wavstart + 2.0

!  NEWKUR.DAT : Trawl through data, using one point extra at both ends

   n = 0
   OPEN(1,FILE =trim(Physics_DataPath)//Trim(SunfileName), err=60, status = 'old')

!mick fix 12/10/2013 - added 1st line; moved and enhanced the 2nd
   nwav = INT(wavnum_finis - wavnum_start) + 3
   allocate ( SUNWAV(nwav),SUNFLUX(nwav),SUNX(nwav),SUNY(nwav),SUN2D(nwav) )

   READ(1,*) ; READ(1,*)
   DO I = 1, 60000
      READ(1,*)INPUT_WAVNUM, INPUT_FLUX
      if ( INPUT_WAVNUM .GE. wavnum_start ) then
         if ( INPUT_WAVNUM .LE. wavnum_finis ) then
            if ( n.eq.0 ) then
               n = n + 1
               sunwav(N)  = INPUT_WAVNUM_PREV
               sunflux(N) = INPUT_FLUX_PREV
            endif
            n = n + 1
            sunwav(N)  = INPUT_WAVNUM
            sunflux(N) = INPUT_FLUX
         else
            n = n + 1
            sunwav(N)  = INPUT_WAVNUM
            sunflux(N) = INPUT_FLUX
            go to 56
         endif
      endif
      INPUT_WAVNUM_PREV = INPUT_WAVNUM
      INPUT_FLUX_PREV   = INPUT_FLUX
   ENDDO
56 continue
   CLOSE(1)

!mick fix 12/10/2013 - added error check
   if (N > NWAV) then
     write(*,*)
     write(*,*) 'NWAV = ',NWAV
     write(*,*) 'N    = ',N
     deallocate ( SUNWAV,SUNFLUX,SUNX,SUNY,SUN2D )
     goto 61
   end if

   NWAV = N

!mick fix 3/10/2015 - converted solar flux to nm space.
!  Convert solar grid and solar flux values just obtained:
!  * Convert solar grid from wavenumbers in cm^-1 to wavelengths in nm and flip array
!  * Convert solar spectral flux avg values from W/(cm^2*cm^-1) to W/(cm^2*nm) and flip array

   do i = 1, NWAV
      ni = NWAV + 1 - i
      SUNX(ni) = 1.0d+07/SUNWAV(I)
      SUNY(ni) = SUNFLUX(I)*(SUNWAV(I)/SUNX(ni))
   enddo

!  Spline

   Call Spline(NWAV,SUNX,SUNY,NWAV,SMALL,SMALL,SUN2D)
   DO w = 1, nlambdas
      X = Lambdas(w)
      CALL SPLINT (NWAV,SUNX,SUNY,SUN2D, NWAV, X, Y)
      SolarSpec(w) = Y
   ENDDO

!  Deallocate

   deallocate ( SUNWAV,SUNFLUX,SUNX,SUNY,SUN2D )

!  Had to reset this, 2/23/16

   FAIL = .false.

!  Normal return

   return

!  Error finishes

60 continue ; close(1)
   message = 'Open file failure for solar spectrum'
   fail = .true. ; return

61 continue
   message = 'Dimension of solar arrays in GEMSTOOL_Solarspec inadequate'
   fail = .true. ; return

!  Finish

   return
end subroutine GEMSTOOL_Solarspec

!  End module

end module GEMSTOOL_Solarspec_m

