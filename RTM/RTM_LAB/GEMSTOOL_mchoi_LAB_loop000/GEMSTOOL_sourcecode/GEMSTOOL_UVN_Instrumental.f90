Module GEMSTOOL_UVN_INSTRUMENTAL_m

!  Modules for Input

   USE GEMSTOOL_Input_types_m

!  Use modules
!  -----------

!  Use numerical Subroutine

   use GEMSTOOL_Numerical_m, ONLY : wtd_mean

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  10/26/16. Satellite data
!   Inst_maxlambdas = maximum channel wavelengths you think you will need.

   integer, parameter :: Inst_maxlambdas = 5001

!  Contains the following subroutine

!   1.  GEMS_Get_RSR                   ! Obtains Instrument RSR

!  everything public

public   :: GEMSTOOL_Instr_Lambdas, GEMSTOOL_Instr_RSR
private  :: fpk, Inst_maxlambdas

contains

subroutine GEMSTOOL_Instr_Lambdas ( &
   Maxwav, GEMSTOOL_Inputs,    & ! Input
   nwav, wav, fail, message)     ! Output

! Subroutine to generate wavelength grid for a given satellite band. 
!  R. Spurr 10/27/16

! INPUT/OUTPUT:
! Inputs      GEMSTOOL input type structure - this subroutine defines the
!             following inputs:
!             nlambda - Number of wavelengths and RSR points
!             ctrwl   - RSR-weighted band central wavelength
!             lambda  - Wavelength(s) (in nm) at which RSR is defined
!             rsr     - The RSR

! OUTPUT:
! wav,nwav    Fine-resolution Instrument wavelength grid
! fail        Fail warning flag
! message     Error message if failure occurs

implicit none

! Input

integer,                intent(in)    :: Maxwav

type(GEMSTOOL_Config_Inputs), intent(inout) :: GEMSTOOL_Inputs

! Output

integer,        intent(out) :: nwav
real(kind=fpk), intent(out) :: wav(Maxwav)

logical,           intent(out) :: fail
character (len=*), intent(out) :: message

!  Local

   integer   :: next, nbin, n, n1, ne
   real(fpk) :: wstart, wfinis, wext, wres, winit, wnum, &
                wstart_ext, wfinis_ext, range, step, step_ext

!  Initialize

   nwav = 0 ; wav = 0.0_fpk

!  Get instrument proxies

   wext       = GEMSTOOL_INPUTS%Instruments%band_Extension
   wstart     = GEMSTOOL_INPUTS%Instruments%band_start
   wfinis     = GEMSTOOL_INPUTS%Instruments%band_finish
   wres       = GEMSTOOL_INPUTS%Instruments%band_fineres

!  Number of wavelengths

   wstart_ext = wstart - wext
   wfinis_ext = wfinis - wext
   range = wfinis - wstart ; step = range / wres ; nbin = INT(step) + 1
   step_ext = wext / wres ; next = INT(step_ext)

   nwav = nbin + 2 * next
   if ( nwav .gt. maxwav) then
      message = 'Number of Satellite wavelengths > Dimension MAXWAV ==> Increase MAXWAV'
      fail    = .true.  ; return
   endif

!  Assign wavelengths

   winit = wstart_ext
   do n = 1, next
      n1 = nwav + 1 - n
      wnum = winit + wres * real(n-1,fpk) ; wav(n1) = 1.0d+07/wnum
   enddo
   winit = wstart
   do ne = 1, nbin
      n = ne + next ; n1 = nwav + 1 - n
      wnum = winit + wres * real(ne-1,fpk) ; wav(n1) = 1.0d+07/wnum
   enddo
   winit = wfinis
   do ne = 1, next
      n = ne + next + nbin; n1 = nwav + 1 - n
      wnum = winit + wres * real(ne,fpk) ; wav(n1) = 1.0d+07/wnum
   enddo

!  debug
!   do n = 1,nwav
!     write(88,*)n,wav(n),1.0d+07/wav(n)
!   enddo
!   pause'wavcheck'

!  Have to reassign the limits to wavelength space. for use in RSR work

   GEMSTOOL_INPUTS%Instruments%band_start  = 1.0d+07/wfinis
   GEMSTOOL_INPUTS%Instruments%band_finish = 1.0d+07/wstart

!  Done

   return
end subroutine GEMSTOOL_Instr_Lambdas

!

subroutine GEMSTOOL_Instr_RSR ( &
   Verbose, Inputs,    & !InOut
   fail, message)        !Output

! Subroutine to read the relative spectral response (RSR) for a given
! satellite band. Return it, and RSR-weighted centre wavelength.
! Use as part of constructing VLIDORT inputs.

! INPUT:
! Verbose     Activates debug output from the subroutine

! INPUT/OUTPUT:
! Inputs      GEMSTOOL input type structure - this subroutine defines the
!             following inputs:
!             nlambda - Number of wavelengths and RSR points
!             ctrwl   - RSR-weighted band central wavelength
!             lambda  - Wavelength(s) (in nm) at which RSR is defined
!             rsr     - The RSR

! OUTPUT:
! fail        Fail warning flag
! message     Error message if failure occurs

implicit none

! Input

logical,                intent(in)    :: Verbose

type(GEMSTOOL_Config_Inputs), intent(inout) :: Inputs

! Output

logical,           intent(out) :: fail
character (len=*), intent(out) :: message

! Local variables

integer :: band,errstat,fileunit,i,ic,j,nheaders,idum
real(kind=fpk) :: rdum,lam,tmplam_in(Inst_maxlambdas)
real(kind=fpk), allocatable, dimension(:)   :: tmprsr_in

integer            :: sensor
character (len=20) :: chsensor
character (len=2)  :: chband
character (len=80) :: rsr_path,fpath,fname

logical :: Quiet, DoOutput

integer        :: nlambda
real(kind=fpk) :: ctrwl
real(kind=fpk) :: lambda(inst_maxlambdas),rsr(inst_maxlambdas)

!logical :: Debug=.true.
logical :: Debug=.false.

!!!!!!!!!!!!!!!!!!!!!!!!!
! Define some variables !
!!!!!!!!!!!!!!!!!!!!!!!!!

rsr_path = trim('InstrumentFiles/RSRs/')

sensor   = Inputs%Instruments%Instrument_choice
chsensor = adjustl(Inputs%Instruments%Instrument_name)

band     = Inputs%Instruments%Instrument_Band
write(chband,'(i2.2)') band

if (Debug) then
  write(*,*)
  write(*,*) 'inside GEMSTOOL_Get_RSR:'
  write(*,*) 'sensor   =  ',sensor
  write(*,*) 'chsensor = |',chsensor,'|'
  write(*,*) 'band     =  ',band
  write(*,*) 'chband   = |',chband,'|'
endif

Quiet     = .not.Verbose

!DoOutput  = .false.
DoOutput  = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variable initialisation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

fail = .false.
message = ' '

ctrwl=0.0
lambda=-999.0_fpk
rsr=-999.0_fpk

fileunit=10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HIMAWARI-8. Filenames vary non-monotonically !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ( sensor .eq. 1 ) then
   ! HIMAWARI-8
   fname=' '

   ! Define filename
   if ( band >= 1  .and. band <= 9  ) fname='Himawari8_Band' // chband(2:2) // '.rsr.txt'
   if ( band >= 10 .and. band <= 16 ) fname='Himawari8_Band' // chband(1:2) // '.rsr.txt'

   fname=adjustl(fname)
   write(*,*)Trim(fname)

   if ( Debug ) write(*,*) '   Reading file = |',trim(fname),'|'

   ! Check for invalid band specification
   if ( fname == ' ' ) then
      message = 'HIMAWARI-8 band not recognised: ' // chband
      fail = .true. ; return
   endif

   ! Read the file
   allocate ( tmprsr_in(inst_maxlambdas) )
   fpath=trim(rsr_path) // 'himawari8/'
   nheaders=4 ! Number of header lines in file
   open ( unit=fileunit,&
          file=trim(fpath) // trim(fname),&
          status='old',&
          action='read' )
   ! Read header lines
   do i=1,nheaders
      read(fileunit,*)
   enddo
   if ( Debug ) write(*,*) '   HIMAWARI-8 band ',band
   i=1 ; ic = 0
   do
     if (mod(i,10).eq.0) then
        read(fileunit,*,iostat=errstat)lam,idum,rdum
     else
        read(fileunit,*,iostat=errstat)lam,rdum,rdum
     endif
     lam = lam * 1000.0d0
     if ( lam .ge. Inputs%Instruments%Band_start .and. &
          lam .le. Inputs%Instruments%Band_finish ) then
        ic = ic + 1 ; tmplam_in(ic) = lam ; tmprsr_in(ic) = rdum
     endif
     if ( errstat == -1 ) exit !EOF
     i=i+1
   enddo
   nlambda=ic
   close (fileunit)

   ! Get centre weighted wavelength
   ! For HIMAWARI, recommendation is only to use points with RSR > 0.01
   j=0
   do i=1,nlambda
     if (tmprsr_in(i) > 0.01_fpk) then
       j=j+1
       lambda(j)=tmplam_in(i)
       rsr(j)=tmprsr_in(i)
       if (Debug) then
         write(*,*) 'i = ',i,' j = ',j,' lambda(j) = ',lambda(j),' rsr(j) = ',rsr(j)
       end if
     endif
   enddo
   nlambda = j
   ctrwl=wtd_mean(nlambda,rsr,lambda)
   deallocate ( tmprsr_in )

   ! Write data to file for later plotting if requested
   if ( DoOutput ) then
      open ( unit=fileunit,&
             file='rsr_' // trim(chsensor) // '_band' // trim(chband) // '.txt',&
             status='replace',&
             action='write' )
      do i=1,nlambda
        write(fileunit,*) lambda(i),rsr(i)
      enddo
      close(fileunit)
   endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unknown sensor, return error !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ( ctrwl == 0.0_fpk ) then
   message = 'Sensor not recognised: ' // trim(chsensor)
   fail = .true. ; return
endif

!!!!!!!!!!!!!!!!!!!!!
! Define the output !
!!!!!!!!!!!!!!!!!!!!!

Inputs%Instruments%rsr_nlambda = nlambda
Inputs%Instruments%rsr_ctrwl   = ctrwl
Inputs%Instruments%rsr_lambda(1:nlambda) = lambda(1:nlambda)
Inputs%Instruments%rsr_vals(1:nlambda)   = rsr(1:nlambda)

end subroutine GEMSTOOL_Instr_RSR

!   End module

end Module GEMSTOOL_UVN_INSTRUMENTAL_m
