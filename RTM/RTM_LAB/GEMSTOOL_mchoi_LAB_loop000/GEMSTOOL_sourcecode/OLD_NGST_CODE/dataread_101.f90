MODULE dataread_101_m

!  Type definitions

  use vlidort_pars, only : fpk

private
public :: read_101prf, read_coe

CONTAINS

  SUBROUTINE read_101prf( fname, p, z, t, x, nlev)

    implicit none

    character*(*), intent(in) :: fname
    integer, intent(out) :: nlev
    integer, parameter :: maxlev=101
    real(fpk),    intent(out), dimension(:) :: p, z, t, x
    integer, parameter :: nhead=3
    integer :: i, iunit

    character(len=100) :: scrap
    integer :: ios

    iunit=9
    open(unit=iunit, file=fname, iostat=ios, &
         status="old", action="read")
    if( ios /= 0 ) then
       write(*,*)fname
       stop "Error trying to open file "
    end if

    ! Read header

    do i=1, nhead
       read(unit=iunit,fmt=*) scrap
       if( ios /= 0 ) then
          write(*,*)fname
          stop 'Error reading through header'
       end if
    end do

    do i = 1, maxlev
       ! Read data
       read(unit=iunit, fmt=*, iostat=ios) p(i), z(i), t(i), x(i)

       ! Stop reading if EOF reached
       if( ios < 0 ) exit

       ! Terminate program on read error.
       if( ios > 0 ) then
          write(*,'(A,I3)') "Error reading "//trim(fname)//" at record", i
          stop
       end if
    end do

    ! Read of profile file is complete. Close file, end.
    nlev = i-1
    close(unit=iunit)


  END SUBROUTINE read_101prf

  SUBROUTINE read_coe( fname, maxwav, wavdat, nwav, coedat)

    ! Allocate coedat array of size (maxwav, norder+1) before 
    ! passing in to this routine. e.g.:

    !    fname='coe.dat'
    !    norder=1
    !    allocate(coedat(maxwav,norder+1))
    !    call read_coe(fname, maxwav, wavdat, nwav, coedat)

    implicit none

    character*(*), intent(in) :: fname
    integer, intent(in) :: maxwav
    integer, intent(out) :: nwav
    real(fpk),    intent(out), dimension(:,:) :: coedat
    real(fpk),    intent(out), dimension(maxwav) :: wavdat
    integer :: unit
    integer, parameter :: nhead=1
    integer :: i

    character(len=100) :: scrap
    integer :: ios

    unit=9
    open(unit=unit, file=fname, iostat=ios, &
         status="old", action="read")
    if( ios /= 0 ) then
       stop 'Error trying to open file '
    end if

    ! Read header

    do i=1, nhead
       read(unit=unit, fmt=*, iostat=ios) scrap
       if( ios /= 0 ) then
          stop 'Error reading through '
       end if
    end do

    do i = 1, maxwav
       ! Read data 
       read(unit=unit, fmt=*, iostat=ios) wavdat(i), coedat(i,:)

       ! Stop reading if EOF reached
       if( ios < 0 ) exit

       ! Terminate program on read error.
       if( ios > 0 ) then
          write(*,'(A,I3)') "Error reading "//trim(fname)//" at record", i
          stop
       end if
    end do

    ! Read of coes file is complete. Close file, end.     
    nwav = i-1
    close(unit=9)

  END SUBROUTINE read_coe

END MODULE dataread_101_m
