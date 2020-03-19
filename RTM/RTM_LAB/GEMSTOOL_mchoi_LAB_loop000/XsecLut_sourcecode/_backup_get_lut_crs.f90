module get_lut_crs_m


INTEGER, PARAMETER :: dp = KIND(1.0D0)

public
private :: dp

contains


SUBROUTINE get_lut_crs &
       ( Lut_path, the_molecule, nlambda, lambda, &
         is_wavenum, nz, ps, ts, &
         crs )

!SUBROUTINE get_lut_crs &
!        ( LUT_path, the_molecule, nwavnums, wavnums(1:nwavnums), &
!          is_wavenum, nlevels, level_press(0:nlayers), level_temps(0:nlayers), &
!          GAS_XSECS(1:nwavnums,0:nlayers,g) )
          

    IMPLICIT NONE

    ! Input arguments

    CHARACTER(LEN=*), INTENT(IN)                       :: Lut_path
    INTEGER, INTENT(IN)                                :: nlambda, nz
    CHARACTER (LEN=6), INTENT(IN)                      :: the_molecule
    LOGICAL, INTENT (IN)                               :: is_wavenum    ! *** T:wavenumber, F: nm ***
!    REAL(KIND=dp), INTENT(IN)                          :: fwhm         ! *** in cm^-1 or nm ***
    REAL(KIND=dp), DIMENSION(nlambda), INTENT(IN)      :: lambda
    REAL(KIND=dp), DIMENSION(nz), INTENT(IN)           :: ps, ts

    integer                                             :: i, j, nwn, wstart_read !mchoi
    REAL(KIND=dp)                                       :: dum, k

!    REAL(KIND=dp), DIMENSION(71)                    :: press
!    REAL(KIND=dp), DIMENSION(17,71)                 :: temp
!    REAL(KIND=dp), DIMENSION(:), allocatable        :: wn
!    REAL(KIND=dp), DIMENSION(:,:,:), allocatable    :: im_crs

!    REAL(KIND=dp), DIMENSION(71)                    :: diff_press


    !  Output arguments
    REAL(KIND=dp), DIMENSION(nlambda, nz), INTENT(OUT) :: crs
    REAL(KIND=dp), DIMENSION(31001, nz) :: crs_read !mchoi; maximum dimension:31001 for O2A/O2B

    !INTEGER, INTENT(OUT)                               :: errstat   ! Error    status
    !CHARACTER*(*), INTENT(OUT)                         :: message   ! output    message

    ! Local variables
    CHARACTER(LEN=2) :: molc
    CHARACTER(LEN=6) :: the_moleculeU, c6

    INTEGER          :: iz, np

    REAL(KIND=dp)    :: wstart, wend
!    CHARACTER(LEN=132)                       :: im_crs_filename, press_filename, temp_filename, wn_filename
!    CHARACTER(LEN=132)                       :: im_c_filename, p_filename, t_filename, w_filename
    CHARACTER(LEN=132)                       :: crs_filename
    CHARACTER(LEN=132)                       :: crs_filename1, crs_filename2, crs_filename3

    the_moleculeU = the_molecule

!    if (the_moleculeU .eq. 'H2O' ) nwn = 140003 
!    if (the_moleculeU .eq. 'CO2' .or. the_moleculeU .eq. 'CH4') nwn = 90002
!    if (the_moleculeU .eq. 'O2' ) nwn = 50001

     
    ! Initialize output
    crs   = 0.0d0
    crs_read = 0.0d0 !mchoi
    
    IF (is_wavenum) THEN
        wstart = lambda(1); wend = lambda(nlambda)
    ELSE
        wstart = 1.0D7/lambda(nlambda); wend = 1.0D7/lambda(1)
    ENDIF

    crs_filename1 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) // '_crs_band1_binary.out'
    crs_filename2 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) // '_crs_band2_binary.out'
    ! crs_filename3 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) // '_crs_band3_binary.out'

    ! Read Xsecs LUT files
    if ( wstart .lt. 13500 ) then 
        crs_filename = crs_filename1
        nwn = 31001 !---12920.00-13230.00; 0.01 wn res.  ;---O2A band
        wstart_read = 12920
    else 
        if ( wstart .gt. 13500  ) then 
            crs_filename = crs_filename2
            nwn = 25001 !--- 14340-14590; 0.01 wn res. ;--O2B band
            wstart_read = 14340
        ! else
        !     crs_filename = crs_filename3
        !     nwn = 14243 !---0.1 wn res. (range is not changed)
        endif
    endif

!    print*, wstart, wend, nlambda, nz
!    print*, crs_filename
! pause
!    open(11, file = crs_filename, status = 'unknown')
!    do iz = 1, nz
!        read(11, 50) crs(1:nlambda,iz)
!    enddo
!    close(11)
!50 format(100000e20.8)
    !print*,nwn,nz

    ! open(11, file = crs_filename, access = 'direct', form = 'unformatted', recl = nlambda*nz*4*dp)
    ! read(11, rec=1) ((crs(i,j),i=1,nlambda),j=1,nz)
    open(11, file = crs_filename, access = 'direct', form = 'unformatted', recl = nwn*nz*4*dp) !mchoi
    read(11, rec=1) ((crs_read(i,j),i=1,nwn),j=1,nz) !mchoi-modified
    close(11)
    
    ! put crs_read to crs for exact wn range (interval is identical as 0.01) !mchoi
    do i=1, nlambda
    crs(i,:)=crs_read(nint(100*(lambda(i)-wstart_read))+1,:)  
    ! print*,crs(1,1),crs_read((100*(lambda(i)-wstart_read))+1,1)
    ! print*,crs(1,1),crs(1.0,1)
    ! pause
    enddo

RETURN

END SUBROUTINE get_lut_crs

!   End module
end module get_lut_crs_m

