! Comments from Kelly's original code
! read logicals and parameters controlling calculations
!
! molnum is the hitran molecule number.
!
! start, step, and npoints define the calculation grid. in order to avoid
! undersampling the grid should be at least as fine as 1/3 of the smallest
! gaussian fwhm of a spectral line, or 0.5550 times the smallest gaussian hw1e.
! this insures that the maximum sampling error is <1.24e-4 of the full-scale
! line-shape. later: add automatic check for undersampling wrt voigt widths
!
! press is the pressure in atmospheres (millibars / 1013.25). temp is the
! temperature in degrees kelvin.
!
! nvoigt is the number of grid points to each side of a spectral line for
! performing the voigt calculation. assuming calculation to <= 1e-4 is desired,
! nvoigt should be the greater of (1) 100 * hwhm_max / step, where hwhm_max is
! the largest lorentzian hwhm of a line and step is the grid spacing;
! (2) 3.035 * hw1e_max / step, where hw1e_max is the largest gaussian hw1e of a
! line.
!
! hw1e is the gaussian slit with at 1/e intensity. hw1e = the full width at
! half-maximum, fwhm, * 0.60056 = fwhm / (2 * sqrt (ln 2)).
!
! nmod gives the option not to write out humongous spectral files by printing
! every nmod^th spectral value

! Notes: 
! 1. Qpower needs to be more accurate
! 2. pressure-induced shift is not considered
! 3. Intensity for different isotopes is weighted by their fraction in the atmosphere
! 4. Is nvoigt enough?
! 5. Need to speed up Voigt calculations?
! 6. Sinc function?

! Updates
! Feb. 2012: add HUMLIK function, which improves the calculation by a factor of 5 compared to voigt
! Feb. 25, 2012: add crsdt
! Aug. 27, 2012: Update Kelly's new partition function

! RT Solutions 12 March 2013
!   * Make into a module
!   * Introduce HITRAN path as variable

! RT Solutions 07 July 2013
!   * Linearized Version of code
!   * Modularize code more, add exception handling

module get_lut_crs_plus_m

INTEGER, PARAMETER :: dp = KIND(1.0D0)

public
private :: dp

contains

SUBROUTINE get_lut_crs_plus &
       ( LUT_path, the_molecule, nlambda, lambda,   &
         is_wavenum, nz, ps, ts, do_dP, do_dT,      &
         crs, crsdt, crsdp)
         !crs, crsdt, crsdp, errstat, message )

IMPLICIT none

! Input arguments

CHARACTER(LEN=*), INTENT(IN)                       :: LUT_path !! Y.Jung
!CHARACTER(LEN=*), INTENT(IN)                       :: Hitran_path  !  Added by R. Spurr
CHARACTER (LEN=6), INTENT(IN)                      :: the_molecule
INTEGER, INTENT(IN)                                :: nlambda, nz
LOGICAL, INTENT (IN)                               :: is_wavenum    ! *** T: wavenumber, F: nm ***
!REAL(KIND=dp), INTENT(IN)                           :: fwhm          ! *** in cm^-1 or nm ***
REAL(KIND=dp), DIMENSION(nlambda), INTENT(IN)       :: lambda
REAL(KIND=dp), DIMENSION(nz), INTENT(IN)            :: ps, ts
LOGICAL, INTENT(IN)                                :: do_dP, do_dT  ! Flags for P and T derivatives

integer                                             :: i, j, nwn, wstart_read !mchoi
REAL(KIND=dp)                                       :: dum, k

!  Output arguments

!INTEGER, INTENT(OUT)                               :: errstat   ! Error status     !! Y.Jung
!CHARACTER*(*), INTENT(OUT)                         :: message   ! output message   !! Y.Jung   
REAL(KIND=dp), DIMENSION(nlambda, nz), INTENT(OUT)  :: crs
REAL(KIND=dp), DIMENSION(nlambda, nz), INTENT(OUT)  :: crsdt
REAL(KIND=dp), DIMENSION(nlambda, nz), INTENT(OUT)  :: crsdp


REAL(KIND=dp), DIMENSION(31001, nz) :: crs_read !mchoi; maximum dimension:31001 for O2A/O2B
REAL(KIND=dp), DIMENSION(31001, nz) :: crsdt_read !mchoi; maximum dimension:31001 for O2A/O2B
REAL(KIND=dp), DIMENSION(31001, nz) :: crsdp_read !mchoi; maximum dimension:31001 for O2A/O2B


! Local arrays
! ------------

! --Original------->

!INTEGER, DIMENSION(maxlines)             :: mol, iso
!REAL(KIND=dp), DIMENSION(maxlines)        :: sigma0, strnth, einstein, alpha, &
!                                            elow, coeff, selfbrdn, pshift
!REAL(KIND=dp), DIMENSION(maxmols, maxiso) :: amu, q296, q
!LOGICAL, DIMENSION(maxmols, maxiso)      :: if_q
!REAL(KIND=dp), DIMENSION(maxpoints)       :: pos, spec, voigtx, v, posnm

!  @@ RTS, commented out not used
!REAL(KIND=dp), DIMENSION(maxpoints)       :: specmod
!REAL(KIND=dp), DIMENSION(maxmols)         :: qpower
!INTEGER          :: j

!  --Linearized------->

!REAL(KIND=dp), DIMENSION(maxmols, maxiso) :: dqdT
!REAL(KIND=dp), DIMENSION(maxpoints)       :: dT_voigtx, dP_voigtx, dT_spec, dP_spec
!REAL(KIND=dp), DIMENSION(maxpoints,2)     :: dX, dv

!  Local variables
!  ---------------

!  @@ RTS, Make sure file string is long enough, add q_file

CHARACTER(LEN=2) :: molc
CHARACTER(LEN=6) :: the_moleculeU, c6

INTEGER          :: iz, np

REAL(KIND=dp)    :: wstart, wend

!CHARACTER(LEN=132)                       :: hitran_filename, q_file
CHARACTER(LEN=132)                       :: crs_filename, crsdt_filename, crsdp_filename        !! Y.Jung
CHARACTER(LEN=132)                       :: crs_filename1, crs_filename2, crs_filename3         !! Y.Jung
CHARACTER(LEN=132)                       :: crsdt_filename1, crsdt_filename2, crsdt_filename3   !! Y.Jung
CHARACTER(LEN=132)                       :: crsdp_filename1, crsdp_filename2, crsdp_filename3   !! Y.Jung

! (Original)

!INTEGER          :: molnum, npoints, nvoigt, ntemp
!REAL(KIND=dp)     :: wstart, wend, step, press, temp, minline_vg, &
!                    minslit_fwhm, min_fwhm, maxline_vg, maxline_hwhm, voigt_extra
!
!INTEGER          :: i, mol_temp, iso_temp, nlines, nvlo, nvhi, idx, MEND, iz
!REAL(KIND=dp)     :: sigma0_temp, strnth_temp, einstein_temp, alpha_temp, &
!                    selfbrdn_temp, elow_temp, coeff_temp, pshift_temp, sigma_temp
!REAL(KIND=dp)     :: vg, voigta, ratio1, ratio2, ratio, vnorm, rt0t, rc2t, rc2t0
!CHARACTER(LEN=2) :: molc
!CHARACTER(LEN=6) :: the_moleculeU, c6

!  Linearized

!INTEGER          :: m, mt, mp, nv, moli, isoi
!REAL(KIND=dp)     :: vg_00, dP_vg, dT_vg, voigta_00, term1, term2, term3, term4, term5, term6, dY(2)
!REAL(KIND=dp)     :: dP_term3, dT_term1, dT_term2, dT_term3, dP_term5
!REAL(KIND=dp)     :: dP_ratio1, dP_ratio2, dP_ratio, dP_vnorm, dP_voigta, dT_rt0t, dT_rc2t
!REAL(KIND=dp)     :: dT_ratio1, dT_ratio, dT_vnorm, dT_voigta, amwt

! Initialize error status

!errstat = 0 ; message = ' '

!  Initialize output

crs   = 0.0d0
crsdp = 0.0d0
crsdt = 0.0d0

crs_read = 0.0d0 !mchoi
crsdp_read = 0.0d0 !mchoi
crsdt_read = 0.0d0 !mchoi

! write(*,*) nlambda
! pause


!  Determine number of derivatives
!   m = 0 (no derivatives), m = 1, Temperature or pressure only, m = 2 (both)

!m = 0
!if ( do_dT ) m = m + 1 ; mt = m
!if ( do_dP ) m = m + 1 ; mp = m

! Determine which molecule and database file

the_moleculeU = the_molecule !StrUpCase(the_molecule)

    IF (is_wavenum) THEN
        wstart = lambda(1); wend = lambda(nlambda)
    ELSE
        wstart = 1.0D7/lambda(nlambda); wend = 1.0D7/lambda(1)
    ENDIF

    crs_filename1 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) //'_crs_band1_binary.out'
    crs_filename2 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) //'_crs_band2_binary.out'
    ! crs_filename3 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) //'_crs_band3_binary.out'
    
    crsdt_filename1 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) //'_dXsecs_dT_band1_binary.out'
    crsdt_filename2 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) //'_dXsecs_dT_band2_binary.out'
    ! crsdt_filename3 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) //'_dXsecs_dT_band3_binary.out'
    
    crsdp_filename1 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) //'_dXsecs_dP_band1_binary.out'
    crsdp_filename2 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) //'_dXsecs_dP_band2_binary.out'
    ! crsdp_filename3 = adjustl(trim(Lut_path)) // adjustl(trim(the_moleculeU)) //'_dXsecs_dP_band3_binary.out'

     ! Read Xsecs LUT files
     if ( wstart .lt. 13500 ) then
         crs_filename = crs_filename1
         crsdt_filename = crsdt_filename1
         crsdp_filename = crsdp_filename1
         nwn = 31001 !---12920.00-13230.00; 0.01 wn res.  ;---O2A band
         wstart_read = 12920
     else
         if ( wstart .gt. 13500  ) then
            crs_filename = crs_filename2
            crsdt_filename = crsdt_filename2
            crsdp_filename = crsdp_filename2
            nwn = 25001 !--- 14340-14590; 0.01 wn res. ;--O2B band
            wstart_read = 14340
        ! else
        !     crs_filename = crs_filename3
        !     crsdt_filename = crsdt_filename3
        !     crsdp_filename = crsdp_filename3
        !     nwn = 14243 !---0.1 wn res. (range is not changed)   
        endif
    endif


    ! open(11, file = crs_filename, access = 'direct', form = 'unformatted', recl = nlambda*nz*4*dp)
    ! read(11, rec=1) ((crs(i,j), i=1,nlambda), j=1,nz)
    open(11, file = crs_filename, access = 'direct', form = 'unformatted', recl = nwn*nz*4*dp)
    read(11, rec=1) ((crs_read(i,j),i=1,nwn),j=1,nz) !mchoi-modified
    close(11)

    ! open(12, file = crsdt_filename, access = 'direct', form = 'unformatted', recl = nlambda*nz*4*dp)
    ! read(12, rec=1) ((crsdt(i,j), i=1,nlambda), j=1,nz)
    open(12, file = crsdt_filename, access = 'direct', form = 'unformatted', recl = nwn*nz*4*dp)
    read(12, rec=1) ((crsdt_read(i,j), i=1,nwn), j=1,nz)
    close(12)

    ! open(13, file = crsdp_filename, access = 'direct', form = 'unformatted', recl = nlambda*nz*4*dp)
    ! read(13, rec=1) ((crsdp(i,j), i=1,nlambda), j=1,nz)
    open(13, file = crsdp_filename, access = 'direct', form = 'unformatted', recl = nwn*nz*4*dp)
    read(13, rec=1) ((crsdp_read(i,j), i=1,nwn), j=1,nz)
    close(13)


    ! put crs_read to crs for exact wn range (interval is identical as 0.01) !mchoi
    do i=1, nlambda
        crs(i,:)=crs_read(nint(100*(lambda(i)-wstart_read))+1,:)  
        crsdt(i,:)=crsdt_read(nint(100*(lambda(i)-wstart_read))+1,:)  
        crsdp(i,:)=crsdp_read(nint(100*(lambda(i)-wstart_read))+1,:)  
        ! print*,crs(1,1),crs_read((100*(lambda(i)-wstart_read))+1,1)
        ! print*,crs(1,1),crs(1.0,1)
        ! pause
    enddo

!    open(11, file = crs_filename, status = 'unknown')
!    do iz = 1, nz
!        read(11, 50) crs(1:nlambda,iz)
!    enddo
!    close(11)

!    open(12, file = crsdt_filename, status = 'unknown')
!    do iz = 1, nz
!        read(12, 50) crsdt(1:nlambda,iz)
!    enddo
!    close(12)

!    open(13, file = crsdp_filename, status = 'unknown')
!    do iz = 1, nz
!        read(13, 50) crsdp(1:nlambda,iz)
!    enddo
!    close(13)

!50 format(100000e20.8)

RETURN

END SUBROUTINE  get_lut_crs_plus

!  End mmdule

END MODULE  get_lut_crs_plus_m
