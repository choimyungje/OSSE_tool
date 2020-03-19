module GEMSTOOL_UVN_Convolution_m

!  New convolution routine, for output. 10/27/16

!  Inclusion of the XSECS_Convolution routine
!    Adapted from CODEC usage by R. Spurr, 10/26/16

!  Use modules
!  -----------

!  Use Tpe structures for inputs

   use GEMSTOOL_Input_Types_m

!  Use Instrument routine

   use GEMSTOOL_UVN_INSTRUMENTAL_m

!  Use numerical Subroutines

   use GEMSTOOL_numerical_m   , only : RSR_to_fine_grid, wtd_mean, asym_gauss_f2c

public
private :: fpk

contains

subroutine GEMSTOOL_Convolution &
   (  MAXWAV, MAXWAV_c, MAX_GEOMETRIES, do_Fluxes, Sun_normalized, & ! GEMSTOOL control
      nwavs, nwavs_C, nstokes, n_geometries, wavs, wavs_C, Inputs, & ! Numbers and wavelengths
      sunspec, STOKES, ACTINIC, REGFLUX,                           & ! Main program, GEMSTOOL Results
      sunspec_C, STOKES_C, ACTINIC_C, REGFLUX_C, DOLP_C )            ! Coarse Grid convolved output

   implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  GEMSTOOL wavelength a dimensions

   integer, intent(in) :: MAXWAV, MAXWAV_C, MAX_GEOMETRIES

!  Flux control flag (GEMSTOOL control = do_SphericalAlbed)

   logical, intent(in) :: do_Fluxes

!  SolarSpectrum control flag, should be set here

   logical, intent(in) :: Sun_Normalized

!  Proxy control inputs

   integer, intent(in) :: nstokes, n_geometries

!  Wavelength inputs

   integer  , intent(in)  :: nwavs, nwavs_C
   real(fpk), intent(in)  :: wavs  (MaxWav)
   real(fpk), intent(in)  :: wavs_C(MaxWav_C)

!  Inputs

   type(GEMSTOOL_Config_Inputs), intent(in) :: INPUTS

!  Fine grid arrays
!     Solar, Stokes Vector, Fluxes, Degrees of linear/circular polarization (DOLP/DOCP)

   real(fpk)   , DIMENSION (MAXWAV)                       :: SUNSPEC
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAV)      :: STOKES
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAV)      :: ACTINIC
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAV)      :: REGFLUX

!  Output
!  ------

   real(fpk)   , DIMENSION (MAXWAV_C)                       :: SUNSPEC_C
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAV_C)      :: STOKES_C
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV_C)        :: DOLP_C
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV_C)        :: DOCP_C
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAV_C)      :: ACTINIC_C
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAV_C)      :: REGFLUX_C

!  Local

   real(fpk)           :: sq2u2,  zero
   integer             :: w, j, nf, nc, nr, nspec, ns

   real(kind=fpk), dimension(:)  , allocatable :: fwave, cwave, hwle, asym, rwave, rsr, frsr
   real(kind=fpk), dimension(:,:), allocatable :: fsun, csun, Vector1, Vector2, Vector1_C, Vector2_C

!  Local flags, revised usage 5/13/16

   logical :: DoGaussianSlit

!  new definitions for Brightness temperature. 5/19/16  Rob add
!   C1 has units W/m^2/sr/um^4 ;  C2 has units K um
!   real(fpk), parameter  :: C1_bt = 1.191042e+08_fpk
!   real(fpk), parameter  :: C2_bt = 1.4387752e+04_fpk

!  Initialization

   zero = 0.0_fpk
   SUNSPEC_C = ZERO
   STOKES_C  = ZERO
   ACTINIC_C = ZERO
   REGFLUX_C = ZERO
   DOLP_C    = ZERO
   DOCP_C    = ZERO

!  CONVOLUTION OF RADIANCES (reflectance * irradiance) WITH SLIT FUNCTIONS
!  =======================================================================

!  Proxies

   nspec = n_geometries
   nf    = nwavs

!  Allocations

   allocate ( fwave(nf),fsun(nf,1),Vector1(nf,nspec),Vector2(nf,nspec) )

!  Assign Fwave & Fsun

   fwave     = Wavs(1:nf) ;  fsun(:,1) = Sunspec(1:nf)

!  Assign Slit Function type.
!       Gaussian  Slit is now defunct. Does not depend on thermal flag. (commented out, 5/13/16)

   DoGaussianSlit =.false.

!  Perform Old-style gaussian slit convolution
!  ===========================================

   if (DoGaussianSlit) then

     !Assign Gaussian slit function
     nc    = nwavs_c

     allocate ( cwave(nc),asym(nc),hwle(nc),csun(nc,1) )
     allocate ( Vector1_C(nc,nspec),Vector2_C(nc,nspec) )

     cwave = Wavs_c(1:nc)
     asym  = zero
     hwle  = Inputs%Instruments%band_slitpars(1)

     !  Convolve Irradiances to the coarse lambda grid

     CALL asym_gauss_f2c ( fwave, fsun, nf, 1, hwle, asym, cwave, csun, nc )
     Sunspec_C(1:nc) = csun(1:nc,1)

     IF ( Sun_Normalized ) then
        do ns = 1, nstokes
           do j = 1, nspec
              Vector1(1:nf,j) = Stokes  (j,ns,1:nf) * Sunspec(1:nf)
           enddo
           CALL asym_gauss_f2c ( fwave, Vector1   (1:nf,1:nspec), nf, nspec, hwle, asym, &
                                 cwave, Vector1_C (1:nc,1:nspec), nc )
           do j = 1, nspec
             Stokes_C(j,ns,1:nc) = Vector1_C(1:nc,j) / Sunspec_C(1:nc) 
           enddo
        enddo
     else
        do ns = 1, nstokes
           do j = 1, nspec
              Vector1(1:nf,j) = Stokes  (j,ns,1:nf)
           enddo
           CALL asym_gauss_f2c ( fwave, Vector1   (1:nf,1:nspec), nf, nspec, hwle, asym, &
                                 cwave, Vector1_C (1:nc,1:nspec), nc )
           do j = 1, nspec
             Stokes_C(j,ns,1:nc) = Vector1_C(1:nc,j)
           enddo
        enddo
     endif

!  Flux Output: Loop over solar_angles and Number of stokes parameters

     if ( do_Fluxes ) then
        IF ( Sun_Normalized ) then
           do ns = 1, nstokes
              do j = 1, nspec
                 Vector1(1:nf,j) = Actinic (j,ns,1:nf) * Sunspec(1:nf)
                 Vector2(1:nf,j) = Regflux (j,ns,1:nf) * Sunspec(1:nf)
              enddo
              CALL asym_gauss_f2c ( fwave, Vector1   (1:nf,1:nspec), nf, nspec, hwle, asym, &
                                    cwave, Vector1_C (1:nc,1:nspec), nc )
              CALL asym_gauss_f2c ( fwave, Vector2   (1:nf,1:nspec), nf, nspec, hwle, asym, &
                                    cwave, Vector2_C (1:nc,1:nspec), nc )
              do j = 1, nspec
                Actinic_C(j,ns,1:nc) = Vector1_C(1:nc,j) / Sunspec_C(1:nc) 
                Regflux_C(j,ns,1:nc) = Vector2_C(1:nc,j) / Sunspec_C(1:nc) 
              enddo
           enddo
        else
           do ns = 1, nstokes
              do j = 1, nspec
                 Vector1(1:nf,j) = Actinic (j,ns,1:nf)
                 Vector2(1:nf,j) = Regflux (j,ns,1:nf)
              enddo
              CALL asym_gauss_f2c ( fwave, Vector1   (1:nf,1:nspec), nf, nspec, hwle, asym, &
                                    cwave, Vector1_C (1:nc,1:nspec), nc )
              CALL asym_gauss_f2c ( fwave, Vector2   (1:nf,1:nspec), nf, nspec, hwle, asym, &
                                    cwave, Vector2_C (1:nc,1:nspec), nc )
              do j = 1, nspec
                Actinic_C(j,ns,1:nc) = Vector1_C(1:nc,j)
                Regflux_C(j,ns,1:nc) = Vector2_C(1:nc,j) 
              enddo
           enddo
        endif
     endif

     deallocate ( cwave,asym,hwle,csun )
     deallocate ( Vector1_C, Vector2_C )

   else

     !Assign instrument empirical slit function
     nr    = Inputs%Instruments%rsr_nlambda ; allocate ( rwave(nr),rsr(nr),frsr(nf) )
     rwave = Inputs%Instruments%rsr_lambda(1:nr)
     rsr   = Inputs%Instruments%rsr_vals(1:nr)

     !Interpolate RSR to the fine spectral grid
     CALL RSR_to_fine_grid ( rwave, rsr, fwave, frsr )

     !Convolve Irradiances to the coarse lambda grid

     Sunspec_C(1) = wtd_mean(nf,frsr,fsun(:,1))

     IF ( Sun_Normalized ) then
        do ns = 1, nstokes
           do j = 1, nspec
              Vector1(1:nf,j)  = Stokes  (j,ns,1:nf) * Sunspec(1:nf)
              Vector1_C(1,j)   = wtd_mean(nf,frsr,Vector1(1:nf,j))
              Stokes_C(j,ns,1) = Vector1_C(1,j) / Sunspec_C(1) 
           enddo
        enddo
     else
        do ns = 1, nstokes
           do j = 1, nspec
              Stokes_C(j,ns,1) = wtd_mean(nf,frsr,Stokes (j,ns,1:nf))
           enddo
        enddo
     endif

!  Flux Output: Loop over solar_angles and Number of stokes parameters

     if ( do_Fluxes ) then
        IF ( Sun_Normalized ) then
           do ns = 1, nstokes
              do j = 1, nspec
                 Vector1(1:nf,j)   = Actinic (j,ns,1:nf) * Sunspec(1:nf)
                 Vector2(1:nf,j)   = Regflux (j,ns,1:nf) * Sunspec(1:nf)
                 Vector1_C(1,j)    = wtd_mean(nf,frsr,Vector1(1:nf,j))
                 Vector2_C(1,j)    = wtd_mean(nf,frsr,Vector2(1:nf,j))
                 Actinic_C(j,ns,1) = Vector1_C(1,j) / Sunspec_C(1) 
                 Regflux_C(j,ns,1) = Vector2_C(1,j) / Sunspec_C(1) 
              enddo
           enddo
        else
           do ns = 1, nstokes
              do j = 1, nspec
                 Actinic_C(j,ns,1) = wtd_mean(nf,frsr,Actinic(j,ns,1:nf))
                 Regflux_C(j,ns,1) = wtd_mean(nf,frsr,Regflux(j,ns,1:nf))
              enddo
           enddo
        endif
     endif

     deallocate ( rwave,rsr,frsr )
     deallocate ( Vector1_C, Vector2_C )

   end if

!  Set DOLP to zero if no polarization (Nstokes = 1), and skip DOLP calculation
!         degree of polarization (only for NSTOKES = 3 or 4)

   IF ( NSTOKES.gt.1 ) THEN
      do w = 1, nc
         do j = 1, nspec
            SQ2U2 = SQRT ( Stokes_C(j,2,w) * Stokes_C(j,2,w) + Stokes_C(j,3,w) * Stokes_C(j,3,w) )
            DOLP_C(J,W) = SQ2U2 / Stokes_C(j,1,w)
            if ( nstokes .eq. 4 ) then
               DOCP_C(J,W) = SQRT( Stokes_C(j,4,w) * Stokes_C(j,4,w) ) / Stokes_C(j,1,w)
            endif
         enddo
      enddo
   ENDIF

!  Final de-allocation

   deallocate ( fwave,fsun, Vector1,Vector2 )

!  End subroutine

   return
end subroutine GEMSTOOL_Convolution

!

subroutine GEMSTOOL_XSECS_CONVOLVE &
          ( MAXLAMBDAS, MAXLAYERS, MAXGASES, Verbose,    & ! INPUT Dimensioning
            NLAYERS, NGASES, WHICH_GASES, INPUTS,        & ! Inputs
            NLAMBDAS, LAMBDAS, SUNSPEC,                  & ! Modified OUTPUT
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL, GAS_XSECS,    & ! Modified OUTPUT
            O3C1_XSECS, O3C2_XSECS, H2O_XSECS, O2_XSECS )  ! Modified OUTPUT

!  New routine, written by R. Spurr, 01 November 2014
!  Solar Spectrum optional R. Spurr, 13 May      2016
!  Effective wavelength calculation by R. Spurr, 31 August 2016.

!  Adapted for GEMSTOOL, 26 October 2016.

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input arguments
!  ---------------

!  Dimensioning

      INTEGER     :: MAXLAMBDAS, MAXLAYERS, MAXGASES

!  Verbose flag

      logical, intent(in) :: Verbose

!  wavelengths, Solar spectrum
 
      INTEGER     :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 
      real(fpk),    dimension ( MAXLAMBDAS ) :: SUNSPEC

!  Layer and gas control

      INTEGER                                  :: NLAYERS
      integer                                  :: NGASES
      character(LEN=4), dimension ( maxgases ) :: WHICH_GASES

!  Inputs

      type(GEMSTOOL_Config_Inputs), intent(inout) :: INPUTS

!  Output arguments
!  ----------------

!  Rayleigh cross-sections and depolarization output

      real(fpk),    dimension ( MAXLAMBDAS ) :: RAYLEIGH_XSEC 
      real(fpk),    dimension ( MAXLAMBDAS ) :: RAYLEIGH_DEPOL

!  Gas cross-sections. ON THE LEVEL

      REAL(fpk),    DIMENSION( MAXLAMBDAS, maxgases ) :: GAS_XSECS  

!  O3 cross-sections temperature-fitting coefficients C1 nad C2

      REAL(fpk),    DIMENSION( MAXLAMBDAS ) :: O3C1_XSECS  
      REAL(fpk),    DIMENSION( MAXLAMBDAS ) :: O3C2_XSECS  

!  H2O and O2

      REAL(fpk),    DIMENSION ( MAXLAMBDAS, 0:maxlayers ) :: H2O_XSECS  
      REAL(fpk),    DIMENSION ( MAXLAMBDAS, 0:maxlayers ) :: O2_XSECS  

!  Local
!  =====

   integer           :: nlevels
   integer           :: nf, nc, nr, g, n, n1, nclambdas
   real(kind=fpk)    :: scalex ( maxgases ), scaleO3_c1c2, clambdas(1)

!  8/31/16. frsrfsun added

   real(kind=fpk), dimension(:)    , allocatable :: fwave, cwave, hwle, asym, rwave, rsr, frsr, frsrfsun
   real(kind=fpk), dimension(:,:)  , allocatable :: fsun, csun, fRay, CRay, fabs1, cabs1, fabs, cabsM
   real(kind=fpk), dimension(:)    , allocatable :: CO3C1, CO3C2
   real(kind=fpk), dimension(:,:)  , allocatable :: CXsec1, CH2O_M, CO2_M

!  Local flags

   logical :: DoGaussianSlit
   logical :: Do_Effective_wavelength

!  CONVOLUTION OF Cross-sections (Using the solar I0-effect) WITH SLIT FUNCTIONS
!  =============================================================================

!  Proxies

   nf = nlambdas

!  8/31/16. Hard-wire the effective wavelength calculation = TRUE for any use of SOLAR SPRECTRUM

   Do_Effective_wavelength = .true.

! allocate  fine grid arrays

   allocate ( fwave(nf), fsun(nf,1), fRay(nf,2), fabs1(nf,1), fabs(nf,nlevels) )

!  Assign Fwave & Fsun

   fwave     = lambdas(1:nf)
   fsun(:,1) = Sunspec(1:nf)

!  Assign Slit Function type.
!       Gaussian  Slit is now defunct. Does not depend on thermal flag. (commented out, 5/13/16)

   DoGaussianSlit=.false.

!  Initialize coarse wavelength here inside the routine

   nclambdas = 1
   clambdas  = 0.0_fpk

!  Perform convolution with Gaussian slit function (assumed symmetrical)
!  =====================================================================

   if (DoGaussianSlit) then

     !Assign Gaussian slit function and define cwave

     nc    = nclambdas ; allocate ( cwave(nc),asym(nc),hwle(nc) )
     cwave(1) = Inputs%Instruments%band_median
     asym     = 0.0_fpk
     hwle     = Inputs%Instruments%band_slitpars(1)

     ! Allocate physical quantities

     allocate ( csun(nc,1), Cray(nc,2), cabs1(nc,1), cabsM(nc,nlevels), CO3C1(nc), CO3C2(nc), &
                CXSec1(nc,ngases), CH2O_M(nc,nlevels), CO2_M(nc,nlevels) )

     ! Convolve Irradiances to the coarse lambda grid. Now optional 5/13/16

     CALL asym_gauss_f2c ( fwave, fsun, nf, 1, hwle, asym, cwave, csun, nc )

     ! Convolve Rayleigh Stuff to the coarse lambda grid

     if (Verbose) write(*,'(4x,A)')'Convolving Rayleigh X-sections and depolarization ratios'
     fray(1:nf,1) = RAYLEIGH_XSEC (1:nf) 
     fray(1:nf,2) = RAYLEIGH_DEPOL(1:nf) 
     CALL asym_gauss_f2c ( fwave, fray, nf, 2, hwle, asym, cwave, Cray, nc )

     !  Solar-I0 convolution of cross sections
     !   - For 5 UV/Vis gases with no temperature-dependence, just do top level and copy rest
     !   - For all other gases, do each level separately (0 through nlayers)
     !  Regular Convolution with no Solar-I0 for the thermal-only cases. 5/13/16

     do g = 1, ngases
        scalex(g) = 1.0d+20 ; IF ( WHICH_GASES(G) == 'O2O2' ) scalex(g) = 1.0d+45
        if (Verbose) write(*,'(4x,A)')'Convolving Trace Gas X-sections (with Solar-I0 effect) for '//which_gases(g)
        IF ( WHICH_GASES(G) == 'O2O2' .or. WHICH_GASES(G) == 'SO2 ' .or. &
             WHICH_GASES(G) == 'HCHO' .or. WHICH_GASES(G) == 'NO2 ' .or. & 
             WHICH_GASES(G) == 'BRO ' ) THEN
           fabs1(1:nf,1) = fsun(1:nf,1) * exp ( - scalex(g) * GAS_XSECS(1:nf,G) ) 
           CALL asym_gauss_f2c ( fwave, fAbs1, nf, 1, hwle, asym, cwave, cabs1, nc )
           CXSec1(1:nc,g) = - Log ( cabs1(1:nc,1)/csun(1:nc,1) ) / scalex(g)
        ELSE IF ( WHICH_GASES(G) == 'O3  ' ) then
           fabs1(1:nf,1) = fsun(1:nf,1) * exp ( - scalex(g) * GAS_XSECS(1:nf,G) ) 
           CALL asym_gauss_f2c ( fwave, fAbs1, nf, 1, hwle, asym, cwave, cabs1, nc )
           CXSec1(1:nc,g) = - Log ( cabs1(1:nc,1)/csun(1:nc,1) ) / scalex(g)
           scaleO3_c1c2 = scalex(g) * 1.0d+05
           fabs1(1:nf,1) = fsun(1:nf,1) * exp ( - scaleO3_c1c2 * O3C1_XSECS(1:nf) ) 
           CALL asym_gauss_f2c ( fwave, fAbs1, nf, 1, hwle, asym, cwave, cabs1, nc )
           CO3C1(1:nc) = - Log ( cabs1(1:nc,1)/csun(1:nc,1) ) / scaleO3_c1c2
           fabs1(1:nf,1) = fsun(1:nf,1) * exp ( - scaleO3_c1c2 * O3C2_XSECS(1:nf) ) 
           CALL asym_gauss_f2c ( fwave, fAbs1, nf, 1, hwle, asym, cwave, cabs1, nc )
           CO3C2(1:nc) = - Log ( cabs1(1:nc,1)/csun(1:nc,1) ) / scaleO3_c1c2
        ELSE IF ( WHICH_GASES(G) == 'H2O ' ) then
           do n = 1, nlevels
              n1 = n - 1 ; fabs(1:nf,n) = fsun(1:nf,1) * exp ( - scalex(g) * H2O_XSECS(1:nf,n1) ) 
           enddo
           CALL asym_gauss_f2c ( fwave, fabs, nf, nlevels, hwle, asym, cwave, cabsM, nc )
           do n = 1, nlevels
              CH2O_M(1:nc,n) = - Log ( cabsM(1:nc,n)/Csun(1:nc,1) ) / scalex(g)
           enddo
        ELSE IF ( WHICH_GASES(G) == 'O2  ' ) then
           do n = 1, nlevels
              n1 = n - 1 ; fabs(1:nf,n) = fsun(1:nf,1) * exp ( - scalex(g) * O2_XSECS(1:nf,n1) ) 
           enddo
           CALL asym_gauss_f2c ( fwave, fabs, nf, nlevels, hwle, asym, cwave, cabsM, nc )
           do n = 1, nlevels
              CO2_M(1:nc,n) = - Log ( cabsM(1:nc,n)/Csun(1:nc,1) ) / scalex(g)
           enddo
        ENDIF
     enddo

!  Re-assignment
!  -------------

!  wavelength grids

     Nlambdas     = nc 
     lambdas      = 0.0_fpk ; lambdas(1:nc)  = clambdas(1:nc)

!  Cross-sections

     GAS_XSECS  = 0.0_fpk ;  GAS_XSECS (1:nc,1:ngases) = CXSec1(1:nc,1:ngases)
     O3C1_XSECS = 0.0_fpk ;  O3C1_XSECS(1:nc) = CO3C1(1:nc)
     O3C2_XSECS = 0.0_fpk ;  O3C2_XSECS(1:nc) = CO3C2(1:nc)
     H2O_XSECS  = 0.0_fpk ;  H2O_XSECS (1:nc,0:nlayers) = CH2O_M(1:nc,1:nlevels)
     O2_XSECS   = 0.0_fpk ;  O2_XSECS  (1:nc,0:nlayers) = CO2_M (1:nc,1:nlevels)

!  Rayleigh

     RAYLEIGH_XSEC   = 0.0_fpk ; RAYLEIGH_XSEC(1:nc)  = CRay (1:nc,1)
     RAYLEIGH_DEPOL  = 0.0_fpk ; RAYLEIGH_DEPOL(1:nc) = CRay (1:nc,2)

!  Re-assign solar spectrum

     SunSpec = 0.0_fpk ; Sunspec(1:nc)  = csun(1:nc,1)

!   deallocate

     deallocate ( csun, Cray, cabs1, cabsM, CXSec1, CO3C1, CO3C2, CH2O_M, CO2_M )
     deallocate ( cwave,asym,hwle )

!  Perform convolution with Empirical Slit function from Instrument RSR
!  ====================================================================

   else

     !Assign instrument empirical slit function. Allocate additional array 8/31/16 for Lambda_eff.

     nr    = Inputs%Instruments%rsr_nlambda ; allocate ( rwave(nr),rsr(nr),frsr(nf), frsrfsun(nf) )
     rwave = Inputs%Instruments%rsr_lambda(1:nr)
     rsr   = Inputs%Instruments%rsr_vals(1:nr)
     nc    = 1

     ! Allocate physical quantities

     allocate ( csun(nc,1), Cray(nc,2), cabs1(nc,1), cabsM(nc,nlevels), CO3C1(nc), CO3C2(nc), &
                CXSec1(nc,ngases), CH2O_M(nc,nlevels), CO2_M(nc,nlevels) )

     !Interpolate RSR to the fine spectral grid

     CALL RSR_to_fine_grid ( rwave, rsr, fwave, frsr )

     !Convolve Solar Spectrum to the coarse lambda grid. Now optional 5/13/16
     ! 8/31/16. Optional calculation of effective wavelength. R. Spurr

     CSun(1,1) = wtd_mean(nf,frsr,fsun(:,1))
     if ( do_effective_wavelength ) then
        do n = 1, nf
           frsrfsun(n) = frsr(n) * fsun(n,1)
        enddo
        clambdas(1) = wtd_mean(nf,frsrfsun,fwave)
     else
        clambdas(1) = Inputs%Instruments%rsr_ctrwl
     endif

!  debug to Fort.90. Response function and Solar Spectrum (fine grid)
!     do n = 1, nf
!        write(90,*)n,frsr(n),1.0d+04*fsun(n,1)
!     enddo

     !Convolve Rayleigh Stuff to the coarse lambda grid
     if (Verbose) write(*,'(4x,A)')'Convolving rayleigh X-sections and depolarization ratios'
     fray(1:nf,1) = RAYLEIGH_XSEC (1:nf) ; CRay(1,1) = wtd_mean(nf,frsr,fray(:,1))
     fray(1:nf,2) = RAYLEIGH_DEPOL(1:nf) ; CRay(1,2) = wtd_mean(nf,frsr,fray(:,2))

     !  Solar-I0 convolution of cross sections
     !   - For 5 UV/Vis gases with no temperature-dependence, just do top level and copy rest
     !   - For all other gases, do each level separately (0 through nlayers)
     !  Regular Convolution with no Solar-I0 for the thermal-only cases. 5/13/16

     do g = 1, ngases
        scalex(g) = 1.0d+20 ; IF ( WHICH_GASES(G) == 'O2O2' ) scalex(g) = 1.0d+45
        if (Verbose) write(*,'(4x,A)')'Convolving Trace Gas X-sections (with Solar-I0 effect) for '//which_gases(g)

        IF ( WHICH_GASES(G) == 'O2O2' .or. WHICH_GASES(G) == 'SO2 ' .or. &
             WHICH_GASES(G) == 'HCHO' .or. WHICH_GASES(G) == 'NO2 ' .or. & 
             WHICH_GASES(G) == 'BRO ' ) THEN
           fabs1(1:nf,1) = fsun(1:nf,1) * exp ( - scalex(g) * GAS_XSECS(1:nf,G) ) 
           Cabs1(1,1) = wtd_mean(nf,frsr,fabs1(:,1))                    
           CXSec1(1:nc,g) = - Log ( cabs1(1:nc,1)/csun(1:nc,1) ) / scalex(g)
        ELSE IF ( WHICH_GASES(G) == 'O3  ' ) then
           fabs1(1:nf,1) = fsun(1:nf,1) * exp ( - scalex(g) * GAS_XSECS(1:nf,G) ) 
           Cabs1(1,1) = wtd_mean(nf,frsr,fabs1(:,1))                    
           CXSec1(1:nc,g) = - Log ( cabs1(1:nc,1)/csun(1:nc,1) ) / scalex(g)
           scaleO3_c1c2 = scalex(g) * 1.0d+05
           fabs1(1:nf,1) = fsun(1:nf,1) * exp ( - scaleO3_c1c2 * O3C1_XSECS(1:nf) ) 
           Cabs1(1,1) = wtd_mean(nf,frsr,fabs1(:,1)) 
           CO3C1(1) = - Log ( cabs1(1,1)/csun(1,1) ) / scaleO3_c1c2
           fabs1(1:nf,1) = fsun(1:nf,1) * exp ( - scaleO3_c1c2 * O3C2_XSECS(1:nf) ) 
           Cabs1(1,1) = wtd_mean(nf,frsr,fabs1(:,1))                       
           CO3C2(1:1) = - Log ( cabs1(1,1)/csun(1,1) ) / scaleO3_c1c2
        ELSE IF ( WHICH_GASES(G) == 'H2O ' ) then
           do n = 1, nlevels
              n1 = n - 1 ; fabs(1:nf,n) = fsun(1:nf,1) * exp ( - scalex(g) * H2O_XSECS(1:nf,n1) ) 
              cabsM(1,n) = wtd_mean(nf,frsr,fabs(1:nf,n))                         
              CH2O_M(1,n) = - Log ( cabsM(1,n)/Csun(1,1) ) / scalex(g)
           enddo
        ELSE IF ( WHICH_GASES(G) == 'O2  ' ) then
           do n = 1, nlevels
              n1 = n - 1 ; fabs(1:nf,n) = fsun(1:nf,1) * exp ( - scalex(g) * O2_XSECS(1:nf,n1) ) 
              cabsM(1,n) = wtd_mean(nf,frsr,fabs(1:nf,n))                         
              CO2_M(1,n) = - Log ( cabsM(1,n)/Csun(1,1) ) / scalex(g)
           enddo
        ENDIF
     enddo

!  Re-assignment
!  -------------

!  wavelength grids

     Nlambdas     = 1 
     lambdas      = 0.0_fpk ; lambdas(1)  = clambdas(1)

!  Cross-sections

     GAS_XSECS  = 0.0_fpk ;  GAS_XSECS (1:nc,1:ngases) = CXSec1(1:nc,1:ngases)
     O3C1_XSECS = 0.0_fpk ;  O3C1_XSECS(1:nc) = CO3C1(1:nc)
     O3C2_XSECS = 0.0_fpk ;  O3C2_XSECS(1:nc) = CO3C2(1:nc)
     H2O_XSECS  = 0.0_fpk ;  H2O_XSECS (1:nc,0:nlayers) = CH2O_M(1:nc,1:nlevels)
     O2_XSECS   = 0.0_fpk ;  O2_XSECS  (1:nc,0:nlayers) = CO2_M (1:nc,1:nlevels)

!  Rayleigh

     RAYLEIGH_XSEC   = 0.0_fpk ; RAYLEIGH_XSEC(1:nc)  = CRay (1:nc,1)
     RAYLEIGH_DEPOL  = 0.0_fpk ; RAYLEIGH_DEPOL(1:nc) = CRay (1:nc,2)

!  Re-assign solar spectrum

     SunSpec = 0.0_fpk ; Sunspec(1:nc)  = csun(1:nc,1)

!   deallocate

     deallocate ( csun, Cray, cabs1, cabsM, CXSec1, CO3C1, CO3C2, CH2O_M, CO2_M )
     deallocate ( rwave,rsr,frsr,frsrfsun )

   end if

!  Deallocate fine grid stuff

   deallocate ( fwave, fsun, fRay, fabs1, fabs )

!  End subroutine

end subroutine GEMSTOOL_XSECS_CONVOLVE

!  End module

end module GEMSTOOL_UVN_Convolution_m
