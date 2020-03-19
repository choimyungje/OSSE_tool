module cross_sections_m

!  Revision. 03 September 2012
!   - GET_O2O2XSEC_1, introduced for O2-O2 absorption (Greenblatt 1990 data)

!  Revision. 18 February 2013
!   - Added several more Cross-sections, from GEOCAPE (SAO) compilation.
!   - separated Cross-section data files in Physics_data

!  Type definitions

   use vlidort_pars, only : fpk

!  Use numerical and dataread_101 Modules

   use dataread_101_m, only : read_coe
   use numerical_m   , only : spline, splint, linear_trawler_ascending

!  Only main routine is public

private :: GET_O3XSEC_1,&
           GET_O3XSEC_2,&
           GET_O3XSEC_3,&
           GET_NO2XSEC_1,  &       ! New 02/18/13
           GET_BROXSEC_1,  &       ! New 02/18/13
           GET_HCHOXSEC_1, &       ! New 02/18/13
           GET_SO2XSEC_1,  &
           GET_SO2XSEC_2,  &       ! New 02/18/13
           GET_O2O2XSEC_1, &       ! New 09/03/12
           GET_O2O2XSEC_2, &       ! New 02/18/13
           RAYLEIGH_FUNCTION

public  :: CROSS_SECTIONS

contains

SUBROUTINE CROSS_SECTIONS                                       &
          ( MAXLAMBDAS, MAXLAYERS, MAXGASES, CO2_PPMV_MIXRATIO, & ! INPUT
            NLAYERS, NLAMBDAS, TEMPERATURES, LAMBDAS,           & ! INPUT
            WHICH_GASES, XSEC_SOURCES, NGASES,                  & ! INPUT
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL, GAS_XSECS,           & ! OUTPUT
            FAIL, MESSAGE )                                       ! OUTPUT

!  Key for Sources
!  ===============

!  List updated, 18 February 2013
!  ------------------------------

!  [SAO compilation = GEOCAPE data set, prepared by K. Chance for GEOCAPE]

!     for O3   : XSEC_SOURCES(G) = 1 --> GOME Flight model
!     for O3   : XSEC_SOURCES(G) = 2 --> coedat    (Daumont)
!     for O3   : XSEC_SOURCES(G) = 3 --> Unrestricted Daumont/Malicet  (SAO)

!     for NO2  : XSEC_SOURCES(G) = 1 --> SAO Compilation

!     for BRO  : XSEC_SOURCES(G) = 1 --> SAO Compilation. BIRA

!     for HCHO : XSEC_SOURCES(G) = 1 --> SAO Compilation. Cantrell

!     for SO2  : XSEC_SOURCES(G) = 1 --> SCIAMACHY (Bogumil)
!     for SO2  : XSEC_SOURCES(G) = 2 --> SAO Compilation. BIRA

!     for O2O2 : XSEC_SOURCES(G) = 1 --> Greenblatt et al., (1990) New 09/03/12
!     for O2O2 : XSEC_SOURCES(G) = 2 --> SAO Compilation. BIRA

!  Input arguments
!  ---------------

!  Dimensioning

      INTEGER     :: MAXLAMBDAS, MAXLAYERS, MAXGASES

!  wavelengths
 
      INTEGER     :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  CO2 mixing ratio

      real(fpk)    :: CO2_PPMV_MIXRATIO

!  Layer control

      INTEGER     :: NLAYERS
      real(fpk),    dimension ( MAXLAYERS ) :: TEMPERATURES 
      
!  Gas control

      integer                                  :: NGASES
      character(LEN=4), dimension ( maxgases ) :: WHICH_GASES
      integer,          dimension ( maxgases ) :: XSEC_SOURCES
      
!  Output arguments
!  ----------------

!  Rayleigh cross-sections and depolarization output

      real(fpk),    dimension ( MAXLAMBDAS ) :: RAYLEIGH_XSEC 
      real(fpk),    dimension ( MAXLAMBDAS ) :: RAYLEIGH_DEPOL

!  Gas cross-sections

      REAL(fpk),    DIMENSION( MAXLAYERS, MAXLAMBDAS, maxgases ) :: GAS_XSECS  

!  status

      LOGICAL       ::        FAIL
      CHARACTER*(*) ::        MESSAGE

!  Local variables
!  ---------------

      integer           ::  G
      CHARACTER(LEN=60) ::  FILENAME

!  Gases

   DO G = 1, NGASES

!  O3 (Ozone)
!  ==========

      IF ( WHICH_GASES(G) == 'O3  ' ) THEN

!  ... GOME FM

         IF ( XSEC_SOURCES(G) == 1 ) then
            FILENAME = 'physics_data/Xsections/gome_o3xsecfm_305365.pre'
            CALL GET_O3XSEC_1                         &  
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          &
            NLAMBDAS, NLAYERS, LAMBDAS, temperatures, &
            gas_xsecs(:,:,G), FAIL, message )

!  ... COE.dat (Fixed wavelengths, no interpolation)

         ELSE IF ( XSEC_SOURCES(G) == 2 ) then 
            CALL GET_O3XSEC_2                         &
          ( MAXLAMBDAS, MAXLAYERS,                    & 
            NLAMBDAS, NLAYERS, LAMBDAS, temperatures, &
            gas_xsecs(:,:,G) )

!  ... Brion/Daumont

         ELSE IF ( XSEC_SOURCES(G) == 3 ) then
            FILENAME ='physics_data/Xsections/o3abs_brion_195_660_vacfinal.dat'
            CALL GET_O3XSEC_3                         &
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          & 
            NLAMBDAS, NLAYERS, LAMBDAS, temperatures, &
            gas_xsecs(:,:,G), FAIL, message )
         ENDIF

!  SO2 (Sulfur Dioxide)
!  ====================

      ELSE IF ( WHICH_GASES(G) == 'SO2 ' ) THEN

!  ... SCIAMACHY flight model (1024 points). Bogumil

         IF ( XSEC_SOURCES(G) == 1 ) then
            FILENAME = 'physics_data/Xsections/so2_.xs'
            CALL GET_SO2XSEC_1                         &
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          & 
            NLAMBDAS, NLAYERS, LAMBDAS,               &
            gas_xsecs(:,:,G), FAIL, message )

!  ... Source unknown BIRA/IASB High Resolution. Bogumil SCIAMACHY ???

         ELSE IF ( XSEC_SOURCES(G) == 2 ) then
            FILENAME = 'physics_data/Xsections/SO2_298_BISA.nm'
            CALL GET_SO2XSEC_2                         &
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          & 
            NLAMBDAS, NLAYERS, LAMBDAS,               &
            gas_xsecs(:,:,G), FAIL, message )
         ENDIF

!  NO2 (Nitrogen Dioxide)
!  ======================

      ELSE IF ( WHICH_GASES(G) == 'NO2 ' ) THEN

!  ... no2r_97.nm

         IF ( XSEC_SOURCES(G) == 1 ) then
            FILENAME = 'physics_data/Xsections/no2r_97.nm'
            CALL GET_NO2XSEC_1                         &
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          & 
            NLAMBDAS, NLAYERS, LAMBDAS,               &
            gas_xsecs(:,:,G), FAIL, message )
         ENDIF

!  BRO (Bromine Monoxide)
!  ======================

      ELSE IF ( WHICH_GASES(G) == 'BRO ' ) THEN

!  ... 228kbro10cm.nm

         IF ( XSEC_SOURCES(G) == 1 ) then
            FILENAME = 'physics_data/Xsections/228kbro10cm.nm'
            CALL GET_BROXSEC_1                         &
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          & 
            NLAMBDAS, NLAYERS, LAMBDAS,               &
            gas_xsecs(:,:,G), FAIL, message )
         ENDIF

!  HCHO (Formaldehyde)
!  ===================

      ELSE IF ( WHICH_GASES(G) == 'HCHO' ) THEN

!  ... h2co_cc.300 (Chris Cantrell, SAO)

         IF ( XSEC_SOURCES(G) == 1 ) then
            FILENAME = 'physics_data/Xsections/h2co_cc.300'
            CALL GET_HCHOXSEC_1                        &
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          & 
            NLAMBDAS, NLAYERS, LAMBDAS,               &
            gas_xsecs(:,:,G), FAIL, message )
         ENDIF

!  O2O2 (Oxygen Dimer)
!  ===================

      ELSE IF ( WHICH_GASES(G) == 'O2O2' ) THEN

!  ... Greenblatt 1990 data

         IF ( XSEC_SOURCES(G) == 1 ) then
            FILENAME = 'physics_data/Xsections/Greenblatt_1990_O2O2.xs'
            CALL GET_O2O2XSEC_1                       &
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          & 
            NLAMBDAS, NLAYERS, LAMBDAS,               &
            gas_xsecs(:,:,G), FAIL, message )

!  ... O4_294K_BISA.dat, high resolution

         ELSE IF ( XSEC_SOURCES(G) == 2 ) then
            FILENAME = 'physics_data/Xsections/O4_294K_BISA.dat'
            CALL GET_O2O2XSEC_2                       &
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          & 
            NLAMBDAS, NLAYERS, LAMBDAS,               &
            gas_xsecs(:,:,G), FAIL, message )
         ENDIF

      ENDIF
      IF ( FAIL ) RETURN
   ENDDO

!  Get the Rayleigh. Input wavelengths in [nm], NOT ANGSTROMS

   CALL RAYLEIGH_FUNCTION                                &
     ( MAXLAMBDAS, CO2_PPMV_MIXRATIO, NLAMBDAS, LAMBDAS, &
       RAYLEIGH_XSEC, RAYLEIGH_DEPOL )

!  FInish

  RETURN
END SUBROUTINE CROSS_SECTIONS

!

SUBROUTINE RAYLEIGH_FUNCTION                       &
          ( FORWARD_MAXLAMBDAS, CO2_PPMV_MIXRATIO, &
            FORWARD_NLAMBDAS,   FORWARD_LAMBDAS,   &
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL )

!  Rayleigh cross sections and depolarization ratios
!     Bodhaine et. al. (1999) formulae
!     Module is stand-alone.
!     Wavelengths in nm (formerly Angstroms)

!  Input arguments
!  ---------------

!  wavelength
 
      INTEGER     :: FORWARD_MAXLAMBDAS, FORWARD_NLAMBDAS
      real(fpk),    dimension ( FORWARD_MAXLAMBDAS ) :: FORWARD_LAMBDAS 

!  CO2 mixing ratio

      real(fpk)    :: CO2_PPMV_MIXRATIO

!  Output arguments
!  ----------------

!  cross-sections and depolarization output

      real(fpk),    dimension ( FORWARD_MAXLAMBDAS ) :: RAYLEIGH_XSEC 
      real(fpk),    dimension ( FORWARD_MAXLAMBDAS ) :: RAYLEIGH_DEPOL

!  Local variables
!  ---------------

      INTEGER      :: W
      real(fpk)    :: MASS_DRYAIR
      real(fpk)    :: NMOL, PI, CONS
      real(fpk)    :: MO2,MN2,MARG,MCO2,MAIR
      real(fpk)    :: FO2,FN2,FARG,FCO2,FAIR
      real(fpk)    :: LAMBDA_C,LAMBDA_M,LPM2,LP2
      real(fpk)    :: N300M1,NCO2M1,NCO2
      real(fpk)    :: NCO2SQ, NSQM1,NSQP2,TERM

!  data statements and parameters
!  ------------------------------

      DATA MO2  / 20.946D0 /
      DATA MN2  / 78.084D0 /
      DATA MARG / 0.934D0 /

      real(fpk),    PARAMETER ::        S0_A = 15.0556D0
      real(fpk),    PARAMETER ::        S0_B = 28.9595D0

      real(fpk),    PARAMETER ::        S1_A = 8060.51D0
      real(fpk),    PARAMETER ::        S1_B = 2.48099D+06
      real(fpk),    PARAMETER ::        S1_C = 132.274D0
      real(fpk),    PARAMETER ::        S1_D = 1.74557D+04
      real(fpk),    PARAMETER ::        S1_E = 39.32957D0

      real(fpk),    PARAMETER ::        S2_A = 0.54D0

      real(fpk),    PARAMETER ::        S3_A = 1.034D0
      real(fpk),    PARAMETER ::        S3_B = 3.17D-04
      real(fpk),    PARAMETER ::        S3_C = 1.096D0
      real(fpk),    PARAMETER ::        S3_D = 1.385D-03
      real(fpk),    PARAMETER ::        S3_E = 1.448D-04

!  Start of code
!  -------------

!  constants

      NMOL = 2.546899D19
      PI   = DATAN(1.0D0)*4.0D0
      CONS = 24.0D0 * PI * PI * PI

!  convert co2

      MCO2 = 1.0D-06 * CO2_PPMV_MIXRATIO

!  mass of dry air: Eq.(17) of BWDS

      MASS_DRYAIR = S0_A * MCO2 + S0_B

!  start loop

      DO W = 1, FORWARD_NLAMBDAS

!  wavelength in micrometers

      LAMBDA_M = 1.0D-03 * FORWARD_LAMBDAS(W)
      LAMBDA_C = 1.0D-07 * FORWARD_LAMBDAS(W)
!      LAMBDA_M = 1.0D-04 * FORWARD_LAMBDAS(W)             ! Angstrom input
!      LAMBDA_C = 1.0D-08 * FORWARD_LAMBDAS(W)             ! Angstrom input
      LPM2     = 1.0D0 / LAMBDA_M / LAMBDA_M

!  step 1: Eq.(18) of BWDS

      N300M1 = S1_A + ( S1_B / ( S1_C - LPM2 ) ) + &
                      ( S1_D / ( S1_E - LPM2 ) )
      N300M1 = N300M1 * 1.0D-08

!  step 2: Eq.(19) of BWDS

      NCO2M1 = N300M1 * ( 1.0D0 + S2_A * ( MCO2  - 0.0003D0 ) )
      NCO2   = NCO2M1 + 1
      NCO2SQ = NCO2 * NCO2

!  step 3: Eqs. (5&6) of BWDS (Bates' results)

      FN2  = S3_A + S3_B * LPM2
      FO2  = S3_C + S3_D * LPM2 + S3_E * LPM2 * LPM2

!  step 4: Eq.(23) of BWDS
!     ---> King factor and depolarization ratio

      FARG = 1.0D0
      FCO2 = 1.15D0
      MAIR = MN2 + MO2 + MARG + MCO2
      FAIR = MN2*FN2 + MO2*FO2 + MARG*FARG + MCO2*FCO2
      FAIR = FAIR / MAIR
      RAYLEIGH_DEPOL(W) = 6.0D0*(FAIR-1.0D0)/(3.0D0+7.0D0*FAIR)

!  step 5: Eq.(22) of BWDS
!     ---> Cross section

      LP2  = LAMBDA_C * LAMBDA_C
      NSQM1 = NCO2SQ - 1.0D0
      NSQP2 = NCO2SQ + 2.0D0
      TERM = NSQM1 / LP2 / NMOL / NSQP2
      RAYLEIGH_XSEC(W) =  CONS * TERM * TERM * FAIR

!  end loop

      ENDDO

!  finish
!  ------

      RETURN
END SUBROUTINE RAYLEIGH_FUNCTION 



SUBROUTINE GET_O3XSEC_1                               &
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          & 
            NLAMBDAS, NLAYERS, LAMBDAS, temperatures, &
            o3_crosssecs, FAIL, message )

!  file detail

      CHARACTER*(*) ::   FILENAME

!  Dimensioning

      INTEGER       :: MAXLAMBDAS, MAXLAYERS

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Layer control

      INTEGER                                :: NLAYERS
      real(fpk),    dimension ( MAXLAYERS )  :: TEMPERATURES 
       
!  Output arguments
!  ----------------

!  Gas cross-sections

      REAL(fpk),    DIMENSION( MAXLAYERS, MAXLAMBDAS ) :: O3_CROSSSECS 

!  status

      LOGICAL       ::        FAIL
      CHARACTER*(*) ::        MESSAGE

!  Local variables
!  ===============

!  Dimensioning

    INTEGER, parameter ::  O3XS_MAX_DATAPOINTS = 1000
    INTEGER, parameter ::  O3XS_MAX_TEMPS      = 5   

!  Cross section data
!  ------------------

    INTEGER   ::       O3XS_N_TEMPS
    real(fpk),   dimension (O3XS_MAX_DATAPOINTS) :: O3XS_WAVEDATA 
    real(fpk),   dimension (O3XS_MAX_DATAPOINTS) :: O3XS_XSECDATA
    real(fpk),   dimension (O3XS_MAX_DATAPOINTS) :: O3XS_XSEC2
    real(fpk),   dimension (O3XS_MAX_TEMPS)      :: O3XS_TEMPDATA 

!  help variables

    INTEGER      :: W, N, K, K1, K2,M, NP
    real(fpk)    :: DF1,DF2

!  variables for splined interpolation

    real(fpk),   dimension(O3XS_MAX_DATAPOINTS,O3XS_MAX_TEMPS) :: XSEC

!  initialise file read

      FAIL = .FALSE.

!  Open file

      OPEN ( 55, FILE = FILENAME, ERR = 900, STATUS = 'OLD' )

!  read and check temperature data first

      READ(55,*)O3XS_N_TEMPS
      IF ( O3XS_N_TEMPS .GT. O3XS_MAX_TEMPS ) GOTO 902
      DO K = 1, O3XS_N_TEMPS
        READ(55,*)O3XS_TEMPDATA(K)
      ENDDO

!  Start Data loop

      DO K = 1, O3XS_N_TEMPS

!  main read section for basic data file

        READ(55,*)NP
        IF ( NP .GT. O3XS_MAX_DATAPOINTS ) GOTO 902
        DO M = 1, NP
          READ(55,*) O3XS_WAVEDATA(M),O3XS_XSECDATA(M)
        ENDDO

!  Spline interpolate immediately

        CALL SPLINE ( O3XS_WAVEDATA,O3XS_XSECDATA,     &
                       NP, 0.0d0, 0.0d0,  O3XS_XSEC2 )
          
        DO W = 1, NLAMBDAS
          IF ( LAMBDAS(W) .le. O3XS_WAVEDATA(NP) .AND. &
               LAMBDAS(W) .ge. O3XS_WAVEDATA(1) ) THEN 
            CALL SPLINT ( O3XS_WAVEDATA,O3XS_XSECDATA, &
                          O3XS_XSEC2,NP,LAMBDAS(W),XSEC(W,K))
          ELSE
            XSEC(W,K) = 0.0d0
          ENDIF
        ENDDO

!  Finish data loop

        ENDDO

!  CLose file

      CLOSE ( 55 )

!  debug
!      do w = 1, nlambdas
!        write(*,'(f8.2,1p5e15.5)')lambdas(w),(xsec(w,k1),k1 = 1,5)
!      enddo
!      pause

!  Set the cross-sections
!  ======================

!  temperature interpolate (linear) clear sky cross sections

      DO N = 1, NLAYERS
        CALL LINEAR_TRAWLER_ASCENDING         &
           ( O3XS_MAX_TEMPS, O3XS_N_TEMPS,    &
             O3XS_TEMPDATA,  temperatures(n), &
             DF1, DF2, K1, K2 )
        DO W = 1, NLAMBDAS
          O3_CROSSSECS(N,W) = DF1 * XSEC(W,K2) + DF2 * XSEC(W,K1)
        ENDDO
      ENDDO

!  Finish

      RETURN

!  Open-file error return

900   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Could not find file of O3 cross-sections to open'
      CLOSE(55)
      RETURN

!  dimensioning error return

902   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Dimensioning insufficient for O3 cross-section data'
      CLOSE(55)
      RETURN

END SUBROUTINE GET_O3XSEC_1




SUBROUTINE GET_O3XSEC_2 ( MAXLAMBDAS, MAXLAYERS,      & 
            NLAMBDAS, NLAYERS, LAMBDAS, temperatures, &
            o3_crosssecs )

!  Dimensioning

      INTEGER       :: MAXLAMBDAS, MAXLAYERS

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Layer control

      INTEGER                                :: NLAYERS
      real(fpk),    dimension ( MAXLAYERS )  :: TEMPERATURES 
       
!  Output arguments
!  ----------------

!  Gas cross-sections

      REAL(fpk),    DIMENSION( MAXLAYERS, MAXLAMBDAS ) :: O3_CROSSSECS 

!  Local stuff
!  -----------

!   real(fpk),    dimension(500,5) :: coedat
   real(fpk),    dimension(:,:), allocatable :: coedat
   real(fpk)    :: t1, t2, sig2
   integer      :: w, i

!  Read the basic data

!   allocate(coedat(maxlambdas,5))
!   call read_coe('physics_data/coe.dat', maxlambdas, lambdas, nlambdas, coedat)
    call read_coe('physics_data/coe.dat', 500, lambdas, nlambdas, coedat)

! Layer cross-sections for ozone

  do w = 1, nlambdas
     do i = 1, nlayers
        t1 = temperatures(i) - 273.15
        t2 =  t1 * t1
        sig2 = coedat(w, 1) + t1 * coedat(w,2) + t2 * coedat(w,3)
        o3_crosssecs(i,w) = sig2
     enddo
  enddo

!  Finish

  return
END SUBROUTINE GET_O3XSEC_2

!
SUBROUTINE GET_O3XSEC_3 ( FILENAME, MAXLAMBDAS, MAXLAYERS,         & 
                          NLAMBDAS, NLAYERS, LAMBDAS, temperatures, &
                          o3_crosssecs, fail, message )

!  Dimensioning

      INTEGER       :: MAXLAMBDAS, MAXLAYERS

!  Filename

      character*(*) :: FILENAME

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Layer control

      INTEGER                                :: NLAYERS
      real(fpk),    dimension ( MAXLAYERS )  :: TEMPERATURES 
       
!  Output arguments
!  ----------------

!  Gas cross-sections

      REAL(fpk),    DIMENSION( MAXLAYERS, MAXLAMBDAS ) :: O3_CROSSSECS 

!  Exceptiopn handling

      logical             :: FAIL
      character*(*)       :: MESSAGE

!  Local stuff
!  -----------

   real(fpk)    :: t1, t2, sig2

   integer, parameter :: maxspline = 1000

   logical       :: reading
   integer       :: w, i, nbuff, ndata, n, maxdata
   real(fpk)     :: lamb1, lamb2, wav, val, c1, c2, conversion, xsec
   real(fpk)     :: x(maxspline), y(maxspline), y2(maxspline)
   real(fpk)     :: yc1(maxspline), yc2(maxspline)

   REAL(fpk),    DIMENSION( MAXLAMBDAS ) :: o3c0_xsecs 
   REAL(fpk),    DIMENSION( MAXLAMBDAS ) :: o3c1_xsecs 
   REAL(fpk),    DIMENSION( MAXLAMBDAS ) :: o3c2_xsecs 

!  initialize exception handling

   fail    = .false.
   message = ' '

!  Initialize output

   o3c0_xsecs = 0.0d0;  o3c1_xsecs = 0.0d0;  o3c2_xsecs = 0.0d0

!  Buffer control, initialization

   nbuff = 0
   ndata = 1
   lamb1 = lambdas(1)        - 0.5d0
   lamb2 = lambdas(nlambdas) + 0.5d0
   maxdata    = 46522
   conversion = 1.0d+20

!  Read the data, check that window is covered

   open(1,file=adjustl(trim(filename)),err=90,status='old')

!  .. First line is a dummy
         read(1,*)
!  .. Read second line, and check window is covered
         read(1,*)wav, val, c1, c2
         if (lamb1 .lt. wav ) then
            fail = .true.
            message = 'O3 Xsec data does not cover input window at lower end'
            return
         endif
!  .. Read more lines, saving to spline buffer
         reading = .true.
         do while (reading .and. ndata .lt. maxdata )
            ndata = ndata + 1
            read(1,*) wav, val, c1, c2
            if ( wav .gt. lamb1 ) then
               if ( wav .lt. lamb2 ) then
                  nbuff = nbuff + 1
                  x(nbuff)   = wav
                  y(nbuff)   = val
                  yc1(nbuff) = c1
                  yc2(nbuff) = c2
               endif
               if ( wav .ge. lamb2 ) reading = .false.
            endif
         enddo
!  .. Check if last data line is present, then window not covered
         if (ndata .eq. maxdata) then
            fail = .true.
            message = 'O3 Xsec data does not cover input window at upper end'
            return
         endif
         close(1)

         if (nbuff .gt. maxspline) then
            fail = .true.
            message = 'need to increase maxspline in GET_O3XSEC_3'
            return
         endif
         
!        write(*,*)nbuff, lamb1, lamb2
!         do n = 1, nbuff
!            write(*,*)x(n),y(n),yc1(n)
!         enddo
!         pause
         
!  Spline the data
         call spline(x,y,nbuff,0.0d0,0.0d0,y2)

         do n = 1, nlambdas
           call splint(x,y,y2,nbuff,lambdas(n),xsec)
           o3c0_xsecs(n) = xsec / conversion
         enddo

         call spline(x,yc1,nbuff,0.0d0,0.0d0,y2)
         do n = 1, nlambdas
           call splint(x,yc1,y2,nbuff,lambdas(n),xsec)
           o3c1_xsecs(n) = xsec / conversion
         enddo

         call spline(x,yc2,nbuff,0.0d0,0.0d0,y2)
         do n = 1, nlambdas
           call splint(x,yc2,y2,nbuff,lambdas(n),xsec)
           o3c2_xsecs(n) = xsec / conversion
         enddo

! Layer cross-sections for ozone

  do w = 1, nlambdas
     do i = 1, nlayers
        t1 = temperatures(i) - 273.15
        t2 =  t1 * t1
        sig2 = o3c0_xsecs(w) + t1 * o3c1_xsecs(w) + t2 * o3c2_xsecs(w)
        o3_crosssecs(i,w) = sig2
     enddo
  enddo

!  normal return

   return

!  error return

90 continue
   fail    = .true.
   message = 'Open failure for Xsec file = '//filename(1:LEN(filename))
   return

!  Finish

  return
END SUBROUTINE GET_O3XSEC_3

!

SUBROUTINE GET_NO2XSEC_1 ( FILENAME, MAXLAMBDAS, MAXLAYERS,  & 
                           NLAMBDAS, NLAYERS, LAMBDAS,       &
                           no2_crosssecs, fail, message )

!  Dimensioning

      INTEGER       :: MAXLAMBDAS, MAXLAYERS

!  Filename

      character*(*) :: FILENAME

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Layer control

      INTEGER                                :: NLAYERS
       
!  Output arguments
!  ----------------

!  Gas cross-sections

      REAL(fpk),    DIMENSION( MAXLAYERS, MAXLAMBDAS ) :: NO2_CROSSSECS 

!  Exceptiopn handling

   logical             :: FAIL
   character*(*)       :: MESSAGE

!  Local stuff
!  -----------

   integer, parameter :: maxspline = 1000
   logical       :: reading
   integer       :: i, nbuff, ndata, n, maxdata
   real(fpk)     :: lamb1, lamb2, wav, val, conversion, xsec
   real(fpk)     :: x(maxspline), y(maxspline), y2(maxspline)

!  initialize exception handling

   fail    = .false.
   message = ' '

!  Buffer control, initialization

   nbuff = 0
   ndata = 1
   lamb1 = lambdas(1)        - 0.5d0
   lamb2 = lambdas(nlambdas) + 0.5d0
   maxdata    = 21485
   conversion = 1.0d+20

!  Read the data, check that window is covered

    open(1,file=adjustl(trim(filename)),err=90,status='old')

!  .. Read First line, and check window is covered
         read(1,*)wav, val
         if (lamb1 .lt. wav ) then
            fail = .true.
            message = 'NO2 Xsec data does not cover input window at lower end'
            return
         endif
!  .. Read more lines, saving to spline buffer
         reading = .true.
         do while (reading .and. ndata .lt. maxdata )
            ndata = ndata + 1
            read(1,*) wav, val
            if ( wav .gt. lamb1 ) then
               if ( wav .lt. lamb2 ) then
                  nbuff = nbuff + 1
                  x(nbuff)   = wav
                  y(nbuff)   = val
               endif
               if ( wav .ge. lamb2 ) reading = .false.
            endif
         enddo
!  .. Check if last data line is present, then window not covered
         if (ndata .eq. maxdata) then
            fail = .true.
            message = 'NO2 Xsec data does not cover input window at upper end'
            return
         endif
         close(1)

         if (nbuff .gt. maxspline) then
            fail = .true.
            message = 'need to increase maxspline in GET_NO2XSEC_1'
            return
         endif
         
!        write(*,*)nbuff, lamb1, lamb2
!         do n = 1, nbuff
!            write(*,*)x(n),y(n),yc1(n)
!         enddo
!         pause
         
!  Spline the data

   call spline(x,y,nbuff,0.0d0,0.0d0,y2)
   do n = 1, nlambdas
      call splint(x,y,y2,nbuff,lambdas(n),xsec)
      no2_crosssecs(1,n) = xsec / conversion
      do i = 2, nlayers
         no2_crosssecs(i,n) = no2_crosssecs(1,n)
      enddo
   enddo

!  normal return

   return

!  error return

90 continue
   fail    = .true.
   message = 'Open failure for Xsec file = '//filename(1:LEN(filename))
   return

!  Finish

   return
END SUBROUTINE GET_NO2XSEC_1

!

SUBROUTINE GET_HCHOXSEC_1 ( FILENAME, MAXLAMBDAS, MAXLAYERS,  & 
                            NLAMBDAS, NLAYERS, LAMBDAS,       &
                            hcho_crosssecs, fail, message )

!  Dimensioning

      INTEGER       :: MAXLAMBDAS, MAXLAYERS

!  Filename

      character*(*) :: FILENAME

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Layer control

      INTEGER                                :: NLAYERS
       
!  Output arguments
!  ----------------

!  Gas cross-sections

      REAL(fpk),    DIMENSION( MAXLAYERS, MAXLAMBDAS ) :: HCHO_CROSSSECS 

!  Exceptiopn handling

   logical             :: FAIL
   character*(*)       :: MESSAGE

!  Local stuff
!  -----------

   integer, parameter :: maxspline = 1000
   logical       :: reading
   integer       :: i, nbuff, ndata, n, maxdata
   real(fpk)     :: lamb1, lamb2, wav, val, conversion, xsec
   real(fpk)     :: x(maxspline), y(maxspline), y2(maxspline)

!  initialize exception handling

   fail    = .false.
   message = ' '

!  Buffer control, initialization

   nbuff = 0
   ndata = 1
   lamb1 = lambdas(1)        - 0.05d0
   lamb2 = lambdas(nlambdas) + 0.05d0
   maxdata    = 15021
   conversion = 1.0d+20

!  Read the data, check that window is covered

    open(1,file=adjustl(trim(filename)),err=90,status='old')

!  .. Read First line, and check window is covered
         read(1,*)wav, val
         if (lamb1 .lt. wav ) then
            fail = .true.
            message = 'HCHO Xsec data does not cover input window at lower end'
            return
         endif
!  .. Read more lines, saving to spline buffer
         reading = .true.
         do while (reading .and. ndata .lt. maxdata )
            ndata = ndata + 1
            read(1,*) wav, val
            if ( wav .gt. lamb1 ) then
               if ( wav .lt. lamb2 ) then
                  nbuff = nbuff + 1
                  x(nbuff)   = wav
                  y(nbuff)   = val
               endif
               if ( wav .ge. lamb2 ) reading = .false.
            endif
         enddo
!  .. Check if last data line is present, then window not covered
         if (ndata .eq. maxdata) then
            fail = .true.
            message = 'HCHO Xsec data does not cover input window at upper end'
            return
         endif
         close(1)

         if (nbuff .gt. maxspline) then
            fail = .true.
            message = 'need to increase maxspline in GET_HCHOXSEC_1'
            return
         endif
         
!        write(*,*)nbuff, lamb1, lamb2
!         do n = 1, nbuff
!            write(*,*)x(n),y(n),yc1(n)
!         enddo
!         pause
         
!  Spline the data

   call spline(x,y,nbuff,0.0d0,0.0d0,y2)
   do n = 1, nlambdas
      call splint(x,y,y2,nbuff,lambdas(n),xsec)
      hcho_crosssecs(1,n) = xsec / conversion
      do i = 2, nlayers
         hcho_crosssecs(i,n) = hcho_crosssecs(1,n)
      enddo
   enddo

!  normal return

   return

!  error return

90 continue
   fail    = .true.
   message = 'Open failure for Xsec file = '//filename(1:LEN(filename))
   return

!  Finish

   return
END SUBROUTINE GET_HCHOXSEC_1

!

SUBROUTINE GET_SO2XSEC_1                              &
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          & 
            NLAMBDAS, NLAYERS, LAMBDAS,               &
            so2_crosssecs, FAIL, message )

!  file detail

      CHARACTER*(*) ::   FILENAME

!  Dimensioning

      INTEGER       :: MAXLAMBDAS, MAXLAYERS

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Layer control

      INTEGER                                :: NLAYERS
       
!  Output arguments
!  ----------------

!  Gas cross-sections

      REAL(fpk),    DIMENSION( MAXLAYERS, MAXLAMBDAS ) :: SO2_CROSSSECS 

!  status

      LOGICAL       ::        FAIL
      CHARACTER*(*) ::        MESSAGE

!  Local variables
!  ===============

!  Dimensioning

      INTEGER, parameter ::   SO2XS_MAX_DATAPOINTS = 1500

!  Cross section data
!  ------------------

    real(fpk),   dimension (SO2XS_MAX_DATAPOINTS) :: SO2XS_WAVEDATA 
    real(fpk),   dimension (SO2XS_MAX_DATAPOINTS) :: SO2XS_XSECDATA
    real(fpk),   dimension (SO2XS_MAX_DATAPOINTS) :: SO2XS_XSEC2

!  help variables

      INTEGER      :: W, I, M, NP
      real(fpk)    :: LOCAL_SO2_CROSSSEC

!  initialise file read

      FAIL = .FALSE.

!  Open file

      OPEN ( 55, FILE = FILENAME, ERR = 900, STATUS = 'OLD' )

!  main read section for basic data file

      READ(55,*)NP
      IF ( NP .GT. SO2XS_MAX_DATAPOINTS ) GOTO 902
      DO M = 1, NP
        READ(55,*) SO2XS_WAVEDATA(M),SO2XS_XSECDATA(M)
      ENDDO

!  Spline interpolate immediately

      CALL SPLINE ( SO2XS_WAVEDATA,SO2XS_XSECDATA,     &
                     NP, 0.0D0,0.0D0, SO2XS_XSEC2 )
          
      DO W = 1, NLAMBDAS
        IF ( LAMBDAS(W).LE. SO2XS_WAVEDATA(NP) .AND.  &
             LAMBDAS(W).GE. SO2XS_WAVEDATA(1) ) THEN 
          CALL SPLINT ( SO2XS_WAVEDATA,SO2XS_XSECDATA,  &
                        SO2XS_XSEC2,NP,LAMBDAS(W),LOCAL_SO2_CROSSSEC )
        ELSE
          LOCAL_SO2_CROSSSEC = 0.0d0
        ENDIF
        DO I = 1, NLAYERS
          SO2_CROSSSECS(I,W) = LOCAL_SO2_CROSSSEC
        ENDDO
      ENDDO

!  CLose file

      CLOSE ( 55 )

!  Finish

      RETURN

!  Open-file error return

900   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Could not find file of SO2 cross-sections to open'
      CLOSE(55)
      RETURN

!  dimensioning error return

902   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Dimensioning insufficient for SO2 cross-section data'
      CLOSE(55)
      RETURN

END SUBROUTINE GET_SO2XSEC_1

!

SUBROUTINE GET_SO2XSEC_2 ( FILENAME, MAXLAMBDAS, MAXLAYERS,  & 
                           NLAMBDAS, NLAYERS, LAMBDAS,       &
                           so2_crosssecs, fail, message )

!  Dimensioning

      INTEGER       :: MAXLAMBDAS, MAXLAYERS

!  Filename

      character*(*) :: FILENAME

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Layer control

      INTEGER                                :: NLAYERS
       
!  Output arguments
!  ----------------

!  Gas cross-sections

      REAL(fpk),    DIMENSION( MAXLAYERS, MAXLAMBDAS ) :: SO2_CROSSSECS 

!  Exceptiopn handling

   logical             :: FAIL
   character*(*)       :: MESSAGE

!  Local stuff
!  -----------

   integer, parameter :: maxspline = 1000
   logical       :: reading
   integer       :: i, nbuff, ndata, n, maxdata
   real(fpk)     :: lamb1, lamb2, wav, val, conversion, xsec
   real(fpk)     :: x(maxspline), y(maxspline), y2(maxspline)

!  initialize exception handling

   fail    = .false.
   message = ' '

!  Buffer control, initialization

   nbuff = 0
   ndata = 1
   lamb1 = lambdas(1)        - 0.05d0
   lamb2 = lambdas(nlambdas) + 0.05d0
   maxdata    = 23450
   conversion = 1.0d+20

!  Read the data, check that window is covered

    open(1,file=adjustl(trim(filename)),err=90,status='old')

!  .. Read First line, and check window is covered
         read(1,*)wav, val
         if (lamb1 .lt. wav ) then
            fail = .true.
            message = 'SO2 Xsec data does not cover input window at lower end'
            return
         endif
!  .. Read more lines, saving to spline buffer
         reading = .true.
         do while (reading .and. ndata .lt. maxdata )
            ndata = ndata + 1
            read(1,*) wav, val
            if ( wav .gt. lamb1 ) then
               if ( wav .lt. lamb2 ) then
                  nbuff = nbuff + 1
                  x(nbuff)   = wav
                  y(nbuff)   = val
               endif
               if ( wav .ge. lamb2 ) reading = .false.
            endif
         enddo
!  .. Check if last data line is present, then window not covered
         if (ndata .eq. maxdata) then
            fail = .true.
            message = 'SO2 Xsec data does not cover input window at upper end'
            return
         endif
         close(1)

         if (nbuff .gt. maxspline) then
            fail = .true.
            message = 'need to increase maxspline in GET_SO2XSEC_2'
            return
         endif
         
!        write(*,*)nbuff, lamb1, lamb2
!         do n = 1, nbuff
!            write(*,*)x(n),y(n),yc1(n)
!         enddo
!         pause
         
!  Spline the data

   call spline(x,y,nbuff,0.0d0,0.0d0,y2)
   do n = 1, nlambdas
      call splint(x,y,y2,nbuff,lambdas(n),xsec)
      so2_crosssecs(1,n) = xsec / conversion
      do i = 2, nlayers
         so2_crosssecs(i,n) = so2_crosssecs(1,n)
      enddo
   enddo

!  normal return

   return

!  error return

90 continue
   fail    = .true.
   message = 'Open failure for Xsec file = '//filename(1:LEN(filename))
   return

!  Finish

   return
END SUBROUTINE GET_SO2XSEC_2

!

SUBROUTINE GET_BROXSEC_1 ( FILENAME, MAXLAMBDAS, MAXLAYERS,  & 
                           NLAMBDAS, NLAYERS, LAMBDAS,       &
                           bro_crosssecs, fail, message )

!  Dimensioning

      INTEGER       :: MAXLAMBDAS, MAXLAYERS

!  Filename

      character*(*) :: FILENAME

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Layer control

      INTEGER                                :: NLAYERS
       
!  Output arguments
!  ----------------

!  Gas cross-sections

      REAL(fpk),    DIMENSION( MAXLAYERS, MAXLAMBDAS ) :: BRO_CROSSSECS 

!  Exceptiopn handling

   logical             :: FAIL
   character*(*)       :: MESSAGE

!  Local stuff
!  -----------

   integer, parameter :: maxspline = 1000
   logical       :: reading
   integer       :: i, nbuff, ndata, n, maxdata
   real(fpk)     :: lamb1, lamb2, wav, val, conversion, xsec
   real(fpk)     :: x(maxspline), y(maxspline), y2(maxspline)

!  initialize exception handling

   fail    = .false.
   message = ' '

!  Buffer control, initialization

   nbuff = 0
   ndata = 1
   lamb1 = lambdas(1)        - 0.05d0
   lamb2 = lambdas(nlambdas) + 0.05d0
   maxdata    = 6094
   conversion = 1.0d+20

!  Read the data, check that window is covered

    open(1,file=adjustl(trim(filename)),err=90,status='old')

!  .. Read First line, and check window is covered
         read(1,*)wav, val
         if (lamb1 .lt. wav ) then
            fail = .true.
            message = 'BRO Xsec data does not cover input window at lower end'
            return
         endif
!  .. Read more lines, saving to spline buffer
         reading = .true.
         do while (reading .and. ndata .lt. maxdata )
            ndata = ndata + 1
            read(1,*) wav, val
            if ( wav .gt. lamb1 ) then
               if ( wav .lt. lamb2 ) then
                  nbuff = nbuff + 1
                  x(nbuff)   = wav
                  y(nbuff)   = val
               endif
               if ( wav .ge. lamb2 ) reading = .false.
            endif
         enddo
!  .. Check if last data line is present, then window not covered
         if (ndata .eq. maxdata) then
            fail = .true.
            message = 'BRO Xsec data does not cover input window at upper end'
            return
         endif
         close(1)

         if (nbuff .gt. maxspline) then
            fail = .true.
            message = 'need to increase maxspline in GET_BROXSEC_1'
            return
         endif
         
!        write(*,*)nbuff, lamb1, lamb2
!         do n = 1, nbuff
!            write(*,*)x(n),y(n),yc1(n)
!         enddo
!         pause
         
!  Spline the data

   call spline(x,y,nbuff,0.0d0,0.0d0,y2)
   do n = 1, nlambdas
      call splint(x,y,y2,nbuff,lambdas(n),xsec)
      if ( xsec.ge.0.0d0 ) bro_crosssecs(1,n) = xsec / conversion
      do i = 2, nlayers
         bro_crosssecs(i,n) = bro_crosssecs(1,n)
      enddo
   enddo

!  normal return

   return

!  error return

90 continue
   fail    = .true.
   message = 'Open failure for Xsec file = '//filename(1:LEN(filename))
   return

!  Finish

   return
END SUBROUTINE GET_BROXSEC_1

!

SUBROUTINE GET_O2O2XSEC_1                             &
          ( FILENAME, MAXLAMBDAS, MAXLAYERS,          & 
            NLAMBDAS, NLAYERS, LAMBDAS,               &
            o2o2_crosssecs, FAIL, message )

!  file detail

      CHARACTER*(*) ::   FILENAME

!  Dimensioning

      INTEGER       :: MAXLAMBDAS, MAXLAYERS

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Layer control

      INTEGER                                :: NLAYERS
       
!  Output arguments
!  ----------------

!  Gas cross-sections

      REAL(fpk),    DIMENSION( MAXLAYERS, MAXLAMBDAS ) :: O2O2_CROSSSECS 

!  status

      LOGICAL       ::        FAIL
      CHARACTER*(*) ::        MESSAGE

!  Local variables
!  ===============

!  Dimensioning

      INTEGER, parameter ::   O2O2XS_MAX_DATAPOINTS = 3090

!  Cross section data
!  ------------------

    real(fpk),   dimension (O2O2XS_MAX_DATAPOINTS) :: O2O2XS_WAVEDATA 
    real(fpk),   dimension (O2O2XS_MAX_DATAPOINTS) :: O2O2XS_XSECDATA
    real(fpk),   dimension (O2O2XS_MAX_DATAPOINTS) :: O2O2XS_XSEC2

!  help variables

      INTEGER      :: W, I, M, NP
      real(fpk)    :: LOCAL_O2O2_CROSSSEC

!  initialise file read

      FAIL = .FALSE.

!  Open file

      OPEN ( 55, FILE = FILENAME, ERR = 900, STATUS = 'OLD' )

!  main read section for basic data file

      READ(55,*)NP
      IF ( NP .GT. O2O2XS_MAX_DATAPOINTS ) GOTO 902
      DO M = 1, NP
         READ(55,*) O2O2XS_WAVEDATA(M),O2O2XS_XSECDATA(M)
         O2O2XS_WAVEDATA(M) = O2O2XS_WAVEDATA(M)/10.0_fpk
      ENDDO

!  Spline interpolate immediately

      CALL SPLINE ( O2O2XS_WAVEDATA,O2O2XS_XSECDATA,     &
                     NP, 0.0D0,0.0D0, O2O2XS_XSEC2 )
          
      DO W = 1, NLAMBDAS
        IF ( LAMBDAS(W).LE. O2O2XS_WAVEDATA(NP) .AND.  &
             LAMBDAS(W).GE. O2O2XS_WAVEDATA(1) ) THEN 
          CALL SPLINT ( O2O2XS_WAVEDATA,O2O2XS_XSECDATA,  &
                        O2O2XS_XSEC2,NP,LAMBDAS(W),LOCAL_O2O2_CROSSSEC )
        ELSE
          LOCAL_O2O2_CROSSSEC = 0.0d0
        ENDIF
        DO I = 1, NLAYERS
          O2O2_CROSSSECS(I,W) = LOCAL_O2O2_CROSSSEC
        ENDDO
      ENDDO

!  CLose file

      CLOSE ( 55 )

!  Finish

      RETURN

!  Open-file error return

900   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Could not find file of O2O2 cross-sections to open'
      CLOSE(55)
      RETURN

!  dimensioning error return

902   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Dimensioning insufficient for O2O2 cross-section data'
      CLOSE(55)
      RETURN

END SUBROUTINE GET_O2O2XSEC_1
!

SUBROUTINE GET_O2O2XSEC_2 ( FILENAME, MAXLAMBDAS, MAXLAYERS,  & 
                            NLAMBDAS, NLAYERS, LAMBDAS,       &
                            o2o2_crosssecs, fail, message )

!  Dimensioning

      INTEGER       :: MAXLAMBDAS, MAXLAYERS

!  Filename

      character*(*) :: FILENAME

!  wavelengths
 
      INTEGER                                :: NLAMBDAS
      real(fpk),    dimension ( MAXLAMBDAS ) :: LAMBDAS 

!  Layer control

      INTEGER                                :: NLAYERS
       
!  Output arguments
!  ----------------

!  Gas cross-sections

      REAL(fpk),    DIMENSION( MAXLAYERS, MAXLAMBDAS ) :: O2O2_CROSSSECS 

!  Exceptiopn handling

   logical             :: FAIL
   character*(*)       :: MESSAGE

!  Local stuff
!  -----------

   integer, parameter :: maxspline = 1000
   logical       :: reading
   integer       :: i, nbuff, ndata, n, maxdata
   real(fpk)     :: lamb1, lamb2, wav, val, conversion, xsec
   real(fpk)     :: x(maxspline), y(maxspline), y2(maxspline)

!  initialize exception handling

   fail    = .false.
   message = ' '

!  Buffer control, initialization

   nbuff = 0
   ndata = 1
   lamb1 = lambdas(1)        - 0.05d0
   lamb2 = lambdas(nlambdas) + 0.05d0
   maxdata    = 12215
   conversion = 1.0d+40

!  Read the data, check that window is covered

    open(1,file=adjustl(trim(filename)),err=90,status='old')

!  .. Read First line, and check window is covered
         read(1,*)wav, val
         if (lamb1 .lt. wav ) then
            fail = .true.
            message = 'O2O2 Xsec data does not cover input window at lower end'
            return
         endif
!  .. Read more lines, saving to spline buffer
         reading = .true.
         do while (reading .and. ndata .lt. maxdata )
            ndata = ndata + 1
            read(1,*) wav, val
            if ( wav .gt. lamb1 ) then
               if ( wav .lt. lamb2 ) then
                  nbuff = nbuff + 1
                  x(nbuff)   = wav
                  y(nbuff)   = val
               endif
               if ( wav .ge. lamb2 ) reading = .false.
            endif
         enddo
!  .. Check if last data line is present, then window not covered
         if (ndata .eq. maxdata) then
            fail = .true.
            message = 'O2O2 Xsec data does not cover input window at upper end'
            return
         endif
         close(1)

         if (nbuff .gt. maxspline) then
            fail = .true.
            message = 'need to increase maxspline in GET_O2O2XSEC_2'
            return
         endif
         
!        write(*,*)nbuff, lamb1, lamb2
!         do n = 1, nbuff
!            write(*,*)x(n),y(n),yc1(n)
!         enddo
!         pause
         
!  Spline the data

   call spline(x,y,nbuff,0.0d0,0.0d0,y2)
   do n = 1, nlambdas
      call splint(x,y,y2,nbuff,lambdas(n),xsec)
      if ( xsec .gt. 0.0d0 ) o2o2_crosssecs(1,n) = xsec / conversion
      do i = 2, nlayers
         o2o2_crosssecs(i,n) = o2o2_crosssecs(1,n)
      enddo
   enddo

!  normal return

   return

!  error return

90 continue
   fail    = .true.
   message = 'Open failure for Xsec file = '//filename(1:LEN(filename))
   return

!  Finish

   return
END SUBROUTINE GET_O2O2XSEC_2

!  End module

End module cross_sections_m

