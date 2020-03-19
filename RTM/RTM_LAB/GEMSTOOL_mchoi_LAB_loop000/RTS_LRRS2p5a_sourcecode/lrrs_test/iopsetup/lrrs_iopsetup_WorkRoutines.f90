
! ###############################################################
! #                                                             #
! #                    THE LIDORT_RRS MODEL                     #
! #                                                             #
! #      (LInearized Discrete Ordinate Radiative Transfer)      #
! #       --         -        -        -         -              #
! #                 (Rotational Raman Scattering)               #
! #                  -          -     -                         #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert J. D. Spurr                           #
! #                                                             #
! #  Address :     RT SOLUTIONS Inc.                            #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email   :     rtsolutions@verizon.net                      #
! #  Website :     www.rtslidort.com                            #
! #                                                             #
! #  Version  #   :  2.5                                        #
! #  Release Date :  March 2017                                 #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  --- History of the model ------------                      #
! #                                                             #
! #  Version 1.0 : 2005, Fortran 77                             #
! #  Version 1.1 : 2007, F77                                    #
! #                                                             #
! #  Version 2.1 : 2009, F77                                    #
! #       * Linearization for Atmospheric-property Jacobians    #
! #       * Single scatter corrections added                    #
! #                                                             #
! #  Version 2.3 : March 2011, Fortran 90                       #
! #       * Simplified Raman-setup procedure                    #
! #       * F90 Version with Type-structure I/O                 #
! #       * Test package developed for installation             #
! #                                                             #
! #  Version 2.5 : March 2017, F90                              #
! #       * Formal BRDF/SLEAVE supplements developed            #
! #       * New test-bed software for testing supplements       #
! #       * Thread-safe Code for OpenMP applications            #
! #       * Complete revision of Taylor-series modules          #
! #       * New User Guide and Review paper                     #
! #                                                             #
! ###############################################################

!    #########################################################
!    #                                                       #
!    #   This Version of LIDORT_RRS comes with a GNU-style   #
!    #   license. Please read the license carefully.         #
!    #                                                       #
!    #########################################################

Module lrrs_iopsetup_WorkRoutines_m

!   Auxiliary routines must be used here.

   use LRRS_Iopsetup_Aux_m, Only : LINEAR_TRAWLER_ASCENDING, SPLINE, SPLINT, LINTP2

!  Contains the following stand-alone routines, all public
!  =======================================================

!    Name                        Called by                   Data Source
!    ----                        ---------                   -----------

! GET_HEIGHTS_A		     O3MONO, O3BIN                  A. Vasilkov
! GET_PTAO3                  O3MONO, O3BIN                  A. Vasilkov

! GET_TV8O3_GM1_Part1        LCSMONO                        TOMSV8 O3 data
! GET_TV8O3_GM1_Part2        LCSMONO                        TOMSV8 O3 data
! GET_TV8O3_GM1_All          LPMONO, LCSBIN, LPBIN          TOMSV8 O3 data

! GET_O3DATA_GOME_Part1      O3MONO                         GOME-1 FM Xsections
! GET_O3XSEC_GOME_Part2      O3MONO                         GOME-1 FM Xsections
! GET_O3XSEC_GOME_All        O3BIN                          GOME-1 FM Xsections

! RAYLEIGH_FUNCTION          All Iopsetup tools             Bodhaine et al. (1999)

! BREON_O3DATA_Part1         LCSMONO, LPMONO                Breon CDM data
! BREON_O3XSEC_Part2         LCSMONO, LPMONO                Breon CDM data
! BREON_O3XSEC_All           LCSBIN,  LPBIN                 Breon CDM data

! GET_O3XSEC_1               Not Used                       Not known

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IMPLICIT NONE

      PRIVATE :: FPK
      PUBLIC

      INTEGER, PARAMETER :: FPK = SELECTED_REAL_KIND(15)

!  TOMS TV8O3 data type structure

      TYPE TV8O3_Data
        REAL ::      PROFILES   (12,18,11,10)
        REAL ::      COVMATRIX  (11,11)
        REAL ::      COLUMNS    (12,18,10)
        INTEGER ::   COLNUMBERS (12,18)
        REAL ::      PRESSURES  (14)
      END TYPE TV8O3_Data

CONTAINS

      SUBROUTINE GET_HEIGHTS_A &
        ( FILENAME, MAX_LAYERS,   & !Inputs
          INPUT_HEIGHTS,          & !Inout
          NLAYERS, FAIL, MESSAGE )  !Outputs

      implicit none
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  File detail

      CHARACTER (LEN=*), INTENT(IN) ::  FILENAME

!  Dimensioning

      INTEGER, INTENT(IN) ::            MAX_LAYERS

!  Output
!  ------

!  Heights and number of layers

      INTEGER, INTENT(OUT) ::           NLAYERS
      REAL(FPK), INTENT(INOUT) ::       INPUT_HEIGHTS ( 0:MAX_LAYERS )

!  Status

      LOGICAL, INTENT(OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE

!  Help variables
!  --------------

      INTEGER ::        K

!  Initialise file read

      FAIL = .FALSE.

!  Open file

      OPEN ( UNIT   = 55, &
             FILE   = Adjustl(Trim(FILENAME)), &
             ERR    = 900, &
             STATUS = 'OLD' )

!  Read height data

      READ(55,*)NLAYERS
      IF ( NLAYERS .GT. MAX_LAYERS ) GOTO 902
      DO K = 0, NLAYERS
        READ(55,*)INPUT_HEIGHTS(K)
      ENDDO

!  Normal return

      CLOSE(55)
      RETURN

!  Open-file error return

900   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Could not find file of input heights to open'
      CLOSE(55)
      RETURN

!  Dimensioning error return

902   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Dimensioning insufficient for heights data'
      CLOSE(55)
      RETURN

      END SUBROUTINE GET_HEIGHTS_A

!  P, T, Air, and O3

      SUBROUTINE GET_PTAO3 &
        ( O3FILENAME, ACTUAL_COLUMN,                 & !Inputs
          MAX_LAYERS, NLAYERS, INPUT_HEIGHTS,        & !Inputs
          LAYER_TEMPERATURES, AIRCOLUMNS, O3UMKEHRS, & !Outputs
          FAIL, MESSAGE )                              !Outputs

      implicit none
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  File detail

      CHARACTER (LEN=*), INTENT(IN)  :: O3FILENAME

!  Actual column amount of ozone in [DU]

      REAL(FPK), INTENT(IN)          :: ACTUAL_COLUMN

!  Dimensioning

      INTEGER, INTENT(IN)            :: MAX_LAYERS

!  Height grid input

      INTEGER, INTENT(IN)            :: NLAYERS
      REAL(FPK), INTENT(IN)          :: INPUT_HEIGHTS ( 0:MAX_LAYERS )

!  Outputs
!  -------

!  Temperatures, air and gas column amounts

      REAL(FPK), INTENT(INOUT)       :: LAYER_TEMPERATURES ( MAX_LAYERS )
      REAL(FPK), INTENT(INOUT)       :: AIRCOLUMNS         ( MAX_LAYERS )

      REAL(FPK), INTENT(OUT)         :: O3UMKEHRS          ( MAX_LAYERS )

!  Status

      LOGICAL, INTENT(OUT)           :: FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE

!  Local variables
!  ---------------

!  Dimensioning

      INTEGER, PARAMETER :: LOCAL_MAXLEVELS  = 61
      INTEGER, PARAMETER :: LOCAL_MAXCOLUMNS = 11

!  Local arrays

      INTEGER   :: LOCAL_NLEVELS
      REAL(FPK) :: LOCAL_HEIGHTS ( LOCAL_MAXLEVELS )
      REAL(FPK) :: LOCAL_TEMPS   ( LOCAL_MAXLEVELS )
      REAL(FPK) :: LOCAL_PRESS   ( LOCAL_MAXLEVELS )
      REAL(FPK) :: LOCAL_O3VMRS  ( LOCAL_MAXLEVELS )

!  USA data arrays

      INTEGER   :: USADATA_NLEVELS, ZONE
      REAL(FPK) :: USADATA_PRESS   ( 46 ), PRESS(6)
      REAL(FPK) :: USADATA_HEIGHTS ( 46 )
      REAL(FPK) :: USADATA_TEMPS   ( 46 ), TEMP(6)

!  O3 data arrays

      INTEGER   :: NCOLUMNS, NO3_NLEVELS
      REAL(FPK) :: O3DATA_HEIGHTS ( LOCAL_MAXLEVELS )
      REAL(FPK) :: O3DATA_VMRS    ( LOCAL_MAXLEVELS )
      REAL(FPK) :: COLUMNS ( LOCAL_MAXCOLUMNS )
      REAL(FPK) :: TEMPO   ( LOCAL_MAXCOLUMNS )

!  Help variables

      REAL(FPK) :: F1, F2, GDENS1, GDENS2, RAT1, RAT2
      REAL(FPK) :: ADENS1, ADENS2, DIFFZ
      CHARACTER (LEN=70) ::     USA_FILENAME
      INTEGER :: DUM, K1, K2, C, N, N1, R, I, K

!  Loschmidt's number (particles/cm3)

      REAL(FPK), PARAMETER :: RHO_STANDARD = 2.68675D+19

!  STP values

      REAL(FPK), PARAMETER :: PZERO = 1013.25D0, TZERO = 273.15D0
      REAL(FPK), PARAMETER :: PTZERO = PZERO / TZERO

      REAL(FPK), PARAMETER :: RHO_ZERO = RHO_STANDARD*TZERO/PZERO
      REAL(FPK), PARAMETER :: CONSTANT = 1.0D+05*RHO_ZERO
      REAL(FPK), PARAMETER :: PPMV = 1.0D-06

!  Initialise file read

      FAIL = .FALSE.
      MESSAGE = ' '

!  Local heights

      LOCAL_NLEVELS = NLAYERS + 1
      DO N = 0, NLAYERS
        LOCAL_HEIGHTS(N+1) = INPUT_HEIGHTS(N)
      ENDDO

!  Choose USA profiles to use, 6 = USA standard

      ZONE = 6
      R = 1

!  Read USA data

      USA_FILENAME = 'lrrs_test/physics_data/ATMOS/pth_usa.dat'
      OPEN ( UNIT   = 55, &
             FILE   = Adjustl(Trim(USA_FILENAME)), &
             ERR    = 900, &
             STATUS = 'OLD' )
      READ(55,*)
      READ(55,*)USADATA_NLEVELS
      DO I = 1, USADATA_NLEVELS
        READ(55,*) USADATA_HEIGHTS(I),(PRESS(K),K=1,6)
          USADATA_PRESS(I) = DLOG(PRESS(ZONE))
      ENDDO
      DO I = 1, USADATA_NLEVELS
        READ(55,*) USADATA_HEIGHTS(I),(TEMP(K),K=1,6)
          USADATA_TEMPS(I) = TEMP(ZONE)
      ENDDO

!  Normal finish

        CLOSE(55)

!  Interpolate to incoming height grid

        CALL LINTP2 &
         ( USADATA_NLEVELS,USADATA_HEIGHTS,USADATA_TEMPS, &
           LOCAL_NLEVELS, LOCAL_HEIGHTS, LOCAL_TEMPS )

        CALL LINTP2 &
         ( USADATA_NLEVELS,USADATA_HEIGHTS,USADATA_PRESS, &
           LOCAL_NLEVELS, LOCAL_HEIGHTS, LOCAL_PRESS )

!  Read ozone file

      R = 2
      OPEN ( UNIT   = 55, &
             FILE   = Adjustl(Trim(O3FILENAME)), &
             ERR    = 900, &
             STATUS = 'OLD' )

      DO K = 1, 6
        READ(55,*)
      ENDDO
      READ(55,*)NCOLUMNS,NO3_NLEVELS
      READ(55,*)(COLUMNS(C),C=1,NCOLUMNS)
      CALL LINEAR_TRAWLER_ASCENDING &
          ( LOCAL_MAXCOLUMNS, NCOLUMNS, COLUMNS, ACTUAL_COLUMN, &
            F1, F2, K1, K2 )
      DO N = 1, NO3_NLEVELS
        N1 = NO3_NLEVELS + 1 - N
        READ(55,*)DUM,(TEMPO(C),C = 1, NCOLUMNS)
        O3DATA_VMRS(N1)    = ( F1*TEMPO(K2) + F2*TEMPO(K1) ) * PPMV
        O3DATA_HEIGHTS(N1) = DBLE(DUM)
      ENDDO

!  Interpolate O3 VMR to incoming height grid

        CALL LINTP2 &
         ( NO3_NLEVELS,O3DATA_HEIGHTS,O3DATA_VMRS, &
           LOCAL_NLEVELS, LOCAL_HEIGHTS, LOCAL_O3VMRS )

!  Set Column Air density and Umkehr amount

      DO N = 1, NLAYERS
        N1 = N + 1
        DIFFZ = LOCAL_HEIGHTS(N) - LOCAL_HEIGHTS(N1)
        RAT1 = DEXP(LOCAL_PRESS(N))  / LOCAL_TEMPS(N)
        RAT2 = DEXP(LOCAL_PRESS(N1)) / LOCAL_TEMPS(N1)
        ADENS1 = CONSTANT *  RAT1
        GDENS1 = ADENS1 * LOCAL_O3VMRS(N)
        ADENS2 = CONSTANT *  RAT2
        GDENS2 = ADENS2 * LOCAL_O3VMRS(N1)
        LAYER_TEMPERATURES(N)= 0.5d0*(LOCAL_TEMPS(N)+LOCAL_TEMPS(N1))
        AIRCOLUMNS(N)        = 0.5d0 * ( ADENS1 + ADENS2 ) * DIFFZ
        O3UMKEHRS(N)         = 0.5d0 * ( GDENS1 + GDENS2 ) * DIFFZ
      ENDDO

!  Normal return

      RETURN

!  Open-file error return

900   CONTINUE
      FAIL = .TRUE.
      IF (R.EQ.1) MESSAGE = 'No file of USA PTH data to open'
      IF (R.EQ.2) MESSAGE = 'No file of O3 VMR data to open'
      CLOSE(55)
      RETURN

!  Dimensioning error return
!      FAIL = .TRUE.
!      MESSAGE = 'not enough dimensions'
!      CLOSE(55)
!      RETURN

      END SUBROUTINE GET_PTAO3

!  Get ozone concentration data: TOMS (Part 1)

      SUBROUTINE GET_TV8O3_GM1_Part1 ( TV8O3NEW_DATA, FAIL, MESSAGE ) !Output

!mick mod 2/27/2017
!  Note: Subroutine GET_TV8O3_GM1 was split into two for this module
!        (a file read portion and an interpolation portion);
!        GET_TV8O3_GM1_Part1 is the file read portion

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  None

!  Output arguments
!  ----------------

!  TOMS data on 4-D grid

      TYPE(TV8O3_Data), INTENT(OUT) :: TV8O3NEW_DATA

!  Status

      LOGICAL, INTENT(OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE

!  Local variables
!  ---------------

      INTEGER ::   C, I, K, M, NC, Z, UNUM
      REAL ::      RCOL, PROF(11)

!  Initialize status

      FAIL = .FALSE.
      MESSAGE = ' '

!  Initialize output

      TV8O3NEW_DATA%PROFILES   = 0.0_fpk
      TV8O3NEW_DATA%COVMATRIX  = 0.0_fpk
      TV8O3NEW_DATA%COLUMNS    = 0.0_fpk
      TV8O3NEW_DATA%COLNUMBERS = 0
      TV8O3NEW_DATA%PRESSURES  = 0.0_fpk

!  Read some data

      UNUM = 99
      OPEN ( UNIT   = UNUM, &
             FILE   = 'lrrs_test/physics_data/ATMOS/tomsv8_ozoneprofiles_new.dat', &
             ERR    = 900, &
             STATUS = 'OLD' )
      READ (UNUM,*)
      READ (UNUM,*)
      READ (UNUM,*,ERR=901)(TV8O3NEW_DATA%PRESSURES(K),K=1,14)
      READ (UNUM,*)
      DO I = 1, 11
         READ (UNUM,*,ERR=901)(TV8O3NEW_DATA%COVMATRIX(I,K),K=1,11)
      ENDDO
      DO M = 1, 12 !Month
         DO Z = 1, 18  !Latitude
            READ(UNUM,*)        !The "Month = M Latitude =  Z" lines
            NC = 0
            DO C = 1, 10 !Pressure level
               READ(UNUM,'(F7.1,11F7.2)')RCOL,(PROF(K),K=11,1,-1)
               IF ( RCOL .NE. 999.0 ) THEN
                  NC = NC + 1
                  TV8O3NEW_DATA%COLUMNS(M,Z,NC) = RCOL
                  DO K = 1, 11
                     TV8O3NEW_DATA%PROFILES(M,Z,K,NC) = PROF(K)
                  ENDDO
               ENDIF
               TV8O3NEW_DATA%COLNUMBERS(M,Z) = NC
            ENDDO
         ENDDO
      ENDDO
      CLOSE (UNUM)

!  Successful file-read

      RETURN

!  Open-file error return

 900  CONTINUE
      FAIL = .TRUE.
      message = 'TOMS V8 Extraction: data file-open failure'
      CLOSE(UNUM)
      RETURN

!  Read-file error return

 901  CONTINUE
      FAIL = .TRUE.
      message = 'TOMS V8 Extraction: data file-read failure'
      CLOSE(UNUM)
      RETURN

!  Finish

      END SUBROUTINE GET_TV8O3_GM1_Part1

!  Get ozone concentration data: TOMS (Part 2)

      SUBROUTINE GET_TV8O3_GM1_Part2 &
          ( ACTUAL_O3COLUMN, LATITUDE, YEAR, MONTH, DAY_OF_MONTH, & !Input
            MAXLAYERS, MAXLEVELS, NLAYERS,                        & !Input
            LEVEL_LOGPRESSURES, SURFACE_PRESS, TV8O3NEW_DATA,     & !Input
            O3UMKEHRS, DC_O3UMKEHRS,                              & !Output
            FAIL, MESSAGE )                                         !Output

!mick mod 2/27/2017
!  Note: Subroutine GET_TV8O3_GM1 was split into two for this module
!        (a file read portion and an interpolation portion);
!        GET_TV8O3_GM1_Part2 is the interpolation portion

      IMPLICIT NONE
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input arguments
!  ---------------

!  Actual O3 column

      REAL(FPK), INTENT(IN) :: ACTUAL_O3COLUMN

!  Time data

      REAL(FPK), INTENT(IN) :: LATITUDE
      INTEGER, INTENT(IN) ::   YEAR, DAY_OF_MONTH

      INTEGER, INTENT(INOUT) :: MONTH

!  Dimensioning

      INTEGER, INTENT(IN) ::   MAXLAYERS, MAXLEVELS

!  Pressure grid input

      INTEGER, INTENT(IN) ::   NLAYERS
      REAL(FPK), INTENT(INOUT) :: LEVEL_LOGPRESSURES ( MAXLEVELS )
      REAL(FPK), INTENT(IN) :: SURFACE_PRESS

!  TOMS data on 4-D grid

      TYPE(TV8O3_Data), INTENT(IN) :: TV8O3NEW_DATA

!  Output arguments
!  ----------------

!  Gas column amounts

      REAL(FPK), INTENT(OUT) :: O3UMKEHRS    ( MAXLAYERS )
      REAL(FPK), INTENT(OUT) :: DC_O3UMKEHRS ( MAXLAYERS )

!  Status

      LOGICAL, INTENT(OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE

!  Local variables
!  ---------------

!  Output mapping coefficients for clear-sky cases
!  --> profile sets which define the mapping domain
!  reference columns spanning the mapping domain

      REAL(FPK) :: CUMO3_MAPPING(14,10)
      INTEGER ::   O3MAPPING_N_REFCOLS
      REAL(FPK) :: O3MAPPING_REFCOLS (10)

!  Interpolation variables

      INTEGER ::   N1, N2, NLEVELS
      REAL(FPK) :: DF1, DF2, GRADN, CGRAD, CGRADP, DIFF, XFAC, LOCAL
      REAL(FPK) :: LOCAL_O31(105), LOCAL_O32(105), SUM
      INTEGER ::   I, NC, C, Z, UNUM
      REAL ::      RCOL, PROF(11)

!  Time/latitude help variables

      LOGICAL ::   LOOP, DO_INTERP
      REAL ::      SLOPE, DIST, O3(11,10,2)
      INTEGER ::   K, N, ZONE, M, REF, TIME,NBDAY(12)
      REAL(FPK) :: COL, DEL, PLAT1
      REAL(FPK) :: PROFILE (14), CUMO3_LNPRESS(14),LATZONES(19)

!  TOMS version 8 profiles (New versions)

      REAL ::      TV8O3NEW_DATA_PROFILES   (12,18,11,10)
      REAL ::      TV8O3NEW_DATA_COVMATRIX  (11, 11)
      REAL ::      TV8O3NEW_DATA_COLUMNS    (12,18,10)
      INTEGER ::   TV8O3NEW_DATA_COLNUMBERS (12,18)
      REAL ::      TV8O3NEW_DATA_PRESSURES  (14)

!  Initialize status

      FAIL = .FALSE.
      MESSAGE = ' '

!mick fix 7/20/2016 - Initialize output

      O3UMKEHRS    = 0.0_FPK
      DC_O3UMKEHRS = 0.0_FPK

!  Define some arrays

      NBDAY = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

      LATZONES = (/ &
           -90.0, -80.0, -70.0, -60.0, -50.0, &
           -40.0, -30.0, -20.0, -10.0,   0.0, &
            10.0,  20.0,  30.0,  40.0,  50.0, &
            60.0,  70.0,  80.0,  90.0 /)

! Pass TOMS input to local variables

      TV8O3NEW_DATA_PROFILES   = TV8O3NEW_DATA%PROFILES
      TV8O3NEW_DATA_COVMATRIX  = TV8O3NEW_DATA%COVMATRIX
      TV8O3NEW_DATA_COLUMNS    = TV8O3NEW_DATA%COLUMNS
      TV8O3NEW_DATA_COLNUMBERS = TV8O3NEW_DATA%COLNUMBERS
      TV8O3NEW_DATA_PRESSURES  = TV8O3NEW_DATA%PRESSURES

!  Do the time/geography interpolation

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
                           O3(K,N,TIME)= &
                                TV8O3NEW_DATA_PROFILES(MONTH,1,K,N)+ &
                                SLOPE*DIST
                        ELSE
                           O3(K,N,TIME)= &
                                TV8O3NEW_DATA_PROFILES(MONTH,1,K,N)
                        ENDIF
                     ENDDO
                  ENDDO

               ELSEIF (LATITUDE .GT. 85) THEN
                  DIST=LATITUDE-85
                  O3MAPPING_N_REFCOLS=TV8O3NEW_DATA_COLNUMBERS(MONTH,18)
                  DO N = 1, O3MAPPING_N_REFCOLS
                     COL=TV8O3NEW_DATA_COLUMNS(MONTH,18,N)
                     M=1
                     DO WHILE(TV8O3NEW_DATA_COLUMNS(MONTH,17,M)/=COL &
                          .AND. &
                          M .LT. TV8O3NEW_DATA_COLNUMBERS(MONTH,17))
                        M = M+1
                     ENDDO

                     DO K=1,11
                        IF (TV8O3NEW_DATA_COLUMNS(MONTH,17,M)==COL) THEN
                           SLOPE=(TV8O3NEW_DATA_PROFILES(MONTH,18,K,N) - &
                                TV8O3NEW_DATA_PROFILES(MONTH,17,K,M))/10
                           O3(K,N,TIME)= &
                                TV8O3NEW_DATA_PROFILES(MONTH,18,K,N)+ &
                                SLOPE*DIST
                        ELSE
                           O3(K,N,TIME)= &
                                TV8O3NEW_DATA_PROFILES(MONTH,18,K,N)
                        ENDIF
                     ENDDO
                  ENDDO
               ELSE             !IF LATITUDE
                  ZONE = 0
                  LOOP = .TRUE.
                  DO WHILE (LOOP)
                     ZONE = ZONE + 1
                     IF ( LATITUDE .GT. LATZONES(ZONE)+5 .AND. &
                          LATITUDE .LE. LATZONES(ZONE+1)+5 ) &
                          LOOP = .FALSE.
                  ENDDO
                  DIST=LATITUDE-(LATZONES(ZONE)+5)

                  IF( TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE) >= &
                       TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1)) THEN
                     REF=ZONE
                     O3MAPPING_N_REFCOLS &
                          =TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE)

                     DO N = 1, O3MAPPING_N_REFCOLS
                        COL=TV8O3NEW_DATA_COLUMNS(MONTH,REF,N)
                        M=1
                        DO WHILE( &
                             TV8O3NEW_DATA_COLUMNS(MONTH,REF+1,M)/=COL &
                             .AND. M .LT. &
                             TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1) &
                             )
                           M = M+1
                        ENDDO

                        DO K=1,11
                           IF (TV8O3NEW_DATA_COLUMNS(MONTH,REF+1,M) &
                                ==COL) THEN
                              SLOPE= &
                              (TV8O3NEW_DATA_PROFILES(MONTH,REF+1,K,M) - &
                              TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N))/10
                              O3(K,N,TIME)= &
                                   TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N) &
                                   + SLOPE*DIST
                           ELSE
                              O3(K,N,TIME)= &
                                   TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N)
                           ENDIF
                        ENDDO
                     ENDDO
                  ELSE
                     REF=ZONE+1
                     O3MAPPING_N_REFCOLS= &
                          TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1)

                     DO N = 1, O3MAPPING_N_REFCOLS
                        COL=TV8O3NEW_DATA_COLUMNS(MONTH,REF,N)
                        M=1
                        DO WHILE( &
                             TV8O3NEW_DATA_COLUMNS(MONTH,REF-1,M)/=COL &
                             .AND. M .LT. &
                             TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE-1) &
                             )
                           M = M+1
                        ENDDO

                        DO K=1,11
                           IF (TV8O3NEW_DATA_COLUMNS(MONTH,REF-1,M) &
                                ==COL) THEN
                              SLOPE= &
                                   (TV8O3NEW_DATA_PROFILES &
                                   (MONTH,REF,K,N) - &
                                   TV8O3NEW_DATA_PROFILES &
                                   (MONTH,REF-1,K,M))/10
                              O3(K,N,TIME)= &
                                   TV8O3NEW_DATA_PROFILES &
                                   (MONTH,REF-1,K,M) &
                                   + SLOPE*DIST
                           ELSE
                              O3(K,N,TIME)= &
                                   TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N)
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
                             (O3(K,N,2)-O3(K,N,1))/NBDAY(12)* &
                             (NBDAY(12)-15+DAY_OF_MONTH)
                     ELSE
                        TV8O3NEW_DATA_PROFILES(MONTH,1,K,N) = O3(K,N,1)+ &
                             (O3(K,N,2)-O3(K,N,1))/NBDAY(MONTH-1)* &
                             (NBDAY(MONTH-1)-15+DAY_OF_MONTH)
                     ENDIF
                  ELSE
                     TV8O3NEW_DATA_PROFILES(MONTH,1,K,N) = O3(K,N,1)+ &
                          (O3(K,N,2)-O3(K,N,1))/ &
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

!  set Log pressure grid

      DO K = 1, 14
        CUMO3_LNPRESS(K) = DLOG(DBLE(TV8O3NEW_DATA_PRESSURES(K)))
      ENDDO
      CUMO3_LNPRESS(14) = DLOG(SURFACE_PRESS)


!  Now interpolate to our incoming column

      CALL LINEAR_TRAWLER_ASCENDING &
            ( 10, O3MAPPING_N_REFCOLS, &
              O3MAPPING_REFCOLS, ACTUAL_O3COLUMN, &
              DF1, DF2, N1, N2 )

!  interpolate Cumulative O3 to incoming pressure grid

      nlevels = nlayers + 1
      CALL LINTP2 &
         ( 14,CUMO3_LNPRESS,CUMO3_MAPPING(1,N1), &
           NLEVELS, LEVEL_LOGPRESSURES, LOCAL_O31 )
      CALL LINTP2 &
         ( 14,CUMO3_LNPRESS,CUMO3_MAPPING(1,N2), &
           NLEVELS, LEVEL_LOGPRESSURES, LOCAL_O32 )

      DIFF = O3MAPPING_REFCOLS(N2) - O3MAPPING_REFCOLS(N1)
      XFAC = ACTUAL_O3COLUMN - O3MAPPING_REFCOLS(N1)
      SUM = 0.0D0
      DO N = 1, NLAYERS
        IF ( N.EQ.1 ) THEN
          CGRAD =  ( LOCAL_O32(2) - LOCAL_O31(2) ) / DIFF
          O3UMKEHRS(N)    = LOCAL_O31(2) + XFAC * CGRAD
          DC_O3UMKEHRS(N) = CGRAD
          CGRADP = CGRAD
        ELSE
          CGRAD = ( LOCAL_O32(N+1) - LOCAL_O31(N+1) ) / DIFF
          GRADN = CGRAD - CGRADP
          LOCAL = LOCAL_O31(N+1) - LOCAL_O31(N)
          O3UMKEHRS(N)    =  LOCAL + XFAC * GRADN
          DC_O3UMKEHRS(N) =  GRADN
          CGRADP = CGRAD
        ENDIF
        SUM = SUM + O3UMKEHRS(N)
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE GET_TV8O3_GM1_Part2

!  Get ozone concentration data TOMS (no partition, used in the Bin setups)

      SUBROUTINE GET_TV8O3_GM1_All &
          ( ACTUAL_O3COLUMN, LATITUDE, YEAR, MONTH, DAY_OF_MONTH, & !Input
            MAXLAYERS, MAXLEVELS, NLAYERS,                        & !Input
            LEVEL_LOGPRESSURES, SURFACE_PRESS,                    & !Input
            O3UMKEHRS, DC_O3UMKEHRS,                              & !Output
            FAIL, MESSAGE )                                         !Output

      IMPLICIT NONE
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input arguments
!  ---------------

!  Actual O3 column

      REAL(FPK), INTENT(IN) :: ACTUAL_O3COLUMN

!     Latitude and month

      REAL(FPK), INTENT(IN) :: LATITUDE
      INTEGER, INTENT(IN) ::   YEAR, DAY_OF_MONTH

      INTEGER, INTENT(INOUT) :: MONTH

!  Dimensioning

      INTEGER, INTENT(IN) ::   MAXLAYERS, MAXLEVELS

!  pressure grid input

      INTEGER, INTENT(IN) ::   NLAYERS
      REAL(FPK), INTENT(INOUT) :: LEVEL_LOGPRESSURES ( MAXLEVELS )
      REAL(FPK), INTENT(IN) :: SURFACE_PRESS

!  Output arguments
!  ----------------

!  gas column amounts

      REAL(FPK), INTENT(OUT) :: O3UMKEHRS    ( MAXLAYERS )
      REAL(FPK), INTENT(OUT) :: DC_O3UMKEHRS ( MAXLAYERS )

!  status

      LOGICAL, INTENT(OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE

!  Local variables
!  ---------------

!  output mapping coefficients for clear-sky cases
!  --> profile sets which define the mapping domain
!  reference columns spanning the mapping domain

      REAL(FPK) :: CUMO3_MAPPING(14,10)
      INTEGER ::   O3MAPPING_N_REFCOLS
      REAL(FPK) :: O3MAPPING_REFCOLS (10)

!  Interpolation variables

      INTEGER ::   N1, N2, NLEVELS
      REAL(FPK) :: DF1, DF2, GRADN, CGRAD, CGRADP, DIFF, XFAC, LOCAL
      REAL(FPK) :: LOCAL_O31(105), LOCAL_O32(105), SUM
      INTEGER ::   I, NC, C, Z, UNUM
      REAL ::      RCOL, PROF(11)

!  Time/latitude help variables

      LOGICAL ::   LOOP, DO_INTERP
      REAL ::      SLOPE, DIST, O3(11,10,2)
      INTEGER ::   K, N, ZONE, M, REF, TIME,NBDAY(12)
      REAL(FPK) :: COL, DEL, PLAT1
      REAL(FPK) :: PROFILE (14), CUMO3_LNPRESS(14),LATZONES(19)

!  output TOMS version 8 profiles (New versions)

      REAL ::      TV8O3NEW_DATA_PROFILES   (12,18,11,10)
      REAL ::      TV8O3NEW_DATA_COVMATRIX  (11, 11)
      REAL ::      TV8O3NEW_DATA_COLUMNS    (12,18,10)
      INTEGER ::   TV8O3NEW_DATA_COLNUMBERS (12,18)
      REAL ::      TV8O3NEW_DATA_PRESSURES  (14)

!  Initialize status

      FAIL = .FALSE.
      MESSAGE = ' '

!mick fix 7/20/2016 - Initialize output

      O3UMKEHRS    = 0.0_FPK
      DC_O3UMKEHRS = 0.0_FPK

!  Define some arrays

      NBDAY = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

      LATZONES = (/ &
           -90.0, -80.0, -70.0, -60.0, -50.0, &
           -40.0, -30.0, -20.0, -10.0,   0.0, &
            10.0,  20.0,  30.0,  40.0,  50.0, &
            60.0,  70.0,  80.0,  90.0 /)

!  Read some data

      UNUM = 99
      OPEN ( UNIT   = UNUM, &
             FILE   = 'lrrs_test/physics_data/ATMOS/tomsv8_ozoneprofiles_new.dat', &
             ERR    = 900, &
             STATUS = 'OLD' )
      READ (UNUM,*)
      READ (UNUM,*)
      READ (UNUM,*,ERR=901)(TV8O3NEW_DATA_PRESSURES(K),K=1,14)
      READ (UNUM,*)
      DO I = 1, 11
         READ (UNUM,*,ERR=901)(TV8O3NEW_DATA_COVMATRIX(I,K),K=1,11)
      ENDDO
      DO M = 1, 12 !Month
         DO Z = 1, 18  !Latitude
            READ(UNUM,*)        !The "Month = M Latitude =  Z" lines
            NC = 0
            DO C = 1, 10 !Pressure level
               READ(UNUM,'(F7.1,11F7.2)')RCOL,(PROF(K),K=11,1,-1)
               IF ( RCOL .NE. 999.0 ) THEN
                  NC = NC + 1
                  TV8O3NEW_DATA_COLUMNS(M,Z,NC) = RCOL
                  DO K = 1, 11
                     TV8O3NEW_DATA_PROFILES(M,Z,K,NC) = PROF(K)
                  ENDDO
               ENDIF
               TV8O3NEW_DATA_COLNUMBERS(M,Z) = NC
            ENDDO
         ENDDO
      ENDDO
      CLOSE (UNUM)

!  Successful file-read

      GO TO 676

!  Open-file error return

 900  CONTINUE
      FAIL = .TRUE.
      message = 'TOMS V8 Extraction: data file-open failure'
      CLOSE(UNUM)
      RETURN

!  Read-file error return

 901  CONTINUE
      FAIL = .TRUE.
      message = 'TOMS V8 Extraction: data file-read failure'
      CLOSE(UNUM)
      RETURN

!  Coninuation point

 676  continue

!  ===============================================

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
                           O3(K,N,TIME)= &
                                TV8O3NEW_DATA_PROFILES(MONTH,1,K,N)+ &
                                SLOPE*DIST
                        ELSE
                           O3(K,N,TIME)= &
                                TV8O3NEW_DATA_PROFILES(MONTH,1,K,N)
                        ENDIF
                     ENDDO
                  ENDDO

               ELSEIF (LATITUDE .GT. 85) THEN
                  DIST=LATITUDE-85
                  O3MAPPING_N_REFCOLS=TV8O3NEW_DATA_COLNUMBERS(MONTH,18)
                  DO N = 1, O3MAPPING_N_REFCOLS
                     COL=TV8O3NEW_DATA_COLUMNS(MONTH,18,N)
                     M=1
                     DO WHILE(TV8O3NEW_DATA_COLUMNS(MONTH,17,M)/=COL &
                          .AND. &
                          M .LT. TV8O3NEW_DATA_COLNUMBERS(MONTH,17))
                        M = M+1
                     ENDDO

                     DO K=1,11
                        IF (TV8O3NEW_DATA_COLUMNS(MONTH,17,M)==COL) THEN
                           SLOPE=(TV8O3NEW_DATA_PROFILES(MONTH,18,K,N) - &
                                TV8O3NEW_DATA_PROFILES(MONTH,17,K,M))/10
                           O3(K,N,TIME)= &
                                TV8O3NEW_DATA_PROFILES(MONTH,18,K,N)+ &
                                SLOPE*DIST
                        ELSE
                           O3(K,N,TIME)= &
                                TV8O3NEW_DATA_PROFILES(MONTH,18,K,N)
                        ENDIF
                     ENDDO
                  ENDDO
               ELSE             !IF LATITUDE
                  ZONE = 0
                  LOOP = .TRUE.
                  DO WHILE (LOOP)
                     ZONE = ZONE + 1
                     IF ( LATITUDE .GT. LATZONES(ZONE)+5 .AND. &
                          LATITUDE .LE. LATZONES(ZONE+1)+5 ) &
                          LOOP = .FALSE.
                  ENDDO
                  DIST=LATITUDE-(LATZONES(ZONE)+5)

                  IF( TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE) >= &
                       TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1)) THEN
                     REF=ZONE
                     O3MAPPING_N_REFCOLS &
                          =TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE)

                     DO N = 1, O3MAPPING_N_REFCOLS
                        COL=TV8O3NEW_DATA_COLUMNS(MONTH,REF,N)
                        M=1
                        DO WHILE( &
                             TV8O3NEW_DATA_COLUMNS(MONTH,REF+1,M)/=COL &
                             .AND. M .LT. &
                             TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1) &
                             )
                           M = M+1
                        ENDDO

                        DO K=1,11
                           IF (TV8O3NEW_DATA_COLUMNS(MONTH,REF+1,M) &
                                ==COL) THEN
                              SLOPE= &
                              (TV8O3NEW_DATA_PROFILES(MONTH,REF+1,K,M) - &
                              TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N))/10
                              O3(K,N,TIME)= &
                                   TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N) &
                                   + SLOPE*DIST
                           ELSE
                              O3(K,N,TIME)= &
                                   TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N)
                           ENDIF
                        ENDDO
                     ENDDO
                  ELSE
                     REF=ZONE+1
                     O3MAPPING_N_REFCOLS= &
                          TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE+1)

                     DO N = 1, O3MAPPING_N_REFCOLS
                        COL=TV8O3NEW_DATA_COLUMNS(MONTH,REF,N)
                        M=1
                        DO WHILE( &
                             TV8O3NEW_DATA_COLUMNS(MONTH,REF-1,M)/=COL &
                             .AND. M .LT. &
                             TV8O3NEW_DATA_COLNUMBERS(MONTH,ZONE-1) &
                             )
                           M = M+1
                        ENDDO

                        DO K=1,11
                           IF (TV8O3NEW_DATA_COLUMNS(MONTH,REF-1,M) &
                                ==COL) THEN
                              SLOPE= &
                                   (TV8O3NEW_DATA_PROFILES &
                                   (MONTH,REF,K,N) - &
                                   TV8O3NEW_DATA_PROFILES &
                                   (MONTH,REF-1,K,M))/10
                              O3(K,N,TIME)= &
                                   TV8O3NEW_DATA_PROFILES &
                                   (MONTH,REF-1,K,M) &
                                   + SLOPE*DIST
                           ELSE
                              O3(K,N,TIME)= &
                                   TV8O3NEW_DATA_PROFILES(MONTH,REF,K,N)
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
                             (O3(K,N,2)-O3(K,N,1))/NBDAY(12)* &
                             (NBDAY(12)-15+DAY_OF_MONTH)
                     ELSE
                        TV8O3NEW_DATA_PROFILES(MONTH,1,K,N) = O3(K,N,1)+ &
                             (O3(K,N,2)-O3(K,N,1))/NBDAY(MONTH-1)* &
                             (NBDAY(MONTH-1)-15+DAY_OF_MONTH)
                     ENDIF
                  ELSE
                     TV8O3NEW_DATA_PROFILES(MONTH,1,K,N) = O3(K,N,1)+ &
                          (O3(K,N,2)-O3(K,N,1))/ &
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

!  set Log pressure grid

      DO K = 1, 14
        CUMO3_LNPRESS(K)     = DLOG(DBLE(TV8O3NEW_DATA_PRESSURES(K)))
      ENDDO
      CUMO3_LNPRESS(14) = DLOG(SURFACE_PRESS)


!  Now interpolate to our incoming column

      CALL LINEAR_TRAWLER_ASCENDING &
            ( 10, O3MAPPING_N_REFCOLS, &
              O3MAPPING_REFCOLS, ACTUAL_O3COLUMN, &
              DF1, DF2, N1, N2 )

!  interpolate Cumulative O3 to incoming pressure grid

      nlevels = nlayers + 1
      CALL LINTP2 &
         ( 14,CUMO3_LNPRESS,CUMO3_MAPPING(1,N1), &
           NLEVELS, LEVEL_LOGPRESSURES, LOCAL_O31 )
      CALL LINTP2 &
         ( 14,CUMO3_LNPRESS,CUMO3_MAPPING(1,N2), &
           NLEVELS, LEVEL_LOGPRESSURES, LOCAL_O32 )

      DIFF = O3MAPPING_REFCOLS(N2) - O3MAPPING_REFCOLS(N1)
      XFAC = ACTUAL_O3COLUMN - O3MAPPING_REFCOLS(N1)
      SUM = 0.0D0
      DO N = 1, NLAYERS
        IF ( N.EQ.1 ) THEN
          CGRAD =  ( LOCAL_O32(2) - LOCAL_O31(2) ) / DIFF
          O3UMKEHRS(N)    = LOCAL_O31(2) + XFAC * CGRAD
          DC_O3UMKEHRS(N) = CGRAD
          CGRADP = CGRAD
        ELSE
          CGRAD = ( LOCAL_O32(N+1) - LOCAL_O31(N+1) ) / DIFF
          GRADN = CGRAD - CGRADP
          LOCAL = LOCAL_O31(N+1) - LOCAL_O31(N)
          O3UMKEHRS(N)    =  LOCAL + XFAC * GRADN
          DC_O3UMKEHRS(N) =  GRADN
          CGRADP = CGRAD
        ENDIF
        SUM = SUM + O3UMKEHRS(N)
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE GET_TV8O3_GM1_All

!  Get Ozone Xsec Data: GOME (Part 1)

      SUBROUTINE GET_O3DATA_GOME_Part1 &
        ( FILENAME,                                    & !Inputs
          O3XS_MAX_DATAPOINTS, O3XS_MAX_TEMPS,         & !Inputs
          O3XS_N_DATAPOINTS, O3XS_N_TEMPS,             & !Outputs
          O3XS_TEMPDATA, O3XS_WAVEDATA, O3XS_XSECDATA, & !Outputs
          FAIL, MESSAGE )                                !Outputs

!mick mod 7/20/2016
!  Note: Subroutine GET_O3XSEC_GOME was split into two for this module
!        (a file read portion and a XSEC calculation portion);
!        GET_O3DATA_GOME_Part1 is the file read portion

      implicit none
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  File detail

      CHARACTER (LEN=*), INTENT(IN)  :: FILENAME

!  Dimensioning

      INTEGER, INTENT(IN)    :: O3XS_MAX_DATAPOINTS, O3XS_MAX_TEMPS

!  Output arguments
!  ----------------

!  Dimensioning

      INTEGER, INTENT(OUT)   :: O3XS_N_DATAPOINTS, O3XS_N_TEMPS

!  Cross section output
                                
      REAL(FPK), INTENT(OUT) :: O3XS_TEMPDATA ( O3XS_MAX_TEMPS ), &
                                O3XS_WAVEDATA ( O3XS_MAX_DATAPOINTS, O3XS_MAX_TEMPS ), &
                                O3XS_XSECDATA ( O3XS_MAX_DATAPOINTS, O3XS_MAX_TEMPS )

!  Status

      LOGICAL, INTENT(OUT)   :: FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE

!  Local variables
!  ===============

!  Help variables

      INTEGER ::   K, M, NP

!  Initialise file read

      FAIL = .FALSE.

!  Open file

      OPEN ( UNIT   = 55, &
             FILE   = Adjustl(Trim(FILENAME)), &
             ERR    = 900, &
             STATUS = 'OLD' )

!  Read and check temperature data first

      READ(55,*)O3XS_N_TEMPS
      IF ( O3XS_N_TEMPS .GT. O3XS_MAX_TEMPS ) GOTO 902
      DO K = 1, O3XS_N_TEMPS
        READ(55,*)O3XS_TEMPDATA(K)
      ENDDO

!  Main read section for basic data file

      DO K = 1, O3XS_N_TEMPS
        READ(55,*)NP
        IF ( NP .GT. O3XS_MAX_DATAPOINTS ) GOTO 902
        DO M = 1, NP
          READ(55,*) O3XS_WAVEDATA(M,K),O3XS_XSECDATA(M,K)
        ENDDO
      ENDDO
      O3XS_N_DATAPOINTS = NP

!  Close file

      CLOSE ( 55 )

!  Finish

      RETURN

!  Open-file error return

900   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Could not find file of O3 cross-sections to open'
      CLOSE(55)
      RETURN

!  Dimensioning error return

902   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Dimensioning insufficient for O3 cross-section data'
      CLOSE(55)
      RETURN

      END SUBROUTINE GET_O3DATA_GOME_Part1

!  Define Ozone X-Sec: GOME (Part 2)

      SUBROUTINE GET_O3XSEC_GOME_Part2 &
        ( O3XS_MAX_DATAPOINTS, O3XS_MAX_TEMPS,         & !Inputs
          O3XS_N_DATAPOINTS, O3XS_N_TEMPS,             & !Inputs
          O3XS_TEMPDATA, O3XS_WAVEDATA, O3XS_XSECDATA, & !Inputs
          MAXLAMBDAS, MAX_LAYERS, NLAMBDAS, NLAYERS,   & !Inputs
          LAMBDAS, TEMPERATURES,                       & !Inputs
          O3_CROSSSECS )                                 !Outputs

!mick mod 7/20/2016
!  Note: Subroutine GET_O3XSEC_GOME was split into two for this module
!        (a file read portion and a XSEC calculation portion);
!        GET_O3XSEC_GOME_Part2 is the XSEC calculation portion

      implicit none
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  Dimensioning

      INTEGER, INTENT(IN) ::    O3XS_MAX_DATAPOINTS, O3XS_MAX_TEMPS

!  Input cross section data

      INTEGER, INTENT(IN) ::    O3XS_N_DATAPOINTS, O3XS_N_TEMPS 
      REAL(FPK), INTENT(IN) ::  O3XS_TEMPDATA ( O3XS_MAX_TEMPS), &
                                O3XS_WAVEDATA ( O3XS_MAX_DATAPOINTS, O3XS_MAX_TEMPS ), &
                                O3XS_XSECDATA ( O3XS_MAX_DATAPOINTS, O3XS_MAX_TEMPS )

!  Dimensioning

      INTEGER, INTENT(IN) ::    MAXLAMBDAS, MAX_LAYERS

!  Input wavelengths

      INTEGER, INTENT(IN) ::    NLAMBDAS
      REAL(FPK), INTENT(IN) ::  LAMBDAS ( MAXLAMBDAS )

!  Input temperatures

      INTEGER, INTENT(IN) ::    NLAYERS
      REAL(FPK), INTENT(IN) ::  TEMPERATURES ( MAX_LAYERS )

!  Output arguments
!  ----------------

!  Cross section output

      REAL(FPK), INTENT(OUT) :: O3_CROSSSECS ( MAXLAMBDAS, MAX_LAYERS )

!  Local variables
!  ===============

!  Cross section data
!  ------------------

      REAL(FPK) ::   O3XS_XSEC2 ( O3XS_MAX_DATAPOINTS )

!  Help variables

      INTEGER ::   W, N, K, K1, K2, NP
      REAL(FPK) :: DF1, DF2

!  Variables for splined interpolation

      REAL(FPK) :: XSEC ( O3XS_MAX_DATAPOINTS, O3XS_MAX_TEMPS )

!  Start data loop

      NP = O3XS_N_DATAPOINTS
      DO K = 1, O3XS_N_TEMPS

!  Spline interpolate immediately

        CALL SPLINE ( O3XS_WAVEDATA(:,K), O3XS_XSECDATA(:,K), &
                      NP, 0.0D0, 0.0D0, O3XS_XSEC2 )

        DO W = 1, NLAMBDAS
          IF ( LAMBDAS(W) .LE. O3XS_WAVEDATA(NP,K) .AND. &
               LAMBDAS(W) .GE. O3XS_WAVEDATA(1,K) ) THEN
            CALL SPLINT ( O3XS_WAVEDATA(:,K), O3XS_XSECDATA(:,K), &
                          O3XS_XSEC2, NP, LAMBDAS(W), XSEC(W,K) )
          ELSE
            XSEC(W,K) = 0.0d0
          ENDIF
        ENDDO

!  Finish data loop

      ENDDO

!  Debug
!      do w = 1, nlambdas
!        write(*,'(f8.2,1p5e15.5)')lambdas(w),(xsec(w,k1),k1 = 1,5)
!      enddo

!  Set the cross-sections
!  ======================

!  Temperature interpolate (linear) clear sky cross sections

      DO N = 1, NLAYERS
        CALL LINEAR_TRAWLER_ASCENDING &
            ( O3XS_MAX_TEMPS, O3XS_N_TEMPS, &
              O3XS_TEMPDATA,  TEMPERATURES(N), &
              DF1, DF2, K1, K2 )
        DO W = 1, NLAMBDAS
          O3_CROSSSECS(W,N) = DF1 * XSEC(W,K2) + DF2 * XSEC(W,K1)
        ENDDO
      ENDDO

!  Finish

      RETURN

      END SUBROUTINE GET_O3XSEC_GOME_Part2

!  Get Ozone Xsecs GOME (no partition, used in the Bin setups)

      SUBROUTINE GET_O3XSEC_GOME_All &
           ( FILENAME, &
             MAXLAMBDAS, MAXLAYERS, NLAMBDAS, NLAYERS, &
             LAMBDAS, TEMPERATURES, &
             O3_CROSSSECS, FAIL, MESSAGE )

      implicit none
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  file detail

      CHARACTER (LEN=*), INTENT(IN) ::    FILENAME

!  Dimensioning

      INTEGER, INTENT(IN) ::          MAXLAMBDAS, MAXLAYERS

!  input wavelengths

      INTEGER, INTENT(IN) ::          NLAMBDAS
      REAL(FPK), INTENT(IN) :: LAMBDAS ( MAXLAMBDAS )

!  input temperatures

      INTEGER, INTENT(IN) ::          NLAYERS
      REAL(FPK), INTENT(IN) :: TEMPERATURES ( MAXLAYERS )

!  Output arguments
!  ----------------

!  Cross section output

      REAL(FPK), INTENT(OUT) :: O3_CROSSSECS ( MAXLAMBDAS, MAXLAYERS )

!  status

      LOGICAL, INTENT(OUT) ::          FAIL
      CHARACTER (LEN=*), INTENT(OUT) ::    MESSAGE

!  Local variables
!  ===============

!  Dimensioning

      INTEGER, PARAMETER :: O3XS_MAX_DATAPOINTS = 1000
      INTEGER, PARAMETER :: O3XS_MAX_TEMPS      = 5

!  Cross section data
!  ------------------

      INTEGER ::          O3XS_N_TEMPS
      REAL(FPK) :: O3XS_WAVEDATA (O3XS_MAX_DATAPOINTS), &
                       O3XS_XSECDATA (O3XS_MAX_DATAPOINTS), &
                       O3XS_XSEC2    (O3XS_MAX_DATAPOINTS), &
                       O3XS_TEMPDATA (O3XS_MAX_TEMPS)

!  help variables

      INTEGER   :: W, N, K, K1, K2, M, NP
      REAL(FPK) :: DF1,DF2

!  variables for splined interpolation

      REAL(FPK) :: XSEC(O3XS_MAX_DATAPOINTS,O3XS_MAX_TEMPS)

!  initialise file read

      FAIL = .FALSE.

!  Open file

      OPEN ( UNIT   = 55, &
             FILE   = Adjustl(Trim(FILENAME)), &
             ERR    = 900, &
             STATUS = 'OLD' )

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

        CALL SPLINE ( O3XS_WAVEDATA,O3XS_XSECDATA, &
                        NP, 0.0D0,0.0D0, O3XS_XSEC2 )

          DO W = 1, NLAMBDAS
            IF ( LAMBDAS(W).LE. O3XS_WAVEDATA(NP) .AND. &
                 LAMBDAS(W).GE. O3XS_WAVEDATA(1) ) THEN
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
        CALL LINEAR_TRAWLER_ASCENDING &
            ( O3XS_MAX_TEMPS, O3XS_N_TEMPS, &
              O3XS_TEMPDATA,  temperatures(n), &
              DF1, DF2, K1, K2 )
        DO W = 1, NLAMBDAS
          O3_CROSSSECS(W,N) = DF1 * XSEC(W,K2) + DF2 * XSEC(W,K1)
        ENDDO
      ENDDO

!  Finish

      RETURN

!  Open-file error return

900   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Could not find file of GOME O3 cross-sections to open'
      CLOSE(55)
      RETURN

!  dimensioning error return

902   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Dimensioning insufficient for GOME O3 X-section data'
      CLOSE(55)
      RETURN

      END SUBROUTINE GET_O3XSEC_GOME_All

!  Get the Breon BDM O3 XSec Data (Part 1)

      SUBROUTINE BREON_O3DATA_Part1 ( &
        MAX_O3ABSDATA, LAMBDA_START, LAMBDA_FINISH, & !Input
        N_O3ABSDATA, O3ABSDATA, FAIL, MESSAGE )          !Output

      implicit none
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input
!  =====

      INTEGER, INTENT(IN)    :: MAX_O3ABSDATA
      REAL(FPK), INTENT(IN)  :: LAMBDA_START, LAMBDA_FINISH

!  Output
!  ======

      INTEGER, INTENT(OUT)   :: N_O3ABSDATA
      REAL(FPK), INTENT(OUT) :: O3ABSDATA (MAX_O3ABSDATA, 4)

      LOGICAL, INTENT(OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE

!  Local
!  =====

      LOGICAL ::   reading
      INTEGER ::   istat, i
      REAL(FPK) :: c0, c1, c2, incr, lamb1, lamb2, wav, wav1, wav2 

      CHARACTER (LEN=90) :: xsec_file_name

!  Initialize status

      fail = .false.
      message = ' '

!  Buffer control, initialization

      i = 0
      lamb1 = LAMBDA_START  - 0.05d0
      lamb2 = LAMBDA_FINISH + 0.05d0

!  Read the data, checking that the desired window is covered

      xsec_file_name = 'lrrs_test/physics_data/XSECS/o3abs_breon_195_660_vacfinal.dat'
      open(unit = 1, file = trim(xsec_file_name), iostat = istat, status = 'old')
      if (istat /= 0) then
        fail    = .true.
        message = 'Open failure for O3 Xsec file'
        return
      endif
!  .. Read header line
      read(1,*)
!  .. Read 1st data line & check that window is covered at the lower end
      read(1,*) wav, c0, c1, c2
      if ( wav .gt. lamb1 ) then
         fail = .true.
         message = 'O3 Xsec data does not cover input window at the lower end'
         return
      elseif ( wav .eq. lamb1 ) then
         i = 1
         O3ABSDATA(i,1) = wav
         O3ABSDATA(i,2) = c0
         O3ABSDATA(i,3) = c1
         O3ABSDATA(i,4) = c2
      endif
!  .. Read more lines
      reading = .true.
      do while (reading)
         read(unit = 1,fmt = *, iostat = istat) wav, c0, c1, c2
         if (istat == -1) then
            wav1 = O3ABSDATA(i-1,1)
            wav2 = O3ABSDATA(i,1)
            incr = wav2 - wav1
            if ( wav2 + incr .gt. lamb2 ) then
               reading = .false.
               N_O3ABSDATA = i
            else
               fail = .true.
               message = 'End of O3 file: O3 Xsec data does not cover input window at the upper end'
               return
            endif
         elseif ( wav .ge. lamb1 ) then
            if ( wav .le. lamb2 ) then
               i = i + 1
               if ( i .gt. MAX_O3ABSDATA ) then
                  fail = .true.
                  message = 'Amount of O3 data too much for O3 data arrays: increase size of MAX_O3ABSDATA ' // &
                            'in subroutine LRRS_INSETUP_LPMONO_M'
                  return
               else
                  O3ABSDATA(i,1) = wav
                  O3ABSDATA(i,2) = c0
                  O3ABSDATA(i,3) = c1
                  O3ABSDATA(i,4) = c2
               endif
            else
              reading = .false.
              N_O3ABSDATA = i
            endif
         endif
      enddo
      close(1)

!  Normal finish

      END SUBROUTINE BREON_O3DATA_Part1

!  Define the Breon BDM O3 XSecs (Part 2)

      SUBROUTINE BREON_O3XSEC_Part2 ( &
        MAX_O3ABSDATA, N_O3ABSDATA, O3ABSDATA,  & !Input
        MAX_POINTS, NLAMBDAS, LAMBDAS,          & !Input
        O3C0, O3C1, O3C2 )                        !Output

      implicit none
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input
!  =====

      INTEGER, INTENT(IN)   :: MAX_O3ABSDATA, N_O3ABSDATA
      REAL(FPK), INTENT(IN) :: O3ABSDATA (MAX_O3ABSDATA, 4)
      INTEGER, INTENT(IN)   :: MAX_POINTS, NLAMBDAS
      REAL(FPK), INTENT(IN) :: LAMBDAS ( MAX_POINTS )

!  Output
!  ======

      REAL(FPK), INTENT(OUT) :: O3C0 ( MAX_POINTS ), &
                                O3C1 ( MAX_POINTS ), &
                                O3C2 ( MAX_POINTS )

!  Local
!  =====

      INTEGER   :: n, nbuff
      REAL(FPK) :: conversion, xsec
      REAL(FPK), DIMENSION(N_O3ABSDATA) :: x, y, y2, yc1, yc2

!  Buffer control, initialization

      conversion = 1.0d+20

      nbuff = N_O3ABSDATA

      x   = O3ABSDATA(1:N_O3ABSDATA,1)
      y   = O3ABSDATA(1:N_O3ABSDATA,2)
      yc1 = O3ABSDATA(1:N_O3ABSDATA,3)
      yc2 = O3ABSDATA(1:N_O3ABSDATA,4)

!  Spline the data

      call spline(x,y,nbuff,0.0d0,0.0d0,y2)
      do n = 1, nlambdas
         call splint(x,y,y2,nbuff,lambdas(n),xsec)
         O3c0(n) = xsec / conversion
      enddo

      call spline(x,yc1,nbuff,0.0d0,0.0d0,y2)
      do n = 1, nlambdas
        call splint(x,yc1,y2,nbuff,lambdas(n),xsec)
        o3c1(n) = xsec / conversion
      enddo

      call spline(x,yc2,nbuff,0.0d0,0.0d0,y2)
      do n = 1, nlambdas
        call splint(x,yc2,y2,nbuff,lambdas(n),xsec)
        o3c2(n) = xsec / conversion
      enddo

      END SUBROUTINE BREON_O3XSEC_Part2

!  Subroutine BREON_O3XSEC_All. Used in LCSbin, LPbin
!     Used to crosscheck BREON_O3XSEC_Part1 and BREON_O3XSEC_Part2

      SUBROUTINE BREON_O3XSEC_All ( &
        MAX_POINTS, NLAMBDAS, LAMBDAS, &
        O3C0, O3C1, O3C2, &
        FAIL, MESSAGE )

      implicit none
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input
!  =====

      INTEGER, INTENT(IN)   :: MAX_POINTS, NLAMBDAS
      REAL(FPK), INTENT(IN) :: LAMBDAS ( MAX_POINTS )

!  Output
!  ======

      REAL(FPK), INTENT(OUT) :: O3C0 ( MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: O3C1 ( MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: O3C2 ( MAX_POINTS )

      LOGICAL, INTENT(OUT) ::           FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE

!  Local
!  =====

      LOGICAL ::   reading
      INTEGER ::   nbuff, ndata, maxdata, nf, n

      INTEGER, PARAMETER :: maxspline = 10000

      REAL(FPK) :: conversion, lamb1, lamb2, wav, val, c1, c2
      REAL(FPK) :: x(maxspline), y(maxspline), y2(maxspline)
      REAL(FPK) :: yc1(maxspline), yc2(maxspline), xsec

      CHARACTER (LEN=90) :: filename, xsec_file_name

!  Initialize status

      fail = .false.
      message = ' '

!  Buffer control, initialization

      nbuff = 0
      ndata = 1
      lamb1 = lambdas(1)        - 0.05d0
      lamb2 = lambdas(nlambdas) + 0.05d0
      maxdata    = 46501
      conversion = 1.0d+20

!  Prepare filename

      xsec_file_name = &
          'lrrs_test/physics_data/XSECS/o3abs_breon_195_660_vacfinal.dat'
      nf = 1
      do while (xsec_file_name(nf:nf).ne. ' ')
        nf = nf + 1
      enddo
      nf = nf - 1
      filename = xsec_file_name(1:nf)

!  Read the data, check that window is covered

      open(1, file = filename, err = 90, status='old')
!  .. First line is a dummy
      read(1,*)
!  .. Read second line, and check window is covered
      read(1,*)wav, val, c1, c2
      if (lamb1 .lt. wav ) then
         fail = .true.
         message = 'O3 Xsec data no cover input window at lower end'
         go to 91
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
         message = 'O3 Xsec data not cover input window at upper end'
         go to 91
      endif
      close(1)

      if (nbuff .gt. maxspline) then
         fail = .true.
         message = 'need to increase maxspline '
         go to 91
      endif

!  Spline the data

      call spline(x,y,nbuff,0.0d0,0.0d0,y2)
      do n = 1, nlambdas
         call splint(x,y,y2,nbuff,lambdas(n),xsec)
         O3c0(n) = xsec / conversion
      enddo

      call spline(x,yc1,nbuff,0.0d0,0.0d0,y2)
      do n = 1, nlambdas
        call splint(x,yc1,y2,nbuff,lambdas(n),xsec)
        o3c1(n) = xsec / conversion
      enddo

      call spline(x,yc2,nbuff,0.0d0,0.0d0,y2)
      do n = 1, nlambdas
        call splint(x,yc2,y2,nbuff,lambdas(n),xsec)
        o3c2(n) = xsec / conversion
      enddo

!  Normal finish

      return

!  error returns

90    continue
      fail    = .true.
      message = 'Open failure for Xsec file '
      return

91    continue
      return

      END SUBROUTINE BREON_O3XSEC_All

!  Get Ozone X-Secs: Other sources

      SUBROUTINE GET_O3XSEC_1 &
        ( FILENAME,                                  & !Inputs
          MAXLAMBDAS, MAX_LAYERS, NLAMBDAS, NLAYERS, & !Inputs
          LAMBDAS, TEMPERATURES,                     & !Inputs
          O3_CROSSSECS, FAIL, MESSAGE )                !Outputs

      implicit none
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  File detail

      CHARACTER (LEN=*), INTENT(IN) ::    FILENAME

!  Dimensioning

      INTEGER, INTENT(IN) ::   MAXLAMBDAS, MAX_LAYERS

!  Input wavelengths

      INTEGER, INTENT(IN) ::   NLAMBDAS
      REAL(FPK), INTENT(IN) :: LAMBDAS ( MAXLAMBDAS )

!  Input temperatures

      INTEGER, INTENT(IN) ::   NLAYERS
      REAL(FPK), INTENT(IN) :: TEMPERATURES ( MAX_LAYERS )

!  Output arguments
!  ----------------

!  Cross section output

      REAL(FPK), INTENT(OUT) :: O3_CROSSSECS ( MAXLAMBDAS, MAX_LAYERS )

!  Status

      LOGICAL, INTENT(OUT)           :: FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE

!  Local variables
!  ===============

!  Dimensioning

      INTEGER, PARAMETER :: O3XS_MAX_DATAPOINTS = 1000
      INTEGER, PARAMETER :: O3XS_MAX_TEMPS      = 5

!  Cross section data
!  ------------------

      INTEGER ::   O3XS_N_TEMPS
      REAL(FPK) :: O3XS_WAVEDATA (O3XS_MAX_DATAPOINTS), &
                   O3XS_XSECDATA (O3XS_MAX_DATAPOINTS), &
                   O3XS_XSEC2    (O3XS_MAX_DATAPOINTS), &
                   O3XS_TEMPDATA (O3XS_MAX_TEMPS)

!  Help variables

      INTEGER ::   W, N, K, K1, K2,M, NP
      REAL(FPK) :: DF1,DF2

!  Variables for splined interpolation

      REAL(FPK) :: XSEC(O3XS_MAX_DATAPOINTS,O3XS_MAX_TEMPS)

!  Initialise file read

      FAIL = .FALSE.

!  Open file

      OPEN ( UNIT   = 55, &
             FILE   = Adjustl(Trim(FILENAME)), &
             ERR    = 900, &
             STATUS = 'OLD' )

!  Read and check temperature data first

      READ(55,*)O3XS_N_TEMPS
      IF ( O3XS_N_TEMPS .GT. O3XS_MAX_TEMPS ) GOTO 902
      DO K = 1, O3XS_N_TEMPS
        READ(55,*)O3XS_TEMPDATA(K)
      ENDDO

!  Start Data loop

      DO K = 1, O3XS_N_TEMPS

!  Main read section for basic data file

        READ(55,*)NP
        IF ( NP .GT. O3XS_MAX_DATAPOINTS ) GOTO 902
        DO M = 1, NP
          READ(55,*) O3XS_WAVEDATA(M),O3XS_XSECDATA(M)
        ENDDO

!  Spline interpolate immediately

        CALL SPLINE ( O3XS_WAVEDATA,O3XS_XSECDATA, &
                      NP, 0.0D0,0.0D0, O3XS_XSEC2 )

        DO W = 1, NLAMBDAS
          IF ( LAMBDAS(W).LE. O3XS_WAVEDATA(NP) .AND. &
               LAMBDAS(W).GE. O3XS_WAVEDATA(1) ) THEN
            CALL SPLINT ( O3XS_WAVEDATA,O3XS_XSECDATA, &
                          O3XS_XSEC2,NP,LAMBDAS(W),XSEC(W,K))
          ELSE
            XSEC(W,K) = 0.0d0
          ENDIF
        ENDDO

!  Finish data loop

        ENDDO

!  Close file

      CLOSE ( 55 )

!  Debug
!      do w = 1, nlambdas
!        write(*,'(f8.2,1p5e15.5)')lambdas(w),(xsec(w,k1),k1 = 1,5)
!      enddo
!      pause

!  Set the cross-sections
!  ======================

!  Temperature interpolate (linear) clear sky cross sections

      DO N = 1, NLAYERS
        CALL LINEAR_TRAWLER_ASCENDING &
            ( O3XS_MAX_TEMPS, O3XS_N_TEMPS, &
              O3XS_TEMPDATA,  TEMPERATURES(N), &
              DF1, DF2, K1, K2 )
        DO W = 1, NLAMBDAS
          O3_CROSSSECS(W,N) = DF1 * XSEC(W,K2) + DF2 * XSEC(W,K1)
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

!  Dimensioning error return

902   CONTINUE
      FAIL = .TRUE.
      MESSAGE = 'Dimensioning insufficient for O3 cross-section data'
      CLOSE(55)
      RETURN

      END SUBROUTINE GET_O3XSEC_1

!  Get Rayleigh Cross-sections and depolarization ratios

      SUBROUTINE RAYLEIGH_FUNCTION &
        ( FORWARD_MAXLAMBDAS, CO2_PPMV_MIXRATIO, & !Inputs
          FORWARD_NLAMBDAS,   FORWARD_LAMBDAS,   & !Inputs
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL )          !Outputs

!  Rayleigh cross sections and depolarization ratios
!     Bodhaine et. al. (1999) formulae

      IMPLICIT NONE
      !integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input arguments
!  ---------------

!  Wavelength

      INTEGER, INTENT(IN)   :: FORWARD_MAXLAMBDAS, FORWARD_NLAMBDAS
      REAL(FPK), INTENT(IN) :: FORWARD_LAMBDAS ( FORWARD_MAXLAMBDAS )

!  CO2 mixing ratio

      REAL(FPK), INTENT(IN) :: CO2_PPMV_MIXRATIO

!  Output arguments
!  ----------------

!  Cross-sections and depolarization output

      REAL(FPK), INTENT(INOUT) :: RAYLEIGH_XSEC  ( FORWARD_MAXLAMBDAS )
      REAL(FPK), INTENT(INOUT) :: RAYLEIGH_DEPOL ( FORWARD_MAXLAMBDAS )

!  Local variables
!  ---------------

      INTEGER ::          W
      REAL(FPK) :: MASS_DRYAIR
      REAL(FPK) :: NMOL, PI, CONS
      REAL(FPK) :: MO2,MN2,MARG,MCO2,MAIR
      REAL(FPK) :: FO2,FN2,FARG,FCO2,FAIR
      REAL(FPK) :: LAMBDA_C,LAMBDA_M,LPM2,LP2
      REAL(FPK) :: N300M1,NCO2M1,NCO2
      REAL(FPK) :: NCO2SQ, NSQM1,NSQP2,TERM

!  Data and parameters
!  -------------------

      REAL(FPK), PARAMETER :: S0_A = 15.0556D0
      REAL(FPK), PARAMETER :: S0_B = 28.9595D0

      REAL(FPK), PARAMETER :: S1_A = 8060.51D0
      REAL(FPK), PARAMETER :: S1_B = 2.48099D+06
      REAL(FPK), PARAMETER :: S1_C = 132.274D0
      REAL(FPK), PARAMETER :: S1_D = 1.74557D+04
      REAL(FPK), PARAMETER :: S1_E = 39.32957D0

      REAL(FPK), PARAMETER :: S2_A = 0.54D0

      REAL(FPK), PARAMETER :: S3_A = 1.034D0
      REAL(FPK), PARAMETER :: S3_B = 3.17D-04
      REAL(FPK), PARAMETER :: S3_C = 1.096D0
      REAL(FPK), PARAMETER :: S3_D = 1.385D-03
      REAL(FPK), PARAMETER :: S3_E = 1.448D-04

      MO2  = 20.946D0
      MN2  = 78.084D0
      MARG =  0.934D0

!  Start of code
!  -------------

!  Constants

      NMOL = 2.546899D19
      PI   = DATAN(1.0D0)*4.0D0
      CONS = 24.0D0 * PI * PI * PI

!  Convert co2

      MCO2 = 1.0D-06 * CO2_PPMV_MIXRATIO

!  Mass of dry air: Eq.(17) of BWDS

      MASS_DRYAIR = S0_A * MCO2 + S0_B

!  Start loop

      DO W = 1, FORWARD_NLAMBDAS

!  Wavelength in micrometers

      LAMBDA_M = 1.0D-03 * FORWARD_LAMBDAS(W)
      LAMBDA_C = 1.0D-07 * FORWARD_LAMBDAS(W)
      LPM2     = 1.0D0 / LAMBDA_M / LAMBDA_M

!  Step 1: Eq.(18) of BWDS

      N300M1 = S1_A + ( S1_B / ( S1_C - LPM2 ) ) + &
                      ( S1_D / ( S1_E - LPM2 ) )
      N300M1 = N300M1 * 1.0D-08

!  Step 2: Eq.(19) of BWDS

      NCO2M1 = N300M1 * ( 1.0D0 + S2_A * ( MCO2  - 0.0003D0 ) )
      NCO2   = NCO2M1 + 1
      NCO2SQ = NCO2 * NCO2

!  Step 3: Eqs. (5&6) of BWDS (Bates' results)

      FN2  = S3_A + S3_B * LPM2
      FO2  = S3_C + S3_D * LPM2 + S3_E * LPM2 * LPM2

!  Step 4: Eq.(23) of BWDS
!     ---> King factor and depolarization ratio

      FARG = 1.0D0
      FCO2 = 1.15D0
      MAIR = MN2 + MO2 + MARG + MCO2
      FAIR = MN2*FN2 + MO2*FO2 + MARG*FARG + MCO2*FCO2
      FAIR = FAIR / MAIR
      RAYLEIGH_DEPOL(W) = 6.0D0*(FAIR-1.0D0)/(3.0D0+7.0D0*FAIR)

!  Step 5: Eq.(22) of BWDS
!     ---> Cross section

      LP2  = LAMBDA_C * LAMBDA_C
      NSQM1 = NCO2SQ - 1.0D0
      NSQP2 = NCO2SQ + 2.0D0
      TERM = NSQM1 / LP2 / NMOL / NSQP2
      RAYLEIGH_XSEC(W) =  CONS * TERM * TERM * FAIR

!  End loop

      ENDDO

!  Finish
!  ------

      RETURN
      END SUBROUTINE RAYLEIGH_FUNCTION

!  End module

End Module lrrs_iopsetup_WorkRoutines_m
