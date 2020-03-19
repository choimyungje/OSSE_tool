
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

      PROGRAM O3BIN_M

!  Used modules

      USE LRRS_PARS_m
      USE LRRS_INPUTS_DEF_m
      USE LRRS_OUTPUTS_DEF_m
      USE LRRS_SUP_INOUT_DEF_m

      USE LRRS_IO_CHECK_m
      USE O3BIN_M_IOPSETUP_NOLIN_m
      USE LRRS_MAIN_MASTER_m

      USE BRDF_SUP_MOD_m
      USE LRRS_SUP_ACCESSORIES_1_m

      IMPLICIT NONE

!  LIDORT-RRS input & output type structures
!  Type structures, augmented by supplemental input, Version 2.5, 9/23/15

      TYPE(LRRS_Fixed_Inputs)       :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs)    :: LRRS_ModIn
      TYPE(LRRS_Sup_InOut)          :: LRRS_Sup
      TYPE(LRRS_Outputs)            :: LRRS_Out

!  BRDF supplement / LRRS BRDF-related inputs consistency check status

      TYPE(LRRS_Exception_Handling) :: LRRS_BRDFCheck_Status

!  BRDF supplement input structure

      TYPE(BRDF_Sup_Inputs)                 :: BRDF_Sup_In

!  BRDF supplement file inputs status structure

      TYPE(BRDF_Input_Exception_Handling)   :: BRDF_Sup_InputStatus

!  LIDORT-RRS related variables
!  ============================

!  Boolean variables

      LOGICAL ::   DO_RRS_OVERALL
      LOGICAL ::   DO_UPWELLING
      LOGICAL ::   DO_DNWELLING
      LOGICAL ::   DO_ENERGY_BALANCING
      LOGICAL ::   DO_BIN_REALIZATION

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)
!    LBOUNDARIES_OUTPUT = layer boundary choices
!    LPARTIALS_OUTPUT = off-boundary level choices

      INTEGER ::   N_LOUTPUT
      INTEGER ::   LBOUNDARIES_OUTPUT ( MAX_LOUTPUT )

!  BINNING: outer/inner wavelength range, and offset for inner range
!       Depends on the choice of solar spectrum

      INTEGER ::   OFFSET_INNER
      INTEGER ::   NPOINTS_INNER
      !INTEGER ::   NPOINTS_OUTER

!  Wavelengths and Fluxes
!  ----------------------

!  These are defined on the outer grid for the binning realization
!  These are defined on 234 points for the monochromatic

!  Depends on the choice of solar spectrum

      REAL(FPK) :: LAMBDAS_RANKED  ( MAX_POINTS )
      REAL(FPK) :: FLUXES_RANKED   ( MAX_POINTS )
      !REAL(FPK) :: ALBEDOS_RANKED   ( MAX_POINTS )

!  ELASTIC FIELD
!  -------------

!  Fourier-summed output

      REAL(FPK) :: ELASTIC_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: ELASTIC_DN &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  RAMAN FIELD
!  -----------

!  Radiance including inelastic scattering

      REAL(FPK) :: RAMAN_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: RAMAN_DN &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  General error Handling

      INTEGER ::   STATUS_INPUTREAD
      INTEGER ::   N_MESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( MAX_MESSAGES )

!  Helper variables
!  ================

!  Start and Finish wavelengths

      REAL(FPK) :: LAMBDA_START, LAMBDA_FINISH

!  Ozone

      LOGICAL ::   No_O3

!  Spectrum of gas optical depth

      REAL(FPK) :: GASVOD ( MAX_POINTS )

!  Surface

      INTEGER   :: NPOINTS_SURF
      REAL(FPK) :: ALBEDO

!  Filling

      REAL(FPK) :: RRSFILLING_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: RRSFILLING_DN &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Instrument control (1 = GOME, 2 = SCIA, 3 = OMI)

      INTEGER ::   INSTRUMENT_INDEX

!  Slit function options

      LOGICAL ::   DO_SLITFUNC_WEIGHTING
      REAL(FPK) :: SLITFUNC_FWHM, SLITFUNC_CUTOFF

!  Error handling for Physics and IOP setup

      LOGICAL ::   FAIL
      !CHARACTER (LEN=70) :: MESSAGE_SUB

!  Other local variables

      INTEGER ::   ATP, CTP, ISZA, N, N_GEOMETRIES, OD, V, WP
      REAL(FPK) :: I_ELA, I_RRS, RATIO

!  File name control

      CHARACTER (LEN=2)   :: C2
      CHARACTER (LEN=256) :: SOLARFILE
      CHARACTER (LEN=80)  :: FILNAM36, FILNAM37

!  Writing to plot file

      LOGICAL :: WRITE_TO_PLOTFILE

!  Timing tests

      REAL :: E1, E2, E3, E4, StdTime

!  Start program

!  Start timing (off for LRRS2.5a comparison)

      !CALL CPU_TIME(E1)

!  INPUT SECTION
!  =============

!  Start and Finish

!     Will calculate 105 points (for OMP testing)
!      LAMBDA_START  = 324.5
!      LAMBDA_FINISH = 336.5

!     Will calculate 52 points (for OMP testing)
      ! LAMBDA_START  = 324.5
      ! LAMBDA_FINISH = 330.5
      ! LAMBDA_START  = 750.0
      ! LAMBDA_FINISH = 753.0
!     Will calculate 8 points
!      LAMBDA_START  = 324.5
!      LAMBDA_FINISH = 325.5
      !----from min/max solar file; option 1 GOME
      ! LAMBDA_START  = 312.3
      ! LAMBDA_FINISH = 405.0

      ! LAMBDA_START  = 309.978+16.0
      ! LAMBDA_FINISH = 377.527-16.0
      
!     Will calculate 2 points (for debugging)
     LAMBDA_START  = 324.5
     LAMBDA_FINISH = 324.8

!     Will calculate 1 point (for debugging)
!      LAMBDA_START  = 324.5
!      LAMBDA_FINISH = 324.65

!  Instrument control (1 = GOME, 2 = SCIA, 3 = OMI)

      INSTRUMENT_INDEX = 1

!  Read status of ozone

      ! No_O3 = .false.
     No_O3 = .true.

!  Surface
!  Note: NPOINTS_SURF should be = 1 or NPOINTS_OUTER

      NPOINTS_SURF = 1
      ALBEDO = 0.05D0

!  Extra output

      !WRITE_TO_PLOTFILE = .TRUE.
      WRITE_TO_PLOTFILE = .FALSE.

!  Input data for LRRS code control (USER-DEFINED INPUT)

      CALL LRRS_READ_INPUTS ( &
        'lrrs_test/lrrs_o3bin_M_tester.cfg', &
         LRRS_FixIn, LRRS_ModIn, &
         STATUS_INPUTREAD, N_MESSAGES, MESSAGES )

!  Exception handling. Go to the finish and write errors.

      IF ( STATUS_INPUTREAD .NE. LRRS_SUCCESS ) THEN
        WRITE(*,*)
        WRITE(*,*) 'FAIL after LRRS_READ_INPUTS'
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = 'input read failed - stop'
        GO TO 5678
      ENDIF

!  Define some local variables

      !If DO_RRS_OVERALL is set, then model will calculate Intensity with RRS
      !If NOT set, then the model will do an Elastic calculation only
      DO_RRS_OVERALL  = LRRS_FixIn%Bool%DO_RRS_OVERALL
      !DO_ELASTIC_ONLY = LRRS_ModIn%MBool%DO_ELASTIC_ONLY

      DO_UPWELLING    = LRRS_FixIn%Bool%DO_UPWELLING
      DO_DNWELLING    = LRRS_FixIn%Bool%DO_DNWELLING

      !These two options are mutually exclusive

      DO_ENERGY_BALANCING = LRRS_ModIn%MBool%DO_ENERGY_BALANCING
      !DO_CABANNES_RAMAN   = LRRS_ModIn%MBool%DO_CABANNES_RAMAN

      DO_BIN_REALIZATION  = LRRS_ModIn%MBool%DO_BIN_REALIZATION

      N_GEOMETRIES        = LRRS_FixIn%UserVal%N_USER_RELAZMS * &
                            LRRS_FixIn%UserVal%N_USER_STREAMS
      N_LOUTPUT           = LRRS_FixIn%UserVal%N_LOUTPUT
      LBOUNDARIES_OUTPUT  = LRRS_FixIn%UserVal%LBOUNDARIES_OUTPUT

! !  This is Binning, no MONO

      IF ( .not. DO_BIN_REALIZATION ) &
        STOP 'This is BINNING - set flag!!'

!  Initialize LRRS standard supplement inputs

      CALL LRRS_Sup_Init ( LRRS_Sup )

!  BRDF input control read

      IF ( LRRS_FixIn%Bool%DO_BRDF_SURFACE ) THEN
        CALL BRDF_INPUTMASTER ( &
          'lrrs_test/lrrs_BRDF_ReadInput_LRRSVal.cfg', & ! Input
          BRDF_Sup_In,         & ! Outputs
          BRDF_Sup_InputStatus ) ! Outputs

        !Exception handling
        IF ( BRDF_Sup_InputStatus%BS_STATUS_INPUTREAD .NE. LRRS_SUCCESS ) then
          OPEN(1,file = 'lrrs_BRDF_ReadInput.log', status = 'unknown')
          WRITE(1,*)' FATAL:   Wrong input from BRDF input file-read'
          WRITE(1,*)'  ------ Here are the messages and actions '
          WRITE(1,'(A,I3)')'    ** Number of messages = ',&
            BRDF_Sup_InputStatus%BS_NINPUTMESSAGES
          DO N = 1, BRDF_Sup_InputStatus%BS_NINPUTMESSAGES
            WRITE(1,'(A,I3,A,A)')'Message # ',N,' : ',&
              adjustl(trim(BRDF_Sup_InputStatus%BS_INPUTMESSAGES(N)))
            WRITE(1,'(A,I3,A,A)')'Action  # ',N,' : ',&
              adjustl(trim(BRDF_Sup_InputStatus%BS_INPUTACTIONS(N)))
          ENDDO
          CLOSE(1)
          STOP 'Read-input fail: Look at file lrrs_BRDF_ReadInput.log'
        ENDIF

        !Check compatibility of BRDF and Main LRRS inputs
        CALL LRRS_BRDF_INPUT_CHECKER ( &
          BRDF_Sup_In,          & ! Inputs
          LRRS_FixIn,           & ! Inputs
          LRRS_ModIn,           & ! Inputs
          LRRS_BRDFCheck_Status ) ! Outputs

        !Exception handling
        IF ( LRRS_BRDFCheck_Status%STATUS_INPUTCHECK .NE. LRRS_SUCCESS ) THEN
          OPEN(1,file = 'lrrs_BRDFcheck.log',status = 'unknown')
          WRITE(1,*)' FATAL: Main and BRDFSup inputs are incompatible'
          WRITE(1,*)'  ------ Here are the messages and actions '
          WRITE(1,'(A,I3)')'    ** Number of messages = ',LRRS_BRDFCheck_Status%NCHECKMESSAGES
          DO N = 1, LRRS_BRDFCheck_Status%NCHECKMESSAGES
            WRITE(1,'(A,I3,A,A)')'Message # ',N,': ',ADJUSTL(TRIM(LRRS_BRDFCheck_Status%CHECKMESSAGES(N)))
            WRITE(1,'(A,I3,A,A)')'Action  # ',N,': ',ADJUSTL(TRIM(LRRS_BRDFCheck_Status%ACTIONS(N)))
          ENDDO
          CLOSE(1)
          STOP 'Checking fail: Look at file lrrs_BRDFcheck.log'
        ENDIF
      ENDIF

!  Obtain input data of solar wavelengths and fluxes (USER FILE-READ)

      IF ( INSTRUMENT_INDEX.EQ.1 ) THEN
        SOLARFILE = 'lrrs_test/physics_data/SOLAR/gome_17296_solarspectrum.txt'
      ELSE IF ( INSTRUMENT_INDEX.EQ.2 ) THEN
        SOLARFILE = 'lrrs_test/physics_data/SOLAR/sciasun_ch2_o3hg.dat_orig'
      ELSE IF ( INSTRUMENT_INDEX.EQ.3 ) THEN
        SOLARFILE = 'lrrs_test/physics_data/SOLAR/solar_omi_uv2.dat'
!        SOLARFILE = 'lrrs_test/physics_data/SOLAR/solar_SOLSPEC.txt'
!        SOLARFILE = 'lrrs_test/physics_data/SOLAR/solar_SOLSPEC_av10nm.txt'
      ELSE
        STOP 'Wrong instrument choice - abort'
      ENDIF

!  Slit function input

!      DO_SLITFUNC_WEIGHTING = .TRUE.
      DO_SLITFUNC_WEIGHTING = .FALSE.

      IF ( INSTRUMENT_INDEX.EQ.1 ) THEN
        SLITFUNC_FWHM         = 0.11d0               ! GOME value = pixe
        SLITFUNC_CUTOFF       = 1.0d-02              ! 1% cutoff
      ELSE IF ( INSTRUMENT_INDEX.EQ.2 ) THEN
        SLITFUNC_FWHM         = 0.11d0               ! SCIA value = pixe
        SLITFUNC_CUTOFF       = 1.0d-02              ! 1% cutoff
      ELSE IF ( INSTRUMENT_INDEX.EQ.3 ) THEN
        SLITFUNC_FWHM         = 0.45d0               ! OMI value = pixel
        SLITFUNC_CUTOFF       = 1.0d-02              ! 1% cutoff
      ENDIF

!  Create PHYSICAL and IOP inputs
!  ==============================

      CALL LRRS_IOPSETUP_O3BIN_M ( &
        SOLARFILE, No_O3,                  & !Inputs
        LAMBDA_START, LAMBDA_FINISH,       & !Inputs
        DO_SLITFUNC_WEIGHTING,             & !Inputs
        SLITFUNC_FWHM, SLITFUNC_CUTOFF,    & !Inputs
        NPOINTS_SURF, ALBEDO, BRDF_Sup_In, & !Inputs
        LRRS_FixIn, LRRS_ModIn, LRRS_Sup,  & !InOut
        GASVOD, FAIL, N_MESSAGES, MESSAGES ) !Outputs

!  Exception handling. Go to the finish and write errors.

      IF ( FAIL ) THEN
        WRITE(*,*)
        WRITE(*,*) 'FAIL after LRRS_IOPSETUP_O3BIN_M'
        GO TO 5678
      ENDIF

!  Define some local variables

      LAMBDAS_RANKED = LRRS_FixIn%Spect%LAMBDAS_RANKED
      FLUXES_RANKED  = LRRS_FixIn%Spect%FLUXES_RANKED
      NPOINTS_INNER  = LRRS_FixIn%Spect%NPOINTS_INNER
      OFFSET_INNER   = LRRS_FixIn%Spect%OFFSET_INNER

!  Call to LRRS master

      LRRS_FixIn%Cont%PROGRESS = 1
      CALL CPU_TIME(E3)
      CALL RAMAN_MASTER (LRRS_FixIn, LRRS_ModIn, LRRS_Sup, LRRS_Out)
      CALL CPU_TIME(E4)
      StdTime = E4 - E3

!  Exception handling. Go to the finish and write errors.

      IF ( LRRS_Out%Status%STATUS_CALCULATION .NE. LRRS_SUCCESS ) THEN
        WRITE(*,*)
        WRITE(*,*) 'FAIL after RAMAN_MASTER'
        N_MESSAGES = LRRS_Out%Status%N_MESSAGES
        MESSAGES   = LRRS_Out%Status%MESSAGES
        GO TO 5678
      END IF

!  Define some local output variables

      ELASTIC_UP = LRRS_Out%Main%ELASTIC_UP
      ELASTIC_DN = LRRS_Out%Main%ELASTIC_DN
      RAMAN_UP   = LRRS_Out%Main%RAMAN_UP
      RAMAN_DN   = LRRS_Out%Main%RAMAN_DN

!  Compute filling (Upwelling only).

      IF ( DO_UPWELLING ) THEN
        DO ATP = 1, NPOINTS_INNER
          DO V = 1, N_GEOMETRIES
            DO OD = 1, N_LOUTPUT
              I_ELA = ELASTIC_UP(OD,V,ATP)
              I_RRS = RAMAN_UP(OD,V,ATP)
              RRSFILLING_UP(OD,V,ATP) = ZERO
              IF ( DABS(I_RRS) .GT. 1.0D-12 ) THEN
                RATIO = I_ELA/I_RRS
                RRSFILLING_UP(OD,V,ATP) = (ONE-RATIO)*100.0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Compute filling (Downwelling only).

      IF ( DO_DNWELLING ) THEN
        DO ATP = 1, NPOINTS_INNER
          DO V = 1, N_GEOMETRIES
            DO OD = 1, N_LOUTPUT
              I_ELA = ELASTIC_DN(OD,V,ATP)
              I_RRS = RAMAN_DN(OD,V,ATP)
              RRSFILLING_DN(OD,V,ATP) = ZERO
              IF ( DABS(I_RRS) .GT. 1.0D-12 ) THEN
                RATIO = I_ELA/I_RRS
                RRSFILLING_DN(OD,V,ATP) = (ONE-RATIO)*100.0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Timing test  (to compare with LRRS2.5a using OpenMP)
!  ===========

      write(*,*)
      write(*,'(11x,a)')                'Timing report'
      write(*,'(4(1x,a))') 'Numwvn','# Threads','Time (sec)'
      write(*,'(4(1x,a))') '------','---------','----------'
      write(*,'(1x,i5,6x,i1,6x,f6.2)') &
        NPOINTS_INNER, 1, StdTime

!  Setup output files

                   !         1         2         3         4         5         6
                   !123456789012345678901234567890123456789012345678901234567890
      IF ( DO_ENERGY_BALANCING ) THEN
        FILNAM36 = 'lrrs_test/results_o3bin_M_tester_EB_SZA**_lamb.resg'
      ELSE
        FILNAM36 = 'lrrs_test/results_o3bin_M_tester_CR_SZA**_lamb.resg'
      ENDIF

      IF (No_O3) FILNAM36(19:20) = 'Ry'

      ISZA = INT(LRRS_FixIn%Cont%SOLAR_ANGLE)
      C2 = '00'
      IF ( ISZA.LT.10 ) WRITE(C2(2:2),'(I1)') ISZA
      IF ( ISZA.GE.10 ) WRITE(C2,'(I2)') ISZA
      FILNAM36(40:41) = C2

      IF (LRRS_FixIn%Bool%DO_BRDF_SURFACE) FILNAM36(43:46) = 'brdf'

      OPEN(36,FILE=TRIM(FILNAM36),STATUS='UNKNOWN')
      IF ( WRITE_TO_PLOTFILE ) THEN
        FILNAM37 = FILNAM36
        FILNAM37(48:51) = 'plot'
        OPEN(37,FILE=TRIM(FILNAM37),STATUS='UNKNOWN')
      ENDIF

!  First line

      WRITE(36,'(/a,a,a/)') &
        ' Wavelength   Solar Flux', &
        '   Geometry   Level/Output', &
        '       Elastic I      Raman  I       Filling %'

!  Write results to file
!  =====================

!  Start wavelength loop

      DO WP = 1, NPOINTS_INNER
        CTP = WP + OFFSET_INNER

        DO V = 1, N_GEOMETRIES
          DO OD = 1, N_LOUTPUT
            WRITE(36,366)LAMBDAS_RANKED(CTP),FLUXES_RANKED(CTP), &
            V,'Upwelling@',DBLE(LBOUNDARIES_OUTPUT(OD)), &
            ELASTIC_UP(OD,V,WP), RAMAN_UP(OD,V,WP), RRSFILLING_UP(OD,V,WP)
          ENDDO

      !     DO OD = 1, N_LOUTPUT
      !       WRITE(36,366)LAMBDAS_RANKED(CTP),FLUXES_RANKED(CTP), &
      !       V,'Dnwelling@',DBLE(LBOUNDARIES_OUTPUT(OD)), &
      !       ELASTIC_DN(OD,V,WP), RAMAN_DN(OD,V,WP), RRSFILLING_DN(OD,V,WP)
      !     ENDDO
      !     WRITE(36,*)' '
        ENDDO
366   FORMAT(1X,F10.5,E14.6,3X,I3,5X,A,F6.2,1P3E15.6)

!  Graphical output

      IF ( WRITE_TO_PLOTFILE ) THEN
        WRITE(37,'(f10.5,1pe14.6,1p500e15.6)') &
          LAMBDAS_RANKED(CTP),FLUXES_RANKED(CTP), &
            (( ELASTIC_UP(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT), &
            ((   RAMAN_UP(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT), &
          ((RRSFILLING_UP(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT), &
            (( ELASTIC_DN(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT), &
            ((   RAMAN_DN(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT), &
          ((RRSFILLING_DN(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT)
      ENDIF

!  End wavelength loop

      ENDDO

!  Display timing (off for LRRS2.5a comparison)

      !CALL CPU_TIME(E2)
      !WRITE(*,*)
      !WRITE(*,'(1X,A,F6.2,A)') 'Test time = ',E2-E1,' sec'

!  Close files and exit

      CLOSE(36)
      IF ( WRITE_TO_PLOTFILE ) CLOSE(37)
      GO TO 901

!  Error section. Write messages to LOGFILE and stop

 5678 CONTINUE
      OPEN( 1, FILE='LRRS_LOGFILE', STATUS = 'UNKNOWN' )
      WRITE(1,'(a,I4)')'Number of messages = ',N_MESSAGES
      DO N = 1, N_MESSAGES
        WRITE(1,'(a)')MESSAGES(N)
      ENDDO
      CLOSE(1)
      WRITE(*,*)
      STOP ' ---> Program failure: Look at LRRS_LOGFILE'

!  Normal stop

901   CONTINUE
      WRITE(*,*)
      STOP ' ---> Program successfully executed'

      END PROGRAM O3BIN_M
