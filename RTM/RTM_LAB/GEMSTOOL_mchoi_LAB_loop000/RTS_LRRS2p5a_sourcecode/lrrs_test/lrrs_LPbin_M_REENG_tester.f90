
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

      PROGRAM LPBIN_M_REENG

!  Used modules

      USE LRRS_PARS_m
      USE LRRS_IO_DEFS_m
      USE LRRS_LIN_IO_DEFS_m

      USE LRRS_IO_CHECK_m
      USE LRRS_MAIN_MASTER_m

      USE LRRS_L_IO_CHECK_m
      USE LPBIN_M_IOPSETUP_m
      USE LRRS_L_MAIN_MASTER_m

      USE BRDF_SUP_MOD_m
      USE LRRS_SUP_ACCESSORIES_1_m

      IMPLICIT NONE

!  LIDORT-RRS input & output type structures

      TYPE(LRRS_Fixed_Inputs)        :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs)     :: LRRS_ModIn
      TYPE(LRRS_Sup_InOut)           :: LRRS_Sup
      TYPE(LRRS_LinInputs)           :: LRRS_LinIn
      TYPE(LRRS_LinSup_InOut)        :: LRRS_LinSup
      TYPE(LRRS_Outputs)             :: LRRS_Out
      TYPE(LRRS_LinOutputs)          :: LRRS_LinOut

!  BRDF supplement / LRRS BRDF-related inputs consistency check status

      TYPE(LRRS_Exception_Handling)  :: LRRS_BRDFCheck_Status

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

!  Number of layers

      INTEGER ::   NLAYERS

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
!    These are defined on the outer grid for the binning realization
!    These are defined on 234 points for the monochromatic
!    Depends on the choice of solar spectrum

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

!  Fourier summed output, Atmospheric Profile Jacobians

      REAL(FPK) :: LP_ELASTIC_UP &
        ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: LP_ELASTIC_DN &
        ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  RAMAN FIELD
!  -----------

!  Radiance including inelastic scattering

      REAL(FPK) :: RAMAN_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: RAMAN_DN &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Fourier summed output, Atmospheric Profile Jacobians

      REAL(FPK) :: LP_RAMAN_UP &
        ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: LP_RAMAN_DN &
        ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  General error handling

      INTEGER ::   STATUS_INPUTREAD
      INTEGER ::   N_MESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( MAX_MESSAGES )

!  Helper variables
!  ================

!  Start and finish wavelengths

      REAL(FPK) :: LAMBDA_START, LAMBDA_FINISH

!  Solar spectrum

      INTEGER, PARAMETER :: MAX_INPUT_POINTS = 1007

      INTEGER   :: N_INPUT_POINTS
      REAL(FPK) :: INPUT_LAMBDAS ( MAX_INPUT_POINTS )
      REAL(FPK) :: INPUT_FLUXES  ( MAX_INPUT_POINTS )

!  Height grid

      REAL(FPK) :: HEIGHT_GRID ( 0:MAX_LAYERS )

!  Temperature-shift [K]

      REAL(FPK) :: TSHIFT

!  Ozone

      REAL(FPK) :: O3COLUMN !total O3 in [DU]

!  Spectrum of gas optical depth

      REAL(FPK) :: GASVOD ( MAX_POINTS )

!  Aerosol information
!   same for all wavelengths.

      LOGICAL ::   DO_AEROSOL
      REAL(FPK) :: AEROSOL_SSALB, AEROSOL_ASYMM, AEROSOL_OPDEP
      INTEGER ::   AEROSOL_TOPLEVEL

!  Surface

      INTEGER   :: NPOINTS_SURF
      REAL(FPK) :: ALBEDO

!  Filling

      REAL(FPK) :: RRSFILLING_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: RRSFILLING_DN &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Saved results

      REAL(FPK) :: BE_UP ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: BE_DN ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: BR_UP ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: BR_DN ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!mick mod 7/20/2016 - changed 1st dim from 4 to MAX_ATMOSWFS
      REAL(FPK) :: SE_UP &
        ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: SE_DN &
        ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: SR_UP &
        ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: SR_DN &
        ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Instrument control (1 = GOME, 2 = SCIA, 3 = OMI)

      INTEGER ::   INSTRUMENT_INDEX

!  Jacobian control

      LOGICAL ::   DO_JAC ( MAX_ATMOSWFS )
      INTEGER ::   JAC, NJACS

!  FD input

      LOGICAL ::   DO_FD
      INTEGER ::   KD, K, ND
      REAL(FPK) :: EPSFAC, EPS

!  Error handling for Physics and IOP setup

      LOGICAL ::   FAIL
      !CHARACTER (LEN=70) :: MESSAGE_SUB

!  Other local variables

      INTEGER ::   ATP, CTP, I, ISZA, N, NW, N_GEOMETRIES, OD, V, W, WP
      REAL(FPK) :: I_ELA, I_RRS, RATIO

!  File name & header control

      CHARACTER (LEN=256) :: SOLARFILE

      CHARACTER (LEN=1)   :: CNJ
      CHARACTER (LEN=2)   :: C2
      CHARACTER (LEN=4)   :: C4
      CHARACTER (LEN=10)  :: JACFORM
      CHARACTER (LEN=56)  :: JACHEADER ( MAX_ATMOSWFS )
      CHARACTER (LEN=80)  :: FILNAM32, FILNAM37

!  Writing to plot file

      LOGICAL :: WRITE_TO_PLOTFILE

!  Timing tests

      REAL :: E1,E2

!  OMP control

      LOGICAL, parameter :: Monitor_CPU = .true.
!     LOGICAL, parameter :: Monitor_CPU = .false.

!      INTEGER, parameter :: NCORES = 1
!      INTEGER, parameter :: NCORES = 2
      INTEGER, parameter :: NCORES = 4

!  OpenMP tests (timing)

      INTEGER :: TIME_DIVIDER
      REAL    :: OMP_E1, OMP_E2, StdTime, LinTime

!  Start program

!  Start timing (off for OpenMP timing)

      !CALL CPU_TIME(E1)

!  INPUT SECTION
!  =============

!  Start and Finish

!     Will calculate 105 points (for OMP testing)
!      LAMBDA_START  = 324.5
!      LAMBDA_FINISH = 336.5

!     Will calculate 52 points (for OMP testing)
      LAMBDA_START  = 324.5
      LAMBDA_FINISH = 330.5

!     Will calculate 8 points
!      LAMBDA_START  = 324.5
!      LAMBDA_FINISH = 325.5

!     Will calculate 2 points (for debugging)
!      LAMBDA_START  = 324.5
!      LAMBDA_FINISH = 324.8

!     Will calculate 1 point (for debugging)
!      LAMBDA_START  = 324.5
!      LAMBDA_FINISH = 324.65

!  Instrument control (1 = GOME, 2 = SCIA, 3 = OMI)

      INSTRUMENT_INDEX = 1

!  Set tshift, total O3, aerosols, & surface albedo

      TSHIFT           = 1.0D0   !in K
      O3COLUMN         = 350.0D0 !in DU

      !DO_AEROSOL       = .FALSE.
      DO_AEROSOL       = .TRUE.

      AEROSOL_TOPLEVEL = 10          ! For 16-layer* atmosphere
!      AEROSOL_TOPLEVEL = 24          ! For 34-layer atmosphere
      AEROSOL_OPDEP    = 0.5D0
      AEROSOL_SSALB    = 0.85D0
      AEROSOL_ASYMM    = 0.75D0

!  Surface
!  Note: NPOINTS_SURF should be = 1 or NPOINTS_OUTER

      NPOINTS_SURF = 1
      ALBEDO       = 0.1D0

!mick mod 7/20/2016 - defined driver array DO_JAC for flexibility in doing different
!                     Jacobians during a given test
!  Jacobian control

      DO_JAC = .FALSE.

!  Jac 1, Ozone profile
!  Jac 2, Aerosol opdep profile
!  Jac 3, Air profile
!  Jac 4, Temperature profile

      DO_JAC(1) = .TRUE.
      DO_JAC(2) = .TRUE. 
      DO_JAC(3) = .TRUE. 
      DO_JAC(4) = .TRUE. 

!  Extra output

      !WRITE_TO_PLOTFILE = .TRUE.
      WRITE_TO_PLOTFILE = .FALSE.

!mick mod 7/20/2016 - added aerosol flag check
!  Check aerosol flags for consistency
      if ( .not.DO_AEROSOL .and. DO_JAC(2) ) then
        WRITE(*,*) 'Driver aerosol flag error'
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = 'Aerosol flag DO_AEROSOL and aerosol Jacobian flag DO_JAC(2) inconsistent'
        GO TO 5678
      endif

!  LRRS input control read (USER-DEFINED INPUT)

      CALL LRRS_READ_INPUTS ( &
        'lrrs_test/lrrs_LPbin_M_REENG_tester.cfg', &
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

!  Initialize LRRS linearized inputs

      CALL LRRS_L_INIT_INPUTS ( LRRS_LinIn )

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

!  This is binning, no MONO

      IF ( .not. DO_BIN_REALIZATION ) &
        STOP 'This is binning - set flag!!'

!  Initialize LRRS standard & linearized supplement inputs

      CALL LRRS_Sup_Init ( LRRS_Sup )
      CALL LRRS_LinSup_Init ( LRRS_LinSup )

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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@  FIRST CALL, BASELINE CALCULATION
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Jac 1, Ozone profile
!  Jac 2, Aerosol opdep profile
!  Jac 3, Air profile
!  Jac 4, Temperature profile

      WRITE(*,*)'*****************************'
      WRITE(*,*)'Doing baseline calculation---'
      WRITE(*,*)'*****************************'

!  Define LRRS linearization control variables

      LRRS_LinIn%Cont%DO_PROFILE_WFS     = .TRUE.
      LRRS_LinIn%Cont%DO_NORMALIZED_WFS  = .TRUE.
      LRRS_LinIn%Cont%DO_AIRPROFILE_WFS  = DO_JAC(3)
      LRRS_LinIn%Cont%DO_TEMPPROFILE_WFS = DO_JAC(4)

!  No FD

      JAC = 999; DO_FD = .FALSE.; ND = 999; EPS = 0.0D0; EPSFAC = 1.0D0

!  Create PHYSICAL and IOP inputs

      CALL LRRS_IOPSETUP_LPBIN_M ( &
        SOLARFILE, DO_JAC, JAC, DO_FD, ND, EPSFAC,    & !Inputs
        LAMBDA_START, LAMBDA_FINISH,                  & !Inputs
        TSHIFT, O3COLUMN,                             & !Inputs
        DO_AEROSOL, AEROSOL_TOPLEVEL,                 & !Inputs
        AEROSOL_SSALB, AEROSOL_ASYMM, AEROSOL_OPDEP,  & !Inputs
        NPOINTS_SURF, ALBEDO, BRDF_Sup_In,            & !Inputs
        LRRS_FixIn, LRRS_ModIn, LRRS_Sup, LRRS_LinIn, & !InOut
        GASVOD, FAIL, N_MESSAGES, MESSAGES )            !Outputs

!  Exception handling. Go to the finish and write errors.

      IF ( FAIL ) THEN
        WRITE(*,*)
        WRITE(*,*) 'FAIL after LRRS_IOPSETUP_LPBIN_M, Baseline call'
        GO TO 5678
      ENDIF

!  Define some local variables

      NLAYERS        = LRRS_FixIn%Cont%NLAYERS
      LAMBDAS_RANKED = LRRS_FixIn%Spect%LAMBDAS_RANKED
      FLUXES_RANKED  = LRRS_FixIn%Spect%FLUXES_RANKED
      NPOINTS_INNER  = LRRS_FixIn%Spect%NPOINTS_INNER
      OFFSET_INNER   = LRRS_FixIn%Spect%OFFSET_INNER

!mick mod 7/20/2016 - define NJACS
      NJACS          = LRRS_LinIn%Cont%LAYER_VARY_NUMBER(1)

!  Heights

      HEIGHT_GRID = LRRS_FixIn%Atmos%HEIGHT_GRID
      LRRS_ModIn%MAtmos%GEOMETRY_SPECHEIGHT = LRRS_FixIn%Atmos%HEIGHT_GRID(NLAYERS)

!  L_Master call

      LRRS_FixIn%Bool%DO_TIMING = Monitor_CPU
      LRRS_FixIn%Cont%PROGRESS  = 1
      LRRS_FixIn%Cont%NTHREADS  = NCORES

      CALL CPU_TIME(OMP_E1)
      CALL L_RAMAN_MASTER &
        ( LRRS_FixIn, LRRS_ModIn, LRRS_Sup, & ! Inputs
          LRRS_LinIn, LRRS_LinSup,          & ! Inputs
          LRRS_Out, LRRS_LinOut )             ! Output
      CALL CPU_TIME(OMP_E2)
      LinTime = OMP_E2 - OMP_E1

!  Exception handling. Go to the finish and write errors.

      IF ( LRRS_Out%Status%STATUS_CALCULATION .NE. LRRS_SUCCESS ) THEN
        WRITE(*,*)
        WRITE(*,*) 'FAIL after L_RAMAN MASTER call'
        N_MESSAGES = LRRS_Out%Status%N_MESSAGES
        MESSAGES   = LRRS_Out%Status%MESSAGES
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = 'FAIL after L_RAMAN MASTER call'
        GO TO 5678
      ENDIF

!  Define some local variables

      ELASTIC_UP    = LRRS_Out%Main%ELASTIC_UP
      ELASTIC_DN    = LRRS_Out%Main%ELASTIC_DN
      RAMAN_UP      = LRRS_Out%Main%RAMAN_UP
      RAMAN_DN      = LRRS_Out%Main%RAMAN_DN

      LP_ELASTIC_UP = LRRS_LinOut%Atmos%LP_ELASTIC_UP
      LP_ELASTIC_DN = LRRS_LinOut%Atmos%LP_ELASTIC_DN
      LP_RAMAN_UP   = LRRS_LinOut%Atmos%LP_RAMAN_UP
      LP_RAMAN_DN   = LRRS_LinOut%Atmos%LP_RAMAN_DN

!  Save Baseline results and compute fillings

      IF ( DO_UPWELLING ) THEN
        DO V = 1, N_GEOMETRIES
          DO OD = 1, N_LOUTPUT
            DO ATP = 1, NPOINTS_INNER
              I_ELA = ELASTIC_UP(OD,V,ATP)
              I_RRS = RAMAN_UP(OD,V,ATP)
              BE_UP(OD,V,ATP) = I_ELA
              BR_UP(OD,V,ATP) = I_RRS
              RRSFILLING_UP(OD,V,ATP) = ZERO
              IF ( DABS(I_RRS) .GT. 1.0D-12 ) THEN
                RATIO = I_ELA/I_RRS
                RRSFILLING_UP(OD,V,ATP) = (ONE-RATIO)*100.0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF ( DO_DNWELLING ) THEN
        DO V = 1, N_GEOMETRIES
          DO OD = 1, N_LOUTPUT
            DO ATP = 1, NPOINTS_INNER
              I_ELA = ELASTIC_DN(OD,V,ATP)
              I_RRS = RAMAN_DN(OD,V,ATP)
              BE_DN(OD,V,ATP) = I_ELA
              BR_DN(OD,V,ATP) = I_RRS
              RRSFILLING_DN(OD,V,ATP) = ZERO
              IF ( DABS(I_RRS) .GT. 1.0D-12 ) THEN
                RATIO = I_ELA/I_RRS
                RRSFILLING_DN(OD,V,ATP) = (ONE-RATIO)*100.0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@  FOUR CALLS, FD PERTURBATIONS 1-4, ALL LAYERS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  FD 1, Ozone    profile perturbation, Turn off linearization flags
!  FD 2, Aerosol  profile perturbation, Turn off linearization flags
!  FD 3, Air Col. profile perturbation, Turn off linearization flags
!  FD 4, Temp.    profile perturbation, Turn off linearization flags

      DO_FD = .TRUE.

      eps = 0.001d0; epsfac = 1.0d0 + eps

!  Define LRRS linearization control variables

      LRRS_LinIn%Cont%DO_PROFILE_WFS     = .FALSE.
      LRRS_LinIn%Cont%DO_NORMALIZED_WFS  = .FALSE.
      LRRS_LinIn%Cont%DO_AIRPROFILE_WFS  = .FALSE.
      LRRS_LinIn%Cont%DO_TEMPPROFILE_WFS = .FALSE.

!  FD loops

      StdTime = 0.0

!mick mod 7/20/2016 - changed loop control for flexibility in doing different Jacobians
!                     during a given test: JAC replaces KD's role as absolute Jacobian 
!                     index; KD now assumes role as relative Jacobian index dependent on
!                     values of array DO_JAC
!                   - MAX_ATMOSWFS introduced for loop limit flexibility

!     DO KD = 1, 4
      KD = 0
      DO JAC = 1, MAX_ATMOSWFS
        IF ( .not.DO_JAC(JAC) ) CYCLE
        KD = KD + 1

        WRITE(*,*)'*****************************'
        WRITE(*,*)'Doing FD calculation # ',KD
        WRITE(*,*)'*****************************'

        DO ND = 1, NLAYERS
          WRITE(*,*) '  Doing FD layer # ',ND

!  Create PHYSICAL and IOP inputs

          CALL LRRS_IOPSETUP_LPBIN_M ( &
            SOLARFILE, DO_JAC, JAC, DO_FD, ND, EPSFAC,    & !Inputs
            LAMBDA_START, LAMBDA_FINISH,                  & !Inputs
            TSHIFT, O3COLUMN,                             & !Inputs
            DO_AEROSOL, AEROSOL_TOPLEVEL,                 & !Inputs
            AEROSOL_SSALB, AEROSOL_ASYMM, AEROSOL_OPDEP,  & !Inputs
            NPOINTS_SURF, ALBEDO, BRDF_Sup_In,            & !Inputs
            LRRS_FixIn, LRRS_ModIn, LRRS_Sup, LRRS_LinIn, & !InOut
            GASVOD, FAIL, N_MESSAGES, MESSAGES )            !Outputs

!  Exception handling. Go to the finish and write errors.

          IF ( FAIL ) THEN
            WRITE(*,*)
            WRITE(*,*) 'FAIL,LRRS_IOPSETUP_LPBIN_M, FD/ND call # ',KD,ND
            GO TO 5678
          ENDIF

!  Define some local variables

          NLAYERS        = LRRS_FixIn%Cont%NLAYERS
          LAMBDAS_RANKED = LRRS_FixIn%Spect%LAMBDAS_RANKED
          FLUXES_RANKED  = LRRS_FixIn%Spect%FLUXES_RANKED
          NPOINTS_INNER  = LRRS_FixIn%Spect%NPOINTS_INNER
          OFFSET_INNER   = LRRS_FixIn%Spect%OFFSET_INNER

!  Geometry Specheight

          LRRS_ModIn%MAtmos%GEOMETRY_SPECHEIGHT = &
            LRRS_FixIn%Atmos%HEIGHT_GRID(NLAYERS)

!  Master call

          LRRS_FixIn%Bool%DO_TIMING = Monitor_CPU
          LRRS_FixIn%Cont%PROGRESS  = 1
          LRRS_FixIn%Cont%NTHREADS  = NCORES

          CALL CPU_TIME(OMP_E1)
          CALL RAMAN_MASTER (LRRS_FixIn, LRRS_ModIn, LRRS_Sup, LRRS_Out)
          CALL CPU_TIME(OMP_E2)
          StdTime = StdTime + (OMP_E2 - OMP_E1)

!  Exception handling. Go to the finish and write errors.

          IF ( LRRS_Out%Status%STATUS_CALCULATION .NE. LRRS_SUCCESS ) THEN
            WRITE(*,*)
            WRITE(*,*) 'FAIL after RAMAN MASTER call, FD/ND ', KD, ND
            N_MESSAGES = LRRS_Out%Status%N_MESSAGES
            MESSAGES   = LRRS_Out%Status%MESSAGES
            GO TO 5678
          ENDIF

!  Define some local variables

          ELASTIC_UP = LRRS_Out%Main%ELASTIC_UP
          ELASTIC_DN = LRRS_Out%Main%ELASTIC_DN
          RAMAN_UP   = LRRS_Out%Main%RAMAN_UP
          RAMAN_DN   = LRRS_Out%Main%RAMAN_DN

!  Save Baseline results

          DO V = 1, N_GEOMETRIES
            DO OD = 1, N_LOUTPUT
              DO ATP = 1, NPOINTS_INNER
                SE_UP(KD,ND,OD,V,ATP) = ELASTIC_UP(OD,V,ATP)
                SR_UP(KD,ND,OD,V,ATP) = RAMAN_UP(OD,V,ATP)
                SE_DN(KD,ND,OD,V,ATP) = ELASTIC_DN(OD,V,ATP)
                SR_DN(KD,ND,OD,V,ATP) = RAMAN_DN(OD,V,ATP)
              ENDDO
            ENDDO
          ENDDO

!  End layer loop

        ENDDO

!  End FD Jacobian loop

      ENDDO

!  Timing test (OpenMP)
!  ====================

!Note: This timing setup assumes the # of OpenMP threads used in the
!      test is <= the # of computational cores available on the
!      machine on which the test is being run and focuses only on the
!      the time to do the RT computations.  It does not include the
!      time to prepare the optical property inputs.

      !IF (OMP_NTHREADS <= N_CORE) THEN
      !  TIME_DIVIDER = OMP_NTHREADS
      !ELSE
        TIME_DIVIDER = NCORES
      !ENDIF

      write(*,*)
      write(*,'(15x,a)')                    'Timing report'
      write(*,'(4(1x,a))') 'Numwvn','# OMP Threads','Time (sec)'
      write(*,'(4(1x,a))') '------','-------------','----------'
      write(*,'(1x,i5,8x,i1,7x,f6.2)') &
        NPOINTS_INNER, NCORES, (LinTime + StdTime)/REAL(TIME_DIVIDER)

!  Setup output files
!  ==================

      ISZA = INT(LRRS_FixIn%Cont%SOLAR_ANGLE)
      C2 = '00'
      IF ( ISZA.LT.10)  WRITE(C2(2:2),'(I1)')ISZA
      IF ( ISZA.GE.10 ) WRITE(C2,'(I2)')ISZA

      IF ( .NOT.LRRS_FixIn%Bool%DO_BRDF_SURFACE ) THEN
        C4 = 'lamb'
      ELSE
        C4 = 'brdf'
      ENDIF

      IF ( DO_ENERGY_BALANCING ) THEN
        FILNAM32 = 'lrrs_test/results_LPAbin_M_tester_EB_SZA' //C2// '_' //C4// '.resg_REENG'
      ELSE
        FILNAM32 = 'lrrs_test/results_LPAbin_M_tester_CR_SZA' //C2// '_' //C4// '.resg_REENG'
      ENDIF

      OPEN(32,FILE=TRIM(FILNAM32),STATUS='REPLACE')
      IF ( WRITE_TO_PLOTFILE ) THEN
        FILNAM37 = FILNAM32
        FILNAM37(49:52) = 'plot'
        OPEN(37,FILE=TRIM(FILNAM37),STATUS='REPLACE')
      ENDIF

!  First lines
!mick mod 7/20/2016 - make headers more flexible

      WRITE(32,'(/3a/)')'Wavelength     Flux  ',&
                        '   Geom   Output-Level', &
                        '      Elastic I      Raman I      Filling %'

      WRITE(cnj,'(i1)') njacs
      jacform='(A,' // cnj // 'A)'
      jacheader(1:njacs) = '     Elastic PJ   Elastic PFD     Raman PJ     Raman PFD'

      WRITE(32,jacform)' Layer #    Midheight ',jacheader(1:njacs)

!  Write results to file
!  =====================

!mick mod 7/20/2016 - use NJACS & implied do loops
      DO W = 1, NPOINTS_INNER
        CTP = W + OFFSET_INNER

!  Start geometry loop

        DO V = 1, N_GEOMETRIES
          DO OD = 1, N_LOUTPUT
            WRITE(32,326)LAMBDAS_RANKED(CTP),FLUXES_RANKED(CTP), &
              V,'Upwelling@',DBLE(LBOUNDARIES_OUTPUT(OD)), &
              BE_UP(OD,V,W),BR_UP(OD,V,W),RRSFILLING_UP(OD,V,W)
            DO K = 1, NLAYERS
              !WRITE(32,327)K,(HEIGHT_GRID(K-1)+HEIGHT_GRID(K))*0.5, &
              !  LP_ELASTIC_UP(1,K,OD,V,W), (SE_UP(1,K,OD,V,W)-BE_UP(OD,V,W))/EPS, &
              !  LP_RAMAN_UP  (1,K,OD,V,W), (SR_UP(1,K,OD,V,W)-BR_UP(OD,V,W))/EPS, &
              !  LP_ELASTIC_UP(2,K,OD,V,W), (SE_UP(2,K,OD,V,W)-BE_UP(OD,V,W))/EPS, &
              !  LP_RAMAN_UP  (2,K,OD,V,W), (SR_UP(2,K,OD,V,W)-BR_UP(OD,V,W))/EPS, &
              !  LP_ELASTIC_UP(3,K,OD,V,W), (SE_UP(3,K,OD,V,W)-BE_UP(OD,V,W))/EPS, &
              !  LP_RAMAN_UP  (3,K,OD,V,W), (SR_UP(3,K,OD,V,W)-BR_UP(OD,V,W))/EPS, &
              !  LP_ELASTIC_UP(4,K,OD,V,W), (SE_UP(4,K,OD,V,W)-BE_UP(OD,V,W))/EPS, &
              !  LP_RAMAN_UP  (4,K,OD,V,W), (SR_UP(4,K,OD,V,W)-BR_UP(OD,V,W))/EPS

              WRITE(32,327)K,(HEIGHT_GRID(K-1)+HEIGHT_GRID(K))*0.5, &
                ( LP_ELASTIC_UP(KD,K,OD,V,W), (SE_UP(KD,K,OD,V,W)-BE_UP(OD,V,W))/EPS, &
                  LP_RAMAN_UP  (KD,K,OD,V,W), (SR_UP(KD,K,OD,V,W)-BR_UP(OD,V,W))/EPS, &
                  KD = 1, NJACS )
            ENDDO
          ENDDO

          DO OD = 1, N_LOUTPUT
            WRITE(32,326)LAMBDAS_RANKED(CTP),FLUXES_RANKED(CTP), &
              V,'Dnwelling@',DBLE(LBOUNDARIES_OUTPUT(OD)), &
              BE_DN(OD,V,W),BR_DN(OD,V,W),RRSFILLING_DN(OD,V,W)
            DO K = 1, NLAYERS
              !WRITE(32,327)K,(HEIGHT_GRID(K-1)+HEIGHT_GRID(K))*0.5, &
              !  LP_ELASTIC_DN(1,K,OD,V,W), (SE_DN(1,K,OD,V,W)-BE_DN(OD,V,W))/EPS, &
              !  LP_RAMAN_DN  (1,K,OD,V,W), (SR_DN(1,K,OD,V,W)-BR_DN(OD,V,W))/EPS, &
              !  LP_ELASTIC_DN(2,K,OD,V,W), (SE_DN(2,K,OD,V,W)-BE_DN(OD,V,W))/EPS, &
              !  LP_RAMAN_DN  (2,K,OD,V,W), (SR_DN(2,K,OD,V,W)-BR_DN(OD,V,W))/EPS, &
              !  LP_ELASTIC_DN(3,K,OD,V,W), (SE_DN(3,K,OD,V,W)-BE_DN(OD,V,W))/EPS, &
              !  LP_RAMAN_DN  (3,K,OD,V,W), (SR_DN(3,K,OD,V,W)-BR_DN(OD,V,W))/EPS, &
              !  LP_ELASTIC_DN(4,K,OD,V,W), (SE_DN(4,K,OD,V,W)-BE_DN(OD,V,W))/EPS, &
              !  LP_RAMAN_DN  (4,K,OD,V,W), (SR_DN(4,K,OD,V,W)-BR_DN(OD,V,W))/EPS

              WRITE(32,327)K,(HEIGHT_GRID(K-1)+HEIGHT_GRID(K))*0.5, &
                ( LP_ELASTIC_DN(KD,K,OD,V,W), (SE_DN(KD,K,OD,V,W)-BE_DN(OD,V,W))/EPS, &
                  LP_RAMAN_DN  (KD,K,OD,V,W), (SR_DN(KD,K,OD,V,W)-BR_DN(OD,V,W))/EPS, &
                  KD = 1, NJACS )
            ENDDO
          ENDDO
          WRITE(32,*)' '

!  End geometry loop

        ENDDO

!  Format

326     format(/1x,F8.4,e14.6,2x,i2,2x,a,f6.2,1p3e14.6/)
327     format(3x,i4,4x,F9.3,4x,1p16e14.6)
!326     format(/1x,F8.4,e11.3,2x,i2,2x,a,f6.2,1p3e13.5/)
!327     format(3x,i4,4x,F9.3,4x,1p16e13.5)
!366     format(1x,F8.4,e11.3,2x,i2,2x,a,f6.2,1p7e13.5)

!  Graphical output

        IF (WRITE_TO_PLOTFILE) THEN
          WRITE(37,'(F10.5,1PE14.6,1P500E15.6)') &
            LAMBDAS_RANKED(CTP),FLUXES_RANKED(CTP), &
            ((        BE_UP(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT), &
            ((        BR_UP(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT), &
            ((RRSFILLING_UP(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT), &
            ((        BE_DN(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT), &
            ((        BR_DN(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT), &
            ((RRSFILLING_DN(OD,V,WP), V=1,N_GEOMETRIES), OD=1,N_LOUTPUT)
        END IF

!  End wavelength loop

      ENDDO

!  Display timing (off for OpenMP timing)

      !CALL CPU_TIME(E2)
      !WRITE(*,*)
      !WRITE(*,'(1X,A,F6.2,A)') 'Test time = ',E2-E1,' sec'

!  Close files and exit

      CLOSE(32)
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

      END PROGRAM LPBIN_M_REENG
