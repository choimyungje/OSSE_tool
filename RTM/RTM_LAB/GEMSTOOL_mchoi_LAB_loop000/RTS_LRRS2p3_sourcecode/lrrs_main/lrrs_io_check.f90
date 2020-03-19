! ###############################################################
! #                                                             #
! #                    THE LIDORT_RRS MODEL                     #
! #                                                             #
! #      (LInearized Discrete Ordinate Radiative Transfer)      #
! #       --         -        -        -         -              #
! #                 (Rotational Raman Scatter)                  #
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
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Version      :  2.3                                        #
! #  Release Date :  March 2011                                 #
! #                                                             #
! ###############################################################

!    #########################################################
!    #                                                       #
!    #   This Version of LIDORT_RRS comes with a GNU-style   #
!    #   license. Please read the license carefully.         #
!    #                                                       #
!    #########################################################

! ###############################################################
! #                                                             #
! #   INPUT SUBROUTINES :                                       #
! #                                                             #
! #     Reading control input from file:                        #
! #                                                             #
! #          1. LRRS_READINPUT                                  #
! #                                                             #
! #     These routines are called at start of lrrs_master:      #
! #                                                             #
! #          2  LRRS_CHECK_INPUTS                               #
! #          3. LRRS_DERIVE_INPUTS                              #
! #                                                             #
! ###############################################################


      MODULE lrrs_io_check

      PRIVATE
      PUBLIC :: LRRS_READINPUT,&
                LRRS_CHECK_INPUTS,&
                LRRS_DERIVE_INPUTS

      CONTAINS

      SUBROUTINE LRRS_READINPUT &
        ( FILNAM, EFILNAM, &
          LRRS_FixIn, LRRS_ModIn, &
          STATUS, NMESS, MESS )

!  Read all control inputs for LRRS
!  This does not read the optical properties!

!  Used modules

      USE LRRS_PARS
      USE LRRS_AUX2
      USE LRRS_INPUTS_DEF

      IMPLICIT NONE

!  Module input arguments (input filename, error filename)
!  -------------------------------------------------------

      CHARACTER (LEN=*), INTENT(IN) ::           FILNAM
      CHARACTER (LEN=*), INTENT(IN) ::           EFILNAM

!  Module output arguments (read-in from file)
!  -------------------------------------------

!  LIDORT-RRS input structures

      TYPE(LRRS_Fixed_Inputs), INTENT(OUT)    :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT(OUT) :: LRRS_ModIn

!  Overall status

      INTEGER, INTENT(OUT) ::                    STATUS
      INTEGER, INTENT(OUT) ::                    NMESS
      CHARACTER (LEN=70), INTENT(OUT) ::         MESS ( MAX_MESSAGES )

!  Local variables
!  ---------------

!  Number of discrete ordinate streams

      INTEGER ::   NSTREAMS

!  Number of fine layers for sscorr outgoing

      INTEGER ::   NFINELAYERS

!  Solar beam input geometry

      REAL(FPK) :: SOLAR_ANGLE

!  User stream variables

      INTEGER ::   N_USER_STREAMS
      REAL(FPK) :: USER_ANGLES  ( MAX_USER_STREAMS )

!  User azimuth variables

      INTEGER ::   N_USER_RELAZMS
      REAL(FPK) :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)
!    LBOUNDARIES_OUTPUT = layer boundary choices
!    LPARTIALS_OUTPUT = off-boundary level choices

      INTEGER ::   N_LOUTPUT
      INTEGER ::   LBOUNDARIES_OUTPUT ( MAX_LOUTPUT )
      REAL(FPK) :: LPARTIALS_OUTPUT   ( MAX_LOUTPUT )

!  Overall accuracy for convergence criterion

      REAL(FPK) :: LRRS_ACCURACY

!  Small number control

      REAL(FPK) :: LRRS_FGSMALL

!  Flux factor

      REAL(FPK) :: FLUX_FACTOR

!  Earth Radius

      REAL(FPK) :: EARTH_RADIUS

!  Lambertian fraction

      REAL(FPK) :: LAMBERTIAN_FRACTION

!  Number of BRDF azimuth streams

      INTEGER ::   NSTREAMS_BRDF

!  BRDF parameters.
!     Index choice of BRDF function
!     For Cox-Munk, first parameter  = Wind Speed
!                 second parameter = refractive index
!                 third parameter  = 0.0d

      INTEGER ::            WHICH_BRDF
      INTEGER ::            BRDF_NPARS
      REAL(FPK) ::          BRDF_PARS ( MAX_BRDF_PARAMETERS)
      REAL(FPK) ::          BRDF_FACTOR
      CHARACTER (LEN=10) :: BRDF_NAMES

!  filenames for output

      CHARACTER (LEN=60) :: DEBUG_FILENAMES(4)

!  CONTROL FLAGS EXPANSION
!  -----------------------

      LOGICAL :: &
         DO_RRS_OVERALL,      DO_DIRECTRRS_ONLY,    DO_ELASTIC_ONLY, &
         DO_CABANNES_RAMAN,   DO_ENERGY_BALANCING,  DO_MSMODE_LRRS, &
         DO_BIN_REALIZATION,  DO_MONO_REALIZATION,  DO_DIRECT_BEAM, &
         DO_SSCORR_OUTGOING,  DO_SSFULL,            DO_MOLECSCAT_ONLY, &
         DO_PLANE_PARALLEL,   DO_DELTAM_SCALING,    DO_DOUBLE_CONVTEST, &
         DO_UPWELLING,        DO_DNWELLING,         DO_USER_STREAMS, &
         DO_NO_AZIMUTH,       DO_LBOUNDARIES,       DO_MVOUT_ONLY, &
         DO_ADDITIONAL_MVOUT, DO_LAMBERTIAN_SURFACE,DO_WATERLEAVING, &
         DO_SHADOW_EFFECT,    DO_GLITTER_DBMS,  DO_LAMBERTIAN_FRACTION, &
         DO_LRRS_WRITESCENARIO, DO_LRRS_WRITERESULTS, &
         DO_LRRS_WRITEINPUT,    DO_LRRS_WRITEFOURIER

!  Helper variables

      CHARACTER (LEN=6), PARAMETER :: PREFIX = 'LRRS -'

      LOGICAL ::            ERROR
      CHARACTER (LEN=80) :: PAR_STR
      !LOGICAL ::            GFINDPAR
      INTEGER ::            I, K, FILUNIT
      !INTEGER ::            LEN_STRING
      CHARACTER (LEN=70) :: MAIL, ACTION
      INTEGER ::            N

!  initialize status

      STATUS = LRRS_SUCCESS

!  initialize error message count

      NMESS = 0
      DO N = 1, MAX_MESSAGES
        MESS(N) = ' '
      ENDDO

!  initialize variables
!  ====================

!  number of discrete ordinates

      NSTREAMS = 0

!  Solar angle (degrees)

      SOLAR_ANGLE = ZERO

!  Output options - viewer angles

      N_USER_STREAMS = 0
      N_USER_RELAZMS = 0
      DO I = 1, MAX_USER_STREAMS
        USER_ANGLES(I) = ZERO
      ENDDO
      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO

!  Output options -  layer output choices

      N_LOUTPUT = 0
      DO I = 1, MAX_LOUTPUT
        LPARTIALS_OUTPUT(I)   = ZERO
        LBOUNDARIES_OUTPUT(I) = 0
      ENDDO

!  Fourier cosine series, accuracy

      LRRS_ACCURACY = ZERO

!  Small number control

      LRRS_FGSMALL  = ZERO

!  Flux factor

      FLUX_FACTOR = ONE

! This variable not used

!      SHIFT_ORDER_CUTOFF = ZERO

!  surface variables

      LAMBERTIAN_FRACTION = ZERO
      NSTREAMS_BRDF = 0
      WHICH_BRDF    = 0
      BRDF_NPARS    = 0
      DO I = 1, MAX_BRDF_PARAMETERS
        BRDF_PARS(I) = ZERO
      ENDDO
      BRDF_FACTOR   = ZERO
      BRDF_NAMES    = ' '

!  Expand control flags
!  ====================

!  Top level control

      DO_RRS_OVERALL       = .FALSE.
      DO_DIRECTRRS_ONLY    = .FALSE.
      DO_ELASTIC_ONLY      = .FALSE.

!  Solution method (options are mutually exclusive)

      DO_CABANNES_RAMAN    = .FALSE.
      DO_ENERGY_BALANCING  = .FALSE.

!  MS mode flag

      DO_MSMODE_LRRS       = .FALSE.

!  Binning or monochromatic realization
!    (Options are mutually exclusive)

      DO_BIN_REALIZATION   = .FALSE.
      DO_MONO_REALIZATION  = .FALSE.

!  Direct beam

      DO_DIRECT_BEAM       = .FALSE.

!  Single scattering flags

      DO_SSCORR_OUTGOING   = .FALSE.
      DO_SSFULL            = .FALSE.

!  Molecular scattering only (no aerosols/clouds)

      DO_MOLECSCAT_ONLY    = .FALSE.

!  Plane parallel

      DO_PLANE_PARALLEL    = .FALSE.

!  Delta-M scaling, only with aerosols

      DO_DELTAM_SCALING    = .FALSE.

!  Double pconvergence test (only with aerosols)

      DO_DOUBLE_CONVTEST   = .FALSE.

!  Upwelling/downwelling

      DO_UPWELLING         = .FALSE.
      DO_DNWELLING         = .FALSE.

!  User angles flag

      DO_USER_STREAMS      = .FALSE.

!  No azimuth flag

      DO_NO_AZIMUTH        = .FALSE.

!  level boundaries flag

      DO_LBOUNDARIES       = .FALSE.

!  Mean-value option

      DO_MVOUT_ONLY        = .FALSE.
      DO_ADDITIONAL_MVOUT  = .FALSE.

!  Surface flags

      DO_LAMBERTIAN_SURFACE  = .FALSE.
      DO_WATERLEAVING        = .FALSE.
      DO_SHADOW_EFFECT       = .FALSE.
      DO_GLITTER_DBMS        = .FALSE.
      DO_LAMBERTIAN_FRACTION = .FALSE.

!  Debug write flags

      DO_LRRS_WRITEINPUT     = .FALSE.
      DO_LRRS_WRITERESULTS   = .FALSE.
      DO_LRRS_WRITESCENARIO  = .FALSE.
      DO_LRRS_WRITEFOURIER   = .FALSE.

!  Read section
!  ============

!  Open file

      FILUNIT = LRRS_INUNIT
      OPEN(LRRS_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  ALGORITHM CONTROL VARIABLES
!  ===========================

!  Overall filling flag

      PAR_STR = 'Do Filling calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_RRS_OVERALL
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Filling calculation for single scatter only

      PAR_STR = 'Do RRS only with primary scatter?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_DIRECTRRS_ONLY
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Elastic only flag

      DO_ELASTIC_ONLY = .not. DO_RRS_OVERALL .and. &
                        .not. DO_DIRECTRRS_ONLY

!  Realization control.  Binning or Monochromatic.
!  options are mutually exclusive

      IF ( .NOT. DO_ELASTIC_ONLY ) THEN
        PAR_STR = 'Do Inelastic scattering with binning realization?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_BIN_REALIZATION
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

        PAR_STR = &
            'Do Inelastic scattering with monochromatic realization?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_MONO_REALIZATION
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  Solution method. Use Cabannes-Raman or Energy-balancing
!  options are mutually exclusive (checked here). One must be true!

      PAR_STR = 'Do Energy-balancing treatment?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_ENERGY_BALANCING
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      DO_CABANNES_RAMAN = .not. DO_ENERGY_BALANCING

!  Pseudo-spherical control

      PAR_STR = 'Plane-parallel treatment of direct beam?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_PLANE_PARALLEL
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Scatterers and phase function control

      PAR_STR='Molecular scattering atmosphere only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_MOLECSCAT_ONLY
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Multiple scatter source function output control
!   Not enabled in this version

      DO_MSMODE_LRRS = .FALSE.

!  Solar beam control

      PAR_STR = 'Include direct beam?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_DIRECT_BEAM
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Outgoing Single scatter correction of whole field
!    ENABLED 29 May 2008.

      PAR_STR = 'Single scatter outgoing correction of all fields?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SSCORR_OUTGOING
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Full single scatter, only with Outgoing..

      IF ( DO_SSCORR_OUTGOING ) THEN
        PAR_STR = 'Single scatter outgoing correction -- finish?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SSFULL
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  Double convergence test

      PAR_STR = 'Perform double convergence test?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_DOUBLE_CONVTEST
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  No azimuth dependence (TRUE means Fourier m = 0 only )

      PAR_STR = 'No azimuth dependence in the solution?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_NO_AZIMUTH
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Delta-M scaling

      PAR_STR = 'Include Delta-M scaling?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_DELTAM_SCALING
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Accuracy input for Fourier series convergence

      PAR_STR = 'Fourier series convergence'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) LRRS_ACCURACY
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Small number control for Taylor series expansions
!    Suggested value = 1.0D-03

      PAR_STR = 'Small number control for Taylor series expansions'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) LRRS_FGSMALL
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Earth radius

      PAR_STR = 'Earth radius'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) EARTH_RADIUS
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Flux factor

      PAR_STR = 'Flux factor'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FLUX_FACTOR
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Number of fine layers for single scatter outgoing correction

      IF ( DO_SSCORR_OUTGOING ) THEN
        PAR_STR = &
            'Number of fine layers for single scatter outgoing'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NFINELAYERS
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  2, GEOMETRY AND OUTPUT VARIABLES
!  ================================

!  Directional output control

      PAR_STR = 'Upwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_UPWELLING
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'Downwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_DNWELLING
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  User-defined Stream angle output control

      PAR_STR = 'User-defined stream angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Flag for output for level boundaries
!    - If not set, then output at off-boundary optical depths

      PAR_STR = 'Output at layer boundaries?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) DO_LBOUNDARIES
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      IF ( .NOT. DO_LBOUNDARIES ) THEN
        MAIL   = 'Partial layer output not enabled in LRRS yet'
        NMESS  = NMESS + 1
        MESS(NMESS) = MAIL(1:LEN_STRING(MAIL))
        GO TO 455
      ENDIF

!  Mean-value output control

      PAR_STR = 'Generate only mean value output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_MVOUT_ONLY
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'Generate mean value output additionally?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_ADDITIONAL_MVOUT
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Number of  discrete ordinate streams
!    -- Number of streams is checked now, to see if exceeds dimensioning

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      IF ( NSTREAMS .GT. MAX_STREAMS ) THEN
        MAIL   = 'Number of half-space streams > maximum dimension'
        CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
        GO TO 455
      ENDIF

!  TOA solar zenith angle input

      PAR_STR = 'TOA solar zenith angle (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) SOLAR_ANGLE
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  User defined stream angles (should be positive from 0 to 90)
!    -- Number of angles is checked now, to see if exceeds dimensioning

      IF ( DO_USER_STREAMS ) THEN
        PAR_STR = 'Number of user-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          MAIL   = 'Number of user streams > maximum dimension'
          CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
          GO TO 455
        ENDIF
        PAR_STR = 'User-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
          ENDDO
        ENDIF
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  Relative azimuth input angles if Azimuth flag set
!    -- Number of angles is checked now, to see if exceeds dimensioning
!    -- Set default to 1 value of 0.0 if the flag is not set

      IF ( .NOT. DO_NO_AZIMUTH ) THEN
        PAR_STR = 'Number of user-defined relative azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
          MAIL   = 'Number of relative azimuths > maximum dimension'
          CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
          GO TO 455
        ENDIF
        PAR_STR = 'User-defined relative azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_RELAZMS
            READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
          ENDDO
        ENDIF
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ELSE
        N_USER_RELAZMS = 1
        USER_RELAZMS(1) = ZERO
      ENDIF

!  User defined off-boundary partial layer output choices
!    -- These are specified as follows
!         LPARTIALS = 17.49 means that output will be in Layer 18
!         and .49 of the way down from the top (in height)
!    -- Number is checked now, to see if exceeds dimensioning

      IF ( .NOT.DO_LBOUNDARIES ) THEN
        PAR_STR = 'Number of partial layer output choices'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) N_LOUTPUT
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        IF ( N_LOUTPUT .GT. MAX_LOUTPUT ) THEN
          MAIL   ='Number of output choices > maximum dimension'
          CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
          GO TO 455
        ENDIF
        PAR_STR = 'Partial layer output choices'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_LOUTPUT
            READ (FILUNIT,*,ERR=998) LPARTIALS_OUTPUT(I)
          ENDDO
        ENDIF
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  User defined LAYER boundary optical depths (INTEGER :: input)
!    -- Number is checked now, to see if exceeds dimensioning

      IF ( DO_LBOUNDARIES ) THEN
        PAR_STR = 'Number of layer boundary output choices'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) N_LOUTPUT
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        IF ( N_LOUTPUT .GT. MAX_LOUTPUT ) THEN
          MAIL   ='Number of output choices > maximum dimension'
          CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
          GO TO 455
        ENDIF
        PAR_STR = 'Layer boundaries for output (0 = TOA)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_LOUTPUT
            READ (FILUNIT,*,ERR=998) LBOUNDARIES_OUTPUT(I)
          ENDDO
        ENDIF
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
      ENDIF

!  4. SURFACE VARIABLES
!  ====================

!  Lambertian surface

      PAR_STR = 'Lambertian surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_LAMBERTIAN_SURFACE
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  BRDF input

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN

!  Lambertian fraction to be included ?

        PAR_STR = 'Do Lambertian component as part of BRDF surface?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_LAMBERTIAN_FRACTION
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Set waterleaving flag

        IF ( DO_LAMBERTIAN_FRACTION ) THEN
          PAR_STR = 'Do water-leaving Lambertian contribution?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
               READ (FILUNIT,*,ERR=998) DO_WATERLEAVING
          CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        ENDIF

!  Set lambertian fraction

        IF ( DO_LAMBERTIAN_FRACTION.AND..NOT.DO_WATERLEAVING ) THEN
          PAR_STR = 'Lambertian fraction as part of BRDF surface'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
               READ (FILUNIT,*,ERR=998) LAMBERTIAN_FRACTION
          CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        ENDIF

!  Direct beam correction. Now automatic, subsumed in the SS calculation
!        PAR_STR = 'Do Direct beam correction for BRDF surface?'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
!     &       READ (FILUNIT,*,ERR=998) DO_DB_CORRECTION
!        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Number of BRDF azimuth streams, check this value

        PAR_STR = 'Number of bidirectional reflectance streams'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) NSTREAMS_BRDF
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

        IF ( NSTREAMS_BRDF .GT. MAX_STREAMS_BRDF ) THEN
          MAIL = 'Number of BRDF streams > maximum dimension.'
          ACTION = ' Re-set input value.'
          STATUS = LRRS_SERIOUS
! mick fix 3/22/11 - replaced call statement with following two lines
          !CALL LRRS_ERROR_TRACE ( MAIL, ACTION, STATUS )
          NMESS = NMESS + 1
          MESS(NMESS) = &
            MAIL(1:LEN_STRING(MAIL))//ACTION(1:LEN_STRING(ACTION))
          RETURN
        ENDIF

!  Main kernel input

        PAR_STR = &
       'Kernel name, index, amplitude, # parameters, parameters'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          READ (FILUNIT,56,ERR=998) &
               BRDF_NAMES, WHICH_BRDF, BRDF_FACTOR, &
               BRDF_NPARS,(BRDF_PARS(K),K=1,3)
        ENDIF
        CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
 56     FORMAT( A10, I2, F6.2, I2, 3F12.6 )

!  Shadowing input (for Cox-Munk type)

        IF ( BRDF_NAMES .EQ. 'Cox-Munk  ' ) THEN
          PAR_STR = 'Do shadow effect for glitter kernels?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_SHADOW_EFFECT
          ENDIF
          CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        ENDIF

!  Multiple reflectance DB correction (for Cox-Munk type)

        IF ( BRDF_NAMES .EQ. 'Cox-Munk  ' ) THEN
          PAR_STR = 'Do multiple reflectance for glitter kernels?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_GLITTER_DBMS
          ENDIF
          CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )
        ENDIF

      ENDIF

!  5. DEBUF FLAGS and FILE NAMES
!  =============================

!  Output write flags

      PAR_STR = 'Input control write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_LRRS_WRITEINPUT
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'Input scenario write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_LRRS_WRITESCENARIO
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'Fourier component output write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_LRRS_WRITEFOURIER
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'Results write?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_LRRS_WRITERESULTS
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

!  Output filenames

      PAR_STR = 'filename for input write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,'(a)',ERR=998) DEBUG_FILENAMES(1)
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'filename for scenario write'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,'(a)',ERR=998) DEBUG_FILENAMES(2)
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'filename for Fourier output'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,'(a)',ERR=998) DEBUG_FILENAMES(3)
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      PAR_STR = 'filename for main output'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,'(a)',ERR=998) DEBUG_FILENAMES(4)
      CALL LRRS_ERROR_TRACE ( ERROR, PAR_STR, NMESS, MESS )

      CLOSE(FILUNIT)

!  Check number of messages at this point

      if ( nmess .gt. 0 ) GO TO 455

!  Define LIDORT-RRS inputs
!  ------------------------

!  Fixed Boolean Inputs

      LRRS_FixIn%Bool%DO_RRS_OVERALL         = DO_RRS_OVERALL
      LRRS_FixIn%Bool%DO_DIRECTRRS_ONLY      = DO_DIRECTRRS_ONLY
      LRRS_FixIn%Bool%DO_DIRECT_BEAM         = DO_DIRECT_BEAM
      LRRS_FixIn%Bool%DO_SSCORR_OUTGOING     = DO_SSCORR_OUTGOING
      LRRS_FixIn%Bool%DO_SSFULL              = DO_SSFULL
      LRRS_FixIn%Bool%DO_MOLECSCAT_ONLY      = DO_MOLECSCAT_ONLY
      LRRS_FixIn%Bool%DO_PLANE_PARALLEL      = DO_PLANE_PARALLEL
      LRRS_FixIn%Bool%DO_DELTAM_SCALING      = DO_DELTAM_SCALING
      LRRS_FixIn%Bool%DO_DOUBLE_CONVTEST     = DO_DOUBLE_CONVTEST
      LRRS_FixIn%Bool%DO_UPWELLING           = DO_UPWELLING
      LRRS_FixIn%Bool%DO_DNWELLING           = DO_DNWELLING
      LRRS_FixIn%Bool%DO_LBOUNDARIES         = DO_LBOUNDARIES
      LRRS_FixIn%Bool%DO_MVOUT_ONLY          = DO_MVOUT_ONLY
      LRRS_FixIn%Bool%DO_ADDITIONAL_MVOUT    = DO_ADDITIONAL_MVOUT
      LRRS_FixIn%Bool%DO_LAMBERTIAN_SURFACE  = DO_LAMBERTIAN_SURFACE
      LRRS_FixIn%Bool%DO_WATERLEAVING        = DO_WATERLEAVING
      LRRS_FixIn%Bool%DO_SHADOW_EFFECT       = DO_SHADOW_EFFECT
      LRRS_FixIn%Bool%DO_GLITTER_DBMS        = DO_GLITTER_DBMS
      LRRS_FixIn%Bool%DO_LAMBERTIAN_FRACTION = DO_LAMBERTIAN_FRACTION
      LRRS_FixIn%Bool%DO_LRRS_WRITEINPUT     = DO_LRRS_WRITEINPUT
      LRRS_FixIn%Bool%DO_LRRS_WRITESCENARIO  = DO_LRRS_WRITESCENARIO
      LRRS_FixIn%Bool%DO_LRRS_WRITEFOURIER   = DO_LRRS_WRITEFOURIER
      LRRS_FixIn%Bool%DO_LRRS_WRITERESULTS   = DO_LRRS_WRITERESULTS

!  Modified Boolean Inputs

      LRRS_ModIn%MBool%DO_ELASTIC_ONLY     = DO_ELASTIC_ONLY
      LRRS_ModIn%MBool%DO_CABANNES_RAMAN   = DO_CABANNES_RAMAN
      LRRS_ModIn%MBool%DO_ENERGY_BALANCING = DO_ENERGY_BALANCING
      LRRS_ModIn%MBool%DO_MSMODE_LRRS      = DO_MSMODE_LRRS
      LRRS_ModIn%MBool%DO_BIN_REALIZATION  = DO_BIN_REALIZATION
      LRRS_ModIn%MBool%DO_MONO_REALIZATION = DO_MONO_REALIZATION
      LRRS_ModIn%MBool%DO_USER_STREAMS     = DO_USER_STREAMS
      LRRS_ModIn%MBool%DO_NO_AZIMUTH       = DO_NO_AZIMUTH

!  Fixed Control Inputs

      !LRRS_FixIn%Cont%PROGRESS        = PROGRESS
      !LRRS_FixIn%Cont%NLAYERS         = NLAYERS
      LRRS_FixIn%Cont%NSTREAMS        = NSTREAMS
      !LRRS_FixIn%Cont%NMOMENTS_INPUT  = NMOMENTS_INPUT
      LRRS_FixIn%Cont%NFINELAYERS     = NFINELAYERS
      LRRS_FixIn%Cont%FLUX_FACTOR     = FLUX_FACTOR
      LRRS_FixIn%Cont%SOLAR_ANGLE     = SOLAR_ANGLE
      LRRS_FixIn%Cont%EARTH_RADIUS    = EARTH_RADIUS
      LRRS_FixIn%Cont%LRRS_ACCURACY   = LRRS_ACCURACY
      LRRS_FixIn%Cont%LRRS_FGSMALL    = LRRS_FGSMALL
      LRRS_FixIn%Cont%DEBUG_FILENAMES = DEBUG_FILENAMES

!  Fixed User-Value Inputs

      LRRS_FixIn%UserVal%N_USER_STREAMS     = N_USER_STREAMS
      LRRS_FixIn%UserVal%USER_ANGLES        = USER_ANGLES
      LRRS_FixIn%UserVal%N_USER_RELAZMS     = N_USER_RELAZMS
      LRRS_FixIn%UserVal%USER_RELAZMS       = USER_RELAZMS
      LRRS_FixIn%UserVal%N_LOUTPUT          = N_LOUTPUT
      LRRS_FixIn%UserVal%LBOUNDARIES_OUTPUT = LBOUNDARIES_OUTPUT
      LRRS_FixIn%UserVal%LPARTIALS_OUTPUT   = LPARTIALS_OUTPUT

!  Fixed Spectral Inputs

      !LRRS_FixIn%Spect%LAMBDAS_RANKED = LAMBDAS_RANKED
      !LRRS_FixIn%Spect%FLUXES_RANKED  = FLUXES_RANKED

      !For monochromatic calculations:
      !LRRS_FixIn%Spect%NPOINTS_MONO   = NPOINTS_MONO
      !LRRS_FixIn%Spect%LAMBDA_EXCIT   = LAMBDA_EXCIT
      !LRRS_FixIn%Spect%W_EXCIT        = W_EXCIT

      !For binning calculations:
      !LRRS_FixIn%Spect%NPOINTS_INNER  = NPOINTS_INNER
      !LRRS_FixIn%Spect%OFFSET_INNER   = OFFSET_INNER
      !LRRS_FixIn%Spect%NPOINTS_OUTER  = NPOINTS_OUTER
      !LRRS_FixIn%Spect%BINLOWER       = BINLOWER
      !LRRS_FixIn%Spect%BINUPPER       = BINUPPER

!  Fixed Atmosphere Inputs

      !LRRS_FixIn%Atmos%HEIGHT_GRID        = HEIGHT_GRID
      !LRRS_FixIn%Atmos%LAYER_TEMPERATURES = LAYER_TEMPERATURES
      !LRRS_FixIn%Atmos%LAYER_AIRCOLUMNS   = LAYER_AIRCOLUMNS
      !LRRS_FixIn%Atmos%RAYLEIGH_XSEC      = RAYLEIGH_XSEC
      !LRRS_FixIn%Atmos%RAYLEIGH_DEPOL     = RAYLEIGH_DEPOL

      !LRRS_FixIn%Atmos%DELTAU_INPUT_UNSCALED      = &
      !                 DELTAU_INPUT_UNSCALED
      !LRRS_FixIn%Atmos%OMEGAMOMS_ELASTIC_UNSCALED = &
      !                 OMEGAMOMS_ELASTIC_UNSCALED

!  Fixed Surface Inputs

      LRRS_FixIn%Surf%LAMBERTIAN_FRACTION = LAMBERTIAN_FRACTION
      LRRS_FixIn%Surf%NSTREAMS_BRDF       = NSTREAMS_BRDF
      LRRS_FixIn%Surf%BRDF_NAMES          = BRDF_NAMES
      LRRS_FixIn%Surf%WHICH_BRDF          = WHICH_BRDF
      LRRS_FixIn%Surf%BRDF_FACTOR         = BRDF_FACTOR
      LRRS_FixIn%Surf%BRDF_NPARS          = BRDF_NPARS
      LRRS_FixIn%Surf%BRDF_PARS           = BRDF_PARS
      !LRRS_FixIn%Surf%ALBEDOS_RANKED      = ALBEDOS_RANKED

!  Normal return

      RETURN

!  Open file error - abort

 300  CONTINUE
      MAIL = 'openfile failure for '//FILNAM(1:LEN_STRING(FILNAM))
      CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
      CLOSE(FILUNIT)
      GO TO 455

!  Line read error - abort immediately

998   CONTINUE
      MAIL = 'read failure for '//PAR_STR(1:LEN_STRING(PAR_STR))
      CALL LRRS_ERROR_TRACE ( .TRUE., MAIL, NMESS, MESS )
      CLOSE(FILUNIT)
      GO TO 455

!  Write collected error messages

455   CONTINUE
      STATUS = LRRS_SERIOUS
      OPEN ( 44, FILE = EFILNAM, STATUS = 'UNKNOWN' )
      write(44,'(a,t40,i5)')'Total number of error messages = ',NMESS
      DO n = 1, NMESS
        write(44,'(a4,i5,a)')'  # ',N,MESS(N)
      ENDDO
      CLOSE(44)

!  Finish

      RETURN
      END SUBROUTINE LRRS_READINPUT

!

      SUBROUTINE LRRS_CHECK_INPUTS &
       ( LRRS_FixIn, LRRS_ModIn, N_USER_STREAMS, N_USER_RELAZMS, N_LOUTPUT, &
         SOLAR_ANGLE, USER_ANGLES, USER_RELAZMS, NLAYERS, HEIGHT_GRID, &
         LBOUNDARIES_OUTPUT, LPARTIALS_OUTPUT, EARTH_RADIUS, WHICH_BRDF, &
         LAMBERTIAN_FRACTION, BRDF_NPARS, BRDF_PARS, BRDF_NAMES, &
         STATUS, N_MESSAGES, MESSAGES )

!  Checking of control inputs

!  Used Modules

      USE LRRS_PARS
      USE LRRS_INPUTS_DEF

      IMPLICIT NONE

!  Input arguments
!  ---------------

      TYPE(LRRS_Fixed_Inputs), INTENT(IN)    :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT(IN) :: LRRS_ModIn

!  Number of layers

      INTEGER, INTENT(IN) ::   NLAYERS

!  Solar beam input geometry

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLE

!  User stream variables

      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_ANGLES  ( MAX_USER_STREAMS )

!  User azimuth variables

      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      REAL(FPK), INTENT(IN) :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  Height grid

      REAL(FPK), INTENT(IN) :: HEIGHT_GRID ( 0: MAX_LAYERS )

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)
!    LBOUNDARIES_OUTPUT = layer boundary choices
!    LPARTIALS_OUTPUT = off-boundary level choices

      INTEGER, INTENT(IN) ::   N_LOUTPUT
      INTEGER, INTENT(IN) ::   LBOUNDARIES_OUTPUT ( MAX_LOUTPUT )
      REAL(FPK), INTENT(IN) :: LPARTIALS_OUTPUT   ( MAX_LOUTPUT )

!  Earth Radius

      REAL(FPK), INTENT(IN) :: EARTH_RADIUS

!  Lambertian fraction

      REAL(FPK), INTENT(IN) :: LAMBERTIAN_FRACTION

!  BRDF parameters.
!     Index choice of BRDF function
!     For Cox-Munk, first parameter  = Wind Speed
!                 second parameter = refractive index
!                 third parameter  = 0.0d

      INTEGER, INTENT(IN) ::      WHICH_BRDF
      INTEGER, INTENT(INOUT) ::   BRDF_NPARS
      REAL(FPK), INTENT(INOUT) :: BRDF_PARS ( MAX_BRDF_PARAMETERS)
      CHARACTER (LEN=10), INTENT(IN) ::  BRDF_NAMES

!  Output arguments
!  ----------------

!  Error handling

      INTEGER, INTENT(OUT) ::  STATUS
      INTEGER, INTENT(OUT) ::  N_MESSAGES
      CHARACTER (LEN=70), INTENT(OUT) :: MESSAGES ( MAX_MESSAGES )

!  Local variables
!  ---------------

      LOGICAL :: &
         DO_RRS_OVERALL,      DO_DIRECTRRS_ONLY,    DO_ELASTIC_ONLY, &
         DO_CABANNES_RAMAN,   DO_ENERGY_BALANCING,  DO_MSMODE_LRRS, &
         DO_BIN_REALIZATION,  DO_MONO_REALIZATION,  DO_DIRECT_BEAM, &
         DO_SSCORR_OUTGOING,  DO_SSFULL,            DO_MOLECSCAT_ONLY, &
         DO_PLANE_PARALLEL,   DO_DELTAM_SCALING,    DO_DOUBLE_CONVTEST, &
         DO_UPWELLING,        DO_DNWELLING,         DO_USER_STREAMS, &
         DO_NO_AZIMUTH,       DO_LBOUNDARIES,       DO_MVOUT_ONLY, &
         DO_ADDITIONAL_MVOUT, DO_LAMBERTIAN_SURFACE,DO_WATERLEAVING, &
         DO_SHADOW_EFFECT,    DO_GLITTER_DBMS,  DO_LAMBERTIAN_FRACTION, &
         DO_LRRS_WRITESCENARIO, DO_LRRS_WRITERESULTS, &
         DO_LRRS_WRITEINPUT,    DO_LRRS_WRITEFOURIER

      INTEGER ::   UTA, N, N1, I, NSTART
      REAL(FPK) :: D1
      CHARACTER (LEN=2) :: C2
      LOGICAL ::   LOOP

!  BRDF Kernel names (check)

      CHARACTER (LEN=10) :: BRDF_CHECK_NAMES ( MAXBRDF_IDX )

!  Define Kernel names
!  -------------------

      BRDF_CHECK_NAMES(1) = 'Lambertian'
      BRDF_CHECK_NAMES(2) = 'Hapke     '
      BRDF_CHECK_NAMES(3) = 'Rahman    '
      BRDF_CHECK_NAMES(4) = 'Cox-Munk  '

!  Initialize output status
!  ------------------------

      STATUS     = LRRS_SUCCESS
! mick fix 3/31/11 - initialise these
      N_MESSAGES = 0
      MESSAGES   = ' '

!  Expand type structure variables
!  -------------------------------

!  Fixed Boolean Inputs

      DO_RRS_OVERALL         = LRRS_FixIn%Bool%DO_RRS_OVERALL
      DO_DIRECTRRS_ONLY      = LRRS_FixIn%Bool%DO_DIRECTRRS_ONLY
      DO_DIRECT_BEAM         = LRRS_FixIn%Bool%DO_DIRECT_BEAM
      DO_SSCORR_OUTGOING     = LRRS_FixIn%Bool%DO_SSCORR_OUTGOING
      DO_SSFULL              = LRRS_FixIn%Bool%DO_SSFULL
      DO_MOLECSCAT_ONLY      = LRRS_FixIn%Bool%DO_MOLECSCAT_ONLY
      DO_PLANE_PARALLEL      = LRRS_FixIn%Bool%DO_PLANE_PARALLEL
      DO_DELTAM_SCALING      = LRRS_FixIn%Bool%DO_DELTAM_SCALING
      DO_DOUBLE_CONVTEST     = LRRS_FixIn%Bool%DO_DOUBLE_CONVTEST
      DO_UPWELLING           = LRRS_FixIn%Bool%DO_UPWELLING
      DO_DNWELLING           = LRRS_FixIn%Bool%DO_DNWELLING
      DO_LBOUNDARIES         = LRRS_FixIn%Bool%DO_LBOUNDARIES
      DO_MVOUT_ONLY          = LRRS_FixIn%Bool%DO_MVOUT_ONLY
      DO_ADDITIONAL_MVOUT    = LRRS_FixIn%Bool%DO_ADDITIONAL_MVOUT
      DO_LAMBERTIAN_SURFACE  = LRRS_FixIn%Bool%DO_LAMBERTIAN_SURFACE
      DO_WATERLEAVING        = LRRS_FixIn%Bool%DO_WATERLEAVING
      DO_SHADOW_EFFECT       = LRRS_FixIn%Bool%DO_SHADOW_EFFECT
      DO_GLITTER_DBMS        = LRRS_FixIn%Bool%DO_GLITTER_DBMS
      DO_LAMBERTIAN_FRACTION = LRRS_FixIn%Bool%DO_LAMBERTIAN_FRACTION
      DO_LRRS_WRITESCENARIO  = LRRS_FixIn%Bool%DO_LRRS_WRITESCENARIO
      DO_LRRS_WRITERESULTS   = LRRS_FixIn%Bool%DO_LRRS_WRITERESULTS
      DO_LRRS_WRITEINPUT     = LRRS_FixIn%Bool%DO_LRRS_WRITEINPUT
      DO_LRRS_WRITEFOURIER   = LRRS_FixIn%Bool%DO_LRRS_WRITEFOURIER

!  Modified Boolean Inputs

      DO_ELASTIC_ONLY     = LRRS_ModIn%MBool%DO_ELASTIC_ONLY
      DO_CABANNES_RAMAN   = LRRS_ModIn%MBool%DO_CABANNES_RAMAN
      DO_ENERGY_BALANCING = LRRS_ModIn%MBool%DO_ENERGY_BALANCING
      DO_MSMODE_LRRS      = LRRS_ModIn%MBool%DO_MSMODE_LRRS
      DO_BIN_REALIZATION  = LRRS_ModIn%MBool%DO_BIN_REALIZATION
      DO_MONO_REALIZATION = LRRS_ModIn%MBool%DO_MONO_REALIZATION
      DO_USER_STREAMS     = LRRS_ModIn%MBool%DO_USER_STREAMS
      DO_NO_AZIMUTH       = LRRS_ModIn%MBool%DO_NO_AZIMUTH

!  Check inputs (both file-read and derived)
!  -----------------------------------------

!  Check consistency of mean value input control

      IF ( DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
              'Bad input: Cannot have both mean-value flags set'
        ENDIF
      ENDIF

      IF ( .NOT.DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          IF ( .NOT. DO_NO_AZIMUTH ) THEN
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
              'Bad input: Mean-value option requires NO azimuth'
          ENDIF
          IF ( DO_USER_STREAMS ) THEN
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
             'Bad input: Mean-value option needs quadratures only'
          ENDIF
        ENDIF
      ENDIF

!  Check directional input

      IF ( .NOT.DO_UPWELLING .AND. .NOT. DO_DNWELLING ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = &
               'Bad input: no directional input is set'
      ENDIF

!  Check Rayleigh-only input against single scatter correction
!    This is disabled now. Wont make any difference with the Nadir SS
!      IF ( DO_RAYLEIGH_ONLY .AND. DO_SSCORRECTION ) THEN
!        N_MESSAGES = N_MESSAGES + 1
!        MESSAGES(N_MESSAGES) =
!     & 'Bad input: SS correction must be turned off, Rayleigh only'
!      ENDIF

!  Check single scattering correction and Do no Azimuth
!    ---WARNING. Do-no-Azimuth Flag turned off

      IF ( DO_SSCORR_OUTGOING ) THEN
        IF ( DO_NO_AZIMUTH ) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
             'Bad input: need azimuth dependence for SS corrections'
        ENDIF
      ENDIF

!  Check azimuth-only conditions
!  =============================

!  Check no-Azimuth flag

      IF ( .NOT.DO_NO_AZIMUTH ) THEN
        IF ( DO_USER_STREAMS .AND. N_USER_STREAMS.EQ.1 ) THEN
          IF ( USER_ANGLES(1) .EQ. ZERO ) THEN
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
          'Bad input: single zenith-sky output requires no azimuth'
          ENDIF
        ENDIF
      ENDIF

!  Check viewing geometry input
!  ============================

!  Check earth radius (Chapman function only)

      IF ( EARTH_RADIUS.LT.6320.0D0 .OR. &
             EARTH_RADIUS.GT.6420.0D0 ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = &
          'Bad input: Earth radius not in range [6320,6420]'
       ENDIF

!  Check solar zenith angle input

      IF ( SOLAR_ANGLE.LT.ZERO .OR.SOLAR_ANGLE.GE.90.0 ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = &
           ' Bad input: Sun angle not in range [0,90)'
       ENDIF

!  Check relative azimuths

      LOOP = .TRUE.
      I = 0
       DO WHILE (LOOP .AND. I.LT.N_USER_RELAZMS)
        I = I + 1
        IF ( USER_RELAZMS(I) .GT. 360.0   .OR. &
               USER_RELAZMS(I) .LT. ZERO ) THEN
          WRITE(C2,'(I2)')I
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
           'Bad input: out-of-range azimuth angle, no. '//C2
          LOOP = .FALSE.
        ENDIF
      ENDDO

!  Check user-defined stream angles (should always be [0,90])

      IF ( DO_USER_STREAMS ) THEN
        LOOP = .TRUE.
        I = 0
         DO WHILE (LOOP .AND. I.LT.N_USER_STREAMS)
          I = I + 1
          IF ( USER_ANGLES(I) .GT. 90.0   .OR. &
                 USER_ANGLES(I) .LT. ZERO ) THEN
            WRITE(C2,'(I2)')I
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
             'Bad input: out-of-range user stream, no. '//C2
            LOOP = .FALSE.
          ENDIF
        ENDDO
      ENDIF

!  Check height grid input (Chapman function only)

      LOOP = .TRUE.
      I = 0
       DO WHILE (LOOP .AND. I.LT.NLAYERS)
        I = I + 1
        IF ( HEIGHT_GRID(I-1).LE.HEIGHT_GRID(I) ) THEN
          WRITE(C2,'(I2)')I
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
      'Bad input: Height-grid not monotonically decreasing; Layer '//C2
          LOOP = .FALSE.
        ENDIF
      ENDDO

!  Check layer output choices
!  ==========================

!  Check layer boundary choices (should always be within atmosphere!)

      IF ( DO_LBOUNDARIES ) THEN

!  Number of choices should not exceed number of levels in atmosphere
!   0 = TOA, 1 = bottom of first layer, NLAYERS = surface)

        IF ( N_LOUTPUT .GT. NLAYERS + 1 ) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
              'Bad input: Too many Layer boundary outputs'
        ENDIF

!  Choices should remain in atmosphere, i.e. no values > NLAYERS

        DO I = 1, N_LOUTPUT
          IF ( LBOUNDARIES_OUTPUT(I) .GT. NLAYERS ) THEN
            WRITE(C2,'(I2)')LBOUNDARIES_OUTPUT(I)
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
              'Bad input: Layer boundary level > NLAYERS: '//C2
            LOOP = .FALSE.
          ENDIF
        ENDDO
      ENDIF

!  Check repetition of layer boundary output choices

      IF ( DO_LBOUNDARIES ) THEN
        UTA = 0
        LOOP = .TRUE.
        DO WHILE (LOOP .AND. UTA .LT. N_LOUTPUT )
          UTA = UTA + 1
          N1 = LBOUNDARIES_OUTPUT(UTA)
          NSTART = 0
          DO N = 1, N_LOUTPUT
            IF ( LBOUNDARIES_OUTPUT(N) .EQ. N1 ) NSTART = NSTART + 1
          ENDDO
          IF ( NSTART .GT. 1 ) THEN
            LOOP = .FALSE.
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
              'Bad input: repetition of optical depth input value'
          ENDIF
        ENDDO
      ENDIF

!  Check off-boundary choices of output levels
!     Cannot be greater than NLAYERS, cannot be negative

      IF ( .NOT.DO_LBOUNDARIES ) THEN
        DO I = 1, N_LOUTPUT
          IF ( LPARTIALS_OUTPUT(I) .GT. DBLE(NLAYERS) ) THEN
            WRITE(C2,'(I2)')I
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
         'Bad input: off-boundary level > DBLE(NLAYERS), entry #'//C2
            LOOP = .FALSE.
          ELSE IF ( LPARTIALS_OUTPUT(I) .LT. ZERO ) THEN
            WRITE(C2,'(I2)')I
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
              'Bad input: off-boundary level < 0.0, entry #'//C2
            LOOP = .FALSE.
          ENDIF
        ENDDO
      ENDIF

!  Check repetition of user-defined optical depth input

      IF ( .NOT.DO_LBOUNDARIES ) THEN
        UTA = 0
        LOOP = .TRUE.
        DO WHILE (LOOP .AND. UTA .LT. N_LOUTPUT )
          UTA = UTA + 1
          D1 = LPARTIALS_OUTPUT(UTA)
          NSTART = 0
          DO N = 1, N_LOUTPUT
            IF ( LPARTIALS_OUTPUT(N) .EQ. D1 ) NSTART = NSTART + 1
          ENDDO
          IF ( NSTART .GT. 1 ) THEN
            LOOP = .FALSE.
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
            'Bad input: repetition of off-boundary output choice'
          ENDIF
        ENDDO
      ENDIF

!  Check surface inputs
!  ====================

!  Check DB correction only for non-Lambertian surface
!    ---WARNING. Flag turned off. NOw disabled. DB term is automatic
!      IF ( DO_DB_CORRECTION ) THEN
!        IF ( DO_LAMBERTIAN_SURFACE ) THEN
!          N_MESSAGES = N_MESSAGES + 1
!          MESSAGES(N_MESSAGES) =
!     &          'Bad input: No DB correction for Lambertian only'
!          DO_DB_CORRECTION = .FALSE.
!        ENDIF
!      ENDIF

!  Bookkeeping for surface kernel Cox-Munk types

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
        IF ( BRDF_NAMES .EQ. 'Cox-Munk  ' ) THEN
         BRDF_NPARS = 3
         IF ( DO_SHADOW_EFFECT ) THEN
            BRDF_PARS(3) = ONE
          ELSE
            BRDF_PARS(3) = ZERO
          ENDIF
        ENDIF
      ENDIF

!  Check list of names (non-Lambertian)

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
        IF ( WHICH_BRDF.GT.MAXBRDF_IDX.OR.WHICH_BRDF.LE.0) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
                'Bad input: BRDF Index not on list of indices'
        ELSE
          IF ( BRDF_NAMES.NE.BRDF_CHECK_NAMES(WHICH_BRDF) ) THEN
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = &
            'Bad input: BRDF kernel name not corresponding'
          ENDIF
        ENDIF
      ENDIF

!  Check water-leaving flag is only operational for Cox-Munk

      IF ( BRDF_NAMES .NE. 'Cox-Munk  '.AND.DO_WATERLEAVING ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = &
            'Bad input: Only Water-leaving with Cox-Munk kernel'
        DO_WATERLEAVING = .FALSE.
      ENDIF

!  Check on Lambertian fraction (if called for the BRDF surface)
!    --- Not required for the Water-leaving option

      IF ( DO_LAMBERTIAN_FRACTION ) THEN
       IF ( .NOT. DO_WATERLEAVING ) THEN
        IF ( LAMBERTIAN_FRACTION.GT.ONE .OR. &
             LAMBERTIAN_FRACTION.LT.ZERO ) THEN
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES) = &
           'Bad input: Lambertian fraction not in range [0,1]'
        ENDIF
       ENDIF
      ENDIF

!  Set status index

      IF ( N_MESSAGES .NE. 0 ) STATUS = LRRS_SERIOUS

!  Finish

      RETURN
      END SUBROUTINE LRRS_CHECK_INPUTS

!

      SUBROUTINE LRRS_DERIVE_INPUTS &
        ( DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, &
          DO_BIN_REALIZATION, DO_MOLECSCAT_ONLY, DO_LBOUNDARIES, &
          NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS, N_LOUTPUT, &
          NLAYERS, NPOINTS_INNER, SOLAR_ANGLE, USER_ANGLES, &
          NMOMENTS, N_OUT_STREAMS, LBOUNDARIES_OUTPUT, LPARTIALS_OUTPUT, &
          FLUX_FACTOR, COS_SZA, LEVELMASK_UP, LEVELMASK_DN, &
          STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, N_CONVTESTS, &
          QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS )

!  Derivation of bookkeeping inputs

!  include file of dimensions and numbers

      USE LRRS_PARS
      USE LRRS_AUX2

      IMPLICIT NONE

!  Input arguments
!  ===============

!  flags

      LOGICAL, INTENT(IN) :: &
         DO_BIN_REALIZATION,  DO_MOLECSCAT_ONLY,    DO_LBOUNDARIES, &
         DO_UPWELLING,        DO_DNWELLING,         DO_USER_STREAMS

!  Number of layers

      INTEGER, INTENT(IN) ::   NLAYERS

!  Number of discrete ordinate streams

      INTEGER, INTENT(IN) ::   NSTREAMS

!  Solar beam input geometry

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLE

!  User stream variables

      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_ANGLES  ( MAX_USER_STREAMS )

!  Number of azimuths

      INTEGER, INTENT(IN) ::   N_USER_RELAZMS

!    N_LOUTPUT = number of level output choices (all)

      INTEGER, INTENT(IN) ::   N_LOUTPUT

!  Number of inner-window points

      INTEGER, INTENT(IN) ::   NPOINTS_INNER

!  Output arguments
!  ================

!  number of moments in the multiple scatter calculation

      INTEGER, INTENT(OUT) ::   NMOMENTS

!  number of output streams (quadrature or user)

      INTEGER, INTENT(OUT) ::   N_OUT_STREAMS

!  number of convergence tests
!    Added by R. Spurr, 18 November 2005

      INTEGER, INTENT(OUT) ::   N_CONVTESTS

!  Cosein of the solar zenith angle

      REAL(FPK), INTENT(OUT) :: COS_SZA

!  User polar cosines

      REAL(FPK), INTENT(OUT) :: USER_STREAMS ( MAX_USER_STREAMS )

!  Quadrature streams and weights

      REAL(FPK), INTENT(OUT) :: QUAD_STREAMS  ( MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: QUAD_WEIGHTS  ( MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: QUAD_STRMWGT  ( MAX_STREAMS )

!    LBOUNDARIES_OUTPUT = layer boundary choices

      INTEGER, INTENT(INOUT) ::   LBOUNDARIES_OUTPUT ( MAX_LOUTPUT )
      REAL(FPK), INTENT(INOUT) :: LPARTIALS_OUTPUT ( MAX_LOUTPUT )

!  output optical depth masks and indices

      INTEGER, INTENT(OUT) ::     LEVELMASK_UP ( MAX_LOUTPUT )
      INTEGER, INTENT(OUT) ::     LEVELMASK_DN ( MAX_LOUTPUT )

!  Layer masks for doing integrated source terms

      LOGICAL, INTENT(OUT) ::     STERM_LAYERMASK_UP ( MAX_LAYERS )
      LOGICAL, INTENT(OUT) ::     STERM_LAYERMASK_DN ( MAX_LAYERS )

!  Flux factor

      REAL(FPK), INTENT(INOUT) :: FLUX_FACTOR

!  Bookkeeping variables not passed
!  --------------------------------

!  output optical depth masks and indices

      LOGICAL ::          PARTIALS_OUTFLAG  (MAX_LOUTPUT)
      INTEGER ::          PARTIALS_OUTINDEX (MAX_LOUTPUT)
      INTEGER ::          PARTIALS_LEVELMASK_UP (MAX_LOUTPUT)
      INTEGER ::          PARTIALS_LEVELMASK_DN (MAX_LOUTPUT)

!  off-grid optical depths (values, masks, indices)

      INTEGER ::          N_PARTIALS
      INTEGER ::          PARTIALS_LAYERIDX (MAX_PARTIALS_LOUTPUT)
      REAL(FPK) ::        PARTIALS_VALUES   (MAX_PARTIALS_LOUTPUT)

!  number of whole layer source terms required

      INTEGER ::          N_LAYERSOURCE_UP
      INTEGER ::          N_LAYERSOURCE_DN
      INTEGER ::          N_ALLLAYERS_UP
      INTEGER ::          N_ALLLAYERS_DN

!  local variables
!  ---------------

      INTEGER ::          UTA, N, N1, I, NT, UT, UT1, UM
      REAL(FPK) ::        REM, DT
      INTEGER ::          N_DIRECTIONS

!  Get derived inputs (for Bookkeeping)
!  ====================================

!  Flux factor is always one

      FLUX_FACTOR    = ONE

!  solar zenith cosine

      COS_SZA = DCOS(SOLAR_ANGLE*DEG_TO_RAD)

!  Number of moments

      NMOMENTS = 2 * NSTREAMS - 1
      IF ( DO_MOLECSCAT_ONLY )  NMOMENTS = 2

!  quadrature

      IF ( NSTREAMS .EQ. 1 ) THEN
        QUAD_STREAMS(1) = DSQRT ( ONE / THREE )
        QUAD_WEIGHTS(1) = ONE
      ELSE
        CALL GAULEG(ZERO,ONE,QUAD_STREAMS,QUAD_WEIGHTS,NSTREAMS)
      ENDIF

      DO I = 1, NSTREAMS
        QUAD_STRMWGT(I) = QUAD_WEIGHTS(I) * QUAD_STREAMS(I)
      ENDDO

!  number of output streams

      IF ( DO_USER_STREAMS ) THEN
        N_OUT_STREAMS = N_USER_STREAMS
      ELSE
        N_OUT_STREAMS = NSTREAMS
      ENDIF

!  Set up the streams and sines

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          USER_STREAMS(UM) = DCOS(USER_ANGLES(UM)*DEG_TO_RAD)
        ENDDO
      ENDIF

!  number of tests to be applied for convergence

      N_DIRECTIONS = 0
      IF ( DO_UPWELLING ) N_DIRECTIONS = N_DIRECTIONS + 1
      IF ( DO_DNWELLING ) N_DIRECTIONS = N_DIRECTIONS + 1
      N_CONVTESTS = N_USER_RELAZMS * N_USER_STREAMS * N_DIRECTIONS
      N_CONVTESTS = N_CONVTESTS * N_LOUTPUT

      IF ( DO_BIN_REALIZATION ) THEN
        N_CONVTESTS = N_CONVTESTS * NPOINTS_INNER
      ENDIF

!        write(*,*)nPOINTS_INNER,N_DIRECTIONS,
!     &    N_USER_RELAZMS,N_USER_STREAMS, N_LOUTPUT

!  Optical depth control (layer boundaries only)
!  =============================================

      IF ( DO_LBOUNDARIES ) THEN

!  Sort out User optical depths
!  ----------------------------

!  assign

        N_PARTIALS = 0
        DO UTA = 1, N_LOUTPUT
          PARTIALS_OUTFLAG(UTA) = .FALSE.
          LPARTIALS_OUTPUT(UTA) = -999.0D0
        ENDDO

!  Sort in ascending order

        IF ( N_LOUTPUT .GT. 1 ) THEN
          CALL IPSORT(N_LOUTPUT,LBOUNDARIES_OUTPUT)
!mick fix 4/4/11 - added if around do loop
          IF (LBOUNDARIES_OUTPUT(1) > LBOUNDARIES_OUTPUT(N_LOUTPUT)) THEN
            DO UTA = 1, N_LOUTPUT/2
              UT =  LBOUNDARIES_OUTPUT(UTA)
              UT1 = LBOUNDARIES_OUTPUT(N_LOUTPUT + 1 - UTA)
              LBOUNDARIES_OUTPUT(UTA) = UT1
              LBOUNDARIES_OUTPUT(N_LOUTPUT + 1 - UTA) = UT
            ENDDO
          END IF
        ENDIF

!  set optical depth mask (for values at layer boundaries)

        DO UTA = 1, N_LOUTPUT
          N1 = LBOUNDARIES_OUTPUT(UTA)
          PARTIALS_LEVELMASK_UP(UTA) = N1
          PARTIALS_LEVELMASK_DN(UTA) = N1
        ENDDO

!  set optical depth mask (for values at layer boundaries)

        DO UTA = 1, N_LOUTPUT
          N1 = LBOUNDARIES_OUTPUT(UTA)
          LEVELMASK_UP(UTA) = N1
          LEVELMASK_DN(UTA) = N1
        ENDDO

!  Partial layer distance-related input

      ELSE

!  Sort in ascending order

        IF ( N_LOUTPUT .GT. 1 ) THEN
          CALL HPSORT(N_LOUTPUT,LPARTIALS_OUTPUT(1:N_LOUTPUT))
        ENDIF

!  mark all optical depths not equal to layer boundary values

        UT = 0
        DO UTA = 1, N_LOUTPUT
          DT = LPARTIALS_OUTPUT(UTA)
          NT = INT(DT)
          REM = DT - DBLE(NT)
          IF ( REM .NE. ZERO ) THEN
            UT = UT + 1
            PARTIALS_OUTFLAG(UTA)      = .TRUE.
            PARTIALS_OUTINDEX(UTA)     = UT
            PARTIALS_LAYERIDX(UT)      = NT + 1
            PARTIALS_LEVELMASK_UP(UTA) = NT + 1
            PARTIALS_LEVELMASK_DN(UTA) = NT
            PARTIALS_VALUES(UT)        = REM
          ELSE
            PARTIALS_OUTFLAG(UTA)      = .FALSE.
            PARTIALS_OUTINDEX(UTA)     = 0
            PARTIALS_LEVELMASK_UP(UTA) = NT
            PARTIALS_LEVELMASK_DN(UTA) = NT
          ENDIF
        ENDDO
        N_PARTIALS = UT

      ENDIF

!      write(*,*)n_offgrid_usertaus
!      do ut = 1, n_offgrid_usertaus
!       write(*,*)ut, PARTIALS_LAYERIDX(UT),
!     &        PARTIALS_VALUES(UT)
!      enddo
!      pause

!  Set masking and number of layer source terms
!  --------------------------------------------

!   .. for upwelling and downwelling

!mick fix 4/6/12 - initialize STERM_LAYERMASK_UP outside if block
      DO N = 1, NLAYERS
        STERM_LAYERMASK_UP(N) = .FALSE.
      ENDDO
      IF ( DO_UPWELLING ) THEN
        !DO N = 1, NLAYERS
        !  STERM_LAYERMASK_UP(N) = .FALSE.
        !ENDDO

!mick fix 4/5/11 - substituted two lines with if block
        !N_LAYERSOURCE_UP = LEVELMASK_UP(1) + 1
        !N_ALLLAYERS_UP   = N_LAYERSOURCE_UP
        IF ( .NOT. PARTIALS_OUTFLAG(1) ) THEN
          N_ALLLAYERS_UP = LEVELMASK_UP(1) + 1
        ELSE
          N_ALLLAYERS_UP = PARTIALS_LAYERIDX(1)
        ENDIF

        DO N = NLAYERS, N_ALLLAYERS_UP, -1
          STERM_LAYERMASK_UP(N) = .TRUE.
        ENDDO
      ENDIF

!mick fix 4/6/12 - initialize STERM_LAYERMASK_DN outside if block
      DO N = 1, NLAYERS
        STERM_LAYERMASK_DN(N) = .FALSE.
      ENDDO
      IF ( DO_DNWELLING ) THEN
        !DO N = 1, NLAYERS
        !  STERM_LAYERMASK_DN(N) = .FALSE.
        !ENDDO

!mick fix 4/5/11 - substituted two lines with if block
        !N_LAYERSOURCE_DN = LEVELMASK_DN(N_LOUTPUT)
        !N_ALLLAYERS_DN   = N_LAYERSOURCE_DN
        IF ( .NOT. PARTIALS_OUTFLAG(N_LOUTPUT) ) THEN
          N_ALLLAYERS_DN = LEVELMASK_DN(N_LOUTPUT)
        ELSE
          N_ALLLAYERS_DN = PARTIALS_LAYERIDX(N_PARTIALS)
        ENDIF

        DO N = 1, N_ALLLAYERS_DN
          STERM_LAYERMASK_DN(N) = .TRUE.
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LRRS_DERIVE_INPUTS

      END MODULE lrrs_io_check

