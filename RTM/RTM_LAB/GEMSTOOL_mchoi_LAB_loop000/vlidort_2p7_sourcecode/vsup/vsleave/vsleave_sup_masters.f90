! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6, 2.7                             #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            VSLEAVE_INPUTMASTER                              #
! #            VSLEAVE_MAINMASTER (master)                      #
! #                                                             #
! ###############################################################

      MODULE vsleave_sup_masters_m

      PRIVATE
      PUBLIC :: VSLEAVE_INPUTMASTER,&
                VSLEAVE_MAINMASTER

      CONTAINS

      SUBROUTINE VSLEAVE_INPUTMASTER ( &
        FILNAM, VSLEAVE_Sup_In, &
        VSLEAVE_Sup_InputStatus )

!  Input routine for VSLEAVE program

!  Version 2.6 Notes
!  -----------------

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 December 2012. 
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 2.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 2.6.
!  For Version 2.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the VBRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the VBRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      USE VLIDORT_PARS
      USE VSLEAVE_FINDPAR_M

      USE vsleave_sup_inputs_def
      USE vsleave_sup_outputs_def

      IMPLICIT NONE

!  Arguments
!  ---------

      CHARACTER (LEN=*), INTENT(IN) :: FILNAM

      TYPE(VSLEAVE_Sup_inputs), INTENT(OUT) :: VSLEAVE_Sup_In

      TYPE(VSLEAVE_Input_Exception_Handling), INTENT(OUT) :: &
        VSLEAVE_Sup_InputStatus

!  Local variables
!  ---------------

!  Inclusion flag (not really necessary, Brian)

      LOGICAL :: DO_SLEAVING

!  Isotropic flag

      LOGICAL :: DO_ISOTROPIC

!  Flo flag

      LOGICAL :: DO_FLUORESCENCE

!   !@@ Solar sources + Observational Geometry flag !@@

      LOGICAL :: DO_SOLAR_SOURCES
      LOGICAL :: DO_USER_OBSGEOMS

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: DO_EXACT
      LOGICAL :: DO_EXACTONLY

!  Geometry and control
!  --------------------

!  Stream angle flag

      LOGICAL :: DO_USER_STREAMS

!  Number of Stokes components

      INTEGER :: NSTOKES

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  Local angle control

      INTEGER :: NBEAMS
      INTEGER :: N_USER_STREAMS
      INTEGER :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: BEAM_SZAS   (MAXBEAMS)
      REAL(fpk) :: USER_RELAZMS(MAX_USER_RELAZMS)
      REAL(fpk) :: USER_ANGLES (MAX_USER_STREAMS)

!  !@@ Local Observational Geometry control and angles

      INTEGER   :: N_USER_OBSGEOMS
      REAL(fpk) :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk) :: WAVELENGTH

!  Changed for Version 2.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Removed, Version 2.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: NSTREAMS_AZQUAD

!  New for Version 2.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude 

!  Input Epoch

      INTEGER   :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk)  :: FL_Amplitude755

!  Flag for using Data Gaussian parameters

      LOGICAL          :: FL_DO_DataGaussian

!  GEMSTOOL Alteration 11/30/16. Additional linear parameterization variables
!   --> Flag for Linear parameterization, linear parameters

      logical   :: FL_Do_LinearParam
      REAL(FPK) :: FL_ScalingConstant, FL_ScalingGradient

!  Exception handling
!  ------------------

!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  local variables
!  ===============

      CHARACTER (LEN=12), PARAMETER :: PREFIX = 'VSLEAVESUP -'

      LOGICAL ::            ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER ::            I, FILUNIT, NM

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Local error handling initialization

      ERROR  = .FALSE.
      NM     = NMESSAGES

!  Open file

      FILUNIT = VLIDORT_INUNIT
      OPEN(VLIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  Initialize Angle control
!  ========================

      DO_SOLAR_SOURCES = .FALSE.  !@@ New line
      DO_USER_OBSGEOMS = .FALSE.  !@@ New line

      DO_USER_STREAMS = .FALSE.
      NSTREAMS = 0
      NSTOKES = 0

      NBEAMS   = 0
      DO I = 1, MAXBEAMS
        BEAM_SZAS(I) = ZERO
      ENDDO
      N_USER_STREAMS = 0
      DO I = 1, MAX_USER_STREAMS
        USER_ANGLES(I) = ZERO
      ENDDO
      N_USER_RELAZMS = 0
      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO

! !@@ Observational Geometry

      N_USER_OBSGEOMS = 0
      DO I = 1, MAX_USER_OBSGEOMS
        USER_OBSGEOMS(I,1:3) = ZERO
      ENDDO

!  Initialize Surface stuff
!  ========================

!  Control flags

      DO_EXACT        = .FALSE.       !@@  New line
      DO_EXACTONLY    = .FALSE.
      DO_ISOTROPIC    = .FALSE.
      DO_SLEAVING     = .FALSE.
      DO_FLUORESCENCE = .FALSE.

!  Fluorescence variables
!     GEMSTOOL Alteration 11/30/16. Additional linear parameterization variables

      FL_LATITUDE   = ZERO
      FL_LONGITUDE  = ZERO
      FL_EPOCH      = 0
      FL_WAVELENGTH = ZERO

      FL_Amplitude755     = ZERO
      FL_DO_DataGaussian  = .false.

      FL_Do_LinearParam  = .false.
      FL_ScalingConstant = ZERO
      FL_ScalingGradient = ZERO

!  Water-leaving variables

      SALINITY   = ZERO
      CHLORCONC  = ZERO
      WAVELENGTH = ZERO

      WINDSPEED  = ZERO
      WINDDIR    = ZERO

      DO_GlintShadow   = .FALSE.
      DO_FoamOption    = .FALSE.
      DO_FacetIsotropy = .FALSE.

!  Geometry and Input Control
!  ==========================

!  !@@ Solar sources is True, always

      DO_SOLAR_SOURCES = .TRUE.

!  user-defined Stream angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  number of Stokes components

      PAR_STR = 'Number of Stokes vector components'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) NSTOKES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR (ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                          ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
        'Number of half-space streams > maximum dimension'
        ACTIONS(NM)  = &
         'Re-set input value or increase MAXSTREAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  Observational Geometry    !@@
!  ======================

!  !@@ New flag, Observational Geometry

      IF ( DO_SOLAR_SOURCES ) THEN
         PAR_STR = 'Do Observation Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_USER_OBSGEOMS
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  !@@ Observational Geometry control
!     ---- check not exceeding dimensioned number

      IF ( DO_USER_OBSGEOMS ) THEN
        PAR_STR = 'Number of Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_OBSGEOMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of ObsGeometry inputs > Maximum dimension'
          ACTIONS(NM)  = 'Re-set input or increase MAX_USER_OBSGEOMS dimension'
          STATUS       = VLIDORT_SERIOUS
          NMESSAGES    = NM
          GO TO 764
        ENDIF
      ENDIF

!  !@@ Observational Geometry control
!     Automatic setting of NBEAMS, N_USER_STREAMS, N_USER_RELAZMS

      IF ( DO_USER_OBSGEOMS ) THEN
         NBEAMS          = N_USER_OBSGEOMS
         N_USER_STREAMS  = N_USER_OBSGEOMS
         N_USER_RELAZMS  = N_USER_OBSGEOMS
         DO_USER_STREAMS = .TRUE.
      ENDIF

!  !@@ Observational Geometry control

      IF ( DO_USER_OBSGEOMS ) THEN
        PAR_STR = 'Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           DO I = 1, N_USER_OBSGEOMS
             READ (FILUNIT,*,ERR=998)&
               USER_OBSGEOMS(I,1), USER_OBSGEOMS(I,2), USER_OBSGEOMS(I,3)
           ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Set angles

      IF ( DO_USER_OBSGEOMS ) THEN
         BEAM_SZAS   (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
         USER_ANGLES (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
         USER_RELAZMS(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)
         GO TO 5667
      ENDIF

!  Solar beams
!  ===========

!  Number of Solar zenith angles

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NBEAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Check not exceeding dimensioned number

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
        'Number of solar zenith angles > maximum dimension'
        ACTIONS(NM)  = &
        'Re-set input value or increase MAXBEAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  TOA solar zenith angle inputs

      PAR_STR = 'Solar zenith angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, NBEAMS
          READ (FILUNIT,*,ERR=998) BEAM_SZAS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Azimuth angles
!  ==============

!  Number of azimuth angles

      PAR_STR = 'Number of user-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Check not exceeding dimensioned number

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
         'Number of relative azimuth angles > maximum dimension'
        ACTIONS(NM)  = &
         'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = VLIDORT_SERIOUS
        NMESSAGES    = NM
        GO TO 764
      ENDIF

! Azimuth angles

      PAR_STR = 'User-defined relative azimuth angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_RELAZMS
          READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  User-defined viewing zenith angles (should be positive)
!  ==================================

      IF ( DO_USER_STREAMS ) THEN

!  Number of user-defined viewing zenith angles

        PAR_STR = 'Number of user-defined viewing zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )

!  Check dimension

        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = &
          'Number of viewing zenith angles > maximum dimension'
          ACTIONS(NM)  = &
          'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  User-defined viewing zenith angles

        PAR_STR = 'User-defined viewing zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )

      ENDIF

!  !@@ Continuation point for Skipping the Lattice-input angles

5667  continue

!  Surface stuff
!  =============

!  SLEAVING input
!  --------------

!  Basic flag

      PAR_STR = 'Do surface-leaving Contributions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SLEAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Skip if not doing surface leaving
!   @@ Rob Fix 04 August 2014

      if ( .not.DO_SLEAVING ) GOTO 652

!  Isotropic flag

      PAR_STR = 'Do Isotropic surface-leaving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) DO_ISOTROPIC
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                         ACTIONS )

!  !@@ Overall-Exact flag

      PAR_STR = 'Do Overall-Exact surface-leaving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
         READ (FILUNIT,*,ERR=998)DO_EXACT
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Exact only flag. Only if above is set (!@@)

      IF ( DO_EXACT ) THEN
        PAR_STR = 'Do Exact-only (no Fourier-term contributions)?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_EXACTONLY
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Basic source

      PAR_STR = 'Do surface-leaving Fluorescence?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_FLUORESCENCE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Continuation point

652   continue

!  Inputs for Water-leaving (Non-Fluorescence case)
!  ------------------------------------------------

      IF ( DO_SLEAVING.and..not.DO_FLUORESCENCE ) THEN

!  salinity, chlorophyll concentration, wavelength

        PAR_STR = 'Ocean water salinity [ppt]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) SALINITY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Chlorophyll concentration in [mg/M]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) CHLORCONC
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Wavelength in [Microns]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) WAVELENGTH
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  New for Version 2.7.   Wind-speed and directions, used for the first time

        PAR_STR = 'Windspeed in [m/s]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) WINDSPEED
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Wind directions (degrees) relative to sun positions'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, NBEAMS
            READ (FILUNIT,*,ERR=998) WINDDIR(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  New for Version 2.7. Control Flag for including Whitecaps

      PAR_STR = 'Do whitecap (foam) calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) Do_FoamOption
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  New for Version 2.7. Glint, Control Flags for Facet Isotropy and Shadowing
!     These are only needed for the non-Isotropic case

      if ( .not. do_Isotropic ) then

         PAR_STR = 'Do glint calculation with facet isotropy?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) Do_FacetIsotropy
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                              ACTIONS )

         PAR_STR = 'Do glint calculation with shadowing?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_GlintShadow
!           READ (FILUNIT,*,ERR=998) Do_FacetIsotropy ! Bug Fix 9/27/14
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                              ACTIONS )
      endif

!  Removed for Version 2.7. --> Quadrature is now internal
!    Non-isotropic input = number of azimuth streams, check this value
!        IF ( .not. DO_ISOTROPIC ) THEN
!          PAR_STR = 'Number of azimuth quadrature streams'
!          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!              READ (FILUNIT,*,ERR=998) NSTREAMS_AZQUAD
!          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!          IF ( NSTREAMS_AZQUAD .GT. MAXSTREAMS_BRDF ) THEN
!            NM = NM + 1
!            MESSAGES(NM) = 'Number of AZQUAD streams > maximum dimension'
!            ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS_BRDF dimension'
!            STATUS = VLIDORT_SERIOUS
!            NMESSAGES = NM
!            GO TO 764
!          ENDIF
!        ENDIF

!  Check Added 9/27/14. Multiple SZAS not allowed with Facet ANISOTROPY

        IF ( .not. Do_FacetIsotropy .and. NBEAMS.gt.1 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Facet Anisotropy (wind direction) not Allowed for NBEAMS > 1'
          ACTIONS(NM)  = 'Either re-set NBEAMS = 1 or re-set Do_FacetIsotropy = .true.'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Inputs for Fluorescence Case
!  ----------------------------

      ELSE IF ( DO_SLEAVING.and.DO_FLUORESCENCE ) THEN

!  Temporary Check 

        IF ( .not. DO_ISOTROPIC ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'DO_ISOTROPIC was set to .FALSE. in fluorescence case'
          ACTIONS(NM)  = 'Tempo! Set DO_ISOTROPIC to .TRUE. if doing fluorescence'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  GEMSTOOL Alteration 11/30/16. Additional linear parameterization variables

        PAR_STR = 'Do Linear Parameterization for Fluorescence?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_DO_LinearParam
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        IF ( FL_DO_LinearParam ) then
           PAR_STR = 'Linear Parameterization Fluorescence: Scaling Constant'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) FL_ScalingConstant
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
           PAR_STR = 'Linear Parameterization Fluorescence: Scaling Gradient'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) FL_ScalingGradient
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        ENDIF

!  Use of Data Gaussians (New, 8 August 2012)
!     IF NOT SET, YOU MUST USE YOUR OWN PARAMETERS
!     GEMSTOOL Alteration 11/30/16. Only with Gaussian parameterization

        IF ( .not. FL_DO_LinearParam ) then
           PAR_STR = 'Do Data Gaussians in Fluorescence?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) FL_DO_DataGaussian
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        ENDIF

!  Amplitude for FS755 (Nominally, this is one).
!     GEMSTOOL Alteration 11/30/16. Only with Gaussian parameterization

        IF ( .not. FL_DO_LinearParam ) then
           PAR_STR = 'Amplitude for Fluorescence model at 755 nm'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) FL_Amplitude755
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        ENDIF

!  Lat/Long, day-of-year, wavelength

        PAR_STR = 'Latitude for Fluorescence model [degs]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_LATITUDE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Longitude for Fluorescence model [degs]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_LONGITUDE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Epoch for Fluorescence model'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_EPOCH(1:6)
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Wavelength for Fluorescence model in [nm]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_WAVELENGTH
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

      ENDIF

!  Successful finish

      CLOSE(FILUNIT)

!mick fix
      NMESSAGES = NM

!  Copy Control inputs

      VSLEAVE_Sup_In%SL_DO_USER_STREAMS  = DO_USER_STREAMS
      VSLEAVE_Sup_In%SL_DO_SLEAVING      = DO_SLEAVING
      VSLEAVE_Sup_In%SL_DO_FLUORESCENCE  = DO_FLUORESCENCE
      VSLEAVE_Sup_In%SL_DO_ISOTROPIC     = DO_ISOTROPIC
      VSLEAVE_Sup_In%SL_DO_EXACT         = DO_EXACT         !@@
      VSLEAVE_Sup_In%SL_DO_EXACTONLY     = DO_EXACTONLY
      VSLEAVE_Sup_In%SL_DO_SOLAR_SOURCES = DO_SOLAR_SOURCES   !@@
      VSLEAVE_Sup_In%SL_DO_USER_OBSGEOMS = DO_USER_OBSGEOMS   !@@

!  Copy Geometry results

      VSLEAVE_Sup_In%SL_NSTOKES           = NSTOKES
      VSLEAVE_Sup_In%SL_NSTREAMS          = NSTREAMS
      VSLEAVE_Sup_In%SL_NBEAMS            = NBEAMS
      VSLEAVE_Sup_In%SL_BEAM_SZAS         = BEAM_SZAS
      VSLEAVE_Sup_In%SL_N_USER_RELAZMS    = N_USER_RELAZMS
      VSLEAVE_Sup_In%SL_USER_RELAZMS      = USER_RELAZMS
      VSLEAVE_Sup_In%SL_N_USER_STREAMS    = N_USER_STREAMS
      VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT = USER_ANGLES
      VSLEAVE_Sup_In%SL_N_USER_OBSGEOMS   = N_USER_OBSGEOMS !@@
      VSLEAVE_Sup_In%SL_USER_OBSGEOMS     = USER_OBSGEOMS   !@@

!  Copy Water-leaving inputs
!  -------------------------

!  Original

      VSLEAVE_Sup_In%SL_SALINITY         = SALINITY
      VSLEAVE_Sup_In%SL_CHLORCONC        = CHLORCONC
      VSLEAVE_Sup_In%SL_WAVELENGTH       = WAVELENGTH

!  Version 2.7 changes

!      VSLEAVE_Sup_In%SL_NSTREAMS_AZQUAD  = NSTREAMS_AZQUAD
      VSLEAVE_Sup_In%SL_WINDSPEED        = WINDSPEED
      VSLEAVE_Sup_In%SL_WINDDIR          = WINDDIR

      VSLEAVE_Sup_In%SL_DO_GlintShadow   = DO_GlintShadow
      VSLEAVE_Sup_In%SL_DO_FoamOption    = DO_FoamOption
      VSLEAVE_Sup_In%SL_DO_FacetIsotropy = DO_FacetIsotropy

!  Copy Fluorescence inputs
!  ------------------------

!  GEMSTOOL Alteration 11/30/16. Additional linear parameterization variables

      VSLEAVE_Sup_In%SL_FL_LATITUDE        = FL_LATITUDE
      VSLEAVE_Sup_In%SL_FL_LONGITUDE       = FL_LONGITUDE
      VSLEAVE_Sup_In%SL_FL_WAVELENGTH      = FL_WAVELENGTH
      VSLEAVE_Sup_In%SL_FL_EPOCH           = FL_EPOCH
      VSLEAVE_Sup_In%SL_FL_Amplitude755    = FL_Amplitude755
      VSLEAVE_Sup_In%SL_FL_DO_DataGaussian = FL_DO_DataGaussian

      VSLEAVE_Sup_In%SL_FL_DO_LinearParam  = FL_DO_LinearParam
      VSLEAVE_Sup_In%SL_FL_ScalingConstant = FL_ScalingConstant
      VSLEAVE_Sup_In%SL_FL_ScalingGradient = FL_ScalingGradient

!  Exception handling

      VSLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD = STATUS
      VSLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES   = NMESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTMESSAGES    = MESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTACTIONS     = ACTIONS

!  Normal return

      RETURN

!  Open file error

300   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
      GO TO 764

!  line read error - abort immediately

998   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'read failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Re-set: Entry is incorrect in input file'
      CLOSE(FILUNIT)
      GO TO 764

!  Final error copying

764   CONTINUE

      VSLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD = STATUS
      VSLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES   = NMESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTMESSAGES    = MESSAGES
      VSLEAVE_Sup_InputStatus%SL_INPUTACTIONS     = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE VSLEAVE_INPUTMASTER

!

      SUBROUTINE VSLEAVE_MAINMASTER ( Datapath, &
        VSLEAVE_Sup_In,         & ! Inputs
        VSLEAVE_Sup_Out )         ! Outputs

!  Prepares the Surface Leaving necessary for VLIDORT.

!  Version 2.6 Notes
!  -----------------

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 December 2012. 
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 2.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 2.6.
!  For Version 2.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the VBRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the VBRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      USE VLIDORT_PARS

      USE vsleave_sup_inputs_def
      USE vsleave_sup_outputs_def

      USE vsleave_sup_aux_m      , only : VSleave_GAULEG
      USE vsleave_sup_routines_m , only : WaterLeaving,              &
                                          get_fluorescence_755,      &
                                          solar_spec_irradiance

      IMPLICIT NONE

!  input

      CHARACTER*(*), intent(in) :: Datapath

!  Input structure
!  ---------------

      TYPE(VSLEAVE_Sup_Inputs), INTENT(IN)   :: VSLEAVE_Sup_In

!  Output structure
!  ----------------

      TYPE(VSLEAVE_Sup_Outputs), INTENT(OUT) :: VSLEAVE_Sup_Out

!  VLIDORT local variables
!  +++++++++++++++++++++++

!  Input arguments
!  ===============

!  Main Boolean flags
!  ------------------

!  Inclusion flag (not really necessary, Brian)

      LOGICAL :: DO_SLEAVING

!  Isotropic flag

      LOGICAL :: DO_ISOTROPIC

!  Flo flag

      LOGICAL :: DO_FLUORESCENCE

!  !@@ Solar sources + Observational Geometry flag !@@

      LOGICAL :: DO_SOLAR_SOURCES
      LOGICAL :: DO_USER_OBSGEOMS

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: DO_EXACT
      LOGICAL :: DO_EXACTONLY

!  Geometry and control
!  --------------------

!  Stream angle flag

      LOGICAL :: DO_USER_STREAMS

!  Number of Stokes components

      INTEGER :: NSTOKES

!  Number of discrete ordinate streams, quadrature
!   Version 2.7, added the qudrature arrays

      INTEGER   :: NSTREAMS
      REAL(fpk) :: STREAMS (MAXSTREAMS)
      REAL(fpk) :: WEIGHTS (MAXSTREAMS)

!  Local angle control

      INTEGER   :: NBEAMS
      INTEGER   :: N_USER_STREAMS
      INTEGER   :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: BEAM_SZAS   (MAXBEAMS)
      REAL(fpk) :: USER_RELAZMS(MAX_USER_RELAZMS)
      REAL(fpk) :: USER_ANGLES (MAX_USER_STREAMS)

!  !@@ Local Observational Geometry control and angles

      INTEGER    :: N_USER_OBSGEOMS
      REAL(fpk)  :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk) :: WAVELENGTH

!  Changed for Version 2.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Removed, Version 2.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: NSTREAMS_AZQUAD

!  New for Version 2.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude

!  Input Epoch

      INTEGER :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk)  :: FL_Amplitude755

!  Flag for using Data Gaussians

      LOGICAL :: FL_DO_DataGaussian

!  GEMSTOOL Alteration 11/30/16. Additional linear parameterization variables
!   --> Flag for Linear parameterization, linear parameters

      logical   :: FL_Do_LinearParam
      REAL(FPK) :: FL_ScalingConstant, FL_ScalingGradient

!  Local functions
!  ===============

!  Fluorescence Isotropic Surface leaving term

   REAL(fpk) :: Fluor_ISOTROPIC

!  Isotropic value. Fast calculation

   REAL(fpk)    :: WLeaving_ISO ( MAXBEAMS )

!  Input solar, output stream angles

   REAL(fpk)    :: WLeaving_SD ( MAXBEAMS, MAXSTREAMS )

!  input solar, output view angles

   REAL(fpk)    :: WLeaving_SV ( MAXBEAMS, MAX_USER_STREAMS )

!  Exact Surface-Leaving term
!   REAL(fpk), dimension ( MAXSTOKES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS ) :: SLTERM_USERANGLES
!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams
!   REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, MAXBEAMS )       :: SLTERM_F_0
!   REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, MAXBEAMS ) :: USER_SLTERM_F_0

!  Other local variables
!  =====================

!  (Version 2.6 code). Water-leaving model
!      REAL :: WAV,CHL,RW,SAL,A,REFR,REFI,N12,RWB,TDS,TDV

!  Fluorescence Gaussian parameters
!     Parameters of the fluorescence Gaussian spectral shape model.
!           Gaussian    A (Wm−2 μm−1 sr−1) Lambda(nm) Sigma(nm)
!              1           1.445           736.8        21.2
!              2           0.868           685.2        9.55

      REAL(FPK) :: FL_DataGAUSSIANS(3,2), FL_GAUSSIANS(3,2)
      data FL_DataGAUSSIANS(1,1) / 1.445d0 /
      data FL_DataGAUSSIANS(2,1) / 736.8d0 /
      data FL_DataGAUSSIANS(3,1) / 21.2d0  /
      data FL_DataGAUSSIANS(1,2) / 0.868d0 /
      data FL_DataGAUSSIANS(2,2) / 685.2d0 /
      data FL_DataGAUSSIANS(3,2) / 9.55d0  /

!  Solar spectral radiance model wavelength

      REAL(FPK) :: ssr_wvl

!  Fluorescence model

      CHARACTER*256 :: Fluofile
      INTEGER   :: IB, K, UM
      REAL(FPK) :: Fs755(MAXBEAMS), FL_SunSpec, FsSum
      REAL(FPK) :: ampli, lamda, sigma, arg, gauss, var, ff
      !REAL(FPK) :: solar_spec_irradiance

      INTEGER, PARAMETER :: LUM = 1   !@@
      INTEGER, PARAMETER :: LUA = 1   !@@

!  Copy from input structure
!  -------------------------

!  Copy Top-level general Control inputs

      DO_USER_STREAMS = VSLEAVE_Sup_In%SL_DO_USER_STREAMS
      DO_SLEAVING     = VSLEAVE_Sup_In%SL_DO_SLEAVING

      DO_EXACT        = VSLEAVE_Sup_In%SL_DO_EXACT          !@@
      DO_EXACTONLY    = VSLEAVE_Sup_In%SL_DO_EXACTONLY
      DO_FLUORESCENCE = VSLEAVE_Sup_In%SL_DO_FLUORESCENCE
      DO_ISOTROPIC    = VSLEAVE_Sup_In%SL_DO_ISOTROPIC

!  !@@ New lines

      DO_SOLAR_SOURCES = VSLEAVE_Sup_In%SL_DO_SOLAR_SOURCES
      DO_USER_OBSGEOMS = VSLEAVE_Sup_In%SL_DO_USER_OBSGEOMS

!  Set number of stokes elements and streams
!   Stream Quadrature new for Version 2.7

      NSTOKES  = VSLEAVE_Sup_In%SL_NSTOKES
      NSTREAMS = VSLEAVE_Sup_In%SL_NSTREAMS
      CALL VSleave_Gauleg ( zero, one, STREAMS(1:nstreams), WEIGHTS(1:nstreams), NSTREAMS )

!   !@@ Observational Geometry + Solar sources Optionalities
!   !@@ Either set from User Observational Geometry
!          Or Copy from Usual lattice input

      IF ( DO_USER_OBSGEOMS ) THEN
        N_USER_OBSGEOMS = VSLEAVE_Sup_In%SL_N_USER_OBSGEOMS
        USER_OBSGEOMS   = VSLEAVE_Sup_In%SL_USER_OBSGEOMS
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS          = N_USER_OBSGEOMS
          N_USER_STREAMS  = N_USER_OBSGEOMS
          N_USER_RELAZMS  = N_USER_OBSGEOMS
          BEAM_SZAS   (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
          USER_ANGLES (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
          USER_RELAZMS(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
          N_USER_STREAMS = N_USER_OBSGEOMS
          USER_ANGLES(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        ENDIF
      ELSE
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS         = VSLEAVE_Sup_In%SL_NBEAMS
          BEAM_SZAS      = VSLEAVE_Sup_In%SL_BEAM_SZAS
          N_USER_RELAZMS = VSLEAVE_Sup_In%SL_N_USER_RELAZMS
          USER_RELAZMS   = VSLEAVE_Sup_In%SL_USER_RELAZMS
          N_USER_STREAMS = VSLEAVE_Sup_In%SL_N_USER_STREAMS
          USER_ANGLES    = VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
          N_USER_STREAMS = VSLEAVE_Sup_In%SL_N_USER_STREAMS
          USER_ANGLES    = VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT
        ENDIF
      ENDIF

!  Copy Water-leaving inputs
!  -------------------------

!  Original

      SALINITY        = VSLEAVE_Sup_In%SL_SALINITY
      CHLORCONC       = VSLEAVE_Sup_In%SL_CHLORCONC
      WAVELENGTH      = VSLEAVE_Sup_In%SL_WAVELENGTH

!  Version 2.7 changes

!      NSTREAMS_AZQUAD = VSLEAVE_Sup_In%SL_NSTREAMS_AZQUAD
      WINDSPEED       = VSLEAVE_Sup_In%SL_WINDSPEED
      WINDDIR         = VSLEAVE_Sup_In%SL_WINDDIR

      DO_GlintShadow   = VSLEAVE_Sup_In%SL_DO_GlintShadow
      DO_FoamOption    = VSLEAVE_Sup_In%SL_DO_FoamOption
      DO_FacetIsotropy = VSLEAVE_Sup_In%SL_DO_FacetIsotropy

!  Copy Fluorescence inputs
!  ------------------------

!  Main variables

      FL_Wavelength   = VSLEAVE_Sup_In%SL_FL_Wavelength
      FL_Latitude     = VSLEAVE_Sup_In%SL_FL_Latitude
      FL_Longitude    = VSLEAVE_Sup_In%SL_FL_Longitude
      FL_Epoch        = VSLEAVE_Sup_In%SL_FL_Epoch
      FL_Amplitude755 = VSLEAVE_Sup_In%SL_FL_Amplitude755
      FL_DO_DataGaussian = VSLEAVE_Sup_In%SL_FL_DO_DataGaussian

!  GEMSTOOL Alteration 11/30/16. Copy linear parameterization variables

      FL_Do_LinearParam  = VSLEAVE_Sup_In%SL_FL_DO_LinearParam
      FL_ScalingConstant = VSLEAVE_Sup_In%SL_FL_ScalingConstant
      FL_ScalingGradient = VSLEAVE_Sup_In%SL_FL_ScalingGradient

!mick fix 8/31/2012 - added outer if block
!  GEMSTOOL Alteration 11/30/16. Gaussian parameterization only, otherwise zero array.

      if (DO_FLUORESCENCE) then
        if ( .not. FL_Do_LinearParam ) then
           if ( FL_DO_DataGaussian ) then
              FL_GAUSSIANS(1:3,1) = FL_DataGAUSSIANS(1:3,1)
              FL_GAUSSIANS(1:3,2) = FL_DataGAUSSIANS(1:3,2)
           else
              FL_GAUSSIANS(1:3,1) = VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,1)
              FL_GAUSSIANS(1:3,2) = VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,2)
           endif
         else
           FL_GAUSSIANS = zero
         endif
      endif

!  Main code
!  ---------

!  Zero the output

      VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC  = ZERO
      VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES = ZERO
      VSLEAVE_Sup_Out%SL_SLTERM_F_0        = ZERO
      VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0   = ZERO

!  Return if not called
!   ROb Fix, 04 August 2014

      IF ( .not. DO_SLEAVING ) RETURN

!  Fluorescence
!  ============

      IF ( DO_FLUORESCENCE ) THEN

!  Temporary - Only Isotropic yet.

        IF ( .not.DO_ISOTROPIC ) &
          Stop 'Non-isotropic not allowed yet if doing fluorescence'

!  F_755 data file

!        Fluofile = 'vlidort_v_test/data/fluor_data_2009_fortran.dat'
        Fluofile = Adjustl(Trim(Datapath))//'fluor_data_2009_fortran.dat'

!  Get solar spectral radiance, in (W m−2 μm−1 sr−1), to normalize data

        !FL_SunSpec = 1.0d0  ! Temporary

        ssr_wvl = FL_Wavelength*1.0d-3 !convert from nm to um
!        FL_SunSpec = solar_spec_irradiance( ssr_wvl )
        FL_SunSpec = solar_spec_irradiance( ssr_wvl, Datapath )

!  factor: After  some fiddling, this is 1.0 (July 30th, 2012)
!    4pi factor is required in DB correction ter,

!         FF = PI4
        FF = 1.0d0
!        FF = 0.0d0

!  For each Solar zenith angle

        DO IB = 1, NBEAMS

 !  Get the F_755 data from the subroutine

          CALL get_fluorescence_755 &
   ( FL_Latitude, FL_Longitude, FL_Epoch, BEAM_SZAS(IB), FluoFile, Fs755(IB) )

!  Apply Gaussians, or linear parameterization.
!    GEMSTOOL Alteration 11/30/16. Linear parameterization introduced as alternative..
!    GEMSTOOL Alteration 11/30/16. For Gaussian option, amplitude scaling done here, not later.

          if ( FL_Do_LinearParam ) then
            FsSum = FL_ScalingConstant + FL_ScalingGradient * ( 1.0d0 - ( FL_Wavelength / 755.0_fpk ) )            
          else
            FsSum = zero
            do k = 1, 2
              ampli = FL_Gaussians(1,k)
              lamda = FL_Gaussians(2,k)
              sigma = FL_Gaussians(3,k)
              var = 0.5d0/sigma/sigma
              arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
              Gauss = zero
              if ( arg.lt.88.0d0 ) gauss = ampli * dexp ( - arg )
              FsSum = FsSum + Gauss
            enddo
            FsSum = FsSum * FL_Amplitude755
          endif

!  Assign output Fluorescence: multiply by Fs755, and normalize to solar spectrum
!   FF is the fudge factor

          Fluor_ISOTROPIC = FF * FsSum * Fs755(IB) / FL_SunSpec
          VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1,IB) = Fluor_ISOTROPIC

!   GEMSTOOL Alteration 11/30/16. Fill in the other variables

          IF ( DO_EXACT .or. DO_EXACTONLY ) then
             VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1,1:n_user_streams,1:n_user_relazms,IB) = Fluor_ISOTROPIC
          ENDIF
          if ( .not. DO_EXACTONLY ) then
             VSLEAVE_Sup_Out%SL_SLTERM_F_0      (0,1,1:nstreams,      ib) = Fluor_ISOTROPIC
             VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0 (0,1,1:n_user_streams,ib) = Fluor_ISOTROPIC
          endif

!  End Beam loop

        ENDDO

      ENDIF

!  WATER-LEAVING
!  =============

      IF ( .not. DO_FLUORESCENCE ) THEN

!  Call the routine for Version 2.7

         Call WaterLeaving &
         ( Maxbeams, Max_User_streams, Maxstreams,                            &
           do_Isotropic, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,     &
           Wavelength, Salinity, ChlorConc, Windspeed, WindDir,               &
           nbeams, n_user_streams, nstreams, beam_szas, user_angles, streams, &
           WLeaving_ISO, WLeaving_SD, WLeaving_SV )

!  Copy to Type structure arrays
!    If Isotropic, VLIDORT takes care of the rest - No, this has been changed, 11/30/16
!   GEMSTOOL Alteration 11/30/16. Fill in the other variables, more carefully

         do ib = 1, nbeams
            VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1,ib) = WLeaving_ISO(ib)
            IF ( DO_EXACT .or. DO_EXACTONLY ) then
               if ( do_Isotropic ) then
                  VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1,1:n_user_streams,1:n_user_relazms,IB) = WLeaving_ISO(ib)
               else
                  do um = 1, n_user_streams
                     VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1,um,1:n_user_relazms,ib) = WLeaving_SV(ib,um)
                  enddo
               endif
            ENDIF
            if ( .not. DO_EXACTONLY ) then
               if ( do_Isotropic ) then
                  VSLEAVE_Sup_Out%SL_SLTERM_F_0      (0,1,1:nstreams,      ib) = WLeaving_ISO(ib)
                  VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0 (0,1,1:n_user_streams,ib) = WLeaving_ISO(ib)
               else
                  do um = 1, n_user_streams
                     VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0 (0,1,um,ib) =  WLeaving_SV(ib,um)
                  enddo
                  do um = 1, nstreams
                     VSLEAVE_Sup_Out%SL_SLTERM_F_0(0,1,um,ib) = WLeaving_SD(ib,um)
                  enddo
               endif
            endif
         enddo

!  Old code, pre 11/30/16
!         VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1,1:NBEAMS) = WLeaving_ISO(1:NBEAMS)
!         if ( .not. do_Isotropic ) then
!            do ib = 1, nbeams
!               do um = 1, n_user_streams
!                  VSLEAVE_Sup_Out%SL_SLTERM_USERANGLES(1,um,1:n_user_relazms,ib) = WLeaving_SV(ib,um)
!                  VSLEAVE_Sup_Out%SL_USER_SLTERM_F_0  (0,1,um,ib)                = WLeaving_SV(ib,um)
!               enddo
!               VSLEAVE_Sup_Out%SL_SLTERM_F_0(0,1,1:nstreams,ib) = WLeaving_SD(ib,1:nstreams)
!            enddo
!         endif

!  Here is the compressed Version 2.6 code ---------------------------------------
!    INDWAT/MORCASIWAT calls . use single precision routine
!        SAL = REAL(SALINITY) ; WAV = REAL(WAVELENGTH)
!        CALL INDWAT(WAV,SAL,refr,refi)
!        CHL = REAL(CHLORCONC) ; WAV = REAL(WAVELENGTH)
!        CALL MORCASIWAT(WAV,CHL,RW,.false.)
!     Perfect Transmittance. Add change in solid angle through surface
!     that accounts for 1/(n12*n12) decrease in directional reflectivity
!        a   = 0.485 ; tds = 1.0 ; tdv = 1.0
!        n12 = refr*refr + refi*refi  ; n12 = sqrt(n12)
!        Rwb=(1.0/(n12*n12))*tds*tdv*Rw/(1-a*Rw)
!        VSLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1,1:NBEAMS)= DBLE(Rwb)

!  end water leaving

      endif

!  Finish

      RETURN
      END SUBROUTINE VSLEAVE_MAINMASTER

      END MODULE vsleave_sup_masters_m

