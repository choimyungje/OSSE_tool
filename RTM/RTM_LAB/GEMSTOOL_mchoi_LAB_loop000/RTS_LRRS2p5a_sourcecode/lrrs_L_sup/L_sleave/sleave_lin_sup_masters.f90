! 
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

!  Developed for LRRS Version 2.5, 9/8/15. R. Spurr, RT SOLUTIONS Inc.
!   Closely follows the LIDORT module with the same name, but
!     (1) No solar angle dimensioning, No observational geometry.
!     (2) Additional wavelength dimensioning (MAX_SLEAVE_POINTS)
!     (3) Additional control for using 1 or all wavelengths

!  Note (9/8/15). The number of wavelengths for SLEAVE has deliberately
!  been left flexible - typically the SLEAVE properties will change very
!  little over a Raman-scattering window (+/- 2 nm in the UV), so to
!  a very good approximation, it is often sufficient to use one point for
!  the calculations of SLEAVE - in which case, MAX_SLEAVE_POINTS will be 1
!  This applies equally to Waterleaving and Fluorescence, where wavelength
!  variation can be important in the Visible and NIR.

!  Now, whenever the SLEAVE supplement is used with LRRS, the choice of 
!  SLEAVE wavelengths is linked to the Raman wavelengths (LAMDAS_RANKED).
!  These wavelengths are now input to the SLEAVE supplement MASTER, and
!  they are not set by hand or by configuration-file read. 

!  However, there is a hard-wired choice (DO_WAV1) which allows you
!  to choose a single wavelength point for calculation. When this flag is
!  set, the number of SLEAVE wavelengths = 1 (Waterleaving or Fluorescence)
!  and the single wavelength is set to the AVERAGE VALUE of LAMBDAS_RANKED.

!  If you want all wavelengths, then control flag DO_Wav1 is False,
!  and number of SLEAVE wavelengths  N_Lambdas_Ranked, and the
!  wavelengths themselves are just copied from LAMBDAS_RANKED
!  A SLEAVE calculation will be done for all Raman-scattered wavelengths!!!

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            SLEAVE_LIN_INPUTMASTER                           #
! #            SLEAVE_LIN_MAINMASTER (master)                   #
! #                                                             #
! ###############################################################

      MODULE sleave_linsup_masters_m

      PRIVATE
      PUBLIC :: SLEAVE_LIN_INPUTMASTER,&
                SLEAVE_LIN_MAINMASTER

      CONTAINS

      SUBROUTINE SLEAVE_LIN_INPUTMASTER ( &
        FILNAM, SLEAVE_Sup_In, SLEAVE_LinSup_In, &
        SLEAVE_Sup_InputStatus )

!  Input routine for SLEAVE program

      USE LRRS_pars_m
      USE SLEAVE_FINDPAR_m

      USE sleave_sup_inputs_def_m
      USE sleave_linsup_inputs_def_m

      USE sleave_sup_outputs_def_m

      IMPLICIT NONE

!  Arguments
!  ---------

      CHARACTER (LEN=*), INTENT(IN) :: FILNAM

      TYPE(SLEAVE_Sup_inputs)   , INTENT(OUT) :: SLEAVE_Sup_In
      TYPE(SLEAVE_LinSup_inputs), INTENT(OUT) :: SLEAVE_LinSup_In

      TYPE(SLEAVE_Input_Exception_Handling), INTENT(OUT) :: &
        SLEAVE_Sup_InputStatus

!  Local variables
!  ===============

!  Main Boolean flags
!  ------------------

!  Inclusion flag (not really necessary, Brian)

      LOGICAL :: DO_SLEAVING

!  Isotropic flag

      LOGICAL :: DO_ISOTROPIC

!  Flo flag

      LOGICAL :: DO_FLUORESCENCE

!  Single- vs multi-wavelength flag

      LOGICAL :: DO_WAV1

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: DO_EXACT
      LOGICAL :: DO_EXACTONLY

!  Linearizations

      LOGICAL :: DO_SL_JACOBIANS
      LOGICAL :: DO_ISO_JACOBIANS
      INTEGER :: N_SLEAVE_WFS

!  Geometry and control
!  --------------------

!  Stream angle flag

      LOGICAL :: DO_USER_STREAMS

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  Local angle control

      INTEGER :: N_USER_STREAMS
      INTEGER :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: BEAM_SZA
      REAL(fpk) :: USER_RELAZMS(MAX_USER_RELAZMS)
      REAL(fpk) :: USER_ANGLES (MAX_USER_STREAMS)

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk) :: WAVELENGTH

!  Changed for Version 3.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR

!  Removed, Version 3.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: NSTREAMS_AZQUAD

!  New for Version 3.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Linearization flags

      LOGICAL :: DO_CHLORCONC_WF
      LOGICAL :: DO_WINDSPEED_WF

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]. Not needed now. See above. 9/9/15
!      REAL(fpk) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude 

!  Input Epoch

      INTEGER   :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk) :: FL_Amplitude755

!  Flag for using Data Gaussian parameters

      LOGICAL   :: FL_DO_DataGaussian

!  Input (non-data) Gaussians

      REAL(fpk) :: FL_InputGAUSSIANS(3,2)

!  Linearization of Gaussians only when the data option is not set

      LOGICAL   :: FL_F755_JACOBIANS
      LOGICAL   :: FL_GAUSS_JACOBIANS(6)

!  Exception handling
!  ------------------

!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  local variables
!  ===============

      CHARACTER (LEN=11), PARAMETER :: PREFIX = 'SLEAVESUP -'

      LOGICAL ::            ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER ::            I, FILUNIT, NM, M, QM

!  Initialize Exception handling

      STATUS = LRRS_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of LRRS Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Local error handling initialization

      ERROR  = .FALSE.
      NM     = NMESSAGES

!  Open file

      FILUNIT = LRRS_INUNIT
      OPEN(LRRS_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  Initialize Angle control
!  ========================

      DO_USER_STREAMS = .FALSE.
      NSTREAMS = 0

      BEAM_SZA = ZERO

      N_USER_STREAMS = 0
      DO I = 1, MAX_USER_STREAMS
        USER_ANGLES(I) = ZERO
      ENDDO
      N_USER_RELAZMS = 0
      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO

!  Initialize Surface stuff
!  ========================

!  Control flags

      DO_EXACT        = .FALSE.  !@@  New line
      DO_EXACTONLY    = .FALSE.
      DO_ISOTROPIC    = .FALSE.
      DO_SLEAVING     = .FALSE.
      DO_FLUORESCENCE = .FALSE.
      DO_WAV1         = .FALSE.  ! New, 9/9/15

      DO_SL_JACOBIANS  = .false.
      DO_ISO_JACOBIANS = .false.
      N_SLEAVE_WFS = 0

!  Fluorescence variables

      FL_LATITUDE   = ZERO
      FL_LONGITUDE  = ZERO
      FL_EPOCH      = 0
      FL_Amplitude755     = ZERO
      FL_DO_DataGaussian  = .false.
!mick fix 7/20/2016 - initialize
      FL_InputGAUSSIANS   = ZERO

      FL_F755_JACOBIANS   = .false.
      FL_GAUSS_JACOBIANS  = .false.

!  Water-leaving variables

      SALINITY   = ZERO
      CHLORCONC  = ZERO

      WINDSPEED  = ZERO
      WINDDIR    = ZERO

      DO_GlintShadow   = .FALSE.
      DO_FoamOption    = .FALSE.
      DO_FacetIsotropy = .FALSE.

      DO_CHLORCONC_WF  = .false.
      DO_WINDSPEED_WF  = .false.

!  Geometry and Input Control
!  ==========================

!  User-defined Stream angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR (ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                          ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAX_STREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
        'Number of half-space streams > MAX_STREAMS'
        ACTIONS(NM)  = &
         'Re-set input value or increase MAX_STREAMS dimension'
        STATUS = LRRS_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  Solar beam
!  ==========

!  BOA solar zenith angle input

      PAR_STR = 'Solar zenith angle (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        READ (FILUNIT,*,ERR=998) BEAM_SZA
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
         'Number of relative azimuth angles > MAX_USER_RELAZMS'
        ACTIONS(NM)  = &
         'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = LRRS_SERIOUS
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
          'Number of viewing zenith angles > MAX_USER_STREAMS'
          ACTIONS(NM)  = &
          'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = LRRS_SERIOUS
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
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Exact only flag. Only if above is set (!@@)

      IF ( DO_EXACT ) THEN
        PAR_STR = 'Do Exact-only (no Fourier-term contributions)?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_EXACTONLY
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )
      ENDIF

!  Basic source

      PAR_STR = 'Do surface-leaving Fluorescence?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_FLUORESCENCE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Single wavelength control

      PAR_STR = 'Do surface-leaving at 1 wavelength?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_WAV1
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Basic linearization flags. Should be true in Isotropic cases

      PAR_STR = 'Do surface-leaving Jacobians?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SL_JACOBIANS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

      IF ( DO_SL_JACOBIANS ) THEN
        PAR_STR = 'Do Isotropic surface-leaving Jacobians?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_ISO_JACOBIANS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )
      ENDIF

!  Local count for Jacobians

      QM = 0

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

!  Wavelength input superceded.....
!        PAR_STR = 'Wavelength in [Microns]'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,*,ERR=998) WAVELENGTH
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
!                             ACTIONS )

!  New for Version 3.7.   Wind-speed and directions, used for the first time

        PAR_STR = 'Windspeed in [m/s]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) WINDSPEED
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )

        PAR_STR = 'Wind direction (degrees) relative to sun position'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          READ (FILUNIT,*,ERR=998) WINDDIR
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )

!  New for Version 3.7. Control Flag for including Whitecaps

        PAR_STR = 'Do whitecap (foam) calculation?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) Do_FoamOption
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )

!  New for Version 3.7. Glint, Control Flags for Facet Isotropy and Shadowing
!     These are only needed for the non-Isotropic case

        if ( .not. do_Isotropic ) then
          PAR_STR = 'Do glint calculation with facet isotropy?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) Do_FacetIsotropy
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                               ACTIONS )

          PAR_STR = 'Do glint calculation with shadowing?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) Do_FacetIsotropy
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                               ACTIONS )
        endif

!  Removed for Version 3.7. --> Quadrature is now internal
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
!            STATUS = LRRS_SERIOUS
!            NMESSAGES = NM
!            GO TO 764
!          ENDIF
!        ENDIF

!  New for Version 3.7. Linearization flags

        IF ( DO_SL_JACOBIANS ) THEN
          PAR_STR = 'Do pigment concentration weighting function?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_CHLORCONC_WF
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          IF ( DO_CHLORCONC_WF ) QM = QM + 1

          PAR_STR = 'Do wind speed weighting function?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_WINDSPEED_WF
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          IF ( DO_WINDSPEED_WF ) QM = QM + 1
        ENDIF

!  Total number of Sleave weighting functions
!    Check not out of bounds

        N_SLEAVE_WFS = QM
        if ( N_SLEAVE_WFS .gt. MAX_SLEAVEWFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Sleave Water-Leaving Jacobians > MAX_SLEAVEWFS'
          ACTIONS(NM)  = 'Re-set input values or increase MAX_SLEAVEWFS dimension'
          STATUS = LRRS_SERIOUS
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
          STATUS = LRRS_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Use of Data Gaussians (New, 8 August 2012)
!    IF NOT SET, YOU MUST USE YOUR OWN PARAMETERS

        PAR_STR = 'Do Data Gaussians in Fluorescence?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_DO_DataGaussian
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Amplitude for FS755 (Nominally, this is one)

        PAR_STR = 'Amplitude for Fluorescence model at 755 nm'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_Amplitude755
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

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

!  Wavelength input superceded 9/9/15
!        PAR_STR = 'Wavelength for Fluorescence model in [nm]'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,*,ERR=998) FL_WAVELENGTH
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
!                           ACTIONS )

!  Linearization inputs for Fluorescence

        IF ( DO_SL_JACOBIANS .and. DO_ISO_JACOBIANS ) THEN
           PAR_STR = 'Do Jacobians for F755 Fluorescence value?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) FL_F755_JACOBIANS
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                              ACTIONS )
           IF ( FL_F755_JACOBIANS ) QM = QM + 1
        ENDIF

!  Gaussian linearizations
!     6 parameters: 1,2,3 for first Gaussian, 4,5,6 for second Gaussian

        IF ( DO_SL_JACOBIANS .and. DO_ISO_JACOBIANS .and..not. FL_DO_DataGaussian ) THEN
           PAR_STR = 'Do Gaussian parameter Jacobians for Fluorescence?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              DO M = 1, 6
                READ (FILUNIT,*,ERR=998) FL_GAUSS_JACOBIANS(m)
                IF ( FL_GAUSS_JACOBIANS(m) ) QM = QM + 1
              ENDDO
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        ENDIF

!  Total number of Sleave weighting functions
!   Check not out of bounds

        N_SLEAVE_WFS = QM
        if ( N_SLEAVE_WFS .gt. MAX_SLEAVEWFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Sleave Fluorescence Jacobians > MAX_SLEAVEWFS'
          ACTIONS(NM)  = 'Re-set input values or increase MAX_SLEAVEWFS dimension'
          STATUS = LRRS_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  End fluorescence parameters

      ENDIF

!  Successful finish

      CLOSE(FILUNIT)

!mick fix
      NMESSAGES = NM

!  Copy General Control inputs

      SLEAVE_Sup_In%SL_DO_USER_STREAMS  = DO_USER_STREAMS
      SLEAVE_Sup_In%SL_DO_SLEAVING      = DO_SLEAVING
      SLEAVE_Sup_In%SL_DO_FLUORESCENCE  = DO_FLUORESCENCE
      SLEAVE_Sup_In%SL_DO_ISOTROPIC     = DO_ISOTROPIC
      SLEAVE_Sup_In%SL_DO_EXACT         = DO_EXACT
      SLEAVE_Sup_In%SL_DO_EXACTONLY     = DO_EXACTONLY
      SLEAVE_Sup_In%SL_DO_WAV1          = DO_WAV1

      SLEAVE_LinSup_In%SL_DO_SL_JACOBIANS  = DO_SL_JACOBIANS
      SLEAVE_LinSup_In%SL_DO_ISO_JACOBIANS = DO_ISO_JACOBIANS

!  Copy Geometry results

      SLEAVE_Sup_In%SL_NSTREAMS          = NSTREAMS
      SLEAVE_Sup_In%SL_BEAM_SZA          = BEAM_SZA
      SLEAVE_Sup_In%SL_N_USER_RELAZMS    = N_USER_RELAZMS
      SLEAVE_Sup_In%SL_USER_RELAZMS      = USER_RELAZMS
      SLEAVE_Sup_In%SL_N_USER_STREAMS    = N_USER_STREAMS
      SLEAVE_Sup_In%SL_USER_ANGLES_INPUT = USER_ANGLES

!  Copy Water-leaving inputs
!  -------------------------

!  Original

      SLEAVE_Sup_In%SL_SALINITY         = SALINITY
      SLEAVE_Sup_In%SL_CHLORCONC        = CHLORCONC
      !SLEAVE_Sup_In%SL_WAVELENGTH       = WAVELENGTH

!  Version 3.7 changes

      SLEAVE_Sup_In%SL_WINDSPEED        = WINDSPEED
      SLEAVE_Sup_In%SL_WINDDIR          = WINDDIR

      SLEAVE_Sup_In%SL_DO_GlintShadow   = DO_GlintShadow
      SLEAVE_Sup_In%SL_DO_FoamOption    = DO_FoamOption
      SLEAVE_Sup_In%SL_DO_FacetIsotropy = DO_FacetIsotropy

      SLEAVE_LinSup_In%SL_DO_CHLORCONC_WF = DO_CHLORCONC_WF
      SLEAVE_LinSup_In%SL_DO_WINDSPEED_WF = DO_WINDSPEED_WF

!  Copy Fluorescence inputs
!  ------------------------

      SLEAVE_Sup_In%SL_FL_LATITUDE        = FL_LATITUDE
      SLEAVE_Sup_In%SL_FL_LONGITUDE       = FL_LONGITUDE
      !SLEAVE_Sup_In%SL_FL_WAVELENGTH      = FL_WAVELENGTH
      SLEAVE_Sup_In%SL_FL_EPOCH           = FL_EPOCH
      SLEAVE_Sup_In%SL_FL_Amplitude755    = FL_Amplitude755
      SLEAVE_Sup_In%SL_FL_DO_DataGaussian = FL_DO_DataGaussian
!mick fix 7/20/2016 - initialize input type structure variable
      SLEAVE_Sup_In%SL_FL_InputGAUSSIANS  = FL_InputGAUSSIANS

      SLEAVE_LinSup_In%SL_FL_F755_JACOBIANS  = FL_F755_JACOBIANS
      SLEAVE_LinSup_In%SL_FL_GAUSS_JACOBIANS = FL_GAUSS_JACOBIANS

!  Total number of weighting functions

      SLEAVE_LinSup_In%SL_N_SLEAVE_WFS = N_SLEAVE_WFS

!  Exception handling

      SLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD = STATUS
      SLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES   = NMESSAGES
      SLEAVE_Sup_InputStatus%SL_INPUTMESSAGES    = MESSAGES
      SLEAVE_Sup_InputStatus%SL_INPUTACTIONS     = ACTIONS

!  Normal return

      RETURN

!  Open file error

300   CONTINUE
      STATUS = LRRS_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
      GO TO 764

!  line read error - abort immediately

998   CONTINUE
      STATUS = LRRS_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'read failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Re-set: Entry is incorrect in input file'
      CLOSE(FILUNIT)
      GO TO 764

!  Final error copying

764   CONTINUE

      SLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD = STATUS
      SLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES   = NMESSAGES
      SLEAVE_Sup_InputStatus%SL_INPUTMESSAGES    = MESSAGES
      SLEAVE_Sup_InputStatus%SL_INPUTACTIONS     = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE SLEAVE_LIN_INPUTMASTER

!

      SUBROUTINE SLEAVE_LIN_MAINMASTER ( &
        N_Lambdas, Lambdas,       & ! Input wavelengths [nm] for LRRS
        SLEAVE_Sup_In,            & ! Inputs
        SLEAVE_LinSup_In,         & ! Inputs
        SLEAVE_Sup_Out,           & ! Outputs
        SLEAVE_LinSup_Out )         ! Outputs

!  Prepares the Surface Leaving (and linearizations) necessary for LRRS.
 
      USE LRRS_pars_m

      USE sleave_sup_inputs_def_m
      USE sleave_linsup_inputs_def_m

      USE sleave_sup_outputs_def_m
      USE sleave_linsup_outputs_def_m

      USE sleave_sup_aux_m          , only : GETQUAD2
      USE sleave_sup_routines_m     , only : WaterLeaving,              &
                                             get_fluorescence_755,      &
                                             solar_spec_irradiance
      USE sleave_lin_sup_routines_m , only : Linearized_WaterLeaving

      IMPLICIT NONE

!  Number of wavelengths used. Introduced 9/8/15 for LRRS
!     = 1, if single-value option is set (dimensioning also 1 - check this!)
!     = NPOINTS_INNER for multiple points (Bin realization)
!     = NPOINTS_MONO  for multiple points (Mono realization)

      INTEGER, intent(in)   :: N_Lambdas

!  Input wavelengths in [nm]. Array introduced 9/8/15 for LRRS.
!    For single-point (Bin realization) , set to average-value of "Lamdas_Ranked" 
!    For single-point (Mono realization), set to Excitation wavelength "Lamdas_Ranked(W_EXCIT)" 
!    For all points, copy LAMDAS_RANKED array (both realizations)
!    - Note units are same as LAMBDAS_RANKED which is in [nm].

      REAL(fpk), intent(in) :: Lambdas ( MAX_POINTS )

!  Input structures
!  ----------------

      TYPE(SLEAVE_Sup_Inputs)   , INTENT(IN)   :: SLEAVE_Sup_In
      TYPE(SLEAVE_LinSup_Inputs), INTENT(IN)   :: SLEAVE_LinSup_In

!  Output structures
!  -----------------

      TYPE(SLEAVE_Sup_Outputs)   , INTENT(OUT) :: SLEAVE_Sup_Out
      TYPE(SLEAVE_LinSup_Outputs), INTENT(OUT) :: SLEAVE_LinSup_Out

!  LRRS local variables
!  ++++++++++++++++++++++

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

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: DO_EXACT
      LOGICAL :: DO_EXACTONLY

!  Do single wavelength Flag. New, 9/9/15

      LOGICAL :: DO_WAV1

!  Linearizations

      LOGICAL :: DO_SL_JACOBIANS
      LOGICAL :: DO_ISO_JACOBIANS

!  Geometry and control
!  --------------------

!  Stream angle flag

      LOGICAL :: DO_USER_STREAMS

!  Number of discrete ordinate streams, quadrature
!   Version 3.7, added the qudrature arrays

      INTEGER   :: NSTREAMS
      REAL(fpk) :: STREAMS (MAX_STREAMS)
      REAL(fpk) :: WEIGHTS (MAX_STREAMS)

!  Local angle control

      INTEGER   :: N_USER_STREAMS
      INTEGER   :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: BEAM_SZA
      REAL(fpk) :: USER_RELAZMS(MAX_USER_RELAZMS)
      REAL(fpk) :: USER_ANGLES (MAX_USER_STREAMS)

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]. Superceded, 9/9/15
!      REAL(fpk) :: WAVELENGTH

!  Changed for Version 3.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR

!  Removed, Version 3.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: NSTREAMS_AZQUAD

!  New for Version 3.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL :: DO_GlintShadow
      LOGICAL :: DO_FoamOption
      LOGICAL :: DO_FacetIsotropy

!  Linearization flags

      LOGICAL :: DO_CHLORCONC_WF
      LOGICAL :: DO_WINDSPEED_WF

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]. Superceded, 9/9/15
!      Double precision :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude

!  Input Epoch

      INTEGER   :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk) :: FL_Amplitude755

!  Flag for using Data Gaussian parameters

      LOGICAL   :: FL_DO_DataGaussian

!  Linearization flags
!     Linearization of Gaussians only when the data option is not set

      LOGICAL   :: FL_F755_JACOBIANS, FL_GAUSS_JACOBIANS(6)

!  Total number of SLEAVE weighting functions
!  ------------------------------------------

      INTEGER   :: N_SLEAVE_WFS

!  Local functions
!  ===============

!  Fluorescence Isotropic Surface leaving term

   REAL(fpk) :: Fluor_ISOTROPIC

!  Isotropic value. Fast calculation

   REAL(fpk)    :: WLeaving_ISO
   REAL(fpk)    :: dC_WLeaving_ISO
   REAL(fpk)    :: dW_WLeaving_ISO

!  Input solar, output stream angles

   REAL(fpk)    :: WLeaving_SD    ( MAX_STREAMS )
   REAL(fpk)    :: dC_WLeaving_SD ( MAX_STREAMS )
   REAL(fpk)    :: dW_WLeaving_SD ( MAX_STREAMS )

!  input solar, output view angles

   REAL(fpk)    :: WLeaving_SV    ( MAX_USER_STREAMS )
   REAL(fpk)    :: dC_WLeaving_SV ( MAX_USER_STREAMS )
   REAL(fpk)    :: dW_WLeaving_SV ( MAX_USER_STREAMS )

!  Exact Surface-Leaving term
!      REAL(fpk), dimension ( MAXSTOKES, MAX_USER_STREAMS, &
!        MAX_USER_RELAZMS, MAXBEAMS ) :: SLTERM_USERANGLES
!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams
!      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAX_STREAMS )       :: SLTERM_F_0
!      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS ) :: USER_SLTERM_F_0

!  Other local variables
!  =====================

!  (Version 3.6 code). Water-leaving model
!      REAL :: WAV,CHL,RW,SAL,A,REFR,REFI,N12,RWB,TDS,TDV

!  Fluorescence Gaussian parameters (Data)
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

!  Solar spectral radiance model

      REAL(FPK) :: ssr_wvl

!  Fluorescence model

      CHARACTER*60 :: Fluofile
      INTEGER   :: K, NSLPARS, M, Q, UM
      REAL(FPK) :: Fs755, FL_SunSpec, FsSum, d_Fssum(6)
      REAL(FPK) :: ampli, lamda, sigma, arg, gauss, var, ff, expl, help1, help2
      REAL(fpk) :: LOCAL_WAVS( MAX_POINTS ), FL_wavelength, wavelength

!  Misc

      INTEGER   :: L,LOCAL_NPOINTS

!  Copy from input structure
!  -------------------------

!  Copy Top-level general Control inputs

      DO_USER_STREAMS  = SLEAVE_Sup_In%SL_DO_USER_STREAMS
      DO_SLEAVING      = SLEAVE_Sup_In%SL_DO_SLEAVING

      DO_EXACT         = SLEAVE_Sup_In%SL_DO_EXACT          !@@
      DO_EXACTONLY     = SLEAVE_Sup_In%SL_DO_EXACTONLY
      DO_FLUORESCENCE  = SLEAVE_Sup_In%SL_DO_FLUORESCENCE
      DO_ISOTROPIC     = SLEAVE_Sup_In%SL_DO_ISOTROPIC

!  Do single wavelength Flag. New, 9/9/15

      DO_WAV1          = SLEAVE_Sup_In%SL_DO_WAV1

      DO_SL_JACOBIANS  = SLEAVE_LinSup_In%SL_DO_SL_JACOBIANS
      DO_ISO_JACOBIANS = SLEAVE_LinSup_In%SL_DO_ISO_JACOBIANS

!  Set number of streams
!   Stream Quadrature new for Version 3.7

      NSTREAMS = SLEAVE_Sup_In%SL_NSTREAMS
      CALL GETQUAD2 ( ZERO, ONE, NSTREAMS, STREAMS, WEIGHTS  )      

!   !@@ Observational Geometry + Solar sources Optionalities
!   !@@ Either set from User Observational Geometry
!          Or Copy from Usual lattice input

      BEAM_SZA       = SLEAVE_Sup_In%SL_BEAM_SZA
      N_USER_RELAZMS = SLEAVE_Sup_In%SL_N_USER_RELAZMS
      USER_RELAZMS   = SLEAVE_Sup_In%SL_USER_RELAZMS
      N_USER_STREAMS = SLEAVE_Sup_In%SL_N_USER_STREAMS
      USER_ANGLES    = SLEAVE_Sup_In%SL_USER_ANGLES_INPUT

!  Copy Water-leaving inputs
!  -------------------------

!  Original

      SALINITY        = SLEAVE_Sup_In%SL_SALINITY
      CHLORCONC       = SLEAVE_Sup_In%SL_CHLORCONC
!      WAVELENGTH      = SLEAVE_Sup_In%SL_WAVELENGTH

      WINDSPEED       = SLEAVE_Sup_In%SL_WINDSPEED
      WINDDIR         = SLEAVE_Sup_In%SL_WINDDIR

      DO_GlintShadow   = SLEAVE_Sup_In%SL_DO_GlintShadow
      DO_FoamOption    = SLEAVE_Sup_In%SL_DO_FoamOption
      DO_FacetIsotropy = SLEAVE_Sup_In%SL_DO_FacetIsotropy

      DO_ChlorConc_WF = SLEAVE_LinSup_In%SL_DO_ChlorConc_WF
      DO_WindSpeed_WF = SLEAVE_LinSup_In%SL_DO_WindSpeed_WF

!  Copy Fluorescence inputs
!  ------------------------

!  Main variables

!      FL_Wavelength   = SLEAVE_Sup_In%SL_FL_Wavelength
      FL_Latitude     = SLEAVE_Sup_In%SL_FL_Latitude
      FL_Longitude    = SLEAVE_Sup_In%SL_FL_Longitude
      FL_Epoch        = SLEAVE_Sup_In%SL_FL_Epoch
      FL_Amplitude755 = SLEAVE_Sup_In%SL_FL_Amplitude755
      FL_DO_DataGaussian = SLEAVE_Sup_In%SL_FL_DO_DataGaussian

!mick fix 8/31/2012 - added outer if block
      if (DO_FLUORESCENCE) then
        if ( FL_DO_DataGaussian ) then
           FL_GAUSSIANS(1:3,1) = FL_DataGAUSSIANS(1:3,1) 
           FL_GAUSSIANS(1:3,2) = FL_DataGAUSSIANS(1:3,2)
        else
           FL_GAUSSIANS(1:3,1) = SLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,1)
           FL_GAUSSIANS(1:3,2) = SLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,2)
        endif
      endif

      FL_F755_JACOBIANS   = SLEAVE_LinSup_In%SL_FL_F755_JACOBIANS
      FL_GAUSS_JACOBIANS  = SLEAVE_LinSup_In%SL_FL_GAUSS_JACOBIANS

!  Total number of SLEAVE weighting functions
!  ------------------------------------------

      N_SLEAVE_WFS = SLEAVE_LinSup_In%SL_N_SLEAVE_WFS

!  Main code
!  ---------

!  Rob Fix. 9/9/15. Initialize the wavelength array
!    Generally only one point, unless control has been set for NewCM
!      - the 0.001 factor is conversion to Microns (Waterleaving only)

      Local_npoints = 1 ; local_wavs = zero
      if ( DO_FLUORESCENCE ) then
         if ( DO_Wav1 ) then
            local_wavs(1) = SUM(lambdas(1:N_lambdas))/DBLE(N_Lambdas)
         else
            Local_npoints = N_Lambdas
            local_wavs(1:N_lambdas) = lambdas(1:N_lambdas)
         endif
      else
         if ( DO_Wav1 ) then
            local_wavs(1) = 0.001_fpk * SUM(lambdas(1:N_lambdas))/DBLE(N_Lambdas)
         else
            Local_npoints = N_Lambdas
            local_wavs(1:N_lambdas) = 0.001_fpk * lambdas(1:N_lambdas)
         endif
      endif

!  Zero the output

      SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC  = ZERO
      SLEAVE_Sup_Out%SL_SLTERM_USERANGLES = ZERO
      SLEAVE_Sup_Out%SL_SLTERM_F_0        = ZERO
      SLEAVE_Sup_Out%SL_USER_SLTERM_F_0   = ZERO

!  Zero the linearized output (8/8/12, only isotropic defined so far)

      SLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC  = ZERO
      SLEAVE_LinSup_Out%SL_LS_SLTERM_USERANGLES = ZERO
      SLEAVE_LinSup_Out%SL_LS_SLTERM_F_0        = ZERO
      SLEAVE_LinSup_Out%SL_LS_USER_SLTERM_F_0   = ZERO

!  Fluorescence
!  ============

      IF ( DO_FLUORESCENCE ) THEN

!  Temporary - Only Isotropic yet.

        IF ( .not.DO_ISOTROPIC ) &
          Stop 'Non-isotropic not allowed yet if doing fluorescence'

!  F_755 data file

        Fluofile = 'lrrs_test/fluorescence_data/fluor_data_2009_fortran.dat'

!  For Given SZA, Get the F_755 data from the subroutine

        CALL get_fluorescence_755 &
          ( FL_Latitude, FL_Longitude, FL_Epoch, BEAM_SZA, FluoFile, Fs755 )

!  Start the LRRS wavelength loop

        DO L = 1, Local_Npoints

!  Get solar spectral irradiance, in (W m−2 μm−1), to normalize data

          !FL_SunSpec = 1.0d0  ! Temporary

          FL_Wavelength = Local_wavs(L)

          ssr_wvl = FL_Wavelength*1.0d-3 !convert from nm to um
          FL_SunSpec = solar_spec_irradiance( ssr_wvl )

!  factor: After some fiddling, this is 1.0 (July 30th, 2012)
!    4pi factor is required in DB correction term

!          FF = PI4
          FF = 1.0d0
!          FF = 0.0d0

!  Zero the running total of weighting functions

          NSLPARS = 0

!  Compute double Gaussian sums

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

!  Assign output Fluorescence (Apply Amplitude)
!  multiply by Fs755, and normalize to solar spectrum
!   FF is the fudge factor

          help2 = FF * FsSum / FL_SunSpec
          Fluor_ISOTROPIC = help2 * Fs755
          SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(L) = FL_Amplitude755 * Fluor_ISOTROPIC

!  Apply Main linearization (w.r.t Fs755 Amplitude)

          if ( DO_ISO_JACOBIANS .and. FL_F755_JACOBIANS ) then
!mick fix 7/20/2016 - changed index of 2nd dim from IB to L
             !NSLPARS = 1 ; SLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(1,IB) = Fluor_ISOTROPIC
             NSLPARS = 1
             SLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(NSLPARS,L) = Fluor_ISOTROPIC
          endif

!  Applyl Gaussian linearizations if flagged.

          if ( DO_ISO_JACOBIANS .and. .not. FL_DO_DataGaussian ) then
             help1 = FL_Amplitude755 * FF * Fs755 / FL_SunSpec;  d_Fssum = zero
             do k = 1, 2
                ampli = FL_Gaussians(1,k)
                lamda = FL_Gaussians(2,k)
                sigma = FL_Gaussians(3,k)
                var = 0.5d0/sigma/sigma
                arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
                m = 3*(k-1);  Gauss = zero
                if ( arg.lt.88.0d0 ) then
                   expl  = dexp ( - arg )
                   gauss = ampli * expl
                   d_Fssum(m+1)  = expl
                   d_Fssum(m+2)  = + gauss * two * ( FL_Wavelength - lamda ) * var
                   d_Fssum(m+3)  = + gauss * two * arg / sigma
                endif
             enddo
             do m = 1, 6
                IF ( FL_GAUSS_JACOBIANS(m) ) then
                   NSLPARS = NSLPARS + 1
                   SLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(NSLPARS,L) = help1 * d_Fssum(m)
                ENDIF
             enddo
          endif

!mick add 7/20/2016 - added consistency check (not active at this time)
!  Consistency check for total number of SLEAVE weighting functions
          !IF ( NSLPARS .ne. N_SLEAVE_WFS ) then
          !   NMESSAGES = NMESSAGES + 1
          !   MESSAGES(NMESSAGES) = 'Fatal - Bookkeeping incorrect for number of SLEAVE WFs'
          !   STATUS = LIDORT_SERIOUS
          !   GO TO 899
          !ENDIF

!         write(*,*)FL_Wavelength,FsSum, FL_SunSpec, SLTERM_ISOTROPIC

!  End Points loop

        ENDDO

      ENDIF

!  WATER-LEAVING
!  =============

      IF ( .not. DO_FLUORESCENCE ) THEN

!  Start the wavelength loop

         do L = 1, Local_npoints

!  wavelength (already in Microns)

            Wavelength = Local_wavs(L)

!  Call the routine

            if ( .not. DO_CHLORCONC_WF .and. .not. DO_WINDSPEED_WF ) then
               Call WaterLeaving &
                ( MAX_USER_STREAMS, MAX_STREAMS,                                     &
                  do_Isotropic, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,     &
                  Wavelength, Salinity, ChlorConc, Windspeed, WindDir,               &
                  n_user_streams, nstreams, beam_sza, user_angles, streams,          &
                  WLeaving_ISO, WLeaving_SD, WLeaving_SV )
            else
               Call Linearized_WaterLeaving &
                ( MAX_USER_STREAMS, MAX_STREAMS,                                     &
                  do_Isotropic, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,     &
                  Wavelength, Salinity, ChlorConc, Windspeed, WindDir,               &
                  n_user_streams, nstreams, beam_sza, user_angles, streams,          &
                  WLeaving_ISO, WLeaving_SD, WLeaving_SV,                            &
                  dC_WLeaving_ISO, dC_WLeaving_SD, dC_WLeaving_SV,                   &
                  dW_WLeaving_ISO, dW_WLeaving_SD, dW_WLeaving_SV )
            endif

!  Copy to Type structure arrays

!    If Isotropic, LRRS takes care of the rest

            SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(L) = WLeaving_ISO
            if ( .not. do_Isotropic ) then
               do um = 1, n_user_streams
                  SLEAVE_Sup_Out%SL_SLTERM_USERANGLES(um,1:n_user_relazms,L) = WLeaving_SV(um)
                  SLEAVE_Sup_Out%SL_USER_SLTERM_F_0  (0,um,L)                = WLeaving_SV(um)
               enddo
               SLEAVE_Sup_Out%SL_SLTERM_F_0(0,1:nstreams,L) = WLeaving_SD(1:nstreams)
            endif

!  Copy to linearized Type structure arrays
!  ----------------------------------------

            Q = 0 ! start Jacobian count

!   Pigment

            if ( DO_CHLORCONC_WF ) then
               Q = Q + 1 ; SLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(Q,L) = dC_WLeaving_ISO
               if ( .not. do_Isotropic ) then
                  do um = 1, n_user_streams
                     SLEAVE_LinSup_Out%SL_LS_SLTERM_USERANGLES(Q,um,1:n_user_relazms,L) = dC_WLeaving_SV(um)
                     SLEAVE_LinSup_Out%SL_LS_USER_SLTERM_F_0  (Q,0,um,L)                = dC_WLeaving_SV(um)
                  enddo
                  SLEAVE_LinSup_Out%SL_LS_SLTERM_F_0(Q,0,1:nstreams,L) = dC_WLeaving_SD(1:nstreams)
               endif
            endif

!   Windspeed

            if ( DO_WINDSPEED_WF ) then
               Q = Q + 1 ; SLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(Q,L) = dW_WLeaving_ISO
               if ( .not. do_Isotropic ) then
                  do um = 1, n_user_streams
                     SLEAVE_LinSup_Out%SL_LS_SLTERM_USERANGLES(Q,um,1:n_user_relazms,L) = dW_WLeaving_SV(um)
                     SLEAVE_LinSup_Out%SL_LS_USER_SLTERM_F_0  (Q,0,um,L)                = dW_WLeaving_SV(um)
                  enddo
                  SLEAVE_LinSup_Out%SL_LS_SLTERM_F_0(Q,0,1:nstreams,L) = dW_WLeaving_SD(1:nstreams)
               endif
            endif

!mick add 7/20/2016 - added consistency check (not active at this time)
!  Consistency check for total number of SLEAVE weighting functions
          !IF ( Q .ne. N_SLEAVE_WFS ) then
          !   NMESSAGES = NMESSAGES + 1
          !   MESSAGES(NMESSAGES) = 'Fatal - Bookkeeping incorrect for number of SLEAVE WFs'
          !   STATUS = LIDORT_SERIOUS
          !   GO TO 899
          !ENDIF

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
!        SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1,L)= DBLE(Rwb)

!  End wavelength loop

         enddo

!  end water leaving

      endif

!  Finish

      RETURN
      END SUBROUTINE SLEAVE_LIN_MAINMASTER

      END MODULE sleave_linsup_masters_m

