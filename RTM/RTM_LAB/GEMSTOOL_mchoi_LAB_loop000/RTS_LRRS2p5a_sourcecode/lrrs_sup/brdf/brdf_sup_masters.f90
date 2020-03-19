
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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            BRDF_INPUTMASTER                                 #
! #            BRDF_MAINMASTER                                  #
! #                                                             #
! ###############################################################

!  Developed for LRRS Version 2.5, 9/8/15. R. Spurr, RT SOLUTIONS Inc.
!   Closely follows the LIDORT module with the same name, but
!     (1) No surface emission, solar angle dimensioning, No observational geometry.
!     (2) Additional wavelength dimensioning (MAX_BRDF_POINTS)
!     (3) Additional control for using 1 or all wavelengths

!  Note (9/8/15). The number of wavelengths for BRDF has deliberately
!  been left flexible - typically the BRDF properties will change very
!  little over a Raman-scattering window (+/- 2 nm in the UV), so to
!  a very good approximation, it is sufficient to use one point for
!  the calculations of BRDF - in which case, MAX_BRDF_POINTS will be 1

!  It turns out that all but one of the BRDF Kernels are wavelength
!  independent, so this consideration only applies to the "NewCM" GLINT
!  BRDF calculation, which has a wavelength input for foam reflectance
!  and refractive index through the ocean water-leaving part.

!  Now, whenever the BRDF supplement is used with LRRS, the choice of 
!  BRDF wavelengths is linked to the Raman wavelengths (LAMDAS_RANKED).
!  These wavelengths are now input to the BRDF supplement MASTER, and
!  they are not set by hand or by configuration-file read. 

!  When NewCM Glint is turned on, and you want just one wavelength, then
!  the control flag DO_NewCM_Wav1 will be set True.  In this case, the
!  number of BRDF "NewCM" wavelengths NewCM_Npoints = 1, and the single
!  wavelength NewCM_Lambdas(1) = AVERAGE VALUE of LAMBDAS_RANKED.

!  When NewCM Glint is turned on, and you want all wavelengths, then
!  control flag DO_NewCM_Wav1 is False, and NewCM_Npoints = N_Lambdas_Ranked
!  and NewCM_Lambdas(1:NewCM_npoints) = LAMBDAS_RANKED(1:NewCM_npoints).
!  A Glint calculation will be done for all Raman-scattered wavelengths!!!

      MODULE brdf_sup_masters_m

! #####################################################################
! #####################################################################

      PRIVATE
      PUBLIC :: BRDF_INPUTMASTER, &
                BRDF_MAINMASTER

      CONTAINS

      SUBROUTINE BRDF_INPUTMASTER ( &
        FILNAM,              & ! Input
        BRDF_Sup_In,         & ! Outputs
        BRDF_Sup_InputStatus ) ! Outputs

!  Input routine for BRDF program

      USE LRRS_pars_m
      USE BRDF_FINDPAR_m

      USE brdf_sup_inputs_def_m
      USE brdf_sup_outputs_def_m

!  Implicit none

      IMPLICIT NONE

!  Module arguments (input filename)

      CHARACTER(LEN=*), intent(in)  :: FILNAM

!  Structures

      TYPE(BRDF_Sup_Inputs), intent(out) :: BRDF_Sup_In

!  Output variables

      TYPE(BRDF_Input_Exception_Handling), intent(out) :: BRDF_Sup_InputStatus

!  LRRS local variables
!  ++++++++++++++++++++++

!  Stream angle flag

      LOGICAL   :: DO_USER_STREAMS

!  BRDF surface flag
!    ---> Really should be true here

      LOGICAL   :: DO_BRDF_SURFACE

!  Number and index-list and names of bidirectional functions

      INTEGER   ::          N_BRDF_KERNELS
      INTEGER   ::          WHICH_BRDF ( MAX_BRDF_KERNELS )
      CHARACTER (LEN=10) :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER   :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(fpk) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Lambertian Surface control

      LOGICAL   :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      REAL(fpk) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 02 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised 12 August 2014; added output option flag

      LOGICAL   :: DO_WSABSA_OUTPUT
      LOGICAL   :: DO_WSA_SCALING
      LOGICAL   :: DO_BSA_SCALING
      REAL(fpk) :: WSA_VALUE, BSA_VALUE

!  Number of azimuth quadrature streams for BRDF

      INTEGER   :: NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL   :: DO_SHADOW_EFFECT

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL   :: DO_EXACT
!      LOGICAL   :: DO_EXACTONLY
      LOGICAL   :: DO_DBONLY

!  Multiple reflectance correction for Glitter kernels

      LOGICAL   :: DO_MSRCORR
      INTEGER   :: MSRCORR_ORDER
      LOGICAL   :: DO_MSRCORR_DBONLY       ! Rob Fix 9/25/14, Variable renamed
      INTEGER   :: MSRCORR_NMUQUAD
      INTEGER   :: MSRCORR_NPHIQUAD

!  Number of discrete ordinate streams

      INTEGER   :: NSTREAMS

!  Solar angle

      REAL(fpk) :: BEAM_SZA

!  Local angle control

      INTEGER   :: N_USER_STREAMS
      INTEGER   :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: USER_RELAZMS      (MAX_USER_RELAZMS)
      REAL(fpk) :: USER_ANGLES       (MAX_USER_STREAMS)

!  NewCM Glitter options (bypasses the usual Kernel system)
!  -------------------------------------------------------

!  Overall flag for this option

      LOGICAL   :: DO_NewCMGLINT

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input wavelength single choice control

      LOGICAL   :: DO_NewCM_Wav1

!  Number of wavelengths used. Introduced 9/8/15 for LRRS
!     = 1, if single-value option is set (dimensioning also 1 - check this!)
!     = NPOINTS_INNER for multiple points (Bin realization)
!     = NPOINTS_MONO  for multiple points (Mono realization)

      INTEGER   :: NewCM_Npoints

!  Input wavelengths in [nm]. Array introduced 9/8/15 for LRRS.
!    For single-point (Bin realization) , set to average-value of "Lamdas_Ranked" 
!    For single-point (Mono realization), set to Excitation wavelength "Lamdas_Ranked(W_EXCIT)" 
!    For all points, copy LAMDAS_RANKED array (both realizations)
!    - Note units are same as LAMBDAS_RANKED which is in [nm].

      REAL(fpk) :: NewCM_Lambdas ( MAX_BRDF_POINTS )

!  Input wavelength in [Microns]. Now taken from LRRS inputs 9/8/15
!      REAL(fpk) :: WAVELENGTH

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR

!  Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Exception handling
!  ------------------

!     Message Length should be at least 120 Characters

      INTEGER        :: STATUS
      INTEGER        :: NMESSAGES
      CHARACTER*120  :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*120  :: ACTIONS (0:MAX_MESSAGES)

!  local variables
!  ===============

      CHARACTER(Len=9), parameter :: PREFIX = 'BRDFSUP -'

      LOGICAL           :: ERROR
      CHARACTER(Len=80) :: PAR_STR
      INTEGER           :: I, K, L, FILUNIT, NM

!  MODIS-style brdfs, some checks on the input

      LOGICAL            :: DO_ROSS, DO_LI, DO_LAMB, DO_CM

!  Check list of Kernel names
!  Rob Extension 12/2/14. BPDF Kernels (replace BREONVEG, BREONSOIL)

      CHARACTER (LEN=10) :: BRDF_CHECK_NAMES ( MAXBRDF_IDX )
      BRDF_CHECK_NAMES = (/ &
                           'Lambertian', &
                           'Ross-thin ', &
                           'Ross-thick', &
                           'Li-sparse ', &
                           'Li-dense  ', &
                           'Hapke     ', &
                           'Roujean   ', &
                           'Rahman    ', &
                           'Cox-Munk  ', &
                           'BPDF-Vegn ', &
                           'BPDF-Soil ', &
                           'BPDF-NDVI ', &
                           'NewCMGlint'/)

!  Initialize Exception handling

      STATUS = LRRS_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of LRRS BRDF Input file'
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

      BEAM_SZA       = zero
      N_USER_STREAMS = 0
      USER_ANGLES    = ZERO
      N_USER_RELAZMS = 0
      USER_RELAZMS   = ZERO

!  Initialize Surface stuff
!  ========================

      NSTREAMS_BRDF  = 0
      N_BRDF_KERNELS = 0
!mick fix 10/29/2015 - initialize N_BRDF_PARAMETERS
      N_BRDF_PARAMETERS = 0

      DO_SHADOW_EFFECT    = .FALSE.

!  Rob Fix 9/25/14. Two variables replaced
!      DO_EXACT            = .FALSE.       !@@  New line
!      DO_EXACTONLY        = .FALSE.
      DO_DBONLY      = .false.

      DO_MSRCORR           = .FALSE.
      MSRCORR_ORDER        = 0
      DO_MSRCORR_DBONLY    = .FALSE.
      MSRCORR_NMUQUAD      = 0
      MSRCORR_NPHIQUAD     = 0

      WHICH_BRDF = 0
      BRDF_NAMES = '          '
      DO K = 1, MAX_BRDF_KERNELS
        LAMBERTIAN_KERNEL_FLAG(K) = .FALSE.
        BRDF_FACTORS(K) = ZERO
        DO L = 1, MAX_BRDF_PARAMETERS
          BRDF_PARAMETERS(K,L) = ZERO
        ENDDO
      ENDDO

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 02-07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised 12 August 2014; added output option flag

      DO_WSABSA_OUTPUT = .false.
      DO_WSA_SCALING = .false.
      DO_BSA_SCALING = .false.
      WSA_VALUE      = zero
      BSA_VALUE      = zero

!  NewCM options

      DO_NewCMGLINT = .false.
      DO_NewCM_Wav1 = .false.
!mick fix 10/29/2015 - initialize BS_NewCM_Npoints & BS_NewCM_Lambdas
      NewCM_Npoints = 0
      NewCM_Lambdas = zero

      SALINITY   = zero
      WINDSPEED  = zero
      WINDDIR    = zero
      DO_GlintShadow   = .false.
      DO_FoamOption    = .false.
      DO_FacetIsotropy = .false.

!  Read Angle stuff
!  ================

!  User-defined viewing zenith angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAX_STREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_STREAMS dimension'
        STATUS = LRRS_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  Solar beam
!  ==========

!  TOA solar zenith angle input

      PAR_STR = 'Solar zenith angle (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          READ (FILUNIT,*,ERR=998) BEAM_SZA
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Azimuth angles
!  ==============

!  Number of angles

      PAR_STR = 'Number of user-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) =  'Number of relative azimuth angles > maximum dimension'
        ACTIONS(NM)  =  'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = LRRS_SERIOUS
        NMESSAGES    = NM
        GO TO 764
      ENDIF

!  Angles

      PAR_STR = 'User-defined relative azimuth angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_RELAZMS
          READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  User defined viewing zenith angles (should be positive)
!  ==================================

      IF ( DO_USER_STREAMS ) THEN

!  Number of angles

        PAR_STR = 'Number of user-defined viewing zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check dimension

        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of viewing zenith angles > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = LRRS_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Angles

        PAR_STR = 'User-defined viewing zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      ENDIF

!  Surface stuff
!  =============

!  BRDF input
!  ----------

!  Basic flag

      PAR_STR = 'Do BRDF surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_BRDF_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  NewCM Flag. New for Version 3.7

      IF ( DO_BRDF_SURFACE ) THEN
        PAR_STR = 'Do NewCM Ocean BRDF reflectance?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_NewCMGLINT
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Only if set

      IF ( DO_BRDF_SURFACE ) THEN

!  Get the NewCM parameters
!  -----------------------

        IF ( DO_NewCMGLINT ) then

!  Flags

          PAR_STR = 'Do NewCM glint shadowing?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_GlintShadow
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'Do NewCM whitecap (foam) reflectance?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_FoamOption
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'Do NewCM facet isotropy?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_FacetIsotropy
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Wavelength, Salinity

          PAR_STR = 'Do NewCM single wavelength option?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_NewCM_Wav1
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'NewCM Ocean water salinity [ppt]'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) SALINITY
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Windspeed and direction

          PAR_STR = 'NewCM Windspeed in [m/s]'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) WINDSPEED
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'NewCM Wind direction (degrees) relative to sun position'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              READ (FILUNIT,*,ERR=998) WINDDIR
          ENDIF
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        ENDIF

!  NewCM input: Set kernel defaults
!     One kernel, scalar only, no surface emission, no scaling

        IF ( DO_NewCMGLINT ) then
           N_BRDF_KERNELS       = 1
           BRDF_NAMES(1)        = 'NewCMGlint'
           WHICH_BRDF(1)        = NewCMGLINT_IDX
           BRDF_FACTORS(1)      = one
           N_BRDF_PARAMETERS(1) = 2
           BRDF_PARAMETERS(1,1) = WINDSPEED
           BRDF_PARAMETERS(1,2) = SALINITY
        endif

!  Skip next section if doing the NewCM kernel

        IF ( DO_NewCMGLINT ) go to 656

!  Basic KERNEL BRDF inputs
!  ------------------------

!  Number of kernels

        PAR_STR = 'Number of BRDF kernels'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_BRDF_KERNELS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check Dimension

        IF ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of BRDF Kernels > maximum dimension (=3)'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_BRDF_KERNELS dimension'
          STATUS = LRRS_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Main kernel input
!   Rob Fix 9/25/14. F8.4 format for the BRDF_FACTORS

        PAR_STR = 'Kernel names, indices, amplitudes, # parameters, parameters'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,56,ERR=998) &
               BRDF_NAMES(I), WHICH_BRDF(I), BRDF_FACTORS(I), &
              N_BRDF_PARAMETERS(I),(BRDF_PARAMETERS(I,K),K=1,3)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
! 56     FORMAT( A10, I2, F6.2, I2, 3F12.6 )
 56     FORMAT( A10, I2, F8.4, I2, 3F12.6 )

!  Check Kernel indices are within bounds. Check BRDF name is on accepted list

        DO K = 1, N_BRDF_KERNELS
          IF ( WHICH_BRDF(K).GT.MAXBRDF_IDX.OR.WHICH_BRDF(K).LE.0) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: BRDF Index not on list of indices'
            ACTIONS(NM)  = 'Re-set input value: Look in LRRS_PARS for correct index'
            STATUS = LRRS_SERIOUS
            NMESSAGES = NM
            GO TO 764
          ELSE
            IF ( BRDF_NAMES(K).NE.BRDF_CHECK_NAMES(WHICH_BRDF(K)) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: BRDF kernel name not one of accepted list'
              ACTIONS(NM)  = 'Re-set input value: Look in LRRS_PARS for correct name'
              STATUS = LRRS_SERIOUS
              NMESSAGES = NM
              GO TO 764
            ENDIF
          ENDIF
        ENDDO

!  **************************************************************
!  Rob Fix 9/25/14. New section on use of MODIS-type kernels

        DO_ROSS = .false. ; do_Li = .false. ; do_Lamb = .false. ; do_CM = .false.
        DO K = 1, N_BRDF_KERNELS
          IF ( WHICH_BRDF(K).eq.ROSSTHIN_IDX.or.WHICH_BRDF(K).eq.ROSSTHICK_IDX ) DO_ROSS = .true.
          IF ( WHICH_BRDF(K).eq.LISPARSE_IDX.or.WHICH_BRDF(K).eq.LIDENSE_IDX   ) DO_LI   = .true.
          IF ( WHICH_BRDF(K).eq.LAMBERTIAN_IDX ) DO_Lamb   = .true.
          IF ( WHICH_BRDF(K).eq.COXMUNK_IDX )    DO_CM   = .true.
        ENDDO
        IF ( ( DO_ROSS .and..not. DO_LI ) .or. ( .not.DO_ROSS .and. DO_LI ) ) then
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Ross-type and no Li-type kernel NOT ALLOWED (or vice-versa)'
          ACTIONS(NM)  = 'Re-set input: Ross/Li kernels must be together for MODIS-based BRDFs'
          STATUS = LRRS_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ELSE IF ( ( DO_ROSS .and. DO_LI ) .and. (.not. DO_Lamb .and. .not. do_CM ) ) then
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Ross/Li kernels must go with Lambertian or Cox-Munk'
          ACTIONS(NM)  = 'Re-set input kernels to include Lambertian or Cox-Munk kernel'
          STATUS = LRRS_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  END section. Rob Fix 9/25/14. Check use of MODIS-type kernels
!  **************************************************************

!  Set the Lambertian kernel flags

        DO I = 1, N_BRDF_KERNELS
          IF ( BRDF_NAMES(I) .EQ. 'Lambertian' ) THEN
            LAMBERTIAN_KERNEL_FLAG(I) = .true.
          ENDIF
        ENDDO

!  Shadowing input (for Cox-Munk type)

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' ) THEN
           PAR_STR = 'Do shadow effect for glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_SHADOW_EFFECT
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Continuation point for skipping Regular kernel input

656     continue

!  General inputs, all types
!  -------------------------

!  Number of BRDF azimuth streams, check this value

        PAR_STR = 'Number of BRDF azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) NSTREAMS_BRDF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of  BRDF streams > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS_BRDF dimension'
          STATUS = LRRS_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  **********************************************************************
!  Rob Fix 9/25/14. New Variable and more explicit character string

! FORMER CODE
!  !@@ Overall-Exact flag
!        PAR_STR = 'Do Overall-Exact kernels?'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!           READ (FILUNIT,*,ERR=998)DO_EXACT
!        ENDIF
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!  Exact only flag. Only if above is set (!@@)
!        IF ( DO_EXACT ) THEN
!          PAR_STR = 'Do Exact-only (no Fourier) kernels?'
!          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!             READ (FILUNIT,*,ERR=998)DO_EXACTONLY
!          ENDIF
!          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!        ENDIF
! END FORMER CODE

!  New flag for turning on the Direct-Bounce only flag ("exact" BRDF)
!   Normally this would be .FALSE., if set, then no multiple-scatter BRDFs will be done

        PAR_STR = 'Do direct-bounce only (no multiple-scatter contributions to BRDF)?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_DBONLY
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!   END Replacement Section Rob Fix 9/25/14
! ********************************************************************************

!  Multiple reflectance correction (for Cox-Munk types)
!  ----------------------------------------------------

!  General flag

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' ) THEN
           PAR_STR = 'Do multiple reflectance for all glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_MSRCORR
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Specific MSRCORR inputs
!   Rob Fix 9/25/14. Variable name changed, the wording made more explicit

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Do multiple reflectance for just the direct-bounce glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_MSRCORR_DBONLY
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  MSRCORR scattering order

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering order for glitter kernels'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_ORDER
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  MSRCORR quadrature orders

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering; Polar quadrature order'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_NMUQUAD
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering; Azimuth quadrature order'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_NPHIQUAD
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Check MSCORR dimensions

        IF ( DO_MSRCORR ) THEN
           IF ( MSRCORR_NMUQUAD .gt. max_msrs_muquad ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR polar quadrature No. > Dimensioning'
              ACTIONS(NM)  = 'Increase value of max_msrs_muquad in LRRS_pars'
              STATUS = LRRS_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
           IF ( MSRCORR_NPHIQUAD .gt. max_msrs_phiquad ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR azimuth quadrature No. > Dimensioning'
              ACTIONS(NM)  = 'Increase value of max_msrs_phiquad in LRRS_pars'
              STATUS = LRRS_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Check on MSRCORR order

        IF ( DO_MSRCORR ) THEN
           IF ( MSRCORR_ORDER .EQ.0 ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR is on, but scattering order = 0'
              ACTIONS(NM)  = 'Turn off MSRCORR flags and proceed with warning'
              DO_MSRCORR = .false. ; DO_MSRCORR_DBONLY = .false.
              STATUS = LRRS_WARNING
              NMESSAGES = NM
           ENDIF
        ENDIF

!  White-Sky and Black-Sky Albedo scalings. New for Version 3.7
!  ============================================================

!  Skip this section for NewCM kernel

        if ( DO_NewCMGLINT ) go to 646

!  WSA and BSA scaling options.
!   Revised, 14-15 April 2014, first introduced 02 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.

!  Output option flag.
!  -------------------

!   Revised 12 August 2014; added output option flag.
!   If you are doing WSA or BSA scaling, this should be set automatically

        PAR_STR = 'Do white-sky and black-sky albedo output?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_WSABSA_OUTPUT
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  White-Sky inputs
!  ----------------

!  White-sky Albedo scaling

        PAR_STR = 'Do white-sky albedo scaling?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_WSA_SCALING
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  WSA value. This could be extracted from a data set.....

        IF ( DO_WSA_SCALING  ) THEN
           PAR_STR = 'White-sky albedo value'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              READ (FILUNIT,*,ERR=998)WSA_VALUE
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)
        ENDIF

!  Check WSA value

        IF ( DO_WSA_SCALING  ) THEN
           IF ( WSA_VALUE .le.zero .or. WSA_VALUE .gt. one ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: White-sky albedo value not in the range [0,1]'
              ACTIONS(NM)  = 'Fix the input'
              STATUS = LRRS_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Black-Sky inputs
!  ----------------

!  Black-sky Albedo scaling.

        PAR_STR = 'Do black-sky albedo scaling?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_BSA_SCALING
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Cannot have BSA and WSA together

        IF ( DO_BSA_SCALING .and. DO_WSA_SCALING ) THEN
           NM = NM + 1
           MESSAGES(NM) = 'Bad input: Cannot apply both Black-sky albedo and White-sky albedo scalings!'
           ACTIONS(NM)  = 'Make a choice of which one you want! '
           STATUS = LRRS_SERIOUS
           NMESSAGES = NM
           GOTO 764
        ENDIF

!  BSA value. This could be extracted from a data set.....
!    WARNING: ONLY ALLOWED ONE VALUE HERE...................

        IF ( DO_BSA_SCALING  ) THEN
           PAR_STR = 'Black-sky albedo value'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              READ (FILUNIT,*,ERR=998)BSA_VALUE
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)
        ENDIF

!  Check BSA value

        IF ( DO_BSA_SCALING  ) THEN
           IF ( BSA_VALUE .le.zero .or. BSA_VALUE .gt. one ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: Black-sky albedo value is not in the range [0,1]'
              ACTIONS(NM)  = 'Fix the input'
              STATUS = LRRS_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  continuation point for avoiding the WSA/BSA albedos

646     continue

!  Checking NewCM Kernel. New for Version 3.7
!    Single kernel, Solar Sources only, No Surface Emission. 
!         Scalar only, No MSR (Multiple-surface reflections), No Scaling

        if ( DO_NewCMGLINT ) then

           IF ( WINDSPEED .le.zero ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: NewCM windspeed value is negative'
              ACTIONS(NM)  = 'Fix the input'
              STATUS = LRRS_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF

        endif

!  End BRDF clause

      ENDIF

!  Successful finish

      CLOSE(FILUNIT)

      NMESSAGES = NM

!  Copy Control inputs

      BRDF_Sup_In%BS_DO_USER_STREAMS     = DO_USER_STREAMS
      BRDF_Sup_In%BS_DO_BRDF_SURFACE     = DO_BRDF_SURFACE

!  Copy Geometry results

      BRDF_Sup_In%BS_NSTREAMS          = NSTREAMS
      BRDF_Sup_In%BS_BEAM_SZA          = BEAM_SZA
      BRDF_Sup_In%BS_N_USER_RELAZMS    = N_USER_RELAZMS
      BRDF_Sup_In%BS_USER_RELAZMS      = USER_RELAZMS
      BRDF_Sup_In%BS_N_USER_STREAMS    = N_USER_STREAMS
      BRDF_Sup_In%BS_USER_ANGLES_INPUT = USER_ANGLES

!  Copy BRDF inputs

      BRDF_Sup_In%BS_N_BRDF_KERNELS         = N_BRDF_KERNELS
      BRDF_Sup_In%BS_BRDF_NAMES             = BRDF_NAMES
      BRDF_Sup_In%BS_WHICH_BRDF             = WHICH_BRDF
      BRDF_Sup_In%BS_N_BRDF_PARAMETERS      = N_BRDF_PARAMETERS
      BRDF_Sup_In%BS_BRDF_PARAMETERS        = BRDF_PARAMETERS
      BRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = LAMBERTIAN_KERNEL_FLAG
      BRDF_Sup_In%BS_BRDF_FACTORS           = BRDF_FACTORS
      BRDF_Sup_In%BS_NSTREAMS_BRDF          = NSTREAMS_BRDF

      BRDF_Sup_In%BS_DO_SHADOW_EFFECT       = DO_SHADOW_EFFECT

!  Rob Fix 9/25/14. Two variables replaced
!      BRDF_Sup_In%BS_DO_EXACT               = DO_EXACT         !@@
!      BRDF_Sup_In%BS_DO_EXACTONLY           = DO_EXACTONLY
      BRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY      = DO_DBONLY

      BRDF_Sup_In%BS_DO_GLITTER_MSRCORR        = DO_MSRCORR
      BRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY = DO_MSRCORR_DBONLY  !  Rob Fix 9/25/14, name change
      BRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER     = MSRCORR_ORDER
      BRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD   = MSRCORR_NMUQUAD
      BRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD  = MSRCORR_NPHIQUAD

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 02-07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised 12 August 2014; added output option flag.

      BRDF_Sup_In%BS_DO_WSABSA_OUTPUT = DO_WSABSA_OUTPUT
      BRDF_Sup_In%BS_DO_WSA_SCALING = DO_WSA_SCALING
      BRDF_Sup_In%BS_DO_BSA_SCALING = DO_BSA_SCALING
      BRDF_Sup_In%BS_WSA_VALUE      = WSA_VALUE
      BRDF_Sup_In%BS_BSA_VALUE      = BSA_VALUE

!  NewCM options
!  9/9/15. Wavelengths will be set later, according to the LRRS model. Zeroed here

      BRDF_Sup_In%BS_DO_NewCMGLINT = DO_NewCMGLINT
      BRDF_Sup_In%BS_SALINITY      = SALINITY
      BRDF_Sup_In%BS_DO_NewCM_Wav1 = DO_NewCM_Wav1
!mick fix 10/29/2015 - modified these two (now initialized above with other inputs)
      !BRDF_Sup_In%BS_NewCM_Npoints = 0
      !BRDF_Sup_In%BS_NewCM_Lambdas = zero
      BRDF_Sup_In%BS_NewCM_Npoints = NewCM_Npoints
      BRDF_Sup_In%BS_NewCM_Lambdas = NewCM_Lambdas

      BRDF_Sup_In%BS_WINDSPEED = WINDSPEED
      BRDF_Sup_In%BS_WINDDIR   = WINDDIR

      BRDF_Sup_In%BS_DO_GlintShadow   = DO_GlintShadow
      BRDF_Sup_In%BS_DO_FoamOption    = DO_FoamOption
      BRDF_Sup_In%BS_DO_FacetIsotropy = DO_FacetIsotropy

!  Exception handling

      BRDF_Sup_InputStatus%BS_STATUS_INPUTREAD = STATUS
      BRDF_Sup_InputStatus%BS_NINPUTMESSAGES   = NMESSAGES
      BRDF_Sup_InputStatus%BS_INPUTMESSAGES    = MESSAGES
      BRDF_Sup_InputStatus%BS_INPUTACTIONS     = ACTIONS

!  Normal return

      RETURN

!  Open file error

300   CONTINUE
      STATUS = LRRS_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//trim(adjustl(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
      GO TO 764

!  Line read error - abort immediately

998   CONTINUE
      STATUS = LRRS_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'read failure for '//trim(adjustl(PAR_STR))
      ACTIONS(NMESSAGES)  = 'Re-set: Entry is incorrect in input file'
      CLOSE(FILUNIT)

!  Final error copying

764   CONTINUE

      BRDF_Sup_InputStatus%BS_STATUS_INPUTREAD = STATUS
      BRDF_Sup_InputStatus%BS_NINPUTMESSAGES   = NMESSAGES
      BRDF_Sup_InputStatus%BS_INPUTMESSAGES    = MESSAGES
      BRDF_Sup_InputStatus%BS_INPUTACTIONS     = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE BRDF_INPUTMASTER

!

      SUBROUTINE BRDF_MAINMASTER ( &
        DO_DEBUG_RESTORATION, & ! Inputs
        NMOMENTS_INPUT,       & ! Inputs
        N_Lambdas, Lambdas,   & ! Input wavelengths [nm] for LRRS
        BRDF_Sup_In,          & ! Inputs
        BRDF_Sup_Out,         & ! Outputs
        BRDF_Sup_OutputStatus )  ! Output Status )

!  Prepares the bidirectional reflectance functions necessary for LRRS.

! #####################################################################
! #####################################################################

      USE LRRS_pars_m

      USE brdf_sup_inputs_def_m
      USE brdf_sup_outputs_def_m

      USE brdf_sup_aux_m, only : GETQUAD2, &
                                 BRDF_QUADRATURE_Gaussian, &
                                 BRDF_QUADRATURE_Trapezoid

      USE brdf_sup_kernels_m
      USE brdf_sup_routines_m

!  Implicit none

      IMPLICIT NONE

!  Input variables
!  ---------------

!  Debug flag for restoration

      LOGICAL, intent(in) :: DO_DEBUG_RESTORATION

!  Input number of moments (only used for restoration debug)

      INTEGER, intent(in) :: NMOMENTS_INPUT

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

      REAL(fpk), intent(in) :: Lambdas ( MAX_BRDF_POINTS )

!  Input structure
!  ---------------

      TYPE(BRDF_Sup_Inputs), intent(in)   :: BRDF_Sup_In

!  Output structure
!  ----------------

      TYPE(BRDF_Sup_Outputs), intent(out) :: BRDF_Sup_Out

!  Exception handling introduced 07 April 2014 for Version 3.7

      TYPE(BRDF_Output_Exception_Handling), INTENT(out) :: BRDF_Sup_OutputStatus

!  LRRS local variables
!  ++++++++++++++++++++++

!  Input arguments
!  ===============

!  Stream angle flag

      LOGICAL    :: DO_USER_STREAMS

!   Number and index-list of bidirectional functions

      INTEGER    :: N_BRDF_KERNELS
      INTEGER    :: WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER    :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(fpk)  :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  BRDF names

      CHARACTER (LEN=10) :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Lambertian Surface control

      LOGICAL    :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      REAL(fpk)  :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 02-07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised 12 August 2014; added output option flag

      LOGICAL   :: DO_WSABSA_OUTPUT
      LOGICAL   :: DO_WSA_SCALING
      LOGICAL   :: DO_BSA_SCALING
      REAL(fpk) :: WSA_VALUE, BSA_VALUE

!  Number of azimuth quadrature streams for BRDF

      INTEGER    :: NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL    :: DO_SHADOW_EFFECT

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL   :: DO_EXACT
!      LOGICAL   :: DO_EXACTONLY
      LOGICAL   :: DO_DBONLY

!  Multiple reflectance correction for Glitter kernels

      LOGICAL    :: DO_MSRCORR
      LOGICAL    :: DO_MSRCORR_DBONLY
      INTEGER    :: MSRCORR_ORDER
      INTEGER    :: N_MUQUAD, N_PHIQUAD

!  Local angle control

      INTEGER    :: NSTREAMS
      INTEGER    :: N_USER_STREAMS
      INTEGER    :: N_USER_RELAZMS

!  Solar angle

      REAL(fpk) :: BEAM_SZA

!  Local angles

      REAL(fpk)  :: USER_RELAZMS (MAX_USER_RELAZMS)
      REAL(fpk)  :: USER_ANGLES  (MAX_USER_STREAMS)

!  NewCM Glitter options (bypasses the usual Kernel system)
!  -------------------------------------------------------

!  Overall flags for this option. 9/9/15 Wav1 flag!

      LOGICAL   :: DO_NewCMGLINT, DO_NewCM_Wav1

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR

!  Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Output arguments
!  ================

!  direct bounce BRDF. Rob Fix 9/25/14, name changed from EXACTDB --> DBOUNCE

      REAL(fpk) :: DBOUNCE_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Fourier components of BRDF, in the following order

!    incident solar direction,    reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar direction,    reflected user streams  ! Truncated DB
!    incident quadrature streams, reflected user streams

      REAL(fpk) :: BRDF_F_0 ( 0:MAX_MOMENTS, MAX_STREAMS )
      REAL(fpk) :: BRDF_F   ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS )
      REAL(fpk) :: USER_BRDF_F_0 ( 0:MAX_MOMENTS, MAX_USER_STREAMS  )
      REAL(fpk) :: USER_BRDF_F   ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS )

!  Local BRDF functions
!  ====================

!  at quadrature (discrete ordinate) angles

      REAL(fpk)  :: BRDFUNC   ( MAX_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF )
      REAL(fpk)  :: BRDFUNC_0 ( MAX_STREAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk)  :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF )
      REAL(fpk)  :: USER_BRDFUNC_0 ( MAX_USER_STREAMS, MAXSTREAMS_BRDF )

!  DB Kernel values

      REAL(fpk)  :: DBKERNEL_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS  )

!  Values for WSA/BSA scaling options. New, Version 3.7

      REAL(fpk)  :: SCALING_BRDFUNC   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      REAL(fpk)  :: SCALING_BRDFUNC_0 ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  WSA/BSA scaling componnets, at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: SCALING_BRDF_F   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING )
      DOUBLE PRECISION :: SCALING_BRDF_F_0 ( MAXSTREAMS_SCALING   )

!  Local angles, and cosine/sines/weights
!  ======================================

!  Azimuths

      REAL(fpk)  :: PHIANG(MAX_USER_RELAZMS)
      REAL(fpk)  :: COSPHI(MAX_USER_RELAZMS)
      REAL(fpk)  :: SINPHI(MAX_USER_RELAZMS)

      REAL(fpk)  :: LPHIANG(MAX_USER_RELAZMS)     ! Local value
      REAL(fpk)  :: LCOSPHI(MAX_USER_RELAZMS)     ! Local value
      REAL(fpk)  :: LSINPHI(MAX_USER_RELAZMS)     ! Local value

!  SZAs

      REAL(fpk)  :: SZASURCOS
      REAL(fpk)  :: SZASURSIN

!  Discrete ordinates

      REAL(fpk)  :: QUAD_STREAMS(MAX_STREAMS)
      REAL(fpk)  :: QUAD_WEIGHTS(MAX_STREAMS)
      REAL(fpk)  :: QUAD_SINES  (MAX_STREAMS)

!  Viewing zenith streams

      REAL(fpk)  :: USER_STREAMS(MAX_USER_STREAMS)
      REAL(fpk)  :: USER_SINES  (MAX_USER_STREAMS)

!  BRDF azimuth quadrature streams

      INTEGER    :: NBRDF_HALF
      REAL(fpk)  :: X_BRDF   ( MAXSTREAMS_BRDF )
      REAL(fpk)  :: CX_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk)  :: LCX_BRDF ( MAXSTREAMS_BRDF )     ! Local value
      REAL(fpk)  :: SX_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk)  :: A_BRDF   ( MAXSTREAMS_BRDF )

!  Azimuth factors

      REAL(fpk)  :: BRDF_AZMFAC(MAXSTREAMS_BRDF)

!  Local arrays for MSR quadrature

      REAL(fpk)  :: X_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk)  :: W_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk)  :: SX_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk)  :: WXX_MUQUAD (MAX_MSRS_MUQUAD)

      REAL(fpk)  :: X_PHIQUAD (MAX_MSRS_PHIQUAD)
      REAL(fpk)  :: W_PHIQUAD (MAX_MSRS_PHIQUAD)

!  Local kernel Fourier components
!  ===============================

!  at quadrature (discrete ordinate) angles

      REAL(fpk)  :: LOCAL_BRDF_F   ( MAX_STREAMS, MAX_STREAMS )
      REAL(fpk)  :: LOCAL_BRDF_F_0 ( MAX_STREAMS )

!  at user-defined stream directions

      REAL(fpk)  :: LOCAL_USER_BRDF_F   ( MAX_USER_STREAMS, MAX_STREAMS )
      REAL(fpk)  :: LOCAL_USER_BRDF_F_0 ( MAX_USER_STREAMS )

!  Exception handling. New code, 09 April 2014. Version 3.7
!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )

!  Other local variables
!  =====================

!  Discrete ordinates (local, for Albedo scaling). Version 3.7

      INTEGER     :: SCAL_NSTREAMS
      REAL(fpk)   :: SCAL_QUAD_STREAMS(MAXSTREAMS_SCALING)
      REAL(fpk)   :: SCAL_QUAD_WEIGHTS(MAXSTREAMS_SCALING)
      REAL(fpk)   :: SCAL_QUAD_SINES  (MAXSTREAMS_SCALING)
      REAL(fpk)   :: SCAL_QUAD_STRMWTS(MAXSTREAMS_SCALING)

!  White-sky and Black-sky albedos. Version 3.7.

      LOGICAL    :: DO_LOCAL_WSA, DO_LOCAL_BSA
      REAL(fpk)  :: WSA_CALC (MAX_BRDF_KERNELS), TOTAL_WSA_CALC
      REAL(fpk)  :: BSA_CALC (MAX_BRDF_KERNELS), TOTAL_BSA_CALC

!  Local NewCM variables. Version 3.7.

      REAL(fpk)  :: Refrac_R, Refrac_I, WC_Reflectance, WC_Lambertian, Wavelength

!  Local check of Albedo, for all regular kernel options. Version 3.7.

      LOGICAL   :: DO_CHECK_ALBEDO

!  help

      INTEGER    :: K, B, I, I1, J, UM, IA, M, L
      INTEGER    :: LOCAL_BRDF_NPARS, NMOMENTS, N_PHIQUAD_HALF, LOCAL_NPOINTS
      REAL(fpk)  :: LOCAL_BRDF_PARS ( MAX_BRDF_PARAMETERS ), LOCAL_WAVS( MAX_BRDF_POINTS )
      REAL(fpk)  :: MUX, DELFAC, HELP_A, HELP_B, XM, SCALING, RM
      LOGICAL    :: ADD_FOURIER, LOCAL_MSR

!  Default, use Gaussian quadrature

      LOGICAL, parameter :: DO_BRDFQUAD_GAUSSIAN = .true.

!  Initialize Exception handling
!  -----------------------------

      STATUS = LRRS_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Execution of LRRS BRDF Sup Master'

!  Copy from input type structure
!  ------------------------------

!  Copy Control inputs

      DO_USER_STREAMS     = BRDF_Sup_In%BS_DO_USER_STREAMS
      !DO_BRDF_SURFACE     = BRDF_Sup_In%BS_DO_BRDF_SURFACE

!  Set number of streams

      NSTREAMS = BRDF_Sup_In%BS_NSTREAMS

!  Copy Geometry results

      BEAM_SZA       = BRDF_Sup_In%BS_BEAM_SZA
      N_USER_RELAZMS = BRDF_Sup_In%BS_N_USER_RELAZMS
      USER_RELAZMS   = BRDF_Sup_In%BS_USER_RELAZMS
      N_USER_STREAMS = BRDF_Sup_In%BS_N_USER_STREAMS
      USER_ANGLES    = BRDF_Sup_In%BS_USER_ANGLES_INPUT

!  Copy BRDF inputs

      N_BRDF_KERNELS         = BRDF_Sup_In%BS_N_BRDF_KERNELS
      BRDF_NAMES             = BRDF_Sup_In%BS_BRDF_NAMES
      WHICH_BRDF             = BRDF_Sup_In%BS_WHICH_BRDF
      N_BRDF_PARAMETERS      = BRDF_Sup_In%BS_N_BRDF_PARAMETERS
      BRDF_PARAMETERS        = BRDF_Sup_In%BS_BRDF_PARAMETERS
      LAMBERTIAN_KERNEL_FLAG = BRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG
      BRDF_FACTORS           = BRDF_Sup_In%BS_BRDF_FACTORS
      NSTREAMS_BRDF          = BRDF_Sup_In%BS_NSTREAMS_BRDF
      DO_SHADOW_EFFECT       = BRDF_Sup_In%BS_DO_SHADOW_EFFECT

!  Rob Fix 9/25/14. Two variables replaced

      DO_DBONLY      = BRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 02-07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised 12 August 2014, Added WSBSA output flag option

      DO_WSABSA_OUTPUT    = BRDF_Sup_In%BS_DO_WSABSA_OUTPUT
      DO_WSA_SCALING      = BRDF_Sup_In%BS_DO_WSA_SCALING
      DO_BSA_SCALING      = BRDF_Sup_In%BS_DO_BSA_SCALING
      WSA_VALUE           = BRDF_Sup_In%BS_WSA_VALUE
      BSA_VALUE           = BRDF_Sup_In%BS_BSA_VALUE

!  NewCM options. 9/9/15. Single wavelength point.

      DO_NewCMGLINT = BRDF_Sup_In%BS_DO_NewCMGLINT
      DO_NewCM_Wav1 = BRDF_Sup_In%BS_DO_NewCM_Wav1

      SALINITY      = BRDF_Sup_In%BS_SALINITY
      WINDSPEED = BRDF_Sup_In%BS_WINDSPEED
      WINDDIR   = BRDF_Sup_In%BS_WINDDIR

      DO_GlintShadow   = BRDF_Sup_In%BS_DO_GlintShadow
      DO_FoamOption    = BRDF_Sup_In%BS_DO_FoamOption
      DO_FacetIsotropy = BRDF_Sup_In%BS_DO_FacetIsotropy 

!  Local check of albedo

      DO_CHECK_ALBEDO = .not. DO_NewCMGLINT

!  Local flags
!   Revised 12 August 2014; added output option flag

      DO_LOCAL_WSA = ( DO_WSA_SCALING .or. DO_WSABSA_OUTPUT ) .or.  DO_CHECK_ALBEDO
      DO_LOCAL_BSA = ( DO_BSA_SCALING .or. DO_WSABSA_OUTPUT )

!  Copy MSR inputs

      DO_MSRCORR          = BRDF_Sup_In%BS_DO_GLITTER_MSRCORR
      DO_MSRCORR_DBONLY   = BRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY
      MSRCORR_ORDER       = BRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER
      N_MUQUAD            = BRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD
      N_PHIQUAD           = BRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD

!  Main code
!  ---------

!  Set up Quadrature streams for output
!    QUAD_STRMWTS dropped for Version 3.7 (now redefined for local WSA/BSA scaling)

      CALL GETQUAD2 ( ZERO, ONE, NSTREAMS, QUAD_STREAMS, QUAD_WEIGHTS )
      DO I = 1, NSTREAMS
        QUAD_SINES(I) = SQRT(ONE-QUAD_STREAMS(I)*QUAD_STREAMS(I))
      enddo

!  Set up Quadrature streams for WSA/BSA Scaling. New code, Version 3.7

      IF ( DO_LOCAL_WSA .or. DO_LOCAL_BSA ) THEN
         SCAL_NSTREAMS = MAXSTREAMS_SCALING
         CALL GETQUAD2 ( ZERO, ONE, SCAL_NSTREAMS, SCAL_QUAD_STREAMS, SCAL_QUAD_WEIGHTS )
         DO I = 1, SCAL_NSTREAMS
            SCAL_QUAD_SINES(I)   = SQRT(ONE-SCAL_QUAD_STREAMS(I)*SCAL_QUAD_STREAMS(I))
            SCAL_QUAD_STRMWTS(I) = SCAL_QUAD_STREAMS(I) * SCAL_QUAD_WEIGHTS(I)
         enddo
      ENDIF

!  Number of Fourier components to calculate

      IF ( DO_DEBUG_RESTORATION ) THEN
        NMOMENTS = NMOMENTS_INPUT
      ELSE
        NMOMENTS = 2 * NSTREAMS - 1
      ENDIF

!  Half number of moments

      NBRDF_HALF = NSTREAMS_BRDF / 2

!    Warning, this should be the BOA angle. OK for the non-refractive case.

      MUX =  COS(BEAM_SZA*DEG_TO_RAD)
      SZASURCOS = MUX
      SZASURSIN = SQRT(1.0_fpk-MUX*MUX)

!  Viewing angles

      DO UM = 1, N_USER_STREAMS
        USER_STREAMS(UM) = COS(USER_ANGLES(UM)*DEG_TO_RAD)
        USER_SINES(UM)   = SQRT(ONE-USER_STREAMS(UM)*USER_STREAMS(UM))
      ENDDO

!   Rob Fix 9/25/14. Removed DO_EXACT; Contribution always required now....

        DO IA = 1, N_USER_RELAZMS
          PHIANG(IA) = USER_RELAZMS(IA)*DEG_TO_RAD
          COSPHI(IA) = COS(PHIANG(IA))
          SINPHI(IA) = SIN(PHIANG(IA))
        ENDDO

!  BRDF quadrature
!  ---------------

!  Save these quantities for efficient coding.
!      Surface emission I/O removed, 9/9/15

      IF ( DO_BRDFQUAD_GAUSSIAN ) then
        CALL BRDF_QUADRATURE_Gaussian &
        ( NSTREAMS_BRDF, NBRDF_HALF,         & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF )   ! Outputs
      ELSE
        CALL BRDF_QUADRATURE_Trapezoid &
        ( NSTREAMS_BRDF,                     & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF )   ! Outputs
      ENDIF

!  Set up the MSR points
!  ---------------------

!  Air to water, Polar quadrature

      IF ( DO_MSRCORR  ) THEN
         CALL GETQUAD2 ( ZERO, ONE, N_MUQUAD, X_MUQUAD, W_MUQUAD )
         DO I = 1, N_MUQUAD
            XM = X_MUQUAD(I)
            SX_MUQUAD(I) = SQRT(ONE-XM*XM)
            WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
         ENDDO
      ENDIF

!  Azimuth quadrature

      IF ( DO_MSRCORR  ) THEN
         N_PHIQUAD_HALF = N_PHIQUAD / 2
         CALL GETQUAD2 ( ZERO, ONE, N_PHIQUAD_HALF, X_PHIQUAD, W_PHIQUAD )
         DO I = 1, N_PHIQUAD_HALF
           I1 = I + N_PHIQUAD_HALF
           X_PHIQUAD(I1) = - X_PHIQUAD(I)
           W_PHIQUAD(I1) =   W_PHIQUAD(I)
         ENDDO
         DO I = 1, N_PHIQUAD
            X_PHIQUAD(I)  = PIE * X_PHIQUAD(I)
         ENDDO
      ENDIF

!  Rob Fix. 9/9/15. Initialize the wavelength array
!    Generally only one point, unless control has been set for NewCM
!      - the 0.001 factor is conversion to Microns.

      Local_npoints = 1 ; local_wavs = zero
      if ( DO_NewCMGLINT ) then
         if ( DO_NewCM_Wav1 ) then
            local_wavs(1) = 0.001_fpk * SUM(lambdas(1:N_lambdas))/DBLE(N_Lambdas)
         else
            Local_npoints = N_Lambdas
            local_wavs(1:N_lambdas) = 0.001_fpk * lambdas(1:N_lambdas)
         endif
      endif

!  Rob Fix. 9/9/15. Initialize the output array

      BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC = zero
      BRDF_Sup_Out%BS_BRDF_F_0        = zero
      BRDF_Sup_Out%BS_BRDF_F          = zero
      BRDF_Sup_Out%BS_USER_BRDF_F_0   = zero
      BRDF_Sup_Out%BS_USER_BRDF_F     = zero
      BRDF_Sup_Out%BS_WSA_CALCULATED  = zero
      BRDF_Sup_Out%BS_BSA_CALCULATED  = zero

!  Start wavelength loop

      do L = 1, Local_npoints

!mick fix 6/29/11 - initialize ALL elements of BRDF arrays

!  Initialise BRDF arrays
!  ----------------------

!   Rob Fix 9/25/14. Direct-bounce name DBOUNCE, replaces EXACTDB

         DBOUNCE_BRDFUNC = ZERO

         BRDF_F_0        = ZERO
         BRDF_F          = ZERO
         USER_BRDF_F_0   = ZERO
         USER_BRDF_F     = ZERO

         SCALING_BRDFUNC   = zero
         SCALING_BRDFUNC_0 = zero

!  Initialize WSA/BSA albedos

         WSA_CALC = zero ; TOTAL_WSA_CALC = zero
         BSA_CALC = zero ; TOTAL_BSA_CALC = zero

!  Fill BRDF arrays
!  ----------------

         DO K = 1, N_BRDF_KERNELS

!  Local variables

           LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
           DO B = 1, MAX_BRDF_PARAMETERS
             LOCAL_BRDF_PARS(B) = BRDF_PARAMETERS(K,B)
           ENDDO

!  Coxmunk shadow flag

           IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
             IF ( DO_SHADOW_EFFECT ) LOCAL_BRDF_PARS(3) = ONE
           ENDIF

!  Local MSRCORR flag

           LOCAL_MSR = .false.
           IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
             LOCAL_MSR = DO_MSRCORR
           ENDIF

!  Get the kernels (all but NewCM)

           IF ( WHICH_BRDF(K) .ne. NewCMGLINT_IDX ) THEN
             CALL BRDF_MAKER &
              ( DO_LOCAL_WSA, DO_LOCAL_BSA, DO_USER_STREAMS,            & ! New line, Version 3.7
                WHICH_BRDF(K), DO_DBONLY, LOCAL_MSR, DO_MSRCORR_DBONLY, & ! Inputs
                MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                     & ! Inputs
                LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, NSTREAMS_BRDF,       & ! Inputs
                NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
                QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,     & ! Inputs
                SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,           & ! Inputs
                SCAL_NSTREAMS, SCAL_QUAD_STREAMS, SCAL_QUAD_SINES,      & ! New line, Version 3.7
                X_BRDF, CX_BRDF, SX_BRDF,                               & ! Inputs
                X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,              & ! Inputs
                X_PHIQUAD, W_PHIQUAD,                                   & ! Inputs
                DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                & ! Outputs
                BRDFUNC_0, USER_BRDFUNC_0,                              & ! Output
                SCALING_BRDFUNC, SCALING_BRDFUNC_0  )                     ! output, New line, Version 3.7 
           ENDIF

!  NewCM Kernel. New for Version 3.7
!  ---------------------------------

!    Single kernel, Solar Sources only, No Surface Emission. 
!         Scalar only, No MSR (Multiple-surface reflections), No Scaling

!  Sequence is (1) Salinity, WhiteCap, Cox-Munk

           IF ( WHICH_BRDF(K) .EQ. NewCMGLINT_IDX ) THEN

!  Reverse angle effect !!!!!
!   NewCM Convention is opposite from LRRS, that is, Phi(6S) = 180 - Phi(L)
!   Once this is realized, NewCM and Regular Cox-Munk will agree perfectly
!          ( Have to turn off Whitecaps in 6S and use Facet Isotropy, make sure RI same)
!   Rob Fix 9/25/14. Removed DO_EXACT

             DO IA = 1, N_USER_RELAZMS
               LPHIANG(IA) = PIE - PHIANG(IA)
               LCOSPHI(IA) = - COSPHI(IA)
               LSINPHI(IA) = + SINPHI(IA)
             ENDDO
             DO I = 1, NSTREAMS_BRDF
               LCX_BRDF(I) = - CX_BRDF(I)
             ENDDO

!  Set local wavelength

             Wavelength = Local_wavs(L)

!  Refractive index. Formerly INDWAT

             Call BRDF_Water_RefracIndex  ( Wavelength, Salinity, Refrac_R, Refrac_I )

!  tempo (fake Cox-Munk)
!          Refrac_R = 1.334d0 ; Refrac_I = 0.0d0
!          Do_FoamOption = .false. ; DO_FacetIsotropy= .true.

!  Foam-reflectance correction.

             WC_Reflectance = zero ; WC_Lambertian = zero
             if ( Do_FoamOption ) then
               call BRDF_WhiteCap_Reflectance &
                 ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian )
             endif

!  Make the reflectance. Includes the WhiteCap term.

             CALL BRDF_NewCM_MAKER &
             ( DO_GlintShadow, DO_FacetIsotropy, WINDSPEED, WINDDIR,            &
               Refrac_R, Refrac_I, WC_Reflectance, WC_Lambertian,               &
               DO_USER_STREAMS, DO_DBONLY,                                      &
               NSTREAMS_BRDF, NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS,         &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,              &
               SZASURCOS, SZASURSIN, LPHIANG, LCOSPHI, LSINPHI,                 &
               X_BRDF, LCX_BRDF, SX_BRDF,                                       &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                         & ! output
               BRDFUNC_0, USER_BRDFUNC_0  )                                       ! output

           ENDIF

!  Compute Exact Direct Beam BRDF
!  ==============================

           DO IA = 1, N_USER_RELAZMS
             DO UM = 1, N_USER_STREAMS
               DBOUNCE_BRDFUNC(UM,IA) = DBOUNCE_BRDFUNC(UM,IA) &
                  + BRDF_FACTORS(K) * DBKERNEL_BRDFUNC(UM,IA)
             ENDDO
           ENDDO

!  Scaling Section. New code, 15 April 2014 for Version 3.7
!  ========================================================

!  Get the requisite Fourier 0 components

           IF ( DO_LOCAL_WSA .or. DO_LOCAL_BSA ) THEN
             CALL SCALING_FOURIER_ZERO &
              ( DO_LOCAL_WSA, DO_LOCAL_BSA, LAMBERTIAN_KERNEL_FLAG(K), &
                SCAL_NSTREAMS, NSTREAMS_BRDF,                          &
                A_BRDF, SCALING_BRDFUNC, SCALING_BRDFUNC_0,            &
                SCALING_BRDF_F, SCALING_BRDF_F_0 )
           ENDIF

!  White-sky Spherical albedo. Code Upgraded for Version 3.7
!  ---------------------------------------------------------

           IF ( DO_LOCAL_WSA ) THEN

!  Only for non-Lambertian kernels (trivially = 1 otherwise)

             WSA_CALC(K) = ONE
             IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
               HELP_A = ZERO
               DO I = 1, SCAL_NSTREAMS
                 HELP_B = DOT_PRODUCT(SCALING_BRDF_F(I,1:SCAL_NSTREAMS),SCAL_QUAD_STRMWTS(1:SCAL_NSTREAMS))
                 HELP_A = HELP_A + HELP_B * SCAL_QUAD_STRMWTS(I)
               ENDDO
               WSA_CALC(K) = HELP_A * FOUR
             ENDIF
             TOTAL_WSA_CALC = TOTAL_WSA_CALC + BRDF_FACTORS(K) * WSA_CALC(K)

!  Perform consistency check on total white-sky spherical albedo
!    -- This is only done after the kernel summation is finished
!    -- If failed, go to 899 and the error output

             IF ( K.eq.N_BRDF_KERNELS ) then
               if ( TOTAL_WSA_CALC .le. zero ) then
                 STATUS = LRRS_SERIOUS ; NMESSAGES = NMESSAGES + 1
                 MESSAGES(NMESSAGES) = 'Fatal error: Total White-sky albedo is Negative; examine BRDF Amplitudes'
               else if ( TOTAL_WSA_CALC .gt. one ) then
                 STATUS = LRRS_SERIOUS ; NMESSAGES = NMESSAGES + 1
                 MESSAGES(NMESSAGES) = 'Fatal error: Total White-sky albedo is > 1; examine BRDF Amplitudes'
               endif
               IF (STATUS.NE.LRRS_success) GO TO 899
             endif

           ENDIF

!  Black-sky Albedo, only for 1 solar beam. Code Upgraded for Version 3.7
!  ---------------------------------------

!  Compute it for non-Lambertian kernels
!     No check necessary, as the WSA is always checked (regardless of whether scaling is applied)

           IF (  DO_LOCAL_BSA ) THEN
             BSA_CALC(K) = ONE
             IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
               BSA_CALC(K) = TWO * DOT_PRODUCT(SCALING_BRDF_F_0(1:SCAL_NSTREAMS),SCAL_QUAD_STRMWTS(1:SCAL_NSTREAMS))
             ENDIF
             TOTAL_BSA_CALC = TOTAL_BSA_CALC + BRDF_FACTORS(K) * BSA_CALC(K)
           ENDIF

!  @@. Skip Fourier section, if Only the Direct bounce term is calculated.
!   Rob Fix 9/25/14. Direct-bounce flag DO_DBONLY, replaces DO_EXACTONLY

           IF ( DO_DBONLY ) go to 676

!  Fourier Work now
!  ================

           DO M = 0, NMOMENTS

!  Fourier addition flag

             ADD_FOURIER = ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) .OR. &
                              (LAMBERTIAN_KERNEL_FLAG(K) .AND. M.EQ.0) )

!  surface reflectance factors, Weighted Azimuth factors

             IF ( M .EQ. 0 ) THEN
               DELFAC   = ONE
               DO I = 1, NSTREAMS_BRDF
                 BRDF_AZMFAC(I) = A_BRDF(I)
               ENDDO
             ELSE
               DELFAC   = TWO ; RM = REAL(M,fpk)
               DO I = 1, NSTREAMS_BRDF
                 BRDF_AZMFAC(I) = A_BRDF(I) * COS ( RM * X_BRDF(I) )
               ENDDO
             ENDIF

!  Call

             CALL BRDF_FOURIER &
             ( DO_USER_STREAMS, LAMBERTIAN_KERNEL_FLAG(K), M, DELFAC,         & ! Inputs
               NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF,                       & ! Inputs
               BRDFUNC, USER_BRDFUNC, BRDFUNC_0, USER_BRDFUNC_0, BRDF_AZMFAC, & ! Inputs
               LOCAL_BRDF_F, LOCAL_BRDF_F_0,                                  & ! Outputs
               LOCAL_USER_BRDF_F, LOCAL_USER_BRDF_F_0 )                         ! Outputs

!  Start Fourier addition

             IF ( ADD_FOURIER ) THEN

!  Kernel combinations (for quadrature-quadrature reflectance)
!   !@@ Code separated

               DO I = 1, NSTREAMS
                 DO J = 1, NSTREAMS
                   BRDF_F(M,I,J) = BRDF_F(M,I,J) &
                      + BRDF_FACTORS(K) * LOCAL_BRDF_F(I,J)
                 ENDDO
               ENDDO

!  Kernel combinations (for Solar-quadrature reflectance)

               DO I = 1, NSTREAMS
                 BRDF_F_0(M,I) = BRDF_F_0(M,I) &
                       + BRDF_FACTORS(K) * LOCAL_BRDF_F_0(I)
               ENDDO

!  Kernel combinations (for Quadrature-to-Userstream reflectance)
!   !@@ Code separated 12/31/12

               IF ( DO_USER_STREAMS ) THEN
                 DO UM = 1, N_USER_STREAMS
                   DO J = 1, NSTREAMS
                     USER_BRDF_F(M,UM,J) = USER_BRDF_F(M,UM,J) &
                        + BRDF_FACTORS(K) * LOCAL_USER_BRDF_F(UM,J)
                   ENDDO
                 ENDDO
               ENDIF

!  Kernel combinations (for Solar-to-Userstream reflectance)
!   !@@ Generally only required for a MS + SS Truncated calculation
!   This is the Truncated Direct bounce calculation. Always be made available.

               IF ( DO_USER_STREAMS ) THEN
                 DO UM = 1, N_USER_STREAMS
                   USER_BRDF_F_0(M,UM) = USER_BRDF_F_0(M,UM) &
                          + BRDF_FACTORS(K) * LOCAL_USER_BRDF_F_0(UM)
                 ENDDO
               ENDIF

!  End Fourier addition and loop

             ENDIF
           ENDDO

!  continuation point for skipping Fourier work. !@@

676        continue

!  End kernel loop

         ENDDO

!  Now perform normalizations and scaling with White-sky or Black-sky albedos. New section, 02-15 April 2014
!  =========================================================================================================

!  SCALING only if flagged.

         IF ( DO_WSA_SCALING .or. DO_BSA_SCALING ) THEN

!  set scaling factor

           if ( DO_WSA_SCALING ) then
             SCALING = WSA_VALUE / TOTAL_WSA_CALC
           else
             SCALING = BSA_VALUE / TOTAL_BSA_CALC
           endif

!  Exact Direct Beam BRDF
           DBOUNCE_BRDFUNC(1:N_USER_STREAMS,1:N_USER_RELAZMS) = &
                 SCALING * DBOUNCE_BRDFUNC(1:N_USER_STREAMS,1:N_USER_RELAZMS)

!  Fourier terms
!     Rob Fix 9/25/14. Change name to DBONLY

           IF ( .not. DO_DBONLY ) THEN
             DO M = 0, NMOMENTS
!        quadrature-quadrature reflectance
!        Solar-quadrature  reflectance
               BRDF_F(M,1:nstreams,1:nstreams) = SCALING * BRDF_F(M,1:nstreams,1:nstreams)
               BRDF_F_0(M,1:nstreams) = SCALING * BRDF_F_0(M,1:nstreams) 
!        Quadrature-to-Userstream reflectance
!        Solar-to-Userstream      reflectance
               IF ( DO_USER_STREAMS ) THEN
                 USER_BRDF_F(M,1:N_USER_STREAMS,1:nstreams) = SCALING * USER_BRDF_F(M,1:N_USER_STREAMS,1:nstreams) 
                 USER_BRDF_F_0(M,1:N_USER_STREAMS)          = SCALING * USER_BRDF_F_0(M,1:N_USER_STREAMS)
               ENDIF
             ENDDO
           ENDIF

!  End WSA/BSA scaling option

         ENDIF

!  Copy to output structure

         BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(:,:,L) = DBOUNCE_BRDFUNC(:,:)

         BRDF_Sup_Out%BS_BRDF_F_0(:,:,L)        = BRDF_F_0(:,:)
         BRDF_Sup_Out%BS_BRDF_F  (:,:,:,L)      = BRDF_F  (:,:,:)
         BRDF_Sup_Out%BS_USER_BRDF_F_0(:,:,L)   = USER_BRDF_F_0(:,:)
         BRDF_Sup_Out%BS_USER_BRDF_F  (:,:,:,L) = USER_BRDF_F  (:,:,:)

!  WSABSA OUTPUT
!  -------------

!  Revision 12 August 2014. Added output option.

         IF ( DO_WSABSA_OUTPUT ) THEN
           BRDF_Sup_Out%BS_WSA_CALCULATED(L) = TOTAL_WSA_CALC
           BRDF_Sup_Out%BS_BSA_CALCULATED(L) = TOTAL_BSA_CALC
         ENDIF

!  End points loop

      enddo

!  Continuation point for Error Finish from Consistency Check of WSA or BSA calculated values

899   continue

!  write Exception handling to output structure

      BRDF_Sup_OutputStatus%BS_STATUS_OUTPUT   = STATUS
      BRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES = NMESSAGES
      BRDF_Sup_OutputStatus%BS_OUTPUTMESSAGES  = MESSAGES

!  Finish

      RETURN
      END SUBROUTINE BRDF_MAINMASTER

      END MODULE brdf_sup_masters_m

