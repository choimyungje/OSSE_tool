
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
! #            LRRS_BRDF_INPUT_CHECKER                          #
! #            LRRS_SLEAVE_INPUT_CHECKER                        #
! #            BRDF_SLEAVE_INPUT_CHECKER                        #
! #                                                             #
! ###############################################################

!  Upgrade to LRRS Version 2.5
!  ---------------------------

!    ** The BRDF and Sleave "Input_Checker" subroutines listed here
!       insure consistency between certain BRDF and SLEAVE supplement
!       inputs and the corresponding LRRS inputs

!    ** A BRDF_SLEAVE_INPUT_CHECKER routine is present to insure 
!       consistency between BRDF and Sleave supplement inputs when the 
!       "New CM" ocean-treatment is in use

      MODULE lrrs_sup_accessories_1_m

      PRIVATE
      PUBLIC :: LRRS_BRDF_INPUT_CHECKER,   &
                LRRS_SLEAVE_INPUT_CHECKER, &
                BRDF_SLEAVE_INPUT_CHECKER

      CONTAINS

      SUBROUTINE LRRS_BRDF_INPUT_CHECKER ( &
        BRDF_Sup_In,          & ! Inputs
        LRRS_FixIn,           & ! Inputs
        LRRS_ModIn,           & ! Inputs
        LRRS_BRDFCheck_Status ) ! Outputs

!  This subroutine checks certain LRRS BRDF supplement inputs with the
!  corresponding LRRS main code inputs for consistency 

!  Modules for use

      USE LRRS_PARS_m
      USE BRDF_Sup_Inputs_def_m
      USE LRRS_Inputs_def_m
      USE LRRS_Outputs_def_m

      IMPLICIT NONE

      TYPE(BRDF_Sup_Inputs), INTENT(IN)           :: BRDF_Sup_In

      TYPE(LRRS_Fixed_Inputs), INTENT (IN)        :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT (IN)     :: LRRS_ModIn

      TYPE(LRRS_Exception_Handling), INTENT(OUT)  :: LRRS_BRDFCheck_Status

!  ---------------
!  Local variables
!  ---------------

!  BRDF supplement inputs
!  ----------------------

!  User stream, BRDF surface and surface emission flags

      LOGICAL ::             BS_DO_USER_STREAMS
      LOGICAL ::             BS_DO_BRDF_SURFACE

!  Number of discrete ordinate streams

      INTEGER ::             BS_NSTREAMS

!  BOA solar zenith angles

      REAL(FPK) ::           BS_BEAM_SZA

!  User-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER ::             BS_N_USER_RELAZMS
      REAL(FPK) ::           BS_USER_RELAZMS (MAX_USER_RELAZMS)

!  User-defined zenith angle input

      INTEGER ::             BS_N_USER_STREAMS
      REAL(FPK) ::           BS_USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  LRRS Main inputs
!  -------------------

!  LRRS_Fixed_Boolean

      LOGICAL ::             DO_MVOUT_ONLY
      LOGICAL ::             DO_BRDF_SURFACE

!  LRRS_Fixed_Control

      INTEGER ::             NSTREAMS
      REAL(FPK) ::           SOLAR_ANGLE

!  LRRS_Fixed_UserValues

      INTEGER ::             N_USER_STREAMS
      REAL(FPK) ::           USER_ANGLES ( MAX_USER_STREAMS )
      INTEGER ::             N_USER_RELAZMS
      REAL(FPK) ::           USER_RELAZMS ( MAX_USER_RELAZMS )

!  LRRS_Modified_Boolean

      LOGICAL ::             DO_USER_STREAMS

!  Exception handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS  ( 0:MAX_MESSAGES )

!  Other

      INTEGER          :: NM, I
      CHARACTER(Len=2) :: C2

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  BRDF Control inputs

      BS_DO_USER_STREAMS     = BRDF_Sup_In%BS_DO_USER_STREAMS
      BS_DO_BRDF_SURFACE     = BRDF_Sup_In%BS_DO_BRDF_SURFACE

!  BRDF Geometry inputs

      BS_NSTREAMS            = BRDF_Sup_In%BS_NSTREAMS
      BS_BEAM_SZA            = BRDF_Sup_In%BS_BEAM_SZA
      BS_N_USER_RELAZMS      = BRDF_Sup_In%BS_N_USER_RELAZMS
      BS_USER_RELAZMS        = BRDF_Sup_In%BS_USER_RELAZMS
      BS_N_USER_STREAMS      = BRDF_Sup_In%BS_N_USER_STREAMS
      BS_USER_ANGLES_INPUT   = BRDF_Sup_In%BS_USER_ANGLES_INPUT

!  LRRS Fixed Boolean inputs

      DO_MVOUT_ONLY          = LRRS_FixIn%Bool%DO_MVOUT_ONLY
      DO_BRDF_SURFACE        = LRRS_FixIn%Bool%DO_BRDF_SURFACE

!  LRRS Fixed Control inputs

      NSTREAMS          = LRRS_FixIn%Cont%NSTREAMS
      SOLAR_ANGLE       = LRRS_FixIn%Cont%SOLAR_ANGLE

!  LRRS Fixed User Value inputs

      N_USER_STREAMS    = LRRS_FixIn%UserVal%N_USER_STREAMS
      USER_ANGLES       = LRRS_FixIn%UserVal%USER_ANGLES

      N_USER_RELAZMS    = LRRS_FixIn%UserVal%N_USER_RELAZMS
      USER_RELAZMS      = LRRS_FixIn%UserVal%USER_RELAZMS

!  LRRS Modified Boolean inputs

      DO_USER_STREAMS   = LRRS_ModIn%MBool%DO_USER_STREAMS

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Initialize output status

      STATUS_INPUTCHECK = LRRS_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES   = 0
      MESSAGES(0) = 'Successful Check of BRDF/MAIN compatibility'
      ACTIONS(0)  = 'No Action required for this Task'

      NM = NMESSAGES

!  Checks

      IF ( BS_DO_BRDF_SURFACE .neqv. DO_BRDF_SURFACE ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'BRDF surface not set for LRRS Main'
        ACTIONS(NM)  = 'DO_LAMBERTIAN_SURFACE should be False!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( BS_DO_USER_STREAMS .neqv. DO_USER_STREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( BS_DO_USER_STREAMS .neqv. (.not.DO_MVOUT_ONLY) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams and DO_MVOUT_ONLY do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( BS_NSTREAMS .ne. NSTREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of discrete ordinates does not agree'
        ACTIONS(NM)  = 'Check NSTREAMS input'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

!  Angles

      IF ( BS_BEAM_SZA .ne. SOLAR_ANGLE ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Solar beam angle does not agree'
        ACTIONS(NM)  = 'Check BS_BEAM_SZA and SOLAR_ANGLE input'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( BS_N_USER_STREAMS .ne. N_USER_STREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check N_USER_STREAMS and N_USER_STREAMS input'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ELSE
        DO I = 1, N_USER_STREAMS
          IF ( BS_USER_ANGLES_INPUT(I) .ne. USER_ANGLES(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'View zenith angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check BS_USER_ANGLES_INPUT & USER_ANGLES input'
            STATUS_INPUTCHECK = LRRS_SERIOUS
          ENDIF
        ENDDO
      ENDIF

      IF ( BS_N_USER_RELAZMS .ne. N_USER_RELAZMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check BS_N_USER_RELAZMS & N_USER_RELAZMS input'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ELSE
        DO I = 1, N_USER_RELAZMS
          if ( BS_USER_RELAZMS(I) .ne.USER_RELAZMS(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Azimuth angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check BS_USER_RELAZMS & USER_RELAZMS input'
            STATUS_INPUTCHECK = LRRS_SERIOUS
          endif
        ENDDO
      ENDIF

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      LRRS_BRDFCheck_Status%STATUS_INPUTCHECK = STATUS_INPUTCHECK
      LRRS_BRDFCheck_Status%NCHECKMESSAGES    = NMESSAGES
      LRRS_BRDFCheck_Status%CHECKMESSAGES     = MESSAGES
      LRRS_BRDFCheck_Status%ACTIONS           = ACTIONS

!  Finish

      END SUBROUTINE LRRS_BRDF_INPUT_CHECKER

!

      SUBROUTINE LRRS_SLEAVE_INPUT_CHECKER ( &
        SLEAVE_Sup_In,          & ! Inputs
        LRRS_FixIn,             & ! Inputs
        LRRS_ModIn,             & ! Inputs
        LRRS_SLEAVECheck_Status ) ! Outputs

!  This subroutine checks certain LRRS SLEAVE supplement inputs with the
!  corresponding LRRS main code inputs for consistency 

!  Modules for use

      USE LRRS_PARS_m
      USE SLEAVE_Sup_Inputs_def_m
      USE LRRS_Inputs_def_m
      USE LRRS_Outputs_def_m

      IMPLICIT NONE

      TYPE(SLEAVE_Sup_Inputs), INTENT(IN)         :: SLEAVE_Sup_In

      TYPE(LRRS_Fixed_Inputs), INTENT (IN)        :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT (IN)     :: LRRS_ModIn

      TYPE(LRRS_Exception_Handling), INTENT(OUT)  :: LRRS_SLEAVECheck_Status

!  ---------------
!  Local variables
!  ---------------

!  SLEAVE supplement inputs
!  -------------------------

!  Surface-leaving control flags

      LOGICAL :: SL_DO_SLEAVING
      LOGICAL :: SL_DO_ISOTROPIC
      LOGICAL :: SL_DO_EXACTONLY
      LOGICAL :: SL_DO_USER_STREAMS

!  Number of discrete ordinate streams

      INTEGER ::      SL_NSTREAMS

!  BOA solar zenith angles

      REAL(FPK) ::    SL_BEAM_SZA

!  User-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER ::      SL_N_USER_RELAZMS
      REAL(FPK) ::    SL_USER_RELAZMS (MAX_USER_RELAZMS)

!  User-defined zenith angle input

      INTEGER ::      SL_N_USER_STREAMS
      REAL(FPK) ::    SL_USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  LRRS Main inputs
!  -------------------

!  LRRS_Fixed_Boolean

      LOGICAL ::      DO_MVOUT_ONLY
      LOGICAL ::      DO_SURFACE_LEAVING
      LOGICAL ::      DO_SL_ISOTROPIC

!  LRRS_Fixed_Control

      INTEGER ::      NSTREAMS
      REAL(FPK) ::    SOLAR_ANGLE

!  LRRS_Fixed_UserValues

      INTEGER ::      N_USER_STREAMS
      REAL(FPK) ::    USER_ANGLES ( MAX_USER_STREAMS )
      INTEGER ::      N_USER_RELAZMS
      REAL(FPK) ::    USER_RELAZMS ( MAX_USER_RELAZMS )

!  LRRS_Modified_Boolean

      LOGICAL ::      DO_USER_STREAMS

!  Exception handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS  ( 0:MAX_MESSAGES )

!  Other

      INTEGER          :: NM, I
      CHARACTER(Len=2) :: C2

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  SLEAVE Control inputs

      SL_DO_SLEAVING         = SLEAVE_Sup_In%SL_DO_SLEAVING
      SL_DO_ISOTROPIC        = SLEAVE_Sup_In%SL_DO_ISOTROPIC
      SL_DO_EXACTONLY        = SLEAVE_Sup_In%SL_DO_EXACTONLY
      SL_DO_USER_STREAMS     = SLEAVE_Sup_In%SL_DO_USER_STREAMS

!  SLEAVE Geometry inputs

      SL_NSTREAMS            = SLEAVE_Sup_In%SL_NSTREAMS
      SL_BEAM_SZA            = SLEAVE_Sup_In%SL_BEAM_SZA
      SL_N_USER_RELAZMS      = SLEAVE_Sup_In%SL_N_USER_RELAZMS
      SL_USER_RELAZMS        = SLEAVE_Sup_In%SL_USER_RELAZMS
      SL_N_USER_STREAMS      = SLEAVE_Sup_In%SL_N_USER_STREAMS
      SL_USER_ANGLES_INPUT   = SLEAVE_Sup_In%SL_USER_ANGLES_INPUT

!  LRRS Fixed Boolean inputs

      DO_MVOUT_ONLY      = LRRS_FixIn%Bool%DO_MVOUT_ONLY
      DO_SURFACE_LEAVING = LRRS_FixIn%Bool%DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC    = LRRS_FixIn%Bool%DO_SL_ISOTROPIC

!  LRRS Fixed Control inputs

      NSTREAMS          = LRRS_FixIn%Cont%NSTREAMS
      SOLAR_ANGLE       = LRRS_FixIn%Cont%SOLAR_ANGLE

!  LRRS Fixed User Value inputs

      N_USER_STREAMS    = LRRS_FixIn%UserVal%N_USER_STREAMS
      USER_ANGLES       = LRRS_FixIn%UserVal%USER_ANGLES

      N_USER_RELAZMS    = LRRS_FixIn%UserVal%N_USER_RELAZMS
      USER_RELAZMS      = LRRS_FixIn%UserVal%USER_RELAZMS

!  LRRS Modified Boolean inputs

      DO_USER_STREAMS   = LRRS_ModIn%MBool%DO_USER_STREAMS

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Initialize output status

      STATUS_INPUTCHECK = LRRS_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES   = 0
      MESSAGES(0) = 'Successful Check of SLEAVE/MAIN compatibility'
      ACTIONS(0)  = 'No Action required for this Task'

      NM = NMESSAGES

!  Checks

      IF ( SL_DO_SLEAVING .neqv. DO_SURFACE_LEAVING ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving control flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( SL_DO_ISOTROPIC .neqv. DO_SL_ISOTROPIC ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving isotropic flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( SL_DO_USER_STREAMS .neqv. DO_USER_STREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( SL_DO_USER_STREAMS .neqv. (.not.DO_MVOUT_ONLY) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams and DO_MVOUT_ONLY not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( SL_NSTREAMS .ne. NSTREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of discrete ordinates does not agree'
        ACTIONS(NM)  = 'Check SL_NSTREAMS and NSTREAMS input'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

!  Angles

      IF ( SL_BEAM_SZA .ne. SOLAR_ANGLE ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Solar beam angle does not agree'
        ACTIONS(NM)  = 'Check SL_BEAM_SZA and SOLAR_ANGLE input'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( SL_N_USER_STREAMS .ne. N_USER_STREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check SL_N_USER_STREAMS and N_USER_STREAMS input'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ELSE
        DO I = 1, N_USER_STREAMS
          if ( SL_USER_ANGLES_INPUT(I) .ne. USER_ANGLES(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'View zenith angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check SL_USER_ANGLES_INPUT & USER_ANGLES input'
            STATUS_INPUTCHECK = LRRS_SERIOUS
          endif
        ENDDO
      ENDIF

      IF ( SL_N_USER_RELAZMS .ne. N_USER_RELAZMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check SL_N_USER_RELAZMS & N_USER_RELAZMS input'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ELSE
        DO I = 1, N_USER_RELAZMS
          if ( SL_USER_RELAZMS(I) .ne. USER_RELAZMS(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Azimuth angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check SL_USER_RELAZMS & USER_RELAZMS input'
            STATUS_INPUTCHECK = LRRS_SERIOUS
          endif
        ENDDO
      ENDIF

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      LRRS_SLEAVECheck_Status%STATUS_INPUTCHECK = STATUS_INPUTCHECK
      LRRS_SLEAVECheck_Status%NCHECKMESSAGES    = NMESSAGES
      LRRS_SLEAVECheck_Status%CHECKMESSAGES     = MESSAGES
      LRRS_SLEAVECheck_Status%ACTIONS           = ACTIONS

!  Finish

      END SUBROUTINE LRRS_SLEAVE_INPUT_CHECKER

!

      SUBROUTINE BRDF_SLEAVE_INPUT_CHECKER ( &
        SLEAVE_Sup_In,             & ! Inputs
        BRDF_Sup_In,               & ! Inputs
        BRDF_SLEAVECheck_Status )    ! Outputs

!  This subroutine checks certain LRRS SLEAVE supplement inputs with the
!  corresponding LRRS BRDF supplement inputs for consistency 


!  Modules for use

      USE LRRS_PARS_m
      USE SLEAVE_Sup_Inputs_def_m
      USE BRDF_Sup_Inputs_def_m
      USE LRRS_Outputs_def_m

      IMPLICIT NONE

      TYPE(SLEAVE_Sup_Inputs), INTENT(IN)         :: SLEAVE_Sup_In
      TYPE(BRDF_Sup_Inputs), INTENT(IN)           :: BRDF_Sup_In

      TYPE(LRRS_Exception_Handling), INTENT(OUT)  :: BRDF_SLEAVECheck_Status

!  ---------------
!  Local variables
!  ---------------

!  SLEAVE supplement inputs
!  -------------------------

!  Surface-leaving control flags #1 (used for gatekeeping checks)

      LOGICAL :: SL_DO_ISOTROPIC
      LOGICAL :: SL_DO_FLUORESCENCE

!  Surface-leaving control flags #2

      LOGICAL :: SL_DO_GlintShadow
      LOGICAL :: SL_DO_FoamOption
      LOGICAL :: SL_DO_FacetIsotropy

!  Salinity

      DOUBLE PRECISION :: SL_SALINITY

!  Wind-speed and directions

      DOUBLE PRECISION :: SL_WINDSPEED
      DOUBLE PRECISION :: SL_WINDDIR

!  BRDF supplement inputs
!  -----------------------

!  Surface-leaving control flags

      LOGICAL :: BS_DO_GlintShadow
      LOGICAL :: BS_DO_FoamOption
      LOGICAL :: BS_DO_FacetIsotropy

!  Salinity

      DOUBLE PRECISION :: BS_SALINITY

!  Wind-speed and directions

      DOUBLE PRECISION :: BS_WINDSPEED
      DOUBLE PRECISION :: BS_WINDDIR

!  Exception handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS  ( 0:MAX_MESSAGES )

!  Other

      INTEGER          :: NM, I
      CHARACTER(Len=2) :: C2

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  SLEAVE Control inputs #1 (used for gatekeeping checks)

      SL_DO_ISOTROPIC     = SLEAVE_Sup_In%SL_DO_ISOTROPIC
      SL_DO_FLUORESCENCE  = SLEAVE_Sup_In%SL_DO_FLUORESCENCE

!  SLEAVE Control inputs #2

      SL_DO_GlintShadow   = SLEAVE_Sup_In%SL_DO_GlintShadow
      SL_DO_FoamOption    = SLEAVE_Sup_In%SL_DO_FoamOption
      SL_DO_FacetIsotropy = SLEAVE_Sup_In%SL_DO_FacetIsotropy

!  SLEAVE Other inputs

      SL_SALINITY         = SLEAVE_Sup_In%SL_SALINITY
      SL_WINDSPEED        = SLEAVE_Sup_In%SL_WINDSPEED
      SL_WINDDIR          = SLEAVE_Sup_In%SL_WINDDIR

!  BRDF Control inputs

      BS_DO_GlintShadow   = BRDF_Sup_In%BS_DO_GlintShadow
      BS_DO_FoamOption    = BRDF_Sup_In%BS_DO_FoamOption
      BS_DO_FacetIsotropy = BRDF_Sup_In%BS_DO_FacetIsotropy

!  BRDF Other inputs

      BS_SALINITY         = BRDF_Sup_In%BS_SALINITY
      BS_WINDSPEED        = BRDF_Sup_In%BS_WINDSPEED
      BS_WINDDIR          = BRDF_Sup_In%BS_WINDDIR

! Debug
!write(*,*)
!write(*,*) 'SL_DO_GlintShadow   = ',SL_DO_GlintShadow
!write(*,*) 'SL_DO_FoamOption    = ',SL_DO_FoamOption
!write(*,*) 'SL_DO_FacetIsotropy = ',SL_DO_FacetIsotropy
!write(*,*) 'SL_SALINITY  = ',SL_SALINITY
!write(*,*) 'SL_WINDSPEED = ',SL_WINDSPEED
!write(*,*) 'SL_WINDDIR   = ',SL_WINDDIR
!write(*,*)
!write(*,*) 'BS_DO_GlintShadow   = ',BS_DO_GlintShadow
!write(*,*) 'BS_DO_FoamOption    = ',BS_DO_FoamOption
!write(*,*) 'BS_DO_FacetIsotropy = ',BS_DO_FacetIsotropy
!write(*,*) 'BS_SALINITY  = ',BS_SALINITY
!write(*,*) 'BS_WINDSPEED = ',BS_WINDSPEED
!write(*,*) 'BS_WINDDIR   = ',BS_WINDDIR

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Initialize output status

      STATUS_INPUTCHECK = LRRS_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES   = 0
      MESSAGES(0) = 'Successful Check of SLEAVE/BRDF compatibility'
      ACTIONS(0)  = 'No Action required for this Task'

      NM = NMESSAGES

!  Checks
!  ------

!  NewCMGLINT/SLEAVE gatekeeper checks
!  (if any fails, don't look at inputs further until correct)

      IF ( SL_DO_ISOTROPIC ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'New Cox-Munk glint BRDF is active and surface-leaving DO_ISOTROPIC is active'
        ACTIONS(NM)  = 'Deactivate surface-leaving isotropy!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( SL_DO_FLUORESCENCE ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'New Cox-Munk glint BRDF is active and surface-leaving DO_FLUORESCENCE is active'
        ACTIONS(NM)  = 'Deactivate surface-leaving fluorescence!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( STATUS_INPUTCHECK == LRRS_SERIOUS ) GOTO 500

!  Control flags

      IF ( SL_DO_GlintShadow.neqv.BS_DO_GlintShadow ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving glint shadow flag does not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( SL_DO_FoamOption.neqv.BS_DO_FoamOption ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving foam option flag does not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( SL_DO_FacetIsotropy.neqv.BS_DO_FacetIsotropy ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving facet isotropy flag does not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

!  Salinity

      IF ( SL_SALINITY .ne. BS_SALINITY) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving salinity does not agree'
        ACTIONS(NM)  = 'Check SL_SALINITY and BS_SALINITY input'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

!  Wind-speed and directions

      IF ( SL_WINDSPEED .ne. BS_WINDSPEED) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving wind-speed does not agree'
        ACTIONS(NM)  = 'Check SL_WINDSPEED and BS_WINDSPEED input'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

      IF ( SL_WINDDIR .ne. BS_WINDDIR ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Wind direction angle does not agree'
        ACTIONS(NM)  = 'Check SL_WINDDIR and BS_WINDDIR input'
        STATUS_INPUTCHECK = LRRS_SERIOUS
      ENDIF

500   CONTINUE

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      BRDF_SLEAVECheck_Status%STATUS_INPUTCHECK = STATUS_INPUTCHECK
      BRDF_SLEAVECheck_Status%NCHECKMESSAGES    = NMESSAGES
      BRDF_SLEAVECheck_Status%CHECKMESSAGES     = MESSAGES
      BRDF_SLEAVECheck_Status%ACTIONS           = ACTIONS

!  Finish

      END SUBROUTINE BRDF_SLEAVE_INPUT_CHECKER

!  Finish Module

      END MODULE lrrs_sup_accessories_1_m
