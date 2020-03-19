
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
! #   INPUT SUBROUTINES :                                       #
! #                                                             #
! #     Initializing inputs:                                    #
! #                                                             #
! #          1. LRRS_L_INIT_INPUTS                              #
! #                                                             #
! #          2. LRRS_LinSup_Init                                #
! #          3. LRRS_BRDF_LinSup_Init                           #
! #          4. LRRS_SLEAVE_LinSup_Init                         #
! #          5. LRRS_SS_LinSup_Init                             #
! #                                                             #
! #     These routines are called at start of the Main          #
! #       linearized LRRS module                                #
! #                                                             #
! #          6. LRRS_L_CHECK_INPUT_DIMS                         #
! #                                                             #
! ###############################################################

!mick mod 7/20/2016

!  This is LRRS Version 2.5. This module added to accompany this version.
!    Main additions/changes:
!    (1) Addition of LRRS_L_INIT_INPUTS subroutine to separately handle the
!        initializing of LRRS linearized input variables
!    (2) Addition of the four subroutines to initialize the LRRS supplement
!        linearized input variables

      MODULE lrrs_L_io_check_m

!      USE LRRS_PARS_m, Only : LDU

      PRIVATE
      PUBLIC :: LRRS_L_INIT_INPUTS,      & 
                LRRS_LinSup_Init,        &
                LRRS_BRDF_LinSup_Init,   &
                LRRS_SLEAVE_LinSup_Init, &
                LRRS_SS_LinSup_Init,     &
                LRRS_L_CHECK_INPUT_DIMS

      CONTAINS

      SUBROUTINE LRRS_L_INIT_INPUTS ( LRRS_LinIn )

!  Module of dimensions and numbers

      USE LRRS_PARS_m
      USE LRRS_LIN_IO_DEFS_m

!  Implicit none

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT-RRS linearized input structures

      TYPE(LRRS_LinInputs), INTENT(OUT) :: LRRS_LinIn

!  Initialize inputs
!  =================

!  LRRS Linearized Control
!  -----------------------

      LRRS_LinIn%Cont%DO_PROFILE_WFS     = .FALSE.
      LRRS_LinIn%Cont%DO_COLUMN_WFS      = .FALSE.
      LRRS_LinIn%Cont%DO_SURFACE_WFS     = .FALSE.
      LRRS_LinIn%Cont%DO_SLEAVE_WFS      = .FALSE.

      LRRS_LinIn%Cont%DO_NORMALIZED_WFS  = .FALSE.
      LRRS_LinIn%Cont%DO_AIRPROFILE_WFS  = .FALSE.
      LRRS_LinIn%Cont%DO_TEMPPROFILE_WFS = .FALSE.
      LRRS_LinIn%Cont%DO_TEMPSHIFT_WF    = .FALSE.

      LRRS_LinIn%Cont%N_TOTALCOLUMN_WFS  = 0

      LRRS_LinIn%Cont%N_SURFACE_WFS      = 0
      LRRS_LinIn%Cont%N_SLEAVE_WFS       = 0

      LRRS_LinIn%Cont%LAYER_VARY_FLAG    = .FALSE.
      LRRS_LinIn%Cont%LAYER_VARY_NUMBER  = 0

!  LRRS Linearized Optical
!  -----------------------

      LRRS_LinIn%Optical%LAYER_AIRCOLUMNS_DT          = ZERO
      LRRS_LinIn%Optical%TEMPERATURES_UNSHIFTED       = ZERO

      LRRS_LinIn%Optical%L_DELTAU_INPUT_UNSCALED      = ZERO
      LRRS_LinIn%Optical%L_OMEGAMOMS_ELASTIC_UNSCALED = ZERO

!  Finish

      END SUBROUTINE LRRS_L_INIT_INPUTS

!

      SUBROUTINE LRRS_LinSup_Init ( LRRS_LinSup )

      USE LRRS_Lin_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LRRS linearized supplement input structure

      TYPE(LRRS_LinSup_InOut), INTENT(INOUT) :: LRRS_LinSup

!  Initialize LRRS linearized supplement inputs
!  ============================================

      CALL LRRS_BRDF_LinSup_Init   ( LRRS_LinSup )
      CALL LRRS_SLEAVE_LinSup_Init ( LRRS_LinSup )
      CALL LRRS_SS_LinSup_Init     ( LRRS_LinSup )

!  Finish

      END SUBROUTINE LRRS_LinSup_Init

!

      SUBROUTINE LRRS_BRDF_LinSup_Init ( LRRS_LinSup )

      USE LRRS_PARS_m, Only : ZERO
      USE LRRS_Lin_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LRRS linearized supplement input structure

      TYPE(LRRS_LinSup_InOut), INTENT(INOUT) :: LRRS_LinSup

!  Initialize LRRS linearized brdf supplement inputs
!  =================================================

      LRRS_LinSup%BRDF%LS_EXACTDB_BRDFUNC = ZERO
      LRRS_LinSup%BRDF%LS_BRDF_F_0        = ZERO
      LRRS_LinSup%BRDF%LS_BRDF_F          = ZERO
      LRRS_LinSup%BRDF%LS_USER_BRDF_F_0   = ZERO
      LRRS_LinSup%BRDF%LS_USER_BRDF_F     = ZERO

!  Finish

      END SUBROUTINE LRRS_BRDF_LinSup_Init

!

      SUBROUTINE LRRS_SLEAVE_LinSup_Init ( LRRS_LinSup )

      USE LRRS_PARS_m, Only : ZERO
      USE LRRS_Lin_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LRRS linearized supplement input structure

      TYPE(LRRS_LinSup_InOut), INTENT(INOUT) :: LRRS_LinSup

!  Initialize LRRS linearized sleave supplement inputs
!  ===================================================

      LRRS_LinSup%SLEAVE%LSSL_SLTERM_ISOTROPIC  = ZERO
      LRRS_LinSup%SLEAVE%LSSL_SLTERM_USERANGLES = ZERO
      LRRS_LinSup%SLEAVE%LSSL_SLTERM_F_0        = ZERO
      LRRS_LinSup%SLEAVE%LSSL_USER_SLTERM_F_0   = ZERO

!  Finish

      END SUBROUTINE LRRS_SLEAVE_LinSup_Init

!

      SUBROUTINE LRRS_SS_LinSup_Init ( LRRS_LinSup )

      USE LRRS_PARS_m, Only : ZERO
      USE LRRS_Lin_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LRRS linearized supplement input structure

      TYPE(LRRS_LinSup_InOut), INTENT(INOUT) :: LRRS_LinSup

!  Initialize LRRS linearized single-scatter supplement inputs
!  ===========================================================

      LRRS_LinSup%SS%Atmos%LC_ELASTIC_SS_UP = ZERO
      LRRS_LinSup%SS%Atmos%LC_ELASTIC_SS_DN = ZERO
      LRRS_LinSup%SS%Atmos%LC_ELASTIC_DB    = ZERO
      LRRS_LinSup%SS%Atmos%LP_ELASTIC_SS_UP = ZERO
      LRRS_LinSup%SS%Atmos%LP_ELASTIC_SS_DN = ZERO
      LRRS_LinSup%SS%Atmos%LP_ELASTIC_DB    = ZERO

      LRRS_LinSup%SS%Atmos%LC_RAMAN_SS_UP   = ZERO
      LRRS_LinSup%SS%Atmos%LC_RAMAN_SS_DN   = ZERO
      LRRS_LinSup%SS%Atmos%LC_RAMAN_DB      = ZERO
      LRRS_LinSup%SS%Atmos%LP_RAMAN_SS_UP   = ZERO
      LRRS_LinSup%SS%Atmos%LP_RAMAN_SS_DN   = ZERO
      LRRS_LinSup%SS%Atmos%LP_RAMAN_DB      = ZERO

      LRRS_LinSup%SS%Surf%LS_ELASTIC_DB     = ZERO
      LRRS_LinSup%SS%Surf%LS_RAMAN_DB       = ZERO

!  Finish

      END SUBROUTINE LRRS_SS_LinSup_Init

!

      SUBROUTINE LRRS_L_CHECK_INPUT_DIMS &
      ( LRRS_LinIn, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

!  module, dimensions and numbers

      USE LRRS_pars_m, only : MAX_ATMOSWFS, MAX_SURFACEWFS, &
                              MAX_SLEAVEWFS, MAX_MESSAGES, &
                              LRRS_SUCCESS, LRRS_SERIOUS

      USE LRRS_LinInputs_def_m

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  LRRS input structures

      TYPE(LRRS_LinInputs), INTENT (IN) :: LRRS_LinIn

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = LRRS_SUCCESS
      NM = NMESSAGES

!  Check LRRS linearized input dimensions
!    against maximum dimensions
!  =========================================

      IF ( LRRS_LinIn%Cont%DO_COLUMN_WFS ) THEN
        IF ( LRRS_LinIn%Cont%N_TOTALCOLUMN_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for column WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_TOTALCOLUMN_WFS or increase MAX_ATMOSWFS in LRRS.PARS'
         STATUS = LRRS_SERIOUS
        ENDIF
      ENDIF

!      IF ( LRRS_LinIn%Cont%DO_PROFILE_WFS ) THEN
!        IF ( LRRS_LinIn%Cont%N_TOTALPROFILE_WFS .GT. MAX_ATMOSWFS ) THEN
!         NM = NM + 1
!         MESSAGES(NM) = &
!             'Bad error: Insuffient dimensioning for profile WFs'
!         ACTIONS(NM)  = &
!             'Action: Decrease N_TOTALPROFILE_WFS or increase MAX_ATMOSWFS in LRRS.PARS'
!         STATUS = LRRS_SERIOUS
!        ENDIF
!      ENDIF

      IF ( LRRS_LinIn%Cont%DO_SURFACE_WFS ) THEN
        IF ( LRRS_LinIn%Cont%N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for Surface WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_SURFACE_WFS or increase MAX_SURFACE_WFS in LRRS.PARS'
         STATUS = LRRS_SERIOUS
        ENDIF
      ENDIF

      IF ( LRRS_LinIn%Cont%DO_SLEAVE_WFS ) THEN
        IF ( LRRS_LinIn%Cont%N_SLEAVE_WFS .GT. MAX_SLEAVEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = &
             'Bad error: Insuffient dimensioning for water surface-leaving (SLEAVE) WFs'
         ACTIONS(NM)  = &
             'Action: Decrease N_SLEAVE_WFS or increase MAX_SLEAVEWFS in LRRS.PARS'
         STATUS = LRRS_SERIOUS
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

      END SUBROUTINE LRRS_L_CHECK_INPUT_DIMS

      END MODULE lrrs_L_io_check_m
