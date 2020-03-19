
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
! #            SET_LRRS_BRDF_INPUTS                             #
! #            SET_LRRS_SLEAVE_INPUTS                           #
! #                                                             #
! ###############################################################

!  Upgrade to LRRS Version 2.5
!  ---------------------------

!    ** Routines to define the main LRRS BRDF/SLEAVE inputs either by
!       just initializing them or by defining them using the corresponding
!       LRRS BRDF/SLEAVE supplement outputs 

      MODULE lrrs_sup_accessories_2_m

      PRIVATE
      PUBLIC :: SET_LRRS_BRDF_INPUTS,      &
                SET_LRRS_SLEAVE_INPUTS

      CONTAINS


      SUBROUTINE SET_LRRS_BRDF_INPUTS ( &
        NPTS_BRDF, NPTS_LRRS,               & !Inputs
        BRDF_Sup_Out, BRDF_LinSup_Out,      & !Inputs
        LRRS_FixIn, LRRS_ModIn, LRRS_LinIn, & !Inputs
        LRRS_Sup, LRRS_LinSup,              & !Outputs
        FAIL, N_MESSAGES, MESSAGES )          !Outputs

!  This subroutine defines the main LRRS BRDF inputs either by just initializing them
!  or by defining them using the corresponding LRRS BRDF supplement outputs 

!  Use Modules

      USE BRDF_Sup_Outputs_def_m
      USE BRDF_LinSup_Outputs_def_m

      USE LRRS_PARS_m

      USE LRRS_IO_DEFS_m
      USE LRRS_LIN_IO_DEFS_m

      USE LRRS_Sup_InOut_def_m
      USE LRRS_Lin_Sup_InOut_def_m

!  Inputs

      INTEGER, INTENT(IN)                    :: NPTS_BRDF, NPTS_LRRS
      TYPE(BRDF_Sup_Outputs), INTENT(IN)     :: BRDF_Sup_Out
      TYPE(BRDF_LinSup_Outputs), INTENT(IN)  :: BRDF_LinSup_Out
      TYPE(LRRS_Fixed_Inputs), INTENT(IN)    :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT(IN) :: LRRS_ModIn
      TYPE(LRRS_LinInputs), INTENT(IN)       :: LRRS_LinIn

!  Outputs

      TYPE(LRRS_Sup_InOut), INTENT(INOUT)    :: LRRS_Sup
      TYPE(LRRS_LinSup_InOut), INTENT(INOUT) :: LRRS_LinSup

!  Error output

      LOGICAL, INTENT(OUT)                  :: FAIL
      INTEGER, INTENT(INOUT)                :: N_MESSAGES
      CHARACTER (LEN=*), INTENT(INOUT)      :: MESSAGES ( MAX_MESSAGES )

!  ---------------
!  Local variables
!  ---------------

      INTEGER :: J, PT
      INTEGER :: NUSERS, NAZIMS, NDISOS, NMOMS, NWFS
      LOGICAL :: DO_USER_STREAMS

!  Start program

!  Check some inputs

      FAIL = .FALSE.
      IF ( (NPTS_BRDF /= 1) .AND. (NPTS_BRDF /= NPTS_LRRS) ) THEN
        FAIL = .TRUE.
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = 'Error: # of BRDF spectral points should be = 1 or ' // &
                               ' # of LRRS spectral points'
        RETURN
      ENDIF

!  Define some local variables

      DO_USER_STREAMS = LRRS_ModIn%MBool%DO_USER_STREAMS

      NDISOS = LRRS_FixIn%Cont%NSTREAMS
      NMOMS  = 2*NDISOS - 1

      IF ( DO_USER_STREAMS ) THEN
        NUSERS = LRRS_FixIn%UserVal%N_USER_STREAMS
        NAZIMS = LRRS_FixIn%UserVal%N_USER_RELAZMS
      ENDIF

!  Initialize LRRS BRDF inputs

      LRRS_Sup%BRDF%EXACTDB_BRDFUNC = ZERO
      LRRS_Sup%BRDF%BRDF_F_0        = ZERO
      LRRS_Sup%BRDF%BRDF_F          = ZERO
      LRRS_Sup%BRDF%USER_BRDF_F_0   = ZERO
      LRRS_Sup%BRDF%USER_BRDF_F     = ZERO
 
      LRRS_LinSup%BRDF%LS_EXACTDB_BRDFUNC = ZERO
      LRRS_LinSup%BRDF%LS_BRDF_F_0        = ZERO
      LRRS_LinSup%BRDF%LS_BRDF_F          = ZERO
      LRRS_LinSup%BRDF%LS_USER_BRDF_F_0   = ZERO
      LRRS_LinSup%BRDF%LS_USER_BRDF_F     = ZERO

!  Set LRRS standard BRDF inputs

      IF (NPTS_BRDF == 1) THEN
        LRRS_Sup%BRDF%DO_BRDF_Wav1 = .TRUE.
      ELSE
        LRRS_Sup%BRDF%DO_BRDF_Wav1 = .FALSE.
      ENDIF

      DO PT=1,NPTS_LRRS
        IF (NPTS_BRDF == 1) THEN
          J = 1
        ELSE
          J = PT
        ENDIF

        LRRS_Sup%BRDF%BRDF_F_0     (0:NMOMS,1:NDISOS,PT)          = &
          BRDF_Sup_Out%BS_BRDF_F_0 (0:NMOMS,1:NDISOS,J)
        LRRS_Sup%BRDF%BRDF_F       (0:NMOMS,1:NDISOS,1:NDISOS,PT) = &
          BRDF_Sup_Out%BS_BRDF_F   (0:NMOMS,1:NDISOS,1:NDISOS,J)

        IF ( DO_USER_STREAMS ) THEN
          LRRS_Sup%BRDF%EXACTDB_BRDFUNC     (1:NUSERS,1:NAZIMS,PT)         = &
            BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC (1:NUSERS,1:NAZIMS,J)
          LRRS_Sup%BRDF%USER_BRDF_F_0       (0:NMOMS,1:NUSERS,PT)          = &
            BRDF_Sup_Out%BS_USER_BRDF_F_0   (0:NMOMS,1:NUSERS,J)
          LRRS_Sup%BRDF%USER_BRDF_F         (0:NMOMS,1:NUSERS,1:NDISOS,PT) = &
            BRDF_Sup_Out%BS_USER_BRDF_F     (0:NMOMS,1:NUSERS,1:NDISOS,J)
        ENDIF
      ENDDO

!  Set LRRS linearized BRDF inputs

      IF ( LRRS_LinIn%Cont%DO_SURFACE_WFS ) THEN
        NWFS = LRRS_LinIn%Cont%N_SURFACE_WFS

        DO PT=1,NPTS_LRRS
          IF (NPTS_BRDF == 1) THEN
            J = 1
          ELSE
            J = PT
          ENDIF

          LRRS_LinSup%BRDF%LS_BRDF_F_0     (1:NWFS,0:NMOMS,1:NDISOS,PT)          = &
            BRDF_LinSup_Out%BS_LS_BRDF_F_0 (1:NWFS,0:NMOMS,1:NDISOS,J)
          LRRS_LinSup%BRDF%LS_BRDF_F       (1:NWFS,0:NMOMS,1:NDISOS,1:NDISOS,PT) = &
            BRDF_LinSup_Out%BS_LS_BRDF_F   (1:NWFS,0:NMOMS,1:NDISOS,1:NDISOS,J)

          IF ( DO_USER_STREAMS ) THEN
            LRRS_LinSup%BRDF%LS_EXACTDB_BRDFUNC     (1:NWFS,1:NUSERS,1:NAZIMS,PT)         = &
              BRDF_LinSup_Out%BS_LS_DBOUNCE_BRDFUNC (1:NWFS,1:NUSERS,1:NAZIMS,J)
            LRRS_LinSup%BRDF%LS_USER_BRDF_F_0       (1:NWFS,0:NMOMS,1:NUSERS,PT)          = &
              BRDF_LinSup_Out%BS_LS_USER_BRDF_F_0   (1:NWFS,0:NMOMS,1:NUSERS,J)
            LRRS_LinSup%BRDF%LS_USER_BRDF_F         (1:NWFS,0:NMOMS,1:NUSERS,1:NDISOS,PT) = &
              BRDF_LinSup_Out%BS_LS_USER_BRDF_F     (1:NWFS,0:NMOMS,1:NUSERS,1:NDISOS,J)
          ENDIF
        ENDDO
      ENDIF

      END SUBROUTINE SET_LRRS_BRDF_INPUTS

!

      SUBROUTINE SET_LRRS_SLEAVE_INPUTS ( &
        NPTS_SLV, NPTS_LRRS,                & !Inputs
        SLEAVE_Sup_Out, SLEAVE_LinSup_Out,  & !Inputs
        LRRS_FixIn, LRRS_ModIn, LRRS_LinIn, & !Inputs
        LRRS_Sup, LRRS_LinSup,              & !Outputs
        FAIL, N_MESSAGES, MESSAGES )          !Outputs

!  This subroutine defines the main LRRS SLEAVE inputs either by just initializing them
!  or by defining them using the corresponding LRRS SLEAVE supplement outputs 

!  Use Modules

      USE SLEAVE_Sup_Outputs_def_m
      USE SLEAVE_LinSup_Outputs_def_m

      USE LRRS_PARS_m

      USE LRRS_IO_DEFS_m
      USE LRRS_LIN_IO_DEFS_m

      USE LRRS_Sup_InOut_def_m
      USE LRRS_Lin_Sup_InOut_def_m

!  Inputs

      INTEGER, INTENT(IN)                     :: NPTS_SLV, NPTS_LRRS
      TYPE(SLEAVE_Sup_Outputs), INTENT(IN)    :: SLEAVE_Sup_Out
      TYPE(SLEAVE_LinSup_Outputs), INTENT(IN) :: SLEAVE_LinSup_Out
      TYPE(LRRS_Fixed_Inputs), INTENT(IN)     :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT(IN)  :: LRRS_ModIn
      TYPE(LRRS_LinInputs), INTENT(IN)        :: LRRS_LinIn

!  Outputs

      TYPE(LRRS_Sup_InOut), INTENT(INOUT)     :: LRRS_Sup
      TYPE(LRRS_LinSup_InOut), INTENT(INOUT)  :: LRRS_LinSup

!  Error output

      LOGICAL, INTENT(OUT)                    :: FAIL
      INTEGER, INTENT(INOUT)                  :: N_MESSAGES
      CHARACTER (LEN=*), INTENT(INOUT)        :: MESSAGES ( MAX_MESSAGES )

!  ---------------
!  Local variables
!  ---------------

      INTEGER :: J, PT
      INTEGER :: NUSERS, NAZIMS, NDISOS, NMOMS, NWFS
      LOGICAL :: DO_USER_STREAMS

!  Start program

!  Check some inputs

      FAIL = .FALSE.
      IF ( (NPTS_SLV /= 1) .AND. (NPTS_SLV /= NPTS_LRRS) ) THEN
        FAIL = .TRUE.
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = 'Error: # of SLEAVE spectral points should be = 1 or ' // &
                               ' # of LRRS spectral points'
        RETURN
      ENDIF

!  Define some local variables

      DO_USER_STREAMS = LRRS_ModIn%MBool%DO_USER_STREAMS

      NDISOS = LRRS_FixIn%Cont%NSTREAMS
      NMOMS  = 2*NDISOS - 1

      IF ( DO_USER_STREAMS ) THEN
        NUSERS = LRRS_FixIn%UserVal%N_USER_STREAMS
        NAZIMS = LRRS_FixIn%UserVal%N_USER_RELAZMS
      ENDIF

!  Initialize LRRS SLEAVE inputs

      LRRS_Sup%SLEAVE%SLTERM_ISOTROPIC  = ZERO
      LRRS_Sup%SLEAVE%SLTERM_USERANGLES = ZERO
      LRRS_Sup%SLEAVE%SLTERM_F_0        = ZERO
      LRRS_Sup%SLEAVE%USER_SLTERM_F_0   = ZERO

      !Note: only isotropic used so far

      LRRS_LinSup%SLEAVE%LSSL_SLTERM_ISOTROPIC  = ZERO
      LRRS_LinSup%SLEAVE%LSSL_SLTERM_USERANGLES = ZERO
      LRRS_LinSup%SLEAVE%LSSL_SLTERM_F_0        = ZERO
      LRRS_LinSup%SLEAVE%LSSL_USER_SLTERM_F_0   = ZERO

!  Set LRRS standard SLEAVE inputs

      IF (NPTS_SLV == 1) THEN
        LRRS_Sup%SLEAVE%DO_SLEAVE_Wav1 = .TRUE.
      ELSE
        LRRS_Sup%SLEAVE%DO_SLEAVE_Wav1 = .FALSE.
      ENDIF

      DO PT=1,NPTS_LRRS
        IF (NPTS_SLV == 1) THEN
          J = 1
        ELSE
          J = PT
        ENDIF

        LRRS_Sup%SLEAVE%SLTERM_ISOTROPIC (PT) = &
          SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC (J)
        LRRS_Sup%SLEAVE%SLTERM_F_0     (0:NMOMS,1:NDISOS,PT) = &
          SLEAVE_Sup_Out%SL_SLTERM_F_0 (0:NMOMS,1:NDISOS,J)

        IF ( DO_USER_STREAMS ) THEN
          LRRS_Sup%SLEAVE%SLTERM_USERANGLES     (1:NUSERS,1:NAZIMS,PT) = &
            SLEAVE_Sup_Out%SL_SLTERM_USERANGLES (1:NUSERS,1:NAZIMS,J)
          LRRS_Sup%SLEAVE%USER_SLTERM_F_0       (0:NMOMS,1:NUSERS,PT)  = &
            SLEAVE_Sup_Out%SL_USER_SLTERM_F_0   (0:NMOMS,1:NUSERS,J)
        ENDIF
      ENDDO

!  Set LRRS linearized SLEAVE inputs

      IF ( LRRS_LinIn%Cont%DO_SLEAVE_WFS ) THEN
        NWFS = LRRS_LinIn%Cont%N_SLEAVE_WFS

        DO PT=1,NPTS_LRRS
          IF (NPTS_SLV == 1) THEN
            J = 1
          ELSE
            J = PT
          ENDIF

          LRRS_LinSup%SLEAVE%LSSL_SLTERM_ISOTROPIC   (1:NWFS,PT) = &
            SLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC (1:NWFS,J)
          LRRS_LinSup%SLEAVE%LSSL_SLTERM_F_0   (1:NWFS,0:NMOMS,1:NDISOS,PT) = &
            SLEAVE_LinSup_Out%SL_LS_SLTERM_F_0 (1:NWFS,0:NMOMS,1:NDISOS,J)

          IF ( DO_USER_STREAMS ) THEN
            LRRS_LinSup%SLEAVE%LSSL_SLTERM_USERANGLES   (1:NWFS,1:NUSERS,1:NAZIMS,PT) = &
              SLEAVE_LinSup_Out%SL_LS_SLTERM_USERANGLES (1:NWFS,1:NUSERS,1:NAZIMS,J)
            LRRS_LinSup%SLEAVE%LSSL_USER_SLTERM_F_0     (1:NWFS,0:NMOMS,1:NUSERS,PT)  = &
              SLEAVE_LinSup_Out%SL_LS_USER_SLTERM_F_0   (1:NWFS,0:NMOMS,1:NUSERS,J)
          ENDIF
        ENDDO
      ENDIF

      END SUBROUTINE SET_LRRS_SLEAVE_INPUTS

!  Finish Module

      END MODULE lrrs_sup_accessories_2_m
