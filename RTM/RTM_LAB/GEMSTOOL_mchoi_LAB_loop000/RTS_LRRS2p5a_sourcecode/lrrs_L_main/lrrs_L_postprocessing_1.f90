
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
! #   LINEARIZATION SUBROUTINES :                               #
! #       ( 1 = Whole-layer, 2 = partial-layer )                #
! #                                                             #
! #     LC  denotes column   linearization                      #
! #     LCH denotes column   linearization with heights         #
! #     LP  denotes profile  linearization                      #
! #     LS  denotes surface  linearization                      #
! #     LPC denotes LC or LP linearization                      #
! #                                                             #
! #   Top level PUBLIC routines                                 #
! #                                                             #
! #             LPC_POSTPROCESSING_1                            #
! #             LS_POSTPROCESSING_1                             #
! #                                                             #
! #   PUBLIC routines                                           #
! #                                                             #
! #             LPC_HOMOGMULT_1                                 #
! #             LPC_TOASOURCE, LPC_BOASOURCE                    #
! #             LS_TOASOURCE , LS_BOASOURCE                     #
! #                                                             #
! #   PRIVATE routines,  called by LPC_POSTPROCESSING_1         #
! #                                                             #
! #             LPC_BEAMMULT_UP_1                               #
! #             LPC_BEAMMULT_DN_1                               #
! #                                                             #
! #             LPC_RAMAN_SOURCETERM_UP_1                       #
! #             LPC_RAMAN_SOURCETERM_DN_1                       #
! #             LS_RAMAN_SOURCETERM_UP_1                        #
! #             LS_RAMAN_SOURCETERM_DN_1                        #
! #                                                             #
! ###############################################################


      MODULE lrrs_L_postprocessing_1_m

!      USE LRRS_PARS_m, Only : LDU

      PRIVATE
      PUBLIC :: LPC_POSTPROCESSING_1, &
                LS_POSTPROCESSING_1,  &
                LPC_HOMOGMULT_1,      &
                LPC_BEAMMULT_UP_1,    &
                LPC_BEAMMULT_DN_1,    &
                LPC_TOASOURCE,        &
                LS_TOASOURCE,         &
                LPC_BOASOURCE,        &
                LS_BOASOURCE

      CONTAINS

      SUBROUTINE LPC_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                                  & ! Inputs
            DO_MSMODE_LIDORT, DO_SSCORRECTION,                           & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                        & ! Inputs
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                      & ! Inputs
            KS, KF, KPARS, NKSTORAGE2, TAYLOR_ORDER, TAYLOR_SMALL,       & ! Inputs
            TRANS_USERM, DELTAUS, USER_STREAMS,                          & ! Inputs
            FLIPPER, CHANGE_EXPONENT, SOLUTION,                          & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                       & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1,        & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,            & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,                              & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,                              & ! Inputs
            L_TRANS_USERM, L_DELTAUS, L_ASOURCE, L_DELTRANS, L_INITRANS, & ! Inputs
            L_QSOURCE_UPOS1, L_QSOURCE_UNEG1, L_ATERM_SAVE,L_BTERM_SAVE, & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,          & ! Inputs
            L_WLAYER_PIST_UP, L_WLAYER_PIST_DN )                           ! Outputs

!  Atmospheric linearization, Source term post processing for one layer

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS, MAX_LAYERS_SQ, MAX_ATMOSWFS

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Up or down

      LOGICAL  , INTENT(IN) :: DO_UPWELLING
      LOGICAL  , INTENT(IN) :: DO_DNWELLING

!  Flags for multiple scattering and single scatter correction

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT
      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION

!  Fourier component

      INTEGER  , INTENT(IN) :: FOURIER

!  Number of layer

      INTEGER  , INTENT(IN) :: N

!  number of streams and solutions

      INTEGER  , INTENT(IN) :: NSTREAMS, N_USER_STREAMS

!  Linearization control numbers

      INTEGER  , INTENT(IN) :: KS, KF, KPARS(0:MAX_LAYERS)
      INTEGER  , INTENT(IN) :: NKSTORAGE2(MAX_LAYERS,0:MAX_LAYERS)

!  layer control

      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN

!  Taylor-series control. Updated, Version 2.5, 10/10/15

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  Control for the particular solutions

      LOGICAL  , INTENT(IN) :: CHANGE_EXPONENT
      LOGICAL  , INTENT(IN) :: FLIPPER
      INTEGER  , INTENT(IN) :: SOLUTION

!  Optical thickness value

      REAL(FPK), INTENT(IN) :: DELTAUS

!  stream directions

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  Layer transmittances in stream directions

      REAL(FPK), INTENT(IN) :: TRANS_USERM  ( MAX_USER_STREAMS )

!  Particular solution attributes

      INTEGER  , INTENT(IN) :: PICUTOFF
      REAL(FPK), INTENT(IN) :: ASOURCE
      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: INITRANS

!  Raman source term vectors

      REAL(FPK), INTENT(IN) :: QSOURCE_UPOS1 ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: QSOURCE_UNEG1 ( MAX_USER_STREAMS )

!  Saved ATERM and BTERM contributions to the Green's function.

      REAL(FPK), INTENT(IN) :: ATERM_SAVE  ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: BTERM_SAVE  ( MAX_STREAMS )

!  Eigensolutions at user-defined stream angles

      REAL(FPK), INTENT(IN) :: &
          U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Eigenvalues and transmittance factors.

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )

!  Whole layer multipliers, homogeneous solutions

      REAL(FPK), INTENT(IN) :: HMULT_1 &
             ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: HMULT_2 &
             ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Particular integral multipliers

      REAL(FPK), INTENT(IN) :: EMULT_UP      ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: EMULT_DN      ( MAX_USER_STREAMS )

!  Saved multiplier arrays

      REAL(FPK), INTENT(IN) :: WSU_SUSAV (MAX_STREAMS,MAX_USER_STREAMS)
      REAL(FPK), INTENT(IN) :: WSU_SDSAV (MAX_STREAMS,MAX_USER_STREAMS)
      REAL(FPK), INTENT(IN) :: WSD_SUSAV (MAX_STREAMS,MAX_USER_STREAMS)
      REAL(FPK), INTENT(IN) :: WSD_SDSAV (MAX_STREAMS,MAX_USER_STREAMS)

!  Linearized inputs

!  Linearizations

      REAL(FPK), INTENT(IN) :: L_DELTAUS     ( MAX_ATMOSWFS )
      REAL(FPK), INTENT(IN) :: L_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS )

      REAL(FPK), INTENT(IN) :: L_ASOURCE      ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_DELTRANS     ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_INITRANS     ( MAX_ATMOSWFS, 0:MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_ATERM_SAVE &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: L_BTERM_SAVE &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS )

      REAL(FPK), INTENT(IN) :: L_KEIGEN ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_HMULT_1 &
          ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_HMULT_2 &
          ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_U_XPOS &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_XNEG &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_QSOURCE_UPOS1 &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: L_QSOURCE_UNEG1 &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS )

!  output
!  ------

!  Layer source terms
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: L_WLAYER_PIST_UP &
          ( MAX_ATMOSWFS, MAX_LAYERS_SQ, MAX_USER_STREAMS )
      REAL(FPK), INTENT(INOUT) :: L_WLAYER_PIST_DN &
          ( MAX_ATMOSWFS, MAX_LAYERS_SQ, MAX_USER_STREAMS )

!  Local variables
!  ---------------

!  Linearized multipliers

      REAL(FPK) :: L_EMULT_UP &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS , MAX_USER_STREAMS )
      REAL(FPK) :: L_EMULT_DN &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS , MAX_USER_STREAMS )

!  Upwelling
!  ---------

      IF ( DO_UPWELLING ) THEN

        CALL LPC_BEAMMULT_UP_1 &
            ( N, N_USER_STREAMS, USER_STREAMS,               & ! Inputs
              TAYLOR_SMALL, TAYLOR_ORDER,                    & ! Inputs
              KS, KF, KPARS, STERM_LAYERMASK_UP,             & ! Inputs
              TRANS_USERM, L_TRANS_USERM,                    & ! Inputs
              DELTAUS, FLIPPER, PICUTOFF, DELTRANS, ASOURCE, & ! Inputs
              EMULT_UP, L_DELTAUS, L_DELTRANS, L_ASOURCE,    & ! Inputs
              L_EMULT_UP )                                     ! Outputs

        CALL LPC_RAMAN_SOURCETERM_UP_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION,                      & ! Inputs
            N, STERM_LAYERMASK_UP, NSTREAMS, N_USER_STREAMS,        & ! Inputs
            TAYLOR_SMALL, TAYLOR_ORDER, KS, KF, KPARS, NKSTORAGE2,  & ! Inputs
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION,  & ! Inputs
            TRANS_USERM, DELTAUS, USER_STREAMS,                     & ! Inputs
            EMULT_UP, ATERM_SAVE, BTERM_SAVE, WSU_SUSAV, WSU_SDSAV, & ! Inputs
            QSOURCE_UPOS1, ASOURCE, DELTRANS, INITRANS,             & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,       & ! Inputs
            L_TRANS_USERM, L_DELTAUS,                               & ! Inputs
            L_EMULT_UP, L_ATERM_SAVE, L_BTERM_SAVE,                 & ! Inputs
            L_QSOURCE_UPOS1, L_ASOURCE, L_DELTRANS, L_INITRANS,     & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,     & ! Inputs
            L_WLAYER_PIST_UP )                                        ! Outputs

      ENDIF

!  Downwelling
!  -----------

      IF ( DO_DNWELLING ) THEN

        CALL LPC_BEAMMULT_DN_1 &
            ( N, N_USER_STREAMS, USER_STREAMS,               & ! Inputs
              TAYLOR_SMALL, TAYLOR_ORDER,                    & ! Inputs
              KS, KF, KPARS, STERM_LAYERMASK_DN,             & ! Inputs
              TRANS_USERM, L_TRANS_USERM,                    & ! Inputs
              DELTAUS, FLIPPER, PICUTOFF, DELTRANS, ASOURCE, & ! Inputs
              EMULT_DN, L_DELTAUS, L_DELTRANS, L_ASOURCE,    & ! Inputs
              L_EMULT_DN )                                     ! Outputs

        CALL LPC_RAMAN_SOURCETERM_DN_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION,                      & ! Inputs
            N, STERM_LAYERMASK_DN, NSTREAMS, N_USER_STREAMS,        & ! Inputs
            TAYLOR_SMALL, TAYLOR_ORDER, KS, KF, KPARS, NKSTORAGE2,  & ! Inputs
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION,  & ! Inputs
            TRANS_USERM, DELTAUS, USER_STREAMS,                     & ! Inputs
            EMULT_DN, ATERM_SAVE, BTERM_SAVE, WSD_SUSAV, WSD_SDSAV, & ! Inputs
            QSOURCE_UNEG1, ASOURCE, DELTRANS, INITRANS,             & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,       & ! Inputs
            L_TRANS_USERM, L_DELTAUS,                               & ! Inputs
            L_EMULT_DN, L_ATERM_SAVE, L_BTERM_SAVE,                 & ! Inputs
            L_QSOURCE_UNEG1, L_ASOURCE, L_DELTRANS, L_INITRANS,     & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,     & ! Inputs
            L_WLAYER_PIST_DN )                                        ! Outputs

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LPC_POSTPROCESSING_1

!

      SUBROUTINE LS_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                 & ! Inputs
            DO_MSMODE_LIDORT, DO_SSCORRECTION, FOURIER, & ! Inputs
            N, NSTREAMS, N_USER_STREAMS, N_SURFACEWFS,  & ! Inputs
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,     & ! Inputs
            SOLUTION, INITRANS, U_XPOS, U_XNEG,         & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,             & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,             & ! Inputs
            LS_ATERM_SAVE, LS_BTERM_SAVE,               & ! Inputs
            LS_QSOURCE_UPOS1, LS_QSOURCE_UNEG1,         & ! Inputs
            LS_WLAYER_PIST_UP, LS_WLAYER_PIST_DN )        ! Outputs

!  Atmospheric linearization, Source term post processing for one layer

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS, MAX_SURFACEWFS

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Up or down

      LOGICAL  , INTENT(IN) :: DO_UPWELLING
      LOGICAL  , INTENT(IN) :: DO_DNWELLING

!  Flags for multiple scattering and single scatter correction

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT
      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION

!  Fourier component

      INTEGER  , INTENT(IN) :: FOURIER

!  Number of layer

      INTEGER  , INTENT(IN) :: N

!  number of streams and solutions

      INTEGER  , INTENT(IN) :: NSTREAMS, N_USER_STREAMS

!  number of Surface Jacobians

      INTEGER  , INTENT(IN) :: N_SURFACEWFS

!  layer control

      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN

!  Control for the particular solutions

      INTEGER  , INTENT(IN) :: SOLUTION

!  Particular solution attributes

      REAL(FPK), INTENT(IN) :: INITRANS

!  Eigensolutions at user-defined stream angles

      REAL(FPK), INTENT(IN) :: &
          U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Particular integral multipliers

      REAL(FPK), INTENT(IN) :: EMULT_UP      ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: EMULT_DN      ( MAX_USER_STREAMS )

!  Saved multiplier arrays

      REAL(FPK), INTENT(IN) :: WSU_SUSAV (MAX_STREAMS,MAX_USER_STREAMS)
      REAL(FPK), INTENT(IN) :: WSU_SDSAV (MAX_STREAMS,MAX_USER_STREAMS)
      REAL(FPK), INTENT(IN) :: WSD_SUSAV (MAX_STREAMS,MAX_USER_STREAMS)
      REAL(FPK), INTENT(IN) :: WSD_SDSAV (MAX_STREAMS,MAX_USER_STREAMS)

!  Linearized inputs
!  -----------------

!  Saved ATERM and BTERM contributions to the Green's function.

      REAL(FPK), INTENT(IN) :: LS_ATERM_SAVE  ( MAX_SURFACEWFS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: LS_BTERM_SAVE  ( MAX_SURFACEWFS, MAX_STREAMS )

!  Raman source vectors

      REAL(FPK), INTENT(IN) :: LS_QSOURCE_UPOS1 ( MAX_SURFACEWFS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: LS_QSOURCE_UNEG1 ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  output
!  ------

!  Linearized Layer source terms
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: LS_WLAYER_PIST_UP (MAX_SURFACEWFS,MAX_LAYERS,MAX_USER_STREAMS)
      REAL(FPK), INTENT(INOUT) :: LS_WLAYER_PIST_DN (MAX_SURFACEWFS,MAX_LAYERS,MAX_USER_STREAMS)

!  Upwelling

      IF ( DO_UPWELLING ) THEN
        CALL LS_RAMAN_SOURCETERM_UP_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION,               & ! Inputs
            N, STERM_LAYERMASK_UP, NSTREAMS, N_USER_STREAMS, & ! Inputs
            N_SURFACEWFS, FOURIER, SOLUTION, INITRANS,       & ! Inputs
            U_XPOS, U_XNEG, EMULT_UP, WSU_SUSAV, WSU_SDSAV,  & ! Inputs
            LS_ATERM_SAVE, LS_BTERM_SAVE, LS_QSOURCE_UPOS1,  & ! Inputs
            LS_WLAYER_PIST_UP )                                ! Output
      ENDIF

!  Downwelling

      IF ( DO_DNWELLING ) THEN
        CALl LS_RAMAN_SOURCETERM_DN_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION,               & ! Inputs
            N, STERM_LAYERMASK_DN, NSTREAMS, N_USER_STREAMS, & ! Inputs
            N_SURFACEWFS, FOURIER, SOLUTION, INITRANS,       & ! Inputs
            U_XPOS, U_XNEG, EMULT_DN, WSD_SUSAV, WSD_SDSAV,  & ! Inputs
            LS_ATERM_SAVE, LS_BTERM_SAVE, LS_QSOURCE_UNEG1,  & ! Inputs
            LS_WLAYER_PIST_DN )                                ! Output
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LS_POSTPROCESSING_1

!

      SUBROUTINE LPC_HOMOGMULT_1 &
          ( NSTREAMS, N, N_USER_STREAMS, DOVARY, NPARS,                 & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, USER_STREAMS, UTRANS, L_UTRANS, & ! Inputs
            KTRANS, L_KEIGEN, L_KTRANS, DELTAUS, L_DELTAUS,             & ! Inputs
            HMULT_1, HMULT_2, ZETA_P, ZETA_M,                           & ! Inputs
            L_HMULT_1, L_HMULT_2 )                                        ! Outputs

!  Whole-layer Linearized homogeneous solution MULTIPLIERS
!      upwelling and downwelling for the complete atmosphere

!   @@@ RobFix 5/5/11. Small numbers analysis:
!                      added FGSMALL, DELTAUS, L_DELTAUS to Argument list.
!   @@@ RobFix 10/10/15. Small numbers analysis updated with Taylor series routines

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m    , Only : FPK, ZERO, ONE, MAX_USER_STREAMS, MAX_ATMOSWFS, &
                                MAX_STREAMS, MAX_2_STREAMS, MAX_LAYERS

      USE lrrs_Taylor_m, Only : Taylor_series_L_1

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  control inputs

      INTEGER  , INTENT(IN) :: N_USER_STREAMS, N, NSTREAMS, NPARS
      LOGICAL  , INTENT(IN) :: DOVARY

!  Taylor-series control. Updated, Version 2.5, 10/10/15

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  optical thickness for use with the small number expansion
!   @@@ RobFix 5/5/11. Small numbers analysis added

      REAL(FPK), INTENT(IN) :: DELTAUS, L_DELTAUS(MAX_ATMOSWFS)

!  User stream input

      REAL(FPK), INTENT(IN) :: UTRANS       ( MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  Eigenproblem input

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  Linearized inputs

      REAL(FPK), INTENT(IN) :: L_KEIGEN ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_UTRANS &
                  ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS )

!  Already calculated multipliers and zeta factors

      REAL(FPK), INTENT(IN) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: &
            ZETA_P ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            ZETA_M ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  output
!mick fix 10/19/2015 - changed intent to inout
      REAL(FPK), INTENT(INOUT) :: L_HMULT_1 &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(INOUT) :: L_HMULT_2 &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Local variables
!  ---------------

      INTEGER   :: UM, AA, Q
      REAL(FPK) :: UDEL, ZDEL, L_UDEL, L_ZDEL
      REAL(FPK) :: L_T2, L_T1, HOM1, HOM2, SM, RHO_M

!  Get the linearized multipliers

      IF ( DOVARY ) THEN
        DO UM = 1, N_USER_STREAMS
          UDEL = UTRANS(UM,N)
          SM  = ONE / USER_STREAMS(UM)
          DO AA = 1, NSTREAMS
            RHO_M = ONE / ZETA_M(AA,UM,N)
            ZDEL  = KTRANS(AA,N)
            DO Q = 1, NPARS
              L_ZDEL = L_KTRANS(Q,AA,N)
              L_UDEL = L_UTRANS(Q,UM,N)
              L_T2 = - ZDEL * L_UDEL - L_ZDEL * UDEL
              L_T1 = L_ZDEL - L_UDEL

              IF ( ABS(RHO_M) .LT. TAYLOR_SMALL ) THEN
                 CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, RHO_M, DELTAUS, L_DELTAUS(Q), L_KEIGEN(Q,AA,N), ZERO, UDEL, SM, HOM1 )
! (2.3 call)     CALL LIMIT_L_GCFUNC ( -RHO_M, DELTAUS, SM, UDEL, 0.0d0, L_KEIGEN(Q,AA,N), L_DELTAUS(Q), HOM1 )

                 L_HMULT_1(Q,AA,UM,N) = SM * HOM1
              ELSE
                 HOM1 =   L_KEIGEN(Q,AA,N)*HMULT_1(AA,UM,N) + SM*L_T1
                 L_HMULT_1(Q,AA,UM,N) = ZETA_M(AA,UM,N) * HOM1
              ENDIF

              HOM2 = - L_KEIGEN(Q,AA,N)*HMULT_2(AA,UM,N) + SM*L_T2
              L_HMULT_2(Q,AA,UM,N) = ZETA_P(AA,UM,N) * HOM2
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LPC_HOMOGMULT_1

!

      SUBROUTINE LPC_BEAMMULT_UP_1 &
            ( N, N_USER_STREAMS, USER_STREAMS, TAYLOR_SMALL,   & ! Inputs
              TAYLOR_ORDER, KS, KF, KPARS, STERM_LAYERMASK_UP, & ! Inputs
              TRANS_USERM, L_TRANS_USERM,                      & ! Inputs
              DELTAUS, FLIPPER, PICUTOFF, DELTRANS, ASOURCE,   & ! Inputs
              EMULT_UP, L_DELTAUS, L_DELTRANS, L_ASOURCE,      & ! Inputs
              L_EMULT_UP )                                       ! Outputs

!  Source term primary multiplier. Upwelling.
!   - Small numbers analysis added 7 December 2006
!   @@@ RobFix 5/5/11. Small numbers analysis added
!      USE LRRS_L_SMALLNOS. Replaced by Taylor series routine

!  Module of dimensions and numbers

      USE LRRS_PARS_m    , Only : FPK, ZERO, ONE, MAX_USER_STREAMS, MAX_LAYERS, MAX_ATMOSWFS

      USE lrrs_Taylor_m, Only : Taylor_series_L_1

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  control INTEGERs

      INTEGER  , INTENT(IN) :: N_USER_STREAMS, N

!  Taylor-series control. Updated, Version 2.5, 10/10/15

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  Linearization control numbers

      INTEGER  , INTENT(IN) :: KS, KF, KPARS(0:MAX_LAYERS)

!  User streams

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP

!  layer transmittances + linearization

      REAL(FPK), INTENT(IN) :: TRANS_USERM   ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: L_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  optical thickness for use with the small number expansion

      REAL(FPK), INTENT(IN) :: DELTAUS
      REAL(FPK), INTENT(IN) :: L_DELTAUS ( MAX_ATMOSWFS )

!  Solution control

      LOGICAL  , INTENT(IN) :: FLIPPER
      INTEGER  , INTENT(IN) :: PICUTOFF

!  Source control and lineaerizations thereof

      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: ASOURCE
      REAL(FPK), INTENT(IN) :: L_DELTRANS ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_ASOURCE  ( MAX_ATMOSWFS, 0:MAX_LAYERS )

!  Multipliers for one solution source

      REAL(FPK), INTENT(IN) :: EMULT_UP ( MAX_USER_STREAMS )

!  output = Linearized multipliers for one solution source
!  -------------------------------------------------

      REAL(FPK), INTENT(OUT) :: L_EMULT_UP &
               ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS )

!  local

      INTEGER   :: UM, Q, K
      REAL(FPK) :: SM, UDEL, SIGMA_M, SIGMA_P, SD, SU
      REAL(FPK) :: L_WDEL, L_UDEL, L_WUDEL, L_DELT
      REAL(FPK) :: L_SIGMA_M, L_SIGMA_P, L_SD, L_SU

!  Only if there is a source term for this layer

      IF ( STERM_LAYERMASK_UP ) THEN

!  Nothing to calculate, outside layer of extinction

        IF ( N .GT. PICUTOFF ) THEN
          DO UM = 1, N_USER_STREAMS
            DO K = KS, KF
              DO Q = 1, KPARS(K)
                L_EMULT_UP(Q,K,UM) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ELSE

!  Start loops

          DO UM = 1, N_USER_STREAMS
            SM   = ONE / USER_STREAMS(UM)
            UDEL = TRANS_USERM(UM)

!  Linearized multiplier
!      K = 0 (column WFs) or K = n (profile WFs) requires non-zero terms
!      K < N (profile WFs), certain cross-layer terms are zero

            IF ( FLIPPER ) THEN
              SIGMA_M = ASOURCE - SM
              SD = EMULT_UP(UM) / SM
              DO K = KS, KF
                DO Q = 1, KPARS(K)
                  L_SIGMA_M = L_ASOURCE(Q,K)
                  L_WDEL    = L_DELTRANS(Q,K)
                  IF (K.EQ.0 .OR. K.EQ.N) THEN
                    L_UDEL = L_TRANS_USERM(Q,UM)
                    L_DELT = l_deltaus(q)
                  ELSE
                    L_UDEL = ZERO
                    L_DELT = ZERO
                  ENDIF
                  IF( ABS(SIGMA_M).LT.TAYLOR_SMALL)THEN
                    !CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, SIGMA_M, DELTAUS, L_DELT, L_ASOURCE(Q,K), ZERO, UDEL, SM, L_SD )
                    CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, -SIGMA_M, DELTAUS, L_DELT, L_ASOURCE(Q,K), ZERO, UDEL, SM, L_SD )
! (2.3 Call)        call limit_l_gcfunc ( sigma_m, deltaus, sm, udel, l_asource(q,k), 0.0d0, l_delt, l_sd )
                  ELSE
                    L_SD = ((L_UDEL-L_WDEL) - SD*L_SIGMA_M) / SIGMA_M
                  ENDIF
                  L_EMULT_UP(Q,K,UM) = SM * L_SD
                ENDDO
              ENDDO
            ELSE
              SIGMA_P = ASOURCE + SM
              SU = EMULT_UP(UM) / SM
              DO K = KS, KF
                DO Q = 1, KPARS(K)
                  L_SIGMA_P = L_ASOURCE(Q,K)
                  IF(K.EQ.0 .OR. K.EQ.N) THEN
                    L_UDEL = L_TRANS_USERM(Q,UM)
                  ELSE
                    L_UDEL = ZERO
                  ENDIF
                  L_WDEL    = L_DELTRANS(Q,K)
                  L_WUDEL   = DELTRANS * L_UDEL + L_WDEL * UDEL
                  L_SU = - ( L_WUDEL + SU * L_SIGMA_P ) / SIGMA_P
                  L_EMULT_UP(Q,K,UM) = SM * L_SU
                ENDDO
              ENDDO
            ENDIF

!  End loops

          ENDDO
        ENDIF
!mick fix 7/20/2016 - added ELSE section with initialization
      ELSE
        DO UM = 1, N_USER_STREAMS
          DO K = KS, KF
            DO Q = 1, KPARS(K)
              L_EMULT_UP(Q,K,UM) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LPC_BEAMMULT_UP_1

!

      SUBROUTINE LPC_BEAMMULT_DN_1 &
            ( N, N_USER_STREAMS, USER_STREAMS, TAYLOR_SMALL,   & ! Inputs
              TAYLOR_ORDER, KS, KF, KPARS, STERM_LAYERMASK_DN, & ! Inputs
              TRANS_USERM, L_TRANS_USERM,                      & ! Inputs
              DELTAUS, FLIPPER, PICUTOFF, DELTRANS, ASOURCE,   & ! Inputs
              EMULT_DN, L_DELTAUS, L_DELTRANS, L_ASOURCE,      & ! Inputs
              L_EMULT_DN )                                       ! Outputs

!  Source term primary multiplier. Downwelling.
!   - Small numbers analysis added 7 December 2006
!   @@@ RobFix 5/5/11. Small numbers analysis revised
!      USE LRRS_L_SMALLNOS. Replaced now by Taylor series routine

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, MAX_USER_STREAMS, MAX_LAYERS, MAX_ATMOSWFS

      USE lrrs_Taylor_m, Only : Taylor_series_L_1

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  control INTEGER ::s

      INTEGER  , INTENT(IN) :: N_USER_STREAMS, N

!  Taylor-series control. Updated, Version 2.5, 10/10/15

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  Linearization control numbers

      INTEGER  , INTENT(IN) :: KS, KF, KPARS(0:MAX_LAYERS)

!  User streams

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN

!  layer transmittances + linearization

      REAL(FPK), INTENT(IN) :: TRANS_USERM   ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: L_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  optical thickness for use with the small number expansion

      REAL(FPK), INTENT(IN) :: DELTAUS
      REAL(FPK), INTENT(IN) :: L_DELTAUS ( MAX_ATMOSWFS )

!  Solution control

      LOGICAL  , INTENT(IN) :: FLIPPER
      INTEGER  , INTENT(IN) :: PICUTOFF

!  Source control and lineaerizations thereof

      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: ASOURCE
      REAL(FPK), INTENT(IN) :: L_DELTRANS ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_ASOURCE  ( MAX_ATMOSWFS, 0:MAX_LAYERS )

!  Multipliers for one solution source

      REAL(FPK), INTENT(IN) :: EMULT_DN ( MAX_USER_STREAMS )

!  output = Linearized multipliers for one solution source
!  -------------------------------------------------

      REAL(FPK), INTENT(OUT) :: L_EMULT_DN &
               ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS )

!  local

      INTEGER   :: UM, Q, K
      REAL(FPK) :: SM, UDEL, SIGMA_M, SIGMA_P, SD, SU
      REAL(FPK) :: L_WDEL, L_UDEL, L_WUDEL, L_DELT
      REAL(FPK) :: L_SIGMA_M, L_SIGMA_P, L_SD, L_SU

      IF ( STERM_LAYERMASK_DN ) THEN

!  Nothing to calculate outside layer of extinction

        IF ( N .GT. PICUTOFF ) THEN
          DO UM = 1, N_USER_STREAMS
            DO K = KS, KF
              DO Q = 1, KPARS(K)
                L_EMULT_DN(Q,K,UM) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ELSE

!  Start loops

          DO UM = 1, N_USER_STREAMS
            SM   = ONE / USER_STREAMS(UM)
            UDEL = TRANS_USERM(UM)

!  Linearized multiplier
!      K = 0 (column WFs) or K = n (profile WFs) requires non-zero terms
!      K < N (profile WFs), certain cross-layer terms are zero

            IF ( FLIPPER ) THEN
              SIGMA_P = ASOURCE + SM
              SU = EMULT_DN(UM) / SM
              DO K = KS, KF
                DO Q = 1, KPARS(K)
                  L_SIGMA_P = L_ASOURCE(Q,K)
                  L_WDEL    = L_DELTRANS(Q,K)
                  IF(K.EQ.0 .OR. K.EQ.N) THEN
                    L_UDEL = L_TRANS_USERM(Q,UM)
                  ELSE
                    L_UDEL = ZERO
                  ENDIF
                  L_WUDEL   = DELTRANS * L_UDEL + L_WDEL * UDEL
                  L_SU = - ( L_WUDEL + SU * L_SIGMA_P ) / SIGMA_P
                  L_EMULT_DN(Q,K,UM) = SM * L_SU
                ENDDO
              ENDDO
            ELSE
              SIGMA_M = ASOURCE - SM
              SD = EMULT_DN(UM) / SM
              DO K = KS, KF
                DO Q = 1, KPARS(K)
                  L_SIGMA_M = L_ASOURCE(Q,K)
                  L_WDEL    = L_DELTRANS(Q,K)
                  IF(K.EQ.0 .OR. K.EQ.N) THEN
                    L_UDEL = L_TRANS_USERM(Q,UM)
                    L_DELT = L_DELTAUS(Q)
                  ELSE
                    L_UDEL = ZERO
                    L_DELT = ZERO
                  ENDIF
                  IF( ABS(SIGMA_M).LT.TAYLOR_SMALL)THEN
                    !CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, SIGMA_M, DELTAUS, L_DELT, L_ASOURCE(Q,K), ZERO, UDEL, SM, L_SD )
                    CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, -SIGMA_M, DELTAUS, L_DELT, L_ASOURCE(Q,K), ZERO, UDEL, SM, L_SD )
!  (2.3 Call)       call limit_l_gcfunc ( sigma_m, deltaus, sm, udel, l_asource(q,k), 0.0d0, l_delt, l_sd )
                  ELSE
                    L_SD = ((L_UDEL-L_WDEL) - SD*L_SIGMA_M) / SIGMA_M
                  ENDIF
                  L_EMULT_DN(Q,K,UM) = SM * L_SD
                ENDDO
              ENDDO
            ENDIF

!  End loops

          ENDDO
        ENDIF
!mick fix 7/20/2016 - added ELSE section with initialization
      ELSE
        DO UM = 1, N_USER_STREAMS
          DO K = KS, KF
            DO Q = 1, KPARS(K)
              L_EMULT_DN(Q,K,UM) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LPC_BEAMMULT_DN_1

!

      SUBROUTINE LPC_TOASOURCE &
          ( N_USER_STREAMS, KPARS, &
            L_TOA_SOURCE )

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m,    Only : fpk, ZERO, MAX_ATMOSWFS,MAX_USER_STREAMS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER  , INTENT(IN)  :: N_USER_STREAMS, KPARS
      REAL(FPK), INTENT(OUT) :: L_TOA_SOURCE(MAX_ATMOSWFS,MAX_USER_STREAMS)

!  local variables

      INTEGER :: UM, Q

!  initialise TOA source function

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, KPARS
          L_TOA_SOURCE(Q,UM) = ZERO
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LPC_TOASOURCE

!

      SUBROUTINE LS_TOASOURCE &
          ( N_USER_STREAMS, N_SURFACEWFS, &
            LS_TOA_SOURCE )

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m,    Only : fpk, ZERO, MAX_SURFACEWFS, MAX_USER_STREAMS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER  , INTENT(IN)  :: N_USER_STREAMS, N_SURFACEWFS
      REAL(FPK), INTENT(OUT) :: LS_TOA_SOURCE(MAX_SURFACEWFS,MAX_USER_STREAMS)

!  local variables

      INTEGER :: UM

!  initialise TOA source function

      DO UM = 1, N_USER_STREAMS
        LS_TOA_SOURCE(1:N_SURFACEWFS,UM) = ZERO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LS_TOASOURCE

!

      SUBROUTINE LPC_BOASOURCE &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,      & ! Inputs
            NSTREAMS, NLAYERS, N_USER_STREAMS, FOURIER, Raman_IDX, Brdf_IDX, & ! Inputs
            K, KPARS, KVARY, QUAD_STRMWGHT, LCON, MCON, XPOS, XNEG, KTRANS,  & ! Inputs
            SURFACE_FACTOR, ALBEDOS_RANKED, USER_BRDF_F, L_USER_DIRECT_BEAM, & ! Inputs
            NCON, PCON, L_XPOS, L_XNEG, L_KTRANS, L_WLOWER,                  & ! Inputs
            L_BOA_SOURCE, L_DIRECT_BOA_SOURCE )                                ! Outputs

!  Bottom of the atmosphere source term

!  This is LRRS Version 2.5. Main changes to this subroutine (from V2.3) are
!    (2) Use of Supplement-derived BRDF inputs and control, for BOASOURCE

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_2_STREAMS, MAX_STREAMS, MAX_USER_STREAMS, &
                              MAX_ATMOSWFS, MAX_MOMENTS, MAX_LAYERS, MAX_POINTS


      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE
      LOGICAL  , INTENT(IN) :: DO_INCLUDE_DIRECTBEAM
      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE

!  Control numbers

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: NLAYERS
      INTEGER  , INTENT(IN) :: N_USER_STREAMS

!  Fourier index

      INTEGER  , INTENT(IN) :: FOURIER

!  Point indices. Version 2.5.
!     Raman_IDX  = APT/CPT for the Albedos_Ranked choice.
!     Brdf_IDX   = 1 (if single-point BRDF in use) or = Raman_Idx (Multiple-point BRDFs)

      INTEGER  , intent(in)  :: Raman_IDX, Brdf_IDX

!  Linearization control

      LOGICAL  , INTENT(IN) :: KVARY
      INTEGER  , INTENT(IN) :: K, KPARS

!  surface multiplier

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Quadrature

      REAL(FPK), INTENT(IN)  :: QUAD_STRMWGHT ( MAX_STREAMS )

!  albedo

      REAL(fpk), intent(in)  :: ALBEDOS_RANKED ( MAX_POINTS )

!  incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: USER_BRDF_F ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Direct beam contribution

      REAL(FPK), INTENT(IN) :: L_USER_DIRECT_BEAM ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  integration constants, homogeneous solutions

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: LCON ( MAX_STREAMS,   MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: MCON ( MAX_STREAMS,   MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Linearizations

      REAL(FPK), INTENT(IN) :: L_WLOWER ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS,   MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: NCON     ( MAX_ATMOSWFS, MAX_STREAMS,   MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: PCON     ( MAX_ATMOSWFS, MAX_STREAMS,   MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XPOS   ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XNEG   ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  output variables, direct and diffuse surface fields.

      REAL(FPK), INTENT(OUT) :: L_BOA_SOURCE        ( MAX_ATMOSWFS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: L_DIRECT_BOA_SOURCE ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER   :: N, I, UM, AA, Q
      REAL(FPK) :: PAR, HOM, REFLEC, LXP, L_LXP, HOM1, HOM2
      REAL(FPK) :: IDOWNSURF ( MAX_ATMOSWFS, MAX_STREAMS )

!  initialise BOA source function

      DO Q = 1, KPARS
        DO UM = 1, N_USER_STREAMS
          L_BOA_SOURCE(Q,UM)        = ZERO
          L_DIRECT_BOA_SOURCE(Q,UM) = ZERO
        ENDDO
      ENDDO

!  Last layer only.

      N = NLAYERS

!  start albedo clause

      IF ( DO_INCLUDE_SURFACE .AND. KVARY ) THEN

!  Get downward intensity at computational angles (beam and homog)
!    And develop reflectance integrand  a(j).x(j).I(-j) (stored in commo

!  If this is also the layer that is varying, extra contributions

        IF ( K .EQ. N .OR. K.EQ.0  ) THEN
          DO Q = 1, KPARS
            DO I = 1, NSTREAMS
              PAR = L_WLOWER(Q,I,N)
              HOM = ZERO
              DO AA = 1, NSTREAMS
                LXP   = LCON(AA,N) * XPOS(I,AA,N)
                L_LXP = NCON(Q,AA,N) *   XPOS(I,AA,N) + &
                        LCON(AA,N)   * L_XPOS(Q,I,AA,N)
                HOM1 = L_LXP *   KTRANS(AA,N) + &
                         LXP * L_KTRANS(Q,AA,N)
                HOM2 =  PCON(Q,AA,N) *   XNEG(I,AA,N) + &
                        MCON(AA,N)   * L_XNEG(Q,I,AA,N)
                HOM = HOM + HOM1 + HOM2
              ENDDO
              IDOWNSURF(Q,I) = QUAD_STRMWGHT(I) * ( PAR + HOM )
            ENDDO
          ENDDO
        ENDIF

!  non-varying layer lower than the varying layer

        IF ( K.LT.N .AND. K.NE.0 ) THEN
          DO Q = 1, KPARS
            DO I = 1, NSTREAMS
              PAR = L_WLOWER(Q,I,N)
              HOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON(Q,AA,N) * XPOS(I,AA,N) * KTRANS(AA,N)
                HOM2 = PCON(Q,AA,N) * XNEG(I,AA,N)
                HOM = HOM + HOM1 + HOM2
              ENDDO
              IDOWNSURF(Q,I) = QUAD_STRMWGHT(I) * ( PAR + HOM )
            ENDDO
          ENDDO
        ENDIF

!  Compute reflectance
!  reflected multiple scatter intensity at user defined-angles

        IF ( DO_BRDF_SURFACE ) THEN
          DO Q = 1, KPARS
            DO UM = 1, N_USER_STREAMS
              REFLEC = DOT_PRODUCT(IDOWNSURF(Q,1:NSTREAMS),USER_BRDF_F(FOURIER,UM,1:NSTREAMS,Brdf_IDX))
              L_BOA_SOURCE(Q,UM) = REFLEC * SURFACE_FACTOR
            ENDDO
          ENDDO
        ELSE 
          IF ( FOURIER .EQ. 0 ) THEN
            DO Q = 1, KPARS
              REFLEC = SUM(IDOWNSURF(Q,1:NSTREAMS)) * SURFACE_FACTOR * ALBEDOS_RANKED(Raman_idx)
              DO UM = 1, N_USER_STREAMS
                L_BOA_SOURCE(Q,UM) = REFLEC
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  Add direct beam if flagged

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          DO Q = 1, KPARS
            DO UM = 1, N_USER_STREAMS
              L_DIRECT_BOA_SOURCE(Q,UM) = L_USER_DIRECT_BEAM(Q,UM)
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LPC_BOASOURCE

!

      SUBROUTINE LS_BOASOURCE &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                         & ! Inputs
            DO_SSCORRECTION, NSTREAMS, NLAYERS, N_USER_STREAMS,          & ! Inputs
            FOURIER, Raman_IDX, Brdf_IDX, N_SURFACEWFS, QUAD_STRMWGHT,   & ! Inputs
            IDOWNSURF, SURFACE_FACTOR, ALBEDOS_RANKED,                   & ! Inputs
            USER_BRDF_F, LS_USER_BRDF_F,  LS_USER_DIRECT_BEAM,           & ! Inputs
            LS_WLOWER, NCON_ALB, PCON_ALB, XPOS, XNEG, KTRANS,           & ! Inputs
            LS_BOA_SOURCE, LS_DIRECT_BOA_SOURCE )                          ! Outputs

!  Bottom of the atmosphere source term

!  This is LRRS Version 2.5. Main changes to this subroutine (from V2.3) are
!    (2) Use of Supplement-derived BRDF inputs and control, for BOASOURCE

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_2_STREAMS, MAX_STREAMS, MAX_USER_STREAMS, &
                              MAX_MOMENTS, MAX_LAYERS, MAX_SURFACEWFS, MAX_POINTS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE
      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE
      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION

!  Control numbers

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: NLAYERS
      INTEGER  , INTENT(IN) :: N_USER_STREAMS

!  Fourier index

      INTEGER  , INTENT(IN) :: FOURIER

!  Point indices. Version 2.5.
!     Raman_IDX  = APT/CPT for the Albedos_Ranked choice.
!     Brdf_IDX   = 1 (if single-point BRDF in use) or = Raman_Idx (Multiple-point BRDFs)

      INTEGER  , intent(in)  :: Raman_IDX, Brdf_IDX

!  Number of surface Jacobians

      INTEGER  , INTENT(IN) :: N_SURFACEWFS

!  Quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT ( MAX_STREAMS )


!  Downwelling radiance integrated

      REAL(FPK), INTENT(IN) :: IDOWNSURF ( MAX_STREAMS )

!  surface multiplier

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  albedo

      REAL(fpk), intent(in)  :: ALBEDOS_RANKED ( MAX_POINTS )

!  incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: USER_BRDF_F ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Fourier components of BRDF, in the following order
!    incident quadrature streams, reflected user streams

      REAL(fpk), dimension ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS ) :: LS_USER_BRDF_F

!  Direct beam contribution

      REAL(FPK), INTENT(IN)  :: LS_USER_DIRECT_BEAM ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  integration constants, homogeneous solutions

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XPOS   ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG   ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Linearizations

      REAL(FPK), INTENT(IN) :: LS_WLOWER ( MAX_SURFACEWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: NCON_ALB  ( MAX_SURFACEWFS, MAX_STREAMS,   MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: PCON_ALB  ( MAX_SURFACEWFS, MAX_STREAMS,   MAX_LAYERS )

!  output variables, Linearized direct and diffuse surface fields.

      REAL(FPK), INTENT(OUT) :: LS_BOA_SOURCE        ( MAX_SURFACEWFS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: LS_DIRECT_BOA_SOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER   :: N, I, UM, AA, M, Q
      REAL(FPK) :: REFLEC, H1, H2, SUMR, SUM1, SUM2
      REAL(FPK) :: LS_IDOWNSURF ( MAX_STREAMS )

!  initialise BOA source functions

      DO UM = 1, N_USER_STREAMS
        LS_BOA_SOURCE(1:N_SURFACEWFS,UM)        = ZERO
        LS_DIRECT_BOA_SOURCE(1:N_SURFACEWFS,UM) = ZERO
      ENDDO

!  Last layer only.

      N = NLAYERS

!  Fourier (for debug only)

      M = FOURIER

!  Get the radiance field at the surface
!    Diffuse and direct contributions. This part is Lambertian-only

      IF ( DO_INCLUDE_SURFACE ) THEN

!  Start loop over surface weighting functions

        DO Q = 1, N_SURFACEWFS

!  First compute derivative of downward intensity Integrand at stream an
!        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

          DO I = 1, NSTREAMS
            SUMR = LS_WLOWER(Q,I,N)
            DO AA = 1, NSTREAMS
              H1  = NCON_ALB(Q,AA,N) * XPOS(I,AA,N) * KTRANS(AA,N)
              H2  = PCON_ALB(Q,AA,N) * XNEG(I,AA,N)
              SUMR = SUMR + H1 + H2
            ENDDO
            LS_IDOWNSURF(I) = SUMR * QUAD_STRMWGHT(I)
          ENDDO

!  integrated reflectance term (Lambertian albedo)
!    --> 2 Contributions from the chain-rule differentiation

          IF ( .not. DO_BRDF_SURFACE .and. M.eq.0 ) THEN
            SUM1 = SUM(IDOWNSURF(1:NSTREAMS))
            SUM2 = SUM(LS_IDOWNSURF(1:NSTREAMS))
            REFLEC = SURFACE_FACTOR * ( SUM1 + ALBEDOS_RANKED(Raman_idx) * SUM2 )
            DO UM = 1, N_USER_STREAMS
               LS_BOA_SOURCE(Q,UM) = REFLEC
            ENDDO
          ENDIF

!  Integrated relfectance (BRDF term)
!    --> 2 Contributions from the chain-rule differentiation

          IF ( DO_BRDF_SURFACE ) THEN
            DO UM = 1, N_USER_STREAMS
              SUM1 = DOT_PRODUCT(   IDOWNSURF(1:NSTREAMS), LS_USER_BRDF_F(Q,M,UM,1:NSTREAMS,BRDF_Idx))
              SUM2 = DOT_PRODUCT(LS_IDOWNSURF(1:NSTREAMS),    USER_BRDF_F(M,UM,1:NSTREAMS,BRDF_Idx))
              REFLEC = ( SUM1 + SUM2 ) * SURFACE_FACTOR
              LS_BOA_SOURCE(Q,UM) = REFLEC
            ENDDO
          ENDIF

!  Add linearization due to Albedo variation of direct beam

          IF ( .NOT. DO_SSCORRECTION ) THEN
            DO UM = 1, N_USER_STREAMS
              LS_DIRECT_BOA_SOURCE(Q,UM) = LS_USER_DIRECT_BEAM (Q,UM)
            ENDDO
          ENDIF

!  End surface-WF do-loop and inclusion clause

        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LS_BOASOURCE

!

      SUBROUTINE LPC_RAMAN_SOURCETERM_UP_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION,                     & ! Inputs
            N, STERM_LAYERMASK_UP, NSTREAMS, N_USER_STREAMS,       & ! Inputs
            TAYLOR_SMALL, TAYLOR_ORDER, KS, KF, KPARS, NKSTORAGE2, & ! Inputs
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION, & ! Inputs
            TRANS_USERM, DELTAUS, USER_STREAMS,                    & ! Inputs
            EMULT_UP, ATERM_SAVE, BTERM_SAVE, SUSAV, SDSAV,        & ! Inputs
            QSOURCE_UPOS1, ASOURCE, DELTRANS, INITRANS,            & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,      & ! Inputs
            L_TRANS_USERM, L_DELTAUS,                              & ! Inputs
            L_EMULT_UP, L_ATERM_SAVE, L_BTERM_SAVE,                & ! Inputs
            L_QSOURCE_UPOS1, L_ASOURCE, L_DELTRANS, L_INITRANS,    & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,    & ! Inputs
            L_WLAYER_PIST_UP )                                       ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m    , Only : FPK, ZERO, ONE, MAX_ATMOSWFS, &
                                  MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS, MAX_LAYERS_SQ

      USE lrrs_Taylor_m, Only : TAYLOR_SERIES_L_2a, TAYLOR_SERIES_L_2b

      IMPLICIT NONE

!  input
!  -----

!  Flags for multiple scattering and single scatter correction to elasti

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT
      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION

!  Fourier component

      INTEGER  , INTENT(IN) :: FOURIER

!  number of streams and solutions

      INTEGER  , INTENT(IN) :: NSTREAMS, N_USER_STREAMS

!  Taylor-series control

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  Linearization control

      INTEGER  , INTENT(IN) :: KS, KF, KPARS(0:MAX_LAYERS)
      INTEGER  , INTENT(IN) :: NKSTORAGE2 ( MAX_LAYERS, 0:MAX_LAYERS )

!  Number of layer, and layer control

      INTEGER  , INTENT(IN) :: N
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP

!  Control for the particular solutions

      LOGICAL  , INTENT(IN) :: CHANGE_EXPONENT
      LOGICAL  , INTENT(IN) :: FLIPPER
      INTEGER  , INTENT(IN) :: SOLUTION

!  Optical thickness values

      REAL(FPK), INTENT(IN) :: DELTAUS

!  stream directions

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  Layer transmittances in stream directions

      REAL(FPK), INTENT(IN) :: TRANS_USERM  ( MAX_USER_STREAMS )

!  Beam attributes to particular solution

      INTEGER  , INTENT(IN) :: PICUTOFF
      REAL(FPK), INTENT(IN) :: ASOURCE
      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: INITRANS

!  Saved ATERM and BTERM contributions to the Green's function.

      REAL(FPK), INTENT(IN) :: ATERM_SAVE  ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: BTERM_SAVE  ( MAX_STREAMS )

!  Eigenvalues and transmittance factors.

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )

!  Whole layer multipliers, homogeneous solutions

      REAL(FPK), INTENT(IN) :: HMULT_1 &
             ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: HMULT_2 &
             ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Eigensolutions at user-defined stream angles

      REAL(FPK), INTENT(IN) :: &
          U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Particular solution attributes

      REAL(FPK), INTENT(IN) :: EMULT_UP      ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: QSOURCE_UPOS1 ( MAX_USER_STREAMS )

!  Saved multiplier arrays

      REAL(FPK), INTENT(IN) :: SUSAV (MAX_STREAMS,MAX_USER_STREAMS)
      REAL(FPK), INTENT(IN) :: SDSAV (MAX_STREAMS,MAX_USER_STREAMS)

!  Linearizations

      REAL(FPK), INTENT(IN) :: L_DELTAUS     ( MAX_ATMOSWFS )
      REAL(FPK), INTENT(IN) :: L_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS )

      REAL(FPK), INTENT(IN) :: L_ASOURCE  ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_DELTRANS ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_INITRANS ( MAX_ATMOSWFS, 0:MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_ATERM_SAVE &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: L_BTERM_SAVE &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS )

      REAL(FPK), INTENT(IN) :: L_KEIGEN ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_HMULT_1 &
          ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_HMULT_2 &
          ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_U_XPOS &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_XNEG &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_EMULT_UP &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: L_QSOURCE_UPOS1 &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS )

!  output
!  ------

!  Linearized Layer source terms
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: L_WLAYER_PIST_UP &
          ( MAX_ATMOSWFS, MAX_LAYERS_SQ, MAX_USER_STREAMS )

!  local variables
!  ---------------

!  help variables

      LOGICAL   :: cedum
      INTEGER   :: AA, UM, Q, M, K, NK2
      REAL(FPK) :: PAR, PAR1, RHOM, RHOP
      REAL(FPK) :: MULT, DMU, UDEL, LT1, LT2, FAC1, FAC2
      REAL(FPK) :: L_PAR, L_PAR1, L_RHOM, L_RHOP
      REAL(FPK) :: L_MULT, L_WDEL, L_SECBAR, L_UDEL

!  Local multiplier arrays

      REAL(FPK) :: SU ( MAX_STREAMS ), SD ( MAX_STREAMS )
      REAL(FPK) :: L_SU    ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_SD    ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_SUSAV ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_SDSAV ( MAX_ATMOSWFS, MAX_STREAMS )

!  Debug

      LOGICAL :: USING_TAYLOR ( MAX_STREAMS )

!  ***************************************************************
!  Single RRS contribution: Greens function multipliers
!  New code programmed with small number analysis and FLIPPER
!   27 September 2006, 29 November 2006. R. Spurr, RT SOLUTIONS.
!  ***************************************************************

!  Initialize (only if SOLUTION = 1)
!mick fix 7/20/2016 - commented out (initializing already done in L_SOURCES_MASTER_1)

      !IF ( SOLUTION .EQ. 1 ) THEN
      !  DO K = KS, KF
      !    NK2 = NKSTORAGE2(N,K)
      !    DO UM = 1, N_USER_STREAMS
      !      DO Q = 1, KPARS(K)
      !        L_WLAYER_PIST_UP(Q,NK2,UM) = ZERO
      !      ENDDO
      !    ENDDO
      !  ENDDO
      !ENDIF

!  Temporary setting

      cedum = CHANGE_EXPONENT

!  Fourier (for debug only)

      M = FOURIER

      IF ( STERM_LAYERMASK_UP ) THEN

!  Start loop over user angles

        DO UM = 1, N_USER_STREAMS

!  Local quantities

          DMU    = ONE/USER_STREAMS(UM)
          UDEL   = TRANS_USERM(UM)
          MULT   = EMULT_UP(UM)

!  Start weighting function loop

          DO K = KS, KF

!  Storage indicator

            NK2 = NKSTORAGE2(N,K)

!  Start parameter loop

            DO Q = 1, KPARS(K)

!  Only require a new multiplier when the exponent changes
!      ----------------------------------------- DISABLE THIS FEATURE
!            IF ( CHANGE_EXPONENT ) THEN

!  set directional multipliers, and their linearizations

              L_UDEL   = L_TRANS_USERM(Q,UM)
              L_MULT   = L_EMULT_UP(Q,K,UM)
              L_SECBAR = L_ASOURCE(Q,K)
              L_WDEL   = L_DELTRANS(Q,K)

!  Zero linearized Greens function multipliers if no solution, and move

              IF ( N .GT. PICUTOFF ) THEN
                DO AA = 1, NSTREAMS
                  L_SDSAV(Q,AA) = ZERO
                  L_SUSAV(Q,AA) = ZERO
                ENDDO
                GO TO 668
              ENDIF

!  Flipper case, Compute linearized multipliers
!    Linearized small numbers routine................10/22/08

              IF ( FLIPPER ) THEN
                IF ( N.EQ.K .OR. K.EQ.0 ) THEN
                  DO AA = 1, NSTREAMS

                    !mick: L(G+-) (flip case)
                    RHOM   = ASOURCE  -   KEIGEN(AA,N)
                    L_RHOM = L_SECBAR - L_KEIGEN(Q,AA,N)

                    IF ( ABS(RHOM) .LT. TAYLOR_SMALL )THEN
                      FAC1 = UDEL ; FAC2 = DELTRANS
                      CALL Taylor_Series_L_2b &
                        ( taylor_order, taylor_small, rhom, keigen(aa,n), deltaus, l_deltaus(q), l_secbar, &
                          l_keigen(q,aa,n), ktrans(aa,n), udel, dmu, asource, fac1, fac2, &
                          l_deltaus(q), l_keigen(q,aa,n), l_susav(q,aa) )
                    ELSE
                      LT1 = L_HMULT_1(Q,AA,UM,N) - L_MULT
                      LT2 = L_RHOM * SUSAV(AA,UM)
                      L_SUSAV(Q,AA)  = ( LT1 - LT2 ) / RHOM
                    ENDIF

                    !mick: L(G++) (flip case)
                    RHOP   =  ASOURCE +   KEIGEN(AA,N)
                    L_RHOP = L_SECBAR + L_KEIGEN(Q,AA,N)
                    LT1 = L_HMULT_2(Q,AA,UM,N) * DELTRANS + &
                            HMULT_2(AA,UM,N)   * L_WDEL
                    LT1 = L_MULT - LT1
                    LT2 = L_RHOP * SDSAV(AA,UM)
                    L_SDSAV(Q,AA)  = ( LT1 - LT2 ) / RHOP
                  ENDDO
                ELSE IF ( N.NE.K .AND. K.GT.0 ) THEN
                  DO AA = 1, NSTREAMS

                    !mick: L(G+-) (flip case)
                    RHOM = ASOURCE - KEIGEN(AA,N)
                    RHOP = ASOURCE + KEIGEN(AA,N)
                    IF ( ABS(RHOM) .LT. TAYLOR_SMALL )THEN
                      FAC1 = UDEL ; FAC2 = DELTRANS
!  Rob 6 feb 2017. Very important Fix
                      if (  abs(L_SECBAR).eq.zero.and.abs(L_mult).eq.zero ) then
                           l_susav(q,aa) = zero
                      else
                        CALL Taylor_Series_L_2b &
                        ( taylor_order, taylor_small, rhom, keigen(aa,n), deltaus, l_deltaus(q), l_secbar, &
                          l_keigen(q,aa,n), ktrans(aa,n), udel, dmu, asource, fac1, fac2, &
                          l_deltaus(q), l_keigen(q,aa,n), l_susav(q,aa) )
                      endif
                    ELSE
                      LT2 = L_SECBAR * SUSAV(AA,UM)
                      L_SUSAV(Q,AA)  = ( - L_MULT - LT2 ) / RHOM
                    ENDIF

                    !mick: L(G++) (flip case)
                    LT1 = L_MULT - HMULT_2(AA,UM,N) * L_WDEL
                    LT2 = L_SECBAR * SDSAV(AA,UM)
                    L_SDSAV(Q,AA)  = ( LT1 - LT2 ) / RHOP
                  ENDDO
                ENDIF
              ENDIF

!  Non Flipper case

              IF ( .NOT.FLIPPER ) THEN
                IF ( N.EQ.K .OR. K.EQ.0 ) THEN
                  DO AA = 1, NSTREAMS

                    !mick: L(G+-)
                    RHOM   =  ASOURCE -   KEIGEN(AA,N)
                    L_RHOM = L_SECBAR - L_KEIGEN(Q,AA,N)

                    IF ( ABS(RHOM) .LT. TAYLOR_SMALL ) THEN
                      FAC1 = ONE ; FAC2 = DELTRANS*UDEL
                      call Taylor_Series_L_2a &
                        ( taylor_order, rhom, keigen(aa,n), deltaus, l_deltaus(q), l_secbar, &
                          l_keigen(q,aa,n), ktrans(aa,n)*udel, dmu, asource, fac1, fac2, &
                          l_deltaus(q), l_keigen(q,aa,n), l_sdsav(q,aa) )
                    ELSE
                      LT1 = L_HMULT_2(Q,AA,UM,N) - L_MULT
                      LT2 = L_RHOM * SDSAV(AA,UM)
                      L_SDSAV(Q,AA)  = ( LT1 - LT2 ) / RHOM
                    ENDIF

                    !mick: L(G++)
                    RHOP   =  ASOURCE +   KEIGEN(AA,N)
                    L_RHOP = L_SECBAR + L_KEIGEN(Q,AA,N)
                    LT1 = L_HMULT_1(Q,AA,UM,N) * DELTRANS + &
                            HMULT_1(AA,UM,N)   * L_WDEL
                    LT1 = L_MULT - LT1
                    LT2 = L_RHOP * SUSAV(AA,UM)
                    L_SUSAV(Q,AA)  = ( LT1 - LT2 ) / RHOP
                  ENDDO
                ELSE IF ( N.NE.K .AND. K.GT.0 ) THEN
                  DO AA = 1, NSTREAMS

                    !mick: L(G+-)
                    RHOM   = ASOURCE -   KEIGEN(AA,N)
                    RHOP   = ASOURCE +   KEIGEN(AA,N)
                    IF ( ABS(RHOM) .LT. TAYLOR_SMALL ) THEN
                      FAC1 = ONE ; FAC2 = DELTRANS*UDEL
!  Rob 6 feb 2017. Very important fix.
                      if (  abs(L_SECBAR).eq.zero.and.abs(L_mult).eq.zero ) then
                           l_sdsav(q,aa) = zero
                      else
                       call Taylor_Series_L_2a &
                        ( taylor_order, rhom, keigen(aa,n), deltaus, zero, l_secbar, &
                          zero, ktrans(aa,n)*udel, dmu, asource, fac1, fac2, &
                          l_deltaus(q), l_keigen(q,aa,n), l_sdsav(q,aa) )
                      endif
                    ELSE
                      LT2 = L_SECBAR * SDSAV(AA,UM)
                      L_SDSAV(Q,AA)  = ( - L_MULT - LT2 ) / RHOM
                    ENDIF

                    !mick: L(G++)
                    LT1 = L_MULT - HMULT_1(AA,UM,N)* L_WDEL
                    LT2 = L_SECBAR * SUSAV(AA,UM)
                    L_SUSAV(Q,AA)  = ( LT1 - LT2 ) / RHOP
                  ENDDO
                ENDIF
              ENDIF

!  Continuation

 668          continue

!  Complete Change-exponent clause
!    Disabled in this version

!            ENDIF

!  End parameter loop

            ENDDO

!  Add the terms independent of optical depth

            DO AA = 1, NSTREAMS
              SD(AA) = SDSAV(AA,UM) * ATERM_SAVE(AA)
              SU(AA) = SUSAV(AA,UM) * BTERM_SAVE(AA)
              DO Q = 1, KPARS(K)
                L_SD(Q,AA) = L_SDSAV(Q,AA)  *   ATERM_SAVE(AA) &
                             + SDSAV(AA,UM) * L_ATERM_SAVE(Q,K,AA)
                L_SU(Q,AA) = L_SUSAV(Q,AA)  *   BTERM_SAVE(AA) &
                             + SUSAV(AA,UM) * L_BTERM_SAVE(Q,K,AA)
              ENDDO
            ENDDO

!  Layer source computation
!  ========================

!  Summation over all streams for the PAR contrribution

            PAR = ZERO
            DO AA = 1, NSTREAMS
              PAR = PAR + U_XPOS(UM,AA,N)*SD(AA)+U_XNEG(UM,AA,N)*SU(AA)
            ENDDO

!  Calculate L_PAR and upgrade the source term
!   For the self-correlated profile WF or the Column WF

            IF ( N.EQ.K .OR. K.EQ.0 ) THEN
              DO Q = 1, KPARS(K)
                L_PAR = ZERO
                DO AA = 1, NSTREAMS
                  L_PAR = L_PAR &
                    + L_U_XPOS(Q,UM,AA,N)*SD(AA) + U_XPOS(UM,AA,N)*L_SD(Q,AA) &
                    + L_U_XNEG(Q,UM,AA,N)*SU(AA) + U_XNEG(UM,AA,N)*L_SU(Q,AA)
                ENDDO
                L_WLAYER_PIST_UP(Q,NK2,UM) = L_WLAYER_PIST_UP(Q,NK2,UM) &
                  + PAR * L_INITRANS(Q,K) + L_PAR * INITRANS
              ENDDO
            ENDIF

!  Note these smart debugs
!       if(k.eq.11.and.q.eq.2)write(65,6)n,solution,flipper,
!     &         par*L_INITRANS(Q,K)+L_PAR*INITRANS
!       if(k.eq.11.and.q.eq.2)write(75,6)n,solution,flipper,
!     &         par*L_INITRANS(Q,K)+L_PAR*INITRANS
! 6        format(2i4,L2,1pe20.10)

!  Calculate L_PAR and upgrade the source term
!   For the self-correlated profile WF or the Column WF

            IF ( N.NE.K .AND. K.GT.0 ) THEN
              DO Q = 1, KPARS(K)
                L_PAR = ZERO
                DO AA = 1, NSTREAMS
                  L_PAR = L_PAR + U_XPOS(UM,AA,N)*L_SD(Q,AA) &
                                + U_XNEG(UM,AA,N)*L_SU(Q,AA)
                ENDDO
                L_WLAYER_PIST_UP(Q,NK2,UM) = L_WLAYER_PIST_UP(Q,NK2,UM) &
                    + PAR * L_INITRANS(Q,K) + L_PAR * INITRANS
              ENDDO
            ENDIF

!  End weighting function loop

          ENDDO

!  End main user stream loop

        ENDDO

!  Add single scatter contribution if flagged.
!  If you are doing an SS correction, avoid this part
!  .. Otherwise, Full radiance mode, add single scatter part
!   2/6/17. Get rid of GOTO 5555, replace by If block

        IF ( .NOT. DO_MSMODE_LIDORT ) THEN
          IF ( .not. DO_SSCORRECTION ) THEN
            DO UM = 1, N_USER_STREAMS
              PAR1 = EMULT_UP(UM) * QSOURCE_UPOS1(UM)
              DO K = KS, KF
                NK2 = NKSTORAGE2(N,K)
                DO Q = 1, KPARS(K)
                  L_PAR1 = L_EMULT_UP(Q,K,UM) *   QSOURCE_UPOS1(UM) &
                          +  EMULT_UP(UM)     * L_QSOURCE_UPOS1(Q,K,UM)
                  L_WLAYER_PIST_UP(Q,NK2,UM) = L_WLAYER_PIST_UP(Q,NK2,UM) &
                    + PAR1 * L_INITRANS(Q,K) + L_PAR1 * INITRANS
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  End layer

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LPC_RAMAN_SOURCETERM_UP_1

!

      SUBROUTINE LPC_RAMAN_SOURCETERM_DN_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION,                     & ! Inputs
            N, STERM_LAYERMASK_DN, NSTREAMS, N_USER_STREAMS,       & ! Inputs
            TAYLOR_SMALL, TAYLOR_ORDER, KS, KF, KPARS, NKSTORAGE2, & ! Inputs
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION, & ! Inputs
            TRANS_USERM, DELTAUS, USER_STREAMS,                    & ! Inputs
            EMULT_DN, ATERM_SAVE, BTERM_SAVE, SUSAV, SDSAV,        & ! Inputs
            QSOURCE_UNEG1, ASOURCE, DELTRANS, INITRANS,            & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,      & ! Inputs
            L_TRANS_USERM, L_DELTAUS,                              & ! Inputs
            L_EMULT_DN, L_ATERM_SAVE, L_BTERM_SAVE,                & ! Inputs
            L_QSOURCE_UNEG1, L_ASOURCE, L_DELTRANS, L_INITRANS,    & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,    & ! Inputs
            L_WLAYER_PIST_DN )                                       ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m    , Only : FPK, ZERO, ONE, MAX_ATMOSWFS, &
                                MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS, MAX_LAYERS_SQ

      USE lrrs_Taylor_m, Only : TAYLOR_SERIES_L_2a, TAYLOR_SERIES_L_2b

      IMPLICIT NONE

!  input
!  -----

!  Flags for multiple scattering and single scatter correction to elastic

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT
      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION

!  Fourier component

      INTEGER  , INTENT(IN) :: FOURIER

!  number of streams and solutions

      INTEGER  , INTENT(IN) :: NSTREAMS, N_USER_STREAMS

!  Taylor-series control

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  Linearization control

      INTEGER  , INTENT(IN) :: KS, KF, KPARS(0:MAX_LAYERS)
      INTEGER  , INTENT(IN) :: NKSTORAGE2 ( MAX_LAYERS, 0:MAX_LAYERS )

!  Number of layer, and layer control

      INTEGER  , INTENT(IN) :: N
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN

!  Control for the particular solutions

      LOGICAL  , INTENT(IN) :: CHANGE_EXPONENT
      LOGICAL  , INTENT(IN) :: FLIPPER
      INTEGER  , INTENT(IN) :: SOLUTION

!  Optical thickness value

      REAL(FPK), INTENT(IN) :: DELTAUS

!  stream directions

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  Layer transmittances in stream directions

      REAL(FPK), INTENT(IN) :: TRANS_USERM  ( MAX_USER_STREAMS )

!  Beam attributes to particular solution

      INTEGER  , INTENT(IN) :: PICUTOFF
      REAL(FPK), INTENT(IN) :: ASOURCE
      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: INITRANS

!  Saved ATERM and BTERM contributions to the Green's function.

      REAL(FPK), INTENT(IN) :: ATERM_SAVE  ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: BTERM_SAVE  ( MAX_STREAMS )

!  Eigenvalues and transmittance factors.

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )

!  Whole layer multipliers, homogeneous solutions

      REAL(FPK), INTENT(IN) :: HMULT_1 &
             ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: HMULT_2 &
             ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Eigensolutions at user-defined stream angles

      REAL(FPK), INTENT(IN) :: &
          U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Particular solution attributes

      REAL(FPK), INTENT(IN) :: EMULT_DN      ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: QSOURCE_UNEG1 ( MAX_USER_STREAMS )

!  Saved multiplier arrays

      REAL(FPK), INTENT(IN) :: SUSAV (MAX_STREAMS,MAX_USER_STREAMS)
      REAL(FPK), INTENT(IN) :: SDSAV (MAX_STREAMS,MAX_USER_STREAMS)

!  Linearizations

      REAL(FPK), INTENT(IN) :: L_DELTAUS     ( MAX_ATMOSWFS )
      REAL(FPK), INTENT(IN) :: L_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS )

      REAL(FPK), INTENT(IN) :: L_ASOURCE      ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_DELTRANS     ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_INITRANS     ( MAX_ATMOSWFS, 0:MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_ATERM_SAVE &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: L_BTERM_SAVE &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS )

      REAL(FPK), INTENT(IN) :: L_KEIGEN ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_HMULT_1 &
          ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_HMULT_2 &
          ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_U_XPOS &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_XNEG &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_EMULT_DN &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS , MAX_USER_STREAMS)
      REAL(FPK), INTENT(IN) :: L_QSOURCE_UNEG1 &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS )

!  output
!  ------

!  Linearized Layer source terms
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: L_WLAYER_PIST_DN &
          ( MAX_ATMOSWFS, MAX_LAYERS_SQ, MAX_USER_STREAMS )

!  local variables
!  ---------------

!  help variables

      LOGICAL :: cedum
      INTEGER :: AA, UM, Q, M, K, NK2
      REAL(FPK) :: PAR, PAR1, RHOM, RHOP
      REAL(FPK) :: MULT, DMU, UDEL, LT1, LT2, FAC1, FAC2
      REAL(FPK) :: L_PAR, L_PAR1, L_RHOM, L_RHOP
      REAL(FPK) :: L_MULT, L_WDEL, L_SECBAR, L_UDEL

!  Local multiplier arrays

      REAL(FPK) :: SU ( MAX_STREAMS ), SD ( MAX_STREAMS )
      REAL(FPK) :: L_SU    ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_SD    ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_SUSAV ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_SDSAV ( MAX_ATMOSWFS, MAX_STREAMS )

!  ***************************************************************
!  Single RRS constribution: Greens function multipliers
!  New code programmed with small number analysis and FLIPPER
!   27 September 2006, 29 November 2006. R. Spurr, RT SOLUTIONS.
!  ***************************************************************

!  Initialize (only if SOLUTION = 1)
!mick fix 7/20/2016 - commented out (initializing already done in L_SOURCES_MASTER_1)

      !IF ( SOLUTION .EQ. 1 ) THEN
      !  DO K = KS, KF
      !    NK2 = NKSTORAGE2(N,K)
      !    DO UM = 1, N_USER_STREAMS
      !      DO Q = 1, KPARS(K)
      !        L_WLAYER_PIST_DN(Q,NK2,UM) = ZERO
      !      ENDDO
      !    ENDDO
      !  ENDDO
      !ENDIF

!  Temporary setting

      cedum = CHANGE_EXPONENT

!  Fourier (for debug only)

      M = FOURIER

      IF ( STERM_LAYERMASK_DN ) THEN

!  Start loop over user angles

        DO UM = 1, N_USER_STREAMS

!  Local quantities

          DMU    = ONE/USER_STREAMS(UM)
          UDEL   = TRANS_USERM(UM)
          MULT   = EMULT_DN(UM)

!  Start weighting function loop

          DO K = KS, KF

!  Set storage indicator

            NK2 = NKSTORAGE2(N,K)

!  Start parameter loop

            DO Q = 1, KPARS(K)

!  Local quantities


!  Only require a new multiplier when the exponent changes
!   -----------------Disabled this feature---------------
!            IF ( CHANGE_EXPONENT ) THEN

!  set directional multipliers, and their linearizations

              L_UDEL   = L_TRANS_USERM(Q,UM)
              L_MULT   = L_EMULT_DN(Q,K,UM)
              L_SECBAR = L_ASOURCE(Q,K)
              L_WDEL   = L_DELTRANS(Q,K)

!  Zero linearized Greens function multipliers if no solution, and move

              IF ( N .GT. PICUTOFF ) THEN
                DO AA = 1, NSTREAMS
                  L_SDSAV(Q,AA) = ZERO
                  L_SUSAV(Q,AA) = ZERO
                ENDDO
                GO TO 668
              ENDIF

!  Flipper case, Compute linearized multipliers
!  Former Calls were.....
!                      call limit_l_ppgmult_wu &
!                       ( rhom, dmu, one, keigen(aa,n), deltaus, ktrans(aa,n), udel, &
!                         l_secbar, l_keigen(q,aa,n), l_deltaus(q), l_susav(q,aa) )
                       !call Taylor_Series_L_2a &
                       ! ( taylor_order, rhom, keigen(aa,n), deltaus, l_deltaus(q), l_secbar, &
                       !   l_keigen(q,aa,n), ktrans(aa,n)*udel, dmu, l_susav(q,aa) )

              IF ( FLIPPER ) THEN
                IF ( N.EQ.K .OR. K.EQ.0 ) THEN
                  DO AA = 1, NSTREAMS

                    !mick: L(G--) (flip case)
                    RHOM   =  ASOURCE -   KEIGEN(AA,N)
                    L_RHOM = L_SECBAR - L_KEIGEN(Q,AA,N)
                    IF ( ABS(RHOM) .LT. TAYLOR_SMALL )THEN
                       FAC1 = ONE ; FAC2 = DELTRANS*UDEL
                       call Taylor_Series_L_2a &
                        ( taylor_order, rhom, keigen(aa,n), deltaus, l_deltaus(q), l_secbar, &
                          l_keigen(q,aa,n), ktrans(aa,n)*udel, dmu, asource, fac1, fac2, &
                          l_deltaus(q), l_keigen(q,aa,n), l_susav(q,aa) )
                    ELSE
                     LT1 = L_HMULT_2(Q,AA,UM,N) - L_MULT
                     LT2 = L_RHOM * SUSAV(AA,UM)
                     L_SUSAV(Q,AA)  = ( LT1 - LT2 ) / RHOM
                    ENDIF

                    !mick: L(G-+) (flip case)
                    RHOP   =  ASOURCE +   KEIGEN(AA,N)
                    L_RHOP = L_SECBAR + L_KEIGEN(Q,AA,N)
                    LT1 = L_HMULT_1(Q,AA,UM,N) * DELTRANS + &
                            HMULT_1(AA,UM,N)   * L_WDEL
                    LT1 = L_MULT - LT1
                    LT2 = L_RHOP * SDSAV(AA,UM)
                    L_SDSAV(Q,AA)  = ( LT1 - LT2 ) / RHOP
                  ENDDO
                ELSE IF ( N.NE.K .AND. K.GT.0 ) THEN
                  DO AA = 1, NSTREAMS

 
                    !mick: L(G--) (flip case)
                    RHOM   = ASOURCE -   KEIGEN(AA,N)
                    RHOP   = ASOURCE +   KEIGEN(AA,N)

                    IF ( ABS(RHOM) .LT. TAYLOR_SMALL )THEN
                       FAC1 = ONE ; FAC2 = DELTRANS*UDEL
!  Rob 6 feb 2017. This is a very important fix...
                      if (  abs(L_SECBAR).eq.zero.and.abs(L_mult).eq.zero ) then
                         l_susav(q,aa) = zero
                      else
                       call Taylor_Series_L_2a &
                        ( taylor_order, rhom, keigen(aa,n), deltaus, zero, l_secbar, &
                          zero, ktrans(aa,n)*udel, dmu, asource, fac1, fac2, &
                          l_deltaus(q), l_keigen(q,aa,n), l_susav(q,aa) )
                      endif
                    ELSE
                     LT2 = L_SECBAR * SUSAV(AA,UM)
                     L_SUSAV(Q,AA)  = ( - L_MULT - LT2 ) / RHOM
                    ENDIF

                    !mick: L(G-+) (flip case)
                    LT1 = L_MULT - HMULT_1(AA,UM,N) * L_WDEL
                    LT2 = L_SECBAR * SDSAV(AA,UM)
                    L_SDSAV(Q,AA)  = ( LT1 - LT2 ) / RHOP
                  ENDDO
                ENDIF
              ENDIF

!  Non Flipper case
!   Former Calls
!                      call limit_l_ppgmult_wd &
!                       ( rhom, dmu, one, keigen(aa,n), deltaus, ktrans(aa,n), udel, &
!                         l_secbar, l_keigen(q,aa,n), l_deltaus(q), l_sdsav(q,aa) )
                      !CALL Taylor_Series_L_2b &
                      !  ( taylor_order, rhom, keigen(aa,n), deltaus, l_deltaus(q), l_secbar, &
                      !    l_keigen(q,aa,n), ktrans(aa,n), udel, dmu, asource, l_sdsav(q,aa) )

              IF ( .NOT.FLIPPER ) THEN
                IF ( N.EQ.K .OR. K.EQ.0 ) THEN
                  DO AA = 1, NSTREAMS

                    !mick: L(G--)
                    RHOM   =  ASOURCE -   KEIGEN(AA,N)
                    L_RHOM = L_SECBAR - L_KEIGEN(Q,AA,N)
                    IF ( ABS(RHOM) .LT. TAYLOR_SMALL ) THEN

                      FAC1 = UDEL ; FAC2 = DELTRANS
                      CALL Taylor_Series_L_2b &
                        ( taylor_order, taylor_small, rhom, keigen(aa,n), deltaus, l_deltaus(q), l_secbar, &
                          l_keigen(q,aa,n), ktrans(aa,n), udel, dmu, asource, fac1, fac2, &
                          l_deltaus(q), l_keigen(q,aa,n), l_sdsav(q,aa) )
                    ELSE
                      LT1 = L_HMULT_1(Q,AA,UM,N) - L_MULT
                      LT2 = L_RHOM * SDSAV(AA,UM)
                      L_SDSAV(Q,AA)  = ( LT1 - LT2 ) / RHOM
                    ENDIF

                    !mick: G-+
                    RHOP   =  ASOURCE +   KEIGEN(AA,N)
                    L_RHOP = L_SECBAR + L_KEIGEN(Q,AA,N)
                    LT1 = L_HMULT_2(Q,AA,UM,N) * DELTRANS + &
                            HMULT_2(AA,UM,N)   * L_WDEL
                    LT1 = L_MULT - LT1
                    LT2 = L_RHOP * SUSAV(AA,UM)
                    L_SUSAV(Q,AA)  = ( LT1 - LT2 ) / RHOP
                  ENDDO
                ELSE IF ( N.NE.K .AND. K.GT.0 ) THEN
                  DO AA = 1, NSTREAMS
                    RHOM   = ASOURCE -   KEIGEN(AA,N)
                    RHOP   = ASOURCE +   KEIGEN(AA,N)

                   !mick: L(G--)
                    IF ( ABS(RHOM) .LT. TAYLOR_SMALL ) THEN
                      FAC1 = UDEL ; FAC2 = DELTRANS
!  Rob 6 feb 2017. Important Fix.
                      if (  abs(L_SECBAR).eq.zero.and.abs(L_mult).eq.zero ) then
                           l_sdsav(q,aa) = zero
                      else
                        CALL Taylor_Series_L_2b &
                        ( taylor_order, taylor_small, rhom, keigen(aa,n), deltaus, zero, l_secbar, &
                          zero, ktrans(aa,n), udel, dmu, asource, fac1, fac2, &
                          l_deltaus(q), l_keigen(q,aa,n), l_sdsav(q,aa) )
                      endif
                    ELSE
                      LT2 = L_SECBAR * SDSAV(AA,UM)
                      L_SDSAV(Q,AA)  = ( - L_MULT - LT2 ) / RHOM
                    ENDIF

                    !mick: L(G-+)
                    LT1 = L_MULT - HMULT_2(AA,UM,N) * L_WDEL
                    LT2 = L_SECBAR * SUSAV(AA,UM)
                    L_SUSAV(Q,AA)  = ( LT1 - LT2 ) / RHOP
                  ENDDO
                ENDIF
              ENDIF

!  Continuation

 668          continue

!  Complete Change-exponent clause
!    DISABLED--------------
!            ENDIF

!  End parameter loop

            ENDDO

!  Add the terms independent of optical depth

            DO AA = 1, NSTREAMS
              SD(AA) = SDSAV(AA,UM) * ATERM_SAVE(AA)
              SU(AA) = SUSAV(AA,UM) * BTERM_SAVE(AA)
              DO Q = 1, KPARS(K)
                L_SD(Q,AA) = L_SDSAV(Q,AA)  *   ATERM_SAVE(AA) &
                             + SDSAV(AA,UM) * L_ATERM_SAVE(Q,K,AA)
                L_SU(Q,AA) = L_SUSAV(Q,AA)  *   BTERM_SAVE(AA) &
                             + SUSAV(AA,UM) * L_BTERM_SAVE(Q,K,AA)
              ENDDO
            ENDDO

!  Layer source computation
!  ========================

!  Summation over all streams for the PAR contribution

            PAR = ZERO
            DO AA = 1, NSTREAMS
              PAR = PAR + U_XNEG(UM,AA,N)*SD(AA)+U_XPOS(UM,AA,N)*SU(AA)
            ENDDO

!  Calculate L_PAR and upgrade the source term
!   For the self-correlated profile WF or the Column WF

            IF ( N.EQ.K .OR. K.EQ.0 ) THEN
              DO Q = 1, KPARS(K)
                L_PAR = ZERO
                DO AA = 1, NSTREAMS
                  L_PAR = L_PAR &
              + L_U_XNEG(Q,UM,AA,N)*SD(AA) + U_XNEG(UM,AA,N)*L_SD(Q,AA) &
              + L_U_XPOS(Q,UM,AA,N)*SU(AA) + U_XPOS(UM,AA,N)*L_SU(Q,AA)
                ENDDO
                L_WLAYER_PIST_DN(Q,NK2,UM) = L_WLAYER_PIST_DN(Q,NK2,UM) &
                    + PAR * L_INITRANS(Q,K) + L_PAR * INITRANS
              ENDDO
            ENDIF

!  Calculate L_PAR and upgrade the source term
!   For the self-correlated profile WF or the Column WF

            IF ( N.NE.K .AND. K.GT.0 ) THEN
              DO Q = 1, KPARS(K)
                L_PAR = ZERO
                DO AA = 1, NSTREAMS
                  L_PAR = L_PAR + U_XNEG(UM,AA,N)*L_SD(Q,AA) &
                                + U_XPOS(UM,AA,N)*L_SU(Q,AA)
                ENDDO
                L_WLAYER_PIST_DN(Q,NK2,UM) = L_WLAYER_PIST_DN(Q,NK2,UM) &
                    + PAR * L_INITRANS(Q,K) + L_PAR * INITRANS
              ENDDO
            ENDIF

!  End weighting function loop

          ENDDO

!  End main user stream loop

        ENDDO

!  Add single scatter contribution if flagged.
!  If you are doing an SS correction, avoid this part
!  .. Otherwise, Full radiance mode, add single scatter part
!   2/6/17. Get rid of GOTO 5555, replace by If block

        IF ( .NOT. DO_MSMODE_LIDORT ) THEN
          IF ( .not. DO_SSCORRECTION ) THEN
            DO UM = 1, N_USER_STREAMS
              PAR1 = EMULT_DN(UM) * QSOURCE_UNEG1(UM)
              DO K = KS, KF
                NK2 = NKSTORAGE2(N,K)
                DO Q = 1, KPARS(K)
                  L_PAR1 = L_EMULT_DN(Q,K,UM) *   QSOURCE_UNEG1(UM) &
                          +  EMULT_DN(UM)     * L_QSOURCE_UNEG1(Q,K,UM)
                  L_WLAYER_PIST_DN(Q,NK2,UM) = L_WLAYER_PIST_DN(Q,NK2,UM) &
                      + PAR1 * L_INITRANS(Q,K) + L_PAR1 * INITRANS
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  End layer

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LPC_RAMAN_SOURCETERM_DN_1

!

      SUBROUTINE LS_RAMAN_SOURCETERM_UP_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION,               & ! Inputs
            N, STERM_LAYERMASK_UP, NSTREAMS, N_USER_STREAMS, & ! Inputs
            N_SURFACEWFS, FOURIER, SOLUTION, INITRANS,       & ! Inputs
            U_XPOS, U_XNEG, EMULT_UP, SUSAV, SDSAV,          & ! Inputs
            LS_ATERM_SAVE, LS_BTERM_SAVE, LS_QSOURCE_UPOS1,  & ! Inputs
            LS_WLAYER_PIST_UP )                                ! Output

!  Source terms for the Surface linearization
!    New routine, 19 November 2008

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS, MAX_SURFACEWFS

      IMPLICIT NONE

!  input
!  -----

!  Flags for multiple scattering and single scatter correction

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT
      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION

!  Fourier component (only to help with debugging)

      INTEGER  , INTENT(IN) :: FOURIER

!  number of streams

      INTEGER  , INTENT(IN) :: NSTREAMS, N_USER_STREAMS

!  Number of Surface weighting functions

      INTEGER  , INTENT(IN) :: N_SURFACEWFS

!  Solution number

      INTEGER  , INTENT(IN) :: SOLUTION

!  Number of layer, and layer control

      INTEGER  , INTENT(IN) :: N
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP

!  Beam attributes to particular solution

      REAL(FPK), INTENT(IN) :: INITRANS

!  Eigensolutions at user-defined stream angles

      REAL(FPK), INTENT(IN) :: &
          U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Particular integral multiplier

      REAL(FPK), INTENT(IN) :: EMULT_UP ( MAX_USER_STREAMS )

!  Saved multiplier arrays

      REAL(FPK), INTENT(IN) :: SUSAV (MAX_STREAMS,MAX_USER_STREAMS)
      REAL(FPK), INTENT(IN) :: SDSAV (MAX_STREAMS,MAX_USER_STREAMS)

!  Linearization source terms

      REAL(FPK), INTENT(IN) :: LS_ATERM_SAVE    ( MAX_SURFACEWFS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: LS_BTERM_SAVE    ( MAX_SURFACEWFS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: LS_QSOURCE_UPOS1 ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  output
!  ------

!  Linearized Layer source terms
!mick fix 7/20/2016 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: LS_WLAYER_PIST_UP (MAX_SURFACEWFS,MAX_LAYERS,MAX_USER_STREAMS)

!  local variables
!  ---------------

!  help variables

      INTEGER   :: AA, UM, NS, M
      REAL(FPK) :: LS_PAR, LS_PAR1, LS_SD, LS_SU

!  Initialize (only if SOLUTION = 1)

!mick fix 7/20/2016 - commented out (initializing already done in L_SOURCES_MASTER_1)
      !IF ( SOLUTION .EQ. 1 ) THEN
      !  DO UM = 1, N_USER_STREAMS
      !    LS_WLAYER_PIST_UP(1:N_SURFACEWFS,N,UM) = ZERO
      !  ENDDO
      !ENDIF

!  Fourier (for debug only)

      M = FOURIER

      IF ( STERM_LAYERMASK_UP ) THEN

        DO NS = 1, N_SURFACEWFS

!  For each user angle, Sum over all streams (LS_PAR), then Upgrade

          DO UM = 1, N_USER_STREAMS
            LS_PAR = ZERO
            DO AA = 1, NSTREAMS
              LS_SD = SDSAV(AA,UM) * LS_ATERM_SAVE(NS,AA)
              LS_SU = SUSAV(AA,UM) * LS_BTERM_SAVE(NS,AA)
              LS_PAR = LS_PAR + U_XPOS(UM,AA,N) * LS_SD &
                              + U_XNEG(UM,AA,N) * LS_SU
            ENDDO
            LS_WLAYER_PIST_UP(NS,N,UM) = LS_WLAYER_PIST_UP(NS,N,UM) + LS_PAR * INITRANS
          ENDDO

!  Add single scatter contribution if flagged.
!  If you are doing an SS correction, avoid this part
!  .. Otherwise, Full radiance mode, add single scatter part
!   2/6/17. Get rid of GOTO 5555, replace by If block

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            IF ( .not.DO_SSCORRECTION ) then
              DO UM = 1, N_USER_STREAMS
                LS_PAR1 = EMULT_UP(UM) * LS_QSOURCE_UPOS1(NS,UM)
                LS_WLAYER_PIST_UP(NS,N,UM) = LS_WLAYER_PIST_UP(NS,N,UM) + LS_PAR1 * INITRANS
              ENDDO
            ENDIF
          ENDIF

!  ENd WF loop

        ENDDO

!  End layer

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LS_RAMAN_SOURCETERM_UP_1

!

      SUBROUTINE LS_RAMAN_SOURCETERM_DN_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION,               & ! Inputs
            N, STERM_LAYERMASK_DN, NSTREAMS, N_USER_STREAMS, & ! Inputs
            N_SURFACEWFS, FOURIER, SOLUTION, INITRANS,       & ! Inputs
            U_XPOS, U_XNEG, EMULT_DN, SUSAV, SDSAV,          & ! Inputs
            LS_ATERM_SAVE, LS_BTERM_SAVE, LS_QSOURCE_UNEG1,  & ! Inputs
            LS_WLAYER_PIST_DN )                                ! Output

!  Source terms for the Surface linearization
!    New routine, 19 November 2008

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS, MAX_SURFACEWFS

      IMPLICIT NONE

!  input
!  -----

!  Flags for multiple scattering and single scatter correction

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT
      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION

!  Fourier component (only to help with debugging)

      INTEGER  , INTENT(IN) :: FOURIER

!  number of streams

      INTEGER  , INTENT(IN) :: NSTREAMS, N_USER_STREAMS

!  Number of Surface weighting functions

      INTEGER  , INTENT(IN) :: N_SURFACEWFS

!  Solution number

      INTEGER  , INTENT(IN) :: SOLUTION

!  Number of layer, and layer control

      INTEGER  , INTENT(IN) :: N
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN

!  Beam attributes to particular solution

      REAL(FPK), INTENT(IN) :: INITRANS

!  Eigensolutions at user-defined stream angles

      REAL(FPK), INTENT(IN) :: &
          U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Particular integral multiplier

      REAL(FPK), INTENT(IN) :: EMULT_DN ( MAX_USER_STREAMS )

!  Saved multiplier arrays

      REAL(FPK), INTENT(IN) :: SUSAV (MAX_STREAMS,MAX_USER_STREAMS)
      REAL(FPK), INTENT(IN) :: SDSAV (MAX_STREAMS,MAX_USER_STREAMS)

!  Linearization source terms

      REAL(FPK), INTENT(IN) :: LS_ATERM_SAVE    ( MAX_SURFACEWFS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: LS_BTERM_SAVE    ( MAX_SURFACEWFS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: LS_QSOURCE_UNEG1 ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  output
!  ------

!  Linearized Layer source terms

!mick fix 7/20/2016 - changed intent to inout
      REAL(FPK), INTENT(INOUT) :: LS_WLAYER_PIST_DN (MAX_SURFACEWFS,MAX_LAYERS,MAX_USER_STREAMS)

!  local variables
!  ---------------

!  help variables

      INTEGER   :: AA, UM, NS, M
      REAL(FPK) :: LS_PAR, LS_PAR1, LS_SD, LS_SU

!  Initialize (only if SOLUTION = 1)
!mick fix 7/20/2016 - commented out (initializing already done in L_SOURCES_MASTER_1)

      !IF ( SOLUTION .EQ. 1 ) THEN
      !  DO UM = 1, N_USER_STREAMS
      !    LS_WLAYER_PIST_DN(1:N_SURFACEWFS,N,UM) = ZERO
      !  ENDDO
      !ENDIF

!  Fourier (for debug only)

      M = FOURIER

      IF ( STERM_LAYERMASK_DN ) THEN

!  For each user angle, Sum over all streams (LS_PAR), then Upgrade

        DO NS = 1, N_SURFACEWFS
          DO UM = 1, N_USER_STREAMS
            LS_PAR = ZERO
            DO AA = 1, NSTREAMS
              LS_SD = SDSAV(AA,UM) * LS_ATERM_SAVE(NS,AA)
              LS_SU = SUSAV(AA,UM) * LS_BTERM_SAVE(NS,AA)
              LS_PAR = LS_PAR + U_XNEG(UM,AA,N) * LS_SD &
                              + U_XPOS(UM,AA,N) * LS_SU
            ENDDO
            LS_WLAYER_PIST_DN(NS,N,UM) = LS_WLAYER_PIST_DN(NS,N,UM) + LS_PAR * INITRANS
          ENDDO

!  Add single scatter contribution if flagged.
!  If you are doing an SS correction, avoid this part
!  .. Otherwise, Full radiance mode, add single scatter part
!   2/6/17. Get rid of GOTO 5555, replace by If block

          IF ( .NOT. DO_MSMODE_LIDORT ) THEN
            IF ( .not. DO_SSCORRECTION ) then
              DO UM = 1, N_USER_STREAMS
                LS_PAR1 = EMULT_DN(UM) * LS_QSOURCE_UNEG1(NS,UM)
                LS_WLAYER_PIST_DN(NS,N,UM) = LS_WLAYER_PIST_DN(NS,N,UM) + LS_PAR1 * INITRANS
              ENDDO
            ENDIF
          ENDIF

!  End WF loop

        ENDDO

!  End layer

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LS_RAMAN_SOURCETERM_DN_1

!  End Module

      END MODULE lrrs_L_postprocessing_1_m

