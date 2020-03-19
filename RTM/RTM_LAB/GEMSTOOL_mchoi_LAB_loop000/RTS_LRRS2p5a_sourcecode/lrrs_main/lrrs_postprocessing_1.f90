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
! #   SUBROUTINES :                                             #
! #                                                             #
! #    Master subroutine(used by Raman Masters)                 #
! #                                                             #
! #             POSTPROCESSING_1                                #
! #                                                             #
! #    Public subroutines (used by Elastic & Raman Masters)     #
! #                                                             #
! #             HOMOGMULT_1                                     #
! #             BEAMMULT_UP_1                                   #
! #             BEAMMULT_DN_1                                   #
! #             TOASOURCE                                       #
! #             BOASOURCE                                       #
! #                                                             #
! #    Private subroutines (called only by POSTPROCESSING_1)    #
! #                                                             #
! #             RAMAN_SOURCETERM_UP_1                           #
! #             RAMAN_SOURCETERM_DN_1                           #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Use of Taylor-series subroutines/module, TAYLOR_ORDER index, in "MULT" routines
!    (2) Use of Supplement-derived BRDF inputs and control, for BOASOURCE
!    (3) Bookkeeping improvements (use of "Only", clearer I/O specifications)

      MODULE lrrs_postprocessing_1_m

!      USE LRRS_PARS_m, Only : SDU

      PRIVATE:: RAMAN_SOURCETERM_UP_1, RAMAN_SOURCETERM_DN_1

      PUBLIC :: POSTPROCESSING_1, HOMOGMULT_1,   &
                BEAMMULT_UP_1,    BEAMMULT_DN_1, &
                TOASOURCE,        BOASOURCE

      CONTAINS

      SUBROUTINE POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LIDORT, DO_SSCORRECTION,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,               & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, SOLUTION,     & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Source term post processing for one layer

!  Use Module of dimensions and numbers
!  ------------------------------------

      USE LRRS_PARS_m, Only : FPK, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS

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

!  layer control

      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN

!  Taylor-series control. Updated, Version 2.5, 9/10/15

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

!      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )

!  Whole layer multipliers, homogeneous solutions

      REAL(FPK), INTENT(IN) :: HMULT_1 &
             ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: HMULT_2 &
             ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  output
!  ------

!  Particular integral multipliers

      REAL(FPK), INTENT(OUT) :: EMULT_UP ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: EMULT_DN ( MAX_USER_STREAMS )

!  Layer source terms

      REAL(FPK), INTENT(INOUT) :: WLAYER_PIST_UP &
        ( MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(INOUT) :: WLAYER_PIST_DN &
        ( MAX_LAYERS, MAX_USER_STREAMS )

!  Saved multiplier arrays

      REAL(FPK), INTENT(OUT) :: WSU_SUSAV ( MAX_STREAMS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: WSU_SDSAV ( MAX_STREAMS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: WSD_SUSAV ( MAX_STREAMS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: WSD_SDSAV ( MAX_STREAMS, MAX_USER_STREAMS )

!  Upwelling
!  ---------

      IF ( DO_UPWELLING ) THEN

        CALL BEAMMULT_UP_1 &
            ( N, N_USER_STREAMS, USER_STREAMS,          & ! Inputs
              TAYLOR_SMALL, TAYLOR_ORDER,               & ! Inputs
              STERM_LAYERMASK_UP, TRANS_USERM, DELTAUS, & ! Inputs
              FLIPPER, PICUTOFF, DELTRANS, ASOURCE,     & ! Inputs
              EMULT_UP )                                  ! Output

        CALL RAMAN_SOURCETERM_UP_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, TAYLOR_SMALL, TAYLOR_ORDER, & ! Inputs
            N, STERM_LAYERMASK_UP, NSTREAMS, N_USER_STREAMS,               & ! Inputs
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION,         & ! Inputs
            TRANS_USERM, DELTAUS, USER_STREAMS, ATERM_SAVE, BTERM_SAVE,    & ! Inputs
            EMULT_UP, QSOURCE_UPOS1, ASOURCE, DELTRANS, INITRANS,          & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,                      & ! Inputs (NO KTRANS)
            WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV )                           ! Outputs

      ENDIF

!  Downwelling
!  -----------

      IF ( DO_DNWELLING ) THEN

        CALL BEAMMULT_DN_1 &
            ( N, N_USER_STREAMS, USER_STREAMS,          & ! Inputs
              TAYLOR_SMALL, TAYLOR_ORDER,               & ! Inputs
              STERM_LAYERMASK_DN, TRANS_USERM, DELTAUS, & ! Inputs
              FLIPPER, PICUTOFF, DELTRANS, ASOURCE,     & ! Inputs
              EMULT_DN )                                  ! Output

        CALL RAMAN_SOURCETERM_DN_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, TAYLOR_SMALL, TAYLOR_ORDER, & ! Inputs
            N, STERM_LAYERMASK_DN, NSTREAMS, N_USER_STREAMS,               & ! Inputs
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION,         & ! Inputs
            TRANS_USERM, DELTAUS, USER_STREAMS, ATERM_SAVE, BTERM_SAVE,    & ! Inputs
            EMULT_DN, QSOURCE_UNEG1, ASOURCE, DELTRANS, INITRANS,          & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,                      & ! Inputs (NO KTRANS)
            WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )                           ! Outputs

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE POSTPROCESSING_1

!

      SUBROUTINE HOMOGMULT_1 &
          ( NSTREAMS, N, N_USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER, & ! Inputs
            UTRANS, USER_STREAMS, KEIGEN, KTRANS, DELTAUS,           & ! Inputs
            HMULT_1, HMULT_2, ZETA_P, ZETA_M )                         ! Outputs

!  Whole-layer INTEGRATED homogeneous solution MULTIPLIERS
!      upwelling and downwelling for the complete atmosphere

!   @@@ RobFix 5/5/11. Small numbers analysis:
!                      added FGSMALL, DELTAUS to Argument list.

!  LRRS Version 2.5. Rob Fix 9/10/15.
!    - Added Taylor_order variable to argument list. FGSMALL = TAYLOR_SMALL
!    - Use Taylor series subroutine and module, replaces older code.

!  Use Module of dimensions and numbers
!  ------------------------------------

      USE LRRS_PARS_m    , Only : FPK, ONE, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS

!  Taylor series module
!  --------------------

      USE lrrs_Taylor_m, Only : Taylor_series_1

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  inputs

      INTEGER  , INTENT(IN) :: N_USER_STREAMS, N, NSTREAMS

!  Taylor-series control

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  optical thickness for use with the small number expansion
!   @@@ RobFix 5/5/11. Small numbers analysis added

      REAL(FPK), INTENT(IN) :: DELTAUS

!  User angles

      REAL(FPK), INTENT(IN) :: UTRANS ( MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  Discrete ordinates

      REAL(FPK), INTENT(IN) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  output

      REAL(FPK), INTENT(INOUT) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(INOUT) :: &
            ZETA_P ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            ZETA_M ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Local variables
!  ---------------

      INTEGER   :: UM, AA
      REAL(FPK) :: UDEL, ZDEL, ZUDEL, SM, THETA_2, THETA_1
      REAL(FPK) :: RHO_P, RHO_M, SD

!  Get the multipliers

      DO UM = 1, N_USER_STREAMS
        UDEL = UTRANS(UM,N)
        SM  = ONE / USER_STREAMS(UM)
        DO AA = 1, NSTREAMS
          RHO_P = SM + KEIGEN(AA,N)
          RHO_M = SM - KEIGEN(AA,N)
          ZETA_P(AA,UM,N) = ONE / RHO_P
          ZETA_M(AA,UM,N) = ONE / RHO_M
          ZDEL    = KTRANS(AA,N)
          ZUDEL   = ZDEL * UDEL
          THETA_2 = ONE - ZUDEL
          THETA_1 = ZDEL - UDEL

!   @@@ RobFix 5/5/11. Small numbers analysis added
!    Rob Fix 9/10/15. Use Taylor series expansion, instead of LIMIT_GCFUNC
!      * Old Call  -->   CALL LIMIT_GCFUNC ( -RHO_M, DELTAUS, UDEL, SD )

          IF ( ABS(RHO_M) .LT. TAYLOR_SMALL ) THEN
            CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, RHO_M, DELTAUS, UDEL, ONE, SD )
          ELSE
            SD = THETA_1 * ZETA_M(AA,UM,N)
          ENDIF
          HMULT_1(AA,UM,N) = SM * SD
          HMULT_2(AA,UM,N) = SM * THETA_2 * ZETA_P(AA,UM,N)
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE HOMOGMULT_1

!

      SUBROUTINE BEAMMULT_UP_1 &
            ( N, N_USER_STREAMS, USER_STREAMS,          & ! Inputs
              TAYLOR_SMALL, TAYLOR_ORDER,               & ! Inputs
              STERM_LAYERMASK_UP, TRANS_USERM, DELTAUS, & ! Inputs
              FLIPPER, PICUTOFF, DELTRANS, ASOURCE,     & ! Inputs
              EMULT_UP )                                  ! Output

!  Source term primary multiplier. Upwelling.
!   - Small numbers analysis added 7 December 2006

!  LRRS Version 2.5. Rob Fix 9/10/15.
!    - Added Taylor_order variable to argument list. FGSMALL = TAYLOR_SMALL
!    - Use Taylor series subroutine and module, replaces older code.

!  Use Module of dimensions and numbers
!  ------------------------------------

      USE LRRS_PARS_m    , Only : FPK, ZERO, ONE, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS

!  Taylor series module
!  --------------------

      USE lrrs_Taylor_m, Only : Taylor_series_1

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  control integers

      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      INTEGER  , INTENT(IN) :: N

!  User streams, layer transmittances

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP
      REAL(FPK), INTENT(IN) :: TRANS_USERM  ( MAX_USER_STREAMS )

!  Taylor-series control

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  optical thickness for use with the small number expansion

      REAL(FPK), INTENT(IN) :: DELTAUS

!  Solution control

      LOGICAL  , INTENT(IN) :: FLIPPER
      INTEGER  , INTENT(IN) :: PICUTOFF
      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: ASOURCE

!  output = multiplier for one solution source

      REAL(FPK), INTENT(OUT) :: EMULT_UP ( MAX_USER_STREAMS )

!  local

      INTEGER   :: UM
      REAL(FPK) :: SM, UDEL, WUDEL, SIGMA_M, SIGMA_P, SD, SU

!    Rob Fix 9/10/15. Use Taylor series expansion, instead of LIMIT_GCFUNC
!      * Old Call  -->   CALL LIMIT_GCFUNC ( SIGMA_M, DELTAUS, UDEL, SD ).
!      * NEW CALL  -->   NOTE MINUS SIGN ON SIGMA_M

!mick fix 7/20/2016 - added ELSE section with initialization

      IF ( STERM_LAYERMASK_UP ) THEN
        IF ( N .GT. PICUTOFF ) THEN
          DO UM = 1, N_USER_STREAMS
            EMULT_UP(UM) = ZERO
          ENDDO
        ELSE
          DO UM = 1, N_USER_STREAMS
            SM   = ONE / USER_STREAMS(UM)
            UDEL = TRANS_USERM(UM)
            IF ( FLIPPER ) THEN
              SIGMA_M = ASOURCE - SM
              IF(ABS(SIGMA_M).LT.TAYLOR_SMALL)THEN
                CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, -SIGMA_M, DELTAUS, UDEL, ONE, SD )
              ELSE
                SD = ( UDEL - DELTRANS ) / SIGMA_M
              ENDIF
              EMULT_UP(UM) = SM * SD
            ELSE
              WUDEL = DELTRANS * UDEL
              SIGMA_P = ASOURCE + SM
              SU = SM * ( ONE - WUDEL ) / SIGMA_P
              EMULT_UP(UM) = SU
            ENDIF
          ENDDO
        ENDIF
      ELSE
        DO UM = 1, N_USER_STREAMS
          EMULT_UP(UM) = ZERO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BEAMMULT_UP_1

!

      SUBROUTINE BEAMMULT_DN_1 &
            ( N, N_USER_STREAMS, USER_STREAMS,          & ! Inputs
              TAYLOR_SMALL, TAYLOR_ORDER,               & ! Inputs
              STERM_LAYERMASK_DN, TRANS_USERM, DELTAUS, & ! Inputs
              FLIPPER, PICUTOFF, DELTRANS, ASOURCE,     & ! Inputs
              EMULT_DN )                                  ! Output

!  Source term primary multiplier. Downwelling.
!   - Small numbers analysis added 7 December 2006

!  LRRS Version 2.5. Rob Fix 9/10/15.
!    - Added Taylor_order variable to argument list. FGSMALL = TAYLOR_SMALL
!    - Use Taylor series subroutine and module, replaces older code.

!  Use Module of dimensions and numbers
!  ------------------------------------

      USE LRRS_PARS_m    , Only : FPK, ZERO, ONE, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS

!  Taylor series module
!  --------------------

      USE lrrs_Taylor_m, Only : Taylor_series_1

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Control integers

      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      INTEGER  , INTENT(IN) :: N

!  User streams, layer transmittances

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN
      REAL(FPK), INTENT(IN) :: TRANS_USERM  ( MAX_USER_STREAMS )

!  Taylor-series control

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  optical thickness for use with the small number expansion

      REAL(FPK), INTENT(IN) :: DELTAUS

!  Solution control

      LOGICAL  , INTENT(IN) :: FLIPPER
      INTEGER  , INTENT(IN) :: PICUTOFF
      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: ASOURCE

!  output = single multiplier for one solution source

      REAL(FPK), INTENT(OUT) :: EMULT_DN ( MAX_USER_STREAMS )

!  local variables

      INTEGER   :: UM
      REAL(FPK) :: SM, UDEL, WUDEL, SIGMA_M, SIGMA_P, SD, SU

!    Rob Fix 9/10/15. Use Taylor series expansion, instead of LIMIT_GCFUNC
!      * Old Call  -->   CALL LIMIT_GCFUNC ( SIGMA_M, DELTAUS, UDEL, SD )

!mick fix 7/20/2016 - added ELSE section with initialization

      IF ( STERM_LAYERMASK_DN ) THEN
        IF ( N .GT. PICUTOFF ) THEN
          DO UM = 1, N_USER_STREAMS
            EMULT_DN(UM) = ZERO
          ENDDO
        ELSE
          DO UM = 1, N_USER_STREAMS
            SM   = ONE / USER_STREAMS(UM)
            UDEL = TRANS_USERM(UM)
            IF ( FLIPPER ) THEN
              WUDEL = DELTRANS * UDEL
              SIGMA_P = ASOURCE + SM
              SU = SM * ( ONE - WUDEL ) / SIGMA_P
              EMULT_DN(UM) = SU
            ELSE
              SIGMA_M = ASOURCE - SM
              IF(ABS(SIGMA_M).LT.TAYLOR_SMALL)THEN
                CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, -SIGMA_M, DELTAUS, UDEL, ONE, SD )
              ELSE
                SD = ( UDEL - DELTRANS ) / SIGMA_M
              ENDIF
              EMULT_DN(UM) = SM * SD
            ENDIF
          ENDDO
        ENDIF
      ELSE
        DO UM = 1, N_USER_STREAMS
          EMULT_DN(UM) = ZERO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BEAMMULT_DN_1

!

      SUBROUTINE TOASOURCE &
          ( N_USER_STREAMS, & ! Input
            TOA_SOURCE )

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_USER_STREAMS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER  , INTENT(IN) :: N_USER_STREAMS

      REAL(FPK), INTENT(OUT) :: TOA_SOURCE(MAX_USER_STREAMS)

!  local variables

      INTEGER :: UM

!  initialise TOA source function

      DO UM = 1, N_USER_STREAMS
        TOA_SOURCE(UM) = ZERO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE TOASOURCE

!

      SUBROUTINE BOASOURCE &
          ( ELASTIC_CALL,                                                    & ! Debug
            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,      & ! Inputs
            NSTREAMS, NLAYERS, N_USER_STREAMS, FOURIER, Raman_IDX, Brdf_IDX, & ! Inputs
            QUAD_STRMWGHT, WLOWER, LCON_XVEC, MCON_XVEC, KTRANS,             & ! Inputs
            SURFACE_FACTOR, ALBEDOS_RANKED, USER_BRDF_F, USER_DIRECT_BEAM,   & ! Inputs
            IDOWNSURF, BOA_SOURCE, DIRECT_BOA_SOURCE )                         ! Outputs

!  Bottom of the atmosphere source term

!  This is LRRS Version 2.5. Main changes to this subroutine (from V2.3) are
!    (2) Use of Supplement-derived BRDF inputs and control, for BOASOURCE

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_2_STREAMS, MAX_STREAMS, MAX_USER_STREAMS, &
                            MAX_MOMENTS, MAX_LAYERS, MAX_POINTS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  debug flag

      LOGICAL  , INTENT(IN) :: ELASTIC_CALL

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

!  Variables for integrating the downwelling surface radiance

      REAL(FPK), INTENT(IN) :: &
          LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  surface multiplier

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Quadrature

      REAL(FPK), INTENT(IN)  :: QUAD_STRMWGHT ( MAX_STREAMS )

!  albedo

      REAL(fpk), intent(in)  :: ALBEDOS_RANKED ( MAX_POINTS )

!  incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: USER_BRDF_F ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Direct beam contribution

      REAL(FPK), INTENT(IN) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS )

!  output variables, direct and diffuse surface fields.

      REAL(FPK), INTENT(OUT) :: IDOWNSURF         ( MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER   :: N, I, UM, AA
      REAL(FPK) :: PAR, HOM, REFLEC

!  initialise BOA source function

      DO UM = 1, N_USER_STREAMS
        BOA_SOURCE(UM)        = ZERO
        DIRECT_BOA_SOURCE(UM) = ZERO
      ENDDO

!  Last layer only.

      N = NLAYERS

!  start albedo clause

      IF ( DO_INCLUDE_SURFACE ) THEN

!  Get downward intensity at computational angles (beam and homog)
!    And develop reflectance integrand  a(j).x(j).I(-j) (stored in commo

        DO I = 1, NSTREAMS
          PAR = WLOWER(I,N)
          HOM = ZERO
          DO AA = 1, NSTREAMS
            HOM = HOM + LCON_XVEC(I,AA,N) * KTRANS(AA,N) + MCON_XVEC(I,AA,N)
          ENDDO
          IDOWNSURF(I) = QUAD_STRMWGHT(I) * ( PAR + HOM )
        ENDDO

!  Compute reflectance
!  reflected multiple scatter intensity at user defined-angles

        IF ( DO_BRDF_SURFACE ) THEN
          DO UM = 1, N_USER_STREAMS
            REFLEC = DOT_PRODUCT(IDOWNSURF(1:NSTREAMS),USER_BRDF_F(FOURIER,UM,1:NSTREAMS,Brdf_IDX)) &
                   * SURFACE_FACTOR
            BOA_SOURCE(UM) = REFLEC
          ENDDO
        ELSE 
          IF ( FOURIER .EQ. 0 ) THEN
            REFLEC = SUM(IDOWNSURF(1:NSTREAMS)) * ALBEDOS_RANKED(Raman_idx) * SURFACE_FACTOR
            DO UM = 1, N_USER_STREAMS
              BOA_SOURCE(UM) = REFLEC
            ENDDO
          ENDIF
        ENDIF

!  Add direct beam if flagged

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          DO UM = 1, N_USER_STREAMS
            DIRECT_BOA_SOURCE(UM) = USER_DIRECT_BEAM(UM)
          ENDDO
        ENDIF

!  End surface clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BOASOURCE

!

      SUBROUTINE RAMAN_SOURCETERM_UP_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, TAYLOR_SMALL, TAYLOR_ORDER, & ! Inputs
            N, STERM_LAYERMASK_UP, NSTREAMS, N_USER_STREAMS,               & ! Inputs
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION,         & ! Inputs
            TRANS_USERM, DELTAUS, USER_STREAMS, ATERM_SAVE, BTERM_SAVE,    & ! Inputs
            EMULT_UP, QSOURCE_UPOS1, ASOURCE, DELTRANS, INITRANS,          & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,                      & ! Inputs (no KTRANS)
            WLAYER_PIST_UP, SUSAV, SDSAV )                                   ! Outputs

!  LRRS Version 2.5. Rob Fix 9/10/15.
!    - Added Taylor_order variable to argument list. FGSMALL = TAYLOR_SMALL
!    - Use Taylor series subroutine and module, replaces older code.

!  Use Module of dimensions and numbers
!  ------------------------------------

      USE LRRS_PARS_m    , Only : FPK, ZERO, ONE, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS

!  Taylor series module
!  --------------------

      USE lrrs_Taylor_m, Only : Taylor_series_2

      IMPLICIT NONE

!  input
!  -----

!  Solution number

      INTEGER  , INTENT(IN) :: SOLUTION

!  number of streams and solutions

      INTEGER  , INTENT(IN) :: NSTREAMS, N_USER_STREAMS

!  Number of layer, and layer control

      INTEGER  , INTENT(IN) :: N
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP

!  Taylor-series control

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  Fourier component

      INTEGER  , INTENT(IN) :: FOURIER

!  Flags for multiple scattering and single scatter correction to elasti

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT
      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION

!  Control for the particular solutions

      LOGICAL  , INTENT(IN) :: CHANGE_EXPONENT
      LOGICAL  , INTENT(IN) :: FLIPPER
      INTEGER  , INTENT(IN) :: PICUTOFF

!  Optical thickness value

      REAL(FPK), INTENT(IN) :: DELTAUS

!  stream directions

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  Layer transmittances in stream directions

      REAL(FPK), INTENT(IN) :: TRANS_USERM  ( MAX_USER_STREAMS )

!  Particular solution attributes

      REAL(FPK), INTENT(IN) :: ASOURCE
      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: INITRANS
      REAL(FPK), INTENT(IN) :: EMULT_UP      ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: QSOURCE_UPOS1 ( MAX_USER_STREAMS )

!  Saved ATERM and BTERM contributions to the Green's function.

      REAL(FPK), INTENT(IN) :: ATERM_SAVE  ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: BTERM_SAVE  ( MAX_STREAMS )

!  Eigensolutions at user-defined stream angles

      REAL(FPK), INTENT(IN) :: &
        U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
        U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Eigenvalues

      REAL(FPK), INTENT(IN) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )
 !     REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS ) ! removed

!  Whole layer multipliers, homogeneous solutions

      REAL(FPK), INTENT(IN) :: HMULT_1 &
        ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: HMULT_2 &
        ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  output
!  ------

!  Layer source terms

      REAL(FPK), INTENT(INOUT) :: WLAYER_PIST_UP &
        ( MAX_LAYERS, MAX_USER_STREAMS )

!  Saved multiplier arrays

      REAL(FPK), INTENT(OUT) :: SUSAV ( MAX_STREAMS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: SDSAV ( MAX_STREAMS, MAX_USER_STREAMS )

!  local variables
!  ---------------

!  help variables

      INTEGER   :: AA, UM, M
      REAL(FPK) :: PAR, RHOM, GP, PAR1, MULT, DMU, UDEL, KVAL, YFAC, FAC1, FAC2

!  Local multiplier arrays

      REAL(FPK) :: SU ( MAX_STREAMS ), SD ( MAX_STREAMS )

!  Version 2.3 upgrades - 
!    * Single RRS constribution: Greens function multipliers
!    * New code programmed with small number analysis and FLIPPER
!    * 27 September 2006, 29 November 2006. R. Spurr, RT SOLUTIONS.

!  Version 2.5 upgrades - 
!    Rob Fix 9/10/15. Use Taylor series expansions, instead of LIMIT_GCFUNC calls
!      * Old Call  -->  call limit_ppgmult_wd ( rhom, dmu, one, kval, deltaus, ktrans(aa,n), udel, susav(aa,um))
!      * Old Call  -->  call limit_ppgmult_wd ( rhom, dmu, one, kval, deltaus, ktrans(aa,n), udel, sdsav(aa,um))

!  Initialize (only if SOLUTION = 1)

      IF ( SOLUTION .EQ. 1 ) THEN
        DO UM = 1, N_USER_STREAMS
          WLAYER_PIST_UP(N,UM) = ZERO
        ENDDO
      ENDIF

!  Fourier (for debug only)

      M = FOURIER

      IF ( STERM_LAYERMASK_UP ) THEN

!  Start loop over user angles

        DO UM = 1, N_USER_STREAMS

!  Local quantities

          DMU    = ONE/USER_STREAMS(UM)
          UDEL   = TRANS_USERM(UM)

!  Only require a new multiplier when the exponent changes

          IF ( CHANGE_EXPONENT ) THEN

!  Set directional multipliers

            MULT = EMULT_UP(UM)

!  Zero the Greens function multipliers if no solution

            IF ( N .GT. PICUTOFF ) THEN

              DO AA = 1, NSTREAMS
                SDSAV(AA,UM) = ZERO
                SUSAV(AA,UM) = ZERO
              ENDDO

            ELSE

              IF ( FLIPPER ) THEN
                DO AA = 1, NSTREAMS
                  KVAL = KEIGEN(AA,N)

                  !mick: G+- (flip case)
                  RHOM = ASOURCE - KVAL
                  IF ( ABS(RHOM) .LT. TAYLOR_SMALL )THEN
                    FAC1 = TRANS_USERM(UM) ; YFAC  = ASOURCE - DMU ; FAC2 = DELTRANS
                    CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, TAYLOR_SMALL, RHOM, YFAC, DELTAUS, FAC1, FAC2, DMU, SUSAV(AA,UM) )
                  ELSE
                    SUSAV(AA,UM)  = ( HMULT_1(AA,UM,N) - MULT ) / RHOM
                  ENDIF

                  !mick: G++ (flip case)
                  GP = ONE / ( ASOURCE + KVAL )
                  SDSAV(AA,UM)  = GP*(MULT-HMULT_2(AA,UM,N)*DELTRANS)
                ENDDO
              ELSE
                DO AA = 1, NSTREAMS
                  KVAL = KEIGEN(AA,N)

                  !mick: G+-
                  RHOM = ASOURCE - KVAL
                  IF ( ABS(RHOM) .LT. TAYLOR_SMALL )THEN
                    FAC1 = ONE ; YFAC  = ASOURCE + DMU ; FAC2 = DELTRANS * TRANS_USERM(UM)
                    CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, TAYLOR_SMALL, RHOM, YFAC, DELTAUS, FAC1, FAC2, DMU, SDSAV(AA,UM) )
                  ELSE
                    SDSAV(AA,UM)  = ( HMULT_2(AA,UM,N) - MULT ) / RHOM
                  ENDIF

                  !mick: G++
                  GP = ONE / ( ASOURCE + KVAL )
                  SUSAV(AA,UM) = GP*(MULT-HMULT_1(AA,UM,N)*DELTRANS)
                ENDDO
              ENDIF

            ENDIF

!  Complete existence clause

          ENDIF

!  Add the terms independent of optical depth

          DO AA = 1, NSTREAMS
            SD(AA) = SDSAV(AA,UM) * ATERM_SAVE(AA)
            SU(AA) = SUSAV(AA,UM) * BTERM_SAVE(AA)
          ENDDO

!  Layer source computation
!  ========================

!  Summation over all streams

          PAR = ZERO
          DO AA = 1, NSTREAMS
            PAR = PAR + U_XPOS(UM,AA,N)*SD(AA) + U_XNEG(UM,AA,N)*SU(AA)
          ENDDO

!  Upgrade contribution

          WLAYER_PIST_UP(N,UM) = WLAYER_PIST_UP(N,UM) + PAR*INITRANS

!  End main user stream loop

        ENDDO

!  Add single scatter contribution if flagged.
!  If you are doing an SS correction, avoid this part
!  .. Otherwise, Full radiance mode, add single scatter part

        IF ( .NOT. DO_MSMODE_LIDORT ) THEN
          IF ( DO_SSCORRECTION ) GO TO 5555
          DO UM = 1, N_USER_STREAMS
            PAR1 = EMULT_UP(UM) * QSOURCE_UPOS1(UM)
            WLAYER_PIST_UP(N,UM) = WLAYER_PIST_UP(N,UM) + PAR1*INITRANS
          ENDDO
 5555     CONTINUE
        ENDIF

!  End layer

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE RAMAN_SOURCETERM_UP_1

!

      SUBROUTINE RAMAN_SOURCETERM_DN_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, TAYLOR_SMALL, TAYLOR_ORDER, & ! Inputs
            N, STERM_LAYERMASK_DN, NSTREAMS, N_USER_STREAMS,               & ! Inputs
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION,         & ! Inputs
            TRANS_USERM, DELTAUS, USER_STREAMS, ATERM_SAVE, BTERM_SAVE,    & ! Inputs
            EMULT_DN, QSOURCE_UNEG1, ASOURCE, DELTRANS, INITRANS,          & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,                      & ! Inputs (no KTRANS)
            WLAYER_PIST_DN, SUSAV, SDSAV )                                   ! Outputs

!  LRRS Version 2.5. Rob Fix 9/10/15.
!    - Added Taylor_order variable to argument list. FGSMALL = TAYLOR_SMALL
!    - Use Taylor series subroutine and module, replaces older code.

!  Use Module of dimensions and numbers
!  ------------------------------------

      USE LRRS_PARS_m    , Only : FPK, ZERO, ONE, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS

!  Taylor series module
!  --------------------

      USE lrrs_Taylor_m, Only : Taylor_series_2

      IMPLICIT NONE

!  arguments
!  ---------

!  Solution number

      INTEGER  , INTENT(IN) :: SOLUTION

!  number of streams and solutions

      INTEGER  , INTENT(IN) :: NSTREAMS, N_USER_STREAMS

!  Number of layer, and layer control

      INTEGER  , INTENT(IN) :: N
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN

!  Taylor-series control

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  Fourier component

      INTEGER  , INTENT(IN) :: FOURIER

!  Flags for multiple scattering and single scatter correction to elasti

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT
      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION

!  Control for the particular solutions

      LOGICAL  , INTENT(IN) :: CHANGE_EXPONENT
      LOGICAL  , INTENT(IN) :: FLIPPER
      INTEGER  , INTENT(IN) :: PICUTOFF

!  Optical thickness value

      REAL(FPK), INTENT(IN) :: DELTAUS

!  stream directions

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  Layer transmittances in stream directions

      REAL(FPK), INTENT(IN) :: TRANS_USERM  ( MAX_USER_STREAMS )

!  Particular solution attributes

      REAL(FPK), INTENT(IN) :: ASOURCE
      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: INITRANS
      REAL(FPK), INTENT(IN) :: EMULT_DN      ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: QSOURCE_UNEG1 ( MAX_USER_STREAMS )

!  Saved ATERM and BTERM contributions to the Green's function.

      REAL(FPK), INTENT(IN) :: ATERM_SAVE  ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: BTERM_SAVE  ( MAX_STREAMS )

!  Eigensolutions at user-defined stream angles

      REAL(FPK), INTENT(IN) :: &
        U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
        U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Eigenvalues and transmittance factors.

!      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )

!  Whole layer multipliers, homogeneous solutions

      REAL(FPK), INTENT(IN) :: HMULT_1 &
        ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: HMULT_2 &
        ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  output
!  ------

!  Layer source terms

      REAL(FPK), INTENT(INOUT) :: WLAYER_PIST_DN &
        ( MAX_LAYERS ,MAX_USER_STREAMS )

!  Saved multiplier arrays

      REAL(FPK), INTENT(OUT) :: SUSAV ( MAX_STREAMS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: SDSAV ( MAX_STREAMS, MAX_USER_STREAMS )

!  local variables
!  ---------------

!  help variables

      INTEGER   :: AA, UM, M
      REAL(FPK) :: PAR, RHOM, GP, MULT, DMU, UDEL, KVAL, YFAC, FAC1, FAC2

!  Local multiplier arrays

      REAL(FPK) :: SU ( MAX_STREAMS ), SD ( MAX_STREAMS )

!  Version 2.3 upgrades - 
!    * Single RRS constribution: Greens function multipliers
!    * New code programmed with small number analysis and FLIPPER
!    * 27 September 2006, 29 November 2006. R. Spurr, RT SOLUTIONS.

!  Version 2.5 upgrades - 
!    Rob Fix 9/10/15. Use Taylor series expansions, instead of LIMIT_GCFUNC calls
!      * Old Call  -->  call limit_ppgmult_wu ( rhom, dmu, one, kval, deltaus, ktrans(aa,n), udel, susav(aa,um))
!      * Old Call  -->  call limit_ppgmult_wd ( rhom, dmu, one, kval, deltaus, ktrans(aa,n), udel, sdsav(aa,um))

!  Initialize (only if SOLUTION = 1)

      IF ( SOLUTION .EQ. 1 ) THEN
        DO UM = 1, N_USER_STREAMS
          WLAYER_PIST_DN(N,UM) = ZERO
        ENDDO
      ENDIF

!  Fourier (for debug only)

      M = FOURIER

      IF ( STERM_LAYERMASK_DN ) THEN

!  Start loop over user angles

        DO UM = 1, N_USER_STREAMS

!  local quantities

          DMU    = ONE/USER_STREAMS(UM)
          UDEL   = TRANS_USERM(UM)

!  Only require a new multiplier when the exponent changes

          IF ( CHANGE_EXPONENT ) THEN

!  set directional multipliers

            MULT   = EMULT_DN(UM)

!  Zero the Greens function multipliers if no solution

            IF ( N .GT. PICUTOFF ) THEN

              DO AA = 1, NSTREAMS
                SDSAV(AA,UM) = ZERO
                SUSAV(AA,UM) = ZERO
              ENDDO

!  Compute multipliers.

            ELSE

              IF ( FLIPPER ) THEN
                DO AA = 1, NSTREAMS
                  KVAL =  KEIGEN(AA,N)

                  !mick: G-- (flip case)
                  RHOM = ASOURCE - KVAL
                  IF ( ABS(RHOM) .LT. TAYLOR_SMALL ) THEN
                    FAC1 = ONE ; YFAC  = ASOURCE + DMU ; FAC2 = DELTRANS * TRANS_USERM(UM)
                    CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, TAYLOR_SMALL, RHOM, YFAC, DELTAUS, FAC1, FAC2, DMU, SUSAV(AA,UM) )
                  ELSE
                    SUSAV(AA,UM)  = ( HMULT_2(AA,UM,N) - MULT ) / RHOM
                  ENDIF

                  !mick: G-+ (flip case)
                  GP = ONE / ( ASOURCE + KVAL )
                  SDSAV(AA,UM) = GP*(MULT-HMULT_1(AA,UM,N)*DELTRANS)
                ENDDO
              ELSE
                DO AA = 1, NSTREAMS
                  KVAL =  KEIGEN(AA,N)

                  !mick: G--
                  RHOM = ASOURCE - KVAL
                  IF ( ABS(RHOM) .LT. TAYLOR_SMALL ) THEN
                    FAC1 = TRANS_USERM(UM) ; YFAC  = ASOURCE - DMU ; FAC2 = DELTRANS
                    CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, TAYLOR_SMALL, RHOM, YFAC, DELTAUS, FAC1, FAC2, DMU, SDSAV(AA,UM) )
                  ELSE
                    SDSAV(AA,UM)  = ( HMULT_1(AA,UM,N) - MULT ) / RHOM
                  ENDIF

                  !mick: G-+
                  GP = ONE / ( ASOURCE + KVAL )
                  SUSAV(AA,UM) = GP*(MULT-HMULT_2(AA,UM,N)*DELTRANS)
                ENDDO
              ENDIF

            ENDIF

!  Complete existence clause

          ENDIF

!  set up multipliers

          DO AA = 1, NSTREAMS
            SD(AA) = SDSAV(AA,UM) * ATERM_SAVE(AA)
            SU(AA) = SUSAV(AA,UM) * BTERM_SAVE(AA)
          ENDDO

!  Layer source computation
!  ========================

!  Summation over all streams

          PAR = ZERO
          DO AA = 1, NSTREAMS
            PAR = PAR + U_XNEG(UM,AA,N)*SD(AA) + U_XPOS(UM,AA,N)*SU(AA)
          ENDDO

!  Upgrade contribution

          WLAYER_PIST_DN(N,UM) = WLAYER_PIST_DN(N,UM) + PAR*INITRANS

!  End main loop over user streams

        ENDDO

!  Add single scatter contribution if flagged.
!  If you are doing an SS correction, avoid this part
!  .. Otherwise, Full radiance mode, add single scatter part

        IF ( .NOT. DO_MSMODE_LIDORT ) THEN
          IF ( DO_SSCORRECTION ) GO TO 5555
          DO UM = 1, N_USER_STREAMS
            PAR = EMULT_DN(UM) * QSOURCE_UNEG1(UM)
            WLAYER_PIST_DN(N,UM) = WLAYER_PIST_DN(N,UM) + PAR*INITRANS
          ENDDO
 5555     CONTINUE
        ENDIF

!  Finish layer

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE RAMAN_SOURCETERM_DN_1

! End module

      END MODULE lrrs_postprocessing_1_m

