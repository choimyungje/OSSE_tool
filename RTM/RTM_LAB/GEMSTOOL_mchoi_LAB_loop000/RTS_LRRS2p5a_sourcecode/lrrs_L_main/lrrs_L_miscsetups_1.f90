
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
! #   SUBROUTINES in this module (5/12/17 for 2p5a)             #
! #   ---------------------------------------------             #
! #                                                             #
! #   Setup SUBROUTINES (formerly in lrrs_L_rtsolutions_1.f90)  #
! #           (1 = Whole-layer, 2 = Part-layer)                 #
! #                                                             #
! #     LC  denotes column   linearization                      #
! #     LCH denotes column   linearization with heights         #
! #     LP  denotes profile  linearization                      #
! #                                                             #
! #             LC_ATTENUATION_SETUP_1                          #
! #             LP_ATTENUATION_SETUP_1                          #
! #                                                             #
! #    Special routines for the LCH linearization               #
! #          (Added, 19 April 2010)                             #
! #                                                             #
! #             CHAPMAN_FUNCTION_PLUS                           #
! #             LCH_ATTENUATION_SETUP_1                         #
! #                                                             #
! #    New Setup routine for the DO_SSCORR_NADIR option         #
! #                                                             #
! #             LPC_BEAM_MULTIPLIERS_1                          #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are

!    1. Introduction of supplement-created BRDF   expressions in DBEAM_SETUP routines
!    2. Introduction of supplement-created SLEAVE expressions in DBEAM_SETUP routines
!    3. Use of new Taylor series subroutine/module in LPC_GREENFUNC_SOLUTION_1
!    4. Bookkeeping improvements (use of "Only", clearer I/O specifications)

!   Items (1) and (2) replace the somewhat ad-hoc SL and BRDF treatments from earlier versions
!     - The new replacements are modeled after the LIDORT Version 3.7 code

!  Regular code was programmed 9/9/15 by R. Spurr. RT SOLUTIONS Inc.
!    Changes there marked by the trope "Rob Fix 9/9/15"

!  Linearized code was programmed 10/10/15 by R. Spurr. RT SOLUTIONS Inc.
!    Changes there marked by the trope "Rob Fix 10/10/15". THIS MODULE.

!   -- Rob mod 5/12/17 for 2p5a. This is a New Module with 5 subroutines
!        4 old subroutines formerly in lrrs_L_rtsolutions_1.f90
!        1 new subroutine  (L_BEAM_MULTIPLIERS_1) for the DO_SSCORR_NADIR option

      MODULE lrrs_L_miscsetups_1_m

      USE LRRS_PARS_m, Only : LDU

      PRIVATE
      PUBLIC :: LC_ATTENUATION_SETUP_1,     &
                LP_ATTENUATION_SETUP_1,     &
                LPC_BEAM_MULTIPLIERS_1,     &
                CHAPMAN_FUNCTION_PLUS,      &
                LCH_ATTENUATION_SETUP_1

      CONTAINS

!

      SUBROUTINE LPC_BEAM_MULTIPLIERS_1 &
        ( DO_UPWELLING, DO_DNWELLING, DO_COLUMN_WFS, DO_PROFILE_WFS,         & ! Inputs
          NPOINTS, NLAYERS, N_USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER,      & ! Inputs
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS, NKSTORAGE2, & ! Inputs
          STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, USER_STREAMS,              & ! Inputs
          DELTAU_VERT, SAVE_TRANS_USERM, L_DELTAU_VERT, L_SAVE_TRANS_USERM,  & ! Inputs
          BEAM_PICUTOFF, BEAM_AVSECANT, BEAM_DTRANS, BEAM_ITRANS,            & ! Inputs
          L_BEAM_ITRANS, L_BEAM_AVSECANT, L_BEAM_DTRANS,                     & ! Inputs
          SAVE_BEAMMULT_UP, L_SAVE_BEAMMULT_UP,                              & ! Output
          SAVE_BEAMMULT_DN, L_SAVE_BEAMMULT_DN )                               ! Output

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, MAX_LAYERS, MAX_LAYERS_NK, &
                              MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_POINTS

!  Dependency

      USE lrrs_postprocessing_1_m  , Only : BEAMMULT_UP_1,     BEAMMULT_DN_1
      USE lrrs_L_postprocessing_1_m, Only : LPC_BEAMMULT_UP_1, LPC_BEAMMULT_DN_1

!  Implicit none

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Flags

      LOGICAL  , INTENT(IN) :: DO_UPWELLING, DO_DNWELLING

!  Linearization control profile Jacobians

      LOGICAL  , INTENT(IN) :: do_profile_wfs
      LOGICAL  , INTENT(IN) :: layer_vary_flag   (max_layers)
      INTEGER  , INTENT(IN) :: layer_vary_number (max_layers)

!  Column linearization control inputs

      LOGICAL  , INTENT(IN) :: do_column_wfs
      INTEGER  , INTENT(IN) :: n_totalcolumn_wfs
      INTEGER  , INTENT(IN) :: NKSTORAGE2(MAX_LAYERS,0:MAX_LAYERS)

!  Numbers and bookkeeping
!  -----------------------

!  Control numbers

      INTEGER  , INTENT(IN) :: NPOINTS, NLAYERS, N_USER_STREAMS

!  layer control

      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP(MAX_LAYERS)
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN(MAX_LAYERS)

!  Taylor-series control. Updated, Version 2.5, 9/10/15

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  stream directions

      REAL(FPK), INTENT(IN) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: DELTAU_VERT   ( MAX_LAYERS, MAX_POINTS )

!  Solar beam transmittances, average secant factors

      INTEGER  , INTENT(IN) :: BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(IN) :: SAVE_TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Linearized inputs

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT   ( MAX_ATMOSWFS, MAX_LAYERS,    MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_BEAM_ITRANS   ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_BEAM_AVSECANT ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_BEAM_DTRANS   ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_SAVE_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Output multipliers
!  ------------------

      REAL(FPK), INTENT(OUT) :: SAVE_BEAMMULT_UP   ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: SAVE_BEAMMULT_DN   ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: L_SAVE_BEAMMULT_UP ( MAX_ATMOSWFS, MAX_USER_STREAMS, 0:MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_SAVE_BEAMMULT_DN ( MAX_ATMOSWFS, MAX_USER_STREAMS, 0:MAX_LAYERS_NK, MAX_POINTS )

!  Local variables
!  ---------------

      LOGICAL, parameter :: FLIPPER = .false.
      INTEGER   :: N, S, K, Q, NK, NU

 !  Linearization control numbers

      INTEGER   :: KS2(MAX_LAYERS),  KF2(MAX_LAYERS)
      INTEGER   :: KPARS(0:MAX_LAYERS)

!  Basic proxies

      REAL(FPK) :: DELTAUS, DELTRANS, ASOURCE, INITRANS
      REAL(FPK) :: TRANS_USERM ( MAX_USER_STREAMS )
      REAL(FPK) :: EMULT_UP    ( MAX_USER_STREAMS )
      REAL(FPK) :: EMULT_DN    ( MAX_USER_STREAMS )

!  Linearized proxies

      REAL(FPK) :: L_DELTAUS      ( MAX_ATMOSWFS )
      REAL(FPK) :: L_TRANS_USERM  ( MAX_ATMOSWFS, MAX_USER_STREAMS )
      REAL(FPK) :: L_ASOURCE      ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK) :: L_DELTRANS     ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK) :: L_INITRANS     ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK) :: L_EMULT_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS , MAX_USER_STREAMS )
      REAL(FPK) :: L_EMULT_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS , MAX_USER_STREAMS )

!  Initialize - No need
!      SAVE_BEAMMULT_UP = zero
!      SAVE_BEAMMULT_DN = zero

!  Monochromatic - 234 points (233 shifts + 1 calculation point)
!  Binned        - all points in outer buffer

!  Atmospheric profile linearization control

      IF ( DO_PROFILE_WFS ) THEN
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            KS2(N)   = 1 ; KF2(N)   = NLAYERS
            KPARS(N) = LAYER_VARY_NUMBER(N)
          ELSE
            KS2(N)   = 1 ; KF2(N)   = 0 ;  KPARS(N) = 0
          ENDIF
        ENDDO
      ENDIF
      IF ( DO_COLUMN_WFS ) THEN
        DO N = 1, NLAYERS
          KS2(N)   = 0 ; KF2(N)   = 0
          KPARS(0) = N_TOTALCOLUMN_WFS
        ENDDO
      ENDIF

!  Proxy 

      NU = N_USER_STREAMS

!  Upwelling
!  ---------

      IF ( DO_UPWELLING ) THEN
        DO S = 1, NPOINTS
          DO N = 1, NLAYERS

!  inputs

            DELTAUS  = DELTAU_VERT(N,S)
            DELTRANS = BEAM_DTRANS(N,S)
            ASOURCE  = BEAM_AVSECANT(N,S)
            INITRANS = BEAM_ITRANS(N,S)
            TRANS_USERM(1:NU) = SAVE_TRANS_USERM (1:NU,N,S)

!  Linearized inputs

            DO K = KS2(N), KF2(N)
              NK = NKSTORAGE2(N,K)
              DO Q = 1, KPARS(K)
                if ( K.eq.N .or. K.eq.0 ) then
                  L_TRANS_USERM(Q,1:NU) = L_SAVE_TRANS_USERM(Q,1:NU,N,S)
                  L_DELTAUS(Q)     = L_DELTAU_VERT(Q,N,S)
                endif
                L_ASOURCE(Q,K)  = L_BEAM_AVSECANT(Q,NK,S)
                L_INITRANS(Q,K) = L_BEAM_ITRANS(Q,NK,S)
                L_DELTRANS(Q,K) = L_BEAM_DTRANS(Q,NK,S)
              ENDDO
            ENDDO

!  Multiplier

            CALL BEAMMULT_UP_1 &
              ( N, N_USER_STREAMS, USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER, & ! Inputs
                STERM_LAYERMASK_UP(N), TRANS_USERM, DELTAUS,   & ! Inputs
                FLIPPER, BEAM_PICUTOFF(S), DELTRANS, ASOURCE,  & ! Inputs
                EMULT_UP )                                       ! Output

!  Linearized multiplier

            CALL LPC_BEAMMULT_UP_1 &
              ( N, N_USER_STREAMS, USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER,            & ! Inputs
                KS2(N), KF2(N), KPARS, STERM_LAYERMASK_UP(N), TRANS_USERM, L_TRANS_USERM, & ! Inputs
                DELTAUS, FLIPPER, BEAM_PICUTOFF(S), DELTRANS, ASOURCE,        & ! Inputs
                EMULT_UP, L_DELTAUS, L_DELTRANS, L_ASOURCE,    & ! Inputs
                L_EMULT_UP )                                     ! Output

!  output

            SAVE_BEAMMULT_UP(1:NU,N,S) = EMULT_UP(1:NU) * INITRANS
            DO K = KS2(N), KF2(N)
              NK = NKSTORAGE2(N,K)
              DO Q = 1, KPARS(K)
                L_SAVE_BEAMMULT_UP(Q,1:NU,NK,S) = EMULT_UP(1:NU) * L_INITRANS(Q,K) + L_EMULT_UP(Q,K,1:NU) * INITRANS
              ENDDO
            ENDDO

!  end layer and point loops

          ENDDO
        ENDDO
      ENDIF

!  Downwelling

      IF ( DO_DNWELLING ) THEN
        DO S = 1, NPOINTS
          DO N = 1, NLAYERS

!  inputs

            DELTAUS  = DELTAU_VERT(N,S)
            DELTRANS = BEAM_DTRANS(N,S)
            ASOURCE  = BEAM_AVSECANT(N,S)
            INITRANS = BEAM_ITRANS(N,S)
            TRANS_USERM(1:NU) = SAVE_TRANS_USERM (1:NU,N,S)

!  Linearized inputs

            DO K = KS2(N), KF2(N)
              NK = NKSTORAGE2(N,K)
              DO Q = 1, KPARS(K)
                if ( K.eq.N .or. K.eq.0 ) then
                  L_TRANS_USERM(Q,1:NU) = L_SAVE_TRANS_USERM(Q,1:NU,N,S)
                  L_DELTAUS(Q)     = L_DELTAU_VERT(Q,N,S)
                endif
                L_ASOURCE(Q,K)  = L_BEAM_AVSECANT(Q,NK,S)
                L_INITRANS(Q,K) = L_BEAM_ITRANS(Q,NK,S)
                L_DELTRANS(Q,K) = L_BEAM_DTRANS(Q,NK,S)
              ENDDO
            ENDDO

!  Multiplier

            CALL BEAMMULT_DN_1 &
              ( N, N_USER_STREAMS, USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER, & ! Inputs
                STERM_LAYERMASK_DN(N), TRANS_USERM, DELTAUS,   & ! Inputs
                FLIPPER, BEAM_PICUTOFF(S), DELTRANS, ASOURCE,  & ! Inputs
                EMULT_DN )                                       ! Output

!  Linearized multiplier

            CALL LPC_BEAMMULT_DN_1 &
              ( N, N_USER_STREAMS, USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER,            & ! Inputs
                KS2(N), KF2(N), KPARS, STERM_LAYERMASK_DN(N), TRANS_USERM, L_TRANS_USERM, & ! Inputs
                DELTAUS, FLIPPER, BEAM_PICUTOFF(S), DELTRANS, ASOURCE,        & ! Inputs
                EMULT_DN, L_DELTAUS, L_DELTRANS, L_ASOURCE,    & ! Inputs
                L_EMULT_DN )                                     ! Output

!  output

            SAVE_BEAMMULT_DN(1:NU,N,S) = EMULT_DN(1:NU) * INITRANS
            DO K = KS2(N), KF2(N)
              NK = NKSTORAGE2(N,K)
              DO Q = 1, KPARS(K)
                L_SAVE_BEAMMULT_DN(Q,1:NU,NK,S) = EMULT_DN(1:NU) * L_INITRANS(Q,K) + L_EMULT_DN(Q,K,1:NU) * INITRANS
              ENDDO
            ENDDO

!  end layer and point loops


          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LPC_BEAM_MULTIPLIERS_1

!

      SUBROUTINE LC_ATTENUATION_SETUP_1 &
          ( DO_PLANE_PARALLEL, NPOINTS, NLAYERS, NPARS,      & ! Inputs
            DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS,   & ! Inputs
            DELTAU_VERT_INPUT, CHAPMAN, BEAM_PICUTOFF,       & ! Inputs
            BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS,         & ! Inputs
            SAVE_TRANS_USERM, L_DELTAU_VERT_INPUT,           & ! Inputs
            L_BEAM_ITRANS, L_BEAM_AVSECANT,                  & ! Outputs 
            L_BEAM_ETRANS, L_BEAM_DTRANS, L_SAVE_TRANS_USERM ) ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m    , Only : FPK, ZERO, ONE, MAX_TAU_UPATH, MAX_USER_STREAMS, MAX_ATMOSWFS, &
                                  MAX_LAYERS, MAX_LAYERS_NK, MAX_POINTS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Control

      LOGICAL  , INTENT(IN) :: DO_PLANE_PARALLEL

!  Control integers

      INTEGER  , INTENT(IN) :: NPOINTS, NLAYERS

!  Linearization control

      INTEGER  , INTENT(IN) :: NPARS

!  Control flag and inputs for user streams

      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS
      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  Chapman factors

      REAL(FPK), INTENT(IN) :: CHAPMAN ( MAX_LAYERS, MAX_LAYERS )

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT   ( MAX_LAYERS, MAX_POINTS )

!  Solar beam transmittances, average secant factors

      INTEGER  , INTENT(IN) :: BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(IN) :: SAVE_TRANS_USERM &
            ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT_INPUT &
                        ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

!  Output arguments
!  ----------------

!  Linearized Solar beam transmittances, average secant factors

      REAL(FPK), INTENT(OUT) :: L_BEAM_ITRANS &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_BEAM_AVSECANT &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_BEAM_ETRANS &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_BEAM_DTRANS &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )

!  Linearized user stream transmittances, whole layers

      REAL(FPK), INTENT(OUT) :: L_SAVE_TRANS_USERM &
            ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Local variables
!  ---------------

      INTEGER   :: N, K, S, UM, Q
      REAL(FPK) :: SPHER, SUMR, SM, L_SPHER, T1, T2
      REAL(FPK) :: TAUSLANT(MAX_LAYERS)
      REAL(FPK) :: L_TAUSLANT(MAX_ATMOSWFS,MAX_LAYERS)

!  Beam attenuation
!  ----------------

!  Start the points loop

      DO S = 1, NPOINTS

!  Tau slant = Beam slant optical thickness

        DO N = 1, NLAYERS
          SUMR = ZERO
          DO K = 1, N
            SUMR = SUMR + DELTAU_VERT_INPUT(K,S) * CHAPMAN(N,K)
          ENDDO
          TAUSLANT(N) = SUMR
          DO Q = 1, NPARS
            SUMR = ZERO
            DO K = 1, N
              SUMR = SUMR + L_DELTAU_VERT_INPUT(Q,K,S) * CHAPMAN(N,K)
            ENDDO
            L_TAUSLANT(Q,N) = SUMR
          ENDDO
        ENDDO

!  First, linearize the initial transmittances

        L_BEAM_ITRANS(:,:,S) = ZERO
        DO N = 2, NLAYERS
          IF ( N.LE.BEAM_PICUTOFF(S) ) THEN
            DO Q = 1, NPARS
              L_BEAM_ITRANS(Q,N,S) = - L_TAUSLANT(Q,N-1) * BEAM_ITRANS(N,S)
            ENDDO
          ENDIF
        ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

        L_BEAM_AVSECANT(:,:,S) = ZERO
        IF ( .NOT. DO_PLANE_PARALLEL ) THEN
          DO N = 2, NLAYERS
            IF  (N.LE.BEAM_PICUTOFF(S) ) THEN
              DO Q = 1, NPARS
                T1 = L_TAUSLANT(Q,N) - L_TAUSLANT(Q,N-1)
                T2 = BEAM_AVSECANT(N,S) * L_DELTAU_VERT_INPUT(Q,N,S)
                L_BEAM_AVSECANT(Q,N,S) = ( T1 - T2 ) / DELTAU_VERT_INPUT(N,S)
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Other variables

        DO N = 1, NLAYERS
          DO Q = 1, NPARS
            L_BEAM_DTRANS(Q,N,S) = - BEAM_DTRANS(N,S) * &
               ( BEAM_AVSECANT(N,S)   * L_DELTAU_VERT_INPUT(Q,N,S) + &
               L_BEAM_AVSECANT(Q,N,S) *   DELTAU_VERT_INPUT(N,S) )
            L_BEAM_ETRANS(Q,N,S) = &
                L_BEAM_DTRANS(Q,N,S) *   BEAM_ITRANS(N,S) + &
                  BEAM_DTRANS(N,S)   * L_BEAM_ITRANS(Q,N,S)
          ENDDO
        ENDDO

!  End points loop

      ENDDO

!  If no user streams, then return

      IF ( .NOT. DO_USER_STREAMS  ) RETURN

!  whole layer transmittances along user streams
!  =============================================

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          SM = ONE / USER_STREAMS(UM)
          DO S = 1, NPOINTS
            DO N = 1, NLAYERS
              SPHER = SM * DELTAU_VERT_INPUT(N,S)
              IF ( SPHER.GT.MAX_TAU_UPATH ) THEN
                DO Q = 1, NPARS
                  L_SAVE_TRANS_USERM(Q,UM,N,S) = ZERO
                ENDDO
              ELSE
                DO Q = 1, NPARS
                  L_SPHER = SM * L_DELTAU_VERT_INPUT(Q,N,S)
                  L_SAVE_TRANS_USERM(Q,UM,N,S) = - L_SPHER * SAVE_TRANS_USERM(UM,N,S)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LC_ATTENUATION_SETUP_1

!

      SUBROUTINE LP_ATTENUATION_SETUP_1 &
          ( DO_PLANE_PARALLEL, NPOINTS, NLAYERS, NKSTORAGE,    & ! Inputs
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                & ! Inputs
            DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS,     & ! Inputs
            DELTAU_VERT_INPUT, CHAPMAN, BEAM_PICUTOFF,         & ! Inputs
            BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS,           & ! Inputs
            SAVE_TRANS_USERM, L_DELTAU_VERT_INPUT,             & ! Inputs
            L_BEAM_ITRANS, L_BEAM_AVSECANT,                    & ! Outputs
            L_BEAM_ETRANS, L_BEAM_DTRANS, L_SAVE_TRANS_USERM )   ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m   , Only : FPK, ZERO, ONE, MAX_TAU_UPATH, MAX_USER_STREAMS, MAX_ATMOSWFS, &
                                 MAX_LAYERS, MAX_LAYERS_NK, MAX_POINTS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Control

      LOGICAL  , INTENT(IN) :: DO_PLANE_PARALLEL

!  Control integers

      INTEGER  , INTENT(IN) :: NPOINTS, NLAYERS

!  Linearization control

      LOGICAL  , INTENT(IN) :: LAYER_VARY_FLAG   ( MAX_LAYERS )
      INTEGER  , INTENT(IN) :: LAYER_VARY_NUMBER ( MAX_LAYERS )
      INTEGER  , INTENT(IN) :: NKSTORAGE ( MAX_LAYERS, 0:MAX_LAYERS )

!  Control flag and inputs for user streams

      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS
      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  Chapman factors

      REAL(FPK), INTENT(IN) :: CHAPMAN ( MAX_LAYERS, MAX_LAYERS )

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT   ( MAX_LAYERS, MAX_POINTS )

!  Solar beam transmittances, average secant factors

      INTEGER  , INTENT(IN) :: BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(IN) :: SAVE_TRANS_USERM &
            ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT_INPUT &
                        ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

!  Output arguments
!  ----------------

!  Linearized Solar beam transmittances, average secant factors

      REAL(FPK), INTENT(OUT) :: L_BEAM_ITRANS &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_BEAM_AVSECANT &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_BEAM_ETRANS &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_BEAM_DTRANS &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )

!  Linearized user stream transmittances, whole layers

      REAL(FPK), INTENT(OUT) :: L_SAVE_TRANS_USERM &
            ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Local variables
!  ---------------

      INTEGER   :: N, K, NK, NK1, NN, S, UM, Q
      REAL(FPK) :: SPHER, SUMR, SM, L_SPHER
      REAL(FPK) :: WDEL, VAR, FAC, DELT

!  Old code
!      REAL(FPK) :: TAUSLANT(0:MAX_LAYERS)
!      REAL(FPK) :: L_TAUSLANT(MAX_ATMOSWFS,0:MAX_LAYERS)

      REAL(FPK) :: TAUSLANT(MAX_LAYERS_NK)
      REAL(FPK) :: L_TAUSLANT(MAX_ATMOSWFS,MAX_LAYERS_NK)

!  Beam attenuation
!  ----------------

!  Start the points loop

      DO S = 1, NPOINTS

!  Tau slant = Beam slant optical thickness

        DO N = 1, NLAYERS
          SUMR = ZERO
          DO K = 1, N
            SUMR = SUMR + DELTAU_VERT_INPUT(K,S) * CHAPMAN(N,K)
          ENDDO
          TAUSLANT(N) = SUMR
        ENDDO

!  Save linearization of slant optical depth

        DO N = 1, NLAYERS
          DO K = 1, N
            NK = NKSTORAGE(N,K)
            IF ( LAYER_VARY_FLAG(K) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_TAUSLANT(Q,NK) = L_DELTAU_VERT_INPUT(Q,K,S) * CHAPMAN(N,K)
              ENDDO
            ELSE
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_TAUSLANT(Q,NK) = ZERO
              ENDDO
            ENDIF
          ENDDO
        ENDDO

!  First, Get the linearized initial transmittances
!  ------------------------------------------------

!  This is modeled after the LIDORT code in lidort_l_miscsetups.f
!    Are we watching for Logartihmic derivative here ??????????

        L_BEAM_ITRANS(:,:,S) = ZERO
        DO N = 1, NLAYERS
          NN = NKSTORAGE(N,N)
          IF ( N.LE.BEAM_PICUTOFF(S) ) THEN
            IF ( N .GT. 1 ) THEN
              DO K = 1, N-1
                NK  = NKSTORAGE(N,K)
                NK1 = NKSTORAGE(N-1,K)
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  L_BEAM_ITRANS(Q,NK,S) = &
                    - L_TAUSLANT(Q,NK1) * BEAM_ITRANS(N,S)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

        L_BEAM_AVSECANT(:,:,S) = ZERO
        IF ( .NOT. DO_PLANE_PARALLEL ) THEN
          DO N = 1, NLAYERS
            NN = NKSTORAGE(N,N)
            IF ( N .GT. 1 ) THEN
              IF  (N.LE.BEAM_PICUTOFF(S) ) THEN
                DELT  = DELTAU_VERT_INPUT(N,S)
                FAC   = ( CHAPMAN(N,N) - BEAM_AVSECANT(N,S) ) / DELT
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_BEAM_AVSECANT(Q,NN,S) = &
                        L_DELTAU_VERT_INPUT(Q,N,S) * FAC
                ENDDO
                DO K = 1, N-1
                  NK  = NKSTORAGE(N,K)
                  FAC   = ( CHAPMAN(N,K) - CHAPMAN(N-1,K) ) / DELT
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    L_BEAM_AVSECANT(Q,NK,S) = &
                        L_DELTAU_VERT_INPUT(Q,K,S) * FAC
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF

!  Layer beam Transmittances
!  -------------------------

!  Start layer loop

        DO N = 1, NLAYERS

!  Help variables

         NN = NKSTORAGE(N,N)
         WDEL  = BEAM_DTRANS(N,S)
         VAR   = - WDEL * DELTAU_VERT_INPUT(N,S)
         FAC   = - WDEL * BEAM_AVSECANT(N,S)

!  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_BEAM_DTRANS(Q,NN,S) = FAC * L_DELTAU_VERT_INPUT(Q,N,S)
            ENDDO
          ELSE
            IF  ( N.LE.BEAM_PICUTOFF(S) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                L_BEAM_DTRANS(Q,NN,S) = L_DELTAU_VERT_INPUT(Q,N,S) * FAC &
                                      + L_BEAM_AVSECANT(Q,NN,S)    * VAR
              ENDDO
              DO K = 1, N-1
                NK = NKSTORAGE(N,K)
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  L_BEAM_DTRANS(Q,NK,S) = L_BEAM_AVSECANT(Q,NK,S) * VAR
                ENDDO
              ENDDO
            ELSE
              DO K = 1, N
                NK = NKSTORAGE(N,K)
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  L_BEAM_DTRANS(Q,NK,S) = ZERO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.BEAM_PICUTOFF(S) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_BEAM_DTRANS(Q,NN,S) = L_DELTAU_VERT_INPUT(Q,N,S) * FAC
            ENDDO
            DO K = 1, N-1
              NK = NKSTORAGE(N,K)
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_BEAM_DTRANS(Q,NK,S) = ZERO
              ENDDO
            ENDDO
          ELSE
            DO K = 1, N
              NK = NKSTORAGE(N,K)
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_BEAM_DTRANS(Q,NK,S) = ZERO
              ENDDO
            ENDDO
          ENDIF

         ENDIF

!  L_BEAM_ETRANS, by Chain rule

         DO K = 1, N
           NK = NKSTORAGE(N,K)
           DO Q = 1, LAYER_VARY_NUMBER(K)
             L_BEAM_ETRANS(Q,NK,S) = &
                L_BEAM_DTRANS(Q,NK,S) *   BEAM_ITRANS(N,S) &
              +   BEAM_DTRANS(N,S)    * L_BEAM_ITRANS(Q,NK,S)
           ENDDO
         ENDDO

!  end layer loop

        ENDDO

!  End points loop

      ENDDO

!  If no user streams, then return

      IF ( .NOT. DO_USER_STREAMS  ) RETURN

!  whole layer transmittances along user streams
!  =============================================

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          SM = ONE / USER_STREAMS(UM)
          DO S = 1, NPOINTS
            DO N = 1, NLAYERS
              SPHER = SM * DELTAU_VERT_INPUT(N,S)
              IF ( SPHER.GT.MAX_TAU_UPATH ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_SAVE_TRANS_USERM(Q,UM,N,S) = ZERO
                ENDDO
              ELSE
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_SPHER = SM * L_DELTAU_VERT_INPUT(Q,N,S)
                  L_SAVE_TRANS_USERM(Q,UM,N,S) = - L_SPHER * &
                                       SAVE_TRANS_USERM(UM,N,S)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LP_ATTENUATION_SETUP_1

!

      SUBROUTINE CHAPMAN_FUNCTION_PLUS &
           ( DO_PLANE_PARALLEL, NLAYERS, NPARS,                 & ! Inputs
             COS_SZA, EARTH_RADIUS, HEIGHT_GRID, L_HEIGHT_GRID, & ! Inputs
             CHAPMAN_FACTORS, L_CHAPMAN_FACTORS )                 ! Outputs

!  This is the Chapman function calculation of the slant path
!  geometrical factors, defined as the ratios between the slant path
!  distances to the corresponding vertical distances.

!  Also calculated are the derivatives of the Chapman factors WRT
!   variations in the height grid itself due to some property.

!  You must specify the Earth_radius and the height grid in order
!  to make this work.

!  This is a straightforward geometrical calculation, and is only
!  valid for a NON-REFRACTIVE atmosphere.

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m   , Only : fpk, ZERO, ONE, TWO, HALF, MAX_LAYERS, MAX_ATMOSWFS

      IMPLICIT NONE

!  Input arguemnts
!  ---------------

      LOGICAL  , INTENT(IN) :: DO_PLANE_PARALLEL
      INTEGER  , INTENT(IN) :: NLAYERS, NPARS
      REAL(FPK), INTENT(IN) :: COS_SZA
      REAL(FPK), INTENT(IN) :: EARTH_RADIUS
      REAL(FPK), INTENT(IN) :: HEIGHT_GRID(0:MAX_LAYERS)
      REAL(FPK), INTENT(IN) :: L_HEIGHT_GRID(MAX_ATMOSWFS,0:MAX_LAYERS)

!  output arguemnts
!  ----------------

      REAL(FPK), INTENT(OUT) :: CHAPMAN_FACTORS(MAX_LAYERS,MAX_LAYERS)
      REAL(FPK), INTENT(OUT) :: L_CHAPMAN_FACTORS(MAX_ATMOSWFS,MAX_LAYERS,MAX_LAYERS)

!  Local variables
!  ---------------

      INTEGER   :: N, M, Q
      REAL(FPK) :: GM_TOA, HELP1, HELP2, HELP3
      REAL(FPK) :: H(0:MAX_LAYERS), DELZ(MAX_LAYERS)
      REAL(FPK) :: STH, CTH, DELS, S1, S0

      REAL(FPK) :: L_HELP1(MAX_ATMOSWFS), L_HELP2(MAX_ATMOSWFS), L_HELP3
      REAL(FPK) :: L_H   (MAX_ATMOSWFS,0:MAX_LAYERS)
      REAL(FPK) :: L_DELZ(MAX_ATMOSWFS,MAX_LAYERS)
      REAL(FPK) :: L_STH(MAX_ATMOSWFS), L_CTH(MAX_ATMOSWFS), L_DELS
      REAL(FPK) :: L_S1 (MAX_ATMOSWFS), L_S0 (MAX_ATMOSWFS)

!  Initialize

      DO N = 1, NLAYERS
        DO M = 1, N
          CHAPMAN_FACTORS(N,M) = zero
          DO Q = 1, NPARS
            L_CHAPMAN_FACTORS(Q,N,M) = zero
          ENDDO
        ENDDO
      ENDDO

!  Prepare spherical attenuation (shell geometry)

      IF ( .NOT.DO_PLANE_PARALLEL ) THEN

        GM_TOA = SQRT ( one - COS_SZA * COS_SZA )
        DO N = 0, NLAYERS
          H(N) = HEIGHT_GRID(N) + EARTH_RADIUS
          DO Q = 1, NPARS
            L_H(Q,N) = L_HEIGHT_GRID(Q,N)
          ENDDO
        ENDDO
        DO N = 1, NLAYERS
          DELZ(N) = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
          DO Q = 1, NPARS
            L_DELZ(Q,N) = L_HEIGHT_GRID(Q,N-1) - L_HEIGHT_GRID(Q,N)
          ENDDO
        ENDDO

        DO N = 1, NLAYERS
          STH = GM_TOA * H(N)/H(0)
          CTH = SQRT ( one - STH * STH )
          DO Q = 1, NPARS
             L_STH(Q) = STH * ( (L_H(Q,N)/H(N)) - (L_H(Q,0)/H(0)) )
             L_CTH(Q) = - STH * L_STH(Q) / CTH
          ENDDO
          S0 = zero
          HELP1 = H(0)*CTH
          HELP2 = -H(0)*H(0)*STH*STH
          DO Q = 1, NPARS
            L_S0(Q) = zero
            L_HELP1(Q) = L_H(Q,0)*CTH + H(0) * L_CTH(Q)
            L_HELP2(Q) = two*HELP2*((L_H(Q,0)/H(0))+(L_STH(Q)/STH))
          ENDDO
          DO M = 1, N
            HELP3 = SQRT(HELP2 + H(M)*H(M))
            S1 = HELP1 - HELP3
            DELS = S1 - S0
            CHAPMAN_FACTORS(N,M) = DELS / DELZ(M)
            DO Q = 1, NPARS
              L_HELP3 = half*(L_HELP2(Q)+two*H(M)*L_H(Q,M))/ HELP3
              L_S1(Q) = L_HELP1(Q) - L_HELP3
              L_DELS  = L_S1(Q) - L_S0(Q)
              L_CHAPMAN_FACTORS(Q,N,M) = CHAPMAN_FACTORS(N,M) * &
                           ( (L_DELS/DELS) - (L_DELZ(Q,M)/DELZ(M)) )
            ENDDO
            S0 = S1
            DO Q = 1, NPARS
              L_S0(Q) = L_S1(Q)
            ENDDO
          ENDDO
        ENDDO

!  Plane parallel

      ELSE

        DO N = 1, NLAYERS
          DO M = 1, N
            CHAPMAN_FACTORS(N,M) = one / COS_SZA
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE CHAPMAN_FUNCTION_PLUS

!

      SUBROUTINE LCH_ATTENUATION_SETUP_1 &
          ( DO_PLANE_PARALLEL, NPOINTS, NLAYERS, NPARS,           & ! Inputs
            DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS,        & ! Inputs
            DELTAU_VERT_INPUT, CHAPMAN, L_CHAPMAN, BEAM_PICUTOFF, & ! Inputs
            BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS,              & ! Inputs
            SAVE_TRANS_USERM, L_DELTAU_VERT_INPUT,                & ! Inputs
            L_BEAM_ITRANS, L_BEAM_AVSECANT,                       & ! Outputs
            L_BEAM_ETRANS, L_BEAM_DTRANS, L_SAVE_TRANS_USERM )      ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m   , Only : fpk, ZERO, ONE, TWO, HALF, MAX_TAU_UPATH, MAX_USER_STREAMS, &
                                 MAX_LAYERS, MAX_LAYERS_NK, MAX_ATMOSWFS, MAX_POINTS


      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Control

      LOGICAL  , INTENT(IN) :: DO_PLANE_PARALLEL

!  Control integers

      INTEGER  , INTENT(IN) :: NPOINTS, NLAYERS

!  Linearization control

      INTEGER  , INTENT(IN) :: NPARS

!  Control flag and inputs for user streams

      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS
      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  Chapman factors

      REAL(FPK), INTENT(IN) :: CHAPMAN ( MAX_LAYERS, MAX_LAYERS )

!  Linearized Input Chapman Factors

      REAL(FPK), INTENT(IN) :: L_CHAPMAN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LAYERS )

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT   ( MAX_LAYERS, MAX_POINTS )

!  Solar beam transmittances, average secant factors

      INTEGER  , INTENT(IN) :: BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(IN) :: SAVE_TRANS_USERM &
            ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT_INPUT &
                        ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

!  Output arguments
!  ----------------

!  Linearized Solar beam transmittances, average secant factors

      REAL(FPK), INTENT(OUT) :: L_BEAM_ITRANS &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_BEAM_AVSECANT &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_BEAM_ETRANS &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_BEAM_DTRANS &
                    ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )

!  Linearized user stream transmittances, whole layers

      REAL(FPK), INTENT(OUT) :: L_SAVE_TRANS_USERM &
            ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Local variables
!  ---------------

      INTEGER   :: N, K, S, UM, Q
      REAL(FPK) :: SPHER, SUMR, SM, L_SPHER, T1, T2
      REAL(FPK) :: TAUSLANT(MAX_LAYERS)
      REAL(FPK) :: L_TAUSLANT(MAX_ATMOSWFS,MAX_LAYERS)

!  Beam attenuation
!  ----------------

!  Start the points loop

      DO S = 1, NPOINTS

!  Tau slant = Beam slant optical thickness

        DO N = 1, NLAYERS
          SUMR = ZERO
          DO K = 1, N
            SUMR = SUMR + DELTAU_VERT_INPUT(K,S) * CHAPMAN(N,K)
          ENDDO
          TAUSLANT(N) = SUMR
          DO Q = 1, NPARS
            SUMR = ZERO
            DO K = 1, N
              SUMR = SUMR + L_DELTAU_VERT_INPUT(Q,K,S) * CHAPMAN(N,K) + &
                            DELTAU_VERT_INPUT(K,S)   * L_CHAPMAN(Q,N,K)
            ENDDO
            L_TAUSLANT(Q,N) = SUMR
          ENDDO
        ENDDO

!  First, linearize the initial transmittances

        L_BEAM_ITRANS(:,:,S) = ZERO
        DO N = 2, NLAYERS
          IF ( N.LE.BEAM_PICUTOFF(S) ) THEN
            DO Q = 1, NPARS
              L_BEAM_ITRANS(Q,N,S) = &
                    - L_TAUSLANT(Q,N-1) * BEAM_ITRANS(N,S)
            ENDDO
          ENDIF
        ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

        L_BEAM_AVSECANT(:,:,S) = ZERO
        IF ( .NOT. DO_PLANE_PARALLEL ) THEN
          N = 1
          DO Q = 1, NPARS
            L_BEAM_AVSECANT(Q,N,S) = L_CHAPMAN(Q,N,N)
          ENDDO
          DO N = 2, NLAYERS
            IF  (N.LE.BEAM_PICUTOFF(S) ) THEN
              DO Q = 1, NPARS
                T1 = L_TAUSLANT(Q,N) - L_TAUSLANT(Q,N-1)
                T2 = BEAM_AVSECANT(N,S) * L_DELTAU_VERT_INPUT(Q,N,S)
                L_BEAM_AVSECANT(Q,N,S) = &
                      (T1 - T2 ) / DELTAU_VERT_INPUT(N,S)
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Other variables

        DO N = 1, NLAYERS
          DO Q = 1, NPARS
            L_BEAM_DTRANS(Q,N,S) = - BEAM_DTRANS(N,S) * &
               ( BEAM_AVSECANT(N,S)   * L_DELTAU_VERT_INPUT(Q,N,S) + &
               L_BEAM_AVSECANT(Q,N,S) *   DELTAU_VERT_INPUT(N,S) )
            L_BEAM_ETRANS(Q,N,S) = &
                L_BEAM_DTRANS(Q,N,S) *   BEAM_ITRANS(N,S) + &
                  BEAM_DTRANS(N,S)   * L_BEAM_ITRANS(Q,N,S)
          ENDDO
        ENDDO

!  End points loop

      ENDDO

!  If no user streams, then return

      IF ( .NOT. DO_USER_STREAMS  ) RETURN

!  whole layer transmittances along user streams
!  =============================================

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          SM = ONE / USER_STREAMS(UM)
          DO S = 1, NPOINTS
            DO N = 1, NLAYERS
              SPHER = SM * DELTAU_VERT_INPUT(N,S)
              IF ( SPHER.GT.MAX_TAU_UPATH ) THEN
                DO Q = 1, NPARS
                  L_SAVE_TRANS_USERM(Q,UM,N,S) = ZERO
                ENDDO
              ELSE
                DO Q = 1, NPARS
                  L_SPHER = SM * L_DELTAU_VERT_INPUT(Q,N,S)
                  L_SAVE_TRANS_USERM(Q,UM,N,S) = - L_SPHER * &
                                       SAVE_TRANS_USERM(UM,N,S)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LCH_ATTENUATION_SETUP_1

!  End Module

      END MODULE lrrs_L_miscsetups_1_m

