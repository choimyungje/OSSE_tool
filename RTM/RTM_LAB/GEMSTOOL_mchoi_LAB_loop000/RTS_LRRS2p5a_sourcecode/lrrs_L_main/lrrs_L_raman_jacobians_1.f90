
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
! #   TOP LEVEL SUBROUTINES :                                   #
! #     (1 = Whole-layer, 2 = part-layer )                      #
! #                                                             #
! #             1. RAMAN_JACOBIAN_UP_1                          #
! #             2. RAMAN_JACOBIAN_DN_1                          #
! #             3. MIFLUX_JACOBIAN_1                            #
! #             4. RAMAN_SURFJAC_UP_1                           #
! #             5. RAMAN_SURFJAC_DN_1                           #
! #             6. MIFLUX_SURFJAC_1                             #
! #                                                             #
! #   Atmospheric Jacobians, quadrature output (called by 1/2)  #
! #                                                             #
! #             QUAD_JACOBIAN_UP_1                              #
! #             QUAD_JACOBIAN_DN_1                              #
! #                                                             #
! ###############################################################

      MODULE lrrs_L_raman_jacobians_1_m

!      USE LRRS_PARS_m, Only : LDU

!  Everything public

      PRIVATE
      PUBLIC :: RAMAN_JACOBIAN_UP_1, &
                RAMAN_JACOBIAN_DN_1, &
                MIFLUX_JACOBIAN_1,   &
                RAMAN_SURFJAC_UP_1,  &
                RAMAN_SURFJAC_DN_1,  &
                MIFLUX_SURFJAC_1,    &
                QUAD_JACOBIAN_UP_1,  &
                QUAD_JACOBIAN_DN_1

      CONTAINS

      SUBROUTINE RAMAN_JACOBIAN_UP_1 &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER,          & ! Inputs
            FOURIER, NLAYERS, NSTREAMS, N_LOUTPUT, N_USER_STREAMS,       & ! Inputs
            LEVELMASK_UP, KVARY, K, KPARS, NKSTORAGE2, TRANS_USERM,      & ! Inputs
            KTRANS, LCON, MCON, XPOS, XNEG, CUMSOURCE_UP, L_TRANS_USERM, & ! Inputs
            L_KTRANS, L_XPOS, L_XNEG, L_WUPPER, L_WLOWER,                & ! Inputs
            NCON, PCON, L_BOA_SOURCE, L_DIRECT_BOA_SOURCE,               & ! Inputs
            L_WLAYER_HOST_UP, L_WLAYER_PIST_UP,                          & ! Inputs
            JACOBIAN_F, QUADJACOBIAN_F )                                   ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_LAYERS_SQ, MAX_STREAMS, &
                              MAX_2_STREAMS, MAX_USER_STREAMS, MAX_LOUTPUT, MAX_ATMOSWFS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS
      LOGICAL  , INTENT(IN) :: DO_INCLUDE_MVOUT

!  FLux multiplier and Fourier component (counter)

      REAL(FPK), INTENT(IN) :: FLUX_MULTIPLIER
      INTEGER  , INTENT(IN) :: FOURIER

!  Integer Control variables

      INTEGER  , INTENT(IN) :: NLAYERS
      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: N_LOUTPUT
      INTEGER  , INTENT(IN) :: N_USER_STREAMS

!  Linearization control

      LOGICAL  , INTENT(IN) :: KVARY
      INTEGER  , INTENT(IN) :: K, KPARS
      INTEGER  , INTENT(IN) :: NKSTORAGE2(MAX_LAYERS,0:MAX_LAYERS)

!  User defined stream angle input (whole layer)

      REAL(FPK), INTENT(IN) :: TRANS_USERM (MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_TRANS_USERM &
                  ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS )

!  layer levels

      INTEGER  , INTENT(IN) :: LEVELMASK_UP ( MAX_LOUTPUT )

!  Linearized Input source terms

      REAL(FPK), INTENT(IN) :: L_BOA_SOURCE &
              ( MAX_ATMOSWFS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: L_DIRECT_BOA_SOURCE &
              ( MAX_ATMOSWFS, MAX_USER_STREAMS )

      REAL(FPK), INTENT(IN) :: L_WLAYER_HOST_UP &
              ( MAX_ATMOSWFS, MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: L_WLAYER_PIST_UP &
              ( MAX_ATMOSWFS, MAX_LAYERS_SQ, MAX_USER_STREAMS )

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  Cumulative source function

      REAL(FPK), INTENT(IN) :: CUMSOURCE_UP ( MAX_USER_STREAMS, 0:MAX_LAYERS)

!  homogeneous vectors and integration constants

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: LCON ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: MCON ( MAX_STREAMS, MAX_LAYERS )

!  Linearized solution terms

      REAL(FPK), INTENT(IN) :: L_WLOWER ( MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_WUPPER ( MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS,MAX_STREAMS,  MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: NCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: PCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_XPOS ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XNEG ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  output
!  ------
!mick fix 10/19/2015 - changed intent to inout

!  post processed field

      REAL(FPK), INTENT(INOUT) :: JACOBIAN_F &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS )

!  Quadrature field

      REAL(FPK), INTENT(INOUT) :: QUADJACOBIAN_F &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

!  local help

      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER   :: UTA, UM, Q, M, NK
      REAL(FPK) :: LST, FMULT

!  local cumulative source terms

      REAL(FPK) :: L_CUMSOURCE_UP ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  If no variation, exit with zeroed results

      IF ( .NOT. KVARY ) THEN
        DO UTA = 1,N_LOUTPUT
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, KPARS
              JACOBIAN_F(Q,K,UTA,UM) = zero
            ENDDO
          ENDDO
        ENDDO
        IF ( DO_INCLUDE_MVOUT ) THEN
          DO UTA = 1, N_LOUTPUT
            DO UM = 1, NSTREAMS
              DO Q = 1, KPARS
                QUADJACOBIAN_F(Q,K,UTA,UM) = zero
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            L_CUMSOURCE_UP(Q,UM) = &
              L_BOA_SOURCE(Q,UM) + L_DIRECT_BOA_SOURCE(Q,UM)
          ENDDO
        ENDDO
      ENDIF

!  Fourier (for debug only)

      M = FOURIER
      FMULT = FLUX_MULTIPLIER

!  initialise cumulative source term loop

      NC = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_LOUTPUT, 1, -1

!  Layer index for given optical depth

        NLEVEL = LEVELMASK_UP(UTA)

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            NK = NKSTORAGE2(N,K)
            DO UM = 1, N_USER_STREAMS
              IF ( N.EQ.K .OR. K.EQ.0 ) THEN
                DO Q = 1, KPARS
                  LST = L_WLAYER_HOST_UP(Q,N,UM) + &
                        L_WLAYER_PIST_UP(Q,NK,UM)
                  L_CUMSOURCE_UP(Q,UM) = LST &
                   + L_TRANS_USERM(Q,UM,N) *   CUMSOURCE_UP(UM,NC-1) &
                   +   TRANS_USERM(UM,N)   * L_CUMSOURCE_UP(Q,UM)
                ENDDO
              ELSE IF ( N.NE.K .AND. K.GT.0 ) THEN
                DO Q = 1, KPARS
                  LST = L_WLAYER_HOST_UP(Q,N,UM) + &
                        L_WLAYER_PIST_UP(Q,NK,UM)
                  L_CUMSOURCE_UP(Q,UM) = LST &
                     + TRANS_USERM(UM,N) * L_CUMSOURCE_UP(Q,UM)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF

!  Ongrid output
!  -------------

!  Quadrature output at layer boundaries, for mean-value calculations

        IF ( DO_INCLUDE_MVOUT ) THEN
          CALL QUAD_JACOBIAN_UP_1 &
                 ( NSTREAMS, NLAYERS, K, KPARS, NLEVEL, UTA, & ! Inputs
                   FMULT, KTRANS, LCON, MCON, XPOS, XNEG,    & ! Inputs
                   L_WLOWER, L_WUPPER, L_KTRANS,             & ! Inputs
                   NCON, PCON, L_XPOS, L_XNEG,               & ! Inputs
                   QUADJACOBIAN_F )                            ! Outputs
        ENDIF

!  User-defined stream output, just set to the cumulative source term

        IF ( DO_USER_STREAMS ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, KPARS
              JACOBIAN_F(Q,K,UTA,UM) = FMULT * L_CUMSOURCE_UP(Q,UM)
            ENDDO
          ENDDO
        ENDIF

!  Check for updating the recursion
!  --------------------------------

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE RAMAN_JACOBIAN_UP_1

!

      SUBROUTINE RAMAN_JACOBIAN_DN_1 &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER,           & ! Inputs
            FOURIER, NSTREAMS, N_LOUTPUT, N_USER_STREAMS,                 & ! Inputs
            LEVELMASK_DN, KVARY, K, KPARS, NKSTORAGE2, TRANS_USERM,       & ! Inputs
            KTRANS, LCON, MCON, XPOS, XNEG, CUMSOURCE_DN, L_TRANS_USERM,  & ! Inputs
            L_KTRANS, L_XPOS, L_XNEG, L_WLOWER, NCON, PCON, L_TOA_SOURCE, & ! Inputs
            L_WLAYER_HOST_DN, L_WLAYER_PIST_DN,                           & ! Inputs
            JACOBIAN_F, QUADJACOBIAN_F )                                    ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_LAYERS_SQ, MAX_STREAMS,  &
                              MAX_2_STREAMS, MAX_USER_STREAMS, MAX_LOUTPUT, MAX_ATMOSWFS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS
      LOGICAL  , INTENT(IN) :: DO_INCLUDE_MVOUT

!  FLux multiplier and Fourier component (counter)

      REAL(FPK), INTENT(IN) :: FLUX_MULTIPLIER
      INTEGER  , INTENT(IN) :: FOURIER

!  Integer Control variables

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: N_LOUTPUT
      INTEGER  , INTENT(IN) :: N_USER_STREAMS

!  Linearization control

      LOGICAL  , INTENT(IN) :: KVARY
      INTEGER  , INTENT(IN) :: K, KPARS
      INTEGER  , INTENT(IN) :: NKSTORAGE2(MAX_LAYERS,0:MAX_LAYERS)

!  User defined stream angle input (whole layer)

      REAL(FPK), INTENT(IN) :: TRANS_USERM (MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_TRANS_USERM &
                  ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS )

!  Partial layer levels

      INTEGER  , INTENT(IN) :: LEVELMASK_DN ( MAX_LOUTPUT )

!  Linearized Input source terms

      REAL(FPK), INTENT(IN) :: L_TOA_SOURCE &
              ( MAX_ATMOSWFS, MAX_USER_STREAMS )

      REAL(FPK), INTENT(IN) :: L_WLAYER_HOST_DN &
              ( MAX_ATMOSWFS, MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: L_WLAYER_PIST_DN &
              ( MAX_ATMOSWFS, MAX_LAYERS_SQ, MAX_USER_STREAMS )

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  Cumulative source function

      REAL(FPK), INTENT(IN) :: CUMSOURCE_DN ( MAX_USER_STREAMS, 0:MAX_LAYERS)

!  homogeneous vectors and integration constants

      REAL(FPK), INTENT(IN) :: &
          XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: &
          LCON ( MAX_STREAMS, MAX_LAYERS ), &
          MCON ( MAX_STREAMS, MAX_LAYERS )

!  Linearized solution terms

      REAL(FPK), INTENT(IN) :: L_WLOWER ( MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS,MAX_STREAMS,  MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: &
            NCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS ), &
            PCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_XPOS &
        ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XNEG &
        ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  output
!  ------
!mick fix 10/19/2015 - changed intent to inout

!  post processed field

      REAL(FPK), INTENT(INOUT) :: JACOBIAN_F &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS )

!  Quadrature field

      REAL(FPK), INTENT(INOUT) :: QUADJACOBIAN_F &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

!  local help

      INTEGER :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER :: UTA, UM, Q, M, NK
      REAL(FPK) :: LST

!  local cumulative source terms

      REAL(FPK) :: L_CUMSOURCE_DN &
             ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  If no variation, exit with zeroed results

      IF ( .NOT. KVARY ) THEN
        DO UTA = 1,N_LOUTPUT
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, KPARS
              JACOBIAN_F(Q,K,UTA,UM) = zero
            ENDDO
          ENDDO
        ENDDO
        IF ( DO_INCLUDE_MVOUT ) THEN
          DO UTA = 1, N_LOUTPUT
            DO UM = 1, NSTREAMS
              DO Q = 1, KPARS
                QUADJACOBIAN_F(Q,K,UTA,UM) = zero
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            L_CUMSOURCE_DN(Q,UM) = L_TOA_SOURCE(Q,UM)
          ENDDO
        ENDDO
      ENDIF

!  Fourier (for debug only)

      M = FOURIER

!  initialise cumulative source term loop

      NC = 0
      NSTART = 1
      NUT_PREV = NSTART - 1
      NUT = 0

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_LOUTPUT

!  Layer index for given optical depth

        NLEVEL = LEVELMASK_DN(UTA)
        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            NC = N
            NK = NKSTORAGE2(N,K)
            DO UM = 1, N_USER_STREAMS
              IF ( N.EQ.K .OR. K.EQ.0 ) THEN
                DO Q = 1, KPARS
                  LST = L_WLAYER_HOST_DN(Q,N,UM) + &
                        L_WLAYER_PIST_DN(Q,NK,UM)
                  L_CUMSOURCE_DN(Q,UM) = LST &
                   + L_TRANS_USERM(Q,UM,N) *   CUMSOURCE_DN(UM,NC-1) &
                   +   TRANS_USERM(UM,N)   * L_CUMSOURCE_DN(Q,UM)
                ENDDO
              ELSE IF ( N.NE.K .AND. K.GT.0 ) THEN
                DO Q = 1, KPARS
                  LST = L_WLAYER_HOST_DN(Q,N,UM) + &
                        L_WLAYER_PIST_DN(Q,NK,UM)
                  L_CUMSOURCE_DN(Q,UM) = LST &
                     + TRANS_USERM(UM,N) * L_CUMSOURCE_DN(Q,UM)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF

!  Ongrid output
!  -------------

!  Quadrature results at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

!  Quadrature output at layer boundaries, for mean-value calculations

        IF ( DO_INCLUDE_MVOUT ) THEN
          CALL QUAD_JACOBIAN_DN_1 &
              ( NSTREAMS, K, KPARS, NLEVEL, UTA, FLUX_MULTIPLIER, & ! Inputs
                KTRANS, LCON, MCON, XPOS, XNEG, L_WLOWER,         & ! Inputs
                L_KTRANS, NCON, PCON, L_XPOS, L_XNEG,             & ! Inputs
                QUADJACOBIAN_F )                                    ! Outputs
        ENDIF

!  User-defined stream output, just set to the cumulative source term

        IF ( DO_USER_STREAMS ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, KPARS
              JACOBIAN_F(Q,K,UTA,UM) = &
                        FLUX_MULTIPLIER * L_CUMSOURCE_DN(Q,UM)
            ENDDO
          ENDDO
        ENDIF

!  Check for updating the recursion
!  --------------------------------

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE RAMAN_JACOBIAN_DN_1

!

      SUBROUTINE MIFLUX_JACOBIAN_1 &
           ( DO_UPWELLING, DO_DNWELLING, NSTREAMS, N_LOUTPUT,    & ! Inputs
             LEVELMASK_DN, KVARY, K, KPARS, NKSTORAGE, COS_SZA,  & ! Inputs
             FLUX_FACTOR, QUAD_WEIGHTS, QUAD_STRMWGHT, PICUTOFF, & ! Inputs
             L_BEAM_ETRANS, QUADJACOBIAN_UP, QUADJACOBIAN_DN,    & ! Inputs
             MEAN_JACOBIAN_UP, FLUX_JACOBIAN_UP,                 & ! Outputs
             MEAN_JACOBIAN_DN, FLUX_JACOBIAN_DN )                  ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, PI2, PI4, HALF, MAX_LAYERS, MAX_STREAMS,  &
                              MAX_LOUTPUT, MAX_LAYERS_NK, MAX_ATMOSWFS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , INTENT(IN) :: DO_UPWELLING, DO_DNWELLING

!  nubmer of streams and user-defined output levels

      INTEGER  , INTENT(IN) :: NSTREAMS, N_LOUTPUT

!  Linearization control

      LOGICAL  , INTENT(IN) :: KVARY
      INTEGER  , INTENT(IN) :: K, KPARS
      INTEGER  , INTENT(IN) :: NKSTORAGE(MAX_LAYERS,0:MAX_LAYERS)

!  Cosine SZA and flux factor

      REAL(FPK), INTENT(IN) :: COS_SZA,  FLUX_FACTOR

!  Control for the user level output

      INTEGER  , INTENT(IN) :: LEVELMASK_DN ( MAX_LOUTPUT )

!  Direct beam transmittances

      REAL(FPK), INTENT(IN) :: L_BEAM_ETRANS ( MAX_ATMOSWFS, MAX_LAYERS_NK )

!  Cutoff for the direct beam

!      REAL(FPK), INTENT(IN) :: PICUTOFF             !  Bug found by JvG, 4/28/1
      INTEGER  , INTENT(IN) :: PICUTOFF

!  Quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS ( MAX_STREAMS )

!  Discrete ordinate Jacobians

      REAL(FPK), INTENT(IN) :: QUADJACOBIAN_UP &
             ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUADJACOBIAN_DN &
             ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS )

!  output
!  ------
!mick fix 10/19/2015 - changed intent to inout

!  Upwelling mean and flux Jacobians

      REAL(FPK), INTENT(INOUT) :: MEAN_JACOBIAN_UP &
                ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT )
      REAL(FPK), INTENT(INOUT) :: FLUX_JACOBIAN_UP &
                ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT )

!  Downwelling mean and flux Jacobians

      REAL(FPK), INTENT(INOUT) :: MEAN_JACOBIAN_DN &
                ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT )
      REAL(FPK), INTENT(INOUT) :: FLUX_JACOBIAN_DN &
                ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT )

!  local variables
!  ---------------

      INTEGER   :: I, UTA, N, Q, NK
      REAL(FPK) :: SUM_MI, SUM_FX, FMU0, A, AX
      REAL(FPK) :: DIRECT_TRANS, DIRECT_FLUX, DIRECT_MEANI

!  mean intensity and flux
!  -----------------------

!  If no results, zero and exit

      IF ( .NOT. KVARY ) THEN
        IF ( DO_UPWELLING ) THEN
          DO UTA = 1, N_LOUTPUT
            DO Q = 1, KPARS
              MEAN_JACOBIAN_UP(Q,K,UTA) = ZERO
              FLUX_JACOBIAN_UP(Q,K,UTA) = ZERO
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_DNWELLING ) THEN
          DO UTA = 1, N_LOUTPUT
            DO Q = 1, KPARS
              MEAN_JACOBIAN_DN(Q,K,UTA) = ZERO
              FLUX_JACOBIAN_DN(Q,K,UTA) = ZERO
            ENDDO
          ENDDO
        ENDIF
        RETURN
      ENDIF

!  Upwelling

!   - integrate the radiance over the half space, using quadrature

      IF ( DO_UPWELLING ) THEN
        DO UTA = 1, N_LOUTPUT
         DO Q = 1, KPARS
          SUM_MI = ZERO
          SUM_FX = ZERO
          DO I = 1, NSTREAMS
            A  = QUAD_WEIGHTS (I)
            AX = QUAD_STRMWGHT(I)
            SUM_MI = SUM_MI + A  * QUADJACOBIAN_UP(Q,K,UTA,I)
            SUM_FX = SUM_FX + AX * QUADJACOBIAN_UP(Q,K,UTA,I)
          ENDDO
          MEAN_JACOBIAN_UP(Q,K,UTA) = SUM_MI * HALF
          FLUX_JACOBIAN_UP(Q,K,UTA) = SUM_FX * PI2
         ENDDO
        ENDDO
      ENDIF

!  Downwelling

      IF ( DO_DNWELLING ) THEN

        DO UTA = 1, N_LOUTPUT
         DO Q = 1, KPARS

!   - integrate the radiance over the half space, using quadrature

          SUM_MI = ZERO
          SUM_FX = ZERO
          DO I = 1, NSTREAMS
            A  = QUAD_WEIGHTS(I)
            AX = QUAD_STRMWGHT(I)
            SUM_MI = SUM_MI + A  * QUADJACOBIAN_DN(Q,K,UTA,I)
            SUM_FX = SUM_FX + AX * QUADJACOBIAN_DN(Q,K,UTA,I)
          ENDDO
          MEAN_JACOBIAN_DN(Q,K,UTA) = SUM_MI * HALF
          FLUX_JACOBIAN_DN(Q,K,UTA) = SUM_FX * PI2

!  For the downward direction, add the direct beam contributions
!  loop over all the output optical depths ( on-grid Values only )

          N  = LEVELMASK_DN(UTA)
          IF ( N .LE. PICUTOFF ) THEN
            FMU0 = COS_SZA * FLUX_FACTOR
            IF ( N .EQ. 0 ) THEN
              DIRECT_TRANS = ZERO
             ELSE
              NK = NKSTORAGE(N,K)
              DIRECT_TRANS = L_BEAM_ETRANS(Q,NK)
            ENDIF
            DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
            MEAN_JACOBIAN_DN(Q,K,UTA) = &
                   MEAN_JACOBIAN_DN(Q,K,UTA) + DIRECT_MEANI
            DIRECT_FLUX  = FMU0 * DIRECT_TRANS
            FLUX_JACOBIAN_DN(Q,K,UTA) = &
                   FLUX_JACOBIAN_DN(Q,K,UTA) + DIRECT_FLUX
          ENDIF

         ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE MIFLUX_JACOBIAN_1

!

      SUBROUTINE RAMAN_SURFJAC_UP_1 &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FMULT,   & ! Inputs
            FOURIER, NLAYERS, NSTREAMS, N_LOUTPUT,      & ! Inputs
            N_USER_STREAMS, N_SURFACEWFS, LEVELMASK_UP, & ! Inputs
            TRANS_USERM, XPOS, XNEG, KTRANS,            & ! Inputs
            LS_WUPPER, LS_WLOWER, NCON_SURF, PCON_SURF, & ! Inputs
            LS_BOA_SOURCE, LS_DIRECT_BOA_SOURCE,        & ! Inputs
            LS_WLAYER_HOST_UP, LS_WLAYER_PIST_UP,       & ! Inputs
            JACOBIAN_F, QUADJACOBIAN_F )                  ! Outputs

!  Version 2.5. Introduced variable number of Surface Jacobians and dimensioning to match.

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_2_STREAMS, MAX_STREAMS,  &
                              MAX_USER_STREAMS, MAX_LOUTPUT, MAX_SURFACEWFS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS
      LOGICAL  , INTENT(IN) :: DO_INCLUDE_MVOUT

!  FLux multiplier and Fourier component (counter)

      REAL(FPK), INTENT(IN) :: FMULT
      INTEGER  , INTENT(IN) :: FOURIER

!  Integer Control variables

      INTEGER  , INTENT(IN) :: NLAYERS
      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: N_LOUTPUT
      INTEGER  , INTENT(IN) :: N_USER_STREAMS

!  Number of surface weighting functions

      INTEGER  , INTENT(IN) :: N_SURFACEWFS

!  User defined stream angle input (whole layer)

      REAL(FPK), INTENT(IN) :: TRANS_USERM (MAX_USER_STREAMS, MAX_LAYERS )

!  Levels of output

      INTEGER  , INTENT(IN) :: LEVELMASK_UP ( MAX_LOUTPUT )

!  Linearized Input source terms

      REAL(FPK), INTENT(IN) :: LS_BOA_SOURCE        ( MAX_SURFACEWFS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: LS_DIRECT_BOA_SOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS )

      REAL(FPK), INTENT(IN) :: LS_WLAYER_HOST_UP &
              ( MAX_SURFACEWFS, MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: LS_WLAYER_PIST_UP &
              ( MAX_SURFACEWFS, MAX_LAYERS, MAX_USER_STREAMS )

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  homogeneous vectors and integration constants

      REAL(FPK), INTENT(IN) :: &
          XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Linearized solution terms

      REAL(FPK), INTENT(IN) :: LS_WLOWER ( MAX_SURFACEWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: LS_WUPPER ( MAX_SURFACEWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: NCON_SURF ( MAX_SURFACEWFS, MAX_STREAMS,   MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: PCON_SURF ( MAX_SURFACEWFS, MAX_STREAMS,   MAX_LAYERS )

!  output
!  ------

!  post processed field

      REAL(FPK), INTENT(OUT) :: JACOBIAN_F ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS )

!  Quadrature field

      REAL(FPK), INTENT(OUT) :: QUADJACOBIAN_F ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

!  local help

      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER   :: UTA, UM, AA, I, I1, NL, M, Q
      REAL(FPK) :: LST, H1, H2, SPAR, SHOM

!  local cumulative source terms

      REAL(FPK) :: LS_CUMSOURCE_UP ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SURFACEWFS
            LS_CUMSOURCE_UP(Q,UM) = LS_BOA_SOURCE(Q,UM) + LS_DIRECT_BOA_SOURCE(Q,UM)
          ENDDO
        ENDDO
      ENDIF

!  Fourier (for debug only)

      M = FOURIER

!  initialise cumulative source term loop

      NC = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_LOUTPUT, 1, -1

!  Layer index for given optical depth

        NLEVEL = LEVELMASK_UP(UTA)

!  Recursion for whole layers

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, N_SURFACEWFS
                LST = LS_WLAYER_HOST_UP(Q,N,UM) + LS_WLAYER_PIST_UP(Q,N,UM)
                LS_CUMSOURCE_UP(Q,UM) = LST + TRANS_USERM(UM,N) * LS_CUMSOURCE_UP(Q,UM)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Ongrid output
!  -------------

!  Quadrature output at layer boundaries, for mean-value calculations

        IF ( DO_INCLUDE_MVOUT ) THEN

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

          NL = NLEVEL
          N = NL + 1

!  For the lowest level

          IF ( NLEVEL .EQ. NLAYERS ) THEN

            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              DO Q = 1, N_SURFACEWFS
                SPAR = LS_WLOWER(Q,I1,NL)
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 = NCON_SURF(Q,AA,NL)*XPOS(I1,AA,NL)*KTRANS(AA,NL)
                  H2 = PCON_SURF(Q,AA,NL)*XNEG(I1,AA,NL)
                  SHOM = SHOM + H1 + H2
                ENDDO
                QUADJACOBIAN_F(Q,UTA,I) = FMULT * ( SHOM + SPAR )
              ENDDO
            ENDDO

!  For other levels in the atmosphere

          ELSE

            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              DO Q = 1, N_SURFACEWFS
                SPAR = LS_WUPPER(Q,I1,N)
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 = NCON_SURF(Q,AA,N)*XPOS(I1,AA,N)
                  H2 = PCON_SURF(Q,AA,N)*XNEG(I1,AA,N)*KTRANS(AA,N)
                  SHOM = SHOM + H1 + H2
                ENDDO
                QUADJACOBIAN_F(Q,UTA,I) = FMULT * ( SHOM + SPAR )
              ENDDO
            ENDDO

          ENDIF

        ENDIF

!  User-defined stream output, just set to the cumulative source term

        IF ( DO_USER_STREAMS ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SURFACEWFS
              JACOBIAN_F(Q,UTA,UM) = FMULT * LS_CUMSOURCE_UP(Q,UM)
            ENDDO
          ENDDO
        ENDIF

!  Check for updating the recursion
!  --------------------------------

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE RAMAN_SURFJAC_UP_1

!

      SUBROUTINE RAMAN_SURFJAC_DN_1 &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER,  & ! Inputs
            FOURIER, NSTREAMS, N_LOUTPUT, N_USER_STREAMS,        & ! Inputs
            N_SURFACEWFS, LEVELMASK_DN, TRANS_USERM,             & ! Inputs
            XPOS, XNEG, KTRANS, LS_WLOWER, NCON_SURF, PCON_SURF, & ! Inputs
            LS_TOA_SOURCE, LS_WLAYER_HOST_DN, LS_WLAYER_PIST_DN, & ! Inputs
            JACOBIAN_F, QUADJACOBIAN_F )                           ! Outputs

!  Version 2.5. Introduced variable number of Surface Jacobians and dimensioning to match.

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_2_STREAMS, MAX_STREAMS,  &
                              MAX_USER_STREAMS, MAX_LOUTPUT, MAX_SURFACEWFS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS
      LOGICAL  , INTENT(IN) :: DO_INCLUDE_MVOUT

!  FLux multiplier and Fourier component (counter)

      REAL(FPK), INTENT(IN) :: FLUX_MULTIPLIER
      INTEGER  , INTENT(IN) :: FOURIER

!  Integer Control variables

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: N_LOUTPUT
      INTEGER  , INTENT(IN) :: N_USER_STREAMS

!  Number of surface weighting functions

      INTEGER  , INTENT(IN) :: N_SURFACEWFS

!  User defined stream angle input (whole layer)

      REAL(FPK), INTENT(IN) :: TRANS_USERM (MAX_USER_STREAMS, MAX_LAYERS )

!  Level input

      INTEGER  , INTENT(IN) :: LEVELMASK_DN ( MAX_LOUTPUT )

!  Linearized Input source terms

      REAL(FPK), INTENT(IN) :: LS_TOA_SOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS )

      REAL(FPK), INTENT(IN) :: LS_WLAYER_HOST_DN ( MAX_SURFACEWFS, MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: LS_WLAYER_PIST_DN ( MAX_SURFACEWFS, MAX_LAYERS, MAX_USER_STREAMS )

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  homogeneous vectors and integration constants

      REAL(FPK), INTENT(IN) :: &
          XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Linearized solution terms

      REAL(FPK), INTENT(IN) :: LS_WLOWER  ( MAX_SURFACEWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: NCON_SURF  ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: PCON_SURF  ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS )

!  output
!  ------

!  post processed field

      REAL(FPK), INTENT(OUT) :: JACOBIAN_F ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS )

!  Quadrature field

      REAL(FPK), INTENT(OUT) :: QUADJACOBIAN_F ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

!  local help

      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER   :: UTA, UM, AA, I, NL, M, Q
      REAL(FPK) :: LST, SHOM, SPAR, H1, H2

!  local cumulative source terms

      REAL(FPK) :: LS_CUMSOURCE_DN ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SURFACEWFS
            LS_CUMSOURCE_DN(Q,UM) = LS_TOA_SOURCE(Q,UM)
          ENDDO
        ENDDO
      ENDIF

!  Fourier (for debug only)

      M = FOURIER

!  initialise cumulative source term loop

      NC = 0
      NSTART = 1
      NUT_PREV = NSTART - 1
      NUT = 0

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_LOUTPUT

!  Layer index

        NLEVEL = LEVELMASK_DN(UTA)

!  Do the recursion for one layer

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            NC = N
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, N_SURFACEWFS
                LST = LS_WLAYER_HOST_DN(Q,N,UM) + LS_WLAYER_PIST_DN(Q,N,UM)
                LS_CUMSOURCE_DN(Q,UM) = LST + TRANS_USERM(UM,N) * LS_CUMSOURCE_DN(Q,UM)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Ongrid output
!  -------------

!  Quadrature results at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

!  Quadrature output at layer boundaries, for mean-value calculations

        IF ( DO_INCLUDE_MVOUT ) THEN

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

          NL = NLEVEL
          N = NL
          IF ( NLEVEL .EQ. 0 ) THEN
            DO I = 1, NSTREAMS
              DO Q = 1, N_SURFACEWFS
                QUADJACOBIAN_F(Q,UTA,I) = ZERO
              ENDDO
            ENDDO
          ELSE
            DO I = 1, NSTREAMS
              DO Q = 1, N_SURFACEWFS
                SPAR = LS_WLOWER(Q,I,N)
                SHOM = ZERO
                  DO AA = 1, NSTREAMS
                  H1 = NCON_SURF(Q,AA,N)*XPOS(I,AA,N)*KTRANS(AA,N)
                  H2 = PCON_SURF(Q,AA,N)*XNEG(I,AA,N)
                  SHOM = SHOM + H1 + H2
                ENDDO
                QUADJACOBIAN_F(Q,UTA,I) = FLUX_MULTIPLIER * ( SHOM + SPAR )
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  User-defined stream output, just set to the cumulative source term

        IF ( DO_USER_STREAMS ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SURFACEWFS
              JACOBIAN_F(Q,UTA,UM) = FLUX_MULTIPLIER * LS_CUMSOURCE_DN(Q,UM)
            ENDDO
          ENDDO
        ENDIF

!  Check for updating the recursion
!  --------------------------------

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE RAMAN_SURFJAC_DN_1

!

      SUBROUTINE MIFLUX_SURFJAC_1 &
           ( DO_UPWELLING, DO_DNWELLING, NSTREAMS, N_LOUTPUT, & ! Inputs
             N_SURFACEWFS, QUAD_WEIGHTS, QUAD_STRMWGHT,       & ! Inputs
             QUADJACOBIAN_UP, QUADJACOBIAN_DN,                & ! Inputs
             MEAN_JACOBIAN_UP, FLUX_JACOBIAN_UP,              & ! Outputs
             MEAN_JACOBIAN_DN, FLUX_JACOBIAN_DN )               ! Outputs

!  Version 2.5. Introduced variable number of Surface Jacobians and dimensioning to match.

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, PI2, HALF, MAX_LAYERS, MAX_STREAMS,  &
                              MAX_LOUTPUT, MAX_SURFACEWFS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , INTENT(IN) :: DO_UPWELLING, DO_DNWELLING

!  nubmer of streams and user-defined output levels

      INTEGER  , INTENT(IN) :: NSTREAMS, N_LOUTPUT

!  Number of surface weighting functions

      INTEGER  , INTENT(IN) :: N_SURFACEWFS

!  Quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS ( MAX_STREAMS )

!  Discrete ordinate Jacobians

      REAL(FPK), INTENT(IN) :: QUADJACOBIAN_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUADJACOBIAN_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_STREAMS )

!  output
!  ------
!mick fix 10/19/2015 - changed intent to inout

!  Upwelling mean and flux Jacobians

      REAL(FPK), INTENT(INOUT) :: MEAN_JACOBIAN_UP ( MAX_SURFACEWFS, MAX_LOUTPUT )
      REAL(FPK), INTENT(INOUT) :: FLUX_JACOBIAN_UP ( MAX_SURFACEWFS, MAX_LOUTPUT )

!  Downwelling mean and flux Jacobians

      REAL(FPK), INTENT(INOUT) :: MEAN_JACOBIAN_DN ( MAX_SURFACEWFS, MAX_LOUTPUT )
      REAL(FPK), INTENT(INOUT) :: FLUX_JACOBIAN_DN ( MAX_SURFACEWFS, MAX_LOUTPUT )

!  local variables
!  ---------------

      INTEGER   :: UTA, Q
      REAL(FPK) :: SUM_MI, SUM_FX

!  mean intensity and flux
!  -----------------------

!  Upwelling

      IF ( DO_UPWELLING ) THEN
        DO UTA = 1, N_LOUTPUT
          DO Q = 1, N_SURFACEWFS
            SUM_MI = DOT_PRODUCT(QUAD_WEIGHTS (1:NSTREAMS),QUADJACOBIAN_UP(Q,UTA,1:NSTREAMS))
            SUM_FX = DOT_PRODUCT(QUAD_STRMWGHT(1:NSTREAMS),QUADJACOBIAN_UP(Q,UTA,1:NSTREAMS))
            MEAN_JACOBIAN_UP(Q,UTA) = SUM_MI * HALF
            FLUX_JACOBIAN_UP(Q,UTA) = SUM_FX * PI2
          ENDDO
        ENDDO
      ENDIF

!  Downwelling

      IF ( DO_DNWELLING ) THEN
        DO UTA = 1, N_LOUTPUT
          DO Q = 1, N_SURFACEWFS
            SUM_MI = DOT_PRODUCT(QUAD_WEIGHTS (1:NSTREAMS),QUADJACOBIAN_DN(Q,UTA,1:NSTREAMS))
            SUM_FX = DOT_PRODUCT(QUAD_STRMWGHT(1:NSTREAMS),QUADJACOBIAN_DN(Q,UTA,1:NSTREAMS))
            MEAN_JACOBIAN_DN(Q,UTA) = SUM_MI * HALF
            FLUX_JACOBIAN_DN(Q,UTA) = SUM_FX * PI2
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE MIFLUX_SURFJAC_1

!

      SUBROUTINE QUAD_JACOBIAN_UP_1 &
           ( NSTREAMS, NLAYERS, K, KPARS,    & ! Inputs
             NL, UTA, FLUX_MULTIPLIER,       & ! Inputs
             KTRANS, LCON, MCON, XPOS, XNEG, & ! Inputs
             L_WLOWER, L_WUPPER, L_KTRANS,   & ! Inputs
             NCON, PCON, L_XPOS, L_XNEG,     & ! Inputs
             QUADJACOBIAN )                    ! Output

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_STREAMS,  MAX_2_STREAMS,  &
                              MAX_LOUTPUT, MAX_ATMOSWFS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  control

      INTEGER  , INTENT(IN) :: NSTREAMS, NLAYERS
      INTEGER  , INTENT(IN) :: NL, UTA

!  Linearization control

      INTEGER  , INTENT(IN) :: K, KPARS

!  Flux factor

      REAL(FPK), INTENT(IN) :: FLUX_MULTIPLIER

!  Eigentransmittances, whole layer

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS,   MAX_LAYERS )

!  integration constants, homogeneous solutions

      REAL(FPK), INTENT(IN) :: LCON ( MAX_STREAMS,   MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: MCON ( MAX_STREAMS,   MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Linearizations

      REAL(FPK), INTENT(IN) :: L_WLOWER ( MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_WUPPER ( MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS,MAX_STREAMS,  MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: NCON     ( MAX_ATMOSWFS,MAX_STREAMS,  MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: PCON     ( MAX_ATMOSWFS,MAX_STREAMS,  MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XPOS &
          ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XNEG &
          ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  output = discrete ordinate Jacobian

!mick fix 10/19/2015 - changed intent to inout
      REAL(FPK), INTENT(INOUT) :: QUADJACOBIAN &
           ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

      INTEGER   :: N, I, I1, AA, Q
      REAL(FPK) :: SPAR, SHOM, H1, H2, H3, H4, H5
      REAL(FPK) :: LXP, L_LXP

!  For those optical depths at layer boundaries
!  --------------------------------------

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the intensity at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling intensity
!  at the bottom of the atmosphere (treated separately).

      N  = NL + 1

!  For the lowest level

      IF ( NL .EQ. NLAYERS ) THEN

!  If this is also the layer that is varying, extra contributions

        IF ( K .EQ. NL .OR. K .EQ. 0 ) THEN
          DO Q = 1, KPARS
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                LXP   = LCON(AA,NL) * XPOS(I1,AA,NL)
                L_LXP = NCON(Q,AA,NL) *   XPOS(I1,AA,NL) + &
                        LCON(AA,NL)   * L_XPOS(Q,I1,AA,NL)
                H1   = L_LXP *   KTRANS(AA,NL) + &
                         LXP * L_KTRANS(Q,AA,NL)
                H2   =  PCON(Q,AA,NL) *   XNEG(I1,AA,NL) + &
                        MCON(AA,NL)   * L_XNEG(Q,I1,AA,NL)
                SHOM = SHOM + H1 + H2
              ENDDO
              SPAR = L_WLOWER(Q,I1,NL)
              QUADJACOBIAN(Q,K,UTA,I) = FLUX_MULTIPLIER * (SPAR+SHOM)
            ENDDO
          ENDDO
        ENDIF

!  non-varying lowest layer

        IF ( K.NE.NL .AND. K.GT.0 ) THEN
          DO Q = 1, KPARS
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                L_LXP = NCON(Q,AA,NL) * XPOS(I1,AA,NL)
                H1 = L_LXP * KTRANS(AA,NL)
                H2 =  PCON(Q,AA,NL) *   XNEG(I1,AA,NL)
                SHOM = SHOM + H1 + H2
              ENDDO
              SPAR = L_WLOWER(Q,I1,NL)
              QUADJACOBIAN(Q,K,UTA,I) = FLUX_MULTIPLIER * (SPAR+SHOM)
            ENDDO
          ENDDO
        ENDIF

!  For other levels in the atmosphere
!  ----------------------------------

      ELSE

!  If this is also the layer that is varying, extra contributions

        IF ( K.EQ.N  .OR. K.EQ.0 ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, KPARS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1 = NCON(Q,AA,N) *   XPOS(I1,AA,N)
                H2 = LCON(AA,N)   * L_XPOS(Q,I1,AA,N)
                H3 = MCON(AA,N)   *   XNEG(I1,AA,N)   * L_KTRANS(Q,AA,N)
                H4 = MCON(AA,N)   * L_XNEG(Q,I1,AA,N) * KTRANS(AA,N)
                H5 = PCON(Q,AA,N) *   XNEG(I1,AA,N)   * KTRANS(AA,N)
                SHOM = SHOM + H1 + H2 + H3 + H4 + H5
              ENDDO
              SPAR = L_WUPPER(Q,I1,N)
              QUADJACOBIAN(Q,K,UTA,I) = FLUX_MULTIPLIER * (SPAR+SHOM)
            ENDDO
          ENDDO

!  non-varying layer lower than the varying layer

        ELSE IF ( K.LT.N .AND. K.NE.0 ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, KPARS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1 = NCON(Q,AA,N) *   XPOS(I1,AA,N)
                H2 = PCON(Q,AA,N) *   XNEG(I1,AA,N) * KTRANS(AA,N)
                SHOM = SHOM + H1 + H2
              ENDDO
              SPAR = L_WUPPER(Q,I1,N)
              QUADJACOBIAN(Q,K,UTA,I) = FLUX_MULTIPLIER * (SPAR+SHOM)
            ENDDO
          ENDDO

!  non-varying layer higher than the varying layer

        ELSE IF ( K.GT.N .AND. K.NE.0 ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, KPARS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1 = NCON(Q,AA,N) *   XPOS(I1,AA,N)
                H2 = PCON(Q,AA,N) *   XNEG(I1,AA,N) * KTRANS(AA,N)
                SHOM = SHOM + H1 + H2
              ENDDO
              QUADJACOBIAN(Q,K,UTA,I) = FLUX_MULTIPLIER * SHOM
            ENDDO
          ENDDO

        ENDIF

!  ENd clause for levels in atmosphere

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUAD_JACOBIAN_UP_1

!

      SUBROUTINE QUAD_JACOBIAN_DN_1 &
           ( NSTREAMS, K, KPARS, NL, UTA, FLUX_MULTIPLIER, & ! Inputs
             KTRANS, LCON, MCON, XPOS, XNEG, L_WLOWER,     & ! Inputs
             L_KTRANS, NCON, PCON, L_XPOS, L_XNEG,         & ! Inputs
             QUADJACOBIAN )                                  ! Output

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_STREAMS, MAX_2_STREAMS,  &
                              MAX_LOUTPUT, MAX_ATMOSWFS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: NL, UTA

!  Linearization control

      INTEGER  , INTENT(IN) :: K, KPARS

!  Flux factor

      REAL(FPK), INTENT(IN) :: FLUX_MULTIPLIER

!  Eigentransmittances, whole layer

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS,   MAX_LAYERS )

!  integration constants, homogeneous solutions

      REAL(FPK), INTENT(IN) :: LCON ( MAX_STREAMS,   MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: MCON ( MAX_STREAMS,   MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Linearizations

      REAL(FPK), INTENT(IN) :: L_WLOWER ( MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS,MAX_STREAMS,  MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: NCON     ( MAX_ATMOSWFS,MAX_STREAMS,  MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: PCON     ( MAX_ATMOSWFS,MAX_STREAMS,  MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XPOS &
          ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XNEG &
          ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  output = discrete ordinate Jacobian

!mick fix 10/19/2015 - changed intent to inout
      REAL(FPK), INTENT(INOUT) :: QUADJACOBIAN &
           ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

      INTEGER   :: I, AA, Q, N
      REAL(FPK) :: SPAR, SHOM, H1, H2, LXP, L_LXP

!  For those optical depths at layer boundaries
!  --------------------------------------

!  Downwelling radiation at TOA ( or N = 0 ) is zero

      IF ( NL .EQ. 0 ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, KPARS
            QUADJACOBIAN(Q,K,UTA,I) = ZERO
          ENDDO
        ENDDO

!  Other levels

      ELSE

        N = NL

!  If this is also the layer that is varying, extra contributions

        IF ( K .EQ. N .OR. K .EQ. 0 ) THEN

          DO Q = 1, KPARS
            DO I = 1, NSTREAMS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                LXP   = LCON(AA,N) * XPOS(I,AA,N)
                L_LXP = NCON(Q,AA,N) *   XPOS(I,AA,N) + &
                        LCON(AA,N)   * L_XPOS(Q,I,AA,N)
                H1   = L_LXP *   KTRANS(AA,N) + &
                         LXP * L_KTRANS(Q,AA,N)
                H2   =  PCON(Q,AA,N) *   XNEG(I,AA,N) + &
                        MCON(AA,N)   * L_XNEG(Q,I,AA,N)
                SHOM = SHOM + H1 + H2
              ENDDO
              SPAR = L_WLOWER(Q,I,N)
              QUADJACOBIAN(Q,K,UTA,I) = FLUX_MULTIPLIER * (SPAR+SHOM)
            ENDDO
          ENDDO

!  non-varying layer lower than the varying layer

        ELSE IF ( K.LT.N .AND. K.NE.0 ) THEN

          DO Q = 1, KPARS
            DO I = 1, NSTREAMS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1   =  NCON(Q,AA,N) *   XPOS(I,AA,N) * KTRANS(AA,N)
                H2   =  PCON(Q,AA,N) *   XNEG(I,AA,N)
                SHOM = SHOM + H1 + H2
              ENDDO
              SPAR = L_WLOWER(Q,I,N)
              QUADJACOBIAN(Q,K,UTA,I) = FLUX_MULTIPLIER * (SPAR+SHOM)
            ENDDO
          ENDDO

!  non-varying layer higher than the varying layer

        ELSE IF ( K.GT.N .AND. K.NE.0 ) THEN

          DO Q = 1, KPARS
            DO I = 1, NSTREAMS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1   =  NCON(Q,AA,N) *   XPOS(I,AA,N) * KTRANS(AA,N)
                H2   =  PCON(Q,AA,N) *   XNEG(I,AA,N)
                SHOM = SHOM + H1 + H2
              ENDDO
              QUADJACOBIAN(Q,K,UTA,I) = FLUX_MULTIPLIER * SHOM
            ENDDO
          ENDDO

        ENDIF

!  End choice of levels

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUAD_JACOBIAN_DN_1

!  End module

      END MODULE lrrs_L_raman_jacobians_1_m


