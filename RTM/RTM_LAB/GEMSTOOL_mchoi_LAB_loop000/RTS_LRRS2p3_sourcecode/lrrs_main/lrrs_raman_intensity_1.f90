! ###############################################################
! #                                                             #
! #                    THE LIDORT_RRS MODEL                     #
! #                                                             #
! #      (LInearized Discrete Ordinate Radiative Transfer)      #
! #       --         -        -        -         -              #
! #                 (Rotational Raman Scatter)                  #
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
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Version      :  2.3                                        #
! #  Release Date :  March 2011                                 #
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
! #             RAMAN_INTENSITY_UP_1                            #
! #             RAMAN_INTENSITY_DN_1                            #
! #             MIFLUX_INTENSITY_1                              #
! #                                                             #
! #   Qudarature field  SUBROUTINES :                           #
! #                                                             #
! #             QUAD_INTENSITY_UP_1                             #
! #             QUAD_INTENSITY_DN_1                             #
! #                                                             #
! ###############################################################


      MODULE lrrs_raman_intensity_1

      PRIVATE
      PUBLIC :: RAMAN_INTENSITY_UP_1,&
                RAMAN_INTENSITY_DN_1,&
                MIFLUX_INTENSITY_1,&
                QUAD_INTENSITY_UP_1,&
                QUAD_INTENSITY_DN_1

      CONTAINS

      SUBROUTINE RAMAN_INTENSITY_UP_1 &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER, &
            FOURIER, NLAYERS, NSTREAMS, N_LOUTPUT, N_USER_STREAMS, &
            LEVELMASK_UP, WUPPER, WLOWER, KTRANS, &
            LCON_XVEC, MCON_XVEC, TRANS_USERM, &
            BOA_SOURCE, DIRECT_BOA_SOURCE, &
            WLAYER_HOST_UP, WLAYER_PIST_UP, &
            CUMSOURCE_UP, INTENSITY_F, QUADINTENS_F )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL, INTENT(IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::          DO_INCLUDE_MVOUT

!  FLux multiplier and Fourier component (counter)

      REAL(FPK), INTENT(IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT(IN) ::          FOURIER

!  Integer Control variables

      INTEGER, INTENT(IN) ::          NLAYERS
      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          N_LOUTPUT
      INTEGER, INTENT(IN) ::          N_USER_STREAMS

!  User defined stream angle input (whole layer/ Partial layers)

      REAL(FPK), INTENT(IN) :: &
          TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS )

!  optical depth indices/flags.

      INTEGER, INTENT(IN) ::          LEVELMASK_UP ( MAX_LOUTPUT )

!  Input source terms

      REAL(FPK), INTENT(IN) :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: WLAYER_HOST_UP ( MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: WLAYER_PIST_UP ( MAX_LAYERS, MAX_USER_STREAMS )

!  Variables for the particular integral (quadrature)

      REAL(FPK), INTENT(IN) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: WUPPER ( MAX_2_STREAMS, MAX_LAYERS )

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  homogeneous vectors with integration constants

      REAL(FPK), INTENT(IN) :: &
          LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  output
!  ------

!  local cumulative source terms

      REAL(FPK), INTENT(OUT) :: CUMSOURCE_UP ( MAX_USER_STREAMS, 0:MAX_LAYERS )

!  post processed field

      REAL(FPK), INTENT(OUT) :: INTENSITY_F ( MAX_LOUTPUT, MAX_USER_STREAMS )

!  Quadrature field

      REAL(FPK), INTENT(OUT) :: QUADINTENS_F ( MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

!  local help

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER ::          UTA, UM, M
      REAL(FPK) :: LST

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          CUMSOURCE_UP(UM,0) = BOA_SOURCE(UM) + DIRECT_BOA_SOURCE(UM)
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
        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO UM = 1, N_USER_STREAMS
              LST = WLAYER_HOST_UP(N,UM) + WLAYER_PIST_UP(N,UM)
              CUMSOURCE_UP(UM,NC) = LST + &
                         TRANS_USERM(UM,N)*CUMSOURCE_UP(UM,NC-1) 
            ENDDO
          ENDDO
        ENDIF

!  Ongrid output only
!  ------------------

!  Quadrature output at layer boundaries, for mean-value calculations

        IF ( DO_INCLUDE_MVOUT ) THEN
          CALL QUAD_INTENSITY_UP_1 &
                 ( NSTREAMS, NLAYERS, &
                   LEVELMASK_UP(UTA), UTA, FLUX_MULTIPLIER, &
                   WLOWER, WUPPER, KTRANS, LCON_XVEC, MCON_XVEC, &
                   QUADINTENS_F )
        ENDIF

!  User-defined stream output, just set to the cumulative source term

        IF ( DO_USER_STREAMS ) THEN
          DO UM = 1, N_USER_STREAMS
            INTENSITY_F(UTA,UM) = &
                      FLUX_MULTIPLIER * CUMSOURCE_UP(UM,NC)
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
      END SUBROUTINE RAMAN_INTENSITY_UP_1

!

      SUBROUTINE RAMAN_INTENSITY_DN_1 &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER, &
            FOURIER, NSTREAMS, N_LOUTPUT, N_USER_STREAMS, &
            LEVELMASK_DN, WLOWER, KTRANS, &
            LCON_XVEC, MCON_XVEC, TRANS_USERM, &
            TOA_SOURCE, WLAYER_HOST_DN, WLAYER_PIST_DN, &
            CUMSOURCE_DN, INTENSITY_F, QUADINTENS_F )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL, INTENT(IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::          DO_INCLUDE_MVOUT

!  FLux multiplier and Fourier component (counter)

      REAL(FPK), INTENT(IN) :: FLUX_MULTIPLIER
      INTEGER, INTENT(IN) ::          FOURIER

!  Integer Control variables

      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          N_LOUTPUT
      INTEGER, INTENT(IN) ::          N_USER_STREAMS

!  User defined stream angle input (whole layer/ Partial layers)

      REAL(FPK), INTENT(IN) :: &
          TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS )

!  optical depth indices

      INTEGER, INTENT(IN) ::          LEVELMASK_DN ( MAX_LOUTPUT )

!  Input source terms

      REAL(FPK), INTENT(IN) :: TOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: WLAYER_HOST_DN ( MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: WLAYER_PIST_DN ( MAX_LAYERS, MAX_USER_STREAMS )

!  Variables for the particular integral (quadrature)

      REAL(FPK), INTENT(IN) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  homogeneous vectors with integration constants

      REAL(FPK), INTENT(IN) :: &
          LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  output
!  ------

!  local cumulative source terms

      REAL(FPK), INTENT(OUT) :: CUMSOURCE_DN ( MAX_USER_STREAMS, 0:MAX_LAYERS )

!  post processed field

      REAL(FPK), INTENT(OUT) :: INTENSITY_F  ( MAX_LOUTPUT, MAX_USER_STREAMS )

!  Quadrature field

      REAL(FPK), INTENT(OUT) :: QUADINTENS_F ( MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NC, UTA, UM, M, NLEVEL
      REAL(FPK) :: LST

!  Initialize recursion
!  --------------------

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          CUMSOURCE_DN(UM,0) = TOA_SOURCE(UM)
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
            DO UM = 1, N_USER_STREAMS
              LST = WLAYER_HOST_DN(N,UM) + WLAYER_PIST_DN(N,UM)
              CUMSOURCE_DN(UM,NC) = LST + &
                        TRANS_USERM(UM,N)*CUMSOURCE_DN(UM,NC-1)
            ENDDO
          ENDDO
        ENDIF

!  Ongrid output only
!  ------------------

!  Quadrature results at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

        IF ( DO_INCLUDE_MVOUT ) THEN
          CALL QUAD_INTENSITY_DN_1 &
                 ( NSTREAMS, LEVELMASK_DN(UTA), UTA, FLUX_MULTIPLIER, &
                   WLOWER, KTRANS, LCON_XVEC, MCON_XVEC, &
                   QUADINTENS_F )
        ENDIF

!  User-defined stream output, just set to the cumulative source term

        IF ( DO_USER_STREAMS ) THEN
          DO UM = 1, N_USER_STREAMS
            INTENSITY_F(UTA,UM) = &
                        FLUX_MULTIPLIER * CUMSOURCE_DN(UM,NC)
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
      END SUBROUTINE RAMAN_INTENSITY_DN_1

!

      SUBROUTINE MIFLUX_INTENSITY_1 &
           ( DO_UPWELLING, DO_DNWELLING, &
             COS_SZA,  FLUX_FACTOR, &
             NSTREAMS, N_LOUTPUT, &
             LEVELMASK_DN, PICUTOFF, BEAM_ETRANS, &
             QUAD_WEIGHTS, QUAD_STRMWGHT, &
             QUADINTENS_UP, QUADINTENS_DN, &
             MEAN_INTENSITY_UP, MEAN_FLUX_UP, &
             MEAN_INTENSITY_DN, MEAN_FLUX_DN )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL, INTENT(IN) ::   DO_UPWELLING, DO_DNWELLING

!  nubmer of streams and user-defined output levels

      INTEGER, INTENT(IN) ::   NSTREAMS, N_LOUTPUT

!  Cosine SZA and flux factor

      REAL(FPK), INTENT(IN) :: COS_SZA, FLUX_FACTOR

!  Control for the user level output

      INTEGER, INTENT(IN) ::   LEVELMASK_DN ( MAX_LOUTPUT )

!  Direct beam transmittances

      REAL(FPK), INTENT(IN) :: BEAM_ETRANS ( MAX_LAYERS )

!  Cutoff for the direct beam

!      REAL(FPK), INTENT(IN) :: PICUTOFF  !  Bug found by JvG, 4/28/1
      INTEGER, INTENT(IN) ::   PICUTOFF

!  Quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS ( MAX_STREAMS )

!  Discrete ordinate radiances

      REAL(FPK), INTENT(IN) :: QUADINTENS_UP( MAX_LOUTPUT, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUADINTENS_DN( MAX_LOUTPUT, MAX_STREAMS )

!  output
!  ------

!  Upwelling mean intensity and mean flux

      REAL(FPK), INTENT(INOUT) :: MEAN_INTENSITY_UP ( MAX_LOUTPUT )
      REAL(FPK), INTENT(INOUT) :: MEAN_FLUX_UP      ( MAX_LOUTPUT )

!  Downwelling mean intensity and mean flux

      REAL(FPK), INTENT(INOUT) :: MEAN_INTENSITY_DN ( MAX_LOUTPUT )
      REAL(FPK), INTENT(INOUT) :: MEAN_FLUX_DN      ( MAX_LOUTPUT )

!  local variables
!  ---------------

      INTEGER ::   I, UTA, N
      REAL(FPK) :: SUM_MI, SUM_FX, FMU0, A, AX
      REAL(FPK) :: DIRECT_TRANS, DIRECT_FLUX, DIRECT_MEANI

!  mean intensity and flux
!  -----------------------

!  Upwelling

!   - integrate the radiance over the half space, using quadrature

      IF ( DO_UPWELLING ) THEN
        DO UTA = 1, N_LOUTPUT
          SUM_MI = ZERO
          SUM_FX = ZERO
          DO I = 1, NSTREAMS
            A  = QUAD_WEIGHTS (I)
            AX = QUAD_STRMWGHT(I)
            SUM_MI = SUM_MI + A  * QUADINTENS_UP(UTA,I)
            SUM_FX = SUM_FX + AX * QUADINTENS_UP(UTA,I)
          ENDDO
          MEAN_INTENSITY_UP(UTA) = SUM_MI * HALF
          MEAN_FLUX_UP(UTA)      = SUM_FX * PI2
        ENDDO
      ENDIF

!  Downwelling

      IF ( DO_DNWELLING ) THEN

        DO UTA = 1, N_LOUTPUT

!   - integrate the radiance over the half space, using quadrature

          SUM_MI = ZERO
          SUM_FX = ZERO
          DO I = 1, NSTREAMS
            A  = QUAD_WEIGHTS(I)
            AX = QUAD_STRMWGHT(I)
            SUM_MI = SUM_MI + A  * QUADINTENS_DN(UTA,I)
            SUM_FX = SUM_FX + AX * QUADINTENS_DN(UTA,I)
          ENDDO
          MEAN_INTENSITY_DN(UTA) = SUM_MI * HALF
          MEAN_FLUX_DN(UTA)      = SUM_FX * PI2

!  For the downward direction, add the direct beam contributions
!  loop over all the output optical depths ( on-grid Values only )

          N = LEVELMASK_DN(UTA)
          IF ( N .LE. PICUTOFF ) THEN
            IF ( N .EQ. 0 ) THEN
              DIRECT_TRANS = ONE
              FMU0 = COS_SZA * FLUX_FACTOR
            ELSE
              DIRECT_TRANS = BEAM_ETRANS(N)
              FMU0 = COS_SZA * FLUX_FACTOR
            ENDIF
            DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
            MEAN_INTENSITY_DN(UTA) = &
                        MEAN_INTENSITY_DN(UTA) + DIRECT_MEANI
            DIRECT_FLUX  = FMU0 * DIRECT_TRANS
            MEAN_FLUX_DN(UTA) = MEAN_FLUX_DN(UTA) + DIRECT_FLUX
          ENDIF

        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE MIFLUX_INTENSITY_1

!

      SUBROUTINE QUAD_INTENSITY_UP_1 &
           ( NSTREAMS, NLAYERS, NLEVEL, UTA, FLUX_MULTIPLIER, &
             WLOWER, WUPPER, KTRANS, LCON_XVEC, MCON_XVEC, &
             QUADINTENS )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  control

      INTEGER, INTENT(IN) ::          NSTREAMS, NLAYERS
      INTEGER, INTENT(IN) ::          NLEVEL, UTA

!  Flux factor

      REAL(FPK), INTENT(IN) :: FLUX_MULTIPLIER

!  particular integrals at layer boundaries

      REAL(FPK), INTENT(IN) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: WUPPER ( MAX_2_STREAMS, MAX_LAYERS )

!  Eigentransmittances, whole layer

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS,   MAX_LAYERS )

!  integration constants multiplied by homogeneous solutions

      REAL(FPK), INTENT(IN) :: &
            LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  output = discrete ordinate radiance

      REAL(FPK), INTENT(INOUT) :: QUADINTENS ( MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, I, I1, AA
      REAL(FPK) :: SPAR, SHOM, HOM1, HOM2

!  For those optical depths at layer boundaries
!  --------------------------------------------

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the intensity at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling intensity
!  at the bottom of the atmosphere (treated separately).

      N = NLEVEL + 1

!  homogeneous and particular solution contributions SHOM and SPAR

!  For the lowest level

      IF ( NLEVEL .EQ. NLAYERS ) THEN

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SHOM = ZERO
          DO AA = 1, NSTREAMS
            HOM1 = LCON_XVEC(I1,AA,NLEVEL) * KTRANS(AA,NLEVEL)
            HOM2 = MCON_XVEC(I1,AA,NLEVEL)
            SHOM = SHOM + HOM1 + HOM2
          ENDDO
          SPAR = WLOWER(I1,NLEVEL)
          QUADINTENS(UTA,I) = FLUX_MULTIPLIER * ( SPAR + SHOM )
        ENDDO

!  For other levels in the atmosphere

      ELSE

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SHOM = ZERO
          DO AA = 1, NSTREAMS
            HOM1 = LCON_XVEC(I1,AA,N)
            HOM2 = MCON_XVEC(I1,AA,N) * KTRANS(AA,N)
            SHOM = SHOM + HOM1 + HOM2
          ENDDO
          SPAR = WUPPER(I1,N)
          QUADINTENS(UTA,I) = FLUX_MULTIPLIER * ( SPAR + SHOM )
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUAD_INTENSITY_UP_1

!

      SUBROUTINE QUAD_INTENSITY_DN_1 &
           ( NSTREAMS, NLEVEL, UTA, FLUX_MULTIPLIER, &
             WLOWER, KTRANS, LCON_XVEC, MCON_XVEC, &
             QUADINTENS )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          NLEVEL, UTA

!  Flux factor

      REAL(FPK), INTENT(IN) :: FLUX_MULTIPLIER

!  Particular integral at lower boundary

      REAL(FPK), INTENT(IN) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )

!  Eigentransmittances, whole layer

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS,   MAX_LAYERS )

!  Integration constant stuff.

      REAL(FPK), INTENT(IN) :: &
            LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  output = discrete ordinate radiances

      REAL(FPK), INTENT(INOUT) :: QUADINTENS ( MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, I, AA
      REAL(FPK) :: SPAR, SHOM, HOM1, HOM2

!  For those optical depths at layer boundaries
!  --------------------------------------------

      N = NLEVEL

!  Downwelling radiation at TOA ( or N = 0 ) is zero

      IF ( NLEVEL .EQ. 0 ) THEN

        DO I = 1, NSTREAMS
          QUADINTENS(UTA,I) = ZERO
        ENDDO

!  Other levels

      ELSE

        DO I = 1, NSTREAMS
          SPAR = WLOWER(I,N)
          SHOM = ZERO
          DO AA = 1, NSTREAMS
            HOM1 = LCON_XVEC(I,AA,N) * KTRANS(AA,N)
            HOM2 = MCON_XVEC(I,AA,N)
            SHOM = SHOM + HOM1 + HOM2
          ENDDO
          QUADINTENS(UTA,I) = FLUX_MULTIPLIER * ( SPAR + SHOM )
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE QUAD_INTENSITY_DN_1


      END MODULE lrrs_raman_intensity_1

