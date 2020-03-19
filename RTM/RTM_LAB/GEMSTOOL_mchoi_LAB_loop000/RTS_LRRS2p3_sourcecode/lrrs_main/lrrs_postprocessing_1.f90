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
! #   SUBROUTINES :                                             #
! #       ( 1 = Whole-layer, 2 = partial-layer )                #
! #                                                             #
! #             POSTPROCESSING_1                                #
! #                                                             #
! #             HOMOGMULT_1                                     #
! #             BEAMMULT_UP_1                                   #
! #             BEAMMULT_DN_1                                   #
! #                                                             #
! #             TOASOURCE                                       #
! #             BOASOURCE                                       #
! #                                                             #
! #             RAMAN_SOURCETERM_UP_1                           #
! #             RAMAN_SOURCETERM_DN_1                           #
! #                                                             #
! ###############################################################


      MODULE lrrs_postprocessing_1

      PRIVATE
      PUBLIC :: POSTPROCESSING_1,&
                HOMOGMULT_1,&
                BEAMMULT_UP_1,&
                BEAMMULT_DN_1,&
                TOASOURCE,&
                BOASOURCE

      CONTAINS

      SUBROUTINE POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING, &
            DO_MSMODE_LIDORT, DO_SSCORRECTION, &
            FOURIER, N, NSTREAMS, N_USER_STREAMS, &
            STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
            FGSMALL, TRANS_USERM, DELTAUS, USER_STREAMS, &
            FLIPPER, CHANGE_EXPONENT, SOLUTION, &
            PICUTOFF, ASOURCE, DELTRANS, INITRANS, &
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, &
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG, &
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV, &
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )

!  Source term post processing for one layer

!  include file of dimensions and numbers

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Up or down

      LOGICAL, INTENT(IN) ::          DO_UPWELLING
      LOGICAL, INTENT(IN) ::          DO_DNWELLING

!  Flags for multiple scattering and single scatter correction

      LOGICAL, INTENT(IN) ::          DO_MSMODE_LIDORT
      LOGICAL, INTENT(IN) ::          DO_SSCORRECTION

!  Fourier component

      INTEGER, INTENT(IN) ::          FOURIER

!  Number of layer

      INTEGER, INTENT(IN) ::          N

!  number of streams and solutions

      INTEGER, INTENT(IN) ::          NSTREAMS, N_USER_STREAMS

!  layer control

      LOGICAL, INTENT(IN) ::          STERM_LAYERMASK_UP
      LOGICAL, INTENT(IN) ::          STERM_LAYERMASK_DN

!  Criterion for small numbers analysis

      REAL(FPK), INTENT(IN) :: FGSMALL

!  Control for the particular solutions

      LOGICAL, INTENT(IN) ::          CHANGE_EXPONENT
      LOGICAL, INTENT(IN) ::          FLIPPER
      INTEGER, INTENT(IN) ::          SOLUTION

!  Optical thickness value

      REAL(FPK), INTENT(IN) :: DELTAUS

!  stream directions

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  Layer transmittances in stream directions

      REAL(FPK), INTENT(IN) :: TRANS_USERM  ( MAX_USER_STREAMS )

!  Particular solution attributes

      INTEGER, INTENT(IN) ::          PICUTOFF
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

!  Local

      INTEGER :: I, UM

!  Upwelling
!  ---------

      IF ( DO_UPWELLING ) THEN

        CALL BEAMMULT_UP_1 &
            ( N, N_USER_STREAMS, USER_STREAMS, FGSMALL, &
              STERM_LAYERMASK_UP, TRANS_USERM, DELTAUS, &
              FLIPPER, PICUTOFF, DELTRANS, ASOURCE, &
              EMULT_UP )

        CALL RAMAN_SOURCETERM_UP_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, FGSMALL, &
            N, STERM_LAYERMASK_UP, NSTREAMS, N_USER_STREAMS, &
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION, &
            TRANS_USERM, DELTAUS, USER_STREAMS, &
            EMULT_UP, ATERM_SAVE, BTERM_SAVE, &
            QSOURCE_UPOS1, ASOURCE, DELTRANS, INITRANS, &
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG, &
            WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV )

      ENDIF

!  Downwelling
!  -----------

      IF ( DO_DNWELLING ) THEN

        CALL BEAMMULT_DN_1 &
            ( N, N_USER_STREAMS, USER_STREAMS, FGSMALL, &
              STERM_LAYERMASK_DN, TRANS_USERM, DELTAUS, &
              FLIPPER, PICUTOFF, DELTRANS, ASOURCE, &
              EMULT_DN )

        CALL RAMAN_SOURCETERM_DN_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, FGSMALL, &
            N, STERM_LAYERMASK_DN, NSTREAMS, N_USER_STREAMS, &
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION, &
            TRANS_USERM, DELTAUS, USER_STREAMS, &
            EMULT_DN, ATERM_SAVE, BTERM_SAVE, &
            QSOURCE_UNEG1, ASOURCE, DELTRANS, INITRANS, &
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG, &
            WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE POSTPROCESSING_1

!

      SUBROUTINE HOMOGMULT_1 &
          ( NSTREAMS, N, N_USER_STREAMS, FGSMALL,          &
            UTRANS, USER_STREAMS, KEIGEN, KTRANS, DELTAUS, &
            HMULT_1, HMULT_2, ZETA_P, ZETA_M )

!  Whole-layer INTEGRATED homogeneous solution MULTIPLIERS
!      upwelling and downwelling for the complete atmosphere

!   @@@ RobFix 5/5/11. Small numbers analysis:
!                      added FGSMALL, DELTAUS to Argument list.

!  include file of dimensions and numbers

      USE LRRS_PARS
      USE LRRS_SMALLNOS

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  inputs

      INTEGER  , INTENT(IN) :: N_USER_STREAMS, N, NSTREAMS

!  Small number control
!   @@@ RobFix 5/5/11. Small numbers analysis added

      REAL(FPK), INTENT(IN) :: FGSMALL

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
!          HMULT_1(AA,UM,N) = SM * THETA_1 * ZETA_M(AA,UM,N)
          IF ( DABS(RHO_M) .LT. FGSMALL ) THEN
             CALL  LIMIT_GCFUNC ( -RHO_M, DELTAUS, UDEL, SD )
          ELSE
            SD = THETA_1 * ZETA_M(AA,UM,N)
          ENDIF
          HMULT_1(AA,UM,N) = SM * SD
!   @@@ RobFix 5/5/11. End of new section

          HMULT_2(AA,UM,N) = SM * THETA_2 * ZETA_P(AA,UM,N)
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE HOMOGMULT_1

!

      SUBROUTINE BEAMMULT_UP_1 &
            ( N, N_USER_STREAMS, USER_STREAMS, FGSMALL, &
              STERM_LAYERMASK_UP, TRANS_USERM, DELTAUS, &
              FLIPPER, PICUTOFF, DELTRANS, ASOURCE, &
              EMULT_UP )

!  Source term primary multiplier. Upwelling.
!   - Small numbers analysis added 7 December 2006

!  include file of dimensions and numbers

      USE LRRS_PARS
      USE LRRS_SMALLNOS

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  control integers

      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N

!  User streams, layer transmittances

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      LOGICAL, INTENT(IN) ::   STERM_LAYERMASK_UP
      REAL(FPK), INTENT(IN) :: TRANS_USERM  ( MAX_USER_STREAMS )

!  Small number control

      REAL(FPK), INTENT(IN) :: FGSMALL

!  optical thickness for use with the small number expansion

      REAL(FPK), INTENT(IN) :: DELTAUS

!  Solution control

      LOGICAL, INTENT(IN) ::   FLIPPER
      INTEGER, INTENT(IN) ::   PICUTOFF
      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: ASOURCE

!  output = multiplier for one solution source

      REAL(FPK), INTENT(OUT) :: EMULT_UP ( MAX_USER_STREAMS )

!  local

      INTEGER ::   UM
      REAL(FPK) :: SM, UDEL, WUDEL, SIGMA_M, SIGMA_P, SD, SU

!  Code

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
              IF(DABS(SIGMA_M).LT.FGSMALL)THEN
                CALL LIMIT_GCFUNC ( SIGMA_M, DELTAUS, UDEL, SD )
!         write(52,*)N,UM,FLIPPER,ABS(SIGMA_M),SD
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
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BEAMMULT_UP_1

!

      SUBROUTINE BEAMMULT_DN_1 &
            ( N, N_USER_STREAMS, USER_STREAMS, FGSMALL, &
              STERM_LAYERMASK_DN, TRANS_USERM, DELTAUS, &
              FLIPPER, PICUTOFF, DELTRANS, ASOURCE, &
              EMULT_DN )

!  Source term primary multiplier. Downwelling.
!   - Small numbers analysis added 7 December 2006

!  include file of dimensions and numbers

      USE LRRS_PARS
      USE LRRS_SMALLNOS

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Control integers

      INTEGER, INTENT(IN) ::          N_USER_STREAMS
      INTEGER, INTENT(IN) ::          N

!  User streams, layer transmittances

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      LOGICAL, INTENT(IN) ::          STERM_LAYERMASK_DN
      REAL(FPK), INTENT(IN) :: TRANS_USERM  ( MAX_USER_STREAMS )

!  Small number control

      REAL(FPK), INTENT(IN) :: FGSMALL

!  optical thickness for use with the small number expansion

      REAL(FPK), INTENT(IN) :: DELTAUS

!  Solution control

      LOGICAL, INTENT(IN) ::          FLIPPER
      INTEGER, INTENT(IN) ::          PICUTOFF
      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: ASOURCE

!  output = single multiplier for one solution source

      REAL(FPK), INTENT(OUT) :: EMULT_DN ( MAX_USER_STREAMS )

!  local variables

      INTEGER ::          UM
      REAL(FPK) :: SM, UDEL, WUDEL, SIGMA_M, SIGMA_P, SD, SU

!  Code

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
              IF(DABS(SIGMA_M).LT.FGSMALL)THEN
                CALL LIMIT_GCFUNC ( SIGMA_M, DELTAUS, UDEL, SD )
!         write(53,*)N,UM,DABS(SIGMA_M),DELTAUS, UDEL,SD
              ELSE
                SD = ( UDEL - DELTRANS ) / SIGMA_M
              ENDIF
              EMULT_DN(UM) = SM * SD
            ENDIF
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BEAMMULT_DN_1
!

      SUBROUTINE TOASOURCE &
          ( N_USER_STREAMS, &
            TOA_SOURCE )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER, INTENT(IN) ::          N_USER_STREAMS

      REAL(FPK), INTENT(OUT) :: TOA_SOURCE(MAX_USER_STREAMS)

!  local variables

      INTEGER ::          UM

!  initialise TOA source function

      DO UM = 1, N_USER_STREAMS
        TOA_SOURCE(UM) = ZERO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE TOASOURCE

!

      SUBROUTINE BOASOURCE &
          ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
            DO_LAMBERTIAN_SURFACE, FOURIER, SURFACE_FACTOR, FL1, FL2, &
            NSTREAMS, NLAYERS, N_USER_STREAMS, &
            WLOWER, LCON_XVEC, MCON_XVEC, KTRANS, &
            QUAD_STRMWGHT, USER_BIREFLEC, USER_DIRECT_BEAM, &
            IDOWNSURF, BOA_SOURCE, DIRECT_BOA_SOURCE )

!  Bottom of the atmosphere source term

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL, INTENT(IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_SURFACE

!  surface multipliers and Fourier index

      REAL(FPK), INTENT(IN) :: SURFACE_FACTOR, FL1, FL2
      INTEGER, INTENT(IN) ::          FOURIER

!  Control numbers

      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          NLAYERS
      INTEGER, INTENT(IN) ::          N_USER_STREAMS

!  Variables for integrating the downwelling surface radiance

      REAL(FPK), INTENT(IN) :: &
          LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  Quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT ( MAX_STREAMS )

!  Bidirectional functions

      REAL(FPK), INTENT(IN) :: USER_BIREFLEC ( MAX_USER_STREAMS, MAX_STREAMS )

!  Direct beam contribution

      REAL(FPK), INTENT(IN) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS )

!  output variables, direct and diffuse surface fields.

      REAL(FPK), INTENT(OUT) :: IDOWNSURF ( MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: BOA_SOURCE ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, J, I, UM, AA
      REAL(FPK) :: PAR, HOM, REFLEC, FRJ, FRAC1, FRAC2

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
            HOM = HOM + LCON_XVEC(I,AA,N) * KTRANS(AA,N) + &
                        MCON_XVEC(I,AA,N)
          ENDDO
          IDOWNSURF(I) = QUAD_STRMWGHT(I) * ( PAR + HOM )
        ENDDO

!  reflected multiple scatter intensity at user defined-angles

        IF ( DO_LAMBERTIAN_SURFACE ) THEN
         IF ( FOURIER .EQ. 0 ) THEN
          FRAC1 = FL1 * SURFACE_FACTOR
          REFLEC = ZERO
          DO J = 1, NSTREAMS
            REFLEC = REFLEC + IDOWNSURF(J)
          ENDDO
          REFLEC = REFLEC * FRAC1
          DO UM = 1, N_USER_STREAMS
            BOA_SOURCE(UM) = REFLEC
          ENDDO
         ENDIF
        ELSE
          FRAC2 = FL2 * SURFACE_FACTOR
          FRAC1 = FL1 * SURFACE_FACTOR
          DO UM = 1, N_USER_STREAMS
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              FRJ = FRAC1 + USER_BIREFLEC(UM,J) * FRAC2
              REFLEC = REFLEC + IDOWNSURF(J) * FRJ
            ENDDO
            BOA_SOURCE(UM) = REFLEC
          ENDDO
        ENDIF

!  Add direct beam if flagged

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          DO UM = 1, N_USER_STREAMS
            DIRECT_BOA_SOURCE(UM) = USER_DIRECT_BEAM(UM)
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BOASOURCE

!

      SUBROUTINE RAMAN_SOURCETERM_UP_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, LRRS_FGSMALL, &
            N, STERM_LAYERMASK_UP, NSTREAMS, N_USER_STREAMS, &
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION, &
            TRANS_USERM, DELTAUS, USER_STREAMS, &
            EMULT_UP, ATERM_SAVE, BTERM_SAVE, &
            QSOURCE_UPOS1, ASOURCE, DELTRANS, INITRANS, &
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG, &
            WLAYER_PIST_UP, SUSAV, SDSAV )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS
      USE LRRS_SMALLNOS

      IMPLICIT NONE

!  input
!  -----

!  Solution number

      INTEGER, INTENT(IN) ::   SOLUTION

!  number of streams and solutions

      INTEGER, INTENT(IN) ::   NSTREAMS, N_USER_STREAMS

!  Number of layer, and layer control

      INTEGER, INTENT(IN) ::   N
      LOGICAL, INTENT(IN) ::   STERM_LAYERMASK_UP

!  Criterion for small numbers analysis

      REAL(FPK), INTENT(IN) :: LRRS_FGSMALL

!  Fourier component

      INTEGER, INTENT(IN) ::   FOURIER

!  Flags for multiple scattering and single scatter correction to elasti

      LOGICAL, INTENT(IN) ::   DO_MSMODE_LIDORT
      LOGICAL, INTENT(IN) ::   DO_SSCORRECTION

!  Control for the particular solutions

      LOGICAL, INTENT(IN) ::   CHANGE_EXPONENT
      LOGICAL, INTENT(IN) ::   FLIPPER
      INTEGER, INTENT(IN) ::   PICUTOFF

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

!  Eigenvalues and transmittance factors.

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )

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

      INTEGER ::   AA, UM, M
      REAL(FPK) :: PAR, RHOM, GP, PAR1, MULT, DMU, UDEL, KVAL

!  Local multiplier arrays

      REAL(FPK) :: SU ( MAX_STREAMS ), SD ( MAX_STREAMS )

!  ***************************************************************
!  Single RRS constribution: Greens function multipliers
!  New code programmed with small number analysis and FLIPPER
!   27 September 2006, 29 November 2006. R. Spurr, RT SOLUTIONS.
!  ***************************************************************

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

!  local quantities

          DMU    = ONE/USER_STREAMS(UM)
          UDEL   = TRANS_USERM(UM)

!  Only require a new multiplier when the exponent changes

          IF ( CHANGE_EXPONENT ) THEN

!  set directional multipliers

            MULT   = EMULT_UP(UM)

!  Zero the Greens function multipliers if no solution

            IF ( N .GT. PICUTOFF ) THEN

              DO AA = 1, NSTREAMS
                SDSAV(AA,UM) = ZERO
                SUSAV(AA,UM) = ZERO
              ENDDO

!  Compute multipliers. Debug for checking (30 November 2006):
!       if ( dabs(rhom).lt.1.0d-06)write(77,*)fourier,n,aa,um,p

            ELSE

              IF ( FLIPPER ) THEN
                DO AA = 1, NSTREAMS
                  KVAL = KEIGEN(AA,N)
                  RHOM = ASOURCE - KVAL
                  IF ( DABS(RHOM) .LT. LRRS_FGSMALL )THEN
                    call limit_ppgmult_wd ( rhom, dmu, one, &
                      kval, deltaus, ktrans(aa,n), udel, susav(aa,um))
!if(RHOm.ne.0.0d0)write(54,'(3i3,1p3e20.12,f12.6)')&
!   n,um,aa,rhom,susav(aa,um),( HMULT_1(AA,UM,N) - MULT ) / RHOM, &
!       abs(1.0d0-(susav(aa,um)*RHOM/( HMULT_1(AA,UM,N) - MULT )))*1000000.0d0
                  ELSE
                   SUSAV(AA,UM)  = ( HMULT_1(AA,UM,N) - MULT ) / RHOM
                  ENDIF
                  GP = ONE / ( ASOURCE + KVAL )
                  SDSAV(AA,UM)  = GP*(MULT-HMULT_2(AA,UM,N)*DELTRANS)
                ENDDO
              ELSE
                DO AA = 1, NSTREAMS
                  KVAL = KEIGEN(AA,N)
                  RHOM = ASOURCE - KVAL
                  IF ( DABS(RHOM) .LT. LRRS_FGSMALL ) THEN
                   call limit_ppgmult_wu ( rhom, dmu, one, &
                      kval, deltaus, ktrans(aa,n), udel, sdsav(aa,um))
!if(RHOm.ne.0.0d0)write(55,'(3i3,1p3e20.12,f12.6)')&
!   n,um,aa,rhom,sdsav(aa,um),( HMULT_2(AA,UM,N) - MULT ) / RHOM, &
!       abs(1.0d0-(sdsav(aa,um)*RHOM/( HMULT_2(AA,UM,N) - MULT )))*1000000.0d0
                  ELSE
                   SDSAV(AA,UM)  = ( HMULT_2(AA,UM,N) - MULT ) / RHOM
                  ENDIF
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

!  Summation over all streamsex

          PAR = ZERO
          DO AA = 1, NSTREAMS
            PAR = PAR + U_XPOS(UM,AA,N)*SD(AA) + U_XNEG(UM,AA,N)*SU(AA)
          ENDDO

!  Upgrade contribution

          WLAYER_PIST_UP(N,UM) = WLAYER_PIST_UP(N,UM) + PAR*INITRANS

!  Smart debug....................
!          if (n.eq.11)write(66,6)n,solution, flipper,par*initrans
!          if (n.gt.11)write(76,6)n,solution, flipper,par*initrans
! 6        format(2i4,L2,1pe20.10)

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
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, LRRS_FGSMALL, &
            N, STERM_LAYERMASK_DN, NSTREAMS, N_USER_STREAMS, &
            FOURIER, FLIPPER, CHANGE_EXPONENT, PICUTOFF, SOLUTION, &
            TRANS_USERM, DELTAUS, USER_STREAMS, &
            EMULT_DN, ATERM_SAVE, BTERM_SAVE, &
            QSOURCE_UNEG1, ASOURCE, DELTRANS, INITRANS, &
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG, &
            WLAYER_PIST_DN, SUSAV, SDSAV )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS
      USE LRRS_SMALLNOS

      IMPLICIT NONE

!  arguments
!  ---------

!  Solution number

      INTEGER, INTENT(IN) ::   SOLUTION

!  number of streams and solutions

      INTEGER, INTENT(IN) ::   NSTREAMS, N_USER_STREAMS

!  Number of layer, and layer control

      INTEGER, INTENT(IN) ::   N
      LOGICAL, INTENT(IN) ::   STERM_LAYERMASK_DN

!  Criterion for small numbers analysis

      REAL(FPK), INTENT(IN) :: LRRS_FGSMALL

!  Fourier component

      INTEGER, INTENT(IN) ::   FOURIER

!  Flags for multiple scattering and single scatter correction to elasti

      LOGICAL, INTENT(IN) ::   DO_MSMODE_LIDORT
      LOGICAL, INTENT(IN) ::   DO_SSCORRECTION

!  Control for the particular solutions

      LOGICAL, INTENT(IN) ::   CHANGE_EXPONENT
      LOGICAL, INTENT(IN) ::   FLIPPER
      INTEGER, INTENT(IN) ::   PICUTOFF

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

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )
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

      INTEGER ::   AA, UM, M
      REAL(FPK) :: PAR, RHOM, GP, MULT, DMU, UDEL, KVAL

!  Local multiplier arrays

      REAL(FPK) :: SU ( MAX_STREAMS ), SD ( MAX_STREAMS )

!  ***************************************************************
!  Loop over all the solutions to get Greens function multipliers
!  New code programmed with small number analysis and FLIPPER
!   27 September 2006, 29 November 2006. R. Spurr, RT SOLUTIONS.
!  ***************************************************************

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
                  RHOM = ASOURCE - KVAL
                  IF ( DABS(RHOM) .LT. LRRS_FGSMALL )THEN
                    call limit_ppgmult_wu ( rhom, dmu, one, &
                      kval, deltaus, ktrans(aa,n), udel, susav(aa,um))
!if(RHOm.ne.0.0d0)write(56,'(3i3,1p3e20.12,f12.6)')&
!   n,um,aa,rhom,susav(aa,um),( HMULT_2(AA,UM,N) - MULT ) / RHOM, &
!       abs(1.0d0-(susav(aa,um)*RHOM/( HMULT_2(AA,UM,N) - MULT )))*1000000.0d0
                  ELSE
                    SUSAV(AA,UM)  = ( HMULT_2(AA,UM,N) - MULT ) / RHOM
                  ENDIF
                  GP = ONE / ( ASOURCE + KVAL )
                  SDSAV(AA,UM)  = GP*(MULT-HMULT_1(AA,UM,N)*DELTRANS)
                ENDDO
              ELSE
                DO AA = 1, NSTREAMS
                  KVAL =  KEIGEN(AA,N)
                  RHOM = ASOURCE - KVAL
                  IF ( DABS(RHOM) .LT. LRRS_FGSMALL ) THEN
                    call limit_ppgmult_wd ( rhom, dmu, one, &
                      kval, deltaus, ktrans(aa,n), udel, sdsav(aa,um))
!if(RHOm.ne.0.0d0)write(57,'(3i3,1p3e20.12,f12.6)')&
!   n,um,aa,rhom,sdsav(aa,um),( HMULT_1(AA,UM,N) - MULT ) / RHOM, &
!       abs(1.0d0-(sdsav(aa,um)*RHOM/( HMULT_1(AA,UM,N) - MULT )))*1000000.0d0
                  ELSE
                    SDSAV(AA,UM)  = ( HMULT_1(AA,UM,N) - MULT ) / RHOM
                  ENDIF
                  GP = ONE / ( ASOURCE + KVAL )
                  SUSAV(AA,UM)  = GP*(MULT-HMULT_2(AA,UM,N)*DELTRANS)
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


      END MODULE lrrs_postprocessing_1

