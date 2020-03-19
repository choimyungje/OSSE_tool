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
! #                                                             #
! #             ELASTIC_INTENSITY_UP_1                          #
! #             ELASTIC_INTENSITY_DN_1                          #
! #                                                             #
! #             ELASTIC_SOURCETERM_UP_1                         #
! #             ELASTIC_SOURCETERM_DN_1                         #
! #                                                             #
! #  These subroutines calculate the upwelling and downwelling  #
! #    post-processed radiance field, using the recursion       #
! #    equation based on transmittance up/down from BOA/TOA     #
! #    and the completion of the whole and partial-layer        #
! #    source functions (STERM_UP, STERM_DN). Source-function   #
! #    integration is based on the use of particular integrals  #
! #    determined by the classical substitution method. The     #
! #    Green's function method is used only for the Inelastic   #
! #    RTE equation.                                            #
! #                                                             #
! ###############################################################


      MODULE lrrs_elastic_intensity_1

      PRIVATE
      PUBLIC :: ELASTIC_INTENSITY_UP_1,&
                ELASTIC_INTENSITY_DN_1

      CONTAINS

      SUBROUTINE ELASTIC_INTENSITY_UP_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, DO_ELASTIC_SCATTERING, &
            DO_QUAD_RESULTS,  DO_INCLUDE_MVOUT,DO_INCLUDE_SURFACE, &
            DO_INCLUDE_DIRECTBEAM, FOURIER, DO_LAMBERTIAN_SURFACE, &
            SURFACE_FACTOR, FL1, FL2, USER_BIREFLEC, FLUXMULT, &
            NSTREAMS, NLAYERS, N_USER_STREAMS, N_LOUTPUT, LEVELMASK_UP, &
            QUAD_STRMWGHT, USER_DIRECT_BEAM, &
            WUPPER, WLOWER, KTRANS, TRANS_USERM, &
            LCON, MCON, LCON_XVEC, MCON_XVEC, &
            U_XPOS,  U_XNEG, U_WPOS1, U_WPOS2, &
            HMULT_1,     HMULT_2,      EMULT_UP, &
            IDOWNSURF,    CUMSOURCE_UP, &
            ELASTIC_F_UP, QUADELASTIC_F_UP )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS
      USE LRRS_POSTPROCESSING_1
      USE LRRS_RAMAN_INTENSITY_1

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL, INTENT(IN) ::          DO_MSMODE_LIDORT
      LOGICAL, INTENT(IN) ::          DO_SSCORRECTION

!  output for quadrature-angle radiances

      LOGICAL, INTENT(IN) ::          DO_QUAD_RESULTS
      LOGICAL, INTENT(IN) ::          DO_INCLUDE_MVOUT

!  Local surface control

      LOGICAL, INTENT(IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_SURFACE

!  elastic scattering flag. If not set for any layers with Fourier > 2,
!  then there will be no layer scattering source terms (source function
!  integration is then trivial)

      LOGICAL, INTENT(IN) ::          DO_ELASTIC_SCATTERING(MAX_LAYERS)

!  Flux multiplier and Fourier component

      REAL(FPK), INTENT(IN) :: FLUXMULT
      INTEGER, INTENT(IN) ::          FOURIER

!  Surface factor = 1 for Fourier 0, otherwise = 0.5

      REAL(FPK), INTENT(IN) :: SURFACE_FACTOR, FL1, FL2

!  Basic control numbers

      INTEGER, INTENT(IN) ::          NSTREAMS, NLAYERS
      INTEGER, INTENT(IN) ::          N_USER_STREAMS
      INTEGER, INTENT(IN) ::          N_LOUTPUT

!  optical depth output control

      INTEGER, INTENT(IN) ::          LEVELMASK_UP    ( MAX_LOUTPUT )

!  Quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT ( MAX_STREAMS )

!  Bidirectional functions

      REAL(FPK), INTENT(IN) :: USER_BIREFLEC ( MAX_USER_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS )

!  Variables for the particular integral (quadrature)

      REAL(FPK), INTENT(IN) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: WUPPER ( MAX_2_STREAMS, MAX_LAYERS )

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  Discrete ordinate solution variables
!  ------------------------------------

!  integration constants

      REAL(FPK), INTENT(IN) :: &
            LCON ( MAX_STREAMS, MAX_LAYERS ), &
            MCON ( MAX_STREAMS, MAX_LAYERS )

!  homogeneous vectors with integration constants

      REAL(FPK), INTENT(IN) :: &
          LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  homogeneous solutions at user angles

      REAL(FPK), INTENT(IN) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  classical particular integral solutions at user angles

      REAL(FPK), INTENT(IN) :: &
            U_WPOS1 ( MAX_USER_STREAMS, MAX_LAYERS ), &
            U_WPOS2 ( MAX_USER_STREAMS, MAX_LAYERS )

!  whole layer transmittances

      REAL(FPK), INTENT(IN) :: TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS )

!  Multiplier input
!  ----------------

!  Whole layer multipliers

      REAL(FPK), INTENT(IN) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: &
            EMULT_UP      ( MAX_USER_STREAMS, MAX_LAYERS )

!  output
!  ------

!  Output diffuse BOA source term

      REAL(FPK), INTENT(OUT) :: IDOWNSURF ( MAX_STREAMS )

!  Cumulative radiance source term

      REAL(FPK), INTENT(OUT) :: CUMSOURCE_UP ( MAX_USER_STREAMS, 0:MAX_LAYERS )

!  elastic scattering Fourier-component fields

      REAL(FPK), INTENT(OUT) :: &
             ELASTIC_F_UP     ( MAX_LOUTPUT, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: &
             QUADELASTIC_F_UP ( MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

!  local help

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER ::          UTA, UM
      LOGICAL ::          SOURCETERM_FLAG

!  local BOA source terms

      REAL(FPK) :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(FPK) :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )

!  Local layer source termes (total and diffuse-only)

      REAL(FPK) :: LAYER_SOURCE      ( MAX_USER_STREAMS )
      REAL(FPK) :: MSCAT_LAYERSOURCE ( MAX_USER_STREAMS )

!  Initialize recursion.
!  ---------------------

!  Get the radiance field at the surface
!    Diffuse and direct contributions. This module is Lambertian-only

      CALL BOASOURCE &
          ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, &
            DO_LAMBERTIAN_SURFACE, FOURIER, SURFACE_FACTOR, FL1, FL2, &
            NSTREAMS, NLAYERS, N_USER_STREAMS, &
            WLOWER, LCON_XVEC, MCON_XVEC, KTRANS, &
            QUAD_STRMWGHT, USER_BIREFLEC, USER_DIRECT_BEAM, &
            IDOWNSURF, BOA_SOURCE, DIRECT_BOA_SOURCE )

!  start the recurrence with this value

      DO UM = 1, N_USER_STREAMS
        CUMSOURCE_UP(UM,0) = BOA_SOURCE(UM) + DIRECT_BOA_SOURCE(UM)
      ENDDO

!  initialise cumulative source term loop

      NC = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1
      NUT = 0

!  loop over all output optical depths
!  -----------------------------------

!  Start loop with the lowest value first

      DO UTA = N_LOUTPUT, 1, -1

!  Level index for given optical depth

        NLEVEL = LEVELMASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        NUT = NLEVEL + 1
        DO N = NSTART, NUT, -1
          NC = NLAYERS + 1 - N
          SOURCETERM_FLAG = DO_ELASTIC_SCATTERING(N)
          CALL ELASTIC_SOURCETERM_UP_1 &
                ( DO_MSMODE_LIDORT, DO_SSCORRECTION, &
                  SOURCETERM_FLAG, N, &
                  NSTREAMS, N_USER_STREAMS, FOURIER, &
                  HMULT_1, HMULT_2, EMULT_UP, &
                  LCON, MCON, U_XPOS, U_XNEG, U_WPOS1, U_WPOS2, &
                  LAYER_SOURCE, MSCAT_LAYERSOURCE )
          DO UM = 1, N_USER_STREAMS
            CUMSOURCE_UP(UM,NC) = LAYER_SOURCE(UM) + &
                         TRANS_USERM(UM,N)*CUMSOURCE_UP(UM,NC-1)
          ENDDO
        ENDDO

!  No offgrid output this module...........!!!!!!!!!!!!!!!!!

!  Ongrid output only
!  ------------------

!  Quadrature output at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

        IF ( DO_INCLUDE_MVOUT .OR. DO_QUAD_RESULTS ) THEN
          CALL QUAD_INTENSITY_UP_1 &
                 ( NSTREAMS, NLAYERS, &
                   NLEVEL, UTA, FLUXMULT, &
                   WLOWER, WUPPER, KTRANS, LCON_XVEC, MCON_XVEC, &
                   QUADELASTIC_F_UP )
        ENDIF

!  User-defined stream output, just set to the cumulative source term

        DO UM = 1, N_USER_STREAMS
          ELASTIC_F_UP(UTA,UM) = FLUXMULT * CUMSOURCE_UP(UM,NC)
        ENDDO

!  Check for updating the recursion
!  --------------------------------

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE ELASTIC_INTENSITY_UP_1

!

      SUBROUTINE ELASTIC_INTENSITY_DN_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, DO_ELASTIC_SCATTERING, &
            DO_QUAD_RESULTS,  DO_INCLUDE_MVOUT, FOURIER, FLUXMULT, &
            NSTREAMS, N_USER_STREAMS, N_LOUTPUT, LEVELMASK_DN, &
            WLOWER, KTRANS, TRANS_USERM, &
            LCON, MCON, LCON_XVEC, MCON_XVEC, &
            U_XPOS,  U_XNEG,  U_WNEG1, U_WNEG2, &
            HMULT_1,     HMULT_2,      EMULT_DN, &
            CUMSOURCE_DN, ELASTIC_F_DN, QUADELASTIC_F_DN )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS
      USE LRRS_POSTPROCESSING_1
      USE LRRS_RAMAN_INTENSITY_1

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL, INTENT(IN) ::          DO_MSMODE_LIDORT
      LOGICAL, INTENT(IN) ::          DO_SSCORRECTION

!  output for quadrature-angle radiances

      LOGICAL, INTENT(IN) ::          DO_QUAD_RESULTS
      LOGICAL, INTENT(IN) ::          DO_INCLUDE_MVOUT

!  elastic scattering flag. If not set for any layers with Fourier > 2,
!  then there will be no layer scattering source terms (source function
!  integration is then trivial)

      LOGICAL, INTENT(IN) ::          DO_ELASTIC_SCATTERING(MAX_LAYERS)

!  Flux factor  and Fourier component

      REAL(FPK), INTENT(IN) :: FLUXMULT
      INTEGER, INTENT(IN) ::          FOURIER

!  Basic control numbers

      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          N_USER_STREAMS
      INTEGER, INTENT(IN) ::          N_LOUTPUT

!  optical depth output control

      INTEGER, INTENT(IN) ::          LEVELMASK_DN     ( MAX_LOUTPUT )

!  Variables for the particular integral (quadrature)

      REAL(FPK), INTENT(IN) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  Multiplier variables, whole layers

      REAL(FPK), INTENT(IN) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: &
            EMULT_DN   ( MAX_USER_STREAMS, MAX_LAYERS )

!  Discrete ordinate solution variables
!  ------------------------------------

!  boundary value integration constants

      REAL(FPK), INTENT(IN) :: &
            LCON ( MAX_STREAMS, MAX_LAYERS ), &
            MCON ( MAX_STREAMS, MAX_LAYERS )

!  homogeneous vectors with integration constants

      REAL(FPK), INTENT(IN) :: &
          LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  User homogeneous solutions

      REAL(FPK), INTENT(IN) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: &
            U_WNEG1 ( MAX_USER_STREAMS, MAX_LAYERS ), &
            U_WNEG2 ( MAX_USER_STREAMS, MAX_LAYERS )

!  whole layer transmittances

      REAL(FPK), INTENT(IN) :: TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS )

!  output
!  ------

!  Cumulative radiance source term

      REAL(FPK), INTENT(OUT) :: CUMSOURCE_DN ( MAX_USER_STREAMS, 0:MAX_LAYERS )

!  elastic scattering Fourier-component fields

      REAL(FPK), INTENT(OUT) :: &
             ELASTIC_F_DN ( MAX_LOUTPUT, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: &
             QUADELASTIC_F_DN ( MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER ::          UTA, UM
      LOGICAL ::          SOURCETERM_FLAG

!  TOA source (usually zero, no diffuse downwelling field at TOA)

      REAL(FPK) :: TOA_SOURCE        ( MAX_USER_STREAMS )

!  Layer source terms, total and diffuse

      REAL(FPK) :: LAYER_SOURCE      ( MAX_USER_STREAMS )
      REAL(FPK) :: MSCAT_LAYERSOURCE ( MAX_USER_STREAMS )

!  Initialize recursion
!  --------------------

!  Get the source term

      CALL TOASOURCE &
          ( N_USER_STREAMS, &
            TOA_SOURCE )

!  Start the recurrence with this term

      DO UM = 1, N_USER_STREAMS
        CUMSOURCE_DN(UM,0) = TOA_SOURCE(UM)
      ENDDO

!  initialise cumulative source term loop

      NC = 0
      NSTART = 1
      NUT_PREV = NSTART - 1
      NUT = 0

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_LOUTPUT

!  Level index for given optical depth

        NLEVEL = LEVELMASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        NUT = NLEVEL
        DO N = NSTART, NUT
          NC = N
          SOURCETERM_FLAG = DO_ELASTIC_SCATTERING(N)
          CALL ELASTIC_SOURCETERM_DN_1 &
                ( DO_MSMODE_LIDORT, DO_SSCORRECTION, &
                  SOURCETERM_FLAG, N, &
                  NSTREAMS, N_USER_STREAMS, FOURIER, &
                  HMULT_1, HMULT_2, EMULT_DN, &
                  LCON, MCON, U_XPOS, U_XNEG, U_WNEG1, U_WNEG2, &
                  LAYER_SOURCE, MSCAT_LAYERSOURCE )

          DO UM = 1, N_USER_STREAMS
            CUMSOURCE_DN(UM,NC) = LAYER_SOURCE(UM) + &
                            TRANS_USERM(UM,N)*CUMSOURCE_DN(UM,NC-1)
          ENDDO
        ENDDO

!  No offgrid output this module............!!!!!!!!!!!!!

!  Ongrid output only
!  ------------------

!  Quadrature results at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

        IF ( DO_INCLUDE_MVOUT .OR. DO_QUAD_RESULTS ) THEN
          CALL QUAD_INTENSITY_DN_1 &
                 ( NSTREAMS, NLEVEL, UTA, FLUXMULT, &
                   WLOWER, KTRANS, LCON_XVEC, MCON_XVEC, &
                   QUADELASTIC_F_DN )
        ENDIF

!  User-defined stream output, just set to the cumulative source term

        DO UM = 1, N_USER_STREAMS
          ELASTIC_F_DN(UTA,UM) = FLUXMULT * CUMSOURCE_DN(UM,NC)
        ENDDO

!  Check for updating the recursion
!  --------------------------------

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE ELASTIC_INTENSITY_DN_1

!

      SUBROUTINE ELASTIC_SOURCETERM_UP_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, &
            SOURCETERM_FLAG, N, &
            NSTREAMS, N_USER_STREAMS, FOURIER, &
            HMULT_1, HMULT_2, EMULT_UP, &
            LCON, MCON, U_XPOS, U_XNEG, U_WPOS1, U_WPOS2, &
            LAYERSOURCE, MSCAT_LAYERSOURCE )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  arguments
!  ---------

!  input flags and control

      INTEGER, INTENT(IN) ::          N
      INTEGER, INTENT(IN) ::          NSTREAMS, N_USER_STREAMS, FOURIER
      LOGICAL, INTENT(IN) ::          DO_MSMODE_LIDORT, DO_SSCORRECTION
      LOGICAL, INTENT(IN) ::          SOURCETERM_FLAG

!  Input multipliers

      REAL(FPK), INTENT(IN) :: &
            EMULT_UP ( MAX_USER_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Input discrete ordinate solution variables

      REAL(FPK), INTENT(IN) :: &
            LCON ( MAX_STREAMS, MAX_LAYERS ), &
            MCON ( MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: &
            U_WPOS1 ( MAX_USER_STREAMS, MAX_LAYERS ), &
            U_WPOS2 ( MAX_USER_STREAMS, MAX_LAYERS )

!  output layer source terms

      REAL(FPK), INTENT(OUT) :: LAYERSOURCE       ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: MSCAT_LAYERSOURCE ( MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER ::   AA, UM, M
      REAL(FPK) :: SHOM, SFOR2, SFOR1
      REAL(FPK) :: LCON_UXVEC ( MAX_USER_STREAMS, MAX_STREAMS )
      REAL(FPK) :: MCON_UXVEC ( MAX_USER_STREAMS, MAX_STREAMS )

!  No layer source term if no scattering in the layer

      IF ( .NOT. SOURCETERM_FLAG ) THEN
        DO UM = 1, N_USER_STREAMS
          LAYERSOURCE(UM)       = ZERO
        ENDDO
        RETURN
      ENDIF

!  Fourier not used

      M = FOURIER

!  Local quantities

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
        ENDDO
      ENDDO

!  Whole layer source function ( Classical solution )

      DO UM = 1, N_USER_STREAMS
        SFOR2 =  EMULT_UP(UM,N) * U_WPOS2(UM,N)
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N) + &
                        MCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N)
       ENDDO
        MSCAT_LAYERSOURCE(UM) = SFOR2 + SHOM
      ENDDO

!  Add single scatter contribution if flagged.

!  .. If operating in Ms-mode only or doing the exact single scatter
!        then copy multiple scatter term
!  .. Otherwise, Full radiance mode, add single scatter part

!  Debug 18 Nov 2005, LCON is NAN for Fourier = 1
!        if ( fourier.eq.1) write(*,*)
!     &     n,do_msmode_lidort,lcon(1,n),u_xpos(1,1,n),u_xneg(1,1,n),
!     &     u_wpos1(1,n),u_wpos2(1,n),emult_up(1,n)

      IF ( DO_MSMODE_LIDORT .OR. DO_SSCORRECTION ) THEN
        DO UM = 1, N_USER_STREAMS
          LAYERSOURCE(UM) = MSCAT_LAYERSOURCE(UM)
        ENDDO
      ELSE
        DO UM = 1, N_USER_STREAMS
          SFOR1 = U_WPOS1(UM,N) * EMULT_UP(UM,N)
          LAYERSOURCE(UM) = MSCAT_LAYERSOURCE(UM) + SFOR1
       ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE ELASTIC_SOURCETERM_UP_1

!

      SUBROUTINE ELASTIC_SOURCETERM_DN_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, &
            SOURCETERM_FLAG, N, &
            NSTREAMS, N_USER_STREAMS, FOURIER, &
            HMULT_1, HMULT_2, EMULT_DN, &
            LCON, MCON, U_XPOS, U_XNEG, U_WNEG1, U_WNEG2, &
            LAYERSOURCE, MSCAT_LAYERSOURCE )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  arguments
!  ---------

!  input flags and control

      INTEGER, INTENT(IN) ::          N
      INTEGER, INTENT(IN) ::          NSTREAMS, N_USER_STREAMS, FOURIER
      LOGICAL, INTENT(IN) ::          DO_MSMODE_LIDORT, DO_SSCORRECTION
      LOGICAL, INTENT(IN) ::          SOURCETERM_FLAG

!  Input Multipliers

      REAL(FPK), INTENT(IN) :: &
            EMULT_DN ( MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  input discrete ordinate solution variables

      REAL(FPK), INTENT(IN) :: &
            LCON ( MAX_STREAMS, MAX_LAYERS ), &
            MCON ( MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: &
            U_WNEG1 ( MAX_USER_STREAMS, MAX_LAYERS ), &
            U_WNEG2 ( MAX_USER_STREAMS, MAX_LAYERS )

!  output layer source terms

      REAL(FPK), INTENT(OUT) :: LAYERSOURCE       ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: MSCAT_LAYERSOURCE ( MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          AA, UM, M
      REAL(FPK) :: SHOM, SFOR2, SFOR1
      REAL(FPK) :: LCON_UXVEC ( MAX_USER_STREAMS, MAX_STREAMS )
      REAL(FPK) :: MCON_UXVEC ( MAX_USER_STREAMS, MAX_STREAMS )

!  No layer source term if no scattering in the layer

      IF ( .NOT. SOURCETERM_FLAG ) THEN
        DO UM = 1, N_USER_STREAMS
          LAYERSOURCE(UM)       = ZERO
        ENDDO
        RETURN
      ENDIF

!  Fourier not used

      M = FOURIER

!  Local quantities

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
        ENDDO
      ENDDO

!  Whole layer source function ( Classical solution )

      DO UM = 1, N_USER_STREAMS
        SFOR2 =  EMULT_DN(UM,N) * U_WNEG2(UM,N)
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N) + &
                        MCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N)
        ENDDO
        MSCAT_LAYERSOURCE(UM) = SFOR2 + SHOM
      ENDDO

!  Add single scatter contribution if flagged.

!  .. If operating in Ms-mode only or doing the exact single scatter
!        then copy multiple scatter term
!  .. Otherwise, Full radiance mode, add single scatter part

      IF ( DO_MSMODE_LIDORT .OR. DO_SSCORRECTION ) THEN
        DO UM = 1, N_USER_STREAMS
          LAYERSOURCE(UM) = MSCAT_LAYERSOURCE(UM)
        ENDDO
      ELSE
        DO UM = 1, N_USER_STREAMS
          SFOR1 = U_WNEG1(UM,N) * EMULT_DN(UM,N)
          LAYERSOURCE(UM) = MSCAT_LAYERSOURCE(UM) + SFOR1
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE ELASTIC_SOURCETERM_DN_1


      END MODULE lrrs_elastic_intensity_1

