
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
! #   SUBROUTINES (PUBLIC):                                     #
! #                                                             #
! #             ELASTIC_JACOBIAN_UP_1                           #
! #             ELASTIC_JACOBIAN_DN_1                           #
! #                                                             #
! #             ELASTIC_SURFJAC_UP_1                            #
! #             ELASTIC_SURFJAC_DN_1                            #
! #                                                             #
! #   SUBROUTINES (PRIVATE, called by ELASTIC_JACOBIAN*):       #
! #                                                             #
! #             L_ELASTIC_LSST_UP                               #
! #             L_ELASTIC_LSST_DN                               #
! #                                                             #
! #  These subroutines calculate the upwelling and downwelling  #
! #    post-processed Jacobian field, using the recursion       #
! #    equation based on transmittance up/down from BOA/TOA     #
! #    and the completion of the whole layer source functions.  #
! #    Integration is based on the use of particular integrals  #
! #    determined by the classical substitution method. The     #
! #    Green's function method is used only for the Inelastic   #
! #    RTE equation.                                            #
! #                                                             #
! #   Version 2.1 with BRDF setup, November 2007.               #
! #   Version 2.2 with Column  Linearization, Nov 2008.         #
! #   Version 2.2 with Profile Linearization, Mar 2009.         #
! #   Version 2.5 Surface supplements, October 2015.            #
! #                                                             #
! ###############################################################


      MODULE lrrs_L_elastic_jacobians_1_m

      USE LRRS_PARS_m, Only : LDU

      PRIVATE :: L_ELASTIC_LSST_UP, &
                 L_ELASTIC_LSST_DN

      PUBLIC  :: ELASTIC_JACOBIAN_UP_1, &
                 ELASTIC_JACOBIAN_DN_1, &
                 ELASTIC_SURFJAC_UP_1,  &
                 ELASTIC_SURFJAC_DN_1

      CONTAINS

      SUBROUTINE ELASTIC_JACOBIAN_UP_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, DO_ELASTIC_SCATTERING, & ! Inputs
            DO_QUAD_RESULTS, DO_INCLUDE_MVOUT, DO_INCLUDE_DIRECTBEAM, & ! Inputs
            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR,      & ! Inputs
            FOURIER, NSTREAMS, NLAYERS, N_USER_STREAMS, N_LOUTPUT,    & ! Inputs
            Raman_Idx, Brdf_Idx, ALBEDOS_RANKED, USER_BRDF_F, FMULT,  & ! Inputs
            QUAD_STRMWGHT, LEVELMASK_UP, KVARY, K, KPARS, NKSTORAGE,  & ! Inputs
            KTRANS, TRANS_USERM, LCON, MCON, XPOS, XNEG,      & ! Inputs
            U_XPOS, U_XNEG, U_WPOS1, U_WPOS2,                 & ! Inputs
            HMULT_1, HMULT_2, EMULT_UP, CUMSOURCE_UP,         & ! Inputs
            L_WUPPER, L_WLOWER, L_KTRANS, L_USER_DIRECT_BEAM, & ! Inputs
            L_TRANS_USERM, NCON, PCON, L_XPOS, L_XNEG,        & ! Inputs
            L_U_XPOS,  L_U_XNEG, L_U_WPOS1, L_U_WPOS2,        & ! Inputs
            L_HMULT_1,     L_HMULT_2,      L_EMULT_UP,        & ! Inputs
            L_ELASTIC_F_UP, L_QUADELASTIC_F_UP )                ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_LAYERS_NK, MAX_STREAMS, MAX_2_STREAMS, MAX_MOMENTS, & 
                              MAX_USER_STREAMS, MAX_LOUTPUT, MAX_POINTS, MAX_POINTS, MAX_ATMOSWFS

!  Use modules

      USE LRRS_L_POSTPROCESSING_1_m  , Only : LPC_BOASOURCE
      USE LRRS_L_RAMAN_JACOBIANS_1_m , Only : QUAD_JACOBIAN_UP_1

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT
      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION

!  Local surface control

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE
      LOGICAL  , INTENT(IN) :: DO_INCLUDE_DIRECTBEAM
      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE

!  control output for quadrature-angle radiances + fluxes

      LOGICAL  , INTENT(IN) :: DO_QUAD_RESULTS
      LOGICAL  , INTENT(IN) :: DO_INCLUDE_MVOUT

!  elastic scattering flag. If not set for any layers with Fourier > 2,
!  then there will be no layer scattering source terms (source function
!  integration is then trivial)

      LOGICAL  , INTENT(IN) :: DO_ELASTIC_SCATTERING(MAX_LAYERS)

!  Flux multiplier and Fourier component

      REAL(FPK), INTENT(IN) :: FMULT
      INTEGER  , INTENT(IN) :: FOURIER

!  Indices, Version 2.5, 9/11/15

      INTEGER  , INTENT(IN) :: Raman_IDX, Brdf_IDX

!  2-deltam0 factor

      REAL(fpk), intent(in) :: SURFACE_FACTOR

!  Albedos

      REAL(fpk), intent(in) :: ALBEDOS_RANKED ( MAX_POINTS)

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in) :: USER_BRDF_F ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_USER_STREAMS, MAX_POINTS )

!  Basic control numbers

      INTEGER  , INTENT(IN) :: NSTREAMS, NLAYERS
      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      INTEGER  , INTENT(IN) :: N_LOUTPUT

!  Linearization control

      LOGICAL  , INTENT(IN) :: KVARY
      INTEGER  , INTENT(IN) :: K, KPARS, NKSTORAGE(MAX_LAYERS,0:MAX_LAYERS)

!  optical depth output control

      INTEGER  , INTENT(IN) :: LEVELMASK_UP ( MAX_LOUTPUT )

!  Quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT ( MAX_STREAMS )

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  Cumulative source function

      REAL(FPK), INTENT(IN) :: CUMSOURCE_UP ( MAX_USER_STREAMS, 0:MAX_LAYERS)

!  Linearizations of Direct Beam, PIs, Eigentrans

      REAL(FPK), INTENT(IN) :: L_USER_DIRECT_BEAM ( MAX_ATMOSWFS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: L_WLOWER ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_WUPPER ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS,   MAX_LAYERS )

!  Discrete ordinate solution variables
!  ------------------------------------

!  integration constants

      REAL(FPK), INTENT(IN) :: LCON ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: MCON ( MAX_STREAMS, MAX_LAYERS )

!  homogeneous vectors

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  homogeneous solutions at user angles

      REAL(FPK), INTENT(IN) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  classical particular integral solutions at user angles

      REAL(FPK), INTENT(IN) :: &
            U_WPOS1 ( MAX_USER_STREAMS, MAX_LAYERS ), &
            U_WPOS2 ( MAX_USER_STREAMS, MAX_LAYERS )

!  whole and partial layer transmittances

      REAL(FPK), INTENT(IN) :: TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS )

!  Whole layer multipliers

      REAL(FPK), INTENT(IN) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: &
            EMULT_UP      ( MAX_USER_STREAMS, MAX_LAYERS )

!  Linearizations
!  --------------

!  integration constants

      REAL(FPK), INTENT(IN) :: &
            NCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS ), &
            PCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  homogeneous vectors

      REAL(FPK), INTENT(IN) :: L_XPOS &
        ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XNEG &
        ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  homogeneous solutions at user angles

      REAL(FPK), INTENT(IN) :: L_U_XPOS &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_XNEG &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  classical particular integral solutions at user angles

      REAL(FPK), INTENT(IN) :: L_U_WPOS1 &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_WPOS2 &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS_NK )

!  whole layer transmittances

      REAL(FPK), INTENT(IN) :: L_TRANS_USERM &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS )

!  Whole layer multipliers

      REAL(FPK), INTENT(IN) :: L_HMULT_1 &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_HMULT_2 &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_EMULT_UP &
            ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS, MAX_LAYERS )

!  output = elastic scattering Fourier-component Jacobians
!  -------------------------------------------------
!mick fix 10/19/2015 - changed intent to inout
      REAL(FPK), INTENT(INOUT) :: L_ELASTIC_F_UP &
            ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS )
      REAL(FPK), INTENT(INOUT) :: L_QUADELASTIC_F_UP &
            ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

!  local help

      INTEGER   :: N, NK, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER   :: UTA, UM, Q
      LOGICAL   :: SOURCETERM_FLAG

!  local BOA source terms

      REAL(FPK) :: L_BOA_SOURCE        ( MAX_ATMOSWFS, MAX_USER_STREAMS )
      REAL(FPK) :: L_DIRECT_BOA_SOURCE ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  Local layer source termes (total and diffuse-only)

      REAL(FPK) :: L_LAYER_SOURCE      ( MAX_ATMOSWFS, MAX_USER_STREAMS )
      REAL(FPK) :: L_MSCAT_LAYERSOURCE ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  Cumulative Jacobian source term

      REAL(FPK) :: L_CUMSOURCE_UP ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  Exit if no variation

      IF ( .NOT. KVARY ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            DO UTA = 1, N_LOUTPUT
              L_ELASTIC_F_UP(Q,K,UTA,UM) = zero
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Initialize recursion.
!  ---------------------

!  Get the radiance field at the surface
!    Diffuse and direct contributions. This module is Lambertian-only

      CALL LPC_BOASOURCE &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,      & ! Inputs
            NSTREAMS, NLAYERS, N_USER_STREAMS, FOURIER, Raman_IDX, Brdf_IDX, & ! Inputs
            K, KPARS, KVARY, QUAD_STRMWGHT, LCON, MCON, XPOS, XNEG, KTRANS,  & ! Inputs
            SURFACE_FACTOR, ALBEDOS_RANKED, USER_BRDF_F, L_USER_DIRECT_BEAM, & ! Inputs
            NCON, PCON, L_XPOS, L_XNEG, L_KTRANS, L_WLOWER,                  & ! Inputs
            L_BOA_SOURCE, L_DIRECT_BOA_SOURCE )                                ! Outputs

!  start the recurrence with this value

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, KPARS
          L_CUMSOURCE_UP(Q,UM) = L_BOA_SOURCE(Q,UM) + L_DIRECT_BOA_SOURCE(Q,UM)
        ENDDO
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

!  Layer index for given optical depth

        NLEVEL = LEVELMASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        NUT = NLEVEL + 1
        DO N = NSTART, NUT, -1

          NK = NKSTORAGE(N,K)
          NC = NLAYERS + 1 - N
          SOURCETERM_FLAG = DO_ELASTIC_SCATTERING(N)

          CALL L_ELASTIC_LSST_UP &
                ( DO_MSMODE_LIDORT, DO_SSCORRECTION, SOURCETERM_FLAG,     & ! Inputs
                  N, K, KPARS, NK, NSTREAMS, N_USER_STREAMS, FOURIER,     & ! Inputs
                  HMULT_1, HMULT_2, EMULT_UP, LCON, MCON, U_XPOS, U_XNEG, & ! Inputs
                  U_WPOS1, U_WPOS2, L_HMULT_1, L_HMULT_2, L_EMULT_UP,     & ! Inputs
                  NCON, PCON, L_U_XPOS, L_U_XNEG, L_U_WPOS1, L_U_WPOS2,   & ! Inputs
                  L_LAYER_SOURCE, L_MSCAT_LAYERSOURCE )                     ! Outputs

          IF ( N.EQ.K .OR. K.EQ.0 ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, KPARS
                L_CUMSOURCE_UP(Q,UM) = L_LAYER_SOURCE(Q,UM) + &
                         TRANS_USERM(UM,N)   * L_CUMSOURCE_UP(Q,UM) + &
                       L_TRANS_USERM(Q,UM,N) *   CUMSOURCE_UP(UM,NC-1)
              ENDDO
            ENDDO
          ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, KPARS
                L_CUMSOURCE_UP(Q,UM) = L_LAYER_SOURCE(Q,UM) + &
                         TRANS_USERM(UM,N) * L_CUMSOURCE_UP(Q,UM)
              ENDDO
            ENDDO
          ENDIF
        ENDDO

!  NO Offgrid output in this module.........................!!!!!!!!

!  Ongrid output
!  -------------

!  Quadrature output at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

        IF ( DO_INCLUDE_MVOUT .OR. DO_QUAD_RESULTS ) THEN
          CALL QUAD_JACOBIAN_UP_1 &
                 ( NSTREAMS, NLAYERS, K, KPARS, NLEVEL, UTA, & ! Inputs
                   FMULT, KTRANS, LCON, MCON, XPOS, XNEG,    & ! Inputs
                   L_WLOWER, L_WUPPER, L_KTRANS,             & ! Inputs
                   NCON, PCON, L_XPOS, L_XNEG,               & ! Inputs
                   L_QUADELASTIC_F_UP )                        ! Outputs
        ENDIF

!  User-defined stream output, just set to the cumulative source term

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            L_ELASTIC_F_UP(Q,K,UTA,UM) = FMULT * L_CUMSOURCE_UP(Q,UM)
          ENDDO
        ENDDO

!  Check for updating the recursion
!  --------------------------------

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Debug
!        write(*,'(i4,1p4e15.6)')K,
!     *       L_ELASTIC_F_UP(1,K,1,1),L_ELASTIC_F_UP(2,K,1,1),
!     *       L_ELASTIC_F_UP(1,K,2,1),L_ELASTIC_F_UP(2,K,2,1)
!      if (k.eq.18)pause'k18 up'

!  Finish

      RETURN
      END SUBROUTINE ELASTIC_JACOBIAN_UP_1

!

      SUBROUTINE ELASTIC_JACOBIAN_DN_1 &
          ( DO_MSMODE_LIDORT, DO_SSCORRECTION, DO_ELASTIC_SCATTERING, & ! Inputs
            DO_QUAD_RESULTS, DO_INCLUDE_MVOUT,                        & ! Inputs
            FOURIER, NSTREAMS, N_USER_STREAMS, N_LOUTPUT,             & ! Inputs
            FMULT,  LEVELMASK_DN, KVARY, K, KPARS, NKSTORAGE,         & ! Inputs
            KTRANS, TRANS_USERM, LCON, MCON, XPOS, XNEG, U_XPOS,      & ! Inputs
            U_XNEG, U_WNEG1, U_WNEG2, HMULT_1, HMULT_2, EMULT_DN,     & ! Inputs
            CUMSOURCE_DN, L_WLOWER, L_KTRANS, L_TRANS_USERM,          & ! Inputs
            NCON, PCON, L_XPOS, L_XNEG, L_U_XPOS, L_U_XNEG,           & ! Inputs
            L_U_WNEG1, L_U_WNEG2, L_HMULT_1, L_HMULT_2, L_EMULT_DN,   & ! Inputs
            L_ELASTIC_F_DN, L_QUADELASTIC_F_DN )                        ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_LAYERS_NK, MAX_STREAMS, &
                              MAX_2_STREAMS, MAX_USER_STREAMS, MAX_LOUTPUT, MAX_ATMOSWFS

!  Use modules

      USE LRRS_L_POSTPROCESSING_1_m  , Only : LPC_TOASOURCE
      USE LRRS_L_RAMAN_JACOBIANS_1_m , Only : QUAD_JACOBIAN_DN_1

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT
      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION

!  control output for quadrature-angle radiances + fluxes

      LOGICAL  , INTENT(IN) :: DO_QUAD_RESULTS
      LOGICAL  , INTENT(IN) :: DO_INCLUDE_MVOUT

!  elastic scattering flag. If not set for any layers with Fourier > 2,
!  then there will be no layer scattering source terms (source function
!  integration is then trivial)

      LOGICAL  , INTENT(IN) :: DO_ELASTIC_SCATTERING(MAX_LAYERS)

!  Flux factor  and Fourier component

      REAL(FPK), INTENT(IN) :: FMULT
      INTEGER  , INTENT(IN) :: FOURIER

!  Basic control numbers

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      INTEGER  , INTENT(IN) :: N_LOUTPUT

!  Linearization control

      LOGICAL  , INTENT(IN) :: KVARY
      INTEGER  , INTENT(IN) :: K, KPARS, NKSTORAGE(MAX_LAYERS,0:MAX_LAYERS)

!  optical depth output control

      INTEGER  , INTENT(IN) :: LEVELMASK_DN ( MAX_LOUTPUT )

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  Multiplier variables, whole layers

      REAL(FPK), INTENT(IN) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: &
            EMULT_DN   ( MAX_USER_STREAMS, MAX_LAYERS )

!  Cumulative source function

      REAL(FPK), INTENT(IN) :: CUMSOURCE_DN ( MAX_USER_STREAMS, 0:MAX_LAYERS)

!  Linearizations of Direct Beam, PIs, Eigenrans

      REAL(FPK), INTENT(IN) :: L_WLOWER ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  Discrete ordinate solution variables
!  ------------------------------------

!  integration constants

      REAL(FPK), INTENT(IN) :: LCON ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: MCON ( MAX_STREAMS, MAX_LAYERS )

!  homogeneous vectors

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  User homogeneous solutions

      REAL(FPK), INTENT(IN) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: &
            U_WNEG1 ( MAX_USER_STREAMS, MAX_LAYERS ), &
            U_WNEG2 ( MAX_USER_STREAMS, MAX_LAYERS )

!  whole layer transmittance

      REAL(FPK), INTENT(IN) :: TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS )

!  Linearizations
!  --------------

!  integration constants

      REAL(FPK), INTENT(IN) :: &
            NCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS ), &
            PCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  homogeneous vectors

      REAL(FPK), INTENT(IN) :: L_XPOS &
        ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XNEG &
        ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  homogeneous solutions at user angles

      REAL(FPK), INTENT(IN) :: L_U_XPOS &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_XNEG &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  classical particular integral solutions at user angles

      REAL(FPK), INTENT(IN) :: L_U_WNEG1 &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_WNEG2 &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS_NK )

!  whole layer transmittances

      REAL(FPK), INTENT(IN) :: L_TRANS_USERM &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS )

!  Whole layer multipliers

      REAL(FPK), INTENT(IN) :: L_HMULT_1 &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_HMULT_2 &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_EMULT_DN &
            ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS, MAX_LAYERS )

!  output = elastic scattering Fourier-component Jacobians
!  -------------------------------------------------
!mick fix 10/19/2015 - changed intent to inout
      REAL(FPK), INTENT(INOUT) :: L_ELASTIC_F_DN &
            ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS )
      REAL(FPK), INTENT(INOUT) :: L_QUADELASTIC_F_DN &
            ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

      INTEGER   :: N, NK, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER   :: UTA, UM, Q
      LOGICAL   :: SOURCETERM_FLAG

!  TOA source (usually zero, no diffuse downwelling field at TOA)

      REAL(FPK) :: L_TOA_SOURCE ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  Local layer source termes (total and diffuse-only)

      REAL(FPK) :: &
          L_LAYER_SOURCE      ( MAX_ATMOSWFS, MAX_USER_STREAMS )
      REAL(FPK) :: &
          L_MSCAT_LAYERSOURCE ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  Cumulative Jacobian source term

      REAL(FPK) :: L_CUMSOURCE_DN( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  Exit if no variation

      IF ( .NOT. KVARY ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            DO UTA = 1, N_LOUTPUT
              L_ELASTIC_F_DN(Q,K,UTA,UM) = zero
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Initialize recursion
!  --------------------

!  Get the source term

      CALL LPC_TOASOURCE &
          ( N_USER_STREAMS, KPARS, L_TOA_SOURCE )

!  Start the recurrence with this term

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, KPARS
          L_CUMSOURCE_DN(Q,UM) = L_TOA_SOURCE(Q,UM)
        ENDDO
      ENDDO

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

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        NUT = NLEVEL
        DO N = NSTART, NUT
          NC = N
          NK = NKSTORAGE(N,K)
          SOURCETERM_FLAG = DO_ELASTIC_SCATTERING(N)

          CALL L_ELASTIC_LSST_DN &
                ( DO_MSMODE_LIDORT, DO_SSCORRECTION, SOURCETERM_FLAG,     & ! Inputs
                  N, K, KPARS, NK, NSTREAMS, N_USER_STREAMS, FOURIER,     & ! Inputs
                  HMULT_1, HMULT_2, EMULT_DN, LCON, MCON, U_XPOS, U_XNEG, & ! Inputs
                  U_WNEG1, U_WNEG2, L_HMULT_1, L_HMULT_2, L_EMULT_DN,     & ! Inputs
                  NCON, PCON, L_U_XPOS, L_U_XNEG, L_U_WNEG1, L_U_WNEG2,   & ! Inputs
                  L_LAYER_SOURCE, L_MSCAT_LAYERSOURCE )                     ! Outputs

          IF ( N.EQ.K .OR. K.EQ.0 ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, KPARS
                L_CUMSOURCE_DN(Q,UM) = L_LAYER_SOURCE(Q,UM) + &
                         TRANS_USERM(UM,N)   * L_CUMSOURCE_DN(Q,UM) + &
                       L_TRANS_USERM(Q,UM,N) *   CUMSOURCE_DN(UM,NC-1)
              ENDDO
            ENDDO
          ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, KPARS
                L_CUMSOURCE_DN(Q,UM) = L_LAYER_SOURCE(Q,UM) + &
                         TRANS_USERM(UM,N) * L_CUMSOURCE_DN(Q,UM)
              ENDDO
            ENDDO
          ENDIF
        ENDDO

!  NO Offgrid output in this module.................!!!!!

!  Ongrid output
!  -------------

!  Quadrature results at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

        IF ( DO_INCLUDE_MVOUT .OR. DO_QUAD_RESULTS  ) THEN
          CALL QUAD_JACOBIAN_DN_1 &
                 ( NSTREAMS, K, KPARS, NLEVEL, UTA, FMULT,   & ! Inputs
                   KTRANS, LCON, MCON, XPOS, XNEG, L_WLOWER, & ! Inputs
                   L_KTRANS, NCON, PCON, L_XPOS, L_XNEG,     & ! Inputs
                   L_QUADELASTIC_F_DN )                        ! Outputs
        ENDIF

!  User-defined stream output, just set to the cumulative source term

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            L_ELASTIC_F_DN(Q,K,UTA,UM) = FMULT * L_CUMSOURCE_DN(Q,UM)
          ENDDO
        ENDDO

!  Check for updating the recursion
!  --------------------------------

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  debug
!        write(*,'(i4,1p4e15.6)')K,
!     *       L_ELASTIC_F_DN(1,K,1,1),L_ELASTIC_F_DN(2,K,1,1),
!     *       L_ELASTIC_F_DN(1,K,2,1),L_ELASTIC_F_DN(2,K,2,1)
!      if (k.eq.18)pause'k18 dn'

!  Finish

      RETURN
      END SUBROUTINE ELASTIC_JACOBIAN_DN_1

!

      SUBROUTINE L_ELASTIC_LSST_UP &
            ( DO_MSMODE_LIDORT, DO_SSCORRECTION, SOURCETERM_FLAG,     & ! Inputs
              N, K, KPARS, NK, NSTREAMS, N_USER_STREAMS, FOURIER,     & ! Inputs
              HMULT_1, HMULT_2, EMULT_UP, LCON, MCON, U_XPOS, U_XNEG, & ! Inputs
              U_WPOS1, U_WPOS2, L_HMULT_1, L_HMULT_2, L_EMULT_UP,     & ! Inputs
              NCON, PCON, L_U_XPOS, L_U_XNEG, L_U_WPOS1, L_U_WPOS2,   & ! Inputs
              L_LAYER_SOURCE, L_MSCAT_LAYERSOURCE )                     ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_LAYERS_NK, MAX_2_STREAMS, MAX_STREAMS, &
                              MAX_USER_STREAMS, MAX_LOUTPUT, MAX_ATMOSWFS

      IMPLICIT NONE

!  arguments
!  ---------

!  input flags and control

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT, DO_SSCORRECTION
      LOGICAL  , INTENT(IN) :: SOURCETERM_FLAG

!  Streams and Fourier

      INTEGER  , INTENT(IN) :: NSTREAMS, N_USER_STREAMS, FOURIER

!  Layer No., linearization control integers

      INTEGER  , INTENT(IN) :: N, K, KPARS, NK

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

!  Linearization inputs
!  --------------------

!  integration constants

      REAL(FPK), INTENT(IN) :: &
            NCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS ), &
            PCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  homogeneous solutions at user angles

      REAL(FPK), INTENT(IN) :: L_U_XPOS &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_XNEG &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  classical particular integral solutions at user angles

      REAL(FPK), INTENT(IN) :: L_U_WPOS1 &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_WPOS2 &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS_NK )

!  Whole layer multipliers

      REAL(FPK), INTENT(IN) :: L_HMULT_1 &
          ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_HMULT_2 &
          ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_EMULT_UP &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS, MAX_LAYERS )

!  output layer source terms
!  -------------------------

      REAL(FPK), INTENT(OUT) :: L_LAYER_SOURCE      ( MAX_ATMOSWFS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: L_MSCAT_LAYERSOURCE ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER   :: AA, UM, Q, M
      REAL(FPK) :: SHOM, S2, S1, LCON_UXVEC, MCON_UXVEC

!  No layer source term if no scattering in the layer

      IF ( .NOT. SOURCETERM_FLAG ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            L_LAYER_SOURCE(Q,UM) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Fourier number (debug only)

      M = FOURIER

!  Homogeneous solutions
!  =====================

!  Special case when N = K, or K = 0 (bulk)

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              LCON_UXVEC = NCON(Q,AA,N) *   U_XPOS(UM,AA,N) &
                         + LCON(AA,N)   * L_U_XPOS(Q,UM,AA,N)
              MCON_UXVEC = PCON(Q,AA,N) *   U_XNEG(UM,AA,N) &
                         + MCON(AA,N)   * L_U_XNEG(Q,UM,AA,N)
              SHOM = SHOM + &
                 LCON_UXVEC * HMULT_2(AA,UM,N) &
               + LCON(AA,N) * U_XPOS(UM,AA,N) * L_HMULT_2(Q,AA,UM,N) &
               + MCON_UXVEC * HMULT_1(AA,UM,N) &
               + MCON(AA,N) * U_XNEG(UM,AA,N) * L_HMULT_1(Q,AA,UM,N)
            ENDDO
            L_MSCAT_LAYERSOURCE(Q,UM) = SHOM
          ENDDO
        ENDDO

!  Other cases when N not equal to K (only variation of Integ-Cons)

      ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              LCON_UXVEC = NCON(Q,AA,N) *   U_XPOS(UM,AA,N)
              MCON_UXVEC = PCON(Q,AA,N) *   U_XNEG(UM,AA,N)
              SHOM = SHOM + LCON_UXVEC * HMULT_2(AA,UM,N) &
                          + MCON_UXVEC * HMULT_1(AA,UM,N)
            ENDDO
            L_MSCAT_LAYERSOURCE(Q,UM) = SHOM
          ENDDO
        ENDDO
      ENDIF

!  Particular and single scatter contributions - classical solution
!  ================================================================

!  add particular solution

      IF ( N.GE.K .OR. K.EQ.0 ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            S2 =  L_EMULT_UP(Q,K,UM,N) *   U_WPOS2(UM,N) &
                  + EMULT_UP(UM,N)     * L_U_WPOS2(Q,UM,NK)
            L_MSCAT_LAYERSOURCE(Q,UM) = L_MSCAT_LAYERSOURCE(Q,UM) + S2
          ENDDO
        ENDDO
      ENDIF

!  Set the term so far

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, KPARS
          L_LAYER_SOURCE(Q,UM) = L_MSCAT_LAYERSOURCE(Q,UM)
        ENDDO
      ENDDO

!  No single scattering contribution if flagged

      IF ( DO_MSMODE_LIDORT .OR. DO_SSCORRECTION ) THEN

!  Single scattering contribution when appropriate
!       Special case when N = K, or K = 0 (bulk)
!        Other cases when N > K

      ELSE
        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, KPARS
              S1 =    U_WPOS1(UM,N)   * L_EMULT_UP(Q,K,UM,N) &
                  + L_U_WPOS1(Q,UM,N) *   EMULT_UP(UM,N)
              L_LAYER_SOURCE(Q,UM) = L_MSCAT_LAYERSOURCE(Q,UM) + S1
            ENDDO
          ENDDO
        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, KPARS
              S1 =    U_WPOS1(UM,N)   * L_EMULT_UP(Q,K,UM,N)
              L_LAYER_SOURCE(Q,UM) = L_MSCAT_LAYERSOURCE(Q,UM) + S1
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE L_ELASTIC_LSST_UP

!

      SUBROUTINE L_ELASTIC_LSST_DN &
            ( DO_MSMODE_LIDORT, DO_SSCORRECTION, SOURCETERM_FLAG,     & ! Inputs
              N, K, KPARS, NK, NSTREAMS, N_USER_STREAMS, FOURIER,     & ! Inputs
              HMULT_1, HMULT_2, EMULT_DN, LCON, MCON, U_XPOS, U_XNEG, & ! Inputs
              U_WNEG1, U_WNEG2, L_HMULT_1, L_HMULT_2, L_EMULT_DN,     & ! Inputs
              NCON, PCON, L_U_XPOS, L_U_XNEG, L_U_WNEG1, L_U_WNEG2,   & ! Inputs
              L_LAYER_SOURCE, L_MSCAT_LAYERSOURCE )                     ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_LAYERS_NK, MAX_STREAMS, MAX_2_STREAMS, &
                              MAX_USER_STREAMS, MAX_LOUTPUT, MAX_ATMOSWFS

      IMPLICIT NONE

!  arguments
!  ---------

!  input flags and control

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LIDORT, DO_SSCORRECTION
      LOGICAL  , INTENT(IN) :: SOURCETERM_FLAG

!  Streams and Fourier

      INTEGER  , INTENT(IN) :: NSTREAMS, N_USER_STREAMS, FOURIER

!  Layer No., linearization control integers

      INTEGER  , INTENT(IN) :: N, K, KPARS, NK

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

!  Linearization inputs
!  --------------------

!  integration constants

      REAL(FPK), INTENT(IN) :: &
            NCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS ), &
            PCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  homogeneous solutions at user angles

      REAL(FPK), INTENT(IN) :: L_U_XPOS &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_XNEG &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  classical particular integral solutions at user angles

      REAL(FPK), INTENT(IN) :: L_U_WNEG1 &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_WNEG2 &
          ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS_NK )

!  Whole layer multipliers

      REAL(FPK), INTENT(IN) :: L_HMULT_1 &
          ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_HMULT_2 &
          ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_EMULT_DN &
          ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS, MAX_LAYERS )

!  output layer source terms
!  -------------------------

      REAL(FPK), INTENT(OUT) :: L_LAYER_SOURCE      ( MAX_ATMOSWFS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: L_MSCAT_LAYERSOURCE ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER   :: AA, UM, Q, M
      REAL(FPK) :: SHOM, S2, S1, LCON_UXVEC, MCON_UXVEC

!  No layer source term if no scattering in the layer

      IF ( .NOT. SOURCETERM_FLAG ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            L_LAYER_SOURCE(Q,UM) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Fourier number (debug only)

      M = FOURIER

!  Homogeneous solutions
!  =====================

!  Special case when N = K, or K = 0 (bulk)

      IF ( N.EQ.K .OR. K.EQ.0 ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              LCON_UXVEC = NCON(Q,AA,N) *   U_XNEG(UM,AA,N) &
                         + LCON(AA,N)   * L_U_XNEG(Q,UM,AA,N)
              MCON_UXVEC = PCON(Q,AA,N) *   U_XPOS(UM,AA,N) &
                         + MCON(AA,N)   * L_U_XPOS(Q,UM,AA,N)
              SHOM = SHOM + &
                 LCON_UXVEC * HMULT_1(AA,UM,N) &
               + LCON(AA,N) * U_XNEG(UM,AA,N) * L_HMULT_1(Q,AA,UM,N) &
               + MCON_UXVEC * HMULT_2(AA,UM,N) &
               + MCON(AA,N) * U_XPOS(UM,AA,N) * L_HMULT_2(Q,AA,UM,N)
            ENDDO
            L_MSCAT_LAYERSOURCE(Q,UM) = SHOM
          ENDDO
        ENDDO

!  Other cases when N not equal to K (only variation of Integ-Cons)

      ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              LCON_UXVEC = NCON(Q,AA,N) *   U_XNEG(UM,AA,N)
              MCON_UXVEC = PCON(Q,AA,N) *   U_XPOS(UM,AA,N)
              SHOM = SHOM + LCON_UXVEC * HMULT_1(AA,UM,N) &
                          + MCON_UXVEC * HMULT_2(AA,UM,N)
            ENDDO
            L_MSCAT_LAYERSOURCE(Q,UM) = SHOM
          ENDDO
        ENDDO
      ENDIF

!  Particular and single scatter contributions - classical solution
!  ================================================================

!  add particular solution

      IF ( N.GE.K .OR. K.EQ.0 ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, KPARS
            S2 =  L_EMULT_DN(Q,K,UM,N) *   U_WNEG2(UM,N) &
                  + EMULT_DN(UM,N)     * L_U_WNEG2(Q,UM,NK)
            L_MSCAT_LAYERSOURCE(Q,UM) = L_MSCAT_LAYERSOURCE(Q,UM) + S2
          ENDDO
        ENDDO
      ENDIF

!  Set contribution so far

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, KPARS
          L_LAYER_SOURCE(Q,UM) = L_MSCAT_LAYERSOURCE(Q,UM)
        ENDDO
      ENDDO

!  No single scattering contribution if flagged

      IF ( DO_MSMODE_LIDORT .OR. DO_SSCORRECTION ) THEN

!  Single scattering contribution when appropriate
!       Special case when N = K, or K = 0 (bulk)
!        Other cases when N > K

      ELSE
        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, KPARS
              S1 =    U_WNEG1(UM,N)   * L_EMULT_DN(Q,K,UM,N) &
                  + L_U_WNEG1(Q,UM,N) *   EMULT_DN(UM,N)
              L_LAYER_SOURCE(Q,UM) = L_MSCAT_LAYERSOURCE(Q,UM) + S1
            ENDDO
          ENDDO
        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, KPARS
              S1 =    U_WNEG1(UM,N)   * L_EMULT_DN(Q,K,UM,N)
              L_LAYER_SOURCE(Q,UM) = L_MSCAT_LAYERSOURCE(Q,UM) + S1
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE L_ELASTIC_LSST_DN

!

      SUBROUTINE ELASTIC_SURFJAC_UP_1 &
          ( DO_SSCORRECTION, DO_QUAD_RESULTS, DO_INCLUDE_MVOUT,        & ! Inputs
            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR,       & ! Inputs
            ALBEDOS_RANKED, Raman_Idx,                                 & ! Inputs
            FOURIER, NSTREAMS, NLAYERS, N_USER_STREAMS, N_LOUTPUT,     & ! Inputs
            N_SURFACEWFS, FLUXMULT, LEVELMASK_UP, Brdf_Idx,            & ! Inputs
            USER_BRDF_F, LS_USER_BRDF_F, LS_USER_DIRECT_BEAM,          & ! Inputs
            TRANS_USERM, QUAD_STRMWGHT, XPOS, XNEG, KTRANS, IDOWNSURF, & ! Inputs
            U_XPOS, U_XNEG, HMULT_1, HMULT_2, NCON_SURF, PCON_SURF,    & ! Inputs
            LS_ELASTIC_F_UP, LS_QUADELASTIC_F_UP )                       ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_2_STREAMS, MAX_STREAMS, MAX_USER_STREAMS, &
                              MAX_MOMENTS, MAX_LOUTPUT, MAX_POINTS, MAX_SURFACEWFS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , INTENT(IN) :: DO_SSCORRECTION
      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE
      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE

!  output for quadrature-angle radiances disabled in this version

      LOGICAL  , INTENT(IN) :: DO_QUAD_RESULTS
      LOGICAL  , INTENT(IN) :: DO_INCLUDE_MVOUT

!  Flux multiplier and Fourier component

      REAL(FPK), INTENT(IN) :: FLUXMULT
      INTEGER  , INTENT(IN) :: FOURIER

!  Basic control numbers

      INTEGER  , INTENT(IN) :: NSTREAMS, NLAYERS
      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      INTEGER  , INTENT(IN) :: N_LOUTPUT

!  Number of Surface Jacobians

      INTEGER  , INTENT(IN) :: N_SURFACEWFS

!  Indices, Version 2.5, 9/11/15

      !INTEGER  , INTENT(IN) :: Brdf_IDX ! Raman_IDX
      INTEGER  , INTENT(IN) :: Brdf_IDX, Raman_IDX

!  2-deltam0 factor

      REAL(fpk), intent(in) :: SURFACE_FACTOR

!  Albedos. Not needed
!      REAL(fpk), intent(in) :: ALBEDOS_RANKED ( MAX_POINTS)
      REAL(fpk), intent(in) :: ALBEDOS_RANKED ( MAX_POINTS)

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

!mick fix 7/20/2016 - swapped stream dimensions
      !REAL(fpk), intent(in) :: USER_BRDF_F    ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_USER_STREAMS, MAX_POINTS )
      !REAL(fpk), intent(in) :: LS_USER_BRDF_F ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_USER_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: USER_BRDF_F    ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LS_USER_BRDF_F ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Direct beam

      REAL(FPK), INTENT(IN) :: LS_USER_DIRECT_BEAM ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  optical depth output control

      INTEGER  , INTENT(IN) :: LEVELMASK_UP ( MAX_LOUTPUT )

!  Quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT ( MAX_STREAMS )

!  diffuse integrated field at surface

      REAL(FPK), INTENT(IN) :: IDOWNSURF ( MAX_STREAMS )

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  Discrete ordinate solution variables
!  ------------------------------------

!  homogeneous vectors

      REAL(FPK), INTENT(IN) :: &
          XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  homogeneous solutions at user angles

      REAL(FPK), INTENT(IN) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  whole and partial layer transmittances

      REAL(FPK), INTENT(IN) :: TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS )

!  Whole layer multipliers

      REAL(FPK), INTENT(IN) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Linearized integration constants

      REAL(FPK), INTENT(IN) :: &
            NCON_SURF ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS ), &
            PCON_SURF ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS )

!  output = elastic scattering Fourier-component Jacobians
!  -------------------------------------------------

      REAL(FPK), INTENT(OUT) :: LS_ELASTIC_F_UP &
            ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: LS_QUADELASTIC_F_UP &
            ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_STREAMS )

!  Local variables
!  ---------------

!  Local help

      REAL(FPK) :: H1, H2, SUMR, SHOM, NCON_HELP, PCON_HELP
      REAL(FPK) :: REFLEC, LS_IDOWNSURF ( MAX_STREAMS )
      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER   :: UTA, UM, I, I1, NL, AA, M, Q

!  Local BOA source terms

      REAL(FPK) :: LS_BOA_SOURCE     ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Local layer source termes (total and diffuse-only)

      REAL(FPK) :: LS_LAYER_SOURCE   ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Cumulative Jacobian source term

      REAL(FPK) :: LS_CUMSOURCE_UP   ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Debug

      REAL(FPK) :: LS_BOA_SOURCE_2ND ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Initialize recursion.
!  ---------------------

!  Initialise

      DO UM = 1, N_USER_STREAMS
        LS_BOA_SOURCE(1:N_SURFACEWFS, UM) = ZERO
      ENDDO

!  Bottom layer

      N = NLAYERS

!  Fourier number

      M = FOURIER

!  Get the radiance field at the surface

      IF ( DO_INCLUDE_SURFACE ) THEN

!  start loop over surface Jacobians

        DO Q = 1, N_SURFACEWFS

!  First compute derivative of downward intensity Integrand at stream angle
!        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

          DO I = 1, NSTREAMS
            SUMR = ZERO
            DO AA = 1, NSTREAMS
              H1  = NCON_SURF(Q,AA,N) * XPOS(I,AA,N) * KTRANS(AA,N)
              H2  = PCON_SURF(Q,AA,N) * XNEG(I,AA,N)
              SUMR = SUMR + H1 + H2
            ENDDO
            LS_IDOWNSURF(I) = SUMR * QUAD_STRMWGHT(I)
          ENDDO

!  Sum the integrand (Lambertian case), same for all user-streams

          IF ( .not. DO_BRDF_SURFACE ) THEN
            IF ( M == 0 ) then
!mick fix 10/19/2015 - added ALBEDOS_RANKED(Raman_Idx) to term
              !REFLEC = SUM(LS_IDOWNSURF(1:NSTREAMS)) * SURFACE_FACTOR
              REFLEC = SUM(LS_IDOWNSURF(1:NSTREAMS)) * ALBEDOS_RANKED(Raman_Idx) * SURFACE_FACTOR
              LS_BOA_SOURCE(Q,1:N_USER_STREAMS) = REFLEC
            ENDIF
          ELSE
            DO UM = 1, N_USER_STREAMS
              REFLEC = DOT_PRODUCT(LS_IDOWNSURF(1:NSTREAMS), USER_BRDF_F(M,UM,1:NSTREAMS,Brdf_Idx)) &
                     * SURFACE_FACTOR
              LS_BOA_SOURCE(Q,UM) = REFLEC
            ENDDO
          ENDIF

!  Add linearization due to Surface variation, diffuse term

          IF ( .not. DO_BRDF_SURFACE ) THEN
            IF ( M == 0 ) then
              REFLEC = SUM(IDOWNSURF(1:NSTREAMS)) * SURFACE_FACTOR
              LS_BOA_SOURCE(Q,1:N_USER_STREAMS) = LS_BOA_SOURCE(Q,1:N_USER_STREAMS) + REFLEC
            ENDIF
          ELSE
            DO UM = 1, N_USER_STREAMS
              REFLEC = DOT_PRODUCT(IDOWNSURF(1:NSTREAMS), LS_USER_BRDF_F(Q,M,UM,1:NSTREAMS,Brdf_Idx)) &
                     * SURFACE_FACTOR
              LS_BOA_SOURCE(Q,UM) = LS_BOA_SOURCE(Q,UM) + REFLEC
            ENDDO
          ENDIF

!  Add linearization due to Albedo variation of direct beam

          IF ( .NOT. DO_SSCORRECTION ) THEN
            DO UM = 1, N_USER_STREAMS
              LS_BOA_SOURCE(Q,UM) = LS_BOA_SOURCE(Q,UM) + LS_USER_DIRECT_BEAM(Q,UM)
            ENDDO
          ENDIF

!  End WF loop

        ENDDO

!  End clause

      ENDIF

!  Initialise cumulative term

      DO UM = 1, N_USER_STREAMS
        LS_CUMSOURCE_UP(1:N_SURFACEWFS,UM) = LS_BOA_SOURCE(1:N_SURFACEWFS,UM)
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

!  Layer index for given optical depth

        NLEVEL = LEVELMASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        NUT = NLEVEL + 1
        DO N = NSTART, NUT, -1
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SURFACEWFS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                NCON_HELP =  NCON_SURF(Q,AA,N) * U_XPOS(UM,AA,N)
                PCON_HELP =  PCON_SURF(Q,AA,N) * U_XNEG(UM,AA,N)
                H1 = NCON_HELP * HMULT_2(AA,UM,N)
                H2 = PCON_HELP * HMULT_1(AA,UM,N)
                SHOM = SHOM + H1 + H2
              ENDDO
              LS_LAYER_SOURCE(Q,UM) = SHOM
              LS_CUMSOURCE_UP(Q,UM) = LS_LAYER_SOURCE(Q,UM) + &
                                      TRANS_USERM(UM,N) * LS_CUMSOURCE_UP(Q,UM)
            ENDDO
          ENDDO
        ENDDO

!  NO Offgrid output in this module.........!!!!!!!!!

!  Ongrid output
!  -------------

!  Quadrature output at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

        IF ( DO_INCLUDE_MVOUT .OR. DO_QUAD_RESULTS ) THEN

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we're
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
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 = NCON_SURF(Q,AA,NL)*XPOS(I1,AA,NL)*KTRANS(AA,NL)
                  H2 = PCON_SURF(Q,AA,NL)*XNEG(I1,AA,NL)
                  SHOM = SHOM + H1 + H2
                ENDDO
                LS_QUADELASTIC_F_UP(Q,UTA,I) = FLUXMULT * SHOM
              ENDDO
            ENDDO

!  For other levels in the atmosphere

          ELSE

            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              DO Q = 1, N_SURFACEWFS
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 = NCON_SURF(Q,AA,N)*XPOS(I1,AA,N)
                  H2 = PCON_SURF(Q,AA,N)*XNEG(I1,AA,N)*KTRANS(AA,N)
                  SHOM = SHOM + H1 + H2
                ENDDO
                LS_QUADELASTIC_F_UP(Q,UTA,I) = FLUXMULT * SHOM
              ENDDO
            ENDDO

          ENDIF

        ENDIF

!  User-defined stream output, just set to the cumulative source term
!mick eff 7/20/2016

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SURFACEWFS
            LS_ELASTIC_F_UP(Q,UTA,UM) = FLUXMULT * LS_CUMSOURCE_UP(Q,UM)
          ENDDO
        ENDDO

!  Check for updating the recursion
!  --------------------------------

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE ELASTIC_SURFJAC_UP_1

!

      SUBROUTINE ELASTIC_SURFJAC_DN_1 &
          ( DO_QUAD_RESULTS, DO_INCLUDE_MVOUT,              & ! Inputs
            FOURIER, NSTREAMS,  N_USER_STREAMS, N_LOUTPUT,  & ! Inputs
            N_SURFACEWFS, FLUXMULT, LEVELMASK_DN,           & ! Inputs
            TRANS_USERM, XPOS, XNEG, KTRANS, U_XPOS,        & ! Inputs
            U_XNEG, HMULT_1, HMULT_2, NCON_SURF, PCON_SURF, & ! Inputs
            LS_ELASTIC_F_DN, LS_QUADELASTIC_F_DN )            ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_2_STREAMS, MAX_STREAMS, MAX_USER_STREAMS, &
                              MAX_LOUTPUT,  MAX_SURFACEWFS

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  output for quadrature-angle radiances

      LOGICAL  , INTENT(IN) :: DO_QUAD_RESULTS

!  local control flags

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_MVOUT

!  Flux multiplier and Fourier component

      REAL(FPK), INTENT(IN) :: FLUXMULT
      INTEGER  , INTENT(IN) :: FOURIER

!  Basic control numbers

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      INTEGER  , INTENT(IN) :: N_LOUTPUT

!  Number of Surface Jacobians

      INTEGER  , INTENT(IN) :: N_SURFACEWFS

!  optical depth output control

      INTEGER  , INTENT(IN) :: LEVELMASK_DN ( MAX_LOUTPUT )

!  Discrete ordinate solution variables
!  ------------------------------------

!  Whole layer Eigensolution transmittances

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!  homogeneous vectors

      REAL(FPK), INTENT(IN) :: &
          XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
          XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  homogeneous solutions at user angles

      REAL(FPK), INTENT(IN) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  whole layer transmittances

      REAL(FPK), INTENT(IN) :: TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS )

!  Whole layer multipliers

      REAL(FPK), INTENT(IN) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Linearized integration constants

      REAL(FPK), INTENT(IN) :: &
            NCON_SURF ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS ), &
            PCON_SURF ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS )

!  output = elastic scattering Fourier-component Jacobians
!  -------------------------------------------------

      REAL(FPK), INTENT(OUT) :: LS_ELASTIC_F_DN      ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: LS_QUADELASTIC_F_DN  ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_STREAMS )

!  local variables
!  ---------------

!  local help

      REAL(FPK) :: H1, H2, SHOM, NCON_HELP, PCON_HELP
      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER   :: UTA, UM, I, NL, AA, M, Q

!  local BOA source terms

      REAL(FPK) :: LS_TOA_SOURCE    ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Local layer source termes (total and diffuse-only)

      REAL(FPK) :: LS_LAYER_SOURCE  ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Cumulative Jacobian source term

      REAL(FPK) :: LS_CUMSOURCE_DN ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Initialize recursion.
!  ---------------------

      DO UM = 1, N_USER_STREAMS
        LS_TOA_SOURCE(1:N_SURFACEWFS,UM)   = ZERO
        LS_CUMSOURCE_DN(1:N_SURFACEWFS,UM) = LS_TOA_SOURCE(1:N_SURFACEWFS,UM)
      ENDDO

!  Fourier number (debug only)

      M = FOURIER

!  initialise cumulative source term loop

      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

!  Start loop with the lowest value first

      DO UTA = 1,  N_LOUTPUT

!  Layer index for given optical depth

        NLEVEL = LEVELMASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term

        NUT = NLEVEL
        DO N = NSTART, NUT
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SURFACEWFS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                NCON_HELP =  NCON_SURF(Q,AA,N) * U_XNEG(UM,AA,N)
                PCON_HELP =  PCON_SURF(Q,AA,N) * U_XPOS(UM,AA,N)
                H1 = NCON_HELP * HMULT_1(AA,UM,N)
                H2 = PCON_HELP * HMULT_2(AA,UM,N)
                SHOM = SHOM + H1 + H2
              ENDDO
              LS_LAYER_SOURCE(Q,UM) = SHOM
              LS_CUMSOURCE_DN(Q,UM) = LS_LAYER_SOURCE(Q,UM) + &
                      TRANS_USERM(UM,N)   * LS_CUMSOURCE_DN(Q,UM)
            ENDDO
          ENDDO
        ENDDO

!  NO Offgrid output in this module...........!!!!!!!!!

!  Ongrid output
!  -------------

!  Quadrature output at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

        IF ( DO_INCLUDE_MVOUT .OR. DO_QUAD_RESULTS ) THEN

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we're
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

          NL = NLEVEL
          N = NL

!  For the highest level, zero. 

          IF ( NLEVEL .EQ. 0 ) THEN
            DO I = 1, NSTREAMS
              DO Q = 1, N_SURFACEWFS
                LS_QUADELASTIC_F_DN(Q,UTA,I) = ZERO
              ENDDO
            ENDDO
          ELSE
            DO I = 1, NSTREAMS
              DO Q = 1, N_SURFACEWFS
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 = NCON_SURF(Q,AA,N)*XPOS(I,AA,N)*KTRANS(AA,N)
                  H2 = PCON_SURF(Q,AA,N)*XNEG(I,AA,N)
                  SHOM = SHOM + H1 + H2
                ENDDO
                LS_QUADELASTIC_F_DN(Q,UTA,I) = FLUXMULT * SHOM
              ENDDO
            ENDDO
          ENDIF

        ENDIF

!  User-defined stream output, just set to the cumulative source term

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SURFACEWFS
            LS_ELASTIC_F_DN(Q,UTA,UM) = FLUXMULT * LS_CUMSOURCE_DN(Q,UM)
          ENDDO
        ENDDO

!  Check for updating the recursion
!  --------------------------------

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE ELASTIC_SURFJAC_DN_1

!  End Module

      END MODULE lrrs_L_elastic_jacobians_1_m

