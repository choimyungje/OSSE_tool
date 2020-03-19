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
! #   Setup SUBROUTINES (formerly in lrrs_rtsolutions_1.f90)    #
! #           (1 = Whole-layer, 2 = Part-layer)                 #
! #                                                             #
! #             CHAPMAN_FUNCTION                                #
! #             ATTENUATION_SETUP_1                             #
! #             LEGENDRE_SETUP                                  #
! #             DBEAM_SETUP                                     #
! #             UDBEAM_SETUP                                    #
! #                                                             #
! #    New Setup routine for the DO_SSCORR_NADIR option         #
! #                                                             #
! #             BEAM_MULTIPLIERS_1                              #
! #                                                             #
! ###############################################################

      MODULE lrrs_miscsetups_1_m

!   -- Rob mod 5/12/17 for 2p5a. This is a New Module with 6 subroutines
!        5 old subroutines formerly in lrrs_rtsolutions_1.f90
!        1 new subroutine  (BEAM_MULTIPLIERS_1) for the DO_SSCORR_NADIR option

!      USE LRRS_PARS_m, Only : SDU

!   Everything public here.

      PUBLIC

      CONTAINS

!

      SUBROUTINE BEAM_MULTIPLIERS_1 &
        ( DO_UPWELLING, DO_DNWELLING, NPOINTS, NLAYERS, & ! Inputs
          N_USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER,   & ! Inputs
          STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,       & ! Inputs
          USER_STREAMS, DELTAU_VERT, SAVE_TRANS_USERM,  & ! Inputs
          BEAM_PICUTOFF, BEAM_AVSECANT, BEAM_DTRANS, BEAM_ITRANS,  & ! Inputs
          SAVE_BEAMMULT_UP, SAVE_BEAMMULT_DN )            ! Output

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, MAX_TAU_SPATH, MAX_TAU_UPATH, &
                              MAX_LAYERS, MAX_USER_STREAMS, MAX_POINTS

!  Dependency

      USE lrrs_postprocessing_1_m, Only : BEAMMULT_UP_1, BEAMMULT_DN_1

!  Implicit none

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Flags

      LOGICAL  , INTENT(IN) :: DO_UPWELLING, DO_DNWELLING

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

!  Output multipliers
!  ------------------

      REAL(FPK), INTENT(OUT) :: SAVE_BEAMMULT_UP ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: SAVE_BEAMMULT_DN ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Local variables
!  ---------------

      INTEGER   :: N, S
      REAL(FPK) :: DELTAUS, DELTRANS, ASOURCE
      REAL(FPK) :: TRANS_USERM (MAX_USER_STREAMS)
      REAL(FPK) :: EMULT_UP    (MAX_USER_STREAMS)
      REAL(FPK) :: EMULT_DN    (MAX_USER_STREAMS)
      LOGICAL, parameter :: FLIPPER = .false.

!  Initialize - No need
!      SAVE_BEAMMULT_UP = zero
!      SAVE_BEAMMULT_DN = zero

!  Monochromatic - 234 points (233 shifts + 1 calculation point)
!  Binned        - all points in outer buffer

!  Upwelling

      IF ( DO_UPWELLING ) THEN
        DO S = 1, NPOINTS
          DO N = 1, NLAYERS
            DELTAUS  = DELTAU_VERT(N,S)
            DELTRANS = BEAM_DTRANS(N,S)
            ASOURCE  = BEAM_AVSECANT(N,S)
            TRANS_USERM(:) = SAVE_TRANS_USERM(:,N,S) 
            CALL BEAMMULT_UP_1 &
              ( N, N_USER_STREAMS, USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER, & ! Inputs
                STERM_LAYERMASK_UP(N), TRANS_USERM, DELTAUS,   & ! Inputs
                FLIPPER, BEAM_PICUTOFF(S), DELTRANS, ASOURCE,  & ! Inputs
                EMULT_UP )                                  ! Output
            SAVE_BEAMMULT_UP(:,N,S) = EMULT_UP(:) * BEAM_ITRANS(N,S)
          ENDDO
        ENDDO
      ENDIF

!  Downwelling

      IF ( DO_DNWELLING ) THEN
        DO S = 1, NPOINTS
          DO N = 1, NLAYERS
            DELTAUS  = DELTAU_VERT(N,S)
            DELTRANS = BEAM_DTRANS(N,S)
            ASOURCE  = BEAM_AVSECANT(N,S)
            TRANS_USERM(:) = SAVE_TRANS_USERM(:,N,S) 
            CALL BEAMMULT_DN_1 &
              ( N, N_USER_STREAMS, USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER, & ! Inputs
                STERM_LAYERMASK_DN(N), TRANS_USERM, DELTAUS,   & ! Inputs
                FLIPPER, BEAM_PICUTOFF(S), DELTRANS, ASOURCE,  & ! Inputs
                EMULT_DN )                                  ! Output
            SAVE_BEAMMULT_DN(:,N,S) = EMULT_DN(:) * BEAM_ITRANS(N,S)
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BEAM_MULTIPLIERS_1

!

      SUBROUTINE ATTENUATION_SETUP_1 &
          ( NPOINTS, NLAYERS, DO_PLANE_PARALLEL,           & ! Input
            DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS, & ! Input
            DELTAU_VERT_INPUT, CHAPMAN, COS_SZA,           & ! Input
            BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT,     & ! Input
            BEAM_ETRANS, BEAM_DTRANS, SAVE_TRANS_USERM )     ! Output

!  Module
!  ------------

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, MAX_TAU_SPATH, MAX_TAU_UPATH, &
                              MAX_LAYERS, MAX_USER_STREAMS, MAX_POINTS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Control INTEGER ::s

      INTEGER  , INTENT(IN) :: NPOINTS, NLAYERS

!  Control flag and inputs for user streams

      LOGICAL  , INTENT(IN) :: DO_PLANE_PARALLEL
      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS
      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_STREAMS  ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: COS_SZA

!  Chapman factors

      REAL(FPK), INTENT(IN) :: CHAPMAN ( MAX_LAYERS, MAX_LAYERS )

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT   ( MAX_LAYERS, MAX_POINTS )

!  Output arguments
!  ----------------

!  Solar beam transmittances, average secant factors

      INTEGER  , INTENT(OUT) :: BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: BEAM_ETRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(OUT) :: SAVE_TRANS_USERM &
                ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Local variables
!  ---------------

      INTEGER   :: N, K, S, UM, CUTOFF
      REAL(FPK) :: SPHER, SUM, LOCAL_LTRANS, TAU_SLANT, ST0, ST1
      REAL(FPK) :: TAUSLANT(0:MAX_LAYERS), DELTA

!  Beam attenuation
!  ----------------

!  Monochromatic - 234 points (233 shifts + 1 calculation point)
!  Binned        - all points in outer buffer

!      IF ( DO_BIN_REALIZATION ) THEN
!        NPOINTS_LOCAL = NPOINTS_OUTER
!      ELSE
!        NPOINTS_LOCAL = NPOINTS_MONO
!      ENDIF

!  start loop over all outer points

      DO S = 1, NPOINTS

!  Tau slant = Beam slant optical thickness

        TAUSLANT(0) = ZERO
        DO N = 1, NLAYERS
          SUM = ZERO
          DO K = 1, N
            SUM = SUM + DELTAU_VERT_INPUT(K,S) * CHAPMAN(N,K)
          ENDDO
          TAUSLANT(N) = SUM
        ENDDO

!  Initial transmittances, and average secants
!   This section lifted from the LIDORT code

        ST0    = ONE
        ST1    = ST0
        CUTOFF = NLAYERS
        DO N = 1, NLAYERS
          DELTA = DELTAU_VERT_INPUT(N,S)
          IF ( N.LE.CUTOFF) THEN
            IF ( TAUSLANT(N) .GT. MAX_TAU_SPATH ) THEN
              CUTOFF = N
            ELSE
              ST1 = EXP ( - TAUSLANT(N) )
            ENDIF
            IF ( DO_PLANE_PARALLEL ) THEN
              BEAM_AVSECANT(N,S) = ONE / COS_SZA
            ELSE
              BEAM_AVSECANT(N,S) = (TAUSLANT(N)-TAUSLANT(N-1))/DELTA
            ENDIF
            BEAM_ITRANS(N,S)   = ST0
            ST0                = ST1
          ELSE
            BEAM_AVSECANT(N,S) = ZERO
            BEAM_ITRANS(N,S)   = ZERO
          ENDIF
        ENDDO
        BEAM_PICUTOFF(S) = CUTOFF

!  Bottom of layer transmittances, ETRANS

        DO N = 1, NLAYERS
          IF ( N.GT.BEAM_PICUTOFF(S) ) THEN
            LOCAL_LTRANS = ZERO
          ELSE
            TAU_SLANT    = DELTAU_VERT_INPUT(N,S) * BEAM_AVSECANT(N,S)
            LOCAL_LTRANS = EXP ( - TAU_SLANT )
          ENDIF
          BEAM_ETRANS(N,S) = BEAM_ITRANS(N,S) * LOCAL_LTRANS
          BEAM_DTRANS(N,S) = LOCAL_LTRANS
        ENDDO

!  End loop over points

      ENDDO

!  whole layer transmittances along user streams
!  =============================================

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO S = 1, NPOINTS
            DO N = 1, NLAYERS
              SPHER = DELTAU_VERT_INPUT(N,S) / USER_STREAMS(UM)
              IF ( SPHER.GT.MAX_TAU_UPATH ) THEN
                SAVE_TRANS_USERM(UM,N,S) = ZERO
              ELSE
                SAVE_TRANS_USERM(UM,N,S) = EXP ( - SPHER )
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE ATTENUATION_SETUP_1

!

      SUBROUTINE LEGENDRE_SETUPS &
            ( FOURIER,  NSTREAMS, NMOMENTS,                  & ! Input
              DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS, & ! Input
              COS_SZA, QUAD_STREAMS, QUAD_WEIGHTS,           & ! Input
              PLM_PXI,     PLM_MXI,     & ! Output
              PLM_PXUI,    PLM_MXUI,    & ! Output
              PLM_00_PXI,  PLM_00_MXI,  & ! Output
              PLM_00_PXUI, PLM_00_MXUI, & ! Output
              PLM_WT_PXI,  PLM_WT_MXI )   ! Output

!  Legendre polynomials and assocaited quantities

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, HALF, MAX_STREAMS, MAX_USER_STREAMS, MAX_MOMENTS 

!  Legendre subroutine

      USE LRRS_AUX2_m, Only : CFPLGARR

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  Fourier number, number of streams and moments

      INTEGER, INTENT(IN) ::    FOURIER
      INTEGER, INTENT(IN) ::    NSTREAMS
      INTEGER, INTENT(IN) ::    NMOMENTS

!  Number of user streams and flag

      INTEGER, INTENT(IN) ::    N_USER_STREAMS
      LOGICAL, INTENT(IN) ::    DO_USER_STREAMS

!  Cosine of the solar zenith angle

      REAL(FPK), INTENT(IN) :: COS_SZA

!  quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STREAMS ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS ( MAX_STREAMS )

!  user stream directions

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  output arguments
!  ----------------

      REAL(FPK), INTENT(OUT) :: PLM_PXI    ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_MXI    ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_00_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_00_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_WT_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_WT_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

      REAL(FPK), INTENT(OUT) :: PLM_PXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_MXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_00_PXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_00_MXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  local variables

      REAL(FPK) :: CFPLG  ( 0:MAX_MOMENTS )
      REAL(FPK) :: PLM_00 ( 0:MAX_MOMENTS )
      INTEGER :: I, UI, L, LPM

!  Legendre polynomials for quadrature streams
!  -------------------------------------

!  .. positive computational angle streams
!  .. negative streams by symmetry relations

      DO I = 1, NSTREAMS
        CALL CFPLGARR ( MAX_MOMENTS, NMOMENTS, FOURIER, QUAD_STREAMS(I), CFPLG)
        DO L = FOURIER, NMOMENTS
          PLM_PXI(I,L) = CFPLG(L)
        ENDDO
      ENDDO
      DO L = FOURIER, NMOMENTS
        LPM = L + FOURIER
        IF (MOD(LPM,2).EQ.0) THEN
          DO I = 1, NSTREAMS
            PLM_MXI(I,L) = PLM_PXI(I,L)
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            PLM_MXI(I,L) = - PLM_PXI(I,L)
          ENDDO
        ENDIF
      ENDDO

!  Legendre polynomials for User-angle streams
!  -------------------------------------

      IF ( DO_USER_STREAMS ) THEN
        DO UI = 1, N_USER_STREAMS
          CALL CFPLGARR ( MAX_MOMENTS, NMOMENTS, FOURIER, USER_STREAMS(UI), CFPLG)
          DO L = FOURIER, NMOMENTS
            PLM_PXUI(UI,L) = CFPLG(L)
          ENDDO
        ENDDO
        DO L = FOURIER, NMOMENTS
          LPM = L + FOURIER
          IF (MOD(LPM,2).EQ.0) THEN
            DO UI = 1, N_USER_STREAMS
              PLM_MXUI(UI,L) = PLM_PXUI(UI,L)
            ENDDO
          ELSE
            DO UI = 1, N_USER_STREAMS
              PLM_MXUI(UI,L) = - PLM_PXUI(UI,L)
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Legendre polynomials for solar streams
!  --------------------------------

!  .. negative solar zenith angle cosine

      CALL CFPLGARR ( MAX_MOMENTS, NMOMENTS, FOURIER, COS_SZA, PLM_00 )

!  Legendre polynomial products
!  ----------------------------

      DO L = FOURIER, NMOMENTS
        DO I = 1, NSTREAMS
          PLM_00_PXI(I,L) = PLM_00(L) * PLM_PXI(I,L)
          PLM_00_MXI(I,L) = PLM_00(L) * PLM_MXI(I,L)
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            PLM_00_PXUI(UI,L) = PLM_00(L) * PLM_MXUI(UI,L)
            PLM_00_MXUI(UI,L) = PLM_00(L) * PLM_PXUI(UI,L)
          ENDDO
        ENDIF
        DO I = 1, NSTREAMS
          PLM_WT_PXI(I,L) = QUAD_WEIGHTS(I) * PLM_PXI(I,L) * HALF
          PLM_WT_MXI(I,L) = QUAD_WEIGHTS(I) * PLM_MXI(I,L) * HALF
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LEGENDRE_SETUPS

!

      SUBROUTINE DBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,             & ! input
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,             & ! input
             FOURIER, Raman_IDX, Brdf_IDX, Sleave_IDX,        & ! Input
             NSTREAMS, DELTA_FACTOR, FLUX_FACTOR,             & ! input
             COS_SZA, BEAM_ETRANS_OPDEP, ALBEDOS_RANKED,      & ! input
             BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0,          & ! input
             ATMOS_ATTN, DIRECT_BEAM )                          ! Output

!  @@@@@@@@@@@@ Rob Fix 9/9/15 @@@@@@@@@@@@@@@@@@@@
!    Some re-writes here, introduction of LIDORT-style BRDF and SLEAVE control and inputs
!   
!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, only : fpk, MAX_STREAMS, MAX_LAYERS, MAX_MOMENTS,      &
                              MAX_POINTS, ZERO, PI4, FOUR

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  surface flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  New Surface-Leaving stuff

      LOGICAL, intent(in)    :: DO_SURFACE_LEAVING
      LOGICAL, intent(in)    :: DO_SL_ISOTROPIC

!  Add SLTERM flag and Isotropic term. @@@ Rob Fix 06 Sep 12. Superceded, Version 2.5
!      LOGICAL  , INTENT(IN)  :: DO_LRRS_SLTERM
!      REAL(FPK), INTENT(IN)  :: LRRS_SLTERM

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER

!  Point indices. Expanded for Version 2.5.
!     Raman_IDX  = APT/CPT for the Albedos_Ranked choice.
!     Brdf_IDX   = 1 (if single-point BRDF in use) or = Raman_Idx (Multiple-point BRDFs)
!     Sleave_IDX = 1 (if single-point SL   in use) or = Raman_Idx (Multiple-point SLs)

      INTEGER  , intent(in)  :: Raman_IDX, Brdf_IDX, Sleave_IDX

!  Computational variable

      INTEGER  , INTENT(IN)  :: NSTREAMS

!  Solar Flux

      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Surface factor

      REAL(fpk), intent(in)  :: DELTA_FACTOR

!  Solar zenith angle cosine

      REAL(FPK), INTENT(IN)  :: COS_SZA

!  Optical thickness value

      REAL(FPK), INTENT(IN)  :: BEAM_ETRANS_OPDEP

!  albedo

      REAL(fpk), intent(in)  :: ALBEDOS_RANKED ( MAX_POINTS )

!  BRDF Fourier components
!    incident solar direction,   reflected quadrature streams

      REAL(fpk), intent(in) :: BRDF_F_0      ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), intent(in) :: SLTERM_ISOTROPIC ( MAX_POINTS )

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams

      REAL(fpk), intent(in) :: SLTERM_F_0 ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )

!  Output Direct beam module
!  -------------------------

!  Direct beam itself

      REAL(FPK), INTENT(OUT) :: DIRECT_BEAM  ( MAX_STREAMS )

!  Attenuation

      REAL(FPK), INTENT(OUT) :: ATMOS_ATTN

!  Local variables
!  ---------------

      REAL(FPK) :: X0_FLUX, REFL_ATTN, HELP, SL
      INTEGER   :: I

!  Initialize
!  ----------

!   Safety first!  Return if there is no reflection.

      ATMOS_ATTN = ZERO ; DIRECT_BEAM = ZERO
      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Attenuation of solar beam
!  -------------------------

!  Bug fixed 23 November 2005, by R. Spurr
!    ( First noticed for the LIDORT codes by Italians!)

!     Old code only worked for FLUX_FACTOR = 1
!      X0_FLUX    = FOUR * COS_SZA * FLUX_FACTOR / DELTA_FACTOR
!    New code: Drop the FLUX_FACTOR in the next statement:

      X0_FLUX    = FOUR * COS_SZA / DELTA_FACTOR
!mick fix 7/20/2016 - added this normalization line (as in LIDORT) - on hold!
      !X0_FLUX    = FLUX_FACTOR * X0_FLUX / PI4
      ATMOS_ATTN = X0_FLUX * BEAM_ETRANS_OPDEP

!  Compute direct beam

      IF ( DO_BRDF_SURFACE ) THEN
        DO I = 1, NSTREAMS
          DIRECT_BEAM(I) = ATMOS_ATTN * BRDF_F_0(FOURIER,I,Brdf_Idx)
        ENDDO
      ELSE
!mick fix 7/20/2016 - removed Fourier condition
        !IF ( FOURIER.eq.0) then
          REFL_ATTN  = ATMOS_ATTN * ALBEDOS_RANKED(Raman_Idx)
          DIRECT_BEAM(1:NSTREAMS) = REFL_ATTN
        !ENDIF
      ENDIF

!  Comments from the initial implementation (2012): ======

!  Add SLTERM flag and Isotropic SLTERM. @@@ RobFix 06 sep 12
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic Fourier = 0 case
!    Flux_Factor  = 1.0 for All calculations
!   Previous kludge code used PI4:
!            HELP = PI4 * FLUX_FACTOR / DELTA_FACTOR --> HELP = PI4
!
!  NOTE TO SELF. Check the PI4 factor 9/10/15.

      IF ( DO_SURFACE_LEAVING ) THEN
        HELP = PI4 * FLUX_FACTOR / DELTA_FACTOR
        IF ( DO_SL_ISOTROPIC .and. FOURIER.EQ.0 ) THEN
          SL = SLTERM_ISOTROPIC(Sleave_idx) * HELP
          DO I = 1, NSTREAMS
            DIRECT_BEAM(I) = DIRECT_BEAM(I) + SL
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            SL = SLTERM_F_0(FOURIER,I,Sleave_idx) * HELP
            DIRECT_BEAM(I) = DIRECT_BEAM(I) + SL
          ENDDO
        ENDIF
      ENDIF

!  end direct beam calculation

      RETURN
      END SUBROUTINE DBEAM_SETUP

!

      SUBROUTINE UDBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,       & ! input
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,       & ! input
             DO_DB_CORRECTION, FOURIER,                 & ! Input
             Raman_IDX, Brdf_IDX, Sleave_IDX,           & ! Input
             N_USER_STREAMS, DELTA_FACTOR, FLUX_FACTOR, & ! input
             ATMOS_ATTN, ALBEDOS_RANKED, USER_BRDF_F_0, & ! input
             SLTERM_ISOTROPIC, USER_SLTERM_F_0,         & ! input
             USER_DIRECT_BEAM )                           ! Output

!  @@@@@@@@@@@@ Rob Fix 9/9/15 @@@@@@@@@@@@@@@@@@@@
!    Some re-writes here, introduction of LIDORT-style BRDF and SLEAVE control and inputs
!   
!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, only : fpk, MAX_USER_STREAMS, MAX_LAYERS, MAX_MOMENTS, &
                              MAX_POINTS, ZERO, PI4

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  surface flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  New Surface-Leaving stuff

      LOGICAL, intent(in)    :: DO_SURFACE_LEAVING
      LOGICAL, intent(in)    :: DO_SL_ISOTROPIC

!  Add SLTERM flag and Isotropic term. @@@ Rob Fix 06 Sep 12. Superceded, Version 2.5
!      LOGICAL  , INTENT(IN)  :: DO_LRRS_SLTERM
!      REAL(FPK), INTENT(IN)  :: LRRS_SLTERM

!  DB correction flag
!   SS now incorporates the exact direct bounce (DB term), at least in Version 2.3 !!!!!

      LOGICAL  , INTENT(IN)  ::   DO_DB_CORRECTION

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER

!  Point indices. Expanded for Version 2.5.
!     Raman_IDX  = APT/CPT for the Albedos_Ranked choice.
!     Brdf_IDX   = 1 (if single-point BRDF in use) or = Raman_Idx (Multiple-point BRDFs)
!     Sleave_IDX = 1 (if single-point SL   in use) or = Raman_Idx (Multiple-point SLs)

      INTEGER  , intent(in)  :: Raman_IDX, Brdf_IDX, Sleave_IDX

!  Computational variable

      INTEGER  , INTENT(IN)  :: N_USER_STREAMS

!  Solar Flux

      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Surface factor

      REAL(fpk), intent(in)  :: DELTA_FACTOR

!  Attenuation input

      REAL(FPK), INTENT(IN)  :: ATMOS_ATTN

!  albedo

      REAL(fpk), intent(in)  :: ALBEDOS_RANKED ( MAX_POINTS )

!  BRDF Fourier components
!    incident solar direction,   reflected user streams

      REAL(fpk), intent(in)  :: USER_BRDF_F_0 ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), intent(in)  :: SLTERM_ISOTROPIC ( MAX_POINTS )

!  Fourier components of Surface-leaving terms:
!    solar direction, SL-transmitted user streams

      REAL(fpk), intent(in)  :: USER_SLTERM_F_0 ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )

!  output Direct beam module
!  -------------------------

      REAL(FPK), INTENT(OUT) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS )

!  Local variables
!  ---------------

      REAL(FPK) :: REFL_ATTN, HELP, SL
      INTEGER   :: I

!  Initialize
!  ----------

!   Safety first!  Return if there is no reflection.

      USER_DIRECT_BEAM = ZERO
      IF ( .NOT. DO_INCLUDE_SURFACE  ) RETURN

!  Compute direct beam

      IF ( DO_BRDF_SURFACE ) THEN
        IF ( .not. DO_DB_CORRECTION ) THEN
          DO I = 1, N_USER_STREAMS
            USER_DIRECT_BEAM(I) = ATMOS_ATTN * USER_BRDF_F_0(FOURIER,I,Brdf_Idx)
          ENDDO
        ENDIF
      ELSE
!mick fix 7/20/2016 - removed Fourier condition
        !IF ( FOURIER.eq.0 .and. .not.DO_DB_CORRECTION ) then
        IF ( .not.DO_DB_CORRECTION ) then
          REFL_ATTN  = ATMOS_ATTN * ALBEDOS_RANKED(Raman_Idx)
          USER_DIRECT_BEAM(1:N_USER_STREAMS) = REFL_ATTN
        ENDIF
      ENDIF

!  Comments from the initial implementation (2012): ======

!  Add SLTERM flag and Isotropic SLTERM. @@@ RobFix 06 sep 12
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic Fourier = 0 case
!    Flux_Factor  = 1.0 for All calculations
!   Previous kludge code used PI4:
!            HELP = PI4 * FLUX_FACTOR / DELTA_FACTOR --> HELP = PI4
!
!  NOTE TO SELF. Check the PI4 factor 9/10/15.

      IF ( DO_SURFACE_LEAVING ) THEN
        HELP = PI4 * FLUX_FACTOR / DELTA_FACTOR
        IF ( DO_SL_ISOTROPIC .and. FOURIER.EQ.0 ) THEN
          SL = SLTERM_ISOTROPIC(Sleave_idx) * HELP
          DO I = 1, N_USER_STREAMS
            USER_DIRECT_BEAM(I) = USER_DIRECT_BEAM(I) + SL
          ENDDO
        ELSE
          DO I = 1, N_USER_STREAMS
            SL = USER_SLTERM_F_0(FOURIER,I,Sleave_idx) * HELP
            USER_DIRECT_BEAM(I) = USER_DIRECT_BEAM(I) + SL
          ENDDO
        ENDIF
      ENDIF

!  end direct beam calculation

      RETURN
      END SUBROUTINE UDBEAM_SETUP

!

      SUBROUTINE CHAPMAN_FUNCTION &
           ( DO_PLANE_PARALLEL, NLAYERS,         & ! Inputs
             COS_SZA, EARTH_RADIUS, HEIGHT_GRID, & ! Inputs
             CHAPMAN_FACTORS )                     ! Outputs

!  This is the Chapman function calculation of the slant path
!  geometrical factors, defined as the ratios between the slant path
!  distances to the corresponding vertical distances.

!  You must specify the Earth_radius and the height grid in order
!  to make this work.

!  This is a straightforward geometrical calculation, and is only
!  valid for a NON-REFRACTIVE atmosphere.

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, MAX_LAYERS

      IMPLICIT NONE

!  Input arguemnts
!  ---------------

      LOGICAL  , INTENT(IN) :: DO_PLANE_PARALLEL
      INTEGER  , INTENT(IN) :: NLAYERS
      REAL(FPK), INTENT(IN) :: COS_SZA
      REAL(FPK), INTENT(IN) :: EARTH_RADIUS
      REAL(FPK), INTENT(IN) :: HEIGHT_GRID(0:MAX_LAYERS)

!  output arguemnts
!  ----------------

      REAL(FPK), INTENT(OUT) :: CHAPMAN_FACTORS(MAX_LAYERS,MAX_LAYERS)

!  Local variables
!  ---------------

      INTEGER   :: N, M
      REAL(FPK) :: GM_TOA, HELP1, HELP2
      REAL(FPK) :: H(0:MAX_LAYERS), DELZ(MAX_LAYERS)
      REAL(FPK) :: STH, CTH, DELS, S1, S0

!  get spherical optical depths
!  ----------------------------

!  Prepare spherical attenuation (shell geometry)

      IF ( .NOT.DO_PLANE_PARALLEL ) THEN

        GM_TOA = SQRT ( ONE - COS_SZA * COS_SZA )
        DO N = 0, NLAYERS
          H(N) = HEIGHT_GRID(N) + EARTH_RADIUS
        ENDDO
        DO N = 1, NLAYERS
          DELZ(N) = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        ENDDO

        DO N = 1, NLAYERS
          STH = GM_TOA * H(N)/H(0)
          CTH = SQRT ( ONE - STH * STH )
          S0 = ZERO
          HELP1 = H(0)*CTH
          HELP2 = -H(0)*H(0)*STH*STH
          DO M = 1, N
            S1 = HELP1 - SQRT(HELP2 + H(M)*H(M))
            DELS = S1 - S0
            CHAPMAN_FACTORS(N,M) = DELS / DELZ(M)
            S0 = S1
          ENDDO
        ENDDO

!  Plane parallel

      ELSE

        DO N = 1, NLAYERS
          DO M = 1, N
            CHAPMAN_FACTORS(N,M) = ONE / COS_SZA
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE CHAPMAN_FUNCTION

! End module

      END MODULE lrrs_miscsetups_1_m
