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
! #             RAMAN_FOURIER_1 (*)                             #
! #                                                             #
! #   Note (*): Version 2.1 with BRDF setup,   November 2007    #
! #   Note    : Version 2.2 with Linearization, November 2008   #
! #                                                             #
! ###############################################################


      MODULE lrrs_fourier_master_1

      PRIVATE
      PUBLIC :: RAMAN_FOURIER_1

      CONTAINS

      SUBROUTINE RAMAN_FOURIER_1 &
       ( DO_RRS_OVERALL, DO_DIRECTRRS_ONLY, DO_MSMODE_LRRS, PROGRESS, &
         DO_ELASTIC_ONLY, DO_SSCORR_OUTGOING, DO_SSFULL, &
         DO_BIN_REALIZATION, DO_ENERGY_BALANCING, &
         DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_DIRECT_BEAM, &
         DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, &
         DO_LAMBERTIAN_SURFACE, DO_WATERLEAVING, &
         DO_LAMBERTIAN_FRACTION, LAMBERTIAN_FRACTION, &
         NPOINTS_OUTER, OFFSET_INNER, NPOINTS_INNER, NPOINTS_MONO, &
         N_RRSBINS, BINMAP, W_EXCIT, N_BRDF, X_BRDF, A_BRDF, &
         BRDFUNC, BRDFUNC_0, USER_BRDFUNC, USER_BRDFUNC_0, &
         NLAYERS, NSTREAMS, NMOMENTS, N_LOUTPUT, N_USER_STREAMS, &
         FOURIER, COS_SZA, LRRS_FGSMALL, FLUX_FACTOR, &
         QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS, &
         STERM_MASK_UP, LEVELMASK_UP, STERM_MASK_DN, LEVELMASK_DN, &
         DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, OMEGAMOMS_CABANNES, &
         OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN, OMEGAMOMS_RRSGAIN, &
         ALBEDOS_RANKED, FLUXES_RANKED, &
         BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, &
         BEAM_DTRANS, BEAM_ETRANS, SAVE_TRANS_USERM, &
         ELASTIC_F_UP, MEAN_ELASTIC_UP, FLUX_ELASTIC_UP, &
         ELASTIC_F_DN, MEAN_ELASTIC_DN, FLUX_ELASTIC_DN, &
         RAMAN_F_UP, MEAN_RAMAN_UP, FLUX_RAMAN_UP, &
         RAMAN_F_DN, MEAN_RAMAN_DN, FLUX_RAMAN_DN, &
         FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE )

!  include file of dimensions and numbers

      USE LRRS_PARS
      USE LRRS_RTSOLUTIONS_1
      USE LRRS_ELASTIC_MASTER_1
      USE LRRS_POSTPROCESSING_1
      USE LRRS_BVPROBLEM
      USE LRRS_AUX2
      USE LRRS_SOURCES_MASTER_1
      USE LRRS_RAMAN_INTENSITY_1

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Progress number
!    If this is zero, no output to screen

      INTEGER, INTENT(IN) ::   PROGRESS

!  elastic control flag

      LOGICAL, INTENT(IN) ::   DO_ELASTIC_ONLY

!  Raman  flags for computing filling

      LOGICAL, INTENT(IN) ::   DO_RRS_OVERALL
      LOGICAL, INTENT(IN) ::   DO_DIRECTRRS_ONLY

!  Multiple scatter mode

      LOGICAL, INTENT(IN) ::   DO_MSMODE_LRRS

!  Overall SS correction flag; Now incorporates the DB term

      LOGICAL, INTENT(IN) ::   DO_SSCORR_OUTGOING
      LOGICAL, INTENT(IN) ::   DO_SSFULL

!  Upwelling and downwelling flags

      LOGICAL, INTENT(IN) ::   DO_UPWELLING
      LOGICAL, INTENT(IN) ::   DO_DNWELLING

!  User stream flag

      LOGICAL, INTENT(IN) ::   DO_USER_STREAMS

!  Binning realization

      LOGICAL, INTENT(IN) ::   DO_BIN_REALIZATION

!  Energy balancing flag

      LOGICAL, INTENT(IN) ::   DO_ENERGY_BALANCING

!  Mean value output

      LOGICAL, INTENT(IN) ::   DO_MVOUT_ONLY
      LOGICAL, INTENT(IN) ::   DO_ADDITIONAL_MVOUT

!  Direct beam

      LOGICAL, INTENT(IN) ::   DO_DIRECT_BEAM

!  Surface control
!  ---------------

!  overall lambertian surface

      LOGICAL, INTENT(IN) ::   DO_LAMBERTIAN_SURFACE

!  Lambertian admixture control

      LOGICAL, INTENT(IN) ::   DO_LAMBERTIAN_FRACTION
      LOGICAL, INTENT(IN) ::   DO_WATERLEAVING
      REAL(FPK), INTENT(IN) :: LAMBERTIAN_FRACTION

!  BRDF control and quadrature

      INTEGER, INTENT(IN) ::   N_BRDF
      REAL(FPK), INTENT(IN) :: X_BRDF ( MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(IN) :: A_BRDF ( MAX_STREAMS_BRDF )

!  BRDFs at quadrature (discrete ordinate) angles

      REAL(FPK), INTENT(IN) :: BRDFUNC &
            ( MAX_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK), INTENT(IN) :: BRDFUNC_0 &
            ( MAX_STREAMS, MAX_STREAMS_BRDF )

!  BRDFs at user-defined stream directions

      REAL(FPK), INTENT(IN) :: USER_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK), INTENT(IN) :: USER_BRDFUNC_0 &
            ( MAX_USER_STREAMS, MAX_STREAMS_BRDF )

!  Numbers
!  -------

!  number of streams

      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS

!  Number of moments

      INTEGER, INTENT(IN) ::   NMOMENTS

!  number of layers

      INTEGER, INTENT(IN) ::   NLAYERS

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)

      INTEGER, INTENT(IN) ::   N_LOUTPUT

!  Fourier number

      INTEGER, INTENT(IN) ::   FOURIER

!  Solar beam input Cosine

      REAL(FPK), INTENT(IN) :: COS_SZA

!  Small number control

      REAL(FPK), INTENT(IN) :: LRRS_FGSMALL

!  Flux factor

      REAL(FPK), INTENT(IN) :: FLUX_FACTOR

!  Streams
!  -------

!  Discrete ordinate quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STREAMS ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUAD_STRMWGT ( MAX_STREAMS )

!  User stream variables

      REAL(FPK), INTENT(IN) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  Number of points
!  ----------------

      INTEGER, INTENT(IN) ::   NPOINTS_OUTER
      INTEGER, INTENT(IN) ::   NPOINTS_INNER
      INTEGER, INTENT(IN) ::   OFFSET_INNER
      INTEGER, INTENT(IN) ::   NPOINTS_MONO

      INTEGER, INTENT(INOUT) :: W_EXCIT

!  Layer masks for doing integrated source terms
!  ---------------------------------------------

      LOGICAL, INTENT(IN) ::   STERM_MASK_UP ( MAX_LAYERS )
      LOGICAL, INTENT(IN) ::   STERM_MASK_DN ( MAX_LAYERS )

      INTEGER, INTENT(IN) ::   LEVELMASK_UP    (MAX_LOUTPUT)
      INTEGER, INTENT(IN) ::   LEVELMASK_DN    (MAX_LOUTPUT)

!  Scaled elastic-scattering optical properties
!  --------------------------------------------

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT ( MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC &
        ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived Cabannes scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_CABANNES &
        ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Scaled Raman-scattering optical properties
!  ------------------------------------------

!  Bin Mapping quantities (internal derivation)

      INTEGER, INTENT(IN) ::   N_RRSBINS  ( MAX_POINTS )
      INTEGER, INTENT(IN) ::   BINMAP ( MAX_BINS, MAX_POINTS )

!  Derived rotational raman scattering (Internal variable). Loss term
!    Single scatter albedo * phase function moments
!    Only required for the Energy-balance approximation.

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSLOSS &
        ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived rotational Raman scattering (Internal variable). BINNING

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSBIN &
        ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Derived rotational Raman scattering (Internal variable). MONO

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSGAIN &
        ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Albedos and fluxes

      REAL(FPK), INTENT(IN) :: ALBEDOS_RANKED ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: FLUXES_RANKED  ( MAX_POINTS )

!  Beam Attenuations
!  -----------------

!  Solar beam transmittances, average secant factors

      INTEGER, INTENT(IN) ::   BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_ETRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(IN) :: SAVE_TRANS_USERM &
        (  MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  OUTPUT ARGUMENTS START HERE
!  ===========================

!  Elastic output at User angles, for one Fourier component
!  --------------------------------------------------------

!  Fourier output

      REAL(FPK), INTENT(OUT) :: ELASTIC_F_UP &
        ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: ELASTIC_F_DN &
        ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  Mean-value output

      REAL(FPK), INTENT(INOUT) :: &
        MEAN_ELASTIC_UP ( MAX_LOUTPUT, MAX_POINTS ), &
        MEAN_ELASTIC_DN ( MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: FLUX_ELASTIC_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: FLUX_ELASTIC_DN ( MAX_LOUTPUT, MAX_POINTS )

!  Raman output at User angles, for one Fourier component
!  ------------------------------------------------------

!  Fourier output

      REAL(FPK), INTENT(OUT) :: RAMAN_F_UP &
        ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: RAMAN_F_DN &
        ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  Mean-value output

      REAL(FPK), INTENT(INOUT) :: &
        MEAN_RAMAN_UP( MAX_LOUTPUT, MAX_POINTS ), &
        MEAN_RAMAN_DN( MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: FLUX_RAMAN_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: FLUX_RAMAN_DN ( MAX_LOUTPUT, MAX_POINTS )

!  module status:
!   Message_sub comes from any of the subroutines called here.
!   Point_trace gives a trace of subroutine where message occurred

      LOGICAL, INTENT(OUT) ::            FAIL
      CHARACTER (LEN=70), INTENT(OUT) :: MESSAGE
      CHARACTER (LEN=70), INTENT(OUT) :: MESSAGE_SUB
      CHARACTER (LEN=70), INTENT(OUT) :: POINT_TRACE

!  LOCAL VARIABLES START HERE
!  ==========================

!  Fourier components of BRDF functions
!  ------------------------------------

!  at quadrature (discrete ordinate) angles

      REAL(FPK) :: BIREFLEC   ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: BIREFLEC_0 ( MAX_STREAMS  )

!  at user-defined stream directions

      REAL(FPK) :: USER_BIREFLEC ( MAX_USER_STREAMS, MAX_STREAMS )
      REAL(FPK) :: USER_BIREFLEC_0 ( MAX_USER_STREAMS  )

!  Legendre polynomial variables
!  -----------------------------

!  Legendre functions at discrete ordinates

      REAL(FPK) :: PLM_PXI    ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK) :: PLM_MXI    ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Product functions, discrete ordinates and solar stream

      REAL(FPK) :: PLM_00_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK) :: PLM_00_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Product of Legendre functions with quadrature weights

      REAL(FPK) :: PLM_WT_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK) :: PLM_WT_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Legendre functions at user defined polar streams

      REAL(FPK) :: PLM_PXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK) :: PLM_MXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  Product functions, User streams and solar stream

      REAL(FPK) :: PLM_00_PXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK) :: PLM_00_MXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  Elastic discrete ordinate solutions at shifted wavelengths
!  ----------------------------------------------------------

!  Eigenvalues

      REAL(FPK) :: L0_KEIGEN &
        ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Eigentransmittances

      REAL(FPK) :: L0_KTRANS &
        ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Particular solution vectors

      REAL(FPK) :: L0_WPARTIC &
        ( MAX_2_STREAMS, MAX_LAYERS, MAX_POINTS )

!  homogeneous solution vectors and integration constants

      REAL(FPK) :: L0_WHOM_XPOS &
        ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L0_WHOM_LCON &
        ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L0_WHOM_MCON &
        ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  output Direct beam module
!  -------------------------

!  Reflected Direct beam at surface (quadrature + User streams)

      REAL(FPK) :: DIRECT_BEAM      ( MAX_STREAMS )
      REAL(FPK) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS )

!  solar beam optical thickness and atmospheric attenuation

!      REAL(FPK) :: SOLAR_BEAM_OPDEP
      REAL(FPK) :: ATMOS_ATTN

!  output arguments from Homogeneous solution modules
!  --------------------------------------------------

!  (Positive) Eigenvalues

      REAL(FPK) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!   Whole layer transmittance factors for +/- eigenvalues

      REAL(FPK) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )

!  Eigenvector solutions

      REAL(FPK) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Saved matrices for eigenvalue computation

      REAL(FPK) :: DAB_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: SAB_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: EIGENMAT_SAVE  ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: EIGENVEC_SAVE  ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: DIFVEC_SAVE    ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: EIGENNORM_SAVE ( MAX_STREAMS, MAX_LAYERS )

!  Reflected homogeneous solutions at ground

      REAL(FPK) :: R2_HOMP(MAX_STREAMS,MAX_STREAMS)
      REAL(FPK) :: R2_HOMM(MAX_STREAMS,MAX_STREAMS)

!  Local arrays
!  ------------

      REAL(FPK) :: LOCAL_DELTAUS ( MAX_LAYERS )
      REAL(FPK) :: LOCAL_TRANS_USERM(MAX_USER_STREAMS,MAX_LAYERS)

!  Multipliers
!  -----------

!  user-stream homogeneous solution vectors

      REAL(FPK) :: &
        U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
        U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Whole layer multipliers (homogeneous)

      REAL(FPK) :: &
        HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
        HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Zeta functions (homogeneous) = 1 / [ 1/mu) +/- k_a ]

      REAL(FPK) :: &
        ZETA_P ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
        ZETA_M ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  output arguments from Beam solution modules
!  -------------------------------------------

!  Particular solutions at upper and lower boundaries
!   (as evaluated for the inelastic field by Green's function method)

!      REAL(FPK) :: RAMAN_WUPPER  ( MAX_LAYERS, MAX_2_STREAMS )
!      REAL(FPK) :: RAMAN_WLOWER  ( MAX_LAYERS, MAX_2_STREAMS )

      REAL(FPK) :: RAMAN_WUPPER  ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK) :: RAMAN_WLOWER  ( MAX_2_STREAMS, MAX_LAYERS )

!  Reflected beam solutions at ground

      REAL(FPK) :: R2_BEAM(MAX_STREAMS)

!  BV problem variables
!  --------------------

!  band matrix and indexingmask, column vector COL2

      REAL(FPK) :: COL2(MAX_TOTAL,1)
      REAL(FPK) :: BANDMAT2(MAX_BANDTOTAL,MAX_TOTAL)
      INTEGER ::   BMAT_ROWMASK(MAX_TOTAL,MAX_TOTAL)

!  LU decomposition PIVOT

      INTEGER ::   IPIVOT ( MAX_TOTAL )

!  Integration constants

      REAL(FPK) :: &
        LCON ( MAX_STREAMS, MAX_LAYERS ), &
        MCON ( MAX_STREAMS, MAX_LAYERS )

!  Integration constants x homogeneous solution vectors

      REAL(FPK) :: &
        LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
        MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Post-processed source terms
!  ---------------------------

      REAL(FPK) :: WLAYER_PIST_UP ( MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK) :: WLAYER_PIST_DN ( MAX_LAYERS, MAX_USER_STREAMS )

      REAL(FPK) :: WLAYER_HOST_UP ( MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK) :: WLAYER_HOST_DN ( MAX_LAYERS, MAX_USER_STREAMS )

      REAL(FPK) :: TOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(FPK) :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(FPK) :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )

      REAL(FPK) :: CUMSOURCE_UP ( MAX_USER_STREAMS, 0:MAX_LAYERS )
      REAL(FPK) :: CUMSOURCE_DN ( MAX_USER_STREAMS, 0:MAX_LAYERS )

!  Local reflected field (not used here) is the diffuse BOA source term

      REAL(FPK) :: IDOWNSURF ( MAX_STREAMS )

!  Help arrays for linearization (not required in this routine)

      REAL(FPK) :: U_HELP_P &
        ( MAX_STREAMS, 0:MAX_MOMENTS, MAX_LAYERS )
      REAL(FPK) :: U_HELP_M &
        ( MAX_STREAMS, 0:MAX_MOMENTS, MAX_LAYERS )

!  Miscellaneous
!  -------------

!  Fourier output (quadrature streams)
!    Use as debug, except for Flux calculations

      REAL(FPK) :: QUADRAMAN_F_UP &
        ( MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

      REAL(FPK) :: QUADRAMAN_F_DN &
        ( MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

!  Surface contributions

      REAL(FPK) :: FL1 ( MAX_POINTS )
      REAL(FPK) :: FL2 ( MAX_POINTS )

!  Log-file output (LU = 6 gives screen output)

      LOGICAL, PARAMETER :: LG = .false.
      INTEGER, PARAMETER :: LU = 8

!  for exception handling

      CHARACTER (LEN=3) ::  C3
      CHARACTER (LEN=2) ::  C2

!  Local integers

      INTEGER ::   NTOTAL, N_SUPDIAG, N_SUBDIAG, N_GAIN_TERMS
      INTEGER ::   NSTR2, N_SOLUTIONS, PROG, M, INFO
      INTEGER ::   LAYER, I, I1, AA, N, C0, APT, CPT, L, UM, W
      INTEGER ::   NPOINTS_TOTAL, NPOINTS_LOCAL

!  Number of RRS solutions

      INTEGER ::   N_LOSS_SOLUTIONS
      INTEGER ::   N_GAIN_SOLUTIONS

!  Monochromatic realization

      LOGICAL ::   DO_MONO_REALIZATION

!  Cabannes-Raman flag

      LOGICAL ::   DO_CABANNES_RAMAN

!  Flags for including surface term

      LOGICAL ::   DO_INCLUDE_SURFACE
      LOGICAL ::   DO_REFLECTED_DIRECTBEAM

!  Flag for mean value output

      LOGICAL ::   DO_INCLUDE_MVOUT

!  Flag for post-processing of elastic field

      LOGICAL ::   DO_ELASTIC_POSTPROCESSING

!  Local flag for CB COPY

      LOGICAL ::   DO_CBCOPY

!  Surface factors and local XT, DT

      REAL(FPK) :: SURFACE_FACTOR, SURFACE_FACTOR_0
      REAL(FPK) :: DELTA_FACTOR, HOM, V1, V2, DETER
      REAL(FPK) :: NORM, T1, T2, LF

!  solar irradianace flux

      REAL(FPK) :: SS_FLUX_MULTIPLIER
      REAL(FPK) :: FLUX_MULTIPLIER

!  Local IOPs

      REAL(FPK) :: LOCAL_OMEGAMOMS ( MAX_LAYERS, 0: MAX_MOMENTS )

!  Overall flag for computing filling

      LOGICAL ::   DO_FILLING

!  ###############
!  INITIAL SECTION
!  ###############

!  set exclusion flags

      DO_MONO_REALIZATION = .not. DO_BIN_REALIZATION
      DO_CABANNES_RAMAN   = .not. DO_ENERGY_BALANCING

!  set overall filling flags

      DO_FILLING = .not.DO_ELASTIC_ONLY.and.DO_RRS_OVERALL

!  Fourier

      M = FOURIER

!  Npoints local

      IF ( DO_BIN_REALIZATION ) THEN
        NPOINTS_TOTAL = NPOINTS_OUTER
        NPOINTS_LOCAL = NPOINTS_INNER
      ELSE
        NPOINTS_TOTAL = NPOINTS_MONO
        NPOINTS_LOCAL = 1
      ENDIF

!  Legendre polynomial setup

      IF ( .NOT. DO_SSFULL ) THEN
        CALL LEGENDRE_SETUPS &
            ( FOURIER,  NSTREAMS, NMOMENTS, &
              DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS, &
              COS_SZA, QUAD_STREAMS, QUAD_WEIGHTS, &
              PLM_PXI,     PLM_MXI, &
              PLM_PXUI,    PLM_MXUI, &
              PLM_00_PXI,  PLM_00_MXI, &
              PLM_00_PXUI, PLM_00_MXUI, &
              PLM_WT_PXI,  PLM_WT_MXI )
      ENDIF

!  Zeroth-order code.
!  ==================

!  Get the Single scatter elastic correction
!    Added, 18 January 2007. RT Solutions Inc.
!      IF ( DO_SSCORRECTION .AND. FOURIER.EQ.0 ) THEN
!        CALL LRRS_SSCORRECTION_BEFORE
!      ENDIF
!  Get the Exact Direct Beam correction
!    Added, 28 November 2007. RT Solutions Inc.
!      IF ( DO_DB_CORRECTION.AND.DO_UPWELLING.AND.FOURIER.EQ.0 ) THEN
!        CALL LRRS_DBCORRECTION_BEFORE
!      ENDIF

!  Surface flag (for inclusion of some kind of reflecting boundary)
!    Fake-Lambertian test: set DO_INCLUDE_SURFACE only for M = 0

      DO_INCLUDE_SURFACE    = .TRUE.
      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        IF ( FOURIER .NE. 0 ) DO_INCLUDE_SURFACE = .FALSE.
      ENDIF

!  setup code
!  ==========

!  Direct beam flag (only if above surface flag has been set)

      DO_REFLECTED_DIRECTBEAM = .FALSE.
      IF ( DO_DIRECT_BEAM ) THEN
        IF ( DO_INCLUDE_SURFACE ) THEN
          DO_REFLECTED_DIRECTBEAM = .TRUE.
        ENDIF
      ENDIF

!  surface reflectance factors

      IF ( FOURIER .EQ. 0 ) THEN
        SURFACE_FACTOR_0 = TWO
        DELTA_FACTOR     = ONE
      ELSE
        SURFACE_FACTOR_0 = ONE
        DELTA_FACTOR     = TWO
      ENDIF

!  Flux multipliers

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4
      FLUX_MULTIPLIER    = SS_FLUX_MULTIPLIER * DELTA_FACTOR

!  set up of associated computational integers

      NSTR2     = 2 * NSTREAMS
      NTOTAL    = NSTR2 * NLAYERS
      IF ( NLAYERS .EQ. 1 ) THEN
        N_SUBDIAG = 2*NSTREAMS - 1
        N_SUPDIAG = 2*NSTREAMS - 1
      ELSE
        N_SUBDIAG = 3*NSTREAMS - 1
        N_SUPDIAG = 3*NSTREAMS - 1
      ENDIF

!  inclusion of mean value output
!   (set local quadrature results flag)

      DO_INCLUDE_MVOUT = .FALSE.
      IF ( DO_MVOUT_ONLY .OR. DO_ADDITIONAL_MVOUT ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_MVOUT = .TRUE.
        ENDIF
      ENDIF

!  Surface stuff
!  -------------

!  BRDF with Lambertian fraction; set up coefficients
!    BRDF Fourier components to be specified in the elastic-field call.
!    Note: extra code for Waterleaving option

!  Zero the arrays

      DO CPT = 1, NPOINTS_TOTAL
        FL1(CPT) = ZERO
        FL2(CPT) = ZERO
      ENDDO

!  lambertian surface

      IF ( DO_LAMBERTIAN_SURFACE .and. FOURIER .eq. 0 ) THEN
        LF = LAMBERTIAN_FRACTION
        DO CPT = 1, NPOINTS_TOTAL
          FL1(CPT) = ALBEDOS_RANKED(CPT)
        ENDDO
      ENDIF

!  BRDF surfaces
!    May include Lambertian fraction and/or waterleaving

      IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
        IF ( DO_LAMBERTIAN_FRACTION ) THEN
          IF ( FOURIER .EQ. 0 ) THEN
            IF ( DO_WATERLEAVING ) THEN
              DO CPT = 1, NPOINTS_TOTAL
                FL1(CPT) = ALBEDOS_RANKED(CPT)
              ENDDO
            ELSE
              DO CPT = 1, NPOINTS_TOTAL
                FL1(CPT) = LF * ALBEDOS_RANKED(CPT)
              ENDDO
            ENDIF
          ENDIF
          IF ( DO_WATERLEAVING ) THEN
            DO CPT = 1, NPOINTS_TOTAL
              FL2(CPT) = ONE
            ENDDO
          ELSE
            DO CPT = 1, NPOINTS_TOTAL
              FL2(CPT) = ONE - LF
            ENDDO
          ENDIF
        ELSE
          DO CPT = 1, NPOINTS_TOTAL
            FL2(CPT) = ONE
          ENDDO
        ENDIF
      ENDIF

!  ###################
!  ELASTIC CALCULATION
!  ###################

!  Post processing flag
!    MVOUT_ONLY added 06 May 2010

      DO_ELASTIC_POSTPROCESSING = .false.
      IF ( DO_USER_STREAMS ) THEN
        IF ( DO_FILLING .OR. DO_ELASTIC_ONLY ) THEN
          DO_ELASTIC_POSTPROCESSING = .true.
        ENDIF
      ELSE
        IF ( DO_MVOUT_ONLY ) DO_ELASTIC_POSTPROCESSING = .true.
      ENDIF

!  This is the coherent light solution with 'Rayleigh' scattering
!  Same for both Cabannes/Raman and Energy Balancing methods

      CALL ELASTIC_MASTER_1 &
       ( DO_ELASTIC_ONLY, DO_ELASTIC_POSTPROCESSING, DO_MSMODE_LRRS, &
         DO_INCLUDE_MVOUT, DO_INCLUDE_SURFACE, DO_SSCORR_OUTGOING, &
         DO_UPWELLING,DO_DNWELLING,DO_USER_STREAMS,DO_BIN_REALIZATION, &
         DO_LAMBERTIAN_SURFACE, FL1, FL2, N_BRDF, X_BRDF, A_BRDF, &
         BRDFUNC, BRDFUNC_0, USER_BRDFUNC, USER_BRDFUNC_0, &
         NLAYERS, NSTREAMS, NMOMENTS, N_LOUTPUT, N_USER_STREAMS, &
         FOURIER, COS_SZA, LRRS_FGSMALL, FLUX_FACTOR, W_EXCIT, &
         NPOINTS_OUTER, OFFSET_INNER, NPOINTS_INNER, NPOINTS_MONO, &
         QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS, &
         STERM_MASK_UP, LEVELMASK_UP, STERM_MASK_DN, LEVELMASK_DN, &
         DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, &
         BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, &
         BEAM_DTRANS,   BEAM_ETRANS, SAVE_TRANS_USERM, &
         PLM_PXI, PLM_MXI, PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI, PLM_MXUI, &
         PLM_00_PXI, PLM_00_MXI, PLM_00_PXUI, PLM_00_MXUI, &
         L0_KEIGEN,     L0_KTRANS,     L0_WPARTIC, &
         L0_WHOM_XPOS,  L0_WHOM_LCON,  L0_WHOM_MCON, &
         BIREFLEC, BIREFLEC_0, USER_BIREFLEC, USER_BIREFLEC_0, &
         ELASTIC_F_UP, MEAN_ELASTIC_UP, FLUX_ELASTIC_UP, &
         ELASTIC_F_DN, MEAN_ELASTIC_DN, FLUX_ELASTIC_DN, &
         FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE )
         
      ! !    write(*,*)NLAYERS, NSTREAMS, NMOMENTS, N_LOUTPUT, N_USER_STREAMS, &
      ! !    FOURIER, COS_SZA, LRRS_FGSMALL
      ! !    pause
      !    write(*,*)DELTAU_VERT_INPUT(1,:)
      ! !    write(*,*)OMEGAMOMS_ELASTIC(1,2,1:100)
      !    pause
   
!  debug, 24   September 2015

!     do w = 1, npoints_outer, 2
!       write(*,*)w,L0_WHOM_LCON(1:2,23,w),L0_WHOM_LCON(7:8,1,w)
!     enddo
!     pause'after    Elastic Master'

!  Return conditions

      IF ( FAIL ) RETURN
      IF ( DO_ELASTIC_ONLY ) RETURN

!  BRDF Fourier code already done in the elastic case
!  --------------------------------------------------

!  Here it is again to remind us.......................!!!!!!!!!!!!!!!!!

!  BRDF Fourier components - coefficient setup
!     Wavelength independent. No Raman scattering.
!   Surface Linearization added 24 November 2008.
!      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
!        IF ( DO_SURFACE_JACOBIANS.AND..NOT.DO_LAMBERT_VARIATION ) THEN
!          CALL LS_LRRS_BRDF_FOURIER
!     I    ( DO_INCLUDE_SURFACE,
!     I      FOURIER, DELFAC )
!        ELSE
!          CALL LRRS_BRDF_FOURIER
!     I    ( DO_INCLUDE_SURFACE,
!     I      FOURIER, DELFAC )
!        ENDIF
!      ENDIF

!  #################################################
!  INELASTIC CALCULATION - MAIN LOOP OVER ALL POINTS
!  #################################################

!  Progress...

!      IF ( NSTREAMS .GT. 10 ) PROG = 5
!      IF ( NSTREAMS .LE. 10 ) PROG = 30
!c      IF ( NSTREAMS .LE.  6 ) PROG = 10
!      IF ( NSTREAMS .LE.  6 ) PROG = 50

      IF ( PROGRESS .EQ. 0 ) THEN
        PROG = 60000
      ELSE
        PROG = PROGRESS
      ENDIF

!  Start points loop
!  -----------------

      DO W = 1, NPOINTS_LOCAL

!  For binning, set counters in inner (APT) and outer (CPT) loops

        IF ( DO_BIN_REALIZATION ) THEN
          APT = W
          CPT = APT + OFFSET_INNER
          IF ( MOD(APT,PROG).EQ.0 ) THEN
            WRITE(*,'(A,I4)')'    --> Doing Calculation point # ',APT
          ENDIF
          W_EXCIT = CPT
        ELSE
          APT = 1
          CPT = W_EXCIT
        ENDIF

!  surface factor

        SURFACE_FACTOR = SURFACE_FACTOR_0

!  Local layer optical depths and layer transmittances

        DO N = 1, NLAYERS
          LOCAL_DELTAUS(N) = DELTAU_VERT_INPUT(N,CPT)
          DO UM = 1, N_USER_STREAMS
            LOCAL_TRANS_USERM(UM,N) = SAVE_TRANS_USERM(UM,N,CPT)
          ENDDO
        ENDDO

!  local IOPs (omega-moments)
!    -Usual elastic-scatering IOPS, unless using Cabannes/Raman method

        DO_CBCOPY =  DO_RRS_OVERALL.AND..NOT.DO_ENERGY_BALANCING
        IF ( DO_CBCOPY ) THEN
         DO N = 1, NLAYERS
          DO L = 0, NMOMENTS
           LOCAL_OMEGAMOMS(N,L) = OMEGAMOMS_CABANNES(N,L,APT)
          ENDDO
         ENDDO
        ELSE
         DO N = 1, NLAYERS
          DO L = 0, NMOMENTS
           LOCAL_OMEGAMOMS(N,L) = OMEGAMOMS_ELASTIC(N,L,CPT)
          ENDDO
         ENDDO
        ENDIF

!  Number of solutions
!  ===================

!  Start with 1 for Elastic solution (no RRS yet)

        N_LOSS_SOLUTIONS = 0
        N_GAIN_SOLUTIONS = 0
        N_SOLUTIONS      = 1

!  Number of RRS gain terms

        IF ( DO_BIN_REALIZATION ) THEN
          N_GAIN_TERMS = N_RRSBINS(APT)
        ELSE
          N_GAIN_TERMS = NPOINTS_MONO - 1
        ENDIF

!  Number of solutions for the Cabannes/Raman solution (no loss terms)

        IF ( DO_CABANNES_RAMAN ) THEN
          IF ( DO_RRS_OVERALL .AND. FOURIER .LE. 2 ) THEN
            N_GAIN_SOLUTIONS = N_GAIN_TERMS
            IF ( .NOT. DO_DIRECTRRS_ONLY ) THEN
                N_GAIN_SOLUTIONS = N_GAIN_SOLUTIONS + &
                                  N_GAIN_TERMS * (2*NSTREAMS+1)
            ENDIF
            N_SOLUTIONS = 1 + N_LOSS_SOLUTIONS + N_GAIN_SOLUTIONS
          ENDIF
        ENDIF

!  Number of solutions for Energy balancing solution (includes loss term

        IF ( DO_ENERGY_BALANCING ) THEN
          IF ( DO_RRS_OVERALL .AND. FOURIER .LE. 2 ) THEN
            N_LOSS_SOLUTIONS = 1
            N_GAIN_SOLUTIONS = N_GAIN_TERMS
            IF ( .NOT. DO_DIRECTRRS_ONLY ) THEN
              N_LOSS_SOLUTIONS = N_LOSS_SOLUTIONS + 2*NSTREAMS+1
              N_GAIN_SOLUTIONS = N_GAIN_SOLUTIONS + &
                              N_GAIN_TERMS*(2*NSTREAMS+1)
            ENDIF
            N_SOLUTIONS = 1 + N_LOSS_SOLUTIONS + N_GAIN_SOLUTIONS
          ENDIF
        ENDIF

!  Reflected Direct beam attenuation
!  =================================

!       DO_SSCORR_OVERALL     ! formerly DO_DB_correction

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  Check that the BEAM_ETRANS(NLAYERS,CPT) is the right value
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        if (lg)write(lu,'(a)')'starting DBEAM_SETUP'
        CALL DBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, &
             NSTREAMS, FL1(CPT), FL2(CPT), BIREFLEC_0, &
             COS_SZA, FOURIER, DELTA_FACTOR, BEAM_ETRANS(NLAYERS,CPT), &
             ATMOS_ATTN, DIRECT_BEAM )

!  RT differential equation solutions (HOMOGENEOUS)
!  ================================================

!  Discrete ordinate solutions
!  ---------------------------

!  Layer loop

        DO LAYER = 1, NLAYERS

          if(lg)write(lu,'(a,i2)')'starting HOMOG_SOLUTION ',LAYER

!  Eigensolver for homogeneous solutions

          CALL HOMOG_SOLUTION &
              ( LAYER, FOURIER, NSTREAMS, NMOMENTS, &
                LOCAL_OMEGAMOMS, LOCAL_DELTAUS, &
                QUAD_STREAMS, QUAD_WEIGHTS, &
                PLM_PXI,     PLM_MXI, &
                XPOS, XNEG, KEIGEN, KTRANS, &
                EIGENMAT_SAVE, EIGENVEC_SAVE, &
                DIFVEC_SAVE, DAB_SAVE, SAB_SAVE, &
                FAIL, MESSAGE_SUB )

!          write(*,'(i4,1p4e20.12)')layer,(keigen(aa,layer),aa=1,4)
!          if (layer.eq.18)pause

!  error handling

          IF ( FAIL ) THEN
            WRITE(C2,'(I2)')LAYER
            MESSAGE = &
                'Fourier master, call to HOMOG_SOLUTION, layer'//C2
            WRITE(C3, '(I3)' ) CPT
            POINT_TRACE = 'Homog. calculation: point number: '//C3
            RETURN
          ENDIF

!  Norms of the solution vectors

          DO AA = 1, NSTREAMS
            NORM = ZERO
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              T1   = XPOS(I,AA,LAYER)  * XPOS(I,AA,LAYER)
              T2   = XPOS(I1,AA,LAYER) * XPOS(I1,AA,LAYER)
              NORM = NORM + QUAD_STRMWGT(I) * ( T1 - T2 )
            ENDDO
            EIGENNORM_SAVE(AA,LAYER) = NORM
          ENDDO

!  end layer loop

        ENDDO

!  HOmogeneous solution user-stream solutions and multipliers
!  ----------------------------------------------------------

!  Start layer loop

        DO N = 1, NLAYERS

!  homogeneous solutions at user-defined streams

          CALL UHOMOG_SOLUTION &
                ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS, &
                  LOCAL_OMEGAMOMS, XPOS, XNEG, &
                  PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI, &
                  U_XPOS, U_XNEG, U_HELP_P(1,0,N), U_HELP_M(1,0,N) )

!  Whole layer multipliers
!   @@@ RobFix 5/5/11. Small numbers analysis: added LRRS_FGSMALL,
!                      DELTAU_VERT_INPUT(N,CPT) to Argument list.

          CALL HOMOGMULT_1 &
            ( NSTREAMS, N, N_USER_STREAMS, LRRS_FGSMALL,       &
              LOCAL_TRANS_USERM, USER_STREAMS, KEIGEN, KTRANS, &
              DELTAU_VERT_INPUT(N,CPT),                        &
              HMULT_1, HMULT_2, ZETA_P, ZETA_M )

!  End layer loop

        ENDDO

!  Boundary Value Matrix
!  =====================

!  Additional setups for the albedo layer

        if(lg)write(lu,'(a)')'starting SURFACE_HOMOG_SETUP'
        CALL SURFACE_HOMOG_SETUP &
         ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, &
           SURFACE_FACTOR, FL1(CPT), FL2(CPT), BIREFLEC, &
           FOURIER, NLAYERS, NSTREAMS, &
           QUAD_STRMWGT, XPOS, XNEG, &
           R2_HOMP, R2_HOMM )

!  initialize compression matrix

        if(lg)write(lu,'(a)')'starting BVPMATRIX_INIT'
        CALL BVPMATRIX_INIT &
           ( NSTREAMS, NLAYERS, &
             NTOTAL, NSTR2, N_SUBDIAG, N_SUPDIAG, &
             BMAT_ROWMASK, BANDMAT2 )

!  set up boundary values matrix in compressed form (the "A" as in AX=B)

        if(lg)write(lu,'(a)')'starting BVPMATRIX_SETUP'
        CALL BVPMATRIX_SETUP &
           (  DO_INCLUDE_SURFACE, &
              NSTREAMS, NLAYERS,  NSTR2, BMAT_ROWMASK, &
              XPOS, XNEG, R2_HOMP, R2_HOMM, KTRANS, &
              BANDMAT2 )

!  LAPACK LU-decomposition for band matrix
!    Only for general case, otherwise this section not required

        IF ( NTOTAL .GT. 2 ) THEN
          CALL DGBTRF &
           ( NTOTAL, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
             BANDMAT2, MAX_BANDTOTAL, IPIVOT, INFO )
          FAIL = ( INFO .NE. 0 )
          IF ( INFO .GT. 0 ) THEN
            WRITE(C3, '(I3)' ) INFO
            MESSAGE_SUB = 'Singular matrix, u(i,i)=0, for i = '//C3
          ELSE IF ( INFO .LT. 0 ) THEN
            WRITE(C3, '(I3)' ) INFO
            MESSAGE_SUB = 'argument i illegal value, for i = '//C3
          ENDIF
          IF ( FAIL ) THEN
            MESSAGE = 'RAMAN_FOURIER_1, Error from DGBTRF call'
            WRITE(C3, '(I3)' ) CPT
            POINT_TRACE = 'point number: '//C3
            RETURN
          ENDIF
        ENDIF

!  Find particular solutions for elastic and all RRS terms
!  =======================================================

!  Get source vectors for wavelength point APT

        if ( lg ) write(lu, '(a)' )'starting SOURCES_MASTER_1'
        CALL SOURCES_MASTER_1 &
       ( FOURIER, APT, DO_ENERGY_BALANCING, &
         DO_RRS_OVERALL, DO_DIRECTRRS_ONLY, DO_MSMODE_LRRS, &
         DO_BIN_REALIZATION, DO_MONO_REALIZATION, DO_SSCORR_OUTGOING, &
         DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, &
         N_RRSBINS, BINMAP, NPOINTS_MONO, OFFSET_INNER, W_EXCIT, &
         NLAYERS, NSTREAMS, NMOMENTS, N_USER_STREAMS, &
         QUAD_WEIGHTS, USER_STREAMS, LRRS_FGSMALL, &
         STERM_MASK_UP, STERM_MASK_DN, FLUXES_RANKED, &
         DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, OMEGAMOMS_CABANNES, &
         OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN, OMEGAMOMS_RRSGAIN, &
         BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, &
         BEAM_DTRANS, SAVE_TRANS_USERM, &
         PLM_PXI, PLM_MXI, PLM_PXUI, PLM_MXUI, &
         PLM_00_PXI, PLM_00_MXI, PLM_00_PXUI, PLM_00_MXUI, &
         PLM_WT_PXI,  PLM_WT_MXI, &
         L0_KEIGEN,     L0_KTRANS,     L0_WPARTIC, &
         L0_WHOM_XPOS,  L0_WHOM_LCON,  L0_WHOM_MCON, &
         XPOS, KEIGEN, KTRANS, U_XPOS, U_XNEG, &
         EIGENNORM_SAVE, HMULT_1, HMULT_2, &
         RAMAN_WUPPER, RAMAN_WLOWER, WLAYER_PIST_UP, WLAYER_PIST_DN, &
         FAIL, MESSAGE_SUB )
      !    write(*,*)RAMAN_WLOWER(:,1)
      !    pause
!        pause'after sources'

!  exception handling

        IF ( FAIL ) THEN
          MESSAGE = 'Error from SOURCES_MASTER_1'
          IF ( DO_BIN_REALIZATION ) THEN
            WRITE(C3, '(I3)' ) CPT
            POINT_TRACE = 'point number: '//C3
          ENDIF
          RETURN
        ENDIF

!  Debug code
!      do n = 1, 23
!         write(67,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i,n),i=1,8)
!         write(67,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i+8,n),i=1,8)
!         write(67,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i,n),i=1,8)
!         write(67,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i+8,n),i=1,8)
!      enddo
!      do n = 1, 23
!         write(*,'(i4,1p8e19.10)') n,WLAYER_PIST_UP(n,1), WLAYER_PIST_DN(n,1)
!      enddo
!     if ( apt.eq.1)pause'end sources'
!      do n = 1, 25
!         write(67,'(i4,1p8e19.10)')n,(WLOWER(i,n),i=1,8)
!         write(67,'(i4,1p8e19.10)')n,(WLOWER(i+8,n),i=1,8)
!         write(67,'(i4,1p8e19.10)')n,(WUPPER(i,n),i=1,8)
!         write(67,'(i4,1p8e19.10)')n,(WUPPER(i+8,n),i=1,8)
!      enddo
!      do n = 1, 25
!         write(67,'(i4,1p8e19.10)')
!     *           n,WLAYER_PIST_UP(n,1), WLAYER_PIST_DN(n,1)
!      enddo
!        if ( apt.eq.1)pause'end sources'

!        write(*,*)'regg',apt,wupper(2,24),wlower(10,25)
!        if ( apt.eq.1) pause

!  ****************
!  Radiation Field:
!  ****************

!  Complete and Solve boundary value problem
!  =========================================

!  Add contributions to the BV problem RHS

        if (lg)write(lu ,'(a)')' **starting SURFACE_BEAM_SETUP'
        CALL SURFACE_BEAM_SETUP &
         ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, &
           SURFACE_FACTOR, FL1(CPT), FL2(CPT), BIREFLEC, &
           FOURIER, NLAYERS, NSTREAMS, &
           QUAD_STRMWGT, RAMAN_WLOWER, &
           R2_BEAM )

!  set up Column COL2 for solution vector (the "B" as in AX=B)

        if (lg)write(lu ,'(a)')' **starting BVPCOLUMN_SETUP'
        CALL BVPCOLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, &
             NSTREAMS, NLAYERS, NTOTAL, NSTR2, &
             DIRECT_BEAM, RAMAN_WUPPER, RAMAN_WLOWER, R2_BEAM, &
             COL2 )

!  Slab problem control
!        IF (NTOTAL.EQ.2) GO TO 750

!  LAPACK substitution using RHS column vector COL2
!     special case, 1 layer, 1 stream

        IF (NTOTAL.GT. 2 ) THEN
          CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAX_BANDTOTAL, IPIVOT, &
              COL2, MAX_TOTAL, INFO )
          IF ( INFO .LT. 0 ) THEN
            FAIL = .TRUE.
            WRITE(C3, '(I3)' ) INFO
            MESSAGE_SUB = 'argument i illegal value, for i = '//C3
            MESSAGE     = 'RAMAN_FOURIER_1, DGBTRS call'
            IF ( DO_BIN_REALIZATION ) THEN
              WRITE(C3, '(I3)' ) CPT
              POINT_TRACE = 'point number: '//C3
            ENDIF
            RETURN
          ENDIF
        ELSE
          DETER = BANDMAT2(1,1) * BANDMAT2(2,2) - &
                  BANDMAT2(1,2) * BANDMAT2(2,1)
          V1 = COL2(1,1)
          V2 = COL2(2,1)
          COL2(1,1) = ( BANDMAT2(2,2)*V1 - BANDMAT2(1,2)*V2 ) / DETER
          COL2(2,1) = ( BANDMAT2(1,1)*V2 - BANDMAT2(2,1)*V1 ) / DETER
        ENDIF

!  Set integration constants
!  =========================

!  Set integration constants LCON and MCON for -/+ eigensolutions, all l

        DO N = 1, NLAYERS
          C0 = (N-1)*NSTR2
          DO I = 1, NSTREAMS
            I1 = I+NSTREAMS
            LCON(I,N) = COL2(C0+I,1)
            MCON(I,N) = COL2(C0+I1,1)
          ENDDO
        ENDDO

!  Debug code, 3 June 2009

!        DO N = 1, NLAYERS
!          DO I = 1, NSTREAMS
!            write(62,5)n,i,lcon(i,n),  mcon(i,n)
!          ENDDO
!        ENDDO
! 5    format(2i5,1p2e25.16)
!        pause'after BVP'

!  Quadrature solutions: multiplied by integration constants

        DO N = 1, NLAYERS
          DO I = 1, NSTR2
            DO AA = 1, NSTREAMS
              LCON_XVEC(I,AA,N) = LCON(AA,N) * XPOS(I,AA,N)
              MCON_XVEC(I,AA,N) = MCON(AA,N) * XNEG(I,AA,N)
            ENDDO
          ENDDO
        ENDDO

!  ***************
!  Post-Processing
!  ***************

!  Get the homgeneous solution layer source terms
!  ==============================================

        DO N = 1, NLAYERS
         IF ( STERM_MASK_UP(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            HOM = ZERO
            DO AA = 1, NSTREAMS
              HOM = HOM + LCON(AA,N)*U_XPOS(UM,AA,N) * HMULT_2(AA,UM,N) &
                        + MCON(AA,N)*U_XNEG(UM,AA,N) * HMULT_1(AA,UM,N)
            ENDDO
            WLAYER_HOST_UP(N,UM) = HOM
          ENDDO
         ENDIF
         IF ( STERM_MASK_DN(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            HOM = ZERO
            DO AA = 1, NSTREAMS
              HOM = HOM + LCON(AA,N)*U_XNEG(UM,AA,N) * HMULT_1(AA,UM,N) &
                        + MCON(AA,N)*U_XPOS(UM,AA,N) * HMULT_2(AA,UM,N)
            ENDDO
            WLAYER_HOST_DN(N,UM) = HOM
          ENDDO
         ENDIF
        ENDDO

!  Get the TOA/BOA source terms
!  ============================

!  Upwelling................................

        IF ( DO_UPWELLING ) THEN

!  User direct beam

          CALL UDBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_SSCORR_OUTGOING, &
             DO_LAMBERTIAN_SURFACE, FL1(CPT), FL2(CPT), USER_BIREFLEC_0, &
             FOURIER, N_USER_STREAMS, ATMOS_ATTN, &
             USER_DIRECT_BEAM )

!  Source term

          CALL BOASOURCE &
          ( DO_INCLUDE_SURFACE, DO_REFLECTED_DIRECTBEAM, &
            DO_LAMBERTIAN_SURFACE,  FOURIER, &
            SURFACE_FACTOR, FL1(CPT), FL2(CPT), &
            NSTREAMS, NLAYERS, N_USER_STREAMS, &
            RAMAN_WLOWER, LCON_XVEC, MCON_XVEC, KTRANS, &
            QUAD_STRMWGT, USER_BIREFLEC, USER_DIRECT_BEAM, &
            IDOWNSURF, BOA_SOURCE, DIRECT_BOA_SOURCE )

        ENDIF

!  Downwelling (no source)..................

        IF ( DO_DNWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            TOA_SOURCE(UM) = ZERO
          ENDDO
        ENDIF

!  quadrature solutions debug----------------------
!   Off-grid solutions check out. 7 December 2006.

!      if ( do_upwelling ) then
!       if ( .not. do_rrs_overall ) then
!        write(39,*)' TOA quadup elastic, Fourier ',FOURIER,APT
!        DO I = 1, NSTREAMS
!         I1 = I+NSTREAMS
!         SUM = WUPPER(I1,1)
!         DO AA = 1, NSTREAMS
!          SUM = SUM + LCON_XVEC(I1,AA,1) +
!     &                  MCON_XVEC(I1,AA,1)*KTRANS(AA,1)
!         ENDDO
!         write(39,*)i,XPOS(I1,1,1),sum/pi4
!        enddo
!       else
!        write(39,*)' TOA quadup recon elastic, Fourier ',FOURIER,APT,CPT
!        DO I = 1, NSTREAMS
!         I1 = I+NSTREAMS
!         SUM = L0_WPARTIC(I1,1,CPT)
!         DO AA = 1, NSTREAMS
!          SUM = SUM + L0_WHOMOG_DN(I1,AA,1,CPT) +
!     &                L0_WHOMOG_UP(I1,AA,1,CPT)*L0_KTRANS(AA,1,CPT)
!         ENDDO
!         write(39,*)i,L0_KTRANS(i,1,CPT),sum
!        enddo
!        write(39,*)' TOA quadup RRS, Fourier ',FOURIER,APT
!        DO I = 1, NSTREAMS
!         I1 = I+NSTREAMS
!         SUM  = WLOWER(I,11)
!         SUM1 = WUTAUS(I,2)
!         DO AA = 1, NSTREAMS
!          SUM = SUM   + LCON_XVEC(I,AA,11)*KTRANS(AA,11) +
!     &                  MCON_XVEC(I,AA,11)
!          SUM1 = SUM1 + LCON_XVEC(I,AA,11)*UTDN_KTRANS(AA,2) +
!     &                  MCON_XVEC(I,AA,11)*UTUP_KTRANS(AA,2)
!         ENDDO
!         write(39,*)i,sum1/pi4,sum/pi4
!        enddo
!       endif
!      endif
!
!      if ( do_dnwelling ) then
!      write(*,*)'dnwelling quadrature BOA layer bottom, Fourier ',FOURIER
!      DO I = 1, NSTREAMS
!        I1 = I+NSTREAMS
!        SUM = WLOWER(I,NLAYERS)
!        DO AA = 1, NSTREAMS
!          SUM = SUM + LCON_XVEC(I,AA,NLAYERS) * KTRANS(AA,NLAYERS)
!     &                + MCON_XVEC(I,AA,NLAYERS)
!        ENDDO
!        write(*,*)i,dacos(quad_streams(i))/deg_to_rad,sum/pi4
!      enddo
!      endif

!  RAMAN output and MIFLUX
!  =======================

!  Upwelling Intensity

        IF ( DO_UPWELLING ) THEN
          if (lg)write(lu ,'(a)')' **starting RAMAN_INTENSITY_UP_1'
          CALL RAMAN_INTENSITY_UP_1 &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER, &
            FOURIER, NLAYERS, NSTREAMS, N_LOUTPUT, N_USER_STREAMS, &
            LEVELMASK_UP, RAMAN_WUPPER, RAMAN_WLOWER, KTRANS, &
            LCON_XVEC, MCON_XVEC, LOCAL_TRANS_USERM, &
            BOA_SOURCE, DIRECT_BOA_SOURCE, &
            WLAYER_HOST_UP, WLAYER_PIST_UP, &
            CUMSOURCE_UP, RAMAN_F_UP(1,1,APT), QUADRAMAN_F_UP(1,1,APT) )

            ! write(*,*)1,RAMAN_F_UP(1,1:3,APT)
            ! write(*,*)2,RAMAN_F_UP(2,1:3,APT)
            ! write(*,*)3,RAMAN_F_UP(3,1:3,APT)
      !      pause!'Fourier_Master, Fourier zero, APT #1'

       
!  Saved results for the Rayleigh planetary problem

!         IF ( DO_RAYLEIGH_ONLY ) THEN
!          DO J = 1, N_LOUTPUT
!           DO K = 1, N_USER_STREAMS
!            IF(M.EQ.0)RAMAN_F0_UP(J,K,APT)=RAMAN_F_UP(J,K,APT)
!            IF(M.EQ.1)RAMAN_F1_UP(J,K,APT)=RAMAN_F_UP(J,K,APT)
!            IF(M.EQ.2)RAMAN_F2_UP(J,K,APT)=RAMAN_F_UP(J,K,APT)
!           ENDDO
!          ENDDO
!         ENDIF

!  End upwelling

        ENDIF

!  debug

!      if ( do_upwelling ) then
!       write(40,*)' TOA user, Fourier ',FOURIER, DO_RRS_OVERALL,APT
!       DO I = 1, N_USER_STREAMS
!        write(*,*)i,FOURIER, &
!           dacos(USER_STREAMS(i))/deg_to_rad,RAMAN_F_UP(1,i,APT)
!       enddo
!       if ( DO_RRS_OVERALL ) then
!        write(40,*)' TOA recon user, Fourier ',FOURIER, APT
!        DO I = 1, N_USER_STREAMS
!         write(40,*)i,FOURIER,
!     &      dacos(USER_STREAMS(i))/deg_to_rad,ELASTIC_F_UP(1,i,APT)
!        enddo
!       endif
!      endif

!  Downwelling Intensity

        IF ( DO_DNWELLING ) THEN
         if (lg)write(lu ,'(a)')' **starting RAMAN_INTENSITY_DN_1'
         CALL RAMAN_INTENSITY_DN_1 &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER, &
            FOURIER, NSTREAMS, N_LOUTPUT, N_USER_STREAMS, &
            LEVELMASK_DN, RAMAN_WLOWER, KTRANS, &
            LCON_XVEC, MCON_XVEC, LOCAL_TRANS_USERM, &
            TOA_SOURCE, WLAYER_HOST_DN, WLAYER_PIST_DN, &
            CUMSOURCE_DN, RAMAN_F_DN(1,1,APT), QUADRAMAN_F_DN(1,1,APT) )

!  Saved results for the Rayleigh planetary problem

!         IF ( DO_RAYLEIGH_ONLY ) THEN
!          DO J = 1, N_LOUTPUT
!           DO K = 1, N_USER_STREAMS
!            IF(M.EQ.0)RAMAN_F0_DN(J,K,APT)=RAMAN_F_DN(J,K,APT)
!            IF(M.EQ.1)RAMAN_F1_DN(J,K,APT)=RAMAN_F_DN(J,K,APT)
!            IF(M.EQ.2)RAMAN_F2_DN(J,K,APT)=RAMAN_F_DN(J,K,APT)
!           ENDDO
!          ENDDO
!         ENDIF

!  end  downwelling

        ENDIF

!  debug

! 5       format(2i5,1p2e25.16)
!          write(*,5)1,APT,RAMAN_F_UP(1,1,APT),RAMAN_F_DN(2,1,APT)

!      if ( do_dnwelling ) then
!      write(*,*)'dnwelling off-quad BOA layer bottom, Fourier ',FOURIER
!      DO I = 1, N_USER_STREAMS
!        write(*,*)i,
!     &    dacos(USER_STREAMS(i))/deg_to_rad,RAMAN_F(NLAYERS+1,i,2)
!      enddo
!      endif

!  mean value output
!  -----------------

!  Check BEAM_ETRANS(1,CPT), formerly BEAM_ETRANS(1,APT)

        IF ( DO_INCLUDE_MVOUT ) THEN
          if (lg)write(lu ,'(a)')' **starting MIFLUX_INTENSITY_1'
          CALL MIFLUX_INTENSITY_1 &
           ( DO_UPWELLING, DO_DNWELLING, &
             COS_SZA,  FLUX_FACTOR, &
             NSTREAMS, N_LOUTPUT, LEVELMASK_DN, &
             BEAM_PICUTOFF(CPT), BEAM_ETRANS(1,CPT), &
             QUAD_WEIGHTS, QUAD_STRMWGT, &
             QUADRAMAN_F_UP(1,1,APT), QUADRAMAN_F_DN(1,1,APT), &
             MEAN_RAMAN_UP(1,APT), FLUX_RAMAN_UP(1,APT), &
             MEAN_RAMAN_DN(1,APT), FLUX_RAMAN_DN(1,APT) )
        ENDIF

!  Debug, 13 March 2008
!            write(*,*)'intensity up TOA',
!     &      apt, MEAN_RAMAN_UP(1,APT),FLUX_RAMAN_UP(1,APT)
!            write(*,*)'intensity up BOA',
!     &      apt, MEAN_RAMAN_UP(2,APT),FLUX_RAMAN_UP(2,APT)
!            write(*,*)'intensity dn TOA',
!     &      apt, MEAN_RAMAN_DN(1,APT),FLUX_RAMAN_DN(1,APT)
!            write(*,*)'intensity dn BOA',
!     &      apt, MEAN_RAMAN_DN(2,APT),FLUX_RAMAN_DN(2,APT)
!        pause

!  Finish point loop

      ENDDO

!      pause 'End Fourier 1, reg'

!  Finish
!  ======

      RETURN
      END SUBROUTINE RAMAN_FOURIER_1


      END MODULE lrrs_fourier_master_1

