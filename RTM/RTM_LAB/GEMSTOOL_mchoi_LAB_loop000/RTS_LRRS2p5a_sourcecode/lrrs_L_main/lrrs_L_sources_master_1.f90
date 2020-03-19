
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
! #   MASTER_MODULE :    L_SOURCES_MASTER_1                     #
! #                        (1 = Whole-layer)                    #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Use of Taylor-series control, TAYLOR_ORDER index, passed to PostProcessing
!    (2) Introduction of variable number of surface Jacobians, dimensioning to match
!    (3) Bookkeeping improvements (use of "Only", clearer I/O specifications)

      MODULE lrrs_L_sources_master_1_m

!      USE LRRS_PARS_m, Only : LDU

      PRIVATE
      PUBLIC :: L_SOURCES_MASTER_1

      CONTAINS

      SUBROUTINE L_SOURCES_MASTER_1 &
       ( FOURIER, ACTUAL_POINT, DO_INCLUDE_SURFACE, DO_ENERGY_BALANCING,      & ! Inputs
         DO_RRS_OVERALL, DO_DIRECTRRS_ONLY, DO_MSMODE_LRRS,                   & ! Inputs
         DO_BIN_REALIZATION, DO_MONO_REALIZATION, DO_SSCORR_OUTGOING,         & ! Inputs
         DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                         & ! Inputs
         DO_COLUMN_WFS, DO_PROFILE_WFS, DO_SURFACE_WFS,                       & ! Inputs
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS, N_SURFACEWFS, & ! Inputs
         N_RRSBINS, BINMAP, NPOINTS_MONO, OFFSET_INNER, W_EXCIT,              & ! Inputs
         NLAYERS, NSTREAMS, NMOMENTS, N_USER_STREAMS,                         & ! Inputs
         QUAD_WEIGHTS, USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER,              & ! Inputs
         STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, FLUXES_RANKED,               & ! Inputs
         DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, OMEGAMOMS_CABANNES,            & ! Inputs
         OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN,  OMEGAMOMS_RRSGAIN,             & ! Inputs
         L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC, L_OMEGAMOMS_CABANNES,      & ! Inputs
         L_OMEGAMOMS_RRSLOSS, L_OMEGAMOMS_RRSBIN,  L_OMEGAMOMS_RRSGAIN,       & ! Inputs
         L_BEAM_ITRANS, L_BEAM_AVSECANT, L_BEAM_DTRANS,                       & ! Inputs
         BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS,              & ! Inputs
         SAVE_TRANS_USERM, L_SAVE_TRANS_USERM,                                & ! Inputs
         PLM_PXI, PLM_MXI, PLM_PXUI, PLM_MXUI,                                & ! Inputs
         PLM_00_PXI, PLM_00_MXI, PLM_00_PXUI, PLM_00_MXUI,                    & ! Inputs
         PLM_WT_PXI,  PLM_WT_MXI,                                             & ! Inputs
         L0_KEIGEN,     L0_KTRANS,     L0_WPARTIC,                            & ! Inputs
         L0_WHOM_XPOS,  L0_WHOM_LCON,  L0_WHOM_MCON,                          & ! Inputs
         L_L0_KEIGEN,      L_L0_KTRANS,     L_L0_WPARTIC,                     & ! Inputs
         L_L0_WHOM_XPOS,   L_L0_WHOM_LCON,  L_L0_WHOM_MCON,                   & ! Inputs
         LS_L0_WHOM_LCON,  LS_L0_WHOM_MCON,                                   & ! Inputs
         XPOS, KEIGEN, KTRANS, U_XPOS, U_XNEG,                                & ! Inputs
         EIGENNORM_SAVE, HMULT_1, HMULT_2,                                    & ! Inputs
         L_XPOS, L_KEIGEN, L_KTRANS, L_U_XPOS, L_U_XNEG,                      & ! Inputs
         L_EIGENNORM_SAVE, L_HMULT_1, L_HMULT_2,                              & ! Inputs
         NKSTORAGE, NKSTORAGE2,                                               & ! Inputs
         RAMAN_WUPPER,   L_RAMAN_WUPPER,   LS_RAMAN_WUPPER,    & ! Outputs
         RAMAN_WLOWER,   L_RAMAN_WLOWER,   LS_RAMAN_WLOWER,    & ! Outputs
         WLAYER_PIST_UP, L_WLAYER_PIST_UP, LS_WLAYER_PIST_UP,  & ! Outputs
         WLAYER_PIST_DN, L_WLAYER_PIST_DN, LS_WLAYER_PIST_DN,  & ! Outputs
         FAIL, MESSAGE )                                         ! Outputs

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, MAX_2_STREAMS, MAX_STREAMS, MAX_USER_STREAMS, &
                              MAX_MOMENTS, MAX_LAYERS, MAX_LAYERS_SQ, MAX_LAYERS_NK,        &
                              MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_BINS, MAX_POINTS

!  Use modules

      USE LRRS_RTSOLUTIONS_1_m     , Only : GREENFUNC_SOLUTION_1
      USE LRRS_POSTPROCESSING_1_m  , Only : POSTPROCESSING_1

      USE LRRS_L_RTSOLUTIONS_1_m   , Only : LPC_GREENFUNC_SOLUTION_1, LS_GREENFUNC_SOLUTION_1
      USE LRRS_L_POSTPROCESSING_1_m, Only : LPC_POSTPROCESSING_1    , LS_POSTPROCESSING_1

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Fourier point, wavelength point

      INTEGER  , INTENT(IN) :: FOURIER
      INTEGER  , INTENT(IN) :: ACTUAL_POINT

!  Flags
!  -----

!  Solution method (options are mutually exclusive)

      LOGICAL  , INTENT(IN) :: DO_ENERGY_BALANCING

!  Surface

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  overall filling flags (directrrs is the solar term only)

      LOGICAL  , INTENT(IN) :: DO_RRS_OVERALL
      LOGICAL  , INTENT(IN) :: DO_DIRECTRRS_ONLY
      LOGICAL  , INTENT(IN) :: DO_MSMODE_LRRS

!  Binning or monochromatic realization
!    (Options are mutually exclusive)

      LOGICAL  , INTENT(IN) :: DO_BIN_REALIZATION
      LOGICAL  , INTENT(IN) :: DO_MONO_REALIZATION

!  SS correction flag; Now incorporates the DB term

      LOGICAL  , INTENT(IN) :: DO_SSCORR_OUTGOING

!  Upwelling and downwelling flags

      LOGICAL  , INTENT(IN) :: DO_UPWELLING
      LOGICAL  , INTENT(IN) :: DO_DNWELLING

!  output for results at user angles (post processing)

      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS

!  Linearization control profile Jacobians

      LOGICAL  , INTENT(IN) :: do_profile_wfs
      LOGICAL  , INTENT(IN) :: layer_vary_flag   (max_layers)
      INTEGER  , INTENT(IN) :: layer_vary_number (max_layers)

!  Column linearization control inputs

      LOGICAL  , INTENT(IN) :: do_column_wfs
      INTEGER  , INTENT(IN) :: n_totalcolumn_wfs

!  Surface linearization

      LOGICAL  , INTENT(IN) :: do_surface_wfs
      INTEGER  , INTENT(IN) :: n_surfacewfs

!  Points
!  ------

!  offset for the inner range
!    Depends on the choice of solar spectrum

      INTEGER  , INTENT(IN) :: OFFSET_INNER

!  Bin Mapping quantities

      INTEGER  , INTENT(IN) :: N_RRSBINS  ( MAX_POINTS )
      INTEGER  , INTENT(IN) :: BINMAP ( MAX_BINS, MAX_POINTS )

!  Number of monochromatic points to be calculated = 234

      INTEGER  , INTENT(IN) :: NPOINTS_MONO

!  W_EXCIT is the index of the excitation wavelength in the Rank.
!     ( Monochromatic only )

      INTEGER  , INTENT(IN) :: W_EXCIT

!  Numbers and bookkeeping
!  -----------------------

!  Number of layers

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of discrete ordinate streams

      INTEGER  , INTENT(IN) :: NSTREAMS

!  Number of moments (elastic scattering)

      INTEGER  , INTENT(IN) :: NMOMENTS

!  User stream variables

      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  Taylor-series control. Updated, Version 2.5, 9/10/15

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  Layer masks for doing integrated source terms

      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP ( MAX_LAYERS )
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN ( MAX_LAYERS )

!  Quadrature weights

      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS  ( MAX_STREAMS )

!  Optical properties
!  ------------------

!  Depends on the choice of solar spectrum

      REAL(FPK), INTENT(IN) :: FLUXES_RANKED   ( MAX_POINTS )

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT   ( MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived Cabannes scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_CABANNES &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived rotational raman scattering (Internal variable). Loss term
!    Single scatter albedo * phase function moments
!    Only required for the Energy-balance approximation.

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSLOSS &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived rotational Raman scattering (binning)

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSBIN &
               ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Derived rotational Raman scattering (Monochromatic). Gain terms
!    Single scatter albedo * phase function moments

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSGAIN &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Linearized Optical Properties
!  -----------------------------

!  Vertical optical thickness values

      REAL(FPK), INTENT(IN) :: &
          L_DELTAU_VERT_INPUT   ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_ELASTIC &
           ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived Cabannes scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_CABANNES &
           ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived rotational raman scattering (Internal variable). Loss term
!    Single scatter albedo * phase function moments
!    Only required for the Energy-balance approximation.

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSLOSS &
               ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived rotational Raman scattering (Internal variable). Binning only

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSBIN &
               ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Derived rotational Raman scattering (Internal variable). Gain terms
!    Single scatter albedo * phase function moments. Monochromatic only.

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSGAIN &
               ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Local setup variables
!  ---------------------

!  Solar beam transmittances, average secant factors

      INTEGER  , INTENT(IN) :: BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(IN) :: SAVE_TRANS_USERM &
            (  MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Solar beam transmittances, average secant factors, LINEARIZED

      REAL(FPK), INTENT(IN) :: L_BEAM_ITRANS &
               ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_BEAM_AVSECANT &
               ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_BEAM_DTRANS &
               ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )

!  user stream transmittances, whole layers, LINEARIZED

      REAL(FPK), INTENT(IN) :: L_SAVE_TRANS_USERM &
            ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Legendre polynomial variables
!  -----------------------------

!  Legendre functions at discrete ordinates

      REAL(FPK), INTENT(IN) :: PLM_PXI    ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_MXI    ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Product functions, discrete ordinates and solar stream

      REAL(FPK), INTENT(IN) :: PLM_00_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_00_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Product of Legendre functions with quadrature weights

      REAL(FPK), INTENT(IN) :: PLM_WT_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_WT_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Legendre functions at user defined polar streams

      REAL(FPK), INTENT(IN) :: PLM_PXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_MXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  Product functions, User streams and solar stream

      REAL(FPK), INTENT(IN) :: PLM_00_PXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_00_MXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  local homogeneous solutions
!  ---------------------------

!  (Positive) Eigenvalues

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )

!   Whole layer transmittance factors for +/- eigenvalues

      REAL(FPK), INTENT(IN) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  user-stream homogeneous solution vectors

      REAL(FPK), INTENT(IN) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Green's function Eigennorms

      REAL(FPK), INTENT(IN) :: EIGENNORM_SAVE ( MAX_STREAMS, MAX_LAYERS )

!  Whole layer multipliers (homogeneous)

      REAL(FPK), INTENT(IN) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  local homogeneous solutions, Linearizations
!  -------------------------------------------

!  (Positive) Eigenvalues

      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!   Whole layer transmittance factors for +/- eigenvalues

      REAL(FPK), INTENT(IN) :: L_KEIGEN ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: L_XPOS &
           ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  user-stream homogeneous solution vectors

      REAL(FPK), INTENT(IN) :: L_U_XPOS &
            ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_U_XNEG &
            ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Whole layer multipliers (homogeneous)

      REAL(FPK), INTENT(IN) :: L_HMULT_1 &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_HMULT_2 &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Eigennorm saved

      REAL(FPK), INTENT(IN) :: L_EIGENNORM_SAVE &
                   ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  Elastic discrete ordinate solutions at shifted wavelengths
!  ----------------------------------------------------------

!  Eigenvalues

      REAL(FPK), INTENT(IN) :: L0_KEIGEN &
              ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Eigentransmittances

      REAL(FPK), INTENT(IN) :: L0_KTRANS &
              ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Particular solution vectors

      REAL(FPK), INTENT(IN) :: L0_WPARTIC &
            ( MAX_2_STREAMS, MAX_LAYERS, MAX_POINTS )

!  homogeneous solution vectors and integration constants

      REAL(FPK), INTENT(IN) :: L0_WHOM_XPOS &
            ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L0_WHOM_LCON &
            ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L0_WHOM_MCON &
            ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Linearized Elastic discrete ordinate solutions at shifted wavelengths
!  ---------------------------------------------------------------------

!  Eigenvalues

      REAL(FPK), INTENT(IN) :: L_L0_KEIGEN &
              ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Eigentransmittances

      REAL(FPK), INTENT(IN) :: L_L0_KTRANS &
              ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Particular solution vectors

      REAL(FPK), INTENT(IN) :: L_L0_WPARTIC &
            ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS_NK, MAX_POINTS )

!  homogeneous solution vectors and integration constants

      REAL(FPK), INTENT(IN) :: L_L0_WHOM_XPOS &
            ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, &
                                     MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_L0_WHOM_LCON &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS_SQ, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_L0_WHOM_MCON &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS_SQ, MAX_POINTS )

!  Surface linearization of integration constants

      REAL(FPK), INTENT(IN) :: LS_L0_WHOM_LCON &
            ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: LS_L0_WHOM_MCON &
            ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Linearization bookkeeping storage

      INTEGER  , INTENT(IN) :: NKSTORAGE  ( MAX_LAYERS, 0:MAX_LAYERS )
      INTEGER  , INTENT(IN) :: NKSTORAGE2 ( MAX_LAYERS, 0:MAX_LAYERS )

!  OUTPUT ARGUMENTS
!  ================

!  output arguments from Beam solution modules
!  -------------------------------------------

!  Particular solutions at upper and lower boundaries
!   (as evaluated for the inelastic field by Green's function method)

      REAL(FPK), INTENT(OUT) :: RAMAN_WUPPER  ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(OUT) :: RAMAN_WLOWER  ( MAX_2_STREAMS, MAX_LAYERS )

!  post processed up and down Raman solutions

      REAL(FPK), INTENT(OUT) :: WLAYER_PIST_UP ( MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: WLAYER_PIST_DN ( MAX_LAYERS, MAX_USER_STREAMS )

!  Particular solutions at upper and lower boundaries
!   (as evaluated for the inelastic field by Green's function method)

      REAL(FPK), INTENT(OUT) :: L_RAMAN_WUPPER &
            ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS_SQ )
      REAL(FPK), INTENT(OUT) :: L_RAMAN_WLOWER &
            ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS_SQ )

!  Particular solutions at upper and lower boundaries
!   (as evaluated for the inelastic field by Green's function method)

      REAL(FPK), INTENT(OUT) :: LS_RAMAN_WUPPER (MAX_SURFACEWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(OUT) :: LS_RAMAN_WLOWER (MAX_SURFACEWFS, MAX_2_STREAMS, MAX_LAYERS )

!  post processed up and down Raman solutions

      REAL(FPK), INTENT(OUT) :: L_WLAYER_PIST_UP &
             ( MAX_ATMOSWFS, MAX_LAYERS_SQ, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: L_WLAYER_PIST_DN &
             ( MAX_ATMOSWFS, MAX_LAYERS_SQ, MAX_USER_STREAMS )

      REAL(FPK), INTENT(OUT) :: LS_WLAYER_PIST_UP &
             ( MAX_SURFACEWFS, MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK), INTENT(OUT) :: LS_WLAYER_PIST_DN &
             ( MAX_SURFACEWFS, MAX_LAYERS, MAX_USER_STREAMS )

!  module exception handling
!  -------------------------

      LOGICAL            , INTENT(OUT) :: FAIL
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE

!  Local variables
!  ===============

!  Source function control

      LOGICAL   :: CHANGE_EXPONENT
      LOGICAL   :: FLIPPER
      INTEGER   :: PICUTOFF
      LOGICAL   :: DO_MULTIPLIER

!  Source function transmittances, average secant factors

      REAL(FPK) :: ASOURCE
      REAL(FPK) :: INITRANS
      REAL(FPK) :: DELTRANS

!  LInearized Source function transmittances, average secant factors

      REAL(FPK) :: L_ASOURCE  ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK) :: L_INITRANS ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK) :: L_DELTRANS ( MAX_ATMOSWFS, 0:MAX_LAYERS )

!  Source function Q-vectors

      REAL(FPK) :: QSOURCE_QPOS  ( MAX_STREAMS )
      REAL(FPK) :: QSOURCE_QNEG  ( MAX_STREAMS )
      REAL(FPK) :: QSOURCE_UPOS1 ( MAX_USER_STREAMS )
      REAL(FPK) :: QSOURCE_UNEG1 ( MAX_USER_STREAMS )

!  Saved local multipliers

      REAL(FPK) :: WSU_SUSAV ( MAX_STREAMS, MAX_USER_STREAMS )
      REAL(FPK) :: WSU_SDSAV ( MAX_STREAMS, MAX_USER_STREAMS )
      REAL(FPK) :: WSD_SUSAV ( MAX_STREAMS, MAX_USER_STREAMS )
      REAL(FPK) :: WSD_SDSAV ( MAX_STREAMS, MAX_USER_STREAMS )

!  Green's function multipliers

      REAL(FPK) :: MULT_M  ( MAX_STREAMS )
      REAL(FPK) :: MULT_P  ( MAX_STREAMS )

!  Source function Beam multipliers, whole layer

      REAL(FPK) :: EMULT_UP ( MAX_USER_STREAMS )
      REAL(FPK) :: EMULT_DN ( MAX_USER_STREAMS )

!  Linearized (Atmospheric) Source function Q-vectors

      REAL(FPK) :: &
        L_QSOURCE_QPOS  ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS ), &
        L_QSOURCE_QNEG  ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS ), &
        L_QSOURCE_UPOS1 ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS ), &
        L_QSOURCE_UNEG1 ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS )

!  Linearized (Surface) function Q-vectors

      REAL(FPK) ::     LS_QSOURCE_QPOS  ( MAX_SURFACEWFS, MAX_STREAMS ), &
                       LS_QSOURCE_QNEG  ( MAX_SURFACEWFS, MAX_STREAMS ), &
                       LS_QSOURCE_UPOS1 ( MAX_SURFACEWFS, MAX_USER_STREAMS ), &
                       LS_QSOURCE_UNEG1 ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Green's function arrays

      REAL(FPK) :: ATERM_SAVE ( MAX_STREAMS )
      REAL(FPK) :: BTERM_SAVE ( MAX_STREAMS )
      REAL(FPK) :: L_ATERM_SAVE ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS )
      REAL(FPK) :: L_BTERM_SAVE ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS )
      REAL(FPK) :: LS_ATERM_SAVE ( MAX_SURFACEWFS, MAX_STREAMS )
      REAL(FPK) :: LS_BTERM_SAVE ( MAX_SURFACEWFS, MAX_STREAMS )

!  Source function Beam multipliers, whole layer
!    Now used locally.......................................
!      REAL(FPK) :: L_EMULT_UP
!     &        ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS )
!      REAL(FPK) :: L_EMULT_DN
!     &        ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS )

!  Local optical depths

      REAL(FPK) :: DELTAUS
      REAL(FPK) :: L_DELTAUS ( MAX_ATMOSWFS )

!  Trans UserM

      REAL(FPK) :: TRANS_USERM   ( MAX_USER_STREAMS )
      REAL(FPK) :: L_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!  Local arrays of omega-moments (IOPs)

      REAL(FPK) :: OMEGA_PHASMOMS     ( MAX_LAYERS, 0: MAX_MOMENTS )
      REAL(FPK) :: OMEGAMOMS_RRSLOCAL ( 0: MAX_MOMENTS )
      REAL(FPK) :: L_OMEGA_PHASMOMS &
                  ( MAX_ATMOSWFS, MAX_LAYERS, 0: MAX_MOMENTS )
      REAL(FPK) :: L_OMEGAMOMS_RRSLOCAL ( MAX_ATMOSWFS, 0: MAX_MOMENTS )

!  Local arrays of optical depth scaling factors  + linearizations

      REAL(FPK) :: DELTA_SCALING
      REAL(FPK) :: L_DELTA_SCALING  ( MAX_ATMOSWFS )

!  Elastic only flag

      LOGICAL   :: DO_ELASTIC_ONLY

!  Flux comparisons

      REAL(FPK) :: FM_SHIFTED
      REAL(FPK) :: FM_SHIFTED_DIRECT
      REAL(FPK) :: FM_SHIFTED_DIFFUSE
      REAL(FPK) :: FM_NONSHIFT
      REAL(FPK) :: FM_NONSHIFT_DIRECT
      REAL(FPK) :: FM_NONSHIFT_DIFFUSE

!  Dummy initial transmittances

      REAL(FPK) :: DUMTRANS, L_DUMTRANS(MAX_ATMOSWFS,0:MAX_LAYERS)

!  Linearization control numbers
!    NK storage arrays are in LRRS_BOOKKEEP.VARS

      INTEGER   :: KS(MAX_LAYERS),  KF(MAX_LAYERS)
      INTEGER   :: KS2(MAX_LAYERS), KF2(MAX_LAYERS)
      INTEGER   :: KPARS(0:MAX_LAYERS)

!  Derived Jacobian flags

      LOGICAL   :: DO_ATMOS_JACOBIANS
      LOGICAL   :: DO_SURFACE_JACOBIANS

!  other local variables

      LOGICAL   :: DO_SSCORR_OVERALL
      INTEGER   :: S, P, Z, ZS, N, I, UM, ID, ID1, SAVCUT, WW
      INTEGER   :: M, J, UI, L, APT, BPT, CPT, Q, K, NK, NC, NC2
      INTEGER   :: NRRS_MOMENTS, NPOINTS_LOCAL, NK2

      REAL(FPK) :: L_HOLD_RRS(MAX_ATMOSWFS,0:MAX_LAYERS,0:2)
      REAL(FPK) :: HOLD_RRS(0:2), LS_HOLD_RRS(MAX_SURFACEWFS,0:2)
      REAL(FPK) :: HOLD_MP(0:2),  L_HOLD_MP,  LS_HOLD_MP
      REAL(FPK) :: MC, LC, L_MC, L_LC,  LS_MC, LS_LC, LP, LM
      REAL(FPK) :: PHIFUNC, SP, SM

!  debug variables

      LOGICAL  , PARAMETER :: DO_LOG_OUTPUT    = .FALSE.
      LOGICAL  , PARAMETER :: DO_DEBUG_SOURCES = .FALSE.
!     LOGICAL  , PARAMETER :: DO_DEBUG_SOURCES = .TRUE.

      INTEGER  , PARAMETER :: LOG_UNIT = 67

!  Factors (from an earlier version)

      INTEGER  , PARAMETER :: Y = -1

      REAL(FPK), PARAMETER :: DECFAC = 1.0D0

!  Start
!  =====

!  initialise output

      FAIL    = .FALSE.
      MESSAGE = ' '

!  Overall single scatter correction flag
!    This must now be set for each PI solution
!     Bug detected and fixed, 23 July 2008
!     For SSCORR_OUTGOING, All SS contributions must be turned off

!  Elastic-only flag:
!    RRS Terms not present for Fourier > 2

      DO_ELASTIC_ONLY = .NOT. DO_RRS_OVERALL .OR. &
                            ( DO_RRS_OVERALL .AND. FOURIER .GT. 2 )

!  flux multiplier;  which is correct factor?
!        Turns out to be DECFAC = 1.0

      IF ( DO_ENERGY_BALANCING ) THEN
        FM_NONSHIFT         = ONE
        FM_NONSHIFT_DIRECT  = FM_NONSHIFT
        FM_NONSHIFT_DIFFUSE = FM_NONSHIFT * DECFAC
      ENDIF

!  Fourier proxy

      M = FOURIER

!  Number of Legendre moments = 2

      IF ( NSTREAMS .EQ. 1 .AND. NLAYERS .EQ. 1 ) THEN
        NRRS_MOMENTS = 0
      ELSE
        NRRS_MOMENTS = 2
      ENDIF

!  calculation point
!  For the binning, we need the inner and outer loops

      WW = W_EXCIT
      IF ( DO_BIN_REALIZATION ) THEN
        APT = ACTUAL_POINT
        CPT = APT + OFFSET_INNER
      ELSE
        APT = ACTUAL_POINT
        CPT = W_EXCIT
      ENDIF

!  Number of points or bins, Raman gain terms

      IF ( DO_BIN_REALIZATION ) THEN
        NPOINTS_LOCAL = N_RRSBINS ( APT )
      ELSE
        NPOINTS_LOCAL = NPOINTS_MONO
      ENDIF

!  Weighting function control

!mick fix 7/20/2016 - moved defining of DO_ATMOS_JACOBIANS from if blocks below to here

      DO_ATMOS_JACOBIANS   =  DO_PROFILE_WFS .OR.  DO_COLUMN_WFS
      DO_SURFACE_JACOBIANS =  DO_SURFACE_WFS .AND. DO_INCLUDE_SURFACE

!  Atmospheric profile linearization control

      IF ( DO_PROFILE_WFS ) THEN
        !DO_ATMOS_JACOBIANS = .TRUE.
        NC  = 0
        NC2 = 0
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            KS(N)    = 1
            KF(N)    = N
            KS2(N)   = 1
            KF2(N)   = NLAYERS
            KPARS(N) = LAYER_VARY_NUMBER(N)
!            DO K = 1, N
!              NC = NC + 1
!              NKSTORAGE(N,K) = NC
!            ENDDO
!            DO K = 1, NLAYERS
!              NC2 = NC2 + 1
!              NKSTORAGE2(N,K) = NC2
!            ENDDO
          ELSE
            KS(N)    = 1
            KF(N)    = 0
            KS2(N)   = 1
            KF2(N)   = 0
            KPARS(N) = 0
          ENDIF
        ENDDO
      ENDIF

!  Atmospheric column linearization control

      IF ( DO_COLUMN_WFS ) THEN
        !DO_ATMOS_JACOBIANS = .TRUE.
        DO N = 1, NLAYERS
          KS(N)    = 0
          KF(N)    = 0
          KS2(N)   = 0
          KF2(N)   = 0
          KPARS(0) = N_TOTALCOLUMN_WFS
        ENDDO
      ENDIF

!  Dummy values of INITRANS

      DUMTRANS = ONE
      IF ( DO_ATMOS_JACOBIANS ) THEN
        DO K = 0, NLAYERS
          DO Q = 1, MAX_ATMOSWFS
            L_DUMTRANS(Q,K) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  local IOPS for Elastic scattering only

      IF ( DO_ELASTIC_ONLY ) THEN
        DO N = 1, NLAYERS
          DO L = 0, NMOMENTS
            OMEGA_PHASMOMS(N,L) = OMEGAMOMS_ELASTIC(N,L,CPT)
          ENDDO
        ENDDO
        IF ( DO_PROFILE_WFS ) THEN
          DO N = 1, NLAYERS
           DO Q = 1, KPARS(N)
            DO L = 0, NMOMENTS
             L_OMEGA_PHASMOMS(Q,N,L) = L_OMEGAMOMS_ELASTIC(Q,N,L,CPT)
            ENDDO
           ENDDO
          ENDDO
        ENDIF
        IF ( DO_COLUMN_WFS ) THEN
          DO N = 1, NLAYERS
           DO Q = 1, KPARS(0)
            DO L = 0, NMOMENTS
             L_OMEGA_PHASMOMS(Q,N,L) = L_OMEGAMOMS_ELASTIC(Q,N,L,CPT)
            ENDDO
           ENDDO
          ENDDO
        ENDIF
      ENDIF

!  local IOPS for INELASTIC scattering (Energy balancing)

      IF ( .NOT.DO_ELASTIC_ONLY ) THEN
        IF ( DO_ENERGY_BALANCING ) THEN
         DO N = 1, NLAYERS
           DO L = 0, NMOMENTS
            OMEGA_PHASMOMS(N,L) = OMEGAMOMS_ELASTIC(N,L,CPT)
           ENDDO
         ENDDO
         IF ( DO_PROFILE_WFS ) THEN
          DO N = 1, NLAYERS
           DO Q = 1, KPARS(N)
            DO L = 0, NMOMENTS
             L_OMEGA_PHASMOMS(Q,N,L) = L_OMEGAMOMS_ELASTIC(Q,N,L,CPT)
            ENDDO
           ENDDO
          ENDDO
         ENDIF
         IF ( DO_COLUMN_WFS ) THEN
          DO N = 1, NLAYERS
           DO Q = 1, KPARS(0)
            DO L = 0, NMOMENTS
             L_OMEGA_PHASMOMS(Q,N,L) = L_OMEGAMOMS_ELASTIC(Q,N,L,CPT)
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDIF
       ENDIF

!  local IOPS for INELASTIC scattering (Cabannes-Raman)
!   RobFix Bug 5/27/11: L_OMEGAMOMS_CABANNES must have APT index

      IF ( .NOT.DO_ELASTIC_ONLY ) THEN
        IF ( .NOT.DO_ENERGY_BALANCING ) THEN
         DO N = 1, NLAYERS
           DO L = 0, NMOMENTS
            OMEGA_PHASMOMS(N,L) = OMEGAMOMS_CABANNES(N,L,APT)
           ENDDO
         ENDDO
         IF ( DO_PROFILE_WFS ) THEN
          DO N = 1, NLAYERS
           DO Q = 1, KPARS(N)
            DO L = 0, NMOMENTS
             L_OMEGA_PHASMOMS(Q,N,L) = L_OMEGAMOMS_CABANNES(Q,N,L,APT)
! Old         L_OMEGA_PHASMOMS(Q,N,L) = L_OMEGAMOMS_CABANNES(Q,N,L,CPT)
            ENDDO
           ENDDO
          ENDDO
         ENDIF
         IF ( DO_COLUMN_WFS ) THEN
          DO N = 1, NLAYERS
           DO Q = 1, KPARS(0)
            DO L = 0, NMOMENTS
             L_OMEGA_PHASMOMS(Q,N,L) = L_OMEGAMOMS_CABANNES(Q,N,L,APT)
! Old         L_OMEGA_PHASMOMS(Q,N,L) = L_OMEGAMOMS_CABANNES(Q,N,L,CPT)
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDIF
      ENDIF

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!         START LAYER LOOP
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      DO N = 1, NLAYERS  !(mick note: this layer loop ends @ ~line 3545!)

!  Local user-stream transmittances

        if ( do_user_streams ) then
          do ui = 1, n_user_streams
            trans_userm(ui) = SAVE_TRANS_USERM(UI,N,CPT)
          enddo
        endif

!  Local user-stream transmittances, linearized

        if ( do_atmos_jacobians ) then
         if ( do_user_streams ) then
          do ui = 1, n_user_streams
           do k = ks(n), kf(n)
            if ( k.eq.n .or. k.eq.0 ) then
             do q = 1, kpars(k)
              l_trans_userm(q,ui) = L_SAVE_TRANS_USERM(Q,UI,N,CPT)
             enddo
            endif
           enddo
          enddo
         endif
        endif

!  Local DELTAUS at CPT

        deltaus = DELTAU_VERT_INPUT(N,CPT)
        if ( do_atmos_jacobians ) then
          do k = ks(n), kf(n)
            if ( k.eq.n .or. k.eq.0 ) then
              do q = 1, kpars(k)
                l_deltaus(q) = L_DELTAU_VERT_INPUT(Q,N,CPT)
              enddo
            endif
          enddo
        endif

!  Initialize solution and exponent numbers
!    P  for the exponent number
!    Z  for the regular solution and Atmos linearization
!    ZS for the surface linearization sources

        P  = 0
        Z  = 0
        ZS = 0

!  Zero the Particular integral
!  ----------------------------

!  (profile linearization)

        if ( do_profile_wfs ) then
          DO K = 1, NLAYERS
            NK2 = NKSTORAGE2(N,K)
            DO Q = 1, KPARS(K)
              DO I = 1, 2*NSTREAMS
                L_RAMAN_WUPPER(Q,I,NK2) = ZERO
                L_RAMAN_WLOWER(Q,I,NK2) = ZERO
              ENDDO
              IF ( DO_UPWELLING ) THEN
                DO UM = 1, N_USER_STREAMS
                  L_WLAYER_PIST_UP(Q,NK2,UM) = ZERO
                ENDDO
              ENDIF
              IF ( DO_DNWELLING ) THEN
                DO UM = 1, N_USER_STREAMS
                  L_WLAYER_PIST_DN(Q,NK2,UM) = ZERO
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        endif

!  (column linearization)

        if ( do_column_wfs ) then
          NK2 = NKSTORAGE2(N,0)
          DO Q = 1, KPARS(0)
            DO I = 1, 2*NSTREAMS
              L_RAMAN_WUPPER(Q,I,NK2) = ZERO
              L_RAMAN_WLOWER(Q,I,NK2) = ZERO
            ENDDO
            IF ( DO_UPWELLING ) THEN
              DO UM = 1, N_USER_STREAMS
                L_WLAYER_PIST_UP(Q,NK2,UM) = ZERO
              ENDDO
            ENDIF
            IF ( DO_DNWELLING ) THEN
              DO UM = 1, N_USER_STREAMS
                L_WLAYER_PIST_DN(Q,NK2,UM) = ZERO
              ENDDO
            ENDIF
          ENDDO
        endif

!  (surface linearization). Updated Version 2.5 10/12/15

        if ( do_surface_wfs ) then
          DO I = 1, 2*NSTREAMS
            LS_RAMAN_WUPPER(1:N_SURFACEWFS,I,N) = ZERO
            LS_RAMAN_WLOWER(1:N_SURFACEWFS,I,N) = ZERO
          ENDDO
          IF ( DO_UPWELLING ) THEN
            DO UM = 1, N_USER_STREAMS
              LS_WLAYER_PIST_UP(1:N_SURFACEWFS,N,UM) = ZERO
            ENDDO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            DO UM = 1, N_USER_STREAMS
              LS_WLAYER_PIST_DN(1:N_SURFACEWFS,N,UM) = ZERO
            ENDDO
          ENDIF
        endif

!   #####################################
!   #####################################
!   #    E L A S T I C   T E R M        #
!   #####################################
!   #####################################

!  Source vector for the Elastic scattering term
!  =============================================

!  increase solution and exponent numbers by 1

        Z = Z + 1
        P = P + 1
        CHANGE_EXPONENT = .TRUE.
        FLIPPER         = .FALSE.
        PICUTOFF        = BEAM_PICUTOFF(CPT)
        DO_MULTIPLIER   = .TRUE.

!  Set single scatter overall flag, set freshly for each

        DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING
!        DO_SSCORR_OVERALL =
!     &     DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING

!  Step 1A. set transmittances and secants for this term
!           Need to re-set this each time.

        ASOURCE  = BEAM_AVSECANT(N,CPT)
        INITRANS = BEAM_ITRANS(N,CPT)
        DELTRANS = BEAM_DTRANS(N,CPT)

!  Step 1B. Atmospheric Linearization quantities

        IF ( DO_ATMOS_JACOBIANS ) THEN
          DO K = KS(N), KF(N)
            NK = NKSTORAGE(N,K)
            DO Q = 1, KPARS(K)
              L_ASOURCE(Q,K)  = L_BEAM_AVSECANT(Q,NK,CPT)
              L_INITRANS(Q,K) = L_BEAM_ITRANS(Q,NK,CPT)
              L_DELTRANS(Q,K) = L_BEAM_DTRANS(Q,NK,CPT)
            ENDDO
          ENDDO
        ENDIF

!  No step 2 here

!  Step 3A. Regular solution source vector, for quadrature streams

        DO I = 1, NSTREAMS
          SP = ZERO
          SM = ZERO
          DO L = FOURIER, NMOMENTS
            SP = SP + OMEGA_PHASMOMS(N,L) * PLM_00_PXI(I,L)
            SM = SM + OMEGA_PHASMOMS(N,L) * PLM_00_MXI(I,L)
          ENDDO
          QSOURCE_QPOS(I) = SP
          QSOURCE_QNEG(I) = SM
        ENDDO

!  Step 3B. Linearizations values of quadrature source vector
!   Safety first Way of Zeroing these arrays for non-specified entries

        IF ( DO_ATMOS_JACOBIANS ) THEN
          DO K = KS(N), KF(N)
            IF ( K.EQ.N .OR. K.EQ.0 ) THEN
              DO Q = 1, KPARS(K)
                DO I = 1, NSTREAMS
                  SP = ZERO
                  SM = ZERO
                  DO L = FOURIER, NMOMENTS
                    SP = SP + L_OMEGA_PHASMOMS(Q,N,L) * PLM_00_PXI(I,L)
                    SM = SM + L_OMEGA_PHASMOMS(Q,N,L) * PLM_00_MXI(I,L)
                  ENDDO
                  L_QSOURCE_QPOS(Q,K,I) = SP
                  L_QSOURCE_QNEG(Q,K,I) = SM
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, KPARS(K)
                DO I = 1, NSTREAMS
                  L_QSOURCE_QPOS(Q,K,I) = ZERO
                  L_QSOURCE_QNEG(Q,K,I) = ZERO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Step 4A. Regular solution source vector for user streams

        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            SP = ZERO
            SM = ZERO
            DO L = FOURIER, NMOMENTS
              SP = SP + OMEGA_PHASMOMS(N,L) * PLM_00_PXUI(UI,L)
              SM = SM + OMEGA_PHASMOMS(N,L) * PLM_00_MXUI(UI,L)
            ENDDO
            QSOURCE_UPOS1(UI) = SP
            QSOURCE_UNEG1(UI) = SM
          ENDDO
        ENDIF

!  Step 4B. Atmospheric Linearization source vector for user streams

        IF ( DO_USER_STREAMS.and.DO_ATMOS_JACOBIANS ) THEN
          DO K = KS(N), KF(N)
            IF ( K.EQ.N .OR. K.EQ.0 ) THEN
              DO Q = 1, KPARS(K)
                DO UI = 1, N_USER_STREAMS
                  SP = ZERO
                  SM = ZERO
                  DO L = FOURIER, NMOMENTS
                    SP = SP + L_OMEGA_PHASMOMS(Q,N,L)*PLM_00_PXUI(UI,L)
                    SM = SM + L_OMEGA_PHASMOMS(Q,N,L)*PLM_00_MXUI(UI,L)
                  ENDDO
                  L_QSOURCE_UPOS1(Q,K,UI) = SP
                  L_QSOURCE_UNEG1(Q,K,UI) = SM
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, KPARS(K)
                DO UI = 1, N_USER_STREAMS
                  L_QSOURCE_UPOS1(Q,K,UI) = ZERO
                  L_QSOURCE_UNEG1(Q,K,UI) = ZERO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  No surface linearization source vectors here (Step 4C)

!  Whole layer Greens function Discrete ordinate solutions

        CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE,             & ! Outputs
              BTERM_SAVE, MULT_M, MULT_P )                          ! Outputs

!  Whole layer Linearized Greens function solutions

        IF ( DO_ATMOS_JACOBIANS ) THEN
          CALL LPC_GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER,           & ! Inputs
              TAYLOR_SMALL, KS(N), KF(N), KPARS,              & ! Inputs
              NKSTORAGE, NKSTORAGE2, NSTREAMS, QUAD_WEIGHTS,  & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,           & ! Inputs
              ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P,         & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS, & ! Inputs
              DELTAUS, QSOURCE_QPOS, QSOURCE_QNEG,            & ! Inputs
              L_XPOS, L_KEIGEN, L_KTRANS, L_EIGENNORM_SAVE,   & ! Inputs
              L_ASOURCE, L_DELTRANS, L_INITRANS,              & ! Inputs
              L_DELTAUS, L_QSOURCE_QPOS, L_QSOURCE_QNEG,      & ! Inputs
              L_RAMAN_WUPPER, L_RAMAN_WLOWER,                 & ! Outputs
              L_ATERM_SAVE, L_BTERM_SAVE )                      ! Outputs
        ENDIF

!  No surface linearization solution, this point

!  Call to the post processing routine

        CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Call to the Linearized post processing routine

        IF ( DO_ATMOS_JACOBIANS ) THEN
          CALL LPC_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                                  & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                           & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                        & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),                & ! Inputs
            KS(N), KF(N), KPARS, NKSTORAGE2, TAYLOR_ORDER,               & ! Inputs
            TAYLOR_SMALL, TRANS_USERM, DELTAUS, USER_STREAMS,            & ! Inputs
            FLIPPER, CHANGE_EXPONENT, Z,                                 & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                       & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1,        & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,            & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,                              & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,                              & ! Inputs
            L_TRANS_USERM, L_DELTAUS, L_ASOURCE, L_DELTRANS, L_INITRANS, & ! Inputs
            L_QSOURCE_UPOS1, L_QSOURCE_UNEG1, L_ATERM_SAVE,L_BTERM_SAVE, & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,          & ! Inputs
            L_WLAYER_PIST_UP, L_WLAYER_PIST_DN )                           ! Outputs
        ENDIF

!  Debug
!      write(68,*)z ; write(69,*)z
!      do n = 1, 25
!         write(68,'(i4,1p8e19.10)')n,WLAYER_PIST_UP(n,1), WLAYER_PIST_DN(n,1)
!         write(69,'(i4,1p8e19.10)')n,L_WLAYER_PIST_UP(1,n,1), L_WLAYER_PIST_DN(1,n,1)
!      enddo
!      pause

!  Finish layer if elastic only

        IF ( DO_ELASTIC_ONLY ) GO TO 5222

!  For the Cabannes-Raman method, Skip the next section

        IF ( .NOT. DO_ENERGY_BALANCING ) GO TO 7676

!   #####################################
!   #####################################
!   #    R R S   L O S S   T E R M S    #
!   #####################################
!   #####################################

!  Source vectors and PI solutions for non-shifted Solar beam term
!  ===============================================================

!  General Procedure in 3 Steps :
!       (1) set INITRANS, ENDTRANS, ASOURCE (transmittances and secants)
!       (2) source vectors for quadrature streams
!       (3) source vectors for user-angle streams

!  increase solution number by 1, P not changed
!    --- Do the exponent again in this linearized version (November 2008

        Z = Z + 1
!      CHANGE_EXPONENT = .FALSE.               ! commented out 11/11/08
        CHANGE_EXPONENT = .TRUE.
        DO_MULTIPLIER   = .FALSE.
        DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING

!  Step 1... set transmittances and secants (once only for FOURIER = 0 )
!     Already done,

!  No step 2

!  Step 3A. set source vector for quadrature streams

        DO I = 1, NSTREAMS
          SP = ZERO
          SM = ZERO
          DO L = FOURIER, NRRS_MOMENTS
            SP = SP + OMEGAMOMS_RRSLOSS(N,L,APT) * PLM_00_PXI(I,L)
            SM = SM + OMEGAMOMS_RRSLOSS(N,L,APT) * PLM_00_MXI(I,L)
          ENDDO
          QSOURCE_QPOS(I) = - FM_NONSHIFT_DIRECT * SP
          QSOURCE_QNEG(I) = - FM_NONSHIFT_DIRECT * SM
        ENDDO

!  Step 3B. set Atmospheric Linearization Source vector for quadrature s

        IF ( DO_ATMOS_JACOBIANS ) THEN
          DO K = KS(N), KF(N)
            IF ( K.EQ.N .OR. K.EQ.0 ) THEN
              DO Q = 1, KPARS(K)
                DO I = 1, NSTREAMS
                  SP = ZERO
                  SM = ZERO
                  DO L = FOURIER, NRRS_MOMENTS
                   SP=SP+L_OMEGAMOMS_RRSLOSS(Q,N,L,APT)*PLM_00_PXI(I,L)
                   SM=SM+L_OMEGAMOMS_RRSLOSS(Q,N,L,APT)*PLM_00_MXI(I,L)
                  ENDDO
                  L_QSOURCE_QPOS(Q,K,I) = - FM_NONSHIFT_DIRECT * SP
                  L_QSOURCE_QNEG(Q,K,I) = - FM_NONSHIFT_DIRECT * SM
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, KPARS(K)
                DO I = 1, NSTREAMS
                  L_QSOURCE_QPOS(Q,K,I) = ZERO
                  L_QSOURCE_QNEG(Q,K,I) = ZERO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  No surface linearization source vectors here (Step 3C)

!  Step 4A. set Regular Solution source vector for user streams

        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            SP = ZERO
            SM = ZERO
            DO L = FOURIER, NRRS_MOMENTS
              SP = SP + OMEGAMOMS_RRSLOSS(N,L,APT)*PLM_00_PXUI(UI,L)
              SM = SM + OMEGAMOMS_RRSLOSS(N,L,APT)*PLM_00_MXUI(UI,L)
            ENDDO
            QSOURCE_UPOS1(UI) = - FM_NONSHIFT_DIRECT * SP
            QSOURCE_UNEG1(UI) = - FM_NONSHIFT_DIRECT * SM
          ENDDO
        ENDIF

!  Step 4B. set Atmospheric Linearization Source vector for user streams

        IF ( DO_USER_STREAMS.and.DO_ATMOS_JACOBIANS ) THEN
          DO K = KS(N), KF(N)
            IF ( K.EQ.N .OR. K.EQ.0 ) THEN
             DO Q = 1, KPARS(K)
              DO UI = 1, N_USER_STREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                SP=SP+L_OMEGAMOMS_RRSLOSS(Q,N,L,APT)*PLM_00_PXUI(UI,L)
                SM=SM+L_OMEGAMOMS_RRSLOSS(Q,N,L,APT)*PLM_00_MXUI(UI,L)
               ENDDO
               L_QSOURCE_UPOS1(Q,K,UI) = - FM_NONSHIFT_DIRECT * SP
               L_QSOURCE_UNEG1(Q,K,UI) = - FM_NONSHIFT_DIRECT * SM
              ENDDO
             ENDDO
            ELSE
             DO Q = 1, KPARS(K)
              DO UI = 1, N_USER_STREAMS
               L_QSOURCE_UPOS1(Q,K,UI) = ZERO
               L_QSOURCE_UNEG1(Q,K,UI) = ZERO
              ENDDO
             ENDDO
            ENDIF
          ENDDO
        ENDIF

!  No surface linearization source vectors here (Step 4C)

!  Whole layer Greens function Discrete ordinate solutions

        CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE,             & ! Outputs
              BTERM_SAVE, MULT_M, MULT_P )                          ! Outputs

!  Whole layer Linearized Greens function solutions

        IF ( DO_ATMOS_JACOBIANS ) THEN
          CALL LPC_GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER,           & ! Inputs
              TAYLOR_SMALL, KS(N), KF(N), KPARS,              & ! Inputs
              NKSTORAGE, NKSTORAGE2, NSTREAMS, QUAD_WEIGHTS,  & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,           & ! Inputs
              ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P,         & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS, & ! Inputs
              DELTAUS, QSOURCE_QPOS, QSOURCE_QNEG,            & ! Inputs
              L_XPOS, L_KEIGEN, L_KTRANS, L_EIGENNORM_SAVE,   & ! Inputs
              L_ASOURCE, L_DELTRANS, L_INITRANS,              & ! Inputs
              L_DELTAUS, L_QSOURCE_QPOS, L_QSOURCE_QNEG,      & ! Inputs
              L_RAMAN_WUPPER, L_RAMAN_WLOWER,                 & ! Outputs
              L_ATERM_SAVE, L_BTERM_SAVE )                      ! Outputs
        ENDIF

!  Debug
!      write(68,*)z
!      do n = 1, 25
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i+8,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i+8,n),i=1,8)
!      enddo
!      pause

!  Call to the post processing routine

        CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Call to the Linearized post processing routine

        IF ( DO_ATMOS_JACOBIANS ) THEN
          CALL LPC_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                                  & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                           & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                        & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),                & ! Inputs
            KS(N), KF(N), KPARS, NKSTORAGE2, TAYLOR_ORDER,               & ! Inputs
            TAYLOR_SMALL, TRANS_USERM, DELTAUS, USER_STREAMS,            & ! Inputs
            FLIPPER, CHANGE_EXPONENT, Z,                                 & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                       & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1,        & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,            & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,                              & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,                              & ! Inputs
            L_TRANS_USERM, L_DELTAUS, L_ASOURCE, L_DELTRANS, L_INITRANS, & ! Inputs
            L_QSOURCE_UPOS1, L_QSOURCE_UNEG1, L_ATERM_SAVE,L_BTERM_SAVE, & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,          & ! Inputs
            L_WLAYER_PIST_UP, L_WLAYER_PIST_DN )                           ! Outputs
        ENDIF

!  debug. 2/6/17
!        q = 1 ; k = 1 ; i = 3 ; nk = nkstorage2(n,k)
!        write(67,225)z,n,WLAYER_PIST_UP(N,1),WLAYER_PIST_DN(N,1)
!        write(67,225)z,n,RAMAN_WUPPER(i,n),RAMAN_WLOWER(i,n)
!        write(67,225)z,n,L_WLAYER_PIST_UP(q,NK,1),L_WLAYER_PIST_DN(q,NK,1)
!        write(67,225)z,n,L_RAMAN_WUPPER(q,i,nk),L_RAMAN_WLOWER(q,i,nk)

!  Sources and PIs for the non-shifted Diffuse-field loss-terms
!  ============================================================

!  Start RAMAN LOSS TERM clause (mick note: this clause ends @ ~line 2315)

        IF ( .NOT. DO_DIRECTRRS_ONLY ) THEN

!  (b1) Due to diffuse-term particular integral of elastic solution
!       -----------------------------------------------------------

!  General Procedure in 4 Steps :
!       (1) set INITRANS, ENDTRANS, ASOURCE (transmittances and secants)
!       (2) integrate solution for each moment = HOLD(L)
!       (3) source vectors for quadrature streams
!       (4) source vectors for user-angle streams

!  Increase number of solutions by 1, P not changed (same exponent)
!       MUST CHANGE EXPONENT for LINEARIZED

          Z = Z + 1
!          CHANGE_EXPONENT = .FALSE.
          CHANGE_EXPONENT = .TRUE.
          DO_MULTIPLIER   = .FALSE.
!          DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING   ! Wrong
          DO_SSCORR_OVERALL = .FALSE.

!  Step 1. (set transmittances and secants for this term)
!    This is the same as previous situation
!  Step 2.... (integrated solution contributions, array HOLD(L) )
!  Step 3.... (set source vector for quadrature streams)
!  Step 4.... (set source vector for user streams)

!  Step 2A. Holding terms, Regular solution

          DO L = FOURIER, NRRS_MOMENTS
            SP = ZERO
            SM = ZERO
            DO ID = 1, NSTREAMS
              ID1 = ID + NSTREAMS
              SP = SP + PLM_WT_PXI(ID,L) * L0_WPARTIC(ID1,N,WW)
              SM = SM + PLM_WT_MXI(ID,L) * L0_WPARTIC(ID,N,WW)
            ENDDO
            HOLD_MP(L)  = SP + SM
            HOLD_RRS(L) = OMEGAMOMS_RRSLOSS(N,L,APT) * HOLD_MP(L)
          ENDDO

!  Step 2B. Holding terms Atmospheric linearization
!      - No zeroing, Extra term if K = 0 (column) or K = N (profile)

          IF ( DO_ATMOS_JACOBIANS ) THEN
           DO K = KS(N), KF(N)
            NK = NKSTORAGE(N,K)
            DO Q = 1, KPARS(K)
             DO L = FOURIER, NRRS_MOMENTS
              SP = ZERO
              SM = ZERO
              DO ID = 1, NSTREAMS
               ID1 = ID + NSTREAMS
               SP = SP + PLM_WT_PXI(ID,L) * L_L0_WPARTIC(Q,ID1,NK,WW)
               SM = SM + PLM_WT_MXI(ID,L) * L_L0_WPARTIC(Q,ID,NK,WW)
              ENDDO
              L_HOLD_MP = SP + SM
              L_HOLD_RRS(Q,K,L) = OMEGAMOMS_RRSLOSS(N,L,APT)*L_HOLD_MP
              IF ( K.EQ.N .OR. K.EQ.0 ) THEN
                L_HOLD_RRS(Q,K,L) = L_HOLD_RRS(Q,K,L) + &
                      L_OMEGAMOMS_RRSLOSS(Q,N,L,APT) * HOLD_MP(L)
              ENDIF
             ENDDO
            ENDDO
           ENDDO
          ENDIF

!  Step 3A. Regular solution Source vector for quadrature streams

          DO I = 1, NSTREAMS
            SP = ZERO
            SM = ZERO
            DO L = FOURIER, NRRS_MOMENTS
              SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
              SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
            ENDDO
            QSOURCE_QPOS(I) = - FM_NONSHIFT_DIFFUSE * SP
            QSOURCE_QNEG(I) = - FM_NONSHIFT_DIFFUSE * SM
          ENDDO

!  Step 3B. Atmospheric Linearization Source vector for quadrature strea

          IF ( DO_ATMOS_JACOBIANS ) THEN
           DO K = KS(N), KF(N)
            DO Q = 1, KPARS(K)
             DO I = 1, NSTREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
               SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXI(I,L)
               SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXI(I,L)
              ENDDO
              L_QSOURCE_QPOS(Q,K,I) = - FM_NONSHIFT_DIFFUSE * SP
              L_QSOURCE_QNEG(Q,K,I) = - FM_NONSHIFT_DIFFUSE * SM
             ENDDO
            ENDDO
           ENDDO
          ENDIF

!  No surface linearization source vectors here (Step 3C)

!  Step 4A. Regular Solution Source vector for user streams

          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
                SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
                SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
              ENDDO
              QSOURCE_UPOS1(UI)  = - FM_NONSHIFT_DIFFUSE * SP
              QSOURCE_UNEG1(UI)  = - FM_NONSHIFT_DIFFUSE * SM
            ENDDO
          ENDIF

!  Step 4B. Atmospheric Linearization Source vector for user streams

          IF ( DO_ATMOS_JACOBIANS .and. DO_USER_STREAMS ) THEN
           DO K = KS(N), KF(N)
            DO Q = 1, KPARS(K)
             DO UI = 1, N_USER_STREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
               SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXUI(UI,L)
               SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXUI(UI,L)
              ENDDO
              L_QSOURCE_UPOS1(Q,K,UI)  = - FM_NONSHIFT_DIFFUSE * SP
              L_QSOURCE_UNEG1(Q,K,UI)  = - FM_NONSHIFT_DIFFUSE * SM
             ENDDO
            ENDDO
           ENDDO
          ENDIF

!  No surface linearization source vectors here (Step 4C)

!  Whole layer Greens function Discrete ordinate solutions

          CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE,             & ! Outputs
              BTERM_SAVE, MULT_M, MULT_P )                          ! Outputs

!  Whole layer Linearized Greens function solutions

          IF ( DO_ATMOS_JACOBIANS ) THEN
            CALL LPC_GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER,           & ! Inputs
              TAYLOR_SMALL, KS(N), KF(N), KPARS,              & ! Inputs
              NKSTORAGE, NKSTORAGE2, NSTREAMS, QUAD_WEIGHTS,  & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,           & ! Inputs
              ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P,         & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS, & ! Inputs
              DELTAUS, QSOURCE_QPOS, QSOURCE_QNEG,            & ! Inputs
              L_XPOS, L_KEIGEN, L_KTRANS, L_EIGENNORM_SAVE,   & ! Inputs
              L_ASOURCE, L_DELTRANS, L_INITRANS,              & ! Inputs
              L_DELTAUS, L_QSOURCE_QPOS, L_QSOURCE_QNEG,      & ! Inputs
              L_RAMAN_WUPPER, L_RAMAN_WLOWER,                 & ! Outputs
              L_ATERM_SAVE, L_BTERM_SAVE )                      ! Outputs
          ENDIF

!      write(68,*)z
!      do n = 1, 25
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i+8,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i+8,n),i=1,8)
!      enddo
!      pause

!  Call to the post processing routine

          CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Call to the Linearized post processing routine

          IF ( DO_ATMOS_JACOBIANS ) THEN
            CALL LPC_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                                  & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                           & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                        & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),                & ! Inputs
            KS(N), KF(N), KPARS, NKSTORAGE2, TAYLOR_ORDER,               & ! Inputs
            TAYLOR_SMALL, TRANS_USERM, DELTAUS, USER_STREAMS,            & ! Inputs
            FLIPPER, CHANGE_EXPONENT, Z,                                 & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                       & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1,        & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,            & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,                              & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,                              & ! Inputs
            L_TRANS_USERM, L_DELTAUS, L_ASOURCE, L_DELTRANS, L_INITRANS, & ! Inputs
            L_QSOURCE_UPOS1, L_QSOURCE_UNEG1, L_ATERM_SAVE,L_BTERM_SAVE, & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,          & ! Inputs
            L_WLAYER_PIST_UP, L_WLAYER_PIST_DN )                           ! Outputs
          ENDIF

!  debug. 2/6/17
!        q = 1 ; k = 1 ; i = 3 ; nk = nkstorage2(n,k)
!        write(67,225)z,n,WLAYER_PIST_UP(N,1),WLAYER_PIST_DN(N,1)
!        write(67,225)z,n,RAMAN_WUPPER(i,n),RAMAN_WLOWER(i,n)
!        write(67,225)z,n,L_WLAYER_PIST_UP(q,NK,1),L_WLAYER_PIST_DN(q,NK,1)
!        write(67,225)z,n,L_RAMAN_WUPPER(q,i,nk),L_RAMAN_WLOWER(q,i,nk)

!  (b2) Due to homogeneous functions of non-shifted elastic solution
!       ------------------------------------------------------------

!  start stream loop  (mick note: this stream loop ends @ ~line 2311)

          DO J = 1, NSTREAMS 

!  homogeneous DN solutions, increase Z and P by 1
!  ===============================================

           IF ( DO_SURFACE_JACOBIANS) ZS = ZS + 1
           Z  = Z + 1
           P  = P + 1
           CHANGE_EXPONENT = .TRUE.
           DO_MULTIPLIER   = .TRUE.
           FLIPPER         = .FALSE.
!           DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING     ! Wrong
           DO_SSCORR_OVERALL = .FALSE.

!  Step 1A. Transmittances

           PICUTOFF = NLAYERS
           INITRANS = ONE
           DELTRANS = L0_KTRANS(J,N,W_EXCIT)
           ASOURCE  = L0_KEIGEN(J,N,W_EXCIT)

!  Step 1B. Linearized Transmittances

           IF ( DO_ATMOS_JACOBIANS ) THEN
            DO K = KS2(N), KF2(N)
             DO Q = 1, KPARS(K)
              L_INITRANS(Q,K) = ZERO
              IF ( K.EQ.N .OR. K.EQ.0 ) THEN
               L_DELTRANS(Q,K) = L_L0_KTRANS(Q,J,N,W_EXCIT)
               L_ASOURCE(Q,K)  = L_L0_KEIGEN(Q,J,N,W_EXCIT)
              ELSE
               L_DELTRANS(Q,K) = ZERO
               L_ASOURCE(Q,K)  = ZERO
              ENDIF
             ENDDO
            ENDDO
           ENDIF

!  Step 2A. Regular Solution Holding terms

           LC = L0_WHOM_LCON(J,N,WW)
           DO L = FOURIER, NRRS_MOMENTS
            SP = ZERO
            SM = ZERO
            DO ID = 1, NSTREAMS
              ID1 = ID + NSTREAMS
              LP = PLM_WT_PXI(ID,L)
              LM = PLM_WT_MXI(ID,L)
              SP = SP + LC * LP * L0_WHOM_XPOS(ID1,J,N,WW)
              SM = SM + LC * LM * L0_WHOM_XPOS(ID,J,N,WW)
            ENDDO
            HOLD_MP(L)  = SP + SM
            HOLD_RRS(L) = OMEGAMOMS_RRSLOSS(N,L,APT) * HOLD_MP(L)
           ENDDO

!  Step 2B. Atmospheric Linearization Holding terms

           IF ( DO_ATMOS_JACOBIANS ) THEN
            DO K = KS2(N), KF2(N)
             NK = NKSTORAGE2(N,K)
             DO Q = 1, KPARS(K)
              L_LC = L_L0_WHOM_LCON(Q,J,NK,WW)
              DO L = FOURIER, NRRS_MOMENTS
               SP = ZERO
               SM = ZERO
               DO ID = 1, NSTREAMS
                ID1 = ID + NSTREAMS
                LP = PLM_WT_PXI(ID,L)
                LM = PLM_WT_MXI(ID,L)
                SP = SP + L_LC * L0_WHOM_XPOS(ID1,J,N,WW) * LP
                SM = SM + L_LC * L0_WHOM_XPOS(ID,J,N,WW)  * LM
               ENDDO
               L_HOLD_MP = SP + SM
               L_HOLD_RRS(Q,K,L) = &
                      +  OMEGAMOMS_RRSLOSS(N,L,APT) * L_HOLD_MP
               IF ( N.EQ.K .OR. K.EQ.0 ) THEN
                SP = ZERO
                SM = ZERO
                DO ID = 1, NSTREAMS
                 ID1 = ID + NSTREAMS
                 LP = PLM_WT_PXI(ID,L)
                 LM = PLM_WT_MXI(ID,L)
                 SP = SP + LC * L_L0_WHOM_XPOS(Q,ID1,J,N,WW) * LP
                 SM = SM + LC * L_L0_WHOM_XPOS(Q,ID,J,N,WW)  * LM
                ENDDO
                L_HOLD_MP = SP + SM
                L_HOLD_RRS(Q,K,L) = L_HOLD_RRS(Q,K,L) &
                      +    OMEGAMOMS_RRSLOSS(N,L,APT)   * L_HOLD_MP &
                      +  L_OMEGAMOMS_RRSLOSS(Q,N,L,APT) *   HOLD_MP(L)
               ENDIF
              ENDDO
             ENDDO
            ENDDO
           ENDIF

!  Step 2C. Surface Linearization Holding terms
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

           IF ( DO_SURFACE_JACOBIANS ) THEN
            DO Q = 1, N_SURFACEWFS
             LS_LC = LS_L0_WHOM_LCON(Q,J,N,WW)
             DO L = FOURIER, NRRS_MOMENTS
              SP = ZERO
              SM = ZERO
              DO ID = 1, NSTREAMS
               ID1 = ID + NSTREAMS
               LP = PLM_WT_PXI(ID,L)
               LM = PLM_WT_MXI(ID,L)
               SP = SP + LP * LS_LC * L0_WHOM_XPOS(ID1,J,N,WW)
               SM = SM + LM * LS_LC * L0_WHOM_XPOS(ID,J,N,WW)
              ENDDO
              LS_HOLD_MP = SP + SM
              LS_HOLD_RRS(Q,L) = OMEGAMOMS_RRSLOSS(N,L,APT) * LS_HOLD_MP
             ENDDO
            ENDDO
           ENDIF

!  Step 3A. Regular Solution Source vector for quadrature streams

           DO I = 1, NSTREAMS
            SP = ZERO
            SM = ZERO
            DO L = FOURIER, NRRS_MOMENTS
             SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
             SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
            ENDDO
            QSOURCE_QPOS(I) = - FM_NONSHIFT_DIFFUSE * SP
            QSOURCE_QNEG(I) = - FM_NONSHIFT_DIFFUSE * SM
           ENDDO

!  Step 3B. Atmospheric Linearization Source vector for quadrature strea

           IF ( DO_ATMOS_JACOBIANS ) THEN
            DO K = KS2(N), KF2(N)
             DO Q = 1, KPARS(K)
              DO I = 1, NSTREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXI(I,L)
                SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXI(I,L)
               ENDDO
               L_QSOURCE_QPOS(Q,K,I) = - FM_NONSHIFT_DIFFUSE * SP
               L_QSOURCE_QNEG(Q,K,I) = - FM_NONSHIFT_DIFFUSE * SM
              ENDDO
             ENDDO
            ENDDO
           ENDIF

!  Step 3C. Surface Linearization Source vector for quadrature streams
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

           IF ( DO_SURFACE_JACOBIANS ) THEN
            DO Q = 1, N_SURFACEWFS
             DO I = 1, NSTREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
               SP = SP + LS_HOLD_RRS(Q,L) * PLM_PXI(I,L)
               SM = SM + LS_HOLD_RRS(Q,L) * PLM_MXI(I,L)
              ENDDO
              LS_QSOURCE_QPOS(Q,I) = - FM_NONSHIFT_DIFFUSE * SP
              LS_QSOURCE_QNEG(Q,I) = - FM_NONSHIFT_DIFFUSE * SM
             ENDDO
            ENDDO
           ENDIF

!  Step 4A. Regular Solution Source vector for user streams

           IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
             SP = ZERO
             SM = ZERO
             DO L = FOURIER, NRRS_MOMENTS
              SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
              SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
             ENDDO
             QSOURCE_UPOS1(UI)  = - FM_NONSHIFT_DIFFUSE * SP
             QSOURCE_UNEG1(UI)  = - FM_NONSHIFT_DIFFUSE * SM
            ENDDO
           ENDIF

!  Step 4B. Atmospheric Linearization Source vector for user streams

           IF ( DO_ATMOS_JACOBIANS .and. DO_USER_STREAMS ) THEN
            DO K = KS2(N), KF2(N)
             DO Q = 1, KPARS(K)
              DO UI = 1, N_USER_STREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXUI(UI,L)
                SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXUI(UI,L)
               ENDDO
               L_QSOURCE_UPOS1(Q,K,UI)  = - FM_NONSHIFT_DIFFUSE * SP
               L_QSOURCE_UNEG1(Q,K,UI)  = - FM_NONSHIFT_DIFFUSE * SM
              ENDDO
             ENDDO
            ENDDO
           ENDIF

!  Step 4C. Surface Linearization Source vector for user streams
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

           IF ( DO_SURFACE_JACOBIANS .and. DO_USER_STREAMS ) THEN
            DO Q = 1, N_SURFACEWFS
             DO UI = 1, N_USER_STREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
               SP = SP + LS_HOLD_RRS(Q,L) * PLM_PXUI(UI,L)
               SM = SM + LS_HOLD_RRS(Q,L) * PLM_MXUI(UI,L)
              ENDDO
              LS_QSOURCE_UPOS1(Q,UI) = - FM_NONSHIFT_DIFFUSE * SP
              LS_QSOURCE_UNEG1(Q,UI) = - FM_NONSHIFT_DIFFUSE * SM
             ENDDO
            ENDDO
           ENDIF

!  Whole layer Greens function Discrete ordinate solutions

           CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE,             & ! Outputs
              BTERM_SAVE, MULT_M, MULT_P )                          ! Outputs

!  Whole layer Linearized Greens function solutions

           IF ( DO_ATMOS_JACOBIANS ) THEN
             CALL LPC_GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER,           & ! Inputs
              TAYLOR_SMALL, KS2(N), KF2(N), KPARS,            & ! Inputs
              NKSTORAGE, NKSTORAGE2, NSTREAMS, QUAD_WEIGHTS,  & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,           & ! Inputs
              ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P,         & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS, & ! Inputs
              DELTAUS, QSOURCE_QPOS, QSOURCE_QNEG,             & ! Inputs
              L_XPOS, L_KEIGEN, L_KTRANS, L_EIGENNORM_SAVE,   & ! Inputs
              L_ASOURCE, L_DELTRANS, L_INITRANS,              & ! Inputs
              L_DELTAUS, L_QSOURCE_QPOS, L_QSOURCE_QNEG,       & ! Inputs
              L_RAMAN_WUPPER, L_RAMAN_WLOWER,                 & ! Outputs
              L_ATERM_SAVE, L_BTERM_SAVE )                      ! Outputs
           ENDIF

!  Whole layer Surface Linearized Greens function solutions

           IF ( DO_SURFACE_JACOBIANS ) THEN
              CALL LS_GREENFUNC_SOLUTION_1 &
                ( M, N, ZS, NSTREAMS, N_SURFACEWFS, QUAD_WEIGHTS,         & ! Inputs
                  FLIPPER, PICUTOFF, INITRANS, XPOS, EIGENNORM_SAVE,      & ! Inputs
                  LS_QSOURCE_QPOS, LS_QSOURCE_QNEG, MULT_M, MULT_P,       & ! Inputs
                  LS_RAMAN_WUPPER, LS_RAMAN_WLOWER, LS_ATERM_SAVE, LS_BTERM_SAVE )    ! Outputs
           ENDIF

!  Debug
!      write(68,*)z, j
!      do n = 1, 25
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i+8,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i+8,n),i=1,8)
!      enddo
!      pause

!  Call to the post processing routine

           CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Call to the Linearized post processing routine

           IF ( DO_ATMOS_JACOBIANS ) THEN
             CALL LPC_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                                  & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                           & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                        & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),                & ! Inputs
            KS2(N), KF2(N), KPARS, NKSTORAGE2, TAYLOR_ORDER,             & ! Inputs
            TAYLOR_SMALL, TRANS_USERM, DELTAUS, USER_STREAMS,            & ! Inputs
            FLIPPER, CHANGE_EXPONENT, Z,                                 & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                       & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1,        & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,            & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,                              & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,                              & ! Inputs
            L_TRANS_USERM, L_DELTAUS, L_ASOURCE, L_DELTRANS, L_INITRANS, & ! Inputs
            L_QSOURCE_UPOS1, L_QSOURCE_UNEG1, L_ATERM_SAVE,L_BTERM_SAVE, & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,          & ! Inputs
            L_WLAYER_PIST_UP, L_WLAYER_PIST_DN )                           ! Outputs
           ENDIF

!  Debug 2/6/17
!        q = 1 ; k = 1 ; i = 3 ; nk = nkstorage2(n,k)
!        write(67,225)z,n,WLAYER_PIST_UP(N,1),WLAYER_PIST_DN(N,1)
!        write(67,225)z,n,RAMAN_WUPPER(i,n),RAMAN_WLOWER(i,n)
!        write(67,225)z,n,L_WLAYER_PIST_UP(q,NK,1),L_WLAYER_PIST_DN(q,NK,1)
!        write(67,225)z,n,L_RAMAN_WUPPER(q,i,nk),L_RAMAN_WLOWER(q,i,nk)

!  Call to the Surface-Linearized post processing routine

           IF ( DO_SURFACE_JACOBIANS ) THEN
             CALL LS_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                   & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL, FOURIER,   & ! Inputs
            N, NSTREAMS, N_USER_STREAMS, N_SURFACEWFS,    & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N), & ! Inputs
            ZS, INITRANS, U_XPOS, U_XNEG,                 & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,               & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,               & ! Inputs
            LS_ATERM_SAVE, LS_BTERM_SAVE,                 & ! Inputs
            LS_QSOURCE_UPOS1, LS_QSOURCE_UNEG1,           & ! Inputs
            LS_WLAYER_PIST_UP, LS_WLAYER_PIST_DN )          ! Outputs
           ENDIF

!  homogeneous UP solutions, increase Z and P by 1
!  ===============================================

           IF ( DO_SURFACE_JACOBIANS) ZS = ZS + 1
           Z = Z + 1
           P = P + 1
           CHANGE_EXPONENT = .TRUE.
           DO_MULTIPLIER   = .TRUE.
           FLIPPER         = .TRUE.
!           DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING         ! Wrong
           DO_SSCORR_OVERALL = .FALSE.

!  Step 1A. Transmittances

           PICUTOFF = NLAYERS
           INITRANS = ONE
           DELTRANS = L0_KTRANS(J,N,W_EXCIT)
           ASOURCE  = L0_KEIGEN(J,N,W_EXCIT)

!  Step 1B. Linearized Transmittances

           IF ( DO_ATMOS_JACOBIANS ) THEN
            DO K = KS2(N), KF2(N)
             DO Q = 1, KPARS(K)
              L_INITRANS(Q,K) = ZERO
              IF ( K.EQ.N .OR. K.EQ.0 ) THEN
               L_DELTRANS(Q,K) = L_L0_KTRANS(Q,J,N,W_EXCIT)
               L_ASOURCE(Q,K)  = L_L0_KEIGEN(Q,J,N,W_EXCIT)
              ELSE
               L_DELTRANS(Q,K) = ZERO
               L_ASOURCE(Q,K)  = ZERO
              ENDIF
             ENDDO
            ENDDO
           ENDIF

!  Step 2A. Regular Solution Holding terms

           MC = L0_WHOM_MCON(J,N,WW)
           DO L = FOURIER, NRRS_MOMENTS
            SP = ZERO
            SM = ZERO
            DO ID = 1, NSTREAMS
              ID1 = ID + NSTREAMS
              LP = PLM_WT_PXI(ID,L)
              LM = PLM_WT_MXI(ID,L)
              SP = SP + MC * LP * L0_WHOM_XPOS(ID,J,N,WW)
              SM = SM + MC * LM * L0_WHOM_XPOS(ID1,J,N,WW)
            ENDDO
            HOLD_MP(L)  = SP + SM
            HOLD_RRS(L) = OMEGAMOMS_RRSLOSS(N,L,APT) * HOLD_MP(L)
           ENDDO

!  Step 2B. Atmospheric Linearization Holding terms

           IF ( DO_ATMOS_JACOBIANS ) THEN
            DO K = KS2(N), KF2(N)
             NK = NKSTORAGE2(N,K)
             DO Q = 1, KPARS(K)
              L_MC = L_L0_WHOM_MCON(Q,J,NK,WW)
              DO L = FOURIER, NRRS_MOMENTS
               SP = ZERO
               SM = ZERO
               DO ID = 1, NSTREAMS
                ID1 = ID + NSTREAMS
                LP = PLM_WT_PXI(ID,L)
                LM = PLM_WT_MXI(ID,L)
                SP = SP + L_MC * L0_WHOM_XPOS(ID,J,N,WW)  * LP
                SM = SM + L_MC * L0_WHOM_XPOS(ID1,J,N,WW) * LM
               ENDDO
               L_HOLD_MP = SP + SM
               L_HOLD_RRS(Q,K,L) = &
                      +  OMEGAMOMS_RRSLOSS(N,L,APT) * L_HOLD_MP
               IF ( N.EQ.K .OR. K.EQ.0 ) THEN
                SP = ZERO
                SM = ZERO
                DO ID = 1, NSTREAMS
                 ID1 = ID + NSTREAMS
                 LP = PLM_WT_PXI(ID,L)
                 LM = PLM_WT_MXI(ID,L)
                 SP = SP + MC * L_L0_WHOM_XPOS(Q,ID,J,N,WW)  * LP
                 SM = SM + MC * L_L0_WHOM_XPOS(Q,ID1,J,N,WW) * LM
                ENDDO
                L_HOLD_MP = SP + SM
                L_HOLD_RRS(Q,K,L) = L_HOLD_RRS(Q,K,L) &
                      +    OMEGAMOMS_RRSLOSS(N,L,APT)   * L_HOLD_MP &
                      +  L_OMEGAMOMS_RRSLOSS(Q,N,L,APT) *   HOLD_MP(L)
               ENDIF
              ENDDO
             ENDDO
            ENDDO
           ENDIF

!  Step 2C. Surface Linearization Holding terms
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

           IF ( DO_SURFACE_JACOBIANS) THEN
            DO Q = 1, N_SURFACEWFS
            LS_MC = LS_L0_WHOM_MCON(Q,J,N,WW)
             DO L = FOURIER, NRRS_MOMENTS
              SP = ZERO
              SM = ZERO
              DO ID = 1, NSTREAMS
               ID1 = ID + NSTREAMS
               LP = PLM_WT_PXI(ID,L)
               LM = PLM_WT_MXI(ID,L)
               SP = SP + LP * LS_MC * L0_WHOM_XPOS(ID,J,N,WW)
               SM = SM + LM * LS_MC * L0_WHOM_XPOS(ID1,J,N,WW)
              ENDDO
              LS_HOLD_MP = SP + SM
              LS_HOLD_RRS(Q,L) = OMEGAMOMS_RRSLOSS(N,L,APT) * LS_HOLD_MP
             ENDDO
            ENDDO
           ENDIF

!  Step 3A. Regular Solution Source vector for quadrature streams

           DO I = 1, NSTREAMS
            SP = ZERO
            SM = ZERO
            DO L = FOURIER, NRRS_MOMENTS
             SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
             SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
            ENDDO
            QSOURCE_QPOS(I) = - FM_NONSHIFT_DIFFUSE * SP
            QSOURCE_QNEG(I) = - FM_NONSHIFT_DIFFUSE * SM
           ENDDO

!  Step 3B. Atmospheric Linearization Source vector for quadrature strea

           IF ( DO_ATMOS_JACOBIANS ) THEN
            DO K = KS2(N), KF2(N)
             DO Q = 1, KPARS(K)
              DO I = 1, NSTREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXI(I,L)
                SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXI(I,L)
               ENDDO
               L_QSOURCE_QPOS(Q,K,I) = - FM_NONSHIFT_DIFFUSE * SP
               L_QSOURCE_QNEG(Q,K,I) = - FM_NONSHIFT_DIFFUSE * SM
              ENDDO
             ENDDO
            ENDDO
           ENDIF

!  Step 3C. Surface Linearization Source vector for quadrature streams
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

           IF ( DO_SURFACE_JACOBIANS ) THEN
            DO Q = 1, N_SURFACEWFS
             DO I = 1, NSTREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
               SP = SP + LS_HOLD_RRS(Q,L) * PLM_PXI(I,L)
               SM = SM + LS_HOLD_RRS(Q,L) * PLM_MXI(I,L)
              ENDDO
              LS_QSOURCE_QPOS(Q,I) = - FM_NONSHIFT_DIFFUSE * SP
              LS_QSOURCE_QNEG(Q,I) = - FM_NONSHIFT_DIFFUSE * SM
             ENDDO
            ENDDO
           ENDIF

!  Step 4A. Regular Solution Source vector for user streams

           IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
             SP = ZERO
             SM = ZERO
             DO L = FOURIER, NRRS_MOMENTS
              SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
              SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
             ENDDO
             QSOURCE_UPOS1(UI)  = - FM_NONSHIFT_DIFFUSE * SP
             QSOURCE_UNEG1(UI)  = - FM_NONSHIFT_DIFFUSE * SM
            ENDDO
           ENDIF

!  Step 4B. Atmospheric Linearization Source vector for user streams

           IF ( DO_ATMOS_JACOBIANS .and. DO_USER_STREAMS ) THEN
            DO K = KS2(N), KF2(N)
             DO Q = 1, KPARS(K)
              DO UI = 1, N_USER_STREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXUI(UI,L)
                SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXUI(UI,L)
               ENDDO
               L_QSOURCE_UPOS1(Q,K,UI)  = - FM_NONSHIFT_DIFFUSE * SP
               L_QSOURCE_UNEG1(Q,K,UI)  = - FM_NONSHIFT_DIFFUSE * SM
              ENDDO
             ENDDO
            ENDDO
           ENDIF

!  Step 4C. Surface Linearization Source vector for user streams
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

           IF ( DO_SURFACE_JACOBIANS .and. DO_USER_STREAMS ) THEN
            DO Q = 1, N_SURFACEWFS
             DO UI = 1, N_USER_STREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
               SP = SP + LS_HOLD_RRS(Q,L) * PLM_PXUI(UI,L)
               SM = SM + LS_HOLD_RRS(Q,L) * PLM_MXUI(UI,L)
              ENDDO
              LS_QSOURCE_UPOS1(Q,UI) = - FM_NONSHIFT_DIFFUSE * SP
              LS_QSOURCE_UNEG1(Q,UI) = - FM_NONSHIFT_DIFFUSE * SM
             ENDDO
            ENDDO
           ENDIF

!  Whole layer Greens function Discrete ordinate solutions

           CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE,             & ! Outputs
              BTERM_SAVE, MULT_M, MULT_P )                          ! Outputs

!  Whole layer Linearized Greens function solutions

           IF ( DO_ATMOS_JACOBIANS ) THEN
             CALL LPC_GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER,           & ! Inputs
              TAYLOR_SMALL, KS2(N), KF2(N), KPARS,            & ! Inputs
              NKSTORAGE, NKSTORAGE2, NSTREAMS, QUAD_WEIGHTS,  & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,           & ! Inputs
              ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P,         & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS, & ! Inputs
              DELTAUS, QSOURCE_QPOS, QSOURCE_QNEG,             & ! Inputs
              L_XPOS, L_KEIGEN, L_KTRANS, L_EIGENNORM_SAVE,   & ! Inputs
              L_ASOURCE, L_DELTRANS, L_INITRANS,              & ! Inputs
              L_DELTAUS, L_QSOURCE_QPOS, L_QSOURCE_QNEG,       & ! Inputs
              L_RAMAN_WUPPER, L_RAMAN_WLOWER,                 & ! Outputs
              L_ATERM_SAVE, L_BTERM_SAVE )                      ! Outputs
           ENDIF

!  Whole layer Surface Linearized Greens function solutions

           IF ( DO_SURFACE_JACOBIANS ) THEN
              CALL LS_GREENFUNC_SOLUTION_1 &
                ( M, N, ZS, NSTREAMS, N_SURFACEWFS, QUAD_WEIGHTS,         & ! Inputs
                  FLIPPER, PICUTOFF, INITRANS, XPOS, EIGENNORM_SAVE,      & ! Inputs
                  LS_QSOURCE_QPOS, LS_QSOURCE_QNEG, MULT_M, MULT_P,       & ! Inputs
                  LS_RAMAN_WUPPER, LS_RAMAN_WLOWER, LS_ATERM_SAVE, LS_BTERM_SAVE )    ! Outputs
           ENDIF

!  Debug
!      write(68,*)z, j
!      do n = 1, 25
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i+8,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i+8,n),i=1,8)
!      enddo
!      pause

!  Call to the post processing routine

           CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Call to the Linearized post processing routine

           IF ( DO_ATMOS_JACOBIANS ) THEN
             CALL LPC_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                                  & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                           & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                        & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),                & ! Inputs
            KS2(N), KF2(N), KPARS, NKSTORAGE2, TAYLOR_ORDER,             & ! Inputs
            TAYLOR_SMALL, TRANS_USERM, DELTAUS, USER_STREAMS,            & ! Inputs
            FLIPPER, CHANGE_EXPONENT, Z,                                 & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                       & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1,        & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,            & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,                              & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,                              & ! Inputs
            L_TRANS_USERM, L_DELTAUS, L_ASOURCE, L_DELTRANS, L_INITRANS, & ! Inputs
            L_QSOURCE_UPOS1, L_QSOURCE_UNEG1, L_ATERM_SAVE,L_BTERM_SAVE, & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,          & ! Inputs
            L_WLAYER_PIST_UP, L_WLAYER_PIST_DN )                           ! Outputs
           ENDIF

!  Call to the Surface-Linearized post processing routine

           IF ( DO_SURFACE_JACOBIANS ) THEN
             CALL LS_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                   & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL, FOURIER,   & ! Inputs
            N, NSTREAMS, N_USER_STREAMS, N_SURFACEWFS,    & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N), & ! Inputs
            ZS, INITRANS, U_XPOS, U_XNEG,                 & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,               & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,               & ! Inputs
            LS_ATERM_SAVE, LS_BTERM_SAVE,                 & ! Inputs
            LS_QSOURCE_UPOS1, LS_QSOURCE_UNEG1,           & ! Inputs
            LS_WLAYER_PIST_UP, LS_WLAYER_PIST_DN )          ! Outputs
           ENDIF

!  debug. 2/6/17
! 225   format(2i4,1p2e25.16)
!        q = 1 ; k = 1 ; i = 3 ; nk = nkstorage2(n,k)
!        write(67,225)z,n,WLAYER_PIST_UP(N,1),WLAYER_PIST_DN(N,1)
!        write(67,225)z,n,RAMAN_WUPPER(i,n),RAMAN_WLOWER(i,n)
!        write(67,225)z,n,L_WLAYER_PIST_UP(q,NK,1),L_WLAYER_PIST_DN(q,NK,1)
!        write(67,225)z,n,L_RAMAN_WUPPER(q,i,nk),L_RAMAN_WLOWER(q,i,nk)

!  End stream loop  (mick note: this stream loop starts @ ~line 1650)

          ENDDO

!  End RAMAN LOSS TERM clause  (mick note: this clause starts @ ~line 1432)

        ENDIF

!     if ( n.eq.2)stop 'end loss'

!  Continuation point for avoiding loss terms.

 7676   CONTINUE

!   #####################################
!   #####################################
!   #    R R S  R A M A N   T E R M S   #
!   #####################################
!   #####################################

!  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  For each transition..   The Raman-shifted contributions
!  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!  Start the transitions loop  (mick note: this transitions loop ends @ ~line 3537)

        DO S = 1, NPOINTS_LOCAL 

!  BPT is the Bin pointer       in the  binning      realization
!  BPT is the wavelength itself in the monochromatic realization

         IF ( DO_BIN_REALIZATION ) THEN
          BPT = BINMAP ( S, APT )
         ELSE
          BPT = S
         ENDIF

!  In the monochromatic case, avoid for the actual point

         IF ( DO_MONO_REALIZATION .AND. S.EQ.CPT ) GO TO 5656

!  copy the OMEGAMOMS array (binned or mono)

         IF ( DO_BIN_REALIZATION ) THEN
           DO L = 0, 2
             OMEGAMOMS_RRSLOCAL(L) = OMEGAMOMS_RRSBIN(N,L,S,APT)
           ENDDO
         ELSE
           DO L = 0, 2
             OMEGAMOMS_RRSLOCAL(L) = OMEGAMOMS_RRSGAIN(N,L,S)
           ENDDO
         ENDIF

!  copy the Linearized OMEGAMOMS array (binned or mono)

         IF ( DO_ATMOS_JACOBIANS ) THEN
           IF ( DO_BIN_REALIZATION ) THEN
             DO L = 0, 2
               DO Q = 1, LAYER_VARY_NUMBER(N)
                 L_OMEGAMOMS_RRSLOCAL(Q,L) = &
                     L_OMEGAMOMS_RRSBIN(Q,N,L,S,APT)
               ENDDO
             ENDDO
           ELSE
             DO L = 0, 2
               DO Q = 1, LAYER_VARY_NUMBER(N)
                 L_OMEGAMOMS_RRSLOCAL(Q,L) = &
                     L_OMEGAMOMS_RRSGAIN(Q,N,L,S)
               ENDDO
             ENDDO
           ENDIF
         ENDIF

!  flux multiplier;  which is correct factor?
!        Turns out to be DECFAC = 1.0

         FM_SHIFTED         = &
              FLUXES_RANKED ( BPT ) / FLUXES_RANKED ( CPT )
         FM_SHIFTED_DIRECT  = FM_SHIFTED
         FM_SHIFTED_DIFFUSE = FM_SHIFTED * DECFAC

!  optical thickness scaling for each shift

         DELTA_SCALING = DELTAU_VERT_INPUT(N,BPT) / &
                              DELTAU_VERT_INPUT(N,CPT)

!  Linearization of optical thickness scaling

         if ( do_atmos_jacobians ) then
           do k = ks(n), kf(n)
             if ( k.eq.n .or. k.eq.0 ) then
               do q = 1, kpars(k)
                 L_DELTA_SCALING(Q) = &
               ( L_DELTAU_VERT_INPUT(Q,N,BPT) - &
                 L_DELTAU_VERT_INPUT(Q,N,CPT)*DELTA_SCALING ) / &
                             DELTAU_VERT_INPUT(N,CPT)
               enddo
             endif
           enddo
         endif

!  Source vectors for the shifted Solar beam term
!  ==============================================

!  General Procedure in 3 Steps :
!       (1) set INITRANS, ENDTRANS, ASOURCE (transmittances and secants)
!       (2) source vectors for quadrature streams
!       (3) source vectors for user-angle streams

!  increase solution and exponent numbers by 1

         Z = Z + 1
         P = P + 1
         CHANGE_EXPONENT = .TRUE.
         DO_MULTIPLIER   = .TRUE.
         FLIPPER         = .FALSE.
         DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING

!  Step 1A. set transmittances and secants
!     Need to re-set this for each Fourier component

         PICUTOFF = BEAM_PICUTOFF(BPT)
         SAVCUT   = PICUTOFF
         ASOURCE  = BEAM_AVSECANT(N,BPT) * DELTA_SCALING
         INITRANS = BEAM_ITRANS(N,BPT)
         DELTRANS = BEAM_DTRANS(N,BPT)

!  Step 1B. Linearizations of transmittances and secants

        IF ( DO_ATMOS_JACOBIANS ) THEN
         DO K = KS(N), KF(N)
          NK = NKSTORAGE(N,K)
          DO Q = 1, KPARS(K)
           L_ASOURCE(Q,K)  = L_BEAM_AVSECANT(Q,NK,BPT) * DELTA_SCALING
           L_INITRANS(Q,K) = L_BEAM_ITRANS(Q,NK,BPT)
           L_DELTRANS(Q,K) = L_BEAM_DTRANS(Q,NK,BPT)
          ENDDO
          IF ( K.EQ.N .OR. K.EQ.0 ) THEN
           DO Q = 1, KPARS(K)
            L_ASOURCE(Q,K) = L_ASOURCE(Q,K) + &
                    BEAM_AVSECANT(N,BPT) * L_DELTA_SCALING(Q)
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  debug. 4 December 2006.
!          write(33,*)s,z,p,delta_scaling,picutoff(p)
!          write(33,*)s,z,p,delta_scaling,picutoff(p)
!        if (s.eq.5.or.s.eq.2)then
!          do n = 1, nlayers
!           write(*,*)S,Z,P,delta_scaling,ASOURCE,INITRANS,DELTRANS
!          enddo
!        endif
!        if (S.eq.npoints_local)pause

!  No Step 2.

!  Step 3A, Regular solution source vector for quadrature streams

         DO I = 1, NSTREAMS
           SP = ZERO
           SM = ZERO
           DO L = FOURIER, NRRS_MOMENTS
             PHIFUNC = OMEGAMOMS_RRSLOCAL(L)
             SP = SP + PHIFUNC * PLM_00_PXI(I,L)
             SM = SM + PHIFUNC * PLM_00_MXI(I,L)
           ENDDO
           QSOURCE_QPOS(I) = + FM_SHIFTED_DIRECT * SP
           QSOURCE_QNEG(I) = + FM_SHIFTED_DIRECT * SM
         ENDDO

!  Step 3B. set Atmospheric Linearization Source vector for quadrature s

         IF ( DO_ATMOS_JACOBIANS ) THEN
          DO K = KS(N), KF(N)
            IF ( K.EQ.N .OR. K.EQ.0 ) THEN
              DO Q = 1, KPARS(K)
                DO I = 1, NSTREAMS
                  SP = ZERO
                  SM = ZERO
                  DO L = FOURIER, NRRS_MOMENTS
                    PHIFUNC = L_OMEGAMOMS_RRSLOCAL(Q,L)
                    SP = SP + PHIFUNC * PLM_00_PXI(I,L)
                    SM = SM + PHIFUNC * PLM_00_MXI(I,L)
                  ENDDO
                  L_QSOURCE_QPOS(Q,K,I) = + FM_SHIFTED_DIRECT * SP
                  L_QSOURCE_QNEG(Q,K,I) = + FM_SHIFTED_DIRECT * SM
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, KPARS(K)
                DO I = 1, NSTREAMS
                  L_QSOURCE_QPOS(Q,K,I) = ZERO
                  L_QSOURCE_QNEG(Q,K,I) = ZERO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
         ENDIF

!  No surface linearization source vectors here (Step 3C)

!  Step 4A, Regular solution source vector for user streams

         IF ( DO_USER_STREAMS ) THEN
           DO UI = 1, N_USER_STREAMS
             SP = ZERO
             SM = ZERO
             DO L = FOURIER, NRRS_MOMENTS
               PHIFUNC = OMEGAMOMS_RRSLOCAL(L)
               SP = SP + PHIFUNC * PLM_00_PXUI(UI,L)
               SM = SM + PHIFUNC * PLM_00_MXUI(UI,L)
             ENDDO
             QSOURCE_UPOS1(UI) = + FM_SHIFTED_DIRECT * SP
             QSOURCE_UNEG1(UI) = + FM_SHIFTED_DIRECT * SM
           ENDDO
         ENDIF

!  Step 4B. set Atmospheric Linearization Source vector for user streams

         IF ( DO_USER_STREAMS.and.DO_ATMOS_JACOBIANS ) THEN
          DO K = KS(N), KF(N)
            IF ( K.EQ.N .OR. K.EQ.0 ) THEN
             DO Q = 1, KPARS(K)
              DO UI = 1, N_USER_STREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                  PHIFUNC = L_OMEGAMOMS_RRSLOCAL(Q,L)
                  SP = SP + PHIFUNC * PLM_00_PXUI(UI,L)
                  SM = SM + PHIFUNC * PLM_00_MXUI(UI,L)
               ENDDO
               L_QSOURCE_UPOS1(Q,K,UI) = + FM_SHIFTED_DIRECT * SP
               L_QSOURCE_UNEG1(Q,K,UI) = + FM_SHIFTED_DIRECT * SM
              ENDDO
             ENDDO
            ELSE
             DO Q = 1, KPARS(K)
              DO UI = 1, N_USER_STREAMS
               L_QSOURCE_UPOS1(Q,K,UI) = ZERO
               L_QSOURCE_UNEG1(Q,K,UI) = ZERO
              ENDDO
             ENDDO
            ENDIF
          ENDDO
         ENDIF

!  No surface linearization source vectors here (Step 4C)

!  Whole layer Greens function Discrete ordinate solutions

         CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE,             & ! Outputs
              BTERM_SAVE, MULT_M, MULT_P )                          ! Outputs

!  Whole layer Linearized Greens function solutions

         IF ( DO_ATMOS_JACOBIANS ) THEN

!mick fix 10/19/2015 - replaced KS2(N) & KF2(N) with KS(N) & KF(N)

           CALL LPC_GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER,           & ! Inputs
              !TAYLOR_SMALL, KS2(N), KF2(N), KPARS,           & ! Inputs
              TAYLOR_SMALL, KS(N), KF(N), KPARS,              & ! Inputs
              NKSTORAGE, NKSTORAGE2, NSTREAMS, QUAD_WEIGHTS,  & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,           & ! Inputs
              ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P,         & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS, & ! Inputs
              DELTAUS, QSOURCE_QPOS, QSOURCE_QNEG,            & ! Inputs
              L_XPOS, L_KEIGEN, L_KTRANS, L_EIGENNORM_SAVE,   & ! Inputs
              L_ASOURCE, L_DELTRANS, L_INITRANS,              & ! Inputs
              L_DELTAUS, L_QSOURCE_QPOS, L_QSOURCE_QNEG,      & ! Inputs
              L_RAMAN_WUPPER, L_RAMAN_WLOWER,                 & ! Outputs
              L_ATERM_SAVE, L_BTERM_SAVE )                      ! Outputs
         ENDIF

!  Debug
!      write(68,*)z
!      do n = 1, 25
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i+8,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i+8,n),i=1,8)
!      enddo
!      pause

!  Call to the post processing routine

         CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Call to the Linearized post processing routine

         IF ( DO_ATMOS_JACOBIANS ) THEN
!mick fix 10/19/2015 - replaced KS2(N) & KF2(N) with KS(N) & KF(N)
           CALL LPC_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                                  & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                           & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                        & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),                & ! Inputs
            !KS2(N), KF2(N), KPARS, NKSTORAGE2, TAYLOR_ORDER,            & ! Inputs
            KS(N), KF(N), KPARS, NKSTORAGE2, TAYLOR_ORDER,               & ! Inputs
            TAYLOR_SMALL, TRANS_USERM, DELTAUS, USER_STREAMS,            & ! Inputs
            FLIPPER, CHANGE_EXPONENT, Z,                                 & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                       & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1,        & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,            & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,                              & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,                              & ! Inputs
            L_TRANS_USERM, L_DELTAUS, L_ASOURCE, L_DELTRANS, L_INITRANS, & ! Inputs
            L_QSOURCE_UPOS1, L_QSOURCE_UNEG1, L_ATERM_SAVE,L_BTERM_SAVE, & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,          & ! Inputs
            L_WLAYER_PIST_UP, L_WLAYER_PIST_DN )                           ! Outputs
         ENDIF

!  Source vectors for the shifted Diffuse-field gain-terms
!  =======================================================

         IF ( .NOT. DO_DIRECTRRS_ONLY ) THEN

!  (b1) Due to particular integral of shifted elastic solution
!       ------------------------------------------------------

!  General Procedure in 4 Steps :
!       (1) set INITRANS, ENDTRANS, ASOURCE (transmittances and secants)
!       (2) integrate solution for each moment = HOLD(L)
!       (3) source vectors for quadrature streams
!       (4) source vectors for user-angle streams

!  increase solution number by 1, do not increase P

           Z = Z + 1
!           CHANGE_EXPONENT = .FALSE.                         ! Safety
           CHANGE_EXPONENT = .TRUE.
           DO_MULTIPLIER   = .FALSE.
!           DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING            ! Wrong
           DO_SSCORR_OVERALL = .FALSE.

!  Step 1.... (set transmittances and secants for this term)
!    This is the same as previous, only need to do this for Fourier 0

!  Step 2A, regular solution holding array

           DO L = FOURIER, NRRS_MOMENTS
             SP = ZERO
             SM = ZERO
             DO ID = 1, NSTREAMS
               ID1 = ID + NSTREAMS
               SP = SP + PLM_WT_PXI(ID,L) * L0_WPARTIC(ID1,N,BPT)
               SM = SM + PLM_WT_MXI(ID,L) * L0_WPARTIC(ID,N,BPT)
             ENDDO
             PHIFUNC     = OMEGAMOMS_RRSLOCAL(L)
             HOLD_MP(L)  = SP + SM
             HOLD_RRS(L) = PHIFUNC * HOLD_MP(L)
           ENDDO

!  Step 2B. Holding terms Atmospheric linearization
!      - No zeroing, Extra term if K = 0 (column) or K = N (profile)

           IF ( DO_ATMOS_JACOBIANS ) THEN
            DO K = KS(N), KF(N)
             NK = NKSTORAGE(N,K)
             DO Q = 1, KPARS(K)
              DO L = FOURIER, NRRS_MOMENTS
               SP = ZERO
               SM = ZERO
               DO ID = 1, NSTREAMS
                ID1 = ID + NSTREAMS
                SP = SP + PLM_WT_PXI(ID,L) * L_L0_WPARTIC(Q,ID1,NK,BPT)
                SM = SM + PLM_WT_MXI(ID,L) * L_L0_WPARTIC(Q,ID,NK,BPT)
               ENDDO
               L_HOLD_MP = SP + SM
               L_HOLD_RRS(Q,K,L) = OMEGAMOMS_RRSLOCAL(L) * L_HOLD_MP
               IF ( K.EQ.N .OR. K.EQ.0 ) THEN
                L_HOLD_RRS(Q,K,L) = L_HOLD_RRS(Q,K,L) + &
                        L_OMEGAMOMS_RRSLOCAL(Q,L) * HOLD_MP(L)
               ENDIF
              ENDDO
             ENDDO
            ENDDO
           ENDIF

!  Step 3A. Regular solution source vector, quadrature streams

           DO I = 1, NSTREAMS
             SP = ZERO
             SM = ZERO
             DO L = FOURIER, NRRS_MOMENTS
               SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
               SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
             ENDDO
             QSOURCE_QPOS(I) = + FM_SHIFTED_DIFFUSE * SP
             QSOURCE_QNEG(I) = + FM_SHIFTED_DIFFUSE * SM
           ENDDO

!  Step 3B. Atmospheric linearization source vector, quadrature streams

           IF ( DO_ATMOS_JACOBIANS ) THEN
            DO K = KS(N), KF(N)
             DO Q = 1, KPARS(K)
              DO I = 1, NSTREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXI(I,L)
                SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXI(I,L)
               ENDDO
               L_QSOURCE_QPOS(Q,K,I) = + FM_SHIFTED_DIFFUSE * SP
               L_QSOURCE_QNEG(Q,K,I) = + FM_SHIFTED_DIFFUSE * SM
              ENDDO
             ENDDO
            ENDDO
           ENDIF

!  No surface linearization source vectors here (Step 3C)

!  Step 4A. Regular solution source vector, user streams

           IF ( DO_USER_STREAMS ) THEN
             DO UI = 1, N_USER_STREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                 SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
                 SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
               ENDDO
               QSOURCE_UPOS1(UI)  = + FM_SHIFTED_DIFFUSE * SP
               QSOURCE_UNEG1(UI)  = + FM_SHIFTED_DIFFUSE * SM
             ENDDO
           ENDIF

!  Step 4B. Atmospheric linearization source vector, user streams

           IF ( DO_ATMOS_JACOBIANS .and. DO_USER_STREAMS ) THEN
            DO K = KS(N), KF(N)
             DO Q = 1, KPARS(K)
              DO UI = 1, N_USER_STREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXUI(UI,L)
                SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXUI(UI,L)
               ENDDO
               L_QSOURCE_UPOS1(Q,K,UI)  = + FM_SHIFTED_DIFFUSE * SP
               L_QSOURCE_UNEG1(Q,K,UI)  = + FM_SHIFTED_DIFFUSE * SM
              ENDDO
             ENDDO
            ENDDO
           ENDIF

!  No surface linearization source vectors here (Step 4C)

!  Whole layer Greens function Discrete ordinate solutions

           CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE,             & ! Outputs
              BTERM_SAVE, MULT_M, MULT_P )                          ! Outputs

!  Whole layer Linearized Greens function solutions

           IF ( DO_ATMOS_JACOBIANS ) THEN
             CALL LPC_GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER,           & ! Inputs
              TAYLOR_SMALL, KS(N), KF(N), KPARS,              & ! Inputs
              NKSTORAGE, NKSTORAGE2, NSTREAMS, QUAD_WEIGHTS,  & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,           & ! Inputs
              ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P,         & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS, & ! Inputs
              DELTAUS, QSOURCE_QPOS, QSOURCE_QNEG,             & ! Inputs
              L_XPOS, L_KEIGEN, L_KTRANS, L_EIGENNORM_SAVE,   & ! Inputs
              L_ASOURCE, L_DELTRANS, L_INITRANS,              & ! Inputs
              L_DELTAUS, L_QSOURCE_QPOS, L_QSOURCE_QNEG,       & ! Inputs
              L_RAMAN_WUPPER, L_RAMAN_WLOWER,                 & ! Outputs
              L_ATERM_SAVE, L_BTERM_SAVE )                      ! Outputs
           ENDIF

!      write(68,*)z
!      do n = 1, 25
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i+8,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i+8,n),i=1,8)
!      enddo
!      pause

!  Call to the post processing routine

           CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Call to the Linearized post processing routine

           IF ( DO_ATMOS_JACOBIANS ) THEN
             CALL LPC_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                                  & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                           & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                        & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),                & ! Inputs
            KS(N), KF(N), KPARS, NKSTORAGE2, TAYLOR_ORDER,               & ! Inputs
            TAYLOR_SMALL, TRANS_USERM, DELTAUS, USER_STREAMS,            & ! Inputs
            FLIPPER, CHANGE_EXPONENT, Z,                                 & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                       & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1,        & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,            & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,                              & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,                              & ! Inputs
            L_TRANS_USERM, L_DELTAUS, L_ASOURCE, L_DELTRANS, L_INITRANS, & ! Inputs
            L_QSOURCE_UPOS1, L_QSOURCE_UNEG1, L_ATERM_SAVE,L_BTERM_SAVE, & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,          & ! Inputs
            L_WLAYER_PIST_UP, L_WLAYER_PIST_DN )                           ! Outputs
           ENDIF

!  (b2) Due to homogeneous functions of shifted elastic solution
!       --------------------------------------------------------

!  start stream loop  (mick note: this stream loop ends @ ~line 3525)

           DO J = 1, NSTREAMS

!  homogeneous DN solutions, increase Z and P by 1
!  ===============================================

            Z = Z + 1
            P = P + 1
            IF ( DO_SURFACE_JACOBIANS ) ZS = ZS + 1
            CHANGE_EXPONENT = .TRUE.
            DO_MULTIPLIER   = .TRUE.
            FLIPPER         = .FALSE.
!            DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING             !  Wrong
            DO_SSCORR_OVERALL = .FALSE.

!  Step 1A. Transmittances and exponent
!            --------Note scaling on the optical thickness !!

!            PICUTOFF = NLAYERS
            PICUTOFF = savcut
            INITRANS = ONE
            DELTRANS = L0_KTRANS(J,N,BPT)
            ASOURCE  = L0_KEIGEN(J,N,BPT) * DELTA_SCALING

!  Step 1B. Linearized Transmittances

            IF ( DO_ATMOS_JACOBIANS ) THEN
             DO K = KS2(N), KF2(N)
              DO Q = 1, KPARS(K)
               L_INITRANS(Q,K) = ZERO
               IF ( K.EQ.N .OR. K.EQ.0 ) THEN
                L_DELTRANS(Q,K) = L_L0_KTRANS(Q,J,N,BPT)
                L_ASOURCE(Q,K)  = &
                         L_L0_KEIGEN(Q,J,N,BPT) * DELTA_SCALING &
                       +   L0_KEIGEN(J,N,BPT)    * L_DELTA_SCALING(Q)
               ELSE
                L_DELTRANS(Q,K) = ZERO
                L_ASOURCE(Q,K)  = ZERO
               ENDIF
              ENDDO
             ENDDO
            ENDIF

!  Step 2A. Regular solution holding terms

            LC = L0_WHOM_LCON(J,N,BPT)
            DO L = FOURIER, NRRS_MOMENTS
              SP = ZERO
              SM = ZERO
              DO ID = 1, NSTREAMS
                ID1 = ID + NSTREAMS
                LP = PLM_WT_PXI(ID,L)
                LM = PLM_WT_MXI(ID,L)
                SP = SP + LC * LP * L0_WHOM_XPOS(ID1,J,N,BPT)
                SM = SM + LC * LM * L0_WHOM_XPOS(ID,J,N,BPT)
              ENDDO
              PHIFUNC     = OMEGAMOMS_RRSLOCAL(L)
              HOLD_MP(L)  = SP + SM
              HOLD_RRS(L) = PHIFUNC * HOLD_MP(L)
            ENDDO

!  Step 2B. Atmospheric Linearization holding terms

            IF ( DO_ATMOS_JACOBIANS ) THEN
             DO K = KS2(N), KF2(N)
              NK = NKSTORAGE2(N,K)
              DO Q = 1, KPARS(K)
               L_LC = L_L0_WHOM_LCON(Q,J,NK,BPT)
               DO L = FOURIER, NRRS_MOMENTS
                SP = ZERO
                SM = ZERO
                DO ID = 1, NSTREAMS
                 ID1 = ID + NSTREAMS
                 LP = PLM_WT_PXI(ID,L)
                 LM = PLM_WT_MXI(ID,L)
                 SP = SP + L_LC * L0_WHOM_XPOS(ID1,J,N,BPT) * LP
                 SM = SM + L_LC * L0_WHOM_XPOS(ID,J,N,BPT)  * LM
                ENDDO
                L_HOLD_MP = SP + SM
                L_HOLD_RRS(Q,K,L) = OMEGAMOMS_RRSLOCAL(L) * L_HOLD_MP
                IF ( N.EQ.K .OR. K.EQ.0 ) THEN
                 SP = ZERO
                 SM = ZERO
                 DO ID = 1, NSTREAMS
                  ID1 = ID + NSTREAMS
                  LP = PLM_WT_PXI(ID,L)
                  LM = PLM_WT_MXI(ID,L)
                  SP = SP + LC * L_L0_WHOM_XPOS(Q,ID1,J,N,BPT) * LP
                  SM = SM + LC * L_L0_WHOM_XPOS(Q,ID,J,N,BPT)  * LM
                 ENDDO
                 L_HOLD_MP = SP + SM
                 L_HOLD_RRS(Q,K,L) = L_HOLD_RRS(Q,K,L) &
                      +    OMEGAMOMS_RRSLOCAL(L)   * L_HOLD_MP &
                      +  L_OMEGAMOMS_RRSLOCAL(Q,L) *   HOLD_MP(L)
                ENDIF
               ENDDO
              ENDDO
             ENDDO
            ENDIF

!  Step 2C. Surface Linearization Holding terms
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

            IF ( DO_SURFACE_JACOBIANS  ) THEN
             DO Q = 1, N_SURFACEWFS
              LS_LC = LS_L0_WHOM_LCON(Q,J,N,BPT)
              DO L = FOURIER, NRRS_MOMENTS
               SP = ZERO
               SM = ZERO
               DO ID = 1, NSTREAMS
                ID1 = ID + NSTREAMS
                LP = PLM_WT_PXI(ID,L)
                LM = PLM_WT_MXI(ID,L)
                SP = SP + LP * LS_LC * L0_WHOM_XPOS(ID1,J,N,BPT)
                SM = SM + LM * LS_LC * L0_WHOM_XPOS(ID,J,N,BPT)
               ENDDO
               LS_HOLD_MP = SP + SM
               LS_HOLD_RRS(Q,L) = OMEGAMOMS_RRSLOCAL(L) * LS_HOLD_MP
              ENDDO
             ENDDO
            ENDIF

!  Step 3A. Regular solution source term quadrature streams

            DO I = 1, NSTREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
                SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
                SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
              ENDDO
              QSOURCE_QPOS(I) = + FM_SHIFTED_DIFFUSE * SP
              QSOURCE_QNEG(I) = + FM_SHIFTED_DIFFUSE * SM
            ENDDO

!  Step 3B. Atmospheric Linearization Source term quadrature streams

            IF ( DO_ATMOS_JACOBIANS ) THEN
             DO K = KS2(N), KF2(N)
              DO Q = 1, KPARS(K)
               DO I = 1, NSTREAMS
                SP = ZERO
                SM = ZERO
                DO L = FOURIER, NRRS_MOMENTS
                 SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXI(I,L)
                 SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXI(I,L)
                ENDDO
                L_QSOURCE_QPOS(Q,K,I) = + FM_SHIFTED_DIFFUSE * SP
                L_QSOURCE_QNEG(Q,K,I) = + FM_SHIFTED_DIFFUSE * SM
               ENDDO
              ENDDO
             ENDDO
            ENDIF

!  Step 3C. Surface Linearization Source vector for quadrature streams
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

            IF ( DO_SURFACE_JACOBIANS  ) THEN
             DO Q = 1, N_SURFACEWFS
              DO I = 1, NSTREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                SP = SP + LS_HOLD_RRS(Q,L) * PLM_PXI(I,L)
                SM = SM + LS_HOLD_RRS(Q,L) * PLM_MXI(I,L)
               ENDDO
               LS_QSOURCE_QPOS(Q,I) = + FM_SHIFTED_DIFFUSE * SP
               LS_QSOURCE_QNEG(Q,I) = + FM_SHIFTED_DIFFUSE * SM
              ENDDO
             ENDDO
            ENDIF

!  Step 4A. Regular solution source vector user streams

            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                SP = ZERO
                SM = ZERO
                DO L = FOURIER, NRRS_MOMENTS
                  SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
                  SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
                ENDDO
                QSOURCE_UPOS1(UI) = + FM_SHIFTED_DIFFUSE * SP
                QSOURCE_UNEG1(UI) = + FM_SHIFTED_DIFFUSE * SM
              ENDDO
            ENDIF

!  Step 4B. Atmospheric Linearization Source vector for user streams

            IF ( DO_ATMOS_JACOBIANS .and. DO_USER_STREAMS ) THEN
             DO K = KS2(N), KF2(N)
              DO Q = 1, KPARS(K)
               DO UI = 1, N_USER_STREAMS
                SP = ZERO
                SM = ZERO
                DO L = FOURIER, NRRS_MOMENTS
                 SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXUI(UI,L)
                 SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXUI(UI,L)
                ENDDO
                L_QSOURCE_UPOS1(Q,K,UI)  = + FM_SHIFTED_DIFFUSE * SP
                L_QSOURCE_UNEG1(Q,K,UI)  = + FM_SHIFTED_DIFFUSE * SM
               ENDDO
              ENDDO
             ENDDO
            ENDIF

!  Step 4C. Surface Linearization Source vector for user streams
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

            IF ( DO_SURFACE_JACOBIANS .and. DO_USER_STREAMS ) THEN
             DO Q = 1, N_SURFACEWFS
              DO UI = 1, N_USER_STREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                SP = SP + LS_HOLD_RRS(Q,L) * PLM_PXUI(UI,L)
                SM = SM + LS_HOLD_RRS(Q,L) * PLM_MXUI(UI,L)
               ENDDO
               LS_QSOURCE_UPOS1(Q,UI) = + FM_SHIFTED_DIFFUSE * SP
               LS_QSOURCE_UNEG1(Q,UI) = + FM_SHIFTED_DIFFUSE * SM
              ENDDO
             ENDDO
            ENDIF

!  Whole layer Greens function Discrete ordinate solutions

            CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE,             & ! Outputs
              BTERM_SAVE, MULT_M, MULT_P )                          ! Outputs

!  Whole layer Linearized Greens function solutions

            IF ( DO_ATMOS_JACOBIANS ) THEN
              CALL LPC_GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER,           & ! Inputs
              TAYLOR_SMALL, KS2(N), KF2(N), KPARS,            & ! Inputs
              NKSTORAGE, NKSTORAGE2, NSTREAMS, QUAD_WEIGHTS,  & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,           & ! Inputs
              ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P,         & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS, & ! Inputs
              DELTAUS, QSOURCE_QPOS, QSOURCE_QNEG,             & ! Inputs
              L_XPOS, L_KEIGEN, L_KTRANS, L_EIGENNORM_SAVE,   & ! Inputs
              L_ASOURCE, L_DELTRANS, L_INITRANS,              & ! Inputs
              L_DELTAUS, L_QSOURCE_QPOS, L_QSOURCE_QNEG,       & ! Inputs
              L_RAMAN_WUPPER, L_RAMAN_WLOWER,                 & ! Outputs
              L_ATERM_SAVE, L_BTERM_SAVE )                      ! Outputs
            ENDIF

!  Whole layer Surface Linearized Greens function solutions

            IF ( DO_SURFACE_JACOBIANS ) THEN
              CALL LS_GREENFUNC_SOLUTION_1 &
                ( M, N, ZS, NSTREAMS, N_SURFACEWFS, QUAD_WEIGHTS,         & ! Inputs
                  FLIPPER, PICUTOFF, INITRANS, XPOS, EIGENNORM_SAVE,      & ! Inputs
                  LS_QSOURCE_QPOS, LS_QSOURCE_QNEG, MULT_M, MULT_P,       & ! Inputs
                  LS_RAMAN_WUPPER, LS_RAMAN_WLOWER, LS_ATERM_SAVE, LS_BTERM_SAVE )    ! Outputs
            ENDIF

!  Debug
!      write(68,*)z, j
!      do n = 1, 25
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i+8,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i+8,n),i=1,8)
!      enddo
!      pause

!  Call to the post processing routine

            CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Call to the Linearized post processing routine

            IF ( DO_ATMOS_JACOBIANS ) THEN
              CALL LPC_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                                  & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                           & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                        & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),                & ! Inputs
            KS2(N), KF2(N), KPARS, NKSTORAGE2, TAYLOR_ORDER,             & ! Inputs
            TAYLOR_SMALL, TRANS_USERM, DELTAUS, USER_STREAMS,            & ! Inputs
            FLIPPER, CHANGE_EXPONENT, Z,                                 & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                       & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1,        & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,            & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,                              & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,                              & ! Inputs
            L_TRANS_USERM, L_DELTAUS, L_ASOURCE, L_DELTRANS, L_INITRANS, & ! Inputs
            L_QSOURCE_UPOS1, L_QSOURCE_UNEG1, L_ATERM_SAVE,L_BTERM_SAVE, & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,          & ! Inputs
            L_WLAYER_PIST_UP, L_WLAYER_PIST_DN )                           ! Outputs
            ENDIF

!  Call to the Surface-Linearized post processing routine

            IF ( DO_SURFACE_JACOBIANS ) THEN
              CALL LS_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                   & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL, FOURIER,   & ! Inputs
            N, NSTREAMS, N_USER_STREAMS, N_SURFACEWFS,    & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N), & ! Inputs
            ZS, INITRANS, U_XPOS, U_XNEG,                 & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,               & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,               & ! Inputs
            LS_ATERM_SAVE, LS_BTERM_SAVE,                 & ! Inputs
            LS_QSOURCE_UPOS1, LS_QSOURCE_UNEG1,           & ! Inputs
            LS_WLAYER_PIST_UP, LS_WLAYER_PIST_DN )          ! Outputs
            ENDIF

!  homogeneous UP solutions, increase Z and P by 1
!  ===============================================

            Z = Z + 1
            P = P + 1
            IF ( DO_SURFACE_JACOBIANS ) ZS = ZS + 1
            CHANGE_EXPONENT = .TRUE.
            DO_MULTIPLIER   = .TRUE.
            FLIPPER         = .TRUE.
!            DO_SSCORR_OVERALL = DO_SSCORR_OUTGOING            ! Wrong
            DO_SSCORR_OVERALL = .FALSE.

!  Step 1A. Transmittances and exponents

!            PICUTOFF = NLAYERS
            PICUTOFF = savcut
            INITRANS = ONE
            DELTRANS = L0_KTRANS(J,N,BPT)
            ASOURCE  = L0_KEIGEN(J,N,BPT) * DELTA_SCALING

!  Step 1B. LInearized Transmittances

            IF ( DO_ATMOS_JACOBIANS ) THEN
             DO K = KS2(N), KF2(N)
              DO Q = 1, KPARS(K)
               L_INITRANS(Q,K) = ZERO
               IF ( K.EQ.N .OR. K.EQ.0 ) THEN
                L_DELTRANS(Q,K) = L_L0_KTRANS(Q,J,N,BPT)
                L_ASOURCE(Q,K)  = &
                         L_L0_KEIGEN(Q,J,N,BPT) * DELTA_SCALING &
                       +   L0_KEIGEN(J,N,BPT)    * L_DELTA_SCALING(Q)
               ELSE
                L_DELTRANS(Q,K) = ZERO
                L_ASOURCE(Q,K)  = ZERO
               ENDIF
              ENDDO
             ENDDO
            ENDIF

!  Step 2A. Regular solution holding terms

            MC = L0_WHOM_MCON(J,N,BPT)
            DO L = FOURIER, NRRS_MOMENTS
              SP = ZERO
              SM = ZERO
              DO ID = 1, NSTREAMS
                ID1 = ID + NSTREAMS
                LP = PLM_WT_PXI(ID,L)
                LM = PLM_WT_MXI(ID,L)
                SP = SP + MC * LP * L0_WHOM_XPOS(ID,J,N,BPT)
                SM = SM + MC * LM * L0_WHOM_XPOS(ID1,J,N,BPT)
              ENDDO
              PHIFUNC     = OMEGAMOMS_RRSLOCAL(L)
              HOLD_MP(L)  = SP + SM
              HOLD_RRS(L) = PHIFUNC * HOLD_MP(L)
            ENDDO

!  Step 2B. Atmospheric Linearization holding terms

            IF ( DO_ATMOS_JACOBIANS ) THEN
             DO K = KS2(N), KF2(N)
              NK = NKSTORAGE2(N,K)
              DO Q = 1, KPARS(K)
               L_MC = L_L0_WHOM_MCON(Q,J,NK,BPT)
               DO L = FOURIER, NRRS_MOMENTS
                SP = ZERO
                SM = ZERO
                DO ID = 1, NSTREAMS
                 ID1 = ID + NSTREAMS
                 LP = PLM_WT_PXI(ID,L)
                 LM = PLM_WT_MXI(ID,L)
                 SP = SP + L_MC * L0_WHOM_XPOS(ID,J,N,BPT) * LP
                 SM = SM + L_MC * L0_WHOM_XPOS(ID1,J,N,BPT)  * LM
                ENDDO
                L_HOLD_MP = SP + SM
                L_HOLD_RRS(Q,K,L) = OMEGAMOMS_RRSLOCAL(L) * L_HOLD_MP
                IF ( N.EQ.K .OR. K.EQ.0 ) THEN
                 SP = ZERO
                 SM = ZERO
                 DO ID = 1, NSTREAMS
                  ID1 = ID + NSTREAMS
                  LP = PLM_WT_PXI(ID,L)
                  LM = PLM_WT_MXI(ID,L)
                  SP = SP + MC * L_L0_WHOM_XPOS(Q,ID,J,N,BPT) * LP
                  SM = SM + MC * L_L0_WHOM_XPOS(Q,ID1,J,N,BPT)  * LM
                 ENDDO
                 L_HOLD_MP = SP + SM
                 L_HOLD_RRS(Q,K,L) = L_HOLD_RRS(Q,K,L) &
                      +    OMEGAMOMS_RRSLOCAL(L)   * L_HOLD_MP &
                      +  L_OMEGAMOMS_RRSLOCAL(Q,L) *   HOLD_MP(L)
                ENDIF
               ENDDO
              ENDDO
             ENDDO
            ENDIF

!  Step 2C. Surface Linearization Holding terms
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

            IF ( DO_SURFACE_JACOBIANS  ) THEN
             DO Q = 1, N_SURFACEWFS
              LS_MC = LS_L0_WHOM_MCON(Q,J,N,BPT)
              DO L = FOURIER, NRRS_MOMENTS
               SP = ZERO
               SM = ZERO
               DO ID = 1, NSTREAMS
                ID1 = ID + NSTREAMS
                LP = PLM_WT_PXI(ID,L)
                LM = PLM_WT_MXI(ID,L)
                SP = SP + LP * LS_MC * L0_WHOM_XPOS(ID,J,N,BPT)
                SM = SM + LM * LS_MC * L0_WHOM_XPOS(ID1,J,N,BPT)
               ENDDO
               LS_HOLD_MP = SP + SM
               LS_HOLD_RRS(Q,L) = OMEGAMOMS_RRSLOCAL(L) * LS_HOLD_MP
              ENDDO
             ENDDO
            ENDIF

!  Step 3A. Regular solution source term quadrature streams

            DO I = 1, NSTREAMS
              SP = ZERO
              SM = ZERO
              DO L = FOURIER, NRRS_MOMENTS
                SP = SP + HOLD_RRS(L) * PLM_PXI(I,L)
                SM = SM + HOLD_RRS(L) * PLM_MXI(I,L)
              ENDDO
              QSOURCE_QPOS(I) = + FM_SHIFTED_DIFFUSE * SP
              QSOURCE_QNEG(I) = + FM_SHIFTED_DIFFUSE * SM
            ENDDO

!  Step 3B. Atmospheric Linearization Source term quadrature streams

            IF ( DO_ATMOS_JACOBIANS ) THEN
             DO K = KS2(N), KF2(N)
              DO Q = 1, KPARS(K)
               DO I = 1, NSTREAMS
                SP = ZERO
                SM = ZERO
                DO L = FOURIER, NRRS_MOMENTS
                 SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXI(I,L)
                 SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXI(I,L)
                ENDDO
                L_QSOURCE_QPOS(Q,K,I) = + FM_SHIFTED_DIFFUSE * SP
                L_QSOURCE_QNEG(Q,K,I) = + FM_SHIFTED_DIFFUSE * SM
               ENDDO
              ENDDO
             ENDDO
            ENDIF

!  Step 3C. Surface Linearization Source vector for quadrature streams
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

            IF ( DO_SURFACE_JACOBIANS ) THEN
             DO Q = 1, N_SURFACEWFS
              DO I = 1, NSTREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                SP = SP + LS_HOLD_RRS(Q,L) * PLM_PXI(I,L)
                SM = SM + LS_HOLD_RRS(Q,L) * PLM_MXI(I,L)
               ENDDO
               LS_QSOURCE_QPOS(Q,I) = + FM_SHIFTED_DIFFUSE * SP
               LS_QSOURCE_QNEG(Q,I) = + FM_SHIFTED_DIFFUSE * SM
              ENDDO
             ENDDO
            ENDIF

!  Step 4A. Regular solution source vector user streams

            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                SP = ZERO
                SM = ZERO
                DO L = FOURIER, NRRS_MOMENTS
                  SP = SP + HOLD_RRS(L) * PLM_PXUI(UI,L)
                  SM = SM + HOLD_RRS(L) * PLM_MXUI(UI,L)
                ENDDO
                QSOURCE_UPOS1(UI) = + FM_SHIFTED_DIFFUSE * SP
                QSOURCE_UNEG1(UI) = + FM_SHIFTED_DIFFUSE * SM
              ENDDO
            ENDIF

!  Step 4B. Atmospheric Linearization Source vector for user streams

            IF ( DO_ATMOS_JACOBIANS .and. DO_USER_STREAMS ) THEN
             DO K = KS2(N), KF2(N)
              DO Q = 1, KPARS(K)
               DO UI = 1, N_USER_STREAMS
                SP = ZERO
                SM = ZERO
                DO L = FOURIER, NRRS_MOMENTS
                 SP = SP + L_HOLD_RRS(Q,K,L) * PLM_PXUI(UI,L)
                 SM = SM + L_HOLD_RRS(Q,K,L) * PLM_MXUI(UI,L)
                ENDDO
                L_QSOURCE_UPOS1(Q,K,UI)  = + FM_SHIFTED_DIFFUSE * SP
                L_QSOURCE_UNEG1(Q,K,UI)  = + FM_SHIFTED_DIFFUSE * SM
               ENDDO
              ENDDO
             ENDDO
            ENDIF

!  Step 4C. Surface Linearization Source vector for user streams
!   Update version 2.5, 10/12/15. Introduce loop over surface WFS

            IF ( DO_SURFACE_JACOBIANS .and. DO_USER_STREAMS  ) THEN
             DO Q = 1, N_SURFACEWFS
              DO UI = 1, N_USER_STREAMS
               SP = ZERO
               SM = ZERO
               DO L = FOURIER, NRRS_MOMENTS
                SP = SP + LS_HOLD_RRS(Q,L) * PLM_PXUI(UI,L)
                SM = SM + LS_HOLD_RRS(Q,L) * PLM_MXUI(UI,L)
               ENDDO
               LS_QSOURCE_UPOS1(Q,UI) = + FM_SHIFTED_DIFFUSE * SP
               LS_QSOURCE_UNEG1(Q,UI) = + FM_SHIFTED_DIFFUSE * SM
              ENDDO
             ENDDO
            ENDIF

!  Whole layer Greens function Discrete ordinate solutions

            CALL GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER, TAYLOR_SMALL, & ! Inputs
              NSTREAMS, QUAD_WEIGHTS, DELTAUS,                    & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,               & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,     & ! Inputs
              QSOURCE_QPOS, QSOURCE_QNEG,                         & ! Inputs
              RAMAN_WUPPER, RAMAN_WLOWER, ATERM_SAVE,             & ! Outputs
              BTERM_SAVE, MULT_M, MULT_P )                          ! Outputs

!  Whole layer Linearized Greens function solutions

            IF ( DO_ATMOS_JACOBIANS ) THEN
              CALL LPC_GREENFUNC_SOLUTION_1 &
            ( M, N, Z, DO_MULTIPLIER, TAYLOR_ORDER,           & ! Inputs
              TAYLOR_SMALL, KS2(N), KF2(N), KPARS,            & ! Inputs
              NKSTORAGE, NKSTORAGE2, NSTREAMS, QUAD_WEIGHTS,  & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,           & ! Inputs
              ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P,         & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS, & ! Inputs
              DELTAUS, QSOURCE_QPOS, QSOURCE_QNEG,             & ! Inputs
              L_XPOS, L_KEIGEN, L_KTRANS, L_EIGENNORM_SAVE,   & ! Inputs
              L_ASOURCE, L_DELTRANS, L_INITRANS,              & ! Inputs
              L_DELTAUS, L_QSOURCE_QPOS, L_QSOURCE_QNEG,       & ! Inputs
              L_RAMAN_WUPPER, L_RAMAN_WLOWER,                 & ! Outputs
              L_ATERM_SAVE, L_BTERM_SAVE )                      ! Outputs
            ENDIF

!  Whole layer Surface Linearized Greens function solutions

            IF ( DO_SURFACE_JACOBIANS ) THEN
              CALL LS_GREENFUNC_SOLUTION_1 &
                ( M, N, ZS, NSTREAMS, N_SURFACEWFS, QUAD_WEIGHTS,         & ! Inputs
                  FLIPPER, PICUTOFF, INITRANS, XPOS, EIGENNORM_SAVE,      & ! Inputs
                  LS_QSOURCE_QPOS, LS_QSOURCE_QNEG, MULT_M, MULT_P,       & ! Inputs
                  LS_RAMAN_WUPPER, LS_RAMAN_WLOWER, LS_ATERM_SAVE, LS_BTERM_SAVE )    ! Outputs
            ENDIF

!  Debug
!      write(68,*)z, j
!      do n = 1, 25
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WLOWER(i+8,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(RAMAN_WUPPER(i+8,n),i=1,8)
!      enddo
!      pause

!  Call to the post processing routine

            CALL POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                           & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                    & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                 & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),         & ! Inputs
            TAYLOR_ORDER, TAYLOR_SMALL, TRANS_USERM, DELTAUS,     & ! Inputs
            USER_STREAMS, FLIPPER, CHANGE_EXPONENT, Z,            & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1, & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, U_XPOS, U_XNEG,             & ! Inputs (No KTRANS)
            EMULT_UP, WLAYER_PIST_UP, WSU_SUSAV, WSU_SDSAV,       & ! Outputs
            EMULT_DN, WLAYER_PIST_DN, WSD_SUSAV, WSD_SDSAV )        ! Outputs

!  Call to the Linearized post processing routine

            IF ( DO_ATMOS_JACOBIANS ) THEN
              CALL LPC_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                                  & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL,                           & ! Inputs
            FOURIER, N, NSTREAMS, N_USER_STREAMS,                        & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N),                & ! Inputs
            KS2(N), KF2(N), KPARS, NKSTORAGE2, TAYLOR_ORDER,             & ! Inputs
            TAYLOR_SMALL, TRANS_USERM, DELTAUS, USER_STREAMS,            & ! Inputs
            FLIPPER, CHANGE_EXPONENT, Z,                                 & ! Inputs
            PICUTOFF, ASOURCE, DELTRANS, INITRANS,                       & ! Inputs
            ATERM_SAVE, BTERM_SAVE, QSOURCE_UPOS1, QSOURCE_UNEG1,        & ! Inputs
            HMULT_1, HMULT_2, KEIGEN, KTRANS, U_XPOS, U_XNEG,            & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,                              & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,                              & ! Inputs
            L_TRANS_USERM, L_DELTAUS, L_ASOURCE, L_DELTRANS, L_INITRANS, & ! Inputs
            L_QSOURCE_UPOS1, L_QSOURCE_UNEG1, L_ATERM_SAVE,L_BTERM_SAVE, & ! Inputs
            L_HMULT_1, L_HMULT_2, L_KEIGEN, L_U_XPOS, L_U_XNEG,          & ! Inputs
            L_WLAYER_PIST_UP, L_WLAYER_PIST_DN )                           ! Outputs
            ENDIF

!  Call to the Surface-Linearized post processing routine

            IF ( DO_SURFACE_JACOBIANS ) THEN
              CALL LS_POSTPROCESSING_1 &
          ( DO_UPWELLING, DO_DNWELLING,                   & ! Inputs
            DO_MSMODE_LRRS, DO_SSCORR_OVERALL, FOURIER,   & ! Inputs
            N, NSTREAMS, N_USER_STREAMS, N_SURFACEWFS,    & ! Inputs
            STERM_LAYERMASK_UP(N), STERM_LAYERMASK_DN(N), & ! Inputs
            ZS, INITRANS, U_XPOS, U_XNEG,                 & ! Inputs
            EMULT_UP, WSU_SUSAV, WSU_SDSAV,               & ! Inputs
            EMULT_DN, WSD_SUSAV, WSD_SDSAV,               & ! Inputs
            LS_ATERM_SAVE, LS_BTERM_SAVE,                 & ! Inputs
            LS_QSOURCE_UPOS1, LS_QSOURCE_UNEG1,           & ! Inputs
            LS_WLAYER_PIST_UP, LS_WLAYER_PIST_DN )          ! Outputs
            ENDIF

!  End stream loop  (mick note: this stream loop starts @ ~line 2870)

           ENDDO

!  End RRS_CONTROL clause

          ENDIF

!  Continuation point for avoiding calculation (monochromatic only)

 5656     CONTINUE

!  End loop over RRS transitions (binned or mono)  (mick note: this transitions loop starts @ ~line 2335)

        ENDDO

!  Elastic-only continuation point

 5222   continue

!  End layer loop  (mick note: this layer loop starts @ ~line 880!)

      ENDDO

!  Debug, 3 June 2009 (profiles). 8 June 2009 (column/surface)
!   Used frequently for debuggin. Version 2.5

      if ( do_debug_sources.and.fourier.eq.0 ) then
 5     format(2i4,1p2e25.16)
       do n = 1, 5
        write(51,5)1,n,WLAYER_PIST_UP(N,1),WLAYER_PIST_DN(N,1)
        do i = 1, 3
          write(41,5)i,n,RAMAN_WUPPER(i,n),RAMAN_WLOWER(i,n)
        enddo
       enddo
       if ( do_atmos_jacobians) then
        q = 1
        k = 1
        do n = 1, 5
         nk = nkstorage2(n,k)
         write(53,5)1,n,L_WLAYER_PIST_UP(q,NK,1),L_WLAYER_PIST_DN(q,NK,1)
         do i = 1, 3
          write(43,5)i,n,L_RAMAN_WUPPER(q,i,nk),L_RAMAN_WLOWER(q,i,nk)
         enddo
        enddo
       endif
       if ( do_surface_jacobians) then
        do n = 1, 8
         write(53,5)1,n,LS_WLAYER_PIST_UP(1,N,1),LS_WLAYER_PIST_DN(1,N,1)
         do i = 1, 3
          write(43,5)i,n,LS_RAMAN_WUPPER(1,i,n),LS_RAMAN_WLOWER(1,i,n)
         enddo
        enddo
       endif
!       pause '4143/5153'
      endif

!  End of module

      RETURN
      END SUBROUTINE L_SOURCES_MASTER_1

!  End module

      END MODULE lrrs_L_sources_master_1_m
