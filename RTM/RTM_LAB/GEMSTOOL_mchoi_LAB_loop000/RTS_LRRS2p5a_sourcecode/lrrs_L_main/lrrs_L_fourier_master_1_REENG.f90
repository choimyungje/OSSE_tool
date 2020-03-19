
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

! ###############################################################
! #                                                             #
! #   SUBROUTINES :                                             #
! #             L_RAMAN_FOURIER_1 (*)                           #
! #                                                             #
! #   Note (*): Version 2.1 with BRDF setup, November 2007      #
! #   Note    : Version 2.2 with Linearization, July 2009       #
! #   Note    : Version 2.3 in F90,              March 2011     #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Use of Supplement-derived BRDF/SLEAVE inputs and control
!    (2) Bookkeeping improvements (use of "Only", clearer I/O specifications)
!    (3) Use of new Taylor series inputs and modules.

      MODULE lrrs_L_fourier_master_1_m

!      USE LRRS_PARS_m, Only : SDU,LDU
!   -- Rob mod 5/12/17 for 2p5a, SSFULL renamed, SSCORR_NADIR added

      PRIVATE
      PUBLIC :: L_RAMAN_FOURIER_1

      CONTAINS

!Rob mod -- 10/10/18. Add OMP Inputs and Outputs

      SUBROUTINE L_RAMAN_FOURIER_1 ( DO_TIMING, NTHREADS,                     & ! OMP Inputs
         DO_RRS_OVERALL, DO_DIRECTRRS_ONLY, DO_MSMODE_LRRS, PROGRESS,         & ! Inputs
         DO_ELASTIC_ONLY, DO_SSCORR_GENERAL, DO_SSCORR_ALONE,                 & ! Inputs
         DO_BIN_REALIZATION, DO_ENERGY_BALANCING, DO_PLANE_PARALLEL,          & ! Inputs
         DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_DIRECT_BEAM,                  & ! Inputs
         DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_BRDF_SURFACE,        & ! Inputs
         DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, do_BRDF_Wav1, do_SLEAVE_Wav1,   & ! Inputs
         DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,        & ! Inputs (Linearized)
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS,               & ! Inputs (Linearized)
         N_SURFACE_WFS, N_SLEAVE_WFS,                                         & ! Inputs (Linearized)
         NPOINTS_OUTER, OFFSET_INNER, NPOINTS_INNER,                          & ! Inputs
         N_RRSBINS, BINMAP, NPOINTS_MONO, W_EXCIT, FLUXES_RANKED,             & ! Inputs
         ALBEDOS_RANKED, BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,        & ! Inputs
         SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                       & ! Inputs
         NLAYERS, NSTREAMS, NMOMENTS, N_LOUTPUT, N_USER_STREAMS,              & ! Inputs
         FOURIER, COS_SZA, TAYLOR_SMALL, TAYLOR_ORDER, FLUX_FACTOR,           & ! Inputs
         QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS,              & ! Inputs
         STERM_MASK_UP, LEVELMASK_UP, STERM_MASK_DN, LEVELMASK_DN,            & ! Inputs
         DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, OMEGAMOMS_CABANNES,            & ! Inputs
         OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN, OMEGAMOMS_RRSGAIN,              & ! Inputs
         BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT,                           & ! Inputs
         BEAM_DTRANS, BEAM_ETRANS, SAVE_TRANS_USERM,                          & ! Inputs
         LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,            & ! Inputs (linearized)
         LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,        & ! Inputs (linearized)
         L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC, L_OMEGAMOMS_CABANNES,      & ! Inputs (linearized)
         L_OMEGAMOMS_RRSLOSS, L_OMEGAMOMS_RRSBIN,  L_OMEGAMOMS_RRSGAIN,       & ! Inputs (linearized)
         NKSTORAGE, NKSTORAGE2, L_BEAM_ITRANS, L_BEAM_AVSECANT,               & ! Inputs (linearized)
         L_BEAM_ETRANS, L_BEAM_DTRANS, L_SAVE_TRANS_USERM,                    & ! Inputs (linearized)
         ELASTIC_F_UP, MEAN_ELASTIC_UP, FLUX_ELASTIC_UP,                      & ! outputs
         ELASTIC_F_DN, MEAN_ELASTIC_DN, FLUX_ELASTIC_DN,                      & ! outputs
         RAMAN_F_UP, MEAN_RAMAN_UP, FLUX_RAMAN_UP,                            & ! outputs
         RAMAN_F_DN, MEAN_RAMAN_DN, FLUX_RAMAN_DN,                            & ! outputs
         L_ELASTIC_F_UP, L_ELASTIC_F_DN, LS_ELASTIC_F_UP, LS_ELASTIC_F_DN,               & ! outputs (linearized)
         L_MEAN_ELASTIC_UP, L_FLUX_ELASTIC_UP, L_MEAN_ELASTIC_DN, L_FLUX_ELASTIC_DN,     & ! outputs (linearized)
         LS_MEAN_ELASTIC_UP, LS_FLUX_ELASTIC_UP, LS_MEAN_ELASTIC_DN, LS_FLUX_ELASTIC_DN, & ! outputs (linearized)
         L_RAMAN_F_UP, L_RAMAN_F_DN, LS_RAMAN_F_UP, LS_RAMAN_F_DN,                       & ! outputs (linearized)
         L_MEAN_RAMAN_UP, L_FLUX_RAMAN_UP, L_MEAN_RAMAN_DN, L_FLUX_RAMAN_DN,             & ! outputs (linearized)
         LS_MEAN_RAMAN_UP, LS_FLUX_RAMAN_UP, LS_MEAN_RAMAN_DN, LS_FLUX_RAMAN_DN,         & ! outputs (linearized)
         FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE,                                        & ! Outputs (Exceptions)
         OMP_Elastic_Time, OMP_Raman_Time )                                                ! OMP Timing Outputs

!  Module of dimensions and numbers
!   -- Rob mod 5/12/17 for 2p5a, introduced LRRS_MISCSETUPS_1_m, reorganized LRRS_RTSOLUTIONS_1_m

      USE LRRS_PARS_m

!  Use Modules

      USE LRRS_RTSOLUTIONS_1_m        , Only : HOMOG_SOLUTION, UHOMOG_SOLUTION
      USE LRRS_MISCSETUPS_1_m         , Only : DBEAM_SETUP, UDBEAM_SETUP, LEGENDRE_SETUPS
      USE LRRS_BVPROBLEM_m            , Only : SURFACE_HOMOG_SETUP, SURFACE_BEAM_SETUP,  &
                                               BVPMATRIX_INIT, BVPMATRIX_SETUP, BVPCOLUMN_SETUP
      USE LRRS_AUX2_m                 , Only : DGBTRF, DGBTRS, DGETRS, DGETRF
      USE LRRS_POSTPROCESSING_1_m     , Only : HOMOGMULT_1, BEAMMULT_UP_1, BEAMMULT_DN_1, BOASOURCE
      USE LRRS_RAMAN_INTENSITY_1_m    , Only : MIFLUX_INTENSITY_1, RAMAN_INTENSITY_UP_1,  RAMAN_INTENSITY_DN_1

      USE LRRS_L_RTSOLUTIONS_1_m      , Only : LPC_HOMOG_SOLUTION, LPC_UHOMOG_SOLUTION
      USE LRRS_L_BVPROBLEM_m          , Only : L_SURFACE_HOMOG_SETUP, L_SURFACE_BEAM_SETUP, LS_SURFACE_BEAM_SETUP, &
                                               LC_BVPCOLUMN_SETUP, LPR_BVPCOLUMN_SETUP, LSR_BVPCOLUMN_SETUP
      USE LRRS_L_POSTPROCESSING_1_m   , Only : LPC_HOMOGMULT_1, LPC_BEAMMULT_UP_1, LPC_BEAMMULT_DN_1, &
                                               LPC_TOASOURCE, LPC_BOASOURCE, LS_TOASOURCE, LS_BOASOURCE
      USE LRRS_L_SOURCES_MASTER_1_m   , Only : L_SOURCES_MASTER_1
      USE LRRS_L_ELASTIC_MASTER_1_m   , Only : L_ELASTIC_MASTER_1
      USE LRRS_L_RAMAN_JACOBIANS_1_m  , Only : MIFLUX_JACOBIAN_1, RAMAN_JACOBIAN_UP_1, RAMAN_JACOBIAN_DN_1, &
                                               MIFLUX_SURFJAC_1,  RAMAN_SURFJAC_UP_1,  RAMAN_SURFJAC_DN_1

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Added 10/10/18. Number of OMP Cores to be used (1, 2, 4, 8...)

      INTEGER  , intent(in) ::  NTHREADS

!  Added 10/10/18. Monitoring flag

      logical  , intent(in) :: DO_TIMING

!  Progress number
!    If this is zero, no output to screen

      INTEGER  , INTENT(IN) :: PROGRESS

!  elastic control flag

      LOGICAL  , INTENT(IN) :: DO_ELASTIC_ONLY

!  Raman  flags for computing filling

      LOGICAL  , INTENT(IN) :: DO_RRS_OVERALL
      LOGICAL  , INTENT(IN) :: DO_DIRECTRRS_ONLY

!  Multiple scatter mode

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LRRS

!  Overall SS correction flag; Now incorporates the DB term
!   -- Rob mod 5/12/17 for 2p5a, SSFULL renamed, SSCORR_GENERAL flag used

      LOGICAL  , INTENT(IN) :: DO_SSCORR_GENERAL
      LOGICAL  , INTENT(IN) :: DO_SSCORR_ALONE

!  Upwelling and downwelling flags

      LOGICAL  , INTENT(IN) :: DO_UPWELLING
      LOGICAL  , INTENT(IN) :: DO_DNWELLING

!  User stream flag

      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS

!  Binning realization

      LOGICAL  , INTENT(IN) :: DO_BIN_REALIZATION

!  Energy balancing flag

      LOGICAL  , INTENT(IN) :: DO_ENERGY_BALANCING

!  Plane parallel

      LOGICAL  , INTENT(IN) ::  DO_PLANE_PARALLEL

!  Mean value output

      LOGICAL  , INTENT(IN) :: DO_MVOUT_ONLY
      LOGICAL  , INTENT(IN) :: DO_ADDITIONAL_MVOUT

!  Direct beam

      LOGICAL  , INTENT(IN) :: DO_DIRECT_BEAM

!  Surface control
!  ---------------

!  New Version 2.5, general BRDF/SLEAVE input

      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE
      LOGICAL  , intent(in) :: DO_SURFACE_LEAVING
      LOGICAL  , intent(in) :: DO_SL_ISOTROPIC

!  Flags for the multiple-point surface options

      LOGICAL  , intent(in) :: DO_BRDF_Wav1
      LOGICAL  , intent(in) :: DO_SLEAVE_Wav1

!  Lambertian albedos. New Version 2.5, 9/11/15

      REAL(fpk), intent(in) :: ALBEDOS_RANKED ( MAX_POINTS )

!  Fourier components of BRDF, in the following order
!    incident solar direction,    reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar direction,    reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), intent(in) :: BRDF_F_0      ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: BRDF_F        ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: USER_BRDF_F_0 ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: USER_BRDF_F   ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

!  FOLLOWING CODE IS COMMENTED OUT, VERSION 2.5.
!  overall lambertian surface
!      LOGICAL  , INTENT(IN) :: DO_LAMBERTIAN_SURFACE
!  Surface leaving prototype for Version 2.4, Now commented out, 9/11/15
!     Add SLTERM flag     . @@@ Rob Fix 06 Sep 12.
!     New Isotropic SLTERM. @@@ Rob Fix 06 sep 12
!      LOGICAL  , INTENT(IN)  :: DO_LRRS_SLTERM
!      REAL(FPK), INTENT(IN)  :: LRRS_SLTERM (MAX_POINTS)
!  Surface contributions for the Version 2.3 code. Now commented out, 9/11/15
!      REAL(FPK), INTENT(IN) :: FL1 ( MAX_POINTS )
!      REAL(FPK), INTENT(IN) :: FL2 ( MAX_POINTS )
!  BRDF control and quadrature for the Version 2.3 code. Now commented out, 9/11/15
!      INTEGER  , INTENT(IN) :: N_BRDF
!      REAL(FPK), INTENT(IN) :: X_BRDF ( MAX_STREAMS_BRDF )
!      REAL(FPK), INTENT(IN) :: A_BRDF ( MAX_STREAMS_BRDF )
!  BRDF Inputs for the Version 2.3 code. Now commented out, 9/11/15
!      REAL(FPK), INTENT(IN) :: BRDFUNC ( MAX_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )
!      REAL(FPK), INTENT(IN) :: BRDFUNC_0  ( MAX_STREAMS, MAX_STREAMS_BRDF )
!      REAL(FPK), INTENT(IN) :: USER_BRDFUNC    ( MAX_USER_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )
!      REAL(FPK), INTENT(IN) :: USER_BRDFUNC_0  ( MAX_USER_STREAMS, MAX_STREAMS_BRDF )

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), intent(in) :: SLTERM_ISOTROPIC ( MAX_POINTS )

!  Fourier components of Surface-leaving terms:
!    solar direction, SL-transmitted quadrature streams
!    solar direction, SL-transmitted user streams

      REAL(fpk), intent(in) ::  SLTERM_F_0      ( 0:MAX_MOMENTS, MAX_STREAMS,      MAX_POINTS )
      REAL(fpk), intent(in) ::  USER_SLTERM_F_0 ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )

!  Linearization control
!  ---------------------

!  Profile/column/surface linearization control inputs

      LOGICAL  , INTENT(IN) :: do_profile_wfs
      LOGICAL  , INTENT(IN) :: do_column_wfs
      LOGICAL  , INTENT(IN) :: do_surface_wfs
      LOGICAL  , INTENT(IN) :: do_sleave_wfs

!  Linearization control

      LOGICAL  , INTENT(IN) :: layer_vary_flag   (max_layers)
      INTEGER  , INTENT(IN) :: layer_vary_number (max_layers)
      INTEGER  , INTENT(IN) :: n_totalcolumn_wfs
      INTEGER  , INTENT(IN) :: n_surface_wfs
      INTEGER  , INTENT(IN) :: n_sleave_wfs

!  Numbers
!  -------

!  number of streams

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: N_USER_STREAMS

!  Number of moments

      INTEGER  , INTENT(IN) :: NMOMENTS

!  number of layers

      INTEGER  , INTENT(IN) :: NLAYERS

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)

      INTEGER  , INTENT(IN) :: N_LOUTPUT

!  Fourier number

      INTEGER  , INTENT(IN) :: FOURIER

!  Solar beam input Cosine

      REAL(FPK), INTENT(IN) :: COS_SZA

!  Small number control. Extended, Version 2.5, 9/11/15

      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL
      INTEGER  , INTENT(IN) :: TAYLOR_ORDER

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

      INTEGER  , INTENT(IN) :: NPOINTS_OUTER
      INTEGER  , INTENT(IN) :: NPOINTS_INNER
      INTEGER  , INTENT(IN) :: OFFSET_INNER
      INTEGER  , INTENT(IN) :: NPOINTS_MONO

      INTEGER  , INTENT(INOUT) :: W_EXCIT

!  fluxes

      REAL(FPK), INTENT(IN) :: FLUXES_RANKED  ( MAX_POINTS )

!  Layer masks for doing integrated source terms
!  ---------------------------------------------

      LOGICAL  , INTENT(IN) :: STERM_MASK_UP ( MAX_LAYERS )
      LOGICAL  , INTENT(IN) :: STERM_MASK_DN ( MAX_LAYERS )

      INTEGER  , INTENT(IN) :: LEVELMASK_UP    (MAX_LOUTPUT)
      INTEGER  , INTENT(IN) :: LEVELMASK_DN    (MAX_LOUTPUT)

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

      INTEGER  , INTENT(IN) :: N_RRSBINS  ( MAX_POINTS )
      INTEGER  , INTENT(IN) :: BINMAP ( MAX_BINS, MAX_POINTS )

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

!  Beam Attenuations
!  -----------------

!  Solar beam transmittances, average secant factors

      INTEGER  , INTENT(IN) :: BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_ETRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(IN) :: SAVE_TRANS_USERM &
        (  MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Surface linearizations (BRDF)

      REAL(fpk), intent(in) :: LS_BRDF_F_0      ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LS_BRDF_F        ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Surface-leaving linearizations

      REAL(fpk), intent(in) :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAX_POINTS )
      REAL(fpk), intent(in) :: LSSL_SLTERM_F_0        ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LSSL_USER_SLTERM_F_0   ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )

!  Scaled elastic-scattering optical properties, LINEARIZATIONS
!  ------------------------------------------------------------

!  Vertical optical thickness values

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_ELASTIC ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived cabannes scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_CABANNES ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  RRS properties, Linearizations

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSLOSS ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSBIN  ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSGAIN ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Other linearizations
!  --------------------

!  Linearizations of Solar beam transmittances, average secant factors

      REAL(FPK), INTENT(IN) :: L_BEAM_ITRANS   ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_BEAM_AVSECANT ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_BEAM_ETRANS   ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_BEAM_DTRANS   ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(IN) :: L_SAVE_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Linearization bookkeeping storage

      INTEGER  , INTENT(IN) :: NKSTORAGE   ( MAX_LAYERS, 0:MAX_LAYERS )
      INTEGER  , INTENT(IN) :: NKSTORAGE2  ( MAX_LAYERS, 0:MAX_LAYERS )

!  OUTPUT ARGUMENTS START HERE
!  ===========================

!  Elastic output at User angles, for one Fourier component
!  --------------------------------------------------------

!  Fourier output

      REAL(FPK), INTENT(OUT) :: ELASTIC_F_UP ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: ELASTIC_F_DN ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  Mean-value output

      REAL(FPK), INTENT(INOUT) :: MEAN_ELASTIC_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: MEAN_ELASTIC_DN ( MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: FLUX_ELASTIC_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: FLUX_ELASTIC_DN ( MAX_LOUTPUT, MAX_POINTS )

!  Fourier output, atmospheric weighting functions

      REAL(FPK), INTENT(OUT) :: L_ELASTIC_F_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: L_ELASTIC_F_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  Fourier output, surface weighting functions

      REAL(FPK), INTENT(OUT) :: LS_ELASTIC_F_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: LS_ELASTIC_F_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  Linearized Mean-value output
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: L_MEAN_ELASTIC_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: L_MEAN_ELASTIC_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: L_FLUX_ELASTIC_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: L_FLUX_ELASTIC_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: LS_MEAN_ELASTIC_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: LS_MEAN_ELASTIC_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: LS_FLUX_ELASTIC_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: LS_FLUX_ELASTIC_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )

!  Raman output at User angles, for one Fourier component
!  ------------------------------------------------------

!  Fourier output

      REAL(FPK), INTENT(OUT) :: RAMAN_F_UP ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: RAMAN_F_DN ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  Regular Mean-value output
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: MEAN_RAMAN_UP( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: MEAN_RAMAN_DN( MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: FLUX_RAMAN_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: FLUX_RAMAN_DN ( MAX_LOUTPUT, MAX_POINTS )

! Fourier output  atmospheric weighting functions, post-processed values

      REAL(FPK), INTENT(OUT) :: L_RAMAN_F_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_RAMAN_F_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  For the surface weighting functions, post-processed values

      REAL(FPK), INTENT(OUT) :: LS_RAMAN_F_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: LS_RAMAN_F_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  Linearized Mean-value output
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: L_MEAN_RAMAN_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: L_MEAN_RAMAN_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: L_FLUX_RAMAN_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: L_FLUX_RAMAN_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: LS_MEAN_RAMAN_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: LS_MEAN_RAMAN_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: LS_FLUX_RAMAN_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: LS_FLUX_RAMAN_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )

!  module status:
!   Message_sub comes from any of the subroutines called here.
!   Point_trace gives a trace of subroutine where message occurred

      LOGICAL            , INTENT(OUT) :: FAIL
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE_SUB
      CHARACTER (LEN=120), INTENT(OUT) :: POINT_TRACE

! Rob mod 10/10/18. Add OMP timing outputs

      REAL, INTENT(INOUT) :: OMP_Elastic_Time, OMP_Raman_Time

!  LOCAL VARIABLES START HERE
!  ==========================

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

!  Eigenvalues and transmittances

      REAL(FPK) :: L0_KEIGEN  ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L0_KTRANS  ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Particular solution vectors

      REAL(FPK) :: L0_WPARTIC ( MAX_2_STREAMS, MAX_LAYERS, MAX_POINTS )

!  homogeneous solution vectors and integration constants

      REAL(FPK) :: L0_WHOM_XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L0_WHOM_LCON ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L0_WHOM_MCON ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Linearized discrete ordinate elastic field solutions
!  ----------------------------------------------------

!  Eigenvalues and transmittancaes

      REAL(FPK) :: L_L0_KEIGEN ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L_L0_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Particular solution vectors

      REAL(FPK) :: L_L0_WPARTIC ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS_NK, MAX_POINTS )

!  homogeneous solution vectors and integration constants

      REAL(FPK) :: L_L0_WHOM_XPOS ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L_L0_WHOM_LCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS_SQ, MAX_POINTS )
      REAL(FPK) :: L_L0_WHOM_MCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS_SQ, MAX_POINTS )

!  Surface linearization of integration constants

      REAL(FPK) :: LS_L0_WHOM_LCON ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: LS_L0_WHOM_MCON ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  output Direct beam module
!  -------------------------

!  Reflected Direct beam at surface (quadrature + User streams)

      REAL(FPK) :: DIRECT_BEAM      ( MAX_STREAMS )
      REAL(FPK) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS )

!  solar beam optical thickness and atmospheric attenuation

!      REAL(FPK) :: SOLAR_BEAM_OPDEP
      REAL(FPK) :: ATMOS_ATTN

!  Reflected Direct beam at surface (quadrature + User streams)

      REAL(FPK) :: L_DIRECT_BEAM      ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_USER_DIRECT_BEAM ( MAX_ATMOSWFS,MAX_USER_STREAMS )

!  Reflected Direct beam at surface (quadrature + User streams)

      REAL(FPK) :: LS_DIRECT_BEAM      ( MAX_SURFACEWFS, MAX_STREAMS )
      REAL(FPK) :: LS_USER_DIRECT_BEAM ( MAX_SURFACEWFS, MAX_USER_STREAMS )

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

!  LINEARIZED output arguments from Homogeneous solution modules
!  -------------------------------------------------------------

!  (Positive) Eigenvalues

      REAL(FPK) :: L_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!   Whole layer transmittance factors for +/- eigenvalues

      REAL(FPK) :: L_KEIGEN ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  Eigenvector solutions

      REAL(FPK) :: L_XPOS ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: L_XNEG ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  LInearized eigenmatrices

      REAL(FPK) :: L_DAB_SAVE ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: L_SAB_SAVE ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: L_EIGENMAT_SAVE ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )

!  Eigennorm saved

      REAL(FPK) :: L_EIGENNORM_SAVE ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  Reflected homogeneous solutions at ground

      REAL(FPK) :: L_R2_HOMP(MAX_ATMOSWFS,MAX_STREAMS,MAX_STREAMS)
      REAL(FPK) :: L_R2_HOMM(MAX_ATMOSWFS,MAX_STREAMS,MAX_STREAMS)

!  user-stream homogeneous solution vectors

      REAL(FPK) :: U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  user-stream homogeneous solution vectors

      REAL(FPK) :: L_U_XPOS ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: L_U_XNEG ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Multipliers
!  -----------

!  Whole layer multipliers (homogeneous)

      REAL(FPK) :: HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Zeta functions (homogeneous) = 1 / [ 1/mu) +/- k_a ]

      REAL(FPK) :: ZETA_P ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: ZETA_M ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Whole layer multipliers (homogeneous)

      REAL(FPK) :: L_HMULT_1 ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: L_HMULT_2 ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  BV problem variables
!  --------------------

!  Compressed band matrix, as in A

      REAL(FPK) :: BANDMAT2(MAX_BANDTOTAL,MAX_TOTAL)
      REAL(FPK) :: SMAT2   (MAX_2_STREAMS, MAX_2_STREAMS)

!  compression indexing for the band matrix

      INTEGER   :: BMAT_ROWMASK(MAX_TOTAL,MAX_TOTAL)

!  LU decomposition PIVOT

      INTEGER   :: IPIVOT ( MAX_TOTAL )
      INTEGER   :: SIPIVOT (MAX_2_STREAMS)

!  Column vectors for Radiances

      REAL(FPK) :: COL2  ( MAX_TOTAL, 1 )
      REAL(FPK) :: SCOL2 ( MAX_2_STREAMS, 1 )

!  Integration constants

      REAL(FPK) :: LCON ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: MCON ( MAX_STREAMS, MAX_LAYERS )

!  Integration constants x homogeneous solution vectors

      REAL(FPK) :: LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Column vectors for Atmospheric weighting functions

      REAL(FPK) :: COL2WF ( MAX_TOTAL, MAX_ATMOSWFS )
      REAL(FPK) :: SCOL2WF ( MAX_2_STREAMS, MAX_ATMOSWFS )

!  Linearized Integration constants (Atmospheric linearization)

      REAL(FPK) :: NCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: PCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  Column vectors for Surface weighting functions

      REAL(FPK) :: COL2WF_SURF  ( MAX_TOTAL, MAX_SURFACEWFS )
      REAL(FPK) :: SCOL2WF_SURF ( MAX_2_STREAMS, MAX_SURFACEWFS )

!  Linearized Integration constants (Surface linearization)

      REAL(FPK) :: NCON_SURF ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: PCON_SURF ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS )

!  Reflected beam solutions at ground

      REAL(FPK) :: R2_BEAM(MAX_STREAMS)

!  output arguments from Beam solution modules
!  -------------------------------------------

!  Particular solutions at upper and lower boundaries
!   (as evaluated for the inelastic field by Green's function method)

!      REAL(FPK) :: RAMAN_WUPPER  ( MAX_LAYERS, MAX_2_STREAMS )
!      REAL(FPK) :: RAMAN_WLOWER  ( MAX_LAYERS, MAX_2_STREAMS )

      REAL(FPK) :: RAMAN_WUPPER  ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK) :: RAMAN_WLOWER  ( MAX_2_STREAMS, MAX_LAYERS )

!  Particular solutions at upper and lower boundaries
!   (as evaluated for the inelastic field by Green's function method)

      REAL(FPK) :: L_RAMAN_WUPPER ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS_SQ )
      REAL(FPK) :: L_RAMAN_WLOWER ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS_SQ )

!  Reflected beam solutions at ground

      REAL(FPK) :: L_R2_BEAM(MAX_ATMOSWFS,MAX_STREAMS)

!  Particular solutions at upper and lower boundaries
!   (as evaluated for the inelastic field by Green's function method)

      REAL(FPK) :: LS_RAMAN_WUPPER ( MAX_SURFACEWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK) :: LS_RAMAN_WLOWER ( MAX_SURFACEWFS, MAX_2_STREAMS, MAX_LAYERS )

!  diffuse beam solutions at ground (includes reflected part)
!    This term is not present in the elastic solution

      REAL(FPK) :: LS_DIFFUSE_BEAM ( MAX_SURFACEWFS, MAX_STREAMS )

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

!  Atmospheric linearizations

      REAL(FPK) :: L_WLAYER_HOST_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK) :: L_WLAYER_HOST_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_USER_STREAMS )

      REAL(FPK) :: L_WLAYER_PIST_UP ( MAX_ATMOSWFS, MAX_LAYERS_SQ, MAX_USER_STREAMS )
      REAL(FPK) :: L_WLAYER_PIST_DN ( MAX_ATMOSWFS, MAX_LAYERS_SQ, MAX_USER_STREAMS )

      REAL(FPK) :: L_TOA_SOURCE        (MAX_ATMOSWFS,MAX_USER_STREAMS)
      REAL(FPK) :: L_BOA_SOURCE        (MAX_ATMOSWFS,MAX_USER_STREAMS)
      REAL(FPK) :: L_DIRECT_BOA_SOURCE (MAX_ATMOSWFS,MAX_USER_STREAMS)

!  Surface linearizations

      REAL(FPK) :: LS_WLAYER_HOST_UP ( MAX_SURFACEWFS, MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK) :: LS_WLAYER_HOST_DN ( MAX_SURFACEWFS, MAX_LAYERS, MAX_USER_STREAMS )

      REAL(FPK) :: LS_WLAYER_PIST_UP ( MAX_SURFACEWFS, MAX_LAYERS, MAX_USER_STREAMS )
      REAL(FPK) :: LS_WLAYER_PIST_DN ( MAX_SURFACEWFS, MAX_LAYERS, MAX_USER_STREAMS )

      REAL(FPK) :: LS_TOA_SOURCE        ( MAX_SURFACEWFS, MAX_USER_STREAMS )
      REAL(FPK) :: LS_BOA_SOURCE        ( MAX_SURFACEWFS, MAX_USER_STREAMS )
      REAL(FPK) :: LS_DIRECT_BOA_SOURCE ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Miscellaneous
!  -------------

!  Help arrays for linearization

      REAL(FPK) :: U_HELP_P ( MAX_STREAMS, 0:MAX_MOMENTS, MAX_LAYERS )
      REAL(FPK) :: U_HELP_M ( MAX_STREAMS, 0:MAX_MOMENTS, MAX_LAYERS )

      REAL(FPK) :: L_WLOWER ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK) :: L_WUPPER ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS )

!  local linearized optical depths and user-stream transmittances

      REAL(FPK) :: LOCAL_DELTAUS ( MAX_LAYERS )
      REAL(FPK) :: LOCAL_TRANS_USERM(MAX_USER_STREAMS,MAX_LAYERS)

      REAL(FPK) :: L_LOCAL_DELTAUS     ( MAX_ATMOSWFS, MAX_LAYERS )
      REAL(FPK) :: L_LOCAL_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS )

!  Local reflected field (not used here) is the diffuse BOA source term

      REAL(FPK) :: IDOWNSURF ( MAX_STREAMS )

!  Fourier output (quadrature streams)
!    Use as debug, except for Flux calculations

      REAL(FPK) :: QUADRAMAN_F_UP &
              ( MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

      REAL(FPK) :: QUADRAMAN_F_DN &
              ( MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

!  For the atmospheric weighting functions, discrete ordinate values

      REAL(FPK) :: L_QUADRAMAN_F_UP &
       ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

      REAL(FPK) :: L_QUADRAMAN_F_DN &
       ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

!  For the surface weighting functions, discrete ordinate values

      REAL(FPK) :: LS_QUADRAMAN_F_UP &
              ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

      REAL(FPK) :: LS_QUADRAMAN_F_DN &
              ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

!  Log-file output (LU = 6 gives screen output)

      LOGICAL  , PARAMETER :: LG = .False.
      INTEGER  , PARAMETER :: LU = 9

!  for exception handling

      CHARACTER (LEN=3) :: C3
      CHARACTER (LEN=2) :: C2

!  Local integers

      INTEGER   :: NTOTAL, N_SUPDIAG, N_SUBDIAG, N_GAIN_TERMS
      INTEGER   :: NSTR2, N_SOLUTIONS, PROG, M, LAYER, INFO
      INTEGER   :: I, I1, AA, N, C0, APT, CPT, L, UM, W, K, UI
      INTEGER   :: NPOINTS_TOTAL, NPOINTS_LOCAL
      INTEGER   :: Raman_IDX, Brdf_IDX, Sleave_IDX

!  Number of RRS solutions

      INTEGER :: N_LOSS_SOLUTIONS
      INTEGER :: N_GAIN_SOLUTIONS

!  Monochromatic realization

      LOGICAL :: DO_MONO_REALIZATION

!  Cabannes-Raman flag

      LOGICAL :: DO_CABANNES_RAMAN

!  Flags for including surface term

      LOGICAL :: DO_INCLUDE_SURFACE
      LOGICAL :: DO_REFLECTED_DIRECTBEAM

!  Flag for mean value output

      LOGICAL :: DO_INCLUDE_MVOUT

!  Flag for post-processing of elastic field

      LOGICAL :: DO_ELASTIC_POSTPROCESSING

!  Local flag for CB COPY

      LOGICAL :: DO_CBCOPY

!  Flag for overall single scatter condition. Commented out 10/10/18
!      LOGICAL :: DO_SSCORR_OVERALL

!  Overall atmospheric Jacobian flag

      LOGICAL :: DO_ATMOS_WFS

!  Surface factors and local XT, DT

      REAL(FPK) :: SURFACE_FACTOR, SURFACE_FACTOR_0
      REAL(FPK) :: DELTA_FACTOR, HOM, NORM, T1, T2
      REAL(FPK) :: H1, L_H1, LS_H1, H2, L_H2, LS_H2

!  solar irradiance flux

      REAL(FPK) :: SS_FLUX_MULTIPLIER
      REAL(FPK) :: FLUX_MULTIPLIER

!  Linearization control

      LOGICAL   :: DO_VARY, KVARY(0:MAX_LAYERS)
      INTEGER   :: Q, NK, KS(0:MAX_LAYERS), KF(0:MAX_LAYERS)
      INTEGER   :: KPARS(0:MAX_LAYERS), NPARS, WFF, WFS
      REAL(FPK) :: REFLEC, L_ATTN, L_AT

!  Local IOPs

      REAL(FPK) :: LOCAL_OMEGAMOMS   ( MAX_LAYERS, 0: MAX_MOMENTS )
      REAL(FPK) :: L_LOCAL_OMEGAMOMS ( MAX_ATMOSWFS, MAX_LAYERS, 0: MAX_MOMENTS )

!  Overall flag for computing filling

      LOGICAL   :: DO_FILLING

!  Rob Fix 26 March 2012. Add new variable DO_AVOID_UDBEAM

      LOGICAL   :: DO_AVOID_UDBEAM

!  Debug

      LOGICAL   :: ELASTIC_CALL=.FALSE.

!  OMP Variables. Added 10/10/18
!  -------------

!  OpenMP tests (general)

      INTEGER       :: OMP_MAXTHREADS
      INTEGER       :: TID, OMP_NTHREADS

!  OpenMP tests (timing)

      REAL          :: omp_e1, omp_e2

!  OpenMP functions

      INTEGER       :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

!  Local shared array, TID identification

      INTEGER       :: LRRS_Tid_SAVE (MAX_POINTS)

!  Local OMP-related variables

      INTEGER      :: W_EXCIT_0
      LOGICAL      :: FAIL_SUB

!  ###############
!  INITIAL SECTION
!  ###############

!  set status

      FAIL            = .FALSE.
      MESSAGE         = ' '
      POINT_TRACE     = ' '

!  Added 10/10/18. Set total number of threads (from NTHREADS input)
!     write(*,*) 'Setting total number of threads'

      OMP_MAXTHREADS = NTHREADS
      CALL OMP_SET_NUM_THREADS(OMP_MAXTHREADS)

!  set exclusion flags

      DO_MONO_REALIZATION = .not. DO_BIN_REALIZATION
      DO_CABANNES_RAMAN   = .not. DO_ENERGY_BALANCING

!  set overall filling flags

      DO_FILLING = .not.DO_ELASTIC_ONLY .and. DO_RRS_OVERALL

!  Overall Linearization flag

      DO_ATMOS_WFS = ( DO_COLUMN_WFS .OR. DO_PROFILE_WFS)

!  Fourier

      M = FOURIER

!mick fix 7/20/2016 - this inelastic-related if block moved after the
!                     return from the elastic-only computation
!  Monochromatic - 234 points (233 shifts + 1 calculation point)
!  Binned        - all points in outer buffer

!      IF ( DO_BIN_REALIZATION ) THEN
!        NPOINTS_TOTAL = NPOINTS_OUTER
!        NPOINTS_LOCAL = NPOINTS_INNER
!      ELSE
!        NPOINTS_TOTAL = NPOINTS_MONO
!        NPOINTS_LOCAL = 1
!      ENDIF

!mick fix 10/19/2015 - added DO_AVOID_UDBEAM lines
!   - Define DO_AVOID_UDBEAM flag
!   - Controls execution of Truncated Postprocessed Direct-Beam Terms
!   -- Rob mod 5/12/17 for 2p5a, Use DO_SSCORR general flag (OUTGOING or NADIR)

      DO_AVOID_UDBEAM = ( DO_SSCORR_GENERAL .or. DO_MSMODE_LRRS )

!  Legendre polynomial setup
!   -- Rob mod 5/12/17 for 2p5a, Use DO_SSCORR_ALONE flag

      IF ( .NOT. DO_SSCORR_ALONE ) THEN
        CALL LEGENDRE_SETUPS &
            ( FOURIER,  NSTREAMS, NMOMENTS,                  & ! Input
              DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS, & ! Input
              COS_SZA, QUAD_STREAMS, QUAD_WEIGHTS,           & ! Input
              PLM_PXI,     PLM_MXI,     & ! Output
              PLM_PXUI,    PLM_MXUI,    & ! Output
              PLM_00_PXI,  PLM_00_MXI,  & ! Output
              PLM_00_PXUI, PLM_00_MXUI, & ! Output
              PLM_WT_PXI,  PLM_WT_MXI )   ! Output
      ENDIF

!  Surface flag (for inclusion of some kind of reflecting boundary)
!    Fake-Lambertian test: set DO_INCLUDE_SURFACE only for M = 0

      DO_INCLUDE_SURFACE    = .TRUE.
      IF ( .not. DO_BRDF_SURFACE ) THEN
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

!mick fix 10/19/2015 - define some local lin control vars
!  note: obtained from L_ELASTIC_MASTER_1
!  Set up the Linearization control

      IF ( DO_PROFILE_WFS ) THEN
        KS(0)    = 0
        KF(0)    = 0
        KPARS(0) = 0
        KVARY(0) = .FALSE.
        DO N = 1, NLAYERS
          KS(N) = 1
          KF(N) = N
          KPARS(N) = LAYER_VARY_NUMBER(N)
          KVARY(N) = LAYER_VARY_FLAG(N)
        ENDDO
      ENDIF

      IF ( DO_COLUMN_WFS ) THEN
        KS(0)    = 0
        KF(0)    = 0
        KPARS(0) = N_TOTALCOLUMN_WFS
        KVARY(0) = .TRUE.
        DO N = 1, NLAYERS
          KS(N)    = 0
          KF(N)    = 0
          KPARS(N) = 0
          KVARY(N) = .FALSE.
        ENDDO
      ENDIF

!  Surface stuff
!  -------------

!  Whole setup with BRDF + Lambertian fraction + Fouriers; Waterleaving, 
!   --> All replaced by supplement-derived inputs Version 2.5

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

!  This is the elastic-scattering solution with 'Rayleigh' scattering
!  Same for both Cabannes/Raman and Energy Balancing methods

!mick fix 10/19/2015 -
!   - Replace DO_SSCORR_GENERAL with DO_AVOID_UDBEAM, in Line # 2
!       of subroutine arguments to L_ELASTIC_MASTER_1. Formerly --
!      CALL L_ELASTIC_MASTER_1 &
!       ( DO_ELASTIC_ONLY, DO_ELASTIC_POSTPROCESSING, DO_MSMODE_LRRS, &
!         DO_INCLUDE_MVOUT, DO_INCLUDE_SURFACE, DO_SSCORR_OVERALL, &
!         ................ etc...

!  Added 10/10/18. OMP inputs and timing output

      CALL L_ELASTIC_MASTER_1 ( DO_TIMING, NTHREADS,                           & ! OMP Input
         DO_ELASTIC_ONLY, DO_ELASTIC_POSTPROCESSING, DO_MSMODE_LRRS,           & ! Inputs
         DO_INCLUDE_MVOUT, DO_INCLUDE_SURFACE, DO_AVOID_UDBEAM,                & ! Inputs
         DO_UPWELLING, DO_DNWELLING, DO_BIN_REALIZATION, DO_PLANE_PARALLEL,    & ! Inputs
         DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,         & ! Inputs
         NKSTORAGE, NKSTORAGE2,                                                & ! Inputs
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS,                & ! Inputs
         N_SURFACE_WFS, N_SLEAVE_WFS,                                          & ! Inputs
         DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_BRDF_Wav1,   & ! Inputs
         ALBEDOS_RANKED, BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,         & ! Inputs
         DO_SLEAVE_Wav1, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,        & ! Inputs
         LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,             & ! Inputs
         LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,         & ! Inputs
         NLAYERS, NSTREAMS, NMOMENTS, N_LOUTPUT, N_USER_STREAMS,               & ! Inputs
         FOURIER, COS_SZA, TAYLOR_SMALL, TAYLOR_ORDER, FLUX_FACTOR,            & ! Inputs
         NPOINTS_OUTER, OFFSET_INNER, NPOINTS_INNER, NPOINTS_MONO, W_EXCIT,    & ! Inputs
         QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS,               & ! Inputs
         STERM_MASK_UP, LEVELMASK_UP, STERM_MASK_DN, LEVELMASK_DN,             & ! Inputs
         DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, SAVE_TRANS_USERM,               & ! Inputs
         BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS, BEAM_ETRANS,  & ! Inputs
         PLM_PXI, PLM_MXI, PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI, PLM_MXUI,         & ! Inputs
         PLM_00_PXI, PLM_00_MXI, PLM_00_PXUI, PLM_00_MXUI,                     & ! Inputs
         L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC, L_SAVE_TRANS_USERM,         & ! Inputs (Linearized)
         L_BEAM_ITRANS, L_BEAM_AVSECANT, L_BEAM_ETRANS, L_BEAM_DTRANS,         & ! Inputs (Linearized)
         L0_KEIGEN, L0_KTRANS, L0_WPARTIC,                                     & ! Outputs
         L0_WHOM_XPOS,  L0_WHOM_LCON,  L0_WHOM_MCON,                           & ! Outputs
         ELASTIC_F_UP, MEAN_ELASTIC_UP, FLUX_ELASTIC_UP,                       & ! Outputs
         ELASTIC_F_DN, MEAN_ELASTIC_DN, FLUX_ELASTIC_DN,                       & ! Outputs
         L_L0_KEIGEN, L_L0_KTRANS, L_L0_WPARTIC, L_L0_WHOM_XPOS,                         & ! Outputs (Linearized)
         L_L0_WHOM_LCON, L_L0_WHOM_MCON, LS_L0_WHOM_LCON,  LS_L0_WHOM_MCON,              & ! Outputs (Linearized)
         L_ELASTIC_F_UP, L_ELASTIC_F_DN, LS_ELASTIC_F_UP, LS_ELASTIC_F_DN,               & ! Outputs (Linearized)
         L_MEAN_ELASTIC_UP,  L_FLUX_ELASTIC_UP,  L_MEAN_ELASTIC_DN,  L_FLUX_ELASTIC_DN,  & ! Outputs (Linearized)
         LS_MEAN_ELASTIC_UP, LS_FLUX_ELASTIC_UP, LS_MEAN_ELASTIC_DN, LS_FLUX_ELASTIC_DN, & ! Outputs (Linearized)
         FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE, OMP_Elastic_Time )                       ! Outputs

!write(*,*)'Baseline',fourier,ELASTIC_F_UP(1:3,1,1)
!write(*,*)'Baseline',fourier,LS_ELASTIC_F_UP(1,1:3,1,1)

!  Return conditions

      IF ( FAIL ) RETURN
      IF ( DO_ELASTIC_ONLY ) RETURN

!  #################################################
!  INELASTIC CALCULATION - MAIN LOOP OVER ALL POINTS
!  #################################################

!  Progress...

!      IF ( NSTREAMS .GT. 10 ) PROG = 5
!      IF ( NSTREAMS .LE. 10 ) PROG = 30
!      IF ( NSTREAMS .LE.  6 ) PROG = 10
!      IF ( NSTREAMS .LE.  6 ) PROG = 50

      IF ( PROGRESS .EQ. 0 ) THEN
        PROG = 60000
      ELSE
        PROG = PROGRESS
      ENDIF

!  Monochromatic - 234 points (233 shifts + 1 calculation point)
!  Binned        - all points in outer buffer

      IF ( DO_BIN_REALIZATION ) THEN
        NPOINTS_TOTAL = NPOINTS_OUTER
        NPOINTS_LOCAL = NPOINTS_INNER
      ELSE
        NPOINTS_TOTAL = NPOINTS_MONO
        NPOINTS_LOCAL = 1
      ENDIF

!  0/7/17. added extra timing for OpenMP

      if ( DO_TIMING ) call cpu_time(omp_e1)

!  @@@@@@@@@@@@@@@@@@@@@
!  Begin parallel region
!  @@@@@@@@@@@@@@@@@@@@@

!  Shared stack for This call
!mick fix 11/28/2018 - added lin vars: DO_ATMOS_WFS, KS, KF, KPARS, KVARY

!  Rob Trial settings 10/10/18. Added Linearized I/O to the original stack.
!     -- All new linearized stuff is marked with the "**!" at the end of the lines.

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (OMP_NTHREADS, DO_TIMING, PROG,                                   & ! Inputs (OMP)
!$OMP     DO_CABANNES_RAMAN, DO_ENERGY_BALANCING,                                & ! Inputs (Module)
!$OMP     DO_MSMODE_LRRS, DO_DIRECTRRS_ONLY, DO_RRS_OVERALL,                     & ! Inputs (Module)
!$OMP     DO_INCLUDE_MVOUT, DO_INCLUDE_SURFACE, DO_SSCORR_GENERAL,               & ! Inputs (Module)
!$OMP     DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                           & ! Inputs (Module)
!$OMP     DO_BIN_REALIZATION, DO_MONO_REALIZATION, DO_BRDF_SURFACE,              & ! Inputs (Module)
!$OMP     DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, do_BRDF_Wav1, do_SLEAVE_Wav1,     & ! Inputs (Module)
!$OMP     ALBEDOS_RANKED, BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,          & ! Inputs (Module)
!$OMP     SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                         & ! Inputs (Module)
!$OMP     NLAYERS, NSTREAMS, NMOMENTS, N_LOUTPUT, N_USER_STREAMS,                & ! Inputs (Module)
!$OMP     FOURIER, COS_SZA, TAYLOR_SMALL, TAYLOR_ORDER, FLUX_FACTOR,             & ! Inputs (Module)
!$OMP     NPOINTS_LOCAL, NSTR2, NTOTAL, N_SUBDIAG, N_SUPDIAG, DELTA_FACTOR,      & ! Inputs (Internal)
!$OMP     FLUX_MULTIPLIER, FLUXES_RANKED, DO_REFLECTED_DIRECTBEAM,               & ! Inputs (Internal)
!$OMP     M, SURFACE_FACTOR_0,                                                   & ! Inputs (Internal)
!$OMP     OFFSET_INNER, NPOINTS_INNER, W_EXCIT, NPOINTS_MONO, N_RRSBINS, BINMAP, & ! Inputs
!$OMP     QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS,                & ! Inputs
!$OMP     STERM_MASK_UP, LEVELMASK_UP, STERM_MASK_DN, LEVELMASK_DN,              & ! Inputs
!$OMP     DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, OMEGAMOMS_CABANNES,              & ! Inputs
!$OMP     OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN, OMEGAMOMS_RRSGAIN,                & ! Inputs
!$OMP     DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,          & ! Inputs (Module Linearized)  **!
!$OMP     LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS,                 & ! Inputs (Module Linearized)  **!
!$OMP     N_SURFACE_WFS, N_SLEAVE_WFS, NKSTORAGE, NKSTORAGE2, DO_PLANE_PARALLEL, & ! Inputs (Module Linearized)  **!
!$OMP     BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS, BEAM_ETRANS,   & ! Inputs
!$OMP     PLM_PXI, PLM_MXI, PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI, PLM_MXUI,          & ! Inputs
!$OMP     PLM_00_PXI, PLM_00_MXI, PLM_00_PXUI, PLM_00_MXUI, SAVE_TRANS_USERM,    & ! Inputs
!$OMP     L0_KEIGEN, L0_KTRANS, L0_WPARTIC,                                      & ! Inputs
!$OMP     L0_WHOM_XPOS,  L0_WHOM_LCON,  L0_WHOM_MCON,                            & ! Inputs
!$OMP     L_BEAM_ITRANS, L_BEAM_AVSECANT, L_BEAM_DTRANS, L_BEAM_ETRANS,          & ! Inputs (Module Linearized)  **!
!$OMP     L_SAVE_TRANS_USERM,                                                    & ! Inputs (Module Linearized)  **!
!$OMP     LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,              & ! Inputs (Module Linearized)  **!
!$OMP     LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,          & ! Inputs (Module Linearized)  **!
!$OMP     L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC, L_OMEGAMOMS_CABANNES,        & ! Inputs (Module Linearized)  **!
!$OMP     L_OMEGAMOMS_RRSLOSS, L_OMEGAMOMS_RRSBIN,  L_OMEGAMOMS_RRSGAIN,         & ! Inputs (Module Linearized)  **!
!$OMP     L_L0_KEIGEN, L_L0_KTRANS, L_L0_WPARTIC, L_L0_WHOM_XPOS,                & ! Inputs (Module Linearized)  **!
!$OMP     L_L0_WHOM_LCON, L_L0_WHOM_MCON, LS_L0_WHOM_LCON,  LS_L0_WHOM_MCON,     & ! Inputs (Module Linearized)  **!
!$OMP     DO_ATMOS_WFS, KS, KF, KPARS, KVARY,                                    & ! Inputs (Internal Linearized)  **!
!$OMP     RAMAN_F_UP, MEAN_RAMAN_UP, FLUX_RAMAN_UP,                              & ! Outputs
!$OMP     RAMAN_F_DN, MEAN_RAMAN_DN, FLUX_RAMAN_DN,                              & ! Outputs
!$OMP     L_RAMAN_F_UP, L_RAMAN_F_DN, LS_RAMAN_F_UP, LS_RAMAN_F_DN,              & ! Outputs (Linearized)  **!
!$OMP     L_MEAN_RAMAN_UP, L_FLUX_RAMAN_UP, L_MEAN_RAMAN_DN, L_FLUX_RAMAN_DN,    & ! Outputs (Linearized)  **!
!$OMP     LS_MEAN_RAMAN_UP, LS_FLUX_RAMAN_UP, LS_MEAN_RAMAN_DN, LS_FLUX_RAMAN_DN,& ! Outputs (Linearized)  **!
!$OMP     LRRS_TID_SAVE, FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE )                 ! Bookkeeping

!  Obtain thread number

   TID = OMP_GET_THREAD_NUM()
   !TID = 0
   write(*,'(1x,a,i1,a)') 'Thread ',TID,' beginning its wavenumber loop'

!  Obtain and display total number of threads and local thread

   IF (TID == 0) THEN
      OMP_NTHREADS = OMP_GET_NUM_THREADS()
      !OMP_NTHREADS = 1
      !write(*,*)
      !write(*,'(1x,a,i1)') 'Total number of threads (inside parallel region) = ', OMP_NTHREADS
      write(*,'(1x,a,i1,a)') 'Running L_Raman_fourier with ',OMP_NTHREADS,' thread(s)'
   END IF

!  Start points loop
!  =================

!$OMP DO

      DO W = 1, NPOINTS_LOCAL

!  For binning, set counters in inner (APT) and outer (CPT) loops
!    Use Local variable W_EXCIT_0 inside parallel region (OMP). Added 10/10/18

        IF ( DO_BIN_REALIZATION ) THEN
          APT = W
          CPT = APT + OFFSET_INNER
          IF ( MOD(APT,PROG).EQ.0 ) THEN
            WRITE(*,'(A,I4)')'    --> Doing Calculation point # ',APT
          ENDIF
          W_EXCIT_0 = CPT
        ELSE
          APT = 1
          W_EXCIT_0 = W_EXCIT
          CPT = W_EXCIT_0
        ENDIF

!  Indices for surface quantities (New, Version 2.5, 9/11/15)

        Raman_IDX   = CPT
        Brdf_Idx    = 1 ; if ( .not. do_BRDF_Wav1   ) Brdf_IDX   = CPT
        Sleave_Idx  = 1 ; if ( .not. do_SLEAVE_Wav1 ) Sleave_IDX = CPT

!  Surface factor

        SURFACE_FACTOR = SURFACE_FACTOR_0

!  Local layer optical depths and layer transmittances

        DO N = 1, NLAYERS
          LOCAL_DELTAUS(N) = DELTAU_VERT_INPUT(N,CPT)
          DO UM = 1, N_USER_STREAMS
            LOCAL_TRANS_USERM(UM,N) = SAVE_TRANS_USERM(UM,N,CPT)
          ENDDO
        ENDDO

!  Linearizations of layer optical depths + transmittances

        IF ( DO_ATMOS_WFS ) THEN
         DO N = 1, NLAYERS
          DO Q = 1, LAYER_VARY_NUMBER(N)
           L_LOCAL_DELTAUS(Q,N) = L_DELTAU_VERT_INPUT(Q,N,CPT)
           DO UM = 1, N_USER_STREAMS
            L_LOCAL_TRANS_USERM(Q,UM,N) = L_SAVE_TRANS_USERM(Q,UM,N,CPT)
           ENDDO
          ENDDO
         ENDDO
        ENDIF

!  Local IOPs (omega-moments).
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

!  Linearized Local IOPs

        IF ( DO_ATMOS_WFS ) THEN
         IF ( DO_CBCOPY ) THEN
          DO N = 1, NLAYERS
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO L = 0, NMOMENTS
             L_LOCAL_OMEGAMOMS(Q,N,L) = L_OMEGAMOMS_CABANNES(Q,N,L,APT)
            ENDDO
           ENDDO
          ENDDO
         ELSE
          DO N = 1, NLAYERS
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO L = 0, NMOMENTS
             L_LOCAL_OMEGAMOMS(Q,N,L) = L_OMEGAMOMS_ELASTIC(Q,N,L,CPT)
            ENDDO
           ENDDO
          ENDDO
         ENDIF
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

!  Number of solutions for Energy balancing solution (includes loss terms)

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

!  Version 2.4: Add SLTERM flag and Isotropic SLTERM. @@@ RobFix 06 sep 12
!  Version 2.5 : replaced BRDF/SLEAVE inputs by supplement-derived variables

        if (lg)write(lu,'(a)')'starting DBEAM_SETUP'
        CALL DBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,               & ! input
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,               & ! input
             FOURIER, Raman_IDX, Brdf_IDX, Sleave_IDX,          & ! Input
             NSTREAMS, DELTA_FACTOR, FLUX_FACTOR,               & ! input
             COS_SZA, BEAM_ETRANS(NLAYERS,CPT), ALBEDOS_RANKED, & ! input
             BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0,            & ! input
             ATMOS_ATTN, DIRECT_BEAM )                            ! Output

!  RT differential equation solutions (HOMOGENEOUS)
!  ================================================

!  Discrete ordinate solutions
!  ---------------------------

!  Layer loop

        DO LAYER = 1, NLAYERS

          if(lg)write(lu,'(a,i2)')'starting HOMOG_SOLUTION ',LAYER

!  Eigensolver for homogeneous solutions

          CALL HOMOG_SOLUTION &
              ( LAYER, FOURIER, NSTREAMS, NMOMENTS,             & ! Input
                LOCAL_OMEGAMOMS, LOCAL_DELTAUS,                 & ! Input
                QUAD_STREAMS, QUAD_WEIGHTS, PLM_PXI, PLM_MXI,   & ! Input
                XPOS, XNEG, KEIGEN, KTRANS, EIGENMAT_SAVE,      & ! Output
                EIGENVEC_SAVE, DIFVEC_SAVE, DAB_SAVE, SAB_SAVE, & ! Output
                FAIL_SUB, MESSAGE_SUB )                           ! Output

!  error handling. Modified for OMP usage 10/10/18

          IF ( FAIL_SUB ) THEN
            WRITE(C2,'(I2)')N   ; MESSAGE     = 'L_Raman_Fourier master, HOMOG_SOLUTION, layer'//C2
            WRITE(C3,'(I3)')CPT ; POINT_TRACE = 'point number: '//C3
            FAIL = .true. ; GOTO 69
          ENDIF

!  Linearization control

          IF ( DO_COLUMN_WFS ) THEN
            NPARS   = KPARS(0)
            DO_VARY = .TRUE.
          ELSE IF ( DO_PROFILE_WFS ) THEN
            NPARS   = KPARS(LAYER)
            DO_VARY = KVARY(LAYER)
          ENDIF

!  Linearized Eigensolutions
!    - Equally valid for profile or column linearization

          IF ( DO_ATMOS_WFS ) THEN

            CALL LPC_HOMOG_SOLUTION &
            ( LAYER, FOURIER, NSTREAMS, NMOMENTS, DO_VARY, NPARS, & ! Input
              L_LOCAL_OMEGAMOMS, LOCAL_DELTAUS, L_LOCAL_DELTAUS,  & ! Input
              QUAD_STREAMS, QUAD_WEIGHTS, PLM_PXI, PLM_MXI,       & ! Input
              KEIGEN, KTRANS, EIGENMAT_SAVE,                      & ! Input
              EIGENVEC_SAVE, DIFVEC_SAVE, DAB_SAVE, SAB_SAVE,     & ! Input
              L_EIGENMAT_SAVE, L_DAB_SAVE, L_SAB_SAVE,            & ! Input
              L_XPOS, L_XNEG, L_KEIGEN, L_KTRANS,                 & ! Output
              FAIL_SUB, MESSAGE_SUB )                               ! Output

!  error handling. Modified for OMP usage 10/10/18

            IF ( FAIL_SUB ) THEN
              WRITE(C2,'(I2)')N   ; MESSAGE     = 'L_Raman_Fourier master, Call to LPC_HOMOG_SOLUTION, layer'//C2
              WRITE(C3,'(I3)')CPT ; POINT_TRACE = 'Linearized Homog. Calculation: point number: '//C3
              FAIL = .true. ; GOTO 69
            ENDIF

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

!  Linearized norms

          IF ( DO_ATMOS_WFS ) THEN
            DO AA = 1, NSTREAMS
              DO Q = 1, NPARS
                NORM = ZERO
                DO I = 1, NSTREAMS
                  I1 = I + NSTREAMS
                  T1 = 2.0d0 * L_XPOS(Q,I,AA,LAYER)  * XPOS(I,AA,LAYER)
                  T2 = 2.0d0 * L_XPOS(Q,I1,AA,LAYER) * XPOS(I1,AA,LAYER)
                  NORM = NORM + QUAD_STRMWGT(I) * ( T1 - T2 )
                ENDDO
                L_EIGENNORM_SAVE(Q,AA,LAYER) = NORM
              ENDDO
            ENDDO
          ENDIF

!  end layer loop

        ENDDO

!  Additional Green function quantities
!  ------------------------------------

!  Norms of the solution vectors

!        CALL HOMOG_SOLUTION_NORMS
!     I      ( NLAYERS, NSTREAMS, QUAD_STRMWGT, XPOS,
!     O        EIGENNORM_SAVE )

!  Linearized norms

!        IF ( DO_ATMOS_WFS ) THEN
!          CALL L_HOMOG_SOLUTION_NORMS
!     I      ( NLAYERS, NSTREAMS, KVARY, KPARS,
!     I        QUAD_STRMWGT, XPOS, L_XPOS,
!     O        L_EIGENNORM_SAVE )
!        ENDIF

!  Homogeneous solution user-stream solutions and multipliers
!  ==========================================================

!  Start layer loop

        DO N = 1, NLAYERS

!  Linearization control

          IF ( DO_COLUMN_WFS ) THEN
            NPARS   = KPARS(0)
            DO_VARY = .TRUE.
          ENDIF
          IF ( DO_PROFILE_WFS ) THEN
            NPARS   = KPARS(N)
            DO_VARY = KVARY(N)
          ENDIF

!  homogeneous solutions at user-defined streams

          CALL UHOMOG_SOLUTION &
                ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS,   & ! Inputs
                  LOCAL_OMEGAMOMS, XPOS, XNEG,                      & ! Inputs
                  PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI,                 & ! Inputs
                  U_XPOS, U_XNEG, U_HELP_P(1,0,N), U_HELP_M(1,0,N) )  ! Output

!  Linearized homogeneous solutions at user-defined streams

          IF ( DO_ATMOS_WFS ) THEN
            CALL LPC_UHOMOG_SOLUTION &
            ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS,     & ! Inputs
              DO_VARY, NPARS, LOCAL_OMEGAMOMS, L_LOCAL_OMEGAMOMS, & ! Inputs
              L_XPOS, L_XNEG, PLM_WT_PXI, PLM_WT_MXI,             & ! Inputs
              PLM_PXUI, U_HELP_P(1,0,N), U_HELP_M(1,0,N),         & ! Inputs
              L_U_XPOS, L_U_XNEG )                                  ! Outputs
          ENDIF

!  Whole layer multipliers
!   @@@ RobFix 5/5/11. Small numbers analysis: added LRRS_FGSMALL,
!                      DELTAU_VERT_INPUT(N,CPT) to Argument list.

          CALL HOMOGMULT_1 &
                ( NSTREAMS, N, N_USER_STREAMS,               & ! Inputs
                  TAYLOR_SMALL, TAYLOR_ORDER,                & ! Inputs
                  LOCAL_TRANS_USERM, USER_STREAMS,           & ! Inputs
                  KEIGEN, KTRANS, DELTAU_VERT_INPUT(N,CPT),  & ! Inputs
                  HMULT_1, HMULT_2, ZETA_P, ZETA_M )           ! Outputs

!  Linearized homogeneous multipliers
!   @@@ RobFix 5/5/11. Small numbers analysis: added LRRS_FGSMALL,
!                      DELTAU_VERT_INPUT(N,CPT) and
!                      L_DELTAU_VERT_INPUT(:,N,CPT) to Argument list.

          IF ( DO_ATMOS_WFS ) THEN
            CALL LPC_HOMOGMULT_1 &
             ( NSTREAMS, N, N_USER_STREAMS, DO_VARY, NPARS,                 & ! Inputs
               TAYLOR_ORDER, TAYLOR_SMALL, USER_STREAMS, LOCAL_TRANS_USERM, & ! Inputs
               L_LOCAL_TRANS_USERM, KTRANS, L_KEIGEN, L_KTRANS,             & ! Inputs
               DELTAU_VERT_INPUT(N,CPT), L_DELTAU_VERT_INPUT(:,N,CPT),      & ! Inputs
               HMULT_1, HMULT_2, ZETA_P, ZETA_M,                            & ! Inputs
               L_HMULT_1, L_HMULT_2 )                                         ! Outputs
          ENDIF

!  End layer loop

        ENDDO

!  Boundary Value Matrix
!  =====================

!  Additional setups for the albedo layer

        if(lg)write(lu,'(a)')'SURFACE_HOMOG_SETUP'
        CALL SURFACE_HOMOG_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR, & ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX,     & ! Inputs
           QUAD_STRMWGT, ALBEDOS_RANKED, BRDF_F, XPOS, XNEG,    & ! Inputs
           R2_HOMP, R2_HOMM )                                     ! Outputs

!  Linearization control for the surface contribution

        IF ( DO_COLUMN_WFS ) THEN
          NPARS   = KPARS(0)
          DO_VARY = .TRUE.
        ENDIF
        IF ( DO_PROFILE_WFS ) THEN
          NPARS   = KPARS(NLAYERS)
          DO_VARY = KVARY(NLAYERS)
        ENDIF

!  Linearized Reflectance of homogeneous solutions

        IF ( DO_ATMOS_WFS ) THEN
          if(lg)write(lu,'(a)')'starting L_SURFACE_HOMOG_SETUP'
          CALL L_SURFACE_HOMOG_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR,   & ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX,       & ! Inputs
           DO_VARY, NPARS, QUAD_STRMWGT, ALBEDOS_RANKED, BRDF_F,  & ! Inputs
           KTRANS, L_KTRANS, XPOS, L_XPOS, L_XNEG,                & ! Inputs
           L_R2_HOMP, L_R2_HOMM )                                   ! Outputs
        ENDIF

!  initialize compression matrix

        if(lg)write(lu,'(a)')'starting BVPMATRIX_INIT'
        CALL BVPMATRIX_INIT &
           ( NSTREAMS, NLAYERS,                   & ! Input
             NTOTAL, NSTR2, N_SUBDIAG, N_SUPDIAG, & ! Input
             BMAT_ROWMASK, BANDMAT2 )               ! Output

!  set up boundary values matrix in compressed form (the "A" as in AX=B)

        if(lg)write(lu,'(a)')'starting BVPMATRIX_SETUP'
        CALL BVPMATRIX_SETUP &
           (  DO_INCLUDE_SURFACE,                     & ! Inputs
              NSTREAMS, NLAYERS, NSTR2, BMAT_ROWMASK, & ! Inputs
              XPOS, XNEG, R2_HOMP, R2_HOMM, KTRANS,   & ! Inputs
              BANDMAT2, SMAT2 )                         ! Output

!  LAPACK LU-decomposition for band matrix

        IF (NLAYERS .GT. 1 ) THEN

          CALL DGBTRF &
            ( NTOTAL, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
              BANDMAT2, MAX_BANDTOTAL, IPIVOT, INFO )

!  (Error tracing). Modified for OMP usage 10/10/18

          IF ( INFO .GT. 0 ) THEN
            WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGBTRF: Singular matrix, u(i,i)=0, for i = '//C2
          ELSE IF ( INFO .LT. 0 ) THEN
            WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGBTRF: argument i illegal value, for i = '//C2
          ENDIF

          IF ( INFO .NE. 0 ) THEN
            WRITE(C3, '(I3)' ) CPT ; POINT_TRACE = 'L_Raman_Fourier Radiances/Jacobians, point number: '//C3
            FAIL = .true. ; GOTO 69
          ENDIF

!  SVD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

!  New for Version 2.5

        ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK LU-decomposition for  matrix

          CALL DGETRF ( NTOTAL, NTOTAL, SMAT2, MAX_2_STREAMS, SIPIVOT, INFO )

!  (Error tracing). Modified for OMP usage 10/10/18

          IF ( INFO .GT. 0 ) THEN
            WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGETRF: argument i illegal value, for i = '//C2
            WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Raman_Fourier Radiances/Jacobians, point number: '//C3
            FAIL = .true. ; GOTO 69
          ENDIF

        ENDIF

!  Find particular solutions for elastic and all RRS terms
!  =======================================================

!  Get source vectors for wavelength point APT
!   -- Rob mod 5/12/17 for 2p5a, Use DO_SSCORR general flag (OUTGOING or NADIR)
!   -- mick fix 11/28/2018 - replace "W_EXCIT" with "W_EXCIT_0" for 2p5a

        if (lg) write(lu, '(a)' )'starting L_SOURCES_MASTER_1'
        CALL L_SOURCES_MASTER_1 &
       ( FOURIER, APT, DO_INCLUDE_SURFACE, DO_ENERGY_BALANCING,               & ! Inputs
         DO_RRS_OVERALL, DO_DIRECTRRS_ONLY, DO_MSMODE_LRRS,                   & ! Inputs
         DO_BIN_REALIZATION, DO_MONO_REALIZATION, DO_SSCORR_GENERAL,          & ! Inputs
         DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                         & ! Inputs
         DO_COLUMN_WFS, DO_PROFILE_WFS, DO_SURFACE_WFS,                       & ! Inputs
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS, N_SURFACE_WFS,& ! Inputs
         N_RRSBINS, BINMAP, NPOINTS_MONO, OFFSET_INNER, W_EXCIT_0,            & ! Inputs
         NLAYERS, NSTREAMS, NMOMENTS, N_USER_STREAMS,                         & ! Inputs
         QUAD_WEIGHTS, USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER,              & ! Inputs
         STERM_MASK_UP, STERM_MASK_DN, FLUXES_RANKED,                         & ! Inputs
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
         FAIL, MESSAGE_SUB )                                     ! Outputs

!  exception handling. Modified for OMP usage 10/10/18

        IF ( FAIL ) THEN
          MESSAGE = 'Error from L_SOURCES_MASTER_1'
          IF ( DO_BIN_REALIZATION ) THEN
            WRITE(C3, '(I3)' ) CPT ; POINT_TRACE = 'L_Raman_Fourier Radiances/Jacobians, point number: '//C3
          ENDIF
          FAIL = .true. ; GOTO 69
        ENDIF

!  Debug code
!      do n = 1, 25
!         write(68,'(i4,1p8e19.10)')n,(WLOWER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(WLOWER(i+8,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(WUPPER(i,n),i=1,8)
!         write(68,'(i4,1p8e19.10)')n,(WUPPER(i+8,n),i=1,8)
!      enddo
!      do n = 1, 25
!         write(69,'(i4,1p8e19.10)')n,(L_WLOWER(1,i,n),i=1,8)
!         write(69,'(i4,1p8e19.10)')n,(L_WLOWER(1,i+8,n),i=1,8)
!         write(69,'(i4,1p8e19.10)')n,(L_WUPPER(1,i,n),i=1,8)
!         write(69,'(i4,1p8e19.10)')n,(L_WUPPER(1,i+8,n),i=1,8)
!      enddo
!      do n = 1, 25
!         write(68,'(i4,1p8e19.10)')
!     *           n,WLAYER_PIST_UP(n,1), WLAYER_PIST_DN(n,1)
!      enddo
!      do n = 1, 25
!         write(69,'(i4,1p8e19.10)')
!     *           n,L_WLAYER_PIST_UP(1,n,1), L_WLAYER_PIST_DN(1,n,1)
!      enddo
!        if ( apt.eq.1) pause'end sources plus'


!  ****************
!  ****************
!  Radiation Field:
!  ****************
!  ****************

!  Complete and Solve boundary value problem
!  =========================================

!  Add contributions to the BV problem RHS

        if (lg)write(lu ,'(a)')' **starting SURFACE_BEAM_SETUP'
        CALL SURFACE_BEAM_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR,& ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX,    & ! Inputs
           QUAD_STRMWGT, ALBEDOS_RANKED, BRDF_F, RAMAN_WLOWER, & ! Inputs
           R2_BEAM )                                             ! Output

!  set up Column COL2 for solution vector (the "B" as in AX=B)
!mick mod 7/20/2016 - added debug var ELASTIC_CALL

        if (lg)write(lu ,'(a)')' **starting BVPCOLUMN_SETUP'
        CALL BVPCOLUMN_SETUP &
           ( ELASTIC_CALL, DO_INCLUDE_SURFACE,                 & ! Inputs
             NSTREAMS, NLAYERS, NTOTAL, NSTR2,                 & ! Inputs
             DIRECT_BEAM, RAMAN_WUPPER, RAMAN_WLOWER, R2_BEAM, & ! Inputs
             COL2, SCOL2 )                                       ! Outputs

!  LAPACK substitution using RHS column vector COL2
!    Single layer case (SCOL2). Version 2.5, 9/11/15
!    Exception handling modifield for OMP 9/7/17.

        IF ( NLAYERS .GT. 1 ) THEN

          CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAX_BANDTOTAL, IPIVOT, &
              COL2, MAX_TOTAL, INFO )

          IF ( INFO .LT. 0 ) THEN
            WRITE(C3, '(I3)' ) INFO
            MESSAGE_SUB = 'argument i illegal value, for i = '//C3
            MESSAGE     = 'L_RAMAN_FOURIER_1, DGBTRS call'
            IF ( DO_BIN_REALIZATION ) THEN
              WRITE(C3,'(I3)')CPT ; POINT_TRACE = 'L_Raman_Fourier Radiances/Jacobians, point number: '//C3
            ENDIF
            FAIL = .true. ; GOTO 69
          ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all l

          DO N = 1, NLAYERS
            C0 = (N-1)*NSTR2
            DO I = 1, NSTREAMS
              I1 = I+NSTREAMS
              LCON(I,N) = COL2(C0+I,1)
              MCON(I,N) = COL2(C0+I1,1)
            ENDDO
          ENDDO

!  special case, 1 layer.
!    More general code, Version 2.5, 9/11/15

        ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2
!    Exception handling modifield for OMP 10/10/18.

          CALL DGETRS ( 'N', NTOTAL, 1, &
                      SMAT2, MAX_2_STREAMS, SIPIVOT, &
                      SCOL2, MAX_2_STREAMS, INFO )

          IF ( INFO .LT. 0 ) THEN
            WRITE(C3, '(I3)' ) INFO ; MESSAGE_SUB = 'argument i illegal value, for i = '//C3
            MESSAGE     = 'L_RAMAN_FOURIER_1, DGETRS call (1 layer)'
            IF ( DO_BIN_REALIZATION ) THEN
              WRITE(C3,'(I3)')CPT ; POINT_TRACE = 'L_Raman_Fourier Radiances/Jacobians, point number: '//C3
            ENDIF
            FAIL = .true. ; GOTO 69
          ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all l

          DO I = 1, NSTREAMS
            I1 = I+NSTREAMS
            LCON(I,1) = SCOL2(I,1)
            MCON(I,1) = SCOL2(I1,1)
          ENDDO
        ENDIF

!  Quadrature solutions: multiplied by integration constants

        DO N = 1, NLAYERS
          DO I = 1, NSTR2
            DO AA = 1, NSTREAMS
              LCON_XVEC(I,AA,N) = LCON(AA,N) * XPOS(I,AA,N)
              MCON_XVEC(I,AA,N) = MCON(AA,N) * XNEG(I,AA,N)
            ENDDO
          ENDDO
        ENDDO

!  Get the homogeneous solution layer source terms
!  ===============================================

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

!  Rob Fix 26 March 2012
!   - Replace DO_SSCORR_GENERAL with DO_AVOID_UDBEAM, in Line # 1
!       of subroutine arguments to UDBEAM_SETUP. Formerly --
!          CALL UDBEAM_SETUP &
!           ( DO_INCLUDE_SURFACE, DO_SSCORR_GENERAL, &
!             ........etc......

!  Earlier code and comments, now gone......9/11/15
!  Add SLTERM flag and Isotropic SLTERM. @@@ RobFix 06 sep 12
!     (Old line)    DO_LAMBERTIAN_SURFACE, FL1(CPT), FL2(CPT), USER_BIREFLEC_0, &
!       CALL UDBEAM_SETUP &
!           ( DO_INCLUDE_SURFACE, DO_SSCORR_GENERAL, &
!             DO_LAMBERTIAN_SURFACE, FL1(CPT), FL2(CPT), USER_BIREFLEC_0, &
!             DO_LRRS_SLTERM, LRRS_SLTERM(CPT), & ! @@@ Rob Fix 06 Sep 12
!             FOURIER, N_USER_STREAMS, ATMOS_ATTN, &
!             USER_DIRECT_BEAM )

!  Add SLTERM flag and Isotropic SLTERM. @@@ RobFix 06 sep 12

          CALL UDBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,       & ! Input
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,       & ! Input
             DO_SSCORR_GENERAL, FOURIER,                & ! Input
             Raman_IDX, Brdf_IDX, Sleave_IDX,           & ! Input
             N_USER_STREAMS, DELTA_FACTOR, FLUX_FACTOR, & ! Input
             ATMOS_ATTN, ALBEDOS_RANKED, USER_BRDF_F_0, & ! Input
             SLTERM_ISOTROPIC, USER_SLTERM_F_0,         & ! Input
             USER_DIRECT_BEAM )                           ! Output

!  Source term
!mick mod 7/20/2016 - added debug var ELASTIC_CALL

          CALL BOASOURCE &
          ( ELASTIC_CALL,                                                    & ! Inputs
            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_REFLECTED_DIRECTBEAM,    & ! Inputs
            NSTREAMS, NLAYERS, N_USER_STREAMS, FOURIER, Raman_IDX, Brdf_IDX, & ! Inputs
            QUAD_STRMWGT, RAMAN_WLOWER, LCON_XVEC, MCON_XVEC, KTRANS,        & ! Inputs
            SURFACE_FACTOR, ALBEDOS_RANKED, USER_BRDF_F, USER_DIRECT_BEAM,   & ! Inputs
            IDOWNSURF, BOA_SOURCE, DIRECT_BOA_SOURCE )                         ! Outputs

        ENDIF

!  Downwelling (no source)..................

        IF ( DO_DNWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            TOA_SOURCE(UM) = ZERO
          ENDDO
        ENDIF

!  RAMAN output and MIFLUX
!  =======================

!  Upwelling Intensity

        IF ( DO_UPWELLING ) THEN
          if (lg)write(lu ,'(a)')' **starting RAMAN_INTENSITY_UP_1'
          CALL RAMAN_INTENSITY_UP_1 &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER,    & ! Inputs
            FOURIER, NLAYERS, NSTREAMS, N_LOUTPUT, N_USER_STREAMS, & ! Inputs
            LEVELMASK_UP, RAMAN_WUPPER, RAMAN_WLOWER, KTRANS,      & ! Inputs
            LCON_XVEC, MCON_XVEC, LOCAL_TRANS_USERM,               & ! Inputs
            BOA_SOURCE, DIRECT_BOA_SOURCE,                         & ! Inputs
            WLAYER_HOST_UP, WLAYER_PIST_UP,                        & ! Inputs
            CUMSOURCE_UP, RAMAN_F_UP(1,1,APT), QUADRAMAN_F_UP(1,1,APT) ) ! Outputs

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

!  Downwelling Intensity

        IF ( DO_DNWELLING ) THEN
          if (lg)write(lu ,'(a)')' **starting RAMAN_INTENSITY_DN_1'
          CALL RAMAN_INTENSITY_DN_1 &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER, & ! Inputs
            FOURIER, NSTREAMS, N_LOUTPUT, N_USER_STREAMS,       & ! Inputs
            LEVELMASK_DN, RAMAN_WLOWER, KTRANS,                 & ! Inputs
            LCON_XVEC, MCON_XVEC, LOCAL_TRANS_USERM,            & ! Inputs
            TOA_SOURCE, WLAYER_HOST_DN, WLAYER_PIST_DN,         & ! Inputs
            CUMSOURCE_DN, RAMAN_F_DN(1,1,APT), QUADRAMAN_F_DN(1,1,APT) ) ! Outputs

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

!  mean value output
!  -----------------

!  Check BEAM_ETRANS(1,CPT), formerly BEAM_ETRANS(1,APT)

        IF ( DO_INCLUDE_MVOUT ) THEN
          if (lg)write(lu ,'(a)')' **starting MIFLUX_INTENSITY_1'
          CALL MIFLUX_INTENSITY_1 &
           ( DO_UPWELLING, DO_DNWELLING,             & ! Inputs
             COS_SZA,  FLUX_FACTOR,                  & ! Inputs
             NSTREAMS, N_LOUTPUT, LEVELMASK_DN,      & ! Inputs
             BEAM_PICUTOFF(CPT), BEAM_ETRANS(1,CPT), & ! Inputs
             QUAD_WEIGHTS, QUAD_STRMWGT,             & ! Inputs
             QUADRAMAN_F_UP(1,1,APT), QUADRAMAN_F_DN(1,1,APT), & ! outputs
             MEAN_RAMAN_UP(1,APT), FLUX_RAMAN_UP(1,APT),       & ! outputs
             MEAN_RAMAN_DN(1,APT), FLUX_RAMAN_DN(1,APT) )        ! outputs
        ENDIF

!  ****************************
!  ****************************
!  Atmospheric Jacobians Field:
!  ****************************
!  ****************************

        IF ( DO_ATMOS_WFS ) THEN

!  Set up the Linearization control

          IF ( DO_PROFILE_WFS ) THEN
            WFS = 1
            WFF = NLAYERS
          ELSE IF ( DO_COLUMN_WFS ) THEN
            WFS = 0
            WFF = 0
          ENDIF

!  Start Major Jacobian loop
!  -------------------------

          DO K = KS(WFS), KF(WFF)

!  Atmospheric Linearization of Direct Beam term

            IF ( DO_INCLUDE_SURFACE .AND. ATMOS_ATTN.GT.ZERO ) THEN
              NK = NKSTORAGE(NLAYERS,K)
              DO I = 1, NSTREAMS
                REFLEC = DIRECT_BEAM(I) / BEAM_ETRANS(NLAYERS,CPT)
                DO Q = 1, KPARS(K)
                  L_ATTN = L_BEAM_ETRANS(Q,NK,CPT)
                  L_DIRECT_BEAM(Q,I) = L_ATTN * REFLEC
                ENDDO
              ENDDO
            ELSE
              DO I = 1, NSTREAMS
                DO Q = 1, KPARS(K)
                  L_DIRECT_BEAM(Q,I) = ZERO
                ENDDO
              ENDDO
            ENDIF

!  Copy linearized BVP particular integrals

            DO N = 1, NLAYERS
              NK = NKSTORAGE2(N,K)
              DO I = 1, NSTR2
                DO Q = 1, KPARS(K)
                  L_WLOWER(Q,I,N) = L_RAMAN_WLOWER(Q,I,NK)
                  L_WUPPER(Q,I,N) = L_RAMAN_WUPPER(Q,I,NK)
                ENDDO
              ENDDO
            ENDDO

!  Linearized  Boundary value problem
!  ----------------------------------

!  Linearized Beam surface reflectance term, for the last layer

            if(lg)write(lu,'(a)')'*starting L_SURFACE_BEAM_SETUP'
            CALL L_SURFACE_BEAM_SETUP &
            ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR,       & ! Inputs
              FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX, KVARY(K), & ! Inputs
              KPARS(K), QUAD_STRMWGT, ALBEDOS_RANKED, BRDF_F, L_WLOWER,  & ! Inputs
              L_R2_BEAM )                                                   ! Output

!  set up Column COL2WF for solution vector (the "B" as in AX=B)

            IF ( K .EQ. 0 ) THEN
              CALL LC_BVPCOLUMN_SETUP &
              ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, KPARS(K), NTOTAL, NSTR2,     & ! Inputs
                LCON, MCON, XPOS, XNEG, L_XPOS, L_XNEG, KTRANS, L_KTRANS,           & ! Inputs
                L_WUPPER, L_WLOWER, L_R2_HOMP, L_R2_HOMM, L_R2_BEAM, L_DIRECT_BEAM, & ! Inputs
                COL2WF, SCOL2WF )                                                     ! Outputs
            ELSE
              CALL LPR_BVPCOLUMN_SETUP &
              ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, K, KPARS(K), & ! Inputs
                NTOTAL, NSTR2, LCON, MCON, XPOS, XNEG, KTRANS,     & ! Inputs
                L_XPOS, L_XNEG,L_KTRANS, L_WUPPER, L_WLOWER,       & ! Inputs
                L_R2_HOMP, L_R2_HOMM, L_R2_BEAM, L_DIRECT_BEAM,    & ! Inputs
                COL2WF, SCOL2WF )                                    ! Outputs
            ENDIF

!  LAPACK substitution using RHS column vector COL2
!    Single layer case (SCOL2). Version 2.5, 9/11/15
!    Exception handling modifield for OMP 10/10/18.

            IF ( NLAYERS .GT. 1 ) THEN

              CALL DGBTRS &
              ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, NPARS, &
                 BANDMAT2, MAX_BANDTOTAL, IPIVOT, &
                 COL2WF, MAX_TOTAL, INFO )

              IF ( INFO .LT. 0 ) THEN
                WRITE(C2, '(I2)' ) INFO ; MESSAGE     = 'DGBTRS: argument i illegal value, for i '//C2
                WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Raman_Fourier Atmospheric Linearization, point number: '//C3
                fail = .true. ; Go To 69
              ENDIF

!  Set linearized integration constants NCON and PCON

              DO Q = 1, KPARS(K)
                DO N = 1, NLAYERS
                  C0 = (N-1)*NSTR2
                  DO I = 1, NSTREAMS
                    I1 = I+NSTREAMS
                    NCON(Q,I,N) = COL2WF(C0+I,Q)
                    PCON(Q,I,N) = COL2WF(C0+I1,Q)
                  ENDDO
                ENDDO
              ENDDO

!    More general code, Version 2.5, 9/11/15

            ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2
!    Exception handling modifield for OMP 10/10/18.

              CALL DGETRS ( 'N', NTOTAL, NPARS, &
                            SMAT2, MAX_2_STREAMS, SIPIVOT, &
                            SCOL2WF, MAX_2_STREAMS, INFO )

              IF ( INFO .LT. 0 ) THEN
                WRITE(C2, '(I2)' ) INFO ; MESSAGE     = 'DGETRS: argument i illegal value, for i '//C2
                WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Raman_Fourier Atmospheric Linearization, point number: '//C3
                fail = .true. ; Go To 69
              ENDIF

!  Set linearized integration constants NCON and PCON

              DO Q = 1, KPARS(K)
                DO I = 1, NSTREAMS
                  I1 = I+NSTREAMS
                  NCON(Q,I,1) = SCOL2WF(I,Q)
                  PCON(Q,I,1) = SCOL2WF(I1,Q)
                ENDDO
              ENDDO

            ENDIF

!  DEbug code, 3 June 2009
!            DO Q = 1, KPARS(K)
!              DO N = 1, NLAYERS
!                DO I = 1, NSTREAMS
!          if (q.eq.1.and.k.eq.16)write(63,5)n,i,ncon(q,i,n),pcon(q,i,n)
!          if (q.eq.1.and.k.eq.16)write(61,5)n,i,lcon(i,n),  mcon(i,n)
!                ENDDO
!              ENDDO
!            ENDDO
! 5    format(2i5,1p2e25.16)
!        if (k.eq.16)pause'after LBVP'

!  Linearizations Homogeneous solution user-solutions & multipliers
!  ================================================================

!  Start layer loop

            DO N = 1, NLAYERS

!  Linearization control

              IF ( DO_COLUMN_WFS ) THEN
                NPARS   = N_TOTALCOLUMN_WFS
                DO_VARY = .TRUE.
              ELSE IF ( DO_PROFILE_WFS ) THEN
                NPARS   = LAYER_VARY_NUMBER(N)
                DO_VARY = LAYER_VARY_FLAG(N)
              ENDIF

!  Linearized homogeneous solutions at user-defined streams

              CALL LPC_UHOMOG_SOLUTION &
              ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS,     & ! Inputs
                DO_VARY, NPARS, LOCAL_OMEGAMOMS, L_LOCAL_OMEGAMOMS, & ! Inputs
                L_XPOS, L_XNEG, PLM_WT_PXI,  PLM_WT_MXI,            & ! Inputs
                PLM_PXUI, U_HELP_P(1,0,N), U_HELP_M(1,0,N),         & ! Inputs
                L_U_XPOS, L_U_XNEG )                                  ! Outputs

!  Linearized Whole layer multipliers
!   @@@ RobFix 5/5/11. Small numbers analysis: added LRRS_FGSMALL,
!                      DELTAU_VERT_INPUT(N,CPT) and
!                      L_DELTAU_VERT_INPUT(:,N,CPT) to Argument list.

              CALL LPC_HOMOGMULT_1 &
               ( NSTREAMS, N, N_USER_STREAMS, DO_VARY, NPARS, TAYLOR_ORDER, & ! Inputs
                 TAYLOR_SMALL, USER_STREAMS, LOCAL_TRANS_USERM,             & ! Inputs
                 L_LOCAL_TRANS_USERM, KTRANS, L_KEIGEN, L_KTRANS,           & ! Inputs
                 DELTAU_VERT_INPUT(N,CPT), L_DELTAU_VERT_INPUT(:,N,CPT),    & ! Inputs
                 HMULT_1, HMULT_2, ZETA_P, ZETA_M,                          & ! Inputs
                 L_HMULT_1, L_HMULT_2 )                                       ! Outputs 

!  End layer loop

            ENDDO

!  Get the Linearized homogeneous solution layer source terms
!  ==========================================================

!  start layer loop

            DO N = 1, NLAYERS

!  Upwelling source terms, for Special case K = N, or K = 0

             IF ( STERM_MASK_UP(N) ) THEN
              IF ( N.EQ.K .OR. K.EQ.0 ) THEN
               DO UM = 1, N_USER_STREAMS
                DO Q = 1, KPARS(K)
                 HOM = ZERO
                 DO AA = 1, NSTREAMS
                  H1   = LCON(AA,N)*U_XPOS(UM,AA,N)
                  L_H1 = NCON(Q,AA,N) *   U_XPOS(UM,AA,N) + &
                         LCON(AA,N)   * L_U_XPOS(Q,UM,AA,N)
                  H2   = MCON(AA,N)*U_XNEG(UM,AA,N)
                  L_H2 = PCON(Q,AA,N) *   U_XNEG(UM,AA,N) + &
                         MCON(AA,N)   * L_U_XNEG(Q,UM,AA,N)
                  HOM = HOM &
                  + H1 * L_HMULT_2(Q,AA,UM,N) + L_H1 * HMULT_2(AA,UM,N) &
                  + H2 * L_HMULT_1(Q,AA,UM,N) + L_H2 * HMULT_1(AA,UM,N)
                 ENDDO
                 L_WLAYER_HOST_UP(Q,N,UM) = HOM
                ENDDO
               ENDDO
              ENDIF
             ENDIF

!  Upwelling source terms, for K /= N, K > 0, Profile linearization

             IF ( STERM_MASK_UP(N) ) THEN
              IF ( N.NE.K .AND. K.NE.0 ) THEN
               DO UM = 1, N_USER_STREAMS
                DO Q = 1, KPARS(K)
                 HOM = ZERO
                 DO AA = 1, NSTREAMS
                  L_H1 = NCON(Q,AA,N) *   U_XPOS(UM,AA,N)
                  L_H2 = PCON(Q,AA,N) *   U_XNEG(UM,AA,N)
                  HOM = HOM + L_H1 * HMULT_2(AA,UM,N) &
                            + L_H2 * HMULT_1(AA,UM,N)
                 ENDDO
                 L_WLAYER_HOST_UP(Q,N,UM) = HOM
                ENDDO
               ENDDO
              ENDIF
             ENDIF

!  downwelling source terms, for Special case K = N, or K = 0

             IF ( STERM_MASK_DN(N) ) THEN
              IF ( N.EQ.K .OR. K.EQ.0 ) THEN
               DO Q = 1, KPARS(K)
                DO UM = 1, N_USER_STREAMS
                 HOM = ZERO
                 DO AA = 1, NSTREAMS
                  H1   = LCON(AA,N)*U_XNEG(UM,AA,N)
                  L_H1 = NCON(Q,AA,N) *   U_XNEG(UM,AA,N) + &
                         LCON(AA,N)   * L_U_XNEG(Q,UM,AA,N)
                  H2   = MCON(AA,N)*U_XPOS(UM,AA,N)
                  L_H2 = PCON(Q,AA,N) *   U_XPOS(UM,AA,N) + &
                         MCON(AA,N)   * L_U_XPOS(Q,UM,AA,N)
                  HOM = HOM &
                  + H1 * L_HMULT_1(Q,AA,UM,N) + L_H1 * HMULT_1(AA,UM,N) &
                  + H2 * L_HMULT_2(Q,AA,UM,N) + L_H2 * HMULT_2(AA,UM,N)
                 ENDDO
                 L_WLAYER_HOST_DN(Q,N,UM) = HOM
                ENDDO
               ENDDO
              ENDIF
             ENDIF

!  downwelling source terms, for K /= N, K > 0, Profile linearization

             IF ( STERM_MASK_DN(N) ) THEN
              IF ( N.NE.K .AND. K.NE.0 ) THEN
               DO Q = 1, KPARS(K)
                DO UM = 1, N_USER_STREAMS
                 HOM = ZERO
                 DO AA = 1, NSTREAMS
                  L_H1 = NCON(Q,AA,N) *   U_XNEG(UM,AA,N)
                  L_H2 = PCON(Q,AA,N) *   U_XPOS(UM,AA,N)
                  HOM = HOM + L_H1 * HMULT_1(AA,UM,N) &
                            + L_H2 * HMULT_2(AA,UM,N)
                 ENDDO
                 L_WLAYER_HOST_DN(Q,N,UM) = HOM
                ENDDO
               ENDDO
              ENDIF
             ENDIF

!  end layer loop

            ENDDO

!  Get the Atmospheric-Linearized TOA/BOA source terms
!  ===================================================

!  Upwelling................................

            IF ( DO_UPWELLING ) THEN

!  Atmospheric Linearization of User direct beam

              DO Q = 1, KPARS(K)
                DO UM = 1, N_USER_STREAMS
                  L_USER_DIRECT_BEAM(Q,UM) = ZERO
                ENDDO
              ENDDO
              IF ( KVARY(K).AND.DO_INCLUDE_SURFACE &
                       .AND.ATMOS_ATTN.NE.ZERO) THEN
                NK = NKSTORAGE(NLAYERS,K)
                DO Q = 1, KPARS(K)
                  L_AT = L_BEAM_ETRANS(Q,NK,CPT) / &
                               BEAM_ETRANS(NLAYERS,CPT)
                  DO UM = 1, N_USER_STREAMS
                    L_USER_DIRECT_BEAM(Q,UM) = L_AT*USER_DIRECT_BEAM(UM)
                  ENDDO
                ENDDO
              ENDIF

!  Linearized BOA term

              CALL LPC_BOASOURCE &
               ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_REFLECTED_DIRECTBEAM,        & ! Inputs
                 NSTREAMS, NLAYERS, N_USER_STREAMS, FOURIER, Raman_IDX, Brdf_IDX,     & ! Inputs
                 K, KPARS(K), KVARY(K), QUAD_STRMWGT, LCON, MCON, XPOS, XNEG, KTRANS, & ! Inputs
                 SURFACE_FACTOR, ALBEDOS_RANKED, USER_BRDF_F, L_USER_DIRECT_BEAM,     & ! Inputs
                 NCON, PCON, L_XPOS, L_XNEG, L_KTRANS, L_WLOWER,                      & ! Inputs
                 L_BOA_SOURCE, L_DIRECT_BOA_SOURCE )                                    ! Outputs

!  End upwelling

            ENDIF

!  Linearized TOA term

            IF ( DO_DNWELLING ) THEN
              CALL LPC_TOASOURCE &
                  ( N_USER_STREAMS, NPARS, L_TOA_SOURCE )
            ENDIF

!  Upwelling Atmospheric Jacobians

            IF ( DO_UPWELLING ) THEN
              if (lg)write(lu ,'(a)')' **starting RAMAN_JACOBIAN_UP_1'
              CALL RAMAN_JACOBIAN_UP_1 &
              ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER,                    & ! Inputs
                FOURIER, NLAYERS, NSTREAMS, N_LOUTPUT, N_USER_STREAMS,                 & ! Inputs
                LEVELMASK_UP, KVARY(K), K, KPARS(K), NKSTORAGE2, LOCAL_TRANS_USERM,    & ! Inputs
                KTRANS, LCON, MCON, XPOS, XNEG, CUMSOURCE_UP, L_LOCAL_TRANS_USERM,     & ! Inputs
                L_KTRANS, L_XPOS, L_XNEG, L_WUPPER, L_WLOWER, NCON, PCON,              & ! Inputs
                L_BOA_SOURCE, L_DIRECT_BOA_SOURCE, L_WLAYER_HOST_UP, L_WLAYER_PIST_UP, & ! Inputs
                L_RAMAN_F_UP(1,0,1,1,APT), L_QUADRAMAN_F_UP(1,0,1,1,APT) )               ! Outputs
            ENDIF

!  Downwelling Atmospheric Jacobians

            IF ( DO_DNWELLING ) THEN
              if (lg)write(lu ,'(a)')' **starting RAMAN_JACOBIAN_DN_1'
              CALL RAMAN_JACOBIAN_DN_1 &
              ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER,                 & ! Inputs
                FOURIER, NSTREAMS, N_LOUTPUT, N_USER_STREAMS,                       & ! Inputs
                LEVELMASK_DN, KVARY(K), K, KPARS(K), NKSTORAGE2, LOCAL_TRANS_USERM, & ! Inputs
                KTRANS, LCON, MCON, XPOS, XNEG, CUMSOURCE_DN, L_LOCAL_TRANS_USERM,  & ! Inputs
                L_KTRANS, L_XPOS, L_XNEG, L_WLOWER, NCON, PCON, L_TOA_SOURCE,       & ! Inputs
                L_WLAYER_HOST_DN, L_WLAYER_PIST_DN,                                 & ! Inputs
                L_RAMAN_F_DN(1,0,1,1,APT), L_QUADRAMAN_F_DN(1,0,1,1,APT) )            ! Outputs
            ENDIF

!  Atmospheric Linearization for Mean/Flux output

            IF ( DO_INCLUDE_MVOUT ) THEN
              if (lg)write(lu ,'(a)')' **starting MIFLUX_JACOBIAN_1'
              CALL MIFLUX_JACOBIAN_1 &
               ( DO_UPWELLING, DO_DNWELLING, NSTREAMS, N_LOUTPUT, LEVELMASK_DN,        & ! Inputs
                 KVARY(K), K, KPARS(K), NKSTORAGE, COS_SZA, FLUX_FACTOR, QUAD_WEIGHTS, & ! Inputs
                 QUAD_STRMWGT, BEAM_PICUTOFF(CPT), L_BEAM_ETRANS(1,1,CPT),             & ! Inputs
                 L_QUADRAMAN_F_UP(1,0,1,1,APT), L_QUADRAMAN_F_DN(1,0,1,1,APT),         & ! Inputs
                 L_MEAN_RAMAN_UP(1,0,1,APT), L_FLUX_RAMAN_UP(1,0,1,APT),               & ! Outputs
                 L_MEAN_RAMAN_DN(1,0,1,APT), L_FLUX_RAMAN_DN(1,0,1,APT) )                ! Outputs
            ENDIF

!  debug, 3 June 2009, 1 pm

!         q = 2
! 5       format(2i5,1p2e25.16)
!         if (k.eq.1 ) then
!          write(*,5)1,APT,RAMAN_F_UP(1,1,APT),RAMAN_F_DN(2,1,APT)
!           write(*,5)1,APT,
!     &     L_RAMAN_F_UP(q,k,1,1,APT),L_RAMAN_F_DN(q,k,2,1,APT)
!         endif

!  End Major Atmospheric Jacobian loop

          ENDDO

!  End Atmospheric Jacobian clause

        ENDIF

!  ***********************
!  ***********************
!  Surface Jacobian Field:
!  ***********************
!  ***********************

        IF ( DO_SURFACE_WFS ) THEN

!  Surface Linearization of Direct beam reflected to quadrature streams
!     Treatment added 21 November 2008. Updated for Version 2.5, 12 October 0215

          LS_DIRECT_BEAM = ZERO
          IF ( DO_INCLUDE_SURFACE ) THEN
            IF ( DO_BRDF_SURFACE ) THEN
              DO Q = 1, N_SURFACE_WFS
                DO I = 1, NSTREAMS
                  LS_DIRECT_BEAM(Q,I) = ATMOS_ATTN * LS_BRDF_F_0(Q,FOURIER,I,BRDF_Idx)
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, N_SURFACE_WFS
                DO I = 1, NSTREAMS
                  LS_DIRECT_BEAM(Q,I) = ATMOS_ATTN
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!  Surface Linearization contributions if flagged
!   (Note - this contribution not present for the Elastic field)

          if(lg)write(lu,'(a)')'*starting LS_SURFACE_BEAM_SETUP'
          CALL LS_SURFACE_BEAM_SETUP &
             ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR, FOURIER,        & ! Inputs
               NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX, N_SURFACE_WFS, QUAD_STRMWGT, & ! Inputs
               ALBEDOS_RANKED, BRDF_F, LS_BRDF_F, RAMAN_WLOWER, LS_RAMAN_WLOWER,    & ! Inputs
               LS_DIFFUSE_BEAM )                                                      ! output

!  Linearized BV Problem for surface Jacobians
!  -------------------------------------------

!  set up Column COL2 for solution vector (the "B" as in AX=B)

          if(lg)write(lu,'(a)')'*starting LSR_BVPCOLUMN_SETUP'
          CALL LSR_BVPCOLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER,       & ! Inputs
             NLAYERS, NSTREAMS, NTOTAL, NSTR2, N_SURFACE_WFS,    & ! Inptus
             SURFACE_FACTOR, Brdf_IDX, QUAD_STRMWGT, LS_BRDF_F,  & ! Inputs
             LCON, MCON, XPOS, XNEG, KTRANS, RAMAN_WLOWER,       & ! Inputs
             R2_HOMP, R2_HOMM, R2_BEAM, LS_RAMAN_WLOWER,         & ! Inputs
             LS_RAMAN_WUPPER, LS_DIRECT_BEAM, LS_DIFFUSE_BEAM,   & ! Inputs
             COL2WF_SURF, SCOL2WF_SURF )                           ! Outputs )

!  General case: LAPACK substitution using RHS column vector COL2WF_SURF
!    Single layer case (SCOL2). Version 2.5, 9/11/15
!    Exception handling modifield for OMP 10/10/18.

          IF ( NLAYERS .GT. 1 ) THEN

            CALL DGBTRS &
            ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SURFACE_WFS, &
               BANDMAT2, MAX_BANDTOTAL, IPIVOT, &
               COL2WF_SURF, MAX_TOTAL, INFO )

            IF ( INFO .LT. 0 ) THEN
              WRITE(C2, '(I2)' ) INFO ; MESSAGE     = 'DGBTRS: argument i illegal value, for i '//C2
              WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Raman_Fourier Surface Linearization, point number: '//C3
              fail = .true. ; Go To 69
            ENDIF

!  Set Linearized integration constants NCON_SURF and PCON_SURF

            DO N = 1, NLAYERS
              C0 = (N-1)*NSTR2
              DO Q = 1, N_SURFACE_WFS
                DO I = 1, NSTREAMS
                  I1 = I+NSTREAMS
                  NCON_SURF(Q,I,N) = COL2WF_SURF(C0+I,Q)
                  PCON_SURF(Q,I,N) = COL2WF_SURF(C0+I1,Q)
                ENDDO
              ENDDO
            ENDDO

!  special case, 1 layer.
!    More general code, Version 2.5, 9/11/15

          ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2WF_SURF
!    Exception handling modifield for OMP 10/10/18.

            CALL DGETRS ( 'N', NTOTAL, N_SURFACE_WFS, &
                       SMAT2, MAX_2_STREAMS, SIPIVOT, &
                       SCOL2WF_SURF, MAX_2_STREAMS, INFO )
 
            IF ( INFO .LT. 0 ) THEN
              WRITE(C2, '(I2)' ) INFO ; MESSAGE     = 'DGETRS: argument i illegal value, for i '//C2
              WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Raman_Fourier Surface Linearization, point number: '//C3
              fail = .true. ; Go To 69
            ENDIF

!  Set Linearized integration constants NCON_SURF and PCON_SURF

            DO Q = 1, N_SURFACE_WFS
              DO I = 1, NSTREAMS
                I1 = I+NSTREAMS
                NCON_SURF(Q,I,1) = SCOL2WF_SURF(I,Q)
                PCON_SURF(Q,I,1) = SCOL2WF_SURF(I1,Q)
              ENDDO
            ENDDO
 
          ENDIF

!  Get the Linearized homogeneous solution layer source terms
!  ==========================================================

          DO N = 1, NLAYERS
            IF ( STERM_MASK_UP(N) ) THEN
              DO Q = 1, N_SURFACE_WFS
                DO UM = 1, N_USER_STREAMS
                  HOM = ZERO
                  DO AA = 1, NSTREAMS
                    LS_H1 = NCON_SURF(Q,AA,N) * U_XPOS(UM,AA,N)
                    LS_H2 = PCON_SURF(Q,AA,N) * U_XNEG(UM,AA,N)
                    HOM = HOM + LS_H1 * HMULT_2(AA,UM,N) &
                              + LS_H2 * HMULT_1(AA,UM,N)
                  ENDDO
                  LS_WLAYER_HOST_UP(Q,N,UM) = HOM
                ENDDO
              ENDDO
            ENDIF
            IF ( STERM_MASK_DN(N) ) THEN
              DO Q = 1, N_SURFACE_WFS
                DO UM = 1, N_USER_STREAMS
                  HOM = ZERO
                  DO AA = 1, NSTREAMS
                    LS_H1 = NCON_SURF(Q,AA,N) * U_XNEG(UM,AA,N)
                    LS_H2 = PCON_SURF(Q,AA,N) * U_XPOS(UM,AA,N)
                    HOM = HOM + LS_H1 * HMULT_1(AA,UM,N) &
                              + LS_H2 * HMULT_2(AA,UM,N)
                  ENDDO
                  LS_WLAYER_HOST_DN(Q,N,UM) = HOM
                ENDDO
              ENDDO
            ENDIF
          ENDDO

!  Get the surface-Linearized BOA source term
!  ==========================================

          IF ( DO_UPWELLING ) THEN

!  Surface Linearization of Direct beam reflected to User streams
!     Treatment added 21 November 2008. Upgraded for Version 2.5

            LS_USER_DIRECT_BEAM = ZERO
!mick fix 10/19/2015 - replaced DO_SSCORR_GENERAL with DO_SSCORR_OUTGOING
!   -- Rob mod 5/12/17 for 2p5a, SSCORR_GENERAL reintroduced

            IF ( .NOT.DO_SSCORR_GENERAL.AND.DO_INCLUDE_SURFACE ) THEN
!            IF ( .NOT.DO_SSCORR_OUTGOING.AND.DO_INCLUDE_SURFACE ) THEN
              IF ( DO_BRDF_SURFACE ) THEN
                DO Q = 1, N_SURFACE_WFS
                  DO UI = 1, N_USER_STREAMS
                    LS_USER_DIRECT_BEAM(Q,UI) = LS_USER_BRDF_F_0(Q,FOURIER,UI,BRDF_Idx) * ATMOS_ATTN
                  ENDDO
                ENDDO
              ELSE
                DO Q = 1, N_SURFACE_WFS
                  DO UI = 1, N_USER_STREAMS
                    LS_USER_DIRECT_BEAM(Q,UI) = ATMOS_ATTN
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

!    This routine is only used by the Raman calculation...

            IF ( DO_USER_STREAMS ) THEN
!mick fix 10/19/2015 - replaced DO_SSCORR_GENERAL with DO_SSCORR_OUTGOING
!   -- Rob mod 5/12/17 for 2p5a, SSCORR_GENERAL reintroduced
              CALL LS_BOASOURCE &
              ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                        & ! Inputs
                DO_SSCORR_GENERAL, NSTREAMS, NLAYERS, N_USER_STREAMS,       & ! Inputs
                FOURIER, Raman_IDX, Brdf_IDX, N_SURFACE_WFS, QUAD_STRMWGT,  & ! Inputs
                IDOWNSURF, SURFACE_FACTOR, ALBEDOS_RANKED,                  & ! Inputs
                USER_BRDF_F, LS_USER_BRDF_F,  LS_USER_DIRECT_BEAM,          & ! Inputs
                LS_RAMAN_WLOWER, NCON_SURF, PCON_SURF, XPOS, XNEG, KTRANS,  & ! Inputs
                LS_BOA_SOURCE, LS_DIRECT_BOA_SOURCE )                         ! Outputs
            ENDIF

!  End upwelling

          ENDIF

!  Linearized downwelling

          IF ( DO_DNWELLING ) THEN
            IF ( DO_USER_STREAMS ) THEN
              CALL LS_TOASOURCE &
                ( N_USER_STREAMS, N_SURFACE_WFS, &
                  LS_TOA_SOURCE )
            ENDIF
          ENDIF

!  Upwelling Surface Jacobians

          IF ( DO_UPWELLING ) THEN
            if (lg)write(lu ,'(a)')' **starting RAMAN_SURFJAC_UP_1'
            CALL RAMAN_SURFJAC_UP_1 &
             ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER,        & ! Inputs
               FOURIER, NLAYERS, NSTREAMS, N_LOUTPUT, N_USER_STREAMS,     & ! Inputs
               N_SURFACE_WFS, LEVELMASK_UP, LOCAL_TRANS_USERM,             & ! Inputs
               XPOS, XNEG, KTRANS, LS_RAMAN_WUPPER, LS_RAMAN_WLOWER,      & ! Inputs
               NCON_SURF, PCON_SURF, LS_BOA_SOURCE, LS_DIRECT_BOA_SOURCE, & ! Inputs
               LS_WLAYER_HOST_UP, LS_WLAYER_PIST_UP,                      & ! Inputs
               LS_RAMAN_F_UP(1,1,1,APT), LS_QUADRAMAN_F_UP(1,1,1,APT) )     ! Outputs
          ENDIF

!  Downwelling Surface Jacobians

          IF ( DO_DNWELLING ) THEN
            if (lg)write(lu ,'(a)')' **starting RAMAN_SURFJAC_DN_1'
            CALL RAMAN_SURFJAC_DN_1 &
             ( DO_USER_STREAMS, DO_INCLUDE_MVOUT, FLUX_MULTIPLIER,        & ! Inputs
               FOURIER, NSTREAMS, N_LOUTPUT, N_USER_STREAMS,              & ! Inputs
               N_SURFACE_WFS, LEVELMASK_DN, LOCAL_TRANS_USERM,            & ! Inputs
               XPOS, XNEG, KTRANS, LS_RAMAN_WLOWER, NCON_SURF, PCON_SURF, & ! Inputs
               LS_TOA_SOURCE, LS_WLAYER_HOST_DN, LS_WLAYER_PIST_DN,       & ! Inputs
               LS_RAMAN_F_DN(1,1,1,APT), LS_QUADRAMAN_F_DN(1,1,1,APT) )     ! Outputs
          ENDIF

!  Surface Linearization for Mean/Flux output

          IF ( DO_INCLUDE_MVOUT ) THEN
            if (lg)write(lu ,'(a)')' **starting MIFLUX_SURFJAC_1'
            CALL MIFLUX_SURFJAC_1 &
              ( DO_UPWELLING, DO_DNWELLING, NSTREAMS, N_LOUTPUT,            & ! Inputs
                N_SURFACE_WFS, QUAD_WEIGHTS, QUAD_STRMWGT,                  & ! Inputs
                LS_QUADRAMAN_F_UP(1,1,1,APT), LS_QUADRAMAN_F_DN(1,1,1,APT), & ! Inputs
                LS_MEAN_RAMAN_UP(1,1,APT), LS_FLUX_RAMAN_UP(1,1,APT),       & ! Outputs
                LS_MEAN_RAMAN_DN(1,1,APT), LS_FLUX_RAMAN_DN(1,1,APT) )        ! Outputs
          ENDIF

!  End surface Jacobians clause

        ENDIF

!  label 69 again defined at bottom of subroutine, as Continuation point Error. OMP usage 10/10/18

69     continue

!  Save OMP thread information

       LRRS_TID_SAVE(w) = TID + 1

!  Finish point loop

      ENDDO

!$OMP END DO

!  End parallel region

!$OMP END PARALLEL

!     @@@@@@@@@@@@@@@@@@@
!     End parallel region
!     @@@@@@@@@@@@@@@@@@@

!      pause 'finished Fourier, lin'

!  Time spent in OpenMP parallel region = OMP_Elastic_Time, to be compared with Ser_Elastic_Time

      if (DO_TIMING) then
        call cpu_time(omp_e2) ; OMP_Raman_Time = (omp_e2 - omp_e1)/REAL(NTHREADS)
      endif

!  Finish
!  ======

      RETURN
      END SUBROUTINE L_RAMAN_FOURIER_1

!  End Module

      END MODULE lrrs_L_fourier_master_1_m
