
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
! #   SUBROUTINE :                                              #
! #             L_ELASTIC_MASTER_1  (master)                    #
! #                                                             #
! #  This module develops complete discrete ordinate elastic    #
! #    scattering solutions, with stored output as required     #
! #    later on for the RRS calculations, and also post-        #
! #    processed elastic-scattering output (if flagged)         #
! #                                                             #
! #  This module is new to Version 2.2 and contains additional  #
! #    software to perform linearization tasks on the elastic   #
! #    solutions. Column/Profile and Surface Jacobians.         #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Use of Supplement-derived BRDF/SLEAVE inputs and control
!    (2) Bookkeeping improvements (use of "Only", clearer I/O specifications)
!    (3) Use of Taylor-series control, TAYLOR_ORDER index, passed to PostProcessing
!    (4) Introduction of variable number of surface Jacobians, dimensioning to match

      MODULE lrrs_L_elastic_master_1_m

!      USE LRRS_PARS_m, Only : LDU

      PRIVATE
      PUBLIC :: L_ELASTIC_MASTER_1

      CONTAINS

!  Added 10/10/18. OMP input and output

      SUBROUTINE L_ELASTIC_MASTER_1 ( DO_TIMING, NTHREADS,                     & ! OMP Input
         DO_ELASTIC_ONLY, DO_ELASTIC_POSTPROCESSING, DO_MSMODE_LRRS,           & ! Inputs
         DO_INCLUDE_MVOUT, DO_INCLUDE_SURFACE, DO_SSCORR_OVERALL,              & ! Inputs
         DO_UPWELLING, DO_DNWELLING, DO_BIN_REALIZATION, DO_PLANE_PARALLEL,    & ! Inputs

         DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,         & ! Inputs
         NKSTORAGE, NKSTORAGE2,                                                & ! Inputs
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS,                & ! Inputs
         N_SURFACE_WFS, N_SLEAVE_WFS,                                          & ! Inputs

         DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, do_BRDF_Wav1,   & ! Inputs
         ALBEDOS_RANKED, BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,         & ! Inputs
         do_SLEAVE_Wav1, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,        & ! Inputs

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
         FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE, OMP_Elastic_Time )                       ! Outputs (Exceptions, OMP Timing)

!  Module of dimensions and numbers

      USE LRRS_PARS_m
!   -- Rob mod 5/12/17 for 2p5a, introduced LRRS_MISCSETUPS_1_m, reorganized LRRS_RTSOLUTIONS_1_m

!  Use Modules (Regular solution)

      USE LRRS_MISCSETUPS_1_m         , Only : DBEAM_SETUP, UDBEAM_SETUP
      USE LRRS_RTSOLUTIONS_1_m        , Only : HOMOG_SOLUTION,  BEAM_SOLUTION_MATRIX, BEAM_SOLUTION_ELASTIC, &
                                               UHOMOG_SOLUTION, UBEAM_SOLUTION_ELASTIC
      USE LRRS_BVPROBLEM_m            , Only : SURFACE_HOMOG_SETUP, SURFACE_BEAM_SETUP,  &
                                               BVPMATRIX_INIT, BVPMATRIX_SETUP, BVPCOLUMN_SETUP
      USE LRRS_AUX2_m                 , Only : DGBTRF, DGBTRS, DGETRS, DGETRF
      USE LRRS_POSTPROCESSING_1_m     , Only : HOMOGMULT_1, BEAMMULT_UP_1, BEAMMULT_DN_1
      USE LRRS_ELASTIC_INTENSITY_1_m  , Only : ELASTIC_INTENSITY_UP_1, ELASTIC_INTENSITY_DN_1
      USE LRRS_RAMAN_INTENSITY_1_m    , Only : MIFLUX_INTENSITY_1

!  Use Modules (Linearized solution)

      USE LRRS_L_RTSOLUTIONS_1_m      , Only : LPC_HOMOG_SOLUTION, LPC_BEAM_SOLUTION_ELASTIC, &
                                               LPC_UHOMOG_SOLUTION,LPC_UBEAM_SOLUTION_ELASTIC
      USE LRRS_L_BVPROBLEM_m          , Only : L_SURFACE_HOMOG_SETUP, L_SURFACE_BEAM_SETUP,  &
                                               LC_BVPCOLUMN_SETUP, LPE_BVPCOLUMN_SETUP, LSE_BVPCOLUMN_SETUP
      USE LRRS_L_POSTPROCESSING_1_m   , Only : LPC_HOMOGMULT_1, LPC_BEAMMULT_UP_1, LPC_BEAMMULT_DN_1
      USE LRRS_L_SLEAVE_m             , Only : LRRS_LSSL_DBSETUPS, LRRS_LSSL_UDBSETUPS, LRRS_LSSL_WFS
      USE LRRS_L_ELASTIC_JACOBIANS_1_m, Only : ELASTIC_JACOBIAN_UP_1, ELASTIC_JACOBIAN_DN_1, &
                                               ELASTIC_SURFJAC_UP_1,  ELASTIC_SURFJAC_DN_1
      USE LRRS_L_RAMAN_JACOBIANS_1_m  , Only : MIFLUX_JACOBIAN_1, MIFLUX_SURFJAC_1

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Top level control
!  -----------------

!  Number of OMP Cores to be used (1, 2, 4, 8...)

      INTEGER  , intent(in) ::  NTHREADS

!  Monitoring flag

      logical  , intent(in) :: DO_TIMING

!  elastic control flags

      LOGICAL  , INTENT(IN) :: DO_ELASTIC_ONLY
      LOGICAL  , INTENT(IN) :: DO_ELASTIC_POSTPROCESSING

!  Multiple scatter mode

      LOGICAL  , INTENT(IN) :: DO_MSMODE_LRRS

!  Output for Flux calculations. Set in the MASTER

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_MVOUT

!  Flag for including surface. set in the MASTER

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  Overall SS correction flag; Now incorporates the DB term

      LOGICAL  , INTENT(IN) :: DO_SSCORR_OVERALL

!  Upwelling and downwelling flags

      LOGICAL  , INTENT(IN) :: DO_UPWELLING
      LOGICAL  , INTENT(IN) :: DO_DNWELLING

!  User stream flag
!      LOGICAL  , INTENT(IN) :: DO_USER_STREAMS

!  Binning realization

      LOGICAL  , INTENT(IN) :: DO_BIN_REALIZATION

!  Plane parallel

      LOGICAL  , INTENT(IN) :: DO_PLANE_PARALLEL

!  Profile/column/surface linearization control inputs

      LOGICAL  , INTENT(IN) :: DO_PROFILE_WFS
      LOGICAL  , INTENT(IN) :: DO_COLUMN_WFS
      LOGICAL  , INTENT(IN) :: DO_SURFACE_WFS
      LOGICAL  , INTENT(IN) :: DO_SLEAVE_WFS

!  Linearization bookkeeping storage

      INTEGER  , INTENT(IN) :: NKSTORAGE   ( MAX_LAYERS, 0:MAX_LAYERS )
      INTEGER  , INTENT(IN) :: NKSTORAGE2  ( MAX_LAYERS, 0:MAX_LAYERS )

!  Linearization control

      LOGICAL  , INTENT(IN) :: LAYER_VARY_FLAG   (MAX_LAYERS)
      INTEGER  , INTENT(IN) :: LAYER_VARY_NUMBER (MAX_LAYERS)
      INTEGER  , INTENT(IN) :: N_TOTALCOLUMN_WFS
      INTEGER  , INTENT(IN) :: N_SURFACE_WFS
      INTEGER  , INTENT(IN) :: N_SLEAVE_WFS

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

      REAL(fpk), intent(in) :: LS_BRDF_F_0      ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LS_BRDF_F        ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

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

      REAL(fpk), intent(in) :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAX_POINTS )
      !REAL(fpk), intent(in) :: LSSL_SLTERM_USERANGLES ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LSSL_SLTERM_F_0        ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LSSL_USER_SLTERM_F_0   ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )

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

!  Small number control. Extended, Version 2.5, 10/12/15

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
      INTEGER  , INTENT(IN) :: W_EXCIT

!  Layer masks for doing integrated source terms
!  ---------------------------------------------

      LOGICAL  , INTENT(IN) :: STERM_MASK_UP ( MAX_LAYERS )
      LOGICAL  , INTENT(IN) :: STERM_MASK_DN ( MAX_LAYERS )

      INTEGER  , INTENT(IN) :: LEVELMASK_UP    (MAX_LOUTPUT)
      INTEGER  , INTENT(IN) :: LEVELMASK_DN    (MAX_LOUTPUT)

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

!  Scaled elastic-scattering optical properties
!  --------------------------------------------

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT   ( MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Vertical optical thickness values

      REAL(FPK), INTENT(IN) :: &
          L_DELTAU_VERT_INPUT   ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_ELASTIC &
           ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Linearized inputs
!  =================

!  Solar beam transmittances, average secant factors

      REAL(FPK), INTENT(IN) :: L_BEAM_ITRANS &
               ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_BEAM_AVSECANT &
               ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_BEAM_ETRANS &
               ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_BEAM_DTRANS &
               ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(IN) :: L_SAVE_TRANS_USERM &
            ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Output arguments
!  ================

!  Local discrete ordinate solutions at shifted wavelengths
!  --------------------------------------------------------

!  Eigenvalues

      REAL(FPK), INTENT(OUT) :: L0_KEIGEN ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Eigentransmittances

      REAL(FPK), INTENT(OUT) :: L0_KTRANS ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Particular solution vectors

      REAL(FPK), INTENT(OUT) :: L0_WPARTIC ( MAX_2_STREAMS, MAX_LAYERS, MAX_POINTS )

!  homogeneous solution vectors and integration constants

      REAL(FPK), INTENT(OUT) :: L0_WHOM_XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L0_WHOM_LCON ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L0_WHOM_MCON ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Linearized Local discrete ordinate solutions at shifted wavelengths
!  -------------------------------------------------------------------

!  Eigenvalues

      REAL(FPK), INTENT(OUT) :: L_L0_KEIGEN ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Eigentransmittances

      REAL(FPK), INTENT(OUT) :: L_L0_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Particular solution vectors

      REAL(FPK), INTENT(OUT) :: L_L0_WPARTIC ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS_NK, MAX_POINTS )

!  homogeneous solution vectors and integration constants

      REAL(FPK), INTENT(OUT) :: L_L0_WHOM_XPOS ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_L0_WHOM_LCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS_SQ, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_L0_WHOM_MCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS_SQ, MAX_POINTS )

!  Surface linearization of integration constants

      REAL(FPK), INTENT(OUT) :: LS_L0_WHOM_LCON ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: LS_L0_WHOM_MCON ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )


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

!  module status

      LOGICAL            , INTENT(OUT) ::  FAIL
      CHARACTER (LEN=120), INTENT(OUT) ::  MESSAGE
      CHARACTER (LEN=120), INTENT(OUT) ::  MESSAGE_SUB
      CHARACTER (LEN=120), INTENT(OUT) ::  POINT_TRACE

!  OMP timing, additional output 10/18/18

      REAL, INTENT(INOUT)              :: OMP_Elastic_Time

!***********************************************************************
!  ARRAYS FOR THE DETERMINATION OF THE ELASTIC DISCRETE ORDINATE SOLUTIO
!     These are based on the LIDORT model with the classical solution
!***********************************************************************

!  Local setups
!  ------------

!  local transmittances along user-stream directions

      REAL(FPK) :: LOCAL_TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS )

!  Homogeneous solutions
!  ---------------------

!  Eigenvector solutions

      REAL(FPK) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Saved matrices for eigenvalue computation

      REAL(FPK) :: DAB ( MAX_STREAMS, MAX_STREAMS)
      REAL(FPK) :: SAB ( MAX_STREAMS, MAX_STREAMS)
      REAL(FPK) :: EIGENMAT ( MAX_STREAMS, MAX_STREAMS)
      REAL(FPK) :: EIGENVEC ( MAX_STREAMS, MAX_STREAMS)
      REAL(FPK) :: DIFVEC   ( MAX_STREAMS, MAX_STREAMS)

!  Reflected solutions at ground

      REAL(FPK) :: R2_HOMP     ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: R2_HOMM     ( MAX_STREAMS, MAX_STREAMS )

!  Help arrays for linearization

      REAL(FPK) :: U_HELP_P  ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK) :: U_HELP_M  ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Particular integrals (for the Chandrasekahar substitution method)
!  --------------------

!  at the layer upper and lower boundaries

      REAL(FPK) :: WUPPER ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )

!  Reflected solutions at ground (direct beam, diffusely reflected beam)

      REAL(FPK) :: DIRECT_BEAM ( MAX_STREAMS )
      REAL(FPK) :: R2_BEAM     ( MAX_STREAMS )

!  Particular solution: LU-decomposed matrix and Pivot

      REAL(FPK) :: QMAT ( MAX_STREAMS, MAX_STREAMS )
      INTEGER   :: QPIVOT    ( MAX_STREAMS )

!  Classical beam solution matrices and vectors
!  Beam solution independent of optical depth (classical solution)
!  Help arrays for linearization

      REAL(FPK) :: QSUMVEC ( MAX_STREAMS )
      REAL(FPK) :: QDIFVEC ( MAX_STREAMS )
      REAL(FPK) :: QVEC    ( MAX_STREAMS )
      REAL(FPK) :: QAUX    ( MAX_STREAMS )

!  Particular solution help variable

      REAL(FPK) :: W_HELP(0:MAX_MOMENTS)

!  Boundary value problem. AX = B
!  ------------------------------

!  Compressed band matrix, as in A

      REAL(FPK) :: BANDMAT2(MAX_BANDTOTAL,MAX_TOTAL)
      REAL(FPK) :: SMAT2   (MAX_2_STREAMS, MAX_2_STREAMS)

!  compression indexing for the band matrix

      INTEGER   :: BMAT_ROWMASK(MAX_TOTAL,MAX_TOTAL)

!  RHS column vector, as in B

      REAL(FPK) :: COL2  ( MAX_TOTAL,1 )
      REAL(FPK) :: SCOL2 ( MAX_2_STREAMS,1 )

!  Pivot for the LU decomposition of the band matrix

      INTEGER   :: IPIVOT ( MAX_TOTAL )
      INTEGER   :: SIPIVOT (MAX_2_STREAMS)

!  Integration constants, as in X

      REAL(FPK) :: LCON ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: MCON ( MAX_STREAMS, MAX_LAYERS )

!  Integration constants multiplied by homogeneous solution vectors


      REAL(FPK) :: LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  post processing solutions (for user-defined polar directions)
!  -------------------------

!  Direct beam transmittance

      REAL(FPK) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS )

!  Homogeneous solutions

      REAL(FPK) :: U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Particular integrals, diffuse term (1), primary scatter term (2)

      REAL(FPK) :: U_WPOS1 ( MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: U_WPOS2 ( MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: U_WNEG1 ( MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: U_WNEG2 ( MAX_USER_STREAMS, MAX_LAYERS )

!  Post processing Multipliers
!  ---------------------------

!  Whole layer multipliers homogeneous solutions


      REAL(FPK) :: HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: ZETA_P  ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: ZETA_M  ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Whole layer solutions, particular integrals

      REAL(FPK) :: EMULT_UP ( MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: EMULT_DN ( MAX_USER_STREAMS, MAX_LAYERS )

!  Cumulative source functions

      REAL(FPK) :: CUMSOURCE_UP ( MAX_USER_STREAMS, 0:MAX_LAYERS)
      REAL(FPK) :: CUMSOURCE_DN ( MAX_USER_STREAMS, 0:MAX_LAYERS)

!  Diffuse BOA source term

      REAL(FPK) :: IDOWNSURF ( MAX_STREAMS )

!  Atmospheric attenuation

      REAL(FPK) :: ATMOS_ATTN

!  Fourier output (quadrature streams)
!    Use as debug, except for Flux calculations

      REAL(FPK) :: QUADELASTIC_F_UP &
              ( MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

      REAL(FPK) :: QUADELASTIC_F_DN &
              ( MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

!  Atmospheric Linearization variables
!  -----------------------------------

!  local linearized  user-stream transmittances

      REAL(FPK) :: L_LOCAL_TRANS_USERM &
           ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS )

!  Linearized Eigenvector solutions

      REAL(FPK) :: L_XPOS ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: L_XNEG ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  LInearized eigenmatrices

      REAL(FPK) :: L_DAB ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: L_SAB ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: L_EIGENMAT ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )

!  Linearized Particular integral at layer boundaries

      REAL(FPK) :: L_WUPPER  ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK) :: L_WLOWER  ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS )

!  Column vectors for weighting functions

      REAL(FPK) :: COL2WF  ( MAX_TOTAL, MAX_ATMOSWFS )
      REAL(FPK) :: SCOL2WF ( MAX_2_STREAMS, MAX_ATMOSWFS )

!  Linearized Integration constants


      REAL(FPK) :: NCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: PCON ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  Linearized Homogeneous solutions at the user defined stream angles

      REAL(FPK) :: L_U_XPOS ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: L_U_XNEG ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Linearized User_angle solution vectors

      REAL(FPK) :: &
            L_U_WPOS1 ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS ), &
            L_U_WPOS2 ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS_NK )

!  Linearized User_angle solution vectors

      REAL(FPK) :: &
            L_U_WNEG1 ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS ), &
            L_U_WNEG2 ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS_NK )

!  Linearized homogeneous solutions multipliers

      REAL(FPK) :: L_HMULT_1 &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: L_HMULT_2 &
            ( MAX_ATMOSWFS, MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Linearized particular integral multipliers

      REAL(FPK) :: L_EMULT_UP &
            ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: L_EMULT_DN &
            ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_USER_STREAMS, MAX_LAYERS )

!  Linearized Reflected homogeneous solutions at ground

      REAL(FPK) :: L_R2_HOMP(MAX_ATMOSWFS,MAX_STREAMS,MAX_STREAMS)
      REAL(FPK) :: L_R2_HOMM(MAX_ATMOSWFS,MAX_STREAMS,MAX_STREAMS)

!  Linearized Reflected beam solutions at ground

      REAL(FPK) :: L_DIRECT_BEAM      ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_R2_BEAM          ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_USER_DIRECT_BEAM ( MAX_ATMOSWFS, MAX_USER_STREAMS )

!  Surface linearization variables
!  -------------------------------

!  Linearized beam solutions at ground

      REAL(FPK) :: LS_DIRECT_BEAM      ( MAX_SURFACEWFS, MAX_STREAMS )
      REAL(FPK) :: LS_USER_DIRECT_BEAM ( MAX_SURFACEWFS, MAX_USER_STREAMS )

      REAL(fpk) :: LSSL_DIRECT_BEAM      ( MAX_SLEAVEWFS, MAX_STREAMS )
      REAL(fpk) :: LSSL_USER_DIRECT_BEAM ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  Column vectors for weighting functions

      REAL(FPK) :: COL2WF_SURF  ( MAX_TOTAL, MAX_SURFACEWFS )
      REAL(FPK) :: SCOL2WF_SURF ( MAX_2_STREAMS, MAX_SURFACEWFS )

!  Linearized Integration constants

      REAL(FPK) :: &
            NCON_SURF ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS ), &
            PCON_SURF ( MAX_SURFACEWFS, MAX_STREAMS, MAX_LAYERS )

!  For the atmospheric weighting functions, discrete ordinate values

      REAL(FPK) :: L_QUADELASTIC_F_UP &
       ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

      REAL(FPK) :: L_QUADELASTIC_F_DN &
       ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

!  For the surface weighting functions, discrete ordinate values

      REAL(FPK) :: LS_QUADELASTIC_F_UP &
              ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

      REAL(FPK) :: LS_QUADELASTIC_F_DN &
              ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

!  Other local help variables
!  --------------------------

!  Original variables

      INTEGER   :: I, I1, AA, N, NC, L, C0, UM, M, K, UI
      INTEGER   :: INFO, CPT, APT, NPOINTS_LOCAL
      INTEGER   :: NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTR2
      INTEGER   :: Raman_IDX, Brdf_IDX, Sleave_IDX
      INTEGER   :: STATUS_SUB
      LOGICAL   :: DO_REFLECTED_DIRECTBEAM
      LOGICAL   :: DO_QUAD_RESULTS, NOFLIP

      REAL(FPK) :: SURFACE_FACTOR, DELFAC, R1, R2
      REAL(FPK) :: FLUX_MULTIPLIER

      LOGICAL   :: DO_ELASTIC_SCATTERING(MAX_LAYERS)

      LOGICAL   :: DO_POST_PROCESSING
      LOGICAL   :: DO_CALCULATION
      LOGICAL   :: FAIL_SUB

      CHARACTER (LEN=2) ::   C2
      CHARACTER (LEN=3) ::   C3

!  Associated with linearization

      LOGICAL   :: DO_ATMOS_WFS

      LOGICAL   :: DO_VARY, KVARY(0:MAX_LAYERS)
      INTEGER   :: Q, NK, KS(0:MAX_LAYERS), KF(0:MAX_LAYERS), LOOP_N
      INTEGER   :: KPARS(0:MAX_LAYERS), NPARS, WFF, WFS
      REAL(FPK) :: REFLEC, L_ATTN, L_AT, AT

!  Local copying for Beam multiplier stuff

      REAL(FPK) :: DELTAUS, DELTRANS, ASOURCE
      INTEGER   :: PICUTOFF

      REAL(FPK) :: L_DELTAUS  ( MAX_ATMOSWFS )
      REAL(FPK) :: L_DELTRANS ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK) :: L_ASOURCE  ( MAX_ATMOSWFS, 0:MAX_LAYERS )

!  Debug

      LOGICAL   :: ELASTIC_CALL=.TRUE.

!  OMP Variables. Added 10/10/18.
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

!  Initialise
!  ----------

!  Added 10/10/18. Set total number of threads (from NTHREADS input)
!     write(*,*) 'Setting total number of threads'

      OMP_MAXTHREADS = NTHREADS
      CALL OMP_SET_NUM_THREADS(OMP_MAXTHREADS)

!  Overall Linearization flag

      DO_ATMOS_WFS = ( DO_COLUMN_WFS .OR. DO_PROFILE_WFS)

!  Set status

      FAIL            = .FALSE.
      MESSAGE         = ' '
      POINT_TRACE     = ' '

!  Local flags

      DO_QUAD_RESULTS = .FALSE.

!  Set number of points
!    - outer number of wavelengths for the binning case
!    - Npoints_mono = 234, monochromatic case

      IF ( DO_BIN_REALIZATION ) THEN
        NPOINTS_LOCAL = NPOINTS_OUTER
      ELSE
        NPOINTS_LOCAL = NPOINTS_MONO
      ENDIF

!  Fourier

      M = FOURIER

!  Start for the actual number of points to be calculated
!Rob fix 10/10/18 - define APT below relative to CPT. OMP usage!

      !APT = 0

!  Set up the elastic scattering flag
!    - always required for Fourier 0, 1 and 2
!    - for Fourier > 2, there may be no scattering (e.g. Rayleigh-only)

      IF ( FOURIER.LT.3 ) THEN
        DO N = 1, NLAYERS
          DO_ELASTIC_SCATTERING(N) = .TRUE.
        ENDDO
      ELSE
        DO N = 1, NLAYERS
          NC = 0
          DO L = FOURIER, 2*NSTREAMS-1
            DO CPT = 1, NPOINTS_LOCAL
              IF ( OMEGAMOMS_ELASTIC(N,L,CPT).NE.ZERO) NC = NC+1
            ENDDO
          ENDDO
          DO_ELASTIC_SCATTERING(N) = (NC .NE. 0)
        ENDDO
      ENDIF

!  Setup INTEGERs for the BV Problem

      NSTR2  = 2 * NSTREAMS
      NTOTAL = NSTR2 * NLAYERS
      IF ( NLAYERS .EQ. 1 ) THEN
        N_SUBDIAG = 2*NSTREAMS - 1
        N_SUPDIAG = 2*NSTREAMS - 1
      ELSE
        N_SUBDIAG = 3*NSTREAMS - 1
        N_SUPDIAG = 3*NSTREAMS - 1
      ENDIF

!  Flux multiplier

       SURFACE_FACTOR = ONE ; DELFAC = TWO
       IF ( FOURIER .EQ. 0 ) THEN
         SURFACE_FACTOR = TWO ; DELFAC = ONE
       ENDIF
       FLUX_MULTIPLIER = FLUX_FACTOR * DELFAC / PI4

!  BRDF Fourier components - coefficient setup
!     Wavelength independent. No Raman scattering.
!   Surface Linearization added 24 November 2008.
!...... COMMENTED OUT for Version 2.5, 9/11/15
!      IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
!        CALL LRRS_BRDF_FOURIER &
!          ( DO_INCLUDE_SURFACE, DO_UPWELLING, DO_USER_STREAMS, &
!            DO_SSCORR_OVERALL, NSTREAMS, N_USER_STREAMS, QUAD_STRMWGT, &
!            N_BRDF, X_BRDF, A_BRDF, FOURIER, DELFAC, &
!            BRDFUNC, BRDFUNC_0, USER_BRDFUNC, USER_BRDFUNC_0, &
!            BIREFLEC, BIREFLEC_0, USER_BIREFLEC, USER_BIREFLEC_0 )
!        IF (DO_SURFACE_WFS.AND..not.DO_LAMBERT_VARIATION)THEN
!           CALL LRRS_BRDF_FOURIER_PLUS &
!           ( DO_INCLUDE_SURFACE, DO_UPWELLING, DO_USER_STREAMS, &
!             DO_SSCORR_OVERALL, NSTREAMS, N_USER_STREAMS, &
!             N_BRDF, X_BRDF, A_BRDF, FOURIER, DELFAC, &
!             LS_BRDFUNC, LS_BRDFUNC_0, LS_USER_BRDFUNC, LS_USER_BRDFUNC_0, &
!             LS_BIREFLEC, LS_BIREFLEC_0, LS_USER_BIREFLEC, LS_USER_BIREFLEC_0 )
!        ENDIF
!      ENDIF

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

!  Following code commented out, JvG 4/28/10
!      IF ( DO_COLUMN_WFS ) THEN
!        DO N = 1, NLAYERS
!          KS(N)    = 0
!          KF(N)    = 0
!          KPARS(0) = N_TOTALCOLUMN_WFS
!          KVARY(0) = .TRUE.
!        ENDDO
!      ENDIF

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

!  Surface reflectance factor R2

      R1 = ONE
      R2 = R1
      IF ( FOURIER .EQ. 0 ) R2 = TWO * R1

!  Direct beam flag (only if above surface flag has been set)

      DO_REFLECTED_DIRECTBEAM = .FALSE.
      IF ( DO_INCLUDE_SURFACE ) THEN
         DO_REFLECTED_DIRECTBEAM = .TRUE.
      ENDIF

!  added extra timing for OpenMP. 9/7/17.

      if ( DO_TIMING ) call cpu_time(omp_e1)

!  @@@@@@@@@@@@@@@@@@@@@
!  Begin parallel region
!  @@@@@@@@@@@@@@@@@@@@@

!  Shared stack for This call
!mick fix 11/28/2018 - added std vars: SURFACE_FACTOR, R1
!                    - added lin vars: DO_ATMOS_WFS, KS, KF, KPARS, KVARY

!  Rob Trial settings 10/10/18. Added Linearized I/O to the original stack.
!     -- All new linearized stuff is marked with the "**!" at the end of the lines.

!$OMP PARALLEL IF(OMP_MAXTHREADS > 1) &
!$OMP   DEFAULT (PRIVATE) &
!$OMP   SHARED (OMP_NTHREADS, DO_TIMING,                                       & ! Inputs (OMP)
!$OMP     DO_ELASTIC_ONLY, DO_ELASTIC_POSTPROCESSING, DO_MSMODE_LRRS,          & ! Inputs (Module)
!$OMP     DO_INCLUDE_MVOUT, DO_INCLUDE_SURFACE, DO_SSCORR_OVERALL,             & ! Inputs (Module)
!$OMP     DO_UPWELLING, DO_DNWELLING, DO_BIN_REALIZATION, DO_BRDF_SURFACE,     & ! Inputs (Module)
!$OMP     DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, do_BRDF_Wav1, do_SLEAVE_Wav1,   & ! Inputs (Module)
!$OMP     ALBEDOS_RANKED, BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,        & ! Inputs (Module)
!$OMP     SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                       & ! Inputs (Module)
!$OMP     NLAYERS, NSTREAMS, NMOMENTS, N_LOUTPUT, N_USER_STREAMS,              & ! Inputs (Module)
!$OMP     NPOINTS_LOCAL, NSTR2, NTOTAL, N_SUBDIAG, N_SUPDIAG, DELFAC,          & ! Inputs (Internal)
!$OMP     FLUX_MULTIPLIER,  DO_ELASTIC_SCATTERING, DO_REFLECTED_DIRECTBEAM,    & ! Inputs (Internal)
!$OMP     DO_QUAD_RESULTS, M, SURFACE_FACTOR, R1, R2,                          & ! Inputs (Internal)
!$OMP     FOURIER, COS_SZA, TAYLOR_SMALL, TAYLOR_ORDER, FLUX_FACTOR,           & ! Inputs (Module)
!$OMP     OFFSET_INNER, NPOINTS_INNER, W_EXCIT,                                & ! Inputs (Module)
!$OMP     QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS,              & ! Inputs (Module)
!$OMP     STERM_MASK_UP, LEVELMASK_UP, STERM_MASK_DN, LEVELMASK_DN,            & ! Inputs (Module)
!$OMP     DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, SAVE_TRANS_USERM,              & ! Inputs (Module)
!$OMP     BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS, BEAM_ETRANS, & ! Inputs (Module)
!$OMP     PLM_PXI, PLM_MXI, PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI, PLM_MXUI,        & ! Inputs (Module)
!$OMP     PLM_00_PXI, PLM_00_MXI, PLM_00_PXUI, PLM_00_MXUI,                    & ! Inputs (Module)
!$OMP     DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,        & ! Inputs (Module Linearized)  **!
!$OMP     NKSTORAGE, NKSTORAGE2, DO_PLANE_PARALLEL,                            & ! Inputs (Module Linearized)  **!
!$OMP     LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS,               & ! Inputs (Module Linearized)  **!
!$OMP     N_SURFACE_WFS, N_SLEAVE_WFS,                                         & ! Inputs (Module Linearized)  **!
!$OMP     LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,            & ! Inputs (Module Linearized)  **!
!$OMP     LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,        & ! Inputs (Module Linearized)  **!
!$OMP     L_BEAM_ITRANS, L_BEAM_AVSECANT, L_BEAM_ETRANS, L_BEAM_DTRANS,        & ! Inputs (Module Linearized)  **!
!$OMP     L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC, L_SAVE_TRANS_USERM,        & ! Inputs (Module Linearized)  **!
!$OMP     DO_ATMOS_WFS, KS, KF, KPARS, KVARY,                                  & ! Inputs (Internal Linearized)  **!
!$OMP     L0_KEIGEN, L0_KTRANS, L0_WPARTIC,                                    & ! Outputs
!$OMP     L0_WHOM_XPOS,  L0_WHOM_LCON,  L0_WHOM_MCON,                          & ! Outputs
!$OMP     ELASTIC_F_UP, MEAN_ELASTIC_UP, FLUX_ELASTIC_UP,                      & ! Outputs
!$OMP     ELASTIC_F_DN, MEAN_ELASTIC_DN, FLUX_ELASTIC_DN,                      & ! Outputs
!$OMP     L_L0_KEIGEN, L_L0_KTRANS, L_L0_WPARTIC, L_L0_WHOM_XPOS,                         & ! Outputs (Linearized)  **!
!$OMP     L_L0_WHOM_LCON, L_L0_WHOM_MCON, LS_L0_WHOM_LCON, LS_L0_WHOM_MCON,               & ! Outputs (Linearized)  **!
!$OMP     L_ELASTIC_F_UP, L_ELASTIC_F_DN, LS_ELASTIC_F_UP, LS_ELASTIC_F_DN,               & ! Outputs (Linearized)  **!
!$OMP     L_MEAN_ELASTIC_UP,  L_FLUX_ELASTIC_UP,  L_MEAN_ELASTIC_DN,  L_FLUX_ELASTIC_DN,  & ! Outputs (Linearized)  **!
!$OMP     LS_MEAN_ELASTIC_UP, LS_FLUX_ELASTIC_UP, LS_MEAN_ELASTIC_DN, LS_FLUX_ELASTIC_DN, & ! Outputs (Linearized)  **!
!$OMP     LRRS_TID_SAVE, FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE )                          ! Bookkeeping

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
      write(*,'(1x,a,i1,a)') 'Running L_Elastic_Fourier with ',OMP_NTHREADS,' thread(s)'
   END IF

!  Start points loop
!  =================

!$OMP DO

      DO CPT = 1, NPOINTS_LOCAL

!  Indices for surface quantities (New, Version 2.5, 9/11/15)

       Raman_IDX   = CPT
       Brdf_Idx    = 1 ; if ( .not. do_BRDF_Wav1   ) Brdf_IDX   = CPT
       Sleave_Idx  = 1 ; if ( .not. do_SLEAVE_Wav1 ) Sleave_IDX = CPT

!  Copy local trans userm

       DO UM = 1, N_USER_STREAMS
        DO N = 1, NLAYERS
          LOCAL_TRANS_USERM(UM,N) = SAVE_TRANS_USERM(UM,N,CPT)
        ENDDO
       ENDDO

!  Linearizations of layer optical depths + transmittances

       IF ( DO_ATMOS_WFS ) THEN
        DO N = 1, NLAYERS
         DO Q = 1, LAYER_VARY_NUMBER(N)
          DO UM = 1, N_USER_STREAMS
           L_LOCAL_TRANS_USERM(Q,UM,N) = L_SAVE_TRANS_USERM(Q,UM,N,CPT)
          ENDDO
         ENDDO
        ENDDO
       ENDIF

!  Local Lambertian surface variation
!       LS_SURFEDO = ZERO
!       IF ( DO_LAMBERT_VARIATION ) LS_SURFEDO = ONE

!  Surface reflectance factor R2
!mick fix 11/28/2018 - removed this 2nd instance of defining R2

!       R2 = R1
!       IF ( FOURIER .EQ. 0 ) R2 = TWO * R1

!  Direct beam reflected to quadrature directions, assumes Flux factor 1
!     BRDF treatment added 2 November 2007.

       CALL DBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,               & ! input
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,               & ! input
             FOURIER, Raman_IDX, Brdf_IDX, Sleave_IDX,          & ! Input
             NSTREAMS, DELFAC, FLUX_FACTOR,                     & ! input
             COS_SZA, BEAM_ETRANS(NLAYERS,CPT), ALBEDOS_RANKED, & ! input
             BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0,            & ! input
             ATMOS_ATTN, DIRECT_BEAM )                            ! Output

!  Earlier code and comments, now gone......9/11/15
!     BRDF treatment added 2 November 2007.
!  Add SLTERM flag and Isotropic SLTERM. @@@ RobFix 06 sep 12
!       CALL DBEAM_SETUP &
!           ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, &
!             DO_LRRS_SLTERM, LRRS_SLTERM(CPT),  & ! @@@ Rob Fix 06 Sep 12 
!             NSTREAMS, FL1(CPT), FL2(CPT), BIREFLEC_0, &
!             COS_SZA, FOURIER, DELFAC, BEAM_ETRANS(NLAYERS,CPT), &
!             ATMOS_ATTN, DIRECT_BEAM )

!  keep track...
!       if ( mod(w,100).eq.0) write(*,*)'-elastic: point #',w

!mick add 10/19/2015 - LRRS_LSSL_DBSETUPS

       IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) THEN
         CALL LRRS_LSSL_DBSETUPS ( &
           DO_SL_ISOTROPIC, DO_REFLECTED_DIRECTBEAM, FOURIER, & ! Input
           Sleave_IDX, NSTREAMS, N_SLEAVE_WFS,                & ! Input
           FLUX_FACTOR, DELFAC,                               & ! Input
           LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,            & ! Input
           LSSL_DIRECT_BEAM )                                   ! Output
       ENDIF

!  ##################################
!  RT differential equation solutions
!  ##################################

!  homogeneous solutions
!  ---------------------

!  Start layer loop

       DO N = 1, NLAYERS

!  Linearization control

        IF ( DO_COLUMN_WFS ) THEN
          NPARS   = KPARS(0)
          DO_VARY = .TRUE.
        ELSE IF ( DO_PROFILE_WFS ) THEN
          NPARS   = KPARS(N)
          DO_VARY = KVARY(N)
        ENDIF

!  Find the solution, by solving the eignproblem.
!    - output eigenvalues L0_KEIGEN, eigentransmittances L0_KTRANS
!       these are required for the RRS calculation, and are stored.
!    -  output saved eigenmatrices and homogeneous solution vectors

        CALL HOMOG_SOLUTION &
            ( N, FOURIER, NSTREAMS, NMOMENTS,                       & ! Input
              OMEGAMOMS_ELASTIC(1,0,CPT), DELTAU_VERT_INPUT(1,CPT), & ! Input
              QUAD_STREAMS, QUAD_WEIGHTS, PLM_PXI, PLM_MXI,         & ! Input
              XPOS, XNEG, L0_KEIGEN(1,1,CPT), L0_KTRANS(1,1,CPT),   & ! Output
              EIGENMAT, EIGENVEC, DIFVEC, DAB, SAB,                 & ! Output
              FAIL_SUB, MESSAGE_SUB )                                 ! Output

!  Exception handling. Modified for OMP usage, 10/10/18.

        IF ( FAIL_SUB ) THEN
          WRITE(C2, '(I2)' ) N   ; MESSAGE = 'L_Elastic_Master, call to HOMOG_SOLUTION, layer'//C2
          WRITE(C3, '(I3)' ) CPT ; POINT_TRACE = 'point number: '//C3
          FAIL = .true. ; GOTO 69
        ENDIF

!  Linearized Eigensolutions
!    - Equally valid for profile or column linearization
!    - Exception handling modified for OMP usage, 10/10/18.

        IF ( DO_ATMOS_WFS ) THEN

          CALL LPC_HOMOG_SOLUTION &
            ( N, FOURIER, NSTREAMS, NMOMENTS, DO_VARY, NPARS,         & ! Input
              L_OMEGAMOMS_ELASTIC(1,1,0,CPT),                         & ! Input
              DELTAU_VERT_INPUT(1,CPT), L_DELTAU_VERT_INPUT(1,1,CPT), & ! Input
              QUAD_STREAMS, QUAD_WEIGHTS, PLM_PXI, PLM_MXI,           & ! Input
              L0_KEIGEN(1,1,CPT),  L0_KTRANS(1,1,CPT),                & ! Input
              EIGENMAT, EIGENVEC, DIFVEC, DAB, SAB,                   & ! Input
              L_EIGENMAT, L_DAB, L_SAB, L_XPOS, L_XNEG,               & ! Output
              L_L0_KEIGEN(1,1,1,CPT), L_L0_KTRANS(1,1,1,CPT),         & ! Output
              FAIL_SUB, MESSAGE_SUB )                                   ! Output

          IF ( FAIL_SUB ) THEN
            WRITE(C2, '(I2)' ) N   ; MESSAGE = 'L_Elastic_Master, call to LPC_HOMOG_SOLUTION, layer'//C2
            WRITE(C3, '(I3)' ) CPT ; POINT_TRACE = 'point number: '//C3
            FAIL = .true. ; GOTO 69
          ENDIF

        ENDIF

!  Store elastic solutions for later use.
!   Save space by storing XPOS, LCON, MCON. April 10th, 2007

        DO AA = 1, NSTREAMS
          DO I = 1, NSTR2
            L0_WHOM_XPOS(I,AA,N,CPT) = XPOS(I,AA,N)
          ENDDO
        ENDDO
        IF ( DO_ATMOS_WFS ) THEN
          DO Q = 1, NPARS
            DO AA = 1, NSTREAMS
              DO I = 1, NSTR2
                L_L0_WHOM_XPOS(Q,I,AA,N,CPT) = L_XPOS(Q,I,AA,N)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Elastic Beam Solution
!  =====================

!  Beam solution Matrix and LU-decomposition
!  -----------------------------------------

!  Linear Algebra setup, AX = B. Creates A in LU-decomposed form.
!     Also create the Linearized Matrix................
!    - Exception handling modified for OMP usage, 10/10/18.

        CALL BEAM_SOLUTION_MATRIX &
              ( N, NSTREAMS, EIGENMAT, &
                BEAM_AVSECANT(N,CPT), BEAM_PICUTOFF(CPT), &
                QMAT, QPIVOT, FAIL_SUB, MESSAGE_SUB )

        IF ( FAIL_SUB ) THEN
          WRITE(C2, '(I2)' ) N   ; MESSAGE = 'L_Elastic_Master, BEAM_SOLUTION_MATRIX, layer'//C2
          WRITE(C3, '(I3)' ) CPT ; POINT_TRACE = 'point number: '//C3
          FAIL = .true. ; GOTO 69
        ENDIF

!  Beam solution solver (Back-substitution)
!  ----------------------------------------

!  Solver for particular integral of the RTE (Elastic)
!   - output = particular integrals at layer boundaries (WUPPER/WLOWER)
!   - output = L0_WPARTIC, the solution source vector for RRS scattering

        CALL BEAM_SOLUTION_ELASTIC &
            ( N, FOURIER,  NMOMENTS, NSTREAMS,                & ! Input
              BEAM_AVSECANT(N,CPT), BEAM_ITRANS(N,CPT),       & ! Input
              QMAT, QPIVOT, SAB, DAB, PLM_00_PXI, PLM_00_MXI, & ! Input
              QUAD_STREAMS, OMEGAMOMS_ELASTIC(1,0,CPT),       & ! Input
              QSUMVEC, QDIFVEC, QVEC, QAUX,                   & ! Output
              L0_WPARTIC(1,1,CPT), FAIL_SUB, MESSAGE_SUB )      ! Output

!  Exception handling for this solution
!    - Exception handling modified for OMP usage, 10/10/18.

        IF ( FAIL_SUB ) THEN
          WRITE(C2, '(I2)' ) N   ; MESSAGE = 'L_Elastic_Master, BEAM_SOLUTION_ELASTIC, layer'//C2
          WRITE(C3, '(I3)' ) CPT ; POINT_TRACE = 'point number: '//C3
          FAIL = .true. ; GOTO 69
        ENDIF

!  Compute Particular integrals at the layer boundaries

        DO I = 1, NSTR2
          WUPPER(I,N) = L0_WPARTIC(I,N,CPT) * BEAM_ITRANS(N,CPT)
          WLOWER(I,N) = L0_WPARTIC(I,N,CPT) * BEAM_ETRANS(N,CPT)
        ENDDO

!  Linearization of the particular integral
!    - output = L_L0_WPARTIC, the solution source vector for RRS scattering.
!    - Exception handling modified for OMP usage, 10/10/18.

        IF ( DO_ATMOS_WFS ) THEN

          CALL LPC_BEAM_SOLUTION_ELASTIC &
            ( N, FOURIER, DO_PLANE_PARALLEL, NMOMENTS, NSTREAMS, & ! Inputs
              KS(N), KF(N), KPARS, KVARY, NKSTORAGE,             & ! Inputs
              QUAD_STREAMS, L_OMEGAMOMS_ELASTIC(1,1,0,CPT),      & ! Inputs
              PLM_00_PXI,  PLM_00_MXI, BEAM_ITRANS(N,CPT),       & ! Inputs
              BEAM_AVSECANT(N,CPT), L_BEAM_AVSECANT(1,1,CPT),    & ! Inputs
              QMAT, QPIVOT, SAB, DAB, QSUMVEC, QDIFVEC,          & ! Inputs
              QVEC, QAUX, L_SAB, L_DAB, L_EIGENMAT,              & ! Inputs
              L_L0_WPARTIC(1,1,1,CPT), FAIL_SUB, MESSAGE_SUB )     ! Outputs

          IF ( FAIL_SUB ) THEN
            WRITE(C2, '(I2)' ) N   ; MESSAGE = 'L_Elastic_Master, LPC_BEAM_SOLUTION_ELASTIC, layer'//C2
            WRITE(C3, '(I3)' ) CPT ; POINT_TRACE = 'point number: '//C3
            FAIL = .true. ; GOTO 69
          ENDIF

!  End linearization clause

        ENDIF

!  End layer loop for all solutions

       ENDDO

!  ######################
!  Boundary Value Problem
!  ######################

!  Matrix Setup section
!  --------------------

!  Additional setups for the lowest layer with surface.
!  Version 2.5, 9/11/15, New setup for the BRDF, using supplement-derived material.

       CALL SURFACE_HOMOG_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, R2,            & ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX,    & ! Inputs
           QUAD_STRMWGT, ALBEDOS_RANKED, BRDF_F, XPOS, XNEG,   & ! Inputs
           R2_HOMP, R2_HOMM )                                    ! Outputs

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

         CALL L_SURFACE_HOMOG_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, R2,              & ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX,      & ! Inputs
           DO_VARY, NPARS, QUAD_STRMWGT, ALBEDOS_RANKED, BRDF_F, & ! Inputs
           L0_KTRANS(1,1,CPT), L_L0_KTRANS(1,1,1,CPT),           & ! Inputs
           XPOS, L_XPOS, L_XNEG,                                 & ! Inputs
           L_R2_HOMP, L_R2_HOMM )                                  ! Outputs

       ENDIF

!  Initialize compression matrix

       CALL BVPMATRIX_INIT &
           ( NSTREAMS, NLAYERS,                   & ! Input
             NTOTAL, NSTR2, N_SUBDIAG, N_SUPDIAG, & ! Input
             BMAT_ROWMASK, BANDMAT2 )               ! Output

!  Set up boundary values matrix in compressed form (the "A" as in AX=B)

       CALL BVPMATRIX_SETUP &
           (  DO_INCLUDE_SURFACE,                               & ! Inputs
              NSTREAMS, NLAYERS, NSTR2, BMAT_ROWMASK,           & ! Inputs
              XPOS, XNEG, R2_HOMP, R2_HOMM, L0_KTRANS(1,1,CPT), & ! Inputs
              BANDMAT2, SMAT2 )                                   ! Output

!  LAPACK LU-decomposition for band matrix

       IF (NLAYERS .GT. 1 ) THEN

        CALL DGBTRF &
           ( NTOTAL, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
             BANDMAT2, MAX_BANDTOTAL, IPIVOT, INFO )

!  Exception handling modified for OMP usage, 10/10/18.

        IF ( INFO .GT. 0 ) THEN
          WRITE(C2, '(I2)' ) INFO
          MESSAGE = 'DGBTRF: Singular matrix, u(i,i)=0, for i = '//C2
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(C2, '(I2)' ) INFO
          MESSAGE = 'DGBTRF: argument i illegal value, for i = '//C2
        ENDIF
        IF ( INFO .NE. 0 ) THEN
          WRITE(C3, '(I3)' ) CPT ; POINT_TRACE = 'L_Elastic_Master Radiances/Jacobians, point number: '//C3
          FAIL = .true. ; GOTO 69
        ENDIF

!  SVD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

!  New for Version 2.5

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK LU-decomposition for  matrix

        CALL DGETRF ( NTOTAL, NTOTAL, SMAT2, MAX_2_STREAMS, SIPIVOT, INFO )

!  Exception handling modified for OMP usage, 10/10/18.

        IF ( INFO .GT. 0 ) THEN
          WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGETRF: argument i illegal value, for i = '//C2
          WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Elastic_Master Radiances/Jacobians, point number: '//C3
          FAIL = .true. ; GOTO 69
        ENDIF

      ENDIF

!  Solution of boundary value problem
!  ----------------------------------

!  Add contributions to the BV problem RHS, for the lowest layer
!  Version 2.5, 9/11/15, New setup for the BRDF, using supplement-derived material.

       CALL SURFACE_BEAM_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, R2,         & ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX, & ! Inputs
           QUAD_STRMWGT, ALBEDOS_RANKED, BRDF_F, WLOWER,    & ! Inputs
           R2_BEAM )                                          ! Output

!  set up Column COL2 for solution vector (the "B" as in AX=B)
!  Additional output for Single layer case (SCOL2). Version 2.5, 9/11/15
!mick mod 7/20/2016 - added debug var ELASTIC_CALL

       CALL BVPCOLUMN_SETUP &
           ( ELASTIC_CALL, DO_INCLUDE_SURFACE,     & ! Inputs
             NSTREAMS, NLAYERS, NTOTAL, NSTR2,     & ! Inputs
             DIRECT_BEAM, WUPPER, WLOWER, R2_BEAM, & ! Inputs
             COL2, SCOL2 )                           ! Outputs

!  LAPACK substitution using RHS column vector COL2
!    Single layer case (SCOL2). Version 2.5, 9/11/15
!    Exception handling modified for OMP usage, 10/10/18.

       IF ( NLAYERS .GT. 1 ) THEN

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAX_BANDTOTAL, IPIVOT, &
              COL2, MAX_TOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGBTRS: argument i illegal value, for i = '//C2
          WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Elastic_Master Radiances/Jacobians, point number: '//C3
          FAIL = .true. ; GOTO 69
        ENDIF

!  special case, 1 layer.
!    More general code, Version 2.5, 9/11/15
!    Exception handling modified for OMP usage, 10/10/18.

       ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS ( 'N', NTOTAL, 1, &
                      SMAT2, MAX_2_STREAMS, SIPIVOT, &
                      SCOL2, MAX_2_STREAMS, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGETRS: argument i illegal value, for i = '//C2
          WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Elastic_Master Radiances/Jacobians, point number: '//C3
          FAIL = .true. ; GOTO 69
        ENDIF

       ENDIF

!  Integration constants to be stored for the RRS calculation
!  ----------------------------------------------------------

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

       DO N = 1, NLAYERS
        C0 = (N-1)*NSTR2
        DO I = 1, NSTREAMS
          I1 = I+NSTREAMS
          LCON(I,N) = COL2(C0+I,1)
          MCON(I,N) = COL2(C0+I1,1)
        ENDDO
       ENDDO

!  Quadrature solutions multiplied by integration constants
!   ( Only required for post processing )

       DO N = 1, NLAYERS
        DO I = 1, NSTR2
          DO AA = 1, NSTREAMS
            LCON_XVEC(I,AA,N) = LCON(AA,N) * XPOS(I,AA,N)
            MCON_XVEC(I,AA,N) = MCON(AA,N) * XNEG(I,AA,N)
          ENDDO
        ENDDO
       ENDDO

!  Copy the above quantities into storage in ELASTICSAVE.VARS
!   Save space by storing XPOS, LCON, MCON. April 10th, 2007

       DO N = 1, NLAYERS
         DO AA = 1, NSTREAMS
           L0_WHOM_LCON(AA,N,CPT) = LCON(AA,N)
           L0_WHOM_MCON(AA,N,CPT) = MCON(AA,N)
         ENDDO
       ENDDO

!  former code............................
!       DO N = 1, NLAYERS
!        DO AA = 1, NSTREAMS
!          AA1 = AA + NSTREAMS
!          DO I = 1, NSTR2
!            L0_WHOMOG_DN(I,AA,N,CPT) = LCON(AA,N) * XPOS(I,AA,N)
!            L0_WHOMOG_UP(I,AA,N,CPT) = MCON(AA,N) * XNEG(I,AA,N)
!          ENDDO
!        ENDDO
!       ENDDO

!  ########################################
!  Calculation of Post Processing Solutions
!  ########################################

!  Only done if the required flags are set

       IF ( DO_ELASTIC_POSTPROCESSING ) THEN

!  Initialize Local post processing flag

         DO_POST_PROCESSING = .FALSE.

!  Only going to be producing post-processed output for:
!    a) 1 point for the monochromatic case, when CPT = W_EXCIT
!    b)  offset < CPT < offset + NPOINTS_INNER + 1, for the binned case

!Rob Fix 10/10/18  - define APT here relative to CPT. OMP Usage   !!!!!
!                  - simplified if structure

         !IF ( DO_ELASTIC_ONLY ) THEN
         !  DO_CALCULATION = .FALSE.
         !  IF ( DO_BIN_REALIZATION ) THEN
         !    IF ( CPT.GT.OFFSET_INNER .AND. &
         !      CPT.LT.OFFSET_INNER + 1 + NPOINTS_INNER ) THEN
         !      DO_CALCULATION = .TRUE.
         !      APT = CPT - OFFSET_INNER
         !    ENDIF
         !  ELSE
         !    IF ( CPT .EQ. W_EXCIT ) THEN
         !      DO_CALCULATION = .TRUE.
         !      APT = 1
         !    ENDIF
         !  ENDIF
         !  IF ( CPT .EQ. W_EXCIT ) THEN
         !    DO_CALCULATION = .TRUE.
         !    APT = 1
         !  ENDIF
         !ELSE
         !  DO_CALCULATION = .FALSE.
         !  IF ( DO_BIN_REALIZATION ) THEN
         !    IF ( CPT.GT.OFFSET_INNER .AND. &
         !      CPT.LT.OFFSET_INNER + 1 + NPOINTS_INNER ) THEN
         !      DO_CALCULATION = .TRUE.
         !      APT = CPT - OFFSET_INNER
         !    ENDIF
         !  ELSE
         !    IF ( CPT .EQ. W_EXCIT ) THEN
         !      DO_CALCULATION = .TRUE.
         !      APT = 1
         !    ENDIF
         !  ENDIF
         !ENDIF

         DO_CALCULATION = .FALSE.
         IF ( DO_BIN_REALIZATION ) THEN
           IF ( CPT .GT. OFFSET_INNER .AND. &
                CPT .LT. OFFSET_INNER + NPOINTS_INNER + 1 ) THEN
             DO_CALCULATION = .TRUE.
             APT = CPT - OFFSET_INNER
           ENDIF
         ELSE
           IF ( CPT .EQ. W_EXCIT ) THEN
             DO_CALCULATION = .TRUE.
             APT = 1
           ENDIF
         ENDIF

         DO_POST_PROCESSING = DO_CALCULATION

       ENDIF

!  Start Post processing solutions
!  -------------------------------

!  Calculation point for the post-processing to be done

       IF ( DO_POST_PROCESSING ) THEN

!  Attenuation of solar beam reflected along user-defined directions
!    Only required for calculations with no SS correction
!  This call updated for Version 2.5, 9/11/15

         CALL UDBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,       & ! Input
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,       & ! Input
             DO_SSCORR_OVERALL, FOURIER,                & ! Input
             Raman_IDX, Brdf_IDX, Sleave_IDX,           & ! Input
             N_USER_STREAMS, DELFAC, FLUX_FACTOR,       & ! Input
             ATMOS_ATTN, ALBEDOS_RANKED, USER_BRDF_F_0, & ! Input
             SLTERM_ISOTROPIC, USER_SLTERM_F_0,         & ! Input
             USER_DIRECT_BEAM )                           ! Output

!  Earlier code and comments, now gone......9/11/15
!  Add SLTERM flag and Isotropic SLTERM. @@@ RobFix 06 sep 12
!     (Old line)    DO_LAMBERTIAN_SURFACE, FL1(CPT), FL2(CPT), USER_BIREFLEC_0, &
!       CALL UDBEAM_SETUP &
!           ( DO_INCLUDE_SURFACE, DO_SSCORR_OVERALL, &
!             DO_LAMBERTIAN_SURFACE, FL1(CPT), FL2(CPT), USER_BIREFLEC_0, &
!             DO_LRRS_SLTERM, LRRS_SLTERM(CPT), & ! @@@ Rob Fix 06 Sep 12
!             FOURIER, N_USER_STREAMS, ATMOS_ATTN, &
!             USER_DIRECT_BEAM )

!mick add 10/19/2015 - LRRS_LSSL_UDBSETUPS
         IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) THEN
           CALL LRRS_LSSL_UDBSETUPS ( &
             DO_SL_ISOTROPIC, DO_REFLECTED_DIRECTBEAM, FOURIER,  & ! Input
             Sleave_IDX, N_USER_STREAMS, N_SLEAVE_WFS,           & ! Input
             FLUX_FACTOR, DELFAC,                                & ! Input
             LSSL_SLTERM_ISOTROPIC, LSSL_USER_SLTERM_F_0,        & ! Input
             LSSL_USER_DIRECT_BEAM )                               ! Output
         ENDIF

!  Post-processed Homogeneous/Particular solutions + All linearizations
!  ====================================================================

!  Individual layer loop

         DO N = 1, NLAYERS

!  Linearization control

           IF ( DO_COLUMN_WFS ) THEN
             NPARS   = N_TOTALCOLUMN_WFS
             DO_VARY = .TRUE.
           ENDIF
           IF ( DO_PROFILE_WFS ) THEN
             NPARS   = LAYER_VARY_NUMBER(N)
             DO_VARY = LAYER_VARY_FLAG(N)
           ENDIF

!  homogeneous solutions at user-defined streams

           CALL UHOMOG_SOLUTION &
                ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS, & ! Inputs
                  OMEGAMOMS_ELASTIC(1,0,CPT), XPOS, XNEG,         & ! Inputs
                  PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI,               & ! Inputs
                  U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )              ! Output

!  Linearized homogeneous solutions at user-defined streams

           IF ( DO_ATMOS_WFS ) THEN
             CALL LPC_UHOMOG_SOLUTION &
              ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS,       & ! Inputs
                DO_VARY, NPARS, OMEGAMOMS_ELASTIC(1,0,CPT),           & ! Inputs
                L_OMEGAMOMS_ELASTIC(1,1,0,CPT), L_XPOS, L_XNEG,       & ! Inputs
                PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI, U_HELP_P, U_HELP_M, & ! Inputs
                L_U_XPOS, L_U_XNEG )                                    ! Outputs
           ENDIF

!  Particular integral solutions at user-defined streams

           CALL UBEAM_SOLUTION_ELASTIC &
            ( DO_UPWELLING, DO_DNWELLING, FOURIER,                    & ! Inputs
              N, STERM_MASK_UP, STERM_MASK_DN,                        & ! Inputs
              NSTREAMS, NMOMENTS, N_USER_STREAMS, BEAM_PICUTOFF(CPT), & ! Inputs
              L0_WPARTIC(1,1,CPT), OMEGAMOMS_ELASTIC(1,0,CPT),        & ! Inputs
              PLM_WT_PXI, PLM_WT_MXI,  PLM_00_PXUI,                   & ! Inputs
              PLM_PXUI,   PLM_00_MXUI, PLM_MXUI,                      & ! Inputs
              W_HELP, U_WPOS1, U_WPOS2, U_WNEG1, U_WNEG2 )              ! Outputs

!  Linearized Particular integral solutions at user-defined streams

           IF ( DO_ATMOS_WFS ) THEN

             CALL LPC_UBEAM_SOLUTION_ELASTIC &
           ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,       & ! Inputs
             FOURIER, N, NSTREAMS, NMOMENTS, N_USER_STREAMS,      & ! Inputs
             KS(N), KF(N), KPARS, KVARY, NKSTORAGE,               & ! Inputs
             STERM_MASK_UP, STERM_MASK_DN, BEAM_PICUTOFF(CPT),    & ! Inputs
             L_L0_WPARTIC(1,1,1,CPT), OMEGAMOMS_ELASTIC(1,0,CPT), & ! Inputs
             L_OMEGAMOMS_ELASTIC(1,1,0,CPT), W_HELP,              & ! Inputs
             PLM_WT_PXI, PLM_WT_MXI, PLM_00_PXUI,                 & ! Inputs
             PLM_PXUI, PLM_00_MXUI, PLM_MXUI,                     & ! Inputs
             L_U_WPOS1, L_U_WPOS2, L_U_WNEG1, L_U_WNEG2 )           ! Outputs

           endif

!  debug. 01 April 2009.

!          write(*,'(i3,1p18e15.5)') n,L_U_WPOS2(1,1,n)

!  End layer loop

         ENDDO

!         pause'end  first pp'

!  Post processing Multipliers + All linearizations
!  ================================================

!  Homogeneous solution multipliers
!  --------------------------------

!  Layer loop

         DO N = 1, NLAYERS

!  Linearization control

           IF ( DO_COLUMN_WFS ) THEN
             NPARS   = N_TOTALCOLUMN_WFS
             DO_VARY = .TRUE.
           ENDIF
           IF ( DO_PROFILE_WFS ) THEN
             NPARS   = LAYER_VARY_NUMBER(N)
             DO_VARY = LAYER_VARY_FLAG(N)
           ENDIF

!  Local copying necessary

           ASOURCE  = BEAM_AVSECANT(N,CPT)
           DELTRANS = BEAM_DTRANS(N,CPT)
           DELTAUS  = DELTAU_VERT_INPUT(N,CPT)
           PICUTOFF = BEAM_PICUTOFF(CPT)

           IF ( DO_ATMOS_WFS ) THEN
!mick fix 10/19/2015 - created LOOP_N & modified loop indexing
             IF (DO_COLUMN_WFS) THEN
               LOOP_N = 0
             ELSE IF (DO_PROFILE_WFS) THEN
               LOOP_N = N
             END IF
             !DO K = KS(N), KF(N)
             DO K = KS(LOOP_N), KF(LOOP_N)
               NK = NKSTORAGE(N,K)
               DO Q = 1, KPARS(K)
                 L_DELTRANS(Q,K) = L_BEAM_DTRANS(Q,NK,CPT)
                 L_ASOURCE(Q,K)  = L_BEAM_AVSECANT(Q,NK,CPT)
               ENDDO
             ENDDO

             DO Q = 1, KPARS(LOOP_N)
               L_DELTAUS(Q) = L_DELTAU_VERT_INPUT(Q,N,CPT)
             ENDDO
           ENDIF

!  homogeneous multipliers for the whole layer
!   @@@ RobFix 5/5/11. Small numbers analysis:
!  Version 2.5, 9/11/15. Expanded Taylor-series control

           CALL HOMOGMULT_1 &
                ( NSTREAMS, N, N_USER_STREAMS,               & ! Inputs
                  TAYLOR_SMALL, TAYLOR_ORDER,                & ! Inputs
                  SAVE_TRANS_USERM(1,1,CPT), USER_STREAMS,   & ! Inputs
                  L0_KEIGEN(1,1,CPT), L0_KTRANS(1,1,CPT),    & ! Inputs
                  DELTAU_VERT_INPUT(N,CPT),                  & ! Inputs
                  HMULT_1, HMULT_2, ZETA_P, ZETA_M )           ! Outputs

!  Linearized homogeneous multipliers
!   @@@ RobFix 5/5/11. Small numbers analysis:
!                      added LRRS_FGSMALL, DELTAUS, L_DELTAUS to Argument list.

           IF ( DO_ATMOS_WFS ) THEN
             CALL LPC_HOMOGMULT_1 &
              ( NSTREAMS, N, N_USER_STREAMS, DO_VARY, NPARS, & ! Inputs
                TAYLOR_ORDER, TAYLOR_SMALL, USER_STREAMS,    & ! Inputs
                LOCAL_TRANS_USERM, L_LOCAL_TRANS_USERM,      & ! Inputs
                L0_KTRANS(1,1,CPT),  L_L0_KEIGEN(1,1,1,CPT), & ! Inputs
                L_L0_KTRANS(1,1,1,CPT), DELTAUS, L_DELTAUS,  & ! Inputs
                HMULT_1, HMULT_2, ZETA_P, ZETA_M,            & ! Inputs
                L_HMULT_1, L_HMULT_2 )                         ! Outputs
           ENDIF

!  multipliers for beam sources
!  ----------------------------

!  Multipliers for the Upwelling solar source functions
!  Version 2.5, 9/11/15. Expanded Taylor-series control

           IF ( DO_UPWELLING ) THEN
             NOFLIP = .FALSE.

             CALL BEAMMULT_UP_1 &
               ( N, N_USER_STREAMS, USER_STREAMS,          & ! Inputs
                 TAYLOR_SMALL, TAYLOR_ORDER,               & ! Inputs
                 STERM_MASK_UP(N), LOCAL_TRANS_USERM(1,N), & ! Inputs
                 DELTAU_VERT_INPUT(N,CPT), NOFLIP,         & ! Inputs
                 BEAM_PICUTOFF(CPT), BEAM_DTRANS(N,CPT),   & ! Inputs
                 BEAM_AVSECANT(N,CPT),                     & ! Inputs
                 EMULT_UP(1,N) )                             ! Output

             IF ( DO_ATMOS_WFS ) THEN
               CALL LPC_BEAMMULT_UP_1 &
               ( N, N_USER_STREAMS, USER_STREAMS, TAYLOR_SMALL,       & ! Inputs
                 TAYLOR_ORDER, KS(N), KF(N), KPARS, STERM_MASK_UP(N), & ! Inputs
                 LOCAL_TRANS_USERM(1,N), L_LOCAL_TRANS_USERM(1,1,N),  & ! Inputs
                 DELTAUS, NOFLIP, PICUTOFF, DELTRANS, ASOURCE,        & ! Inputs
                 EMULT_UP(1,N), L_DELTAUS, L_DELTRANS, L_ASOURCE,     & ! Inputs
                 L_EMULT_UP(1,0,1,N) )                                  ! Outputs
             ENDIF

           ENDIF

!  Multipliers for the Downwelling solar source functions
!  Version 2.5, 9/11/15. Expanded Taylor-series control

           IF ( DO_DNWELLING ) THEN
             NOFLIP = .FALSE.
             CALL BEAMMULT_DN_1 &
               ( N, N_USER_STREAMS, USER_STREAMS,          & ! Inputs
                 TAYLOR_SMALL, TAYLOR_ORDER,               & ! Inputs
                 STERM_MASK_DN(N), LOCAL_TRANS_USERM(1,N), & ! Inputs
                 DELTAU_VERT_INPUT(N,CPT), NOFLIP,         & ! Inputs
                 BEAM_PICUTOFF(CPT), BEAM_DTRANS(N,CPT),   & ! Inputs
                 BEAM_AVSECANT(N,CPT),                     & ! Inputs
                 EMULT_DN(1,N) )                             ! Output

             IF ( DO_ATMOS_WFS ) THEN
               CALL LPC_BEAMMULT_DN_1 &
               ( N, N_USER_STREAMS, USER_STREAMS, TAYLOR_SMALL,       & ! Inputs
                 TAYLOR_ORDER, KS(N), KF(N), KPARS, STERM_MASK_DN(N), & ! Inputs
                 LOCAL_TRANS_USERM(1,N), L_LOCAL_TRANS_USERM(1,1,N),  & ! Inputs
                 DELTAUS, NOFLIP, PICUTOFF, DELTRANS, ASOURCE,        & ! Inputs
                 EMULT_DN(1,N), L_DELTAUS, L_DELTRANS, L_ASOURCE,     & ! Inputs
                 L_EMULT_DN(1,0,1,N))                                   ! Outputs
             ENDIF

           ENDIF

!  End layer loop

         ENDDO

!         pause

!  Apply Initial Transmittances, Upwelling solutions
!     01 April 2009, This is necessary here.
!       * EMULT_UP nor U_WPOS1/U_WPOS2 not multiplied by ITRANS

         IF ( DO_UPWELLING ) THEN
           DO N = 1, NLAYERS
             AT = BEAM_ITRANS(N,CPT)
             IF ( DO_ATMOS_WFS ) THEN
               DO K = KS(N), KF(N)
                 NK = NKSTORAGE(N,K)
                 DO Q = 1, KPARS(K)
                   L_AT = L_BEAM_ITRANS(Q,NK,CPT)
                   DO UM = 1, N_USER_STREAMS
                     L_EMULT_UP(Q,K,UM,N) = AT   * L_EMULT_UP(Q,K,UM,N) &
                                          + L_AT *   EMULT_UP(UM,N)
                   ENDDO
                 ENDDO
               ENDDO
             ENDIF
             DO UM = 1, N_USER_STREAMS
               EMULT_UP(UM,N) = AT * EMULT_UP(UM,N)
             ENDDO
           ENDDO
         ENDIF

!  Apply Initial Transmittances, Downwelling solutions
!     01 April 2009, This is necessary here.
!       * EMULT_DN nor U_WNEG1/U_WNEG2 not multiplied by ITRANS

         IF ( DO_DNWELLING ) THEN
           DO N = 1, NLAYERS
             AT = BEAM_ITRANS(N,CPT)
             IF ( DO_ATMOS_WFS ) THEN
               DO K = KS(N), KF(N)
                 NK = NKSTORAGE(N,K)
                 DO Q = 1, KPARS(K)
                   L_AT = L_BEAM_ITRANS(Q,NK,CPT)
                   DO UM = 1, N_USER_STREAMS
                     L_EMULT_DN(Q,K,UM,N) = AT   * L_EMULT_DN(Q,K,UM,N) &
                                          + L_AT *   EMULT_DN(UM,N)
                   ENDDO
                 ENDDO
               ENDDO
             ENDIF
             DO UM = 1, N_USER_STREAMS
               EMULT_DN(UM,N) = AT * EMULT_DN(UM,N)
             ENDDO
           ENDDO
         ENDIF

!  End section for determining post-processing solutions

       ENDIF

!  ######################################
!  Post processing of the Radiance Fields
!  ######################################

!  Start post-processing clause

       IF ( DO_POST_PROCESSING ) THEN

!  Calculation counter is increased by 1
!mick fix 9/7/2017, 10/10/18- define APT above relative to CPT
 
         !APT = APT + 1

!  Upwelling intensity
!     Added Supplement-derived BRDF inputs, Version 2.5, 9/11/15

         IF ( DO_UPWELLING ) THEN
           CALL ELASTIC_INTENSITY_UP_1 &
          ( DO_MSMODE_LRRS, DO_SSCORR_OVERALL, DO_ELASTIC_SCATTERING,                 & ! Inputs
            DO_QUAD_RESULTS, DO_INCLUDE_MVOUT, DO_INCLUDE_SURFACE,                    & ! Inputs
            DO_REFLECTED_DIRECTBEAM, DO_BRDF_SURFACE, FOURIER, Raman_IDX, Brdf_IDX,   & ! Inputs
            NSTREAMS, NLAYERS, N_USER_STREAMS, N_LOUTPUT, LEVELMASK_UP,               & ! Inputs
            FLUX_MULTIPLIER, QUAD_STRMWGT, SAVE_TRANS_USERM(1,1,CPT),                 & ! Inputs
            R2, ALBEDOS_RANKED, USER_BRDF_F, USER_DIRECT_BEAM,                        & ! Inputs
            LCON, MCON, LCON_XVEC, MCON_XVEC, WUPPER, WLOWER, L0_KTRANS(1,1,CPT),     & ! Inputs
            U_XPOS,  U_XNEG, U_WPOS1, U_WPOS2, HMULT_1, HMULT_2,  EMULT_UP,           & ! Inputs
            IDOWNSURF, CUMSOURCE_UP, ELASTIC_F_UP(1,1,APT), QUADELASTIC_F_UP(1,1,APT) ) ! Outputs
         ENDIF

!  Saved results for the Rayleigh planetary problem

!         IF ( DO_UPWELLING ) THEN
!           IF ( DO_RAYLEIGH_ONLY ) THEN
!             DO J = 1, N_LOUTPUT
!               DO K = 1, N_USER_STREAMS
!                IF(M.EQ.0)ELASTIC_F0_UP(J,K,APT)=ELASTIC_F_UP(J,K,APT)
!                IF(M.EQ.1)ELASTIC_F1_UP(J,K,APT)=ELASTIC_F_UP(J,K,APT)
!                IF(M.EQ.2)ELASTIC_F2_UP(J,K,APT)=ELASTIC_F_UP(J,K,APT)
!               ENDDO
!             ENDDO
!           ENDIF
!         ENDIF

!  Downwelling intensity

         IF ( DO_DNWELLING ) THEN
           CALL ELASTIC_INTENSITY_DN_1 &
          ( DO_MSMODE_LRRS, DO_SSCORR_OVERALL, DO_ELASTIC_SCATTERING,      & ! Inputs
            DO_QUAD_RESULTS, DO_INCLUDE_MVOUT, FOURIER, FLUX_MULTIPLIER,   & ! Inputs
            NSTREAMS, N_USER_STREAMS, N_LOUTPUT, LEVELMASK_DN,             & ! Inputs
            WLOWER, L0_KTRANS(1,1,CPT), SAVE_TRANS_USERM(1,1,CPT),         & ! Inputs
            LCON, MCON, LCON_XVEC, MCON_XVEC,                              & ! Inputs
            U_XPOS, U_XNEG, U_WNEG1, U_WNEG2, HMULT_1, HMULT_2, EMULT_DN,  & ! Inputs
            CUMSOURCE_DN, ELASTIC_F_DN(1,1,APT), QUADELASTIC_F_DN(1,1,APT) ) ! Outputs
         ENDIF

!  Saved results for the Rayleigh planetary problem

!         IF ( DO_DNWELLING ) THEN
!           IF ( DO_RAYLEIGH_ONLY ) THEN
!             DO J = 1, N_LOUTPUT
!               DO K = 1, N_USER_STREAMS
!                IF(M.EQ.0)ELASTIC_F0_DN(J,K,APT)=ELASTIC_F_DN(J,K,APT)
!                IF(M.EQ.1)ELASTIC_F1_DN(J,K,APT)=ELASTIC_F_DN(J,K,APT)
!                IF(M.EQ.2)ELASTIC_F2_DN(J,K,APT)=ELASTIC_F_DN(J,K,APT)
!               ENDDO
!             ENDDO
!           ENDIF
!         ENDIF

!  Post processing of the Mean-Intensity & Flux Fields (Azimuth-integrated)

         IF ( DO_INCLUDE_MVOUT ) THEN
           CALL MIFLUX_INTENSITY_1 &
               ( DO_UPWELLING, DO_DNWELLING, COS_SZA, FLUX_FACTOR,      & ! Inputs
                 NSTREAMS, N_LOUTPUT, LEVELMASK_DN, BEAM_PICUTOFF(CPT), & ! Inputs
                 BEAM_ETRANS(1,CPT), QUAD_WEIGHTS, QUAD_STRMWGT,        & ! Inputs
                 QUADELASTIC_F_UP(1,1,APT), QUADELASTIC_F_DN(1,1,APT),  & ! Inputs
                 MEAN_ELASTIC_UP(1,APT), FLUX_ELASTIC_UP(1,APT),        & ! Outputs
                 MEAN_ELASTIC_DN(1,APT), FLUX_ELASTIC_DN(1,APT) )         ! Outputs
         ENDIF

!  End post processing clause for Intensity field

       ENDIF

!  #############################
!  Atmospheric Jacobian Section
!  #############################

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

!  Form L_WUPPER and L_WLOWER
!  --------------------------

!  Remark, should we zero this first ?

           DO N = 1, NLAYERS
             IF ( K.LE.N ) THEN
               NK = NKSTORAGE(N,K)
               DO Q = 1, KPARS(K)
                 DO I = 1, NSTR2
                   L_WUPPER(Q,I,N) = &
                    L_L0_WPARTIC(Q,I,NK,CPT) *   BEAM_ITRANS(N,CPT) + &
                      L0_WPARTIC(I,N,CPT)    * L_BEAM_ITRANS(Q,NK,CPT)
                   L_WLOWER(Q,I,N) = &
                    L_L0_WPARTIC(Q,I,NK,CPT) *   BEAM_ETRANS(N,CPT) + &
                      L0_WPARTIC(I,N,CPT)    * L_BEAM_ETRANS(Q,NK,CPT)
                 ENDDO
               ENDDO
             ELSE
               DO Q = 1, KPARS(K)
                 DO I = 1, NSTR2
                   L_WUPPER(Q,I,N) = ZERO
                   L_WLOWER(Q,I,N) = ZERO
                 ENDDO
               ENDDO
             ENDIF
           ENDDO

!  Linearized Boundary value problem
!  ---------------------------------

!  Linearized Beam surface reflectance term, for the last layer

           CALL L_SURFACE_BEAM_SETUP &
            ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, R2,                   & ! Inputs
              FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX, KVARY(K), & ! Inputs
              KPARS(K), QUAD_STRMWGT, ALBEDOS_RANKED, BRDF_F, L_WLOWER,  & ! Inputs
              L_R2_BEAM )                                                  ! Output

!  set up Column COL2WF for solution vector (the "B" as in AX=B)

           IF ( K .EQ. 0 ) THEN
             CALL LC_BVPCOLUMN_SETUP &
             ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, KPARS(K), NTOTAL, NSTR2,                     & ! Inputs
               LCON, MCON, XPOS, XNEG, L_XPOS, L_XNEG, L0_KTRANS(1,1,CPT), L_L0_KTRANS(1,1,1,CPT), & ! Inputs
               L_WUPPER, L_WLOWER, L_R2_HOMP, L_R2_HOMM, L_R2_BEAM, L_DIRECT_BEAM,                 & ! Inputs
               COL2WF, SCOL2WF )                                                                    ! Outputs
           ELSE
             CALL LPE_BVPCOLUMN_SETUP &
             ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, K, KPARS(K),         & ! Inputs
               NTOTAL, NSTR2, LCON, MCON, XPOS, XNEG, L0_KTRANS(1,1,CPT),  & ! Inputs
               L_XPOS, L_XNEG, L_L0_KTRANS(1,1,1,CPT), L_WUPPER, L_WLOWER, & ! Inputs
               L_R2_HOMP, L_R2_HOMM, L_R2_BEAM, L_DIRECT_BEAM,             & ! Inputs
               COL2WF, SCOL2WF)                                              ! Outputs
           ENDIF

!  LAPACK substitution using RHS column vector COL2
!    Single layer case (SCOL2). Version 2.5, 9/11/15
!    Exception handling modified for OMP usage, 10/10/18.

           IF ( NLAYERS .GT. 1 ) THEN

             CALL DGBTRS &
             ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, NPARS, &
                BANDMAT2, MAX_BANDTOTAL, IPIVOT, &
                COL2WF, MAX_TOTAL, INFO )

             IF ( INFO .LT. 0 ) THEN
               WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGBTRS: argument i illegal value, for i '//C2
               WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Elastic_Master Atmospheric Linearization, point number: '//C3
               FAIL = .TRUE. ; go to 69
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
!    Exception handling modified for OMP usage, 10/10/18.

             CALL DGETRS ( 'N', NTOTAL, NPARS, &
                      SMAT2, MAX_2_STREAMS, SIPIVOT, &
                      SCOL2WF, MAX_2_STREAMS, INFO )

             IF ( INFO .LT. 0 ) THEN
               WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGETRS: argument i illegal value, for i '//C2
               WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Elastic_Master Atmospheric Linearization, point number: '//C3
               FAIL = .TRUE. ; go to 69
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

!  Copy the above quantities into the output arrays

          DO Q = 1, KPARS(K)
             DO N = 1, NLAYERS
               NK = NKSTORAGE2(N,K)
               DO AA = 1, NSTREAMS
                 L_L0_WHOM_LCON(Q,AA,NK,CPT) = NCON(Q,AA,N)
                 L_L0_WHOM_MCON(Q,AA,NK,CPT) = PCON(Q,AA,N)
               ENDDO
             ENDDO
           ENDDO

!  Post processing for the Atmospheric Jacobians (elastic)
!  -------------------------------------------------------

           IF ( DO_POST_PROCESSING ) THEN

!  Atmospheric Linearization of User direct beam

             DO Q = 1, NPARS
               DO UM = 1, N_USER_STREAMS
                 L_USER_DIRECT_BEAM(Q,UM) = ZERO
               ENDDO
             ENDDO
             IF ( KVARY(K) .AND. DO_INCLUDE_SURFACE &
                           .AND. ATMOS_ATTN.NE.ZERO ) THEN
              NK = NKSTORAGE(NLAYERS,K)
              DO Q = 1, KPARS(K)
               L_AT = L_BEAM_ETRANS(Q,NK,CPT) / BEAM_ETRANS(NLAYERS,CPT)
               DO UM = 1, N_USER_STREAMS
                L_USER_DIRECT_BEAM(Q,UM) = L_AT*USER_DIRECT_BEAM(UM)
               ENDDO
              ENDDO
             ENDIF

!  Upwelling Jacobians
!   Version 2.5, upgrade for BRDF treatment and

             IF ( DO_UPWELLING ) THEN
               CALL ELASTIC_JACOBIAN_UP_1 &
          ( DO_MSMODE_LRRS, DO_SSCORR_OVERALL, DO_ELASTIC_SCATTERING,      & ! Inputs
            DO_QUAD_RESULTS, DO_INCLUDE_MVOUT, DO_REFLECTED_DIRECTBEAM,    & ! Inputs
            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, R2, FOURIER,              & ! Inputs
            NSTREAMS, NLAYERS, N_USER_STREAMS, N_LOUTPUT, Raman_Idx,       & ! Inputs
            Brdf_Idx, ALBEDOS_RANKED, USER_BRDF_F, FLUX_MULTIPLIER,        & ! Inputs
            QUAD_STRMWGT, LEVELMASK_UP, KVARY(K), K, KPARS(K), NKSTORAGE,  & ! Inputs
            L0_KTRANS(1,1,CPT), LOCAL_TRANS_USERM, LCON, MCON, XPOS, XNEG, & ! Inputs
            U_XPOS, U_XNEG, U_WPOS1, U_WPOS2, HMULT_1, HMULT_2, EMULT_UP,  & ! Inputs
            CUMSOURCE_UP, L_WUPPER, L_WLOWER, L_L0_KTRANS(1,1,1,CPT),      & ! Inputs
            L_USER_DIRECT_BEAM, L_LOCAL_TRANS_USERM, NCON, PCON,           & ! Inputs
            L_XPOS, L_XNEG, L_U_XPOS, L_U_XNEG,                            & ! Inputs
            L_U_WPOS1, L_U_WPOS2, L_HMULT_1, L_HMULT_2, L_EMULT_UP,        & ! Inputs
            L_ELASTIC_F_UP(1,0,1,1,APT), L_QUADELASTIC_F_UP(1,0,1,1,APT) )   ! Outputs
             ENDIF

!  Downwelling Jacobians

             IF ( DO_DNWELLING ) THEN
               CALL ELASTIC_JACOBIAN_DN_1 &
          ( DO_MSMODE_LRRS, DO_SSCORR_OVERALL, DO_ELASTIC_SCATTERING,    & ! Inputs
            DO_QUAD_RESULTS, DO_INCLUDE_MVOUT, FOURIER,                  & ! Inputs
            NSTREAMS, N_USER_STREAMS, N_LOUTPUT, FLUX_MULTIPLIER,        & ! Inputs
            LEVELMASK_DN, KVARY(K), K, KPARS(K), NKSTORAGE,              & ! Inputs
            L0_KTRANS(1,1,CPT), LOCAL_TRANS_USERM, LCON, MCON,           & ! Inputs
            XPOS, XNEG, U_XPOS, U_XNEG, U_WNEG1, U_WNEG2,                & ! Inputs
            HMULT_1, HMULT_2, EMULT_DN, CUMSOURCE_DN,                    & ! Inputs
            L_WLOWER, L_L0_KTRANS(1,1,1,CPT), L_LOCAL_TRANS_USERM,       & ! Inputs
            NCON, PCON, L_XPOS, L_XNEG, L_U_XPOS, L_U_XNEG,              & ! Inputs
            L_U_WNEG1, L_U_WNEG2, L_HMULT_1, L_HMULT_2, L_EMULT_DN,      & ! Inputs
            L_ELASTIC_F_DN(1,0,1,1,APT), L_QUADELASTIC_F_DN(1,0,1,1,APT) ) ! Outputs
             ENDIF

!  Linearization for Mean/Flux output

             IF ( DO_INCLUDE_MVOUT ) THEN
               CALL MIFLUX_JACOBIAN_1 &
            ( DO_UPWELLING, DO_DNWELLING, NSTREAMS, N_LOUTPUT, LEVELMASK_DN,          & ! Inputs
              KVARY(K), K, KPARS(K), NKSTORAGE, COS_SZA, FLUX_FACTOR,                 & ! Inputs
              QUAD_WEIGHTS, QUAD_STRMWGT, BEAM_PICUTOFF(CPT), L_BEAM_ETRANS(1,1,CPT), & ! Inputs
              L_QUADELASTIC_F_UP(1,0,1,1,APT), L_QUADELASTIC_F_DN(1,0,1,1,APT),       & ! Inputs
              L_MEAN_ELASTIC_UP(1,0,1,APT), L_FLUX_ELASTIC_UP(1,0,1,APT),             & ! Outputs
              L_MEAN_ELASTIC_DN(1,0,1,APT), L_FLUX_ELASTIC_DN(1,0,1,APT) )              ! Outputs
             ENDIF

!  End post processing

           ENDIF

!  End Major Atmospheric Jacobian loop

         ENDDO

!  End Atmospheric Jacobian clause

       ENDIF

!  ########################
!  Surface Jacobian Section
!  ########################

!  Restricted to the Lambertian albedo for now

       IF ( DO_SURFACE_WFS ) THEN

!  Surface Linearization of Direct beam reflected to quadrature streams
!     Treatment added 21 November 2008. Updated for Version 2.5, 12 October 0215

         LS_DIRECT_BEAM = ZERO
!mick question 7/20/2016 - add ".NOT.DO_SSCORR_OVERALL" condition here?
!                          (the LS_USER_DIRECT_BEAM if block has it)
         IF ( DO_INCLUDE_SURFACE ) THEN
           IF ( DO_BRDF_SURFACE ) THEN
             DO Q = 1, N_SURFACE_WFS
               DO I = 1, NSTREAMS
                 LS_DIRECT_BEAM(Q,I) = ATMOS_ATTN * LS_BRDF_F_0(Q,FOURIER,I,BRDF_Idx)
               ENDDO
             ENDDO
           ELSE
!mick fix 7/20/2016 - added Fourier condition
             IF ( FOURIER.eq.0 ) then
               DO Q = 1, N_SURFACE_WFS
                 DO I = 1, NSTREAMS
                   LS_DIRECT_BEAM(Q,I) = ATMOS_ATTN
                 ENDDO
               ENDDO
             ENDIF
           ENDIF
         ENDIF

!  Linearized BV Problem for surface Jacobians
!  -------------------------------------------

!  Set up Column COL2 for solution vector (the "B" as in AX=B)

         CALL LSE_BVPCOLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER,       & ! Inputs
             NLAYERS, NSTREAMS, NTOTAL, NSTR2, N_SURFACE_WFS,    & ! Inputs
             R2, Brdf_IDX, QUAD_STRMWGT, LS_BRDF_F,              & ! Inputs
             LCON, MCON, XPOS, XNEG, L0_KTRANS(1,1,CPT), WLOWER, & ! Inputs
             R2_HOMP, R2_HOMM, R2_BEAM, LS_DIRECT_BEAM,          & ! Inputs
             COL2WF_SURF, SCOL2WF_SURF )                           ! Outputs

!  General case: LAPACK substitution using RHS column vector COL2WF_SURF
!    Single layer case (SCOL2). Version 2.5, 9/11/15
!    Exception handling modified for OMP usage, 10/10/18.

         IF ( NLAYERS .GT. 1 ) THEN

           CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SURFACE_WFS, &
              BANDMAT2, MAX_BANDTOTAL, IPIVOT, &
              COL2WF_SURF, MAX_TOTAL, INFO )

           IF ( INFO .LT. 0 ) THEN
             WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGBTRS: argument i illegal value, for i = '//C2
             WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Elastic_Master, Surface Linearization, point number: '//C3
             FAIL = .TRUE. ; go to 69
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
!    Exception handling modified for OMP usage, 10/10/18.

           CALL DGETRS ( 'N', NTOTAL, N_SURFACE_WFS, &
                      SMAT2, MAX_2_STREAMS, SIPIVOT, &
                      SCOL2WF_SURF, MAX_2_STREAMS, INFO )

           IF ( INFO .LT. 0 ) THEN
             WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGETRS: argument i illegal value, for i = '//C2
             WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'L_Elastic_Master, Surface Linearization, point number: '//C3
             FAIL = .TRUE. ; go to 69
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

!  Copy the above quantities into the output arrays.
!    November 19th, 2008

         DO N = 1, NLAYERS
           DO Q = 1, N_SURFACE_WFS
             DO AA = 1, NSTREAMS
               LS_L0_WHOM_LCON(Q,AA,N,CPT) = NCON_SURF(Q,AA,N)
               LS_L0_WHOM_MCON(Q,AA,N,CPT) = PCON_SURF(Q,AA,N)
             ENDDO
           ENDDO
         ENDDO

!  Post processing for Surface Jacobians
!  -------------------------------------

         IF ( DO_POST_PROCESSING ) THEN

!  Surface Linearization of Direct beam reflected to User streams
!     Treatment added 21 November 2008. Upgraded for Version 2.5

           LS_USER_DIRECT_BEAM = ZERO
           IF ( .NOT.DO_SSCORR_OVERALL .AND. DO_INCLUDE_SURFACE ) THEN
             IF ( DO_BRDF_SURFACE ) THEN
               DO Q = 1, N_SURFACE_WFS
                 DO UI = 1, N_USER_STREAMS
                   LS_USER_DIRECT_BEAM(Q,UI) = ATMOS_ATTN * LS_USER_BRDF_F_0(Q,FOURIER,UI,BRDF_Idx)
                 ENDDO
               ENDDO
             ELSE
!mick fix 7/20/2016 - added Fourier condition
               IF ( FOURIER.eq.0 ) THEN
                 DO Q = 1, N_SURFACE_WFS
                   DO UI = 1, N_USER_STREAMS
                     LS_USER_DIRECT_BEAM(Q,UI) = ATMOS_ATTN
                   ENDDO
                 ENDDO
               ENDIF
             ENDIF
           ENDIF

!write(*,*)DO_SSCORR_OVERALL, DO_BRDF_SURFACE, FOURIER, LS_USER_DIRECT_BEAM(1,1:2)

!  Upwelling surface Jacobians

           IF ( DO_UPWELLING ) THEN
             CALL ELASTIC_SURFJAC_UP_1 &
             ( DO_SSCORR_OVERALL, DO_QUAD_RESULTS, DO_INCLUDE_MVOUT,              & ! Inputs
               DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, R2, ALBEDOS_RANKED, Raman_Idx,& ! Inputs
               FOURIER, NSTREAMS, NLAYERS, N_USER_STREAMS, N_LOUTPUT,             & ! Inputs
               N_SURFACE_WFS, FLUX_MULTIPLIER, LEVELMASK_UP,                      & ! Inputs
               Brdf_Idx, USER_BRDF_F, LS_USER_BRDF_F, LS_USER_DIRECT_BEAM,        & ! Inputs
               LOCAL_TRANS_USERM, QUAD_STRMWGT, XPOS, XNEG, L0_KTRANS(1,1,CPT),   & ! Inputs
               IDOWNSURF, U_XPOS, U_XNEG, HMULT_1, HMULT_2, NCON_SURF, PCON_SURF, & ! Inputs
               LS_ELASTIC_F_UP(1,1,1,APT), LS_QUADELASTIC_F_UP(1,1,1,APT) )         ! Outputs
           ENDIF

!  Downwelling surface Jacobians

           IF ( DO_DNWELLING ) THEN
             CALL ELASTIC_SURFJAC_DN_1 &
             ( DO_QUAD_RESULTS, DO_INCLUDE_MVOUT,                         & ! Inputs
               FOURIER, NSTREAMS, N_USER_STREAMS, N_LOUTPUT,              & ! Inputs
               N_SURFACE_WFS, FLUX_MULTIPLIER, LEVELMASK_DN,              & ! Inputs
               LOCAL_TRANS_USERM, XPOS, XNEG, L0_KTRANS(1,1,CPT),         & ! Inputs
               U_XPOS, U_XNEG, HMULT_1, HMULT_2, NCON_SURF, PCON_SURF,    & ! Inputs
               LS_ELASTIC_F_DN(1,1,1,APT), LS_QUADELASTIC_F_DN(1,1,1,APT) ) ! Outputs
           ENDIF

!  Surface Linearization for Mean/Flux output

           IF ( DO_INCLUDE_MVOUT ) THEN
             CALL MIFLUX_SURFJAC_1 &
             ( DO_UPWELLING, DO_DNWELLING, NSTREAMS, N_LOUTPUT,                & ! Inputs
               N_SURFACE_WFS, QUAD_WEIGHTS, QUAD_STRMWGT,                      & ! Inputs
               LS_QUADELASTIC_F_UP(1,1,1,APT), LS_QUADELASTIC_F_DN(1,1,1,APT), & ! Inputs
               LS_MEAN_ELASTIC_UP(1,1,APT), LS_FLUX_ELASTIC_UP(1,1,APT),       & ! Outputs
               LS_MEAN_ELASTIC_DN(1,1,APT), LS_FLUX_ELASTIC_DN(1,1,APT) )        ! Outputs
           ENDIF

!  End post processing for surface Jacobians

         ENDIF

!  End surface Jacobian clause

       ENDIF

!mick add 10/19/2015 - Sleave Jac section

!  ################################
!  Surface-leaving Jacobian Section
!  ################################

       IF ( DO_SLEAVE_WFS ) THEN

!  Main call

         CALL LRRS_LSSL_WFS ( &
           DO_REFLECTED_DIRECTBEAM, DO_INCLUDE_MVOUT,                       & ! Inputs
           DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_SL_ISOTROPIC,    & ! Inputs
           FOURIER, Raman_IDX, Brdf_IDX, Sleave_IDX,                        & ! Inputs
           NLAYERS, NSTREAMS, N_USER_STREAMS, N_LOUTPUT,                    & ! Inputs
           NSTR2, N_SUBDIAG, N_SUPDIAG, NTOTAL, N_SURFACE_WFS, N_SLEAVE_WFS,& ! Inputs
           LEVELMASK_UP, LEVELMASK_DN, ALBEDOS_RANKED, USER_BRDF_F,         & ! Inputs
           SURFACE_FACTOR, FLUX_MULTIPLIER, QUAD_WEIGHTS, QUAD_STRMWGT,     & ! Inputs
           BANDMAT2, IPIVOT, SMAT2, SIPIVOT,                                & ! Inputs
           L0_KTRANS(1,1,CPT), XPOS, XNEG, HMULT_1, HMULT_2,                & ! Inputs
           LOCAL_TRANS_USERM, U_XPOS, U_XNEG,                               & ! Inputs
           LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM,                         & ! Inputs
           LS_ELASTIC_F_UP(1,1,1,APT),  LS_ELASTIC_F_DN(1,1,1,APT),         & ! Outputs
           LS_MEAN_ELASTIC_UP(1,1,APT), LS_MEAN_ELASTIC_DN(1,1,APT),        & ! Outputs
           LS_FLUX_ELASTIC_UP(1,1,APT), LS_FLUX_ELASTIC_DN(1,1,APT),        & ! Outputs
           STATUS_SUB, MESSAGE_SUB, POINT_TRACE )                             ! Outputs

!  Exception handling modified for OMP usage, 10/10/18

         IF ( STATUS_SUB /= LRRS_SUCCESS ) THEN
           MESSAGE = 'L_Elastic_Master, LRRS_LSSL_WFS'
           WRITE(C3, '(I3)' ) CPT ; POINT_TRACE = 'point number: ' // C3 // TRIM(POINT_TRACE)
           FAIL = .TRUE. ; go to 69
         ENDIF

!  End Sleave WFs clause

       ENDIF

!  label 69 again defined at bottom of subroutine, as Continuation point Error

69     continue

!  Save OMP thread information

       LRRS_TID_SAVE(CPT) = TID + 1

!  End of main wavelength loop

      ENDDO

!$OMP END DO

!  End parallel region

!$OMP END PARALLEL

!     @@@@@@@@@@@@@@@@@@@
!     End parallel region
!     @@@@@@@@@@@@@@@@@@@

!  Time spent in OpenMP parallel region = OMP_Elastic_Time, to be compared with Ser_Elastic_Time

      if ( DO_TIMING ) then
        !write(*,*)
        !do CPT = 1, NPOINTS_LOCAL
        !  write(*,*) 'CPT = ',CPT,'LRRS_TID_SAVE(CPT) = ',LRRS_TID_SAVE(CPT)
        !enddo
        call cpu_time(omp_e2) ; OMP_Elastic_Time = (omp_e2 - omp_e1)/REAL(NTHREADS)
      endif

!  Finish

      RETURN
      END SUBROUTINE L_ELASTIC_MASTER_1

!  End Module

      END MODULE lrrs_L_elastic_master_1_m

