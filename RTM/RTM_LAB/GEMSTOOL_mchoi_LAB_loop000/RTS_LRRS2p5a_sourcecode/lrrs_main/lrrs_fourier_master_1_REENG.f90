
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
! #             RAMAN_FOURIER_1 (*)                             #
! #                                                             #
! #   Note (*): Version 2.1 with BRDF setup,   November 2007    #
! #   Note    : Version 2.2 with Linearization, November 2008   #
! #   Note    : Version 2.3 in F90,              March 2011     #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Use of Supplement-derived BRDF/SLEAVE inputs and control
!    (2) Bookkeeping improvements (use of "Only", clearer I/O specifications)
!    (3) Use of new Taylor series inputs and modules.

      MODULE lrrs_fourier_master_1_m

!      USE LRRS_PARS_m, Only : SDU
!   -- Rob mod 5/12/17 for 2p5a, SSFULL renamed, SSCORR_NADIR added

      PRIVATE
      PUBLIC :: RAMAN_FOURIER_1

      CONTAINS

!Rob/Mick mod -- 9/17. Add OMP Inputs and Outputs

      SUBROUTINE RAMAN_FOURIER_1 ( DO_TIMING, NTHREADS,                       & ! OMP Inputs
         DO_RRS_OVERALL, DO_DIRECTRRS_ONLY, DO_MSMODE_LRRS, PROGRESS,         & ! Inputs
         DO_ELASTIC_ONLY, DO_SSCORR_GENERAL, DO_SSCORR_ALONE,                 & ! Inputs
         DO_BIN_REALIZATION, DO_ENERGY_BALANCING,                             & ! Inputs
         DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_DIRECT_BEAM,                  & ! Inputs
         DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_BRDF_SURFACE,        & ! Inputs
         DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, do_BRDF_Wav1, do_SLEAVE_Wav1,   & ! Inputs
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
         ELASTIC_F_UP, MEAN_ELASTIC_UP, FLUX_ELASTIC_UP,                      & ! Outputs
         ELASTIC_F_DN, MEAN_ELASTIC_DN, FLUX_ELASTIC_DN,                      & ! Outputs
         RAMAN_F_UP, MEAN_RAMAN_UP, FLUX_RAMAN_UP,                            & ! Outputs
         RAMAN_F_DN, MEAN_RAMAN_DN, FLUX_RAMAN_DN,                            & ! Outputs
         FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE,                             & ! Outputs
         OMP_Elastic_Time, OMP_Raman_Time )                                     ! OMP Timing Outputs

!  Module of dimensions and numbers

      USE LRRS_PARS_m

!  Use Modules
!   -- Rob mod 5/12/17 for 2p5a, introduced LRRS_MISCSETUPS_1_m, reorganized LRRS_RTSOLUTIONS_1_m

      USE LRRS_RTSOLUTIONS_1_m        , Only : HOMOG_SOLUTION, UHOMOG_SOLUTION
      USE LRRS_MISCSETUPS_1_m         , Only : DBEAM_SETUP, UDBEAM_SETUP, LEGENDRE_SETUPS
      USE LRRS_BVPROBLEM_m            , Only : SURFACE_HOMOG_SETUP, SURFACE_BEAM_SETUP,  &
                                               BVPMATRIX_INIT, BVPMATRIX_SETUP, BVPCOLUMN_SETUP
      USE LRRS_AUX2_m                 , Only : DGBTRF, DGBTRS, DGETRS, DGETRF
      USE LRRS_POSTPROCESSING_1_m     , Only : HOMOGMULT_1, BEAMMULT_UP_1, BEAMMULT_DN_1, BOASOURCE
      USE LRRS_ELASTIC_MASTER_1_m     , Only : ELASTIC_MASTER_1
      USE LRRS_SOURCES_MASTER_1_m     , Only : SOURCES_MASTER_1
      USE LRRS_RAMAN_INTENSITY_1_m    , Only : MIFLUX_INTENSITY_1, RAMAN_INTENSITY_UP_1,  RAMAN_INTENSITY_DN_1

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Number of OMP Cores to be used (1, 2, 4, 8...)

      INTEGER  , intent(in) ::  NTHREADS

!  Monitoring flag

      logical  , intent(in) :: DO_TIMING

!  Progress number
!    If this is zero, no output to screen

      INTEGER  , INTENT(IN) :: PROGRESS

!  elastic control flag

      LOGICAL  , INTENT(IN) :: DO_ELASTIC_ONLY

!  Raman flags for computing filling

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

!  Raman output at User angles, for one Fourier component
!  ------------------------------------------------------

!  Fourier output

      REAL(FPK), INTENT(OUT) :: RAMAN_F_UP ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: RAMAN_F_DN ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  Mean-value output

      REAL(FPK), INTENT(INOUT) :: MEAN_RAMAN_UP( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: MEAN_RAMAN_DN( MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: FLUX_RAMAN_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK), INTENT(INOUT) :: FLUX_RAMAN_DN ( MAX_LOUTPUT, MAX_POINTS )

!  module status:
!   Message_sub comes from any of the subroutines called here.
!   Point_trace gives a trace of subroutine where message occurred

      LOGICAL            , INTENT(OUT) :: FAIL
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE_SUB
      CHARACTER (LEN=120), INTENT(OUT) :: POINT_TRACE

!  Mick mod 9/7/17. Add OMP timing outputs

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

      REAL(FPK) :: RAMAN_WUPPER  ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK) :: RAMAN_WLOWER  ( MAX_2_STREAMS, MAX_LAYERS )

!  Reflected beam solutions at ground

      REAL(FPK) :: R2_BEAM(MAX_STREAMS)

!  Boundary value problem. AX = B
!  ------------------------------

!  Compressed band matrix or Regular matrix  as in A

      REAL(FPK) :: BANDMAT2(MAX_BANDTOTAL,MAX_TOTAL)
      REAL(FPK) :: SMAT2   (MAX_2_STREAMS, MAX_2_STREAMS)

!  compression indexing for the band matrix

      INTEGER   :: BMAT_ROWMASK(MAX_TOTAL,MAX_TOTAL)

!  RHS column vectors, as in B

      REAL(FPK) :: COL2  ( MAX_TOTAL,1 )
      REAL(FPK) :: SCOL2 ( MAX_2_STREAMS,1 )

!  Pivot for the LU decomposition of the band matrix

      INTEGER   :: IPIVOT  ( MAX_TOTAL )
      INTEGER   :: SIPIVOT (MAX_2_STREAMS)

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

      REAL(FPK) :: U_HELP_P ( MAX_STREAMS, 0:MAX_MOMENTS, MAX_LAYERS )
      REAL(FPK) :: U_HELP_M ( MAX_STREAMS, 0:MAX_MOMENTS, MAX_LAYERS )

!  Miscellaneous
!  -------------

!  Fourier output (quadrature streams)
!    Use as debug, except for Flux calculations

      REAL(FPK) :: QUADRAMAN_F_UP ( MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )
      REAL(FPK) :: QUADRAMAN_F_DN ( MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

!  Log-file output (LU = 6 gives screen output)

      LOGICAL  , PARAMETER :: LG = .false.
      INTEGER  , PARAMETER :: LU = 8

!  for exception handling

      CHARACTER (LEN=3) ::  C3
      CHARACTER (LEN=2) ::  C2

!  Local integers

      INTEGER   :: NTOTAL, N_SUPDIAG, N_SUBDIAG, N_GAIN_TERMS
      INTEGER   :: NSTR2, N_SOLUTIONS, PROG, M, INFO
      INTEGER   :: LAYER, I, I1, AA, N, C0, APT, CPT, L, UM, W
      INTEGER   :: NPOINTS_TOTAL, NPOINTS_LOCAL
      INTEGER   :: Raman_IDX, Brdf_IDX, Sleave_IDX

!  Number of RRS solutions

      INTEGER   :: N_LOSS_SOLUTIONS
      INTEGER   :: N_GAIN_SOLUTIONS

!  Monochromatic realizationresults_o3bin_M_tester_CR_SZA45_lamb.resg

      LOGICAL   :: DO_MONO_REALIZATION

!  Cabannes-Raman flag

      LOGICAL   :: DO_CABANNES_RAMAN

!  Flags for including surface term

      LOGICAL   :: DO_INCLUDE_SURFACE
      LOGICAL   :: DO_REFLECTED_DIRECTBEAM

!  Flag for mean value output

      LOGICAL   :: DO_INCLUDE_MVOUT

!  Flag for post-processing of elastic field

      LOGICAL   :: DO_ELASTIC_POSTPROCESSING

!  Local flag for CB COPY

      LOGICAL   :: DO_CBCOPY

!  Surface factors and local XT, DT

      REAL(FPK) :: SURFACE_FACTOR, SURFACE_FACTOR_0
      REAL(FPK) :: DELTA_FACTOR, HOM, NORM, T1, T2

!  solar irradianace flux

      REAL(FPK) :: SS_FLUX_MULTIPLIER
      REAL(FPK) :: FLUX_MULTIPLIER

!  Local IOPs

      REAL(FPK) :: LOCAL_OMEGAMOMS ( MAX_LAYERS, 0: MAX_MOMENTS )

!  Overall flag for computing filling

      LOGICAL   :: DO_FILLING

!  Rob Fix 26 March 2012. Add new variable DO_AVOID_UDBEAresults_o3bin_M_tester_CR_SZA45_lamb.resgM

      LOGICAL   :: DO_AVOID_UDBEAM

!  Debug

      LOGICAL   :: ELASTIC_CALL=.FALSE.

!  OMP Variables. Added September 2017
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

!  Added September 2017. Set total number of threads (from NTHREADS input)
!     write(*,*) 'Setting total number of threads'

      OMP_MAXTHREADS = NTHREADS
      CALL OMP_SET_NUM_THREADS(OMP_MAXTHREADS)

!  set exclusion flags

      DO_MONO_REALIZATION = .not. DO_BIN_REALIZATION
      DO_CABANNES_RAMAN   = .not. DO_ENERGY_BALANCING

!  set overall filling flags

      DO_FILLING = .not.DO_ELASTIC_ONLY.and.DO_RRS_OVERALL

!  Fourier

      M = FOURIER

!mick fix 7/20/2016 - this inelastic-related if block moved after the
!                     return from the elastic-only computation
!  Monochromatic - all points
!  Binned        - all points in outer buffer

!      IF ( DO_BIN_REALIZATION ) THEN
!        NPOINTS_TOTAL = NPOINTS_OUTER
!        NPOINTS_LOCAL = NPOINTS_INNER
!      ELSE
!        NPOINTS_TOTAL = NPOINTS_MONO
!        NPOINTS_LOCAL = 1
!      ENDIF

!  Rob Fix 26 March 2012
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

!  This is the coherent light solution with 'Rayleigh' scattering
!  Same for both Cabannes/Raman and Energy Balancing methods

!  Rob Fix 26 March 2012 
!   - Replace DO_SSCORR_GENERAL with DO_AVOID_UDBEAM, in Line # 2
!       of subroutine arguments to ELASTIC_MASTER_1. Formerly --
!      CALL ELASTIC_MASTER_1 &
!       ( DO_ELASTIC_ONLY, DO_ELASTIC_POSTPROCESSING, DO_MSMODE_LRRS, &
!         DO_INCLUDE_MVOUT, DO_INCLUDE_SURFACE, DO_SSCORR_GENERAL, &
!         ................ etc...

!  Add SLTERM flag and Isotropic SLTERM. @@@ RobFix 06 sep 12

!  Added September 2017. OMP inputs and timing output

      CALL ELASTIC_MASTER_1 ( DO_TIMING, NTHREADS,                            & ! OMP Input
         DO_ELASTIC_ONLY, DO_ELASTIC_POSTPROCESSING, DO_MSMODE_LRRS,          & ! Inputs
         DO_INCLUDE_MVOUT, DO_INCLUDE_SURFACE, DO_AVOID_UDBEAM,               & ! Inputs
         DO_UPWELLING, DO_DNWELLING, DO_BIN_REALIZATION, DO_BRDF_SURFACE,     & ! Inputs
         DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, do_BRDF_Wav1, do_SLEAVE_Wav1,   & ! Inputs
         ALBEDOS_RANKED, BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,        & ! Inputs
         SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                       & ! Inputs
         NLAYERS, NSTREAMS, NMOMENTS, N_LOUTPUT, N_USER_STREAMS,              & ! Inputs
         FOURIER, COS_SZA, TAYLOR_SMALL, TAYLOR_ORDER, FLUX_FACTOR,           & ! Inputs
         NPOINTS_OUTER, OFFSET_INNER, NPOINTS_INNER, NPOINTS_MONO, W_EXCIT,   & ! Inputs
         QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS,              & ! Inputs
         STERM_MASK_UP, LEVELMASK_UP, STERM_MASK_DN, LEVELMASK_DN,            & ! Inputs
         DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, SAVE_TRANS_USERM,              & ! Inputs
         BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS, BEAM_ETRANS, & ! Inputs
         PLM_PXI, PLM_MXI, PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI, PLM_MXUI,        & ! Inputs
         PLM_00_PXI, PLM_00_MXI, PLM_00_PXUI, PLM_00_MXUI,                    & ! Inputs
         L0_KEIGEN, L0_KTRANS, L0_WPARTIC,                                    & ! Outputs
         L0_WHOM_XPOS,  L0_WHOM_LCON,  L0_WHOM_MCON,                          & ! Outputs
         ELASTIC_F_UP, MEAN_ELASTIC_UP, FLUX_ELASTIC_UP,                      & ! Outputs
         ELASTIC_F_DN, MEAN_ELASTIC_DN, FLUX_ELASTIC_DN,                      & ! Outputs
         FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE, OMP_Elastic_Time )            ! Outputs

!write(*,*)'Pert1',fourier,ELASTIC_F_UP(1:3,1,1)
!write(*,*)'Pert1',fourier,ELASTIC_F_UP(1:3,1,6) ; stop

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
!mick fix 9/7/2017 - added M, SURFACE_FACTOR_0

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
!$OMP     NPOINTS_LOCAL, NSTR2, NTOTAL, N_SUBDIAG, N_SUPDIAG, DELTA_FACTOR,      & ! Inputs (Internal)
!$OMP     FOURIER, COS_SZA, TAYLOR_SMALL, TAYLOR_ORDER, FLUX_FACTOR,             & ! Inputs
!$OMP     FLUX_MULTIPLIER, FLUXES_RANKED, DO_REFLECTED_DIRECTBEAM,               & ! Inputs (Internal)
!$OMP     M, SURFACE_FACTOR_0,                                                   & ! Inputs (Internal)
!$OMP     OFFSET_INNER, NPOINTS_INNER, W_EXCIT, NPOINTS_MONO, N_RRSBINS, BINMAP, & ! Inputs
!$OMP     QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS,                & ! Inputs
!$OMP     STERM_MASK_UP, LEVELMASK_UP, STERM_MASK_DN, LEVELMASK_DN,              & ! Inputs
!$OMP     DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, OMEGAMOMS_CABANNES,              & ! Inputs
!$OMP     OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN, OMEGAMOMS_RRSGAIN,                & ! Inputs
!$OMP     BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS, BEAM_ETRANS,   & ! Inputs
!$OMP     PLM_PXI, PLM_MXI, PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI, PLM_MXUI,          & ! Inputs
!$OMP     PLM_00_PXI, PLM_00_MXI, PLM_00_PXUI, PLM_00_MXUI, SAVE_TRANS_USERM,    & ! Inputs
!$OMP     L0_KEIGEN, L0_KTRANS, L0_WPARTIC,                                      & ! Inputs
!$OMP     L0_WHOM_XPOS,  L0_WHOM_LCON,  L0_WHOM_MCON,                            & ! Inputs
!$OMP     RAMAN_F_UP, MEAN_RAMAN_UP, FLUX_RAMAN_UP,                              & ! Outputs
!$OMP     RAMAN_F_DN, MEAN_RAMAN_DN, FLUX_RAMAN_DN,                              & ! Outputs
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
      write(*,'(1x,a,i1,a)') 'Running Raman_Fourier with ',OMP_NTHREADS,' thread(s)'
   END IF

!  Start points loop
!  =================

!$OMP DO

      DO W = 1, NPOINTS_LOCAL

!  For binning, set counters in inner (APT) and outer (CPT) loops
!    Use Local variable W_EXCIT_0 inside Parallel region (OMP). 9/7/17.

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

!  Local IOPs (omega-moments)
!    -Usual elastic-scattering IOPS, unless using Cabannes/Raman method

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

!  error handling. Modified for OMP usage 9/7/17

          IF ( FAIL_SUB ) THEN
            WRITE(C2,'(I2)')N   ; MESSAGE     = 'Raman_Fourier master, HOMOG_SOLUTION, layer'//C2
            WRITE(C3,'(I3)')CPT ; POINT_TRACE = 'point number: '//C3
            FAIL = .true. ; GOTO 69
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
                ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS,   & ! Inputs
                  LOCAL_OMEGAMOMS, XPOS, XNEG,                      & ! Inputs
                  PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI,                 & ! Inputs
                  U_XPOS, U_XNEG, U_HELP_P(1,0,N), U_HELP_M(1,0,N) )  ! Output

!  Whole layer multipliers
!   @@@ RobFix 5/5/11. Small numbers analysis: added LRRS_FGSMALL,
!                      DELTAU_VERT_INPUT(N,CPT) to Argument list.

          CALL HOMOGMULT_1 &
                ( NSTREAMS, N, N_USER_STREAMS,               & ! Inputs
                  TAYLOR_SMALL, TAYLOR_ORDER,                & ! Inputs
                  LOCAL_TRANS_USERM, USER_STREAMS,           & ! Inputs
                  KEIGEN, KTRANS, DELTAU_VERT_INPUT(N,CPT),  & ! Inputs
                  HMULT_1, HMULT_2, ZETA_P, ZETA_M )           ! Outputs

!  End layer loop

        ENDDO

!  Boundary Value Matrix
!  =====================

!  Additional setups for the albedo layer

        if(lg)write(lu,'(a)')'starting SURFACE_HOMOG_SETUP'
        CALL SURFACE_HOMOG_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR, & ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX,     & ! Inputs
           QUAD_STRMWGT, ALBEDOS_RANKED, BRDF_F, XPOS, XNEG,    & ! Inputs
           R2_HOMP, R2_HOMM )                                     ! Outputs

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

!  (Error tracing). Modified for OMP usage 9/7/17

        IF ( INFO .GT. 0 ) THEN
          WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGBTRF: Singular matrix, u(i,i)=0, for i = '//C2
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(C2, '(I2)' ) INFO ; MESSAGE = 'DGBTRF: argument i illegal value, for i = '//C2
        ENDIF

        IF ( INFO .ne. 0 ) THEN
          WRITE(C3,'(I3)')CPT ; POINT_TRACE = 'Raman_Fourier Radiances, point number: '//C3
          FAIL = .true. ; GOTO 69
        ENDIF

!  SVD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

!  New for Version 2.5

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK LU-decomposition for  matrix

        CALL DGETRF ( NTOTAL, NTOTAL, SMAT2, MAX_2_STREAMS, SIPIVOT, INFO )

!  (Error tracing). Modified for OMP usage 9/7/17

        IF ( INFO .GT. 0 ) THEN
          WRITE(C2, '(I2)' ) INFO ; MESSAGE     = 'DGBTRF: argument i illegal value, for i = '//C2
          WRITE(C3, '(I3)' ) CPT  ; POINT_TRACE = 'Raman_Fourier Radiances, point number: '//C3
          FAIL = .true. ; GOTO 69
        ENDIF

      ENDIF

!  Find particular solutions for elastic and all RRS terms
!  =======================================================

!  Get source vectors for wavelength point APT
!   -- Rob mod 5/12/17 for 2p5a, Use DO_SSCORR general flag (OUTGOING or NADIR)

        if ( lg ) write(lu, '(a)' )'starting SOURCES_MASTER_1'
        CALL SOURCES_MASTER_1 &
       ( FOURIER, APT, DO_ENERGY_BALANCING,                                & ! Inputs
         DO_RRS_OVERALL, DO_DIRECTRRS_ONLY, DO_MSMODE_LRRS,                & ! Inputs
         DO_BIN_REALIZATION, DO_MONO_REALIZATION, DO_SSCORR_GENERAL,       & ! Inputs
         DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                      & ! Inputs
         N_RRSBINS, BINMAP, NPOINTS_MONO, OFFSET_INNER, W_EXCIT_0,           & ! Inputs
         NLAYERS, NSTREAMS, NMOMENTS, N_USER_STREAMS,                      & ! Inputs
         QUAD_WEIGHTS, USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER,           & ! Inputs
         STERM_MASK_UP, STERM_MASK_DN, FLUXES_RANKED,                      & ! Inputs
         DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, OMEGAMOMS_CABANNES,         & ! Inputs
         OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN, OMEGAMOMS_RRSGAIN,           & ! Inputs
         BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT,                        & ! Inputs
         BEAM_DTRANS, SAVE_TRANS_USERM, PLM_PXI, PLM_MXI,                  & ! Inputs
         PLM_PXUI, PLM_MXUI, PLM_00_PXI, PLM_00_MXI,                       & ! Inputs
         PLM_00_PXUI, PLM_00_MXUI, PLM_WT_PXI,  PLM_WT_MXI,                & ! Inputs
         L0_KEIGEN,     L0_KTRANS,     L0_WPARTIC,                         & ! Inputs
         L0_WHOM_XPOS,  L0_WHOM_LCON,  L0_WHOM_MCON,                       & ! Inputs
         XPOS, KEIGEN, KTRANS, U_XPOS, U_XNEG,                             & ! Inputs
         EIGENNORM_SAVE, HMULT_1, HMULT_2,                                 & ! Inputs
         RAMAN_WUPPER, RAMAN_WLOWER, WLAYER_PIST_UP, WLAYER_PIST_DN,       & ! Output
         FAIL_SUB, MESSAGE_SUB )                                             ! Output

!  exception handling. Modified for OMP usage 9/7/17

        IF ( FAIL_SUB ) THEN
          MESSAGE = 'Error from SOURCES_MASTER_1'
          IF ( DO_BIN_REALIZATION ) THEN
            WRITE(C3,'(I3)')CPT ; POINT_TRACE = 'Raman_Fourier Radiances, point number: '//C3
          ENDIF
          FAIL = .true. ; GOTO 69
        ENDIF

!  ****************
!  Radiation Field:
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
            WRITE(C3, '(I3)' ) INFO ; MESSAGE_SUB = 'argument i illegal value, for i = '//C3
            MESSAGE     = 'RAMAN_FOURIER_1, DGBTRS call'
            IF ( DO_BIN_REALIZATION ) THEN
              WRITE(C3,'(I3)')CPT ; POINT_TRACE = 'Raman_Fourier Radiances, point number: '//C3
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
!    Exception handling modifield for OMP 9/7/17.

          CALL DGETRS ( 'N', NTOTAL, 1, &
                      SMAT2, MAX_2_STREAMS, SIPIVOT, &
                      SCOL2, MAX_2_STREAMS, INFO )

          IF ( INFO .LT. 0 ) THEN
            WRITE(C3, '(I3)' ) INFO ; MESSAGE_SUB = 'argument i illegal value, for i = '//C3
            MESSAGE     = 'RAMAN_FOURIER_1, DGETRS call (1 layer)'
            IF ( DO_BIN_REALIZATION ) THEN
              WRITE(C3,'(I3)')CPT ; POINT_TRACE = 'Raman_Fourier Radiances, point number: '//C3
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

!  label 69 again defined at bottom of subroutine, as Continuation point Error

69     continue

!  Save OMP thread information

       LRRS_TID_SAVE(w) = TID + 1

!write(*,*)'Pert1',fourier,APT,RAMAN_F_UP(1:3,1,APT)

!  Finish point loop

      ENDDO

!$OMP END DO

!  End parallel region

!$OMP END PARALLEL

!     @@@@@@@@@@@@@@@@@@@
!     End parallel region
!     @@@@@@@@@@@@@@@@@@@

!      pause 'End Fourier 1, reg'

!  Time spent in OpenMP parallel region = OMP_Elastic_Time, to be compared with Ser_Elastic_Time

      if (DO_TIMING) then
        call cpu_time(omp_e2) ; OMP_Raman_Time = (omp_e2 - omp_e1)/REAL(NTHREADS)
      endif

!  Finish
!  ======

      RETURN
      END SUBROUTINE RAMAN_FOURIER_1


      END MODULE lrrs_fourier_master_1_m

