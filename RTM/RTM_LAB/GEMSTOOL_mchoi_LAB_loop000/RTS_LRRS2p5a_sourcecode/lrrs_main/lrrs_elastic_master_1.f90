
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
! #             ELASTIC_MASTER  (master)                        #
! #                                                             #
! #  This module develops complete discrete ordinate elastic    #
! #    scattering solutions, with stored output as required     #
! #    later on for the RRS calculations, and also post-        #
! #    processed elastic-scattering radiances (if flagged)      #
! #                                                             #
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
! #   Note (*): Version 2.1 with BRDF setup, November 2007.     #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Use of Supplement-derived BRDF/SLEAVE inputs and control
!    (2) Bookkeeping improvements (use of "Only", clearer I/O specifications)
!    (3) Use of new Taylor series inputs and modules.

      MODULE lrrs_elastic_master_1_m

!      USE LRRS_PARS_m, Only : SDU

      PRIVATE
      PUBLIC :: ELASTIC_MASTER_1

      CONTAINS

      SUBROUTINE ELASTIC_MASTER_1 &
       ( DO_ELASTIC_ONLY, DO_ELASTIC_POSTPROCESSING, DO_MSMODE_LRRS,          & ! Inputs
         DO_INCLUDE_MVOUT, DO_INCLUDE_SURFACE, DO_SSCORR_OVERALL,             & ! Inputs
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
         FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE )                              ! Outputs

!  Module of dimensions and numbers

      USE LRRS_PARS_m
!   -- Rob mod 5/12/17 for 2p5a, introduced LRRS_MISCSETUPS_1_m, reorganized LRRS_RTSOLUTIONS_1_m

!  Use Modules

      USE LRRS_MISCSETUPS_1_m         , Only : DBEAM_SETUP, UDBEAM_SETUP
      USE LRRS_RTSOLUTIONS_1_m        , Only : HOMOG_SOLUTION,  BEAM_SOLUTION_MATRIX, BEAM_SOLUTION_ELASTIC, &
                                               UHOMOG_SOLUTION, UBEAM_SOLUTION_ELASTIC
      USE LRRS_BVPROBLEM_m            , Only : SURFACE_HOMOG_SETUP, SURFACE_BEAM_SETUP,  &
                                               BVPMATRIX_INIT, BVPMATRIX_SETUP, BVPCOLUMN_SETUP
      USE LRRS_AUX2_m                 , Only : DGBTRF, DGBTRS, DGETRS, DGETRF
      USE LRRS_POSTPROCESSING_1_m     , Only : HOMOGMULT_1, BEAMMULT_UP_1, BEAMMULT_DN_1
      USE LRRS_ELASTIC_INTENSITY_1_m  , Only : ELASTIC_INTENSITY_UP_1, ELASTIC_INTENSITY_DN_1
      USE LRRS_RAMAN_INTENSITY_1_m    , Only : MIFLUX_INTENSITY_1

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Top level control
!  -----------------

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

      REAL(fpk) :: BRDF_F_0      ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk) :: BRDF_F        ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk) :: USER_BRDF_F_0 ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )
      REAL(fpk) :: USER_BRDF_F   ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

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

!  module status
!  -------------

      LOGICAL            , INTENT(OUT) :: FAIL
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE_SUB
      CHARACTER (LEN=120), INTENT(OUT) :: POINT_TRACE

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

!  Help arrays (for use with linearization)

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

      REAL(FPK) :: QMAT  ( MAX_STREAMS, MAX_STREAMS )
      INTEGER   :: QPIVOT( MAX_STREAMS )

!  Classical beam solution matrices and vectors
!  Beam solution independent of optical depth (classical solution)
!  Help arrays for linearization

      REAL(FPK) :: QSUMVEC ( MAX_STREAMS )
      REAL(FPK) :: QDIFVEC ( MAX_STREAMS )
      REAL(FPK) :: QVEC    ( MAX_STREAMS )
      REAL(FPK) :: QAUX    ( MAX_STREAMS )

!  Particular solution help variable

      REAL(FPK) :: W_HELP ( 0:MAX_MOMENTS )

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

!  Integration constants, as in X

      REAL(FPK) :: LCON ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK) :: MCON ( MAX_STREAMS, MAX_LAYERS )

!  Integration constants multiplied by homogeneous solution vectors

      REAL(FPK) :: &
            LCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            MCON_XVEC ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  post processing solutions (for user-defined polar directions)
!  -------------------------

!  Direct beam transmittance

      REAL(FPK) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS )

!  Homogeneous solutions

      REAL(FPK) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Particular integrals, diffuse term (1), primary scatter term (2)

      REAL(FPK) :: &
            U_WPOS1 ( MAX_USER_STREAMS, MAX_LAYERS ), &
            U_WPOS2 ( MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: &
            U_WNEG1 ( MAX_USER_STREAMS, MAX_LAYERS ), &
            U_WNEG2 ( MAX_USER_STREAMS, MAX_LAYERS )

!  Post processing Multipliers
!  ---------------------------

!  Whole layer multipliers homogeneous solutions

      REAL(FPK) :: &
            HMULT_1 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            HMULT_2 ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: &
            ZETA_P ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS ), &
            ZETA_M ( MAX_STREAMS, MAX_USER_STREAMS, MAX_LAYERS )

!  Whole layer solutions, particular integrals

      REAL(FPK) :: EMULT_UP ( MAX_USER_STREAMS, MAX_LAYERS )
      REAL(FPK) :: EMULT_DN ( MAX_USER_STREAMS, MAX_LAYERS )

!  Cumulative radiance source terms

      REAL(FPK) :: CUMSOURCE_UP ( MAX_USER_STREAMS, 0:MAX_LAYERS )
      REAL(FPK) :: CUMSOURCE_DN ( MAX_USER_STREAMS, 0:MAX_LAYERS )

!  Output diffuse downwelling field term
!     Wrong dimensioning, bug JvG found March 2010

!      REAL(FPK) :: IDOWNSURF ( MAX_USER_STREAMS )
      REAL(FPK) :: IDOWNSURF ( MAX_STREAMS )

!  Atmospheric attenuation

      REAL(FPK) :: ATMOS_ATTN

!  Fourier output (quadrature streams)
!    Use as debug, except for Flux calculations

      REAL(FPK) :: QUADELASTIC_F_UP &
              ( MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

      REAL(FPK) :: QUADELASTIC_F_DN &
              ( MAX_LOUTPUT, MAX_STREAMS, MAX_POINTS )

!  Other local help variables
!  --------------------------

      INTEGER   :: I, I1, AA, N, NC, L, C0, M, UM
      INTEGER   :: INFO, CPT, APT, NPOINTS_LOCAL
      INTEGER   :: NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTR2
      INTEGER   :: Raman_IDX, Brdf_IDX, Sleave_IDX
      LOGICAL   :: DO_REFLECTED_DIRECTBEAM
      LOGICAL   :: DO_QUAD_RESULTS, NOFLIP

      REAL(FPK) :: DELFAC, R1, R2
      REAL(FPK) :: FLUX_MULTIPLIER, AT

      LOGICAL   :: DO_ELASTIC_SCATTERING(MAX_LAYERS)

      LOGICAL   :: DO_POST_PROCESSING
      LOGICAL   :: DO_CALCULATION

      CHARACTER (LEN=2) :: C2
      CHARACTER (LEN=3) :: C3

!  Debug

      LOGICAL   :: ELASTIC_CALL=.TRUE.

!  initialise
!  ----------

!  set status

      FAIL            = .FALSE.
      MESSAGE         = ' '
      POINT_TRACE     = ' '

!  Local flags

      DO_QUAD_RESULTS = .FALSE.

!  Monochromatic - all points
!  Binned        - all points in outer buffer

      IF ( DO_BIN_REALIZATION ) THEN
        NPOINTS_LOCAL = NPOINTS_OUTER
      ELSE
        NPOINTS_LOCAL = NPOINTS_MONO
      ENDIF

!  Fourier

      M = FOURIER

!  Start for the actual number of points to be calculated

      APT  = 0

!  set up the elastic scattering flag
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

!  setup INTEGERs for the BV Problem

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

      DELFAC = TWO
      IF ( FOURIER .EQ. 0 ) DELFAC = ONE
      FLUX_MULTIPLIER       = FLUX_FACTOR * DELFAC / PI4

!  Create BRDF Fourier component stuff here
!...... COMMENTED OUt for Version 2.5, 9/11/15
!      IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
!        CALL LRRS_BRDF_FOURIER &
!          ( DO_INCLUDE_SURFACE, DO_UPWELLING, DO_USER_STREAMS, &
!            DO_SSCORR_OVERALL, NSTREAMS, N_USER_STREAMS, QUAD_STRMWGT, &
!            N_BRDF, X_BRDF, A_BRDF, FOURIER, DELFAC, &
!            BRDFUNC, BRDFUNC_0, USER_BRDFUNC, USER_BRDFUNC_0, &
!            BIREFLEC, BIREFLEC_0, USER_BIREFLEC, USER_BIREFLEC_0 )
!      ENDIF

!  Surface reflectance factor R2

      R1 = ONE
      R2 = R1
      IF ( FOURIER .EQ. 0 ) R2 = TWO * R1

!  Direct beam flag (only if above surface flag has been set)

      DO_REFLECTED_DIRECTBEAM = .FALSE.
      IF ( DO_INCLUDE_SURFACE ) THEN
         DO_REFLECTED_DIRECTBEAM = .TRUE.
      ENDIF

!  Start points loop
!  =================

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

!  Direct beam reflected to quadrature directions, assumes Flux factor 1

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

!  ##################################
!  RT differential equation solutions
!  ##################################

!  homogeneous solutions
!  ---------------------

!  Start layer loop

       DO N = 1, NLAYERS

!  Find the solution, by solving the eignproblem.
!    - output eigenvalues L0_KEIGEN, eigentransmittances L0_KTRANS
!       these are required for the RRS calculation.
!    -  output saved eigenmatrices and homogeneous soltuion vectors

        CALL HOMOG_SOLUTION &
            ( N, FOURIER, NSTREAMS, NMOMENTS,                       & ! Input
              OMEGAMOMS_ELASTIC(1,0,CPT), DELTAU_VERT_INPUT(1,CPT), & ! Input
              QUAD_STREAMS, QUAD_WEIGHTS, PLM_PXI, PLM_MXI,         & ! Input
              XPOS, XNEG, L0_KEIGEN(1,1,CPT), L0_KTRANS(1,1,CPT),   & ! Output
              EIGENMAT, EIGENVEC, DIFVEC, DAB, SAB,                 & ! Output
              FAIL, MESSAGE_SUB )                                     ! Output

!  error handling

        IF ( FAIL ) THEN
          WRITE(C2,'(I2)')N
          MESSAGE = 'Elastic master, HOMOG_SOLUTION, layer'//C2
          WRITE(C3, '(I3)' ) CPT
          POINT_TRACE = 'point number: '//C3
          RETURN
        ENDIF

!  Store elastic solutions for later use
!   Save space by storing XPOS, LCON, MCON. April 10th, 2007

        DO AA = 1, NSTREAMS
          DO I = 1, NSTR2
            L0_WHOM_XPOS(I,AA,N,CPT) = XPOS(I,AA,N)
          ENDDO
        ENDDO

!  Elastic Beam Solution
!  =====================

!  Beam solution Matrix and LU-decomposition
!  -----------------------------------------

!  Linear Algebra setup, AX = B. Creates A in LU-decomposed form.
!     ALso create the Linearized Matrix................

        CALL BEAM_SOLUTION_MATRIX &
              ( N, NSTREAMS, EIGENMAT,                    & ! Input
                BEAM_AVSECANT(N,CPT), BEAM_PICUTOFF(CPT), & ! Input
                QMAT, QPIVOT, FAIL, MESSAGE_SUB )           ! Output

        IF ( FAIL ) THEN
          WRITE(C2,'(I2)')N
          MESSAGE = 'L_Elastic master, BEAM_SOLUTION_MATRIX, layer'//C2
          WRITE(C3, '(I3)' ) CPT
          POINT_TRACE = 'point number: '//C3
          RETURN
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
              L0_WPARTIC(1,1,CPT), FAIL, MESSAGE_SUB )          ! Output

!  Exception handling for this solution

        IF ( FAIL ) THEN
          WRITE(C2,'(I2)')N
          MESSAGE = 'Elastic master, BEAM_SOLUTION_ELASTIC, layer'//C2
          WRITE(C3, '(I3)' ) CPT
          POINT_TRACE = 'point number: '//C3
          RETURN
        ENDIF

!  Compute particular integrals at the layer boundaries

        DO I = 1, NSTR2
          WUPPER(I,N) = L0_WPARTIC(I,N,CPT) * BEAM_ITRANS(N,CPT)
          WLOWER(I,N) = L0_WPARTIC(I,N,CPT) * BEAM_ETRANS(N,CPT)
        ENDDO

!  End layer loop for all solutions

       ENDDO

!  ######################
!  Boundary Value Problem
!  ######################

!  Additional setups for the lowest layer with surface.
!  Version 2.5, 9/11/15, New setup for the BRDF, using supplement-derived material.

       CALL SURFACE_HOMOG_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, R2,            & ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX,    & ! Inputs
           QUAD_STRMWGT, ALBEDOS_RANKED, BRDF_F, XPOS, XNEG,   & ! Inputs
           R2_HOMP, R2_HOMM )                                    ! Outputs

!  initialize compression matrix

       CALL BVPMATRIX_INIT &
           ( NSTREAMS, NLAYERS,                   & ! Input
             NTOTAL, NSTR2, N_SUBDIAG, N_SUPDIAG, & ! Input
             BMAT_ROWMASK, BANDMAT2 )               ! Output

!  set up boundary values matrix in compressed form (the "A" as in AX=B)

       CALL BVPMATRIX_SETUP &
           (  DO_INCLUDE_SURFACE,                               & ! Inputs
              NSTREAMS, NLAYERS, NSTR2, BMAT_ROWMASK,           & ! Inputs
              XPOS, XNEG, R2_HOMP, R2_HOMM, L0_KTRANS(1,1,CPT), & ! Inputs
              BANDMAT2, SMAT2 )                                   ! Output

!  SVD the BVP matrix: With compression (multilayers)
!  --------------------------------------------------

!  LAPACK LU-decomposition for band matrix

       IF (NLAYERS .GT. 1 ) THEN

        CALL DGBTRF &
           ( NTOTAL, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
             BANDMAT2, MAX_BANDTOTAL, IPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(C2, '(I2)' ) INFO
          MESSAGE = 'DGBTRF: Singular matrix, u(i,i)=0, for i = '//C2
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(C2, '(I2)' ) INFO
          MESSAGE = 'DGBTRF: argument i illegal value, for i = '//C2
        ENDIF
        IF ( INFO .NE. 0 ) THEN
          FAIL = .TRUE.
          WRITE(C3, '(I3)' ) CPT
          POINT_TRACE = 'Elastic Radiances, point number: '//C3
          RETURN
        ENDIF

!  SVD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

!  New for Version 2.5

       ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK LU-decomposition for  matrix

        CALL DGETRF ( NTOTAL, NTOTAL, SMAT2, MAX_2_STREAMS, SIPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(C2, '(I2)' ) INFO
          MESSAGE = 'DGETRF: argument i illegal value, for i = '//C2
          FAIL = .TRUE.
          WRITE(C3, '(I3)' ) CPT
          POINT_TRACE = 'Elastic Radiances, point number: '//C3
          RETURN
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

       IF ( NLAYERS .GT. 1 ) THEN

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAX_BANDTOTAL, IPIVOT, &
              COL2, MAX_TOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          FAIL = .TRUE.
          WRITE(C2, '(I2)' ) INFO
          MESSAGE = 'DGBTRS: argument i illegal value, for i = '//C2
          WRITE(C3, '(I3)' ) CPT
          POINT_TRACE = 'Elastic Radiances, point number: '//C3
          RETURN
        ENDIF

!  Integration constants to be stored for the RRS calculation
!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

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

         CALL DGETRS ( 'N', NTOTAL, 1, &
                      SMAT2, MAX_2_STREAMS, SIPIVOT, &
                      SCOL2, MAX_2_STREAMS, INFO )

         IF ( INFO .LT. 0 ) THEN
           FAIL = .TRUE.
           WRITE(C2, '(I2)' ) INFO
           MESSAGE = 'DGETRS: argument i illegal value, for i = '//C2
           WRITE(C3, '(I3)' ) CPT
           POINT_TRACE = 'Elastic Radiances, point number: '//C3
           RETURN
         ENDIF

!  Integration constants to be stored for the RRS calculation
!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

         DO I = 1, NSTREAMS
           I1 = I+NSTREAMS
           LCON(I,1) = SCOL2(I,1)
           MCON(I,1) = SCOL2(I1,1)
         ENDDO

       ENDIF

!  Quadrature solutions: multiplied by integration constants
!   ( Only required for post processing )

       DO N = 1, NLAYERS
         DO I = 1, NSTR2
           DO AA = 1, NSTREAMS
             LCON_XVEC(I,AA,N) = LCON(AA,N) * XPOS(I,AA,N)
             MCON_XVEC(I,AA,N) = MCON(AA,N) * XNEG(I,AA,N)
           ENDDO
         ENDDO
       ENDDO

!  Copy the above quantities into the output arrays L0_WHOMOG
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
!            L0_WHOMOG_DN(I,AA,N,CPT) = LCON_XVEC(I,AA,N)
!            L0_WHOMOG_UP(I,AA,N,CPT) = MCON_XVEC(I,AA,N)
!          ENDDO
!        ENDDO
!       ENDDO

!  Make  discrete ordinate field. Debug.
!    --Check functioning of multipliers. 4 December 2006. OK
!       DO I = 1, NSTREAMS
!         HELP = WUPPER(I+NSTREAMS,1)
!         DO AA = 1, NSTREAMS
!           HELP = HELP + LCON_XVEC(I+NSTREAMS,AA,1) +
!     &           L0_KTRANS(AA,1,cpt) * MCON_XVEC(I+NSTREAMS,AA,1)
!         ENDDO
!         if (cpt.eq.115)write(*,*)'up',cpt,I,quad_streams(i),help/pi4
!         if(i.eq.6)write(15,*)'up',cpt,I,quad_streams(i),help/pi4
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

         IF ( DO_ELASTIC_ONLY ) THEN
           DO_CALCULATION = .FALSE.
           IF ( DO_BIN_REALIZATION ) THEN
             IF ( CPT.GT.OFFSET_INNER .AND. &
               CPT.LT.OFFSET_INNER + 1 + NPOINTS_INNER ) THEN
               DO_CALCULATION = .TRUE.
             ENDIF
           ELSE
             IF ( CPT .EQ. W_EXCIT ) THEN
               DO_CALCULATION = .TRUE.
             ENDIF
           ENDIF
           IF ( CPT .EQ. W_EXCIT ) THEN
              DO_CALCULATION = .TRUE.
           ENDIF
         ELSE
           DO_CALCULATION = .FALSE.
           IF ( DO_BIN_REALIZATION ) THEN
             IF ( CPT.GT.OFFSET_INNER .AND. &
               CPT.LT.OFFSET_INNER + 1 + NPOINTS_INNER ) THEN
               DO_CALCULATION = .TRUE.
             ENDIF
           ELSE
             IF ( CPT .EQ. W_EXCIT ) THEN
               DO_CALCULATION = .TRUE.
             ENDIF
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

!  Post-processed Homogeneous/Particular solutions + All linearizations
!  ====================================================================

!  Individual layer loop

         DO N = 1, NLAYERS

!  homogeneous solutions at user-defined streams

           CALL UHOMOG_SOLUTION &
                ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS, & ! Inputs
                  OMEGAMOMS_ELASTIC(1,0,CPT), XPOS, XNEG,         & ! Inputs
                  PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI,               & ! Inputs
                  U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )              ! Output

!  user classical solutions

           CALL UBEAM_SOLUTION_ELASTIC &
            ( DO_UPWELLING, DO_DNWELLING, FOURIER,                    & ! Inputs
              N, STERM_MASK_UP, STERM_MASK_DN,                        & ! Inputs
              NSTREAMS, NMOMENTS, N_USER_STREAMS, BEAM_PICUTOFF(CPT), & ! Inputs
              L0_WPARTIC(1,1,CPT), OMEGAMOMS_ELASTIC(1,0,CPT),        & ! Inputs
              PLM_WT_PXI, PLM_WT_MXI,  PLM_00_PXUI,                   & ! Inputs
              PLM_PXUI,   PLM_00_MXUI, PLM_MXUI,                      & ! Inputs
              W_HELP, U_WPOS1, U_WPOS2, U_WNEG1, U_WNEG2 )              ! Outputs

!  End layer loop

         ENDDO

!  Post processing Multipliers
!  ===========================

!  layer loop

         DO N = 1, NLAYERS

!  Homogeneous solution multipliers
!   @@@ RobFix 5/5/11. Small numbers analysis: added LRRS_FGSMALL,
!                      DELTAU_VERT_INPUT(N,CPT) to Argument list.
!  Version 2.5, 9/11/15. Expanded Taylor-series control

           CALL HOMOGMULT_1 &
                ( NSTREAMS, N, N_USER_STREAMS,               & ! Inputs
                  TAYLOR_SMALL, TAYLOR_ORDER,                & ! Inputs
                  SAVE_TRANS_USERM(1,1,CPT), USER_STREAMS,   & ! Inputs
                  L0_KEIGEN(1,1,CPT), L0_KTRANS(1,1,CPT),    & ! Inputs
                  DELTAU_VERT_INPUT(N,CPT),                  & ! Inputs
                  HMULT_1, HMULT_2, ZETA_P, ZETA_M )           ! Outputs

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
           ENDIF

!  End layer loop

         ENDDO

!  Apply Initial Transmittances, Upwelling solutions
!     01 June 2009, This is necessary here.
!       * EMULT_UP nor U_WPOS1/U_WPOS2 not multiplied by ITRANS

         IF ( DO_UPWELLING ) THEN
           DO N = 1, NLAYERS
             AT = BEAM_ITRANS(N,CPT)
             DO UM = 1, N_USER_STREAMS
               EMULT_UP(UM,N) = AT * EMULT_UP(UM,N)
             ENDDO
           ENDDO
         ENDIF

!  Apply Initial Transmittances, Downwelling solutions
!     01 June 2009, This is necessary here.
!       * EMULT_DN nor U_WNEG1/U_WNEG2 not multiplied by ITRANS

         IF ( DO_DNWELLING ) THEN
           DO N = 1, NLAYERS
             AT = BEAM_ITRANS(N,CPT)
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

         APT = APT + 1

!  Post processing of the Upwelling intensity
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

!  Post processing of the Downwelling intensity

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

!  Miflux

         IF ( DO_INCLUDE_MVOUT ) THEN
           CALL MIFLUX_INTENSITY_1 &
               ( DO_UPWELLING, DO_DNWELLING, COS_SZA, FLUX_FACTOR,      & ! Inputs
                 NSTREAMS, N_LOUTPUT, LEVELMASK_DN, BEAM_PICUTOFF(CPT), & ! Inputs
                 BEAM_ETRANS(1,CPT), QUAD_WEIGHTS, QUAD_STRMWGT,        & ! Inputs
                 QUADELASTIC_F_UP(1,1,APT), QUADELASTIC_F_DN(1,1,APT),  & ! Inputs
                 MEAN_ELASTIC_UP(1,APT), FLUX_ELASTIC_UP(1,APT),        & ! Outputs
                 MEAN_ELASTIC_DN(1,APT), FLUX_ELASTIC_DN(1,APT) )         ! Outputs
         ENDIF

!  debug, 13 March 2008
!            write(*,*)'elastic up TOA',
!     &      apt, MEAN_ELASTIC_UP(1,APT),FLUX_ELASTIC_UP(1,APT)
!            write(*,*)'elastic up BOA',
!     &      apt, MEAN_ELASTIC_UP(2,APT),FLUX_ELASTIC_UP(2,APT)
!            write(*,*)'elastic dn TOA',
!     &      apt, MEAN_ELASTIC_DN(1,APT),FLUX_ELASTIC_DN(1,APT)
!            write(*,*)'elastic dn BOA',
!     &      apt, MEAN_ELASTIC_DN(2,APT),FLUX_ELASTIC_DN(2,APT)

!  End post processing clause

       ENDIF

!  Finish points loop

      ENDDO

!      pause
!  Finish

      RETURN
      END SUBROUTINE ELASTIC_MASTER_1

!  End Module

      END MODULE lrrs_elastic_master_1_m

