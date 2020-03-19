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
! #             ELASTIC_MASTER  (master)                        #
! #                                                             #
! #  This module develops complete discrete ordinate elastic    #
! #    scattering solutions, with stored output as required     #
! #    later on for the RRS calculations, and also post-        #
! #    processed elastic-scattering radiances (if flagged)      #
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
! #   Note (*): Version 2.1 with BRDF setup, November 2007.     #
! #                                                             #
! ###############################################################


      MODULE lrrs_elastic_master_1

      PRIVATE
      PUBLIC :: ELASTIC_MASTER_1

      CONTAINS

      SUBROUTINE ELASTIC_MASTER_1 &
       ( DO_ELASTIC_ONLY, DO_ELASTIC_POSTPROCESSING, DO_MSMODE_LRRS, &
         DO_INCLUDE_MVOUT, DO_INCLUDE_SURFACE, DO_SSCORR_OVERALL, &
         DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_BIN_REALIZ, &
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

!  include files
!  -------------

!  include file of dimensions and numbers

      USE LRRS_PARS
      USE LRRS_BRDF
      USE LRRS_RTSOLUTIONS_1
      USE LRRS_BVPROBLEM
      USE LRRS_AUX2
      USE LRRS_POSTPROCESSING_1
      USE LRRS_ELASTIC_INTENSITY_1
      USE LRRS_RAMAN_INTENSITY_1

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Top level control
!  -----------------

!  elastic control flags

      LOGICAL, INTENT(IN) ::          DO_ELASTIC_ONLY
      LOGICAL, INTENT(IN) ::          DO_ELASTIC_POSTPROCESSING

!  Multiple scatter mode

      LOGICAL, INTENT(IN) ::          DO_MSMODE_LRRS

!  Output for Flux calculations. Set in the MASTER

      LOGICAL, INTENT(IN) ::          DO_INCLUDE_MVOUT

!  Flag for including surface. set in the MASTER

      LOGICAL, INTENT(IN) ::          DO_INCLUDE_SURFACE

!  Overall SS correction flag; Now incorporates the DB term

      LOGICAL, INTENT(IN) ::          DO_SSCORR_OVERALL

!  Upwelling and downwelling flags

      LOGICAL, INTENT(IN) ::          DO_UPWELLING
      LOGICAL, INTENT(IN) ::          DO_DNWELLING

!  User stream flag

      LOGICAL, INTENT(IN) ::          DO_USER_STREAMS

!  Binning realization

      LOGICAL, INTENT(IN) ::          DO_BIN_REALIZ

!  Surface control
!  ---------------

!  overall lambertian surface

      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_SURFACE

!  Surface contributions

      REAL(FPK), INTENT(IN) :: FL1 ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: FL2 ( MAX_POINTS )

!  BRDF control and quadrature

      INTEGER, INTENT(IN) ::          N_BRDF
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

      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          N_USER_STREAMS

!  Number of moments

      INTEGER, INTENT(IN) ::          NMOMENTS

!  number of layers

      INTEGER, INTENT(IN) ::          NLAYERS

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)

      INTEGER, INTENT(IN) ::          N_LOUTPUT

!  Fourier number

      INTEGER, INTENT(IN) ::          FOURIER

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

      INTEGER, INTENT(IN) ::          NPOINTS_OUTER
      INTEGER, INTENT(IN) ::          NPOINTS_INNER
      INTEGER, INTENT(IN) ::          OFFSET_INNER
      INTEGER, INTENT(IN) ::          NPOINTS_MONO
      INTEGER, INTENT(IN) ::          W_EXCIT

!  Layer masks for doing integrated source terms
!  ---------------------------------------------

      LOGICAL, INTENT(IN) ::          STERM_MASK_UP ( MAX_LAYERS )
      LOGICAL, INTENT(IN) ::          STERM_MASK_DN ( MAX_LAYERS )

      INTEGER, INTENT(IN) ::          LEVELMASK_UP    (MAX_LOUTPUT)
      INTEGER, INTENT(IN) ::          LEVELMASK_DN    (MAX_LOUTPUT)

!  Beam Attenuations
!  -----------------

!  Solar beam transmittances, average secant factors

      INTEGER, INTENT(IN) ::          BEAM_PICUTOFF ( MAX_POINTS )
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

      REAL(FPK), INTENT(OUT) :: L0_KEIGEN &
              ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Eigentransmittances

      REAL(FPK), INTENT(OUT) :: L0_KTRANS &
              ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Particular solution vectors

      REAL(FPK), INTENT(OUT) :: L0_WPARTIC &
            ( MAX_2_STREAMS, MAX_LAYERS, MAX_POINTS )

!  homogeneous solution vectors and integration constants

      REAL(FPK), INTENT(OUT) :: L0_WHOM_XPOS &
            ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L0_WHOM_LCON &
            ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L0_WHOM_MCON &
            ( MAX_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Fourier components of BRDF functions
!  ------------------------------------

!  at quadrature (discrete ordinate) angles

      REAL(FPK), INTENT(OUT) :: BIREFLEC   ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: BIREFLEC_0 ( MAX_STREAMS  )

!  at user-defined stream directions

      REAL(FPK), INTENT(OUT) :: USER_BIREFLEC ( MAX_USER_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: USER_BIREFLEC_0 ( MAX_USER_STREAMS  )

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

!  module status
!  -------------

      LOGICAL, INTENT(OUT) ::          FAIL
      CHARACTER (LEN=70), INTENT(OUT) ::     MESSAGE
      CHARACTER (LEN=70), INTENT(OUT) ::     MESSAGE_SUB
      CHARACTER (LEN=70), INTENT(OUT) ::     POINT_TRACE

!***********************************************************************
!  ARRAYS FOR THE DETERMINATION OF THE ELASTIC DISCRETE ORDINATE SOLUTIO
!     These are based on the LIDORT model with the classical solution
!***********************************************************************

!  Local setups
!  ------------

!  local transmittances along user-stream directions

      REAL(FPK) :: &
          LOCAL_TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS )

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
!  at the user-taus partial layer values (WUTAUS)

      REAL(FPK) :: WUPPER ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )

!  Reflected solutions at ground (direct beam, diffusely reflected beam)

      REAL(FPK) :: DIRECT_BEAM ( MAX_STREAMS )
      REAL(FPK) :: R2_BEAM     ( MAX_STREAMS )

!  Particular solution: LU-decomposed matrix and Pivot

      REAL(FPK) :: QMAT ( MAX_STREAMS, MAX_STREAMS )
      INTEGER ::          QPIVOT    ( MAX_STREAMS )

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

!  Compressed band matrix, as in A

      REAL(FPK) :: BANDMAT2(MAX_BANDTOTAL,MAX_TOTAL)

!  compression indexing for the band matrix

      INTEGER ::          BMAT_ROWMASK(MAX_TOTAL,MAX_TOTAL)

!  RHS column vector, as in B

      REAL(FPK) :: COL2 ( MAX_TOTAL,1 )

!  Pivot for the LU decomposition of the band matrix

      INTEGER ::          IPIVOT ( MAX_TOTAL )

!  Integration constants, as in X

      REAL(FPK) :: &
            LCON ( MAX_STREAMS, MAX_LAYERS ), &
            MCON ( MAX_STREAMS, MAX_LAYERS )

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

      INTEGER ::          I, I1, AA, N, NC, L, C0, M, UM
      INTEGER ::          INFO, CPT, APT, NPOINTS_LOCAL
      INTEGER ::          NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTR2

      LOGICAL ::          DO_REFLECTED_DIRECTBEAM
      LOGICAL ::          DO_QUAD_RESULTS, NOFLIP

      REAL(FPK) :: DELFAC, R1, R2, DETER, V1, V2
      REAL(FPK) :: FLUX_MULTIPLIER, AT

      LOGICAL ::          DO_ELASTIC_SCATTERING(MAX_LAYERS)

      LOGICAL ::          DO_POST_PROCESSING
      LOGICAL ::          DO_CALCULATION

      CHARACTER (LEN=2) ::      C2
      CHARACTER (LEN=3) ::      C3

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

      IF ( DO_BIN_REALIZ ) THEN
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

!  debug

      !  DO N = 1, NLAYERS
      !    DO L = FOURIER, 2*NSTREAMS-1
      !      DO CPT = 1, NPOINTS_LOCAL
      !        write(88,*)n,l,cpt,OMEGAMOMS_ELASTIC(N,L,CPT)
      !      ENDDO
      !    ENDDO
      !  ENDDO
      !  pause

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

      IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
        CALL LRRS_BRDF_FOURIER &
          ( DO_INCLUDE_SURFACE, DO_UPWELLING, DO_USER_STREAMS, &
            DO_SSCORR_OVERALL, NSTREAMS, N_USER_STREAMS, QUAD_STRMWGT, &
            N_BRDF, X_BRDF, A_BRDF, FOURIER, DELFAC, &
            BRDFUNC, BRDFUNC_0, USER_BRDFUNC, USER_BRDFUNC_0, &
            BIREFLEC, BIREFLEC_0, USER_BIREFLEC, USER_BIREFLEC_0 )
      ENDIF

!  albedo reflectance factor R2

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

!  Copy local trans userm

       DO UM = 1, N_USER_STREAMS
        DO N = 1, NLAYERS
          LOCAL_TRANS_USERM(UM,N) = SAVE_TRANS_USERM(UM,N,CPT)
        ENDDO
       ENDDO

!  Direct beam reflected to quadrature directions, assumes Flux factor 1
!     BRDF treatment added 2 November 2007.

       CALL DBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, &
             NSTREAMS, FL1(CPT), FL2(CPT), BIREFLEC_0, &
             COS_SZA, FOURIER, DELFAC, BEAM_ETRANS(NLAYERS,CPT), &
             ATMOS_ATTN, DIRECT_BEAM )

!         if ( cpt.eq.4.and.FOURIER.eq.0) &
!           write(*,*)'Calculation',DIRECT_BEAM(1)

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
!       these are required for the RRS calculation, and are stored
!       in common blocks in ELASTIC_SAVE.VARS.
!    -  output saved eigenmatrices and homogeneous soltuion vectors

        CALL HOMOG_SOLUTION &
            ( N, FOURIER, NSTREAMS, NMOMENTS, &
              OMEGAMOMS_ELASTIC(1,0,CPT), DELTAU_VERT_INPUT(1,CPT), &
              QUAD_STREAMS, QUAD_WEIGHTS, &
              PLM_PXI,     PLM_MXI, &
              XPOS, XNEG, L0_KEIGEN(1,1,CPT), L0_KTRANS(1,1,CPT), &
              EIGENMAT, EIGENVEC, DIFVEC, DAB, SAB, &
              FAIL, MESSAGE_SUB )

!  error handling

        IF ( FAIL ) THEN
          WRITE(C2,'(I2)')N
          MESSAGE = 'Elastic master, HOMOGENEOUS_SOLUTION, layer'//C2
          WRITE(C3, '(I3)' ) CPT
          POINT_TRACE = 'point number: '//C3
          RETURN
        ENDIF

!  Store elastic solutions for later use
!   These arrays are stored in commons in ELASTICSAVE.VARS
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
              ( N, NSTREAMS, EIGENMAT, &
                BEAM_AVSECANT(N,CPT), BEAM_PICUTOFF(CPT), &
                QMAT, QPIVOT, FAIL, MESSAGE_SUB )
        IF ( FAIL ) THEN
          WRITE(C2,'(I2)')N
          MESSAGE = &
              'L_Elastic master, BEAM_SOLUTION_MATRIX, layer'//C2
          WRITE(C3, '(I3)' ) CPT
          POINT_TRACE = 'point number: '//C3
          RETURN
        ENDIF

!  Beam solution solver (Back-substitution)
!  ----------------------------------------

!   - output = particular integrals at layer boundaries (WUPPER/WLOWER)
!   - output = L0_WPARTIC, the solution source vector for RRS scattering
!         This latter quantity is stored in commons (ELASTIC_SAVE.VARS)

        CALL BEAM_SOLUTION_ELASTIC &
            ( N, FOURIER,  NMOMENTS, NSTREAMS, &
              BEAM_AVSECANT(N,CPT), BEAM_ITRANS(N,CPT), &
              QMAT, QPIVOT, SAB, DAB, &
              PLM_00_PXI,  PLM_00_MXI, &
              QUAD_STREAMS, OMEGAMOMS_ELASTIC(1,0,CPT), &
              QSUMVEC, QDIFVEC, QVEC, QAUX, &
              L0_WPARTIC(1,1,CPT), FAIL, MESSAGE_SUB )

!  Exception handling for this solution

        IF ( FAIL ) THEN
          WRITE(C2,'(I2)')N
          MESSAGE = 'Elastic master, BEAM_SOLUTION_ELASTIC, layer'//C2
          WRITE(C3, '(I3)' ) CPT
          POINT_TRACE = 'point number: '//C3
          RETURN
        ENDIF

!  COmpute Particular integrals at the layer boundaries

        DO I = 1, NSTR2 !=2 * NSTREAMS
          WUPPER(I,N) = L0_WPARTIC(I,N,CPT) * BEAM_ITRANS(N,CPT)
          WLOWER(I,N) = L0_WPARTIC(I,N,CPT) * BEAM_ETRANS(N,CPT)
        ENDDO

!  End layer loop for all solutions

       ENDDO

!      if ( fourier.eq.0)
!     &  call wupperlower_write ( max_layers, max_2_streams,
!     &        nlayers, nstreams, cpt, 1, 13, wupper, wlower )

!  ######################
!  Boundary Value Problem
!  ######################

!  Additional setups for the lowest layer with surface.

       CALL SURFACE_HOMOG_SETUP &
         ( DO_INCLUDE_SURFACE, &
           DO_LAMBERTIAN_SURFACE, &
           R2, FL1(CPT), FL2(CPT), BIREFLEC, &
           FOURIER, NLAYERS, NSTREAMS, &
           QUAD_STRMWGT, XPOS, XNEG, &
           R2_HOMP, R2_HOMM )

!  initialize compression matrix

       CALL BVPMATRIX_INIT &
           ( NSTREAMS, NLAYERS, &
             NTOTAL, NSTR2, N_SUBDIAG, N_SUPDIAG, &
             BMAT_ROWMASK, BANDMAT2 )

!  set up boundary values matrix in compressed form (the "A" as in AX=B)

       CALL BVPMATRIX_SETUP &
           (  DO_INCLUDE_SURFACE, &
              NSTREAMS, NLAYERS, NSTR2, BMAT_ROWMASK, &
              XPOS, XNEG, R2_HOMP, R2_HOMM, L0_KTRANS(1,1,CPT), &
              BANDMAT2 )

!  LAPACK LU-decomposition for band matrix

       IF (NTOTAL.GT. 2 ) THEN

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
          POINT_TRACE = 'point number: '//C3
          RETURN
        ENDIF

       ENDIF

!  More debug, 18 November. Problem solved !!!!!

!        if ( fourier.le.1)  then
!         do n = nlayers-10, nlayers
!           write(*,*)cpt,n,(wupper(i,n),i=1,3)
!         enddo
!        endif

!  Solution of boundary value problem
!  ----------------------------------

!  Add contributions to the BV problem RHS, for the lowest layer

       CALL SURFACE_BEAM_SETUP &
         ( DO_INCLUDE_SURFACE, &
           DO_LAMBERTIAN_SURFACE, &
           R2, FL1(CPT), FL2(CPT), BIREFLEC, &
           FOURIER, NLAYERS, NSTREAMS, &
           QUAD_STRMWGT, WLOWER, &
           R2_BEAM )

!       write(35,*)cpt, fourier
!       do i = 1, nstreams
!         write(35,*)i,r2_beam(i),r2_homm(i,1),r2_homp(i,2)
!       enddo
!       if ( cpt.eq.20.and.fourier.eq.0)pause

!  set up Column COL2 for solution vector (the "B" as in AX=B)

       CALL BVPCOLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, &
             NSTREAMS, NLAYERS, NTOTAL, NSTR2, &
             DIRECT_BEAM, WUPPER, WLOWER, R2_BEAM, &
             COL2 )

!        do i = 1, ntotal
!           if ( cpt.eq.1)write(45,*)i,col2(i,1)
!        enddo

!  LAPACK substitution using RHS column vector COL2

       IF ( NTOTAL.GT. 2 ) THEN

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

!  special case, 1 layer, 1 stream

       ELSE

        DETER = BANDMAT2(1,1) * BANDMAT2(2,2) - &
                BANDMAT2(1,2) * BANDMAT2(2,1)
        V1 = COL2(1,1)
        V2 = COL2(2,1)
        COL2(1,1) = ( BANDMAT2(2,2)*V1 - BANDMAT2(1,2)*V2 ) / DETER
        COL2(2,1) = ( BANDMAT2(1,1)*V2 - BANDMAT2(2,1)*V1 ) / DETER

       ENDIF

!  Integration constants to be stored for the RRS calculation
!  ----------------------------------------------------------

!  Set integration constants LCON and MCON for -/+ eigensolutions, all l

       DO N = 1, NLAYERS
        C0 = (N-1)*NSTR2
        DO I = 1, NSTREAMS
          I1 = I+NSTREAMS
          LCON(I,N) = COL2(C0+I,1)
          MCON(I,N) = COL2(C0+I1,1)
        ENDDO
!        if ( cpt.eq.1) write(*,*) N, LCON(1:2,N) 
       ENDDO

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
           DO_CALCULATION = .TRUE.
         ELSE
           DO_CALCULATION = .FALSE.
           IF ( DO_BIN_REALIZ ) THEN
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

         CALL UDBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_SSCORR_OVERALL, &
             DO_LAMBERTIAN_SURFACE, FL1(CPT), FL2(CPT), USER_BIREFLEC_0, &
             FOURIER, N_USER_STREAMS, ATMOS_ATTN, &
             USER_DIRECT_BEAM )

!  Post-processed Homogeneous/Particular solutions + All linearizations
!  ====================================================================

!  Individual layer loop

         DO N = 1, NLAYERS

!  homogeneous solutions at user-defined streams

           CALL UHOMOG_SOLUTION &
                ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS, &
                  OMEGAMOMS_ELASTIC(1,0,CPT), XPOS, XNEG, &
                  PLM_WT_PXI, PLM_WT_MXI, PLM_PXUI, &
                  U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )

!  user classical solutions

           CALL UBEAM_SOLUTION_ELASTIC &
            ( DO_UPWELLING, DO_DNWELLING, FOURIER, &
              N, STERM_MASK_UP, STERM_MASK_DN, &
              NSTREAMS, NMOMENTS, N_USER_STREAMS, BEAM_PICUTOFF(CPT), &
              L0_WPARTIC(1,1,CPT), OMEGAMOMS_ELASTIC(1,0,CPT), &
              PLM_WT_PXI, PLM_WT_MXI,  PLM_00_PXUI, &
              PLM_PXUI,   PLM_00_MXUI, PLM_MXUI, &
              W_HELP, U_WPOS1, U_WPOS2, U_WNEG1, U_WNEG2 )

!  End layer loop

         ENDDO

!  Post processing Multipliers
!  ===========================

!  layer loop

         DO N = 1, NLAYERS

!  Homogeneous solution multipliers
!   @@@ RobFix 5/5/11. Small numbers analysis: added LRRS_FGSMALL,
!                      DELTAU_VERT_INPUT(N,CPT) to Argument list.

           CALL HOMOGMULT_1 &
                ( NSTREAMS, N, N_USER_STREAMS, LRRS_FGSMALL, &
                  SAVE_TRANS_USERM(1,1,CPT), USER_STREAMS,   &
                  L0_KEIGEN(1,1,CPT), L0_KTRANS(1,1,CPT),    &
                  DELTAU_VERT_INPUT(N,CPT),                  &
                  HMULT_1, HMULT_2, ZETA_P, ZETA_M )

!  Multipliers for the Upwelling solar source functions

           IF ( DO_UPWELLING ) THEN
             NOFLIP = .FALSE.
             CALL BEAMMULT_UP_1 &
               ( N, N_USER_STREAMS, USER_STREAMS, LRRS_FGSMALL, &
                 STERM_MASK_UP(N), LOCAL_TRANS_USERM(1,N), &
                 DELTAU_VERT_INPUT(N,CPT), &
                 NOFLIP, BEAM_PICUTOFF(CPT), BEAM_DTRANS(N,CPT), &
                 BEAM_AVSECANT(N,CPT), &
                 EMULT_UP(1,N) )
           ENDIF

!  Multipliers for the Downwelling solar source functions

           IF ( DO_DNWELLING ) THEN
             NOFLIP = .FALSE.
             CALL BEAMMULT_DN_1 &
               ( N, N_USER_STREAMS, USER_STREAMS, LRRS_FGSMALL, &
                 STERM_MASK_DN(N), LOCAL_TRANS_USERM(1,N), &
                 DELTAU_VERT_INPUT(N,CPT), &
                 NOFLIP, BEAM_PICUTOFF(CPT), BEAM_DTRANS(N,CPT), &
                 BEAM_AVSECANT(N,CPT), &
                 EMULT_DN(1,N) )
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

!  debug
!         write(*,*)'hello'
!           DO N = 1, NLAYERS
!             DO UM = 1, N_USER_STREAMS
!               write(87,'(2i3,1p3e15.6)')n,um,emult_up(um,n)
!           enddo
!          enddo
!         pause'hello'

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

         IF ( DO_UPWELLING ) THEN

           CALL ELASTIC_INTENSITY_UP_1 &
           ( DO_MSMODE_LRRS, DO_SSCORR_OVERALL, DO_ELASTIC_SCATTERING, &
             DO_QUAD_RESULTS, DO_INCLUDE_MVOUT, &
             DO_INCLUDE_SURFACE, DO_REFLECTED_DIRECTBEAM, FOURIER, &
             DO_LAMBERTIAN_SURFACE, R2, FL1(CPT), FL2(CPT), &
             USER_BIREFLEC, FLUX_MULTIPLIER, &
             NSTREAMS, NLAYERS, N_USER_STREAMS, N_LOUTPUT, LEVELMASK_UP, &
             QUAD_STRMWGT, USER_DIRECT_BEAM, WUPPER, WLOWER, &
             L0_KTRANS(1,1,CPT), SAVE_TRANS_USERM(1,1,CPT), &
             LCON, MCON, LCON_XVEC, MCON_XVEC, &
             U_XPOS, U_XNEG, U_WPOS1, U_WPOS2, &
             HMULT_1, HMULT_2, EMULT_UP, &
             IDOWNSURF, CUMSOURCE_UP, &
             ELASTIC_F_UP(1,1,APT), &
             QUADELASTIC_F_UP(1,1,APT) )

!         if ( fourier.eq.0)write(*,*)apt,lcon(4,23)
!      write(*,*)'PER',APT,FOURIER,ELASTIC_F_UP(1,1,APT), IDOWNSURF(2)
!           if ( fourier.eq.0)write(*,*)APT,USER_DIRECT_BEAM(1)

!           write(*,'(a,2i4,1pe17.8)')
!     *         'elastic up',fourier,apt, ELASTIC_F_UP(1,1,APT)
!           pause

!       write(*,*)'fou',APT,ELASTIC_F_UP(1,1,APT)
!         write(44,66)APT, U_XPOS(1,5,23),U_WPOS1(1,23),U_WPOS2(1,23)
! 66      format(i4,1p3e20.11)
!if ( APT.eq.8)stop

!         write(*,'(2i3,1p4e15.7)')fourier,apt,do_include_mvout,
!     *     ELASTIC_F_UP(1,1,APT), QUADELASTIC_F_UP(1,1,APT),
!     *     ELASTIC_F_UP(2,1,APT), QUADELASTIC_F_UP(2,1,APT)

!  Saved results for the Rayleigh planetary problem

!          IF ( DO_RAYLEIGH_ONLY ) THEN
!            DO J = 1, N_LOUTPUT
!             DO K = 1, N_USER_STREAMS
!               IF(M.EQ.0)ELASTIC_F0_UP(J,K,APT)=ELASTIC_F_UP(J,K,APT)
!               IF(M.EQ.1)ELASTIC_F1_UP(J,K,APT)=ELASTIC_F_UP(J,K,APT)
!               IF(M.EQ.2)ELASTIC_F2_UP(J,K,APT)=ELASTIC_F_UP(J,K,APT)
!             ENDDO
!            ENDDO
!           ENDIF

         ENDIF

!  Post processing of the Downwelling intensity

         IF ( DO_DNWELLING ) THEN
           CALL ELASTIC_INTENSITY_DN_1 &
           ( DO_MSMODE_LRRS, DO_SSCORR_OVERALL, DO_ELASTIC_SCATTERING, &
             DO_QUAD_RESULTS, DO_INCLUDE_MVOUT, &
             FOURIER, FLUX_MULTIPLIER, NSTREAMS, &
             N_USER_STREAMS, N_LOUTPUT, LEVELMASK_DN, &
             WLOWER, L0_KTRANS(1,1,CPT), &
             SAVE_TRANS_USERM(1,1,CPT), &
             LCON, MCON, LCON_XVEC, MCON_XVEC, &
             U_XPOS, U_XNEG, U_WNEG1, U_WNEG2, &
             HMULT_1, HMULT_2, EMULT_DN, &
             CUMSOURCE_DN, &
             ELASTIC_F_DN(1,1,APT), &
             QUADELASTIC_F_DN(1,1,APT) )

!             if (fourier.eq.2.and.apt.eq.1)pause'm1 apt 1 elasreg dn'

!  Saved results for the Rayleigh planetary problem

!           IF ( DO_RAYLEIGH_ONLY ) THEN
!             DO J = 1, N_LOUTPUT
!               DO K = 1, N_USER_STREAMS
!                IF(M.EQ.0)ELASTIC_F0_DN(J,K,APT)=ELASTIC_F_DN(J,K,APT)
!                IF(M.EQ.1)ELASTIC_F1_DN(J,K,APT)=ELASTIC_F_DN(J,K,APT)
!                IF(M.EQ.2)ELASTIC_F2_DN(J,K,APT)=ELASTIC_F_DN(J,K,APT)
!               ENDDO
!              ENDDO
!           ENDIF

!  End downwelling

         ENDIF

!  Miflux

         IF ( DO_INCLUDE_MVOUT ) THEN
           CALL MIFLUX_INTENSITY_1 &
               ( DO_UPWELLING, DO_DNWELLING, &
                 COS_SZA,  FLUX_FACTOR, &
                 NSTREAMS, N_LOUTPUT, LEVELMASK_DN, &
                 BEAM_PICUTOFF(CPT), BEAM_ETRANS(1,CPT), &
                 QUAD_WEIGHTS, QUAD_STRMWGT, &
                 QUADELASTIC_F_UP(1,1,APT), &
                 QUADELASTIC_F_DN(1,1,APT), &
                 MEAN_ELASTIC_UP(1,APT), FLUX_ELASTIC_UP(1,APT), &
                 MEAN_ELASTIC_DN(1,APT), FLUX_ELASTIC_DN(1,APT) )
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


      END MODULE lrrs_elastic_master_1

