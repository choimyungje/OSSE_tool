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
! #   SUBROUTINES :                                             #
! #                                                             #
! #         lrrs_raman_ops_bin_plus   (with linearization)      #
! #         lrrs_raman_ops_mono_plus  (with linearization)      #
! #                                                             #
! #  These subroutines calculate the linearizations of the      #
! #  Raman optical property inputs.                             #
! #                                                             #
! #  The calculation is based on the elastic IOP inputs which   #
! #  are assumed known, and the routines adjust these inputs    #
! #  and create new Raman-only IOP inputs by using the Raman    #
! #  cross-sections and the Air column density in each layer.   #
! #                                                             #
! #  Monochromatic Routines added 8 February 2011               #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Bookkeeping improvements (use of "Only", clearer I/O specifications)

!  -- Rob mod 5/12/17 for 2p5a, separate Bin/Mono RAMAN_OPS routines into 2 types: -
!       1. For unscaled quantities, Now including new phase function inputs
!       2. For scaled quantities,   Same as original subroutines.

      MODULE lrrs_l_generate_ramanops_m

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, LDU

      PRIVATE
      PUBLIC :: LRRS_RAMAN_OPS_BIN_PLUS_1,  LRRS_RAMAN_OPS_BIN_PLUS_2, &
                LRRS_RAMAN_OPS_MONO_PLUS_1, LRRS_RAMAN_OPS_MONO_PLUS_2

      CONTAINS

!  Linearized version

      SUBROUTINE LRRS_RAMAN_OPS_BIN_PLUS_1 &
        ( MAX_LAYERS, MAX_POINTS, MAX_GEOMETRIES, MAX_BINS, MAX_VARS,                & ! Input dimensions
          DO_UPWELLING, DO_DNWELLING, DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,        & ! Input flags
          DO_NORMALIZED_WFS, DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF, & ! Input flags
          NLAYERS, NGEOMETRIES, NPOINTS_INNER, OFFSET_INNER, NVARS,                  & ! Input integers
          COSSCAT_UP, COSSCAT_DN, TEMPERATURES, TEMPERATURES_UNSHIFTED,              & ! Inputs
          LAYER_AIRCOLUMNS, LAYER_AIRCOLUMNS_DT,                                     & ! Inputs
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, DELTAU_VERT, L_DELTAU_VERT,                 & ! Inputs Optical
          OMEGAPHASFUNC_ELASTIC_UP, L_OMEGAPHASFUNC_ELASTIC_UP,                      & ! Inputs Optical
          OMEGAPHASFUNC_ELASTIC_DN, L_OMEGAPHASFUNC_ELASTIC_DN,                      & ! Inputs Optical
          N_RRSBINS, BINXSEC, RRSXSEC_OUT, BINXSEC_DT, RRSXSEC_OUT_DT,               & ! Inputs Raman Spec
          OMEGAPHASFUNC_CABANNES_UP, L_OMEGAPHASFUNC_CABANNES_UP,                    & ! Outputs
          OMEGAPHASFUNC_CABANNES_DN, L_OMEGAPHASFUNC_CABANNES_DN,                    & ! Outputs
          OMEGAMOMS_RRSLOSS,  L_OMEGAMOMS_RRSLOSS,                                   & ! Output
          OMEGAMOMS_RRSBIN,   L_OMEGAMOMS_RRSBIN )                                     ! Outputs

       IMPLICIT NONE

!  Notes
!  =====

!  Flags
!  -----

!  There are separate flags for Air and Temperature PROFILE Jacobians.
!  There is a separate flag for Temperature-Shift   COLUMN  Jacobian.
!  There is also a flag for use of normalized derivatives (PROFILE only).

!  IOPS FOR THE LRRS MODEL, BINNING REALIZATION
!  ********************************************

!  Given a set of optical properties defined for an OUTER window with
!  total number of points NPOINTS_OUTER and wavelengths LAMBDAS_RANKED,
!  determine the necessary optical properties for the RRS model over
!  an INNER window which has total number of points NPOINTS_INNER
!  and starts after the offset value OFFSET_INNER. Each wavelength
!  is at the center of a bin defined by arrays BINLOWER and BINUPPER,
!  which are specified for the outer window.
!
!  Here is the indexing system for this subroutine:
!    n  is the layer index,            n  = 1, .... NLAYERS
!    v  is the Geometry index,         v  = 0, .... NGEOMETRIES
!    wo is the outer wavelength index, wo = 1, ...  NPOINTS_TOTAL
!    wi is the inner wavelength index, wi = 1, ...  NPOINTS_INNER
!    b  is the inner grid bin index,   b  = 1, ...  NBINS(wi)
!    q  is the Weight function index,  q  = 1, .... NVARS
!
!  The control inputs for the windows are then
!
!      NPOINTS_INNER  = number of points in inner window
!      OFFSET_INNER   = offset for start of inner window
!
!  The wavelengths and bin-limits are (these are not required here)
!
!      NPOINTS_OUTER      = number of points in outer window
!      LAMBDAS_RANKED(wo) = wavelengths (bin-centers) for the outer grid
!      BINLOWER(wo)       = Bin lower limit wavelengths for the outer gr
!      BINUPPER(wo)       = Bin upper limit wavelengths for the outer gr
!
!  The control integers for layers and geometries are
!
!      NLAYERS     = actual number of atmospheric layers
!      NGEOMETRIES = actual number of phase function geometries
!
!  There are scattering angle inputs to go with the phase function inputs
!
!      COSSCAT_UP (v) = Upwelling scattering angle cosines
!      COSSCAT_DN (v) = Dnwelling scattering angle cosines
!!
!  The already-computed optical properties are total quantities defined
!  for all elastic scatterers for the total number of layers NLAYERS.
!  The optical properties are :
!
!      DELTAU_VERT(n,wo)                = total optical depth for extinction.
!      OMEGAPHASFUNC_ELASTIC_UP(n,v,wo) = SSA x Elastic upwelling   phase function
!      OMEGAPHASFUNC_ELASTIC_DN(n,v,wo) = SSA x Elastic downwelling phase function
!      L_OMEGAPHASFUNC_ELASTIC_UP(q,n,v,wo) = Linearized (SSA x Elastic upwelling   phase function)
!      L_OMEGAPHASFUNC_ELASTIC_DN(q,n,v,wo) = Linearized (SSA x Elastic downwelling phase function)

!  Also required as input for this process are:
!
!      TEMPERATURES(n)           = Layer Temperature Profile, [K]
!      TEMPERATURES_UNSHIFTED(n) = Unshifted Temperature Profile, [K]

!      LAYER_AIRCOLUMNS(n)    = Air column-density Profile, [mol/cm^2]
!      LAYER_AIRCOLUMNS_dT(n) = Air Profile T-derivative,   [mol/cm^2/K]

!      RAYLEIGH_XSEC(wo)   = Rayleigh cross-sections on outer grid [cm^2/mol]
!      RAYLEIGH_DEPOL(wo)  = Rayleigh depolarization ratios on outer grid
!
!  The binning assignment is
!
!      NRRSBINS(wi) = number of bins for inner wavelength wi
!      BINMAP(b,wi) = Bin map        for inner wavelength   wi
!
!  The scattering outputs are:
!
!      OMEGAPHASFUNC_CABANNES_UP(n,v,wi) = SSA x Cabannes-adjusted upwelling   phase function
!      OMEGAPHASFUNC_CABANNES_DN(n,v,wi) = SSA x Cabannes-adjusted downwelling phase function
!      OMEGAMOMS_RRSBIN(n,l,b,wi) = SSA x RRS phase moments. L = 0, 1, 2
!                                   Raman Gain terms for each bin.
!      OMEGAMOMS_RRSLOSS(n,l,wi)  = SSA x RRS phase moments. L = 0, 1, 2.
!                                   Raman loss term for the inner point

!      L_OMEGAPHASFUNC_CABANNES_UP(n,v,wi) = Linearized SSA x Cabannes-adjusted upwelling   phasefunction
!      L_OMEGAPHASFUNC_CABANNES_DN(n,v,wi) = Linearized SSA x Cabannes-adjusted downwelling phasefunction
!      L_OMEGAMOMS_RRSBIN(q,n,l,b,wi) = Linearized Raman Gain terms for each bin
!      L_OMEGAMOMS_RRSLOSS(q,n,l,wi)  = Linearized Raman loss term for the inner point

!  Input arguments
!  ===============

!  Maximum number of spectral points and bins

      INTEGER  , INTENT(IN) :: MAX_POINTS, MAX_BINS

!  Maximum number of layers and geometries
!   -- Rob mod 5/12/17 for 2p5a, add MAX_GEOMETRIES, remove max_moments

      INTEGER  , INTENT(IN) :: MAX_LAYERS, MAX_GEOMETRIES

!  Maximum number of weighting functions

      INTEGER  , INTENT(IN) :: MAX_VARS

!  Flags
!  -----

!  Directional flags
!   -- Rob mod 5/12/17 for 2p5a, added

      LOGICAL  , INTENT(IN) :: DO_UPWELLING, DO_DNWELLING

!  RRS Method

      LOGICAL, INTENT(IN) :: DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Temperature and Air Profile Linearization flags

      LOGICAL, INTENT(IN) :: DO_TEMPPROFILE_WFS
      LOGICAL, INTENT(IN) :: DO_AIRPROFILE_WFS

!  Temperature shift Linearization flag (Column WF only)

      LOGICAL, INTENT(IN) :: DO_TEMPSHIFT_WF

!  use of normalized weighting functions (Profile WFs only)
!   T-shift Jacobian is automatically normalized

      LOGICAL, INTENT(IN) :: DO_NORMALIZED_WFS

!  Numbers
!  -------

!  Layering

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of Geometries
!   -- Rob mod 5/12/17 for 2p5a, add NGEOMETRIES, remove NMOMENTS

      INTEGER  , INTENT(IN) :: NGEOMETRIES

!  Number of points between the limiting wavelengths

      INTEGER  , INTENT(IN) :: NPOINTS_INNER

!  outer wavelength offset

      INTEGER  , INTENT(IN) :: OFFSET_INNER

!  Number of weighting functions

      INTEGER  , INTENT(IN) :: NVARS ( MAX_LAYERS )

!  Scattering angles
!   -- Rob mod 5/12/17 for 2p5a, added

      REAL(FPK), INTENT(IN) :: COSSCAT_UP ( MAX_GEOMETRIES )
      REAL(FPK), INTENT(IN) :: COSSCAT_DN ( MAX_GEOMETRIES )

!  Physical Properties
!  -------------------

!  Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  Layer temperatures (actual and unshifted)

      REAL(FPK), INTENT(IN) :: TEMPERATURES ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: TEMPERATURES_UNSHIFTED ( MAX_LAYERS )

!  Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )

!  Temperature derivatives of Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS_dT ( MAX_LAYERS )

!  Vertical optical depths

      REAL(FPK), INTENT(IN) :: DELTAU_VERT ( MAX_LAYERS, MAX_POINTS )

!  Input elastic scattering phase function products
!   -- Rob mod 5/12/17 for 2p5a, Replaces the phasmoms input (OMEGAMOMS_ELASTIC)

      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Linearization of optical properties
!   -- Rob mod 5/12/17 for 2p5a, Replaces the phasmoms input (OMEGAMOMS_ELASTIC)

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT ( MAX_VARS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_UP ( MAX_VARS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_DN ( MAX_VARS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Bin Mapping numbers

      INTEGER  , INTENT(IN) :: N_RRSBINS  ( MAX_POINTS )

!  Binned Gain term cross-sections

      REAL(FPK), INTENT(IN) :: BINXSEC ( MAX_LAYERS, MAX_POINTS, MAX_BINS )

!  Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT ( MAX_LAYERS, MAX_POINTS )

!  T-derivatives of Binned Gain term cross-sections

      REAL(FPK), INTENT(IN) :: BINXSEC_dT ( MAX_LAYERS,MAX_POINTS,MAX_BINS )

!  T-derivatives Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT_dT ( MAX_LAYERS, MAX_POINTS )

!  Output arguments
!  ================

!  Output Cabannes-adjusted scattering phase function products
!   -- Rob mod 5/12/17 for 2p5a, Replaces the phasmoms output (OMEGAMOMS_CABANNES)
!                                Only calculated for the Cabannes-Raman option

      REAL(FPK), INTENT(OUT) :: OMEGAPHASFUNC_CABANNES_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAPHASFUNC_CABANNES_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSBIN  ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Linearizations

      REAL(FPK), INTENT(OUT) :: L_OMEGAPHASFUNC_CABANNES_UP ( MAX_VARS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAPHASFUNC_CABANNES_DN ( MAX_VARS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSLOSS ( MAX_VARS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSBIN  ( MAX_VARS, MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Local variables
!  ===============

      INTEGER   :: WI, WO, N, B, Q, NQ, V
      REAL(FPK) :: DIF0(MAX_LAYERS), DIF2(MAX_LAYERS), ADJUST, RAMANPHASFUNC, NORM, NORM2
      REAL(FPK) :: L_DIF0(MAX_VARS,MAX_LAYERS), L_DIF2(MAX_VARS,MAX_LAYERS), L_ADJUST, RATIO_1
      REAL(FPK) :: SIGRAY, DEPOL, LOSSM0(MAX_LAYERS), LOSSM2(MAX_LAYERS)
      REAL(FPK) :: TSHIFT, LOSSM0_dT, LOSSM2_dT
      REAL(FPK) :: RAYMOM, CABMOM, RRS_PHASMOM2, AIRCOL, L_RAMANPHASFUNC
      REAL(FPK) :: TERM_1, TERM_2, TERM_3
      REAL(FPK) :: MU, LEG2_UP(MAX_GEOMETRIES), LEG2_DN(MAX_GEOMETRIES)

!  Initialize

      OMEGAMOMS_RRSBIN   = zero
      OMEGAMOMS_RRSLOSS  = zero
      OMEGAPHASFUNC_CABANNES_UP = zero
      OMEGAPHASFUNC_CABANNES_DN = zero
      L_OMEGAMOMS_RRSBIN   = zero
      L_OMEGAPHASFUNC_CABANNES_UP = zero
      L_OMEGAPHASFUNC_CABANNES_DN = zero

!  Legendre calculations
!   -- Rob mod 5/12/17 for 2p5a, need Legendre polynomials for Cabannes adjustment

      IF ( DO_CABANNES_RAMAN ) THEN
         DO V = 1, NGEOMETRIES
           MU = COSSCAT_UP(V) ; LEG2_UP(V) = 1.5_fpk * MU * MU - 0.5_fpk
           MU = COSSCAT_DN(V) ; LEG2_DN(V) = 1.5_fpk * MU * MU - 0.5_fpk
         ENDDO
      ENDIF

!  Optical properties on Inner wavelength loop
!  -------------------------------------------

!  Start the loop

      DO WI = 1, NPOINTS_INNER

!  Use offset to get point index in the outer loop

        WO = WI + OFFSET_INNER

!  Cabannes adjustment
!  -------------------

        IF ( DO_CABANNES_RAMAN ) THEN

!  Make Cabannes cross-sections------------------------
!    Essentially adjusting the Rayleigh scattering part to make sure tha
!    the elastic component is just the due to Cabannes scattering
!    Only requires the first 0,1,2 moments to be adjusted.

          DEPOL  = RAYLEIGH_DEPOL(WO)
          SIGRAY = RAYLEIGH_XSEC(WO)

!  DIF0 = RRS Loss-scattering optical depth / total optical depth
!  DIF2 = scatter-weighted RRS Loss / total OD, for moment L = 2

!      Rayleigh and Cabannes Legendre(2) moments

          CABMOM = 0.25D0 * &
             ( 8.0D0 - 9.0D0*DEPOL ) / ( 4.0D0 - 3.0D0*DEPOL )
          RAYMOM = ( one - DEPOL ) / ( 2.0D0 + DEPOL )

!   -- Rob mod 5/12/17 for 2p5a, Adjust the phase functions for Cabannes scattering

!             DIF0 = RRS Loss-scattering optical depth / total optical depth
!             DIF2 = scatter-weighted RRS Loss / total OD, for moment L = 2

          DO N = 1, NLAYERS
            AIRCOL = LAYER_AIRCOLUMNS(N)
            LOSSM0(N) = RRSXSEC_OUT(N,WI)
            LOSSM2(N) = SIGRAY*RAYMOM - (SIGRAY - LOSSM0(N)) * CABMOM
            ADJUST = AIRCOL / DELTAU_VERT(N,WO)
            DIF0(N)   = LOSSM0(N) * ADJUST
            DIF2(N)   = LOSSM2(N) * ADJUST
          ENDDO

!  Now make the adjustments

          IF ( DO_UPWELLING ) THEN
            DO V = 1, NGEOMETRIES
              DO N = 1, NLAYERS
                RAMANPHASFUNC = DIF0(N) + DIF2(N) * LEG2_UP(V)
                OMEGAPHASFUNC_CABANNES_UP(N,V,WI) = OMEGAPHASFUNC_ELASTIC_UP(N,V,WO) - RAMANPHASFUNC
              ENDDO
            ENDDO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            DO V = 1, NGEOMETRIES
              DO N = 1, NLAYERS
                RAMANPHASFUNC = DIF0(N) + DIF2(N) * LEG2_DN(V)
                OMEGAPHASFUNC_CABANNES_DN(N,V,WI) = OMEGAPHASFUNC_ELASTIC_DN(N,V,WO) - RAMANPHASFUNC
              ENDDO
            ENDDO
          ENDIF

!  Here is the original adjustment of the moments
!     Cabannes L = 0, adjust OMEGAMOMS by DIF0
!     Cabannes L = 2, adjust OMEGAMOMS by DIF2
!     Cabannes L = 1,3,4.... no adjustment (just equal to elastic values)
!          DO N = 1, NLAYERS
!           OMEGAMOMS_CABANNES(N,0,WI) = OMEGAMOMS_ELASTIC(N,0,WO) - DIF0(N)
!           OMEGAMOMS_CABANNES(N,1,WI) = OMEGAMOMS_ELASTIC(N,1,WO)
!           OMEGAMOMS_CABANNES(N,2,WI) = OMEGAMOMS_ELASTIC(N,2,WO) - DIF2(N)
!           OMEGAMOMS_CABANNES (N,3:NMOMENTS,WI) = OMEGAMOMS_ELASTIC(N,3:NMOMENTS,WO)
!          ENDDO

!  Linearized Adjustment
!  ---------------------

          DO N = 1, NLAYERS

!  Counting: Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

            NQ = NVARS(N)
            IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
            IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
            IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type). Derivatives can be normalized or un-normalized.

            DO Q = 1, NQ
              RATIO_1 = - L_DELTAU_VERT(Q,N,WO) / DELTAU_VERT(N,WO)
              L_DIF0(Q,N)   = DIF0(N) * RATIO_1
              L_DIF2(Q,N)   = DIF2(N) * RATIO_1
            ENDDO
            Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG

            IF ( DO_AIRPROFILE_WFS ) THEN
              Q = Q + 1
              NORM = one ;  IF ( DO_NORMALIZED_WFS) NORM = LAYER_AIRCOLUMNS(N)
              TERM_1 = L_DELTAU_VERT(Q,N,WO)  * DIF0(N)
              L_DIF0(Q,N) = (LOSSM0(N) * NORM - TERM_1 ) / DELTAU_VERT(N,WO)
              TERM_1 = L_DELTAU_VERT(Q,N,WO)  * DIF2(N)
              L_DIF2(Q,N) = (LOSSM2(N) * NORM - TERM_1 ) / DELTAU_VERT(N,WO)
            ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

            IF ( DO_TEMPPROFILE_WFS ) THEN
              Q = Q + 1
              NORM = RRSXSEC_OUT_dT(N,WI) ; IF ( DO_NORMALIZED_WFS ) NORM = NORM * TEMPERATURES(N)
              NORM2 = NORM * CABMOM
              TERM_1 = NORM * LAYER_AIRCOLUMNS(N)
              TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM0(N)
              TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF0(N)
              IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
              L_DIF0(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
              TERM_1 = NORM2 * LAYER_AIRCOLUMNS(N)
              TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM2(N)
              TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF2(N)
              IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
              L_DIF2(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
            ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here.

            IF ( DO_TEMPSHIFT_WF ) THEN
              Q = Q + 1
              TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
              LOSSM0_dT = RRSXSEC_OUT_dT(N,WI)
              LOSSM2_dT = LOSSM0_dT * CABMOM
              TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM0_dT
              TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM0(N)
              TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF0(N)
              L_DIF0(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
              TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM2_dT
              TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM2(N)
              TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF2(N)
              L_DIF2(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
            ENDIF

!  End saved layer loop

          ENDDO

!  Now make the adjustments

          IF ( DO_UPWELLING ) THEN
            DO V = 1, NGEOMETRIES
              DO N = 1, NLAYERS
                DO Q = 1, NVARS(N)
                  L_RAMANPHASFUNC = L_DIF0(Q,N) + L_DIF2(Q,N) * LEG2_UP(V)
                  L_OMEGAPHASFUNC_CABANNES_UP(Q,N,V,WI) = L_OMEGAPHASFUNC_ELASTIC_UP(Q,N,V,WO) - L_RAMANPHASFUNC
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            DO V = 1, NGEOMETRIES
              DO N = 1, NLAYERS
                DO Q = 1, NVARS(N)
                  L_RAMANPHASFUNC = L_DIF0(Q,N) + L_DIF2(Q,N) * LEG2_DN(V)
                  L_OMEGAPHASFUNC_CABANNES_DN(Q,N,V,WI) = L_OMEGAPHASFUNC_ELASTIC_DN(Q,N,V,WO) - L_RAMANPHASFUNC
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  Here is the original adjustment of the moments
!     Cabannes L = 0, adjust OMEGAMOMS by DIF0
!     Cabannes L = 2, adjust OMEGAMOMS by DIF2
!     Cabannes L = 1,3,4.... no adjustment (just equal to elastic values)
!          DO N = 1, NLAYERS
!             DO Q = 1, NVARS(N)
!               L_OMEGAMOMS_CABANNES(Q,N,0,WI) = L_OMEGAMOMS_ELASTIC(Q,N,0,WO) - L_DIF0
!               L_OMEGAMOMS_CABANNES(Q,N,1,WI) = L_OMEGAMOMS_ELASTIC(Q,N,1,WO)
!               L_OMEGAMOMS_CABANNES(Q,N,2,WI) = L_OMEGAMOMS_ELASTIC(Q,N,2,WO) - L_DIF2
!               L_OMEGAMOMS_CABANNES(Q,N,3:NMOMENTS,WI) = L_OMEGAMOMS_ELASTIC(Q,N,3:NMOMENTS,WO)
!             ENDDO
!           ENDDO

!  End Cabannes-Raman treatment

        ENDIF

!  Loss term (Energy balancing mode only)
!  --------------------------------------

        IF ( DO_ENERGY_BALANCING ) THEN

!  Raman phase function moment

          RRS_PHASMOM2 = 0.05D0

!  Start layer loop

          DO N = 1, NLAYERS

!  Loss term

           AIRCOL = LAYER_AIRCOLUMNS(N)
           ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,WO)
           LOSSM0(N) = RRSXSEC_OUT(N,WI)
           DIF0(N)   = LOSSM0(N) * ADJUST
           OMEGAMOMS_RRSLOSS(N,0,WI) = DIF0(N)
           OMEGAMOMS_RRSLOSS(N,1,WI) = ZERO
           OMEGAMOMS_RRSLOSS(N,2,WI) = DIF0(N) * RRS_PHASMOM2

!  Counting
!    Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

           NQ = NVARS(N)
           IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
           IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
           IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type)
!    Derivatives can be normalized or un-normalized.

           DO Q = 1, NQ
             RATIO_1 = - L_DELTAU_VERT(Q,N,WO) / DELTAU_VERT(N,WO)
             L_DIF0(Q,N)  = DIF0(N) * RATIO_1
             L_OMEGAMOMS_RRSLOSS(Q,N,0,WI) = L_DIF0(Q,N)
             L_OMEGAMOMS_RRSLOSS(Q,N,2,WI) = L_DIF0(Q,N) * RRS_PHASMOM2
           ENDDO
           Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_AIRPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = ONE
             IF ( DO_NORMALIZED_WFS) NORM = AIRCOL
             TERM_1 = L_DELTAU_VERT(Q,N,WO)  * DIF0(N)
             L_DIF0(Q,N) = (LOSSM0(N) * NORM - TERM_1 ) / DELTAU_VERT(N,WO)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,WI) = L_DIF0(Q,N)
             L_OMEGAMOMS_RRSLOSS(Q,N,2,WI) = L_DIF0(Q,N) * RRS_PHASMOM2
           ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_TEMPPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = RRSXSEC_OUT_dT(N,WI)
             TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF0(N)
             TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM0(N)
             TERM_1 = NORM * LAYER_AIRCOLUMNS(N)
             IF ( DO_NORMALIZED_WFS ) TERM_1 = TERM_1 * TEMPERATURES(N)
             IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
             L_DIF0(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,WI) = L_DIF0(Q,N)
             L_OMEGAMOMS_RRSLOSS(Q,N,2,WI) = L_DIF0(Q,N) * RRS_PHASMOM2
           ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here.

           IF ( DO_TEMPSHIFT_WF ) THEN
             Q = Q + 1
             TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
             LOSSM0_dT = RRSXSEC_OUT_dT(N,WI)
             TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM0_dT
             TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM0(N)
             TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF0(N)
             L_DIF0(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,WI) = L_DIF0(Q,N)
             L_OMEGAMOMS_RRSLOSS(Q,N,2,WI) = L_DIF0(Q,N) * RRS_PHASMOM2
           ENDIF

!  First moment linearization is zero, always

           DO Q = 1, NVARS(N)
              L_OMEGAMOMS_RRSLOSS(Q,N,1,WI) = ZERO
           ENDDO

!  End layer loop

          ENDDO

!  End Energy-balancing loss-term clause

        ENDIF

!  Gain terms (both implementations)
!  ---------------------------------

!  Inelastic scattering: RRS input---------------
!  Calculate optical property input for the Gain RRS source terms
!  The second moment for RRS scattering = 1/20.
!  The adjustment = aircolumn / total optical thickness

        RRS_PHASMOM2 = 0.05D0

!  Start layer loop

        DO N = 1, NLAYERS

!  Bin optical properties

          AIRCOL = LAYER_AIRCOLUMNS(N)
          ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,WO)
          DO B = 1, N_RRSBINS(WI)
           OMEGAMOMS_RRSBIN(N,0,B,WI) = ADJUST * BINXSEC(N,WI,B)
           OMEGAMOMS_RRSBIN(N,1,B,WI) = ZERO
           OMEGAMOMS_RRSBIN(N,2,B,WI) = &
                     OMEGAMOMS_RRSBIN(N,0,B,WI) * RRS_PHASMOM2
          ENDDO

!  Counting
!    Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

          NQ = NVARS(N)
          IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
          IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
          IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type)
!    Derivatives can be normalized or un-normalized.

          DO Q = 1, NQ
            RATIO_1  = - L_DELTAU_VERT(Q,N,WO) / DELTAU_VERT(N,WO)
            L_ADJUST = ADJUST * RATIO_1
            DO B = 1, N_RRSBINS(WI)
             L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) = L_ADJUST*BINXSEC(N,WI,B)
             L_OMEGAMOMS_RRSBIN(Q,N,2,B,WI) = &
                   L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) * RRS_PHASMOM2
            ENDDO
          ENDDO
          Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_AIRPROFILE_WFS ) THEN
            Q = Q + 1
            NORM = ONE
            IF ( DO_NORMALIZED_WFS) NORM = AIRCOL
            TERM_1 = L_DELTAU_VERT(Q,N,WO) * ADJUST
            L_ADJUST = ( NORM - TERM_1 ) / DELTAU_VERT(N,WO)
            DO B = 1, N_RRSBINS(WI)
             L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) = L_ADJUST*BINXSEC(N,WI,B)
             L_OMEGAMOMS_RRSBIN(Q,N,2,B,WI) = &
                   L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) * RRS_PHASMOM2
            ENDDO
          ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

          IF ( DO_TEMPPROFILE_WFS ) THEN
            Q = Q + 1
            TERM_3   = L_DELTAU_VERT(Q,N,WO) * ADJUST
            TERM_2   = LAYER_AIRCOLUMNS_dT(N)
            IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
            L_ADJUST = ( TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
            NORM = ADJUST
            IF ( DO_NORMALIZED_WFS ) NORM = NORM * TEMPERATURES(N)
            DO B = 1, N_RRSBINS(WI)
             L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) = &
                L_ADJUST*BINXSEC(N,WI,B) + NORM * BINXSEC_dT(N,WI,B)
             L_OMEGAMOMS_RRSBIN(Q,N,2,B,WI) = &
                L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) * RRS_PHASMOM2
            ENDDO
          ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here.

          IF ( DO_TEMPSHIFT_WF ) THEN
            Q = Q + 1
            TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
            TERM_1 = TSHIFT * ADJUST
            TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N)
            TERM_3 = L_DELTAU_VERT(Q,N,WO) * ADJUST
            L_ADJUST = ( TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
            DO B = 1, N_RRSBINS(WI)
             L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) = &
                L_ADJUST*BINXSEC(N,WI,B) + TERM_1 * BINXSEC_dT(N,WI,B)
             L_OMEGAMOMS_RRSBIN(Q,N,2,B,WI) = &
                L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) * RRS_PHASMOM2
            ENDDO
          ENDIF

!  All first moment entries are zero

          DO Q = 1, NVARS(N)
            DO B = 1, N_RRSBINS(WI)
              L_OMEGAMOMS_RRSBIN(Q,N,1,B,WI) = ZERO
            ENDDO
          ENDDO

!  End layer loop

        ENDDO

!  End loop over all inner wavelengths

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_BIN_PLUS_1

!

      SUBROUTINE LRRS_RAMAN_OPS_BIN_PLUS_2 &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, MAX_BINS, MAX_VARS,   & ! Inputs
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, DO_NORMALIZED_WFS, & ! Inputs
          DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF,    & ! Inputs
          NLAYERS, NMOMENTS, NPOINTS_INNER, OFFSET_INNER, NVARS,     & ! Inputs
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL,                             & ! Inputs
          TEMPERATURES, TEMPERATURES_UNSHIFTED,                      & ! Inputs
          LAYER_AIRCOLUMNS, LAYER_AIRCOLUMNS_DT,                     & ! Inputs
          DELTAU_VERT,     OMEGAMOMS_ELASTIC,                        & ! Inputs
          L_DELTAU_VERT, L_OMEGAMOMS_ELASTIC,                        & ! Inputs
          N_RRSBINS, BINXSEC, RRSXSEC_OUT,                           & ! Inputs
          BINXSEC_DT, RRSXSEC_OUT_DT,                                & ! Inputs
          OMEGAMOMS_CABANNES, L_OMEGAMOMS_CABANNES,  & ! Outputs
          OMEGAMOMS_RRSLOSS,  L_OMEGAMOMS_RRSLOSS,   & ! Outputs
          OMEGAMOMS_RRSBIN,   L_OMEGAMOMS_RRSBIN )     ! Outputs

       IMPLICIT NONE

!  Purpose
!  =======

!  This is a generic routine for generating ROPS and their linearizations.

!  Notes
!  =====

!  Flags
!  -----

!  There are separate flags for Air and Temperature PROFILE Jacobians.
!  There is a separate flag for Temperature-Shift   COLUMN  Jacobian.
!  There is also a flag for use of normalized derivatives (PROFILE only).

!  IOPS FOR THE LRRS MODEL, BINNING REALIZATION
!  ********************************************

!  Given a set of optical properties defined for an OUTER window with
!  total number of points NPOINTS_OUTER and wavelengths LAMBDAS_RANKED,
!  determine the necessary optical properties for the RRS model over
!  an INNER window which has total number of points NPOINTS_INNER
!  and starts after the offset value OFFSET_INNER. Each wavelength
!  is at the center of a bin defined by arrays BINLOWER and BINUPPER,
!  which are specified for the outer window.
!
!  Here is the indexing system for this subroutine:
!    n  is the layer index,            n  = 1, .... NLAYERS
!    wo is the outer wavelength index, wo = 1, ...  NPOINTS_TOTAL
!    wi is the inner wavelength index, wi = 1, ...  NPOINTS_INNER
!    b  is the inner grid bin index,   b  = 1, ...  NBINS(wi)
!    l  is the Legendre moment index,  l  = 0, .... NMOMENTS_TOTAL
!    q  is the Weight function index,  q  = 1, .... NVARS
!
!  The control inputs for the windows are then
!
!      NPOINTS_INNER  = number of points in inner window
!      OFFSET_INNER   = offset for start of inner window
!
!  The wavelengths and bin-limits are (these are not required here)
!
!      NPOINTS_OUTER      = number of points in outer window
!      LAMBDAS_RANKED(wo) = wavelengths (bin-centers) for the outer grid
!      BINLOWER(wo)       = Bin lower limit wavelengths for the outer gr
!      BINUPPER(wo)       = Bin upper limit wavelengths for the outer gr
!
!  The control INTEGERs for layers and moments are
!
!      NLAYERS  = actual number of atmospheric layers
!      NMOMENTS = actual number of phase function moments
!
!  The already-computed optical properties are total quantities defined
!  for all elastic scatterers for the total number of layers NLAYERS.
!  The optical properties are :
!
!      DELTAU_VERT(n,wo)         = total optical depth for extinction.
!      OMEGAMOMS_ELASTIC(n,l,wo) = SSA x Elastic phasemoments
!                                   (l = 0, 1, ...NMOMENTS_TOTAL)
!    n  is the layer index,            n  = 1, .... NLAYERS
!    wo is the outer wavelength index, wo = 1, ...  NPOINTS_TOTAL
!    wi is the inner wavelength index, wi = 1, ...  NPOINTS_INNER
!    b  is the inner grid bin index,   b  = 1, ...  NBINS(wi)
!    l  is the Legendre moment index,  l  = 0, .... NMOMENTS_TOTAL
!
!  Also required as input for this process are:
!
!      TEMPERATURES(n)           = Layer Temperature Profile, [K]
!      TEMPERATURES_UNSHIFTED(n) = Unshifted Temperature Profile, [K]

!      LAYER_AIRCOLUMNS(n)    = Air column-density Profile, [mol/cm^2]
!      LAYER_AIRCOLUMNS_dT(n) = Air Profile T-derivative,   [mol/cm^2/K]

!      RAYLEIGH_XSEC(wo)   = Rayleigh cross-sections on outer grid [cm^2/mol]
!      RAYLEIGH_DEPOL(wo)  = Rayleigh depolarization ratios on outer grid
!
!  The binning assignment is
!
!      NRRSBINS(wi) = number of bins for inner wavelength wi
!      BINMAP(b,wi) = Bin map        for inner wavelength   wi
!
!  The scattering outputs are:
!
!      OMEGAMOMS_CABANNES(n,l,wi) = SSA x Cabannes-adjusted phasemoments
!                                   (l = 0, 1, ...NMOMENTS_TOTAL)
!      OMEGAMOMS_RRSBIN(n,l,b,wi) = SSA x RRS phasemoments. L = 0, 1, 2
!                                   Raman Gain terms for each bin.
!      OMEGAMOMS_RRSLOSS(n,l,wi)  = SSA x RRS phasemoments. L = 0, 1, 2.
!                                   Raman loss term for the iner point

!  Input arguments
!  ===============

!  Maximum number of spectral points and bins

      INTEGER  , INTENT(IN) :: MAX_POINTS, MAX_BINS

!  Maximum number of layers and input moments

      INTEGER  , INTENT(IN) :: MAX_LAYERS, MAX_MOMENTS

!  Maximum number of weighting functions

      INTEGER  , INTENT(IN) :: MAX_VARS

!  Flags
!  -----

!  RRS Method

      LOGICAL, INTENT(IN) :: DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Temperature and Air Profile Linearization flags

      LOGICAL, INTENT(IN) :: DO_TEMPPROFILE_WFS
      LOGICAL, INTENT(IN) :: DO_AIRPROFILE_WFS

!  Temperature shift Linearization flag (Column WF only)

      LOGICAL, INTENT(IN) :: DO_TEMPSHIFT_WF

!  use of normalized weighting functions (Profile WFs only)
!   T-shift Jacobian is automatically normalized

      LOGICAL, INTENT(IN) :: DO_NORMALIZED_WFS

!  Numbers
!  -------

!  Layering

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of elastic moments

      INTEGER  , INTENT(IN) :: NMOMENTS

!  Number of points between the limiting wavelengths

      INTEGER  , INTENT(IN) :: NPOINTS_INNER

!  outer wavelength offset

      INTEGER  , INTENT(IN) :: OFFSET_INNER

!  Number of weighting functions

      INTEGER  , INTENT(IN) :: NVARS ( MAX_LAYERS )

!  Physical Properties
!  -------------------

!  Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  Layer temperatures (actual and unshifted)

      REAL(FPK), INTENT(IN) :: TEMPERATURES ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: TEMPERATURES_UNSHIFTED ( MAX_LAYERS )

!  Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )

!  Temperature derivatives of Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS_dT ( MAX_LAYERS )

!  Vertical optical depths

      REAL(FPK), INTENT(IN) :: DELTAU_VERT ( MAX_LAYERS, MAX_POINTS )

!  Input elastic scattering law

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Linearization of optical properties

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT &
                ( MAX_VARS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_ELASTIC &
           ( MAX_VARS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Bin Mapping numbers

      INTEGER  , INTENT(IN) :: N_RRSBINS  ( MAX_POINTS )

!  Binned Gain term cross-sections

      REAL(FPK), INTENT(IN) :: BINXSEC ( MAX_LAYERS, MAX_POINTS, MAX_BINS )

!  Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT ( MAX_LAYERS, MAX_POINTS )

!  T-derivatives of Binned Gain term cross-sections

      REAL(FPK), INTENT(IN) :: BINXSEC_dT ( MAX_LAYERS,MAX_POINTS,MAX_BINS )

!  T-derivatives Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT_dT ( MAX_LAYERS, MAX_POINTS )

!  Output arguments
!  ================

!  Scattering quantities

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSBIN ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Linearizations

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_CABANNES ( MAX_VARS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSLOSS ( MAX_VARS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSBIN ( MAX_VARS, MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Local variables
!  ===============

      INTEGER   :: WI, WO, N, L, B, Q, NQ
      REAL(FPK) :: DIF0, DIF2, ADJUST, NORM, NORM2
      REAL(FPK) :: L_DIF0, L_DIF2, L_ADJUST, RATIO_1
      REAL(FPK) :: SIGRAY, DEPOL, LOSSM0, LOSSM2
      REAL(FPK) :: TSHIFT, LOSSM0_dT, LOSSM2_dT
      REAL(FPK) :: RAYMOM, CABMOM, RRS_PHASMOM2, AIRCOL
      REAL(FPK) :: TERM_1, TERM_2, TERM_3

!  Initialize

      OMEGAMOMS_RRSBIN   = zero
      OMEGAMOMS_RRSLOSS  = zero
      OMEGAMOMS_CABANNES = zero
      L_OMEGAMOMS_RRSBIN   = zero
      L_OMEGAMOMS_RRSLOSS  = zero
      L_OMEGAMOMS_CABANNES = zero

!  Optical properties on Inner wavelength loop
!  -------------------------------------------

!  Start the loop

      DO WI = 1, NPOINTS_INNER

!  Use offset to get point index in the outer loop

        WO = WI + OFFSET_INNER

!  Cabannes adjustment
!  -------------------

        IF ( DO_CABANNES_RAMAN ) THEN

!  Make Cabannes cross-sections------------------------
!    Essentially adjusting the Rayleigh scattering part to make sure tha
!    the elastic component is just the due to Cabannes scattering
!    Only requires the first 0,1,2 moments to be adjusted.

          DEPOL  = RAYLEIGH_DEPOL(WO)
          SIGRAY = RAYLEIGH_XSEC(WO)

!  DIF0 = RRS Loss-scattering optical depth / total optical depth
!  DIF2 = scatter-weighted RRS Loss / total OD, for moment L = 2

!      Rayleigh and Cabannes Legendre(2) moments

          CABMOM = 0.25D0 * &
             ( 8.0D0 - 9.0D0*DEPOL ) / ( 4.0D0 - 3.0D0*DEPOL )
          RAYMOM = ( one - DEPOL ) / ( 2.0D0 + DEPOL )

!  Cabannes L = 0, adjust OMEGAMOMS by DIF0
!  Cabannes L = 2, adjust OMEGAMOMS by DIF2
!  Cabannes L = 1,3,4.. no adjustment (just equal to elastic values)

          DO N = 1, NLAYERS

!  Basic quantities

           AIRCOL = LAYER_AIRCOLUMNS(N)
           LOSSM0 = RRSXSEC_OUT(N,WI)
           LOSSM2 = SIGRAY*RAYMOM - (SIGRAY - LOSSM0) * CABMOM
           ADJUST = AIRCOL / DELTAU_VERT(N,WO)
           DIF0   = LOSSM0 * ADJUST
           DIF2   = LOSSM2 * ADJUST
           OMEGAMOMS_CABANNES(N,0,WI) = OMEGAMOMS_ELASTIC(N,0,WO) - DIF0
           OMEGAMOMS_CABANNES(N,1,WI) = OMEGAMOMS_ELASTIC(N,1,WO)
           OMEGAMOMS_CABANNES(N,2,WI) = OMEGAMOMS_ELASTIC(N,2,WO) - DIF2
           DO L = 3, NMOMENTS
            OMEGAMOMS_CABANNES (N,L,WI) = OMEGAMOMS_ELASTIC(N,L,WO)
           ENDDO

!  Counting
!    Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

           NQ = NVARS(N)
           IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
           IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
           IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type)
!    Derivatives can be normalized or un-normalized.

           DO Q = 1, NQ
             RATIO_1 = - L_DELTAU_VERT(Q,N,WO) / DELTAU_VERT(N,WO)
             L_DIF0   = DIF0 * RATIO_1
             L_DIF2   = DIF2 * RATIO_1
             L_OMEGAMOMS_CABANNES(Q,N,0,WI) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,0,WO) - L_DIF0
             L_OMEGAMOMS_CABANNES(Q,N,2,WI) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,2,WO) - L_DIF2
           ENDDO
           Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_AIRPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = one
             IF ( DO_NORMALIZED_WFS) NORM = LAYER_AIRCOLUMNS(N)
             TERM_1 = L_DELTAU_VERT(Q,N,WO)  * DIF0
             L_DIF0 = (LOSSM0 * NORM - TERM_1 ) / DELTAU_VERT(N,WO)
             TERM_1 = L_DELTAU_VERT(Q,N,WO)  * DIF2
             L_DIF2 = (LOSSM2 * NORM - TERM_1 ) / DELTAU_VERT(N,WO)
             L_OMEGAMOMS_CABANNES(Q,N,0,WI) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,0,WO) - L_DIF0
             L_OMEGAMOMS_CABANNES(Q,N,2,WI) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,2,WO) - L_DIF2
           ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_TEMPPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = RRSXSEC_OUT_dT(N,WI)
             IF ( DO_NORMALIZED_WFS ) NORM = NORM * TEMPERATURES(N)
             NORM2 = NORM * CABMOM
             TERM_1 = NORM * LAYER_AIRCOLUMNS(N)
             TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM0 
             TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF0
             IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
             L_DIF0 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
             TERM_1 = NORM2 * LAYER_AIRCOLUMNS(N)
             TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM2
             TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF2
             IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
             L_DIF2 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
             L_OMEGAMOMS_CABANNES(Q,N,0,WI) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,0,WO) - L_DIF0
             L_OMEGAMOMS_CABANNES(Q,N,2,WI) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,2,WO) - L_DIF2
           ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here.

           IF ( DO_TEMPSHIFT_WF ) THEN
             Q = Q + 1
             TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
             LOSSM0_dT = RRSXSEC_OUT_dT(N,WI)
             LOSSM2_dT = LOSSM0_dT * CABMOM
             TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM0_dT
             TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM0
             TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF0
             L_DIF0 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
             TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM2_dT
             TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM2
             TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF2
             L_DIF2 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
             L_OMEGAMOMS_CABANNES(Q,N,0,WI) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,0,WO) - L_DIF0
             L_OMEGAMOMS_CABANNES(Q,N,2,WI) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,2,WO) - L_DIF2
           ENDIF

!  1st moment, Elastic terms are all just copies

           DO Q = 1, NVARS(N)
             L_OMEGAMOMS_CABANNES(Q,N,1,WI) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,1,WO)
             DO L = 3, NMOMENTS
              L_OMEGAMOMS_CABANNES(Q,N,L,WI) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,L,WO)
             ENDDO
           ENDDO

!  End layer loop

          ENDDO

!  End Cabannes-Raman treatment

        ENDIF

!  Loss term (Energy balancing mode only)
!  --------------------------------------

        IF ( DO_ENERGY_BALANCING ) THEN

!  Raman phase function moment

          RRS_PHASMOM2 = 0.05D0

!  Start layer loop

          DO N = 1, NLAYERS

!  Loss term

           AIRCOL = LAYER_AIRCOLUMNS(N)
           ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,WO)
           LOSSM0 = RRSXSEC_OUT(N,WI)
           DIF0   = LOSSM0 * ADJUST
           OMEGAMOMS_RRSLOSS(N,0,WI) = DIF0
           OMEGAMOMS_RRSLOSS(N,1,WI) = ZERO
           OMEGAMOMS_RRSLOSS(N,2,WI) = DIF0 * RRS_PHASMOM2

!  Counting
!    Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

           NQ = NVARS(N)
           IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
           IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
           IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type)
!    Derivatives can be normalized or un-normalized.

           DO Q = 1, NQ
             RATIO_1 = - L_DELTAU_VERT(Q,N,WO) / DELTAU_VERT(N,WO)
             L_DIF0   = DIF0 * RATIO_1
             L_OMEGAMOMS_RRSLOSS(Q,N,0,WI) = L_DIF0
             L_OMEGAMOMS_RRSLOSS(Q,N,2,WI) = L_DIF0 * RRS_PHASMOM2
           ENDDO
           Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_AIRPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = ONE
             IF ( DO_NORMALIZED_WFS) NORM = AIRCOL
             TERM_1 = L_DELTAU_VERT(Q,N,WO)  * DIF0
             L_DIF0 = (LOSSM0 * NORM - TERM_1 ) / DELTAU_VERT(N,WO)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,WI) = L_DIF0
             L_OMEGAMOMS_RRSLOSS(Q,N,2,WI) = L_DIF0 * RRS_PHASMOM2
           ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_TEMPPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = RRSXSEC_OUT_dT(N,WI)
             TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF0
             TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM0
             TERM_1 = NORM * LAYER_AIRCOLUMNS(N)
             IF ( DO_NORMALIZED_WFS ) TERM_1 = TERM_1 * TEMPERATURES(N)
             IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
             L_DIF0 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,WI) = L_DIF0
             L_OMEGAMOMS_RRSLOSS(Q,N,2,WI) = L_DIF0 * RRS_PHASMOM2
           ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here.

           IF ( DO_TEMPSHIFT_WF ) THEN
             Q = Q + 1
             TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
             LOSSM0_dT = RRSXSEC_OUT_dT(N,WI)
             TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM0_dT
             TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM0
             TERM_3 = L_DELTAU_VERT(Q,N,WO)  * DIF0
             L_DIF0 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,WI) = L_DIF0
             L_OMEGAMOMS_RRSLOSS(Q,N,2,WI) = L_DIF0 * RRS_PHASMOM2
           ENDIF

!  First moment linearization is zero, always

           DO Q = 1, NVARS(N)
              L_OMEGAMOMS_RRSLOSS(Q,N,1,WI) = ZERO
           ENDDO

!  End layer loop

          ENDDO

!  End Energy-balancing loss-term clause

        ENDIF

!  Gain terms (both implementations)
!  ---------------------------------

!  Inelastic scattering: RRS input---------------
!  Calculate optical property input for the Gain RRS source terms
!  The second moment for RRS scattering = 1/20.
!  The adjustment = aircolumn / total optical thickness

        RRS_PHASMOM2 = 0.05D0

!  Start layer loop

        DO N = 1, NLAYERS

!  Bin optical properties

          AIRCOL = LAYER_AIRCOLUMNS(N)
          ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,WO)
          DO B = 1, N_RRSBINS(WI)
           OMEGAMOMS_RRSBIN(N,0,B,WI) = ADJUST * BINXSEC(N,WI,B)
           OMEGAMOMS_RRSBIN(N,1,B,WI) = ZERO
           OMEGAMOMS_RRSBIN(N,2,B,WI) = &
                     OMEGAMOMS_RRSBIN(N,0,B,WI) * RRS_PHASMOM2
          ENDDO

!  Counting
!    Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

          NQ = NVARS(N)
          IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
          IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
          IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type)
!    Derivatives can be normalized or un-normalized.

          DO Q = 1, NQ
            RATIO_1  = - L_DELTAU_VERT(Q,N,WO) / DELTAU_VERT(N,WO)
            L_ADJUST = ADJUST * RATIO_1
            DO B = 1, N_RRSBINS(WI)
             L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) = L_ADJUST*BINXSEC(N,WI,B)
             L_OMEGAMOMS_RRSBIN(Q,N,2,B,WI) = &
                   L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) * RRS_PHASMOM2
            ENDDO
          ENDDO
          Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_AIRPROFILE_WFS ) THEN
            Q = Q + 1
            NORM = ONE
            IF ( DO_NORMALIZED_WFS) NORM = LAYER_AIRCOLUMNS(N)
            TERM_1 = L_DELTAU_VERT(Q,N,WO) * ADJUST
            L_ADJUST = ( NORM - TERM_1 ) / DELTAU_VERT(N,WO)
            DO B = 1, N_RRSBINS(WI)
             L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) = L_ADJUST*BINXSEC(N,WI,B)
             L_OMEGAMOMS_RRSBIN(Q,N,2,B,WI) = &
                   L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) * RRS_PHASMOM2
            ENDDO
          ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

          IF ( DO_TEMPPROFILE_WFS ) THEN
            Q = Q + 1
            TERM_3   = L_DELTAU_VERT(Q,N,WO) * ADJUST
            TERM_2   = LAYER_AIRCOLUMNS_dT(N)
            IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
            L_ADJUST = ( TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
            NORM = ADJUST
            IF ( DO_NORMALIZED_WFS ) NORM = NORM * TEMPERATURES(N)
            DO B = 1, N_RRSBINS(WI)
             L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) = &
                L_ADJUST*BINXSEC(N,WI,B) + NORM * BINXSEC_dT(N,WI,B)
             L_OMEGAMOMS_RRSBIN(Q,N,2,B,WI) = &
                L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) * RRS_PHASMOM2
            ENDDO
          ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here.

          IF ( DO_TEMPSHIFT_WF ) THEN
            Q = Q + 1
            TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
            TERM_1 = TSHIFT * ADJUST
            TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N)
            TERM_3 = L_DELTAU_VERT(Q,N,WO) * ADJUST
            L_ADJUST = ( TERM_2 - TERM_3 ) / DELTAU_VERT(N,WO)
            DO B = 1, N_RRSBINS(WI)
             L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) = &
                L_ADJUST*BINXSEC(N,WI,B) + TERM_1 * BINXSEC_dT(N,WI,B)
             L_OMEGAMOMS_RRSBIN(Q,N,2,B,WI) = &
                L_OMEGAMOMS_RRSBIN(Q,N,0,B,WI) * RRS_PHASMOM2
            ENDDO
          ENDIF

!  All first moment entries are zero

          DO Q = 1, NVARS(N)
            DO B = 1, N_RRSBINS(WI)
              L_OMEGAMOMS_RRSBIN(Q,N,1,B,WI) = ZERO
            ENDDO
          ENDDO

!  End layer loop

        ENDDO

!  End loop over all inner wavelengths

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_BIN_PLUS_2

!  Linearized MONO code

      SUBROUTINE LRRS_RAMAN_OPS_MONO_PLUS_1 &
        ( MAX_LAYERS, MAX_POINTS, MAX_GEOMETRIES, MAX_VARS,                          & ! Input dimensions
          DO_UPWELLING, DO_DNWELLING, DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,        & ! Input flags
          DO_NORMALIZED_WFS, DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF, & ! Input flags
          NLAYERS, NGEOMETRIES, NPOINTS_MONO, W_EXCIT, NVARS,                        & ! Inputs
          COSSCAT_UP, COSSCAT_DN, TEMPERATURES, TEMPERATURES_UNSHIFTED,              & ! Inputs
          LAYER_AIRCOLUMNS, LAYER_AIRCOLUMNS_DT,                                     & ! Inputs
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, DELTAU_VERT, L_DELTAU_VERT,                 & ! Inputs Optical
          OMEGAPHASFUNC_ELASTIC_UP, L_OMEGAPHASFUNC_ELASTIC_UP,                      & ! Inputs Optical
          OMEGAPHASFUNC_ELASTIC_DN, L_OMEGAPHASFUNC_ELASTIC_DN,                      & ! Inputs Optical
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO, RRSXSEC_RANKED_DT, RRSXSEC_OUT_MONO_DT,  & ! Inputs
          OMEGAPHASFUNC_CABANNES_UP, L_OMEGAPHASFUNC_CABANNES_UP,                    & ! Outputs
          OMEGAPHASFUNC_CABANNES_DN, L_OMEGAPHASFUNC_CABANNES_DN,                    & ! Outputs
          OMEGAMOMS_RRSLOSS,  L_OMEGAMOMS_RRSLOSS,                                   & ! Output
          OMEGAMOMS_RRSGAIN,  L_OMEGAMOMS_RRSGAIN )                                    ! Outputs

       IMPLICIT NONE

!  Notes
!  =====

!  IOPS FOR THE LRRS MODEL, MONOCHROMATIC REALIZATION
!  **************************************************

!  The  binning arrays are not used here.

!  IOPs are created for a complete calculation at one excitation wavelen
!  indexed by W_EXCIT in the ranked array of wavelengths.
!     W_EXCIT = 115 out of a total of 234 points = NPOINTS_MONO.
!  There are 233 Raman shifted wavelengths.

!  Given a set of optical properties defined for all 234 points which
!  contribute to the inelastic scattering for the W_EXCIT wavelength,
!  determine the necessary optical properties for the RRS model over
!  this window which has total number of points NPOINTS_MONO.

!  The source terms OMEGAMOMS_RRSGAIN (for light scattered into the
!  excitation wavelength), and OMEGAMOMS_RRSLOSS (light inelastically
!  scattered out), the latter term only considered for the Energy Balanc
!  Method.

!  Here is the indexing system for this subroutine
!    n  is the layer index,       n  = 1, .... NLAYERS
!    w  is the wavelength index,  w  = 1, ...  NPOINTS_MONO
!    V  is the geometry index,    v  = 0, .... NGEOMETRIES
!    q  is the Wtfunction index,  q  = 1, .... NVARS
!
!  The control integers for layers and geometries are
!
!      NLAYERS     = actual number of atmospheric layers
!      NGEOMETRIES = actual number of phase function geometries
!
!  There are scattering angle inputs to go with the phase function inputs
!
!      COSSCAT_UP (v) = Upwelling scattering angle cosines
!      COSSCAT_DN (v) = Dnwelling scattering angle cosines
!
!  The already-computed optical properties are total quantities defined
!  for all elastic scatterers for the total number of layers NLAYERS.
!  The optical properties are :
!
!      DELTAU_VERT(n,w)                = total optical depth for extinction.
!      OMEGAPHASFUNC_ELASTIC_UP(n,v,w) = SSA x Elastic upwelling   phase function
!      OMEGAPHASFUNC_ELASTIC_DN(n,v,w) = SSA x Elastic downwelling phase function
!
!  Also required as input for this process are:
!
!      TEMPERATURES(n)     = layer average temperatures, deg K
!      LAYER_AIRCOLUMNS(n) = column air density, [mol/cm^2]
!      RAYLEIGH_XSEC(w)    = Rayleigh cross-sections on complete grid [c
!      RAYLEIGH_DEPOL(w)   = Rayleigh depolarization ratios on complete
!
!  The scattering outputs are:
!
!      OMEGAPHASFUNC_CABANNES_UP(n,v,w) = SSA x Cabannes-adjusted upwelling   phase function
!      OMEGAPHASFUNC_CABANNES_DN(n,v,w) = SSA x Cabannes-adjusted downwelling phase function
!      OMEGAMOMS_RRSGAIN(n,l,w)  = SSA x RRS phasemoments. L = 0, 1, 2
!                                   Raman Gain terms for each bin.
!      OMEGAMOMS_RRSLOSS(n,l,1)  = SSA x RRS phasemoments. L = 0, 1, 2.
!                                   Raman loss term for the iner point

!      L_OMEGAPHASFUNC_CABANNES_UP(n,v,wi) = Linearized SSA x Cabannes-adjusted upwelling   phasefunction
!      L_OMEGAPHASFUNC_CABANNES_DN(n,v,wi) = Linearized SSA x Cabannes-adjusted downwelling phasefunction
!      L_OMEGAMOMS_RRSBIN(q,n,l,b,wi) = Linearized Raman Gain terms for each bin
!      L_OMEGAMOMS_RRSLOSS(q,n,l,wi)  = Linearized Raman loss term for the inner point


!  Input arguments
!  ===============

!  Maximum number of spectral points

      INTEGER  , INTENT(IN) :: MAX_POINTS

!  Maximum number of layers and geometries, weighting functions
!   -- Rob mod 5/12/17 for 2p5a, add MAX_GEOMETRIES, remove max_moments

      INTEGER  , INTENT(IN) :: MAX_LAYERS, MAX_GEOMETRIES, MAX_VARS

!  flags
!  -----

!  Directional flags
!   -- Rob mod 5/12/17 for 2p5a, added

      LOGICAL  , INTENT(IN) :: DO_UPWELLING, DO_DNWELLING

!  Method

      LOGICAL, INTENT(IN) :: DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Temperature and Air Profile Linearization flags

      LOGICAL, INTENT(IN) :: DO_TEMPPROFILE_WFS
      LOGICAL, INTENT(IN) :: DO_AIRPROFILE_WFS

!  Temperature shift Linearization flag (Column WF only)

      LOGICAL, INTENT(IN) :: DO_TEMPSHIFT_WF

!  use of normalized weighting functions (Profile WFs only)
!   T-shift Jacobian is automatically normalized

      LOGICAL, INTENT(IN) :: DO_NORMALIZED_WFS

!  Numbers
!  -------

!  Layering

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of Geometries
!   -- Rob mod 5/12/17 for 2p5a, add NGEOMETRIES, remove NMOMENTS

      INTEGER  , INTENT(IN) :: NGEOMETRIES

!  Number of points nd excitation position

      INTEGER  , INTENT(IN) :: NPOINTS_MONO, W_EXCIT

!  Number of weighting functions

      INTEGER  , INTENT(IN) :: NVARS ( MAX_LAYERS )

!  Scattering angles
!   -- Rob mod 5/12/17 for 2p5a, added

      REAL(FPK), INTENT(IN) :: COSSCAT_UP ( MAX_GEOMETRIES )
      REAL(FPK), INTENT(IN) :: COSSCAT_DN ( MAX_GEOMETRIES )

!  Layer temperatures (actual and unshifted)

      REAL(FPK), INTENT(IN) :: TEMPERATURES           ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: TEMPERATURES_UNSHIFTED ( MAX_LAYERS )

!  Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )

!  Temperature derivatives of Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS_dT ( MAX_LAYERS )

!  Vertical optical depths

      REAL(FPK), INTENT(IN) :: DELTAU_VERT ( MAX_LAYERS, MAX_POINTS )

!  Input elastic scattering phase function products
!   -- Rob mod 5/12/17 for 2p5a, Replaces the phasmoms input (OMEGAMOMS_ELASTIC)

      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Linearization of optical properties
!   -- Rob mod 5/12/17 for 2p5a, Replaces the phasmoms input (OMEGAMOMS_ELASTIC)

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT ( MAX_VARS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_UP ( MAX_VARS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_DN ( MAX_VARS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!   Gain term cross-sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_RANKED ( MAX_LAYERS, MAX_POINTS )

!  Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT_MONO ( MAX_LAYERS )

!  T-derivatives of Gain term cross-sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_RANKED_dT ( MAX_LAYERS,MAX_POINTS )

!  T-derivatives Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT_MONO_dT ( MAX_LAYERS )

!  Output arguments
!  ================

!  Output Cabannes-adjusted scattering phase function products
!   -- Rob mod 5/12/17 for 2p5a, Replaces the phasmoms output (OMEGAMOMS_CABANNES)
!                                Only calculated for the Cabannes-Raman option

      REAL(FPK), INTENT(OUT) :: OMEGAPHASFUNC_CABANNES_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAPHASFUNC_CABANNES_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSGAIN ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Linearizations

      REAL(FPK), INTENT(OUT) :: L_OMEGAPHASFUNC_CABANNES_UP ( MAX_VARS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAPHASFUNC_CABANNES_DN ( MAX_VARS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSLOSS ( MAX_VARS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSGAIN ( MAX_VARS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Local variables
!  ===============

      INTEGER   :: N, Q, NQ, V, W, W1
      REAL(FPK) :: DIF0(MAX_LAYERS), DIF2(MAX_LAYERS), ADJUST, RAMANPHASFUNC, NORM, NORM2
      REAL(FPK) :: L_DIF0(MAX_VARS,MAX_LAYERS), L_DIF2(MAX_VARS,MAX_LAYERS), L_ADJUST, RATIO_1
      REAL(FPK) :: SIGRAY, DEPOL, LOSSM0(MAX_LAYERS), LOSSM2(MAX_LAYERS)
      REAL(FPK) :: TSHIFT, LOSSM0_dT, LOSSM2_dT
      REAL(FPK) :: RAYMOM, CABMOM, RRS_PHASMOM2, AIRCOL, L_RAMANPHASFUNC
      REAL(FPK) :: TERM_1, TERM_2, TERM_3
      REAL(FPK) :: MU, LEG2_UP(MAX_GEOMETRIES), LEG2_DN(MAX_GEOMETRIES)

!  Initialize

      OMEGAMOMS_RRSGAIN  = zero
      OMEGAMOMS_RRSLOSS  = zero
      OMEGAPHASFUNC_CABANNES_UP = zero
      OMEGAPHASFUNC_CABANNES_DN = zero
      L_OMEGAMOMS_RRSGAIN  = zero
      L_OMEGAMOMS_RRSLOSS  = zero
      L_OMEGAPHASFUNC_CABANNES_UP = zero
      L_OMEGAPHASFUNC_CABANNES_DN = zero

!  Legendre calculations
!   -- Rob mod 5/12/17 for 2p5a, need Legendre polynomials for Cabannes adjustment

      IF ( DO_CABANNES_RAMAN ) THEN
         DO V = 1, NGEOMETRIES
           MU = COSSCAT_UP(V) ; LEG2_UP(V) = 1.5_fpk * MU * MU - 0.5_fpk
           MU = COSSCAT_DN(V) ; LEG2_DN(V) = 1.5_fpk * MU * MU - 0.5_fpk
         ENDDO
      ENDIF

!  Cabannes adjustment
!  -------------------

      IF ( DO_CABANNES_RAMAN ) THEN

!    Essentially adjusting the Rayleigh scattering part to make sure tha
!    the elastic component is just the due to Cabannes scattering
!    Only requires the first 0,1,2 moments to be adjusted.

!  We only need to do this for the excitation wavelength.

!  Rayleigh and Cabannes  Legendre(2) moments

         W1 = W_EXCIT
         DEPOL  = RAYLEIGH_DEPOL(W1)
         SIGRAY = RAYLEIGH_XSEC (W1)
         CABMOM = 0.25D0 * &
             ( 8.0D0 - 9.0D0*DEPOL ) / ( 4.0D0 - 3.0D0*DEPOL )
         RAYMOM = ( one - DEPOL ) / ( 2.0D0 + DEPOL )

!   -- Rob mod 5/12/17 for 2p5a, Adjust the phase functions for Cabannes scattering

!             DIF0 = RRS Loss-scattering optical depth / total optical depth
!             DIF2 = scatter-weighted RRS Loss / total OD, for moment L = 2

         DO N = 1, NLAYERS
           AIRCOL = LAYER_AIRCOLUMNS(N)
           LOSSM0(N) = RRSXSEC_OUT_MONO(N)
           LOSSM2(N) = SIGRAY*RAYMOM - (SIGRAY - LOSSM0(N)) * CABMOM
           ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,W1)
           DIF0(N)   = LOSSM0(N) * ADJUST
           DIF2(N)   = LOSSM2(N) * ADJUST
         ENDDO

!  Now make the adjustments

         IF ( DO_UPWELLING ) THEN
           DO V = 1, NGEOMETRIES
             DO N = 1, NLAYERS
               RAMANPHASFUNC = DIF0(N) + DIF2(N) * LEG2_UP(V)
               OMEGAPHASFUNC_CABANNES_UP(N,V,1) = OMEGAPHASFUNC_ELASTIC_UP(N,V,W1) - RAMANPHASFUNC
             ENDDO
           ENDDO
         ENDIF
         IF ( DO_DNWELLING ) THEN
           DO V = 1, NGEOMETRIES
             DO N = 1, NLAYERS
               RAMANPHASFUNC = DIF0(N) + DIF2(N) * LEG2_DN(V)
               OMEGAPHASFUNC_CABANNES_DN(N,V,1) = OMEGAPHASFUNC_ELASTIC_DN(N,V,W1) - RAMANPHASFUNC
             ENDDO
           ENDDO
         ENDIF

!  Here is the original adjustment of the moments
!     Cabannes L = 0, adjust OMEGAMOMS by DIF0
!     Cabannes L = 2, adjust OMEGAMOMS by DIF2
!     Cabannes L = 1,3,4.... no adjustment (just equal to elastic values)
!         DO N = 1, NLAYERS
!           OMEGAMOMS_CABANNES(N,0,1) = OMEGAMOMS_ELASTIC(N,0,W1) - DIF0(N)
!           OMEGAMOMS_CABANNES(N,1,1) = OMEGAMOMS_ELASTIC(N,1,W1)
!           OMEGAMOMS_CABANNES(N,2,1) = OMEGAMOMS_ELASTIC(N,2,W1) - DIF2(N)
!           OMEGAMOMS_CABANNES (N,3:NMOMENTS,1) = OMEGAMOMS_ELASTIC(N,3:NMOMENTS,WO)
!         ENDDO

!  Start layer loop

         DO N = 1, NLAYERS

!  Counting
!    Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

           NQ = NVARS(N)
           IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
           IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
           IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type)
!    Derivatives can be normalized or un-normalized.

           DO Q = 1, NQ
             RATIO_1 = - L_DELTAU_VERT(Q,N,W1) / DELTAU_VERT(N,W1)
             L_DIF0(Q,N)   = DIF0(N) * RATIO_1
             L_DIF2(Q,N)   = DIF2(N) * RATIO_1
           ENDDO
           Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_AIRPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = one
             IF ( DO_NORMALIZED_WFS) NORM = LAYER_AIRCOLUMNS(N)
             TERM_1 = L_DELTAU_VERT(Q,N,W1)  * DIF0(N)
             L_DIF0(Q,N) = (LOSSM0(N) * NORM - TERM_1 ) / DELTAU_VERT(N,W1)
             TERM_1 = L_DELTAU_VERT(Q,N,W1)  * DIF2(N)
             L_DIF2(Q,N) = (LOSSM2(N) * NORM - TERM_1 ) / DELTAU_VERT(N,W1)
           ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_TEMPPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = RRSXSEC_OUT_MONO_dT(N)
             IF ( DO_NORMALIZED_WFS ) NORM = NORM * TEMPERATURES(N)
             NORM2 = NORM * CABMOM
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF0(N)
             TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM0(N)
             TERM_1 = NORM * LAYER_AIRCOLUMNS(N)
             L_DIF0(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF2(N)
             TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM2(N)
             TERM_1 = NORM2 * LAYER_AIRCOLUMNS(N)
             L_DIF2(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
           ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here. 

           IF ( DO_TEMPSHIFT_WF ) THEN
             Q = Q + 1
             TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
             LOSSM0_dT = RRSXSEC_OUT_MONO_dT(N)
             LOSSM2_dT = LOSSM0_dT * CABMOM
             TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM0_dT
             TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM0(N)
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF0(N)
             L_DIF0(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
             TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM2_dT
             TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM2(N)
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF2(N)
             L_DIF2(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
           ENDIF

!  End saved layer loop

         ENDDO

!  Now make the adjustments

         IF ( DO_UPWELLING ) THEN
           DO V = 1, NGEOMETRIES
             DO N = 1, NLAYERS
               DO Q = 1, NVARS(N)
                 L_RAMANPHASFUNC = L_DIF0(Q,N) + L_DIF2(Q,N) * LEG2_UP(V)
                 L_OMEGAPHASFUNC_CABANNES_UP(Q,N,V,1) = L_OMEGAPHASFUNC_ELASTIC_UP(Q,N,V,W1) - L_RAMANPHASFUNC
               ENDDO
             ENDDO
           ENDDO
         ENDIF
         IF ( DO_DNWELLING ) THEN
           DO V = 1, NGEOMETRIES
             DO N = 1, NLAYERS
               DO Q = 1, NVARS(N)
                 L_RAMANPHASFUNC = L_DIF0(Q,N) + L_DIF2(Q,N) * LEG2_DN(V)
                 L_OMEGAPHASFUNC_CABANNES_DN(Q,N,V,1) = L_OMEGAPHASFUNC_ELASTIC_DN(Q,N,V,W1) - L_RAMANPHASFUNC
               ENDDO
             ENDDO
           ENDDO
         ENDIF

!  Here is the original adjustment of the moments
!     Cabannes L = 0, adjust OMEGAMOMS by DIF0
!     Cabannes L = 2, adjust OMEGAMOMS by DIF2
!     Cabannes L = 1,3,4.... no adjustment (just equal to elastic values)
!         DO N = 1, NLAYERS
!            DO Q = 1, NVARS(N)
!               L_OMEGAMOMS_CABANNES(Q,N,0,WI) = L_OMEGAMOMS_ELASTIC(Q,N,0,WO) - L_DIF0
!               L_OMEGAMOMS_CABANNES(Q,N,1,WI) = L_OMEGAMOMS_ELASTIC(Q,N,1,WO)
!               L_OMEGAMOMS_CABANNES(Q,N,2,WI) = L_OMEGAMOMS_ELASTIC(Q,N,2,WO) - L_DIF2
!               L_OMEGAMOMS_CABANNES(Q,N,3:NMOMENTS,WI) = L_OMEGAMOMS_ELASTIC(Q,N,3:NMOMENTS,WO)
!            ENDDO
!         ENDDO

!  End Cabannes-Raman treatment

      ENDIF

!  Loss adjustment
!  ---------------

!  In the energy balancing mode, need the loss terms
!    Only required at the excitation wavelength in the monochromatic cas

      IF ( DO_ENERGY_BALANCING ) THEN

         RRS_PHASMOM2 = 0.05D0
         W1 = W_EXCIT

!  Start layer loop

         DO N = 1, NLAYERS

!    Only required at the excitation wavelength in the monochromatic cas

           ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,W1)
           LOSSM0(N) = RRSXSEC_OUT_MONO(N)
           DIF0(N)   = LOSSM0(N) * ADJUST
           OMEGAMOMS_RRSLOSS(N,0,1) = DIF0(N)
           OMEGAMOMS_RRSLOSS(N,1,1) = zero
           OMEGAMOMS_RRSLOSS(N,2,1) = DIF0(N) * RRS_PHASMOM2

!  Counting
!    Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

           NQ = NVARS(N)
           IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
           IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
           IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type)
!    Derivatives can be normalized or un-normalized.

           DO Q = 1, NQ
             RATIO_1 = - L_DELTAU_VERT(Q,N,W1) / DELTAU_VERT(N,W1)
             L_DIF0(Q,N)   = DIF0(N) * RATIO_1
             L_OMEGAMOMS_RRSLOSS(Q,N,0,1) = L_DIF0(Q,N)
             L_OMEGAMOMS_RRSLOSS(Q,N,2,1) = L_DIF0(Q,N) * RRS_PHASMOM2
           ENDDO
           Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG
!    Robfix 06/01/11 @@@@. AIRCOL not defined.

           IF ( DO_AIRPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = one
! Wrong      IF ( DO_NORMALIZED_WFS) NORM = AIRCOL
             IF ( DO_NORMALIZED_WFS) NORM = LAYER_AIRCOLUMNS(N)
             TERM_1 = L_DELTAU_VERT(Q,N,W1)  * DIF0(N)
             L_DIF0(Q,N) = (LOSSM0(N) * NORM - TERM_1 ) / DELTAU_VERT(N,W1)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,1) = L_DIF0(Q,N)
             L_OMEGAMOMS_RRSLOSS(Q,N,2,1) = L_DIF0(Q,N) * RRS_PHASMOM2
           ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_TEMPPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = RRSXSEC_OUT_MONO_dT(N)
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF0(N)
             TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM0(N)
             TERM_1 = NORM * LAYER_AIRCOLUMNS(N)
             IF ( DO_NORMALIZED_WFS ) TERM_1 = TERM_1 * TEMPERATURES(N)
             IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
             L_DIF0(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,1) = L_DIF0(Q,N)
             L_OMEGAMOMS_RRSLOSS(Q,N,2,1) = L_DIF0(Q,N) * RRS_PHASMOM2
           ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here.

           IF ( DO_TEMPSHIFT_WF ) THEN
             Q = Q + 1
             TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
             LOSSM0_dT = RRSXSEC_OUT_MONO_dT(N)
!mick fix 10/19/2015 - line not needed for EB
             !LOSSM2_dT = LOSSM0_dT * CABMOM
             TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM0_dT
             TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM0(N)
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF0(N)
             L_DIF0(Q,N) = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,1) = L_DIF0(Q,N)
             L_OMEGAMOMS_RRSLOSS(Q,N,2,1) = L_DIF0(Q,N) * RRS_PHASMOM2
           ENDIF

!  First moment linearization is zero, always

           DO Q = 1, NVARS(N)
              L_OMEGAMOMS_RRSLOSS(Q,N,1,1) = zero
           ENDDO

!  End layer loop

         ENDDO

!  End loss term

      ENDIF

!  Gain terms for Inelastic scattering: RRS input
!  ----------------------------------------------

!  Calculate optical property input for the Gain RRS source terms
!  The second moment for RRS scattering = 1/20.

      RRS_PHASMOM2 = 0.05D0
      DO W = 1, NPOINTS_MONO

!  Start layer loop

        DO N = 1, NLAYERS

!  The adjustment = aircolumn / total optical thickness
!  The gain term for the excitation wavelength is set to zero

          AIRCOL = LAYER_AIRCOLUMNS(N)
          ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,W1)
          OMEGAMOMS_RRSGAIN(N,0,W) = ADJUST * RRSXSEC_RANKED(N,W)
          OMEGAMOMS_RRSGAIN(N,1,W) = zero
          OMEGAMOMS_RRSGAIN(N,2,W) = OMEGAMOMS_RRSGAIN(N,0,W) * RRS_PHASMOM2

!  Counting
!    Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

          NQ = NVARS(N)
          IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
          IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
          IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type)
!    Derivatives can be normalized or un-normalized.

          DO Q = 1, NQ
            RATIO_1  = - L_DELTAU_VERT(Q,N,W1) / DELTAU_VERT(N,W1)
            L_ADJUST = ADJUST * RATIO_1
            L_OMEGAMOMS_RRSGAIN(Q,N,0,W) = L_ADJUST*RRSXSEC_RANKED(N,W)
            L_OMEGAMOMS_RRSGAIN(Q,N,2,W) = L_OMEGAMOMS_RRSGAIN(Q,N,0,W) * RRS_PHASMOM2
          ENDDO
          Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG

          IF ( DO_AIRPROFILE_WFS ) THEN
            Q = Q + 1
            NORM = one
            IF ( DO_NORMALIZED_WFS) NORM = AIRCOL
            TERM_1 = L_DELTAU_VERT(Q,N,W1) * ADJUST
            L_ADJUST = ( NORM - TERM_1 ) / DELTAU_VERT(N,W1)
            L_OMEGAMOMS_RRSGAIN(Q,N,0,W) = L_ADJUST*RRSXSEC_RANKED(N,W)
            L_OMEGAMOMS_RRSGAIN(Q,N,2,W) = L_OMEGAMOMS_RRSGAIN(Q,N,0,W) * RRS_PHASMOM2
          ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_TEMPPROFILE_WFS ) THEN
            Q = Q + 1
            TERM_3   = L_DELTAU_VERT(Q,N,W1) * ADJUST
            TERM_2   = LAYER_AIRCOLUMNS_dT(N)
            IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
            L_ADJUST = ( TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
            NORM = ADJUST
            IF ( DO_NORMALIZED_WFS ) NORM = NORM * TEMPERATURES(N)
            L_OMEGAMOMS_RRSGAIN(Q,N,0,W) = L_ADJUST*RRSXSEC_RANKED(N,W) &
                                         + NORM * RRSXSEC_RANKED_dT(N,W)
            L_OMEGAMOMS_RRSGAIN(Q,N,2,W) = L_OMEGAMOMS_RRSGAIN(Q,N,0,W) * RRS_PHASMOM2
          ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here.

           IF ( DO_TEMPSHIFT_WF ) THEN
            Q = Q + 1
            TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
            TERM_3   = L_DELTAU_VERT(Q,N,W1) * ADJUST
            TERM_2   = TSHIFT * LAYER_AIRCOLUMNS_dT(N)
            L_ADJUST = ( TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
            NORM = TSHIFT * ADJUST
            L_OMEGAMOMS_RRSGAIN(Q,N,0,W) = L_ADJUST*RRSXSEC_RANKED(N,W) &
                                         + NORM * RRSXSEC_RANKED_dT(N,W)
            L_OMEGAMOMS_RRSGAIN(Q,N,2,W) = L_OMEGAMOMS_RRSGAIN(Q,N,0,W) * RRS_PHASMOM2
           ENDIF

!  All first moment entries are zero

          DO Q = 1, NVARS(N)
            L_OMEGAMOMS_RRSGAIN(Q,N,1,W) = zero
          ENDDO

!  End layer loop

        ENDDO

!  End loop over points

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_MONO_PLUS_1

!

      SUBROUTINE LRRS_RAMAN_OPS_MONO_PLUS_2 &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, MAX_VARS,             & ! Inputs
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, DO_NORMALIZED_WFS, & ! Inputs
          DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF,    & ! Inputs
          NLAYERS, NMOMENTS_INPUT, NPOINTS_MONO, W_EXCIT, NVARS,     & ! Inputs
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL,                             & ! Inputs
          TEMPERATURES, TEMPERATURES_UNSHIFTED,                      & ! Inputs
          LAYER_AIRCOLUMNS, LAYER_AIRCOLUMNS_DT,                     & ! Inputs
          DELTAU_VERT,     OMEGAMOMS_ELASTIC,                        & ! Inputs
          L_DELTAU_VERT, L_OMEGAMOMS_ELASTIC,                        & ! Inputs
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO,                          & ! Inputs
          RRSXSEC_RANKED_DT, RRSXSEC_OUT_MONO_DT,                    & ! Inputs
          OMEGAMOMS_CABANNES, L_OMEGAMOMS_CABANNES,                  & ! Outputs
          OMEGAMOMS_RRSLOSS,  L_OMEGAMOMS_RRSLOSS,                   & ! Outputs
          OMEGAMOMS_RRSGAIN,  L_OMEGAMOMS_RRSGAIN )                    ! Outputs

       IMPLICIT NONE

!  Notes
!  =====

!  IOPS FOR THE LRRS MODEL, MONOCHROMATIC REALIZATION
!  **************************************************

!  The  binning arrays are not used here.

!  IOPs are created for a complete calculation at one excitation wavelen
!  indexed by W_EXCIT in the ranked array of wavelengths.
!     W_EXCIT = 115 out of a total of 234 points = NPOINTS_MONO.
!  There are 233 Raman shifted wavelengths.

!  Given a set of optical properties defined for all 234 points which
!  contribute to the inelastic scattering for the W_EXCIT wavelength,
!  determine the necessary optical properties for the RRS model over
!  this window which has total number of points NPOINTS_MONO.

!  The source terms OMEGAMOMS_RRSGAIN (for light scattered into the
!  excitation wavelength), and OMEGAMOMS_RRSLOSS (light inelastically
!  scattered out), the latter term only considered for the Energy-Balancing Method.

!  Here is the indexing system for this subroutine
!    n  is the layer index,       n  = 1, .... NLAYERS
!    w  is the wavelength index,  w  = 1, ...  NPOINTS_MONO
!    l  is the Leg. moment index, l  = 0, .... NMOMENTS_INPUT
!
!  The control INTEGERs for layers and moments are
!
!      NLAYERS        = actual number of atmospheric layers
!      NMOMENTS_INPUT = actual number of phase function moments
!
!  The already-computed optical properties are total quantities defined
!  for all elastic scatterers for the total number of layers NLAYERS.
!  The optical properties are :
!
!      DELTAU_VERT_INPUT(n,w)   = total optical depth for extinction.
!      OMEGAMOMS_ELASTIC(n,l,w) = single scatter albedo X phase-moment.

!    n  is the layer index,       n  = 1, .... NLAYERS
!    w  is the wavelength index,  w  = 1, ...  NPOINTS_MONO
!    l  is the Leg. moment index, l  = 0, .... NMOMENTS_INPUT
!    q  is the Wtfunction index,  q  = 1, .... NVARS
!
!  Also required as input for this process are:
!
!      TEMPERATURES(n)     = layer average temperatures, deg K
!      LAYER_AIRCOLUMNS(n) = column air density, [mol/cm^2]
!      RAYLEIGH_XSEC(w)    = Rayleigh cross-sections on complete grid [c
!      RAYLEIGH_DEPOL(w)   = Rayleigh depolarization ratios on complete
!
!  The scattering outputs are:
!
!      OMEGAMOMS_CABANNES(n,l,1) = SSA x Cabannes-adjusted phasemoments
!                                   (l = 0, 1, ...NMOMENTS_TOTAL)
!      OMEGAMOMS_RRSGAIN(n,l,w)  = SSA x RRS phasemoments. L = 0, 1, 2
!                                   Raman Gain terms for each bin.
!      OMEGAMOMS_RRSLOSS(n,l,1)  = SSA x RRS phasemoments. L = 0, 1, 2.
!                                   Raman loss term for the iner point

!  Input arguments
!  ===============

!  Maximum number of spectral points

      INTEGER  , INTENT(IN) :: MAX_POINTS

!  Maximum number of layers, input moments, weighting functions

      INTEGER  , INTENT(IN) :: MAX_LAYERS, MAX_MOMENTS, MAX_VARS

!  flags
!  -----

!  Method

      LOGICAL, INTENT(IN) :: DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Temperature and Air Profile Linearization flags

      LOGICAL, INTENT(IN) :: DO_TEMPPROFILE_WFS
      LOGICAL, INTENT(IN) :: DO_AIRPROFILE_WFS

!  Temperature shift Linearization flag (Column WF only)

      LOGICAL, INTENT(IN) :: DO_TEMPSHIFT_WF

!  use of normalized weighting functions (Profile WFs only)
!   T-shift Jacobian is automatically normalized

      LOGICAL, INTENT(IN) :: DO_NORMALIZED_WFS

!  Numbers
!  -------

!  Layering

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of elastic moments

      INTEGER  , INTENT(IN) :: NMOMENTS_INPUT

!  Number of points nd excitation position

      INTEGER  , INTENT(IN) :: NPOINTS_MONO, W_EXCIT

!  Number of weighting functions

      INTEGER  , INTENT(IN) :: NVARS ( MAX_LAYERS )

!  Layer temperatures (actual and unshifted)

      REAL(FPK), INTENT(IN) :: TEMPERATURES           ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: TEMPERATURES_UNSHIFTED ( MAX_LAYERS )

!  Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )

!  Temperature derivatives of Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS_dT ( MAX_LAYERS )

!  Vertical optical depths

      REAL(FPK), INTENT(IN) :: DELTAU_VERT ( MAX_LAYERS, MAX_POINTS )

!  Input elastic scattering law

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Linearization of optical properties

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT       ( MAX_VARS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_ELASTIC ( MAX_VARS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!   Gain term cross-sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_RANKED ( MAX_LAYERS, MAX_POINTS )

!  Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT_MONO ( MAX_LAYERS )

!  T-derivatives of Gain term cross-sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_RANKED_dT ( MAX_LAYERS,MAX_POINTS )

!  T-derivatives Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT_MONO_dT ( MAX_LAYERS )

!  Output arguments
!  ================

!  Scattering quantities

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSGAIN ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Linearizations

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_CABANNES ( MAX_VARS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSLOSS ( MAX_VARS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSGAIN ( MAX_VARS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Local variables
!  ===============

      INTEGER   :: W1, W, N, L, Q, NQ
      REAL(FPK) :: DIF0, DIF2, ADJUST, NORM, NORM2
      REAL(FPK) :: L_DIF0, L_DIF2, L_ADJUST, RATIO_1
      REAL(FPK) :: SIGRAY, DEPOL, LOSSM0, LOSSM2
      REAL(FPK) :: TSHIFT, LOSSM0_dT, LOSSM2_dT
      REAL(FPK) :: RAYMOM, CABMOM, RRS_PHASMOM2, AIRCOL
      REAL(FPK) :: TERM_1, TERM_2, TERM_3

!  Initialize

      OMEGAMOMS_RRSGAIN  = zero
      OMEGAMOMS_RRSLOSS  = zero
      OMEGAMOMS_CABANNES = zero
      L_OMEGAMOMS_RRSGAIN  = zero
      L_OMEGAMOMS_RRSLOSS  = zero
      L_OMEGAMOMS_CABANNES = zero

!  Cabannes adjustment
!  -------------------

      IF ( DO_CABANNES_RAMAN ) THEN

!    Essentially adjusting the Rayleigh scattering part to make sure tha
!    the elastic component is just the due to Cabannes scattering
!    Only requires the first 0,1,2 moments to be adjusted.

!  We only need to do this for the excitation wavelength.

!  Rayleigh and Cabannes  Legendre(2) moments

         W1 = W_EXCIT
         DEPOL  = RAYLEIGH_DEPOL(W1)
         SIGRAY = RAYLEIGH_XSEC (W1)
         CABMOM = 0.25D0 * &
             ( 8.0D0 - 9.0D0*DEPOL ) / ( 4.0D0 - 3.0D0*DEPOL )
         RAYMOM = ( one - DEPOL ) / ( 2.0D0 + DEPOL )

!  Start layer loop

         DO N = 1, NLAYERS

!  DIF0 = RRS Loss-scattering optical depth / total optical depth
!  DIF2 = scatter-weighted RRS Loss / total OD, for moment L = 2
!  Cabannes L = 0, adjust OMEGAMOMS by DIF0
!  Cabannes L = 1, no adjustment (just equal to elastic value)
!  Cabannes L = 2, adjust OMEGAMOMS by DIF2
!  Cabannes L = 3,4.. no adjustment (just equal to elastic values)

           AIRCOL = LAYER_AIRCOLUMNS(N)
           LOSSM0 = RRSXSEC_OUT_MONO(N)
           LOSSM2 = SIGRAY*RAYMOM - (SIGRAY - LOSSM0) * CABMOM
           ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,W1)
           DIF0   = LOSSM0 * ADJUST
           DIF2   = LOSSM2 * ADJUST
           OMEGAMOMS_CABANNES(N,0,1) = OMEGAMOMS_ELASTIC(N,0,W1) - DIF0
           OMEGAMOMS_CABANNES(N,1,1) = OMEGAMOMS_ELASTIC(N,1,W1)
           OMEGAMOMS_CABANNES(N,2,1) = OMEGAMOMS_ELASTIC(N,2,W1) - DIF2
           DO L = 3, NMOMENTS_INPUT
             OMEGAMOMS_CABANNES (N,L,1) = OMEGAMOMS_ELASTIC(N,L,W1)
           ENDDO

!  Counting
!    Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

           NQ = NVARS(N)
           IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
           IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
           IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type)
!    Derivatives can be normalized or un-normalized.

           DO Q = 1, NQ
             RATIO_1 = - L_DELTAU_VERT(Q,N,W1) / DELTAU_VERT(N,W1)
             L_DIF0   = DIF0 * RATIO_1
             L_DIF2   = DIF2 * RATIO_1
             L_OMEGAMOMS_CABANNES(Q,N,0,1) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,0,W1) - L_DIF0
             L_OMEGAMOMS_CABANNES(Q,N,2,1) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,2,W1) - L_DIF2
           ENDDO
           Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_AIRPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = one
             IF ( DO_NORMALIZED_WFS) NORM = AIRCOL
             TERM_1 = L_DELTAU_VERT(Q,N,W1)  * DIF0
             L_DIF0 = (LOSSM0 * NORM - TERM_1 ) / DELTAU_VERT(N,W1)
             TERM_1 = L_DELTAU_VERT(Q,N,W1)  * DIF2
             L_DIF2 = (LOSSM2 * NORM - TERM_1 ) / DELTAU_VERT(N,W1)
             L_OMEGAMOMS_CABANNES(Q,N,0,1) = &
                  L_OMEGAMOMS_ELASTIC(Q,N,0,W1) - L_DIF0
             L_OMEGAMOMS_CABANNES(Q,N,2,1) = &
                  L_OMEGAMOMS_ELASTIC(Q,N,2,W1) - L_DIF2
           ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_TEMPPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = RRSXSEC_OUT_MONO_dT(N)
             IF ( DO_NORMALIZED_WFS ) NORM = NORM * TEMPERATURES(N)
             NORM2 = NORM * CABMOM
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF0
             TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM0
             TERM_1 = NORM * LAYER_AIRCOLUMNS(N)
             L_DIF0 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF2
             TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM2
             TERM_1 = NORM2 * LAYER_AIRCOLUMNS(N)
             L_DIF2 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
             L_OMEGAMOMS_CABANNES(Q,N,0,1) = &
                  L_OMEGAMOMS_ELASTIC(Q,N,0,W1) - L_DIF0
             L_OMEGAMOMS_CABANNES(Q,N,2,1) = &
                  L_OMEGAMOMS_ELASTIC(Q,N,2,W1) - L_DIF2
           ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here. 

           IF ( DO_TEMPSHIFT_WF ) THEN
             Q = Q + 1
             TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
             LOSSM0_dT = RRSXSEC_OUT_MONO_dT(N)
             LOSSM2_dT = LOSSM0_dT * CABMOM
             TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM0_dT
             TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM0
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF0
             L_DIF0 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
             TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM2_dT
             TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM2
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF2
             L_DIF2 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
             L_OMEGAMOMS_CABANNES(Q,N,0,1) = &
                  L_OMEGAMOMS_ELASTIC(Q,N,0,W1) - L_DIF0
             L_OMEGAMOMS_CABANNES(Q,N,2,1) = &
                  L_OMEGAMOMS_ELASTIC(Q,N,2,W1) - L_DIF2
           ENDIF

!  1st moment, Elastic terms are all just copies

           DO Q = 1, NVARS(N)
             L_OMEGAMOMS_CABANNES(Q,N,1,1) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,1,W1)
             DO L = 3, NMOMENTS_INPUT
              L_OMEGAMOMS_CABANNES(Q,N,L,1) = &
                   L_OMEGAMOMS_ELASTIC(Q,N,L,W1)
             ENDDO
           ENDDO

!  End layer loop

         ENDDO

!  End Cabannes-Raman treatment

      ENDIF

!  Loss adjustment
!  ---------------

!  In the energy balancing mode, need the loss terms
!    Only required at the excitation wavelength in the monochromatic cas

      IF ( DO_ENERGY_BALANCING ) THEN

         RRS_PHASMOM2 = 0.05D0
         W1 = W_EXCIT

!  Start layer loop

         DO N = 1, NLAYERS

!    Only required at the excitation wavelength in the monochromatic cas

           ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,W1)
           LOSSM0 = RRSXSEC_OUT_MONO(N)
           DIF0   = LOSSM0 * ADJUST
           OMEGAMOMS_RRSLOSS(N,0,1) = DIF0
           OMEGAMOMS_RRSLOSS(N,1,1) = zero
           OMEGAMOMS_RRSLOSS(N,2,1) = DIF0 * RRS_PHASMOM2

!  Counting
!    Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

           NQ = NVARS(N)
           IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
           IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
           IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type)
!    Derivatives can be normalized or un-normalized.

           DO Q = 1, NQ
             RATIO_1 = - L_DELTAU_VERT(Q,N,W1) / DELTAU_VERT(N,W1)
             L_DIF0   = DIF0 * RATIO_1
             L_OMEGAMOMS_RRSLOSS(Q,N,0,1) = L_DIF0
             L_OMEGAMOMS_RRSLOSS(Q,N,2,1) = L_DIF0 * RRS_PHASMOM2
           ENDDO
           Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG
!    Robfix 06/01/11 @@@@. AIRCOL not defined.

           IF ( DO_AIRPROFILE_WFS ) THEN
             Q = Q + 1

             NORM = one
! Wrong      IF ( DO_NORMALIZED_WFS) NORM = AIRCOL
             IF ( DO_NORMALIZED_WFS) NORM = LAYER_AIRCOLUMNS(N)
             TERM_1 = L_DELTAU_VERT(Q,N,W1)  * DIF0
             L_DIF0 = (LOSSM0 * NORM - TERM_1 ) / DELTAU_VERT(N,W1)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,1) = L_DIF0
             L_OMEGAMOMS_RRSLOSS(Q,N,2,1) = L_DIF0 * RRS_PHASMOM2
           ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_TEMPPROFILE_WFS ) THEN
             Q = Q + 1
             NORM = RRSXSEC_OUT_MONO_dT(N)
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF0
             TERM_2 = LAYER_AIRCOLUMNS_dT(N) * LOSSM0
             TERM_1 = NORM * LAYER_AIRCOLUMNS(N)
             IF ( DO_NORMALIZED_WFS ) TERM_1 = TERM_1 * TEMPERATURES(N)
             IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
             L_DIF0 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,1) = L_DIF0
             L_OMEGAMOMS_RRSLOSS(Q,N,2,1) = L_DIF0 * RRS_PHASMOM2
           ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here.

           IF ( DO_TEMPSHIFT_WF ) THEN
             Q = Q + 1
             TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
             LOSSM0_dT = RRSXSEC_OUT_MONO_dT(N)
!mick fix 10/19/2015 - line not needed for EB
             !LOSSM2_dT = LOSSM0_dT * CABMOM
             TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM0_dT
             TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM0
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF0
             L_DIF0 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,1) = L_DIF0
             L_OMEGAMOMS_RRSLOSS(Q,N,2,1) = L_DIF0 * RRS_PHASMOM2
           ENDIF

!  First moment linearization is zero, always

           DO Q = 1, NVARS(N)
              L_OMEGAMOMS_RRSLOSS(Q,N,1,1) = zero
           ENDDO

!  End layer loop

         ENDDO

!  End loss term

      ENDIF

!  Gain terms for Inelastic scattering: RRS input
!  ----------------------------------------------

!  Calculate optical property input for the Gain RRS source terms
!  The second moment for RRS scattering = 1/20.

      RRS_PHASMOM2 = 0.05D0
      DO W = 1, NPOINTS_MONO

!  Start layer loop

        DO N = 1, NLAYERS

!  The adjustment = aircolumn / total optical thickness
!  The gain term for the excitation wavelength is set to zero

          AIRCOL = LAYER_AIRCOLUMNS(N)
          ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,W1)
          OMEGAMOMS_RRSGAIN(N,0,W) = ADJUST * RRSXSEC_RANKED(N,W)
          OMEGAMOMS_RRSGAIN(N,1,W) = zero
          OMEGAMOMS_RRSGAIN(N,2,W) = &
                   OMEGAMOMS_RRSGAIN(N,0,W) * RRS_PHASMOM2

!  Counting
!    Always do "Regular" weighting functions first. Then ---
!      For profiles, ORDER IS IMPORTANT: Air first, then Temperature
!      For Columns,  T-shift follows, if flagged

          NQ = NVARS(N)
          IF ( DO_AIRPROFILE_WFS  ) NQ = NQ - 1
          IF ( DO_TEMPPROFILE_WFS ) NQ = NQ - 1
          IF ( DO_TEMPSHIFT_WF )    NQ = NQ - 1

!  Linearized values (Regular Type)
!    Derivatives can be normalized or un-normalized.

          DO Q = 1, NQ
            RATIO_1  = - L_DELTAU_VERT(Q,N,W1) / DELTAU_VERT(N,W1)
            L_ADJUST = ADJUST * RATIO_1
            L_OMEGAMOMS_RRSGAIN(Q,N,0,W) = L_ADJUST*RRSXSEC_RANKED(N,W)
            L_OMEGAMOMS_RRSGAIN(Q,N,2,W) = &
                   L_OMEGAMOMS_RRSGAIN(Q,N,0,W) * RRS_PHASMOM2
          ENDDO
          Q = NQ

!  Linearization, w.r.t. Profile of Layer Air columns.
!    Derivatives can be normalized or un-normalized. Need to use FLAG

          IF ( DO_AIRPROFILE_WFS ) THEN
            Q = Q + 1
            NORM = one
            IF ( DO_NORMALIZED_WFS) NORM = AIRCOL
            TERM_1 = L_DELTAU_VERT(Q,N,W1) * ADJUST
            L_ADJUST = ( NORM - TERM_1 ) / DELTAU_VERT(N,W1)
            L_OMEGAMOMS_RRSGAIN(Q,N,0,W) = L_ADJUST*RRSXSEC_RANKED(N,W)
            L_OMEGAMOMS_RRSGAIN(Q,N,2,W) = &
                   L_OMEGAMOMS_RRSGAIN(Q,N,0,W) * RRS_PHASMOM2
          ENDIF

!  Linearization, w.r.t. Temperature Profile
!    Derivatives can be normalized or un-normalized. Need to use FLAG

           IF ( DO_TEMPPROFILE_WFS ) THEN
            Q = Q + 1
            TERM_3   = L_DELTAU_VERT(Q,N,W1) * ADJUST
            TERM_2   = LAYER_AIRCOLUMNS_dT(N)
            IF ( DO_NORMALIZED_WFS ) TERM_2 = TERM_2 * TEMPERATURES(N)
            L_ADJUST = ( TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
            NORM = ADJUST
            IF ( DO_NORMALIZED_WFS ) NORM = NORM * TEMPERATURES(N)
            L_OMEGAMOMS_RRSGAIN(Q,N,0,W) = L_ADJUST*RRSXSEC_RANKED(N,W) &
                                         + NORM * RRSXSEC_RANKED_dT(N,W)
            L_OMEGAMOMS_RRSGAIN(Q,N,2,W) = &
                  L_OMEGAMOMS_RRSGAIN(Q,N,0,W) * RRS_PHASMOM2
          ENDIF

!  Linearization, w.r.t. Temperature Shift (Column WF)
!    Derivatives Must be normalized here.

           IF ( DO_TEMPSHIFT_WF ) THEN
            Q = Q + 1
            TSHIFT = TEMPERATURES(N) - TEMPERATURES_UNSHIFTED(N)
            TERM_3   = L_DELTAU_VERT(Q,N,W1) * ADJUST
            TERM_2   = TSHIFT * LAYER_AIRCOLUMNS_dT(N)
            L_ADJUST = ( TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
            NORM = TSHIFT * ADJUST
            L_OMEGAMOMS_RRSGAIN(Q,N,0,W) = L_ADJUST*RRSXSEC_RANKED(N,W) &
                                         + NORM * RRSXSEC_RANKED_dT(N,W)
            L_OMEGAMOMS_RRSGAIN(Q,N,2,W) = &
                  L_OMEGAMOMS_RRSGAIN(Q,N,0,W) * RRS_PHASMOM2
           ENDIF

!  All first moment entries are zero

          DO Q = 1, NVARS(N)
            L_OMEGAMOMS_RRSGAIN(Q,N,1,W) = zero
          ENDDO

!  End layer loop

        ENDDO

!  End loop over points

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_MONO_PLUS_2

!   End module

      END MODULE lrrs_l_generate_ramanops_m

