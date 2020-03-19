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
! #         lrrs_raman_ops_bin        ( no  linearization)      #
! #         lrrs_raman_ops_mono       ( no  linearization)      #
! #                                                             #
! #  These subroutines calculate the Raman optical property     #
! #  inputs.                                                    #
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

      MODULE lrrs_generate_ramanops_m

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE

      PRIVATE
      PUBLIC :: LRRS_RAMAN_OPS_BIN_1,  LRRS_RAMAN_OPS_BIN_2,  &
                LRRS_RAMAN_OPS_MONO_1, LRRS_RAMAN_OPS_MONO_2

      CONTAINS

      SUBROUTINE LRRS_RAMAN_OPS_BIN_1 &
        ( MAX_LAYERS, MAX_POINTS, MAX_GEOMETRIES, MAX_BINS,                   & ! Input dimensions
          DO_UPWELLING, DO_DNWELLING, DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, & ! Input flags
          NLAYERS, NGEOMETRIES, NPOINTS_INNER, OFFSET_INNER,                  & ! Input integers
          COSSCAT_UP, COSSCAT_DN, LAYER_AIRCOLUMNS,              & ! Input 
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, DELTAU_VERT,            & ! Inputs Optical
          OMEGAPHASFUNC_ELASTIC_UP, OMEGAPHASFUNC_ELASTIC_DN,    & ! Inputs Optical
          N_RRSBINS, BINXSEC, RRSXSEC_OUT,                       & ! Inputs Raman Spec.
          OMEGAPHASFUNC_CABANNES_UP, OMEGAPHASFUNC_CABANNES_DN,  & ! Outputs Adjusted Cabannes
          OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN )                    ! Outputs raman Loss//Gain

       IMPLICIT NONE

!  Notes
!  =====

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
!
!  The already-computed optical properties are total quantities defined
!  for all elastic scatterers for the total number of layers NLAYERS.
!  The optical properties are :
!
!      DELTAU_VERT(n,wo)                = total optical depth for extinction.
!      OMEGAPHASFUNC_ELASTIC_UP(n,v,wo) = SSA x Elastic upwelling   phase function
!      OMEGAPHASFUNC_ELASTIC_DN(n,v,wo) = SSA x Elastic downwelling phase function

!  Also required as input for this process are:
!
!      LAYER_AIRCOLUMNS(n) = column air density, [mol/cm^2]
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

!  Input arguments
!  ===============

!  Maximum number of spectral points and bins

      INTEGER  , INTENT(IN) :: MAX_POINTS, MAX_BINS

!  Maximum number of layers and geometries
!   -- Rob mod 5/12/17 for 2p5a, add MAX_GEOMETRIES, remove max_moments

      INTEGER  , INTENT(IN) :: MAX_LAYERS, MAX_GEOMETRIES

!  Directional flags
!   -- Rob mod 5/12/17 for 2p5a, added

      LOGICAL  , INTENT(IN) :: DO_UPWELLING, DO_DNWELLING

!  LRRS Calculation Method

      LOGICAL  , INTENT(IN) :: DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Layering

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of Geometries
!   -- Rob mod 5/12/17 for 2p5a, add NGEOMETRIES, remove NMOMENTS

      INTEGER  , INTENT(IN) :: NGEOMETRIES

!  Number of points between the limiting wavelengths

      INTEGER  , INTENT(IN) :: NPOINTS_INNER

!  outer wavelength offset

      INTEGER  , INTENT(IN) :: OFFSET_INNER

!  Scattering angles
!   -- Rob mod 5/12/17 for 2p5a, added

      REAL(FPK), INTENT(IN) :: COSSCAT_UP ( MAX_GEOMETRIES )
      REAL(FPK), INTENT(IN) :: COSSCAT_DN ( MAX_GEOMETRIES )

!  Vertical optical depths

      REAL(FPK), INTENT(IN) :: DELTAU_VERT ( MAX_LAYERS, MAX_POINTS )

!  Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )

!  Input elastic scattering phase function products
!   -- Rob mod 5/12/17 for 2p5a, Replaces the phasmoms input (OMEGAMOMS_ELASTIC)

      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  Bin Mapping numbers

      INTEGER  , INTENT(IN) :: N_RRSBINS  ( MAX_POINTS )

!  Binned Gain term cross-sections

      REAL(FPK), INTENT(IN) :: BINXSEC ( MAX_LAYERS, MAX_POINTS, MAX_BINS )

!  Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT ( MAX_LAYERS, MAX_POINTS )

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

!  Local variables
!  ===============

      INTEGER   :: WI, WO, N, V, B
      REAL(FPK) :: DIF0(MAX_LAYERS), DIF2(MAX_LAYERS), ADJUST, RRS_PHASMOM2,RAMANPHASFUNC
      REAL(FPK) :: SIGRAY, SIGCAB, RAYMOM, CABMOM, DEPOL
      REAL(FPK) :: MU, LEG2_UP(MAX_GEOMETRIES), LEG2_DN(MAX_GEOMETRIES)

!  Initialize
!   -- Rob mod 5/12/17 for 2p5a, replaces the phasmoms output (OMEGAMOMS_CABANNES)

      OMEGAMOMS_RRSBIN   = zero
      OMEGAMOMS_RRSLOSS  = zero
      OMEGAPHASFUNC_CABANNES_UP = zero
      OMEGAPHASFUNC_CABANNES_DN = zero

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
!    Essentially adjusting the Rayleigh scattering part to make sure that
!    the elastic component is just that due to Cabannes scattering
!    Only requires the first 0,1,2 moments to be calculated

          DEPOL  = RAYLEIGH_DEPOL(WO)
          SIGRAY = RAYLEIGH_XSEC(WO)

!      Rayleigh and Cabannes Legendre(2) moments

          CABMOM = 0.25D0 * &
             ( 8.0D0 - 9.0D0*DEPOL ) / ( 4.0D0 - 3.0D0*DEPOL )
          RAYMOM = ( one - DEPOL ) / ( 2.0D0 + DEPOL )

!   -- Rob mod 5/12/17 for 2p5a, Adjust the phase functions for Cabannes scattering

!             DIF0 = RRS Loss-scattering optical depth / total optical depth
!             DIF2 = scatter-weighted RRS Loss / total OD, for moment L = 2

          DO N = 1, NLAYERS
            SIGCAB = SIGRAY - RRSXSEC_OUT(N,WI)
            ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,WO)
            DIF0(N)   = RRSXSEC_OUT(N,WI) * ADJUST
            DIF2(N)   = ( SIGRAY*RAYMOM - SIGCAB*CABMOM ) * ADJUST
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

!  End Cabannes-Raman treatment

        ENDIF

!  Loss term (Energy balancing mode only)
!  --------------------------------------

!  In the energy balancing mode, need the loss term IOPs-------------

        IF ( DO_ENERGY_BALANCING ) THEN
          RRS_PHASMOM2 = 0.05D0
          DO N = 1, NLAYERS
            ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,WO)
            DIF0(N)   = RRSXSEC_OUT(N,WI) * ADJUST
            OMEGAMOMS_RRSLOSS(N,0,WI) = DIF0(N)
            OMEGAMOMS_RRSLOSS(N,1,WI) = ZERO
            OMEGAMOMS_RRSLOSS(N,2,WI) = DIF0(N) * RRS_PHASMOM2
          ENDDO
        ENDIF

!  Gain terms (both implementations)
!  ---------------------------------

!  Inelastic scattering: RRS input---------------
!  Calculate optical property input for the Gain RRS source terms
!  The second moment for RRS scattering = 1/20.
!  The adjustment = aircolumn / total optical thickness

        RRS_PHASMOM2 = 0.05D0
        DO B = 1, N_RRSBINS(WI)
          DO N = 1, NLAYERS
            OMEGAMOMS_RRSBIN(N,0,B,WI) = LAYER_AIRCOLUMNS(N) * &
                      BINXSEC(N,WI,B) / DELTAU_VERT(N,WO)
            OMEGAMOMS_RRSBIN(N,1,B,WI) = ZERO
            OMEGAMOMS_RRSBIN(N,2,B,WI) = &
                      OMEGAMOMS_RRSBIN(N,0,B,WI) * RRS_PHASMOM2
          ENDDO
        ENDDO

!  End loop over all inner wavelengths

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_BIN_1

!

      SUBROUTINE LRRS_RAMAN_OPS_BIN_2 &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, MAX_BINS,            & ! Inputs
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,                   & ! Inputs
          NLAYERS, NMOMENTS, NPOINTS_INNER, OFFSET_INNER,           & ! Inputs
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL,                            & ! Inputs
          LAYER_AIRCOLUMNS, DELTAU_VERT, OMEGAMOMS_ELASTIC,         & ! Inputs
          N_RRSBINS, BINXSEC, RRSXSEC_OUT,                          & ! Inputs
          OMEGAMOMS_CABANNES, OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN )   ! Outputs

       IMPLICIT NONE

!  Notes
!  =====

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
!
!  Also required as input for this process are:
!
!      LAYER_AIRCOLUMNS(n) = column air density, [mol/cm^2]
!      RAYLEIGH_XSEC(wo)   = Rayleigh cross-sections on outer grid [cm^2/mol]
!      RAYLEIGH_DEPOL(wo)  = Rayleigh depolarization ratios on outer grib
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
!      OMEGAMOMS_RRSBIN(n,l,b,wi) = SSA x RRS phase moments. L = 0, 1, 2
!                                   Raman Gain terms for each bin.
!      OMEGAMOMS_RRSLOSS(n,l,wi)  = SSA x RRS phase moments. L = 0, 1, 2.
!                                   Raman loss term for the inner point

!  Input arguments
!  ===============

!  Maximum number of spectral points and bins

      INTEGER  , INTENT(IN) :: MAX_POINTS, MAX_BINS

!  Maximum number of layers and input moments

      INTEGER  , INTENT(IN) :: MAX_LAYERS, MAX_MOMENTS

!  Layering

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of elastic moments

      INTEGER  , INTENT(IN) :: NMOMENTS

!  Number of points between the limiting wavelengths

      INTEGER  , INTENT(IN) :: NPOINTS_INNER

!  outer wavelength offset

      INTEGER  , INTENT(IN) :: OFFSET_INNER

!  Method

      LOGICAL, INTENT(IN) :: DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Vertical optical depths

      REAL(FPK), INTENT(IN) :: DELTAU_VERT ( MAX_LAYERS, MAX_POINTS )

!  Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )

!  Input elastic scattering law

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  Bin Mapping numbers

      INTEGER  , INTENT(IN) :: N_RRSBINS  ( MAX_POINTS )

!  Binned Gain term cross-sections

      REAL(FPK), INTENT(IN) :: BINXSEC ( MAX_LAYERS, MAX_POINTS, MAX_BINS )

!  Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT ( MAX_LAYERS, MAX_POINTS )

!  Output arguments
!  ================

!  Scattering quantities

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSBIN  ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Local variables
!  ===============

      INTEGER   :: WI, WO, N, L, B
      REAL(FPK) :: DIF0, DIF2, ADJUST, RRS_PHASMOM2
      REAL(FPK) :: SIGRAY, SIGCAB, RAYMOM, CABMOM, DEPOL

!  Initialize

      OMEGAMOMS_RRSBIN   = zero
      OMEGAMOMS_RRSLOSS  = zero
      OMEGAMOMS_CABANNES = zero

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
!    Essentially adjusting the Rayleigh scattering part to make sure that
!    the elastic component is just that due to Cabannes scattering
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
           SIGCAB = SIGRAY - RRSXSEC_OUT(N,WI)
           ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,WO)
           DIF0   = RRSXSEC_OUT(N,WI) * ADJUST
           DIF2   = ( SIGRAY*RAYMOM - SIGCAB*CABMOM ) * ADJUST
           OMEGAMOMS_CABANNES(N,0,WI) = OMEGAMOMS_ELASTIC(N,0,WO) - DIF0
           OMEGAMOMS_CABANNES(N,1,WI) = OMEGAMOMS_ELASTIC(N,1,WO)
           OMEGAMOMS_CABANNES(N,2,WI) = OMEGAMOMS_ELASTIC(N,2,WO) - DIF2
           DO L = 3, NMOMENTS
            OMEGAMOMS_CABANNES (N,L,WI) = OMEGAMOMS_ELASTIC(N,L,WO)
           ENDDO
          ENDDO

        ENDIF

!  Loss term (Energy balancing mode only)
!  --------------------------------------

!  In the energy balancing mode, need the loss term IOPs-------------

        IF ( DO_ENERGY_BALANCING ) THEN
          RRS_PHASMOM2 = 0.05D0
          DO N = 1, NLAYERS
            ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,WO)
            DIF0   = RRSXSEC_OUT(N,WI) * ADJUST
            OMEGAMOMS_RRSLOSS(N,0,WI) = DIF0
            OMEGAMOMS_RRSLOSS(N,1,WI) = ZERO
            OMEGAMOMS_RRSLOSS(N,2,WI) = DIF0 * RRS_PHASMOM2
          ENDDO
        ENDIF

!  Gain terms (both implementations)
!  ---------------------------------

!  Inelastic scattering: RRS input---------------
!  Calculate optical property input for the Gain RRS source terms
!  The second moment for RRS scattering = 1/20.
!  The adjustment = aircolumn / total optical thickness

        RRS_PHASMOM2 = 0.05D0
        DO B = 1, N_RRSBINS(WI)
          DO N = 1, NLAYERS
            OMEGAMOMS_RRSBIN(N,0,B,WI) = LAYER_AIRCOLUMNS(N) * &
                      BINXSEC(N,WI,B) / DELTAU_VERT(N,WO)
            OMEGAMOMS_RRSBIN(N,1,B,WI) = ZERO
            OMEGAMOMS_RRSBIN(N,2,B,WI) = &
                      OMEGAMOMS_RRSBIN(N,0,B,WI) * RRS_PHASMOM2
          ENDDO
        ENDDO

!  End loop over all inner wavelengths

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_BIN_2

!

      SUBROUTINE LRRS_RAMAN_OPS_MONO_1 &
        ( MAX_LAYERS, MAX_POINTS, MAX_GEOMETRIES,                             & ! Input dimensions
          DO_UPWELLING, DO_DNWELLING, DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, & ! Input flags
          NLAYERS, NGEOMETRIES, NPOINTS_MONO, W_EXCIT,           & ! Input integers
          COSSCAT_UP, COSSCAT_DN, LAYER_AIRCOLUMNS,              & ! Input 
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, DELTAU_VERT,            & ! Inputs Optical
          OMEGAPHASFUNC_ELASTIC_UP, OMEGAPHASFUNC_ELASTIC_DN,    & ! Inputs Optical
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO,                      & ! Inputs Raman Spec.
          OMEGAPHASFUNC_CABANNES_UP, OMEGAPHASFUNC_CABANNES_DN,  & ! Outputs Adjusted Cabannes
          OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSGAIN )                   ! Outputs raman Loss//Gain

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
!
!  Given a set of optical properties defined for all 234 points which
!  contribute to the inelastic scattering for the W_EXCIT wavelength,
!  determine the necessary optical properties for the RRS model over
!  this window which has total number of points NPOINTS_MONO.
!
!  The source terms OMEGAMOMS_RRSGAIN (for light scattered into the
!  excitation wavelength), and OMEGAMOMS_RRSLOSS (light inelastically
!  scattered out), the latter term only considered for the Energy Balanc
!  Method.
!
!  Here is the indexing system for this subroutine
!    n  is the layer index,       n  = 1, .... NLAYERS
!    w  is the wavelength index,  w  = 1, ...  NPOINTS_MONO
!    V  is the geometry index,    v  = 0, .... NGEOMETRIES
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

!  Also required as input for this process are:
!
!      LAYER_AIRCOLUMNS(n) = column air density, [mol/cm^2]
!      RAYLEIGH_XSEC(w)    = Rayleigh cross-sections on complete grid [cm^2/mol]
!      RAYLEIGH_DEPOL(w)   = Rayleigh depolarization ratios on complete grid
!
!  The scattering outputs are:
!
!      OMEGAPHASFUNC_CABANNES_UP(n,v,w) = SSA x Cabannes-adjusted upwelling   phase function
!      OMEGAPHASFUNC_CABANNES_DN(n,v,w) = SSA x Cabannes-adjusted downwelling phase function
!      OMEGAMOMS_RRSGAIN(n,l,w)  = SSA x RRS phase moments. L = 0, 1, 2
!                                   Raman Gain terms for each bin.
!      OMEGAMOMS_RRSLOSS(n,l,1)  = SSA x RRS phase moments. L = 0, 1, 2.
!                                   Raman loss term for the inner point

!  Input arguments
!  ===============

!  Maximum number of spectral points

      INTEGER  , INTENT(IN) :: MAX_POINTS

!  Maximum number of layers and geometries
!   -- Rob mod 5/12/17 for 2p5a, add MAX_GEOMETRIES, remove max_moments

      INTEGER  , INTENT(IN) :: MAX_LAYERS, MAX_GEOMETRIES

!  Directional flags
!   -- Rob mod 5/12/17 for 2p5a, added

      LOGICAL  , INTENT(IN) :: DO_UPWELLING, DO_DNWELLING

!  LRRS Calculation Method

      LOGICAL  , INTENT(IN) :: DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Layering

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of Geometries
!   -- Rob mod 5/12/17 for 2p5a, add NGEOMETRIES, remove NMOMENTS

      INTEGER  , INTENT(IN) :: NGEOMETRIES

!  Number of points nd excitation position

      INTEGER  , INTENT(IN) :: NPOINTS_MONO, W_EXCIT

!  Scattering angles
!   -- Rob mod 5/12/17 for 2p5a, added

      REAL(FPK), INTENT(IN) :: COSSCAT_UP ( MAX_GEOMETRIES )
      REAL(FPK), INTENT(IN) :: COSSCAT_DN ( MAX_GEOMETRIES )

!  Vertical optical depths

      REAL(FPK), INTENT(IN) :: DELTAU_VERT ( MAX_LAYERS, MAX_POINTS )

!  Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )

!  Input elastic scattering phase function products
!   -- Rob mod 5/12/17 for 2p5a, Replaces the phasmoms input (OMEGAMOMS_ELASTIC)

      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  Gain term cross-sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_RANKED ( MAX_LAYERS, MAX_POINTS )

!  Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT_MONO ( MAX_LAYERS )

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

!  Local variables
!  ===============

      INTEGER   :: W, W1, N, V
      REAL(FPK) :: DIF0(MAX_LAYERS), DIF2(MAX_LAYERS), ADJUST, RRS_PHASMOM2, RAMANPHASFUNC
      REAL(FPK) :: SIGRAY, SIGCAB, RAYMOM, CABMOM, DEPOL
      REAL(FPK) :: MU, LEG2_UP(MAX_GEOMETRIES), LEG2_DN(MAX_GEOMETRIES)

!  Initialize
!   -- Rob mod 5/12/17 for 2p5a, replaces the phasmoms output (OMEGAMOMS_CABANNES)

      OMEGAMOMS_RRSGAIN  = zero
      OMEGAMOMS_RRSLOSS  = zero
      OMEGAPHASFUNC_CABANNES_UP = zero
      OMEGAPHASFUNC_CABANNES_DN = zero

!  Legendre calculations
!   -- Rob mod 5/12/17 for 2p5a, need Legendre polynomials for Cabannes adjustment

      IF ( DO_CABANNES_RAMAN ) THEN
         DO V = 1, NGEOMETRIES
           MU = COSSCAT_UP(V) ; LEG2_UP(V) = 1.5_fpk * MU * MU - 0.5_fpk
           MU = COSSCAT_DN(V) ; LEG2_DN(V) = 1.5_fpk * MU * MU - 0.5_fpk
         ENDDO
      ENDIF

      W1 = W_EXCIT

!  Cabannes adjustment
!  -------------------

      IF ( DO_CABANNES_RAMAN ) THEN

!    Essentially adjusting the Rayleigh scattering part to make sure that
!    the elastic component is just that due to Cabannes scattering
!    Only requires the first 0,1,2 moments to be adjusted.

!  We only need to do this for the excitation wavelength.

!  Rayleigh and Cabannes  Legendre(2) moments

        DEPOL  = RAYLEIGH_DEPOL(W1)
        SIGRAY = RAYLEIGH_XSEC (W1)
        CABMOM = 0.25D0 * &
             ( 8.0D0 - 9.0D0*DEPOL ) / ( 4.0D0 - 3.0D0*DEPOL )
        RAYMOM = ( one - DEPOL ) / ( 2.0D0 + DEPOL )

!   -- Rob mod 5/12/17 for 2p5a, Adjust the phase functions for Cabannes scattering

!             DIF0 = RRS Loss-scattering optical depth / total optical depth
!             DIF2 = scatter-weighted RRS Loss / total OD, for moment L = 2

        DO N = 1, NLAYERS
          SIGCAB = SIGRAY - RRSXSEC_OUT_MONO(N)
          ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,W1)
          DIF0(N)   = RRSXSEC_OUT_MONO(N) * ADJUST
          DIF2(N)   = ( SIGRAY*RAYMOM - SIGCAB*CABMOM ) * ADJUST
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
!        DO N = 1, NLAYERS
!          OMEGAMOMS_CABANNES(N,0,1) = OMEGAMOMS_ELASTIC(N,0,W1) - DIF0(N)
!          OMEGAMOMS_CABANNES(N,1,1) = OMEGAMOMS_ELASTIC(N,1,W1)
!          OMEGAMOMS_CABANNES(N,2,1) = OMEGAMOMS_ELASTIC(N,2,W1) - DIF2(N)
!          OMEGAMOMS_CABANNES (N,3:NMOMENTS,1) = OMEGAMOMS_ELASTIC(N,3:NMOMENTS,WO)
!        ENDDO

      ENDIF

!  Loss adjustment
!  ---------------

!  In the energy balancing mode, need the loss terms
!    Only required at the excitation wavelength in the monochromatic case

      IF ( DO_ENERGY_BALANCING ) THEN
        RRS_PHASMOM2 = 0.05D0
        DO N = 1, NLAYERS
          ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,W1)
          DIF0(N)   = RRSXSEC_OUT_MONO(N) * ADJUST
          OMEGAMOMS_RRSLOSS(N,0,1) = DIF0(N)
          OMEGAMOMS_RRSLOSS(N,1,1) = ZERO
          OMEGAMOMS_RRSLOSS(N,2,1) = DIF0(N) * RRS_PHASMOM2
        ENDDO
      ENDIF

!  Gain terms for Inelastic scattering: RRS input
!  ----------------------------------------------

!  Calculate optical property input for the Gain RRS source terms
!  The second moment for RRS scattering = 1/20.
!  The adjustment = aircolumn / total optical thickness
!  The gain term for the excitation wavelength is set to zero

      RRS_PHASMOM2 = 0.05D0
      DO W = 1, NPOINTS_MONO
        DO N = 1, NLAYERS
          OMEGAMOMS_RRSGAIN(N,0,W) = LAYER_AIRCOLUMNS(N) * &
                    RRSXSEC_RANKED(N,W) / DELTAU_VERT(N,W1)
          OMEGAMOMS_RRSGAIN(N,1,W) = ZERO
          OMEGAMOMS_RRSGAIN(N,2,W) = OMEGAMOMS_RRSGAIN(N,0,W) * RRS_PHASMOM2
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_MONO_1

!

      SUBROUTINE LRRS_RAMAN_OPS_MONO_2 &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS,                      & ! Inputs
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,                   & ! Inputs
          NLAYERS, NMOMENTS, NPOINTS_MONO, W_EXCIT,                 & ! Inputs
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL,                            & ! Inputs
          LAYER_AIRCOLUMNS, DELTAU_VERT, OMEGAMOMS_ELASTIC,         & ! Inputs
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO,                         & ! Inputs
          OMEGAMOMS_CABANNES, OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSGAIN )  ! Outputs

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
!
!  Here is the indexing system for this subroutine
!    n  is the layer index,       n  = 1, .... NLAYERS
!    w  is the wavelength index,  w  = 1, ...  NPOINTS_MONO
!    l  is the Leg. moment index, l  = 0, .... NMOMENTS_INPUT
!
!  The control integers for layers and moments are
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
!
!  Also required as input for this process are:
!
!      TEMPERATURES(n)           = Layer Temperature Profile, [K]
!      TEMPERATURES_UNSHIFTED(n) = Unshifted Temperature Profile, [K]
!
!      LAYER_AIRCOLUMNS(n)    = Air column-density Profile, [mol/cm^2]
!      LAYER_AIRCOLUMNS_dT(n) = Air Profile T-derivative,   [mol/cm^2/K]
!
!      RAYLEIGH_XSEC(w)    = Rayleigh cross-sections on complete
!                            grid [cm^2/mol]
!      RAYLEIGH_DEPOL(w)   = Rayleigh depolarization ratios on complete grid
!
!  The scattering outputs are:
!
!      OMEGAMOMS_CABANNES(n,l,1) = SSA x Cabannes-adjusted phasemoments
!                                   (l = 0, 1, ...NMOMENTS_TOTAL)
!      OMEGAMOMS_RRSGAIN(n,l,w)  = SSA x RRS phase moments. L = 0, 1, 2
!                                   Raman Gain terms for each bin.
!      OMEGAMOMS_RRSLOSS(n,l,1)  = SSA x RRS phase moments. L = 0, 1, 2.
!                                   Raman loss term for the inner point

!  Input arguments
!  ===============

!  Maximum number of spectral points

      INTEGER  , INTENT(IN) :: MAX_POINTS

!  Maximum number of layers and input moments

      INTEGER  , INTENT(IN) :: MAX_LAYERS, MAX_MOMENTS

!  Layering

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of elastic moments

      INTEGER  , INTENT(IN) :: NMOMENTS

!  Number of points nd excitation position

      INTEGER  , INTENT(IN) :: NPOINTS_MONO, W_EXCIT

!  Method

      LOGICAL  , INTENT(IN) :: DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )

!  Vertical optical depths

      REAL(FPK), INTENT(IN) :: DELTAU_VERT ( MAX_LAYERS, MAX_POINTS )

!  Input elastic scattering law

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Gain term cross-sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_RANKED ( MAX_LAYERS, MAX_POINTS )

!  Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT_MONO ( MAX_LAYERS )

!  Output arguments
!  ================

!  Scattering quantities

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSGAIN ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Local variables
!  ===============

      INTEGER   :: W, W1, N, L
      REAL(FPK) :: DIF0, DIF2, ADJUST, RRS_PHASMOM2
      REAL(FPK) :: SIGRAY, SIGCAB, RAYMOM, CABMOM, DEPOL

!  Initialize

      OMEGAMOMS_RRSGAIN  = zero
      OMEGAMOMS_RRSLOSS  = zero
      OMEGAMOMS_CABANNES = zero

      W1 = W_EXCIT

!  Cabannes adjustment
!  -------------------

      IF ( DO_CABANNES_RAMAN ) THEN

!    Essentially adjusting the Rayleigh scattering part to make sure that
!    the elastic component is just that due to Cabannes scattering
!    Only requires the first 0,1,2 moments to be adjusted.

!  We only need to do this for the excitation wavelength.

!  Rayleigh and Cabannes  Legendre(2) moments

        DEPOL  = RAYLEIGH_DEPOL(W1)
        SIGRAY = RAYLEIGH_XSEC (W1)
        CABMOM = 0.25D0 * &
             ( 8.0D0 - 9.0D0*DEPOL ) / ( 4.0D0 - 3.0D0*DEPOL )
        RAYMOM = ( one - DEPOL ) / ( 2.0D0 + DEPOL )

!  DIF0 = RRS Loss-scattering optical depth / total optical depth
!  DIF2 = scatter-weighted RRS Loss / total OD, for moment L = 2
!  Cabannes L = 0, adjust OMEGAMOMS by DIF0
!  Cabannes L = 1, no adjustment (just equal to elastic value)
!  Cabannes L = 2, adjust OMEGAMOMS by DIF2
!  Cabannes L = 3,4.. no adjustment (just equal to elastic values)

        DO N = 1, NLAYERS
          SIGCAB = SIGRAY - RRSXSEC_OUT_MONO(N)
          ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,W1)
          DIF0   = RRSXSEC_OUT_MONO(N) * ADJUST
          DIF2   = ( SIGRAY*RAYMOM - SIGCAB*CABMOM ) * ADJUST
          OMEGAMOMS_CABANNES(N,0,1) = OMEGAMOMS_ELASTIC(N,0,W1) - DIF0
          OMEGAMOMS_CABANNES(N,1,1) = OMEGAMOMS_ELASTIC(N,1,W1)
          OMEGAMOMS_CABANNES(N,2,1) = OMEGAMOMS_ELASTIC(N,2,W1) - DIF2
          DO L = 3, NMOMENTS
            OMEGAMOMS_CABANNES (N,L,1) = OMEGAMOMS_ELASTIC(N,L,W1)
          ENDDO
        ENDDO

      ENDIF

!  Loss adjustment
!  ---------------

!  In the energy balancing mode, need the loss terms
!    Only required at the excitation wavelength in the monochromatic case

      IF ( DO_ENERGY_BALANCING ) THEN
        RRS_PHASMOM2 = 0.05D0
        DO N = 1, NLAYERS
          ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,W1)
          DIF0   = RRSXSEC_OUT_MONO(N) * ADJUST
          OMEGAMOMS_RRSLOSS(N,0,1) = DIF0
          OMEGAMOMS_RRSLOSS(N,1,1) = ZERO
          OMEGAMOMS_RRSLOSS(N,2,1) = DIF0 * RRS_PHASMOM2
        ENDDO
      ENDIF

!  Gain terms for Inelastic scattering: RRS input
!  ----------------------------------------------

!  Calculate optical property input for the Gain RRS source terms
!  The second moment for RRS scattering = 1/20.
!  The adjustment = aircolumn / total optical thickness
!  The gain term for the excitation wavelength is set to zero

      RRS_PHASMOM2 = 0.05D0
      DO W = 1, NPOINTS_MONO
        DO N = 1, NLAYERS
          OMEGAMOMS_RRSGAIN(N,0,W) = LAYER_AIRCOLUMNS(N) * &
                    RRSXSEC_RANKED(N,W) / DELTAU_VERT(N,W1)
          OMEGAMOMS_RRSGAIN(N,1,W) = ZERO
          OMEGAMOMS_RRSGAIN(N,2,W) = &
                    OMEGAMOMS_RRSGAIN(N,0,W) * RRS_PHASMOM2
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_MONO_2

!   End module

      END MODULE lrrs_generate_ramanops_m

