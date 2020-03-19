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
! #                                                             #
! #         lrrs_raman_ops_bin       ( no  linearization)       #
! #         lrrs_raman_ops_bin_plus  (with linearization)       #
! #                                                             #
! #         lrrs_raman_ops_mono       ( no  linearization)      #
! #         lrrs_raman_ops_mono_plus  (with linearization)      #
! #                                                             #
! #  These subroutines calculate the Raman optical property     #
! #  inputs and their linearizations. The calculation is based  #
! #  on the elastic IOP inputs which are assumed known, and the #
! #  routines adjust these inputs and create new Raman-only     #
! #  IOP inputs by using the Raman cross-sections and the Air   #
! #  column density in each layer.                              #
! #                                                             #
! #  Monochromatic Routines added 8 February 2011               #
! #                                                             #
! ###############################################################


      MODULE lrrs_generate_ramanops

      USE LRRS_PARS

      PRIVATE
      PUBLIC :: LRRS_RAMAN_OPS_BIN,&
                LRRS_RAMAN_OPS_BIN_PLUS,&
                LRRS_RAMAN_OPS_MONO,&
                LRRS_RAMAN_OPS_MONO_PLUS

      CONTAINS

      SUBROUTINE LRRS_RAMAN_OPS_BIN &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, MAX_BINS, &
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, &
          NLAYERS, NMOMENTS, NPOINTS_INNER, OFFSET_INNER, &
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, &
          LAYER_AIRCOLUMNS, DELTAU_VERT, OMEGAMOMS_ELASTIC, &
          N_RRSBINS, BINXSEC, RRSXSEC_OUT, &
          OMEGAMOMS_CABANNES, OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN )

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
!  The control INTEGER ::s for layers and moments are
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
!      LAYER_AIRCOLUMNS(n) = column air density, [mol/cm^2]
!      RAYLEIGH_XSEC(wo)   = Rayleigh cross-sections on outer grid [cm^2
!      RAYLEIGH_DEPOL(wo)  = Rayleigh depolarization ratios on outer gri
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

      INTEGER, INTENT(IN) ::          MAX_POINTS, MAX_BINS

!  Maximum number of layers and input moments

      INTEGER, INTENT(IN) ::          MAX_LAYERS, MAX_MOMENTS

!  Layering

      INTEGER, INTENT(IN) ::          NLAYERS

!  Number of elastic moments

      INTEGER, INTENT(IN) ::          NMOMENTS

!  Number of points between the limiting wavelengths

      INTEGER, INTENT(IN) ::          NPOINTS_INNER

!  outer wavelength offset

      INTEGER, INTENT(IN) ::          OFFSET_INNER

!  Method

      LOGICAL, INTENT(IN) ::          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Vertical optical depths

      REAL(FPK), INTENT(IN) :: DELTAU_VERT ( MAX_LAYERS, MAX_POINTS )

!  Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )

!  Input elastic scattering law

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  Bin Mapping numbers

      INTEGER, INTENT(IN) ::          N_RRSBINS  ( MAX_POINTS )

!  Binned Gain term cross-sections

      REAL(FPK), INTENT(IN) :: BINXSEC ( MAX_LAYERS, MAX_POINTS, MAX_BINS )

!  Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT ( MAX_LAYERS, MAX_POINTS )

!  Output arguments
!  ================

!  Scattering quantities

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS &
               ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSBIN &
               ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Local variables
!  ===============

      INTEGER ::   WI, WO, N, L, B
      REAL(FPK) :: DIF0, DIF2, ADJUST, RRS_PHASMOM2
      REAL(FPK) :: SIGRAY, SIGCAB, RAYMOM, CABMOM, DEPOL

!  Initialize

      OMEGAMOMS_RRSBIN   = 0.0d0
      OMEGAMOMS_RRSLOSS  = 0.0d0
      OMEGAMOMS_CABANNES = 0.0d0

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
          RAYMOM = ( 1.0D0 - DEPOL ) / ( 2.0D0 + DEPOL )

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
!if ( n.eq.27.and.wi.eq.1) write(*,*)WI,OMEGAMOMS_CABANNES(N,0,WI)
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
            OMEGAMOMS_RRSLOSS(N,1,WI) = 0.0d0
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
            OMEGAMOMS_RRSBIN(N,1,B,WI) = 0.0D0
            OMEGAMOMS_RRSBIN(N,2,B,WI) = &
                   OMEGAMOMS_RRSBIN(N,0,B,WI) * RRS_PHASMOM2
          ENDDO
        ENDDO

!  End loop over all inner wavelengths

      enddo

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_BIN

!  Linearized version

      SUBROUTINE LRRS_RAMAN_OPS_BIN_PLUS &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, MAX_BINS, MAX_VARS, &
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, DO_NORMALIZED_WFS, &
          DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF, &
          NLAYERS, NMOMENTS, NPOINTS_INNER, OFFSET_INNER, NVARS, &
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, &
          TEMPERATURES, TEMPERATURES_UNSHIFTED, &
          LAYER_AIRCOLUMNS, LAYER_AIRCOLUMNS_DT, &
          DELTAU_VERT,     OMEGAMOMS_ELASTIC, &
          L_DELTAU_VERT, L_OMEGAMOMS_ELASTIC, &
          N_RRSBINS, BINXSEC, RRSXSEC_OUT, &
          BINXSEC_DT, RRSXSEC_OUT_DT, &
          OMEGAMOMS_CABANNES, L_OMEGAMOMS_CABANNES, &
          OMEGAMOMS_RRSLOSS,  L_OMEGAMOMS_RRSLOSS, &
          OMEGAMOMS_RRSBIN,   L_OMEGAMOMS_RRSBIN )

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
!  The control INTEGER ::s for layers and moments are
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

      INTEGER, INTENT(IN) ::          MAX_POINTS, MAX_BINS

!  Maximum number of layers and input moments

      INTEGER, INTENT(IN) ::          MAX_LAYERS, MAX_MOMENTS

!  Maximum number of weighting functions

      INTEGER, INTENT(IN) ::          MAX_VARS

!  Flags
!  -----

!  RRS Method

      LOGICAL, INTENT(IN) ::          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Temperature and Air Profile Linearization flags

      LOGICAL, INTENT(IN) ::          DO_TEMPPROFILE_WFS
      LOGICAL, INTENT(IN) ::          DO_AIRPROFILE_WFS

!  Temperature shift Linearization flag (Column WF only)

      LOGICAL, INTENT(IN) ::          DO_TEMPSHIFT_WF

!  use of normalized weighting functions (Profile WFs only)
!   T-shift Jacobian is automatically normalized

      LOGICAL, INTENT(IN) ::          DO_NORMALIZED_WFS

!  Numbers
!  -------

!  Layering

      INTEGER, INTENT(IN) ::          NLAYERS

!  Number of elastic moments

      INTEGER, INTENT(IN) ::          NMOMENTS

!  Number of points between the limiting wavelengths

      INTEGER, INTENT(IN) ::          NPOINTS_INNER

!  outer wavelength offset

      INTEGER, INTENT(IN) ::          OFFSET_INNER

!  Number of weighting functions

      INTEGER, INTENT(IN) ::          NVARS ( MAX_LAYERS )

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

      INTEGER, INTENT(IN) ::          N_RRSBINS  ( MAX_POINTS )

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

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS &
               ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSBIN &
               ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Linearizations

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_CABANNES &
           ( MAX_VARS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSLOSS &
               ( MAX_VARS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSBIN &
               ( MAX_VARS, MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Local variables
!  ===============

      INTEGER ::          WI, WO, N, L, B, Q, NQ
      REAL(FPK) :: DIF0, DIF2, ADJUST, NORM, NORM2
      REAL(FPK) :: L_DIF0, L_DIF2, L_ADJUST, RATIO_1
      REAL(FPK) :: SIGRAY, DEPOL, LOSSM0, LOSSM2
      REAL(FPK) :: TSHIFT, LOSSM0_dT, LOSSM2_dT
      REAL(FPK) :: RAYMOM, CABMOM, RRS_PHASMOM2, AIRCOL
      REAL(FPK) :: TERM_1, TERM_2, TERM_3

!  Initialize

      OMEGAMOMS_RRSBIN   = 0.0d0
      OMEGAMOMS_RRSLOSS  = 0.0d0
      OMEGAMOMS_CABANNES = 0.0d0
      L_OMEGAMOMS_RRSBIN   = 0.0d0
      L_OMEGAMOMS_RRSLOSS  = 0.0d0
      L_OMEGAMOMS_CABANNES = 0.0d0

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
          RAYMOM = ( 1.0D0 - DEPOL ) / ( 2.0D0 + DEPOL )

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
             NORM = 1.0d0
             IF ( DO_NORMALIZED_WFS) NORM = AIRCOL
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

!           if ( n.eq.27.and.wi.eq.1) write(*,*)WI,1,&
!              OMEGAMOMS_CABANNES(N,0,WI),L_OMEGAMOMS_CABANNES(1,N,0,WI)
!           if ( n.eq.27.and.wi.eq.1) write(*,*)WI,2,&
!              OMEGAMOMS_CABANNES(N,0,WI),L_OMEGAMOMS_CABANNES(2,N,0,WI)

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
           OMEGAMOMS_RRSLOSS(N,1,WI) = 0.0d0
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
             NORM = 1.0d0
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
              L_OMEGAMOMS_RRSLOSS(Q,N,1,WI) = 0.0d0
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
           OMEGAMOMS_RRSBIN(N,1,B,WI) = 0.0D0
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
            NORM = 1.0d0
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
              L_OMEGAMOMS_RRSBIN(Q,N,1,B,WI) = 0.0D0
            ENDDO
          ENDDO

!  End layer loop

        ENDDO

!  End loop over all inner wavelengths

      enddo

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_BIN_PLUS


      SUBROUTINE LRRS_RAMAN_OPS_MONO &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, &
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, &
          NLAYERS, NMOMENTS_INPUT, NPOINTS_MONO, W_EXCIT, &
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, &
          LAYER_AIRCOLUMNS, DELTAU_VERT, OMEGAMOMS_ELASTIC, &
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO, &
          OMEGAMOMS_CABANNES, OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSGAIN )

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

!  The control INTEGER ::s for layers and moments are
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
!      OMEGAMOMS_RRSGAIN(n,l,w)  = SSA x RRS phasemoments. L = 0, 1, 2
!                                   Raman Gain terms for each bin.
!      OMEGAMOMS_RRSLOSS(n,l,1)  = SSA x RRS phasemoments. L = 0, 1, 2.
!                                   Raman loss term for the iner point

!  Input arguments
!  ===============

!  Maximum number of spectral points

      INTEGER, INTENT(IN) ::          MAX_POINTS

!  Maximum number of layers and input moments

      INTEGER, INTENT(IN) ::          MAX_LAYERS, MAX_MOMENTS

!  Layering

      INTEGER, INTENT(IN) ::          NLAYERS

!  Number of elastic moments

      INTEGER, INTENT(IN) ::          NMOMENTS_INPUT

!  Number of points nd excitation position

      INTEGER, INTENT(IN) ::          NPOINTS_MONO, W_EXCIT

!  Method

      LOGICAL, INTENT(IN) ::          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  Layer air columns

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )

!  Vertical optical depths

      REAL(FPK), INTENT(IN) :: DELTAU_VERT ( MAX_LAYERS, MAX_POINTS )

!  Input elastic scattering law

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Gain term cross-sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_RANKED ( MAX_LAYERS, MAX_POINTS )

!  Loss term cross sections

      REAL(FPK), INTENT(IN) :: RRSXSEC_OUT_MONO ( MAX_LAYERS )

!  Output arguments
!  ================

!  Scattering quantities

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS &
               ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSGAIN &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Local variables
!  ===============

      INTEGER ::          W, W1, N, L
      REAL(FPK) :: DIF0, DIF2, ADJUST, RRS_PHASMOM2
      REAL(FPK) :: SIGRAY, SIGCAB, RAYMOM, CABMOM, DEPOL

!  Initialize

      OMEGAMOMS_RRSGAIN  = 0.0d0
      OMEGAMOMS_RRSLOSS  = 0.0d0
      OMEGAMOMS_CABANNES = 0.0d0

      W1 = W_EXCIT

!  Cabannes adjustment
!  -------------------

      IF ( DO_CABANNES_RAMAN ) THEN

!    Essentially adjusting the Rayleigh scattering part to make sure tha
!    the elastic component is just the due to Cabannes scattering
!    Only requires the first 0,1,2 moments to be adjusted.

!  We only need to do this for the excitation wavelength.

!  Rayleigh and Cabannes  Legendre(2) moments

        DEPOL  = RAYLEIGH_DEPOL(W1)
        SIGRAY = RAYLEIGH_XSEC (W1)
        CABMOM = 0.25D0 * &
             ( 8.0D0 - 9.0D0*DEPOL ) / ( 4.0D0 - 3.0D0*DEPOL )
        RAYMOM = ( 1.0D0 - DEPOL ) / ( 2.0D0 + DEPOL )

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
          DO L = 3, NMOMENTS_INPUT
            OMEGAMOMS_CABANNES (N,L,1) = OMEGAMOMS_ELASTIC(N,L,W1)
          ENDDO
        ENDDO

      ENDIF

!  Loss adjustment
!  ---------------

!  In the energy balancing mode, need the loss terms
!    Only required at the excitation wavelength in the monochromatic cas

      IF ( DO_ENERGY_BALANCING ) THEN
        RRS_PHASMOM2 = 0.05D0
        DO N = 1, NLAYERS
          ADJUST = LAYER_AIRCOLUMNS(N) / DELTAU_VERT(N,W1)
          DIF0   = RRSXSEC_OUT_MONO(N) * ADJUST
          OMEGAMOMS_RRSLOSS(N,0,1) = DIF0
          OMEGAMOMS_RRSLOSS(N,1,1) = 0.0d0
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
          OMEGAMOMS_RRSGAIN(N,1,W) = 0.0D0
          OMEGAMOMS_RRSGAIN(N,2,W) = &
                   OMEGAMOMS_RRSGAIN(N,0,W) * RRS_PHASMOM2
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_MONO

!  Linearized code

      SUBROUTINE LRRS_RAMAN_OPS_MONO_PLUS &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, MAX_VARS, &
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, DO_NORMALIZED_WFS, &
          DO_AIRPROFILE_WFS, DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF, &
          NLAYERS, NMOMENTS_INPUT, NPOINTS_MONO, W_EXCIT, NVARS, &
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, &
          TEMPERATURES, TEMPERATURES_UNSHIFTED, &
          LAYER_AIRCOLUMNS, LAYER_AIRCOLUMNS_DT, &
          DELTAU_VERT,     OMEGAMOMS_ELASTIC, &
          L_DELTAU_VERT, L_OMEGAMOMS_ELASTIC, &
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO, &
          RRSXSEC_RANKED_DT, RRSXSEC_OUT_MONO_DT, &
          OMEGAMOMS_CABANNES, L_OMEGAMOMS_CABANNES, &
          OMEGAMOMS_RRSLOSS,  L_OMEGAMOMS_RRSLOSS, &
          OMEGAMOMS_RRSGAIN,  L_OMEGAMOMS_RRSGAIN )

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

!  The control INTEGER ::s for layers and moments are
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

      INTEGER, INTENT(IN) ::          MAX_POINTS

!  Maximum number of layers, input moments, weighting functions

      INTEGER, INTENT(IN) ::          MAX_LAYERS, MAX_MOMENTS, MAX_VARS

!  flags
!  -----

!  Method

      LOGICAL, INTENT(IN) ::          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN

!  Temperature and Air Profile Linearization flags

      LOGICAL, INTENT(IN) ::          DO_TEMPPROFILE_WFS
      LOGICAL, INTENT(IN) ::          DO_AIRPROFILE_WFS

!  Temperature shift Linearization flag (Column WF only)

      LOGICAL, INTENT(IN) ::          DO_TEMPSHIFT_WF

!  use of normalized weighting functions (Profile WFs only)
!   T-shift Jacobian is automatically normalized

      LOGICAL, INTENT(IN) ::          DO_NORMALIZED_WFS

!  Numbers
!  -------

!  Layering

      INTEGER, INTENT(IN) ::          NLAYERS

!  Number of elastic moments

      INTEGER, INTENT(IN) ::          NMOMENTS_INPUT

!  Number of points nd excitation position

      INTEGER, INTENT(IN) ::          NPOINTS_MONO, W_EXCIT

!  Number of weighting functions

      INTEGER, INTENT(IN) ::          NVARS ( MAX_LAYERS )

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

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Linearization of optical properties

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT &
                ( MAX_VARS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_ELASTIC &
           ( MAX_VARS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

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

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Loss and gain terms

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS &
               ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSGAIN &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Linearizations

      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_CABANNES &
               ( MAX_VARS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSLOSS &
               ( MAX_VARS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_RRSGAIN &
               ( MAX_VARS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Local variables
!  ===============

      INTEGER ::          W1, W, N, L, Q, NQ
      REAL(FPK) :: DIF0, DIF2, ADJUST, NORM, NORM2
      REAL(FPK) :: L_DIF0, L_DIF2, L_ADJUST, RATIO_1
      REAL(FPK) :: SIGRAY, DEPOL, LOSSM0, LOSSM2
      REAL(FPK) :: TSHIFT, LOSSM0_dT, LOSSM2_dT
      REAL(FPK) :: RAYMOM, CABMOM, RRS_PHASMOM2, AIRCOL
      REAL(FPK) :: TERM_1, TERM_2, TERM_3

!  Initialize

      OMEGAMOMS_RRSGAIN  = 0.0d0
      OMEGAMOMS_RRSLOSS  = 0.0d0
      OMEGAMOMS_CABANNES = 0.0d0
      L_OMEGAMOMS_RRSGAIN  = 0.0d0
      L_OMEGAMOMS_RRSLOSS  = 0.0d0
      L_OMEGAMOMS_CABANNES = 0.0d0

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
         RAYMOM = ( 1.0D0 - DEPOL ) / ( 2.0D0 + DEPOL )

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
             NORM = 1.0d0
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
           OMEGAMOMS_RRSLOSS(N,1,1) = 0.0d0
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
             NORM = 1.0d0
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
             LOSSM2_dT = LOSSM0_dT * CABMOM
             TERM_1 = TSHIFT * LAYER_AIRCOLUMNS(N)    * LOSSM0_dT
             TERM_2 = TSHIFT * LAYER_AIRCOLUMNS_dT(N) * LOSSM0
             TERM_3 = L_DELTAU_VERT(Q,N,W1)  * DIF0
             L_DIF0 = (TERM_1 + TERM_2 - TERM_3 ) / DELTAU_VERT(N,W1)
             L_OMEGAMOMS_RRSLOSS(Q,N,0,1) = L_DIF0
             L_OMEGAMOMS_RRSLOSS(Q,N,2,1) = L_DIF0 * RRS_PHASMOM2
           ENDIF

!  First moment linearization is zero, always

           DO Q = 1, NVARS(N)
              L_OMEGAMOMS_RRSLOSS(Q,N,1,1) = 0.0d0
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
          OMEGAMOMS_RRSGAIN(N,1,W) = 0.0D0
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
            NORM = 1.0d0
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
            L_OMEGAMOMS_RRSGAIN(Q,N,1,W) = 0.0D0
          ENDDO

!  End layer loop

        ENDDO

!  End loop over points

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_RAMAN_OPS_MONO_PLUS


      END MODULE lrrs_generate_ramanops

