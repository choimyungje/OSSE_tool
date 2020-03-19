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
! #             SETUP_MASTER_BIN  (master)                      #
! #             SETUP_MASTER_MONO (master)                      #
! #                                                             #
! #  The BIN module generates complete LRRS set-up, 3 tasks     #
! #     1. Delta-M scaling of Elastic OPs                       #
! #     2. Raman Cross-sections (Spectroscopy)                  #
! #     3. Raman Optical Properties                             #
! #                                                             #
! #  The MONO module finishes complete LRRS set-up, 3 tasks     #
! #     1. Delta-M scaling of Elastic OPs                       #
! #     2. Raman Optical Properties                             #
! #         [Raman Cross-sections (Spectroscopy) already done   #
! #                                                             #
! #    Monochromatic routine added 8 february 2011              #
! #                                                             #
! ###############################################################


      MODULE lrrs_setup_master

      PRIVATE
      PUBLIC :: SETUP_MASTER_BIN,&
                SETUP_MASTER_MONO

      CONTAINS

      SUBROUTINE SETUP_MASTER_BIN &
        ( DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, DO_DELTAM_SCALING, &
          NLAYERS, NSTREAMS, NMOMENTS_INPUT, &
          NPOINTS_INNER, OFFSET_INNER, NPOINTS_OUTER, &
          LAYER_TEMPERATURES, LAYER_AIRCOLUMNS, &
          LAMBDAS_RANKED, FLUXES_RANKED, BINLOWER, BINUPPER, &
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, &
          DELTAU_INPUT_UNSCALED, OMEGAMOMS_ELASTIC_UNSCALED, &
            NMOMENTS, TRUNC_FACTORS, &
            DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, &
            N_RRSBINS, BINMAP, RINGSPEC_1, &
            OMEGAMOMS_CABANNES_UNSCALED, OMEGAMOMS_CABANNES, &
            OMEGAMOMS_RRSLOSS_UNSCALED,  OMEGAMOMS_RRSLOSS, &
            OMEGAMOMS_RRSBIN_UNSCALED,   OMEGAMOMS_RRSBIN, &
            FAIL, MESSAGE, ACTION )

!  LRRS include files
!  ==================

      USE LRRS_PARS
      USE LRRS_DELTAMSCALING
      USE LRRS_RAMAN_SPECTROSCOPY
      USE LRRS_GENERATE_RAMANOPS

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Flags
!  -----

!  deltam scaling flag (introduced July 24th, 2006)

      LOGICAL, INTENT(IN) ::          DO_DELTAM_SCALING

!  Solution method (options are mutually exclusive)

      LOGICAL, INTENT(IN) ::          DO_CABANNES_RAMAN
      LOGICAL, INTENT(IN) ::          DO_ENERGY_BALANCING

!  Scattering discretization
!  -------------------------

!  number of streams

      INTEGER, INTENT(IN) ::          NSTREAMS

!  Atmospheric quantities
!  ----------------------

!  number of layers

      INTEGER, INTENT(IN) ::          NLAYERS

!  Input layer temperatures, must be in deg K

      REAL(FPK), INTENT(IN) :: LAYER_TEMPERATURES ( MAX_LAYERS )

!  Input layer Air columns, should be in mol/cm^2 or [DU]

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS ( MAX_LAYERS )

!  Wavelengths/Fluxes are defined on  outer grid for the binning realiza

      REAL(FPK), INTENT(IN) :: LAMBDAS_RANKED ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: FLUXES_RANKED  ( MAX_POINTS )

!  Basic Rayleigh data
!     Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  For the binning realization
!  ---------------------------

!  outer/inner wavelength range, and offset for the inner range
!    Depends on the choice of solar spectrum

      INTEGER, INTENT(IN) ::          OFFSET_INNER
      INTEGER, INTENT(IN) ::          NPOINTS_INNER
      INTEGER, INTENT(IN) ::          NPOINTS_OUTER

!  Bin upper and lower limits

      REAL(FPK), INTENT(IN) :: BINLOWER ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BINUPPER ( MAX_POINTS )

!  Unscaled elastic scattering properties
!  ---------------------------------------

!  These must be defined on the outer wavelength grid (binning)

!  Number of input phase function Legendre moments
!    THIS IS A BASIC INPUT THAT USER MUST PROVIDE

      INTEGER, INTENT(IN) ::          NMOMENTS_INPUT

!  Unscaled quantities, Elastic input

      REAL(FPK), INTENT(IN) :: DELTAU_INPUT_UNSCALED ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

!  Output arguments
!  ================

!  Scaled elastic-scattering optical properties
!  --------------------------------------------

!    These must be defined on the outer wavelength grid (binning)

!  Number of moments ( = 2 * NSTREAMS )

      INTEGER, INTENT(OUT) ::          NMOMENTS

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(OUT) :: DELTAU_VERT_INPUT   ( MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_ELASTIC &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Truncation factors (zero if no deltam scaling)

      REAL(FPK), INTENT(OUT) :: TRUNC_FACTORS         ( MAX_LAYERS, MAX_POINTS )

!  RRS spectroscopic quantities
!  ----------------------------

!  Bin Mapping quantities (internal derivation)

      INTEGER, INTENT(OUT) ::          N_RRSBINS  ( MAX_POINTS )
      INTEGER, INTENT(OUT) ::          BINMAP ( MAX_BINS, MAX_POINTS )

!  Chance/Spurr Ring Spectrum. DIAGNOSTIC ONLY.

      REAL(FPK), INTENT(OUT) :: RINGSPEC_1 ( MAX_POINTS )

!  Unscaled Raman quantities
!  -------------------------

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSBIN_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Scaled Raman-scattering optical properties
!  ------------------------------------------

!  Derived Cabannes scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived rotational raman scattering (Internal variable). Loss term
!    Single scatter albedo * phase function moments
!    Only required for the Energy-balance approximation.

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived rotational Raman scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSBIN &
               ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@
!  Exception handling for MAX_BINS Dimension
      LOGICAL, INTENT(OUT)   :: FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE, ACTION
!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@

!  Local Arrays
!  ============

!  shift dimensions

      INTEGER, PARAMETER :: MAXN2 = 48, MAXO2 = 185

!  Positions

      REAL(FPK) :: O2POS (MAXO2)
      REAL(FPK) :: N2POS (MAXN2)

!  Gammas

      REAL(FPK) :: GAMMA_O2(3)
      REAL(FPK) :: GAMMA_N2(3)

!  Basic coefficients for O2 and N2

      REAL(FPK) :: RRS_N2REF ( MAX_LAYERS, MAXN2 )
      REAL(FPK) :: RRS_O2REF ( MAX_LAYERS, MAXO2 )

!  T-derivatives of Basic coefficients for O2 and N2
!   ******** NOT REQUIRED HERE

      REAL(FPK) :: LINRRS_N2REF ( MAX_LAYERS, MAXN2 )
      REAL(FPK) :: LINRRS_O2REF ( MAX_LAYERS, MAXO2 )

!  Loss term cross-section and number of shifts = N_Stokes + N_Antistoke

      REAL(FPK) :: RRSXSEC_OUT ( MAX_LAYERS, MAX_POINTS )

!  Binned Gain term cross-sections

      REAL(FPK) :: BINXSEC ( MAX_LAYERS, MAX_POINTS, MAX_BINS )

!  Local flag. No derivatives here.

      LOGICAL, PARAMETER :: DO_WF = .false.

!  debug

!      INTEGER ::    n, l, cpt

!  DELTA-M Scaling
!  ===============

!  REMARK: In Version 2.2, this is automatic, controlled only
!          by the DO_DELTAM_SCALING flag.

      CALL LRRS_DELTAM_SCALING_2P2 &
         ( MAX_LAYERS, MAX_MOMENTS_INPUT, MAX_MOMENTS, &
           MAX_POINTS, DO_DELTAM_SCALING, &
           NPOINTS_OUTER, NLAYERS, NSTREAMS, NMOMENTS_INPUT, &
           DELTAU_INPUT_UNSCALED, OMEGAMOMS_ELASTIC_UNSCALED, &
           NMOMENTS, TRUNC_FACTORS, &
           DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC )

!      write(*,*)DELTAU_VERT_INPUT(27,1)
!      write(*,*)OMEGAMOMS_ELASTIC_UNSCALED(26,5,1)
!  debug
!      write(*,*)nmoments_input, nmoments
!      do n = 30, 30
!        write(*,*)n,deltau_input_unscaled(n,1),deltau_vert_input(n,1)
!        do l = 0, 16
!          write(*,*)OMEGAMOMS_ELASTIC_UNSCALED(n,l,1),
!     &      OMEGAMOMS_ELASTIC(n,l,1)
!        enddo
!      enddo
!      pause'after dmscaling'

!  debug
!      DO CPT = 1, NPOINTS_OUTER
!        DO N = 1, NLAYERS
!      write(86,*)cpt,n, DELTAU_VERT_INPUT(n,cpt), &
!           OMEGAMOMS_ELASTIC(N,0,CPT),OMEGAMOMS_ELASTIC(N,1,CPT), &
!           OMEGAMOMS_ELASTIC(N,2,CPT),OMEGAMOMS_ELASTIC(N,5,CPT)
!          ENDDO
!        ENDDO
!        pause

!  Generation of Raman-scattering Optical Properties (ROPS)
!  ========================================================

!  REMARK: THIS is an Automatic procedure in VERSION 2.2

!  THERE ARE 2 MAIN STEPS
!     (1) SPECTROSCOPY - generate the Raman cross-sections, optional T-d
!     (2) MAKE ROPS    - calculate the additional OPs

!  STEP 1. Raman Cross-sections generator
!  --------------------------------------

!  Get the Raman spectroscopy

      CALL RRS_BASIC_SPECTROSCOPY &
         ( MAX_LAYERS, NLAYERS, DO_WF, LAYER_TEMPERATURES, MAXN2, MAXO2, &
           N2POS, O2POS, GAMMA_N2, GAMMA_O2, &
           RRS_N2REF, RRS_O2REF, LINRRS_N2REF, LINRRS_O2REF )

!  Given outer and inner windows and wavelengths + bins and solar fluxes
!  this routine generates the following Raman cross-sections:

!   Each wavelength w,                 the number of bins (N_RRSBINS)
!   Each wavelength w, bin b,          the binning map BINMAP(w,b)
!   Each wavelength w, bin b, layer n, Raman gain X-section BINXSEC(w,b,
!   Each wavelength w, layer n,        Raman loss X-section RRSXSEC_OUT(

!  RINGSPEC_1 is the Fraunhofer "Ring spectrum" for this input field.
!    --This is the standard SAO calculation, output as a diagnostic here

      call LRRS_RAMAN_XSECS_BIN &
        ( MAX_LAYERS, MAX_POINTS, MAX_BINS, &
          NLAYERS, NPOINTS_INNER, NPOINTS_OUTER, OFFSET_INNER, &
          BINLOWER, BINUPPER, LAMBDAS_RANKED, FLUXES_RANKED, &
          MAXN2, MAXO2, RRS_N2REF, RRS_O2REF, &
          N2POS, O2POS, GAMMA_N2, GAMMA_O2, &
          N_RRSBINS, BINMAP, BINXSEC, RRSXSEC_OUT, RINGSPEC_1, &
          FAIL, MESSAGE, ACTION )

!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@
!  Exception handling for MAX_BINS Dimension. Exit immediately
          IF ( FAIL) RETURN
!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@

!  STEP 2. Make the Raman Optical Depths
!  -------------------------------------

!  Formerly, this was the responsibility of the User
!  Now taken care of automatically with these two routines:

!          lrrs_raman_ops_bin_plus (with linearizations);
!          lrrs_raman_ops_bin      (no linearization)

!  These routines are completely stand-alone.  In order to make them work
!  the user must supply the following:

!      Dimensioning inputs
!      LRRS method (Energy-balancing vs. Cabannes). First is usual.
!      Number of layers, Inner window points and offset, # of bins
!      Number of phase function expansion coefficients (Elastic)
!      Rayleigh X-secs(**) + depolarizations, Layer Air columns(**)
!      Elastic optical properties as defined above
!      Raman Loss (RRSXSEC_OUT) and binned gain (BINXSEC) cross-sections

!  (**) Note - cross sections must be [cm^2/mol], Air columns [mol/cm^2]

      call lrrs_raman_ops_bin &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS_INPUT, MAX_BINS, &
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, &
          NLAYERS, NMOMENTS_INPUT, NPOINTS_INNER, OFFSET_INNER, &
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, LAYER_AIRCOLUMNS, &
          DELTAU_INPUT_UNSCALED, OMEGAMOMS_ELASTIC_UNSCALED, &
          N_RRSBINS, BINXSEC, RRSXSEC_OUT, &
          OMEGAMOMS_CABANNES_UNSCALED, &
          OMEGAMOMS_RRSLOSS_UNSCALED, &
          OMEGAMOMS_RRSBIN_UNSCALED )

!  Debug
!      write(*,*)n_RRSBINS(1)
!      do n = 1, nlayers
!        write(65,'(i4,1p3e22.12)')n,OMEGAMOMS_RRSBIN_UNSCALED(n,0,24,1)
!      enddo
!      pause

!  Call the Raman IOP generator, scaled

      call lrrs_raman_ops_bin &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, MAX_BINS, &
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, &
          NLAYERS, NMOMENTS, NPOINTS_INNER, OFFSET_INNER, &
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, LAYER_AIRCOLUMNS, &
          DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, &
          N_RRSBINS, BINXSEC, RRSXSEC_OUT, &
          OMEGAMOMS_CABANNES, OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN )

!  debug
!        DO N = 1, NLAYERS
!          DO L = 0, NMOMENTS
!            DO CPT = 1, NPOINTS_OUTER
!              write(87,*)n,l,cpt,OMEGAMOMS_ELASTIC(N,L,CPT)
!            ENDDO
!          ENDDO
!        ENDDO
!        pause

!  Finish

      return
      END SUBROUTINE SETUP_MASTER_BIN

!

      SUBROUTINE SETUP_MASTER_MONO &
        ( DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, DO_DELTAM_SCALING, &
          NLAYERS, NSTREAMS, NMOMENTS_INPUT, &
          LAMBDA_EXCIT, NPOINTS_MONO, W_EXCIT, &
          LAYER_TEMPERATURES, LAYER_AIRCOLUMNS, &
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, &
          DELTAU_INPUT_UNSCALED, OMEGAMOMS_ELASTIC_UNSCALED, &
            NMOMENTS, TRUNC_FACTORS, &
            DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, &
            OMEGAMOMS_CABANNES_UNSCALED, OMEGAMOMS_CABANNES, &
            OMEGAMOMS_RRSLOSS_UNSCALED,  OMEGAMOMS_RRSLOSS, &
            OMEGAMOMS_RRSGAIN_UNSCALED,  OMEGAMOMS_RRSGAIN )

!  LRRS include files
!  ==================

      USE LRRS_PARS
      USE LRRS_DELTAMSCALING
      USE LRRS_RAMAN_SPECTROSCOPY
      USE LRRS_GENERATE_RAMANOPS

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Flags
!  -----

!  deltam scaling flag (introduced July 24th, 2006)

      LOGICAL, INTENT(IN) ::          DO_DELTAM_SCALING

!  Solution method (options are mutually exclusive)

      LOGICAL, INTENT(IN) ::          DO_CABANNES_RAMAN
      LOGICAL, INTENT(IN) ::          DO_ENERGY_BALANCING

!  Scattering discretization
!  -------------------------

!  number of streams

      INTEGER, INTENT(IN) ::          NSTREAMS

!  Atmospheric quantities
!  ----------------------

!  number of layers

      INTEGER, INTENT(IN) ::          NLAYERS

!  Input layer temperatures, must be in deg K

      REAL(FPK), INTENT(IN) :: LAYER_TEMPERATURES ( MAX_LAYERS )

!  Input layer Air columns, should be in mol/cm^2 or [DU]

      REAL(FPK), INTENT(IN) :: LAYER_AIRCOLUMNS ( MAX_LAYERS )

!  MONOCHROMATIC Wavelengths
!  -------------------------

!  Excitation wavelength

      REAL(FPK), INTENT(IN) :: LAMBDA_EXCIT

!  Ranked wavelengths and NPOINTS_MONO = 234 points
!  W_EXCIT = position of excitation wavelength in the ranked set
!    NOTE: These are inputs, and should be calculated in the Iopsetup
!          routine RRS_LAMDASRANKED_MONO.

      INTEGER, INTENT(IN) ::          NPOINTS_MONO, W_EXCIT
!      REAL(FPK), INTENT(IN) :: LAMBDAS_RANKED  ( MAX_POINTS )

!  Basic Rayleigh data
!     Rayleigh Cross-sections and depolarization ratios

      REAL(FPK), INTENT(IN) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  Unscaled elastic scattering properties
!  ---------------------------------------

!  These must be defined on the outer wavelength grid (binning)

!  Number of input phase function Legendre moments
!    THIS IS A BASIC INPUT THAT USER MUST PROVIDE

      INTEGER, INTENT(IN) ::          NMOMENTS_INPUT

!  Unscaled quantities, Elastic input

      REAL(FPK), INTENT(IN) :: DELTAU_INPUT_UNSCALED ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

!  Output arguments
!  ================

!  Scaled elastic-scattering optical properties
!  --------------------------------------------

!    These must be defined on the outer wavelength grid (binning)

!  Number of moments ( = 2 * NSTREAMS )

      INTEGER, INTENT(OUT) ::          NMOMENTS

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(OUT) :: DELTAU_VERT_INPUT   ( MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_ELASTIC &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Truncation factors (zero if no deltam scaling)

      REAL(FPK), INTENT(OUT) :: TRUNC_FACTORS         ( MAX_LAYERS, MAX_POINTS )

!  Unscaled Raman quantities
!  -------------------------

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSGAIN_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Scaled Raman-scattering optical properties
!  ------------------------------------------

!  Derived Cabannes scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_CABANNES &
           ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived rotational raman scattering (Internal variable). Loss term
!    Single scatter albedo * phase function moments
!    Only required for the Energy-balance approximation.

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSLOSS &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived rotational Raman scattering (Internal variable)

      REAL(FPK), INTENT(OUT) :: OMEGAMOMS_RRSGAIN &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Local Arrays
!  ============

!  shift dimensions

      INTEGER, PARAMETER :: MAXN2 = 48, MAXO2 = 185

!  Positions

      REAL(FPK) :: O2POS (MAXO2)
      REAL(FPK) :: N2POS (MAXN2)

!  Gammas

      REAL(FPK) :: GAMMA_O2(3)
      REAL(FPK) :: GAMMA_N2(3)

!  Basic coefficients for O2 and N2

      REAL(FPK) :: RRS_N2REF ( MAX_LAYERS, MAXN2 )
      REAL(FPK) :: RRS_O2REF ( MAX_LAYERS, MAXO2 )

!  T-derivatives of Basic coefficients for O2 and N2
!   ******** NOT REQUIRED HERE

      REAL(FPK) :: LINRRS_N2REF ( MAX_LAYERS, MAXN2 )
      REAL(FPK) :: LINRRS_O2REF ( MAX_LAYERS, MAXO2 )

!  Loss term cross-section

      REAL(FPK) :: RRSXSEC_OUT_MONO ( MAX_LAYERS )

!  Ranked Gain term cross-sections

      REAL(FPK) :: RRSXSEC_RANKED( MAX_LAYERS, MAX_POINTS )

!  Local flag. No derivatives here.

      LOGICAL, PARAMETER :: DO_WF = .false.

!  DELTA-M Scaling
!  ===============

!  REMARK: In Version 2.2, this is automatic, controlled only
!          by the DO_DELTAM_SCALING flag.
            
      CALL LRRS_DELTAM_SCALING_2P2 &
         ( MAX_LAYERS, MAX_MOMENTS_INPUT, MAX_MOMENTS, &
           MAX_POINTS, DO_DELTAM_SCALING, &
           NPOINTS_MONO, NLAYERS, NSTREAMS, NMOMENTS_INPUT, &
           DELTAU_INPUT_UNSCALED, OMEGAMOMS_ELASTIC_UNSCALED, &
           NMOMENTS, TRUNC_FACTORS, &
           DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC )

      !      write(*,*)DO_DELTAM_SCALING
      !      write(*,*)DELTAU_INPUT_UNSCALED(1,:)
      !      write(*,*)
      !      write(*,*)DELTAU_VERT_INPUT(1,:)
      !      write(*,*)
      !      write(*,*)OMEGAMOMS_ELASTIC_UNSCALED(1,2,:)
      !      write(*,*)
      !      write(*,*)OMEGAMOMS_ELASTIC(1,2,:)
      !      write(*,*)
      ! !      write(*,*)OMEGAMOMS_ELASTIC(1,1:20,1)
      !      pause
!  Generation of Raman-scattering Optical Properties (ROPS)
!  ========================================================

!  STEP 1. Raman Cross-sections generator
!  --------------------------------------

!  REMARK: THIS is an Automatic procedure in VERSION 2.2

!  Get the Raman spectroscopy

      CALL RRS_BASIC_SPECTROSCOPY &
         ( MAX_LAYERS, NLAYERS, DO_WF, LAYER_TEMPERATURES, MAXN2, MAXO2, &
           N2POS, O2POS, GAMMA_N2, GAMMA_O2, &
           RRS_N2REF, RRS_O2REF, LINRRS_N2REF, LINRRS_O2REF )
           
!  Given outer and inner windows and wavelengths + bins and solar fluxes
!  this routine generates the following Raman cross-sections:

!   Wavelength LAMBDA_EXCIT, --> Ranked wavelengths, and excitation poin
!                           (Ranked) Raman gain X-sections RRSXSEC_RANKE
!                                    Raman loss X-section RRSXSEC_OUT_MO

      call LRRS_RAMAN_XSECS_MONO &
        ( MAX_LAYERS, MAX_POINTS, NLAYERS, LAMBDA_EXCIT, &
          MAXN2, MAXO2, RRS_N2REF, RRS_O2REF, &
          N2POS, O2POS, GAMMA_N2, GAMMA_O2, &
          RRSXSEC_OUT_MONO, RRSXSEC_RANKED ) 
          ! RRSXSEC_OUT_MONO=LossTerm // RRSXSEC_RANKED=GainTerm
      !     write(*,*)RRSXSEC_RANKED(1,:)
      !     pause
!@@@@@@@@@ RT Solut
!  STEP 2. Make the Raman Optical Depths
!  -------------------------------------

!  Formerly, this was the responsibility of the User
!  Now taken care of automatically with these two routines:

!          lrrs_raman_ops_mono_plus (with linearizations);
!          lrrs_raman_ops_mono      (no linearization)

!  These routines are completely stand-alone.  In order to make them wor
!  the user must supply the following:

!      Dimensioning inputs
!      LRRS method (Energy-balancing vs. Cabannes). First is usual.
!      Number of layers, Inner window points and excitation position
!      Number of phase function expansion coefficients (Elastic)
!      Rayleigh X-secs(**) + depolarizations, Layer Air columns(**)
!      Elastic optical properties as defined above
!      Loss (RRSXSEC_OUT_MONO) and gain (RRSXSEC_RANKED) cross-sections(

!  (**) Note - cross sections must be [cm^2/mol], Air columns [mol/cm^2]

      call lrrs_raman_ops_mono &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS_INPUT, &
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, &
          NLAYERS, NMOMENTS_INPUT, NPOINTS_MONO, W_EXCIT, &
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, LAYER_AIRCOLUMNS, &
          DELTAU_INPUT_UNSCALED, OMEGAMOMS_ELASTIC_UNSCALED, &
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO, &
          OMEGAMOMS_CABANNES_UNSCALED, &
          OMEGAMOMS_RRSLOSS_UNSCALED, &
          OMEGAMOMS_RRSGAIN_UNSCALED )

      ! !     write(*,*)OMEGAMOMS_CABANNES_UNSCALED(1,0:5,1)
      !      write(*,*)OMEGAMOMS_RRSLOSS_UNSCALED(1,0,1)
      !      write(*,*)sum(OMEGAMOMS_RRSGAIN_UNSCALED(1,0,:))
      !     pause
!  Call the Raman IOP generator, scaled

      call lrrs_raman_ops_mono &
        ( MAX_LAYERS, MAX_POINTS, MAX_MOMENTS, &
          DO_ENERGY_BALANCING, DO_CABANNES_RAMAN, &
          NLAYERS, NMOMENTS, NPOINTS_MONO, W_EXCIT, &
          RAYLEIGH_XSEC, RAYLEIGH_DEPOL, LAYER_AIRCOLUMNS, &
          DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, &
          RRSXSEC_RANKED, RRSXSEC_OUT_MONO, &
          OMEGAMOMS_CABANNES, &
          OMEGAMOMS_RRSLOSS, &
          OMEGAMOMS_RRSGAIN )

!  debug
!        DO N = 1, NLAYERS
!          DO L = 0, NMOMENTS
!            DO CPT = 1, NPOINTS_MONO
!              write(87,*)n,l,cpt,OMEGAMOMS_ELASTIC(N,L,CPT)
!            ENDDO
!          ENDDO
!        ENDDO
!        pause

!  Finish

      return
      END SUBROUTINE SETUP_MASTER_MONO

      END MODULE lrrs_setup_master

