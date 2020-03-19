
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
! # Subroutines in this Module (Suffix '_1' no partials)        #
! #                                                             #
! #             LRRS_SSCORR_GENERAL_BIN_1  (master)             #
! #             LRRS_SSCORR_GENERAL_MONO_1 (master)             #
! #                                                             #
! #                 outgoing_attenuations_1                     #
! #                                                             #
! #                 outgoing_integration_bin_up_1               #
! #                 outgoing_integration_bin_dn_1               #
! #                                                             #
! #                 outgoing_integration_mono_up_1              #
! #                 outgoing_integration_mono_dn_1              #
! #                                                             #
! #   Single scatter exact calculations for outgoing LOS path   #
! #   Inelastic and Elastic output.                             #
! #                                                             #
! #   Programmed by R. Spurr, RT Solutions Inc.                 #
! #                                                             #
! #    First  Draft, April    2008                              #
! #    Second Draft, February 2011 Added Monochromatic routines #
! #    Third  Draft, March    2011 Use adjusted geometries      #
! #                                                             #
! ###############################################################

!  THIS IS THE LEVEL_BOUNDARY Version (No partials)

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Supplement-generated BRDF/SLEAVE inputs have replaced Ad-hoc inputs 
!    (2) Geometry routines have been separated to lrrs_geometry_1.f90, Now used here.
!    (3) Bookkeeping improvements (use of "Only", clearer I/O specifications)

!   -- Rob mod 5/12/17 for 2p5a, Switch to Phase-function inputs, Implement in code instead of Legendre expansions
!   -- Rob mod 5/12/17 for 2p5a, Incoporate SAVE_BEAM quantities for SSCORR_NADIR implementation
!   -- Rob mod 5/12/17 for 2p5a, Change subroutine names to reflect more general treatment.

      MODULE lrrs_corrections_1_m

!      USE LRRS_PARS_m, Only : SDU

      USE lrrs_geometry_1_m , Only : OUTGOING_SPHERGEOM_FINE_UP_1, &
                                     OUTGOING_SPHERGEOM_FINE_DN_1

      PRIVATE  :: outgoing_integration_bin_up_1 , outgoing_integration_bin_dn_1 , &
                  outgoing_integration_mono_up_1, outgoing_integration_mono_dn_1, &
                  outgoing_attenuations_1

      PUBLIC   :: LRRS_SSCORR_GENERAL_BIN_1, &
                  LRRS_SSCORR_GENERAL_MONO_1


      CONTAINS

      SUBROUTINE LRRS_SSCORR_GENERAL_BIN_1 &
         ( DO_FULL_RAMAN, DO_UPWELLING, DO_DNWELLING,                     & ! Input Mode flags
           DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,                           & ! Input Mode flags 2p5a
           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_DMSCALING,          & ! Input Mode flags
           DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,          & ! Input Surface flags
           DO_BRDF_Wav1, DO_Sleave_Wav1,                                  & ! Input Supplement flags
           NLAYERS, NFINELAYERS, N_USER_STREAMS, N_USER_RELAZMS,          & ! Input Control Integers
           N_LOUTPUT, LEVELMASK_UP, LEVELMASK_DN,                         & ! Input levels for output
           NPOINTS_INNER, OFFSET_INNER, NPOINTS_OUTER, N_RRSBINS, BINMAP, & ! Input Raman bin control
           SOLAR_ANGLE, COSSCAT_NADIR_UP, COSSCAT_NADIR_DN,               & ! Input Geometry
           SOLAR_ANGLES_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST,  & ! Input Geometry
           FLUXES_RANKED, ALBEDOS_RANKED, EARTH_RADIUS, HEIGHT_GRID,      & ! Input Fluxes/heights
           EXACTDB_BRDFUNC,  SLTERM_ISOTROPIC, SLTERM_USERANGLES,         & ! Input Surface stuff
           TRUNC_FACTORS, DELTAU_VERT_INPUT, OMEGA_UNSCALED,              & ! Input Optical Elastic
           OMEGAPHASFUNC_ELASTIC_UP,   OMEGAPHASFUNC_ELASTIC_DN,          & ! Input Optical Elastic
           OMEGAPHASFUNC_CABANNES_UP,  OMEGAPHASFUNC_CABANNES_DN,         & ! Input Optical Cabannes
           OMEGAMOMS_RRSLOSS_UNSCALED, OMEGAMOMS_RRSBIN_UNSCALED,         & ! Input Optical Raman
           SAVE_TRANS_USERM, SAVE_BEAMMULT_UP, SAVE_BEAMMULT_DN,          & ! Input Optical for NADIR SSCORR
           ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN,        & ! Outputs
           FAIL, MESSAGE)                                                   ! Outputs

!  Single scatter exact calculation for the outgoing LOS

!   Programmed by R. Spurr, RT Solutions Inc.
!    First Draft,  for LRRS Version 2.1, April 23rd 2008
!    Second draft, for LRRS Version 2.2, 23 July 2009

!   -- Rob mod 5/12/17 for 2p5a, Added DO_SSCORR_NADIR flag, plus precalculated multipliers/transmittances
!   -- Rob mod 5/12/17 for 2p5a, Using PHASFUNC input now, removed NMOMENTS_INPUT.

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, FOUR, PI4, DEG_TO_RAD,              &
                              MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES, &
                              MAX_LAYERS, MAX_FINE_LAYERS, MAX_LOUTPUT, MAX_POINTS, MAX_BINS

      IMPLICIT NONE

!  Input
!  -----

!  Flags, General

      LOGICAL  , INTENT(IN) :: DO_FULL_RAMAN
      LOGICAL  , INTENT(IN) :: DO_UPWELLING
      LOGICAL  , INTENT(IN) :: DO_DNWELLING

      LOGICAL  , INTENT(IN) :: DO_SSCORR_NADIR
      LOGICAL  , INTENT(IN) :: DO_SSCORR_OUTGOING

      LOGICAL  , INTENT(IN) :: DO_DMSCALING
      LOGICAL  , INTENT(IN) :: DO_CABANNES_RAMAN
      LOGICAL  , INTENT(IN) :: DO_ENERGY_BALANCING

!  New Version 2.5, Surface and supplement BRDF/SLEAVE Flags
!       Use of multiple or single point surface inputs. (Wav1 flags) 

      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE
      LOGICAL  , intent(in) :: DO_SURFACE_LEAVING
      LOGICAL  , intent(in) :: DO_SL_ISOTROPIC
      LOGICAL  , intent(in) :: DO_BRDF_Wav1
      LOGICAL  , intent(in) :: DO_SLEAVE_Wav1

!  Old Version 2.3 Code --------- commented out
!      LOGICAL  , INTENT(IN) :: DO_LAMBERTIAN_SURFACE
!      LOGICAL  , INTENT(IN) :: DO_WATERLEAVING
!      LOGICAL  , INTENT(IN) :: DO_LAMBERTIAN_FRACTION
!  Old Version 2.4 Code --------- commented out
!  Add SLTERM flag     . @@@ Rob Fix 06 Sep 12.
!  New Isotropic SLTERM. @@@ Rob Fix 06 sep 12
!      LOGICAL  , INTENT(IN)  :: DO_LRRS_SLTERM
!      REAL(FPK), INTENT(IN)  :: LRRS_SLTERM (MAX_POINTS)

!  number of layers

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of input phase function Legendre moments. THIS IS A BASIC INPUT THAT USER MUST PROVIDE
!   -- Rob mod 5/12/17 for 2p5a, No Longer required....!!!!!!!!
!      INTEGER  , INTENT(IN) :: NMOMENTS_INPUT

!  Number of fine layers for sscorr outgoing

      INTEGER  , INTENT(IN) :: NFINELAYERS

!  Los Control

      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      INTEGER  , INTENT(IN) :: N_USER_RELAZMS

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)

      INTEGER  , INTENT(IN) :: N_LOUTPUT

!  Layer masks for doing integrated source terms

      INTEGER  , INTENT(IN) :: LEVELMASK_UP (MAX_LOUTPUT)
      INTEGER  , INTENT(IN) :: LEVELMASK_DN (MAX_LOUTPUT)

!  outer/inner wavelength range, and offset for the inner range
!    Depends on the choice of solar spectrum

      INTEGER  , INTENT(IN) :: OFFSET_INNER
      INTEGER  , INTENT(IN) :: NPOINTS_INNER
      INTEGER  , INTENT(IN) :: NPOINTS_OUTER

!  Bin Mapping quantities (internal derivation)

      INTEGER  , INTENT(IN) :: N_RRSBINS  ( MAX_POINTS )
      INTEGER  , INTENT(IN) :: BINMAP ( MAX_BINS, MAX_POINTS )

!  Solar beam input geometry

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLES_ADJUST ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  User VZA angles
!   @@@ RobFix 5/5/11, Wrong dimensioning. Fixed 5/5/11

!      REAL(FPK), INTENT(IN) :: USER_ANGLES_ADJUST ( MAX_USER_RELAZMS )
      REAL(FPK), INTENT(IN) :: USER_ANGLES_ADJUST ( MAX_USER_STREAMS )

!  User azimuth angles

      REAL(FPK), INTENT(IN) :: USER_RELAZMS_ADJUST ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  original solar angle for NADIR correction

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLE

!  NADIR-view scattering angle cosines

      REAL(FPK), INTENT(IN) :: COSSCAT_NADIR_UP ( MAX_GEOMETRIES )
      REAL(FPK), INTENT(IN) :: COSSCAT_NADIR_DN ( MAX_GEOMETRIES )

!  solar fluxes

      REAL(FPK), INTENT(IN) :: FLUXES_RANKED   ( MAX_POINTS )

!  Lambertian albedos. New Version 2.5, 9/11/15

      REAL(fpk), intent(in) :: ALBEDOS_RANKED ( MAX_POINTS )

!  Earth Radius

      REAL(FPK), INTENT(IN) :: EARTH_RADIUS

!  Height grid

      REAL(FPK), INTENT(IN) :: HEIGHT_GRID ( 0:MAX_LAYERS )

!  These must be defined on the outer wavelength grid (binning)
!  Unscaled quantities, Elastic input
!   -- Rob mod 5/12/17 for 2p5a, Only need OMEGA_UNSCALED, plus PHASFUNC products

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: TRUNC_FACTORS     ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGA_UNSCALED    ( MAX_LAYERS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_UP  ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_DN  ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_CABANNES_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_CABANNES_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSLOSS_UNSCALED ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSBIN_UNSCALED  ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Exact (direct bounce) BRDF. Version 2.5 9/11/15

      REAL(fpk), intent(in) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

!  Isotropic Surface leaving term (if flag set). New Version 2.5, 9/11/15

      REAL(fpk), intent(in) :: SLTERM_ISOTROPIC ( MAX_POINTS )

!  Exact Surface-Leaving term. New Version 2.5, 9/11/15

      REAL(fpk), intent(in) :: SLTERM_USERANGLES (MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

!  NADIR correction - user stream transmittances, whole layers
!   -- Rob mod 5/12/17 for 2p5a

      REAL(FPK), INTENT(IN) :: SAVE_TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  NADIR correction - pre-calculated multipliers
!   -- Rob mod 5/12/17 for 2p5a

      REAL(FPK), INTENT(IN) :: SAVE_BEAMMULT_UP ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: SAVE_BEAMMULT_DN ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Output
!  ------

!  single scatter results

      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_UP (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_DN (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_UP (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_DN (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Exception handling

      LOGICAL          , INTENT(OUT) :: fail
      CHARACTER (LEN=*), INTENT(OUT) :: message

!  Geometry routine inputs and outputs
!  -----------------------------------

!  control

      LOGICAL   :: do_fine
      REAL(FPK) :: alpha_boa, theta_boa, phi_boa

!  main outputs (geometry)

      INTEGER   :: ntraverse(0:max_layers)
      REAL(FPK) :: sunpaths(0:max_layers,max_layers)
      REAL(FPK) :: radii   (0:max_layers)
      REAL(FPK) :: alpha_all  (0:max_layers)

!  Fine level output (geometry)

      INTEGER   :: ntraverse_fine(max_layers,max_fine_layers)
      REAL(FPK) :: sunpaths_fine (max_layers,max_layers,max_fine_layers)
      REAL(FPK) :: radii_fine    (max_layers,max_fine_layers)
      REAL(FPK) :: alpha_fine    (max_layers,max_fine_layers)

!  Other (incidental) geometrical output

      REAL(FPK) :: lospaths(max_layers)
      REAL(FPK) :: theta_all  (0:max_layers)
      REAL(FPK) :: phi_all    (0:max_layers)
      REAL(FPK) :: cosscat_up (0:max_layers)
      REAL(FPK) :: cosscat_dn (0:max_layers)

!  Extinction

      REAL(FPK) :: extinction ( max_layers, max_points )

!  Saved Legendre polynomials. 
!   -- Rob mod 5/12/17 for 2p5a, Only need 0:2 now, drop MAX_MOMENTS_INPUT

      REAL(FPK) :: SS_PLEG_UP(MAX_LAYERS,0:2), SS_PLEG_DN(MAX_LAYERS,0:2)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(FPK) :: TMS(MAX_LAYERS,MAX_POINTS)

!  Local truncation factors for additional DELTAM scaling. DISABLED
!      REAL(FPK) :: SSFDEL ( MAXLAYERS )

!  Exact Phase function calculations

      REAL(FPK) :: ESCAT_UP(MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ESCAT_DN(MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ISCAT_UP(MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ISCAT_DN(MAX_LAYERS,MAX_POINTS)

!  Cumulative single scatter source terms

      REAL(FPK) :: ESS_CUMSOURCE_UP(0:MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ESS_CUMSOURCE_DN(0:MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ISS_CUMSOURCE_UP(0:MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ISS_CUMSOURCE_DN(0:MAX_LAYERS,MAX_POINTS)

!  Outgoing sphericity stuff
!  -------------------------

!  Multipliers and LOS transmittances

      REAL(FPK) :: &
         UP_IMULTIPLIERS(MAX_LAYERS,MAX_BINS,MAX_POINTS), &
         DN_IMULTIPLIERS(MAX_LAYERS,MAX_BINS,MAX_POINTS), &
         UP_EMULTIPLIERS(MAX_LAYERS,MAX_POINTS), &
         DN_EMULTIPLIERS(MAX_LAYERS,MAX_POINTS), &
         UP_LOSTRANS(MAX_LAYERS,MAX_POINTS), &
         DN_LOSTRANS(MAX_LAYERS,MAX_POINTS)

!  Solar beam attenuations to BOA (required for exact DB calculation)

      REAL(FPK) :: attn      ( 0:max_layers, max_points )
      REAL(FPK) :: attn_fine ( max_layers, max_fine_layers, max_points )

!  Local variables
!  ---------------

!  Indices

      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, W, B
      INTEGER   :: UTA, UM, IA, NC, V, WR, WB
      INTEGER   :: Sleave_IDX, Brdf_IDX

!  Help variables (REAL(FPK) ::)
!   -- Rob mod 5/12/17 for 2p5a, removed DFL1, DFL2

      REAL(FPK) :: boa_attn, ssflux, x0_boa, x0_boa_4, ratio, Brdf_Exact
      REAL(FPK) :: SOURCE, HELP, COSSCAT, CONTRIB
      REAL(FPK) :: PHAS_E, PHAS_LR, PHAS_LG(MAX_BINS), PHAS_CB
      REAL(FPK) :: FACTOR, CAB, GAIN, LOSS
      REAL(FPK) :: REFLEC(MAX_POINTS), Sleave_Iso, Sleave_User

!  Local surface flag

      LOGICAL  , PARAMETER :: DO_INCLUDE_SURFACE = .true.

!  Set up operations
!  -----------------

!  Ss flux scaling factor

      SSFLUX = one / PI4

!  Fine layering flag

      do_fine = .true.

!   Number of fine layers now a control input (7/25/08)

!  Floating point numbers for Legendre polynomials
!   -- Rob mod 5/12/17 for 2p5a, removed DFL1, DFL2
!      DO L = 2, NMOMENTS_INPUT
!        HELP = DBLE(L) ; DF1(L) = DBLE(2*L-1)/HELP ; DF2(L) = DBLE(L-1)/HELP
!      ENDDO

!  Create TMS factors, these get stored
!    Delta-M Scaling introduced April 2005.
!   -- Rob mod 5/12/17 for 2p5a, Use OMEGA_UNSCALED

      DO N = 1, NLAYERS
        DO W = 1, NPOINTS_INNER
          TMS(N,W) = ONE
          IF ( DO_DMSCALING ) THEN
            WR = W + OFFSET_INNER
            HELP = TRUNC_FACTORS(N,WR) * OMEGA_UNSCALED(N,WR)
            TMS(N,W) = TMS(N,W) / (ONE - HELP)
          ENDIF
        ENDDO
      ENDDO

!  Create extinctions (using scaled optical depths)

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        DO W = 1, NPOINTS_OUTER
          EXTINCTION(N,W) = DELTAU_VERT_INPUT(N,W) / HELP
        ENDDO
      ENDDO

!  Additional Delta-M scaling. DISABLED
!  ------------------------------------

!  New section. R. Spurr, 07 September 2007.
!   TMS gets modified by (1-F). Save the truncation factor.
!   Phase function moments are modified later on.
!      IF ( DO_SSCORR_TRUNCATION ) THEN
!        NM1  = NMOMENTS_INPUT
!        DNM1 = DBLE(2*NM1+1)
!        DO N = 1, NLAYERS
!          BETA = OMEGAMOMS_ELASTIC  (N, NM1, W_EXCIT ) / OMEGA_LOCAL(N)
!          SSFDEL(N) = BETA / DNM1
!          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
!        ENDDO
!      ENDIF

!  Start the main loop over all viewing geometries
!  -----------------------------------------------

!    Outgoing correction, use the adjusted values of the angles -------> March 2011
!   -- Rob mod 5/12/17 for 2p5a, re-introduce the SSCORR_NADIR computation

      DO UM = 1, N_USER_STREAMS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  UPWELLING NADIR CORRECTION: multipliers, transmittances
!  =======================================================

!  Upwelling Multipliers and Transmittances for the SSCORR_NADIR case
!     Does not use the ADJUSTED geometrical angles, Same for each azimuth.
!   -- Rob mod 5/12/17 for 2p5a, SSCORR_NADIR values are already pre-calculated

        IF ( DO_SSCORR_NADIR .and.DO_UPWELLING ) THEN
          DO N = 1, NLAYERS
            DO W = 1, NPOINTS_INNER
              WR = W + OFFSET_INNER
              UP_LOSTRANS(N,W) = SAVE_TRANS_USERM(UM,N,W)
              UP_EMULTIPLIERS(N,W) = SAVE_BEAMMULT_UP(UM,N,W)
              DO B = 1, N_RRSBINS(W)
                WB = BINMAP(B,W)
                UP_IMULTIPLIERS(N,B,W) = SAVE_BEAMMULT_UP(UM,N,WB)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  DOWNWELLING NADIR CORRECTION: multipliers, transmittances
!  =========================================================

!  Downwelling Multipliers and Transmittances for the SSCORR_NADIR case
!     Does not use the ADJUSTED geometrical angles, Same for each azimuth.
!   -- Rob mod 5/12/17 for 2p5a, SSCORR_NADIR values are already pre-calculated

        IF ( DO_SSCORR_NADIR .and.DO_DNWELLING ) THEN
          DO N = 1, NLAYERS
            DO W = 1, NPOINTS_INNER
              WR = W + OFFSET_INNER
              DN_LOSTRANS(N,W) = SAVE_TRANS_USERM(UM,N,W)
              DN_EMULTIPLIERS(N,W) = SAVE_BEAMMULT_DN(UM,N,W)
              DO B = 1, N_RRSBINS(W)
                WB = BINMAP(B,W)
                DN_IMULTIPLIERS(N,B,W) = SAVE_BEAMMULT_DN(UM,N,WB)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  Start the Azimuth loop
!  ----------------------

        DO IA = 1, N_USER_RELAZMS

!  Adjusted geometry for the Outgoing correction

          ALPHA_BOA = USER_ANGLES_ADJUST(UM)
          THETA_BOA = SOLAR_ANGLES_ADJUST(UM,IA)
          PHI_BOA   = USER_RELAZMS_ADJUST(UM,IA)
          V = N_USER_RELAZMS * (UM-1) + IA

!  UPWELLING OUTGOING CORRECTION: multipliers, transmittances
!  ==========================================================

          IF ( DO_UPWELLING .and. DO_SSCORR_OUTGOING ) THEN

!  Call to geometry routine, path distances etc....

            call outgoing_sphergeom_fine_up_1 &
        ( max_layers, max_fine_layers, do_fine, nlayers, nfinelayers, & ! Inputs
          height_grid, earth_radius, alpha_boa, theta_boa, phi_boa,   & ! Inputs
          sunpaths,      radii,      ntraverse,      alpha_all,       & ! Inputs
          sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine,      & ! Inputs
          lospaths, theta_all, phi_all, cosscat_up,                   & ! Inputs
          fail, message )                                               ! Outputs

!  Exception handling

            if ( fail ) return

!  Get the attenuations

            call outgoing_attenuations_1 &
        ( max_layers, max_fine_layers, max_points, & ! Inputs
          npoints_outer, nlayers, nfinelayers,     & ! Inputs
          extinction, sunpaths, ntraverse,         & ! Inputs
          sunpaths_fine, ntraverse_fine,           & ! Inputs
          attn, attn_fine )                          ! Outputs

!  Multipliers, transmittances

            call outgoing_integration_bin_up_1 &
        ( max_layers, max_fine_layers, max_points, max_bins, & ! Inputs
          npoints_inner, offset_inner, n_rrsbins, binmap,    & ! Inputs
          nlayers, nfinelayers, extinction,                  & ! Inputs
          radii, alpha_all, alpha_fine, attn, attn_fine,     & ! Inputs
          up_emultipliers, up_imultipliers, up_lostrans )     ! Outputs

!   End Outgoing correction multipliers etc.

          endif

!  UPWELLING SINGLE SCATTER CORRECTION:  Source terms
!  ==================================================

          IF ( DO_UPWELLING ) THEN

!   Cosine scatter angles
!   -- Rob mod 5/12/17 for 2p5a, Careful of selection.

            if ( DO_SSCORR_NADIR ) THEN
              COSSCAT = COSSCAT_NADIR_UP(V)
            else
              COSSCAT = COSSCAT_UP(NLAYERS)
            endif

!  legendre polynomials
!   -- Rob mod 5/12/17 for 2p5a, Only need #2 = (3x^2-1)/2

            SS_PLEG_UP(1,0) = ONE
            SS_PLEG_UP(1,1) = COSSCAT
            SS_PLEG_UP(1,2) = 1.5_fpk * COSSCAT * COSSCAT - 0.5_fpk
!            DO L = 2, NMOMENTS_INPUT
!               SS_PLEG_UP(1,L) = DF1(L) * SS_PLEG_UP(1,L-1) * COSSCAT - DF2(L) * SS_PLEG_UP(1,L-2)
!            ENDDO

!  Elastic scattering Phase functions (multiplied by TMS factor)
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_E with Phase function input
!   -- PHAS_E = DOT_PRODUCT(OMEGAMOMS_ELASTIC_UNSCALED(N,0:NMOMENTS_INPUT,WR),SS_PLEG_UP(1,0:NMOMENTS_INPUT))

            DO N = 1, NLAYERS
              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
                PHAS_E = OMEGAPHASFUNC_ELASTIC_UP(N,V,WR) * TMS(N,W)
                ESCAT_UP(N,W) = PHAS_E * UP_EMULTIPLIERS(N,W)
              ENDDO
            ENDDO

!  Add inelastic scattering Phase functions if flagged

            IF ( DO_FULL_RAMAN ) THEN

!  Energy balancing method
!  -----------------------

             IF ( DO_ENERGY_BALANCING ) THEN

!  Start loops

              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
                DO N = 1, NLAYERS

!  Loss term

                  PHAS_LR = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,W) + &
                            OMEGAMOMS_RRSLOSS_UNSCALED(N,2,W) * SS_PLEG_UP(1,2)
                  PHAS_LR = PHAS_LR * TMS(N,W)
                  LOSS = PHAS_LR * UP_EMULTIPLIERS(N,W)

!  Gain term

                  GAIN = ZERO
                  DO B = 1, N_RRSBINS(W)
                    WB = BINMAP(B,W)
                    RATIO   = FLUXES_RANKED(WB)/FLUXES_RANKED(WR)
                    CONTRIB = OMEGAMOMS_RRSBIN_UNSCALED(N,0,B,W) + &
                              OMEGAMOMS_RRSBIN_UNSCALED(N,2,B,W) * SS_PLEG_UP(1,2)
                    PHAS_LG(B) = CONTRIB * RATIO
                    GAIN = GAIN + PHAS_LG(B) * UP_IMULTIPLIERS(N,B,W)
                  ENDDO
                  GAIN = GAIN * TMS(N,W)

!  Whole layer contributions

                  ISCAT_UP(N,W) = ESCAT_UP(N,W) + GAIN - LOSS

!  End loops

                ENDDO
              ENDDO

!  Cabannes Raman method
!  ---------------------

             ELSE IF ( DO_CABANNES_RAMAN ) THEN

!  Start loops

               DO N = 1, NLAYERS
                 DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER

!  Cabannes term
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_CB with Phase function input
!   -- PHAS_CB = DOT_PRODUCT(OMEGAMOMS_CABANNES_UNSCALED(N,0:NMOMENTS_INPUT,W),SS_PLEG_UP(1,0:NMOMENTS_INPUT))

                  PHAS_CB = OMEGAPHASFUNC_CABANNES_UP(N,V,W) * TMS(N,W)
                  CAB = PHAS_CB * UP_EMULTIPLIERS(N,W)

!  Gain term

                  GAIN = ZERO
                  DO B = 1, N_RRSBINS(W)
                    WB = BINMAP(B,W)
                    RATIO   = FLUXES_RANKED(WB)/FLUXES_RANKED(WR)
                    CONTRIB = OMEGAMOMS_RRSBIN_UNSCALED(N,0,B,W) + &
                              OMEGAMOMS_RRSBIN_UNSCALED(N,2,B,W) * SS_PLEG_UP(1,2)
                    PHAS_LG(B)  = CONTRIB * RATIO
                    GAIN = GAIN + PHAS_LG(B) * UP_IMULTIPLIERS(N,B,W)
                  ENDDO
                  GAIN = GAIN * TMS(N,W)

!  Whole layer contribution

                  ISCAT_UP(N,W) = CAB + GAIN

!  End loops

                ENDDO
              ENDDO

!  End inelastic scattering options

             ENDIF
            ENDIF

!  End upwelling clause for Source terms

          ENDIF

!  DOWNWELLING OUTGOING CORRECTION: multipliers, transmittances
!  ============================================================

          IF ( DO_DNWELLING .and. DO_SSCORR_OUTGOING ) THEN

!  Call to geometry routine, path distances etc....

            call outgoing_sphergeom_fine_dn_1 &
        ( max_layers, max_fine_layers, do_fine, nlayers, nfinelayers, & ! Inputs
          height_grid, earth_radius, alpha_boa, theta_boa, phi_boa,   & ! Inputs
          sunpaths,      radii,      ntraverse,      alpha_all,       & ! Inputs
          sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine,      & ! Inputs
          lospaths, theta_all, phi_all, cosscat_dn,                   & ! Inputs
          fail, message )                                               ! Outputs

!  Exception handling

            if ( fail ) return

!  Get the attenuations

            call outgoing_attenuations_1 &
        ( max_layers, max_fine_layers, max_points, & ! Inputs
          npoints_outer, nlayers, nfinelayers,     & ! Inputs
          extinction, sunpaths, ntraverse,         & ! Inputs
          sunpaths_fine, ntraverse_fine,           & ! Inputs
          attn, attn_fine )                          ! Outputs

!  Multipliers and transmittances

            call outgoing_integration_bin_dn_1 &
        ( max_layers, max_fine_layers, max_points, max_bins, & ! Inputs
          npoints_inner, offset_inner, n_rrsbins, binmap,    & ! Inputs
          nlayers, nfinelayers, extinction,                  & ! Inputs
          radii, alpha_all, alpha_fine, attn, attn_fine,     & ! Inputs
          dn_emultipliers, dn_imultipliers, dn_lostrans )      ! Outputs

!   End Outgoing correction multipliers etc.

          endif

!  DOWNWELLING SINGLE SCATTER CORRECTION:  Source terms
!  ====================================================

          IF ( DO_DNWELLING ) THEN

!   Cosine scatter angles
!   -- Rob mod 5/12/17 for 2p5a, Careful of selection.

            if ( DO_SSCORR_NADIR ) THEN
              COSSCAT = COSSCAT_NADIR_DN(V)
            else
              COSSCAT = COSSCAT_DN(NLAYERS)
            endif

!  Legendre polynomials
!   -- Rob mod 5/12/17 for 2p5a, Only need #2 = (3x^2-1)/2

            SS_PLEG_DN(1,0) = ONE
            SS_PLEG_DN(1,1) = COSSCAT
            SS_PLEG_DN(1,2) = 1.5_fpk * COSSCAT * COSSCAT - 0.5_fpk
!            DO L = 2, NMOMENTS_INPUT
!               SS_PLEG_DN(1,L) = DF1(L) * SS_PLEG_DN(1,L-1) * COSSCAT - DF2(L) * SS_PLEG_DN(1,L-2)
!            ENDDO

!  Elastic scattering Phase functions (multiplied by TMS factor)
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_E with Phase function input
!   -- PHAS_E = DOT_PRODUCT(OMEGAMOMS_ELASTIC_UNSCALED(N,0:NMOMENTS_INPUT,WR),SS_PLEG_DN(1,0:NMOMENTS_INPUT))

            DO N = 1, NLAYERS
              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
                PHAS_E = OMEGAPHASFUNC_ELASTIC_DN(N,V,WR) * TMS(N,W)
                ESCAT_DN(N,W) = PHAS_E * DN_EMULTIPLIERS(N,W)
              ENDDO
            ENDDO

!  Add inelastic scattering Phase functions if flagged

            IF ( DO_FULL_RAMAN ) THEN

!  Energy balancing method
!  -----------------------

             IF ( DO_ENERGY_BALANCING ) THEN

!  Start loops

              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
                DO N = 1, NLAYERS

!  Loss term

                  PHAS_LR = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,W) + &
                            OMEGAMOMS_RRSLOSS_UNSCALED(N,2,W) * SS_PLEG_DN(1,2)
                  PHAS_LR = PHAS_LR * TMS(N,W)
                  LOSS = PHAS_LR * DN_EMULTIPLIERS(N,W)

!  Gain term

                  GAIN = ZERO
                  DO B = 1, N_RRSBINS(W)
                    WB = BINMAP(B,W)
                    RATIO   = FLUXES_RANKED(WB)/FLUXES_RANKED(WR)
                    CONTRIB = OMEGAMOMS_RRSBIN_UNSCALED(N,0,B,W) + &
                              OMEGAMOMS_RRSBIN_UNSCALED(N,2,B,W) * SS_PLEG_DN(1,2)
                    PHAS_LG(B) = CONTRIB * RATIO
                    GAIN = GAIN + PHAS_LG(B) * DN_IMULTIPLIERS(N,B,W)
                  ENDDO
                  GAIN = GAIN * TMS(N,W)

!  Whole layer contributions

                ISCAT_DN(N,W) = ESCAT_DN(N,W) + GAIN - LOSS

!  End loops

                ENDDO
              ENDDO

!  Cabannes Raman method
!  ---------------------

             ELSE IF ( DO_CABANNES_RAMAN ) THEN

!  Start loops

               DO N = 1, NLAYERS
                 DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER

!  Cabannes term
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_CB with Phase function input
!   -- PHAS_CB = DOT_PRODUCT(OMEGAMOMS_CABANNES_UNSCALED(N,0:NMOMENTS_INPUT,W),SS_PLEG_DN(1,0:NMOMENTS_INPUT))

                  PHAS_CB = OMEGAPHASFUNC_CABANNES_DN(N,V,W) * TMS(N,W)
                  CAB = PHAS_CB * DN_EMULTIPLIERS(N,W)

!  Gain term

                  GAIN = ZERO
                  DO B = 1, N_RRSBINS(W)
                    WB = BINMAP(B,W)
                    RATIO   = FLUXES_RANKED(WB)/FLUXES_RANKED(WR)
                    CONTRIB = OMEGAMOMS_RRSBIN_UNSCALED(N,0,B,W) + &
                              OMEGAMOMS_RRSBIN_UNSCALED(N,2,B,W) * SS_PLEG_DN(1,2)
                    PHAS_LG(B)  = CONTRIB * RATIO
                    GAIN = GAIN + PHAS_LG(B) * DN_IMULTIPLIERS(N,B,W)
                  ENDDO
                  GAIN = GAIN * TMS(N,W)

!  Whole layer contribution

                  ISCAT_DN(N,W) = CAB + GAIN

!  End loops

                ENDDO
              ENDDO

!  End inelastic scattering options

             ENDIF
            ENDIF

!  End Downwelling clause for Source terms

          ENDIF

!  Recurrence relation for the UPWELLING intensity
!  ===============================================

          IF ( DO_UPWELLING ) THEN

!  This section rewritten for Version 2.5

            IF ( DO_INCLUDE_SURFACE ) THEN
              IF ( DO_BRDF_SURFACE ) THEN
                IF ( DO_BRDF_Wav1 ) then
                  Brdf_IDX = 1 ; BRDF_exact = EXACTDB_BRDFUNC(UM,IA,Brdf_IDX)
                  DO W = 1, NPOINTS_INNER
                    REFLEC(W) = BRDF_exact
                  ENDDO
                ELSE
                  DO W = 1, NPOINTS_INNER
                    WR = W + OFFSET_INNER ; Brdf_IDX = WR
                    REFLEC(W) = EXACTDB_BRDFUNC(UM,IA,Brdf_IDX)
                  ENDDO
                ENDIF
              ELSE
                DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER
                  REFLEC(W) = ALBEDOS_RANKED(WR)
                ENDDO
              ENDIF
            ENDIF

!  initialize cumulative source terms = F.A.mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux, A = albedo
!    Only require the Stokes total intensity component.
!    This contribution is purely elastic.
!    This line Replaced: X0_BOA = cos(SOLAR_ANGLE * DEG_TO_RAD)  . reintroduced for 2p5a
!     Introduction of the Adjusted Gometry value, 18 March 2011
!    
            NC =  0
            IF ( DO_INCLUDE_SURFACE ) THEN
              IF ( DO_SSCORR_NADIR ) then
                X0_BOA = cos(SOLAR_ANGLE * DEG_TO_RAD)
              ELSE
                X0_BOA = cos(SOLAR_ANGLES_ADJUST(UM,IA) * DEG_TO_RAD)
              ENDIF
              X0_BOA_4 = FOUR * X0_BOA
              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
                FACTOR = REFLEC(W)
                BOA_ATTN = X0_BOA_4 * ATTN(NLAYERS,WR)
                ESS_CUMSOURCE_UP(NC,W) = FACTOR * BOA_ATTN
              ENDDO
            ELSE
              DO W = 1, NPOINTS_INNER
                ESS_CUMSOURCE_UP(NC,W) = ZERO
              ENDDO
            ENDIF

!  Version 2.4 code commented out. 9/11/15
!  Add SLTERM flag     . @@@ Rob Fix 06 Sep 12.
!  New Isotropic SLTERM. @@@ Rob Fix 06 sep 12
!   Multiply terms by 4pi to get correct normalization. 7/30/12
!            IF ( DO_LRRS_SLTERM ) THEN
!              DO W = 1, NPOINTS_INNER
!                WR = W + OFFSET_INNER
!                ESS_CUMSOURCE_UP(NC,W) = &
!                     ESS_CUMSOURCE_UP(NC,W) + LRRS_SLTERM(WR)*PI4
!              ENDDO
!            ENDIF

!  New Surface-Leaving stuff. Version 2.5, 9/11/15
!    Multiply terms by 4pi to get correct normalization.

            IF ( DO_SURFACE_LEAVING ) THEN
              IF ( DO_SL_ISOTROPIC ) THEN
                IF ( DO_SLEAVE_Wav1 ) then
                  Sleave_IDX = 1 ; Sleave_Iso = SLTERM_ISOTROPIC(Sleave_IDX) * PI4
                  DO W = 1, NPOINTS_INNER
                    ESS_CUMSOURCE_UP(NC,W) = ESS_CUMSOURCE_UP(NC,W) + Sleave_Iso
                  ENDDO
                ELSE
                  DO W = 1, NPOINTS_INNER
                    WR = W + OFFSET_INNER ; Sleave_IDX = WR
                    ESS_CUMSOURCE_UP(NC,W) = ESS_CUMSOURCE_UP(NC,W) + SLTERM_ISOTROPIC(Sleave_IDX)*PI4
                  ENDDO
                ENDIF
              ELSE
                IF ( DO_SLEAVE_Wav1 ) then
                  Sleave_IDX = 1 ; Sleave_User = SLTERM_USERANGLES(UM,IA,Sleave_IDX) * PI4
                  DO W = 1, NPOINTS_INNER
                    ESS_CUMSOURCE_UP(NC,W) = ESS_CUMSOURCE_UP(NC,W) + Sleave_User
                  ENDDO
                ELSE
                  DO W = 1, NPOINTS_INNER
                    WR = W + OFFSET_INNER ; Sleave_IDX = WR
                    ESS_CUMSOURCE_UP(NC,W) = ESS_CUMSOURCE_UP(NC,W) + SLTERM_USERANGLES(UM,IA,Sleave_IDX)*PI4
                  ENDDO
                ENDIF
              ENDIF
            ENDIF

!  Duplicate the result for the inelastic case

            IF ( DO_FULL_RAMAN ) THEN
              DO W = 1, NPOINTS_INNER
                ISS_CUMSOURCE_UP(NC,W) = ESS_CUMSOURCE_UP(NC,W)
              ENDDO
            ENDIF

!  initialize optical depth loop

            NSTART = NLAYERS
            NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

            DO UTA = N_LOUTPUT, 1, -1

!  Layer index for given optical depth

              NLEVEL = LEVELMASK_UP(UTA)
              NUT    = NLEVEL + 1

!  Cumulative single scatter source terms :
!      For loop over layers working upwards to level NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier
!  Multiplier using new integration scheme

              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N
                DO W = 1, NPOINTS_INNER
                  ESS_CUMSOURCE_UP(NC,W) =  ESCAT_UP(N,W) + &
                       UP_LOSTRANS(N,W) * ESS_CUMSOURCE_UP(NC-1,W)
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO W = 1, NPOINTS_INNER
                    ISS_CUMSOURCE_UP(NC,W) = ISCAT_UP(N,W) + &
                       UP_LOSTRANS(N,W) * ISS_CUMSOURCE_UP(NC-1,W)
                  ENDDO
                ENDIF
              ENDDO

!  Ongrid output--------
!     Set final cumulative source and single scatter intensity

              DO W = 1, NPOINTS_INNER
                SOURCE = ESS_CUMSOURCE_UP(NC,W)
                ELASTIC_SS_UP(UTA,V,W) = SSFLUX * SOURCE
              ENDDO
              IF ( DO_FULL_RAMAN ) THEN
                DO W = 1, NPOINTS_INNER
                  SOURCE = ISS_CUMSOURCE_UP(NC,W)
                  RAMAN_SS_UP(UTA,V,W) = SSFLUX * SOURCE
                ENDDO
              ENDIF

!  Check for updating the recursion

              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
              NUT_PREV = NUT

!  end optical depth loop and Upwelling clause

            ENDDO
          ENDIF

!  Recurrence relation for the DOWNWELLING intensity
!  =================================================

          IF ( DO_DNWELLING ) THEN

!  initialize cumulative source terms, and optical depth loop

            NC =  0
            DO W = 1, NPOINTS_INNER
              ESS_CUMSOURCE_DN(NC,W) = ZERO
            ENDDO
            IF ( DO_FULL_RAMAN ) THEN
              DO W = 1, NPOINTS_INNER
                ISS_CUMSOURCE_DN(NC,W) = ZERO
              ENDDO
            ENDIF

!  Start loop over optical depths

            NSTART = 1
            NUT_PREV = NSTART - 1

!  Main loop over all output optical depths

            DO UTA = 1, N_LOUTPUT

!  Layer index for given optical depth

              NLEVEL = LEVELMASK_DN(UTA)
              NUT = NLEVEL

!  Cumulative single scatter source terms :
!      For loop over layers working downwards to NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier
!      Multiplier by new integration method

              DO N = NSTART, NUT
                NC = N
                DO W = 1, NPOINTS_INNER
                  ESS_CUMSOURCE_DN(NC,W) = ESCAT_DN(N,W) + &
                       DN_LOSTRANS(N,W) * ESS_CUMSOURCE_DN(NC-1,W)
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO W = 1, NPOINTS_INNER
                    ISS_CUMSOURCE_DN(NC,W) = ISCAT_DN(N,W)+ &
                       DN_LOSTRANS(N,W) * ISS_CUMSOURCE_DN(NC-1,W)
                  ENDDO
                ENDIF
              ENDDO

!  Ongrid output :
!     Set final cumulative source and correct Stokes vector

              DO W = 1, NPOINTS_INNER
                SOURCE = ESS_CUMSOURCE_DN(NC,W)
                ELASTIC_SS_DN(UTA,V,W) = SSFLUX * SOURCE
              ENDDO
              IF ( DO_FULL_RAMAN ) THEN
                DO W = 1, NPOINTS_INNER
                  SOURCE = ISS_CUMSOURCE_DN(NC,W)
                  RAMAN_SS_DN(UTA,V,W) = SSFLUX * SOURCE
                ENDDO
              ENDIF

!  Check for updating the recursion

              IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
              NUT_PREV = NUT

!  End optical depth loop and Downwelling clause

            ENDDO
          ENDIF

!  Finish geometry loop

        enddo
      enddo

!  Finish

      RETURN
      END SUBROUTINE LRRS_SSCORR_GENERAL_BIN_1

!

      SUBROUTINE LRRS_SSCORR_GENERAL_MONO_1 &
         ( DO_FULL_RAMAN, DO_UPWELLING, DO_DNWELLING,                    & ! Inputs
           DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,                          & ! Inputs
           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_DMSCALING,         & ! Inputs
           DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,         & ! Inputs
           DO_BRDF_Wav1, DO_Sleave_Wav1,                                 & ! Inputs
           NLAYERS, NFINELAYERS, N_USER_STREAMS, N_USER_RELAZMS,         & ! Input Control Integers
           N_LOUTPUT, LEVELMASK_UP, LEVELMASK_DN,                        & ! Input levels for output
           NPOINTS_MONO, W_EXCIT,                                        & ! Inputs Raman Spec.
           SOLAR_ANGLE, COSSCAT_NADIR_UP, COSSCAT_NADIR_DN,              & ! Input Geometry
           SOLAR_ANGLES_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST, & ! Input Geometry
           FLUXES_RANKED, ALBEDOS_RANKED, EARTH_RADIUS, HEIGHT_GRID,     & ! Input Fluxes/heights
           EXACTDB_BRDFUNC, SLTERM_ISOTROPIC, SLTERM_USERANGLES,         & ! Input Surface stuff
           TRUNC_FACTORS, DELTAU_VERT_INPUT, OMEGA_UNSCALED,             & ! Input Optical Elastic
           OMEGAPHASFUNC_ELASTIC_UP,   OMEGAPHASFUNC_ELASTIC_DN,         & ! Input Optical Elastic
           OMEGAPHASFUNC_CABANNES_UP,  OMEGAPHASFUNC_CABANNES_DN,        & ! Input Optical Cabannes
           OMEGAMOMS_RRSLOSS_UNSCALED, OMEGAMOMS_RRSGAIN_UNSCALED,       & ! Input Optical Raman
           SAVE_TRANS_USERM, SAVE_BEAMMULT_UP, SAVE_BEAMMULT_DN,         & ! Input Optical for NADIR SSCORR
           ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN,       & ! Outputs
           FAIL, MESSAGE)                                                  ! Outputs

!  Single scatter exact calculation for the outgoing LOS

!   Programmed by R. Spurr, RT Solutions Inc.
!    First Draft,  for LRRS Version 2.1, April 23rd 2008
!    Second draft, for LRRS Version 2.2, 23 July 2009

!   -- Rob mod 5/12/17 for 2p5a, Added DO_SSCORR_NADIR flag, plus precalculated multipliers/transmittances
!   -- Rob mod 5/12/17 for 2p5a, Using PHASFUNC input now, removed NMOMENTS_INPUT.

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, FOUR, PI4, DEG_TO_RAD,              &
                              MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES, &
                              MAX_LAYERS, MAX_FINE_LAYERS, MAX_LOUTPUT, MAX_POINTS

      IMPLICIT NONE

!  Input
!  -----

!  Flags, General

      LOGICAL  , INTENT(IN) :: DO_FULL_RAMAN
      LOGICAL  , INTENT(IN) :: DO_UPWELLING
      LOGICAL  , INTENT(IN) :: DO_DNWELLING

      LOGICAL  , INTENT(IN) :: DO_SSCORR_NADIR
      LOGICAL  , INTENT(IN) :: DO_SSCORR_OUTGOING

      LOGICAL  , INTENT(IN) :: DO_DMSCALING
      LOGICAL  , INTENT(IN) :: DO_CABANNES_RAMAN
      LOGICAL  , INTENT(IN) :: DO_ENERGY_BALANCING

!  New Version 2.5, Surface and supplement BRDF/SLEAVE Flags
!       Use of multiple or single point surface inputs. (Wav1 flags) 

      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE
      LOGICAL  , intent(in) :: DO_SURFACE_LEAVING
      LOGICAL  , intent(in) :: DO_SL_ISOTROPIC
      LOGICAL  , intent(in) :: DO_BRDF_Wav1
      LOGICAL  , intent(in) :: DO_SLEAVE_Wav1

!  Old Version 2.3 Code --------- commented out
!      LOGICAL  , INTENT(IN) :: DO_LAMBERTIAN_SURFACE
!      LOGICAL  , INTENT(IN) :: DO_WATERLEAVING
!      LOGICAL  , INTENT(IN) :: DO_LAMBERTIAN_FRACTION
!  Old Version 2.4 Code --------- commented out
!  Add SLTERM flag     . @@@ Rob Fix 06 Sep 12.
!  New Isotropic SLTERM. @@@ Rob Fix 06 sep 12
!      LOGICAL  , INTENT(IN)  :: DO_LRRS_SLTERM
!      REAL(FPK), INTENT(IN)  :: LRRS_SLTERM (MAX_POINTS)

!  number of layers

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of input phase function Legendre moments. THIS IS A BASIC INPUT THAT USER MUST PROVIDE
!   -- Rob mod 5/12/17 for 2p5a, No Longer required....!!!!!!!!
!      INTEGER  , INTENT(IN) :: NMOMENTS_INPUT

!  Number of fine layers for sscorr outgoing

      INTEGER  , INTENT(IN) :: NFINELAYERS

!  Los Control

      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      INTEGER  , INTENT(IN) :: N_USER_RELAZMS

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)

      INTEGER  , INTENT(IN) :: N_LOUTPUT

!  Layer masks for doing integrated source terms

      INTEGER  , INTENT(IN) :: LEVELMASK_UP (MAX_LOUTPUT)
      INTEGER  , INTENT(IN) :: LEVELMASK_DN (MAX_LOUTPUT)

!  Number of points, excitation index

      INTEGER  , INTENT(IN) :: NPOINTS_MONO, W_EXCIT

!  original solar angle for NADIR correction

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLE

!  NADIR-view scattering angle cosines

      REAL(FPK), INTENT(IN) :: COSSCAT_NADIR_UP ( MAX_GEOMETRIES )
      REAL(FPK), INTENT(IN) :: COSSCAT_NADIR_DN ( MAX_GEOMETRIES )

!  Solar beam input geometry

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLES_ADJUST ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  User VZA angles
!   @@@ RobFix 5/5/11, Wrong dimensioning. Fixed 5/5/11

!      REAL(FPK), INTENT(IN) :: USER_ANGLES_ADJUST ( MAX_USER_RELAZMS )
      REAL(FPK), INTENT(IN) :: USER_ANGLES_ADJUST ( MAX_USER_STREAMS )

!  User azimuth angles

      REAL(FPK), INTENT(IN) :: USER_RELAZMS_ADJUST ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Solar fluxes

      REAL(FPK), INTENT(IN) :: FLUXES_RANKED   ( MAX_POINTS )

!  Lambertian albedos. New Version 2.5, 9/11/15

      REAL(fpk), intent(in) :: ALBEDOS_RANKED ( MAX_POINTS )

!  Earth Radius

      REAL(FPK), INTENT(IN) :: EARTH_RADIUS

!  Height grid

      REAL(FPK), INTENT(IN) :: HEIGHT_GRID ( 0:MAX_LAYERS )

!  Exact (direct bounce) BRDF. Version 2.5 9/11/15

      REAL(fpk), intent(in) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

!  Isotropic Surface leaving term (if flag set). New Version 2.5, 9/11/15

      REAL(fpk), intent(in) :: SLTERM_ISOTROPIC ( MAX_POINTS )

!  Exact Surface-Leaving term. New Version 2.5, 9/11/15

      REAL(fpk), intent(in) :: SLTERM_USERANGLES (MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

!  These must be defined on the outer wavelength grid (binning)
!  Unscaled quantities, Elastic input

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: TRUNC_FACTORS     ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGA_UNSCALED    ( MAX_LAYERS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_UP  ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_ELASTIC_DN  ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_CABANNES_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAPHASFUNC_CABANNES_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSLOSS_UNSCALED ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSGAIN_UNSCALED ( MAX_LAYERS, 0:2, MAX_POINTS )

!  NADIR correction - user stream transmittances, whole layers

      REAL(FPK), INTENT(IN) :: SAVE_TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  NADIR correction - pre-calculated multipliers

      REAL(FPK), INTENT(IN) :: SAVE_BEAMMULT_UP ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: SAVE_BEAMMULT_DN ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Output
!  ------

!  Single scatter results

      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_UP (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_DN (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_UP   (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_DN   (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Exception handling

      LOGICAL          , INTENT(OUT) :: fail
      CHARACTER (LEN=*), INTENT(OUT) :: message

!  Geometry routine inputs and outputs
!  -----------------------------------

!  Control

      LOGICAL   :: do_fine
      REAL(FPK) :: alpha_boa, theta_boa, phi_boa

!  Main outputs (geometry)

      INTEGER   :: ntraverse(0:max_layers)
      REAL(FPK) :: sunpaths(0:max_layers,max_layers)
      REAL(FPK) :: radii   (0:max_layers)
      REAL(FPK) :: alpha_all  (0:max_layers)

!  Fine level output (geometry)

      INTEGER   :: ntraverse_fine(max_layers,max_fine_layers)
      REAL(FPK) :: sunpaths_fine (max_layers,max_layers,max_fine_layers)
      REAL(FPK) :: radii_fine    (max_layers,max_fine_layers)
      REAL(FPK) :: alpha_fine    (max_layers,max_fine_layers)

!  Other (incidental) geometrical output

      REAL(FPK) :: lospaths(max_layers)
      REAL(FPK) :: theta_all  (0:max_layers)
      REAL(FPK) :: phi_all    (0:max_layers)
      REAL(FPK) :: cosscat_up (0:max_layers)
      REAL(FPK) :: cosscat_dn (0:max_layers)

!  Extinction

      REAL(FPK) :: extinction ( max_layers, max_points )

!  Saved Legendre polynomials. 
!   -- Rob mod 5/12/17 for 2p5a, Only need 0:2 now, drop MAX_MOMENTS_INPUT

      REAL(FPK) :: SS_PLEG_UP(MAX_LAYERS,0:2), SS_PLEG_DN(MAX_LAYERS,0:2)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(FPK) :: TMS(MAX_LAYERS)

!  Local truncation factors for additional DELTAM scaling.  DISABLED
!      REAL(FPK) :: SSFDEL ( MAXLAYERS )

!  Exact Phase function calculations

      REAL(FPK) :: ESCAT_UP(MAX_LAYERS)
      REAL(FPK) :: ESCAT_DN(MAX_LAYERS)
      REAL(FPK) :: ISCAT_UP(MAX_LAYERS)
      REAL(FPK) :: ISCAT_DN(MAX_LAYERS)

!  Cumulative single scatter source terms

      REAL(FPK) :: ESS_CUMSOURCE_UP(0:MAX_LAYERS)
      REAL(FPK) :: ESS_CUMSOURCE_DN(0:MAX_LAYERS)
      REAL(FPK) :: ISS_CUMSOURCE_UP(0:MAX_LAYERS)
      REAL(FPK) :: ISS_CUMSOURCE_DN(0:MAX_LAYERS)

!  Outgoing sphericity stuff
!  -------------------------

!  Multipliers and LOS Transmittances

      REAL(FPK) :: &
         UP_EMULTIPLIERS(MAX_LAYERS), &
         DN_EMULTIPLIERS(MAX_LAYERS), &
         UP_IMULTIPLIERS(MAX_LAYERS,MAX_POINTS), &
         DN_IMULTIPLIERS(MAX_LAYERS,MAX_POINTS), &
         UP_LOSTRANS(MAX_LAYERS), &
         DN_LOSTRANS(MAX_LAYERS)

!  Solar beam attenuations to BOA (required for exact DB calculation)

      REAL(FPK) :: attn      ( 0:max_layers, max_points )
      REAL(FPK) :: attn_fine ( max_layers, max_fine_layers, max_points )

!  Local variables
!  ---------------

!  Indices

      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, W
      INTEGER   :: UTA, UM, IA, NC, V, WR
      INTEGER   :: Sleave_IDX, Brdf_IDX

!  Help variables (REAL(FPK) ::)
!   -- Rob mod 5/12/17 for 2p5a, removed DFL1, DFL2

      REAL(FPK) :: x0_boa, X0_BOA_4, boa_attn, omega_local, ssflux, ratio
      REAL(FPK) :: SOURCE, HELP, LSOURCE, COSSCAT, CONTRIB, REFLEC
      REAL(FPK) :: PHAS_E, PHAS_LR, PHAS_LG(MAX_POINTS), PHAS_CB
      REAL(FPK) :: FACTOR, CAB, GAIN, LOSS

!  Local surface flag

      LOGICAL, PARAMETER :: DO_INCLUDE_SURFACE = .true.

!  Set up operations
!  -----------------

!  Set flux scaling

      SSFLUX = one / PI4

!  Fine layering flag

      do_fine = .true.

!  Number of fine layers now a control input (7/25/08)

!  Floating point numbers for Legendre polynomials
!   -- Rob mod 5/12/17 for 2p5a, removed DFL1, DFL2
!      DO L = 2, NMOMENTS_INPUT
!        HELP = DBLE(L) ; DF1(L) = DBLE(2*L-1)/HELP ; DF2(L) = DBLE(L-1)/HELP
!      ENDDO

!  Create TMS factors, these get stored
!    Delta-M Scaling introduced April 2005.
!   -- Rob mod 5/12/17 for 2p5a, Use OMEGA_UNSCALED

      DO N = 1, NLAYERS
        TMS(N) = ONE
        IF ( DO_DMSCALING ) THEN
          OMEGA_LOCAL = OMEGA_UNSCALED (N,W_EXCIT )
          HELP = TRUNC_FACTORS(N,W_EXCIT) * OMEGA_LOCAL
          TMS(N) = TMS(N) / (ONE - HELP)
        ENDIF
      ENDDO

!  Create extinctions

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        WR = 0
        DO W = 1, NPOINTS_MONO
          EXTINCTION(N,W) = DELTAU_VERT_INPUT(N,W) / HELP
        ENDDO
      ENDDO

!  Additional Delta-M scaling. DISABLED
!  ------------------------------------

!  New section. R. Spurr, 07 September 2007.
!   TMS gets modified by (1-F). Save the truncation factor.
!   Phase function moments are modified later on.
!      IF ( DO_SSCORR_TRUNCATION ) THEN
!        NM1  = NMOMENTS_INPUT
!        DNM1 = DBLE(2*NM1+1)
!        DO N = 1, NLAYERS
!          BETA = OMEGAMOMS_ELASTIC  (N, NM1, W_EXCIT ) / OMEGA_LOCAL(N)
!          SSFDEL(N) = BETA / DNM1
!          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
!        ENDDO
!      ENDIF

!  Start the main loop over all viewing geometries
!  -----------------------------------------------

!    Outgoing correction, use the adjusted values of the angles -------> March 2011
!   -- Rob mod 5/12/17 for 2p5a, re-introduce the SSCORR_NADIR computation

      DO UM = 1, N_USER_STREAMS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  UPWELLING NADIR CORRECTION: multipliers, transmittances
!  =======================================================

!  Upwelling Multipliers and Transmittances for the SSCORR_NADIR case
!     Does not use the ADJUSTED geometrical angles, Same for each azimuth.
!   -- Rob mod 5/12/17 for 2p5a, SSCORR_NADIR values are already pre-calculated

        IF ( DO_SSCORR_NADIR .and.DO_UPWELLING ) THEN
          DO N = 1, NLAYERS
            UP_LOSTRANS(N) = SAVE_TRANS_USERM(UM,N,W_EXCIT)
            UP_EMULTIPLIERS(N) = SAVE_BEAMMULT_UP(UM,N,W_EXCIT)
            DO W = 1, NPOINTS_MONO
              UP_IMULTIPLIERS(N,W) = SAVE_BEAMMULT_UP(UM,N,W)
            ENDDO
          ENDDO
        ENDIF

!  DOWNWELLING NADIR CORRECTION: multipliers, transmittances
!  =========================================================

!  Downwelling Multipliers and Transmittances for the SSCORR_NADIR case
!     Does not use the ADJUSTED geometrical angles, Same for each azimuth.
!   -- Rob mod 5/12/17 for 2p5a, SSCORR_NADIR values are already pre-calculated

        IF ( DO_SSCORR_NADIR .and.DO_DNWELLING ) THEN
          DO N = 1, NLAYERS
            DN_LOSTRANS(N) = SAVE_TRANS_USERM(UM,N,W_EXCIT)
            DN_EMULTIPLIERS(N) = SAVE_BEAMMULT_DN(UM,N,W_EXCIT)
            DO W = 1, NPOINTS_MONO
              DN_IMULTIPLIERS(N,W) = SAVE_BEAMMULT_DN(UM,N,W)
            ENDDO
          ENDDO
        ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  Start the Azimuth loop
!  ----------------------

        DO IA = 1, N_USER_RELAZMS

!  Adjusted geometry for the Outgoing correction

          ALPHA_BOA = USER_ANGLES_ADJUST(UM)
          THETA_BOA = SOLAR_ANGLES_ADJUST(UM,IA)
          PHI_BOA   = USER_RELAZMS_ADJUST(UM,IA)
          V = N_USER_RELAZMS * (UM-1) + IA

!  UPWELLING OUTGOING CORRECTION: multipliers, transmittances
!  ==========================================================

          IF ( DO_UPWELLING .and. DO_SSCORR_OUTGOING ) THEN

!  Call to geometry routine, path distances etc....

            call outgoing_sphergeom_fine_up_1 &
        ( max_layers, max_fine_layers, do_fine, nlayers, nfinelayers, & ! Inputs
          height_grid, earth_radius, alpha_boa, theta_boa, phi_boa,   & ! Inputs
          sunpaths,      radii,      ntraverse,      alpha_all,       & ! Inputs
          sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine,      & ! Inputs
          lospaths, theta_all, phi_all, cosscat_up,                   & ! Inputs
          fail, message )                                               ! Outputs

!  Exception handling

            if ( fail ) return

!  Get the attenuations

            call outgoing_attenuations_1 &
        ( max_layers, max_fine_layers, max_points, & ! Inputs
          npoints_mono, nlayers, nfinelayers,      & ! Inputs
          extinction, sunpaths, ntraverse,         & ! Inputs
          sunpaths_fine, ntraverse_fine,           & ! Inputs
          attn, attn_fine )                          ! Outputs

!  Multipliers, transmittances

            call outgoing_integration_mono_up_1 &
        ( max_layers, max_fine_layers,                   & ! Inputs
          max_points, npoints_mono, w_excit,             & ! Inputs
          nlayers, nfinelayers, extinction,              & ! Inputs
          radii, alpha_all, alpha_fine, attn, attn_fine, & ! Inputs
          up_emultipliers, up_imultipliers, up_lostrans )  ! Outputs

!   End Outgoing correction multipliers etc.

          endif

!  UPWELLING SINGLE SCATTER CORRECTION:  Source terms
!  ==================================================

          IF ( DO_UPWELLING ) THEN

!   Cosine scatter angles
!   -- Rob mod 5/12/17 for 2p5a, Careful of selection.

            if ( DO_SSCORR_NADIR ) THEN
              COSSCAT = COSSCAT_NADIR_UP(V)
            else
              COSSCAT = COSSCAT_UP(NLAYERS)
            endif

!  legendre polynomials
!   -- Rob mod 5/12/17 for 2p5a, Only need #2 = (3x^2-1)/2

            SS_PLEG_UP(1,0) = ONE
            SS_PLEG_UP(1,1) = COSSCAT
            SS_PLEG_UP(1,2) = 1.5_fpk * COSSCAT * COSSCAT - 0.5_fpk
!            DO L = 2, NMOMENTS_INPUT
!               SS_PLEG_UP(1,L) = DF1(L) * SS_PLEG_UP(1,L-1) * COSSCAT - DF2(L) * SS_PLEG_UP(1,L-2)
!            ENDDO

!  Elastic scattering Phase functions (multiplied by TMS factor)
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_E with Phase function input
!   -- PHAS_E = DOT_PRODUCT(OMEGAMOMS_ELASTIC_UNSCALED(N,0:NMOMENTS_INPUT,W_EXCIT),SS_PLEG_UP(1,0:NMOMENTS_INPUT)

            DO N = 1, NLAYERS
              PHAS_E = OMEGAPHASFUNC_ELASTIC_UP(N,V,W_EXCIT) * TMS(N)
              ESCAT_UP(N) = PHAS_E * UP_EMULTIPLIERS(N)
            ENDDO

!  Add inelastic scattering Phase functions if flagged

            IF ( DO_FULL_RAMAN ) THEN

!  Energy Balancing method
!  -----------------------

             IF ( DO_ENERGY_BALANCING ) THEN

!  start layer loop

              DO N = 1, NLAYERS

!  Loss term

                PHAS_LR = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,1) + &
                          OMEGAMOMS_RRSLOSS_UNSCALED(N,2,1) * SS_PLEG_UP(1,2)
                LOSS = PHAS_LR * UP_EMULTIPLIERS(N)
! @@@ Robfix 06/01/11. Add TMS
                LOSS = LOSS * TMS(N)

!  Gain term

                GAIN = ZERO
                DO W = 1, NPOINTS_MONO
                  RATIO   = FLUXES_RANKED(W)/FLUXES_RANKED(W_EXCIT)
                  CONTRIB = OMEGAMOMS_RRSGAIN_UNSCALED(N,0,W) + &
                            OMEGAMOMS_RRSGAIN_UNSCALED(N,2,W) * SS_PLEG_UP(1,2)
                  PHAS_LG(W) = CONTRIB * RATIO
                  GAIN = GAIN + PHAS_LG(W) * UP_IMULTIPLIERS(N,W)
                ENDDO
! @@@ Robfix 06/01/11. Add TMS
                GAIN = GAIN * TMS(N)

!  Whole layer contributions

                ISCAT_UP(N) = ESCAT_UP(N) + GAIN - LOSS

!  End loop

              ENDDO

!  Cabannes Raman method
!  ---------------------

             ELSE IF ( DO_CABANNES_RAMAN ) THEN

!  Start loop

              DO N = 1, NLAYERS

!  Cabannes term
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_CB with Phase function input
!   -- PHAS_CB = DOT_PRODUCT(OMEGAMOMS_CABANNES_UNSCALED(N,0:NMOMENTS_INPUT,W),SS_PLEG_UP(1,0:NMOMENTS_INPUT)

                PHAS_CB = OMEGAPHASFUNC_CABANNES_UP(N,V,1) * TMS(N)
                CAB = PHAS_CB * UP_EMULTIPLIERS(N)

!  Gain term

                GAIN = ZERO
                DO W = 1, NPOINTS_MONO
                  RATIO   = FLUXES_RANKED(W)/FLUXES_RANKED(W_EXCIT)
                  CONTRIB = OMEGAMOMS_RRSGAIN_UNSCALED(N,0,W) + &
                            OMEGAMOMS_RRSGAIN_UNSCALED(N,2,W) * SS_PLEG_UP(1,2)
                  PHAS_LG(W)  = CONTRIB * RATIO
                  GAIN = GAIN + PHAS_LG(W) * UP_IMULTIPLIERS(N,W)
                ENDDO
                GAIN = GAIN * TMS(N)

!  Whole layer contribution

                ISCAT_UP(N) = CAB + GAIN

!  End loop

              ENDDO

!  End inelastic scattering options

             ENDIF
            ENDIF

!  End upwelling clause

          ENDIF

!  DOWNWELLING OUTGOING CORRECTION: multipliers, transmittances
!  ============================================================

          IF ( DO_DNWELLING .and. DO_SSCORR_OUTGOING ) THEN

!  Call to geometry routine, path distances etc....

            call outgoing_sphergeom_fine_dn_1 &
        ( max_layers, max_fine_layers, do_fine, nlayers, nfinelayers, & ! Inputs
          height_grid, earth_radius, alpha_boa, theta_boa, phi_boa,   & ! Inputs
          sunpaths,      radii,      ntraverse,      alpha_all,       & ! Inputs
          sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine,      & ! Inputs
          lospaths, theta_all, phi_all, cosscat_dn,                   & ! Inputs
          fail, message )                                               ! Outputs

!  Exception handling

            if ( fail ) return

!  Get the attenuations

            call outgoing_attenuations_1 &
        ( max_layers, max_fine_layers, max_points, & ! Inputs
          npoints_mono, nlayers, nfinelayers,      & ! Inputs
          extinction, sunpaths, ntraverse,         & ! Inputs
          sunpaths_fine, ntraverse_fine,           & ! Inputs
          attn, attn_fine )                          ! Outputs

!  Multipliers and transmittances

            call outgoing_integration_mono_dn_1 &
        ( max_layers, max_fine_layers,                   & ! Inputs
          max_points, npoints_mono, w_excit,             & ! Inputs
          nlayers, nfinelayers, extinction,              & ! Inputs
          radii, alpha_all, alpha_fine, attn, attn_fine, & ! Inputs
          dn_emultipliers, dn_imultipliers, dn_lostrans )  ! Outputs

!   End Outgoing correction multipliers etc.

          endif

!  DOWNWELLING SINGLE SCATTER CORRECTION:  Source terms
!  ====================================================

          IF ( DO_DNWELLING ) THEN

!   Cosine scatter angles
!   -- Rob mod 5/12/17 for 2p5a, Careful of selection.

            if ( DO_SSCORR_NADIR ) THEN
              COSSCAT = COSSCAT_NADIR_DN(V)
            else
              COSSCAT = COSSCAT_DN(NLAYERS)
            endif

!  Legendre polynomials
!   -- Rob mod 5/12/17 for 2p5a, Only need #2 = (3x^2-1)/2

            SS_PLEG_DN(1,0) = ONE
            SS_PLEG_DN(1,1) = COSSCAT
            SS_PLEG_DN(1,2) = 1.5_fpk * COSSCAT * COSSCAT - 0.5_fpk
!            DO L = 2, NMOMENTS_INPUT
!               SS_PLEG_DN(1,L) = DF1(L) * SS_PLEG_DN(1,L-1) * COSSCAT - DF2(L) * SS_PLEG_DN(1,L-2)
!            ENDDO

!  Elastic scattering Phase functions (multiplied by TMS factor)
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_E with Phase function input
!   -- PHAS_E = DOT_PRODUCT(OMEGAMOMS_ELASTIC_UNSCALED(N,0:NMOMENTS_INPUT,W_EXCIT),SS_PLEG_DN(1,0:NMOMENTS_INPUT)

            DO N = 1, NLAYERS
              PHAS_E = OMEGAPHASFUNC_ELASTIC_DN(N,V,W_EXCIT) * TMS(N)
              ESCAT_DN(N) = PHAS_E * DN_EMULTIPLIERS(N)
            ENDDO

!  Add inelastic scattering Phase functions if flagged

            IF ( DO_FULL_RAMAN ) THEN

!  Energy balancing method
!  -----------------------

             IF ( DO_ENERGY_BALANCING ) THEN

!  Start loop

              DO N = 1, NLAYERS

!  Loss term

                PHAS_LR = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,1) + &
                          OMEGAMOMS_RRSLOSS_UNSCALED(N,2,1) * SS_PLEG_DN(1,2)
                LOSS = PHAS_LR * DN_EMULTIPLIERS(N)
                LOSS = LOSS * TMS(N)

!  Gain term

                GAIN = ZERO
                DO W = 1, NPOINTS_MONO
                  RATIO   = FLUXES_RANKED(W)/FLUXES_RANKED(W_EXCIT)
                  CONTRIB = OMEGAMOMS_RRSGAIN_UNSCALED(N,0,W) + &
                            OMEGAMOMS_RRSGAIN_UNSCALED(N,2,W) * SS_PLEG_DN(1,2)
                  PHAS_LG(W) = CONTRIB * RATIO
                  GAIN = GAIN + PHAS_LG(W) * DN_IMULTIPLIERS(N,W)
                ENDDO
                GAIN = GAIN * TMS(N)

!  Whole layer contributions

                ISCAT_DN(N) = ESCAT_DN(N) + GAIN - LOSS

!  End loop

              ENDDO

!  Cabannes Raman method
!  ---------------------

             ELSE IF ( DO_CABANNES_RAMAN ) THEN

!  Start loop

              DO N = 1, NLAYERS

!  Cabannes term
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_CB with Phase function input
!   -- PHAS_CB = DOT_PRODUCT(OMEGAMOMS_CABANNES_UNSCALED(N,0:NMOMENTS_INPUT,W),SS_PLEG_UP(1,0:NMOMENTS_INPUT)

                PHAS_CB = OMEGAPHASFUNC_CABANNES_DN(N,V,1) * TMS(N)
                CAB = PHAS_CB * DN_EMULTIPLIERS(N)

!  Gain term

                GAIN = ZERO
                DO W = 1, NPOINTS_MONO
                  RATIO   = FLUXES_RANKED(W)/FLUXES_RANKED(W_EXCIT)
                  CONTRIB = OMEGAMOMS_RRSGAIN_UNSCALED(N,0,W) + &
                            OMEGAMOMS_RRSGAIN_UNSCALED(N,2,W) * SS_PLEG_DN(1,2)
                  PHAS_LG(W)  = CONTRIB * RATIO
                  GAIN = GAIN + PHAS_LG(W) * DN_IMULTIPLIERS(N,W)
                ENDDO
                GAIN = GAIN * TMS(N)

!  Whole layer contribution

                ISCAT_DN(N) = CAB + GAIN

!  End loop

              ENDDO

!  End inelastic scattering options

             ENDIF
            ENDIF

!  End Downwelling clause for source terms

          ENDIF

!  Recurrence relation for the UPWELLING intensity
!  ===============================================

          IF ( DO_UPWELLING ) THEN

!   This section of code revised for the Version 2.5 supplement-derived BRDF.

            IF ( DO_INCLUDE_SURFACE ) THEN
              IF ( DO_BRDF_SURFACE ) THEN
                Brdf_IDX = 1 ; if ( DO_BRDF_Wav1 ) Brdf_IDX = W_EXCIT
                REFLEC = EXACTDB_BRDFUNC(UM,IA,Brdf_IDX)
              ELSE
                REFLEC = ALBEDOS_RANKED(W_EXCIT)
              ENDIF
            ENDIF

!  initialize cumulative source term = F.A.mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux, A = albedo
!    Only require the Stokes total intensity component
!    This contribution is purely elastic.
!    This line Replaced: X0_BOA = cos(SOLAR_ANGLE * DEG_TO_RAD). reintroduced for Version 2.5a
!     Introduction of the Adjusted Gometry value, 18 March 2011

            NC =  0
            IF ( DO_INCLUDE_SURFACE ) THEN
              IF ( DO_SSCORR_NADIR ) then
                X0_BOA = cos(SOLAR_ANGLE * DEG_TO_RAD)
              ELSE
                X0_BOA = cos(SOLAR_ANGLES_ADJUST(UM,IA) * DEG_TO_RAD)
              ENDIF
              X0_BOA_4 = FOUR * X0_BOA
              FACTOR = REFLEC
              BOA_ATTN = X0_BOA_4 * ATTN(NLAYERS,W_EXCIT)
              ESS_CUMSOURCE_UP(NC) = FACTOR * BOA_ATTN
            ELSE
              ESS_CUMSOURCE_UP(NC) = ZERO
            ENDIF

!  Version 2.4 code commented out................... 9/11/15
!  Add SLTERM flag     . @@@ Rob Fix 06 Sep 12.
!  New Isotropic SLTERM. @@@ Rob Fix 06 sep 12
!   Multiply term by 4pi to get correct normalization. 7/30/12
!            IF ( DO_LRRS_SLTERM ) THEN
!                ESS_CUMSOURCE_UP(NC) = ESS_CUMSOURCE_UP(NC) + LRRS_SLTERM(W_EXCIT)*PI4
!            ENDIF

!  New Surface-Leaving stuff for Version 2.5.. 9/11/15
!    Multiply terms by 4pi to get correct normalization. 7/30/12

            IF ( DO_SURFACE_LEAVING ) THEN
              Sleave_IDX = 1 ; if ( DO_SLEAVE_Wav1 ) Sleave_IDX = W_EXCIT
              IF ( DO_SL_ISOTROPIC ) THEN
                ESS_CUMSOURCE_UP(NC) = ESS_CUMSOURCE_UP(NC) + SLTERM_ISOTROPIC(Sleave_IDX)*PI4
              ELSE
                ESS_CUMSOURCE_UP(NC) = ESS_CUMSOURCE_UP(NC) + SLTERM_USERANGLES(UM,IA,Sleave_IDX)*PI4
              ENDIF
            ENDIF

!  Copy for the Raman

            IF ( DO_FULL_RAMAN ) THEN
              ISS_CUMSOURCE_UP(NC) = ESS_CUMSOURCE_UP(NC)
            ENDIF

!  initialize optical depth loop

            NSTART = NLAYERS
            NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

            DO UTA = N_LOUTPUT, 1, -1

!  Layer index for given optical depth

              NLEVEL = LEVELMASK_UP(UTA)
              NUT    = NLEVEL + 1

!  Cumulative single scatter source terms :
!      For loop over layers working upwards to level NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier
!  Multiplier using new integration scheme

              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N
                LSOURCE =  ESCAT_UP(N)
                ESS_CUMSOURCE_UP(NC) = ESCAT_UP(N)+ &
                       UP_LOSTRANS(N) * ESS_CUMSOURCE_UP(NC-1)
                IF ( DO_FULL_RAMAN ) THEN
                  ISS_CUMSOURCE_UP(NC) = ISCAT_UP(N) + &
                       UP_LOSTRANS(N) * ISS_CUMSOURCE_UP(NC-1)
                ENDIF
              ENDDO

!  Ongrid output--------
!     Set final cumulative source and single scatter intensity

              SOURCE = ESS_CUMSOURCE_UP(NC)
              ELASTIC_SS_UP(UTA,V,1) = SSFLUX * SOURCE
              IF ( DO_FULL_RAMAN ) THEN
                SOURCE = ISS_CUMSOURCE_UP(NC)
                RAMAN_SS_UP(UTA,V,1) = SSFLUX * SOURCE
              ENDIF

!  Check for updating the recursion

              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
              NUT_PREV = NUT

!  end optical depth loop and Upwelling clause

            ENDDO
          ENDIF

!  Recurrence relation for the DOWNWELLING intensity
!  =================================================

          IF ( DO_DNWELLING ) THEN

!  initialize cumulative source term

            NC =  0
            ESS_CUMSOURCE_DN(NC) = ZERO
            IF ( DO_FULL_RAMAN ) THEN
              ISS_CUMSOURCE_DN(NC) = ESS_CUMSOURCE_DN(NC)
            ENDIF

!  initialise optical depth loop

            NSTART = 1
            NUT_PREV = NSTART - 1

!  Main loop over all output optical depths

            DO UTA = 1, N_LOUTPUT

!  Layer index for given optical depth

              NLEVEL = LEVELMASK_DN(UTA)
              NUT = NLEVEL

!  Cumulative single scatter source terms :
!      For loop over layers working downwards to NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier
!      Multiplier by new integration method

              DO N = NSTART, NUT
                NC = N
                LSOURCE =  ESCAT_DN(N)
                ESS_CUMSOURCE_DN(NC) = ESCAT_DN(N)+ &
                       DN_LOSTRANS(N) * ESS_CUMSOURCE_DN(NC-1)
                IF ( DO_FULL_RAMAN ) THEN
                  ISS_CUMSOURCE_DN(NC) = ISCAT_DN(N) + &
                       DN_LOSTRANS(N) * ISS_CUMSOURCE_DN(NC-1)
                ENDIF
              ENDDO

!  Ongrid output :
!     Set final cumulative source and correct Stokes vector

              SOURCE = ESS_CUMSOURCE_DN(NC)
              ELASTIC_SS_DN(UTA,V,1) = SSFLUX * SOURCE
              IF ( DO_FULL_RAMAN ) THEN
                SOURCE = ISS_CUMSOURCE_DN(NC)
                RAMAN_SS_DN(UTA,V,1) = SSFLUX * SOURCE
              ENDIF

!  Check for updating the recursion

              IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
              NUT_PREV = NUT

!  end optical depth loop and Downwelling clause

            ENDDO
          ENDIF

!  Finish geometry loop

        enddo
      enddo

!  Finish

      RETURN
      END SUBROUTINE LRRS_SSCORR_GENERAL_MONO_1

!

      SUBROUTINE OUTGOING_ATTENUATIONS_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS,  & ! Inputs
          N_RAMAN_POINTS, NLAYERS, NFINELAYERS, & ! Inputs
          EXTINCTION, SUNPATHS, NTRAVERSE,      & ! Inputs
          SUNPATHS_FINE, NTRAVERSE_FINE,        & ! Inputs
          ATTN, ATTN_FINE )                       ! Outputs

      USE LRRS_PARS_m, Only : FPK, ZERO

!  Does attenuations, Whole layers only. No partials.

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER  , INTENT(IN) :: maxlayers, maxfinelayers, maxpoints

!  control

      INTEGER  , INTENT(IN) :: nfinelayers, nlayers, n_raman_points

!  Whole layers

      INTEGER  , INTENT(IN) :: ntraverse(0:maxlayers)
      REAL(FPK), INTENT(IN) :: sunpaths(0:maxlayers,maxlayers)

!  Fine level

      INTEGER  , INTENT(IN) :: ntraverse_fine(maxlayers,maxfinelayers)
      REAL(FPK), INTENT(IN) :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)

!  Extinction

      REAL(FPK), INTENT(IN) :: extinction (maxlayers, maxpoints)

!  outputs
!  -------

      REAL(FPK), INTENT(OUT) :: attn      ( 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(OUT) :: attn_fine ( maxlayers, maxfinelayers, maxpoints )

!  help variables
!  --------------

      INTEGER   :: w, n, j, k
      REAL(FPK) :: tau

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0_FPK

!  attenuation functions, whole layers

      do w = 1, n_raman_points
       do n = 0, nlayers
        tau     = zero
        attn(n,w) = zero
        do k = 1, ntraverse(n)
         tau = tau + sunpaths(n,k) * extinction(k,w)
        enddo
        if ( tau .le. local_cutoff ) attn(n,w) = exp(-tau)
       enddo
      enddo

!  Fine grid attenuations

      do n = 1, nlayers
       do j = 1, nfinelayers
        do w = 1, n_raman_points
         tau            = zero
         attn_fine(n,j,w) = zero
         do k = 1, ntraverse_fine(n,j)
          tau = tau + sunpaths_fine(n,k,j) * extinction(k,w)
         enddo
         if ( tau .le. local_cutoff ) attn_fine(n,j,w) = exp(-tau)
        enddo
       enddo
      enddo

!  Finish

      return
      END SUBROUTINE OUTGOING_ATTENUATIONS_1

!

      SUBROUTINE OUTGOING_INTEGRATION_BIN_UP_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXBINS,  & ! Inputs
          NPOINTS_INNER, OFFSET_INNER, NBINS, BINMAP,    & ! Inputs
          NLAYERS, NFINELAYERS, EXTINCTION,              & ! Inputs
          RADII, ALPHA_ALL, ALPHA_FINE, ATTN, ATTN_FINE, & ! Inputs
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS )           ! Outputs

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, HALF

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER  , INTENT(IN) :: maxlayers, maxfinelayers
      INTEGER  , INTENT(IN) :: maxpoints, maxbins

!  control

      INTEGER  , INTENT(IN) :: nfinelayers, nlayers
      INTEGER  , INTENT(IN) :: npoints_inner, offset_inner
      INTEGER  , INTENT(IN) :: nbins  ( maxpoints )
      INTEGER  , INTENT(IN) :: binmap ( maxbins, maxpoints )

!  radii

      REAL(FPK), INTENT(IN) :: radii   (0:maxlayers)

!  line of sight angles

      REAL(FPK), INTENT(IN) :: alpha_all  (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_fine    (maxlayers,maxfinelayers)

!  Extinction

      REAL(FPK), INTENT(IN) :: extinction (maxlayers, maxpoints)

!  Attenuations

      REAL(FPK), INTENT(IN) :: attn      ( 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: attn_fine ( maxlayers, maxfinelayers, maxpoints )

!  outputs
!  -------

      REAL(FPK), INTENT(OUT) :: imultipliers   (maxlayers,   maxbins, maxpoints)
      REAL(FPK), INTENT(OUT) :: emultipliers   (maxlayers,   maxpoints)
      REAL(FPK), INTENT(OUT) :: lostrans      (maxlayers,   maxpoints)

!  local arrays
!  ------------

!  Local geoemetry arrays

      REAL(FPK) :: csq_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: cot_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: tranj     ( maxfinelayers)

!  help variables
!  --------------

      INTEGER   :: w, b, wb, wr, n, j
      REAL(FPK) :: sumr, salpha, calpha, dfine1, raycon, kn
      REAL(FPK) :: csq_1, csqt_1, csq_2, cot_1, cot_2
      REAL(FPK) :: func_1, func_2, tran_1, step, func, skn

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0_FPK

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        do w = 1, npoints_inner
          lostrans(n,w)    = zero
          emultipliers(n,w) = zero
          do b = 1, nbins(w)
            imultipliers(n,b,w) = zero
          enddo
        enddo
      enddo

!  Ray constant at TOA

      salpha = sin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work up from the bottom of the atmosphere
!  =========================================

!  initialise

      n = nlayers
      salpha = sin(alpha_all(n))
      calpha = cos(alpha_all(n))
      csq_1 = one / salpha / salpha
      cot_1 = calpha / salpha

!  save some quantities

      do n = nlayers, 1, -1
        do j = 1, nfinelayers
          calpha = cos(alpha_fine(n,j))
          salpha = sin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = one / salpha / salpha
        enddo
      enddo

!  Start layer loop

      do n = nlayers, 1, -1

!  Save some quantities

        salpha = sin(alpha_all(n-1))
        calpha = cos(alpha_all(n-1))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  Trapezium-rule integrated sourceterm multipliers + transmittances
!  -----------------------------------------------------------------

        do w = 1, npoints_inner

!  set up

          wr = w + offset_inner
          kn = raycon * extinction(n,wr)
          skn = step * kn
          tran_1 = exp ( - kn * ( cot_2 - cot_1 ) )
          csqt_1 = csq_1 * tran_1
          do j = 1, nfinelayers
            tranj(j) = exp ( - kn * ( cot_2 - cot_fine(n,j) ) )
          enddo

!  Line of sight transmittance factor

          lostrans(n,w)    = tran_1

!  Elastic scattering multiplier

          func_2 = attn(n-1,wr) *  csq_2
          func_1 = attn(n,wr)   * csqt_1
          sumr = half * ( func_1 + func_2 )
          do j = 1, nfinelayers
            func = attn_fine(n,j,wr) * tranj(j) * csq_fine(n,j)
            sumr  = sumr + func
          enddo
          emultipliers(n,w) = sumr * skn

!  Inelastic scattering multipliers (by bins)

          do b = 1, nbins(w)
            wb = binmap(b,w)
            func_2 = attn(n-1,wb) *  csq_2
            func_1 = attn(n,wb)   * csqt_1
            sumr = half * ( func_1 + func_2 )
            do j = 1, nfinelayers
              func = attn_fine(n,j,wb) * tranj(j) * csq_fine(n,j)
              sumr  = sumr + func
            enddo
            imultipliers(n,b,w) = sumr * skn
          enddo

!  End points loop

        enddo

!  update geometry

        csq_1 = csq_2
        cot_1 = cot_2

!  Finish layer

      ENDDO

!  finish

      RETURN
      END SUBROUTINE OUTGOING_INTEGRATION_BIN_UP_1

!

      SUBROUTINE OUTGOING_INTEGRATION_BIN_DN_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXBINS,  & ! Inputs
          NPOINTS_INNER, OFFSET_INNER, NBINS, BINMAP,    & ! Inputs
          NLAYERS, NFINELAYERS, EXTINCTION,              & ! Inputs
          RADII, ALPHA_ALL, ALPHA_FINE, ATTN, ATTN_FINE, & ! Inputs
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS )           ! Outputs

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, HALF

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER  , INTENT(IN) :: maxlayers, maxfinelayers
      INTEGER  , INTENT(IN) :: maxpoints, maxbins

!  control

      INTEGER  , INTENT(IN) :: nfinelayers, nlayers
      INTEGER  , INTENT(IN) :: npoints_inner, offset_inner
      INTEGER  , INTENT(IN) :: nbins  ( maxpoints )
      INTEGER  , INTENT(IN) :: binmap ( maxbins, maxpoints )

!  radii

      REAL(FPK), INTENT(IN) :: radii   (0:maxlayers)

!  line of sight angles

      REAL(FPK), INTENT(IN) :: alpha_all  (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_fine (maxlayers,maxfinelayers)

!  Extinction

      REAL(FPK), INTENT(IN) :: extinction (maxlayers, maxpoints)

!  Attenuations

      REAL(FPK), INTENT(IN) :: attn      ( 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: attn_fine &
               ( maxlayers, maxfinelayers, maxpoints )

!  outputs
!  -------

      REAL(FPK), INTENT(OUT) :: imultipliers   (maxlayers,   maxbins, maxpoints)
      REAL(FPK), INTENT(OUT) :: emultipliers   (maxlayers,   maxpoints)
      REAL(FPK), INTENT(OUT) :: lostrans      (maxlayers,   maxpoints)

!  local arrays
!  ------------

!  Local geometry arrays

      REAL(FPK) :: csq_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: cot_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: tranj     ( maxfinelayers)

!  help variables
!  --------------

      INTEGER   :: w, wb, b, wr, n, j
      REAL(FPK) :: sumr, salpha, calpha, dfine1, raycon, kn
      REAL(FPK) :: csq_1, csqt_1, csq_2, cot_1, cot_2
      REAL(FPK) :: func_1, func_2, tran_1, step, func, skn

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0_FPK

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        do w = 1, npoints_inner
          lostrans(n,w)    = zero
          emultipliers(n,w) = zero
          do b = 1, nbins(w)
            imultipliers(n,b,w) = zero
          enddo
        enddo
      enddo

!  Ray constant at TOA

      salpha = sin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work Down from the top of the atmosphere
!  ========================================

!  initialise

      n = 0
      salpha = sin(alpha_all(n))
      calpha = cos(alpha_all(n))
      csq_1 = one / salpha / salpha
      cot_1 = calpha / salpha

!  Save some quantities

      do n = 1, nlayers
        do j = nfinelayers, 1, -1
          calpha = cos(alpha_fine(n,j))
          salpha = sin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = one / salpha / salpha
        enddo
      enddo

!  Start layer loop

      do n = 1, nlayers

!  Save some quantities

        salpha = sin(alpha_all(n))
        calpha = cos(alpha_all(n))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  Trapezium-rule integrated sourceterm multiplier + transmittance

        do w = 1, npoints_inner

!  set up

          wr = w + offset_inner
          kn = raycon * extinction(n,wr)
          skn = step * kn
          tran_1 = exp ( - kn * ( cot_1 - cot_2 ) )
          csqt_1 = csq_1 * tran_1
          do j = 1, nfinelayers
            tranj(j) = exp ( - kn * ( cot_fine(n,j) - cot_2 ) )
          enddo

!  Line of sight transmittance factor

          lostrans(n,w)    = tran_1

!  Elastic multipliers

          func_1 = attn(n-1,wr) * csqt_1
          func_2 = attn(n,wr)   * csq_2
          sumr = half * ( func_1 + func_2 )
          do j = nfinelayers, 1, -1
            func = attn_fine(n,j,wr) * tranj(j) * csq_fine(n,j)
            sumr = sumr + func
          enddo
          emultipliers(n,w) = sumr * step * kn

!  Inelastic multipliers

          do b = 1, nbins(w)
            wb = binmap(b,w)
            func_1 = attn(n-1,wb) * csqt_1
            func_2 = attn(n,wb)   * csq_2
            sumr = half * ( func_1 + func_2 )
            do j = nfinelayers, 1, -1
              func = attn_fine(n,j,wb) * tranj(j) * csq_fine(n,j)
              sumr = sumr + func
            enddo
            imultipliers(n,b,w) = sumr * step * kn
          enddo

!  End points loop

        enddo

!  update geometry

        csq_1 = csq_2
        cot_1 = cot_2

!  Finish layer

      ENDDO

!  finish

      RETURN
      END SUBROUTINE OUTGOING_INTEGRATION_BIN_DN_1

!

      SUBROUTINE OUTGOING_INTEGRATION_MONO_UP_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS,           & ! Inputs
          NPOINTS_MONO, W_EXCIT,                         & ! Inputs
          NLAYERS, NFINELAYERS, EXTINCTION,              & ! Inputs
          RADII, ALPHA_ALL, ALPHA_FINE, ATTN, ATTN_FINE, & ! Inputs
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS )           ! Outputs

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, HALF

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER  , INTENT(IN) :: maxlayers, maxfinelayers
      INTEGER  , INTENT(IN) :: maxpoints

!  control

      INTEGER  , INTENT(IN) :: nfinelayers, nlayers
      INTEGER  , INTENT(IN) :: npoints_mono, w_excit

!  radii

      REAL(FPK), INTENT(IN) :: radii   (0:maxlayers)

!  line of sight angles

      REAL(FPK), INTENT(IN) :: alpha_all  (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_fine    (maxlayers,maxfinelayers)

!  Extinction

      REAL(FPK), INTENT(IN) :: extinction (maxlayers, maxpoints)

!  Attenuations

      REAL(FPK), INTENT(IN) :: attn      ( 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: attn_fine &
               ( maxlayers, maxfinelayers, maxpoints )

!  outputs
!  -------

      REAL(FPK), INTENT(OUT) :: imultipliers   (maxlayers,   maxpoints)
      REAL(FPK), INTENT(OUT) :: emultipliers   (maxlayers)
      REAL(FPK), INTENT(OUT) :: lostrans      (maxlayers)

!  local arrays
!  ------------

!  Local geoemetry arrays

      REAL(FPK) :: csq_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: cot_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: tranj     ( maxfinelayers)

!  help variables
!  --------------

      INTEGER   :: w, n, j
      REAL(FPK) :: sumr, salpha, calpha, dfine1, raycon, kn
      REAL(FPK) :: csq_1, csqt_1, csq_2, cot_1, cot_2
      REAL(FPK) :: func_1, func_2, tran_1, step, func, skn

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0_FPK

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        lostrans(n)     = zero
        emultipliers(n) = zero
        do w = 1, npoints_mono
          imultipliers(n,w) = zero
        enddo
      enddo

!  Ray constant at TOA

      salpha = sin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work up from the bottom of the atmosphere
!  =========================================

!  initialise

      n = nlayers
      salpha = sin(alpha_all(n))
      calpha = cos(alpha_all(n))
      csq_1 = one / salpha / salpha
      cot_1 = calpha / salpha

!  save some quantities

      do n = nlayers, 1, -1
        do j = 1, nfinelayers
          calpha = cos(alpha_fine(n,j))
          salpha = sin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = one / salpha / salpha
        enddo
      enddo

!  Start layer loop

      do n = nlayers, 1, -1

!  Save some quantities

        salpha = sin(alpha_all(n-1))
        calpha = cos(alpha_all(n-1))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  Trapezium-rule integrated sourceterm multipliers + transmittances
!  -----------------------------------------------------------------

!  set up

        kn = raycon * extinction(n,w_excit)
        skn = step * kn
        tran_1 = exp ( - kn * ( cot_2 - cot_1 ) )
        csqt_1 = csq_1 * tran_1
        do j = 1, nfinelayers
          tranj(j) = exp ( - kn * ( cot_2 - cot_fine(n,j) ) )
        enddo

!  Line of sight transmittance factor

        lostrans(n)    = tran_1

!  Elastic scattering multiplier

        func_2 = attn(n-1,w_excit) *  csq_2
        func_1 = attn(n,w_excit)   * csqt_1
        sumr = half * ( func_1 + func_2 )
        do j = 1, nfinelayers
          func = attn_fine(n,j,w_excit) * tranj(j) * csq_fine(n,j)
          sumr  = sumr + func
        enddo
        emultipliers(n) = sumr * skn

!  Inelastic scattering multipliers (by bins)

        do w = 1, npoints_mono
          func_2 = attn(n-1,w) *  csq_2
          func_1 = attn(n,w)   * csqt_1
          sumr = half * ( func_1 + func_2 )
          do j = 1, nfinelayers
            func = attn_fine(n,j,w) * tranj(j) * csq_fine(n,j)
            sumr  = sumr + func
          enddo
          imultipliers(n,w) = sumr * skn
        enddo

!  update geometry

        csq_1 = csq_2
        cot_1 = cot_2

!  Finish layer

      ENDDO

!  finish

      RETURN
      END SUBROUTINE OUTGOING_INTEGRATION_MONO_UP_1

!

      SUBROUTINE OUTGOING_INTEGRATION_MONO_DN_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS,           & ! Inputs
          NPOINTS_MONO, W_EXCIT,                         & ! Inputs
          NLAYERS, NFINELAYERS, EXTINCTION,              & ! Inputs
          RADII, ALPHA_ALL, ALPHA_FINE, ATTN, ATTN_FINE, & ! Inputs
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS )           ! Outputs

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, HALF

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER  , INTENT(IN) :: maxlayers, maxfinelayers
      INTEGER  , INTENT(IN) :: maxpoints

!  control

      INTEGER  , INTENT(IN) :: nfinelayers, nlayers
      INTEGER  , INTENT(IN) :: npoints_mono, w_excit

!  radii

      REAL(FPK), INTENT(IN) :: radii   (0:maxlayers)

!  line of sight angles

      REAL(FPK), INTENT(IN) :: alpha_all  (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_fine (maxlayers,maxfinelayers)

!  Extinction

      REAL(FPK), INTENT(IN) :: extinction (maxlayers, maxpoints)

!  Attenuations

      REAL(FPK), INTENT(IN) :: attn      ( 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: attn_fine &
               ( maxlayers, maxfinelayers, maxpoints )

!  outputs
!  -------

      REAL(FPK), INTENT(OUT) :: imultipliers   (maxlayers,   maxpoints)
      REAL(FPK), INTENT(OUT) :: emultipliers   (maxlayers)
      REAL(FPK), INTENT(OUT) :: lostrans      (maxlayers)

!  local arrays
!  ------------

!  Local geometry arrays

      REAL(FPK) :: csq_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: cot_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: tranj     ( maxfinelayers)

!  help variables
!  --------------

      INTEGER   :: w, n, j
      REAL(FPK) :: sumr, salpha, calpha, dfine1, raycon, kn
      REAL(FPK) :: csq_1, csqt_1, csq_2, cot_1, cot_2
      REAL(FPK) :: func_1, func_2, tran_1, step, func, skn

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0_FPK

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        lostrans(n)     = zero
        emultipliers(n) = zero
        do w = 1, npoints_mono
          imultipliers(n,w) = zero
        enddo
      enddo

!  Ray constant at TOA

      salpha = sin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work Down from the top of the atmosphere
!  ========================================

!  initialise

      n = 0
      salpha = sin(alpha_all(n))
      calpha = cos(alpha_all(n))
      csq_1 = one / salpha / salpha
      cot_1 = calpha / salpha

!  Save some quantities

      do n = 1, nlayers
        do j = nfinelayers, 1, -1
          calpha = cos(alpha_fine(n,j))
          salpha = sin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = one / salpha / salpha
        enddo
      enddo

!  Start layer loop

      do n = 1, nlayers

!  Save some quantities

        salpha = sin(alpha_all(n))
        calpha = cos(alpha_all(n))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  Trapezium-rule integrated sourceterm multiplier + transmittance
!  ---------------------------------------------------------------

!  set up

        kn = raycon * extinction(n,w_excit)
        skn = step * kn
        tran_1 = exp ( - kn * ( cot_1 - cot_2 ) )
        csqt_1 = csq_1 * tran_1
        do j = 1, nfinelayers
          tranj(j) = exp ( - kn * ( cot_fine(n,j) - cot_2 ) )
        enddo

!  Line of sight transmittance factor

        lostrans(n)    = tran_1

!  Elastic multipliers

        func_1 = attn(n-1,w_excit) * csqt_1
        func_2 = attn(n,w_excit)   * csq_2
        sumr = half * ( func_1 + func_2 )
        do j = nfinelayers, 1, -1
          func = attn_fine(n,j,w_excit) * tranj(j) * csq_fine(n,j)
          sumr = sumr + func
        enddo
        emultipliers(n) = sumr * step * kn

!  Inelastic multipliers

        do w = 1, npoints_mono
          func_1 = attn(n-1,w) * csqt_1
          func_2 = attn(n,w)   * csq_2
          sumr = half * ( func_1 + func_2 )
          do j = nfinelayers, 1, -1
            func = attn_fine(n,j,w) * tranj(j) * csq_fine(n,j)
            sumr = sumr + func
          enddo
          imultipliers(n,w) = sumr * step * kn
        enddo

!  Check and alter code

        emultipliers(n) = imultipliers(n,w_excit)

!  update geometry

        csq_1 = csq_2
        cot_1 = cot_2

!  Finish layer

      ENDDO

!  finish

      RETURN
      END SUBROUTINE OUTGOING_INTEGRATION_MONO_DN_1

!

      END MODULE lrrs_corrections_1_m

