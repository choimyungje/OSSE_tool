
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
! # Subroutines in this Module                                  #
! #                                                             #
! #             LRRS_SSCORR_GENERAL_BIN_PLUS_1  (master)        #
! #             LRRS_SSCORR_GENERAL_MONO_PLUS_1 (master)        #
! #                                                             #
! #                 og_integration_bin_up_plus_1                #
! #                 og_integration_bin_dn_plus_1                #
! #                                                             #
! #                 og_integration_mono_up_plus_1               #
! #                 og_integration_mono_dn_plus_1               #
! #                                                             #
! #                 og_attenuations_plus_1                      #
! #                                                             #
! #    Single scatter exact calculations for outgoing LOS path  #
! #   Inelastic and Elastic output.                             #
! #                                                             #
! #  Programmed by R. Spurr, RT Solutions Inc.                  #
! #                                                             #
! #    First Draft,  April    2008                              #
! #    Second Draft, October  2008 with Column  Jacobians       #
! #    Third Draft,  November 2008 with Surface Jacobians extra #
! #    Fourth Draft, July     2009 with Profile Jacobians extra #
! #    Fifth Draft,  February 2011 Added Monochromatic stuff    #
! #    Sixth  Draft, March    2011 Use adjusted geometries      #
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

      MODULE lrrs_corrections_plus_1_m

!      USE LRRS_PARS_m, Only : LDU

      USE lrrs_geometry_1_m   , Only : OUTGOING_SPHERGEOM_FINE_UP_1, &
                                       OUTGOING_SPHERGEOM_FINE_DN_1

      PRIVATE  :: OG_INTEGRATION_BIN_UP_PLUS_1 , OG_INTEGRATION_BIN_DN_PLUS_1 , &
                  OG_INTEGRATION_MONO_UP_PLUS_1, OG_INTEGRATION_MONO_DN_PLUS_1, &
                  og_attenuations_plus_1

      PUBLIC   :: LRRS_SSCORR_GENERAL_BIN_PLUS_1,&
                  LRRS_SSCORR_GENERAL_MONO_PLUS_1

      CONTAINS

      SUBROUTINE LRRS_SSCORR_GENERAL_BIN_PLUS_1 &
         ( DO_FULL_RAMAN, DO_UPWELLING, DO_DNWELLING,                     & ! Input Mode flags
           DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,                           & ! Input Mode flags 2p5a
           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_DMSCALING,          & ! Input Mode flags
           DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,          & ! Input Surface flags
           DO_BRDF_Wav1, DO_Sleave_Wav1,                                  & ! Input Supplement flags
           DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,  & ! Inputs Lin.Control
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMNWFS,               & ! Inputs Lin.Control
           NKSTORAGE2, N_SURFACE_WFS, N_SLEAVE_WFS,                       & ! InputsLin.Control
           NLAYERS, NFINELAYERS, N_USER_STREAMS, N_USER_RELAZMS,          & ! Input Control Integers
           N_LOUTPUT, LEVELMASK_UP, LEVELMASK_DN,                         & ! Input levels for output
           NPOINTS_INNER, OFFSET_INNER, NPOINTS_OUTER, N_RRSBINS, BINMAP, & ! Input Raman bin control
           SOLAR_ANGLE, COSSCAT_NADIR_UP, COSSCAT_NADIR_DN,               & ! Input Geometry 2p5a
           SOLAR_ANGLES_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST,  & ! Input Geometry
           FLUXES_RANKED, ALBEDOS_RANKED, EARTH_RADIUS, HEIGHT_GRID,      & ! Input Fluxes/heights
           EXACTDB_BRDFUNC,  SLTERM_ISOTROPIC, SLTERM_USERANGLES,         & ! Input Surface stuff
           TRUNC_FACTORS, DELTAU_VERT_INPUT, OMEGA_UNSCALED,              & ! Input Optical Elastic
           OMEGAPHASFUNC_ELASTIC_UP,   OMEGAPHASFUNC_ELASTIC_DN,          & ! Input Optical Elastic
           OMEGAPHASFUNC_CABANNES_UP,  OMEGAPHASFUNC_CABANNES_DN,         & ! Input Optical Cabannes
           OMEGAMOMS_RRSLOSS_UNSCALED, OMEGAMOMS_RRSBIN_UNSCALED,         & ! Input Optical Raman
           SAVE_TRANS_USERM, SAVE_BEAMMULT_UP, SAVE_BEAMMULT_DN,          & ! Input Optical for NADIR SSCORR 2p5a
           LS_EXACTDB_BRDFUNC, LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES, & ! Input Lin Surface
           L_DELTAU_VERT_INPUT, L_TRUNC_FACTORS, L_OMEGA_UNSCALED,            & ! Input Lin Optical
           L_OMEGAPHASFUNC_ELASTIC_UP,   L_OMEGAPHASFUNC_ELASTIC_DN,          & ! Input Lin Optical Elastic  2p5a
           L_OMEGAPHASFUNC_CABANNES_UP,  L_OMEGAPHASFUNC_CABANNES_DN,         & ! Input Lin Optical Cabannes 2p5a
           L_OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSBIN_UNSCALED,         & ! Input Lin Optical Raman
           L_SAVE_TRANS_USERM, L_SAVE_BEAMMULT_UP, L_SAVE_BEAMMULT_DN,    & ! Input Optical for NADIR SSCORR 2p5a
           ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN,        & ! Outputs
           LC_ELASTIC_SS_UP, LP_ELASTIC_SS_UP, LS_ELASTIC_SS_UP,          & ! Outputs
           LC_ELASTIC_SS_DN, LP_ELASTIC_SS_DN,                            & ! Outputs
           LC_RAMAN_SS_UP,   LP_RAMAN_SS_UP,   LS_RAMAN_SS_UP,            & ! Outputs
           LC_RAMAN_SS_DN,   LP_RAMAN_SS_DN,                              & ! Outputs
           FAIL, MESSAGE)                                                   ! Outputs

!  Single scatter exact calculation for the outgoing LOS

!   Programmed by R. Spurr, RT Solutions Inc.
!    First Draft,  for LRRS Version 2.1, April 23rd 2008
!    Second draft, for LRRS Version 2.2, 23 July 2009

!   -- Rob mod 5/12/17 for 2p5a, Added DO_SSCORR_NADIR flag, plus precalculated multipliers/transmittances
!   -- Rob mod 5/12/17 for 2p5a, Using PHASFUNC input now, removed NMOMENTS_INPUT.

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, FOUR, PI4, DEG_TO_RAD,                &
                              MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES,   &
                              MAX_LAYERS, MAX_FINE_LAYERS, MAX_LAYERS_NK,           &
                              MAX_LOUTPUT, MAX_POINTS, MAX_BINS, MAX_ATMOSWFS,      &
                              MAX_SURFACEWFS, MAX_SLEAVEWFS

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

!  Profile WFS control inputs

      LOGICAL  , INTENT(IN) :: DO_PROFILE_WFS
      LOGICAL  , INTENT(IN) :: LAYER_VARY_FLAG   (MAX_LAYERS)
      INTEGER  , INTENT(IN) :: LAYER_VARY_NUMBER (MAX_LAYERS)

!  Linearization bookkeeping

      INTEGER  , INTENT(IN) :: NKSTORAGE2(MAX_LAYERS,0:MAX_LAYERS)

!  Column WFS control inputs

      LOGICAL  , INTENT(IN) :: DO_COLUMN_WFS
      INTEGER  , INTENT(IN) :: N_COLUMNWFS

!  Surface WFS

      LOGICAL  , INTENT(IN) :: DO_SURFACE_WFS
      INTEGER  , INTENT(IN) :: N_SURFACE_WFS

!  Surface-leaving WFS

      LOGICAL  , INTENT(IN) :: DO_SLEAVE_WFS
      INTEGER  , INTENT(IN) :: N_SLEAVE_WFS

!  Number of layers

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

!  WFSs of all optical properties
!  ------------------------------

!  Atmos
!   -- Rob mod 5/12/17 for 2p5a, Only need OMEGA_UNSCALED, plus PHASFUNC products

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_TRUNC_FACTORS     ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGA_UNSCALED    ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_UP  ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_DN  ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_CABANNES_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_CABANNES_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSLOSS_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSBIN_UNSCALED  ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Surface BRDF

      REAL(FPK), INTENT(IN) :: LS_EXACTDB_BRDFUNC ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

!  Surface-leaving

      REAL(fpk), INTENT(IN) :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAX_POINTS )
      REAL(fpk), INTENT(IN) :: LSSL_SLTERM_USERANGLES ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

!  NADIR correction - Linearized user stream transmittances, whole layers
!   -- Rob mod 5/12/17 for 2p5a

      REAL(FPK), INTENT(IN) :: L_SAVE_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  NADIR correction - Linearized pre-calculated multipliers
!   -- Rob mod 5/12/17 for 2p5a

      REAL(FPK), INTENT(IN) :: L_SAVE_BEAMMULT_UP ( MAX_ATMOSWFS, MAX_USER_STREAMS, 0:MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_SAVE_BEAMMULT_DN ( MAX_ATMOSWFS, MAX_USER_STREAMS, 0:MAX_LAYERS_NK, MAX_POINTS )

!  Output
!  ------

!mick fix 10/19/2015 - changed intent to inout

!  Single scatter Elastic and Raman Intensity results

      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_UP (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_DN (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_UP (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_DN (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Single scatter Elastic jacobian results

      REAL(FPK), INTENT(INOUT) :: LC_ELASTIC_SS_UP (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: LC_ELASTIC_SS_DN (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LP_ELASTIC_SS_UP (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: LP_ELASTIC_SS_DN (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LS_ELASTIC_SS_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Single scatter Raman jacobian results

      REAL(FPK), INTENT(INOUT) :: LC_RAMAN_SS_UP (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: LC_RAMAN_SS_DN (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LP_RAMAN_SS_UP (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: LP_RAMAN_SS_DN (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LS_RAMAN_SS_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

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

      REAL(FPK) :: extinction   ( max_layers, max_points )
      REAL(FPK) :: L_extinction ( MAX_ATMOSWFS, max_layers, max_points )

!  Saved Legendre polynomials. 
!   -- Rob mod 5/12/17 for 2p5a, Only need 0:2 now, drop MAX_MOMENTS_INPUT

      REAL(FPK) :: SS_PLEG_UP(MAX_LAYERS,0:2), SS_PLEG_DN(MAX_LAYERS,0:2)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(FPK) :: TMS(MAX_LAYERS,MAX_POINTS)

!  Linearized TMS factor

      REAL(FPK) :: L_TMS(MAX_ATMOSWFS,MAX_LAYERS,MAX_POINTS)

!  Local truncation factors for additional DELTAM scaling

!      REAL(FPK) :: SSFDEL ( MAXLAYERS )

!  Exact Phase function calculations

      REAL(FPK) :: ESCAT_UP(MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ESCAT_DN(MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: CSCAT_UP(MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: CSCAT_DN(MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: GSCAT_UP(MAX_LAYERS,MAX_BINS,MAX_POINTS)
      REAL(FPK) :: GSCAT_DN(MAX_LAYERS,MAX_BINS,MAX_POINTS)

!  Linearized Exact Phase function calculations

      REAL(FPK) :: L_ESCAT_UP(MAX_ATMOSWFS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_ESCAT_DN(MAX_ATMOSWFS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_CSCAT_UP(MAX_ATMOSWFS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_CSCAT_DN(MAX_ATMOSWFS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_GSCAT_UP (MAX_ATMOSWFS,MAX_LAYERS,MAX_BINS,MAX_POINTS)
      REAL(FPK) :: L_GSCAT_DN (MAX_ATMOSWFS,MAX_LAYERS,MAX_BINS,MAX_POINTS)

!  Cumulative single scatter source terms

      REAL(FPK) :: ESS_CUMSOURCE_UP(0:MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ESS_CUMSOURCE_DN(0:MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ISS_CUMSOURCE_UP(0:MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ISS_CUMSOURCE_DN(0:MAX_LAYERS,MAX_POINTS)

!  Linearized Cumulative single scatter source terms

      REAL(FPK) :: L_ESS_CUMSOURCE(MAX_ATMOSWFS,MAX_POINTS)
      REAL(FPK) :: L_ISS_CUMSOURCE(MAX_ATMOSWFS,MAX_POINTS)

!  Linearized DB sources

      REAL(FPK) :: LS_ESS_CUMSOURCE_UP(MAX_SURFACEWFS,MAX_POINTS)

!  Outgoing sphericity stuff
!  -------------------------

!  Whole layer Elastic and Inelastic multipliers

      REAL(FPK) :: &
         UP_IMULT(MAX_LAYERS,MAX_BINS,MAX_POINTS), &
         DN_IMULT(MAX_LAYERS,MAX_BINS,MAX_POINTS), &
         UP_EMULT(MAX_LAYERS,MAX_POINTS), &
         DN_EMULT(MAX_LAYERS,MAX_POINTS), &
         UP_LOSTRANS(MAX_LAYERS,MAX_POINTS), &
         DN_LOSTRANS(MAX_LAYERS,MAX_POINTS)

!  Linearized Whole layer Elastic and Inelastic multipliers

      REAL(FPK) :: L_UP_IMULT (MAX_ATMOSWFS,0:MAX_LAYERS,MAX_LAYERS,MAX_BINS,MAX_POINTS)
      REAL(FPK) :: L_DN_IMULT (MAX_ATMOSWFS,0:MAX_LAYERS,MAX_LAYERS,MAX_BINS,MAX_POINTS)

      REAL(FPK) :: L_UP_EMULT (MAX_ATMOSWFS,0:MAX_LAYERS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_DN_EMULT (MAX_ATMOSWFS,0:MAX_LAYERS,MAX_LAYERS,MAX_POINTS)

      REAL(FPK) :: &
         L_UP_LOSTRANS(MAX_ATMOSWFS,MAX_LAYERS,MAX_POINTS), &
         L_DN_LOSTRANS(MAX_ATMOSWFS,MAX_LAYERS,MAX_POINTS)

!  Solar beam attenuations

      REAL(FPK) :: attn      ( 0:max_layers, max_points )
      REAL(FPK) :: attn_fine ( max_layers, max_fine_layers, max_points )

      REAL(FPK) :: l_attn ( MAX_ATMOSWFS, 0:max_layers, 0:max_layers, max_points )
      REAL(FPK) :: l_attn_fine ( MAX_ATMOSWFS, 0:max_layers, max_layers, max_fine_layers, max_points )

!  Local variables
!  ---------------

!  Indices

      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, W, B, NK
      INTEGER   :: UTA, UM, IA, NC, V, WR, WB, Q, Q1, K, KS
      INTEGER   :: Sleave_IDX, Brdf_IDX

!  Help variables
!   -- Rob mod 5/12/17 for 2p5a, removed DFL1, DFL2

      REAL(FPK) :: boa_attn, ssflux, x0_boa, x0_boa_4, ratio, Brdf_Exact, &
                   Sleave_Iso, Sleave_User, L_Sleave_Iso, L_Sleave_User
      REAL(FPK) :: GSOURCE, CSOURCE, ESOURCE, ISOURCE, SOURCE
      REAL(FPK) :: HELP, COSSCAT, CONTRIB, FACTOR, OMEGA_LOCAL, BS, T0, T2, P2
      REAL(FPK) :: PHAS_E, PHAS_CB, PHAS_LR, PHAS_LG
      REAL(FPK) :: REFLEC(MAX_POINTS), LS_REFLEC(MAX_SURFACEWFS,MAX_POINTS)

!  WFS help variables

      REAL(FPK) :: L_OF, LSS, LCSS, LGSS, LBSS
      REAL(FPK) :: L_PHAS_E, L_PHAS_CB, L_PHAS_LR, L_PHAS_LG
      REAL(FPK) :: L_ESOURCE, L_ISOURCE, L_SSCORRECTION, L_OMEGA_LOCAL
      REAL(FPK) :: L_BOA_ATTN, L_LOSS, L_CSOURCE, L_GSOURCE

!  Bookkeeping

      LOGICAL   :: do_atmos_wfs
      INTEGER   :: NTYPE_VARY, K_PARAMETERS
      LOGICAL   :: DO_RTSOL_VARY(max_layers)
      INTEGER   :: NPARAMS_VARY(max_layers)
      INTEGER   :: KINDEX(max_layers)

!  Sleave variables

      REAL(fpk)  :: LSSL_EXACTDB_SOURCE ( MAX_SLEAVEWFS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(fpk)  :: LSSL_DB_CUMSOURCE   ( MAX_GEOMETRIES, MAX_POINTS )

!  Local surface flag

      LOGICAL  , PARAMETER :: DO_INCLUDE_SURFACE = .true.

!  Local value

      REAL(FPK), PARAMETER :: ls_albedos_ranked = one

!  Set up operations
!  -----------------

!  SS flux scaling factor

      SSFLUX = one / PI4

!  Fine layering flag

      do_fine = .true.

!  Number of fine layers now a control input (7/25/08)

!  General atmospheric WFS term

      DO_ATMOS_WFS = ( DO_PROFILE_WFS .OR. DO_COLUMN_WFS )

!  Parameter control
!  -----------------

      IF ( DO_PROFILE_WFS ) THEN
        NTYPE_VARY = NLAYERS
        DO N = 1, NLAYERS
         DO_RTSOL_VARY(n) = LAYER_VARY_FLAG(n)
         NPARAMS_VARY(n)  = LAYER_VARY_NUMBER(n)
         KINDEX(N) = N
        ENDDO
      ELSE IF ( DO_COLUMN_WFS ) THEN
        NTYPE_VARY = 1
        DO_RTSOL_VARY(1) = .TRUE.
        NPARAMS_VARY(1)  = N_COLUMNWFS
        KINDEX(1) = 0
      ENDIF

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

!  Create Linearized TMS factor for each layer
!   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)
!      SHOULD ONLY BE USED WITH THE ENERGY BALANCING
!   -- Rob mod 5/12/17 for 2p5a, Use L_OMEGA_UNSCALED and OMEGA_UNSCALED

      IF ( DO_ATMOS_WFS ) THEN
        DO W = 1, NPOINTS_INNER
          WR = W + OFFSET_INNER
          DO K = 1, NLAYERS
            K_PARAMETERS = LAYER_VARY_NUMBER(K)
            DO Q = 1, K_PARAMETERS
              IF ( DO_DMSCALING ) THEN
                OMEGA_LOCAL    = OMEGA_UNSCALED(K,WR)
                L_OMEGA_LOCAL  = L_OMEGA_UNSCALED(Q,K,WR)
                L_OF = L_OMEGA_LOCAL * TRUNC_FACTORS(K,WR) + L_TRUNC_FACTORS(Q,K,WR) * OMEGA_LOCAL
                L_TMS(Q,K,W) = TMS(K,W) * TMS(K,W) * L_OF
              ELSE
                L_TMS(Q,K,W) = zero
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Create extinctions (using scaled optical depths)

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        DO W = 1, NPOINTS_OUTER
          EXTINCTION(N,W) = DELTAU_VERT_INPUT(N,W) / HELP
        ENDDO
      ENDDO

!  Linearized extinctions

      IF ( DO_ATMOS_WFS ) THEN
       DO K = 1, NLAYERS
         HELP = HEIGHT_GRID(K-1) - HEIGHT_GRID(K)
         K_PARAMETERS = LAYER_VARY_NUMBER(K)
         DO Q = 1, K_PARAMETERS
           DO W = 1, NPOINTS_OUTER
             L_EXTINCTION(Q,K,W) = L_DELTAU_VERT_INPUT(Q,K,W) / HELP
           ENDDO
         ENDDO
       ENDDO
      ENDIF

!  Additional Delta-M scaling
!  --------------------------

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
        ALPHA_BOA = USER_ANGLES_ADJUST(UM)

!  UPWELLING NADIR CORRECTION: multipliers, transmittances
!  =======================================================

!  Upwelling Multipliers and Transmittances for the SSCORR_NADIR case
!     Does not use the ADJUSTED geometrical angles, Same for each azimuth.
!   -- Rob mod 5/12/17 for 2p5a, SSCORR_NADIR values are already pre-calculated

        IF ( DO_SSCORR_NADIR .and.DO_UPWELLING ) THEN

          DO W = 1, NPOINTS_INNER
            WR = W + OFFSET_INNER

!  Multipliers

            DO N = 1, NLAYERS
              UP_LOSTRANS(N,W) = SAVE_TRANS_USERM(UM,N,W)
              UP_EMULT(N,W)    = SAVE_BEAMMULT_UP(UM,N,W)
              DO B = 1, N_RRSBINS(W)
                WB = BINMAP(B,W)
                UP_IMULT(N,B,W) = SAVE_BEAMMULT_UP(UM,N,WB)
              ENDDO
            ENDDO

!  Linearized Multipliers

            DO N = 1, NLAYERS
              DO K = 1, NTYPE_VARY
                IF ( DO_RTSOL_VARY(K) ) THEN
                  K_PARAMETERS = NPARAMS_VARY(K) ; KS = KINDEX(K); NK = NKSTORAGE2(N,KS)
                  IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                    DO Q = 1, K_PARAMETERS
                      L_UP_LOSTRANS(Q,N,W) = L_SAVE_TRANS_USERM(Q,UM,N,W)
                    ENDDO
                  ENDIF
                  DO Q = 1, K_PARAMETERS
                    L_UP_EMULT(Q,KS,N,W) = L_SAVE_BEAMMULT_UP(Q,UM,NK,W)
                    DO B = 1, N_RRSBINS(W)
                      WB = BINMAP(B,W)
                      L_UP_IMULT(Q,KS,N,B,W) = L_SAVE_BEAMMULT_UP(Q,UM,NK,WB)
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO
            ENDDO

!  End upwelling Nadir multipliers

          ENDDO
        ENDIF

!  DOWNWELLING NADIR CORRECTION: multipliers, transmittances
!  =========================================================

!  Downwelling Multipliers and Transmittances for the SSCORR_NADIR case
!     Does not use the ADJUSTED geometrical angles, Same for each azimuth.
!   -- Rob mod 5/12/17 for 2p5a, SSCORR_NADIR values are already pre-calculated

        IF ( DO_SSCORR_NADIR .and.DO_DNWELLING ) THEN

          DO W = 1, NPOINTS_INNER
            WR = W + OFFSET_INNER

!  Multipliers

            DO N = 1, NLAYERS
              DN_LOSTRANS(N,W) = SAVE_TRANS_USERM(UM,N,W)
              DN_EMULT(N,W)    = SAVE_BEAMMULT_DN(UM,N,W)
              DO B = 1, N_RRSBINS(W)
                WB = BINMAP(B,W)
                DN_IMULT(N,B,W) = SAVE_BEAMMULT_DN(UM,N,WB)
              ENDDO
            ENDDO

!  Linearized Multipliers

            DO N = 1, NLAYERS
              DO K = 1, NTYPE_VARY
                IF ( DO_RTSOL_VARY(K) ) THEN
                  K_PARAMETERS = NPARAMS_VARY(K) ; KS = KINDEX(K) ; NK = NKSTORAGE2(N,KS)
                  IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                    DO Q = 1, K_PARAMETERS
                      L_DN_LOSTRANS(Q,N,W) = L_SAVE_TRANS_USERM(Q,UM,N,W)
                    ENDDO
                  ENDIF
                  DO Q = 1, K_PARAMETERS
                    L_DN_EMULT(Q,KS,N,W) = L_SAVE_BEAMMULT_DN(Q,UM,NK,W)
                    DO B = 1, N_RRSBINS(W)
                      WB = BINMAP(B,W)
                      L_DN_IMULT(Q,KS,N,B,W) = L_SAVE_BEAMMULT_DN(Q,UM,NK,WB)
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO
            ENDDO

!  End Downwelling Nadir multipliers

          ENDDO
        ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  Start the Azimuth loop
!  ----------------------

        DO IA = 1, N_USER_RELAZMS

!  Adjusted geometry for the Outgoing correction
!mick mod 11/28/2018 - moved defining of ALPHA_BOA further up

          !ALPHA_BOA = USER_ANGLES_ADJUST(UM)
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

            call og_attenuations_plus_1 &
        ( max_layers, max_fine_layers, max_points, MAX_ATMOSWFS, & ! Inputs
          npoints_outer, nlayers, nfinelayers, do_profile_wfs,   & ! Inputs
          layer_vary_flag, layer_vary_number, do_column_wfs,     & ! Inputs
          N_COLUMNWFS, extinction, L_extinction,                 & ! Inputs
          sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,    & ! Inputs
          attn, attn_fine, l_attn, l_attn_fine )                   ! Outputs

!  Multipliers, transmittances

            call og_integration_bin_up_plus_1 &
        ( max_layers, max_fine_layers, max_points, max_bins,        & ! Inputs
          MAX_ATMOSWFS, do_full_raman, npoints_inner, offset_inner, & ! Inputs
          n_rrsbins, binmap, nlayers, nfinelayers,                  & ! Inputs
          do_profile_wfs, layer_vary_flag, layer_vary_number,       & ! Inputs
          do_column_wfs,  n_columnwfs,                              & ! Inputs
          extinction, l_extinction, radii, alpha_all,               & ! Inputs
          alpha_fine, attn, attn_fine, l_attn, l_attn_fine,         & ! Inputs
          up_emult, up_imult, up_lostrans,                          & ! Inputs
          l_up_emult, l_up_imult, l_up_lostrans )                     ! Outputs

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
            P2 = SS_PLEG_UP(1,2)

!  Elastic scatting Phase functions (multiplied by TMS factor)
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_E with Phase function input
!   -- PHAS_E   = DOT_PRODUCT(OMEGAMOMS_ELASTIC_UNSCALED(N,0:NMOMENTS_INPUT,WR),SS_PLEG_UP(1,0:NMOMENTS_INPUT))
!   -- L_PHAS_E = DOT_PRODUCT(L_OMEGAMOMS_ELASTIC_UNSCALED(Q,N,0:NMOMENTS_INPUT,WR),SS_PLEG_UP(1,0:NMOMENTS_INPUT))

            DO N = 1, NLAYERS
              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
                PHAS_E = OMEGAPHASFUNC_ELASTIC_UP(N,V,WR)
                ESCAT_UP(N,W) = PHAS_E * TMS(N,W)
                IF ( LAYER_VARY_FLAG(N) ) THEN
                  DO Q = 1, LAYER_VARY_NUMBER(N)
                    L_PHAS_E =  L_OMEGAPHASFUNC_ELASTIC_UP(Q,N,V,WR)
                    L_ESCAT_UP(Q,N,W) = L_PHAS_E * TMS(N,W) + PHAS_E * L_TMS(Q,N,W)
                  ENDDO
                ENDIF
              ENDDO
            ENDDO

!  Add inelastic scattering source functions if flagged

            IF ( DO_FULL_RAMAN ) THEN

!  Energy-balancing method, CSCAT_UP = Elastic - Loss
! ---------------------------------------------------

              IF ( DO_ENERGY_BALANCING ) THEN
                DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER
                  DO N = 1, NLAYERS
                    PHAS_LR = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,W) + &
                              OMEGAMOMS_RRSLOSS_UNSCALED(N,2,W) * P2
                    CSCAT_UP(N,W) = ESCAT_UP(N,W) - PHAS_LR * TMS(N,W)
                    IF ( LAYER_VARY_FLAG(N) ) THEN
                      DO Q = 1, LAYER_VARY_NUMBER(N)
                        L_PHAS_LR = L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,0,W) + &
                                    L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,2,W) * P2
                        L_LOSS = PHAS_LR * L_TMS(Q,N,W) + L_PHAS_LR * TMS(N,W)
                        L_CSCAT_UP(Q,N,W) = L_ESCAT_UP(Q,N,W) - L_LOSS
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF

!  Cabannes-Raman method, CSCAT_UP = Cabannes (direct)
!  ---------------------------------------------------

!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_CB with Phase function input
!   -- PHAS_CB = DOT_PRODUCT(OMEGAMOMS_CABANNES_UNSCALED(N,0:NMOMENTS_INPUT,W),SS_PLEG_UP(1,0:NMOMENTS_INPUT))

              IF ( DO_CABANNES_RAMAN ) THEN
                DO N = 1, NLAYERS
                  DO W = 1, NPOINTS_INNER
                    WR = W + OFFSET_INNER
                    PHAS_CB = OMEGAPHASFUNC_CABANNES_UP(N,V,W)
                    CSCAT_UP(N,W) = PHAS_CB * TMS(N,W)
                    IF ( LAYER_VARY_FLAG(N) ) THEN
                      DO Q = 1, LAYER_VARY_NUMBER(N)
                        L_PHAS_CB = L_OMEGAPHASFUNC_CABANNES_UP(Q,N,V,W)
                        L_CSCAT_UP(Q,N,W) = L_PHAS_CB * TMS(N,W) + PHAS_CB * L_TMS(Q,N,W)
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF

!  Gain term (Same for both realizations)
!  --------------------------------------
!mick fix 11/28/2018 - added defining of WR here

              DO N = 1, NLAYERS
                DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER
                  DO B = 1, N_RRSBINS(W)
                    WB = BINMAP(B,W)
                    RATIO   = FLUXES_RANKED(WB)/FLUXES_RANKED(WR)
                    CONTRIB = OMEGAMOMS_RRSBIN_UNSCALED(N,0,B,W) + &
                              OMEGAMOMS_RRSBIN_UNSCALED(N,2,B,W) * P2
                    PHAS_LG  =  CONTRIB * RATIO
                    GSCAT_UP(N,B,W) = PHAS_LG * TMS(N,W)
                    IF ( LAYER_VARY_FLAG(N) ) THEN
                      DO Q = 1, LAYER_VARY_NUMBER(N)
                        T0 = L_OMEGAMOMS_RRSBIN_UNSCALED(Q,N,0,B,W)
                        T2 = L_OMEGAMOMS_RRSBIN_UNSCALED(Q,N,2,B,W) * P2
                        L_PHAS_LG = ( T0 + T2 ) * RATIO
                        L_GSCAT_UP(Q,N,B,W) = L_PHAS_LG * TMS(N,W) + PHAS_LG * L_TMS(Q,N,W)
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

!  End inelastic scattering option

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

            call og_attenuations_plus_1 &
        ( max_layers, max_fine_layers, max_points, MAX_ATMOSWFS, & ! Inputs
          npoints_outer, nlayers, nfinelayers, do_profile_wfs,   & ! Inputs
          layer_vary_flag, layer_vary_number, do_column_wfs,     & ! Inputs
          N_COLUMNWFS, extinction, L_extinction,                 & ! Inputs
          sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,    & ! Inputs
          attn, attn_fine, l_attn, l_attn_fine )                   ! Outputs

!  Multipliers and transmittances

            call og_integration_bin_dn_plus_1 &
        ( max_layers, max_fine_layers, max_points, max_bins,        & ! Inputs
          MAX_ATMOSWFS, do_full_raman, npoints_inner, offset_inner, & ! Inputs
          n_rrsbins, binmap, nlayers, nfinelayers,                  & ! Inputs
          do_profile_wfs, layer_vary_flag, layer_vary_number,       & ! Inputs
          do_column_wfs,  N_COLUMNWFS,                              & ! Inputs
          extinction, l_extinction, radii, alpha_all,               & ! Inputs
          alpha_fine, attn, attn_fine, l_attn, l_attn_fine,         & ! Inputs
          dn_emult, dn_imult, dn_lostrans,                          & ! Inputs
          l_dn_emult, l_dn_imult, l_dn_lostrans )                     ! Outputs

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
            P2 = SS_PLEG_DN(1,2)

!  Elastic scattering Phase functions (multiplied by TMS factor)
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_E with Phase function input
!   -- PHAS_E   = DOT_PRODUCT(OMEGAMOMS_ELASTIC_UNSCALED  (N,0:NMOMENTS_INPUT,WR),  SS_PLEG_DN(1,0:NMOMENTS_INPUT))
!   -- L_PHAS_E = DOT_PRODUCT(L_OMEGAMOMS_ELASTIC_UNSCALED(Q,N,0:NMOMENTS_INPUT,WR),SS_PLEG_DN(1,0:NMOMENTS_INPUT))

            DO N = 1, NLAYERS
              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
                PHAS_E = OMEGAPHASFUNC_ELASTIC_DN(N,V,WR)
                ESCAT_DN(N,W) = PHAS_E * TMS(N,W)
                IF ( LAYER_VARY_FLAG(N) ) THEN
                  DO Q = 1, LAYER_VARY_NUMBER(N)
                    L_PHAS_E =  L_OMEGAPHASFUNC_ELASTIC_DN(Q,N,V,WR)
                    L_ESCAT_DN(Q,N,W) = L_PHAS_E * TMS(N,W) + PHAS_E * L_TMS(Q,N,W)
                  ENDDO
                ENDIF
              ENDDO
            ENDDO

!  Add inelastic scatting Phase functions if flagged

            IF ( DO_FULL_RAMAN ) THEN

!  Energy-balancing method, CSCAT_DN = Elastic - Loss
! ---------------------------------------------------

              IF ( DO_ENERGY_BALANCING ) THEN
                DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER
                  DO N = 1, NLAYERS
                    PHAS_LR = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,W) + &
                              OMEGAMOMS_RRSLOSS_UNSCALED(N,2,W) * P2
                    CSCAT_DN(N,W) = ESCAT_DN(N,W) - PHAS_LR * TMS(N,W)
                    IF ( LAYER_VARY_FLAG(N) ) THEN
                      DO Q = 1, LAYER_VARY_NUMBER(N)
                        L_PHAS_LR = L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,0,W) + &
                                    L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,2,W) * P2
                        L_LOSS = PHAS_LR * L_TMS(Q,N,W) + L_PHAS_LR * TMS(N,W)
                        L_CSCAT_DN(Q,N,W) = L_ESCAT_DN(Q,N,W) - L_LOSS
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF

!  Cabannes-Raman method, CSCAT_DN = Cabannes (direct)
!  ---------------------------------------------------

!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_CB with Phase function input
!   -- PHAS_CB = DOT_PRODUCT(OMEGAMOMS_CABANNES_UNSCALED(N,0:NMOMENTS_INPUT,W),SS_PLEG_DN(1,0:NMOMENTS_INPUT))

              IF ( DO_CABANNES_RAMAN ) THEN
                DO N = 1, NLAYERS
                  DO W = 1, NPOINTS_INNER
                    WR = W + OFFSET_INNER
                    PHAS_CB = OMEGAPHASFUNC_CABANNES_DN(N,V,W)
                    CSCAT_DN(N,W) = PHAS_CB * TMS(N,W)
                    IF ( LAYER_VARY_FLAG(N) ) THEN
                      DO Q = 1, LAYER_VARY_NUMBER(N)
                        L_PHAS_CB = L_OMEGAPHASFUNC_CABANNES_DN(Q,N,V,W)
                        L_CSCAT_DN(Q,N,W) = L_PHAS_CB * TMS(N,W) + PHAS_CB * L_TMS(Q,N,W)
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF

!  Gain term (Same for both realizations)
!  --------------------------------------
!mick fix 11/28/2018 - added defining of WR here

              DO N = 1, NLAYERS
                DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER
                  DO B = 1, N_RRSBINS(W)
                    WB = BINMAP(B,W)
                    RATIO   = FLUXES_RANKED(WB)/FLUXES_RANKED(WR)
                    CONTRIB = OMEGAMOMS_RRSBIN_UNSCALED(N,0,B,W) + &
                              OMEGAMOMS_RRSBIN_UNSCALED(N,2,B,W) * P2
                    PHAS_LG  =  CONTRIB * RATIO
                    GSCAT_DN(N,B,W) = PHAS_LG * TMS(N,W)
                    IF ( LAYER_VARY_FLAG(N) ) THEN
                      DO Q = 1, LAYER_VARY_NUMBER(N)
                        T0 = L_OMEGAMOMS_RRSBIN_UNSCALED(Q,N,0,B,W)
                        T2 = L_OMEGAMOMS_RRSBIN_UNSCALED(Q,N,2,B,W) * P2
                        L_PHAS_LG = ( T0 + T2 ) * RATIO
                        L_GSCAT_DN(Q,N,B,W) = L_PHAS_LG * TMS(N,W) + PHAS_LG * L_TMS(Q,N,W)
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

!  End inelastic scattering option

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

!  Initialize cumulative source terms = F.A.mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux, A = albedo
!    Only require the Stokes total intensity component.
!    This contribution is purely elastic.
!    This line Replaced: X0_BOA = cos(SOLAR_ANGLE * DEG_TO_RAD)
!    Introduction of the Adjusted Gometry value, 18 March 2011

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
!  Multiply terms by 4pi to get correct normalization. 7/30/12
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

!  Initialize optical depth loop

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
                  ESOURCE = ESCAT_UP(N,W) * UP_EMULT(N,W)
                  ESS_CUMSOURCE_UP(NC,W) =  ESOURCE + &
                       UP_LOSTRANS(N,W) * ESS_CUMSOURCE_UP(NC-1,W)
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO W = 1, NPOINTS_INNER
                    CSOURCE = CSCAT_UP(N,W) * UP_EMULT(N,W)
                    GSOURCE = zero
                    DO B = 1, N_RRSBINS(W)
                      BS = GSCAT_UP(N,B,W) * UP_IMULT(N,B,W)
                      GSOURCE = GSOURCE + BS
                    ENDDO
                    ISOURCE = CSOURCE + GSOURCE
                    ISS_CUMSOURCE_UP(NC,W) = ISOURCE + &
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

!  End optical depth loop and Upwelling clause

            ENDDO
          ENDIF

!  Recurrence relation for the DOWNWELLING intensity
!  =================================================

          IF ( DO_DNWELLING ) THEN

!  Initialize cumulative source terms, and optical depth loop

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
                  ESOURCE = ESCAT_DN(N,W) * DN_EMULT(N,W)
                  ESS_CUMSOURCE_DN(NC,W) = ESOURCE + &
                       DN_LOSTRANS(N,W) * ESS_CUMSOURCE_DN(NC-1,W)
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO W = 1, NPOINTS_INNER
                    CSOURCE = CSCAT_DN(N,W) * DN_EMULT(N,W)
                    GSOURCE = zero
                    DO B = 1, N_RRSBINS(W)
                      BS = GSCAT_DN (N,B,W) * DN_IMULT(N,B,W)
                      GSOURCE = GSOURCE + BS
                    ENDDO
                    ISOURCE = CSOURCE + GSOURCE
                    ISS_CUMSOURCE_DN(NC,W) = ISOURCE + &
                       DN_LOSTRANS(N,W) * ISS_CUMSOURCE_DN(NC-1,W)
                  ENDDO
                ENDIF
              ENDDO


!  Ongrid output:
!    Set final cumulative source and correct Stokes vector

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

!  Recurrence relation for the UPWELLING Atmsopheric Jacobians
!  ===========================================================

          IF ( DO_UPWELLING .AND. DO_ATMOS_WFS ) THEN

!  Start the main layer variation loop
!  -----------------------------------

           DO K = 1, NTYPE_VARY
            IF ( DO_RTSOL_VARY(K) ) THEN
             K_PARAMETERS = NPARAMS_VARY(K)
             KS = KINDEX(K)

!  Initialize cumulative source terms = F.A.mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux, A = albedo
!    Only require the Stokes total intensity component.
!    This contribution is purely elastic.

             IF ( DO_INCLUDE_SURFACE ) THEN
              DO Q = 1, K_PARAMETERS
                DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER
                  FACTOR = REFLEC(W)
                  L_BOA_ATTN = X0_BOA_4 * L_ATTN(Q,KS,NLAYERS,WR)
                  L_ESS_CUMSOURCE(Q,W) = FACTOR * L_BOA_ATTN
                ENDDO
              ENDDO
             ELSE
              DO Q = 1, K_PARAMETERS
                DO W = 1, NPOINTS_INNER
                  L_ESS_CUMSOURCE(Q,W) = zero
                ENDDO
              ENDDO
             ENDIF

!  Duplicate the result for the inelastic case

             IF ( DO_FULL_RAMAN ) THEN
              DO Q = 1, K_PARAMETERS
                DO W = 1, NPOINTS_INNER
                 L_ISS_CUMSOURCE(Q,W) = L_ESS_CUMSOURCE(Q,W)
                ENDDO
              ENDDO
             ENDIF

!  Initialize optical depth loop

             NSTART = NLAYERS
             NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

             DO UTA = N_LOUTPUT, 1, -1

!  Layer index for given optical depth

              NLEVEL = LEVELMASK_UP(UTA)
              NUT    = NLEVEL + 1

!  Cumulative single scatter source terms :
!      For loop over layers working upwards to level NUT,

              DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N

!  Elastic terms

               IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                DO W = 1, NPOINTS_INNER
                  DO Q = 1, K_PARAMETERS
                    L_ESOURCE = ESCAT_UP(N,W)   * L_UP_EMULT(Q,KS,N,W) &
                            + L_ESCAT_UP(Q,N,W) *   UP_EMULT(N,W)
                    L_ESS_CUMSOURCE(Q,W) = L_ESOURCE &
                    +   UP_LOSTRANS(N,W)   * L_ESS_CUMSOURCE(Q,W) &
                    + L_UP_LOSTRANS(Q,N,W) *   ESS_CUMSOURCE_UP(NC-1,W)
                  ENDDO
                ENDDO
               ELSE
                DO W = 1, NPOINTS_INNER
                 DO Q = 1, K_PARAMETERS
                  L_ESOURCE = ESCAT_UP(N,W) * L_UP_EMULT(Q,KS,N,W)
                  L_ESS_CUMSOURCE(Q,W) = L_ESOURCE &
                  +   UP_LOSTRANS(N,W)   * L_ESS_CUMSOURCE(Q,W)
                 ENDDO
                ENDDO
               ENDIF

!  Raman terms for Inelastic scattering

               IF ( DO_FULL_RAMAN ) THEN
                IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                  DO W = 1, NPOINTS_INNER
                    DO Q = 1, K_PARAMETERS
                      L_CSOURCE = CSCAT_UP(N,W)   * L_UP_EMULT(Q,KS,N,W) &
                              + L_CSCAT_UP(Q,N,W) *   UP_EMULT(N,W)
                      L_GSOURCE = zero
                      DO B = 1, N_RRSBINS(W)
                        LBSS = GSCAT_UP(N,B,W) * L_UP_IMULT(Q,KS,N,B,W) &
                           + L_GSCAT_UP(Q,N,B,W) * UP_IMULT(N,B,W)
                        L_GSOURCE = L_GSOURCE + LBSS
                      ENDDO
                      L_ISOURCE  = L_CSOURCE + L_GSOURCE
                      L_ISS_CUMSOURCE(Q,W) = L_ISOURCE &
                     +   UP_LOSTRANS(N,W)   * L_ISS_CUMSOURCE(Q,W) &
                     + L_UP_LOSTRANS(Q,N,W) *   ISS_CUMSOURCE_UP(NC-1,W)
                    ENDDO
                  ENDDO
                ELSE
                  DO W = 1, NPOINTS_INNER
                    DO Q = 1, K_PARAMETERS
                      LCSS =  CSCAT_UP(N,W) *   L_UP_EMULT(Q,KS,N,W)
                      LGSS = zero
                      DO B = 1, N_RRSBINS(W)
                        LBSS =  GSCAT_UP(N,B,W) * L_UP_IMULT(Q,KS,N,B,W)
                        LGSS = LGSS + LBSS
                      ENDDO
                      LSS = LCSS + LGSS
                      L_ISS_CUMSOURCE(Q,W) = LSS &
                           + UP_LOSTRANS(N,W) * L_ISS_CUMSOURCE(Q,W)
                    ENDDO
                  ENDDO
                ENDIF
               ENDIF

!  End layer source term loop

              ENDDO

!  Ongrid output--------
!     Set final cumulative source and single scatter intensity

!  If N = K (the layer that is varying)
!    add WFS of additional partial layer source term =
!        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)
!  Variations when N > K
!    add WFS of additional partial layer source term =
!         Exact_Scat * L_Multiplier(k)

!  Profile WFS

              IF ( DO_PROFILE_WFS ) THEN
                DO W = 1, NPOINTS_INNER
                 DO Q = 1, K_PARAMETERS
                  L_ESOURCE = L_ESS_CUMSOURCE(Q,W)
                  L_SSCORRECTION = SSFLUX * L_ESOURCE
                  LP_ELASTIC_SS_UP(Q,K,UTA,V,W) = L_SSCORRECTION
                 ENDDO
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                 DO W = 1, NPOINTS_INNER
                  DO Q = 1, K_PARAMETERS
                   L_ISOURCE = L_ISS_CUMSOURCE(Q,W)
                   L_SSCORRECTION = SSFLUX * L_ISOURCE
                   LP_RAMAN_SS_UP(Q,K,UTA,V,W) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDIF
              ENDIF

!  Column WFS

              IF ( DO_COLUMN_WFS ) THEN
                DO W = 1, NPOINTS_INNER
                  DO Q = 1, K_PARAMETERS
                    L_ESOURCE = L_ESS_CUMSOURCE(Q,W)
                    L_SSCORRECTION = SSFLUX * L_ESOURCE
                    LC_ELASTIC_SS_UP(Q,UTA,V,W) = L_SSCORRECTION
                  ENDDO
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO W = 1, NPOINTS_INNER
                    DO Q = 1, K_PARAMETERS
                      L_ISOURCE = L_ISS_CUMSOURCE(Q,W)
                      L_SSCORRECTION = SSFLUX * L_ISOURCE
                      LC_RAMAN_SS_UP(Q,UTA,V,W) = L_SSCORRECTION
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF

!  Check for updating the recursion

              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
              NUT_PREV = NUT

!  End optical depth loop

             ENDDO

!  End loop over varying layers

            ENDIF
           ENDDO

!  End Upwelling Jacobian clause

          ENDIF

!  Recurrence relation for the DOWNWELLING Atmospheric Jacobians
!  =============================================================

          IF ( DO_DNWELLING .AND. DO_ATMOS_WFS ) THEN

!  Start the main layer variation loop
!  -----------------------------------

           DO K = 1, NTYPE_VARY
            IF ( DO_RTSOL_VARY(K) ) THEN
             K_PARAMETERS = NPARAMS_VARY(K)
             KS = KINDEX(K)

!  Initialise elastic source term

             DO Q = 1, K_PARAMETERS
              DO W = 1, NPOINTS_INNER
                L_ESS_CUMSOURCE(Q,W) = zero
              ENDDO
             ENDDO

!  Duplicate the result for the inelastic case

             IF ( DO_FULL_RAMAN ) THEN
              DO Q = 1, K_PARAMETERS
                DO W = 1, NPOINTS_INNER
                  L_ISS_CUMSOURCE(Q,W) = L_ESS_CUMSOURCE(Q,W)
                ENDDO
              ENDDO
             ENDIF

!  Initialize optical depth loop

             NSTART = 1
             NUT_PREV = NSTART - 1

!  Main loop over all output optical depths

             DO UTA = 1, N_LOUTPUT

!  Layer index for given optical depth

              NLEVEL = LEVELMASK_DN(UTA)
              NUT    = NLEVEL

!  Cumulative single scatter source terms :
!      For loop over layers working upwards to level NUT,

              DO N = NSTART, NUT
               NC = N

!  If N = K (profile) or KS = 0 (column), there are extra WFSs
!  If N > K, N < K, transmittance + some multiplication WFSs

!  Elastic terms

               IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                DO W = 1, NPOINTS_INNER
                 DO Q = 1, K_PARAMETERS
                  LSS = ESCAT_DN(N,W)   * L_DN_EMULT(Q,KS,N,W) &
                    + L_ESCAT_DN(Q,N,W) *   DN_EMULT(N,W)
                  L_ESS_CUMSOURCE(Q,W) = LSS &
                  +   DN_LOSTRANS(N,W)   * L_ESS_CUMSOURCE(Q,W) &
                  + L_DN_LOSTRANS(Q,N,W) *   ESS_CUMSOURCE_DN(NC-1,W)
                 ENDDO
                ENDDO
               ELSE
                DO W = 1, NPOINTS_INNER
                 DO Q = 1, K_PARAMETERS
                  LSS = ESCAT_DN(N,W)   * L_DN_EMULT(Q,KS,N,W)
                  L_ESS_CUMSOURCE(Q,W) = LSS &
                  +   DN_LOSTRANS(N,W)   * L_ESS_CUMSOURCE(Q,W)
                 ENDDO
                ENDDO
               ENDIF

!  Raman terms for Inelastic scattering

               IF ( DO_FULL_RAMAN ) THEN
                IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                 DO W = 1, NPOINTS_INNER
                  DO Q = 1, K_PARAMETERS
                   LCSS = CSCAT_DN(N,W)   * L_DN_EMULT(Q,KS,N,W) &
                      + L_CSCAT_DN(Q,N,W) *   DN_EMULT(N,W)
                   LGSS = zero
                   DO B = 1, N_RRSBINS(W)
                    LBSS = GSCAT_DN(N,B,W)   * L_DN_IMULT(Q,KS,N,B,W) &
                       + L_GSCAT_DN(Q,N,B,W) *   DN_IMULT(N,B,W)
                    LGSS = LGSS + LBSS
                   ENDDO
                   LSS = LCSS + LGSS
                   L_ISS_CUMSOURCE(Q,W) = LSS &
                  +   DN_LOSTRANS(N,W)   * L_ISS_CUMSOURCE(Q,W) &
                  + L_DN_LOSTRANS(Q,N,W) *   ISS_CUMSOURCE_DN(NC-1,W)
                  ENDDO
                 ENDDO
                ELSE
                 DO W = 1, NPOINTS_INNER
                  DO Q = 1, K_PARAMETERS
                   LCSS =  CSCAT_DN(N,W) *   L_DN_EMULT(Q,KS,N,W)
                   LGSS = zero
                   DO B = 1, N_RRSBINS(W)
                    LBSS =  GSCAT_DN(N,B,W) * L_DN_IMULT(Q,KS,N,B,W)
                    LGSS = LGSS + LBSS
                   ENDDO
                   LSS = LCSS + LGSS
                   L_ISS_CUMSOURCE(Q,W) = LSS &
                           + DN_LOSTRANS(N,W) * L_ISS_CUMSOURCE(Q,W)
                  ENDDO
                 ENDDO
                ENDIF
               ENDIF

!  End layer source term loop

              ENDDO

!  Ongrid output--------
!     Set final cumulative source and single scatter intensity


!  Profile WFS

              IF ( DO_PROFILE_WFS ) THEN
                DO W = 1, NPOINTS_INNER
                 DO Q = 1, K_PARAMETERS
                  L_ESOURCE = L_ESS_CUMSOURCE(Q,W)
                  L_SSCORRECTION = SSFLUX * L_ESOURCE
                  LP_ELASTIC_SS_DN(Q,K,UTA,V,W) = L_SSCORRECTION
                 ENDDO
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                 DO W = 1, NPOINTS_INNER
                  DO Q = 1, K_PARAMETERS
                   L_ISOURCE = L_ISS_CUMSOURCE(Q,W)
                   L_SSCORRECTION = SSFLUX * L_ISOURCE
                   LP_RAMAN_SS_DN(Q,K,UTA,V,W) = L_SSCORRECTION
                  ENDDO
                 ENDDO
                ENDIF
              ENDIF

!  Column WFS

              IF ( DO_COLUMN_WFS ) THEN
                DO W = 1, NPOINTS_INNER
                  DO Q = 1, K_PARAMETERS
                    L_ESOURCE = L_ESS_CUMSOURCE(Q,W)
                    L_SSCORRECTION = SSFLUX * L_ESOURCE
                    LC_ELASTIC_SS_DN(Q,UTA,V,W) = L_SSCORRECTION
                  ENDDO
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO W = 1, NPOINTS_INNER
                    DO Q = 1, K_PARAMETERS
                      L_ISOURCE = L_ISS_CUMSOURCE(Q,W)
                      L_SSCORRECTION = SSFLUX * L_ISOURCE
                      LC_RAMAN_SS_DN(Q,UTA,V,W) = L_SSCORRECTION
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF

!  Check for updating the recursion

              IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
              NUT_PREV = NUT

!  End optical depth loop

             ENDDO

!  End loop over varying layers

            ENDIF
           ENDDO

!  End Downwelling Jacobian clause

          ENDIF

!  Recurrence relation for the UPWELLING Surface Jacobian
!  ======================================================

!  This contribution is purely Elastic

          IF ( DO_SURFACE_WFS .and. DO_UPWELLING ) THEN

!  Finish if no surface (zero output and go to the next point).

            IF ( .NOT. DO_INCLUDE_SURFACE ) THEN
              DO W = 1, NPOINTS_INNER
                DO UTA = N_LOUTPUT, 1, -1
                  LS_ELASTIC_SS_UP(1:N_SURFACE_WFS,UTA,V,W) = ZERO
                  LS_RAMAN_SS_UP(1:N_SURFACE_WFS,UTA,V,W)   = ZERO
                ENDDO
              ENDDO
            ENDIF

!  If there is a surface.....
!  This section rewritten for Version 2.5
!    Surface Jacobians not normalized, in line with current LIDORT practice....
!mick fix 7/20/2016 - replaced REFLEC with LS_REFLEC in two lines below

            IF ( DO_INCLUDE_SURFACE ) THEN

              IF ( DO_BRDF_SURFACE ) THEN
                IF ( DO_BRDF_Wav1 ) then
                  DO Q = 1, N_SURFACE_WFS
                    Brdf_IDX = 1 ; BRDF_exact = LS_EXACTDB_BRDFUNC(Q,UM,IA,Brdf_IDX)
                    DO W = 1, NPOINTS_INNER
                      LS_REFLEC(Q,W) = BRDF_exact
                    ENDDO
                  ENDDO
                ELSE
                  DO Q = 1, N_SURFACE_WFS
                    DO W = 1, NPOINTS_INNER
                      WR = W + OFFSET_INNER ; Brdf_IDX = WR
                      !REFLEC(W) = LS_EXACTDB_BRDFUNC(Q,UM,IA,Brdf_IDX)
                      LS_REFLEC(Q,W) = LS_EXACTDB_BRDFUNC(Q,UM,IA,Brdf_IDX)
                    ENDDO
                  ENDDO
                ENDIF
              ELSE
                DO Q = 1, N_SURFACE_WFS
                  DO W = 1, NPOINTS_INNER
                    !WR = W + OFFSET_INNER
                    !REFLEC(W) = LS_ALBEDOS_RANKED
                    LS_REFLEC(Q,W) = LS_ALBEDOS_RANKED
                  ENDDO
                ENDDO
              ENDIF

!  Initialize cumulative source terms = F.L(A).mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux.
!    L(A) = Linearized surface term (e.g. albedo)
!
              NC =  0
              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
                BOA_ATTN = X0_BOA_4 * ATTN(NLAYERS,WR)
                DO Q = 1, N_SURFACE_WFS
                  LS_ESS_CUMSOURCE_UP(Q,W) = LS_REFLEC(Q,W) * BOA_ATTN
                ENDDO
              ENDDO

!  Initialize optical depth loop

              NSTART = NLAYERS
              NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

              DO UTA = N_LOUTPUT, 1, -1
                NLEVEL = LEVELMASK_UP(UTA)
                NUT    = NLEVEL + 1

!  Tranmsittance

                DO N = NSTART, NUT, -1
                  !NC = NLAYERS + 1 - N
                  DO W = 1, NPOINTS_INNER
                    DO Q = 1, N_SURFACE_WFS
                      LS_ESS_CUMSOURCE_UP(Q,W) = UP_LOSTRANS(N,W) * LS_ESS_CUMSOURCE_UP(Q,W)
                    ENDDO
                  ENDDO
                ENDDO

!  Set output

                DO W = 1, NPOINTS_INNER
                  DO Q = 1, N_SURFACE_WFS
                    SOURCE = LS_ESS_CUMSOURCE_UP(Q,W)
                    LS_ELASTIC_SS_UP(Q,UTA,V,W) = SSFLUX * SOURCE
                    IF ( DO_FULL_RAMAN ) LS_RAMAN_SS_UP(Q,UTA,V,W) = SSFLUX * SOURCE
                  ENDDO
                ENDDO

!  Check for updating the recursion

                IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                NUT_PREV = NUT

!  End optical depth loop

              ENDDO

!  End inclusion of surface

            ENDIF

!  End Upwelling clause (surface Jacobian)

          ENDIF

!mick add 11/14/2015 - Sleave Jac section

!  Recurrence relation for the UPWELLING Sleave Jacobian
!  =====================================================

         IF ( DO_SLEAVE_WFS .and. DO_UPWELLING ) THEN

!  Finish if no surface (Zero output and go to the next point).

            IF ( .NOT. DO_INCLUDE_SURFACE ) THEN
              DO W = 1, NPOINTS_INNER
                DO UTA = N_LOUTPUT, 1, -1
                  DO Q = 1, N_SLEAVE_WFS
                    Q1 = Q + N_SURFACE_WFS
                    LS_ELASTIC_SS_UP(Q1,UTA,V,W) = ZERO
                    LS_RAMAN_SS_UP(Q1,UTA,V,W)   = ZERO
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  Get the sleaving term linearization

            IF ( DO_INCLUDE_SURFACE ) THEN
              DO Q = 1, N_SLEAVE_WFS
                IF ( DO_SL_ISOTROPIC ) THEN
                  IF ( DO_SLEAVE_Wav1 ) then
                    Sleave_IDX = 1 ; L_Sleave_Iso = LSSL_SLTERM_ISOTROPIC(Q,Sleave_IDX)*PI4
                    DO W = 1, NPOINTS_INNER
                      LSSL_EXACTDB_SOURCE(Q,V,W) = L_Sleave_Iso
                    ENDDO
                  ELSE
                    DO W = 1, NPOINTS_INNER
                      WR = W + OFFSET_INNER ; Sleave_IDX = WR
                      LSSL_EXACTDB_SOURCE(Q,V,W) = LSSL_SLTERM_ISOTROPIC(Q,Sleave_IDX)*PI4
                    ENDDO
                  ENDIF
                ELSE
                  IF ( DO_SLEAVE_Wav1 ) then
                    Sleave_IDX = 1 ; L_Sleave_User = LSSL_SLTERM_USERANGLES(Q,UM,IA,Sleave_IDX)*PI4
                    DO W = 1, NPOINTS_INNER
                      LSSL_EXACTDB_SOURCE(Q,V,W) = L_Sleave_User
                    ENDDO
                  ELSE
                    DO W = 1, NPOINTS_INNER
                      WR = W + OFFSET_INNER ; Sleave_IDX = WR
                      LSSL_EXACTDB_SOURCE(Q,V,W) = LSSL_SLTERM_USERANGLES(Q,UM,IA,Sleave_IDX)*PI4
                    ENDDO
                  ENDIF
                ENDIF
              ENDDO

!  ====================================================
!  2. TRANSMITTANCE OF SLEAVE TERM: UPWELLING RECURSION
!  ====================================================

              DO Q = 1, N_SLEAVE_WFS

!  Offset

                Q1 = Q + N_SURFACE_WFS

!  Initialize cumulative source term

                NC =  0
                DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER
                  !orig LIDORT code line: 
                  !LSSL_DB_CUMSOURCE(V) = LSSL_EXACTDB_SOURCE(Q,V)*FLUXMULT                      !FLUXMULT=FLUX_FACTOR / PI4
                  LSSL_DB_CUMSOURCE(V,W) = LSSL_EXACTDB_SOURCE(Q,V,W)*FLUXES_RANKED(WR)*SSFLUX   !SSFLUX=1.0/PI4
                ENDDO

!  Initialize optical depth loop

                NSTART = NLAYERS
                NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

                DO UTA = N_LOUTPUT, 1, -1

!  Layer index for given optical depth

                  NLEVEL = LEVELMASK_UP(UTA)
                  NUT    = NLEVEL + 1

!  Cumulative layer transmittance :
!    Loop over layers working upwards to level nut

                  DO N = NSTART, NUT, -1
                    !NC = NLAYERS + 1 - N
                    DO W = 1, NPOINTS_INNER
                      LSSL_DB_CUMSOURCE(V,W) = UP_LOSTRANS(N,W) * LSSL_DB_CUMSOURCE(V,W)
                    ENDDO
                  ENDDO

!  Set output
                  DO W = 1, NPOINTS_INNER
                    LS_ELASTIC_SS_UP(Q1,UTA,V,W) = LSSL_DB_CUMSOURCE(V,W)
                    IF ( DO_FULL_RAMAN ) LS_RAMAN_SS_UP(Q1,UTA,V,W) = LSSL_DB_CUMSOURCE(V,W)
                  ENDDO

!  Check for updating the recursion

                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                  NUT_PREV = NUT

!  End optical depth loop

                enddo

!  End loop over all sleave weighting functions

              ENDDO

!  End inclusion of surface

            ENDIF

!  End Upwelling clause (Sleave Jacobian)

         ENDIF 

!  Finish geometry loops (note: these loops start @ ~line 640)

        ENDDO
      ENDDO

!stop'corrplus'

!  Finish

      RETURN
      END SUBROUTINE LRRS_SSCORR_GENERAL_BIN_PLUS_1

!

      SUBROUTINE LRRS_SSCORR_GENERAL_MONO_PLUS_1 &
         ( DO_FULL_RAMAN, DO_UPWELLING, DO_DNWELLING,                     & ! Inputs
           DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,                           & ! Inputs
           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_DMSCALING,          & ! Inputs
           DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,          & ! Inputs
           DO_BRDF_Wav1, DO_Sleave_Wav1,                                  & ! Inputs
           DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,  & ! Inputs
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMNWFS,               & ! Inputs
           NKSTORAGE2, N_SURFACE_WFS, N_SLEAVE_WFS,                       & ! Inputs
           NLAYERS, NFINELAYERS, N_USER_STREAMS, N_USER_RELAZMS,          & ! Input Control Integers
           N_LOUTPUT, LEVELMASK_UP, LEVELMASK_DN,                         & ! Input levels for output
           NPOINTS_MONO, W_EXCIT,                                         & ! Inputs Raman Spec.
           SOLAR_ANGLE, COSSCAT_NADIR_UP, COSSCAT_NADIR_DN,               & ! Input Geometry
           SOLAR_ANGLES_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST,  & ! Input Geometry
           FLUXES_RANKED, ALBEDOS_RANKED, EARTH_RADIUS, HEIGHT_GRID,      & ! Input Fluxes
           EXACTDB_BRDFUNC, SLTERM_ISOTROPIC, SLTERM_USERANGLES,          & ! Input   Surface
           TRUNC_FACTORS, DELTAU_VERT_INPUT, OMEGA_UNSCALED,              & ! Input Optical
           OMEGAPHASFUNC_ELASTIC_UP,   OMEGAPHASFUNC_ELASTIC_DN,          & ! Input Optical Elastic
           OMEGAPHASFUNC_CABANNES_UP,  OMEGAPHASFUNC_CABANNES_DN,         & ! Input Optical Cabannes
           OMEGAMOMS_RRSLOSS_UNSCALED, OMEGAMOMS_RRSGAIN_UNSCALED,        & ! Input Optical Raman
           SAVE_TRANS_USERM, SAVE_BEAMMULT_UP, SAVE_BEAMMULT_DN,          & ! Input Optical for NADIR SSCORR
           LS_EXACTDB_BRDFUNC, LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES, & ! Input Lin Surface
           L_DELTAU_VERT_INPUT, L_TRUNC_FACTORS, L_OMEGA_UNSCALED,            & ! Input Lin Optical
           L_OMEGAPHASFUNC_ELASTIC_UP,   L_OMEGAPHASFUNC_ELASTIC_DN,          & ! Input Lin Optical Elastic
           L_OMEGAPHASFUNC_CABANNES_UP,  L_OMEGAPHASFUNC_CABANNES_DN,         & ! Input Lin Optical Cabannes
           L_OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSGAIN_UNSCALED,        & ! Input Lin Optical Raman
           L_SAVE_TRANS_USERM, L_SAVE_BEAMMULT_UP, L_SAVE_BEAMMULT_DN,    & ! Input Optical for NADIR SSCORR 2p5a
           ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN,        & ! Outputs
           LC_ELASTIC_SS_UP, LP_ELASTIC_SS_UP, LS_ELASTIC_SS_UP,          & ! Outputs
           LC_ELASTIC_SS_DN, LP_ELASTIC_SS_DN,                            & ! Outputs
           LC_RAMAN_SS_UP,   LP_RAMAN_SS_UP,   LS_RAMAN_SS_UP,            & ! Outputs
           LC_RAMAN_SS_DN,   LP_RAMAN_SS_DN,                              & ! Outputs
           FAIL, MESSAGE)                                                   ! Outputs

!  Single scatter exact calculation for the outgoing LOS

!   Programmed by R. Spurr, RT Solutions Inc.
!    First Draft,  for LRRS Version 2.1, April 23rd 2008
!    Second draft, for LRRS Version 2.2, 23 July 2009

!   -- Rob mod 5/12/17 for 2p5a, Added DO_SSCORR_NADIR flag, plus precalculated multipliers/transmittances
!   -- Rob mod 5/12/17 for 2p5a, Using PHASFUNC input now, removed NMOMENTS_INPUT.

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, FOUR, PI4, DEG_TO_RAD,                &
                              MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES,   &
                              MAX_LAYERS, MAX_FINE_LAYERS, MAX_LAYERS_NK,           &
                              MAX_LOUTPUT, MAX_POINTS, MAX_BINS, MAX_ATMOSWFS,      &
                              MAX_SURFACEWFS, MAX_SLEAVEWFS

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

!  Profile WFS control inputs

      LOGICAL  , INTENT(IN) :: DO_PROFILE_WFS
      LOGICAL  , INTENT(IN) :: LAYER_VARY_FLAG   (MAX_LAYERS)
      INTEGER  , INTENT(IN) :: LAYER_VARY_NUMBER (MAX_LAYERS)

!  Linearization bookkeeping

      INTEGER  , INTENT(IN) :: NKSTORAGE2(MAX_LAYERS,0:MAX_LAYERS)

!  Column WFS control inputs

      LOGICAL  , INTENT(IN) :: DO_COLUMN_WFS
      INTEGER  , INTENT(IN) :: N_COLUMNWFS

!  Surface WFS

      LOGICAL  , INTENT(IN) :: DO_SURFACE_WFS
      INTEGER  , INTENT(IN) :: N_SURFACE_WFS

!  Surface-leaving WFS

      LOGICAL  , INTENT(IN) :: DO_SLEAVE_WFS
      INTEGER  , INTENT(IN) :: N_SLEAVE_WFS

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

!  solar fluxes

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

!  WFSs of all optical properties
!  ------------------------------

!  Atmos

      REAL(FPK), INTENT(IN) :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_TRUNC_FACTORS     ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGA_UNSCALED    ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_UP  ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_ELASTIC_DN  ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_CABANNES_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAPHASFUNC_CABANNES_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSLOSS_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSGAIN_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Surface BRDF

      REAL(FPK), INTENT(IN) :: LS_EXACTDB_BRDFUNC ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

!  Surface-leaving

      REAL(fpk), INTENT(IN) :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAX_POINTS )
      REAL(fpk), INTENT(IN) :: LSSL_SLTERM_USERANGLES ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

!  NADIR correction - Linearized user stream transmittances, whole layers
!   -- Rob mod 5/12/17 for 2p5a

      REAL(FPK), INTENT(IN) :: L_SAVE_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  NADIR correction - Linearized pre-calculated multipliers
!   -- Rob mod 5/12/17 for 2p5a

      REAL(FPK), INTENT(IN) :: L_SAVE_BEAMMULT_UP ( MAX_ATMOSWFS, MAX_USER_STREAMS, 0:MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: L_SAVE_BEAMMULT_DN ( MAX_ATMOSWFS, MAX_USER_STREAMS, 0:MAX_LAYERS_NK, MAX_POINTS )

!  Output
!  ------

!  single scatter Elastic and Raman Intensity results

      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_UP (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_DN (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_UP   (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_DN   (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  single scatter Elastic Jacobian results

      REAL(FPK), INTENT(INOUT) :: LC_ELASTIC_SS_UP (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: LC_ELASTIC_SS_DN (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LP_ELASTIC_SS_UP (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: LP_ELASTIC_SS_DN (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LS_ELASTIC_SS_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  single scatter Raman Jacobian results

      REAL(FPK), INTENT(INOUT) :: LC_RAMAN_SS_UP (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: LC_RAMAN_SS_DN (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LP_RAMAN_SS_UP (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(INOUT) :: LP_RAMAN_SS_DN (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LS_RAMAN_SS_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

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

      REAL(FPK) :: extinction   ( max_layers, max_points )
      REAL(FPK) :: L_extinction ( MAX_ATMOSWFS, max_layers, max_points )

!  Saved Legendre polynomials
!   -- Rob mod 5/12/17 for 2p5a, Only need 0:2 now, drop MAX_MOMENTS_INPUT

      REAL(FPK) :: SS_PLEG_UP(MAX_LAYERS,0:2), SS_PLEG_DN(MAX_LAYERS,0:2)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(FPK) :: TMS(MAX_LAYERS)

!  Linearized TMS factor

      REAL(FPK) :: L_TMS(MAX_ATMOSWFS,MAX_LAYERS)

!  Local truncation factors for additional DELTAM scaling

!      REAL(FPK) :: SSFDEL ( MAXLAYERS )

!  Exact Phase function calculations

      REAL(FPK) :: ESCAT_UP(MAX_LAYERS)
      REAL(FPK) :: ESCAT_DN(MAX_LAYERS)
      REAL(FPK) :: CSCAT_UP(MAX_LAYERS)
      REAL(FPK) :: CSCAT_DN(MAX_LAYERS)
      REAL(FPK) :: GSCAT_UP(MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: GSCAT_DN(MAX_LAYERS,MAX_POINTS)

!  Linearized Exact Phase function calculations

      REAL(FPK) :: L_ESCAT_UP(MAX_ATMOSWFS,MAX_LAYERS)
      REAL(FPK) :: L_ESCAT_DN(MAX_ATMOSWFS,MAX_LAYERS)
      REAL(FPK) :: L_CSCAT_UP(MAX_ATMOSWFS,MAX_LAYERS)
      REAL(FPK) :: L_CSCAT_DN(MAX_ATMOSWFS,MAX_LAYERS)
      REAL(FPK) :: L_GSCAT_UP(MAX_ATMOSWFS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_GSCAT_DN(MAX_ATMOSWFS,MAX_LAYERS,MAX_POINTS)

!  Cumulative single scatter source terms

      REAL(FPK) :: ESS_CUMSOURCE_UP(0:MAX_LAYERS)
      REAL(FPK) :: ESS_CUMSOURCE_DN(0:MAX_LAYERS)
      REAL(FPK) :: ISS_CUMSOURCE_UP(0:MAX_LAYERS)
      REAL(FPK) :: ISS_CUMSOURCE_DN(0:MAX_LAYERS)

!  Linearized Cumulative single scatter source terms

      REAL(FPK) :: L_ESS_CUMSOURCE(MAX_ATMOSWFS)
      REAL(FPK) :: L_ISS_CUMSOURCE(MAX_ATMOSWFS)

!  Linearized Cumulative DB source terms

      REAL(FPK) :: LS_ESS_CUMSOURCE_UP(MAX_SURFACEWFS)

!  Outgoing sphericity stuff
!  -------------------------

!  Whole layer Elastic and Inelastic multipliers

      REAL(FPK) :: &
         UP_IMULT(MAX_LAYERS,MAX_POINTS), &
         DN_IMULT(MAX_LAYERS,MAX_POINTS), &
         UP_EMULT(MAX_LAYERS), &
         DN_EMULT(MAX_LAYERS), &
         UP_LOSTRANS(MAX_LAYERS), &
         DN_LOSTRANS(MAX_LAYERS)

!  Linearized Whole layer Elastic and Inelastic multipliers

      REAL(FPK) :: L_UP_IMULT (MAX_ATMOSWFS,0:MAX_LAYERS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_DN_IMULT (MAX_ATMOSWFS,0:MAX_LAYERS,MAX_LAYERS,MAX_POINTS)

      REAL(FPK) :: L_UP_EMULT (MAX_ATMOSWFS,0:MAX_LAYERS,MAX_LAYERS)
      REAL(FPK) :: L_DN_EMULT (MAX_ATMOSWFS,0:MAX_LAYERS,MAX_LAYERS)

      REAL(FPK) :: &
         L_UP_LOSTRANS(MAX_ATMOSWFS,MAX_LAYERS), &
         L_DN_LOSTRANS(MAX_ATMOSWFS,MAX_LAYERS)

!  Solar beam attenuations

      REAL(FPK) :: attn      ( 0:max_layers, max_points )
      REAL(FPK) :: attn_fine ( max_layers, max_fine_layers, max_points )

      REAL(FPK) :: l_attn      ( MAX_ATMOSWFS, 0:max_layers, 0:max_layers, max_points )
      REAL(FPK) :: l_attn_fine ( MAX_ATMOSWFS, 0:max_layers, max_layers, max_fine_layers, max_points )

!  local variables
!  ---------------

!  Indices

      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, W
      INTEGER   :: UTA, UM, IA, NC, V, Q, Q1, K, KS, NK
      INTEGER   :: Sleave_IDX, Brdf_IDX

!  help variables (REAL(FPK) ::)
!   -- Rob mod 5/12/17 for 2p5a, removed DFL1, DFL2

      REAL(FPK) :: boa_attn, ssflux, x0_boa, x0_boa_4, ratio
      REAL(FPK) :: GSOURCE, CSOURCE, ESOURCE, ISOURCE, SOURCE
      REAL(FPK) :: HELP, COSSCAT, CONTRIB, FACTOR, OMEGA_LOCAL, T0, T2, P2
      REAL(FPK) :: PHAS_E, PHAS_LR, PHAS_LG, PHAS_CB
      REAL(FPK) :: REFLEC, LS_REFLEC(MAX_SURFACEWFS)

!  WFS help variables

      REAL(FPK) :: L_OF, L_LOSS, LSS, LCSS, LGSS, LBSS
      REAL(FPK) :: L_ESOURCE, L_ISOURCE, L_SSCORRECTION, L_OMEGA_LOCAL
      REAL(FPK) :: L_PHAS_E, L_PHAS_LR, L_PHAS_LG, L_PHAS_CB
      REAL(FPK) :: L_BOA_ATTN, L_CSOURCE, L_GSOURCE

!  Bookkeeping

      LOGICAL   :: do_atmos_wfs
      INTEGER   :: NTYPE_VARY, K_PARAMETERS
      LOGICAL   :: DO_RTSOL_VARY(max_layers)
      INTEGER   :: NPARAMS_VARY(max_layers)
      INTEGER   :: KINDEX(max_layers)

!  Sleave variables

      REAL(fpk)  :: LSSL_EXACTDB_SOURCE ( MAX_SLEAVEWFS, MAX_GEOMETRIES )
      REAL(fpk)  :: LSSL_DB_CUMSOURCE   ( MAX_GEOMETRIES )

!  Local surface flag

      LOGICAL  , PARAMETER :: DO_INCLUDE_SURFACE = .true.

!  Local value

      REAL(FPK), PARAMETER :: ls_albedos_ranked = one

!  Set up operations
!  -----------------

!  Ss flux scaling factor

      SSFLUX = one / PI4

!  Fine layering flag

      do_fine = .true.

!  Number of fine layers now a control input (7/25/08)

!  Floating point numbers for Legendre polynomials
!   -- Rob mod 5/12/17 for 2p5a, removed DFL1, DFL2
!      DO L = 2, NMOMENTS_INPUT
!        HELP = DBLE(L) ; DF1(L) = DBLE(2*L-1)/HELP ; DF2(L) = DBLE(L-1)/HELP
!      ENDDO

!  General atmospheric WFS term

      do_atmos_wfs = ( do_profile_wfs .or. do_column_wfs )

!  Parameter control
!  -----------------

      IF ( DO_PROFILE_WFS ) THEN
        NTYPE_VARY = NLAYERS
        DO N = 1, NLAYERS
         DO_RTSOL_VARY(n) = LAYER_VARY_FLAG(n)
         NPARAMS_VARY(n)  = LAYER_VARY_NUMBER(n)
         KINDEX(N) = N
        ENDDO
      ELSE IF ( DO_COLUMN_WFS ) THEN
        NTYPE_VARY = 1
        DO_RTSOL_VARY(1) = .TRUE.
        NPARAMS_VARY(1)  = N_COLUMNWFS
        KINDEX(1) = 0
      ENDIF

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

!  Create Linearized TMS factor for each layer
!   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)
!      SHOULD ONLY BE USED WITH THE ENERGY BALANCING
!   -- Rob mod 5/12/17 for 2p5a, Use L_OMEGA_UNSCALED and OMEGA_UNSCALED

      IF ( DO_ATMOS_WFS ) THEN
        DO K = 1, NLAYERS
          K_PARAMETERS = LAYER_VARY_NUMBER(K)
          DO Q = 1, K_PARAMETERS
            IF ( DO_DMSCALING ) THEN
              OMEGA_LOCAL    = OMEGA_UNSCALED(K,W_EXCIT)
              L_OMEGA_LOCAL  = L_OMEGA_UNSCALED(Q,K,W_EXCIT)
              L_OF = L_OMEGA_LOCAL * TRUNC_FACTORS(K,W_EXCIT) + L_TRUNC_FACTORS(Q,K,W_EXCIT) * OMEGA_LOCAL
              L_TMS(Q,K) = TMS(K) * TMS(K) * L_OF
            ELSE
              L_TMS(Q,K) = zero
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  Create extinctions (using scaled optical depths)

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        DO W = 1, NPOINTS_MONO
          EXTINCTION(N,W) = DELTAU_VERT_INPUT(N,W) / HELP
        ENDDO
      ENDDO

!  Linearized extinctions

      IF ( DO_ATMOS_WFS ) THEN
       DO K = 1, NLAYERS
         HELP = HEIGHT_GRID(K-1) - HEIGHT_GRID(K)
         K_PARAMETERS = LAYER_VARY_NUMBER(K)
         DO Q = 1, K_PARAMETERS
           DO W = 1, NPOINTS_MONO
             L_EXTINCTION(Q,K,W) = L_DELTAU_VERT_INPUT(Q,K,W) / HELP
           ENDDO
         ENDDO
       ENDDO
      ENDIF

!  Additional Delta-M scaling
!  --------------------------

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

!  Multipliers

          DO N = 1, NLAYERS
            UP_LOSTRANS(N) = SAVE_TRANS_USERM(UM,N,W_EXCIT)
            UP_EMULT(N)    = SAVE_BEAMMULT_UP(UM,N,W_EXCIT)
            DO W = 1, NPOINTS_MONO
              UP_IMULT(N,W) = SAVE_BEAMMULT_UP(UM,N,W)
            ENDDO
          ENDDO

!  Linearized Multipliers

          DO N = 1, NLAYERS
            DO K = 1, NTYPE_VARY
              IF ( DO_RTSOL_VARY(K) ) THEN
                K_PARAMETERS = NPARAMS_VARY(K) ; KS = KINDEX(K) ; NK = NKSTORAGE2(N,KS)
                IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                  DO Q = 1, K_PARAMETERS
                    L_UP_LOSTRANS(Q,N) = L_SAVE_TRANS_USERM(Q,UM,N,W_EXCIT)
                  ENDDO
                ENDIF
                DO Q = 1, K_PARAMETERS
                  L_UP_EMULT(Q,KS,N) = L_SAVE_BEAMMULT_UP(Q,UM,NK,W_EXCIT)
                  DO W = 1, NPOINTS_MONO
                     L_UP_IMULT(Q,KS,N,W) = L_SAVE_BEAMMULT_UP(Q,UM,NK,W)
                   ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDDO

!  End upwelling Nadir multipliers

        ENDIF

!  DOWNWELLING NADIR CORRECTION: multipliers, transmittances
!  =========================================================

!  Downwelling Multipliers and Transmittances for the SSCORR_NADIR case
!     Does not use the ADJUSTED geometrical angles, Same for each azimuth.
!   -- Rob mod 5/12/17 for 2p5a, SSCORR_NADIR values are already pre-calculated

        IF ( DO_SSCORR_NADIR .and.DO_DNWELLING ) THEN

!  Multipliers

          DO N = 1, NLAYERS
            DN_LOSTRANS(N) = SAVE_TRANS_USERM(UM,N,W_EXCIT)
            DN_EMULT(N) = SAVE_BEAMMULT_DN(UM,N,W_EXCIT)
            DO W = 1, NPOINTS_MONO
              DN_IMULT(N,W) = SAVE_BEAMMULT_DN(UM,N,W)
            ENDDO
          ENDDO

!  Linearized Multipliers

          DO N = 1, NLAYERS
            DO K = 1, NTYPE_VARY
              IF ( DO_RTSOL_VARY(K) ) THEN
                K_PARAMETERS = NPARAMS_VARY(K) ; KS = KINDEX(K) ; NK = NKSTORAGE2(N,KS)
                IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                  DO Q = 1, K_PARAMETERS
                    L_DN_LOSTRANS(Q,N) = L_SAVE_TRANS_USERM(Q,UM,N,W_EXCIT)
                  ENDDO
                ENDIF
                DO Q = 1, K_PARAMETERS
                  L_DN_EMULT(Q,KS,N) = L_SAVE_BEAMMULT_DN(Q,UM,NK,W_EXCIT)
                  DO W = 1, NPOINTS_MONO
                     L_DN_IMULT(Q,KS,N,W) = L_SAVE_BEAMMULT_DN(Q,UM,NK,W)
                   ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDDO

!  End Downwelling Nadir multipliers

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

            call og_attenuations_plus_1 &
        ( max_layers, max_fine_layers, max_points, MAX_ATMOSWFS, & ! Inputs
          npoints_mono, nlayers, nfinelayers, do_profile_wfs,    & ! Inputs
          layer_vary_flag, layer_vary_number, do_column_wfs,     & ! Inputs
          n_columnwfs, extinction, L_extinction,                 & ! Inputs
          sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,    & ! Inputs
          attn, attn_fine, l_attn, l_attn_fine )                   ! Outputs

!  Multipliers, transmittances

            call og_integration_mono_up_plus_1 &
        ( max_layers, max_fine_layers, max_points, max_atmoswfs,       & ! Inputs
          do_full_raman,  npoints_mono, w_excit, nlayers, nfinelayers, & ! Inputs
          do_profile_wfs, layer_vary_flag, layer_vary_number,          & ! Inputs
          do_column_wfs,  n_columnwfs,                                 & ! Inputs
          extinction, l_extinction, radii, alpha_all, alpha_fine,      & ! Inputs
          attn, attn_fine, l_attn, l_attn_fine,                        & ! Inputs
          up_emult, up_imult, up_lostrans,                             & ! Inputs
          l_up_emult, l_up_imult, l_up_lostrans )                        ! Outputs

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
            P2 = SS_PLEG_UP(1,2)

!  Elastic scatting Phase functions (multiplied by TMS factor)
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_E with Phase function input
!   -- PHAS_E   = DOT_PRODUCT(OMEGAMOMS_ELASTIC_UNSCALED(N,0:NMOMENTS_INPUT,W1),SS_PLEG_UP(1,0:NMOMENTS_INPUT))
!   -- L_PHAS_E = DOT_PRODUCT(L_OMEGAMOMS_ELASTIC_UNSCALED(Q,N,0:NMOMENTS_INPUT,W1),SS_PLEG_UP(1,0:NMOMENTS_INPUT))

            DO N = 1, NLAYERS
              PHAS_E = OMEGAPHASFUNC_ELASTIC_UP(N,V,W_EXCIT)
              ESCAT_UP(N) = PHAS_E * TMS(N)
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_PHAS_E =  L_OMEGAPHASFUNC_ELASTIC_UP(Q,N,V,W_EXCIT)
                  L_ESCAT_UP(Q,N) = L_PHAS_E * TMS(N) + PHAS_E * L_TMS(Q,N)
                ENDDO
              ENDIF
            ENDDO

!  Add inelastic scatting Phase functions if flagged

            IF ( DO_FULL_RAMAN ) THEN

!  Energy-balancing method, CSCAT_UP = Elastic - Loss
! ---------------------------------------------------

              IF ( DO_ENERGY_BALANCING ) THEN
                DO N = 1, NLAYERS
                  PHAS_LR = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,1) + &
                            OMEGAMOMS_RRSLOSS_UNSCALED(N,2,1) * P2
                  CSCAT_UP(N) = ESCAT_UP(N) - PHAS_LR * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      L_PHAS_LR = L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,0,1) + &
                                  L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,2,1) * P2
                      L_LOSS = PHAS_LR * L_TMS(Q,N) + L_PHAS_LR * TMS(N)
                      L_CSCAT_UP(Q,N) = L_ESCAT_UP(Q,N) - L_LOSS
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF

!  Cabannes-Raman method, CSCAT_UP = Cabannes (direct)
!  ---------------------------------------------------

!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_CB with Phase function input
!   -- PHAS_CB = DOT_PRODUCT(OMEGAMOMS_CABANNES_UNSCALED(N,0:NMOMENTS_INPUT,W_EXCIT),SS_PLEG_UP(1,0:NMOMENTS_INPUT))

              IF ( DO_CABANNES_RAMAN ) THEN
                DO N = 1, NLAYERS
                  PHAS_CB = OMEGAPHASFUNC_CABANNES_UP(N,V,1)
                  CSCAT_UP(N) = PHAS_CB * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      L_PHAS_CB = L_OMEGAPHASFUNC_CABANNES_UP(Q,N,V,1)
                      L_CSCAT_UP(Q,N) = L_PHAS_CB * TMS(N) + PHAS_CB * L_TMS(Q,N)
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF

!  Gain term (Same for both realizations)
!  --------------------------------------

              DO N = 1, NLAYERS
                DO W = 1, NPOINTS_MONO
                  RATIO   = FLUXES_RANKED(W)/FLUXES_RANKED(W_EXCIT)
                  CONTRIB = OMEGAMOMS_RRSGAIN_UNSCALED(N,0,W) + &
                            OMEGAMOMS_RRSGAIN_UNSCALED(N,2,W) * P2
                  PHAS_LG  =  CONTRIB * RATIO
                  GSCAT_UP(N,W) = PHAS_LG * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      T0 = L_OMEGAMOMS_RRSGAIN_UNSCALED(Q,N,0,W)
                      T2 = L_OMEGAMOMS_RRSGAIN_UNSCALED(Q,N,2,W) * P2
                      L_PHAS_LG = ( T0 + T2 ) * RATIO
                      L_GSCAT_UP(Q,N,W) = L_PHAS_LG * TMS(N) + PHAS_LG * L_TMS(Q,N)
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO

!  End inelastic scattering option

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

            call og_attenuations_plus_1 &
        ( max_layers, max_fine_layers, max_points, MAX_ATMOSWFS, & ! Inputs
          npoints_mono, nlayers, nfinelayers, do_profile_wfs,    & ! Inputs
          layer_vary_flag, layer_vary_number, do_column_wfs,     & ! Inputs
          N_COLUMNWFS, extinction, L_extinction,                 & ! Inputs
          sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,    & ! Inputs
          attn, attn_fine, l_attn, l_attn_fine )                   ! Outputs

!  Multipliers, transmittances

            call og_integration_mono_dn_plus_1 &
        ( max_layers, max_fine_layers, max_points, MAX_ATMOSWFS,       & ! Inputs
          do_full_raman,  npoints_mono, w_excit, nlayers, nfinelayers, & ! Inputs
          do_profile_wfs, layer_vary_flag, layer_vary_number,          & ! Inputs
          do_column_wfs,  N_COLUMNWFS,                                 & ! Inputs
          extinction, l_extinction, radii, alpha_all, alpha_fine,      & ! Inputs
          attn, attn_fine, l_attn, l_attn_fine,                        & ! Inputs
          dn_emult, dn_imult, dn_lostrans,                             & ! Inputs
          l_dn_emult, l_dn_imult, l_dn_lostrans )

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
            P2 = SS_PLEG_DN(1,2)

!  Elastic scatting Phase functions (multiplied by TMS factor)
!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_E with Phase function input
!   -- PHAS_E   = DOT_PRODUCT(OMEGAMOMS_ELASTIC_UNSCALED(N,0:NMOMENTS_INPUT,W1),SS_PLEG_UP(1,0:NMOMENTS_INPUT))
!   -- L_PHAS_E = DOT_PRODUCT(L_OMEGAMOMS_ELASTIC_UNSCALED(Q,N,0:NMOMENTS_INPUT,W1),SS_PLEG_UP(1,0:NMOMENTS_INPUT))

            DO N = 1, NLAYERS
              PHAS_E = OMEGAPHASFUNC_ELASTIC_DN(N,V,W_EXCIT)
              ESCAT_DN(N) = PHAS_E * TMS(N)
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_PHAS_E =  L_OMEGAPHASFUNC_ELASTIC_DN(Q,N,V,W_EXCIT)
                  L_ESCAT_DN(Q,N) = L_PHAS_E * TMS(N) + PHAS_E * L_TMS(Q,N)
                ENDDO
              ENDIF
            ENDDO

!  Add inelastic scatting Phase functions if flagged
!  -------------------------------------------------

            IF ( DO_FULL_RAMAN ) THEN

!  Energy-balancing method, CSCAT_DN = Elastic - Loss
! ---------------------------------------------------

              IF ( DO_ENERGY_BALANCING ) THEN
                DO N = 1, NLAYERS
                  PHAS_LR = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,1) + &
                            OMEGAMOMS_RRSLOSS_UNSCALED(N,2,1) * P2
                  CSCAT_DN(N) = ESCAT_DN(N) - PHAS_LR * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      L_PHAS_LR = L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,0,1) + &
                                  L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,2,1) * P2
                      L_LOSS = PHAS_LR * L_TMS(Q,N) + L_PHAS_LR * TMS(N)
                      L_CSCAT_DN(Q,N) = L_ESCAT_DN(Q,N) - L_LOSS
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF

!  Cabannes-Raman method, CSCAT_UP = Cabannes (direct)
!  ---------------------------------------------------

!   -- Rob mod 5/12/17 for 2p5a, replace PHAS_CB with Phase function input
!   -- PHAS_CB = DOT_PRODUCT(OMEGAMOMS_CABANNES_UNSCALED(N,0:NMOMENTS_INPUT,W_EXCIT),SS_PLEG_UP(1,0:NMOMENTS_INPUT))

              IF ( DO_CABANNES_RAMAN ) THEN
                DO N = 1, NLAYERS
                  PHAS_CB = OMEGAPHASFUNC_CABANNES_DN(N,V,1)
                  CSCAT_DN(N) = PHAS_CB * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      L_PHAS_CB = L_OMEGAPHASFUNC_CABANNES_DN(Q,N,V,1)
                      L_CSCAT_DN(Q,N) = L_PHAS_CB * TMS(N) + PHAS_CB * L_TMS(Q,N)
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF

!  Gain term (Same for both realizations)
!  --------------------------------------

              DO N = 1, NLAYERS
                DO W = 1, NPOINTS_MONO
                  RATIO   = FLUXES_RANKED(W)/FLUXES_RANKED(W_EXCIT)
                  CONTRIB = OMEGAMOMS_RRSGAIN_UNSCALED(N,0,W) + &
                            OMEGAMOMS_RRSGAIN_UNSCALED(N,2,W) * P2
                  PHAS_LG  =  CONTRIB * RATIO
                  GSCAT_DN(N,W) = PHAS_LG * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      T0 = L_OMEGAMOMS_RRSGAIN_UNSCALED(Q,N,0,W)
                      T2 = L_OMEGAMOMS_RRSGAIN_UNSCALED(Q,N,2,W) * P2
                      L_PHAS_LG = ( T0 + T2 ) * RATIO
                      L_GSCAT_DN(Q,N,W) = L_PHAS_LG * TMS(N) + PHAS_LG * L_TMS(Q,N)
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO

!  End inelastic scattering option

            ENDIF

!  End Downwelling clause

          ENDIF

!  Recurrence relation for the UPWELLING intensity
!  ===============================================

          IF ( DO_UPWELLING ) THEN

!  This section rewritten for Version 2.5

            IF ( DO_INCLUDE_SURFACE ) THEN
              IF ( DO_BRDF_SURFACE ) THEN
                Brdf_IDX = 1 ; if ( DO_BRDF_Wav1 ) Brdf_IDX = W_EXCIT
                REFLEC = EXACTDB_BRDFUNC(UM,IA,Brdf_IDX)
              ELSE
                REFLEC = ALBEDOS_RANKED(W_EXCIT)
              ENDIF
            ENDIF

!  initialize cumulative source terms = F.A.mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux, A = albedo
!    Only require the total intensity component.
!    This contribution is purely elastic.
!    This line Replaced: X0_BOA = cos(SOLAR_ANGLE * DEG_TO_RAD)
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

!  Duplicate the result for the inelastic case

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
                ESOURCE =  ESCAT_UP(N) * UP_EMULT(N)
                ESS_CUMSOURCE_UP(NC) = ESOURCE + &
                       UP_LOSTRANS(N) * ESS_CUMSOURCE_UP(NC-1)
                IF ( DO_FULL_RAMAN ) THEN
                  CSOURCE = CSCAT_UP(N) * UP_EMULT(N)
                  GSOURCE = ZERO
                  DO W = 1, NPOINTS_MONO
                    GSOURCE = GSOURCE + GSCAT_UP(N,W) * UP_IMULT(N,W)
                  ENDDO
                  ISOURCE = CSOURCE + GSOURCE
                  ISS_CUMSOURCE_UP(NC) = ISOURCE + &
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
                ESOURCE =  ESCAT_DN(N) * DN_EMULT(N)
                ESS_CUMSOURCE_DN(NC) = ESOURCE + &
                       DN_LOSTRANS(N) * ESS_CUMSOURCE_DN(NC-1)
                IF ( DO_FULL_RAMAN ) THEN
                  CSOURCE = CSCAT_DN(N) * DN_EMULT(N)
                  GSOURCE = ZERO
                  DO W = 1, NPOINTS_MONO
                    GSOURCE = GSOURCE + GSCAT_DN(N,W) * DN_IMULT(N,W)
                  ENDDO
                  ISOURCE = CSOURCE + GSOURCE
                  ISS_CUMSOURCE_DN(NC) = ISOURCE + &
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

!  Recurrence relation for the UPWELLING Atmospheric Jacobians
!  ===========================================================

          IF ( DO_UPWELLING .AND. DO_ATMOS_WFS ) THEN

!  Start the main layer variation loop
!  -----------------------------------

           DO K = 1, NTYPE_VARY
            IF ( DO_RTSOL_VARY(K) ) THEN
             K_PARAMETERS = NPARAMS_VARY(K)
             KS = KINDEX(K)

!  initialize cumulative source terms = F.A.mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux, A = albedo
!    Only require the Stokes total intensity component.
!    This contribution is purely elastic.

             IF ( DO_INCLUDE_SURFACE ) THEN
              DO Q = 1, K_PARAMETERS
                FACTOR = REFLEC
                L_BOA_ATTN = X0_BOA_4 * L_ATTN(Q,KS,NLAYERS,W_EXCIT)
                L_ESS_CUMSOURCE(Q) = FACTOR * L_BOA_ATTN
              ENDDO
             ELSE
              DO Q = 1, K_PARAMETERS
                L_ESS_CUMSOURCE(Q) = zero
              ENDDO
             ENDIF

!  Duplicate the result for the inelastic case

             IF ( DO_FULL_RAMAN ) THEN
               DO Q = 1, K_PARAMETERS
                 L_ISS_CUMSOURCE(Q) = L_ESS_CUMSOURCE(Q)
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

              DO N = NSTART, NUT, -1
               NC = NLAYERS + 1 - N

!  Elastic terms

               IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                 DO Q = 1, K_PARAMETERS
                   L_ESOURCE = ESCAT_UP(N)   * L_UP_EMULT(Q,KS,N) &
                           + L_ESCAT_UP(Q,N) *   UP_EMULT(N)
                   L_ESS_CUMSOURCE(Q) = L_ESOURCE &
                   +   UP_LOSTRANS(N)   * L_ESS_CUMSOURCE(Q) &
                   + L_UP_LOSTRANS(Q,N) *   ESS_CUMSOURCE_UP(NC-1)
                 ENDDO
               ELSE
                 DO Q = 1, K_PARAMETERS
                   L_ESOURCE = ESCAT_UP(N) * L_UP_EMULT(Q,KS,N)
                   L_ESS_CUMSOURCE(Q) = L_ESOURCE &
                   +   UP_LOSTRANS(N)   * L_ESS_CUMSOURCE(Q)
                 ENDDO
               ENDIF

!  Raman terms for Inelastic scattering

               IF ( DO_FULL_RAMAN ) THEN
                 IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                   DO Q = 1, K_PARAMETERS
                     L_CSOURCE = CSCAT_UP(N)   * L_UP_EMULT(Q,KS,N) &
                             + L_CSCAT_UP(Q,N) *   UP_EMULT(N)
                     L_GSOURCE = zero
                     DO W = 1, NPOINTS_MONO
                       LBSS = GSCAT_UP(N,W) * L_UP_IMULT(Q,KS,N,W) &
                           + L_GSCAT_UP(Q,N,W) * UP_IMULT(N,W)
                       L_GSOURCE = L_GSOURCE + LBSS
                     ENDDO
                     L_ISOURCE  = L_CSOURCE + L_GSOURCE
                     L_ISS_CUMSOURCE(Q) = L_ISOURCE &
                     +   UP_LOSTRANS(N)   * L_ISS_CUMSOURCE(Q) &
                     + L_UP_LOSTRANS(Q,N) *   ISS_CUMSOURCE_UP(NC-1)
                   ENDDO
                 ELSE
                   DO Q = 1, K_PARAMETERS
                     LCSS =  CSCAT_UP(N) *   L_UP_EMULT(Q,KS,N)
                     LGSS = zero
                     DO W = 1, NPOINTS_MONO
                       LBSS =  GSCAT_UP(N,W) * L_UP_IMULT(Q,KS,N,W)
                       LGSS = LGSS + LBSS
                     ENDDO
                     LSS = LCSS + LGSS
                     L_ISS_CUMSOURCE(Q) = LSS &
                           + UP_LOSTRANS(N) * L_ISS_CUMSOURCE(Q)
                   ENDDO
                 ENDIF
               ENDIF

!  End layer source term loop

              ENDDO

!  Ongrid output--------
!     Set final cumulative source and single scatter intensity

!  If N = K (the layer that is varying)
!    add WFS of additional partial layer source term =
!        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)
!  Variations when N > K
!    add WFS of additional partial layer source term =
!         Exact_Scat * L_Multiplier(k)

!  Profile WFS

              IF ( DO_PROFILE_WFS ) THEN
                DO Q = 1, K_PARAMETERS
                  L_ESOURCE = L_ESS_CUMSOURCE(Q)
                  L_SSCORRECTION = SSFLUX * L_ESOURCE
                  LP_ELASTIC_SS_UP(Q,K,UTA,V,1) = L_SSCORRECTION
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO Q = 1, K_PARAMETERS
                    L_ISOURCE = L_ISS_CUMSOURCE(Q)
                    L_SSCORRECTION = SSFLUX * L_ISOURCE
                    LP_RAMAN_SS_UP(Q,K,UTA,V,1) = L_SSCORRECTION
                  ENDDO
                ENDIF
              ENDIF

!  Column WFS

              IF ( DO_COLUMN_WFS ) THEN
                DO Q = 1, K_PARAMETERS
                  L_ESOURCE = L_ESS_CUMSOURCE(Q)
                  L_SSCORRECTION = SSFLUX * L_ESOURCE
                  LC_ELASTIC_SS_UP(Q,UTA,V,1) = L_SSCORRECTION
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO Q = 1, K_PARAMETERS
                    L_ISOURCE = L_ISS_CUMSOURCE(Q)
                    L_SSCORRECTION = SSFLUX * L_ISOURCE
                    LC_RAMAN_SS_UP(Q,UTA,V,1) = L_SSCORRECTION
                  ENDDO
                ENDIF
              ENDIF

!  Check for updating the recursion

              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
              NUT_PREV = NUT

!  end optical depth loop

             ENDDO

!  end loop over varying layers

            ENDIF
           ENDDO

!  end Upwelling Jacobian clause

          ENDIF

!  Recurrence relation for the DOWNWELLING Atmospheric Jacobians
!  =============================================================

          IF ( DO_DNWELLING .AND. DO_ATMOS_WFS ) THEN

!  Start the main layer variation loop
!  -----------------------------------

           DO K = 1, NTYPE_VARY
            IF ( DO_RTSOL_VARY(K) ) THEN
             K_PARAMETERS = NPARAMS_VARY(K)
             KS = KINDEX(K)

!  Initialise elastic source term

             DO Q = 1, K_PARAMETERS
               L_ESS_CUMSOURCE(Q) = zero
             ENDDO

!  Duplicate the result for the inelastic case

             IF ( DO_FULL_RAMAN ) THEN
               DO Q = 1, K_PARAMETERS
                 L_ISS_CUMSOURCE(Q) = L_ESS_CUMSOURCE(Q)
               ENDDO
             ENDIF

!  initialize optical depth loop

             NSTART = 1
             NUT_PREV = NSTART - 1

!  Main loop over all output optical depths

             DO UTA = 1, N_LOUTPUT

!  Layer index for given optical depth

              NLEVEL = LEVELMASK_DN(UTA)
              NUT    = NLEVEL

!  Cumulative single scatter source terms :
!      For loop over layers working upwards to level NUT,

              DO N = NSTART, NUT
               NC = N

!  If N = K (profile) or KS = 0 (column), there are extra WFSs
!  If N > K, N < K, transmittance + some multiplication WFSs

!  Elastic terms

               IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                 DO Q = 1, K_PARAMETERS
                  LSS = ESCAT_DN(N)   * L_DN_EMULT(Q,KS,N) &
                    + L_ESCAT_DN(Q,N) *   DN_EMULT(N)
                  L_ESS_CUMSOURCE(Q) = LSS &
                  +   DN_LOSTRANS(N)   * L_ESS_CUMSOURCE(Q) &
                  + L_DN_LOSTRANS(Q,N) *   ESS_CUMSOURCE_DN(NC-1)
                 ENDDO
               ELSE
                 DO Q = 1, K_PARAMETERS
                  LSS = ESCAT_DN(N)   * L_DN_EMULT(Q,KS,N)
                  L_ESS_CUMSOURCE(Q) = LSS &
                  +   DN_LOSTRANS(N)   * L_ESS_CUMSOURCE(Q)
                 ENDDO
               ENDIF

!  Raman terms for Inelastic scattering

               IF ( DO_FULL_RAMAN ) THEN
                 IF ( N.EQ.K .OR. KS.EQ.0 ) THEN
                   DO Q = 1, K_PARAMETERS
                     LCSS = CSCAT_DN(N)   * L_DN_EMULT(Q,KS,N) &
                        + L_CSCAT_DN(Q,N) *   DN_EMULT(N)
                     LGSS = zero
                     DO W = 1, NPOINTS_MONO
                       LBSS = GSCAT_DN(N,W)   * L_DN_IMULT(Q,KS,N,W) &
                          + L_GSCAT_DN(Q,N,W) *   DN_IMULT(N,W)
                       LGSS = LGSS + LBSS
                     ENDDO
                     LSS = LCSS + LGSS
                     L_ISS_CUMSOURCE(Q) = LSS &
                    +   DN_LOSTRANS(N)   * L_ISS_CUMSOURCE(Q) &
                    + L_DN_LOSTRANS(Q,N) *   ISS_CUMSOURCE_DN(NC-1)
                   ENDDO
                 ELSE
                   DO Q = 1, K_PARAMETERS
                     LCSS =  CSCAT_DN(N) *   L_DN_EMULT(Q,KS,N)
                     LGSS = zero
                     DO W = 1, NPOINTS_MONO
                       LBSS =  GSCAT_DN(N,W) * L_DN_IMULT(Q,KS,N,W)
                       LGSS = LGSS + LBSS
                     ENDDO
                     LSS = LCSS + LGSS
                     L_ISS_CUMSOURCE(Q) = LSS &
                           + DN_LOSTRANS(N) * L_ISS_CUMSOURCE(Q)
                   ENDDO
                 ENDIF
               ENDIF

!  End layer source term loop

              ENDDO

!  Ongrid output--------
!     Set final cumulative source and single scatter intensity

!  Profile WFS

              IF ( DO_PROFILE_WFS ) THEN
                DO Q = 1, K_PARAMETERS
                  L_ESOURCE = L_ESS_CUMSOURCE(Q)
                  L_SSCORRECTION = SSFLUX * L_ESOURCE
                  LP_ELASTIC_SS_DN(Q,K,UTA,V,1) = L_SSCORRECTION
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO Q = 1, K_PARAMETERS
                    L_ISOURCE = L_ISS_CUMSOURCE(Q)
                    L_SSCORRECTION = SSFLUX * L_ISOURCE
                    LP_RAMAN_SS_DN(Q,K,UTA,V,1) = L_SSCORRECTION
                  ENDDO
                ENDIF
              ENDIF

!  Column WFS

              IF ( DO_COLUMN_WFS ) THEN
                DO Q = 1, K_PARAMETERS
                  L_ESOURCE = L_ESS_CUMSOURCE(Q)
                  L_SSCORRECTION = SSFLUX * L_ESOURCE
                  LC_ELASTIC_SS_DN(Q,UTA,V,1) = L_SSCORRECTION
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO Q = 1, K_PARAMETERS
                    L_ISOURCE = L_ISS_CUMSOURCE(Q)
                    L_SSCORRECTION = SSFLUX * L_ISOURCE
                    LC_RAMAN_SS_DN(Q,UTA,V,1) = L_SSCORRECTION
                  ENDDO
                ENDIF
              ENDIF

!  Check for updating the recursion

              IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
              NUT_PREV = NUT

!  end optical depth loop

             ENDDO

!  end loop over varying layers

            ENDIF
           ENDDO

!  end Downwelling Jacobian clause

          ENDIF

!  Recurrence relation for the UPWELLING Surface Jacobian
!  ======================================================

!  This contribution is purely Elastic

          IF ( DO_SURFACE_WFS .and. DO_UPWELLING ) THEN

!  Finish if no surface (Zero output and go to the next point.

            IF ( .NOT. DO_INCLUDE_SURFACE ) THEN
              DO UTA = N_LOUTPUT, 1, -1
                LS_ELASTIC_SS_UP(1:N_SURFACE_WFS,UTA,V,1) = ZERO
                LS_RAMAN_SS_UP(1:N_SURFACE_WFS,UTA,V,1)   = ZERO
              ENDDO
            ENDIF

!  If there is a surface.....

            IF ( DO_INCLUDE_SURFACE ) THEN

!   This section of code revised for the Version 2.5 supplement-derived BRDF.

              IF ( DO_BRDF_SURFACE ) THEN
                Brdf_IDX = 1 ; if ( DO_BRDF_Wav1 ) Brdf_IDX = W_EXCIT
                LS_REFLEC(1:N_SURFACE_WFS) = LS_EXACTDB_BRDFUNC(1:N_SURFACE_WFS,UM,IA,Brdf_IDX)
              ELSE
                LS_REFLEC(1:N_SURFACE_WFS) = LS_ALBEDOS_RANKED
              ENDIF

!  initialize cumulative source terms = F.L(A).mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux.
!    L(A) = Linearized surface term (e.g. albedo)
!
              NC =  0
              BOA_ATTN = X0_BOA_4 * ATTN(NLAYERS,W_EXCIT)
              LS_ESS_CUMSOURCE_UP(1:N_SURFACE_WFS) = LS_REFLEC(1:N_SURFACE_WFS) * BOA_ATTN

!  initialize optical depth loop

              NSTART = NLAYERS
              NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

              DO UTA = N_LOUTPUT, 1, -1
                NLEVEL = LEVELMASK_UP(UTA)
                NUT    = NLEVEL + 1

!  Cumulative single scatter source terms :
!      For loop over layers working upwards to level NUT.
!  Rob fix 2/7/17. No Transmittance in the original !!!

                DO N = NSTART, NUT, -1
                  NC = NLAYERS + 1 - N
!                  LS_ESS_CUMSOURCE_UP(1:N_SURFACE_WFS) = LS_ESS_CUMSOURCE_UP(1:N_SURFACE_WFS)
                  LS_ESS_CUMSOURCE_UP(1:N_SURFACE_WFS) = UP_LOSTRANS(N) * LS_ESS_CUMSOURCE_UP(1:N_SURFACE_WFS)
                ENDDO

!  Set output

                LS_ELASTIC_SS_UP(1:N_SURFACE_WFS,UTA,V,1) = SSFLUX * LS_ESS_CUMSOURCE_UP(1:N_SURFACE_WFS)
                IF  ( DO_FULL_RAMAN ) LS_RAMAN_SS_UP(1:N_SURFACE_WFS,UTA,V,1) = SSFLUX * LS_ESS_CUMSOURCE_UP(1:N_SURFACE_WFS)

!  Check for updating the recursion

                IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                NUT_PREV = NUT

!  end optical depth loop

              ENDDO

!  No surface clause

            ENDIF

!  end Upwelling clause (surface Jacobian)

         ENDIF

!mick add 11/14/2015 - Sleave Jac section

!  Recurrence relation for the UPWELLING Sleave Jacobian
!  =====================================================

         IF ( DO_SLEAVE_WFS .and. DO_UPWELLING ) THEN

!  Finish if no surface (Zero output and go to the next point.

            IF ( .NOT. DO_INCLUDE_SURFACE ) THEN
              DO UTA = N_LOUTPUT, 1, -1
                DO Q = 1, N_SLEAVE_WFS
                  Q1 = Q + N_SURFACE_WFS
                  LS_ELASTIC_SS_UP(Q1,UTA,V,1) = ZERO
                  LS_RAMAN_SS_UP(Q1,UTA,V,1)   = ZERO
                ENDDO
              ENDDO
            ENDIF

!  Get the sleaving term linearization

            IF ( DO_INCLUDE_SURFACE ) THEN
              DO Q = 1, N_SLEAVE_WFS
                Sleave_IDX = 1 ; if ( DO_Sleave_Wav1 ) Sleave_IDX = W_EXCIT
                IF ( DO_SL_ISOTROPIC ) THEN
                  LSSL_EXACTDB_SOURCE(Q,V) = LSSL_SLTERM_ISOTROPIC(Q,Sleave_IDX)*PI4
                ELSE
                  LSSL_EXACTDB_SOURCE(Q,V) = LSSL_SLTERM_USERANGLES(Q,UM,IA,Sleave_IDX)*PI4
                ENDIF
              ENDDO

!  ====================================================
!  2. TRANSMITTANCE OF SLEAVE TERM: UPWELLING RECURSION
!  ====================================================

              DO Q = 1, N_SLEAVE_WFS

!  Offset

                Q1 = Q + N_SURFACE_WFS

!  Initialize cumulative source term

                NC =  0
                !orig LIDORT code line: 
                !LSSL_DB_CUMSOURCE(V) = LSSL_EXACTDB_SOURCE(Q,V)*FLUXMULT                     !FLUXMULT=FLUX_FACTOR / PI4
                LSSL_DB_CUMSOURCE(V) = LSSL_EXACTDB_SOURCE(Q,V)*FLUXES_RANKED(W_EXCIT)*SSFLUX   !SSFLUX=1.0/PI4

!  Initialize optical depth loop

                NSTART = NLAYERS
                NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

                DO UTA = N_LOUTPUT, 1, -1

!  Layer index for given optical depth

                  NLEVEL = LEVELMASK_UP(UTA)
                  NUT    = NLEVEL + 1

!  Cumulative layer transmittance :
!    Loop over layers working upwards to level nut

                  DO N = NSTART, NUT, -1
                    NC = NLAYERS + 1 - N
                    LSSL_DB_CUMSOURCE(V) = UP_LOSTRANS(N) * LSSL_DB_CUMSOURCE(V)
                  ENDDO

!  Set output

                  LS_ELASTIC_SS_UP(Q1,UTA,V,1) = LSSL_DB_CUMSOURCE(V)
                  IF ( DO_FULL_RAMAN ) LS_RAMAN_SS_UP(Q1,UTA,V,1) = LSSL_DB_CUMSOURCE(V)

!  Check for updating the recursion

                  IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                  NUT_PREV = NUT

!  End optical depth loop

                enddo

!  End loop over all sleave weighting functions

              ENDDO

!  No surface clause

            ENDIF

!  End Upwelling clause (Sleave Jacobian)

         ENDIF 

!  Finish geometry loop

        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_SSCORR_GENERAL_MONO_PLUS_1

!

      SUBROUTINE OG_INTEGRATION_BIN_UP_PLUS_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXBINS, MAXVARS,  & ! Inputs
          DO_RAMAN, NPOINTS_INNER, OFFSET_INNER,                  & ! Inputs
          NBINS, BINMAP, NLAYERS, NFINELAYERS,                    & ! Inputs
          DO_PROFILE_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,     & ! Inputs
          DO_COLUMN_WFS, N_COLUMNWFS,                             & ! Inputs
          EXTINCTION, L_EXTINCTION, RADII, ALPHA_ALL, ALPHA_FINE, & ! Inputs
          ATTN, ATTN_FINE, L_ATTN, L_ATTN_FINE,                   & ! Inputs
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS,                   & ! Inputs
          L_EMULTIPLIERS, L_IMULTIPLIERS, L_LOSTRANS )              ! Outputs

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007 (Removed for version 2.3 on)

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, HALF

      IMPLICIT NONE

!  Inputs
!  ------

!  Dimensioning

      INTEGER  , INTENT(IN) :: maxlayers, maxfinelayers
      INTEGER  , INTENT(IN) :: maxpoints, maxbins, maxvars

!  Control

      LOGICAL  , INTENT(IN) :: do_raman
      INTEGER  , INTENT(IN) :: nfinelayers, nlayers
      INTEGER  , INTENT(IN) :: npoints_inner, offset_inner
      INTEGER  , INTENT(IN) :: nbins  ( maxpoints )
      INTEGER  , INTENT(IN) :: binmap ( maxbins, maxpoints )

!  Profile WFS control inputs

      LOGICAL  , INTENT(IN) :: do_profile_WFS
      LOGICAL  , INTENT(IN) :: layer_vary_flag   (maxlayers)
      INTEGER  , INTENT(IN) :: layer_vary_number (maxlayers)

!  Column WFS control inputs

      LOGICAL  , INTENT(IN) :: do_column_WFS
      INTEGER  , INTENT(IN) :: N_COLUMNWFS

!  Radii

      REAL(FPK), INTENT(IN) :: radii   (0:maxlayers)

!  Line of sight angles

      REAL(FPK), INTENT(IN) :: alpha_all  (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_fine    (maxlayers,maxfinelayers)

!  Extinction inputs

      REAL(FPK), INTENT(IN) :: extinction   (maxlayers,maxpoints)
      REAL(FPK), INTENT(IN) :: L_extinction (maxvars,maxlayers,maxpoints)

!  Attenuations

      REAL(FPK), INTENT(IN) :: attn      ( 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: attn_fine &
               ( maxlayers, maxfinelayers, maxpoints )

!  Linearized attenuations

      REAL(FPK), INTENT(IN) :: l_attn &
         ( maxvars, 0:maxlayers, 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: l_attn_fine &
         ( maxvars, 0:maxlayers, maxlayers, maxfinelayers, maxpoints )

!  Outputs
!  -------

!  Regular quantities

      REAL(FPK), INTENT(OUT) :: imultipliers   (maxlayers,   maxbins, maxpoints)
      REAL(FPK), INTENT(OUT) :: emultipliers   (maxlayers,   maxpoints)
      REAL(FPK), INTENT(OUT) :: lostrans      (maxlayers,   maxpoints)

!  Linearized quantities

      REAL(FPK), INTENT(OUT) :: l_imultipliers &
            ( maxvars, 0:maxlayers, maxlayers, maxbins, maxpoints)
      REAL(FPK), INTENT(OUT) :: l_emultipliers &
            ( maxvars, 0:maxlayers, maxlayers,   maxpoints)
      REAL(FPK), INTENT(OUT) :: l_lostrans (maxvars, maxlayers, maxpoints)

!  Local arrays
!  ------------

!  Local geometry arrays

      REAL(FPK) :: csq_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: cot_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: tranj     ( maxfinelayers)

!  Help variables
!  --------------

      INTEGER :: n, j, q, k
      INTEGER :: w, b, wb, wr

      REAL(FPK) :: argm_1, argj
      REAL(FPK) :: l_emult, l_imult
      REAL(FPK) :: l_func_1, l_func_2, l_tran_1, l_tran
      REAL(FPK) :: l_sum, l_func, l_kn

      REAL(FPK) :: esum,   rsum  (maxpoints)
      REAL(FPK) :: salpha, calpha, dfine1, raycon, kn
      REAL(FPK) :: csq_1, csqt_1, csq_2, cot_1, cot_2
      REAL(FPK) :: func_1, func_2, tran_1, step, func, skn

!  Local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0D0

!  Initialise output
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

!  Whole layer WFSs

      if ( do_profile_WFS ) then
        do n = 1, nlayers
         do w = 1, npoints_inner
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_emultipliers(q,k,n,w) = zero
                do b = 1, nbins(w)
                  l_imultipliers(q,k,n,b,w) = zero
                enddo
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans(q,n,w) = zero
            enddo
          endif
         enddo
        enddo
      endif

      if ( do_column_WFS ) then
        do n = 1, nlayers
         do w = 1, npoints_inner
          do q = 1, N_COLUMNWFS
            L_emultipliers(q,0,n,w) = zero
            L_lostrans(q,n,w)        = zero
            do b = 1, nbins(w)
              l_imultipliers(q,0,n,b,w) = zero
            enddo
          enddo
         enddo
        enddo
      endif

!  Ray constant at TOA

      salpha = sin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work up from the bottom of the atmosphere
!  =========================================

!  Initialise

      n = nlayers
      salpha = sin(alpha_all(n))
      calpha = cos(alpha_all(n))
      csq_1 = one / salpha / salpha
      cot_1 = calpha / salpha

!  Save some quantities

      do n = nlayers, 1, -1
        do j = 1, nfinelayers
          calpha = cos(alpha_fine(n,j))
          salpha = sin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = one / salpha / salpha
        enddo
      enddo

!  Start layer loop
!  ================

      do n = nlayers, 1, -1

!  Save some quantities

        salpha = sin(alpha_all(n-1))
        calpha = cos(alpha_all(n-1))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  Start spectral loop

        do w = 1, npoints_inner

!  Set up

          wr = w + offset_inner
          kn = raycon * extinction(n,wr)
          skn = step * kn
          argm_1 = cot_2 - cot_1
          tran_1 = dexp ( - kn * argm_1 )
          csqt_1 = csq_1 * tran_1
          do j = 1, nfinelayers
            tranj(j) = dexp ( - kn * ( cot_2 - cot_fine(n,j) ) )
          enddo

!  Intensity multipliers + transmittance
!  -------------------------------------

!  Line of sight transmittance factor

          lostrans(n,w)    = tran_1

!  Elastic scattering multiplier

          func_2 = attn(n-1,wr) *  csq_2
          func_1 = attn(n,wr)   * csqt_1
          esum = half * ( func_1 + func_2 )
          do j = 1, nfinelayers
            func = attn_fine(n,j,wr) * tranj(j) * csq_fine(n,j)
            esum  = esum + func
          enddo
          emultipliers(n,w) = esum * skn

!  Inelastic scattering multipliers (by bins)

          do b = 1, nbins(w)
            wb = binmap(b,w)
            func_2 = attn(n-1,wb) *  csq_2
            func_1 = attn(n,wb)   * csqt_1
            rsum(b) = half * ( func_1 + func_2 )
            do j = 1, nfinelayers
              func = attn_fine(n,j,wb) * tranj(j) * csq_fine(n,j)
              rsum(b)  = rsum(b) + func
            enddo
            imultipliers(n,b,w) = rsum(b) * skn
          enddo

!  Profile WFSs
!  ----------------------

          if ( do_profile_WFS ) then

!  Elastic

           do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
             do q = 1, layer_vary_number(k)
              if ( k.eq.n ) then
               l_kn = raycon * l_extinction(q,n,wr)
               l_func_2 = l_attn(q,n,n-1,wr) * csq_2
               l_tran_1 = -l_kn * argm_1 * tran_1
               l_func_1 = csq_1 * ( l_attn(q,n,n,wr) *   tran_1 &
                                    + attn(n,wr)     * l_tran_1 )
               l_sum = half * ( l_func_1 + l_func_2 )
               do j = 1, nfinelayers
                argj = cot_2 - cot_fine(n,j)
                l_tran = - l_kn * argj * tranj(j)
                l_func =  ( l_attn_fine(q,n,n,j,wr) * tranj(j) + &
                              attn_fine(n,j,wr)       * l_tran )
                l_sum = l_sum + csq_fine(n,j) * l_func
               enddo
               l_emult =  l_kn * esum + kn * l_sum
               l_lostrans(q,n,w)       = l_tran_1
               l_emultipliers(q,n,n,w) = step * l_emult
              else
               l_func_2 = l_attn(q,k,n-1,wr) * csq_2
               l_func_1 = csq_1 * l_attn(q,k,n,wr) * tran_1
               l_sum = half * ( l_func_1 + l_func_2 )
               do j = 1, nfinelayers
                l_func = l_attn_fine(q,k,n,j,wr) * tranj(j)
                l_sum  = l_sum + csq_fine(n,j) * l_func
               enddo
               l_emultipliers(q,k,n,w) = step * kn * l_sum
              endif
             enddo
            endif
           enddo

!  Inelastic

           if ( do_raman ) then
            do k = 1, nlayers
             if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
               do b = 1, nbins(w)
                wb = binmap(b,w)
                if ( k.eq.n) then
                 l_kn = raycon * l_extinction(q,n,wr)
                 l_func_2 = l_attn(q,n,n-1,wb) * csq_2
                 l_tran_1 = -l_kn * argm_1 * tran_1
                 l_func_1 = csq_1 * ( l_attn(q,n,n,wb) *   tran_1 &
                                      + attn(n,wb)     * l_tran_1 )
                 l_sum = half * ( l_func_1 + l_func_2 )
                 do j = 1, nfinelayers
                  argj = cot_2 - cot_fine(n,j)
                  l_tran = - l_kn * argj * tranj(j)
                  l_func =  ( l_attn_fine(q,n,n,j,wb) * tranj(j) &
                            +   attn_fine(n,j,wb)       * l_tran )
                  l_sum  = l_sum + csq_fine(n,j) * l_func
                 enddo
                 l_imult = l_kn * rsum(b) + kn * l_sum
                 l_imultipliers(q,n,n,b,w) = step * l_imult
                else
                 l_func_2 = l_attn(q,k,n-1,wb) * csq_2
                 l_func_1 = csq_1 * l_attn(q,k,n,wb) * tran_1
                 l_sum = half * ( l_func_1 + l_func_2 )
                 do j = 1, nfinelayers
                  l_func = l_attn_fine(q,k,n,j,wb) * tranj(j)
                  l_sum  = l_sum + csq_fine(n,j) * l_func
                 enddo
                 l_imultipliers(q,k,n,b,w) = step * kn * l_sum
                endif
               enddo
              enddo
             endif
            enddo
           endif

!  End profile WFS clause

          endif

!  Column WFSs
!  ---------------------

          if ( do_column_WFS ) then


!  Elastic

           do q = 1, N_COLUMNWFS
            l_kn = raycon * l_extinction(q,n,wr)
            l_func_2 = l_attn(q,0,n-1,wr) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(q,0,n,wr) *   tran_1 &
                                 + attn(n,wr)     * l_tran_1 )
            l_sum = half * ( l_func_1 + l_func_2 )
            do j = 1, nfinelayers
             argj = cot_2 - cot_fine(n,j)
             l_tran = - l_kn * argj * tranj(j)
             l_func = ( l_attn_fine(q,0,n,j,wr) * tranj(j) &
                      +   attn_fine(n,j,wr)      * l_tran )
             l_sum = l_sum + csq_fine(n,j) * l_func
            enddo
            l_lostrans(q,n,w)        = l_tran_1
            l_emult                  = l_kn * esum + kn * l_sum
            l_emultipliers(q,0,n,w) = step * l_emult
           enddo

!  Inelastic

           if ( do_raman ) then
            do q = 1, N_COLUMNWFS
             l_kn = raycon * l_extinction(q,n,wr)
             do b = 1, nbins(w)
              wb = binmap(b,w)
              l_func_2 = l_attn(q,0,n-1,wb) * csq_2
              l_tran_1 = -l_kn * argm_1 * tran_1
              l_func_1 = csq_1 * ( l_attn(q,0,n,wb) *   tran_1 &
                                   + attn(n,wb)      * l_tran_1 )
              l_sum = half * ( l_func_1 + l_func_2 )
              do j = 1, nfinelayers
               argj = cot_2 - cot_fine(n,j)
               l_tran = - l_kn * argj * tranj(j)
               l_func = ( l_attn_fine(q,0,n,j,wb) * tranj(j) &
                        +   attn_fine(n,j,wb)      * l_tran )
               l_sum = l_sum + csq_fine(n,j) * l_func
              enddo
              l_imult                  = l_kn * rsum(b) + kn * l_sum
              l_imultipliers(q,0,n,b,w) = step * l_imult
             enddo
            enddo
           endif

!  End column WFS clause

          endif

!  End spectral points loop

        enddo

!  Update geometry

        csq_1 = csq_2
        cot_1 = cot_2

!  Finish layer

      ENDDO

!  finish

      RETURN
      END SUBROUTINE OG_INTEGRATION_BIN_UP_PLUS_1

!

      SUBROUTINE OG_INTEGRATION_BIN_DN_PLUS_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXBINS, MAXVARS,  & ! Inputs
          DO_RAMAN, NPOINTS_INNER, OFFSET_INNER,                  & ! Inputs
          NBINS, BINMAP, NLAYERS, NFINELAYERS,                    & ! Inputs
          DO_PROFILE_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,     & ! Inputs
          DO_COLUMN_WFS,N_COLUMNWFS,                             & ! Inputs
          EXTINCTION, L_EXTINCTION, RADII, ALPHA_ALL, ALPHA_FINE, & ! Inputs
          ATTN, ATTN_FINE, L_ATTN, L_ATTN_FINE,                   & ! Inputs
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS,                   & ! Inputs
          L_EMULTIPLIERS, L_IMULTIPLIERS, L_LOSTRANS )              ! Outputs

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007 (Removed for version 2.3 on)

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, HALF

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER  , INTENT(IN) :: maxlayers, maxfinelayers
      INTEGER  , INTENT(IN) :: maxpoints, maxbins, maxvars

!  control

      LOGICAL  , INTENT(IN) :: do_raman
      INTEGER  , INTENT(IN) :: nfinelayers, nlayers
      INTEGER  , INTENT(IN) :: npoints_inner, offset_inner
      INTEGER  , INTENT(IN) :: nbins  ( maxpoints )
      INTEGER  , INTENT(IN) :: binmap ( maxbins, maxpoints )

!  Profile WFS control inputs

      LOGICAL  , INTENT(IN) :: do_profile_WFS
      LOGICAL  , INTENT(IN) :: layer_vary_flag   (maxlayers)
      INTEGER  , INTENT(IN) :: layer_vary_number (maxlayers)

!  Column WFS control inputs

      LOGICAL  , INTENT(IN) :: do_column_WFS
      INTEGER  , INTENT(IN) :: N_COLUMNWFS

!  radii

      REAL(FPK), INTENT(IN) :: radii   (0:maxlayers)

!  line of sight angles

      REAL(FPK), INTENT(IN) :: alpha_all  (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_fine    (maxlayers,maxfinelayers)

!  Extinction inputs

      REAL(FPK), INTENT(IN) :: extinction   (maxlayers,maxpoints)
      REAL(FPK), INTENT(IN) :: L_extinction (maxvars,maxlayers,maxpoints)

!  Attenuations

      REAL(FPK), INTENT(IN) :: attn      ( 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: attn_fine &
               ( maxlayers, maxfinelayers, maxpoints )

!  Linearized attenuations

      REAL(FPK), INTENT(IN) :: l_attn &
         ( maxvars, 0:maxlayers, 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: l_attn_fine &
         ( maxvars, 0:maxlayers, maxlayers, maxfinelayers, maxpoints )

!  outputs
!  -------

!  Regular quantities

      REAL(FPK), INTENT(OUT) :: imultipliers   (maxlayers, maxbins, maxpoints)
      REAL(FPK), INTENT(OUT) :: emultipliers   (maxlayers, maxpoints)
      REAL(FPK), INTENT(OUT) :: lostrans       (maxlayers, maxpoints)

!  Linearized quantities

      REAL(FPK), INTENT(OUT) :: l_imultipliers &
               ( maxvars, 0:maxlayers, maxlayers, maxbins, maxpoints)
      REAL(FPK), INTENT(OUT) :: l_emultipliers &
               ( maxvars, 0:maxlayers, maxlayers, maxpoints)
      REAL(FPK), INTENT(OUT) :: l_lostrans (maxvars, maxlayers, maxpoints)

!  local arrays
!  ------------

!  Local geometry arrays

      REAL(FPK) :: csq_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: cot_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: tranj     ( maxfinelayers)

!  help variables
!  --------------

      INTEGER :: n, j, q, k
      INTEGER :: w, b, wb, wr

      REAL(FPK) :: argm_1, argj
      REAL(FPK) :: l_emult, l_imult
      REAL(FPK) :: l_func_1, l_func_2, l_tran_1, l_tran
      REAL(FPK) :: l_sum, l_func, l_kn

      REAL(FPK) :: esum,   rsum  (maxpoints)
      REAL(FPK) :: salpha, calpha, dfine1, raycon, kn
      REAL(FPK) :: csq_1, csqt_1, csq_2, cot_1, cot_2
      REAL(FPK) :: func_1, func_2, tran_1, step, func, skn

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0D0

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

!  Whole layer WFSs

      if ( do_profile_WFS ) then
        do n = 1, nlayers
         do w = 1, npoints_inner
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_emultipliers(q,k,n,w) = zero
                do b = 1, nbins(w)
                  l_imultipliers(q,k,n,b,w) = zero
                enddo
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans(q,n,w) = zero
            enddo
          endif
         enddo
        enddo
      endif

      if ( do_column_WFS ) then
        do n = 1, nlayers
         do w = 1, npoints_inner
          do q = 1, N_COLUMNWFS
            L_emultipliers(q,0,n,w) = zero
            L_lostrans(q,n,w)        = zero
            do b = 1, nbins(w)
              l_imultipliers(q,0,n,b,w) = zero
            enddo
          enddo
         enddo
        enddo
      endif

!  Ray constant at TOA

      salpha = sin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work down from the top of the atmosphere
!  ========================================

!  initialise

      n = 0
      salpha = sin(alpha_all(n))
      calpha = cos(alpha_all(n))
      csq_1 = one / salpha / salpha
      cot_1 = calpha / salpha

!  save some quantities

      do n = 1, nlayers
        do j = 1, nfinelayers
          calpha = cos(alpha_fine(n,j))
          salpha = sin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = one / salpha / salpha
        enddo
      enddo

!  Start layer loop
!  ================

      do n = 1, nlayers

!  Save some quantities

        salpha = sin(alpha_all(n))
        calpha = cos(alpha_all(n))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  Start spectral loop

        do w = 1, npoints_inner

!  set up

          wr = w + offset_inner
          kn = raycon * extinction(n,wr)
          skn = step * kn
          argm_1 =  cot_1 - cot_2
          tran_1 = dexp ( - kn * argm_1  )
          csqt_1 = csq_1 * tran_1
          do j = nfinelayers, 1, -1
            tranj(j) = dexp ( - kn * ( cot_fine(n,j) -cot_2 ) )
          enddo

!  Intensity multipliers + transmittance
!  -------------------------------------

!  Line of sight transmittance factor

          lostrans(n,w)    = tran_1

!  Elastic scattering multiplier

          func_2 = attn(n,wr)   * csq_2
          func_1 = attn(n-1,wr) * csqt_1
          esum = half * ( func_1 + func_2 )
          do j = 1, nfinelayers
            func = attn_fine(n,j,wr) * tranj(j) * csq_fine(n,j)
            esum  = esum + func
          enddo
          emultipliers(n,w) = esum * skn

!  Inelastic scattering multipliers (by bins)

          do b = 1, nbins(w)
            wb = binmap(b,w)
            func_2 = attn(n,wb)   *  csq_2
            func_1 = attn(n-1,wb) * csqt_1
            rsum(b) = half * ( func_1 + func_2 )
            do j = 1, nfinelayers
              func = attn_fine(n,j,wb) * tranj(j) * csq_fine(n,j)
              rsum(b) = rsum(b) + func
            enddo
            imultipliers(n,b,w) = rsum(b) * skn
          enddo

!  Profile WFSs
!  ----------------------

          if ( do_profile_WFS ) then

!  Elastic

           do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
             do q = 1, layer_vary_number(k)
              if ( k.eq.n) then
               l_kn = raycon * l_extinction(q,n,wr)
               l_func_2 = l_attn(q,n,n,wr) * csq_2
               l_tran_1 = -l_kn * argm_1 * tran_1
               l_func_1 = csq_1 * ( l_attn(q,n,n-1,wr) *   tran_1 &
                                    + attn(n-1,wr)     * l_tran_1 )
               l_sum = half * ( l_func_1 + l_func_2 )
               do j = nfinelayers, 1, -1
                argj = cot_fine(n,j) - cot_2
                l_tran = - l_kn * argj * tranj(j)
                l_func =  ( l_attn_fine(q,n,n,j,wr) * tranj(j) + &
                              attn_fine(n,j,wr)       * l_tran )
                l_sum = l_sum + csq_fine(n,j) * l_func
               enddo
               l_emult =  l_kn * esum + kn * l_sum
               l_lostrans(q,n,w)       = l_tran_1
               l_emultipliers(q,n,n,w) = step * l_emult
              else
               l_func_1 = l_attn(q,k,n-1,wr) * csqt_1
               l_func_2 = l_attn(q,k,n,wr)   * csq_2
               l_sum = half * ( l_func_1 + l_func_2 )
               do j = 1, nfinelayers
                l_func = l_attn_fine(q,k,n,j,wr) * tranj(j)
                l_sum  = l_sum + csq_fine(n,j) * l_func
               enddo
               l_emultipliers(q,k,n,w) = step * kn * l_sum
              endif
             enddo
            endif
           enddo

!  Inelastic

           if ( do_raman ) then
            do k = 1, nlayers
             if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
               do b = 1, nbins(w)
                wb = binmap(b,w)
                if ( k.eq.n) then
                 l_kn = raycon * l_extinction(q,n,wr)
                 l_func_2 = l_attn(q,n,n,wb) * csq_2
                 l_tran_1 = -l_kn * argm_1 * tran_1
                 l_func_1 = csq_1 * ( l_attn(q,n,n-1,wb) *   tran_1 &
                                      + attn(n-1,wb)     * l_tran_1 )
                 l_sum = half * ( l_func_1 + l_func_2 )
                 do j = nfinelayers, 1, -1
                  argj = cot_fine(n,j) - cot_2
                  l_tran = - l_kn * argj * tranj(j)
                  l_func =  ( l_attn_fine(q,n,n,j,wb) * tranj(j) &
                            +   attn_fine(n,j,wb)       * l_tran )
                  l_sum  = l_sum + csq_fine(n,j) * l_func
                 enddo
                 l_imult = l_kn * rsum(b) + kn * l_sum
                 l_imultipliers(q,n,n,b,w) = step * l_imult
                else
                 l_func_2 = l_attn(q,k,n,wb)   * csq_2
                 l_func_1 = l_attn(q,k,n-1,wb) * csqt_1
                 l_sum = half * ( l_func_1 + l_func_2 )
                 do j = 1, nfinelayers
                  l_func = l_attn_fine(q,k,n,j,wb) * tranj(j)
                  l_sum  = l_sum + csq_fine(n,j) * l_func
                 enddo
                 l_imultipliers(q,k,n,b,w) = step * kn * l_sum
                endif
               enddo
              enddo
             endif
            enddo
           endif

!  End profile WFS clause

          endif

!  Column WFSs
!  ---------------------

          if ( do_column_WFS ) then

!  Elastic

           do q = 1, N_COLUMNWFS
            l_kn = raycon * l_extinction(q,n,wr)
            l_func_2 = l_attn(q,0,n,wr) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(q,0,n-1,wr) *   tran_1 &
                                 + attn(n-1,wr)      * l_tran_1 )
            l_sum = half * ( l_func_1 + l_func_2 )
            do j = nfinelayers, 1, -1
             argj = cot_fine(n,j) - cot_2
             l_tran = - l_kn * argj * tranj(j)
             l_func = ( l_attn_fine(q,0,n,j,wr) * tranj(j) &
                      +   attn_fine(n,j,wr)      * l_tran )
             l_sum = l_sum + csq_fine(n,j) * l_func
            enddo
            l_lostrans(q,n,w)        = l_tran_1
            l_emult                  = l_kn * esum + kn * l_sum
            l_emultipliers(q,0,n,w)  = step * l_emult
           enddo

!  Inelastic

           if ( do_raman ) then
            do q = 1, N_COLUMNWFS
             l_kn = raycon * l_extinction(q,n,wr)
             do b = 1, nbins(w)
              wb = binmap(b,w)
              l_func_2 = l_attn(q,0,n,wb) * csq_2
              l_tran_1 = -l_kn * argm_1 * tran_1
              l_func_1 = csq_1 * ( l_attn(q,0,n-1,wb) *   tran_1 &
                                   + attn(n-1,wb)      * l_tran_1 )
              l_sum = half * ( l_func_1 + l_func_2 )
              do j = nfinelayers, 1, -1
               argj = cot_fine(n,j) - cot_2
               l_tran = - l_kn * argj * tranj(j)
               l_func = ( l_attn_fine(q,0,n,j,wb) * tranj(j) &
                        +   attn_fine(n,j,wb)      * l_tran )
               l_sum = l_sum + csq_fine(n,j) * l_func
              enddo
              l_imult                  = l_kn * rsum(b) + kn * l_sum
              l_imultipliers(q,0,n,b,w) = step * l_imult
             enddo
            enddo
           endif

!  end column WFS clause

          endif

!  End spectral points loop

        enddo

!  update geometry

        csq_1 = csq_2
        cot_1 = cot_2

!  Finish layer

      ENDDO

!  finish

      RETURN
      END SUBROUTINE OG_INTEGRATION_BIN_DN_PLUS_1

!

      SUBROUTINE OG_INTEGRATION_MONO_UP_PLUS_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXVARS,            & ! Inputs
          DO_RAMAN, NPOINTS_MONO, W_EXCIT, NLAYERS, NFINELAYERS,   & ! Inputs
          DO_PROFILE_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,      & ! Inputs
          DO_COLUMN_WFS, N_COLUMNWFS,                             & ! Inputs
          EXTINCTION, L_EXTINCTION, RADII, ALPHA_ALL, ALPHA_FINE,  & ! Inputs
          ATTN, ATTN_FINE, L_ATTN, L_ATTN_FINE,                    & ! Inputs
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS,                    & ! Inputs
          L_EMULTIPLIERS, L_IMULTIPLIERS, L_LOSTRANS )               ! Outputs

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007 (Removed for version 2.3 on)

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, HALF

       IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER  , INTENT(IN) :: maxlayers, maxfinelayers
      INTEGER  , INTENT(IN) :: maxpoints, maxvars

!  control

      LOGICAL  , INTENT(IN) :: do_raman
      INTEGER  , INTENT(IN) :: nfinelayers, nlayers
      INTEGER  , INTENT(IN) :: npoints_mono, w_excit

!  Profile WFS control inputs

      LOGICAL  , INTENT(IN) :: do_profile_WFS
      LOGICAL  , INTENT(IN) :: layer_vary_flag   (maxlayers)
      INTEGER  , INTENT(IN) :: layer_vary_number (maxlayers)

!  Column WFS control inputs

      LOGICAL  , INTENT(IN) :: do_column_WFS
      INTEGER  , INTENT(IN) :: N_COLUMNWFS

!  radii

      REAL(FPK), INTENT(IN) :: radii   (0:maxlayers)

!  line of sight angles

      REAL(FPK), INTENT(IN) :: alpha_all  (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_fine    (maxlayers,maxfinelayers)

!  Extinction inputs

      REAL(FPK), INTENT(IN) :: extinction   (maxlayers,maxpoints)
      REAL(FPK), INTENT(IN) :: L_extinction (maxvars,maxlayers,maxpoints)

!  Attenuations

      REAL(FPK), INTENT(IN) :: attn      ( 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: attn_fine &
               ( maxlayers, maxfinelayers, maxpoints )

!  Linearized attenuations

      REAL(FPK), INTENT(IN) :: l_attn &
         ( maxvars, 0:maxlayers, 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: l_attn_fine &
         ( maxvars, 0:maxlayers, maxlayers, maxfinelayers, maxpoints )

!  outputs
!  -------

!  Regular quantities

      REAL(FPK), INTENT(OUT) :: imultipliers   (maxlayers, maxpoints)
      REAL(FPK), INTENT(OUT) :: emultipliers   (maxlayers)
      REAL(FPK), INTENT(OUT) :: lostrans       (maxlayers)

!  Linearized quantities

      REAL(FPK), INTENT(OUT) :: l_imultipliers &
            ( maxvars, 0:maxlayers, maxlayers, maxpoints )
      REAL(FPK), INTENT(OUT) :: l_emultipliers &
            ( maxvars, 0:maxlayers, maxlayers )
      REAL(FPK), INTENT(OUT) :: l_lostrans (maxvars, maxlayers )

!  local arrays
!  ------------

!  Local geometry arrays

      REAL(FPK) :: csq_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: cot_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: tranj     ( maxfinelayers)

!  help variables
!  --------------

      INTEGER   :: n, j, w, q, k

      REAL(FPK) :: argm_1, argj
      REAL(FPK) :: l_emult, l_imult
      REAL(FPK) :: l_func_1, l_func_2, l_tran_1, l_tran
      REAL(FPK) :: l_sum, l_func, l_kn

      REAL(FPK) :: esum,   rsum  (maxpoints)
      REAL(FPK) :: salpha, calpha, dfine1, raycon, kn
      REAL(FPK) :: csq_1, csqt_1, csq_2, cot_1, cot_2
      REAL(FPK) :: func_1, func_2, tran_1, step, func, skn

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0D0

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        lostrans(n)    = zero
        emultipliers(n) = zero
        do w = 1, npoints_mono
          imultipliers(n,w) = zero
        enddo
      enddo

!  Whole layer WFSs

      if ( do_profile_WFS ) then
        do n = 1, nlayers
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_emultipliers(q,k,n) = zero
                do w = 1, npoints_mono
                  l_imultipliers(q,k,n,w) = zero
                enddo
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans(q,n) = zero
            enddo
          endif
        enddo
      endif

      if ( do_column_WFS ) then
        do n = 1, nlayers
          do q = 1, N_COLUMNWFS
            L_emultipliers(q,0,n) = zero
            L_lostrans(q,n)        = zero
            do w = 1, npoints_mono
              l_imultipliers(q,0,n,w) = zero
            enddo
          enddo
        enddo
      endif

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
!  ================

      do n = nlayers, 1, -1

!  Save some quantities

        salpha = sin(alpha_all(n-1))
        calpha = cos(alpha_all(n-1))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  set up

        kn = raycon * extinction(n,w_excit)
        skn = step * kn
        argm_1 = cot_2 - cot_1
        tran_1 = dexp ( - kn * argm_1 )
        csqt_1 = csq_1 * tran_1
        do j = 1, nfinelayers
          tranj(j) = dexp ( - kn * ( cot_2 - cot_fine(n,j) ) )
        enddo

!  Intensity multipliers + transmittance
!  -------------------------------------

!  Line of sight transmittance factor

        lostrans(n)    = tran_1

!  Elastic scattering multiplier

        func_2 = attn(n-1,w_excit) *  csq_2
        func_1 = attn(n,w_excit)   * csqt_1
        esum = half * ( func_1 + func_2 )
        do j = 1, nfinelayers
          func = attn_fine(n,j,w_excit) * tranj(j) * csq_fine(n,j)
          esum  = esum + func
        enddo
        emultipliers(n) = esum * skn

!  Inelastic scattering multipliers (by bins)

        do w = 1, npoints_mono
          func_2 = attn(n-1,w) *  csq_2
          func_1 = attn(n,w)   * csqt_1
          rsum(w) = half * ( func_1 + func_2 )
          do j = 1, nfinelayers
            func = attn_fine(n,j,w) * tranj(j) * csq_fine(n,j)
            rsum(w)  = rsum(w) + func
          enddo
          imultipliers(n,w) = rsum(w) * skn
        enddo

!  Profile WFSs
!  ----------------------

        if ( do_profile_WFS ) then

!  Elastic

          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
             do q = 1, layer_vary_number(k)
              if ( k.eq.n ) then
               l_kn = raycon * l_extinction(q,n,w_excit)
               l_func_2 = l_attn(q,n,n-1,w_excit) * csq_2
               l_tran_1 = -l_kn * argm_1 * tran_1
               l_func_1 = csq_1 * ( l_attn(q,n,n,w_excit) *   tran_1 &
                                    + attn(n,w_excit)     * l_tran_1 )
               l_sum = half * ( l_func_1 + l_func_2 )
               do j = 1, nfinelayers
                argj = cot_2 - cot_fine(n,j)
                l_tran = - l_kn * argj * tranj(j)
                l_func =  ( l_attn_fine(q,n,n,j,w_excit) * tranj(j) + &
                              attn_fine(n,j,w_excit)       * l_tran )
                l_sum = l_sum + csq_fine(n,j) * l_func
               enddo
               l_emult =  l_kn * esum + kn * l_sum
               l_lostrans(q,n)       = l_tran_1
               l_emultipliers(q,n,n) = step * l_emult
              else
               l_func_2 = l_attn(q,k,n-1,w_excit) * csq_2
               l_func_1 = csq_1 * l_attn(q,k,n,w_excit) * tran_1
               l_sum = half * ( l_func_1 + l_func_2 )
               do j = 1, nfinelayers
                l_func = l_attn_fine(q,k,n,j,w_excit) * tranj(j)
                l_sum  = l_sum + csq_fine(n,j) * l_func
               enddo
               l_emultipliers(q,k,n) = step * kn * l_sum
              endif
             enddo
            endif
          enddo

!  Inelastic

          if ( do_raman ) then
            do k = 1, nlayers
             if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
               do w = 1, npoints_mono
                if ( k.eq.n) then
                 l_kn = raycon * l_extinction(q,n,w_excit)
                 l_func_2 = l_attn(q,n,n-1,w) * csq_2
                 l_tran_1 = -l_kn * argm_1 * tran_1
                 l_func_1 = csq_1 * ( l_attn(q,n,n,w) *   tran_1 &
                                      + attn(n,w)     * l_tran_1 )
                 l_sum = half * ( l_func_1 + l_func_2 )
                 do j = 1, nfinelayers
                  argj = cot_2 - cot_fine(n,j)
                  l_tran = - l_kn * argj * tranj(j)
                  l_func =  ( l_attn_fine(q,n,n,j,w) * tranj(j) &
                            +   attn_fine(n,j,w)       * l_tran )
                  l_sum  = l_sum + csq_fine(n,j) * l_func
                 enddo
                 l_imult = l_kn * rsum(w) + kn * l_sum
                 l_imultipliers(q,n,n,w) = step * l_imult
                else
                 l_func_2 = l_attn(q,k,n-1,w) * csq_2
                 l_func_1 = csq_1 * l_attn(q,k,n,w) * tran_1
                 l_sum = half * ( l_func_1 + l_func_2 )
                 do j = 1, nfinelayers
                  l_func = l_attn_fine(q,k,n,j,w) * tranj(j)
                  l_sum  = l_sum + csq_fine(n,j) * l_func
                 enddo
                 l_imultipliers(q,k,n,w) = step * kn * l_sum
                endif
               enddo
              enddo
             endif
            enddo
          endif

!  End profile WFS clause

        endif

!  Column WFSs
!  ---------------------

        if ( do_column_WFS ) then

!  Elastic

          do q = 1, N_COLUMNWFS
            l_kn = raycon * l_extinction(q,n,w_excit)
            l_func_2 = l_attn(q,0,n-1,w_excit) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(q,0,n,w_excit) *   tran_1 &
                                 + attn(n,w_excit)     * l_tran_1 )
            l_sum = half * ( l_func_1 + l_func_2 )
            do j = 1, nfinelayers
             argj = cot_2 - cot_fine(n,j)
             l_tran = - l_kn * argj * tranj(j)
             l_func = ( l_attn_fine(q,0,n,j,w_excit) * tranj(j) &
                      +   attn_fine(n,j,w_excit)      * l_tran )
             l_sum = l_sum + csq_fine(n,j) * l_func
            enddo
            l_lostrans(q,n)        = l_tran_1
            l_emult                = l_kn * esum + kn * l_sum
            l_emultipliers(q,0,n) = step * l_emult
          enddo

!  Inelastic

          if ( do_raman ) then
            do q = 1, N_COLUMNWFS
             l_kn = raycon * l_extinction(q,n,w_excit)
             do w = 1, npoints_mono
              l_func_2 = l_attn(q,0,n-1,w) * csq_2
              l_tran_1 = -l_kn * argm_1 * tran_1
              l_func_1 = csq_1 * ( l_attn(q,0,n,w) *   tran_1 &
                                   + attn(n,w)      * l_tran_1 )
              l_sum = half * ( l_func_1 + l_func_2 )
              do j = 1, nfinelayers
               argj = cot_2 - cot_fine(n,j)
               l_tran = - l_kn * argj * tranj(j)
               l_func = ( l_attn_fine(q,0,n,j,w) * tranj(j) &
                        +   attn_fine(n,j,w)      * l_tran )
               l_sum = l_sum + csq_fine(n,j) * l_func
              enddo
              l_imult                  = l_kn * rsum(w) + kn * l_sum
              l_imultipliers(q,0,n,w) = step * l_imult
             enddo
            enddo
          endif

!  end column WFS clause

        endif

!  update geometry

        csq_1 = csq_2
        cot_1 = cot_2

!  Finish layer

      ENDDO

!  finish

      RETURN
      END SUBROUTINE OG_INTEGRATION_MONO_UP_PLUS_1

!

      SUBROUTINE OG_INTEGRATION_MONO_DN_PLUS_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXVARS,           & ! Inputs
          DO_RAMAN, NPOINTS_MONO, W_EXCIT, NLAYERS, NFINELAYERS,  & ! Inputs
          DO_PROFILE_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,     & ! Inputs
          DO_COLUMN_WFS,N_COLUMNWFS,                             & ! Inputs
          EXTINCTION, L_EXTINCTION, RADII, ALPHA_ALL, ALPHA_FINE, & ! Inputs
          ATTN, ATTN_FINE, L_ATTN, L_ATTN_FINE,                   & ! Inputs
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS,                   & ! Inputs
          L_EMULTIPLIERS, L_IMULTIPLIERS, L_LOSTRANS )              ! Output

!  Does the optical depth integration over layers, with linearized output
!  Partial layer integration added September 2007 (Removed for version 2.3 on)

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, HALF

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER  , INTENT(IN) :: maxlayers, maxfinelayers
      INTEGER  , INTENT(IN) :: maxpoints, maxvars

!  control

      LOGICAL  , INTENT(IN) :: do_raman
      INTEGER  , INTENT(IN) :: nfinelayers, nlayers
      INTEGER  , INTENT(IN) :: npoints_mono, w_excit

!  Profile WFS control inputs

      LOGICAL  , INTENT(IN) :: do_profile_WFS
      LOGICAL  , INTENT(IN) :: layer_vary_flag   (maxlayers)
      INTEGER  , INTENT(IN) :: layer_vary_number (maxlayers)

!  Column WFS control inputs

      LOGICAL  , INTENT(IN) :: do_column_WFS
      INTEGER  , INTENT(IN) :: N_COLUMNWFS

!  radii

      REAL(FPK), INTENT(IN) :: radii   (0:maxlayers)

!  line of sight angles

      REAL(FPK), INTENT(IN) :: alpha_all  (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_fine    (maxlayers,maxfinelayers)

!  Extinction inputs

      REAL(FPK), INTENT(IN) :: extinction   (maxlayers,maxpoints)
      REAL(FPK), INTENT(IN) :: L_extinction (maxvars,maxlayers,maxpoints)

!  Attenuations

      REAL(FPK), INTENT(IN) :: attn      ( 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: attn_fine &
               ( maxlayers, maxfinelayers, maxpoints )

!  Linearized attenuations

      REAL(FPK), INTENT(IN) :: l_attn &
         ( maxvars, 0:maxlayers, 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(IN) :: l_attn_fine &
         ( maxvars, 0:maxlayers, maxlayers, maxfinelayers, maxpoints )

!  outputs
!  -------

!  Regular quantities

      REAL(FPK), INTENT(OUT) :: imultipliers   (maxlayers, maxpoints)
      REAL(FPK), INTENT(OUT) :: emultipliers   (maxlayers )
      REAL(FPK), INTENT(OUT) :: lostrans       (maxlayers )

!  Linearized quantities

      REAL(FPK), INTENT(OUT) :: l_imultipliers &
               ( maxvars, 0:maxlayers, maxlayers, maxpoints)
      REAL(FPK), INTENT(OUT) :: l_emultipliers &
               ( maxvars, 0:maxlayers, maxlayers )
      REAL(FPK), INTENT(OUT) :: l_lostrans (maxvars, maxlayers )

!  local arrays
!  ------------

!  Local geometry arrays

      REAL(FPK) :: csq_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: cot_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: tranj     ( maxfinelayers)

!  help variables
!  --------------

      INTEGER   :: n, j, q, k, w

      REAL(FPK) :: argm_1, argj
      REAL(FPK) :: l_emult, l_imult
      REAL(FPK) :: l_func_1, l_func_2, l_tran_1, l_tran
      REAL(FPK) :: l_sum, l_func, l_kn

      REAL(FPK) :: esum,   rsum  (maxpoints)
      REAL(FPK) :: salpha, calpha, dfine1, raycon, kn
      REAL(FPK) :: csq_1, csqt_1, csq_2, cot_1, cot_2
      REAL(FPK) :: func_1, func_2, tran_1, step, func, skn

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0D0

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        lostrans(n)    = zero
        emultipliers(n) = zero
        do w = 1, npoints_mono
          imultipliers(n,w) = zero
        enddo
      enddo

!  Whole layer WFSs

      if ( do_profile_WFS ) then
        do n = 1, nlayers
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_emultipliers(q,k,n) = zero
                do w = 1, npoints_mono
                  l_imultipliers(q,k,n,w) = zero
                enddo
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans(q,n) = zero
            enddo
          endif
        enddo
      endif

      if ( do_column_WFS ) then
        do n = 1, nlayers
          do q = 1, N_COLUMNWFS
            L_emultipliers(q,0,n) = zero
            L_lostrans(q,n)        = zero
            do w = 1, npoints_mono
              l_imultipliers(q,0,n,w) = zero
            enddo
          enddo
        enddo
      endif

!  Ray constant at TOA

      salpha = sin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work down from the top of the atmosphere
!  ========================================

!  initialise

      n = 0
      salpha = sin(alpha_all(n))
      calpha = cos(alpha_all(n))
      csq_1 = one / salpha / salpha
      cot_1 = calpha / salpha

!  save some quantities

      do n = 1, nlayers
        do j = 1, nfinelayers
          calpha = cos(alpha_fine(n,j))
          salpha = sin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = one / salpha / salpha
        enddo
      enddo

!  Start layer loop
!  ================

      do n = 1, nlayers

!  Save some quantities

        salpha = sin(alpha_all(n))
        calpha = cos(alpha_all(n))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  set up

        kn = raycon * extinction(n,w_excit)
        skn = step * kn
        argm_1 =  cot_1 - cot_2
        tran_1 = dexp ( - kn * argm_1  )
        csqt_1 = csq_1 * tran_1
        do j = nfinelayers, 1, -1
          tranj(j) = dexp ( - kn * ( cot_fine(n,j) -cot_2 ) )
        enddo

!  Intensity multipliers + transmittance
!  -------------------------------------

!  Line of sight transmittance factor

        lostrans(n)    = tran_1

!  Elastic scattering multiplier

        func_2 = attn(n,w_excit)   * csq_2
        func_1 = attn(n-1,w_excit) * csqt_1
        esum = half * ( func_1 + func_2 )
        do j = 1, nfinelayers
          func = attn_fine(n,j,w_excit) * tranj(j) * csq_fine(n,j)
          esum  = esum + func
        enddo
        emultipliers(n) = esum * skn

!  Inelastic scattering multipliers (by bins)

        do w = 1, npoints_mono
          func_2 = attn(n,w)   *  csq_2
          func_1 = attn(n-1,w) * csqt_1
          rsum(w) = half * ( func_1 + func_2 )
          do j = 1, nfinelayers
            func = attn_fine(n,j,w) * tranj(j) * csq_fine(n,j)
            rsum(w) = rsum(w) + func
          enddo
          imultipliers(n,w) = rsum(w) * skn
        enddo

!  Profile WFSs
!  ----------------------

        if ( do_profile_WFS ) then

!  Elastic

          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
             do q = 1, layer_vary_number(k)
              if ( k.eq.n) then
               l_kn = raycon * l_extinction(q,n,w_excit)
               l_func_2 = l_attn(q,n,n,w_excit) * csq_2
               l_tran_1 = -l_kn * argm_1 * tran_1
               l_func_1 = csq_1 * ( l_attn(q,n,n-1,w_excit) *   tran_1 &
                                    + attn(n-1,w_excit)     * l_tran_1 )
               l_sum = half * ( l_func_1 + l_func_2 )
               do j = nfinelayers, 1, -1
                argj = cot_fine(n,j) - cot_2
                l_tran = - l_kn * argj * tranj(j)
                l_func =  ( l_attn_fine(q,n,n,j,w_excit) * tranj(j) + &
                              attn_fine(n,j,w_excit)       * l_tran )
                l_sum = l_sum + csq_fine(n,j) * l_func
               enddo
               l_emult =  l_kn * esum + kn * l_sum
               l_lostrans(q,n)       = l_tran_1
               l_emultipliers(q,n,n) = step * l_emult
              else
               l_func_1 = l_attn(q,k,n-1,w_excit) * csqt_1
               l_func_2 = l_attn(q,k,n,w_excit)   * csq_2
               l_sum = half * ( l_func_1 + l_func_2 )
               do j = 1, nfinelayers
                l_func = l_attn_fine(q,k,n,j,w_excit) * tranj(j)
                l_sum  = l_sum + csq_fine(n,j) * l_func
               enddo
               l_emultipliers(q,k,n) = step * kn * l_sum
              endif
             enddo
            endif
          enddo

!  Inelastic

          if ( do_raman ) then
            do k = 1, nlayers
             if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
               do w = 1, npoints_mono
                if ( k.eq.n) then
                 l_kn = raycon * l_extinction(q,n,w_excit)
                 l_func_2 = l_attn(q,n,n,w) * csq_2
                 l_tran_1 = -l_kn * argm_1 * tran_1
                 l_func_1 = csq_1 * ( l_attn(q,n,n-1,w) *   tran_1 &
                                      + attn(n-1,w)     * l_tran_1 )
                 l_sum = half * ( l_func_1 + l_func_2 )
                 do j = nfinelayers, 1, -1
                  argj = cot_fine(n,j) - cot_2
                  l_tran = - l_kn * argj * tranj(j)
                  l_func =  ( l_attn_fine(q,n,n,j,w) * tranj(j) &
                            +   attn_fine(n,j,w)       * l_tran )
                  l_sum  = l_sum + csq_fine(n,j) * l_func
                 enddo
                 l_imult = l_kn * rsum(w) + kn * l_sum
                 l_imultipliers(q,n,n,w) = step * l_imult
                else
                 l_func_2 = l_attn(q,k,n,w)   * csq_2
                 l_func_1 = l_attn(q,k,n-1,w) * csqt_1
                 l_sum = half * ( l_func_1 + l_func_2 )
                 do j = 1, nfinelayers
                  l_func = l_attn_fine(q,k,n,j,w) * tranj(j)
                  l_sum  = l_sum + csq_fine(n,j) * l_func
                 enddo
                 l_imultipliers(q,k,n,w) = step * kn * l_sum
                endif
               enddo
              enddo
             endif
            enddo
          endif

!  End profile WFS clause

        endif

!  Column WFSs
!  ---------------------

        if ( do_column_WFS ) then

!  Elastic

          do q = 1, N_COLUMNWFS
            l_kn = raycon * l_extinction(q,n,w_excit)
            l_func_2 = l_attn(q,0,n,w_excit) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(q,0,n-1,w_excit) *   tran_1 &
                                 + attn(n-1,w_excit)      * l_tran_1 )
            l_sum = half * ( l_func_1 + l_func_2 )
            do j = nfinelayers, 1, -1
             argj = cot_fine(n,j) - cot_2
             l_tran = - l_kn * argj * tranj(j)
             l_func = ( l_attn_fine(q,0,n,j,w_excit) * tranj(j) &
                      +   attn_fine(n,j,w_excit)      * l_tran )
             l_sum = l_sum + csq_fine(n,j) * l_func
            enddo
            l_lostrans(q,n)        = l_tran_1
            l_emult                = l_kn * esum + kn * l_sum
            l_emultipliers(q,0,n)  = step * l_emult
          enddo

!  Inelastic

          if ( do_raman ) then
            do q = 1, N_COLUMNWFS
             l_kn = raycon * l_extinction(q,n,w_excit)
             do w = 1, npoints_mono
              l_func_2 = l_attn(q,0,n,w) * csq_2
              l_tran_1 = -l_kn * argm_1 * tran_1
              l_func_1 = csq_1 * ( l_attn(q,0,n-1,w) *   tran_1 &
                                   + attn(n-1,w)      * l_tran_1 )
              l_sum = half * ( l_func_1 + l_func_2 )
              do j = nfinelayers, 1, -1
               argj = cot_fine(n,j) - cot_2
               l_tran = - l_kn * argj * tranj(j)
               l_func = ( l_attn_fine(q,0,n,j,w) * tranj(j) &
                        +   attn_fine(n,j,w)      * l_tran )
               l_sum = l_sum + csq_fine(n,j) * l_func
              enddo
              l_imult                  = l_kn * rsum(w) + kn * l_sum
              l_imultipliers(q,0,n,w) = step * l_imult
             enddo
            enddo
          endif

!  end column WFS clause

        endif

!  update geometry

        csq_1 = csq_2
        cot_1 = cot_2

!  Finish layer

      ENDDO

!  finish

      RETURN
      END SUBROUTINE OG_INTEGRATION_MONO_DN_PLUS_1

!

      SUBROUTINE OG_ATTENUATIONS_PLUS_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXVARS,       & ! Inputs
          N_RAMAN_POINTS, NLAYERS, NFINELAYERS,               & ! Inputs
          DO_PROFILE_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, & ! Inputs
          DO_COLUMN_WFS, N_COLUMNWFS,                        & ! Inputs
          EXTINCTION, L_EXTINCTION,                           & ! Inputs
          SUNPATHS, NTRAVERSE, SUNPATHS_FINE, NTRAVERSE_FINE, & ! Inputs
          ATTN, ATTN_FINE, L_ATTN, L_ATTN_FINE )                ! Outputs

!  Does attenuations with linearization

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, HALF

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER  , INTENT(IN) :: maxlayers, maxfinelayers
      INTEGER  , INTENT(IN) :: maxpoints, maxvars

!  control

      INTEGER  , INTENT(IN) :: nfinelayers, nlayers, n_raman_points

!  Profile WFS control inputs

      LOGICAL  , INTENT(IN) :: do_profile_WFS
      LOGICAL  , INTENT(IN) :: layer_vary_flag   (maxlayers)
      INTEGER  , INTENT(IN) :: layer_vary_number (maxlayers)

!  Column WFS control inputs

      LOGICAL  , INTENT(IN) :: do_column_WFS
      INTEGER  , INTENT(IN) :: N_COLUMNWFS

!  Whole layers

      INTEGER  , INTENT(IN) :: ntraverse(0:maxlayers)
      REAL(FPK), INTENT(IN) :: sunpaths(0:maxlayers,maxlayers)

!  Fine level

      INTEGER  , INTENT(IN) :: ntraverse_fine(maxlayers,maxfinelayers)
      REAL(FPK), INTENT(IN) :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)

!  Extinction inputs

      REAL(FPK), INTENT(IN) :: extinction   (maxlayers,maxpoints)
      REAL(FPK), INTENT(IN) :: L_extinction (maxvars,maxlayers,maxpoints)

!  outputs
!  -------

      REAL(FPK), INTENT(OUT) :: attn      ( 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(OUT) :: attn_fine &
               ( maxlayers, maxfinelayers, maxpoints )

      REAL(FPK), INTENT(OUT) :: l_attn &
         ( maxvars, 0:maxlayers, 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(OUT) :: l_attn_fine &
         ( maxvars, 0:maxlayers, maxlayers, maxfinelayers, maxpoints )

!  help variables
!  --------------

      INTEGER   :: w, n, j, k, q
      REAL(FPK) :: tau, l_iop, sumr

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0D0

!  attenuation functions, whole layers

      do w = 1, n_raman_points
       do n = 0, nlayers
        tau     = zero
        attn(n,w) = zero
        do k = 1, ntraverse(n)
         tau = tau + sunpaths(n,k) * extinction(k,w)
        enddo
        if ( tau .le. local_cutoff ) attn(n,w) = dexp(-tau)
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
         if ( tau .le. local_cutoff ) attn_fine(n,j,w) = dexp(-tau)
        enddo
       enddo
      enddo

!  Linearized attenuations
!  -----------------------

!  Linearized attenuation factors. Profile WFS.
!     Note the extra zeroing...............bug, 01 November 2007

      if ( do_profile_WFS ) then
       do w = 1, n_raman_points

        do n = 0, nlayers
!mick chg 10/29/2015 - changed loop extent to initialize everything
          !do k = n+1, nlayers
          do k = 1, nlayers
            do q = 1, layer_vary_number(k)
              l_attn(q,k,n,w) = zero
            enddo
          enddo

          do k = 1, ntraverse(n)
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                l_iop = - l_extinction(q,k,w)
                l_attn(q,k,n,w) = sunpaths(n,k) * attn(n,w) * l_iop
              enddo
            endif
          enddo
        enddo

        do n = 1, nlayers
         do j = 1, nfinelayers
!mick fix 10/19/2015 - added loops to initialize everything
          do k = 1, nlayers
            do q = 1, layer_vary_number(k)
              l_attn_fine(q,k,n,j,w) = zero
            enddo
          enddo

          do k = 1, ntraverse_fine(n,j)
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                l_iop = - l_extinction(q,k,w) * attn_fine(n,j,w)
                l_attn_fine(q,k,n,j,w) = sunpaths_fine(n,k,j) * l_iop
              enddo
            endif
          enddo
         enddo
        enddo

       enddo !Raman pt loop
      endif

!  Linearized attenuation factors. Column WFS

      if ( do_column_WFS ) then
       do w = 1, n_raman_points
        do n = 0, nlayers
          do q = 1, N_COLUMNWFS
            sumr = zero
            do k = 1, ntraverse(n)
              l_iop = - l_extinction(q,k,w)
              sumr = sumr + sunpaths(n,k) * l_iop
            enddo
            l_attn(q,0,n,w) = sumr * attn(n,w)
          enddo
        enddo
        do n = 1, nlayers
         do j = 1, nfinelayers
          do q = 1, N_COLUMNWFS
            sumr = zero
            do k = 1, ntraverse_fine(n,j)
              l_iop = - l_extinction(q,k,w)
              sumr = sumr + sunpaths_fine(n,k,j) * l_iop
            enddo
            l_attn_fine(q,0,n,j,w) = sumr * attn_fine(n,j,w)
          enddo
         enddo
        enddo
       enddo
      endif

!  Finish

      return
      END SUBROUTINE OG_ATTENUATIONS_PLUS_1

!  End module

      END MODULE lrrs_corrections_plus_1_m
