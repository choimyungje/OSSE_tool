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
! # Subroutines in this Module                                  #
! #                                                             #
! #             LRRS_SSCORR_OG_BIN_PLUS_1  (master)             #
! #             LRRS_SSCORR_OG_MONO_PLUS_1  (master)            #
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


      MODULE lrrs_corrections_plus_1

      PRIVATE
      PUBLIC :: LRRS_SSCORR_OG_BIN_PLUS_1,&
                LRRS_SSCORR_OG_MONO_PLUS_1

      CONTAINS

      SUBROUTINE LRRS_SSCORR_OG_BIN_PLUS_1 &
         ( DO_FULL_RAMAN, DO_UPWELLING, DO_DNWELLING, &
           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_WATERLEAVING, &
           DO_LAMBERTIAN_SURFACE, DO_LAMBERTIAN_FRACTION, DO_DMSCALING, &
           DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, &
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS, &
           NLAYERS, NFINELAYERS, NMOMENTS_INPUT, &
           NPOINTS_INNER, OFFSET_INNER, NPOINTS_OUTER, &
           N_LOUTPUT, LEVELMASK_UP, LEVELMASK_DN, &
           N_USER_STREAMS, N_USER_RELAZMS, &
           SOLAR_ANGLES_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST, &
           LAMBERTIAN_FRACTION, EXACTDB_BRDFUNC, LS_EXACTDB_BRDFUNC, &
           ALBEDOS_RANKED, FLUXES_RANKED, EARTH_RADIUS, HEIGHT_GRID, &
             DELTAU_VERT_INPUT,   OMEGAMOMS_ELASTIC_UNSCALED, &
           L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC_UNSCALED, &
           TRUNC_FACTORS, L_TRUNC_FACTORS, N_RRSBINS, BINMAP, &
           OMEGAMOMS_CABANNES_UNSCALED,L_OMEGAMOMS_CABANNES_UNSCALED, &
           OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSLOSS_UNSCALED, &
           OMEGAMOMS_RRSBIN_UNSCALED,  L_OMEGAMOMS_RRSBIN_UNSCALED, &
           ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN, &
           LC_ELASTIC_SS_UP, LP_ELASTIC_SS_UP, LS_ELASTIC_SS_UP, &
           LC_ELASTIC_SS_DN, LP_ELASTIC_SS_DN, &
           LC_RAMAN_SS_UP,   LP_RAMAN_SS_UP,   LS_RAMAN_SS_UP, &
           LC_RAMAN_SS_DN,   LP_RAMAN_SS_DN, &
           FAIL, MESSAGE)

!  include file of dimensions and numbers

      USE LRRS_PARS
      USE LRRS_CORRECTIONS_1

      IMPLICIT NONE

!  Input
!  -----

!  Flags

      LOGICAL, INTENT(IN) ::          DO_CABANNES_RAMAN
      LOGICAL, INTENT(IN) ::          DO_ENERGY_BALANCING
      LOGICAL, INTENT(IN) ::          DO_UPWELLING
      LOGICAL, INTENT(IN) ::          DO_DNWELLING
      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT(IN) ::          DO_WATERLEAVING
      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_FRACTION
      LOGICAL, INTENT(IN) ::          DO_DMSCALING
      LOGICAL, INTENT(IN) ::          DO_FULL_RAMAN

!  Profile linearization control inputs

      LOGICAL, INTENT(IN) ::          do_profile_wfs
      LOGICAL, INTENT(IN) ::          layer_vary_flag   (max_layers)
      INTEGER, INTENT(IN) ::          layer_vary_number (max_layers)

!  Column linearization control inputs

      LOGICAL, INTENT(IN) ::          do_column_wfs
      INTEGER, INTENT(IN) ::          n_totalcolumn_wfs

!  Surface linearization

      LOGICAL, INTENT(IN) ::          do_surface_wfs

!  number of layers

      INTEGER, INTENT(IN) ::          NLAYERS

!  Number of input phase function Legendre moments
!    THIS IS A BASIC INPUT THAT USER MUST PROVIDE

      INTEGER, INTENT(IN) ::          NMOMENTS_INPUT

!  Number of fine layers for sscorr outgoing

      INTEGER, INTENT(IN) ::          NFINELAYERS

!  Solar beam input geometry

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLES_ADJUST &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  User stream variables

      INTEGER, INTENT(IN) ::          N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_ANGLES_ADJUST ( MAX_USER_STREAMS )

!  User azimuth variables

      INTEGER, INTENT(IN) ::          N_USER_RELAZMS
      REAL(FPK), INTENT(IN) :: USER_RELAZMS_ADJUST &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)

      INTEGER, INTENT(IN) ::          N_LOUTPUT

!  outer/inner wavelength range, and offset for the inner range
!    Depends on the choice of solar spectrum

      INTEGER, INTENT(IN) ::          OFFSET_INNER
      INTEGER, INTENT(IN) ::          NPOINTS_INNER
      INTEGER, INTENT(IN) ::          NPOINTS_OUTER

!  Layer masks for doing integrated source terms

      INTEGER, INTENT(IN) ::          LEVELMASK_UP    (MAX_LOUTPUT)
      INTEGER, INTENT(IN) ::          LEVELMASK_DN    (MAX_LOUTPUT)

!  Earth Radius

      REAL(FPK), INTENT(IN) :: EARTH_RADIUS

!  Height grid

      REAL(FPK), INTENT(IN) :: HEIGHT_GRID ( 0:MAX_LAYERS )

!  Scaled optical depths and truncation factors
!    (If no scaling, trunc factors are zero)

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: TRUNC_FACTORS     ( MAX_LAYERS, MAX_POINTS )

!  Unscaled quantities, Elastic input

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

!  Raman and Cabannes input optical properties (UNSCALED)

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_CABANNES_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSLOSS_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSBIN_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Linearizations of all optical properties

      REAL(FPK), INTENT(IN) :: &
          L_DELTAU_VERT_INPUT ( MAX_PARS, MAX_LAYERS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: &
          L_TRUNC_FACTORS       ( MAX_PARS, MAX_LAYERS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_ELASTIC_UNSCALED &
           ( MAX_PARS, MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_CABANNES_UNSCALED &
           ( MAX_PARS, MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSLOSS_UNSCALED &
               ( MAX_PARS, MAX_LAYERS, 0:2, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSBIN_UNSCALED &
               ( MAX_PARS, MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Bin Mapping quantities (internal derivation)

      INTEGER, INTENT(IN) ::          N_RRSBINS  ( MAX_POINTS )
      INTEGER, INTENT(IN) ::          BINMAP ( MAX_BINS, MAX_POINTS )

!  Exact DB values

      REAL(FPK), INTENT(IN) :: EXACTDB_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

      REAL(FPK), INTENT(IN) :: LS_EXACTDB_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Lambertian fraction

      REAL(FPK), INTENT(IN) :: LAMBERTIAN_FRACTION

!  Albedos and solar fluxes

      REAL(FPK), INTENT(IN) :: FLUXES_RANKED   ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: ALBEDOS_RANKED  ( MAX_POINTS )

!  Output
!  ------

!  single scatter Elastic and Raman Intensity results

      REAL(FPK), INTENT(OUT) :: ELASTIC_SS_UP &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: ELASTIC_SS_DN &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: RAMAN_SS_UP &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: RAMAN_SS_DN &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  single scatter Elastic jacobian results

      REAL(FPK), INTENT(OUT) :: LC_ELASTIC_SS_UP &
       (MAX_PARS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LC_ELASTIC_SS_DN &
       (MAX_PARS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LP_ELASTIC_SS_UP &
       (MAX_PARS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LP_ELASTIC_SS_DN &
       (MAX_PARS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LS_ELASTIC_SS_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  single scatter Raman jacobian results

      REAL(FPK), INTENT(OUT) :: LC_RAMAN_SS_UP &
       (MAX_PARS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LC_RAMAN_SS_DN &
       (MAX_PARS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LP_RAMAN_SS_UP &
       (MAX_PARS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LP_RAMAN_SS_DN &
       (MAX_PARS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LS_RAMAN_SS_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Exception handling

      LOGICAL, INTENT(OUT) ::          fail
      CHARACTER (LEN=*), INTENT(OUT) ::    message

!  Geometry routine inputs and outputs
!  -----------------------------------

!  control

      LOGICAL ::          do_fine
      REAL(FPK) :: alpha_boa, theta_boa, phi_boa

!  main outputs (geometry)

      INTEGER ::          ntraverse(0:max_layers)
      REAL(FPK) :: sunpaths(0:max_layers,max_layers)
      REAL(FPK) :: radii   (0:max_layers)
      REAL(FPK) :: alpha_all  (0:max_layers)

!  Fine level output (geometry)

      INTEGER ::          ntraverse_fine(max_layers,max_fine_layers)
      REAL(FPK) :: &
            sunpaths_fine (max_layers,max_layers,max_fine_layers)
      REAL(FPK) :: &
            radii_fine    (max_layers,max_fine_layers)
      REAL(FPK) :: alpha_fine    (max_layers,max_fine_layers)

!  Other (incidental) geometrical output

      REAL(FPK) :: lospaths(max_layers)
      REAL(FPK) :: theta_all  (0:max_layers)
      REAL(FPK) :: phi_all    (0:max_layers)
      REAL(FPK) :: cosscat_up (0:max_layers)
      REAL(FPK) :: cosscat_dn (0:max_layers)

!  Extinction

      REAL(FPK) :: extinction   ( max_layers, max_points )
      REAL(FPK) :: L_extinction ( max_pars, max_layers, max_points )

!  Saved Legendre polynomials

      REAL(FPK) :: &
         SS_PLEG_UP(MAX_LAYERS,0:MAX_MOMENTS_INPUT), &
         SS_PLEG_DN(MAX_LAYERS,0:MAX_MOMENTS_INPUT)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(FPK) :: TMS(MAX_LAYERS,MAX_POINTS)

!  Linearized TMS factor

      REAL(FPK) :: L_TMS(MAX_PARS,MAX_LAYERS,MAX_POINTS)

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

      REAL(FPK) :: L_ESCAT_UP(MAX_PARS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_ESCAT_DN(MAX_PARS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_CSCAT_UP(MAX_PARS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_CSCAT_DN(MAX_PARS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_GSCAT_UP &
              (MAX_PARS,MAX_LAYERS,MAX_BINS,MAX_POINTS)
      REAL(FPK) :: L_GSCAT_DN &
              (MAX_PARS,MAX_LAYERS,MAX_BINS,MAX_POINTS)

!  Cumulative single scatter source terms

      REAL(FPK) :: ESS_CUMSOURCE_UP(0:MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ESS_CUMSOURCE_DN(0:MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ISS_CUMSOURCE_UP(0:MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: ISS_CUMSOURCE_DN(0:MAX_LAYERS,MAX_POINTS)

!  Linearized Cumulative single scatter source terms

      REAL(FPK) :: L_ESS_CUMSOURCE(MAX_PARS,MAX_POINTS)
      REAL(FPK) :: L_ISS_CUMSOURCE(MAX_PARS,MAX_POINTS)

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

      REAL(FPK) :: L_UP_IMULT &
              (MAX_PARS,0:MAX_LAYERS,MAX_LAYERS,MAX_BINS,MAX_POINTS)
      REAL(FPK) :: L_DN_IMULT &
              (MAX_PARS,0:MAX_LAYERS,MAX_LAYERS,MAX_BINS,MAX_POINTS)

      REAL(FPK) :: L_UP_EMULT &
              (MAX_PARS,0:MAX_LAYERS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_DN_EMULT &
              (MAX_PARS,0:MAX_LAYERS,MAX_LAYERS,MAX_POINTS)

      REAL(FPK) :: &
         L_UP_LOSTRANS(MAX_PARS,MAX_LAYERS,MAX_POINTS), &
         L_DN_LOSTRANS(MAX_PARS,MAX_LAYERS,MAX_POINTS)

!  Solar beam attenuations

      REAL(FPK) :: attn      ( 0:max_layers, max_points )
      REAL(FPK) :: attn_fine &
               ( max_layers, max_fine_layers, max_points )

      REAL(FPK) :: l_attn ( max_pars, 0:max_layers, &
                              0:max_layers, max_points )
      REAL(FPK) :: l_attn_fine ( max_pars, 0:max_layers, &
                              max_layers, max_fine_layers, max_points )

!  local variables
!  ---------------

!  Indices

      INTEGER ::       N, NUT, NSTART, NUT_PREV, NLEVEL, L, W, B
      INTEGER ::       UTA, UM, IA, NC, V, WR, WB, Q, K, KS

!  help variables (REAL(FPK) ::)

      REAL(FPK) :: boa_attn, ssflux, x0_boa, x0_boa_4, ratio
      REAL(FPK) :: GSOURCE, CSOURCE, ESOURCE, ISOURCE, SOURCE
      REAL(FPK) :: HELP, COSSCAT, CONTRIB, PHAS
      REAL(FPK) :: PHAS_E, BS, BETA, T0, T2, P2
      REAL(FPK) :: DF1(0:MAX_MOMENTS_INPUT)
      REAL(FPK) :: DF2(0:MAX_MOMENTS_INPUT), OMEGA_LOCAL
      REAL(FPK) :: REFLEC(MAX_POINTS), LS_REFLEC(MAX_POINTS)
      REAL(FPK) :: BRDF, FACTOR, FL1, FL2

!  Linearization help variables

      REAL(FPK) :: L_O, L_F, L_OF, LSS, LCSS, LGSS, LBSS
      REAL(FPK) :: L_PHAS_E, L_ESOURCE, L_ISOURCE, L_SSCORRECTION
      REAL(FPK) :: L_BOA_ATTN, L_LOSS, L_PHAS, L_CSOURCE, L_GSOURCE

!  Bookkeeping

      LOGICAL ::          do_atmos_wfs
      INTEGER ::          NTYPE_VARY, K_PARAMETERS
      LOGICAL ::          DO_RTSOL_VARY(max_layers)
      INTEGER ::          NPARAMS_VARY(max_layers)
      INTEGER ::          KINDEX(max_layers)

!  Local surface flag

      LOGICAL, PARAMETER :: DO_INCLUDE_SURFACE = .true.

!  Local value

      REAL(FPK), PARAMETER :: ls_albedos_ranked = 1.0d0

! #######################################################
! #######################################################
!          Set up operations
! #######################################################
! #######################################################

!  Ss flux scaling factor

      SSFLUX = one / PI4

!  Fine layering flag

      do_fine = .true.

!   Number of fine layers now a control input (7/25/08)

!  General atmospheric linearization term

      do_atmos_wfs = (do_profile_wfs.or. do_column_wfs )

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
        NPARAMS_VARY(1)  = N_TOTALCOLUMN_WFS
        KINDEX(1) = 0
      ENDIF

!  Floating point numbers for Legendre polynomials

      DO L = 2, NMOMENTS_INPUT
        HELP = DBLE(L)
        DF1(L) = DBLE(2*L-1)/HELP
        DF2(L) = DBLE(L-1)/HELP
      ENDDO

!  Create TMS factors, these get stored
!    Delta-M Scaling introduced April 2005.

      DO N = 1, NLAYERS
        DO W = 1, NPOINTS_INNER
          TMS(N,W) = ONE
          IF ( DO_DMSCALING ) THEN
            WR = W + OFFSET_INNER
            OMEGA_LOCAL = OMEGAMOMS_ELASTIC_UNSCALED(N,0,WR)
            HELP = TRUNC_FACTORS(N,WR) * OMEGA_LOCAL
            TMS(N,W) = TMS(N,W) / (ONE - HELP)
          ENDIF
        ENDDO
      ENDDO

!  Create Linearized TMS factor for each layer
!   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)
!C      SHOULD ONLY BE USED WITH THE ENERGY BALANCING
!  Only if required

      IF ( DO_ATMOS_WFS ) THEN
        DO W = 1, NPOINTS_INNER
          WR = W + OFFSET_INNER
          DO K = 1, NLAYERS
            K_PARAMETERS = LAYER_VARY_NUMBER(K)
            DO Q = 1, K_PARAMETERS
              IF ( DO_DMSCALING ) THEN
                OMEGA_LOCAL = OMEGAMOMS_ELASTIC_UNSCALED(K,0,WR)
                L_O  = L_OMEGAMOMS_ELASTIC_UNSCALED(Q,K,0,WR)
                L_F  = L_TRUNC_FACTORS(Q,K,WR)
                L_OF = L_O * TRUNC_FACTORS(K,WR) + L_F * OMEGA_LOCAL
                L_TMS(Q,K,W) = TMS(K,W) * TMS(K,W) * L_OF
              ELSE
                L_TMS(Q,K,W) = 0.0d0
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

!  Source term calculations
!  ========================

! ######################################################################
! ######################################################################
! ######################################################################
!    Start the main loop over all viewing geometries
!     Use the adjusted values of the angles -------> March 2011
! ######################################################################
! ######################################################################
! ######################################################################

      DO UM = 1, N_USER_STREAMS
        ALPHA_BOA = USER_ANGLES_ADJUST(UM)
        DO IA = 1, N_USER_RELAZMS
          THETA_BOA = SOLAR_ANGLES_ADJUST(UM,IA)
          PHI_BOA   = USER_RELAZMS_ADJUST(UM,IA)
          V = N_USER_RELAZMS * (UM-1) + IA

!  *****************************************************************
!  *****************************************************************
!  Upwelling Source terms: Phase matrix, multipliers, transmittances
!  *****************************************************************
!  *****************************************************************

          IF ( DO_UPWELLING ) THEN

!  Call to geometry routine, path distances etc....

            call outgoing_sphergeom_fine_up_1 &
        ( max_layers, max_fine_layers, &
          do_fine, nlayers, nfinelayers, &
          height_grid, earth_radius, &
          alpha_boa, theta_boa, phi_boa, &
          sunpaths,      radii,      ntraverse,      alpha_all, &
          sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine, &
          lospaths, theta_all, phi_all, cosscat_up, &
          fail, message )

!  Exception handling

            if ( fail ) return

!  debug geometry
!            do n = 0, nlayers
!              write(45,'(i4,101f10.5)')n,(sunpaths(n,v),v=1,nlayers)
!            enddo
!            pause

!  Get the attenuations

            call og_attenuations_plus_1 &
        ( max_layers, max_fine_layers, max_points, max_pars, &
          npoints_outer, nlayers, nfinelayers, &
          do_profile_wfs, layer_vary_flag, layer_vary_number, &
          do_column_wfs, n_totalcolumn_wfs, &
          extinction, L_extinction, &
          sunpaths, ntraverse, sunpaths_fine, ntraverse_fine, &
          attn, attn_fine, l_attn, l_attn_fine )

!         do w = 1, npoints_outer
!           write(56,*)w,attn(23,w),attn(1,w)
!         enddo
!         pause

!  Multipliers, transmittances

            call og_integration_bin_up_plus_1 &
        ( max_layers, max_fine_layers, max_points, max_bins, &
          max_pars, do_full_raman, npoints_inner, offset_inner, &
          n_rrsbins, binmap, nlayers, nfinelayers, &
          do_profile_wfs, layer_vary_flag, layer_vary_number, &
          do_column_wfs,  n_totalcolumn_wfs, &
          extinction, l_extinction, radii, alpha_all, alpha_fine, &
          attn, attn_fine, l_attn, l_attn_fine, &
          up_emult, up_imult, up_lostrans, &
          l_up_emult, l_up_imult, l_up_lostrans )

!  Debug
!              w = 1
!              do n = 1, nlayers
!                write(29,'(2i5,1p4e18.10)')n,w,
!     &           up_lostrans(n,w),up_emult(n,w),
!     &               l_up_lostrans(2,n,w),l_up_emult(2,33,n,w)
!               do b = 1, n_rrsbins(w)
!                write(29,'(2i5,1pe18.10)')
!     &                 w,binmap(b,w),up_imult(n,b,w)
!               enddo
!              enddo
!              pause

!  Debug: Multipliers (all layers), Upwelling OLD WAY
!              call outgoing_integration_old_up
!     i           ( nlayers, extinction, deltau_vert,
!     i             lospaths, sunpaths, radii, ntraverse, alpha_all,
!     o             up_multipliers(1,v), up_lostrans(1,v) )
!              write(46,*)v
!              do n = 1, nlayers
!                write(46,*)up_lostrans(n,v),up_multipliers(n,v)
!              enddo

!  legendre polynomials

            COSSCAT = COSSCAT_UP(NLAYERS)
            SS_PLEG_UP(1,0) = 1.0d0
            SS_PLEG_UP(1,1) = COSSCAT
            DO L = 2, NMOMENTS_INPUT
              SS_PLEG_UP(1,L) = &
                 DF1(L) * SS_PLEG_UP(1,L-1) * COSSCAT  - &
                 DF2(L) * SS_PLEG_UP(1,L-2)
            ENDDO
            P2 = SS_PLEG_UP(1,2)

!  Start the wavelength loop

            DO W = 1, NPOINTS_INNER
              WR = W + OFFSET_INNER

!  Elastic scattering
!  ==================

              DO N = 1, NLAYERS
                PHAS_E = 0.0d0
                DO L = 0, NMOMENTS_INPUT
                  BETA   = OMEGAMOMS_ELASTIC_UNSCALED(N,L,WR)
                  PHAS_E = PHAS_E + BETA * SS_PLEG_UP(1,L)
                ENDDO
                ESCAT_UP(N,W) = PHAS_E * TMS(N,W)
                IF ( LAYER_VARY_FLAG(N) ) THEN
                 DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_PHAS_E = 0.0d0
                  DO L = 0, NMOMENTS_INPUT
                   BETA = L_OMEGAMOMS_ELASTIC_UNSCALED(Q,N,L,WR)
                   L_PHAS_E = L_PHAS_E + BETA * SS_PLEG_UP(1,L)
                  ENDDO
                  L_ESCAT_UP(Q,N,W) = L_PHAS_E *   TMS(N,W) + &
                                        PHAS_E * L_TMS(Q,N,W)
                ENDDO
               ENDIF
              ENDDO

!  Add inelastic scattering source functions if flagged
!  ====================================================

              IF ( DO_FULL_RAMAN ) THEN

!  Cabannes Term
!  -------------

!    Either by Energy Balancing or by Cabannes cross-sections directly

                IF ( DO_ENERGY_BALANCING ) THEN
                  DO N = 1, NLAYERS
                   PHAS = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,W) + &
                          OMEGAMOMS_RRSLOSS_UNSCALED(N,2,W) * P2
                   CSCAT_UP(N,W) = ESCAT_UP(N,W) - PHAS * TMS(N,W)
                   IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                     T0 =  L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,0,W)
                     T2 =  L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,2,W) * P2
                     L_PHAS = T0 + T2
                     L_LOSS = PHAS * L_TMS(Q,N,W) + L_PHAS * TMS(N,W)
                     L_CSCAT_UP(Q,N,W) = L_ESCAT_UP(Q,N,W) - L_LOSS
                    ENDDO
                   ENDIF
                  ENDDO
                ELSE IF ( DO_CABANNES_RAMAN ) THEN
                  DO N = 1, NLAYERS
                   PHAS = 0.0d0
                   DO L = 0, NMOMENTS_INPUT
                    BETA   = OMEGAMOMS_CABANNES_UNSCALED(N,L,W)
                    PHAS = PHAS + BETA * SS_PLEG_UP(1,L)
                   ENDDO
                   CSCAT_UP(N,W) = PHAS * TMS(N,W)
                   IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                     L_PHAS = 0.0d0
                     DO L = 0, NMOMENTS_INPUT
                      BETA   = L_OMEGAMOMS_CABANNES_UNSCALED(Q,N,L,W)
                      L_PHAS = L_PHAS + BETA * SS_PLEG_UP(1,L)
                     ENDDO
                     L_CSCAT_UP(Q,N,W) = L_PHAS *   TMS(N,W) + &
                                           PHAS * L_TMS(Q,N,W)
                    ENDDO
                   ENDIF
                  ENDDO
                ENDIF

!  Gain term (Same for both realizations)
!  --------------------------------------

                DO N = 1, NLAYERS
                  DO B = 1, N_RRSBINS(W)
                    WB = BINMAP(B,W)
                    RATIO   = FLUXES_RANKED(WB)/FLUXES_RANKED(WR)
                    CONTRIB = OMEGAMOMS_RRSBIN_UNSCALED(N,0,B,W) + &
                              OMEGAMOMS_RRSBIN_UNSCALED(N,2,B,W) * P2
                    PHAS    =  CONTRIB * RATIO
                    GSCAT_UP(N,B,W) = PHAS * TMS(N,W)
                    IF ( LAYER_VARY_FLAG(N) ) THEN
                     DO Q = 1, LAYER_VARY_NUMBER(N)
                      T0 = L_OMEGAMOMS_RRSBIN_UNSCALED(Q,N,0,B,W)
                      T2 = L_OMEGAMOMS_RRSBIN_UNSCALED(Q,N,2,B,W) * P2
                      L_PHAS = ( T0 + T2 ) * RATIO
                      L_GSCAT_UP(Q,N,B,W) = L_PHAS *   TMS(N,W) &
                                            + PHAS * L_TMS(Q,N,W)
                     ENDDO
                    ENDIF
                  ENDDO
                ENDDO

!  End inelastic scattering option

              ENDIF

!  End the wavelength loop

            ENDDO

!  End upwelling clause

          ENDIF

!  *******************************************************************
!  *******************************************************************
!  Downwelling Source terms: Phase matrix, multipliers, transmittances
!  *******************************************************************
!  *******************************************************************

          IF ( DO_DNWELLING ) THEN

!  Call to geometry routine, path distances etc....

            call outgoing_sphergeom_fine_dn_1 &
        ( max_layers, max_fine_layers, &
          do_fine, nlayers, nfinelayers, &
          height_grid, earth_radius, &
          alpha_boa, theta_boa, phi_boa, &
          sunpaths,      radii,      ntraverse,      alpha_all, &
          sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine, &
          lospaths, theta_all, phi_all, cosscat_dn, &
          fail, message )

!  Exception handling

            if ( fail ) return

!  debug geometry
!            do n = 0, nlayers
!              write(45,'(i4,101f10.5)')n,(sunpaths(n,v),v=1,nlayers)
!            enddo
!            pause

!  Get the attenuations

            call og_attenuations_plus_1 &
        ( max_layers, max_fine_layers, max_points, max_pars, &
          npoints_outer, nlayers, nfinelayers, &
          do_profile_wfs, layer_vary_flag, layer_vary_number, &
          do_column_wfs, n_totalcolumn_wfs, &
          extinction, L_extinction, &
          sunpaths, ntraverse, sunpaths_fine, ntraverse_fine, &
          attn, attn_fine, l_attn, l_attn_fine )

!         do w = 1, npoints_outer
!           write(56,*)w,attn(23,w),l_attn(1,23,w)
!         enddo
!         pause

!  Multipliers and transmittances

            call og_integration_bin_dn_plus_1 &
        ( max_layers, max_fine_layers, max_points, max_bins, &
          max_pars, do_full_raman, npoints_inner, offset_inner, &
          n_rrsbins, binmap, nlayers, nfinelayers, &
          do_profile_wfs, layer_vary_flag, layer_vary_number, &
          do_column_wfs,  n_totalcolumn_wfs, &
          extinction, l_extinction, radii, alpha_all, alpha_fine, &
          attn, attn_fine, l_attn, l_attn_fine, &
          dn_emult, dn_imult, dn_lostrans, &
          l_dn_emult, l_dn_imult, l_dn_lostrans )

!  Debug
!              write(65,*)v
!              do n = 1, nlayers
!                write(65,'(1p2e18.10)')
!     &           dn_lostrans(n,v),dn_multipliers(n,v)
!              enddo

!  Debug: Multipliers (all layers), Down welling OLD WAY
!              call outgoing_integration_old_dn
!     i           ( nlayers, extinction, deltau_vert,
!     i             lospaths, sunpaths, radii, ntraverse, alpha_all,
!     o             dn_multipliers(1,v), dn_lostrans(1,v) )
!              write(66,*)v
!              do n = 1, nlayers
!                write(66,*)dn_lostrans(n,v),dn_multipliers(n,v)
!              enddo

!  Legendre polynomials

            COSSCAT = COSSCAT_DN(NLAYERS)
            SS_PLEG_DN(1,0) = 1.0d0
            SS_PLEG_DN(1,1) = COSSCAT
            DO L = 2, NMOMENTS_INPUT
              SS_PLEG_DN(1,L) = &
                 DF1(L) * SS_PLEG_DN(1,L-1) * COSSCAT  - &
                 DF2(L) * SS_PLEG_DN(1,L-2)
            ENDDO
            P2 = SS_PLEG_DN(1,2)

!  Start the wavelength loop

            DO W = 1, NPOINTS_INNER
              WR = W + OFFSET_INNER

!  Elastic scatting
!  ================

              DO N = 1, NLAYERS
                PHAS_E = 0.0d0
                DO L = 0, NMOMENTS_INPUT
                  BETA   = OMEGAMOMS_ELASTIC_UNSCALED(N,L,WR)
                  PHAS_E = PHAS_E + BETA * SS_PLEG_DN(1,L)
                ENDDO
                ESCAT_DN(N,W) = PHAS_E * TMS(N,W)
                IF ( LAYER_VARY_FLAG(N) ) THEN
                 DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_PHAS_E = 0.0d0
                  DO L = 0, NMOMENTS_INPUT
                    BETA = L_OMEGAMOMS_ELASTIC_UNSCALED(Q,N,L,WR)
                    L_PHAS_E = L_PHAS_E + BETA * SS_PLEG_DN(1,L)
                  ENDDO
                  L_ESCAT_DN(Q,N,W) = L_PHAS_E *   TMS(N,W) + &
                                        PHAS_E * L_TMS(Q,N,W)
                 ENDDO
                ENDIF
              ENDDO

!  Add inelastic scattering source functions if flagged
!  ====================================================

              IF ( DO_FULL_RAMAN ) THEN

!  Cabannes Term
!  -------------

!    Either by Energy Balancing or by Cabannes cross-sections directly

                IF ( DO_ENERGY_BALANCING ) THEN
                  DO N = 1, NLAYERS
                   PHAS = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,W) + &
                          OMEGAMOMS_RRSLOSS_UNSCALED(N,2,W) * P2
                   CSCAT_DN(N,W) = ESCAT_DN(N,W) - PHAS * TMS(N,W)
                   IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                     T0 =  L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,0,W)
                     T2 =  L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,2,W) * P2
                     L_PHAS = T0 + T2
                     L_LOSS = PHAS * L_TMS(Q,N,W) + L_PHAS * TMS(N,W)
                     L_CSCAT_DN(Q,N,W) = L_ESCAT_DN(Q,N,W) - L_LOSS
                    ENDDO
                   ENDIF
                  ENDDO
                ELSE IF ( DO_CABANNES_RAMAN ) THEN
                  DO N = 1, NLAYERS
                   PHAS = 0.0d0
                   DO L = 0, NMOMENTS_INPUT
                    BETA   = OMEGAMOMS_CABANNES_UNSCALED(N,L,W)
                    PHAS = PHAS + BETA * SS_PLEG_DN(1,L)
                   ENDDO
                   CSCAT_DN(N,W) = PHAS * TMS(N,W)
                   IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                     L_PHAS = 0.0d0
                     DO L = 0, NMOMENTS_INPUT
                      BETA   = L_OMEGAMOMS_CABANNES_UNSCALED(Q,N,L,W)
                      L_PHAS = L_PHAS + BETA * SS_PLEG_DN(1,L)
                     ENDDO
                     L_CSCAT_DN(Q,N,W) = L_PHAS *   TMS(N,W) + &
                                           PHAS * L_TMS(Q,N,W)
                    ENDDO
                   ENDIF
                  ENDDO
                ENDIF

!  Gain term (Same for both realizations)
!  --------------------------------------

                DO N = 1, NLAYERS
                  DO B = 1, N_RRSBINS(W)
                    WB = BINMAP(B,W)
                    RATIO   = FLUXES_RANKED(WB)/FLUXES_RANKED(WR)
                    CONTRIB = OMEGAMOMS_RRSBIN_UNSCALED(N,0,B,W) + &
                              OMEGAMOMS_RRSBIN_UNSCALED(N,2,B,W) * P2
                    PHAS    =  CONTRIB * RATIO
                    GSCAT_DN(N,B,W) = PHAS * TMS(N,W)
                    IF ( LAYER_VARY_FLAG(N) ) THEN
                     DO Q = 1, LAYER_VARY_NUMBER(N)
                      T0 = L_OMEGAMOMS_RRSBIN_UNSCALED(Q,N,0,B,W)
                      T2 = L_OMEGAMOMS_RRSBIN_UNSCALED(Q,N,2,B,W) * P2
                      L_PHAS = ( T0 + T2 ) * RATIO
                      L_GSCAT_DN(Q,N,B,W) = L_PHAS *   TMS(N,W) &
                                            + PHAS * L_TMS(Q,N,W)
                     ENDDO
                    ENDIF
                  ENDDO
                ENDDO

!  End inelastic scattering option

              ENDIF

!  End the wavelength loop

            ENDDO

!  End Downwelling clause

          ENDIF

!  *******************************************************************
!  *******************************************************************
!  Recurrence relation for the UPWELLING intensity
!  *******************************************************************
!  *******************************************************************

          IF ( DO_UPWELLING ) THEN

!  Full BRDF treatment.

            IF ( DO_INCLUDE_SURFACE ) THEN
              IF ( DO_LAMBERTIAN_SURFACE ) THEN
                DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER
                  REFLEC(W) = ALBEDOS_RANKED(WR)
                ENDDO
              ELSE
                BRDF = EXACTDB_BRDFUNC(UM,IA)
                IF ( DO_LAMBERTIAN_FRACTION ) THEN
                  DO W = 1, NPOINTS_INNER
                    WR = W + OFFSET_INNER
                    IF ( DO_WATERLEAVING ) THEN
                      FL1 = ALBEDOS_RANKED(WR)
                      FL2 = ONE
                    ELSE
                      FL1 = LAMBERTIAN_FRACTION * ALBEDOS_RANKED(WR)
                      FL2 = ONE - LAMBERTIAN_FRACTION
                    ENDIF
                    REFLEC(W) = FL1 + FL2 * BRDF
                  ENDDO
                ELSE
                  DO W = 1, NPOINTS_INNER
                    REFLEC(W) = BRDF
                  ENDDO
                ENDIF
              ENDIF
            ENDIF

!  initialize cumulative source terms = F.A.mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux, A = albedo
!    Only require the Stokes total intensity component.
!    This contribution is purely elastic.
!    This line Replaced: X0_BOA = DCOS(SOLAR_ANGLE * DEG_TO_RAD)
!     Introduction of the Adjusted Gometry value, 18 March 2011

            NC =  0
            IF ( DO_INCLUDE_SURFACE ) THEN
              X0_BOA   = DCOS(SOLAR_ANGLES_ADJUST(UM,IA) * DEG_TO_RAD)
              X0_BOA_4 = 4.0d0 * X0_BOA
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
                  ESOURCE = ESCAT_UP(N,W) * UP_EMULT(N,W)
                  ESS_CUMSOURCE_UP(NC,W) =  ESOURCE + &
                       UP_LOSTRANS(N,W) * ESS_CUMSOURCE_UP(NC-1,W)
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO W = 1, NPOINTS_INNER
                    CSOURCE = CSCAT_UP(N,W) * UP_EMULT(N,W)
                    GSOURCE = 0.0d0
                    DO B = 1, N_RRSBINS(W)
                      BS = GSCAT_UP (N,B,W) * UP_IMULT(N,B,W)
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

!  end optical depth loop and Upwelling clause

            ENDDO
          ENDIF

!  *******************************************************************
!  *******************************************************************
!  Recurrence relation for the DOWNWELLING intensity
!  *******************************************************************
!  *******************************************************************

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
                  ESOURCE = ESCAT_DN(N,W) * DN_EMULT(N,W)
                  ESS_CUMSOURCE_DN(NC,W) = ESOURCE + &
                       DN_LOSTRANS(N,W) * ESS_CUMSOURCE_DN(NC-1,W)
                ENDDO
                IF ( DO_FULL_RAMAN ) THEN
                  DO W = 1, NPOINTS_INNER
                    CSOURCE = CSCAT_DN(N,W) * DN_EMULT(N,W)
                    GSOURCE = 0.0d0
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

!  end optical depth loop and Downwelling clause

            ENDDO
          ENDIF

!  *******************************************************************
!  *******************************************************************
!  Recurrence relation for the UPWELLING Atmsopheric Jacobians
!  *******************************************************************
!  *******************************************************************

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
                  L_ESS_CUMSOURCE(Q,W) = 0.0d0
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
                      L_GSOURCE = 0.0d0
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
                      LGSS = 0.0d0
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
!    add Linearization of additional partial layer source term =
!        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)
!  Variations when N > K
!    add Linearization of additional partial layer source term =
!         Exact_Scat * L_Multiplier(k)

!  Profile linearization

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

!  Column linearization

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

!  end optical depth loop

             ENDDO

!  end loop over varying layers

            ENDIF
           ENDDO

!  end Upwelling Jacobian clause

          ENDIF

!  *******************************************************************
!  *******************************************************************
!  Recurrence relation for the DOWNWELLING Atmospheric Jacobians
!  *******************************************************************
!  *******************************************************************

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
                L_ESS_CUMSOURCE(Q,W) = 0.0d0
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

!  If N = K (profile) or KS = 0 (column), there are extra linearizations
!  If N > K, N < K, transmittance + some multiplication linearizations

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
                   LGSS = 0.0d0
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
                   LGSS = 0.0d0
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


!  Profile linearization

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

!  Column linearization

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

!  end optical depth loop

             ENDDO

!  end loop over varying layers

            ENDIF
           ENDDO

!  end Downwelling Jacobian clause

          ENDIF

!  *******************************************************************
!  *******************************************************************
!  Recurrence relation for the UPWELLING Surface Jacobian
!  *******************************************************************
!  *******************************************************************

!  This contribution is purely Elastic

          IF ( DO_SURFACE_WFS .and. DO_UPWELLING ) THEN

!  Finish if no surface (Zero output and go to the next point.

            IF ( .NOT. DO_INCLUDE_SURFACE ) THEN
              DO W = 1, NPOINTS_INNER
                DO UTA = N_LOUTPUT, 1, -1
                  LS_ELASTIC_SS_UP(UTA,V,W) = ZERO
                  LS_RAMAN_SS_UP(UTA,V,W)   = ZERO
                ENDDO
              ENDDO
            ENDIF

!  If there is a surface.....

            IF ( DO_INCLUDE_SURFACE ) THEN

!  Lambertian

              IF ( DO_LAMBERTIAN_SURFACE ) THEN
                DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER
                  LS_REFLEC(W) = LS_ALBEDOS_RANKED * ALBEDOS_RANKED(WR)
                ENDDO
              ENDIF

!  Full BRDF treatment.

              IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
                BRDF = LS_EXACTDB_BRDFUNC(UM,IA)
                IF ( DO_LAMBERTIAN_FRACTION ) THEN
                 DO W = 1, NPOINTS_INNER
                  WR = W + OFFSET_INNER
                  IF ( DO_WATERLEAVING ) THEN
                    FL1 = LS_ALBEDOS_RANKED
                    FL2 = ONE
                  ELSE
                    FL1 = LAMBERTIAN_FRACTION * &
                           LS_ALBEDOS_RANKED * ALBEDOS_RANKED(WR)
                    FL2 = ONE - LAMBERTIAN_FRACTION
                  ENDIF
                  LS_REFLEC(W) = FL1 + FL2 * BRDF
                 ENDDO
                ELSE
                 DO W = 1, NPOINTS_INNER
                  LS_REFLEC(W) = BRDF
                 ENDDO
                ENDIF
              ENDIF

!  initialize cumulative source terms = F.L(A).mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux.
!    L(A) = Linearized surface term (e.g. albedo)
!
              NC =  0
              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
                FACTOR = LS_REFLEC(W)
                BOA_ATTN = X0_BOA_4 * ATTN(NLAYERS,WR)
                ESS_CUMSOURCE_UP(NC,W) = FACTOR * BOA_ATTN
              ENDDO

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
                  DO W = 1, NPOINTS_INNER
                    ESS_CUMSOURCE_UP(NC,W) = &
                         UP_LOSTRANS(N,W) * ESS_CUMSOURCE_UP(NC-1,W)
                  ENDDO
                ENDDO

!  Set output

                DO W = 1, NPOINTS_INNER
                  SOURCE = ESS_CUMSOURCE_UP(NC,W)
                  LS_ELASTIC_SS_UP(UTA,V,W) = SSFLUX * SOURCE
                  IF ( DO_FULL_RAMAN ) &
                       LS_RAMAN_SS_UP(UTA,V,W) = SSFLUX * SOURCE
                ENDDO

!  Check for updating the recursion

                IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                NUT_PREV = NUT

!  end optical depth loop

              ENDDO

!  End inclusion of surface

            ENDIF

!  end Upwelling clause (surface Jacobian)

          ENDIF

!  Finish geometry loop

        enddo
      enddo

!  Finish

      RETURN
      END SUBROUTINE LRRS_SSCORR_OG_BIN_PLUS_1

!

      SUBROUTINE LRRS_SSCORR_OG_MONO_PLUS_1 &
         ( DO_FULL_RAMAN, DO_UPWELLING, DO_DNWELLING, &
           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_WATERLEAVING, &
           DO_LAMBERTIAN_SURFACE, DO_LAMBERTIAN_FRACTION, DO_DMSCALING, &
           DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, &
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS, &
           NLAYERS, NFINELAYERS, NMOMENTS_INPUT, &
           NPOINTS_MONO, W_EXCIT, &
           N_LOUTPUT, LEVELMASK_UP, LEVELMASK_DN, &
           N_USER_STREAMS, N_USER_RELAZMS, &
           SOLAR_ANGLES_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST, &
           LAMBERTIAN_FRACTION, EXACTDB_BRDFUNC, LS_EXACTDB_BRDFUNC, &
           ALBEDOS_RANKED, FLUXES_RANKED, EARTH_RADIUS, HEIGHT_GRID, &
             DELTAU_VERT_INPUT,   OMEGAMOMS_ELASTIC_UNSCALED, &
           L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC_UNSCALED, &
           TRUNC_FACTORS, L_TRUNC_FACTORS, &
           OMEGAMOMS_CABANNES_UNSCALED,L_OMEGAMOMS_CABANNES_UNSCALED, &
           OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSLOSS_UNSCALED, &
           OMEGAMOMS_RRSGAIN_UNSCALED, L_OMEGAMOMS_RRSGAIN_UNSCALED, &
           ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN, &
           LC_ELASTIC_SS_UP, LP_ELASTIC_SS_UP, LS_ELASTIC_SS_UP, &
           LC_ELASTIC_SS_DN, LP_ELASTIC_SS_DN, &
           LC_RAMAN_SS_UP,   LP_RAMAN_SS_UP,   LS_RAMAN_SS_UP, &
           LC_RAMAN_SS_DN,   LP_RAMAN_SS_DN, &
           FAIL, MESSAGE)

!  include file of dimensions and numbers

      USE LRRS_PARS
      USE LRRS_CORRECTIONS_1

      IMPLICIT NONE

!  Input
!  -----

!  Flags

      LOGICAL, INTENT(IN) ::          DO_CABANNES_RAMAN
      LOGICAL, INTENT(IN) ::          DO_ENERGY_BALANCING
      LOGICAL, INTENT(IN) ::          DO_UPWELLING
      LOGICAL, INTENT(IN) ::          DO_DNWELLING
      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT(IN) ::          DO_WATERLEAVING
      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_FRACTION
      LOGICAL, INTENT(IN) ::          DO_DMSCALING
      LOGICAL, INTENT(IN) ::          DO_FULL_RAMAN

!  Profile linearization control inputs

      LOGICAL, INTENT(IN) ::          do_profile_wfs
      LOGICAL, INTENT(IN) ::          layer_vary_flag   (max_layers)
      INTEGER, INTENT(IN) ::          layer_vary_number (max_layers)

!  Column linearization control inputs

      LOGICAL, INTENT(IN) ::          do_column_wfs
      INTEGER, INTENT(IN) ::          n_totalcolumn_wfs

!  Surface linearization

      LOGICAL, INTENT(IN) ::          do_surface_wfs

!  number of layers

      INTEGER, INTENT(IN) ::          NLAYERS

!  Number of input phase function Legendre moments
!    THIS IS A BASIC INPUT THAT USER MUST PROVIDE

      INTEGER, INTENT(IN) ::          NMOMENTS_INPUT

!  Number of fine layers for sscorr outgoing

      INTEGER, INTENT(IN) ::          NFINELAYERS

!  Solar beam input geometry

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLES_ADJUST &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  User stream variables

      INTEGER, INTENT(IN) ::          N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_ANGLES_ADJUST ( MAX_USER_STREAMS )

!  User azimuth variables

      INTEGER, INTENT(IN) ::          N_USER_RELAZMS
      REAL(FPK), INTENT(IN) :: USER_RELAZMS_ADJUST &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)

      INTEGER, INTENT(IN) ::          N_LOUTPUT

!  Number of points, excitation index

      INTEGER, INTENT(IN) ::          NPOINTS_MONO, W_EXCIT

!  Layer masks for doing integrated source terms

      INTEGER, INTENT(IN) ::          LEVELMASK_UP    (MAX_LOUTPUT)
      INTEGER, INTENT(IN) ::          LEVELMASK_DN    (MAX_LOUTPUT)

!  Earth Radius

      REAL(FPK), INTENT(IN) :: EARTH_RADIUS

!  Height grid

      REAL(FPK), INTENT(IN) :: HEIGHT_GRID ( 0:MAX_LAYERS )

!  Scaled optical depths and truncation factors
!    (If no scaling, trunc factors are zero)

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: TRUNC_FACTORS     ( MAX_LAYERS, MAX_POINTS )

!  Unscaled quantities, Elastic input

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

!  Raman and Cabannes input optical properties (UNSCALED)

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_CABANNES_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSLOSS_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSGAIN_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Linearizations of all optical properties

      REAL(FPK), INTENT(IN) :: &
          L_DELTAU_VERT_INPUT ( MAX_PARS, MAX_LAYERS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: &
          L_TRUNC_FACTORS       ( MAX_PARS, MAX_LAYERS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_ELASTIC_UNSCALED &
           ( MAX_PARS, MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_CABANNES_UNSCALED &
           ( MAX_PARS, MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSLOSS_UNSCALED &
               ( MAX_PARS, MAX_LAYERS, 0:2, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_RRSGAIN_UNSCALED &
               ( MAX_PARS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Exact DB values

      REAL(FPK), INTENT(IN) :: EXACTDB_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

      REAL(FPK), INTENT(IN) :: LS_EXACTDB_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Lambertian fraction

      REAL(FPK), INTENT(IN) :: LAMBERTIAN_FRACTION

!  Albedos and solar fluxes

      REAL(FPK), INTENT(IN) :: FLUXES_RANKED   ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: ALBEDOS_RANKED  ( MAX_POINTS )

!  Output
!  ------

!  single scatter Elastic and Raman Intensity results

      REAL(FPK), INTENT(OUT) :: ELASTIC_SS_UP &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: ELASTIC_SS_DN &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: RAMAN_SS_UP &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: RAMAN_SS_DN &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  single scatter Elastic jacobian results

      REAL(FPK), INTENT(OUT) :: LC_ELASTIC_SS_UP &
       (MAX_PARS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LC_ELASTIC_SS_DN &
       (MAX_PARS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LP_ELASTIC_SS_UP &
       (MAX_PARS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LP_ELASTIC_SS_DN &
       (MAX_PARS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LS_ELASTIC_SS_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  single scatter Raman jacobian results

      REAL(FPK), INTENT(OUT) :: LC_RAMAN_SS_UP &
       (MAX_PARS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LC_RAMAN_SS_DN &
       (MAX_PARS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LP_RAMAN_SS_UP &
       (MAX_PARS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LP_RAMAN_SS_DN &
       (MAX_PARS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(OUT) :: LS_RAMAN_SS_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Exception handling

      LOGICAL, INTENT(OUT) ::          fail
      CHARACTER (LEN=*), INTENT(OUT) ::    message

!  Geometry routine inputs and outputs
!  -----------------------------------

!  control

      LOGICAL ::          do_fine
      REAL(FPK) :: alpha_boa, theta_boa, phi_boa

!  main outputs (geometry)

      INTEGER ::          ntraverse(0:max_layers)
      REAL(FPK) :: sunpaths(0:max_layers,max_layers)
      REAL(FPK) :: radii   (0:max_layers)
      REAL(FPK) :: alpha_all  (0:max_layers)

!  Fine level output (geometry)

      INTEGER ::          ntraverse_fine(max_layers,max_fine_layers)
      REAL(FPK) :: &
            sunpaths_fine (max_layers,max_layers,max_fine_layers)
      REAL(FPK) :: &
            radii_fine    (max_layers,max_fine_layers)
      REAL(FPK) :: alpha_fine    (max_layers,max_fine_layers)


!  Other (incidental) geometrical output

      REAL(FPK) :: lospaths(max_layers)
      REAL(FPK) :: theta_all  (0:max_layers)
      REAL(FPK) :: phi_all    (0:max_layers)
      REAL(FPK) :: cosscat_up (0:max_layers)
      REAL(FPK) :: cosscat_dn (0:max_layers)

!  Extinction

      REAL(FPK) :: extinction   ( max_layers, max_points )
      REAL(FPK) :: L_extinction ( max_pars, max_layers, max_points )

!  Saved Legendre polynomials

      REAL(FPK) :: &
         SS_PLEG_UP(MAX_LAYERS,0:MAX_MOMENTS_INPUT), &
         SS_PLEG_DN(MAX_LAYERS,0:MAX_MOMENTS_INPUT)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(FPK) :: TMS(MAX_LAYERS)

!  Linearized TMS factor

      REAL(FPK) :: L_TMS(MAX_PARS,MAX_LAYERS)

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

      REAL(FPK) :: L_ESCAT_UP(MAX_PARS,MAX_LAYERS)
      REAL(FPK) :: L_ESCAT_DN(MAX_PARS,MAX_LAYERS)
      REAL(FPK) :: L_CSCAT_UP(MAX_PARS,MAX_LAYERS)
      REAL(FPK) :: L_CSCAT_DN(MAX_PARS,MAX_LAYERS)
      REAL(FPK) :: L_GSCAT_UP(MAX_PARS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_GSCAT_DN(MAX_PARS,MAX_LAYERS,MAX_POINTS)

!  Cumulative single scatter source terms

      REAL(FPK) :: ESS_CUMSOURCE_UP(0:MAX_LAYERS)
      REAL(FPK) :: ESS_CUMSOURCE_DN(0:MAX_LAYERS)
      REAL(FPK) :: ISS_CUMSOURCE_UP(0:MAX_LAYERS)
      REAL(FPK) :: ISS_CUMSOURCE_DN(0:MAX_LAYERS)

!  Linearized Cumulative single scatter source terms

      REAL(FPK) :: L_ESS_CUMSOURCE(MAX_PARS)
      REAL(FPK) :: L_ISS_CUMSOURCE(MAX_PARS)

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

      REAL(FPK) :: L_UP_IMULT &
              (MAX_PARS,0:MAX_LAYERS,MAX_LAYERS,MAX_POINTS)
      REAL(FPK) :: L_DN_IMULT &
              (MAX_PARS,0:MAX_LAYERS,MAX_LAYERS,MAX_POINTS)

      REAL(FPK) :: L_UP_EMULT &
              (MAX_PARS,0:MAX_LAYERS,MAX_LAYERS)
      REAL(FPK) :: L_DN_EMULT &
              (MAX_PARS,0:MAX_LAYERS,MAX_LAYERS)

      REAL(FPK) :: &
         L_UP_LOSTRANS(MAX_PARS,MAX_LAYERS), &
         L_DN_LOSTRANS(MAX_PARS,MAX_LAYERS)

!  Solar beam attenuations

      REAL(FPK) :: attn      ( 0:max_layers, max_points )
      REAL(FPK) :: attn_fine &
               ( max_layers, max_fine_layers, max_points )

      REAL(FPK) :: l_attn ( max_pars, 0:max_layers, &
                              0:max_layers, max_points )
      REAL(FPK) :: l_attn_fine ( max_pars, 0:max_layers, &
                              max_layers, max_fine_layers, max_points )

!  local variables
!  ---------------

!  Indices

      INTEGER ::       N, NUT, NSTART, NUT_PREV, NLEVEL, L, W
      INTEGER ::       UTA, UM, IA, NC, V, Q, K, KS

!  help variables (REAL(FPK) ::)

      REAL(FPK) :: boa_attn, ssflux, x0_boa, x0_boa_4, ratio
      REAL(FPK) :: GSOURCE, CSOURCE, ESOURCE, ISOURCE, SOURCE
      REAL(FPK) :: HELP, COSSCAT, CONTRIB, PHAS
      REAL(FPK) :: PHAS_E, BETA, T0, T2, P2
      REAL(FPK) :: DF1(0:MAX_MOMENTS_INPUT)
      REAL(FPK) :: DF2(0:MAX_MOMENTS_INPUT), OMEGA_LOCAL
      REAL(FPK) :: REFLEC, LS_REFLEC
      REAL(FPK) :: BRDF, FACTOR, FL1, FL2

!  Linearization help variables

      REAL(FPK) :: L_O, L_F, L_OF, L_LOSS, LSS, LCSS, LGSS, LBSS
      REAL(FPK) :: L_PHAS_E, L_ESOURCE, L_ISOURCE, L_SSCORRECTION
      REAL(FPK) :: L_BOA_ATTN, L_PHAS, L_CSOURCE, L_GSOURCE

!  Bookkeeping

      LOGICAL ::          do_atmos_wfs
      INTEGER ::          NTYPE_VARY, K_PARAMETERS
      LOGICAL ::          DO_RTSOL_VARY(max_layers)
      INTEGER ::          NPARAMS_VARY(max_layers)
      INTEGER ::          KINDEX(max_layers)

!  Local surface flag

      LOGICAL, PARAMETER :: DO_INCLUDE_SURFACE = .true.

!  Local value

      REAL(FPK), PARAMETER :: ls_albedos_ranked = 1.0d0

! #######################################################
! #######################################################
!          Set up operations
! #######################################################
! #######################################################

!  Ss flux scaling factor

      SSFLUX = one / PI4

!  Fine layering flag

      do_fine = .true.

!   Number of fine layers now a control input (7/25/08)

!  General atmospheric linearization term

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
        NPARAMS_VARY(1)  = N_TOTALCOLUMN_WFS
        KINDEX(1) = 0
      ENDIF

!  Floating point numbers for Legendre polynomials

      DO L = 2, NMOMENTS_INPUT
        HELP = DBLE(L)
        DF1(L) = DBLE(2*L-1)/HELP
        DF2(L) = DBLE(L-1)/HELP
      ENDDO

!  Create TMS factors, these get stored
!    Delta-M Scaling introduced April 2005.

      DO N = 1, NLAYERS
        TMS(N) = ONE
        IF ( DO_DMSCALING ) THEN
          OMEGA_LOCAL = OMEGAMOMS_ELASTIC_UNSCALED(N,0,W_EXCIT)
          HELP = TRUNC_FACTORS(N,W_EXCIT) * OMEGA_LOCAL
          TMS(N) = TMS(N) / (ONE - HELP)
        ENDIF
      ENDDO

!  Create Linearized TMS factor for each layer
!   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)
!C      SHOULD ONLY BE USED WITH THE ENERGY BALANCING
!  Only if required

      IF ( DO_ATMOS_WFS ) THEN
        DO K = 1, NLAYERS
          K_PARAMETERS = LAYER_VARY_NUMBER(K)
          DO Q = 1, K_PARAMETERS
            IF ( DO_DMSCALING ) THEN
              OMEGA_LOCAL = OMEGAMOMS_ELASTIC_UNSCALED(K,0,W_EXCIT)
              L_O  = L_OMEGAMOMS_ELASTIC_UNSCALED(Q,K,0,W_EXCIT)
              L_F  = L_TRUNC_FACTORS(Q,K,W_EXCIT)
              L_OF = L_O * TRUNC_FACTORS(K,W_EXCIT) + L_F * OMEGA_LOCAL
              L_TMS(Q,K) = TMS(K) * TMS(K) * L_OF
            ELSE
              L_TMS(Q,K) = 0.0d0
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

!  Source term calculations
!  ========================

! ######################################################################
! ######################################################################
! ######################################################################
!  Start the main loop over all viewing geometries
!   Use the adjusted values of the angles -------> March 2011
! ######################################################################
! ######################################################################
! ######################################################################

      DO UM = 1, N_USER_STREAMS
        ALPHA_BOA = USER_ANGLES_ADJUST(UM)
        DO IA = 1, N_USER_RELAZMS
          THETA_BOA = SOLAR_ANGLES_ADJUST(UM,IA)
          PHI_BOA   = USER_RELAZMS_ADJUST(UM,IA)
          V = N_USER_RELAZMS * (UM-1) + IA

!  *****************************************************************
!  *****************************************************************
!  Upwelling Source terms: Phase matrix, multipliers, transmittances
!  *****************************************************************
!  *****************************************************************

          IF ( DO_UPWELLING ) THEN

!  Call to geometry routine, path distances etc....

            call outgoing_sphergeom_fine_up_1 &
        ( max_layers, max_fine_layers, do_fine, nlayers, nfinelayers, &
          height_grid, earth_radius, alpha_boa, theta_boa, phi_boa, &
          sunpaths,      radii,      ntraverse,      alpha_all, &
          sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine, &
          lospaths, theta_all, phi_all, cosscat_up, &
          fail, message )

!  Exception handling

            if ( fail ) return

!  debug geometry
!            do n = 0, nlayers
!              write(45,'(i4,101f10.5)')n,(sunpaths(n,v),v=1,nlayers)
!            enddo
!            pause

!  Get the attenuations

            call og_attenuations_plus_1 &
        ( max_layers, max_fine_layers, max_points, max_pars, &
          npoints_mono, nlayers, nfinelayers, &
          do_profile_wfs, layer_vary_flag, layer_vary_number, &
          do_column_wfs, n_totalcolumn_wfs, &
          extinction, L_extinction, &
          sunpaths, ntraverse, sunpaths_fine, ntraverse_fine, &
          attn, attn_fine, l_attn, l_attn_fine )

!         do w = 1, npoints_outer
!           write(56,*)w,attn(23,w),attn(1,w)
!         enddo
!         pause

!  Multipliers, transmittances

            call og_integration_mono_up_plus_1 &
        ( max_layers, max_fine_layers, max_points, max_pars, &
          do_full_raman,  npoints_mono, w_excit, nlayers, nfinelayers, &
          do_profile_wfs, layer_vary_flag, layer_vary_number, &
          do_column_wfs,  n_totalcolumn_wfs, &
          extinction, l_extinction, radii, alpha_all, alpha_fine, &
          attn, attn_fine, l_attn, l_attn_fine, &
          up_emult, up_imult, up_lostrans, &
          l_up_emult, l_up_imult, l_up_lostrans )

!  Debug
!              do n = 1, nlayers
!                write(29,'(i5,1p6e18.10)')n,
!     &           up_lostrans(n),up_emult(n), up_imult(n,45),
!     &    l_up_lostrans(1,n),l_up_emult(1,2,n),l_up_imult(1,2,n,45)
!              enddo
!              pause

!  legendre polynomials

            COSSCAT = COSSCAT_UP(NLAYERS)
            SS_PLEG_UP(1,0) = 1.0d0
            SS_PLEG_UP(1,1) = COSSCAT
            DO L = 2, NMOMENTS_INPUT
              SS_PLEG_UP(1,L) = &
                 DF1(L) * SS_PLEG_UP(1,L-1) * COSSCAT  - &
                 DF2(L) * SS_PLEG_UP(1,L-2)
            ENDDO
            P2 = SS_PLEG_UP(1,2)

!  Elastic scatting Phase functions (multiplied by TMS factor)
!  --------------------------------

            DO N = 1, NLAYERS
              PHAS_E = ZERO
              DO L = 0, NMOMENTS_INPUT
                BETA   = OMEGAMOMS_ELASTIC_UNSCALED(N,L,W_EXCIT)
                PHAS_E = PHAS_E + BETA * SS_PLEG_UP(1,L)
              ENDDO
              ESCAT_UP(N) = PHAS_E * TMS(N)
              IF ( LAYER_VARY_FLAG(N) ) THEN
               DO Q = 1, LAYER_VARY_NUMBER(N)
                L_PHAS_E = 0.0d0
                DO L = 0, NMOMENTS_INPUT
                  BETA = L_OMEGAMOMS_ELASTIC_UNSCALED(Q,N,L,W_EXCIT)
                  L_PHAS_E = L_PHAS_E + BETA * SS_PLEG_UP(1,L)
                ENDDO
                L_ESCAT_UP(Q,N) = L_PHAS_E *   TMS(N) + &
                                    PHAS_E * L_TMS(Q,N)
               ENDDO
             ENDIF
            ENDDO

!  Add inelastic scatting Phase functions if flagged
!  -------------------------------------------------

            IF ( DO_FULL_RAMAN ) THEN

!  Elastic Term (eiher by energy balancing or Cabannes-Raman)
!  ------------

              IF ( DO_ENERGY_BALANCING ) THEN
                DO N = 1, NLAYERS
                  PHAS = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,1) + &
                         OMEGAMOMS_RRSLOSS_UNSCALED(N,2,1) * P2
                  CSCAT_UP(N) = ESCAT_UP(N) - PHAS * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      T0 =  L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,0,1)
                      T2 =  L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,2,1) * P2
                      L_PHAS = T0 + T2
                      L_LOSS = PHAS * L_TMS(Q,N) + L_PHAS * TMS(N)
                      L_CSCAT_UP(Q,N) = L_ESCAT_UP(Q,N) - L_LOSS
                    ENDDO
                  ENDIF
                ENDDO
              ELSE IF ( DO_CABANNES_RAMAN ) THEN
                DO N = 1, NLAYERS
                  PHAS = ZERO
                  DO L = 0, NMOMENTS_INPUT
                    BETA    = OMEGAMOMS_CABANNES_UNSCALED(N,L,1)
                    PHAS = PHAS + BETA * SS_PLEG_UP(1,L)
                  ENDDO
                  CSCAT_UP(N) = PHAS * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      L_PHAS = ZERO
                      DO L = 0, NMOMENTS_INPUT
                        BETA   = L_OMEGAMOMS_CABANNES_UNSCALED(Q,N,L,1)
                        L_PHAS = L_PHAS + BETA * SS_PLEG_UP(1,L)
                      ENDDO
                      L_CSCAT_UP(Q,N) = L_PHAS *   TMS(N) + &
                                          PHAS * L_TMS(Q,N)
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF

!  gain term (Common to both approaches)

              DO N = 1, NLAYERS
                DO W = 1, NPOINTS_MONO
                  RATIO   = FLUXES_RANKED(W)/FLUXES_RANKED(W_EXCIT)
                  CONTRIB = OMEGAMOMS_RRSGAIN_UNSCALED(N,0,W) + &
                            OMEGAMOMS_RRSGAIN_UNSCALED(N,2,W) * P2
                  PHAS    =  CONTRIB * RATIO
                  GSCAT_UP(N,W) = PHAS * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      T0 = L_OMEGAMOMS_RRSGAIN_UNSCALED(Q,N,0,W)
                      T2 = L_OMEGAMOMS_RRSGAIN_UNSCALED(Q,N,2,W) * P2
                      L_PHAS = ( T0 + T2 ) * RATIO
                      L_GSCAT_UP(Q,N,W) = L_PHAS *   TMS(N) &
                                          + PHAS * L_TMS(Q,N)
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO

!  End inelastic scattering option

            ENDIF

!  End upwelling clause

          ENDIF

!  *******************************************************************
!  *******************************************************************
!  Downwelling Source terms: Phase matrix, multipliers, transmittances
!  *******************************************************************
!  *******************************************************************

          IF ( DO_DNWELLING ) THEN

!  Call to geometry routine, path distances etc....

            call outgoing_sphergeom_fine_dn_1 &
        ( max_layers, max_fine_layers, do_fine, nlayers, nfinelayers, &
          height_grid, earth_radius, alpha_boa, theta_boa, phi_boa, &
          sunpaths,      radii,      ntraverse,      alpha_all, &
          sunpaths_fine, radii_fine, ntraverse_fine, alpha_fine, &
          lospaths, theta_all, phi_all, cosscat_dn, &
          fail, message )

!  Exception handling

            if ( fail ) return

!  debug geometry
!            do n = 0, nlayers
!              write(45,'(i4,101f10.5)')n,(sunpaths(n,v),v=1,nlayers)
!            enddo
!            pause

!  Get the attenuations

            call og_attenuations_plus_1 &
        ( max_layers, max_fine_layers, max_points, max_pars, &
          npoints_mono, nlayers, nfinelayers, &
          do_profile_wfs, layer_vary_flag, layer_vary_number, &
          do_column_wfs, n_totalcolumn_wfs, &
          extinction, L_extinction, &
          sunpaths, ntraverse, sunpaths_fine, ntraverse_fine, &
          attn, attn_fine, l_attn, l_attn_fine )

!         do w = 1, npoints_outer
!           write(56,*)w,attn(23,w),attn(1,w)
!         enddo
!         pause

!  Multipliers, transmittances

            call og_integration_mono_dn_plus_1 &
        ( max_layers, max_fine_layers, max_points, max_pars, &
          do_full_raman,  npoints_mono, w_excit, nlayers, nfinelayers, &
          do_profile_wfs, layer_vary_flag, layer_vary_number, &
          do_column_wfs,  n_totalcolumn_wfs, &
          extinction, l_extinction, radii, alpha_all, alpha_fine, &
          attn, attn_fine, l_attn, l_attn_fine, &
          dn_emult, dn_imult, dn_lostrans, &
          l_dn_emult, l_dn_imult, l_dn_lostrans )

!  Debug
!              do n = 1, nlayers
!                write(30,'(i5,1p6e18.10)')n,
!     &           dn_lostrans(n),dn_emult(n), dn_imult(n),
!     &      l_dn_lostrans(2,n),l_dn_emult(2,33,n),l_dn_imult(2,33,n)
!              enddo
!              pause

!  Legendre polynomials

            COSSCAT = COSSCAT_DN(NLAYERS)
            SS_PLEG_DN(1,0) = ONE
            SS_PLEG_DN(1,1) = COSSCAT
            DO L = 2, NMOMENTS_INPUT
            SS_PLEG_DN(1,L) = &
               DF1(L) * SS_PLEG_DN(1,L-1) * COSSCAT  - &
               DF2(L) * SS_PLEG_DN(1,L-2)
            ENDDO
            P2 = SS_PLEG_DN(1,2)

!  Elastic scatting Phase functions (multiplied by TMS factor)
!  --------------------------------

            DO N = 1, NLAYERS
              PHAS_E = ZERO
              DO L = 0, NMOMENTS_INPUT
                BETA   = OMEGAMOMS_ELASTIC_UNSCALED(N,L,W_EXCIT)
                PHAS_E = PHAS_E + BETA * SS_PLEG_DN(1,L)
              ENDDO
              ESCAT_DN(N) = PHAS_E * TMS(N)
              IF ( LAYER_VARY_FLAG(N) ) THEN
               DO Q = 1, LAYER_VARY_NUMBER(N)
                L_PHAS_E = 0.0d0
                DO L = 0, NMOMENTS_INPUT
                  BETA = L_OMEGAMOMS_ELASTIC_UNSCALED(Q,N,L,W_EXCIT)
                  L_PHAS_E = L_PHAS_E + BETA * SS_PLEG_DN(1,L)
                ENDDO
                L_ESCAT_DN(Q,N) = L_PHAS_E *   TMS(N) + &
                                    PHAS_E * L_TMS(Q,N)
               ENDDO
             ENDIF
            ENDDO

!  Add inelastic scatting Phase functions if flagged
!  -------------------------------------------------

            IF ( DO_FULL_RAMAN ) THEN

!  Elastic Term (eiher by energy balancing or Cabannes-Raman)
!  ------------

              IF ( DO_ENERGY_BALANCING ) THEN
                DO N = 1, NLAYERS
                  PHAS = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,1) + &
                         OMEGAMOMS_RRSLOSS_UNSCALED(N,2,1) * P2
                  CSCAT_DN(N) = ESCAT_DN(N) - PHAS * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      T0 =  L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,0,1)
                      T2 =  L_OMEGAMOMS_RRSLOSS_UNSCALED(Q,N,2,1) * P2
                      L_PHAS = T0 + T2
                      L_LOSS = PHAS * L_TMS(Q,N) + L_PHAS * TMS(N)
                      L_CSCAT_DN(Q,N) = L_ESCAT_DN(Q,N) - L_LOSS
                    ENDDO
                  ENDIF
                ENDDO
              ELSE IF ( DO_CABANNES_RAMAN ) THEN
                DO N = 1, NLAYERS
                  PHAS = ZERO
                  DO L = 0, NMOMENTS_INPUT
                    BETA    = OMEGAMOMS_CABANNES_UNSCALED(N,L,1)
                    PHAS = PHAS + BETA * SS_PLEG_DN(1,L)
                  ENDDO
                  CSCAT_DN(N) = PHAS * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      L_PHAS = ZERO
                      DO L = 0, NMOMENTS_INPUT
                        BETA   = L_OMEGAMOMS_CABANNES_UNSCALED(Q,N,L,1)
                        L_PHAS = L_PHAS + BETA * SS_PLEG_DN(1,L)
                      ENDDO
                      L_CSCAT_DN(Q,N) = L_PHAS *   TMS(N) + &
                                          PHAS * L_TMS(Q,N)
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF

!  gain term (Common to both approaches)

              DO N = 1, NLAYERS
                DO W = 1, NPOINTS_MONO
                  RATIO   = FLUXES_RANKED(W)/FLUXES_RANKED(W_EXCIT)
                  CONTRIB = OMEGAMOMS_RRSGAIN_UNSCALED(N,0,W) + &
                            OMEGAMOMS_RRSGAIN_UNSCALED(N,2,W) * P2
                  PHAS    =  CONTRIB * RATIO
                  GSCAT_DN(N,W) = PHAS * TMS(N)
                  IF ( LAYER_VARY_FLAG(N) ) THEN
                    DO Q = 1, LAYER_VARY_NUMBER(N)
                      T0 = L_OMEGAMOMS_RRSGAIN_UNSCALED(Q,N,0,W)
                      T2 = L_OMEGAMOMS_RRSGAIN_UNSCALED(Q,N,2,W) * P2
                      L_PHAS = ( T0 + T2 ) * RATIO
                      L_GSCAT_DN(Q,N,W) = L_PHAS *   TMS(N) &
                                          + PHAS * L_TMS(Q,N)
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO

!  End inelastic scattering option

            ENDIF

!  End Downwelling clause

          ENDIF

!  *******************************************************************
!  *******************************************************************
!  Recurrence relation for the UPWELLING intensity
!  *******************************************************************
!  *******************************************************************

          IF ( DO_UPWELLING ) THEN

!  Full BRDF treatment.

            IF ( DO_INCLUDE_SURFACE ) THEN
              IF ( DO_LAMBERTIAN_SURFACE ) THEN
                REFLEC = ALBEDOS_RANKED(W_EXCIT)
              ELSE
                BRDF = EXACTDB_BRDFUNC(UM,IA)
!                BRDF = BRDF_FACTOR * EXACTDB_BRDFUNC(UM,IA)
                IF ( DO_LAMBERTIAN_FRACTION ) THEN
                  IF ( DO_WATERLEAVING ) THEN
                    FL1 = ALBEDOS_RANKED(W_EXCIT)
                    FL2 = ONE
                  ELSE
                    FL1 = LAMBERTIAN_FRACTION * ALBEDOS_RANKED(W_EXCIT)
                    FL2 = ONE - LAMBERTIAN_FRACTION
                  ENDIF
                  REFLEC = FL1 + FL2 * BRDF
                ELSE
                  REFLEC = BRDF
                ENDIF
              ENDIF
            ENDIF

!  initialize cumulative source terms = F.A.mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux, A = albedo
!    Only require the total intensity component.
!    This contribution is purely elastic.
!    This line Replaced: X0_BOA = DCOS(SOLAR_ANGLE * DEG_TO_RAD)
!     Introduction of the Adjusted Gometry value, 18 March 2011

            NC =  0
            IF ( DO_INCLUDE_SURFACE ) THEN
              X0_BOA   = DCOS(SOLAR_ANGLES_ADJUST(UM,IA) * DEG_TO_RAD)
              X0_BOA_4 = 4.0d0 * X0_BOA
              FACTOR = REFLEC
              BOA_ATTN = X0_BOA_4 * ATTN(NLAYERS,W_EXCIT)
              ESS_CUMSOURCE_UP(NC) = FACTOR * BOA_ATTN
            ELSE
              ESS_CUMSOURCE_UP(NC) = ZERO
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

!  *******************************************************************
!  *******************************************************************
!  Recurrence relation for the UPWELLING Atmospheric Jacobians
!  *******************************************************************
!  *******************************************************************

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
                L_ESS_CUMSOURCE(Q) = 0.0d0
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
                     L_GSOURCE = 0.0d0
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
                     LGSS = 0.0d0
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
!    add Linearization of additional partial layer source term =
!        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)
!  Variations when N > K
!    add Linearization of additional partial layer source term =
!         Exact_Scat * L_Multiplier(k)

!  Profile linearization

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

!  Column linearization

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

!  *******************************************************************
!  *******************************************************************
!  Recurrence relation for the DOWNWELLING Atmospheric Jacobians
!  *******************************************************************
!  *******************************************************************

          IF ( DO_DNWELLING .AND. DO_ATMOS_WFS ) THEN

!  Start the main layer variation loop
!  -----------------------------------

           DO K = 1, NTYPE_VARY
            IF ( DO_RTSOL_VARY(K) ) THEN
             K_PARAMETERS = NPARAMS_VARY(K)
             KS = KINDEX(K)

!  Initialise elastic source term

             DO Q = 1, K_PARAMETERS
               L_ESS_CUMSOURCE(Q) = 0.0d0
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

!  If N = K (profile) or KS = 0 (column), there are extra linearizations
!  If N > K, N < K, transmittance + some multiplication linearizations

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
                     LGSS = 0.0d0
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
                     LGSS = 0.0d0
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

!  Profile linearization

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

!  Column linearization

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

!  *******************************************************************
!  *******************************************************************
!  Recurrence relation for the UPWELLING Surface Jacobian
!  *******************************************************************
!  *******************************************************************

!  This contribution is purely Elastic

          IF ( DO_SURFACE_WFS .and. DO_UPWELLING ) THEN

!  Finish if no surface (Zero output and go to the next point.

            IF ( .NOT. DO_INCLUDE_SURFACE ) THEN
              DO UTA = N_LOUTPUT, 1, -1
                LS_ELASTIC_SS_UP(UTA,V,1) = ZERO
                LS_RAMAN_SS_UP(UTA,V,1)   = ZERO
              ENDDO
            ENDIF

!  If there is a surface.....

            IF ( DO_INCLUDE_SURFACE ) THEN

!  Lambertian

              IF ( DO_LAMBERTIAN_SURFACE ) THEN
                LS_REFLEC = LS_ALBEDOS_RANKED * ALBEDOS_RANKED(W_EXCIT)
              ENDIF

!  Full BRDF treatment.

              IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
                BRDF = LS_EXACTDB_BRDFUNC(UM,IA)
                IF ( DO_LAMBERTIAN_FRACTION ) THEN
                  IF ( DO_WATERLEAVING ) THEN
                    FL1 = LS_ALBEDOS_RANKED*ALBEDOS_RANKED(W_EXCIT)
                    FL2 = ONE
                  ELSE
                    FL1 = LAMBERTIAN_FRACTION * &
                           LS_ALBEDOS_RANKED * ALBEDOS_RANKED(W_EXCIT)
                    FL2 = ONE - LAMBERTIAN_FRACTION
                  ENDIF
                  LS_REFLEC = FL1 + FL2 * BRDF
                ELSE
                  LS_REFLEC = BRDF
                ENDIF
              ENDIF

!  initialize cumulative source terms = F.L(A).mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux.
!    L(A) = Linearized surface term (e.g. albedo)
!
              NC =  0
              FACTOR = LS_REFLEC
              BOA_ATTN = X0_BOA_4 * ATTN(NLAYERS,W_EXCIT)
              ESS_CUMSOURCE_UP(NC) = FACTOR * BOA_ATTN

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
                    ESS_CUMSOURCE_UP(NC) = &
                         UP_LOSTRANS(N) * ESS_CUMSOURCE_UP(NC-1)
                ENDDO

!  Set output

                SOURCE = ESS_CUMSOURCE_UP(NC)
                LS_ELASTIC_SS_UP(UTA,V,1) = SSFLUX * SOURCE
                IF  ( DO_FULL_RAMAN ) &
                      LS_RAMAN_SS_UP(UTA,V,1) = SSFLUX*SOURCE

!  Check for updating the recursion

                IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                NUT_PREV = NUT

!  end optical depth loop

              ENDDO

!  No surface clause

            ENDIF

!  end Upwelling clause (surface Jacobian)

         ENDIF

!  Finish geometry loop

        enddo
      enddo

!  Finish

      RETURN
      END SUBROUTINE LRRS_SSCORR_OG_MONO_PLUS_1

!

      SUBROUTINE OG_INTEGRATION_BIN_UP_PLUS_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXBINS, MAXVARS, &
          DO_RAMAN, NPOINTS_INNER, OFFSET_INNER, &
          NBINS, BINMAP, NLAYERS, NFINELAYERS, &
          DO_PROFILE_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
          DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
          EXTINCTION, L_EXTINCTION, RADII, ALPHA_ALL, ALPHA_FINE, &
          ATTN, ATTN_FINE, L_ATTN, L_ATTN_FINE, &
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS, &
          L_EMULTIPLIERS, L_IMULTIPLIERS, L_LOSTRANS )

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

      USE LRRS_PARS

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER, INTENT(IN) ::          maxlayers, maxfinelayers
      INTEGER, INTENT(IN) ::          maxpoints, maxbins, maxvars

!  control

      LOGICAL, INTENT(IN) ::          do_raman
      INTEGER, INTENT(IN) ::          nfinelayers, nlayers
      INTEGER, INTENT(IN) ::          npoints_inner, offset_inner
      INTEGER, INTENT(IN) ::          nbins  ( maxpoints )
      INTEGER, INTENT(IN) ::          binmap ( maxbins, maxpoints )

!  Profile linearization control inputs

      LOGICAL, INTENT(IN) ::          do_profile_linearization
      LOGICAL, INTENT(IN) ::          layer_vary_flag   (maxlayers)
      INTEGER, INTENT(IN) ::          layer_vary_number (maxlayers)

!  Column linearization control inputs

      LOGICAL, INTENT(IN) ::          do_column_linearization
      INTEGER, INTENT(IN) ::          n_totalcolumn_wfs

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

      REAL(FPK), INTENT(OUT) :: imultipliers   (maxlayers,   maxbins, maxpoints)
      REAL(FPK), INTENT(OUT) :: emultipliers   (maxlayers,   maxpoints)
      REAL(FPK), INTENT(OUT) :: lostrans      (maxlayers,   maxpoints)

!  Linearized quantities

      REAL(FPK), INTENT(OUT) :: l_imultipliers &
            ( maxvars, 0:maxlayers, maxlayers, maxbins, maxpoints)
      REAL(FPK), INTENT(OUT) :: l_emultipliers &
            ( maxvars, 0:maxlayers, maxlayers,   maxpoints)
      REAL(FPK), INTENT(OUT) :: l_lostrans (maxvars, maxlayers, maxpoints)

!  local arrays
!  ------------

!  Local geometry arrays

      REAL(FPK) :: csq_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: cot_fine ( maxlayers, maxfinelayers)
      REAL(FPK) :: tranj     ( maxfinelayers)

!  help variables
!  --------------

      INTEGER ::          n, j, q, k
      INTEGER ::          w, b, wb, wr

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
          lostrans(n,w)    = 0.0d0
          emultipliers(n,w) = 0.0d0
          do b = 1, nbins(w)
            imultipliers(n,b,w) = 0.0d0
          enddo
        enddo
      enddo

!  Whole layer linearizations

      if ( do_profile_linearization ) then
        do n = 1, nlayers
         do w = 1, npoints_inner
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_emultipliers(q,k,n,w) = 0.0d0
                do b = 1, nbins(w)
                  l_imultipliers(q,k,n,b,w) = 0.0d0
                enddo
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans(q,n,w) = 0.0d0
            enddo
          endif
         enddo
        enddo
      endif

      if ( do_column_linearization ) then
        do n = 1, nlayers
         do w = 1, npoints_inner
          do q = 1, n_totalcolumn_wfs
            L_emultipliers(q,0,n,w) = 0.0d0
            L_lostrans(q,n,w)        = 0.0d0
            do b = 1, nbins(w)
              l_imultipliers(q,0,n,b,w) = 0.0d0
            enddo
          enddo
         enddo
        enddo
      endif

!  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work up from the bottom of the atmosphere
!  =========================================

!  initialise

      n = nlayers
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

!  save some quantities

      do n = nlayers, 1, -1
        do j = 1, nfinelayers
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = 1.0d0 / salpha / salpha
        enddo
      enddo

!  Start layer loop
!  ================

      do n = nlayers, 1, -1

!  Save some quantities

        salpha = dsin(alpha_all(n-1))
        calpha = dcos(alpha_all(n-1))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  Start spectral loop

        do w = 1, npoints_inner

!  set up

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
          esum = 0.5d0 * ( func_1 + func_2 )
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
            rsum(b) = 0.5d0 * ( func_1 + func_2 )
            do j = 1, nfinelayers
              func = attn_fine(n,j,wb) * tranj(j) * csq_fine(n,j)
              rsum(b)  = rsum(b) + func
            enddo
            imultipliers(n,b,w) = rsum(b) * skn
          enddo

!  Profile linearizations
!  ----------------------

          if ( do_profile_linearization ) then

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
               l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
               l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
                 l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
                 l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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

!  End profile linearization clause

          endif

!  Column Linearizations
!  ---------------------

          if ( do_column_linearization ) then


!  Elastic

           do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(q,n,wr)
            l_func_2 = l_attn(q,0,n-1,wr) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(q,0,n,wr) *   tran_1 &
                                 + attn(n,wr)     * l_tran_1 )
            l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
            do q = 1, n_totalcolumn_wfs
             l_kn = raycon * l_extinction(q,n,wr)
             do b = 1, nbins(w)
              wb = binmap(b,w)
              l_func_2 = l_attn(q,0,n-1,wb) * csq_2
              l_tran_1 = -l_kn * argm_1 * tran_1
              l_func_1 = csq_1 * ( l_attn(q,0,n,wb) *   tran_1 &
                                   + attn(n,wb)      * l_tran_1 )
              l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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

!  end column linearization clause

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
      END SUBROUTINE OG_INTEGRATION_BIN_UP_PLUS_1

!

      SUBROUTINE OG_INTEGRATION_BIN_DN_PLUS_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXBINS, MAXVARS, &
          DO_RAMAN, NPOINTS_INNER, OFFSET_INNER, &
          NBINS, BINMAP, NLAYERS, NFINELAYERS, &
          DO_PROFILE_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
          DO_COLUMN_LINEARIZATION,N_TOTALCOLUMN_WFS, &
          EXTINCTION, L_EXTINCTION, RADII, ALPHA_ALL, ALPHA_FINE, &
          ATTN, ATTN_FINE, L_ATTN, L_ATTN_FINE, &
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS, &
          L_EMULTIPLIERS, L_IMULTIPLIERS, L_LOSTRANS )

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

      USE LRRS_PARS

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER, INTENT(IN) ::          maxlayers, maxfinelayers
      INTEGER, INTENT(IN) ::          maxpoints, maxbins, maxvars

!  control

      LOGICAL, INTENT(IN) ::          do_raman
      INTEGER, INTENT(IN) ::          nfinelayers, nlayers
      INTEGER, INTENT(IN) ::          npoints_inner, offset_inner
      INTEGER, INTENT(IN) ::          nbins  ( maxpoints )
      INTEGER, INTENT(IN) ::          binmap ( maxbins, maxpoints )

!  Profile linearization control inputs

      LOGICAL, INTENT(IN) ::          do_profile_linearization
      LOGICAL, INTENT(IN) ::          layer_vary_flag   (maxlayers)
      INTEGER, INTENT(IN) ::          layer_vary_number (maxlayers)

!  Column linearization control inputs

      LOGICAL, INTENT(IN) ::          do_column_linearization
      INTEGER, INTENT(IN) ::          n_totalcolumn_wfs

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

      INTEGER ::          n, j, q, k
      INTEGER ::          w, b, wb, wr

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
          lostrans(n,w)    = 0.0d0
          emultipliers(n,w) = 0.0d0
          do b = 1, nbins(w)
            imultipliers(n,b,w) = 0.0d0
          enddo
        enddo
      enddo

!  Whole layer linearizations

      if ( do_profile_linearization ) then
        do n = 1, nlayers
         do w = 1, npoints_inner
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_emultipliers(q,k,n,w) = 0.0d0
                do b = 1, nbins(w)
                  l_imultipliers(q,k,n,b,w) = 0.0d0
                enddo
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans(q,n,w) = 0.0d0
            enddo
          endif
         enddo
        enddo
      endif

      if ( do_column_linearization ) then
        do n = 1, nlayers
         do w = 1, npoints_inner
          do q = 1, n_totalcolumn_wfs
            L_emultipliers(q,0,n,w) = 0.0d0
            L_lostrans(q,n,w)        = 0.0d0
            do b = 1, nbins(w)
              l_imultipliers(q,0,n,b,w) = 0.0d0
            enddo
          enddo
         enddo
        enddo
      endif

!  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work down from the top of the atmosphere
!  ========================================

!  initialise

      n = 0
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

!  save some quantities

      do n = 1, nlayers
        do j = 1, nfinelayers
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = 1.0d0 / salpha / salpha
        enddo
      enddo

!  Start layer loop
!  ================

      do n = 1, nlayers

!  Save some quantities

        salpha = dsin(alpha_all(n))
        calpha = dcos(alpha_all(n))
        csq_2 = 1.0d0 / salpha / salpha
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
          esum = 0.5d0 * ( func_1 + func_2 )
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
            rsum(b) = 0.5d0 * ( func_1 + func_2 )
            do j = 1, nfinelayers
              func = attn_fine(n,j,wb) * tranj(j) * csq_fine(n,j)
              rsum(b) = rsum(b) + func
            enddo
            imultipliers(n,b,w) = rsum(b) * skn
          enddo

!  Profile linearizations
!  ----------------------

          if ( do_profile_linearization ) then

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
               l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
               l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
                 l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
                 l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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

!  End profile linearization clause

          endif

!  Column Linearizations
!  ---------------------

          if ( do_column_linearization ) then

!  Elastic

           do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(q,n,wr)
            l_func_2 = l_attn(q,0,n,wr) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(q,0,n-1,wr) *   tran_1 &
                                 + attn(n-1,wr)      * l_tran_1 )
            l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
            do q = 1, n_totalcolumn_wfs
             l_kn = raycon * l_extinction(q,n,wr)
             do b = 1, nbins(w)
              wb = binmap(b,w)
              l_func_2 = l_attn(q,0,n,wb) * csq_2
              l_tran_1 = -l_kn * argm_1 * tran_1
              l_func_1 = csq_1 * ( l_attn(q,0,n-1,wb) *   tran_1 &
                                   + attn(n-1,wb)      * l_tran_1 )
              l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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

!  end column linearization clause

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
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXVARS, &
          DO_RAMAN, NPOINTS_MONO, W_EXCIT, NLAYERS, NFINELAYERS, &
          DO_PROFILE_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
          DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
          EXTINCTION, L_EXTINCTION, RADII, ALPHA_ALL, ALPHA_FINE, &
          ATTN, ATTN_FINE, L_ATTN, L_ATTN_FINE, &
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS, &
          L_EMULTIPLIERS, L_IMULTIPLIERS, L_LOSTRANS )

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007 (not included here)

       USE LRRS_PARS

       IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER, INTENT(IN) ::          maxlayers, maxfinelayers
      INTEGER, INTENT(IN) ::          maxpoints, maxvars

!  control

      LOGICAL, INTENT(IN) ::          do_raman
      INTEGER, INTENT(IN) ::          nfinelayers, nlayers
      INTEGER, INTENT(IN) ::          npoints_mono, w_excit

!  Profile linearization control inputs

      LOGICAL, INTENT(IN) ::          do_profile_linearization
      LOGICAL, INTENT(IN) ::          layer_vary_flag   (maxlayers)
      INTEGER, INTENT(IN) ::          layer_vary_number (maxlayers)

!  Column linearization control inputs

      LOGICAL, INTENT(IN) ::          do_column_linearization
      INTEGER, INTENT(IN) ::          n_totalcolumn_wfs

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

      INTEGER ::          n, j, w, q, k

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
        lostrans(n)    = 0.0d0
        emultipliers(n) = 0.0d0
        do w = 1, npoints_mono
          imultipliers(n,w) = 0.0d0
        enddo
      enddo

!  Whole layer linearizations

      if ( do_profile_linearization ) then
        do n = 1, nlayers
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_emultipliers(q,k,n) = 0.0d0
                do w = 1, npoints_mono
                  l_imultipliers(q,k,n,w) = 0.0d0
                enddo
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans(q,n) = 0.0d0
            enddo
          endif
        enddo
      endif

      if ( do_column_linearization ) then
        do n = 1, nlayers
          do q = 1, n_totalcolumn_wfs
            L_emultipliers(q,0,n) = 0.0d0
            L_lostrans(q,n)        = 0.0d0
            do w = 1, npoints_mono
              l_imultipliers(q,0,n,w) = 0.0d0
            enddo
          enddo
        enddo
      endif

!  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work up from the bottom of the atmosphere
!  =========================================

!  initialise

      n = nlayers
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

!  save some quantities

      do n = nlayers, 1, -1
        do j = 1, nfinelayers
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = 1.0d0 / salpha / salpha
        enddo
      enddo

!  Start layer loop
!  ================

      do n = nlayers, 1, -1

!  Save some quantities

        salpha = dsin(alpha_all(n-1))
        calpha = dcos(alpha_all(n-1))
        csq_2 = 1.0d0 / salpha / salpha
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
        esum = 0.5d0 * ( func_1 + func_2 )
        do j = 1, nfinelayers
          func = attn_fine(n,j,w_excit) * tranj(j) * csq_fine(n,j)
          esum  = esum + func
        enddo
        emultipliers(n) = esum * skn

!  Inelastic scattering multipliers (by bins)

        do w = 1, npoints_mono
          func_2 = attn(n-1,w) *  csq_2
          func_1 = attn(n,w)   * csqt_1
          rsum(w) = 0.5d0 * ( func_1 + func_2 )
          do j = 1, nfinelayers
            func = attn_fine(n,j,w) * tranj(j) * csq_fine(n,j)
            rsum(w)  = rsum(w) + func
          enddo
          imultipliers(n,w) = rsum(w) * skn
        enddo

!  Profile linearizations
!  ----------------------

        if ( do_profile_linearization ) then

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
               l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
               l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
                 l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
                 l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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

!  End profile linearization clause

        endif

!  Column Linearizations
!  ---------------------

        if ( do_column_linearization ) then

!  Elastic

          do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(q,n,w_excit)
            l_func_2 = l_attn(q,0,n-1,w_excit) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(q,0,n,w_excit) *   tran_1 &
                                 + attn(n,w_excit)     * l_tran_1 )
            l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
            do q = 1, n_totalcolumn_wfs
             l_kn = raycon * l_extinction(q,n,w_excit)
             do w = 1, npoints_mono
              l_func_2 = l_attn(q,0,n-1,w) * csq_2
              l_tran_1 = -l_kn * argm_1 * tran_1
              l_func_1 = csq_1 * ( l_attn(q,0,n,w) *   tran_1 &
                                   + attn(n,w)      * l_tran_1 )
              l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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

!  end column linearization clause

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
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXVARS, &
          DO_RAMAN, NPOINTS_MONO, W_EXCIT, NLAYERS, NFINELAYERS, &
          DO_PROFILE_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
          DO_COLUMN_LINEARIZATION,N_TOTALCOLUMN_WFS, &
          EXTINCTION, L_EXTINCTION, RADII, ALPHA_ALL, ALPHA_FINE, &
          ATTN, ATTN_FINE, L_ATTN, L_ATTN_FINE, &
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS, &
          L_EMULTIPLIERS, L_IMULTIPLIERS, L_LOSTRANS )

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

       USE LRRS_PARS

       IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER, INTENT(IN) ::          maxlayers, maxfinelayers
      INTEGER, INTENT(IN) ::          maxpoints, maxvars

!  control

      LOGICAL, INTENT(IN) ::          do_raman
      INTEGER, INTENT(IN) ::          nfinelayers, nlayers
      INTEGER, INTENT(IN) ::          npoints_mono, w_excit

!  Profile linearization control inputs

      LOGICAL, INTENT(IN) ::          do_profile_linearization
      LOGICAL, INTENT(IN) ::          layer_vary_flag   (maxlayers)
      INTEGER, INTENT(IN) ::          layer_vary_number (maxlayers)

!  Column linearization control inputs

      LOGICAL, INTENT(IN) ::          do_column_linearization
      INTEGER, INTENT(IN) ::          n_totalcolumn_wfs

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

      INTEGER ::          n, j, q, k, w

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
        lostrans(n)    = 0.0d0
        emultipliers(n) = 0.0d0
        do w = 1, npoints_mono
          imultipliers(n,w) = 0.0d0
        enddo
      enddo

!  Whole layer linearizations

      if ( do_profile_linearization ) then
        do n = 1, nlayers
          do k = 1, nlayers
            if ( layer_vary_flag(k) ) then
              do q = 1, layer_vary_number(k)
                L_emultipliers(q,k,n) = 0.0d0
                do w = 1, npoints_mono
                  l_imultipliers(q,k,n,w) = 0.0d0
                enddo
              enddo
            endif
          enddo
          if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
              L_lostrans(q,n) = 0.0d0
            enddo
          endif
        enddo
      endif

      if ( do_column_linearization ) then
        do n = 1, nlayers
          do q = 1, n_totalcolumn_wfs
            L_emultipliers(q,0,n) = 0.0d0
            L_lostrans(q,n)        = 0.0d0
            do w = 1, npoints_mono
              l_imultipliers(q,0,n,w) = 0.0d0
            enddo
          enddo
        enddo
      endif

!  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work down from the top of the atmosphere
!  ========================================

!  initialise

      n = 0
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

!  save some quantities

      do n = 1, nlayers
        do j = 1, nfinelayers
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = 1.0d0 / salpha / salpha
        enddo
      enddo

!  Start layer loop
!  ================

      do n = 1, nlayers

!  Save some quantities

        salpha = dsin(alpha_all(n))
        calpha = dcos(alpha_all(n))
        csq_2 = 1.0d0 / salpha / salpha
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
        esum = 0.5d0 * ( func_1 + func_2 )
        do j = 1, nfinelayers
          func = attn_fine(n,j,w_excit) * tranj(j) * csq_fine(n,j)
          esum  = esum + func
        enddo
        emultipliers(n) = esum * skn

!  Inelastic scattering multipliers (by bins)

        do w = 1, npoints_mono
          func_2 = attn(n,w)   *  csq_2
          func_1 = attn(n-1,w) * csqt_1
          rsum(w) = 0.5d0 * ( func_1 + func_2 )
          do j = 1, nfinelayers
            func = attn_fine(n,j,w) * tranj(j) * csq_fine(n,j)
            rsum(w) = rsum(w) + func
          enddo
          imultipliers(n,w) = rsum(w) * skn
        enddo

!  Profile linearizations
!  ----------------------

        if ( do_profile_linearization ) then

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
               l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
               l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
                 l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
                 l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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

!  End profile linearization clause

        endif

!  Column Linearizations
!  ---------------------

        if ( do_column_linearization ) then

!  Elastic

          do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(q,n,w_excit)
            l_func_2 = l_attn(q,0,n,w_excit) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(q,0,n-1,w_excit) *   tran_1 &
                                 + attn(n-1,w_excit)      * l_tran_1 )
            l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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
            do q = 1, n_totalcolumn_wfs
             l_kn = raycon * l_extinction(q,n,w_excit)
             do w = 1, npoints_mono
              l_func_2 = l_attn(q,0,n,w) * csq_2
              l_tran_1 = -l_kn * argm_1 * tran_1
              l_func_1 = csq_1 * ( l_attn(q,0,n-1,w) *   tran_1 &
                                   + attn(n-1,w)      * l_tran_1 )
              l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
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

!  end column linearization clause

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
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXVARS, &
          N_RAMAN_POINTS, NLAYERS, NFINELAYERS, &
          DO_PROFILE_LINEARIZATION, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
          DO_COLUMN_LINEARIZATION, N_TOTALCOLUMN_WFS, &
          EXTINCTION, L_EXTINCTION, &
          SUNPATHS, NTRAVERSE, SUNPATHS_FINE, NTRAVERSE_FINE, &
          ATTN, ATTN_FINE, L_ATTN, L_ATTN_FINE )

!  Does attenuations

      USE LRRS_PARS

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER, INTENT(IN) ::          maxlayers, maxfinelayers
      INTEGER, INTENT(IN) ::          maxpoints, maxvars

!  control

      INTEGER, INTENT(IN) ::          nfinelayers, nlayers, n_raman_points

!  Profile linearization control inputs

      LOGICAL, INTENT(IN) ::          do_profile_linearization
      LOGICAL, INTENT(IN) ::          layer_vary_flag   (maxlayers)
      INTEGER, INTENT(IN) ::          layer_vary_number (maxlayers)

!  Column linearization control inputs

      LOGICAL, INTENT(IN) ::          do_column_linearization
      INTEGER, INTENT(IN) ::          n_totalcolumn_wfs

!  Whole layers

      INTEGER, INTENT(IN) ::          ntraverse(0:maxlayers)
      REAL(FPK), INTENT(IN) :: sunpaths(0:maxlayers,maxlayers)

!  Fine level

      INTEGER, INTENT(IN) ::          ntraverse_fine(maxlayers,maxfinelayers)
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

      INTEGER ::          w, n, j, k, q
      REAL(FPK) :: tau, l_iop, sum

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0D0

!  attenuation functions, whole layers

      do w = 1, n_raman_points
       do n = 0, nlayers
        tau     = 0.0d0
        attn(n,w) = 0.0d0
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
         tau            = 0.0d0
         attn_fine(n,j,w) = 0.0d0
         do k = 1, ntraverse_fine(n,j)
          tau = tau + sunpaths_fine(n,k,j) * extinction(k,w)
         enddo
         if ( tau .le. local_cutoff ) attn_fine(n,j,w) = dexp(-tau)
        enddo
       enddo
      enddo

!  Linearized attenuations
!  -----------------------

!  Linearized attenuation factors. Profile linearization.
!     Note the extra zeroing...............bug, 01 November 2007

      if ( do_profile_linearization ) then
       do w = 1, n_raman_points
        do n = 0, nlayers
          do k = n+1, nlayers
            do q = 1, layer_vary_number(k)
              l_attn(q,k,n,w) = 0.0d0
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
       enddo
      endif

!  Linearized attenuation factors. Column linearization

      if ( do_column_linearization ) then
       do w = 1, n_raman_points
        do n = 0, nlayers
          do q = 1, n_totalcolumn_wfs
            sum = 0.0d0
            do k = 1, ntraverse(n)
              l_iop = - l_extinction(q,k,w)
              sum = sum + sunpaths(n,k) * l_iop
            enddo
            l_attn(q,0,n,w) = sum * attn(n,w)
          enddo
        enddo
        do n = 1, nlayers
         do j = 1, nfinelayers
          do q = 1, n_totalcolumn_wfs
            sum = 0.0d0
            do k = 1, ntraverse_fine(n,j)
              l_iop = - l_extinction(q,k,w)
              sum = sum + sunpaths_fine(n,k,j) * l_iop
            enddo
            l_attn_fine(q,0,n,j,w) = sum * attn_fine(n,j,w)
          enddo
         enddo
        enddo
       enddo
      endif

!  Finish

      return
      END SUBROUTINE OG_ATTENUATIONS_PLUS_1

      END MODULE lrrs_corrections_plus_1
