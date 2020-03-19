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
! # Subroutines in this Module (Suffix '_1' no partials)        #
! #                                                             #
! #             LRRS_SSCORR_OUTGOING_BIN_1  (master)            #
! #             LRRS_SSCORR_OUTGOING_MONO_1 (master)            #
! #                                                             #
! #                 outgoing_sphergeom_fine_up_1                #
! #                 outgoing_sphergeom_fine_dn_1                #
! #                                                             #
! #                 outgoing_attenuations_1                     #
! #                                                             #
! #                 outgoing_integration_bin_up_1               #
! #                 outgoing_integration_bin_dn_1               #
! #                                                             #
! #                 outgoing_integration_mono_up_1              #
! #                 outgoing_integration_mono_dn_1              #
! #                                                             #
! #                 multi_outgoing_adjustgeom                   #
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


      MODULE lrrs_corrections_1

      USE LRRS_PARS

      PRIVATE
      PUBLIC :: LRRS_SSCORR_OUTGOING_BIN_1,&
                LRRS_SSCORR_OUTGOING_MONO_1,&
                OUTGOING_SPHERGEOM_FINE_UP_1,&
                OUTGOING_SPHERGEOM_FINE_DN_1,&
                MULTI_OUTGOING_ADJUSTGEOM

      CONTAINS

      SUBROUTINE LRRS_SSCORR_OUTGOING_BIN_1 &
         ( DO_FULL_RAMAN, DO_UPWELLING, DO_DNWELLING, &
           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_WATERLEAVING, &
           DO_LAMBERTIAN_SURFACE, DO_LAMBERTIAN_FRACTION, DO_DMSCALING, &
           NLAYERS, NFINELAYERS, NMOMENTS_INPUT, &
           NPOINTS_INNER, OFFSET_INNER, NPOINTS_OUTER, &
           N_LOUTPUT, LEVELMASK_UP, LEVELMASK_DN, &
           N_USER_STREAMS, N_USER_RELAZMS, &
           SOLAR_ANGLES_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST, &
           LAMBERTIAN_FRACTION, EXACTDB_BRDFUNC, ALBEDOS_RANKED, &
           FLUXES_RANKED, EARTH_RADIUS, HEIGHT_GRID, &
           TRUNC_FACTORS, DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC_UNSCALED, &
           N_RRSBINS, BINMAP, &
           OMEGAMOMS_CABANNES_UNSCALED, &
           OMEGAMOMS_RRSLOSS_UNSCALED, &
           OMEGAMOMS_RRSBIN_UNSCALED, &
           ELASTIC_SS_UP, ELASTIC_SS_DN, &
           RAMAN_SS_UP,   RAMAN_SS_DN, &
           FAIL, MESSAGE)

!  Single scatter exact calculation for the outgoing LOS

!   Programmed by R. Spurr, RT Solutions Inc.
!    First Draft,  for LRRS Version 2.1, April 23rd 2008
!    Second draft, for LRRS Version 2.2, 23 July 2009

!  Include files
!  -------------

!  include file of dimensions and numbers

      USE LRRS_PARS

      IMPLICIT NONE

!  Input
!  -----

!  Flags

      LOGICAL, INTENT(IN) ::   DO_CABANNES_RAMAN
      LOGICAL, INTENT(IN) ::   DO_ENERGY_BALANCING
      LOGICAL, INTENT(IN) ::   DO_UPWELLING
      LOGICAL, INTENT(IN) ::   DO_DNWELLING
      LOGICAL, INTENT(IN) ::   DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT(IN) ::   DO_WATERLEAVING
      LOGICAL, INTENT(IN) ::   DO_LAMBERTIAN_FRACTION
      LOGICAL, INTENT(IN) ::   DO_DMSCALING
      LOGICAL, INTENT(IN) ::   DO_FULL_RAMAN

!  number of layers

      INTEGER, INTENT(IN) ::   NLAYERS

!  Number of input phase function Legendre moments
!    THIS IS A BASIC INPUT THAT USER MUST PROVIDE

      INTEGER, INTENT(IN) ::   NMOMENTS_INPUT

!  Number of fine layers for sscorr outgoing

      INTEGER, INTENT(IN) ::   NFINELAYERS

!  Solar beam input geometry

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLES_ADJUST &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  User stream variables
!   @@@ RobFix 5/5/11, Wrong dimensioning. Fixed 5/5/11

      INTEGER, INTENT(IN) ::   N_USER_STREAMS
!      REAL(FPK), INTENT(IN) :: USER_ANGLES_ADJUST ( MAX_USER_RELAZMS )
      REAL(FPK), INTENT(IN) :: USER_ANGLES_ADJUST ( MAX_USER_STREAMS )

!  User azimuth variables

      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      REAL(FPK), INTENT(IN) :: USER_RELAZMS_ADJUST &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)

      INTEGER, INTENT(IN) ::   N_LOUTPUT

!  outer/inner wavelength range, and offset for the inner range
!    Depends on the choice of solar spectrum

      INTEGER, INTENT(IN) ::   OFFSET_INNER
      INTEGER, INTENT(IN) ::   NPOINTS_INNER
      INTEGER, INTENT(IN) ::   NPOINTS_OUTER

!  Layer masks for doing integrated source terms

      INTEGER, INTENT(IN) ::   LEVELMASK_UP (MAX_LOUTPUT)
      INTEGER, INTENT(IN) ::   LEVELMASK_DN (MAX_LOUTPUT)

!  Earth Radius

      REAL(FPK), INTENT(IN) :: EARTH_RADIUS

!  Height grid

      REAL(FPK), INTENT(IN) :: HEIGHT_GRID ( 0:MAX_LAYERS )

!  These must be defined on the outer wavelength grid (binning)
!  Unscaled quantities, Elastic input

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: TRUNC_FACTORS     ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_CABANNES_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSLOSS_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSBIN_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Bin Mapping quantities (internal derivation)

      INTEGER, INTENT(IN) ::   N_RRSBINS  ( MAX_POINTS )
      INTEGER, INTENT(IN) ::   BINMAP ( MAX_BINS, MAX_POINTS )

!  Exact DB values

      REAL(FPK), INTENT(IN) :: EXACTDB_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Lambertian fraction

      REAL(FPK), INTENT(IN) :: LAMBERTIAN_FRACTION

!  Albedos and solar fluxes

      REAL(FPK), INTENT(IN) :: FLUXES_RANKED   ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: ALBEDOS_RANKED  ( MAX_POINTS )

!  Output
!  ------

!  single scatter results

      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_UP &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_DN &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_UP &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_DN &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Exception handling

      LOGICAL, INTENT(OUT) ::           fail
      CHARACTER (LEN=*), INTENT(OUT) :: message

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

      REAL(FPK) :: extinction ( max_layers, max_points )

!  Saved Legendre polynomials

      REAL(FPK) :: &
         SS_PLEG_UP(MAX_LAYERS,0:MAX_MOMENTS_INPUT), &
         SS_PLEG_DN(MAX_LAYERS,0:MAX_MOMENTS_INPUT)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(FPK) :: TMS(MAX_LAYERS,MAX_POINTS)

!  Local truncation factors for additional DELTAM scaling

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

      REAL(FPK) :: &
         UP_IMULTIPLIERS(MAX_LAYERS,MAX_BINS,MAX_POINTS), &
         DN_IMULTIPLIERS(MAX_LAYERS,MAX_BINS,MAX_POINTS), &
         UP_EMULTIPLIERS(MAX_LAYERS,MAX_POINTS), &
         DN_EMULTIPLIERS(MAX_LAYERS,MAX_POINTS), &
         UP_LOSTRANS(MAX_LAYERS,MAX_POINTS), &
         DN_LOSTRANS(MAX_LAYERS,MAX_POINTS)

!  Solar beam attenuations
!         to BOA (required for exact DB calculation)

      REAL(FPK) :: attn      ( 0:max_layers, max_points )
      REAL(FPK) :: attn_fine &
               ( max_layers, max_fine_layers, max_points )

!  local variables
!  ---------------

!  Indices

      INTEGER ::       N, NUT, NSTART, NUT_PREV, NLEVEL, L, W, B
      INTEGER ::       UTA, UM, IA, NC, V, WR, WB

!  help variables (REAL(FPK) ::)

      REAL(FPK) :: boa_attn, ssflux, x0_boa, ratio
      REAL(FPK) :: SOURCE, HELP, COSSCAT, CONTRIB
      REAL(FPK) :: PHAS_E, PHAS_LR, PHAS_LG(MAX_BINS), PHAS_CB
      REAL(FPK) :: FACTOR, BETA, CAB, GAIN, LOSS
      REAL(FPK) :: DF1(0:MAX_MOMENTS_INPUT)
      REAL(FPK) :: DF2(0:MAX_MOMENTS_INPUT), OMEGA_LOCAL
      REAL(FPK) :: REFLEC(MAX_POINTS), FL1, FL2, BRDF

!  Local surface flag

      LOGICAL, PARAMETER :: DO_INCLUDE_SURFACE = .true.

!  Set up operations
!  -----------------

!  Ss flux scaling factor

      SSFLUX = one / PI4

!  Fine layering flag

      do_fine = .true.

!   Number of fine layers now a control input (7/25/08)

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
            OMEGA_LOCAL = OMEGAMOMS_ELASTIC_UNSCALED(N,0,WR )
            HELP = TRUNC_FACTORS(N,WR) * OMEGA_LOCAL
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

!  Start the main loop over all viewing geometries
!   Use the adjusted values of the angles -------> March 2011

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

            call outgoing_attenuations_1 &
        ( max_layers, max_fine_layers, max_points, &
          npoints_outer, nlayers, nfinelayers, &
          extinction, sunpaths, ntraverse, &
          sunpaths_fine, ntraverse_fine, &
          attn, attn_fine )

!         do w = 1, npoints_outer
!           write(55,*)w,attn(23,w),attn(1,w)
!         enddo
!         pause

!  Multipliers, transmittances

            call outgoing_integration_bin_up_1 &
        ( max_layers, max_fine_layers, max_points, max_bins, &
          npoints_inner, offset_inner, n_rrsbins, binmap, &
          nlayers, nfinelayers, extinction, &
          radii, alpha_all, alpha_fine, attn, attn_fine, &
          up_emultipliers,  up_imultipliers, up_lostrans )

!  Debug
!              w = 1
!              do n = 1, nlayers
!                write(28,'(2i5,1p2e18.10)')n,w,
!     &           up_lostrans(n,w),up_emultipliers(n,w)
!               do b = 1, n_rrsbins(w)
!                write(44,'(2i5,1pe18.10)')
!     &    w,binmap(b,w),up_imultipliers(n,b,w)
!               enddo
!              enddo
!            pause

!             write(*,*)npoints_inner, npoints_outer, offset_inner
!              n = 3
!              do w = 1, npoints_inner
!                write(45,'(2i5,1p2e18.10)')n,w,
!     &           up_lostrans(n,w),up_emultipliers(n,w)
!               do b = 1, n_rrsbins(w)
!                write(45,'(2i5,1pe18.10)')
!     &    w,binmap(b,w),up_imultipliers(n,b,w)
!               enddo
!              enddo
!             enddo
!            pause

!  legendre polynomials

            COSSCAT = COSSCAT_UP(NLAYERS)
            SS_PLEG_UP(1,0) = ONE
            SS_PLEG_UP(1,1) = COSSCAT
            DO L = 2, NMOMENTS_INPUT
              SS_PLEG_UP(1,L) = &
               DF1(L) * SS_PLEG_UP(1,L-1) * COSSCAT  - &
               DF2(L) * SS_PLEG_UP(1,L-2)
            ENDDO

!  Elastic scatting Phase functions (multiplied by TMS factor)
!  ===========================================================

!  Layer/point loops, whole layer contributions

            DO N = 1, NLAYERS
              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
                PHAS_E = ZERO
                DO L = 0, NMOMENTS_INPUT
                  BETA   = OMEGAMOMS_ELASTIC_UNSCALED(N,L,WR)
                  PHAS_E = PHAS_E + BETA * SS_PLEG_UP(1,L)
                ENDDO
                PHAS_E = PHAS_E * TMS(N,W)
                ESCAT_UP(N,W) = PHAS_E * UP_EMULTIPLIERS(N,W)
              ENDDO
            ENDDO

!  Add inelastic scatting Phase functions if flagged
!  =================================================

            IF ( DO_FULL_RAMAN ) THEN

!  Energy balancing
!  ----------------

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

                  PHAS_CB = ZERO
                  DO L = 0, NMOMENTS_INPUT
                    BETA    = OMEGAMOMS_CABANNES_UNSCALED(N,L,W)
                    PHAS_CB = PHAS_CB + BETA * SS_PLEG_UP(1,L)
                  ENDDO
                  PHAS_CB = PHAS_CB * TMS(N,W)
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

            call outgoing_attenuations_1 &
        ( max_layers, max_fine_layers, max_points, &
          npoints_outer, nlayers, nfinelayers, &
          extinction, sunpaths, ntraverse, &
          sunpaths_fine, ntraverse_fine, &
          attn, attn_fine )

!         do w = 1, npoints_outer
!           write(55,*)w,attn(23,w),attn(1,w)
!         enddo
!         pause

!  Multipliers and transmittances

            call outgoing_integration_bin_dn_1 &
        ( max_layers, max_fine_layers, max_points, max_bins, &
          npoints_inner, offset_inner, n_rrsbins, binmap, &
          nlayers, nfinelayers, extinction, &
          radii, alpha_all, alpha_fine, attn, attn_fine, &
          dn_emultipliers, dn_imultipliers, dn_lostrans )

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
            SS_PLEG_DN(1,0) = ONE
            SS_PLEG_DN(1,1) = COSSCAT
            DO L = 2, NMOMENTS_INPUT
            SS_PLEG_DN(1,L) = &
               DF1(L) * SS_PLEG_DN(1,L-1) * COSSCAT  - &
               DF2(L) * SS_PLEG_DN(1,L-2)
            ENDDO

!  Elastic scatting Phase functions (multiplied by TMS factor)
!  ===========================================================

!  layer/point loops, whole layer contribuions

            DO N = 1, NLAYERS
              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
                PHAS_E = ZERO
                DO L = 0, NMOMENTS_INPUT
                  BETA   = OMEGAMOMS_ELASTIC_UNSCALED(N,L,WR)
                  PHAS_E = PHAS_E + BETA * SS_PLEG_DN(1,L)
                ENDDO
                PHAS_E = PHAS_E * TMS(N,W)
                ESCAT_DN(N,W) = PHAS_E * DN_EMULTIPLIERS(N,W)
              ENDDO
            ENDDO

!  Add inelastic scatting Phase functions if flagged
!  =================================================

            IF ( DO_FULL_RAMAN ) THEN

!  Energy balancing
!  ----------------

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

                  PHAS_CB = ZERO
                  DO L = 0, NMOMENTS_INPUT
                    BETA    = OMEGAMOMS_CABANNES_UNSCALED(N,L,W)
                    PHAS_CB = PHAS_CB + BETA * SS_PLEG_DN(1,L)
                  ENDDO
                  PHAS_CB = PHAS_CB * TMS(N,W)
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

!  End Downwelling clause

          ENDIF

!  Recurrence relation for the UPWELLING intensity
!  ===============================================

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
!              BRDF = BRDF_FACTOR * EXACTDB_BRDFUNC(UM,IA)
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
              X0_BOA = DCOS(SOLAR_ANGLES_ADJUST(UM,IA) * DEG_TO_RAD)
              DO W = 1, NPOINTS_INNER
                WR = W + OFFSET_INNER
!                FACTOR = ALBEDOS_RANKED(WR)
                FACTOR = REFLEC(W)
                BOA_ATTN = 4.0d0 * X0_BOA * ATTN(NLAYERS,WR)
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

!  end optical depth loop and Downwelling clause

            ENDDO
          ENDIF

!  Finish geometry loop

        enddo
      enddo

!  Finish

      RETURN
      END SUBROUTINE LRRS_SSCORR_OUTGOING_BIN_1

!

      SUBROUTINE LRRS_SSCORR_OUTGOING_MONO_1 &
         ( DO_FULL_RAMAN, DO_UPWELLING, DO_DNWELLING, &
           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_WATERLEAVING, &
           DO_LAMBERTIAN_SURFACE, DO_LAMBERTIAN_FRACTION, DO_DMSCALING, &
           NLAYERS, NFINELAYERS, NMOMENTS_INPUT, &
           NPOINTS_MONO, W_EXCIT, &
           N_LOUTPUT, LEVELMASK_UP, LEVELMASK_DN, &
           N_USER_STREAMS, N_USER_RELAZMS, &
           SOLAR_ANGLES_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST, &
           LAMBERTIAN_FRACTION, EXACTDB_BRDFUNC, ALBEDOS_RANKED, &
           FLUXES_RANKED, EARTH_RADIUS, HEIGHT_GRID, &
           TRUNC_FACTORS, DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC_UNSCALED, &
           OMEGAMOMS_CABANNES_UNSCALED, &
           OMEGAMOMS_RRSLOSS_UNSCALED, &
           OMEGAMOMS_RRSGAIN_UNSCALED, &
           ELASTIC_SS_UP, ELASTIC_SS_DN, &
           RAMAN_SS_UP,   RAMAN_SS_DN, &
           FAIL, MESSAGE)

!  Single scatter exact calculation for the outgoing LOS

!   Programmed by R. Spurr, RT Solutions Inc.
!    First Draft,  for LRRS Version 2.1, April 23rd 2008
!    Second draft, for LRRS Version 2.2, 23 July 2009

!  Include files
!  -------------

!  include file of dimensions and numbers

      USE LRRS_PARS

      IMPLICIT NONE

!  Input
!  -----

!  Flags

      LOGICAL, INTENT(IN) ::   DO_CABANNES_RAMAN
      LOGICAL, INTENT(IN) ::   DO_ENERGY_BALANCING
      LOGICAL, INTENT(IN) ::   DO_UPWELLING
      LOGICAL, INTENT(IN) ::   DO_DNWELLING
      LOGICAL, INTENT(IN) ::   DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT(IN) ::   DO_WATERLEAVING
      LOGICAL, INTENT(IN) ::   DO_LAMBERTIAN_FRACTION
      LOGICAL, INTENT(IN) ::   DO_DMSCALING
      LOGICAL, INTENT(IN) ::   DO_FULL_RAMAN

!  number of layers

      INTEGER, INTENT(IN) ::   NLAYERS

!  Number of input phase function Legendre moments
!    THIS IS A BASIC INPUT THAT USER MUST PROVIDE

      INTEGER, INTENT(IN) ::   NMOMENTS_INPUT

!  Number of fine layers for sscorr outgoing

      INTEGER, INTENT(IN) ::   NFINELAYERS

!  Solar beam input geometry

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLES_ADJUST &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  User stream variables

      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_ANGLES_ADJUST ( MAX_USER_STREAMS )

!  User azimuth variables

      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      REAL(FPK), INTENT(IN) :: USER_RELAZMS_ADJUST &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)

      INTEGER, INTENT(IN) ::   N_LOUTPUT

!  Number of points, excitation index

      INTEGER, INTENT(IN) ::   NPOINTS_MONO, W_EXCIT

!  Layer masks for doing integrated source terms

      INTEGER, INTENT(IN) ::   LEVELMASK_UP (MAX_LOUTPUT)
      INTEGER, INTENT(IN) ::   LEVELMASK_DN (MAX_LOUTPUT)

!  Earth Radius

      REAL(FPK), INTENT(IN) :: EARTH_RADIUS

!  Height grid

      REAL(FPK), INTENT(IN) :: HEIGHT_GRID ( 0:MAX_LAYERS )

!  These must be defined on the outer wavelength grid (binning)
!  Unscaled quantities, Elastic input

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: TRUNC_FACTORS     ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) :: OMEGAMOMS_ELASTIC_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_CABANNES_UNSCALED &
           ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSLOSS_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: OMEGAMOMS_RRSGAIN_UNSCALED &
               ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Exact DB values

      REAL(FPK), INTENT(IN) :: EXACTDB_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Lambertian fraction

      REAL(FPK), INTENT(IN) :: LAMBERTIAN_FRACTION

!  Albedos and solar fluxes

      REAL(FPK), INTENT(IN) :: FLUXES_RANKED   ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: ALBEDOS_RANKED  ( MAX_POINTS )

!  Output
!  ------

!  single scatter results

      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_UP &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: ELASTIC_SS_DN &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_UP &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: RAMAN_SS_DN &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Exception handling

      LOGICAL, INTENT(OUT) ::           fail
      CHARACTER (LEN=*), INTENT(OUT) :: message

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

      REAL(FPK) :: extinction ( max_layers, max_points )

!  Saved Legendre polynomials

      REAL(FPK) :: &
         SS_PLEG_UP(MAX_LAYERS,0:MAX_MOMENTS_INPUT), &
         SS_PLEG_DN(MAX_LAYERS,0:MAX_MOMENTS_INPUT)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(FPK) :: TMS(MAX_LAYERS)

!  Local truncation factors for additional DELTAM scaling

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

      REAL(FPK) :: &
         UP_EMULTIPLIERS(MAX_LAYERS), &
         DN_EMULTIPLIERS(MAX_LAYERS), &
         UP_IMULTIPLIERS(MAX_LAYERS,MAX_POINTS), &
         DN_IMULTIPLIERS(MAX_LAYERS,MAX_POINTS), &
         UP_LOSTRANS(MAX_LAYERS), &
         DN_LOSTRANS(MAX_LAYERS)

!  Solar beam attenuations
!         to BOA (required for exact DB calculation)

      REAL(FPK) :: attn      ( 0:max_layers, max_points )
      REAL(FPK) :: attn_fine &
               ( max_layers, max_fine_layers, max_points )

!  local variables
!  ---------------

!  Indices

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, L, W
      INTEGER ::          UTA, UM, IA, NC, V, WR

!  help variables (REAL(FPK) ::)

      REAL(FPK) :: x0_boa, boa_attn, omega_local, ssflux, ratio
      REAL(FPK) :: SOURCE, HELP, LSOURCE, COSSCAT, CONTRIB
      REAL(FPK) :: PHAS_E, PHAS_LR, PHAS_LG(MAX_POINTS), PHAS_CB
      REAL(FPK) :: FACTOR, BETA, CAB, GAIN, LOSS
      REAL(FPK) :: DF1(0:MAX_MOMENTS_INPUT)
      REAL(FPK) :: DF2(0:MAX_MOMENTS_INPUT)
      REAL(FPK) :: REFLEC, FL1, FL2, BRDF

!  Local surface flag

      LOGICAL, PARAMETER :: DO_INCLUDE_SURFACE = .true.

!  Set up operations
!  -----------------

!  set flux scaling

      SSFLUX = one / PI4

!  Fine layering flag

      do_fine = .true.

!   Number of fine layers now a control input (7/25/08)

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
          OMEGA_LOCAL = OMEGAMOMS_ELASTIC_UNSCALED (N, 0, W_EXCIT )
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

!  Start the main loop over all viewing geometries
!   Use the adjusted values of the angles -------> March 2011

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

            call outgoing_attenuations_1 &
        ( max_layers, max_fine_layers, max_points, &
          npoints_mono, nlayers, nfinelayers, &
          extinction, sunpaths, ntraverse, &
          sunpaths_fine, ntraverse_fine, &
          attn, attn_fine )

!  Multipliers, transmittances

            call outgoing_integration_mono_up_1 &
        ( max_layers, max_fine_layers, &
          max_points, npoints_mono, w_excit, &
          nlayers, nfinelayers, extinction, &
          radii, alpha_all, alpha_fine, attn, attn_fine, &
          up_emultipliers, up_imultipliers, up_lostrans )

!  Debug
!              do n = 1, nlayers
!                write(45,'(1p2e18.10)')
!     &           up_lostrans(n),up_emultipliers(n)
!              enddo
!              pause

!  legendre polynomials

            COSSCAT = COSSCAT_UP(NLAYERS)
            SS_PLEG_UP(1,0) = ONE
            SS_PLEG_UP(1,1) = COSSCAT
            DO L = 2, NMOMENTS_INPUT
              SS_PLEG_UP(1,L) = &
               DF1(L) * SS_PLEG_UP(1,L-1) * COSSCAT  - &
               DF2(L) * SS_PLEG_UP(1,L-2)
            ENDDO

!  Elastic scatting Phase functions (multiplied by TMS factor)
!  ===========================================================

!  Start loop

            DO N = 1, NLAYERS
              PHAS_E = ZERO
              DO L = 0, NMOMENTS_INPUT
                BETA   = OMEGAMOMS_ELASTIC_UNSCALED(N,L,W_EXCIT)
                PHAS_E = PHAS_E + BETA * SS_PLEG_UP(1,L)
              ENDDO
              PHAS_E = PHAS_E * TMS(N)
              ESCAT_UP(N) = PHAS_E * UP_EMULTIPLIERS(N)
            ENDDO

!  Add inelastic scatting Phase functions if flagged
!  -------------------------------------------------

            IF ( DO_FULL_RAMAN ) THEN

!  Energy Balancing
!  ----------------

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

                PHAS_CB = ZERO
                DO L = 0, NMOMENTS_INPUT
                  BETA    = OMEGAMOMS_CABANNES_UNSCALED(N,L,1)
                  PHAS_CB = PHAS_CB + BETA * SS_PLEG_UP(1,L)
                ENDDO
                CAB = PHAS_CB * UP_EMULTIPLIERS(N)
! @@@ Robfix 06/01/11. Add TMS
                CAB = CAB * TMS(N)

!  Gain term

                GAIN = ZERO
                DO W = 1, NPOINTS_MONO
                  RATIO   = FLUXES_RANKED(W)/FLUXES_RANKED(W_EXCIT)
                  CONTRIB = OMEGAMOMS_RRSGAIN_UNSCALED(N,0,W) + &
                    OMEGAMOMS_RRSGAIN_UNSCALED(N,2,W) * SS_PLEG_UP(1,2)
                  PHAS_LG(W)  = CONTRIB * RATIO
                  GAIN = GAIN + PHAS_LG(W) * UP_IMULTIPLIERS(N,W)
                ENDDO
! @@@ Robfix 06/01/11. Add TMS
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

            call outgoing_attenuations_1 &
        ( max_layers, max_fine_layers, max_points, &
          npoints_mono, nlayers, nfinelayers, &
          extinction, sunpaths, ntraverse, &
          sunpaths_fine, ntraverse_fine, &
          attn, attn_fine )

!  Multipliers and transmittances

            call outgoing_integration_mono_dn_1 &
        ( max_layers, max_fine_layers, &
          max_points, npoints_mono, w_excit, &
          nlayers, nfinelayers, extinction, &
          radii, alpha_all, alpha_fine, attn, attn_fine, &
          dn_emultipliers, dn_imultipliers, dn_lostrans )

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
            SS_PLEG_DN(1,0) = ONE
            SS_PLEG_DN(1,1) = COSSCAT
            DO L = 2, NMOMENTS_INPUT
            SS_PLEG_DN(1,L) = &
               DF1(L) * SS_PLEG_DN(1,L-1) * COSSCAT  - &
               DF2(L) * SS_PLEG_DN(1,L-2)
            ENDDO

!  Elastic scatting Phase functions (multiplied by TMS factor)
!  ===========================================================

!  Start loop

            DO N = 1, NLAYERS
              PHAS_E = ZERO
              DO L = 0, NMOMENTS_INPUT
                BETA   = OMEGAMOMS_ELASTIC_UNSCALED(N,L,W_EXCIT)
                PHAS_E = PHAS_E + BETA * SS_PLEG_DN(1,L)
              ENDDO
              PHAS_E = PHAS_E * TMS(N)
              ESCAT_DN(N) = PHAS_E * DN_EMULTIPLIERS(N)
            ENDDO

!  Add inelastic scatting Phase functions if flagged
!  =================================================

            IF ( DO_FULL_RAMAN ) THEN

!  Energy balancing
!  ----------------

             IF ( DO_ENERGY_BALANCING ) THEN

!  Start loop

              DO N = 1, NLAYERS

!  Loss term

                PHAS_LR = OMEGAMOMS_RRSLOSS_UNSCALED(N,0,1) + &
                    OMEGAMOMS_RRSLOSS_UNSCALED(N,2,1) * SS_PLEG_DN(1,2)
                LOSS = PHAS_LR * DN_EMULTIPLIERS(N)
! @@@ Robfix 06/01/11. Add TMS
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
! @@@ Robfix 06/01/11. Add TMS
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

                PHAS_CB = ZERO
                DO L = 0, NMOMENTS_INPUT
                  BETA    = OMEGAMOMS_CABANNES_UNSCALED(N,L,1)
                  PHAS_CB = PHAS_CB + BETA * SS_PLEG_DN(1,L)
                ENDDO
                CAB = PHAS_CB * DN_EMULTIPLIERS(N)
! @@@ Robfix 06/01/11. Add TMS
                CAB = CAB * TMS(N)

!  Gain term

                GAIN = ZERO
                DO W = 1, NPOINTS_MONO
                  RATIO   = FLUXES_RANKED(W)/FLUXES_RANKED(W_EXCIT)
                  CONTRIB = OMEGAMOMS_RRSGAIN_UNSCALED(N,0,W) + &
                   OMEGAMOMS_RRSGAIN_UNSCALED(N,2,W) * SS_PLEG_DN(1,2)
                  PHAS_LG(W)  = CONTRIB * RATIO
                  GAIN = GAIN + PHAS_LG(W) * DN_IMULTIPLIERS(N,W)
                ENDDO
! @@@ Robfix 06/01/11. Add TMS
                GAIN = GAIN * TMS(N)

!  Whole layer contribution

                ISCAT_DN(N) = CAB + GAIN

!  End loop

              ENDDO

!  End inelastic scattering options

             ENDIF
            ENDIF

!  End Downwelling clause

          ENDIF

!  Recurrence relation for the UPWELLING intensity
!  ===============================================

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

!  initialize cumulative source term = F.A.mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux, A = albedo
!    Only require the Stokes total intensity component
!    This contribution is purely elastic.
!    This line Replaced: X0_BOA = DCOS(SOLAR_ANGLE * DEG_TO_RAD)
!     Introduction of the Adjusted Gometry value, 18 March 2011

            NC =  0
            IF ( DO_INCLUDE_SURFACE ) THEN
              X0_BOA = DCOS(SOLAR_ANGLES_ADJUST(UM,IA) * DEG_TO_RAD)
              FACTOR = REFLEC
              BOA_ATTN = 4.0d0 * X0_BOA * ATTN(NLAYERS,W_EXCIT)
              ESS_CUMSOURCE_UP(NC) = FACTOR * BOA_ATTN
            ELSE
              ESS_CUMSOURCE_UP(NC) = ZERO
            ENDIF
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
      END SUBROUTINE LRRS_SSCORR_OUTGOING_MONO_1

!

      SUBROUTINE OUTGOING_SPHERGEOM_FINE_UP_1 &
        ( MAXLAYERS, MAXFINE, DO_FINE, NLAYERS, NFINE, &
          HEIGHTS, ERADIUS, ALPHA_BOA, THETA_BOA, PHI_BOA, &
          SUNPATHS,      RADII,      NTRAVERSE,      ALPHA_ALL, &
          SUNPATHS_FINE, RADII_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
          LOSPATHS, THETA_ALL, PHI_ALL, COSSCAT_UP, &
          FAIL, MESSAGE )

      IMPLICIT NONE

!  Completely stand-alone geometry routine for the outgoing correction
!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

!  This routine has the fine gridding treatment
!  Version 2.1. September 2007, Partial layer geometries added
!  Version 2.2. July 2009, This version without partials

!  Separate Downwelling and Upwelling routines, 15 February 2011

!  inputs

      INTEGER, INTENT(IN) ::          maxlayers, maxfine
      INTEGER, INTENT(IN) ::          nlayers, nfine
      LOGICAL, INTENT(IN) ::          do_fine
      REAL(FPK), INTENT(IN) :: eradius, heights (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_boa, theta_boa

!  modified inputs

      REAL(FPK), INTENT(INOUT) :: phi_boa

!  main outputs (geometry)

      INTEGER, INTENT(OUT) ::          ntraverse   (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: sunpaths   (0:maxlayers,maxlayers)
      REAL(FPK), INTENT(OUT) :: radii      (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: alpha_all  (0:maxlayers)

!  Fine level output (geometry)

      INTEGER, INTENT(OUT) ::          ntraverse_fine(maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: sunpaths_fine (maxlayers,maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: radii_fine    (maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: alpha_fine    (maxlayers,maxfine)

!  Other (incidental) geometrical output

      REAL(FPK), INTENT(OUT) :: lospaths(maxlayers)
      REAL(FPK), INTENT(OUT) :: theta_all  (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: phi_all    (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: cosscat_up (0:maxlayers)

!  Status output

      LOGICAL, INTENT(OUT) ::          fail
      CHARACTER (LEN=*), INTENT(OUT) ::    message

!  Local

      LOGICAL ::          direct_sun
      INTEGER ::          n, k, krad, n1
      REAL(FPK) :: deg_to_rad, ex, ey, ez, px, py, pz
      REAL(FPK) :: salpha_boa, calpha_boa, sphi_boa, pie, pi2
      REAL(FPK) :: stheta_boa, ctheta_boa, cphi_boa
      REAL(FPK) :: ksi, cksi, sksi, xicum, tangr, fac
      REAL(FPK) :: ctheta, stheta, calpha, salpha, cphi
      REAL(FPK) :: b, sth0, th0, ks1, sth1, th1

!  Local arrays associated with fine grid output

      INTEGER ::          j

      INTEGER, parameter   ::          maxlocalfine = 20

      LOGICAL ::          direct_sunf(maxlocalfine)
      REAL(FPK) :: difz, dfine1, saf, xicum0, path
      REAL(FPK) :: thetaf(maxlocalfine), xicumf, difa
      REAL(FPK) :: cthetaf(maxlocalfine)
      REAL(FPK) :: sthetaf(maxlocalfine)
      REAL(FPK) :: ksif(maxlocalfine)

!  Initialise output

      fail = .false.
      message = ' '

!  check range of inputs

      if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
        message = 'boa LOS angle outside range [0,90])'
        fail    = .true.
        return
      endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 10 January 2011------------------------- V. Natraj
!      if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
      if ( phi_boa.lt.0.0d0 )   phi_boa = 360.0d0 + phi_boa
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr

!  IMPORTANT, DO NOT LIMIT THE AZIMUTH TO < 180 degrees

!      if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa    ! OLD
       if ( phi_boa.gt.360.0d0 ) phi_boa = phi_boa - 360.0d0    ! NEW

!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
        message = 'boa SZA angle outside range [0,90])'
        fail    = .true.
        return
      endif
      if ( do_fine ) then
       if ( nfine.gt.maxlocalfine ) then
         message = 'local finelayer dimensioning insufficient'
         fail    = .true.
         return
       endif
      endif

!  zero the sun paths
!  Initialize number of layers traversed  (nominal conditions)

      do n = 0, nlayers
        ntraverse(n) = n
        do k = 1, nlayers
         sunpaths(n,k) = 0.0d0
        enddo
      enddo

!  Zero the fine data paths

      if ( do_fine ) then
       dfine1 = dble(nfine) + 1
       do n = 1, nlayers
        do j = 1, nfine
         ntraverse_fine(n,j) = n
         do k = 1, nlayers
          sunpaths_fine(n,k,j) = 0.0d0
         enddo
        enddo
       enddo
      endif

!  start at BOA

      pie        = dacos(-1.0d0)
      deg_to_rad = pie / 180.0d0
      pi2        = 2.0d0 * pie

      alpha_all(nlayers) = alpha_boa * deg_to_rad
      theta_all(nlayers) = theta_boa * deg_to_rad
      phi_all(nlayers)   = phi_boa   * deg_to_rad

!  Cosine of scattering angle at boa

      salpha_boa = dsin(alpha_all(nlayers))
      calpha_boa = dcos(alpha_all(nlayers))
      stheta_boa = dsin(theta_all(nlayers))
      ctheta_boa = dcos(theta_all(nlayers))
      cphi_boa   = dcos(phi_all(nlayers))
      sphi_boa   = dsin(phi_all(nlayers))

      cosscat_up (nlayers) = - calpha_boa * ctheta_boa + &
                               salpha_boa * stheta_boa * cphi_boa

!  Radii
!  -----

!  layer levels

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Fine levels

      if ( do_fine ) then
        do n = 1, nlayers
          difz = (radii(n-1)-radii(n))/dfine1
          do j = 1, nfine
            radii_fine(n,j) = radii(n) + difz * dble(j)
          enddo
        enddo
      endif

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit.
!    (This is the same as the regular pseudo-spherical )

      if ( salpha_boa.eq.0.0d0 ) then

!  WHOLE LAYER and FINE divisions
!  ------------------------------

!  Start layer loop, working upwards

        do n = nlayers,1,-1

!  set main output.

          alpha_all(n-1)   = alpha_all(n)
          theta_all(n-1)   = theta_all(n)
          phi_all(n-1)     = phi_all(n)
          cosscat_up(n-1) = cosscat_up(n)
          lospaths(n) = radii(n-1)-radii(n)
          if ( do_fine ) then
            do j = 1, nfine
              alpha_fine(n,j) = 0.0d0
            enddo
          endif

!  Overhead sun

          if (stheta_boa.eq.0.0d0 ) then
            do k = n, 1, -1
              sunpaths(n,k) = radii(k-1)-radii(k)
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                do k = n - 1, 1, -1
                  sunpaths_fine(n,k,j) = radii(k-1)-radii(k)
                enddo
                sunpaths_fine(n,n,j) = radii(n-1)-radii_fine(n,j)
              enddo
            endif
          endif

!  Non-overhead sun
!  Main output of solar paths
!  Solar path distances for fine output

          if (stheta_boa.gt.0.0d0 ) then
            sth0 = stheta_boa
            th0  = theta_all(n)
            do k = n, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                sth0 = stheta_boa
                th0  = theta_all(n)
                sth1 = sth0*radii_fine(n,j)/radii(n-1)
                th1  = dasin(sth1)
                ks1  = th0-th1
                sunpaths_fine(n,n,j) = dsin(ks1)*radii_fine(n,j)/sth1
                sth0 = sth1
                th0  = th1
                do k = n-1, 1, -1
                  sth1 = sth0*radii(k)/radii(k-1)
                  th1  = dasin(sth1)
                  ks1  = th0-th1
                  sunpaths_fine(n,k,j) = dsin(ks1)*radii(k)/sth1
                  sth0 = sth1
                  th0  = th1
                enddo
              enddo
            endif
          endif

!  End main layer loop

        enddo

!  Return, as everything now done

        return

!  end regular pseudo-spherical clause, LOS is zero

      endif

!  Outgoing spehricity geometry
!  ============================

!  define Unit solar vector at BOA

      ex = - stheta_boa * cphi_boa
      ey = - stheta_boa * sphi_boa
      ez = - ctheta_boa

!  Sun paths, boa geometry, always directly illuminated

      if ( stheta_boa.eq.0.0d0 ) then
        do k = nlayers, 1, -1
          sunpaths(nlayers,k) = radii(k-1)-radii(k)
        enddo
      else
        sth0 = stheta_boa
        th0  = theta_all(nlayers)
        do k = nlayers, 1, -1
          sth1 = sth0*radii(k)/radii(k-1)
          th1  = dasin(sth1)
          ks1  = th0-th1
          sunpaths(nlayers,k) = dsin(ks1)*radii(k)/sth1
          sth0 = sth1
          th0  = th1
        enddo
      endif

!  Check single illumination
!      --Not required, now we have the tangent point treatment
!      if (stheta_boa.gt.0.0d0 ) then
!        xicum = dasin(radii(nlayers)*salpha_boa/radii(0))
!        xicum = alpha_all(nlayers)-xicum
!        px = - radii(0) * dsin(xicum)
!        py = 0.0d0
!        pz =   radii(0) * dcos(xicum)
!        b = ex*px + ey*py + ez*pz
!        ctheta = -b/radii(0)
!        if ( ctheta.le.0.0d0 ) then
!          write(*,*)'limit value = ',90.0d0-xicum/deg_to_rad
!        endif
!      endif

!  initialise los cumulative angle

      xicum  = 0.0d0

!  set TOA direct illumination flag

      direct_sun = .true.
      if ( do_fine ) then
        do j = 1, nfine
          direct_sunf(j) = .true.
        enddo
      endif

!  Start loop over positions (layer upper boundaries)

      do n = nlayers - 1, 0, -1

!  Next level up

        n1 = n + 1

!  Los angles at level boundaries

        salpha = radii(nlayers) * salpha_boa / radii(n)
        alpha_all(n)  = dasin(salpha)
        calpha = dcos(alpha_all(n))

!  Lospaths

        ksi = alpha_all(n1) - alpha_all(n)
        sksi = dsin(ksi)
        cksi = dcos(ksi)
        lospaths(n1) = sksi * radii(n1) / salpha
        xicum0 = xicum
        xicum  = xicum + ksi

!  Fine grid lospath output (angle and radius)
!    Locally save the earth-center angle ksif

        if ( do_fine ) then
          difa = (alpha_all(n1)-alpha_all(n))/dfine1
          do j = 1, nfine
            alpha_fine(n1,j) = alpha_all(n1) - difa * dble(j)
            saf = dsin(alpha_fine(n1,j))
            radii_fine(n1,j) = salpha_boa * radii(nlayers) / saf
            ksif(j) = alpha_all(n1) - alpha_fine(n1,j)
          enddo
        endif

!  Sun angles for the Direct Nadir case

        if (stheta_boa.eq.0.0d0 ) then
         theta_all(n) = xicum
         ctheta = dcos(theta_all(n))
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             thetaf(j)  = xicum0 + ksif(j)
             cthetaf(j) = dcos(thetaf(j))
             sthetaf(j) = dsqrt(1.0d0-ctheta*ctheta)
           enddo
         endif
        endif

!  Sun angles for the general case
!    Local save of angles, cosines, sines and  illumination flags

        if (stheta_boa.gt.0.0d0 ) then
         px = - radii(n) * dsin(xicum)
         py = 0.0d0
         pz =   radii(n) * dcos(xicum)
         b = ex*px + ey*py + ez*pz
         ctheta = -b/radii(n)
         direct_sun = (direct_sun.and.ctheta.ge.0.d0)
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         theta_all(n) = dacos(ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             xicumf  = xicum0 + ksif(j)
             px = - radii_fine(n1,j) * dsin(xicumf)
             py = 0.0d0
             pz =   radii_fine(n1,j) * dcos(xicumf)
             b  = ex*px + ey*py + ez*pz
             cthetaf(j) = -b/radii_fine(n1,j)
             direct_sunf(j) = (direct_sunf(j).and.cthetaf(j).ge.0.d0)
             sthetaf(j) = dsqrt(1.0d0-cthetaf(j)*cthetaf(j))
             thetaf(j)  = dacos(cthetaf(j))
           enddo
         endif
        endif

!  Unit vector f2(i) perpendicular to OP but in plane of path
!  projection of f2(i) on solar path gives the relative azimuth at P
!        f2x = dsin(xicum)
!        f2y = 0.0d0
!        f2z = dcos(xicum)
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        if ( cphi.gt.1.0d0)  cphi = 1.0d0
!        if ( cphi.lt.-1.0d0) cphi = -1.0d0
! ********************************************* Apparently not correct

!  Fix phi by using constancy of scatter angle
!  Only for the scattering up directions..................

        cosscat_up(n) = cosscat_up(n+1)
        if (stheta_boa.eq.0.0d0 ) then
          phi_all(n)     = phi_all(n+1)
        else
         cphi = (cosscat_up(n)+calpha*ctheta)/stheta/salpha
         if ( cphi.gt.1.0d0) cphi = 1.0d0
         if ( cphi.lt.-1.0d0) cphi = -1.0d0
         phi_all(n)     = dacos(cphi)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr

!  IF AZIMUTH > 180 degrees. Must ensure consistency, add this line.

         if ( phi_boa.gt.180.0d0 ) phi_all(n) = pi2 - phi_all(n)

 !  B G CORRECTED 01 October 2010------------------------- RTS, R. Spurr
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        endif

!  Sun paths, Direct sun at layer top
!  ==================================

!   Means that the SZA at layer top is < 90.
!    ===> SZA < 90 for all fine points in layer beneath

        if ( direct_sun ) then

!  Work up from level n to TOA
!    Layer top calculation gets left out at TOA

         if ( n .gt. 0 ) then
          sth0 = stheta
          th0  = theta_all(n)
          do k = n, 1, -1
           sth1 = sth0*radii(k)/radii(k-1)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
           sth0 = sth1
           th0  = th1
          enddo
         endif

! DBG         write(*,*)'regular',n1,(sunpaths(n1,k),k=n1,1,-1)

!  Fine grid calculation is always required
!  ----------------------------------------

!   Start at the grid point on LOS path and work upwards,
!     first across partial layer to upper boundary, then continue
!     upwards to TOA on a whole layer basis. Sine rule.

         if ( do_fine ) then
          do j = 1, nfine
           sth0 = sthetaf(j)
           th0  = thetaf(j)
           sth1 = sth0*radii_fine(n1,j)/radii(n)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sunpaths_fine(n1,n1,j) = dsin(ks1)*radii_fine(n1,j)/sth1
           sth0 = sth1
           th0  = th1
           do k = n, 1, -1
            sth1 = sth0*radii(k)/radii(k-1)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k)/sth1
            sth0 = sth1
            th0  = th1
           enddo
!  DBG           write(*,*)'regular',n1,j,(sunpaths_fine(n1,k,j),k=n1,1,
          enddo

!  DBG         write(*,*)'regular',n,(sunpaths(n,k),k=n,1,-1)

         endif

!  Complete direct sun computations

        endif

!  Sun paths, Not direct sun , with tangent point
!  ==============================================

!  Although layer top has a tangent point, not all of the fine-grid
!   points will have a tangent point.

        if (.not.direct_sun ) then

!  First do the layer-top calculation.
!  -----------------------------------

!  TANGR = tangent point radius.

         tangr = stheta*radii(n)

!  ntraverse(n) is the number of layers traversed by ray.

         krad = nlayers
         do while (tangr.gt.radii(krad))
          krad = krad - 1
         enddo
         ntraverse(n) = krad + 1

!  Start at the TOA angles

         sth0 = tangr/radii(0)
         th0 = dasin(sth0)

!  Work downwards from TOA (sine rule) to level immediately above
!  the tangent layer. Don't forget to double the path length for
!  any layers which are traversed twice.

         do k = 1, krad
          sth1 = radii(k-1)*sth0/radii(k)
          th1 = dasin(sth1)
          ks1 = th1-th0
          fac = 1.0d0
          if ( k.gt.n) fac = 2.0d0
          sunpaths(n,k) = fac*dsin(ks1)*radii(k-1)/sth1
          sth0 = sth1
          th0  = th1
         enddo

!  Tangent layer path length. Twice again. The following check is good.
!  check       write(*,*)tangr/dtan(th1),radii(krad)*dcos(th1)

         sunpaths(n,krad+1)=2.0d0*radii(krad)*dcos(th1)

! DBG         write(*,*)'--tangent',n1,(sunpaths(n1,k),k=ntraverse(n1),1

!  Fine layer development (slightly different)
!  ---------------------

         if ( do_fine ) then
          do j = 1, nfine

!  If there is no tangent point, repeat calculation above.

           if ( direct_sunf(j) ) then
            sth0 = sthetaf(j)
            th0  = thetaf(j)
            sth1 = sth0*radii_fine(n1,j)/radii(n)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,n1,j) = dsin(ks1)*radii_fine(n1,j)/sth1
            sth0 = sth1
            th0  = th1
            do k = n, 1, -1
             sth1 = sth0*radii(k)/radii(k-1)
             th1  = dasin(sth1)
             ks1  = th0-th1
             sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k)/sth1
             sth0 = sth1
             th0  = th1
            enddo

! DBG           write(*,*)'regular',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Fine grid Calculation with tangent point
!  ----------------------------------------

           else

!  Local tangent radius and number of layers traversed

            tangr = sthetaf(j)*radii_fine(n1,j)
            krad = nlayers
            do while (tangr.gt.radii(krad))
             krad = krad - 1
            enddo
            ntraverse_fine(n1,j) = krad + 1

!  Start again at TOA

            sth0 = tangr/radii(0)
            th0  = dasin(sth0)

!  Work down to the level n

            do k = 1, n
             sth1 = radii(k-1)*sth0/radii(k)
             th1 = dasin(sth1)
             ks1 = th1-th0
             sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k-1)/sth1
             sth0 = sth1
             th0  = th1
            enddo

!  In layer below level n, Work down to level of fine grid radius
!     (single contribution)

            sth1 = radii(n)*sth0/radii_fine(n1,j)
            th1 = dasin(sth1)
            ks1 = th1-th0
            path = dsin(ks1)*radii(n)/sth1
            sth0 = sth1
            th0  = th1

!  In layer below level n, Complete the path down to the tangent point
!    (double contribution). Finish.

            if ( krad.eq.n ) then

              path = path + 2.0d0*radii_fine(n1,j)*dcos(th1)
              sunpaths_fine(n1,krad+1,j) = path

!  Check    write(*,*)'1',tangr/dtan(th1),radii_fine(n1,j)*dcos(th1)

!  In layer below level n, Need to go down further.
!    (from now on, we need double contributions)
!     --- Continue the path to the bottom of this layer
!     --- Continue down by whole-layer steps until reach level above Tan
!     --- Complete the path down to the tangent point

            else
              sth1 = radii_fine(n1,j)*sth0/radii(n1)
              th1 = dasin(sth1)
              ks1 = th1-th0
              path = path + 2.0d0 * dsin(ks1)*radii_fine(n1,j)/sth1
              sunpaths_fine(n1,n1,j) = path
              sth0 = sth1
              th0  = th1
              do k = n1 + 1, krad
               sth1 = radii(k-1)*sth0/radii(k)
               th1 = dasin(sth1)
               ks1 = th1-th0
               sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k-1)/sth1
               sth0 = sth1
               th0  = th1
              enddo
              sunpaths_fine(n1,krad+1,j)=2.0d0*radii(krad)*dcos(th1)
!  Check 2    write(*,*)'2',tangr/dtan(th1),radii(krad)*dcos(th1)
            endif

! DBG           write(*,*)'tangent',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Complete tangent point clause

           endif

!  Complete fine-grid loop

          enddo

!  DBG        write(*,*)'tangent',n,(sunpaths(n,k),k=ntraverse(n),1,-1)

!  Complete tangent point calculation

         endif

        endif

!  End layer loop

      enddo

!  Finish

      return
      END SUBROUTINE OUTGOING_SPHERGEOM_FINE_UP_1

!

      SUBROUTINE OUTGOING_SPHERGEOM_FINE_DN_1 &
        ( MAXLAYERS, MAXFINE, DO_FINE, NLAYERS, NFINE, &
          HEIGHTS, ERADIUS, ALPHA_BOA, THETA_BOA, PHI_BOA, &
          SUNPATHS,      RADII,      NTRAVERSE,      ALPHA_ALL, &
          SUNPATHS_FINE, RADII_FINE, NTRAVERSE_FINE, ALPHA_FINE, &
          LOSPATHS, THETA_ALL, PHI_ALL, COSSCAT_DN, &
          FAIL, MESSAGE )

      IMPLICIT NONE

!  Completely stand-alone geometry routine for the outgoing correction
!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

!  This routine has the fine gridding treatment
!  Version 2.1. September 2007, Partial layer geometries added
!  Version 2.2. July 2009, This version without partials

!  Separate Downwelling and Upwelling routines, 15 February 2011

!  inputs

      INTEGER, INTENT(IN) ::          maxlayers, maxfine
      INTEGER, INTENT(IN) ::          nlayers, nfine
      LOGICAL, INTENT(IN) ::          do_fine
      REAL(FPK), INTENT(IN) :: eradius, heights (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_boa, theta_boa

!  modified inputs

      REAL(FPK), INTENT(INOUT) :: phi_boa

!  main outputs (geometry)

      INTEGER, INTENT(OUT) ::          ntraverse   (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: sunpaths   (0:maxlayers,maxlayers)
      REAL(FPK), INTENT(OUT) :: radii      (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: alpha_all  (0:maxlayers)

!  Fine level output (geometry)

      INTEGER, INTENT(OUT) ::          ntraverse_fine(maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: sunpaths_fine (maxlayers,maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: radii_fine    (maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: alpha_fine    (maxlayers,maxfine)

!  Other (incidental) geometrical output

      REAL(FPK), INTENT(OUT) :: lospaths(maxlayers)
      REAL(FPK), INTENT(OUT) :: theta_all  (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: phi_all    (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: cosscat_dn (0:maxlayers)

!  Status output

      LOGICAL, INTENT(OUT) ::          fail
      CHARACTER (LEN=*), INTENT(OUT) ::    message

!  Local

      LOGICAL ::          direct_sun
      INTEGER ::          n, k, krad, n1
      REAL(FPK) :: deg_to_rad, ex, ey, ez, px, py, pz
      REAL(FPK) :: salpha_boa, calpha_boa, sphi_boa, pie, pi2
      REAL(FPK) :: stheta_boa, ctheta_boa, cphi_boa
      REAL(FPK) :: ksi, cksi, sksi, xicum, tangr, fac
      REAL(FPK) :: ctheta, stheta, calpha, salpha, cphi
      REAL(FPK) :: b, sth0, th0, ks1, sth1, th1

!  Local arrays associated with fine grid output

      INTEGER ::          j

      INTEGER, parameter   ::          maxlocalfine = 20

      LOGICAL ::          direct_sunf(maxlocalfine)
      REAL(FPK) :: difz, dfine1, saf, xicum0, path
      REAL(FPK) :: thetaf(maxlocalfine), xicumf, difa
      REAL(FPK) :: cthetaf(maxlocalfine)
      REAL(FPK) :: sthetaf(maxlocalfine)
      REAL(FPK) :: ksif(maxlocalfine)

!  Initialise output

      fail = .false.
      message = ' '

!  check range of inputs

      if ( alpha_boa.ge.90.0d0.or.alpha_boa.lt.0.0d0 ) then
        message = 'boa LOS angle outside range [0,90])'
        fail    = .true.
        return
      endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 10 January 2011------------------------- V. Natraj
!      if ( phi_boa.lt.0.0d0 )   phi_boa = - phi_boa
      if ( phi_boa.lt.0.0d0 )   phi_boa = 360.0d0 + phi_boa
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr

!  IMPORTANT, DO NOT LIMIT THE AZIMUTH TO < 180 degrees

!      if ( phi_boa.gt.180.0d0 ) phi_boa = 360.0d0 - phi_boa    ! OLD
       if ( phi_boa.gt.360.0d0 ) phi_boa = phi_boa - 360.0d0    ! NEW

!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
        message = 'boa SZA angle outside range [0,90])'
        fail    = .true.
        return
      endif
      if ( do_fine ) then
       if ( nfine.gt.maxlocalfine ) then
         message = 'local finelayer dimensioning insufficient'
         fail    = .true.
         return
       endif
      endif

!  zero the sun paths
!  Initialize number of layers traversed  (nominal conditions)

      do n = 0, nlayers
        ntraverse(n) = n
        do k = 1, nlayers
         sunpaths(n,k) = 0.0d0
        enddo
      enddo

!  Zero the fine data paths

      if ( do_fine ) then
       dfine1 = dble(nfine) + 1
       do n = 1, nlayers
        do j = 1, nfine
         ntraverse_fine(n,j) = n
         do k = 1, nlayers
          sunpaths_fine(n,k,j) = 0.0d0
         enddo
        enddo
       enddo
      endif

!  start at BOA

      pie        = dacos(-1.0d0)
      deg_to_rad = pie / 180.0d0
      pi2        = 2.0d0 * pie

      alpha_all(nlayers) = alpha_boa * deg_to_rad
      theta_all(nlayers) = theta_boa * deg_to_rad
      phi_all(nlayers)   = phi_boa   * deg_to_rad

!  Cosine of scattering angle at boa

      salpha_boa = dsin(alpha_all(nlayers))
      calpha_boa = dcos(alpha_all(nlayers))
      stheta_boa = dsin(theta_all(nlayers))
      ctheta_boa = dcos(theta_all(nlayers))
      cphi_boa   = dcos(phi_all(nlayers))
      sphi_boa   = dsin(phi_all(nlayers))

      cosscat_dn (nlayers) = + calpha_boa * ctheta_boa + &
                               salpha_boa * stheta_boa * cphi_boa

!  Radii
!  -----

!  layer levels

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Fine levels

      if ( do_fine ) then
        do n = 1, nlayers
          difz = (radii(n-1)-radii(n))/dfine1
          do j = 1, nfine
            radii_fine(n,j) = radii(n) + difz * dble(j)
          enddo
        enddo
      endif

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit.
!    (This is the same as the regular pseudo-spherical )

      if ( salpha_boa.eq.0.0d0 ) then

!  WHOLE LAYER and FINE divisions
!  ------------------------------

!  Start layer loop, working upwards

        do n = nlayers,1,-1

!  set main output.

          alpha_all(n-1)   = alpha_all(n)
          theta_all(n-1)   = theta_all(n)
          phi_all(n-1)     = phi_all(n)
          cosscat_dn(n-1) = cosscat_dn(n)
          lospaths(n) = radii(n-1)-radii(n)
          if ( do_fine ) then
            do j = 1, nfine
              alpha_fine(n,j) = 0.0d0
            enddo
          endif

!  Overhead sun

          if (stheta_boa.eq.0.0d0 ) then
            do k = n, 1, -1
              sunpaths(n,k) = radii(k-1)-radii(k)
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                do k = n - 1, 1, -1
                  sunpaths_fine(n,k,j) = radii(k-1)-radii(k)
                enddo
                sunpaths_fine(n,n,j) = radii(n-1)-radii_fine(n,j)
              enddo
            endif
          endif

!  Non-overhead sun
!  Main output of solar paths
!  Solar path distances for fine output

          if (stheta_boa.gt.0.0d0 ) then
            sth0 = stheta_boa
            th0  = theta_all(n)
            do k = n, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = dasin(sth1)
              ks1  = th0-th1
              sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                sth0 = stheta_boa
                th0  = theta_all(n)
                sth1 = sth0*radii_fine(n,j)/radii(n-1)
                th1  = dasin(sth1)
                ks1  = th0-th1
                sunpaths_fine(n,n,j) = dsin(ks1)*radii_fine(n,j)/sth1
                sth0 = sth1
                th0  = th1
                do k = n-1, 1, -1
                  sth1 = sth0*radii(k)/radii(k-1)
                  th1  = dasin(sth1)
                  ks1  = th0-th1
                  sunpaths_fine(n,k,j) = dsin(ks1)*radii(k)/sth1
                  sth0 = sth1
                  th0  = th1
                enddo
              enddo
            endif
          endif

!  End main layer loop

        enddo

!  Return, as everything now done

        return

!  end regular pseudo-spherical clause, LOS is zero

      endif

!  Outgoing spehricity geometry
!  ============================

!  define Unit solar vector at BOA
!   Note change of sign for ex here.

      ex = + stheta_boa * cphi_boa
      ey = - stheta_boa * sphi_boa
      ez = - ctheta_boa

!  Sun paths, boa geometry, always directly illuminated

      if ( stheta_boa.eq.0.0d0 ) then
        do k = nlayers, 1, -1
          sunpaths(nlayers,k) = radii(k-1)-radii(k)
        enddo
      else
        sth0 = stheta_boa
        th0  = theta_all(nlayers)
        do k = nlayers, 1, -1
          sth1 = sth0*radii(k)/radii(k-1)
          th1  = dasin(sth1)
          ks1  = th0-th1
          sunpaths(nlayers,k) = dsin(ks1)*radii(k)/sth1
          sth0 = sth1
          th0  = th1
        enddo
      endif

!  Check single illumination
!      --Not required, now we have the tangent point treatment
!      if (stheta_boa.gt.0.0d0 ) then
!        xicum = dasin(radii(nlayers)*salpha_boa/radii(0))
!        xicum = alpha_all(nlayers)-xicum
!        px = - radii(0) * dsin(xicum)
!        py = 0.0d0
!        pz =   radii(0) * dcos(xicum)
!        b = ex*px + ey*py + ez*pz
!        ctheta = -b/radii(0)
!        if ( ctheta.le.0.0d0 ) then
!          write(*,*)'limit value = ',90.0d0-xicum/deg_to_rad
!        endif
!      endif

!  initialise los cumulative angle

      xicum  = 0.0d0

!  set TOA direct illumination flag

      direct_sun = .true.
      if ( do_fine ) then
        do j = 1, nfine
          direct_sunf(j) = .true.
        enddo
      endif

!  Start loop over positions (layer upper boundaries)

      do n = nlayers - 1, 0, -1

!  Next level up

        n1 = n + 1

!  Los angles at level boundaries

        salpha = radii(nlayers) * salpha_boa / radii(n)
        alpha_all(n)  = dasin(salpha)
        calpha = dcos(alpha_all(n))

!  Lospaths

        ksi = alpha_all(n1) - alpha_all(n)
        sksi = dsin(ksi)
        cksi = dcos(ksi)
        lospaths(n1) = sksi * radii(n1) / salpha
        xicum0 = xicum
        xicum  = xicum + ksi

!  Fine grid lospath output (angle and radius)
!    Locally save the earth-center angle ksif

        if ( do_fine ) then
          difa = (alpha_all(n1)-alpha_all(n))/dfine1
          do j = 1, nfine
            alpha_fine(n1,j) = alpha_all(n1) - difa * dble(j)
            saf = dsin(alpha_fine(n1,j))
            radii_fine(n1,j) = salpha_boa * radii(nlayers) / saf
            ksif(j) = alpha_all(n1) - alpha_fine(n1,j)
          enddo
        endif

!  Sun angles for the Direct Nadir case

        if (stheta_boa.eq.0.0d0 ) then
         theta_all(n) = xicum
         ctheta = dcos(theta_all(n))
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             thetaf(j)  = xicum0 + ksif(j)
             cthetaf(j) = dcos(thetaf(j))
             sthetaf(j) = dsqrt(1.0d0-ctheta*ctheta)
           enddo
         endif
        endif

!  Sun angles for the general case
!    Local save of angles, cosines, sines and  illumination flags

        if (stheta_boa.gt.0.0d0 ) then
         px = - radii(n) * dsin(xicum)
         py = 0.0d0
         pz =   radii(n) * dcos(xicum)
         b = ex*px + ey*py + ez*pz
         ctheta = -b/radii(n)
         direct_sun = (direct_sun.and.ctheta.ge.0.d0)
         stheta = dsqrt(1.0d0-ctheta*ctheta)
         theta_all(n) = dacos(ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             xicumf  = xicum0 + ksif(j)
             px = - radii_fine(n1,j) * dsin(xicumf)
             py = 0.0d0
             pz =   radii_fine(n1,j) * dcos(xicumf)
             b  = ex*px + ey*py + ez*pz
             cthetaf(j) = -b/radii_fine(n1,j)
             direct_sunf(j) = (direct_sunf(j).and.cthetaf(j).ge.0.d0)
             sthetaf(j) = dsqrt(1.0d0-cthetaf(j)*cthetaf(j))
             thetaf(j)  = dacos(cthetaf(j))
           enddo
         endif
        endif

!  Unit vector f2(i) perpendicular to OP but in plane of path
!  projection of f2(i) on solar path gives the relative azimuth at P
!        f2x = dsin(xicum)
!        f2y = 0.0d0
!        f2z = dcos(xicum)
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        if ( cphi.gt.1.0d0)  cphi = 1.0d0
!        if ( cphi.lt.-1.0d0) cphi = -1.0d0
! ********************************************* Apparently not correct

!  Fix phi by using constancy of scatter angle
!  Only for the scattering up directions..................

        cosscat_dn(n) = cosscat_dn(n+1)
        if (stheta_boa.eq.0.0d0 ) then
          phi_all(n)     = phi_all(n+1)
        else
         cphi = (cosscat_dn(n)-calpha*ctheta)/stheta/salpha
         if ( cphi.gt.1.0d0) cphi = 1.0d0
         if ( cphi.lt.-1.0d0) cphi = -1.0d0
         phi_all(n)     = dacos(cphi)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr

!  IF AZIMUTH > 180 degrees. Must ensure consistency, add this line.

         if ( phi_boa.gt.180.0d0 ) phi_all(n) = pi2 - phi_all(n)

 !  B G CORRECTED 01 October 2010------------------------- RTS, R. Spurr
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        endif

!  Sun paths, Direct sun at layer top
!  ==================================

!   Means that the SZA at layer top is < 90.
!    ===> SZA < 90 for all fine points in layer beneath

        if ( direct_sun ) then

!  Work up from level n to TOA
!    Layer top calculation gets left out at TOA

         if ( n .gt. 0 ) then
          sth0 = stheta
          th0  = theta_all(n)
          do k = n, 1, -1
           sth1 = sth0*radii(k)/radii(k-1)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sunpaths(n,k) = dsin(ks1)*radii(k)/sth1
           sth0 = sth1
           th0  = th1
          enddo
         endif

! DBG         write(*,*)'regular',n1,(sunpaths(n1,k),k=n1,1,-1)

!  Fine grid calculation is always required
!  ----------------------------------------

!   Start at the grid point on LOS path and work upwards,
!     first across partial layer to upper boundary, then continue
!     upwards to TOA on a whole layer basis. Sine rule.

         if ( do_fine ) then
          do j = 1, nfine
           sth0 = sthetaf(j)
           th0  = thetaf(j)
           sth1 = sth0*radii_fine(n1,j)/radii(n)
           th1  = dasin(sth1)
           ks1  = th0-th1
           sunpaths_fine(n1,n1,j) = dsin(ks1)*radii_fine(n1,j)/sth1
           sth0 = sth1
           th0  = th1
           do k = n, 1, -1
            sth1 = sth0*radii(k)/radii(k-1)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k)/sth1
            sth0 = sth1
            th0  = th1
           enddo
!  DBG           write(*,*)'regular',n1,j,(sunpaths_fine(n1,k,j),k=n1,1,
          enddo

!  DBG         write(*,*)'regular',n,(sunpaths(n,k),k=n,1,-1)

         endif

!  Complete direct sun computations

        endif

!  Sun paths, Not direct sun , with tangent point
!  ==============================================

!  Although layer top has a tangent point, not all of the fine-grid
!   points will have a tangent point.

        if (.not.direct_sun ) then

!  First do the layer-top calculation.
!  -----------------------------------

!  TANGR = tangent point radius.

         tangr = stheta*radii(n)

!  ntraverse(n) is the number of layers traversed by ray.

         krad = nlayers
         do while (tangr.gt.radii(krad))
          krad = krad - 1
         enddo
         ntraverse(n) = krad + 1

!  Start at the TOA angles

         sth0 = tangr/radii(0)
         th0 = dasin(sth0)

!  Work downwards from TOA (sine rule) to level immediately above
!  the tangent layer. Don't forget to double the path length for
!  any layers which are traversed twice.

         do k = 1, krad
          sth1 = radii(k-1)*sth0/radii(k)
          th1 = dasin(sth1)
          ks1 = th1-th0
          fac = 1.0d0
          if ( k.gt.n) fac = 2.0d0
          sunpaths(n,k) = fac*dsin(ks1)*radii(k-1)/sth1
          sth0 = sth1
          th0  = th1
         enddo

!  Tangent layer path length. Twice again. The following check is good.
!  check       write(*,*)tangr/dtan(th1),radii(krad)*dcos(th1)

         sunpaths(n,krad+1)=2.0d0*radii(krad)*dcos(th1)

! DBG         write(*,*)'--tangent',n1,(sunpaths(n1,k),k=ntraverse(n1),1

!  Fine layer development (slightly different)
!  ---------------------

         if ( do_fine ) then
          do j = 1, nfine

!  If there is no tangent point, repeat calculation above.

           if ( direct_sunf(j) ) then
            sth0 = sthetaf(j)
            th0  = thetaf(j)
            sth1 = sth0*radii_fine(n1,j)/radii(n)
            th1  = dasin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,n1,j) = dsin(ks1)*radii_fine(n1,j)/sth1
            sth0 = sth1
            th0  = th1
            do k = n, 1, -1
             sth1 = sth0*radii(k)/radii(k-1)
             th1  = dasin(sth1)
             ks1  = th0-th1
             sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k)/sth1
             sth0 = sth1
             th0  = th1
            enddo

! DBG           write(*,*)'regular',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Fine grid Calculation with tangent point
!  ----------------------------------------

           else

!  Local tangent radius and number of layers traversed

            tangr = sthetaf(j)*radii_fine(n1,j)
            krad = nlayers
            do while (tangr.gt.radii(krad))
             krad = krad - 1
            enddo
            ntraverse_fine(n1,j) = krad + 1

!  Start again at TOA

            sth0 = tangr/radii(0)
            th0  = dasin(sth0)

!  Work down to the level n

            do k = 1, n
             sth1 = radii(k-1)*sth0/radii(k)
             th1 = dasin(sth1)
             ks1 = th1-th0
             sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k-1)/sth1
             sth0 = sth1
             th0  = th1
            enddo

!  In layer below level n, Work down to level of fine grid radius
!     (single contribution)

            sth1 = radii(n)*sth0/radii_fine(n1,j)
            th1 = dasin(sth1)
            ks1 = th1-th0
            path = dsin(ks1)*radii(n)/sth1
            sth0 = sth1
            th0  = th1

!  In layer below level n, Complete the path down to the tangent point
!    (double contribution). Finish.

            if ( krad.eq.n ) then

              path = path + 2.0d0*radii_fine(n1,j)*dcos(th1)
              sunpaths_fine(n1,krad+1,j) = path

!  Check    write(*,*)'1',tangr/dtan(th1),radii_fine(n1,j)*dcos(th1)

!  In layer below level n, Need to go down further.
!    (from now on, we need double contributions)
!     --- Continue the path to the bottom of this layer
!     --- Continue down by whole-layer steps until reach level above Tan
!     --- Complete the path down to the tangent point

            else
              sth1 = radii_fine(n1,j)*sth0/radii(n1)
              th1 = dasin(sth1)
              ks1 = th1-th0
              path = path + 2.0d0 * dsin(ks1)*radii_fine(n1,j)/sth1
              sunpaths_fine(n1,n1,j) = path
              sth0 = sth1
              th0  = th1
              do k = n1 + 1, krad
               sth1 = radii(k-1)*sth0/radii(k)
               th1 = dasin(sth1)
               ks1 = th1-th0
               sunpaths_fine(n1,k,j) = dsin(ks1)*radii(k-1)/sth1
               sth0 = sth1
               th0  = th1
              enddo
              sunpaths_fine(n1,krad+1,j)=2.0d0*radii(krad)*dcos(th1)
!  Check 2    write(*,*)'2',tangr/dtan(th1),radii(krad)*dcos(th1)
            endif

! DBG           write(*,*)'tangent',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Complete tangent point clause

           endif

!  Complete fine-grid loop

          enddo

!  DBG        write(*,*)'tangent',n,(sunpaths(n,k),k=ntraverse(n),1,-1)

!  Complete tangent point calculation

         endif

        endif

!  End layer loop

      enddo

!  Finish

      return
      END SUBROUTINE OUTGOING_SPHERGEOM_FINE_DN_1

!

      SUBROUTINE OUTGOING_ATTENUATIONS_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, &
          N_RAMAN_POINTS, NLAYERS, NFINELAYERS, &
          EXTINCTION, SUNPATHS, NTRAVERSE, &
          SUNPATHS_FINE, NTRAVERSE_FINE, &
          ATTN, ATTN_FINE )

!  Does attenuations, Whole layers only. No partials.

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER, INTENT(IN) ::          maxlayers, maxfinelayers, maxpoints

!  control

      INTEGER, INTENT(IN) ::          nfinelayers, nlayers, n_raman_points

!  Whole layers

      INTEGER, INTENT(IN) ::          ntraverse(0:maxlayers)
      REAL(FPK), INTENT(IN) :: sunpaths(0:maxlayers,maxlayers)

!  Fine level

      INTEGER, INTENT(IN) ::          ntraverse_fine(maxlayers,maxfinelayers)
      REAL(FPK), INTENT(IN) :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)

!  Extinction

      REAL(FPK), INTENT(IN) :: extinction (maxlayers, maxpoints)

!  outputs
!  -------

      REAL(FPK), INTENT(OUT) :: attn      ( 0:maxlayers, maxpoints )
      REAL(FPK), INTENT(OUT) :: attn_fine &
               ( maxlayers, maxfinelayers, maxpoints )

!  help variables
!  --------------

      INTEGER ::          w, n, j, k
      REAL(FPK) :: tau

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

!  Finish

      return
      END SUBROUTINE OUTGOING_ATTENUATIONS_1

!

      SUBROUTINE OUTGOING_INTEGRATION_BIN_UP_1 &
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXBINS, &
          NPOINTS_INNER, OFFSET_INNER, NBINS, BINMAP, &
          NLAYERS, NFINELAYERS, EXTINCTION, &
          RADII, ALPHA_ALL, ALPHA_FINE, ATTN, ATTN_FINE, &
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS )

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

      IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER, INTENT(IN) ::          maxlayers, maxfinelayers
      INTEGER, INTENT(IN) ::          maxpoints, maxbins

!  control

      INTEGER, INTENT(IN) ::          nfinelayers, nlayers
      INTEGER, INTENT(IN) ::          npoints_inner, offset_inner
      INTEGER, INTENT(IN) ::          nbins  ( maxpoints )
      INTEGER, INTENT(IN) ::          binmap ( maxbins, maxpoints )

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

      INTEGER ::          w, b, wb, wr, n, j
      REAL(FPK) :: sum, salpha, calpha, dfine1, raycon, kn
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

      do n = nlayers, 1, -1

!  Save some quantities

        salpha = dsin(alpha_all(n-1))
        calpha = dcos(alpha_all(n-1))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  Trapezium-rule integrated sourceterm multipliers + transmittances
!  -----------------------------------------------------------------

        do w = 1, npoints_inner

!  set up

          wr = w + offset_inner
          kn = raycon * extinction(n,wr)
          skn = step * kn
          tran_1 = dexp ( - kn * ( cot_2 - cot_1 ) )
          csqt_1 = csq_1 * tran_1
          do j = 1, nfinelayers
            tranj(j) = dexp ( - kn * ( cot_2 - cot_fine(n,j) ) )
          enddo

!  Line of sight transmittance factor

          lostrans(n,w)    = tran_1

!  Elastic scattering multiplier

          func_2 = attn(n-1,wr) *  csq_2
          func_1 = attn(n,wr)   * csqt_1
          sum = 0.5d0 * ( func_1 + func_2 )
          do j = 1, nfinelayers
            func = attn_fine(n,j,wr) * tranj(j) * csq_fine(n,j)
            sum  = sum + func
          enddo
          emultipliers(n,w) = sum * skn

!  Inelastic scattering multipliers (by bins)

          do b = 1, nbins(w)
            wb = binmap(b,w)
            func_2 = attn(n-1,wb) *  csq_2
            func_1 = attn(n,wb)   * csqt_1
            sum = 0.5d0 * ( func_1 + func_2 )
            do j = 1, nfinelayers
              func = attn_fine(n,j,wb) * tranj(j) * csq_fine(n,j)
              sum  = sum + func
            enddo
            imultipliers(n,b,w) = sum * skn
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
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, MAXBINS, &
          NPOINTS_INNER, OFFSET_INNER, NBINS, BINMAP, &
          NLAYERS, NFINELAYERS, EXTINCTION, &
          RADII, ALPHA_ALL, ALPHA_FINE, ATTN, ATTN_FINE, &
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS )

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

       IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER, INTENT(IN) ::          maxlayers, maxfinelayers
      INTEGER, INTENT(IN) ::          maxpoints, maxbins

!  control

      INTEGER, INTENT(IN) ::          nfinelayers, nlayers
      INTEGER, INTENT(IN) ::          npoints_inner, offset_inner
      INTEGER, INTENT(IN) ::          nbins  ( maxpoints )
      INTEGER, INTENT(IN) ::          binmap ( maxbins, maxpoints )

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

      INTEGER ::          w, wb, b, wr, n, j
      REAL(FPK) :: sum, salpha, calpha, dfine1, raycon, kn
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

!  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work Down from the top of the atmosphere
!  ========================================

!  initialise

      n = 0
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

!  Save some quantities

      do n = 1, nlayers
        do j = nfinelayers, 1, -1
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = 1.0d0 / salpha / salpha
        enddo
      enddo

!  Start layer loop

      do n = 1, nlayers

!  Save some quantities

        salpha = dsin(alpha_all(n))
        calpha = dcos(alpha_all(n))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  Trapezium-rule integrated sourceterm multiplier + transmittance

        do w = 1, npoints_inner

!  set up

          wr = w + offset_inner
          kn = raycon * extinction(n,wr)
          skn = step * kn
          tran_1 = dexp ( - kn * ( cot_1 - cot_2 ) )
          csqt_1 = csq_1 * tran_1
          do j = 1, nfinelayers
            tranj(j) = dexp ( - kn * ( cot_fine(n,j) - cot_2 ) )
          enddo

!  Line of sight transmittance factor

          lostrans(n,w)    = tran_1

!  Elastic multipliers

          func_1 = attn(n-1,wr) * csqt_1
          func_2 = attn(n,wr)   * csq_2
          sum = 0.5d0 * ( func_1 + func_2 )
          do j = nfinelayers, 1, -1
            func = attn_fine(n,j,wr) * tranj(j) * csq_fine(n,j)
            sum = sum + func
          enddo
          emultipliers(n,w) = sum * step * kn

!  Inelastic multipliers

          do b = 1, nbins(w)
            wb = binmap(b,w)
            func_1 = attn(n-1,wb) * csqt_1
            func_2 = attn(n,wb)   * csq_2
            sum = 0.5d0 * ( func_1 + func_2 )
            do j = nfinelayers, 1, -1
              func = attn_fine(n,j,wb) * tranj(j) * csq_fine(n,j)
              sum = sum + func
            enddo
            imultipliers(n,b,w) = sum * step * kn
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
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, &
          NPOINTS_MONO, W_EXCIT, &
          NLAYERS, NFINELAYERS, EXTINCTION, &
          RADII, ALPHA_ALL, ALPHA_FINE, ATTN, ATTN_FINE, &
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS )

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

       IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER, INTENT(IN) ::          maxlayers, maxfinelayers
      INTEGER, INTENT(IN) ::          maxpoints

!  control

      INTEGER, INTENT(IN) ::          nfinelayers, nlayers
      INTEGER, INTENT(IN) ::          npoints_mono, w_excit

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

      INTEGER ::          w, n, j
      REAL(FPK) :: sum, salpha, calpha, dfine1, raycon, kn
      REAL(FPK) :: csq_1, csqt_1, csq_2, cot_1, cot_2
      REAL(FPK) :: func_1, func_2, tran_1, step, func, skn

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0D0

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        lostrans(n)     = 0.0d0
        emultipliers(n) = 0.0d0
        do w = 1, npoints_mono
          imultipliers(n,w) = 0.0d0
        enddo
      enddo

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

      do n = nlayers, 1, -1

!  Save some quantities

        salpha = dsin(alpha_all(n-1))
        calpha = dcos(alpha_all(n-1))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  Trapezium-rule integrated sourceterm multipliers + transmittances
!  -----------------------------------------------------------------

!  set up

        kn = raycon * extinction(n,w_excit)
        skn = step * kn
        tran_1 = dexp ( - kn * ( cot_2 - cot_1 ) )
        csqt_1 = csq_1 * tran_1
        do j = 1, nfinelayers
          tranj(j) = dexp ( - kn * ( cot_2 - cot_fine(n,j) ) )
        enddo

!  Line of sight transmittance factor

        lostrans(n)    = tran_1

!  Elastic scattering multiplier

        func_2 = attn(n-1,w_excit) *  csq_2
        func_1 = attn(n,w_excit)   * csqt_1
        sum = 0.5d0 * ( func_1 + func_2 )
        do j = 1, nfinelayers
          func = attn_fine(n,j,w_excit) * tranj(j) * csq_fine(n,j)
          sum  = sum + func
        enddo
        emultipliers(n) = sum * skn

!  Inelastic scattering multipliers (by bins)

        do w = 1, npoints_mono
          func_2 = attn(n-1,w) *  csq_2
          func_1 = attn(n,w)   * csqt_1
          sum = 0.5d0 * ( func_1 + func_2 )
          do j = 1, nfinelayers
            func = attn_fine(n,j,w) * tranj(j) * csq_fine(n,j)
            sum  = sum + func
          enddo
          imultipliers(n,w) = sum * skn
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
        ( MAXLAYERS, MAXFINELAYERS, MAXPOINTS, &
          NPOINTS_MONO, W_EXCIT, &
          NLAYERS, NFINELAYERS, EXTINCTION, &
          RADII, ALPHA_ALL, ALPHA_FINE, ATTN, ATTN_FINE, &
          EMULTIPLIERS, IMULTIPLIERS, LOSTRANS )

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

       IMPLICIT NONE

!  inputs
!  ------

!  dimensioning

      INTEGER, INTENT(IN) ::          maxlayers, maxfinelayers
      INTEGER, INTENT(IN) ::          maxpoints

!  control

      INTEGER, INTENT(IN) ::          nfinelayers, nlayers
      INTEGER, INTENT(IN) ::          npoints_mono, w_excit

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

      INTEGER ::          w, n, j
      REAL(FPK) :: sum, salpha, calpha, dfine1, raycon, kn
      REAL(FPK) :: csq_1, csqt_1, csq_2, cot_1, cot_2
      REAL(FPK) :: func_1, func_2, tran_1, step, func, skn

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in VLIDORT)

      REAL(FPK), PARAMETER   :: LOCAL_CUTOFF = 32.0D0

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        lostrans(n)     = 0.0d0
        emultipliers(n) = 0.0d0
        do w = 1, npoints_mono
          imultipliers(n,w) = 0.0d0
        enddo
      enddo

!  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  Work Down from the top of the atmosphere
!  ========================================

!  initialise

      n = 0
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

!  Save some quantities

      do n = 1, nlayers
        do j = nfinelayers, 1, -1
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(n,j) = calpha / salpha
          csq_fine(n,j) = 1.0d0 / salpha / salpha
        enddo
      enddo

!  Start layer loop

      do n = 1, nlayers

!  Save some quantities

        salpha = dsin(alpha_all(n))
        calpha = dcos(alpha_all(n))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

!  Trapezium-rule integrated sourceterm multiplier + transmittance
!  ---------------------------------------------------------------

!  set up

        kn = raycon * extinction(n,w_excit)
        skn = step * kn
        tran_1 = dexp ( - kn * ( cot_1 - cot_2 ) )
        csqt_1 = csq_1 * tran_1
        do j = 1, nfinelayers
          tranj(j) = dexp ( - kn * ( cot_fine(n,j) - cot_2 ) )
        enddo

!  Line of sight transmittance factor

        lostrans(n)    = tran_1

!  Elastic multipliers

        func_1 = attn(n-1,w_excit) * csqt_1
        func_2 = attn(n,w_excit)   * csq_2
        sum = 0.5d0 * ( func_1 + func_2 )
        do j = nfinelayers, 1, -1
          func = attn_fine(n,j,w_excit) * tranj(j) * csq_fine(n,j)
          sum = sum + func
        enddo
        emultipliers(n) = sum * step * kn

!  Inelastic multipliers

        do w = 1, npoints_mono
          func_1 = attn(n-1,w) * csqt_1
          func_2 = attn(n,w)   * csq_2
          sum = 0.5d0 * ( func_1 + func_2 )
          do j = nfinelayers, 1, -1
            func = attn_fine(n,j,w) * tranj(j) * csq_fine(n,j)
            sum = sum + func
          enddo
          imultipliers(n,w) = sum * step * kn
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

      SUBROUTINE MULTI_OUTGOING_ADJUSTGEOM &
         ( MAX_VZA, MAX_AZM, N_VZA, N_AZM, &
           HSURFACE, ERADIUS, DO_ADJUST_SURFACE, &
           ALPHA_BOA, THETA_BOA, PHI_BOA, &
           ALPHA_SSA, THETA_SSA, PHI_SSA, FAIL, MESSAGE, TRACE )

! stand-alone geometry routine for adjusting the outgoing correction
!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height of new surface, earth radius.

       IMPLICIT NONE

!  Height grid here is artificial

!  inputs

      INTEGER, INTENT(IN) ::   max_vza, max_azm
      INTEGER, INTENT(IN) ::   n_vza, n_azm
      REAL(FPK), INTENT(IN) :: eradius, hsurface
      LOGICAL, INTENT(IN) ::   do_adjust_surface
      REAL(FPK), INTENT(IN) :: alpha_boa (max_vza)
      REAL(FPK), INTENT(IN) :: theta_boa

!  modified inputs

      REAL(FPK), INTENT(INOUT) :: phi_boa   (max_azm)

!  outputs

      REAL(FPK), INTENT(OUT) :: alpha_ssa (max_vza)
      REAL(FPK), INTENT(OUT) :: theta_ssa (max_vza,max_azm)
      REAL(FPK), INTENT(OUT) :: phi_ssa   (max_vza,max_azm)
      LOGICAL, INTENT(OUT) ::   fail
      CHARACTER (LEN=*), INTENT(OUT) :: message, trace

!  Local

      INTEGER ::   j, k
      REAL(FPK) :: deg_to_rad, ex, ey, ez, px, py, pz, pie, pi2
      REAL(FPK) :: salpha_boa, calpha_boa, sphi_boa
      REAL(FPK) :: stheta_boa, ctheta_boa, cphi_boa
      REAL(FPK) :: ksi, cksi, sksi, xicum, cosscat_up
      REAL(FPK) :: phi_all, alpha_all, theta_all
      REAL(FPK) :: ctheta, stheta, calpha, salpha, cphi
      REAL(FPK) :: b,rssa
      CHARACTER (LEN=2) :: c2

!  Initialise output

      fail = .false.
      message = ' '
      trace   = ' '

!  check range of inputs

      do j = 1, n_vza
       if ( alpha_boa(j).ge.90.0d0.or.alpha_boa(j).lt.0.0d0 ) then
        write(c2,'(I2)')J
        message = 'boa LOS angle outside range [0,90])'
        trace   = 'Change Boa Los angle, number '//C2
        fail    = .true.
        return
       endif
      enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr (1)
!  BUG CORRECTED 10 January 2011------------------------- V. Natraj     (2)

!  IMPORTANT, DO NOT LIMIT THE AZIMUTH TO < 180 degrees

      do k = 1, n_azm
!        if ( phi_boa(k).lt.0.0d0   ) phi_boa(k) = - phi_boa(k)          (2)
        if ( phi_boa(k).lt.0.0d0   ) phi_boa(k) = 360.0d0 + phi_boa(k)
!        if ( phi_boa(k).gt.180.0d0 ) phi_boa(k) = 360.0d0 - phi_boa(k)  (1)
        if ( phi_boa(k).gt.360.0d0 ) phi_boa(k) = phi_boa(k) - 360.0d0
      enddo

!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if ( theta_boa.ge.90.0d0.or.theta_boa.lt.0.0d0 ) then
         message = 'boa SZA angle outside range [0,90])'
         trace   = 'Change Boa SZA angle'
         fail    = .true.
         return
      endif

!  No adjustment, just copy and exit

      if ( .not. do_adjust_surface ) then
         do j = 1, n_vza
            alpha_ssa(j)   = alpha_boa(j)
            do k = 1, n_azm
               theta_ssa(j,k) = theta_boa
               phi_ssa(j,k)   = phi_boa(k)
            enddo
         enddo
         return
      endif

!  conversion

      pie = dacos(-1.0d0)
      pi2 = 2.0d0 * pie
      deg_to_rad = dacos(-1.0d0) / 180.0d0

!  Radius of surface

      rssa = hsurface + eradius

!  Start VZA loop

      do j = 1, n_vza

        alpha_all = alpha_boa(j) * deg_to_rad
        salpha_boa = dsin(alpha_all)
        calpha_boa = dcos(alpha_all)

!  Special case. Direct nadir viewing. Compute everything and Exit.
!    (This is the same as the regular pseudo-spherical )

        if ( salpha_boa.eq.0.0d0 ) then
           alpha_ssa(j)   = alpha_boa(j)
           do k = 1, n_azm
              theta_ssa(j,k) = theta_boa
              phi_ssa(j,k)   = phi_boa(k)
           enddo
           go to 567
        endif

!  Los angle

        salpha       = eradius * salpha_boa / rssa
        alpha_ssa(j) = dasin(salpha)
        calpha       = dcos(alpha_ssa(j))

!  Lospaths

        ksi  = alpha_all - alpha_ssa(j)
        sksi = dsin(ksi)
        cksi = dcos(ksi)
        xicum = ksi

!  output angle in degrees

        alpha_ssa(j) = alpha_ssa(j) / deg_to_rad

!  Los vector

        px = - rssa * dsin(xicum)
        py = 0.0d0
        pz =   rssa * dcos(xicum)

        theta_all = theta_boa * deg_to_rad
        stheta_boa = dsin(theta_all)
        ctheta_boa = dcos(theta_all)

!  Start azimuth loop

        do k = 1, n_azm

          phi_all   = phi_boa(k)   * deg_to_rad
          cphi_boa  = dcos(phi_all)
          sphi_boa  = dsin(phi_all)

!  define Unit solar vector

          ex = - stheta_boa * cphi_boa
          ey = - stheta_boa * sphi_boa
          ez = - ctheta_boa

!  Sun angle

          b = ex*px + ey*py + ez*pz
          ctheta = -b/rssa
          stheta = dsqrt(1.0d0-ctheta*ctheta)
          theta_ssa(j,k) = dacos(ctheta)/deg_to_rad
          if ( ctheta.lt.0.0d0 ) then
            write(c2,'(I2)')J
            message = 'LOS-path SZA angle outside range [0,90])'
            trace   = 'Check inputs for LOS angle '//c2
            fail    = .true.
            return
          endif

!  scattering angle

          cosscat_up  = - calpha_boa * ctheta_boa + &
                          salpha_boa * stheta_boa * cphi_boa

!  Fix phi by using constancy of scatter angle

          if ( phi_boa(k).eq.180.0d0 ) then
            phi_ssa(j,k) = phi_all
          else if ( phi_boa(k) .eq. 0.0d0 ) then
            phi_ssa(j,k) = 0.0d0
          else
            cphi = (cosscat_up+calpha*ctheta)/stheta/salpha
            if ( cphi.gt.1.0d0) cphi = 1.0d0
            if ( cphi.lt.-1.0d0) cphi = -1.0d0
            phi_ssa(j,k) = dacos(cphi)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 05 October 2010----------------- V. Natraj, R. Spurr
!  IF AZIMUTH > 180 degrees. Must ensure consistency, add this line.
            if ( phi_boa(k).gt.180.0d0 ) then
               phi_ssa(j,k)= pi2 - phi_ssa(j,k)
            endif
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

          endif
          phi_ssa(j,k) = phi_ssa(j,k) / deg_to_rad

!  End Azimuth loop

        enddo

!  Continuation point

 567    continue

!  End los loop

      enddo

!  Finish

      return
      END SUBROUTINE MULTI_OUTGOING_ADJUSTGEOM


      END MODULE lrrs_corrections_1

