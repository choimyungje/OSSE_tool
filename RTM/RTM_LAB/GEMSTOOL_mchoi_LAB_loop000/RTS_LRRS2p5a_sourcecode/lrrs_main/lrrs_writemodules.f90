
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
! #              LRRS_WRITE_STD_INPUT                           #
! #              LRRS_WRITE_SUP_BRDF_INPUT                      #
! #              LRRS_WRITE_SUP_SLEAVE_INPUT                    #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5.  Changes to this routine are---
!    (1) Write out additional Taylor-seris and BRDF/SLEAVE inputs.
!    (2) Add BRDF & SLEAVE subroutines

      MODULE lrrs_writemodules_m

      USE LRRS_PARS_m, Only : SDU

      PRIVATE
      PUBLIC :: LRRS_WRITE_STD_INPUT, &
                LRRS_WRITE_SUP_BRDF_INPUT, &
                LRRS_WRITE_SUP_SLEAVE_INPUT

      CONTAINS

      SUBROUTINE LRRS_WRITE_STD_INPUT ( &
        DO_RRS_OVERALL,DO_DIRECTRRS_ONLY,DO_DIRECT_BEAM,DO_SSCORR_OUTGOING,DO_SSCORR_NADIR, &
        DO_SSCORR_ALONE,DO_MOLECSCAT_ONLY,DO_PLANE_PARALLEL,DO_DELTAM_SCALING,              &
        DO_DOUBLE_CONVTEST,DO_UPWELLING,DO_DNWELLING,DO_LBOUNDARIES,           &
        DO_MVOUT_ONLY,DO_ADDITIONAL_MVOUT,DO_BRDF_SURFACE,DO_SURFACE_LEAVING,  &
        DO_SL_ISOTROPIC,DO_LRRS_WRITEINPUT,DO_LRRS_WRITESCENARIO,              &
        DO_LRRS_WRITEFOURIER,DO_LRRS_WRITERESULTS,                             &
        PROGRESS,NLAYERS,NSTREAMS,NFINELAYERS,FLUX_FACTOR,                     &
        SOLAR_ANGLE,EARTH_RADIUS,LRRS_ACCURACY,TAYLOR_ORDER,TAYLOR_SMALL,      &
        DEBUG_FILENAMES,N_USER_STREAMS,USER_ANGLES,N_USER_RELAZMS,USER_RELAZMS,&
        N_LOUTPUT,LPARTIALS_OUTPUT,LBOUNDARIES_OUTPUT,LAMBDAS_RANKED,          &
        FLUXES_RANKED,NPOINTS_MONO,LAMBDA_EXCIT,W_EXCIT,NPOINTS_INNER,         &
        OFFSET_INNER,NPOINTS_OUTER,BINLOWER,BINUPPER,HEIGHT_GRID,              &
        LAYER_TEMPERATURES,LAYER_AIRCOLUMNS,RAYLEIGH_XSEC,RAYLEIGH_DEPOL,      &
        DELTAU_INPUT_UNSCALED,OMEGAMOMS_ELASTIC_UNSCALED,ALBEDOS_RANKED,       &
        OMEGAPHASFUNC_ELASTIC_UP, OMEGAPHASFUNC_ELASTIC_DN,                    &
        DO_ELASTIC_ONLY,DO_CABANNES_RAMAN,DO_ENERGY_BALANCING,DO_MSMODE_LRRS,  &
        DO_BIN_REALIZATION,DO_MONO_REALIZATION,DO_USER_STREAMS,DO_NO_AZIMUTH,  &
        GEOMETRY_SPECHEIGHT,NPOINTS )

!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT, Added PHASFUNC inputs
!   -- Rob mod 5/12/17 for 2p5a, renamed SSFULL --> SSCORR_ALONE, Added SSCORR_NADIR

      USE LRRS_PARS_m

      IMPLICIT NONE

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

!  Fixed Boolean Inputs

      LOGICAL  , INTENT(IN) ::          DO_RRS_OVERALL
      LOGICAL  , INTENT(IN) ::          DO_DIRECTRRS_ONLY
      LOGICAL  , INTENT(IN) ::          DO_DIRECT_BEAM
      LOGICAL  , INTENT(IN) ::          DO_SSCORR_NADIR
      LOGICAL  , INTENT(IN) ::          DO_SSCORR_OUTGOING
      LOGICAL  , INTENT(IN) ::          DO_SSCORR_ALONE
      LOGICAL  , INTENT(IN) ::          DO_MOLECSCAT_ONLY
      LOGICAL  , INTENT(IN) ::          DO_PLANE_PARALLEL
      LOGICAL  , INTENT(IN) ::          DO_DELTAM_SCALING
      LOGICAL  , INTENT(IN) ::          DO_DOUBLE_CONVTEST
      LOGICAL  , INTENT(IN) ::          DO_UPWELLING
      LOGICAL  , INTENT(IN) ::          DO_DNWELLING
      LOGICAL  , INTENT(IN) ::          DO_LBOUNDARIES
      LOGICAL  , INTENT(IN) ::          DO_MVOUT_ONLY
      LOGICAL  , INTENT(IN) ::          DO_ADDITIONAL_MVOUT

!  Next 5 commented out (formerly Version 2.3)
!      LOGICAL  , INTENT(IN) ::            DO_LAMBERTIAN_SURFACE
!      LOGICAL  , INTENT(IN) ::            DO_WATERLEAVING
!      LOGICAL  , INTENT(IN) ::            DO_SHADOW_EFFECT
!      LOGICAL  , INTENT(IN) ::            DO_GLITTER_DBMS
!      LOGICAL  , INTENT(IN) ::            DO_LAMBERTIAN_FRACTION

!  Next 3 are new for Version 2.5
      LOGICAL  , INTENT(IN) ::          DO_BRDF_SURFACE
      LOGICAL  , INTENT(IN) ::          DO_SURFACE_LEAVING
      LOGICAL  , INTENT(IN) ::          DO_SL_ISOTROPIC

      LOGICAL  , INTENT(IN) ::          DO_LRRS_WRITEINPUT
      LOGICAL  , INTENT(IN) ::          DO_LRRS_WRITESCENARIO
      LOGICAL  , INTENT(IN) ::          DO_LRRS_WRITEFOURIER
      LOGICAL  , INTENT(IN) ::          DO_LRRS_WRITERESULTS

!  Fixed Control Inputs. Taylor -series inputs have changed, version 2.5, 9/11/15
!      INTEGER  , INTENT(IN) ::          NMOMENTS_INPUT

      INTEGER  , INTENT(IN) ::          PROGRESS
      INTEGER  , INTENT(IN) ::          NLAYERS
      INTEGER  , INTENT(IN) ::          NSTREAMS
      INTEGER  , INTENT(IN) ::          NFINELAYERS
      REAL(FPK), INTENT(IN) ::          FLUX_FACTOR
      REAL(FPK), INTENT(IN) ::          SOLAR_ANGLE
      REAL(FPK), INTENT(IN) ::          EARTH_RADIUS
      REAL(FPK), INTENT(IN) ::          LRRS_ACCURACY
      INTEGER  , INTENT(IN) ::          TAYLOR_ORDER
      REAL(FPK), INTENT(IN) ::          TAYLOR_SMALL
      CHARACTER (LEN=60), INTENT(IN) :: DEBUG_FILENAMES(4)

!  Fixed User-Value Inputs

      INTEGER  , INTENT(IN) ::          N_USER_STREAMS
      REAL(FPK), INTENT(IN) ::          USER_ANGLES ( MAX_USER_STREAMS )
      INTEGER  , INTENT(IN) ::          N_USER_RELAZMS
      REAL(FPK), INTENT(IN) ::          USER_RELAZMS ( MAX_USER_RELAZMS )
      INTEGER  , INTENT(IN) ::          N_LOUTPUT
      REAL(FPK), INTENT(IN) ::          LPARTIALS_OUTPUT ( MAX_LOUTPUT ) !not active yet (7/20/2016)
      INTEGER  , INTENT(INOUT) ::       LBOUNDARIES_OUTPUT ( MAX_LOUTPUT )

!  Fixed Spectral Inputs

      REAL(FPK), INTENT(IN) ::          LAMBDAS_RANKED ( MAX_POINTS )
      REAL(FPK), INTENT(IN) ::          FLUXES_RANKED ( MAX_POINTS )

      !For monochromatic calculations:
      INTEGER  , INTENT(IN) ::          NPOINTS_MONO
      REAL(FPK), INTENT(IN) ::          LAMBDA_EXCIT
      INTEGER  , INTENT(IN) ::          W_EXCIT

      !For binning calculations:
      INTEGER  , INTENT(IN) ::          NPOINTS_INNER
      INTEGER  , INTENT(IN) ::          OFFSET_INNER
      INTEGER  , INTENT(IN) ::          NPOINTS_OUTER
      REAL(FPK), INTENT(IN) ::          BINLOWER ( MAX_POINTS )
      REAL(FPK), INTENT(IN) ::          BINUPPER ( MAX_POINTS )

!  Fixed Atmosphere Inputs

      REAL(FPK), INTENT(IN) ::          HEIGHT_GRID ( 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) ::          LAYER_TEMPERATURES ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) ::          LAYER_AIRCOLUMNS ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) ::          RAYLEIGH_XSEC ( MAX_POINTS )
      REAL(FPK), INTENT(IN) ::          RAYLEIGH_DEPOL ( MAX_POINTS )

!  -- Rob mod 5/12/17 for 2p5a, changed moment dimension in OMEGAMOMS_ELASTIC_UNSCALED

      REAL(FPK), INTENT(IN) ::          DELTAU_INPUT_UNSCALED ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) ::          OMEGAMOMS_ELASTIC_UNSCALED ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  -- Rob mod 5/12/17 for 2p5a, added these inputs for Write-up

      REAL(FPK), INTENT(IN) ::          OMEGAPHASFUNC_ELASTIC_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) ::          OMEGAPHASFUNC_ELASTIC_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Fixed Surface Inputs. Only the albedo survives, Version 2.5, 9/11/15

!      REAL(FPK), INTENT(IN) ::          LAMBERTIAN_FRACTION
!      INTEGER  , INTENT(IN) ::            NSTREAMS_BRDF
!      INTEGER  , INTENT(IN) ::            WHICH_BRDF
!      INTEGER  , INTENT(IN) ::            BRDF_NPARS
!      REAL(FPK), INTENT(IN) ::          BRDF_PARS ( MAX_BRDF_PARAMETERS)
!      REAL(FPK), INTENT(IN) ::          BRDF_FACTOR
!      CHARACTER (LEN=10), INTENT(IN) :: BRDF_NAMES

      REAL(FPK), INTENT(IN) ::          ALBEDOS_RANKED ( MAX_POINTS )

!  --------------------------
!  Standard Inputs - Variable
!  --------------------------

!  Modified Boolean Inputs

      LOGICAL  , INTENT(IN) ::          DO_ELASTIC_ONLY
      LOGICAL  , INTENT(IN) ::          DO_CABANNES_RAMAN
      LOGICAL  , INTENT(IN) ::          DO_ENERGY_BALANCING
      LOGICAL  , INTENT(IN) ::          DO_MSMODE_LRRS
      LOGICAL  , INTENT(IN) ::          DO_BIN_REALIZATION
      LOGICAL  , INTENT(IN) ::          DO_MONO_REALIZATION
      LOGICAL  , INTENT(IN) ::          DO_USER_STREAMS
      LOGICAL  , INTENT(IN) ::          DO_NO_AZIMUTH

!  Modified Atmosphere Inputs

      REAL(FPK), INTENT(IN) ::          GEOMETRY_SPECHEIGHT

!  LRRS derived input

      INTEGER  , INTENT(IN) ::          NPOINTS

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: I,LAY,MOM,PT,ULEV,URA,UVA,NG,V

!  Open output file

      OUTUNIT = 101
      OPEN (OUTUNIT,file = 'LRRS_WRITE_STD_INPUT.dbg',status = 'replace')

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Fixed'
      WRITE(OUTUNIT,'(A)') '-----------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL) 'DO_RRS_OVERALL         = ',DO_RRS_OVERALL
      WRITE(OUTUNIT,DWFL) 'DO_DIRECTRRS_ONLY      = ',DO_DIRECTRRS_ONLY
      WRITE(OUTUNIT,DWFL) 'DO_DIRECT_BEAM         = ',DO_DIRECT_BEAM
      WRITE(OUTUNIT,DWFL) 'DO_SSCORR_NADIR        = ',DO_SSCORR_NADIR
      WRITE(OUTUNIT,DWFL) 'DO_SSCORR_OUTGOING     = ',DO_SSCORR_OUTGOING
      WRITE(OUTUNIT,DWFL) 'DO_SSCORR_ALONE        = ',DO_SSCORR_ALONE
      WRITE(OUTUNIT,DWFL) 'DO_MOLECSCAT_ONLY      = ',DO_MOLECSCAT_ONLY
      WRITE(OUTUNIT,DWFL) 'DO_PLANE_PARALLEL      = ',DO_PLANE_PARALLEL
      WRITE(OUTUNIT,DWFL) 'DO_DELTAM_SCALING      = ',DO_DELTAM_SCALING
      WRITE(OUTUNIT,DWFL) 'DO_DOUBLE_CONVTEST     = ',DO_DOUBLE_CONVTEST
      WRITE(OUTUNIT,DWFL) 'DO_UPWELLING           = ',DO_UPWELLING
      WRITE(OUTUNIT,DWFL) 'DO_DNWELLING           = ',DO_DNWELLING
      WRITE(OUTUNIT,DWFL) 'DO_LBOUNDARIES         = ',DO_LBOUNDARIES
      WRITE(OUTUNIT,DWFL) 'DO_MVOUT_ONLY          = ',DO_MVOUT_ONLY
      WRITE(OUTUNIT,DWFL) 'DO_ADDITIONAL_MVOUT    = ',DO_ADDITIONAL_MVOUT

      WRITE(OUTUNIT,DWFL) 'DO_BRDF_SURFACE        = ',DO_BRDF_SURFACE
      WRITE(OUTUNIT,DWFL) 'DO_SURFACE_LEAVING     = ',DO_SURFACE_LEAVING
      WRITE(OUTUNIT,DWFL) 'DO_SL_ISOTROPIC        = ',DO_SL_ISOTROPIC

!      WRITE(OUTUNIT,DWFL) 'DO_LAMBERTIAN_SURFACE  = ',DO_LAMBERTIAN_SURFACE
!      WRITE(OUTUNIT,DWFL) 'DO_WATERLEAVING        = ',DO_WATERLEAVING
!      WRITE(OUTUNIT,DWFL) 'DO_SHADOW_EFFECT       = ',DO_SHADOW_EFFECT
!      WRITE(OUTUNIT,DWFL) 'DO_GLITTER_DBMS        = ',DO_GLITTER_DBMS
!      WRITE(OUTUNIT,DWFL) 'DO_LAMBERTIAN_FRACTION = ',DO_LAMBERTIAN_FRACTION

      WRITE(OUTUNIT,DWFL) 'DO_LRRS_WRITESCENARIO  = ',DO_LRRS_WRITESCENARIO
      WRITE(OUTUNIT,DWFL) 'DO_LRRS_WRITERESULTS   = ',DO_LRRS_WRITERESULTS
      WRITE(OUTUNIT,DWFL) 'DO_LRRS_WRITEINPUT     = ',DO_LRRS_WRITEINPUT
      WRITE(OUTUNIT,DWFL) 'DO_LRRS_WRITEFOURIER   = ',DO_LRRS_WRITEFOURIER

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI) 'PROGRESS          = ',PROGRESS
      WRITE(OUTUNIT,DWFI) 'NLAYERS           = ',NLAYERS
      WRITE(OUTUNIT,DWFI) 'NSTREAMS          = ',NSTREAMS
!      WRITE(OUTUNIT,DWFI) 'NMOMENTS_INPUT    = ',NMOMENTS_INPUT   -- Rob mod 5/2/17 for 2p5a, removed
      WRITE(OUTUNIT,DWFI) 'NFINELAYERS       = ',NFINELAYERS
      WRITE(OUTUNIT,DWFR) 'FLUX_FACTOR       = ',FLUX_FACTOR
      WRITE(OUTUNIT,DWFR) 'LRRS_ACCURACY     = ',LRRS_ACCURACY
      WRITE(OUTUNIT,DWFI) 'TAYLOR_ORDER      = ',TAYLOR_ORDER
      WRITE(OUTUNIT,DWFR) 'TAYLOR_SMALL      = ',TAYLOR_SMALL
      WRITE(OUTUNIT,DWFR) 'EARTH_RADIUS      = ',EARTH_RADIUS
      WRITE(OUTUNIT,DWFR) 'SOLAR_ANGLE       = ',SOLAR_ANGLE
      WRITE(OUTUNIT,*)
      DO I=1,4
        WRITE(OUTUNIT,DWFC1) 'I = ',I,&
          ' DEBUG_FILENAMES(I) = ',DEBUG_FILENAMES(I)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI) 'N_USER_STREAMS  = ',N_USER_STREAMS
      DO UVA=1,N_USER_STREAMS
        WRITE(OUTUNIT,DWFR1)  'UVA = ',UVA,&
          ' USER_ANGLES(UVA)  = ',USER_ANGLES(UVA)
      END DO
      WRITE(OUTUNIT,DWFI) 'N_USER_RELAZMS  = ',N_USER_RELAZMS
      DO URA=1,N_USER_RELAZMS
        WRITE(OUTUNIT,DWFR1)  'URA = ',URA,&
          ' USER_RELAZMS(URA) = ',USER_RELAZMS(URA)
      END DO
      WRITE(OUTUNIT,DWFI) 'N_LOUTPUT       = ',N_LOUTPUT
      DO ULEV=1,N_LOUTPUT
        WRITE(OUTUNIT,DWFR1)  'ULEV = ',ULEV,&
          ' LPARTIALS_OUTPUT(ULEV)   = ',LPARTIALS_OUTPUT(ULEV)
      END DO
      DO ULEV=1,N_LOUTPUT
        WRITE(OUTUNIT,DWFI1)  'ULEV = ',ULEV,&
          ' LBOUNDARIES_OUTPUT(ULEV) = ',LBOUNDARIES_OUTPUT(ULEV)
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        WRITE(OUTUNIT,DWFR1)  'PT = ',PT,&
          ' LAMBDAS_RANKED(PT) = ',LAMBDAS_RANKED(PT)
      END DO
      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        WRITE(OUTUNIT,DWFR1)  'PT = ',PT,&
          ' FLUXES_RANKED(PT)  = ',FLUXES_RANKED(PT)
      END DO

      IF (DO_MONO_REALIZATION) THEN
        !For monochromatic calculations:
        WRITE(OUTUNIT,*)
        WRITE(OUTUNIT,DWFI) 'NPOINTS_MONO  = ',NPOINTS_MONO
        WRITE(OUTUNIT,DWFR) 'LAMBDA_EXCIT  = ',LAMBDA_EXCIT
        WRITE(OUTUNIT,DWFI) 'W_EXCIT       = ',W_EXCIT
      END IF

      IF (DO_BIN_REALIZATION) THEN
        !For binning calculations:
        WRITE(OUTUNIT,*)
        WRITE(OUTUNIT,DWFI) 'NPOINTS_INNER = ',NPOINTS_INNER
        WRITE(OUTUNIT,DWFI) 'OFFSET_INNER  = ',OFFSET_INNER
        WRITE(OUTUNIT,DWFI) 'NPOINTS_OUTER = ',NPOINTS_OUTER
        WRITE(OUTUNIT,*)
        DO PT=1,NPOINTS
          WRITE(OUTUNIT,DWFR1)  'PT = ',PT,&
            ' BINLOWER(PT) = ',BINLOWER(PT)
        END DO
        WRITE(OUTUNIT,*)
        DO PT=1,NPOINTS
          WRITE(OUTUNIT,DWFR1)  'PT = ',PT,&
            ' BINUPPER(PT) = ',BINUPPER(PT)
        END DO
      END IF

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' HEIGHT_GRID(LAY)        = ',HEIGHT_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' LAYER_TEMPERATURES(LAY) = ',LAYER_TEMPERATURES(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' LAYER_AIRCOLUMNS(LAY)   = ',LAYER_AIRCOLUMNS(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        WRITE(OUTUNIT,DWFR1)  'PT = ',PT,&
          ' RAYLEIGH_XSEC(PT)  = ',RAYLEIGH_XSEC(PT)
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        WRITE(OUTUNIT,DWFR1)  'PT = ',PT,&
          ' RAYLEIGH_DEPOL(PT) = ',RAYLEIGH_DEPOL(PT)
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO LAY=1,NLAYERS
          WRITE(OUTUNIT,DWFR2)  'PT = ',PT,' LAY = ',LAY,&
            ' DELTAU_INPUT_UNSCALED(LAY,PT) = ',&
              DELTAU_INPUT_UNSCALED(LAY,PT)
        END DO
      END DO

      !THE 1ST & 2ND DIMENSIONS OF OMEGAMOMS_ELASTIC_UNSCALED ARE
      !"IN THE WRONG ORDER"??
!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT, output all moments to MAX_MOMENTS instead

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO LAY=1,NLAYERS
          DO MOM=0,MAX_MOMENTS      !  Rob mod 5/12/17 for 2p5a
            WRITE(OUTUNIT,DWFR3)  'PT = ',PT,' LAY = ',LAY,' MOM = ',MOM,&
              ' OMEGAMOMS_ELASTIC_UNSCALED(LAY,MOM,PT) = ',&
                OMEGAMOMS_ELASTIC_UNSCALED(LAY,MOM,PT)
          END DO
        END DO
      END DO

!  Write the OMEGAPHASFUNC inputs
!   -- Rob mod 5/12/17 for 2p5a

      NG = N_USER_STREAMS * N_USER_RELAZMS
      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO LAY=1,NLAYERS
          DO V=1,NG
            WRITE(OUTUNIT,DWFR3)  'PT = ',PT,' LAY = ',LAY,' GEOM = ',V,&
              ' OMEGAPHASFUNC_ELASTIC_UP(LAY,GEOM,PT) = ',&
                OMEGAPHASFUNC_ELASTIC_UP(LAY,V,PT) 
          END DO
        END DO
      END DO
      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO LAY=1,NLAYERS
          DO V=1,NG
            WRITE(OUTUNIT,DWFR3)  'PT = ',PT,' LAY = ',LAY,' GEOM = ',V,&
              ' OMEGAPHASFUNC_ELASTIC_DN(LAY,GEOM,PT) = ',&
                OMEGAPHASFUNC_ELASTIC_DN(LAY,V,PT) 
          END DO
        END DO
      END DO

!  Only the albedo survives, Version 2.5.

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        WRITE(OUTUNIT,DWFR1)  'PT = ',PT,&
          ' ALBEDOS_RANKED(PT) = ',ALBEDOS_RANKED(PT)
      END DO

!      WRITE(OUTUNIT,*) 'LAMBERTIAN_FRACTION = ',LAMBERTIAN_FRACTION
!      WRITE(OUTUNIT,*) 'NSTREAMS_BRDF       = ',NSTREAMS_BRDF
!      WRITE(OUTUNIT,*) 'WHICH_BRDF          = ',WHICH_BRDF
!      WRITE(OUTUNIT,*) 'BRDF_NPARS          = ',BRDF_NPARS
!      DO PAR=1,MAX_BRDF_PARAMETERS
!        WRITE(OUTUNIT,*)  'PAR = ',PAR,&
!          ' BRDF_PARS(PAR)     = ',BRDF_PARS(PAR)
!      END DO
!      WRITE(OUTUNIT,*) 'BRDF_NAMES          = ',BRDF_NAMES

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '--------------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Variable'
      WRITE(OUTUNIT,'(A)') '--------------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL) 'DO_ELASTIC_ONLY     = ',DO_ELASTIC_ONLY
      WRITE(OUTUNIT,DWFL) 'DO_CABANNES_RAMAN   = ',DO_CABANNES_RAMAN
      WRITE(OUTUNIT,DWFL) 'DO_ENERGY_BALANCING = ',DO_ENERGY_BALANCING
      WRITE(OUTUNIT,DWFL) 'DO_MSMODE_LRRS      = ',DO_MSMODE_LRRS
      WRITE(OUTUNIT,DWFL) 'DO_BIN_REALIZATION  = ',DO_BIN_REALIZATION
      WRITE(OUTUNIT,DWFL) 'DO_MONO_REALIZATION = ',DO_MONO_REALIZATION
      WRITE(OUTUNIT,DWFL) 'DO_USER_STREAMS     = ',DO_USER_STREAMS
      WRITE(OUTUNIT,DWFL) 'DO_NO_AZIMUTH       = ',DO_NO_AZIMUTH

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) 'GEOMETRY_SPECHEIGHT = ',GEOMETRY_SPECHEIGHT

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LRRS_WRITE_STD_INPUT

!

      SUBROUTINE LRRS_WRITE_SUP_BRDF_INPUT ( &
        NSTREAMS,N_USER_STREAMS,N_USER_RELAZMS,NPOINTS,&
        EXACTDB_BRDFUNC,BRDF_F_0,BRDF_F,USER_BRDF_F_0,USER_BRDF_F )

      USE LRRS_PARS_m

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   NPOINTS

      REAL(fpk), INTENT(IN) ::   EXACTDB_BRDFUNC &
          ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )
      REAL(fpk), INTENT(IN) ::   BRDF_F_0 &
          ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), INTENT(IN) ::   BRDF_F &
          ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), INTENT(IN) ::   USER_BRDF_F_0 &
          ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )
      REAL(fpk), INTENT(IN) ::   USER_BRDF_F &
          ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: MOM,PT,STRM,STRMI,STRMJ,USTRM,URA
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 102
      OPEN (OUTUNIT,file = 'LRRS_WRITE_SUP_BRDF_INPUT.dbg',&
            status = 'replace')

!  Define local variables

      NMOMENTS = 2*NSTREAMS

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '----------------------'
      WRITE(OUTUNIT,'(A)') 'BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '----------------------'

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_STREAMS
              WRITE(OUTUNIT,DWFR3) &
                'PT = ',PT,' URA = ',URA,' USTRM = ',USTRM,&
                ' EXACTDB_BRDFUNC(USTRM,URA,PT) = ',&
                  EXACTDB_BRDFUNC(USTRM,URA,PT)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO STRM=1,NSTREAMS
          DO MOM=0,NMOMENTS
            WRITE(OUTUNIT,DWFR3) &
              'PT = ',PT,' STRM = ',STRM,' MOM = ',MOM,&
              ' BRDF_F_0(MOM,STRM,PT) = ',BRDF_F_0(MOM,STRM,PT)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO STRMJ=1,NSTREAMS
          DO STRMI=1,NSTREAMS
            DO MOM=0,NMOMENTS
              WRITE(OUTUNIT,DWFR4) &
                'PT = ',PT,' STRMJ = ',STRMJ,' STRMI = ',STRMI,' MOM = ',MOM,&
                ' BRDF_F(MOM,STRMI,STRMJ,PT) = ',BRDF_F(MOM,STRMI,STRMJ,PT)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO USTRM=1,N_USER_STREAMS
          DO MOM=0,NMOMENTS
            WRITE(OUTUNIT,DWFR3) &
              'PT = ',PT,' USTRM = ',USTRM,' MOM = ',MOM,&
              ' USER_BRDF_F_0(MOM,USTRM,PT) = ',&
                USER_BRDF_F_0(MOM,USTRM,PT)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO STRM=1,NSTREAMS
          DO USTRM=1,N_USER_STREAMS
            DO MOM=0,NMOMENTS
              WRITE(OUTUNIT,DWFR4) &
                'PT = ',PT,' STRM = ',STRM,' USTRM = ',USTRM,' MOM = ',MOM,&
                ' USER_BRDF_F(MOM,USTRM,STRM,PT) = ',&
                  USER_BRDF_F(MOM,USTRM,STRM,PT)
            END DO
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LRRS_WRITE_SUP_BRDF_INPUT

!

      SUBROUTINE LRRS_WRITE_SUP_SLEAVE_INPUT ( &
        NSTREAMS,N_USER_STREAMS,N_USER_RELAZMS,NPOINTS,&
        SLTERM_ISOTROPIC,SLTERM_USERANGLES,SLTERM_F_0,USER_SLTERM_F_0)

      USE LRRS_PARS_m

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   NPOINTS

      REAL(fpk), INTENT(IN) ::   SLTERM_ISOTROPIC  ( MAX_POINTS )
      REAL(fpk), INTENT(IN) ::   SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

      REAL(fpk), INTENT(IN) ::   SLTERM_F_0      ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), INTENT(IN) ::   USER_SLTERM_F_0 ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: MOM,PT,STRM,USTRM,URA
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'LRRS_WRITE_SUP_SLEAVE_INPUT.dbg',status = 'replace')

!  Define local variable

      NMOMENTS = 2*NSTREAMS

!  Write all surface-leaving input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '------------------------'
      WRITE(OUTUNIT,'(A)') 'SLEAVE Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '------------------------'

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        WRITE(OUTUNIT,DWFR1) &
          'PT = ',PT,&
          ' SLTERM_ISOTROPIC(PT) = ',SLTERM_ISOTROPIC(PT)
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_STREAMS
            WRITE(OUTUNIT,DWFR3) &
              'PT = ',PT,' URA = ',URA,' USTRM = ',USTRM,&
              ' SLTERM_USERANGLES(USTRM,URA,PT) = ',&
                SLTERM_USERANGLES(USTRM,URA,PT)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO STRM=1,NSTREAMS
          DO MOM=0,NMOMENTS
            WRITE(OUTUNIT,DWFR3) &
              'PT = ',PT,' STRM = ',STRM,' MOM = ',MOM,&
              ' SLTERM_F_0(MOM,STRM,PT) = ',SLTERM_F_0(MOM,STRM,PT)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO USTRM=1,N_USER_STREAMS
          DO MOM=0,NMOMENTS
            WRITE(OUTUNIT,DWFR3) &
              'PT = ',PT,' USTRM = ',USTRM,' MOM = ',MOM,&
              ' USER_SLTERM_F_0(MOM,USTRM,PT) = ',&
                USER_SLTERM_F_0(MOM,USTRM,PT)
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LRRS_WRITE_SUP_SLEAVE_INPUT

      END MODULE lrrs_writemodules_m
