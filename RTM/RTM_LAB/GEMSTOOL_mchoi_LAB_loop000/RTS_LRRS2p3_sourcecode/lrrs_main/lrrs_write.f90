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
! #   SUBROUTINE :                                              #
! #             LRRS_WRITEINPUT_ALL                             #
! #                                                             #
! ###############################################################

      MODULE lrrs_write

      PRIVATE
      PUBLIC :: LRRS_WRITEINPUT_ALL

      CONTAINS

      SUBROUTINE LRRS_WRITEINPUT_ALL ( &
        DO_RRS_OVERALL,DO_DIRECTRRS_ONLY,DO_DIRECT_BEAM,&
        DO_SSCORR_OUTGOING,DO_SSFULL,DO_MOLECSCAT_ONLY,&
        DO_PLANE_PARALLEL,DO_DELTAM_SCALING,DO_DOUBLE_CONVTEST,&
        DO_UPWELLING,DO_DNWELLING,DO_LBOUNDARIES,DO_MVOUT_ONLY,&
        DO_ADDITIONAL_MVOUT,DO_LAMBERTIAN_SURFACE,DO_WATERLEAVING,&
        DO_SHADOW_EFFECT,DO_GLITTER_DBMS,DO_LAMBERTIAN_FRACTION,&
        DO_LRRS_WRITESCENARIO,DO_LRRS_WRITERESULTS,DO_LRRS_WRITEINPUT,&
        DO_LRRS_WRITEFOURIER,DO_ELASTIC_ONLY,DO_CABANNES_RAMAN,&
        DO_ENERGY_BALANCING,DO_MSMODE_LRRS,DO_BIN_REALIZATION,&
        DO_MONO_REALIZATION,DO_USER_STREAMS,DO_NO_AZIMUTH,&
        PROGRESS,NLAYERS,NSTREAMS,NMOMENTS_INPUT,NFINELAYERS,&
        FLUX_FACTOR,LRRS_ACCURACY,LRRS_FGSMALL,EARTH_RADIUS,&
        SOLAR_ANGLE,DEBUG_FILENAMES,N_USER_STREAMS,USER_ANGLES,&
        N_USER_RELAZMS,USER_RELAZMS,N_LOUTPUT,LBOUNDARIES_OUTPUT,&
        LPARTIALS_OUTPUT,LAMBDAS_RANKED,FLUXES_RANKED,NPOINTS_MONO,&
        LAMBDA_EXCIT,W_EXCIT,NPOINTS_INNER,OFFSET_INNER,NPOINTS_OUTER,&
        BINLOWER,BINUPPER,HEIGHT_GRID,LAYER_TEMPERATURES,&
        LAYER_AIRCOLUMNS,RAYLEIGH_XSEC,RAYLEIGH_DEPOL,&
        DELTAU_INPUT_UNSCALED,OMEGAMOMS_ELASTIC_UNSCALED,&
        LAMBERTIAN_FRACTION,NSTREAMS_BRDF,WHICH_BRDF,BRDF_NPARS,&
        BRDF_PARS,BRDF_FACTOR,BRDF_NAMES,ALBEDOS_RANKED)
        !74 INPUTS

      USE LRRS_PARS

      IMPLICIT NONE

!  Fixed Boolean Inputs

      LOGICAL, INTENT(IN) ::            DO_RRS_OVERALL
      LOGICAL, INTENT(IN) ::            DO_DIRECTRRS_ONLY
      LOGICAL, INTENT(IN) ::            DO_DIRECT_BEAM
      LOGICAL, INTENT(IN) ::            DO_SSCORR_OUTGOING
      LOGICAL, INTENT(IN) ::            DO_SSFULL
      LOGICAL, INTENT(IN) ::            DO_MOLECSCAT_ONLY
      LOGICAL, INTENT(IN) ::            DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN) ::            DO_DELTAM_SCALING
      LOGICAL, INTENT(IN) ::            DO_DOUBLE_CONVTEST
      LOGICAL, INTENT(IN) ::            DO_UPWELLING
      LOGICAL, INTENT(IN) ::            DO_DNWELLING
      LOGICAL, INTENT(IN) ::            DO_LBOUNDARIES
      LOGICAL, INTENT(IN) ::            DO_MVOUT_ONLY
      LOGICAL, INTENT(IN) ::            DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT(IN) ::            DO_LAMBERTIAN_SURFACE
      LOGICAL, INTENT(IN) ::            DO_WATERLEAVING
      LOGICAL, INTENT(IN) ::            DO_SHADOW_EFFECT
      LOGICAL, INTENT(IN) ::            DO_GLITTER_DBMS
      LOGICAL, INTENT(IN) ::            DO_LAMBERTIAN_FRACTION
      LOGICAL, INTENT(IN) ::            DO_LRRS_WRITESCENARIO
      LOGICAL, INTENT(IN) ::            DO_LRRS_WRITERESULTS
      LOGICAL, INTENT(IN) ::            DO_LRRS_WRITEINPUT
      LOGICAL, INTENT(IN) ::            DO_LRRS_WRITEFOURIER

!  Modified Boolean Inputs

      LOGICAL, INTENT(IN) ::            DO_ELASTIC_ONLY
      LOGICAL, INTENT(IN) ::            DO_CABANNES_RAMAN
      LOGICAL, INTENT(IN) ::            DO_ENERGY_BALANCING
      LOGICAL, INTENT(IN) ::            DO_MSMODE_LRRS
      LOGICAL, INTENT(IN) ::            DO_BIN_REALIZATION
      LOGICAL, INTENT(IN) ::            DO_MONO_REALIZATION
      LOGICAL, INTENT(IN) ::            DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::            DO_NO_AZIMUTH

!  Fixed Control Inputs

      INTEGER, INTENT(IN) ::            PROGRESS
      INTEGER, INTENT(IN) ::            NLAYERS
      INTEGER, INTENT(IN) ::            NSTREAMS
      INTEGER, INTENT(IN) ::            NMOMENTS_INPUT
      INTEGER, INTENT(IN) ::            NFINELAYERS
      REAL(FPK), INTENT(IN) ::          FLUX_FACTOR
      REAL(FPK), INTENT(IN) ::          LRRS_ACCURACY
      REAL(FPK), INTENT(IN) ::          LRRS_FGSMALL
      REAL(FPK), INTENT(IN) ::          EARTH_RADIUS
      REAL(FPK), INTENT(IN) ::          SOLAR_ANGLE
      CHARACTER (LEN=60), INTENT(IN) :: DEBUG_FILENAMES(4)

!  Fixed User-Value Inputs

      INTEGER, INTENT(IN) ::            N_USER_STREAMS
      REAL(FPK), INTENT(IN) ::          USER_ANGLES ( MAX_USER_STREAMS )
      INTEGER, INTENT(IN) ::            N_USER_RELAZMS
      REAL(FPK), INTENT(IN) ::          USER_RELAZMS ( MAX_USER_RELAZMS )
      INTEGER, INTENT(IN) ::            N_LOUTPUT
      INTEGER, INTENT(INOUT) ::         LBOUNDARIES_OUTPUT ( MAX_LOUTPUT )
      REAL(FPK), INTENT(IN) ::          LPARTIALS_OUTPUT ( MAX_LOUTPUT )

!  Fixed Spectral Inputs

      REAL(FPK), INTENT(IN) ::          LAMBDAS_RANKED ( MAX_POINTS )
      REAL(FPK), INTENT(IN) ::          FLUXES_RANKED ( MAX_POINTS )

      !For monochromatic calculations:
      INTEGER, INTENT(IN) ::            NPOINTS_MONO
      REAL(FPK), INTENT(IN) ::          LAMBDA_EXCIT
      INTEGER, INTENT(IN) ::            W_EXCIT

      !For binning calculations:
      INTEGER, INTENT(IN) ::            NPOINTS_INNER
      INTEGER, INTENT(IN) ::            OFFSET_INNER
      INTEGER, INTENT(IN) ::            NPOINTS_OUTER
      REAL(FPK), INTENT(IN) ::          BINLOWER ( MAX_POINTS )
      REAL(FPK), INTENT(IN) ::          BINUPPER ( MAX_POINTS )

!  Fixed Atmosphere Inputs

      REAL(FPK), INTENT(IN) ::          HEIGHT_GRID ( 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) ::          LAYER_TEMPERATURES ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) ::          LAYER_AIRCOLUMNS ( MAX_LAYERS )
      REAL(FPK), INTENT(IN) ::          RAYLEIGH_XSEC ( MAX_POINTS )
      REAL(FPK), INTENT(IN) ::          RAYLEIGH_DEPOL ( MAX_POINTS )

      REAL(FPK), INTENT(IN) ::          DELTAU_INPUT_UNSCALED &
                                          ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(IN) ::          OMEGAMOMS_ELASTIC_UNSCALED &
                                          ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, &
                                            MAX_POINTS )

!  Fixed Surface Inputs

      REAL(FPK), INTENT(IN) ::          LAMBERTIAN_FRACTION
      INTEGER, INTENT(IN) ::            NSTREAMS_BRDF
      INTEGER, INTENT(IN) ::            WHICH_BRDF
      INTEGER, INTENT(IN) ::            BRDF_NPARS
      REAL(FPK), INTENT(IN) ::          BRDF_PARS ( MAX_BRDF_PARAMETERS)
      REAL(FPK), INTENT(IN) ::          BRDF_FACTOR
      CHARACTER (LEN=10), INTENT(IN) :: BRDF_NAMES
      REAL(FPK), INTENT(IN) ::          ALBEDOS_RANKED ( MAX_POINTS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: I,LAY,LEV,MOM,PAR,PT,ULEV,URA,UVA

      INTEGER :: LRRS_CALL=0
      CHARACTER (LEN=2)  :: LC_CHAR
      CHARACTER (LEN=26) :: FILENAME

!  Open output file

      OUTUNIT = 100
      LRRS_CALL = LRRS_CALL + 1
      WRITE(LC_CHAR,'(I2.2)') LRRS_CALL
      FILENAME = 'LRRS_WRITEINPUT_ALL_' // LC_CHAR // '.txt'
      WRITE(*,*) 'FILENAME = |',FILENAME,'|'
      OPEN (OUTUNIT,file = FILENAME,status = 'unknown')

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) '-----------------------'
      WRITE(OUTUNIT,*) 'Fixed Boolean Inputs'
      WRITE(OUTUNIT,*) '-----------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) 'DO_RRS_OVERALL         = ',DO_RRS_OVERALL
      WRITE(OUTUNIT,*) 'DO_DIRECTRRS_ONLY      = ',DO_DIRECTRRS_ONLY
      WRITE(OUTUNIT,*) 'DO_DIRECT_BEAM         = ',DO_DIRECT_BEAM
      WRITE(OUTUNIT,*) 'DO_SSCORR_OUTGOING     = ',DO_SSCORR_OUTGOING
      WRITE(OUTUNIT,*) 'DO_SSFULL              = ',DO_SSFULL
      WRITE(OUTUNIT,*) 'DO_MOLECSCAT_ONLY      = ',DO_MOLECSCAT_ONLY
      WRITE(OUTUNIT,*) 'DO_PLANE_PARALLEL      = ',DO_PLANE_PARALLEL
      WRITE(OUTUNIT,*) 'DO_DELTAM_SCALING      = ',DO_DELTAM_SCALING
      WRITE(OUTUNIT,*) 'DO_DOUBLE_CONVTEST     = ',DO_DOUBLE_CONVTEST
      WRITE(OUTUNIT,*) 'DO_UPWELLING           = ',DO_UPWELLING
      WRITE(OUTUNIT,*) 'DO_DNWELLING           = ',DO_DNWELLING
      WRITE(OUTUNIT,*) 'DO_LBOUNDARIES         = ',DO_LBOUNDARIES
      WRITE(OUTUNIT,*) 'DO_MVOUT_ONLY          = ',DO_MVOUT_ONLY
      WRITE(OUTUNIT,*) 'DO_ADDITIONAL_MVOUT    = ',DO_ADDITIONAL_MVOUT
      WRITE(OUTUNIT,*) 'DO_LAMBERTIAN_SURFACE  = ',DO_LAMBERTIAN_SURFACE
      WRITE(OUTUNIT,*) 'DO_WATERLEAVING        = ',DO_WATERLEAVING
      WRITE(OUTUNIT,*) 'DO_SHADOW_EFFECT       = ',DO_SHADOW_EFFECT
      WRITE(OUTUNIT,*) 'DO_GLITTER_DBMS        = ',DO_GLITTER_DBMS
      WRITE(OUTUNIT,*) 'DO_LAMBERTIAN_FRACTION = ',DO_LAMBERTIAN_FRACTION
      WRITE(OUTUNIT,*) 'DO_LRRS_WRITESCENARIO  = ',DO_LRRS_WRITESCENARIO
      WRITE(OUTUNIT,*) 'DO_LRRS_WRITERESULTS   = ',DO_LRRS_WRITERESULTS
      WRITE(OUTUNIT,*) 'DO_LRRS_WRITEINPUT     = ',DO_LRRS_WRITEINPUT
      WRITE(OUTUNIT,*) 'DO_LRRS_WRITEFOURIER   = ',DO_LRRS_WRITEFOURIER

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) '-----------------------'
      WRITE(OUTUNIT,*) 'Modified Boolean Inputs'
      WRITE(OUTUNIT,*) '-----------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) 'DO_ELASTIC_ONLY     = ',DO_ELASTIC_ONLY
      WRITE(OUTUNIT,*) 'DO_CABANNES_RAMAN   = ',DO_CABANNES_RAMAN
      WRITE(OUTUNIT,*) 'DO_ENERGY_BALANCING = ',DO_ENERGY_BALANCING
      WRITE(OUTUNIT,*) 'DO_MSMODE_LRRS      = ',DO_MSMODE_LRRS
      WRITE(OUTUNIT,*) 'DO_BIN_REALIZATION  = ',DO_BIN_REALIZATION
      WRITE(OUTUNIT,*) 'DO_MONO_REALIZATION = ',DO_MONO_REALIZATION
      WRITE(OUTUNIT,*) 'DO_USER_STREAMS     = ',DO_USER_STREAMS
      WRITE(OUTUNIT,*) 'DO_NO_AZIMUTH       = ',DO_NO_AZIMUTH

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) '-----------------------'
      WRITE(OUTUNIT,*) 'Fixed Control Inputs'
      WRITE(OUTUNIT,*) '-----------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) 'PROGRESS       = ',PROGRESS
      WRITE(OUTUNIT,*) 'NLAYERS        = ',NLAYERS
      WRITE(OUTUNIT,*) 'NSTREAMS       = ',NSTREAMS
      WRITE(OUTUNIT,*) 'NMOMENTS_INPUT = ',NMOMENTS_INPUT
      WRITE(OUTUNIT,*) 'NFINELAYERS    = ',NFINELAYERS
      WRITE(OUTUNIT,*) 'FLUX_FACTOR    = ',FLUX_FACTOR
      WRITE(OUTUNIT,*) 'LRRS_ACCURACY  = ',LRRS_ACCURACY
      WRITE(OUTUNIT,*) 'LRRS_FGSMALL   = ',LRRS_FGSMALL
      WRITE(OUTUNIT,*) 'EARTH_RADIUS   = ',EARTH_RADIUS
      WRITE(OUTUNIT,*) 'SOLAR_ANGLE    = ',SOLAR_ANGLE
      DO I=1,4
        WRITE(OUTUNIT,*) 'I = ',I,&
          'DEBUG_FILENAMES(I) = ',DEBUG_FILENAMES(I)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) '-----------------------'
      WRITE(OUTUNIT,*) 'Fixed User-Value Inputs'
      WRITE(OUTUNIT,*) '-----------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) 'N_USER_STREAMS  = ',N_USER_STREAMS
      DO UVA=1,MAX_USER_STREAMS
        WRITE(OUTUNIT,*)  'UVA = ',UVA,&
          ' USER_ANGLES(UVA)  = ',USER_ANGLES(UVA)
      END DO
      WRITE(OUTUNIT,*) 'N_USER_RELAZMS  = ',N_USER_RELAZMS
      DO URA=1,MAX_USER_RELAZMS
        WRITE(OUTUNIT,*)  'URA = ',URA,&
          ' USER_RELAZMS(URA) = ',USER_RELAZMS(URA)
      END DO
      WRITE(OUTUNIT,*) 'N_LOUTPUT       = ',N_LOUTPUT
      DO ULEV=1,MAX_LOUTPUT
        WRITE(OUTUNIT,*)  'ULEV = ',ULEV,&
          ' LBOUNDARIES_OUTPUT(ULEV) = ',LBOUNDARIES_OUTPUT(ULEV)
      END DO
      DO ULEV=1,MAX_LOUTPUT
        WRITE(OUTUNIT,*)  'ULEV = ',ULEV,&
          ' LPARTIALS_OUTPUT(ULEV)   = ',LPARTIALS_OUTPUT(ULEV)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) '-----------------------'
      WRITE(OUTUNIT,*) 'Fixed Spectral Inputs'
      WRITE(OUTUNIT,*) '-----------------------'

      WRITE(OUTUNIT,*)
      DO PT=1,MAX_POINTS
        WRITE(OUTUNIT,*)  'PT = ',PT,&
          ' LAMBDAS_RANKED(PT) = ',LAMBDAS_RANKED(PT)
      END DO
      DO PT=1,MAX_POINTS
        WRITE(OUTUNIT,*)  'PT = ',PT,&
          ' FLUXES_RANKED(PT)  = ',FLUXES_RANKED(PT)
      END DO

      IF (DO_MONO_REALIZATION) THEN
        !For monochromatic calculations:
        WRITE(OUTUNIT,*)
        WRITE(OUTUNIT,*) 'NPOINTS_MONO  = ',NPOINTS_MONO
        WRITE(OUTUNIT,*) 'LAMBDA_EXCIT  = ',LAMBDA_EXCIT
        WRITE(OUTUNIT,*) 'W_EXCIT       = ',W_EXCIT
      END IF

      IF (DO_BIN_REALIZATION) THEN
        !For binning calculations:
        WRITE(OUTUNIT,*)
        WRITE(OUTUNIT,*) 'NPOINTS_INNER = ',NPOINTS_INNER
        WRITE(OUTUNIT,*) 'OFFSET_INNER  = ',OFFSET_INNER
        WRITE(OUTUNIT,*) 'NPOINTS_OUTER = ',NPOINTS_OUTER
        DO PT=1,MAX_POINTS
          WRITE(OUTUNIT,*)  'PT = ',PT,&
            ' BINLOWER(PT) = ',BINLOWER(PT)
        END DO
        DO PT=1,MAX_POINTS
          WRITE(OUTUNIT,*)  'PT = ',PT,&
            ' BINUPPER(PT) = ',BINUPPER(PT)
        END DO
      END IF

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) '-----------------------'
      WRITE(OUTUNIT,*) 'Fixed Atmosphere Inputs'
      WRITE(OUTUNIT,*) '-----------------------'

      WRITE(OUTUNIT,*)
      DO LAY=0,MAX_LAYERS
        WRITE(OUTUNIT,*)  'LAY = ',LAY,&
          ' HEIGHT_GRID(LAY)        = ',HEIGHT_GRID(LAY)
      END DO
      DO LAY=1,MAX_LAYERS
        WRITE(OUTUNIT,*)  'LAY = ',LAY,&
          ' LAYER_TEMPERATURES(LAY) = ',LAYER_TEMPERATURES(LAY)
      END DO
      DO LAY=1,MAX_LAYERS
        WRITE(OUTUNIT,*)  'LAY = ',LAY,&
          ' LAYER_AIRCOLUMNS(LAY)   = ',LAYER_AIRCOLUMNS(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,MAX_POINTS
        WRITE(OUTUNIT,*)  'PT = ',PT,&
          ' RAYLEIGH_XSEC(PT)  = ',RAYLEIGH_XSEC(PT)
      END DO
      DO PT=1,MAX_POINTS
        WRITE(OUTUNIT,*)  'PT = ',PT,&
          ' RAYLEIGH_DEPOL(PT) = ',RAYLEIGH_DEPOL(PT)
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,MAX_POINTS
        DO LAY=1,MAX_LAYERS
          WRITE(OUTUNIT,*)  'PT = ',PT,' LAY = ',LAY,&
            ' DELTAU_INPUT_UNSCALED(LAY,PT) = ',&
              DELTAU_INPUT_UNSCALED(LAY,PT)
        END DO
      END DO
      !THE DIMENSIONS OF OMEGAMOMS_ELASTIC_UNSCALED ARE
      !"IN THE WRONG ORDER"??
      DO PT=1,MAX_POINTS
        DO LAY=1,MAX_LAYERS
          DO MOM=0,MAX_MOMENTS_INPUT
            WRITE(OUTUNIT,*)  'PT = ',PT,' LAY = ',LAY,' MOM = ',MOM,&
              ' OMEGAMOMS_ELASTIC_UNSCALED(LAY,MOM,PT) = ',&
                OMEGAMOMS_ELASTIC_UNSCALED(LAY,MOM,PT)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) '-----------------------'
      WRITE(OUTUNIT,*) 'Fixed Surface Inputs'
      WRITE(OUTUNIT,*) '-----------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,*) 'LAMBERTIAN_FRACTION = ',LAMBERTIAN_FRACTION
      WRITE(OUTUNIT,*) 'NSTREAMS_BRDF       = ',NSTREAMS_BRDF
      WRITE(OUTUNIT,*) 'WHICH_BRDF          = ',WHICH_BRDF
      WRITE(OUTUNIT,*) 'BRDF_NPARS          = ',BRDF_NPARS
      DO PAR=1,MAX_BRDF_PARAMETERS
        WRITE(OUTUNIT,*)  'PAR = ',PAR,&
          ' BRDF_PARS(PAR)     = ',BRDF_PARS(PAR)
      END DO
      WRITE(OUTUNIT,*) 'BRDF_NAMES          = ',BRDF_NAMES
      DO PT=1,MAX_POINTS
        WRITE(OUTUNIT,*)  'PT = ',PT,&
          ' ALBEDOS_RANKED(PT) = ',ALBEDOS_RANKED(PT)
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LRRS_WRITEINPUT_ALL

      END MODULE lrrs_write
