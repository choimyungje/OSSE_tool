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
! #             RAMAN_MASTER (master)                           #
! #                                                             #
! ###############################################################

      MODULE lrrs_main_master

      PRIVATE
      PUBLIC :: RAMAN_MASTER

      CONTAINS

      SUBROUTINE RAMAN_MASTER (LRRS_FixIn, LRRS_ModIn, LRRS_Out)

!  include file of dimensions and numbers

      USE LRRS_PARS
      USE LRRS_INPUTS_DEF
      USE LRRS_OUTPUTS_DEF
      USE LRRS_IO_CHECK
      USE LRRS_SETUP_MASTER
      USE LRRS_RTSOLUTIONS_1
      USE LRRS_BRDF
      USE LRRS_CORRECTIONS_1
      USE LRRS_FOURIER_MASTER_1
      USE LRRS_CONVERGE
      USE LRRS_WRITE

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

      TYPE(LRRS_Fixed_Inputs), INTENT(IN)       :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT(INOUT) :: LRRS_ModIn
      TYPE(LRRS_Outputs), INTENT(OUT)           :: LRRS_Out

!  Local model input-related variables
!  ===================================

!  Control variables

      LOGICAL :: &
        DO_RRS_OVERALL,      DO_DIRECTRRS_ONLY,    DO_ELASTIC_ONLY, &
        DO_CABANNES_RAMAN,   DO_ENERGY_BALANCING,  DO_MSMODE_LRRS, &
        DO_BIN_REALIZATION,  DO_MONO_REALIZATION,  DO_DIRECT_BEAM, &
        DO_SSCORR_OUTGOING,  DO_SSFULL,            DO_MOLECSCAT_ONLY, &
        DO_PLANE_PARALLEL,   DO_DELTAM_SCALING,    DO_DOUBLE_CONVTEST, &
        DO_UPWELLING,        DO_DNWELLING,         DO_USER_STREAMS, &
        DO_NO_AZIMUTH,       DO_LBOUNDARIES,       DO_MVOUT_ONLY, &
        DO_ADDITIONAL_MVOUT, DO_LAMBERTIAN_SURFACE,DO_WATERLEAVING, &
        DO_SHADOW_EFFECT,    DO_GLITTER_DBMS,  DO_LAMBERTIAN_FRACTION, &
        DO_LRRS_WRITESCENARIO, DO_LRRS_WRITERESULTS, &
        DO_LRRS_WRITEINPUT,    DO_LRRS_WRITEFOURIER

!  Number of discrete ordinate streams

      INTEGER ::   NSTREAMS

!  number of layers

      INTEGER ::   NLAYERS

!  Number of input phase function Legendre moments
!    THIS IS A BASIC INPUT THAT USER MUST PROVIDE

      INTEGER ::   NMOMENTS_INPUT

!  Number of fine layers for sscorr outgoing

      INTEGER ::   NFINELAYERS

!  Solar beam input geometry

      REAL(FPK) :: SOLAR_ANGLE

!  User stream variables

      INTEGER ::   N_USER_STREAMS
      REAL(FPK) :: USER_ANGLES  ( MAX_USER_STREAMS )

!  User azimuth variables

      INTEGER ::   N_USER_RELAZMS
      REAL(FPK) :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)
!    LBOUNDARIES_OUTPUT = layer boundary choices
!    LPARTIALS_OUTPUT = off-boundary level choices

      INTEGER ::   N_LOUTPUT

      INTEGER ::   LBOUNDARIES_OUTPUT ( MAX_LOUTPUT )
      REAL(FPK) :: LPARTIALS_OUTPUT   ( MAX_LOUTPUT )

!  Progress control (If set to 0, no screen output)

      INTEGER ::   PROGRESS

!  BINNING: outer/inner wavelength range, and offset for inner range
!       Depends on the choice of solar spectrum

      INTEGER ::   OFFSET_INNER
      INTEGER ::   NPOINTS_INNER
      INTEGER ::   NPOINTS_OUTER

!  MONOCHROMATIC: Excitation wavelength and index. NPOINTS_MONO = 234

      INTEGER ::   NPOINTS_MONO
      REAL(FPK) :: LAMBDA_EXCIT
      INTEGER ::   W_EXCIT

!  Overall accuracy for convergence criterion

      REAL(FPK) :: LRRS_ACCURACY

!  Small number control

      REAL(FPK) :: LRRS_FGSMALL

!  Flux factor

      REAL(FPK) :: FLUX_FACTOR

!  Earth Radius

      REAL(FPK) :: EARTH_RADIUS

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007 (LIDORT)
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 18 March  2011 (LRRS)

!    This is only required when the outgoing sphericity correction is
!    in operation. Otherwise, the regular pseudo-spherical correction
!    (wiht or without an exact single-scatter correction) uses the same
!    set of angles all the up the nadir from the bottom of atmosphere.

!     This height is normally set equal to the height at the lowest
!     level of the atmosphere: GEOMETRY_SPECHEIGHT = HEIGHT_GRID(NLAYERS)
!     In this case, no adjustment to the geometrical inputs is needed
!     for the outgoing sphericity correction.

!     If there is a situation GEOMETRY_SPECHEIGHT < HEIGHT_GRID(NLAYERS),
!     then an adjustment to the geometrical inputs is needed for the
!     outgoing sphericity correction. This adjustment is internal and
!     the model output will still be given at the geometrical angles
!     as specified by the user, even though these angles may not be the
!     ones at which the calculations were done. This situation will occur
!     when we are given a BOA geometry but we want to make a calculation
!     for a reflecting surface (such as a cloud-top) which is above the
!     BOA level. In this case, GEOMETRY_SPECHEIGHT = 0.0, and the lowest
!     height HEIGHT_GRID(NLAYERS) = cloud-top.

!     This height cannot be greater than HEIGHT_GRID(NLAYERS). If this is
!     the case, this height will be set equal to HEIGHT_GRID(NLAYERS), and
!     the calculation will go through without the adjustment. A warning
!     about this incorrect input choice will be sent to LOGFILE.

      REAL(FPK) :: GEOMETRY_SPECHEIGHT

!  Lambertian fraction

      REAL(FPK) :: LAMBERTIAN_FRACTION

!  Number of BRDF azimuth streams

      INTEGER ::   NSTREAMS_BRDF

!  BRDF parameters.
!     Index choice of BRDF function
!     For Cox-Munk, first parameter  = Wind Speed
!                 second parameter = refractive index
!                 third parameter  = 0.0d

      INTEGER ::   WHICH_BRDF
      INTEGER ::   BRDF_NPARS
      REAL(FPK) :: BRDF_PARS ( MAX_BRDF_PARAMETERS)
      REAL(FPK) :: BRDF_FACTOR
      CHARACTER (LEN=10) :: BRDF_NAMES

!  filenames for output

      CHARACTER (LEN=60) :: DEBUG_FILENAMES(4)

!  Atmospheric quantities
!  ----------------------

!  Height grid

      REAL(FPK) :: HEIGHT_GRID ( 0:MAX_LAYERS )

!  Input layer temperatures, must be in deg K

      REAL(FPK) :: LAYER_TEMPERATURES ( MAX_LAYERS )

!  Input layer Air columns, should be in mol/cm^2 or [DU]

      REAL(FPK) :: LAYER_AIRCOLUMNS ( MAX_LAYERS )

!  Wavelengths and Fluxes
!  ----------------------

!  These are defined on the outer grid for the binning realization
!  These are defined on 234 points for the monochromatic

!  Depends on the choice of solar spectrum

      REAL(FPK) :: LAMBDAS_RANKED  ( MAX_POINTS )
      REAL(FPK) :: FLUXES_RANKED   ( MAX_POINTS )
      REAL(FPK) :: ALBEDOS_RANKED  ( MAX_POINTS )

!  Basic Rayleigh data
!     Rayleigh Cross-sections and depolarization ratios

      REAL(FPK) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  For the binning realization
!  ---------------------------

!  Bin upper and lower limits

      REAL(FPK) :: BINLOWER ( MAX_POINTS )
      REAL(FPK) :: BINUPPER ( MAX_POINTS )

!  Unscaled elastic scattering properties
!  ---------------------------------------

!  These must be defined on the outer wavelength grid (binning)

!  Unscaled quantities, Elastic input

      REAL(FPK) :: DELTAU_INPUT_UNSCALED &
        ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: OMEGAMOMS_ELASTIC_UNSCALED &
        ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

!  ROUTINE OUTPUT
!  ==============

!  ELASTIC FIELD
!  -------------

!   Fourier-summed output

      REAL(FPK) :: ELASTIC_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK) :: ELASTIC_DN &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  single scatter results

      REAL(FPK) :: ELASTIC_SS_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK) :: ELASTIC_SS_DN &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Mean-value output

      REAL(FPK) :: &
        MEAN_ELASTIC_UP ( MAX_LOUTPUT, MAX_POINTS ), &
        MEAN_ELASTIC_DN ( MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK) :: FLUX_ELASTIC_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: FLUX_ELASTIC_DN ( MAX_LOUTPUT, MAX_POINTS )

!  RAMAN FIELD
!  -----------

!  Radiance including inelastic scattering

      REAL(FPK) :: RAMAN_UP &
             ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS)

      REAL(FPK) :: RAMAN_DN &
             ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS)

!  single scatter results

      REAL(FPK) :: RAMAN_SS_UP &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK) :: RAMAN_SS_DN &
        ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Mean-value output

      REAL(FPK) :: MEAN_RAMAN_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: MEAN_RAMAN_DN ( MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK) :: FLUX_RAMAN_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: FLUX_RAMAN_DN ( MAX_LOUTPUT, MAX_POINTS )

!  Number of Fourier terms used

      INTEGER ::   FOURIER_SAVED

!  Error handling
!  --------------

      INTEGER ::   STATUS
      INTEGER ::   N_MESSAGES
      CHARACTER (LEN=70) :: MESSAGES ( MAX_MESSAGES)

!  Local variables
!  ===============

!  bookkeping variables
!  --------------------

!  Number of moments

      INTEGER ::   NMOMENTS

!  Number of convergence tests

      INTEGER ::   N_CONVTESTS

!  number of output streams (quadrature and user)

      INTEGER ::   N_OUT_STREAMS

!  Solar beam input Cosine

      REAL(FPK) :: COS_SZA

!  Discrete ordinate quadrature

      REAL(FPK) :: QUAD_STREAMS ( MAX_STREAMS )
      REAL(FPK) :: QUAD_WEIGHTS ( MAX_STREAMS )
      REAL(FPK) :: QUAD_STRMWGT ( MAX_STREAMS )

!  User stream variables

      REAL(FPK) :: USER_STREAMS ( MAX_USER_STREAMS )

!  Number of Geometries, Geometry offsetting
!    --Added by R. Spurr, 04 August 2009

      INTEGER ::   N_GEOMETRIES
      INTEGER ::   GEOM_OFFSETS ( MAX_USER_STREAMS )

!  Layer masks for doing integrated source terms

      LOGICAL ::   STERM_MASK_UP ( MAX_LAYERS )
      LOGICAL ::   STERM_MASK_DN ( MAX_LAYERS )

      INTEGER ::   LEVELMASK_UP (MAX_LOUTPUT)
      INTEGER ::   LEVELMASK_DN (MAX_LOUTPUT)

!  Fourier output
!  --------------

!  Elastic field

      REAL(FPK) :: ELASTIC_F_UP &
        ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

      REAL(FPK) :: ELASTIC_F_DN &
        ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  Raman Fourier output, not saved

      REAL(FPK) :: RAMAN_F_UP &
        ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK) :: RAMAN_F_DN &
        ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  BRDF quantities
!  ---------------

!  BRDF control and quadrature

      REAL(FPK) :: X_BRDF ( MAX_STREAMS_BRDF )
      REAL(FPK) :: A_BRDF ( MAX_STREAMS_BRDF )

!  BRDFs at quadrature (discrete ordinate) angles

      REAL(FPK) :: BRDFUNC &
        ( MAX_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK) :: BRDFUNC_0 &
        ( MAX_STREAMS, MAX_STREAMS_BRDF )

!  BRDFs at user-defined stream directions

      REAL(FPK) :: USER_BRDFUNC &
        ( MAX_USER_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK) :: USER_BRDFUNC_0 &
        ( MAX_USER_STREAMS, MAX_STREAMS_BRDF )

!  Exact DB values

      REAL(FPK) :: EXACTDB_BRDFUNC &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Beam Attenuations
!  -----------------

!  Solar beam transmittances, average secant factors

      INTEGER ::   BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: BEAM_ETRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK) :: SAVE_TRANS_USERM &
        (  MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  setup variables
!  ---------------

!  Chapman factors

      REAL(FPK) :: CHAPMAN ( MAX_LAYERS, MAX_LAYERS )

!  Basic input quantities for elastic scattering

      REAL(FPK) :: DELTAU_VERT_INPUT ( MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK) :: OMEGAMOMS_ELASTIC &
        ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Truncation factors (zero if no deltam scaling)

      REAL(FPK) :: TRUNC_FACTORS         ( MAX_LAYERS, MAX_POINTS )

!  Bin Mapping quantities (internal derivation)

      INTEGER ::   N_RRSBINS  ( MAX_POINTS )
      INTEGER ::   BINMAP ( MAX_BINS, MAX_POINTS )

!  Chance/Spurr Ring Spectrum. DIAGNOSTIC ONLY.

      REAL(FPK) :: RINGSPEC_1 ( MAX_POINTS )

!  Unscaled Raman optical properties

      REAL(FPK) :: OMEGAMOMS_CABANNES_UNSCALED &
        ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK) :: OMEGAMOMS_RRSLOSS_UNSCALED &
        ( MAX_LAYERS, 0:2, MAX_POINTS )

      REAL(FPK) :: OMEGAMOMS_RRSBIN_UNSCALED &
        ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

      REAL(FPK) :: OMEGAMOMS_RRSGAIN_UNSCALED &
        ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived Cabannes scattering (Internal variable)

      REAL(FPK) :: OMEGAMOMS_CABANNES &
        ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived rotational raman scattering (Internal variable). Loss term
!    Single scatter albedo * phase function moments
!    Only required for the Energy-balance approximation.

      REAL(FPK) :: OMEGAMOMS_RRSLOSS &
        ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived rotational Raman scattering (Internal variable). BIN

      REAL(FPK) :: OMEGAMOMS_RRSBIN &
        ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )

!  Derived rotational Raman scattering (Internal variable). MONO

      REAL(FPK) :: OMEGAMOMS_RRSGAIN &
        ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Help variables
!  --------------

!  Fourier number

      INTEGER ::   FOURIER
      INTEGER ::   N_FOURIER_COMPONENTS

!  Local flags

      LOGICAL ::   ADJUST_SURFACE
      LOGICAL ::   LOCAL_DO_NO_AZIMUTH
      LOGICAL ::   SAVE_DO_NO_AZIMUTH
      LOGICAL ::   LOCAL_ITERATION

!  Modified eradius

      REAL(FPK) :: MODIFIED_ERADIUS

!  Adjusted geometries. New for Version 2.3, March 18, 2011.

      REAL(FPK) :: &
        USER_ANGLES_ADJUST  ( MAX_USER_STREAMS ), &
        SOLAR_ANGLES_ADJUST ( MAX_USER_STREAMS, MAX_USER_RELAZMS ), &
        USER_RELAZMS_ADJUST ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Other variables

      INTEGER ::   UA, UM, LOCAL_PROGRESS
      INTEGER ::   NPOINTS_LOCAL, NPOINTS_TOTAL
      INTEGER ::   FUNIT, RUNIT, IUNIT, SUNIT
      INTEGER ::   TESTCONV, LOCAL_N_USERAZM

!  Fourier cosine arguments

      REAL(FPK) :: AZM_ARGUMENT, DFC
      REAL(FPK) :: AZMFAC ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Exception handling variables

      CHARACTER (LEN=3)  :: C3
      CHARACTER (LEN=70) :: MESSAGE
      CHARACTER (LEN=70) :: MESSAGE_SUB
      CHARACTER (LEN=70) :: POINT_TRACE
      LOGICAL ::   FAIL

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  Fixed Boolean Inputs ;--T/F?

      DO_RRS_OVERALL         = LRRS_FixIn%Bool%DO_RRS_OVERALL
      DO_DIRECTRRS_ONLY      = LRRS_FixIn%Bool%DO_DIRECTRRS_ONLY
      DO_DIRECT_BEAM         = LRRS_FixIn%Bool%DO_DIRECT_BEAM
      DO_SSCORR_OUTGOING     = LRRS_FixIn%Bool%DO_SSCORR_OUTGOING
      DO_SSFULL              = LRRS_FixIn%Bool%DO_SSFULL
      DO_MOLECSCAT_ONLY      = LRRS_FixIn%Bool%DO_MOLECSCAT_ONLY
      DO_PLANE_PARALLEL      = LRRS_FixIn%Bool%DO_PLANE_PARALLEL
      DO_DELTAM_SCALING      = LRRS_FixIn%Bool%DO_DELTAM_SCALING
      DO_DOUBLE_CONVTEST     = LRRS_FixIn%Bool%DO_DOUBLE_CONVTEST
      DO_UPWELLING           = LRRS_FixIn%Bool%DO_UPWELLING
      DO_DNWELLING           = LRRS_FixIn%Bool%DO_DNWELLING
      DO_LBOUNDARIES         = LRRS_FixIn%Bool%DO_LBOUNDARIES
      DO_MVOUT_ONLY          = LRRS_FixIn%Bool%DO_MVOUT_ONLY
      DO_ADDITIONAL_MVOUT    = LRRS_FixIn%Bool%DO_ADDITIONAL_MVOUT
      DO_LAMBERTIAN_SURFACE  = LRRS_FixIn%Bool%DO_LAMBERTIAN_SURFACE
      DO_WATERLEAVING        = LRRS_FixIn%Bool%DO_WATERLEAVING
      DO_SHADOW_EFFECT       = LRRS_FixIn%Bool%DO_SHADOW_EFFECT
      DO_GLITTER_DBMS        = LRRS_FixIn%Bool%DO_GLITTER_DBMS
      DO_LAMBERTIAN_FRACTION = LRRS_FixIn%Bool%DO_LAMBERTIAN_FRACTION
      DO_LRRS_WRITEINPUT     = LRRS_FixIn%Bool%DO_LRRS_WRITEINPUT
      DO_LRRS_WRITESCENARIO  = LRRS_FixIn%Bool%DO_LRRS_WRITESCENARIO
      DO_LRRS_WRITEFOURIER   = LRRS_FixIn%Bool%DO_LRRS_WRITEFOURIER
      DO_LRRS_WRITERESULTS   = LRRS_FixIn%Bool%DO_LRRS_WRITERESULTS

!  Modified Boolean Inputs ;--T/F?

      DO_ELASTIC_ONLY     = LRRS_ModIn%MBool%DO_ELASTIC_ONLY
      DO_CABANNES_RAMAN   = LRRS_ModIn%MBool%DO_CABANNES_RAMAN
      DO_ENERGY_BALANCING = LRRS_ModIn%MBool%DO_ENERGY_BALANCING
      DO_MSMODE_LRRS      = LRRS_ModIn%MBool%DO_MSMODE_LRRS
      DO_BIN_REALIZATION  = LRRS_ModIn%MBool%DO_BIN_REALIZATION
      DO_MONO_REALIZATION = LRRS_ModIn%MBool%DO_MONO_REALIZATION
      DO_USER_STREAMS     = LRRS_ModIn%MBool%DO_USER_STREAMS
      DO_NO_AZIMUTH       = LRRS_ModIn%MBool%DO_NO_AZIMUTH

!  Fixed Control Inputs

      PROGRESS        = LRRS_FixIn%Cont%PROGRESS
      NLAYERS         = LRRS_FixIn%Cont%NLAYERS
      NSTREAMS        = LRRS_FixIn%Cont%NSTREAMS
      NMOMENTS_INPUT  = LRRS_FixIn%Cont%NMOMENTS_INPUT
      NFINELAYERS     = LRRS_FixIn%Cont%NFINELAYERS
      FLUX_FACTOR     = LRRS_FixIn%Cont%FLUX_FACTOR
      SOLAR_ANGLE     = LRRS_FixIn%Cont%SOLAR_ANGLE
      EARTH_RADIUS    = LRRS_FixIn%Cont%EARTH_RADIUS
      LRRS_ACCURACY   = LRRS_FixIn%Cont%LRRS_ACCURACY
      LRRS_FGSMALL    = LRRS_FixIn%Cont%LRRS_FGSMALL
      DEBUG_FILENAMES = LRRS_FixIn%Cont%DEBUG_FILENAMES

!  Fixed User-Value Inputs

      N_USER_STREAMS     = LRRS_FixIn%UserVal%N_USER_STREAMS
      USER_ANGLES        = LRRS_FixIn%UserVal%USER_ANGLES
      N_USER_RELAZMS     = LRRS_FixIn%UserVal%N_USER_RELAZMS
      USER_RELAZMS       = LRRS_FixIn%UserVal%USER_RELAZMS
      N_LOUTPUT          = LRRS_FixIn%UserVal%N_LOUTPUT
      LBOUNDARIES_OUTPUT = LRRS_FixIn%UserVal%LBOUNDARIES_OUTPUT
      LPARTIALS_OUTPUT   = LRRS_FixIn%UserVal%LPARTIALS_OUTPUT

!  Fixed Spectral Inputs
      NPOINTS_MONO   = LRRS_FixIn%Spect%NPOINTS_MONO
      LAMBDAS_RANKED = LRRS_FixIn%Spect%LAMBDAS_RANKED
      FLUXES_RANKED  = LRRS_FixIn%Spect%FLUXES_RANKED
      ! write(*,*) NPOINTS_MONO
      ! write(*,*) LAMBDAS_RANKED(1:2), LAMBDAS_RANKED(NPOINTS_MONO-1:NPOINTS_MONO)
      ! write(*,*) FLUXES_RANKED(1:2), FLUXES_RANKED(NPOINTS_MONO-1:NPOINTS_MONO)

      IF (DO_MONO_REALIZATION) THEN
        !For monochromatic calculations:
        NPOINTS_MONO   = LRRS_FixIn%Spect%NPOINTS_MONO
        LAMBDA_EXCIT   = LRRS_FixIn%Spect%LAMBDA_EXCIT
        W_EXCIT        = LRRS_FixIn%Spect%W_EXCIT
      !   write(*,*) NPOINTS_MONO, LAMBDA_EXCIT, W_EXCIT
      !   pause
      END IF

      IF (DO_BIN_REALIZATION) THEN
        !For binning calculations:
        NPOINTS_INNER  = LRRS_FixIn%Spect%NPOINTS_INNER
        OFFSET_INNER   = LRRS_FixIn%Spect%OFFSET_INNER
        NPOINTS_OUTER  = LRRS_FixIn%Spect%NPOINTS_OUTER
        BINLOWER       = LRRS_FixIn%Spect%BINLOWER
        BINUPPER       = LRRS_FixIn%Spect%BINUPPER
      END IF

!  Fixed Atmosphere Inputs

      HEIGHT_GRID        = LRRS_FixIn%Atmos%HEIGHT_GRID
      LAYER_TEMPERATURES = LRRS_FixIn%Atmos%LAYER_TEMPERATURES
      LAYER_AIRCOLUMNS   = LRRS_FixIn%Atmos%LAYER_AIRCOLUMNS
      RAYLEIGH_XSEC      = LRRS_FixIn%Atmos%RAYLEIGH_XSEC
      RAYLEIGH_DEPOL     = LRRS_FixIn%Atmos%RAYLEIGH_DEPOL

      DELTAU_INPUT_UNSCALED      = &
        LRRS_FixIn%Atmos%DELTAU_INPUT_UNSCALED
      OMEGAMOMS_ELASTIC_UNSCALED = &
        LRRS_FixIn%Atmos%OMEGAMOMS_ELASTIC_UNSCALED

!  Modified Atmosphere Inputs

      GEOMETRY_SPECHEIGHT = LRRS_ModIn%MAtmos%GEOMETRY_SPECHEIGHT

!  Fixed Surface Inputs

      LAMBERTIAN_FRACTION = LRRS_FixIn%Surf%LAMBERTIAN_FRACTION
      NSTREAMS_BRDF       = LRRS_FixIn%Surf%NSTREAMS_BRDF
      BRDF_NAMES          = LRRS_FixIn%Surf%BRDF_NAMES
      WHICH_BRDF          = LRRS_FixIn%Surf%WHICH_BRDF
      !Note: BRDF_FACTOR not used at present
      BRDF_FACTOR         = LRRS_FixIn%Surf%BRDF_FACTOR
      BRDF_NPARS          = LRRS_FixIn%Surf%BRDF_NPARS
      BRDF_PARS           = LRRS_FixIn%Surf%BRDF_PARS
      ALBEDOS_RANKED      = LRRS_FixIn%Surf%ALBEDOS_RANKED

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Write input to file (debug only)

!      CALL LRRS_WRITEINPUT_ALL ( &
!        DO_RRS_OVERALL,DO_DIRECTRRS_ONLY,DO_DIRECT_BEAM,&
!        DO_SSCORR_OUTGOING,DO_SSFULL,DO_MOLECSCAT_ONLY,&
!        DO_PLANE_PARALLEL,DO_DELTAM_SCALING,DO_DOUBLE_CONVTEST,&
!        DO_UPWELLING,DO_DNWELLING,DO_LBOUNDARIES,DO_MVOUT_ONLY,&
!        DO_ADDITIONAL_MVOUT,DO_LAMBERTIAN_SURFACE,DO_WATERLEAVING,&
!        DO_SHADOW_EFFECT,DO_GLITTER_DBMS,DO_LAMBERTIAN_FRACTION,&
!        DO_LRRS_WRITESCENARIO,DO_LRRS_WRITERESULTS,DO_LRRS_WRITEINPUT,&
!        DO_LRRS_WRITEFOURIER,DO_ELASTIC_ONLY,DO_CABANNES_RAMAN,&
!        DO_ENERGY_BALANCING,DO_MSMODE_LRRS,DO_BIN_REALIZATION,&
!        DO_MONO_REALIZATION,DO_USER_STREAMS,DO_NO_AZIMUTH,&
!        PROGRESS,NLAYERS,NSTREAMS,NMOMENTS_INPUT,NFINELAYERS,&
!        FLUX_FACTOR,LRRS_ACCURACY,LRRS_FGSMALL,EARTH_RADIUS,&
!        SOLAR_ANGLE,DEBUG_FILENAMES,N_USER_STREAMS,USER_ANGLES,&
!        N_USER_RELAZMS,USER_RELAZMS,N_LOUTPUT,LBOUNDARIES_OUTPUT,&
!        LPARTIALS_OUTPUT,LAMBDAS_RANKED,FLUXES_RANKED,NPOINTS_MONO,&
!        LAMBDA_EXCIT,W_EXCIT,NPOINTS_INNER,OFFSET_INNER,NPOINTS_OUTER,&
!        BINLOWER,BINUPPER,HEIGHT_GRID,LAYER_TEMPERATURES,&
!        LAYER_AIRCOLUMNS,RAYLEIGH_XSEC,RAYLEIGH_DEPOL,&
!        DELTAU_INPUT_UNSCALED,OMEGAMOMS_ELASTIC_UNSCALED,&
!        LAMBERTIAN_FRACTION,NSTREAMS_BRDF,WHICH_BRDF,BRDF_NPARS,&
!        BRDF_PARS,BRDF_FACTOR,BRDF_NAMES,ALBEDOS_RANKED)

!  Initialize outputs
!  ==================

!  Main Outputs

      !Elastic Field
      ELASTIC_UP      = ZERO
      ELASTIC_DN      = ZERO
      ELASTIC_SS_UP   = ZERO
      ELASTIC_SS_DN   = ZERO
      MEAN_ELASTIC_UP = ZERO
      MEAN_ELASTIC_DN = ZERO
      FLUX_ELASTIC_UP = ZERO
      FLUX_ELASTIC_DN = ZERO

      !Raman Field
      RAMAN_UP        = ZERO
      RAMAN_DN        = ZERO
      RAMAN_SS_UP     = ZERO
      RAMAN_SS_DN     = ZERO
      MEAN_RAMAN_UP   = ZERO
      MEAN_RAMAN_DN   = ZERO
      FLUX_RAMAN_UP   = ZERO
      FLUX_RAMAN_DN   = ZERO

      !Number of Fourier terms used
      FOURIER_SAVED   = ZERO

!  Status Outputs

      STATUS     = LRRS_SUCCESS
      N_MESSAGES = 0
      MESSAGES   = ' '

!  Internal settings
!  =================

!  Elastic-only flag

      DO_ELASTIC_ONLY = .NOT.DO_RRS_OVERALL.AND..NOT.DO_DIRECTRRS_ONLY

!  exclusion flags

      DO_CABANNES_RAMAN = .not. DO_ENERGY_BALANCING
      IF ( .NOT. DO_ELASTIC_ONLY ) THEN
        DO_MONO_REALIZATION = .not.DO_BIN_REALIZATION
      ENDIF

!  Check inputs. Exit if any failure
!  =================================

      CALL LRRS_CHECK_INPUTS &
       ( LRRS_FixIn, LRRS_ModIn, N_USER_STREAMS, N_USER_RELAZMS, N_LOUTPUT, &
         SOLAR_ANGLE, USER_ANGLES, USER_RELAZMS, NLAYERS, HEIGHT_GRID, &
         LBOUNDARIES_OUTPUT, LPARTIALS_OUTPUT, EARTH_RADIUS, WHICH_BRDF, &
         LAMBERTIAN_FRACTION, BRDF_NPARS, BRDF_PARS, BRDF_NAMES, &
         STATUS, N_MESSAGES, MESSAGES )

      IF ( STATUS .NE. LRRS_SUCCESS ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = 'Failure from Input checking, Main'
        LRRS_Out%Status%STATUS_CALCULATION = STATUS
        LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
        LRRS_Out%Status%MESSAGES           = MESSAGES
        RETURN
      ENDIF

!  If there's no azimuth dependence, just do one value in azimuth loop

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_N_USERAZM = 1
        N_USER_RELAZMS  = 1
       ELSE
        LOCAL_N_USERAZM = N_USER_RELAZMS
      ENDIF

!  Set Number of Geometries, and the offsetting

      N_GEOMETRIES = LOCAL_N_USERAZM * N_USER_STREAMS
      DO UM = 1, N_USER_STREAMS
        GEOM_OFFSETS(UM) = LOCAL_N_USERAZM * (UM - 1)
      ENDDO

!  Geometry adjustment. New section 18 March 2011
!  ==============================================

!  Adjust surface condition

      ADJUST_SURFACE = .FALSE.
      IF ( DO_SSCORR_OUTGOING ) THEN
        IF ( HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) THEN
         ADJUST_SURFACE = .TRUE.
        ENDIF
      ENDIF

!  Perform adjustment

      MODIFIED_ERADIUS = EARTH_RADIUS + GEOMETRY_SPECHEIGHT
      CALL MULTI_OUTGOING_ADJUSTGEOM &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS, &
          N_USER_STREAMS,   N_USER_RELAZMS, &
          HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE, &
          USER_ANGLES,  SOLAR_ANGLE, USER_RELAZMS, &
          USER_ANGLES_ADJUST, SOLAR_ANGLES_ADJUST, USER_RELAZMS_ADJUST, &
          FAIL, MESSAGE_SUB, MESSAGE )

!  Exception handling

      IF ( FAIL ) THEN
        STATUS = LRRS_SERIOUS
        MESSAGES(N_MESSAGES+1) = MESSAGE_SUB
        MESSAGES(N_MESSAGES+2) = MESSAGE
        MESSAGES(N_MESSAGES+3) = ' Failure in multi_outgoing_adjustgeom'
        N_MESSAGES = N_MESSAGES + 3
        LRRS_Out%Status%STATUS_CALCULATION = STATUS
        LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
        LRRS_Out%Status%MESSAGES           = MESSAGES
        RETURN
      ENDIF

!  Derive inputs - bookkeeping!
!  ============================

      CALL LRRS_DERIVE_INPUTS &
        ( DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, &
          DO_BIN_REALIZATION, DO_MOLECSCAT_ONLY, DO_LBOUNDARIES, &
          NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS, N_LOUTPUT, &
          NLAYERS, NPOINTS_INNER, SOLAR_ANGLE, USER_ANGLES, &
          NMOMENTS, N_OUT_STREAMS, LBOUNDARIES_OUTPUT, LPARTIALS_OUTPUT, &
          FLUX_FACTOR, COS_SZA, LEVELMASK_UP, LEVELMASK_DN, &
          STERM_MASK_UP, STERM_MASK_DN, N_CONVTESTS, &
          QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS )

!  Setups for Raman problem
!  ========================

!    - delta-m scaling
!    - Raman spectroscopy
!    - Raman optical depths

      IF ( DO_BIN_REALIZATION ) THEN
        CALL SETUP_MASTER_BIN &
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
          FAIL, MESSAGE, MESSAGE_SUB )
      ELSE
        CALL SETUP_MASTER_MONO &
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
            
            ! write(*,*)OMEGAMOMS_RRSLOSS(1,2,:)
            ! write(*,*)
            ! write(*,*)OMEGAMOMS_RRSGAIN(1,2,:)
            ! pause
            ! write(*,*)TRUNC_FACTORS
            
      ENDIF

!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@
!  Exception handling for MAX_BINS Dimension
      IF ( DO_BIN_REALIZATION ) THEN
         IF  ( FAIL ) THEN
            STATUS = LRRS_SERIOUS
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = MESSAGE
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = MESSAGE_SUB
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = 'Failure from SETUP_MASTER_BIN'
            LRRS_Out%Status%STATUS_CALCULATION = STATUS
            LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
            LRRS_Out%Status%MESSAGES           = MESSAGES
            RETURN
         ENDIF
      ENDIF
!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@

!  Chapman function calculation of CHAPMAN_FACTORS

      CALL CHAPMAN_FUNCTION &
           ( DO_PLANE_PARALLEL, NLAYERS, &
             COS_SZA, EARTH_RADIUS, HEIGHT_GRID, &
             CHAPMAN )
            !  write(*,*)CHAPMAN(1:NLAYERS,1)
            ! pause
!  BRDF input functions (non-Lambertian)

      IF ( .NOT. DO_LAMBERTIAN_SURFACE ) THEN
        CALL LRRS_BRDF_MASTER &
          ( WHICH_BRDF, DO_GLITTER_DBMS, BRDF_PARS, &
            DO_UPWELLING, DO_USER_STREAMS, &
            DO_SSCORR_OUTGOING, DO_SSFULL, &
            NSTREAMS_BRDF,  SOLAR_ANGLE, &
            NSTREAMS,       QUAD_STREAMS, &
            N_USER_STREAMS, USER_STREAMS, &
            N_USER_RELAZMS, USER_RELAZMS, &
            X_BRDF, A_BRDF, EXACTDB_BRDFUNC, &
            BRDFUNC, BRDFUNC_0, &
            USER_BRDFUNC, USER_BRDFUNC_0 )
      ENDIF

!  Write input variables (if flagged)

      IF ( DO_LRRS_WRITEINPUT ) THEN
        IUNIT = LRRS_INUNIT
        OPEN(IUNIT,FILE=DEBUG_FILENAMES(1),STATUS='UNKNOWN')
!        CALL LRRS_WRITEINPUT ( IUNIT )
        CLOSE(IUNIT)
      ENDIF

!  Write scenario variables (if flagged)

      IF ( DO_LRRS_WRITESCENARIO ) THEN
        SUNIT = LRRS_SCENUNIT
        OPEN(SUNIT,FILE=DEBUG_FILENAMES(2),STATUS='UNKNOWN')
!        CALL LRRS_WRITESCENARIO ( SUNIT )
        CLOSE(SUNIT)
      ENDIF

!  Get the Single scatter elastic correction. Now includes the DB term
!    DB correction is automatic now, even for the Lambertian surface
!      28 May 2008,  RT Solutions Inc.

      IF ( DO_SSCORR_OUTGOING ) THEN
        IF ( DO_BIN_REALIZATION ) THEN
          CALL LRRS_SSCORR_OUTGOING_BIN_1 &
         ( DO_RRS_OVERALL, DO_UPWELLING, DO_DNWELLING, &
           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_WATERLEAVING, &
           DO_LAMBERTIAN_SURFACE, DO_LAMBERTIAN_FRACTION, &
           DO_DELTAM_SCALING, NLAYERS, NFINELAYERS, NMOMENTS_INPUT, &
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
           FAIL, MESSAGE_SUB )

!  debug 27 may 2011
!         write(*,*)raman_ss_up(1,1,1)
!         pause

          IF ( FAIL ) THEN
            STATUS = LRRS_SERIOUS
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = MESSAGE_SUB
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = 'Failure from LRRS_SSCORR_OUTGOING_BIN_1'
            LRRS_Out%Status%STATUS_CALCULATION = STATUS
            LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
            LRRS_Out%Status%MESSAGES           = MESSAGES
            RETURN
          ENDIF

        ELSE IF ( DO_MONO_REALIZATION ) THEN
          CALL LRRS_SSCORR_OUTGOING_MONO_1 &
         ( DO_RRS_OVERALL, DO_UPWELLING, DO_DNWELLING, &
           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_WATERLEAVING, &
           DO_LAMBERTIAN_SURFACE, DO_LAMBERTIAN_FRACTION, &
           DO_DELTAM_SCALING, NLAYERS, NFINELAYERS, NMOMENTS_INPUT, &
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
           FAIL, MESSAGE_SUB)

      !      write(*,*)OMEGAMOMS_RRSGAIN_UNSCALED(1,2,:)
      !      write(*,*)DELTAU_VERT_INPUT(1,:)
      !      write(*,*)
           
      !      write(*,*),ELASTIC_SS_UP(1,1,:)
      !      write(*,*)
      !      write(*,*),RAMAN_SS_UP(1,1,:)
           
      !      pause
           !--------note mchoi: no problem untill here.

          IF ( FAIL ) THEN
            STATUS = LRRS_SERIOUS
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = MESSAGE_SUB
            N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = 'Failure from LRRS_SSCORR_OUTGOING_MONO_1'
            LRRS_Out%Status%STATUS_CALCULATION = STATUS
            LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
            LRRS_Out%Status%MESSAGES           = MESSAGES
            RETURN
          ENDIF
        ENDIF
      ENDIF

!  Npoints local

      IF ( DO_BIN_REALIZATION ) THEN
        NPOINTS_TOTAL = NPOINTS_OUTER
        NPOINTS_LOCAL = NPOINTS_INNER
      ELSE
        NPOINTS_TOTAL = NPOINTS_MONO
        NPOINTS_LOCAL = 1
      ENDIF

!  Get the Beam average secants, transmittances

      IF ( .NOT. DO_SSFULL ) THEN
        CALL ATTENUATION_SETUP_1 &
          ( NPOINTS_TOTAL, NLAYERS, &
            DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS, &
            DELTAU_VERT_INPUT, CHAPMAN, &
            BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, &
            BEAM_ETRANS, BEAM_DTRANS, SAVE_TRANS_USERM )
            ! write(*,*)BEAM_PICUTOFF(:)
            ! write(*,*)BEAM_ETRANS(1,:)
            ! write(*,*)
            ! write(*,*)SAVE_TRANS_USERM(1,1,:)
            ! pause
      ENDIF

!  Debug, 24 September 2015
!      do um = 1, npoints_inner
!        write(*,'(i4,1p6e15.6)')UM,RAMAN_SS_UP(1:3,3:4,UM)
!      enddo
!      pause

!  Initialise Fourier loop
!  =======================

!  Set Number of Fourier terms (NMOMENTS = Maximum).
!    ( Starting from 0 = Fundamental )

      SAVE_DO_NO_AZIMUTH  = DO_NO_AZIMUTH
      LOCAL_DO_NO_AZIMUTH = DO_NO_AZIMUTH

!  No azimuth dependency for following cases

      IF ( DO_MVOUT_ONLY  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

!  set Fourier number (2 for Rayleigh only)

      IF ( LOCAL_DO_NO_AZIMUTH .OR. DO_SSFULL ) THEN
        N_FOURIER_COMPONENTS = 0
      ELSE
        IF ( DO_MOLECSCAT_ONLY  ) THEN
          N_FOURIER_COMPONENTS = 2
        ELSE
          N_FOURIER_COMPONENTS = NMOMENTS
        ENDIF
      ENDIF

!  re-set no-azimuth flag

      DO_NO_AZIMUTH = LOCAL_DO_NO_AZIMUTH

!  Fourier loop
!  ============

!  pre-set counters and flag

      LOCAL_ITERATION = .TRUE.
      FOURIER  = -1
      TESTCONV = 0

!  start of loop

      DO WHILE ( LOCAL_ITERATION .AND. &
                    FOURIER.LT.N_FOURIER_COMPONENTS )

!  Go to final settings, is only the SS is required

        IF ( DO_SSFULL ) GO TO 7443

!  Fourier counter

        FOURIER = FOURIER + 1

!  Local progress counter

        IF ( FOURIER .GT. 2 ) THEN
          LOCAL_PROGRESS = 10000
        ELSE
          LOCAL_PROGRESS = PROGRESS
        ENDIF

!  azimuth cosine factors
!    Adjusted Azimuths, now. 18 March 2011, in line with LIDORT 3.5

        IF ( FOURIER .GT. 0 ) THEN
          DFC = DBLE(FOURIER)
          DO UA = 1, LOCAL_N_USERAZM
            DO UM = 1, N_USER_STREAMS
              AZM_ARGUMENT = USER_RELAZMS_ADJUST(UM,UA) * DFC
              AZMFAC(UM,UA)   = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
            ENDDO
          ENDDO
        ENDIF

!  Progress - could get rid of this

!        IF (DO_ELASTIC_ONLY) THEN
!          write(*,*)' ..Elastic solution: Fourier component',FOURIER
!        ELSE
      !     write(*,*)' ..RRS solution: Fourier component',FOURIER
!        ENDIF

!  Main call to LRRS Fourier module

        CALL RAMAN_FOURIER_1 &
       ( DO_RRS_OVERALL, DO_DIRECTRRS_ONLY, DO_MSMODE_LRRS, &
         LOCAL_PROGRESS, DO_ELASTIC_ONLY, DO_SSCORR_OUTGOING, DO_SSFULL, &
         DO_BIN_REALIZATION, DO_ENERGY_BALANCING, &
         DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_DIRECT_BEAM, &
         DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, &
         DO_LAMBERTIAN_SURFACE, DO_WATERLEAVING, &
         DO_LAMBERTIAN_FRACTION, LAMBERTIAN_FRACTION, &
         NPOINTS_OUTER, OFFSET_INNER, NPOINTS_INNER, NPOINTS_MONO, &
         N_RRSBINS, BINMAP, W_EXCIT, NSTREAMS_BRDF, X_BRDF, A_BRDF, &
         BRDFUNC, BRDFUNC_0, USER_BRDFUNC, USER_BRDFUNC_0, &
         NLAYERS, NSTREAMS, NMOMENTS, N_LOUTPUT, N_USER_STREAMS, &
         FOURIER, COS_SZA, LRRS_FGSMALL, FLUX_FACTOR, &
         QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS, &
         STERM_MASK_UP, LEVELMASK_UP, STERM_MASK_DN, LEVELMASK_DN, &
         DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, OMEGAMOMS_CABANNES, &
         OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN, OMEGAMOMS_RRSGAIN, &
         ALBEDOS_RANKED, FLUXES_RANKED, &
         BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, &
         BEAM_DTRANS, BEAM_ETRANS, SAVE_TRANS_USERM, &
         ELASTIC_F_UP, MEAN_ELASTIC_UP, FLUX_ELASTIC_UP, &
         ELASTIC_F_DN, MEAN_ELASTIC_DN, FLUX_ELASTIC_DN, &
         RAMAN_F_UP, MEAN_RAMAN_UP, FLUX_RAMAN_UP, &
         RAMAN_F_DN, MEAN_RAMAN_DN, FLUX_RAMAN_DN, &
         FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE )

      !    write(*,*)ELASTIC_F_UP(1,1,1)
      !    write(*,*)RAMAN_F_UP(1,1,:)
      !    pause
         
      !    write(*,*) sum(RAMAN_F_UP) ! why.. all 0?
      !    write(*,*)MAX_USER_STREAMS !1
      !    write(*,*)SAVE_TRANS_USERM(1,1,:) !filled as '0.999~~~' (not identical value per points)
      !    write(*,*)BEAM_ETRANS(1,:)! filled as '0.999~~~' (not identical value per points)
      !    write(*,*)BEAM_DTRANS(1,:)! filled as '0.999~~~' (not identical value per points)
      !    write(*,*)BEAM_AVSECANT(1,:) ! filled as '1.413~~'
            ! write(*,*)BEAM_ITRANS(1, :) ! same; filled as '1'
      !    write(*,*)BEAM_PICUTOFF ! filled with nlayers 24
      !    write(*,*)STERM_MASK_UP ! different but from nlayers
      !    write(*,*)
      !    write(*,*)QUAD_WEIGHTS ! similar
      !    write(*,*)QUAD_STRMWGT ! different
      !    write(*,*)QUAD_WEIGHTS ! same
      !    write(*,*)QUAD_STREAMS ! same
      !    write(*,*)NSTREAMS, NMOMENTS ! 4, 2; same
      !    write(*,*)FLUXES_RANKED(:)
      !      write(*,*)OMEGAMOMS_RRSLOSS(1,2,:)
      !      write(*,*)OMEGAMOMS_RRSGAIN(1,2,:)
      !    write(*,*)ELASTIC_F_UP(1,1,:)
      !    write(*,*)
      !    write(*,*)RAMAN_F_UP(1,1,:)
      !    write(*,*)N_USER_STREAMS ! = 1
            ! write(*,*)DO_MVOUT_ONLY ! = F
      !       write(*,*)ELASTIC_F_UP(1,1,:)
      !    pause
      !    write(*,*) NPOINTS_OUTER, OFFSET_INNER, NPOINTS_INNER, NPOINTS_MONO, &
      !    N_RRSBINS
      !    pause

!  error handling. Add to collection of messages

        C3 = '   '
        IF ( FAIL ) THEN
          STATUS = LRRS_SERIOUS
          IF (FOURIER.LT.10) WRITE(C3(3:3),'(I1)') FOURIER
          IF (FOURIER.GE.10) WRITE(C3(2:3),'(I2)') FOURIER
          MESSAGES(N_MESSAGES+1)='Error from RAMAN_FOURIER_1, m = '//C3
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES+1) = MESSAGE
          N_MESSAGES = N_MESSAGES + 1
          MESSAGES(N_MESSAGES+1) = MESSAGE_SUB
          N_MESSAGES = N_MESSAGES + 1
          IF ( DO_BIN_REALIZATION ) THEN
            MESSAGES(N_MESSAGES+1) = POINT_TRACE
            N_MESSAGES = N_MESSAGES + 1
          ENDIF
          LRRS_Out%Status%STATUS_CALCULATION = STATUS
          LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
          LRRS_Out%Status%MESSAGES           = MESSAGES
          RETURN
        ENDIF

!  Continuation point for SSFULL

 7443   continue

!  Convergence examination

        IF ( .not. DO_MVOUT_ONLY  ) THEN
          CALL RAMAN_CONVERGE &
       ( DO_UPWELLING, DO_DNWELLING, DO_DOUBLE_CONVTEST, &
         DO_SSCORR_OUTGOING, DO_SSFULL, DO_MOLECSCAT_ONLY, &
         DO_NO_AZIMUTH, DO_RRS_OVERALL, DO_ELASTIC_ONLY, FOURIER, &
         NPOINTS_LOCAL, N_LOUTPUT, N_OUT_STREAMS, LOCAL_N_USERAZM, &
         N_CONVTESTS, N_GEOMETRIES, GEOM_OFFSETS, AZMFAC, LRRS_ACCURACY, &
         ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN, &
         ELASTIC_F_UP, ELASTIC_F_DN, RAMAN_F_UP, RAMAN_F_DN, &
         ELASTIC_UP, ELASTIC_DN, RAMAN_UP, RAMAN_DN, &
         FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )
        ENDIF

!  Exit if SSFULL

        if ( DO_SSFULL ) return

!  Fourier output
!  --------------

!  Open file if Fourier = 0
!  Write Standard Fourier output
!  Close file if iteration has finished

        IF ( DO_LRRS_WRITEFOURIER ) THEN
          FUNIT = LRRS_FUNIT
          IF ( FOURIER .EQ. 0 ) THEN
            OPEN(FUNIT,FILE=DEBUG_FILENAMES(3),STATUS='UNKNOWN')
          ENDIF
!          CALL LRRS_WRITEFOURIER ( FUNIT, FOURIER )
          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
        ENDIF

!  end iteration loop

      ENDDO

!  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

!  Major result output
!  ===================

!  Standard output
!  - - - - - - - -

      IF ( DO_LRRS_WRITERESULTS ) THEN
        RUNIT = LRRS_RESUNIT
        OPEN(RUNIT,FILE=DEBUG_FILENAMES(4),STATUS='UNKNOWN')
!        CALL LRRS_WRITERESULTS ( RUNIT )
        CLOSE(RUNIT)
      ENDIF

!  ========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (IN/OUT variables)
!  ========================================================

!  Modified Boolean Inputs

      LRRS_ModIn%MBool%DO_ELASTIC_ONLY     = DO_ELASTIC_ONLY
      LRRS_ModIn%MBool%DO_CABANNES_RAMAN   = DO_CABANNES_RAMAN
      LRRS_ModIn%MBool%DO_ENERGY_BALANCING = DO_ENERGY_BALANCING
      LRRS_ModIn%MBool%DO_MSMODE_LRRS      = DO_MSMODE_LRRS
      LRRS_ModIn%MBool%DO_BIN_REALIZATION  = DO_BIN_REALIZATION
      LRRS_ModIn%MBool%DO_MONO_REALIZATION = DO_MONO_REALIZATION
      LRRS_ModIn%MBool%DO_USER_STREAMS     = DO_USER_STREAMS
      LRRS_ModIn%MBool%DO_NO_AZIMUTH       = DO_NO_AZIMUTH

!  ==========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (pure OUT variables)
!  ==========================================================

!  Main Outputs

      !Elastic Field
      LRRS_Out%Main%ELASTIC_UP      = ELASTIC_UP
      LRRS_Out%Main%ELASTIC_DN      = ELASTIC_DN
      LRRS_Out%Main%ELASTIC_SS_UP   = ELASTIC_SS_UP
      LRRS_Out%Main%ELASTIC_SS_DN   = ELASTIC_SS_DN
      LRRS_Out%Main%MEAN_ELASTIC_UP = MEAN_ELASTIC_UP
      LRRS_Out%Main%MEAN_ELASTIC_DN = MEAN_ELASTIC_DN
      LRRS_Out%Main%FLUX_ELASTIC_UP = FLUX_ELASTIC_UP
      LRRS_Out%Main%FLUX_ELASTIC_DN = FLUX_ELASTIC_DN

      !Raman Field
      LRRS_Out%Main%RAMAN_UP        = RAMAN_UP
      LRRS_Out%Main%RAMAN_DN        = RAMAN_DN
      LRRS_Out%Main%RAMAN_SS_UP     = RAMAN_SS_UP
      LRRS_Out%Main%RAMAN_SS_DN     = RAMAN_SS_DN
      LRRS_Out%Main%MEAN_RAMAN_UP   = MEAN_RAMAN_UP
      LRRS_Out%Main%MEAN_RAMAN_DN   = MEAN_RAMAN_DN
      LRRS_Out%Main%FLUX_RAMAN_UP   = FLUX_RAMAN_UP
      LRRS_Out%Main%FLUX_RAMAN_DN   = FLUX_RAMAN_DN

      !Number of Fourier terms used
      LRRS_Out%Main%FOURIER_SAVED   = FOURIER_SAVED

!  Status Outputs

      LRRS_Out%Status%STATUS_CALCULATION = STATUS
      LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
      LRRS_Out%Status%MESSAGES           = MESSAGES

!  ===================================
!  END COPY LOCAL VARIABLES TO OUTPUTS
!  ===================================

!  Finish

      RETURN
      END SUBROUTINE RAMAN_MASTER


      END MODULE lrrs_main_master

