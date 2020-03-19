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

      MODULE LRRS_Inputs_Def

!  This Module contains the following LIDORT_RRS Input Structures,
!  with Intents :

!            LRRS_Fixed_Inputs    Intent(In)
!         LRRS_Modified_Inputs    Intent(InOut)

      USE LRRS_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LRRS_Fixed_Boolean


!  Top level control

      LOGICAL :: DO_RRS_OVERALL
      LOGICAL :: DO_DIRECTRRS_ONLY

!  Direct beam flag

      LOGICAL :: DO_DIRECT_BEAM

!  Single scattering flags

      LOGICAL :: DO_SSCORR_OUTGOING
      LOGICAL :: DO_SSFULL

!  Molecular scattering only (no aerosols/clouds)

      LOGICAL :: DO_MOLECSCAT_ONLY

!  Plane parallel flag

      LOGICAL :: DO_PLANE_PARALLEL

!  Delta-M scaling (only with aerosols)

      LOGICAL :: DO_DELTAM_SCALING

!  Double convergence test (only with aerosols)

      LOGICAL :: DO_DOUBLE_CONVTEST

!  Upwelling/downwelling flags

      LOGICAL :: DO_UPWELLING
      LOGICAL :: DO_DNWELLING

!  Level boundaries flag

      LOGICAL :: DO_LBOUNDARIES

!  Mean-value output options

      LOGICAL :: DO_MVOUT_ONLY
      LOGICAL :: DO_ADDITIONAL_MVOUT

!  Surface flags

      LOGICAL :: DO_LAMBERTIAN_SURFACE
      LOGICAL :: DO_WATERLEAVING
      LOGICAL :: DO_SHADOW_EFFECT
      LOGICAL :: DO_GLITTER_DBMS
      LOGICAL :: DO_LAMBERTIAN_FRACTION

!  Debug write flags

      LOGICAL :: DO_LRRS_WRITEINPUT
      LOGICAL :: DO_LRRS_WRITESCENARIO
      LOGICAL :: DO_LRRS_WRITEFOURIER
      LOGICAL :: DO_LRRS_WRITERESULTS


      END TYPE LRRS_Fixed_Boolean

! #####################################################################
! #####################################################################

      TYPE LRRS_Fixed_Control


!  Progress control (If set to 0, no screen output)

      INTEGER :: PROGRESS

!  Number of layers

      INTEGER :: NLAYERS

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  number of moments in the multiple scatter calculation

      INTEGER :: NMOMENTS_INPUT

!  Number of fine layers for sscorr outgoing

      INTEGER :: NFINELAYERS

!  Flux factor

      REAL(FPK) :: FLUX_FACTOR

!  Solar beam input geometry

      REAL(FPK) :: SOLAR_ANGLE

!  Earth Radius

      REAL(FPK) :: EARTH_RADIUS

!  Overall accuracy for convergence criterion

      REAL(FPK) :: LRRS_ACCURACY

!  Small number control

      REAL(FPK) :: LRRS_FGSMALL

!  filenames for output

      CHARACTER (LEN=60), DIMENSION(4) :: DEBUG_FILENAMES


      END TYPE LRRS_Fixed_Control

! #####################################################################
! #####################################################################

      TYPE LRRS_Fixed_UserValues


!  User stream variables

      INTEGER :: N_USER_STREAMS
      REAL(FPK), DIMENSION (MAX_USER_STREAMS) :: USER_ANGLES

!  User azimuth variables

      INTEGER :: N_USER_RELAZMS
      REAL(FPK), DIMENSION (MAX_USER_RELAZMS) :: USER_RELAZMS

!  Number of level output choices (all)

      INTEGER :: N_LOUTPUT

!  Off-boundary level choices

      REAL(FPK), DIMENSION (MAX_LOUTPUT) :: LPARTIALS_OUTPUT

!  Layer boundary choices

      INTEGER, DIMENSION (MAX_LOUTPUT) :: LBOUNDARIES_OUTPUT


      END TYPE LRRS_Fixed_UserValues

! #####################################################################
! #####################################################################

      TYPE LRRS_Fixed_Spectral


!  Ranked wavelengths

      REAL(FPK), DIMENSION (MAX_POINTS) :: LAMBDAS_RANKED

!  Ranked fluxes

      REAL(FPK), DIMENSION (MAX_POINTS) :: FLUXES_RANKED


!  ### For monochromatic calculations ###

!  Number of monochromatic points

      INTEGER :: NPOINTS_MONO

!  Excitation wavelength

      REAL(FPK) :: LAMBDA_EXCIT

!  Excitation index

      INTEGER :: W_EXCIT


!  ### For binning calculations ###

!  Number of inner-window points

      INTEGER :: NPOINTS_INNER

!  Offset for inner-window range

      INTEGER :: OFFSET_INNER

!  Number of outer-window points

      INTEGER :: NPOINTS_OUTER

!  Bin upper and lower limits

      REAL(FPK), DIMENSION (MAX_POINTS) :: BINLOWER
      REAL(FPK), DIMENSION (MAX_POINTS) :: BINUPPER


      END TYPE LRRS_Fixed_Spectral

! #####################################################################
! #####################################################################

      TYPE LRRS_Fixed_Atmosphere


!  Height grid

      REAL(FPK), DIMENSION (0:MAX_LAYERS) :: HEIGHT_GRID

!  Input layer temperatures, must be in deg K

      REAL(FPK), DIMENSION (MAX_LAYERS) :: LAYER_TEMPERATURES

!  Input layer Air columns, should be in mol/cm^2 or [DU]

      REAL(FPK), DIMENSION (MAX_LAYERS) :: LAYER_AIRCOLUMNS

!  Rayleigh cross-sections

      REAL(FPK), DIMENSION (MAX_POINTS) :: RAYLEIGH_XSEC

!  Rayleigh depolarization ratios

      REAL(FPK), DIMENSION (MAX_POINTS) :: RAYLEIGH_DEPOL

!  Unscaled optical depth (elastic input)

      REAL(FPK), DIMENSION (MAX_LAYERS, MAX_POINTS) :: DELTAU_INPUT_UNSCALED

!  Product of unscaled single-scatter albedo and phase function
!    moments (elastic input)

      REAL(FPK), DIMENSION (MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS) :: &
        OMEGAMOMS_ELASTIC_UNSCALED


      END TYPE LRRS_Fixed_Atmosphere

! #####################################################################
! #####################################################################

      TYPE LRRS_Fixed_Surface


!  Fraction of surface treated as Lambertian

      REAL(FPK) :: LAMBERTIAN_FRACTION

!  Number of BRDF azimuth streams

      INTEGER :: NSTREAMS_BRDF

!  BRDF name

      CHARACTER (LEN=10) :: BRDF_NAMES

!  LRRS index of BRDF

      INTEGER :: WHICH_BRDF

!  BRDF weighting factor

      REAL(FPK) :: BRDF_FACTOR

!  Number of parameters associated with BRDF

      INTEGER :: BRDF_NPARS

!  BRDF parameters

      REAL(FPK), DIMENSION (MAX_BRDF_PARAMETERS) :: BRDF_PARS

!  Ranked surface albedos

      REAL(FPK), DIMENSION (MAX_POINTS) :: ALBEDOS_RANKED


      END TYPE LRRS_Fixed_Surface

! #####################################################################
! #####################################################################

      TYPE LRRS_Fixed_Inputs


      TYPE(LRRS_Fixed_Boolean)    :: Bool
      TYPE(LRRS_Fixed_Control)    :: Cont
      TYPE(LRRS_Fixed_UserValues) :: UserVal
      TYPE(LRRS_Fixed_Spectral)   :: Spect
      TYPE(LRRS_Fixed_Atmosphere) :: Atmos
      TYPE(LRRS_Fixed_Surface)    :: Surf


      END TYPE LRRS_Fixed_Inputs

! #####################################################################
! #####################################################################

      TYPE LRRS_Modified_Boolean


!  Top level control

      LOGICAL :: DO_ELASTIC_ONLY

!  Solution method (options are mutually exclusive)

      LOGICAL :: DO_CABANNES_RAMAN
      LOGICAL :: DO_ENERGY_BALANCING

!  MS mode flag

      LOGICAL :: DO_MSMODE_LRRS

!  Binning or monochromatic realization
!    (options are mutually exclusive)

      LOGICAL :: DO_BIN_REALIZATION
      LOGICAL :: DO_MONO_REALIZATION

!  User angles flag

      LOGICAL :: DO_USER_STREAMS

!  No azimuth flag

      LOGICAL :: DO_NO_AZIMUTH


      END TYPE LRRS_Modified_Boolean

! #####################################################################
! #####################################################################

      TYPE LRRS_Modified_Atmosphere


!  A height used when using the outgoing sphericity correction

      REAL(FPK) :: GEOMETRY_SPECHEIGHT


      END TYPE LRRS_Modified_Atmosphere

! #####################################################################
! #####################################################################

      TYPE LRRS_Modified_Inputs


      TYPE(LRRS_Modified_Boolean)    :: MBool
      TYPE(LRRS_Modified_Atmosphere) :: MAtmos


      END TYPE LRRS_Modified_Inputs

! #####################################################################
! #####################################################################

      PRIVATE
      PUBLIC :: LRRS_Fixed_Boolean,&
                LRRS_Fixed_Control,&
                LRRS_Fixed_UserValues,&
                LRRS_Fixed_Spectral,&
                LRRS_Fixed_Atmosphere,&
                LRRS_Fixed_Surface,&
                LRRS_Fixed_Inputs, &
                LRRS_Modified_Boolean,&
                LRRS_Modified_Atmosphere,&
                LRRS_Modified_Inputs

      END MODULE LRRS_Inputs_Def
