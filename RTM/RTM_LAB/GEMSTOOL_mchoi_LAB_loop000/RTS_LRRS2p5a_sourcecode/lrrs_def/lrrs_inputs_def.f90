
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

      MODULE LRRS_Inputs_Def_m

!  This Module contains the following LIDORT_RRS Input Structures,
!  with Intents :

!            LRRS_Fixed_Inputs    Intent(In)
!         LRRS_Modified_Inputs    Intent(InOut)

      USE LRRS_PARS_m

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
!   -- Rob mod 5/12/17 for 2p5a. NADIR flag is new, SSFULL flag renamed
!   -- Rob mod 5/12/17 for 2p5a. NADIR fand OUTGOING options are exclusive.
 
      LOGICAL :: DO_SSCORR_NADIR
      LOGICAL :: DO_SSCORR_OUTGOING
      LOGICAL :: DO_SSCORR_ALONE         ! formerly DO_SSFULL

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

!  Surface flags. removed for LRRS Version 2.5
!    Control is now in the new BRDF/SLEAVE supplements
!      LOGICAL :: DO_LAMBERTIAN_SURFACE
!      LOGICAL :: DO_WATERLEAVING
!      LOGICAL :: DO_SHADOW_EFFECT
!      LOGICAL :: DO_GLITTER_DBMS
!      LOGICAL :: DO_LAMBERTIAN_FRACTION

!  Surface control. Added, 9/9/15 for LRRS Version 2.5
!   - Modeled after the LIDORT 3.7 code

      LOGICAL :: DO_BRDF_SURFACE

!  Surface leaving Control. Added, 9/9/15 for LRRS Version 2.5
!   - Modeled after the LIDORT 3.7 code

      LOGICAL :: DO_SURFACE_LEAVING
      LOGICAL :: DO_SL_ISOTROPIC

!  Debug write flags

      LOGICAL :: DO_LRRS_WRITEINPUT
      LOGICAL :: DO_LRRS_WRITESCENARIO
      LOGICAL :: DO_LRRS_WRITEFOURIER
      LOGICAL :: DO_LRRS_WRITERESULTS

!  New SLTERM flag. @@@ RobFix 06 sep 12
!    - Replaced by new supplement controls, 9/9/15 for LRRS Version 2.5
!      LOGICAL :: DO_LRRS_SLTERM

!  Perform/display internal timing (primarily for OpenMP)
!New 06Nov17
      LOGICAL :: DO_TIMING

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
!   -- Rob mod 5/12/17 for 2p5a, Dont need this now
!      INTEGER :: NMOMENTS_INPUT

!  Number of fine layers for sscorr outgoing

      INTEGER   :: NFINELAYERS
      
!  Number of parallel threads for elastic & Raman computations
!  New 06Nov17 (primarily for OpenMP)

      INTEGER   :: NTHREADS
      
      
!  Flux factor

      REAL(FPK) :: FLUX_FACTOR

!  Solar beam input geometry

      REAL(FPK) :: SOLAR_ANGLE

!  Earth Radius

      REAL(FPK) :: EARTH_RADIUS

!  Overall accuracy for convergence criterion

      REAL(FPK) :: LRRS_ACCURACY

!  Small Number control.
!   Taylor ordering parameter, added 9/9/15. Should be set to 2 or 3

      INTEGER   :: LRRS_TAYLOR_ORDER
      REAL(FPK) :: LRRS_TAYLOR_SMALL

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

!  New Isotropic SLTERM. @@@ RobFix 06 sep 12
!      REAL(FPK), DIMENSION (MAX_POINTS) :: LRRS_SLTERM

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

!  Product of unscaled single-scatter albedo and expansion coefficients (elastic input)
!    -- Rob Mod 5/12/17 for 2p5a, change dimension to MAX_MOMENTS

      REAL(FPK), DIMENSION (MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS) :: OMEGAMOMS_ELASTIC_UNSCALED

!  Product of unscaled single-scatter albedo and Phase functions (elastic input)
!   alternative to use of the above expansion coefficient products in single-scatter codes.
!    -- Rob Mod 5/12/17, introduced for 2p5a

      REAL(fpk), dimension ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS ) :: OMEGAPHASFUNC_ELASTIC_UP
      REAL(fpk), dimension ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS ) :: OMEGAPHASFUNC_ELASTIC_DN


      END TYPE LRRS_Fixed_Atmosphere

! #####################################################################
! #####################################################################

      TYPE LRRS_Fixed_Surface

!  All the BRDF stuff has been removed for LRRS Version 2.5. 9/9/15
!    All BRDF control now in the BRDF supplement....
!  Fraction of surface treated as Lambertian
!      REAL(FPK) :: LAMBERTIAN_FRACTION
!  Number of BRDF azimuth streams
!      INTEGER :: NSTREAMS_BRDF
!  BRDF name
!      CHARACTER (LEN=10) :: BRDF_NAMES
!  LRRS index of BRDF
!      INTEGER :: WHICH_BRDF
!  BRDF weighting factor
!      REAL(FPK) :: BRDF_FACTOR
!  Number of parameters associated with BRDF
!      INTEGER :: BRDF_NPARS
!  BRDF parameters
!      REAL(FPK), DIMENSION (MAX_BRDF_PARAMETERS) :: BRDF_PARS

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

      END MODULE LRRS_Inputs_Def_m
