
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
! #   SUBROUTINE :                                              #
! #             L_RAMAN_MASTER (master)                         #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Use of Supplement-derived BRDF/SLEAVE inputs and control
!    (2) Bookkeeping improvements (use of "Only", clearer I/O specifications)
!    (3) Use of new Taylor series inputs and modules.

!   -- Rob mod 5/12/17 for 2p5a, Using PHASFUNC input now.

      MODULE lrrs_L_main_master_m

!      USE LRRS_PARS_m, Only : LDU

      PRIVATE
      PUBLIC :: L_RAMAN_MASTER

      CONTAINS

      SUBROUTINE L_RAMAN_MASTER &
        ( LRRS_FixIn, LRRS_ModIn, LRRS_Sup, & ! Inputs
          LRRS_LinIn, LRRS_LinSup,          & ! Inputs
          LRRS_Out, LRRS_LinOut )             ! Output

!  Module of dimensions and numbers

      USE LRRS_PARS_m

!  Type structure definitions

      USE LRRS_INPUTS_DEF_m
      USE LRRS_OUTPUTS_DEF_m
      USE LRRS_SUP_INOUT_DEF_m

      USE LRRS_LININPUTS_DEF_m
      USE LRRS_LINOUTPUTS_DEF_m
      USE LRRS_LIN_SUP_INOUT_DEF_m

!  Sourcecode modules, regular
!   -- Rob mod 5/12/17 for 2p5a, Introduced MISCSETUPS, Renamed SETUP_MASTER
!   -- Rob mod 5/12/17 for 2p5a, Renamed SSCORR subroutines, with DO_SSCORR_NADIR flag and facility.

      USE LRRS_IO_CHECK_m        , Only : LRRS_CHECK_INPUTS, LRRS_CHECK_INPUT_DIMS, LRRS_DERIVE_INPUTS 
      USE LRRS_MISCSETUPS_1_m    , Only : CHAPMAN_FUNCTION, ATTENUATION_SETUP_1, BEAM_MULTIPLIERS_1
      USE LRRS_GEOMETRY_1_m      , Only : MULTI_OUTGOING_ADJUSTGEOM
      USE LRRS_CONVERGE_m
      USE LRRS_WRITEMODULES_m

!      USE LRRS_WRITE_m.            Disabled Version 2.5.

!  Sourcecode modules, linearized

      USE LRRS_L_IO_CHECK_m          , Only : LRRS_L_CHECK_INPUT_DIMS
      USE LRRS_L_MISCSETUPS_1_m      , Only : LP_ATTENUATION_SETUP_1, LC_ATTENUATION_SETUP_1, LPC_BEAM_MULTIPLIERS_1
      USE LRRS_CORRECTIONS_PLUS_1_m  , Only : LRRS_SSCORR_GENERAL_BIN_PLUS_1, LRRS_SSCORR_GENERAL_MONO_PLUS_1

      USE LRRS_L_FOURIER_MASTER_1_m
      USE LRRS_L_RRSOPTICAL_MASTER_m  !  Formerly LRRS_L_SETUP_MASTER_m
      USE LRRS_L_CONVERGE_m
      USE LRRS_L_WRITEMODULES_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

      TYPE(LRRS_Fixed_Inputs), INTENT(IN)        :: LRRS_FixIn
      TYPE(LRRS_Modified_Inputs), INTENT(INOUT)  :: LRRS_ModIn
      TYPE(LRRS_Sup_InOut), INTENT (INOUT)       :: LRRS_Sup
      TYPE(LRRS_Outputs), INTENT(INOUT)          :: LRRS_Out

      TYPE(LRRS_LinInputs), INTENT(IN)           :: LRRS_LinIn
      TYPE(LRRS_LinSup_InOut), INTENT (INOUT)    :: LRRS_LinSup
      TYPE(LRRS_LinOutputs), INTENT(INOUT)       :: LRRS_LinOut

!  Local model input-related variables
!  ===================================

!  Control variables
!    ** BRDF flags removed, general surface flags added for Version 2.5. 10/12/15.
!   -- Rob mod 5/12/17 for 2p5a, added SSCORR_NADIR, renamed SSFULL --> SSCORR_ALONE.
!   -- 10/10/18. Added DO_TIMING flag to this list (for OMP)

      LOGICAL :: &
         DO_RRS_OVERALL,      DO_DIRECTRRS_ONLY,    DO_ELASTIC_ONLY, &
         DO_CABANNES_RAMAN,   DO_ENERGY_BALANCING,  DO_MSMODE_LRRS, &
         DO_BIN_REALIZATION,  DO_MONO_REALIZATION,  DO_DIRECT_BEAM, &
         DO_SSCORR_NADIR,     DO_SSCORR_OUTGOING,   DO_SSCORR_ALONE, DO_MOLECSCAT_ONLY, &
         DO_PLANE_PARALLEL,   DO_DELTAM_SCALING,    DO_DOUBLE_CONVTEST, &
         DO_UPWELLING,        DO_DNWELLING,         DO_USER_STREAMS, &
         DO_NO_AZIMUTH,       DO_LBOUNDARIES,       DO_MVOUT_ONLY, &
         DO_ADDITIONAL_MVOUT, DO_BRDF_SURFACE,      DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, &
         DO_LRRS_WRITESCENARIO, DO_LRRS_WRITERESULTS, &
         DO_LRRS_WRITEINPUT,    DO_LRRS_WRITEFOURIER, DO_TIMING

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  number of layers

      INTEGER :: NLAYERS

!  Number of input phase function Legendre moments  THIS IS A BASIC INPUT THAT USER MUST PROVIDE
!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT !!!!!!!!!!!!!!!
!      INTEGER :: NMOMENTS_INPUT

!  Number of fine layers for sscorr outgoing

      INTEGER :: NFINELAYERS

!  Number of parallel threads for elastic & Raman computations
!  (primarily for OpenMP). Introduced 10/10/18 for Linearized code.

      INTEGER   :: NTHREADS

!  Solar beam input geometry

      REAL(FPK) :: SOLAR_ANGLE

!  User stream variables

      INTEGER :: N_USER_STREAMS
      REAL(FPK) :: USER_ANGLES  ( MAX_USER_STREAMS )

!  User azimuth variables

      INTEGER :: N_USER_RELAZMS
      REAL(FPK) :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  Level/layer output control
!    N_LOUTPUT = number of level output choices (all)
!    LBOUNDARIES_OUTPUT = layer boundary choices
!    LPARTIALS_OUTPUT = off-boundary level choices

      INTEGER   :: N_LOUTPUT
      REAL(FPK) :: LPARTIALS_OUTPUT   ( MAX_LOUTPUT )
      INTEGER   :: LBOUNDARIES_OUTPUT ( MAX_LOUTPUT )

!  Progress control (If set to 0, no screen output)

      INTEGER :: PROGRESS

!  BINNING: outer/inner wavelength range, and offset for inner range
!    Depends on the choice of solar spectrum

      INTEGER :: OFFSET_INNER
      INTEGER :: NPOINTS_INNER
      INTEGER :: NPOINTS_OUTER

!  MONOCHROMATIC: Excitation wavelength and index. NPOINTS_MONO = 234

      INTEGER   :: NPOINTS_MONO
      REAL(FPK) :: LAMBDA_EXCIT
      INTEGER   :: W_EXCIT

!  Surface leaving prototype for Version 2.4, Now commented out, 9/11/15
!     Add SLTERM flag     . @@@ Rob Fix 06 Sep 12.
!     New Isotropic SLTERM. @@@ Rob Fix 06 sep 12
!      LOGICAL    :: DO_LRRS_SLTERM
!      REAL(FPK)  :: LRRS_SLTERM (MAX_POINTS)

!  Overall accuracy for convergence criterion

      REAL(FPK) :: LRRS_ACCURACY

!  Small number control. Expanded for Version 2.5, 9/11/15

      INTEGER   :: TAYLOR_ORDER
      REAL(FPK) :: TAYLOR_SMALL

!  Flux factor

      REAL(FPK) :: FLUX_FACTOR

!  Earth Radius

      REAL(FPK) :: EARTH_RADIUS

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007 (LIDORT)
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 18 March  2011 (LRRS)

!    This is only required when the outgoing sphericity correction is
!    in operation. Otherwise, the regular pseudo-spherical correction
!    (with or without an exact single-scatter correction) uses the same
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

!  FOLLOWING CODE IS COMMENTED OUT, VERSION 2.5.
!     Lambertian fraction, Number of BRDF azimuth streams
!  BRDF parameters. Index choice of BRDF function
!     For Cox-Munk, first parameter  = Wind Speed
!                 second parameter = refractive index
!                 third parameter  = 0.0d
!      REAL(FPK) :: LAMBERTIAN_FRACTION
!      INTEGER :: NSTREAMS_BRDF
!      INTEGER :: WHICH_BRDF
!      INTEGER :: BRDF_NPARS
!      REAL(FPK) :: BRDF_PARS ( MAX_BRDF_PARAMETERS)
!      REAL(FPK) :: BRDF_FACTOR
!      CHARACTER (LEN=10) :: BRDF_NAMES

!  filenames for output

      CHARACTER (LEN=60) :: DEBUG_FILENAMES(4)

!  Atmospheric quantities
!  ----------------------

!  Height grid

      REAL(FPK) :: HEIGHT_GRID ( 0:MAX_LAYERS )

!  Linearized multilayer Height inputs. ******* New, 4.19.2010 ********
!      REAL(FPK) :: LH_HEIGHT_GRID ( MAX_ATMOSWFS, 0:MAX_LAYERS )

!  Input layer temperatures, must be in deg K

      REAL(FPK) :: LAYER_TEMPERATURES ( MAX_LAYERS )

!  Input layer Air columns, should be in mol/cm^2 or [DU]

      REAL(FPK) :: LAYER_AIRCOLUMNS    ( MAX_LAYERS )

!  Wavelengths and Fluxes
!  ----------------------

!  These are defined on the outer grid for the binning realization
!  These are defined on 234 points for the monochromatic

!  Depends on the choice of solar spectrum

      REAL(FPK) :: LAMBDAS_RANKED  ( MAX_POINTS )
      REAL(FPK) :: FLUXES_RANKED   ( MAX_POINTS )

!  Surface optical
!  ---------------

!  Flags for the multiple-point surface options

      LOGICAL :: DO_BRDF_Wav1
      LOGICAL :: DO_SLEAVE_Wav1

!  Lambertian albedos. New Version 2.5, 9/11/15

      REAL(fpk) :: ALBEDOS_RANKED ( MAX_POINTS )

!  Exact (direct bounce) BRDF. Version 2.5 9/11/15

      REAL(fpk) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

!  Fourier components of BRDF, in the following order
!    incident solar direction,    reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar direction,    reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk) :: BRDF_F_0      ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk) :: BRDF_F        ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk) :: USER_BRDF_F_0 ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )
      REAL(fpk) :: USER_BRDF_F   ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Isotropic Surface leaving term (if flag set). New Version 2.5, 9/11/15
!  Exact Surface-Leaving term. New Version 2.5, 9/11/15

      REAL(fpk) :: SLTERM_ISOTROPIC ( MAX_POINTS )
      REAL(fpk) :: SLTERM_USERANGLES (MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

!  Fourier components of Surface-leaving terms:
!    solar direction, SL-transmitted quadrature streams
!    solar direction, SL-transmitted user streams

      REAL(fpk) ::  SLTERM_F_0      ( 0:MAX_MOMENTS, MAX_STREAMS,      MAX_POINTS )
      REAL(fpk) ::  USER_SLTERM_F_0 ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )

!  Basic Rayleigh data
!  -------------------

!     Rayleigh Cross-sections and depolarization ratios

      REAL(FPK) :: RAYLEIGH_XSEC  ( MAX_POINTS )
      REAL(FPK) :: RAYLEIGH_DEPOL ( MAX_POINTS )

!  For the binning realization
!  ---------------------------

!  Bin upper and lower limits

      REAL(FPK) :: BINLOWER ( MAX_POINTS )
      REAL(FPK) :: BINUPPER ( MAX_POINTS )

!  Unscaled Optical properties, elastic
!  ------------------------------------

!  These must be defined on the outer wavelength grid (binning)
!  Unscaled quantities, Elastic input
!   -- Rob mod 5/12/17 for 2p5a, OMEGAMOMS_ELASTIC_UNSCALED replacing dimension MAX_MOMENTS_INPUT

      REAL(FPK) :: DELTAU_INPUT_UNSCALED      ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: OMEGAMOMS_ELASTIC_UNSCALED ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!   -- Rob mod 5/12/17 for 2p5a, Introducing the phase function product input

      REAL(FPK) :: OMEGAPHASFUNC_ELASTIC_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: OMEGAPHASFUNC_ELASTIC_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Linearized input-related variables
!  ==================================

!  Linearization control
!  ---------------------

!  Profile/column/surface linearization control inputs

      LOGICAL :: DO_PROFILE_WFS
      LOGICAL :: DO_COLUMN_WFS
      LOGICAL :: DO_SURFACE_WFS
      LOGICAL :: DO_SLEAVE_WFS

!  Flag for generating normalized weighting functions

      LOGICAL :: DO_NORMALIZED_WFS

!  Flag for generating air column weighting functions

      LOGICAL :: DO_AIRPROFILE_WFS

!  Flag for generating temperature profile weighting functions

      LOGICAL :: DO_TEMPPROFILE_WFS

!  Flag for generating temperature shift weighting functions

      LOGICAL :: DO_TEMPSHIFT_WF

!  H-level Linearization control ******* New, 4.19.2010 ********
!      LOGICAL :: DO_HLEVEL_LINEARIZATION
!      LOGICAL :: LHMASK_COLUMN_WFS ( MAX_ATMOSWFS )

!  Linearization atmosphere control
!  --------------------------------

!  Total number of column weighting functions

      INTEGER :: N_TOTALCOLUMN_WFS

!  Control for atmospheric linearizations, layer by layer

      LOGICAL :: LAYER_VARY_FLAG (MAX_LAYERS)
      INTEGER :: LAYER_VARY_NUMBER (MAX_LAYERS)

!  Input layer air columns, should be in mol/cm^2 or [DU]

      REAL(FPK) :: LAYER_AIRCOLUMNS_DT (MAX_LAYERS)

!  Input layer temperatures (unshifted), must be in deg K

      REAL(FPK) :: TEMPERATURES_UNSHIFTED (MAX_LAYERS)

!  Linearized unscaled optical properties, elastic input
!   -- Rob mod 5/12/17 for 2p5a, OMEGAMOMS_ELASTIC_UNSCALED replacing dimension MAX_MOMENTS_INPUT

      REAL(FPK) :: L_DELTAU_INPUT_UNSCALED      ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L_OMEGAMOMS_ELASTIC_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!   -- Rob mod 5/12/17 for 2p5a, Introducing the phase function product input

      REAL(FPK) :: L_OMEGAPHASFUNC_ELASTIC_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: L_OMEGAPHASFUNC_ELASTIC_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Linearized BRDF
!  ---------------

!  Number of surface Jacobians

      INTEGER   :: N_SURFACE_WFS

!  Linearized Exact (direct bounce) BRDF (same all threads)

      REAL(fpk) :: LS_EXACTDB_BRDFUNC ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

!  Linearized Fourier components of BRDF, in the following order
!    incident solar direction,    reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar direction,    reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk) :: LS_BRDF_F_0      ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk) :: LS_BRDF_F        ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )
      REAL(fpk) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Linearized SLEAVE
!  -----------------

      INTEGER   :: N_SLEAVE_WFS

      REAL(fpk) :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAX_POINTS )
      REAL(fpk) :: LSSL_SLTERM_USERANGLES ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )
      REAL(fpk) :: LSSL_SLTERM_F_0        ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk) :: LSSL_USER_SLTERM_F_0   ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )

!  ROUTINE OUTPUT
!  ==============

!  ELASTIC FIELD
!  -------------

!   Fourier-summed output

      REAL(FPK) :: ELASTIC_UP ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: ELASTIC_DN ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  single scatter results

      REAL(FPK) :: ELASTIC_SS_UP ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: ELASTIC_SS_DN ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Mean-value output

      REAL(FPK) :: MEAN_ELASTIC_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: MEAN_ELASTIC_DN ( MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK) :: FLUX_ELASTIC_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: FLUX_ELASTIC_DN ( MAX_LOUTPUT, MAX_POINTS )

!  Fourier summed output, Atmospheric Column Jacobians

      REAL(FPK) :: LC_ELASTIC_UP ( MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: LC_ELASTIC_DN ( MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Fourier summed output, Atmospheric Profile Jacobians

      REAL(FPK) :: LP_ELASTIC_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: LP_ELASTIC_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Fourier summed output, Surface Jacobians

      REAL(FPK) :: LS_ELASTIC_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: LS_ELASTIC_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Linearized Single scatter results

      REAL(FPK) :: LC_ELASTIC_SS_UP ( MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: LC_ELASTIC_SS_DN ( MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK) :: LP_ELASTIC_SS_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: LP_ELASTIC_SS_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK) :: LS_ELASTIC_SS_UP  ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Linearized Mean-value output

      REAL(FPK) :: L_MEAN_ELASTIC_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: L_MEAN_ELASTIC_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK) :: L_FLUX_ELASTIC_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: L_FLUX_ELASTIC_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK) :: LS_MEAN_ELASTIC_UP ( MAX_SURFACEWFS, MAX_LOUTPUT,MAX_POINTS )
      REAL(FPK) :: LS_MEAN_ELASTIC_DN ( MAX_SURFACEWFS, MAX_LOUTPUT,MAX_POINTS )
      REAL(FPK) :: LS_FLUX_ELASTIC_UP ( MAX_SURFACEWFS, MAX_LOUTPUT,MAX_POINTS )
      REAL(FPK) :: LS_FLUX_ELASTIC_DN ( MAX_SURFACEWFS, MAX_LOUTPUT,MAX_POINTS )

!  RAMAN FIELD
!  -----------

!  Radiance including inelastic scattering

      REAL(FPK) :: RAMAN_UP ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS)
      REAL(FPK) :: RAMAN_DN ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS)

!  single scatter results

      REAL(FPK) :: RAMAN_SS_UP ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: RAMAN_SS_DN ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Mean-value output

      REAL(FPK) :: MEAN_RAMAN_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: MEAN_RAMAN_DN ( MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK) :: FLUX_RAMAN_UP ( MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: FLUX_RAMAN_DN ( MAX_LOUTPUT, MAX_POINTS )

!  Fourier summed output, Atmospheric Column Jacobians

      REAL(FPK) :: LC_RAMAN_UP ( MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: LC_RAMAN_DN ( MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Fourier summed output, Atmospheric Profile Jacobians

      REAL(FPK) :: LP_RAMAN_UP ( MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK) :: LP_RAMAN_DN ( MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Fourier summed output, Surface Jacobians

      REAL(FPK) :: LS_RAMAN_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: LS_RAMAN_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Linearized Single scatter results

      REAL(FPK) :: LC_RAMAN_SS_UP ( MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: LC_RAMAN_SS_DN ( MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK) :: LP_RAMAN_SS_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: LP_RAMAN_SS_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK) :: LS_RAMAN_SS_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Linearized Mean-value output

      REAL(FPK) :: L_MEAN_RAMAN_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: L_MEAN_RAMAN_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK) :: L_FLUX_RAMAN_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: L_FLUX_RAMAN_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS )

      REAL(FPK) :: LS_MEAN_RAMAN_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: LS_MEAN_RAMAN_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: LS_FLUX_RAMAN_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )
      REAL(FPK) :: LS_FLUX_RAMAN_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS )

!  Number of Fourier terms used

      INTEGER :: FOURIER_SAVED

!  Error handling
!  --------------

!  Inputs

      INTEGER             :: STATUS_INPUTCHECK
      INTEGER             :: NCHECKMESSAGES
      CHARACTER(LEN=120)  :: CHECKMESSAGES (0:MAX_MESSAGES)
      CHARACTER(LEN=120)  :: ACTIONS (0:MAX_MESSAGES)

!  Main

      INTEGER :: STATUS
      INTEGER :: N_MESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( MAX_MESSAGES)

!  Local variables
!  ===============

!  bookkeping variables
!  --------------------

!  Number of moments

      INTEGER   :: NMOMENTS

!  Number of convergence tests

      INTEGER   :: N_CONVTESTS

!  number of output streams (quadrature and user)

      INTEGER   :: N_OUT_STREAMS

!  Solar beam input Cosine

      REAL(FPK) :: COS_SZA

!  Discrete ordinate quadrature

      REAL(FPK) :: QUAD_STREAMS ( MAX_STREAMS )
      REAL(FPK) :: QUAD_WEIGHTS ( MAX_STREAMS )
      REAL(FPK) :: QUAD_STRMWGT ( MAX_STREAMS )

!  User stream variables

      REAL(FPK) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  Number of Geometries, Geometry offsetting
!    --Added by R. Spurr, 04 August 2009

      INTEGER :: N_GEOMETRIES
      INTEGER :: GEOM_OFFSETS ( MAX_USER_STREAMS )

!  Layer masks for doing integrated source terms

      LOGICAL :: STERM_MASK_UP ( MAX_LAYERS )
      LOGICAL :: STERM_MASK_DN ( MAX_LAYERS )

      INTEGER :: LEVELMASK_UP ( MAX_LOUTPUT )
      INTEGER :: LEVELMASK_DN ( MAX_LOUTPUT )

!  Fourier output
!  --------------

!  Elastic field
!  -------------

!  Radiances

      REAL(FPK) :: ELASTIC_F_UP ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK) :: ELASTIC_F_DN ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  For the atmospheric weighting functions, post-processed values

      REAL(FPK) :: L_ELASTIC_F_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK) :: L_ELASTIC_F_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  For the surface weighting functions, post-processed values

      REAL(FPK) :: LS_ELASTIC_F_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK) :: LS_ELASTIC_F_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  RAMAN FIELD
!  -----------

!  Raman Fourier output, not saved

      REAL(FPK) :: RAMAN_F_UP ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK) :: RAMAN_F_DN ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  For the atmospheric weighting functions, post-processed values

      REAL(FPK) :: L_RAMAN_F_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK) :: L_RAMAN_F_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  For the surface weighting functions, post-processed values

      REAL(FPK) :: LS_RAMAN_F_UP ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )
      REAL(FPK) :: LS_RAMAN_F_DN ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  Miscellaneous Setups (to do with Beam)
!  --------------------------------------

!  Chapman factors

      REAL(FPK) :: CHAPMAN ( MAX_LAYERS, MAX_LAYERS )

!  LH-Linearizated Chapman factors  ****** New 4.19.2010. *******
!      REAL(FPK) :: L_CHAPMAN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_LAYERS )

!  Solar beam transmittances, average secant factors

      INTEGER ::   BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: BEAM_ETRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )

!  User stream transmittances, whole layers

      REAL(FPK) :: SAVE_TRANS_USERM ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Beam Multipliers for use in the SSCORR_NADIR single-scattering correction
!   -- Rob mod 5/12/17 for 2p5a, Newly introduced as output from BEAM_MULTIPLIERS_1

      REAL(FPK) :: SAVE_BEAMMULT_UP ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: SAVE_BEAMMULT_DN ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Solar beam transmittances, average secant factors

      REAL(FPK) :: L_BEAM_ITRANS   ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK) :: L_BEAM_AVSECANT ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK) :: L_BEAM_ETRANS   ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK) :: L_BEAM_DTRANS   ( MAX_ATMOSWFS, MAX_LAYERS_NK, MAX_POINTS )

!  User stream transmittances, whole layers

      REAL(FPK) :: L_SAVE_TRANS_USERM ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Beam Multipliers for use in the SSCORR_NADIR single-scattering correction
!   -- Rob mod 5/12/17 for 2p5a, Newly introduced as output from L_BEAM_MULTIPLIERS_1

      REAL(FPK) :: L_SAVE_BEAMMULT_UP ( MAX_ATMOSWFS, MAX_USER_STREAMS, 0:MAX_LAYERS_NK, MAX_POINTS )
      REAL(FPK) :: L_SAVE_BEAMMULT_DN ( MAX_ATMOSWFS, MAX_USER_STREAMS, 0:MAX_LAYERS_NK, MAX_POINTS )

!  Optical setup variables, including RRS
!  --------------------------------------

!  Basic input quantities for elastic scattering
!   -- Rob mod 5/12/17 for 2p5a, OMEGA_UNSCALED added

      REAL(FPK) :: DELTAU_VERT_INPUT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: OMEGA_UNSCALED    ( MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK) :: OMEGAMOMS_ELASTIC ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Truncation factors (zero if no deltam scaling)

      REAL(FPK) :: TRUNC_FACTORS ( MAX_LAYERS, MAX_POINTS )

!  Linearized values 
!   -- Rob mod 5/12/17 for 2p5a, L_OMEGA_UNSCALED added

      REAL(FPK) :: L_TRUNC_FACTORS     ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(FPK) :: L_OMEGA_UNSCALED    ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )

!  Derived Elastic scattering (Internal variable)

      REAL(FPK) :: L_OMEGAMOMS_ELASTIC ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Bin Mapping quantities (internal derivation)

      INTEGER   :: N_RRSBINS  ( MAX_POINTS )
      INTEGER   :: BINMAP ( MAX_BINS, MAX_POINTS )

!  Chance/Spurr Ring Spectrum. DIAGNOSTIC ONLY.

      REAL(FPK) :: RINGSPEC_1 ( MAX_POINTS )

!  Unscaled Raman optical properties
!   -- Rob mod 5/12/17 for 2p5a, Introducing the phase function product instead of OMEGAMOMS_CABANNES_UNSCALED
!      REAL(FPK) :: OMEGAMOMS_CABANNES_UNSCALED ( MAX_LAYERS, 0:MAX_MOMENTS_INPUT, MAX_POINTS )

      REAL(FPK) :: OMEGAPHASFUNC_CABANNES_UP ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: OMEGAPHASFUNC_CABANNES_DN ( MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK) :: OMEGAMOMS_RRSLOSS_UNSCALED  ( MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK) :: OMEGAMOMS_RRSBIN_UNSCALED   ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )
      REAL(FPK) :: OMEGAMOMS_RRSGAIN_UNSCALED  ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Linearized

      REAL(FPK) :: L_OMEGAPHASFUNC_CABANNES_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: L_OMEGAPHASFUNC_CABANNES_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK) :: L_OMEGAMOMS_RRSLOSS_UNSCALED  ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK) :: L_OMEGAMOMS_RRSBIN_UNSCALED   ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )
      REAL(FPK) :: L_OMEGAMOMS_RRSGAIN_UNSCALED  ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived Cabannes scattering (Internal variable)

      REAL(FPK) :: OMEGAMOMS_CABANNES ( MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  Derived rotational raman scattering (Internal variable). Loss term
!    Only required for the Energy-balance approximation.

      REAL(FPK) :: OMEGAMOMS_RRSLOSS ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Derived rotational Raman scattering (Internal variable). BIN and MONO

      REAL(FPK) :: OMEGAMOMS_RRSBIN   ( MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )
      REAL(FPK) :: OMEGAMOMS_RRSGAIN  ( MAX_LAYERS, 0:2, MAX_POINTS )

!  Linearizations

      REAL(FPK) :: L_OMEGAMOMS_CABANNES ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )
      REAL(FPK) :: L_OMEGAMOMS_RRSLOSS  ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )
      REAL(FPK) :: L_OMEGAMOMS_RRSBIN   ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_BINS, MAX_POINTS )
      REAL(FPK) :: L_OMEGAMOMS_RRSGAIN  ( MAX_ATMOSWFS, MAX_LAYERS, 0:2, MAX_POINTS )

!  Help variables
!  --------------

!  general SSCORR flag
!   -- Rob mod 5/12/17 for 2p5a, SSCORR_GENERAL = SSCORR_NADIR or SSCORR_OUTGOING

      LOGICAL :: DO_SSCORR_GENERAL

!  Fourier number

      INTEGER :: FOURIER
      INTEGER :: N_FOURIER_COMPONENTS

!  Local flags

      LOGICAL :: ADJUST_SURFACE
      LOGICAL :: LOCAL_DO_NO_AZIMUTH
      LOGICAL :: SAVE_DO_NO_AZIMUTH
      LOGICAL :: LOCAL_ITERATION

!  Linearization bookkeeping storage

      INTEGER :: NKSTORAGE   ( MAX_LAYERS, 0:MAX_LAYERS )
      INTEGER :: NKSTORAGE2  ( MAX_LAYERS, 0:MAX_LAYERS )

!  Total number of surface Jacobians (both surface and surface-leaving)

      INTEGER :: N_TOTALSURFACE_WFS

!  Modified eradius

      REAL(FPK) :: MODIFIED_ERADIUS

!  Adjusted geometries. New for Version 2.3, March 18, 2011.

      REAL(FPK) :: USER_ANGLES_ADJUST  ( MAX_USER_STREAMS )
      REAL(FPK) :: SOLAR_ANGLES_ADJUST ( MAX_USER_STREAMS, MAX_USER_RELAZMS )
      REAL(FPK) :: USER_RELAZMS_ADJUST ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Scattering angles
!   -- Rob mod 5/12/17 for 2p5a, added

      REAL(FPK)  :: COSSCAT_UP ( MAX_GEOMETRIES )
      REAL(FPK)  :: COSSCAT_DN ( MAX_GEOMETRIES )

!  Other variables

      INTEGER   :: I, UA, UM, V, N, K, NK, N2, LOCAL_PROGRESS
      INTEGER   :: NPOINTS_LOCAL, NPOINTS_TOTAL
      INTEGER   :: FUNIT, RUNIT, IUNIT, SUNIT
      INTEGER   :: TESTCONV, LOCAL_N_USERAZM
      INTEGER   :: STATUS_SUB
      LOGICAL   :: DO_ATMOS_WFS

!      LOGICAL   :: DO_GRIDLEVEL_WFS   ! Disabled 10 March 2011

!  Geometry variables
!   -- Rob mod 5/12/17 for 2p5a, added

      REAL(FPK) :: ALPHA_BOA, THETA_BOA, PHI_BOA
      REAL(FPK) :: SALPHA_BOA, STHETA_BOA, CALPHA_BOA, CTHETA_BOA, CPHI_BOA

!  Fourier cosine arguments

      REAL(FPK) :: AZM_ARGUMENT, DFC
      REAL(FPK) :: AZMFAC ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Exception handling variables

      CHARACTER (LEN=3)   :: C3
      CHARACTER (LEN=120) :: MESSAGE
      CHARACTER (LEN=120) :: MESSAGE_SUB
      CHARACTER (LEN=120) :: POINT_TRACE
      LOGICAL             :: FAIL

!  Test variables

      LOGICAL         :: DO_DEBUG_INPUT=.FALSE.
!      LOGICAL         :: DO_DEBUG_INPUT=.TRUE.

!  OMP Variables. Added 10/10/18
!  -------------

!  Fourier times

      REAL      :: OMP_Elastic_Time (0:MAX_MOMENTS), OMP_Raman_Time (0:MAX_MOMENTS)

!mick mod 7/20/2016 - added initializations and input dimensions checking

!  Initialize some variables
!  -------------------------

!  Main status

      STATUS_INPUTCHECK  = LRRS_SUCCESS
      STATUS             = LRRS_SUCCESS

!  Input checks

      NCHECKMESSAGES = 0
      CHECKMESSAGES  = ' '
      ACTIONS        = ' '

!  Check input dimensions
!  ----------------------

!  Regular

      CALL LRRS_CHECK_INPUT_DIMS &
        ( LRRS_FixIn, LRRS_ModIn,                            & !Inputs
          STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS ) !Outputs

!  Exception handling

      IF ( STATUS_SUB .EQ. LRRS_SERIOUS ) THEN
        STATUS_INPUTCHECK = LRRS_SERIOUS
        LRRS_Out%Status%STATUS_INPUTCHECK = STATUS_INPUTCHECK
        LRRS_Out%Status%NCHECKMESSAGES    = NCHECKMESSAGES
        LRRS_Out%Status%CHECKMESSAGES     = CHECKMESSAGES
        LRRS_Out%Status%ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  Linearized

      CALL LRRS_L_CHECK_INPUT_DIMS &
        ( LRRS_LinIn,                                        & !Inputs
          STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS ) !Outputs

!  Exception handling

      IF ( STATUS_SUB .EQ. LRRS_SERIOUS ) THEN
        STATUS_INPUTCHECK = LRRS_SERIOUS
        LRRS_Out%Status%STATUS_INPUTCHECK = STATUS_INPUTCHECK
        LRRS_Out%Status%NCHECKMESSAGES    = NCHECKMESSAGES
        LRRS_Out%Status%CHECKMESSAGES     = CHECKMESSAGES
        LRRS_Out%Status%ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  Standard Inputs
!  ---------------

!  Fixed Boolean Inputs
!   -- Rob mod 5/12/17 for 2p5a, added SSCORR_NADIR, renamed SSFULL --> SSCORR_ALONE

      DO_RRS_OVERALL         = LRRS_FixIn%Bool%DO_RRS_OVERALL
      DO_DIRECTRRS_ONLY      = LRRS_FixIn%Bool%DO_DIRECTRRS_ONLY
      DO_DIRECT_BEAM         = LRRS_FixIn%Bool%DO_DIRECT_BEAM

      DO_SSCORR_NADIR        = LRRS_FixIn%Bool%DO_SSCORR_NADIR
      DO_SSCORR_OUTGOING     = LRRS_FixIn%Bool%DO_SSCORR_OUTGOING
      DO_SSCORR_ALONE        = LRRS_FixIn%Bool%DO_SSCORR_ALONE

      DO_MOLECSCAT_ONLY      = LRRS_FixIn%Bool%DO_MOLECSCAT_ONLY
      DO_PLANE_PARALLEL      = LRRS_FixIn%Bool%DO_PLANE_PARALLEL
      DO_DELTAM_SCALING      = LRRS_FixIn%Bool%DO_DELTAM_SCALING
      DO_DOUBLE_CONVTEST     = LRRS_FixIn%Bool%DO_DOUBLE_CONVTEST
      DO_UPWELLING           = LRRS_FixIn%Bool%DO_UPWELLING
      DO_DNWELLING           = LRRS_FixIn%Bool%DO_DNWELLING
      DO_LBOUNDARIES         = LRRS_FixIn%Bool%DO_LBOUNDARIES
      DO_MVOUT_ONLY          = LRRS_FixIn%Bool%DO_MVOUT_ONLY
      DO_ADDITIONAL_MVOUT    = LRRS_FixIn%Bool%DO_ADDITIONAL_MVOUT

!  general SSCORR flag
!   -- Rob mod 5/12/17 for 2p5a, SSCORR_GENERAL = SSCORR_NADIR or SSCORR_OUTGOING

      DO_SSCORR_GENERAL      = DO_SSCORR_NADIR .or. DO_SSCORR_OUTGOING

!  Surface variables removed. Version 2.5 9/11/15
!      DO_LAMBERTIAN_SURFACE  = LRRS_FixIn%Bool%DO_LAMBERTIAN_SURFACE
!      DO_WATERLEAVING        = LRRS_FixIn%Bool%DO_WATERLEAVING
!      DO_SHADOW_EFFECT       = LRRS_FixIn%Bool%DO_SHADOW_EFFECT
!      DO_GLITTER_DBMS        = LRRS_FixIn%Bool%DO_GLITTER_DBMS
!      DO_LAMBERTIAN_FRACTION = LRRS_FixIn%Bool%DO_LAMBERTIAN_FRACTION

!     Surface variables added, Version 2.5
      DO_BRDF_SURFACE    = LRRS_FixIn%Bool%DO_BRDF_SURFACE
      DO_SURFACE_LEAVING = LRRS_FixIn%Bool%DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC    = LRRS_FixIn%Bool%DO_SL_ISOTROPIC

      DO_LRRS_WRITEINPUT     = LRRS_FixIn%Bool%DO_LRRS_WRITEINPUT
      DO_LRRS_WRITESCENARIO  = LRRS_FixIn%Bool%DO_LRRS_WRITESCENARIO
      DO_LRRS_WRITEFOURIER   = LRRS_FixIn%Bool%DO_LRRS_WRITEFOURIER
      DO_LRRS_WRITERESULTS   = LRRS_FixIn%Bool%DO_LRRS_WRITERESULTS

!  Timing variable added, Version 2.5a with OMP. 10/10/18
      DO_TIMING              = LRRS_FixIn%Bool%DO_TIMING

!  Modified Boolean Inputs

      DO_ELASTIC_ONLY     = LRRS_ModIn%MBool%DO_ELASTIC_ONLY
      DO_CABANNES_RAMAN   = LRRS_ModIn%MBool%DO_CABANNES_RAMAN
      DO_ENERGY_BALANCING = LRRS_ModIn%MBool%DO_ENERGY_BALANCING
      DO_MSMODE_LRRS      = LRRS_ModIn%MBool%DO_MSMODE_LRRS
      DO_BIN_REALIZATION  = LRRS_ModIn%MBool%DO_BIN_REALIZATION
      DO_MONO_REALIZATION = LRRS_ModIn%MBool%DO_MONO_REALIZATION
      DO_USER_STREAMS     = LRRS_ModIn%MBool%DO_USER_STREAMS
      DO_NO_AZIMUTH       = LRRS_ModIn%MBool%DO_NO_AZIMUTH

!  Fixed Control Inputs
!   -- Rob mod 5/12/17  for 2p5a, removed NMOMENTS_INPUT
!   -- Rob mod 10/10/17 for 2p5a, added NTHREADS

      PROGRESS        = LRRS_FixIn%Cont%PROGRESS
      NLAYERS         = LRRS_FixIn%Cont%NLAYERS
      NSTREAMS        = LRRS_FixIn%Cont%NSTREAMS
!      NMOMENTS_INPUT  = LRRS_FixIn%Cont%NMOMENTS_INPUT
      NFINELAYERS     = LRRS_FixIn%Cont%NFINELAYERS
      NTHREADS        = LRRS_FixIn%Cont%NTHREADS
      FLUX_FACTOR     = LRRS_FixIn%Cont%FLUX_FACTOR
      SOLAR_ANGLE     = LRRS_FixIn%Cont%SOLAR_ANGLE
      EARTH_RADIUS    = LRRS_FixIn%Cont%EARTH_RADIUS
      LRRS_ACCURACY   = LRRS_FixIn%Cont%LRRS_ACCURACY

!  Taylor variables new for version 2.5........

      TAYLOR_ORDER    = LRRS_FixIn%Cont%LRRS_TAYLOR_ORDER
      TAYLOR_SMALL    = LRRS_FixIn%Cont%LRRS_TAYLOR_SMALL

      DEBUG_FILENAMES = LRRS_FixIn%Cont%DEBUG_FILENAMES

!  Fixed User-Value Inputs

      N_USER_STREAMS     = LRRS_FixIn%UserVal%N_USER_STREAMS
      USER_ANGLES        = LRRS_FixIn%UserVal%USER_ANGLES
      N_USER_RELAZMS     = LRRS_FixIn%UserVal%N_USER_RELAZMS
      USER_RELAZMS       = LRRS_FixIn%UserVal%USER_RELAZMS
      N_LOUTPUT          = LRRS_FixIn%UserVal%N_LOUTPUT
      LBOUNDARIES_OUTPUT = LRRS_FixIn%UserVal%LBOUNDARIES_OUTPUT  !not active yet (7/20/2016)
      LPARTIALS_OUTPUT   = LRRS_FixIn%UserVal%LPARTIALS_OUTPUT

!  Fixed Spectral Inputs

      LAMBDAS_RANKED = LRRS_FixIn%Spect%LAMBDAS_RANKED
      FLUXES_RANKED  = LRRS_FixIn%Spect%FLUXES_RANKED

      IF ( DO_MONO_REALIZATION ) THEN
        !For monochromatic calculations:
        NPOINTS_MONO   = LRRS_FixIn%Spect%NPOINTS_MONO
        LAMBDA_EXCIT   = LRRS_FixIn%Spect%LAMBDA_EXCIT
        W_EXCIT        = LRRS_FixIn%Spect%W_EXCIT
      END IF

      IF ( DO_BIN_REALIZATION ) THEN
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

      DELTAU_INPUT_UNSCALED      = LRRS_FixIn%Atmos%DELTAU_INPUT_UNSCALED
      OMEGAMOMS_ELASTIC_UNSCALED = LRRS_FixIn%Atmos%OMEGAMOMS_ELASTIC_UNSCALED

!   -- Rob mod 5/12/17 for 2p5a, Setting local OMEGA_UNSCLAED

      OMEGA_UNSCALED(:,:)        = OMEGAMOMS_ELASTIC_UNSCALED(:,0,:)

!   -- Rob mod 5/12/17 for 2p5a, Copy Phase function product input

      OMEGAPHASFUNC_ELASTIC_UP = LRRS_FixIn%Atmos%OMEGAPHASFUNC_ELASTIC_UP
      OMEGAPHASFUNC_ELASTIC_DN = LRRS_FixIn%Atmos%OMEGAPHASFUNC_ELASTIC_DN

!  Modified Atmosphere Inputs

      GEOMETRY_SPECHEIGHT = LRRS_ModIn%MAtmos%GEOMETRY_SPECHEIGHT

!  Fixed Surface Inputs. Removed version 2.5, 9/11/15
!   This kind of input is now found in the LRRS BRDF supplement.
!      LAMBERTIAN_FRACTION = LRRS_FixIn%Surf%LAMBERTIAN_FRACTION
!      NSTREAMS_BRDF       = LRRS_FixIn%Surf%NSTREAMS_BRDF
!      BRDF_NAMES          = LRRS_FixIn%Surf%BRDF_NAMES
!      WHICH_BRDF          = LRRS_FixIn%Surf%WHICH_BRDF
!      BRDF_FACTOR         = LRRS_FixIn%Surf%BRDF_FACTOR
!      BRDF_NPARS          = LRRS_FixIn%Surf%BRDF_NPARS
!      BRDF_PARS           = LRRS_FixIn%Surf%BRDF_PARS

!  Fixed surface inputs

      ALBEDOS_RANKED      = LRRS_FixIn%Surf%ALBEDOS_RANKED

!  BRDF inputs from the supplement-generated structure. New, Version 2.5, 9/23/15

      DO_BRDF_Wav1 = LRRS_Sup%BRDF%DO_BRDF_Wav1

      EXACTDB_BRDFUNC(1:N_USER_STREAMS,1:N_USER_RELAZMS,:)        = &
                     LRRS_Sup%BRDF%EXACTDB_BRDFUNC(1:N_USER_STREAMS,1:N_USER_RELAZMS,:)

      BRDF_F_0 (0:MAX_MOMENTS,1:NSTREAMS,:)                       = & 
                     LRRS_Sup%BRDF%BRDF_F_0(0:MAX_MOMENTS,1:NSTREAMS,:)
      BRDF_F   (0:MAX_MOMENTS,1:NSTREAMS,1:NSTREAMS,:)            = & 
                     LRRS_Sup%BRDF%BRDF_F(0:MAX_MOMENTS,1:NSTREAMS,1:NSTREAMS,:)
      USER_BRDF_F_0 (0:MAX_MOMENTS,1:N_USER_STREAMS,:)            = &
                     LRRS_Sup%BRDF%USER_BRDF_F_0(0:MAX_MOMENTS,1:N_USER_STREAMS,:)
      USER_BRDF_F   (0:MAX_MOMENTS,1:N_USER_STREAMS,1:NSTREAMS,:) = &
                     LRRS_Sup%BRDF%USER_BRDF_F(0:MAX_MOMENTS,1:N_USER_STREAMS,1:NSTREAMS,:)

!  Fixed Surface-leaving inputs. Removed version, 2.5 9/11/15
!  Add SLTERM flag     . @@@ Rob Fix 06 Sep 12.
!  New Isotropic SLTERM. @@@ Rob Fix 06 sep 12
!      DO_LRRS_SLTERM = LRRS_FixIn%Bool%DO_LRRS_SLTERM
!      LRRS_SLTERM    = LRRS_FixIn%Spect%LRRS_SLTERM

!  Surface-leaving inputs from the supplement-generated structure. New, Version 2.5, 9/23/15

      DO_SLEAVE_Wav1 = LRRS_Sup%SLEAVE%DO_SLEAVE_Wav1

      SLTERM_ISOTROPIC(:)                                    = LRRS_Sup%SLEAVE%SLTERM_ISOTROPIC(:)
      SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,:) = &
                 LRRS_Sup%SLEAVE%SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,:)

      SLTERM_F_0(0:MAX_MOMENTS,1:NSTREAMS,:)            = &
                 LRRS_Sup%SLEAVE%SLTERM_F_0(0:MAX_MOMENTS,1:NSTREAMS,:)
      USER_SLTERM_F_0(0:MAX_MOMENTS,1:N_USER_STREAMS,:) = &
                 LRRS_Sup%SLEAVE%USER_SLTERM_F_0(0:MAX_MOMENTS,1:N_USER_STREAMS,:)

!  Linearized Inputs
!  -----------------

!  Boolean Inputs

      DO_PROFILE_WFS           = LRRS_LinIn%Cont%DO_PROFILE_WFS
      DO_COLUMN_WFS            = LRRS_LinIn%Cont%DO_COLUMN_WFS
      DO_SURFACE_WFS           = LRRS_LinIn%Cont%DO_SURFACE_WFS
      DO_SLEAVE_WFS            = LRRS_LinIn%Cont%DO_SLEAVE_WFS

      DO_NORMALIZED_WFS        = LRRS_LinIn%Cont%DO_NORMALIZED_WFS
      DO_AIRPROFILE_WFS        = LRRS_LinIn%Cont%DO_AIRPROFILE_WFS
      DO_TEMPPROFILE_WFS       = LRRS_LinIn%Cont%DO_TEMPPROFILE_WFS
      DO_TEMPSHIFT_WF          = LRRS_LinIn%Cont%DO_TEMPSHIFT_WF

!  Linearized Atmosphere Inputs

      N_TOTALCOLUMN_WFS        = LRRS_LinIn%Cont%N_TOTALCOLUMN_WFS

      LAYER_VARY_FLAG          = LRRS_LinIn%Cont%LAYER_VARY_FLAG
      LAYER_VARY_NUMBER        = LRRS_LinIn%Cont%LAYER_VARY_NUMBER

      LAYER_AIRCOLUMNS_DT      = LRRS_LinIn%Optical%LAYER_AIRCOLUMNS_DT
      TEMPERATURES_UNSHIFTED   = LRRS_LinIn%Optical%TEMPERATURES_UNSHIFTED

      L_DELTAU_INPUT_UNSCALED      = LRRS_LinIn%Optical%L_DELTAU_INPUT_UNSCALED
      L_OMEGAMOMS_ELASTIC_UNSCALED = LRRS_LinIn%Optical%L_OMEGAMOMS_ELASTIC_UNSCALED

!   -- Rob mod 5/12/17 for 2p5a, Setting local L_OMEGA_UNSCLAED

      L_OMEGA_UNSCALED(:,:,:)      = L_OMEGAMOMS_ELASTIC_UNSCALED(:,:,0,:)

!   -- Rob mod 5/12/17 for 2p5a, Copy Phase function product input

      L_OMEGAPHASFUNC_ELASTIC_UP = LRRS_LinIn%Optical%L_OMEGAPHASFUNC_ELASTIC_UP
      L_OMEGAPHASFUNC_ELASTIC_DN = LRRS_LinIn%Optical%L_OMEGAPHASFUNC_ELASTIC_DN
  
!  Linearized Surface Inputs.  Upgraded, Version 2.5

!      BRDFPAR_DERIV_INDEX = LRRS_LinIn%Surf%BRDFPAR_DERIV_INDEX
      N_SURFACE_WFS = LRRS_LinIn%Cont%N_SURFACE_WFS

!  BRDF

      LS_EXACTDB_BRDFUNC(1:MAX_SURFACEWFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,:)        = &
        LRRS_LinSup%BRDF%LS_EXACTDB_BRDFUNC(1:MAX_SURFACEWFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,:)

      LS_BRDF_F_0 (1:MAX_SURFACEWFS,0:MAX_MOMENTS,1:NSTREAMS,:)                       = & 
        LRRS_LinSup%BRDF%LS_BRDF_F_0(1:MAX_SURFACEWFS,0:MAX_MOMENTS,1:NSTREAMS,:)
      LS_BRDF_F   (1:MAX_SURFACEWFS,0:MAX_MOMENTS,1:NSTREAMS,1:NSTREAMS,:)            = & 
        LRRS_LinSup%BRDF%LS_BRDF_F(1:MAX_SURFACEWFS,0:MAX_MOMENTS,1:NSTREAMS,1:NSTREAMS,:)
      LS_USER_BRDF_F_0 (1:MAX_SURFACEWFS,0:MAX_MOMENTS,1:N_USER_STREAMS,:)            = &
        LRRS_LinSup%BRDF%LS_USER_BRDF_F_0(1:MAX_SURFACEWFS,0:MAX_MOMENTS,1:N_USER_STREAMS,:)
      LS_USER_BRDF_F   (1:MAX_SURFACEWFS,0:MAX_MOMENTS,1:N_USER_STREAMS,1:NSTREAMS,:) = &
        LRRS_LinSup%BRDF%LS_USER_BRDF_F(1:MAX_SURFACEWFS,0:MAX_MOMENTS,1:N_USER_STREAMS,1:NSTREAMS,:)

!  SLEAVE

      N_SLEAVE_WFS  = LRRS_LinIn%Cont%N_SLEAVE_WFS

      LSSL_SLTERM_ISOTROPIC(1:N_SLEAVE_WFS,:)                                    = &
        LRRS_LinSup%SLEAVE%LSSL_SLTERM_ISOTROPIC(1:N_SLEAVE_WFS,:)
      LSSL_SLTERM_USERANGLES(1:N_SLEAVE_WFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,:) = &
        LRRS_LinSup%SLEAVE%LSSL_SLTERM_USERANGLES(1:N_SLEAVE_WFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,:)
      LSSL_SLTERM_F_0(1:N_SLEAVE_WFS,0:MAX_MOMENTS,1:NSTREAMS,:)                 = &
        LRRS_LinSup%SLEAVE%LSSL_SLTERM_F_0(1:N_SLEAVE_WFS,0:MAX_MOMENTS,1:NSTREAMS,:)
      LSSL_USER_SLTERM_F_0(1:N_SLEAVE_WFS,0:MAX_MOMENTS,1:N_USER_STREAMS,:)      = &
        LRRS_LinSup%SLEAVE%LSSL_USER_SLTERM_F_0(1:N_SLEAVE_WFS,0:MAX_MOMENTS,1:N_USER_STREAMS,:)

      N_TOTALSURFACE_WFS = N_SURFACE_WFS + N_SLEAVE_WFS

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Define local values of Npoints

      IF ( DO_BIN_REALIZATION ) THEN
        NPOINTS_TOTAL = NPOINTS_OUTER
        NPOINTS_LOCAL = NPOINTS_INNER
      ELSE
        NPOINTS_TOTAL = NPOINTS_MONO
        NPOINTS_LOCAL = 1
      ENDIF

!  LRRS input debug

      IF (DO_DEBUG_INPUT) THEN
        CALL LRRS_DEBUG_INPUT_MASTER()
        IF (DO_COLUMN_WFS .OR. DO_PROFILE_WFS .OR. DO_SURFACE_WFS) &
          CALL LRRS_DEBUG_LIN_INPUT_MASTER()
      END IF

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
      FOURIER_SAVED   = 0

!  Status Outputs

      STATUS     = LRRS_SUCCESS
      N_MESSAGES = 0
      MESSAGES   = ' '

!Rob add 10/10/18 - initialize OpenMP timing
      IF ( DO_TIMING ) THEN
        OMP_Elastic_Time = 0.0
        OMP_Raman_Time   = 0.0
      ENDIF

!  Linearized Atmosphere Outputs

      !Elastic Field
      LC_ELASTIC_UP     = ZERO
      LC_ELASTIC_DN     = ZERO
      LP_ELASTIC_UP     = ZERO
      LP_ELASTIC_DN     = ZERO

      LC_ELASTIC_SS_UP  = ZERO
      LC_ELASTIC_SS_DN  = ZERO
      LP_ELASTIC_SS_UP  = ZERO
      LP_ELASTIC_SS_DN  = ZERO

      L_MEAN_ELASTIC_UP = ZERO
      L_MEAN_ELASTIC_DN = ZERO
      L_FLUX_ELASTIC_UP = ZERO
      L_FLUX_ELASTIC_DN = ZERO

      !Raman Field
      LC_RAMAN_UP       = ZERO
      LC_RAMAN_DN       = ZERO
      LP_RAMAN_UP       = ZERO
      LP_RAMAN_DN       = ZERO

      LC_RAMAN_SS_UP    = ZERO
      LC_RAMAN_SS_DN    = ZERO
      LP_RAMAN_SS_UP    = ZERO
      LP_RAMAN_SS_DN    = ZERO

      L_MEAN_RAMAN_UP   = ZERO
      L_MEAN_RAMAN_DN   = ZERO
      L_FLUX_RAMAN_UP   = ZERO
      L_FLUX_RAMAN_DN   = ZERO

!  Linearized Surface Outputs

      !Elastic Field
      LS_ELASTIC_UP      = ZERO
      LS_ELASTIC_DN      = ZERO

      LS_ELASTIC_SS_UP   = ZERO

      LS_MEAN_ELASTIC_UP = ZERO
      LS_MEAN_ELASTIC_DN = ZERO
      LS_FLUX_ELASTIC_UP = ZERO
      LS_FLUX_ELASTIC_DN = ZERO

      !Raman Field
      LS_RAMAN_UP        = ZERO
      LS_RAMAN_DN        = ZERO

      LS_RAMAN_SS_UP     = ZERO

      LS_MEAN_RAMAN_UP   = ZERO
      LS_MEAN_RAMAN_DN   = ZERO
      LS_FLUX_RAMAN_UP   = ZERO
      LS_FLUX_RAMAN_DN   = ZERO

!  Internal settings
!  =================

!mick fix 7/20/2016 - commented out the defining (and possible re-defining) of
!                     the three following local control flags w/o user knowledge
!                   - added consistency checks for these flags in LRRS_CHECK_INPUTS below

!  Elastic-only flag

!      DO_ELASTIC_ONLY = .NOT.DO_RRS_OVERALL .AND. .NOT.DO_DIRECTRRS_ONLY

!  Exclusion flags

!      DO_CABANNES_RAMAN = .NOT.DO_ENERGY_BALANCING
!      IF ( .NOT.DO_ELASTIC_ONLY ) THEN
!        DO_MONO_REALIZATION = .NOT.DO_BIN_REALIZATION
!      ENDIF

!  Check inputs. Exit if any failure
!  =================================

!  Version 2.5. removal of specific BRDF information

      CALL LRRS_CHECK_INPUTS &
       ( LRRS_FixIn, LRRS_ModIn, N_USER_STREAMS, N_USER_RELAZMS, N_LOUTPUT, & ! Inputs
         SOLAR_ANGLE, USER_ANGLES, USER_RELAZMS, NLAYERS, HEIGHT_GRID,      & ! Inputs
         LBOUNDARIES_OUTPUT, LPARTIALS_OUTPUT, EARTH_RADIUS,                & ! Inputs
         STATUS, N_MESSAGES, MESSAGES )

      IF ( STATUS .NE. LRRS_SUCCESS ) THEN
        N_MESSAGES = N_MESSAGES + 1
        MESSAGES(N_MESSAGES) = 'Failure from Input checking, L_main'
        LRRS_Out%Status%STATUS_CALCULATION = STATUS
        LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
        LRRS_Out%Status%MESSAGES           = MESSAGES
        RETURN
      ENDIF

!  Additional check. ******** New 4.19.2010 **********
!    Height linearizations. DISABLED 10 MARCH 2011.
!      IF ( DO_HLEVEL_LINEARIZATION ) THEN
!        IF ( .not. DO_COLUMN_WFS ) THEN
!          n_messages = n_messages + 1
!          messages(n_messages) =
!     &        'Height linearization requires column WF flag'
!          STATUS = LRRS_SERIOUS
!          RETURN
!        ENDIF
!        IF ( NFINELAYERS .NE. 0 ) THEN
!          n_messages = n_messages + 1
!          messages(n_messages) =
!     &           'H-linearization: Reset NLINELAYERS = 0'
!          STATUS = LRRS_WARNING
!          NFINELAYERS = 0
!        ENDIF
!      ENDIF

!  if there's no azimuth dependence, just do one value in azimuth loop

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_N_USERAZM = 1
        N_USER_RELAZMS  = 1
       ELSE
        LOCAL_N_USERAZM = N_USER_RELAZMS
      ENDIF

!  Set Number of Geometries, and the offsetting

      N_GEOMETRIES = LOCAL_N_USERAZM * N_USER_STREAMS
      DO UM = 1, N_USER_STREAMS
        GEOM_OFFSETS(UM) = LOCAL_N_USERAZM * ( UM - 1)
      ENDDO

!  Geometry adjustment. New section 18 March 2011
!  ==============================================

!  Adjust surface condition

      ADJUST_SURFACE = .FALSE.
      IF ( DO_SSCORR_OUTGOING ) THEN
        IF (HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) THEN
         ADJUST_SURFACE = .TRUE.
        ENDIF
      ENDIF

!  Perform adjustment

      MODIFIED_ERADIUS = EARTH_RADIUS + GEOMETRY_SPECHEIGHT
      CALL MULTI_OUTGOING_ADJUSTGEOM &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS,                     & ! Inputs
          N_USER_STREAMS,   N_USER_RELAZMS,                       & ! Inputs
          HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE, & ! Inputs
          USER_ANGLES, SOLAR_ANGLE, USER_RELAZMS,                 & ! Inputs
          USER_ANGLES_ADJUST, SOLAR_ANGLES_ADJUST, USER_RELAZMS_ADJUST, & ! Outputs
          FAIL, MESSAGE_SUB, MESSAGE )                                    ! Outputs

!   -- Rob mod 5/12/17 for 2p5a, Need to generate local Cosscat variables for the Setups
!       ** COSSCAT variables here will be same as those used internally in SSCORR routines
!       ** Same as those in External supplement for PHASFUNC_ELASTIC products  ??????????????????

      if ( DO_SSCORR_OUTGOING ) then
        DO UM = 1, N_USER_STREAMS
          ALPHA_BOA = USER_ANGLES_ADJUST(UM)
          DO UA = 1, N_USER_RELAZMS
            THETA_BOA = SOLAR_ANGLES_ADJUST(UM,UA)
            PHI_BOA   = USER_RELAZMS_ADJUST(UM,UA)
            V = N_USER_RELAZMS * (UM-1) + UA
            stheta_boa = sin(theta_boa * deg_to_rad) ; ctheta_boa = sqrt(one-stheta_boa*stheta_boa)
            salpha_boa = sin(alpha_boa * deg_to_rad) ; calpha_boa = sqrt(one-salpha_boa*salpha_boa)
            cphi_boa   = cos(phi_boa * deg_to_rad)
            cosscat_up (v) = - calpha_boa * ctheta_boa + salpha_boa * stheta_boa * cphi_boa
            cosscat_dn (v) = + calpha_boa * ctheta_boa + salpha_boa * stheta_boa * cphi_boa
          ENDDO
        ENDDO
      ENDIF

!   -- Rob mod 5/12/17 for 2p5a, Need to generate local Cosscat variables for the Setups
!       ** COSSCAT variables here will be same as those used internally in SSCORR routines
!       ** Same as those in External supplement for PHASFUNC_ELASTIC products  ??????????????????

      if ( DO_SSCORR_NADIR ) then
        DO UM = 1, N_USER_STREAMS
         ALPHA_BOA = USER_ANGLES(UM)
         DO UA = 1, N_USER_RELAZMS
            THETA_BOA = SOLAR_ANGLE
            PHI_BOA   = USER_RELAZMS(UA)
            V = N_USER_RELAZMS * (UM-1) + UA
            stheta_boa = sin(theta_boa * deg_to_rad) ; ctheta_boa = sqrt(one-stheta_boa*stheta_boa)
            salpha_boa = sin(alpha_boa * deg_to_rad) ; calpha_boa = sqrt(one-salpha_boa*salpha_boa)
            cphi_boa   = cos(phi_boa * deg_to_rad)
            cosscat_up (v) = - calpha_boa * ctheta_boa + salpha_boa * stheta_boa * cphi_boa
            cosscat_dn (v) = + calpha_boa * ctheta_boa + salpha_boa * stheta_boa * cphi_boa
          ENDDO
        ENDDO
      ENDIF

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
        ( DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                   & ! Inputs
          DO_BIN_REALIZATION, DO_MOLECSCAT_ONLY, DO_LBOUNDARIES,         & ! Inputs
          NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS, N_LOUTPUT,           & ! Inputs
          NLAYERS, NPOINTS_INNER, SOLAR_ANGLE, USER_ANGLES,              & ! Inputs
          NMOMENTS, N_OUT_STREAMS, LBOUNDARIES_OUTPUT, LPARTIALS_OUTPUT, & ! Inputs
          FLUX_FACTOR, COS_SZA, LEVELMASK_UP, LEVELMASK_DN,              & ! Outputs
          STERM_MASK_UP, STERM_MASK_DN, N_CONVTESTS,                     & ! Outputs
          QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS )         ! Outputs

!  Bookkeeping code for linearization
!  ==================================

!mick chg 7/20/2016 - set linearization to 0 (leave for now)

      nkstorage  = 0
      nkstorage2 = 0

      if ( do_column_wfs ) then
        do n = 1, nlayers
          nkstorage(n,0)  = n
          nkstorage2(n,0) = n
        enddo
      else if ( do_profile_wfs ) then
        do n = 1, nlayers
          nk = n*(n-1)/2
          do k = 1, n
            nkstorage(n,k) = nk + k
          enddo
          n2 = nlayers*(n-1)
          do k = 1, nlayers
            nkstorage2(n,k) = n2 + k
          enddo
        enddo
      endif

!  Overall linearization flag for atmosphere

      do_atmos_wfs = ( do_profile_wfs .or. do_column_wfs )

!  Grid level flag  . New May 2010. Omitted for now..........
!      DO_GRIDLEVEL_WFS = DO_HLEVEL_LINEARIZATION

!  Setups for Raman problem
!  ========================

!    - delta-m scaling with option for linearization
!    - Raman spectroscopy
!    - Raman optical depths with option for linearization

!   -- Rob mod 5/12/17 for 2p5a, remove N_MOMENTS_INPUT, add Geometry information (flags, cosines)
!   -- Rob mod 5/12/17 for 2p5a, add PHASE FUNCTION Product input (Elastic) and output (Cabannes)
!    - New derivatives option, 7 May 2010. DISABLED.

      IF ( DO_BIN_REALIZATION ) THEN

        CALL L_RRSOPTICAL_MASTER_BIN &
        ( DO_DELTAM_SCALING, DO_ELASTIC_ONLY, DO_UPWELLING,          & ! Inputs
          DO_DNWELLING, DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,      & ! Inputs
          DO_PROFILE_WFS, DO_COLUMN_WFS, DO_AIRPROFILE_WFS,          & ! Inputs
          DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF, DO_NORMALIZED_WFS,    & ! Inputs
          NLAYERS, NSTREAMS, N_GEOMETRIES,                           & ! Inputs
          NPOINTS_INNER, OFFSET_INNER, NPOINTS_OUTER,                & ! Inputs
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                        & ! Inputs
          LAMBDAS_RANKED, FLUXES_RANKED, BINLOWER, BINUPPER,         & ! Inputs
          COSSCAT_UP, COSSCAT_DN, LAYER_TEMPERATURES,                & ! Inputs
          LAYER_AIRCOLUMNS, RAYLEIGH_XSEC, RAYLEIGH_DEPOL,           & ! Inputs
          TEMPERATURES_UNSHIFTED, LAYER_AIRCOLUMNS_dT,               & ! Inputs
          DELTAU_INPUT_UNSCALED,      L_DELTAU_INPUT_UNSCALED,       & ! Inputs
          OMEGAMOMS_ELASTIC_UNSCALED, L_OMEGAMOMS_ELASTIC_UNSCALED,  & ! Inputs
          OMEGAPHASFUNC_ELASTIC_UP,   L_OMEGAPHASFUNC_ELASTIC_UP,    & ! Inputs
          OMEGAPHASFUNC_ELASTIC_DN,   L_OMEGAPHASFUNC_ELASTIC_DN,    & ! Inputs
          N_RRSBINS, BINMAP, RINGSPEC_1,                             & ! Outputs
          NMOMENTS, TRUNC_FACTORS,    L_TRUNC_FACTORS,               & ! Outputs
          DELTAU_VERT_INPUT,          L_DELTAU_VERT_INPUT,           & ! Outputs
          OMEGAMOMS_ELASTIC,          L_OMEGAMOMS_ELASTIC,           & ! Outputs
          OMEGAMOMS_CABANNES,         L_OMEGAMOMS_CABANNES,          & ! Outputs
          OMEGAPHASFUNC_CABANNES_UP,  L_OMEGAPHASFUNC_CABANNES_UP,   & ! Inputs
          OMEGAPHASFUNC_CABANNES_DN,  L_OMEGAPHASFUNC_CABANNES_DN,   & ! Inputs
          OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSLOSS_UNSCALED,  & ! Outputs
          OMEGAMOMS_RRSBIN_UNSCALED,  L_OMEGAMOMS_RRSBIN_UNSCALED,   & ! Outputs
          OMEGAMOMS_RRSLOSS,          L_OMEGAMOMS_RRSLOSS,           & ! Outputs
          OMEGAMOMS_RRSBIN,           L_OMEGAMOMS_RRSBIN,            & ! Outputs
          FAIL, MESSAGE, MESSAGE_SUB )                                 ! Outputs

      ELSE

        CALL L_RRSOPTICAL_MASTER_MONO &
        ( DO_DELTAM_SCALING, DO_ELASTIC_ONLY, DO_UPWELLING,          & ! Inputs
          DO_DNWELLING, DO_ENERGY_BALANCING, DO_CABANNES_RAMAN,      & ! Inputs
          DO_PROFILE_WFS, DO_COLUMN_WFS, DO_AIRPROFILE_WFS,          & ! Inputs
          DO_TEMPPROFILE_WFS, DO_TEMPSHIFT_WF, DO_NORMALIZED_WFS,    & ! Inputs
          NLAYERS, NSTREAMS, N_GEOMETRIES,                           & ! Inputs
          LAMBDA_EXCIT, NPOINTS_MONO, W_EXCIT,                       & ! Inputs
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                        & ! Inputs
          COSSCAT_UP, COSSCAT_DN, LAYER_TEMPERATURES,                & ! Inputs
          TEMPERATURES_UNSHIFTED, LAYER_AIRCOLUMNS,                  & ! Inputs
          LAYER_AIRCOLUMNS_dT, RAYLEIGH_XSEC, RAYLEIGH_DEPOL,        & ! Inputs
          DELTAU_INPUT_UNSCALED,      L_DELTAU_INPUT_UNSCALED,       & ! Inputs
          OMEGAMOMS_ELASTIC_UNSCALED, L_OMEGAMOMS_ELASTIC_UNSCALED,  & ! Inputs
          OMEGAPHASFUNC_ELASTIC_UP,   L_OMEGAPHASFUNC_ELASTIC_UP,    & ! Inputs
          OMEGAPHASFUNC_ELASTIC_DN,   L_OMEGAPHASFUNC_ELASTIC_DN,    & ! Inputs
          NMOMENTS, TRUNC_FACTORS,    L_TRUNC_FACTORS,               & ! Outputs
          DELTAU_VERT_INPUT,          L_DELTAU_VERT_INPUT,           & ! Outputs
          OMEGAMOMS_ELASTIC,          L_OMEGAMOMS_ELASTIC,           & ! Outputs
          OMEGAMOMS_CABANNES,         L_OMEGAMOMS_CABANNES,          & ! Outputs
          OMEGAPHASFUNC_CABANNES_UP,  L_OMEGAPHASFUNC_CABANNES_UP,   & ! Inputs
          OMEGAPHASFUNC_CABANNES_DN,  L_OMEGAPHASFUNC_CABANNES_DN,   & ! Inputs
          OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSLOSS_UNSCALED,  & ! Outputs
          OMEGAMOMS_RRSGAIN_UNSCALED, L_OMEGAMOMS_RRSGAIN_UNSCALED,  & ! Outputs
          OMEGAMOMS_RRSLOSS,          L_OMEGAMOMS_RRSLOSS,           & ! Outputs
          OMEGAMOMS_RRSGAIN,          L_OMEGAMOMS_RRSGAIN )            ! Outputs

      ENDIF

!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@
!  Exception handling for MAX_BINS Dimension
      IF ( DO_BIN_REALIZATION ) THEN
         IF  ( FAIL ) THEN
            STATUS = LRRS_SERIOUS              ; N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = MESSAGE     ; N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = MESSAGE_SUB ; N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = 'Failure from RRSOPTICAL_MASTER_BIN'
            LRRS_Out%Status%STATUS_CALCULATION = STATUS
            LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
            LRRS_Out%Status%MESSAGES           = MESSAGES
            RETURN
         ENDIF
      ENDIF
!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@

!  Chapman function calculation of CHAPMAN_FACTORS

      CALL CHAPMAN_FUNCTION &
           ( DO_PLANE_PARALLEL, NLAYERS,         & ! Inputs
             COS_SZA, EARTH_RADIUS, HEIGHT_GRID, & ! Inputs
             CHAPMAN )                             ! Outputs

!  Disabled code (alternative)
!      IF ( DO_COLUMN_WFS .and. DO_HLEVEL_LINEARIZATION ) THEN
!        CALL CHAPMAN_FUNCTION_PLUS
!     I     ( DO_PLANE_PARALLEL, NLAYERS, N_TOTALCOLUMN_WFS,
!     I       COS_SZA, EARTH_RADIUS, HEIGHT_GRID, LH_HEIGHT_GRID,
!     O       CHAPMAN, L_CHAPMAN )
!      ELSE
!        CALL CHAPMAN_FUNCTION
!     I     ( DO_PLANE_PARALLEL, NLAYERS,
!     I       COS_SZA, EARTH_RADIUS, HEIGHT_GRID,
!     O       CHAPMAN )
!      ENDIF

!  Get the Beam average secants, beam transmittances and view-angle transmittances
!   -- Rob mod 5/12/17 for 2p5a, This routine moved, was formerly just before FOURIER loop
!   -- Rob mod 5/12/17 for 2p5a, Added DO_SSCORR_NADIR flag,renamed SSFULL

!      IF (.NOT. DO_SSFULL ) THEN
      IF ( DO_SSCORR_NADIR .or. .NOT.DO_SSCORR_ALONE ) THEN
        CALL ATTENUATION_SETUP_1 &
          ( NPOINTS_TOTAL, NLAYERS, DO_PLANE_PARALLEL,     & ! Input
            DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS, & ! Input
            DELTAU_VERT_INPUT, CHAPMAN, COS_SZA,           & ! Input
            BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT,     & ! Output
            BEAM_ETRANS, BEAM_DTRANS, SAVE_TRANS_USERM )     ! Output
      ENDIF

!  Profile WFS
!  -----------

!   -- Rob mod 5/12/17 for 2p5a, This routine moved, was formerly just before FOURIER loop
!   -- Rob mod 5/12/17 for 2p5a, Added DO_SSCORR_NADIR flag, renamed SSFULL

!      IF (.NOT. DO_SSFULL ) THEN
      IF ( DO_SSCORR_NADIR .or. .NOT.DO_SSCORR_ALONE ) THEN
        IF ( DO_PROFILE_WFS ) THEN
          CALL LP_ATTENUATION_SETUP_1 &
          ( DO_PLANE_PARALLEL, NPOINTS_TOTAL, NLAYERS, NKSTORAGE, & ! Inputs
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                   & ! Inputs
            DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS,        & ! Inputs
            DELTAU_VERT_INPUT, CHAPMAN, BEAM_PICUTOFF,            & ! Inputs
            BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS,              & ! Inputs
            SAVE_TRANS_USERM, L_DELTAU_VERT_INPUT,                & ! Inputs
            L_BEAM_ITRANS, L_BEAM_AVSECANT,                       & ! Outputs
            L_BEAM_ETRANS, L_BEAM_DTRANS, L_SAVE_TRANS_USERM )      ! Outputs
        ENDIF
      ENDIF

!  Column WFS
!  ----------

!   -- Rob mod 5/12/17 for 2p5a, This routine moved, was formerly just before FOURIER loop
!   -- Rob mod 5/12/17 for 2p5a, Added DO_SSCORR_NADIR flag, renamed SSFULL

!      IF (.NOT. DO_SSFULL ) THEN
      IF ( DO_SSCORR_NADIR .or. .NOT.DO_SSCORR_ALONE ) THEN
        IF ( DO_COLUMN_WFS ) THEN
          CALL LC_ATTENUATION_SETUP_1 &
          ( DO_PLANE_PARALLEL, NPOINTS_TOTAL, NLAYERS, N_TOTALCOLUMN_WFS, & ! Inputs
            DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS,                & ! Inputs
            DELTAU_VERT_INPUT, CHAPMAN, BEAM_PICUTOFF,                    & ! Inputs
            BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS,                      & ! Inputs
            SAVE_TRANS_USERM, L_DELTAU_VERT_INPUT,                        & ! Inputs
            L_BEAM_ITRANS, L_BEAM_AVSECANT,                               & ! Outputs 
            L_BEAM_ETRANS, L_BEAM_DTRANS, L_SAVE_TRANS_USERM )              ! Outputs
        ENDIF
      ENDIF

!  Get the Linearized Beam Multipliers
!   -- Rob mod 5/12/17 for 2p5a, New routine for the DO_SSCORR_NADIR option

      IF ( DO_SSCORR_NADIR .and. DO_USER_STREAMS ) THEN
        IF ( DO_ATMOS_WFS ) THEN
          CALL LPC_BEAM_MULTIPLIERS_1 &
           ( DO_UPWELLING, DO_DNWELLING, DO_COLUMN_WFS, DO_PROFILE_WFS,         & ! Inputs
             NPOINTS_TOTAL, NLAYERS, N_USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER,& ! Inputs
             LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS, NKSTORAGE2, & ! Inputs
             STERM_MASK_UP, STERM_MASK_DN, USER_STREAMS,                        & ! Inputs
             DELTAU_VERT_INPUT,   SAVE_TRANS_USERM,                             & ! Inputs
             L_DELTAU_VERT_INPUT, L_SAVE_TRANS_USERM,                           & ! Inputs
             BEAM_PICUTOFF, BEAM_AVSECANT, BEAM_DTRANS, BEAM_ITRANS,            & ! Inputs
             L_BEAM_ITRANS, L_BEAM_AVSECANT, L_BEAM_DTRANS,                     & ! Inputs
             SAVE_BEAMMULT_UP, L_SAVE_BEAMMULT_UP,                              & ! Output
             SAVE_BEAMMULT_DN, L_SAVE_BEAMMULT_DN )                               ! Output
        ELSE
          CALL BEAM_MULTIPLIERS_1 &
           ( DO_UPWELLING, DO_DNWELLING, NPOINTS_TOTAL, NLAYERS, & ! Inputs
             N_USER_STREAMS, TAYLOR_SMALL, TAYLOR_ORDER,         & ! Inputs
             STERM_MASK_UP, STERM_MASK_DN,                       & ! Inputs
             USER_STREAMS, DELTAU_VERT_INPUT, SAVE_TRANS_USERM,  & ! Inputs
             BEAM_PICUTOFF, BEAM_AVSECANT, BEAM_DTRANS, BEAM_ITRANS, & ! Inputs
             SAVE_BEAMMULT_UP, SAVE_BEAMMULT_DN )                  ! Output
        ENDIF
      ENDIF

!  Disable LH linearization. Old code suggested
!      IF (.NOT. DO_SSFULL ) THEN
!        IF ( DO_COLUMN_WFS ) THEN
!          IF ( DO_HLEVEL_LINEARIZATION ) THEN
!            CALL LCH_ATTENUATION_SETUP_1 &
!               ( DO_PLANE_PARALLEL, NPOINTS_TOTAL,NLAYERS, N_TOTALCOLUMN_WFS, &
!                 DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS, &
!                 DELTAU_VERT_INPUT, CHAPMAN, L_CHAPMAN, BEAM_PICUTOFF, &
!                 BEAM_ITRANS, BEAM_AVSECANT, BEAM_DTRANS, &
!                 SAVE_TRANS_USERM, L_DELTAU_VERT_INPUT,&
!                 L_BEAM_ITRANS, L_BEAM_AVSECANT,        &
!                 L_BEAM_ETRANS, L_BEAM_DTRANS, L_SAVE_TRANS_USERM )
!         ENDIF
!       ENDIF

!  BRDF input functions (non-Lambertian)
!   CODE REMOVED, VERSION 2.5. Replaced by Supplement-derived input

!  Write input variables (if flagged). Disabled Version 2.5.
!   ........... linearized routine not yet available

      IF ( DO_LRRS_WRITEINPUT ) THEN
        IUNIT = LRRS_INUNIT
        OPEN(IUNIT,FILE=DEBUG_FILENAMES(1),STATUS='UNKNOWN')
!        CALL LRRS_WRITEINPUT ( IUNIT )
!        IF ( DO_ATMOS_WFS ) CALL L_LRRS_WRITEINPUT ( IUNIT )
        CLOSE(IUNIT)
      ENDIF

!  Write scenario variables (if flagged). Disabled Version 2.5.
!   ........... linearized routine not yet available

      IF ( DO_LRRS_WRITESCENARIO ) THEN
        SUNIT = LRRS_SCENUNIT
        OPEN(SUNIT,FILE=DEBUG_FILENAMES(2),STATUS='UNKNOWN')
!        CALL LRRS_WRITESCENARIO ( SUNIT )
!        IF ( DO_ATMOS_WFS ) CALL L_LRRS_WRITESCENARIO ( SUNIT )
        CLOSE(SUNIT)
      ENDIF

!  Get the Single scatter elastic correction. Now includes the DB term
!    DB correction is automatic now, even for the Lambertian surface
!      28 May 2008,  RT Solutions Inc.

!   -- Rob mod 5/12/17 for 2p5a, Switch to Phase-function inputs
!   -- Rob mod 5/12/17 for 2p5a. Use the SSCORR_GENERAL FLAG, extra inputs for NADIR SSCORR

      IF ( DO_SSCORR_GENERAL ) THEN

!  Bin realization, 2 possibilities

!         1 = Normal with no height linearization
!         2 = Height/Column linearization. No fine-layering. DISABLED

        IF ( DO_BIN_REALIZATION ) THEN

!  THIS IS POSSIBILITY 1
!   Note 17 October 2015. Master has not been compiled with this routine....

!   -- Rob mod 5/12/17 for 2p5a, Switch to Phase-function inputs

          CALL LRRS_SSCORR_GENERAL_BIN_PLUS_1 &
         ( DO_RRS_OVERALL, DO_UPWELLING, DO_DNWELLING,                    & ! Input Mode flags
           DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,                           & ! Input Mode flags
           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_DELTAM_SCALING,     & ! Input Mode flags
           DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,          & ! Input Surface flags
           DO_BRDF_Wav1, DO_Sleave_Wav1,                                  & ! Input Supplement flags
           DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,  & ! Inputs Lin.Control
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS,         & ! Inputs Lin.Control
           NKSTORAGE2, N_SURFACE_WFS, N_SLEAVE_WFS,                       & ! Inputs Lin.Control
           NLAYERS, NFINELAYERS, N_USER_STREAMS, N_USER_RELAZMS,          & ! Input Control Integers
           N_LOUTPUT, LEVELMASK_UP, LEVELMASK_DN,                         & ! Input levels for output
           NPOINTS_INNER, OFFSET_INNER, NPOINTS_OUTER, N_RRSBINS, BINMAP, & ! Input Raman bin control
           SOLAR_ANGLE, COSSCAT_UP, COSSCAT_DN,                           & ! Input Geometry
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
           L_OMEGAPHASFUNC_ELASTIC_UP,   L_OMEGAPHASFUNC_ELASTIC_DN,          & ! Input Lin Optical Elastic
           L_OMEGAPHASFUNC_CABANNES_UP,  L_OMEGAPHASFUNC_CABANNES_DN,         & ! Input Lin Optical Cabannes
           L_OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSBIN_UNSCALED,         & ! Input Lin Optical Raman
           L_SAVE_TRANS_USERM, L_SAVE_BEAMMULT_UP, L_SAVE_BEAMMULT_DN,    & ! Input Optical for NADIR SSCORR 2p5a
           ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN,        & ! Outputs
           LC_ELASTIC_SS_UP, LP_ELASTIC_SS_UP, LS_ELASTIC_SS_UP,          & ! Outputs
           LC_ELASTIC_SS_DN, LP_ELASTIC_SS_DN,                            & ! Outputs
           LC_RAMAN_SS_UP,   LP_RAMAN_SS_UP,   LS_RAMAN_SS_UP,            & ! Outputs
           LC_RAMAN_SS_DN,   LP_RAMAN_SS_DN,                              & ! Outputs
           FAIL, MESSAGE_SUB)                                               ! Outputs

          IF ( FAIL ) THEN
            STATUS = LRRS_SERIOUS              ; N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = MESSAGE_SUB ; N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = 'Failure, LRRS_SSCORR_GENERAL_BIN_PLUS_1'
            LRRS_Out%Status%STATUS_CALCULATION = STATUS
            LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
            LRRS_Out%Status%MESSAGES           = MESSAGES
            RETURN
          ENDIF

!  THIS IS POSSIBILITY 2 (DISABLED, 10 March 2011). F77 CODE
!     REMARK, NEEDS THE ADJUSTED ANGLES (18 March 2011)
!          if (DO_COLUMN_WFS.and.DO_HLEVEL_LINEARIZATION) THEN
!            CALL LRRS_SSCORR_OGNFLH_BIN_PLUS_1
!     I   ( DO_RRS_OVERALL, DO_UPWELLING, DO_DNWELLING,
!     I     DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_WATERLEAVING,
!     I     DO_LAMBERTIAN_SURFACE, DO_LAMBERTIAN_FRACTION, DO_DELTAM_SCALING,
!     I     DO_COLUMN_WFS, DO_SURFACE_WFS, LHMASK_COLUMN_WFS,
!     I     NLAYERS, NMOMENTS_INPUT, N_TOTALCOLUMN_WFS,
!     I     N_USER_STREAMS, N_USER_RELAZMS, N_LOUTPUT,
!     I     NPOINTS_INNER, OFFSET_INNER, NPOINTS_OUTER, SOLAR_ANGLE,
!     I     USER_ANGLES, USER_RELAZMS, LEVELMASK_UP, LEVELMASK_DN,
!     I     LAMBERTIAN_FRACTION, EXACTDB_BRDFUNC, LS_EXACTDB_BRDFUNC,
!     I     FLUXES_RANKED, ALBEDOS_RANKED,
!     I      EARTH_RADIUS, HEIGHT_GRID, LH_HEIGHT_GRID,
!     I       DELTAU_VERT_INPUT,   OMEGAMOMS_ELASTIC_UNSCALED,
!     I     L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC_UNSCALED,
!     I     TRUNC_FACTORS, L_TRUNC_FACTORS, N_RRSBINS, BINMAP,
!     I     OMEGAMOMS_CABANNES_UNSCALED,L_OMEGAMOMS_CABANNES_UNSCALED,
!     I     OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSLOSS_UNSCALED,
!     I     OMEGAMOMS_RRSBIN_UNSCALED,  L_OMEGAMOMS_RRSBIN_UNSCALED,
!     O     ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN,
!     O     LC_ELASTIC_SS_UP, LC_ELASTIC_SS_DN, LS_ELASTIC_SS_UP,
!     O     LC_RAMAN_SS_UP,   LC_RAMAN_SS_DN,   LS_RAMAN_SS_UP,
!     O     FAIL, MESSAGE_SUB)
!            IF ( FAIL ) THEN
!              MESSAGE = 'Failure from LRRS_SSCORR_OGNFLH_BIN_PLUS_1'
!              RETURN
!            ENDIF
!         endif

!  Monochromatic realization, 2 possibilities

!         1 = Normal with no height linearization
!         2 = Height/Column linearization. No fine-layering. PLACEHOLDER

        ELSE IF ( DO_MONO_REALIZATION ) THEN

!  THIS IS POSSIBILITY 1
!   Note 17 October 2015. Master has not been compiled with this routine....

          CALL LRRS_SSCORR_GENERAL_MONO_PLUS_1 &
          ( DO_RRS_OVERALL, DO_UPWELLING, DO_DNWELLING,                   & ! Inputs
            DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,                          & ! Input Mode flags
            DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_DELTAM_SCALING,    & ! Inputs
            DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,         & ! Inputs
            DO_BRDF_Wav1, DO_Sleave_Wav1,                                 & ! Inputs
            DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,      & ! Inputs
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS,             & ! Inputs
            NKSTORAGE2, N_SURFACE_WFS, N_SLEAVE_WFS,                           & ! Inputs Lin.Control
            NLAYERS, NFINELAYERS, N_USER_STREAMS, N_USER_RELAZMS,              & ! Input Control Integers
            N_LOUTPUT, LEVELMASK_UP, LEVELMASK_DN,                             & ! Input levels for output
            NPOINTS_MONO, W_EXCIT,                                             & ! Inputs Raman Spec.
            SOLAR_ANGLE, COSSCAT_UP, COSSCAT_DN,                               & ! Input Geometry
            SOLAR_ANGLES_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST,      & ! Input Geometry
            FLUXES_RANKED, ALBEDOS_RANKED, EARTH_RADIUS, HEIGHT_GRID,          & ! Input Fluxes
            EXACTDB_BRDFUNC, SLTERM_ISOTROPIC, SLTERM_USERANGLES,              & ! Input   Surface
            TRUNC_FACTORS, DELTAU_VERT_INPUT, OMEGA_UNSCALED,                  & ! Input Optical
            OMEGAPHASFUNC_ELASTIC_UP,   OMEGAPHASFUNC_ELASTIC_DN,              & ! Input Optical Elastic
            OMEGAPHASFUNC_CABANNES_UP,  OMEGAPHASFUNC_CABANNES_DN,             & ! Input Optical Cabannes
            OMEGAMOMS_RRSLOSS_UNSCALED, OMEGAMOMS_RRSGAIN_UNSCALED,            & ! Input Optical Raman
            SAVE_TRANS_USERM, SAVE_BEAMMULT_UP, SAVE_BEAMMULT_DN,              & ! Input Optical for NADIR SSCORR 2p5a
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
            FAIL, MESSAGE_SUB)                                               ! Outputs

!  Old Vesion 2.3 code
!          CALL LRRS_SSCORR_OG_MONO_PLUS_1 &
!         ( DO_RRS_OVERALL, DO_UPWELLING, DO_DNWELLING, &
!           DO_CABANNES_RAMAN, DO_ENERGY_BALANCING, DO_WATERLEAVING, &
!           DO_LAMBERTIAN_SURFACE, DO_LAMBERTIAN_FRACTION, DO_DELTAM_SCALING, &
!           DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, &
!           LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS, &
!           NLAYERS, NFINELAYERS, NMOMENTS_INPUT, &
!           NPOINTS_MONO, W_EXCIT, &
!           N_LOUTPUT, LEVELMASK_UP, LEVELMASK_DN, &
!           N_USER_STREAMS, N_USER_RELAZMS, &
!           SOLAR_ANGLES_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST, &
!           LAMBERTIAN_FRACTION, EXACTDB_BRDFUNC, LS_EXACTDB_BRDFUNC, &
!           ALBEDOS_RANKED, FLUXES_RANKED, EARTH_RADIUS, HEIGHT_GRID, &
!             DELTAU_VERT_INPUT,   OMEGAMOMS_ELASTIC_UNSCALED, &
!           L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC_UNSCALED, &
!           TRUNC_FACTORS, L_TRUNC_FACTORS, &
!           OMEGAMOMS_CABANNES_UNSCALED,L_OMEGAMOMS_CABANNES_UNSCALED, &
!           OMEGAMOMS_RRSLOSS_UNSCALED, L_OMEGAMOMS_RRSLOSS_UNSCALED, &
!           OMEGAMOMS_RRSGAIN_UNSCALED, L_OMEGAMOMS_RRSGAIN_UNSCALED, &
!           ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN, &
!           LC_ELASTIC_SS_UP, LP_ELASTIC_SS_UP, LS_ELASTIC_SS_UP, &
!           LC_ELASTIC_SS_DN, LP_ELASTIC_SS_DN, &
!           LC_RAMAN_SS_UP,   LP_RAMAN_SS_UP,   LS_RAMAN_SS_UP, &
!           LC_RAMAN_SS_DN,   LP_RAMAN_SS_DN, &
!           FAIL, MESSAGE_SUB)

          IF ( FAIL ) THEN
            STATUS = LRRS_SERIOUS              ; N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = MESSAGE_SUB ; N_MESSAGES = N_MESSAGES + 1
            MESSAGES(N_MESSAGES) = 'Failure from LRRS_SSCORR_GENERAL_MONO_PLUS_1'
            LRRS_Out%Status%STATUS_CALCULATION = STATUS
            LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
            LRRS_Out%Status%MESSAGES           = MESSAGES
            RETURN
          ENDIF

!  THIS IS POSSIBILITY 2 (PLACEHOLDER, 10 March 2011)
!     REMARK, NEEDS THE ADJUSTED ANGLES (18 March 2011)
!            ***********  P L A C E H O L D E R  ***************

        ENDIF
      ENDIF

!mick chg 7/20/2016 - move this block after input arg passing
!  Npoints local

      !IF ( DO_BIN_REALIZATION ) THEN
      !  NPOINTS_TOTAL = NPOINTS_OUTER
      !  NPOINTS_LOCAL = NPOINTS_INNER
      !ELSE
      !  NPOINTS_TOTAL = NPOINTS_MONO
      !  NPOINTS_LOCAL = 1
      !ENDIF

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
!   -- Rob mod 5/12/17 for 2p5a. Renamed SSCORR_ALONE flag

      IF ( LOCAL_DO_NO_AZIMUTH .OR. DO_SSCORR_ALONE ) THEN
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
!   -- Rob mod 5/12/17 for 2p5a. Renamed SSCORR_ALONE flag

        IF ( DO_SSCORR_ALONE ) GO TO 7443

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
              AZM_ARGUMENT  = USER_RELAZMS_ADJUST(UM,UA) * DFC
              AZMFAC(UM,UA) = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
            ENDDO
          ENDDO
        ENDIF

!  Progress - could get rid of this

        IF (.NOT.DO_RRS_OVERALL) THEN
          write(*,*)' ..Elastic solution: Fourier component',FOURIER
        ELSE
          write(*,*)' ..RRS solution: Fourier component',FOURIER
        ENDIF

!  Main call to LRRS Fourier module
!    (1) Use of Supplement-derived BRDF/SLEAVE inputs and control
!    (2) Use of new Taylor series inputs and modules.
!   -- Rob mod 5/12/17  for 2p5a. Use SSCORR_GENERAL and SSCORR_ALONE flags
!   -- Rob mod 10/10/18 for 2p5a. Add DO_TIMING, NTHREADS inputs, OMP Timing Outputs

        CALL L_RAMAN_FOURIER_1 ( DO_TIMING, NTHREADS,                         & ! OMP Input
         DO_RRS_OVERALL, DO_DIRECTRRS_ONLY, DO_MSMODE_LRRS, PROGRESS,         & ! Inputs
         DO_ELASTIC_ONLY, DO_SSCORR_GENERAL, DO_SSCORR_ALONE,                 & ! Inputs
         DO_BIN_REALIZATION, DO_ENERGY_BALANCING, DO_PLANE_PARALLEL,          & ! Inputs
         DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_DIRECT_BEAM,                  & ! Inputs
         DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_BRDF_SURFACE,        & ! Inputs
         DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_BRDF_WAV1, DO_SLEAVE_WAV1,   & ! Inputs
         DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,        & ! Inputs (Linearized)
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS,               & ! Inputs (Linearized)
         N_SURFACE_WFS, N_SLEAVE_WFS,                                         & ! Inputs (Linearized)
         NPOINTS_OUTER, OFFSET_INNER, NPOINTS_INNER,                          & ! Inputs
         N_RRSBINS, BINMAP, NPOINTS_MONO, W_EXCIT, FLUXES_RANKED,             & ! Inputs
         ALBEDOS_RANKED, BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,        & ! Inputs
         SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                       & ! Inputs
         NLAYERS, NSTREAMS, NMOMENTS, N_LOUTPUT, N_USER_STREAMS,              & ! Inputs
         FOURIER, COS_SZA, TAYLOR_SMALL, TAYLOR_ORDER, FLUX_FACTOR,           & ! Inputs
         QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWGT, USER_STREAMS,              & ! Inputs
         STERM_MASK_UP, LEVELMASK_UP, STERM_MASK_DN, LEVELMASK_DN,            & ! Inputs
         DELTAU_VERT_INPUT, OMEGAMOMS_ELASTIC, OMEGAMOMS_CABANNES,            & ! Inputs
         OMEGAMOMS_RRSLOSS, OMEGAMOMS_RRSBIN,  OMEGAMOMS_RRSGAIN,             & ! Inputs
         BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT,                           & ! Inputs
         BEAM_DTRANS, BEAM_ETRANS, SAVE_TRANS_USERM,                          & ! Inputs
         LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,            & ! Inputs (linearized)
         LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,        & ! Inputs (linearized)
         L_DELTAU_VERT_INPUT, L_OMEGAMOMS_ELASTIC, L_OMEGAMOMS_CABANNES,      & ! Inputs (linearized)
         L_OMEGAMOMS_RRSLOSS, L_OMEGAMOMS_RRSBIN,  L_OMEGAMOMS_RRSGAIN,       & ! Inputs (linearized)
         NKSTORAGE, NKSTORAGE2, L_BEAM_ITRANS, L_BEAM_AVSECANT,               & ! Inputs (linearized)
         L_BEAM_ETRANS, L_BEAM_DTRANS, L_SAVE_TRANS_USERM,                    & ! Inputs (linearized)
         ELASTIC_F_UP, MEAN_ELASTIC_UP, FLUX_ELASTIC_UP,                      & ! outputs
         ELASTIC_F_DN, MEAN_ELASTIC_DN, FLUX_ELASTIC_DN,                      & ! outputs
         RAMAN_F_UP, MEAN_RAMAN_UP, FLUX_RAMAN_UP,                            & ! outputs
         RAMAN_F_DN, MEAN_RAMAN_DN, FLUX_RAMAN_DN,                            & ! outputs
         L_ELASTIC_F_UP, L_ELASTIC_F_DN, LS_ELASTIC_F_UP, LS_ELASTIC_F_DN,               & ! outputs (linearized)
         L_MEAN_ELASTIC_UP, L_FLUX_ELASTIC_UP, L_MEAN_ELASTIC_DN, L_FLUX_ELASTIC_DN,     & ! outputs (linearized)
         LS_MEAN_ELASTIC_UP, LS_FLUX_ELASTIC_UP, LS_MEAN_ELASTIC_DN, LS_FLUX_ELASTIC_DN, & ! outputs (linearized)
         L_RAMAN_F_UP, L_RAMAN_F_DN, LS_RAMAN_F_UP, LS_RAMAN_F_DN,                       & ! outputs (linearized)
         L_MEAN_RAMAN_UP, L_FLUX_RAMAN_UP, L_MEAN_RAMAN_DN, L_FLUX_RAMAN_DN,             & ! outputs (linearized)
         LS_MEAN_RAMAN_UP, LS_FLUX_RAMAN_UP, LS_MEAN_RAMAN_DN, LS_FLUX_RAMAN_DN,         & ! outputs (linearized)
         FAIL, MESSAGE_SUB, MESSAGE, POINT_TRACE,                                        & ! outputs (Exceptions)
         OMP_Elastic_Time(fourier), OMP_Raman_Time(fourier) )                              ! OMP Timing outputs

! Error handling. Add to collection of messages

        C3 = '   '
        IF ( FAIL ) THEN
          STATUS = LRRS_SERIOUS
          IF (FOURIER.LT.10)WRITE(C3(3:3),'(I1)')FOURIER
          IF (FOURIER.GE.10)WRITE(C3(2:3),'(I2)')FOURIER
          MESSAGES(N_MESSAGES+1)= 'Error from L_RAMAN_FOURIER_1, m = '//C3
          N_MESSAGES = N_MESSAGES + 1 ; MESSAGES(N_MESSAGES+1) = MESSAGE
          N_MESSAGES = N_MESSAGES + 1 ; MESSAGES(N_MESSAGES+1) = MESSAGE_SUB
          N_MESSAGES = N_MESSAGES + 1
          IF ( DO_BIN_REALIZATION ) THEN
            MESSAGES(N_MESSAGES+1) = POINT_TRACE ; N_MESSAGES = N_MESSAGES + 1
          ENDIF
          LRRS_Out%Status%STATUS_CALCULATION = STATUS
          LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
          LRRS_Out%Status%MESSAGES           = MESSAGES
          RETURN
        ENDIF

!  Continuation point for SSCORR_ALONE

 7443   continue

!  Convergence examination
!   -- Rob mod 5/12/17 for 2p5a. Use SSCORR_GENERAL and SSCORR_ALONE flags

        IF ( .not. DO_MVOUT_ONLY  ) THEN
          CALL RAMAN_CONVERGE &
           ( DO_UPWELLING, DO_DNWELLING, DO_DOUBLE_CONVTEST,                 & ! Inputs
             DO_SSCORR_GENERAL, DO_SSCORR_ALONE, DO_MOLECSCAT_ONLY,          & ! Inputs
             DO_NO_AZIMUTH, DO_RRS_OVERALL, DO_ELASTIC_ONLY, FOURIER,        & ! Inputs
             NPOINTS_LOCAL, N_LOUTPUT, N_OUT_STREAMS, LOCAL_N_USERAZM,       & ! Inputs
             N_CONVTESTS, N_GEOMETRIES, GEOM_OFFSETS, AZMFAC, LRRS_ACCURACY, & ! Inputs
             ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN,         & ! Inputs
             ELASTIC_F_UP, ELASTIC_F_DN, RAMAN_F_UP, RAMAN_F_DN,             & ! Inputs
             ELASTIC_UP, ELASTIC_DN, RAMAN_UP, RAMAN_DN,                     & ! Outputs
             FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )                        ! Outputs
        ENDIF

!  Linearizations

        IF ( .not. DO_MVOUT_ONLY  ) THEN
          IF ( DO_PROFILE_WFS .OR. DO_COLUMN_WFS .OR. &
               DO_SURFACE_WFS .OR. DO_SLEAVE_WFS )THEN
            CALL L_RAMAN_CONVERGE &
            ( DO_UPWELLING, DO_DNWELLING, DO_SSCORR_GENERAL, DO_SSCORR_ALONE,    & ! Inputs
              DO_NO_AZIMUTH, DO_RRS_OVERALL, DO_ELASTIC_ONLY, FOURIER,           & ! Inputs
              DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,      & ! Inputs
              LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS, & ! Inputs
              NPOINTS_LOCAL, N_LOUTPUT, N_OUT_STREAMS, LOCAL_N_USERAZM,          & ! Inputs
              NLAYERS, N_GEOMETRIES, GEOM_OFFSETS, AZMFAC,                       & ! Inputs
              LC_ELASTIC_SS_UP, LC_ELASTIC_SS_DN, LP_ELASTIC_SS_UP, LP_ELASTIC_SS_DN, LS_ELASTIC_SS_UP, & ! E-SS Inputs
              LC_RAMAN_SS_UP,   LC_RAMAN_SS_DN,   LP_RAMAN_SS_UP,   LP_RAMAN_SS_DN,   LS_RAMAN_SS_UP,   & ! R-SS Inputs
              L_ELASTIC_F_UP,   L_ELASTIC_F_DN,   LS_ELASTIC_F_UP,  LS_ELASTIC_F_DN,                    & ! E-Fr Inputs
              L_RAMAN_F_UP,     L_RAMAN_F_DN,     LS_RAMAN_F_UP,    LS_RAMAN_F_DN,                      & ! R-Fr Inputs
              LC_ELASTIC_UP,    LC_ELASTIC_DN,    LP_ELASTIC_UP,    LP_ELASTIC_DN,    LS_ELASTIC_UP,    LS_ELASTIC_DN, & ! E Outputs
              LC_RAMAN_UP,      LC_RAMAN_DN,      LP_RAMAN_UP,      LP_RAMAN_DN,      LS_RAMAN_UP,      LS_RAMAN_DN,   & ! R Outputs
              LOCAL_ITERATION ) !  Diagnostic output
          ENDIF
        ENDIF

!  Exit if SSFULL
!   -- Rob mod 5/12/17 for 2p5a. Renamed SSCORR_ALONE flag

        if ( DO_SSCORR_ALONE ) return

!  Fourier output
!  --------------

!  Open file if Fourier = 0
!  Write Standard Fourier output. Disabled Version 2.5.
!  Close file if iteration has finished
!   ........... linearized routine not yet available

        IF ( DO_LRRS_WRITEFOURIER ) THEN
          FUNIT = LRRS_FUNIT
          IF ( FOURIER .EQ. 0 ) THEN
            OPEN(FUNIT,FILE=DEBUG_FILENAMES(3),STATUS='UNKNOWN')
          ENDIF
!          CALL LRRS_WRITEFOURIER ( FUNIT, FOURIER )
!        IF ( DO_ATMOS_WFS ) CALL L_LRRS_WRITEFOURIER ( FUNIT )
          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
        ENDIF

!  end iteration loop

      ENDDO

!  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

!  Major result output
!  ===================

!  Standard output. Disabled for Version 2.5
!   ........... linearized routine not yet available

      IF ( DO_LRRS_WRITERESULTS ) THEN
        RUNIT = LRRS_RESUNIT
        OPEN(RUNIT,FILE=DEBUG_FILENAMES(4),STATUS='UNKNOWN')
!        CALL LRRS_WRITERESULTS ( RUNIT )
!        IF ( DO_ATMOS_WFS ) CALL L_LRRS_WRITERESULTS ( RUNIT )
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

      LRRS_Out%Status%STATUS_INPUTCHECK               = STATUS_INPUTCHECK
      LRRS_Out%Status%NCHECKMESSAGES                  = NCHECKMESSAGES
      LRRS_Out%Status%CHECKMESSAGES(0:NCHECKMESSAGES) = CHECKMESSAGES(0:NCHECKMESSAGES)
      LRRS_Out%Status%ACTIONS(0:NCHECKMESSAGES)       = ACTIONS(0:NCHECKMESSAGES)

      LRRS_Out%Status%STATUS_CALCULATION = STATUS
      LRRS_Out%Status%N_MESSAGES         = N_MESSAGES
      LRRS_Out%Status%MESSAGES           = MESSAGES

!  Linearized Atmosphere Outputs

      !Elastic Field
      LRRS_LinOut%Atmos%LC_ELASTIC_UP     = LC_ELASTIC_UP
      LRRS_LinOut%Atmos%LC_ELASTIC_DN     = LC_ELASTIC_DN
      LRRS_LinOut%Atmos%LP_ELASTIC_UP     = LP_ELASTIC_UP
      LRRS_LinOut%Atmos%LP_ELASTIC_DN     = LP_ELASTIC_DN

      LRRS_LinOut%Atmos%LC_ELASTIC_SS_UP  = LC_ELASTIC_SS_UP
      LRRS_LinOut%Atmos%LC_ELASTIC_SS_DN  = LC_ELASTIC_SS_DN
      LRRS_LinOut%Atmos%LP_ELASTIC_SS_UP  = LP_ELASTIC_SS_UP
      LRRS_LinOut%Atmos%LP_ELASTIC_SS_DN  = LP_ELASTIC_SS_DN

      LRRS_LinOut%Atmos%LC_MEAN_ELASTIC_UP(:,:,:)   = L_MEAN_ELASTIC_UP(:,0,:,:)
      LRRS_LinOut%Atmos%LC_MEAN_ELASTIC_DN(:,:,:)   = L_MEAN_ELASTIC_DN(:,0,:,:)
      LRRS_LinOut%Atmos%LC_FLUX_ELASTIC_UP(:,:,:)   = L_FLUX_ELASTIC_UP(:,0,:,:)
      LRRS_LinOut%Atmos%LC_FLUX_ELASTIC_DN(:,:,:)   = L_FLUX_ELASTIC_DN(:,0,:,:)

      LRRS_LinOut%Atmos%LP_MEAN_ELASTIC_UP(:,:,:,:)   = L_MEAN_ELASTIC_UP(:,1:MAX_LAYERS,:,:)
      LRRS_LinOut%Atmos%LP_MEAN_ELASTIC_DN(:,:,:,:)   = L_MEAN_ELASTIC_DN(:,1:MAX_LAYERS,:,:)
      LRRS_LinOut%Atmos%LP_FLUX_ELASTIC_UP(:,:,:,:)   = L_FLUX_ELASTIC_UP(:,1:MAX_LAYERS,:,:)
      LRRS_LinOut%Atmos%LP_FLUX_ELASTIC_DN(:,:,:,:)   = L_FLUX_ELASTIC_DN(:,1:MAX_LAYERS,:,:)

      !Raman Field
      LRRS_LinOut%Atmos%LC_RAMAN_UP       = LC_RAMAN_UP
      LRRS_LinOut%Atmos%LC_RAMAN_DN       = LC_RAMAN_DN
      LRRS_LinOut%Atmos%LP_RAMAN_UP       = LP_RAMAN_UP
      LRRS_LinOut%Atmos%LP_RAMAN_DN       = LP_RAMAN_DN

      LRRS_LinOut%Atmos%LC_RAMAN_SS_UP    = LC_RAMAN_SS_UP
      LRRS_LinOut%Atmos%LC_RAMAN_SS_DN    = LC_RAMAN_SS_DN
      LRRS_LinOut%Atmos%LP_RAMAN_SS_UP    = LP_RAMAN_SS_UP
      LRRS_LinOut%Atmos%LP_RAMAN_SS_DN    = LP_RAMAN_SS_DN

      LRRS_LinOut%Atmos%LC_MEAN_RAMAN_UP(:,:,:)   = L_MEAN_RAMAN_UP(:,0,:,:)
      LRRS_LinOut%Atmos%LC_MEAN_RAMAN_DN(:,:,:)   = L_MEAN_RAMAN_DN(:,0,:,:)
      LRRS_LinOut%Atmos%LC_FLUX_RAMAN_UP(:,:,:)   = L_FLUX_RAMAN_UP(:,0,:,:)
      LRRS_LinOut%Atmos%LC_FLUX_RAMAN_DN(:,:,:)   = L_FLUX_RAMAN_DN(:,0,:,:)

      LRRS_LinOut%Atmos%LP_MEAN_RAMAN_UP(:,:,:,:)   = L_MEAN_RAMAN_UP(:,1:MAX_LAYERS,:,:)
      LRRS_LinOut%Atmos%LP_MEAN_RAMAN_DN(:,:,:,:)   = L_MEAN_RAMAN_DN(:,1:MAX_LAYERS,:,:)
      LRRS_LinOut%Atmos%LP_FLUX_RAMAN_UP(:,:,:,:)   = L_FLUX_RAMAN_UP(:,1:MAX_LAYERS,:,:)
      LRRS_LinOut%Atmos%LP_FLUX_RAMAN_DN(:,:,:,:)   = L_FLUX_RAMAN_DN(:,1:MAX_LAYERS,:,:)

!  Linearized Surface Outputs

      !Elastic Field
      LRRS_LinOut%Surf%LS_ELASTIC_UP      = LS_ELASTIC_UP
      LRRS_LinOut%Surf%LS_ELASTIC_DN      = LS_ELASTIC_DN

!      LRRS_LinOut%Surf%LS_ELASTIC_SS_UP   = LS_ELASTIC_SS_UP

      LRRS_LinOut%Surf%LS_MEAN_ELASTIC_UP = LS_MEAN_ELASTIC_UP
      LRRS_LinOut%Surf%LS_MEAN_ELASTIC_DN = LS_MEAN_ELASTIC_DN
      LRRS_LinOut%Surf%LS_FLUX_ELASTIC_UP = LS_FLUX_ELASTIC_UP
      LRRS_LinOut%Surf%LS_FLUX_ELASTIC_DN = LS_FLUX_ELASTIC_DN

      !Raman Field
      LRRS_LinOut%Surf%LS_RAMAN_UP        = LS_RAMAN_UP
      LRRS_LinOut%Surf%LS_RAMAN_DN        = LS_RAMAN_DN

!      LRRS_LinOut%Surf%LS_RAMAN_SS_UP     = LS_RAMAN_SS_UP

      LRRS_LinOut%Surf%LS_MEAN_RAMAN_UP   = LS_MEAN_RAMAN_UP
      LRRS_LinOut%Surf%LS_MEAN_RAMAN_DN   = LS_MEAN_RAMAN_DN
      LRRS_LinOut%Surf%LS_FLUX_RAMAN_UP   = LS_FLUX_RAMAN_UP
      LRRS_LinOut%Surf%LS_FLUX_RAMAN_DN   = LS_FLUX_RAMAN_DN

!  ===================================
!  END COPY LOCAL VARIABLES TO OUTPUTS
!  ===================================

!Rob add - display OpenMP timing. 10/10/18
      IF ( DO_TIMING ) THEN
        WRITE(*,*)
        WRITE(*,'(7x,3(2x,a))') 'Fourier','Elastic Time','Raman Time'
        WRITE(*,'(7x,3(2x,a))') '-------','------------','----------'
        DO I=0,FOURIER
          WRITE(*,'(12x,i2,7x,f7.3,6x,f7.3)') I,OMP_Elastic_Time(I),OMP_Raman_Time(I)
          !WRITE(*,*) I,OMP_Elastic_Time(I),OMP_Raman_Time(I)
        ENDDO
      ENDIF

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE LRRS_DEBUG_INPUT_MASTER()

!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT, added PHASFUNC inputs
!   -- Rob mod 5/12/17 for 2p5a. Use SSCORR_NADIR and SSCORR_ALONE flags

        CALL LRRS_WRITE_STD_INPUT ( &
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
          GEOMETRY_SPECHEIGHT,NPOINTS_TOTAL )

      IF (DO_BRDF_SURFACE) THEN
        CALL LRRS_WRITE_SUP_BRDF_INPUT ( &
          NSTREAMS,N_USER_STREAMS,N_USER_RELAZMS,NPOINTS_TOTAL,&
          EXACTDB_BRDFUNC,BRDF_F_0,BRDF_F,USER_BRDF_F_0,USER_BRDF_F )
      END IF

      IF (DO_SURFACE_LEAVING) THEN
        CALL LRRS_WRITE_SUP_SLEAVE_INPUT ( &
          NSTREAMS,N_USER_STREAMS,N_USER_RELAZMS,NPOINTS_TOTAL,&
          SLTERM_ISOTROPIC,SLTERM_USERANGLES,SLTERM_F_0,USER_SLTERM_F_0 )
      END IF

      END SUBROUTINE LRRS_DEBUG_INPUT_MASTER

      SUBROUTINE LRRS_DEBUG_LIN_INPUT_MASTER()

!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT, added PHASFUNC inputs

      CALL LRRS_WRITE_LIN_INPUT ( &
        NLAYERS,NPOINTS_TOTAL,N_USER_STREAMS,N_USER_RELAZMS,&
        DO_COLUMN_WFS,DO_PROFILE_WFS,DO_SURFACE_WFS,DO_SLEAVE_WFS,&
        DO_NORMALIZED_WFS,DO_AIRPROFILE_WFS,DO_TEMPPROFILE_WFS,DO_TEMPSHIFT_WF,&
        LAYER_VARY_FLAG,LAYER_VARY_NUMBER,&
        N_TOTALCOLUMN_WFS,N_SURFACE_WFS,N_SLEAVE_WFS,&
        LAYER_AIRCOLUMNS_DT,TEMPERATURES_UNSHIFTED,&
        L_DELTAU_INPUT_UNSCALED,L_OMEGAMOMS_ELASTIC_UNSCALED, &
        L_OMEGAPHASFUNC_ELASTIC_UP, L_OMEGAPHASFUNC_ELASTIC_DN )

      IF (DO_BRDF_SURFACE .AND. DO_SURFACE_WFS) THEN
        CALL LRRS_WRITE_LIN_SUP_BRDF_INPUT ( &
          NSTREAMS,N_USER_STREAMS,N_USER_RELAZMS,NPOINTS_TOTAL,N_SURFACE_WFS,&
          LS_EXACTDB_BRDFUNC,LS_BRDF_F_0,LS_BRDF_F,LS_USER_BRDF_F_0,LS_USER_BRDF_F )
      END IF

      IF (DO_SURFACE_LEAVING .AND. DO_SURFACE_WFS .AND. DO_SLEAVE_WFS) THEN
        CALL LRRS_WRITE_LIN_SUP_SLEAVE_INPUT ( &
          NSTREAMS,N_USER_STREAMS,N_USER_RELAZMS,NPOINTS_TOTAL,N_SLEAVE_WFS,&
          LSSL_SLTERM_ISOTROPIC,LSSL_SLTERM_USERANGLES,&
          LSSL_SLTERM_F_0,LSSL_USER_SLTERM_F_0 )
      END IF

      END SUBROUTINE LRRS_DEBUG_LIN_INPUT_MASTER

      END SUBROUTINE L_RAMAN_MASTER

!  End module

      END MODULE lrrs_L_main_master_m

