
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

      MODULE LRRS_LinInputs_Def_m

!    New for LRRS Version 2.5. 10/9/15, R. Spurr, RT Solutions Inc.
!       Modeled after the LIDORT Version 3.7 file.

!  This Module contains the following LRRS Input Structures,
!  with Intents :

!  This Module contains the following LRRS Input Structures, with Intents :

!          LRRS_LinControl    nested in LRRS_LinInputs
!          LRRS_LinOptical    nested in LRRS_LinInputs
!          LRRS_LinInputs     Intent(In)

      USE LRRS_PARS_m

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LRRS_LinControl

!  Profile/column/surface linearization control inputs

      LOGICAL :: DO_PROFILE_WFS
      LOGICAL :: DO_COLUMN_WFS
      LOGICAL :: DO_SURFACE_WFS
      LOGICAL :: DO_SLEAVE_WFS

!  Flag for generating normalized weighting functions

      LOGICAL :: DO_NORMALIZED_WFS

!  Flag for generating air column-density weighting functions

      LOGICAL :: DO_AIRPROFILE_WFS

!  Flag for generating temperature profile weighting functions

      LOGICAL :: DO_TEMPPROFILE_WFS

!  Flag for generating temperature shift weighting functions

      LOGICAL :: DO_TEMPSHIFT_WF

!  Total number of column weighting functions

      INTEGER :: N_TOTALCOLUMN_WFS

!  Total number of Surface/Sleave Jacobians
!     New for Version 2.5, modeled after LIDORT

      INTEGER  :: N_SURFACE_WFS
      INTEGER  :: N_SLEAVE_WFS

!  Control for atmospheric linearizations, layer by layer

      LOGICAL, DIMENSION (MAX_LAYERS) :: LAYER_VARY_FLAG
      INTEGER, DIMENSION (MAX_LAYERS) :: LAYER_VARY_NUMBER

      END TYPE LRRS_LinControl

! #####################################################################
! #####################################################################

      TYPE LRRS_LinOptical

!  Input layer air columns, should be in mol/cm^2 or [DU]

      REAL(FPK), DIMENSION (MAX_LAYERS) :: LAYER_AIRCOLUMNS_DT

!  Input layer temperatures (unshifted), must be in deg K

      REAL(FPK), DIMENSION (MAX_LAYERS) :: TEMPERATURES_UNSHIFTED

!  Linearized unscaled optical properties, elastic input
!    -- Rob Mod 5/12/17 for 2p5a, change dimension to MAX_MOMENTS

      REAL(FPK), DIMENSION ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS)                :: L_DELTAU_INPUT_UNSCALED
      REAL(FPK), DIMENSION ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS) :: L_OMEGAMOMS_ELASTIC_UNSCALED

!  Product of unscaled single-scatter albedo and linearized Phase functions (elastic input)
!   alternative to use of the above expansion coefficient products in single-scatter codes.
!    -- Rob Mod 5/12/17, introduced for 2p5a

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS ) :: L_OMEGAPHASFUNC_ELASTIC_UP
      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS ) :: L_OMEGAPHASFUNC_ELASTIC_DN

      END TYPE LRRS_LinOptical

! #####################################################################
! #####################################################################

      TYPE LRRS_LinInputs

      TYPE(LRRS_LinControl)    :: Cont
      TYPE(LRRS_LinOptical)    :: Optical

      END TYPE LRRS_LinInputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LRRS_LinControl,&
                LRRS_LinOptical,&
                LRRS_LinInputs

      END MODULE LRRS_LinInputs_Def_m

