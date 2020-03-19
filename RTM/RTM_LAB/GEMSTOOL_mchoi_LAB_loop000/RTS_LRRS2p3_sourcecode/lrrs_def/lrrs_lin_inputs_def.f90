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

      MODULE LRRS_Linearized_Inputs_Def

!  This Module contains the following LIDORT_RRS Input Structures,
!  with Intents :

!          LRRS_Linearized_Inputs    Intent(In)

      USE LRRS_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LRRS_Lin_Boolean


!  Profile/column/surface linearization control inputs

      LOGICAL :: DO_PROFILE_WFS
      LOGICAL :: DO_COLUMN_WFS
      LOGICAL :: DO_SURFACE_WFS
      LOGICAL :: DO_LAMBERTIAN_WF

!  Flag for generating normalized weighting functions

      LOGICAL :: DO_NORMALIZED_WFS

!  Flag for generating air column weighting functions

      LOGICAL :: DO_AIRPROFILE_WFS

!  Flag for generating temperature profile weighting functions

      LOGICAL :: DO_TEMPPROFILE_WFS

!  Flag for generating temperature shift weighting functions

      LOGICAL :: DO_TEMPSHIFT_WF


      END TYPE LRRS_Lin_Boolean

! #####################################################################
! #####################################################################

      TYPE LRRS_Lin_Atmos


!  Total number of column weighting functions

      INTEGER :: N_TOTALCOLUMN_WFS

!  Control for atmospheric linearizations, layer by layer

      LOGICAL, DIMENSION (MAX_LAYERS) :: LAYER_VARY_FLAG
      INTEGER, DIMENSION (MAX_LAYERS) :: LAYER_VARY_NUMBER

!  Input layer air columns, should be in mol/cm^2 or [DU]

      REAL(FPK), DIMENSION (MAX_LAYERS) :: LAYER_AIRCOLUMNS_DT

!  Input layer temperatures (unshifted), must be in deg K

      REAL(FPK), DIMENSION (MAX_LAYERS) :: TEMPERATURES_UNSHIFTED

!  Linearized unscaled optical properties, elastic input

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LAYERS, MAX_POINTS) &
        :: L_DELTAU_INPUT_UNSCALED

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LAYERS, 0:MAX_MOMENTS_INPUT, &
        MAX_POINTS) :: L_OMEGAMOMS_ELASTIC_UNSCALED


      END TYPE LRRS_Lin_Atmos

! #####################################################################
! #####################################################################

      TYPE LRRS_Lin_Surf


!  Index of BRDF parameter for which surface weighting function is desired

      INTEGER :: BRDFPAR_DERIV_INDEX


      END TYPE LRRS_Lin_Surf

! #####################################################################
! #####################################################################

      TYPE LRRS_Linearized_Inputs


      TYPE(LRRS_Lin_Boolean) :: Bool
      TYPE(LRRS_Lin_Atmos)   :: Atmos
      TYPE(LRRS_Lin_Surf)    :: Surf


      END TYPE LRRS_Linearized_Inputs

! #####################################################################
! #####################################################################

      PRIVATE
      PUBLIC :: LRRS_Lin_Boolean,&
                LRRS_Lin_Atmos,&
                LRRS_Lin_Surf,&
                LRRS_Linearized_Inputs

      END MODULE LRRS_Linearized_Inputs_Def
