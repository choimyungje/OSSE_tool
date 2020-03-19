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

      MODULE LRRS_Linearized_Outputs_Def

!  This Module contains the following LIDORT_RRS Output Structures,
!  with Intents :

!        LRRS_Linearized_Outputs    Intent(Out)

      USE LRRS_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LRRS_Lin_Atmos_Outputs


!  ELASTIC FIELD
!  -------------

!  Atmospheric Column Jacobians

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LOUTPUT, MAX_GEOMETRIES, &
        MAX_POINTS) :: LC_ELASTIC_UP

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LOUTPUT, MAX_GEOMETRIES, &
        MAX_POINTS) :: LC_ELASTIC_DN

!  Atmospheric Profile Jacobians

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LAYERS, MAX_LOUTPUT, &
        MAX_GEOMETRIES, MAX_POINTS) :: LP_ELASTIC_UP

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LAYERS, MAX_LOUTPUT, &
        MAX_GEOMETRIES, MAX_POINTS) :: LP_ELASTIC_DN

!  Single scatter results

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LOUTPUT, MAX_GEOMETRIES, &
        MAX_POINTS) :: LC_ELASTIC_SS_UP

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LOUTPUT, MAX_GEOMETRIES, &
        MAX_POINTS) :: LC_ELASTIC_SS_DN

      REAL(FPK), DIMENSION ( MAX_PARS, MAX_LAYERS, MAX_LOUTPUT, &
        MAX_GEOMETRIES, MAX_POINTS) :: LP_ELASTIC_SS_UP

      REAL(FPK), DIMENSION ( MAX_PARS, MAX_LAYERS, MAX_LOUTPUT, &
        MAX_GEOMETRIES, MAX_POINTS) :: LP_ELASTIC_SS_DN

!  Mean-value output

      REAL(FPK), DIMENSION (MAX_PARS, 0:MAX_LAYERS, MAX_LOUTPUT, &
        MAX_POINTS) :: L_MEAN_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_PARS, 0:MAX_LAYERS, MAX_LOUTPUT, &
        MAX_POINTS) :: L_MEAN_ELASTIC_DN

      REAL(FPK), DIMENSION (MAX_PARS, 0:MAX_LAYERS, MAX_LOUTPUT, &
        MAX_POINTS) :: L_FLUX_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_PARS, 0:MAX_LAYERS, MAX_LOUTPUT, &
        MAX_POINTS) :: L_FLUX_ELASTIC_DN


!  RAMAN FIELD
!  -----------

!  Atmospheric Column Jacobians

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LOUTPUT, MAX_GEOMETRIES, &
        MAX_POINTS) :: LC_RAMAN_UP

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LOUTPUT, MAX_GEOMETRIES, &
        MAX_POINTS) :: LC_RAMAN_DN

!  Atmospheric Profile Jacobians

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LAYERS, MAX_LOUTPUT, &
        MAX_GEOMETRIES, MAX_POINTS) :: LP_RAMAN_UP

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LAYERS, MAX_LOUTPUT, &
        MAX_GEOMETRIES, MAX_POINTS) :: LP_RAMAN_DN

!  Single scatter results

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LOUTPUT, MAX_GEOMETRIES, &
        MAX_POINTS) :: LC_RAMAN_SS_UP

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LOUTPUT, MAX_GEOMETRIES, &
        MAX_POINTS) :: LC_RAMAN_SS_DN

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LAYERS, MAX_LOUTPUT, &
        MAX_GEOMETRIES, MAX_POINTS) :: LP_RAMAN_SS_UP

      REAL(FPK), DIMENSION (MAX_PARS, MAX_LAYERS, MAX_LOUTPUT, &
        MAX_GEOMETRIES, MAX_POINTS) :: LP_RAMAN_SS_DN

!  Mean-value output

      REAL(FPK), DIMENSION (MAX_PARS, 0:MAX_LAYERS, MAX_LOUTPUT, &
        MAX_POINTS) :: L_MEAN_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_PARS, 0:MAX_LAYERS, MAX_LOUTPUT, &
        MAX_POINTS) :: L_MEAN_RAMAN_DN

      REAL(FPK), DIMENSION (MAX_PARS, 0:MAX_LAYERS, MAX_LOUTPUT, &
        MAX_POINTS) :: L_FLUX_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_PARS, 0:MAX_LAYERS, MAX_LOUTPUT, &
        MAX_POINTS) :: L_FLUX_RAMAN_DN


      END TYPE LRRS_Lin_Atmos_Outputs

! #####################################################################
! #####################################################################

      TYPE LRRS_Lin_Surf_Outputs


!  ELASTIC FIELD
!  -------------

!  Surface Jacobians

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) &
        :: LS_ELASTIC_UP

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) &
        :: LS_ELASTIC_DN

!  Single scatter results

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) &
        :: LS_ELASTIC_SS_UP

!  Mean-value output

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: LS_MEAN_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: LS_MEAN_ELASTIC_DN
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: LS_FLUX_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: LS_FLUX_ELASTIC_DN


!  RAMAN FIELD
!  -----------

!  Surface Jacobians

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) &
        :: LS_RAMAN_UP

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) &
        :: LS_RAMAN_DN

!  Single scatter results

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) &
        :: LS_RAMAN_SS_UP

!  Mean-value output

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: LS_MEAN_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: LS_MEAN_RAMAN_DN
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: LS_FLUX_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: LS_FLUX_RAMAN_DN


      END TYPE LRRS_Lin_Surf_Outputs

! #####################################################################
! #####################################################################

      TYPE LRRS_Linearized_Outputs


      TYPE(LRRS_Lin_Atmos_Outputs) :: Atmos
      TYPE(LRRS_Lin_Surf_Outputs)  :: Surf


      END TYPE LRRS_Linearized_Outputs

! #####################################################################
! #####################################################################

      PRIVATE
      PUBLIC :: LRRS_Lin_Atmos_Outputs,&
                LRRS_Lin_Surf_Outputs,&
                LRRS_Linearized_Outputs

      END MODULE LRRS_Linearized_Outputs_Def
