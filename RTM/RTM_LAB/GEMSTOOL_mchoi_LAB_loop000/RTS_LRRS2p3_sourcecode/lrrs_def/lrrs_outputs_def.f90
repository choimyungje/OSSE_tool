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

      MODULE LRRS_Outputs_Def

!  This Module contains the following LIDORT_RRS Output Structures,
!  with Intents :

!               LRRS_Outputs    Intent(Out)

      USE LRRS_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LRRS_Main_Outputs


!  ELASTIC FIELD
!  -------------

!   Radiances including only elastic scattering

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: &
        ELASTIC_UP

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: &
        ELASTIC_DN

!  Single scatter results

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: &
        ELASTIC_SS_UP

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: &
        ELASTIC_SS_DN

!  Mean-value output

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: &
        MEAN_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: &
        MEAN_ELASTIC_DN

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: &
        FLUX_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: &
        FLUX_ELASTIC_DN


!  RAMAN FIELD
!  -----------

!  Radiance including both elastic and inelastic scattering

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: &
        RAMAN_UP

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: &
        RAMAN_DN

!  Single scatter results

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: &
        RAMAN_SS_UP

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: &
        RAMAN_SS_DN

!  Mean-value output

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: MEAN_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: MEAN_RAMAN_DN

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: FLUX_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: FLUX_RAMAN_DN

!  Number of Fourier terms used

      INTEGER :: FOURIER_SAVED


      END TYPE LRRS_Main_Outputs

! #####################################################################
! #####################################################################

      TYPE LRRS_Exception_Handling


!  LRRS calculation status code

      INTEGER :: STATUS_CALCULATION

!  Number of error messages

      INTEGER :: N_MESSAGES

!  Error message(s)

      CHARACTER (LEN=70), DIMENSION (MAX_MESSAGES) :: MESSAGES


      END TYPE LRRS_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE LRRS_Outputs


      TYPE(LRRS_Main_Outputs)       :: Main
      TYPE(LRRS_Exception_Handling) :: Status


      END TYPE LRRS_Outputs

! #####################################################################
! #####################################################################

      PRIVATE
      PUBLIC :: LRRS_Main_Outputs,&
                LRRS_Exception_Handling,&
                LRRS_Outputs

      END MODULE LRRS_Outputs_Def
