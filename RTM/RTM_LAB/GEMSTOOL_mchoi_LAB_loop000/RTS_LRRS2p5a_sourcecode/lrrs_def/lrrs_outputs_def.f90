
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

      MODULE LRRS_Outputs_Def_m

!  9/9/15 LRRS Version 2.5. Unchanged from Version 2.3

!  This Module contains the following LIDORT_RRS Output Structures,
!  with Intents :

!               LRRS_Outputs    Intent(Out)

      USE LRRS_PARS_m

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LRRS_Main_Outputs


!  ELASTIC FIELD
!  -------------

!   Radiances including only elastic scattering

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: ELASTIC_DN

!  Single scatter results

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: ELASTIC_SS_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: ELASTIC_SS_DN

!  Mean-value output

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: MEAN_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: MEAN_ELASTIC_DN

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: FLUX_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_POINTS) :: FLUX_ELASTIC_DN

! Rob fix 13/7/16, 28/11/16, 23/1/17.Output Lospaths/Sunpaths
!   -- 23/1/17. Los paths with geometry

      real(fpk) :: lospaths_up_out(max_user_streams, max_layers)
      real(fpk) :: sunpaths_up_out(max_layers)
      real(fpk) :: lospaths_dn_out(max_user_streams, max_layers)
      real(fpk) :: sunpaths_dn_out(max_layers)

!  RAMAN FIELD
!  -----------

!  Radiance including both elastic and inelastic scattering

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: RAMAN_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: RAMAN_DN

!  Single scatter results

      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: RAMAN_SS_UP
      REAL(FPK), DIMENSION (MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: RAMAN_SS_DN

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


!  Exception handling for Input Checking
!     Message lengths should be at least 120 Characters

      INTEGER :: STATUS_INPUTCHECK
      INTEGER :: NCHECKMESSAGES

      CHARACTER(Len=120), DIMENSION(0:MAX_MESSAGES) :: CHECKMESSAGES
      CHARACTER(Len=120), DIMENSION(0:MAX_MESSAGES) :: ACTIONS

!  Exception handling for Model Calculation

      INTEGER :: STATUS_CALCULATION

      INTEGER :: N_MESSAGES
      CHARACTER (Len=120), DIMENSION (MAX_MESSAGES) :: MESSAGES

      !CHARACTER(Len=120)  :: MESSAGE, TRACE_1, TRACE_2, TRACE_3


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

      END MODULE LRRS_Outputs_Def_m
