
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

    MODULE LRRS_LinOutputs_Def_m

!    New for LRRS Version 2.5. 10/9/15, R. Spurr, RT Solutions Inc.
!       Modeled after the LIDORT Version 3.7 file. 

!  This module contains the following LRRS output structures,
!  with intents :

!         LRRS_LinAtmos      nested in LRRS_LinOutputs
!         LRRS_LinSurf       nested in LRRS_LinOutputs
!         LRRS_LinOutputs    Intent(Out)

      USE LRRS_PARS_m, only : fpk, MAX_GEOMETRIES, MAX_USER_STREAMS, MAX_LOUTPUT, &
                              MAX_POINTS, MAX_LAYERS, MAX_ATMOSWFS, MAX_SURFACEWFS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LRRS_LinAtmos


!  ELASTIC FIELD
!  -------------

!  Atmospheric Column Jacobians

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LC_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LC_ELASTIC_DN

!  Atmospheric Profile Jacobians

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LP_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LP_ELASTIC_DN

!  SS output

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LC_ELASTIC_SS_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LC_ELASTIC_SS_DN
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LP_ELASTIC_SS_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LP_ELASTIC_SS_DN

!  Mean-value output

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS) :: LP_MEAN_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS) :: LP_MEAN_ELASTIC_DN

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS) :: LP_FLUX_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS) :: LP_FLUX_ELASTIC_DN

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_POINTS) :: LC_MEAN_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_POINTS) :: LC_MEAN_ELASTIC_DN

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_POINTS) :: LC_FLUX_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_POINTS) :: LC_FLUX_ELASTIC_DN

!  RAMAN FIELD
!  -----------

!  Atmospheric Column Jacobians

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LC_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LC_RAMAN_DN

!  Atmospheric Profile Jacobians

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LP_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LP_RAMAN_DN

!  SS output

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LC_RAMAN_SS_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LC_RAMAN_SS_DN
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LP_RAMAN_SS_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LP_RAMAN_SS_DN

!  Mean-value output

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS) :: LP_MEAN_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS) :: LP_MEAN_RAMAN_DN

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS) :: LP_FLUX_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LAYERS, MAX_LOUTPUT, MAX_POINTS) :: LP_FLUX_RAMAN_DN

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_POINTS) :: LC_MEAN_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_POINTS) :: LC_MEAN_RAMAN_DN

      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_POINTS) :: LC_FLUX_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_ATMOSWFS, MAX_LOUTPUT, MAX_POINTS) :: LC_FLUX_RAMAN_DN

      END TYPE LRRS_LinAtmos

! #####################################################################
! #####################################################################

      TYPE LRRS_LinSurf

!  ELASTIC FIELD
!  -------------

!  Surface Jacobians

      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LS_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LS_ELASTIC_DN

!  Mean-value output

      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS) :: LS_MEAN_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS) :: LS_MEAN_ELASTIC_DN
      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS) :: LS_FLUX_ELASTIC_UP
      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS) :: LS_FLUX_ELASTIC_DN

!  RAMAN FIELD
!  -----------

!  Surface Jacobians

      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LS_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS) :: LS_RAMAN_DN

!  Mean-value output

      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS) :: LS_MEAN_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS) :: LS_MEAN_RAMAN_DN
      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS) :: LS_FLUX_RAMAN_UP
      REAL(FPK), DIMENSION (MAX_SURFACEWFS, MAX_LOUTPUT, MAX_POINTS) :: LS_FLUX_RAMAN_DN

      END TYPE LRRS_LinSurf

! #####################################################################
! #####################################################################

      TYPE LRRS_LinOutputs

      TYPE(LRRS_LinAtmos)  :: Atmos
      TYPE(LRRS_LinSurf)   :: Surf

      END TYPE LRRS_LinOutputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LRRS_LinAtmos,&
                LRRS_LinSurf,&
                LRRS_LinOutputs

      end module LRRS_LinOutputs_def_m

