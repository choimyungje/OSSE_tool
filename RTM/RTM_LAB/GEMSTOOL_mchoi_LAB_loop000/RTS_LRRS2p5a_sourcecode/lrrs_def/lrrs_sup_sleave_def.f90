
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

      MODULE LRRS_Sup_SLEAVE_def_m

!   Newly Introduced to LRRS Version 2.5. 9/9/15, R. Spurr, RT Solutions Inc.
!       Modeled after the LIDORT Version 3.7 file. Differences:--

!         (1) No solar angle dimensioning
!         (2) Additional wavelength dimensioning (MAX_SLEAVE_POINTS)

!  This module contains the following structures:
!     LRRS_Sup_SLEAVE    Intent(In)  for LRRS
!                        Intent(Out) for LRRS SleaveSup

      USE LRRS_PARS_m, only : fpk, MAX_USER_STREAMS, MAX_USER_RELAZMS, &
                              MAX_MOMENTS, MAX_STREAMS, MAX_POINTS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LRRS_Sup_SLEAVE

!  Use of Single-Wavelength flag

      LOGICAL :: DO_SLEAVE_Wav1

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), dimension ( MAX_POINTS ) :: SLTERM_ISOTROPIC

!  Exact Surface-Leaving term

      REAL(fpk), dimension ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS ) :: SLTERM_USERANGLES

!  Fourier components of Surface-leaving terms:
!    solar direction, SL-transmitted quadrature streams
!    solar direction, SL-transmitted user streams

      REAL(fpk), dimension ( 0:MAX_MOMENTS, MAX_STREAMS,      MAX_POINTS ) :: SLTERM_F_0
      REAL(fpk), dimension ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS ) :: USER_SLTERM_F_0

      END TYPE LRRS_Sup_SLEAVE

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LRRS_Sup_SLEAVE

      END MODULE LRRS_Sup_SLEAVE_def_m
