
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

      module LRRS_Sup_SS_def_m

!   Newly Introduced to LRRS Version 2.5. 9/9/15, R. Spurr, RT Solutions Inc.
!       Modeled after the LIDORT Version 3.7 file. Differences:--

!         (1) No solar angle dimensioning
!         (2) Additional wavelength dimensioning (MAX_POINTS)
!         (3) Raman and Elastic fields

!  This module contains the following structures,
!  with intents :

!      LRRS_Sup_SS_def  Intent(InOut)

      USE LRRS_PARS_m, only : fpk, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS

      implicit none

! #####################################################################
! #####################################################################

      type LRRS_Sup_SS

!  Elastic SS Intensity Results at all angles, levels and points

      REAL(FPK) :: ELASTIC_SS_UP ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: ELASTIC_SS_DN ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Elastic DB Intensity Results at all angles, levels and points

      REAL(FPK) :: ELASTIC_DB_UP ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Raman SS Intensity Results at all angles, levels and points

      REAL(FPK) :: RAMAN_SS_UP ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK) :: RAMAN_SS_DN ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Raman DB Intensity Results at all angles, levels and points

      REAL(FPK) :: RAMAN_DB_UP ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

      end type LRRS_Sup_SS

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LRRS_Sup_SS

      end module LRRS_Sup_SS_def_m

