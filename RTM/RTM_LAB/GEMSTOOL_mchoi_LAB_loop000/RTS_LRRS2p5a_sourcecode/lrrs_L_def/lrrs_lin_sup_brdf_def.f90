
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

     module LRRS_Lin_Sup_BRDF_def_m

!   Newly Introduced to LRRS Version 2.5. 10/9/15, R. Spurr, RT Solutions Inc.
!       Modeled after the LIDORT Version 3.7 file. Differences:--

!         (1) No solar angle dimensioning
!         (2) Additional wavelength dimensioning (MAX_POINTS)
!         (3) Emissivity removed

!  This module contains the following Structures

!      LRRS_Lin_Sup_BRDF      Intent(In)  for LRRS,
!                             Intent(Out) for LRRS LinBRDFSup

      use LRRS_PARS_m, only : fpk, MAX_MOMENTS, MAX_STREAMS, MAX_POINTS, &
                              MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS

      implicit none

! #####################################################################
! #####################################################################

      type LRRS_LinSup_BRDF

!  Linearized Exact (direct bounce) BRDF (same all threads)

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS ) :: LS_EXACTDB_BRDFUNC

!  Fourier components of BRDF, in the following order
!    incident solar direction,    reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar direction,    reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), dimension ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )                   :: LS_BRDF_F_0
      REAL(fpk), dimension ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )      :: LS_BRDF_F
      REAL(fpk), dimension ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )              :: LS_USER_BRDF_F_0
      REAL(fpk), dimension ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS ) :: LS_USER_BRDF_F

      end type LRRS_LinSup_BRDF

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LRRS_LinSup_BRDF

      end module LRRS_Lin_Sup_BRDF_def_m

