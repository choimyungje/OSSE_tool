
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

!  Developed for LRRS Version 2.5, 9/8/15. R. Spurr, RT SOLUTIONS Inc.
!   Closely follows the LIDORT module with the same name, but
!     (1) No solar angle dimensioning, No observational geometry.
!     (2) Additional wavelength dimensioning (MAX_SLEAVE_POINTS)
!     (3) Additional control for using 1 or all wavelengths

!  Note (9/8/15). The number of wavelengths for SLEAVE has deliberately
!  been left flexible - typically the SLEAVE properties will change very
!  little over a Raman-scattering window (+/- 2 nm in the UV), so to
!  a very good approximation, it is often sufficient to use one point for
!  the calculations of SLEAVE - in which case, MAX_SLEAVE_POINTS will be 1
!  This applies equally to Waterleaving and Fluorescence, where wavelength
!  variation can be important in the Visible and NIR.

!  Now, whenever the SLEAVE supplement is used with LRRS, the choice of 
!  SLEAVE wavelengths is linked to the Raman wavelengths (LAMDAS_RANKED).
!  These wavelengths are now input to the SLEAVE supplement MASTER, and
!  they are not set by hand or by configuration-file read. 

!  However, there is a hard-wired choice (DO_WAV1) which allows you
!  to choose a single wavelength point for calculation. When this flag is
!  set, the number of SLEAVE wavelengths = 1 (Waterleaving or Fluorescence)
!  and the single wavelength is set to the AVERAGE VALUE of LAMBDAS_RANKED.

!  If you want all wavelengths, then control flag DO_Wav1 is False,
!  and number of SLEAVE wavelengths  N_Lambdas_Ranked, and the
!  wavelengths themselves are just copied from LAMBDAS_RANKED
!  A SLEAVE calculation will be done for all Raman-scattered wavelengths!!!

module SLEAVE_LinSup_Outputs_def_m

!  This module contains the following structures:

!  SLEAVE_LinSup_Outputs - Intent(In)  for LRRS,
!                          Intent(Out) for SLEAVE_LinSup

use LRRS_PARS_m, Only : fpk, MAX_SLEAVEWFS, MAX_MOMENTS, MAX_STREAMS, &
                        MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_SLEAVE_POINTS

implicit none

! #####################################################################
! #####################################################################

type SLEAVE_LinSup_Outputs

!  Isotropic Surface leaving term (if flag set), derivatives

      REAL(fpk), dimension ( MAX_SLEAVEWFS, MAX_SLEAVE_POINTS ) :: SL_LS_SLTERM_ISOTROPIC

!  Exact Surface-Leaving term

      REAL(fpk), dimension ( MAX_SLEAVEWFS, MAX_USER_STREAMS, &
        MAX_USER_RELAZMS, MAX_SLEAVE_POINTS ) :: SL_LS_SLTERM_USERANGLES

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk), dimension ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_STREAMS, &
        MAX_SLEAVE_POINTS ) :: SL_LS_SLTERM_F_0
      REAL(fpk), dimension ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, &
        MAX_SLEAVE_POINTS ) :: SL_LS_USER_SLTERM_F_0

end type SLEAVE_LinSup_Outputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: SLEAVE_LinSup_Outputs

end module SLEAVE_LinSup_Outputs_def_m
