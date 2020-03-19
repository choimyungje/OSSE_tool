
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

      module BRDF_Sup_Outputs_def_m

!  Developed for LRRS Version 2.5, 9/8/15. R. Spurr, RT SOLUTIONS Inc.
!   Closely follows the LIDORT module with the same name, but
!     (1) No surface emission, solar angle dimensioning, No observational geometry.
!     (2) Additional wavelength dimensioning (MAX_BRDF_POINTS)
!     (3) Additional control for using 1 or all wavelengths

!  Note (9/8/15). The number of wavelengths for BRDF has deliberately
!  been left flexible - typically the BRDF properties will change very
!  little over a Raman-scattering window (+/- 2 nm in the UV), so to
!  a very good approximation, it is sufficient to use one point for
!  the calculations of BRDF - in which case, MAX_BRDF_POINTS will be 1

!  It turns out that all but one of the BRDF Kernels are wavelength
!  independent, so this consideration only applies to the "NewCM" GLINT
!  BRDF calculation, which has a wavelength input for foam reflectance
!  and refractive index through the ocean water-leaving part.

!  Now, whenever the BRDF supplement is used with LRRS, the choice of 
!  BRDF wavelengths is linked to the Raman wavelengths (LAMDAS_RANKED).
!  These wavelengths are now input to the BRDF supplement MASTER, and
!  they are not set by hand or by configuration-file read. 

!  When NewCM Glint is turned on, and you want just one wavelength, then
!  the control flag DO_NewCM_Wav1 will be set True.  In this case, the
!  number of BRDF "NewCM" wavelengths NewCM_Npoints = 1, and the single
!  wavelength NewCM_Lambdas(1) = AVERAGE VALUE of LAMBDAS_RANKED.

!  When NewCM Glint is turned on, and you want all wavelengths, then
!  control flag DO_NewCM_Wav1 is False, and NewCM_Npoints = N_Lambdas_Ranked
!  and NewCM_Lambdas(1:NewCM_npoints) = LAMBDAS_RANKED(1:NewCM_npoints).
!  A Glint calculation will be done for all Raman-scattered wavelengths!!!

! #####################################################################
! #####################################################################

!  This module contains the following structures:

!  BRDF_Sup_Outputs - Intent(In) for LRRS,
!                     Intent(Out) for BRDF_Sup

      use LRRS_PARS_m, only : fpk, MAX_MOMENTS, MAX_STREAMS, MAX_BRDF_POINTS, &
                              MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_MESSAGES

      implicit none

! #####################################################################
! #####################################################################

      type BRDF_Sup_Outputs

!  Exact (direct bounce) BRDF
!  mick fix 12/29/2014 - Name changed from EXACTDB --> DBOUNCE

      REAL(fpk), dimension ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_BRDF_POINTS ) :: BS_DBOUNCE_BRDFUNC

!  Fourier components of BRDF, in the following order
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), dimension ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_BRDF_POINTS )                   :: BS_BRDF_F_0
      REAL(fpk), dimension ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_BRDF_POINTS )      :: BS_BRDF_F
      REAL(fpk), dimension ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_BRDF_POINTS )              :: BS_USER_BRDF_F_0
      REAL(fpk), dimension ( 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_BRDF_POINTS ) :: BS_USER_BRDF_F

!  Revision 12 august 2014. Add WSA and BSA output (one SZA only for the latter)
!  11/27/14. Added the Kernel output

      REAL(fpk) :: BS_WSA_CALCULATED(MAX_BRDF_POINTS), BS_WSA_KERNELS(3, MAX_BRDF_POINTS)
      REAL(fpk) :: BS_BSA_CALCULATED(MAX_BRDF_POINTS), BS_BSA_KERNELS(3, MAX_BRDF_POINTS)

      end type BRDF_Sup_Outputs

! #####################################################################
! #####################################################################

      TYPE BRDF_Input_Exception_Handling

!  Exception handling for Input Checking settings. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: BS_STATUS_INPUTREAD
      INTEGER      :: BS_NINPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_INPUTMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_INPUTACTIONS

      END TYPE BRDF_Input_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE BRDF_Output_Exception_Handling

!  Exception handling for Output. New code, 02 April 2014
!     Message Length should be at least 120 Characters

      INTEGER      :: BS_STATUS_OUTPUT
      INTEGER      :: BS_NOUTPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_OUTPUTMESSAGES

      END TYPE BRDF_Output_Exception_Handling

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: BRDF_Sup_Outputs,              &
                BRDF_Input_Exception_Handling, &
                BRDF_Output_Exception_Handling

      end module BRDF_Sup_Outputs_def_m
