
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

      module BRDF_LinSup_Inputs_def_m

!  Version 3.6 Notes
!  -----------------

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!     New OG inputs are :
!       Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control
!     Added Overall Exact flag for better control

! #####################################################################
! #####################################################################

!  BRDF Upgrades for Version 3.7
!  -----------------------------

!  A. White-sky and Black-sky scaling options
!  ==========================================

!  WSA and BSA scaling options.
!   first introduced 14 April 2014, Revised, 17-22 April 2014
!      WSA = White-sky albedo. BSA = Black-sky albedo.

!  These options are mutually exclusive. If either is set, the BRDF code
!  will automatically perform an albedo calculation (either WS or BS) for
!  the (1,1) component of the complete 3-kernel BRDF, then normalize the
!  entire BRDF with this albedo before scaling up with the externally chosen
!  WSA or BSA (from input).

!  Additional exception handling has been introduced to make sure that
!  the spherical or planar albedos for a complete 3-kernel BRDF are in
!  the [0,1] range. Otherwise, the WSA/BSA scaling makes no sense. It is
!  still the case that some of the MODIS-tpe kernels give NEGATIVE albedos.

!  The albedo scaling process has been linearized for all existing surface
!  property Jacobians available to the BRDF linearized supplement. In
!  addition, it is also possible to derive single Jacobians of the BRDFs
!  with respect to the WS or BS albedo - this is separate from (and orthogonal
!  to) the usual kernel derivatives.

!  B. Alternative Cox-Munk Glint Reflectance
!  =========================================

!  In conjunction with new Water-leaving code developed for the SLEAVE
!  supplement, we have given the BRDF supplement a new option to 
!  return the (scalar) Cox-Munk glint reflectance, based on code originally
!  written for the 6S code.

!  Developed and tested by R. Spurr, 21-29  April 2014
!  Based in part on Modified-6S code by A. Sayer (NASA-GSFC).
!  Validated against Modified-6S OCEABRDF.F code, 24-28 April 2014.

!  The new glint option depends on Windspeed/direction, with refractive
!  indices now computed using salinity and wavelength. There is now an 
!  optional correction for (Foam) Whitecaps (Foam). These choices come from
!  the 6S formulation.

!  Need to make sure that the wind input information is the same as that
!  use for glint calculations in the SLEAVE supplement when the Glitter
!  kernels are in use. Also, the Foam correction applied here in the
!  surface-leaving code should also be applied in the SLEAVE system..

!  Choosing this new glint option bypasses the normal kernel inputs (and
!  also the WSA/BSA inputs). Instead, a single-kernel glint reflectance
!  is calculated with amplitude 1.0, based on a separate set of dedicated
!  inputs. This option does not apply for surface emission - the solar
!  sources flag must be turned on. There is only 1 Jacobian available - 
!  with respect to the windspeed. 

!  Note that the use of the facet isotropy flag is recommended for
!  multi-beam runs. This is because the wind-direction is a function of
!  the solar angle, and including this wind-direction in the glint
!  calculations for BRDF will only work if there is just one SZA. This
!  condition is checked for internally.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! #####################################################################
! #####################################################################

!  This module contains the following structures:

!  BRDF_LinSup_Inputs - Intent(In) for BRDF_LinSup

      use LRRS_PARS_m, only : fpk, MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_BRDF_POINTS, &
                              MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS

      implicit none

! #####################################################################
! #####################################################################

      type BRDF_LinSup_Inputs

!  Linearizaion material
!  ---------------------

!  Flags for WF of bidirectional function parameters and factors

      LOGICAL, dimension ( MAX_BRDF_KERNELS ) :: &
        BS_DO_KERNEL_FACTOR_WFS
      LOGICAL, dimension ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS ) :: &
        BS_DO_KERNEL_PARAMS_WFS

!  Derived quantity (tells you when to do BRDF derivatives)

      LOGICAL, dimension ( MAX_BRDF_KERNELS )  :: BS_DO_KPARAMS_DERIVS

!  Number of surface weighting functions

      INTEGER    :: BS_N_SURFACE_WFS
      INTEGER    :: BS_N_KERNEL_FACTOR_WFS
      INTEGER    :: BS_N_KERNEL_PARAMS_WFS

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!      Warning: BSA only works if NBEAMS = 1, since BSA is dependent on solar angle.

      LOGICAL    :: BS_DO_BSAVALUE_WF 
      LOGICAL    :: BS_DO_WSAVALUE_WF

!  NewCM option. Only linearization w.r.t Wind speed. Version 3.7

      LOGICAL    :: BS_DO_WINDSPEED_WF

      end type BRDF_LinSup_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: BRDF_LinSup_Inputs

      end module BRDF_LinSup_Inputs_def_m

