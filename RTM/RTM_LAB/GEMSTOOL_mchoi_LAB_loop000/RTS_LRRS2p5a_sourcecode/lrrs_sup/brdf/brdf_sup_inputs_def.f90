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

      module BRDF_Sup_Inputs_def_m

! #####################################################################
! #####################################################################

!  This module contains the following structures:

!  BRDF_Sup_Inputs - Intent(In) for BRDF_Sup

      use LRRS_PARS_m, only : fpk, MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_BRDF_POINTS, &
                              MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS

      implicit none

! #####################################################################
! #####################################################################

      type BRDF_Sup_Inputs

!  Stream angle flag

      LOGICAL :: BS_DO_USER_STREAMS

!  BRDF surface flag

      LOGICAL :: BS_DO_BRDF_SURFACE

!  Number of discrete ordinate streams

      INTEGER :: BS_NSTREAMS

!  Solar angle

      REAL(fpk) :: BS_BEAM_SZA

!  user-defined relative azimuths

      INTEGER                                 :: BS_N_USER_RELAZMS
      REAL(fpk), dimension (MAX_USER_RELAZMS) :: BS_USER_RELAZMS

!  User-defined zenith angle input

      INTEGER                                 :: BS_N_USER_STREAMS
      REAL(fpk), dimension (MAX_USER_STREAMS) :: BS_USER_ANGLES_INPUT

!  BRDF-specific inputs
!  --------------------

!   Number and index-list of bidirectional functions

      INTEGER                                            :: BS_N_BRDF_KERNELS
      CHARACTER (LEN=10), dimension ( MAX_BRDF_KERNELS ) :: BS_BRDF_NAMES
      INTEGER, dimension ( MAX_BRDF_KERNELS )            :: BS_WHICH_BRDF

!  Parameters required for Kernel families

      INTEGER  , dimension ( MAX_BRDF_KERNELS )                      :: BS_N_BRDF_PARAMETERS
      REAL(fpk), dimension ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS ) :: BS_BRDF_PARAMETERS

!  Lambertian Surface control

      LOGICAL, dimension ( MAX_BRDF_KERNELS )   :: BS_LAMBERTIAN_KERNEL_FLAG

!  Input kernel amplitude factors

      REAL(fpk), dimension ( MAX_BRDF_KERNELS ) :: BS_BRDF_FACTORS

!  Number of azimuth quadrature streams for BRDF

      INTEGER :: BS_NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL :: BS_DO_SHADOW_EFFECT

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!  Rob Fix 9/25/14. Two variables replaced

!      LOGICAL :: BS_DO_EXACT
!      LOGICAL :: BS_DO_EXACTONLY
      LOGICAL :: BS_DO_DIRECTBOUNCE_ONLY
!      LOGICAL :: BS_DO_DIRECTBOUNCE_TRUNCATED

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 02 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised, 12 August 2014, added option to output calculated WSA/BSA

      LOGICAL   :: BS_DO_WSABSA_OUTPUT
      LOGICAL   :: BS_DO_WSA_SCALING
      LOGICAL   :: BS_DO_BSA_SCALING
      REAL(fpk) :: BS_WSA_VALUE, BS_BSA_VALUE

!  New Cox-Munk Glint options (bypasses the usual Kernel system)
!  -------------------------------------------------------------

!  Overall flag for this option

      LOGICAL   :: BS_DO_NewCMGLINT

!  Input Salinity in [ppt]

      REAL(fpk) :: BS_SALINITY

!  Flag for using a single wavelength. Introduced 9/8/15 for LRRS

      LOGICAL   :: BS_DO_NewCM_Wav1

!  Number of wavelengths used. Introduced 9/8/15 for LRRS
!     = 1, if single-value option is set (dimensioning also 1 - check this!)
!     = NPOINTS_INNER for multiple points (Bin realization)
!     = NPOINTS_MONO  for multiple points (Mono realization)

      INTEGER   :: BS_NewCM_Npoints

!  Input wavelengths in [nm]. Array introduced 9/8/15 for LRRS.
!    For single-point (Bin realization) , set to average-value of "Lamdas_Ranked" 
!    For single-point (Mono realization), set to Excitation wavelength "Lamdas_Ranked(W_EXCIT)" 
!    For all points, copy LAMDAS_RANKED array (both realizations)
!    - Note units are same as LAMBDAS_RANKED which is in [nm].

      REAL(fpk) :: BS_NewCM_Lambdas ( MAX_BRDF_POINTS )

!  Input Wind speed in m/s, and azimuth directions relative to Sun position

      REAL(fpk) :: BS_WINDSPEED, BS_WINDDIR

!  Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: BS_DO_GlintShadow
      LOGICAL   :: BS_DO_FoamOption
      LOGICAL   :: BS_DO_FacetIsotropy

!  Multiple-scattering Glitter options
!  -----------------------------------

!  Multiple reflectance corrections for GLITTER kernels (All of them!)

      LOGICAL :: BS_DO_GLITTER_MSRCORR

!  Multiple reflectance correction for exact-term Glitter kernels only
!  mick fix 12/29/2014 - Name changed from EXACTONLY --> DBONLY

      LOGICAL :: BS_DO_GLITTER_MSRCORR_DBONLY

!  Correction order for the Multiple reflectance computations
!    ( = 0, no correction, 1, 2, 3   ec.)
!  Warning; using S > 0 can increase CPU dramatically

      INTEGER :: BS_GLITTER_MSRCORR_ORDER

!  Quadrature orders for MSRCORR

      INTEGER :: BS_GLITTER_MSRCORR_NMUQUAD
      INTEGER :: BS_GLITTER_MSRCORR_NPHIQUAD

      end type BRDF_Sup_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: BRDF_Sup_Inputs

      end module BRDF_Sup_Inputs_def_m

