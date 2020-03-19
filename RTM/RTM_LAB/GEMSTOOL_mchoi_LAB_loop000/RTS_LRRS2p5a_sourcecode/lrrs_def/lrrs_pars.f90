
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

!  ======================================================
!  This is the module of TYPES, DIMENSIONS, and CONSTANTS
!  ======================================================

      MODULE lrrs_pars_m

      IMPLICIT NONE

!  Real number type definitions

      INTEGER, PARAMETER :: LRRS_SPKIND = SELECTED_REAL_KIND(6)
      INTEGER, PARAMETER :: LRRS_DPKIND = SELECTED_REAL_KIND(15)
      INTEGER, PARAMETER :: FPK = LRRS_DPKIND

!  DIMENSIONS
!  ==========

!  Basic: Physics Dimensioning
!  ---------------------------

!  Maximum number of spectral points
!  MONO: Should always be set to (at least) number of shifts + 1
!  BIN:  Can be less

!      INTEGER, PARAMETER :: MAX_POINTS = 234
!     INTEGER, PARAMETER :: MAX_POINTS = 100
      INTEGER, PARAMETER :: MAX_POINTS = 1500
!      INTEGER, PARAMETER :: MAX_POINTS = 600

!  Maximum number of layers

      INTEGER, PARAMETER :: MAX_LAYERS = 36

!  Maximum number of fine layers

      INTEGER, PARAMETER :: MAX_FINE_LAYERS = 4 !2

!  Maximum number of MS phase function moments
!    Can set this to twice the number of streams

      INTEGER, PARAMETER :: MAX_MOMENTS = 20
!      INTEGER, PARAMETER :: MAX_MOMENTS = 432

!  Maximum number of INPUT phase function moments
!     Should be at least 2 * MAXSTREAMS
!      INTEGER, PARAMETER :: MAX_MOMENTS_INPUT = 60
!      INTEGER, PARAMETER :: MAX_MOMENTS_INPUT = 301

!  Maximum number of RRS bins and shifts

!      INTEGER, PARAMETER :: MAX_BINS = 80
!      INTEGER, PARAMETER :: MAX_BINS = 50
!      INTEGER, PARAMETER :: MAX_BINS = 234
      INTEGER, PARAMETER :: MAX_BINS = 1500
      INTEGER, PARAMETER :: MAX_SHIFTS = 233

!  Basic: RT Dimensioning
!  ----------------------

!  Maximum number of discrete ordinates

      INTEGER, PARAMETER :: MAX_STREAMS = 10

!  Maximum numbers of off-boundary and total output levels

      INTEGER, PARAMETER :: MAX_PARTIALS_LOUTPUT = 1
      INTEGER, PARAMETER :: MAX_LOUTPUT = 5

!  Maximum numbers of user defined zeniths and azimuths

      INTEGER, PARAMETER :: MAX_USER_STREAMS = 3
      INTEGER, PARAMETER :: MAX_USER_RELAZMS = 3

!  Upwelling and downwelling

      INTEGER, PARAMETER :: MAX_DIRECTIONS = 2

!  Maximum number of Terms for Taylor series expansions. Introduced, Version 2.5 9/9/15
!    If you are retaining contributions of order EPS^n, then you need at least n+2 Taylor terms

      INTEGER, PARAMETER :: MAX_TAYLOR_TERMS = 7

!  Surface BRDF dimensioning
!  -------------------------

!  The next 5 definitions are new for version 2.5 (Rob Fix 9/8/15)
!    -- Taken straight from the LIDORT 3.7 code
!    -- Only required for the BRDF supplement

!  Maximum number of BRDF kernels

      INTEGER, PARAMETER :: MAX_BRDF_KERNELS = 3

!  Maximum number of BRDF parameters per kernel

      INTEGER, PARAMETER :: MAX_BRDF_PARAMETERS = 3

!  Maximum number of azimuth-quadrature streams for BRDF Fourier.

!      INTEGER, PARAMETER :: MAXSTREAMS_BRDF = 2
      INTEGER, PARAMETER :: MAXSTREAMS_BRDF = 101   ! best

!  Maximum numbers for the MSR quadratures

      INTEGER, PARAMETER :: MAX_MSRS_MUQUAD  = 50
      INTEGER, PARAMETER :: MAX_MSRS_PHIQUAD = 100

!  Number of quadrature streams for internal WSA/BSA scaling
!    New, Version 3.7. User does not need to know this value.

      INTEGER, parameter :: MAXSTREAMS_SCALING = 24

!  Maximum Number of points for the BRDF calculation
!     = 1, if single-value option is to be used 
!     = MAX_POINTS for multiple points (Bin/Mono realizations)

      INTEGER, parameter :: MAX_BRDF_POINTS = 1
!      INTEGER, parameter :: MAX_BRDF_POINTS = MAX_POINTS

!  Maximum Number of points for the SLEAVE calculation
!     = 1, if single-value option is to be used 
!     = MAX_POINTS for multiple points (Bin/Mono realizations)

      INTEGER, parameter :: MAX_SLEAVE_POINTS = 1
!      INTEGER, parameter :: MAX_SLEAVE_POINTS = MAX_POINTS

!  Weighting functions
!  -------------------

!  Maximum number of profile/column weighting functions

   INTEGER, PARAMETER :: MAX_ATMOSWFS = 4

!  Maximum number of surface property weighting functions

   INTEGER, PARAMETER :: MAX_SURFACEWFS = 7

!  Maximum number of surface-leaving weighting functions

   INTEGER, PARAMETER :: MAX_SLEAVEWFS = 2

!  Messages dimensioning
!  ---------------------

      INTEGER, PARAMETER :: MAX_MESSAGES = 25

!  Derived Dimensioning
!  --------------------

!  NK storage (for linearization arrays)

      INTEGER, PARAMETER :: MAX_LAYERS_NK = &
                            MAX_LAYERS * ( MAX_LAYERS + 3 ) / 2
      INTEGER, PARAMETER :: MAX_LAYERS_SQ = MAX_LAYERS * MAX_LAYERS

!  Max Geometries

      INTEGER, PARAMETER :: MAX_GEOMETRIES = &
                            MAX_USER_STREAMS * MAX_USER_RELAZMS

!  Half the number of BRDF azimuth quadratures. Introduced 2.5, 9/8/15

      INTEGER, PARAMETER :: MAXSTHALF_BRDF = MAXSTREAMS_BRDF / 2

!  Stream numbers

      INTEGER, PARAMETER :: MAX_2_STREAMS   = 2*MAX_STREAMS
      INTEGER, PARAMETER :: MAX_OUT_STREAMS = MAX_USER_STREAMS+MAX_STREAMS
      INTEGER, PARAMETER :: MAX_STREAMS_P1  = MAX_STREAMS + 1

!  Maximum Number of additional sourceterms in the RTS

      INTEGER, PARAMETER :: MAX_EXPONENTS = &
                            ( MAX_2_STREAMS + 1 ) * MAX_BINS + 1
      INTEGER, PARAMETER :: MAX_SOLUTIONS = &
                            ( MAX_2_STREAMS + 2 ) * MAX_BINS + 1

!  Dimensioning for the boundary value problem

      INTEGER, PARAMETER :: MAX_TOTAL     = MAX_LAYERS * MAX_2_STREAMS
      INTEGER, PARAMETER :: MAX_BANDTOTAL = 9 * MAX_STREAMS - 2

!  CONSTANTS
!  =========

!  Version numbers

      CHARACTER (LEN=8), PARAMETER :: LRRS_VERSION_NUMBER = 'LRRS_2.5'

!  File i/o unit numbers
!  ---------------------

      INTEGER, PARAMETER :: LRRS_INUNIT   = 21
      INTEGER, PARAMETER :: LRRS_SCENUNIT = 21
      INTEGER, PARAMETER :: LRRS_FUNIT    = 23
      INTEGER, PARAMETER :: LRRS_RESUNIT  = 24
      INTEGER, PARAMETER :: LRRS_ERRUNIT  = 25
      INTEGER, PARAMETER :: LRRS_DBGUNIT  = 71

!  Special file debug units

      INTEGER :: SDU,LDU

!  Format constants (Internal debug only)
!  ================

!  This section upgraded to conform with LIDORT Version 3.7, 9/8/15

      CHARACTER (LEN=*), PARAMETER :: FMT_HEADING = '( / T6, ''-----> '', A, /)'
      CHARACTER (LEN=*), PARAMETER :: FMT_INTEGER = '(T6, A, T58, I10)'
      CHARACTER (LEN=*), PARAMETER :: FMT_REAL    = '(T6, A, T58, 1PG14.6)'
      CHARACTER (LEN=*), PARAMETER :: FMT_CHAR    = '(T6, A, T48, A20)'
      CHARACTER (LEN=*), PARAMETER :: FMT_SECTION = '( / T6, ''-----> '', A, /)'

!  Debug write format (DWF) constants
!  ==================================

!  This section upgraded to conform with LIDORT Version 3.7, 9/8/15

      CHARACTER (LEN=*), PARAMETER :: DWFL  = '(A,L1)'
      CHARACTER (LEN=*), PARAMETER :: DWFL1 = '(A,I3,A,L1)'
      CHARACTER (LEN=*), PARAMETER :: DWFL2 = '(2(A,I3),A,L1)'

      CHARACTER (LEN=*), PARAMETER :: DWFI  = '(A,I5)'
      CHARACTER (LEN=*), PARAMETER :: DWFI1 = '(A,I3,A,I5)'
      CHARACTER (LEN=*), PARAMETER :: DWFI2 = '(2(A,I3),A,I5)'

      CHARACTER (LEN=*), PARAMETER :: DWFR  = '(A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR1 = '(A,I3,A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR2 = '(2(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR3 = '(3(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR4 = '(4(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR5 = '(5(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR6 = '(6(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR7 = '(7(A,I3),A,ES13.6E2)'

      CHARACTER (LEN=*), PARAMETER :: DWFR1_3 = '(A,I3,3(A,ES13.6E2))'

      CHARACTER (LEN=*), PARAMETER :: DWFC  = '(2A)'
      CHARACTER (LEN=*), PARAMETER :: DWFC1 = '(A,I3,2A)'
      CHARACTER (LEN=*), PARAMETER :: DWFC2 = '(2(A,I3),2A)'

!  Numbers
!  -------

!  This section upgraded to conform with LIDORT Version 3.7, 9/8/15

      real(fpk), PARAMETER :: ONE = 1.0_fpk, ZERO  = 0.0_fpk, ONEP5 = 1.5_fpk
      real(fpk), PARAMETER :: TWO = 2.0_fpk, THREE = 3.0_fpk, FOUR  = 4.0_fpk
      real(fpk), PARAMETER :: QUARTER = 0.25_fpk, HALF = 0.5_fpk

      real(fpk), PARAMETER :: MINUS_ONE = - ONE
      real(fpk), PARAMETER :: MINUS_TWO = - TWO

      real(fpk), PARAMETER :: PIE = ACOS(MINUS_ONE)
      real(fpk), PARAMETER :: DEG_TO_RAD = PIE/180.0_fpk

      real(fpk), PARAMETER :: PI2  = TWO  * PIE
      real(fpk), PARAMETER :: PI4  = FOUR * PIE
      real(fpk), PARAMETER :: PIO2 = HALF * PIE
      real(fpk), PARAMETER :: PIO4 = QUARTER * PIE

      real(fpk), PARAMETER :: EPS3 = 0.001_fpk
      real(fpk), PARAMETER :: EPS4 = 0.0001_fpk
      real(fpk), PARAMETER :: EPS5 = 0.00001_fpk

!  Rob fix 5/6/13 - Taylor series limiting values

      !real(fpk), PARAMETER :: TAYLOR_SMALL = 0.001_fpk

   !real(fpk), PARAMETER :: TAYLOR_SMALL = 0.0005_fpk
   !real(fpk), PARAMETER :: TAYLOR_SMALL = 0.0010_fpk
   !real(fpk), PARAMETER :: TAYLOR_LARGE = 10000.0_fpk (not used now)

      real(fpk), PARAMETER :: SMALLNUM = 0.000000001_fpk

!  Maximum negative exponential argument

      real(fpk), PARAMETER :: BIGEXP = 32.0_fpk

!  Control for Using L'Hopital's Rule

      real(fpk), PARAMETER :: HOPITAL_TOLERANCE = EPS5

!  Control for limits of single scatter albedo

      real(fpk), PARAMETER :: OMEGA_SMALLNUM = 0.00000001_fpk

!  Control for limits of extinction optical depth along solar path
!  Control for limits of extinction optical depth along USER paths
!  Control for limits of extinction optical depth along QUADRATURE paths

      real(fpk), PARAMETER :: MAX_TAU_SPATH = BIGEXP
      real(fpk), PARAMETER :: MAX_TAU_UPATH = BIGEXP
      real(fpk), PARAMETER :: MAX_TAU_QPATH = BIGEXP

!  Error indices
!  -------------

      INTEGER, PARAMETER :: LRRS_SERIOUS = 2
      INTEGER, PARAMETER :: LRRS_WARNING = 1
      INTEGER, PARAMETER :: LRRS_SUCCESS = 0

!  Directional indices
!  -------------------

      INTEGER, PARAMETER :: UPIDX = 1
      INTEGER, PARAMETER :: DNIDX = 2

!  Surface Type indices
!  --------------------

!  Upgraded in conformance with the latest LIDORT and VLIDORT indices

!  These refer to the BRDF kernel functions currently included.

   INTEGER, PARAMETER :: LAMBERTIAN_IDX  = 1
   INTEGER, PARAMETER :: ROSSTHIN_IDX    = 2
   INTEGER, PARAMETER :: ROSSTHICK_IDX   = 3
   INTEGER, PARAMETER :: LISPARSE_IDX    = 4
   INTEGER, PARAMETER :: LIDENSE_IDX     = 5
   INTEGER, PARAMETER :: HAPKE_IDX       = 6
   INTEGER, PARAMETER :: ROUJEAN_IDX     = 7
   INTEGER, PARAMETER :: RAHMAN_IDX      = 8
   INTEGER, PARAMETER :: COXMUNK_IDX     = 9

!  New BPDF functions, matching those in LIDORT 3.7 and VLIDORT 2.7

   INTEGER, PARAMETER :: BPDFSOIL_IDX         = 10
   INTEGER, PARAMETER :: BPDFVEGN_IDX         = 11
   INTEGER, PARAMETER :: BPDFNDVI_IDX         = 12

!  New Cox-Munk function for LIDORT 3.7, VLIDORT 2.7

   INTEGER, PARAMETER :: NewCMGLINT_IDX       = 13
   INTEGER, PARAMETER :: MAXBRDF_IDX = NewCMGLINT_IDX

!  End of Module.

      END MODULE lrrs_pars_m

