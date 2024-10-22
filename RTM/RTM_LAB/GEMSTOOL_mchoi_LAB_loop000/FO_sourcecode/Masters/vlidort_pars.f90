! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5                                       #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

!  This is VLIDORT.PARS.

!  File of constants for VLIDORT model.

      MODULE vlidort_pars

      IMPLICIT NONE

!  Real number type definitions

      INTEGER, PARAMETER :: VLIDORT_SPKIND = SELECTED_REAL_KIND(6)
      INTEGER, PARAMETER :: VLIDORT_DPKIND = SELECTED_REAL_KIND(15)
      INTEGER, PARAMETER :: FPK = VLIDORT_DPKIND

!  Version number
!  ==============

      CHARACTER (LEN=5), PARAMETER :: VLIDORT_VERSION_NUMBER = '2.5'

!  File i/o unit numbers
!  ======================


      INTEGER, PARAMETER ::       VLIDORT_INUNIT   = 21
      INTEGER, PARAMETER ::       VLIDORT_SCENUNIT = 22
      INTEGER, PARAMETER ::       VLIDORT_FUNIT    = 23
      INTEGER, PARAMETER ::       VLIDORT_RESUNIT  = 24
      INTEGER, PARAMETER ::       VLIDORT_ERRUNIT  = 25

!  Basic dimensions
!  ================

!  Computational dimensioning
!  --------------------------

!  Number of computational streams in the half-space

      INTEGER, PARAMETER :: MAXSTREAMS = 30

!  Maximum number of computational layers

      INTEGER, PARAMETER :: MAXLAYERS = 130

!  Maximum number of fine layers used in single scattering corrections

      INTEGER, PARAMETER :: MAXFINELAYERS = 1

!  Maximum number of input moments.
!    (Use full range for exact single scatter calculations)

      INTEGER, PARAMETER :: MAXMOMENTS_INPUT = 500

!  Max number of thermal coefficients
!   --------- New for Version 2.4RT -----------------

      INTEGER, PARAMETER :: MAX_THERMAL_COEFFS = 2

!  Geometrical and output parameters
!  ---------------------------------

!  Maximum number of solar zenith angles

      INTEGER, PARAMETER :: MAX_SZANGLES = 6 

!  maximum number of user-defined viewing zenith angles

      INTEGER, PARAMETER :: MAX_USER_VZANGLES = 6

!  maximum number of user-defined output relative azimuth angles

      INTEGER, PARAMETER :: MAX_USER_RELAZMS = 6

!  Maximum number of output optical depths

      INTEGER, PARAMETER :: MAX_USER_LEVELS = 2

!  Maximum number of output optical depths away from layer boundaries
!   This must be less than or equal to the previous entry

      INTEGER, PARAMETER :: MAX_PARTLAYERS = 1

!  Fixed parameters
!  ----------------

!  Number of Stokes parameters

      INTEGER, PARAMETER :: MAXSTOKES = 4

!  Two directions (Up and Down)

      INTEGER, PARAMETER :: MAX_DIRECTIONS = 2

!  Surface BRDF dimensioning
!  -------------------------

!  Maximum number of BRDF kernels

      INTEGER, PARAMETER :: MAX_BRDF_KERNELS = 3

!  Maximum number of BRDF parameters per kernel

      INTEGER, PARAMETER :: MAX_BRDF_PARAMETERS = 3

!  Maximum number of azimuth-quadrature streams for BGRDF Fourier.

      INTEGER, PARAMETER :: MAXSTREAMS_BRDF = 2
!      INTEGER, PARAMETER :: MAXSTREAMS_BRDF = 101    ! best

!  Weighting functions
!  -------------------

!  Maximum number of profile/column weighting functions

!      INTEGER, PARAMETER :: MAX_ATMOSWFS = 15
!      INTEGER, PARAMETER :: MAX_ATMOSWFS = 22
      INTEGER, PARAMETER :: MAX_ATMOSWFS = 1

!  Maximum number of surface property weighting functions

      INTEGER, PARAMETER :: MAX_SURFACEWFS = 1

!  Maximum number of error messages

      INTEGER, PARAMETER :: MAX_MESSAGES = 25

!  Derived dimensions
!  ==================

!  Copy Beam dimensioning

      INTEGER, PARAMETER :: MAXBEAMS = MAX_SZANGLES

!  Copy viewing zenith angle dimension

      INTEGER, PARAMETER :: MAX_USER_STREAMS = MAX_USER_VZANGLES

!  Maximum possible geometries

      INTEGER, PARAMETER :: MAX_GEOMETRIES = &
                            MAX_USER_VZANGLES*MAX_USER_RELAZMS*MAX_SZANGLES

!  All streams

      INTEGER, PARAMETER :: MAX_ALLSTRMS = MAX_USER_STREAMS + MAXSTREAMS

!  All streams for the Legendre PI-matrix setup.
!   Refractive Goemetry setting: Watch out for Kill. Memory Hog
!      INTEGER, PARAMETER :: MAX_ALLSTRMS_P1 = &
!                            MAX_ALLSTRMS + MAXBEAMS*MAXLAYERS

!  All streams for the Legendre PI-matrix setup.
!   Straightline setting: This setting should avoid dimensioning error

      INTEGER, PARAMETER :: MAX_ALLSTRMS_P1 = MAX_ALLSTRMS + MAXBEAMS

!  Maximum number of moments in the diffuse field calculation
!   This is always 2*MAXSTREAMS, in case we need DELTA-M

      INTEGER, PARAMETER :: MAXMOMENTS = 2*MAXSTREAMS

!  Maximum number of Fourier components = 2*MAXSTREAMS - 1

      INTEGER, PARAMETER :: MAXFOURIER = 2*MAXSTREAMS - 1

!  Number of Stokes streams squared

      INTEGER, PARAMETER :: MAXSTOKES_SQ = 16

!  Half the number of BRDF azimuth quadratures

      INTEGER, PARAMETER :: MAXSTHALF_BRDF = MAXSTREAMS_BRDF / 2

!  Other derived dimensions

      INTEGER, PARAMETER :: MAXSTREAMS_2   = 2*MAXSTREAMS
      INTEGER, PARAMETER :: MAXSTREAMS_21  = 2*MAXSTREAMS - 1
      INTEGER, PARAMETER :: MAXSTREAMS_P1  = MAXSTREAMS + 1
      INTEGER, PARAMETER :: MAXSTREAMS_P2  = MAXSTREAMS + 2

      INTEGER, PARAMETER :: MAXSTRMSTKS    = MAXSTREAMS * MAXSTOKES
      INTEGER, PARAMETER :: MAXSTRMSTKS_2  = 2*MAXSTRMSTKS
      INTEGER, PARAMETER :: MAXSTRMSTKS_21 = 2*MAXSTRMSTKS - 1

      INTEGER, PARAMETER :: MAXSTRMSTKS_P1 = MAXSTRMSTKS + 1
      INTEGER, PARAMETER :: MAXSTRMSTKS_P2 = MAXSTRMSTKS + 2
      INTEGER, PARAMETER :: MAXSTRMSTKS_P4 = MAXSTRMSTKS + 4

      INTEGER, PARAMETER :: MAX_USTRMSTKS = MAX_USER_STREAMS * MAXSTOKES

!  Maximum number of eigenvalues

      INTEGER, PARAMETER :: MAXEVALUES = MAXSTRMSTKS

!  For the BVP problem

      INTEGER, PARAMETER :: MAXTOTAL = MAXLAYERS*MAXSTRMSTKS_2
      INTEGER, PARAMETER :: MAXBANDTOTAL = 9*MAXSTRMSTKS - 2

!  Not so far used

      INTEGER, PARAMETER :: MAX_PSOLS = 2
      INTEGER, PARAMETER :: MAX_SCATPSOLS = MAX_PSOLS

!  Format constants
!  ================

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_HEADING = '( / T6, ''-----> '', A, /)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_INTEGER = '(T6, A, T58, I10)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_REAL    = '(T6, A, T58, 1PG14.6)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_CHAR    = '(T6, A, T48, A20)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_SECTION = '( / T6, ''-----> '', A, /)'
!     $                  '( // T6
!     $                  ''----------------------------------------'',
!     $                  ''-----------------------------------'',
!     $                      / T6 A,
!     $                      / T6
!     $                  ''----------------------------------------'',
!     $                  ''-----------------------------------'',
!     $                          / )' )

!  numbers
!  =======

      DOUBLE PRECISION, PARAMETER :: &
        ONE = 1.0D0, ZERO = 0.0D0,  ONEP5 = 1.5D0
      DOUBLE PRECISION, PARAMETER :: &
        TWO = 2.0D0, THREE = 3.0D0, FOUR = 4.0D0
      DOUBLE PRECISION, PARAMETER :: &
        QUARTER = 0.25D0, HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER :: &
        MINUS_ONE = -ONE
      DOUBLE PRECISION, PARAMETER :: &
        MINUS_TWO = -TWO
      DOUBLE PRECISION, PARAMETER :: &
        DEG_TO_RAD = 1.7453292519943D-02
      DOUBLE PRECISION, PARAMETER :: &
        PIE = 180.0D0*DEG_TO_RAD
      DOUBLE PRECISION, PARAMETER :: &
        PI2 = 2.0D0 * PIE
      DOUBLE PRECISION, PARAMETER :: &
        PI4 = 4.0D0 * PIE
      DOUBLE PRECISION, PARAMETER :: &
        PIO2 = HALF * PIE
      DOUBLE PRECISION, PARAMETER :: &
        PIO4 = QUARTER * PIE
      DOUBLE PRECISION, PARAMETER :: &
        EPS3 = 0.001D0
      DOUBLE PRECISION, PARAMETER :: &
        EPS4 = 0.0001D0
      DOUBLE PRECISION, PARAMETER :: &
        EPS5 = 0.00001D0
      DOUBLE PRECISION, PARAMETER :: &
        SMALLNUM = 1.0D-15
      DOUBLE PRECISION, PARAMETER :: &
        BIGEXP = 32.0D0

!  Control for Using L'Hopital's Rule
!   Changed, January 2009 for Version 2.4..........

      DOUBLE PRECISION, PARAMETER :: HOPITAL_TOLERANCE = EPS3
!      DOUBLE PRECISION, PARAMETER :: HOPITAL_TOLERANCE = EPS5

!  Control for limits of single scatter albedo

      DOUBLE PRECISION, PARAMETER :: OMEGA_SMALLNUM = 1.0D-15

!  Control for limits of extinction optical depth along solar path

      DOUBLE PRECISION, PARAMETER :: MAX_TAU_SPATH = 32.0D0

!  Control for limits of extinction optical depth along USER paths

      DOUBLE PRECISION, PARAMETER :: MAX_TAU_UPATH = 32.0D0

!  Control for limits of extinction optical depth along QUADRATURE paths

      DOUBLE PRECISION, PARAMETER :: MAX_TAU_QPATH = 32.0D0

!  error indices
!  =============

      INTEGER, PARAMETER :: VLIDORT_SERIOUS  = 4
      INTEGER, PARAMETER :: VLIDORT_WARNING  = 3
      INTEGER, PARAMETER :: VLIDORT_INFO     = 2
      INTEGER, PARAMETER :: VLIDORT_DEBUG    = 1
      INTEGER, PARAMETER :: VLIDORT_SUCCESS  = 0

!  directional indices

      INTEGER, PARAMETER :: UPIDX  = 1
      INTEGER, PARAMETER :: DNIDX  = 2

!  surface indices
!  ---------------

!  These refer to the BRDF kernel functions currently included.

      INTEGER, PARAMETER :: LAMBERTIAN_IDX       = 1
      INTEGER, PARAMETER :: ROSSTHIN_IDX         = 2
      INTEGER, PARAMETER :: ROSSTHICK_IDX        = 3
      INTEGER, PARAMETER :: LISPARSE_IDX         = 4
      INTEGER, PARAMETER :: LIDENSE_IDX          = 5
      INTEGER, PARAMETER :: HAPKE_IDX            = 6
      INTEGER, PARAMETER :: ROUJEAN_IDX          = 7
      INTEGER, PARAMETER :: RAHMAN_IDX           = 8
      INTEGER, PARAMETER :: COXMUNK_IDX          = 9
      INTEGER, PARAMETER :: GISSCOXMUNK_IDX      = 10
      INTEGER, PARAMETER :: GISSCOXMUNK_CRI_IDX  = 11

!  New for Version 2.4RTC

      INTEGER, PARAMETER :: BPDF2009_IDX         = 12

      INTEGER, PARAMETER :: MAXBRDF_IDX = BPDF2009_IDX

!  End of file.

      END MODULE vlidort_pars
