
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

! ###############################################################
! #                                                             #
! #   SUBROUTINE :                                              #
! #              RAMAN_CONVERGE  (*)                            #
! #                                                             #
! #      (*) Version 2.1 contains Exact DB correction (BRDFs)   #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Bookkeeping improvements (use of "Only", clearer I/O specifications)

!   -- Rob mod 5/12/17 for 2p5a, DO_SSFULL renamed --> DO_SSCORR_ALONE
!   -- Rob mod 5/12/17 for 2p5a, Use DO_SSCORR_GENERAL (OUTGOING or NADIR)

      MODULE lrrs_converge_m

!      USE LRRS_PARS_m, Only : SDU

      PRIVATE
      PUBLIC :: RAMAN_CONVERGE

      CONTAINS

      SUBROUTINE RAMAN_CONVERGE &
       ( DO_UPWELLING, DO_DNWELLING, DO_DOUBLE_CONVTEST,            & ! Inputs
         DO_SSCORR_GENERAL, DO_SSCORR_ALONE, DO_RAYLEIGH_ONLY,      & ! Inputs
         DO_NO_AZIMUTH, DO_FILLING, DO_ELASTIC_ONLY, FOURIER,       & ! Inputs
         NPOINTS_LOCAL, N_LOUTPUT, N_OUT_STREAMS, LOCAL_N_USERAZM,  & ! Inputs
         N_CONVTESTS, N_GEOMETRIES, GEOM_OFFSETS, AZMFAC, ACCURACY, & ! Inputs
         ELASTIC_SS_UP, ELASTIC_SS_DN, RAMAN_SS_UP, RAMAN_SS_DN,    & ! Inputs
         ELASTIC_F_UP, ELASTIC_F_DN, RAMAN_F_UP, RAMAN_F_DN,        & ! Inputs
         ELASTIC_UP, ELASTIC_DN, RAMAN_UP, RAMAN_DN,                & ! Outputs
         FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )                   ! Outputs

!  include file of dimensions and numbers 

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LOUTPUT, MAX_USER_STREAMS, &
                              MAX_USER_RELAZMS, MAX_GEOMETRIES, MAX_POINTS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Upwelling and downwelling flags

      LOGICAL  , INTENT(IN) :: DO_UPWELLING
      LOGICAL  , INTENT(IN) :: DO_DNWELLING

!  Convergence test control

      LOGICAL  , INTENT(IN) :: DO_DOUBLE_CONVTEST

!  Do Rayleigh only flag - should really be "Do Molecular Only"

      LOGICAL  , INTENT(IN) :: DO_RAYLEIGH_ONLY

!  Azimuth only control

      LOGICAL  , INTENT(IN) :: DO_NO_AZIMUTH

!  SS correction flags; Now incorporates the DB term
!   -- Rob mod 5/12/17 for 2p5a, Use DO_SSCORR_GENERAL (nadir or outgoing), rename SSFULL

      LOGICAL  , INTENT(IN) :: DO_SSCORR_GENERAL
      LOGICAL  , INTENT(IN) :: DO_SSCORR_ALONE

!  Filling flag..................
!  If set,  Fourier sums for elastic AND inelastic solutions are done

      LOGICAL  , INTENT(IN) :: DO_FILLING

!  Elastic only flag

      LOGICAL  , INTENT(IN) :: DO_ELASTIC_ONLY

!  Number of points

      INTEGER  , INTENT(IN) :: NPOINTS_LOCAL

!  Number of output streams

      INTEGER  , INTENT(IN) :: N_OUT_STREAMS

!  Number of output levels

      INTEGER  , INTENT(IN) :: N_LOUTPUT

!  Local number of azimuths

      INTEGER  , INTENT(IN) :: LOCAL_N_USERAZM

!  Fourier index

      INTEGER  , INTENT(IN) :: FOURIER

!  Number of convergence tests

      INTEGER  , INTENT(IN) :: N_CONVTESTS

!  Number of Geometries, Geometry offsetting
!    --Added by R. Spurr, 04 August 2009

      INTEGER  , INTENT(IN) :: N_GEOMETRIES
      INTEGER  , INTENT(IN) :: GEOM_OFFSETS ( MAX_USER_STREAMS )

!  Azimuth angles. Adjusted factors, 18 March 2011.

      REAL(FPK), INTENT(IN) :: AZMFAC( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Overall accuracy for convergence criterion

      REAL(FPK), INTENT(IN) :: ACCURACY

!  Input single scatter fields
!  ---------------------------

      REAL(FPK), INTENT(IN) :: ELASTIC_SS_UP &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(IN) :: ELASTIC_SS_DN &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(IN) :: RAMAN_SS_UP &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(IN) :: RAMAN_SS_DN &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Input Fourier component fields
!  ------------------------------

      REAL(FPK), INTENT(IN) :: ELASTIC_F_UP &
              ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: ELASTIC_F_DN &
              ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: RAMAN_F_UP &
              ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: RAMAN_F_DN &
              ( MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  Modified/output variables
!  -------------------------

!  Fourier summed output: Radiance purely elastic scattering

      REAL(FPK), INTENT(INOUT) :: ELASTIC_UP &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: ELASTIC_DN &
       (MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Fourier summed output: Radiance including inelastic scattering

      REAL(FPK), INTENT(INOUT) :: RAMAN_UP &
             ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: RAMAN_DN &
             ( MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS)

!  convergence counter and flag

      LOGICAL  , INTENT(INOUT) :: LOCAL_ITERATION
      INTEGER  , INTENT(INOUT) :: TESTCONV

!  Record of the number of Fourier terms used

      INTEGER  , INTENT(INOUT) :: FOURIER_SAVED

!  local variables
!  ---------------

      INTEGER   :: COUNT(MAX_POINTS), COUNT_TOTAL
      INTEGER   :: S, I, UT, UA, V
      REAL(FPK) :: TNEW, TOLD, TAZM

!  Single scattering ONLY (SSCORR_ALONE is flagged). Copy and return
!  -----------------------------------------------------------------

      IF ( DO_SSCORR_GENERAL .AND. DO_SSCORR_ALONE ) THEN
        IF ( DO_UPWELLING ) THEN
          DO UT = 1, N_LOUTPUT
            DO V = 1, N_GEOMETRIES
              DO S = 1, NPOINTS_LOCAL
                ELASTIC_UP(UT,V,S) = ELASTIC_SS_UP(UT,V,S)
              ENDDO
              IF ( .NOT. DO_ELASTIC_ONLY ) THEN
                DO S = 1, NPOINTS_LOCAL
                  RAMAN_UP(UT,V,S)   = RAMAN_SS_UP(UT,V,S)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_DNWELLING ) THEN
          DO UT = 1, N_LOUTPUT
            DO V = 1, N_GEOMETRIES
              DO S = 1, NPOINTS_LOCAL
                ELASTIC_DN(UT,V,S) = ELASTIC_SS_DN(UT,V,S)
              ENDDO
              IF ( .NOT. DO_ELASTIC_ONLY ) THEN
                DO S = 1, NPOINTS_LOCAL
                  RAMAN_DN(UT,V,S)   = RAMAN_SS_DN(UT,V,S)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        RETURN
      ENDIF

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Diffuse term addition
!  ---------------------

!  Skip if elastic only

        IF ( DO_ELASTIC_ONLY ) GO TO 677

!  Copy Fourier component at all output angles and layer output choices

        IF ( DO_UPWELLING ) THEN
          DO UT = 1, N_LOUTPUT
            DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = GEOM_OFFSETS(I) + UA
               DO S = 1, NPOINTS_LOCAL
                RAMAN_UP(UT,V,S) = RAMAN_F_UP(UT,I,S)
               ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        IF ( DO_DNWELLING ) THEN
          DO UT = 1, N_LOUTPUT
            DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = GEOM_OFFSETS(I) + UA
               DO S = 1, NPOINTS_LOCAL
                RAMAN_DN(UT,V,S) = RAMAN_F_DN(UT,I,S)
               ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Continuation point for Elastic-only

 677    CONTINUE

        IF ( DO_FILLING .OR. DO_ELASTIC_ONLY ) THEN
         IF ( DO_UPWELLING ) THEN
          DO UT = 1, N_LOUTPUT
            DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = GEOM_OFFSETS(I) + UA
               DO S = 1, NPOINTS_LOCAL
                ELASTIC_UP(UT,V,S) = ELASTIC_F_UP(UT,I,S)
               ENDDO
              ENDDO
            ENDDO
          ENDDO
         ENDIF
         IF ( DO_DNWELLING ) THEN
          DO UT = 1, N_LOUTPUT
            DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = GEOM_OFFSETS(I) + UA
               DO S = 1, NPOINTS_LOCAL
                ELASTIC_DN(UT,V,S) = ELASTIC_F_DN(UT,I,S)
               ENDDO
              ENDDO
            ENDDO
          ENDDO
         ENDIF
        ENDIF

!  Add the single scatter component if flagged
!  -------------------------------------------

!  Skip if elastic only

        IF ( DO_ELASTIC_ONLY ) GO TO 678

!  OLD STYLE NADIR CORRECTION IS REMOVED.......

!    If set,convergence RADIANCE = INELASTIC_DIFFUSE + ELASTIC_SSEXACT
!        IF ( DO_SSCORR_NADIR ) THEN
!         IF ( DO_UPWELLING ) THEN
!           DO UT = 1, N_LOUTPUT
!             DO V = 1, N_GEOMETRIES
!               DO S = 1, NPOINTS_LOCAL
!                 RAMAN_UP(UT,V,S) =
!     &             RAMAN_UP(UT,V,S) + ELASTIC_SS_UP(UT,V,S)
!               ENDDO
!             ENDDO
!           ENDDO
!         ENDIF
!         IF ( DO_DNWELLING ) THEN
!           DO UT = 1, N_LOUTPUT
!             DO V = 1, N_GEOMETRIES
!               DO S = 1, NPOINTS_LOCAL
!                 RAMAN_DN(UT,V,S) =
!     &           RAMAN_DN(UT,V,S) + ELASTIC_SS_DN(UT,V,S)
!               ENDDO
!             ENDDO
!           ENDDO
!          ENDIF
!        ENDIF

!  GENERAL (new-style) correction:
!  --------------------------------

!    If set,convergence RADIANCE = INELASTIC_DIFFUSE + INELASTIC_SSEXACT

        IF ( DO_SSCORR_GENERAL ) THEN
          IF ( DO_UPWELLING ) THEN
            DO UT = 1, N_LOUTPUT
              DO V = 1, N_GEOMETRIES
                DO S = 1, NPOINTS_LOCAL
                  RAMAN_UP(UT,V,S) = &
                     RAMAN_UP(UT,V,S) + RAMAN_SS_UP(UT,V,S)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            DO UT = 1, N_LOUTPUT
              DO V = 1, N_GEOMETRIES
                DO S = 1, NPOINTS_LOCAL
                  RAMAN_DN(UT,V,S) = &
                     RAMAN_DN(UT,V,S) + RAMAN_SS_DN(UT,V,S)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  Continuation point for Elastic-only

 678    continue

!    If set,convergence RADIANCE = ELASTIC_DIFFUSE + ELASTIC_SSEXACT
!     Only if filling flag is set, or elastic-only

!        IF ( DO_SSCORR_GENERAL .or. DO_SSCORR_NADIR ) THEN

        IF (  DO_SSCORR_GENERAL ) THEN
          IF ( DO_FILLING .OR. DO_ELASTIC_ONLY ) THEN
            IF ( DO_UPWELLING ) THEN
              DO UT = 1, N_LOUTPUT
                DO V = 1, N_GEOMETRIES
                  DO S = 1, NPOINTS_LOCAL
                    ELASTIC_UP(UT,V,S) = &
                      ELASTIC_UP(UT,V,S) + ELASTIC_SS_UP(UT,V,S)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
            IF ( DO_DNWELLING ) THEN
              DO UT = 1, N_LOUTPUT
                DO V = 1, N_GEOMETRIES
                  DO S = 1, NPOINTS_LOCAL
                    ELASTIC_DN(UT,V,S) = &
                      ELASTIC_DN(UT,V,S) + ELASTIC_SS_DN(UT,V,S)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF

! *********************************************************************
!           DB TERM NOW INCORPORATED in THE SS TERM UPWELLING
! *********************************************************************
!  Direct Beam correction, only with Upwelling BRDF
!   New code, 28 November 2007. Add direct beam component if flagged
!        IF ( DO_DB_CORRECTION.AND.DO_UPWELLING ) THEN
!    If set, convergence INELASTIC_RADIANCE =
!       INELASTIC_DIFFUSE + ELASTIC_SSEXACT + ELASTIC_DBEXACT
!          DO UT = 1, N_LOUTPUT
!            DO V = 1, N_GEOMETRIES
!              DO S = 1, NPOINTS_LOCAL
!                RAMAN_UP(UT,V,S) =
!     &           RAMAN_UP(UT,V,S) + ELASTIC_DB(UT,V,S)
!              ENDDO
!            ENDDO
!          ENDDO
!     Elastic Field, Only if filling flag is set
!    If set, convergence ELASTIC_RADIANCE =
!       ELASTIC_DIFFUSE + ELASTIC_SSEXACT + ELASTIC_DBEXACT
!         IF ( DO_FILLING ) THEN
!           DO UT = 1, N_LOUTPUT
!             DO V = 1, N_GEOMETRIES
!               DO S = 1, NPOINTS_LOCAL
!                 ELASTIC_UP(UT,V,S) =
!     &            ELASTIC_UP(UT,V,S) + ELASTIC_DB(UT,V,S)
!               ENDDO
!             ENDDO
!           ENDDO
!         ENDIF
!  End direct beam correction
!        ENDIF
! *********************************************************************
! *********************************************************************

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  For Rayleigh atmosphere
!  -----------------------

!   skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY ) THEN

!  Skip if elastic only

          IF ( DO_ELASTIC_ONLY ) GO TO 679

!  For each azimuth, add Fourier component
!     - for direction, user output choice, out stream

          IF ( DO_UPWELLING ) THEN
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = RAMAN_UP(UT,V,S)
                  TAZM = AZMFAC(I,UA)*RAMAN_F_UP(UT,I,S)
                  RAMAN_UP(UT,V,S) = TOLD + TAZM
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = RAMAN_DN(UT,V,S)
                  TAZM = AZMFAC(I,UA)*RAMAN_F_DN(UT,I,S)
                  RAMAN_DN(UT,V,S) = TOLD + TAZM
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  Continuation point for Elastic-only

 679      continue

!  Repeat for the elastic only calculation

          IF ( DO_FILLING .OR. DO_ELASTIC_ONLY ) THEN
           IF ( DO_UPWELLING ) THEN
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = ELASTIC_UP(UT,V,S)
                  TAZM = AZMFAC(I,UA)*ELASTIC_F_UP(UT,I,S)
                  ELASTIC_UP(UT,V,S) = TOLD + TAZM
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDIF
           IF ( DO_DNWELLING ) THEN
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = ELASTIC_DN(UT,V,S)
                  TAZM = AZMFAC(I,UA)*ELASTIC_F_DN(UT,I,S)
                  ELASTIC_DN(UT,V,S) = TOLD + TAZM
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDIF
          ENDIF

!  Examine convergence on INTENSITY only
!  -------------------------------------

!  DO NOT EXAMINE CONVERGENCE ON ELASTIC radiances
!    ---------- Except when elastic-only flag is set !!!!!

!  convergence test applied to ALL directions AND
!                              ALL stream values (except near zenith) AN
!                              ALL azimuths taken together
!                              ALL user optical depths
!                      (New!)  ALL wavealength points

        ELSE

!  Initialise count

          DO S = 1, NPOINTS_LOCAL
            COUNT(S) = 0
          ENDDO

!  Skip if elastic-only flag is set

          IF ( DO_ELASTIC_ONLY ) GO TO 680

!  Count number of occasions Fourier term addition is below accuracy lev

          DO UA = 1, LOCAL_N_USERAZM
            IF ( DO_UPWELLING ) THEN
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = RAMAN_UP(UT,V,S)
                  TAZM = AZMFAC(I,UA)*RAMAN_F_UP(UT,I,S)
                  TNEW = TOLD + TAZM
                  IF ( TAZM .NE. ZERO ) THEN
                    IF ( ABS(TAZM/TNEW) .LT. ACCURACY ) THEN
                      COUNT(S) = COUNT(S) + 1
                    ENDIF
                  ELSE
                    COUNT(S)  = COUNT(S) + 1
                  ENDIF
                  RAMAN_UP(UT,V,S) = TNEW
                  IF ( DO_FILLING ) THEN
                    TOLD = ELASTIC_UP(UT,V,S)
                    TAZM = AZMFAC(I,UA)*ELASTIC_F_UP(UT,I,S)
                    ELASTIC_UP(UT,V,S) = TOLD + TAZM
                  ENDIF
                 ENDDO
                ENDDO
              ENDDO
            ENDIF
            IF ( DO_DNWELLING ) THEN
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = RAMAN_DN(UT,V,S)
                  TAZM = AZMFAC(I,UA)*RAMAN_F_DN(UT,I,S)
                  TNEW = TOLD + TAZM
                  IF ( TAZM .NE. ZERO ) THEN
                    IF ( ABS(TAZM/TNEW) .LT. ACCURACY ) THEN
                      COUNT(S) = COUNT(S) + 1
                    ENDIF
                  ELSE
                    COUNT(S)  = COUNT(S) + 1
                  ENDIF
                  RAMAN_DN(UT,V,S) = TNEW
                  IF ( DO_FILLING ) THEN
                    TOLD = ELASTIC_DN(UT,V,S)
                    TAZM = AZMFAC(I,UA)*ELASTIC_F_DN(UT,I,S)
                    ELASTIC_DN(UT,V,S) = TOLD + TAZM
                  ENDIF
                 ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO

!  Continuation point for elastic-only

 680      continue

          IF ( DO_ELASTIC_ONLY ) THEN
           IF ( DO_UPWELLING ) THEN
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = ELASTIC_UP(UT,V,S)
                  TAZM = AZMFAC(I,UA)*ELASTIC_F_UP(UT,I,S)
                  TNEW = TOLD + TAZM
                  IF ( TAZM .NE. ZERO ) THEN
                    IF ( ABS(TAZM/TNEW) .LT. ACCURACY ) THEN
                      COUNT(S) = COUNT(S) + 1
                    ENDIF
                  ELSE
                    COUNT(S)  = COUNT(S) + 1
                  ENDIF
                  ELASTIC_UP(UT,V,S) = TOLD + TAZM
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDIF
           IF ( DO_DNWELLING ) THEN
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = ELASTIC_DN(UT,V,S)
                  TAZM = AZMFAC(I,UA)*ELASTIC_F_DN(UT,I,S)
                  TNEW = TOLD + TAZM
                  IF ( TAZM .NE. ZERO ) THEN
                    IF ( ABS(TAZM/TNEW) .LT. ACCURACY ) THEN
                      COUNT(S) = COUNT(S) + 1
                    ENDIF
                  ELSE
                    COUNT(S)  = COUNT(S) + 1
                  ENDIF
                  ELASTIC_DN(UT,V,S) = TOLD + TAZM
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
           ENDIF
          ENDIF

!  Total number of passed criteria for all wavelengths

          COUNT_TOTAL = 0
          DO S = 1, NPOINTS_LOCAL
            COUNT_TOTAL = COUNT_TOTAL + COUNT(S)
          ENDDO

!  set convergence counter TESTCONV (convergence requires twice)

          IF ( COUNT_TOTAL .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
              LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED = FOURIER
            ENDIF
          ELSE
            TESTCONV = 0
          ENDIF

!  end convergence clause

        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER  .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED   = FOURIER
          ENDIF
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE RAMAN_CONVERGE

!  End module

      END MODULE lrrs_converge_m

