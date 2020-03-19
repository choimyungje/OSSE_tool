
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
! #               L_RAMAN_CONVERGE                              #
! #                                                             #
! #     New to Version 2.2 - Routines with linearizations       #
! #                                                             #
! ###############################################################


!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Bookkeeping improvements (use of "Only", clearer I/O specifications)
!    (2) Introduction of N_SURFACE_WFS & N_SLEAVE_WFS, extra dimensioning
!        for Surface Jacobians

!   -- Rob mod 5/12/17 for 2p5a, DO_SSFULL renamed --> DO_SSCORR_ALONE
!   -- Rob mod 5/12/17 for 2p5a, Use DO_SSCORR_GENERAL (OUTGOING or NADIR)

      MODULE lrrs_L_converge_m

!      USE LRRS_PARS_m, Only : LDU

      PRIVATE
      PUBLIC :: L_RAMAN_CONVERGE

      CONTAINS

      SUBROUTINE L_RAMAN_CONVERGE &
       ( DO_UPWELLING, DO_DNWELLING, DO_SSCORR_GENERAL, DO_SSCORR_ALONE,    & ! Inputs
         DO_NO_AZIMUTH, DO_FILLING, DO_ELASTIC_ONLY, FOURIER,               & ! Inputs
         DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,      & ! Inputs
         LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS, & ! Inputs
         NPOINTS_LOCAL, N_LOUTPUT, N_OUT_STREAMS, LOCAL_N_USERAZM,          & ! Inputs
         NLAYERS, N_GEOMETRIES, GEOM_OFFSETS, AZMFAC,                       & ! Inputs
         LC_ELASTIC_SS_UP, LC_ELASTIC_SS_DN, LP_ELASTIC_SS_UP, LP_ELASTIC_SS_DN, LS_ELASTIC_SS_UP, & ! E-SS Inputs
         LC_RAMAN_SS_UP,   LC_RAMAN_SS_DN,   LP_RAMAN_SS_UP,   LP_RAMAN_SS_DN,   LS_RAMAN_SS_UP,   & ! R-SS Inputs
         L_ELASTIC_F_UP,   L_ELASTIC_F_DN,   LS_ELASTIC_F_UP,  LS_ELASTIC_F_DN,                    & ! E-Fr Inputs
         L_RAMAN_F_UP,     L_RAMAN_F_DN,     LS_RAMAN_F_UP,    LS_RAMAN_F_DN,                      & ! R-Fr Inputs
         LC_ELASTIC_UP,    LC_ELASTIC_DN,    LP_ELASTIC_UP,    LP_ELASTIC_DN,    LS_ELASTIC_UP,    LS_ELASTIC_DN, & ! E Outputs
         LC_RAMAN_UP,      LC_RAMAN_DN,      LP_RAMAN_UP,      LP_RAMAN_DN,      LS_RAMAN_UP,      LS_RAMAN_DN,   & ! R Outputs
         LOCAL_ITERATION ) !  Diagnostic output

!  include file of dimensions and numbers 

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES, &
                              MAX_LAYERS, MAX_POINTS, MAX_ATMOSWFS, MAX_SURFACEWFS

      IMPLICIT NONE

!  Input arguments (Control)
!  =========================

!  Upwelling and downwelling flags

      LOGICAL, INTENT(IN) :: DO_UPWELLING
      LOGICAL, INTENT(IN) :: DO_DNWELLING

!  Azimuth only control

      LOGICAL, INTENT(IN) :: DO_NO_AZIMUTH

!  SS correction flags; Now incorporates the DB term
!   -- Rob mod 5/12/17 for 2p5a, Use DO_SSCORR_GENERAL (nadir or outgoing), rename SSFULL

      LOGICAL  , INTENT(IN) :: DO_SSCORR_GENERAL
      LOGICAL  , INTENT(IN) :: DO_SSCORR_ALONE

!  Filling flag..................
!  If set,  Fourier sums for elastic AND inelastic solutions are done

      LOGICAL, INTENT(IN) :: DO_FILLING

!  Elastic only flag

      LOGICAL, INTENT(IN) :: DO_ELASTIC_ONLY

!  Profile/column/surface linearization control inputs

      LOGICAL, INTENT(IN) :: do_profile_wfs
      LOGICAL, INTENT(IN) :: do_column_wfs
      LOGICAL, INTENT(IN) :: do_surface_wfs
      LOGICAL, INTENT(IN) :: do_sleave_wfs

!  Linearization control

!      LOGICAL, INTENT(IN) :: layer_vary_flag   (max_layers)
      INTEGER, INTENT(IN) :: layer_vary_number (max_layers)
      INTEGER, INTENT(IN) :: n_totalcolumn_wfs
      INTEGER, INTENT(IN) :: n_surface_wfs
      INTEGER, INTENT(IN) :: n_sleave_wfs

!  Number of points

      INTEGER, INTENT(IN) :: NPOINTS_LOCAL

!  Number of output streams

      INTEGER, INTENT(IN) :: N_OUT_STREAMS

!  Number of output levels

      INTEGER, INTENT(IN) :: N_LOUTPUT

!  Number of layers

      INTEGER, INTENT(IN) :: NLAYERS

!  Local number of azimuths

      INTEGER, INTENT(IN) :: LOCAL_N_USERAZM

!  Fourier index

      INTEGER, INTENT(IN) :: FOURIER

!  Number of Geometries, Geometry offsetting
!    --Added by R. Spurr, 04 August 2009

      INTEGER, INTENT(IN) :: N_GEOMETRIES
      INTEGER, INTENT(IN) :: GEOM_OFFSETS ( MAX_USER_STREAMS )

!  Azimuth angles. Adjusted factors, 18 March 2011.

      REAL(FPK), INTENT(IN) :: AZMFAC ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Input arguments (Elastic fields)
!  ================================

!  Fourier output of post-processed field
!  --------------------------------------

!  For the atmospheric weighting functions, post-processed values

      REAL(FPK), INTENT(IN) :: L_ELASTIC_F_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, &
                     MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_ELASTIC_F_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, &
                  MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  For the surface weighting functions, post-processed values

      REAL(FPK), INTENT(IN) :: LS_ELASTIC_F_UP &
            ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: LS_ELASTIC_F_DN &
            ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  single scatter results, elastic Jacobian fields
!  -----------------------------------------------

      REAL(FPK), INTENT(IN) :: LC_ELASTIC_SS_UP &
       (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(IN) :: LC_ELASTIC_SS_DN &
       (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(IN) :: LP_ELASTIC_SS_UP &
       (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(IN) :: LP_ELASTIC_SS_DN &
       (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(IN) :: LS_ELASTIC_SS_UP &
        ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Input arguments (Raman fields)
!  ==============================

!  Fourier output of post-processed field
!  --------------------------------------

!  For the atmospheric weighting functions, post-processed values

      REAL(FPK), INTENT(IN) :: L_RAMAN_F_UP ( MAX_ATMOSWFS, 0:MAX_LAYERS, &
                     MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: L_RAMAN_F_DN ( MAX_ATMOSWFS, 0:MAX_LAYERS, &
                  MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  For the surface weighting functions, post-processed values

      REAL(FPK), INTENT(IN) :: LS_RAMAN_F_UP &
            ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

      REAL(FPK), INTENT(IN) :: LS_RAMAN_F_DN &
            ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_USER_STREAMS, MAX_POINTS )

!  single scatter results, Raman Jacobian fields
!  ---------------------------------------------

      REAL(FPK), INTENT(IN) :: LC_RAMAN_SS_UP &
       (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(IN) :: LC_RAMAN_SS_DN &
       (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(IN) :: LP_RAMAN_SS_UP &
       (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)
      REAL(FPK), INTENT(IN) :: LP_RAMAN_SS_DN &
       (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(IN) :: LS_RAMAN_SS_UP &
        ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Output Arguments (Elastic fields)
!  =================================
!mick fix 10/19/2015 - changed intent to inout

!  Fourier summed output, Atmospheric Column Jacobians

      REAL(FPK), INTENT(INOUT) :: LC_ELASTIC_UP &
       (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LC_ELASTIC_DN &
       (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Fourier summed output, Atmospheric Profile Jacobians

      REAL(FPK), INTENT(INOUT) :: LP_ELASTIC_UP &
       (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LP_ELASTIC_DN &
       (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Fourier summed output, Surface Jacobians

      REAL(FPK), INTENT(INOUT) :: LS_ELASTIC_UP &
        ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: LS_ELASTIC_DN &
        ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  Output Arguments (Raman fields)
!  ===============================
!mick fix 10/19/2015 - changed intent to inout

!  Fourier summed output, Atmospheric Column Jacobians

      REAL(FPK), INTENT(INOUT) :: LC_RAMAN_UP &
       (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LC_RAMAN_DN &
       (MAX_ATMOSWFS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Fourier summed output, Atmospheric Profile Jacobians

      REAL(FPK), INTENT(INOUT) :: LP_RAMAN_UP &
       (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

      REAL(FPK), INTENT(INOUT) :: LP_RAMAN_DN &
       (MAX_ATMOSWFS,MAX_LAYERS,MAX_LOUTPUT,MAX_GEOMETRIES,MAX_POINTS)

!  Fourier summed output, Surface Jacobians

      REAL(FPK), INTENT(INOUT) :: LS_RAMAN_UP &
        ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

      REAL(FPK), INTENT(INOUT) :: LS_RAMAN_DN &
        ( MAX_SURFACEWFS, MAX_LOUTPUT, MAX_GEOMETRIES, MAX_POINTS )

!  modified/output variables
!  -------------------------

!  convergence flag

      LOGICAL, INTENT(INOUT) :: LOCAL_ITERATION

!  local variables
!  ---------------

      INTEGER   :: S, I, UT, T, UA, V, Q, K
      REAL(FPK) :: TOLD, TAZM

!  Single scattering ONLY (SSCORR_ALONE is flagged). Copy and return
!  -----------------------------------------------------------------

!  Copy Column Jacobians for SSCORR_ALONE

      IF ( DO_COLUMN_WFS ) THEN
        IF ( DO_SSCORR_GENERAL .AND. DO_SSCORR_ALONE ) THEN
          IF ( DO_UPWELLING ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO T = 1, N_LOUTPUT
                DO V = 1, N_GEOMETRIES
                  DO S = 1, NPOINTS_LOCAL
                    LC_ELASTIC_UP(Q,T,V,S) = LC_ELASTIC_SS_UP(Q,T,V,S)
                  ENDDO
                  IF ( .NOT. DO_ELASTIC_ONLY ) THEN
                    DO S = 1, NPOINTS_LOCAL
                      LC_RAMAN_UP(Q,T,V,S) = LC_RAMAN_SS_UP  (Q,T,V,S)
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO T = 1, N_LOUTPUT
                DO V = 1, N_GEOMETRIES
                  DO S = 1, NPOINTS_LOCAL
                    LC_ELASTIC_DN(Q,T,V,S) = LC_ELASTIC_SS_DN(Q,T,V,S)
                  ENDDO
                  IF ( .NOT. DO_ELASTIC_ONLY ) THEN
                    DO S = 1, NPOINTS_LOCAL
                      LC_RAMAN_DN(Q,T,V,S) = LC_RAMAN_SS_DN  (Q,T,V,S)
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!  Copy Profile Jacobians for SSCORR_ALONE

      IF ( DO_PROFILE_WFS ) THEN
        IF ( DO_SSCORR_GENERAL .AND. DO_SSCORR_ALONE ) THEN
          IF ( DO_UPWELLING ) THEN
            DO K = 1, NLAYERS
              DO Q = 1, LAYER_VARY_NUMBER(K)
                DO T = 1, N_LOUTPUT
                  DO V = 1, N_GEOMETRIES
                    DO S = 1, NPOINTS_LOCAL
                      LP_ELASTIC_UP(Q,K,T,V,S) = LP_ELASTIC_SS_UP(Q,K,T,V,S)
                    ENDDO
                    IF ( .NOT. DO_ELASTIC_ONLY ) THEN
                      DO S = 1, NPOINTS_LOCAL
                        LP_RAMAN_UP(Q,K,T,V,S) = LP_RAMAN_SS_UP(Q,K,T,V,S)
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            DO K = 1, NLAYERS
              DO Q = 1, LAYER_VARY_NUMBER(K)
                DO T = 1, N_LOUTPUT
                  DO V = 1, N_GEOMETRIES
                    DO S = 1, NPOINTS_LOCAL
                      LP_ELASTIC_DN(Q,K,T,V,S) = LP_ELASTIC_SS_DN(Q,K,T,V,S)
                    ENDDO
                    IF ( .NOT. DO_ELASTIC_ONLY ) THEN
                      DO S = 1, NPOINTS_LOCAL
                        LP_RAMAN_DN(Q,K,T,V,S) = LP_RAMAN_SS_DN(Q,K,T,V,S)
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!  Copy Surface Jacobians for SSCORR_ALONE

      IF ( DO_SURFACE_WFS ) THEN
        IF ( DO_SSCORR_GENERAL .AND. DO_SSCORR_ALONE ) THEN
          IF ( DO_UPWELLING ) THEN
            DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO T = 1, N_LOUTPUT
                DO V = 1, N_GEOMETRIES
                  DO S = 1, NPOINTS_LOCAL
                    LS_ELASTIC_UP(Q,T,V,S) = LS_ELASTIC_SS_UP(Q,T,V,S)
                  ENDDO
                  IF ( .NOT. DO_ELASTIC_ONLY ) THEN
                    DO S = 1, NPOINTS_LOCAL
                      LS_RAMAN_UP(Q,T,V,S) = LS_RAMAN_SS_UP(Q,T,V,S)
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO T = 1, N_LOUTPUT
                DO V = 1, N_GEOMETRIES
                  DO S = 1, NPOINTS_LOCAL
                    LS_ELASTIC_DN(Q,T,V,S) = ZERO
                  ENDDO
                  IF ( .NOT. DO_ELASTIC_ONLY ) THEN
                    DO S = 1, NPOINTS_LOCAL
                      LS_RAMAN_DN(Q,T,V,S)   = ZERO
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!  Exit (All finished) if DO_SSCORR_GENERAL .AND. DO_SSCORR_ALONE

      IF ( DO_SSCORR_GENERAL .AND. DO_SSCORR_ALONE ) RETURN

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Diffuse term addition
!  ---------------------

!  Skip if elastic only

        IF ( DO_ELASTIC_ONLY ) GO TO 677

!  Column Jacobians

        IF ( DO_COLUMN_WFS .AND. DO_UPWELLING ) THEN
         DO Q = 1, N_TOTALCOLUMN_WFS
          DO T = 1, N_LOUTPUT
           DO I = 1, N_OUT_STREAMS
            DO UA = 1, LOCAL_N_USERAZM
             V = GEOM_OFFSETS(I) + UA
             DO S = 1, NPOINTS_LOCAL
              LC_RAMAN_UP(Q,T,V,S) = L_RAMAN_F_UP(Q,0,T,I,S)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF

        IF ( DO_COLUMN_WFS .AND. DO_DNWELLING ) THEN
         DO Q = 1, N_TOTALCOLUMN_WFS
          DO T = 1, N_LOUTPUT
           DO I = 1, N_OUT_STREAMS
            DO UA = 1, LOCAL_N_USERAZM
             V = GEOM_OFFSETS(I) + UA
             DO S = 1, NPOINTS_LOCAL
              LC_RAMAN_DN(Q,T,V,S) = L_RAMAN_F_DN(Q,0,T,I,S)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF

!  Profile Jacobians

        IF ( DO_PROFILE_WFS .AND. DO_UPWELLING ) THEN
         DO K = 1, NLAYERS
          DO Q = 1, LAYER_VARY_NUMBER(K)
           DO T = 1, N_LOUTPUT
            DO I = 1, N_OUT_STREAMS
             DO UA = 1, LOCAL_N_USERAZM
              V = GEOM_OFFSETS(I) + UA
              DO S = 1, NPOINTS_LOCAL
               LP_RAMAN_UP(Q,K,T,V,S) = L_RAMAN_F_UP(Q,K,T,I,S)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF

        IF ( DO_PROFILE_WFS .AND. DO_DNWELLING ) THEN
         DO K = 1, NLAYERS
          DO Q = 1, LAYER_VARY_NUMBER(K)
           DO T = 1, N_LOUTPUT
            DO I = 1, N_OUT_STREAMS
             DO UA = 1, LOCAL_N_USERAZM
              V = GEOM_OFFSETS(I) + UA
              DO S = 1, NPOINTS_LOCAL
               LP_RAMAN_DN(Q,K,T,V,S) = L_RAMAN_F_DN(Q,K,T,I,S)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF

!  Surface Jacobians

        IF ( DO_SURFACE_WFS .AND. DO_UPWELLING ) THEN
          DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
            DO T = 1, N_LOUTPUT
              DO I = 1, N_OUT_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = GEOM_OFFSETS(I) + UA
                  DO S = 1, NPOINTS_LOCAL
                    LS_RAMAN_UP(Q,T,V,S) = LS_RAMAN_F_UP(Q,T,I,S)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        IF ( DO_SURFACE_WFS .AND. DO_DNWELLING ) THEN
          DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
            DO T = 1, N_LOUTPUT
              DO I = 1, N_OUT_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = GEOM_OFFSETS(I) + UA
                  DO S = 1, NPOINTS_LOCAL
                    LS_RAMAN_DN(Q,T,V,S) = LS_RAMAN_F_DN(Q,T,I,S)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Repeat for the elastic only calculation (if FILLING is flagged)
!  ---------------------------------------------------------------

!  Continuation point for Elastic-only

 677    CONTINUE

!  Skip if not flagged (no filling, no calculation)

        IF ( DO_FILLING .OR. DO_ELASTIC_ONLY ) THEN

!  Column Jacobians

         IF ( DO_COLUMN_WFS .AND. DO_UPWELLING )THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
           DO T = 1, N_LOUTPUT
            DO I = 1, N_OUT_STREAMS
             DO UA = 1, LOCAL_N_USERAZM
              V = GEOM_OFFSETS(I) + UA
              DO S = 1, NPOINTS_LOCAL
               LC_ELASTIC_UP(Q,T,V,S) = L_ELASTIC_F_UP(Q,0,T,I,S)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF

         IF ( DO_COLUMN_WFS .AND. DO_DNWELLING ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
           DO T = 1, N_LOUTPUT
            DO I = 1, N_OUT_STREAMS
             DO UA = 1, LOCAL_N_USERAZM
              V = GEOM_OFFSETS(I) + UA
              DO S = 1, NPOINTS_LOCAL
               LC_ELASTIC_DN(Q,T,V,S) = L_ELASTIC_F_DN(Q,0,T,I,S)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF

!  Profile Jacobians

         IF ( DO_PROFILE_WFS .AND. DO_UPWELLING ) THEN
          DO K = 1, NLAYERS
           DO Q = 1, LAYER_VARY_NUMBER(K)
            DO T = 1, N_LOUTPUT
             DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = GEOM_OFFSETS(I) + UA
               DO S = 1, NPOINTS_LOCAL
                LP_ELASTIC_UP(Q,K,T,V,S) = L_ELASTIC_F_UP(Q,K,T,I,S)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF

         IF ( DO_PROFILE_WFS .AND. DO_DNWELLING ) THEN
          DO K = 1, NLAYERS
           DO Q = 1, LAYER_VARY_NUMBER(K)
            DO T = 1, N_LOUTPUT
             DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = GEOM_OFFSETS(I) + UA
               DO S = 1, NPOINTS_LOCAL
                LP_ELASTIC_DN(Q,K,T,V,S) = L_ELASTIC_F_DN(Q,K,T,I,S)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF

!  Surface Jacobians

         IF ( DO_SURFACE_WFS .AND. DO_UPWELLING ) THEN
           DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
            DO T = 1, N_LOUTPUT
             DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = GEOM_OFFSETS(I) + UA
               DO S = 1, NPOINTS_LOCAL
                LS_ELASTIC_UP(Q,T,V,S) = LS_ELASTIC_F_UP(Q,T,I,S)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
         ENDIF

         IF ( DO_SURFACE_WFS .AND. DO_DNWELLING ) THEN
           DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
            DO T = 1, N_LOUTPUT
             DO I = 1, N_OUT_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
               V = GEOM_OFFSETS(I) + UA
               DO S = 1, NPOINTS_LOCAL
                LS_ELASTIC_DN(Q,T,V,S) = LS_ELASTIC_F_DN(Q,T,I,S)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
           ENDDO
         ENDIF

!  End elastic contribution

        ENDIF

!  Add the single scatter component if flagged.
!  -------------------------------------------

!  Skip if elastic only

        IF ( DO_ELASTIC_ONLY ) GO TO 678

!  OLD STYLE NADIR CORRECTION IS REMOVED.......

!        IF ( DO_SSCORR_NADIR ) THEN
!         IF ( DO_COLUMN_WFS .AND. DO_UPWELLING ) THEN
!           DO Q = 1, N_TOTALCOLUMN_WFS
!             DO T = 1, N_LOUTPUT
!               DO V = 1, N_GEOMETRIES
!                 DO S = 1, NPOINTS_LOCAL
!                   LC_RAMAN_UP(Q,T,V,S) =
!     &             LC_RAMAN_UP(Q,T,V,S) + LC_ELASTIC_SS_UP(Q,T,V,S)
!                 ENDDO
!               ENDDO
!             ENDDO
!           ENDDO
!         ENDIF
!         IF ( DO_COLUMN_WFS .AND. DO_DNWELLING ) THEN
!           DO Q = 1, N_TOTALCOLUMN_WFS
!             DO T = 1, N_LOUTPUT
!               DO V = 1, N_GEOMETRIES
!                 DO S = 1, NPOINTS_LOCAL
!                   LC_RAMAN_DN(Q,T,V,S) = LC_RAMAN_DN(Q,T,V,S)
!     &                        + LC_ELASTIC_SS_DN(Q,T,V,S)
!                 ENDDO
!               ENDDO
!             ENDDO
!           ENDDO
!         ENDIF
!         IF ( DO_PROFILE_WFS .AND. DO_UPWELLING ) THEN
!           DO K = 1, NLAYERS
!             DO Q = 1, LAYER_VARY_NUMBER(K)
!               DO T = 1, N_LOUTPUT
!                 DO V = 1, N_GEOMETRIES
!                   DO S = 1, NPOINTS_LOCAL
!                     LP_RAMAN_UP(Q,K,T,V,S) = LP_RAMAN_UP(Q,K,T,V,S)
!     &                          + LP_ELASTIC_SS_UP(Q,K,T,V,S)
!                   ENDDO
!                 ENDDO
!               ENDDO
!             ENDDO
!           ENDDO
!         ENDIF
!         IF ( DO_PROFILE_WFS .AND. DO_DNWELLING ) THEN
!           DO K = 1, NLAYERS
!             DO Q = 1, LAYER_VARY_NUMBER(K)
!               DO T = 1, N_LOUTPUT
!                 DO V = 1, N_GEOMETRIES
!                   DO S = 1, NPOINTS_LOCAL
!                     LP_RAMAN_DN(Q,K,T,V,S) = LP_RAMAN_DN(Q,K,T,V,S)
!     &                          + LP_ELASTIC_SS_DN(Q,K,T,V,S)
!                   ENDDO
!                 ENDDO
!               ENDDO
!             ENDDO
!           ENDDO
!         ENDIF
!         IF ( DO_SURFACE_WFS .AND. DO_UPWELLING ) THEN
!           DO T = 1, N_LOUTPUT
!             DO V = 1, N_GEOMETRIES
!               DO S = 1, NPOINTS_LOCAL
!                 LS_RAMAN_UP(T,V,S) =
!     &            LS_RAMAN_UP(T,V,S) + LS_ELASTIC_SS_UP(T,V,S)
!               ENDDO
!             ENDDO
!           ENDDO
!         ENDIF
!         IF ( DO_SURFACE_WFS .AND. DO_DNWELLING ) THEN
!           DO T = 1, N_LOUTPUT
!             DO V = 1, N_GEOMETRIES
!               DO S = 1, NPOINTS_LOCAL
!                 LS_RAMAN_DN(T,V,S) =
!     &            LS_RAMAN_DN(T,V,S) + LS_ELASTIC_SS_DN(T,V,S)
!               ENDDO
!             ENDDO
!           ENDDO
!         ENDIF
!        ENDIF

!  OUTGOING (new-style) correction:
!  --------------------------------

        IF ( DO_SSCORR_GENERAL ) THEN

!  Converged column Jacobians

         IF ( DO_COLUMN_WFS .AND. DO_UPWELLING ) THEN
           DO Q = 1, N_TOTALCOLUMN_WFS
             DO T = 1, N_LOUTPUT
               DO V = 1, N_GEOMETRIES
                 DO S = 1, NPOINTS_LOCAL
                   LC_RAMAN_UP(Q,T,V,S) = &
                     LC_RAMAN_UP(Q,T,V,S) + LC_RAMAN_SS_UP(Q,T,V,S)
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDIF
         IF ( DO_COLUMN_WFS .AND. DO_DNWELLING ) THEN
           DO Q = 1, N_TOTALCOLUMN_WFS
             DO T = 1, N_LOUTPUT
               DO V = 1, N_GEOMETRIES
                 DO S = 1, NPOINTS_LOCAL
                   LC_RAMAN_DN(Q,T,V,S) = &
                     LC_RAMAN_DN(Q,T,V,S) + LC_RAMAN_SS_DN(Q,T,V,S)
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDIF

!  Converged profile Jacobians

         IF ( DO_PROFILE_WFS .AND. DO_UPWELLING ) THEN
           DO K = 1, NLAYERS
             DO Q = 1, LAYER_VARY_NUMBER(K)
               DO T = 1, N_LOUTPUT
                 DO V = 1, N_GEOMETRIES
                   DO S = 1, NPOINTS_LOCAL
                     LP_RAMAN_UP(Q,K,T,V,S) = &
                       LP_RAMAN_UP(Q,K,T,V,S) + LP_RAMAN_SS_UP(Q,K,T,V,S)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDIF
         IF ( DO_PROFILE_WFS .AND. DO_DNWELLING ) THEN
           DO K = 1, NLAYERS
             DO Q = 1, LAYER_VARY_NUMBER(K)
               DO T = 1, N_LOUTPUT
                 DO V = 1, N_GEOMETRIES
                   DO S = 1, NPOINTS_LOCAL
                     LP_RAMAN_DN(Q,K,T,V,S) = &
                       LP_RAMAN_DN(Q,K,T,V,S) + LP_RAMAN_SS_DN(Q,K,T,V,S)
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDIF

!  Converged Surface Jacobians (Upwelling only)

         IF ( DO_SURFACE_WFS .AND. DO_UPWELLING ) THEN
           DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
             DO T = 1, N_LOUTPUT
               DO V = 1, N_GEOMETRIES
                 DO S = 1, NPOINTS_LOCAL
                   LS_RAMAN_UP(Q,T,V,S) = &
                     LS_RAMAN_UP(Q,T,V,S) + LS_RAMAN_SS_UP(Q,T,V,S)
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDIF

!  End SSCORR General correction clause

        ENDIF

!  Continuation point for Elastic-only

 678    continue

!  ELASTIC correction (same for both SS methods):
!  ----------------------------------------------

!   Only if filling flag is set, or elastic-only

!        IF ( DO_SSCORR_GENERAL .or. DO_SSCORR_NADIR ) THEN

        IF (  DO_SSCORR_GENERAL ) THEN
         IF ( DO_FILLING .OR. DO_ELASTIC_ONLY ) THEN

!  Converged column Jacobians

          IF ( DO_COLUMN_WFS .AND. DO_UPWELLING ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO T = 1, N_LOUTPUT
                DO V = 1, N_GEOMETRIES
                  DO S = 1, NPOINTS_LOCAL
                    LC_ELASTIC_UP(Q,T,V,S) = &
                      LC_ELASTIC_UP(Q,T,V,S) + LC_ELASTIC_SS_UP(Q,T,V,S)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF ( DO_COLUMN_WFS .AND. DO_DNWELLING ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO T = 1, N_LOUTPUT
                DO V = 1, N_GEOMETRIES
                  DO S = 1, NPOINTS_LOCAL
                    LC_ELASTIC_DN(Q,T,V,S) = &
                      LC_ELASTIC_DN(Q,T,V,S) + LC_ELASTIC_SS_DN(Q,T,V,S)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  Converged profile Jacobians

          IF ( DO_PROFILE_WFS .AND. DO_UPWELLING ) THEN
            DO K = 1, NLAYERS
              DO Q = 1, LAYER_VARY_NUMBER(K)
                DO T = 1, N_LOUTPUT
                  DO V = 1, N_GEOMETRIES
                    DO S = 1, NPOINTS_LOCAL
                      LP_ELASTIC_UP(Q,K,T,V,S) = &
                        LP_ELASTIC_UP(Q,K,T,V,S) + LP_ELASTIC_SS_UP(Q,K,T,V,S)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF ( DO_PROFILE_WFS .AND. DO_DNWELLING ) THEN
            DO K = 1, NLAYERS
              DO Q = 1, LAYER_VARY_NUMBER(K)
                DO T = 1, N_LOUTPUT
                  DO V = 1, N_GEOMETRIES
                    DO S = 1, NPOINTS_LOCAL
                      LP_ELASTIC_DN(Q,K,T,V,S) = &
                        LP_ELASTIC_DN(Q,K,T,V,S)+LP_ELASTIC_SS_DN(Q,K,T,V,S)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  Surface Jacobians (upwelling only)

          IF ( DO_SURFACE_WFS .AND. DO_UPWELLING ) THEN
            DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
              DO T = 1, N_LOUTPUT
                DO V = 1, N_GEOMETRIES
                  DO S = 1, NPOINTS_LOCAL
                    LS_ELASTIC_UP(Q,T,V,S) = &
                      LS_ELASTIC_UP(Q,T,V,S) + LS_ELASTIC_SS_UP(Q,T,V,S)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  End Filling/Elastic-only and SS correction clauses

         ENDIF
        ENDIF

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  Skip if elastic only

        IF ( DO_ELASTIC_ONLY ) GO TO 679

!  For each azimuth, add Fourier component
!     - for direction, user output choice, out stream

!  Column Jacobians

        IF ( DO_COLUMN_WFS .AND. DO_UPWELLING ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = LC_RAMAN_UP(Q,UT,V,S)
                  TAZM = AZMFAC(I,UA)*L_RAMAN_F_UP(Q,0,UT,I,S)
                  LC_RAMAN_UP(Q,UT,V,S) = TOLD + TAZM
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_COLUMN_WFS .AND. DO_DNWELLING ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = LC_RAMAN_DN(Q,UT,V,S)
                  TAZM = AZMFAC(I,UA)*L_RAMAN_F_DN(Q,0,UT,I,S)
                  LC_RAMAN_DN(Q,UT,V,S) = TOLD + TAZM
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Profile Jacobians

        IF ( DO_PROFILE_WFS .AND. DO_UPWELLING ) THEN
          DO K = 1, NLAYERS
            DO Q = 1, LAYER_VARY_NUMBER(K)
             DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
               DO I = 1, N_OUT_STREAMS
                V = GEOM_OFFSETS(I) + UA
                DO S = 1, NPOINTS_LOCAL
                 TOLD = LP_RAMAN_UP(Q,K,UT,V,S)
                 TAZM = AZMFAC(I,UA)*L_RAMAN_F_UP(Q,K,UT,I,S)
                 LP_RAMAN_UP(Q,K,UT,V,S) = TOLD + TAZM
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_PROFILE_WFS .AND. DO_DNWELLING ) THEN
          DO K = 1, NLAYERS
            DO Q = 1, LAYER_VARY_NUMBER(K)
             DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
               DO I = 1, N_OUT_STREAMS
                V = GEOM_OFFSETS(I) + UA
                DO S = 1, NPOINTS_LOCAL
                 TOLD = LP_RAMAN_DN(Q,K,UT,V,S)
                 TAZM = AZMFAC(I,UA)*L_RAMAN_F_DN(Q,K,UT,I,S)
                 LP_RAMAN_DN(Q,K,UT,V,S) = TOLD + TAZM
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Surface Jacobians

        IF ( DO_SURFACE_WFS .AND. DO_UPWELLING ) THEN
          DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                  V = GEOM_OFFSETS(I) + UA
                  DO S = 1, NPOINTS_LOCAL
                    TOLD = LS_RAMAN_UP(Q,UT,V,S)
                    TAZM = AZMFAC(I,UA)*LS_RAMAN_F_UP(Q,UT,I,S)
                    LS_RAMAN_UP(Q,UT,V,S) = TOLD + TAZM
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_SURFACE_WFS .AND. DO_DNWELLING ) THEN
          DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                  V = GEOM_OFFSETS(I) + UA
                  DO S = 1, NPOINTS_LOCAL
                    TOLD = LS_RAMAN_DN(Q,UT,V,S)
                    TAZM = AZMFAC(I,UA)*LS_RAMAN_F_DN(Q,UT,I,S)
                    LS_RAMAN_DN(Q,UT,V,S) = TOLD + TAZM
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Repeat for the elastic only calculation (if FILLING is flagged)
!  ---------------------------------------

!  Elastic-only continuation point

 679    continue

!  Skip if no filling or calculation

        IF ( .NOT.DO_FILLING.AND..NOT.DO_ELASTIC_ONLY ) GO TO 680

!  Column Jacobians

        IF ( DO_COLUMN_WFS .AND. DO_UPWELLING ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = LC_ELASTIC_UP(Q,UT,V,S)
                  TAZM = AZMFAC(I,UA)*L_ELASTIC_F_UP(Q,0,UT,I,S)
                  LC_ELASTIC_UP(Q,UT,V,S) = TOLD + TAZM
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_COLUMN_WFS .AND. DO_DNWELLING ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                 V = GEOM_OFFSETS(I) + UA
                 DO S = 1, NPOINTS_LOCAL
                  TOLD = LC_ELASTIC_DN(Q,UT,V,S)
                  TAZM = AZMFAC(I,UA)*L_ELASTIC_F_DN(Q,0,UT,I,S)
                  LC_ELASTIC_DN(Q,UT,V,S) = TOLD + TAZM
                 ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Profile Jacobians

        IF ( DO_PROFILE_WFS .AND. DO_UPWELLING ) THEN
          DO K = 1, NLAYERS
            DO Q = 1, LAYER_VARY_NUMBER(K)
             DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
               DO I = 1, N_OUT_STREAMS
                V = GEOM_OFFSETS(I) + UA
                DO S = 1, NPOINTS_LOCAL
                 TOLD = LP_ELASTIC_UP(Q,K,UT,V,S)
                 TAZM = AZMFAC(I,UA)*L_ELASTIC_F_UP(Q,K,UT,I,S)
                 LP_ELASTIC_UP(Q,K,UT,V,S) = TOLD + TAZM
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_PROFILE_WFS .AND. DO_DNWELLING ) THEN
          DO K = 1, NLAYERS
            DO Q = 1, LAYER_VARY_NUMBER(K)
             DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
               DO I = 1, N_OUT_STREAMS
                V = GEOM_OFFSETS(I) + UA
                DO S = 1, NPOINTS_LOCAL
                 TOLD = LP_ELASTIC_DN(Q,K,UT,V,S)
                 TAZM = AZMFAC(I,UA)*L_ELASTIC_F_DN(Q,K,UT,I,S)
                 LP_ELASTIC_DN(Q,K,UT,V,S) = TOLD + TAZM
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Surface Jacobians

        IF ( DO_SURFACE_WFS .AND. DO_UPWELLING ) THEN
          DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                  V = GEOM_OFFSETS(I) + UA
                  DO S = 1, NPOINTS_LOCAL
                    TOLD = LS_ELASTIC_UP(Q,UT,V,S)
                    TAZM = AZMFAC(I,UA)*LS_ELASTIC_F_UP(Q,UT,I,S)
                    LS_ELASTIC_UP(Q,UT,V,S) = TOLD + TAZM
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_SURFACE_WFS .AND. DO_DNWELLING ) THEN
          DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
            DO UA = 1, LOCAL_N_USERAZM
              DO UT = 1, N_LOUTPUT
                DO I = 1, N_OUT_STREAMS
                  V = GEOM_OFFSETS(I) + UA
                  DO S = 1, NPOINTS_LOCAL
                    TOLD = LS_ELASTIC_DN(Q,UT,V,S)
                    TAZM = AZMFAC(I,UA)*LS_ELASTIC_F_DN(Q,UT,I,S)
                    LS_ELASTIC_DN(Q,UT,V,S) = TOLD + TAZM
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Continuation point for skipping elastic calculation

 680    CONTINUE

!  end convergence clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE L_RAMAN_CONVERGE

!  End module

      END MODULE lrrs_L_converge_m

