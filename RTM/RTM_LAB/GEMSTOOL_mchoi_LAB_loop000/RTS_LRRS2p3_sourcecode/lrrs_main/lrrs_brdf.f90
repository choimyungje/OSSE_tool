! ###############################################################
! #                                                             #
! #                    THE LIDORT_RRS MODEL                     #
! #                                                             #
! #      (LInearized Discrete Ordinate Radiative Transfer)      #
! #       --         -        -        -         -              #
! #                 (Rotational Raman Scatter)                  #
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
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Version      :  2.3                                        #
! #  Release Date :  March 2011                                 #
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
! # Subroutines in this Module (all new for  Version 2.1)       #
! #                                                             #
! #            LRRS_BRDF_MASTER (master), calling               #
! #              LRRS_BRDF_MAKER                                #
! #              LRRS_BRDF_FUNCTION                             #
! #                                                             #
! #            LRRS_BRDF_FOURIER (called by FOURIER)            #
! #                                                             #
! # Choice of BRDFs:                                            #
! #                                                             #
! #            HAPKE_FUNCTION                                   #
! #            RAHMAN_FUNCTION                                  #
! #            COXMUNK_FUNCTION                                 #
! #            COXMUNK_FUNCTION_DB                              #
! #                                                             #
! #   Note (*): Version 2.1 with BRDF setup, November 2007.     #
! #                                                             #
! ###############################################################


      MODULE lrrs_brdf

      USE LRRS_PARS

      PRIVATE
      PUBLIC  :: LRRS_BRDF_MASTER,   &
                 LRRS_BRDF_MAKER,    &
                 LRRS_BRDF_FUNCTION, &
                 LRRS_BRDF_FOURIER,  &
                 HAPKE_FUNCTION,     &
                 RAHMAN_FUNCTION,    &
                 COXMUNK_FUNCTION,   &
                 COXMUNK_FUNCTION_DB,&
                 DERFC_E

!  @@@ RobFix 5/5/11. Drop the External calls.
!    Need now the new routine, LRRS_BRDF_FUNCTION, modeled after LIDORT.

      CONTAINS

      SUBROUTINE LRRS_BRDF_MASTER &
          ( WHICH_BRDF, DO_GLITTER_DBMS, BRDF_PARS, &
               DO_UPWELLING, DO_USER_STREAMS, &
               DO_SSCORR_OUTGOING, DO_SSFULL, &
               NSTREAMS_BRDF,  SOLAR_ANGLE, &
               NSTREAMS,       QUAD_STREAMS, &
               N_USER_STREAMS, USER_STREAMS, &
               N_USER_RELAZMS, USER_RELAZMS, &
               X_BRDF, A_BRDF, EXACTDB_BRDFUNC, &
               BRDFUNC, BRDFUNC_0, &
               USER_BRDFUNC, USER_BRDFUNC_0 )

!  Master routine for creation of BRDFs
!  @@@ RobFix 5/5/11. Drop the External calls

!  include file of dimensions and numbers

      USE LRRS_PARS
      USE LRRS_AUX2

      IMPLICIT NONE

!  Input arguments
!  ===============

!  BRDF parameters

      INTEGER, INTENT(IN) ::          WHICH_BRDF
      LOGICAL, INTENT(IN) ::          DO_GLITTER_DBMS
      REAL(FPK), INTENT(IN) :: BRDF_PARS(3)

!  Flags

      LOGICAL, INTENT(IN) ::          DO_UPWELLING
      LOGICAL, INTENT(IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT(IN) ::          DO_SSFULL

!  Geometrical variables

      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          NSTREAMS_BRDF
      INTEGER, INTENT(IN) ::          N_USER_STREAMS
      INTEGER, INTENT(IN) ::          N_USER_RELAZMS

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLE
      REAL(FPK), INTENT(IN) :: QUAD_STREAMS ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  Output arguments (BRDFs)
!  ========================

!  BRDF quadratures

      REAL(FPK), INTENT(OUT) :: X_BRDF  ( MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(OUT) :: A_BRDF  ( MAX_STREAMS_BRDF )

!  at quadrature (discrete ordinate) angles

      REAL(FPK), INTENT(OUT) :: BRDFUNC &
            ( MAX_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK), INTENT(OUT) :: BRDFUNC_0 &
            ( MAX_STREAMS, MAX_STREAMS_BRDF )

!  at user-defined stream directions

      REAL(FPK), INTENT(OUT) :: USER_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK), INTENT(OUT) :: USER_BRDFUNC_0 &
            ( MAX_USER_STREAMS, MAX_STREAMS_BRDF )

!  Exact DB values

      REAL(FPK), INTENT(OUT) :: EXACTDB_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Local definitions
!  =================

!  BRDF functions, External subroutine
!  Single Kernel functions. No 3-kernel implementation here.

!      EXTERNAL       HAPKE_FUNCTION
!      EXTERNAL       RAHMAN_FUNCTION
!      EXTERNAL       COXMUNK_FUNCTION
!      EXTERNAL       COXMUNK_FUNCTION_DB

!  Hapke old uses exact DISORT code
!      EXTERNAL       HAPKE_FUNCTION_OLD

!  other variables

      INTEGER ::          I, I1, NBRDF_HALF
      REAL(FPK) :: CX_BRDF ( MAX_STREAMS_BRDF )
      REAL(FPK) :: SX_BRDF ( MAX_STREAMS_BRDF )

!  BRDF quadrature (Gauss-Legendre)
!  ---------------

!  Save these quantities for efficient coding

      NBRDF_HALF = NSTREAMS_BRDF / 2
      CALL GAULEG ( ZERO, ONE, X_BRDF, A_BRDF, NBRDF_HALF )
      DO I = 1, NBRDF_HALF
        I1 = I + NBRDF_HALF
        X_BRDF(I1) = - X_BRDF(I)
        A_BRDF(I1) =   A_BRDF(I)
      ENDDO
      DO I = 1, NSTREAMS_BRDF
        X_BRDF(I) = PIE * X_BRDF(I)
        CX_BRDF(I) = DCOS ( X_BRDF(I) )
        SX_BRDF(I) = DSIN ( X_BRDF(I) )
      ENDDO

!  Fill BRDF arrays
!  ----------------

!  @@@ RobFix 5/5/11. Drop the External calls

      CALL LRRS_BRDF_MAKER &
             ( WHICH_BRDF, BRDF_PARS, &
               DO_UPWELLING, DO_USER_STREAMS, DO_GLITTER_DBMS, &
               DO_SSCORR_OUTGOING, DO_SSFULL, &
               NSTREAMS_BRDF,  X_BRDF, CX_BRDF, SOLAR_ANGLE, &
               NSTREAMS,       QUAD_STREAMS, &
               N_USER_STREAMS, USER_STREAMS, &
               N_USER_RELAZMS, USER_RELAZMS, &
               EXACTDB_BRDFUNC, BRDFUNC, BRDFUNC_0, &
               USER_BRDFUNC, USER_BRDFUNC_0 )

!  Fake Lambertian test with albedo 0.05
!        exactdb_brdfunc(1,1) = 0.05d0
!        do j = 1, 8
!         do k = 1, 50
!          brdfunc_0(j,k) = 0.05d0
!          user_brdfunc_0(1,k) = 0.05d0
!          do i = 1, 8
!           brdfunc(i,j,k) = 0.05d0
!          enddo
!          user_brdfunc(1,j,k) = 0.05d0
!         enddo
!        enddo
!       write(*,'(i2,1p8e12.4)')i,(brdfunc(i,j,1),j=1,8)

!  Finish

      RETURN
      END SUBROUTINE LRRS_BRDF_MASTER

!

      SUBROUTINE LRRS_BRDF_MAKER &
             ( WHICH_BRDF, BRDF_PARS, &
               DO_UPWELLING, DO_USER_STREAMS, DO_GLITTER_DBMS, &
               DO_SSCORR_OUTGOING, DO_SSFULL, &
               NSTREAMS_BRDF,  X_BRDF, CX_BRDF, SOLAR_ANGLE, &
               NSTREAMS,       QUAD_STREAMS, &
               N_USER_STREAMS, USER_STREAMS, &
               N_USER_RELAZMS, USER_RELAZMS, &
               EXACTDB_BRDFUNC, BRDFUNC, BRDFUNC_0, &
               USER_BRDFUNC, USER_BRDFUNC_0 )

!  Prepares the bidirectional reflectance scatter matrices
!  @@@ RobFix 5/5/11. Drop the External calls

!  include file of dimensions and numbers

      USE LRRS_PARS

      IMPLICIT NONE

!  Input arguments
!  ===============

!  @@@ RobFix 5/5/11. Drop external calls. Use index

      INTEGER, INTENT(IN) ::   WHICH_BRDF

!  BRDF functions (external calls). BRDF parameters
!      EXTERNAL         BRDF_FUNCTION
!      EXTERNAL         BRDF_FUNCTION_DB

      REAL(FPK), INTENT(IN) :: BRDF_PARS(3)

!  Flags

      LOGICAL, INTENT(IN) ::          DO_UPWELLING
      LOGICAL, INTENT(IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::          DO_GLITTER_DBMS
      LOGICAL, INTENT(IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT(IN) ::          DO_SSFULL

!  Geometrical variables

      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          NSTREAMS_BRDF
      INTEGER, INTENT(IN) ::          N_USER_STREAMS
      INTEGER, INTENT(IN) ::          N_USER_RELAZMS

      REAL(FPK), INTENT(IN) :: SOLAR_ANGLE
      REAL(FPK), INTENT(IN) :: QUAD_STREAMS ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(FPK), INTENT(IN) :: USER_RELAZMS ( MAX_USER_RELAZMS )
      REAL(FPK), INTENT(IN) :: X_BRDF  ( MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(IN) :: CX_BRDF ( MAX_STREAMS_BRDF )

!  Output arguments (BRDFs)
!  ========================

!  at quadrature (discrete ordinate) angles

      REAL(FPK), INTENT(OUT) :: BRDFUNC &
            ( MAX_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK), INTENT(OUT) :: BRDFUNC_0 &
            ( MAX_STREAMS, MAX_STREAMS_BRDF )

!  at user-defined stream directions

      REAL(FPK), INTENT(OUT) :: USER_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK), INTENT(OUT) :: USER_BRDFUNC_0 &
            ( MAX_USER_STREAMS, MAX_STREAMS_BRDF )

!  Exact DB values

      REAL(FPK), INTENT(OUT) :: EXACTDB_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  local variables
!  --------------

      LOGICAL ::          DO_DB_CORRECTION
      INTEGER ::          I, UI, J, K
      REAL(FPK) :: MUX, SZASURCOS,SZASURSIN
      REAL(FPK) :: USER_SINES ( MAX_USER_STREAMS )
      REAL(FPK) :: QUAD_SINES ( MAX_STREAMS )
      REAL(FPK) :: PHIANG, COSPHI

!  Local Do Db correction flag

!       DO_DB_CORRECTION = DO_UPWELLING .AND.
!     &    ( DO_SSCORR_OUTGOING .OR. DO_SSCORR_NADIR )

       DO_DB_CORRECTION = DO_UPWELLING .AND. DO_SSCORR_OUTGOING

!  Solar cosine/sine

      MUX =  DCOS(SOLAR_ANGLE * DEG_TO_RAD)
      SZASURCOS = MUX
      SZASURSIN = DSQRT(ONE-MUX*MUX)

!  quadrature sines

      DO I = 1, NSTREAMS
        QUAD_SINES(I)    = DSQRT(ONE-QUAD_STREAMS(I)*QUAD_STREAMS(I))
      ENDDO

!  quadrature sines

      DO I = 1, N_USER_STREAMS
        USER_SINES(I)    = DSQRT(ONE-USER_STREAMS(I)*USER_STREAMS(I))
      ENDDO

!  Exact DB calculation
!  --------------------

      IF ( DO_DB_CORRECTION ) THEN
        IF ( WHICH_BRDF .eq. COXMUNK_IDX .and. DO_GLITTER_DBMS ) THEN
          DO K = 1, N_USER_RELAZMS
            PHIANG = USER_RELAZMS(K)
            COSPHI = DCOS(PHIANG*DEG_TO_RAD)
            DO UI = 1, N_USER_STREAMS
              CALL COXMUNK_FUNCTION_DB &
                ( MAX_BRDF_PARAMETERS, BRDF_PARS, &
                  SZASURCOS, SZASURSIN, &
                  USER_STREAMS(UI), USER_SINES(UI), &
                  PHIANG, COSPHI, &
                  EXACTDB_BRDFUNC(UI,K) )
            ENDDO
          ENDDO
        ELSE
          DO K = 1, N_USER_RELAZMS
            PHIANG = USER_RELAZMS(K)
            COSPHI = DCOS(PHIANG*DEG_TO_RAD)
            DO UI = 1, N_USER_STREAMS
              CALL LRRS_BRDF_FUNCTION &
                ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_PARS, &
                  SZASURCOS, SZASURSIN, USER_STREAMS(UI),     &
                  USER_SINES(UI),       PHIANG, COSPHI,       &
                  EXACTDB_BRDFUNC(UI,K) )
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Exit for SSFULL calculation

      IF ( DO_SSFULL ) RETURN

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam

      DO I = 1, NSTREAMS
        DO K = 1, NSTREAMS_BRDF
          CALL LRRS_BRDF_FUNCTION &
              ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_PARS, &
                SZASURCOS, SZASURSIN, &
                QUAD_STREAMS(I), QUAD_SINES(I), &
                X_BRDF(K), CX_BRDF(K), &
                BRDFUNC_0(I,K) )
        ENDDO
      ENDDO

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL LRRS_BRDF_FUNCTION &
                ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_PARS, &
                  QUAD_STREAMS(J), QUAD_SINES(J), &
                  QUAD_STREAMS(I), QUAD_SINES(I), &
                  X_BRDF(K), CX_BRDF(K), &
                  BRDFUNC(I,J,K) )
          ENDDO
        ENDDO
      ENDDO

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam. Not required if DB_CORRECTION on

        IF ( .NOT. DO_DB_CORRECTION ) THEN
          DO UI = 1, N_USER_STREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL LRRS_BRDF_FUNCTION &
                ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_PARS, &
                  SZASURCOS, SZASURSIN, &
                  USER_STREAMS(UI), USER_SINES(UI), &
                  X_BRDF(K), CX_BRDF(K), &
                  USER_BRDFUNC_0(UI,K) )
            ENDDO
          ENDDO
        ENDIF

!  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL LRRS_BRDF_FUNCTION &
                  ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_PARS, &
                    QUAD_STREAMS(J), QUAD_SINES(J), &
                    USER_STREAMS(UI), USER_SINES(UI), &
                    X_BRDF(K), CX_BRDF(K), &
                    USER_BRDFUNC(UI,J,K) )
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  debug

!      DO I = 1, NSTREAMS
!        write(86,'(i3,1p10e11.3)')
!     &        I,BRDFUNC_0(I,1),USER_BRDFUNC(1,I,2),
!     &          (BRDFUNC(I,J,1),J = 1, NSTREAMS)
!      ENDDO
!      pause

!  Finish

      RETURN
      END SUBROUTINE LRRS_BRDF_MAKER

!

      SUBROUTINE LRRS_BRDF_FUNCTION     &
           ( MAXPARS, WHICH_BRDF, PARS, &
             XJ, SXJ, XI, SXI, PHI, CPHI, KERNEL )

!  @@@ RobFix 5/5/11. New subroutine

!  module, dimensions and numbers

      USE LRRS_PARS, only :  RAHMAN_IDX, HAPKE_IDX, COXMUNK_IDX

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS
      INTEGER  , intent(in)  :: WHICH_BRDF
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI
      REAL(fpk), intent(out) :: KERNEL

!  Trawl through

      IF ( WHICH_BRDF .EQ. RAHMAN_IDX ) THEN
        CALL RAHMAN_FUNCTION ( MAXPARS, PARS, &
                               XJ, SXJ, XI, SXI, PHI, CPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. HAPKE_IDX ) THEN
        CALL HAPKE_FUNCTION  ( MAXPARS, PARS, &
                               XJ, SXJ, XI, SXI, PHI, CPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. COXMUNK_IDX ) THEN
        CALL COXMUNK_FUNCTION ( MAXPARS, PARS, &
                               XJ, SXJ, XI, SXI, PHI, CPHI, KERNEL )
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LRRS_BRDF_FUNCTION

!

      SUBROUTINE LRRS_BRDF_FOURIER &
          ( DO_INCLUDE_SURFACE, DO_UPWELLING, DO_USER_STREAMS, &
            DO_SSCORR, NSTREAMS, N_USER_STREAMS, QUAD_STRMWGT, &
            NSTREAMS_BRDF, X_BRDF, A_BRDF, FOURIER, DELFAC, &
            BRDFUNC, BRDFUNC_0, USER_BRDFUNC, USER_BRDFUNC_0, &
            BIREFLEC, BIREFLEC_0, USER_BIREFLEC, USER_BIREFLEC_0 )

!  Prepares Fourier components of the bidirectional reflectance function

!  include file of dimensions and numbers

      USE LRRS_PARS

      IMPLICIT NONE

!  Input arguments (flags, numbers and BRDFs)
!  ==========================================

!  Flags

      LOGICAL, INTENT(IN) ::          DO_UPWELLING
      LOGICAL, INTENT(IN) ::          DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::          DO_SSCORR

!  Fourier values

      INTEGER, INTENT(IN) ::          FOURIER
      REAL(FPK), INTENT(IN) :: DELFAC

!  Geometrical variables

      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          NSTREAMS_BRDF
      INTEGER, INTENT(IN) ::          N_USER_STREAMS

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGT ( MAX_STREAMS )

      REAL(FPK), INTENT(IN) :: X_BRDF ( MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(IN) :: A_BRDF ( MAX_STREAMS_BRDF )

!  at quadrature (discrete ordinate) angles

      REAL(FPK), INTENT(IN) :: BRDFUNC &
            ( MAX_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK), INTENT(IN) :: BRDFUNC_0 &
            ( MAX_STREAMS, MAX_STREAMS_BRDF )

!  at user-defined stream directions

      REAL(FPK), INTENT(IN) :: USER_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK), INTENT(IN) :: USER_BRDFUNC_0 &
            ( MAX_USER_STREAMS, MAX_STREAMS_BRDF )

!  Output arguments (BRDF Fourier components)
!  ==========================================

!  at quadrature (discrete ordinate) angles

      REAL(FPK), INTENT(OUT) :: BIREFLEC   ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: BIREFLEC_0 ( MAX_STREAMS  )

!  at user-defined stream directions

      REAL(FPK), INTENT(OUT) :: USER_BIREFLEC ( MAX_USER_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: USER_BIREFLEC_0 ( MAX_USER_STREAMS  )

!  local variables
!  ===============

!  Help

      LOGICAL ::          DO_DB_CORRECTION
      INTEGER ::          I, UI, J, K
      REAL(FPK) :: SUM, HELP, HELP_A

!  BRDF Azimuth factors

      REAL(FPK) :: BRDF_AZMFAC(MAX_STREAMS_BRDF)

!  Spherical albedo

      REAL(FPK) :: SPHERICAL_ALBEDO

!  Code
!  ====

!  Local Do Db correction flag

!       DO_DB_CORRECTION = DO_UPWELLING .AND.
!     &    ( DO_SSCORR_OUTGOING .OR. DO_SSCORR_NADIR )
       DO_DB_CORRECTION = DO_UPWELLING .AND. DO_SSCORR

!  Weighted azimuth factor
!  ( Results are stored in commons )

      IF ( FOURIER .NE. 0 ) THEN
        DO K = 1, NSTREAMS_BRDF
          BRDF_AZMFAC(K) = A_BRDF(K) * DCOS ( FOURIER * X_BRDF(K) )
        ENDDO
      ELSE
        DO K = 1, NSTREAMS_BRDF
          BRDF_AZMFAC(K) = A_BRDF(K)
        ENDDO
      ENDIF

!  surface factor

      HELP = HALF * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)

      DO I = 1, NSTREAMS
        SUM = ZERO
        DO K = 1, NSTREAMS_BRDF
          SUM  = SUM + BRDFUNC_0(I,K)*BRDF_AZMFAC(K)
        ENDDO
        BIREFLEC_0(I) = SUM * HELP
      ENDDO

!  incident quadrature directions (surface multiple reflections)

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            SUM = ZERO
            DO K = 1, NSTREAMS_BRDF
              SUM  = SUM + BRDFUNC(I,J,K) * BRDF_AZMFAC(K)
            ENDDO
             BIREFLEC(I,J) = SUM * HELP
          ENDDO
        ENDDO
      ENDIF

!  debug information

!      IF ( DO_DEBUG_WRITE ) THEN
!        WRITE(555,'(A)')'BRDF_1 Fourier 0 quad values'
!        IF ( FOURIER .EQ. 0 ) THEN
!          DO I = 1, NSTREAMS
!          WRITE(555,'(1PE12.5,3x,1P10E12.5)')
!     &     BIREFLEC_0(I,1),(BIREFLEC(I,J),J=1,NSTREAMS)
!         ENDDO
!        ENDIF
!      ENDIF

!  albedo check, always calculate the spherical albedo.
!   (Plane albedo calculations are commented out)

      IF ( FOURIER .EQ. 0 ) THEN
!  ............................... Plane albedo calculations
!          SUM = ZERO
!          DO I = 1, NSTREAMS
!            SUM = SUM + BIREFLEC_0(I) * QUAD_STRMWGHT(I)
!          ENDDO
!          SUM = SUM * TWO
!          write(*,*)0,BRDF_FACTOR, sum
!          DO J = 1, NSTREAMS
!            SUM = ZERO
!            DO I = 1, NSTREAMS
!              SUM = SUM + BIREFLEC(I,J) * QUAD_STRMWGHT(I)
!            ENDDO
!            SUM = SUM * TWO
!            write(*,*)j,sum
!          ENDDO
!  ............................... Spherical albedo calculation
        HELP_A = ZERO
        DO I = 1, NSTREAMS
          SUM = ZERO
          DO J = 1, NSTREAMS
            SUM = SUM + BIREFLEC(I,J) * QUAD_STRMWGT(J)
          ENDDO
          HELP_A = HELP_A + SUM * QUAD_STRMWGT(I)
        ENDDO
        SPHERICAL_ALBEDO = HELP_A*FOUR
!        write(*,*)(BRDF_PARS(1)-0.003d0)/0.00512d0,spherical_albedo
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam (direct beam reflection).
!        Not needed if DB_CORRECTION is in operation

        IF ( .NOT. DO_DB_CORRECTION ) THEN
          DO UI = 1, N_USER_STREAMS
            SUM = ZERO
            DO K = 1, NSTREAMS_BRDF
              SUM = SUM + USER_BRDFUNC_0(UI,K)*BRDF_AZMFAC(K)
            ENDDO
            USER_BIREFLEC_0(UI) = SUM * HELP
          ENDDO
        ENDIF

!  incident quadrature directions (surface multiple reflections)

        IF ( DO_INCLUDE_SURFACE ) THEN
          DO UI = 1, N_USER_STREAMS
            DO J = 1, NSTREAMS
              SUM = ZERO
              DO K = 1, NSTREAMS_BRDF
                SUM = SUM + USER_BRDFUNC(UI,J,K)*BRDF_AZMFAC(K)
              ENDDO
              USER_BIREFLEC(UI,J) = SUM * HELP
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LRRS_BRDF_FOURIER

!

      SUBROUTINE HAPKE_FUNCTION &
             ( MAXPARS, PARS, &
               XJ, SXJ, XI, SXI, PHI, CPHI, &
               HAPKE_KERNEL )

!  include file of constants

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER, INTENT(IN) ::          MAXPARS
      REAL(FPK), INTENT(IN) :: PARS ( MAXPARS )
      REAL(FPK), INTENT(IN) :: XI, SXI, XJ, SXJ, PHI, CPHI
      REAL(FPK), INTENT(OUT) :: HAPKE_KERNEL

!  Hapke Kernel function.
!    - New version, Fresh Coding
!    - Old version uses DISORT code; for validation.

!  input variables:

!    XI, SXI  : Cosine/Sine of angle of reflection (positive)
!    XJ, SXJ  : Cosine/Sine of angle of incidence (positive)
!    XPHI     : Difference of azimuth angles of incidence and reflection
!    PARS(1)  : single scattering albedo in Hapke's BDR model
!    PARS(2)  : angular width parameter of opposition effect in Hapke's
!    PARS(3)  : Empirical hot spot multiplier

!  local variables
!    B0_EMPIR : empirical factor to account for the finite size of
!               particles in Hapke's BDR model
!    B_HOT    : term that accounts for the opposition effect
!               (retroreflectance, hot spot) in Hapke's BDR model
!    CTHETA   : cosine of phase angle in Hapke's BDR model
!    GAMMA    : albedo factor in Hapke's BDR model
!    PHASE    : scattering phase function in Hapke's BDR model
!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model

!  local variables

      REAL(FPK) :: CTHETA, THETA, PHASE
      REAL(FPK) :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      REAL(FPK) :: SSALBEDO, GAMMA, REFLEC, FUNC
      REAL(FPK) :: HELP_J, TERM_J, HELP_I, TERM_I
      REAL(FPK) :: XPHI, CKPHI

!  Initialise

      HAPKE_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!       CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      IF ( CTHETA .GT. ONE ) CTHETA = ONE
      THETA  = DACOS( CTHETA )
      PHASE  = ONE + HALF * CTHETA

!  hot spot parameterization

      HOTSPOT  = PARS(2)
      B0_EMPIR = PARS(3)
      HELP_HOT = HOTSPOT + DTAN ( HALF * THETA )
      B_HOT    = B0_EMPIR * HOTSPOT / HELP_HOT

!  Albedo parameterization

      SSALBEDO = PARS(1)
      GAMMA    = DSQRT ( ONE - SSALBEDO )
      HELP_J   = TWO * XJ
      TERM_J   = ( ONE + HELP_J ) / ( ONE + HELP_J * GAMMA )
      HELP_I   = TWO * XI
      TERM_I   = ( ONE + HELP_I ) / ( ONE + HELP_I * GAMMA )

!  Function

      REFLEC   = SSALBEDO * QUARTER / ( XI + XJ )
      FUNC     = ( ONE + B_HOT ) * PHASE + TERM_J * TERM_I - ONE
      HAPKE_KERNEL = REFLEC * FUNC

!  Finish

      RETURN
      END SUBROUTINE HAPKE_FUNCTION

!

      SUBROUTINE RAHMAN_FUNCTION &
             ( MAXPARS, PARS, &
               XJ, SXJ, XI, SXI, PHI, CPHI, &
               RAHMAN_KERNEL )

!  Revision. 24 October 2007.
!  --------------------------

!    * Limiting cases and hotspot evaluation.
!    * Revision based on the DISORT_2 code
!    *  In Disort, this kernel is known as the RPV^ BRDF.

!     The RPV reference is:
!       Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere
!       Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable
!       With NOAA Advanced Very High Resolution Radiometer Data,
!       J. Geophys. Res., 98, 20791-20801.

!  The hotspot should occur when XI = XJ and PHI = 180.

!  include file of constants

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER, INTENT(IN) ::          MAXPARS
      REAL(FPK), INTENT(IN) :: PARS ( MAXPARS )
      REAL(FPK), INTENT(IN) :: XI, SXI, XJ, SXJ, PHI, CPHI
      REAL(FPK), INTENT(OUT) :: RAHMAN_KERNEL

!  local variables

      REAL(FPK) :: T_INC, T_REF, DT1, DT2
      REAL(FPK) :: CXI, DELTA, K1_SQ, FACT
      REAL(FPK) :: GEOM, PHASE, RFAC, K0, K1, K2
      REAL(FPK) :: XPHI, CKPHI, HSPOT, UPPER_LIMIT

      REAL(FPK), PARAMETER :: SMALL = 1.0d-04

!  Initial section
!  ---------------

!  Initialise output

      RAHMAN_KERNEL = ZERO

!  Limiting case, formerly
!      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

!  Limiting case, revised

      IF ( XJ.LT.SMALL ) RETURN

!  Azimuth convettion

      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  Hot Spot
!  --------

!  Value of hot spot

      FACT = K0 * ( TWO - K0 )
      FACT = FACT * ( ONE - K1 ) / ( ONE + K1 ) / ( ONE + K1 )
      GEOM = ( TWO * XJ * XJ * XJ ) ** ( K2 - ONE )
      HSPOT = FACT * GEOM

!  Upper limit ( 5 times hotspot value ). Follwing comments inserted.
!     This function needs more checking; some constraints are
!     required to avoid albedos larger than 1; in particular,
!     the BDREF is limited to 5 times the hotspot value to
!     avoid extremely large values at low polar angles

      UPPER_LIMIT = 5.0d0 * HSPOT

!  hot spot value

      IF ( DABS(PHI) .LT. SMALL .AND. XI.EQ.XJ ) THEN
        RAHMAN_KERNEL = HSPOT
        RETURN
      ENDIF

!  Use upper limit value at edges (low incidence or reflection)

      IF ( XI.LT.SMALL .OR. XJ.LT.SMALL ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
        RETURN
      ENDIF

!  Main section
!  ------------

!  geometrical angle xi

      CXI = XI * XJ + SXI * SXJ * CKPHI
      IF ( CXI .GT. ONE ) CXI = ONE

!  Phase function

      K1_SQ = K1 * K1
      FACT  = ( ONE + K1_SQ + TWO * K1 * CXI ) ** ONEP5
      PHASE = ( ONE - K1_SQ ) / FACT

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      RFAC = ( ONE - K0 ) / ( ONE + DELTA )

!  Geom factor and kernel

      GEOM = ( XI * XJ * ( XI + XJ ) ) ** ( K2 - ONE)
      RAHMAN_KERNEL = K0 * PHASE * ( ONE + RFAC ) * GEOM

!  Check upper limit not exceeded

      IF ( RAHMAN_KERNEL .GT. UPPER_LIMIT ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE RAHMAN_FUNCTION

!

      SUBROUTINE COXMUNK_FUNCTION &
             ( MAXPARS, PARS, &
               XJ, SXJ, XI, SXI, PHI, CPHI, &
               COXMUNK_KERNEL )

!  include file of constants

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER, INTENT(IN) ::           MAXPARS
      REAL(FPK), INTENT(IN)  :: PARS ( MAXPARS )
      REAL(FPK), INTENT(IN)  :: XI, SXI, XJ, SXJ, PHI, CPHI
      REAL(FPK), INTENT(OUT) :: COXMUNK_KERNEL

!  Critical exponent taken out

      REAL(FPK), PARAMETER  :: CRITEXP = 88.0D0

!  Local variables

      REAL(FPK) :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      REAL(FPK) :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      REAL(FPK) :: XPHI, CKPHI
      REAL(FPK) :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(FPK) :: SHADOWI, SHADOWR, SHADOW

!  Shadow variables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      COXMUNK_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Kernel

!  ..Scatter angles

! old   Z = - XI * XJ + SXI * SXJ * CKPHI
! old   IF ( Z .LT. MINUS_ONE) Z = MINUS_ONE
! old   Z1 = DACOS(-Z)
! old   Z2 = DCOS(Z1*HALF)

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

!  .. Fresnel coefficients

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(2) * Z2
      H2 = DSQRT ( PARS(2) + Z2_SQ_M1 )
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      XMP = HALF * ( RP*RP + RL*RL )

!  Coxmunk Function

      A = TWO * Z2
      B = ( XI + XJ ) / A
      IF ( B .GT. ONE ) B = ONE
      A = PIO2 - DASIN(B)
      TA = DTAN(A)
      ARGUMENT = TA * TA  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB  / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
        COXMUNK_KERNEL = XMP * FAC1 * FAC2 / XJ
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT = XI/DSQRT(ONE-XXI)
      T1   = DEXP(-DCOT*DCOT*S2)
!      T2   = DERFC(DCOT*S3)
      T2   = DERFC_E(DCOT*S3)
      SHADOWI = HALF*(S1*T1/DCOT-T2)

      XXJ  = XJ*XJ
      DCOT = XJ/DSQRT(ONE-XXJ)
      T1   = DEXP(-DCOT*DCOT*S2)
!      T2   = DERFC(DCOT*S3)
      T2   = DERFC_E(DCOT*S3)
      SHADOWR = HALF*(S1*T1/DCOT-T2)

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW

!     Finish

      RETURN
      END SUBROUTINE COXMUNK_FUNCTION

!

      SUBROUTINE COXMUNK_FUNCTION_DB &
             ( MAXPARS, PARS, &
               XJ, SXJ, XI, SXI, PHI, CPHI, &
               COXMUNK_KERNEL )

!  include file of constants

      USE LRRS_PARS
      USE LRRS_AUX2

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER, INTENT(IN) ::          MAXPARS
      REAL(FPK), INTENT(IN) :: PARS ( MAXPARS )
      REAL(FPK), INTENT(IN) :: XI, SXI, XJ, SXJ, PHI, CPHI
      REAL(FPK), INTENT(OUT) :: COXMUNK_KERNEL

!  local variables
!  ---------------

!  help variables

      INTEGER ::          n, k, i, i1, N_phiquad_HALF
      REAL(FPK) :: XM, SXM, sum_pr
      REAL(FPK) :: sumr, w_p, reflec_0, reflec_1
      REAL(FPK) :: phi_sub1, cphi_sub1
      REAL(FPK) :: phi_sub2, cphi_sub2

!  Local quadrature stuff

      INTEGER, PARAMETER :: max_msrs_muquad  = 40, &
                            max_msrs_phiquad = 50
      INTEGER ::            n_muquad, n_phiquad

!  arrays

      REAL(FPK) :: X_MUQUAD (max_msrs_muquad)
      REAL(FPK) :: W_MUQUAD (max_msrs_muquad)
      REAL(FPK) :: SX_MUQUAD (max_msrs_muquad)
      REAL(FPK) :: WXX_MUQUAD(max_msrs_muquad)

      REAL(FPK) :: X_PHIQUAD (max_msrs_phiquad)
      REAL(FPK) :: W_PHIQUAD (max_msrs_phiquad)

      REAL(FPK) :: R0_QUAD_IN  (max_msrs_muquad,max_msrs_phiquad)
      REAL(FPK) :: R0_OUT_QUAD (max_msrs_muquad,max_msrs_phiquad)

!  Safety first zeroing

      REFLEC_0 = ZERO
      REFLEC_1 = ZERO
      COXMUNK_KERNEL = ZERO

!  Air to water, Polar quadrature

      n_muquad = 20
      CALL GAULEG ( ZERO, ONE, X_muquad, W_muquad, n_muquad )
      DO I = 1, N_MUQUAD
        XM = X_MUQUAD(I)
        SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
        WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
      ENDDO

!  Azimuth quadrature

      n_phiquad = 40
      N_phiquad_HALF = N_PHIQUAD / 2
      CALL GAULEG ( ZERO, ONE, X_PHIQUAD, W_PHIQUAD, N_PHIQUAD_HALF )
      DO I = 1, N_PHIQUAD_HALF
        I1 = I + N_PHIQUAD_HALF
        X_PHIQUAD(I1) = - X_PHIQUAD(I)
        W_PHIQUAD(I1) =   W_PHIQUAD(I)
      ENDDO
      DO I = 1, N_PHIQUAD
        X_PHIQUAD(I)  = PIE * X_PHIQUAD(I)
      ENDDO

!  Single scattering (zero order), Phi is in degrees here!

      CALL COXMUNK_FUNCTION &
             ( MAXPARS, PARS, &
               XJ, SXJ, XI, SXI, PHI, CPHI, &
               REFLEC_0 )

!  Quadrature output for first order R/T calculations

      DO K = 1, n_muquad
        XM  = X_MUQUAD(K)
        SXM = SX_MUQUAD(K)
        DO N = 1, N_PHIQUAD
          PHI_SUB1  = X_PHIQUAD(N)
          CPHI_SUB1 = DCOS(PHI_SUB1)
          PHI_SUB2  = PHI*DEG_TO_RAD - X_PHIQUAD(N)
          CPHI_SUB2 = DCOS(PHI_SUB2)
          CALL COXMUNK_FUNCTION &
             ( MAXPARS, PARS, &
               XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, &
               R0_OUT_QUAD(K,N) )
          CALL COXMUNK_FUNCTION &
             ( MAXPARS, PARS, &
               XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, &
               R0_QUAD_IN(K,N) )
        ENDDO
      ENDDO

!  compute the next order

      SUMR = ZERO
      DO K = 1, n_muquad
        SUM_PR = ZERO
        DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          SUM_PR = SUM_PR + W_P * R0_QUAD_IN(K,N) * R0_OUT_QUAD(K,N)
        ENDDO
        SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
      ENDDO
      REFLEC_1 = SUMR

!  Compute total

      COXMUNK_KERNEL = REFLEC_0 + REFLEC_1

!      write(34,'(1p6e14.5)')reflec_0, reflec_1, coxmunk_kernel,
!     &   dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_FUNCTION_DB

!

      FUNCTION DERFC_E(X)

      REAL(FPK), INTENT(IN) :: X
      REAL(FPK) :: DERFC_E

! Returns the complementary error function erfc(x) with fractional error
! everywhere less than 1.2 * 10^7.

      REAL(FPK) :: T,Z

      Z = DABS(X)
      T = 1.D0/(1.D0+0.5D0*Z)
      DERFC_E = T*DEXP(-Z*Z-1.26551223D0+T*(1.00002368D0+T*(.37409196D0+ &
              T*(.09678418D0+T*(-.18628806D0+T*(.27886807D0+T* &
              (-1.13520398D0+T*(1.48851587D0+T*(-.82215223D0+T* &
              .17087277D0)))))))))
      IF (X .LT. 0.D0) DERFC_E = 2.D0-DERFC_E

      RETURN
      END FUNCTION DERFC_E

      END MODULE lrrs_brdf

