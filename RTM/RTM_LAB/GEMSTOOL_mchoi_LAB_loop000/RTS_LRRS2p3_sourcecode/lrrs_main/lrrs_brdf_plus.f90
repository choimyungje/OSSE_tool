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
! #            LRRS_BRDF_MASTER_PLUS (master), calling          #
! #              LRRS_BRDF_MAKER_PLUS                           #
! #              LRRS_BRDF_FUNCTION_PLUS                        #
! #                                                             #
! #            LRRS_BRDF_FOURIER_PLUS (called by FOURIER_PLUS)  #
! #                                                             #
! # Choice of BRDFs:                                            #
! #                                                             #
! #            HAPKE_FUNCTION_PLUS                              #
! #            RAHMAN_FUNCTION_PLUS                             #
! #            COXMUNK_FUNCTION_PLUS                            #
! #            COXMUNK_FUNCTION_DB_PLUS                         #
! #                                                             #
! #   Note (*): Version 2.1 with BRDF setup, November 2007.     #
! #                                                             #
! ###############################################################


      MODULE lrrs_brdf_plus

      USE LRRS_PARS

      PRIVATE
      PUBLIC :: LRRS_BRDF_MASTER_PLUS ,  &
                LRRS_BRDF_MAKER_PLUS,    &
                LRRS_BRDF_FUNCTION_PLUS, &
                LRRS_BRDF_FOURIER_PLUS,  &
                HAPKE_FUNCTION_PLUS,     &
                RAHMAN_FUNCTION_PLUS,    &
                COXMUNK_FUNCTION_PLUS,   &
                COXMUNK_FUNCTION_PLUS_DB

!  @@@ RobFix 5/5/11. Drop the External calls.
!    Need new routine, LRRS_BRDF_FUNCTION_PLUS, modeled after LIDORT.

      CONTAINS

      SUBROUTINE LRRS_BRDF_MASTER_PLUS &
          ( WHICH_BRDF, DO_GLITTER_DBMS, BRDF_PARS, &
            BRDFPAR_DERIV_INDEX, &
            DO_UPWELLING, DO_USER_STREAMS, &
            DO_LAMBERT_VARIATION, &
            DO_SSCORR_OUTGOING, DO_SSFULL, &
            NSTREAMS_BRDF,  SOLAR_ANGLE, &
            NSTREAMS,       QUAD_STREAMS, &
            N_USER_STREAMS, USER_STREAMS, &
            N_USER_RELAZMS, USER_RELAZMS, &
            X_BRDF, A_BRDF, &
            EXACTDB_BRDFUNC, LS_EXACTDB_BRDFUNC, &
            BRDFUNC,         LS_BRDFUNC, &
            BRDFUNC_0,       LS_BRDFUNC_0, &
            USER_BRDFUNC,    LS_USER_BRDFUNC, &
            USER_BRDFUNC_0,  LS_USER_BRDFUNC_0 )

!  Master routine for creation of BRDFs
!  @@@ RobFix 5/5/11. Drop the External calls

!  include file of dimensions and numbers

      USE LRRS_PARS
      USE LRRS_AUX2
      USE LRRS_BRDF

      IMPLICIT NONE

!  Input arguments
!  ===============

!  BRDF parameters

      INTEGER, INTENT(IN) ::          WHICH_BRDF
      LOGICAL, INTENT(IN) ::          DO_GLITTER_DBMS
      REAL(FPK), INTENT(IN) ::        BRDF_PARS(3)
      INTEGER, INTENT(IN) ::          BRDFPAR_DERIV_INDEX

!  Flags

      LOGICAL, INTENT(IN) ::          DO_UPWELLING
      LOGICAL, INTENT(IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT(IN) ::          DO_SSFULL

!  Lambertian variation flag

      LOGICAL, INTENT(IN) ::          DO_LAMBERT_VARIATION

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

!  Linearizations

      REAL(FPK), INTENT(OUT) :: LS_BRDFUNC &
            ( MAX_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(OUT) :: LS_BRDFUNC_0 &
            ( MAX_STREAMS, MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(OUT) :: LS_USER_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(OUT) :: LS_USER_BRDFUNC_0 &
            ( MAX_USER_STREAMS, MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(OUT) :: LS_EXACTDB_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Local definitions
!  =================

!  BRDF functions
!  Single Kernel functions. No 3-kernel implementation here.

!  ordinary BRDF without derivatives

      !EXTERNAL       HAPKE_FUNCTION
      !EXTERNAL       RAHMAN_FUNCTION
      !EXTERNAL       COXMUNK_FUNCTION
      !EXTERNAL       COXMUNK_FUNCTION_DB

!  BRDFs with derivatives

      !EXTERNAL       HAPKE_FUNCTION_PLUS
      !EXTERNAL       RAHMAN_FUNCTION_PLUS
      !EXTERNAL       COXMUNK_FUNCTION_PLUS
      !EXTERNAL       COXMUNK_FUNCTION_PLUS_DB

!  Hapke old uses exact DISORT code
!      EXTERNAL       HAPKE_FUNCTION_OLD

!  other variables

      INTEGER ::          I, I1, NBRDF_HALF
      REAL(FPK) :: CX_BRDF ( MAX_STREAMS_BRDF )
      REAL(FPK) :: SX_BRDF ( MAX_STREAMS_BRDF )

!  BRDF quadrature
!  ---------------

!  Save these quantities for efficient coding

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

      IF ( .NOT. DO_LAMBERT_VARIATION ) THEN
         CALL LRRS_BRDF_MAKER_PLUS &
             ( WHICH_BRDF, BRDF_PARS, &
               DO_UPWELLING, DO_USER_STREAMS, DO_GLITTER_DBMS,     &
               BRDFPAR_DERIV_INDEX, DO_SSCORR_OUTGOING, DO_SSFULL, &
               NSTREAMS_BRDF,  X_BRDF, CX_BRDF, SOLAR_ANGLE, &
               NSTREAMS,       QUAD_STREAMS, &
               N_USER_STREAMS, USER_STREAMS, &
               N_USER_RELAZMS, USER_RELAZMS, &
               EXACTDB_BRDFUNC, LS_EXACTDB_BRDFUNC, &
               BRDFUNC,         LS_BRDFUNC, &
               BRDFUNC_0,       LS_BRDFUNC_0, &
               USER_BRDFUNC,    LS_USER_BRDFUNC, &
               USER_BRDFUNC_0,  LS_USER_BRDFUNC_0 )
      ELSE
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
      ENDIF

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
      END SUBROUTINE LRRS_BRDF_MASTER_PLUS

!

      SUBROUTINE LRRS_BRDF_MAKER_PLUS &
             ( WHICH_BRDF, BRDF_PARS, &
               DO_UPWELLING, DO_USER_STREAMS, DO_GLITTER_DBMS,      &
               BRDFPAR_DERIV_INDEX,  DO_SSCORR_OUTGOING, DO_SSFULL, &
               NSTREAMS_BRDF,  X_BRDF, CX_BRDF, SOLAR_ANGLE, &
               NSTREAMS,       QUAD_STREAMS, &
               N_USER_STREAMS, USER_STREAMS, &
               N_USER_RELAZMS, USER_RELAZMS, &
               EXACTDB_BRDFUNC, LS_EXACTDB_BRDFUNC, &
               BRDFUNC,         LS_BRDFUNC, &
               BRDFUNC_0,       LS_BRDFUNC_0, &
               USER_BRDFUNC,    LS_USER_BRDFUNC, &
               USER_BRDFUNC_0,  LS_USER_BRDFUNC_0 )

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
!      EXTERNAL         BRDF_FUNCTION_PLUS
!      EXTERNAL         BRDF_FUNCTION_DB_PLUS

!  parameters and derivative index

      REAL(FPK), INTENT(IN) :: BRDF_PARS(3)
      INTEGER, INTENT(IN)   :: BRDFPAR_DERIV_INDEX

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

!  Linearizations

      REAL(FPK), INTENT(OUT) :: LS_BRDFUNC &
            ( MAX_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(OUT) :: LS_BRDFUNC_0 &
            ( MAX_STREAMS, MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(OUT) :: LS_USER_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(OUT) :: LS_USER_BRDFUNC_0 &
            ( MAX_USER_STREAMS, MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(OUT) :: LS_EXACTDB_BRDFUNC &
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
              CALL COXMUNK_FUNCTION_PLUS_DB &
                ( MAX_BRDF_PARAMETERS, &
                  BRDF_PARS, BRDFPAR_DERIV_INDEX, &
                  SZASURCOS, SZASURSIN, &
                  USER_STREAMS(UI), USER_SINES(UI), &
                  PHIANG, COSPHI, &
                  EXACTDB_BRDFUNC(UI,K), &
                  LS_EXACTDB_BRDFUNC(UI,K) )
            ENDDO
          ENDDO
        ELSE
          DO K = 1, N_USER_RELAZMS
            PHIANG = USER_RELAZMS(K)
            COSPHI = DCOS(PHIANG*DEG_TO_RAD)
            DO UI = 1, N_USER_STREAMS
              CALL LRRS_BRDF_FUNCTION_PLUS &
                ( MAX_BRDF_PARAMETERS, WHICH_BRDF, &
                  BRDF_PARS, BRDFPAR_DERIV_INDEX,  &
                  SZASURCOS, SZASURSIN, &
                  USER_STREAMS(UI), USER_SINES(UI), &
                  PHIANG, COSPHI, &
                  EXACTDB_BRDFUNC(UI,K), &
                  LS_EXACTDB_BRDFUNC(UI,K) )
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
          CALL LRRS_BRDF_FUNCTION_PLUS &
              ( MAX_BRDF_PARAMETERS, WHICH_BRDF, &
                BRDF_PARS, BRDFPAR_DERIV_INDEX, &
                SZASURCOS, SZASURSIN, &
                QUAD_STREAMS(I), QUAD_SINES(I), &
                X_BRDF(K), CX_BRDF(K), &
                BRDFUNC_0(I,K), &
                LS_BRDFUNC_0(I,K) )
        ENDDO
      ENDDO

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL LRRS_BRDF_FUNCTION_PLUS &
                ( MAX_BRDF_PARAMETERS, WHICH_BRDF, &
                  BRDF_PARS, BRDFPAR_DERIV_INDEX, &
                  QUAD_STREAMS(J), QUAD_SINES(J), &
                  QUAD_STREAMS(I), QUAD_SINES(I), &
                  X_BRDF(K), CX_BRDF(K), &
                  BRDFUNC(I,J,K), &
                  LS_BRDFUNC(I,J,K) )
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
              CALL LRRS_BRDF_FUNCTION_PLUS &
                ( MAX_BRDF_PARAMETERS, WHICH_BRDF, &
                  BRDF_PARS, BRDFPAR_DERIV_INDEX, &
                  SZASURCOS, SZASURSIN, &
                  USER_STREAMS(UI), USER_SINES(UI), &
                  X_BRDF(K), CX_BRDF(K), &
                  USER_BRDFUNC_0(UI,K), &
                  LS_USER_BRDFUNC_0(UI,K) )
            ENDDO
          ENDDO
        ENDIF

!  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL LRRS_BRDF_FUNCTION_PLUS &
                  ( MAX_BRDF_PARAMETERS, WHICH_BRDF, &
                    BRDF_PARS, BRDFPAR_DERIV_INDEX, &
                    QUAD_STREAMS(J), QUAD_SINES(J), &
                    USER_STREAMS(UI), USER_SINES(UI), &
                    X_BRDF(K), CX_BRDF(K), &
                    USER_BRDFUNC(UI,J,K), &
                    LS_USER_BRDFUNC(UI,J,K) )
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
      END SUBROUTINE LRRS_BRDF_MAKER_PLUS

!

      SUBROUTINE LRRS_BRDF_FUNCTION_PLUS           &
           ( MAXPARS, WHICH_BRDF, PARS, PAR_INDEX, &
             XJ, SXJ, XI, SXI, PHI, CPHI, KERNEL, DKERNEL )

!  @@@ RobFix 5/5/11. New subroutine

!  module, dimensions and numbers

      USE LRRS_PARS, only :  RAHMAN_IDX, HAPKE_IDX, COXMUNK_IDX

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS
      INTEGER  , intent(in)  :: WHICH_BRDF
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      INTEGER  , intent(in)  :: PAR_INDEX
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI
      REAL(fpk), intent(out) :: KERNEL
      REAL(fpk), intent(out) :: DKERNEL

!  Trawl through

      IF ( WHICH_BRDF .EQ. RAHMAN_IDX ) THEN
        CALL RAHMAN_FUNCTION_PLUS         &
             ( MAXPARS, PARS, PAR_INDEX,  &
              XJ, SXJ, XI, SXI, PHI, CPHI, KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. HAPKE_IDX ) THEN
        CALL HAPKE_FUNCTION_PLUS          &
             ( MAXPARS, PARS, PAR_INDEX,  &
              XJ, SXJ, XI, SXI, PHI, CPHI, KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. COXMUNK_IDX ) THEN
        CALL COXMUNK_FUNCTION_PLUS        &
             ( MAXPARS, PARS, PAR_INDEX,  &
              XJ, SXJ, XI, SXI, PHI, CPHI, KERNEL, DKERNEL )
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LRRS_BRDF_FUNCTION_PLUS

!

      SUBROUTINE LRRS_BRDF_FOURIER_PLUS &
          ( DO_INCLUDE_SURFACE, DO_UPWELLING, DO_USER_STREAMS, &
            DO_SSCORR, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, &
            X_BRDF, A_BRDF, FOURIER, DELFAC, &
            LS_BRDFUNC,        LS_BRDFUNC_0, &
            LS_USER_BRDFUNC,   LS_USER_BRDFUNC_0, &
            LS_BIREFLEC,       LS_BIREFLEC_0, &
            LS_USER_BIREFLEC,  LS_USER_BIREFLEC_0 )

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

      REAL(FPK), INTENT(IN) :: X_BRDF ( MAX_STREAMS_BRDF )
      REAL(FPK), INTENT(IN) :: A_BRDF ( MAX_STREAMS_BRDF )

!  at quadrature (discrete ordinate) angles

      REAL(FPK), INTENT(IN) :: LS_BRDFUNC &
            ( MAX_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK), INTENT(IN) :: LS_BRDFUNC_0 &
            ( MAX_STREAMS, MAX_STREAMS_BRDF )

!  at user-defined stream directions

      REAL(FPK), INTENT(IN) :: LS_USER_BRDFUNC &
            ( MAX_USER_STREAMS, MAX_STREAMS, MAX_STREAMS_BRDF )

      REAL(FPK), INTENT(IN) :: LS_USER_BRDFUNC_0 &
            ( MAX_USER_STREAMS, MAX_STREAMS_BRDF )

!  Output arguments (BRDF Fourier components)
!  ==========================================

!  at quadrature (discrete ordinate) angles

      REAL(FPK), INTENT(OUT) :: LS_BIREFLEC   ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: LS_BIREFLEC_0 ( MAX_STREAMS  )

!  at user-defined stream directions

      REAL(FPK), INTENT(OUT) :: LS_USER_BIREFLEC (MAX_USER_STREAMS,MAX_STREAMS)
      REAL(FPK), INTENT(OUT) :: LS_USER_BIREFLEC_0 ( MAX_USER_STREAMS  )

!  local variables
!  ===============

      LOGICAL ::             DO_DB_CORRECTION
      INTEGER ::             I, UI, J, K
      REAL(FPK) ::    SUM, HELP

!  BRDF Azimuth factors

      REAL(FPK) :: BRDF_AZMFAC(MAX_STREAMS_BRDF)

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
          SUM  = SUM + LS_BRDFUNC_0(I,K)*BRDF_AZMFAC(K)
        ENDDO
        LS_BIREFLEC_0(I) = SUM * HELP
      ENDDO

!  incident quadrature directions (surface multiple reflections)

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            SUM = ZERO
            DO K = 1, NSTREAMS_BRDF
              SUM  = SUM + LS_BRDFUNC(I,J,K) * BRDF_AZMFAC(K)
            ENDDO
             LS_BIREFLEC(I,J) = SUM * HELP
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

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam (direct beam reflection).
!        Not needed if DB_CORRECTION is in operation

        IF ( .NOT. DO_DB_CORRECTION ) THEN
          DO UI = 1, N_USER_STREAMS
            SUM = ZERO
            DO K = 1, NSTREAMS_BRDF
              SUM = SUM + LS_USER_BRDFUNC_0(UI,K)*BRDF_AZMFAC(K)
            ENDDO
            LS_USER_BIREFLEC_0(UI) = SUM * HELP
          ENDDO
        ENDIF

!  incident quadrature directions (surface multiple reflections)

        IF ( DO_INCLUDE_SURFACE ) THEN
          DO UI = 1, N_USER_STREAMS
            DO J = 1, NSTREAMS
              SUM = ZERO
              DO K = 1, NSTREAMS_BRDF
                SUM = SUM + LS_USER_BRDFUNC(UI,J,K)*BRDF_AZMFAC(K)
              ENDDO
              LS_USER_BIREFLEC(UI,J) = SUM * HELP
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LRRS_BRDF_FOURIER_PLUS

!

      SUBROUTINE HAPKE_FUNCTION_PLUS &
             ( MAXPARS, PARS, PAR_INDEX, &
               XJ, SXJ, XI, SXI, PHI, CPHI, &
               HAPKE_KERNEL, HAPKE_DERIVATIVE )

!  include file of constants amd dimensions

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER, INTENT(IN) ::          MAXPARS
      REAL(FPK), INTENT(IN) :: PARS ( MAXPARS )
      INTEGER, INTENT(IN) ::          PAR_INDEX
      REAL(FPK), INTENT(IN) :: XI, SXI, XJ, SXJ, PHI, CPHI
      REAL(FPK), INTENT(OUT) :: HAPKE_KERNEL
      REAL(FPK), INTENT(OUT) :: HAPKE_DERIVATIVE

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

!  Local variables

      REAL(FPK) :: CTHETA, THETA, PHASE
      REAL(FPK) :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      REAL(FPK) :: SSALBEDO, GAMMA, REFLEC, FUNC
      REAL(FPK) :: HELP_J, GHELP_J, TERM_J
      REAL(FPK) :: HELP_I, GHELP_I, TERM_I
      REAL(FPK) :: TI_TJ, DT1, DT2
      REAL(FPK) :: XPHI, CKPHI

!  Initialise

      HAPKE_KERNEL     = ZERO
      HAPKE_DERIVATIVE = ZERO

!  Switch azimuth

      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!        CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
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
      GHELP_J  = ( ONE + HELP_J * GAMMA )
      TERM_J   = ( ONE + HELP_J ) / GHELP_J
      HELP_I   = TWO * XI
      GHELP_I  = ( ONE + HELP_I * GAMMA )
      TERM_I   = ( ONE + HELP_I ) / GHELP_I
      TI_TJ    = TERM_J * TERM_I

!  Function

      REFLEC   = SSALBEDO * QUARTER / ( XI + XJ )
      FUNC     = ( ONE + B_HOT ) * PHASE + TI_TJ - ONE
      HAPKE_KERNEL = REFLEC * FUNC

!  ssalbedo derivative

      IF ( PAR_INDEX .EQ. 1 ) THEN
        DT1 = HAPKE_KERNEL / SSALBEDO
        DT2 = ( HELP_J / GHELP_J ) + ( HELP_I / GHELP_I )
        DT2 = DT2 * TI_TJ * HALF / GAMMA
        HAPKE_DERIVATIVE = DT1 + DT2 * REFLEC
      ENDIF

!  Hotspot  derivative

      IF ( PAR_INDEX .EQ. 2 ) THEN
        DT1 = B_HOT * ( B0_EMPIR - B_HOT ) / B0_EMPIR / HOTSPOT
        HAPKE_DERIVATIVE = DT1 * REFLEC * PHASE
      ENDIF

!  empirical factor derivative

      IF ( PAR_INDEX .EQ. 3 ) THEN
        DT1 = B_HOT / B0_EMPIR
        HAPKE_DERIVATIVE = DT1 * REFLEC * PHASE
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE HAPKE_FUNCTION_PLUS

!

      SUBROUTINE RAHMAN_FUNCTION_PLUS &
             ( MAXPARS, PARS, PAR_INDEX, &
               XJ, SXJ, XI, SXI, PHI, CPHI, &
               RAHMAN_KERNEL, RAHMAN_DERIVATIVE )

!  include file of constants amd dimensions

      USE LRRS_PARS

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER, INTENT(IN) ::          MAXPARS
      REAL(FPK), INTENT(IN) :: PARS ( MAXPARS )
      INTEGER, INTENT(IN) ::          PAR_INDEX
      REAL(FPK), INTENT(IN) :: XI, SXI, XJ, SXJ, PHI, CPHI
      REAL(FPK), INTENT(OUT) :: RAHMAN_KERNEL
      REAL(FPK), INTENT(OUT) :: RAHMAN_DERIVATIVE

!  local variables

      REAL(FPK) :: T_INC, T_REF, DT1, DT2
      REAL(FPK) :: CXI, DELTA, K1_SQ, FACT, K0, K1, K2
      REAL(FPK) :: HELPM, HELPR, HELPG, D_HELPM, D_FACT
      REAL(FPK) :: GEOM, PHASE, RFAC, RFAC1, D_K0, D_K1, D_K2
      REAL(FPK) :: XPHI, CKPHI, HSPOT, UPPER_LIMIT

      REAL(FPK), PARAMETER :: SMALL = 1.0d-10

!  Initialise

      RAHMAN_KERNEL = ZERO
      RAHMAN_DERIVATIVE = ZERO

      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  Hot Spot
!  --------

!  Value of hot spot

      FACT = K0 * ( 2.0d0 - K0 )
      FACT = FACT * ( 1.0d0 - K1 ) / ( 1.0d0 + K1 ) / ( 1.0d0 + K1 )
      GEOM = ( 2.0d0 * XJ * XJ * XJ ) ** ( K2 - 1.0d0 )
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
      HELPM = ONE - K1_SQ
      FACT  = ONE + K1_SQ + TWO * K1 * CXI
      PHASE = HELPM / ( FACT ** ONEP5 )

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      HELPR = ONE / ( ONE + DELTA )
      RFAC  = ( ONE - K0 ) * HELPR
      RFAC1 = ONE + RFAC

!  Geom factor and kernel

      HELPG = XI * XJ * ( XI + XJ )
      GEOM  = HELPG ** ( K2 - ONE)
      RAHMAN_KERNEL = K0 * PHASE * RFAC1 * GEOM

!  K0 derivative

      IF ( PAR_INDEX .EQ. 1 ) THEN
        D_K0   = ( ONE / K0 ) - ( HELPR / RFAC1 )
        RAHMAN_DERIVATIVE = RAHMAN_KERNEL * D_K0
      ENDIF

!  Phase function derivative

      IF ( PAR_INDEX .EQ. 2 ) THEN
        D_FACT  =   TWO * K1 + TWO * CXI
        D_HELPM = - TWO * K1
        D_K1    = ( D_HELPM / HELPM ) - ONEP5 * ( D_FACT / FACT )
        RAHMAN_DERIVATIVE = RAHMAN_KERNEL * D_K1
      ENDIF

!  K2 derivative

      IF ( PAR_INDEX .EQ. 3 ) THEN
        D_K2 = DLOG ( HELPG )
        RAHMAN_DERIVATIVE = RAHMAN_KERNEL * D_K2
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE RAHMAN_FUNCTION_PLUS

!

      SUBROUTINE COXMUNK_FUNCTION_PLUS &
             ( MAXPARS, PARS, PAR_INDEX, &
               XJ, SXJ, XI, SXI, PHI, CPHI, &
               COXMUNK_KERNEL, COXMUNK_DERIVATIVE )

!  include file of constants amd dimensions

      USE LRRS_PARS
      USE LRRS_BRDF

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER, INTENT(IN) ::           MAXPARS
      REAL(FPK), INTENT(IN)  :: PARS ( MAXPARS )
      INTEGER, INTENT(IN) ::           PAR_INDEX
      REAL(FPK), INTENT(IN)  :: XI, SXI, XJ, SXJ, PHI, CPHI
      REAL(FPK), INTENT(OUT) :: COXMUNK_KERNEL
      REAL(FPK), INTENT(OUT) :: COXMUNK_DERIVATIVE

!  Critical exponent taken out

      REAL(FPK), PARAMETER  :: CRITEXP = 88.0D0

!  Local variables

      REAL(FPK) :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      REAL(FPK) :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      REAL(FPK) :: S1, S2, S3, XXI, XXJ
      REAL(FPK) :: T1_I, T2_I, DCOT_I
      REAL(FPK) :: T1_R, T2_R, DCOT_R
      REAL(FPK) :: SHADOWI, SHADOWR, SHADOW
      REAL(FPK) :: XPHI, CKPHI

      REAL(FPK) :: H1H2, H2Z2, TA_SQ, DFAC2, DH1, DH2, DRP, DRL
      REAL(FPK) :: D_S1, D_S2, D_T1, D_T2
      REAL(FPK) :: D_SHADOWI, D_SHADOWR, D_SHADOW

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      COXMUNK_KERNEL     = ZERO
      COXMUNK_DERIVATIVE = ZERO

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
      H1H2 = H1 + H2
      RP = ( H1 - H2 ) / H1H2
      H2Z2 = Z2 + H2
      RL = ( Z2 - H2 ) / H2Z2
      XMP = HALF * ( RP*RP + RL*RL )

!  Coxmunk Function

      A = TWO * Z2
      B = ( XI + XJ ) / A
      IF ( B .GT. ONE ) B = ONE
      A  = PIO2 - DASIN(B)
      TA = DTAN(A)
      TA_SQ = TA * TA
      ARGUMENT = TA_SQ  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
        COXMUNK_KERNEL = XMP * FAC1 * FAC2 / XJ
      ENDIF

!  inverse slope-squared derivative

      IF ( PAR_INDEX .EQ. 1 ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DFAC2 = ( PARS(1) - TA_SQ ) / PARS(1) / PARS(1)
          COXMUNK_DERIVATIVE = - COXMUNK_KERNEL * DFAC2
        ENDIF
      ENDIF

!  square refractive index derivative
!  --This section of code was formerly at the end of routine
!    -- Now moved here before the shadowing option
!    -- otherwise derivative will not get done
!         Bug found by V. Natraj in VLIDORT. 02 February 2007.

      IF ( PAR_INDEX .EQ. 2 ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DH1 = Z2
          DH2 = HALF / H2
          DRP = ( DH1 * ( ONE - RP ) - DH2 * ( ONE + RP ) ) / H1H2
          DRL =  - DH2 * ( ONE + RL ) / H2Z2
          DFAC2 = ( RP*DRP + RL*DRL ) / XMP
          COXMUNK_DERIVATIVE = COXMUNK_KERNEL * DFAC2
        ENDIF
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code
!  -----------

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT_I = XI/DSQRT(ONE-XXI)
      T1_I   = DEXP(-DCOT_I*DCOT_I*S2)
!      T2_I   = DERFC(DCOT_I*S3)
      T2_I   = DERFC_E(DCOT_I*S3)
      SHADOWI = HALF * ( S1*T1_I/DCOT_I - T2_I )

      XXJ  = XJ*XJ
      DCOT_R = XJ/DSQRT(ONE-XXJ)
      T1_R   = DEXP(-DCOT_R*DCOT_R*S2)
!      T2_R   = DERFC(DCOT_R*S3)
      T2_R   = DERFC_E(DCOT_R*S3)
      SHADOWR = HALF * ( S1*T1_R/DCOT_R - T2_R )

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW

!  Update Scalar derivatives
!  -------------------------

!  add the shadow derivative to inverse slope-squared derivative

      IF ( PAR_INDEX .EQ. 1 ) THEN
        D_S1 = HALF / PIE / S1
        D_S2 = - S2 * S2
        D_T1 = - T1_I * DCOT_I * DCOT_I * D_S2
        D_T2 = S2 * S2 * DCOT_I * S1 * T1_I
        D_SHADOWI = HALF * ( D_S1*T1_I/DCOT_I + S1*D_T1/DCOT_I - D_T2 )
        D_T1 = - T1_R * DCOT_R * DCOT_R * D_S2
        D_T2 = S2 * S2 * DCOT_R * S1 * T1_R
        D_SHADOWR = HALF * ( D_S1*T1_R/DCOT_R + S1*D_T1/DCOT_R - D_T2 )
        D_SHADOW = - SHADOW * SHADOW * ( D_SHADOWI + D_SHADOWR )
        COXMUNK_DERIVATIVE = &
               COXMUNK_DERIVATIVE * SHADOW+ &
               COXMUNK_KERNEL * D_SHADOW / SHADOW
      ENDIF

!  Refractive index derivative, update

      IF ( PAR_INDEX .EQ. 2 ) THEN
        COXMUNK_DERIVATIVE = COXMUNK_DERIVATIVE * SHADOW
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_FUNCTION_PLUS

!

      SUBROUTINE COXMUNK_FUNCTION_PLUS_DB &
             ( MAXPARS, PARS, PAR_INDEX, &
               XJ, SXJ, XI, SXI, PHI, CPHI, &
               COXMUNK_KERNEL, COXMUNK_DERIVATIVE )

!  include file of constants

      USE LRRS_PARS
      USE LRRS_AUX2

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER, INTENT(IN) ::          MAXPARS
      REAL(FPK), INTENT(IN) :: PARS ( MAXPARS )
      INTEGER, INTENT(IN) ::          PAR_INDEX
      REAL(FPK), INTENT(IN) :: XI, SXI, XJ, SXJ, PHI, CPHI
      REAL(FPK), INTENT(OUT) :: COXMUNK_KERNEL
      REAL(FPK), INTENT(OUT) :: COXMUNK_DERIVATIVE

!  local variables
!  ---------------

!  help variables

      INTEGER ::          n, k, i, i1, N_phiquad_HALF
      REAL(FPK) :: XM, SXM, sum_pr, pr
      REAL(FPK) :: sumr, w_p, reflec_0, reflec_1
      REAL(FPK) :: d_reflec_0, d_reflec_1
      REAL(FPK) :: phi_sub1, cphi_sub1
      REAL(FPK) :: phi_sub2, cphi_sub2

!  Local quadrature stuff

      INTEGER, PARAMETER :: max_msrs_muquad= 40, &
                            max_msrs_phiquad= 100
      INTEGER ::          n_muquad, n_phiquad

!  arrays

      REAL(FPK) :: X_MUQUAD (max_msrs_muquad)
      REAL(FPK) :: W_MUQUAD (max_msrs_muquad)
      REAL(FPK) :: SX_MUQUAD (max_msrs_muquad)
      REAL(FPK) :: WXX_MUQUAD(max_msrs_muquad)

      REAL(FPK) :: X_PHIQUAD (max_msrs_phiquad)
      REAL(FPK) :: W_PHIQUAD (max_msrs_phiquad)

      REAL(FPK) :: R0_QUAD_IN  (max_msrs_muquad,max_msrs_phiquad)
      REAL(FPK) :: R0_OUT_QUAD (max_msrs_muquad,max_msrs_phiquad)
      REAL(FPK) :: &
             D_R0_QUAD_IN  (max_msrs_muquad,max_msrs_phiquad), &
             D_R0_OUT_QUAD (max_msrs_muquad,max_msrs_phiquad)

!  Safety first zeroing

      REFLEC_0 = ZERO
      REFLEC_1 = ZERO

      COXMUNK_KERNEL     = ZERO
      COXMUNK_DERIVATIVE = ZERO

      D_REFLEC_0 = ZERO
      D_REFLEC_1 = ZERO

!  Air to water, Polar quadrature

      n_muquad = 40
      CALL GAULEG ( ZERO, ONE, X_muquad, W_muquad, n_muquad )
      DO I = 1, N_MUQUAD
        XM = X_MUQUAD(I)
        SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
        WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
      ENDDO

!  Azimuth quadrature

      n_phiquad = 100
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

      CALL COXMUNK_FUNCTION_PLUS &
             ( MAXPARS, PARS, PAR_INDEX, &
               XJ, SXJ, XI, SXI, PHI, CPHI, &
               REFLEC_0, D_REFLEC_0 )

!  Quadrature output for first order R/T calculations

      DO K = 1, n_muquad
        XM  = X_MUQUAD(K)
        SXM = SX_MUQUAD(K)
        DO N = 1, N_PHIQUAD
          PHI_SUB1  = X_PHIQUAD(N)
          CPHI_SUB1 = DCOS(PHI_SUB1)
          PHI_SUB2  = PHI*DEG_TO_RAD - X_PHIQUAD(N)
          CPHI_SUB2 = DCOS(PHI_SUB2)
          CALL COXMUNK_FUNCTION_PLUS &
             ( MAXPARS, PARS, PAR_INDEX, &
               XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, &
               R0_OUT_QUAD(K,N), D_R0_OUT_QUAD(K,N) )
! mick fix 4/22/11 - changed call statement
!          CALL COXMUNK_FUNCTION_PLUS &
!             ( MAXPARS, PARS, PAR_INDEX, &
!               XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, &
!               R0_QUAD_IN(K,N), D_R0_OUT_QUAD(K,N) )
          CALL COXMUNK_FUNCTION_PLUS &
             ( MAXPARS, PARS, PAR_INDEX, &
               XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, &
               R0_QUAD_IN(K,N), D_R0_QUAD_IN(K,N) )
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

!  Derivative

      SUMR = ZERO
      DO K = 1, n_muquad
        PR = ZERO
        DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          PR = PR + W_P *   R0_QUAD_IN(K,N) * D_R0_OUT_QUAD(K,N) &
                  + W_P * D_R0_QUAD_IN(K,N) *   R0_OUT_QUAD(K,N)
        ENDDO
        SUMR = SUMR + PR * WXX_MUQUAD(K)
      ENDDO
      D_REFLEC_1 = SUMR
      COXMUNK_DERIVATIVE = D_REFLEC_0 + D_REFLEC_1

!      write(34,'(1p6e14.5)')reflec_0, reflec_1, pars(1),
!     &   dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_FUNCTION_PLUS_DB

      END MODULE lrrs_brdf_plus

