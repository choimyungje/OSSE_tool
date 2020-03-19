
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
! #   Master SUBROUTINES (Public, these are stand-alone):       #
! #                                                             #
! #            LINLRRS_RAMAN_XSECS_BIN                          #
! #            LINLRRS_RAMAN_XSECS_MONO                         #
! #                                                             #
! #   SUBROUTINES (Private) :                                   #
! #                                                             #
! #                LINRRS_LOSSTERM_XSEC                         #
! #                LINRRS_GAINTERM_XSEC_BIN                     #
! #                LINRRS_GAINTERM_XSEC_MONO                    #
! #                                                             #
! #  Monochromatic Routines added 8 February 2011               #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Bookkeeping improvements (use of "Only", clearer I/O specifications)

      MODULE lrrs_L_raman_spectroscopy_m

      USE LRRS_PARS_m, Only : FPK, zero, One, half, two, PIE, MAX_BINS, MAX_MESSAGES, LDU

      PRIVATE :: LINRRS_LOSSTERM_XSEC,      &
                 LINRRS_GAINTERM_XSEC_BIN,  &
                 LINRRS_GAINTERM_XSEC_MONO

      PUBLIC  :: LINLRRS_RAMAN_XSECS_BIN,&
                 LINLRRS_RAMAN_XSECS_MONO

      CONTAINS

      SUBROUTINE LINLRRS_RAMAN_XSECS_BIN &
        ( MAX_LAYERS, MAX_POINTS, MAX_BINS,                             & ! Inputs
          NLAYERS, NPOINTS_INNER, NPOINTS_OUTER, OFFSET_INNER,          & ! Inputs
          BINLOWER, BINUPPER, LAMBDAS_RANKED, FLUXES_RANKED,            & ! Inputs
          N_RRS_N2REF, N_RRS_O2REF, RRS_N2REF, RRS_O2REF,               & ! Inputs
          LINRRS_N2REF, LINRRS_O2REF, N2POS, O2POS, GAMMA_N2, GAMMA_O2, & ! Inputs
          N_RRSBINS, BINMAP, BINXSEC, RRSXSEC_OUT, RINGSPEC_1,          & ! Output
          BINXSEC_dT, RRSXSEC_OUT_dT,                                   & ! Output
          FAIL, MESSAGE, ACTION )                                         ! Output

!  Purpose:
!  ========

!  This is stand alone routine for generation of Raman Loss and Gain
!  Cross-sections in each layer, plus the temperature derivative of
!  these cross-sections.

!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@
!  Exception handling for MAX_BINS Dimension

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Maximum number of spectral points and bins

      INTEGER  , INTENT(IN) :: MAX_POINTS, MAX_BINS

!  Maximum number of layers

      INTEGER  , INTENT(IN) :: MAX_LAYERS

!  Layering

      INTEGER  , INTENT(IN) :: NLAYERS

!  Number of points between the limiting wavelengths

      INTEGER  , INTENT(IN) :: NPOINTS_INNER

!  outer wavelength and offset

      INTEGER  , INTENT(IN) :: NPOINTS_OUTER, OFFSET_INNER

!  Wavelengths and Fluxes
!  These are defined on the outer grid for the binning realization

      REAL(FPK), INTENT(IN) :: LAMBDAS_RANKED  ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: FLUXES_RANKED   ( MAX_POINTS )

!  Wavelength bins

      REAL(FPK), INTENT(IN) :: BINLOWER ( MAX_POINTS )
      REAL(FPK), INTENT(IN) :: BINUPPER ( MAX_POINTS )

!  Dimensions

      INTEGER  , INTENT(IN) :: N_RRS_O2REF
      INTEGER  , INTENT(IN) :: N_RRS_N2REF

!  Basic coefficients from the previous module (see above)

      REAL(FPK), INTENT(IN) :: RRS_N2REF ( MAX_LAYERS, N_RRS_N2REF )
      REAL(FPK), INTENT(IN) :: RRS_O2REF ( MAX_LAYERS, N_RRS_O2REF )


!  Linearized Basic coefficients from the Spectroscopy module

      REAL(FPK), INTENT(IN) :: LINRRS_N2REF ( MAX_LAYERS, N_RRS_N2REF )
      REAL(FPK), INTENT(IN) :: LINRRS_O2REF ( MAX_LAYERS, N_RRS_O2REF )

!  Positions

      REAL(FPK), INTENT(IN) :: O2POS (N_RRS_O2REF)
      REAL(FPK), INTENT(IN) :: N2POS (N_RRS_N2REF)

!  Gammas

      REAL(FPK), INTENT(IN) :: GAMMA_O2(3)
      REAL(FPK), INTENT(IN) :: GAMMA_N2(3)

!  Output arguments
!  ================

!  Bin Mapping quantities (internal derivation)

      INTEGER  , INTENT(OUT) :: N_RRSBINS  ( MAX_POINTS )
      INTEGER  , INTENT(OUT) :: BINMAP ( MAX_BINS, MAX_POINTS )

!  Binned Gain term cross-sections

      REAL(FPK), INTENT(OUT) :: BINXSEC ( MAX_LAYERS, MAX_POINTS, MAX_BINS )

!  Loss term cross sections

      REAL(FPK), INTENT(OUT) :: RRSXSEC_OUT ( MAX_LAYERS, MAX_POINTS )

!  Chance/Spurr Ring Spectrum. DIAGNOSTIC ONLY.

      REAL(FPK), INTENT(OUT) :: RINGSPEC_1 ( MAX_POINTS )

!  T-derivative of Binned Gain term cross-sections

      REAL(FPK), INTENT(OUT) :: BINXSEC_dT (MAX_LAYERS,MAX_POINTS,MAX_BINS)

!  T-derivative of Loss term cross sections

      REAL(FPK), INTENT(OUT) :: RRSXSEC_OUT_dT ( MAX_LAYERS, MAX_POINTS )

!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@
!  Exception handling for MAX_BINS Dimension

      LOGICAL          , INTENT(OUT)   :: FAIL
      CHARACTER (LEN=*), INTENT(OUT)   :: MESSAGE, ACTION

!  Local variables
!  ===============

!  1. For the RRS spectroscic arrays and associated quantitites
!  ------------------------------------------------------------

!  Shift numbers

      INTEGER, PARAMETER :: MAX_SHIFTS = 233

!  number of shifts = N_Stokes + N_Antistokes

      INTEGER   :: NSHIFTS, N_STOKES, N_ANTISTOKES

!  Gain term cross-sections and wavelengths (shifted)

      REAL(FPK) :: LAMBDAS_SHIFTED    ( MAX_SHIFTS )
      REAL(FPK) :: RRSXSEC_SHIFTED    ( MAX_LAYERS, MAX_SHIFTS )

!  Stokes and Antistokes cross sections (debug only)

      REAL(FPK) :: RRSXSEC_STOKES     ( MAX_LAYERS, MAX_SHIFTS )
      REAL(FPK) :: RRSXSEC_ANTISTOKES ( MAX_LAYERS, MAX_SHIFTS )

!  T-derivatives of Raman cross-sections (all)

      REAL(FPK) :: LINRRSXSEC_SHIFTED ( MAX_LAYERS, MAX_SHIFTS )

!  2. Help variables
!  -----------------

!  C3 and BMAX added, 14 April 2011

      LOGICAL   :: LOOP
      INTEGER   :: I, N, B, WI, WO, BMAX
      INTEGER   :: K, IK, J
      REAL(FPK) :: CONTRIB, RMS
      CHARACTER (LEN=3) :: C3

!  Initial section
!  ===============

!  Zero the binned maps and cross-sections

      NSHIFTS = 233
      DO K = 1, NPOINTS_INNER
       N_RRSBINS(K) = 0
       DO B = 1, MAX_BINS
        BINMAP(B,K) = 0
        DO N = 1, NLAYERS
         BINXSEC(N,K,B)    = zero
         BINXSEC_dT(N,K,B) = zero
        ENDDO
       ENDDO
      ENDDO

!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@
!  Initialize exception handling

      FAIL = .false.
      message = ' '
      action  = ' '
      BMAX = -9999

!  Gain cross-sections (binned)
!  ----------------------------

!  Start outer window loop

      DO I = 1, NPOINTS_OUTER

!  Get the RRS cross sections from point I in this window

        CALL LINRRS_GAINTERM_XSEC_BIN &
           ( MAX_LAYERS, MAX_SHIFTS,                  & ! Inputs
             NLAYERS, N_RRS_N2REF, N_RRS_O2REF,       & ! Inputs
             N2POS, O2POS, GAMMA_N2, GAMMA_O2,        & ! Inputs
             LAMBDAS_RANKED(I), RRS_N2REF, RRS_O2REF, & ! Inputs
             LINRRS_N2REF, LINRRS_O2REF,              & ! Inputs
             N_STOKES,     RRSXSEC_STOKES,            & ! Outputs
             N_ANTISTOKES, RRSXSEC_ANTISTOKES,        & ! Outputs
             LAMBDAS_SHIFTED, RRSXSEC_SHIFTED,        & ! Outputs
             LINRRSXSEC_SHIFTED )                       ! Outputs

!  Start loop over points in inner window
!    If there any contributions, add them to this bin, assign map.

        DO K = 1, NPOINTS_INNER

          IK = K + OFFSET_INNER
          J = 0
          LOOP = .TRUE.

!    Check if any shifts are within the bin of the inner point IK

          DO WHILE ( LOOP .AND. J.LT.NSHIFTS )
            J = J + 1
            IF ( LAMBDAS_SHIFTED(J).GT.BINLOWER(IK).AND. &
                 LAMBDAS_SHIFTED(J).LT.BINUPPER(IK) ) LOOP = .FALSE.
          ENDDO

!    If there any contributions add them to this bin, assign map.

          IF ( .NOT. LOOP ) THEN

! @@@@@@@ RT Solutions, 14 April 2011
!  New code, checks maximum value of B
!      only assigns map and cross-sections if this value is not exceeded

            B = N_RRSBINS(K) + 1
            N_RRSBINS(K) = B
            BMAX = MAX( B, BMAX )
            IF ( B .LE. MAX_BINS ) THEN
               BINMAP(B,K) = I
               DO J = 1, NSHIFTS
                 IF ( LAMBDAS_SHIFTED(J).GT.BINLOWER(IK).AND. &
                      LAMBDAS_SHIFTED(J).LT.BINUPPER(IK) )  THEN
                   DO N = 1, NLAYERS
                     CONTRIB = RRSXSEC_SHIFTED(N,J)
                     BINXSEC(N,K,B) = BINXSEC(N,K,B) + CONTRIB
                     CONTRIB = LINRRSXSEC_SHIFTED(N,J)
                     BINXSEC_dT(N,K,B) = BINXSEC_dT(N,K,B) + CONTRIB
                   ENDDO
                 ENDIF
               ENDDO
            ENDIF

!  Old code
!            B = N_RRSBINS(K) + 1
!            BINMAP(B,K) = I
!            N_RRSBINS(K) = B
!            DO J = 1, NSHIFTS
!              IF ( LAMBDAS_SHIFTED(J).GT.BINLOWER(IK).AND. &
!                   LAMBDAS_SHIFTED(J).LT.BINUPPER(IK) )  THEN
!                DO N = 1, NLAYERS
!                  CONTRIB = RRSXSEC_SHIFTED(N,J)
!                  BINXSEC(N,K,B) = BINXSEC(N,K,B) + CONTRIB
!                 CONTRIB = LINRRSXSEC_SHIFTED(N,J)
!                  BINXSEC_dT(N,K,B) = BINXSEC_dT(N,K,B) + CONTRIB
!               ENDDO
!              ENDIF
!            ENDDO

!  Finish contributions

          ENDIF

!  Finish inner loop

        ENDDO

!  Finish outer loop for Gain cross-sections

      ENDDO

! @@@@@@@ RT Solutions, 14 April 2011
!  Exception handling, when Binning dimension exceeded

      IF ( BMAX.GT.MAX_BINS ) THEN
         FAIL    = .true.
         write(C3,'(I3)')BMAX
         MESSAGE = ' Dimension MAX_BINS exceeded in Subroutine LINLRRS_RAMAN_XSECS_BIN'
         ACTION  = ' Increase MAX_BINS in lrrs_pars.f90 to at least '//c3
         RETURN
      ENDIF

!  debug and ring spec
!  -------------------

      do k = 1, npoints_inner
       ringspec_1(k) = zero
       ik = k + offset_inner
       do b = 1, n_rrsbins(k)
         ringspec_1(k) = ringspec_1(k) + &
           fluxes_ranked(binmap(b,k))*binxsec(1,k,b)
       enddo
        ringspec_1(k)  = ringspec_1(k) /  fluxes_ranked(ik)
      enddo
      rms = zero
      do k = 1, npoints_inner
        rms = rms + 1.0d+20 * ringspec_1(k) * ringspec_1(k)
      enddo

      rms = 1.0d-10*dsqrt(rms/dble(npoints_inner))
      do k = 1, npoints_inner
       ringspec_1(k) = ringspec_1(k) / rms
      enddo

!  Loss term RRS
!  -------------

!  For each wavelength in the inner loop, compute the total cross-sectio
!  for RRS light scattered OUT OF THIS WAVELENGTH. This is the Loss term
!    Requires temperatures as input, and the basic reference values

      DO WI = 1, NPOINTS_INNER
        WO = WI + OFFSET_INNER
        CALL LINRRS_LOSSTERM_XSEC &
            ( MAX_LAYERS, NLAYERS, N_RRS_N2REF, N_RRS_O2REF,        & ! Inputs
              N2POS, O2POS, GAMMA_N2, GAMMA_O2, LAMBDAS_RANKED(WO), & ! Inputs
              RRS_N2REF, RRS_O2REF, LINRRS_N2REF, LINRRS_O2REF,     & ! Inputs
              RRSXSEC_OUT(1,WI), RRSXSEC_OUT_dT(1,WI), NSHIFTS )      ! Outputs
      enddo

!  Finish routine

      return
      END SUBROUTINE LINLRRS_RAMAN_XSECS_BIN

!

      SUBROUTINE LINLRRS_RAMAN_XSECS_MONO &
        ( MAX_LAYERS, MAX_POINTS, NLAYERS, LAMBDA_EXCIT,  & ! Inputs
          N_RRS_N2REF, N_RRS_O2REF, RRS_N2REF, RRS_O2REF, & ! Inputs
          LINRRS_N2REF, LINRRS_O2REF,                     & ! Inputs
          N2POS, O2POS, GAMMA_N2, GAMMA_O2,               & ! Inputs
          RRSXSEC_OUT_MONO,    RRSXSEC_RANKED,            & ! Outputs
          RRSXSEC_OUT_MONO_dT, RRSXSEC_RANKED_dT )          ! Outputs

!  This is stand alone.

!  Routine is called Inside LRRS_L_Setups_Master module; the only inputs
!  the monochromatic wavelength LAMBDA_EXCIT and spectroscopy
!    Ranked Raman cross-sections are the output + Loss term.
!    Ranked Wavelengths are generated automatically but not output here.
!     (Already generated in the setup routine RRS_LAMDASRANKED_MONO).

      USE LRRS_AUX2_m, Only : RQSORT_IDX

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Maximum number of spectral points and layers

      INTEGER  , INTENT(IN) :: MAX_POINTS
      INTEGER  , INTENT(IN) :: MAX_LAYERS

!  Actual number

      INTEGER  , INTENT(IN) :: NLAYERS

!  Excitation wavelength

      REAL(FPK), INTENT(IN) :: LAMBDA_EXCIT

!  Dimensions

      INTEGER  , INTENT(IN) :: N_RRS_O2REF
      INTEGER  , INTENT(IN) :: N_RRS_N2REF

!  Basic coefficients from the Spectroscopy module

      REAL(FPK), INTENT(IN) :: RRS_N2REF ( MAX_LAYERS, N_RRS_N2REF )
      REAL(FPK), INTENT(IN) :: RRS_O2REF ( MAX_LAYERS, N_RRS_O2REF )

!  Linearized Basic coefficients from the Spectroscopy module

      REAL(FPK), INTENT(IN) :: LINRRS_N2REF ( MAX_LAYERS, N_RRS_N2REF )
      REAL(FPK), INTENT(IN) :: LINRRS_O2REF ( MAX_LAYERS, N_RRS_O2REF )

!  Positions

      REAL(FPK), INTENT(IN) :: O2POS (N_RRS_O2REF)
      REAL(FPK), INTENT(IN) :: N2POS (N_RRS_N2REF)

!  Gammas

      REAL(FPK), INTENT(IN) :: GAMMA_O2(3)
      REAL(FPK), INTENT(IN) :: GAMMA_N2(3)

!  Output arguments
!  ================

!  Gain term cross-sections

      REAL(FPK), INTENT(OUT) :: RRSXSEC_RANKED ( MAX_LAYERS, MAX_POINTS )

!  Loss term cross sections

      REAL(FPK), INTENT(OUT) :: RRSXSEC_OUT_MONO    ( MAX_LAYERS )

!  Temperature linearizations of cross-sections

      REAL(FPK), INTENT(OUT) :: RRSXSEC_RANKED_dT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: RRSXSEC_OUT_MONO_dT    ( MAX_LAYERS )

!  Local variables
!  ===============

!  Shift numbers

      INTEGER, PARAMETER :: MAX_SHIFTS   = 233
      INTEGER, PARAMETER :: MAX_SHIFTS_1 = 234

!  Number of points
!  number of shifts = N_Stokes + N_Antistokes

      INTEGER   :: NPOINTS_MONO
      INTEGER   :: NSHIFTS

!  Gain term cross-sections and wavelengths (shifted) of excitatoion)
!   (Includes separately the STokes and AntiStokes cross-sections)

      REAL(FPK) :: LAMBDAS_SHIFTED    ( MAX_SHIFTS )
      REAL(FPK) :: RRSXSEC_SHIFTED    ( MAX_LAYERS, MAX_SHIFTS )

!  T-derivatives of Raman cross-sections (all)

      REAL(FPK) :: RRSXSEC_SHIFTED_dT  ( MAX_LAYERS,MAX_SHIFTS )

!  Ranking

      INTEGER   :: INDEX_ALL           ( MAX_SHIFTS_1 )
      REAL(FPK) :: LAMBDAS_ALL         ( MAX_SHIFTS_1 )
      REAL(FPK) :: RRSXSEC_ALL1        ( MAX_LAYERS, MAX_SHIFTS_1 )
      REAL(FPK) :: RRSXSEC_ALL2        ( MAX_LAYERS, MAX_SHIFTS_1 )

!  Help variables

      INTEGER   :: I, W, N

!  Initial section
!  ===============

!  Number of points is always the same

      NPOINTS_MONO = MAX_SHIFTS_1

!  Zero the cross-sections and T-linearizations

      RRSXSEC_OUT_MONO = zero
      RRSXSEC_RANKED   = zero

      RRSXSEC_OUT_MONO_dT = zero
      RRSXSEC_RANKED_dT   = zero

!  Loss term RRS
!  -------------

!  For each wavelength in the inner loop, compute the total cross-sectio
!  for RRS light scattered OUT OF THIS WAVELENGTH. This is the Loss term
!    Requires temperatures as input and the basic reference values

      CALL LINRRS_LOSSTERM_XSEC &
            ( MAX_LAYERS, NLAYERS, N_RRS_N2REF, N_RRS_O2REF,    & ! Inputs
              N2POS, O2POS, GAMMA_N2, GAMMA_O2, LAMBDA_EXCIT,   & ! Inputs
              RRS_N2REF, RRS_O2REF, LINRRS_N2REF, LINRRS_O2REF, & ! Inputs
              RRSXSEC_OUT_MONO, RRSXSEC_OUT_MONO_dT, NSHIFTS )    ! Outputs

!  Gain term RRS (Mono
!  -------------------

!  find shifted wavelengths and cross sections for the RRS Gain terms =
!    ineslatic scattering of light INTO the wavelength LAMBDAS(wo)
!       Not sorted or ranked

! mick fix 4/22/11 - replaced call statement
!      CALL LINRRS_GAINTERM_XSEC_MONO &
!           ( MAX_LAYERS, MAX_SHIFTS, NLAYERS, LAMBDA_EXCIT, &
!             RRS_N2REF, RRS_O2REF, LINRRS_N2REF, LINRRS_O2REF, &
!             LAMBDAS_SHIFTED, RRSXSEC_SHIFTED, RRSXSEC_SHIFTED_dT )

      CALL LINRRS_GAINTERM_XSEC_MONO &
           ( MAX_LAYERS, MAX_SHIFTS, NLAYERS, N_RRS_N2REF, N_RRS_O2REF, & ! Inputs
             N2POS, O2POS, GAMMA_N2, GAMMA_O2, LAMBDA_EXCIT,            & ! Inputs
             RRS_N2REF, RRS_O2REF, LINRRS_N2REF, LINRRS_O2REF,          & ! Inputs
             LAMBDAS_SHIFTED, RRSXSEC_SHIFTED, RRSXSEC_SHIFTED_dT )       ! Outputs

!  Ranking procedure
!  =================

!  Copy to holding arrays and add Excitation wavelength
!  Number of points is always the same

      NPOINTS_MONO = MAX_SHIFTS_1
      DO W = 1, NSHIFTS
        LAMBDAS_ALL(W) = LAMBDAS_SHIFTED(W)
        DO N = 1, NLAYERS
          RRSXSEC_ALL1(N,W) = RRSXSEC_SHIFTED(N,W)
          RRSXSEC_ALL2(N,W) = RRSXSEC_SHIFTED_dT(N,W)
        ENDDO
      ENDDO
      LAMBDAS_ALL(NPOINTS_MONO) = LAMBDA_EXCIT
      DO N = 1, NLAYERS
        RRSXSEC_ALL1(N,NPOINTS_MONO) = zero
        RRSXSEC_ALL2(N,NPOINTS_MONO) = zero
      ENDDO

!  Use RQSORT_IDX to rank the cross sections

      DO I=1,NPOINTS_MONO
        INDEX_ALL(I) = I
      ENDDO
      CALL RQSORT_IDX(NPOINTS_MONO, LAMBDAS_ALL, NPOINTS_MONO, INDEX_ALL)
      DO W = 1, NPOINTS_MONO
        DO N = 1, NLAYERS
          RRSXSEC_RANKED(N,W)    = RRSXSEC_ALL1(N,INDEX_ALL(W))
          RRSXSEC_RANKED_dT(N,W) = RRSXSEC_ALL2(N,INDEX_ALL(W))
        ENDDO
      ENDDO

!  Finish routine

      return
      END SUBROUTINE LINLRRS_RAMAN_XSECS_MONO

!

      SUBROUTINE LINRRS_LOSSTERM_XSEC &
            ( MAXLAYERS, NLAYERS, N_RRS_N2REF, N_RRS_O2REF,       & ! Inputs
              N2POS, O2POS, GAMMA_N2, GAMMA_O2, WAVELENGTH_INPUT, & ! Inputs
              RRS_N2REF, RRS_O2REF, LINRRS_N2REF, LINRRS_O2REF,   & ! Inputs
              RRSXSEC_OUT, LINRRSXSEC_OUT, N_RRSOUT )               ! Outputs

!  Loss term cross section
!  =======================

!  Given a wavelength WAVELENGTH_INPUT, and a bunch of temperatures
!  the compute the complete  cross-section of energy scattered out
!  of this wavelength. Requires as input the basic transition factors
!  from the previous module. Calculates internally the wavelength shifts
!  out from WAVELENGTH_INPUT, and includes lambda^-4 factor in result.

!  Get the derivatives of these cross-sections w.r.t Temperatures
!   ...Coded by R. Spurr, 29 July 2008

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Number of layers

      INTEGER  , INTENT(IN) :: MAXLAYERS, NLAYERS

!  Dimensions

      INTEGER  , INTENT(IN) :: N_RRS_O2REF
      INTEGER  , INTENT(IN) :: N_RRS_N2REF

!  Positions

      REAL(FPK), INTENT(IN) :: O2POS (N_RRS_O2REF)
      REAL(FPK), INTENT(IN) :: N2POS (N_RRS_N2REF)

!  Gammas

      REAL(FPK), INTENT(IN) :: GAMMA_O2(3)
      REAL(FPK), INTENT(IN) :: GAMMA_N2(3)

!  Input wavelength

      REAL(FPK), INTENT(IN) :: WAVELENGTH_INPUT

!  Basic coefficients from the previous module (see above)

      REAL(FPK), INTENT(IN) :: RRS_N2REF ( MAXLAYERS, N_RRS_N2REF )
      REAL(FPK), INTENT(IN) :: RRS_O2REF ( MAXLAYERS, N_RRS_O2REF )

!  T-derivatives of Basic coefficients for O2 and N2

      REAL(FPK), INTENT(IN) :: LINRRS_N2REF ( MAXLAYERS, N_RRS_N2REF )
      REAL(FPK), INTENT(IN) :: LINRRS_O2REF ( MAXLAYERS, N_RRS_O2REF )

!  Output arguments
!  ----------------

!  Number of Raman shifts

      INTEGER  , INTENT(OUT) :: N_RRSOUT

!  Total Loss-term cross sections, one per layer

      REAL(FPK), INTENT(OUT) :: RRSXSEC_OUT ( MAXLAYERS )

!  T-derivative of Total Loss-term cross sections, one per layer

      REAL(FPK), INTENT(OUT) :: LINRRSXSEC_OUT ( MAXLAYERS )

!  Local variables
!  ---------------

      INTEGER   :: N, J
      REAL(FPK) :: SIGMA_0, SIGMA_4, SUMR, N2XSEC, O2XSEC
      REAL(FPK) :: SIGMA, SIGSQ, DRV, D_N2XSEC, D_O2XSEC
      REAL(FPK) :: GAMMAN2,    GAMMAO2
      REAL(FPK) :: GAMMAN2_SQ, GAMMAO2_SQ
      REAL(FPK) :: HELPN2 ( 48 )
      REAL(FPK) :: HELPO2 ( 185 )

      REAL(FPK), PARAMETER :: DCONV24 = 1.0D-24

!  wavelength is in nm

      SIGMA = 1.0D+3 / WAVELENGTH_INPUT
      SIGSQ = SIGMA * SIGMA
      GAMMAN2 = GAMMA_N2(1) + GAMMA_N2(2) / (GAMMA_N2(3) - SIGSQ)
      GAMMAO2 = GAMMA_O2(1) + GAMMA_O2(2) / (GAMMA_O2(3) - SIGSQ)
      GAMMAN2_SQ = GAMMAN2 * GAMMAN2
      GAMMAO2_SQ = GAMMAO2 * GAMMAO2
      SIGMA_0 = 1.0d+07 / WAVELENGTH_INPUT

!  wavelength shifts for N2

      DO J = 1, N_RRS_N2REF
        SIGMA   = SIGMA_0 - N2POS(J)
        SIGMA_4 = SIGMA * SIGMA * SIGMA * SIGMA
        HELPN2(J) = GAMMAN2_SQ * SIGMA_4
      ENDDO

!  Wavelength shifts for O2

      DO J = 1, N_RRS_O2REF
        SIGMA = SIGMA_0 - O2POS(J)
        SIGMA_4 = SIGMA * SIGMA * SIGMA * SIGMA
        HELPO2(J) = GAMMAO2_SQ * SIGMA_4
      ENDDO

!  Generate Loss-term cross sections

      DO N = 1, NLAYERS
        SUMR = zero
        DRV  = zero
        DO J = 1, N_RRS_N2REF
          N2XSEC    = RRS_N2REF   (N,J) * HELPN2(J)
          D_N2XSEC  = LINRRS_N2REF(N,J) * HELPN2(J)
          SUMR = SUMR + N2XSEC
          DRV  = DRV  + D_N2XSEC
        ENDDO
        DO J = 1, N_RRS_O2REF
          O2XSEC    = RRS_O2REF   (N,J) * HELPO2(J)
          D_O2XSEC  = LINRRS_O2REF(N,J) * HELPO2(J)
          SUMR = SUMR + O2XSEC
          DRV  = DRV  + D_O2XSEC
        ENDDO
        RRSXSEC_OUT(N)    = SUMR
        LINRRSXSEC_OUT(N) = DRV
      ENDDO

!   we need to multiply by d-48 to get the correct
!   units for the absolute cross sections in cm2: We calculate gamman2
!   and gammao2 in d-30 m3 and then calculate n2xsec and o2xsec values
!   as (gamma**2) * (sigprime**4), where sigprime is in cm-1.

      DO N = 1, NLAYERS
        RRSXSEC_OUT(N)    = RRSXSEC_OUT(N)    * DCONV24 * DCONV24
        LINRRSXSEC_OUT(N) = LINRRSXSEC_OUT(N) * DCONV24 * DCONV24
      ENDDO

!  number of shifts

      N_RRSOUT = N_RRS_N2REF + N_RRS_O2REF

!  Finish

      RETURN
      END SUBROUTINE LINRRS_LOSSTERM_XSEC

!

      SUBROUTINE LINRRS_GAINTERM_XSEC_BIN &
           ( MAXLAYERS, MAXSHIFTS,                   & ! Inputs
             NLAYERS, N_RRS_N2REF, N_RRS_O2REF,      & ! Inputs
             N2POS, O2POS, GAMMA_N2, GAMMA_O2,       & ! Inputs
             WAVELENGTH_INPUT, RRS_N2REF, RRS_O2REF, & ! Inputs
             LINRRS_N2REF, LINRRS_O2REF,             & ! Inputs
             N_STOKES,     RRSXSEC_STOKES,           & ! Inputs
             N_ANTISTOKES, RRSXSEC_ANTISTOKES,       & ! Inputs
             WAVELENGTHS_SHIFTED, RRSXSEC_SHIFTED,   & ! Outputs
             LINRRSXSEC_SHIFTED )                      ! Outputs

!  Gain term cross-sections
!  ========================

!  Given a wavelength WAVELENGTH_INPUT, and a bunch of temperatures
!  the compute all cross-section of energy scattered INTO this
!  wavelength from the RRS shifted excitation wavelengths. Requires as
!  input the basic transition factors from previous module.
!  Calculates and outputs the shifted wavelengths WAVELENGTHS_SHIFTED
!  into WAVELENGTH_INPUT, and includes lambda^-4 factor in results of
!  cross sections RRSXSEC_SHIFTED.

!  Module also outputs the Stokes and Antistokes cross-sections, but
!  these are only for information.

      IMPLICIT NONE

!  Input arguments
!  ---------------

      INTEGER  , INTENT(IN) :: MAXLAYERS, MAXSHIFTS

!  actual number of layers

      INTEGER  , INTENT(IN) :: NLAYERS

!  Dimensions

      INTEGER  , INTENT(IN) :: N_RRS_O2REF
      INTEGER  , INTENT(IN) :: N_RRS_N2REF

!  Positions

      REAL(FPK), INTENT(IN) :: O2POS (N_RRS_O2REF)
      REAL(FPK), INTENT(IN) :: N2POS (N_RRS_N2REF)

!  Gammas

      REAL(FPK), INTENT(IN) :: GAMMA_O2(3)
      REAL(FPK), INTENT(IN) :: GAMMA_N2(3)

!  Excitation wavelength

      REAL(FPK), INTENT(IN) :: WAVELENGTH_INPUT

!  Basic coefficients from the previous module (see above)

      REAL(FPK), INTENT(IN) :: RRS_N2REF ( MAXLAYERS, N_RRS_N2REF )
      REAL(FPK), INTENT(IN) :: RRS_O2REF ( MAXLAYERS, N_RRS_O2REF )

!  T-derivatives of Basic coefficients for O2 and N2

      REAL(FPK), INTENT(IN) :: LINRRS_N2REF ( MAXLAYERS, N_RRS_N2REF )
      REAL(FPK), INTENT(IN) :: LINRRS_O2REF ( MAXLAYERS, N_RRS_O2REF )

!  Output arguments
!  ----------------

!  Raman shifted wavelengths

      REAL(FPK), INTENT(OUT) :: WAVELENGTHS_SHIFTED ( MAXSHIFTS )

!  Raman cross-sections (all)

      REAL(FPK), INTENT(OUT) :: RRSXSEC_SHIFTED     ( MAXLAYERS, MAXSHIFTS )

!  T-derivatives of Raman cross-sections (all)

      REAL(FPK), INTENT(OUT) :: LINRRSXSEC_SHIFTED  ( MAXLAYERS, MAXSHIFTS )

!  Stokes and Antistokes cross sections (debug only)

      REAL(FPK), INTENT(OUT) :: RRSXSEC_STOKES      ( MAXLAYERS, MAXSHIFTS )
      REAL(FPK), INTENT(OUT) :: RRSXSEC_ANTISTOKES  ( MAXLAYERS, MAXSHIFTS )

!  Local variables
!  ---------------

      INTEGER   :: N, J, S, NSHIFTS, N_STOKES, N_ANTISTOKES, S1, S2
      REAL(FPK) :: ANTISTOKES_SHIFTED ( 233 )
      REAL(FPK) :: STOKES_SHIFTED     ( 233 )
      REAL(FPK) :: SIGMA, SIGSQ, SIGMA_0, SIGMA_2, SIGMA_4
      REAL(FPK) :: GAMMAN2, GAMMAO2, GAMMAN2_SQ, GAMMAO2_SQ
      REAL(FPK) :: HELPN2, HELPO2

      REAL(FPK), PARAMETER  :: D24 = 1.0D-24

!  baseline values of GAMMA

      SIGMA_0 = 1.0D+07 / WAVELENGTH_INPUT
      SIGSQ   = 1.0D+06 / WAVELENGTH_INPUT / WAVELENGTH_INPUT
      GAMMAN2 = GAMMA_N2(1) + GAMMA_N2(2) / (GAMMA_N2(3) - SIGSQ)
      GAMMAO2 = GAMMA_O2(1) + GAMMA_O2(2) / (GAMMA_O2(3) - SIGSQ)
      GAMMAN2_SQ = GAMMAN2 * GAMMAN2
      GAMMAO2_SQ = GAMMAO2 * GAMMAO2

!  initialize counts

      S  = 0
      S1 = 0
      S2 = 0

!  Shifted wavelengths and cross sections (N2)

      DO J = 1, N_RRS_N2REF
        SIGMA   = SIGMA_0 - N2POS(J)
        SIGMA_2 = SIGMA * SIGMA
        SIGMA_4 = SIGMA_2 * SIGMA_2
        HELPN2  = GAMMAN2_SQ * SIGMA_4
        S = S + 1
        WAVELENGTHS_SHIFTED(S) = 1.0D+7 / SIGMA
        IF ( WAVELENGTHS_SHIFTED(S) .LT. WAVELENGTH_INPUT ) THEN
          S1 = S1 + 1
          ANTISTOKES_SHIFTED(S1)   = WAVELENGTHS_SHIFTED(S)
          DO N = 1, NLAYERS
            RRSXSEC_ANTISTOKES(N,S1) = RRS_N2REF(N,J) * HELPN2
          ENDDO
        ELSE IF ( WAVELENGTHS_SHIFTED(S) .GT. WAVELENGTH_INPUT ) THEN
          S2 = S2 + 1
          STOKES_SHIFTED(S2)   = WAVELENGTHS_SHIFTED(S)
          DO N = 1, NLAYERS
            RRSXSEC_STOKES(N,S2) = RRS_N2REF(N,J) * HELPN2
          ENDDO
        ENDIF
        DO N = 1, NLAYERS
          RRSXSEC_SHIFTED   (N,S) = RRS_N2REF   (N,J) * HELPN2
          LINRRSXSEC_SHIFTED(N,S) = LINRRS_N2REF(N,J) * HELPN2
        ENDDO
      ENDDO

!  Shifted wavelengths and cross sections (O2)

      DO J = 1, N_RRS_O2REF
        SIGMA = SIGMA_0 - O2POS(J)
        SIGMA_2 = SIGMA * SIGMA
        SIGMA_4 = SIGMA_2 * SIGMA_2
        HELPO2  = GAMMAO2_SQ * SIGMA_4
        S = S + 1
        WAVELENGTHS_SHIFTED(S) = 1.0D+7 / SIGMA
        IF ( WAVELENGTHS_SHIFTED(S) .LT. WAVELENGTH_INPUT ) THEN
          S1 = S1 + 1
          ANTISTOKES_SHIFTED(S1)   = WAVELENGTHS_SHIFTED(S)
          DO N = 1, NLAYERS
            RRSXSEC_ANTISTOKES(N,S1) = RRS_O2REF(N,J) * HELPO2
          ENDDO
        ELSE IF ( WAVELENGTHS_SHIFTED(S) .GT. WAVELENGTH_INPUT ) THEN
          S2 = S2 + 1
          STOKES_SHIFTED(S2)   = WAVELENGTHS_SHIFTED(S)
          DO N = 1, NLAYERS
            RRSXSEC_STOKES(N,S2) = RRS_O2REF(N,J) * HELPO2
          ENDDO
        ENDIF
        DO N = 1, NLAYERS
          RRSXSEC_SHIFTED   (N,S) =    RRS_O2REF(N,J) * HELPO2
          LINRRSXSEC_SHIFTED(N,S) = LINRRS_O2REF(N,J) * HELPO2
        ENDDO
      ENDDO

!  Total number of shifts

      NSHIFTS      = N_RRS_N2REF + N_RRS_O2REF
      N_ANTISTOKES = S1
      N_STOKES     = S2

!   we need to multiply by d-48 to get the correct
!   units for the absolute cross sections in cm2: We calculate gamman2
!   and gammao2 in d-30 m3 and then calculate n2xsec and o2xsec values
!   as (gamma**2) * (sigprime**4), where sigprime is in cm-1.

      DO S = 1, NSHIFTS
        DO N = 1, NLAYERS
          RRSXSEC_SHIFTED   (N,S) = RRSXSEC_SHIFTED   (N,S) * D24 * D24
          LINRRSXSEC_SHIFTED(N,S) = LINRRSXSEC_SHIFTED(N,S) * D24 * D24
        ENDDO
      ENDDO
      DO S = 1, N_STOKES
        DO N = 1, NLAYERS
          RRSXSEC_STOKES(N,S) = RRSXSEC_STOKES(N,S)*D24*D24
        ENDDO
      ENDDO
      DO S = 1, N_ANTISTOKES
        DO N = 1, NLAYERS
          RRSXSEC_ANTISTOKES(N,S) = RRSXSEC_ANTISTOKES(N,S)*D24*D24
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LINRRS_GAINTERM_XSEC_BIN

!

      SUBROUTINE LINRRS_GAINTERM_XSEC_MONO &
           ( MAXLAYERS, MAXSHIFTS,                   & ! Inputs
             NLAYERS, N_RRS_N2REF, N_RRS_O2REF,      & ! Inputs
             N2POS, O2POS, GAMMA_N2, GAMMA_O2,       & ! Inputs
             WAVELENGTH_INPUT, RRS_N2REF, RRS_O2REF, & ! Inputs
             LINRRS_N2REF, LINRRS_O2REF,             & ! Inputs
             WAVELENGTHS_SHIFTED, RRSXSEC_SHIFTED,   & ! Inputs
             LINRRSXSEC_SHIFTED )                      ! Outputs

!  Gain term cross-sections
!  ========================

!  Given a wavelength WAVELENGTH_INPUT, and a bunch of temperatures
!  the compute all cross-section of energy scattered INTO this
!  wavelength from the RRS shifted excitation wavelengths. Requires as
!  input the basic transition factors from previous module.
!  Calculates and outputs the shifted wavelengths WAVELENGTHS_SHIFTED
!  into WAVELENGTH_INPUT, and includes lambda^-4 factor in results of
!  cross sections RRSXSEC_SHIFTED.

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Dimensioning parameters

      INTEGER  , INTENT(IN) :: MAXLAYERS, MAXSHIFTS

!  actual number of layers

      INTEGER  , INTENT(IN) :: NLAYERS

!  Dimensions

      INTEGER  , INTENT(IN) :: N_RRS_O2REF
      INTEGER  , INTENT(IN) :: N_RRS_N2REF

!  Positions

      REAL(FPK), INTENT(IN) :: O2POS (N_RRS_O2REF)
      REAL(FPK), INTENT(IN) :: N2POS (N_RRS_N2REF)

!  Gammas

      REAL(FPK), INTENT(IN) :: GAMMA_O2(3)
      REAL(FPK), INTENT(IN) :: GAMMA_N2(3)

!  Excitation wavelength

      REAL(FPK), INTENT(IN) :: WAVELENGTH_INPUT

!  Basic coefficients from the previous module (see above)

      REAL(FPK), INTENT(IN) :: RRS_N2REF ( MAXLAYERS, N_RRS_N2REF )
      REAL(FPK), INTENT(IN) :: RRS_O2REF ( MAXLAYERS, N_RRS_O2REF )

!  T-derivatives of Basic coefficients for O2 and N2

      REAL(FPK), INTENT(IN) :: LINRRS_N2REF ( MAXLAYERS, N_RRS_N2REF )
      REAL(FPK), INTENT(IN) :: LINRRS_O2REF ( MAXLAYERS, N_RRS_O2REF )

!  Output arguments
!  ----------------

!  Raman shifted wavelengths

      REAL(FPK), INTENT(OUT) :: WAVELENGTHS_SHIFTED ( MAXSHIFTS )

!  Raman cross-sections (all)

      REAL(FPK), INTENT(OUT) :: RRSXSEC_SHIFTED     ( MAXLAYERS, MAXSHIFTS )

!  T-derivatives of Raman cross-sections (all)

      REAL(FPK), INTENT(OUT) :: LINRRSXSEC_SHIFTED  ( MAXLAYERS, MAXSHIFTS )

!  Local variables
!  ---------------

      INTEGER   :: N, J, S, NSHIFTS
      REAL(FPK) :: SIGMA, SIGSQ, SIGMA_0, SIGMA_2, SIGMA_4
      REAL(FPK) :: GAMMAN2, GAMMAO2, GAMMAN2_SQ, GAMMAO2_SQ
      REAL(FPK) :: HELPN2, HELPO2, LAMDA

      REAL(FPK), PARAMETER :: D24 = 1.0D-24

!  baseline values of GAMMA

      SIGMA_0 = 1.0D+07 / WAVELENGTH_INPUT
      SIGMA_2 = SIGMA_0 * SIGMA_0
      SIGMA_4 = SIGMA_2 * SIGMA_2

!  initialize counts

      S  = 0

!  Shifted wavelengths and cross sections (N2)

      DO J = 1, N_RRS_N2REF
        S = S + 1
        SIGMA = SIGMA_0 + N2POS(J)
        LAMDA = 1.0d+07 / SIGMA
        SIGSQ = 1.0d+06 / LAMDA / LAMDA
        GAMMAN2 = GAMMA_N2(1) + GAMMA_N2(2) / (GAMMA_N2(3) - SIGSQ)
        GAMMAN2_SQ = GAMMAN2 * GAMMAN2
        HELPN2  = GAMMAN2_SQ * SIGMA_4
        WAVELENGTHS_SHIFTED(S) = 1.0D+7 / SIGMA
        DO N = 1, NLAYERS
          RRSXSEC_SHIFTED   (N,S) = RRS_N2REF   (N,J) * HELPN2
          LINRRSXSEC_SHIFTED(N,S) = LINRRS_N2REF(N,J) * HELPN2
        ENDDO
      ENDDO

!  Shifted wavelengths and cross sections (O2)

      DO J = 1, N_RRS_O2REF
        S = S + 1
        SIGMA = SIGMA_0 + O2POS(J)
        LAMDA = 1.0d+07 / SIGMA
        SIGSQ = 1.0d+06 / LAMDA / LAMDA
        GAMMAO2 = GAMMA_O2(1) + GAMMA_O2(2) / (GAMMA_O2(3) - SIGSQ)
        GAMMAO2_SQ = GAMMAO2 * GAMMAO2
        HELPO2  = GAMMAO2_SQ * SIGMA_4
        WAVELENGTHS_SHIFTED(S) = 1.0D+7 / SIGMA
        DO N = 1, NLAYERS
          RRSXSEC_SHIFTED   (N,S) = RRS_O2REF   (N,J) * HELPO2
          LINRRSXSEC_SHIFTED(N,S) = LINRRS_O2REF(N,J) * HELPO2
        ENDDO
      ENDDO

!  Total number of shifts

      NSHIFTS      = N_RRS_N2REF + N_RRS_O2REF

!   we need to multiply by d-48 to get the correct
!   units for the absolute cross sections in cm2: We calculate gamman2
!   and gammao2 in d-30 m3 and then calculate n2xsec and o2xsec values
!   as (gamma**2) * (sigprime**4), where sigprime is in cm-1.

      DO S = 1, NSHIFTS
        DO N = 1, NLAYERS
          RRSXSEC_SHIFTED   (N,S) = RRSXSEC_SHIFTED   (N,S) * D24 * D24
          LINRRSXSEC_SHIFTED(N,S) = LINRRSXSEC_SHIFTED(N,S) * D24 * D24
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LINRRS_GAINTERM_XSEC_MONO

!  End Module

      END MODULE lrrs_L_raman_spectroscopy_m

