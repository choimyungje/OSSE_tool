
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
! #            LRRS_BUFFERING_BIN                               #
! #            LRRS_BUFFERING_MONO                              #
! #            LRRS_RAMAN_XSECS_BIN                             #
! #            LRRS_RAMAN_XSECS_MONO                            #
! #                                                             #
! #   Other SUBROUTINES called by the RAMAN_XSEC modules:       #
! #                                                             #
! #                RRS_BASIC_SPECTROSCOPY   (Public)            #
! #                RRS_LOSSTERM_XSEC        (Private)           #
! #                RRS_GAINTERM_XSEC_BIN    (Private)           #
! #                RRS_GAINTERM_XSEC_MONO   (Private)           #
! #                                                             #
! #  Utility Routines (Public)                                  #
! #                                                             #
! #            RRS_LAMDASRANKED_BIN                             #
! #            RRS_LAMDASRANKED_MONO                            #
! #                                                             #
! #            RRSBIN_REMAPPING                                 #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Bookkeeping improvements (use of "Only", clearer I/O specifications)

      MODULE lrrs_raman_spectroscopy_m

      USE LRRS_PARS_m, Only : FPK, zero, One, half, two, PIE, MAX_BINS, MAX_MESSAGES, SDU

      PRIVATE :: RRS_LOSSTERM_XSEC,      &
                 RRS_GAINTERM_XSEC_BIN,  &
                 RRS_GAINTERM_XSEC_MONO

      PUBLIC  :: LRRS_BUFFERING_BIN,     &
                 LRRS_BUFFERING_MONO,    &
                 LRRS_RAMAN_XSECS_BIN,   &
                 LRRS_RAMAN_XSECS_MONO,  &
                 RRS_BASIC_SPECTROSCOPY, &
                 RRS_LAMDASRANKED_BIN,   &
                 RRS_LAMDASRANKED_MONO,  &
                 RRSBIN_REMAPPING

      CONTAINS

      SUBROUTINE LRRS_BUFFERING_BIN &
        ( MAX_INPUT_LAMBDAS, MAX_POINTS, CRIT_2, DO_PURE_ELASTIC, & ! Inputs
          N_INPUT_LAMBDAS, INPUT_LAMBDAS, INPUT_FLUXES,           & ! Inputs
          LAMBDA_START, LAMBDA_FINISH,                            & ! Inputs
          NPOINTS_INNER, NPOINTS_TOTAL, OFFSET_INNER,             & ! Outputs
          LAMBDAS, FLUXES, BINLOWER, BINUPPER,                    & ! Outputs
          FAIL, N_MESSAGES, MESSAGES )                              ! Outputs

!  Stand-alone routine, no include files

!  Perform buffering
!  -----------------

!   given the starting and ending wavelengths, and
!   given a set of wavelengths and fluxes from data, then:
!      - Set the inner and outer windows and the inner window offset
!      - Set the wavelengths and fluxes on the outer window
!      - Get the binning limits for the outer window

      IMPLICIT NONE

!  input arguments
!  ---------------

!  dimensions

      INTEGER  , INTENT(IN) :: max_input_lambdas
      INTEGER  , INTENT(IN) :: max_points

!  Critical step at the end

      REAL(FPK), INTENT(IN) :: crit_2

!  pure elastic flag

      LOGICAL  , INTENT(IN) :: do_pure_elastic

!  Data input of wavelengths and fluxes

      INTEGER  , INTENT(IN) :: n_input_lambdas
      REAL(FPK), INTENT(IN) :: input_lambdas ( max_input_lambdas )
      REAL(FPK), INTENT(IN) :: input_fluxes  ( max_input_lambdas )

!  Starting and ending wavelengths of inner window

      REAL(FPK), INTENT(IN) :: lambda_start
      REAL(FPK), INTENT(IN) :: lambda_finish

!  output arguments
!  ----------------

!  inner and outer buffer sizes + inner window offset

      INTEGER  , INTENT(OUT) :: npoints_inner
      INTEGER  , INTENT(OUT) :: npoints_total
      INTEGER  , INTENT(OUT) :: offset_inner

!  Wavelengths and fluxes defined on the OUTER window

      REAL(FPK), INTENT(INOUT) :: lambdas  ( max_points )
      REAL(FPK), INTENT(INOUT) :: fluxes   ( max_points )

!  Binning limits defined for the OUTER window

      REAL(FPK), INTENT(INOUT) :: binlower ( max_points )
      REAL(FPK), INTENT(INOUT) :: binupper ( max_points )

!  status flag

      LOGICAL          , INTENT(OUT) ::  fail
      INTEGER          , INTENT(OUT) ::  n_messages
      CHARACTER (LEN=*), INTENT(OUT) :: messages ( max_messages )

!  local variables
!  ---------------

      INTEGER :: n_inner,     n_outer
      INTEGER :: imark_outer, imark_inner
      INTEGER :: ioff_outer,  ioff_inner

      REAL(FPK) :: sig_0, sig_t, num, nup
      REAL(FPK) :: lambda_start_shifted
      REAL(FPK) :: lambda_finish_shifted
      REAL(FPK) :: lam1, lam2, lam3
      LOGICAL :: loop_outer, loop_inner
      INTEGER :: q, i, nip

!  Start of code
!  -------------

!  set status

      fail = .false.
      n_messages = 0
      messages   = ' '

      nip = n_input_lambdas

!  Go to simpler code for pure elastic calculation

      if ( do_pure_elastic ) go to 6666

!  find extreme wavelength limits, determined by the highest Raman
!   Stokes and Anti-Stokes shifts (194.30 and 196.83 wavenumbers)

      num = - 194.3015d0
      sig_0  = 1.0d+7 / lambda_start
      sig_t  = sig_0 - num
      lambda_start_shifted = 1.0d+7 / sig_t
      lambda_start_shifted = lambda_start_shifted - crit_2

      nup = + 196.8269d0
      sig_0  = 1.0d+7 / lambda_finish
      sig_t  = sig_0 - nup
      lambda_finish_shifted = 1.0d+7 / sig_t
      lambda_finish_shifted = lambda_finish_shifted + crit_2

!  find range of calculation points

      loop_outer = .true.
      loop_inner = .true.
      imark_outer = 0
      imark_inner = 0
      i = 0
      n_inner = 0
      n_outer = 0
      do while ( loop_outer .or. loop_inner )
        i = i + 1
        if ( input_lambdas(i) .ge. lambda_start_shifted ) then
          if ( input_lambdas(i) .le. lambda_finish_shifted ) then
            n_outer = n_outer + 1
            if ( input_lambdas(i) .ge. lambda_start ) then
              if ( input_lambdas(i) .le. lambda_finish ) then
                n_inner = n_inner + 1
              else
                if ( loop_inner ) imark_inner = i - 1
                loop_inner = .false.
              endif
            endif
          else
            if ( loop_outer ) imark_outer = i - 1
            loop_outer  = .false.
          endif
        endif
      enddo
      ioff_outer = imark_outer - n_outer
      ioff_inner = imark_inner - n_inner

!  get overall size of wavelength grid

      npoints_inner = n_inner
      offset_inner  = ioff_inner - ioff_outer

!  Reset outer bin number and offset
!    ( Check also that total number of points does not exceed limits)

      ioff_outer    = ioff_outer - 1
      npoints_total = n_outer + 2
      offset_inner  = offset_inner + 1

!mick mod 6/6/2011 - additional buffer error handling
      if ( ( offset_inner*2 + 1 .gt. max_bins ) .or. &
           ( npoints_total .gt. max_points ) ) then
        fail=.true.

        !For reference:
        !write(*,*) 'offset_inner   = ',offset_inner
        !write(*,*) '# of RRS bins  = ',offset_inner*2 + 1
        !write(*,*) 'npoints_inner  = ',npoints_inner
        !write(*,*) 'npoints_total  = ',npoints_total
        !write(*,*) 'max_bins       = ',max_bins
        !write(*,*) 'max_points     = ',max_points

        N_MESSAGES = N_MESSAGES + 1
        write(MESSAGES(N_MESSAGES),'(a)') &
          '  Some Raman spectroscopy array dimensions inadequate'
        N_MESSAGES = N_MESSAGES + 1
        write(MESSAGES(N_MESSAGES),'(a,i4)') &
          '  # of RRS bins  = ',offset_inner*2 + 1
        N_MESSAGES = N_MESSAGES + 1
        write(MESSAGES(N_MESSAGES),'(a,i4)') &
          '  npoints_inner  = ',npoints_inner
        N_MESSAGES = N_MESSAGES + 1
        write(MESSAGES(N_MESSAGES),'(a,i4)') &
          '  max_bins       = ',max_bins
        N_MESSAGES = N_MESSAGES + 1
        write(MESSAGES(N_MESSAGES),'(a,i4)') &
          '  max_points     = ',max_points
        N_MESSAGES = N_MESSAGES + 1
        write(MESSAGES(N_MESSAGES),'(a)') &
          '  Inside file lrrs_pars.f90_NPTS_100, insure:'
        N_MESSAGES = N_MESSAGES + 1
        write(MESSAGES(N_MESSAGES),'(a)') &
          '    max_bins   > # of RRS bins'
        N_MESSAGES = N_MESSAGES + 1
        write(MESSAGES(N_MESSAGES),'(a)') &
          '    max_points > # of RRS bins + npoints_inner'

        return
      endif

!  Set bin limits

      lam1 = input_lambdas(ioff_outer)
      lam2 = input_lambdas(ioff_outer+1)
      do q = 1, npoints_total
        lam3 = input_lambdas(q+1+ioff_outer)
        lambdas(q)  = lam2
        fluxes(q)   = input_fluxes(q+ioff_outer)
        binlower(q) = half * (lam1 + lam2 )
        binupper(q) = half * (lam3 + lam2 )
        lam1 = lam2
        lam2 = lam3
      enddo

!  debug
!      do q = 1, npoints_total
!        write(*,'(i5,3f12.5,1pe15.5)')
!     &     q,lambdas(q),binlower(q),binupper(q),fluxes(q)
!      enddo
!      do i = 1, npoints_inner
!        q = i + offset_inner
!        write(*,'(i5,3f12.5,1pe15.5)')
!     &     q,lambdas(q),binlower(q),binupper(q),fluxes(q)
!      enddo

!  Normal return

      return

!  continuation point for pure elastic binning

6666  continue

!  find range of calculation points

      loop_inner = .true.
      i = 0
      n_inner = 0
      imark_inner = 0
      do while ( loop_inner )
        i = i + 1
        if ( input_lambdas(i) .ge. lambda_start ) then
          if ( input_lambdas(i) .le. lambda_finish ) then
            n_inner = n_inner + 1
          else
            if ( loop_inner ) imark_inner = i - 1
            loop_inner = .false.
          endif
        endif
      enddo
      ioff_inner = imark_inner - n_inner

!  get overall size of wavelength grid

      npoints_inner = n_inner
      offset_inner  = 0
      npoints_total = n_inner

!  Set bin limits

      lam1 = input_lambdas(ioff_inner)
      lam2 = input_lambdas(ioff_inner+1)
      do q = 1, npoints_total
        lam3 = input_lambdas(q+1+ioff_inner)
        lambdas(q)  = lam2
        fluxes(q)   = input_fluxes(q+ioff_inner)
        binlower(q) = half * (lam1 + lam2 )
        binupper(q) = half * (lam3 + lam2 )
        lam1 = lam2
        lam2 = lam3
      enddo

!  finish

      return
      END SUBROUTINE LRRS_BUFFERING_BIN

!

      SUBROUTINE LRRS_BUFFERING_MONO &
           ( MAX_SUN_LAMBDAS, DO_PURE_ELASTIC, & ! Inputs
             SUN_LAMBDAS, SUN_FLUXES,          & ! Inputs
             LAMBDA_START, LAMBDA_FINISH,      & ! Inputs
             N_USED_LAMBDAS, USED_LAMBDAS, USED_FLUXES ) ! Outputs

!  Stand-alone routine, no include files

!  Perform buffering
!  -----------------

!   given the starting and ending wavelengths, and
!   given a set of wavelengths and fluxes from data, then:
!      - set an outer window for splining fluxes

      IMPLICIT NONE

!  input arguments
!  ---------------

!  dimensions

      INTEGER  , INTENT(IN) :: max_sun_lambdas

!  pure elastic flag

      LOGICAL  , INTENT(IN) :: do_pure_elastic

!  Data input of wavelengths and fluxes

      REAL(FPK), INTENT(IN) :: sun_lambdas ( max_sun_lambdas )
      REAL(FPK), INTENT(IN) :: sun_fluxes  ( max_sun_lambdas )

!  Starting and ending wavelengths of inner window

      REAL(FPK), INTENT(IN) :: lambda_start
      REAL(FPK), INTENT(IN) :: lambda_finish

!  output arguments
!  ----------------

!  Wavelengths and fluxes defined on the OUTER window

      INTEGER  , INTENT(OUT) :: n_used_lambdas
      REAL(FPK), INTENT(OUT) :: used_lambdas ( max_sun_lambdas )
      REAL(FPK), INTENT(OUT) :: used_fluxes  ( max_sun_lambdas )

!  local variables
!  ---------------

      REAL(FPK) :: sig_0, sig_t, num, nup
      REAL(FPK) :: lambda_start_shifted
      REAL(FPK) :: lambda_finish_shifted
      LOGICAL :: loop_outer
      INTEGER :: i, ic

!  Start of code
!  -------------

!  initialize output

      n_used_lambdas = 0 ; used_lambdas = zero ; used_fluxes = zero

!  set wavelengths for buffering. Elastic is trivial.
!     * find extreme wavelength limits, determined by the highest Raman
!       Stokes and Anti-Stokes shifts (194.30 and 196.83 wavenumbers)

      if ( do_pure_elastic ) then
         lambda_start_shifted  = lambda_start
         lambda_finish_shifted = lambda_finish
      else
         num = - 194.3015d0
         sig_0  = 1.0d+7 / lambda_start
         sig_t  = sig_0 - num
         lambda_start_shifted = 1.0d+7 / sig_t
         nup = + 196.8269d0
         sig_0  = 1.0d+7 / lambda_finish
         sig_t  = sig_0 - nup
         lambda_finish_shifted = 1.0d+7 / sig_t
      endif

!  find range of calculation points

      loop_outer = .true.
      ic = 0  ; i  = 0
      do while ( loop_outer )
        i = i + 1
        if ( sun_lambdas(i) .ge. lambda_start_shifted ) then
          if ( sun_lambdas(i) .le. lambda_finish_shifted ) then
            ic = ic + 1
            used_lambdas(ic) = sun_lambdas(i)
            used_fluxes(ic)  = sun_fluxes(i)
          else
            loop_outer  = .false.
          endif
        endif
      enddo
      n_used_lambdas = ic

!  finish

      return
      END SUBROUTINE LRRS_BUFFERING_MONO

!

      SUBROUTINE LRRS_RAMAN_XSECS_BIN &
        ( MAX_LAYERS, MAX_POINTS, MAX_BINS,                    & ! Inputs
          NLAYERS, NPOINTS_INNER, NPOINTS_OUTER, OFFSET_INNER, & ! Inputs
          BINLOWER, BINUPPER, LAMBDAS_RANKED, FLUXES_RANKED,   & ! Inputs
          N_RRS_N2REF, N_RRS_O2REF, RRS_N2REF, RRS_O2REF,      & ! Inputs
          N2POS, O2POS, GAMMA_N2, GAMMA_O2,                    & ! Inputs
          N_RRSBINS, BINMAP, BINXSEC, RRSXSEC_OUT, RINGSPEC_1, & ! Outputs
          FAIL, MESSAGE, ACTION )                                ! Outputs

       IMPLICIT NONE

!  Comments still required
!  =======================

!  This is stand alone.

!  Routine is called INSIDE the LRRS Master module; wavelength buffering
!    is done first, so that ranked LAMBDAS/FLUXES are inputs, along with
!    the inner and outer point numbers and the binning.

!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@
!  Exception handling for MAX_BINS Dimension

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

!  Chance Ring Spectrum

      REAL(FPK), INTENT(OUT) :: RINGSPEC_1 ( MAX_POINTS )

!@@@@@@@@@ RT Solutions, 14 April 2011 @@@@@@@@@
!  Exception handling for MAX_BINS Dimension

      LOGICAL          , INTENT(OUT)   :: FAIL
      CHARACTER (LEN=*), INTENT(OUT)   :: MESSAGE, ACTION

!  Local variables
!  ===============

!  1. For the RRS spectroscic arrays and associated quantitites
!  ------------------------------------------------------------

!  Shift numbers

      INTEGER  , PARAMETER :: MAX_SHIFTS = 233

!  number of shifts = N_Stokes + N_Antistokes

      INTEGER   :: NSHIFTS, N_STOKES, N_ANTISTOKES

!  Gain term cross-sections and wavelengths (shifted) of excitatoion)
!   (Includes separately the Stokes and AntiStokes cross-sections)

      REAL(FPK) :: LAMBDAS_SHIFTED    ( MAX_SHIFTS )
      REAL(FPK) :: RRSXSEC_SHIFTED    ( MAX_LAYERS, MAX_SHIFTS )
      REAL(FPK) :: RRSXSEC_STOKES     ( MAX_LAYERS, MAX_SHIFTS )
      REAL(FPK) :: RRSXSEC_ANTISTOKES ( MAX_LAYERS, MAX_SHIFTS )

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
         BINXSEC(N,K,B) = zero
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

        CALL RRS_GAINTERM_XSEC_BIN &
           ( MAX_LAYERS, MAX_SHIFTS, NLAYERS, N_RRS_N2REF, N_RRS_O2REF, & ! Inputs
             N2POS, O2POS, GAMMA_N2, GAMMA_O2,                          & ! Inputs
             LAMBDAS_RANKED(I), RRS_N2REF, RRS_O2REF,                   & ! Inputs
             N_STOKES,     RRSXSEC_STOKES,                              & ! Outputs
             N_ANTISTOKES, RRSXSEC_ANTISTOKES,                          & ! Outputs
             LAMBDAS_SHIFTED, RRSXSEC_SHIFTED )                           ! Outputs

!  debug

!        write(89,'(i5,1p233e16.6)')I,(RRSXSEC_SHIFTED(1,J),J=1,233)
!        write(90,*)'----------- ',I,LAMBDAS_RANKED(I)
!        do j = 1, 233
!          write(90,'(i5,f12.6)')J,LAMBDAS_SHIFTED(J)
!        enddo

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
!                ENDDO
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
         MESSAGE = ' Dimension MAX_BINS exceeded in Subroutine LRRS_RAMAN_XSECS_BIN'
         ACTION  = ' Increase MAX_BINS in lrrs_pars.f90 to at least '//c3
         RETURN
      ENDIF

!  debug and ring spec
!  -------------------

!      write(*,*)lambda_start, lambda_finish, sigma
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
      rms = 1.0e-10*dsqrt(rms/dble(npoints_inner))
      do k = 1, npoints_inner
       ringspec_1(k) = ringspec_1(k) / rms
      enddo

!  Loss term RRS
!  -------------

!  For each wavelength in the inner loop, compute the total cross-section
!  for RRS light scattered OUT OF THIS WAVELENGTH. This is the Loss term
!    Requires temperatures as input and the basic reference values

      DO WI = 1, NPOINTS_INNER
        WO = WI + OFFSET_INNER
        CALL RRS_LOSSTERM_XSEC &
            ( MAX_LAYERS, NLAYERS, N_RRS_N2REF, N_RRS_O2REF, & ! Inputs
              N2POS, O2POS, GAMMA_N2, GAMMA_O2,              & ! Inputs
              LAMBDAS_RANKED(WO), RRS_N2REF, RRS_O2REF,      & ! Inputs
              RRSXSEC_OUT(1,WI), NSHIFTS )                     ! Outputs
      enddo

!  Finish routine

      return
      END SUBROUTINE LRRS_RAMAN_XSECS_BIN

!

      SUBROUTINE LRRS_RAMAN_XSECS_MONO &
        ( MAX_LAYERS, MAX_POINTS, NLAYERS, LAMBDA_EXCIT,  & ! Inputs
          N_RRS_N2REF, N_RRS_O2REF, RRS_N2REF, RRS_O2REF, & ! Inputs
          N2POS, O2POS, GAMMA_N2, GAMMA_O2,               & ! Inputs
          RRSXSEC_OUT_MONO, RRSXSEC_RANKED )                ! Outputs

!  This is stand alone.

!  Routine is called Inside LRRS_Setups_Master module; the only inputs are
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

!  Positions

      REAL(FPK), INTENT(IN) :: O2POS (N_RRS_O2REF)
      REAL(FPK), INTENT(IN) :: N2POS (N_RRS_N2REF)

!  Gammas

      REAL(FPK), INTENT(IN) :: GAMMA_O2(3)
      REAL(FPK), INTENT(IN) :: GAMMA_N2(3)

!  Output arguments
!  ================

!  Gain term cross-sections

      REAL(FPK), intent(inout) :: RRSXSEC_RANKED ( MAX_LAYERS, MAX_POINTS )

!  Loss term cross sections

      REAL(FPK), intent(inout) :: RRSXSEC_OUT_MONO    ( MAX_LAYERS )

!  Local variables
!  ===============

!  Shift numbers

      INTEGER  , PARAMETER  :: MAX_SHIFTS   = 233
      INTEGER  , PARAMETER  :: MAX_SHIFTS_1 = 234

!  Number of points
!  number of shifts = N_Stokes + N_Antistokes

      INTEGER   :: NPOINTS_MONO
      INTEGER   :: NSHIFTS

!  Gain term cross-sections and wavelengths (shifted) of excitatoion)
!   (Includes separately the STokes and AntiStokes cross-sections)

      REAL(FPK) :: LAMBDAS_SHIFTED    ( MAX_SHIFTS )
      REAL(FPK) :: RRSXSEC_SHIFTED    ( MAX_LAYERS, MAX_SHIFTS )

!  Ranking

      INTEGER   :: INDEX_ALL          ( MAX_SHIFTS_1 )
      REAL(FPK) :: LAMBDAS_ALL        ( MAX_SHIFTS_1 )
      REAL(FPK) :: RRSXSEC_ALL        ( MAX_LAYERS, MAX_SHIFTS_1 )

!  Help variables

      INTEGER   :: I, W, N

!  Initial section
!  ===============

!  Zero the Ranked cross-sections

      RRSXSEC_OUT_MONO = zero
      RRSXSEC_RANKED   = zero

!  Loss term RRS
!  -------------

!  For each wavelength in the inner loop, compute the total cross-section
!  for RRS light scattered OUT OF THIS WAVELENGTH. This is the Loss term
!    Requires temperatures as input and the basic reference values

      CALL RRS_LOSSTERM_XSEC &
            ( MAX_LAYERS, NLAYERS, N_RRS_N2REF, N_RRS_O2REF, & ! Inputs
              N2POS, O2POS, GAMMA_N2, GAMMA_O2,              & ! Inputs
              LAMBDA_EXCIT, RRS_N2REF, RRS_O2REF,            & ! Inputs
              RRSXSEC_OUT_MONO, NSHIFTS )                      ! Output

!  Gain term RRS (Mono)
!  --------------------

!  find shifted wavelengths and cross sections for the RRS Gain terms =
!    ineslatic scattering of light INTO the wavelength LAMBDAS(wo)

      CALL RRS_GAINTERM_XSEC_MONO &
           ( MAX_LAYERS, MAX_SHIFTS, NLAYERS, N_RRS_N2REF, N_RRS_O2REF, & ! Inputs
             N2POS, O2POS, GAMMA_N2, GAMMA_O2,                          & ! Inputs
             LAMBDA_EXCIT,  RRS_N2REF, RRS_O2REF,                       & ! Inputs
             LAMBDAS_SHIFTED, RRSXSEC_SHIFTED )                           ! Output

!  Ranking procedure
!  =================

!  Copy to holding arrays and add Excitation wavelength
!  Number of points is always the same

      NPOINTS_MONO = MAX_SHIFTS_1
      DO W = 1, NSHIFTS
        LAMBDAS_ALL(W) = LAMBDAS_SHIFTED(W)
        DO N = 1, NLAYERS
          RRSXSEC_ALL(N,W) = RRSXSEC_SHIFTED(N,W)
        ENDDO
      ENDDO
      LAMBDAS_ALL(NPOINTS_MONO) = LAMBDA_EXCIT
      DO N = 1, NLAYERS
        RRSXSEC_ALL(N,NPOINTS_MONO) = zero
      ENDDO

!  Use RQSORT_IDX to rank the cross sections

      DO I=1,NPOINTS_MONO
        INDEX_ALL(I) = I
      ENDDO
      CALL RQSORT_IDX(NPOINTS_MONO, LAMBDAS_ALL, NPOINTS_MONO, INDEX_ALL)
      DO W = 1, NPOINTS_MONO
        DO N = 1, NLAYERS
          RRSXSEC_RANKED(N,W) = RRSXSEC_ALL(N,INDEX_ALL(W))
        ENDDO
      ENDDO

!  Finish routine

      return
      END SUBROUTINE LRRS_RAMAN_XSECS_MONO

!

      SUBROUTINE RRS_BASIC_SPECTROSCOPY &
         ( MAXLAYERS, NLAYERS, DO_TEMPERATURE_WF,             & ! Inputs
           TEMPERATURES, M_RRS_N2REF, M_RRS_O2REF,            & ! Inputs
           N2POS, O2POS, GAMMA_N2, GAMMA_O2,                  & ! Outputs
           RRS_N2REF, RRS_O2REF, LINRRS_N2REF, LINRRS_O2REF )   ! Outputs

!  Purpose
!  -------

!  Get the basic transition coefficients for RRS:
!     Spin degeneracy factors
!     Partition functions (temperature dependent)
!     Placzek-Teller Coefficients
!     Mixing-ratios, 256.pi^5/27
!   ... all these multiplied together

!  Version 2.2M.
!   *** Will also generate optional T-derivatives

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Temperature derivative flag

      LOGICAL  , INTENT(IN) :: DO_TEMPERATURE_WF

!  shift dimensions

      INTEGER  , INTENT(IN) :: M_RRS_N2REF
      INTEGER  , INTENT(IN) :: M_RRS_O2REF

!  Layering and local layer temperatures (one per layer)

      INTEGER  , INTENT(IN) :: MAXLAYERS, NLAYERS
      REAL(FPK), INTENT(IN) :: TEMPERATURES ( MAXLAYERS )

!  Output arguments
!  ----------------

!  Positions

      REAL(FPK), INTENT(OUT) :: O2POS (M_RRS_O2REF)
      REAL(FPK), INTENT(OUT) :: N2POS (M_RRS_N2REF)

!  Gammas

      REAL(FPK), INTENT(OUT) :: GAMMA_O2(3)
      REAL(FPK), INTENT(OUT) :: GAMMA_N2(3)

!  Basic coefficients for O2 and N2

      REAL(FPK), INTENT(OUT) :: RRS_N2REF ( MAXLAYERS, M_RRS_N2REF )
      REAL(FPK), INTENT(OUT) :: RRS_O2REF ( MAXLAYERS, M_RRS_O2REF )

!  T-derivatives of Basic coefficients for O2 and N2

      REAL(FPK), INTENT(OUT) :: LINRRS_N2REF ( MAXLAYERS, M_RRS_N2REF )
      REAL(FPK), INTENT(OUT) :: LINRRS_O2REF ( MAXLAYERS, M_RRS_O2REF )

!  Local variables

      REAL(FPK) :: C2, C0, EMULT, TEMP, QN2, QO2, QN2_I, QO2_I
      REAL(FPK) :: D_EMULT, D_QN2, D_QO2, D_LOGQN2, D_LOGQO2
      REAL(FPK) :: N2FRAC, O2FRAC, PREFIX_N2, PREFIX_O2, PCPREF
      REAL(FPK) :: D_N2FRAC, D_O2FRAC, D_TERM, D_QN2_I, D_QO2_I
      INTEGER   :: I, N

!  Spectroscopy stuff
!  ------------------

!  single numbers

      REAL(FPK) :: GAMMA_N2_A
      REAL(FPK) :: GAMMA_N2_B
      REAL(FPK) :: GAMMA_N2_C
      REAL(FPK) :: MIXRATIO_N2

      REAL(FPK) :: GAMMA_O2_A
      REAL(FPK) :: GAMMA_O2_B
      REAL(FPK) :: GAMMA_O2_C
      REAL(FPK) :: MIXRATIO_O2

!     Statistical parameters, for partition function calculation

      REAL(FPK) :: N2STAT_1 (31)
      REAL(FPK) :: N2STAT_2 (31)
      REAL(FPK) :: N2STAT_3 (31)
      REAL(FPK) :: O2STAT_1 (54)
      REAL(FPK) :: O2STAT_2 (54)

!     Rotational Raman line-specific parameters

      REAL(FPK) :: N2TERM (48)
      REAL(FPK) :: N2PLACTEL (48)
      REAL(FPK) :: N2DEG (48)
      REAL(FPK) :: N2NUC (48)
      REAL(FPK) :: RN2POS (48)

      REAL(FPK) :: O2TERM (185)
      REAL(FPK) :: O2PLACTEL (185)
!      REAL(FPK) :: O2PLACTEL_KC99 (185)
      REAL(FPK) :: O2DEG (185)
      REAL(FPK) :: RO2POS (185)

!mick fix 7/20/2016 - initialized all data below to double precision (i.e. "fpk")

!  N2 spectroscopy
!  ===============

!  Gamma constants

      DATA  GAMMA_N2_A / -0.601466_fpk /
      DATA  GAMMA_N2_B /   238.557_fpk /
      DATA  GAMMA_N2_C /   186.099_fpk /

!  VMR constants

      DATA  MIXRATIO_N2 /   0.7905_fpk /
!      DATA  MIXRATIO_N2 / 0.7809  /

!  number of transitions
!      DATA  N_RRS_N2REF / 48  /

!  Transitions (wavenumber shifts)

      DATA RN2POS / &
        -194.3015_fpk, -186.4226_fpk, -178.5374_fpk, -170.6459_fpk, -162.7484_fpk, &
        -154.8453_fpk, -146.9368_fpk, -139.0233_fpk, -131.1049_fpk, -123.1819_fpk, &
        -115.2547_fpk, -107.3235_fpk,  -99.3886_fpk,  -91.4502_fpk,  -83.5086_fpk, &
         -75.5642_fpk,  -67.6171_fpk,  -59.6676_fpk,  -51.7162_fpk,  -43.7629_fpk, &
         -35.8080_fpk,  -27.8519_fpk,  -19.8950_fpk,  -11.9373_fpk,   11.9373_fpk, &
          19.8950_fpk,   27.8519_fpk,   35.8080_fpk,   43.7629_fpk,   51.7162_fpk, &
          59.6676_fpk,   67.6171_fpk,   75.5642_fpk,   83.5086_fpk,   91.4502_fpk, &
          99.3886_fpk,  107.3235_fpk,  115.2547_fpk,  123.1819_fpk,  131.1049_fpk, &
         139.0233_fpk,  146.9368_fpk,  154.8453_fpk,  162.7484_fpk,  170.6459_fpk, &
         178.5374_fpk,  186.4226_fpk,  194.3015_fpk /

!  N2 statistics

      DATA N2STAT_1 / &
         0._fpk,  1._fpk,  2._fpk,  3._fpk,  4._fpk,  5._fpk,  6._fpk,  7._fpk,  8._fpk,  9._fpk, &
        10._fpk, 11._fpk, 12._fpk, 13._fpk, 14._fpk, 15._fpk, 16._fpk, 17._fpk, 18._fpk, 19._fpk, &
        20._fpk, 21._fpk, 22._fpk, 23._fpk, 24._fpk, 25._fpk, 26._fpk, 27._fpk, 28._fpk, 29._fpk, &
        30._fpk /

      DATA N2STAT_2 / &
        6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, &
        6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, &
        6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, &
        6._fpk /

      DATA N2STAT_3 / &
           0.0000_fpk,    3.9791_fpk,   11.9373_fpk,   23.8741_fpk,   39.7892_fpk, &
          59.6821_fpk,   83.5521_fpk,  111.3983_fpk,  143.2197_fpk,  179.0154_fpk, &
         218.7839_fpk,  262.5240_fpk,  310.2341_fpk,  361.9126_fpk,  417.5576_fpk, &
         477.1673_fpk,  540.7395_fpk,  608.2722_fpk,  679.7628_fpk,  755.2090_fpk, &
         834.6081_fpk,  917.9574_fpk, 1005.2540_fpk, 1096.4948_fpk, 1191.6766_fpk, &
        1290.7963_fpk, 1393.8503_fpk, 1500.8350_fpk, 1611.7467_fpk, 1726.5816_fpk, &
        1845.3358_fpk /

!  N2 Term Energies

      DATA N2TERM / &
        1290.7963_fpk, 1191.6766_fpk, 1096.4948_fpk, 1005.2540_fpk,  917.9574_fpk,  834.6081_fpk, &
         755.2090_fpk,  679.7628_fpk,  608.2722_fpk,  540.7395_fpk,  477.1673_fpk,  417.5576_fpk, &
         361.9126_fpk,  310.2341_fpk,  262.5240_fpk,  218.7839_fpk,  179.0154_fpk,  143.2197_fpk, &
         111.3983_fpk,   83.5521_fpk,   59.6821_fpk,   39.7892_fpk,   23.8741_fpk,   11.9373_fpk, &
           0.0000_fpk,    3.9791_fpk,   11.9373_fpk,   23.8741_fpk,   39.7892_fpk,   59.6821_fpk, &
          83.5521_fpk,  111.3983_fpk,  143.2197_fpk,  179.0154_fpk,  218.7839_fpk,  262.5240_fpk, &
         310.2341_fpk,  361.9126_fpk,  417.5576_fpk,  477.1673_fpk,  540.7395_fpk,  608.2722_fpk, &
         679.7628_fpk,  755.2090_fpk,  834.6081_fpk,  917.9574_fpk, 1005.2540_fpk, 1096.4948_fpk /

!  N2 Placzek-Teller coefficiens

      DATA N2PLACTEL / &
        3.601e-1_fpk, 3.595e-1_fpk, 3.589e-1_fpk, 3.581e-1_fpk, 3.573e-1_fpk, &
        3.565e-1_fpk, 3.555e-1_fpk, 3.544e-1_fpk, 3.532e-1_fpk, 3.519e-1_fpk, &
        3.504e-1_fpk, 3.487e-1_fpk, 3.467e-1_fpk, 3.443e-1_fpk, 3.416e-1_fpk, &
        3.383e-1_fpk, 3.344e-1_fpk, 3.294e-1_fpk, 3.231e-1_fpk, 3.147e-1_fpk, &
        3.030e-1_fpk, 2.857e-1_fpk, 2.571e-1_fpk, 2.000e-1_fpk, 1.000e+0_fpk, &
        6.000e-1_fpk, 5.143e-1_fpk, 4.762e-1_fpk, 4.545e-1_fpk, 4.406e-1_fpk, &
        4.308e-1_fpk, 4.235e-1_fpk, 4.180e-1_fpk, 4.135e-1_fpk, 4.099e-1_fpk, &
        4.070e-1_fpk, 4.044e-1_fpk, 4.023e-1_fpk, 4.004e-1_fpk, 3.988e-1_fpk, &
        3.974e-1_fpk, 3.961e-1_fpk, 3.950e-1_fpk, 3.940e-1_fpk, 3.931e-1_fpk, &
        3.922e-1_fpk, 3.915e-1_fpk, 3.908e-1_fpk /

!  N2 degeneracy labels

      DATA N2DEG / &
        25._fpk, 24._fpk, 23._fpk, 22._fpk, 21._fpk, 20._fpk, &
        19._fpk, 18._fpk, 17._fpk, 16._fpk, 15._fpk, 14._fpk, &
        13._fpk, 12._fpk, 11._fpk, 10._fpk,  9._fpk,  8._fpk, &
         7._fpk,  6._fpk,  5._fpk,  4._fpk,  3._fpk,  2._fpk, &
         0._fpk,  1._fpk,  2._fpk,  3._fpk,  4._fpk,  5._fpk, &
         6._fpk,  7._fpk,  8._fpk,  9._fpk, 10._fpk, 11._fpk, &
        12._fpk, 13._fpk, 14._fpk, 15._fpk, 16._fpk, 17._fpk, &
        18._fpk, 19._fpk, 20._fpk, 21._fpk, 22._fpk, 23._fpk /

!  N2 Nuclear spin degeneracy (?)

      DATA N2NUC / &
        3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, &
        3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, &
        3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, &
        3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, &
        6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, &
        6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, &
        6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk, &
        6._fpk, 3._fpk, 6._fpk, 3._fpk, 6._fpk, 3._fpk /

!  O2 spectroscopy
!  ===============

!  Gamma constants

      DATA  GAMMA_O2_A /  0.07149_fpk /
      DATA  GAMMA_O2_B /  45.9364_fpk /
      DATA  GAMMA_O2_C /  48.2716_fpk /

!  VMR constants

      DATA  MIXRATIO_O2 /  0.2095_fpk /

!  number of transitions
!      DATA  N_RRS_O2REF / 185  /

!  Transitions

      DATA RO2POS / &
        -185.5861_fpk, -185.5690_fpk, -185.5512_fpk, -174.3154_fpk, -174.2980_fpk, &
        -174.2802_fpk, -163.0162_fpk, -162.9980_fpk, -162.9809_fpk, -151.6906_fpk, &
        -151.6730_fpk, -151.6551_fpk, -140.3405_fpk, -140.3230_fpk, -140.3046_fpk, &
        -128.9677_fpk, -128.9490_fpk, -128.9314_fpk, -117.5740_fpk, -117.5560_fpk, &
        -117.5372_fpk, -108.2632_fpk, -106.1613_fpk, -106.1430_fpk, -106.1240_fpk, &
        -104.3008_fpk,  -96.8136_fpk,  -94.7318_fpk,  -94.7120_fpk,  -94.6934_fpk, &
         -92.8515_fpk,  -85.3489_fpk,  -83.2871_fpk,  -83.2671_fpk,  -83.2473_fpk, &
         -81.3868_fpk,  -73.8694_fpk,  -71.8294_fpk,  -71.8080_fpk,  -71.7876_fpk, &
         -69.9077_fpk,  -62.3771_fpk,  -60.3611_fpk,  -60.3373_fpk,  -60.3157_fpk, &
         -58.4156_fpk,  -50.8730_fpk,  -48.8851_fpk,  -48.8571_fpk,  -48.8332_fpk, &
         -46.9116_fpk,  -40.1158_fpk,  -39.3565_fpk,  -37.4069_fpk,  -37.3688_fpk, &
         -37.3406_fpk,  -35.3953_fpk,  -27.8241_fpk,  -25.9472_fpk,  -25.8745_fpk, &
         -25.8364_fpk,  -25.8125_fpk,  -23.8629_fpk,  -16.2529_fpk,  -16.2529_fpk, &
         -14.3761_fpk,  -14.3033_fpk,  -14.1686_fpk,  -12.2918_fpk,   -2.1392_fpk, &
          -2.1202_fpk,   -2.1016_fpk,   -2.0843_fpk,   -2.0843_fpk,   -2.0818_fpk, &
          -2.0614_fpk,   -2.0398_fpk,   -2.0159_fpk,   -2.0116_fpk,   -1.9877_fpk, &
          -1.9735_fpk,   -1.9496_fpk,   -1.9455_fpk,   -1.9217_fpk,   -1.9003_fpk, &
          -1.8803_fpk,   -1.8768_fpk,   -1.8605_fpk,   -1.8422_fpk,   -0.1347_fpk, &
           0.1347_fpk,    1.8422_fpk,    1.8605_fpk,    1.8768_fpk,    1.8803_fpk, &
           1.9003_fpk,    1.9217_fpk,    1.9455_fpk,    1.9496_fpk,    1.9735_fpk, &
           1.9877_fpk,    2.0116_fpk,    2.0159_fpk,    2.0398_fpk,    2.0614_fpk, &
           2.0818_fpk,    2.0843_fpk,    2.0843_fpk,    2.1016_fpk,    2.1202_fpk, &
           2.1392_fpk,   12.2918_fpk,   14.1686_fpk,   14.3033_fpk,   14.3761_fpk, &
          16.2529_fpk,   16.2529_fpk,   23.8629_fpk,   25.8125_fpk,   25.8364_fpk, &
          25.8745_fpk,   25.9472_fpk,   27.8241_fpk,   35.3953_fpk,   37.3406_fpk, &
          37.3688_fpk,   37.4069_fpk,   39.3565_fpk,   40.1158_fpk,   46.9116_fpk, &
          48.8332_fpk,   48.8571_fpk,   48.8851_fpk,   50.8730_fpk,   58.4156_fpk, &
          60.3157_fpk,   60.3373_fpk,   60.3611_fpk,   62.3771_fpk,   69.9077_fpk, &
          71.7876_fpk,   71.8080_fpk,   71.8294_fpk,   73.8694_fpk,   81.3868_fpk, &
          83.2473_fpk,   83.2671_fpk,   83.2871_fpk,   85.3489_fpk,   92.8515_fpk, &
          94.6934_fpk,   94.7120_fpk,   94.7318_fpk,   96.8136_fpk,  104.3008_fpk, &
         106.1240_fpk,  106.1430_fpk,  106.1613_fpk,  108.2632_fpk,  115.7318_fpk, &
         117.5372_fpk,  117.5560_fpk,  117.5740_fpk,  119.6952_fpk,  128.9314_fpk, &
         128.9490_fpk,  128.9677_fpk,  140.3046_fpk,  140.3230_fpk,  140.3405_fpk, &
         151.6551_fpk,  151.6730_fpk,  151.6906_fpk,  162.9809_fpk,  162.9980_fpk, &
         163.0162_fpk,  174.2802_fpk,  174.2980_fpk,  174.3154_fpk,  185.5512_fpk, &
         185.5690_fpk,  185.5861_fpk,  196.7919_fpk,  196.8100_fpk,  196.8269_fpk /

!  O2 statistics

      DATA O2STAT_1 / &
         0._fpk,  2._fpk,  1._fpk,  2._fpk,  4._fpk,  3._fpk, &
         4._fpk,  6._fpk,  5._fpk,  8._fpk,  6._fpk,  7._fpk, &
        10._fpk,  8._fpk,  9._fpk, 12._fpk, 10._fpk, 11._fpk, &
        14._fpk, 12._fpk, 13._fpk, 16._fpk, 14._fpk, 15._fpk, &
        18._fpk, 16._fpk, 17._fpk, 20._fpk, 18._fpk, 19._fpk, &
        22._fpk, 20._fpk, 21._fpk, 24._fpk, 22._fpk, 23._fpk, &
        26._fpk, 24._fpk, 25._fpk, 28._fpk, 26._fpk, 27._fpk, &
        30._fpk, 28._fpk, 29._fpk, 32._fpk, 30._fpk, 31._fpk, &
        34._fpk, 32._fpk, 33._fpk, 36._fpk, 34._fpk, 35._fpk /

      DATA O2STAT_2 / &
          0.0000_fpk,    2.0843_fpk,    3.9611_fpk,   16.2529_fpk,   16.3876_fpk,   18.3372_fpk, &
         42.2001_fpk,   42.2240_fpk,   44.2117_fpk,   79.5646_fpk,   79.6070_fpk,   81.5805_fpk, &
        128.3978_fpk,  128.4921_fpk,  130.4376_fpk,  188.7135_fpk,  188.8532_fpk,  190.7749_fpk, &
        260.5011_fpk,  260.6826_fpk,  262.5829_fpk,  343.7484_fpk,  343.9697_fpk,  345.8500_fpk, &
        438.4418_fpk,  438.7015_fpk,  440.5620_fpk,  544.5658_fpk,  544.8628_fpk,  546.7050_fpk, &
        662.1030_fpk,  662.4368_fpk,  664.2610_fpk,  791.0344_fpk,  791.4045_fpk,  793.2100_fpk, &
        931.3390_fpk,  931.7450_fpk,  933.5330_fpk, 1082.9941_fpk, 1083.4356_fpk, 1085.2060_fpk, & 
       1245.9750_fpk, 1246.4518_fpk, 1248.2040_fpk, 1420.2552_fpk, 1420.7672_fpk, 1422.5020_fpk, &
       1605.8064_fpk, 1606.3533_fpk, 1608.0710_fpk, 1802.5983_fpk, 1803.1802_fpk, 1804.8810_fpk /

!  O2 Term Energies

      DATA O2TERM / &
        1606.3533_fpk, 1608.0710_fpk, 1605.8064_fpk, 1420.7672_fpk, 1422.5020_fpk, &
        1420.2552_fpk, 1246.4518_fpk, 1248.2040_fpk, 1245.9750_fpk, 1083.4356_fpk, &
        1085.2060_fpk, 1082.9941_fpk,  931.7450_fpk,  933.5330_fpk,  931.3390_fpk, &
         791.4045_fpk,  793.2100_fpk,  791.0344_fpk,  662.4368_fpk,  664.2610_fpk, &
         662.1030_fpk,  546.7050_fpk,  544.8628_fpk,  546.7050_fpk,  544.5658_fpk, &
         544.8628_fpk,  440.5620_fpk,  438.7015_fpk,  440.5620_fpk,  438.4418_fpk, &
         438.7015_fpk,  345.8500_fpk,  343.9697_fpk,  345.8500_fpk,  343.7484_fpk, &
         343.9697_fpk,  262.5829_fpk,  260.6826_fpk,  262.5829_fpk,  260.5011_fpk, &
         260.6826_fpk,  190.7749_fpk,  188.8532_fpk,  190.7749_fpk,  188.7135_fpk, &
         188.8532_fpk,  130.4376_fpk,  128.4921_fpk,  130.4376_fpk,  128.3978_fpk, &
         128.4921_fpk,   42.2001_fpk,   81.5805_fpk,   79.6070_fpk,   81.5805_fpk, &
          79.5646_fpk,   79.6070_fpk,   44.2117_fpk,   42.2001_fpk,   44.2117_fpk, &
          42.2240_fpk,   42.2001_fpk,   42.2001_fpk,   16.2529_fpk,   18.3372_fpk, &
          18.3372_fpk,   16.3876_fpk,   16.2529_fpk,   16.2529_fpk,  546.7050_fpk, &
         440.5620_fpk,  345.8500_fpk,    2.0843_fpk,   18.3372_fpk,  262.5829_fpk, &
         190.7749_fpk,  130.4376_fpk,   81.5805_fpk,   44.2117_fpk,   44.2117_fpk, &
          81.5805_fpk,   18.3372_fpk,  130.4376_fpk,  190.7749_fpk,  262.5829_fpk, &
         345.8500_fpk,    3.9611_fpk,  440.5620_fpk,  546.7050_fpk,   16.3876_fpk, &
          16.2529_fpk,  544.8628_fpk,  438.7015_fpk,    2.0843_fpk,  343.9697_fpk, &
         260.6826_fpk,  188.8532_fpk,  128.4921_fpk,   16.3876_fpk,   79.6070_fpk, &
          42.2240_fpk,   42.2001_fpk,   79.5646_fpk,  128.3978_fpk,  188.7135_fpk, &
         260.5011_fpk,    0.0000_fpk,   16.2529_fpk,  343.7484_fpk,  438.4418_fpk, &
         544.5658_fpk,    3.9611_fpk,    2.0843_fpk,    2.0843_fpk,    3.9611_fpk, &
           0.0000_fpk,    2.0843_fpk,   18.3372_fpk,   16.3876_fpk,   16.3876_fpk, &
          18.3372_fpk,   16.2529_fpk,   16.3876_fpk,   44.2117_fpk,   42.2240_fpk, &
          44.2117_fpk,   42.2001_fpk,   42.2240_fpk,    2.0843_fpk,   81.5805_fpk, &
          79.5646_fpk,   81.5805_fpk,   79.6070_fpk,   79.5646_fpk,  130.4376_fpk, &
         128.3978_fpk,  130.4376_fpk,  128.4921_fpk,  128.3978_fpk,  190.7749_fpk, &
         188.7135_fpk,  190.7749_fpk,  188.8532_fpk,  188.7135_fpk,  262.5829_fpk, &
         260.5011_fpk,  262.5829_fpk,  260.6826_fpk,  260.5011_fpk,  345.8500_fpk, &
         343.7484_fpk,  345.8500_fpk,  343.9697_fpk,  343.7484_fpk,  440.5620_fpk, &
         438.4418_fpk,  440.5620_fpk,  438.7015_fpk,  438.4418_fpk,  546.7050_fpk, &
         544.5658_fpk,  546.7050_fpk,  544.8628_fpk,  544.5658_fpk,  662.1030_fpk, &
         664.2610_fpk,  662.4368_fpk,  791.0344_fpk,  793.2100_fpk,  791.4045_fpk, &
         931.3390_fpk,  933.5330_fpk,  931.7450_fpk, 1082.9941_fpk, 1085.2060_fpk, &
        1083.4356_fpk, 1245.9750_fpk, 1248.2040_fpk, 1246.4518_fpk, 1420.2552_fpk, &
        1422.5020_fpk, 1420.7672_fpk, 1605.8064_fpk, 1608.0710_fpk, 1606.3533_fpk /

!  O2 Placzek-teller coefficients

!      data o2plactel_kc99 /
!     * 3.630d-1, 3.630d-1, 3.637d-1, 3.622d-1, 3.622d-1,
!     * 3.630d-1, 3.613d-1, 3.613d-1, 3.622d-1, 3.602d-1,
!     * 3.602d-1, 3.612d-1, 3.589d-1, 3.589d-1, 3.601d-1,
!     * 3.574d-1, 3.574d-1, 3.589d-1, 3.556d-1, 3.556d-1,
!     * 3.573d-1, 2.079d-3, 3.533d-1, 3.534d-1, 3.555d-1,
!     * 2.191d-3, 2.597d-3, 3.505d-1, 3.506d-1, 3.532d-1,
!     * 2.755d-3, 3.337d-3, 3.468d-1, 3.471d-1, 3.504d-1,
!     * 3.567d-3, 4.444d-3, 3.418d-1, 3.422d-1, 3.467d-1,
!     * 4.800d-3, 6.211d-3, 3.348d-1, 3.354d-1, 3.416d-1,
!     * 6.803d-3, 9.288d-3, 3.220d-1, 3.251d-1, 3.344d-1,
!     * 1.127d-2, 1.062d-3, 1.387d-2, 3.013d-1, 3.077d-1,
!     * 3.223d-1, 1.979d-2, 2.613d-2, 2.544d-1, 2.727d-1,
!     * 3.020d-1, 1.476d-3, 4.342d-2, 9.234d-2, 6.596d-2,
!     * 1.714d-1, 2.571d-1, 1.843d-2, 1.615d-1, 1.970d-3,
!     * 2.445d-3, 3.116d-3, 1.077d-1, 7.690d-2, 4.105d-3,
!     * 5.652d-3, 8.271d-3, 1.222d-2, 2.842d-2, 2.207d-2,
!     * 1.470d-2, 5.132d-2, 8.256d-3, 5.647d-3, 4.103d-3,
!     * 3.115d-3, 2.308d-1, 2.445d-3, 1.970d-3, 2.116d-3,
!     * 3.810d-3, 2.076d-3, 2.593d-3, 1.385d-1, 3.329d-3,
!     * 4.431d-3, 6.184d-3, 9.227d-3, 3.991d-2, 1.696d-2,
!     * 1.867d-2, 3.474d-2, 1.078d-2, 7.483d-3, 5.200d-3,
!     * 3.822d-3, 5.383d-1, 1.077d-1, 2.927d-3, 2.313d-3,
!     * 1.874d-3, 2.692d-1, 1.843d-2, 4.628d-1, 4.000d-1,
!     * 4.617d-1, 9.234d-2, 5.583d-2, 1.476d-3, 4.362d-1,
!     * 4.286d-1, 4.579d-1, 3.193d-2, 2.339d-2, 4.214d-1,
!     * 4.196d-1, 4.352d-1, 1.600d-2, 1.911d-3, 1.278d-2,
!     * 4.130d-1, 4.118d-1, 4.210d-1, 1.038d-2, 7.519d-3,
!     * 4.067d-1, 4.060d-1, 4.135d-1, 6.803d-3, 5.217d-3,
!     * 4.021d-1, 4.017d-1, 4.070d-1, 4.800d-3, 3.831d-3,
!     * 3.987d-1, 3.985d-1, 4.023d-1, 3.567d-3, 2.933d-3,
!     * 3.961d-1, 3.959d-1, 3.988d-1, 2.755d-3, 2.317d-3,
!     * 3.939d-1, 3.938d-1, 3.961d-1, 2.191d-3, 1.876d-3,
!     * 3.922d-1, 3.921d-1, 3.940d-1, 1.785d-3, 3.908d-1,
!     * 3.907d-1, 3.922d-1, 3.895d-1, 3.895d-1, 3.908d-1,
!     * 3.885d-1, 3.885d-1, 3.896d-1, 3.876d-1, 3.876d-1,
!     * 3.885d-1, 3.868d-1, 3.868d-1, 3.876d-1, 3.861d-1,
!     * 3.861d-1, 3.868d-1, 3.855d-1, 3.855d-1, 3.861d-1 /

      DATA O2PLACTEL / &
       3.630e-1_fpk, 3.630e-1_fpk, 3.637e-1_fpk, 3.622e-1_fpk, 3.622e-1_fpk, &
       3.630e-1_fpk, 3.612e-1_fpk, 3.613e-1_fpk, 3.622e-1_fpk, 3.602e-1_fpk, &
       3.602e-1_fpk, 3.612e-1_fpk, 3.589e-1_fpk, 3.589e-1_fpk, 3.602e-1_fpk, &
       3.574e-1_fpk, 3.574e-1_fpk, 3.589e-1_fpk, 3.555e-1_fpk, 3.556e-1_fpk, &
       3.574e-1_fpk, 2.003e-3_fpk, 3.533e-1_fpk, 3.534e-1_fpk, 3.555e-1_fpk, &
       2.275e-3_fpk, 2.492e-3_fpk, 3.504e-1_fpk, 3.506e-1_fpk, 3.533e-1_fpk, &
       2.874e-3_fpk, 3.184e-3_fpk, 3.467e-1_fpk, 3.471e-1_fpk, 3.505e-1_fpk, &
       3.743e-3_fpk, 4.209e-3_fpk, 3.417e-1_fpk, 3.422e-1_fpk, 3.468e-1_fpk, &
       5.076e-3_fpk, 5.822e-3_fpk, 3.345e-1_fpk, 3.354e-1_fpk, 3.418e-1_fpk, &
       7.271e-3_fpk, 8.577e-3_fpk, 3.233e-1_fpk, 3.251e-1_fpk, 3.347e-1_fpk, &
       1.127e-2_fpk, 1.030e-3_fpk, 1.387e-2_fpk, 3.037e-1_fpk, 3.077e-1_fpk, &
       3.236e-1_fpk, 1.979e-2_fpk, 2.613e-2_fpk, 2.599e-1_fpk, 2.727e-1_fpk, &
       3.045e-1_fpk, 1.458e-3_fpk, 4.342e-2_fpk, 9.234e-2_fpk, 6.596e-2_fpk, &
       1.714e-1_fpk, 2.627e-1_fpk, 1.775e-2_fpk, 1.615e-1_fpk, 1.905e-3_fpk, &
       2.356e-3_fpk, 2.989e-3_fpk, 1.077e-1_fpk, 7.690e-2_fpk, 3.916e-3_fpk, &
       5.353e-3_fpk, 7.753e-3_fpk, 1.222e-2_fpk, 2.842e-2_fpk, 2.207e-2_fpk, &
       1.470e-2_fpk, 5.132e-2_fpk, 8.967e-3_fpk, 6.036e-3_fpk, 4.338e-3_fpk, &
       3.268e-3_fpk, 2.308e-1_fpk, 2.550e-3_fpk, 2.045e-3_fpk, 2.073e-4_fpk, &
       3.731e-4_fpk, 2.156e-3_fpk, 2.705e-3_fpk, 1.385e-1_fpk, 3.493e-3_fpk, &
       4.685e-3_fpk, 6.610e-3_fpk, 1.002e-2_fpk, 3.991e-2_fpk, 1.696e-2_fpk, &
       1.867e-2_fpk, 3.474e-2_fpk, 1.078e-2_fpk, 7.014e-3_fpk, 4.924e-3_fpk, &
       3.646e-3_fpk, 5.383e-1_fpk, 1.077e-1_fpk, 2.808e-3_fpk, 2.229e-3_fpk, &
       1.812e-3_fpk, 2.692e-1_fpk, 1.775e-2_fpk, 4.729e-1_fpk, 4.000e-1_fpk, &
       4.617e-1_fpk, 9.234e-2_fpk, 5.583e-2_fpk, 1.458e-3_fpk, 4.398e-1_fpk, &
       4.286e-1_fpk, 4.678e-1_fpk, 3.193e-2_fpk, 2.339e-2_fpk, 4.232e-1_fpk, &
       4.196e-1_fpk, 4.387e-1_fpk, 1.600e-2_fpk, 1.854e-3_fpk, 1.278e-2_fpk, &
       4.134e-1_fpk, 4.118e-1_fpk, 4.228e-1_fpk, 9.586e-3_fpk, 8.037e-3_fpk, &
       4.069e-1_fpk, 4.060e-1_fpk, 4.132e-1_fpk, 6.377e-3_fpk, 5.517e-3_fpk, &
       4.022e-1_fpk, 4.017e-1_fpk, 4.068e-1_fpk, 4.546e-3_fpk, 4.020e-3_fpk, &
       3.988e-1_fpk, 3.985e-1_fpk, 4.022e-1_fpk, 3.403e-3_fpk, 3.059e-3_fpk, &
       3.961e-1_fpk, 3.959e-1_fpk, 3.988e-1_fpk, 2.643e-3_fpk, 2.405e-3_fpk, &
       3.940e-1_fpk, 3.938e-1_fpk, 3.961e-1_fpk, 2.112e-3_fpk, 1.941e-3_fpk, &
       3.922e-1_fpk, 3.921e-1_fpk, 3.940e-1_fpk, 1.726e-3_fpk, 3.908e-1_fpk, &
       3.907e-1_fpk, 3.922e-1_fpk, 3.896e-1_fpk, 3.895e-1_fpk, 3.908e-1_fpk, &
       3.884e-1_fpk, 3.885e-1_fpk, 3.896e-1_fpk, 3.876e-1_fpk, 3.876e-1_fpk, &
       3.884e-1_fpk, 3.868e-1_fpk, 3.868e-1_fpk, 3.876e-1_fpk, 3.861e-1_fpk, &
       3.861e-1_fpk, 3.868e-1_fpk, 3.855e-1_fpk, 3.855e-1_fpk, 3.861e-1_fpk /

!  O2 degeneracy

      DATA O2DEG / &
       32._fpk, 33._fpk, 34._fpk, 30._fpk, 31._fpk, &
       32._fpk, 28._fpk, 29._fpk, 30._fpk, 26._fpk, &
       27._fpk, 28._fpk, 24._fpk, 25._fpk, 26._fpk, &
       22._fpk, 23._fpk, 24._fpk, 20._fpk, 21._fpk, &
       22._fpk, 19._fpk, 18._fpk, 19._fpk, 20._fpk, &
       18._fpk, 17._fpk, 16._fpk, 17._fpk, 18._fpk, &
       16._fpk, 15._fpk, 14._fpk, 15._fpk, 16._fpk, &
       14._fpk, 13._fpk, 12._fpk, 13._fpk, 14._fpk, &
       12._fpk, 11._fpk, 10._fpk, 11._fpk, 12._fpk, &
       10._fpk,  9._fpk,  8._fpk,  9._fpk, 10._fpk, &
        8._fpk,  4._fpk,  7._fpk,  6._fpk,  7._fpk, &
        8._fpk,  6._fpk,  5._fpk,  4._fpk,  5._fpk, &
        6._fpk,  4._fpk,  4._fpk,  2._fpk,  3._fpk, &
        3._fpk,  4._fpk,  2._fpk,  2._fpk, 19._fpk, &
       17._fpk, 15._fpk,  2._fpk,  3._fpk, 13._fpk, &
       11._fpk,  9._fpk,  7._fpk,  5._fpk,  5._fpk, &
        7._fpk,  3._fpk,  9._fpk, 11._fpk, 13._fpk, &
       15._fpk,  1._fpk, 17._fpk, 19._fpk,  4._fpk, &
        2._fpk, 18._fpk, 16._fpk,  2._fpk, 14._fpk, &
       12._fpk, 10._fpk,  8._fpk,  4._fpk,  6._fpk, &
        6._fpk,  4._fpk,  8._fpk, 10._fpk, 12._fpk, &
       14._fpk,  0._fpk,  2._fpk, 16._fpk, 18._fpk, &
       20._fpk,  1._fpk,  2._fpk,  2._fpk,  1._fpk, &
        0._fpk,  2._fpk,  3._fpk,  4._fpk,  4._fpk, &
        3._fpk,  2._fpk,  4._fpk,  5._fpk,  6._fpk, &
        5._fpk,  4._fpk,  6._fpk,  2._fpk,  7._fpk, &
        8._fpk,  7._fpk,  6._fpk,  8._fpk,  9._fpk, &
       10._fpk,  9._fpk,  8._fpk, 10._fpk, 11._fpk, &
       12._fpk, 11._fpk, 10._fpk, 12._fpk, 13._fpk, &
       14._fpk, 13._fpk, 12._fpk, 14._fpk, 15._fpk, &
       16._fpk, 15._fpk, 14._fpk, 16._fpk, 17._fpk, &
       18._fpk, 17._fpk, 16._fpk, 18._fpk, 19._fpk, &
       20._fpk, 19._fpk, 18._fpk, 20._fpk, 22._fpk, &
       21._fpk, 20._fpk, 24._fpk, 23._fpk, 22._fpk, &
       26._fpk, 25._fpk, 24._fpk, 28._fpk, 27._fpk, &
       26._fpk, 30._fpk, 29._fpk, 28._fpk, 32._fpk, &
       31._fpk, 30._fpk, 34._fpk, 33._fpk, 32._fpk /

!  some numbers

      C2   = 1.438769D0

!  the constants factor

      C0   = 256.0D0 * (PIE**5.0D0) / 27.0D0

!  set mixing ratio prefix

      PREFIX_N2 = C0 * MIXRATIO_N2
      PREFIX_O2 = C0 * MIXRATIO_O2

!  start layer loop

      DO N = 1, NLAYERS

        TEMP = TEMPERATURES(N)
        EMULT = - C2 / TEMP
        D_EMULT = - EMULT / TEMP

!  N2 partition function

        IF ( DO_TEMPERATURE_WF ) THEN
          QN2   = ZERO
          D_QN2 = ZERO
          DO I = 1, 31
            QN2_I =  (TWO * N2STAT_1 (I) + ONE) * N2STAT_2 (I) * &
                     EXP (EMULT * N2STAT_3 (I))
            D_QN2_I = QN2_I * D_EMULT * N2STAT_3 (I)
            QN2   = QN2   + QN2_I
            D_QN2 = D_QN2 + D_QN2_I
          END DO
          D_LOGQN2 = D_QN2 / QN2
        ELSE
          QN2 = ZERO
          DO I = 1, 31
            QN2 = QN2 + (TWO * N2STAT_1 (I) + ONE) * N2STAT_2 (I) * &
              EXP (EMULT * N2STAT_3 (I))
          END DO
        ENDIF

!  N2 Calculate population fractions for rotational Raman lines and
!     the cross sections, except for gamma**2/lambda**4

        IF ( DO_TEMPERATURE_WF ) THEN
          DO I = 1, 48
            PCPREF =  PREFIX_N2 * N2PLACTEL (I)
            D_TERM =  D_EMULT   * N2TERM    (I)
            N2FRAC = (TWO * N2DEG (I) + ONE) * N2NUC (I) * &
                        EXP (EMULT * N2TERM (I)) / QN2
            D_N2FRAC = N2FRAC * ( D_TERM - D_LOGQN2 )
            RRS_N2REF    (N,I) = PCPREF * N2FRAC
            LINRRS_N2REF (N,I) = PCPREF * D_N2FRAC
          END DO
        ELSE
          DO I = 1, 48
            N2FRAC = (TWO * N2DEG (I) + ONE) * N2NUC (I) * &
              EXP (EMULT * N2TERM (I)) / QN2
            RRS_N2REF (N,I) = PREFIX_N2 * N2FRAC * N2PLACTEL (I)
          END DO
        ENDIF

!  O2 Calculate partition function

        IF ( DO_TEMPERATURE_WF ) THEN
          QO2   = ZERO
          D_QO2 = ZERO
          DO I = 1, 54
            QO2_I = (TWO * o2stat_1 (I) + ONE) * &
                         EXP (EMULT * O2STAT_2 (I))
            D_QO2_I = QO2_I * D_EMULT * O2STAT_2 (I)
            QO2   = QO2   + QO2_I
            D_QO2 = D_QO2 + D_QO2_I
          END DO
          D_LOGQO2 = D_QO2 / QO2
        ELSE
          QO2 = ZERO
          DO I = 1, 54
            QO2 = QO2 + (TWO * o2stat_1 (I) + ONE) * &
                EXP (EMULT * O2STAT_2 (I))
          END DO
        ENDIF

!  O2 Calculate population fractions for rotational Raman lines and
!     the cross sections, except for gamma**2/lambda**4

        IF ( DO_TEMPERATURE_WF ) THEN
          DO I = 1, 185
            PCPREF =  PREFIX_O2 * O2PLACTEL (I)
            D_TERM =  D_EMULT   * O2TERM    (I)
            O2FRAC = (TWO * O2DEG (I) + ONE) * &
                EXP (EMULT * O2TERM (I)) / QO2
            D_O2FRAC = O2FRAC * ( D_TERM - D_LOGQO2 )
            RRS_O2REF    (N,I) = PCPREF * O2FRAC
            LINRRS_O2REF (N,I) = PCPREF * D_O2FRAC
          END DO
        ELSE
          DO I = 1, 185
            O2FRAC = (TWO * O2DEG (I) + ONE) * &
                EXP (EMULT * O2TERM (I)) / QO2
            RRS_O2REF (N,I) = PREFIX_O2 * O2FRAC * O2PLACTEL (I)
          END DO
        ENDIF

!  End layer loop

      ENDDO

!  Set the spectroscopy output

      GAMMA_O2(1) = GAMMA_O2_A
      GAMMA_O2(2) = GAMMA_O2_B
      GAMMA_O2(3) = GAMMA_O2_C
      GAMMA_N2(1) = GAMMA_N2_A
      GAMMA_N2(2) = GAMMA_N2_B
      GAMMA_N2(3) = GAMMA_N2_C
      DO I = 1, 48
        N2POS(I) = RN2POS(I)
      ENDDO
      DO I = 1, 185
        O2POS(I) = RO2POS(I)
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE RRS_BASIC_SPECTROSCOPY

!

      SUBROUTINE RRS_LOSSTERM_XSEC &
            ( MAXLAYERS, NLAYERS, N_RRS_N2REF, N_RRS_O2REF, & ! Inputs
              N2POS, O2POS, GAMMA_N2, GAMMA_O2,             & ! Inputs
              WAVELENGTH_INPUT, RRS_N2REF, RRS_O2REF,       & ! Inputs
              RRSXSEC_OUT, N_RRSOUT )                         ! Outputs

!  Loss term cross section
!  =======================

!  Given a wavelength WAVELENGTH_INPUT, and a bunch of temperatures
!  the compute the complete  cross-section of energy scattered out
!  of this wavelength. Requires as input the basic transition factors
!  from the previous module. Calculates internally the wavelength shifts
!  out from WAVELENGTH_INPUT, and includes lambda^-4 factor in result.

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

!  Output arguments
!  ----------------

!  Number of Raman shifts

      INTEGER  , INTENT(OUT) :: N_RRSOUT

!  Total Loss-term cross sections, one per layer

      REAL(FPK), INTENT(OUT) :: RRSXSEC_OUT ( MAXLAYERS )

!  Local variables
!  ---------------

      INTEGER   :: N, J
      REAL(FPK) :: SIGMA_0, SIGMA_4, SUMR, N2XSEC, O2XSEC
      REAL(FPK) :: SIGMA, SIGSQ
      REAL(FPK) :: GAMMAN2,    GAMMAO2
      REAL(FPK) :: GAMMAN2_SQ, GAMMAO2_SQ
      REAL(FPK) :: HELPN2 ( 48 )
      REAL(FPK) :: HELPO2 ( 185 )

      REAL(FPK), PARAMETER   :: DCONV24 = 1.0D-24

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
        DO J = 1, N_RRS_N2REF
          N2XSEC  = RRS_N2REF(N,J) * HELPN2(J)
          SUMR = SUMR + N2XSEC
        ENDDO
        DO J = 1, N_RRS_O2REF
          O2XSEC  = RRS_O2REF(N,J) * HELPO2(J)
          SUMR = SUMR + O2XSEC
        ENDDO
        RRSXSEC_OUT(N) = SUMR
      ENDDO

!   we need to multiply by d-48 to get the correct
!   units for the absolute cross sections in cm2: We calculate gamman2
!   and gammao2 in d-30 m3 and then calculate n2xsec and o2xsec values
!   as (gamma**2) * (sigprime**4), where sigprime is in cm-1.

      DO N = 1, NLAYERS
        RRSXSEC_OUT(N) = RRSXSEC_OUT(N) * DCONV24 * DCONV24
      ENDDO

!  number of shifts

      N_RRSOUT = N_RRS_N2REF + N_RRS_O2REF

!  Finish

      RETURN
      END SUBROUTINE RRS_LOSSTERM_XSEC

!

      SUBROUTINE RRS_GAINTERM_XSEC_BIN &
           ( MAXLAYERS, MAXSHIFTS, NLAYERS, N_RRS_N2REF, N_RRS_O2REF, & ! Inputs
             N2POS, O2POS, GAMMA_N2, GAMMA_O2,                        & ! Inputs
             WAVELENGTH_INPUT, RRS_N2REF, RRS_O2REF,                  & ! Inputs
             N_STOKES,     RRSXSEC_STOKES,                            & ! Inputs
             N_ANTISTOKES, RRSXSEC_ANTISTOKES,                        & ! Inputs
             WAVELENGTHS_SHIFTED, RRSXSEC_SHIFTED )                     ! Outputs       

!  Gain term cross-sections
!  ========================

!  Given a wavelength WAVELENGTH_INPUT,
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

!  Output arguments
!  ----------------

!  Number of Stokes/Anti-Stokes lines

      INTEGER  , intent(OUT) :: N_STOKES, N_ANTISTOKES

!  Raman shifted wavelengths

      REAL(FPK), INTENT(OUT) :: WAVELENGTHS_SHIFTED ( MAXSHIFTS )

!  Raman cross-sections (all)

      REAL(FPK), INTENT(OUT) :: RRSXSEC_SHIFTED     ( MAXLAYERS, MAXSHIFTS )

!  Stokes and Antistokes cross sections (debug only)

      REAL(FPK), INTENT(OUT) :: RRSXSEC_STOKES      ( MAXLAYERS, MAXSHIFTS )
      REAL(FPK), INTENT(OUT) :: RRSXSEC_ANTISTOKES  ( MAXLAYERS, MAXSHIFTS )

!  Local variables
!  ---------------

      INTEGER   :: N, J, S, NSHIFTS, S1, S2
      REAL(FPK) :: ANTISTOKES_SHIFTED ( 233 )
      REAL(FPK) :: STOKES_SHIFTED     ( 233 )
      REAL(FPK) :: SIGMA, SIGSQ, SIGMA_0, SIGMA_2, SIGMA_4
      REAL(FPK) :: GAMMAN2, GAMMAO2, GAMMAN2_SQ, GAMMAO2_SQ
      REAL(FPK) :: HELPN2, HELPO2

      REAL(FPK), PARAMETER :: DCONV24 = 1.0D-24

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
          RRSXSEC_SHIFTED(N,S) = RRS_N2REF(N,J) * HELPN2
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
          RRSXSEC_SHIFTED(N,S) = RRS_O2REF(N,J) * HELPO2
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
          RRSXSEC_SHIFTED(N,S) = RRSXSEC_SHIFTED(N,S)*DCONV24*DCONV24
        ENDDO
      ENDDO
      DO S = 1, N_STOKES
        DO N = 1, NLAYERS
          RRSXSEC_STOKES(N,S) = RRSXSEC_STOKES(N,S)*DCONV24*DCONV24
        ENDDO
      ENDDO
      DO S = 1, N_ANTISTOKES
        DO N = 1, NLAYERS
          RRSXSEC_ANTISTOKES(N,S) = &
                  RRSXSEC_ANTISTOKES(N,S)*DCONV24*DCONV24
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE RRS_GAINTERM_XSEC_BIN

!

      SUBROUTINE RRS_GAINTERM_XSEC_MONO &
           ( MAXLAYERS, MAXSHIFTS, NLAYERS, N_RRS_N2REF, N_RRS_O2REF, & ! Inputs
             N2POS, O2POS, GAMMA_N2, GAMMA_O2,                        & ! Inputs
             WAVELENGTH_INPUT, RRS_N2REF, RRS_O2REF,                  & ! Inputs
             WAVELENGTHS_SHIFTED, RRSXSEC_SHIFTED )                     ! Outputs

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

!  Output arguments
!  ----------------

!  Raman shifted wavelengths

      REAL(FPK), INTENT(OUT) :: WAVELENGTHS_SHIFTED ( MAXSHIFTS )

!  Raman cross-sections (all)

      REAL(FPK), INTENT(OUT) :: RRSXSEC_SHIFTED     ( MAXLAYERS, MAXSHIFTS )

!  Local variables
!  ---------------

      REAL(FPK) :: SIGMA, SIGSQ, SIGMA_0, SIGMA_2, SIGMA_4
      REAL(FPK) :: GAMMAN2, GAMMAO2, GAMMAN2_SQ, GAMMAO2_SQ
      REAL(FPK) :: HELPN2, HELPO2, LAMDA

      REAL(FPK), PARAMETER :: DCONV24 = 1.0D-24

      INTEGER   :: N, J, S, NSHIFTS

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
          RRSXSEC_SHIFTED(N,S) = RRS_N2REF(N,J) * HELPN2
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
          RRSXSEC_SHIFTED(N,S) = RRS_O2REF(N,J) * HELPO2
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
          RRSXSEC_SHIFTED(N,S) = RRSXSEC_SHIFTED(N,S)*DCONV24*DCONV24
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE RRS_GAINTERM_XSEC_MONO

!

      SUBROUTINE RRS_LAMDASRANKED_BIN &
           ( LAMBDA_EXCIT,     & ! Input
             LAMBDAS_SHIFTED )   ! Modified Output

!  Given a wavelength LAMBDA_EXCIT, Calculates Raman-shifted wavelengths
!  Output is ranked in ascending order

      USE LRRS_AUX2_m, Only : RQSORT_IDX

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Excitation wavelength

      REAL(FPK), INTENT(IN) :: LAMBDA_EXCIT

!  Output arguments
!  ----------------

!  Ranked wavelengths

      REAL(FPK), INTENT(OUT) :: LAMBDAS_SHIFTED ( 233 )

!  Local variables
!  ---------------

      INTEGER, PARAMETER :: N_RRS_N2REF = 48, &
                            N_RRS_O2REF = 185

!  Positions

      REAL(FPK) :: O2POS1 (N_RRS_O2REF)
      REAL(FPK) :: N2POS1 (N_RRS_N2REF)

!  help variables

      REAL(FPK) :: SIGMA, SIGMA_0, LAMBDAS_ALL ( 233 )
      INTEGER   :: I, W, J, S, INDEX_ALL ( 233 )
      INTEGER   :: NSHIFT

!  N2 Spectroscopy Transitions (wavenumber shifts)

      N2POS1 = (/ &
       -194.3015, -186.4226, -178.5374, -170.6459, -162.7484, &
       -154.8453, -146.9368, -139.0233, -131.1049, -123.1819, &
       -115.2547, -107.3235, -99.3886, -91.4502, -83.5086, &
       -75.5642, -67.6171, -59.6676, -51.7162, -43.7629, &
       -35.8080, -27.8519, -19.8950, -11.9373, 11.9373, &
       19.8950, 27.8519, 35.8080, 43.7629, 51.7162, &
       59.6676, 67.6171, 75.5642, 83.5086, 91.4502, &
       99.3886, 107.3235, 115.2547, 123.1819, 131.1049, &
       139.0233, 146.9368, 154.8453, 162.7484, 170.6459, &
       178.5374, 186.4226, 194.3015 /)

!  O2 spectroscopy Transitions

      O2POS1 = (/ &
       -185.5861, -185.5690, -185.5512, -174.3154, -174.2980, &
       -174.2802, -163.0162, -162.9980, -162.9809, -151.6906, &
       -151.6730, -151.6551, -140.3405, -140.3230, -140.3046, &
       -128.9677, -128.9490, -128.9314, -117.5740, -117.5560, &
       -117.5372, -108.2632, -106.1613, -106.1430, -106.1240, &
       -104.3008, -96.8136, -94.7318, -94.7120, -94.6934, &
       -92.8515, -85.3489, -83.2871, -83.2671, -83.2473, &
       -81.3868, -73.8694, -71.8294, -71.8080, -71.7876, &
       -69.9077, -62.3771, -60.3611, -60.3373, -60.3157, &
       -58.4156, -50.8730, -48.8851, -48.8571, -48.8332, &
       -46.9116, -40.1158, -39.3565, -37.4069, -37.3688, &
       -37.3406, -35.3953, -27.8241, -25.9472, -25.8745, &
       -25.8364, -25.8125, -23.8629, -16.2529, -16.2529, &
       -14.3761, -14.3033, -14.1686, -12.2918, -2.1392, &
       -2.1202, -2.1016, -2.0843, -2.0843, -2.0818, &
       -2.0614, -2.0398, -2.0159, -2.0116, -1.9877, &
       -1.9735, -1.9496, -1.9455, -1.9217, -1.9003, &
       -1.8803, -1.8768, -1.8605, -1.8422, -0.1347, &
       0.1347, 1.8422, 1.8605, 1.8768, 1.8803, &
       1.9003, 1.9217, 1.9455, 1.9496, 1.9735, &
       1.9877, 2.0116, 2.0159, 2.0398, 2.0614, &
       2.0818, 2.0843, 2.0843, 2.1016, 2.1202, &
       2.1392, 12.2918, 14.1686, 14.3033, 14.3761, &
       16.2529, 16.2529, 23.8629, 25.8125, 25.8364, &
       25.8745, 25.9472, 27.8241, 35.3953, 37.3406, &
       37.3688, 37.4069, 39.3565, 40.1158, 46.9116, &
       48.8332, 48.8571, 48.8851, 50.8730, 58.4156, &
       60.3157, 60.3373, 60.3611, 62.3771, 69.9077, &
       71.7876, 71.8080, 71.8294, 73.8694, 81.3868, &
       83.2473, 83.2671, 83.2871, 85.3489, 92.8515, &
       94.6934, 94.7120, 94.7318, 96.8136, 104.3008, &
       106.1240, 106.1430, 106.1613, 108.2632, 115.7318, &
       117.5372, 117.5560, 117.5740, 119.6952, 128.9314, &
       128.9490, 128.9677, 140.3046, 140.3230, 140.3405, &
       151.6551, 151.6730, 151.6906, 162.9809, 162.9980, &
       163.0162, 174.2802, 174.2980, 174.3154, 185.5512, &
       185.5690, 185.5861, 196.7919, 196.8100, 196.8269 /)

!  Initialize

      LAMBDAS_SHIFTED = zero

!  baseline values of GAMMA

      SIGMA_0 = 1.0D+07 / LAMBDA_EXCIT

!  initialize counts

      S  = 0

!  Shifted wavelengths (N2)

      DO J = 1, N_RRS_N2REF
        S = S + 1
        SIGMA = SIGMA_0 - N2POS1(J)
        LAMBDAS_ALL(S) = 1.0D+7 / SIGMA
      ENDDO

!  Shifted wavelengths (O2)

      DO J = 1, N_RRS_O2REF
        S = S + 1
        SIGMA = SIGMA_0 - O2POS1(J)
        LAMBDAS_ALL(S) = 1.0D+7 / SIGMA
      ENDDO

!  Total number of shifts

      NSHIFT = N_RRS_N2REF + N_RRS_O2REF

!  Use RQSORT_IDX to rank the wavelengths

      DO I=1,NSHIFT
        INDEX_ALL(I) = I
      ENDDO
      CALL RQSORT_IDX(NSHIFT, LAMBDAS_ALL, NSHIFT, INDEX_ALL)
      DO W = 1, NSHIFT
        LAMBDAS_SHIFTED(W) = LAMBDAS_ALL(W)
        LAMBDAS_SHIFTED(W) = LAMBDAS_ALL(INDEX_ALL(W))
      ENDDO

!  Finish routine

      RETURN
      END SUBROUTINE RRS_LAMDASRANKED_BIN

!

      SUBROUTINE RRS_LAMDASRANKED_MONO &
           ( MAX_POINTS, LAMBDA_EXCIT,               & ! Inputs
             LAMBDAS_RANKED, NPOINTS_MONO, W_EXCIT )   ! Outputs

!  Given a wavelength WAVELENGTH_INPUT, Calculates Raman-shifted wavelen
!  and forms the Ranked wavelengths for the Monochromatic case.

      USE LRRS_AUX2_m, Only : RQSORT_IDX

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Dimensioning parameters

      INTEGER  , INTENT(IN) :: MAX_POINTS

!  Excitation wavelength

      REAL(FPK), INTENT(IN) :: LAMBDA_EXCIT

!  Output arguments
!  ----------------

!  Ranked wavelengths

      REAL(FPK), INTENT(OUT) :: LAMBDAS_RANKED ( MAX_POINTS )

!  Index for excitation wavelength and number of MONO points (234)

      INTEGER  , INTENT(OUT) :: NPOINTS_MONO, W_EXCIT

!  Local variables
!  ---------------

      INTEGER  , PARAMETER :: N_RRS_N2REF = 48
      INTEGER  , PARAMETER :: N_RRS_O2REF = 185

!  Positions

      REAL(FPK) :: O2POS (N_RRS_O2REF)
      REAL(FPK) :: N2POS (N_RRS_N2REF)

!  help variables

      REAL(FPK) :: SIGMA, SIGMA_0, LAMBDAS_ALL ( MAX_POINTS )
      INTEGER   :: I, W, J, S, NSHIFTS, INDEX_ALL ( MAX_POINTS )

!  N2 Spectroscopy Transitions (wavenumber shifts)

      N2POS = (/ &
       -194.3015, -186.4226, -178.5374, -170.6459, -162.7484, &
       -154.8453, -146.9368, -139.0233, -131.1049, -123.1819, &
       -115.2547, -107.3235, -99.3886, -91.4502, -83.5086, &
       -75.5642, -67.6171, -59.6676, -51.7162, -43.7629, &
       -35.8080, -27.8519, -19.8950, -11.9373, 11.9373, &
       19.8950, 27.8519, 35.8080, 43.7629, 51.7162, &
       59.6676, 67.6171, 75.5642, 83.5086, 91.4502, &
       99.3886, 107.3235, 115.2547, 123.1819, 131.1049, &
       139.0233, 146.9368, 154.8453, 162.7484, 170.6459, &
       178.5374, 186.4226, 194.3015 /)

!  O2 spectroscopy Transitions

      O2POS = (/ &
       -185.5861, -185.5690, -185.5512, -174.3154, -174.2980, &
       -174.2802, -163.0162, -162.9980, -162.9809, -151.6906, &
       -151.6730, -151.6551, -140.3405, -140.3230, -140.3046, &
       -128.9677, -128.9490, -128.9314, -117.5740, -117.5560, &
       -117.5372, -108.2632, -106.1613, -106.1430, -106.1240, &
       -104.3008, -96.8136, -94.7318, -94.7120, -94.6934, &
       -92.8515, -85.3489, -83.2871, -83.2671, -83.2473, &
       -81.3868, -73.8694, -71.8294, -71.8080, -71.7876, &
       -69.9077, -62.3771, -60.3611, -60.3373, -60.3157, &
       -58.4156, -50.8730, -48.8851, -48.8571, -48.8332, &
       -46.9116, -40.1158, -39.3565, -37.4069, -37.3688, &
       -37.3406, -35.3953, -27.8241, -25.9472, -25.8745, &
       -25.8364, -25.8125, -23.8629, -16.2529, -16.2529, &
       -14.3761, -14.3033, -14.1686, -12.2918, -2.1392, &
       -2.1202, -2.1016, -2.0843, -2.0843, -2.0818, &
       -2.0614, -2.0398, -2.0159, -2.0116, -1.9877, &
       -1.9735, -1.9496, -1.9455, -1.9217, -1.9003, &
       -1.8803, -1.8768, -1.8605, -1.8422, -0.1347, &
       0.1347, 1.8422, 1.8605, 1.8768, 1.8803, &
       1.9003, 1.9217, 1.9455, 1.9496, 1.9735, &
       1.9877, 2.0116, 2.0159, 2.0398, 2.0614, &
       2.0818, 2.0843, 2.0843, 2.1016, 2.1202, &
       2.1392, 12.2918, 14.1686, 14.3033, 14.3761, &
       16.2529, 16.2529, 23.8629, 25.8125, 25.8364, &
       25.8745, 25.9472, 27.8241, 35.3953, 37.3406, &
       37.3688, 37.4069, 39.3565, 40.1158, 46.9116, &
       48.8332, 48.8571, 48.8851, 50.8730, 58.4156, &
       60.3157, 60.3373, 60.3611, 62.3771, 69.9077, &
       71.7876, 71.8080, 71.8294, 73.8694, 81.3868, &
       83.2473, 83.2671, 83.2871, 85.3489, 92.8515, &
       94.6934, 94.7120, 94.7318, 96.8136, 104.3008, &
       106.1240, 106.1430, 106.1613, 108.2632, 115.7318, &
       117.5372, 117.5560, 117.5740, 119.6952, 128.9314, &
       128.9490, 128.9677, 140.3046, 140.3230, 140.3405, &
       151.6551, 151.6730, 151.6906, 162.9809, 162.9980, &
       163.0162, 174.2802, 174.2980, 174.3154, 185.5512, &
       185.5690, 185.5861, 196.7919, 196.8100, 196.8269 /)

!  Initialize

      LAMBDAS_RANKED = zero
      NPOINTS_MONO   = 0
      W_EXCIT        = 0

!  baseline values of GAMMA

      SIGMA_0 = 1.0D+07 / LAMBDA_EXCIT

!  initialize counts

      S  = 0

!  Shifted wavelengths (N2)

      DO J = 1, N_RRS_N2REF
        S = S + 1
        SIGMA = SIGMA_0 + N2POS(J)
        LAMBDAS_ALL(S) = 1.0D+7 / SIGMA
      ENDDO

!  Shifted wavelengths (O2)

      DO J = 1, N_RRS_O2REF
        S = S + 1
        SIGMA = SIGMA_0 + O2POS(J)
        LAMBDAS_ALL(S) = 1.0D+7 / SIGMA
      ENDDO

!  Total number of shifts

      NSHIFTS = N_RRS_N2REF + N_RRS_O2REF
      NPOINTS_MONO = NSHIFTS + 1
      LAMBDAS_ALL(NPOINTS_MONO) = LAMBDA_EXCIT

!  Use RQSORT_IDX to rank the wavelengths
!  Find the position of the excitation wavelength in the rank

      DO I=1,NPOINTS_MONO
        INDEX_ALL(I) = I
      ENDDO
      CALL RQSORT_IDX(NPOINTS_MONO, LAMBDAS_ALL, NPOINTS_MONO, INDEX_ALL)

      DO W = 1, NPOINTS_MONO
        LAMBDAS_RANKED(W) = LAMBDAS_ALL(INDEX_ALL(W))
        IF ( LAMBDAS_RANKED(W) .EQ. LAMBDA_EXCIT) W_EXCIT = W
      ENDDO

!  Finish routine

      RETURN
      END SUBROUTINE RRS_LAMDASRANKED_MONO

!

      SUBROUTINE RRSBIN_REMAPPING &
         ( MAX_POINTS_0, MAX_POINTS,                            & ! Inputs
           NPOINTS_INNER_0, NPOINTS_OUTER_0, OFFSET_INNER_0,    & ! Inputs
           LAMBDAS_0, FLUXES_0, BINLOWER_0, BINUPPER_0,         & ! Inputs
           NPOINTS_INNER, NPOINTS_OUTER, OFFSET_INNER,          & ! Outputs
           LAMBDAS, FLUXES, BINLOWER, BINUPPER, FAIL, MESSAGE )   ! Outputs

!  Developed 28 March 2011, RT SOLUTIONS Inc.

!  Given a Binned solar spectrum for an outer range of RRS transitions,
!  This routine provides another binned solar spectrum as output, this
!  one with all zero bins removed.

!  In order to do this remapping, one has to pick out the empty bins in
!  in the original spectrum and omit them.

!  Original binned spectrum is input, with the "_0" extension

!  In this manner, we can keep the number of points and bins to a minimum.
!  The number of bins need be no greater than 233 for a single point.

       IMPLICIT NONE

!  Inputs
!  ======

!  Dimensioning

      INTEGER  , INTENT(IN) :: MAX_POINTS_0, MAX_POINTS

!  Input binned solar spectrum
!    - inner and outer buffer sizes + inner window offset
!    - Wavelengths and fluxes defined on the OUTER window
!    - Binning limits defined for the OUTER window

      INTEGER  , INTENT(IN) :: NPOINTS_INNER_0
      INTEGER  , INTENT(IN) :: NPOINTS_OUTER_0
      INTEGER  , INTENT(IN) :: OFFSET_INNER_0
      REAL(FPK), INTENT(IN) :: LAMBDAS_0  ( MAX_POINTS_0 )
      REAL(FPK), INTENT(IN) :: FLUXES_0   ( MAX_POINTS_0 )
      REAL(FPK), INTENT(IN) :: BINLOWER_0 ( MAX_POINTS_0 )
      REAL(FPK), INTENT(IN) :: BINUPPER_0 ( MAX_POINTS_0 )

!  Output Re-binned solar spectrum

      INTEGER  , INTENT(OUT) :: NPOINTS_INNER
      INTEGER  , INTENT(OUT) :: NPOINTS_OUTER
      INTEGER  , INTENT(OUT) :: OFFSET_INNER
      REAL(FPK), INTENT(OUT) :: LAMBDAS  ( MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: FLUXES   ( MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: BINLOWER ( MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: BINUPPER ( MAX_POINTS )

!  Errors

      LOGICAL          , INTENT(OUT) ::  FAIL
      CHARACTER (LEN=*), INTENT(OUT) :: MESSAGE

!  Local variables
!  ---------------

      INTEGER :: N, K, IK, NC, NCIN, J, B
      LOGICAL :: LOOP
      LOGICAL :: HITA(MAX_POINTS_0), INBUFFA(MAX_POINTS_0)
      REAL(FPK) :: LAMBDAS_SHIFTED ( 233 )
      INTEGER :: N_RRSBINS  ( MAX_POINTS )
      INTEGER :: BINMAP ( 233, MAX_POINTS )

!  Initialize
!  ----------

      NPOINTS_OUTER = 0
      NPOINTS_INNER = 0
      OFFSET_INNER  = 0
      LAMBDAS  = zero; FLUXES   = zero
      BINLOWER = zero; BINUPPER = zero

      FAIL = .FALSE. ; MESSAGE = ' '

!  Remapping
!  ---------

      N_RRSBINS = 0
      BINMAP    = 0

!  Start counting

      NC   = 0
      NCIN = 0

!      write(*,*)NPOINTS_INNER_0, OFFSET_INNER_0, NPOINTS_OUTER_0

!  Start outer window loop

      DO N = 1, NPOINTS_OUTER_0

!  Get the NSHIFT RRS shifted wavelengths from point I in this window

        CALL RRS_LAMDASRANKED_BIN &
           ( LAMBDAS_0(N),      & ! Input
             LAMBDAS_SHIFTED )    ! Output 

!  Start loop over points in inner window
!    Check if there any contributions or if this is an inner point
!    Compute remapped values

        DO K = 1, NPOINTS_INNER_0
          IK = K + OFFSET_INNER_0
          J = 0
          LOOP = .TRUE.
          DO WHILE ( LOOP .AND. J.LT.233 )
            J = J + 1
            IF ( LAMBDAS_SHIFTED(J).GT.BINLOWER_0(IK).AND. &
                 LAMBDAS_SHIFTED(J).LT.BINUPPER_0(IK) ) LOOP = .FALSE.
          ENDDO
          IF ( .NOT. LOOP ) THEN
            B = N_RRSBINS(K) + 1
            BINMAP(B,K) = N
            N_RRSBINS(K) = B
          ENDIF
        ENDDO

!  Finish outer loop

      ENDDO

!  Again

      HITA = .FALSE.; INBUFFA = .FALSE.
      DO N = 1, NPOINTS_OUTER_0
        DO K  = 1, NPOINTS_INNER_0
          IK = K + OFFSET_INNER_0
          DO B = 1, N_RRSBINS(K)
            IF ( BINMAP(B,K) .EQ. N .OR. IK .EQ. N ) THEN
               INBUFFA(N) = INBUFFA(N) .OR. ( IK.EQ.N )
!               IF ( IK.EQ.N )INBUFF(N) = .TRUE.
               HITA(N) = .TRUE.
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!  Remap

      NC = 0
      NCIN = 0
      DO N = 1, NPOINTS_OUTER_0
         IF ( HITA(N) ) THEN
           NC = NC + 1
           IF ( NC.GT. MAX_POINTS ) THEN
              FAIL = .TRUE.
              MESSAGE = 'MAX_POINTS DIMENSION INSUFFICIENT'
              RETURN
           ENDIF

!####  Rob Fix 3/29/11
           IF ( INBUFFA(N) ) THEN
               NCIN = NCIN + 1
               IF ( NCIN.EQ.1 ) OFFSET_INNER = NC - 1
           ENDIF
!####  End Rob Fix 3/29/11

           LAMBDAS(NC) = LAMBDAS_0(N)
           FLUXES(NC)  = FLUXES_0(N)
           BINLOWER(NC)= BINLOWER_0(N)
           BINUPPER(NC)= BINUPPER_0(N)
        ENDIF
      ENDDO
      NPOINTS_OUTER = NC
      NPOINTS_INNER = NCIN

!!####  Rob Fix 3/29/11, remove write statement
!      write(*,*)'lllllll',nc,ncin,offset_inner

!  Done

      RETURN
      END SUBROUTINE RRSBIN_REMAPPING

!  End Module

      END MODULE lrrs_raman_spectroscopy_m

