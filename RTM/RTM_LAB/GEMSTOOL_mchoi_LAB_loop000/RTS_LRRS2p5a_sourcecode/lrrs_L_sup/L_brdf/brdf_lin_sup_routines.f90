
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
! # Subroutines in this Module                                  #
! #                                                             #
! #              BRDF_LIN_MAKER, calling                        #
! #                BRDF_FUNCTION_PLUS                           #
! #              LIN_SCALING_FOURIER_ZERO (New, Version 3.7)    #
! #              BRDF_LIN_FOURIER                               #
! #                                                             #
! # New Cox-Munk Subroutine in this Module  (Version 3.7)       #
! #                                                             #
! #              BRDF_LIN_NewCM_MAKER                           #
! #                                                             #
! ###############################################################

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      MODULE brdf_LinSup_routines_m

      PRIVATE :: BRDF_FUNCTION_PLUS
      PUBLIC  :: BRDF_LIN_MAKER, BRDF_LIN_NewCM_MAKER, &
                 BRDF_LIN_FOURIER, LIN_SCALING_FOURIER_ZERO

      CONTAINS

      SUBROUTINE BRDF_LIN_MAKER                                           &
           ( DO_LOCAL_WSA, DO_LOCAL_BSA, DO_WSA_SCALING,                  & ! New line, Version 3.7
             WHICH_BRDF, DO_DBONLY,                                       & ! Inputs !@@
             DO_MSRCORR, DO_MSRCORR_DBONLY,                               & ! Inputs
             MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                          & ! Inputs
             BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                          & ! Inputs
             DO_USER_STREAMS,                                             & ! Inputs
             NSTREAMS_BRDF, NSTREAMS,                                     & ! Inputs
             N_USER_STREAMS, N_USER_RELAZMS,                              & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,          & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                & ! Inputs
             SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,  & ! New line, Version 3.7
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                & ! Inputs
             X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,                   & ! Inputs
             X_PHIQUAD, W_PHIQUAD,                                        & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                     & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0,                                   & ! Outputs
             SCALING_BRDFUNC, SCALING_BRDFUNC_0,                          & ! output, New line, Version 3.7
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,               & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0,                               & ! output
             D_SCALING_BRDFUNC, D_SCALING_BRDFUNC_0 )                       ! output, New line, Version 3.7

!  Prepares the bidirectional reflectance scatter matrices

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!     Observation-Geometry input control.         (DO_USER_OBSGEOMS)
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)

!  Rob Fix 9/25/14. Changed names with EXACTDB --->. DBKERNEL
!  Rob Fix 9/25/14. Two variables DO_EXACT, DO_EXACTONLY have been replaced
!                   with just the one variable DO_DBONLY

!  Rob Extension 12/2/14. BPDF Kernels (replace BREONVEG, BREONSOIL)

!  module, dimensions and numbers

      USE LRRS_pars_m, only : fpk, COXMUNK_IDX, MAX_BRDF_PARAMETERS, &
                              MAX_USER_RELAZMS, MAXSTREAMS_SCALING, &
                              MAX_STREAMS, MAX_USER_STREAMS, &
                              MAXSTREAMS_BRDF, MAXSTHALF_BRDF, &
                              MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD

      USE brdf_LinSup_kernels_m, only : COXMUNK_FUNCTION_MSR_PLUS, &
                                        BRDF_Generalized_Glint_plus

      implicit none

!  Input arguments
!  ===============

!  White-sky and Black-sky albedo scaling flags. New Version 3.7

      LOGICAL , intent(in)   :: DO_LOCAL_WSA           ! Required only for Regular                SCALING_BRDFUNCs
      LOGICAL , intent(in)   :: DO_LOCAL_BSA           ! Required both for Regular and linearized SCALING_BRDFUNCs
      LOGICAL , intent(in)   :: DO_WSA_SCALING         ! Required only for linearized             SCALING_BRDFUNCs

!  Which BRDF index

      INTEGER  , intent(in)  :: WHICH_BRDF

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!      LOGICAL  , intent(in)  :: DO_EXACT        ! Direct Bounce contribution is included
!      LOGICAL  , intent(in)  :: DO_EXACTONLY    ! Only Direct Bounce calculation (no Multiple scatter)
      LOGICAL  , intent(in)  :: DO_DBONLY        ! Only Direct Bounce calculation

!  Multiple reflectance correction for Glitter kernels

      LOGICAL  , intent(in)  :: DO_MSRCORR
      LOGICAL  , intent(in)  :: DO_MSRCORR_DBONLY ! Name change, Rob Fix 9/25/14
      INTEGER  , intent(in)  :: MSRCORR_ORDER
      INTEGER  , intent(in)  :: N_MUQUAD, N_PHIQUAD

!  Local number of parameters and local parameter array

      INTEGER  , intent(in)  :: BRDF_NPARS
      REAL(fpk), intent(in)  :: BRDF_PARS ( MAX_BRDF_PARAMETERS )
      LOGICAL  , intent(in) ::  BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

!  Local flags

      LOGICAL  , intent(in) :: DO_USER_STREAMS

!  Local angle control

      INTEGER  , intent(in) :: NSTREAMS
      INTEGER  , intent(in) :: N_USER_STREAMS
      INTEGER  , intent(in) :: N_USER_RELAZMS

!  Local angles

      REAL(fpk), intent(in) ::  PHIANG(MAX_USER_RELAZMS)
      REAL(fpk), intent(in) ::  COSPHI(MAX_USER_RELAZMS)
      REAL(fpk), intent(in) ::  SINPHI(MAX_USER_RELAZMS)

      REAL(fpk), intent(in) ::  SZASURCOS
      REAL(fpk), intent(in) ::  SZASURSIN

      REAL(fpk), intent(in) ::  QUAD_STREAMS(MAX_STREAMS)
      REAL(fpk), intent(in) ::  QUAD_SINES  (MAX_STREAMS)

      REAL(fpk), intent(in) ::  USER_STREAMS(MAX_USER_STREAMS)
      REAL(fpk), intent(in) ::  USER_SINES  (MAX_USER_STREAMS)

!  Discrete ordinates (local, for Albedo scaling). New Version 3.7

      INTEGER  , intent(in)  :: SCALING_NSTREAMS
      REAL(fpk), intent(in)  :: SCALING_QUAD_STREAMS(MAXSTREAMS_SCALING)
      REAL(fpk), intent(in)  :: SCALING_QUAD_SINES  (MAXSTREAMS_SCALING)

!  azimuth quadrature streams for BRDF

      INTEGER  , intent(in) ::  NSTREAMS_BRDF
      REAL(fpk), intent(in) ::  X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in) ::  CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in) ::  SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in) ::  CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(fpk), intent(in) ::  SXE_BRDF ( MAXSTHALF_BRDF )

!  Local arrays for MSR quadrature

      REAL(fpk), intent(in)  :: X_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk), intent(in)  :: W_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk), intent(in)  :: SX_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk), intent(in)  :: WXX_MUQUAD (MAX_MSRS_MUQUAD)

      REAL(fpk), intent(in)  :: X_PHIQUAD (MAX_MSRS_PHIQUAD)
      REAL(fpk), intent(in)  :: W_PHIQUAD (MAX_MSRS_PHIQUAD)

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(out) :: BRDFUNC   ( MAX_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: BRDFUNC_0 ( MAX_STREAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk), intent(out) :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: USER_BRDFUNC_0 ( MAX_USER_STREAMS, MAXSTREAMS_BRDF )

!  Direct Bounce (DB) values

      REAL(fpk), intent(out) :: DBKERNEL_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Output for WSA/BSA scaling options. New, Version 3.7

      REAL(fpk), intent(out) :: SCALING_BRDFUNC   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: SCALING_BRDFUNC_0 ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  Output Linearizations of BRDF functions (parameter derivatives)
!  ===============================================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(out) :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAX_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAX_STREAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk), intent(out) :: D_USER_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: D_USER_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS_BRDF )

!  Direct Bounce (DB) values

      REAL(fpk), intent(out) :: D_DBKERNEL_BRDFUNC ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Values for WSA/BSA scaling options. New, Version 3.7

      REAL(fpk), intent(out) ::  D_SCALING_BRDFUNC   &
           ( MAX_BRDF_PARAMETERS, MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) ::  D_SCALING_BRDFUNC_0 &
           ( MAX_BRDF_PARAMETERS, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      LOGICAL       :: MSRFLAG
      INTEGER       :: Q, I, UI, J, K, KE, IB
      REAL(fpk)     :: DFUNC ( MAX_BRDF_PARAMETERS )
      REAL(fpk)     :: KERNEL, D_KERNEL ( MAX_BRDF_PARAMETERS )

!  Rob Fix 9/25/14. DBFLAG has been renamed to MSRFLAG

     MSRFLAG = ( WHICH_BRDF .eq. COXMUNK_IDX ) .and. &
               ( DO_MSRCORR .or. DO_MSRCORR_DBONLY )

!  Direct-bounce calculation
!  -------------------------

!  CoxMunk special

      IF ( MSRFLAG ) THEN
         DO K = 1, N_USER_RELAZMS
            DO UI = 1, N_USER_STREAMS
               CALL COXMUNK_FUNCTION_MSR_PLUS &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, & ! Inputs
                MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                       & ! Inputs
                SZASURCOS, SZASURSIN, USER_STREAMS(UI),                   & ! Inputs
                USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),          & ! Inputs
                X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,                & ! Inputs
                X_PHIQUAD, W_PHIQUAD,                                     & ! Inputs
                DBKERNEL_BRDFUNC(UI,K), DFUNC )                             ! Output
              D_DBKERNEL_BRDFUNC(:,UI,K)  = DFUNC(:)
            ENDDO
         ENDDO

!  All other BRDFs

      ELSE
         DO K = 1, N_USER_RELAZMS
            DO UI = 1, N_USER_STREAMS
               CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                        & ! Inputs
                 BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                     & ! Inputs
                 SZASURCOS, SZASURSIN, USER_STREAMS(UI),                 & ! Inputs
                 USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),        & ! Inputs
                 DBKERNEL_BRDFUNC(UI,K), DFUNC )                         ! Output
                D_DBKERNEL_BRDFUNC(:,UI,K)  = DFUNC(:)
            ENDDO
         ENDDO
      ENDIF

!  SCALING OPTIONS (New Section, Version 3.7)
!  ------------------------------------------

!  White-sky albedo, scaling.
!     Use Local "Scaling_streams", both incident and outgoing

      IF ( DO_LOCAL_WSA .or. DO_WSA_SCALING ) THEN
         DO I = 1, SCALING_NSTREAMS
            DO J = 1, SCALING_NSTREAMS
               DO K = 1, NSTREAMS_BRDF
                  IF ( MSRFLAG ) THEN
                     CALL COXMUNK_FUNCTION_MSR_PLUS &
                      ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, & ! Inputs
                        MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                      & ! Inputs
                        SCALING_QUAD_STREAMS(J), SCALING_QUAD_SINES(J),          &
                        SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I),          &         
                        X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                       &
                        X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,               & ! Inputs
                        X_PHIQUAD, W_PHIQUAD,                                    & ! Inputs
                        KERNEL, D_KERNEL )                                         ! Output
                  ELSE
                     CALL BRDF_FUNCTION_PLUS &
                      ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                         & ! Inputs
                        BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                      & ! Inputs
                        SCALING_QUAD_STREAMS(J), SCALING_QUAD_SINES(J),          &
                        SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I),          &         
                        X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                       &
                        KERNEL, D_KERNEL )                                         ! Output
                  ENDIF
                  SCALING_BRDFUNC(I,J,K) = KERNEL
                  if ( DO_WSA_SCALING ) D_SCALING_BRDFUNC(:,I,J,K) = D_KERNEL(:)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  Black-sky albedo, scaling
!     Use Local "Scaling_streams" for outgoing, solar beam for incoming (IB = 1)

      IF ( DO_LOCAL_BSA ) THEN
         DO I = 1, SCALING_NSTREAMS
            DO K = 1, NSTREAMS_BRDF
               IF ( MSRFLAG ) THEN
                  CALL COXMUNK_FUNCTION_MSR_PLUS &
                      ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, & ! Inputs
                        MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                      & ! Inputs
                        SZASURCOS, SZASURSIN,                                    & ! Inputs
                        SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I),          & ! Inputs       
                        X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                       & ! Inputs
                        X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,               & ! Inputs
                        X_PHIQUAD, W_PHIQUAD,                                    & ! Inputs
                        KERNEL, D_KERNEL )                                         ! Output
               ELSE
                  CALL BRDF_FUNCTION_PLUS &
                      ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                         & ! Inputs
                        BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                      & ! Inputs
                        SZASURCOS, SZASURSIN,                                    & ! Inputs
                        SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I),          & ! Inputs         
                        X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                       & ! Inputs
                        KERNEL, D_KERNEL )                                         ! Output
               ENDIF
               SCALING_BRDFUNC_0(I,K) = KERNEL
               D_SCALING_BRDFUNC_0(:,I,K) = D_KERNEL(:)
            ENDDO
         ENDDO
      ENDIF

!  Return if the Direct-bounce BRDF is all that is required (scaled or not!)

      IF ( DO_DBONLY ) RETURN

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam
!    !@@  Solar Optionality. 12/31/12

      DO I = 1, NSTREAMS
        DO K = 1, NSTREAMS_BRDF
          CALL BRDF_FUNCTION_PLUS &
           ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                        & ! Inputs
             BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                     & ! Inputs
             SZASURCOS, SZASURSIN, QUAD_STREAMS(I),                  &
             QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
             BRDFUNC_0(I,K), DFUNC )
          DO Q = 1, BRDF_NPARS
            D_BRDFUNC_0(Q,I,K)  = DFUNC(Q)
          ENDDO
        ENDDO
      ENDDO

!  Incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                        & ! Inputs
                 BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                     & ! Inputs
                 QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),        &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                 BRDFUNC(I,J,K), DFUNC )
            DO Q = 1, BRDF_NPARS
              D_BRDFUNC(Q,I,J,K)  = DFUNC(Q)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam, Outgoing User-stream

         DO UI = 1, N_USER_STREAMS
            DO K = 1, NSTREAMS_BRDF
               CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                        & ! Inputs
                 BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                     & ! Inputs
                 SZASURCOS, SZASURSIN, USER_STREAMS(UI),                 &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                 USER_BRDFUNC_0(UI,K), DFUNC )
               DO Q = 1, BRDF_NPARS
                  D_USER_BRDFUNC_0(Q,UI,K)  = DFUNC(Q)
               ENDDO
            ENDDO
         ENDDO

!  Incident quadrature directions, Outgoing User-stream

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION_PLUS &
                 ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                        & ! Inputs
                   BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                     & ! Inputs
                   QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),       &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                   USER_BRDFUNC(UI,J,K), DFUNC )
              DO Q = 1, BRDF_NPARS
                D_USER_BRDFUNC(Q,UI,J,K)  = DFUNC(Q)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_LIN_MAKER

!

      SUBROUTINE BRDF_FUNCTION_PLUS  &
      ( MAXPARS, WHICH_BRDF, NPARS, PARS, DERIVS, &
        XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
        KERNEL, DKERNEL )

!  Rob Extension 12/2/14. BPDF Kernels (replace BREONVEG, BREONSOIL)

!  module, dimensions and numbers

      USE LRRS_pars_m, only : fpk, LISPARSE_IDX, LIDENSE_IDX,       &
                              RAHMAN_IDX, HAPKE_IDX, COXMUNK_IDX,   &
                              BPDFVEGN_IDX,   BPDFSOIL_IDX, BPDFNDVI_IDX
      USE brdf_LinSup_kernels_m

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: WHICH_BRDF
      INTEGER  , intent(in)  :: MAXPARS, NPARS
      LOGICAL  , intent(in)  :: DERIVS ( MAXPARS )
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: KERNEL
      REAL(fpk), intent(out) :: DKERNEL ( MAXPARS )

!  Trawl through

      IF ( WHICH_BRDF .EQ. LISPARSE_IDX ) THEN
        CALL LISPARSE_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LIDENSE_IDX ) THEN
        CALL LIDENSE_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. RAHMAN_IDX ) THEN
        CALL RAHMAN_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. HAPKE_IDX ) THEN
        CALL HAPKE_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. COXMUNK_IDX ) THEN
        CALL COXMUNK_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFVEGN_IDX ) THEN
        CALL BPDFVEGN_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFSOIL_IDX ) THEN
        CALL BPDFSOIL_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFNDVI_IDX ) THEN
        CALL BPDFNDVI_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_FUNCTION_PLUS

!

      SUBROUTINE LIN_SCALING_FOURIER_ZERO &
            ( DO_LOCAL_WSA, DO_LOCAL_BSA, LAMBERTIAN_FLAG,              &
              BRDF_NPARS, BRDF_DERIVS, SCALING_NSTREAMS, NSTREAMS_BRDF, &
              A_BRDF, D_SCALING_BRDFUNC, D_SCALING_BRDFUNC_0,           &
              D_SCALING_BRDF_F, D_SCALING_BRDF_F_0 )

!  include file of dimensions and numbers

      USE LRRS_pars_m, Only : fpk, zero, half, &
          MAX_BRDF_PARAMETERS, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF

      IMPLICIT NONE

!  This is a new routine for developing Fourier = 0 components for WSA/BSA computations.
!   Installed, 17 April 2014 for Version 3.7

!  Input arguments
!  ===============

!  Local flags

      LOGICAL  , intent(in) :: DO_LOCAL_WSA, DO_LOCAL_BSA

!  Control

      LOGICAL  , intent(in) :: LAMBERTIAN_FLAG

!  Local numbers

      INTEGER  , intent(in) :: SCALING_NSTREAMS, NSTREAMS_BRDF

!  Azimuth weights

      REAL(fpk), intent(in) :: A_BRDF ( MAXSTREAMS_BRDF )

!  linearization Control

      INTEGER  , intent(in) :: BRDF_NPARS
      LOGICAL  , intent(in) :: BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

!  Input for WSA/BSA scaling options. New, Version 3.7

      REAL(fpk), intent(in)  :: D_SCALING_BRDFUNC &
          ( MAX_BRDF_PARAMETERS, MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: D_SCALING_BRDFUNC_0 &
          ( MAX_BRDF_PARAMETERS, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  Output: Derivative-kernel Fourier components
!  ============================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(out)  :: D_SCALING_BRDF_F &
       ( MAX_BRDF_PARAMETERS, MAXSTREAMS_SCALING, MAXSTREAMS_SCALING )
      REAL(fpk), intent(out)  :: D_SCALING_BRDF_F_0 &
       ( MAX_BRDF_PARAMETERS, MAXSTREAMS_SCALING  )

!  local variables
!  ===============

      INTEGER   :: I, J, K, W
      real(fpk) :: SUM

!  Zeroing

      D_SCALING_BRDF_F        = ZERO
      D_SCALING_BRDF_F_0      = ZERO

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)
!    !@@ Solar Optionality, added 12/31/12

!  BSA: Incident Solar beam

      IF ( DO_LOCAL_BSA ) THEN
         IF ( .NOT. LAMBERTIAN_FLAG ) THEN
            DO W = 1, BRDF_NPARS
               IF ( BRDF_DERIVS(W) ) THEN
                  DO I = 1, SCALING_NSTREAMS
                     SUM = ZERO
                     DO K = 1, NSTREAMS_BRDF
                        SUM  = SUM + D_SCALING_BRDFUNC_0(W,I,K)*A_BRDF(K)
                     ENDDO
                     D_SCALING_BRDF_F_0(W,I) = SUM * HALF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDIF

!  WSA: incident quadrature directions

      if ( DO_LOCAL_WSA ) THEN
         IF ( .NOT. LAMBERTIAN_FLAG ) THEN
            DO W = 1, BRDF_NPARS
               IF ( BRDF_DERIVS(W) ) THEN
                  DO I = 1, SCALING_NSTREAMS
                     DO J = 1, SCALING_NSTREAMS
                        SUM = ZERO
                        DO K = 1, NSTREAMS_BRDF
                           SUM  = SUM + D_SCALING_BRDFUNC(W,I,J,K)*A_BRDF(K)
                        ENDDO
                        D_SCALING_BRDF_F(W,I,J) = SUM * HALF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LIN_SCALING_FOURIER_ZERO

!

      SUBROUTINE BRDF_Lin_NewCM_MAKER &
         ( DO_GlintShadow, DO_FacetIsotropy, WINDSPEED, WINDDIR, Refrac_R, Refrac_I, &
           WC_Reflectance, WC_Lambertian, DWC_Reflectance, DWC_Lambertian,           &
           DO_USER_STREAMS, DO_DBONLY,                                               &
           NSTREAMS_BRDF, NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS,                  &
           QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,                       &
           SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, X_BRDF, CX_BRDF, SX_BRDF,   &
           DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, BRDFUNC_0, USER_BRDFUNC_0,       &
           D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_USER_BRDFUNC_0 )     

!  include file of dimensions and numbers

      USE LRRS_pars_m,        only : MAX_USER_RELAZMS, MAX_USER_STREAMS, &
                                     MAX_STREAMS, MAXSTREAMS_BRDF,   &
                                     MAX_BRDF_PARAMETERS, fpk, zero, one, DEG_TO_RAD

      USE brdf_LinSup_kernels_m, only : BRDF_Generalized_Glint_plus

      IMPLICIT NONE

!  Prepares the bidirectional reflectance scatter matrices

!  Input arguments
!  ===============

!  NewCM Glitter options (bypasses the usual Kernel system)
!  --------------------------------------------------------

!  Flags for glint shadowing, Facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FacetIsotropy

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk):: WINDSPEED, WINDDIR

!  Refractive Index

      REAL(fpk) :: Refrac_R, Refrac_I

!  Whitecap correction (Zero if not flagged), linearization w.r.t Windspeed

      REAL(fpk) :: WC_Reflectance, WC_Lambertian
      REAL(fpk) :: DWC_Reflectance, DWC_Lambertian

!  Local flags

      LOGICAL ::   DO_USER_STREAMS

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL    :: DO_EXACT
!      LOGICAL   :: DO_EXACTONLY
      LOGICAL ::   DO_DBONLY

!  Number of Azimuth quadrature streams

      INTEGER ::   NSTREAMS_BRDF

!  Local angle control

      INTEGER ::   NSTREAMS
      INTEGER ::   N_USER_STREAMS
      INTEGER ::   N_USER_RELAZMS

!  Local angles

      REAL(fpk) :: PHIANG(MAX_USER_RELAZMS)
      REAL(fpk) :: COSPHI(MAX_USER_RELAZMS)
      REAL(fpk) :: SINPHI(MAX_USER_RELAZMS)

      REAL(fpk) :: SZASURCOS
      REAL(fpk) :: SZASURSIN

      REAL(fpk) :: QUAD_STREAMS(MAX_STREAMS)
      REAL(fpk) :: QUAD_SINES  (MAX_STREAMS)

      REAL(fpk) :: USER_STREAMS(MAX_USER_STREAMS)
      REAL(fpk) :: USER_SINES  (MAX_USER_STREAMS)

!  azimuth quadrature streams for BRDF

      REAL(fpk) :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk) :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk) :: SX_BRDF ( MAXSTREAMS_BRDF )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      REAL(fpk) :: BRDFUNC   ( MAX_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF )
      REAL(fpk) :: BRDFUNC_0 ( MAX_STREAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk) :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF )
      REAL(fpk) :: USER_BRDFUNC_0 ( MAX_USER_STREAMS, MAXSTREAMS_BRDF )

!  Exact DB values

      REAL(fpk) :: DBKERNEL_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Output Linearizations of BRDF functions (parameter derivatives)
!  ===============================================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk) :: D_BRDFUNC & 
                 ( MAX_BRDF_PARAMETERS, MAX_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF )
      REAL(fpk) :: D_BRDFUNC_0 &
                 ( MAX_BRDF_PARAMETERS, MAX_STREAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk) :: D_USER_BRDFUNC &
                 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF )
      REAL(fpk) :: D_USER_BRDFUNC_0 &
                 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS_BRDF )

!  Exact DB values

      REAL(fpk) :: D_DBKERNEL_BRDFUNC &
                 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  local variables
!  ---------------

      LOGICAL   :: DO_COEFFS, Local_Isotropy
      INTEGER   :: I, UI, J, K, IB
      REAL(fpk) :: PHI_W, CPHI_W, SPHI_W
      REAL(fpk) :: SUNGLINT_COEFFS(7), DSUNGLINT_COEFFS(7)
      REAL(fpk) :: WC_correction, DWC_correction, KERNEL, DKERNEL

!   Wind-direction and coefficient set-up

      DO_COEFFS = .true.
      PHI_W = zero ; CPHI_W = one ; SPHI_W = zero
      Local_Isotropy = DO_FacetIsotropy 
      if ( .not.Local_Isotropy ) then
            PHI_W  = WINDDIR
            CPHI_W = cos(WINDDIR * deg_to_rad) 
            SPHI_W = sin(WINDDIR * deg_to_rad)
      endif

!  Whitecap correction to glint

      WC_correction  = one - WC_Lambertian
      DWC_correction = - DWC_Lambertian

!  Direct Bounce calculation
!  -------------------------

      DO K = 1, N_USER_RELAZMS
         DO UI = 1, N_USER_STREAMS
            CALL BRDF_Generalized_Glint_plus &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,       &
               REFRAC_R, REFRAC_I, WINDSPEED,                   &
               PHI_W, CPHI_W, SPHI_W,                           &
               SZASURCOS, SZASURSIN, USER_STREAMS(UI),          &
               USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K), &
               SUNGLINT_COEFFS, DSUNGLINT_COEFFS, KERNEL, DKERNEL )
            DBKERNEL_BRDFUNC  (UI,K)   = WC_Reflectance  +  WC_correction * KERNEL
            D_DBKERNEL_BRDFUNC(1,UI,K) = DWC_Reflectance + DWC_correction * KERNEL &
                                                         +  WC_correction * DKERNEL
         ENDDO
      ENDDO

!      pause'after direct bounce'

!  Return if this is all you require

      IF ( DO_DBONLY ) RETURN

!  Incident Solar beam
!  ===================

!  Quadrature outgoing directions

      DO I = 1, NSTREAMS
         DO K = 1, NSTREAMS_BRDF
            CALL BRDF_Generalized_Glint_plus &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               REFRAC_R, REFRAC_I, WINDSPEED,                    &
               PHI_W, CPHI_W, SPHI_W,                            &
               SZASURCOS, SZASURSIN, QUAD_STREAMS(I),            &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, DSUNGLINT_COEFFS, KERNEL, DKERNEL )
            BRDFUNC_0  (I,K)   = WC_Reflectance  + WC_correction  * KERNEL
            D_BRDFUNC_0(1,I,K) = DWC_Reflectance + DWC_correction * KERNEL &
                                                 +  WC_correction * DKERNEL
         ENDDO
      ENDDO

!  User-streams outgoing directions
!   This is the "Truncated" Direct Bounce calculation

      IF ( DO_USER_STREAMS ) THEN
         DO UI = 1, N_USER_STREAMS
            DO K = 1, NSTREAMS_BRDF
               CALL BRDF_Generalized_Glint_plus &
                ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                  REFRAC_R, REFRAC_I, WINDSPEED,                     &
                  PHI_W, CPHI_W, SPHI_W,                             &
                  SZASURCOS, SZASURSIN, USER_STREAMS(UI),            &
                  USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                  SUNGLINT_COEFFS, DSUNGLINT_COEFFS, KERNEL, DKERNEL )
               USER_BRDFUNC_0  (UI,K)   = WC_Reflectance  + WC_correction  * KERNEL
               D_USER_BRDFUNC_0(1,UI,K) = DWC_Reflectance + DWC_correction * KERNEL &
                                                          +  WC_correction * DKERNEL
            ENDDO
         ENDDO
      ENDIF

!  Incident quadrature directions (MULTIPLE SCATTERING)
!  ==============================

!  Outgoing quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_Generalized_Glint_plus &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               REFRAC_R, REFRAC_I, WINDSPEED,                    &
               PHI_W, CPHI_W, SPHI_W,                            &
               QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),  &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, DSUNGLINT_COEFFS, KERNEL, DKERNEL )
            BRDFUNC  (I,J,K)   = WC_Reflectance  + WC_correction  * KERNEL
            D_BRDFUNC(1,I,J,K) = DWC_Reflectance + DWC_correction * KERNEL &
                                                 +  WC_correction * DKERNEL
          ENDDO
        ENDDO
      ENDDO

!  User stream outgoing directions

      IF ( DO_USER_STREAMS ) THEN
        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_Generalized_Glint_plus &
               ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                 REFRAC_R, REFRAC_I, WINDSPEED,                     &
                 PHI_W, CPHI_W, SPHI_W,                             &
                 QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),  &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 SUNGLINT_COEFFS, DSUNGLINT_COEFFS, KERNEL, DKERNEL )
              USER_BRDFUNC  (UI,J,K)   = WC_Reflectance  + WC_correction  * KERNEL
              D_USER_BRDFUNC(1,UI,J,K) = DWC_Reflectance + DWC_correction * KERNEL &
                                                         +  WC_correction * DKERNEL
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_Lin_NewCM_MAKER

!

      SUBROUTINE BRDF_LIN_FOURIER                                       &
         ( DO_USER_STREAMS, LOCAL_BRDF_NPARS, LOCAL_BRDF_DERIVS,        & ! Inputs
           LAMBERTIAN_FLAG, FACTOR, M, DELFAC,                          & ! Inputs
           NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,         & ! Inputs
           D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_USER_BRDFUNC_0,    & ! Inputs
           BRDF_AZMFAC, A_BRDF, BAX_BRDF,                               & ! Inputs
           D_LOCAL_BRDF_F,      D_LOCAL_BRDF_F_0,                       & ! Outputs
           D_LOCAL_USER_BRDF_F, D_LOCAL_USER_BRDF_F_0 )                   ! Outputs

!  Prepares Fourier component of the bidirectional reflectance functions

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!     Observation-Geometry input control.         (DO_USER_OBSGEOMS)
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)

!  module, dimensions and numbers

      USE LRRS_pars_m, only : fpk, ZERO, ONE, HALF, &
                              MAX_STREAMS, MAX_USER_STREAMS, &
                              MAXSTREAMS_BRDF, MAXSTHALF_BRDF, &
                              MAX_BRDF_PARAMETERS

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Control

      LOGICAL  , intent(in)  :: LAMBERTIAN_FLAG
      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      REAL(fpk), intent(in)  :: DELFAC, FACTOR
      INTEGER  , intent(in)  :: M, LOCAL_BRDF_NPARS
      LOGICAL  , intent(in)  :: LOCAL_BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

!  Local numbers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NSTREAMS_BRDF, NBRDF_HALF

!  Azimuth cosines and weights

      REAL(fpk), intent(in)  :: BRDF_AZMFAC ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: A_BRDF      ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: BAX_BRDF    ( MAXSTHALF_BRDF  )

!  Local Linearizations of BRDF functions (parameter derivatives)
!  ==============================================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(in)  :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAX_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF ) 
      REAL(fpk), intent(in)  :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAX_STREAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk), intent(in)  :: D_USER_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAX_STREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: D_USER_BRDFUNC_0 ( MAX_BRDF_PARAMETERS ,MAX_USER_STREAMS, MAXSTREAMS_BRDF )

!  Output: Derivative-kernel Fourier components
!  ============================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(out) :: D_LOCAL_BRDF_F   ( MAX_BRDF_PARAMETERS, MAX_STREAMS, MAX_STREAMS )
      REAL(fpk), intent(out) :: D_LOCAL_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAX_STREAMS )

!  at user-defined stream directions

      REAL(fpk), intent(out) :: D_LOCAL_USER_BRDF_F   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAX_STREAMS )
      REAL(fpk), intent(out) :: D_LOCAL_USER_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS)

!  local variables
!  ===============

      INTEGER      :: I, UI, J, K, KPHI, IB, Q
      REAL(fpk)    :: SUM, REFL, HELP

!  surface factor

      HELP = HALF * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)

      IF ( .not.LAMBERTIAN_FLAG ) THEN
        DO Q = 1, LOCAL_BRDF_NPARS
          IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
            DO I = 1, NSTREAMS
              SUM = ZERO
              DO K = 1, NSTREAMS_BRDF
                SUM  = SUM + D_BRDFUNC_0(Q,I,K)*BRDF_AZMFAC(K)
              ENDDO
              D_LOCAL_BRDF_F_0(Q,I) = SUM * HELP
            ENDDO
          ENDIF
        ENDDO
      ELSE
        D_LOCAL_BRDF_F_0 = ZERO
      ENDIF

!  Incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO Q = 1, LOCAL_BRDF_NPARS
          IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
            DO I = 1, NSTREAMS
              DO J = 1, NSTREAMS
                SUM = ZERO
                DO K = 1, NSTREAMS_BRDF
                  SUM  = SUM + D_BRDFUNC(Q,I,J,K) * BRDF_AZMFAC(K)
                ENDDO
                D_LOCAL_BRDF_F(Q,I,J) = SUM * HELP
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam (direct beam reflections)

        IF ( .not. LAMBERTIAN_FLAG ) THEN
           DO Q = 1, LOCAL_BRDF_NPARS
              IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
                 DO UI = 1, N_USER_STREAMS
                    SUM = ZERO
                    DO K = 1, NSTREAMS_BRDF
                      SUM = SUM + D_USER_BRDFUNC_0(Q,UI,K)*BRDF_AZMFAC(K)
                    ENDDO
                    D_LOCAL_USER_BRDF_F_0(Q,UI) = SUM * HELP
                 ENDDO
              ENDIF
           ENDDO
        ELSE
          D_LOCAL_USER_BRDF_F_0 = zero
        ENDIF

!  Incident quadrature directions (surface multiple reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO Q = 1, LOCAL_BRDF_NPARS
            IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
              DO UI = 1, N_USER_STREAMS
                DO J = 1, NSTREAMS
                  SUM = ZERO
                  DO K = 1, NSTREAMS_BRDF
                    SUM = SUM + D_USER_BRDFUNC(Q,UI,J,K)*BRDF_AZMFAC(K)
                  ENDDO
                  D_LOCAL_USER_BRDF_F(Q,UI,J) = SUM * HELP
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_LIN_FOURIER

!  End module

      END MODULE brdf_LinSup_routines_m

