
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
! #   SUBROUTINES (Everything public):                          #
! #                                                             #
! #             L_SURFACE_HOMOG_SETUP                           #
! #             L_SURFACE_BEAM_SETUP                            #
! #             LS_SURFACE_BEAM_SETUP                           #
! #                                                             #
! #             LC_BVPCOLUMN_SETUP                              #
! #             LPE_BVPCOLUMN_SETUP    (Elastic)                #
! #             LPR_BVPCOLUMN_SETUP    (Raman)                  #
! #             LSE_BVPCOLUMN_SETUP    (Elastic)                #
! #             LSR_BVPCOLUMN_SETUP    (Raman)                  #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are
!    (1) Bookkeeping improvements (use of "Only", clearer I/O specifications)
!    (2) Specific surface input replaced by general supplement-derived BRDF stuff

      MODULE lrrs_L_bvproblem_m

!      USE LRRS_PARS_m, Only : LDU

      PRIVATE
      PUBLIC :: L_SURFACE_HOMOG_SETUP, &
                L_SURFACE_BEAM_SETUP,  &
                LS_SURFACE_BEAM_SETUP, &
                LC_BVPCOLUMN_SETUP,    &
                LPE_BVPCOLUMN_SETUP,   &
                LPR_BVPCOLUMN_SETUP,   &
                LSE_BVPCOLUMN_SETUP,   &
                LSR_BVPCOLUMN_SETUP

      CONTAINS

      SUBROUTINE L_SURFACE_HOMOG_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR,  & ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX,      & ! Inputs
           DOVARY, NPARS, QUAD_STRMWGHT, ALBEDOS_RANKED, BRDF_F, & ! Inputs
           KTRANS, L_KTRANS, XPOS, L_XPOS, L_XNEG,               & ! Inputs
           L_R2_HOMP, L_R2_HOMM )                                  ! Outputs

!  Additional sums for the final surface-reflecting layer

!   Version 2.0 Lambertian only for the RRS model
!   Version 2.1 BRDF treatment introduced. 3 November 2007 by R. Spurr
!   Version 2.2 Treatment standardized. 21 November 2008 by R. Spurr
!   Version 2.5 New BRDF Treatment using supplement material, 10/10/15 by R. Spurr

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_MOMENTS, MAX_STREAMS, MAX_2_STREAMS, &
                              MAX_LAYERS, MAX_ATMOSWFS, MAX_POINTS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  inclusion of a surface

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  BRDF flag (replaces Lamabertian flag, Version 2.5, 9/11/15)

      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE

!  Fourier number, number of layers and streams

      INTEGER  , INTENT(IN) :: FOURIER, NLAYERS, NSTREAMS

!  Indices, Version 2.5, 9/11/15

      INTEGER  , INTENT(IN) :: Raman_IDX, Brdf_IDX

!  2-deltam0 factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Linearization control

      LOGICAL  , INTENT(IN) :: DOVARY
      INTEGER  , INTENT(IN) :: NPARS

!  quadrature inputs

      REAL(FPK), INTENT(IN)  :: QUAD_STRMWGHT(MAX_STREAMS)

!  Albedos

      REAL(fpk), intent(in)  :: ALBEDOS_RANKED ( MAX_POINTS)

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS (MAX_2_STREAMS,MAX_STREAMS,MAX_LAYERS)

      REAL(FPK), INTENT(IN) :: L_XPOS (MAX_ATMOSWFS,MAX_2_STREAMS,MAX_STREAMS,MAX_LAYERS)
      REAL(FPK), INTENT(IN) :: L_XNEG (MAX_ATMOSWFS,MAX_2_STREAMS,MAX_STREAMS,MAX_LAYERS)

!  Whole-layer Eigenstream transmittances

      REAL(FPK), INTENT(IN) :: KTRANS   (MAX_STREAMS,MAX_LAYERS)
      REAL(FPK), INTENT(IN) :: L_KTRANS (MAX_ATMOSWFS,MAX_STREAMS,MAX_LAYERS)

!  Output arguments
!  ----------------

!  Linearized Reflected homogeneous solutions at ground

      REAL(FPK), INTENT(OUT) :: L_R2_HOMP(MAX_ATMOSWFS,MAX_STREAMS,MAX_STREAMS)
      REAL(FPK), INTENT(OUT) :: L_R2_HOMM(MAX_ATMOSWFS,MAX_STREAMS,MAX_STREAMS)

!  local variables
!  ---------------

      REAL(FPK) :: H_1, H_2, FRJ, HSP_U, HSM_U, FACTOR
      INTEGER   :: I, AA, J, N, Q

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        DO AA = 1, NSTREAMS
          L_R2_HOMP(1:NPARS,I,AA) = ZERO
          L_R2_HOMM(1:NPARS,I,AA) = ZERO
        ENDDO
      ENDDO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE .or. .not. DOVARY ) RETURN

!  Layer

      N = NLAYERS

!  For Lambertian reflectance, all streams are the same

      IF ( .not. DO_BRDF_SURFACE ) THEN
        FACTOR = SURFACE_FACTOR * ALBEDOS_RANKED(Raman_IDX)
        IF ( FOURIER .EQ. 0 ) THEN
           DO Q = 1, NPARS
            DO AA = 1, NSTREAMS
              H_1 = ZERO ; H_2 = ZERO
              DO J = 1, NSTREAMS
                HSP_U = L_XPOS(Q,J,AA,N) *   KTRANS(AA,N) + &
                          XPOS(J,AA,N)   * L_KTRANS(Q,AA,N)
                HSM_U = L_XNEG(Q,J,AA,N)
                H_1 = H_1 + QUAD_STRMWGHT(J) * HSP_U
                H_2 = H_2 + QUAD_STRMWGHT(J) * HSM_U
              ENDDO
              L_R2_HOMP(Q,1:NSTREAMS,AA) = H_1 * FACTOR
              L_R2_HOMM(Q,1:NSTREAMS,AA) = H_2 * FACTOR
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  For BRDF surface

      IF ( DO_BRDF_SURFACE ) THEN
        DO Q = 1, NPARS
          DO AA = 1, NSTREAMS
            DO I = 1, NSTREAMS
              H_1 = ZERO ; H_2 = ZERO
              DO J = 1, NSTREAMS
                FRJ = QUAD_STRMWGHT(J) * BRDF_F(FOURIER,I,J,Brdf_IDX)
                HSP_U = L_XPOS(Q,J,AA,N) *   KTRANS(AA,N) + &
                          XPOS(J,AA,N)   * L_KTRANS(Q,AA,N)
                HSM_U = L_XNEG(Q,J,AA,N)
                H_1 = H_1 + FRJ * HSP_U
                H_2 = H_2 + FRJ * HSM_U
              ENDDO
              L_R2_HOMP(Q,I,AA) = H_1 * SURFACE_FACTOR
              L_R2_HOMM(Q,I,AA) = H_2 * SURFACE_FACTOR
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE L_SURFACE_HOMOG_SETUP

!

      SUBROUTINE L_SURFACE_BEAM_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR,     & ! Inputs
           FOURIER, NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX, DOVARY, & ! Inputs
           NVARY, QUAD_STRMWGHT, ALBEDOS_RANKED, BRDF_F, L_WLOWER,  & ! Inputs
           L_R2_BEAM )                                                ! Output

!  Linearized reflectance of BOA downwelling field (Atmospheric Jacobian)

!  Version 2.0 Lambertian only for the RRS model
!  Version 2.1 BRDF treatment introduced. 3 November 2007 by R. Spurr
!  Version 2.2 Treatment standardized. 21 November 2008 by R. Spurr
!  Version 2.5 Surface input from supplement-derived BRDFs, 10/10/15, R. Spurr

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_MOMENTS, MAX_STREAMS, MAX_2_STREAMS, &
                              MAX_LAYERS, MAX_ATMOSWFS, MAX_POINTS, MAX_POINTS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  inclusion of surface

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  BRDF flag (replaces Lamabertian flag, Version 2.5, 9/11/15)

      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE

!  Fourier number, number of layers and streams

      INTEGER  , INTENT(IN) :: FOURIER, NLAYERS, NSTREAMS

!  Indices, Version 2.5, 9/11/15

      INTEGER  , INTENT(IN) :: Raman_IDX, Brdf_IDX

!  2-deltam0 factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Linearization control

      LOGICAL  , INTENT(IN) :: DOVARY
      INTEGER  , INTENT(IN) :: NVARY

!  quadrature inputs

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT(MAX_STREAMS)

!  Albedos

      REAL(fpk), intent(in)  :: ALBEDOS_RANKED ( MAX_POINTS)

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Linearized particular solution at lower boundary

      REAL(FPK), INTENT(IN)  :: L_WLOWER(MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS)
!      REAL(FPK), INTENT(IN) :: L_WLOWER(MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS_NK) !

!  Output arguments
!  ----------------

!  Linearized Reflected beam solution at ground

      REAL(FPK), INTENT(OUT) :: L_R2_BEAM(MAX_ATMOSWFS,MAX_STREAMS)

!  local variables
!  ---------------

      INTEGER   :: I, J, Q, N
      REAL(FPK) :: SUMR, FACTOR, REFL_B

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        L_R2_BEAM(1:NVARY,I) = ZERO
      ENDDO

!  Return with Zeroed values if flags not set

      IF ( .NOT. DO_INCLUDE_SURFACE .or. .not. DOVARY) RETURN

!  Surface layer

      N = NLAYERS

!  For Lambertian reflectance, all streams are the same

      IF ( .not. DO_BRDF_SURFACE ) THEN
        FACTOR = SURFACE_FACTOR * ALBEDOS_RANKED(Raman_IDX)
        IF ( FOURIER .EQ. 0 ) THEN
          DO Q = 1, NVARY
            SUMR = FACTOR * DOT_PRODUCT(QUAD_STRMWGHT(1:NSTREAMS),L_WLOWER(Q,1:NSTREAMS,N))
            L_R2_BEAM(Q,1:NSTREAMS) = SUMR
          ENDDO
        ENDIF
      ENDIF

!  For bidirectional reflecting surface

      IF ( DO_BRDF_SURFACE ) THEN
        DO Q = 1, NVARY
          DO I = 1, NSTREAMS
            REFL_B = ZERO
            DO J = 1, NSTREAMS
              REFL_B = REFL_B + BRDF_F(FOURIER,I,J,Brdf_IDX) * QUAD_STRMWGHT(J) * L_WLOWER(Q,J,N)
            ENDDO
            L_R2_BEAM(Q,I) = REFL_B * SURFACE_FACTOR
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE L_SURFACE_BEAM_SETUP

!

      SUBROUTINE LS_SURFACE_BEAM_SETUP &
         ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, SURFACE_FACTOR, FOURIER,  & ! Inputs
           NLAYERS, NSTREAMS, Raman_IDX, Brdf_IDX, NSVARY, QUAD_STRMWGHT, & ! Inputs
           ALBEDOS_RANKED, BRDF_F, LS_BRDF_F, WLOWER, LS_WLOWER,          & ! Inputs
           LS_R2_BEAM )                                           ! output

!  Linearized BOA downwelling field (w.r.t. surface albedo)

!  Additional sums for the final surface-reflecting layer

!   Version 2.0 Lambertian only for the RRS model
!   Version 2.1 BRDF treatment introduced. 3 November 2007 by R. Spurr
!   Version 2.2 Treatment standardized. 21 November 2008 by R. Spurr
!     Earlier Versions: Not considering the variation of BIREFLEC here.

!   Version 2.5 Treatment extended to include BRDF variations
!               Use of more general supplement-derived BRDF inputs.
!               Bug found with Lambertian case, 2/2/17..

!  This term is not present for the elastic field.
!  This routine is NOT CALLED in the elastic treatment.

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_MOMENTS, MAX_STREAMS, MAX_2_STREAMS, &
                              MAX_LAYERS, MAX_SURFACEWFS, MAX_POINTS


      IMPLICIT NONE

!  Input arguments
!  ---------------

!  inclusion of a surface

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  BRDF flag (replaces Lamabertian flag, Version 2.5, 9/11/15)

      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE

!  Fourier number, number of layers, streams

      INTEGER  , INTENT(IN) :: FOURIER, NLAYERS, NSTREAMS

!  Indices, Version 2.5, 9/11/15

      INTEGER  , INTENT(IN) :: Brdf_IDX, Raman_IDX 

!  Number of surface Jacobians

      INTEGEr  , INTENT(IN) :: NSVARY

!  2-deltam0 factor

      REAL(fpk), intent(in) :: SURFACE_FACTOR

!  quadrature inputs

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT(MAX_STREAMS)

!  Albedos. Nt needed. Returned. 2/2/17

      REAL(fpk), intent(in) :: ALBEDOS_RANKED ( MAX_POINTS)

!  Fourier components of BRDF, in the following order
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in) :: BRDF_F    ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LS_BRDF_F ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )

!  particular solution at lower boundary

      REAL(FPK), INTENT(IN) :: WLOWER(MAX_2_STREAMS,MAX_LAYERS)

!  Linearized particular solution at lower boundary

      REAL(FPK), INTENT(IN) :: LS_WLOWER(MAX_SURFACEWFS,MAX_2_STREAMS,MAX_LAYERS)

!  Output arguments
!  ----------------

!  Linearized beam solution at ground, including reflected part

      REAL(FPK), INTENT(OUT) :: LS_R2_BEAM(MAX_SURFACEWFS,MAX_STREAMS)

!  local variables
!  ---------------

      INTEGER   :: I, I1, J, LAY, Q
      REAL(FPK) :: REFL_B, FACTOR

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        LS_R2_BEAM(1:NSVARY,I) = ZERO
      ENDDO

!  Return with Zeroed values if surface flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Surface layer

      LAY = NLAYERS

!  For Lambertian reflectance, all streams are the same
!  ----------------------------------------------------

!  2/2/17. debugged treatment.

!  Bad old code
!      IF ( .not. DO_BRDF_SURFACE ) THEN
!        FACTOR = SURFACE_FACTOR
!        IF ( FOURIER .EQ. 0 ) THEN
!          DO Q = 1, NSVARY
!            REFL_B = FACTOR * DOT_PRODUCT(QUAD_STRMWGHT(1:NSTREAMS),LS_WLOWER(Q,1:NSTREAMS,LAY))      !   Wrong
!            DO I = 1, NSTREAMS
!              I1 = I + NSTREAMS
!              LS_R2_BEAM(Q,I) = LS_WLOWER(Q,I1,LAY) - REFL_B * SURFACE_FACTOR
!            ENDDO
!            LS_R2_BEAM(Q,1:NSTREAMS) = REFL_B   !     Wrong
!          ENDDO
!        ENDIF
!      ENDIF

      IF ( .not. DO_BRDF_SURFACE ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO Q = 1, NSVARY
            REFL_B = ZERO
            DO J = 1, NSTREAMS
              REFL_B = REFL_B + ( WLOWER(J,LAY) + ALBEDOS_RANKED(Raman_IDX) * LS_WLOWER(Q,J,LAY) ) * QUAD_STRMWGHT(J)
            ENDDO
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              LS_R2_BEAM(Q,I) = LS_WLOWER(Q,I1,LAY) - REFL_B * SURFACE_FACTOR
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  For bidirectional reflecting surface
!  ------------------------------------

      IF ( DO_BRDF_SURFACE ) THEN
        DO Q = 1, NSVARY
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            REFL_B = ZERO
            DO J = 1, NSTREAMS
              REFL_B = REFL_B + ( LS_BRDF_F(Q,FOURIER,I,J,Brdf_IDX) *    WLOWER(J,LAY) &
                                   + BRDF_F(FOURIER,I,J,Brdf_IDX)   * LS_WLOWER(Q,J,LAY) ) * QUAD_STRMWGHT(J)
            ENDDO
            LS_R2_BEAM(Q,I) = LS_WLOWER(Q,I1,LAY) - REFL_B * SURFACE_FACTOR
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LS_SURFACE_BEAM_SETUP

!

      SUBROUTINE LC_BVPCOLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, KPARS, NTOTAL, NSTR2,    & ! Inputs
             LCON, MCON, XPOS, XNEG, L_XPOS, L_XNEG, KTRANS, L_KTRANS,       & ! Inputs
             L_WUPPER, L_WLOWER, L_R2_HOMP, L_R2_HOMM, L_R2_BEAM, L_DIRBEAM, & ! Inputs
             COL2WF, SCOL2WF )                                                 ! Outputs

!  Version 2.5, 10/11/15. Added SCOL2. Output declared "InOut" for OpenMP

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_TOTAL, MAX_ATMOSWFS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Inclusion of an albedo term

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  Computational variables

      INTEGER  , INTENT(IN) :: NSTREAMS, NLAYERS
      INTEGER  , INTENT(IN) :: NTOTAL, NSTR2

!  Linearization control

      INTEGER  , INTENT(IN) :: KPARS

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Linearized Eigenvector solutions

      REAL(FPK), INTENT(IN) :: L_XPOS ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XNEG ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Eigentransmittances and their linearizations

      REAL(FPK), INTENT(IN) :: KTRANS   ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  Reflected solutions at ground

      REAL(FPK), INTENT(IN) :: L_R2_HOMP ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: L_R2_HOMM ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )

!  Constants of integration

      REAL(FPK), INTENT(IN) :: LCON ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: MCON ( MAX_STREAMS, MAX_LAYERS )

!  Direct beam

      REAL(FPK), INTENT(IN) :: L_DIRBEAM(MAX_ATMOSWFS,MAX_STREAMS)

!  Beam solutions at upper and lower layer boundaries

      REAL(FPK), INTENT(IN) :: L_WUPPER(MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS)
      REAL(FPK), INTENT(IN) :: L_WLOWER(MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS)

!  Surface-reflected beam solution

      REAL(FPK), INTENT(IN) :: L_R2_BEAM(MAX_ATMOSWFS,MAX_STREAMS)

!  output arguments
!  ----------------

!  Column vectors

      REAL(FPK), INTENT(INOUT) :: COL2WF (MAX_TOTAL,MAX_ATMOSWFS)
      REAL(FPK), INTENT(INOUT) :: SCOL2WF(MAX_2_STREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      INTEGER   :: I,N,N1,I1,C0,CM,Q,AA
      REAL(FPK) :: COL2_D (MAX_TOTAL,MAX_ATMOSWFS)
      REAL(FPK) :: CPOS, CNEG, L_HOM, L_BEAM, L_HOMD, L_HOMU

!  zero column vector

      DO I = 1, NTOTAL
        DO Q = 1, KPARS
          COL2WF(I,Q) = ZERO
          COL2_D(I,Q) = ZERO
        ENDDO
      ENDDO

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------

      N  = 1

!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

      DO I = 1, NSTREAMS
        DO Q = 1, KPARS
          L_BEAM = - L_WUPPER(Q,I,N)
          L_HOM  = ZERO
          DO AA = 1, NSTREAMS
            CPOS = L_XPOS(Q,I,AA,N)
            CNEG =   KTRANS(AA,N)   * L_XNEG(Q,I,AA,N) + &
                   L_KTRANS(Q,AA,N) *   XNEG(I,AA,N)
            L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
          ENDDO
          COL2_D(I,Q) = L_BEAM - L_HOM
        ENDDO
      ENDDO

!  intermediate boundary conditions
!  --------------------------------

      DO N = 1, NLAYERS - 1

!  N1 is the layer below, C0 is the offset

        N1 = N + 1
        C0 = N*NSTR2 - NSTREAMS

!  .. 2 contributions to L_BEAM, from variations L_WUPPER and L_WLOWER
!  .. 2 contributions to L_HOM,  from variations above and below

        DO I = 1, NSTR2
          CM = C0 + I
          DO Q = 1, KPARS
            L_HOMD = ZERO
            L_HOMU = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(Q,I,AA,N1)
              CNEG =   KTRANS(AA,N1)   * L_XNEG(Q,I,AA,N1) + &
                     L_KTRANS(Q,AA,N1) *   XNEG(I,AA,N1)
              L_HOMU = L_HOMU + LCON(AA,N1) * CPOS + MCON(AA,N1) * CNEG
              CNEG = L_XNEG(Q,I,AA,N)
              CPOS =   KTRANS(AA,N)   * L_XPOS(Q,I,AA,N) + &
                     L_KTRANS(Q,AA,N) *   XPOS(I,AA,N)
              L_HOMD = L_HOMD + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            L_BEAM       = L_WUPPER(Q,I,N1) - L_WLOWER(Q,I,N)
            L_HOM        = L_HOMU - L_HOMD
            COL2_D(CM,Q) = L_BEAM + L_HOM
          ENDDO
        ENDDO

!  End layer loop

      ENDDO

!  lowest (surface) reflecting boundary (diffuse radiation terms only)
!  -------------------------------------------------------------------

      N  = NLAYERS
      C0 = (N-1)*NSTR2 + NSTREAMS

!  with non-zero surface reflectance, need integrated downward reflectan

      IF ( DO_INCLUDE_SURFACE ) THEN

        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, KPARS
            L_BEAM = L_WLOWER(Q,I1,N) - L_R2_BEAM(Q,I)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS =   KTRANS(AA,N)   * L_XPOS(Q,I1,AA,N) + &
                     L_KTRANS(Q,AA,N) *   XPOS(I1,AA,N)
              CPOS =        CPOS       - L_R2_HOMP(Q,I,AA)
              CNEG = L_XNEG(Q,I1,AA,N) - L_R2_HOMM(Q,I,AA)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COL2_D(CM,Q) = - L_BEAM - L_HOM
          ENDDO
        ENDDO

!  no surface, similar code excluding integrated reflectance

      ELSE

        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, KPARS
            L_BEAM = L_WLOWER(Q,I1,N)
            L_HOM    = ZERO
            DO AA = 1, NSTREAMS
              CPOS =   KTRANS(AA,N)   * L_XPOS(Q,I1,AA,N) + &
                     L_KTRANS(Q,AA,N) *   XPOS(I1,AA,N)
              CNEG = L_XNEG(Q,I1,AA,N)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COL2_D(CM,Q) = - L_BEAM - L_HOM
          ENDDO
        ENDDO

      ENDIF

!  Add direct beam solution (only to final level)
!  ----------------------------------------

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          DO Q = 1, KPARS
            COL2_D(CM,Q) = COL2_D(CM,Q) + L_DIRBEAM(Q,I)
          ENDDO
        ENDDO
      ENDIF

!  Add to existing solution
!   Note 9/11/15. Why do we need this ???????

      DO I = 1, NTOTAL
        DO Q = 1, KPARS
          COL2WF(I,Q) = COL2WF(I,Q) + COL2_D(I,Q)
        ENDDO
      ENDDO

!  Return, full case

      IF ( NLAYERS .NE. 1 ) RETURN

!  If Nlayers = 1, special case. Just Copy to output

      DO I = 1, NTOTAL
        SCOL2WF(I,1:KPARS) = COL2WF(I,1:KPARS)
      ENDDO

!  finish

      RETURN
      END SUBROUTINE LC_BVPCOLUMN_SETUP

!

      SUBROUTINE LPE_BVPCOLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, K, KPARS, & ! Inputs
             NTOTAL, NSTR2, LCON, MCON, XPOS, XNEG, KTRANS,   & ! Inputs
             L_XPOS, L_XNEG, L_KTRANS, L_WUPPER, L_WLOWER,    & ! Inputs
             L_R2_HOMP, L_R2_HOMM, L_R2_BEAM, L_DIRECT_BEAM,  & ! Inputs
             COL2WF, SCOL2WF)                                   ! Outputs

!  Subroutine for the Profile linearization of the BVP column vector
!    Elastic sources only.............................

!  Version 2.5, 10/11/15. Added SCOL2. Output declared "InOut" for OpenMP

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_TOTAL, MAX_ATMOSWFS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Inclusion of an albedo term

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  Computational variables

      INTEGER  , INTENT(IN) :: NSTREAMS, NLAYERS
      INTEGER  , INTENT(IN) :: NTOTAL, NSTR2

!  Linearization control

      INTEGER  , INTENT(IN) :: K, KPARS

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Linearized Eigenvector solutions

      REAL(FPK), INTENT(IN) :: L_XPOS &
              ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XNEG &
              ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Eigentransmittances and their linearizations

      REAL(FPK), INTENT(IN) :: KTRANS             ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  Reflected solutions at ground

      REAL(FPK), INTENT(IN) :: L_R2_HOMP ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: L_R2_HOMM ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )

!  Constants of integration

      REAL(FPK), INTENT(IN) :: LCON ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: MCON ( MAX_STREAMS, MAX_LAYERS )

!  Direct beam

      REAL(FPK), INTENT(IN) :: L_DIRECT_BEAM(MAX_ATMOSWFS,MAX_STREAMS)

!  Beam solutions at upper and lower layer boundaries

      REAL(FPK), INTENT(IN) :: L_WUPPER(MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS)
      REAL(FPK), INTENT(IN) :: L_WLOWER(MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS)

!  Surface-reflected beam solution

      REAL(FPK), INTENT(IN) :: L_R2_BEAM(MAX_ATMOSWFS,MAX_STREAMS)

!  output arguments
!  ----------------

!  Column vector

      REAL(FPK), INTENT(INOUT) :: COL2WF (MAX_TOTAL,MAX_ATMOSWFS)
      REAL(FPK), INTENT(INOUT) :: SCOL2WF(MAX_2_STREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      LOGICAL   :: MODIFIED_BCL3, MODIFIED_BCL4
      LOGICAL   :: REGULAR_BCL3,  REGULAR_BCL4
      INTEGER   :: I,N,N1,I1,C0,CM,Q,AA
      REAL(FPK) :: CPOS, CNEG, L_HOM, L_BEAM

!  Boundary condition flags for special cases

      MODIFIED_BCL3 = ( K .EQ. 1 )
      MODIFIED_BCL4 = ( K .EQ. NLAYERS )

!  zero column vector

      DO I = 1, NTOTAL
        DO Q = 1, KPARS
          COL2WF(I,Q) = ZERO
        ENDDO
      ENDDO

!  complete boundary condition flags

      REGULAR_BCL3 = .NOT.MODIFIED_BCL3
      REGULAR_BCL4 = .NOT.MODIFIED_BCL4

!  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
!  ------------------------------------------------------------

!    If this layer is the one that is varied, use MODIFIED_BCL3 (BCL3M)
!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations
!     Otherwise  zero.

      N  = 1
      IF ( MODIFIED_BCL3 ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, KPARS
            L_BEAM = - L_WUPPER(Q,I,N)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(Q,I,AA,N)
              CNEG =   KTRANS(AA,N)   * L_XNEG(Q,I,AA,N) + &
                     L_KTRANS(Q,AA,N) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            COL2WF(I,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO
      ENDIF

!  BCL2 Intermediate levels between top layer and varying layer
!  ------------------------------------------------------

!  YOU COULD GET RID OF THIS CODE...............
!  [not required if top layer is varying, case MODIFIED_BCL3 above]
!  .. nothing varying in these layers
!      IF ( REGULAR_BCL3 ) THEN
!        DO N = 2, K - 1
!          N1 = N - 1
!          C0  = N1*NSTR2 - NSTREAMS
!          DO I = 1, NSTR2
!            CM = C0 + I
!            DO Q = 1, KPARS
!              COL2WF(CM,Q) = ZERO
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDIF

!  BCL3 - regular upper boundary condition for layer that is varying
!  -----------------------------------------------------------

      IF ( REGULAR_BCL3 ) THEN

        N   = K
        N1  = N - 1
        C0  = N1*NSTR2 - NSTREAMS

!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO I = 1, NSTR2
          CM = C0 + I
          DO Q = 1, KPARS
            L_BEAM = + L_WUPPER(Q,I,N)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(Q,I,AA,N)
              CNEG =   KTRANS(AA,N)   * L_XNEG(Q,I,AA,N) + &
                     L_KTRANS(Q,AA,N) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            COL2WF(CM,Q) = L_BEAM + L_HOM
          ENDDO
        ENDDO

      ENDIF

!  BCL4 - LOWER boundary condition for varying layer
!  -------------------------------------------

!   special case when layer-to-vary = last (albedo) layer is treated
!   separately below under MODIFIED BCL4.
!  .. 2 contributions to WVAR from beam solution variations BEAM_V and B
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

      IF ( REGULAR_BCL4 ) THEN
        N   = K
        N1  = N + 1
        C0 = N*NSTR2 - NSTREAMS
        DO I = 1, NSTR2
          CM = C0 + I
          DO Q = 1, KPARS
            L_BEAM = L_WUPPER(Q,I,N1) - L_WLOWER(Q,I,N)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CNEG = L_XNEG(Q,I,AA,N)
              CPOS =   KTRANS(AA,N)   * L_XPOS(Q,I,AA,N) + &
                     L_KTRANS(Q,AA,N) *   XPOS(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COL2WF(CM,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO
      ENDIF

!  BCL5 - Intermediate boundary conditions between varying layer & final
!  ---------------------------------------------------------------

!  .. contributions from beam solution (direct assign). No homog. variat

      IF ( REGULAR_BCL4 ) THEN
        DO N = K + 1, NLAYERS - 1
          N1  = N + 1
          C0  = N*NSTR2 - NSTREAMS
          DO I = 1, NSTR2
            CM = C0 + I
            DO Q = 1, KPARS
              COL2WF(CM,Q) = L_WUPPER(Q,I,N1) - L_WLOWER(Q,I,N)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Final layer - use BCL6 or BCL4M (last layer is varying)
!  -------------------------------------------------

      N  = NLAYERS

!  Modified BCL4M Component loop

      IF ( MODIFIED_BCL4 ) THEN
        C0 = (N-1)*NSTR2 + NSTREAMS
        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, KPARS
            L_BEAM = L_WLOWER(Q,I1,N) - L_R2_BEAM(Q,I)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = KTRANS(AA,N)   * L_XPOS(Q,I1,AA,N) + &
                   L_KTRANS(Q,AA,N) *   XPOS(I1,AA,N)
              CPOS =        CPOS       - L_R2_HOMP(Q,I,AA)
              CNEG = L_XNEG(Q,I1,AA,N) - L_R2_HOMM(Q,I,AA)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COL2WF(CM,Q) = - L_BEAM - L_HOM
         ENDDO
        ENDDO

!  ordinary BCL6 Component loop

      ELSE

!  Compute the solution.

        C0 = (N-1)*NSTR2 + NSTREAMS
          DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, KPARS
            L_BEAM = L_WLOWER(Q,I1,N) - L_R2_BEAM(Q,I)
            COL2WF(CM,Q) = - L_BEAM
          ENDDO
        ENDDO

      ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------

!  Don't think we need this code......................
!      IF ( DO_INCLUDE_SURFACE ) THEN
!        DO I = 1, NSTREAMS
!          CM = C0 + I
!          FAC = - DIRECT_BEAM(I) * DELTAU_SLANT(N,K)
!          DO Q = 1, KPARS
!            L_BEAM = L_DELTAU_VERT(Q,K) * FAC
!            COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
!          ENDDO
!        ENDDO
!      ENDIF

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          DO Q = 1, KPARS
            COL2WF(CM,Q) = COL2WF(CM,Q) + L_DIRECT_BEAM(Q,I)
          ENDDO
        ENDDO
      ENDIF

!  Return, full case

      IF ( NLAYERS .NE. 1 ) RETURN

!  If Nlayers = 1, special case. Just Copy to output

      DO I = 1, NTOTAL
        SCOL2WF(I,1:KPARS) = COL2WF(I,1:KPARS)
      ENDDO

!  finish

      RETURN
      END SUBROUTINE LPE_BVPCOLUMN_SETUP

!

      SUBROUTINE LPR_BVPCOLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, K, KPARS, & ! Inputs
             NTOTAL, NSTR2, LCON, MCON, XPOS, XNEG, KTRANS,   & ! Inputs
             L_XPOS, L_XNEG,L_KTRANS, L_WUPPER, L_WLOWER,     & ! Inputs
             L_R2_HOMP, L_R2_HOMM, L_R2_BEAM, L_DIRECT_BEAM,  & ! Inputs
             COL2WF, SCOL2WF )                                  ! Outputs

!  Subroutine for the Profile linearization of the BVP column vector
!    Raman sources only.............................

!  Version 2.5, 10/11/15. Added SCOL2. Output declared "InOut" for OpenMP

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_TOTAL, MAX_ATMOSWFS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Inclusion of an albedo term

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  Computational variables

      INTEGER  , INTENT(IN) :: NSTREAMS, NLAYERS
      INTEGER  , INTENT(IN) :: NTOTAL, NSTR2

!  Linearization control

      INTEGER  , INTENT(IN) :: K, KPARS

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Linearized Eigenvector solutions

      REAL(FPK), INTENT(IN) :: L_XPOS &
              ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XNEG &
              ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Eigentransmittances and their linearizations

      REAL(FPK), INTENT(IN) :: KTRANS             ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  Reflected solutions at ground

      REAL(FPK), INTENT(IN) :: L_R2_HOMP ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: L_R2_HOMM ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )

!  Constants of integration

      REAL(FPK), INTENT(IN) :: &
            LCON ( MAX_STREAMS, MAX_LAYERS ), &
            MCON ( MAX_STREAMS, MAX_LAYERS )

!  Direct beam

      REAL(FPK), INTENT(IN) :: L_DIRECT_BEAM(MAX_ATMOSWFS,MAX_STREAMS)

!  Beam solutions at upper and lower layer boundaries

      REAL(FPK), INTENT(IN) :: L_WUPPER(MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS)
      REAL(FPK), INTENT(IN) :: L_WLOWER(MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS)

!  Surface-reflected beam solution

      REAL(FPK), INTENT(IN) :: L_R2_BEAM(MAX_ATMOSWFS,MAX_STREAMS)

!  output arguments
!  ----------------

!  Column vector

      REAL(FPK), INTENT(INOUT) :: COL2WF  (MAX_TOTAL,MAX_ATMOSWFS)
      REAL(FPK), INTENT(INOUT) :: SCOL2WF (MAX_2_STREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      LOGICAL   :: MODIFIED_BCL3, MODIFIED_BCL4
      LOGICAL   :: REGULAR_BCL3,  REGULAR_BCL4
      INTEGER   :: I,N,N1,I1,C0,CM,Q,AA
      REAL(FPK) :: CPOS, CNEG, L_HOM, L_BEAM

!  Boundary condition flags for special cases

      MODIFIED_BCL3 = ( K .EQ. 1 )
      MODIFIED_BCL4 = ( K .EQ. NLAYERS )

!  zero column vector

      DO I = 1, NTOTAL
        COL2WF(I,1:KPARS) = ZERO
      ENDDO

!  complete boundary condition flags

      REGULAR_BCL3 = .NOT.MODIFIED_BCL3
      REGULAR_BCL4 = .NOT.MODIFIED_BCL4

!  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
!  ------------------------------------------------------------

!    If this layer is the one that is varied, use MODIFIED_BCL3 (BCL3M)
!  .. contribution WVAR from beam solution variations
!  .. contribution HVAR homogeneous (eigenvalue) solution variations
!     Otherwise Particular integral contributions only (NON-ZERO HERE)

      N  = 1
      IF ( MODIFIED_BCL3 ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, KPARS
            L_BEAM = - L_WUPPER(Q,I,N)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(Q,I,AA,N)
              CNEG =   KTRANS(AA,N)   * L_XNEG(Q,I,AA,N) + &
                     L_KTRANS(Q,AA,N) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            COL2WF(I,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO
      ELSE
        DO I = 1, NSTREAMS
          DO Q = 1, KPARS
            COL2WF(I,Q) = - L_WUPPER(Q,I,N)
          ENDDO
        ENDDO
      ENDIF

!  BCL2 Intermediate levels between top layer and varying layer
!  ------------------------------------------------------

!  Only the particular integrals are varying (both layers)

      IF ( REGULAR_BCL3 ) THEN
        DO N = 2, K - 1
          N1 = N - 1
          C0  = N1*NSTR2 - NSTREAMS
          DO I = 1, NSTR2
            CM = C0 + I
            DO Q = 1, KPARS
              COL2WF(CM,Q) = + L_WUPPER(Q,I,N) - L_WLOWER(Q,I,N1)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  BCL3 - regular upper boundary condition for layer that is varying
!  -----------------------------------------------------------

      IF ( REGULAR_BCL3 ) THEN

        N   = K
        N1  = N - 1
        C0  = N1*NSTR2 - NSTREAMS

!  .. contribution WVAR from beam solution variations (both layers)
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO I = 1, NSTR2
          CM = C0 + I
          DO Q = 1, KPARS
            L_BEAM = + L_WUPPER(Q,I,N) - L_WLOWER(Q,I,N1)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = L_XPOS(Q,I,AA,N)
              CNEG =   KTRANS(AA,N)   * L_XNEG(Q,I,AA,N) + &
                     L_KTRANS(Q,AA,N) *   XNEG(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
            ENDDO
            COL2WF(CM,Q) = L_BEAM + L_HOM
          ENDDO
        ENDDO

      ENDIF

!  BCL4 - LOWER boundary condition for varying layer
!  -------------------------------------------

!   special case when layer-to-vary = last (albedo) layer is treated
!   separately below under MODIFIED BCL4.
!  .. 2 contributions to WVAR from beam solution variations BEAM_V and B
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

      IF ( REGULAR_BCL4 ) THEN
        N   = K
        N1  = N + 1
        C0 = N*NSTR2 - NSTREAMS
        DO I = 1, NSTR2
          CM = C0 + I
          DO Q = 1, KPARS
            L_BEAM = L_WUPPER(Q,I,N1) - L_WLOWER(Q,I,N)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CNEG = L_XNEG(Q,I,AA,N)
              CPOS =   KTRANS(AA,N)   * L_XPOS(Q,I,AA,N) + &
                     L_KTRANS(Q,AA,N) *   XPOS(I,AA,N)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COL2WF(CM,Q) = L_BEAM - L_HOM
          ENDDO
        ENDDO
      ENDIF

!  BCL5 - Intermediate boundary conditions between varying layer & final
!  ---------------------------------------------------------------

!  .. contributions from beam solution (direct assign). No homog. variat

      IF ( REGULAR_BCL4 ) THEN
        DO N = K + 1, NLAYERS - 1
          N1  = N + 1
          C0  = N*NSTR2 - NSTREAMS
          DO I = 1, NSTR2
            CM = C0 + I
            DO Q = 1, KPARS
              COL2WF(CM,Q) = L_WUPPER(Q,I,N1) - L_WLOWER(Q,I,N)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Final layer - use BCL6 or BCL4M (last layer is varying)
!  -------------------------------------------------

      N  = NLAYERS

!  Modified BCL4M Component loop

      IF ( MODIFIED_BCL4 ) THEN
        C0 = (N-1)*NSTR2 + NSTREAMS
        DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, KPARS
            L_BEAM = L_WLOWER(Q,I1,N) - L_R2_BEAM(Q,I)
            L_HOM  = ZERO
            DO AA = 1, NSTREAMS
              CPOS = KTRANS(AA,N)   * L_XPOS(Q,I1,AA,N) + &
                   L_KTRANS(Q,AA,N) *   XPOS(I1,AA,N)
              CPOS =        CPOS       - L_R2_HOMP(Q,I,AA)
              CNEG = L_XNEG(Q,I1,AA,N) - L_R2_HOMM(Q,I,AA)
              L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
            ENDDO
            COL2WF(CM,Q) = - L_BEAM - L_HOM
         ENDDO
        ENDDO

!  ordinary BCL6 Component loop

      ELSE

!  Compute the solution.

        C0 = (N-1)*NSTR2 + NSTREAMS
          DO I = 1, NSTREAMS
          CM = C0 + I
          I1 = I + NSTREAMS
          DO Q = 1, KPARS
            L_BEAM = L_WLOWER(Q,I1,N) - L_R2_BEAM(Q,I)
            COL2WF(CM,Q) = - L_BEAM
          ENDDO
        ENDDO

      ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------

!  Don't think we need this code......................
!      IF ( DO_INCLUDE_SURFACE ) THEN
!        DO I = 1, NSTREAMS
!          CM = C0 + I
!          FAC = - DIRECT_BEAM(I) * DELTAU_SLANT(N,K)
!          DO Q = 1, KPARS
!            L_BEAM = L_DELTAU_VERT(Q,K) * FAC
!            COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
!          ENDDO
!        ENDDO
!      ENDIF

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO I = 1, NSTREAMS
          CM = C0 + I
          DO Q = 1, KPARS
            COL2WF(CM,Q) = COL2WF(CM,Q) + L_DIRECT_BEAM(Q,I)
          ENDDO
        ENDDO
      ENDIF

!  Return, full case

      IF ( NLAYERS .NE. 1 ) RETURN

!  If Nlayers = 1, special case. Just Copy to output

      DO I = 1, NTOTAL
        SCOL2WF(I,1:KPARS) = COL2WF(I,1:KPARS)
      ENDDO

!  finish

      RETURN
      END SUBROUTINE LPR_BVPCOLUMN_SETUP

!

      SUBROUTINE LSE_BVPCOLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                & ! Inputs
             FOURIER, NLAYERS, NSTREAMS, NTOTAL, NSTR2, NSVARY,  & ! Inptus
             SURFACE_FACTOR, Brdf_IDX, QUAD_STRMWGHT, LS_BRDF_F, & ! Inputs
             LCON, MCON, XPOS, XNEG, KTRANS, WLOWER,             & ! Inputs
             R2_HOMP, R2_HOMM, R2_BEAM, LS_DIRECT_BEAM,          & ! Inputs
             COL2WF_SURF, SCOL2WF_SURF )                           ! Outputs

!  19 November 2008. First programming.
!  21 November 2008, Extend to more general linearization
!  08 June 2009,     Final form for the elastic solution

!  Only applies to the Elastic solution

!  Version 2.5, 10/11/15. Added SCOL2. Output declared "InOut" for OpenMP

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, PI4, &
                              MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_TOTAL, &
                              MAX_MOMENTS, MAX_POINTS, MAX_SURFACEWFS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Inclusion of an albedo term

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  BRDF flag (replaces Lamabertian flag, Version 2.5, 9/11/15)

      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE

!  Fourier number, number of layers and streams

      INTEGER  , INTENT(IN) :: FOURIER, NLAYERS, NSTREAMS
      INTEGER  , INTENT(IN) :: NTOTAL, NSTR2

!  Number of surface weighting functions

      INTEGER  , INTENT(IN) :: NSVARY

!  Indices, Version 2.5, 9/11/15

      INTEGER  , INTENT(IN) :: Brdf_IDX ! Raman_IDX

!  2-deltam0 factor

      REAL(fpk), INTENT(IN) :: SURFACE_FACTOR

!  Albedos. Not needed
!      REAL(fpk), intent(in) :: ALBEDOS_RANKED ( MAX_POINTS)

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), INTENT(IN) :: LS_BRDF_F ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, &
                                           MAX_STREAMS, MAX_POINTS )

!  quadrature inputs

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT ( MAX_STREAMS )

!  Constants of integration

      REAL(FPK), INTENT(IN) :: LCON ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: MCON ( MAX_STREAMS, MAX_LAYERS )

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Eigentransmittances

      REAL(FPK), INTENT(IN) :: KTRANS  ( MAX_STREAMS, MAX_LAYERS )

!  particular integral

      REAL(FPK), INTENT(IN) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )

!  Reflected solutions at ground

      REAL(FPK), INTENT(IN) :: R2_HOMP ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: R2_HOMM ( MAX_STREAMS, MAX_STREAMS )

!  Surface-reflected beam solution

      REAL(FPK), INTENT(IN) :: R2_BEAM ( MAX_STREAMS )

!  Surface-linearized Direct beam

      REAL(FPK), INTENT(IN) :: LS_DIRECT_BEAM ( MAX_SURFACEWFS, MAX_STREAMS )

!  Output arguments
!  ----------------

!  Column vector

      REAL(FPK), INTENT(INOUT) :: COL2WF_SURF  ( MAX_TOTAL,     MAX_SURFACEWFS )
      REAL(FPK), INTENT(INOUT) :: SCOL2WF_SURF ( MAX_2_STREAMS, MAX_SURFACEWFS )

!  Local variables
!  ---------------

      INTEGER   :: I, J, M, N, Q, AA, C0, CM
      REAL(FPK) :: HELP, REFL_B, REFL_P, REFL_M, FACTOR, FRAC1, F1
      REAL(FPK) :: H_XPOS   ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: H_XNEG   ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: H_WLOWER ( MAX_STREAMS )

!  Zero column vector

      DO I = 1, NTOTAL
        COL2WF_SURF(I,1:NSVARY) = ZERO
      ENDDO

!  Return if no surface

      if ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Ground level boundary condition
!  -------------------------------

!  Initialise

      M  = FOURIER
      N  = NLAYERS
      C0 = (N-1)*NSTR2 + NSTREAMS

!mick mod 7/20/2016 - moved these calculations from if blocks below
!                     (applied to both Lambertian & BRDF cases)
!  Compute some intermediate quantities

      DO J = 1, NSTREAMS
        H_WLOWER(J) = QUAD_STRMWGHT(J) * WLOWER(J,N)
        DO AA = 1, NSTREAMS
          H_XPOS(J,AA) = QUAD_STRMWGHT(J) * XPOS(J,AA,N)
          H_XNEG(J,AA) = QUAD_STRMWGHT(J) * XNEG(J,AA,N)
        ENDDO
      ENDDO

!  Diffuse scatter contributions (Lambertian case)

      IF ( .not. DO_BRDF_SURFACE ) then
        IF ( M .eq. 0 ) then

!         FRAC1 = R2 * FL1                !@@@@@@@ ROB debug 4/1/11. Wrong
!         FRAC1 = 1.0d0                   !@@@@@@@ ROB debug 4/1/11. Right

!mick note 7/20/2016 - R2_BEAM here = LIDORT's REFL_B  * SURFACE_FACTOR * ALBEDOS_RANKED(Raman_IDX)
!                        (see lrrs_bvproblem.f90 ~line 266)
!mick note 7/20/2016 - R2_HOMP here = LIDORT's REFL_HP * SURFACE_FACTOR * ALBEDOS_RANKED(Raman_IDX) 
!                      R2_HOMM here = LIDORT's REFL_HM * SURFACE_FACTOR * ALBEDOS_RANKED(Raman_IDX)
!                        (see lrrs_bvproblem.f90 ~line 129)
!mick fix 7/20/2016  - replaced this stream loop with new code below
          !DO I = 1, NSTREAMS
          !  CM = C0 + I
          !  HELP = R2_BEAM(I)
          !  DO AA = 1, NSTREAMS
          !    HELP = HELP + LCON(AA,N) * R2_HOMP(I,AA) * KTRANS(AA,N) &
          !                + MCON(AA,N) * R2_HOMM(I,AA)
          !  ENDDO
          !  FACTOR = SURFACE_FACTOR * HELP
          !  COL2WF_SURF(CM,1:NSVARY) = FACTOR
          !ENDDO

!mick fix 7/20/2016 - new code (now more LIDORT-like & similar to BRDF code below)

          DO I = 1, NSTREAMS
            CM = C0 + I
            REFL_B = SUM(H_WLOWER(1:NSTREAMS))
            HELP   = REFL_B
            DO AA = 1, NSTREAMS
              REFL_P = SUM(H_XPOS(1:NSTREAMS,AA))
              REFL_M = SUM(H_XNEG(1:NSTREAMS,AA))
              HELP = HELP + LCON(AA,N) * REFL_P * KTRANS(AA,N) &
                          + MCON(AA,N) * REFL_M
            ENDDO
            COL2WF_SURF(CM,1:NSVARY) = SURFACE_FACTOR * HELP

          ENDDO
        ENDIF
      ENDIF

!  Diffuse scatter contributions (BRDF case)

      IF ( DO_BRDF_SURFACE ) then
!mick fix 7/20/2016 - swapped indices (1:NSTREAMS,I) ---> (I,1:NSTREAMS) on LS_BRDF_F array.
!                     1:NSTREAMS is the usual J index where J is incident.
        DO I = 1, NSTREAMS
          CM = C0 + I
          DO Q = 1, NSVARY
            !REFL_B = DOT_PRODUCT(H_WLOWER(1:NSTREAMS),LS_BRDF_F(Q,M,1:NSTREAMS,I,BRDF_idx))
            REFL_B = DOT_PRODUCT(H_WLOWER(1:NSTREAMS),LS_BRDF_F(Q,M,I,1:NSTREAMS,BRDF_idx))
            HELP = REFL_B

            DO AA = 1, NSTREAMS
              !REFL_P = DOT_PRODUCT(H_XPOS(1:NSTREAMS,AA),LS_BRDF_F(Q,M,1:NSTREAMS,I,BRDF_idx))
              !REFL_M = DOT_PRODUCT(H_XNEG(1:NSTREAMS,AA),LS_BRDF_F(Q,M,1:NSTREAMS,I,BRDF_idx))
              REFL_P = DOT_PRODUCT(H_XPOS(1:NSTREAMS,AA),LS_BRDF_F(Q,M,I,1:NSTREAMS,BRDF_idx))
              REFL_M = DOT_PRODUCT(H_XNEG(1:NSTREAMS,AA),LS_BRDF_F(Q,M,I,1:NSTREAMS,BRDF_idx))

              HELP = HELP + LCON(AA,N) * REFL_P * KTRANS(AA,N) &
                          + MCON(AA,N) * REFL_M
            ENDDO
            COL2WF_SURF(CM,Q) = SURFACE_FACTOR * HELP
          ENDDO
        ENDDO
      ENDIF

!  Add direct beam variation

      DO I = 1, NSTREAMS
        CM = C0 + I
        COL2WF_SURF(CM,1:NSVARY) = COL2WF_SURF(CM,1:NSVARY) + LS_DIRECT_BEAM(1:NSVARY,I)
      ENDDO

!  Return, full case

      IF ( NLAYERS .NE. 1 ) RETURN

!  If Nlayers = 1, special case. Just Copy to output

      DO I = 1, NTOTAL
        SCOL2WF_SURF(I,1:NSVARY) = COL2WF_SURF(I,1:NSVARY)
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSE_BVPCOLUMN_SETUP

!

      SUBROUTINE LSR_BVPCOLUMN_SETUP &
           ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                & ! Inputs
             FOURIER, NLAYERS, NSTREAMS, NTOTAL, NSTR2, NSVARY,  & ! Inptus
             SURFACE_FACTOR, Brdf_IDX, QUAD_STRMWGHT, LS_BRDF_F, & ! Inputs
             LCON, MCON, XPOS, XNEG, KTRANS, WLOWER,             & ! Inputs
             R2_HOMP, R2_HOMM, R2_BEAM, LS_WLOWER, LS_WUPPER,    & ! Inputs
             LS_DIRECT_BEAM, LS_DIFFUSE_BEAM,                    & ! Inputs
             COL2WF_SURF, SCOL2WF_SURF )                           ! Outputs

!  19 November 2008. First programming.
!  21 November 2008, Extend to more general linearization

!  Version 2.5, 10/11/15. Added SCOL2. Output declared "InOut" for OpenMP
!               02/02/17. Slight alteration of treatment of Diffuse Beam terms

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS, MAX_TOTAL, &
                              MAX_MOMENTS, MAX_POINTS, MAX_SURFACEWFS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Inclusion of an albedo term

      LOGICAL  , INTENT(IN) :: DO_INCLUDE_SURFACE

!  BRDF flag (replaces Lamabertian flag, Version 2.5, 9/11/15)

      LOGICAL  , INTENT(IN) :: DO_BRDF_SURFACE

!  Fourier number, number of layers and streams

      INTEGER  , INTENT(IN) :: FOURIER, NLAYERS, NSTREAMS
      INTEGER  , INTENT(IN) :: NTOTAL, NSTR2

!  Number of surface weighting functions

      INTEGER  , INTENT(IN) :: NSVARY

!  Indices, Version 2.5, 9/11/15

      INTEGER  , INTENT(IN) :: Brdf_IDX ! Raman_idx

!  2-deltam0 factor

      REAL(fpk), intent(in) :: SURFACE_FACTOR

!  Albedos. Not needed
!      REAL(fpk), intent(in) :: ALBEDOS_RANKED ( MAX_POINTS)

!  Fourier components of BRDF, in the following order (same all threads)
!    incident quadrature streams, reflected quadrature streams

!      REAL(fpk), intent(in) :: BRDF_F    ( 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), intent(in) :: LS_BRDF_F ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Quadrature inputs

      REAL(FPK), INTENT(IN) :: QUAD_STRMWGHT(MAX_STREAMS)

!  Constants of integration

      REAL(FPK), INTENT(IN) :: LCON ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: MCON ( MAX_STREAMS, MAX_LAYERS )

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Eigentransmittances

      REAL(FPK), INTENT(IN) :: KTRANS  ( MAX_STREAMS, MAX_LAYERS )

!  Particular integral

      REAL(FPK), INTENT(IN) :: WLOWER ( MAX_2_STREAMS, MAX_LAYERS )

!  Reflected solutions at ground

      REAL(FPK), INTENT(IN) :: R2_HOMP ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: R2_HOMM ( MAX_STREAMS, MAX_STREAMS )

!  Surface-reflected beam solution

      REAL(FPK), INTENT(IN) :: R2_BEAM(MAX_STREAMS)

!  Linearized particular solutions at layer boundaries

      REAL(FPK), INTENT(IN) :: LS_WLOWER(MAX_SURFACEWFS,MAX_2_STREAMS,MAX_LAYERS)
      REAL(FPK), INTENT(IN) :: LS_WUPPER(MAX_SURFACEWFS,MAX_2_STREAMS,MAX_LAYERS)

!  Surface-linearized Direct beam

      REAL(FPK), INTENT(IN) :: LS_DIRECT_BEAM(MAX_SURFACEWFS,MAX_STREAMS)

!  Surface-linearized beam solution at surface (including reflected part
!   (Not required in the elastic case)

      REAL(FPK), INTENT(IN) :: LS_DIFFUSE_BEAM(MAX_SURFACEWFS,MAX_STREAMS)

!  Output arguments
!  ----------------

!  Column vector

      REAL(FPK), INTENT(INOUT) :: COL2WF_SURF (MAX_TOTAL,MAX_SURFACEWFS)
      REAL(FPK), INTENT(INOUT) :: SCOL2WF_SURF(MAX_2_STREAMS,MAX_SURFACEWFS)

!  Local variables
!  ---------------

      INTEGER   :: I, J, M, N, Q, AA, C0, CM, N1
      REAL(FPK) :: HELP, REFL_B, REFL_P, REFL_M
      REAL(FPK) :: H_XPOS   ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: H_XNEG   ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK) :: H_WLOWER ( MAX_STREAMS ), FACTOR

!  Zero column vector

      DO I = 1, NTOTAL
        COL2WF_SURF(I,1:NSVARY) = ZERO
      ENDDO

!  Return if no surface

      if ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------

!  .. 1 contribution to LS_BEAM, from variations LS_WUPPER

!mick fix 7/20/2016 - turned off for LS calc here

      N  = 1
      DO I = 1, NSTREAMS
        COL2WF_SURF(I,1:NSVARY) = - LS_WUPPER(1:NSVARY,I,N)
      ENDDO

!  Intermediate boundary conditions
!  --------------------------------

!  N1 is the layer below, C0 is the offset
!  .. 2 contributions to LS_BEAM, from variations LS_WUPPER and LS_WLOWER

!mick fix 7/20/2016 - turned off for LS calc here

      DO N = 1, NLAYERS - 1
        N1 = N + 1
        C0 = N*NSTR2 - NSTREAMS
        DO I = 1, NSTR2
          CM = C0 + I
          COL2WF_SURF(CM,1:NSVARY) = LS_WUPPER(1:NSVARY,I,N1) - LS_WLOWER(1:NSVARY,I,N)
        ENDDO
      ENDDO

!  Ground level boundary condition
!  -------------------------------

!  Initialise

      M  = FOURIER
      N  = NLAYERS
      C0 = (N-1)*NSTR2 + NSTREAMS

!mick mod 7/20/2016 - moved these calculations from if blocks below
!                     (applied to both Lambertian & BRDF cases)
!  Compute some intermediate quantities

      DO J = 1, NSTREAMS
        H_WLOWER(J) = QUAD_STRMWGHT(J) * WLOWER(J,N)
        DO AA = 1, NSTREAMS
          H_XPOS(J,AA) = QUAD_STRMWGHT(J) * XPOS(J,AA,N)
          H_XNEG(J,AA) = QUAD_STRMWGHT(J) * XNEG(J,AA,N)
        ENDDO
      ENDDO

!  Diffuse scatter contributions (Lambertian case)

      IF ( .not. DO_BRDF_SURFACE ) then
        IF ( M .eq. 0 ) then

!  Old Code now commented out, 2/2/17 (see below)
!        FRAC1 = R2 * FL1                !@@@@@@@ ROB debug 4/1/11. Wrong
!        FRAC1 = 1.0d0                   !@@@@@@@ ROB debug 4/1/11. Right
!mick note 7/20/2016 - R2_BEAM here = LIDORT's REFL_B  * SURFACE_FACTOR * ALBEDOS_RANKED(Raman_IDX)
!                        (see lrrs_bvproblem.f90 ~line 266)
!mick note 7/20/2016 - R2_HOMP here = LIDORT's REFL_HP * SURFACE_FACTOR * ALBEDOS_RANKED(Raman_IDX) 
!                      R2_HOMM here = LIDORT's REFL_HM * SURFACE_FACTOR * ALBEDOS_RANKED(Raman_IDX)
!                        (see lrrs_bvproblem.f90 ~line 129)
!mick fix 7/20/2016  - replaced this stream loop with new code below
          !DO I = 1, NSTREAMS
          !  CM = C0 + I
          !  HELP = R2_BEAM(I)
          !  DO AA = 1, NSTREAMS
          !    HELP = HELP + LCON(AA,N) * R2_HOMP(I,AA) * KTRANS(AA,N) &
          !                + MCON(AA,N) * R2_HOMM(I,AA)
          !  ENDDO
          !  FACTOR = SURFACE_FACTOR * HELP
          !  COL2WF_SURF(CM,1:NSVARY) = FACTOR
          !ENDDO

!mick fix 7/20/2016 - new code (now more LIDORT-like & similar to BRDF code below)
!  Rob Fix, 2/2/17. Removed Beam contribution, now treated in LS_DIFFUSE_BEAM

          DO I = 1, NSTREAMS
            CM = C0 + I
            HELP   = ZERO       !  HELP = SUM(H_WLOWER(1:NSTREAMS)) ! This is incorrect
            DO AA = 1, NSTREAMS
              REFL_P = SUM(H_XPOS(1:NSTREAMS,AA))
              REFL_M = SUM(H_XNEG(1:NSTREAMS,AA))
              HELP = HELP + LCON(AA,N) * REFL_P * KTRANS(AA,N) &
                          + MCON(AA,N) * REFL_M
            ENDDO
            COL2WF_SURF(CM,1:NSVARY) = SURFACE_FACTOR * HELP
          ENDDO

        ENDIF
      ENDIF

!  Diffuse scatter contributions (BRDF case)
!  Rob Fix, 2/2/17. Removed Beam contribution, now treated in LS_DIFFUSE_BEAM
            !REFL_B = DOT_PRODUCT(H_WLOWER(1:NSTREAMS),LS_BRDF_F(Q,M,1:NSTREAMS,I,BRDF_idx))
!            REFL_B = DOT_PRODUCT(H_WLOWER(1:NSTREAMS),LS_BRDF_F(Q,M,I,1:NSTREAMS,BRDF_idx))
!            HELP = REFL_B

!  mick fix 7/20/2016 - swapped indices (1:NSTREAMS,I) ---> (I,1:NSTREAMS) on LS_BRDF_F array.
!                     1:NSTREAMS is the usual J index where J is incident.
            !REFL_B = DOT_PRODUCT(H_WLOWER(1:NSTREAMS),LS_BRDF_F(Q,M,1:NSTREAMS,I,BRDF_idx))
!            REFL_B = DOT_PRODUCT(H_WLOWER(1:NSTREAMS),LS_BRDF_F(Q,M,I,1:NSTREAMS,BRDF_idx))

      IF ( DO_BRDF_SURFACE ) then
        DO I = 1, NSTREAMS
          CM = C0 + I
          DO Q = 1, NSVARY
            HELP   = ZERO
            DO AA = 1, NSTREAMS
              REFL_P = DOT_PRODUCT(H_XPOS(1:NSTREAMS,AA),LS_BRDF_F(Q,M,I,1:NSTREAMS,BRDF_idx))
              REFL_M = DOT_PRODUCT(H_XNEG(1:NSTREAMS,AA),LS_BRDF_F(Q,M,I,1:NSTREAMS,BRDF_idx))
              HELP = HELP + LCON(AA,N) * REFL_P * KTRANS(AA,N) &
                          + MCON(AA,N) * REFL_M
            ENDDO
            COL2WF_SURF(CM,Q) = SURFACE_FACTOR * HELP
          ENDDO
        ENDDO
      ENDIF

!  Add the contribution from LS_DIFFUSE_BEAM
!   Rob 1/24/17. Appears to be a wrong sign in front of LS_DIFFUSE

      DO I = 1, NSTREAMS
        CM = C0 + I
        COL2WF_SURF(CM,1:NSVARY) = COL2WF_SURF(CM,1:NSVARY) - LS_DIFFUSE_BEAM(1:NSVARY,I)
      ENDDO

!  Add direct beam variation

      DO I = 1, NSTREAMS
        CM = C0 + I
        COL2WF_SURF(CM,1:NSVARY) = COL2WF_SURF(CM,1:NSVARY) + LS_DIRECT_BEAM(1:NSVARY,I)

      ENDDO

!  Return, full case

      IF ( NLAYERS .NE. 1 ) RETURN

!  If Nlayers = 1, special case. Just Copy to output

      DO I = 1, NTOTAL
        SCOL2WF_SURF(I,1:NSVARY) = COL2WF_SURF(I,1:NSVARY)
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSR_BVPCOLUMN_SETUP

!  End Module

      END MODULE lrrs_L_bvproblem_m

