
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
! #   SUBROUTINES :                                             #
! #                                                             #
! #    Homogeneous solutions (elastic and inelastic Eqs.)       #
! #                                                             #
! #             HOMOG_SOLUTION                                  #
! #             UHOMOG_SOLUTION                                 #
! #                                                             #
! #    Classical solutions for elastic solar source field       #
! #                                                             #
! #             BEAM_SOLUTION_MATRIX                            #
! #             BEAM_SOLUTION_ELASTIC                           #
! #             UBEAM_SOLUTION_ELASTIC                          #
! #                                                             #
! #    Green's function solutions for Raman source fields       #
! #           (1 = Whole-layer, 2 = Part-layer)                 #
! #                                                             #
! #             GREENFUNC_SOLUTION_1                            #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are

!    1. Introduction of supplement-created BRDF   expressions in DBEAM_SETUP routines
!    2. Introduction of supplement-created SLEAVE expressions in DBEAM_SETUP routines
!    3. Use of new Taylor series subroutine/module in GREENFUNC_SOLUTION_1
!    4. Bookkeeping improvements (use of "Only", clearer I/O specifications)

!   Items (1) and (2) replace the somewhat ad-hoc SL and BRDF treatments from earlier versions
!     - The new replacements are modeled after the LIDORT Version 3.7 code

!  Programmed, 9/9/15 by R. Spurr. RT SOLUTIONS Inc.
!          Changes marked by the trope "Rob Fix 9/9/15"

!    -- Rob mod 5/12/17 for 2p5a. These 5 setup subroutines moved to new module lrrs_miscsetups.f90
!              CHAPMAN_FUNCTION
!              ATTENUATION_SETUP_1
!              LEGENDRE_SETUP
!              DBEAM_SETUP
!              UDBEAM_SETUP

      MODULE lrrs_rtsolutions_1_m

!      USE LRRS_PARS_m, Only : SDU

      PRIVATE
      PUBLIC :: HOMOG_SOLUTION,        &
                UHOMOG_SOLUTION,       &
                BEAM_SOLUTION_MATRIX,  &
                BEAM_SOLUTION_ELASTIC, &
                UBEAM_SOLUTION_ELASTIC,&
                GREENFUNC_SOLUTION_1

      CONTAINS

      SUBROUTINE HOMOG_SOLUTION &
            ( LAYER, FOURIER, NSTREAMS, NMOMENTS, & ! Input
              OMEGA_PHASMOMS, DELTA_TAU,          & ! Input
              QUAD_STREAMS, QUAD_WEIGHTS,         & ! Input
              PLM_PXI, PLM_MXI,                   & ! Input
              XPOS, XNEG, KEIGEN, KTRANS,         & ! Output
              EIGENMAT_SAVE, EIGENVEC_SAVE,       & ! Output
              DIFVEC_SAVE, DAB_SAVE, SAB_SAVE,    & ! Output
              FAIL, MESSAGE )                       ! Output

!  Numerical solution of Eigenproblem.

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, MAX_STREAMS, MAX_2_STREAMS, MAX_LAYERS, MAX_MOMENTS, &
                              ZERO, ONE, HALF, MAX_TAU_QPATH

!  Module with the Eigenvalue routine

      USE LRRS_AUX2_m, Only : ASYMTX

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  Given layer index and Fourier number

      INTEGER  , INTENT(IN) :: LAYER
      INTEGER  , INTENT(IN) :: FOURIER

!  Number of streams and moments

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: NMOMENTS

!  IOP inputs

      REAL(FPK), INTENT(IN) :: OMEGA_PHASMOMS ( MAX_LAYERS, 0: MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: DELTA_TAU      ( MAX_LAYERS )

!  quadrature inputs

      REAL(FPK), INTENT(IN) :: QUAD_STREAMS ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS ( MAX_STREAMS )

!  Legendre polynomials input

      REAL(FPK), INTENT(IN) :: PLM_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

!  output module arguments
!  -----------------------

!  Eigenvalues + transmittance factors

      REAL(FPK), INTENT(INOUT) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(INOUT) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )

!  Eigenvector solutions

      REAL(FPK), INTENT(INOUT) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(INOUT) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Saved matrices for eigenvalue computation

      REAL(FPK), INTENT(OUT) :: DAB_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: SAB_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: EIGENMAT_SAVE  ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: EIGENVEC_SAVE  ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: DIFVEC_SAVE    ( MAX_STREAMS, MAX_STREAMS )

!  output status

      LOGICAL            , INTENT(OUT) :: FAIL
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE

!  Local variables
!  ---------------

!  local matrix for eigenvalue computation

      REAL(FPK) :: EIGENMAT ( MAX_STREAMS, MAX_STREAMS )

!  (output from Eigenpackage module ASYMTX)

      REAL(FPK) :: KSQ      ( MAX_STREAMS )
!      REAL(FPK) :: WK       ( 2*MAX_2_STREAMS )
      REAL(FPK) :: WK       ( MAX_2_STREAMS )
      REAL(FPK) :: EVEC     ( MAX_STREAMS,MAX_STREAMS )
      INTEGER   :: IER
      LOGICAL   :: ASYMTX_FAILURE

!  Miscellaneous local variables

      INTEGER   :: I, J, I1, L, AA, K
      REAL(FPK) :: DP, DM, SUML, FAC, KVAL
      REAL(FPK) :: NORM, XINV, HELP

!  local variables for error tracing

      CHARACTER (LEN=3) :: CN, CI

!  initialise status
!  -----------------

      FAIL    = .FALSE.
      MESSAGE = ' '

!  Construct Eigenmatrix
!  ---------------------

!  zero the Eigenmatrix

      EIGENMAT = zero
      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          EIGENMAT(I,J) = ZERO
        ENDDO
      ENDDO

!  Develop Sum and Difference matrices

      DO I = 1, NSTREAMS
        XINV = ONE/QUAD_STREAMS(I)
        DO J = 1, NSTREAMS
          FAC = XINV * HALF * QUAD_WEIGHTS(J)
          DP = ZERO
          DM = ZERO
          DO L = FOURIER, NMOMENTS
            DP = DP + PLM_PXI(I,L)*PLM_PXI(J,L)*OMEGA_PHASMOMS(LAYER,L)
            DM = DM + PLM_PXI(I,L)*PLM_MXI(J,L)*OMEGA_PHASMOMS(LAYER,L)
          ENDDO
          SAB_SAVE(I,J) = FAC * ( DP + DM )
          DAB_SAVE(I,J) = FAC * ( DP - DM )
        ENDDO
        SAB_SAVE(I,I) = SAB_SAVE(I,I) - XINV
        DAB_SAVE(I,I) = DAB_SAVE(I,I) - XINV
      ENDDO

!  Compute Eigenmatrix

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          SUML = ZERO
          DO K = 1, NSTREAMS
            SUML = SUML + DAB_SAVE(I,K) * SAB_SAVE(K,J)
          ENDDO
          EIGENMAT(I,J) = SUML
        ENDDO
      ENDDO

!  save Eigenmatrix (original is destroyed by ASMTYX)

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          EIGENMAT_SAVE(I,J) = EIGENMAT(I,J)
        ENDDO
      ENDDO
      
!  Eigensolution package
!  ---------------------

!  Call to ASYMTX

      CALL  ASYMTX &
               ( EIGENMAT, NSTREAMS, MAX_STREAMS, MAX_STREAMS, &
                 EVEC, KSQ, IER, WK, &
                 MESSAGE, ASYMTX_FAILURE )

!  error tracing

      IF ( ASYMTX_FAILURE ) THEN
        WRITE(CN,'(I3)')LAYER
        MESSAGE ='ASYMTX error, Layer='//CN
      ENDIF
      IF ( IER.GT.0 ) THEN
        WRITE(CI,'(I3)')IER
        WRITE(CN,'(I3)')LAYER
        MESSAGE = 'eigenvalue '//CI//'has not converged, Layer='//CN
      ENDIF
      FAIL = ( IER.GT.0 .OR. ASYMTX_FAILURE )
      IF ( FAIL) RETURN

!  Find solution eigenvectors XPOS, XNEG for all eigenvalues
!  ---------------------------------------------------

      DO AA = 1, NSTREAMS

!  Store positive values in output array for each layer

        KVAL = SQRT(KSQ(AA))
        KEIGEN(AA,LAYER) = KVAL

!  Normalize eigenvectors to 1

        NORM = ZERO
        DO I = 1, NSTREAMS
          NORM = NORM + EVEC(I,AA)*EVEC(I,AA)
        ENDDO
        NORM = SQRT(NORM)

!  Find normalized eigenvector EIGENVEC_SAVE

        DO I = 1, NSTREAMS
          EIGENVEC_SAVE(I,AA) = EVEC(I,AA)/NORM
        ENDDO

!  Find difference eigenvector DIFVEC (Siewert's notation)

        DO I = 1, NSTREAMS
          SUML = ZERO
          DO K = 1, NSTREAMS
            SUML = SUML - SAB_SAVE(I,K) * EIGENVEC_SAVE(K,AA)
          ENDDO
          DIFVEC_SAVE(I,AA) = SUML / KVAL
        ENDDO

!  assign original evectors; first N are "DOWN", last N are "UP" (stream

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          XPOS(I,AA,LAYER)  = &
                  HALF * ( EIGENVEC_SAVE(I,AA) + DIFVEC_SAVE(I,AA) )
          XPOS(I1,AA,LAYER) = &
                  HALF * ( EIGENVEC_SAVE(I,AA) - DIFVEC_SAVE(I,AA) )
        ENDDO

!  Use symmetry properties to set -ve eigenvectors

        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          XNEG(I1,AA,LAYER) = XPOS(I,AA,LAYER)
          XNEG(I,AA,LAYER)  = XPOS(I1,AA,LAYER)
        ENDDO

!  End eigenstream loop

      ENDDO

!  Eigenstream transmittance factors for this layer
!     No longer required : KTRANS(AA,LAYER) = SMALLNUM.

      DO AA = 1, NSTREAMS
        HELP = KEIGEN(AA,LAYER)*DELTA_TAU(LAYER)
        IF ( HELP .GT. MAX_TAU_QPATH ) THEN
          KTRANS(AA,LAYER) = ZERO
        ELSE
          KTRANS(AA,LAYER) = DEXP(-HELP)
        ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE HOMOG_SOLUTION

!

      SUBROUTINE UHOMOG_SOLUTION &
            ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS, & ! Input
              OMEGA_PHASMOMS, XPOS, XNEG,                     & ! Input
              PLM_WT_PXI,  PLM_WT_MXI, PLM_PXUI,              & ! Input
              U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )              ! Output

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_USER_STREAMS, MAX_STREAMS, &
                              MAX_2_STREAMS, MAX_LAYERS, MAX_MOMENTS

      IMPLICIT NONE 

!  input arguments
!  ---------------

!  Computational variable control

      INTEGER  , INTENT(IN) :: FOURIER, N
      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: NMOMENTS

!  Number of user streams

      INTEGER  , INTENT(IN) :: N_USER_STREAMS

!  optical properties

      REAL(FPK), INTENT(IN) :: OMEGA_PHASMOMS  ( MAX_LAYERS, 0:MAX_MOMENTS )

!  Quadrature solution vectors

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: XNEG ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Legendre functions

      REAL(FPK), INTENT(IN) :: PLM_WT_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_WT_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_PXUI   ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  output arguments
!  ----------------

!  Homogeneous solutions at the user defined stream angles

      REAL(FPK), INTENT(INOUT) :: U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(INOUT) :: U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Help arrays (for use with linearization)

      REAL(FPK), INTENT(OUT) :: U_HELP_P  ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: U_HELP_M  ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Local variables
!  ---------------

      INTEGER   :: UM, J, J1, L, AA
      REAL(FPK) :: SUM_NEG, SUM_POS, POS1, POS2, NEG1, NEG2, ULP

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------

!  For each eigenvector

      DO AA = 1, NSTREAMS

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors

        DO L = FOURIER, NMOMENTS
          SUM_POS = ZERO
          SUM_NEG = ZERO
          DO  J = 1, NSTREAMS
            J1 = J + NSTREAMS
            POS1 = XPOS(J1,AA,N) * PLM_WT_PXI(J,L)
            POS2 = XPOS(J,AA,N)  * PLM_WT_MXI(J,L)
            NEG1 = XNEG(J1,AA,N) * PLM_WT_PXI(J,L)
            NEG2 = XNEG(J,AA,N)  * PLM_WT_MXI(J,L)
            SUM_POS = SUM_POS + POS1 + POS2
            SUM_NEG = SUM_NEG + NEG1 + NEG2
          ENDDO
          U_HELP_P(AA,L) = SUM_POS
          U_HELP_M(AA,L) = SUM_NEG
        ENDDO

!  Now sum over all harmonic contributions at each user-defined stream

        DO UM = 1, N_USER_STREAMS
          SUM_POS = ZERO
          SUM_NEG = ZERO
          DO L = FOURIER, NMOMENTS
            ULP = PLM_PXUI(UM,L) * OMEGA_PHASMOMS(N,L)
            SUM_POS = SUM_POS + U_HELP_P(AA,L) * ULP
            SUM_NEG = SUM_NEG + U_HELP_M(AA,L) * ULP
          ENDDO
          U_XPOS(UM,AA,N) = SUM_POS
          U_XNEG(UM,AA,N) = SUM_NEG
        ENDDO

!  end eigenvector loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE UHOMOG_SOLUTION

!

      SUBROUTINE BEAM_SOLUTION_MATRIX &
              ( LAYER, NSTREAMS, EIGENMAT_SAVE,  & ! Input
                LOCAL_ASOURCE, LOCAL_PICUTOFF,   & ! Input
                QMAT_SAVE, QPIVOT, FAIL, MESSAGE ) ! Output

!  This is the classical Chandrasekhar beam solution.
!  ( plane parallel or average secant only)
!  Linear Matrix algebra.

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, MAX_STREAMS

!  LAPACK routine

      USE LRRS_AUX2_m, Only : DGETRF

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN) :: LAYER
      INTEGER, INTENT(IN) :: NSTREAMS

!  Beam solution cutoff layer

      INTEGER, INTENT(IN) :: LOCAL_PICUTOFF

!  Beam solution average secant

      REAL(FPK), INTENT(IN) :: LOCAL_ASOURCE

!  Saved matrix from eigenvalue computation

      REAL(FPK), INTENT(IN) :: EIGENMAT_SAVE ( MAX_STREAMS, MAX_STREAMS )

!  output module arguments
!  -----------------------

!  Particular solution: LU-decomposed matrix and Pivot

      REAL(FPK), INTENT(OUT) :: QMAT_SAVE ( MAX_STREAMS, MAX_STREAMS )
      INTEGER  , INTENT(OUT) :: QPIVOT    ( MAX_STREAMS )

!  output status

      LOGICAL            , INTENT(OUT) :: FAIL
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE

!  Local variables
!  ---------------

      INTEGER   :: I, J, INFO
      REAL(FPK) :: INV_COS_SUNZENSQ

!  No particular solution beyond the cutoff layer.
!  -----------------------------------------------

      FAIL    = .FALSE.
      MESSAGE = ' '
      IF ( LAYER .GT. LOCAL_PICUTOFF ) RETURN

!  set local values

      INV_COS_SUNZENSQ = LOCAL_ASOURCE * LOCAL_ASOURCE

!  solution matrix for the reduced problem
!  ( matrix should be saved in the LU decomposition form )

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          QMAT_SAVE(I,J) = EIGENMAT_SAVE(I,J)
        ENDDO
        QMAT_SAVE(I,I) = QMAT_SAVE(I,I) - INV_COS_SUNZENSQ
      ENDDO

!  L-U decomposition of the solution matrix

      CALL DGETRF &
           ( NSTREAMS, NSTREAMS, QMAT_SAVE, MAX_STREAMS, QPIVOT, INFO )

!  error handling

      IF ( INFO .NE. 0 ) THEN
        MESSAGE = 'Beam solution LU decomposition: DGETRF failed'
        FAIL    = .TRUE.
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BEAM_SOLUTION_MATRIX

!

      SUBROUTINE BEAM_SOLUTION_ELASTIC &
            ( LAYER, FOURIER, NMOMENTS, NSTREAMS,    & ! Input
              LOCAL_ASOURCE, LOCAL_INITRANS,         & ! Input
              QMAT_SAVE, QPIVOT, SAB_SAVE, DAB_SAVE, & ! Input
              PLM_00_PXI,  PLM_00_MXI,               & ! Input
              QUAD_STREAMS, OMEGA_PHASMOMS,          & ! Input
              QSUMVEC_SAVE, QDIFVEC_SAVE,            & ! Output
              QVEC_SAVE, QAUX_SAVE,                  & ! Output
              WVECTOR, FAIL, MESSAGE )                 ! Output

!  This is the classical Chandrasekhar beam solution.
!  ( plane parallel or average secant only)
!       Linear Matrix algebra solution

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, HALF, MAX_STREAMS, MAX_2_STREAMS, MAX_LAYERS, MAX_MOMENTS

!  Lapack routine

      USE LRRS_AUX2_m, Only : DGETRS

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  layer index, Fourier number

      INTEGER  , INTENT(IN) :: LAYER
      INTEGER  , INTENT(IN) :: FOURIER

!  Given number of streams and moments

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: NMOMENTS

!  Beam solution inputs
!  --------------------

!  average secant

      REAL(FPK), INTENT(IN) :: LOCAL_ASOURCE

!  transmittances to layer top

      REAL(FPK), INTENT(IN) :: LOCAL_INITRANS

!  Quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STREAMS ( MAX_STREAMS )

!  IOP inputs

      REAL(FPK), INTENT(IN) :: OMEGA_PHASMOMS  ( MAX_LAYERS, 0:MAX_MOMENTS )

!  Saved matrices from eigenvalue computation

      REAL(FPK), INTENT(IN) :: DAB_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: SAB_SAVE ( MAX_STREAMS, MAX_STREAMS )

!  Particular solution: LU-decomposed matrix and Pivot

      REAL(FPK), INTENT(IN) :: QMAT_SAVE ( MAX_STREAMS, MAX_STREAMS )
      INTEGER  , INTENT(IN) :: QPIVOT    ( MAX_STREAMS )

!  legendre polynomial functions

      REAL(FPK), INTENT(IN) :: PLM_00_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_00_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

!  output module arguments
!  -----------------------

!  Particular integral solution  vector

      REAL(FPK), INTENT(INOUT) :: WVECTOR ( MAX_2_STREAMS, MAX_LAYERS )

!  Classical beam solution matrices and vectors
!   Beam solution independent of optical depth (classical solution)
!    [These are required for the linearization]

      REAL(FPK), INTENT(OUT) :: QSUMVEC_SAVE ( MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: QDIFVEC_SAVE ( MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: QVEC_SAVE    ( MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: QAUX_SAVE    ( MAX_STREAMS )

!  output status

      LOGICAL            , INTENT(OUT) :: FAIL
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE

!  Local variables
!  ---------------

      INTEGER   :: I, J, L, INFO, I1, NSTR2
      REAL(FPK) :: TP, TM, EXPFAC, HELP

!  initialise status

      FAIL    = .FALSE.
      MESSAGE = ' '

!  set local values

      EXPFAC = LOCAL_ASOURCE
      NSTR2 = 2 * NSTREAMS

!  No particular solution beyond the cutoff layer.
!  -----------------------------------------------

!  Zero the particular solution and return
!    Zeroing added 18 November 2005 by R. Spurr
!  former code wasn't right...............
!      IF ( LAYER .GT. LOCAL_PICUTOFF ) THEN
!        DO I = 1, 2*NSTREAMS
!          WUPPER(I,LAYER) = ZERO
!          WLOWER(I,LAYER) = ZERO
!            WVECTOR(I,LAYER)= ZERO
!        ENDDO
!          RETURN
!        ENDIF

!  Set the cutoff using the initial transmittance

      IF (LOCAL_INITRANS.EQ.ZERO) THEN
        DO I = 1, NSTR2
          WVECTOR(I,LAYER)= ZERO
        ENDDO
        RETURN
      ENDIF

!  Calculation
!  -----------

!  Set up sum and difference vectors for Beam source terms
!  ( sum vector may be required again in linearization )

      DO I = 1, NSTREAMS
        TP = ZERO
        TM = ZERO
        DO L = FOURIER, NMOMENTS
          TP = TP + OMEGA_PHASMOMS(LAYER,L) * PLM_00_PXI(I,L)
          TM = TM + OMEGA_PHASMOMS(LAYER,L) * PLM_00_MXI(I,L)
        ENDDO
        QSUMVEC_SAVE(I) =  ( TP + TM ) / QUAD_STREAMS(I)
        QDIFVEC_SAVE(I) =  ( TP - TM ) / QUAD_STREAMS(I)
      ENDDO

!  RHS vector for the reduced problem
!  ( this vector will be the answer after the linear algebra solution,
!    and may be needed again if there is linearization )

      DO I = 1, NSTREAMS
        HELP = ZERO
        DO J = 1, NSTREAMS
          HELP = HELP - DAB_SAVE(I,J)*QSUMVEC_SAVE(J)
        ENDDO
        QVEC_SAVE(I) = HELP + QDIFVEC_SAVE(I) * EXPFAC
      ENDDO

!  Solution of reduced problem by back-substitution

      CALL DGETRS &
          ('N',NSTREAMS,1,QMAT_SAVE,MAX_STREAMS,QPIVOT, &
           QVEC_SAVE,MAX_STREAMS,INFO)

      IF ( INFO .NE. 0 ) THEN
        MESSAGE='Classical Beam solution back-substitution (DGETRS)'
        FAIL   = .TRUE.
        RETURN
      ENDIF

!  Assigning beam particular integral solution vector

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        HELP = ZERO
        DO J = 1, NSTREAMS
          HELP = HELP - SAB_SAVE(I,J)*QVEC_SAVE(J)
        ENDDO
        QAUX_SAVE(I) = ( HELP - QSUMVEC_SAVE(I) ) / EXPFAC
        WVECTOR(I,LAYER)  = HALF * ( QVEC_SAVE(I) + QAUX_SAVE(I) )
        WVECTOR(I1,LAYER) = HALF * ( QVEC_SAVE(I) - QAUX_SAVE(I) )
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE BEAM_SOLUTION_ELASTIC

!

      SUBROUTINE UBEAM_SOLUTION_ELASTIC &
                ( DO_UPWELLING, DO_DNWELLING, FOURIER,       & ! Inputs
                  N, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Inputs
                  NSTREAMS, NMOMENTS, N_USER_STREAMS,        & ! Inputs
                  PICUTOFF, WPARTIC, OMEGAMOMS,              & ! Inputs
                  PLM_WT_PXI, PLM_WT_MXI,                    & ! Inputs
                  PLM_00_PXUI, PLM_PXUI,                     & ! Inputs
                  PLM_00_MXUI, PLM_MXUI,                     & ! Inputs
                  W_HELP, U_WPOS1, U_WPOS2, U_WNEG1, U_WNEG2 ) ! Outputs

!  Particular integral solutions at user angles

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_USER_STREAMS, &
                              MAX_STREAMS, MAX_2_STREAMS, MAX_MOMENTS

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Directional control

      LOGICAL, INTENT(IN) :: DO_UPWELLING
      LOGICAL, INTENT(IN) :: DO_DNWELLING

!  Computational control

      INTEGER  , INTENT(IN) :: N
      INTEGER  , INTENT(IN) :: FOURIER
      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: NMOMENTS

!  Layer existence

      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP ( MAX_LAYERS )
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN ( MAX_LAYERS )

!  number of user-defined stream angles

      INTEGER  , INTENT(IN) :: N_USER_STREAMS

!  Beam solution cutoff layer

      INTEGER  , INTENT(IN) :: PICUTOFF

!  Beam solution vector

      REAL(FPK), INTENT(IN) :: WPARTIC   ( MAX_2_STREAMS, MAX_LAYERS )

!  Initial transmittance of beam to top of layer. NOT required now.....................................
!      REAL(FPK), INTENT(IN) :: INITRANS  ( MAX_LAYERS )

!  IOP input

      REAL(FPK), INTENT(IN) :: OMEGAMOMS ( MAX_LAYERS, 0:MAX_MOMENTS )

!  Legendre functions

      REAL(FPK), INTENT(IN) :: PLM_WT_PXI  ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_WT_MXI  ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_00_PXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_PXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_00_MXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_MXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  subroutine output arguments
!  ---------------------------

!  User-defined stream angles beam solution vectors

      REAL(FPK), INTENT(INOUT) :: &
            U_WPOS1 ( MAX_USER_STREAMS, MAX_LAYERS ), &
            U_WPOS2 ( MAX_USER_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(INOUT) :: &
            U_WNEG1 ( MAX_USER_STREAMS, MAX_LAYERS ), &
            U_WNEG2 ( MAX_USER_STREAMS, MAX_LAYERS )

!  Help array (for linearization only)

      REAL(FPK), INTENT(OUT) :: W_HELP(0:MAX_MOMENTS)

!  Local variables
!  ---------------

      INTEGER   :: UM, J, J1, L
      REAL(FPK) :: SUML, POS1, POS2, NEG1, HELP(0:MAX_MOMENTS)

!  No particular solution beyond the cutoff layer, Exit if so.

      IF ( N .GT. PICUTOFF ) THEN
        IF ( DO_UPWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            U_WPOS1(UM,N) = ZERO
            U_WPOS2(UM,N) = ZERO
          ENDDO
        ENDIF
        IF ( DO_DNWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            U_WNEG1(UM,N) = ZERO
            U_WNEG2(UM,N) = ZERO
          ENDDO
        ENDIF
        RETURN
      ENDIF

!  For each moment & particulate, do inner sum over computational angles
!  (required for both directions)

      DO L = FOURIER, NMOMENTS
        SUML = ZERO
        DO  J = 1, NSTREAMS
          J1 = J + NSTREAMS
          POS1 = WPARTIC(J1,N) * PLM_WT_PXI(J,L)
          POS2 = WPARTIC(J,N)  * PLM_WT_MXI(J,L)
          SUML = SUML + POS1 + POS2
        ENDDO
        W_HELP(L) = SUML
        HELP(L)   = SUML * OMEGAMOMS(N,L)
      ENDDO

!  Upwelling, sum over all harmonic contributions, each user-defined str

      IF ( DO_UPWELLING ) THEN
        IF ( STERM_LAYERMASK_UP(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            POS1 = ZERO
            DO L = FOURIER, NMOMENTS
              POS1 = POS1 + OMEGAMOMS(N,L) * PLM_00_PXUI(UM,L)
            ENDDO
!            U_WPOS1(UM,N) = POS1 * INITRANS(N)
            U_WPOS1(UM,N) = POS1
            POS1 = ZERO
            DO L = FOURIER, NMOMENTS
              POS1 = POS1 + HELP(L) * PLM_PXUI(UM,L)
            ENDDO
!            U_WPOS2(UM,N) = POS1 * INITRANS(N)
            U_WPOS2(UM,N) = POS1
          ENDDO
        ENDIF
      ENDIF

!  Downwelling, sum over all harmonic contributions, each user-defined s

      IF ( DO_DNWELLING ) THEN
        IF ( STERM_LAYERMASK_DN(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            NEG1 = ZERO
            DO L = FOURIER, NMOMENTS
              NEG1 = NEG1 + OMEGAMOMS(N,L) * PLM_00_MXUI(UM,L)
            ENDDO
!            U_WNEG1(UM,N) = NEG1 * INITRANS(N)
            U_WNEG1(UM,N) = NEG1
            NEG1 = ZERO
            DO L = FOURIER, NMOMENTS
              NEG1 = NEG1 + HELP(L) * PLM_MXUI(UM,L)
            ENDDO
!            U_WNEG2(UM,N) = NEG1 * INITRANS(N)
            U_WNEG2(UM,N) = NEG1
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE UBEAM_SOLUTION_ELASTIC

!

      SUBROUTINE GREENFUNC_SOLUTION_1 &
         ( FOURIER, LAYER, SOLUTION, DO_MULTIPLIER,               & ! Inputs
           TAYLOR_ORDER, TAYLOR_SMALL,                            & ! Inputs
           NSTREAMS, QUAD_WEIGHTS, DELTAS,                        & ! Inputs
           XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,                  & ! Inputs
           PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS,        & ! Inputs
           QSOURCE_QPOS, QSOURCE_QNEG,                            & ! Inputs
           WUPPER, WLOWER, ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P ) ! Output

!  This is the Green's function solution, one layer one contribution

!  LRRS Version 2.5. Rob Fix 9/10/15.
!    - Added Taylor_order variable to argument list
!    - Use Taylor series subroutine and module, replaces older code.

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m  , Only : FPK, ZERO, ONE, PI4, UPIDX, DNIDX, &
                                MAX_STREAMS, MAX_2_STREAMS, MAX_LAYERS

!  Taylor series module

      USE lrrs_Taylor_m, Only : Taylor_series_1

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  Fourier Index

      INTEGER  , INTENT(IN) :: FOURIER

!  Given layer, solution index

      INTEGER  , INTENT(IN) :: LAYER
      INTEGER  , INTENT(IN) :: SOLUTION

!  Control for calculating Green's function multiplier


      LOGICAL  , INTENT(IN) :: DO_MULTIPLIER

!  Small number limit, Taylor order (Taylor series expansions)

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  quadrature inputs

      INTEGER  , INTENT(IN) :: NSTREAMS
      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS ( MAX_STREAMS )

!  source vector inputs
!  --------------------

!  layer cutoff

      INTEGER, INTENT(IN) ::   PICUTOFF

!  Flipper for using -x instead of x as inteegration argument

      LOGICAL, INTENT(IN) ::   FLIPPER

!  Average secant, initial transmittance

      REAL(FPK), INTENT(IN) :: ASOURCE
      REAL(FPK), INTENT(IN) :: INITRANS

!  layer transmittance

      REAL(FPK), INTENT(IN) :: DELTRANS

!  optical thickness (only required for Taylor series expansions)

      REAL(FPK), INTENT(IN) :: DELTAS

!  Source vectors in quadrature directions

      REAL(FPK), INTENT(IN) :: QSOURCE_QPOS  ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QSOURCE_QNEG  ( MAX_STREAMS )

!  homogeneous solution inputs
!  ---------------------------

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Eigenvalues and whole-layer transmittance factors

      REAL(FPK), INTENT(IN) :: KTRANS    ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KEIGEN    ( MAX_STREAMS, MAX_LAYERS )

!  Norm values

      REAL(FPK), INTENT(IN) :: EIGENNORM_SAVE ( MAX_STREAMS, MAX_LAYERS)

!  output module arguments
!  -----------------------

!  Total RRS-summed Particular integral at layer boundaries

      REAL(FPK), INTENT(INOUT) :: WUPPER  ( MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(INOUT) :: WLOWER  ( MAX_2_STREAMS, MAX_LAYERS )

!  Saved ATERM and BTERM  vectors

      REAL(FPK), INTENT(OUT) :: ATERM_SAVE ( MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: BTERM_SAVE ( MAX_STREAMS )

!  output MULTIPLIERS

      REAL(FPK), INTENT(INOUT) :: MULT_M ( MAX_STREAMS )
      REAL(FPK), INTENT(INOUT) :: MULT_P ( MAX_STREAMS )

!  local variables
!  ===============

      INTEGER ::   AA, I, I1, N, M
      REAL(FPK) :: TPA, TMA, SUM_LA, SUM_LB !,F1
      REAL(FPK) :: S_P_U, S_P_L, S_M_U, S_M_L
      REAL(FPK) :: RHO_M, RHO_P, ZDEL, ZWDEL

!   RRS Particular integral at layer boundaries (one term)

      REAL(FPK) :: WLEVELS ( MAX_2_STREAMS, 2 )

!  multipliers

      REAL(FPK) :: CFUNC    ( MAX_STREAMS )
      REAL(FPK) :: DFUNC    ( MAX_STREAMS )
      REAL(FPK) :: GFUNC_UP ( MAX_STREAMS )
      REAL(FPK) :: GFUNC_DN ( MAX_STREAMS )

!  Zero the solution (initialise boundary values for first solution)

      IF ( SOLUTION .EQ. 1 ) THEN
        DO I = 1, 2 * NSTREAMS
          WUPPER(I,LAYER) = ZERO
          WLOWER(I,LAYER) = ZERO
        ENDDO
      ENDIF

!  Layer

      N = LAYER

!  Fourier (debug only)

      M  = FOURIER

!  No particular solution beyond the cutoff layer.
!    ( Zero the boundary layer values and exit )
!mick fix 10/11/2016 - define ATERM_SAVE & BTERM_SAVE for this case

      IF ( N .GT. PICUTOFF ) THEN
        DO AA = 1, NSTREAMS
          ATERM_SAVE(AA) = ZERO
          BTERM_SAVE(AA) = ZERO
        ENDDO
        DO I = 1, 2 * NSTREAMS
          WLEVELS(I,UPIDX) = ZERO
          WLEVELS(I,DNIDX) = ZERO
        ENDDO
        GO TO 4567
      ENDIF

!  Optical depth integrations for the discrete ordinate solution
!  =============================================================

!  Avoid this multiplier calculation if already done

!   ---> Allows for degeneracy (old code, this was treated separately)
!   ---> Allows for the Flipper formalism (now fully worked out)
!   ---> Note the small numbers expansion (from Fluoresence model)

      IF ( DO_MULTIPLIER ) THEN
        DO AA = 1, NSTREAMS
          ZDEL  = KTRANS(AA,N)
          ZWDEL = ZDEL * DELTRANS
          RHO_M = ASOURCE - KEIGEN(AA,N)
          RHO_P = ASOURCE + KEIGEN(AA,N)

!  Set the downwelling multiplier
!    Rob Fix 9/10/15. Use Taylor series expansion, instead of LIMIT_GCFUNC
!      * Old Call  -->   CALL LIMIT_GCFUNC ( RHO_M, DELTAS, ZDEL, MULT_M(AA) )
!      * Debugged again. 1/15/17

          IF ( ABS(RHO_M).LT.TAYLOR_SMALL ) THEN
            CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, RHO_M, DELTAS, DELTRANS, ONE, MULT_M(AA) )
          ELSE
            MULT_M(AA) = ( ZDEL - DELTRANS ) / RHO_M
          ENDIF
          MULT_P(AA) = ( ONE - ZWDEL ) / RHO_P

!  End stream cases and multiplier existence

        ENDDO
      ENDIF

!  Distinguish between the flipper cases
!   (multipliers are reversed if -x rather than +x is tau-variable)

      IF ( .NOT.FLIPPER ) THEN
        DO AA = 1, NSTREAMS
          CFUNC(AA)  = MULT_M(AA) * INITRANS
          DFUNC(AA)  = MULT_P(AA) * INITRANS
        ENDDO
      ELSE
        DO AA = 1, NSTREAMS
          CFUNC(AA)  = MULT_P(AA) * INITRANS
          DFUNC(AA)  = MULT_M(AA) * INITRANS
        ENDDO
      ENDIF

!  For each eigenstream, get the terms ATERM_SAVE and BTERM_SAVE
!  =============================================================

      DO AA = 1, NSTREAMS
        SUM_LA = ZERO
        SUM_LB = ZERO
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          TPA = QSOURCE_QPOS(I) * XPOS(I,AA,N)  + &
                QSOURCE_QNEG(I) * XPOS(I1,AA,N)
          TMA = QSOURCE_QNEG(I) * XPOS(I,AA,N)  + &
                QSOURCE_QPOS(I) * XPOS(I1,AA,N)
          SUM_LA  = SUM_LA + TPA * QUAD_WEIGHTS(I)
          SUM_LB  = SUM_LB + TMA * QUAD_WEIGHTS(I)
        ENDDO
        ATERM_SAVE(AA) = SUM_LA / EIGENNORM_SAVE(AA,N)
        BTERM_SAVE(AA) = SUM_LB / EIGENNORM_SAVE(AA,N)
      ENDDO

!  Local Green function multipliers
!  ================================

!  For each eigenstream
!    (Initrans multiplication taken care of in source function definition above)

      DO AA = 1, NSTREAMS
        GFUNC_DN(AA) = CFUNC(AA) * ATERM_SAVE(AA)
        GFUNC_UP(AA) = DFUNC(AA) * BTERM_SAVE(AA)
      ENDDO

!  Set particular integral from Green function expansion
!  =====================================================

!  particular integral WLEVELS at lower and upper boundaries

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        S_P_U = ZERO
        S_P_L = ZERO
        S_M_U = ZERO
        S_M_L = ZERO
        DO AA = 1, NSTREAMS
          S_P_U = S_P_U + GFUNC_UP(AA)*XPOS(I1,AA,N)
          S_M_U = S_M_U + GFUNC_UP(AA)*XPOS(I,AA,N)
          S_P_L = S_P_L + GFUNC_DN(AA)*XPOS(I,AA,N)
          S_M_L = S_M_L + GFUNC_DN(AA)*XPOS(I1,AA,N)
        ENDDO
        WLEVELS(I,UPIDX)  = S_P_U
        WLEVELS(I1,UPIDX) = S_M_U
        WLEVELS(I1,DNIDX) = S_M_L
        WLEVELS(I,DNIDX)  = S_P_L
      ENDDO

!  Continuation point for avoiding this calculation

 4567 CONTINUE

!  Add single RRS particular integral to the RRS-summed total P.I.


      DO I = 1, 2 * NSTREAMS
        WUPPER(I,N) = WUPPER(I,N) + WLEVELS(I,UPIDX)
        WLOWER(I,N) = WLOWER(I,N) + WLEVELS(I,DNIDX)
      ENDDO

!  Debug
!  Note the Debug Write -> Very potent, gets all sources
!      i   = 4 ; deb = .true.
!      do ss = 1, nss
!        if ( solution.eq.solmask(ss)) deb = .false.
!      enddo
!  profiles
!      kfd = 3
!      if ((n.ge.kfd).or.((n.lt.kfd).and.deb)) &
!               write(56,'(3i5,1p2e24.14)')solution,n,i,WLEVELS(I,UPIDX),WLEVELS(I,DNIDX)
!  column
!      kfd = 0
!      write(56,'(3i5,1p2e24.14)')solution,n,i,WLEVELS(I,UPIDX),WLEVELS(I,DNIDX)

!  Finish

      RETURN
      END SUBROUTINE GREENFUNC_SOLUTION_1

!  Finish Module

      END MODULE lrrs_rtsolutions_1_m
