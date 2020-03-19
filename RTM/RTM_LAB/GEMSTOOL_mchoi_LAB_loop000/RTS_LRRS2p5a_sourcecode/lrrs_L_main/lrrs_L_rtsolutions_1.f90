
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
! #   LINEARIZATION SUBROUTINES :                               #
! #                                                             #
! #     LC  denotes column   linearization                      #
! #     LCH denotes column   linearization with heights         #
! #     LP  denotes profile  linearization                      #
! #     LS  denotes surface  linearization                      #
! #     LPC denotes LC or LP linearization                      #
! #                                                             #
! #    Homogeneous solutions (elastic and inelastic Eqs.)       #
! #                                                             #
! #             LPC_HOMOG_SOLUTION                              #
! #             LPC_UHOMOG_SOLUTION                             #
! #                                                             #
! #    Classical solutions for elastic solar source field       #
! #                                                             #
! #             LPC_BEAM_SOLUTION_ELASTIC                       #
! #             LPC_UBEAM_SOLUTION_ELASTIC                      #
! #                                                             #
! #    Green's function solutions for Raman source fields       #
! #           (1 = Whole-layer, 2 = Part-layer)                 #
! #                                                             #
! #             LPC_GREENFUNC_SOLUTION_1                        #
! #             LS_GREENFUNC_SOLUTION_1                         #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. Main changes to this module (from V2.3) are

!    1. Introduction of supplement-created BRDF   expressions in DBEAM_SETUP routines
!    2. Introduction of supplement-created SLEAVE expressions in DBEAM_SETUP routines
!    3. Use of new Taylor series subroutine/module in LPC_GREENFUNC_SOLUTION_1
!    4. Bookkeeping improvements (use of "Only", clearer I/O specifications)

!   Items (1) and (2) replace the somewhat ad-hoc SL and BRDF treatments from earlier versions
!     - The new replacements are modeled after the LIDORT Version 3.7 code

!  Regular code was programmed 9/9/15 by R. Spurr. RT SOLUTIONS Inc.
!    Changes there marked by the trope "Rob Fix 9/9/15"

!  Linearized code was programmed 10/10/15 by R. Spurr. RT SOLUTIONS Inc.
!    Changes there marked by the trope "Rob Fix 10/10/15". THIS MODULE.

!    -- Rob mod 5/12/17 for 2p5a. These 4 setup subroutines moved to new module lrrs_L_miscsetups.f90
!              CHAPMAN_FUNCTION_PLUS
!              LC_ATTENUATION_SETUP_1
!              LP_ATTENUATION_SETUP_1
!              LCH_ATTENUATION_SETUP_1

      MODULE lrrs_L_rtsolutions_1_m

      USE LRRS_PARS_m, Only : LDU

      PRIVATE
      PUBLIC :: LPC_HOMOG_SOLUTION,         &
                LPC_UHOMOG_SOLUTION,        &
                LPC_BEAM_SOLUTION_ELASTIC,  &
                LPC_UBEAM_SOLUTION_ELASTIC, &
                LPC_GREENFUNC_SOLUTION_1,   &
                LS_GREENFUNC_SOLUTION_1

      CONTAINS

      SUBROUTINE LPC_HOMOG_SOLUTION &
            ( LAYER, FOURIER, NSTREAMS, NMOMENTS, DO_VARY, NPARS, & ! Input
              L_OMEGA_PHASMOMS, DELTA_TAU, L_DELTA_TAU,           & ! Input
              QUAD_STREAMS, QUAD_WEIGHTS, PLM_PXI, PLM_MXI,       & ! Input
              KEIGEN, KTRANS, EIGENMAT_SAVE,                      & ! Input
              EIGENVEC_SAVE, DIFVEC_SAVE, DAB_SAVE, SAB_SAVE,     & ! Input
              L_EIGENMAT, L_DAB_SAVE, L_SAB_SAVE,                 & ! Output
              L_XPOS, L_XNEG, L_KEIGEN, L_KTRANS,                 & ! Output
              FAIL, MESSAGE )                                       ! Output

!  This is the linearization of the Eigenproblem:
!     Generates linearized Eigenvalues          L_KEIGEN
!                          Eigensolutions       L_XPOS/L_XNEG
!                          Eigen transmittances L_KTRANS

!     Other outputs are needed for the Beam solution linearizations:
!        L_EIGENMAT,    L_DAB_SAVE, L_SAB_SAVE

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, MAX_STREAMS, MAX_2_STREAMS, MAX_STREAMS_P1, MAX_LAYERS, MAX_MOMENTS, &
                              MAX_ATMOSWFS, ZERO, ONE, HALF, TWO, MAX_TAU_QPATH

!  Linear algebra routines

      USE LRRS_AUX2_m, Only : DGETRF, DGETRS

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  Given layer index and Fourier number

      INTEGER  , INTENT(IN) :: LAYER
      INTEGER  , INTENT(IN) :: FOURIER

!  Number of streams and moments

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: NMOMENTS

!  Linearization control

      LOGICAL  , INTENT(IN) :: DO_VARY
      INTEGER  , INTENT(IN) :: NPARS

!  Linearized IOP inputs

      REAL(FPK), INTENT(IN) :: L_OMEGA_PHASMOMS ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: L_DELTA_TAU      ( MAX_ATMOSWFS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: DELTA_TAU        ( MAX_LAYERS )

!  quadrature inputs

      REAL(FPK), INTENT(IN) :: QUAD_STREAMS ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS ( MAX_STREAMS )

!  Legendre polynomials input

      REAL(FPK), INTENT(IN) :: PLM_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Eigenvalues + transmittance factors

      REAL(FPK), INTENT(IN) :: KTRANS ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KEIGEN ( MAX_STREAMS, MAX_LAYERS )

!  Saved matrices from eigenvalue computation

      REAL(FPK), INTENT(IN) :: DAB_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: SAB_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: EIGENMAT_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: EIGENVEC_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: DIFVEC_SAVE   ( MAX_STREAMS, MAX_STREAMS )

!  output module arguments
!  -----------------------

!  Linearized Eigenvalues + transmittance factors
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: L_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(INOUT) :: L_KEIGEN ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

!  Linearized Eigenvector solutions
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: L_XPOS &
              ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(INOUT) :: L_XNEG &
              ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Linearized eigenmatrices

      REAL(FPK), INTENT(OUT) :: L_DAB_SAVE ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: L_SAB_SAVE ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: L_EIGENMAT ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )

!  output status

      LOGICAL            , INTENT(OUT) :: FAIL
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE

!  Local variables
!  ---------------

!  local matrices for computation

      REAL(FPK) :: &
            HMAT(MAX_STREAMS_P1,MAX_STREAMS_P1), &
            HVEC(MAX_STREAMS_P1,MAX_ATMOSWFS), &
            L_DIFVEC(MAX_STREAMS), L_EIGENVEC(MAX_STREAMS)

!  extra output/input for linear-algebra solver (LAPACK)

      INTEGER   :: IPIV(MAX_STREAMS_P1), INFO

!  Miscellaneous local variables

      INTEGER   :: I, J, I1, L, AA, K, N, Q, NSTRM_P1
      REAL(FPK) :: FAC, KVAL, DP, DM, SUMR, L_KVAL, KSQD
      REAL(FPK) :: XINV, HELP, LTQ

!  initialise status
!  -----------------

      FAIL    = .FALSE.
      MESSAGE = ' '

!  Layer

      N = LAYER

!  Zero entries if layer not varying, and return

      IF ( .NOT. DO_VARY ) THEN
        DO Q = 1, NPARS
          DO AA = 1, NSTREAMS
            L_KEIGEN(Q,AA,N) = ZERO ; L_KTRANS(Q,AA,N) = ZERO
            DO I = 1, 2*NSTREAMS
              L_XPOS(Q,I,AA,N) = ZERO ; L_XNEG(Q,I,AA,N)  = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Construct Linearized Eigenmatrix
!  --------------------------------

!  Develop Linearized Sum and Difference matrices

      DO I = 1, NSTREAMS
        XINV = ONE/QUAD_STREAMS(I)
        DO J = 1, NSTREAMS
          FAC = XINV * HALF * QUAD_WEIGHTS(J)
          DO Q = 1, NPARS
           DP = ZERO
           DM = ZERO
           DO L = FOURIER, NMOMENTS
            DP = DP + PLM_PXI(I,L)*PLM_PXI(J,L)*L_OMEGA_PHASMOMS(Q,N,L)
            DM = DM + PLM_PXI(I,L)*PLM_MXI(J,L)*L_OMEGA_PHASMOMS(Q,N,L)
           ENDDO
           L_SAB_SAVE(Q,I,J) = FAC * ( DP + DM )
           L_DAB_SAVE(Q,I,J) = FAC * ( DP - DM )
          ENDDO
        ENDDO
      ENDDO

!  Compute linearized Eigenmatrix

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO Q = 1, NPARS
           SUMR = ZERO
           DO K = 1, NSTREAMS
            SUMR = SUMR + L_DAB_SAVE(Q,I,K) *   SAB_SAVE(K,J) &
                      +   DAB_SAVE(I,K)   * L_SAB_SAVE(Q,K,J)
           ENDDO
           L_EIGENMAT(Q,I,J) = SUMR
          ENDDO
        ENDDO
      ENDDO

!  Linearization matrix and vector setup
!  -------------------------------------

!  Start eigenstream loop

      DO AA = 1, NSTREAMS

        KVAL = KEIGEN(AA,N)
        KSQD = KVAL * KVAL

!  initialise solution matrix HMAT (important to do this!)

        DO I = 1, MAX_STREAMS_P1
          DO J = 1, MAX_STREAMS_P1
            HMAT(I,J) = ZERO
          ENDDO
        ENDDO

!  Determine solution matrix HMAT (A in the matrix equation A.x = B)

        NSTRM_P1 = NSTREAMS + 1
        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            HMAT(I,J+1) = - EIGENMAT_SAVE(I,J)
          ENDDO
          HMAT(I,I+1) = HMAT(I,I+1) + KSQD
          HMAT(I,1)   = TWO * KEIGEN(AA,N) * EIGENVEC_SAVE(I,AA)
          HMAT(NSTRM_P1,I+1) = EIGENVEC_SAVE(I,AA)
        ENDDO
        HMAT(NSTRM_P1,1) = ZERO

!  solution column vectors (B in the matrix equation A.x = B).
!    (One for each parameter to be varied)

        DO Q = 1, NPARS
          DO I = 1, NSTREAMS
            HELP = ZERO
            DO J = 1, NSTREAMS
              HELP = HELP + L_EIGENMAT(Q,I,J) * EIGENVEC_SAVE(J,AA)
            ENDDO
            HVEC(I,Q) = HELP
          ENDDO
          HVEC(NSTRM_P1,Q) = ZERO
        ENDDO

!  Solve matrix system for linearization factors
!  ---------------------------------------------

!   .. LU-decomposition of the matrix HMAT using DGETRF

        CALL DGETRF (NSTRM_P1,NSTRM_P1,HMAT,MAX_STREAMS_P1,IPIV,INFO)
        IF ( INFO .NE. 0 ) THEN
          Message =  ' Homogeneous Linearization DGETRF'
          FAIL = .TRUE. ; RETURN
        ENDIF

!   .. Back substitution for column vectors (one for each parameter vary

        CALL DGETRS ('N',NSTRM_P1,NPARS,HMAT, &
                         MAX_STREAMS_P1,IPIV,HVEC,MAX_STREAMS_P1,INFO)

        IF ( INFO .NE. 0 ) THEN
          Message =  ' Homogeneous Linearization DGETRS'
          FAIL = .TRUE. ; RETURN
        ENDIF

!  Assign linearized eigensolution vectors
!  ---------------------------------------

!  Start loop over varying parameters

        DO Q = 1, NPARS

!  assign linearization for eigenvalue K

          L_KVAL = HVEC(1,Q)
          L_KEIGEN(Q,AA,N) = L_KVAL

!  linearization of the actual eigenvector

          DO I = 1, NSTREAMS
            L_EIGENVEC(I) = HVEC(I+1,Q)
          ENDDO

!  linearized difference vector

          DO I = 1, NSTREAMS
            SUMR = ZERO
            DO K = 1, NSTREAMS
              SUMR = SUMR - L_SAB_SAVE(Q,I,K) *   EIGENVEC_SAVE(K,AA) &
                          -   SAB_SAVE(I,K)   * L_EIGENVEC(K)
            ENDDO
            HELP = ( SUMR - L_KVAL * DIFVEC_SAVE(I,AA) ) / KVAL
            L_DIFVEC(I) = HELP
          ENDDO

!  assign linearization for 'positive' homogeneous solution vectors

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            L_XPOS(Q,I1,AA,N) = HALF * ( L_EIGENVEC(I) - L_DIFVEC(I) )
            L_XPOS(Q,I,AA,N)  = HALF * ( L_EIGENVEC(I) + L_DIFVEC(I) )
          ENDDO

!  symmetry for linearized 'negative' homogeneous solution vectors

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            L_XNEG(Q,I1,AA,N) = L_XPOS(Q,I,AA,N)
            L_XNEG(Q,I,AA,N)  = L_XPOS(Q,I1,AA,N)
          ENDDO

!  End parameter loop

        ENDDO

!  End eigenstream loop

      ENDDO

!  Linearized Eigentransmittances
!  ------------------------------

      DO AA = 1, NSTREAMS
        DO Q = 1, NPARS
          LTQ = L_KEIGEN(Q,AA,N) *   DELTA_TAU(N) &
                + KEIGEN(AA,N)   * L_DELTA_TAU(Q,N)
          L_KTRANS(Q,AA,N) = - LTQ * KTRANS(AA,N)
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LPC_HOMOG_SOLUTION

!

      SUBROUTINE LPC_UHOMOG_SOLUTION &
            ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS,   & ! Inputs
              DO_VARY, NPARS, OMEGA_PHASMOMS, L_OMEGA_PHASMOMS, & ! Inputs
              L_XPOS, L_XNEG, PLM_WT_PXI,  PLM_WT_MXI,          & ! Inputs
              PLM_PXUI, U_HELP_P, U_HELP_M,                     & ! Inputs
              L_U_XPOS, L_U_XNEG )                                ! Outputs

!  Module of dimensions and numbers
!  --------------------------------

!  Module of dimensions and numbers

      USE LRRS_PARS_m, Only : FPK, MAX_STREAMS, MAX_2_STREAMS, MAX_USER_STREAMS, MAX_LAYERS,  &
                              MAX_MOMENTS, MAX_ATMOSWFS, ZERO, ONE, HALF, MAX_TAU_QPATH

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Linearization control

      LOGICAL  , INTENT(IN) :: DO_VARY
      INTEGER  , INTENT(IN) :: NPARS

!  Computational variables

      INTEGER  , INTENT(IN) :: FOURIER, N
      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: NMOMENTS

!  Number of user streams

      INTEGER  , INTENT(IN) :: N_USER_STREAMS

!  optical properties and linearizations

      REAL(FPK), INTENT(IN) :: OMEGA_PHASMOMS   ( MAX_LAYERS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: L_OMEGA_PHASMOMS ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS )

!  linearization of Quadrature solution vectors

      REAL(FPK), INTENT(IN) :: L_XPOS ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_XNEG ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Legendre functions

      REAL(FPK), INTENT(IN) :: PLM_WT_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_WT_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_PXUI   ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  Help matrices

      REAL(FPK), INTENT(IN) :: U_HELP_P  ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: U_HELP_M  ( MAX_STREAMS, 0:MAX_MOMENTS )

!  output arguments
!  ----------------

!  Linearized Homogeneous solutions at the user defined stream angles
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: L_U_XPOS &
            ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(INOUT) :: L_U_XNEG &
            ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Local variables
!  ---------------

      INTEGER   :: UM, J, J1, L, AA, Q
      REAL(FPK) :: SUM_NEG, SUM_POS
      REAL(FPK) :: POS1, POS2, NEG1, NEG2, ULP, L_ULP
      REAL(FPK) :: L_U_HELP_P  ( MAX_ATMOSWFS, 0:MAX_MOMENTS )
      REAL(FPK) :: L_U_HELP_M  ( MAX_ATMOSWFS, 0:MAX_MOMENTS )

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

!  If there is no variation, zero solutions and exit

      IF ( .NOT. DO_VARY ) THEN
        DO AA = 1, NSTREAMS
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, NPARS
              L_U_XPOS(Q,UM,AA,N) = ZERO
              L_U_XNEG(Q,UM,AA,N) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For each eigenvector

      DO AA = 1, NSTREAMS

!  For each moment, do inner sum over computational angles
!  for the positive and negative linearized eigenvectors

        DO Q = 1, NPARS
          DO L = FOURIER, NMOMENTS
            SUM_POS = ZERO
            SUM_NEG = ZERO
            DO  J = 1, NSTREAMS
              J1 = J + NSTREAMS
              POS1 = L_XPOS(Q,J1,AA,N) * PLM_WT_PXI(J,L)
              POS2 = L_XPOS(Q,J,AA,N)  * PLM_WT_MXI(J,L)
              NEG1 = L_XNEG(Q,J1,AA,N) * PLM_WT_PXI(J,L)
              NEG2 = L_XNEG(Q,J,AA,N)  * PLM_WT_MXI(J,L)
              SUM_POS = SUM_POS + POS1 + POS2
              SUM_NEG = SUM_NEG + NEG1 + NEG2
            ENDDO
            L_U_HELP_P(Q,L) = SUM_POS
            L_U_HELP_M(Q,L) = SUM_NEG
          ENDDO
        ENDDO

!  Now sum over all harmonic contributions at each user-defined stream

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, NPARS
            SUM_POS = ZERO
            SUM_NEG = ZERO
            DO L = FOURIER, NMOMENTS
              ULP     = PLM_PXUI(UM,L) *   OMEGA_PHASMOMS(N,L)
              L_ULP   = PLM_PXUI(UM,L) * L_OMEGA_PHASMOMS(Q,N,L)
              SUM_POS = SUM_POS +   U_HELP_P(AA,L) * L_ULP + &
                                  L_U_HELP_P(Q,L)  *   ULP
              SUM_NEG = SUM_NEG +   U_HELP_M(AA,L) * L_ULP + &
                                  L_U_HELP_M(Q,L)  *   ULP
            ENDDO
            L_U_XPOS(Q,UM,AA,N) = SUM_POS
            L_U_XNEG(Q,UM,AA,N) = SUM_NEG
          ENDDO
        ENDDO

!  end eigenvector loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LPC_UHOMOG_SOLUTION

!

      SUBROUTINE LPC_BEAM_SOLUTION_ELASTIC &
            ( LAYER, FOURIER, DO_PLANE_PARALLEL, NMOMENTS, NSTREAMS, & ! Inputs
              KS, KF, KPARS, KVARY, NKSTORAGE,                       & ! Inputs
              QUAD_STREAMS, L_OMEGA_PHASMOMS,                        & ! Inputs
              PLM_00_PXI,  PLM_00_MXI, LOCAL_INITRANS,               & ! Inputs
              LOCAL_ASOURCE,    L_LOCAL_ASOURCE,                     & ! Inputs
              QMAT_SAVE, QPIVOT, SAB_SAVE, DAB_SAVE,                 & ! Inputs
              QSUMVEC_SAVE, QDIFVEC_SAVE, QVEC_SAVE, QAUX_SAVE,      & ! Inputs
              L_SAB_SAVE, L_DAB_SAVE, L_EIGENMAT,                    & ! Inputs
              L_WVECTOR, FAIL, MESSAGE )                               ! Outputs

!  This is the Linearization of the classical Chandrasekhar beam solutio
!     Plane parallel or average secant.

!  Output is the single matrix L_WVECTOR of linearized beam solution.

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, ONE, HALF, TWO, MAX_STREAMS, MAX_2_STREAMS, &
                              MAX_LAYERS, MAX_LAYERS_NK, MAX_MOMENTS, MAX_ATMOSWFS

!  Linear-alglebra routine

      USE LRRS_AUX2_m, Only : DGETRS

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  layer index, Fourier number

      INTEGER  , INTENT(IN) :: LAYER
      INTEGER  , INTENT(IN) :: FOURIER

!  Plane parallel flag

      LOGICAL  , INTENT(IN) :: DO_PLANE_PARALLEL

!  Given number of streams and moments

      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: NMOMENTS

!  Linearization control

      INTEGER  , INTENT(IN) :: KS, KF
      LOGICAL  , INTENT(IN) :: KVARY  ( 0:MAX_LAYERS )
      INTEGER  , INTENT(IN) :: KPARS  ( 0:MAX_LAYERS )
      INTEGER  , INTENT(IN) :: NKSTORAGE ( MAX_LAYERS,0:MAX_LAYERS )
!      INTEGER  , INTENT(IN) :: NKSTORAGE ( MAX_LAYERS,0:MAX_LAYERS_NK ) ! Bug JvG 4/2

!  Beam solution inputs
!  --------------------

!  average secant

      REAL(FPK), INTENT(IN) :: LOCAL_ASOURCE
      REAL(FPK), INTENT(IN) :: L_LOCAL_ASOURCE ( MAX_ATMOSWFS, MAX_LAYERS_NK )

!  Local initial transmittance

      REAL(FPK), INTENT(IN) :: LOCAL_INITRANS

!  Quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STREAMS ( MAX_STREAMS )

!  Linearized IOP inputs

      REAL(FPK), INTENT(IN) :: L_OMEGA_PHASMOMS &
                          ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS )

!  Saved matrices from eigenvalue computation

      REAL(FPK), INTENT(IN) :: DAB_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: SAB_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: L_DAB_SAVE ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: L_SAB_SAVE ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: L_EIGENMAT ( MAX_ATMOSWFS, MAX_STREAMS, MAX_STREAMS )

!  Particular solution: LU-decomposed matrix and Pivot

      REAL(FPK), INTENT(IN) :: QMAT_SAVE ( MAX_STREAMS, MAX_STREAMS )
      INTEGER  , INTENT(IN) :: QPIVOT    ( MAX_STREAMS )

!  legendre polynomial functions

      REAL(FPK), INTENT(IN) :: PLM_00_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_00_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Classical beam solution matrices and vectors
!  Beam solution independent of optical depth (classical solution)

      REAL(FPK), INTENT(IN) :: QSUMVEC_SAVE ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QDIFVEC_SAVE ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QVEC_SAVE    ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QAUX_SAVE    ( MAX_STREAMS )

!  output module arguments
!  -----------------------

!  Linearized Particular integral solution  vector
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: L_WVECTOR &
                ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS_NK )

!  output status

      LOGICAL            , INTENT(OUT) :: FAIL
      CHARACTER (LEN=120), INTENT(OUT) :: MESSAGE

!  Local variables
!  ---------------

!  Saved vector

      REAL(FPK) :: QSAVE(MAX_STREAMS)

!  linearization arrays

      REAL(FPK) :: L_QVEC   (MAX_STREAMS,MAX_ATMOSWFS)
      REAL(FPK) :: L_QSUMVEC(MAX_STREAMS,MAX_ATMOSWFS)
      REAL(FPK) :: L_QDIFVEC(MAX_STREAMS,MAX_ATMOSWFS)

!  help variables

      LOGICAL   :: DO_VARY
      INTEGER   :: I, J, L, INFO, I1, NSTR2, Q, N
      INTEGER   :: K, NK, NPARS, LVARY
      REAL(FPK) :: TP, TM, TM1, TM2, TM3, XINV, HELP, WSUMR, WDIF

!  initialise status

      FAIL    = .FALSE.
      MESSAGE = ' '

!  Local Integers

      NSTR2 = 2 * NSTREAMS
      N = LAYER

!  Linearization calculation
!  -------------------------

!  set up driving vector (using saved results)
!  This must be done regardless of whether layer N is varying or not.

      DO I = 1, NSTREAMS
        QSAVE(I) = QDIFVEC_SAVE(I) + TWO * LOCAL_ASOURCE * QVEC_SAVE(I)
      ENDDO

!  Linearization control

      IF ( KS.EQ.0 ) THEN
        NPARS   = KPARS(KS)
        LVARY   = 0
        DO_VARY = .TRUE.
      ELSE
        NPARS   = KPARS(N)
        LVARY   = N
        DO_VARY = KVARY(N)
      ENDIF

!  No particular solution beyond the cutoff layer.
!  Skip this section if there is nothing varying
!    [ Zero the boundary layer values and start Part 2 ]

      IF ( LOCAL_INITRANS .EQ. ZERO .OR. .NOT.DO_VARY) THEN
        NK = NKSTORAGE(N,LVARY)
        DO I = 1, NSTR2
          DO Q = 1, NPARS
            L_WVECTOR(Q,I,NK) = ZERO
          ENDDO
        ENDDO
        GO TO 2222
      ENDIF

!  For each varying parameter

      DO Q = 1, NPARS

!  Set up linearized sum and difference vectors for Beam source terms

        DO I = 1, NSTREAMS
          XINV = ONE / QUAD_STREAMS(I)
          TP = ZERO
          TM = ZERO
          DO L = FOURIER, NMOMENTS
            TP = TP + L_OMEGA_PHASMOMS(Q,N,L) * PLM_00_PXI(I,L)
            TM = TM + L_OMEGA_PHASMOMS(Q,N,L) * PLM_00_MXI(I,L)
          ENDDO
          L_QSUMVEC(I,Q) =  ( TP + TM ) * XINV
          L_QDIFVEC(I,Q) =  ( TP - TM ) * XINV
        ENDDO

!   setup linearized RHS vector
!  ( use results from the original solution )

        DO I = 1, NSTREAMS
          HELP = ZERO
          DO J = 1, NSTREAMS
            TM1 = L_EIGENMAT(Q,I,J) * QVEC_SAVE(J)
            TM2 =   DAB_SAVE(I,J)   * L_QSUMVEC(J,Q)
            TM3 = L_DAB_SAVE(Q,I,J) *   QSUMVEC_SAVE(J)
            HELP = HELP - TM1 - TM2 - TM3
          ENDDO
          L_QVEC(I,Q)  = HELP + L_QDIFVEC(I,Q) * LOCAL_ASOURCE
        ENDDO

      ENDDO

!  additional terms for the quasi-spherical case
!  ( layers greater than one )

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        IF ( N.GT.1 ) THEN
          NK = NKSTORAGE(N,LVARY)
          DO Q = 1, NPARS
            DO I = 1, NSTREAMS
              HELP = QSAVE(I) * L_LOCAL_ASOURCE(Q,NK)
              L_QVEC(I,Q) = L_QVEC(I,Q) + HELP
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Solve problem by back substitution for all of the RHS vectors
!  ( uses L_U decomposition of QMAT_SAVE from original solution )

      CALL DGETRS &
              ('N',NSTREAMS,NPARS,QMAT_SAVE, &
                MAX_STREAMS,QPIVOT,L_QVEC,MAX_STREAMS,INFO)

      IF ( INFO .NE. 0 ) THEN
        FAIL    = .TRUE.
        MESSAGE = 'Elastic Beam P.I. linearization (N,N) (DGETRS)'
        RETURN
      ENDIF

!  assign solutions for the quasi-spherical case, N > 1

      IF ( .NOT. DO_PLANE_PARALLEL .AND. N.GT.1 ) THEN

        NK = NKSTORAGE(N,LVARY)
        DO I = 1, NSTREAMS
          TM3 = - QAUX_SAVE(I) / LOCAL_ASOURCE
          DO Q = 1, NPARS
            I1 = I + NSTREAMS
            HELP = ZERO
            DO J = 1, NSTREAMS
              TM1 =   SAB_SAVE(I,J)   * L_QVEC(J,Q)
              TM2 = L_SAB_SAVE(Q,I,J) *   QVEC_SAVE(J)
              HELP = HELP - TM1 - TM2
            ENDDO
            WSUMR = L_QVEC(I,Q)
            TM2 = ( HELP - L_QSUMVEC(I,Q) ) / LOCAL_ASOURCE
            WDIF = L_LOCAL_ASOURCE(Q,NK) * TM3 + TM2
            L_WVECTOR(Q,I,NK)  = HALF * ( WSUMR + WDIF )
            L_WVECTOR(Q,I1,NK) = HALF * ( WSUMR - WDIF )
          ENDDO
        ENDDO

!  assign solutions for plane/parallel & quasi-spherical case N = 1

      ELSE

        NK = NKSTORAGE(N,LVARY)
        DO I = 1, NSTREAMS
          DO Q = 1, NPARS
            I1 = I + NSTREAMS
            HELP = ZERO
            DO J = 1, NSTREAMS
              TM1 =   SAB_SAVE(I,J)   * L_QVEC(J,Q)
              TM2 = L_SAB_SAVE(Q,I,J) *   QVEC_SAVE(J)
              HELP = HELP - TM1 - TM2
            ENDDO
            WSUMR = L_QVEC(I,Q)
            WDIF = ( HELP - L_QSUMVEC(I,Q) ) / LOCAL_ASOURCE
            L_WVECTOR(Q,I,NK)  = HALF * ( WSUMR + WDIF )
            L_WVECTOR(Q,I1,NK) = HALF * ( WSUMR - WDIF )
          ENDDO
        ENDDO

      ENDIF

!  Part 2.
!  =======

!  Continuation point

 2222 CONTINUE

!  Only for the pseudo-spherical case, profile weighting functions
!  Also not required for the column weighting functions (K = 0)

      IF ( DO_PLANE_PARALLEL )  RETURN
      IF ( KS .EQ. 0 )          RETURN

!  No particular solution beyond the cutoff layer.
!  No solution if layer is inactive
!    [ Zero the boundary layer values and exit ]

      IF ( LOCAL_INITRANS.EQ.ZERO .or. .NOT.DO_VARY ) THEN
        DO K = KS, KF - 1
          NK = NKSTORAGE(N,K)
          DO I = 1, NSTR2
            DO Q = 1, KPARS(K)
              L_WVECTOR(Q,I,NK) = ZERO
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Loop over layers K above N

      DO K = 1, KF - 1

!  If there is a varying layer

        IF ( KVARY(K) ) THEN

!  number of varying parameters for this layer

          NK = NKSTORAGE(N,K)

!  Set up vector

          DO I = 1, NSTREAMS
            DO Q = 1, KPARS(K)
              L_QVEC(I,Q) = QSAVE(I) * L_LOCAL_ASOURCE(Q,NK)
            ENDDO
          ENDDO

!  Solve problem by back substitution for all of the RHS vectors
!  ( uses L_U decomposition of QMAT_SAVE from original solution )

          CALL DGETRS &
              ('N',NSTREAMS,KPARS(K),QMAT_SAVE, &
                MAX_STREAMS,QPIVOT,L_QVEC,MAX_STREAMS,INFO)

          IF ( INFO .NE. 0 ) THEN
            FAIL    = .TRUE.
            MESSAGE = 'Beam P.I. linearization (N,K) (DGETRS)'
            RETURN
          ENDIF

!  assign linearized solutions for layer N due to variations in layer K

          DO I = 1, NSTREAMS
            TM3 = - QAUX_SAVE(I) / LOCAL_ASOURCE
            DO Q = 1, KPARS(K)
              I1 = I + NSTREAMS
              HELP = ZERO
              DO J = 1, NSTREAMS
                HELP = HELP - SAB_SAVE(I,J) * L_QVEC(J,Q)
              ENDDO
              WSUMR = L_QVEC(I,Q)
              WDIF = L_LOCAL_ASOURCE(Q,NK) * TM3 + (HELP/LOCAL_ASOURCE)
              L_WVECTOR(Q,I,NK)  = HALF * ( WSUMR + WDIF )
              L_WVECTOR(Q,I1,NK) = HALF * ( WSUMR - WDIF )
            ENDDO
          ENDDO

!  end K-layer loop

        ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LPC_BEAM_SOLUTION_ELASTIC

!

      SUBROUTINE LPC_UBEAM_SOLUTION_ELASTIC &
           ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,       & ! Inputs
             FOURIER, N, NSTREAMS, NMOMENTS, N_USER_STREAMS,      & ! Inputs
             KS, KF, KPARS, KVARY, NKSTORAGE,                     & ! Inputs
             STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,              & ! Inputs
             PICUTOFF, L_WPARTIC, OMEGAMOMS, L_OMEGAMOMS, W_HELP, & ! Inputs
             PLM_WT_PXI, PLM_WT_MXI, PLM_00_PXUI,                 & ! Inputs
             PLM_PXUI, PLM_00_MXUI, PLM_MXUI,                     & ! Inputs
             L_U_WPOS1, L_U_WPOS2, L_U_WNEG1, L_U_WNEG2 )           ! Outputs

!  Linearized Solar-beam solutions (Elastic field) at user angles

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m, Only : FPK, ZERO, MAX_LAYERS, MAX_LAYERS_NK, MAX_USER_STREAMS, &
                              MAX_STREAMS, MAX_2_STREAMS, MAX_MOMENTS, MAX_ATMOSWFS

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Directional flags

      LOGICAL  , INTENT(IN) :: DO_UPWELLING
      LOGICAL  , INTENT(IN) :: DO_DNWELLING

!  Plane-parallel flag

      LOGICAL  , INTENT(IN) :: DO_PLANE_PARALLEL

!  Computational control

      INTEGER  , INTENT(IN) :: N
      INTEGER  , INTENT(IN) :: FOURIER
      INTEGER  , INTENT(IN) :: NSTREAMS
      INTEGER  , INTENT(IN) :: NMOMENTS

!  Linearization control

      INTEGER  , INTENT(IN) :: KS, KF
      LOGICAL  , INTENT(IN) :: KVARY  ( 0:MAX_LAYERS )
      INTEGER  , INTENT(IN) :: KPARS  ( 0:MAX_LAYERS )
      INTEGER  , INTENT(IN) :: NKSTORAGE ( MAX_LAYERS,0:MAX_LAYERS )
!      INTEGER  , INTENT(IN) :: NKSTORAGE ( MAX_LAYERS,0:MAX_LAYERS_NK ) ! Bug JvG 4/2

!  Layer existence

      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_UP ( MAX_LAYERS )
      LOGICAL  , INTENT(IN) :: STERM_LAYERMASK_DN ( MAX_LAYERS )

!  number of user-defined stream angles

      INTEGER  , INTENT(IN) :: N_USER_STREAMS

!  Beam solution cutoff layer

      INTEGER  , INTENT(IN) :: PICUTOFF

!  IOP input

      REAL(FPK), INTENT(IN) :: OMEGAMOMS ( MAX_LAYERS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS ( MAX_ATMOSWFS,MAX_LAYERS,0:MAX_MOMENTS)

!  Linearized paprticular solution vector

      REAL(FPK), INTENT(IN) :: L_WPARTIC ( MAX_ATMOSWFS,MAX_2_STREAMS,MAX_LAYERS_NK)

!  Legendre functions

      REAL(FPK), INTENT(IN) :: PLM_WT_PXI  ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_WT_MXI  ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_00_PXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_PXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_00_MXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(IN) :: PLM_MXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  Help array

      REAL(FPK), INTENT(IN) :: W_HELP(0:MAX_MOMENTS)

!  subroutine output arguments
!  ---------------------------

!  Linearized User_angle solution vectors
!mick fix 10/19/2015 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: &
            L_U_WPOS1 ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS ), &
            L_U_WPOS2 ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS_NK )

      REAL(FPK), INTENT(INOUT) :: &
            L_U_WNEG1 ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS ), &
            L_U_WNEG2 ( MAX_ATMOSWFS, MAX_USER_STREAMS, MAX_LAYERS_NK )

!  Local variables
!  ---------------

      LOGICAL   :: DO_VARY
      INTEGER   :: UM, J, J1, L, Q, LVARY, NPARS, NN, NK, K
      REAL(FPK) :: SUMR, L_SUMR, POS1, POS2, L_POS, L_NEG
      REAL(FPK) :: L_W_HELP(MAX_ATMOSWFS,0:MAX_MOMENTS)
      REAL(FPK) :: HELP    (0:MAX_MOMENTS)

!  Linearization control

      IF ( KS.EQ.0 ) THEN
        NPARS   = KPARS(KS)
        LVARY   = 0
        DO_VARY = .TRUE.
      ELSE
        NPARS   = KPARS(N)
        LVARY   = N
        DO_VARY = KVARY(N)
      ENDIF

!  Initial storage

      NN = NKSTORAGE(N,LVARY)

!  Check existence of first linearization
!  If no solution or no variation, zero output and go to Part 2

      IF ( .NOT. DO_VARY .OR. N.GT.PICUTOFF ) THEN
        IF ( DO_UPWELLING ) THEN
          DO Q = 1, NPARS
            DO UM = 1, N_USER_STREAMS
              L_U_WPOS1(Q,UM,N)  = ZERO
              L_U_WPOS2(Q,UM,NN) = ZERO
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_DNWELLING ) THEN
          DO Q = 1, NPARS
            DO UM = 1, N_USER_STREAMS
              L_U_WNEG1(Q,UM,N)  = ZERO
              L_U_WNEG2(Q,UM,NN) = ZERO
            ENDDO
          ENDDO
        ENDIF
        GO TO 2222
      ENDIF

!  For each moment & particulate, do inner sum over computational angles
!  (required for both directions)

      DO L = FOURIER, NMOMENTS
        DO Q = 1, NPARS
          L_SUMR = ZERO
          DO  J = 1, NSTREAMS
            J1 = J + NSTREAMS
            POS1 = L_WPARTIC(Q,J1,NN) * PLM_WT_PXI(J,L)
            POS2 = L_WPARTIC(Q,J,NN)  * PLM_WT_MXI(J,L)
            L_SUMR = L_SUMR + POS1 + POS2
          ENDDO
          L_W_HELP(Q,L) = W_HELP(L) * L_OMEGAMOMS(Q,N,L) + &
                             L_SUMR *   OMEGAMOMS(N,L)
        ENDDO
      ENDDO

!  Upwelling: Sum over harmonic contributions at each user-defined strea
!   Two solutions: POS1 = single scatter, POS2 = multiple scatter

      IF ( DO_UPWELLING ) THEN
        IF ( STERM_LAYERMASK_UP(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, NPARS
              L_POS = ZERO
              DO L = FOURIER, NMOMENTS
                L_POS = L_POS + L_OMEGAMOMS(Q,N,L) * PLM_00_PXUI(UM,L)
              ENDDO
              L_U_WPOS1(Q,UM,N) = L_POS
              L_POS = ZERO
              DO L = FOURIER, NMOMENTS
                L_POS = L_POS + L_W_HELP(Q,L) * PLM_PXUI(UM,L)
              ENDDO
              L_U_WPOS2(Q,UM,NN) = L_POS
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Downwelling: Sum over harmonic contributions at each user-defined str
!   Two solutions: NEG1 = single scatter, NEG2 = multiple scatter

      IF ( DO_DNWELLING ) THEN
        IF ( STERM_LAYERMASK_DN(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, NPARS
              L_NEG = ZERO
              DO L = FOURIER, NMOMENTS
                L_NEG = L_NEG + L_OMEGAMOMS(Q,N,L) * PLM_00_MXUI(UM,L)
              ENDDO
              L_U_WNEG1(Q,UM,N) = L_NEG
              L_NEG = ZERO
              DO L = FOURIER, NMOMENTS
                L_NEG = L_NEG + L_W_HELP(Q,L) * PLM_MXUI(UM,L)
              ENDDO
              L_U_WNEG2(Q,UM,NN) = L_NEG
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Continuation point

 2222 CONTINUE

!  Only these extra variations when the following conditions
!  are not satisfied (pseudo-spherical, layer > 1)
!  Profile linearization only - Skip otherwise

      IF ( N .EQ. 1 )          RETURN
      IF ( DO_PLANE_PARALLEL ) RETURN
      IF ( KS .EQ. 0 )         RETURN

!  If no solution, zero output and skip to next layer (exit)

      IF ( N .GT. PICUTOFF ) THEN
        IF ( DO_UPWELLING ) THEN
          DO K = KS, KF - 1
            NK    = NKSTORAGE(N,K)
            DO Q = 1, KPARS(K)
              DO UM = 1, N_USER_STREAMS
                L_U_WPOS2(Q,UM,NK) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_DNWELLING ) THEN
          DO K = KS, KF - 1
            NK    = NKSTORAGE(N,K)
            DO Q = 1, KPARS(K)
              DO UM = 1, N_USER_STREAMS
                L_U_WNEG2(Q,UM,NK) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        RETURN
      ENDIF

!  start loop over all layers above N

      DO K = 1, N - 1

!  Storage index

        NK = NKSTORAGE(N,K)

!  only do if layer K has some variation

        IF ( KVARY(K) ) THEN
          DO Q = 1, KPARS(K)

!  For each moment do inner sum over computational angles
!  taking care to use chain rule on the linearization

            DO L = FOURIER, NMOMENTS
              SUMR = ZERO
              DO  J = 1, NSTREAMS
                J1 = J + NSTREAMS
                POS1 = L_WPARTIC(Q,J1,NK) * PLM_WT_PXI(J,L)
                POS2 = L_WPARTIC(Q,J,NK)  * PLM_WT_MXI(J,L)
                SUMR = SUMR + POS1 + POS2
              ENDDO
              HELP(L) =  SUMR * OMEGAMOMS(N,L)
            ENDDO

!  Now sum over all harmonic contributions at each user-defined stream

            IF ( DO_UPWELLING ) THEN
              IF ( STERM_LAYERMASK_UP(N) ) THEN
                 DO UM = 1, N_USER_STREAMS
                   L_POS = ZERO
                   DO L = FOURIER, NMOMENTS
                     L_POS = L_POS + HELP(L) * PLM_PXUI(UM,L)
                   ENDDO
                   L_U_WPOS2(Q,UM,NK) = L_POS
                 ENDDO
              ENDIF
            ENDIF
            IF ( DO_DNWELLING ) THEN
              IF ( STERM_LAYERMASK_DN(N) ) THEN
                 DO UM = 1, N_USER_STREAMS
                   L_NEG = ZERO
                   DO L = FOURIER, NMOMENTS
                     L_NEG = L_NEG + HELP(L) * PLM_MXUI(UM,L)
                   ENDDO
                   L_U_WNEG2(Q,UM,NK) = L_NEG
                 ENDDO
              ENDIF
            ENDIF

!  end parameter and K-layer loops

          ENDDO
        ENDIF
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LPC_UBEAM_SOLUTION_ELASTIC

!

      SUBROUTINE LPC_GREENFUNC_SOLUTION_1 &
            ( FOURIER, LAYER, SOLUTION, DO_MULTIPLIER,        & ! Inputs
              TAYLOR_ORDER, TAYLOR_SMALL, KS, KF, KPARS,      & ! Inputs
              NKSTORAGE, NKSTORAGE2, NSTREAMS, QUAD_WEIGHTS,  & ! Inputs
              XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE,           & ! Inputs
              ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P,         & ! Inputs
              PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS, & ! Inputs
              DELTAS, QSOURCE_QPOS, QSOURCE_QNEG,             & ! Inputs
              L_XPOS, L_KEIGEN, L_KTRANS, L_EIGENNORM_SAVE,   & ! Inputs
              L_ASOURCE, L_DELTRANS, L_INITRANS,              & ! Inputs
              L_DELTAS, L_QSOURCE_QPOS, L_QSOURCE_QNEG,       & ! Inputs
              L_RAMAN_WUPPER, L_RAMAN_WLOWER,                 & ! Outputs
              L_ATERM_SAVE, L_BTERM_SAVE )                      ! Outputs

!  This is the Green's function solution in one layer, for one contribution
!   to the series of particular integrals

!  Module of dimensions and numbers
!  ---------------------------------

      USE LRRS_PARS_m   , Only : FPK, ZERO, ONE, UPIDX, DNIDX, MAX_ATMOSWFS, &
                                 MAX_STREAMS, MAX_2_STREAMS, MAX_LAYERS, MAX_LAYERS_SQ

!  Taylor serioes module
!  ---------------------

      USE lrrs_Taylor_m, Only : Taylor_series_L_1

      IMPLICIT NONE

!  INPUT ARGUMENTS
!  ===============

!  input Control arguments
!  -----------------------

!  Fourier index

      INTEGER  , INTENT(IN) :: FOURIER

!  Given layer, solution index

      INTEGER  , INTENT(IN) :: LAYER
      INTEGER  , INTENT(IN) :: SOLUTION

!  Control for calculating Green's function multiplier

      LOGICAL  , INTENT(INOUT) :: DO_MULTIPLIER

!  Small number limit, Taylor order (Taylor series expansions)

      INTEGER  , INTENT(IN) :: TAYLOR_ORDER
      REAL(FPK), INTENT(IN) :: TAYLOR_SMALL

!  quadrature inputs

      INTEGER  , INTENT(IN) :: NSTREAMS
      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS ( MAX_STREAMS )

!  Number of weighting function control

      INTEGER  , INTENT(IN) :: KS, KF, KPARS(0:MAX_LAYERS)
      INTEGER  , INTENT(IN) :: NKSTORAGE (MAX_LAYERS,0:MAX_LAYERS)
      INTEGER  , INTENT(IN) :: NKSTORAGE2(MAX_LAYERS,0:MAX_LAYERS)

!  source vector inputs
!  --------------------

!  layer cutoff

      INTEGER  , INTENT(IN) :: PICUTOFF

!  Flipper for using -x instead of x as integration argument

      LOGICAL  , INTENT(IN) :: FLIPPER

!  Beam: Average secant, layer transmittance,Transmittance to start

      REAL(FPK), INTENT(IN) :: ASOURCE
      REAL(FPK), INTENT(IN) :: DELTRANS
      REAL(FPK), INTENT(IN) :: INITRANS

!  optical thickness (only required for Taylor series expansions)

      REAL(FPK), INTENT(IN) :: DELTAS

!  Source vectors in quadrature directions

      REAL(FPK), INTENT(IN) :: QSOURCE_QPOS  ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QSOURCE_QNEG  ( MAX_STREAMS )

!  Atmospheric Linearization Source vector inputs
!  ----------------------------------------------

      REAL(FPK), INTENT(IN) :: L_ASOURCE  ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_DELTRANS ( MAX_ATMOSWFS, 0:MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_INITRANS ( MAX_ATMOSWFS, 0:MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_DELTAS   ( MAX_ATMOSWFS )

      REAL(FPK), INTENT(IN) :: L_QSOURCE_QPOS &
                ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: L_QSOURCE_QNEG &
                ( MAX_ATMOSWFS, 0:MAX_LAYERS, MAX_STREAMS )

!  homogeneous solution inputs
!  ---------------------------

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Eigenvalues and whole-layer transmittance factors

      REAL(FPK), INTENT(IN) :: KTRANS    ( MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: KEIGEN    ( MAX_STREAMS, MAX_LAYERS )

!  Norm values

      REAL(FPK), INTENT(IN) :: EIGENNORM_SAVE ( MAX_STREAMS, MAX_LAYERS)

!  Green's function variables
!  --------------------------

!  multipliers, from the regular calculation

      REAL(FPK), INTENT(IN) :: MULT_M ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: MULT_P ( MAX_STREAMS )

!  Saved ATERM and BTERM vectors, from the regular calculation

      REAL(FPK), INTENT(IN) :: ATERM_SAVE ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: BTERM_SAVE ( MAX_STREAMS )

!  Linearized homogeneous solution inputs
!  --------------------------------------

      REAL(FPK), INTENT(IN) :: L_XPOS &
              ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_KTRANS ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(IN) :: L_KEIGEN ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS )

      REAL(FPK), INTENT(IN) :: L_EIGENNORM_SAVE ( MAX_ATMOSWFS, MAX_STREAMS, MAX_LAYERS)

!  output arguments
!  ================
!mick fix 10/19/2015 - changed intent on both pairs to inout

!  Linearized solution variables
!  -----------------------------

      REAL(FPK), INTENT(INOUT) :: L_RAMAN_WUPPER &
                ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS_SQ )
      REAL(FPK), INTENT(INOUT) :: L_RAMAN_WLOWER &
                ( MAX_ATMOSWFS, MAX_2_STREAMS, MAX_LAYERS_SQ )

!  Saved ATERM and BTERM  vectors

      REAL(FPK), INTENT(INOUT) :: L_ATERM_SAVE (MAX_ATMOSWFS,0:MAX_LAYERS,MAX_STREAMS)
      REAL(FPK), INTENT(INOUT) :: L_BTERM_SAVE (MAX_ATMOSWFS,0:MAX_LAYERS,MAX_STREAMS)

!  local Variables
!  ===============

!  Multiplier arrays

      REAL(FPK) :: GFUNC_UP   ( MAX_STREAMS )
      REAL(FPK) :: GFUNC_DN   ( MAX_STREAMS )
      REAL(FPK) :: CFUNC      ( MAX_STREAMS )
      REAL(FPK) :: DFUNC      ( MAX_STREAMS )
      REAL(FPK) :: L_CFUNC    ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_DFUNC    ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_GFUNC_UP ( MAX_ATMOSWFS, MAX_STREAMS )
      REAL(FPK) :: L_GFUNC_DN ( MAX_ATMOSWFS, MAX_STREAMS )

!   RRS Particular integral at layer boundaries (one term)

      REAL(FPK) :: L_WLEVELS ( MAX_ATMOSWFS, MAX_2_STREAMS, 2 )

!  local variables

      INTEGER   :: AA, I, I1, N, Q, M, K, NK, NK2, NSTR2
      REAL(FPK) :: TPA, TMA, SUM_LA, SUM_LB
      REAL(FPK) :: S_P_U, S_P_L, S_M_U, S_M_L
      REAL(FPK) :: RHO_M, RHO_P, NORM, ZDEL

      REAL(FPK) :: L_RHO_M, L_RHO_P, L_MULT_M, L_MULT_P
      REAL(FPK) :: L_ZDEL, L_ZWDEL, L_TRAN, L_NORM, L_EVAL, L_DELT

!  Layer variables

      N     = LAYER
      NSTR2 = 2 * NSTREAMS

!  Temporary setting

      DO_MULTIPLIER = .TRUE.

!  Fourier (for debug only)

      M = FOURIER

!  Start loop over all weighting functions

      DO K = KS, KF

!  Storage value

        NK  = NKSTORAGE(N,K)
        NK2 = NKSTORAGE2(N,K)

!  Zero the TOTAL particular integrals at layer boudaries
!     if this is the FIRST solution
!mick fix 7/20/2016 - commented out (initializing already done in L_SOURCES_MASTER_1)

        !IF ( SOLUTION .EQ. 1 ) THEN
        !  DO I = 1, NSTR2
        !    DO Q = 1, KPARS(K)
        !      L_RAMAN_WUPPER(Q,I,NK2) = ZERO
        !      L_RAMAN_WLOWER(Q,I,NK2) = ZERO
        !    ENDDO
        !  ENDDO
        !ENDIF

!  No particular solution beyond the cutoff layer.
!    [ Zero the boundary layer values and skip to next WF ]

        IF ( N .GT. PICUTOFF ) THEN
          DO I = 1, NSTR2
            DO Q = 1, KPARS(K)
              L_WLEVELS(Q,I,UPIDX) = ZERO
              L_WLEVELS(Q,I,DNIDX) = ZERO
            ENDDO
          ENDDO
          GO TO 4567
        ENDIF

!  For each eigenstream, get Linearized ATERM_SAVE and BTERM_SAVE
!  ==============================================================

!  Linearize ATERM_SAVE and BTERM_SAVE, straightforward chain rule.
!   - For K = 0 (column) or K = N (profile), no cross-layer derivatives

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO Q = 1, KPARS(K)
            DO AA = 1, NSTREAMS
              SUM_LA = ZERO
              SUM_LB = ZERO
              DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                TPA = L_QSOURCE_QPOS(Q,K,I) *   XPOS(I,AA,N)    + &
                        QSOURCE_QPOS(I)     * L_XPOS(Q,I,AA,N)  + &
                      L_QSOURCE_QNEG(Q,K,I) *   XPOS(I1,AA,N)   + &
                        QSOURCE_QNEG(I)     * L_XPOS(Q,I1,AA,N)
                TMA = L_QSOURCE_QNEG(Q,K,I) *   XPOS(I,AA,N)    + &
                        QSOURCE_QNEG(I)     * L_XPOS(Q,I,AA,N)  + &
                      L_QSOURCE_QPOS(Q,K,I) *   XPOS(I1,AA,N)   + &
                        QSOURCE_QPOS(I)     * L_XPOS(Q,I1,AA,N)
                SUM_LA  = SUM_LA + TPA * QUAD_WEIGHTS(I)
                SUM_LB  = SUM_LB + TMA * QUAD_WEIGHTS(I)
              ENDDO
              NORM   =   EIGENNORM_SAVE(AA,N)
              L_NORM = L_EIGENNORM_SAVE(Q,AA,N)
              L_ATERM_SAVE(Q,K,AA) = &
                 ( SUM_LA - L_NORM * ATERM_SAVE(AA) ) / NORM
              L_BTERM_SAVE(Q,K,AA) = &
                 ( SUM_LB - L_NORM * BTERM_SAVE(AA) ) / NORM
            ENDDO
          ENDDO
        ENDIF

!  Linearize ATERM_SAVE and BTERM_SAVE, straightforward chain rule.
!   - For K > 0 and K =/ N (profiles), cross-layer derivatives

        IF ( N.NE.K .AND. K.GT.0 ) THEN
          DO Q = 1, KPARS(K)
            DO AA = 1, NSTREAMS
              SUM_LA = ZERO
              SUM_LB = ZERO
              DO I = 1, NSTREAMS
                I1 = I + NSTREAMS
                TPA = L_QSOURCE_QPOS(Q,K,I) *   XPOS(I,AA,N)  + &
                      L_QSOURCE_QNEG(Q,K,I) *   XPOS(I1,AA,N)
                TMA = L_QSOURCE_QNEG(Q,K,I) *   XPOS(I,AA,N)  + &
                      L_QSOURCE_QPOS(Q,K,I) *   XPOS(I1,AA,N)
                SUM_LA  = SUM_LA + TPA * QUAD_WEIGHTS(I)
                SUM_LB  = SUM_LB + TMA * QUAD_WEIGHTS(I)
              ENDDO
              NORM   =   EIGENNORM_SAVE(AA,N)
              L_ATERM_SAVE(Q,K,AA) = SUM_LA / NORM
              L_BTERM_SAVE(Q,K,AA) = SUM_LB / NORM
            ENDDO
          ENDDO
        ENDIF

!  Optical depth integrations, Linearized discrete ordinate solution
!  =================================================================

!  Avoid this multiplier calculation if already done
!     ------------ DISABLE THIS FEATURE FOR NOW.

!   ---> Allows for degeneracy (old code, this was treated separately)
!   ---> Allows for the Flipper formalism (now fully worked out)
!   ---> The small numbers expansion (from Fluorescence model)

!  ------------------------------------- Feature disabled
!      IF ( DO_MULTIPLIER ) THEN
!  ------------------------------------- Feature disabled

        DO AA = 1, NSTREAMS

!  Always require these

          RHO_M = ASOURCE - KEIGEN(AA,N)
          RHO_P = ASOURCE + KEIGEN(AA,N)
          ZDEL  = KTRANS(AA,N)

!  Set linearized solution multipliers (if flagged)
!  ------------------------------------------------

!  start parameter loop

           DO Q = 1, KPARS(K)

!  Linearized transmittance factors and RHO

            IF ( N.EQ.K .OR. K.EQ.0 ) THEN
              L_ZDEL  = L_KTRANS(Q,AA,N)
              L_ZWDEL = L_ZDEL*DELTRANS + ZDEL*L_DELTRANS(Q,K)
              L_EVAL  = L_KEIGEN(Q,AA,N)
              L_DELT  = l_deltas(q)
              L_RHO_M = L_ASOURCE(Q,K) - L_EVAL
              L_RHO_P = L_ASOURCE(Q,K) + L_EVAL
            ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN
              L_ZDEL  = ZERO
              L_ZWDEL = ZDEL * L_DELTRANS(Q,K)
              L_EVAL  = ZERO
              L_DELT  = ZERO
              L_RHO_M = L_ASOURCE(Q,K)
              L_RHO_P = L_ASOURCE(Q,K)
            ENDIF

!  Set the downwelling multiplier (small-number analysis may apply)

! (2.3) call limit_l_gcfunc ( rho_m, deltas, keigen(aa,n), zdel, l_asource(q,k), l_eval, l_delt, l_mult_m )
! (2.5) mick fix 7/20/2016 - modified three inputs from the following
              !CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, -RHO_M, DELTAS, L_DELT, L_EVAL, L_ASOURCE(Q,K), ZDEL, ONE, L_MULT_M )

            IF ( ABS(RHO_M).LT.TAYLOR_SMALL ) THEN
              CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, RHO_M, DELTAS, L_DELT, L_EVAL, L_ASOURCE(Q,K), DELTRANS, ASOURCE, L_MULT_M )
            ELSE
              L_MULT_M = ( L_ZDEL - L_DELTRANS(Q,K) - MULT_M(AA) * L_RHO_M ) / RHO_M
            ENDIF

!  Set the Upwelling multiplier

            L_MULT_P = - ( L_ZWDEL + MULT_P(AA) * L_RHO_P ) / RHO_P

!  Distinguish between the flipper cases

            L_TRAN = L_INITRANS(Q,K)
            IF ( .NOT.FLIPPER ) THEN
              L_CFUNC(Q,AA)  = L_MULT_M * INITRANS + MULT_M(AA) * L_TRAN
              L_DFUNC(Q,AA)  = L_MULT_P * INITRANS + MULT_P(AA) * L_TRAN
            ELSE
              L_CFUNC(Q,AA)  = L_MULT_P * INITRANS + MULT_P(AA) * L_TRAN
              L_DFUNC(Q,AA)  = L_MULT_M * INITRANS + MULT_M(AA) * L_TRAN
            ENDIF

!  End parameter Q-loop and eigenstream AA-loop

          ENDDO
        ENDDO

!   Set local CFUNC, DFUNC Distinguish between the flipper cases
!         (multipliers are reversed if -x rather than +x is tau-variable

        DO AA = 1, NSTREAMS
          IF ( .NOT.FLIPPER ) THEN
            CFUNC(AA)  = MULT_M(AA) * INITRANS
            DFUNC(AA)  = MULT_P(AA) * INITRANS
          ELSE
            CFUNC(AA)  = MULT_P(AA) * INITRANS
            DFUNC(AA)  = MULT_M(AA) * INITRANS
          ENDIF
        ENDDO

!  Local linearized Green function multipliers
!  ===========================================

!  For Column WFS (K=0) or self-correlated profile WFs (K=N),
!      Set local Green's function multiplierS GFUNC_UP, GFUNC_DN

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO AA = 1, NSTREAMS
            GFUNC_DN(AA) = CFUNC(AA) * ATERM_SAVE(AA)
            GFUNC_UP(AA) = DFUNC(AA) * BTERM_SAVE(AA)
          ENDDO
        ENDIF

!  All WFS, same calculation

        DO AA = 1, NSTREAMS
          DO Q = 1, KPARS(K)
            L_GFUNC_DN(Q,AA) = L_CFUNC(Q,AA) *   ATERM_SAVE(AA) + &
                                 CFUNC(AA)   * L_ATERM_SAVE(Q,K,AA)
            L_GFUNC_UP(Q,AA) = L_DFUNC(Q,AA) *   BTERM_SAVE(AA) + &
                                 DFUNC(AA)   * L_BTERM_SAVE(Q,K,AA)
          ENDDO
        ENDDO

!  Set Linearized particular integral from Green function expansion
!  ================================================================

!  Column WFS (K=0) or self-correlated profile WFs (K=N)
!   Must include variation of Eigensolutions

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO Q = 1, KPARS(K)
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              S_P_U = ZERO
              S_P_L = ZERO
              S_M_U = ZERO
              S_M_L = ZERO
              DO AA = 1, NSTREAMS
                S_P_U = S_P_U + L_GFUNC_UP(Q,AA) *   XPOS(I1,AA,N) &
                              +   GFUNC_UP(AA)   * L_XPOS(Q,I1,AA,N)
                S_M_U = S_M_U + L_GFUNC_UP(Q,AA) *   XPOS(I,AA,N) &
                              +   GFUNC_UP(AA)   * L_XPOS(Q,I,AA,N)
                S_P_L = S_P_L + L_GFUNC_DN(Q,AA) *   XPOS(I,AA,N) &
                              +   GFUNC_DN(AA)   * L_XPOS(Q,I,AA,N)
                S_M_L = S_M_L + L_GFUNC_DN(Q,AA) *   XPOS(I1,AA,N) &
                              +   GFUNC_DN(AA)   * L_XPOS(Q,I1,AA,N)
              ENDDO
              L_WLEVELS(Q,I,UPIDX)  = S_P_U
              L_WLEVELS(Q,I1,UPIDX) = S_M_U
              L_WLEVELS(Q,I1,DNIDX) = S_M_L
              L_WLEVELS(Q,I,DNIDX)  = S_P_L
            ENDDO
          ENDDO
        ENDIF

!  Cross-layer profile WFs (K /= N, K > 0 )

        IF ( N.NE.K .AND. K.GT.0 ) THEN
          DO Q = 1, KPARS(K)
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              S_P_U = ZERO
              S_P_L = ZERO
              S_M_U = ZERO
              S_M_L = ZERO
              DO AA = 1, NSTREAMS
                S_P_U = S_P_U + L_GFUNC_UP(Q,AA) *   XPOS(I1,AA,N)
                S_M_U = S_M_U + L_GFUNC_UP(Q,AA) *   XPOS(I,AA,N)
                S_P_L = S_P_L + L_GFUNC_DN(Q,AA) *   XPOS(I,AA,N)
                S_M_L = S_M_L + L_GFUNC_DN(Q,AA) *   XPOS(I1,AA,N)
              ENDDO
              L_WLEVELS(Q,I,UPIDX)  = S_P_U
              L_WLEVELS(Q,I1,UPIDX) = S_M_U
              L_WLEVELS(Q,I1,DNIDX) = S_M_L
              L_WLEVELS(Q,I,DNIDX)  = S_P_L
            ENDDO
          ENDDO

        ENDIF

!  Continuation point for avoiding this WF calculation

 4567   CONTINUE

!  Add single RRS particular integral to the RRS-summed total P.I.
!  ---------------------------------------------------------------

        DO Q = 1, KPARS(K)
          DO I = 1, 2 * NSTREAMS
            L_RAMAN_WUPPER(Q,I,NK2) = &
                L_RAMAN_WUPPER(Q,I,NK2) + L_WLEVELS(Q,I,UPIDX)
            L_RAMAN_WLOWER(Q,I,NK2) = &
                L_RAMAN_WLOWER(Q,I,NK2) + L_WLEVELS(Q,I,DNIDX)
          ENDDO
         ENDDO

!  Note the Debug Write -> Very potent, gets all sources
!    output when k = KFD

!        i   = 4 ;  q   = 1
!        if (k.eq.0)write(55,'(3i5,1p2e24.14)')solution,n,i,
!     &             L_WLEVELS(Q,I,UPIDX),L_WLEVELS(Q,I,DNIDX)

!  End Jacobian loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LPC_GREENFUNC_SOLUTION_1

!

      SUBROUTINE LS_GREENFUNC_SOLUTION_1 &
            ( FOURIER, LAYER, SOLUTION, NSTREAMS, N_SURFACEWFS, QUAD_WEIGHTS, & ! Inputs
              FLIPPER, PICUTOFF, INITRANS, XPOS, EIGENNORM_SAVE,              & ! Inputs
              LS_QSOURCE_QPOS, LS_QSOURCE_QNEG, MULT_M, MULT_P,               & ! Inputs
              LS_RAMAN_WUPPER, LS_RAMAN_WLOWER, LS_ATERM, LS_BTERM )            ! Outputs

!  Addition of a number of surface weighting functions
!   Rob Fix 10/10/15

!  This is the Green's function solution, one layer contribution
!    - Solution for the albedo linearization

!  Module of dimensions and numbers
!  --------------------------------

      USE LRRS_PARS_m    , Only : FPK, ZERO, ONE, UPIDX, DNIDX, MAX_SURFACEWFS, &
                                  MAX_STREAMS, MAX_2_STREAMS, MAX_LAYERS

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  Fourier Index

      INTEGER  , INTENT(IN) :: FOURIER

!  Given layer, solution index

      INTEGER  , INTENT(IN) :: LAYER
      INTEGER  , INTENT(IN) :: SOLUTION

!  quadrature inputs

      INTEGER  , INTENT(IN) :: NSTREAMS
      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS ( MAX_STREAMS )

!  Number of weighting functions

      INTEGER  , INTENT(IN) :: N_SURFACEWFS

!  Cut off layer

      INTEGER  , INTENT(IN) :: PICUTOFF

!  Flipper for using -x instead of x as integration argument

      LOGICAL  , INTENT(IN) :: FLIPPER

!  Initial transmittance

      REAL(FPK), INTENT(IN) :: INITRANS

!  source vector inputs
!  --------------------

!  SAVED MULTIPLIERS

      REAL(FPK), INTENT(IN) :: MULT_M  ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: MULT_P  ( MAX_STREAMS )

!  Source vectors in quadrature directions

      REAL(FPK), INTENT(IN) :: LS_QSOURCE_QPOS  ( MAX_SURFACEWFS, MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: LS_QSOURCE_QNEG  ( MAX_SURFACEWFS, MAX_STREAMS )

!  homogeneous solution inputs
!  ---------------------------

!  Eigenvector solutions

      REAL(FPK), INTENT(IN) :: XPOS ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Norm values

      REAL(FPK), INTENT(IN) :: EIGENNORM_SAVE ( MAX_STREAMS, MAX_LAYERS)

!  output module arguments
!  -----------------------

!  Total RRS-summed Particular integral at layer boundaries
!mick fix 7/20/2016 - changed intent to inout

      REAL(FPK), INTENT(INOUT) :: LS_RAMAN_WUPPER  ( MAX_SURFACEWFS, MAX_2_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(INOUT) :: LS_RAMAN_WLOWER  ( MAX_SURFACEWFS, MAX_2_STREAMS, MAX_LAYERS )

!  Surface-linearized ATERM and BTERM vectors

      REAL(FPK), INTENT(OUT) :: LS_ATERM ( MAX_SURFACEWFS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: LS_BTERM ( MAX_SURFACEWFS, MAX_STREAMS )

!  local variables
!  ===============

      INTEGER   :: AA, I, I1, N, M, NSTR2, NS
      REAL(FPK) :: TPA, TMA, SUM_LA, SUM_LB
      REAL(FPK) :: S_P_U, S_P_L, S_M_U, S_M_L, CFUNC, DFUNC

      REAL(FPK) :: GFUNC_UP ( MAX_SURFACEWFS, MAX_STREAMS )
      REAL(FPK) :: GFUNC_DN ( MAX_SURFACEWFS, MAX_STREAMS )

!   RRS Particular integral at layer boundaries (one term)

      REAL(FPK) :: WLEVELS ( MAX_SURFACEWFS, MAX_2_STREAMS, 2 )

!  Zero the solution initialise boundary value

      N = LAYER
      NSTR2 = 2 * NSTREAMS

!mick fix 7/20/2016 - commented out (initializing already done in L_SOURCES_MASTER_1)
      !IF ( SOLUTION .EQ. 1 ) THEN
      !  DO I = 1, NSTR2
      !    LS_RAMAN_WUPPER(1:N_SURFACEWFS,I,LAYER) = ZERO
      !    LS_RAMAN_WLOWER(1:N_SURFACEWFS,I,LAYER) = ZERO
      !  ENDDO
      !ENDIF

!  Fourier (for debug only)

      M = FOURIER

!  No particular solution beyond the cutoff layer.
!    [ Zero the boundary layer values and exit )

      IF ( N .GT. PICUTOFF ) THEN
        DO I = 1, NSTR2
          WLEVELS(1:N_SURFACEWFS,I,UPIDX) = ZERO
          WLEVELS(1:N_SURFACEWFS,I,DNIDX) = ZERO
        ENDDO
        GO TO 4567
      ENDIF

!  For each eigenstream, get the terms ATERM_SAVE and BTERM_SAVE

      DO NS = 1, N_SURFACEWFS
        DO AA = 1, NSTREAMS
          SUM_LA = ZERO
          SUM_LB = ZERO
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            TPA = LS_QSOURCE_QPOS(NS,I) * XPOS(I,AA,N)  + &
                  LS_QSOURCE_QNEG(NS,I) * XPOS(I1,AA,N)
            TMA = LS_QSOURCE_QNEG(NS,I) * XPOS(I,AA,N)  + &
                  LS_QSOURCE_QPOS(NS,I) * XPOS(I1,AA,N)
            SUM_LA  = SUM_LA + TPA * QUAD_WEIGHTS(I)
            SUM_LB  = SUM_LB + TMA * QUAD_WEIGHTS(I)
          ENDDO
          LS_ATERM(NS,AA) = SUM_LA / EIGENNORM_SAVE(AA,N)
          LS_BTERM(NS,AA) = SUM_LB / EIGENNORM_SAVE(AA,N)
        ENDDO
      ENDDO

!  Local Green function multipliers

      NS = N_SURFACEWFS
      DO AA = 1, NSTREAMS
        IF ( .NOT.FLIPPER ) THEN
          CFUNC  = MULT_M(AA) * INITRANS
          DFUNC  = MULT_P(AA) * INITRANS
        ELSE
          CFUNC  = MULT_P(AA) * INITRANS
          DFUNC  = MULT_M(AA) * INITRANS
        ENDIF
        GFUNC_DN(1:NS,AA) = CFUNC * LS_ATERM(1:NS,AA)
        GFUNC_UP(1:NS,AA) = DFUNC * LS_BTERM(1:NS,AA)
      ENDDO

!  Set particular integral from Green function expansion
!    particular integral WLEVELS at lower and upper boundaries

      DO NS = 1, N_SURFACEWFS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          S_P_U = ZERO ; S_P_L = ZERO ; S_M_U = ZERO ; S_M_L = ZERO
          DO AA = 1, NSTREAMS
            S_P_U = S_P_U + GFUNC_UP(NS,AA)*XPOS(I1,AA,N)
            S_M_U = S_M_U + GFUNC_UP(NS,AA)*XPOS(I,AA,N)
            S_P_L = S_P_L + GFUNC_DN(NS,AA)*XPOS(I,AA,N)
            S_M_L = S_M_L + GFUNC_DN(NS,AA)*XPOS(I1,AA,N)
          ENDDO
          WLEVELS(NS,I,UPIDX)  = S_P_U
          WLEVELS(NS,I1,UPIDX) = S_M_U
          WLEVELS(NS,I1,DNIDX) = S_M_L
          WLEVELS(NS,I,DNIDX)  = S_P_L
        ENDDO
      ENDDO

!  Continuation point for avoiding this calculation

 4567 CONTINUE

!  Add single RRS particular integral to the RRS-summed total P.I.
 
      DO NS = 1, N_SURFACEWFS
        DO I = 1, NSTR2
          LS_RAMAN_WUPPER(NS,I,N) = LS_RAMAN_WUPPER(NS,I,N) + WLEVELS(NS,I,UPIDX)
          LS_RAMAN_WLOWER(NS,I,N) = LS_RAMAN_WLOWER(NS,I,N) + WLEVELS(NS,I,DNIDX)
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LS_GREENFUNC_SOLUTION_1

!  End Module

      END MODULE lrrs_L_rtsolutions_1_m

