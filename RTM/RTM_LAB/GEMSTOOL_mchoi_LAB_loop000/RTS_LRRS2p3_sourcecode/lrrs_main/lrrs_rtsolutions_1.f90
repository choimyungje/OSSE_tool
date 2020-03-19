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
! #    Other setup SUBROUTINES :                                #
! #           (1 = Whole-layer, 2 = Part-layer)                 #
! #                                                             #
! #             ATTENUATION_SETUP_1                             #
! #             LEGENDRE_SETUP                                  #
! #             DBEAM_SETUP                                     #
! #             UDBEAM_SETUP                                    #
! #                                                             #
! #             CHAPMAN_FUNCTION                                #
! #                                                             #
! ###############################################################


      MODULE lrrs_rtsolutions_1

      PRIVATE
      PUBLIC :: HOMOG_SOLUTION,&
                UHOMOG_SOLUTION,&
                BEAMSOLUTION_MATRIX,&
                BEAM_SOLUTION_MATRIX,&
                BEAM_SOLUTION_ELASTIC,&
                UBEAM_SOLUTION_ELASTIC,&
                GREENFUNC_SOLUTION_1,&
                ATTENUATION_SETUP_1,&
                LEGENDRE_SETUPS,&
                DBEAM_SETUP,&
                UDBEAM_SETUP,&
                CHAPMAN_FUNCTION

      CONTAINS

      SUBROUTINE HOMOG_SOLUTION &
            ( LAYER, FOURIER, NSTREAMS, NMOMENTS, &
              OMEGA_PHASMOMS, DELTA_TAU, &
              QUAD_STREAMS, QUAD_WEIGHTS, &
              PLM_PXI,     PLM_MXI, &
              XPOS, XNEG, KEIGEN, KTRANS, &
              EIGENMAT_SAVE, EIGENVEC_SAVE, &
              DIFVEC_SAVE, DAB_SAVE, SAB_SAVE, &
              FAIL, MESSAGE )

!  Numerical solution of Eigenproblem.

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS
      USE LRRS_AUX2

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  Given layer index and Fourier number

      INTEGER, INTENT(IN) ::          LAYER
      INTEGER, INTENT(IN) ::          FOURIER

!  Number of streams and moments

      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          NMOMENTS

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

      REAL(FPK), INTENT(INOUT) :: XPOS &
        ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )
      REAL(FPK), INTENT(INOUT) :: XNEG &
        ( MAX_2_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Saved matrices for eigenvalue computation

      REAL(FPK), INTENT(OUT) :: DAB_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: SAB_SAVE ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: EIGENMAT_SAVE  ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: EIGENVEC_SAVE  ( MAX_STREAMS, MAX_STREAMS )
      REAL(FPK), INTENT(OUT) :: DIFVEC_SAVE    ( MAX_STREAMS, MAX_STREAMS )

!  output status

      LOGICAL, INTENT(OUT) ::          FAIL
      CHARACTER (LEN=70), INTENT(OUT) ::     MESSAGE

!  Local variables
!  ---------------

!  local matrix for eigenvalue computation

      REAL(FPK) :: EIGENMAT ( MAX_STREAMS, MAX_STREAMS )

!  (output from Eigenpackage module ASYMTX)

      REAL(FPK) :: KSQ      ( MAX_STREAMS )
      REAL(FPK) :: WK       ( MAX_2_STREAMS )
      REAL(FPK) :: EVEC     ( MAX_STREAMS,MAX_STREAMS )
      INTEGER ::          IER
      LOGICAL ::          ASYMTX_FAILURE

!  Miscellaneous local variables

      INTEGER ::          I, J, I1, L, AA, K
      REAL(FPK) :: DP, DM, SUM, FAC, KVAL
      REAL(FPK) :: NORM, XINV, HELP

!  local variables for error tracing

      CHARACTER (LEN=3) ::      CN, CI

!  initialise status
!  -----------------

      FAIL    = .FALSE.
      MESSAGE = ' '

!  Construct Eigenmatrix
!  ---------------------

!  zero the Eigenmatrix

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
          SUM = ZERO
          DO K = 1, NSTREAMS
            SUM = SUM + DAB_SAVE(I,K) * SAB_SAVE(K,J)
          ENDDO
          EIGENMAT(I,J) = SUM
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

!      DO AA = 1, NSTREAMS
!        write(*,*)evec(1:2,aa),KSQ(AA)
!      enddo
!       if ( layer.eq.23)stop'fisrt'

!  Find solution eigenvectors XPOS, XNEG for all eigenvalues
!  ---------------------------------------------------------

      DO AA = 1, NSTREAMS

!  Store positive values in output array for each layer

        KVAL = DSQRT(KSQ(AA))
        KEIGEN(AA,LAYER) = KVAL

!  Normalize eigenvectors to 1

        NORM = ZERO
        DO I = 1, NSTREAMS
          NORM = NORM + EVEC(I,AA)*EVEC(I,AA)
        ENDDO
        NORM = DSQRT(NORM)

!  Find normalized eigenvector EIGENVEC_SAVE

        DO I = 1, NSTREAMS
          EIGENVEC_SAVE(I,AA) = EVEC(I,AA)/NORM
        ENDDO

!  Find difference eigenvector DIFVEC (Siewert's notation)

        DO I = 1, NSTREAMS
          SUM = ZERO
          DO K = 1, NSTREAMS
            SUM = SUM - SAB_SAVE(I,K) * EIGENVEC_SAVE(K,AA)
          ENDDO
          DIFVEC_SAVE(I,AA) = SUM / KVAL
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
            ( N, FOURIER, NSTREAMS, NMOMENTS, N_USER_STREAMS, &
              OMEGA_PHASMOMS, XPOS, XNEG, &
              PLM_WT_PXI,  PLM_WT_MXI, PLM_PXUI, &
              U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Computational variable control

      INTEGER, INTENT(IN) ::          FOURIER, N
      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          NMOMENTS

!  Number of user streams

      INTEGER, INTENT(IN) ::          N_USER_STREAMS

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

      REAL(FPK), INTENT(INOUT) :: &
            U_XPOS ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS ), &
            U_XNEG ( MAX_USER_STREAMS, MAX_STREAMS, MAX_LAYERS )

!  Help arrays (for use with linearization)

      REAL(FPK), INTENT(OUT) :: U_HELP_P  ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: U_HELP_M  ( MAX_STREAMS, 0:MAX_MOMENTS )

!  Local variables
!  ---------------

      INTEGER ::   UM, J, J1, L, AA
      REAL(FPK) :: SUM_NEG, SUM_POS, POS1, POS2, NEG1, NEG2, ULP

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

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

      SUBROUTINE BEAMSOLUTION_MATRIX &
              ( LAYER, NSTREAMS, EIGENMAT_SAVE, &
                LOCAL_ASOURCE, LOCAL_PICUTOFF, &
                QMAT_SAVE, QPIVOT, &
                FAIL, MESSAGE )

!  This is the classical Chandrasekhar beam solution.
!  ( plane parallel or average secant only)
!  Linear Matrix algebra.

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS
      USE LRRS_AUX2

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN) ::          LAYER
      INTEGER, INTENT(IN) ::          NSTREAMS

!  Beam solution cutoff layer

      INTEGER, INTENT(IN) ::          LOCAL_PICUTOFF

!  Beam solution average secant

      REAL(FPK), INTENT(IN) :: LOCAL_ASOURCE

!  Saved matrix from eigenvalue computation

      REAL(FPK), INTENT(IN) :: EIGENMAT_SAVE ( MAX_STREAMS, MAX_STREAMS )

!  output module arguments
!  -----------------------

!  Particular solution: LU-decomposed matrix and Pivot

      REAL(FPK), INTENT(OUT) :: QMAT_SAVE ( MAX_STREAMS, MAX_STREAMS )
      INTEGER, INTENT(OUT) ::          QPIVOT    ( MAX_STREAMS )

!  output status

      LOGICAL, INTENT(OUT) ::          FAIL
      CHARACTER (LEN=70), INTENT(OUT) ::     MESSAGE

!  Local variables
!  ---------------

      INTEGER ::          I, J, INFO
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
      END SUBROUTINE BEAMSOLUTION_MATRIX

!

      SUBROUTINE BEAM_SOLUTION_MATRIX &
              ( LAYER, NSTREAMS, EIGENMAT_SAVE, &
                LOCAL_ASOURCE, LOCAL_PICUTOFF, &
                QMAT_SAVE, QPIVOT, &
                FAIL, MESSAGE )

!  This is the classical Chandrasekhar beam solution.
!  ( plane parallel or average secant only)
!  Linear Matrix algebra.

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS
      USE LRRS_AUX2

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN) ::          LAYER
      INTEGER, INTENT(IN) ::          NSTREAMS

!  Beam solution cutoff layer

      INTEGER, INTENT(IN) ::          LOCAL_PICUTOFF

!  Beam solution average secant

      REAL(FPK), INTENT(IN) :: LOCAL_ASOURCE

!  Saved matrix from eigenvalue computation

      REAL(FPK), INTENT(IN) :: EIGENMAT_SAVE ( MAX_STREAMS, MAX_STREAMS )

!  output module arguments
!  -----------------------

!  Particular solution: LU-decomposed matrix and Pivot

      REAL(FPK), INTENT(OUT) :: QMAT_SAVE ( MAX_STREAMS, MAX_STREAMS )
      INTEGER, INTENT(OUT) ::          QPIVOT    ( MAX_STREAMS )

!  output status

      LOGICAL, INTENT(OUT) ::          FAIL
      CHARACTER (LEN=70), INTENT(OUT) ::     MESSAGE

!  Local variables
!  ---------------

      INTEGER ::          I, J, INFO
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
            ( LAYER, FOURIER, NMOMENTS, NSTREAMS, &
              LOCAL_ASOURCE, LOCAL_INITRANS, &
              QMAT_SAVE, QPIVOT, SAB_SAVE, DAB_SAVE, &
              PLM_00_PXI,  PLM_00_MXI, &
              QUAD_STREAMS, OMEGA_PHASMOMS, &
              QSUMVEC_SAVE, QDIFVEC_SAVE, QVEC_SAVE, QAUX_SAVE, &
              WVECTOR, FAIL, MESSAGE )

!  This is the classical Chandrasekhar beam solution.
!  ( plane parallel or average secant only)
!       Linear Matrix algebra solution

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS
      USE LRRS_AUX2

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  layer index, Fourier number

      INTEGER, INTENT(IN) ::          LAYER
      INTEGER, INTENT(IN) ::          FOURIER

!  Given number of streams and moments

      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          NMOMENTS

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
      INTEGER, INTENT(IN) ::          QPIVOT    ( MAX_STREAMS )

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

      LOGICAL, INTENT(OUT) ::          FAIL
      CHARACTER (LEN=70), INTENT(OUT) ::     MESSAGE

!  Local variables
!  ---------------

      INTEGER ::          I, J, L, INFO, I1, NSTR2
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
                ( DO_UPWELLING, DO_DNWELLING, FOURIER, &
                  N, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, &
                  NSTREAMS, NMOMENTS, N_USER_STREAMS, &
                  PICUTOFF, WPARTIC, OMEGAMOMS, &
                  PLM_WT_PXI, PLM_WT_MXI, &
                  PLM_00_PXUI, PLM_PXUI, &
                  PLM_00_MXUI, PLM_MXUI, &
                  W_HELP, U_WPOS1, U_WPOS2, U_WNEG1, U_WNEG2 )

!  Particular integral solutions at user angles

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Directional control

      LOGICAL, INTENT(IN) ::          DO_UPWELLING
      LOGICAL, INTENT(IN) ::          DO_DNWELLING

!  Computational control

      INTEGER, INTENT(IN) ::          N
      INTEGER, INTENT(IN) ::          FOURIER
      INTEGER, INTENT(IN) ::          NSTREAMS
      INTEGER, INTENT(IN) ::          NMOMENTS

!  Layer existence

      LOGICAL, INTENT(IN) ::          STERM_LAYERMASK_UP ( MAX_LAYERS )
      LOGICAL, INTENT(IN) ::          STERM_LAYERMASK_DN ( MAX_LAYERS )

!  number of user-defined stream angles

      INTEGER, INTENT(IN) ::          N_USER_STREAMS

!  Beam solution cutoff layer

      INTEGER, INTENT(IN) ::          PICUTOFF

!  Beam solution vector

      REAL(FPK), INTENT(IN) :: WPARTIC   ( MAX_2_STREAMS, MAX_LAYERS )

!  Initial transmittance of beam to top of layer
!    NOT required now.....................................
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

      INTEGER ::          UM, J, J1, L
      REAL(FPK) :: SUM, POS1, POS2, NEG1, HELP(0:MAX_MOMENTS)

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
        SUM = ZERO
        DO  J = 1, NSTREAMS
          J1 = J + NSTREAMS
          POS1 = WPARTIC(J1,N) * PLM_WT_PXI(J,L)
          POS2 = WPARTIC(J,N)  * PLM_WT_MXI(J,L)
          SUM = SUM + POS1 + POS2
        ENDDO
        W_HELP(L) = SUM
        HELP(L)   = SUM * OMEGAMOMS(N,L)
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
         ( FOURIER, LAYER, SOLUTION, DO_MULTIPLIER, LRRS_FGSMALL, &
           NSTREAMS, QUAD_WEIGHTS, DELTAS, &
           XPOS, KEIGEN, KTRANS, EIGENNORM_SAVE, &
           PICUTOFF, FLIPPER, ASOURCE, DELTRANS, INITRANS, &
           QSOURCE_QPOS, QSOURCE_QNEG, &
           WUPPER, WLOWER, &
           ATERM_SAVE, BTERM_SAVE, MULT_M, MULT_P )

!  This is the Green's function solution, one layer one contribution

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS
      USE LRRS_SMALLNOS

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  Fourier Index

      INTEGER, INTENT(IN) ::   FOURIER

!  Given layer, solution index

      INTEGER, INTENT(IN) ::   LAYER
      INTEGER, INTENT(IN) ::   SOLUTION

!  Control for calculating Green's function multiplier

      LOGICAL, INTENT(IN) ::   DO_MULTIPLIER

!  Small number limit (Taylor series expansions)

      REAL(FPK), INTENT(IN) :: LRRS_FGSMALL

!  quadrature inputs

      INTEGER, INTENT(IN) ::   NSTREAMS
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
      REAL(FPK) :: TPA, TMA, SUM_LA, SUM_LB
      REAL(FPK) :: S_P_U, S_P_L, S_M_U, S_M_L
      REAL(FPK) :: RHO_M, RHO_P, ZDEL, ZWDEL

!  debug checking, for 37 Raman bins, 8 half-space streams. 03 June 2009
!      INTEGER ::          ss, nss, solmask(77), kfd
!      LOGICAL ::          deb
!      data nss / 77 /
!      data solmask /
!     &    1,   2,   3,  20,  21,  38,  39,  56,  57,
!     &   74,  75,  92,  93, 110, 111, 128, 129, 146, 147, 164, 165,
!     &  182, 183, 200, 201, 218, 219, 236, 237, 254, 255, 272, 273,
!     &  290, 291, 308, 309, 326, 327, 344, 345, 362, 363, 380, 381,
!     &  398, 399, 416, 417, 434, 435, 452, 453, 470, 471, 488, 489,
!     &  506, 507, 524, 525, 542, 543, 560, 561, 578, 579, 596, 597,
!     &  614, 615, 632, 633, 650, 651, 668, 669 /

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

      IF ( N .GT. PICUTOFF ) THEN
!mick fix 10/11/2016 - define ATERM_SAVE & BTERM_SAVE for this case
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
!     Decide If the small-numbers analysis applies
!         Old code   POWER = DELTA*RHO_M
!            SLIM = ZDEL * DELTA * (ONE-HALF*POWER+POWER*POWER/6.0d0)
!          if (n.eq.16)write(44,*)n,p,aa,rho_m,asource(n,p),keigen(aa,n)

          IF ( DABS(RHO_M).LT.LRRS_FGSMALL ) THEN
            CALL LIMIT_GCFUNC ( RHO_M, DELTAS, ZDEL, MULT_M(AA) )
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
!    (Initrans multiplication taken care of in Source function definitio

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
        WLEVELS(I,UPIDX) = S_P_U
        WLEVELS(I1,UPIDX) = S_M_U
        WLEVELS(I1,DNIDX) = S_M_L
        WLEVELS(I,DNIDX) = S_P_L
      ENDDO

!  Continuation point for avoiding this calculation

 4567 CONTINUE

!  Add single RRS particular integral to the RRS-summed total P.I.

      DO I = 1, 2 * NSTREAMS
        WUPPER(I,N) = WUPPER(I,N) + WLEVELS(I,UPIDX)
        WLOWER(I,N) = WLOWER(I,N) + WLEVELS(I,DNIDX)
      ENDDO

!  Note the Debug Write -> Very potent, gets all sources

!      i   = 4
!      deb = .true.
!      do ss = 1, nss
!        if ( solution.eq.solmask(ss)) deb = .false.
!      enddo
!  profiles
!      kfd = 3
!      if ((n.ge.kfd).or.((n.lt.kfd).and.deb))
!     &           write(56,'(3i5,1p2e24.14)')solution,n,i,
!     &              WLEVELS(I,UPIDX),WLEVELS(I,DNIDX)
!  column
!      kfd = 0
!      write(56,'(3i5,1p2e24.14)')solution,n,i,
!     &              WLEVELS(I,UPIDX),WLEVELS(I,DNIDX)

!  Finish

      RETURN
      END SUBROUTINE GREENFUNC_SOLUTION_1

!

      SUBROUTINE ATTENUATION_SETUP_1 &
          ( NPOINTS, NLAYERS, &
            DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS, &
            DELTAU_VERT_INPUT, CHAPMAN, &
            BEAM_PICUTOFF, BEAM_ITRANS, BEAM_AVSECANT, &
            BEAM_ETRANS, BEAM_DTRANS, SAVE_TRANS_USERM )

!  Include file
!  ------------

!  include file of dimensions and numbers

      USE LRRS_PARS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Control INTEGER ::s

      INTEGER, INTENT(IN) ::          NPOINTS, NLAYERS

!  Control flag and inputs for user streams

      LOGICAL, INTENT(IN) ::          DO_USER_STREAMS
      INTEGER, INTENT(IN) ::          N_USER_STREAMS
      REAL(FPK), INTENT(IN) :: USER_STREAMS  ( MAX_USER_STREAMS )

!  Chapman factors

      REAL(FPK), INTENT(IN) :: CHAPMAN ( MAX_LAYERS, MAX_LAYERS )

!  Basic input quantities for elastic scattering

      REAL(FPK), INTENT(IN) :: DELTAU_VERT_INPUT   ( MAX_LAYERS, MAX_POINTS )

!  Output arguments
!  ----------------

!  Solar beam transmittances, average secant factors

      INTEGER, INTENT(OUT) ::          BEAM_PICUTOFF ( MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: BEAM_ITRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: BEAM_AVSECANT ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: BEAM_ETRANS   ( MAX_LAYERS, MAX_POINTS )
      REAL(FPK), INTENT(OUT) :: BEAM_DTRANS   ( MAX_LAYERS, MAX_POINTS )

!  user stream transmittances, whole layers

      REAL(FPK), INTENT(OUT) :: SAVE_TRANS_USERM &
                ( MAX_USER_STREAMS, MAX_LAYERS, MAX_POINTS )

!  Local variables
!  ---------------

      INTEGER ::          N, K, S, UM, CUTOFF
      REAL(FPK) :: SPHER, SUM, LOCAL_LTRANS, TAU_SLANT, ST0, ST1
      REAL(FPK) :: TAUSLANT(0:MAX_LAYERS), DELTA

!  Beam attenuation
!  ----------------

!  Monochromatic - 234 points (233 shifts + 1 calculation point)
!  Binned        - all points in outer buffer

!      IF ( DO_BIN_REALIZATION ) THEN
!        NPOINTS_LOCAL = NPOINTS_OUTER
!      ELSE
!        NPOINTS_LOCAL = NPOINTS_MONO
!      ENDIF

!  start loop over all outer points

      DO S = 1, NPOINTS

!  Tau slant = Beam slant optical thickness

        TAUSLANT(0) = ZERO
        DO N = 1, NLAYERS
          SUM = ZERO
          DO K = 1, N
            SUM = SUM + DELTAU_VERT_INPUT(K,S) * CHAPMAN(N,K)
          ENDDO
          TAUSLANT(N) = SUM
        ENDDO

!  Initial transmittances, and average secants
!   This section lifted from the LIDORT code

        ST0    = ONE
        ST1    = ST0
        CUTOFF = NLAYERS
        DO N = 1, NLAYERS
          DELTA = DELTAU_VERT_INPUT(N,S)
          IF ( N.LE.CUTOFF) THEN
            IF ( TAUSLANT(N) .GT. MAX_TAU_SPATH ) THEN
              CUTOFF = N
            ELSE
              ST1 = DEXP ( - TAUSLANT(N) )
            ENDIF
            BEAM_AVSECANT(N,S) = (TAUSLANT(N)-TAUSLANT(N-1))/DELTA
            BEAM_ITRANS(N,S)   = ST0
            ST0                = ST1
          ELSE
            BEAM_AVSECANT(N,S) = ZERO
            BEAM_ITRANS(N,S)   = ZERO
          ENDIF
        ENDDO
        BEAM_PICUTOFF(S) = CUTOFF

!  Bottom of layer transmittances, ETRANS

        DO N = 1, NLAYERS
          IF ( N.GT.BEAM_PICUTOFF(S) ) THEN
            LOCAL_LTRANS = ZERO
          ELSE
            TAU_SLANT    = DELTAU_VERT_INPUT(N,S) * BEAM_AVSECANT(N,S)
            LOCAL_LTRANS = DEXP ( - TAU_SLANT )
          ENDIF
          BEAM_ETRANS(N,S) = BEAM_ITRANS(N,S) * LOCAL_LTRANS
          BEAM_DTRANS(N,S) = LOCAL_LTRANS
        ENDDO

!  End loop over points

      ENDDO

!  whole layer transmittances along user streams
!  =============================================

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO S = 1, NPOINTS
            DO N = 1, NLAYERS
              SPHER = DELTAU_VERT_INPUT(N,S) / USER_STREAMS(UM)
              IF ( SPHER.GT.MAX_TAU_UPATH ) THEN
                SAVE_TRANS_USERM(UM,N,S) = ZERO
              ELSE
                SAVE_TRANS_USERM(UM,N,S) = DEXP ( - SPHER )
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE ATTENUATION_SETUP_1

!

      SUBROUTINE LEGENDRE_SETUPS &
            ( FOURIER,  NSTREAMS, NMOMENTS, &
              DO_USER_STREAMS, N_USER_STREAMS, USER_STREAMS, &
              COS_SZA, QUAD_STREAMS, QUAD_WEIGHTS, &
              PLM_PXI,     PLM_MXI, &
              PLM_PXUI,    PLM_MXUI, &
              PLM_00_PXI,  PLM_00_MXI, &
              PLM_00_PXUI, PLM_00_MXUI, &
              PLM_WT_PXI,  PLM_WT_MXI )

!  Legendre polynomials and assocaited quantities

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS
      USE LRRS_AUX2

      IMPLICIT NONE

!  input module arguments
!  ----------------------

!  Fourier number, number of streams and moments

      INTEGER, INTENT(IN) ::             FOURIER
      INTEGER, INTENT(IN) ::             NSTREAMS
      INTEGER, INTENT(IN) ::             NMOMENTS

!  Number of user streams and flag

      INTEGER, INTENT(IN) ::             N_USER_STREAMS
      LOGICAL, INTENT(IN) ::             DO_USER_STREAMS

!  Cosine of the solar zenith angle

      REAL(FPK), INTENT(IN) :: COS_SZA

!  quadrature

      REAL(FPK), INTENT(IN) :: QUAD_STREAMS ( MAX_STREAMS )
      REAL(FPK), INTENT(IN) :: QUAD_WEIGHTS ( MAX_STREAMS )

!  user stream directions

      REAL(FPK), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )

!  output arguments
!  ----------------

      REAL(FPK), INTENT(OUT) :: PLM_PXI    ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_MXI    ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_00_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_00_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_WT_PXI ( MAX_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_WT_MXI ( MAX_STREAMS, 0:MAX_MOMENTS )

      REAL(FPK), INTENT(OUT) :: PLM_PXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_MXUI    ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_00_PXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )
      REAL(FPK), INTENT(OUT) :: PLM_00_MXUI ( MAX_USER_STREAMS, 0:MAX_MOMENTS )

!  local variables

      REAL(FPK) :: CFPLG  ( 0:MAX_MOMENTS )
      REAL(FPK) :: PLM_00 ( 0:MAX_MOMENTS )
      INTEGER ::          I, UI, L, LPM

!  Legendre polynomials for quadrature streams
!  -------------------------------------------

!  .. positive computational angle streams
!  .. negative streams by symmetry relations

      DO I = 1, NSTREAMS
        CALL CFPLGARR &
            ( MAX_MOMENTS, NMOMENTS, FOURIER, QUAD_STREAMS(I), &
              CFPLG)
        DO L = FOURIER, NMOMENTS
          PLM_PXI(I,L) = CFPLG(L)
        ENDDO
      ENDDO
      DO L = FOURIER, NMOMENTS
        LPM = L + FOURIER
        IF (MOD(LPM,2).EQ.0) THEN
          DO I = 1, NSTREAMS
            PLM_MXI(I,L) = PLM_PXI(I,L)
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            PLM_MXI(I,L) = - PLM_PXI(I,L)
          ENDDO
        ENDIF
      ENDDO

!  Legendre polynomials for User-angle streams
!  -------------------------------------------

      IF ( DO_USER_STREAMS ) THEN
        DO UI = 1, N_USER_STREAMS
          CALL CFPLGARR &
            ( MAX_MOMENTS, NMOMENTS, FOURIER, USER_STREAMS(UI), &
              CFPLG)
          DO L = FOURIER, NMOMENTS
            PLM_PXUI(UI,L) = CFPLG(L)
          ENDDO
        ENDDO
        DO L = FOURIER, NMOMENTS
          LPM = L + FOURIER
          IF (MOD(LPM,2).EQ.0) THEN
            DO UI = 1, N_USER_STREAMS
              PLM_MXUI(UI,L) = PLM_PXUI(UI,L)
            ENDDO
          ELSE
            DO UI = 1, N_USER_STREAMS
              PLM_MXUI(UI,L) = - PLM_PXUI(UI,L)
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Legendre polynomials for solar streams
!  --------------------------------------

!  .. negative solar zenith angle cosine

      CALL CFPLGARR &
           ( MAX_MOMENTS, NMOMENTS, FOURIER, COS_SZA, &
             PLM_00 )

!  Legendre polynomial products
!  ----------------------------

      DO L = FOURIER, NMOMENTS
        DO I = 1, NSTREAMS
          PLM_00_PXI(I,L) = PLM_00(L) * PLM_PXI(I,L)
          PLM_00_MXI(I,L) = PLM_00(L) * PLM_MXI(I,L)
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            PLM_00_PXUI(UI,L) = PLM_00(L) * PLM_MXUI(UI,L)
            PLM_00_MXUI(UI,L) = PLM_00(L) * PLM_PXUI(UI,L)
          ENDDO
        ENDIF
        DO I = 1, NSTREAMS
          PLM_WT_PXI(I,L) = QUAD_WEIGHTS(I) * PLM_PXI(I,L) * HALF
          PLM_WT_MXI(I,L) = QUAD_WEIGHTS(I) * PLM_MXI(I,L) * HALF
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LEGENDRE_SETUPS

!

      SUBROUTINE DBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, &
             NSTREAMS, FL1, FL2, BIREFLEC_0, &
             COS_SZA, FOURIER, DELTA_FACTOR, BEAM_ETRANS_OPDEP, &
             ATMOS_ATTN, DIRECT_BEAM )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Inclusion flag

      LOGICAL, INTENT(IN) ::          DO_INCLUDE_SURFACE

!  Lambertian flag

      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_SURFACE

!  Computational variable

      INTEGER, INTENT(IN) ::          NSTREAMS

!  fractions

      REAL(FPK), INTENT(IN) :: FL1, FL2

!  BRDF inputs

      REAL(FPK), INTENT(IN) :: BIREFLEC_0      ( MAX_STREAMS )

!  Solar zenith angle cosine

      REAL(FPK), INTENT(IN) :: COS_SZA

!  Fourier number

      INTEGER, INTENT(IN) ::          FOURIER

!  2-deltam0 factor

      REAL(FPK), INTENT(IN) :: DELTA_FACTOR

!  Optical thickness value

      REAL(FPK), INTENT(IN) :: BEAM_ETRANS_OPDEP

!  output Direct beam module
!  -------------------------

!  Direct beam itself

      REAL(FPK), INTENT(OUT) :: DIRECT_BEAM      ( MAX_STREAMS )

!  Attenuation

      REAL(FPK), INTENT(OUT) :: ATMOS_ATTN

!  Local variables
!  ---------------

      REAL(FPK) :: X0_FLUX, REFL_ATTN, FRAC1, FRAC2
      INTEGER ::          I

!  Initialize
!  ----------

!   Safety first!  Return if there is no reflection.

      DO I = 1, NSTREAMS
        DIRECT_BEAM(I) = ZERO
      ENDDO
      ATMOS_ATTN = ZERO                       ! JvG 4/28/10
      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Attenuation of solar beam
!  -------------------------

!  Bug fixed 23 November 2005, by R. Spurr
!    ( First noticed for the LIDORT codes by Italians!)

!     Old code only worked for FLUX_FACTOR = 1
!      X0_FLUX    = FOUR * COS_SZA * FLUX_FACTOR / DELTA_FACTOR
!    New code: Drop the FLUX_FACTOR in the next statement:

      X0_FLUX    = FOUR * COS_SZA / DELTA_FACTOR
      ATMOS_ATTN = X0_FLUX * BEAM_ETRANS_OPDEP

!  For the Lambertian or BRDF case
!    - reflected along the quadrature directions
!    - reflected along the user-defined directions

!   Fractional components are FL1 and FL2
!     BRDF case may have Lambertian fraction
!     BRDF case may have waterleaving option

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        REFL_ATTN = FL1 * ATMOS_ATTN
        IF ( FOURIER .EQ. 0 ) THEN
          DO I = 1, NSTREAMS
            DIRECT_BEAM(I) = REFL_ATTN
          ENDDO
        ENDIF
      ELSE
        FRAC1 = FL1 * ATMOS_ATTN
        FRAC2 = FL2 * ATMOS_ATTN
        DO I = 1, NSTREAMS
          DIRECT_BEAM(I) = FRAC1 + FRAC2 * BIREFLEC_0(I)
        ENDDO
      ENDIF

!  end direct beam calculation

      RETURN
      END SUBROUTINE DBEAM_SETUP

!

      SUBROUTINE UDBEAM_SETUP &
           ( DO_INCLUDE_SURFACE, DO_SS_CORRECTION, &
             DO_LAMBERTIAN_SURFACE, FL1, FL2, USER_BIREFLEC_0, &
             FOURIER, N_USER_STREAMS, ATMOS_ATTN, &
             USER_DIRECT_BEAM )

!  include file of dimensions and numbers
!  --------------------------------------

      USE LRRS_PARS

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Inclusion flag

      LOGICAL, INTENT(IN) ::          DO_INCLUDE_SURFACE

!  SS correction flag (formerly DB correction flag)
!   SS now incorporates the exact direct bounce (DB term)

      LOGICAL, INTENT(IN) ::          DO_SS_CORRECTION

!  Lambertian flag

      LOGICAL, INTENT(IN) ::          DO_LAMBERTIAN_SURFACE

!  fractions

      REAL(FPK), INTENT(IN) :: FL1, FL2

!  BRDF inputs

      REAL(FPK), INTENT(IN) :: USER_BIREFLEC_0 ( MAX_USER_STREAMS )

!  Fourier number

      INTEGER, INTENT(IN) ::          FOURIER

!  User stream angles input

      INTEGER, INTENT(IN) ::          N_USER_STREAMS

!  Attenuation input

      REAL(FPK), INTENT(IN) :: ATMOS_ATTN

!  output Direct beam module
!  -------------------------

      REAL(FPK), INTENT(OUT) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS )

!  Local variables
!  ---------------

      REAL(FPK) :: REFL_ATTN, FRAC1, FRAC2
      INTEGER ::          UI

!  Initialize
!  ----------

!   Safety first!  Return if there is no reflection.

      DO UI = 1, N_USER_STREAMS
        USER_DIRECT_BEAM(UI) = ZERO
      ENDDO
      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!   Fractional components are FL1 and FL2
!     BRDF case may have Lambertian fraction
!     BRDF case may have waterleaving option

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        REFL_ATTN = FL1 * ATMOS_ATTN
        IF ( FOURIER .EQ. 0 ) THEN
          IF ( .NOT.DO_SS_CORRECTION ) THEN
           DO UI = 1, N_USER_STREAMS
              USER_DIRECT_BEAM(UI) = REFL_ATTN
            ENDDO
          ENDIF
        ENDIF
      ELSE
        FRAC1 = FL1 * ATMOS_ATTN
        FRAC2 = FL2 * ATMOS_ATTN
        IF ( .NOT.DO_SS_CORRECTION ) THEN
          DO UI = 1, N_USER_STREAMS
            USER_DIRECT_BEAM(UI) = FRAC1 + FRAC2*USER_BIREFLEC_0(UI)
          ENDDO
        ENDIF
      ENDIF

!  end direct beam calculation

      RETURN
      END SUBROUTINE UDBEAM_SETUP

!

      SUBROUTINE CHAPMAN_FUNCTION &
           ( DO_PLANE_PARALLEL, NLAYERS, &
             COS_SZA, EARTH_RADIUS, HEIGHT_GRID, &
             CHAPMAN_FACTORS )

!  This is the Chapman function calculation of the slant path
!  geometrical factors, defined as the ratios between the slant path
!  distances to the corresponding vertical distances.

!  You must specify the Earth_radius and the height grid in order
!  to make this work.

!  This is a straightforward geometrical calculation, and is only
!  valid for a NON-REFRACTIVE atmosphere.

!  include file of dimensions and numbers

      USE LRRS_PARS

      IMPLICIT NONE

!  Input arguemnts
!  ---------------

      LOGICAL, INTENT(IN) ::          DO_PLANE_PARALLEL
      INTEGER, INTENT(IN) ::          NLAYERS
      REAL(FPK), INTENT(IN) :: COS_SZA
      REAL(FPK), INTENT(IN) :: EARTH_RADIUS
      REAL(FPK), INTENT(IN) :: HEIGHT_GRID(0:MAX_LAYERS)

!  output arguemnts
!  ----------------

      REAL(FPK), INTENT(OUT) :: CHAPMAN_FACTORS(MAX_LAYERS,MAX_LAYERS)

!  Local variables
!  ---------------

      INTEGER ::          N, M
      REAL(FPK) :: GM_TOA, HELP1, HELP2
      REAL(FPK) :: H(0:MAX_LAYERS), DELZ(MAX_LAYERS)
      REAL(FPK) :: STH, CTH, DELS, S1, S0

!  get spherical optical depths
!  ----------------------------

!  Prepare spherical attenuation (shell geometry)

      IF ( .NOT.DO_PLANE_PARALLEL ) THEN

        GM_TOA = DSQRT ( 1.0D0 - COS_SZA * COS_SZA )
        DO N = 0, NLAYERS
          H(N) = HEIGHT_GRID(N) + EARTH_RADIUS
        ENDDO
        DO N = 1, NLAYERS
          DELZ(N) = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        ENDDO

        DO N = 1, NLAYERS
          STH = GM_TOA * H(N)/H(0)
          CTH = DSQRT ( ONE - STH * STH )
          S0 = ZERO
          HELP1 = H(0)*CTH
          HELP2 = -H(0)*H(0)*STH*STH
          DO M = 1, N
            S1 = HELP1 - DSQRT(HELP2 + H(M)*H(M))
            DELS = S1 - S0
            CHAPMAN_FACTORS(N,M) = DELS / DELZ(M)
            S0 = S1
          ENDDO
        ENDDO

!  Plane parallel

      ELSE

        DO N = 1, NLAYERS
          DO M = 1, N
            CHAPMAN_FACTORS(N,M) = ONE / COS_SZA
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE CHAPMAN_FUNCTION


      END MODULE lrrs_rtsolutions_1

