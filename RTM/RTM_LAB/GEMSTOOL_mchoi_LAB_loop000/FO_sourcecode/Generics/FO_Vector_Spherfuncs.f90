module FO_Vector_spherfuncs_m

public

contains

SUBROUTINE FO_Vector_spherfuncs &
        ( STARTER, MAXMOMENTS, NMOMENTS, NSTOKES, DTR, & ! Inputs
          DO_SUNLIGHT, THETA, ALPHA, PHI, MUIN, VSIGN, & ! Inputs
          ROTATIONS, GSHELP, GENSPHER )                  ! Outputs

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs. [Starter flag may be re-set]

      LOGICAL  , intent(inout) :: STARTER
      LOGICAL  , intent(in)    :: DO_SUNLIGHT
      INTEGER  , intent(in)    :: MAXMOMENTS, NMOMENTS, NSTOKES
      REAL(fpk), intent(in)    :: DTR, THETA, ALPHA, PHI, MUIN, VSIGN

!  Outputs

      REAL(fpk), intent(inout) :: ROTATIONS(4)
      REAL(fpk), intent(inout) :: GENSPHER(0:MAXMOMENTS,4)
      REAL(fpk), intent(inout) :: GSHELP(7,MAXMOMENTS)

!  Local

      integer   :: L, LNEW, LOLD, ITMP
      REAL(fpk) :: Dl, Dl1, Dl21, DL41, DLL1
      REAL(fpk) :: SQL4, SQL41, QROOT6, PERLL4, WFACT
      REAL(fpk) :: CTHETA, STHETA, CALPHA, SALPHA, CPHI
      REAL(fpk) :: THETA_R, ALPHA_R, PHI_R, C1, S1, C2, S2
      REAL(fpk) :: SINSCAT, HELP_SINSCAT, MUP1, MUM1, MUS
      REAL(fpk) :: CSIG1, CSIG2, SSIG1, SSIG2, CSIG1_2, CSIG2_2
      REAL(fpk) :: P00(2), P02(2), P2P2(2), P2M2(2)

!  Parameter numbers

      REAL(fpk), parameter :: ZERO  = 0.0_fpk
      REAL(fpk), parameter :: ONE   = 1.0_fpk
      REAL(fpk), parameter :: TWO   = 2.0_fpk
      REAL(fpk), parameter :: FOUR  = 4.0_fpk
      REAL(fpk), parameter :: SMALL = 1.0e-12_fpk

!  Section 1: Generate help array, only do this once

      IF ( STARTER ) THEN
         GSHELP = ZERO ; GENSPHER = ZERO ; ROTATIONS = ZERO
         DO L = 2, NMOMENTS
            DL = DBLE(L) ; DL1 = DBLE(L-1) ; DL21 = DBLE(2*L-1)
            GSHELP(1,L) = DL21 / DL
            GSHELP(2,L) = DL1  / DL
            if ( nstokes .gt. 1 ) then
               if ( L .EQ. 2 ) THEN
                  SQL41 = ZERO
               ELSE
                  SQL4  = SQL41 ; DL41  = DL*DL - FOUR
                  SQL41 = DSQRT(DL41)
                  GSHELP(3,L)  = DL21 / SQL41
                  GSHELP(4,L)  = SQL4 / SQL41
                  if ( .not. do_sunlight ) then
                     DLL1 = DL * DL1  ;  PERLL4 = ONE / DL1 / DL41
                     GSHELP(5,L) = DL21 * DLL1 * PERLL4
                     GSHELP(6,L) = DL21 * FOUR * PERLL4
                     GSHELP(7,L) = DL  * ( DL1*DL1 - FOUR) * PERLL4
                  endif
               ENDIF
            endif
         ENDDO
         STARTER = .false.
      ENDIF

!  Section 2. Rotation angles
!  --------------------------

     mus = muin

     THETA_R = DTR * THETA ; CTHETA = COS(THETA_R) ; STHETA = SIN(THETA_R)
     ALPHA_R = DTR * ALPHA ; CALPHA = COS(ALPHA_R) ; SALPHA = SIN(ALPHA_R)
     PHI_R   = DTR * PHI   ; CPHI   = COS(PHI_R)

!  THIS IS VALID ONLY FOR NON-REFRACTING ATMOSPHER
!  VSIGN = -1 FOR UPWELLING, +1 FOR DOWNWELLING
!  COSINE SIGMA 1 AND 2. H/VDM, EQS. (99)-(101)

!  A. SAFETY
!    WATCH FOR SIN^2(SCATTER ANGLE) LESS THAN ZERO (MACHINE PRECISION)
!    R. SPURR, 16 JANUARY 2006, RT SOLUTIONS INC.

      HELP_SINSCAT = ( ONE - MUS * MUS )
      IF ( HELP_SINSCAT.LE.ZERO ) THEN
        SINSCAT = SMALL
      ELSE
        SINSCAT = SQRT ( HELP_SINSCAT )
      ENDIF

!  B. NECESSARY LIMIT ANALYSES - HOVENIER LIMITS.
!     R. SPURR AND V. NATRAJ, 17 JANUARY 2006

      IF ( DABS(SINSCAT) .LE. SMALL ) THEN
        CSIG1 = ZERO
        CSIG2 = ZERO
      ELSE
        IF ( STHETA .EQ. ZERO ) THEN
          CSIG1 = -  CPHI
        ELSE
          CSIG1 = ( + VSIGN * CALPHA + CTHETA * MUS ) / SINSCAT / STHETA
!          CSIG1 = ( + CALPHA + CTHETA * MUS ) / SINSCAT / STHETA
        ENDIF
        IF ( SALPHA .EQ. ZERO ) THEN
          CSIG2 = -  CPHI
        ELSE
          CSIG2 = ( - CTHETA - VSIGN * CALPHA * MUS ) / SINSCAT / SALPHA
!          CSIG2 = ( - CTHETA - CALPHA * MUS ) / SINSCAT / SALPHA
        ENDIF
      ENDIF
!      write(*,*)'hello',CSIG1,CSIG2

!  DEBUG, 06 OCTOBER 2010. FOUND DOWNWELLING BUG, NEED SEPARATE GEOMETRY
!      WRITE(7,'(I4,2F6.1,1P6E15.7)')N,VSIGN,PHI/DEG_TO_RAD,
!     &      SALPHA, CALPHA, STHETA, CTHETA, CSIG1,CSIG2

!  THESE LINES ARE NECESSARY TO AVOID BAD VALUES

      IF ( CSIG1 .GT. ONE  ) CSIG1 = ONE
      IF ( CSIG1 .LT. -ONE ) CSIG1 = -ONE
      IF ( CSIG2 .GT. ONE  ) CSIG2 = ONE
      IF ( CSIG2 .LT. -ONE ) CSIG2 = -ONE

!  OUTPUT, H/VDM, EQS. (89)-(94)

      CSIG1_2 = TWO * CSIG1
      CSIG2_2 = TWO * CSIG2
      IF ( ABS(CSIG1-ONE) .LT. SMALL ) THEN
        SSIG1 = ZERO
      ELSE
        SSIG1 = SQRT ( ONE - CSIG1 * CSIG1 )
      ENDIF
      IF ( ABS(CSIG2-ONE).LT. SMALL ) THEN
        SSIG2 = ZERO
      ELSE
        SSIG2 = SQRT ( ONE - CSIG2 * CSIG2 )
      ENDIF

!  RELAZM IN [180,360), NEED SIGN REVERSAL FOR S1/S2. H/VDM, EQS. 94-95.

      C1 = CSIG1_2 * CSIG1 - ONE
      C2 = CSIG2_2 * CSIG2 - ONE

!  PHI IS INPUT IN DEGREES

      IF ( PHI .LE. 180.0_fpk ) THEN
        S1 = CSIG1_2 * SSIG1
        S2 = CSIG2_2 * SSIG2
      ELSE
        S1 = -CSIG1_2 * SSIG1
        S2 = -CSIG2_2 * SSIG2
      ENDIF

!  OUTPUT

      ROTATIONS(1) = C1
      ROTATIONS(2) = S1
      ROTATIONS(3) = C2
      ROTATIONS(4) = S2

!  Section 3a. Legendre only: Nstokes = 1 (scalar case). Return when done
!  -----------------------------------------------------

      IF ( NSTOKES .eq. 1 ) THEN
         GENSPHER(0,1) = ONE
         GENSPHER(1,1) = MUS
         DO L = 2, NMOMENTS
            GENSPHER(L,1) = GSHELP(1,L)*GENSPHER(L-1,1)*MUS - GSHELP(2,L)*GENSPHER(L-2,1)
         ENDDO
         RETURN
      ENDIF

!  Section 3b. : Nstokes > 1 (vector case). General case
!  -----------------------------------------------------

!  START LOOP OVER THE COEFFICIENT INDEX L
!  FIRST UPDATE GENERALIZED SPHERICAL FUNCTIONS, THEN CALCULATE COEFS.
!  LOLD AND LNEW ARE POINTER-LIKE INDICES USED IN RECURRENCE

      LNEW = 1
      LOLD = 2
      QROOT6 = - 0.25_fpk * SQRT(6.0_fpk)

      DO L = 0, NMOMENTS

!  FIRST MOMENT

        IF ( L .EQ. 0 ) THEN
          P00(LOLD)  = ONE  ; P00(LNEW)  = ZERO
          P02(LOLD)  = ZERO ; P02(LNEW)  = ZERO
          P2P2(LOLD) = ZERO ; P2P2(LNEW) = ZERO
          P2M2(LOLD) = ZERO ; P2M2(LNEW) = ZERO
        ELSE IF ( L.eq.1 ) then
          P00(LOLD) = MUS
        ELSE
          P00(LOLD) = GSHELP(1,L) * MUS * P00(LNEW) - GSHELP(2,L)*P00(LOLD)
        END IF
        GENSPHER(L,1) = P00(LOLD)

!  Other moments

        IF ( L .EQ. 2 ) THEN

! ADDING PAPER EQ. (78)
! SQL4 CONTAINS THE FACTOR DSQRT((L+1)*(L+1)-4) NEEDED IN
! THE RECURRENCE EQS. (81) AND (82)

          MUM1 = ONE - MUS
          MUP1 = ONE + MUS
          P02(LOLD) = QROOT6 * MUM1 * MUP1
          P02(LNEW) = ZERO
          SQL41     = ZERO
          GENSPHER(L,2) = P02(LOLD)

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L = 2

          if ( .not. do_sunlight ) then
             P2P2(LOLD)= MUP1 * MUP1 / FOUR
             P2M2(LOLD)= MUM1 * MUM1 / FOUR
             GENSPHER(L,3) = P2P2(LOLD)
             GENSPHER(L,4) = P2M2(LOLD)
          endif

        ELSE IF ( L .GT. 2) THEN

! ADDING PAPER EQ. (82) WITH M=0

          P02(LOLD) = GSHELP(3,L) * MUS * P02(LNEW) - GSHELP(4,L) * P02(LOLD)
          GENSPHER(L,2) = P02(LOLD)

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L > 2

          IF ( .not. do_sunlight ) then
             WFACT = GSHELP(5,L) * MUS - GSHELP(6,L)
             P2P2(LOLD) = WFACT * P2P2(LNEW) - GSHELP(7,L) * P2P2(LOLD)
             WFACT = GSHELP(5,L) * MUS + GSHELP(6,L)
             P2M2(LOLD) = WFACT * P2M2(LNEW) - GSHELP(7,L) * P2M2(LOLD)
             GENSPHER(L,3) = P2P2(LOLD)
             GENSPHER(L,4) = P2M2(LOLD)
           ENDIF

        END IF

! SWITCH INDICES SO THAT LNEW INDICATES THE FUNCTION WITH
! THE PRESENT INDEX VALUE L; PREVENTS SWAPPING OF ENTIRE ARRAYS.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

!  END MOMENTS

      END DO

!  Finish

      RETURN
END SUBROUTINE FO_Vector_spherfuncs

!  Finish

end module FO_Vector_spherfuncs_m

