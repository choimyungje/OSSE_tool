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
! #  This is brdf_sup_aux.f90.  Utility routines.               #
! #  The subroutines are listed below with their                #
! #  source of origin (order of appearance).                    #
! #                                                             #
! #      GETQUAD2                    M. Christi, 2017           #
! #      DERFC_E                     V. Natraj, 2005            #
! #      BRDF_Fresnel_Complex        R. Spurr, 2014 (Vers. 3.7) #
! #      BRDF_QUADRATURE_Gaussian    R. Spurr, 2004             #
! #      BRDF_QUADRATURE_Trapezoid*  R. Spurr, 2004 (not used)  #
! #                                                             #
! ###############################################################

      MODULE brdf_sup_aux_m

!  Developed for LRRS Version 2.5, 9/8/15. R. Spurr, RT SOLUTIONS Inc.
!   Closely follows the LIDORT modul with the same name, but
!     (1) No surface emission, solar angle dimensioning, No observational geometry.
!     (2) Additional wavelength dimensioning (MAX_BRDF_POINTS)
!     (3) Additional control for using 1 or all wavelengths

      USE LRRS_pars_m, only : FPK, ZERO, ONE, TWO, HALF, QUARTER, PIE, &
                              MAXSTREAMS_BRDF, MAXSTHALF_BRDF

!  Everything public here

      private
      public :: GETQUAD2, &
                DERFC_E, BRDF_Fresnel_Complex, &
                BRDF_QUADRATURE_Gaussian, &
                BRDF_QUADRATURE_Trapezoid

      CONTAINS

      SUBROUTINE GETQUAD2(A,B,N,ROOTS,WGTS)

!  Computes N roots and weights for Gauss-Legendre quadrature on the interval (a,b)

      IMPLICIT NONE

!  Limits of interval

      REAL(FPK), INTENT(IN)  :: A, B

!  Dimension

      INTEGER, INTENT(IN) :: N

!  Quadrature roots and weights

      REAL(FPK), INTENT(OUT) :: ROOTS(N), WGTS(N)

!  Local variables

      INTEGER   :: I, M, N2, NM1
      REAL(FPK) :: IR, MR, NR
      REAL(FPK) :: MIDPT, SFAC
      REAL(FPK) :: DLP_DX, LP, LPM1, LPM2, X, XOLD, XX

!  Threshold for Newton's Method

      REAL(FPK), PARAMETER :: QEPS = 1.0D-13

!  Since roots are symmetric about zero on the interval (-1,1), split the interval
!  in half and only work on the lower half of the interval (-1,0).

      N2 = INT((N + 1)/2)
      NR = REAL(N,FPK)

!  Define the shift [midpoint of (a,b)] and scale factor to later move roots from
!  the interval (-1,1) to the interval (a,b)

      MIDPT = HALF*(B + A)
      SFAC  = HALF*(B - A)

      DO M = 1, N2

!  Find current root of the related Nth order Legendre Polynomial on (-1,0) by Newton's
!  Method using two Legendre Polynomial recurrence relations (e.g. see Abramowitz &
!  Stegan (1972))

         !Define starting point [ after Tricomi (1950) ]
         MR = REAL(M,FPK)
         XX = PIE*(MR - QUARTER)/(NR + HALF)
         X  = (ONE - (NR - ONE)/(8.0_FPK*NR**3) &
             - ONE/(384.0_FPK*NR**4)*(39.0_FPK - 28.0_FPK/SIN(XX)**2))*COS(XX)

         !Use Newton's Method
         DO 
            LPM1 = ZERO ; LP = ONE
            DO I = 1, N
               IR = REAL(I,FPK) ; LPM2 = LPM1 ; LPM1 = LP
               LP = ((TWO*IR - ONE)*X*LPM1 - (IR - ONE)*LPM2)/IR
            ENDDO
            DLP_DX = NR*(X*LP - LPM1)/(X**2 - ONE)
            XOLD = X ; X = XOLD - LP/DLP_DX
            IF (ABS(X-XOLD) <= QEPS) EXIT
         ENDDO

!  Shift and scale the current root (and its symmetric counterpart) from the interval (-1,1)
!  to the interval (a,b).  Define their related weights (e.g. see Abramowitz & Stegan (1972)).
!  Note:
!  If (1) N is even or (2) N is odd and M /= N2, then ROOTS(M) and ROOTS(NM1) are unique.
!  If N is odd and M = N2, then M = NM1 and ROOTS(M) = ROOTS(NM1) are one and the same root.

         !On interval lower half: (a,midpt)
         ROOTS(M)   = MIDPT - SFAC*X
         WGTS(M)    = (TWO*SFAC)/((ONE - X**2)*DLP_DX**2)

         !On interval upper half: (midpt,b)
         NM1 = N - M + 1
         ROOTS(NM1) = MIDPT + SFAC*X
         WGTS (NM1) = WGTS(M)

      ENDDO

      END SUBROUTINE GETQUAD2

!

      REAL(fpk) function derfc_e(x)

      IMPLICIT NONE

      REAL(fpk) :: x

! Returns the complementary error function erfc(x) with fractional error
! everywhere less than 1.2 * 10^7.

      REAL(fpk) :: t,z

      z = dabs(x)
      t = 1.d0/(1.d0+0.5d0*z)
      derfc_e = t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+ &
              t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t* &
              (-1.13520398d0+t*(1.48851587d0+t*(-.82215223d0+t* &
              .17087277d0)))))))))
      if (x .lt. 0.d0) derfc_e = 2.d0-derfc_e

      return
      END function derfc_e

!

Subroutine BRDF_Fresnel_Complex ( MR, MI, COSCHI, FP )

!  Renamed for the BRDF supplement.
!    (Same routine occurs also in the SLEAVE suite)
!    Adapated from SixS code, this is essentially Born/Wolf computation

   implicit none
  
!  Arguments

   real(fpk), intent(in)  :: MR, MI, COSCHI
   real(fpk), intent(out) :: FP

!  Local

   real(fpk) :: MRSQ, MISQ, MSQ, MRMI2, SINCHI_SQ, AA, A1, A2, B1, B2
   real(fpk) :: U, V, VSQ, CMU, CPU, RR2
   real(fpk) :: B1MU, B1PU, B2MV, B2PV, RL2

!  Calculation of FP, Complex RI

   IF ( MI.eq.zero) goto 67

   MRSQ = MR * MR ; MISQ = MI * MI
   MSQ   = MRSQ - MISQ
   MRMI2 = two * MR * MI

   SINCHI_SQ = one - COSCHI * COSCHI 
   AA = MSQ - SINCHI_SQ
   A1 = abs(AA)
   A2 = SQRT ( AA*AA + MRMI2 * MRMI2 )

   U = sqrt(half*abs(A1+A2))
   V = sqrt(half*abs(-A1+A2))
   VSQ = V * V
   CMU = ( COSCHI - U ) ; CPU = ( COSCHI + U )
   RR2 = ( CMU*CMU + VSQ ) / ( CPU*CPU + VSQ )

   B1 = MSQ   * COSCHI
   B2 = MRMI2 * COSCHI
   B1MU = B1 - U ; B1PU = B1 + U 
   B2PV = B2 + V ; B2MV = B2 - V 

   RL2 = ( B1MU*B1MU + B2PV*B2PV ) / ( B1PU*B1PU + B2MV*B2MV )
   FP = half * ( RR2 + RL2 )
   return

!  Calculation of FP. Real RI

67 continue
   MSQ = MR * MR
   SINCHI_SQ = one - COSCHI * COSCHI 
   U = sqrt(abs(MSQ - SINCHI_SQ))
   CMU = ( COSCHI - U ) ; CPU = ( COSCHI + U )
   RR2  = CMU*CMU / ( CPU*CPU )
   B1   = MSQ * COSCHI
   B1MU = B1 - U ; B1PU = B1 + U 
   RL2  = B1MU*B1MU / ( B1PU*B1PU )
   FP   = half * ( RR2 + RL2 )

!  Finish

   return
end subroutine BRDF_Fresnel_Complex


!

SUBROUTINE BRDF_QUADRATURE_GAUSSIAN &
        ( NSTREAMS_BRDF, NBRDF_HALF,         & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF )   ! Outputs

      IMPLICIT NONE

!  Input
!  =====


!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS_BRDF, NBRDF_HALF

!  OUTPUT
!  ======

!  azimuth quadrature streams for BRDF

      REAL(fpk), intent(out) :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: A_BRDF  ( MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER   :: I, I1

!  BRDF quadrature (Gauss-Legendre)
!  ---------------

!  Save these quantities for efficient coding

      CALL GETQUAD2 ( ZERO, ONE, NBRDF_HALF, X_BRDF, A_BRDF )
      DO I = 1, NBRDF_HALF
        I1 = I + NBRDF_HALF
        X_BRDF(I1) = - X_BRDF(I)
        A_BRDF(I1) =   A_BRDF(I)
      ENDDO
      DO I = 1, NSTREAMS_BRDF
        X_BRDF(I) = PIE * X_BRDF(I)
        CX_BRDF(I) = COS ( X_BRDF(I) )
        SX_BRDF(I) = SIN ( X_BRDF(I) )
      ENDDO

!  Finish

      RETURN
END SUBROUTINE BRDF_QUADRATURE_GAUSSIAN

!

SUBROUTINE BRDF_QUADRATURE_TRAPEZOID &
        ( NSTREAMS_BRDF,                     & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF )   ! Outputs

      IMPLICIT NONE

!  Input
!  =====

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS_BRDF

!  OUTPUT
!  ======

!  azimuth quadrature streams for BRDF

      REAL(fpk), intent(out) :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: A_BRDF  ( MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER    :: I, I1
      REAL(fpk)  :: DF1, DEL

!  BRDF quadrature (Trapezium)
!  ---------------

!  Save these quantities for efficient coding

      DF1 = DBLE(NSTREAMS_BRDF - 1 )
      DEL = TWO * PIE / DF1
      DO I = 1, NSTREAMS_BRDF
        I1 = I - 1
        X_BRDF(I) = DBLE(I1) * DEL - PIE
        X_BRDF(I) = DBLE(I1) * DEL
        CX_BRDF(I) = DCOS ( X_BRDF(I) )
        SX_BRDF(I) = DSIN ( X_BRDF(I) )
      ENDDO
      DO I = 2, NSTREAMS_BRDF - 1
        A_BRDF(I)  = DEL / PIE
      ENDDO
      A_BRDF(1)              = DEL * HALF / PIE
      A_BRDF(NSTREAMS_BRDF)  = DEL * HALF / PIE

!  Finish

      RETURN
END SUBROUTINE BRDF_QUADRATURE_TRAPEZOID

!  End module

      END MODULE brdf_sup_aux_m

