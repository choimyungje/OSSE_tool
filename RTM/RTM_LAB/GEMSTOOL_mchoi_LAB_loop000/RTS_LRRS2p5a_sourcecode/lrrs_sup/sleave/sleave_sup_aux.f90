
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
! #  This is sleave_sup_aux.f90.  Utility routines.             #
! #  The subroutines are listed below with their                #
! #  source of origin (order of appearance).                    #
! #                                                             #
! #      GETQUAD2              M. Christi, 2017                 #
! #      HOMEGROWN_ERRFUNC     R. Spurr (modified from          #
! #                               DERFC_E by V. Natraj, 2005)   #
! #                                                             #
! ###############################################################

      MODULE sleave_sup_aux_m

      USE LRRS_pars_m, only : FPK, ZERO, ONE, TWO, HALF, QUARTER, PIE

      PRIVATE
      PUBLIC :: GETQUAD2, &
                HOMEGROWN_ERRFUNC

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

      SUBROUTINE HOMEGROWN_ERRFUNC ( X )

      IMPLICIT NONE

      REAL(FPK), intent(inout) :: x

! Returns the complementary error function erfc(x) with fractional error
! everywhere less than 1.2 * 10^7. 

      REAL(FPK) :: t,z,xin

      z = abs(x) ; xin = x
      t = 1.d0/(1.d0+0.5d0*z)
      x = t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+ &
              t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t* &
              (-1.13520398d0+t*(1.48851587d0+t*(-.82215223d0+t* &
              .17087277d0)))))))))
      if ( xin .lt. 0.d0 ) x = 2.d0 - x

      return
      END SUBROUTINE HOMEGROWN_ERRFUNC

!  End module

      END MODULE sleave_sup_aux_m

