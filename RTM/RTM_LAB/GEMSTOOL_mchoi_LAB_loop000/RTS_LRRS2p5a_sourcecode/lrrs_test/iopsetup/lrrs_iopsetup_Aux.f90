
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

      MODULE LRRS_Iopsetup_Aux_m

!  Module of Auxiliary routines for the LRRS IOPSETUPs
!  ---------------------------------------------------

!  contains the following stand-alone routines (all public)
      
      ! LINEAR_TRAWLER_ASCENDING    Linear Interpolation
      ! spline                      (based on Numerical Recipes)
      ! splint                      (based on Numerical Recipes)
      ! LINTP2                      Linear Interpolation, old code
      ! CHANGE_RESOLUTION           Making coarser grids (RTS)

!  Note the exception handling is not great in some of these routines !!!

private :: fpk
public

contains

      SUBROUTINE LINEAR_TRAWLER_ASCENDING &
        ( MAXGRID, NGRID, GRID, X, &
          F1, F2, K1, K2 )

      implicit none
      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

      INTEGER, INTENT(IN)    :: MAXGRID, NGRID
      REAL(FPK), INTENT(IN)  :: GRID(MAXGRID), X

      INTEGER, INTENT(OUT)   :: K1, K2
      REAL(FPK), INTENT(OUT) :: F1, F2

      LOGICAL :: LOOP

      IF ( X .LT. GRID(2) ) THEN
        K1 = 1
      ELSE IF ( X .GE. GRID(NGRID-1) ) THEN
        K1 = NGRID - 1
      ELSE
        LOOP = .TRUE.
        K1 = 0
        DO WHILE (LOOP)
          K1 = K1 + 1
          IF ( X .GE. GRID(K1) .AND. X .LT. GRID(K1+1) ) LOOP=.FALSE.
        ENDDO
      ENDIF

      K2 = K1 + 1
      F1 = ( X - GRID(K1) ) / ( GRID(K2) - GRID(K1) )
      F2 = 1.0D0 - F1

      RETURN
      END SUBROUTINE LINEAR_TRAWLER_ASCENDING

!

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)

      implicit none
      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

      INTEGER, INTENT(IN) ::           n
      REAL(FPK), INTENT(IN)  :: yp1,ypn,x(n),y(n)
      REAL(FPK), INTENT(OUT) :: y2(n)

      INTEGER, PARAMETER  ::           NMAX = 6094
      INTEGER ::                       i,k
      REAL(FPK) ::              p,qn,sig,un,u(NMAX)

!  Check dimension

      if (n.gt.NMAX) then
        write(*,*)'Error in spline routine: too small NMAX =',NMAX
        stop
      endif

      if (yp1.gt..99e30) then
        y2(1)=0.0d0
        u(1)=0.0d0
      else
        y2(1)=-0.5d0
        u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.0d0*( (y(i+1)-y(i)) / (x(i+1)-x(i)) - &
                      (y(i)-y(i-1)) / (x(i)-x(i-1)) &
                    ) / (x(i+1)-x(i-1)) - sig*u(i-1) &
             )/p
      enddo

      if (ypn.gt..99d30) then
        qn=0.0d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo

      return
      END SUBROUTINE spline

!

      SUBROUTINE splint(xa,ya,y2a,n,x,y)

      implicit none
      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

      INTEGER, INTENT(IN)    :: n
      REAL(FPK), INTENT(IN)  :: x,xa(n),y2a(n),ya(n)
      REAL(FPK), INTENT(OUT) :: y

      INTEGER   ::              k,khi,klo
      REAL(FPK) ::              a,b,h

      klo=1
      khi=n

 1    if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x)then
            khi=k
          else
            klo=k
          endif
       goto 1
       endif
       h=xa(khi)-xa(klo)
       if (h.eq.0.0d0) stop 'splint'
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+ &
         ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0

      return
      END SUBROUTINE splint

!

      SUBROUTINE LINTP2 (N1,X1,Y1,N2,X2,Y2)

      implicit none
      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

      INTEGER, INTENT(IN) ::             N1,N2
      REAL(FPK), INTENT(IN) ::    X1(N1),Y1(N1)
      REAL(FPK), INTENT(INOUT) :: X2(N2)
      REAL(FPK), INTENT(OUT) ::   Y2(N2)

      INTEGER ::                       INV, N2P, J, I2, I1, J2, I
      REAL(FPK) ::              X, YP, XP, S

      IF(N1.LE.0.OR.N2.LE.0.OR.N1.GT.32767.OR.N2.GT.32767)GOTO 9
!
! --- Checking the condition for inverting arrays ---
!
      N2P=0
      INV=0
      IF (N1.EQ.1.OR.N2.EQ.1) GOTO 500
      IF ((X1(N1)-X1(1))*(X2(N2)-X2(1))) 300,300,500
!
! --- Inversion of new grid ---
!
      N2P=0
  300 CONTINUE
      INV=1
      N2P=N2/2
      DO 301 J=1,N2P
      XP=X2(J)
      X2(J)=X2(N2-J+1)
      X2(N2-J+1)=XP
  301 CONTINUE
!
! --- Main block ---
!
  500 IF (N1.EQ.1) GOTO 7
      S=DSIGN(1.0D0,X1(N1)-X1(1))
      I2=1
      I1=2
    1 IF((X2(I2)-X1(1))*S) 2,2,3
    2 Y2(I2)=Y1(1)
      I2=I2+1
      IF(I2.GT.N2) GOTO 999
      GOTO 1
    3 IF((X2(I2)-X1(I1))*S) 4,4,5
    4 X=(X2(I2)-X1(I1-1))/(X1(I1)-X1(I1-1))
      Y2(I2)=Y1(I1)*X+Y1(I1-1)*(1.0D0-X)
      I2=I2+1
      IF(I2.GT.N2) GOTO 999
      GOTO 3
    5 I1=I1+1
      IF(I1.LE.N1) GOTO 3
      DO 6 J2=I2,N2
    6 Y2(J2)=Y1(N1)
      GOTO 999
    7 DO 8 I=1,N2
    8 Y2(I)=Y1(1)
  999 CONTINUE
!
! --- Checking the condition for back inversion ---
!
      IF (INV.NE.1) GOTO 1000
!
! --- Back inversion ---
!

      DO 302 J=1,N2P
      XP=X2(J)
      YP=Y2(J)
      X2(J)=X2(N2-J+1)
      Y2(J)=Y2(N2-J+1)
      X2(N2-J+1)=XP
      Y2(N2-J+1)=YP
  302 CONTINUE
!
! --- Exit block ---
!
 1000 RETURN
    9 PRINT 100,N1,N2
  100 FORMAT(10X,'  Error in subroutine LINTP2 !',2I15)
      STOP '  Error in subroutine LINTP2 !'

      END SUBROUTINE LINTP2

!

      SUBROUTINE CHANGE_RESOLUTION ( &
        MAX_LOCAL_LAYERS, NRES, FINER_RES, HALF_RES, QTR_RES, &
        DATA_NHEIGHTS, DATA_HEIGHTS, DATA_TEMPS )

!  Change resolution, either By Introducing layers (NRES, FINER_RES)
!                     or by making Half or Quarter (HALF_RES, QTR_RES)

      IMPLICIT NONE
      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input

      INTEGER, INTENT(IN) :: MAX_LOCAL_LAYERS, NRES
      LOGICAL, INTENT(IN) :: FINER_RES, HALF_RES, QTR_RES

!  Modified input

      INTEGER, INTENT(INOUT) ::   DATA_NHEIGHTS
      REAL(FPK), INTENT(INOUT) :: DATA_HEIGHTS (MAX_LOCAL_LAYERS)
      REAL(FPK), INTENT(INOUT) :: DATA_TEMPS   (MAX_LOCAL_LAYERS)

!  Local

      INTEGER ::   P, N, NC
      REAL(FPK) :: DHITES (MAX_LOCAL_LAYERS)
      REAL(FPK) :: DTEMPS (MAX_LOCAL_LAYERS)
      REAL(FPK) :: f1(4), f2(4)

!  Finer layering

      IF ( FINER_RES ) THEN
        DO P = 1, NRES - 1
          F1(P) = DBLE(NRES-P) / DBLE(NRES)
          F2(P) = DBLE(P)      / DBLE(NRES)
        ENDDO
        NC = 1
        DHITES(1) = DATA_HEIGHTS(1)
        DTEMPS(1) = DATA_TEMPS  (1)
        DO N = 2, DATA_NHEIGHTS
          DO P = 1, NRES - 1
           DHITES(NC+P) = F1(P)*DATA_HEIGHTS(N-1)+F2(P)*DATA_HEIGHTS(N)
           DTEMPS(NC+P) = F1(P)*DATA_TEMPS(N-1)  +F2(P)*DATA_TEMPS(N)
          ENDDO
          DHITES(NC+NRES) = DATA_HEIGHTS(N)
          DTEMPS(NC+NRES) = DATA_TEMPS(N)
          NC = NC + NRES
        ENDDO
        DATA_NHEIGHTS = NC
        DO N = 1, DATA_NHEIGHTS
          DATA_HEIGHTS(N) = DHITES(N)
          DATA_TEMPS(N)   = DTEMPS(N)
        ENDDO
      ENDIF

!  Keep one out of every 2 layers

      IF ( HALF_RES ) THEN
        NC = 1
        DHITES(1) = DATA_HEIGHTS(1)
        DTEMPS(1) = DATA_TEMPS  (1)
        DO N = 3, DATA_NHEIGHTS, 2
          NC = NC + 1
          DHITES(NC) = DATA_HEIGHTS(N)
          DTEMPS(NC) = DATA_TEMPS(N)
        ENDDO
        DATA_NHEIGHTS = NC
        DO N = 1, DATA_NHEIGHTS
          DATA_HEIGHTS(N) = DHITES(N)
          DATA_TEMPS(N)   = DTEMPS(N)
        ENDDO
      ENDIF

!  Keep one out of every 4 layers

      IF ( QTR_RES ) THEN
        NC = 1
        DHITES(1) = DATA_HEIGHTS(1)
        DTEMPS(1) = DATA_TEMPS  (1)
        DO N = 3, DATA_NHEIGHTS, 4
          NC = NC + 1
          DHITES(NC) = DATA_HEIGHTS(N)
          DTEMPS(NC) = DATA_TEMPS(N)
        ENDDO
        DATA_NHEIGHTS = NC
        DO N = 1, DATA_NHEIGHTS
          DATA_HEIGHTS(N) = DHITES(N)
          DATA_TEMPS(N)   = DTEMPS(N)
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE CHANGE_RESOLUTION

!  End module

      END MODULE LRRS_Iopsetup_Aux_m
