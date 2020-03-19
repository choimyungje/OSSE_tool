MODULE GEMSTOOL_numerical_m

!  Everything public
!   Double precision.

implicit none

!  Precision

integer, parameter :: fpk = SELECTED_REAL_KIND(15)

private :: fpk
public  :: linear_trawler_ascending,&
           spline,&
           splint,&
           lintp2,&
           makechar3

CONTAINS


  SUBROUTINE LINEAR_TRAWLER_ASCENDING &
         ( MAXGRID, NGRID, GRID, X,   &
           F1, F2, K1, K2 )


      INTEGER          :: MAXGRID, NGRID
      real(kind=fpk)        :: GRID(MAXGRID), X
      INTEGER          :: K1, K2
      real(kind=fpk)        :: F1, F2
      LOGICAL          :: LOOP

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


  SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER            ::  n
      REAL(kind=fpk)          :: yp1,ypn,x(n),y(n),y2(n)
      INTEGER, PARAMETER :: NMAX=6094
      INTEGER            :: i,k
      REAL(kind=fpk)          :: p,qn,sig,un,u(NMAX)

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
                     (y(i)-y(i-1)) / (x(i)-x(i-1))   &
                   ) / (x(i+1)-x(i-1)) - sig*u(i-1)  &
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
  END SUBROUTINE SPLINE


  SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER            ::  n
      REAL(kind=fpk)          :: x, y, xa(n),ya(n),y2a(n)
      INTEGER            ::  k,khi,klo
      REAL(kind=fpk)          :: a,b,h

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
!       if (h.eq.0.0d0) pause
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+ &
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0

      return
  END SUBROUTINE SPLINT


  SUBROUTINE LINTP2 (N1,X1,Y1,N2,X2,Y2)

      INTEGER        ::     N1,N2
      REAL(kind=fpk)      :: X1(N1),Y1(N1),X2(N2),Y2(N2)
      INTEGER        :: INV, N2P, J, I2, I1, J2, I
      REAL(kind=fpk)     :: X, YP, XP, S

      IF(N1.LE.0.OR.N2.LE.0.OR.N1.GT.32767.OR.N2.GT.32767)GOTO 9
!
! --- Checking the condition for inverting arrays ---
!
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

  subroutine makechar3(w,c3)
  integer, intent(in)          :: w
  character(LEN=3),intent(out) :: c3
  c3 = '000'
  if ( w.lt.10 ) then
     write(c3(3:3),'(I1)')W
  else if ( w.gt.99 ) then
     write(c3(1:3),'(I3)')W
  else
     write(c3(2:3),'(I2)')W
  endif
  end subroutine makechar3

END MODULE GEMSTOOL_numerical_m
    
