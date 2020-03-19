Module SMCConv_Interp_m

private 
public  SMCConv_Interp, SMCConv_Interp_2, SMC_GAULEG, SMC_spline, SMC_splint

contains

subroutine SMCConv_Interp &
    ( Max_InAngles, Max_OutAngles, do_Quadrature, dtr, &
      N_OutAngles, N_InAngles, InAngles, InScatMat,    &
      OutAngles, OutCosines, OutWeights, OutScatMat )

!  Stand-alone routine to interpret Scattering Matrix onto Quadrature grid

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensioning

   Integer, intent(in) :: Max_InAngles,  Max_OutAngles

!  Logical for Quadrature

   Logical, intent(in) :: do_Quadrature

!  Angle numbers

   Integer, intent(in) :: N_InAngles,  N_OutAngles

!  Input angles and scattering matrix entries

   Real(dpk), intent(in) :: dtr, InAngles ( Max_InAngles )
   Real(dpk), intent(in) :: InScatMat ( Max_InAngles, 6 )

!  Outputs
!  -------

!   Quadrature

   Real(dpk), intent(inout) :: OutAngles ( Max_OutAngles )
   Real(dpk), intent(inout) :: OutCosines ( Max_OutAngles )
   Real(dpk), intent(inout) :: OutWeights ( Max_OutAngles )

!  Interpolated scattering matrix entries

   Real(dpk), intent(out) :: OutScatMat ( Max_OutAngles, 6 )

!  Local
!  -----

   real(dpk), parameter :: zero = 0.0_dpk, one = 1.0_dpk, mone = -1.0_dpk
   logical    :: loop
   integer    :: n, m, k, k1, k2, kstart
   real(dpk)  :: InCosines(Max_InAngles), x, f1, f2, x1, x2, w1, w2

!  Start of Code
!  -------------

!  initialize

   OutScatMat = zero
   if ( do_Quadrature ) then
      OutAngles  = zero
      OutCosines = zero
      OutWeights = zero
   endif

!  Develop Quadrature

   if ( do_Quadrature ) then
      Call SMC_GAULEG ( mone, one, OutCosines, OutWeights, N_OutAngles, Max_OutAngles )
      do k1 = 1, N_OutAngles / 2
         k2 = N_OutAngles + 1 - k1 
         x1 = OutCosines(k1) ; x2 = OutCosines(k2)
         w1 = OutWeights(k1) ; w2 = OutWeights(k2)
         OutCosines(k1) = x2 ; OutCosines(k2) = x1
         OutWeights(k1) = w2 ; OutWeights(k2) = w1
      enddo
      do k = 1, N_OutAngles
         OutAngles(k) = acos ( OutCosines(k) ) / dtr
      enddo
   endif

!  Cosines of input angles

   do k = 1, N_InAngles
      InCosines(k) = cos ( InAngles(k) * dtr )
   enddo

!  interpolate

   kstart = 0
   do n = 1, N_OutAngles
      x = OutCosines(n)
      k = kstart ; loop = .true.
      do while (loop .and.k.lt.N_InAngles)
         k = k + 1 
         loop = (x.lt.InCosines(k))
      enddo
      kstart = k - 1 ; k1 = k - 1 ; k2 = k
      f1 = ( InCosines(k2) - x ) / ( InCosines(k2) - InCosines(k1) ) ; f2 = one - f1
      do m = 1, 6
         OutScatmat(n,m) = f1 * InScatmat(k1,m) + f2 * InScatmat(k2,m)
      enddo
   enddo

!  Done

   return
end subroutine SMCConv_Interp

!
subroutine SMCConv_Interp_2 &
    ( Max_InAngles, Max_OutAngles, do_Quadrature, dtr, &
      N_OutAngles, N_InAngles, InAngles, InScatMat,    &
      OutAngles, OutCosines, OutWeights, OutScatMat )

!  Stand-alone routine to interpret Scattering Matrix onto Quadrature grid
!  This one based on Splines

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensioning

   Integer, intent(in) :: Max_InAngles,  Max_OutAngles

!  Logical for Quadrature

   Logical, intent(in) :: do_Quadrature

!  Angle numbers

   Integer, intent(in) :: N_InAngles,  N_OutAngles

!  Input angles and scattering matrix entries

   Real(dpk), intent(in) :: dtr, InAngles ( Max_InAngles )
   Real(dpk), intent(in) :: InScatMat ( Max_InAngles, 6 )

!  Outputs
!  -------

!   Quadrature

   Real(dpk), intent(inout) :: OutAngles  ( Max_OutAngles )
   Real(dpk), intent(inout) :: OutCosines ( Max_OutAngles )
   Real(dpk), intent(inout) :: OutWeights ( Max_OutAngles )

!  Interpolated scattering matrix entries

   Real(dpk), intent(out) :: OutScatMat ( Max_OutAngles, 6 )

!  Local
!  -----

   real(dpk), parameter :: zero = 0.0_dpk, one = 1.0_dpk, mone = -1.0_dpk
   integer    :: m, k, k1, k2
   real(dpk)  :: InCosines(Max_InAngles), Local_InScatMat(Max_InAngles,6)
   real(dpk)  :: x1,x2,w1,w2,x,y,y2(Max_InAngles,6),yp1,ypn

!  Start of Code
!  -------------

!  initialize

   OutScatMat = zero
   if ( do_Quadrature ) then
      OutAngles  = zero
      OutCosines = zero
      OutWeights = zero
   endif

!  Develop Quadrature

   if ( do_Quadrature ) then
      Call SMC_GAULEG ( mone, one, OutCosines, OutWeights, N_OutAngles, Max_OutAngles )
      do k1 = 1, N_OutAngles / 2
         k2 = N_OutAngles + 1 - k1 
         x1 = OutCosines(k1) ; x2 = OutCosines(k2)
         w1 = OutWeights(k1) ; w2 = OutWeights(k2)
         OutCosines(k1) = x2 ; OutCosines(k2) = x1
         OutWeights(k1) = w2 ; OutWeights(k2) = w1
      enddo
      do k = 1, N_OutAngles
         OutAngles(k) = acos ( OutCosines(k) ) / dtr
      enddo
   endif

!  Cosines of input angles.
!    -- Reverse directions for the Splining (Monotonically increasing)

   do k = 1, N_InAngles
     k1 = N_InAngles + 1 - k 
     InCosines(k1) = cos ( InAngles(k) * dtr )
     Local_Inscatmat(k1,1:6) = Inscatmat(k,1:6)  
   enddo

!  Develop Splines
!    - Set the End-point gradient (YPN) = input gradient at forward-peak
!    - This make a HUGE difference to the accuracy

    do m = 1, 6
       yp1 = zero
       ypn = ( Local_InScatMat(N_InAngles,m) -  Local_InScatMat(N_InAngles-1,m) ) / &
              ( InCosines(N_InAngles) -  InCosines(N_InAngles-1) )
       Call SMC_SPLINE(Max_InAngles,InCosines,Local_InScatMat(:,m),N_InAngles,yp1,ypn,y2(:,m))
    enddo

!  interpolate

   do k = 1, N_OutAngles
      k1 = N_OutAngles + 1 - k 
      x = OutCosines(k1)
      do m = 1, 6
         Call SMC_SPLINT(Max_InAngles,InCosines,Local_InScatMat(:,m),y2(:,m),N_InAngles,x,y)
         OutScatMat(k1,m) = y
      enddo
   enddo

!  Done

   return
end subroutine SMCConv_Interp_2

!

SUBROUTINE SMC_GAULEG(X1,X2,X,W,N,NMAX)

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Inputs

   INTEGER  , intent(in)  :: N, NMAX
   Real(dpk), intent(in)  :: X1, X2

!  Outputs

   Real(dpk)  :: X(NMAX),W(NMAX)

!  Local

   INTEGER    :: I, M, J
   Real(dpk)  :: XM,XL,P1,P2,P3,PP,Z,Z1,DN,DJ, PIE
   Real(dpk), PARAMETER :: EPS=3.0E-14_dpk
   Real(dpk), PARAMETER :: half = 0.5_dpk, quarter = 0.25_dpk
   Real(dpk), PARAMETER :: one  = 1.0_dpk, two     = 2.0_dpk

!  Code

   M  = (N+1)/2 ; DN = real(N,dpk) ; PIE = ACOS(-one)
   XM = half*(X2+X1)
   XL = half*(X2-X1)
   DO I = 1, M
      Z = COS ( pie * ( real(I,dpk) - quarter ) / ( real(N,dpk) + half ) )
1     CONTINUE
      P1 = one ; P2 = 0.0_dpk
      DO J = 1, N
         P3 = P2 ; P2 = P1 ; DJ = real(J,dpk)
         P1 = ( (two*DJ-one) * Z * P2 - (DJ-one) * P3) / DJ
      ENDDO
      PP = DN * (Z*P1-P2) / (Z*Z-one)
      Z1 = Z ; Z = Z1 - P1 / PP
      IF (ABS(Z-Z1).GT.EPS) GO TO 1
      X(I)     = XM - XL * Z
      X(N+1-I) = XM + XL * Z
      W(I)     = two * XL / ( (one-Z*Z)*PP*PP )
      W(N+1-I) = W(I)
   ENDDO
   RETURN
END SUBROUTINE SMC_GAULEG


SUBROUTINE SMC_SPLINE(nmax,x,y,n,yp1,ypn,y2)

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input

   INTEGER        :: nmax !max # of x & y points
   INTEGER        :: n !actual # of x & y points
   REAL(kind=dpk) :: x(nmax),y(nmax) !pts of orig func
   REAL(kind=dpk) :: yp1,ypn  !1st derivs of orig func at endpts

!  output

   REAL(kind=dpk) :: y2(nmax) !2nd derivs of the spline interp func

!  local

   INTEGER        :: i,k
   REAL(kind=dpk) :: p,qn,sig,un,u(nmax)

  if (yp1.gt..99e30) then
    y2(1) = 0.0_dpk
    u(1)  = 0.0_dpk
  else
    y2(1)=- 0.5_dpk
    u(1)=(3.0_dpk/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif

  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.0_dpk
    y2(i)=(sig-1.)/p
    u(i)=(6.0_dpk*( (y(i+1)-y(i)) / (x(i+1)-x(i)) - &
                  (y(i)-y(i-1)) / (x(i)-x(i-1))   &
                ) / (x(i+1)-x(i-1)) - sig*u(i-1)  &
         )/p
  enddo

  if (ypn.gt..99d30) then
    qn=0.0_dpk
    un=0.0_dpk
  else
    qn=0.5_dpk
    un=(3.0_dpk/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif

  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo

END SUBROUTINE SMC_SPLINE


SUBROUTINE SMC_SPLINT(nmax,xa,ya,y2a,n,x,y)

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input

   INTEGER        :: nmax !max # of x & y points
   INTEGER        :: n !actual # of x & y points
   REAL(kind=dpk) :: xa(nmax),ya(nmax) !pts of orig func
   REAL(kind=dpk) :: y2a(nmax) !2nd derivs of the spline interp func
   REAL(kind=dpk) :: x !x value for which interp y value is desired

!  output

   REAL(kind=dpk) :: y !interp output value

!  local

   INTEGER        :: k,khi,klo
   REAL(kind=dpk) :: a,b,h

   klo=1
   khi=n

1  continue
   if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
         khi=k
      else
         klo=k
      endif
      goto 1
   endif
   h=xa(khi)-xa(klo)
!  if (h.eq.0.0_dpk) pause
   a=(xa(khi)-x)/h
   b=(x-xa(klo))/h
   y=a*ya(klo)+b*ya(khi)+ &
    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dpk

END SUBROUTINE SMC_SPLINT

!  End module

End Module SMCConv_Interp_m
