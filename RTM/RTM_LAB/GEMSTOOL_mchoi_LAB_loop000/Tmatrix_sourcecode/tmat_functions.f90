module tmat_functions

 use tmat_parameters, only : fpk => tmat_fpkind, NPNG1, NPNG2, NPN1 

!  LIST OF SUBROUTINES
!  ===================

!    Eqv_radius_chebyshev
!    Eqv_radius_cylinder
!    Eqv_radius_spheroids

!    tmatrix_constants

!    tmatrix_vary             ----- calling ----
!       spherical_coords
!       spheroid_coords
!       chebyshev_coords
!       cylinder_coords
!       Bessel_functions_master  ----- calling ----
!          Bessel_rj
!          Bessel_ry
!          Bessel_Cspher

!    GAULEG

!  Top level routines are PUBLIC here, along with others as required

private
public        tmatrix_constants, tmatrix_vary,    &
              GAULEG_right,                       &
              Eqv_radius_chebyshev,               &
              Eqv_radius_cylinder,                &
              Eqv_radius_spheroids,               &
              spheroid_coords,  spherical_coords, &
              chebyshev_coords, cylinder_coords

contains

subroutine tmatrix_constants            &
     ( ngauss, nmax, np, aspect,        & ! inputs
       x, w, an, ann, s, ss )             ! outputs

   use tmat_parameters, only : NPNG1, NPNG2, NPN1

   implicit none

!  input arguments
!  ---------------

!  Quadrature numbers

   integer, intent(in) :: ngauss
   integer, intent(in) :: nmax

!  Particle shape parameterization
!             For spheroids NP=-1 and ASPECT is the ratio of the          
!                 horizontal to rotational axes.  ASPECT is larger than   
!                 1 for oblate spheroids and smaller than 1 for       
!                 prolate spheroids.                                   
!             For cylinders NP=-2 and ASPECT is the ratio of the          
!                 diameter to the length.                              
!             For Chebyshev particles NP must be positive and 
!                 is the degree of the Chebyshev polynomial, while     
!                 ASPECT is the deformation parameter                     

   integer  , intent(in) :: np
   real(fpk), intent(in) :: aspect

!  Output arguments
!  ----------------

!  Quadrature

   REAL(fpk), intent(out) :: X(NPNG2), W(NPNG2)
   REAL(fpk), intent(out) :: S(NPNG2), SS(NPNG2)

!  Constants

   REAL(fpk), intent(out) :: AN(NPN1), ANN(NPN1,NPN1)

!  Local variables
!  ---------------

!  Help arrays

   REAL(fpk) :: X1(NPNG1),W1(NPNG1), DD(NPN1)
   REAL(fpk) :: X2(NPNG1),W2(NPNG1)

!  other variables

   INTEGER   :: N1, N, NN, NG, NG1, NG2, I
   REAL(fpk) :: D, DDD, XX, Y

!  Zero the output

   X    = 0.0d0
   W    = 0.0d0
   S    = 0.0d0
   SS   = 0.0d0
   AN   = 0.0d0
   ANN  = 0.0d0

!  Set up the local constants

   DO N = 1, NMAX
      NN = N*(N+1)
      AN(N) = DBLE(NN)
      D = DSQRT(DBLE(2*N+1)/DBLE(NN))
      DD(N) = D
      DO N1 = 1,N
         DDD = D*DD(N1)*0.5D0
         ANN(N,N1) = DDD
         ANN(N1,N) = DDD
      ENDDO
   ENDDO

!  Twice the number of quadratures

   NG=2*NGAUSS

!  Cylinders. Double quadrature

   if ( np .eq. -2 ) then

      NG1 = DBLE(NGAUSS)/2.0D0
      NG2 = NGAUSS - NG1
      XX  = -DCOS(DATAN(aspect))
!      CALL GAULEG_wrong (-1.0d0,1.0d0,X1,W1,NG1)
!      CALL GAULEG_wrong (-1.0d0,1.0d0,X2,W2,NG2)
      CALL GAULEG_right (NG1,0,0,X1,W1)
      CALL GAULEG_right (NG2,0,0,X2,W2)

      DO I = 1,NG1
         W(I)=0.5D0*(XX+1.0D0)*W1(I)
         X(I)=0.5D0*(XX+1.0D0)*X1(I)+0.5D0*(XX-1.0D0)
      ENDDO

      DO I = 1,NG2
         W(I+NG1)=-0.5D0*XX*W2(I)
         X(I+NG1)=-0.5D0*XX*X2(I)+0.5D0*XX
      ENDDO

      DO I = 1,NGAUSS
         W(NG-I+1)=W(I)
         X(NG-I+1)=-X(I)
      ENDDO

!  Spheroids, Chebyshevs (Single quadrature, whole range)

   else

!      CALL GAULEG_wrong (-1.0d0,1.0d0,X,W,NG)
      CALL GAULEG_right (NG,0,0,X,W)

   endif

!  S vectors (same for all types)

   DO I = 1,NGAUSS
      Y = X(I)
      Y = 1D0/(1D0-Y*Y)
      SS(I)      = Y
      SS(NG-I+1) = Y
      Y = DSQRT(Y)
      S(I)       = Y
      S(NG-I+1)  = Y
   ENDDO

!  Finish

   return
end subroutine tmatrix_constants

subroutine tmatrix_vary                          &
     ( n_real, n_imag, radius, PI, x,            & ! inputs
       aspect, np, ngauss, nmax,                 & ! inputs
       R, DR, DDR, DRR, DRI,                     & ! outputs
       j_bess,  y_bess,  jr_bess,  ji_bess,      & ! Outputs
       dj_bess, dy_bess, djr_bess, dji_bess )      ! Outputs
     
   use tmat_parameters, only : NPNG1, NPNG2, NPN1

   implicit none

!  input arguments
!  ---------------

   integer  , intent(in)    :: np, ngauss, nmax
   real(fpk), intent(in)    :: n_real, n_imag, radius, aspect, pi
   real(fpk), intent(inout) :: X(NPNG2)

!  Output arguments
!  ----------------

!  R and DR

   real(fpk), intent(out) :: R(NPNG2), DR(NPNG2)
   real(fpk), intent(out) :: DDR(NPNG2), DRR(NPNG2), DRI(NPNG2)

!  Bessel Master output

   real(fpk), intent(out) :: J_BESS  (NPNG2,NPN1)
   real(fpk), intent(out) :: Y_BESS  (NPNG2,NPN1)
   real(fpk), intent(out) :: JR_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: JI_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DJ_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DY_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DJR_BESS(NPNG2,NPN1)
   real(fpk), intent(out) :: DJI_BESS(NPNG2,NPN1)

!  Local variables
!  ---------------

!  Local arrays

   real(fpk) :: Z(NPNG2),ZR(NPNG2),ZI(NPNG2)

!  Other lcoal variables

   integer   :: i, ng, NNMAX1, NNMAX2
   real(fpk) :: Modulus, PRR, PRI
   real(fpk) :: TA, TB, V, VV, V1, V2

!  Constants
!  ---------

   NG = NGAUSS*2

   Modulus  = n_real*n_real + n_imag*n_imag
   PRR = + n_real / Modulus
   PRI=  - n_imag / Modulus

!  Coordinates
!  -----------

   if ( np .eq. -1 ) then       ! Spheroids
      call spheroid_coords                &
        ( ng, ngauss, radius, aspect, x,  & ! inputs
          r, dr )                           ! Outputs
   else if ( np .eq. -2 ) then  !  Cylinders
      call cylinder_coords                &
        ( ng, ngauss, radius, aspect, x,  & ! inputs
          r, dr )                           ! Outputs
   else if ( np .gt. 0 ) then   !  Chebyshev
      call chebyshev_coords               &
        ( ng, radius, aspect, np, x,      & ! inputs
          r, dr )                           ! Outputs
   else if ( np .eq. 0 ) then   !  Mie
      call spherical_coords                &
        ( ng, radius,                      & ! inputs
          r, dr )                            ! Outputs
   endif

!  loop for DDR, DRR, DRI, ZR, ZI
!  ------------------------------

   TA=0D0
   DO I = 1, NG
        VV=DSQRT(R(I))
        V=VV*PI
        TA=MAX(TA,V)
        VV=1D0/V
        DDR(I)=VV
        DRR(I) = PRR*VV
        DRI(I) = PRI*VV
        V1 = V * n_real
        V2 = V * n_imag
        Z(I)  = V
        ZR(I) = V1
        ZI(I) = V2
   ENDDO

!  Maximum numbers for exapnsions

   TB = TA*DSQRT(Modulus)
   TB = DMAX1 ( TB,DBLE(NMAX) )

   NNMAX1 = 1.2D0*DSQRT(DMAX1(TA,DBLE(NMAX)))+3D0
   NNMAX2 = (TB+4D0*(TB**0.33333D0)+1.2D0*DSQRT(TB))
   NNMAX2 = NNMAX2-NMAX+5

!  Bessel functions

   call Bessel_functions_master                  &
     ( ng, nmax, nnmax1, nnmax2, z, zr, zi,      & ! Inputs
        j_bess,  y_bess,  jr_bess,  ji_bess,     & ! Outputs
       dj_bess, dy_bess, djr_bess, dji_bess )      ! Outputs

!  Finish

   return
end subroutine tmatrix_vary


subroutine Eqv_radius_chebyshev  ( norder, deformation, radius )

   implicit none

!  Input

   integer  , intent(in) ::  norder
   REAL(fpk), intent(in)  :: deformation

!  Output

   REAL(fpk), intent(inout)  :: radius

!  Local

   integer    :: i, ngauss
   REAL(fpk)  :: X(60),W(60)
   real(fpk)  :: dn,e,e2,en,s,v,xi,dx,dxn,ds,dsn,dcn,a,a2,ens,rs,rv

!  Constants

   E  = deformation
   DN = DBLE(norder)
   E2 = E*E
   EN = E*DN

!  Quadrature

   ngauss = 60
!   CALL GAULEG_wrong (-1.0d0,1.0d0,X,W,ngauss)
   CALL GAULEG_right (NGAUSS,0,0,X,W)

!  Summation for surface-area and volume

   S = 0.0D0
   V = 0.0D0
   DO i = 1, ngauss
      XI = X(I)
      DX = DACOS(XI)
      DXN = DN*DX
      DS  = DSIN(DX)
      DSN = DSIN(DXN)
      DCN = DCOS(DXN)
      A  = 1.0D0 + E * DCN
      A2 = A * A
      ENS = EN * DSN
      S = S+W(I) * A * DSQRT(A2+ENS*ENS)
      V = V+W(I) * ( DS*A + XI*ENS ) * DS * A2
   ENDDO

!  final answer for the radius

   RS  = DSQRT(S*0.5D0)
   RV  = (V*3.0D0/4.0D0)**(1.0D0/3.0D0)
   radius = RV/RS

!  Finish

   RETURN
end subroutine Eqv_radius_chebyshev

subroutine Eqv_radius_cylinder ( Diamlength_ratio, radius )

   implicit none

!  input argument

   real(fpk), intent(in) :: Diamlength_ratio

!  Output argument (equivalent radius)

   real(fpk), intent(inout) :: radius

!  formula

   radius = (1.5D0/Diamlength_ratio)**(1.0D0/3.0D0)
   radius = radius / DSQRT((Diamlength_ratio+2.0D0)/(2.0D0*Diamlength_ratio))

!  Finish

   return   
end subroutine Eqv_radius_cylinder

subroutine Eqv_radius_spheroids ( aspect_ratio, radius )

   implicit none

!  input argument

   real(fpk), intent(in) :: aspect_ratio

!  Output argument (equivalent radius)

   real(fpk), intent(inout) :: radius

!  Local variables

   real(fpk) :: sinphi, phi, r1, r2, ratiosq

!  Prolate spheroids

   if ( aspect_ratio .lt. 1.0d0 ) then
      sinphi = dsqrt(1.0d0 - aspect_ratio * aspect_ratio )
      phi    = dasin(sinphi)
      r1 = 0.5d0 * aspect_ratio ** ( 2.0d0 / 3.0d0 )
      r2 = 1.0d0 + ( phi / sinphi / aspect_ratio )
      radius = 1.0d0 / dsqrt(r1*r2)
   endif

!  Oblate spheroids

   if ( aspect_ratio .gt. 1.0d0 ) then
      ratiosq = aspect_ratio * aspect_ratio
      sinphi = dsqrt ( ratiosq - 1.0d0 ) / aspect_ratio
      phi    = dlog((1.0d0+sinphi)/(1.0d0-sinphi))
      r1 = 0.25d0 * aspect_ratio ** ( 2.0d0 / 3.0d0 )
      r2 = 2.0d0 + ( phi / ratiosq / sinphi )
      radius = 1.0d0 / dsqrt(r1*r2)
   endif

!  Sphere

   if ( aspect_ratio .eq. 1.0d0 ) then
      radius = 1.0d0
   endif

!  Finish

   return
end subroutine Eqv_radius_spheroids

subroutine spheroid_coords       &
     ( ng, ngauss, rev, eps, x,  & ! inputs
       r, dr )                     ! Outputs

!  Coordinates for the spheroidal particles

   implicit none

!  inputs

   integer  , intent(in)  :: ng, ngauss
   real(fpk), intent(in)  :: rev, eps
   real(fpk), intent(in)  :: x(ng)

!  outputs

   real(fpk), intent(out) :: r (ng)
   real(fpk), intent(out) :: dr(ng)

!  Local

   integer   :: i, i1
   real(fpk) :: a, aa, ee, ee1, c, cc, ss, s, rr

!  initialize

   R = 0.0d0   ; DR = 0.0d0

!  Constants for this rotuine

   A   = REV * EPS**(1.0D0/3.0D0)
   AA  = A*A
   EE  = EPS*EPS
   EE1 = EE - 1.0D0

!  Loop over quadrautre radii

   DO I = 1, NGAUSS
      I1 = NG-I+1
      C  = X(I)
      CC = C*C
      SS = 1.0D0-CC
      S  = DSQRT(SS)
      RR = 1.0D0/(SS+EE*CC)
      R (I) = AA*RR        ; R (I1) =   R (I)
      DR(I) = RR*C*S*EE1   ; DR(I1) = - DR(I)
   ENDDO

!  Finish

   return
end subroutine spheroid_coords

subroutine chebyshev_coords          &
     ( ng, rev, eps, np, x,          & ! inputs
       r, dr )                         ! Outputs

!  Coordinates for the Chebyshev particles

   implicit none

!  inputs

   integer  , intent(in)  :: ng, np
   real(fpk), intent(in)  :: rev, eps
   real(fpk), intent(in)  :: x(ng)

!  outputs

   real(fpk), intent(out) :: r (ng)
   real(fpk), intent(out) :: dr(ng)

!  Local

   integer   :: i
   real(fpk) :: a, dn, dnp, dn4, ep, r0, xi, ri

!  initialize

   R = 0.0d0   ; DR = 0.0d0

!  Constants

   DNP = DBLE(NP)
   DN  = DNP*DNP
   DN4 = DN*4.0D0
   EP  = EPS*EPS
   A = 1.0D0 + 1.5D0 * EP * (DN4-2.0D0)/(DN4-1.0D0)

!  Is this really an integer ?

!   I = DNP + 0.1D0 ) * 0.5D0 )                ! Original code
   I = INT ( ( DNP + 0.1D0 ) * 0.5D0 )
   I = 2*I
   IF (I.EQ.NP) THEN
      A = A - 3.0D0*EPS*(1.0D0+0.25D0*EP)/ &
                 (DN-1.0D0)-0.25D0*EP*EPS/(9.0D0*DN-1.0D0)
   ENDIF

   R0 = REV * A ** (-1.0D0/3.0D0)

!  Quadrature loop

   DO I = 1, NG
      XI = DACOS(X(I))*DNP
      RI = R0 * (1D0+EPS*DCOS(XI))
      R(I)  = RI*RI
      DR(I) = -R0*EPS*DNP*DSIN(XI)/RI
   ENDDO

!  Finish

   return
end subroutine chebyshev_coords

subroutine spherical_coords          &
     ( ng, rev,                      & ! inputs
       r, dr )                         ! Outputs

!  Coordinates for Spherical particles (Mie code)

   implicit none

!  inputs

   integer  , intent(in)  :: ng
   real(fpk), intent(in)  :: rev

!  outputs

   real(fpk), intent(out) :: r (ng)
   real(fpk), intent(out) :: dr(ng)

!  Local

   integer      :: i
   real(fpk) :: r0

   R = 0.0d0   ; DR = 0.0d0
   R0 = REV * REV
   DO I = 1, NG
      R(I)  = R0
   ENDDO

!  Finish

   return
end subroutine spherical_coords

subroutine cylinder_coords       &
     ( ng, ngauss, rev, eps, x,  & ! inputs
       r, dr )                     ! Outputs

!  Coordinates for the Cylindrical particles

   implicit none

!  inputs

   integer  , intent(in)  :: ng, ngauss
   real(fpk), intent(in)  :: rev, eps
   real(fpk), intent(in)  :: x(ng)

!  outputs

   real(fpk), intent(out) :: r (ng)
   real(fpk), intent(out) :: dr(ng)

!  Local

   integer   :: i, i1
   real(fpk) :: a, h, co, si, rad, rthet

!  initialize

   R = 0.0d0   ; DR = 0.0d0

!  Constants

   H = REV*( (2D0/(3D0*EPS*EPS))**(1D0/3D0) )
   A = H*EPS

!  Quadrature loop

   DO I = 1, NGAUSS
      I1 = NG-I+1
      CO = -X(I)
      SI = DSQRT(1.0D0-CO*CO)
      IF (SI/CO.GT.A/H) THEN
         RAD   = A/SI
         RTHET = -A*CO/(SI*SI)
      ELSE
         RAD   = H/CO
         RTHET = H*SI/(CO*CO)
      ENDIF
      R(I)  =  RAD*RAD       ; R(I1)  =    R(I)
      DR(I) = -RTHET/RAD     ; DR(I1) = - DR(I)
   ENDDO

!  Finish

   return
end subroutine cylinder_coords

subroutine Bessel_functions_master               &
     ( ng, nmax, nnmax1, nnmax2, x, xr, xi,      & ! Inputs
        j_bess,  y_bess,  jr_bess,  ji_bess,     & ! Outputs
       dj_bess, dy_bess, djr_bess, dji_bess )      ! Outputs

   use tmat_parameters, only : NPNG2, NPN1

   implicit none

!  Input arguments

   integer  , intent(in)  :: ng, nmax, nnmax1, nnmax2
   real(fpk), intent(in)  :: X(NG),XR(NG),XI(NG)

!  output arguments

   real(fpk), intent(out) :: J_BESS  (NPNG2,NPN1)
   real(fpk), intent(out) :: Y_BESS  (NPNG2,NPN1)
   real(fpk), intent(out) :: JR_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: JI_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DJ_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DY_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DJR_BESS(NPNG2,NPN1)
   real(fpk), intent(out) :: DJI_BESS(NPNG2,NPN1)

!  Local arrays

   real(fpk) :: AJ (NPN1), AY (NPN1), AJR (NPN1), AJI (NPN1)
   real(fpk) :: ADJ(NPN1), ADY(NPN1), ADJR(NPN1), ADJI(NPN1)

!  Other local variables

   integer   :: I, N
   real(fpk) :: XX,YR,YI

!  Main loop
 
   DO I = 1, NG

!  Argument

      XX = X(I)

!  J and Y Bessel functions

      CALL Bessel_rj(NMAX,NNMAX1,XX,AJ,ADJ)
      CALL Bessel_ry(NMAX,XX,AY,ADY)

!  Complex spherical Bessel functions

      YR=XR(I)
      YI=XI(I)
      CALL Bessel_Cspher(NMAX,NNMAX2,YR,YI,AJR,AJI,ADJR,ADJI)

!  Copy to output arrays (this could be skipped or tidied up)

      DO N = 1, NMAX
         J_BESS(I,N)   = AJ(N)
         Y_BESS(I,N)   = AY(N)
         JR_BESS(I,N)  = AJR(N)
         JI_BESS(I,N)  = AJI(N)
         DJ_BESS(I,N)  = ADJ(N)
         DY_BESS(I,N)  = ADY(N)
         DJR_BESS(I,N) = ADJR(N)
         DJI_BESS(I,N) = ADJI(N)
      ENDDO

!  End quadrature loop

   ENDDO

!  FInish

   return
end subroutine Bessel_functions_master

subroutine Bessel_rj    &
     ( NMAX, NNMAX, X,  & ! Inputs
       Y, U )             ! Outputs

   implicit none

!  Input arguments

   integer  , intent(in)  :: nmax, nnmax
   real(fpk), intent(in)  :: X

!  output arguments

   real(fpk), intent(out) :: U  (NMAX)
   real(fpk), intent(out) :: Y  (NMAX)

!  Local arrays

   real(fpk) :: Z(800)

!  Other local variables

   integer   :: L, L1, I1, I
   real(fpk) :: XX, Z0, Y0, Y1, YI1, YI

!  initialize

   L    = NMAX+NNMAX 
   XX   = 1.0D0/X
   Z(L) = 1.0D0/(DBLE(2*L+1)*XX)
   L1   = L - 1

!  others

   DO I = 1, L1
      I1 = L - I
      Z(I1) = 1.0D0/(DBLE(2*I1+1)*XX-Z(I1+1))
   ENDDO

!  Fill up first value

   Z0 = 1.0D0 / (XX-Z(1))
   Y0 = Z0 * DCOS(X) * XX
   Y1 = Y0 * Z(1)
   U(1) = Y0 - Y1*XX
   Y(1) = Y1

!  Remainder

   DO I = 2, NMAX
      YI1  = Y(I-1)
      YI   = YI1*Z(I)
      U(I) = YI1 - DBLE(I)*YI*XX
      Y(I) = YI
   ENDDO

!  Finish

   return
end subroutine Bessel_rj
 
subroutine Bessel_ry    &
     ( NMAX, X,         & ! Inputs
       Y, V )             ! Outputs

   implicit none

!  Input arguments

   integer  , intent(in)  :: nmax
   real(fpk), intent(in)  :: X

!  output arguments

   real(fpk), intent(out) :: V  (NMAX)
   real(fpk), intent(out) :: Y  (NMAX)

!  Local  variables

   integer   :: I, NMAX1
   real(fpk) :: C, S, X1, X2, X3, Y1

!  initialize

   C  = DCOS(X)
   S  = DSIN(X)
   X1 = 1.0D0/X
   X2 = X1*X1
   X3 = X2*X1
   Y1 = -C*X2-S*X1

   NMAX1 = NMAX-1

!  Fill y

   Y(1) = Y1
   Y(2) = (-3D0*X3+X1)*C-3D0*X2*S
   DO I = 2, NMAX1
      Y(I+1) = DBLE(2*I+1)*X1*Y(I)-Y(I-1)
   ENDDO

!  Fill V

   V(1) = -X1*(C+Y1)
   DO I = 2, NMAX
      V(I) = Y(I-1)-DBLE(I)*X1*Y(I)
   ENDDO

!  Finish

   return
end subroutine Bessel_ry
 
subroutine Bessel_Cspher      &
       ( NMAX, NNMAX, XR, XI, & ! Inputs
         YR, YI, UR, UI )       ! Outputs

!**********************************************************************
!                                                                     *
!   CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND       *
!   J=JR+I*JI OF COMPLEX ARGUMENT X=XR+I*XI OF ORDERS FROM 1 TO NMAX  *
!   BY USING BACKWARD RECURSION. PARAMETR NNMAX DETERMINES NUMERICAL  *
!   ACCURACY. U=UR+I*UI - FUNCTION (1/X)(D/DX)(X*J(X))                *
!                                                                     *
!**********************************************************************

   use tmat_parameters, only : NPN1

   implicit none

!  Input arguments

   integer  , intent(in)  :: nmax, nnmax
   real(fpk), intent(in)  :: XR, XI

!  output arguments

   real(fpk), intent(out) :: YR(NPN1),YI(NPN1)
   real(fpk), intent(out) :: UR(NPN1),UI(NPN1)

!  Local arrays

   real(fpk) :: CYR(NPN1),CYI(NPN1)
   real(fpk) :: CUR(NPN1),CUI(NPN1)
   real(fpk) :: CZR(1200),CZI(1200)

!  Other local variables

   integer   :: L, L1, I, I1
   real(fpk) :: XRXI, CXXR, CXXI, QF, QI, AR, AI, ARI
   real(fpk) :: CZ0R, CZ0I, CR, CI, CU1R, CU1I, CUIR, CUII
   real(fpk) :: CY0R, CY0I, CY1R, CY1I, CYIR, CYII, CYI1R, CYI1I

!  Code

   L=NMAX+NNMAX
   XRXI=1D0/(XR*XR+XI*XI)
   CXXR=XR*XRXI
   CXXI=-XI*XRXI 
   !QF=1D0/DFLOAT(2*L+1)
   QF=1D0/REAL(2*L+1,fpk)
   CZR(L)=XR*QF
   CZI(L)=XI*QF

   L1=L-1
   DO I=1,L1
      I1=L-I
      QF=DBLE(2*I1+1)
      AR=QF*CXXR-CZR(I1+1)
      AI=QF*CXXI-CZI(I1+1)
      ARI=1D0/(AR*AR+AI*AI)
      CZR(I1)=AR*ARI
      CZI(I1)=-AI*ARI
   ENDDO

   AR=CXXR-CZR(1)
   AI=CXXI-CZI(1)
   ARI=1D0/(AR*AR+AI*AI)
   CZ0R=AR*ARI
   CZ0I=-AI*ARI
   CR=DCOS(XR)*DCOSH(XI)
   CI=-DSIN(XR)*DSINH(XI)
   AR=CZ0R*CR-CZ0I*CI
   AI=CZ0I*CR+CZ0R*CI
   CY0R=AR*CXXR-AI*CXXI
   CY0I=AI*CXXR+AR*CXXI
   CY1R=CY0R*CZR(1)-CY0I*CZI(1)
   CY1I=CY0I*CZR(1)+CY0R*CZI(1)
   AR=CY1R*CXXR-CY1I*CXXI
   AI=CY1I*CXXR+CY1R*CXXI
   CU1R=CY0R-AR
   CU1I=CY0I-AI

   CYR(1)=CY1R
   CYI(1)=CY1I
   CUR(1)=CU1R
   CUI(1)=CU1I
   YR(1)=CY1R
   YI(1)=CY1I
   UR(1)=CU1R
   UI(1)=CU1I

   DO I=2,NMAX
      QI=DBLE(I)
      CYI1R=CYR(I-1)
      CYI1I=CYI(I-1)
      CYIR=CYI1R*CZR(I)-CYI1I*CZI(I)
      CYII=CYI1I*CZR(I)+CYI1R*CZI(I)
      AR=CYIR*CXXR-CYII*CXXI
      AI=CYII*CXXR+CYIR*CXXI
      CUIR=CYI1R-QI*AR
      CUII=CYI1I-QI*AI
      CYR(I)=CYIR
      CYI(I)=CYII
      CUR(I)=CUIR
      CUI(I)=CUII
      YR(I)=CYIR
      YI(I)=CYII
      UR(I)=CUIR
      UI(I)=CUII
   ENDDO   

!  FInish

   RETURN
end subroutine Bessel_Cspher
 
SUBROUTINE GAULEG_wrong(X1,X2,X,W,N)

!  30 December 2010. Does not work with [-1,1] and N Odd.

      implicit none

      INTEGER  , intent(in)  :: N
      REAL(fpk), intent(in)  :: X1, X2
      REAL(fpk), intent(out) :: X(N),W(N)

      INTEGER     :: I, M, J
      REAL(fpk)   :: XM,XL,P1,P2,P3,PP,Z,Z1,PIE
      REAL(fpk), PARAMETER :: EPS = 3.0D-16

      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      PIE = DACOS(-1.0d0)

      DO I=1,M
            Z=DCOS(PIE*(I-.25D0)/(N+.5D0))
!            Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
            Z1 = 0.0d0
            DO WHILE (DABS(Z-Z1).GT.EPS)
                  P1=1.D0
                  P2=0.D0
                  DO J=1,N
                        P3=P2
                        P2=P1
                        P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
                  ENDDO
                  PP=N*(Z*P1-P2)/(Z*Z-1.D0)
                  Z1=Z
                  Z=Z1-P1/PP
            ENDDO
            X(I)=XM-XL*Z
            X(N+1-I)=XM+XL*Z
            W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
            W(N+1-I)=W(I)
      ENDDO
      RETURN
END SUBROUTINE GAULEG_wrong


SUBROUTINE GAULEG_right ( N,IND1,IND2,Z,W )

      implicit none

      INTEGER  , intent(in)  :: N
      INTEGER  , intent(in)  :: IND1,IND2
      REAL(fpk), intent(out) :: Z(N),W(N)

      INTEGER        :: I, J, K, IND, M, NITER
      REAL(fpk)   :: A, B, C, F, X, CHECK, PA, PB, PC, DJ, ZZ

      A=1D0
      B=2D0
      C=3D0
      IND=MOD(N,2)
      K=N/2+IND
      F=DBLE(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.CHECK*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE

      IF(IND2.NE.0) THEN
         PRINT 1100,N
         DO 105 I=1,K
            ZZ=-Z(I)
            PRINT 1200,I,ZZ,I,W(I)
   105   CONTINUE
      ENDIF

 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER')
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)

      IF(IND1.NE.0) THEN
         DO 120 I=1,N
            Z(I)=(A+Z(I))/B
  120    CONTINUE
      ENDIF

      RETURN
END SUBROUTINE GAULEG_right

!  End module

end module tmat_functions
