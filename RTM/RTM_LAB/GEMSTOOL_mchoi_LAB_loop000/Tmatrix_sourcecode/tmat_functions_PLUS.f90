module tmat_functions_plus

  use tmat_parameters, only : fpk => tmat_fpkind, NPNG1, NPNG2, NPN1
  use tmat_functions,  only : GAULEG_right, spheroid_coords, spherical_coords, &
                              chebyshev_coords, cylinder_coords

!  LIST OF SUBROUTINES
!  ===================

!    tmatrix_constants_plus

!    tmatrix_vary_plus             ----- calling ----
!       spheroid_coords_plus
!       chebyshev_coords_plus
!       cylinder_coords_plus
!       Bessel_functions_master_plus  ----- calling ----
!          Bessel_rj_plus
!          Bessel_ry_plus
!          Bessel_Cspher_plus

!    Eqv_radius_chebyshev_plus
!    Eqv_radius_cylinder_plus
!    Eqv_radius_spheroids_plus

!  Top level routines are PUBLIC here, along with others as required

private
public        tmatrix_constants_plus, tmatrix_vary_plus,  &
              Eqv_radius_chebyshev_plus,                  &
              Eqv_radius_cylinder_plus,                   &
              Eqv_radius_spheroids_plus

contains

subroutine tmatrix_constants_plus                        &
     ( ngauss, nmax, np, aspect,                         & ! inputs
       x, w, Le_x, Le_w, an, ann, s, ss, Le_s, Le_ss )     ! outputs

!  Only linearizable thing here for Cylinders !!!

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
   REAL(fpk), intent(out) :: Le_X(NPNG2), Le_W(NPNG2)
   REAL(fpk), intent(out) :: Le_S(NPNG2), Le_SS(NPNG2)

!  Constants

   REAL(fpk), intent(out) :: AN(NPN1), ANN(NPN1,NPN1)

!  Local variables
!  ---------------

!  Help arrays

   REAL(fpk) :: X1(NPNG1),W1(NPNG1), DD(NPN1)
   REAL(fpk) :: X2(NPNG1),W2(NPNG1)

!  other variables

   INTEGER   :: N1, N, NN, NG, NG1, NG2, I
   REAL(fpk) :: D, DDD, XX, Y, YS, YSS, Le_Y, Le_YS, Le_YSS, EPS21, L_XX

!  Zero the output

   X    = 0.0d0 ; Le_X  = 0.0d0
   W    = 0.0d0 ; Le_W  = 0.0d0
   S    = 0.0d0 ; Le_S  = 0.0d0
   SS   = 0.0d0 ; Le_SS = 0.0d0
   AN   = 0.0d0 ; ANN   = 0.0d0

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
      eps21 = aspect * aspect + 1.0d0
      L_XX = dsqrt(1.0d0 - XX*XX) * aspect / eps21

!      CALL GAULEG_wrong (-1.0d0,1.0d0,X1,W1,NG1)
!      CALL GAULEG_wrong (-1.0d0,1.0d0,X2,W2,NG2)
      CALL GAULEG_right (NG1,0,0,X1,W1)
      CALL GAULEG_right (NG2,0,0,X2,W2)

      DO I = 1,NG1
         W(I)=0.5D0*(XX+1.0D0)*W1(I)
         X(I)=0.5D0*(XX+1.0D0)*X1(I)+0.5D0*(XX-1.0D0)
         Le_W(I)=0.5D0*L_XX*W1(I)
         Le_X(I)=0.5D0*L_XX*X1(I)+0.5D0*L_XX
      ENDDO

      DO I = 1,NG2
         W(I+NG1)=-0.5D0*XX*W2(I)
         X(I+NG1)=-0.5D0*XX*X2(I)+0.5D0*XX
         Le_W(I+NG1)=-0.5D0*L_XX*W2(I)
         Le_X(I+NG1)=-0.5D0*L_XX*X2(I)+0.5D0*L_XX
      ENDDO

      DO I = 1,NGAUSS
         W(NG-I+1)=W(I)
         X(NG-I+1)=-X(I)
         Le_W(NG-I+1)=Le_W(I)
         Le_X(NG-I+1)=-Le_X(I)
      ENDDO

!  Spheroids, Chebyshevs (Single quadrature, whole range)

   else

!      CALL GAULEG_wrong (-1.0d0,1.0d0,X,W,NG)
      CALL GAULEG_right (NG,0,0,X,W)

   endif

!  S vectors (same for all types)

!   DO I = 1,NGAUSS
!      Y = X(I)
!      Y = 1D0/(1D0-Y*Y)
!      SS(I)      = Y
!      SS(NG-I+1) = Y
!      Y = DSQRT(Y)
!      S(I)       = Y
!      S(NG-I+1)  = Y
!   ENDDO

!  S and SS vectors
!  ================

!  Cylinders - Fix the Linearized S, SS vectors
!    ***** New code 30 June 2011

   if ( np .eq. -2 ) then
      DO I = 1,NGAUSS
         Y     = X(I)        ; Le_Y  = Le_X(I)
         YSS = 1D0/(1D0-Y*Y) ; YS = DSQRT(YSS)
         Le_YSS = 2.0d0 * Y * Le_Y * YSS * YSS
         Le_YS  = Le_YSS * 0.5d0 / YS
         SS(I)      = YSS      ; SS(NG-I+1)    = YSS
         S (I)      = YS       ; S (NG-I+1)    = YS
         Le_SS(I)   = Le_YSS   ; Le_SS(NG-I+1) = Le_YSS
         Le_S (I)   = Le_YS    ; Le_S (NG-I+1) = Le_YS
      ENDDO
   else
      DO I = 1,NGAUSS
         Y  = X(I)  ;  YSS = 1D0/(1D0-Y*Y) ; YS = DSQRT(YSS)
         SS(I)      = YSS   ; SS(NG-I+1) = YSS
         S (I)      = YS    ; S (NG-I+1) = YS
      ENDDO
   endif

!  Finish

   return
end subroutine tmatrix_constants_plus


subroutine tmatrix_vary_plus                             &
     ( np, ngauss, nmax, n_real, n_imag, aspect,         & ! inputs
       Do_LinearEps, PI, radius, x, Le_radius, Le_x,     & ! inputs
       R, DR, DDR, DRR, DRI,                             & ! outputs
       j_bess,  y_bess,  jr_bess,  ji_bess,              & ! Outputs
       dj_bess, dy_bess, djr_bess, dji_bess,             &  ! Outputs
       Le_R, Le_DR, L_DDR, L_DRR, L_DRI, NLIN,           & ! outputs
       L_j_bess,  L_y_bess,  L_jr_bess,  L_ji_bess,      & ! Outputs
       L_dj_bess, L_dy_bess, L_djr_bess, L_dji_bess)       ! Outputs

   use tmat_parameters, only : NPNG1, NPNG2, NPN1

   implicit none

!  input arguments
!  ---------------

   logical  , intent(in)  :: Do_LinearEps
   integer  , intent(in)  :: np, ngauss, nmax
   real(fpk), intent(in)  :: n_real, n_imag
   real(fpk), intent(in)  :: Le_radius,radius, aspect, pi
   real(fpk), intent(in)  :: X(NPNG2), Le_X(NPNG2)

!  Output arguments
!  ----------------

!  number of linearizations

   integer  , intent(inout) :: nlin

!  R and DR

   real(fpk), intent(out) :: R(NPNG2), DR(NPNG2)
   real(fpk), intent(out) :: DDR(NPNG2), DRR(NPNG2), DRI(NPNG2)

   real(fpk), intent(out) :: Le_R(NPNG2), Le_DR(NPNG2)
   real(fpk), intent(out) :: L_DDR(3,NPNG2), L_DRR(3,NPNG2), L_DRI(3,NPNG2)

!  Bessel Master output

   real(fpk), intent(out) :: J_BESS  (NPNG2,NPN1)
   real(fpk), intent(out) :: Y_BESS  (NPNG2,NPN1)
   real(fpk), intent(out) :: JR_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: JI_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DJ_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DY_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DJR_BESS(NPNG2,NPN1)
   real(fpk), intent(out) :: DJI_BESS(NPNG2,NPN1)

   real(fpk), intent(out) :: L_J_BESS  (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_Y_BESS  (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_JR_BESS (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_JI_BESS (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_DJ_BESS (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_DY_BESS (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_DJR_BESS(3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_DJI_BESS(3,NPNG2,NPN1)

!  Local variables
!  ---------------

!  Local arrays

   real(fpk) :: Z(NPNG2),ZR(NPNG2),ZI(NPNG2)
   real(fpk) :: L_Z(3,NPNG2),L_ZR(3,NPNG2),L_ZI(3,NPNG2)

!  Other lcoal variables

   integer   :: i, ng, NNMAX1, NNMAX2
   real(fpk) :: Modulus, PRR, PRI, Ln_PRR(2),Ln_PRI(2)
   real(fpk) :: TA, TB, V, VV, Vinv, V1, V2
   real(fpk) :: Le_V, Le_VV, Le_Vinv

!  Constants
!  ---------

   NG = NGAUSS*2

   Modulus  = n_real*n_real + n_imag*n_imag
   PRR = + n_real / Modulus
   PRI=  - n_imag / Modulus

   Ln_PRR(1) =   ( 1.0d0 - 2.0d0 * n_real * PRR ) * n_real / Modulus
   Ln_PRR(2) =   (       - 2.0d0 * n_imag * PRR ) * n_imag / Modulus
   Ln_PRI(1) = - (       + 2.0d0 * n_real * PRI ) * n_real / Modulus
   Ln_PRI(2) = - ( 1.0d0 + 2.0d0 * n_imag * PRI ) * n_imag / Modulus

!  derivative counting

   NLIN = 2
   if ( np .ne. 0 ) then
      if ( do_LinearEps ) nlin = 3
   endif

!  Coordinates
!  -----------

   if ( np .eq. -1 ) then       ! Spheroids
      if ( do_LinearEps ) then
         call  spheroid_coords_plus                    &
          ( ng, ngauss, radius, Le_radius, aspect, x,  & ! inputs
            r, dr, Le_r, Le_dr )                         ! Outputs
      else
         Le_r = 0.0d0 ; Le_dr = 0.0d0
         call  spheroid_coords                    &
          ( ng, ngauss, radius, aspect, x,        & ! inputs
            r, dr )                                 ! Outputs
      endif
   else if ( np .eq. -2 ) then  !  Cylinders
      if ( do_LinearEps ) then
         call  cylinder_coords_plus                    &
          ( ng, ngauss, radius, Le_radius, aspect, x, Le_x, & ! inputs
            r, dr, Le_r, Le_dr )                              ! Outputs
      else
         Le_r = 0.0d0 ; Le_dr = 0.0d0
         call  cylinder_coords                   &
          ( ng, ngauss, radius, aspect, x,       & ! inputs
            r, dr )                                ! Outputs
      endif
   else if ( np .gt. 0 ) then   !  Chebyshev
      if ( do_LinearEps ) then
         call  chebyshev_coords_plus               &
          ( ng, radius, Le_radius, aspect, np, x,  & ! inputs
            r, dr, Le_r, Le_dr )                     ! Outputs
      else
         Le_r = 0.0d0 ; Le_dr = 0.0d0
         call  chebyshev_coords                  &
          ( ng, radius, aspect, np, x,           & ! inputs
            r, dr )                                ! Outputs
      endif
   else if ( np .eq. 0 ) then   !  Mie
      Le_r = 0.0d0 ; Le_dr = 0.0d0
      call  spherical_coords ( ng, radius, r, dr ) 
   endif

!  loop for DDR, DRR, DRI, ZR, ZI
!  ------------------------------

!  Initialise derivatives

   L_DDR = 0.0d0 ; L_DRR = 0.0d0 ; L_DRI = 0.0d0  
   L_Z   = 0.0d0 ; L_ZR  = 0.0d0 ; L_ZI  = 0.0d0

!  Start
  
   TA=0D0
   DO I = 1, NG

!  Functions

        VV=DSQRT(R(I))
        V=VV*PI
        TA=MAX(TA,V)
        VInv = 1D0/V
        V1 = V * n_real
        V2 = V * n_imag
        Z(I)  = V
        ZR(I) = V1
        ZI(I) = V2
        DDR(I) = Vinv
        DRR(I) = PRR*Vinv
        DRI(I) = PRI*Vinv

!  Linearizations w.r.t n_real/n_imag

        L_ZR(1,I) = V1
        L_ZI(2,I) = V2
        L_DRR(1,I) = Ln_PRR(1)*Vinv
        L_DRR(2,I) = Ln_PRR(2)*Vinv
        L_DRI(1,I) = Ln_PRI(1)*Vinv
        L_DRI(2,I) = Ln_PRI(2)*Vinv

!  Linearization w.r.t eps

        if ( do_LinearEps ) then
           Le_VV   = 0.5d0 * Le_R(I) / VV
           Le_V    = Le_VV * PI
           Le_Vinv = - Vinv * Vinv * Le_V
           L_DDR(3,I) = Le_Vinv
           L_DRR(3,I) = PRR * Le_Vinv
           L_DRI(3,I) = PRI * Le_Vinv
           L_Z(3,I)   = Le_V
           L_ZR(3,I)  = Le_V * n_real
           L_ZI(3,I)  = Le_V * n_imag
        endif

   ENDDO

!  Maximum numbers for exapnsions

   TB = TA*DSQRT(Modulus)
   TB = DMAX1 ( TB,DBLE(NMAX) )

   NNMAX1 = 1.2D0*DSQRT(DMAX1(TA,DBLE(NMAX)))+3D0
   NNMAX2 = (TB+4D0*(TB**0.33333D0)+1.2D0*DSQRT(TB))
   NNMAX2 = NNMAX2-NMAX+5

!  Bessel functions

   call Bessel_functions_master_plus                   &
       ( ng, nmax, nnmax1, nnmax2, nlin,                  & ! Inputs
         z, zr, zi, L_z, L_zr, L_zi,                      & ! Inputs
         j_bess,  y_bess,  jr_bess,  ji_bess,             & ! Outputs
        dj_bess, dy_bess, djr_bess, dji_bess,             & ! Outputs
         L_j_bess,  L_y_bess,  L_jr_bess,  L_ji_bess,     & ! Outputs
        L_dj_bess, L_dy_bess, L_djr_bess, L_dji_bess)       ! Outputs

!  Finish

   return
end subroutine tmatrix_vary_plus


subroutine Bessel_functions_master_plus                   &
       ( ng, nmax, nnmax1, nnmax2, nlin,                  & ! Inputs
         x, xr, xi, L_x, L_xr, L_xi,                      & ! Inputs
         j_bess,  y_bess,  jr_bess,  ji_bess,             & ! Outputs
        dj_bess, dy_bess, djr_bess, dji_bess,             & ! Outputs
         L_j_bess,  L_y_bess,  L_jr_bess,  L_ji_bess,     & ! Outputs
        L_dj_bess, L_dy_bess, L_djr_bess, L_dji_bess)       ! Outputs

   use tmat_parameters, only : NPNG2, NPN1

   implicit none

!  Input arguments

   integer  , intent(in)  :: ng, nmax, nnmax1, nnmax2,nlin
   real(fpk), intent(in)  :: X(NG),XR(NG),XI(NG)
   real(fpk), intent(in)  :: L_X(3,NG),L_XR(3,NG),L_XI(3,NG)

!  output arguments

   real(fpk), intent(out) :: J_BESS  (NPNG2,NPN1)
   real(fpk), intent(out) :: Y_BESS  (NPNG2,NPN1)
   real(fpk), intent(out) :: JR_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: JI_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DJ_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DY_BESS (NPNG2,NPN1)
   real(fpk), intent(out) :: DJR_BESS(NPNG2,NPN1)
   real(fpk), intent(out) :: DJI_BESS(NPNG2,NPN1)

   real(fpk), intent(out) :: L_J_BESS  (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_Y_BESS  (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_JR_BESS (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_JI_BESS (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_DJ_BESS (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_DY_BESS (3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_DJR_BESS(3,NPNG2,NPN1)
   real(fpk), intent(out) :: L_DJI_BESS(3,NPNG2,NPN1)

!  Local arrays

   real(fpk) :: AJ (NPN1), AY (NPN1), AJR (NPN1), AJI (NPN1)
   real(fpk) :: ADJ(NPN1), ADY(NPN1), ADJR(NPN1), ADJI(NPN1)
   real(fpk) :: L_AJ  (3,NPN1),  L_AY  (3,NPN1)
   real(fpk) :: L_AJR (3,NPN1),  L_AJI (3,NPN1)
   real(fpk) :: L_ADJ (3,NPN1),  L_ADY (3,NPN1)
   real(fpk) :: L_ADJR(3,NPN1),  L_ADJI(3,NPN1)

!  Other local variables

   integer   :: I, N, Q
   real(fpk) :: XX,YR,YI,L_XX(3),L_YR(3),L_YI(3)

!  Main loop
 
   DO I = 1, NG

!  Arguments

      XX = X(I)  ; L_XX(1:NLIN) = L_X(1:NLIN,I)
      YR = XR(I) ; L_YR(1:NLIN) = L_XR(1:NLIN,I)
      YI = XI(I) ; L_YI(1:NLIN) = L_XI(1:NLIN,I)

!  J and Y Bessel functions

      call Bessel_rj_plus             &
     ( NMAX, NNMAX1, NLIN, XX, L_XX,  & ! Inputs
       AJ, ADJ, L_AJ, L_ADJ )           ! Outputs

      call Bessel_ry_plus         &
     ( NMAX, NLIN, XX, L_XX,      & ! Inputs
       AY, ADY, L_AY, L_ADY )       ! Outputs

!  Complex spherical Bessel functions

      call Bessel_Cspher_plus                   &
     ( NMAX, NNMAX2, NLIN, YR, YI, L_YR, L_YI,  & ! Inputs
       AJR, AJI, ADJR, ADJI,                    & ! Outputs
       L_AJR, L_AJI, L_ADJR, L_ADJI )             ! Outputs

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
         DO Q = 1, NLIN
            L_J_BESS(Q,I,N)   = L_AJ(Q,N)
            L_Y_BESS(Q,I,N)   = L_AY(Q,N)
            L_JR_BESS(Q,I,N)  = L_AJR(Q,N)
            L_JI_BESS(Q,I,N)  = L_AJI(Q,N)
            L_DJ_BESS(Q,I,N)  = L_ADJ(Q,N)
            L_DY_BESS(Q,I,N)  = L_ADY(Q,N)
            L_DJR_BESS(Q,I,N) = L_ADJR(Q,N)
            L_DJI_BESS(Q,I,N) = L_ADJI(Q,N)
         ENDDO
      ENDDO

!  End quadrature loop

   ENDDO

!  FInish

   return
end subroutine Bessel_functions_master_plus


subroutine spheroid_coords_plus                &
     ( ng, ngauss, radius, Le_radius, eps, x,  & ! inputs
       r, dr, Le_r, Le_dr )                      ! Outputs

!  Coordinates for the spheroidal particles

!    Includes the EPS derivatives of R and DR as output
!    Input radius may have an EPS linearization (Equiv. Surf. Area)
!    Input x-values do not depend on EPS linearization

   implicit none

!  inputs

   integer  , intent(in)  :: ng, ngauss
   real(fpk), intent(in)  :: radius, Le_radius, eps
   real(fpk), intent(in)  :: x(ng)

!  outputs

   real(fpk), intent(out) :: r (ng)
   real(fpk), intent(out) :: dr(ng)
   real(fpk), intent(out) :: Le_r (ng)
   real(fpk), intent(out) :: Le_dr(ng)

!  Local

   integer   :: i, i1
   real(fpk) :: third, a, aa, ee, ee1, eps3, c, cc, ss, s, cs, rr
   real(fpk) :: L_a, L_aa, L_ee, L_ee1, L_eps3, L_rr

!  initialize

   R    = 0.0d0   ; DR    = 0.0d0
   Le_R = 0.0d0   ; Le_DR = 0.0d0

!  Constants for this rotuine

   THIRD = (1.0D0/3.0D0)
   EPS3  = EPS ** THIRD

   A   = radius * EPS3
   AA  = A * A
   EE  = EPS * EPS
   EE1 = EE - 1.0D0

!  Linearized constants. Normalized derivative throughout !!!!

   L_EPS3 = THIRD  * EPS3
   L_A    = radius * L_EPS3
   if ( Le_radius .ne. 0.0d0 ) L_A = L_A + Le_radius * eps3
   L_AA  = 2.0d0 * A * L_A
   L_EE  = 2.0d0 * EE
   L_EE1 = L_EE

!  Loop over quadrature radii

   DO I = 1, NGAUSS
      I1 = NG-I+1
      C  = X(I)
      CC = C*C
      SS = 1.0D0-CC
      S  = DSQRT(SS)
      CS = C * S
      RR = 1.0D0/(SS+EE*CC)
      L_RR = - RR * RR * CC * L_EE
      R (I) = AA*RR
      DR(I) = RR*CS*EE1
      Le_R (I) = L_AA * RR + AA * L_RR
      Le_DR(I) = CS * ( L_RR * EE1 + RR * L_EE1 )
      R   (I1) = R   (I)  ; DR   (I1) = - DR   (I)
      Le_R(I1) = Le_R(I)  ; Le_DR(I1) = - Le_DR(I)
   ENDDO

!  Finish

   return
end subroutine spheroid_coords_plus


subroutine chebyshev_coords_plus           &
     ( ng, radius, Le_radius, eps, np, x,  & ! inputs
       r, dr, Le_r, Le_dr )                  ! Outputs

!  Coordinates for the Chebyshev particles

!    Includes the EPS derivatives of R and DR as output
!    Input radius may have an EPS linearization (Equiv. Surf. Area)
!    Input x-values do not depend on EPS linearization

   implicit none

!  inputs

   integer  , intent(in)  :: ng, np
   real(fpk), intent(in)  :: radius, Le_radius, eps
   real(fpk), intent(in)  :: x(ng)

!  outputs

   real(fpk), intent(out) :: r (ng)
   real(fpk), intent(out) :: dr(ng)
   real(fpk), intent(out) :: Le_r (ng)
   real(fpk), intent(out) :: Le_dr(ng)

!  Local

   integer   :: i
   real(fpk) :: a, dn, dnp, dn4, ep, e3, r0, xi, ri, fac, third
   real(fpk) :: l_a, l_ep, s1, s2, t1, L_r0, L_ri, L_t1, dci, dsi

!  initialize

   R = 0.0d0   ; DR = 0.0d0
   Le_R = 0.0d0   ; Le_DR = 0.0d0

!  Constants

   DNP = DBLE(NP)
   DN  = DNP*DNP
   DN4 = DN*4.0D0
   EP  = EPS*EPS
   FAC = 1.5d0 * (DN4-2.0D0)/(DN4-1.0D0)
   A = 1.0D0 + EP * FAC

   L_EP = 2.0d0 * EP
   L_A  = L_EP * FAC

!  Is this really an integer ?

!   I = ( ( DNP + 0.1D0 ) * 0.5D0 )                ! Original code
   I = INT ( ( DNP + 0.1D0 ) * 0.5D0 )
   I = 2*I
   IF (I.EQ.NP) THEN
      E3 = EP * EPS
      S1 = 3.0D0 / (DN-1.0D0)
      S2 = 0.25D0/(9.0D0*DN-1.0D0)
      A   = A - S1 * EPS*(1.0D0+0.25D0*EP) - S2 * E3
      L_A = L_A - S1 * EPS * (1.0D0+0.75D0*EP) - 3.0d0 * S2 * E3
   ENDIF

   third = (-1.0D0/3.0D0)
   R0   = radius * A ** third
   L_R0 = R0 * ( (Le_radius/radius) + third * ( L_A / A ) )

!  Quadrature loop

   DO I = 1, NG
      XI = DACOS(X(I))*DNP
      DCI = DCOS(XI)
      DSI = DSIN(XI)
      S1 = 1D0 + EPS * DCI
      RI = R0 * S1
      L_RI = L_R0 * S1 + R0 * EPS * DCI
      S2 = DNP*DSI
      T1 = R0 * EPS * S2
      L_T1 = S2 * EPS * ( L_R0 + R0 )
      R(I)  = RI*RI
      Le_R(I) = 2.0 * RI * L_RI
      DR(I) = - T1 / RI
      Le_DR(I) = DR(I)  * ( (L_T1/T1) - (L_RI/RI) )
   ENDDO

!  Finish

   return
end subroutine chebyshev_coords_plus


subroutine spherical_coords_plus       &
     ( ng, radius, Le_radius,          & ! inputs
       r, dr, Le_r, Le_dr )              ! Outputs

!  Coordinates for Spherical particles (Mie code)

   implicit none

!  inputs

   integer  , intent(in)  :: ng
   real(fpk), intent(in)  :: radius, Le_radius

!  outputs

   real(fpk), intent(out) :: r (ng)
   real(fpk), intent(out) :: dr(ng)
   real(fpk), intent(out) :: Le_r (ng)
   real(fpk), intent(out) :: Le_dr(ng)

!  Local

   integer   :: i
   real(fpk) :: r0, L_r0

   R    = 0.0d0   ; DR    = 0.0d0
   Le_R = 0.0d0   ; Le_DR = 0.0d0

   R0 = Radius * Radius
   L_R0 = Le_Radius * Radius * 2.0d0

   DO I = 1, NG
      R(I)     = R0
      Le_R(I)  = L_R0
   ENDDO

!  Finish

   return
end subroutine spherical_coords_plus



subroutine cylinder_coords_plus                      &
     ( ng, ngauss, radius, Le_radius, eps, x, Le_x,  & ! inputs
       r, dr, Le_r, Le_dr )                            ! Outputs

!  Coordinates for the Cylindrical particles

!    Includes the EPS derivatives of R and DR as output
!    Input radius may have an EPS linearization (Equiv. Surf. Area)
!    Input x-values do depend on EPS linearization

   implicit none

!  inputs

   integer  , intent(in)  :: ng, ngauss
   real(fpk), intent(in)  :: radius, Le_radius, eps
   real(fpk), intent(in)  :: x(ng), Le_x(ng)

!  outputs

   real(fpk), intent(out) :: r (ng)
   real(fpk), intent(out) :: dr(ng)
   real(fpk), intent(out) :: Le_r (ng)
   real(fpk), intent(out) :: Le_dr(ng)

!  Local

   integer   :: i, i1
   real(fpk) :: a, h, co, si, rad, rthet, c13, c23, e23, e23p
   real(fpk) :: L_a, L_h, L_co, L_si, aco, hsi, L_rad, L_rthet

!  initialize

   R = 0.0d0    ; DR = 0.0d0
   Le_R = 0.0d0 ; Le_DR = 0.0d0

!  Constants

   C13 = 1.0d0 / 3.0d0
   C23 = 2.0d0 * C13
   E23 = C23 / EPS / EPS
   E23P = E23 ** C13

   H = Radius * E23P
   A = H*EPS
   L_H = ( Le_Radius - C23 * radius )  * E23P
   L_A = EPS * ( H + L_H )

!  Quadrature loop

   DO I = 1, NGAUSS
      I1 = NG-I+1
      CO = -X(I)
      SI = DSQRT(1.0D0-CO*CO)
      L_CO = - Le_x(I)
      L_SI = + CO * Le_x(I) / SI
      IF (SI/CO.GT.A/H) THEN
         ACO = A * CO
         RAD   = A/SI
         RTHET = -ACO/SI/SI
         L_RAD   = ( L_A - RAD * L_SI ) / SI
         L_RTHET = RTHET * ( ((L_A*CO+A*L_CO)/ACO) - (2.0*L_SI/SI) )
      ELSE
         HSI = H * SI
         RAD   = H/CO
         RTHET = HSI/(CO*CO)
         L_RAD = ( L_H - RAD * L_CO ) / CO
         L_RTHET = RTHET * ( ((L_H*SI+H*L_SI)/HSI) - (2.0*L_CO/CO) )
      ENDIF
      R(I)  =  RAD*RAD       ; R(I1)  =   R(I)
      DR(I) = -RTHET/RAD     ; DR(I1) = - DR(I)
      Le_R(I)  =  2.0d0*L_RAD*RAD                ; Le_R(I1)  =   Le_R(I)
      Le_DR(I) = - (L_RTHET+DR(I)*L_RAD)/RAD     ; Le_DR(I1) = - Le_DR(I)
   ENDDO

!  Finish

   return
end subroutine cylinder_coords_plus


subroutine Bessel_rj_plus          &
     ( NMAX, NNMAX, NLIN, X, L_X,  & ! Inputs
       Y, U, L_Y, L_U )              ! Outputs

   implicit none

!  Input arguments

   integer  , intent(in)  :: nmax, nnmax, nlin
   real(fpk), intent(in)  :: X
   real(fpk), intent(in)  :: L_X(3)

!  output arguments

   real(fpk), intent(out) :: U  (NMAX)
   real(fpk), intent(out) :: Y  (NMAX)
   real(fpk), intent(out) :: L_U  (3,NMAX)
   real(fpk), intent(out) :: L_Y  (3,NMAX)

!  Local arrays

   real(fpk) :: Z(800)
   real(fpk) :: L_Z(3,800)

!  Other local variables

   integer   :: L, L1, I1, I, Q
   real(fpk) :: XX, Z0, Y0, Y1, YI1, YI, DL, DI, MZSQ, DCX, DSX
   real(fpk) :: L_XX(3), L_Z0, L_Y0, L_Y1, L_YI1, L_YI

!  initialize

   L    = NMAX+NNMAX 
   L1   = L - 1
   DL   = 1.0d0/DBLE(2*L+1)
   XX   = 1.0D0/X
   Z(L) = X * DL
   DO Q = 1, NLIN
     L_Z(Q,L) = L_X(Q) * DL
     L_XX(Q)  = - XX * XX * L_X(Q)
   ENDDO

!  others

   DO I = 1, L1
      I1 = L - I
      DI = DBLE(2*I1+1)
      Z(I1) = 1.0D0/(DI*XX-Z(I1+1))
      MZSQ = - Z(I1) * Z(I1)
      DO Q = 1, NLIN
         L_Z(Q,I1) = MZSQ * ( DI * L_XX(Q) - L_Z(Q,I1+1) )
      ENDDO
   ENDDO

!  Fill up first value

   Z0 = 1.0D0 / (XX-Z(1))
   DCX = DCOS(X)
   DSX = DSIN(X)
   MZSQ = - Z0 * Z0

   Y0 = Z0 * DCX * XX
   Y1 = Y0 * Z(1)
   U(1) = Y0 - Y1*XX
   Y(1) = Y1

   DO Q = 1, NLIN
      L_Z0 = MZSQ * ( L_XX(Q) - L_Z(Q,1) )
      L_Y0 = L_Z0 * DCX * XX + Z0 * ( - DSX * L_X(Q) * XX + DCX * L_XX(Q) )
      L_Y1 = Y0 * L_Z(Q,1) + L_Y0 * Z(1)
      L_U(Q,1) = L_Y0 - L_Y1 * XX - Y1 * L_XX(Q)
      L_Y(Q,1) = L_Y1
   ENDDO

!  Remainder

   DO I = 2, NMAX
      DI = DBLE(I)
      YI1  = Y(I-1)
      YI   = YI1*Z(I)
      U(I) = YI1 - DI * YI * XX
      Y(I) = YI
      DO Q = 1, NLIN
         L_YI1 = L_Y(Q,I-1)
         L_YI  = L_YI1 * Z(I) + YI1 * L_Z(Q,I)
         L_U(Q,I) = L_YI1 - DI * ( L_YI * XX + YI * L_XX(Q) )
         L_Y(Q,I) = L_YI
      ENDDO
   ENDDO

!  Finish

   return
end subroutine Bessel_rj_plus
 

subroutine Bessel_ry_plus          &
     ( NMAX, NLIN, X, L_X,         & ! Inputs
       Y, V, L_Y, L_V )              ! Outputs

   implicit none

!  Input arguments

   integer  , intent(in)  :: nmax, nlin
   real(fpk), intent(in)  :: X
   real(fpk), intent(in)  :: L_X(3)

!  output arguments

   real(fpk), intent(out) :: V  (NMAX)
   real(fpk), intent(out) :: Y  (NMAX)
   real(fpk), intent(out) :: L_V  (3,NMAX)
   real(fpk), intent(out) :: L_Y  (3,NMAX)

!  Local  variables

   integer   :: I, NMAX1, Q
   real(fpk) :: C, S, X1, X2, X3, Y1, V1C, Y2C, Y2S, DI
   real(fpk) :: L_C(3), L_S(3), L_X1(3), L_X2(3), L_X3(3), L_Y1(3)

!  initialize

   C  = DCOS(X)
   S  = DSIN(X)
   X1 = 1.0D0/X
   X2 = X1*X1
   X3 = X2*X1
   Y1 = -C*X2-S*X1
   DO Q = 1, NLIN
      L_C(Q) = - S * L_X(Q)
      L_S(Q) = + C * L_X(Q)
      L_X1(Q) = - X2 * L_X(Q)
      L_X2(Q) = 2.0d0 * X1 * L_X1(Q)
      L_X3(Q) = 3.0d0 * X2 * L_X1(Q)
      L_Y1(Q) = - ( L_C(Q)*X2 + C*L_X2(Q) + L_S(Q)*X1 + S*L_X1(Q) )
   ENDDO 

   NMAX1 = NMAX-1

!  Fill y

   Y(1) = Y1
   Y2C  = -3D0 * X3 + X1
   Y2S  =  3D0 * X2
   Y(2) = Y2C*C - Y2S*S
   DO Q = 1, NLIN
      L_Y(Q,1) = L_Y1(Q)
      L_Y(Q,2) = L_C(Q) * Y2C + C * ( L_X1(Q)-3D0*L_X3(Q) ) - &
                 L_S(Q) * Y2S - S * 3D0 * L_X2(Q)
   ENDDO

   DO I = 2, NMAX1
      DI = DBLE(2*I+1)
      Y(I+1) = DI*X1*Y(I)-Y(I-1)
      DO Q = 1, NLIN
         L_Y(Q,I+1) = DI * ( L_X1(Q)*Y(I) + X1*L_Y(Q,I) ) - L_Y(Q,I-1)
      ENDDO
   ENDDO

!  Fill V

   V1C  = C + Y1
   V(1) = - X1 * V1C
   DO Q = 1, NLIN
      L_V(Q,1) = - L_X1(Q) * V1C - X1 * ( L_C(Q) + L_Y1(Q) )
   ENDDO

   DO I = 2, NMAX
      DI = DBLE(I)
      V(I) = Y(I-1) - DI * X1 * Y(I)
      DO Q = 1, NLIN
         L_V(Q,I) = L_Y(Q,I-1) - DI * ( L_X1(Q)*Y(I) + X1*L_Y(Q,I) ) 
      ENDDO
   ENDDO

!  Finish

   return
end subroutine Bessel_ry_plus
 
 
subroutine Bessel_Cspher_plus                  &
     ( NMAX, NNMAX, NLIN, XR, XI, L_XR, L_XI,  & ! Inputs
       YR, YI, UR, UI, L_YR, L_YI, L_UR, L_UI )  ! Outputs

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

   integer  , intent(in)  :: nmax, nnmax, nlin
   real(fpk), intent(in)  :: XR, XI
   real(fpk), intent(in)  :: L_XR(3), L_XI(3)

!  output arguments

   real(fpk), intent(out) :: YR(NPN1),YI(NPN1)
   real(fpk), intent(out) :: UR(NPN1),UI(NPN1)
   real(fpk), intent(out) :: L_YR(3,NPN1),L_YI(3,NPN1)
   real(fpk), intent(out) :: L_UR(3,NPN1),L_UI(3,NPN1)

!  Local arrays

   real(fpk) :: CYR(NPN1),CYI(NPN1)
!   real(fpk) :: CUR(NPN1),CUI(NPN1)          !  Why do you need these ??
   real(fpk) :: L_CYR(3,NPN1),L_CYI(3,NPN1)
!   real(fpk) :: L_CUR(3,NPN1),L_CUI(3,NPN1)          !  Why do you need these ??

   real(fpk) :: CZR(1200),CZI(1200)
   real(fpk) :: L_CZR(3,1200),L_CZI(3,1200)

!  Other local variables

   integer   :: L, L1, I, I1, Q
   real(fpk) :: XRXI, CXXR, CXXI, QF, QI, AR, AI, ARI, MSQ
   real(fpk) :: AR0, AI0, AR2, AI2, AR3, AI3
   real(fpk) :: COSXR, SINXR, COSHXI, SINHXI
   real(fpk) :: CZ0R, CZ0I, CR, CI, CU1R, CU1I, CUIR, CUII
   real(fpk) :: CY0R, CY0I, CY1R, CY1I, CYIR, CYII, CYI1R, CYI1I

   real(fpk) :: L_XRXI, L_CXXR(3), L_CXXI(3), L_AR, L_AI, L_ARI, L_MSQ
   real(fpk) :: L_AR0, L_AI0, L_AR2, L_AI2, L_AR3, L_AI3
   real(fpk) :: L_COSXR, L_SINXR, L_COSHXI, L_SINHXI
   real(fpk) :: L_CZ0R, L_CZ0I, L_CR, L_CI, L_CU1R, L_CU1I, L_CUIR, L_CUII
   real(fpk) :: L_CY0R, L_CY0I, L_CY1R, L_CY1I, L_CYIR, L_CYII, L_CYI1R, L_CYI1I

!  Initialize CZ

   L=NMAX+NNMAX

   MSQ  = XR*XR+XI*XI
   XRXI = 1D0/MSQ
   CXXR = XR*XRXI
   CXXI = -XI*XRXI 
   QF = 1D0/DBLE(2*L+1)
   CZR(L) = XR*QF
   CZI(L) = XI*QF

!  Initialize L_CZ

   DO Q = 1, NLIN
     L_MSQ  = 2D0 * ( XR * L_XR(Q) + XI * L_XI(Q) )
     L_XRXI = - XRXI * XRXI * L_MSQ
     L_CXXR(Q)  =   XR * L_XRXI + L_XR(Q) * XRXI
     L_CXXI(Q)  = - XI * L_XRXI - L_XI(Q) * XRXI
     L_CZR(Q,L) = L_XR(Q) * QF
     L_CZI(Q,L) = L_XI(Q) * QF
   ENDDO

!  Downward recursion for CZ and L_CZ

   L1=L-1
   DO I=1,L1
      I1=L-I
      QF=DBLE(2*I1+1)
      AR = QF*CXXR-CZR(I1+1)
      AI = QF*CXXI-CZI(I1+1)
      MSQ = AR*AR+AI*AI
      ARI = 1D0 / MSQ
      CZR(I1) =  AR*ARI
      CZI(I1) = -AI*ARI
      DO Q = 1, NLIN
         L_AR = QF*L_CXXR(Q) - L_CZR(Q,I1+1)
         L_AI = QF*L_CXXI(Q) - L_CZI(Q,I1+1)
         L_MSQ  = 2D0 * ( AR*L_AR+AI*L_AI )
         L_ARI = - ARI * ARI * L_MSQ
         L_CZR(Q,I1) =   L_AR*ARI + AR * L_ARI
         L_CZI(Q,I1) = - L_AI*ARI - AI * L_ARI
      ENDDO
   ENDDO

!  More setups and their linearizations

   AR0 = CXXR-CZR(1)
   AI0 = CXXI-CZI(1)
   MSQ = AR0*AR0+AI0*AI0
   ARI = 1D0 / MSQ

   CZ0R =  AR0*ARI
   CZ0I = -AI0*ARI

   COSXR  = DCOS(XR)
   SINXR  = DSIN(XR)
   COSHXI = DCOSH(XI)
   SINHXI = DSINH(XI)

   CR = + COSXR * COSHXI
   CI = - SINXR * SINHXI

   AR2=CZ0R*CR-CZ0I*CI
   AI2=CZ0I*CR+CZ0R*CI

   CY0R=AR2*CXXR-AI2*CXXI
   CY0I=AI2*CXXR+AR2*CXXI

   CY1R=CY0R*CZR(1)-CY0I*CZI(1)
   CY1I=CY0I*CZR(1)+CY0R*CZI(1)

   AR3=CY1R*CXXR-CY1I*CXXI
   AI3=CY1I*CXXR+CY1R*CXXI

   CU1R=CY0R-AR3
   CU1I=CY0I-AI3

   CYR(1)=CY1R     !  ;  CUR(1)=CU1R
   CYI(1)=CY1I     !  ;  CUI(1)=CU1I
   YR(1)=CY1R
   YI(1)=CY1I
   UR(1)=CU1R
   UI(1)=CU1I

!  Linearizations of initial function values

   DO Q = 1, NLIN

      L_AR0 = L_CXXR(Q) - L_CZR(Q,1)
      L_AI0 = L_CXXI(Q) - L_CZI(Q,1)
      L_MSQ = 2.0d0 * ( L_AR0*AR0 + L_AI0*AI0 )
      L_ARI = - ARI * ARI * L_MSQ

      L_CZ0R =   L_AR0*ARI + AR0*L_ARI
      L_CZ0I = - L_AI0*ARI - AI0*L_ARI
      
      L_COSXR  = - SINXR * L_XR(Q)
      L_SINXR  = + COSXR * L_XR(Q)
      L_COSHXI = + SINHXI * L_XI(Q)
      L_SINHXI = + COSHXI * L_XI(Q)

      L_CR = + L_COSXR * COSHXI + COSXR * L_COSHXI
      L_CI = - L_SINXR * SINHXI - SINXR * L_SINHXI

      L_AR2 = L_CZ0R * CR - L_CZ0I * CI + CZ0R * L_CR - CZ0I * L_CI
      L_AI2 = L_CZ0I * CR + L_CZ0R * CI + CZ0I * L_CR + CZ0R * L_CI

      L_CY0R = L_AR2 * CXXR - L_AI2 * CXXI + AR2 * L_CXXR(Q) - AI2 * L_CXXI(Q)
      L_CY0I = L_AI2 * CXXR + L_AR2 * CXXI + AI2 * L_CXXR(Q) + AR2 * L_CXXI(Q)

      L_CY1R = L_CY0R *   CZR(1)   - L_CY0I *   CZI(1) + &
                 CY0R * L_CZR(Q,1) -   CY0I * L_CZI(Q,1)
      L_CY1I = L_CY0I *   CZR(1)   + L_CY0R *   CZI(1) + &
                 CY0I * L_CZR(Q,1) +   CY0R * L_CZI(Q,1)

      L_AR3 = L_CY1R * CXXR - L_CY1I * CXXI + CY1R * L_CXXR(Q) - CY1I * L_CXXI(Q)
      L_AI3 = L_CY1I * CXXR + L_CY1R * CXXI + CY1I * L_CXXR(Q) + CY1R * L_CXXI(Q)

      L_CU1R = L_CY0R - L_AR3
      L_CU1I = L_CY0I - L_AI3

      L_YR(Q,1) = L_CY1R
      L_YI(Q,1) = L_CY1I
      L_UR(Q,1) = L_CU1R
      L_UI(Q,1) = L_CU1I    
      L_CYR(Q,1) = L_CY1R
      L_CYI(Q,1) = L_CY1I
!      L_CUR(Q,1) = L_CU1R
!      L_CUI(Q,1) = L_CU1I    

   ENDDO

!  Recursion upwards for Y and U, and their linearizations

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
      YR(I)=CYIR
      YI(I)=CYII
      UR(I)=CUIR
      UI(I)=CUII
      CYR(I)=CYIR       ! ; CUR(I)=CUIR
      CYI(I)=CYII       ! ; CUI(I)=CUII
      DO Q = 1, NLIN
         L_CYI1R = L_CYR(Q,I-1)
         L_CYI1I = L_CYI(Q,I-1)
         L_CYIR  = L_CYI1R *   CZR(I)   - L_CYI1I *   CZI(I) + & 
                     CYI1R * L_CZR(Q,I) -   CYI1I * L_CZI(Q,I)
         L_CYII  = L_CYI1I *   CZR(I)   + L_CYI1R *   CZI(I) + &
                     CYI1I * L_CZR(Q,I) +   CYI1R * L_CZI(Q,I)
         L_AR = L_CYIR *   CXXR    - L_CYII *   CXXI + &
                  CYIR * L_CXXR(Q) -   CYII * L_CXXI(Q)
         L_AI = L_CYII *   CXXR    + L_CYIR *   CXXI + &
                  CYII * L_CXXR(Q) +   CYIR * L_CXXI(Q)
         L_CUIR = L_CYI1R - QI * L_AR
         L_CUII = L_CYI1I - QI * L_AI
         L_CYR(Q,I) = L_CYIR      !   ;      L_CUR(Q,I) = L_CUIR
         L_CYI(Q,I) = L_CYII      !   ;      L_CUI(Q,I) = L_CUII
         L_YR(Q,I) = L_CYIR
         L_YI(Q,I) = L_CYII
         L_UR(Q,I) = L_CUIR
         L_UI(Q,I) = L_CUII
      ENDDO
   ENDDO   

!  FInish

   RETURN
end subroutine Bessel_Cspher_plus

         
subroutine Eqv_radius_chebyshev_plus ( norder, deformation, radius, Le_radius )

   implicit none

!  input argument

   integer  , intent(in) ::  norder
   REAL(fpk), intent(in)  :: deformation

!  Output argument (equivalent radius)

   real(fpk), intent(inout) :: radius, Le_radius

!  Local variables

   integer    :: i, ngauss
   REAL(fpk)  :: X(60),W(60)
   real(fpk)  :: dn,xi,dx,dxn,ds,dsn,dcn,p13,cv
   real(fpk)  :: e,e2,en,s,v,a,a2,ens,ens2,t1,t2,rs,rv
   real(fpk)  :: L_e2,L_en,L_s,L_v,L_a,L_a2,L_ens,L_ens2,L_t1,L_t2,L_rs,L_rv

!  Constants

   E  = deformation
   DN = DBLE(norder)
   E2 = E*E
   EN = E*DN

   L_EN = EN
   L_E2 = 2.0d0 * E2

!  Quadrature

   ngauss = 60
!   CALL GAULEG_wrong (-1.0d0,1.0d0,X,W,ngauss)
   CALL GAULEG_right (NGAUSS,0,0,X,W)

!  Summation for surface-area and volume

   S = 0.0D0
   V = 0.0D0
   L_S = 0.0d0
   L_V = 0.0d0

   DO i = 1, ngauss

      XI = X(I)
      DX = DACOS(XI)
      DXN = DN*DX
      DS  = DSIN(DX)
      DSN = DSIN(DXN)
      DCN = DCOS(DXN)

      A   = 1.0D0 + E * DCN
      A2 = A * A
      ENS  = EN * DSN
      ENS2 = ENS * ENS

      L_A    = A - 1.0d0
      L_A2   = 2.0d0 * A * L_A
      L_ENS  = L_EN * DSN
      L_ENS2 = 2.0d0 * ENS * L_ENS

      T1 = DSQRT(A2+ENS2)
      T2 = DS*A + XI*ENS
      S = S + W(I) * A * T1
      V = V + W(I) * DS * A2 * T2

      L_T1 = 0.5d0 * ( L_A2 + L_ENS2 ) / T1
      L_T2 = DS * L_A + XI * L_ENS
      L_S = L_S + W(I) *      ( L_A  * T1 + A  * L_T1 )
      L_V = L_V + W(I) * DS * ( L_A2 * T2 + A2 * L_T2 )
   ENDDO

!  final answer for the radius

   P13 = 1.0D0/3.0D0
   CV  = 0.75d0  ** P13
   RS  = DSQRT(S*0.5D0)
   RV  = CV * V ** P13

   L_RS = L_S * 0.5d0 * RS / S
   L_RV = L_V * P13   * RV / V

   radius = RV/RS
   Le_radius = ( L_RV - radius * L_RS ) / RS

!  Finish

   RETURN
end subroutine Eqv_radius_chebyshev_plus


subroutine Eqv_radius_cylinder_plus ( Diamlength_ratio, radius, Le_radius )

   implicit none

!  input argument

   real(fpk), intent(in) :: Diamlength_ratio

!  Output argument (equivalent radius)

   real(fpk), intent(inout) :: radius, Le_radius

!  Local

   real(fpk) :: p13, c23, r1, r2, dinv, L_r1, L_r2
 
!  formula

   p13 = - 1.0D0/3.0D0
   c23 = ( - 2.0d0 * p13 ) ** p13
   r1 = c23 * Diamlength_ratio ** p13
   L_r1 = p13 * r1

   dinv = 1.0d0 / Diamlength_ratio
   r2 = 0.5d0 * dinv * (Diamlength_ratio+2.0D0)
   L_r2  = - dinv
   radius    = r1 / DSQRT(r2)
   Le_radius = radius * ( (L_r1/r1) - (0.5d0*L_r2/r2) )

!  Finish

   return   
end subroutine Eqv_radius_cylinder_plus


subroutine Eqv_radius_spheroids_plus ( aspect_ratio, radius, Le_radius )

   implicit none

!  input argument

   real(fpk), intent(in) :: aspect_ratio

!  Output argument (equivalent radius)

   real(fpk), intent(inout) :: radius, Le_radius

!  Local variables

   real(fpk) :: sinphi, phi, cphi, sinas, r1, r2, ratiosq, p23, radsq
   real(fpk) :: L_sinphi, L_phi, L_r1, L_r2, L_sinas

!  Prolate spheroids

   if ( aspect_ratio .lt. 1.0d0 ) then

      ratiosq = aspect_ratio * aspect_ratio
      p23  = 2.0d0 / 3.0d0
      r1   = 0.5d0 * aspect_ratio ** p23
      L_r1 =  p23 * r1

      sinphi = dsqrt(1.0d0 - aspect_ratio * aspect_ratio )
      phi    = dasin(sinphi)
      L_sinphi = - ratiosq / sinphi
      L_phi    = L_sinphi / aspect_ratio

      sinas = sinphi * aspect_ratio
      L_sinas = (L_sinphi + sinphi ) * aspect_ratio
      r2 = 1.0d0 + ( phi / sinas )
      L_r2 = ( L_phi - (r2 - 1.0d0) * L_sinas ) / sinas

      radius = 1.0d0 / dsqrt(r1*r2)
      radsq  = radius * radius
      Le_radius = -0.5d0 * radsq * radius * (L_r1*r2 + L_r2*r1)
      
   endif

!  Oblate spheroids

   if ( aspect_ratio .gt. 1.0d0 ) then

      ratiosq = aspect_ratio * aspect_ratio
      p23  = 2.0d0 / 3.0d0
      r1   = 0.25d0 * aspect_ratio ** p23
      L_r1 =  p23 * r1

      cphi   = dsqrt ( ratiosq - 1.0d0 )
      sinphi = cphi / aspect_ratio
      L_sinphi = 1.0d0 / cphi / aspect_ratio

      phi    = dlog((1.0d0+sinphi)/(1.0d0-sinphi))
      L_phi  =  L_sinphi * 2.0d0 * ratiosq

      sinas = ratiosq * sinphi
      L_sinas = ( L_sinphi + 2.0d0 * sinphi ) * ratiosq
      r2 = 2.0d0 + ( phi / sinas )
      L_r2 = ( L_phi - (r2 - 2.0d0) * L_sinas ) / sinas

      radius = 1.0d0 / dsqrt(r1*r2)
      radsq  = radius * radius
      Le_radius = -0.5d0 * radsq * radius * (L_r1*r2 + L_r2*r1)

   endif

      sinphi = dsqrt ( ratiosq - 1.0d0 ) / aspect_ratio
      phi    = dlog((1.0d0+sinphi)/(1.0d0-sinphi))
      r1 = 0.25d0 * aspect_ratio ** ( 2.0d0 / 3.0d0 )
      r2 = 2.0d0 + ( phi / ratiosq / sinphi )

!  Sphere

   if ( aspect_ratio .eq. 1.0d0 ) then
      radius    = 1.0d0
      Le_radius = 0.0d0
   endif

!  Finish

   return
end subroutine Eqv_radius_spheroids_plus

!  End module

end module tmat_functions_plus
