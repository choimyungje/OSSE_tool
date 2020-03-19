module SMCConv_Develop_m

public

contains

subroutine SMCConv_Develop_Mk2 &
   ( max_Angles, nangles, Cutoff, &
    cosines, weights, Fmat, Expcoeffs, Ncoeffs )

!  Stand-alone routine to develop coefficients from Scattering Matrix on Quadrature grid
!  Based on the Meerhoff Mie code (as found in RTSMie package), and adapted

!************************************************************************
!*  Calculate the expansion coefficients of the scattering matrix in    *
!*  generalized spherical functions by numerical integration over the   *
!*  scattering angle.                                   *
!************************************************************************

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input

   INTEGER          , INTENT (IN) :: max_Angles
   INTEGER          , INTENT (IN) :: nangles
 
   REAL    (KIND=dpk), INTENT (IN) :: Cutoff
   REAL    (KIND=dpk), INTENT (IN) :: cosines(max_Angles)
   REAL    (KIND=dpk), INTENT (IN) :: weights(max_Angles)
   REAL    (KIND=dpk), INTENT (IN) :: FMAT(max_Angles,6)

!  output. Initialized

   INTEGER           , INTENT (OUT) :: ncoeffs
   REAL    (KIND=dpk), INTENT (OUT) :: expcoeffs(0:max_Angles,6)

!  local variables

   REAL    (KIND=dpk) :: P00(max_Angles,2)
   REAL    (KIND=dpk) :: P02(max_Angles,2)
   REAL    (KIND=dpk) :: P22(max_Angles,2)
   REAL    (KIND=dpk) :: P2m2(max_Angles,2)
   REAL    (KIND=dpk) :: fmatw(6,max_Angles)

   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk
   real(dpk), parameter :: d_half  = 0.5_dpk, d_two  = 2.0_dpk
   real(dpk), parameter :: d_three = 3.0_dpk, d_four = 4.0_dpk

   INTEGER            :: i, j, l, lnew, lold, itmp
   INTEGER            :: index_11, index_12, index_22, index_33, index_34, index_44 
   REAL    (KIND=dpk) :: dl, dl2, qroot6, fac1, fac2, fac3, fl, &
                         sql4, sql41, twol1, tmp1, tmp2, denom, &
                         alfap, alfam, ps, mom1, Sig1

!  Initialization

  qroot6 = -0.25_dpk*SQRT(6.0_dpk)
  expcoeffs = d_zero
  ncoeffs = 0

!  Old indexing

!  ps = 1.0d0
!  index_11 = 1
!  index_12 = 2
!  index_22 = 3
!  index_33 = 4
!  index_34 = 5
!  index_44 = 6

!  New indexing consistent with Tmatrix output

  ps = -d_one
  index_11 = 1
  index_22 = 2
  index_33 = 3
  index_44 = 4
  index_12 = 5
  index_34 = 6

!  Multiply the scattering matrix F with the weights w for all angles  *
!  We do this here because otherwise it should be done for each l      *

  DO i = 1, 6
    DO j = 1, nangles
      fmatw(i,j) = weights(j)*FMAT(j,i)
    END DO
  END DO

!  Start loop over the coefficient index l                             *
!  first update generalized spherical functions, then calculate coefs. *
!  lold and lnew are pointer-like indices used in recurrence           *

  lnew = 1
  lold = 2
  mom1 = 1.0d0 ; l = -1
  Sig1 = 1.0d0
  DO while ( Sig1 .eq. 1.0d0 .and. mom1.gt.cutoff .and. l.lt. nangles )
    l = l + 1
!   do l = 0, ncoeffs

    IF (l == 0) THEN

      dl   = d_zero
      DO  i=1, nangles
        P00(i,lold) = d_one
        P00(i,lnew) = d_zero
        P02(i,lold) = d_zero
        P22(i,lold) = d_zero
        P2m2(i,lold)= d_zero
        P02(i,lnew) = d_zero
        P22(i,lnew) = d_zero
        P2m2(i,lnew)= d_zero
      END DO

    ELSE

      dl   = DBLE(l)
      dl2  = dl * dl
      fac1 = (d_two*dl-d_one)/dl
      fac2 = (dl-d_one)/dl
      DO  i=1, nangles
        P00(i,lold) = fac1*cosines(i)*P00(i,lnew) - fac2*P00(i,lold)
      END DO

    ENDIF

    IF (l == 2) THEN

      DO  i=1, nangles
        P02(i,lold) = qroot6*(d_one-cosines(i)*cosines(i))
        P22(i,lold) = 0.25_dpk*(d_one+cosines(i))*(d_one+cosines(i))
        P2m2(i,lold)= 0.25_dpk*(d_one-cosines(i))*(d_one-cosines(i))
        P02(i,lnew) = d_zero
        P22(i,lnew) = d_zero
        P2m2(i,lnew)= d_zero
      END DO
      sql41 = d_zero

    ELSE IF (l > 2) THEN

      sql4  = sql41
      sql41 = dsqrt(dl2-d_four)
      twol1 = d_two*dl - d_one
      tmp1  = twol1/sql41
      tmp2  = sql4/sql41
      denom = (dl-d_one)*(dl2-d_four)
      fac1  = twol1*(dl-d_one)*dl/denom
      fac2  = d_four*twol1/denom
      fac3  = dl*((dl-d_one)*(dl-d_one)-d_four)/denom
      DO i=1, nangles
        P02(i,lold) = tmp1*cosines(i)*P02(i,lnew)         - tmp2*P02(i,lold)
        P22(i,lold) = (fac1*cosines(i)-fac2)*P22(i,lnew)  - fac3*P22(i,lold)
        P2m2(i,lold)= (fac1*cosines(i)+fac2)*P2m2(i,lnew) - fac3*P2m2(i,lold)
      END DO

    END IF

    itmp = lnew
    lnew = lold
    lold = itmp
    alfap = d_zero
    alfam = d_zero

    fl = dl+d_half
    do i=1, nangles
      expcoeffs(L,index_11) = expcoeffs(L,index_11) + P00(i,lnew)*fmatw(1,i)
      alfap = alfap + P22(i,lnew)  * (fmatw(2,i)+fmatw(3,i))
      alfam = alfam + P2m2(i,lnew) * (fmatw(2,i)-fmatw(3,i))
      expcoeffs(L,index_44) = expcoeffs(L,index_44) + P00(i,lnew)*fmatw(4,i)
      expcoeffs(L,index_12) = expcoeffs(L,index_12) + P02(i,lnew)*fmatw(5,i)
      expcoeffs(L,index_34) = expcoeffs(L,index_34) + P02(i,lnew)*fmatw(6,i)
    END DO
    expcoeffs(L,index_11) =  fl*expcoeffs(L,index_11)
    expcoeffs(L,index_22) =  fl*d_half*(alfap+alfam)
    expcoeffs(L,index_33) =  fl*d_half*(alfap-alfam)
    expcoeffs(L,index_44) =  fl*expcoeffs(L,index_44)
    expcoeffs(L,index_12) =  fl*expcoeffs(L,index_12)
    expcoeffs(L,index_34) =  ps * fl*expcoeffs(L,index_34)
    if ( l.gt.0 ) mom1 = abs(expcoeffs(L,index_11))
    if ( l.gt.0 ) Sig1 = expcoeffs(L,index_11)/mom1

  END DO

!  Phase function normalization

  expcoeffs(0,index_11)        = d_one
  expcoeffs(L,index_11)        = d_zero
  ncoeffs = L - 1

!  finish

   return
end subroutine SMCConv_Develop_Mk2

!

subroutine SMCConv_Develop &
   ( max_Angles, ncoeffs, nangles, &
    cosines, weights, Fmat, Expcoeffs )

!  Stand-alone routine to develop coefficients from Scattering Matrix on Quadrature grid
!  Based on the Meerhoff Mie code (as found in RTSMie package), and adapted

!************************************************************************
!*  Calculate the expansion coefficients of the scattering matrix in    *
!*  generalized spherical functions by numerical integration over the   *
!*  scattering angle.                                   *
!************************************************************************

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input

   INTEGER          , INTENT (IN) :: max_Angles
   INTEGER          , INTENT (IN) :: ncoeffs, nangles
 
   REAL    (KIND=dpk), INTENT (IN) :: cosines(max_Angles)
   REAL    (KIND=dpk), INTENT (IN) :: weights(max_Angles)
   REAL    (KIND=dpk), INTENT (IN) :: FMAT(max_Angles,6)

!  output. Initialized

   REAL    (KIND=dpk), INTENT (OUT) :: expcoeffs(0:max_Angles,6)

!  local variables

   REAL    (KIND=dpk) :: P00(max_Angles,2)
   REAL    (KIND=dpk) :: P02(max_Angles,2)
   REAL    (KIND=dpk) :: P22(max_Angles,2)
   REAL    (KIND=dpk) :: P2m2(max_Angles,2)
   REAL    (KIND=dpk) :: fmatw(6,max_Angles)

   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk
   real(dpk), parameter :: d_half  = 0.5_dpk, d_two  = 2.0_dpk
   real(dpk), parameter :: d_three = 3.0_dpk, d_four = 4.0_dpk

   INTEGER            :: i, j, l, lnew, lold, itmp
   INTEGER            :: index_11, index_12, index_22, index_33, index_34, index_44 
   REAL    (KIND=dpk) :: dl, dl2, qroot6, fac1, fac2, fac3, fl, &
                         sql4, sql41, twol1, tmp1, tmp2, denom, &
                         alfap, alfam, ps

!  Initialization

  qroot6 = -0.25_dpk*SQRT(6.0_dpk)
  expcoeffs = d_zero

!  Old indexing

!  ps = 1.0d0
!  index_11 = 1
!  index_12 = 2
!  index_22 = 3
!  index_33 = 4
!  index_34 = 5
!  index_44 = 6

!  New indexing consistent with Tmatrix output

  ps = -d_one
  index_11 = 1
  index_22 = 2
  index_33 = 3
  index_44 = 4
  index_12 = 5
  index_34 = 6

!  Multiply the scattering matrix F with the weights w for all angles  *
!  We do this here because otherwise it should be done for each l      *

  DO i = 1, 6
    DO j = 1, nangles
      fmatw(i,j) = weights(j)*FMAT(j,i)
    END DO
  END DO

!  Start loop over the coefficient index l                             *
!  first update generalized spherical functions, then calculate coefs. *
!  lold and lnew are pointer-like indices used in recurrence           *

  lnew = 1
  lold = 2

  DO l = 0, ncoeffs

    IF (l == 0) THEN

      dl   = d_zero
      DO  i=1, nangles
        P00(i,lold) = d_one
        P00(i,lnew) = d_zero
        P02(i,lold) = d_zero
        P22(i,lold) = d_zero
        P2m2(i,lold)= d_zero
        P02(i,lnew) = d_zero
        P22(i,lnew) = d_zero
        P2m2(i,lnew)= d_zero
      END DO

    ELSE

      dl   = DBLE(l)
      dl2  = dl * dl
      fac1 = (d_two*dl-d_one)/dl
      fac2 = (dl-d_one)/dl
      DO  i=1, nangles
        P00(i,lold) = fac1*cosines(i)*P00(i,lnew) - fac2*P00(i,lold)
      END DO

    ENDIF

    IF (l == 2) THEN

      DO  i=1, nangles
        P02(i,lold) = qroot6*(d_one-cosines(i)*cosines(i))
        P22(i,lold) = 0.25_dpk*(d_one+cosines(i))*(d_one+cosines(i))
        P2m2(i,lold)= 0.25_dpk*(d_one-cosines(i))*(d_one-cosines(i))
        P02(i,lnew) = d_zero
        P22(i,lnew) = d_zero
        P2m2(i,lnew)= d_zero
      END DO
      sql41 = d_zero

    ELSE IF (l > 2) THEN

      sql4  = sql41
      sql41 = dsqrt(dl2-d_four)
      twol1 = d_two*dl - d_one
      tmp1  = twol1/sql41
      tmp2  = sql4/sql41
      denom = (dl-d_one)*(dl2-d_four)
      fac1  = twol1*(dl-d_one)*dl/denom
      fac2  = d_four*twol1/denom
      fac3  = dl*((dl-d_one)*(dl-d_one)-d_four)/denom
      DO i=1, nangles
        P02(i,lold) = tmp1*cosines(i)*P02(i,lnew)         - tmp2*P02(i,lold)
        P22(i,lold) = (fac1*cosines(i)-fac2)*P22(i,lnew)  - fac3*P22(i,lold)
        P2m2(i,lold)= (fac1*cosines(i)+fac2)*P2m2(i,lnew) - fac3*P2m2(i,lold)
      END DO

    END IF

    itmp = lnew
    lnew = lold
    lold = itmp
    alfap = d_zero
    alfam = d_zero

    fl = dl+d_half
    do i=1, nangles
      expcoeffs(L,index_11) = expcoeffs(L,index_11) + P00(i,lnew)*fmatw(1,i)
      alfap = alfap + P22(i,lnew)  * (fmatw(2,i)+fmatw(3,i))
      alfam = alfam + P2m2(i,lnew) * (fmatw(2,i)-fmatw(3,i))
      expcoeffs(L,index_44) = expcoeffs(L,index_44) + P00(i,lnew)*fmatw(4,i)
      expcoeffs(L,index_12) = expcoeffs(L,index_12) + P02(i,lnew)*fmatw(5,i)
      expcoeffs(L,index_34) = expcoeffs(L,index_34) + P02(i,lnew)*fmatw(6,i)
    END DO
    expcoeffs(L,index_11) =  fl*expcoeffs(L,index_11)
    expcoeffs(L,index_22) =  fl*d_half*(alfap+alfam)
    expcoeffs(L,index_33) =  fl*d_half*(alfap-alfam)
    expcoeffs(L,index_44) =  fl*expcoeffs(L,index_44)
    expcoeffs(L,index_12) =  fl*expcoeffs(L,index_12)
    expcoeffs(L,index_34) =  ps * fl*expcoeffs(L,index_34)
  END DO

!  Phase function normalization

  expcoeffs(0,index_11)        = d_one

!  finish

   return
end subroutine SMCConv_Develop

!  End module

End Module SMCConv_Develop_m
