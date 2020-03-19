 module tmat_distributions

   USE tmat_parameters,    ONLY :  fpk => tmat_fpkind, &
                                   d_zero, d_one, d_two

!  New distribution software

!     DISTRB_bulk_plus
!     DISTRB_bulk
!     DISTRB_new_plus
!     DISTRB_new

!  Old sitributions

!     distr, power, zeroin

!  Private routines

!     sizedist_nod
!     gammafunction
!     rminmax

private gammafunction, rminmax

public  DISTRB_bulk_plus, DISTRB_bulk, &
        DISTRB_new_plus,  DISTRB_new,  &
        DISTRB, POWER, ZEROIN

contains

subroutine DISTRB_bulk_plus                          &
      ( MAXK, XG1, NK, pars, nr, nr_derivs,          & ! Inputs
        WG1, WG1_derivs, tmat_dist, LPSD_tmat_dist )   ! Outputs

!  modules

   USE tmat_parameters,    ONLY : d_zero, d_one, d_two

   IMPLICIT NONE

! subroutine input arguments

   INTEGER          , INTENT (IN)     ::  maxk, nk
   REAL     (KIND=fpk), INTENT (IN)     ::  XG1(maxk)
   REAL     (KIND=fpk), INTENT (IN)     ::  pars(3)
   REAL     (KIND=fpk), INTENT (IN)     ::  nr(maxk)
   REAL     (KIND=fpk), INTENT (IN)     ::  nr_derivs(maxk,3)

!  subroutine modified output

   REAL     (KIND=fpk), INTENT (INOUT)  ::  WG1(maxk)
   REAL     (KIND=fpk), INTENT (INOUT)  ::  WG1_derivs(maxk,3)

!  Subroutine output (Distribution parameters)
!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(KIND=fpk), intent(out) :: Tmat_dist (5)
   real(KIND=fpk), intent(out) :: LPSD_Tmat_dist (5,3)

!  Local variables

   INTEGER        :: i, kd
   REAL (KIND=fpk)  :: d_pi, quad, quadr2, quadr3, quadr4, factor
   REAL (KIND=fpk)  :: quad_d, quadr2_d, quadr3_d, quadr4_d, awt
   REAL (KIND=fpk)  :: ndens, ndens_d(3), gxsec, reff, volume, veff
   REAL (KIND=fpk)  :: gxsec_d(3), reff_d(3), volume_d(3), veff_d(3)
   REAL (KIND=fpk)  :: J2_d(3), J3_d(3), J4_d(3), J2, J3, J4

!  Initialise

   Tmat_dist      = d_zero
   LPSD_Tmat_dist = d_zero

!  Normalized density and derivatives

   ndens = d_zero ;  ndens_d = d_zero
   DO i = 1, nk
      ndens  = ndens + wg1(i) * nr(i)
      DO kd = 1, 3
         quad_d   = nr_derivs(i,kd) * wg1(i)
         ndens_d(kd) = ndens_d(kd) + quad_d
      ENDDO
   ENDDO

!************************************
!  Do not normalize WG1 until the end
!************************************

!  Find moments of PSD
!    Number density, geometric cross-section, 3rd and 4th powers

   J2   = d_zero ; J3   = d_zero ; J4   = d_zero
   J2_d = d_zero ; J3_d = d_zero ; J4_d = d_zero
   DO i = 1, nk
      quad   = nr(i)  * wg1(i)
      quadr2 = quad   * xg1(i) * xg1(i)
      quadr3 = quadr2 * xg1(i)
      quadr4 = quadr3 * xg1(i)
      J2 = J2 + quadr2
      J3 = J3 + quadr3
      J4 = J4 + quadr4
      DO kd = 1, 3
         quad_d   = nr_derivs(i,kd) * wg1(i)
         quadr2_d = quad_d   * xg1(i) * xg1(i)
         quadr3_d = quadr2_d * xg1(i)
         quadr4_d = quadr3_d * xg1(i)
         J2_d(kd) = J2_d(kd) + quadr2_d
         J3_d(kd) = J3_d(kd) + quadr3_d
         J4_d(kd) = J4_d(kd) + quadr4_d
      ENDDO
   ENDDO

!  geometric cross-section and derivatives

   d_pi    = dacos(-d_one)
   factor  = d_pi
   gxsec   = J2 / ndens
   DO kd = 1, 3
      gxsec_d(kd) = factor * ( J2_d(kd) - gxsec * ndens_d(kd) ) / ndens
   END DO
   gxsec   = factor * gxsec
 
!  VOLUME

   factor  = d_pi * 4.0d0 / 3.0d0
   volume = J3 / ndens
   DO kd = 1, 3
      volume_d(kd) = factor * ( J3_d(kd) - volume * ndens_d(kd) ) / ndens
   END DO
   volume = factor * volume

!  REFF

   reff = J3 / J2
   DO kd = 1, 3
      reff_d(kd) = ( J3_d(kd) - reff * J2_d(kd) ) / J2
   END DO

!  VEFF

   factor = 1.0d0 / J3 / J3
   veff  = ( J4 * J2 * factor ) - d_one
   DO kd = 1, 3
      veff_d(kd) = factor * ( J4_d(kd) * J2 + J2_d(kd) * J4 )  &
                    - d_two * J3_d(kd) * ( veff + d_one ) / J3
   END DO

!  Assign output

  tmat_dist(1) = ndens
  tmat_dist(2) = gxsec
  tmat_dist(3) = volume
  tmat_dist(4) = reff
  tmat_dist(5) = veff

   DO kd = 1, 3
      LPSD_tmat_dist(1,kd) = ndens_d(kd)   
      LPSD_tmat_dist(2,kd) = gxsec_d(kd)   
      LPSD_tmat_dist(3,kd) = volume_d(kd)   
      LPSD_tmat_dist(4,kd) = reff_d(kd)   
      LPSD_tmat_dist(5,kd) = veff_d(kd)
   END DO

!  Normalize derivatives

   DO kd = 1, 3
      do i = 1, 5
         LPSD_tmat_dist(i,kd) = pars(kd) * LPSD_tmat_dist(i,kd)
      enddo
   enddo

!  Scaling and Density normalization of weights

   DO i = 1, nk
      awt    = wg1(i)
      WG1(i) = awt * nr(i) / ndens
      do kd = 1, 3
         quad_d   = nr_derivs(i,kd) * awt
         WG1_derivs(i,kd) = ( quad_d - WG1(i) * ndens_d(kd) ) / ndens
         WG1_derivs(i,kd) = pars(kd) * WG1_derivs(i,kd)
      enddo
   ENDDO

!  Finish

  RETURN
END SUBROUTINE DISTRB_bulk_plus

subroutine DISTRB_bulk      &
      ( MAXK, XG1, NK, nr,  & ! Inputs
        WG1, tmat_dist )      ! Outputs

!  !  modules

   USE tmat_parameters,    ONLY : d_zero, d_one

   IMPLICIT NONE

! subroutine input arguments

   INTEGER          , INTENT (IN)     ::  maxk, nk
   REAL     (KIND=fpk), INTENT (IN)     ::  XG1(maxk)
   REAL     (KIND=fpk), INTENT (IN)     ::  nr(maxk)

!  subroutine modified output

   REAL     (KIND=fpk), INTENT (INOUT)  ::  WG1(maxk)

!  Subroutine output (Distribution parameters)
!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(KIND=fpk), intent(out) :: Tmat_dist (5)

!  Local variables

   INTEGER        :: i
   REAL (KIND=fpk)  :: d_pi, factor, J2, J3, J4
   REAL (KIND=fpk)  :: quad, quadr2, quadr3, quadr4
   REAL (KIND=fpk)  :: ndens, gxsec, reff, volume, veff

!  Initialise

   Tmat_dist      = d_zero

!************************************
!  Do not normalize WG1 until the end
!************************************

!  Find density and moments of PSD

   ndens = d_zero ; J2   = d_zero ; J3   = d_zero ; J4   = d_zero
   DO i = 1, nk
      quad   = nr(i)  * wg1(i)
      quadr2 = quad   * xg1(i) * xg1(i)
      quadr3 = quadr2 * xg1(i)
      quadr4 = quadr3 * xg1(i)
      ndens  = ndens + quad
      J2     = J2 + quadr2
      J3     = J3 + quadr3
      J4     = J4 + quadr4
   ENDDO

!  geometric cross-section and derivatives

   d_pi    = dacos(-d_one)
   factor  = d_pi
   gxsec   = J2 / ndens
   gxsec   = factor * gxsec
 
!  VOLUME

   factor  = d_pi * 4.0d0 / 3.0d0
   volume = J3 / ndens
   volume = factor * volume

!  REFF

   reff = J3 / J2

!  VEFF

   factor = 1.0d0 / J3 / J3
   veff  = J4 * J2 * factor - d_one

!  Assign Distribution output

   tmat_dist(1) = ndens
   tmat_dist(2) = gxsec
   tmat_dist(3) = volume
   tmat_dist(4) = reff
   tmat_dist(5) = veff

!  Compound and Normalization by Ndens

   DO i = 1, nk
      wg1(i) = wg1(i) * nr(i) / ndens
   ENDDO

!  Finish

  RETURN
END SUBROUTINE DISTRB_bulk


SUBROUTINE DISTRB_new_plus                          &
    ( max_distpoints, idis, par, radius, numradius, & ! Inputs
      nwithr, nwithr_d, message, faild )              ! Outputs

!************************************************************************
!*  Calculate the size distribution n(r) for the numr radius values     *
!*  contained in array r and return the results through the array nwithr*
!*  The size distributions are normalized such that the integral over   *
!*  all r is equal to one.                                              *
!************************************************************************

!  modules

   USE tmat_parameters,    ONLY : d_zero, d_half, d_one, d_two, d_three

   IMPLICIT NONE

! subroutine input arguments

   INTEGER          , INTENT (IN)  :: max_distpoints
   REAL     (KIND=fpk), INTENT (IN)  :: par(3)
   INTEGER          , INTENT (IN)  :: idis, numradius
   REAL     (KIND=fpk), INTENT (IN)  :: radius(max_distpoints)

!  Subroutine output

   REAL    (KIND=fpk), INTENT (OUT) :: nwithr  (max_distpoints)
   REAL    (KIND=fpk), INTENT (OUT) :: nwithr_d(max_distpoints,3)
   CHARACTER*(*)   , INTENT (OUT) :: message
   LOGICAL         , INTENT (OUT) :: faild

! local variables

  INTEGER        :: i
  LOGICAL        :: deriv(3)
  REAL  (KIND=fpk) :: pi,r,logr,root2p
  REAL  (KIND=fpk) :: alpha,alpha1,b,logb,arg1,arg2,arg,argsq,r3
  REAL  (KIND=fpk) :: b1,b2,b11,b13,b22,b23,logb1,logb2,rc
  REAL  (KIND=fpk) :: logrg,logsi,logsi_inv,fac_d1,gamma,gamma1,rg
  REAL  (KIND=fpk) :: rmin,rmax,fac1,fac2,aperg
  REAL  (KIND=fpk) :: alpha2, fac_d2a
  REAL  (KIND=fpk) :: n1, n2, n1_d1, n1_d3, n2_d2, n2_d3

!  redundant variables
!  REAL  (KIND=dp) :: sigfac, logC1_d2

  REAL  (KIND=fpk) :: C, logC, logC_d1, logC_d2, logC_d3
  REAL  (KIND=fpk) :: logC1, logC2, logC1_d1, logC1_d3, logC2_d2, logC2_d3

  REAL    (KIND=fpk)  :: gammln, dgammln
  CHARACTER*70      :: message_gamma
  LOGICAL           :: fail
  character*1       :: cdis

!  Initialize

   faild = .FALSE.
   message = ' '
   nwithr(1:numradius)       = d_zero
   nwithr_d(1:numradius,1:3) = d_zero

!  Set local derivs flags (this could be an input)

   deriv = .false.
   deriv(1) = .true.
   if ( idis .ne. 6 ) then
      deriv(2) = .true.
      if ( idis==3 .or. idis==7 .or. idis==8 ) deriv(3) = .true.
   endif

!  Check

   if (idis == 0 ) RETURN
   IF ( IDIS > 8 ) THEN
      faild = .TRUE.
      message = 'illegal index in sizedis'
      RETURN
   END IF

!  setup

   pi     = dacos(-1.d0)
   root2p = dsqrt(pi+pi)

!  IDIS = 1 : TWO-PARAMETER GAMMA with alpha and b given

   IF ( idis == 1 ) THEN

      deriv(1) = .true.
      deriv(2) = .true.
      alpha  = par(1)
      b      = par(2)
      alpha1 = alpha + d_one
      logb   = LOG(b)
      CALL gammafunction ( alpha1, deriv(1), gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240
      logC  = alpha1*logb - gammln

      IF ( deriv(1) .and. deriv(2)  ) then
        logC_d1 = logb - dgammln
        logC_d2 = alpha1/b
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          nwithr(i) = EXP ( arg1 - b*r )
          nwithr_d(i,2) = ( logC_d2 - r )    * nwithr(i)
          nwithr_d(i,1) = ( logC_d1 + logr ) * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r  = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          nwithr(i) = EXP ( arg1 - b*r )
        END DO
      END IF

!  IDIS = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given

   ELSE IF ( idis == 2 ) THEN

      alpha  = d_one/par(2) - d_three
      b      = d_one/(par(1)*par(2))
      alpha1 = alpha + d_one
      logb   = LOG(b)
      CALL gammafunction ( alpha1, deriv(2), gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = alpha1*logb - gammln

      IF ( deriv(1) .and. deriv(2) ) then
        b1      = b / par(1)
        b2      = b / par(2) 
        logC_d1 = - alpha1 / par(1)
        logC_d2 = ( dgammln - logb - alpha1*par(2) ) / par(2) / par(2)
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          nwithr(i) = EXP ( arg1 - b*r )
          nwithr_d(i,1) = ( logC_d1 + b1*r  )  * nwithr(i)
          nwithr_d(i,2) = ( logC_d2 - logr/par(2)/par(2) + b2*r )  * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          nwithr(i) = EXP ( arg1 - b*r )
        END DO
      END IF

!  IDIS = 3 : BIMODAL GAMMA with equal mode weights

   ELSE IF ( idis == 3 ) THEN

      alpha  = d_one/par(3) - d_three
      b1     = d_one/(par(1)*par(3))
      b2     = d_one/(par(2)*par(3))
      alpha1 = alpha + d_one
      CALL gammafunction ( alpha1, deriv(3), gammln, dgammln, fail, message_gamma )
      logb1 = LOG(b1)
      logb2 = LOG(b2)
      logC1 = alpha1*logb1 - gammln
      logC2 = alpha1*logb2 - gammln

      IF ( deriv(1) .and. deriv(2) .and. deriv(3) ) then
        b11      = b1 / par(1)
        b13      = b1 / par(3)
        b22      = b2 / par(2)
        b23      = b2 / par(3)
        logC1_d1 = - alpha1 / par(1)
        logC1_d3 = ( dgammln - logb1 - alpha1*par(3) ) / par(3) / par(3)
        logC2_d2 = - alpha1 / par(2)
        logC2_d3 = ( dgammln - logb2 - alpha1*par(3) ) / par(3) / par(3)
        DO i = 1, numradius
          r  = radius(i)
          logr = LOG(r)
          arg1 = logC1 + alpha*logr
          arg2 = logC2 + alpha*logr
          n1   = EXP(arg1 - b1*r)
          n2   = EXP(arg2 - b2*r)
          nwithr(i) = d_half * ( n1 + n2 )
          n1_d1 = ( logC1_d1 + b11*r  )                 * n1
          n1_d3 = ( logC1_d3 - logr/par(3)/par(3) + b13*r ) * n1
          n2_d2 = ( logC2_d2 + b22*r  )                 * n2
          n2_d3 = ( logC2_d3 - logr/par(3)/par(3) + b23*r ) * n2
          nwithr_d(i,1) = d_half * n1_d1
          nwithr_d(i,2) = d_half * n2_d2
          nwithr_d(i,3) = d_half * ( n1_d3 + n2_d3 ) 
        END DO
      ELSE
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC1 + alpha*logr
          arg2 = logC2 + alpha*logr
          n1   = EXP(arg1 - b1*r)
          n2   = EXP(arg2 - b2*r)
          nwithr(i) = d_half * ( n1 + n2 )
        END DO
      END IF
  
!  4  LOG-NORMAL with rg and sigma given

   ELSE IF ( idis == 4 ) THEN

      logrg = dlog(par(1))
      logsi = dabs(dlog(par(2)))
      logsi_inv = d_one / logsi
      C      = logsi_inv / root2p
      IF ( deriv(1) .and. deriv(2) ) then
        logC_d2 = - logsi_inv / par(2)
        fac_d1  =   logsi_inv / par(1)
        DO i = 1, numradius
          r     = radius(i)
          logr  = LOG(r)
          arg   = ( logr - logrg ) / logsi
          argsq = arg * arg
          nwithr(i) = C * dexp( - d_half * argsq ) / r
          nwithr_d(i,1) = arg * fac_d1 * nwithr(i)
          nwithr_d(i,2) = logC_d2 * ( d_one - argsq ) * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r     = radius(i)
          logr  = LOG(r)
          arg   = ( logr - logrg ) / logsi
          argsq = arg * arg
          nwithr(i) = C * dexp( - d_half * argsq ) / r
        END DO
      END IF

!  5 LOG-NORMAL with reff and veff given                               *

   ELSE IF ( idis == 5 ) THEN

      alpha1 = d_one + par(2)
      alpha2 = dlog(alpha1)
      rg     = par(1)/(d_one+par(2))**2.5d0
      logrg  = dlog(rg)
      logsi  = dsqrt(alpha2)
      logsi_inv = d_one / logsi
      C         = logsi_inv / root2p
      IF ( deriv(1) .and. deriv(2) ) then
        logC_d2 = - d_half / alpha2 / alpha1
        fac_d1  =   logsi_inv / par(1)
        fac_d2a  =  - 2.5d0 * logsi_inv / alpha1
        DO i = 1, numradius
          r     = radius(i)
          logr  = LOG(r)
          arg   = ( logr - logrg ) / logsi
          argsq = arg * arg
          nwithr(i) = C * dexp( - d_half * argsq ) / r
          nwithr_d(i,1) = arg * fac_d1 * nwithr(i)
          nwithr_d(i,2) = ( arg * fac_d2a + logC_d2*(d_one-argsq) ) * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r     = radius(i)
          logr  = LOG(r)
          arg   = ( logr - logrg ) / logsi
          argsq = arg * arg
          nwithr(i) = C * dexp( - d_half * argsq ) / r
        END DO
      END IF

!  6 POWER LAW                               *

   ELSE IF ( idis == 6 ) THEN

      alpha = par(1)
      rmin  = par(2)
      rmax  = par(3)
      alpha1 = alpha - d_one
      fac1 = (rmax/rmin)**alpha1
      fac2 = d_one / ( fac1 - d_one )
      C = alpha1 * rmax**alpha1 * fac2
      IF ( deriv(1) .and. deriv(2) .and. deriv(3) ) then
        logC_d1 = (d_one/alpha1) + LOG(par(3)) - fac1 * fac2 * LOG(par(3)/par(2))
        DO i = 1, numradius
          r     = radius(i)
          if ( (r < rmax) .and. (r > rmin) ) then
            nwithr(i)    = C*r**(-alpha)
            nwithr_d(i,1) = ( logC_d1 - log(r) ) * nwithr(i)
          else
            nwithr(i)    = d_zero
            nwithr_d(i,1) = d_zero
          endif
        END DO
      ELSE
        DO i = 1, numradius
          r     = radius(i)
          if ( (r < rmax) .and. (r > rmin) ) then
            nwithr(i)    = C*r**(-alpha)
          else
            nwithr(i)    = d_zero
          endif
        END DO
      END IF

!  7 MODIFIED GAMMA with alpha, rc and gamma given

   ELSE IF ( idis == 7 ) THEN

      alpha = par(1)
      rc    = par(2)
      gamma = par(3)
      b     = alpha / (gamma*rc**gamma)
      logb  = LOG(b)
      alpha1 = alpha + d_one
      gamma1 = d_one / gamma
      aperg = alpha1/gamma
      CALL gammafunction ( aperg, deriv(1), gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      IF ( deriv(1) .and. deriv(2) .and. deriv(3) ) then
        logC_d1 = ( logb - dgammln ) * gamma1 + aperg/par(1)
        logC_d2 = - aperg * gamma / par(2)
        logC_d3 = gamma1 - aperg * ( logb - dgammln ) * gamma1 - aperg * (gamma1 + LOG(par(2)) ) 
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          r3   = b * r ** gamma
          nwithr(i) = EXP ( arg1 - r3 )
          nwithr_d(i,1) = ( logC_d1 + logr - r3/par(1) ) * nwithr(i)
          nwithr_d(i,2) = ( logC_d2 + r3*gamma/par(2) )  * nwithr(i)
          nwithr_d(i,3) = ( logC_d3 + r3*(gamma1+LOG(par(2))-logr) ) * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          r3   = b*r ** gamma
          nwithr(i) = EXP ( arg1 - r3 )
        END DO
      END IF

!  8 MODIFIED GAMMA with alpha, b and gamma given

   ELSE IF ( idis == 8 ) THEN

      alpha = par(1)
      b     = par(2)
      gamma = par(3)
      alpha1 = alpha + d_one
      gamma1 = d_one / gamma
      logb = LOG(b)
      aperg = alpha1/gamma
      CALL gammafunction ( aperg, deriv(1), gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      IF ( deriv(1) .and. deriv(2) .and. deriv(3) ) then
        b1      = b / par(1)
        b2      = b / par(2)
        logC_d1 = ( logb - dgammln ) * gamma1
        logC_d2 = aperg / b
        logC_d3 = gamma1 - aperg * logC_d1
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          r3   = r ** gamma
          nwithr(i) = EXP ( arg1 - b*r3 )
          nwithr_d(i,1) = ( logC_d1 + logr )       * nwithr(i)
          nwithr_d(i,2) = ( logC_d2 - r3 )         * nwithr(i)
          nwithr_d(i,3) = ( logC_d3 - b*logr*r3 )  * nwithr(i)
        END DO
      ELSE
        DO i = 1, numradius
          r    = radius(i)
          logr = LOG(r)
          arg1 = logC + alpha*logr
          r3   = r ** gamma
          nwithr(i) = EXP ( arg1 - b*r3 )
        END DO
      END IF

   END IF

!  normal return

   RETURN

!  special return

240  CONTINUE
   faild = .TRUE.
   write(cdis,'(I1)')idis
   message = message_gamma(1:LEN(message_gamma))//', distribution : '//cdis
   RETURN

END SUBROUTINE DISTRB_new_plus


SUBROUTINE DISTRB_new                               &
    ( max_distpoints, idis, par, radius, numradius, &
      nwithr, message, faild )

!************************************************************************
!*  Calculate the size distribution n(r) for the numr radius values     *
!*  contained in array r and return the results through the array nwithr*
!*  The size distributions are normalized such that the integral over   *
!*  all r is equal to one.                                              *
!************************************************************************

!  modules

   USE tmat_parameters,    ONLY : d_zero, d_half, d_one, d_two, d_three

   IMPLICIT NONE

! subroutine input arguments

   INTEGER          , INTENT (IN)  :: max_distpoints
   REAL     (KIND=fpk), INTENT (IN)  :: par(3)
   INTEGER          , INTENT (IN)  :: idis, numradius
   REAL     (KIND=fpk), INTENT (IN)  :: radius(max_distpoints)

!  Subroutine output

   REAL    (KIND=fpk), INTENT (OUT) :: nwithr(max_distpoints)
   CHARACTER*(*)   , INTENT (OUT) :: message
   LOGICAL         , INTENT (OUT) :: faild

! local variables

   INTEGER        :: i
   REAL  (KIND=fpk) :: pi,r,logr,root2p
   REAL  (KIND=fpk) :: alpha,alpha1,b,logb,arg1,arg2,arg,argsq,r3
   REAL  (KIND=fpk) :: b1,b2,logb1,logb2,rc
   REAL  (KIND=fpk) :: logrg,logsi,logsi_inv,gamma,gamma1,rg
   REAL  (KIND=fpk) :: rmin,rmax,fac1,fac2,aperg, alpha2
   REAL  (KIND=fpk) :: n1, n2, C, logC, logC1, logC2

   REAL    (KIND=fpk)  :: gammln, dummy
   CHARACTER*70      :: message_gamma
   LOGICAL           :: fail
   character*1       :: cdis

!  Initialise

   faild = .FALSE.
   message = ' '
   nwithr(1:numradius) = d_zero

!  Check index

   if (idis == 0 ) RETURN
   IF ( IDIS > 8 ) THEN
      faild = .TRUE.
      message = 'illegal index in sizedis'
      RETURN
   END IF

!  setup

   pi     = dacos(-1.d0)
   root2p = dsqrt(pi+pi)

!  IDIS = 1 : TWO-PARAMETER GAMMA with alpha and b given

   IF ( idis == 1 ) THEN

      alpha  = par(1)
      b      = par(2)
      alpha1 = alpha + d_one
      logb   = LOG(b)
      CALL gammafunction ( alpha1, .false., gammln, dummy, fail, message_gamma )
      IF ( fail ) go to 240
      logC  = alpha1*logb - gammln

      DO i = 1, numradius
         r  = radius(i)
         logr = LOG(r)
         arg1 = logC + alpha*logr
         nwithr(i) = EXP ( arg1 - b*r )
      END DO

!  IDIS = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given

   ELSE IF ( idis == 2 ) THEN

      alpha  = d_one/par(2) - d_three
      b      = d_one/(par(1)*par(2))
      alpha1 = alpha + d_one
      logb   = LOG(b)
      CALL gammafunction ( alpha1, .false., gammln, dummy, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = alpha1*logb - gammln

      DO i = 1, numradius
         r    = radius(i)
         logr = LOG(r)
         arg1 = logC + alpha*logr
         nwithr(i) = EXP ( arg1 - b*r )
      END DO

!  IDIS = 3 : BIMODAL GAMMA with equal mode weights

   ELSE IF ( idis == 3 ) THEN

      alpha  = d_one/par(3) - d_three
      b1     = d_one/(par(1)*par(3))
      b2     = d_one/(par(2)*par(3))
      alpha1 = alpha + d_one
      CALL gammafunction ( alpha1, .false., gammln, dummy, fail, message_gamma )
      logb1 = LOG(b1)
      logb2 = LOG(b2)
      logC1 = alpha1*logb1 - gammln
      logC2 = alpha1*logb2 - gammln

      DO i = 1, numradius
         r    = radius(i)
         logr = LOG(r)
         arg1 = logC1 + alpha*logr
         arg2 = logC2 + alpha*logr
         n1   = EXP(arg1 - b1*r)
         n2   = EXP(arg2 - b2*r)
         nwithr(i) = d_half * ( n1 + n2 )
      END DO
  
!  4  LOG-NORMAL with rg and sigma given

   ELSE IF ( idis == 4 ) THEN

      logrg = dlog(par(1))
      logsi = dabs(dlog(par(2)))
      logsi_inv = d_one / logsi
      C      = logsi_inv / root2p
      DO i = 1, numradius
         r     = radius(i)
         logr  = LOG(r)
         arg   = ( logr - logrg ) / logsi
         argsq = arg * arg
         nwithr(i) = C * dexp( - d_half * argsq ) / r
      END DO

!  5 LOG-NORMAL with reff and veff given                               *

   ELSE IF ( idis == 5 ) THEN

      alpha1 = d_one + par(2)
      alpha2 = dlog(alpha1)
      rg     = par(1)/(d_one+par(2))**2.5d0
      logrg  = dlog(rg)
      logsi  = dsqrt(alpha2)
      logsi_inv = d_one / logsi
      C         = logsi_inv / root2p
      DO i = 1, numradius
         r     = radius(i)
         logr  = LOG(r)
         arg   = ( logr - logrg ) / logsi
         argsq = arg * arg
         nwithr(i) = C * dexp( - d_half * argsq ) / r
      END DO

!  6 POWER LAW                               *

   ELSE IF ( idis == 6 ) THEN

      alpha = par(1)
      rmin  = par(2)
      rmax  = par(3)
      alpha1 = alpha - d_one
      fac1 = (rmax/rmin)**alpha1
      fac2 = d_one / ( fac1 - d_one )
      C = alpha1 * rmax**alpha1 * fac2
      DO i = 1, numradius
         r     = radius(i)
         if ( (r < rmax) .and. (r > rmin) ) then
            nwithr(i)    = C*r**(-alpha)
         else
            nwithr(i)    = d_zero
         endif
      END DO

!  7 MODIFIED GAMMA with alpha, rc and gamma given

   ELSE IF ( idis == 7 ) THEN

      alpha = par(1)
      rc    = par(2)
      gamma = par(3)
      b     = alpha / (gamma*rc**gamma)
      logb  = LOG(b)
      alpha1 = alpha + d_one
      gamma1 = d_one / gamma
      aperg = alpha1/gamma
      CALL gammafunction ( aperg, .false., gammln, dummy, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      DO i = 1, numradius
         r    = radius(i)
         logr = LOG(r)
         arg1 = logC + alpha*logr
         r3   = b*r ** gamma
         nwithr(i) = EXP ( arg1 - r3 )
      END DO

!  8 MODIFIED GAMMA with alpha, b and gamma given

   ELSE IF ( idis == 8 ) THEN

      alpha = par(1)
      b     = par(2)
      gamma = par(3)
      alpha1 = alpha + d_one
      gamma1 = d_one / gamma
      logb = LOG(b)
      aperg = alpha1/gamma
      CALL gammafunction ( aperg, .false., gammln, dummy, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      DO i = 1, numradius
         r    = radius(i)
         logr = LOG(r)
         arg1 = logC + alpha*logr
         r3   = r ** gamma
         nwithr(i) = EXP ( arg1 - b*r3 )
      END DO

   END IF

!  normal return

   RETURN

!  special return

240  CONTINUE
   faild = .TRUE.
   write(cdis,'(I1)')idis
   message = message_gamma(1:LEN(message_gamma))//', distribution : '//cdis
   RETURN

END SUBROUTINE DISTRB_new


subroutine DISTRB                                  &
    ( NNK, YY, NDISTR, R1, DO_OUTPUT, AA, BB, GAM, & ! Inputs
      WY, NDENS, GXSEC, VOLUME, REFF, VEFF)          ! Outputs
   
   implicit none

!  input arguments
!  ---------------

!  Number of points and radii

   integer     , intent(in)  ::  NNK
   real(KIND=fpk), intent(in)  ::  YY(NNK)

!  Distribution index and parameters

   integer     , intent(in)  :: NDISTR
   real(KIND=fpk), intent(in)  :: AA,BB,GAM

!  Minimum radius

   real(KIND=fpk), intent(in)  :: R1

!  output flag

   logical     , intent(in)  :: do_output

!  Output arguments
!  ----------------

!  PSD function

   real(KIND=fpk), intent(inout) ::  WY(NNK)

!  Density, Geometric Cross-section and Volume

   real(KIND=fpk), intent(out) ::  NDENS, GXSEC, VOLUME

!  Effective Radius and volume

   real(KIND=fpk), intent(out) ::  REFF, VEFF

!  Local variables
!  ---------------

   integer      :: I
   real(KIND=fpk) :: X, XI, G, SUM, Y, A2, DB, DA, B2, DAB, D_PI
 
!  NDISTR = 1, Modified Gamma

   IF ( NDISTR .eq. 1 ) THEN
                                              
      A2=AA/GAM                                                                 
      DB=1D0/BB
      DO 50 I=1,NNK                                                             
         X=YY(I)                                                             
         Y=X**AA                                                                
         X=X*DB
         Y=Y*DEXP(-A2*(X**GAM))                                                 
         WY(I)=WY(I)*Y                                                       
   50 CONTINUE                                                                  

!   NDISTR = 2, Log normal

   ELSE IF ( NDISTR .eq. 2 ) THEN
  
      DA=1D0/AA                                                                 
      DO 150 I=1,NNK                                                            
         X=YY(I)                                                                
         Y=DLOG(X*DA)                                                          
         Y=DEXP(-Y*Y*0.5D0/BB)/X                                             
         WY(I)=WY(I)*Y                                                          
  150 CONTINUE            
                                                      
!   NDISTR = 3, Power law

   ELSE IF ( NDISTR .eq. 3 ) THEN

      DO 250 I=1,NNK                                                            
         X=YY(I)                                                                
         WY(I)=WY(I)/(X*X*X)                                                 
  250 CONTINUE
                                                                  
!   NDISTR = 4, Gamma

   ELSE IF ( NDISTR .eq. 4 ) THEN

      B2=(1D0-3D0*BB)/BB                                                        
      DAB=1D0/(AA*BB)                                                          
      DO 350 I=1,NNK                                                            
         X=YY(I)                                                                
         X=(X**B2)*DEXP(-X*DAB)                                                 
         WY(I)=WY(I)*X                                                       
  350 CONTINUE  
                                                                
!   NDISTR = 5, Modified Power law

   ELSE IF ( NDISTR .eq. 5 ) THEN

      DO 370 I=1,NNK
         X=YY(I)
         IF (X.LE.R1) WY(I)=WY(I)
         IF (X.GT.R1) WY(I)=WY(I)*(X/R1)**BB
  370 CONTINUE

   ENDIF

!  Print out

      if ( do_output .and. ndistr .eq. 1 ) then
       PRINT 1001,AA,BB,GAM    
 1001  FORMAT('MODIFIED GAMMA DISTRIBUTION, alpha=',F6.4,'  r_c=', &
       F6.4,'  gamma=',F6.4)
      endif
      if ( do_output .and. NDISTR .eq. 2 ) then
       PRINT 1002,AA,BB    
 1002  FORMAT('LOG-NORMAL DISTRIBUTION, r_g=',F8.4, &
            '  [ln(sigma_g)]**2=', F6.4)        
      endif
      if ( do_output .and. NDISTR .eq. 3) then
       PRINT 1003                                                               
 1003  FORMAT('POWER LAW DISTRIBUTION OF HANSEN & TRAVIS 1974')
      endif
      if ( do_output .and. ndistr .eq. 4) then
       PRINT 1004,AA,BB                                                         
 1004  FORMAT ('GAMMA DISTRIBUTION,  a=',F6.3,'  b=',F6.4)
      endif
      if ( do_output .and. ndistr .eq. 5 ) then
       PRINT 1005,BB
 1005  FORMAT ('MODIFIED POWER LAW DISTRIBUTION,  alpha=',D10.4)
      endif

!  Normalize, compute Geometric Cross-section G

      d_pi = dacos(-1.0d0)
      SUM=0D0
      DO 450 I=1,NNK
         SUM=SUM+WY(I)
  450 CONTINUE
      NDENS = SUM

      SUM=1D0/SUM
      DO 500 I=1,NNK
         WY(I)=WY(I)*SUM
  500 CONTINUE
      G=0D0
      DO 550 I=1,NNK
         X=YY(I)
         G=G+X*X*WY(I)
  550 CONTINUE
      GXSEC  = G * d_pi

!  Compute REFF and VEFF

      REFF=0D0
      DO 600 I=1,NNK
         X=YY(I)
         REFF=REFF+X*X*X*WY(I)
  600 CONTINUE
      VOLUME = (4.0d0/3.0d0) * reff * d_pi
      REFF=REFF/G
      VEFF=0D0
      DO 650 I=1,NNK
         X=YY(I)
         XI=X-REFF
         VEFF=VEFF+XI*XI*X*X*WY(I)
  650 CONTINUE
      VEFF=VEFF/(G*REFF*REFF)

!  Finish

   return
end subroutine distrb

 
!  COMPUTATION OF R1 AND R2 FOR A POWER LAW SIZE DISTRIBUTION WITH
!  EFFECTIVE RADIUS A AND EFFECTIVE VARIANCE B
 
      SUBROUTINE POWER (A,B,R1,R2)

!  Old
!      IMPLICIT REAL*8 (A-H,O-Z)

!  New
      real(KIND=fpk), intent(in)  :: A, B
      real(KIND=fpk), intent(out) :: R1, R2
      real(KIND=fpk) :: AA, BB, AX, BX

      AA=A
      BB=B
      AX=1D-5
      BX=A-1D-5
      CALL ZEROIN (AX,BX,0D0,AA,BB,R1)
      R2=(1D0+B)*2D0*A-R1
      RETURN
      END SUBROUTINE POWER
  
      SUBROUTINE ZEROIN (AX,BX,TOL,APAR,BPAR,ZEROIN_OUT)

!  Old
!      IMPLICIT REAL*8 (A-H,O-Z)

!  New
      real(KIND=fpk), intent(in)  :: AX, BX, TOL, APAR, BPAR
      real(KIND=fpk), intent(out) :: ZEROIN_OUT
      real(KIND=fpk) :: EPS, TOL1, A, B, FA, FB, C, FC, D, E, XM, S, P, Q, R

      EPS=1D0
   10 EPS=0.5D0*EPS
      TOL1=1D0+EPS
      IF (TOL1.GT.1D0) GO TO 10
      A=AX     ! 15
      B=BX
      FA=(((1D0+BPAR)*2D0*APAR-A)-A)/DLOG(((1D0+BPAR)*2D0*APAR-A)/A)-APAR
      FB=(((1D0+BPAR)*2D0*APAR-B)-B)/DLOG(((1D0+BPAR)*2D0*APAR-B)/B)-APAR
   20 C=A
      FC=FA
      D=B-A
      E=D
   30 IF (DABS(FC).GE.DABS(FB)) GO TO 40
      A=B      ! 35
      B=C
      C=A
      FA=FB
      FB=FC
      FC=FA
   40 TOL1=2D0*EPS*DABS(B)+0.5D0*TOL
      XM=0.5D0*(C-B)
      IF (DABS(XM).LE.TOL1) GO TO 90
      IF (FB.EQ.0D0) GO TO 90                ! 44
      IF (DABS(E).LT.TOL1) GO TO 70          ! 45
      IF (DABS(FA).LE.DABS(FB)) GO TO 70     ! 46
      IF (A.NE.C) GO TO 50                   ! 47
      S=FB/FA                                ! 48
      P=2D0*XM*S
      Q=1D0-S
      GO TO 60
   50 Q=FA/FC
      R=FB/FC
      S=FB/FA
      P=S*(2D0*XM*Q*(Q-R)-(B-A)*(R-1D0))
      Q=(Q-1D0)*(R-1D0)*(S-1D0)
   60 IF (P.GT.0D0) Q=-Q
      P=DABS(P)
      IF ((2D0*P).GE.(3D0*XM*Q-DABS(TOL1*Q))) GO TO 70
      IF (P.GE.DABS(0.5D0*E*Q)) GO TO 70     ! 64
      E=D                                    ! 65
      D=P/Q
      GO TO 80
   70 D=XM
      E=D
   80 A=B
      FA=FB
      IF (DABS(D).GT.TOL1) B=B+D
      IF (DABS(D).LE.TOL1) B=B+DSIGN(TOL1,XM)
      FB = (((1D0+BPAR)*2D0*APAR-B)-B)/DLOG(((1D0+BPAR)*2D0*APAR-B)/B)-APAR
      IF ((FB*(FC/DABS(FC))).GT.0D0) GO TO 20
      GO TO 30                                 ! 85
   90 ZEROIN_OUT=B
      RETURN
      END SUBROUTINE ZEROIN
  


SUBROUTINE sizedis_nod &
    ( idis, par, numradius, radius, nwithr, message, failnod )

!************************************************************************
!*  Calculate the size distribution n(r) for the numr radius values     *
!*  contained in array r and return the results through the array nwithr*
!*  The size distributions are normalized such that the integral over   *
!*  all r is equal to one.                                              *
!************************************************************************

!  modules

  USE tmat_parameters,    ONLY : d_zero, d_half, d_one, d_two, d_three

  IMPLICIT NONE

!* subroutine arguments

  REAL    (KIND=fpk), INTENT (IN)  :: par(3)
  INTEGER         , INTENT (IN)  :: idis, numradius
  REAL    (KIND=fpk), DIMENSION (numradius), INTENT (IN)  :: radius

  REAL    (KIND=fpk), DIMENSION (numradius), INTENT (INOUT) :: nwithr
  CHARACTER*(*)   , INTENT (OUT) :: message
  LOGICAL         , INTENT (OUT) :: failnod

!* local variables

  LOGICAL        :: deriv
  INTEGER        :: i
  REAL  (KIND=fpk) :: pi,r,logr,root2p
  REAL  (KIND=fpk) :: alpha,alpha1,b,logb,arg1,arg2,arg,argsq,r3

  REAL  (KIND=fpk) :: b1,b2,logb1,logb2,rc, gammln, dgammln
  REAL  (KIND=fpk) :: logrg,logsi,logsi_inv,gamma,gamma1,rg
  REAL  (KIND=fpk) :: rmin,rmax,fac1,fac2,aperg,alpha2
  REAL  (KIND=fpk) :: n1, n2, C, logC, logC1, logC2

  CHARACTER*70   :: message_gamma
  LOGICAL        :: fail
  character*1    :: cdis

!  check

  failnod = .FALSE.
  if (idis == 0 ) RETURN
  IF ( IDIS > 8 ) THEN
     failnod = .TRUE.
     message = 'illegal index in sizedis'
     RETURN
  END IF

!  setup

  pi     = dacos(-1.d0)
  root2p = dsqrt(pi+pi)
  deriv  = .FALSE.

!  IDIS = 1 : TWO-PARAMETER GAMMA with alpha and b given

  IF ( idis == 1 ) THEN

      alpha  = par(1)
      b      = par(2)
      alpha1 = alpha + d_one
      logb   = LOG(b)
      CALL gammafunction ( alpha1, deriv, gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240
      logC  = alpha1*logb - gammln

      DO i = 1, numradius
        r  = radius(i)
        logr = LOG(r)
        arg1 = logC + alpha*logr
        nwithr(i) = EXP ( arg1 - b*r )
      END DO

!  IDIS = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given

  ELSE IF ( idis == 2 ) THEN

      alpha  = d_one/par(2) - d_three
      b      = d_one/(par(1)*par(2))
      alpha1 = alpha + d_one
      logb   = LOG(b)
      CALL gammafunction ( alpha1, deriv, gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = alpha1*logb - gammln

      DO i = 1, numradius
        r    = radius(i)
        logr = LOG(r)
        arg1 = logC + alpha*logr
        nwithr(i) = EXP ( arg1 - b*r )
      END DO

!  IDIS = 3 : BIMODAL GAMMA with equal mode weights

  ELSE IF ( idis == 3 ) THEN

      alpha  = d_one/par(3) - d_three
      b1     = d_one/(par(1)*par(3))
      b2     = d_one/(par(2)*par(3))
      alpha1 = alpha + d_one
      CALL gammafunction ( alpha1, deriv, gammln, dgammln, fail, message_gamma )
      logb1 = LOG(b1)
      logb2 = LOG(b2)
      logC1 = alpha1*logb1 - gammln
      logC2 = alpha1*logb2 - gammln

      DO i = 1, numradius
        r    = radius(i)
        logr = LOG(r)
        arg1 = logC1 + alpha*logr
        arg2 = logC2 + alpha*logr
        n1   = EXP(arg1 - b1*r)
        n2   = EXP(arg2 - b2*r)
        nwithr(i) = d_half * ( n1 + n2 )
      END DO
  
!  4  LOG-NORMAL with rg and sigma given

  ELSE IF ( idis == 4 ) THEN

      logrg = dlog(par(1))
      logsi = dabs(dlog(par(2)))
      logsi_inv = d_one / logsi
      C      = logsi_inv / root2p
      DO i = 1, numradius
        r     = radius(i)
        logr  = LOG(r)
        arg   = ( logr - logrg ) / logsi
        argsq = arg * arg
        nwithr(i) = C * dexp( - d_half * argsq ) / r
      END DO

!  5 LOG-NORMAL with reff and veff given                               *

  ELSE IF ( idis == 5 ) THEN

      alpha1 = d_one + par(2)
      alpha2 = dlog(alpha1)
      rg     = par(1)/(d_one+par(2))**2.5d0
      logrg  = dlog(rg)
      logsi  = dsqrt(alpha2)
      logsi_inv = d_one / logsi
      C         = logsi_inv / root2p
      DO i = 1, numradius
        r     = radius(i)
        logr  = LOG(r)
        arg   = ( logr - logrg ) / logsi
        argsq = arg * arg
        nwithr(i) = C * dexp( - d_half * argsq ) / r
      END DO

!  6 POWER LAW                               *

  ELSE IF ( idis == 6 ) THEN

      alpha = par(1)
      rmin  = par(2)
      rmax  = par(3)
      alpha1 = alpha - d_one
      fac1 = (rmax/rmin)**alpha1
      fac2 = d_one / ( fac1 - d_one )
      C = alpha1 * rmax**alpha1 * fac2
      DO i = 1, numradius
        r     = radius(i)
        if ( (r < rmax) .and. (r > rmin) ) then
          nwithr(i)    = C*r**(-alpha)
        else
          nwithr(i)    = d_zero
        endif
      END DO

!  7 MODIFIED GAMMA with alpha, rc and gamma given

  ELSE IF ( idis == 7 ) THEN

      alpha = par(1)
      rc    = par(2)
      gamma = par(3)
      b     = alpha / (gamma*rc**gamma)
      logb  = LOG(b)
      alpha1 = alpha + d_one
      gamma1 = d_one / gamma
      aperg = alpha1/gamma
      CALL gammafunction ( aperg, deriv, gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      DO i = 1, numradius
        r    = radius(i)
        logr = LOG(r)
        arg1 = logC + alpha*logr
        r3   = b*r ** gamma
        nwithr(i) = EXP ( arg1 - r3 )
      END DO

!  8 MODIFIED GAMMA with alpha, b and gamma given

  ELSE IF ( idis == 8 ) THEN

      alpha = par(1)
      b     = par(2)
      gamma = par(3)
      alpha1 = alpha + d_one
      gamma1 = d_one / gamma
      logb = LOG(b)
      aperg = alpha1/gamma
      CALL gammafunction ( aperg, deriv, gammln, dgammln, fail, message_gamma )
      IF ( fail ) go to 240      
      logC  = dlog(gamma) + aperg*logb - gammln
      DO i = 1, numradius
        r    = radius(i)
        logr = LOG(r)
        arg1 = logC + alpha*logr
        r3   = r ** gamma
        nwithr(i) = EXP ( arg1 - b*r3 )
      END DO

  END IF

!  normal return

  RETURN

!  special return

240  CONTINUE
     failnod = .TRUE.
  write(cdis,'(I1)')idis
  message = message_gamma(1:LEN(message_gamma))//', distribution : '//cdis
  RETURN

END SUBROUTINE sizedis_nod


SUBROUTINE rminmax( idis, par, cutoff, rmin, rmax, message, fail )

  USE tmat_parameters,    ONLY : d_zero, d_half, d_one, d_two, d_three

  IMPLICIT NONE

!  subroutine arguments

  REAL    (KIND=fpk), INTENT (IN)  :: par(3), cutoff
  INTEGER         , INTENT (IN)  :: idis
  CHARACTER*(*)   , INTENT (OUT) :: message
  REAL    (KIND=fpk), INTENT (OUT) :: rmin, rmax
  LOGICAL         , INTENT (OUT) :: fail

!  local variables

  REAL    (KIND=fpk)            :: r(1), nwithr(1),ref,rref,sef,r0,r1,eps
  LOGICAL                     :: failnod

!************************************************************************
!*  Find the integration bounds rmin and rmax for the integration over  *
!*  a size distribution. These bounds are chosen such that the size     *
!*  distribution falls below the user specified cutoff. It is essential *
!*  that the size distribution is normalized such that the integral     *
!*  over all r is equal to one !                                        *
!*  This is programmed rather clumsy and will in the future be changed  *
!************************************************************************

  fail = .FALSE.
  eps = 1.0d-10

  IF (idis == 0) THEN
     rmin= par(1)
     rmax= par(1)
     RETURN
  ELSE IF ( idis == 1)  THEN
     sef  = d_one/SQRT(par(2)+d_three)
     ref  = d_one/(sef*sef*par(2))
     rref = ref
  ELSE IF ( idis == 2)  THEN
     ref = par(1)
     sef = SQRT(par(2))
     rref= ref
  ELSE IF ( idis == 3)  THEN
     sef = SQRT(par(3))
     ref = MAX(par(1),par(2))+sef
     rref= MAX(par(1),par(2))
  ELSE IF ( idis == 4)  THEN
     sef = SQRT(EXP(LOG(par(2))**d_two)-d_one)
     ref = par(1)*(d_one+sef*sef)**0.4d0
     rref= ref
  ELSE IF ( idis == 5)  THEN
     ref = par(1)
     sef = SQRT(ref)
     rref= ref
  ELSE IF ( idis == 6)  THEN
     rmin= par(2)
     rmax= par(3)
     RETURN
  ELSE IF ( idis == 7)  THEN
     ref = par(2)
     sef = d_two*ref
     rref= d_half*ref
  ELSE IF ( idis == 8)  THEN
     ref = (par(1)/(par(2)*par(3)))**par(3)
     sef = d_two*ref
     rref= d_half*ref
  END IF

!************************************************************************
!*  search for a value of r such that the size distribution
!*  is less than the cutoff. Start the search at ref+sef which          *
!*  guarantees that such a value will be found on the TAIL of the       *
!*  distribution.                                                       *
!************************************************************************

  r(1) = ref+sef
  r0   = ref
200 CONTINUE
  CALL sizedis_nod( idis, par, 1, r, nwithr, message, failnod )
  IF ( failnod ) GO TO 899
  IF ( nwithr(1) > cutoff) THEN
     r0   = r(1)
     r(1) = d_two*r(1)
     goto 200
  END IF
  r1 = r(1)

!************************************************************************
!*  Now the size distribution assumes the cutoff value somewhere        *
!*  between r0 and r1  Use bisection to find the corresponding r        *
!************************************************************************

300 CONTINUE
  r(1) = d_half*(r0+r1)
  CALL sizedis_nod( idis, par, 1, r, nwithr, message, failnod )
  IF ( failnod ) GO TO 899
  IF ( nwithr(1) > cutoff) THEN
     r0   = r(1)
  ELSE
     r1 = r(1)
  END IF
  IF ((r1-r0) > eps) GOTO 300
  rmax = d_half*(r0+r1)

!************************************************************************
!*  Search for a value of r on the low end of the size distribution     *
!*  such that the distribution falls below the cutoff. There is no      *
!*  guarantee that such a value exists, so use an extra test to see if  *
!*  the search comes very near to r = 0                                 *
!************************************************************************

  r1 = rref
  r0 = d_zero
400 CONTINUE
  r(1) = d_half*r1
  CALL sizedis_nod( idis, par, 1, r, nwithr, message, failnod )
  IF ( failnod ) GO TO 899
  IF ( nwithr(1) > cutoff) THEN
     r1 = r(1)
     IF (r1 > eps) GOTO 400
  ELSE
     r0 = r(1)
  END IF

!************************************************************************
!*  Possibly the size distribution goes through cutoff between r0       *
!*  and r1 try to find the exact value of r where this happens by       *
!*  bisection.                                                          *
!*  In case there is no solution, the algorithm will terminate soon.    *
!************************************************************************

500 CONTINUE
  r(1) = d_half*(r0+r1)
  CALL sizedis_nod( idis, par, 1, r, nwithr, message, failnod )
  IF ( failnod ) GO TO 899
  IF ( nwithr(1) > cutoff) THEN
     r1 = r(1)
  ELSE
     r0 = r(1)
  END IF
  IF ( (r1-r0) > eps) GOTO 500
  IF (r1 <= eps) THEN
     rmin = d_zero
  ELSE
     rmin = d_half*(r0+r1)
  END IF

! normal return

  RETURN

!  error return

899  CONTINUE
  fail = .TRUE.
  RETURN
END SUBROUTINE rminmax

SUBROUTINE gammafunction &
     ( xarg, do_derivative, gammln, dgammln, gammafail, message )

!************************************************************************
!*  Return the value of the natural logarithm of the gamma function.    *
!*  and its derivative (the Digamma function)                           *
!*  The argument xarg must be real and positive.                        *
!*  This function is documented in :                                    *
!*                                                                      *
!*  W.H. Press et al. 1986, 'Numerical Recipes' Cambridge Univ. Pr.     *
!*  page 157 (ISBN 0-521-30811)                                         *
!*                                                                      *
!*  When the argument xarg is between zero and one, the relation (6.1.4)*
!*  on page 156 of the book by Press is used.                           *
!*                                         V.L. Dolman April 18 1989    *
!************************************************************************

!  modules

  USE tmat_parameters,    ONLY : d_zero, d_half, d_one, d_two

  IMPLICIT NONE

!* subroutine arguments

  REAL    (KIND=fpk), INTENT (IN)  :: xarg
  LOGICAL         , INTENT (IN)  :: do_derivative
  REAL    (KIND=fpk), INTENT (OUT) :: gammln, dgammln
  CHARACTER*(*)   , INTENT (OUT) :: message
  LOGICAL         , INTENT (OUT) :: gammafail

!* local parameters and data

  REAL    (KIND=fpk), PARAMETER :: gammaf_eps = 1.d-10
  REAL    (KIND=fpk), PARAMETER :: gammaf_fpf = 5.5D0
  REAL    (KIND=fpk) :: cof(6),stp, c0
  DATA cof /    76.18009172947146D0,    &
               -86.50532032941677D0,    &
                24.01409824083091D0,    &
                -1.231739572450155D0,   &
                 0.1208650973866179D-2, &
                -0.5395239384953D-5      /
  DATA c0  / 1.000000000190015D0  /
  DATA stp / 2.5066282746310005D0 /

!* local variables

  INTEGER        :: j
  REAL  (KIND=fpk) :: x,xx,xxx,tmp,x1,x2,logtmp,pi,ser,dser,gtmp,dgtmp,pix

!*  initialize output

  message   = ' '
  gammln    = d_zero
  dgammln   = d_zero
  gammafail = .FALSE.

!* check for bad input

  IF (xarg <= d_zero) THEN
     message = ' gammafunction: called with negative argument xarg'
     gammafail = .TRUE.
     RETURN
  END IF
  IF (ABS(xarg-d_one) < gammaf_eps) THEN
     message = ' gammafunction: argument too close to one '
     gammafail = .TRUE.
     RETURN
  END IF

!*  set up

  pi = 4.0d0 * ATAN(d_one)
  IF (xarg .ge. d_one) THEN
     xxx = xarg
  ELSE
     xxx = xarg + d_two
  END IF

!* Numerical Recipes stuff

  xx = xxx - d_one
  x1 = xx  + gammaf_fpf
  x2 = xx  + d_half

  logtmp = LOG(x1)
  tmp = x2*logtmp-x1

  ser = c0
  x = xx
  DO j =1, 6
     x = x + d_one
     ser = ser+cof(j)/x
  END DO
  gtmp = tmp + LOG(stp*ser)

!* derivative of gammln

  IF ( do_derivative ) THEN
     dser = d_zero
     x = xx
     DO j = 1, 6
        x = x+d_one
        dser = dser+cof(j)/x/x
     END DO
     dgtmp = logtmp - (5.0d0/x1) - (dser/ser)
  END IF

!* assign output

  IF ( do_derivative ) THEN
     IF  (xarg > d_one) THEN
        gammln  = gtmp
        dgammln = dgtmp
     ELSE
        pix = pi*(d_one-xarg)
        gammln  = LOG(pix/SIN(pix))-gtmp
        dgammln = - dgtmp - (pi*COS(pix)/SIN(pix)) - (d_one/(d_one-xarg))
     END IF
  ELSE
     IF  (xarg > d_one) THEN
        gammln  = gtmp
     ELSE
        pix = pi*(d_one-xarg)
        gammln  = LOG(pix/SIN(pix))-gtmp
     END IF
  END IF

  RETURN
end SUBROUTINE gammafunction

!  End module

end module tmat_distributions
