MODULE RTSMie_distributions_m

use RTSMie_parameters_m

private
public sizedis, sizedis_plus, Mie_gauleg, rminmax

!  Contains the following routines for the Extended calculation

!     sizedist_plus

!  Contains the following routines for the Regular calculation

!     sizedist
!     sizedist_nod
!     gammafunction
!     Mie_gauleg
!     rminmax

contains

SUBROUTINE sizedis_plus                                         &
    ( max_Mie_distpoints, idis, par, Gderiv, radius, numradius, &
      nwithr, nwithr_d, message, faild )

!************************************************************************
!*  Calculate the size distribution n(r) for the numr radius values     *
!*  contained in array r and return the results through the array nwithr*
!*  The size distributions are normalized such that the integral over   *
!*  all r is equal to one.                                              *
!************************************************************************
!  modules

  USE RTSMie_parameters_m, ONLY : dp, d_zero, d_half, d_one, d_two, d_three

  IMPLICIT NONE

!* subroutine arguments

  INTEGER          , INTENT (IN)  :: max_Mie_distpoints
  REAL    (KIND=dp), INTENT (IN)  :: par(3)
  LOGICAL          , INTENT (IN)  :: Gderiv
  INTEGER          , INTENT (IN)  :: idis, numradius
  CHARACTER*(*)    , INTENT (OUT) :: message
  LOGICAL          , INTENT (OUT) :: faild

  REAL    (KIND=dp), DIMENSION (max_Mie_distpoints),   INTENT (IN)  :: radius
  REAL    (KIND=dp), DIMENSION (max_Mie_distpoints),   INTENT (OUT) :: nwithr
  REAL    (KIND=dp), DIMENSION (max_Mie_distpoints,3), INTENT (OUT) :: nwithr_d

!* local variables

  LOGICAL         :: deriv(3)
  INTEGER         :: i
  REAL  (KIND=dp) :: pi,r,logr,root2p
  REAL  (KIND=dp) :: alpha,alpha1,b,logb,arg1,arg2,arg,argsq,r3
  REAL  (KIND=dp) :: b1,b2,b11,b13,b22,b23,logb1,logb2,rc
  REAL  (KIND=dp) :: logrg,logsi,logsi_inv,fac_d1,gamma,gamma1,rg
  REAL  (KIND=dp) :: rmin,rmax,fac1,fac2,aperg
  REAL  (KIND=dp) :: alpha2, fac_d2a
  REAL  (KIND=dp) :: n1, n2, n1_d1, n1_d3, n2_d2, n2_d3

!  redundant variables
!  REAL  (KIND=dp) :: sigfac, logC1_d2

  REAL  (KIND=dp) :: C, logC, logC_d1, logC_d2, logC_d3
  REAL  (KIND=dp) :: logC1, logC2, logC1_d1, logC1_d3, logC2_d2, logC2_d3

  REAL    (KIND=dp) :: gammln, dgammln
  CHARACTER*70      :: message_gamma
  LOGICAL           :: fail
  character*1       :: cdis

!  check

   faild = .FALSE. ; message = ' '
   if (idis == 0 ) RETURN
   IF ( IDIS > 8 ) THEN
      faild = .TRUE.
      message = 'illegal index in sizedis'
      RETURN
   END IF

!  setup

   pi     = dacos(-1.d0)
   root2p = dsqrt(pi+pi)

!mick fix 6/19/2014 - initialize output

   nwithr   = d_zero
   nwithr_d = d_zero

!  Derivative flags

   deriv = .false.
   if ( Gderiv ) then
      if ( idis.eq.1 .or. idis.eq.2 .or. idis.eq.4 .or. idis.eq.5 ) deriv(1:2) = .true.
      if ( idis.eq.3 .or. idis.eq.6 .or. idis.eq.7 .or. idis.eq.8 ) deriv(1:3) = .true.
   endif

!  IDIS = 1 : TWO-PARAMETER GAMMA with alpha and b given

  IF ( idis == 1 ) THEN

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
      rg     = par(1)/(d_one+par(2))**2.5_dp
      logrg  = dlog(rg)
      logsi  = dsqrt(alpha2)
      logsi_inv = d_one / logsi
      C         = logsi_inv / root2p
      IF ( deriv(1) .and. deriv(2) ) then
        logC_d2 = - d_half / alpha2 / alpha1
        fac_d1  =   logsi_inv / par(1)
        fac_d2a  =  - 2.5_dp * logsi_inv / alpha1
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

END SUBROUTINE sizedis_plus


SUBROUTINE sizedis                                      &
    ( max_Mie_distpoints, idis, par, radius, numradius, &
      nwithr, message, faild )

!************************************************************************
!*  Calculate the size distribution n(r) for the numr radius values     *
!*  contained in array r and return the results through the array nwithr*
!*  The size distributions are normalized such that the integral over   *
!*  all r is equal to one.                                              *
!************************************************************************
!  modules

  USE RTSMie_parameters_m, ONLY : dp, d_zero, d_half, d_one, d_two, d_three

  IMPLICIT NONE

!* subroutine arguments

  INTEGER          , INTENT (IN)  :: max_Mie_distpoints
  REAL    (KIND=dp), INTENT (IN)  :: par(3)
  INTEGER          , INTENT (IN)  :: idis, numradius
  CHARACTER*(*)    , INTENT (OUT) :: message
  LOGICAL          , INTENT (OUT) :: faild

  REAL    (KIND=dp), DIMENSION (max_Mie_distpoints),   INTENT (IN)  :: radius
  REAL    (KIND=dp), DIMENSION (max_Mie_distpoints),   INTENT (OUT) :: nwithr

!* local variables

  INTEGER         :: i
  REAL  (KIND=dp) :: pi,r,logr,root2p
  REAL  (KIND=dp) :: alpha,alpha1,b,logb,arg1,arg2,arg,argsq,r3
  REAL  (KIND=dp) :: b1,b2,logb1,logb2,rc
  REAL  (KIND=dp) :: logrg,logsi,logsi_inv,gamma,gamma1,rg
  REAL  (KIND=dp) :: R1,R2,fac1,fac2,aperg, alpha2
  REAL  (KIND=dp) :: n1, n2, C, logC, logC1, logC2

  REAL    (KIND=dp) :: gammln, dummy
  CHARACTER*70      :: message_gamma
  LOGICAL           :: fail
  character*1       :: cdis

!  check

      faild = .FALSE.

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
      rg     = par(1)/(d_one+par(2))**2.5_dp
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
      R1  = par(2)
      R2  = par(3)
      alpha1 = alpha - d_one
      fac1 = (R2/R1)**alpha1
      fac2 = d_one / ( fac1 - d_one )
      C = alpha1 * R2**alpha1 * fac2
      DO i = 1, numradius
         r     = radius(i)
         if ( (r < R2) .and. (r > R1) ) then
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

END SUBROUTINE sizedis



SUBROUTINE sizedis_nod &
    ( idis, par, numradius, radius, nwithr, message, failnod )

!************************************************************************
!*  Calculate the size distribution n(r) for the numr radius values     *
!*  contained in array r and return the results through the array nwithr*
!*  The size distributions are normalized such that the integral over   *
!*  all r is equal to one.                                              *
!************************************************************************
!  modules

  USE RTSMie_parameters_m, ONLY : dp, d_zero, d_half, d_one, d_two, d_three

  IMPLICIT NONE

!* subroutine arguments

  REAL    (KIND=dp), INTENT (IN)  :: par(3)
  INTEGER          , INTENT (IN)  :: idis, numradius
  CHARACTER*(*)    , INTENT (OUT) :: message
  LOGICAL          , INTENT (OUT) :: failnod

  REAL    (KIND=dp), DIMENSION (numradius), INTENT (IN)  :: radius
  REAL    (KIND=dp), DIMENSION (numradius), INTENT (OUT) :: nwithr

!* local variables

  LOGICAL         :: deriv
  INTEGER         :: i
  REAL  (KIND=dp) :: pi,r,logr,root2p
  REAL  (KIND=dp) :: alpha,alpha1,b,logb,arg1,arg2,arg,argsq,r3

  REAL  (KIND=dp) :: b1,b2,logb1,logb2,rc, gammln, dgammln
  REAL  (KIND=dp) :: logrg,logsi,logsi_inv,gamma,gamma1,rg
  REAL  (KIND=dp) :: R1,R2,fac1,fac2,aperg,alpha2
  REAL  (KIND=dp) :: n1, n2, C, logC, logC1, logC2

  CHARACTER*70      :: message_gamma
  LOGICAL           :: fail
  character*1       :: cdis

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
!         write(*,*)i,logr,arg1
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
      rg     = par(1)/(d_one+par(2))**2.5_dp
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
      R1  = par(2)
      R2  = par(3)
      alpha1 = alpha - d_one
      fac1 = (R2/R1)**alpha1
      fac2 = d_one / ( fac1 - d_one )
      C = alpha1 * R2**alpha1 * fac2
      DO i = 1, numradius
        r     = radius(i)
        if ( (r < R2) .and. (r > R1) ) then
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


SUBROUTINE rminmax( idis, par, cutoff, R1, R2, message, fail )

  USE RTSMie_parameters_m, ONLY : dp, d_zero, d_half, d_one, d_two, d_three

  IMPLICIT NONE

!  subroutine arguments

  REAL    (KIND=dp), INTENT (IN)  :: par(3), cutoff
  INTEGER          , INTENT (IN)  :: idis
  CHARACTER*(*)    , INTENT (OUT) :: message
  REAL    (KIND=dp), INTENT (OUT) :: R1, R2
  LOGICAL          , INTENT (OUT) :: fail

!  local variables

  REAL    (KIND=dp)            :: r(1), nwithr(1),ref,rref,sef,r0,r11,eps
  LOGICAL                      :: failnod

!************************************************************************
!*  Find the integration bounds R1 and R2 for the integration over  *
!*  a size distribution. These bounds are chosen such that the size     *
!*  distribution falls below the user specified cutoff. It is essential *
!*  that the size distribution is normalized such that the integral     *
!*  over all r is equal to one !                                        *
!*  This is programmed rather clumsy and will in the future be changed  *
!************************************************************************

  fail = .FALSE.
  eps = 1.0E-10_dp

  IF (idis == 0) THEN
     R1= par(1)
     R2= par(1)
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
     ref = par(1)*(d_one+sef*sef)**0.4_dp
     rref= ref
  ELSE IF ( idis == 5)  THEN
     ref = par(1)
     sef = SQRT(ref)
     rref= ref
  ELSE IF ( idis == 6)  THEN
     R1= par(2)
     R2= par(3)
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
  r11 = r(1)

!************************************************************************
!*  Now the size distribution assumes the cutoff value somewhere        *
!*  between r0 and r1  Use bisection to find the corresponding r        *
!************************************************************************

300 CONTINUE
  r(1) = d_half*(r0+r11)
  CALL sizedis_nod( idis, par, 1, r, nwithr, message, failnod )
  IF ( failnod ) GO TO 899
  IF ( nwithr(1) > cutoff) THEN
     r0   = r(1)
  ELSE
     r11 = r(1)
  END IF
  IF ((r11-r0) > eps) GOTO 300
  R2 = d_half*(r0+r11)

!************************************************************************
!*  Search for a value of r on the low end of the size distribution     *
!*  such that the distribution falls below the cutoff. There is no      *
!*  guarantee that such a value exists, so use an extra test to see if  *
!*  the search comes very near to r = 0                                 *
!************************************************************************

  r11 = rref
  r0 = d_zero
400 CONTINUE
  r(1) = d_half*(r11+R0)           ! Lost R11, bug 16 January 2013
  CALL sizedis_nod( idis, par, 1, r, nwithr, message, failnod )
  IF ( failnod ) GO TO 899
  IF ( nwithr(1) > cutoff) THEN
     r11 = r(1)
     IF (r11 > eps) GOTO 400
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
  r(1) = d_half*(r0+r11)
  CALL sizedis_nod( idis, par, 1, r, nwithr, message, failnod )
  IF ( failnod ) GO TO 899
  IF ( nwithr(1) > cutoff) THEN
     r11 = r(1)
  ELSE
     r0 = r(1)
  END IF
  IF ( (r11-r0) > eps) GOTO 500
  IF (r11 <= eps) THEN
     R1 = d_zero
  ELSE
     R1 = d_half*(r0+r11)
  END IF

! normal return

  RETURN

!  error return

899  CONTINUE
  fail = .TRUE.
  RETURN
END SUBROUTINE rminmax

!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************

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

  USE RTSMie_parameters_m, ONLY : dp, d_zero, d_half, d_one, d_two

  IMPLICIT NONE

!* subroutine arguments

  REAL    (KIND=dp), INTENT (IN)  :: xarg
  LOGICAL          , INTENT (IN)  :: do_derivative
  REAL    (KIND=dp), INTENT (OUT) :: gammln, dgammln
  CHARACTER*(*)    , INTENT (OUT) :: message
  LOGICAL          , INTENT (OUT) :: gammafail

!* local parameters and data

  REAL    (KIND=dp), PARAMETER :: gammaf_eps = 1.d-10
  REAL    (KIND=dp), PARAMETER :: gammaf_fpf = 5.5D0
  REAL    (KIND=dp) :: cof(6),stp, c0
  DATA cof /    76.18009172947146D0,    &
               -86.50532032941677D0,    &
                24.01409824083091D0,    &
                -1.231739572450155D0,   &
                 0.1208650973866179D-2, &
                -0.5395239384953D-5      /
  DATA c0  / 1.000000000190015D0  /
  DATA stp / 2.5066282746310005D0 /

!* local variables

  INTEGER         :: j
  REAL  (KIND=dp) :: x,xx,xxx,tmp,x1,x2,logtmp,pi,ser,dser,gtmp,dgtmp,pix

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

  pi = 4.0_dp * ATAN(d_one)
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
     dgtmp = logtmp - (5.0_dp/x1) - (dser/ser)
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

!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************

SUBROUTINE mie_gauleg(maxn,n,x1,x2,x,w)

  USE RTSMie_parameters_m, ONLY : dp

  IMPLICIT NONE

!* subroutine arguments

  INTEGER          , INTENT (IN)  :: maxn, n
  REAL    (KIND=dp), INTENT (IN)  :: x1,x2
  REAL    (KIND=dp), INTENT (OUT) :: x(maxn),w(maxn)

  INTEGER            :: i, m, j
  REAL    (KIND=dp)  ::  xm,xl,p1,p2,p3,pp,z,z1,eps

  eps=3.0e-14_dp
  m=(n+1)/2
  xm=0.5_dp*(x2+x1)
  xl=0.5_dp*(x2-x1)

  DO i=1,m
     z=COS(3.1415926540_dp*(DBLE(i)-0.250_dp)/(n+0.50_dp))
1  CONTINUE
     p1=1.0_dp
     p2=0.0_dp

     DO j=1,n
        p3=p2
        p2=p1
        p1=((2.0_dp*DBLE(j)-1.0_dp)*z*p2-(DBLE(j)-1.0_dp)*p3)/DBLE(j)
     END DO

     pp=n*(z*p1-p2)/(z*z-1.0_dp)
     z1=z
     z=z1-p1/pp
     if(dabs(z-z1) > eps) GO TO 1
     x(i)=xm-xl*z
     x(n+1-i)=xm+xl*z
     w(i)=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
     w(n+1-i)=w(i)
  END DO

  RETURN
END SUBROUTINE mie_gauleg

END MODULE RTSMie_distributions_m
