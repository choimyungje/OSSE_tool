MODULE PY_distributions_m

private  gammafunction
public   sizedis

!  Contains the following routines for the Regular calculation
!     sizedist
!     gammafunction

contains

SUBROUTINE sizedis                                      &
    ( max_Mie_distpoints, idis, par, radius, numradius, &
      nwithr, message, faild )

!************************************************************************
!*  Calculate the size distribution n(r) for the numr radius values     *
!*  contained in array r and return the results through the array nwithr*
!*  The size distributions are normalized such that the integral over   *
!*  all r is equal to one.                                              *
!************************************************************************

  IMPLICIT NONE

   INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

! Numbers such as 1, 2, pie, etc.......

   REAL (KIND=dp), PARAMETER :: d_zero  = 0.0_dp
   REAL (KIND=dp), PARAMETER :: d_one   = 1.0_dp
   REAL (KIND=dp), PARAMETER :: d_two   = 2.0_dp
   REAL (KIND=dp), PARAMETER :: d_three = 3.0_dp
   REAL (KIND=dp), PARAMETER :: d_half  = 0.5_dp

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

!

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

  IMPLICIT NONE

   INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

! Numbers such as 1, 2, pie, etc.......

   REAL (KIND=dp), PARAMETER :: d_zero  = 0.0_dp
   REAL (KIND=dp), PARAMETER :: d_one   = 1.0_dp
   REAL (KIND=dp), PARAMETER :: d_two   = 2.0_dp
   REAL (KIND=dp), PARAMETER :: d_half  = 0.5_dp

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

!

END MODULE PY_distributions_m
