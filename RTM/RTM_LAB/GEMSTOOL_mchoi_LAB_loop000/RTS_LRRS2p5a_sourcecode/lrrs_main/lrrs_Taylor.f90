
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
! #    Taylor series, small number expansions:                  #
! #                                                             #
! #             TAYLOR_SERIES_1                                 #
! #             TAYLOR_SERIES_2                                 #
! #             TAYLOR_SERIES_L_1                               #
! #             TAYLOR_SERIES_L_2a                              #
! #             TAYLOR_SERIES_L_2b                              #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5, Taylor series module.
!     - Replaces the Former LRRS_SMALLNOS.f90 Module (version 2.3)

!    New for this version, but mainly a copy of LIDORT Version 3.7 module of the same name.
!     - as in LIDORT, fixed-parameter "Taylor_small" is set in "PARS" file.
!     - Recommended values are 3 tfor the Taylor order, 1.0d-05 for Taylor_small

!     - First programmed, 9/10/15 by R. Spurr. RT SOLUTIONS Inc.
!     - Debugging done  , December 2016 (M. Christi), January 2017 (R. Spurr)
!     - This is "Clean Version", see "lrrs_Taylor_final_debugged.f90" for debug details

MODULE lrrs_Taylor_m

   USE LRRS_PARS_m, only : max_Taylor_terms, fpk, zero, one, two

PUBLIC

CONTAINS

subroutine Taylor_Series_1 &
          ( order, eps, delta, udel, sm, mult )

!  Good for homogeneous & particular solution multipliers.
!      Small number expansion to any order up to 4.

!  Note: In LIDORT, this subroutine is applied to quantities of the form
!                [exp(-b*DELTA) - exp(-a*DELTA)]
!            q = -------------------------------
!                            (a-b)
!        where (a-b) has become small to the point of causing instability.
!        Note the positions of "a" and "b" in the numerator are swapped
!        from their positions in the denominator.
!
!        Using the above form for "q", the I/O for the subroutine is as follows:
!        * ORDER = order of last term in the series
!        * EPS   = (a-b)
!        * DELTA = optical thickness (whole or partial layer)
!        * UDEL  = exp(-a*DELTA) (usually UDEL or WDEL)
!        * SM    = a (note: "a" currently applied outside subroutine; thus, SM=1 currently)
!        * MULT  = q

   IMPLICIT NONE

!  arguments

   INTEGER  , INTENT(IN)  :: order
   REAL(FPK), INTENT(IN)  :: eps, delta, udel, sm
   REAL(FPK), INTENT(OUT) :: mult

!  local declarations

   INTEGER   :: mterms, m
   REAL(FPK) :: power, d(0:max_Taylor_terms)

!  Zero output

   mult = zero

!  exp(De) expansion coefficients

   mterms = order + 1 ; d(0) = one
   do m = 1, mterms
      d(m) = delta * d(m-1) / real(m,fpk)
   enddo

!  evaluate multiplier

   mult = d(1) ; power = one
   do m = 2, mterms
      power = power * eps ; mult = mult + power * d(m)
   enddo
   mult = mult * udel * sm

!  Equivalent to the following, for order = 3 (highest power of eps)
!      power   = delta*eps ;  power2  = power*power ; power3 = power2 * power
!      mult = udel * sm * delta *(one + half*power + power2/6.0_fpk + power3/24.0_fpk)
!power = delta*eps
!mult3 = udel * sm * delta *(one + 0.5_fpk*power + (power**2.0_fpk)/6.0_fpk + (power**3.0_fpk)/24.0_fpk)
!mult4 = udel * sm * delta *(one + 0.5_fpk*power + (power**2.0_fpk)/6.0_fpk + (power**3.0_fpk)/24.0_fpk+(power**4.0_fpk)/120.0_fpk)

!  Finish

   return
end subroutine Taylor_Series_1

!

subroutine Taylor_Series_2 &
       ( order, small, eps, y, delta, fac1, fac2, sm, mult )

!  Post-processing Green's function multipliers.
!    Small number expansion to any order

   implicit none

!  Arguments

   integer  , intent(in)  :: order
   real(fpk), intent(in)  :: small, eps, y, delta, fac1, fac2, sm
   real(fpk), intent(out) :: mult

!  Local

   integer    :: mterms, m, j
   real(fpk)  :: d(0:max_Taylor_terms), ac(0:max_Taylor_terms), cc(0:max_Taylor_terms), COF_0
   real(fpk)  :: Term_1(0:max_Taylor_terms),Term_2
   real(fpk)  :: Y1, power, power2

!  Zero output

   mult = zero

!  number of terms. One more if you have a double small-number 

   mterms = order + 1
   if ( abs(y) .lt. small ) mterms = mterms+1

!  Delta powers

   d(0) = one
   do m = 1, mterms
      d(m) = delta * d(m-1) / real(m,fpk)
   enddo

!  Double small number expansion......Return exit

   if ( abs(y) .lt. small ) then
      power = one ; power2 = one
      mult = d(2)
      do m = 3, mterms
         power  = power  * ( eps - y )
         power2 = power - y * power2
         mult = mult + d(m) * power2       
      enddo
      mult = mult * fac1 * sm
      return
   endif

!  Develop ( 1 - x/s )^-1  for small x

   y1 = one / y
   ac(0) = one
   do m = 1, mterms
      ac(m) = y1 * ac(m-1)
   enddo

!  Develop ( 1 - x/s )^-1 . exp(Dx), for small x

   cc(0) = one
   do m = 1, mterms
      cc(m) = zero
      do j = 0, m
         cc(m) = cc(m) + ac(j)*d(m-j)
      enddo
   enddo

!  Term 1. Use AC and CC

   do m = 0, mterms
      Term_1(m) = fac1 * ac(m) - fac2 * cc(m)
   enddo

!  Term 2, then Check 0 coefficient is zero

   Term_2 = fac1 - fac2
   COF_0  = Term_1(0) - Term_2

!  Final answer

   power = one ; mult = Term_1(1)
   do m = 2, mterms
      power = eps * power
      mult = mult + power * Term_1(m)
   enddo
   mult = mult * sm * y1

return

!  Done

   return
end subroutine Taylor_Series_2

!

subroutine Taylor_Series_L_1 ( order, eps, delta, ddot, kdot, Ldot, uterm, sm, L_mult )

!  Small number expansion for derivatives of Series 1 quantities
!     sm is required, but result is NOT SCALED

!   L_HMULT --> sm = user-secant   ,    Ldot = zero
!   L_EMULT --> sm = user_secant   ,    Ldot = zero
!   L_GMULT --> sm = average_secant,    kdot/Ldot present

   implicit none

!  arguments

   INTEGER  , INTENT(IN)  :: order
   REAL(FPK), INTENT(IN)  :: eps, delta, ddot, kdot, Ldot, uterm, sm
   REAL(FPK), INTENT(OUT) :: L_mult

!  local declarations

   INTEGER   :: mterms, m, m1, m2
   REAL(FPK) :: power, d(0:max_Taylor_terms), series1, series2, series3

!  Zero output

   L_mult = zero

!  exp(De) expansion coefficients

   mterms = order + 2 ; d(0) = one
   do m = 1, mterms
      d(m) = delta * d(m-1) / real(m,fpk)
   enddo

!  Develop  series
!  Series 3 absent for HMULT/EMULT, only present for GMULT

   power = one
   series1 = d(0) - d(1)*sm
   series2 = d(2) - d(1)*delta
   series3 = d(2)
   do m = 1, order
      m1 = m + 1 ; m2 = m1 + 1
      power = power * eps
      series1 = series1 + power * (d(m)  - d(m1)*sm   )
      series2 = series2 + power * (d(m2) - d(m1)*delta)
      series3 = series3 + power *  d(m2)
   enddo

!  Final

   L_mult = ( ddot*series1 - Ldot*series3 + kdot*series2 ) * uterm

return

!  done

end subroutine Taylor_Series_L_1

!

subroutine Taylor_Series_L_2a &
       ( order, eps, y, delta, deldot, lamdot, kdot, wudel, sm, &
         lam, fac1, fac2, deldot2, kdot2, L_mult )

!  For L(G+-)  &  L(G--)(flip case)
!     rewritten December 2016 by M. Christi

   implicit none

!  Arguments

   integer  , intent(in)  :: order
   REAL(FPK), intent(in)  :: eps, y, delta, deldot, lamdot, kdot, wudel, sm, &
                             lam, fac1, fac2, deldot2, kdot2
   REAL(FPK), intent(out) :: L_mult

!  Local

   integer   :: mterms1, mterms2, m, j

   real(fpk) :: c1,c2,c3,sigma,series1,series2,power
   real(fpk) :: a(0:max_Taylor_terms),b(0:max_Taylor_terms),&
                g(0:max_Taylor_terms),fact(0:max_Taylor_terms+1)

   real(fpk) :: dc1dx,dc2dx,dc3dx,depsdx,dsigdx
   real(fpk) :: dadx(0:max_Taylor_terms),dbdx(0:max_Taylor_terms),&
                dgdx(0:max_Taylor_terms),ABC(0:max_Taylor_terms),&
                dABCdx(0:max_Taylor_terms)

!  Zero it

   L_mult = zero

!  Number of terms

   mterms2 = order + 2
   mterms1 = order + 1  ! Need this for same order of accuracy as Radiances

!mick check 
!write(*,*)
!write(*,*) 'From Taylor_Series_L_2a'
!write(*,*) 'inputs:'
!write(*,*) 'order      = ',order   !R & M: taylor_order
!write(*,*) 'eps        = ',eps     !R    : rhom
!write(*,*) 'y (keigen) = ',y       !R    : keigen(aa,n)
!write(*,*) 'delta      = ',delta   !R & M: deltaus
!write(*,*) 'deldot     = ',deldot  !R    : l_deltaus(q) OR zero
!write(*,*) 'lamdot     = ',lamdot  !R & M: l_secbar
!write(*,*) 'kdot       = ',kdot    !R    : l_keigen(q,aa,n) OR zero
!write(*,*) 'wudel      = ',wudel   !R    : ktrans(aa,n)*udel
!write(*,*) 'sm         = ',sm      !R & M: dmu
!write(*,*) 'lam        = ',lam     !M    : asource
!write(*,*) 'fac1       = ',fac1    !M    : fac1
!write(*,*) 'fac2       = ',fac2    !M    : fac2
!write(*,*) 'deldot2    = ',deldot2 !M    : l_deltaus(q)
!write(*,*) 'kdot2      = ',kdot2   !M    : l_keigen(q,aa,n)

   depsdx = lamdot - kdot2
   sigma  = one/(lam + sm)
   dsigdx = -sigma**2*lamdot

   c1 = fac1 - fac2 !1 - Wdel*Udel
   c2 = fac2        !    Wdel*Udel
   dc2dx = -c2*(lamdot*delta + (lam + sm)*deldot2)
   dc1dx = -dc2dx

   fact(0) = one
   do m = 1, mterms2+1
      fact(m) = fact(m-1)*real(m,fpk)
   enddo

   a(0) = one   ; dadx(0) = zero
   a(1) = sigma ; dadx(1) = dsigdx 
   b(0) = delta ; dbdx(0) = deldot2
   g(0) = delta ; dgdx(0) = deldot2
   do m = 1, mterms2
      a(m) = sigma**m
      b(m) = delta**(m+1)/fact(m+1)
      dadx(m) = a(m-1)*dsigdx*real(m,fpk)
      dbdx(m) = b(m-1)*deldot2
      g(m) = zero ; dgdx(m) = zero
      do j = 0, m
         g(m) = g(m) + a(j)*b(m-j)
         dgdx(m) = dgdx(m) + (dadx(j)*b(m-j) + a(j)*dbdx(m-j))
      enddo
   enddo

   c3 = a(mterms1)*eps**mterms1
   dc3dx = dadx(mterms1)*eps**mterms1 + a(mterms1)*real(mterms1,fpk)*eps**(mterms1-1)*depsdx

   do m = 1, mterms2
      !ABC(m-1) = AA(m-1) - BB(m-1) - CC(m-1) for m=0,mterms1
      ABC(m-1) = c1*a(m) - c2*g(m-1) - c2*c3*b(m-1)
      !Derivative of ABC(m-1) = dAAdx(m-1) - dBBdx(m-1) - dCCdx(m-1) for m=1,mterms1 
      dABCdx(m-1) = (dc1dx*a(m) + c1*dadx(m)) - (dc2dx*g(m-1) + c2*dgdx(m-1)) &
                  - (dc2dx*c3 + c2*dc3dx)*b(m-1) + c2*c3*dbdx(m-1)
   enddo

   !for m=0
   power = one
   series1 = ABC(0)
   series2 = dABCdx(0) + ABC(1)*depsdx

   !for m=1,mterms1-1
   do m = 1, mterms1
      power = power*eps
      series1 = series1 + ABC(m)*power !S
      if (m<mterms1) series2 = series2 + (dABCdx(m) + ABC(m+1)*real(m+1,fpk)*depsdx)*power !dSdx
   enddo
   series2 = series2 + dABCdx(mterms1)*power*eps

!  Final answer

   L_mult = sm*(dsigdx*series1 + sigma*series2)

!  Done

   return
end subroutine Taylor_Series_L_2a

!

subroutine Taylor_Series_L_2b &
       ( order, small, eps, y, delta, deldot, lamdot, kdot, wdel, udel, sm, &
         lam, fac1, fac2, deldot2, kdot2, L_mult )

!  For L(G+-)(flip case)  &  L(G--)
!     rewritten December 2016 by M. Christi

   implicit none

!  Arguments

   integer  , intent(in)  :: order
   REAL(FPK), intent(in)  :: small, eps, y, delta, deldot, lamdot, kdot, wdel, udel, sm, &
                             lam, fac1, fac2, deldot2, kdot2
   REAL(FPK), intent(out) :: L_mult

!  Local

   integer   :: mterms1, mterms2, m, j

   real(fpk) :: c1,c2,c3,sigma,series1,series2,power
   real(fpk) :: a(0:max_Taylor_terms),b(0:max_Taylor_terms),&
                g(0:max_Taylor_terms),fact(0:max_Taylor_terms+1)

   real(fpk) :: dc1dx,dc2dx,dc3dx,depsdx,dsigdx
   real(fpk) :: dadx(0:max_Taylor_terms),dbdx(0:max_Taylor_terms),&
                dgdx(0:max_Taylor_terms),ABC(0:max_Taylor_terms),&
                dABCdx(0:max_Taylor_terms)

!  Zero it

   L_mult = zero

!  number of terms

   mterms2 = order + 2
   mterms1 = order + 1  ! Need this for same order of accuracy as Radiances


   fact   = zero
   sigma  = zero
   L_mult = zero
   depsdx = zero
   dsigdx = zero
   c1 = zero ; c2 = zero ; c3 = zero
   dc1dx = zero ; dc2dx = zero ; dc3dx = zero
   a = zero ; dadx = zero ; b = zero ; dbdx = zero ; g = zero ; dgdx = zero
   series1 = zero ; series2 = zero
   ABC = zero ; dABCdx = zero


   sigma = one/(lam - sm)
   c2 = fac2 !Wdel
   dc2dx = -c2*(lamdot*delta + lam*deldot2)

   fact(0) = one
   do m = 1, mterms2+1
      fact(m) = fact(m-1)*real(m,fpk)
   enddo

   if ( abs(sigma) .lt. small ) then
      !special case: triple small #

      a(0) = one ; dadx(0) = zero
      g(0) = one ; dgdx(0) = zero
      do m = 1, mterms2
         a(m) = delta**m/fact(m+1)
         dadx(m) = a(m-1)*real(m,fpk)*deldot2/real(m+1,fpk)

         g(m) = zero ; dgdx(m) = zero
         do j = 0, m
            g(m) = g(m) + a(j)*a(m-j)
            dgdx(m) = dgdx(m) + (dadx(j)*a(m-j) + a(j)*dadx(m-j))
         enddo
      enddo

      !for m=0
      series1 = g(0)
      !series2 = dgdx(0) + g(1)*depsdx
      series2 = g(1)*depsdx !=zero when eps=0? (I think yes).
      !for m=1,mterms1-1
      power = one
      do m = 1, mterms1
         power = power*eps
         series1 = series1 + g(m)*power !S
         if (m<mterms1) series2 = series2 + (dgdx(m) + g(m+1)*real(m+1,fpk)*depsdx)*power !dSdx
      enddo
      !for m=mterms1
      series2 = series2 + dgdx(mterms1)*power*eps

      L_mult = sm*delta*((dc2dx*delta + two*c2*deldot2)*series1 + c2*delta*series2)

   else

      !normal case: single small #
      depsdx = lamdot - kdot2
      dsigdx = -sigma**2*lamdot

      c1 = fac1 - fac2 !Udel - Wdel
      dc1dx = -fac1*sm*deldot2 - dc2dx

      a(0) = one   ; dadx(0) = zero
      a(1) = sigma ; dadx(1) = dsigdx 
      b(0) = delta ; dbdx(0) = deldot2
      g(0) = delta ; dgdx(0) = deldot2
      do m = 1, mterms2
         a(m) = sigma**m
         b(m) = delta**(m+1)/fact(m+1)
         dadx(m) = a(m-1)*dsigdx*real(m,fpk)
         dbdx(m) = b(m-1)*deldot2

         g(m) = zero ; dgdx(m) = zero
         do j = 0, m
            g(m) = g(m) + a(j)*b(m-j)
            dgdx(m) = dgdx(m) + (dadx(j)*b(m-j) + a(j)*dbdx(m-j))
         enddo
      enddo

      c3 = a(mterms1)*eps**mterms1
      dc3dx = dadx(mterms1)*eps**mterms1 + a(mterms1)*real(mterms1,fpk)*eps**(mterms1-1)*depsdx

      do m = 1, mterms2
         !ABC(m-1) = AA(m-1) - BB(m-1) - CC(m-1) for m=0,mterms1
         ABC(m-1) = c1*a(m) - c2*g(m-1) - c2*c3*b(m-1)
         !Derivative of ABC(m-1) = dAAdx(m-1) - dBBdx(m-1) - dCCdx(m-1) for m=0,mterms1 
         dABCdx(m-1) = (dc1dx*a(m) + c1*dadx(m)) - (dc2dx*g(m-1) + c2*dgdx(m-1)) &
                     - (dc2dx*c3 + c2*dc3dx)*b(m-1) + c2*c3*dbdx(m-1)
      enddo

      !for m=0
      power = one
      series1 = ABC(0)
      series2 = dABCdx(0) + ABC(1)*depsdx

      !for m=1,mterms1-1
      do m = 1, mterms1
         power = power*eps
         series1 = series1 + ABC(m)*power !S
         if (m<mterms1) series2 = series2 + (dABCdx(m) + ABC(m+1)*real(m+1,fpk)*depsdx)*power !dSdx
      enddo
      !for m=mterms1
      series2 = series2 + dABCdx(mterms1)*power*eps

      L_mult = sm*(dsigdx*series1 + sigma*series2)
   endif

!   Done

   return
end subroutine Taylor_Series_L_2b

!  Finish module

END MODULE lrrs_Taylor_m
