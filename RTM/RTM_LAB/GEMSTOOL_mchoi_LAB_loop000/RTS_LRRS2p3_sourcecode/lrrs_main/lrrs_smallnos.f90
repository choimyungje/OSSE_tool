! ###############################################################
! #                                                             #
! #                    THE LIDORT_RRS MODEL                     #
! #                                                             #
! #      (LInearized Discrete Ordinate Radiative Transfer)      #
! #       --         -        -        -         -              #
! #                 (Rotational Raman Scatter)                  #
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
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Version      :  2.3                                        #
! #  Release Date :  March 2011                                 #
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
! #             limit_gcfunc                                    #
! #             limit_gcfunc_ut                                 #
! #             limit_ppgmult_wu                                #
! #             limit_ppgmult_wd                                #
! #             limit_ppgmult_pu                                #
! #             limit_ppgmult_pd                                #
! #                                                             #
! ###############################################################

!  @@@ RobFix 5/5/11.
!    LIMIT_PPGMULT_WD, LIMIT_PPGMULT_PD
!       double small-number treatment to zeroth order

      MODULE lrrs_smallnos

      USE LRRS_PARS

      PRIVATE
      PUBLIC :: LIMIT_GCFUNC,&
                LIMIT_GCFUNC_UT,&
                LIMIT_PPGMULT_WD,&
                LIMIT_PPGMULT_WU,&
                LIMIT_PPGMULT_PD,&
                LIMIT_PPGMULT_PU

      CONTAINS

      SUBROUTINE LIMIT_GCFUNC &
          ( EPS, DELTA, ZDEL, CFUNC )

!  Green's function solution multiplier CFUNC (whole layer)
!    Small number expansion to second order

      IMPLICIT NONE

!  input arguments

      REAL(FPK), INTENT(IN) :: eps, delta, zdel

!  output arguments

      REAL(FPK), INTENT(OUT) :: cfunc

!  local declarations

      REAL(FPK) :: power, power2

!  initialise

      cfunc = 0.0d0

!  evaluate

      power = delta * eps
      power2 = power * power
      cfunc = zdel*delta*(1.0d0-0.5d0*power+power2/6.0d0)

!  Finish

      return
      END SUBROUTINE LIMIT_GCFUNC

!  ####################################################################

      SUBROUTINE LIMIT_GCFUNC_UT &
          ( EPS, TAUX, ZXDN, CFUNC_UT )

!  Green's function solution multiplier CFUNC_UT (partial layer)
!    Small number expansion to second order

      IMPLICIT NONE

!  input arguments

      REAL(FPK), INTENT(IN) :: eps, taux, zxdn

!  output arguments

      REAL(FPK), INTENT(OUT) :: cfunc_ut

!  local declarations

      REAL(FPK) :: power, power2

!  initialise

      cfunc_ut = 0.0d0

!  evaluate

      power = taux * eps
      power2 = power * power
      cfunc_ut = zxdn*taux*(1.0d0-0.5d0*power+power2/6.0d0)

!  Finish

      return
      END SUBROUTINE LIMIT_GCFUNC_UT

!  ####################################################################

      SUBROUTINE LIMIT_PPGMULT_WD &
          ( EPS, SM, CONS, KALPHA, DELTA, ZDEL, UDEL, GMULT)

!  Post-processing Green's function multiplier, whole layer downwelling
!    Small number expansion to second order

      IMPLICIT NONE

!  input arguments

      REAL(FPK), INTENT(IN) :: eps, sm, cons, kalpha, delta, zdel, udel

!  output arguments

      REAL(FPK), INTENT(OUT) :: gmult

!  local declarations

      REAL(FPK) :: delta2, rm, smrm, h1, g1,term0,term1,term2

!  initialise

      gmult = 0.0d0

!  @@@ RobFix 5/5/11. Double small-number clause.
      if ( dabs(sm - kalpha) .lt. 1.0d-03 ) then
        delta2 = delta*delta
        gmult = 0.5d0 * udel * delta2
        return
      endif
!  @@@ End RobFix 5/5/11.

!  evaluate

      delta2 = delta*delta
      rm   = 1.0d0 / ( sm - kalpha )
      smrm = sm * rm
      h1 = smrm * ( zdel - udel)
      g1 = smrm * zdel * delta
      term0 =  g1 - h1*rm
      term1 = rm*term0 - g1*delta*0.5d0
      term2 = rm*term1 + g1*delta2/6.0d0
      gmult = cons * ( term0 + eps*term1 + eps*eps*term2)

!  Finish

      return
      END SUBROUTINE LIMIT_PPGMULT_WD

!  ####################################################################

      SUBROUTINE LIMIT_PPGMULT_WU &
          ( EPS, SM, CONS, KALPHA, DELTA, ZDEL, UDEL, GMULT)

!  Post-processing Green's function multiplier, whole layer upwelling
!    Small number expansion to second order

      IMPLICIT NONE

!  input arguments

      REAL(FPK), INTENT(IN) :: eps, sm, cons, kalpha, delta, zdel, udel

!  output arguments

      REAL(FPK), INTENT(OUT) :: gmult

!  local declarations

      REAL(FPK) :: delta2,rp,smrp,h2,g2,zu,term0,term1,term2

!  initialise

      gmult = 0.0d0

!  evaluate

      delta2 = delta*delta
      rp   = 1.0d0 / ( sm + kalpha )
      smrp = sm * rp
      zu = zdel * udel
      h2 = smrp * ( 1.0d0 - zu)
      g2 = smrp * zu * delta
      term0 =  g2 - h2*rp
      term1 = -rp*term0 - g2*delta*0.5d0
      term2 = -rp*term1 + g2*delta2/6.0d0
      gmult = - cons * ( term0 + eps*term1 + eps*eps*term2)

!  Finish

      return
      END SUBROUTINE LIMIT_PPGMULT_WU

!  ####################################################################

      SUBROUTINE LIMIT_PPGMULT_PD &
          ( EPS, SM, CONS, KALPHA, TAUX, ZXDN, UXDN, GMULT)

!  Post-processing Green's function multiplier, part-layer downwelling
!    Small number expansion to second order

      IMPLICIT NONE

!  input arguments

      REAL(FPK), INTENT(IN) :: eps, sm, cons, kalpha, taux, zxdn, uxdn

!  output arguments

      REAL(FPK), INTENT(OUT) :: gmult

!  local declarations

      REAL(FPK) :: taux2, rm, smrm, h1, g1,term0,term1,term2

!  initialise

      gmult = 0.0d0

!  evaluate
!  @@@ RobFix 5/5/11. Double small number clause.

      taux2 = taux*taux
      if ( dabs(sm - kalpha) .lt. 1.0d-03 ) then
        gmult = 0.5d0 * uxdn * taux2
      else
        rm   = 1.0d0 / ( sm - kalpha )
        smrm = sm * rm
        h1 = smrm * ( zxdn - uxdn)
        g1 = smrm * zxdn * taux
        term0 =  g1 - h1*rm
        term1 = rm*term0 - g1*taux*0.5d0
        term2 = rm*term1 + g1*taux2/6.0d0
        gmult = cons * ( term0 + eps*term1 + eps*eps*term2)
      endif

!  Finish

      return
      END SUBROUTINE LIMIT_PPGMULT_PD

!  ####################################################################

      SUBROUTINE LIMIT_PPGMULT_PU &
        ( EPS, SM, CONS, KALPHA, DELTA, TAUX, ZDEL, ZXDN, UXUP, GMULT)

!  Post-processing Green's function multiplier, part-layer upwelling
!    Small number expansion to second order

      IMPLICIT NONE

!  input arguments

      REAL(FPK), INTENT(IN) :: eps, sm, cons, kalpha, delta, taux
      REAL(FPK), INTENT(IN) :: zdel, zxdn, uxup

!  output arguments

      REAL(FPK), INTENT(OUT) :: gmult

!  local declarations

      REAL(FPK) :: delta2,delta3,taux2,taux3,rp,rp2,rp3,smrp
      REAL(FPK) :: zu,c0,c1,c2,c3,term0,term1,term2

!  initialise

      gmult = 0.0d0

!  evaluate

      delta2 = delta*delta
      delta3 = delta2*delta
      taux2 = taux  * taux
      taux3 = taux2 * taux
      rp   = 1.0d0 / ( sm + kalpha )
      rp2 = rp * rp
      rp3 = rp2 * rp
      smrp = sm * rp
      zu = zdel * uxup
      c0 = zxdn - zu
      c1 =   - taux  * zxdn + zu * delta
      c2 = (   taux2 * zxdn - zu * delta2 ) * 0.5d0
      c3 = ( - taux3 * zxdn + zu * delta3 ) / 6.0d0
      term0 = c1 - c0*rp
      term1 = c2 - c1*rp + c0*rp2
      term2 = c3 - c2*rp + c1*rp2 - c0*rp3
      gmult = - cons * smrp * ( term0 + eps*term1 + eps*eps*term2)

!  Finish

      return
      END SUBROUTINE LIMIT_PPGMULT_PU


      END MODULE lrrs_smallnos

