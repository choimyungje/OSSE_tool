MODULE Rayleigh_function_m

!  Rayleigh Cross-sections and depolarizaiton ratios
!   Given separate module, 21 October 2013

!  stand-alone routine

public

contains

SUBROUTINE RAYLEIGH_FUNCTION                       &
          ( FORWARD_MAXLAMBDAS, CO2_PPMV_MIXRATIO, &
            FORWARD_NLAMBDAS,   FORWARD_LAMBDAS,   &
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL )

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Rayleigh cross sections and depolarization ratios
!     Bodhaine et. al. (1999) formulae
!     Module is stand-alone.
!     Wavelengths in nm (formerly Angstroms)

!  Input arguments
!  ---------------

!  wavelength
 
      INTEGER     :: FORWARD_MAXLAMBDAS, FORWARD_NLAMBDAS
      real(fpk),    dimension ( FORWARD_MAXLAMBDAS ) :: FORWARD_LAMBDAS 

!  CO2 mixing ratio

      real(fpk)    :: CO2_PPMV_MIXRATIO

!  Output arguments
!  ----------------

!  cross-sections and depolarization output

      real(fpk),    dimension ( FORWARD_MAXLAMBDAS ) :: RAYLEIGH_XSEC 
      real(fpk),    dimension ( FORWARD_MAXLAMBDAS ) :: RAYLEIGH_DEPOL

!  Local variables
!  ---------------

      INTEGER      :: W
      real(fpk)    :: MASS_DRYAIR
      real(fpk)    :: NMOL, PI, CONS
      real(fpk)    :: MO2,MN2,MARG,MCO2,MAIR
      real(fpk)    :: FO2,FN2,FARG,FCO2,FAIR
      real(fpk)    :: LAMBDA_C,LAMBDA_M,LPM2,LP2
      real(fpk)    :: N300M1,NCO2M1,NCO2
      real(fpk)    :: NCO2SQ, NSQM1,NSQP2,TERM

!  data statements and parameters
!  ------------------------------

      DATA MO2  / 20.946D0 /
      DATA MN2  / 78.084D0 /
      DATA MARG / 0.934D0 /

      real(fpk),    PARAMETER ::        S0_A = 15.0556D0
      real(fpk),    PARAMETER ::        S0_B = 28.9595D0

      real(fpk),    PARAMETER ::        S1_A = 8060.51D0
      real(fpk),    PARAMETER ::        S1_B = 2.48099D+06
      real(fpk),    PARAMETER ::        S1_C = 132.274D0
      real(fpk),    PARAMETER ::        S1_D = 1.74557D+04
      real(fpk),    PARAMETER ::        S1_E = 39.32957D0

      real(fpk),    PARAMETER ::        S2_A = 0.54D0

      real(fpk),    PARAMETER ::        S3_A = 1.034D0
      real(fpk),    PARAMETER ::        S3_B = 3.17D-04
      real(fpk),    PARAMETER ::        S3_C = 1.096D0
      real(fpk),    PARAMETER ::        S3_D = 1.385D-03
      real(fpk),    PARAMETER ::        S3_E = 1.448D-04

!  Start of code
!  -------------

!  constants

      NMOL = 2.546899D19
      PI   = DATAN(1.0D0)*4.0D0
      CONS = 24.0D0 * PI * PI * PI

!  convert co2

      MCO2 = 1.0D-06 * CO2_PPMV_MIXRATIO

!  mass of dry air: Eq.(17) of BWDS

      MASS_DRYAIR = S0_A * MCO2 + S0_B

!  start loop

      DO W = 1, FORWARD_NLAMBDAS

!  wavelength in micrometers

      LAMBDA_M = 1.0D-03 * FORWARD_LAMBDAS(W)
      LAMBDA_C = 1.0D-07 * FORWARD_LAMBDAS(W)
!      LAMBDA_M = 1.0D-04 * FORWARD_LAMBDAS(W)             ! Angstrom input
!      LAMBDA_C = 1.0D-08 * FORWARD_LAMBDAS(W)             ! Angstrom input
      LPM2     = 1.0D0 / LAMBDA_M / LAMBDA_M

!  step 1: Eq.(18) of BWDS

      N300M1 = S1_A + ( S1_B / ( S1_C - LPM2 ) ) + &
                      ( S1_D / ( S1_E - LPM2 ) )
      N300M1 = N300M1 * 1.0D-08

!  step 2: Eq.(19) of BWDS

      NCO2M1 = N300M1 * ( 1.0D0 + S2_A * ( MCO2  - 0.0003D0 ) )
      NCO2   = NCO2M1 + 1
      NCO2SQ = NCO2 * NCO2

!  step 3: Eqs. (5&6) of BWDS (Bates' results)

      FN2  = S3_A + S3_B * LPM2
      FO2  = S3_C + S3_D * LPM2 + S3_E * LPM2 * LPM2

!  step 4: Eq.(23) of BWDS
!     ---> King factor and depolarization ratio

      FARG = 1.0D0
      FCO2 = 1.15D0
      MAIR = MN2 + MO2 + MARG + MCO2
      FAIR = MN2*FN2 + MO2*FO2 + MARG*FARG + MCO2*FCO2
      FAIR = FAIR / MAIR
      RAYLEIGH_DEPOL(W) = 6.0D0*(FAIR-1.0D0)/(3.0D0+7.0D0*FAIR)

!  step 5: Eq.(22) of BWDS
!     ---> Cross section

      LP2  = LAMBDA_C * LAMBDA_C
      NSQM1 = NCO2SQ - 1.0D0
      NSQP2 = NCO2SQ + 2.0D0
      TERM = NSQM1 / LP2 / NMOL / NSQP2
      RAYLEIGH_XSEC(W) =  CONS * TERM * TERM * FAIR

!  end loop

      ENDDO

!  finish
!  ------

      RETURN
END SUBROUTINE RAYLEIGH_FUNCTION

!  End module

END MODULE Rayleigh_function_m

