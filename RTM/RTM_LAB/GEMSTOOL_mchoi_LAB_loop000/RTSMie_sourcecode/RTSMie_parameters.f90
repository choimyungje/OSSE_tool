
!  Contains the following module
!       MODULE Mie_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE RTSMie_parameters_m

  IMPLICIT NONE

!  Dimensioning input for Mie code
!  ===============================

!  INTEGER, PARAMETER :: max_Mie_angles     = 200
!  INTEGER, PARAMETER :: max_Mie_sizes      = 20

!  This is an extreme case (5441). Use with Caution 
!    Former use of max_Mie_points caused some Segmentation errors
!    Needed lots of memory for arrays polyplus, polyminus
!    Interduced allocatable arrays inside Mie code to alleviate this.
!    RT SOLUTIONS Inc, R. Spurr, 25 March 2011

  INTEGER, PARAMETER :: max_Mie_angles     = 4000
  INTEGER, PARAMETER :: max_Mie_sizes      = 24

!  Define KIND variables for single and double precision

   INTEGER, PARAMETER :: sprec  = KIND( 1.0 )
   INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

! Numbers such as 1, 2, pie, etc.......

   REAL (KIND=dp), PARAMETER :: d_zero  = 0.0_dp
   REAL (KIND=dp), PARAMETER :: d_one   = 1.0_dp
   REAL (KIND=dp), PARAMETER :: d_two   = 2.0_dp
   REAL (KIND=dp), PARAMETER :: d_three = 3.0_dp
   REAL (KIND=dp), PARAMETER :: d_four  = 4.0_dp
   REAL (KIND=dp), PARAMETER :: d_half  = 0.5_dp

!  COMPLEX (KIND=dp), PARAMETER :: c_zero = ( d_zero, d_zero )
!  COMPLEX (KIND=dp), PARAMETER :: c_one =  ( d_one,  d_zero )
!  COMPLEX (KIND=dp), PARAMETER :: c_i =    ( d_zero,  d_one )

!  Mie version

   CHARACTER (LEN=10), PARAMETER :: Mie_f90_version = 'F90_Mie_V1'

END MODULE RTSMie_parameters_m
