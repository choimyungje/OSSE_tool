Module aerosol_loading_routines_m

!  Notes: 18 February 2013
!     profiles_lidar routine added
!      - Read LIDAR Extinction [km-1] and height [km] profiles from FILE
!      - Parcel the entire LIDAR profile into the output array
!      - Ignores z_upperlimit, z_lowerlimit
!      - Only one offline test so far (19 february, 2013)

!  Type definitions

   use vlidort_pars, only : fpk

private
public :: profiles_uniform,&
          profiles_expone,&
          exp1_profile_only,&
          exp1_profile_plus,&
          profiles_gdfone,&
          gdf1_profile_only,&
          gdf1_profile_plus, &
          profiles_lidar

contains

subroutine profiles_uniform                                 &
          ( maxlayers, nlayers, heightgrid, do_derivatives, & ! Inputs
            z_upperlimit, z_lowerlimit, total_column,       & ! Inputs
            profile, d_profile_dcol,                        & ! output
            fail, message, action )

      implicit none

!  Inputs
!  ======

!  Height grid in [km]

      integer, intent(in) :: maxlayers
      integer, intent(in) :: nlayers
      real(fpk),    intent(in) ::  heightgrid(0:maxlayers)

!  height levels controlling profile shape. All in [km].

!     z_upperlimit = uppermost height value at start  of Uniform slab
!     z_lowerlimit = lowest    height value at bottom of Uniform slab

!  Notes:
!    (1) Must have  z_u >/= z_l. This is checked internally
!    (2) Must have  z_l >/= heightgrid(nlayers). Also checked internally

      real(fpk),    intent(in) ::  z_upperlimit, z_lowerlimit

!  Total column between upper and lower limits

      real(fpk),    intent(in) ::  total_column

!  Control for calculating derivatives

      logical, intent(in) :: do_derivatives

!  Output
!  ======

!  Umkehr profile

      real(fpk),    intent(out) ::  profile ( maxlayers )

!  Umkehr profile derivatives  (w.r.t. total column )

      real(fpk),    intent(out) :: d_profile_dcol ( maxlayers )

!  Exception handling

      logical      , intent(out) ::  fail
      character*(*), intent(out) ::  message, action

!  Local variables
!  ===============

!  help variables

      integer          :: n, nupper, nlower
      real(fpk)        :: z1,z2,density, scaling, za,zb

!  Initial checks
!  ==============

!  Initialise exception handling

      fail = .false.
      message = ' '
      action  = ' '

!  Check input

      if ( z_lowerlimit .lt. heightgrid(nlayers) ) then
        fail = .true.
        message = 'Bad input: Uniform slab Lower limit is below ground'
        action  = ' increase z_lowerlimit > BOA, heightgrid(nlayers)'
        return
      endif

      if ( z_upperlimit .gt. heightgrid(0)  ) then
        fail = .true.
        message = 'Bad input: Uniform slab upperlimit is above TOA'
        action  = ' decrease z_upperlimit < TOA, heightgrid(0)'
        return
      endif

      if ( z_lowerlimit .gt. z_upperlimit ) then
        fail = .true.
        message = 'Bad input:  Uniform slab lower > upper limit'
        action  = ' increase z_upperlimit > z_Lowerlimit'
        return
      endif

!  Find intercept layers
!  ---------------------

      n = 0
      do while ( heightgrid(n).ge.z_upperlimit  )
        n = n + 1
      enddo
      nupper = n

      do while (heightgrid(n).gt.z_lowerlimit)
        n = n + 1
      enddo
      nlower = n

!  Zero output 

      do n = 1, nlayers
        profile(n) = 0.0d0
      enddo
      if ( do_derivatives ) then
        do n = 1, nlayers
          d_profile_dcol(n) = 0.0d0
        enddo
      endif

!  Basic settings

      z1 = z_upperlimit
      z2 = z_lowerlimit
      density = total_column / ( z1 - z2 )
 
!  construct profile Umkehrs
!  -------------------------

!  Layer containing z_upperlimit

      n = nupper
      za = z1
      if ( nlower .eq. nupper ) then
        scaling = 1.0d0
      else
        zb = heightgrid(nupper)
        scaling = ( z1 - zb ) / ( z1 - z2 )
      endif
      profile(n) = total_column * scaling
      if ( do_derivatives) d_profile_dcol(n) = scaling
 
!  Intermediate Layers

      do n = nupper+1,nlower-1
        za = heightgrid(n-1)
        zb = heightgrid(n)
        scaling = ( za - zb ) / ( z1 - z2 )
        profile(n)        = total_column * scaling
        if ( do_derivatives) d_profile_dcol(n) = scaling
      enddo

!  Layer containing z_lowerlimit
!   ( Will not be done if nupper = nlower)

      if ( nupper .lt. nlower ) then
        n = nlower
        za = heightgrid(n-1)
        zb = z2
        scaling = ( za - zb ) / ( z1 - z2 )
        profile(n)        = total_column * scaling
        if ( do_derivatives) d_profile_dcol(n) = scaling
      endif

!  Debug Check integrated profile
!  ------------------------------

!      sum = 0.0d0
!      sumd = 0.0d0
!      do n = nupper, nlower
!        sum = sum + profile(n)
!        if(do_derivatives)sumd = sumd + d_profile_dcol(n)
!      enddo
!      write(*,*)sum, total_column, sumd

!  Finish

      return
end subroutine profiles_uniform

!

subroutine profiles_expone                                         &
          ( maxlayers, nlayers, heightgrid, do_derivatives,        & ! Inputs
            z_upperlimit, z_lowerlimit, relaxation, total_column,  & ! Inputs
            profile, d_profile_dcol, d_profile_drxn,               & ! Output
            fail, message, action )

      implicit none

!  Inputs
!  ======

!  Height grid in [km]

      integer, intent(in) :: maxlayers
      integer, intent(in) :: nlayers
      real(fpk),    intent(in) ::  heightgrid(0:maxlayers)

!  height levels controlling profile shape. All in [km].

!     z_upperlimit = uppermost height value at start  of EXP profile
!     z_lowerlimit = lowest    height value at bottom of EXP profile

!  Notes:
!    (1) Must have  z_u >/= z_l. This is checked internally
!    (2) Must have  z_l >/= heightgrid(nlayers). Also checked internally

      real(fpk),    intent(in) ::  z_upperlimit, z_lowerlimit

!  Relaxation constant for the Exponential function

      real(fpk),    intent(in) ::  Relaxation

!  Total column between upper and lower limits

      real(fpk),    intent(in) ::  total_column

!  Control for calculating derivatives

      logical, intent(in) :: do_derivatives

!  Output
!  ======

!  Umkehr profile 

      real(fpk),    intent(out) ::  profile ( maxlayers )

!  Umkehr profile derivatives
!    (w.r.t. total column and relaxation parameter)

      real(fpk),    intent(out) :: d_profile_dcol ( maxlayers )
      real(fpk),    intent(out) :: d_profile_drxn ( maxlayers )

!  Exception handling

      logical      , intent(out) ::  fail
      character*(*), intent(out) ::  message, action

!  Local variables
!  ===============

!  help variables

      integer          :: n, nupper, nlower
      real(fpk)    :: k,b,ed,zd,kd,za,zb,term,omega,z1,z2,db_dk

!  Debug
!      logical          do_debug_profile
!      integer          m,msav,mdebug
!      parameter        ( mdebug = 1001 )
!      real(fpk)    pp(mdebug),zp(mdebug), sum, sumd, z, q

!  Initial checks
!  ==============

!  set debug output

!      do_debug_profile = .false.

!  Initialise exception handling

      fail = .false.
      message = ' '
      action  = ' '

!  Check input

      if ( z_lowerlimit .lt. heightgrid(nlayers) ) then
        fail = .true.
        message = 'Bad input: GDF Lower limit is below ground'
        action  = ' increase z_lowerlimit > heightgrid(nlayers)'
        return
      endif

      if ( z_lowerlimit .gt. z_upperlimit ) then
        fail = .true.
        message = 'Bad input: GDF transition height above upper limit'
        action  = ' increase z_upperlimit > z_lowerlimit'
        return
      endif

!  Find intercept layers
!  ---------------------

      n = 0
      do while ( heightgrid(n).ge.z_upperlimit  )
        n = n + 1
      enddo
      nupper = n

      do while (heightgrid(n).gt.z_lowerlimit)
        n = n + 1
      enddo
      nlower = n

!  Zero output

      do n = 1, nlayers
        profile(n) = 0.0d0
        if ( do_derivatives ) then
          d_profile_dcol(n) = 0.0d0
          d_profile_drxn(n) = 0.0d0
        endif
      enddo

!  set profile parameters

      k  = relaxation
      z1 = z_upperlimit
      z2 = z_lowerlimit
      omega = total_column
      zd = z1 - z2
      kd = k * zd
      ed = dexp ( - kd)
      term = kd - 1.0d0 + ed
      b  = omega * k / term

!  Set profile derivative parameters

      db_dk = ( omega - b * zd * (1.0 - ed ) ) / term
 
!  construct profile Umkehrs
!  -------------------------

!  Layer containing z_upperlimit. Partly EXP

      n = nupper
      za = z1
      zb = heightgrid(nupper)
      if ( do_derivatives ) then
        call exp1_profile_plus (za,zb,z1,b,k,omega,db_dk, &
              profile(n), d_profile_dcol(n), d_profile_drxn(n))
      else
        call exp1_profile_only(za,zb,z1,b,k,profile(n))
      endif

!  Intermediate Layers

      do n = nupper+1,nlower-1
        za = heightgrid(n-1)
        zb = heightgrid(n)
        if ( do_derivatives ) then
          call exp1_profile_plus (za,zb,z1,b,k,omega,db_dk, &
              profile(n), d_profile_dcol(n), d_profile_drxn(n))
        else
          call exp1_profile_only(za,zb,z1,b,k,profile(n))
        endif
      enddo

!  Layer containing z_lowerlimit
!   ( Will not be done if nupper = nlower)

      if ( nupper .lt. nlower ) then
        n = nlower
        za = heightgrid(n-1)
        zb = z2
        if ( do_derivatives ) then
          call exp1_profile_plus (za,zb,z1,b,k,omega,db_dk, &
             profile(n), d_profile_dcol(n), d_profile_drxn(n))
        else
          call exp1_profile_only(za,zb,z1,b,k,profile(n))
        endif
      endif

!  debug Check integrated profile
!  ------------------------------

!      sum = 0.0d0
!      sumd = 0.0d0
!      do n = nupper, nlower
!        sum = sum + profile(n)
!        if(do_derivatives)sumd = sumd + d_profile_dcol(n)
!      enddo

!  Debug profile itself
!  --------------------

!      if ( do_debug_profile ) then
!        m = 0
!        z = z1
!        do while (z.ge.z2)
!          z = z - 0.01
!          m = m + 1
!          if(m.gt.mdebug)stop'debug dimension needs to be increased'
!          zp(m) = z
!          q = dexp ( - k * (z1 - z ) )
!          pp(m) = b * (z1-z) + (b*q/k)
!        enddo
!        msav = m
!        do m = 1, msav
!          write(81,'(13f10.4)')zp(m),pp(m)
!        enddo
!      endif

!  Finish

      return
end subroutine profiles_expone

!

subroutine exp1_profile_only(za,zb,z1,b,k,prof)

!  Profile assignation

!  Arguments

      real(fpk),    intent(in)  :: za,zb,z1,b,k
      real(fpk),    intent(out) :: prof

!  Local variables

      real(fpk)    :: da,db,bok,ea,eb

!  Code

      prof = 0.0d0
      da = z1 - za
      db = z1 - zb
      bok = b / k
      ea = dexp(-k*da)
      eb = dexp(-k*db)
      prof = b * ( za - zb ) - bok * ( ea - eb )

!  Finish

      return
end subroutine exp1_profile_only

!

subroutine exp1_profile_plus ( za,zb,z1,b,k,omega,db_dk, &
                               prof,deriv1,deriv2)

!  Profile and derivative assignations

!  Arguments

      real(fpk),    intent(in)  :: za,zb,z1,b,k,omega,db_dk
      real(fpk),    intent(out) :: prof,deriv1,deriv2

!  Local variables

      real(fpk)    :: da,db,ea,eb,bok,term1,term2,term3

!  initialise

      prof   = 0.0d0
      deriv1 = 0.0d0
      deriv2 = 0.0d0

!  Code for profile

      da = z1 - za
      db = z1 - zb
      bok = b / k
      ea = dexp(-k*da)
      eb = dexp(-k*db)
      prof = b * ( za - zb ) - bok * ( ea - eb )

!  Derivatives

      deriv1 = prof / omega
      term1  = prof * db_dk / b
      term2  = (ea - eb ) * bok / k
      term3  = bok * ( da*ea - db*eb )
      deriv2 = term1 + term2 + term3

!  Finish

      return
end subroutine exp1_profile_plus

!

subroutine profiles_gdfone                                                     &
         ( maxlayers, nlayers, heightgrid, do_derivatives,                     &
           z_upperlimit, z_peakheight, z_lowerlimit, half_width, total_column, &
           profile, d_profile_dcol, d_profile_dpkh, d_profile_dhfw,            &
           fail, message, action )

      implicit none

!  Inputs
!  ======

!  Height grid in [km]

      integer, intent(in) ::  maxlayers
      integer, intent(in) ::  nlayers
      real(fpk),    intent(in) :: heightgrid(0:maxlayers)

!  height levels controlling profile shape. All in [km].

!     z_upperlimit = uppermost height value at start  of GDF profile
!     z_lowerlimit = lowest    height value at bottom of GDF profile
!     z_peakheight = height at which peak value of GDF profile occurs

!  Notes:
!    (1) Must have  z_u > z_p >/= z_l. This is checked internally
!    (2) Must have  z_l > heightgrid(nlayers). Also checked internally

      real(fpk),    intent(in) :: z_upperlimit, z_lowerlimit,  z_peakheight

!  Half width

      real(fpk),    intent(in) :: half_width

!  Total column between upper and lower limits

      real(fpk),    intent(in) :: total_column

!  Control for calculating derivatives

      logical, intent(in)   :: do_derivatives

!  Output
!  ======

!  Umkehr profile 

      real(fpk),    intent(out) :: profile ( maxlayers )

!  Umkehr profile derivatives
!    (w.r.t. total column and peak height and half-width)

      real(fpk),    intent(out)  :: d_profile_dcol ( maxlayers )
      real(fpk),    intent(out)  :: d_profile_dpkh ( maxlayers )
      real(fpk),    intent(out)  :: d_profile_dhfw ( maxlayers )

!  Exception handling

      logical      , intent(out)     :: fail
      character*(*), intent(out)     :: message, action

!  Local variables
!  ===============

!  help variables
     
      real(fpk)    :: omega,z1,z0,z2,h
      integer      :: n, nupper, nlower
      real(fpk)    :: a,zt1,zt2,za,zb
      real(fpk)    :: a_d,a_p,a_w,term,q1,q2,r1,r2,rsq1,rsq2

!  Debug

!      logical          do_debug_profile
!      integer          m,msav,mdebug
!      parameter        ( mdebug = 1001 )
!      real(fpk)    pp(mdebug),zp(mdebug), sum, sumd, z, q, r
     
!  Initial checks
!  ==============

!  set debug output

!      do_debug_profile = .true.

!  Initialise exception handling

      fail = .false.
      message = ' '
      action  = ' '

!  Check input

      if ( z_lowerlimit .lt. heightgrid(nlayers) ) then
        fail = .true.
        message = 'Bad input: GDF Lower limit is below ground'
        action  = ' increase z_lowerlimit > heightgrid(nlayers)'
        return
      endif

      if ( z_lowerlimit .gt. z_peakheight ) then
        fail = .true.
        message = 'Bad input: GDF Lower limit is above peak height'
        action  = ' increase z_peakheight > z_lowerlimit'
        return
      endif

      if ( z_peakheight .gt. z_upperlimit ) then
        fail = .true.
        message = 'Bad input: GDF transition height above upper limit'
        action  = ' increase z_upperlimit > z_peakheight'
        return
      endif

!  Find intercept layers
!  ---------------------

      n = 0
      do while ( heightgrid(n).ge.z_upperlimit  )
        n = n + 1
      enddo
      nupper = n

      do while (heightgrid(n).gt.z_lowerlimit)
        n = n + 1
      enddo
      nlower = n

!  Zero output 

      profile = 0.0d0
      if ( do_derivatives ) then
        d_profile_dcol = 0.0d0
        d_profile_dpkh = 0.0d0
        d_profile_dhfw = 0.0d0
      endif

!  Basic settings

      h  = half_width
      z1 = z_upperlimit
      z0 = z_peakheight
      z2 = z_lowerlimit
      omega = total_column

!  set profile parameters

      zt1 = z1 - z0
      zt2 = z0 - z2
      q1 = dexp ( - h * zt1)
      q2 = dexp ( - h * zt2)
      r1 = 1.0d0 / ( 1.0d0 + q1 )
      r2 = 1.0d0 / ( 1.0d0 + q2 )
      term = ( r1 + r2 - 1.0d0 )
      a  = omega * h / term

!  Set profile derivative parameters

      if ( do_derivatives ) then
        rsq1 = r1 * r1 * q1
        rsq2 = r2 * r2 * q2
        a_d = a / omega
        a_p  = - a * h * ( rsq2 - rsq1 ) / term
        a_w = ( omega - a * ( rsq2*zt2 + rsq1*zt1 ) ) / term
       endif

!  construct profile Umkehrs
!  -------------------------

!  Layer containing z_upperlimit. Partly EXP

      n = nupper
      za = z1
      zb = heightgrid(nupper)
      if ( do_derivatives ) then
        call gdf1_profile_plus (za,zb,z0,h,a,a_d,a_p,a_w,                           & ! Inputs
              profile(n), d_profile_dcol(n),  d_profile_dpkh(n),d_profile_dhfw(n))    ! Output
      else
        call gdf1_profile_only(za,zb,z0,h,a,profile(n))
      endif

!  Intermediate Layers

      do n = nupper+1,nlower-1
        za = heightgrid(n-1)
        zb = heightgrid(n)
        if ( do_derivatives ) then
          call gdf1_profile_plus (za,zb,z0,h,a,a_d,a_p,a_w,                         & ! Inputs
              profile(n), d_profile_dcol(n),  d_profile_dpkh(n),d_profile_dhfw(n))    ! Output
        else
          call gdf1_profile_only(za,zb,z0,h,a,profile(n))
        endif
      enddo

!  Layer containing z_lowerlimit
!   ( Will not be done if nupper = nlower)

      if ( nupper .lt. nlower ) then
        n = nlower
        za = heightgrid(n-1)
        zb = z2
        if ( do_derivatives ) then
          call gdf1_profile_plus (za,zb,z0,h,a,a_d,a_p,a_w,                         & ! Inputs
              profile(n), d_profile_dcol(n),  d_profile_dpkh(n),d_profile_dhfw(n))    ! Output
        else
          call gdf1_profile_only(za,zb,z0,h,a,profile(n))
        endif
      endif

!  Debug Check integrated profile
!  ------------------------------

!      sum = 0.0d0
!      sumd = 0.0d0
!      do n = nupper, nlower
!        sum = sum + profile(n)
!        if(do_derivatives)sumd = sumd + d_profile_dcol(n)
!      enddo
!      write(*,*)sum, omega, sumd

!  Debug profile itself
!  --------------------

!      if ( do_debug_profile ) then
!        m = 0
!        z = z1
!        do while (z.ge.z2)
!         z = z - 0.01
!         m = m + 1
!         if(m.gt.mdebug)stop'debug dimension needs to be increased'
!          zp(m) = z
!          q = dexp ( - h * dabs(( z - z0 )) )
!          r = 1.0d0 + q
!          pp(m) = a * q / r / r
!        enddo
!        msav = m
!        do m = 1, msav
!          write(81,'(13f10.4)')zp(m),pp(m)
!        enddo
!      endif

!  Finish

      return
end subroutine profiles_gdfone

!

subroutine gdf1_profile_only(za,zb,z0,h,a,prof)

!  Profile assignation

!  Arguments

      real(fpk),    intent(in)   :: za,zb,z0,h,a
      real(fpk),    intent(out)  :: prof

!  Local variables

      real(fpk)     :: fac,ra,rb

!  initialise

      prof = 0.0d0

      fac = a / h
      ra = 1.0d0 / ( 1.0d0 + dexp ( - h * dabs(za-z0)) ) 
      rb = 1.0d0 / ( 1.0d0 + dexp ( - h * dabs(zb-z0)) )
      if ( za.gt.z0 ) then
        if ( zb .ge. z0 ) then
          prof = fac * ( ra - rb )
        else 
          prof = fac * ( ra + rb - 1.0d0 )
        endif
      else
        prof = fac * ( rb - ra )
      endif

!  Finish

      return
end subroutine gdf1_profile_only

!

subroutine gdf1_profile_plus ( za,zb,z0,h,a,a_d,a_p,a_w,   & ! Input
                              prof,deriv1,deriv2,deriv3)     ! Output

!  Profile and derivative assignation

!  Arguments

      real(fpk),    intent(in)   :: za,zb,z0,h,a,a_d,a_p,a_w
      real(fpk),    intent(out)  :: prof,deriv1,deriv2,deriv3

!  Local variables

      real(fpk)    :: fac,qa,qb,ra,rb,mza,mzb,rsqa,rsqb,rab,fac_d
      real(fpk)    :: fac_p,ra_p,rb_p,rab_p
      real(fpk)    :: fac_w,ra_w,rb_w,rab_w

!  initialise

      prof   = 0.0d0
      deriv1 = 0.0d0
      deriv2 = 0.0d0
      deriv3 = 0.0d0

!  Type 1

      fac = a / h
      mza =  dabs(za-z0)
      mzb =  dabs(zb-z0)
      qa = dexp ( - h * mza)
      qb = dexp ( - h * mzb)
      ra = 1.0d0 / ( 1.0d0 + qa ) 
      rb = 1.0d0 / ( 1.0d0 + qb )
      fac_d = a_d/h
      fac_p = a_p/h
      fac_w = (a_w/h) - (fac/h)
      rsqa =  qa * ra * ra
      rsqb =  qb * rb * rb
      ra_p = rsqa * h
      rb_p = rsqb * h
      ra_w = rsqa * mza
      rb_w = rsqb * mzb
      if ( za.gt.z0 ) then
        if ( zb .ge. z0 ) then
          rab   =   ra   - rb
          rab_p = - ra_p + rb_p
          rab_w =   ra_w - rb_w
        else 
          rab   =   ra   + rb   - 1.0d0
          rab_p = - ra_p + rb_p
          rab_w =   ra_w + rb_w
        endif
      else
        rab    = rb   - ra
        rab_p  = rb_p - ra_p
        rab_w  = rb_w - ra_w
      endif
      prof   = fac   * rab
      deriv1 = fac_d * rab
      deriv2 = fac_p * rab + fac * rab_p
      deriv3 = fac_w * rab + fac * rab_w

!  Finish

      return
end subroutine gdf1_profile_plus

!

subroutine profiles_lidar &
      ( maxlayers, nlayers, filename, heightgrid,  & ! Inputs
        profile, d_profile_dtau, d_profile_dpars,  & ! Outputs
        fail, message, action )                      ! Exception Handling

   implicit none
   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Height grid in [km]

   integer, intent(in)   ::  maxlayers
   integer, intent(in)   ::  nlayers
   real(fpk), intent(in) :: heightgrid(0:maxlayers)

!  Filename fo Lidar sources

   character*(*), intent(in) :: filename

!  Output
!  ======

!  profile 

   real(fpk),    intent(out)   :: profile ( maxlayers )

!  profile derivatives. Dummies here

   real(fpk),    intent(out)   :: d_profile_dtau  ( maxlayers )
   real(fpk),    intent(out)   :: d_profile_dpars ( maxlayers, 2 )

!  Exception handling

   logical      , intent(out)  :: fail
   character*(*), intent(out)  :: message, action

!  Local variables
!  ===============

!  Lidar data

   integer, parameter :: max_lidar_levels = 10000
   integer    :: n_lidar_levels
   real(fpk)  :: lidar_levels  ( 0:max_lidar_levels )
   real(fpk)  :: lidar_extincs ( 0:max_lidar_levels )
   real(fpk)  :: lidar_taus    ( max_lidar_levels )

!  Other variables

   integer   :: n, n1, nz, nupper, nlower, nmax, nstart
   logical   :: loop
   real(fpk) :: maxl, ext, dif, tau, fac, tau_div1, tau_div2, x, x1, y, y1

!  Zero output

   profile         = 0.0_fpk
   d_profile_dtau  = 0.0_fpk
   d_profile_dpars = 0.0_fpk

   fail = .false.
   message = ' ' ; action = ' '

!   n_lidar_levels = 9001 ; write(14,'(I5)')n_lidar_levels
!   do n = 1, n_lidar_levels
!      write(14,'(1p2e17.7)')10.0777d0-0.001*real(n),0.00004d0+0.0001d0*real(n)
!   enddo
!   pause

!  Read file

   open(1,file=trim(filename),status='old')
   read(1,*)n_lidar_levels ; maxl = -999.99_fpk
   do n = 1, n_lidar_levels
      n1 = n - 1
      read(1,*)lidar_levels(n1),lidar_extincs(n1)
      maxl = max(maxl,lidar_levels(n1)) ; if (maxl.eq.lidar_levels(n1))nmax=n1
   enddo
   close(1)

!  reverse order if necessary

   n_lidar_levels = n_lidar_levels - 1
   if ( nmax.ne.0) then
      do n = 0, n_lidar_levels
         n1 = n_lidar_levels - n
         x  = lidar_levels(n)  ; y  = lidar_extincs(n)
         x1 = lidar_levels(n1) ; y1 = lidar_extincs(n1)
         lidar_levels(n)  = x1 ; lidar_extincs(n)  = y1
         lidar_levels(n1) = x  ; lidar_extincs(n1) = y
      enddo
   endif

!  develop tau

   do n = 0, n_lidar_levels-1
      n1 = n + 1
      ext = 0.5_fpk * ( lidar_extincs(n) + lidar_extincs(n1) )
      dif = lidar_levels(n) - lidar_levels(n1)
      lidar_taus(n1) = ext * dif
   enddo

!  Find intercept layers
!  ---------------------

   n = 0
   do while ( heightgrid(n).ge.lidar_levels(0)  )
      n = n + 1
   enddo
   nupper = n

   do while (heightgrid(n).ge.lidar_levels(n_lidar_levels))
      n = n + 1
   enddo
   nlower = n

!   write(*,*)nupper,nlower ; pause'2'

!  Parceling for each layer
!  ------------------------

!  Top layer

   nz = nupper ; nstart = 0
   n = nstart ; loop = .true. ; tau = 0.0_fpk
   do while ( loop .and.(n.ge.nstart) )
      n = n + 1
      loop = (lidar_levels(n).gt.heightgrid(nz))
      if ( loop ) tau = tau + lidar_taus(n)
   enddo
   n1 = n - 1 ; dif = lidar_levels(n1) - lidar_levels(n)
   fac = ( lidar_levels(n1) - heightgrid(nz)  ) / dif
   tau_div1 = fac*lidar_taus(n) ; tau_div2 = (1.0d0-fac)*lidar_taus(n)
   profile(nz) = tau + tau_div1

!  Intermediate

   do nz = nupper+1, nlower -1
      nstart = n ; n = nstart ; loop = .true. ; tau = tau_div2
      do while ( loop .and.(n.ge.nstart) )
         n = n + 1
         loop = (lidar_levels(n).gt.heightgrid(nz))
         if ( loop ) tau = tau + lidar_taus(n)
      enddo
      n1 = n - 1 ; dif = lidar_levels(n1) - lidar_levels(n)
      fac = ( lidar_levels(n1) - heightgrid(nz)  ) / dif
      tau_div1 = fac*lidar_taus(n) ; tau_div2 = (1.0d0-fac)*lidar_taus(n)
      profile(nz) = tau + tau_div1
   enddo

!  Lowest

   nstart = n ; n = nstart ; loop = .true. ; tau = tau_div2
   do while ( loop .and.(n.ge.nstart) )
      n = n + 1
      loop = (lidar_levels(n).gt.heightgrid(nz))
      if ( loop ) tau = tau + lidar_taus(n)
   enddo
   n1 = n - 1 ; dif = lidar_levels(n1) - lidar_levels(n)
   fac = ( lidar_levels(n1) - heightgrid(nz)  ) / dif
   tau_div1 = fac*lidar_taus(n) ; tau_div2 = (1.0d0-fac)*lidar_taus(n)
   profile(nz) = tau + tau_div1

!  Check: these two numbers should be the same

!   write(*,*)'total LIDAR    tau =  ',sum(lidar_taus(1:n_lidar_levels))
!   write(*,*)'total Parceled tau =  ',sum(profile(1:nlayers))

!  Finish

   return
end subroutine profiles_lidar

!  End module
 
End Module aerosol_loading_routines_m

