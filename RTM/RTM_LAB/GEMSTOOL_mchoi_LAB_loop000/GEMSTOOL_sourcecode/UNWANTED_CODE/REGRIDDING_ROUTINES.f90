Module Regridding_Routines_m

!  Type definitions

   use vlidort_pars, only : fpk

!  All routines are Public

private
public :: PTH_REGRIDDING_NOAER,&
          PTH_REGRIDDING_GENERAL,&
          PTH_REGRIDDING_ALT_NOAER,&
          PTH_REGRIDDING_ALT_GENERAL

contains

!  PURE INPUTS

!maxlayers, maxgases (dimensioning)

!ngases = 1
!do_gases(1) = .true.
!cloud_method   = 2 (mie scattering cloud)
!do_clearsky, Flag for doing a clear sky calculation, with (Method = 1/2) or without (Method = 0) regridding
!cloud_opdep.  Cloud optical thickness at one wavelength, GOING to be PARCELLED (in cparcel)

!surface_press
!ctop_pressures
!cbot_pressure

!  MODIFIED INPUTS

!nlayers_total. Originally N layers, come out (usually) with N+2
!press   = Level pressures
!heights = level heights
!temps   = LAYER temperatures
!shapegrid = LAYER Temperature-shift functions, set these to 0.0 throughout (not useful)

!profiles(gas =1) = ozone profile
!profiles_dc = derivative of the profile w.r.t its total column. Set this to 0.0 (not needed)

!airdens = air column desities, each layer
!ds_airdens = derivative of airdens w.r.t. T-shift. Set to 0.0d0 (Not required)

!  PURE OUTPUTS

!nlayers_cloud. Only for Method 1
!cloudflag(n) = (Method 2) .true. if n is the cloud layer. SAME as the aerosol FLAG.
!cparcel(n)   = (Method 2) loading. Same as AEROSOL_LOADING (layer optical thickness values)
!                        in your master routine. 


      subroutine PTH_REGRIDDING_NOAER                       &
       ( maxlayers, maxgases, ngases, do_gases,             &
         cloud_method, do_clearsky, cloud_opdep,            &
         surface_press, ctop_pressure, cbot_pressure,       &
         nlayers_total, press, heights, temps, shapegrid,   &      ! MODIFIED INPUTS
         profiles, dc_profiles, airdens, ds_airdens,        &      ! MODIFIED INPUTS
         nlayers_cloud, cloudflag, cparcel )                       !  NEW OUTPUTS

      implicit none

!  dimensioning

      integer          maxlayers, maxgases

!  gas control

      integer          ngases
      logical          do_gases ( maxgases)

!  control inputs (index, surface and cloud information)

      logical          do_clearsky
      integer          cloud_method
      real(fpk)        surface_press
      real(fpk)        cloud_opdep
      real(fpk)        ctop_pressure,cbot_pressure

!  Modified Grid outputs

      integer          nlayers_total
      real(fpk)        press     (0:maxlayers)
      real(fpk)        temps     (  maxlayers)
      real(fpk)        heights   (0:maxlayers)
      real(fpk)        shapegrid (  maxlayers)

      real(fpk)        profiles    (maxlayers,maxgases)
      real(fpk)        dc_profiles (maxlayers,maxgases)
      real(fpk)        airdens     (maxlayers)
      real(fpk)        ds_airdens  (maxlayers)

!  Completely new outputs for cloud information

      integer          nlayers_cloud
      logical          cloudflag (maxlayers)
      real(fpk)        cparcel   (maxlayers)

!  status

      character*70    message
      logical         do_fail

!  local variables
!  ---------------

!  holding arrays 

      integer          maxlayers_hold, maxgases_hold, nlayers_hold
      parameter        ( maxlayers_hold = 101 )
      parameter        ( maxgases_hold  = 2  )

      real(fpk)        hold_press  (0:maxlayers_hold)
      real(fpk)        hold_temps  (  maxlayers_hold)
      real(fpk)        hold_shape  (  maxlayers_hold)

      real(fpk)        hold_heights(0:maxlayers_hold)
      real(fpk)        hold_logp   (0:maxlayers_hold)
      real(fpk)        hold_airdens(  maxlayers_hold)
      real(fpk)        hold_dsair  (  maxlayers_hold)

      real(fpk)        hold_pr      (maxlayers_hold,maxgases_hold)
      real(fpk)        hold_prderiv (maxlayers_hold,maxgases_hold)

      real(fpk)        pr_t  (maxgases_hold), ad_t
      real(fpk)        prd_t (maxgases_hold), add_t
      real(fpk)        pr_b  (maxgases_hold), ad_b
      real(fpk)        prd_b (maxgases_hold), add_b

!  other help variables

      integer          n, nc, n1,n2, nc1, nc2, g
      integer          nlayers_s, nlayers_t, nlayers_b
      logical          loop
      real(fpk)        gradh, dlp, dlpb, dlpt, diffh, fac_t, fac_b
      real(fpk)        hgt_b, hgt_t, tmp_b, tmp_t, shp_t, shp_b
      real(fpk)        log_ptop, log_pbot

!  Initialise (important!)
!  ----------

      do_fail    = .false.
      message    = ' '

!  initialise nlayers_cloud and cloudflag and cparcel
!       (the only new variables to be set)

      nlayers_cloud = 0
      do n = 1, maxlayers
        cloudflag(n) = .false.
        cparcel(n)   = 0.0d0
      enddo

!  Nothing to do if Clearsky, Cloud_method = 0

      if ( cloud_method .eq. 0 ) return

!  Checking input
!  --------------

!  Check the cloud-top pressure is not out of bounds.

      if ( ctop_pressure .gt. surface_press ) then
        do_fail = .true.
        message = 'Abort: Ctop pressure > surface pressure'
        write(*,*)message
        stop
      endif

!  Check the cloud bottom pressure (clouds as layers only)

      if ( cloud_method.eq.2 ) then
        if ( cbot_pressure .gt. surface_press ) then
          do_fail = .true.
          message = 'Abort: Cbot pressure > surface pressure'
          write(*,*)message
          stop
        endif
      endif

!  This code is an automatic calculation of bottom pressure
!        maxdrop = 0.999*(surface_press - ctop_pressure)
!        pdiff = 12.5d0 * cloud_opdep
!        pdiff = min(pdiff,maxdrop)
!        cbot_pressure = ctop_pressure + pdiff

!  copy to holding arrays
!  ----------------------

      nlayers_hold = nlayers_total
      nlayers_s    = nlayers_hold
      do n = 0, nlayers_hold
        hold_heights(n) = heights(n)
        hold_press(n)   = press(n)
        hold_logp(n)    = dlog(hold_press(n))
        if ( n .gt. 0 ) then
          hold_airdens(n) = airdens(n)
          hold_dsair(n)   = ds_airdens(n)
          hold_temps(n) = temps(n)
          hold_shape(n) = shapegrid(n)
        endif
      enddo
      do g = 1, ngases
        if ( do_gases(g) ) then
          do n = 1, nlayers_hold
            hold_pr(n,g)      = profiles(n,g)
            hold_prderiv(n,g) = dc_profiles(n,g)
          enddo
        endif
      enddo

!  Find the layers containing cloud levels
!  ---------------------------------------

!  Find layer containing cloud top
!    interpolate to cloud-top pressure

      loop = .true.
      n = 0
      do while (loop)
        n = n + 1
        if ( ctop_pressure .lt. hold_press(n) ) loop = .false.
      enddo
      nlayers_t = n

      log_ptop = dlog(ctop_pressure)
      dlp  = hold_logp(n) - hold_logp(n-1)
      dlpt = log_ptop     - hold_logp(n-1)
      gradh = (hold_heights(n) - hold_heights(n-1))/dlp
      hgt_t = hold_heights(n-1) + gradh *dlpt
      tmp_t = hold_temps(n)
      shp_t = hold_shape(n)
      fac_t = dlpt / dlp
      ad_t  = fac_t * hold_airdens(n)
      add_t = fac_t * hold_dsair(n)
      do g = 1, ngases
        if ( do_gases(g) ) then
          pr_t(g)  = fac_t * hold_pr(n,g)
          prd_t(g) = fac_t * hold_prderiv(n,g)
        endif
      enddo

!  These two lines are for temperatures on levels
!      gradt = (hold_temps(n) - hold_temps(n-1))/dlp
!      tmp_t = hold_temps(n-1) + gradt *dlpt

!  Find layer with cloud bottom (clouds-as-layers only)
!   interpolate to cloud - bottom pressure

      if ( cloud_method .eq. 2 ) then

        loop = .true.
        n = 0
        do while (loop)
          n = n + 1
          if ( cbot_pressure .lt. hold_press(n) ) loop = .false.
        enddo

        nlayers_b = n
        log_pbot = dlog(cbot_pressure)
        dlp  = hold_logp(n) - hold_logp(n-1)
        dlpb = log_pbot     - hold_logp(n-1)
        gradh = (hold_heights(n) - hold_heights(n-1))/dlp
        hgt_b = hold_heights(n-1) + gradh *dlpb
        tmp_b = hold_temps(n)
        shp_b = hold_shape(n)
        fac_b = dlpb / dlp
        ad_b  = fac_b * hold_airdens(n)
        add_b = fac_b * hold_dsair(n)
        do g = 1, ngases
          if ( do_gases(g) ) then
            pr_b(g)  = fac_b * hold_pr(n,g)
            prd_b(g) = fac_b * hold_prderiv(n,g)
          endif
        enddo

!  These two lines are for temperatures on levels
!        gradt = (hold_temps(n) - hold_temps(n-1))/dlp
!        tmp_b = hold_temps(n-1) + gradt *dlpb

!  Assign variables for the parceling

         diffh = cloud_opdep / ( hgt_t - hgt_b ) 
         
      endif

!  debug output

!      write(*,*)nlayers_t, hgt_t,tmp_t, o3_t, ctop_pressure
!      if(do_cloud_layering)write(*,*)&
!            nlayers_b, hgt_b,tmp_b, o3_b, cbot_pressure
!      write(*,*)' '
!      pause

!  Clouds as reflecting boundaries
!  -------------------------------

!  skip this if not appropriate

      if ( cloud_method .eq. 2 ) go to 566

!  total grid to cloud top

      nlayers_cloud = nlayers_t
      heights(0) = hold_heights(0)
      press(0)   = hold_press(0)
      do n = 1, nlayers_t - 1
        heights(n)    = hold_heights(n)
        temps(n)      = hold_temps(n)
        shapegrid(n)  = hold_shape(n)
        press(n)      = hold_press(n)
        airdens(n)    = hold_airdens(n)
        ds_airdens(n) = hold_dsair(n)
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = hold_pr(n,g)
            dc_profiles(n,g) = hold_prderiv(n,g)
          endif
        enddo
      enddo
      n = nlayers_t
      heights(n)    = hgt_t
      temps(n)      = tmp_t
      shapegrid(n)  = shp_t
      press(n)      = ctop_pressure
      airdens(n)    = ad_t
      ds_airdens(n) = add_t
      do g = 1, ngases
        if ( do_gases(g) ) then
          profiles(n,g)    = pr_t(g)
          dc_profiles(n,g) = prd_t(g)
        endif
      enddo
      nlayers_total = nlayers_t

!  total grid: remaining layers to surface. Clearsky

      if ( do_clearsky ) then
        nc  = nlayers_t
        nc1 = nc + 1
        heights(nc1)    = hold_heights(nc)
        temps(nc1)      = hold_temps(nc)
        shapegrid(nc1)  = hold_shape(nc)
        press(nc1)      = hold_press(nc)
        airdens(nc1)    = hold_airdens(nc) - ad_t
        ds_airdens(nc1) = hold_dsair(nc) - add_t
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc1,g)    = hold_pr(nc,g) - pr_t(g)
            dc_profiles(nc1,g) = hold_prderiv(nc,g) - prd_t(g)
          endif
        enddo
        do n = nc1, nlayers_s
          n1 = n + 1
          heights(n1)    = hold_heights(n)
          temps(n1)      = hold_temps(n)
          shapegrid(n1)  = hold_shape(n)
          press(n1)      = hold_press(n)
          airdens(n1)    = hold_airdens(n)
          ds_airdens(n1) = hold_dsair(n)
         do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n1,g)    = hold_pr(n,g)
              dc_profiles(n1,g) = hold_prderiv(n,g)
            endif
          enddo
        enddo
        nlayers_total = nlayers_s + 1
      endif

!  finish

      return

!  Continuation point for clouds as layers

 566  continue

!  clouds as layers grid
!  ---------------------

!  Assign layers down to level immediately above cloud top

      heights(0) = hold_heights(0)
      press(0)   = hold_press(0)
      do n = 1, nlayers_t - 1
        heights(n)    = hold_heights(n)
        temps(n)      = hold_temps(n)
        shapegrid(n)  = hold_shape(n)
        press(n)      = hold_press(n)
        airdens(n)    = hold_airdens(n)
        ds_airdens(n) = hold_dsair(n)
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = hold_pr(n,g)
            dc_profiles(n,g) = hold_prderiv(n,g)
          endif
        enddo
      enddo

!  Assign reduced layer above cloud top

      n = nlayers_t
      heights(n)    = hgt_t
      temps(n)      = tmp_t
      shapegrid(n)  = shp_t
      press(n)      = ctop_pressure
      airdens(n)    = ad_t
      ds_airdens(n) = add_t
      do g = 1, ngases
        if ( do_gases(g) ) then
          profiles(n,g)    = pr_t(g)
          dc_profiles(n,g) = prd_t(g)
        endif
      enddo

!  If the cloud is within one grid layer

      if (nlayers_b .eq. nlayers_t ) then

!....assign layer to bottom cloud level. THE CLOUD LAYER.

        n = n + 1
        heights(n)    = hgt_b
        temps(n)      = tmp_b
        shapegrid(n)  = shp_b
        press(n)      = cbot_pressure
        airdens(n)    = ad_b - ad_t
        ds_airdens(n) = add_b - add_t
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = pr_b(g)  - pr_t(g)
            dc_profiles(n,g) = prd_b(g) - prd_t(g)
          endif
        enddo
        if ( .not. do_clearsky ) then
          cloudflag(n) = .true.
          cparcel(n)   = cloud_opdep
        endif

!....additional layers between cloud-bottom and surface

        nc  = nlayers_b
        nc2 = nc + 2
        heights(nc2)    = hold_heights(nc)
        temps(nc2)      = hold_temps(nc)
        shapegrid(nc2)  = hold_shape(nc)
        press(nc2)      = hold_press(nc)
        airdens(nc2)    = hold_airdens(nc) - ad_b
        ds_airdens(nc2) = hold_dsair(nc)   - add_b
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc2,g)    = hold_pr(nc,g) - pr_b(g)
            dc_profiles(nc2,g) = hold_prderiv(nc,g) - prd_b(g)
          endif
        enddo
        do n = nc+1, nlayers_s
          n2 = n + 2
          heights(n2)    = hold_heights(n)
          shapegrid(n2)  = hold_shape(n)
          temps(n2)      = hold_temps(n)
          press(n2)      = hold_press(n)
          airdens(n2)    = hold_airdens(n)
          ds_airdens(n2) = hold_dsair(n)
          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n2,g)    = hold_pr(n,g)
              dc_profiles(n2,g) = hold_prderiv(n,g)
            endif
          enddo
        enddo
        nlayers_total = nlayers_s + 2

!  If the cloud straddles a number of layers

      else

!.....assign first layer below the cloud top

        nc  = nlayers_t
        nc1 = nc + 1
        heights(nc1)    = hold_heights(nc)
        temps(nc1)      = hold_temps(nc)
        shapegrid(nc1)  = hold_shape(nc)
        press(nc1)      = hold_press(nc)
        airdens(nc1)    = hold_airdens(nc) - ad_t
        ds_airdens(nc1) = hold_dsair(nc) - add_t
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc1,g)    = hold_pr(nc,g) - pr_t(g)
            dc_profiles(nc1,g) = hold_prderiv(nc,g) - prd_t(g)
          endif
        enddo
        if ( .not. do_clearsky ) then
          cloudflag(nc1) = .true.
          cparcel(nc1)   = ( hgt_t - hold_heights(nc) ) * diffh
        endif

!....assign intervening layers to level immediately above cloud bottom

        do n = nc1, nlayers_b - 1
          n1 = n + 1
          heights(n1)    = hold_heights(n)
          shapegrid(n1)  = hold_shape(n)
          temps(n1)      = hold_temps(n)
          press(n1)      = hold_press(n)
          airdens(n1)    = hold_airdens(n)
          ds_airdens(n1) = hold_dsair(n)
          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n1,g)    = hold_pr(n,g)
              dc_profiles(n1,g) = hold_prderiv(n,g)
            endif
          enddo
          if ( .not. do_clearsky ) then
            cloudflag(n1) = .true.
            cparcel(n1)   = (hold_heights(n-1)-hold_heights(n)) * diffh
          endif
        enddo

!....assign final cloud layer to cloud bottom level

        n = nlayers_b + 1
        heights(n)    = hgt_b
        temps(n)      = tmp_b
        shapegrid(n)  = shp_b
        press(n)      = cbot_pressure
        airdens(n)    = ad_b
        ds_airdens(n) = add_b
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = pr_b(g)
            dc_profiles(n,g) = prd_b(g)
          endif
        enddo
        if ( .not. do_clearsky ) then
          cloudflag(n) = .true.
          cparcel(n)   = ( hold_heights(nlayers_b-1) - hgt_b ) * diffh
        endif

!  ...more than one additional intervening layers to the surface

        nc = nlayers_b
        nc2 = nc + 2
        heights(nc2)    = hold_heights(nc)
        temps(nc2)      = hold_temps(nc)
        shapegrid(nc2)  = hold_shape(nc)
        press(nc2)      = hold_press(nc)
        airdens(nc2)    = hold_airdens(nc) - ad_b
        ds_airdens(nc2) = hold_dsair(nc)   - add_b
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc2,g)    = hold_pr(nc,g) - pr_b(g)
            dc_profiles(nc2,g) = hold_prderiv(nc,g) - prd_b(g)
          endif
        enddo
        do n = nc+1, nlayers_s
          n2 = n + 2
          heights(n2)    = hold_heights(n)
          temps(n2)      = hold_temps(n)
          shapegrid(n2)  = hold_shape(n)
          press(n2)      = hold_press(n)
          airdens(n2)    = hold_airdens(n)
          ds_airdens(n2) = hold_dsair(n)
          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n2,g)    = hold_pr(n,g)
              dc_profiles(n2,g) = hold_prderiv(n,g)
            endif
          enddo
        enddo
        nlayers_total = nlayers_s + 2

      endif

!  cloud grid

      nlayers_cloud = nlayers_total
       
!  Finish

      return
      end subroutine PTH_REGRIDDING_NOAER


!

      subroutine PTH_REGRIDDING_GENERAL                     &
       ( maxlayers, maxgases, maxaerwfs, ngases, do_gases,  &
         do_aerosols, do_aer_Jacobians, n_aerosol_wfs,      &
         cloud_method, do_clearsky, cloud_opdep,            &
         surface_press, ctop_pressure, cbot_pressure,       &
         nlayers_total, press, heights, temps, shapegrid,   &      ! MODIFIED INPUTS
         profiles, dc_profiles, airdens, ds_airdens,        &      ! MODIFIED INPUTS
         aerflag, aertau_unscaled, L_aertau_unscaled,       &      ! MODIFIED INPUTS
         nlayers_cloud, cloudflag, cparcel )                       !  NEW OUTPUTS

      implicit none

!  dimensioning

      integer          maxlayers, maxgases,  maxaerwfs 

!  gas control

      integer          ngases
      logical          do_gases ( maxgases)

!  Aerosol control

      logical          do_aerosols
      logical          do_aer_Jacobians
      integer          n_aerosol_wfs

!  control inputs (index, surface and cloud information)

      logical          do_clearsky
      integer          cloud_method
      real(fpk)        surface_press
      real(fpk)        cloud_opdep
      real(fpk)        ctop_pressure,cbot_pressure

!  Modified Grid outputs

      integer          nlayers_total
      real(fpk)        press     (0:maxlayers)
      real(fpk)        temps     (  maxlayers)
      real(fpk)        heights   (0:maxlayers)
      real(fpk)        shapegrid (  maxlayers)

      real(fpk)        profiles    (maxlayers,maxgases)
      real(fpk)        dc_profiles (maxlayers,maxgases)
      real(fpk)        airdens     (maxlayers)
      real(fpk)        ds_airdens  (maxlayers)

      LOGICAL          aerflag              (maxlayers)
      real(fpk)        aertau_unscaled      (maxlayers)
      real(fpk)        L_aertau_unscaled    (maxaerwfs,maxlayers)
 
!  Completely new outputs for cloud information

      integer          nlayers_cloud
      logical          cloudflag (maxlayers)
      real(fpk)        cparcel   (maxlayers)

!  status

      character*70    message
      logical         do_fail

!  local variables
!  ---------------

!  holding arrays 

      integer          maxlayers_hold, maxgases_hold, maxaerwfs_hold, nlayers_hold
      parameter        ( maxlayers_hold = 101 )
      parameter        ( maxgases_hold   = 2  )
      parameter        ( maxaerwfs_hold  = 15  )

      real(fpk)        hold_press  (0:maxlayers_hold)
      real(fpk)        hold_temps  (  maxlayers_hold)
      real(fpk)        hold_shape  (  maxlayers_hold)

      real(fpk)        hold_heights(0:maxlayers_hold)
      real(fpk)        hold_logp   (0:maxlayers_hold)
      real(fpk)        hold_airdens(  maxlayers_hold)
      real(fpk)        hold_dsair  (  maxlayers_hold)

      real(fpk)        hold_pr      (maxlayers_hold,maxgases_hold)
      real(fpk)        hold_prderiv (maxlayers_hold,maxgases_hold)

      logical          hold_flag     (maxlayers_hold)
      real(fpk)        hold_aer      (maxlayers_hold)
      real(fpk)        hold_aerderiv (maxaerwfs_hold,maxlayers_hold)

      real(fpk)        pr_t  (maxgases_hold), ad_t
      real(fpk)        prd_t (maxgases_hold), add_t
      real(fpk)        pr_b  (maxgases_hold), ad_b
      real(fpk)        prd_b (maxgases_hold), add_b

      real(fpk)        aer_t, aerd_t(maxaerwfs_hold)
      real(fpk)        aer_b, aerd_b(maxaerwfs_hold)

!  other help variables

      integer          n, nc, n1,n2, nc1, nc2, g, q
      integer          nlayers_s, nlayers_t, nlayers_b
      logical          loop, flag_b, flag_t
      real(fpk)        gradh, dlp, dlpb, dlpt, diffh, fac_t, fac_b
      real(fpk)        hgt_b, hgt_t, tmp_b, tmp_t, shp_t, shp_b
      real(fpk)        log_ptop, log_pbot

!  Initialise (important!)
!  ----------

      do_fail    = .false.
      message    = ' '

!  initialise nlayers_cloud and cloudflag and cparcel
!       (the only new variables to be set)

      nlayers_cloud = 0
      do n = 1, maxlayers
        cloudflag(n) = .false.
        cparcel(n)   = 0.0d0
      enddo

!  Nothing to do if Clearsky, Cloud_method = 0

      if ( cloud_method .eq. 0 ) return

!  Checking input
!  --------------

!  Check the cloud-top pressure is not out of bounds.

      if ( ctop_pressure .gt. surface_press ) then
        do_fail = .true.
        message = 'Abort: Ctop pressure > surface pressure'
        write(*,*)message
        stop
      endif

!  Check the cloud bottom pressure (clouds as layers only)

      if ( cloud_method.eq.2 ) then
        if ( cbot_pressure .gt. surface_press ) then
          do_fail = .true.
          message = 'Abort: Cbot pressure > surface pressure'
          write(*,*)message
          stop
        endif
      endif

!  This code is an automatic calculation of bottom pressure
!        maxdrop = 0.999*(surface_press - ctop_pressure)
!        pdiff = 12.5d0 * cloud_opdep
!        pdiff = min(pdiff,maxdrop)
!        cbot_pressure = ctop_pressure + pdiff

!  copy to holding arrays
!  ----------------------

      nlayers_hold = nlayers_total
      nlayers_s    = nlayers_hold
      do n = 0, nlayers_hold
        hold_heights(n) = heights(n)
        hold_press(n)   = press(n)
        hold_logp(n)    = dlog(hold_press(n))
        if ( n .gt. 0 ) then
          hold_airdens(n) = airdens(n)
          hold_dsair(n)   = ds_airdens(n)
          hold_temps(n) = temps(n)
          hold_shape(n) = shapegrid(n)
        endif
      enddo

      do g = 1, ngases
        if ( do_gases(g) ) then
          do n = 1, nlayers_hold
            hold_pr(n,g)      = profiles(n,g)
            hold_prderiv(n,g) = dc_profiles(n,g)
          enddo
        endif
      enddo

      if ( do_aerosols ) then
         do n = 1, nlayers_hold
            hold_flag(n) = aerflag(n)
            hold_aer(n)  = aertau_unscaled(n)
            if ( do_aer_Jacobians ) then
               do q = 1, n_aerosol_wfs
                  hold_aerderiv(q,n) = L_aertau_unscaled(q,n)
               enddo
            endif
         enddo
      endif

!  Find the layers containing cloud levels
!  ---------------------------------------

!  Find layer containing cloud top
!    interpolate to cloud-top pressure

      loop = .true.
      n = 0
      do while (loop)
        n = n + 1
        if ( ctop_pressure .lt. hold_press(n) ) loop = .false.
      enddo
      nlayers_t = n

      log_ptop = dlog(ctop_pressure)
      dlp  = hold_logp(n) - hold_logp(n-1)
      dlpt = log_ptop     - hold_logp(n-1)
      gradh = (hold_heights(n) - hold_heights(n-1))/dlp
      hgt_t = hold_heights(n-1) + gradh *dlpt
      tmp_t = hold_temps(n)
      shp_t = hold_shape(n)
      fac_t = dlpt / dlp
      ad_t  = fac_t * hold_airdens(n)
      add_t = fac_t * hold_dsair(n)

      do g = 1, ngases
        if ( do_gases(g) ) then
          pr_t(g)  = fac_t * hold_pr(n,g)
          prd_t(g) = fac_t * hold_prderiv(n,g)
        endif
      enddo

      if ( do_aerosols ) then
         flag_t = hold_flag(n)
         aer_t  = fac_t * hold_aer(n)
         if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              aerd_t(q) = fac_t * hold_aerderiv(q,n)
            enddo
         endif
      endif

!  These two lines are for temperatures on levels
!      gradt = (hold_temps(n) - hold_temps(n-1))/dlp
!      tmp_t = hold_temps(n-1) + gradt *dlpt

!  Find layer with cloud bottom (clouds-as-layers only)
!   interpolate to cloud - bottom pressure

      if ( cloud_method .eq. 2 ) then

        loop = .true.
        n = 0
        do while (loop)
          n = n + 1
          if ( cbot_pressure .lt. hold_press(n) ) loop = .false.
        enddo

        nlayers_b = n
        log_pbot = dlog(cbot_pressure)
        dlp  = hold_logp(n) - hold_logp(n-1)
        dlpb = log_pbot     - hold_logp(n-1)
        gradh = (hold_heights(n) - hold_heights(n-1))/dlp
        hgt_b = hold_heights(n-1) + gradh *dlpb
        tmp_b = hold_temps(n)
        shp_b = hold_shape(n)
        fac_b = dlpb / dlp
        ad_b  = fac_b * hold_airdens(n)
        add_b = fac_b * hold_dsair(n)
        do g = 1, ngases
          if ( do_gases(g) ) then
            pr_b(g)  = fac_b * hold_pr(n,g)
            prd_b(g) = fac_b * hold_prderiv(n,g)
          endif
        enddo

        if ( do_aerosols ) then
           flag_b = hold_flag(n)
           aer_b  = fac_b * hold_aer(n)
           if ( do_aer_Jacobians ) then
              do q = 1, n_aerosol_wfs
                aerd_b(q) = fac_b * hold_aerderiv(q,n)
              enddo
           endif
        endif

!  These two lines are for temperatures on levels
!        gradt = (hold_temps(n) - hold_temps(n-1))/dlp
!        tmp_b = hold_temps(n-1) + gradt *dlpb

!  Assign variables for the parceling

         diffh = cloud_opdep / ( hgt_t - hgt_b ) 
         
      endif

!  debug output

!      write(*,*)nlayers_t, hgt_t,tmp_t, o3_t, ctop_pressure
!      if(do_cloud_layering)write(*,*)&
!            nlayers_b, hgt_b,tmp_b, o3_b, cbot_pressure
!      write(*,*)' '
!      pause

!  Clouds as reflecting boundaries
!  -------------------------------

!  skip this if not appropriate

      if ( cloud_method .eq. 2 ) go to 566

!  total grid to cloud top

      nlayers_cloud = nlayers_t
      heights(0) = hold_heights(0)
      press(0)   = hold_press(0)
      do n = 1, nlayers_t - 1
        heights(n)    = hold_heights(n)
        temps(n)      = hold_temps(n)
        shapegrid(n)  = hold_shape(n)
        press(n)      = hold_press(n)
        airdens(n)    = hold_airdens(n)
        ds_airdens(n) = hold_dsair(n)
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = hold_pr(n,g)
            dc_profiles(n,g) = hold_prderiv(n,g)
          endif
        enddo
        if (do_aerosols ) then
          aerflag(n)         = hold_flag(n)
          aertau_unscaled(n) = hold_aer(n)
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,n) = hold_aerderiv(q,n)
            enddo
          endif
        endif
      enddo

!  reduced layer above cloud-top

      n = nlayers_t
      heights(n)    = hgt_t
      temps(n)      = tmp_t
      shapegrid(n)  = shp_t
      press(n)      = ctop_pressure
      airdens(n)    = ad_t
      ds_airdens(n) = add_t

      do g = 1, ngases
        if ( do_gases(g) ) then
          profiles(n,g)    = pr_t(g)
          dc_profiles(n,g) = prd_t(g)
        endif
      enddo

      if (do_aerosols ) then
        aerflag(n)         = flag_t
        aertau_unscaled(n) = aer_t
        if ( do_aer_Jacobians ) then
          do q = 1, n_aerosol_wfs
            L_aertau_unscaled(q,n) = aerd_t(q)
          enddo
        endif
      endif

      nlayers_total = nlayers_t

!   Clear-sky situation, fill in rest of atmosphere

      if ( do_clearsky ) then

!  Partial layer immediately below cloud-top

        nc  = nlayers_t
        nc1 = nc + 1

        heights(nc1)    = hold_heights(nc)
        temps(nc1)      = hold_temps(nc)
        shapegrid(nc1)  = hold_shape(nc)
        press(nc1)      = hold_press(nc)
        airdens(nc1)    = hold_airdens(nc) - ad_t
        ds_airdens(nc1) = hold_dsair(nc) - add_t

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc1,g)    = hold_pr(nc,g) - pr_t(g)
            dc_profiles(nc1,g) = hold_prderiv(nc,g) - prd_t(g)
          endif
        enddo

        if (do_aerosols ) then
         aerflag(nc1)         = hold_flag(nc)
         aertau_unscaled(nc1) = hold_aer(nc) - aer_t
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,nc1) = hold_aerderiv(q,nc) - aerd_t(q)
            enddo
          endif
        endif

!  remaining layers to surface

        do n = nc1, nlayers_s

          n1 = n + 1
          heights(n1)    = hold_heights(n)
          temps(n1)      = hold_temps(n)
          shapegrid(n1)  = hold_shape(n)
          press(n1)      = hold_press(n)
          airdens(n1)    = hold_airdens(n)
          ds_airdens(n1) = hold_dsair(n)

          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n1,g)    = hold_pr(n,g)
              dc_profiles(n1,g) = hold_prderiv(n,g)
            endif
          enddo

          if ( do_aerosols ) then
            aerflag(n1)         = hold_flag(n)
            aertau_unscaled(n1) = hold_aer(n)
            if ( do_aer_Jacobians ) then
              do q = 1, n_aerosol_wfs
                L_aertau_unscaled(q,n1) =  hold_aerderiv(q,n) 
              enddo
            endif
          endif

        enddo

        nlayers_total = nlayers_s + 1
 
     endif

!  finish

      return

!  Continuation point for clouds as layers

 566  continue

!  clouds as layers grid
!  ---------------------

!  Assign layers down to level immediately above cloud top

      heights(0) = hold_heights(0)
      press(0)   = hold_press(0)

      do n = 1, nlayers_t - 1

        heights(n)    = hold_heights(n)
        temps(n)      = hold_temps(n)
        shapegrid(n)  = hold_shape(n)
        press(n)      = hold_press(n)
        airdens(n)    = hold_airdens(n)
        ds_airdens(n) = hold_dsair(n)

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = hold_pr(n,g)
            dc_profiles(n,g) = hold_prderiv(n,g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(n)         = hold_flag(n)
          aertau_unscaled(n) = hold_aer(n)
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,n) = hold_aerderiv(q,n)
            enddo
          endif
        endif

      enddo

!  Assign reduced layer above cloud top

      n = nlayers_t
      heights(n)    = hgt_t
      temps(n)      = tmp_t
      shapegrid(n)  = shp_t
      press(n)      = ctop_pressure
      airdens(n)    = ad_t
      ds_airdens(n) = add_t

      do g = 1, ngases
        if ( do_gases(g) ) then
          profiles(n,g)    = pr_t(g)
          dc_profiles(n,g) = prd_t(g)
        endif
      enddo

      if (do_aerosols ) then
        aerflag(n)         = flag_t
        aertau_unscaled(n) = aer_t
        if ( do_aer_Jacobians ) then
          do q = 1, n_aerosol_wfs
            L_aertau_unscaled(q,n) = aerd_t(q)
          enddo
        endif
      endif

!  If the cloud is within one grid layer
!  =====================================

      if ( nlayers_b .eq. nlayers_t ) then

!....assign layer to bottom cloud level. THE CLOUD LAYER.

        n = n + 1
        heights(n)    = hgt_b
        temps(n)      = tmp_b
        shapegrid(n)  = shp_b
        press(n)      = cbot_pressure
        airdens(n)    = ad_b - ad_t
        ds_airdens(n) = add_b - add_t

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = pr_b(g)  - pr_t(g)
            dc_profiles(n,g) = prd_b(g) - prd_t(g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(n)         = flag_t
          aertau_unscaled(n) = aer_b - aer_t
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,n) = aerd_b(q) - aerd_t(q)
            enddo
          endif
        endif

        if ( .not. do_clearsky ) then
          cloudflag(n) = .true.
          cparcel(n)   = cloud_opdep
        endif

!....additional layer between cloud-bottom and next lowest grid level

        nc  = nlayers_b
        nc2 = nc + 2
        heights(nc2)    = hold_heights(nc)
        temps(nc2)      = hold_temps(nc)
        shapegrid(nc2)  = hold_shape(nc)
        press(nc2)      = hold_press(nc)
        airdens(nc2)    = hold_airdens(nc) - ad_b
        ds_airdens(nc2) = hold_dsair(nc)   - add_b

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc2,g)    = hold_pr(nc,g) - pr_b(g)
            dc_profiles(nc2,g) = hold_prderiv(nc,g) - prd_b(g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(nc2)         = hold_flag(nc)
          aertau_unscaled(nc2) = hold_aer(nc) - aer_b
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,nc2) = hold_aerderiv(q,nc) - aerd_b(q)
            enddo
          endif
        endif

!....additional layers to surface, follwing original grid

        do n = nc+1, nlayers_s

          n2 = n + 2
          heights(n2)    = hold_heights(n)
          shapegrid(n2)  = hold_shape(n)
          temps(n2)      = hold_temps(n)
          press(n2)      = hold_press(n)
          airdens(n2)    = hold_airdens(n)
          ds_airdens(n2) = hold_dsair(n)

          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n2,g)    = hold_pr(n,g)
              dc_profiles(n2,g) = hold_prderiv(n,g)
            endif
          enddo

          if ( do_aerosols ) then
            aerflag(n2)         = hold_flag(n)
            aertau_unscaled(n2) = hold_aer(n)
            if ( do_aer_Jacobians ) then
              do q = 1, n_aerosol_wfs
                L_aertau_unscaled(q,n2) =  hold_aerderiv(q,n) 
              enddo
            endif
          endif

        enddo

        nlayers_total = nlayers_s + 2

!  If the cloud straddles a number of layers
!  =========================================

      else

!.....assign first layer below the cloud top

        nc  = nlayers_t
        nc1 = nc + 1
        heights(nc1)    = hold_heights(nc)
        temps(nc1)      = hold_temps(nc)
        shapegrid(nc1)  = hold_shape(nc)
        press(nc1)      = hold_press(nc)
        airdens(nc1)    = hold_airdens(nc) - ad_t
        ds_airdens(nc1) = hold_dsair(nc) - add_t

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc1,g)    = hold_pr(nc,g) - pr_t(g)
            dc_profiles(nc1,g) = hold_prderiv(nc,g) - prd_t(g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(nc1)         = hold_flag(nc)
          aertau_unscaled(nc1) = hold_aer(nc) - aer_t
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,nc1) = hold_aerderiv(q,nc) - aerd_t(q)
            enddo
          endif
        endif

        if ( .not. do_clearsky ) then
          cloudflag(nc1) = .true.
          cparcel(nc1)   = ( hgt_t - hold_heights(nc) ) * diffh
        endif

!....assign intervening layers to level immediately above cloud bottom

        do n = nc1, nlayers_b - 1

          n1 = n + 1
          heights(n1)    = hold_heights(n)
          shapegrid(n1)  = hold_shape(n)
          temps(n1)      = hold_temps(n)
          press(n1)      = hold_press(n)
          airdens(n1)    = hold_airdens(n)
          ds_airdens(n1) = hold_dsair(n)

          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n1,g)    = hold_pr(n,g)
              dc_profiles(n1,g) = hold_prderiv(n,g)
            endif
          enddo

          if ( do_aerosols ) then
            aerflag(n1)         = hold_flag(n)
            aertau_unscaled(n1) = hold_aer(n)
            if ( do_aer_Jacobians ) then
              do q = 1, n_aerosol_wfs
                L_aertau_unscaled(q,n1) =  hold_aerderiv(q,n) 
              enddo
            endif
          endif

          if ( .not. do_clearsky ) then
            cloudflag(n1) = .true.
            cparcel(n1)   = (hold_heights(n-1)-hold_heights(n)) * diffh
          endif

        enddo

!....assign final cloud layer to cloud bottom level

        n = nlayers_b + 1
        heights(n)    = hgt_b
        temps(n)      = tmp_b
        shapegrid(n)  = shp_b
        press(n)      = cbot_pressure
        airdens(n)    = ad_b
        ds_airdens(n) = add_b

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = pr_b(g)
            dc_profiles(n,g) = prd_b(g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(n)         = flag_b
          aertau_unscaled(n) = aer_b
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,n) = aerd_b(q)
            enddo
          endif
        endif

        if ( .not. do_clearsky ) then
          cloudflag(n) = .true.
          cparcel(n)   = ( hold_heights(nlayers_b-1) - hgt_b ) * diffh
        endif

!  ...additional intervening layer to the next grid level

        nc = nlayers_b
        nc2 = nc + 2
        heights(nc2)    = hold_heights(nc)
        temps(nc2)      = hold_temps(nc)
        shapegrid(nc2)  = hold_shape(nc)
        press(nc2)      = hold_press(nc)
        airdens(nc2)    = hold_airdens(nc) - ad_b
        ds_airdens(nc2) = hold_dsair(nc)   - add_b

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc2,g)    = hold_pr(nc,g) - pr_b(g)
            dc_profiles(nc2,g) = hold_prderiv(nc,g) - prd_b(g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(nc2)         = hold_flag(nc)
          aertau_unscaled(nc2) = hold_aer(nc) - aer_b
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,nc2) = hold_aerderiv(q,nc) - aerd_b(q)
            enddo
          endif
        endif

!  ...additional gridded layers to surface

        do n = nc+1, nlayers_s

          n2 = n + 2
          heights(n2)    = hold_heights(n)
          temps(n2)      = hold_temps(n)
          shapegrid(n2)  = hold_shape(n)
          press(n2)      = hold_press(n)
          airdens(n2)    = hold_airdens(n)
          ds_airdens(n2) = hold_dsair(n)

          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n2,g)    = hold_pr(n,g)
              dc_profiles(n2,g) = hold_prderiv(n,g)
            endif
          enddo

          if ( do_aerosols ) then
            aerflag(n2)         = hold_flag(n)
            aertau_unscaled(n2) = hold_aer(n)
            if ( do_aer_Jacobians ) then
              do q = 1, n_aerosol_wfs
                L_aertau_unscaled(q,n2) =  hold_aerderiv(q,n) 
              enddo
            endif
          endif

        enddo

        nlayers_total = nlayers_s + 2

      endif

!  cloud grid

      nlayers_cloud = nlayers_total

!  Finish

      return
      end subroutine PTH_REGRIDDING_GENERAL

!  ALTERNATIVE Height-based Regridding routines

!   ** PURE (FIXED) INPTUS

!maxlayers, maxgases (dimensioning)
!ngases = 1
!do_gases(1) = .true.
!cloud_method   = 2 (mie scattering cloud)
!do_clearsky, Flag for doing a clear sky calculation, with (Method = 1/2) or without (Method = 0) regridding
!cloud_opdep.  Cloud optical thickness at one wavelength, GOING to be PARCELLED (in cparcel)
!surface_press


!  ** CLOUD INPUTS

! IF ( ALTERNATIVE )       then CTOP/CBOT Heights   are Given, Pressures Calculated
! IF ( .not. ALTERNATIVE ) then CTOP/CBOT Pressures are Given, Heights   Calculated
!      ctop_height, ctop_pressure
!      cbot_height, cbot_pressure

!  ** MODIFIED (REGRIDDED) INPUTS

!nlayers_total. Originally N layers, come out (usually) with N+2
!press   = Level pressures
!heights = level heights
!temps   = LAYER temperatures
!shapegrid = LAYER Temperature-shift functions, set these to 0.0 throughout (not useful)

!profiles(gas =1) = ozone profile
!profiles_dc = derivative of the profile w.r.t its total column. Set this to 0.0 (not needed)

!airdens = air column desities, each layer
!ds_airdens = derivative of airdens w.r.t. T-shift. Set to 0.0d0 (Not required)

!  PURE OUTPUTS

!nlayers_cloud. Only for Method 1
!cloudflag(n) = (Method 2) .true. if n is the cloud layer. SAME as the aerosol FLAG.
!cparcel(n)   = (Method 2) loading. Same as AEROSOL_LOADING (layer optical thickness values)
!                        in your master routine.

      subroutine PTH_REGRIDDING_ALT_NOAER                        &
       ( maxlayers, maxgases, ngases, do_gases, surface_press,   & ! PURE INPUTS
         alternative, cloud_method, do_clearsky, cloud_opdep,    & ! PURE INPUTS
         ctop_height, cbot_height, ctop_pressure, cbot_pressure, & ! MODIFIED/PURE INPUTS
         nlayers_total, press, heights, temps, shapegrid,        & ! MODIFIED INPUTS
         profiles, dc_profiles, airdens, ds_airdens,             & ! MODIFIED INPUTS
         nlayers_cloud, cloudflag, cparcel,                      & ! NEW OUTPUTS
         do_fail, message )                                        ! STATUS

      implicit none

!  dimensioning

      integer, intent(in)          :: maxlayers, maxgases

!  gas control

      integer, intent(in)          :: ngases
      logical, intent(in)          :: do_gases ( maxgases)

!  control inputs and cloud information

      logical, intent(in)          :: do_clearsky
      integer, intent(in)          :: cloud_method
      real(fpk)       , intent(in) :: cloud_opdep
      real(fpk)       , intent(in) :: surface_press
      logical, intent(in)          :: alternative

!  Cloud boundaries (in or out)

      real(fpk)       , intent(inout) :: ctop_height,   cbot_height
      real(fpk)       , intent(inout) :: ctop_pressure, cbot_pressure

!  Modified Grid outputs

      integer         , intent(inout) :: nlayers_total
      real(fpk)       , intent(inout) :: press     (0:maxlayers)
      real(fpk)       , intent(inout) :: temps     (  maxlayers)
      real(fpk)       , intent(inout) :: heights   (0:maxlayers)
      real(fpk)       , intent(inout) :: shapegrid (  maxlayers)

      real(fpk)       , intent(inout) :: profiles    (maxlayers,maxgases)
      real(fpk)       , intent(inout) :: dc_profiles (maxlayers,maxgases)
      real(fpk)       , intent(inout) :: airdens     (maxlayers)
      real(fpk)       , intent(inout) :: ds_airdens  (maxlayers)

!  Completely new outputs for cloud information

      integer         , intent(out) :: nlayers_cloud
      logical         , intent(out) :: cloudflag (maxlayers)
      real(fpk)       , intent(out) :: cparcel   (maxlayers)

!  status

      character*(*), intent(out) :: message
      logical      , intent(out) :: do_fail

!  local variables
!  ---------------

!  holding arrays 

      integer          maxlayers_hold, maxgases_hold, nlayers_hold
      parameter        ( maxlayers_hold = 201 )
      parameter        ( maxgases_hold  = 2  )

      real(fpk)        hold_press  (0:maxlayers_hold)
      real(fpk)        hold_temps  (  maxlayers_hold)
      real(fpk)        hold_shape  (  maxlayers_hold)

      real(fpk)        hold_heights(0:maxlayers_hold)
      real(fpk)        hold_logp   (0:maxlayers_hold)
      real(fpk)        hold_airdens(  maxlayers_hold)
      real(fpk)        hold_dsair  (  maxlayers_hold)

      real(fpk)        hold_pr      (maxlayers_hold,maxgases_hold)
      real(fpk)        hold_prderiv (maxlayers_hold,maxgases_hold)

      real(fpk)        pr_t  (maxgases_hold), ad_t
      real(fpk)        prd_t (maxgases_hold), add_t
      real(fpk)        pr_b  (maxgases_hold), ad_b
      real(fpk)        prd_b (maxgases_hold), add_b

!  other help variables

      integer          n, nc, n1,n2, nc1, nc2, g
      integer          nlayers_s, nlayers_t, nlayers_b
      logical          loop
      real(fpk)        gradh, dlp, dlpb, dlpt, diffh, fac_t, fac_b
      real(fpk)        hgt_b, hgt_t, tmp_b, tmp_t, shp_t, shp_b
      real(fpk)        log_ptop, log_pbot, prs_t, prs_b

!  Initialise (important!)
!  ----------

      do_fail    = .false.
      message    = ' '

!  initialise nlayers_cloud and cloudflag and cparcel
!       (the only new variables to be set)

      nlayers_cloud = 0
      do n = 1, maxlayers
        cloudflag(n) = .false.
        cparcel(n)   = 0.0d0
      enddo

!  Nothing to do if Clearsky, Cloud_method = 0

      if ( cloud_method .eq. 0 ) return

!  Checking inputs
!  ---------------

!  Check the cloud-top pressure is not out of bounds.

      if ( .not. alternative ) then
        if ( ctop_pressure .gt. surface_press ) then
          do_fail = .true.
          message = 'Abort: Ctop pressure > surface pressure'
          return
        endif
      endif

!  Check the cloud bottom pressure (clouds as layers only)

      if ( .not. alternative ) then
        if ( cloud_method.eq.2 ) then
          if ( cbot_pressure .gt. surface_press ) then
            do_fail = .true.
            message = 'Abort: Cbot pressure > surface pressure'
            return
          endif
        endif
      endif

!  Check the cloud-top height is not out of bounds.

      if ( alternative ) then
        if ( ctop_height .lt. heights(nlayers_total)  ) then
          do_fail = .true.
          message = 'Abort: Ctop height below surface height'
          return
        endif
      endif

!  Check the cloud bottom height (clouds as layers only)

      if ( alternative ) then
        if ( cloud_method.eq.2 ) then
          if ( cbot_height .lt. heights(nlayers_total) ) then
            do_fail = .true.
            message = 'Abort: Cbot height below surface height'
            return
          endif
        endif
      endif

!  copy to holding arrays
!  ----------------------

      nlayers_hold = nlayers_total
      nlayers_s    = nlayers_hold
      do n = 0, nlayers_hold
        hold_heights(n) = heights(n)
        hold_press(n)   = press(n)
        hold_logp(n)    = dlog(hold_press(n))
        if ( n .gt. 0 ) then
          hold_airdens(n) = airdens(n)
          hold_dsair(n)   = ds_airdens(n)
          hold_temps(n) = temps(n)
          hold_shape(n) = shapegrid(n)
        endif
      enddo
      do g = 1, ngases
        if ( do_gases(g) ) then
          do n = 1, nlayers_hold
            hold_pr(n,g)      = profiles(n,g)
            hold_prderiv(n,g) = dc_profiles(n,g)
          enddo
        endif
      enddo

!  Find the layers containing cloud levels
!  ---------------------------------------

!  Find layer containing cloud top
!    interpolate to cloud-top pressure

      loop = .true.
      n = 0
      do while (loop)
        n = n + 1
        if (     alternative.and.(ctop_height   .gt. hold_heights(n))) loop = .false.
        if (.not.alternative.and.(ctop_pressure .lt. hold_press(n))  ) loop = .false.
      enddo
      nlayers_t = n

      dlp   = hold_logp(n) - hold_logp(n-1)
      gradh = (hold_heights(n) - hold_heights(n-1))/dlp
      if ( alternative ) then
         hgt_t = ctop_height
         dlpt  = ( hgt_t - hold_heights(n-1) ) / gradh
         log_ptop = dlpt + hold_logp(n-1)
         prs_t = dexp(log_ptop)
         ctop_pressure = prs_t
      else
         log_ptop = dlog(ctop_pressure)
         dlpt = log_ptop     - hold_logp(n-1)
         hgt_t = hold_heights(n-1) + gradh *dlpt
         prs_t = ctop_pressure
         ctop_height = hgt_t
      endif

      tmp_t = hold_temps(n)
      shp_t = hold_shape(n)
      fac_t = dlpt / dlp
      ad_t  = fac_t * hold_airdens(n)
      add_t = fac_t * hold_dsair(n)

      do g = 1, ngases
        if ( do_gases(g) ) then
          pr_t(g)  = fac_t * hold_pr(n,g)
          prd_t(g) = fac_t * hold_prderiv(n,g)
        endif
      enddo

!  These two lines are for temperatures on levels
!      gradt = (hold_temps(n) - hold_temps(n-1))/dlp
!      tmp_t = hold_temps(n-1) + gradt *dlpt

!  Find layer with cloud bottom (clouds-as-layers only)
!   interpolate to cloud - bottom pressure

      if ( cloud_method .eq. 2 ) then

        loop = .true.
        n = 0
        do while (loop)
          n = n + 1
          if (     alternative.and.(cbot_height   .gt. hold_heights(n))) loop = .false.
          if (.not.alternative.and.(cbot_pressure .lt. hold_press(n))  ) loop = .false.
        enddo
        nlayers_b = n

        dlp   = hold_logp(n) - hold_logp(n-1)
        gradh = (hold_heights(n) - hold_heights(n-1))/dlp
        if ( alternative ) then
           hgt_b = cbot_height
           dlpb  = ( hgt_b - hold_heights(n-1) ) / gradh
           log_pbot = dlpb + hold_logp(n-1)
           prs_b = dexp(log_pbot)
           cbot_pressure = prs_b
        else
           log_pbot = dlog(cbot_pressure)
           dlpb = log_pbot     - hold_logp(n-1)
           hgt_b = hold_heights(n-1) + gradh *dlpb
           prs_b = cbot_pressure
           cbot_height = hgt_b
        endif

        tmp_b = hold_temps(n)
        shp_b = hold_shape(n)
        fac_b = dlpb / dlp
        ad_b  = fac_b * hold_airdens(n)
        add_b = fac_b * hold_dsair(n)
        do g = 1, ngases
          if ( do_gases(g) ) then
            pr_b(g)  = fac_b * hold_pr(n,g)
            prd_b(g) = fac_b * hold_prderiv(n,g)
          endif
        enddo

!  These two lines are for temperatures on levels
!        gradt = (hold_temps(n) - hold_temps(n-1))/dlp
!        tmp_b = hold_temps(n-1) + gradt *dlpb

!  Assign variables for the parceling

         diffh = cloud_opdep / ( hgt_t - hgt_b ) 
         
      endif

!  debug output

!      write(*,*)nlayers_t, hgt_t,tmp_t, o3_t, ctop_pressure
!      if(do_cloud_layering)write(*,*)&
!            nlayers_b, hgt_b,tmp_b, o3_b, cbot_pressure
!      write(*,*)' '
!      pause

!  Clouds as reflecting boundaries
!  -------------------------------

!  skip this if not appropriate

      if ( cloud_method .eq. 2 ) go to 566

!  total grid to cloud top

      nlayers_cloud = nlayers_t
      heights(0) = hold_heights(0)
      press(0)   = hold_press(0)
      do n = 1, nlayers_t - 1
        heights(n)    = hold_heights(n)
        temps(n)      = hold_temps(n)
        shapegrid(n)  = hold_shape(n)
        press(n)      = hold_press(n)
        airdens(n)    = hold_airdens(n)
        ds_airdens(n) = hold_dsair(n)
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = hold_pr(n,g)
            dc_profiles(n,g) = hold_prderiv(n,g)
          endif
        enddo
      enddo
      n = nlayers_t
      heights(n)    = hgt_t
      temps(n)      = tmp_t
      shapegrid(n)  = shp_t
      press(n)      = prs_t
      airdens(n)    = ad_t
      ds_airdens(n) = add_t
      do g = 1, ngases
        if ( do_gases(g) ) then
          profiles(n,g)    = pr_t(g)
          dc_profiles(n,g) = prd_t(g)
        endif
      enddo
      nlayers_total = nlayers_t

!  total grid: remaining layers to surface. Clearsky

      if ( do_clearsky ) then
        nc  = nlayers_t
        nc1 = nc + 1
        heights(nc1)    = hold_heights(nc)
        temps(nc1)      = hold_temps(nc)
        shapegrid(nc1)  = hold_shape(nc)
        press(nc1)      = hold_press(nc)
        airdens(nc1)    = hold_airdens(nc) - ad_t
        ds_airdens(nc1) = hold_dsair(nc) - add_t
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc1,g)    = hold_pr(nc,g) - pr_t(g)
            dc_profiles(nc1,g) = hold_prderiv(nc,g) - prd_t(g)
          endif
        enddo
        do n = nc1, nlayers_s
          n1 = n + 1
          heights(n1)    = hold_heights(n)
          temps(n1)      = hold_temps(n)
          shapegrid(n1)  = hold_shape(n)
          press(n1)      = hold_press(n)
          airdens(n1)    = hold_airdens(n)
          ds_airdens(n1) = hold_dsair(n)
         do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n1,g)    = hold_pr(n,g)
              dc_profiles(n1,g) = hold_prderiv(n,g)
            endif
          enddo
        enddo
        nlayers_total = nlayers_s + 1
      endif

!  finish

      return

!  Continuation point for clouds as layers

 566  continue

!  clouds as layers grid
!  ---------------------

!  Assign layers down to level immediately above cloud top

      heights(0) = hold_heights(0)
      press(0)   = hold_press(0)
      do n = 1, nlayers_t - 1
        heights(n)    = hold_heights(n)
        temps(n)      = hold_temps(n)
        shapegrid(n)  = hold_shape(n)
        press(n)      = hold_press(n)
        airdens(n)    = hold_airdens(n)
        ds_airdens(n) = hold_dsair(n)
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = hold_pr(n,g)
            dc_profiles(n,g) = hold_prderiv(n,g)
          endif
        enddo
      enddo

!  Assign reduced layer above cloud top

      n = nlayers_t
      heights(n)    = hgt_t
      temps(n)      = tmp_t
      shapegrid(n)  = shp_t
      press(n)      = prs_t
      airdens(n)    = ad_t
      ds_airdens(n) = add_t
      do g = 1, ngases
        if ( do_gases(g) ) then
          profiles(n,g)    = pr_t(g)
          dc_profiles(n,g) = prd_t(g)
        endif
      enddo

!  If the cloud is within one grid layer

      if (nlayers_b .eq. nlayers_t ) then

!....assign layer to bottom cloud level. THE CLOUD LAYER.

        n = n + 1
        heights(n)    = hgt_b
        temps(n)      = tmp_b
        shapegrid(n)  = shp_b
        press(n)      = prs_b
        airdens(n)    = ad_b - ad_t
        ds_airdens(n) = add_b - add_t
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = pr_b(g)  - pr_t(g)
            dc_profiles(n,g) = prd_b(g) - prd_t(g)
          endif
        enddo
        if ( .not. do_clearsky ) then
          cloudflag(n) = .true.
          cparcel(n)   = cloud_opdep
        endif

!....additional layers between cloud-bottom and surface

        nc  = nlayers_b
        nc2 = nc + 2
        heights(nc2)    = hold_heights(nc)
        temps(nc2)      = hold_temps(nc)
        shapegrid(nc2)  = hold_shape(nc)
        press(nc2)      = hold_press(nc)
        airdens(nc2)    = hold_airdens(nc) - ad_b
        ds_airdens(nc2) = hold_dsair(nc)   - add_b
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc2,g)    = hold_pr(nc,g) - pr_b(g)
            dc_profiles(nc2,g) = hold_prderiv(nc,g) - prd_b(g)
          endif
        enddo
        do n = nc+1, nlayers_s
          n2 = n + 2
          heights(n2)    = hold_heights(n)
          shapegrid(n2)  = hold_shape(n)
          temps(n2)      = hold_temps(n)
          press(n2)      = hold_press(n)
          airdens(n2)    = hold_airdens(n)
          ds_airdens(n2) = hold_dsair(n)
          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n2,g)    = hold_pr(n,g)
              dc_profiles(n2,g) = hold_prderiv(n,g)
            endif
          enddo
        enddo
        nlayers_total = nlayers_s + 2

!  If the cloud straddles a number of layers

      else

!.....assign first layer below the cloud top

        nc  = nlayers_t
        nc1 = nc + 1
        heights(nc1)    = hold_heights(nc)
        temps(nc1)      = hold_temps(nc)
        shapegrid(nc1)  = hold_shape(nc)
        press(nc1)      = hold_press(nc)
        airdens(nc1)    = hold_airdens(nc) - ad_t
        ds_airdens(nc1) = hold_dsair(nc) - add_t
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc1,g)    = hold_pr(nc,g) - pr_t(g)
            dc_profiles(nc1,g) = hold_prderiv(nc,g) - prd_t(g)
          endif
        enddo
        if ( .not. do_clearsky ) then
          cloudflag(nc1) = .true.
          cparcel(nc1)   = ( hgt_t - hold_heights(nc) ) * diffh
        endif

!....assign intervening layers to level immediately above cloud bottom

        do n = nc1, nlayers_b - 1
          n1 = n + 1
          heights(n1)    = hold_heights(n)
          shapegrid(n1)  = hold_shape(n)
          temps(n1)      = hold_temps(n)
          press(n1)      = hold_press(n)
          airdens(n1)    = hold_airdens(n)
          ds_airdens(n1) = hold_dsair(n)
          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n1,g)    = hold_pr(n,g)
              dc_profiles(n1,g) = hold_prderiv(n,g)
            endif
          enddo
          if ( .not. do_clearsky ) then
            cloudflag(n1) = .true.
            cparcel(n1)   = (hold_heights(n-1)-hold_heights(n)) * diffh
          endif
        enddo

!....assign final cloud layer to cloud bottom level

        n = nlayers_b + 1
        heights(n)    = hgt_b
        temps(n)      = tmp_b
        shapegrid(n)  = shp_b
        press(n)      = prs_b
        airdens(n)    = ad_b
        ds_airdens(n) = add_b
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = pr_b(g)
            dc_profiles(n,g) = prd_b(g)
          endif
        enddo
        if ( .not. do_clearsky ) then
          cloudflag(n) = .true.
          cparcel(n)   = ( hold_heights(nlayers_b-1) - hgt_b ) * diffh
        endif

!  ...more than one additional intervening layers to the surface

        nc = nlayers_b
        nc2 = nc + 2
        heights(nc2)    = hold_heights(nc)
        temps(nc2)      = hold_temps(nc)
        shapegrid(nc2)  = hold_shape(nc)
        press(nc2)      = hold_press(nc)
        airdens(nc2)    = hold_airdens(nc) - ad_b
        ds_airdens(nc2) = hold_dsair(nc)   - add_b
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc2,g)    = hold_pr(nc,g) - pr_b(g)
            dc_profiles(nc2,g) = hold_prderiv(nc,g) - prd_b(g)
          endif
        enddo
        do n = nc+1, nlayers_s
          n2 = n + 2
          heights(n2)    = hold_heights(n)
          temps(n2)      = hold_temps(n)
          shapegrid(n2)  = hold_shape(n)
          press(n2)      = hold_press(n)
          airdens(n2)    = hold_airdens(n)
          ds_airdens(n2) = hold_dsair(n)
          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n2,g)    = hold_pr(n,g)
              dc_profiles(n2,g) = hold_prderiv(n,g)
            endif
          enddo
        enddo
        nlayers_total = nlayers_s + 2

      endif

!  cloud grid

      nlayers_cloud = nlayers_total

!      write(*,*) nlayers_total,nlayers_s
!      do n = 90,nlayers_total
!        write(*,*)n, press(n),heights(n)
!      enddo
!      pause

!  Finish

      return
      end subroutine PTH_REGRIDDING_ALT_NOAER


!

      subroutine PTH_REGRIDDING_ALT_GENERAL                        &
       ( maxlayers, maxgases, maxaerwfs,                         & ! PURE INPUTS
         do_gases, do_aerosols, do_aer_Jacobians,                & ! PURE INPUTS
         ngases, surface_press,  n_aerosol_wfs,                  & ! PURE INPUTS
         alternative, cloud_method, do_clearsky, cloud_opdep,    & ! PURE INPUTS
         ctop_height, cbot_height, ctop_pressure, cbot_pressure, & ! MODIFIED/PURE INPUTS
         nlayers_total, press, heights, temps, shapegrid,        & ! MODIFIED INPUTS
         profiles, dc_profiles, airdens, ds_airdens,             & ! MODIFIED INPUTS
         aerflag, aertau_unscaled, L_aertau_unscaled,            & ! MODIFIED INPUTS
         nlayers_cloud, cloudflag, cparcel,                      & !  NEW OUTPUTS
         do_fail, message )                                        !  STATUS

      implicit none

!  dimensioning

      integer, intent(in)          :: maxlayers, maxgases, maxaerwfs

!  gas control

      integer, intent(in)          :: ngases
      logical, intent(in)          :: do_gases ( maxgases)

!  Aerosol control

      logical, intent(in)          :: do_aerosols
      logical, intent(in)          :: do_aer_Jacobians
      integer, intent(in)          :: n_aerosol_wfs

!  control inputs and cloud information

      logical, intent(in)          :: do_clearsky
      integer, intent(in)          :: cloud_method
      real(fpk)       , intent(in) :: cloud_opdep
      real(fpk)       , intent(in) :: surface_press
      logical, intent(in)          :: alternative

!  Cloud boundaries (in or out)

      real(fpk)       , intent(inout) :: ctop_height,   cbot_height
      real(fpk)       , intent(inout) :: ctop_pressure, cbot_pressure

!  Modified Grid outputs

      integer         , intent(inout) :: nlayers_total
      real(fpk)       , intent(inout) :: press     (0:maxlayers)
      real(fpk)       , intent(inout) :: temps     (  maxlayers)
      real(fpk)       , intent(inout) :: heights   (0:maxlayers)
      real(fpk)       , intent(inout) :: shapegrid (  maxlayers)

      real(fpk)       , intent(inout) :: profiles    (maxlayers,maxgases)
      real(fpk)       , intent(inout) :: dc_profiles (maxlayers,maxgases)
      real(fpk)       , intent(inout) :: airdens     (maxlayers)
      real(fpk)       , intent(inout) :: ds_airdens  (maxlayers)

      LOGICAL         , intent(inout) :: aerflag              (maxlayers)
      real(fpk)       , intent(inout) :: aertau_unscaled      (maxlayers)
      real(fpk)       , intent(inout) :: L_aertau_unscaled    (maxaerwfs,maxlayers)
 
!  Completely new outputs for cloud information

      integer         , intent(out) :: nlayers_cloud
      logical         , intent(out) :: cloudflag (maxlayers)
      real(fpk)       , intent(out) :: cparcel   (maxlayers)

!  status

      character*(*), intent(out) :: message
      logical      , intent(out) :: do_fail

!  local variables
!  ---------------

!  holding arrays 

      integer          maxlayers_hold, maxgases_hold, maxaerwfs_hold, nlayers_hold
      parameter        ( maxlayers_hold = 201 )
      parameter        ( maxgases_hold   = 2  )
      parameter        ( maxaerwfs_hold  = 15  )

      real(fpk)        hold_press  (0:maxlayers_hold)
      real(fpk)        hold_temps  (  maxlayers_hold)
      real(fpk)        hold_shape  (  maxlayers_hold)

      real(fpk)        hold_heights(0:maxlayers_hold)
      real(fpk)        hold_logp   (0:maxlayers_hold)
      real(fpk)        hold_airdens(  maxlayers_hold)
      real(fpk)        hold_dsair  (  maxlayers_hold)

      real(fpk)        hold_pr      (maxlayers_hold,maxgases_hold)
      real(fpk)        hold_prderiv (maxlayers_hold,maxgases_hold)

      logical          hold_flag     (maxlayers_hold)
      real(fpk)        hold_aer      (maxlayers_hold)
      real(fpk)        hold_aerderiv (maxaerwfs_hold,maxlayers_hold)

      real(fpk)        pr_t  (maxgases_hold), ad_t
      real(fpk)        prd_t (maxgases_hold), add_t
      real(fpk)        pr_b  (maxgases_hold), ad_b
      real(fpk)        prd_b (maxgases_hold), add_b

      real(fpk)        aer_t, aerd_t(maxaerwfs_hold)
      real(fpk)        aer_b, aerd_b(maxaerwfs_hold)

!  other help variables

      integer          n, nc, n1,n2, nc1, nc2, g, q
      integer          nlayers_s, nlayers_t, nlayers_b
      logical          loop, flag_b, flag_t
      real(fpk)        gradh, dlp, dlpb, dlpt, diffh, fac_t, fac_b
      real(fpk)        hgt_b, hgt_t, tmp_b, tmp_t, shp_t, shp_b
      real(fpk)        log_ptop, log_pbot, prs_b, prs_t

!  Initialise (important!)
!  ----------

      do_fail    = .false.
      message    = ' '

!  initialise nlayers_cloud and cloudflag and cparcel
!       (the only new variables to be set)

      nlayers_cloud = 0
      do n = 1, maxlayers
        cloudflag(n) = .false.
        cparcel(n)   = 0.0d0
      enddo

!  Nothing to do if Clearsky, Cloud_method = 0

      if ( cloud_method .eq. 0 ) return

!  Checking input
!  --------------

!  Check the cloud-top pressure is not out of bounds.

      if ( .not. alternative ) then
        if ( ctop_pressure .gt. surface_press ) then
          do_fail = .true.
          message = 'Abort: Ctop pressure > surface pressure'
          return
        endif
      endif

!  Check the cloud bottom pressure (clouds as layers only)

      if ( .not. alternative ) then
        if ( cloud_method.eq.2 ) then
          if ( cbot_pressure .gt. surface_press ) then
            do_fail = .true.
            message = 'Abort: Cbot pressure > surface pressure'
            return
          endif
        endif
      endif

!  Check the cloud-top height is not out of bounds.

      if ( alternative ) then
        if ( ctop_height .lt. heights(nlayers_total)  ) then
          do_fail = .true.
          message = 'Abort: Ctop height below surface height'
          return
        endif
      endif

!  Check the cloud bottom height (clouds as layers only)

      if ( alternative ) then
        if ( cloud_method.eq.2 ) then
          if ( cbot_height .lt. heights(nlayers_total) ) then
            do_fail = .true.
            message = 'Abort: Cbot height below surface height'
            return
          endif
        endif
      endif

!  copy to holding arrays
!  ----------------------

      nlayers_hold = nlayers_total
      nlayers_s    = nlayers_hold
      do n = 0, nlayers_hold
        hold_heights(n) = heights(n)
        hold_press(n)   = press(n)
        hold_logp(n)    = dlog(hold_press(n))
        if ( n .gt. 0 ) then
          hold_airdens(n) = airdens(n)
          hold_dsair(n)   = ds_airdens(n)
          hold_temps(n) = temps(n)
          hold_shape(n) = shapegrid(n)
        endif
      enddo

      do g = 1, ngases
        if ( do_gases(g) ) then
          do n = 1, nlayers_hold
            hold_pr(n,g)      = profiles(n,g)
            hold_prderiv(n,g) = dc_profiles(n,g)
          enddo
        endif
      enddo

      if ( do_aerosols ) then
         do n = 1, nlayers_hold
            hold_flag(n) = aerflag(n)
            hold_aer(n)  = aertau_unscaled(n)
            if ( do_aer_Jacobians ) then
               do q = 1, n_aerosol_wfs
                  hold_aerderiv(q,n) = L_aertau_unscaled(q,n)
               enddo
            endif
         enddo
      endif

!  Find the layers containing cloud levels
!  ---------------------------------------

!  Find layer containing cloud top
!    interpolate to cloud-top pressure

      loop = .true.
      n = 0
      do while (loop)
        n = n + 1
        if (     alternative.and.(ctop_height   .gt. hold_heights(n))) loop = .false.
        if (.not.alternative.and.(ctop_pressure .lt. hold_press(n))  ) loop = .false.
      enddo
      nlayers_t = n

      dlp   = hold_logp(n) - hold_logp(n-1)
      gradh = (hold_heights(n) - hold_heights(n-1))/dlp
      if ( alternative ) then
         hgt_t = ctop_height
         dlpt  = ( hgt_t - hold_heights(n-1) ) / gradh
         log_ptop = dlpt + hold_logp(n-1)
         prs_t = dexp(log_ptop)
         ctop_pressure = prs_t
      else
         log_ptop = dlog(ctop_pressure)
         dlpt = log_ptop     - hold_logp(n-1)
         hgt_t = hold_heights(n-1) + gradh *dlpt
         prs_t = ctop_pressure
         ctop_height = hgt_t
      endif

      tmp_t = hold_temps(n)
      shp_t = hold_shape(n)
      fac_t = dlpt / dlp
      ad_t  = fac_t * hold_airdens(n)
      add_t = fac_t * hold_dsair(n)

      do g = 1, ngases
        if ( do_gases(g) ) then
          pr_t(g)  = fac_t * hold_pr(n,g)
          prd_t(g) = fac_t * hold_prderiv(n,g)
        endif
      enddo

      if ( do_aerosols ) then
         flag_t = hold_flag(n)
         aer_t  = fac_t * hold_aer(n)
         if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              aerd_t(q) = fac_t * hold_aerderiv(q,n)
            enddo
         endif
      endif

!  These two lines are for temperatures on levels
!      gradt = (hold_temps(n) - hold_temps(n-1))/dlp
!      tmp_t = hold_temps(n-1) + gradt *dlpt

!  Find layer with cloud bottom (clouds-as-layers only)
!   interpolate to cloud - bottom pressure

      if ( cloud_method .eq. 2 ) then

        loop = .true.
        n = 0
        do while (loop)
          n = n + 1
          if (     alternative.and.(cbot_height   .gt. hold_heights(n))) loop = .false.
          if (.not.alternative.and.(cbot_pressure .lt. hold_press(n))  ) loop = .false.
        enddo
        nlayers_b = n

        dlp   = hold_logp(n) - hold_logp(n-1)
        gradh = (hold_heights(n) - hold_heights(n-1))/dlp
        if ( alternative ) then
           hgt_b = cbot_height
           dlpb  = ( hgt_b - hold_heights(n-1) ) / gradh
           log_pbot = dlpb + hold_logp(n-1)
           prs_b = dexp(log_pbot)
           cbot_pressure = prs_b
        else
           log_pbot = dlog(cbot_pressure)
           dlpb = log_pbot     - hold_logp(n-1)
           hgt_b = hold_heights(n-1) + gradh *dlpb
           prs_b = cbot_pressure
           cbot_height = hgt_b
        endif

        tmp_b = hold_temps(n)
        shp_b = hold_shape(n)
        fac_b = dlpb / dlp
        ad_b  = fac_b * hold_airdens(n)
        add_b = fac_b * hold_dsair(n)
        do g = 1, ngases
          if ( do_gases(g) ) then
            pr_b(g)  = fac_b * hold_pr(n,g)
            prd_b(g) = fac_b * hold_prderiv(n,g)
          endif
        enddo

        if ( do_aerosols ) then
           flag_b = hold_flag(n)
           aer_b  = fac_b * hold_aer(n)
           if ( do_aer_Jacobians ) then
              do q = 1, n_aerosol_wfs
                aerd_b(q) = fac_b * hold_aerderiv(q,n)
              enddo
           endif
        endif

!  These two lines are for temperatures on levels
!        gradt = (hold_temps(n) - hold_temps(n-1))/dlp
!        tmp_b = hold_temps(n-1) + gradt *dlpb

!  Assign variables for the parceling

         diffh = cloud_opdep / ( hgt_t - hgt_b ) 
         
      endif

!  debug output

!      write(*,*)nlayers_t, hgt_t,tmp_t, o3_t, ctop_pressure
!      if(do_cloud_layering)write(*,*)&
!            nlayers_b, hgt_b,tmp_b, o3_b, cbot_pressure
!      write(*,*)' '
!      pause

!  Clouds as reflecting boundaries
!  -------------------------------

!  skip this if not appropriate

      if ( cloud_method .eq. 2 ) go to 566

!  total grid to cloud top

      nlayers_cloud = nlayers_t
      heights(0) = hold_heights(0)
      press(0)   = hold_press(0)
      do n = 1, nlayers_t - 1
        heights(n)    = hold_heights(n)
        temps(n)      = hold_temps(n)
        shapegrid(n)  = hold_shape(n)
        press(n)      = hold_press(n)
        airdens(n)    = hold_airdens(n)
        ds_airdens(n) = hold_dsair(n)
        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = hold_pr(n,g)
            dc_profiles(n,g) = hold_prderiv(n,g)
          endif
        enddo
        if (do_aerosols ) then
          aerflag(n)         = hold_flag(n)
          aertau_unscaled(n) = hold_aer(n)
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,n) = hold_aerderiv(q,n)
            enddo
          endif
        endif
      enddo

!  reduced layer above cloud-top

      n = nlayers_t
      heights(n)    = hgt_t
      temps(n)      = tmp_t
      shapegrid(n)  = shp_t
      press(n)      = prs_t
      airdens(n)    = ad_t
      ds_airdens(n) = add_t

      do g = 1, ngases
        if ( do_gases(g) ) then
          profiles(n,g)    = pr_t(g)
          dc_profiles(n,g) = prd_t(g)
        endif
      enddo

      if (do_aerosols ) then
        aerflag(n)         = flag_t
        aertau_unscaled(n) = aer_t
        if ( do_aer_Jacobians ) then
          do q = 1, n_aerosol_wfs
            L_aertau_unscaled(q,n) = aerd_t(q)
          enddo
        endif
      endif

      nlayers_total = nlayers_t

!   Clear-sky situation, fill in rest of atmosphere

      if ( do_clearsky ) then

!  Partial layer immediately below cloud-top

        nc  = nlayers_t
        nc1 = nc + 1

        heights(nc1)    = hold_heights(nc)
        temps(nc1)      = hold_temps(nc)
        shapegrid(nc1)  = hold_shape(nc)
        press(nc1)      = hold_press(nc)
        airdens(nc1)    = hold_airdens(nc) - ad_t
        ds_airdens(nc1) = hold_dsair(nc) - add_t

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc1,g)    = hold_pr(nc,g) - pr_t(g)
            dc_profiles(nc1,g) = hold_prderiv(nc,g) - prd_t(g)
          endif
        enddo

        if (do_aerosols ) then
         aerflag(nc1)         = hold_flag(nc)
         aertau_unscaled(nc1) = hold_aer(nc) - aer_t
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,nc1) = hold_aerderiv(q,nc) - aerd_t(q)
            enddo
          endif
        endif

!  remaining layers to surface

        do n = nc1, nlayers_s

          n1 = n + 1
          heights(n1)    = hold_heights(n)
          temps(n1)      = hold_temps(n)
          shapegrid(n1)  = hold_shape(n)
          press(n1)      = hold_press(n)
          airdens(n1)    = hold_airdens(n)
          ds_airdens(n1) = hold_dsair(n)

          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n1,g)    = hold_pr(n,g)
              dc_profiles(n1,g) = hold_prderiv(n,g)
            endif
          enddo

          if ( do_aerosols ) then
            aerflag(n1)         = hold_flag(n)
            aertau_unscaled(n1) = hold_aer(n)
            if ( do_aer_Jacobians ) then
              do q = 1, n_aerosol_wfs
                L_aertau_unscaled(q,n1) =  hold_aerderiv(q,n) 
              enddo
            endif
          endif

        enddo

        nlayers_total = nlayers_s + 1
 
     endif

!  finish

      return

!  Continuation point for clouds as layers

 566  continue

!  clouds as layers grid
!  ---------------------

!  Assign layers down to level immediately above cloud top

      heights(0) = hold_heights(0)
      press(0)   = hold_press(0)

      do n = 1, nlayers_t - 1

        heights(n)    = hold_heights(n)
        temps(n)      = hold_temps(n)
        shapegrid(n)  = hold_shape(n)
        press(n)      = hold_press(n)
        airdens(n)    = hold_airdens(n)
        ds_airdens(n) = hold_dsair(n)

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = hold_pr(n,g)
            dc_profiles(n,g) = hold_prderiv(n,g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(n)         = hold_flag(n)
          aertau_unscaled(n) = hold_aer(n)
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,n) = hold_aerderiv(q,n)
            enddo
          endif
        endif

      enddo

!  Assign reduced layer above cloud top

      n = nlayers_t
      heights(n)    = hgt_t
      temps(n)      = tmp_t
      shapegrid(n)  = shp_t
      press(n)      = prs_t
      airdens(n)    = ad_t
      ds_airdens(n) = add_t

      do g = 1, ngases
        if ( do_gases(g) ) then
          profiles(n,g)    = pr_t(g)
          dc_profiles(n,g) = prd_t(g)
        endif
      enddo

      if (do_aerosols ) then
        aerflag(n)         = flag_t
        aertau_unscaled(n) = aer_t
        if ( do_aer_Jacobians ) then
          do q = 1, n_aerosol_wfs
            L_aertau_unscaled(q,n) = aerd_t(q)
          enddo
        endif
      endif

!  If the cloud is within one grid layer
!  =====================================

      if ( nlayers_b .eq. nlayers_t ) then

!....assign layer to bottom cloud level. THE CLOUD LAYER.

        n = n + 1
        heights(n)    = hgt_b
        temps(n)      = tmp_b
        shapegrid(n)  = shp_b
        press(n)      = prs_b
        airdens(n)    = ad_b - ad_t
        ds_airdens(n) = add_b - add_t

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = pr_b(g)  - pr_t(g)
            dc_profiles(n,g) = prd_b(g) - prd_t(g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(n)         = flag_t
          aertau_unscaled(n) = aer_b - aer_t
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,n) = aerd_b(q) - aerd_t(q)
            enddo
          endif
        endif

        if ( .not. do_clearsky ) then
          cloudflag(n) = .true.
          cparcel(n)   = cloud_opdep
        endif

!....additional layer between cloud-bottom and next lowest grid level

        nc  = nlayers_b
        nc2 = nc + 2
        heights(nc2)    = hold_heights(nc)
        temps(nc2)      = hold_temps(nc)
        shapegrid(nc2)  = hold_shape(nc)
        press(nc2)      = hold_press(nc)
        airdens(nc2)    = hold_airdens(nc) - ad_b
        ds_airdens(nc2) = hold_dsair(nc)   - add_b

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc2,g)    = hold_pr(nc,g) - pr_b(g)
            dc_profiles(nc2,g) = hold_prderiv(nc,g) - prd_b(g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(nc2)         = hold_flag(nc)
          aertau_unscaled(nc2) = hold_aer(nc) - aer_b
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,nc2) = hold_aerderiv(q,nc) - aerd_b(q)
            enddo
          endif
        endif

!....additional layers to surface, follwing original grid

        do n = nc+1, nlayers_s

          n2 = n + 2
          heights(n2)    = hold_heights(n)
          shapegrid(n2)  = hold_shape(n)
          temps(n2)      = hold_temps(n)
          press(n2)      = hold_press(n)
          airdens(n2)    = hold_airdens(n)
          ds_airdens(n2) = hold_dsair(n)

          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n2,g)    = hold_pr(n,g)
              dc_profiles(n2,g) = hold_prderiv(n,g)
            endif
          enddo

          if ( do_aerosols ) then
            aerflag(n2)         = hold_flag(n)
            aertau_unscaled(n2) = hold_aer(n)
            if ( do_aer_Jacobians ) then
              do q = 1, n_aerosol_wfs
                L_aertau_unscaled(q,n2) =  hold_aerderiv(q,n) 
              enddo
            endif
          endif

        enddo

        nlayers_total = nlayers_s + 2

!  If the cloud straddles a number of layers
!  =========================================

      else

!.....assign first layer below the cloud top

        nc  = nlayers_t
        nc1 = nc + 1
        heights(nc1)    = hold_heights(nc)
        temps(nc1)      = hold_temps(nc)
        shapegrid(nc1)  = hold_shape(nc)
        press(nc1)      = hold_press(nc)
        airdens(nc1)    = hold_airdens(nc) - ad_t
        ds_airdens(nc1) = hold_dsair(nc) - add_t

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc1,g)    = hold_pr(nc,g) - pr_t(g)
            dc_profiles(nc1,g) = hold_prderiv(nc,g) - prd_t(g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(nc1)         = hold_flag(nc)
          aertau_unscaled(nc1) = hold_aer(nc) - aer_t
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,nc1) = hold_aerderiv(q,nc) - aerd_t(q)
            enddo
          endif
        endif

        if ( .not. do_clearsky ) then
          cloudflag(nc1) = .true.
          cparcel(nc1)   = ( hgt_t - hold_heights(nc) ) * diffh
        endif

!....assign intervening layers to level immediately above cloud bottom

        do n = nc1, nlayers_b - 1

          n1 = n + 1
          heights(n1)    = hold_heights(n)
          shapegrid(n1)  = hold_shape(n)
          temps(n1)      = hold_temps(n)
          press(n1)      = hold_press(n)
          airdens(n1)    = hold_airdens(n)
          ds_airdens(n1) = hold_dsair(n)

          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n1,g)    = hold_pr(n,g)
              dc_profiles(n1,g) = hold_prderiv(n,g)
            endif
          enddo

          if ( do_aerosols ) then
            aerflag(n1)         = hold_flag(n)
            aertau_unscaled(n1) = hold_aer(n)
            if ( do_aer_Jacobians ) then
              do q = 1, n_aerosol_wfs
                L_aertau_unscaled(q,n1) =  hold_aerderiv(q,n) 
              enddo
            endif
          endif

          if ( .not. do_clearsky ) then
            cloudflag(n1) = .true.
            cparcel(n1)   = (hold_heights(n-1)-hold_heights(n)) * diffh
          endif

        enddo

!....assign final cloud layer to cloud bottom level

        n = nlayers_b + 1
        heights(n)    = hgt_b
        temps(n)      = tmp_b
        shapegrid(n)  = shp_b
        press(n)      = prs_b
        airdens(n)    = ad_b
        ds_airdens(n) = add_b

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(n,g)    = pr_b(g)
            dc_profiles(n,g) = prd_b(g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(n)         = flag_b
          aertau_unscaled(n) = aer_b
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,n) = aerd_b(q)
            enddo
          endif
        endif

        if ( .not. do_clearsky ) then
          cloudflag(n) = .true.
          cparcel(n)   = ( hold_heights(nlayers_b-1) - hgt_b ) * diffh
        endif

!  ...additional intervening layer to the next grid level

        nc = nlayers_b
        nc2 = nc + 2
        heights(nc2)    = hold_heights(nc)
        temps(nc2)      = hold_temps(nc)
        shapegrid(nc2)  = hold_shape(nc)
        press(nc2)      = hold_press(nc)
        airdens(nc2)    = hold_airdens(nc) - ad_b
        ds_airdens(nc2) = hold_dsair(nc)   - add_b

        do g = 1, ngases
          if ( do_gases(g) ) then
            profiles(nc2,g)    = hold_pr(nc,g) - pr_b(g)
            dc_profiles(nc2,g) = hold_prderiv(nc,g) - prd_b(g)
          endif
        enddo

        if (do_aerosols ) then
          aerflag(nc2)         = hold_flag(nc)
          aertau_unscaled(nc2) = hold_aer(nc) - aer_b
          if ( do_aer_Jacobians ) then
            do q = 1, n_aerosol_wfs
              L_aertau_unscaled(q,nc2) = hold_aerderiv(q,nc) - aerd_b(q)
            enddo
          endif
        endif

!  ...additional gridded layers to surface

        do n = nc+1, nlayers_s

          n2 = n + 2
          heights(n2)    = hold_heights(n)
          temps(n2)      = hold_temps(n)
          shapegrid(n2)  = hold_shape(n)
          press(n2)      = hold_press(n)
          airdens(n2)    = hold_airdens(n)
          ds_airdens(n2) = hold_dsair(n)

          do g = 1, ngases
            if ( do_gases(g) ) then
              profiles(n2,g)    = hold_pr(n,g)
              dc_profiles(n2,g) = hold_prderiv(n,g)
            endif
          enddo

          if ( do_aerosols ) then
            aerflag(n2)         = hold_flag(n)
            aertau_unscaled(n2) = hold_aer(n)
            if ( do_aer_Jacobians ) then
              do q = 1, n_aerosol_wfs
                L_aertau_unscaled(q,n2) =  hold_aerderiv(q,n) 
              enddo
            endif
          endif

        enddo

        nlayers_total = nlayers_s + 2

      endif

!  cloud grid

      nlayers_cloud = nlayers_total
       
!  Finish

      return
      end subroutine PTH_REGRIDDING_ALT_GENERAL

!  End module

End Module Regridding_Routines_m
