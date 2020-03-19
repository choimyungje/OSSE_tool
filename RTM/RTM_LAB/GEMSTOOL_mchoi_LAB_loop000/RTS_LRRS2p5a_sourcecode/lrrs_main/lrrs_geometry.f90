
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
! # Subroutines in this Module                                  #
! #                                                             #
! #            OUTGOING_SPHERGEOM_FINE_UP                       #
! #            OUTGOING_SPHERGEOM_FINE_DN                       #
! #            MULTI_OUTGOING_ADJUSTGEOM                        #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5. This is a New module (from V2.3)
!    (1) Three Geometry routines removed from original CORRECTION_1 module, put in here
!    (2) Bookkeeping improvements (use of "Only", clearer I/O specifications)

      MODULE lrrs_geometry_1_m

!  Parameter types

   USE LRRS_PARS_m, only : fpk, zero, one, two, pie, pi2, deg_to_rad, SDU

!private
public

contains

      SUBROUTINE OUTGOING_SPHERGEOM_FINE_UP_1 &
        ( MAXLAYERS, MAXFINE, DO_FINE, NLAYERS, NFINE,           & ! Inputs 
          HEIGHTS, ERADIUS, ALPHA_BOA, THETA_BOA, PHI_BOA,       & ! Inputs
          SUNPATHS,      RADII,      NTRAVERSE,      ALPHA_ALL,  & ! Output
          SUNPATHS_FINE, RADII_FINE, NTRAVERSE_FINE, ALPHA_FINE, & ! Output
          LOSPATHS, THETA_ALL, PHI_ALL, COSSCAT_UP,              & ! Output
          FAIL, MESSAGE )                                          ! Output

      IMPLICIT NONE

!  Completely stand-alone geometry routine for the outgoing correction
!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

!  This routine has the fine gridding treatment
!  Version 2.1. September 2007, Partial layer geometries added
!  Version 2.2. July 2009, This version without partials

!  Separate Downwelling and Upwelling routines, 15 February 2011

!  inputs

      INTEGER  , INTENT(IN) :: maxlayers, maxfine
      INTEGER  , INTENT(IN) :: nlayers, nfine
      LOGICAL  , INTENT(IN) :: do_fine
      REAL(FPK), INTENT(IN) :: eradius, heights (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_boa, theta_boa

!  modified inputs

      REAL(FPK), INTENT(INOUT) :: phi_boa

!  main outputs (geometry)

      INTEGER  , INTENT(OUT) :: ntraverse  (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: sunpaths   (0:maxlayers,maxlayers)
      REAL(FPK), INTENT(OUT) :: radii      (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: alpha_all  (0:maxlayers)

!  Fine level output (geometry)

      INTEGER  , INTENT(OUT) :: ntraverse_fine(maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: sunpaths_fine (maxlayers,maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: radii_fine    (maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: alpha_fine    (maxlayers,maxfine)

!  Other (incidental) geometrical output

      REAL(FPK), INTENT(OUT) :: lospaths   (maxlayers)
      REAL(FPK), INTENT(OUT) :: theta_all  (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: phi_all    (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: cosscat_up (0:maxlayers)

!  Status output

      LOGICAL          , INTENT(OUT) :: fail
      CHARACTER (LEN=*), INTENT(OUT) :: message

!  Local

      LOGICAL   :: direct_sun
      INTEGER   :: n, k, krad, n1
      REAL(FPK) :: ex, ey, ez, px, py, pz
      REAL(FPK) :: salpha_boa, calpha_boa, sphi_boa
      REAL(FPK) :: stheta_boa, ctheta_boa, cphi_boa
      REAL(FPK) :: ksi, cksi, sksi, xicum, tangr, fac
      REAL(FPK) :: ctheta, stheta, calpha, salpha, cphi
      REAL(FPK) :: b, sth0, th0, ks1, sth1, th1

!  Local arrays associated with fine grid output

      INTEGER   :: j

      INTEGER, parameter :: maxlocalfine = 20

      LOGICAL   :: direct_sunf(maxlocalfine)
      REAL(FPK) :: difz, dfine1, saf, xicum0, path
      REAL(FPK) :: thetaf(maxlocalfine), xicumf, difa
      REAL(FPK) :: cthetaf(maxlocalfine)
      REAL(FPK) :: sthetaf(maxlocalfine)
      REAL(FPK) :: ksif(maxlocalfine)

!  Initialise output

      fail = .false.
      message = ' '

!  Check range of inputs

      if ( alpha_boa.ge.90.0_fpk.or.alpha_boa.lt.zero ) then
        message = 'boa LOS angle outside range [0,90])'
        fail    = .true.
        return
      endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 10 January 2011------------------------- V. Natraj
!      if ( phi_boa.lt.zero )   phi_boa = - phi_boa
      if ( phi_boa.lt.zero )   phi_boa = 360.0_fpk + phi_boa
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr

!  IMPORTANT, DO NOT LIMIT THE AZIMUTH TO < 180 degrees

!      if ( phi_boa.gt.180.0_fpk ) phi_boa = 360.0_fpk - phi_boa    ! OLD
       if ( phi_boa.gt.360.0_fpk ) phi_boa = phi_boa - 360.0_fpk    ! NEW

!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if ( theta_boa.ge.90.0_fpk.or.theta_boa.lt.zero ) then
        message = 'boa SZA angle outside range [0,90])'
        fail    = .true.
        return
      endif
      if ( do_fine ) then
       if ( nfine.gt.maxlocalfine ) then
         message = 'local finelayer dimensioning insufficient'
         fail    = .true.
         return
       endif
      endif

!  Zero the sun paths
!  Initialize number of layers traversed  (nominal conditions)

      do n = 0, nlayers
        ntraverse(n) = n
        do k = 1, nlayers
         sunpaths(n,k) = zero
        enddo
      enddo

!  Zero the fine data paths

      if ( do_fine ) then
       dfine1 = dble(nfine) + 1
       do n = 1, nlayers
        do j = 1, nfine
         ntraverse_fine(n,j) = n
         do k = 1, nlayers
          sunpaths_fine(n,k,j) = zero
         enddo
        enddo
       enddo
      endif

!  Start at BOA

      alpha_all(nlayers) = alpha_boa * deg_to_rad
      theta_all(nlayers) = theta_boa * deg_to_rad
      phi_all(nlayers)   = phi_boa   * deg_to_rad

!  Cosine of scattering angle at boa

      salpha_boa = sin(alpha_all(nlayers))
      calpha_boa = cos(alpha_all(nlayers))
      stheta_boa = sin(theta_all(nlayers))
      ctheta_boa = cos(theta_all(nlayers))
      cphi_boa   = cos(phi_all(nlayers))
      sphi_boa   = sin(phi_all(nlayers))

      cosscat_up (nlayers) = - calpha_boa * ctheta_boa + &
                               salpha_boa * stheta_boa * cphi_boa

!  Radii
!  -----

!  Layer levels

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Fine levels

      if ( do_fine ) then
        do n = 1, nlayers
          difz = (radii(n-1)-radii(n))/dfine1
          do j = 1, nfine
            radii_fine(n,j) = radii(n) + difz * dble(j)
          enddo
        enddo
      endif

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit.
!    (This is the same as the regular pseudo-spherical )

      if ( salpha_boa.eq.zero ) then

!  WHOLE LAYER and FINE divisions
!  ------------------------------

!  Start layer loop, working upwards

        do n = nlayers,1,-1

!  Set main output.

          alpha_all(n-1)   = alpha_all(n)
          theta_all(n-1)   = theta_all(n)
          phi_all(n-1)     = phi_all(n)
          cosscat_up(n-1) = cosscat_up(n)
          lospaths(n) = radii(n-1)-radii(n)
          if ( do_fine ) then
            do j = 1, nfine
              alpha_fine(n,j) = zero
            enddo
          endif

!  Overhead sun

          if (stheta_boa.eq.zero ) then
            do k = n, 1, -1
              sunpaths(n,k) = radii(k-1)-radii(k)
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                do k = n - 1, 1, -1
                  sunpaths_fine(n,k,j) = radii(k-1)-radii(k)
                enddo
                sunpaths_fine(n,n,j) = radii(n-1)-radii_fine(n,j)
              enddo
            endif
          endif

!  Non-overhead sun
!  Main output of solar paths
!  Solar path distances for fine output

          if (stheta_boa.gt.zero ) then
            sth0 = stheta_boa
            th0  = theta_all(n)
            do k = n, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = asin(sth1)
              ks1  = th0-th1
              sunpaths(n,k) = sin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                sth0 = stheta_boa
                th0  = theta_all(n)
                sth1 = sth0*radii_fine(n,j)/radii(n-1)
                th1  = asin(sth1)
                ks1  = th0-th1
                sunpaths_fine(n,n,j) = sin(ks1)*radii_fine(n,j)/sth1
                sth0 = sth1
                th0  = th1
                do k = n-1, 1, -1
                  sth1 = sth0*radii(k)/radii(k-1)
                  th1  = asin(sth1)
                  ks1  = th0-th1
                  sunpaths_fine(n,k,j) = sin(ks1)*radii(k)/sth1
                  sth0 = sth1
                  th0  = th1
                enddo
              enddo
            endif
          endif

!  End main layer loop

        enddo

!  Return, as everything now done

        return

!  End regular pseudo-spherical clause, LOS is zero

      endif

!  Outgoing spehricity geometry
!  ============================

!  Define Unit solar vector at BOA

      ex = - stheta_boa * cphi_boa
      ey = - stheta_boa * sphi_boa
      ez = - ctheta_boa

!  Sun paths, boa geometry, always directly illuminated

      if ( stheta_boa.eq.zero ) then
        do k = nlayers, 1, -1
          sunpaths(nlayers,k) = radii(k-1)-radii(k)
        enddo
      else
        sth0 = stheta_boa
        th0  = theta_all(nlayers)
        do k = nlayers, 1, -1
          sth1 = sth0*radii(k)/radii(k-1)
          th1  = asin(sth1)
          ks1  = th0-th1
          sunpaths(nlayers,k) = sin(ks1)*radii(k)/sth1
          sth0 = sth1
          th0  = th1
        enddo
      endif

!  Check single illumination
!      --Not required, now we have the tangent point treatment
!      if (stheta_boa.gt.zero ) then
!        xicum = asin(radii(nlayers)*salpha_boa/radii(0))
!        xicum = alpha_all(nlayers)-xicum
!        px = - radii(0) * sin(xicum)
!        py = zero
!        pz =   radii(0) * cos(xicum)
!        b = ex*px + ey*py + ez*pz
!        ctheta = -b/radii(0)
!        if ( ctheta.le.zero ) then
!          write(*,*)'limit value = ',90.0_fpk-xicum/deg_to_rad
!        endif
!      endif

!  Initialise los cumulative angle

      xicum  = zero

!  Set TOA direct illumination flag

      direct_sun = .true.
      if ( do_fine ) then
        do j = 1, nfine
          direct_sunf(j) = .true.
        enddo
      endif

!  Start loop over positions (layer upper boundaries)

      do n = nlayers - 1, 0, -1

!  Next level up

        n1 = n + 1

!  Los angles at level boundaries

        salpha = radii(nlayers) * salpha_boa / radii(n)
        alpha_all(n)  = asin(salpha)
        calpha = cos(alpha_all(n))

!  Lospaths

        ksi = alpha_all(n1) - alpha_all(n)
        sksi = sin(ksi)
        cksi = cos(ksi)
        lospaths(n1) = sksi * radii(n1) / salpha
        xicum0 = xicum
        xicum  = xicum + ksi

!  Fine grid lospath output (angle and radius)
!    Locally save the earth-center angle ksif

        if ( do_fine ) then
          difa = (alpha_all(n1)-alpha_all(n))/dfine1
          do j = 1, nfine
            alpha_fine(n1,j) = alpha_all(n1) - difa * dble(j)
            saf = sin(alpha_fine(n1,j))
            radii_fine(n1,j) = salpha_boa * radii(nlayers) / saf
            ksif(j) = alpha_all(n1) - alpha_fine(n1,j)
          enddo
        endif

!  Sun angles for the Direct Nadir case

        if (stheta_boa.eq.zero ) then
         theta_all(n) = xicum
         ctheta = cos(theta_all(n))
         stheta = sqrt(one-ctheta*ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             thetaf(j)  = xicum0 + ksif(j)
             cthetaf(j) = cos(thetaf(j))
             sthetaf(j) = sqrt(one-ctheta*ctheta)
           enddo
         endif
        endif

!  Sun angles for the general case
!    Local save of angles, cosines, sines and  illumination flags

        if (stheta_boa.gt.zero ) then
         px = - radii(n) * sin(xicum)
         py = zero
         pz =   radii(n) * cos(xicum)
         b = ex*px + ey*py + ez*pz
         ctheta = -b/radii(n)
         direct_sun = (direct_sun.and.ctheta.ge.0.d0)
         stheta = sqrt(one-ctheta*ctheta)
         theta_all(n) = acos(ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             xicumf  = xicum0 + ksif(j)
             px = - radii_fine(n1,j) * sin(xicumf)
             py = zero
             pz =   radii_fine(n1,j) * cos(xicumf)
             b  = ex*px + ey*py + ez*pz
             cthetaf(j) = -b/radii_fine(n1,j)
             direct_sunf(j) = (direct_sunf(j).and.cthetaf(j).ge.0.d0)
             sthetaf(j) = sqrt(one-cthetaf(j)*cthetaf(j))
             thetaf(j)  = acos(cthetaf(j))
           enddo
         endif
        endif

!  Unit vector f2(i) perpendicular to OP but in plane of path
!  projection of f2(i) on solar path gives the relative azimuth at P
!        f2x = sin(xicum)
!        f2y = zero
!        f2z = cos(xicum)
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        if ( cphi.gt.one)  cphi = one
!        if ( cphi.lt.-one) cphi = -one
! ********************************************* Apparently not correct

!  Fix phi by using constancy of scatter angle
!  Only for the scattering up directions..................

        cosscat_up(n) = cosscat_up(n+1)
        if (stheta_boa.eq.zero ) then
          phi_all(n)     = phi_all(n+1)
        else
         cphi = (cosscat_up(n)+calpha*ctheta)/stheta/salpha
         if ( cphi.gt.one) cphi = one
         if ( cphi.lt.-one) cphi = -one
         phi_all(n)     = acos(cphi)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr

!  IF AZIMUTH > 180 degrees. Must ensure consistency, add this line.

         if ( phi_boa.gt.180.0_fpk ) phi_all(n) = pi2 - phi_all(n)

 !  B G CORRECTED 01 October 2010------------------------- RTS, R. Spurr
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        endif

!  Sun paths, Direct sun at layer top
!  ==================================

!   Means that the SZA at layer top is < 90.
!    ===> SZA < 90 for all fine points in layer beneath

        if ( direct_sun ) then

!  Work up from level n to TOA
!    Layer top calculation gets left out at TOA

         if ( n .gt. 0 ) then
          sth0 = stheta
          th0  = theta_all(n)
          do k = n, 1, -1
           sth1 = sth0*radii(k)/radii(k-1)
           th1  = asin(sth1)
           ks1  = th0-th1
           sunpaths(n,k) = sin(ks1)*radii(k)/sth1
           sth0 = sth1
           th0  = th1
          enddo
         endif

! DBG         write(*,*)'regular',n1,(sunpaths(n1,k),k=n1,1,-1)

!  Fine grid calculation is always required
!  ----------------------------------------

!  Start at the grid point on LOS path and work upwards,
!     first across partial layer to upper boundary, then continue
!     upwards to TOA on a whole layer basis. Sine rule.

         if ( do_fine ) then
          do j = 1, nfine
           sth0 = sthetaf(j)
           th0  = thetaf(j)
           sth1 = sth0*radii_fine(n1,j)/radii(n)
           th1  = asin(sth1)
           ks1  = th0-th1
           sunpaths_fine(n1,n1,j) = sin(ks1)*radii_fine(n1,j)/sth1
           sth0 = sth1
           th0  = th1
           do k = n, 1, -1
            sth1 = sth0*radii(k)/radii(k-1)
            th1  = asin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,k,j) = sin(ks1)*radii(k)/sth1
            sth0 = sth1
            th0  = th1
           enddo
!  DBG           write(*,*)'regular',n1,j,(sunpaths_fine(n1,k,j),k=n1,1,
          enddo

!  DBG         write(*,*)'regular',n,(sunpaths(n,k),k=n,1,-1)

         endif

!  Complete direct sun computations

        endif

!  Sun paths, Not direct sun , with tangent point
!  ==============================================

!  Although layer top has a tangent point, not all of the fine-grid
!   points will have a tangent point.

        if (.not.direct_sun ) then

!  First do the layer-top calculation.
!  -----------------------------------

!  TANGR = tangent point radius.

         tangr = stheta*radii(n)

!  ntraverse(n) is the number of layers traversed by ray.

         krad = nlayers
         do while (tangr.gt.radii(krad))
          krad = krad - 1
         enddo
         ntraverse(n) = krad + 1

!  Start at the TOA angles

         sth0 = tangr/radii(0)
         th0 = asin(sth0)

!  Work downwards from TOA (sine rule) to level immediately above
!  the tangent layer. Don't forget to double the path length for
!  any layers which are traversed twice.

         do k = 1, krad
          sth1 = radii(k-1)*sth0/radii(k)
          th1 = asin(sth1)
          ks1 = th1-th0
          fac = one
          if ( k.gt.n) fac = two
          sunpaths(n,k) = fac*sin(ks1)*radii(k-1)/sth1
          sth0 = sth1
          th0  = th1
         enddo

!  Tangent layer path length. Twice again. The following check is good.
!  check       write(*,*)tangr/dtan(th1),radii(krad)*cos(th1)

         sunpaths(n,krad+1)=two*radii(krad)*cos(th1)

!  DBG         write(*,*)'--tangent',n1,(sunpaths(n1,k),k=ntraverse(n1),1

!  Fine layer development (slightly different)
!  ---------------------

         if ( do_fine ) then
          do j = 1, nfine

!  If there is no tangent point, repeat calculation above.

           if ( direct_sunf(j) ) then
            sth0 = sthetaf(j)
            th0  = thetaf(j)
            sth1 = sth0*radii_fine(n1,j)/radii(n)
            th1  = asin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,n1,j) = sin(ks1)*radii_fine(n1,j)/sth1
            sth0 = sth1
            th0  = th1
            do k = n, 1, -1
             sth1 = sth0*radii(k)/radii(k-1)
             th1  = asin(sth1)
             ks1  = th0-th1
             sunpaths_fine(n1,k,j) = sin(ks1)*radii(k)/sth1
             sth0 = sth1
             th0  = th1
            enddo

! DBG           write(*,*)'regular',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Fine grid Calculation with tangent point
!  ----------------------------------------

           else

!  Local tangent radius and number of layers traversed

            tangr = sthetaf(j)*radii_fine(n1,j)
            krad = nlayers
            do while (tangr.gt.radii(krad))
             krad = krad - 1
            enddo
            ntraverse_fine(n1,j) = krad + 1

!  Start again at TOA

            sth0 = tangr/radii(0)
            th0  = asin(sth0)

!  Work down to the level n

            do k = 1, n
             sth1 = radii(k-1)*sth0/radii(k)
             th1 = asin(sth1)
             ks1 = th1-th0
             sunpaths_fine(n1,k,j) = sin(ks1)*radii(k-1)/sth1
             sth0 = sth1
             th0  = th1
            enddo

!  In layer below level n, Work down to level of fine grid radius
!     (single contribution)

            sth1 = radii(n)*sth0/radii_fine(n1,j)
            th1 = asin(sth1)
            ks1 = th1-th0
            path = sin(ks1)*radii(n)/sth1
            sth0 = sth1
            th0  = th1

!  In layer below level n, Complete the path down to the tangent point
!    (double contribution). Finish.

            if ( krad.eq.n ) then

              path = path + two*radii_fine(n1,j)*cos(th1)
              sunpaths_fine(n1,krad+1,j) = path

!  Check    write(*,*)'1',tangr/dtan(th1),radii_fine(n1,j)*cos(th1)

!  In layer below level n, Need to go down further.
!    (from now on, we need double contributions)
!     --- Continue the path to the bottom of this layer
!     --- Continue down by whole-layer steps until reach level above Tan
!     --- Complete the path down to the tangent point

            else
              sth1 = radii_fine(n1,j)*sth0/radii(n1)
              th1 = asin(sth1)
              ks1 = th1-th0
              path = path + two * sin(ks1)*radii_fine(n1,j)/sth1
              sunpaths_fine(n1,n1,j) = path
              sth0 = sth1
              th0  = th1
              do k = n1 + 1, krad
               sth1 = radii(k-1)*sth0/radii(k)
               th1 = asin(sth1)
               ks1 = th1-th0
               sunpaths_fine(n1,k,j) = sin(ks1)*radii(k-1)/sth1
               sth0 = sth1
               th0  = th1
              enddo
              sunpaths_fine(n1,krad+1,j)=two*radii(krad)*cos(th1)
!  Check 2    write(*,*)'2',tangr/dtan(th1),radii(krad)*cos(th1)
            endif

!  DBG           write(*,*)'tangent',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Complete tangent point clause

           endif

!  Complete fine-grid loop

          enddo

!  DBG        write(*,*)'tangent',n,(sunpaths(n,k),k=ntraverse(n),1,-1)

!  Complete tangent point calculation

         endif

        endif

!  End layer loop

      enddo

!  Finish

      return
      END SUBROUTINE OUTGOING_SPHERGEOM_FINE_UP_1

!

      SUBROUTINE OUTGOING_SPHERGEOM_FINE_DN_1 &
        ( MAXLAYERS, MAXFINE, DO_FINE, NLAYERS, NFINE,           & ! Input
          HEIGHTS, ERADIUS, ALPHA_BOA, THETA_BOA, PHI_BOA,       & ! Input
          SUNPATHS,      RADII,      NTRAVERSE,      ALPHA_ALL,  & ! Output
          SUNPATHS_FINE, RADII_FINE, NTRAVERSE_FINE, ALPHA_FINE, & ! Output
          LOSPATHS, THETA_ALL, PHI_ALL, COSSCAT_DN,              & ! Output
          FAIL, MESSAGE )                                          ! Output

      IMPLICIT NONE

!  Completely stand-alone geometry routine for the outgoing correction
!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

!  This routine has the fine gridding treatment
!  Version 2.1. September 2007, Partial layer geometries added
!  Version 2.2. July 2009, This version without partials

!  Separate Downwelling and Upwelling routines, 15 February 2011

!  Inputs

      INTEGER  , INTENT(IN) :: maxlayers, maxfine
      INTEGER  , INTENT(IN) :: nlayers, nfine
      LOGICAL  , INTENT(IN) :: do_fine
      REAL(FPK), INTENT(IN) :: eradius, heights (0:maxlayers)
      REAL(FPK), INTENT(IN) :: alpha_boa, theta_boa

!  Modified inputs

      REAL(FPK), INTENT(INOUT) :: phi_boa

!  Main outputs (geometry)

      INTEGER  , INTENT(OUT) :: ntraverse  (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: sunpaths   (0:maxlayers,maxlayers)
      REAL(FPK), INTENT(OUT) :: radii      (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: alpha_all  (0:maxlayers)

!  Fine level output (geometry)

      INTEGER  , INTENT(OUT) :: ntraverse_fine(maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: sunpaths_fine (maxlayers,maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: radii_fine    (maxlayers,maxfine)
      REAL(FPK), INTENT(OUT) :: alpha_fine    (maxlayers,maxfine)

!  Other (incidental) geometrical output

      REAL(FPK), INTENT(OUT) :: lospaths   (maxlayers)
      REAL(FPK), INTENT(OUT) :: theta_all  (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: phi_all    (0:maxlayers)
      REAL(FPK), INTENT(OUT) :: cosscat_dn (0:maxlayers)

!  Status output

      LOGICAL          , INTENT(OUT) :: fail
      CHARACTER (LEN=*), INTENT(OUT) :: message

!  Local

      LOGICAL   :: direct_sun
      INTEGER   :: n, k, krad, n1
      REAL(FPK) :: ex, ey, ez, px, py, pz
      REAL(FPK) :: salpha_boa, calpha_boa, sphi_boa
      REAL(FPK) :: stheta_boa, ctheta_boa, cphi_boa
      REAL(FPK) :: ksi, cksi, sksi, xicum, tangr, fac
      REAL(FPK) :: ctheta, stheta, calpha, salpha, cphi
      REAL(FPK) :: b, sth0, th0, ks1, sth1, th1

!  Local arrays associated with fine grid output

      INTEGER   :: j

      INTEGER, parameter :: maxlocalfine = 20

      LOGICAL   :: direct_sunf(maxlocalfine)
      REAL(FPK) :: difz, dfine1, saf, xicum0, path
      REAL(FPK) :: thetaf(maxlocalfine), xicumf, difa
      REAL(FPK) :: cthetaf(maxlocalfine)
      REAL(FPK) :: sthetaf(maxlocalfine)
      REAL(FPK) :: ksif(maxlocalfine)

!  Initialise output

      fail = .false.
      message = ' '

!  Check range of inputs

      if ( alpha_boa.ge.90.0_fpk.or.alpha_boa.lt.zero ) then
        message = 'boa LOS angle outside range [0,90])'
        fail    = .true.
        return
      endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 10 January 2011------------------------- V. Natraj
!      if ( phi_boa.lt.zero )   phi_boa = - phi_boa
      if ( phi_boa.lt.zero )   phi_boa = 360.0_fpk + phi_boa
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr

!  IMPORTANT, DO NOT LIMIT THE AZIMUTH TO < 180 degrees

!      if ( phi_boa.gt.180.0_fpk ) phi_boa = 360.0_fpk - phi_boa    ! OLD
       if ( phi_boa.gt.360.0_fpk ) phi_boa = phi_boa - 360.0_fpk    ! NEW

!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if ( theta_boa.ge.90.0_fpk.or.theta_boa.lt.zero ) then
        message = 'boa SZA angle outside range [0,90])'
        fail    = .true.
        return
      endif
      if ( do_fine ) then
       if ( nfine.gt.maxlocalfine ) then
         message = 'local finelayer dimensioning insufficient'
         fail    = .true.
         return
       endif
      endif

!  Zero the sun paths
!  Initialize number of layers traversed  (nominal conditions)

      do n = 0, nlayers
        ntraverse(n) = n
        do k = 1, nlayers
         sunpaths(n,k) = zero
        enddo
      enddo

!  Zero the fine data paths

      if ( do_fine ) then
       dfine1 = dble(nfine) + 1
       do n = 1, nlayers
        do j = 1, nfine
         ntraverse_fine(n,j) = n
         do k = 1, nlayers
          sunpaths_fine(n,k,j) = zero
         enddo
        enddo
       enddo
      endif

!  Start at BOA

      alpha_all(nlayers) = alpha_boa * deg_to_rad
      theta_all(nlayers) = theta_boa * deg_to_rad
      phi_all(nlayers)   = phi_boa   * deg_to_rad

!  Cosine of scattering angle at boa

      salpha_boa = sin(alpha_all(nlayers))
      calpha_boa = cos(alpha_all(nlayers))
      stheta_boa = sin(theta_all(nlayers))
      ctheta_boa = cos(theta_all(nlayers))
      cphi_boa   = cos(phi_all(nlayers))
      sphi_boa   = sin(phi_all(nlayers))

      cosscat_dn (nlayers) = + calpha_boa * ctheta_boa + &
                               salpha_boa * stheta_boa * cphi_boa

!  Radii
!  -----

!  Layer levels

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Fine levels

      if ( do_fine ) then
        do n = 1, nlayers
          difz = (radii(n-1)-radii(n))/dfine1
          do j = 1, nfine
            radii_fine(n,j) = radii(n) + difz * dble(j)
          enddo
        enddo
      endif

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit.
!    (This is the same as the regular pseudo-spherical )

      if ( salpha_boa.eq.zero ) then

!  WHOLE LAYER and FINE divisions
!  ------------------------------

!  Start layer loop, working upwards

        do n = nlayers,1,-1

!  Set main output.

          alpha_all(n-1)   = alpha_all(n)
          theta_all(n-1)   = theta_all(n)
          phi_all(n-1)     = phi_all(n)
          cosscat_dn(n-1) = cosscat_dn(n)
          lospaths(n) = radii(n-1)-radii(n)
          if ( do_fine ) then
            do j = 1, nfine
              alpha_fine(n,j) = zero
            enddo
          endif

!  Overhead sun

          if (stheta_boa.eq.zero ) then
            do k = n, 1, -1
              sunpaths(n,k) = radii(k-1)-radii(k)
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                do k = n - 1, 1, -1
                  sunpaths_fine(n,k,j) = radii(k-1)-radii(k)
                enddo
                sunpaths_fine(n,n,j) = radii(n-1)-radii_fine(n,j)
              enddo
            endif
          endif

!  Non-overhead sun
!  Main output of solar paths
!  Solar path distances for fine output

          if (stheta_boa.gt.zero ) then
            sth0 = stheta_boa
            th0  = theta_all(n)
            do k = n, 1, -1
              sth1 = sth0*radii(k)/radii(k-1)
              th1  = asin(sth1)
              ks1  = th0-th1
              sunpaths(n,k) = sin(ks1)*radii(k)/sth1
              sth0 = sth1
              th0  = th1
            enddo
            if ( do_fine ) then
              do j = 1, nfine
                sth0 = stheta_boa
                th0  = theta_all(n)
                sth1 = sth0*radii_fine(n,j)/radii(n-1)
                th1  = asin(sth1)
                ks1  = th0-th1
                sunpaths_fine(n,n,j) = sin(ks1)*radii_fine(n,j)/sth1
                sth0 = sth1
                th0  = th1
                do k = n-1, 1, -1
                  sth1 = sth0*radii(k)/radii(k-1)
                  th1  = asin(sth1)
                  ks1  = th0-th1
                  sunpaths_fine(n,k,j) = sin(ks1)*radii(k)/sth1
                  sth0 = sth1
                  th0  = th1
                enddo
              enddo
            endif
          endif

!  End main layer loop

        enddo

!  Return, as everything now done

        return

!  End regular pseudo-spherical clause, LOS is zero

      endif

!  Outgoing sphericity geometry
!  ============================

!  Define Unit solar vector at BOA
!     Note change of sign for ex here.

      ex = + stheta_boa * cphi_boa
      ey = - stheta_boa * sphi_boa
      ez = - ctheta_boa

!  Sun paths, boa geometry, always directly illuminated

      if ( stheta_boa.eq.zero ) then
        do k = nlayers, 1, -1
          sunpaths(nlayers,k) = radii(k-1)-radii(k)
        enddo
      else
        sth0 = stheta_boa
        th0  = theta_all(nlayers)
        do k = nlayers, 1, -1
          sth1 = sth0*radii(k)/radii(k-1)
          th1  = asin(sth1)
          ks1  = th0-th1
          sunpaths(nlayers,k) = sin(ks1)*radii(k)/sth1
          sth0 = sth1
          th0  = th1
        enddo
      endif

!  Check single illumination
!      --Not required, now we have the tangent point treatment
!      if (stheta_boa.gt.zero ) then
!        xicum = asin(radii(nlayers)*salpha_boa/radii(0))
!        xicum = alpha_all(nlayers)-xicum
!        px = - radii(0) * sin(xicum)
!        py = zero
!        pz =   radii(0) * cos(xicum)
!        b = ex*px + ey*py + ez*pz
!        ctheta = -b/radii(0)
!        if ( ctheta.le.zero ) then
!          write(*,*)'limit value = ',90.0_fpk-xicum/deg_to_rad
!        endif
!      endif

!  Initialise los cumulative angle

      xicum  = zero

!  Set TOA direct illumination flag

      direct_sun = .true.
      if ( do_fine ) then
        do j = 1, nfine
          direct_sunf(j) = .true.
        enddo
      endif

!  Start loop over positions (layer upper boundaries)

      do n = nlayers - 1, 0, -1

!  Next level up

        n1 = n + 1

!  Los angles at level boundaries

        salpha = radii(nlayers) * salpha_boa / radii(n)
        alpha_all(n)  = asin(salpha)
        calpha = cos(alpha_all(n))

!  Lospaths

        ksi = alpha_all(n1) - alpha_all(n)
        sksi = sin(ksi)
        cksi = cos(ksi)
        lospaths(n1) = sksi * radii(n1) / salpha
        xicum0 = xicum
        xicum  = xicum + ksi

!  Fine grid lospath output (angle and radius)
!    Locally save the earth-center angle ksif

        if ( do_fine ) then
          difa = (alpha_all(n1)-alpha_all(n))/dfine1
          do j = 1, nfine
            alpha_fine(n1,j) = alpha_all(n1) - difa * dble(j)
            saf = sin(alpha_fine(n1,j))
            radii_fine(n1,j) = salpha_boa * radii(nlayers) / saf
            ksif(j) = alpha_all(n1) - alpha_fine(n1,j)
          enddo
        endif

!  Sun angles for the Direct Nadir case

        if (stheta_boa.eq.zero ) then
         theta_all(n) = xicum
         ctheta = cos(theta_all(n))
         stheta = sqrt(one-ctheta*ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             thetaf(j)  = xicum0 + ksif(j)
             cthetaf(j) = cos(thetaf(j))
             sthetaf(j) = sqrt(one-ctheta*ctheta)
           enddo
         endif
        endif

!  Sun angles for the general case
!    Local save of angles, cosines, sines and  illumination flags

        if (stheta_boa.gt.zero ) then
         px = - radii(n) * sin(xicum)
         py = zero
         pz =   radii(n) * cos(xicum)
         b = ex*px + ey*py + ez*pz
         ctheta = -b/radii(n)
         direct_sun = (direct_sun.and.ctheta.ge.0.d0)
         stheta = sqrt(one-ctheta*ctheta)
         theta_all(n) = acos(ctheta)
         if ( do_fine ) then
           do j = 1, nfine
             xicumf  = xicum0 + ksif(j)
             px = - radii_fine(n1,j) * sin(xicumf)
             py = zero
             pz =   radii_fine(n1,j) * cos(xicumf)
             b  = ex*px + ey*py + ez*pz
             cthetaf(j) = -b/radii_fine(n1,j)
             direct_sunf(j) = (direct_sunf(j).and.cthetaf(j).ge.0.d0)
             sthetaf(j) = sqrt(one-cthetaf(j)*cthetaf(j))
             thetaf(j)  = acos(cthetaf(j))
           enddo
         endif
        endif

!  Unit vector f2(i) perpendicular to OP but in plane of path
!  projection of f2(i) on solar path gives the relative azimuth at P
!        f2x = sin(xicum)
!        f2y = zero
!        f2z = cos(xicum)
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        cphi = - (ex*f2x + ey*f2y + ez*f2z ) / stheta
!        if ( cphi.gt.one)  cphi = one
!        if ( cphi.lt.-one) cphi = -one
! ********************************************* Apparently not correct

!  Fix phi by using constancy of scatter angle
!  Only for the scattering up directions..................

        cosscat_dn(n) = cosscat_dn(n+1)
        if (stheta_boa.eq.zero ) then
          phi_all(n)     = phi_all(n+1)
        else
         cphi = (cosscat_dn(n)-calpha*ctheta)/stheta/salpha
         if ( cphi.gt.one) cphi = one
         if ( cphi.lt.-one) cphi = -one
         phi_all(n)     = acos(cphi)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr

!  IF AZIMUTH > 180 degrees. Must ensure consistency, add this line.

         if ( phi_boa.gt.180.0_fpk ) phi_all(n) = pi2 - phi_all(n)

 !  B G CORRECTED 01 October 2010------------------------- RTS, R. Spurr
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        endif

!  Sun paths, Direct sun at layer top
!  ==================================

!   Means that the SZA at layer top is < 90.
!    ===> SZA < 90 for all fine points in layer beneath

        if ( direct_sun ) then

!  Work up from level n to TOA
!    Layer top calculation gets left out at TOA

         if ( n .gt. 0 ) then
          sth0 = stheta
          th0  = theta_all(n)
          do k = n, 1, -1
           sth1 = sth0*radii(k)/radii(k-1)
           th1  = asin(sth1)
           ks1  = th0-th1
           sunpaths(n,k) = sin(ks1)*radii(k)/sth1
           sth0 = sth1
           th0  = th1
          enddo
         endif

! DBG         write(*,*)'regular',n1,(sunpaths(n1,k),k=n1,1,-1)

!  Fine grid calculation is always required
!  ----------------------------------------

!   Start at the grid point on LOS path and work upwards,
!     first across partial layer to upper boundary, then continue
!     upwards to TOA on a whole layer basis. Sine rule.

         if ( do_fine ) then
          do j = 1, nfine
           sth0 = sthetaf(j)
           th0  = thetaf(j)
           sth1 = sth0*radii_fine(n1,j)/radii(n)
           th1  = asin(sth1)
           ks1  = th0-th1
           sunpaths_fine(n1,n1,j) = sin(ks1)*radii_fine(n1,j)/sth1
           sth0 = sth1
           th0  = th1
           do k = n, 1, -1
            sth1 = sth0*radii(k)/radii(k-1)
            th1  = asin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,k,j) = sin(ks1)*radii(k)/sth1
            sth0 = sth1
            th0  = th1
           enddo
!  DBG           write(*,*)'regular',n1,j,(sunpaths_fine(n1,k,j),k=n1,1,
          enddo

!  DBG         write(*,*)'regular',n,(sunpaths(n,k),k=n,1,-1)

         endif

!  Complete direct sun computations

        endif

!  Sun paths, Not direct sun , with tangent point
!  ==============================================

!  Although layer top has a tangent point, not all of the fine-grid
!    points will have a tangent point.

        if (.not.direct_sun ) then

!  First do the layer-top calculation.
!  -----------------------------------

!  TANGR = tangent point radius.

         tangr = stheta*radii(n)

!  ntraverse(n) is the number of layers traversed by ray.

         krad = nlayers
         do while (tangr.gt.radii(krad))
          krad = krad - 1
         enddo
         ntraverse(n) = krad + 1

!  Start at the TOA angles

         sth0 = tangr/radii(0)
         th0 = asin(sth0)

!  Work downwards from TOA (sine rule) to level immediately above
!  the tangent layer. Don't forget to double the path length for
!  any layers which are traversed twice.

         do k = 1, krad
          sth1 = radii(k-1)*sth0/radii(k)
          th1 = asin(sth1)
          ks1 = th1-th0
          fac = one
          if ( k.gt.n) fac = two
          sunpaths(n,k) = fac*sin(ks1)*radii(k-1)/sth1
          sth0 = sth1
          th0  = th1
         enddo

!  Tangent layer path length. Twice again. The following check is good.
!  check       write(*,*)tangr/dtan(th1),radii(krad)*cos(th1)

         sunpaths(n,krad+1)=two*radii(krad)*cos(th1)

!  DBG         write(*,*)'--tangent',n1,(sunpaths(n1,k),k=ntraverse(n1),1

!  Fine layer development (slightly different)
!  ---------------------

         if ( do_fine ) then
          do j = 1, nfine

!  If there is no tangent point, repeat calculation above.

           if ( direct_sunf(j) ) then
            sth0 = sthetaf(j)
            th0  = thetaf(j)
            sth1 = sth0*radii_fine(n1,j)/radii(n)
            th1  = asin(sth1)
            ks1  = th0-th1
            sunpaths_fine(n1,n1,j) = sin(ks1)*radii_fine(n1,j)/sth1
            sth0 = sth1
            th0  = th1
            do k = n, 1, -1
             sth1 = sth0*radii(k)/radii(k-1)
             th1  = asin(sth1)
             ks1  = th0-th1
             sunpaths_fine(n1,k,j) = sin(ks1)*radii(k)/sth1
             sth0 = sth1
             th0  = th1
            enddo

!  DBG           write(*,*)'regular',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Fine grid Calculation with tangent point
!  ----------------------------------------

           else

!  Local tangent radius and number of layers traversed

            tangr = sthetaf(j)*radii_fine(n1,j)
            krad = nlayers
            do while (tangr.gt.radii(krad))
             krad = krad - 1
            enddo
            ntraverse_fine(n1,j) = krad + 1

!  Start again at TOA

            sth0 = tangr/radii(0)
            th0  = asin(sth0)

!  Work down to the level n

            do k = 1, n
             sth1 = radii(k-1)*sth0/radii(k)
             th1 = asin(sth1)
             ks1 = th1-th0
             sunpaths_fine(n1,k,j) = sin(ks1)*radii(k-1)/sth1
             sth0 = sth1
             th0  = th1
            enddo

!  In layer below level n, Work down to level of fine grid radius
!     (single contribution)

            sth1 = radii(n)*sth0/radii_fine(n1,j)
            th1 = asin(sth1)
            ks1 = th1-th0
            path = sin(ks1)*radii(n)/sth1
            sth0 = sth1
            th0  = th1

!  In layer below level n, Complete the path down to the tangent point
!    (double contribution). Finish.

            if ( krad.eq.n ) then

              path = path + two*radii_fine(n1,j)*cos(th1)
              sunpaths_fine(n1,krad+1,j) = path

!  Check    write(*,*)'1',tangr/dtan(th1),radii_fine(n1,j)*cos(th1)

!  In layer below level n, Need to go down further.
!    (from now on, we need double contributions)
!     --- Continue the path to the bottom of this layer
!     --- Continue down by whole-layer steps until reach level above Tan
!     --- Complete the path down to the tangent point

            else
              sth1 = radii_fine(n1,j)*sth0/radii(n1)
              th1 = asin(sth1)
              ks1 = th1-th0
              path = path + two * sin(ks1)*radii_fine(n1,j)/sth1
              sunpaths_fine(n1,n1,j) = path
              sth0 = sth1
              th0  = th1
              do k = n1 + 1, krad
               sth1 = radii(k-1)*sth0/radii(k)
               th1 = asin(sth1)
               ks1 = th1-th0
               sunpaths_fine(n1,k,j) = sin(ks1)*radii(k-1)/sth1
               sth0 = sth1
               th0  = th1
              enddo
              sunpaths_fine(n1,krad+1,j)=two*radii(krad)*cos(th1)
!  Check 2    write(*,*)'2',tangr/dtan(th1),radii(krad)*cos(th1)
            endif

!  DBG           write(*,*)'tangent',n1,j,
!     &       (sunpaths_fine(n1,k,j),k=ntraverse_fine(n1,j),1,-1)

!  Complete tangent point clause

           endif

!  Complete fine-grid loop

          enddo

!  DBG        write(*,*)'tangent',n,(sunpaths(n,k),k=ntraverse(n),1,-1)

!  Complete tangent point calculation

         endif

        endif

!  End layer loop

      enddo

!  Finish

      return
      END SUBROUTINE OUTGOING_SPHERGEOM_FINE_DN_1

!

      SUBROUTINE MULTI_OUTGOING_ADJUSTGEOM &
         ( MAX_VZA, MAX_AZM, N_VZA, N_AZM,       & ! Inputs
           HSURFACE, ERADIUS, DO_ADJUST_SURFACE, & ! Inputs
           ALPHA_BOA, THETA_BOA, PHI_BOA,        & ! Inputs
           ALPHA_SSA, THETA_SSA, PHI_SSA,        & ! Output
           FAIL, MESSAGE, TRACE )                  ! Output

!  Stand-alone geometry routine for adjusting the outgoing correction
!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height of new surface, earth radius.

       IMPLICIT NONE

!  Height grid here is artificial

!  Inputs

      INTEGER  , INTENT(IN) :: max_vza, max_azm
      INTEGER  , INTENT(IN) :: n_vza, n_azm
      REAL(FPK), INTENT(IN) :: eradius, hsurface
      LOGICAL  , INTENT(IN) :: do_adjust_surface
      REAL(FPK), INTENT(IN) :: alpha_boa (max_vza)
      REAL(FPK), INTENT(IN) :: theta_boa

!  Modified inputs

      REAL(FPK), INTENT(INOUT) :: phi_boa   (max_azm)

!  Outputs

      REAL(FPK), INTENT(OUT) :: alpha_ssa (max_vza)
      REAL(FPK), INTENT(OUT) :: theta_ssa (max_vza,max_azm)
      REAL(FPK), INTENT(OUT) :: phi_ssa   (max_vza,max_azm)
      LOGICAL  , INTENT(OUT) :: fail
      CHARACTER (LEN=*), INTENT(OUT) :: message, trace

!  Local

      INTEGER   :: j, k
      REAL(FPK) :: deg_to_rad, ex, ey, ez, px, py, pz, pie, pi2
      REAL(FPK) :: salpha_boa, calpha_boa, sphi_boa
      REAL(FPK) :: stheta_boa, ctheta_boa, cphi_boa
      REAL(FPK) :: ksi, cksi, sksi, xicum, cosscat_up
      REAL(FPK) :: phi_all, alpha_all, theta_all
      REAL(FPK) :: ctheta, stheta, calpha, salpha, cphi
      REAL(FPK) :: b,rssa
      CHARACTER (LEN=2) :: c2

!  Initialise output

      fail = .false.
      message = ' '
      trace   = ' '

!  Check range of inputs

      do j = 1, n_vza
       if ( alpha_boa(j).ge.90.0_fpk.or.alpha_boa(j).lt.zero ) then
        write(c2,'(I2)')J
        message = 'boa LOS angle outside range [0,90])'
        trace   = 'Change Boa Los angle, number '//C2
        fail    = .true.
        return
       endif
      enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr (1)
!  BUG CORRECTED 10 January 2011------------------------- V. Natraj     (2)

!  IMPORTANT, DO NOT LIMIT THE AZIMUTH TO < 180 degrees

      do k = 1, n_azm
!        if ( phi_boa(k).lt.zero   ) phi_boa(k) = - phi_boa(k)          (2)
        if ( phi_boa(k).lt.zero   ) phi_boa(k) = 360.0_fpk + phi_boa(k)
!        if ( phi_boa(k).gt.180.0_fpk ) phi_boa(k) = 360.0_fpk - phi_boa(k)  (1)
        if ( phi_boa(k).gt.360.0_fpk ) phi_boa(k) = phi_boa(k) - 360.0_fpk
      enddo

!  BUG CORRECTED 01 October 2010------------------------- RTS, R. Spurr
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if ( theta_boa.ge.90.0_fpk.or.theta_boa.lt.zero ) then
         message = 'boa SZA angle outside range [0,90])'
         trace   = 'Change Boa SZA angle'
         fail    = .true.
         return
      endif

!  No adjustment, just copy and exit

      if ( .not. do_adjust_surface ) then
         do j = 1, n_vza
            alpha_ssa(j)   = alpha_boa(j)
            do k = 1, n_azm
               theta_ssa(j,k) = theta_boa
               phi_ssa(j,k)   = phi_boa(k)
            enddo
         enddo
         return
      endif

!  Conversion

      pie = acos(-one)
      pi2 = two * pie
      deg_to_rad = acos(-one) / 180.0_fpk

!  Radius of surface

      rssa = hsurface + eradius

!  Start VZA loop

      do j = 1, n_vza

        alpha_all = alpha_boa(j) * deg_to_rad
        salpha_boa = sin(alpha_all)
        calpha_boa = cos(alpha_all)

!  Special case. Direct nadir viewing. Compute everything and Exit.
!    (This is the same as the regular pseudo-spherical )

        if ( salpha_boa.eq.zero ) then
           alpha_ssa(j)   = alpha_boa(j)
           do k = 1, n_azm
              theta_ssa(j,k) = theta_boa
              phi_ssa(j,k)   = phi_boa(k)
           enddo
           go to 567
        endif

!  Los angle

        salpha       = eradius * salpha_boa / rssa
        alpha_ssa(j) = asin(salpha)
        calpha       = cos(alpha_ssa(j))

!  Lospaths

        ksi  = alpha_all - alpha_ssa(j)
        sksi = sin(ksi)
        cksi = cos(ksi)
        xicum = ksi

!  Output angle in degrees

        alpha_ssa(j) = alpha_ssa(j) / deg_to_rad

!  Los vector

        px = - rssa * sin(xicum)
        py = zero
        pz =   rssa * cos(xicum)

        theta_all = theta_boa * deg_to_rad
        stheta_boa = sin(theta_all)
        ctheta_boa = cos(theta_all)

!  Start azimuth loop

        do k = 1, n_azm

          phi_all   = phi_boa(k)   * deg_to_rad
          cphi_boa  = cos(phi_all)
          sphi_boa  = sin(phi_all)

!  Define Unit solar vector

          ex = - stheta_boa * cphi_boa
          ey = - stheta_boa * sphi_boa
          ez = - ctheta_boa

!  Sun angle

          b = ex*px + ey*py + ez*pz
          ctheta = -b/rssa
          stheta = sqrt(one-ctheta*ctheta)
          theta_ssa(j,k) = acos(ctheta)/deg_to_rad
          if ( ctheta.lt.zero ) then
            write(c2,'(I2)')J
            message = 'LOS-path SZA angle outside range [0,90])'
            trace   = 'Check inputs for LOS angle '//c2
            fail    = .true.
            return
          endif

!  Scattering angle

          cosscat_up  = - calpha_boa * ctheta_boa + &
                          salpha_boa * stheta_boa * cphi_boa

!  Fix phi by using constancy of scatter angle

          if ( phi_boa(k).eq.180.0_fpk ) then
            phi_ssa(j,k) = phi_all
          else if ( phi_boa(k) .eq. zero ) then
            phi_ssa(j,k) = zero
          else
            cphi = (cosscat_up+calpha*ctheta)/stheta/salpha
            if ( cphi.gt.one) cphi = one
            if ( cphi.lt.-one) cphi = -one
            phi_ssa(j,k) = acos(cphi)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  BUG CORRECTED 05 October 2010----------------- V. Natraj, R. Spurr
!  IF AZIMUTH > 180 degrees. Must ensure consistency, add this line.
            if ( phi_boa(k).gt.180.0_fpk ) then
               phi_ssa(j,k)= pi2 - phi_ssa(j,k)
            endif
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

          endif
          phi_ssa(j,k) = phi_ssa(j,k) / deg_to_rad

!  End Azimuth loop

        enddo

!  Continuation point

 567    continue

!  End los loop

      enddo

!  Finish

      return
      END SUBROUTINE MULTI_OUTGOING_ADJUSTGEOM

!  Finish module

      END MODULE lrrs_geometry_1_m


