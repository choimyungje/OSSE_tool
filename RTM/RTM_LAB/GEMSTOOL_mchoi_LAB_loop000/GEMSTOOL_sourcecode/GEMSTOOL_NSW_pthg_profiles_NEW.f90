module GEMSTOOL_NSW_pthgas_profiles_new_m


!  11/16/18: MyungjeChoi
!  Line 210: nlevel 20 -->51-->101-->41
!  Line 161-162: PZERO, TZERO change    
!  Line 159: zsur change (meter unit)
!  Line 209: Lat change 34.201691

!  NSW-specific Type Structures. 12 August 2013

!  Uses data provided on 12 August 2013 from Y. Jung.

!  Major new code, 21 October 2013
!    Use fine-layering schemes


private
public :: GEMSTOOL_NSW_PTHGAS_PROFILES_NEW

contains

subroutine GEMSTOOL_NSW_PTHGAS_PROFILES_NEW &
           ( Maxlay, MaxFinelay, MaxGases, ngases, which_gases,                      & ! Input
             do_fine_divisions, do_fine_resolutions, Tshift, FDGas, FDLev, FDeps,    & ! Input
             PhysicsDataPath, PTH_profile_name, GAS_profile_names,                   & ! Input
             nlevels, level_heights, level_temps, level_press, level_logpress,       & ! Output (LEVELS)
             nlayers, aircolumns, daircolumns_dS, daircolumn_dP, Tshape, gascolumns, & ! Output (LAYERS). Air/Gas densities
             level_vmrs, levelgas, dLevelGas_dV, dLevelGas_dS, dLevelGas_dP,         & ! Output (LEVELS)
             nfinediv,  fineres,  heights_fine, temps_fine, press_fine, Gasfine,     & ! Fine-grid output
             dGasfine_dV,  dGasfine_dS, dGasfine_dP, dTfine_dS,                      & ! Fine-grid output
             fail, message )                                                           ! output

!             gascolumns, dgascolumns_dS, dgascolumn_dP, dgascolumns_dv,             &  ! Output (LAYERS). Gas densities

  implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  INPUTS
!  ======

!  dimensioning

   integer, intent(in)       :: Maxlay, MaxFinelay, MaxGases

!  Fine-layering control

   logical, intent(in)       :: do_fine_divisions, do_fine_resolutions

!  Gas control

   integer                                 , intent(in) :: NGASES
   character(LEN=4), dimension ( MaxGases ), intent(in) :: WHICH_GASES

!  temperature shift

   real(fpk), intent(in) :: Tshift

!  Indices for Finite Difference testing, Perturbation FDeps
!    FDGas Must be 1, 2, 3.... up to ngases
!    FDLev Must be 0, 1, 2, 3 ... up to nlayers

   integer  , intent(in) :: FDGas, FDLev
   real(fpk), intent(in) :: FDeps

!  File paths and file names

   character*(*), intent(in) :: PhysicsDataPath, PTH_profile_name, GAS_profile_names ( Maxgases)

!  OUTPUTS
!  =======

    !      Trace gas Layer column amounts and cross-sections
    !      Units should be consistent, e.g:

    !      gas X-sections   are in cm^2/mol  or inverse [DU]
    !      gas columns      are in mol/cm^2  or [DU]

    !      O2O2 X-sections  are in cm^5/mol^2    ! ALWAYS
    !      O2O2 columns     are in mol^2/cm^5    ! ALWAYS

!  Levels output

    integer                                      , intent(out)     :: NLEVELS
    REAL(fpk)   , DIMENSION( 0:maxlay )          , intent(out)     :: LEVEL_TEMPS
    REAL(fpk)   , DIMENSION( 0:maxlay )          , intent(out)     :: LEVEL_HEIGHTS
    REAL(fpk)   , DIMENSION( 0:maxlay )          , intent(out)     :: LEVEL_PRESS
    REAL(fpk)   , DIMENSION( 0:maxlay )          , intent(out)     :: LEVEL_LOGPRESS
    REAL(fpk)   , DIMENSION( 0:maxlay, maxgases ), intent(out)     :: LEVEL_VMRS
    REAL(fpk)   , DIMENSION( maxlay, maxgases )  , intent(out)     :: LEVELGAS

!  Fine-layer output

    integer     , DIMENSION( maxlay )            , intent(out)     :: NFINEDIV
    REAL(fpk)   , DIMENSION( maxlay )            , intent(out)     :: FINERES
    REAL(fpk)   , DIMENSION( maxlay, maxfinelay ), intent(out)     :: HEIGHTS_FINE
    REAL(fpk)   , DIMENSION( maxlay, maxfinelay ), intent(out)     :: TEMPS_FINE
    REAL(fpk)   , DIMENSION( maxlay, maxfinelay ), intent(out)     :: PRESS_FINE
    REAL(fpk)   , DIMENSION( maxlay, maxgases, maxfinelay ), intent(out) ::    Gasfine

!  Layer outputs, Air densities

    integer                                       , intent(out)    :: NLAYERS
    REAL(fpk)   , DIMENSION( maxlay )             , intent(out)    :: AIRCOLUMNS
    REAL(fpk)   , DIMENSION( maxlay )             , intent(out)    :: dAIRCOLUMNS_dS
    REAL(fpk)                                     , intent(out)    :: dAIRCOLUMN_dP

!  Layer outputs, Trace Gas densities

    REAL(fpk)   , DIMENSION( maxlay, maxgases )   , intent(out)    :: GASCOLUMNS

!  Derivative output for gas and temperatures

    REAL(fpk)   , DIMENSION( 0:maxlay, maxgases, 2 )            , intent(out)    :: dLevelGas_dV
    REAL(fpk)   , DIMENSION( 0:maxlay, maxgases, maxfinelay, 2 ), intent(out)    :: dGasfine_dV

    REAL(fpk)   , DIMENSION( 0:maxlay, maxgases, 2 )            , intent(out)    :: dLevelGas_dS
    REAL(fpk)   , DIMENSION( maxgases)                          , intent(out)    :: dLevelGas_dP
    REAL(fpk)   , DIMENSION( 0:maxlay, maxgases, maxfinelay, 2 ), intent(out)    :: dGasfine_dS
    REAL(fpk)   , DIMENSION(           maxgases, maxfinelay )   , intent(out)    :: dGasfine_dP

    REAL(fpk)   , DIMENSION( maxlay, maxfinelay, 2 )            , intent(out)    :: dTfine_dS


!    REAL(fpk)   , DIMENSION( maxlay, maxgases, 2 ), intent(out)    :: dGASCOLUMNS_dS
!    REAL(fpk)   , DIMENSION( maxlay, maxgases, 2 ), intent(out)    :: dGASCOLUMNS_dv
!    REAL(fpk)   , DIMENSION( maxgases )           , intent(out)    :: dGASCOLUMN_dP

!  Shape function (needed for later)

    REAL(fpk)   , DIMENSION ( 0:maxlay )          , intent(out)    :: TSHAPE

!  Excpetion handling

    logical      , intent(out) :: fail
    character*(*), intent(out) :: message

!  Local variables
!  ---------------

!  help

    integer            :: n, n1, g, k, ndum, j
    real(fpk)          :: constant, reftemp, refpress, dum, gasdum(6), zsur, avit, ccon
    character*256      :: file_pth, file_gas, file_fine

    real(fpk)  :: rhg1(maxgases), rhg2(maxgases), gas2(maxgases), gascol(maxgases), alpha, beta
    real(fpk)  :: rho1, rho2, air2, aircol, hght1, hght2, temp1, temp2, pres1, pres2, logp2, hdiff
    real(fpk)  :: diffh(maxlay), vmr1(maxgases), vmr2(maxgases)
    real(fpk)  :: drho1_dS(2), drho2_dS(2), dair2_dS(2), daircol_dS(2), dtemp1_dS(2), dtemp2_dS(2)
    real(fpk)  :: drho1_dP, drho2_dP, dair2_dP, daircol_dP

    real(fpk), PARAMETER :: ZERO = 0.0_fpk

!  Loschmidt's number (particles/cm3), STP values

    real(fpk)   , PARAMETER :: RHO_STANDARD = 2.68675D+19
    real(fpk)   , PARAMETER :: PZERO = 1015.0811D0 !1013.25D0
    real(fpk)   , PARAMETER :: TZERO = 289.92776D0 !!273.15D0

!  Zero the output
!  ===============

   fail = .false. ; message = ' '

   nlevels = 0 ; nlayers = 0

   level_heights = 0.0_fpk
   level_press   = 0.0_fpk
   level_temps   = 0.0_fpk
   level_logpress= 0.0_fpk
   level_vmrs    = 0.0_fpk
   levelgas      = 0.0_fpk

   aircolumns     = 0.0_fpk ; gascolumns     = 0.0_fpk
   daircolumns_dS = 0.0_fpk ; daircolumn_dP  = 0.0_fpk 

   dLevelGas_dS = 0.0_fpk
   dLevelGas_dV = 0.0_fpk
   dLevelGas_dP = 0.0_fpk

!  Fine level stuff

   heights_fine = zero
   temps_fine   = zero
   press_fine   = zero
   gasfine      = zero

   fineres  = zero
   nfinediv = 0   

   dGasfine_dV = zero
   dGasfine_dS = zero
   dGasfine_dP = zero
   dTfine_dS   = zero

!  Uniform T-shift profile shape

   Tshape = 1.0_fpk

!  PTH profiles
!  ============

!  Open PTH file and read the data.
!      Assume file = Initial_NSW_atmosphere_12aug13.dat, then we get VMRs at the same Pth levels/

!  Set the number of layers, check dimensions immediately]

   nlevels = 25 ; nlayers = nlevels - 1 
   !    nlevels = 20 ; nlayers = nlevels - 1 
   if ( nlayers .gt. maxlay ) then
      message = 'NLAYERS > Dimension MAXLAYERS ==> Increase VLIDORT parameter MAXLAYERS'
      fail = .true. ; return
   endif

!  Now read the file

   file_pth = adjustl(trim(PhysicsDataPath))//'PTH_PROFILES/'//adjustl(trim(Pth_profile_name))
   open(1,file=trim(file_pth),err=90,status='old')
   do k = 1, 12
      read(1,*)
   enddo
   do n = 1, nlevels
      n1 = n - 1 ; read(1,*)refpress, reftemp, dum,dum !,dum,dum,dum,dum
      level_press(n1) = refpress * 1.0d-02
      level_temps(n1) = reftemp + Tshift *  Tshape(n1)
   enddo

!  heights by Hydrostatic Eqn. (includes TOA)
!    For now, zsur = surface height in [m] ------Topography data set ????

   zsur = 0.347d0 !0.0d0
   Level_heights(nlayers) = zsur / 1000.0d0
   ccon = - 9.81d0 * 28.9d0 / 8314.0d0 * 500.0d0
   do n = nlayers, 1, -1
      avit = (1.0d0/level_temps(n-1))+(1.0d0/level_temps(n))
      level_heights(n-1) = level_heights(n) - log(Level_press(n)/Level_press(n-1))/avit/ccon
   enddo

! !  read from file

!    open(1,file='NSW_heightgrid_1.dat',status='old')
!    do n = 0, nlayers
!       read(1,*)n1,level_heights(n)
!    enddo
!    close(1)

!  Open GAS files and read the data

!    @@@@ 7/25/13. Assume for now that GAS profiles use same level scheme as PTH. UVN profile (Discover AQ)
!                  Assume everything is VMR (not ppmv)
!                  Assume only 5 gases available ---> NO BRO

!    @@@@ 8/12/13. Assume for now that GAS profiles use same level scheme as PTH. NSW profile (Initial 20-level)
!                  Assume everything is VMR (not ppmv)
!                  Assume 4 gases available ---> O2, H2O, CH4, CO2

   level_vmrs = 0.0d0
   do g = 1, ngases
      file_gas = adjustl(trim(PhysicsDataPath))//'GAS_PROFILES/'//adjustl(trim(Gas_profile_names(g)))
      open(1,file=trim(file_gas),err=91,status='old')
      do k = 1, 12
         read(1,*)
      enddo
      do n = 1, nlevels
            n1 = n - 1 ; read(1,*)dum,dum,(gasdum(k),k=1,2) ! 6 gas entries here ;---order is changed O2<->CO2
            if ( which_gases(g).eq.'H2O ')level_vmrs(n1,g) = gasdum(1)
            if ( which_gases(g).eq.'O2  ')level_vmrs(n1,g) = gasdum(2)
            ! if ( which_gases(g).eq.'CO2 ')level_vmrs(n1,g) = gasdum(3)
            ! if ( which_gases(g).eq.'CH4 ')level_vmrs(n1,g) = gasdum(4)
            if (FDgas.eq.g.and.FDLev.eq.n1) level_vmrs(n1,g) = level_vmrs(n1,g) * ( 1.0d0 + FDeps )
      enddo
      close(1)
   enddo

!  Number of layers = number of levels MINUS 1

   nlayers = nlevels - 1

!  Gas constant

   CONSTANT  = 1.0D+05 * 0.5d0 * RHO_STANDARD * TZERO / PZERO

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ NEW SECTION 21 OCTOBER 2013 @@@@@

!  height differentials

   do n = 1, nlayers
      diffh(n) = level_heights(n-1) - level_heights(n)
   enddo

!  Read the fine-resolution or fine-division data

   if ( do_fine_divisions ) then
      file_fine = adjustl(trim(PhysicsDataPath))//'PTH_PROFILES/Initial_NSW_atmosphere_FineDivisions.dat'
      open(1,file=trim(file_fine),err=90,status='old')
      do n = 1, nlayers
         read(1,*)ndum,nfinediv(n)
         fineres(n) = diffh(n) / real(nfinediv(n),fpk)
      enddo
   else if ( do_fine_resolutions ) then
      file_fine = adjustl(trim(PhysicsDataPath))//'PTH_PROFILES/Initial_NSW_atmosphere_FineResolutions.dat'
      open(1,file=trim(file_fine),err=90,status='old')
      do n = 1, nlayers
         read(1,*)ndum,fineres(n)
         nfinediv(n) = nint(diffh(n)/fineres(n)) - 1
      enddo
   endif

!  Find the fine level heights

   do n = 1, nlayers
      do j = 1, Nfinediv(n)
         n1 = n - 1
         heights_fine(n,j)  = Level_heights(n1) - fineres(n) * real(j)
      enddo
   enddo

!  Log pressures

   Level_logpress(0:nlayers) = log (level_press(0:nlayers))

!  Start layer loop

   do n = 1, nlayers

!  Upper level values of PTHG

      n1 = n - 1
      aircol = zero ; gascol(1:ngases) = zero ; daircol_dS = zero

      temp1 = level_temps(n1)
      pres1 = level_press(n1)
      hght1 = level_heights(n1)

      dTemp1_dS(1) = Tshape(n1) ; dTemp1_dS(2) = zero
      rho1  = constant * pres1 / temp1
      drho1_dS(1:2)  = -rho1 * dTemp1_dS(1:2) / temp1

      do g = 1, ngases
         vmr1(g) = level_vmrs(n1,g)
         rhg1(g) = rho1 * vmr1(g)
         Levelgas(n1,g) = rhg1(g)
         dLevelGas_dV(n1,g,1) = rho1
!         Atmos%dLevelGas_dV(n1,gasvmr_mask(v),2) = zero
         dLevelGas_dS(n1,g,1) = vmr1(g) * drho1_dS(1)
      enddo

!  Surface pressure

      if ( n.eq.nlayers ) drho1_dP = zero

!  Fine layer computations

      if ( Nfinediv(n) .gt. 0 ) then
         do j = 1, Nfinediv(n)

!  Fine-level values of PTHG

            hght2 = Heights_fine(n,j)
            alpha = - ( Level_heights(n) - Heights_fine(n,j) ) / diffh(n)
            beta  = 1.0_fpk - alpha
            temp2 = alpha * level_temps(n1)   + beta * level_temps(n)

            dTemp2_dS(1) = alpha * Tshape(n1) ; dTemp2_dS(2) = beta * Tshape(n)
 
            Temps_fine(n,j)  = temp2
            dTfine_dS(n,j,:) = dtemp2_dS(:)

            logp2        = alpha * Level_logpress (n1)   + beta * Level_logpress (n)
            pres2        = exp(logp2)
            press_fine(n,j)  = pres2

            rho2          = constant * pres2 / temp2 
            drho2_dS(1:2) = -rho2 * dTemp2_dS(1:2) / temp2

            do g = 1, ngases
               vmr2(g) = alpha * Level_vmrs(n1,g) + beta * Level_vmrs (n,g)
               rhg2(g) = rho2 * vmr2(g)
               Gasfine(n,g,j) = rhg2(g)
               dGasFine_dV(n,g,j,1) = rho2 * alpha
               dGasFine_dV(n,g,j,2) = rho2 * beta
               dGasFine_dS(n,g,j,:) = vmr2(g) * drho2_dS(:)
               if(n.eq.nlayers)dGasFine_dP(g,j) = vmr2(g) * drho2_dP
            enddo

!  surface pressure

            if ( n.eq.nlayers ) drho2_dP = constant * beta / temp2 

!  Compute intermediate columns

            hdiff      = hght1 - hght2
            air2       = hdiff * (rho1 + rho2)
            aircol     = aircol + air2

            dair2_dS(1:2)   = hdiff * (drho1_dS(1:2) + drho2_dS(1:2)) 
            daircol_dS(1:2) = daircol_dS(1:2) + dair2_dS(1:2)

            gas2(1:ngases)   = hdiff * (rhg1(1:ngases) + rhg2(1:ngases))
            gascol(1:ngases) = gascol(1:ngases) + gas2(1:ngases)

            if ( n.eq.nlayers ) then
                dair2_dP   = hdiff * (drho1_dP + drho2_dP)
                daircol_dP = daircol_dP + dair2_dP
            endif

 !  update

            hght1 = hght2 ; temp1 = temp2 ; pres1 = pres2 ; vmr1(1:ngases) = vmr2(1:ngases)
            rho1 = rho2   ; drho1_dS(:) = drho2_dS(:)
            rhg1(1:ngases) = rhg2(1:ngases)
            if ( n.eq.nlayers)drho1_dP = drho2_dP

!  End fine-level computations

         enddo
      endif

!  Lower Level values of PTHG

      temp2 = Level_temps(n)
      pres2 = Level_press(n)
      hght2 = Level_heights(n)

      dTemp2_dS(1) = zero ; dTemp2_dS(2) = Tshape(n)

      rho2        = constant * pres2 / temp2
      drho2_dS(:) = - rho2 * dTemp2_dS(:) / temp2

      if ( n.eq.nlayers )drho2_dP = constant / temp2

      do g = 1, ngases
         vmr2(g) = Level_vmrs(n,g)
         rhg2(g) = rho2 * vmr2(g)
         Levelgas(n,g) = rhg2(g)
         dLevelGas_dV(n,g,2) = rho2
         dLevelGas_dS(n,g,2) = vmr2(g) * drho2_dS(2)
         if (n.eq.nlayers)dLevelGas_dP(g) = vmr2(g) * drho2_dP
      enddo

!  Add remaining columns

      hdiff      = hght1 - hght2
      air2       = hdiff * (rho1 + rho2)
      dair2_dS   = hdiff * (drho1_dS + drho2_dS) 
      aircol     = aircol + air2
      daircol_dS = daircol_dS + dair2_dS
      if (n.eq.nlayers) then
         dair2_dP   = hdiff * drho2_dP
         daircol_dP = daircol_dP + dair2_dP
      endif

      gas2(1:ngases)   = hdiff * (rhg1(1:ngases) + rhg2(1:ngases))
      gascol(1:ngases) = gascol(1:ngases) + gas2(1:ngases) 

!  Final layer column values

      aircolumns(n)     = aircol
      daircolumns_dS(n) = daircol_dS(1) + daircol_dS(2)
      if (n.eq.nlayers) daircolumn_dP = daircol_dP
      gascolumns(n,1:ngases) = gascol(1:ngases)

!  End layer loop

   enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END NEW SECTION  @@@@@


!  normal Finish

   return

!  error finish

90 continue
   fail    = .true.
   message = 'PTH file not found' 
   return
91 continue
   fail    = .true.
   message = 'GAS file not found : '//which_gases(g) 
   return

!  End

end subroutine GEMSTOOL_NSW_PTHGAS_PROFILES_NEW

!  End module

end module GEMSTOOL_NSW_pthgas_profiles_new_m


