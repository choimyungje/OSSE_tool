module GEMSTOOL_NSW_ODSETTER_m

!  Subroutines to calculate trace-gas layer absorption optical depths, based on fine-layering scheme.

!    GEMSTOOL_NSW_ODSETTER        : Trace-gas absorption optical depths; NO DERIVATIVES
!    GEMSTOOL_NSW_ODSETTER_PLUS_1 : Trace-gas absorption optical depths; VMR derivatives alone
!    GEMSTOOL_NSW_ODSETTER_PLUS_2 : Trace-gas absorption optical depths; VMR, T-shift and Surf-P derivatives

!  R. Spurr, 21 October 2013

public

contains

subroutine GEMSTOOL_NSW_ODSETTER &
     ( maxlayers, maxfinelayers, maxwav, maxgases, & ! Input dimensions
       nlayers, nfinediv, ngases, do_gases, w,  &  ! Input control
       level_heights, levelgas, gas_xsecs,         & ! Input (LEVELS)
       heights_fine, Gas_fine,  Gasxsecs_fine,     & ! Fine-grid input
       odabs, odabstot )                             ! Output

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input arguments
!  ---------------

!  Dimensioning

   integer  , intent(in) :: maxlayers, maxfinelayers, maxwav, maxgases

!  Layer control

   INTEGER  , intent(in)  :: NLAYERS
   INTEGER  , intent(in)  :: NFINEDIV ( MAXLAYERS)


!  Gas control

   integer  , intent(in)  :: NGASES
   logical  , intent(in)  :: do_GASES( maxgases )

!  Spectral index

   integer  , intent(in)  :: w

!  Level Heights and Gas values + cross-sections

   real(fpk), intent(in)  :: LEVEL_HEIGHTS ( 0:MAXLAYERS) 
   real(fpk), intent(in)  :: LEVELGAS      ( 0:MAXLAYERS, maxgases )
   REAL(fpk), intent(in)  :: GAS_XSECS     ( MAXWAV, 0:maxlayers, maxgases ) 

!  Heights and Gas values, cross-sections: FINE-LEVEL values

   real(fpk), intent(in)  :: HEIGHTS_FINE  ( MAXLAYERS, maxfinelayers ) 
   real(fpk), intent(in)  :: GAS_FINE      ( MAXLAYERS, maxfinelayers, maxgases )
   REAL(fpk), intent(in)  :: GASXSECS_FINE ( MAXWAV, maxlayers, maxfinelayers, maxgases ) 

!  Output
!  ------

!  Gas absoprtion optical depths

   real(fpk), intent(out) :: odabs    ( MAXLAYERS, maxgases )
   real(fpk), intent(out) :: odabstot ( MAXLAYERS )

!  Local variables
!  ---------------

   integer   :: n, n1, j, g
   real(fpk) :: hght1, hght2, od1, od2, hdiff, tau

!  Initialize

   odabs    = 0.0_fpk
   odabstot = 0.0_fpk

!  Start layer and gas loops

   do n = 1, nlayers
      n1 = n - 1
      do g = 1, ngases
         if ( do_gases(g) ) then

!  Initialize

            tau = 0.0_fpk

!  Upper boundary values

            hght1 = level_heights(n1)
            od1   = gas_xsecs(w,n1,g) * LevelGas(n1,g)

!  Loop over fine layer subdivisions

            do j = 1, nfinediv(n)
               hght2 = heights_fine(n,j)
               od2   = Gasxsecs_fine(w,n,j,g) * Gas_Fine(n,g,j)
               hdiff = hght1 - hght2
               tau   = tau + hdiff * (od1 + od2)
               hght1 = hght2 ; od1 = od2
            enddo

!  Lower boundary values

            hght2 = level_heights(n)
            od2   = gas_xsecs(w,n,g) * LevelGas(n,g)
            hdiff = hght1 - hght2
            tau = tau + hdiff * ( od1 + od2 )      

!  Set output

            odabs(n,g) = tau

!  End gas and layer loops

         endif
      enddo
   enddo

!  total molecular absorption

   do n = 1, nlayers
      odabstot(n) = SUM(odabs(n,1:ngases))
   enddo

!  Rayleigh scattering absorption (commented out)
!  ==============================

!   do n = 1, nlayers
!      odair(n) = aircolumns(n)  * Rayleigh_xsecs(w)
!   enddo

!  End subroutine

   return
end subroutine GEMSTOOL_NSW_ODSETTER


subroutine GEMSTOOL_NSW_ODSETTER_PLUS_1 &
     ( maxlayers, maxfinelayers, maxwav, maxgases,           & ! Input dimensions
       nlayers, nfinediv, ngases, do_gases, w,               & ! Input control
       level_heights, levelgas, dLevelGas_dV, gas_xsecs,     & ! Input (LEVELS)
       heights_fine,  Gas_fine, dGasFine_dV,  Gasxsecs_fine, & ! Fine-grid input
       odabs, dodabs_dV, odabstot )                             ! Output

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input arguments
!  ---------------

!  Dimensioning

   integer  , intent(in) :: maxlayers, maxfinelayers, maxwav, maxgases

!  Layer control

   INTEGER  , intent(in)  :: NLAYERS
   INTEGER  , intent(in)  :: NFINEDIV ( MAXLAYERS)


!  Gas control

   integer  , intent(in)  :: NGASES
   logical  , intent(in)  :: do_GASES( maxgases )

!  Spectral index

   integer  , intent(in)  :: w

!  Level Heights and Gas values + cross-sections

   real(fpk), intent(in)  :: LEVEL_HEIGHTS ( 0:MAXLAYERS) 
   real(fpk), intent(in)  :: LEVELGAS      ( 0:MAXLAYERS, maxgases )
   real(fpk), intent(in)  :: dLEVELGAS_dV  ( 0:MAXLAYERS, maxgases, 2 )
   REAL(fpk), intent(in)  :: GAS_XSECS     ( MAXWAV, 0:maxlayers, maxgases ) 

!  Heights and Gas values, cross-sections: FINE-LEVEL values

   real(fpk), intent(in)  :: HEIGHTS_FINE  ( MAXLAYERS, maxfinelayers ) 
   real(fpk), intent(in)  :: GAS_FINE      ( MAXLAYERS, maxfinelayers, maxgases )
   real(fpk), intent(in)  :: dGASFINE_dV   ( MAXLAYERS, maxfinelayers, maxgases, 2 )
   REAL(fpk), intent(in)  :: GASXSECS_FINE ( MAXWAV, maxlayers, maxfinelayers, maxgases ) 

!  Output
!  ------

!  Gas absoprtion optical depths

   real(fpk), intent(out) :: odabs     ( MAXLAYERS, maxgases )
   real(fpk), intent(out) :: dodabs_dV ( MAXLAYERS, maxgases, 2 )
   real(fpk), intent(out) :: odabstot  ( MAXLAYERS )

!  Local variables
!  ---------------

   integer   :: n, n1, j, g
   real(fpk) :: hght1, hght2, od1, od2, dod1_dV(2), dod2_dV(2), hdiff, tau, dtau_dV(2)

!  Initialize

   odabs     = 0.0_fpk
   dodabs_dV = 0.0_fpk
   odabstot  = 0.0_fpk

!  Start layer and gas loops

   do n = 1, nlayers
      n1 = n - 1
      do g = 1, ngases
         if ( do_gases(g) ) then

!  Initialize

            tau = 0.0_fpk ; dtau_dV = 0.0_fpk

!  Upper boundary values

            hght1      = level_heights(n1)
            od1        = gas_xsecs(w,n1,g) * LevelGas(n1,g)
            dod1_dV(1) = gas_xsecs(w,n1,g) * dLevelGas_dV(n1,g,1) ; dod1_dV(2)   = 0.0_fpk

!  Loop over fine layer subdivisions

            do j = 1, nfinediv(n)

               hght2        = heights_fine(n,j)
               od2          = Gasxsecs_fine(w,n,j,g) * Gas_Fine(n,g,j)
               dod2_dV(1:2) = Gasxsecs_fine(w,n,j,g) * dGasFine_dV(n,g,j,1:2)

               hdiff        = hght1 - hght2
               tau          = tau + hdiff * (od1 + od2)
               dtau_dV(1:2) = dtau_dV(1:2)  + hdiff * (dod1_dV(1:2) + dod2_dV(1:2))

               hght1 = hght2 ; od1 = od2 ; dod1_dV = dod2_dV

            enddo

!  Lower boundary values

            hght2        = level_heights(n)
            od2          = gas_xsecs(w,n,g) * LevelGas(n,g)
            dod2_dV(1)   = 0.0_fpk ; dod2_dV(2)  = gas_xsecs(w,n,g) * dLevelGas_dV(n,g,2)

            hdiff        = hght1 - hght2       
            tau          = Tau + hdiff * ( od1 + od2 )      
            dtau_dV(1:2) = dtau_dV(1:2)  + hdiff * (dod1_dV(1:2) + dod2_dV(1:2))

!  Set output

            odabs(n,g)         = tau
            dodabs_dV(n,g,1:2) = dtau_dV(1:2) 

!  End gas and layer loops

         endif
      enddo
   enddo

!  total molecular absorption

   do n = 1, nlayers
      odabstot(n) = SUM(odabs(n,1:ngases))
   enddo

!  Rayleigh scattering absorption (commented out)
!  ==============================

!   do n = 1, nlayers
!      odair(n) = aircolumns(n)  * Rayleigh_xsecs(w)
!   enddo

!  End subroutine

   return
end subroutine GEMSTOOL_NSW_ODSETTER_PLUS_1

!subroutine GEMSTOOL_NSW_ODSETTER_PLUS_2 &
!     ( maxlayers, maxfinelayers, maxwav, maxgases, & ! Input dimensions
!       nlayers, nfinediv, ngases, which_gases, w,  &  ! Input control
!       level_heights, level_temps, level_press, level_logpress,  & ! Input (LEVELS)
!       levelgas, dLevelGas_dV, gas_xsecs,                 & ! Input (LEVELS)
!       heights_fine, temps_fine, press_fine, Gas_fine,     & ! Fine-grid input
!       dGasfine_dV,  dGasfine_dS, dGasfine_dP, dTfine_dS,  & ! Fine-grid input
!       Gasxsecs_fine, dGasxsecsfine_dP, dGasxsecsfine_dT,  & ! Fine-grid input
!       odabs, dodabs_dV, dodabs_dS, dodabs_dP, odabstot )     ! Output
!  End subroutine
!   return
!end subroutine GEMSTOOL_NSW_ODSETTER_PLUS_2

!  End module

end module GEMSTOOL_NSW_ODSETTER_m

