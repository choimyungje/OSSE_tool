
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

!  Developed for LRRS Version 2.5, 9/8/15. R. Spurr, RT SOLUTIONS Inc.
!   Closely follows the LIDORT module with the same name, but
!     (1) No solar angle dimensioning, No observational geometry.
!     (2) Additional wavelength dimensioning (MAX_SLEAVE_POINTS)
!     (3) Additional control for using 1 or all wavelengths

!  Note (9/8/15). The number of wavelengths for SLEAVE has deliberately
!  been left flexible - typically the SLEAVE properties will change very
!  little over a Raman-scattering window (+/- 2 nm in the UV), so to
!  a very good approximation, it is often sufficient to use one point for
!  the calculations of SLEAVE - in which case, MAX_SLEAVE_POINTS will be 1
!  This applies equally to Waterleaving and Fluorescence, where wavelength
!  variation can be important in the Visible and NIR.

!  Now, whenever the SLEAVE supplement is used with LRRS, the choice of 
!  SLEAVE wavelengths is linked to the Raman wavelengths (LAMDAS_RANKED).
!  These wavelengths are now input to the SLEAVE supplement MASTER, and
!  they are not set by hand or by configuration-file read. 

!  However, there is a hard-wired choice (DO_WAV1) which allows you
!  to choose a single wavelength point for calculation. When this flag is
!  set, the number of SLEAVE wavelengths = 1 (Waterleaving or Fluorescence)
!  and the single wavelength is set to the AVERAGE VALUE of LAMBDAS_RANKED.

!  If you want all wavelengths, then control flag DO_Wav1 is False,
!  and number of SLEAVE wavelengths  N_Lambdas_Ranked, and the
!  wavelengths themselves are just copied from LAMBDAS_RANKED
!  A SLEAVE calculation will be done for all Raman-scattered wavelengths!!!

! ###############################################################
! #                                                             #
! # Water-Leaving Subroutines in this Module                    #
! #                                                             #
! #         WaterLeaving (Top-level)                            #
! #           - Water_RefracIndex (formerly INDWAT)             #
! #           - Ocean_Reflectance (formerly MORCASIWAT)         #
! #           - WhiteCap_Reflectance                            #
! #           - Water_Transmittance_Quads                       #
! #           - Water_Transmittance                             #
! #               * GENERAL_SUNGLINT                            #
! #                  --> Fresnel_Complex                        #
! #                                                             #
! # Fluorescence Subroutines in this Module                     #
! #                                                             #
! #              get_fluorescence_755                           #
! #              average_solar_cosine                           #
! #              solar_spec_irradiance                          #
! #                                                             #
! ###############################################################

      MODULE sleave_sup_routines_m

      use sleave_sup_aux_m

      PRIVATE
      PUBLIC :: WaterLeaving,              &
                Water_RefracIndex,         &
                Water_Transmittance_Quads, &
                Fresnel_Complex,           &
                get_fluorescence_755,      &
                solar_spec_irradiance

      CONTAINS

subroutine WaterLeaving &
     ( Maxvzas, MAX_STREAMS,                                         &
       do_Isotropy, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy, &
       Wavelength, Salinity, PigmentConc, Windspeed, WindSunAngle,   &
       nvzas, nstreams, sza, vzas, streams,                          &
       WLeaving_ISO, WLeaving_SD, WLeaving_SV )

! mick fix 12/28/2014 - Using normalization correction by A Sayer, 04 Nov 2014
! Apply a normalisation factor of 1/(pi/mu0) to output water-leaving reflectance, to
! bring things in line with results from e.g. 6S simulations and expected behaviour.
! Think this is a subtlety related to reflectance vs. normalised radiance treatment,
! although it is very obvious if you don't do it. Correction applied at end of the
! subroutine.

   IMPLICIT NONE
   INTEGER  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  This is a Stand-alone subroutine.

!  inputs
!  ======

!  Dimensioning

   integer  , intent(in)   :: Maxvzas, MAX_STREAMS

!  Logical flags
!  -------------

!  Isotropic (Fast Calculation) option, assumes all transmittances = 1

   Logical  , intent(in)   :: do_Isotropy

!  Optional inclusion of Foam term

   Logical  , intent(in)   :: Do_FoamOption

!  Optional inclusion of Shadow term for Glitter (Air-water only?)

   Logical  , intent(in)   :: Do_GlintShadow

!  Flag for using Isotropic Facet distribution

   LOGICAL  , intent(in)   :: Do_FacetIsotropy

!  Physical
!  --------

!  Wavelength in Micrometers

   real(fpk), intent(in)   :: Wavelength

!  Salinity

   real(fpk), intent(in)   :: Salinity

!  Pigment concentration

   real(fpk), intent(in)   :: PigmentConc

!  Windspeed m/s

   REAL(fpk), intent(in)    :: WINDSPEED

!  Azimuth angle in Radians

   REAL(fpk), intent(in)    :: WindSunAngle

!  Sun, viewing and stream angles
!  ------------------------------

   integer  , intent(in) :: nvzas, nstreams
   real(fpk), intent(in) :: sza
   real(fpk), intent(in) :: vzas   (Maxvzas)
   real(fpk), intent(in) :: streams(MAX_STREAMS)

!  OUTPUT
!  ======

!  Isotropic value. Fast calculation

   REAL(fpk), intent(out)    :: WLeaving_ISO

!  Input solar, output stream angles

   REAL(fpk), intent(out)    :: WLeaving_SD ( MAX_STREAMS )

!  input solar, output view angles

   REAL(fpk), intent(out)    :: WLeaving_SV ( Maxvzas )

!  remark: Still no Azimuth dependence here.....

!  LOCAL
!  =====

!  Transmittance Quadratures
!  -------------------------

   integer, parameter   :: Max_PolarQuads=24, Max_AzimQuads=48
   real(fpk) :: PolarQuads    (Max_PolarQuads)   ! Radians
   real(fpk) :: CosPolarQuads (Max_PolarQuads)
   real(fpk) :: SinPolarQuads (Max_PolarQuads)
   real(fpk) :: PolarWeights  (Max_PolarQuads)

   real(fpk) :: AzimQuads    (Max_AzimQuads)     ! Radians
   real(fpk) :: CosAzimQuads (Max_AzimQuads)
   real(fpk) :: SinAzimQuads (Max_AzimQuads)
   real(fpk) :: AzimWeights  (Max_AzimQuads)

!  Glitter control
!  ---------------

   logical   :: do_coeffs, local_do_shadow
   REAL(fpk) :: SUNGLINT_COEFFS(7), Refrac_R, Refrac_I, Refrac_sq
   real(fpk) :: phi_w, cphi_w, sphi_w, local_refrac_R, local_refrac_I

!  other variables
!  ---------------

!  Help

   logical   :: noWL
   integer   :: I
   real(fpk) :: dtr, pi, Albedo, Const, Rwprime, incident_angle, Local_Sine
   real(fpk) :: Foam_correction, WC_Reflectance, WC_Lambertian
   real(fpk) :: Ocean_Reflec, SZA_cosine

!  transmittance help variables

   logical   :: do_transmittances
   real(fpk) :: Trans_Norm, Trans_solar
   real(fpk) :: Trans_stream(MAX_STREAMS), Trans_Viewing(Maxvzas)

!  Parameters

   real(fpk), parameter :: zero = 0.0_fpk
   real(fpk), parameter :: one  = 1.0_fpk

!  Initial Setup
!  -------------

!  conversions

   pi = acos(-1.0_fpk)
   dtr = pi / 180.0_fpk
   SZA_cosine = cos(sza*dtr)

!  Zero the output

   WLeaving_ISO = zero
   WLeaving_SD  = zero
   WLeaving_SV  = zero

!  Refractive index. Formerly INDWAT

   Call  Water_RefracIndex  ( Wavelength, Salinity, Refrac_R, Refrac_I )
   Refrac_sq = Refrac_R * Refrac_R + Refrac_I * Refrac_I

!  Ocean Leaving ; Exit if no contribution. Formerly MORCASIWAT

   call  Ocean_Reflectance &
     ( SZA_cosine, Wavelength, PigmentConc, noWL, Ocean_Reflec )
   if ( noWL ) return

!  Foam-reflectance correction.

   Foam_correction = one
   if ( Do_FoamOption ) then
      call WhiteCap_Reflectance &
         ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian )
      Foam_correction = one - WC_Reflectance
   endif

!  Initial debug...........
!   write(*,*)Ocean_Reflec
!   write(*,*)WC_Reflectance,WC_Lambertian, Refrac_sq
!   pause'after first routines'

!  set Coeffs flag, initialize local shadow flag

   do_coeffs       = .true.
   local_do_shadow = do_GlintShadow

!  Initialize Transmittances to 1.0

   trans_norm    = one
   Trans_solar   = one 
   Trans_viewing = one
   Trans_stream  = one

!  Set Transmittances flag. Not required fort he Fast calculation option

   Trans_solar = one ; trans_norm = one
   do_transmittances = .false.
   if ( .not. Do_Isotropy ) then
      if ( Ocean_Reflec .gt.0.0001 ) do_transmittances =.true.
   endif

!  Get quadratures

   if ( do_transmittances ) then
      call Water_Transmittance_Quads &
       ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
         PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Output
         AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Output
         TRANS_NORM )
   endif

!   write(*,*) trans_norm ; pause'after quads'

!  Skip if no transmittance calculation (values are initialized to 1.0)

   if ( .not. do_transmittances ) go to 67

!  Downward Solar transmittances
!     Set shadow flag for this

   if ( Ocean_Reflec.gt.0.0001 ) then
      phi_w = WindSunAngle
      cphi_w = cos(phi_w*dtr)
      sphi_w = sin(phi_w*dtr)
      local_do_shadow = do_GlintShadow
      call Water_Transmittance &
         ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
           PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
           AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
           do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
           sza, REFRAC_R, REFRAC_I,                                & ! Input
           WINDSPEED, SUNGLINT_COEFFS, PHI_W, CPHI_W, SPHI_W,      & ! Input
           TRANS_NORM, TRANS_SOLAR )
   endif

!  Upward transmittances into View directions
!     Passing from water to air, use Snell's Law.  no absorption
!     Local shadow flag turned off here.

   local_do_shadow  = .false.
   local_refrac_R   = one / refrac_R
   local_refrac_I   = zero
   if ( Ocean_Reflec.gt.0.0001 ) then
      do i = 1, nvzas
         incident_angle = asin(sin(vzas(i) * dtr)/refrac_R)/dtr
         call Water_Transmittance &
            ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
              PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
              AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
              do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
              incident_angle, local_refrac_R, local_refrac_I,         & ! Input
              WINDSPEED, SUNGLINT_COEFFS, PHI_W, CPHI_W, SPHI_W,      & ! Input
              TRANS_NORM, TRANS_VIEWING(I) )
      enddo
   endif

!  Upward transmittances into stream directions
!     Passing from water to air, use Snell's Law.  no absorption

   local_do_shadow  = .false.
   local_refrac_R   = one / refrac_R
   local_refrac_I   = zero
   if ( Ocean_Reflec.gt.0.0001 ) then
      do i = 1, nstreams
         local_sine = sqrt ( one - streams(i) * streams(i) )
         incident_angle = asin(local_sine/refrac_R)/dtr
         call Water_Transmittance &
            ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
              PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
              AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
              do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
              incident_angle, local_refrac_R, local_refrac_I,         & ! Input
              WINDSPEED, SUNGLINT_COEFFS, PHI_W, CPHI_W, SPHI_W,      & ! Input
              TRANS_NORM, TRANS_STREAM(I) )
      enddo
   endif

!  continuation point for avoiding transmittance calculations

67 continue

!  Final computation

   Albedo  = 0.485_fpk
   Rwprime = Ocean_Reflec / ( one - Albedo * Ocean_Reflec  )

   Const = ( one / refrac_sq ) * trans_solar * Rwprime
   Const = Const * Foam_correction
   WLeaving_Iso = Const

   do i = 1, nstreams
      WLeaving_SD(I) = Const * Trans_Stream(I)
   enddo
   do i = 1, nvzas
      WLeaving_SV(I) = Const * Trans_Viewing(I)
   enddo

   ! Correction of normalisation factor: divide by (pi/mu0). A Sayer 04 Nov 2014.
   ! Have kept this outside the above loops so it is more obvious.
   WLeaving_Iso = WLeaving_Iso/(pi/cos(sza*dtr))
   do i = 1, nstreams
      WLeaving_SD(I) = WLeaving_SD(I)/(pi/cos(sza*dtr))
   enddo
   do i = 1, nvzas
      WLeaving_SV(I) = WLeaving_SV(I)/(pi/cos(sza*dtr))
   enddo

! Finish

   return
end subroutine WaterLeaving

!
subroutine WhiteCap_Reflectance &
    ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian )

!  Stand-alone routine for computing the WhiteCap Reflectance
!   Based on 6S code, as updated by A. Sayer (2011)

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!    (Wind speed in [m/s], Wavelength in Microns)

   real(fpk), intent(in)  :: WindSpeed
   real(fpk), intent(in)  :: Wavelength

!  output

   real(fpk), intent(out) :: WC_Reflectance
   real(fpk), intent(out) :: WC_Lambertian

!  Data
!  ----

!  Single precision

   real :: Effective_WCRef(39)

! effective reflectance of the whitecaps (Koepke, 1984)
! These are the original values - superseded, A Sayer 05 Jul 2011.
!      data Effective_WCRef/ &
!     0.220,0.220,0.220,0.220,0.220,0.220,0.215,0.210,0.200,0.190,&
!     0.175,0.155,0.130,0.080,0.100,0.105,0.100,0.080,0.045,0.055,&
!     0.065,0.060,0.055,0.040,0.000,0.000,0.000,0.000,0.000,0.000,&
!     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

! effective reflectance of the whitecaps (Frouin et al, 1996)
! Assume linear trends between the node points they give
! This is the spectral shape

      data Effective_WCRef/ &
     1.000,1.000,1.000,1.000,0.950,0.900,0.700,0.550,0.500,0.450,&
     0.400,0.350,0.300,0.250,0.200,0.150,0.100,0.050,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

!  Local variables
!  ---------------

!  Single precision in the original code

   integer :: iwl, iref
   real    :: Wlb, WLP, Ref(39), wspd, wl, Ref_i, Rwc

!  Initialize

   WC_Reflectance = 0.0_fpk
   WC_Lambertian  = 0.0_fpk

!  Single precision inputs in the original

   wspd = real(WindSpeed)
   wl   = real(Wavelength)

!  Scale data for value of 0.22 in the midvisible.

   DO iref = 1,39
      Ref(iref) = 0.22 * Effective_WCRef(iref)
   ENDDO

!  COMPUTE WHITECAPS REFLECTANCE (LAMBERTIAN)

   Wlb    = 0.0
   IF (wspd .le. 9.25) THEN
       Wlb = 0.01*((3.18e-03)*((wspd-3.7)**3.0))
   ELSE IF (wspd .gt. 9.25) THEN
       Wlb = 0.01*((4.82e-04)*((wspd+1.8)**3.0))
   END IF

! Original whitecap calculation - superseded, A. Sayer 05 Jul 2011.
!      W=2.95e-06*(wspd**3.52)

!  Find data point, Linearly interpolate

   iwl   = 1+int((wl-0.2)/0.1)
   wlp   = 0.5+(iwl-1)*0.1
   Ref_i = Ref(iwl+1) + ( wl-wlp)/0.1*(Ref(iwl)-Ref(iwl+1))
   Rwc   = Wlb*Ref_i

!  Final values

   WC_Lambertian  = real(Wlb,fpk)
   WC_Reflectance = real(Rwc,fpk)

!  Finish

   return
end subroutine WhiteCap_Reflectance

Subroutine Fresnel_Complex ( MR, MI, COSCHI, FP )

  implicit none
  integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Arguments

  real(fpk), intent(in)  :: MR, MI, COSCHI
  real(fpk), intent(out) :: FP

!  Local

  real(fpk) :: MRSQ, MISQ, MSQ, MRMI2, SINCHI_SQ, AA, A1, A2, B1, B2
  real(fpk) :: U, V, VSQ, CMU, CPU, RR2
  real(fpk) :: B1MU, B1PU, B2MV, B2PV, RL2

!  Calculation of FP, Complex RI

   IF ( MI.eq.0.0_fpk) goto 67

   MRSQ = MR * MR ; MISQ = MI * MI
   MSQ   = MRSQ - MISQ
   MRMI2 = 2.0_fpk * MR * MI

   SINCHI_SQ = 1.0_fpk - COSCHI * COSCHI 
   AA = MSQ - SINCHI_SQ
   A1 = abs(AA)
   A2 = SQRT ( AA*AA + MRMI2 * MRMI2 )

   U = sqrt(0.5_fpk*abs(A1+A2))
   V = sqrt(0.5_fpk*abs(-A1+A2))
   VSQ = V * V
   CMU = ( COSCHI - U ) ; CPU = ( COSCHI + U )
   RR2 = ( CMU*CMU + VSQ ) / ( CPU*CPU + VSQ )

   B1 = MSQ * COSCHI
   B2 = MRMI2 * COSCHI
   B1MU = B1 - U ; B1PU = B1 + U 
   B2PV = B2 + V ; B2MV = B2 - V 

   RL2 = ( B1MU*B1MU + B2PV*B2PV ) / ( B1PU*B1PU + B2MV*B2MV )
   FP = 0.5_fpk * ( RR2 + RL2 )
   return

!  Calculation of FP. Real RI

67 continue
   MSQ = MR * MR
   SINCHI_SQ = 1.0_fpk - COSCHI * COSCHI 
   U = sqrt(abs(MSQ - SINCHI_SQ))
   CMU = ( COSCHI - U ) ; CPU = ( COSCHI + U )
   RR2 = CMU*CMU / ( CPU*CPU )
   B1 = MSQ * COSCHI
   B1MU = B1 - U ; B1PU = B1 + U 
   RL2 = B1MU*B1MU / ( B1PU*B1PU )
   FP = 0.5_fpk * ( RR2 + RL2 )

!  Finish

   return
end subroutine Fresnel_Complex

subroutine Ocean_Reflectance &
       ( SZA_cosine, Wavelength, PigmentConc, noWL, Ocean_Reflec )

!  THIS IS FORMERLY CALLED "MORCASIWAT", as modified by A. Sayer for 6S

! mick fix 12/28/2014 - Using updates by A Sayer November 03 2014:
! Extended functionality down to 200 nm. Achieved by:
! - Extended data arrays down to 200 nm (another 40 elements).
! - Changed logic check for contribution to 0.2-0.9 microns from 0.4-0.9 microns, and started table
!   lookup calculation from 0.2 microns instead of 0.4 microns.
! Note, this is based on a simple extension of the published optical model for vis wavelengths.
!   Possible that other scatterers/absorbers
! which are neglected in this model may be important at UV wavelengths.
! Do linear interpolation of optical property LUTs, rather than nearest neighbour, to remove discontinuities. Achieved by:
! - Replicated final element of LUTs to avoid potential for extrapolation errors.
! - Replace nint() call with floor() call to correctly get lower bound
! - Define variable dwl, fractional distance along the 5 nm LUT grid
! - Implement the interpolation using dwl rather than direct lookup of nearest value.
! Also:
! - Corrected Prieur and Sathyendranath, Limnol. Oceanogr. reference year to 1981 instead of 1983.
! - Corrected typo in water scattering coefficient at 410 nm: 0.0068 was written instead of 0.0061.
!   Removes artificial spike at 410 nm in calculated reflectance.

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input/Output
!  ------------


   real(fpk), intent(in)   :: SZA_cosine

   real(fpk), intent(in)   :: Wavelength
   real(fpk), intent(in)   :: PigmentConc

   logical  , intent(out)  :: noWL
   real(fpk), intent(out)  :: Ocean_Reflec

!  Local
!  -----

!      subroutine morcasiwat(wl,C,R2,mu_sol)
! Rewritten, beginning 07 July 2011, Andrew Sayer
! Now extends underlight calculations out to 400-900 nm.
! and more modern formulations,
! but still based on Case 1 principle.

! Spectral diffuse attenuation coefficient of Case I Waters as Predicted
! by MOREL within the spectral range 400-700nm (1988, Journal of Geophysical
! Research, Vol.93, No C9, pp 10749-10768)
!
! input parameters:	wl wavelength (IN MICROMETERS)
!			C  pigment concentration
!                       mu_sol : cosine of solar zenith angle
! output parameter:	R2  reflectance of water below the surface

! Tabulated absorption coefficient for water, scattering coefficient for water,
! and absorption coefficient for chlorophyll-a, tabulated from 200 nm to 900 nm
! in 5 nm increments
      real(fpk) water_abs_coeff(142),water_scat_coeff(142),abs_coeff_chl(142)
! Input/output parameters
      real(fpk) mu_sol,r2,C,wl
! Absorption/scattering terms, and parameters for calculation of f
      real(fpk) a_wat,b_wat,b_wat_all,a_chl,f,a_tot,b_tot,a_ph,a_cdom,v,bp,bbp,eta,dwl
! Wavelength index for tables
      integer iwl, J

! Smith and Baker, AO(20) 1981, table 1, out to 800 nm. This has 10 nm increments so
! linearly interpolate between these.
! Hale & Qurry, AO(12) 1973, table 1, for 805-900 nm. This has 25 nm increments so
! linearly interpolate between these. Provided as extinction coefficient
! so convert using a=4*pi*k/lambda (note lambda in m for units m^{-1})
      data water_abs_coeff/ &
       3.0700_fpk,2.5300_fpk,1.9900_fpk,1.6500_fpk,1.3100_fpk,&
       1.1185_fpk,0.9270_fpk,0.8235_fpk,0.7200_fpk,0.6395_fpk,&
       0.5590_fpk,0.5080_fpk,0.4570_fpk,0.4150_fpk,0.3730_fpk,&
       0.3305_fpk,0.2880_fpk,0.2515_fpk,0.2150_fpk,0.1780_fpk,&
       0.1410_fpk,0.1230_fpk,0.1050_fpk,0.0907_fpk,0.0844_fpk,&
       0.0761_fpk,0.0678_fpk,0.0620_fpk,0.0561_fpk,0.0512_fpk,&
       0.0463_fpk,0.0421_fpk,0.0379_fpk,0.0340_fpk,0.0300_fpk,&
       0.0260_fpk,0.0220_fpk,0.0206_fpk,0.0191_fpk,0.0181_fpk,&
       0.0171_fpk,0.0166_fpk,0.0162_fpk,0.0158_fpk,0.0153_fpk,&
       0.0149_fpk,0.0144_fpk,0.0144_fpk,0.0145_fpk,0.0145_fpk,&
       0.0145_fpk,0.0150_fpk,0.0156_fpk,0.0156_fpk,0.0156_fpk,&
       0.0166_fpk,0.0176_fpk,0.0186_fpk,0.0196_fpk,0.0227_fpk,&
       0.0257_fpk,0.0307_fpk,0.0357_fpk,0.0417_fpk,0.0477_fpk,&
       0.0492_fpk,0.0507_fpk,0.0532_fpk,0.0558_fpk,0.0598_fpk,&
       0.0638_fpk,0.0673_fpk,0.0708_fpk,0.0753_fpk,0.0799_fpk,&
       0.0940_fpk,0.1080_fpk,0.1330_fpk,0.1570_fpk,0.2005_fpk,&
       0.2440_fpk,0.2660_fpk,0.2890_fpk,0.2990_fpk,0.3090_fpk,&
       0.3145_fpk,0.3190_fpk,0.3245_fpk,0.3290_fpk,0.3390_fpk,&
       0.3490_fpk,0.3740_fpk,0.4000_fpk,0.4150_fpk,0.4300_fpk,&
       0.4400_fpk,0.4500_fpk,0.4750_fpk,0.5000_fpk,0.5750_fpk,&
       0.6500_fpk,0.7445_fpk,0.8390_fpk,1.0040_fpk,1.1690_fpk,&
       1.4840_fpk,1.7990_fpk,2.0895_fpk,2.3800_fpk,2.4250_fpk,&
       2.4700_fpk,2.5100_fpk,2.5500_fpk,2.5300_fpk,2.5100_fpk,&
       2.4350_fpk,2.3600_fpk,2.2600_fpk,2.1600_fpk,2.1150_fpk,&
       2.0700_fpk,2.2104_fpk,2.3509_fpk,2.4913_fpk,2.6318_fpk,&
       2.7722_fpk,3.0841_fpk,3.3960_fpk,3.7079_fpk,4.0198_fpk,&
       4.3317_fpk,4.5884_fpk,4.8451_fpk,5.1019_fpk,5.3586_fpk,&
       5.6153_fpk,5.8495_fpk,6.0836_fpk,6.3177_fpk,6.5518_fpk,&
       6.7858_fpk,6.7858_fpk/

! Smith and Baker, 1981 out to 800 nm.
! This is again at 10 nm increments, interpolated to 5 nm.
! Set to 0.0003 out to 870 nm and 0.0002 after this
! based on Morel's power law (lambda^-4.32)
! and the fact that results not very sensitive to
! the exact value here.
      data water_scat_coeff/&
       0.1510_fpk,0.1350_fpk,0.1190_fpk,0.1093_fpk,0.0995_fpk,&
       0.0908_fpk,0.0820_fpk,0.0753_fpk,0.0685_fpk,0.0630_fpk,&
       0.0575_fpk,0.0530_fpk,0.0485_fpk,0.0450_fpk,0.0415_fpk,&
       0.0384_fpk,0.0353_fpk,0.0329_fpk,0.0305_fpk,0.0284_fpk,&
       0.0262_fpk,0.0246_fpk,0.0229_fpk,0.0215_fpk,0.0200_fpk,&
       0.0188_fpk,0.0175_fpk,0.0164_fpk,0.0153_fpk,0.0144_fpk,&
       0.0134_fpk,0.0127_fpk,0.0120_fpk,0.0113_fpk,0.0106_fpk,&
       0.0100_fpk,0.0094_fpk,0.0089_fpk,0.0084_fpk,0.0080_fpk,&
       0.0076_fpk,0.0072_fpk,0.0068_fpk,0.0064_fpk,0.0061_fpk,&
       0.0058_fpk,0.0055_fpk,0.0052_fpk,0.0049_fpk,0.0047_fpk,&
       0.0045_fpk,0.0043_fpk,0.0041_fpk,0.0039_fpk,0.0037_fpk,&
       0.0036_fpk,0.0034_fpk,0.0033_fpk,0.0031_fpk,0.0030_fpk,&
       0.0029_fpk,0.0028_fpk,0.0026_fpk,0.0025_fpk,0.0024_fpk,&
       0.0023_fpk,0.0022_fpk,0.0022_fpk,0.0021_fpk,0.0020_fpk,&
       0.0019_fpk,0.0019_fpk,0.0018_fpk,0.0018_fpk,0.0017_fpk,&
       0.0017_fpk,0.0016_fpk,0.0016_fpk,0.0015_fpk,0.0015_fpk,&
       0.0014_fpk,0.0014_fpk,0.0013_fpk,0.0013_fpk,0.0012_fpk,&
       0.0012_fpk,0.0011_fpk,0.0011_fpk,0.0010_fpk,0.0010_fpk,&
       0.0010_fpk,0.0009_fpk,0.0008_fpk,0.0008_fpk,0.0008_fpk,&
       0.0008_fpk,0.0007_fpk,0.0007_fpk,0.0007_fpk,0.0007_fpk,&
       0.0007_fpk,0.0007_fpk,0.0007_fpk,0.0007_fpk,0.0006_fpk,&
       0.0006_fpk,0.0006_fpk,0.0006_fpk,0.0006_fpk,0.0006_fpk,&
       0.0005_fpk,0.0005_fpk,0.0005_fpk,0.0005_fpk,0.0005_fpk,&
       0.0005_fpk,0.0004_fpk,0.0004_fpk,0.0004_fpk,0.0004_fpk,&
       0.0004_fpk,0.0003_fpk,0.0003_fpk,0.0003_fpk,0.0003_fpk,&
       0.0003_fpk,0.0003_fpk,0.0003_fpk,0.0003_fpk,0.0003_fpk,&
       0.0003_fpk,0.0003_fpk,0.0003_fpk,0.0003_fpk,0.0003_fpk,&
       0.0002_fpk,0.0002_fpk,0.0002_fpk,0.0002_fpk,0.0002_fpk,&
       0.0002_fpk,0.0002_fpk/

! Prieur and Sathyendranath, Limnol. Oceanogr. 26, 1983 table 2
! assumed zero after 700 nm: no data given and water absorption becomes very dominating here
      data abs_coeff_chl/&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.053_fpk,0.123_fpk,0.195_fpk,0.264_fpk,&
       0.335_fpk,0.405_fpk,0.476_fpk,0.546_fpk,0.617_fpk,&
       0.687_fpk,0.781_fpk,0.828_fpk,0.883_fpk,0.913_fpk,&
       0.939_fpk,0.973_fpk,1.001_fpk,1.000_fpk,0.971_fpk,&
       0.944_fpk,0.928_fpk,0.917_fpk,0.902_fpk,0.870_fpk,&
       0.839_fpk,0.798_fpk,0.773_fpk,0.750_fpk,0.717_fpk,&
       0.668_fpk,0.645_fpk,0.618_fpk,0.582_fpk,0.528_fpk,&
       0.504_fpk,0.474_fpk,0.444_fpk,0.416_fpk,0.384_fpk,&
       0.357_fpk,0.321_fpk,0.294_fpk,0.273_fpk,0.276_fpk,&
       0.268_fpk,0.291_fpk,0.274_fpk,0.282_fpk,0.249_fpk,&
       0.236_fpk,0.279_fpk,0.252_fpk,0.268_fpk,0.276_fpk,&
       0.299_fpk,0.317_fpk,0.333_fpk,0.334_fpk,0.326_fpk,&
       0.356_fpk,0.389_fpk,0.441_fpk,0.534_fpk,0.595_fpk,&
       0.544_fpk,0.502_fpk,0.420_fpk,0.329_fpk,0.262_fpk,&
       0.215_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,0.000_fpk,&
       0.000_fpk,0.000_fpk/

!  Initialize

      wl     = Wavelength
      C      = PigmentConc
      noWL = .false.  ; Ocean_Reflec = 0.0_fpk

! If wavelength out of range, no need to calculate underlight

      if (wl.lt.0.200_fpk.or.wl.gt.0.900_fpk)then
        noWL = .true. ; goto 60
      endif

! Extract tabulated values for this wavelength
      iwl=1+floor((wl-0.200_fpk)/0.005_fpk)
! 03 Nov 2014 A Sayer now linear interpolation rather than nearest neighbour
!      a_wat=water_abs_coeff(iwl)
!      b_wat_all=water_scat_coeff(iwl)
!      a_chl=abs_coeff_chl(iwl)
      dwl=(wl-0.200_fpk)/0.005_fpk-floor((wl-0.200_fpk)/0.005_fpk)
      a_wat=water_abs_coeff(iwl)+dwl*(water_abs_coeff(iwl+1)-water_abs_coeff(iwl))
      b_wat_all=water_scat_coeff(iwl)+dwl*(water_scat_coeff(iwl+1)-water_scat_coeff(iwl))
      a_chl=abs_coeff_chl(iwl)+dwl*(abs_coeff_chl(iwl+1)-abs_coeff_chl(iwl))

! Morel and Maritorena, 2001
      a_ph=0.06_fpk*a_chl*(C**0.65_fpk)

! Equations 2 and 4a from Morel and Gentili, RSE 113, 2009
      a_cdom=0.0524_fpk*(C**0.63_fpk)*exp(-0.018_fpk*(wl*1000._fpk-412._fpk))

      a_tot=a_wat + a_ph + a_cdom

      b_wat=0.5_fpk*b_wat_all

! Morel and Maritorena, 2001 (also earlier work by Loisel and Morel)
! exponent for spectral dependence of scattering
      if (C .le. 2._fpk) then
        v=0.5_fpk*(log10(C)-0.3_fpk)
      endif
      if (C .gt. 2._fpk) then
        v=0._fpk
      endif

      bp=0.416_fpk*(C**0.766_fpk)

      bbp=0.002_fpk+0.01_fpk*(0.5_fpk-0.25_fpk*log10(C))*((wl/0.55_fpk)**v)

      b_tot=b_wat + bbp*bp
      eta=b_wat/b_tot

      !write(*,*)
      !write(*,*)bbp
      !write(*,*)b_tot
      !write(*,*)eta

! Morel and Gentili, 1991.
      mu_sol = SZA_cosine
      f=0.6279_fpk - (0.2227_fpk*eta) - (0.0513_fpk*eta*eta) + (0.2465_fpk*eta -0.3119_fpk)*mu_sol

! R=f*b(a+b) : constant multiplied by ratio of backscatter to total extinction + backscatter

      R2=f*b_tot/(a_tot+b_tot)
      Ocean_Reflec = R2

      !write(*,*)
      !write(*,*)Ocean_Reflec

!  continuation point

 60  continue

!  Finish

      return
end subroutine Ocean_Reflectance

subroutine Water_Transmittance_Quads &
    ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
      PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Output
      AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & 
      TRANS_NORM )

   IMPLICIT NONE
   INTEGER  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   integer, intent(in)   :: Max_PolarQuads, Max_AzimQuads

!  Output Quadratures
!  ------------------

   real(fpk), intent(out) :: PolarQuads    (Max_PolarQuads)   ! Radians
   real(fpk), intent(out) :: CosPolarQuads (Max_PolarQuads)
   real(fpk), intent(out) :: SinPolarQuads (Max_PolarQuads)
   real(fpk), intent(out) :: PolarWeights  (Max_PolarQuads)

   real(fpk), intent(out) :: AzimQuads    (Max_AzimQuads)     ! Radians
   real(fpk), intent(out) :: CosAzimQuads (Max_AzimQuads)
   real(fpk), intent(out) :: SinAzimQuads (Max_AzimQuads)
   real(fpk), intent(out) :: AzimWeights  (Max_AzimQuads)

!  Pre-computed Norm

   REAL(fpk), intent(out)  :: TRANS_NORM

!  Local

   integer   :: I, K
   real(fpk) :: zero, pi, pi2, pio2, weight

!  setups

   pi = acos(-1.0_fpk) ; pi2 = 2.0_fpk * pi ; pio2 = 0.5_fpk * pi
   zero = 0.0_fpk

!  Gaussian-quadrature calls

   CALL GETQUAD2 ( zero, pio2, Max_PolarQuads, PolarQuads, PolarWeights )
   CALL GETQUAD2 ( zero, pi2,  Max_AzimQuads,  AzimQuads,  AzimWeights  )

!  Develop ancillaries

   DO I = 1, Max_PolarQuads
      CosPolarQuads(i) = cos(PolarQuads(I))
      SinPolarQuads(i) = Sin(PolarQuads(I))
      PolarWeights (i) = PolarWeights(i) * CosPolarQuads(i) * SinPolarQuads(i)
   ENDDO
   do k = 1, Max_AzimQuads
      CosAzimQuads(k) = cos(AzimQuads(k))
      SinAzimQuads(k) = Sin(AzimQuads(k))
   ENDDO

!  Normalization

   trans_norm = zero
   do k = 1, Max_AzimQuads
      do i = 1, Max_PolarQuads
         weight = PolarWeights(i) * AzimWeights(k)
         trans_norm = trans_norm + weight
      enddo
   enddo

!  Done

   return
end subroutine Water_Transmittance_Quads

subroutine Water_Transmittance &
    ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
      PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
      AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
      DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,                     & ! Input
      INCIDENT_ANGLE, REFRAC_R, REFRAC_I,                     & ! Input
      WINDSPEED, SUNGLINT_COEFFS, PHI_W, CPHI_W, SPHI_W,      & ! Input
      TRANS_NORM, TRANS )

      IMPLICIT NONE
      INTEGER  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   integer, intent(in)   :: Max_PolarQuads, Max_AzimQuads

!  Quadratures

   real(fpk), intent(in) :: PolarQuads    (Max_PolarQuads)   ! Radians
   real(fpk), intent(in) :: CosPolarQuads (Max_PolarQuads)
   real(fpk), intent(in) :: SinPolarQuads (Max_PolarQuads)
   real(fpk), intent(in) :: PolarWeights  (Max_PolarQuads)

   real(fpk), intent(in) :: AzimQuads    (Max_AzimQuads)     ! Radians
   real(fpk), intent(in) :: CosAzimQuads (Max_AzimQuads)
   real(fpk), intent(in) :: SinAzimQuads (Max_AzimQuads)
   real(fpk), intent(in) :: AzimWeights  (Max_AzimQuads)

!  Windspeed, coefficients
!  -----------------------

!  Flag for Calculating Cox-Munk Coefficients
!     Only needs to be done once, so intent(inout)

      LOGICAL  , intent(inout) :: DO_COEFFS

!  Windspeed m/s

      REAL(fpk), intent(in)    :: WINDSPEED

!  Azimuth between Sun and Wind directions. angle in Radians + Cosine/sine

      REAL(fpk), intent(in)    :: PHI_W, CPHI_W, SPHI_W

!  Cox-Munk Coefficients. Intent(inout).

      REAL(fpk), intent(inout) :: SUNGLINT_COEFFS(7)

!  Other inputs
!  ------------

!  Flag for using Isotropic Facet distribution

      LOGICAL  , intent(in)    :: DO_ISOTROPIC

!  Flag for including Shadow effect

      LOGICAL  , intent(in)    :: DO_SHADOW

!  incident angle in degrees

      REAL(fpk), intent(in)    :: INCIDENT_ANGLE

!  Real and imaginary parts of refractive index

      REAL(fpk), intent(in)    :: REFRAC_R
      REAL(fpk), intent(in)    :: REFRAC_I

!  Pre-computed Norm

      REAL(fpk), intent(in)    :: TRANS_NORM

!  Output
!  ======

      REAL(fpk), intent(out)   :: TRANS

!  Local
!  =====

      integer   :: i, k
      real(fpk) :: dtr, xj, sxj, xi, sxi, phi, cphi, sphi, weight
      real(fpk) :: SUNGLINT_REFLEC

!  Computation setup

      TRANS = 0.0_fpk
      DTR   = ACOS(-1.0d0) / 180.0_fpk
      XJ  = COS ( INCIDENT_ANGLE * DTR )
      SXJ = SQRT ( 1.0_fpk - XJ * XJ )

!  Loops

      do k = 1, Max_AzimQuads
         PHI  = AzimQuads(K)/dtr
         CPHI = CosAzimQuads(K)
         SPHI = SinAzimQuads(K)
         do i = 1, Max_PolarQuads
            XI  = CosPolarQuads(I)
            SXI = SinPolarQuads(I)
            Call GENERAL_SUNGLINT &
             ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
               REFRAC_R, REFRAC_I, WINDSPEED,         &
               PHI_W, CPHI_W, SPHI_W,                 &
               XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
               SUNGLINT_COEFFS, SUNGLINT_REFLEC )
            WEIGHT = PolarWeights(I) * AzimWeights(k)
            TRANS = TRANS + SUNGLINT_REFLEC * WEIGHT
         enddo
      enddo
      TRANS = 1.0_fpk - (TRANS/TRANS_NORM)

!  done

      RETURN
end subroutine Water_Transmittance

!

subroutine Water_RefracIndex &
     ( Wavelength, Salinity, Refrac_R, Refrac_I )

!  THIS IS FORMERLY CALLED "INDWAT"

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input/Output
!  ------------

   real(fpk), intent(in)   :: Wavelength
   real(fpk), intent(in)   :: Salinity
   real(fpk), intent(out)  :: Refrac_R
   real(fpk), intent(out)  :: Refrac_I

!  Local
!  -----

!subroutine indwat(wl,xsal,nr,ni)
!
! input parameters:  wl=wavelength (in micrometers)
!                    xsal=salinity (in ppt), if xsal<0 then 34.3ppt by default
! output parameters: nr=index of refraction of sea water
!                    ni=extinction coefficient of sea water

       real twl(62),tnr(62),tni(62)
       real nr,ni,wl,xwl,yr,yi,nrc,nic,xsal
       integer i

! Indices of refraction for pure water from Hale and Querry, 
! Applied Optics, March 1973, Vol. 12,  No. 3, pp. 555-563
       data twl/&
        0.250,0.275,0.300,0.325,0.345,0.375,0.400,0.425,0.445,0.475,&
        0.500,0.525,0.550,0.575,0.600,0.625,0.650,0.675,0.700,0.725,&
        0.750,0.775,0.800,0.825,0.850,0.875,0.900,0.925,0.950,0.975,&
        1.000,1.200,1.400,1.600,1.800,2.000,2.200,2.400,2.600,2.650,&
        2.700,2.750,2.800,2.850,2.900,2.950,3.000,3.050,3.100,3.150,&
        3.200,3.250,3.300,3.350,3.400,3.450,3.500,3.600,3.700,3.800,&
        3.900,4.000/
        data tnr/&
        1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,1.336,&
        1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.331,1.330,&
        1.330,1.330,1.329,1.329,1.329,1.328,1.328,1.328,1.327,1.327,&
        1.327,1.324,1.321,1.317,1.312,1.306,1.296,1.279,1.242,1.219,&
        1.188,1.157,1.142,1.149,1.201,1.292,1.371,1.426,1.467,1.483,&
        1.478,1.467,1.450,1.432,1.420,1.410,1.400,1.385,1.374,1.364,&
        1.357,1.351/
        data tni/&
        3.35E-08,2.35E-08,1.60E-08,1.08E-08,6.50E-09,&
        3.50E-09,1.86E-09,1.30E-09,1.02E-09,9.35E-10,&
        1.00E-09,1.32E-09,1.96E-09,3.60E-09,1.09E-08,&
        1.39E-08,1.64E-08,2.23E-08,3.35E-08,9.15E-08,&
        1.56E-07,1.48E-07,1.25E-07,1.82E-07,2.93E-07,&
        3.91E-07,4.86E-07,1.06E-06,2.93E-06,3.48E-06,&
        2.89E-06,9.89E-06,1.38E-04,8.55E-05,1.15E-04,&
        1.10E-03,2.89E-04,9.56E-04,3.17E-03,6.70E-03,&
        1.90E-02,5.90E-02,1.15E-01,1.85E-01,2.68E-01,&
        2.98E-01,2.72E-01,2.40E-01,1.92E-01,1.35E-01,&
        9.24E-02,6.10E-02,3.68E-02,2.61E-02,1.95E-02,&
        1.32E-02,9.40E-03,5.15E-03,3.60E-03,3.40E-03,&
        3.80E-03,4.60E-03/

!  Assign input

      wl   = real(WAVELENGTH)
      xsal = real(SALINITY)
      Refrac_R = 0.0_fpk
      Refrac_I = 0.0_fpk

!  Find wavelength point for interpolation

        i=2
 10     if (wl.lt.twl(i)) goto 20
        if (i.lt.62) then
           i=i+1
           goto 10
         endif

!  Interpolate

 20     xwl=twl(i)-twl(i-1)
        yr=tnr(i)-tnr(i-1)
        yi=tni(i)-tni(i-1)
        nr=tnr(i-1)+(wl-twl(i-1))*yr/xwl
        ni=tni(i-1)+(wl-twl(i-1))*yi/xwl
!
! Correction to be applied to the index of refraction and to the extinction 
! coefficients of the pure water to obtain the ocean water one (see for 
! example Friedman). By default, a typical sea water is assumed 
! (Salinity=34.3ppt, Chlorinity=19ppt) as reported by Sverdrup. 
! In that case there is no correction for the extinction coefficient between 
! 0.25 and 4 microns. For the index of refraction, a correction of +0.006 
! has to be applied (McLellan). For a chlorinity of 19.0ppt the correction 
! is a linear function of the salt concentration. Then, in 6S users are able 
! to enter the salt concentration (in ppt).
! REFERENCES:
! Friedman D., Applied Optics, 1969, Vol.8, No.10, pp.2073-2078.
! McLellan H.J., Elements of physical Oceanography, Pergamon Press, Inc.,
!        New-York, 1965, p 129.
! Sverdrup H.V. et al., The Oceans (Prentice-Hall, Inc., Englewood Cliffs,
!        N.J., 1942, p 173.

        nrc=0.006
        nic=0.000
        nr=nr+nrc*(xsal/34.3)
        ni=ni+nic*(xsal/34.3)

!  Assign output

    REFRAC_R = real(nr,fpk)
    REFRAC_I = real(ni,fpk)

    return
end subroutine Water_RefracIndex

!

      SUBROUTINE GENERAL_SUNGLINT &
         ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
           REFRAC_R, REFRAC_I, WINDSPEED,         &
           PHI_W, CPHI_W, SPHI_W,                 &
           XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
           SUNGLINT_COEFFS, SUNGLINT_REFLEC )

      implicit none

      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Subroutine Input arguments
!  --------------------------

!  Flag for using Isotropic Facet distribution

      LOGICAL  , intent(in)    :: DO_ISOTROPIC

!  Flag for including Shadow effect

      LOGICAL  , intent(in)    :: DO_SHADOW

!  Flag for Calculating Cox-Munk Coefficients
!     Only needs to be done once, so intent(inout)

      LOGICAL  , intent(inout) :: DO_COEFFS

!  Real and imaginary parts of refractive index

      REAL(fpk), intent(in)    :: REFRAC_R
      REAL(fpk), intent(in)    :: REFRAC_I

!  Windspeed m/s

      REAL(fpk), intent(in)    :: WINDSPEED

!  Azimuth between Sun and Wind directions. angle in Radians + Cosine/sine

      REAL(fpk), intent(in)    :: PHI_W, CPHI_W, SPHI_W

!  Incident and reflected ddirections: sines/cosines. Relative azimuth (angle in radians)

      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SPHI

!  Subroutine output arguments
!  ---------------------------

!   Glitter reflectance

      REAL(fpk), intent(out)   :: SUNGLINT_REFLEC

!  Cox-Munk Coefficients. Intent(inout).

      REAL(fpk), intent(inout) :: SUNGLINT_COEFFS(7)

!  Local arguments
!  ---------------

!  parameters from LIDORT

   real(fpk), PARAMETER :: ONE = 1.0_fpk, ZERO  = 0.0_fpk, ONEP5 = 1.5_fpk
   real(fpk), PARAMETER :: TWO = 2.0_fpk, THREE = 3.0_fpk, FOUR  = 4.0_fpk
   real(fpk), PARAMETER :: six = two * three, twentyfour = six * four
   real(fpk), PARAMETER :: QUARTER = 0.25_fpk, HALF = 0.5_fpk

   real(fpk), PARAMETER :: MINUS_ONE = - ONE
   real(fpk), PARAMETER :: MINUS_TWO = - TWO

   real(fpk), PARAMETER :: PIE = ACOS(MINUS_ONE)
   real(fpk), PARAMETER :: DEG_TO_RAD = PIE/180.0_fpk

   real(fpk), PARAMETER :: PI2  = TWO  * PIE
   real(fpk), PARAMETER :: PI4  = FOUR * PIE
   real(fpk), PARAMETER :: PIO2 = HALF * PIE
   real(fpk), PARAMETER :: PIO4 = QUARTER * PIE

   REAL(fpk), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      REAL(fpk)  :: B, ZX, ZY, Z, Z1, Z2, XMP
      REAL(fpk)  :: TILT, TANTILT, TANTILT_SQ, COSTILT
      REAL(fpk)  :: ARGUMENT, PROB, FAC2, COEFF, VAR, WSigC, WSigU
      REAL(fpk)  :: XE, XN, XE_sq, XN_sq, XE_sq_1, XN_sq_1
      REAL(fpk)  :: XPHI, CKPHI, SKPHI, XPHI_W, CKPHI_W, SKPHI_W
      REAL(fpk)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(fpk)  :: SHADOWI, SHADOWR, SHADOW

!  Initialise output

      SUNGLINT_REFLEC = ZERO

!  COmpute coefficients, according to 6S formulation

      IF ( DO_COEFFS ) THEN
         SUNGLINT_COEFFS = zero
         IF ( DO_ISOTROPIC ) THEN
            SUNGLINT_COEFFS(1) = 0.003_fpk + 0.00512_fpk * WINDSPEED
         ELSE
            SUNGLINT_COEFFS(1) = 0.003_fpk + 0.00192_fpk * WINDSPEED ! sigmaC
            SUNGLINT_COEFFS(2) =             0.00316_fpk * WINDSPEED ! sigmaU
            SUNGLINT_COEFFS(3) = 0.010_fpk - 0.00860_fpk * WINDSPEED ! C21
            SUNGLINT_COEFFS(4) = 0.040_fpk - 0.03300_fpk * WINDSPEED ! C03
            SUNGLINT_COEFFS(5) = 0.400_fpk                           ! C40
            SUNGLINT_COEFFS(6) = 0.230_fpk                           ! C04
            SUNGLINT_COEFFS(7) = 0.120_fpk                           ! C22
         ENDIF
         DO_COEFFS = .false.
      ENDIF

!  Local angles

      XPHI   = PIE - PHI       ! Not used
!     CKPHI  = - CPHI          ! Original, not correct.

      CKPHI  = + CPHI
      SKPHI  = + SPHI

      XPHI_W  = PHI_W
      CKPHI_W = CPHI_W
      SKPHI_W = SPHI_W

!  Tilt angle

      B  = ONE / ( XI + XJ )
      ZX = - SXI * SKPHI * B
      ZY = ( SXJ + SXI * CKPHI ) * B
      TANTILT_SQ = ZX * ZX + ZY * ZY
      TANTILT    = SQRT ( TANTILT_SQ )
      TILT       = ATAN(TANTILT)
      COSTILT    = COS(TILT)

!  Scatter angle

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)

!  Fresnel
!  -------

       CALL Fresnel_Complex ( REFRAC_R, REFRAC_I, Z2, XMP )

!  Anisotropic
!  -----------

      IF ( .not. DO_ISOTROPIC ) THEN

!  Variance

         WSigC = ONE / Sqrt(SUNGLINT_COEFFS(1))
         WSigU = ONE / Sqrt(SUNGLINT_COEFFS(2))
         VAR   = WSigC * WSigU * HALF ; VAR = ONE / VAR

!  angles

         XE = (  CKPHI_W * ZX + SKPHI_W * ZY ) * WSigC ; XE_sq = XE * XE ; XE_sq_1 = xe_sq - one
         XN = ( -SKPHI_W * ZX + CKPHI_W * ZY ) * WSigU ; XN_sq = XN * XN ; XN_sq_1 = xn_sq - one

!  GC Coefficient

         Coeff = ONE - SUNGLINT_COEFFS(3) *      XE_sq_1      * XN * half &
                     - SUNGLINT_COEFFS(4) * ( XN_sq - three ) * XN / six  &
                     + SUNGLINT_COEFFS(5) * ( XE_sq * XE_sq - six * XE_sq + three ) / twentyfour &
                     + SUNGLINT_COEFFS(6) * ( XN_sq * XN_sq - six * XN_sq + three ) / twentyfour &
                     + SUNGLINT_COEFFS(7) * XE_sq_1 * XN_sq_1 / four

!  Probability and finish

         ARGUMENT = ( XE_sq  + XN_sq ) * HALF
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            PROB = COEFF * EXP ( - ARGUMENT ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC = XMP * PROB * FAC2
         ENDIF

      ENDIF

!  Isotropic
!  ---------

      IF ( DO_ISOTROPIC ) THEN

!  Compute Probability and finish

         COEFF = ONE
         VAR   = SUNGLINT_COEFFS(1)
         ARGUMENT = TANTILT_SQ / VAR
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            PROB = COEFF * EXP ( - ARGUMENT ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC = XMP * PROB * FAC2
         ENDIF

      ENDIF

!  No Shadow code if not flagged

      IF ( .not. DO_SHADOW  ) RETURN

!  Shadow code

      S1 = SQRT ( VAR / PIE )
      S3 = ONE / ( SQRT(VAR) )
      S2 = S3 * S3

      XXI  = XI*XI
      DCOT = XI / SQRT ( ONE - XXI )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function
      SHADOWI = HALF * ( S1 * T1 / DCOT - T2 )

      XXJ  = XJ*XJ
      DCOT = XJ / SQRT ( ONE - XXJ )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function
      SHADOWR = HALF * ( S1 * T1 / DCOT - T2 )

      SHADOW = ONE / ( ONE + SHADOWI + SHADOWR )
      SUNGLINT_REFLEC = SUNGLINT_REFLEC * SHADOW

!     Finish

      RETURN
      END SUBROUTINE GENERAL_SUNGLINT

!

subroutine get_fluorescence_755(lat, lon, epoch, sza, fluo_file, Fs755)

!  Function from Chris O'Dell, 12 July 12
!    Minor reprogramming for the SL supplement by R. Spurr, 12 July 2012

! This subroutine calculates the fluorescence intensity at 755 nm in W/m^2/um/sr,
! as a function of day of year (given via "epoch"), lat, lon, and solar zenith angle (sza).

      implicit none

!  I/O

      double precision, intent(in)       :: lat, lon  ! latitude & longitude of desired location
      integer, dimension(:), intent(in)  :: epoch     ! 6-7 element array with Year/month/day/hour/min/sec/msec
      double precision, intent(in)       :: sza       ! Solar zenith angle in degrees
      character(LEN=*), intent(in)       :: fluo_file ! file containing fluorescence climatology
      double precision, intent(out)      :: Fs755     ! fluorescence at 755 nm in W/m^2/um/sr

!  local variables

      integer, parameter :: NLAT_FLUO_FILE = 291
      integer, parameter :: NLON_FLUO_FILE = 720
      double precision, dimension(NLAT_FLUO_FILE) :: lat_data
      double precision, dimension(NLON_FLUO_FILE) :: lon_data

      real, SAVE, dimension(NLON_FLUO_FILE, NLAT_FLUO_FILE, 12) :: fluo_data
      logical, SAVE :: Fluor_Data_Loaded=.FALSE.

      integer, parameter, dimension(12) :: DAYS_IN_MONTH = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

      integer :: i, j, mon1, mon2, loc1(1)
      real :: scene_lon, fmon
      double precision :: doy, Fs_corr,  avmu, pi, d2r

      integer :: FUNIT, ios, lunit
      logical :: Assigned, verbose

!New variables
      integer :: k,kmax,mon
      real, dimension(NLON_FLUO_FILE*NLAT_FLUO_FILE*12) :: fluo_data_temp

      logical :: use_nag_compiler=.false.
      !logical :: use_nag_compiler=.true.

!  initialize

      Fs755 = 0.0d0
      PI = acos(-1.0d0) ; D2R = PI / 180.0d0
      verbose = .false. ; lunit = 4555

!  Grids

      do i = 1, NLAT_FLUO_FILE
          lat_data(i) = 90. - (i-1)*0.5
      enddo
      do j = 1, NLON_FLUO_FILE
          lon_data(j) = -180.0 + 0.5 * (j-1)
      enddo

      if (.NOT. Fluor_Data_Loaded) then
!       Select the next available unit number.
        FUNIT=1
        INQUIRE(UNIT=FUNIT,OPENED=Assigned)
        DO WHILE (Assigned)
           FUNIT=FUNIT+1
           INQUIRE(UNIT=FUNIT,OPENED=Assigned)
        END DO
        if (verbose) write(lunit,*) 'Fluo File = ' // trim(fluo_file)
        open(FUNIT, file=trim(fluo_file),&
             form='UNFORMATTED', status='OLD', IOSTAT=ios)

        if (ios /=0) then
           print *, 'Error Opening Fluorescence file ' // trim(fluo_file)
           STOP
        endif

        if (.not. use_nag_compiler) then
          !original read
          read(FUNIT, IOSTAT=ios) fluo_data
        else
          !modified read section
          read(FUNIT, IOSTAT=ios) fluo_data_temp

          !prepare to read from array "fluo_data_temp"
          !starting at position 2 (not 1!) since NAG reads
          !the binary file record header as a data point
          k=1
          kmax=NLON_FLUO_FILE*NLAT_FLUO_FILE*12
          do mon=1,12
            do i=1,NLAT_FLUO_FILE
              do j=1,NLON_FLUO_FILE
                k = k + 1
                if (k <= kmax) then
                  fluo_data(j,i,mon) = fluo_data_temp(k)
                else
                  fluo_data(j,i,mon) = 0.0
                end if
              enddo
            enddo
          enddo
        endif

        if (ios /=0) then
           print *, 'Error Reading Fluorescence file ' // trim(fluo_file)
           STOP
        endif
        close(FUNIT)
        Fluor_Data_Loaded = .TRUE.
      endif

      ! find closest lat and lon points
      scene_lon = real(lon)
      if (scene_lon >= 179.75) scene_lon = scene_lon - 360.0
      loc1 = minloc( abs(scene_lon-lon_data) )
      j = loc1(1)
      loc1 = minloc(abs(lat-lat_data))
      i = loc1(1)

      ! Do an interpolation in month.  Assume data file contains month days in middle of each month
      mon1 = epoch(2)
      ! this quantity is 0.5 in the middle of the month
      fmon = (epoch(3)-0.5) / days_in_month(mon1)
      ! This quantity is 0.0 in the middle of the month,
      !                 -0.5 at the beginning of the month, and
      !                 +0.5 at the end
      fmon = fmon - 0.5
      if (fmon < 0.) then
         mon1 = mon1-1
         fmon = fmon + 1.
      endif
      mon2 = mon1 + 1
      if (mon1==0) mon1=12
      if (mon2==13) mon2=1

      if (mon1<1 .OR. mon1 >12 .OR. &
          mon2<1 .OR. mon2>12 .OR. &
          i<1 .OR. i>NLAT_FLUO_FILE .OR. &
          j<1 .OR. j>NLON_FLUO_FILE) then
          print *, 'BAD PIXEL ATTEMPT IN FLUORESCENCE MODULE.'
          write(*, "('mon1=', i3, '; mon2=', i3, '; i=', i5, '; j=', i5)") mon1, mon2, i ,j
          write(*, "('epoch = ', 7i8)") epoch
          write(*, "('year, month, day = ', 3i8)") epoch(1:3)
          write(*, "('hr, min, sec = ', 3i8)") epoch(4:6)
          write(*, "('Lat, Lon, Sza = ', 3f12.4)") lat, lon, sza
      endif

      ! Calculate daily-averaged 755nm Fluorescence for this location & DOY.
      fs_corr = fluo_data(j,i,mon1) * (1.-fmon) + fluo_data(j,i,mon2) * fmon
      ! (NOTE: Fs_Corr currently has units of W/m2/um/sr)

      ! Now, convert from daily average Fluoresence to instantaneous fluorescence
      doy = epoch(3) + (epoch(4) + epoch(5)/60. + epoch(6)/3600.)/24.
      do i = 1, epoch(2)-1
        doy = doy + days_in_month(i)
      enddo
      avmu = average_solar_cosine(lat, doy, pi, d2r)

      Fs755 = Fs_corr / avmu * cos(sza * D2R)

      END subroutine get_fluorescence_755

!

      FUNCTION average_solar_cosine(lat, doy, pi, d2r) result(avmu)
        implicit none
        double precision, INTENT(IN) :: lat, doy, pi, d2r
        double precision             :: avmu

        real, dimension(0:3), parameter :: cn = (/ 0.006918, -0.399912, -0.006758, -0.002697 /)
        real, dimension(0:3), parameter :: dn = (/ 0., 0.070257, 0.000907, 0.000148 /)

        double precision :: dec, t, H, cH
        integer :: i

        t = 2.d0*PI*(doy-1.d0)/365.d0
        ! solar declination in radians
        dec = 0.d0
        do i = 0,3
           dec = dec + cn(i) * cos(i*t) + dn(i)*sin(i*t)
        enddo
        cH = - tan(lat*D2R) * tan(dec)

        ! H = length of solar half-day in hours
        if (cH .LT. 1.0) then ! there is some sun
            if (cH .LT. -1.0) then
                ! sun is always up
                H = PI
            else
                ! sun rises and sets like a normal place
                H = abs(acos(cH))
            endif
        else
            ! there is no sun at all
            H = 0.d0
        endif

        avmu = 1.d0/PI * (sin(lat*D2R)*sin(dec)*H + cos(lat*D2R)*cos(dec)*sin(H))

      END FUNCTION average_solar_cosine

!

      FUNCTION solar_spec_irradiance(wavelength) result(ssi)

!Reads a solar spectral irradiance file and returns the solar
!spectral irradiance for an air mass of zero in
!units of W m^-2 μm^-1 for the wavelength input

!Function by Mick Christi - 16 July 2012

      implicit none

      !Input variable
      double precision, intent(in)   :: wavelength !in um

      !Output variable
      double precision               :: ssi

      !Local variables

      !Number of data for solar data arrays
      integer, parameter             :: maxfiledata = 24000

      !Regular help variables
      integer                        :: i,numfiledata,solar_file,obs_period
      double precision               :: w,slope,data_src,norm_factor,&
                                        solar_integ_irad_in
      double precision, dimension(3) :: solar_spec_irad_in
      logical                        :: normalize

      !Saved help variables
      double precision, dimension(maxfiledata) :: &
                                        wvl = -1.0d0,&
                                        solar_spec_irad = -1.0d0

!Start program

!Obtain solar data if necessary

      !solar_file = 1  1985 Wehrli Standard Extraterrestrial Spectral
      !                Solar Irradiance
      !                (TSI = 1367 W/m2)
      !solar_file = 2  Solar Spectral Irradiance Reference Spectra for
      !                Whole Heliosphere Interval (WHI) 2008
      !                (TSI = ?; however, may be normalized by user)
      solar_file = 2

      if (solar_file == 1) then
        numfiledata = 920
      else if (solar_file == 2) then
        numfiledata = 24000
      endif

      if (wvl(1) < 0.0d0) then
        if (solar_file == 1) then
          !1985 Wehrli Standard Extraterrestrial Spectral Solar Irradiance file
          open(unit=50,file='lrrs_test/fluorescence_data/wehrli85.dat',&
               status='old',action='read')

          !Skip file header
          do i=1,5
            read(50,*)
          enddo

          !Read solar data and convert:
          !(1) wavelength from nm to um
          !(2) spectral irradiance data from W m^-2 nm^-1 to W m^-2 μm^-1
          do i=1,numfiledata
            read(50,*) wvl(i),solar_spec_irad_in(1),solar_integ_irad_in
            wvl(i) = wvl(i)*1.0d-3
            solar_spec_irad(i) = solar_spec_irad_in(1)*1.0d3
          enddo
          close(50)
        else if (solar_file == 2) then
          !Solar Spectral Irradiance Reference Spectra for
          !Whole Heliosphere Interval (WHI) 2008 file
          open(unit=50,file='lrrs_test/fluorescence_data/ref_solar_irradiance_whi-2008_ver2.dat',&
               status='old',action='read')

          !Skip file header
          do i=1,144
            read(50,*)
          enddo

          !Note: observation period for this data set -
          ! obs_period = 1  Moderately low solar activity with sunspot darkening.
          !                 TSI = 1360.696 W/m^2
          ! obs_period = 2  Moderately low solar activity with faculae brightening.
          !                 TSI = 1360.944 W/m^2
          ! obs_period = 3  Very close to solar cycle minimum condition.
          !                 TSI = 1360.840 W/m^2
          obs_period = 3

          !Create normalization factor to adjust radiances to correspond
          !to a TSI of 1366.1 W/m^2
          if (obs_period == 1) then
            norm_factor = 1366.1d0/1360.696d0
          else if (obs_period == 2) then
            norm_factor = 1366.1d0/1360.944d0
          else if (obs_period == 3) then
            norm_factor = 1366.1d0/1360.840d0
          end if
          !write(*,*) 'norm factor = ',norm_factor
          !read(*,*)

          !Read solar data and convert:
          !(1) wavelength from nm to um
          !(2) spectral irradiance data from W m^-2 nm^-1 to W m^-2 μm^-1
          do i=1,numfiledata
            read(50,'(F8.2,3E12.4,F4.0)') wvl(i),solar_spec_irad_in(1:3),data_src
            wvl(i) = wvl(i)*1.0d-3
            solar_spec_irad(i) = solar_spec_irad_in(obs_period)*1.0d3

            !Actually use normalization factor if desired
            normalize = .true. !.false.
            if (normalize) then
              solar_spec_irad(i) = solar_spec_irad(i)*norm_factor
            end if
          enddo
          close(50)
        endif
      end if

!Find indices of wavelengths in spectral irradiance data surrounding
!the input wavelength

      w = wavelength

      !Handle special cases
      if (w < wvl(1)) then
        write(*,*) 'Error in FUNCTION solar_spec_irradiance: input wavelength ' // &
                   'is less than minimum wavelength in solar data file'
        write(*,*) 'input wavelength is: ',w
        write(*,*) 'min   wavelength is: ',wvl(1)
        stop
      end if
      if (w > wvl(numfiledata)) then
        write(*,*) 'Error in FUNCTION solar_spec_irradiance: input wavelength ' // &
                   'is greater than maximum wavelength in solar file'
        write(*,*) 'input wavelength is: ',w
        write(*,*) 'max   wavelength is: ',wvl(numfiledata)
        stop
      end if

      i = 1
      do
        if ( (w >= wvl(i)) .and. (w < wvl(i+1)) ) exit
        i = i + 1
      end do

!Linearly interpolate solar spectral irradiance to the input wavelength

      slope = (solar_spec_irad(i+1) - solar_spec_irad(i)) / &
              (wvl(i+1) - wvl(i))
      ssi = slope*(w - wvl(i)) + solar_spec_irad(i)

      END FUNCTION solar_spec_irradiance

!  End module

      END MODULE sleave_sup_routines_m

