module GEMSTOOL_NSW_Fluorescence_m

!  Rob, created 10/18/16. Only designed to work with NSW tool.
!       updated 11/30/16. Optional use of linear parameterization

! This Module is designed to calculats the sun-normalized fluorescence radiance 
! as a function of the time, latitude, longitude, and solar zenith angle.

!  The 755-Amplitude Scaling derivative is a trivial extra piece of output
!     GEMSTOOL_SIF          = SIFScaling * SLEAVE_output
!     dGEMSTOOLSIF_dScaling = SLEAVE_output (non-normalized)

! Routine calls the VLIDORT SLEAVE supplement, operating with Fluorescence option
!  -- inputs are provided by variables in the GEMSTOOL Inputs

!  Six variables are required for the Time:
!     Year/month/day/hour/min/sec/

use GEMSTOOL_Input_Types_m

!  Rob, add the SLEAVE type structures 

USE VSLEAVE_SUP_MOD

public
private fpk

contains

subroutine GEMSTOOL_Fluorescence &
      ( do_wavnums, Wav, GEMSTOOL_INPUTS, SIF_datapath, &
        GEMSTOOL_SLEAVE_Results )

! Formerly  ---> Gemstool_SIF, dGEMSTOOLSIF_dScaling )
     
   implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Flag for using wavenumbers

   logical  , intent(in) :: do_wavnums

!  Wavelength/wavenumber input

   real(fpk), intent(in) :: Wav

! GEOMTOOL inputs, structure. Intent(in) here

   TYPE(GEMSTOOL_Config_Inputs), INTENT(IN) :: GEMSTOOL_INPUTS

!  Path to the Fluorescence data. set in the Wrapper.

   character*(*), intent(in) :: SIF_datapath

!  VSLEAVE supplement output structure, added 10/25/16

   TYPE(VSLEAVE_Sup_Outputs), intent(inout)              :: GEMSTOOL_SLEAVE_Results

!  SIF output
!   real(fpk), intent(out) :: Gemstool_SIF          (Maxwav,MaxGeoms)
!   real(fpk), intent(out) :: dGEMSTOOLSIF_dScaling (Maxwav,MaxGeoms)

!  Local
!  =====

!  VSLEAVE supplement structures (from VLIDORT)
!  --------------------------------------------

!  VSLEAVE supplement input structures

   TYPE(VSLEAVE_Sup_Inputs)               :: VSLEAVE_Sup_In

!  other variables

   integer   :: g, ngeoms
   real(fpk), parameter :: fzero = 0.0_fpk

!  initialize output
!   Gemstool_SIF          = fzero
!   dGEMSTOOLSIF_dScaling = fzero

!  proxy

   ngeoms = GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries

!  General inputs
!  --------------

!  Following flags ar always required

   VSLEAVE_Sup_In%SL_DO_SOLAR_SOURCES = .TRUE.
   VSLEAVE_Sup_In%SL_DO_USER_OBSGEOMS = .TRUE.
   VSLEAVE_Sup_In%SL_DO_USER_STREAMS  = .TRUE.

!  Stokes and streams

   VSLEAVE_Sup_In%SL_NSTREAMS = GEMSTOOL_INPUTS%RTControl%NVlidort_nstreams
   VSLEAVE_Sup_In%SL_NSTOKES  = GEMSTOOL_INPUTS%RTControl%NVlidort_nstokes

!  Geometry - copy GEMSTOOL inputs
!  --------

   VSLEAVE_Sup_In%SL_NBEAMS          = ngeoms
   VSLEAVE_Sup_In%SL_N_USER_STREAMS  = ngeoms
   VSLEAVE_Sup_In%SL_N_USER_RELAZMS  = ngeoms
   VSLEAVE_Sup_In%SL_N_USER_OBSGEOMS = ngeoms
   do g =  1, ngeoms
      VSLEAVE_Sup_In%SL_BEAM_SZAS(g)         = GEMSTOOL_INPUTS%Geometry%GEMS_szas(g)
      VSLEAVE_Sup_In%SL_USER_ANGLES_INPUT(g) = GEMSTOOL_INPUTS%Geometry%GEMS_vzas(g)
      VSLEAVE_Sup_In%SL_USER_RELAZMS(g)      = GEMSTOOL_INPUTS%Geometry%GEMS_azms(g)
      VSLEAVE_Sup_In%SL_USER_OBSGEOMS(g,1)   = GEMSTOOL_INPUTS%Geometry%GEMS_szas(g)
      VSLEAVE_Sup_In%SL_USER_OBSGEOMS(g,2)   = GEMSTOOL_INPUTS%Geometry%GEMS_vzas(g)
      VSLEAVE_Sup_In%SL_USER_OBSGEOMS(g,3)   = GEMSTOOL_INPUTS%Geometry%GEMS_azms(g)
   ENDDO
!mick fix 1/24/2017 - initialize these elements
   VSLEAVE_Sup_In%SL_USER_OBSGEOMS(ngeoms+1:,:) = FZERO

!  Water-leaving variables.
!  -----------------------

!  Not required, just initialized

   VSLEAVE_Sup_In%SL_SALINITY   = FZERO
   VSLEAVE_Sup_In%SL_CHLORCONC  = FZERO
   VSLEAVE_Sup_In%SL_WAVELENGTH = FZERO
   VSLEAVE_Sup_In%SL_WINDSPEED  = FZERO
   VSLEAVE_Sup_In%SL_WINDDIR    = FZERO
   VSLEAVE_Sup_In%SL_DO_GlintShadow   = .FALSE.
   VSLEAVE_Sup_In%SL_DO_FoamOption    = .FALSE.
   VSLEAVE_Sup_In%SL_DO_FacetIsotropy = .FALSE.

!  Fluorescence variables
!  ----------------------

!  Set Fluorescence control
!    11/30/16. Add control for Exact-only approximation, linear parameterization

   VSLEAVE_Sup_In%SL_DO_ISOTROPIC    = .TRUE.
   VSLEAVE_Sup_In%SL_DO_SLEAVING     = .TRUE.
   VSLEAVE_Sup_In%SL_DO_FLUORESCENCE = .TRUE.
   VSLEAVE_Sup_In%SL_DO_EXACT        = .TRUE.
   VSLEAVE_Sup_In%SL_DO_EXACTONLY      = GEMSTOOL_INPUTS%Sleave%DO_SIF_ExactOnly
   VSLEAVE_Sup_In%SL_FL_DO_LinearParam = GEMSTOOL_INPUTS%Sleave%DO_SIF_LinearParam

!  Gaussians (automatic here, not used if lineari parameterization is in force)

   VSLEAVE_Sup_In%SL_FL_DO_DataGaussian  = .true.
   VSLEAVE_Sup_In%SL_FL_InputGAUSSIANS   = fzero  ! Will be set internally

!  Zero the GEOMSTOOL dependencies

   VSLEAVE_Sup_In%SL_FL_ScalingConstant    = 0.0d0
   VSLEAVE_Sup_In%SL_FL_ScalingGradient    = 0.0d0
   VSLEAVE_Sup_In%SL_FL_Amplitude755       = 0.0d0

!  Lat/Lon and Epoch (Year/month/day/hour/min/sec)

   VSLEAVE_Sup_In%SL_FL_LATITUDE   = GEMSTOOL_INPUTS%TimePos%latitude
   VSLEAVE_Sup_In%SL_FL_LONGITUDE  = GEMSTOOL_INPUTS%TimePos%longitude

   VSLEAVE_Sup_In%SL_FL_EPOCH(1)      = GEMSTOOL_INPUTS%TimePos%Year
   VSLEAVE_Sup_In%SL_FL_EPOCH(2)      = GEMSTOOL_INPUTS%TimePos%Month
   VSLEAVE_Sup_In%SL_FL_EPOCH(3)      = GEMSTOOL_INPUTS%TimePos%Day_of_Month
   VSLEAVE_Sup_In%SL_FL_EPOCH(4)      = GEMSTOOL_INPUTS%TimePos%Hour
   VSLEAVE_Sup_In%SL_FL_EPOCH(5)      = GEMSTOOL_INPUTS%TimePos%Minute
   VSLEAVE_Sup_In%SL_FL_EPOCH(6)      = GEMSTOOL_INPUTS%TimePos%Second

!  Wavelength setting, MUST be [nm]

   if ( do_wavnums ) then
      VSLEAVE_Sup_In%SL_FL_WAVELENGTH = 1.0d7/ wav
   else
      VSLEAVE_Sup_In%SL_FL_WAVELENGTH = wav
   endif

!   11/30/16. Add variables for Linear and gaussian parameterizations
!             - Linear parameterization only for the basic SIF_755 value

   if ( VSLEAVE_Sup_In%SL_FL_DO_LinearParam ) then
!      VSLEAVE_Sup_In%SL_FL_ScalingConstant    = GEMSTOOL_INPUTS%Sleave%SIF755_Scaling_Constant
!      VSLEAVE_Sup_In%SL_FL_ScalingGradient    = GEMSTOOL_INPUTS%Sleave%SIF755_Scaling_Gradient
      VSLEAVE_Sup_In%SL_FL_ScalingConstant    = 1.0d0
      VSLEAVE_Sup_In%SL_FL_ScalingGradient    = 0.0d0
   else
      VSLEAVE_Sup_In%SL_FL_Amplitude755     = GEMSTOOL_INPUTS%Sleave%SIF755_Amplitude
   endif

!  Regular VSLEAVE supplement call -
!    the output will be used for the VLIDORT calculation

   CALL VSLEAVE_MAINMASTER ( SIF_datapath, &
        VSLEAVE_Sup_In,         & ! Inputs
        GEMSTOOL_SLEAVE_Results ) ! Outputs

!  Former wavelength loop
!  ----------------------

!   do w = 1, nwav
!      if ( do_wavnums ) then
!         VSLEAVE_Sup_In%SL_FL_WAVELENGTH = 1.0d7/ wav(w)
!      else
!         VSLEAVE_Sup_In%SL_FL_WAVELENGTH = wav(w)
!      endif
!      CALL VSLEAVE_MAINMASTER ( SIF_datapath, &
!        VSLEAVE_Sup_In,            & ! Inputs
!        GEMSTOOL_SLEAVE_Results )    ! Outputs
!      do g = 1, ngeoms
!         Gemstool_SIF(w,g)  = VSGEMSTOOL_SLEAVE_Results%SL_SLTERM_ISOTROPIC(1,g)
!         if ( GEMSTOOL_INPUTS%LinControl%do_SIFScaling_Jacobian ) then
!            dGEMSTOOLSIF_dScaling(w,g)  = Gemstool_SIF(w,g)/GEMSTOOL_INPUTS%Sleave%SIF755_Amplitude
!         endif
!      enddo
!   enddo

!  Done

   return
end subroutine GEMSTOOL_Fluorescence

!  End module

end module GEMSTOOL_NSW_Fluorescence_m
