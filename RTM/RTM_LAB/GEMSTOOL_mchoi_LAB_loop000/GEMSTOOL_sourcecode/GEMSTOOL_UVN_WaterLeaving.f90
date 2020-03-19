module GEMSTOOL_UVN_WaterLeaving_m

!  Rob, created 10/25/16. This is based on the VLIDORT combination of
!       Water-leaving with Ocean Glitter (the "NewCM" kernel

!  Only designed to work with the UVN tool.

!  We will be calling both the VSLEAVE and VBRDF supplements in sequence,
!   and for each wavelength

! Routine calls the VLIDORT VBRDF supplement, operating with Kernel choices
!  -- inputs are provided by variables in the GEMSTOOL Inputs

use GEMSTOOL_Input_Types_m

!  Add the BRDF and SLEAVE type structures (contains everything you need)

USE VBRDF_SUP_MOD
USE VSLEAVE_SUP_MOD

public
private fpk

contains

subroutine GEMSTOOL_WaterLeaving &
      ( GEMSTOOL_INPUTS, WaveLength,                    & ! inputs
        GEMSTOOL_BRDF_Results, GEMSTOOL_SLEAVE_Results, & ! outputs
        fail, nmessages, messages )
     
   implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

! GEOMTOOL inputs, structure. Intent(in) here

   TYPE(GEMSTOOL_Config_Inputs), INTENT(IN) :: GEMSTOOL_INPUTS

!  Wavelength input in [nm]

   real(fpk), intent(in) :: WaveLength

!  VBRDF and VSLEAVE output

   TYPE(VBRDF_Sup_Outputs)  , INTENT(INOUT) :: GEMSTOOL_BRDF_Results
   TYPE(VSLEAVE_Sup_Outputs), INTENT(INOUT) :: GEMSTOOL_SLEAVE_Results

!  Exception handling

   logical      , intent(out) :: fail
   integer      , intent(out) :: nmessages
   character*(*), intent(out) :: messages(3)

!  Local
!  =====

!  VBRDF supplement input structure

   TYPE(VBRDF_Sup_Inputs)                :: VBRDF_Sup_In

!  VBRDF supplement exception handling structure

   TYPE(VBRDF_Output_Exception_Handling) :: VBRDF_Sup_OutputStatus

!  VSLEAVE supplement input structure

   TYPE(VSLEAVE_Sup_Inputs)              :: VSLEAVE_Sup_In

!  other variables

   logical   :: do_restoration
   integer   :: m, k, g, ngeoms, nmoments_input
   real(fpk), parameter :: fzero = 0.0_fpk
   CHARACTER (LEN=10) :: CheckList ( 16 )
   CheckList = (/  'Lambertian', 'Ross-thin ', 'Ross-thick', 'Li-sparse ', &
                   'Li-dense  ', 'Hapke     ', 'Roujean   ', 'Rahman    ', &
                   'Cox-Munk  ', 'GissCoxMnk', 'GCMcomplex', 'BPDF-Vegn ', &
                   'BPDF-Soil ', 'BPDF-NDVI ', 'NewCMGlint', 'NewGCMGlit'/)

!  Initialize output exception handling

   fail       = .false.
   nmessages  = 0
   messages   = ' '

!  Proxy

   ngeoms = GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries

!  PART 1. WATERLEAVING CALCULATION
!  ================================

!  General inputs
!  --------------

!  Following flags are always required

   VSLEAVE_Sup_In%SL_DO_SLEAVING      = .TRUE.
   VSLEAVE_Sup_In%SL_DO_ISOTROPIC     = .FALSE.  ! Simplest - could be true.
   VSLEAVE_Sup_In%SL_DO_EXACT         = .TRUE.
   VSLEAVE_Sup_In%SL_DO_EXACTONLY     = .FALSE.
   VSLEAVE_Sup_In%SL_DO_FLUORESCENCE  = .FALSE.

   VSLEAVE_Sup_In%SL_DO_SOLAR_SOURCES = .TRUE.
   VSLEAVE_Sup_In%SL_DO_USER_OBSGEOMS = .TRUE.
   VSLEAVE_Sup_In%SL_DO_USER_STREAMS  = .TRUE.

!  Stokes and streams

   VSLEAVE_Sup_In%SL_NSTREAMS = GEMSTOOL_INPUTS%RTControl%NVlidort_nstreams
   VSLEAVE_Sup_In%SL_NSTOKES  = GEMSTOOL_INPUTS%RTControl%NVlidort_nstokes

!  Geometry - copy GEMSTOOL inputs
!  -------------------------------

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
!mick fix 1/24/2016 - initialize remaining input array elements
   VSLEAVE_Sup_In%SL_USER_OBSGEOMS(ngeoms+1:,:) = FZERO

!  Water-leaving Inputs from GEMSTOOL
!  ----------------------------------

!  Input Salinity in [ppt]
!  Input Chlorophyll concentration in [mg/M]
!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

   VSLEAVE_Sup_In%SL_SALINITY   = GEMSTOOL_INPUTS%Sleave%Salinity
   VSLEAVE_Sup_In%SL_CHLORCONC  = GEMSTOOL_INPUTS%Sleave%ChlorConc
   VSLEAVE_Sup_In%SL_WINDSPEED  = GEMSTOOL_INPUTS%Sleave%WindSpeed
   VSLEAVE_Sup_In%SL_WINDDIR    = FZERO

!    Flags for glint shadowing, Foam Correction, facet Isotropy

   VSLEAVE_Sup_In%SL_DO_GlintShadow   = GEMSTOOL_INPUTS%Sleave%DO_GlintShadow
   VSLEAVE_Sup_In%SL_DO_FoamOption    = GEMSTOOL_INPUTS%Sleave%DO_FoamOption
   VSLEAVE_Sup_In%SL_DO_FacetIsotropy = GEMSTOOL_INPUTS%Sleave%DO_FacetIsotropy
   if ( .not.GEMSTOOL_INPUTS%Sleave%DO_FacetIsotropy ) then
      VSLEAVE_Sup_In%SL_WINDDIR(1) = GEMSTOOL_INPUTS%Sleave%WindDir
   endif

!  Wavelength in [Microns]. Always needed
!  -----------------------

   VSLEAVE_Sup_In%SL_WAVELENGTH = Wavelength * 0.001d0  !convert from nm to um

!  SLEAVE calculation
!  ------------------

   CALL VSLEAVE_MAINMASTER ( 'DummyPath', &
        VSLEAVE_Sup_In,            & ! Inputs
        GEMSTOOL_SLEAVE_Results )    ! Outputs

!  If no Glitter calculation, return

   if ( .not. GEMSTOOL_INPUTS%Sleave%DO_OceanColor ) return

!  PART 2. GLITTER (GLINT) CALCULATION
!  ===================================

!  General inputs
!  --------------

!  Following flags are always required

   VBRDF_Sup_In%BS_DO_BRDF_SURFACE  = .TRUE.
   VBRDF_Sup_In%BS_DO_SOLAR_SOURCES = .TRUE.
   VBRDF_Sup_In%BS_DO_USER_OBSGEOMS = .TRUE.
   VBRDF_Sup_In%BS_DO_USER_STREAMS  = .TRUE.

   VBRDF_Sup_In%BS_DO_SURFACE_EMISSION = .FALSE.

!  Stokes and streams

   VBRDF_Sup_In%BS_NSTREAMS = GEMSTOOL_INPUTS%RTControl%NVlidort_nstreams
   VBRDF_Sup_In%BS_NSTOKES  = GEMSTOOL_INPUTS%RTControl%NVlidort_nstokes

!  Geometry - copy GEMSTOOL inputs
!  -------------------------------

   VBRDF_Sup_In%BS_NBEAMS          = ngeoms
   VBRDF_Sup_In%BS_N_USER_STREAMS  = ngeoms
   VBRDF_Sup_In%BS_N_USER_RELAZMS  = ngeoms
   VBRDF_Sup_In%BS_N_USER_OBSGEOMS = ngeoms
   do g =  1, ngeoms
      VBRDF_Sup_In%BS_BEAM_SZAS(g)         = GEMSTOOL_INPUTS%Geometry%GEMS_szas(g)
      VBRDF_Sup_In%BS_USER_ANGLES_INPUT(g) = GEMSTOOL_INPUTS%Geometry%GEMS_vzas(g)
      VBRDF_Sup_In%BS_USER_RELAZMS(g)      = GEMSTOOL_INPUTS%Geometry%GEMS_azms(g)
      VBRDF_Sup_In%BS_USER_OBSGEOMS(g,1)   = GEMSTOOL_INPUTS%Geometry%GEMS_szas(g)
      VBRDF_Sup_In%BS_USER_OBSGEOMS(g,2)   = GEMSTOOL_INPUTS%Geometry%GEMS_vzas(g)
      VBRDF_Sup_In%BS_USER_OBSGEOMS(g,3)   = GEMSTOOL_INPUTS%Geometry%GEMS_azms(g)
   ENDDO
!mick fix 1/24/2016 - initialize remaining input array elements
   VBRDF_Sup_In%BS_USER_OBSGEOMS(ngeoms+1:,:) = FZERO

!  SURFACE Inputs from GEMSTOOL
!  ----------------------------

!  Initialize for Water-Leaving, NewCM (Number 15 in the above check-list)
!   Parameters are not required for this kernel

   VBRDF_Sup_In%BS_N_BRDF_KERNELS         = 1
   VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = .false.
   VBRDF_Sup_In%BS_BRDF_PARAMETERS        = fzero

   VBRDF_Sup_In%BS_N_BRDF_PARAMETERS = 0   ;  VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1) = 3
   VBRDF_Sup_In%BS_WHICH_BRDF        = 0   ;  VBRDF_Sup_In%BS_WHICH_BRDF(1)        = 15
   VBRDF_Sup_In%BS_BRDF_NAMES        = ' ' ;  VBRDF_Sup_In%BS_BRDF_NAMES(1)        = 'NewCMGlint'
   VBRDF_Sup_In%BS_BRDF_FACTORS    = fzero ;  VBRDF_Sup_In%BS_BRDF_FACTORS(1) = 1.0d0

!  Other default values. Not required for this kernel

   VBRDF_Sup_In%BS_NSTREAMS_BRDF    = 100
   VBRDF_Sup_In%BS_DO_SHADOW_EFFECT = .false.
   VBRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY = .false.
   VBRDF_Sup_In%BS_DO_WSABSA_OUTPUT = .false.
   VBRDF_Sup_In%BS_DO_WSA_SCALING   = .false.
   VBRDF_Sup_In%BS_DO_BSA_SCALING   = .false.
   VBRDF_Sup_In%BS_WSA_VALUE = FZERO
   VBRDF_Sup_In%BS_BSA_VALUE = FZERO

!  Options with water leaving. These are the smae as the SLEAVE settings

   VBRDF_Sup_In%BS_DO_NewCMGLINT  = .true.
   VBRDF_Sup_In%BS_SALINITY       = GEMSTOOL_INPUTS%Sleave%Salinity
   VBRDF_Sup_In%BS_WINDSPEED      = GEMSTOOL_INPUTS%Sleave%WindSpeed
   VBRDF_Sup_In%BS_WINDDIR        = FZERO
   VBRDF_Sup_In%BS_DO_GlintShadow   = GEMSTOOL_INPUTS%Sleave%DO_GlintShadow
   VBRDF_Sup_In%BS_DO_FoamOption    = GEMSTOOL_INPUTS%Sleave%DO_FoamOption
   VBRDF_Sup_In%BS_DO_FacetIsotropy = GEMSTOOL_INPUTS%Sleave%DO_FacetIsotropy
   if ( .not. GEMSTOOL_INPUTS%Sleave%DO_FacetIsotropy ) then
      VBRDF_Sup_In%BS_WINDDIR(1)     = GEMSTOOL_INPUTS%Sleave%WindDir
   endif

!  Multiple-scattering Glitter options, Turn them off

   VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR        = .false.
   VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY = .false.
   VBRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER     = 0
   VBRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD   = 0
   VBRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD  = 0

!  Wavelength in [Microns]. Always needed
!  -----------------------

   VBRDF_Sup_In%BS_WAVELENGTH = Wavelength * 0.001d0

!  BRDF calculation
!  ----------------

!  Local dummy inputs

   do_Restoration  =  .false.
   NMOMENTS_INPUT  = 2 * VBRDF_Sup_In%BS_NSTREAMS

!  Call vector BRDF (VBRDF) supplement

    CALL VBRDF_MAINMASTER (  &
           DO_RESTORATION, NMOMENTS_INPUT, & ! Inputs
           VBRDF_Sup_In,            & ! Inputs
           GEMSTOOL_BRDF_Results,   & ! Outputs
           VBRDF_Sup_OutputStatus ) ! Output Status

!  Exception handling. Up to 3 messages allowed

    if ( VBRDF_Sup_OutputStatus%BS_STATUS_OUTPUT .ne. 0 ) then
       nmessages = min(VBRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES,3)
       do m = 1, nmessages
          messages(m) = Trim(VBRDF_Sup_OutputStatus%BS_OUTPUTMESSAGES(m))
       enddo
       fail = .true. ; return
    endif

!  Quick debug

!  write(*,44)'Waterleaving @ wavelength ',Wavelength,&
!             GEMSTOOL_SLEAVE_Results%SL_SLTERM_ISOTROPIC(1,1),&
!             GEMSTOOL_SLEAVE_Results%SL_SLTERM_F_0(0,1,1:8,1)
!  write(*,44)'NewCM BRDF   @ wavelength ',Wavelength,&
!             GEMSTOOL_BRDF_Results%BS_DBOUNCE_BRDFUNC(1,1,1,1),&
!             GEMSTOOL_BRDF_Results%BS_BRDF_F_0(0,1,1:8,1)
!44 format(A,f10.4,1p9e12.4)

!pause

!  Debug write-up

!  PLACEHOLDER

!  Done

   return
end subroutine GEMSTOOL_WaterLeaving

!  End module

end module GEMSTOOL_UVN_WaterLeaving_m
