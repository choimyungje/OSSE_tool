module GEMSTOOL_BRDFSurface_m

!  Rob, created 10/24/16. Very simple version, first attempt.

! This Module is designed to calculate the BRDF reflectances required by VLIDORT 
!   -- Basically, we fill us the inputs for the VBRDF supplement, then call it
!   -- Results are output to a local type structure which is modeled after
!   -- the VLIDORT VBRDF structure

! Routine calls the VLIDORT VBRDF supplement, operating with Kernel choices
!  -- inputs are provided by variables in the GEMSTOOL Inputs

use GEMSTOOL_Input_Types_m

!  Rob, add the BRDF type structures 

USE VBRDF_SUP_MOD

public
private fpk

contains

subroutine GEMSTOOL_BRDFSurface &
      ( GEMSTOOL_INPUTS,  &
        GEMSTOOL_BRDF_Results, fail, nmessages, messages )
     
   implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

! GEOMTOOL inputs, structure. Intent(in) here

   TYPE(GEMSTOOL_Config_Inputs), INTENT(IN) :: GEMSTOOL_INPUTS

!  VBRDF output

   TYPE(VBRDF_Sup_Outputs), INTENT(INOUT) :: GEMSTOOL_BRDF_Results

!  Exception handling

   logical      , intent(out) :: fail
   integer      , intent(out) :: nmessages
   character*(*), intent(out) :: messages(3)

!  Local
!  =====

!  VBRDF supplement input structure

   TYPE(VBRDF_Sup_Inputs)               :: VBRDF_Sup_In

!  VBRDF supplement exception handling structure

   TYPE(VBRDF_Output_Exception_Handling) :: VBRDF_Sup_OutputStatus

!  other variables

   logical   :: do_restoration
   integer   :: m, k, g, ngeoms, nmoments_input
   real(fpk), parameter :: fzero = 0.0_fpk
   CHARACTER (LEN=10) :: CheckList ( 16 )
   CheckList = (/  'Lambertian', 'Ross-thin ', 'Ross-thick', 'Li-sparse ', &
                   'Li-dense  ', 'Hapke     ', 'Roujean   ', 'Rahman    ', &
                   'Cox-Munk  ', 'GissCoxMnk', 'GCMcomplex', 'BPDF-Vegn ', &
                   'BPDF-Soil ', 'BPDF-NDVI ', 'NewCMGlint',  'NewGCMGlit'/)

!  Initialize output exception handling

   fail       = .false.
   nmessages  = 0
   messages   = ' '

!  proxy

   ngeoms = GEMSTOOL_INPUTS%Geometry%N_GEMS_geometries

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

!  SURFACE Inputs
!  ==============

!  From GEMSTOOL
!  -------------

!  initialize

   VBRDF_Sup_In%BS_N_BRDF_PARAMETERS = 0
   VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = .false.
   VBRDF_Sup_In%BS_WHICH_BRDF = 0  
   VBRDF_Sup_In%BS_BRDF_NAMES = ' '
   VBRDF_Sup_In%BS_BRDF_FACTORS   = fzero
   VBRDF_Sup_In%BS_BRDF_PARAMETERS = fzero

!  Copy from GEMSTOOL inputs

   VBRDF_Sup_In%BS_N_BRDF_KERNELS    = GEMSTOOL_INPUTS%BRDF%n_kernels
   do g = 1, VBRDF_Sup_In%BS_N_BRDF_KERNELS
      VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(g) = GEMSTOOL_INPUTS%BRDF%Kernel_n_parameters(g)
      VBRDF_Sup_In%BS_WHICH_BRDF(g)        = GEMSTOOL_INPUTS%BRDF%Kernel_indices(g)
      VBRDF_Sup_In%BS_BRDF_NAMES(g)        = CheckList(GEMSTOOL_INPUTS%BRDF%Kernel_indices(g))
      VBRDF_Sup_In%BS_BRDF_FACTORS(g)      = GEMSTOOL_INPUTS%BRDF%Kernel_factors(g)
      if ( GEMSTOOL_INPUTS%BRDF%Kernel_indices(g) .eq. 1 ) VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG(g) = .true.
      do k = 1, VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(g)
        VBRDF_Sup_In%BS_BRDF_PARAMETERS(g,k)  = GEMSTOOL_INPUTS%BRDF%Kernel_parameters(g,k)
      enddo
   enddo

!  Other default values
!  --------------------

   VBRDF_Sup_In%BS_NSTREAMS_BRDF    = 100
   VBRDF_Sup_In%BS_DO_SHADOW_EFFECT = .true.

   VBRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY = .false.

   VBRDF_Sup_In%BS_DO_WSABSA_OUTPUT = .false.
   VBRDF_Sup_In%BS_DO_WSA_SCALING   = .false.
   VBRDF_Sup_In%BS_DO_BSA_SCALING   = .false.
   VBRDF_Sup_In%BS_WSA_VALUE = FZERO
   VBRDF_Sup_In%BS_BSA_VALUE = FZERO

!  options with water leaving (Not currently implemented)

   VBRDF_Sup_In%BS_DO_NewCMGLINT  = .false.
   VBRDF_Sup_In%BS_SALINITY       = FZERO
   VBRDF_Sup_In%BS_WAVELENGTH     = FZERO
   VBRDF_Sup_In%BS_WINDSPEED      = FZERO
   VBRDF_Sup_In%BS_WINDDIR        = FZERO
   VBRDF_Sup_In%BS_DO_GlintShadow   = .false.
   VBRDF_Sup_In%BS_DO_FoamOption    = .false.
   VBRDF_Sup_In%BS_DO_FacetIsotropy = .false.

!   VBRDF_Sup_In%BS_CHLORCONC      = FZERO
!   VBRDF_Sup_In%BS_DO_NewGCMGLINT = .false.

!  Multiple-scattering Glitter options

   VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR        = .false.
   VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY = .false.
   VBRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER     = 0
   VBRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD   = 0
   VBRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD  = 0

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

!  Debug
!   open(55, file = 'BRDF_model_input.out')
!   write(55,*)'Directbounce',GEMSTOOL_BRDF_Results%BS_DBOUNCE_BRDFUNC(1,1,1,1)
!   do m = 0, NMOMENTS_INPUT-1
!      write(55,44)m,GEMSTOOL_BRDF_Results%BS_BRDF_F(m,1,1,1:VBRDF_Sup_In%BS_NSTREAMS)
!   enddo
!44 format(i4,1p8e12.4)
!   close(55) 

!  Done

   return
end subroutine GEMSTOOL_BRDFSurface

!  End module

end module GEMSTOOL_BRDFSurface_m
