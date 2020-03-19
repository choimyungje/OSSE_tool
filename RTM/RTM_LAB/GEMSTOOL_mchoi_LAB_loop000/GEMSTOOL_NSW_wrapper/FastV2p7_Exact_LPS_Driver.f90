PROGRAM FastV2p7_Exact_LPS_Driver

!  Original TEMPO notes ---------
!     First and Second Versions, August 2014. Not Successful
!     Third Version, Upgrade to VLIDORT 2.7, 30 June 2015
!     Linearization, Upgrade to VLIDORT 2.7, 20 July 2015

!  Notes for the FASTPCA TOOL--------
!     * This Driver adapted  14 December  2015  (TEMPO trials)
!     * This Driver upgraded 29 December  2015  (Geocape Leveraging)
!     * This Driver upgraded 26 January   2016  (Full testing of system)
!     * This Driver upgraded 18 February  2016  (IPA and clouds, Bookkeeping)
!     * This Driver upgraded 30 March     2016  (Convolution/Raman/IPA)
!     * This Driver upgraded 4-6 May      2016  (RingSF/Albcoeffs/ILS/Outputs)
!     * This Driver upgraded 22 September 2016  (Optional use of PGE-TES atmospheric system)

!  ORIGINAL HISTORY NOTES have been omitted (see text file for details)

!  VLIDORT pars file (for checking dimensions)

   use vlidort_pars, only : MAXSTREAMS, MAXLAYERS, MAXMOMENTS_INPUT

!  Case Selector and Bookkeeping Modules

   use FastV2p7_Case_Selector_m
   use FastV2p7_MiscBkkeep_m , only : mkchar, WriteTiming_Exact

!  Create properties module

   use FastV2p7_CreateProps_Plus_Master_m

!  Convolution modules
!       AsymGauss tool first introduced 3/29/16
!       AsymGauss tool updated          5/5/16
!       OMIILS    tool first introduced 5/5/16

   use AsymGauss_Convolution_Tool_m
   use OMIILS_Convolution_Tool_m

!  Sioris Raman Modules. Introduced 3/30/16

   use Sioris_Raman_Specdata_m
   use Sioris_Raman_m

!  LRRS Raman Module. Introduced 4/1/16

   use FastV2p7_LRRS2p3_Master_m

!  Master Module

   use FastV2p7_Exact_LPS_Master_m

   implicit none

!  Constants and EXTERNAL DIMENSIONING
!  ===================================

!  precision

   integer, parameter :: dpm = SELECTED_REAL_KIND(15)
   integer, parameter :: spm = SELECTED_REAL_KIND(6)

!  Number of points, layers, moments
!    Aerosol moments, No bigger than 2*MAXSTREAMS = MAXMOMENTS + 1

!  Original Geocape settings
!  -------------------------

!   integer, parameter :: E_ndat      = 45000
!   integer, parameter :: E_nlayers   = 114      ! GEOCAPE
!   integer, parameter :: E_ngases    = 5        ! 5 gases
!   integer, parameter :: E_nmoms_all = 32       ! Aerosols

!  Tempo settings (for testing)
!  ----------------------------

!   integer, parameter :: E_ndat      = 11001    ! TEMPO, 540-650 nm @ 0.01 nm
!   integer, parameter :: E_nlayers   = 47       ! TEMPO May/Aug 2014
!   integer, parameter :: E_nmoms_all = 32       ! Aerosols
!   integer, parameter :: E_ngases    = 5        ! 5 gases (O3,NO2,H2O,O4,O2)

!  Tool settings (2/19/16 update)
!  ------------------------------

!  Fine-grid and general settings

   integer, parameter :: E_ndat      = 5201     ! MW1 maximum size
   integer, parameter :: E_nlayers   = 64       ! AIRS-OMI value.        Same value in LRRS_pars and VLIDORT_pars
   integer, parameter :: E_nlevels   = 65       ! AIRS-OMI value (always 1 more, used for Level Jacobians)
   integer, parameter :: E_nmoms_all = 32       ! Aerosols, Leveraged.   Same value in LRRS_pars and VLIDORT_pars
   integer, parameter :: E_ngases    = 5        ! 5 gases (O3,NO2,H2O,O4,O2)
   integer, parameter :: E_nscenes   = 3        ! Clear/Cloud/Total

!  Coarse settings. Added 3/28/16. Same value in LRRS_pars (4/2/16)

   integer, parameter :: E_ndat_c      = 200     ! 559 points in OMI Channel 2 overall, ~186 MW2 points

!mick test 2/26/2016
   !INTEGER, parameter :: E_npars     = 1 !for testing with 1 gas
   INTEGER, parameter :: E_npars     = 6 !for testing with 1 gas + 5 aero
   INTEGER, parameter :: E_nspars    = 1

!  OMI ILS settings, Added 5/4/16

   integer, parameter :: E_ndat_ILS   = 505

!  Geometry settings
!  -----------------

! Observational geometry only
!    Same values in LRRS_pars and VLIDORT_pars

   integer, parameter :: E_ngeoms    = 1
   integer, parameter :: E_nszas     = 1
   integer, parameter :: E_nvzas     = 1
   integer, parameter :: E_nazms     = 1

!  Numbers

   real(kind=dpm), parameter :: pi = 3.14159265358979323846_dpm
   real(kind=dpm), parameter :: dtor = pi/180.0_dpm
   real(kind=dpm), parameter :: epsilon = 1.0e-08_dpm

!  Input Control
!  =============

!  High-level methods for Convolution and Raman. Added 3/30/16, OMIILS added 5/4/16

   logical       :: OMIILS_Convolution_Method
   logical       :: Original_Convolution_Method
   logical       :: Sioris_Raman_Method

!  AIRS_OMI Data Path

   character*256 :: AIRSOMI_DataPath

!  Geocape Path (feature introduced 5/15/15)

   character*256 :: Geocape_Path

!  Aerosol OD Reference Wavelength Control

   integer :: Select_index

!  AIRS_OMI Scene Index

   integer :: AIRSOMI_SceneIndex

!  OMI Pixel index introduced 5/4/16. Initially set to the AIRS_OMI scene index

   integer :: OMI_PixelIndex

!  Input choices for 'Configfiles/FastV2p7_Main.cfg'
!  ------------------------------------------------

!  Optional TES usage flag. 9/22/16. RT Solutions

   logical :: do_PGE_TES

!  Boolean flag (aerosol inclusion)

   logical :: do_aerosols

!  No gas flag

   logical :: do_Continuum

!  OMI Window (MW1, MW2, MW3) index. 3/29/16
!     --OMI Channel (1 or 2). Value = 1 for MW1, = 2 for MW2/MW3

   integer :: OMI_Channel
   integer :: OMI_window

!  Base window limits. New 3/31/16

   real(kind=dpm) :: wstart_base, wfinis_base

!  Spectral window parameters

   integer        :: ndat
   real(kind=dpm) :: wstart, wresol

!  Gas control,

   integer          :: ngases
!mick fix 2/26/2016 - modified length of which_gases
   !character(Len=4) :: which_gases (E_ngases)
   character(Len=6) :: which_gases (E_ngases)

!  Input choices from 'Configfiles/FastV2p7_Lineariz.cfg'. OR HARDWIRED
!  --------------------------------------------------------------------

!  Linearization control

   Logical   :: which_gaspars (E_ngases), do_GasAbs_WF, do_AOD_WF(5)
   LOGICAL   :: DO_SURFACE_WFS
   LOGICAL   :: DO_PROFILE_WFS
   LOGICAL   :: LFVARY(E_nlayers)
   INTEGER   :: LNVARY(E_nlayers)
   INTEGER   :: NPARS, NSPARS

!  Input choices for 'Configfiles/FastV2p7_RT.cfg'
!  -----------------------------------------------

!  Enhanced sphericity flag

   logical :: do_enhanced_ps

!  Observational-geometry flag

   logical :: do_obsgeoms

!  Number of VLIDORT discrete ordinates and Stokes components

   integer :: nstreams, nstr2, nmoms_2ns, nstokes

!  Optional control for including the Exact 2S calculation
!    REALLY FOR DEBUGGING and ANALYSIS

   logical :: do_optional_2stream

!  2STREAM BVP INdex (0 = LAPACK,1 = PentaDiag) and inverse flag

   logical :: DO_PDINVERSE
   INTEGER :: BVPINDEX

!  Earth radius

   real(kind=dpm)  :: eradius

!  Input choices for 'Configfiles/FastV2p7_ObsGeoms.cfg' OR 'Configfiles/FastV2p7_LatticeGeoms.cfg'
!  ------------------------------------------------------------------------------------------------

!    ngeoms              = Number of geometries
!    nszas, nvzas, nazms = Number of angles

   integer :: ngeoms
   integer :: nszas
   integer :: nvzas
   integer :: nazms

!  Geometry. Scattering angle added 3/30/16 (for Raman calculation)

   real(kind=dpm)  :: sza_boa(E_nszas), vza_boa(E_nvzas), azm_boa(E_nazms)
   real(kind=dpm)  :: scatang_boa(E_ngeoms)

!  Optical Variables, I/O from Create_properties Master
!  ====================================================

!  Control (output). Ncloud_layers added, 2/17/16

   integer        :: nlayers, nlevels, ncloud_layers, nmoms_all

!  input wavelengths

   real(kind=dpm) :: lambdas(E_ndat)       ! input
   real(kind=dpm) :: wavnums(E_ndat)       ! input, GEOCAPE Only

!  whether to use hitran or not. NOT REQUIRED HERE

   logical :: use_hitran

!  Output from the create_properties routine. All Fine-grid

   real(kind=dpm) :: albedo(E_ndat)
   real(kind=dpm) :: depol (E_ndat)
   real(kind=dpm) :: RayXsec (E_ndat)
   real(kind=dpm) :: taug  (E_nlayers,E_ndat)
   real(kind=dpm) :: taudp (E_nlayers,E_ndat)
   real(kind=dpm) :: omega (E_nlayers,E_ndat)
   real(kind=dpm) :: heights(0:E_nlayers)
   real(kind=dpm) :: fr(E_nlayers,E_ndat)
   real(kind=dpm) :: fa(5,E_nlayers,E_ndat)

   real(kind=dpm) :: solar_flux (E_ndat)                ! GeoCape output
   real(kind=dpm) :: aergreeks(2,0:E_nmoms_all,5,6)     ! GeoCape output

!  Local finegrid array

   real(kind=dpm) :: omegamoms_f (E_nlayers,E_ndat)

!  You might need these again !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   real(kind=dpm) :: taur  (E_nlayers,E_ndat)
!   real(kind=dpm) :: taua  (E_nlayers,E_ndat)
!   real(kind=dpm) :: fa(E_nlayers,E_ndat)                          ! Tempo output
!   real(kind=dpm) :: gascols(E_ngases)                             ! Tempo output
!   logical        :: aerflags(E_nlayers,E_ndat)                    ! Tempo output
!   real(kind=dpm) :: aercoeffs(0:E_nmoms_all,E_nlayers,6,E_ndat)   ! Tempo output

   real(kind=dpm) :: tausp (E_nlayers,E_ndat,E_ngases)         ! Individual Gases
   real(kind=dpm) :: tauO3 (E_nlayers,E_ndat)                 ! Just 1 gas = O3....

!  Linearized.  First dimension is always Linearization parameters.

   real(kind=dpm) :: L_taudp(E_npars,E_nlayers,E_ndat)
   real(kind=dpm) :: L_omega(E_npars,E_nlayers,E_ndat)
   real(kind=dpm) :: L_fr   (E_npars,E_nlayers,E_ndat)
   real(kind=dpm) :: L_fa   (E_npars,5,E_nlayers,E_ndat)

!  Linearization of Gas absorption w.r.t Log(VMR)
!     Additional argument, 3/28/16

   real(kind=dpm)  :: dtaug_dLogv(E_nlayers,2,E_ndat,E_ngases)

!  Other outputs
!  -------------

!  Layer Temperatures and Air Columns (needed for Sioris-Raman calculations)
!    ---Added output, 3/30/16.  Aircolumns [mol/cm/cm], temperatures in [K]

   real(kind=dpm)    :: aircolumns  (E_nlayers)
   real(kind=dpm)    :: temperatures(E_nlayers)

!  Angles (1=SZA,2=EAZ,3=VZA,4=SCA)

   real(kind=dpm) :: Angles(2,4)

!  CLoud information

   real(kind=dpm) :: CldFraction, CldTopPressure, CldTopAlbedo

!  Geographical target attribute

   real(kind=dpm) :: EarthRadius

!  OMI specifications, Ring scaling factor, albedo coefficients
!    ---Added output, 5/4/16. 

   real(kind=dpm) :: Ring_Scaling
   real(kind=dpm) :: Albedo_0, Albedo_1, Lambda_Ref_MW2

!  Intensity/Profile-Jacobians and Status, output results
!  ======================================================

!  Variables for Exact intensity
!    Scene 1 = Clearsky, 2 = cloudy sky, 3 = IPA total

   real(kind=dpm) :: intensity_Exact   (E_ndat,E_ngeoms,E_nscenes)
   real(kind=dpm) :: intensity_LD_exact(E_ndat,E_ngeoms,E_nscenes)
   real(kind=dpm) :: intensity_2S_Exact(E_ndat,E_ngeoms,E_nscenes)
   real(kind=dpm) :: intensity_FO_Exact(E_ndat,E_ngeoms,E_nscenes)

!  Variables for Exact Jacobians
!    Scene 1 = Clearsky, 2 = cloudy sky, 3 = IPA total

   REAL(kind=dpm) :: LP_Jacobians_exact    (E_npars,E_nlayers,E_ndat,E_ngeoms,E_nscenes)
   REAL(kind=dpm) :: LP_Jacobians_LD_exact (E_npars,E_nlayers,E_ndat,E_ngeoms,E_nscenes)
   REAL(kind=dpm) :: LP_Jacobians_2S_exact (E_npars,E_nlayers,E_ndat,E_ngeoms,E_nscenes)
   REAL(kind=dpm) :: LP_Jacobians_FO_exact (E_npars,E_nlayers,E_ndat,E_ngeoms,E_nscenes)

   REAL(kind=dpm) :: LS_Jacobians_exact    (E_nspars,E_ndat,E_ngeoms,E_nscenes)
   REAL(kind=dpm) :: LS_Jacobians_LD_exact (E_nspars,E_ndat,E_ngeoms,E_nscenes)
   REAL(kind=dpm) :: LS_Jacobians_2S_exact (E_nspars,E_ndat,E_ngeoms,E_nscenes)
   REAL(kind=dpm) :: LS_Jacobians_FO_exact (E_nspars,E_ndat,E_ngeoms,E_nscenes)

!  Variables for Exact Trace-gas Log(VMR) LEVEL Jacobians
!    Introduced 3/29/16. Only work with the final result.

   REAL(kind=dpm) :: LP_LogV_Jacobians_exact  (E_ngases,E_nlevels,E_ndat,E_ngeoms,E_nscenes)

!  Variables for Exact Scattering Weights
!   12/30/15 Status uncertain
!    Scene 1 = Clearsky, 2 = cloudy sky, 3 = IPA total

   real(kind=dpm) :: Scatwts_Exact      (E_nlayers,E_ndat,E_ngeoms,E_nscenes)
   real(kind=dpm) :: Scatwts_LD_Exact   (E_nlayers,E_ndat,E_ngeoms,E_nscenes)
   real(kind=dpm) :: Scatwts_2S_Exact   (E_nlayers,E_ndat,E_ngeoms,E_nscenes)
   real(kind=dpm) :: Scatwts_FO_Exact   (E_nlayers,E_ndat,E_ngeoms,E_nscenes)

!  Coarse-grid variables. Introduced, 3/29/16. Linearized variables 3/30/16

   integer        :: ndat_c
   real(kind=dpm) :: lambdas_c   (E_ndat_c)
   real(kind=dpm) :: solarflux_c (E_ndat_c)
   real(kind=dpm) :: IExact_c    (E_ndat_c,E_ngeoms,E_nscenes)

   real(kind=dpm) :: asym_c      (E_ndat_c)     ! For Convolution
   real(kind=dpm) :: hwhm_c      (E_ndat_c)     ! For Convolution

   REAL(kind=dpm) :: LogV_JExact_c  (E_ngases,E_nlevels,E_ndat_c,E_ngeoms,E_nscenes)
   REAL(kind=dpm) :: LS_JExact_c    (E_nspars,E_ndat_c,E_ngeoms,E_nscenes)

   real(kind=dpm) :: taudp_c     (E_ndat_c,E_nlayers)                  !  For Sioris-Raman
   real(kind=dpm) :: deltau_c    (E_nlayers,E_ndat_c)                  !  For LRRS
   real(kind=dpm) :: omegamoms_c (E_nlayers, 0:E_nmoms_all, E_ndat_c)  !  For LRRS
   real(kind=dpm) :: Depol_c     (E_ndat_c), RayXsec_c   (E_ndat_c)    !  For LRRS

!  Auxiliary albedo coefficient Jacobians. Added, 5/5/16

   real(kind=dpm) :: JExact_AlbCoeffs_C (2,E_ndat_c,E_ngeoms,E_nscenes)

!  ILS convolution variables and arrays. Introduced 5/4/16

   integer        :: ndat_ILS
   integer        :: ILS_offsets (E_ndat_c)
   integer        :: ILS_npoints (E_ndat_c)
   real(kind=dpm) :: ILS_wavs    (E_ndat_ILS,E_ndat_c)
   real(kind=dpm) :: ILS_values  (E_ndat_ILS,E_ndat_c)
   logical        :: Fail_ILS
   character*100  :: Message_ILS

!  OMI-specific arrays on the coarse grid. Linearized variables 3/30/16, Updated 5/6/16
!  --------------------------------------

!  Filling and Flux

   real(kind=dpm) :: OMI_solarflux_c (E_ndat_c)
   real(kind=dpm) :: OMI_filling_c   (E_ndat_c,E_ngeoms,E_nscenes)

!  Main Output. Formalized, 5/6/16. Added LFJ (cloud fraction Jacobian), 9/21/16.

   real(kind=dpm) :: OMI_RadExact_c  (E_ndat_c,E_ngeoms)
   REAL(kind=dpm) :: OMI_LPJExact_c  (E_ngases,E_nlevels,E_ndat_c,E_ngeoms)
   REAL(kind=dpm) :: OMI_LAJExact_c  (2,E_ndat_c,E_ngeoms)
   REAL(kind=dpm) :: OMI_LRJExact_c  (E_ndat_c,E_ngeoms)
   REAL(kind=dpm) :: OMI_LFJExact_c  (E_ndat_c,E_ngeoms)

!  Raman variables
!  ---------------

!  Sioris-Raman variables. NOTE the Sioris-Raman code is NOT LINEARIZED

   integer        :: SR_nlayers
   real(kind=dpm) :: SR_albedo
   logical        :: SR_do_upwelling, Fail_SRaman
   character*100  :: Message_SRaman
   TYPE(Sioris_Raman_Specdata)      :: RamanData

!  LRRS Variables. Not using the Linearized LRRS code (YET!!!). 2/4/16

   integer        :: LRRS_nlayers, Nmoments_LRRS
   real(kind=dpm) :: LRRS_Albedo   (E_ndat_c)

   real(kind=dpm) :: Elastic_c    (E_ndat_c)
   real(kind=dpm) :: Raman_c      (E_ndat_c)
   real(kind=dpm) :: Filling_SS_c (E_ndat_c)

   logical        :: Fail_LRRS_Raman
   integer        :: n_LRRS_messages
   character*100  :: LRRS_Messages(100)

!  Exception handling variables
!  ----------------------------

   integer, parameter :: b_max_messages   = 25
   logical       :: b_fail
   integer       :: b_nmessages
   character*100 :: b_messages (b_max_messages)

!  LOCAL variables
!  ===============

!  Local Raman-code settings. 2/18/16

   integer            :: local_nlayers
   real(kind=dpm)     :: local_albedo(E_ndat)
   logical            :: local_do_surf_wfs

!  Help variables, some new 2/18/16, 3/29/16, 3/30/16

   logical            :: do_Dump1, do_Create1
   integer            :: i,j,k,l,n,n1,m,v,w, idum, ndum, FDIndex, q, qs, g
   integer            :: ndatr, nstep, nwav_ch1, nwav_ch2
   character(len=1)   :: c1, c1d
   character(len=2)   :: c2d
   character(len=3)   :: c3d
   character(len=4)   :: c4, c4d
   character(len=5)   :: Transept
   character(len=6)   :: Cstokes
   real(kind=dpm)     :: taugas, tauray, tautot, tauaer, outgastau(E_ndat)
   real(kind=dpm)     :: beta2, kn1, kn1_gasabs, kn2, kn2_gasabs, Filling, albfac
   real(kind=dpm)     :: cx0, cx1, sx0, sx1, lam, sol

!  Component arrays. 5/6/16 updated

   real(kind=dpm)     :: Rad_Component (E_nscenes)
   real(kind=dpm)     :: LPJ_COmponent (E_ngases,E_nlevels,E_nscenes)
   real(kind=dpm)     :: LAJ_Component (2,E_nscenes)
   real(kind=dpm)     :: LRJ_Component (E_nscenes)

!  output filenames. Revised 4/1/16. Upgraded 5/6/16

   character(len=120) :: Fileheader0, Fileheader, Dumpfile, Results_Path, Fileheader_save(E_ngeoms)
   character(len=150) :: OMIOutput_Timing(E_ngeoms), OMIOutput_AllResults(E_ngeoms),  OMIOutput_LPJContour(E_ngeoms)
   character(len=150) :: fileOutput1_FineGrid(E_ngeoms), fileOutput2_FineGrid(E_ngeoms)
   character(len=150) :: fileOutput3_FineGrid(E_ngeoms), fileOutput4_FineGrid(E_ngeoms)

!  VLIDORT initialize flag

   logical        :: DO_VLidort_initialize

!  Scene bookkeeping

   integer        :: nscenes, iscene, iout, output_index(2)
   logical        :: do_clearsky(2)
   character*5    :: Charscene(2)
   real(kind=dpm) :: ClrFraction
   real(kind=dpm), parameter :: CldFraction_Threshold_Lower = 0.05_dpm
   real(kind=dpm), parameter :: CldFraction_Threshold_Upper = 0.95_dpm

!  Timing variables
!     ---- Raman/Convolution/IPA added 3/30/16

   logical        :: Monitor_CPU
   real(kind=spm) :: e1,e2,e3
   real(kind=spm) :: CreateProptime, Writetime, Overalltime
   real(kind=spm) :: Exacttimes(9,2), Exactruntime(2)
   real(kind=spm) :: Convolutiontime, RamanFilltime, IPACalctime, PostProcTime

!  Local debugging flags. Revised 5/6/16
!  ---------------------

!  General Debug flag

   logical, parameter :: do_debug_ExactTool =.true.

!  short test flag (for debugging). 4/1/16.
!    If set, Avoids Convolution and Raman

   logical, parameter :: do_shorttest   = .false.
   !logical, parameter :: do_shorttest  = .true.

!  Dump/Create Flags

   !logical, parameter :: do_Dump   = .false.
   logical, parameter :: do_Dump   = .true.
   !logical, parameter :: do_Create = .false.
   logical, parameter :: do_Create = .true.

!  Skip fine calculation
!    If set, Avboids the Finegrid RT calculation

   logical, parameter :: do_Skip_Fine = .false.
!   logical, parameter :: do_Skip_Fine = .true.

!  output Contour plot of profile Jacobians
!    If set, will output profile Jacobians in contour-plot form (GNUPLOT)

   logical, parameter :: do_LPJ_ContourPlot = .true.

!  output fine Grid results
!    If set, will output 4 files with fine-grid Radiances + 3 Jacobians

   logical, parameter :: do_Output_FineGrid = .true.

!  0. INITIAL SECTION
!  ==================

!  Set high-level Convolution and Raman methods. Added 3/30/16, 5/4/16

!   OMIILS_Convolution_Method = .true.
   OMIILS_Convolution_Method = .false.

!   Original_Convolution_Method = .false.
   Original_Convolution_Method = .true.

!   Sioris_Raman_Method         = .true.
   Sioris_Raman_Method         = .false.

!  Results path

   Results_Path = 'Results/'

!  Data paths

   Geocape_Path     = '../../GEOCAPE_DataBase/'
   AIRSOMI_DataPath = '../../AIRSOMI_DataBase/'

!  Basic selection
!    3/29/16. OMI_Channel, now dependent on window selection

!   OMI_Channel      = 1
   AIRSOMI_SceneIndex = 106
   Transept = 'Tr000' ; call mkchar(3,c1d,c2d,Transept(3:5),c4d,AIRSOMI_SceneIndex)

!  OMI pixel index. 5/4/16, initially set to AIRSOMI_SceneIndex
!    Really needs to come from someplace else................................*!*!*!*!*!*!*!*!

   OMI_PixelIndex = AIRSOMI_SceneIndex

!  FD Index 

   FDIndex = 0

!  initialize exception handling

   b_fail      = .false.
   b_nmessages = 0
   b_messages  = ' '

!  VLIDORT Needs to be Initialised

   DO_VLidort_initialize = .false.

!  Initializing the CPU times. Initialization has been re-introduced here 2/18/16
!     ---- Raman/Convolution/IPA added 3/30/16

   CreatePropTime  = 0.0
   Exacttimes      = 0.0
   ExactRuntime    = 0.0

   Convolutiontime = 0.0
   RamanFilltime   = 0.0
   IPACalctime     = 0.0

   Writetime       = 0.0
   Overalltime     = 0.0

!  Monitor CPU flag

   Monitor_CPU = .true.

!  initialize some inputs

   nlayers = 0 ; nmoms_all = 0

!  Open and Read configuration files
!  ---------------------------------

!     File 1  = Main     control
!     File 2  = RT Model control
!     File 3A = Geometry control (Observational)
!     File 3B = Geometry control (Lattice)

!  File 1
!  ------

!  Tempo example

!   if ( do_TempoCase ) then
!      OPEN(1,file='Configfiles/FastV2p7_Main.cfg',status='old')
!      read(1,*)do_Continuum
!      read(1,*)do_aerosols
!      read(1,*)nlayers
!      read(1,*)nmoms_all
!      read(1,*)ndat, wstart, wresol
!      if ( do_Continuum ) then
!         ngases = 0 ; which_gases = '    '
!      else
!         read(1,*)ngases
!         read(1,*)(which_gases(g),g=1,ngases)
!      endif 
!      Close(1)
!   endif

!  Top Level input file
!  9/22/16. PGE TES option added

   OPEN(1,file='Configfiles/FastV2p7_Main.cfg',status='old')
   read(1,*)do_PGE_TES
   read(1,*)do_Continuum
   read(1,*)do_aerosols
   read(1,*)OMI_Window    ! MW1, MW2 or MW3
   Close(1)

!  Set OMI Channel. 3/29/16

   if ( OMI_Window .eq. 1 ) then
      OMI_Channel = 1
   else
      OMI_Channel = 2
   endif
 
!  Case selection.
!    Revised 3/30/16 to output base-window limits. Short-test 4/1/16.
   
   write(*,*)
   write(*,'(1x,a,i2.2)') 'Doing OMI_Window = ',OMI_Window
   CALL FastV2p7_Case_Selector &
      ( E_ngases, OMI_Window, do_shorttest,             & ! Input
        wstart_base, wfinis_base, ndat, wstart, wresol, & ! Window output
        ngases, which_gases, Select_index, use_hitran )   ! Other output

!  NOT REQUIRED
  ! if ( do_Continuum ) then
  !    ngases = 0 ; which_gases = '      '
  !    ndat = 0 ! Need to change this 
  !    wstart = 0.d0 ! Need to change this
  !    wresol = 0.d0 ! Need to change this
  ! endif

!  File 2
!  ------

   OPEN(1,file='Configfiles/FastV2p7_RT.cfg',status='old')
   read(1,*)do_enhanced_ps
   read(1,*)do_obsgeoms
   read(1,*)nstreams
   read(1,*)nstokes
   read(1,*)do_optional_2stream
   read(1,*)eradius
   close(1)

!  Other quantities

   nstr2 = 2 * nstreams
   nmoms_2ns = nstr2

!  other 2-stream options

   DO_PDINVERSE = .false.
   BVPINDEX     = 1

!  File 3A
!  -------

   if ( do_obsgeoms ) then
      OPEN(1,file='Configfiles/FastV2p7_ObsGeoms.cfg',status='old')
      read(1,*)ngeoms
      if ( ngeoms .gt. E_ngeoms ) go to 53         ! Skip if out of bounds
      do n = 1, ngeoms
         read(1,*)sza_boa(n), vza_boa(n), azm_boa(n)
      enddo
      nszas = ngeoms ; nvzas  = ngeoms ; nazms = ngeoms !  Necessary otherwise crashes
      close(1)
   endif

!  File 3B
!  -------

   if ( .not.do_obsgeoms ) then
      OPEN(1,file='Configfiles/FastV2p7_LatticeGeoms.cfg',status='old')
      read(1,*)nszas
      if ( nszas .gt. E_nszas ) go to 53         ! Skip if out of bounds
      do n = 1, nszas
         read(1,*)sza_boa(n)
      enddo
      read(1,*)nvzas
      if ( nvzas .gt. E_nvzas ) go to 53         ! Skip if out of bounds
      do n = 1, nvzas
         read(1,*)vza_boa(n)
      enddo
      read(1,*)nazms
      if ( nazms .gt. E_nazms ) go to 53         ! Skip if out of bounds
      do n = 1, nazms
         read(1,*)azm_boa(n)
      enddo
      close(1)
      ngeoms = nszas * nvzas * nazms
   endif

!  Linearization Control
!  ---------------------

!  Hardwiring or from 'Configfiles/FastV2p7_LinCont.cfg'

!  initialize

   which_gaspars = .false.
   do_GasAbs_WF  = .false.
   do_AOD_WF     = .false.
   do_profile_wfs = .false. ; npars  = 0 ; LFVARY = .false. ; LNVARY = 0
   do_surface_wfs = .false. ; nspars = 0
   
!  up to 5 gases, but only want ozone
!  12/29/15, Allow for possibility of aerosol Jacobians

   which_gaspars(1) = .true.
   do_GasAbs_WF     = .true.  !  These together
   if ( do_aerosols ) do_AOD_WF = .true.
   do_profile_wfs   = do_GasAbs_WF
   do_surface_wfs = .true. ; if ( do_surface_wfs ) nspars = 1

!  CHECK SECTION
!  -------------

!  Check Configuration inputs against given dimensions

53 continue

   if ( do_obsgeoms ) then
      if ( ngeoms .gt. E_ngeoms ) then
         b_fail = .true. ;  b_nmessages =  b_nmessages + 1
         b_messages(b_nmessages) = 'Input ngeoms > Allowed dimension' 
      endif
   else
      if ( nszas .gt. E_nszas ) then
         b_fail = .true. ;  b_nmessages =  b_nmessages + 1
         b_messages(b_nmessages) = 'Input nszas > Allowed dimension' 
      endif
      if ( nvzas .gt. E_nvzas ) then
         b_fail = .true. ;  b_nmessages =  b_nmessages + 1
         b_messages(b_nmessages) = 'Input nvzas > Allowed dimension' 
      endif
      if ( nazms .gt. E_nazms ) then
         b_fail = .true. ;  b_nmessages =  b_nmessages + 1
         b_messages(b_nmessages) = 'Input nazms > Allowed dimension' 
      endif
   endif

   if ( ndat .gt. E_NDAT ) then
      b_fail = .true. ;  b_nmessages =  b_nmessages + 1
      b_messages(b_nmessages) = 'Input ndat > Allowed dimension' 
   endif

   if ( ngases .gt. E_NGASES ) then
      b_fail = .true. ;  b_nmessages =  b_nmessages + 1
      b_messages(b_nmessages) = 'Input ngases > Allowed dimension' 
   endif

   if ( nlayers .gt. E_NLAYERS ) then
      b_fail = .true. ;  b_nmessages =  b_nmessages + 1
      b_messages(b_nmessages) = 'Input nlayers > Allowed dimension' 
   endif

!  Check Configuration inputs against VLIDORT dimensions

   if ( nstreams .gt. MAXSTREAMS ) then
      b_fail = .true. ;  b_nmessages =  b_nmessages + 1
      b_messages(b_nmessages) = 'Input nstreams > VLIDORT MAXSTREAMS dimension' 
   endif

   if ( nlayers .gt. MAXLAYERS ) then
      b_fail = .true. ;  b_nmessages =  b_nmessages + 1
      b_messages(b_nmessages) = 'Input nlayers > VLIDORT MAXLAYERS dimension' 
   endif

   if ( do_obsgeoms ) then
      if ( ngeoms .gt. MAX_GEOMETRIES ) then
         b_fail = .true. ;  b_nmessages =  b_nmessages + 1
         b_messages(b_nmessages) = 'Input ngeoms > VLIDORT MAX_GEOMETRIES dimension' 
      endif
   else
      if ( nszas .gt. MAXBEAMS ) then
         b_fail = .true. ;  b_nmessages =  b_nmessages + 1
         b_messages(b_nmessages) = 'Input nszas > VLIDORT MAXBEAMS dimension' 
      endif
      if ( nvzas .gt. MAX_USER_STREAMS ) then
         b_fail = .true. ;  b_nmessages =  b_nmessages + 1
         b_messages(b_nmessages) = 'Input nvzas > VLIDORT MAX_USER_STREAMS dimension' 
      endif
      if ( nazms .gt. MAX_USER_RELAZMS ) then
         b_fail = .true. ;  b_nmessages =  b_nmessages + 1
         b_messages(b_nmessages) = 'Input nazms > VLIDORT MAX_USER_RELAZMS dimension' 
      endif
   endif

!  Errors, finish

   if ( b_fail ) go to 69

!  Scattering angle calculation added, 3/30/16. Needed for Raman calculation

   if ( do_obsgeoms ) then
      do n = 1, ngeoms
         cx0 = cos(dtor*sza_boa(n)) ; sx0 = sin(dtor*sza_boa(n))
         cx1 = cos(dtor*vza_boa(n)) ; sx1 = sin(dtor*vza_boa(n))
         scatang_boa(n) = acos ( cx0 * cx1 - sx0 * sx1 * cos(dtor*azm_boa(n)) ) / dtor
      enddo
   else
      do i = 1, nszas
         cx0 = cos(dtor*sza_boa(i)) ; sx0 = sin(dtor*sza_boa(i))
         do j = 1, nvzas
            cx1 = cos(dtor*vza_boa(j)) ; sx1 = sin(dtor*vza_boa(j))
            do k = 1, nazms
               n = nvzas * nazms * ( i - 1 ) + nazms * ( j - 1 ) + k
               scatang_boa(n) = acos ( cx0 * cx1 - sx0 * sx1 * cos(dtor*azm_boa(k)) ) / dtor
            enddo
         enddo
      enddo
   endif

!  Initial CPU time setting
!  ------------------------

   if ( Monitor_CPU ) call cpu_time(e1)

!  1. Create optical properties.
!     ==========================

!  Wavelengths/Wavenumbers

   do n = 1, ndat
      n1 = ndat+1-n
      if (.not. use_hitran) then
         lambdas(n)  = wstart + real(n-1,kind=dpm) * wresol
         wavnums(n1) = 1.0d+07/lambdas(n)
      else
         wavnums(n)  = wstart + real(n-1,kind=dpm) * wresol
         lambdas(n1) = 1.0d+07/wavnums(n)
      endif
   enddo

!  Get OMI Coarse wavelengths and Solar Spectrum
!   - Use the base window limits to collect OMI points
!   - OMI_SolarSpectrum_Ch*.dat file extracted from Averaged data ("omisol_avg..."

   if ( OMI_Window .eq. 1 ) then
      open(23,file= 'OMI_SolarSpectrum_Ch1.dat', status = 'old')
      read(23,*)nwav_ch1
      w = 0
      do i = 1, nwav_ch1
         read(23,*) idum, lam, sol
         if ( lam.ge.wstart_base.and.lam.le.wfinis_base ) then
            w = w + 1 ; lambdas_c(w) = lam ;  OMI_solarflux_c(w) = sol
         endif
      enddo
      ndat_c = w
      close(23)
   else
      open(23,file= 'OMI_SolarSpectrum_Ch2.dat', status = 'old')
      read(23,*)nwav_ch2
      w = 0
      do i = 1, nwav_ch2
         read(23,*) idum, lam, sol
         if ( lam.ge.wstart_base.and.lam.le.wfinis_base ) then
            w = w + 1 ; lambdas_c(w) = lam ;  OMI_solarflux_c(w) = sol
         endif
      enddo
      ndat_c = w
      close(23)
   endif

!  Actual creation of optical properties

   if ( Monitor_CPU ) call cpu_time(e2)

!  Fd testing. Must create anew if FDINDEX not zero
!   -- That way, do not disturb the original Dumps

   nstep = 1
   FDINDEX = 0
   if (FDINDEX.ne.0 ) then
      do_Dump1 = .true. ; do_Create1 = .true.
   else
      do_Dump1 = do_Dump ; do_Create1 = do_Create
   endif

!   write(*,*)'Just before creation: Dump/Create Flags = ', do_Dump1, do_Create1 ; pause
!  Call to create properties. (Ncloud_layers added, 2/17/16)
!  New arguments, 29 March 2016
!    ** derivative of Gas absorption w.r.t. Log(VMR)
!    ** input OMI_Window (MW1, MW2, MW3), for setting the albedo
!    ** Temperatures/Aircolumns outputs added 3/30/16
!    ** Rayleigh Xsec output added 4/1/16
!    ** RingSF+Albedo output added 5/4/16

!  New argument, 22 September 2016. do_PGE_TES

   if ( do_Dump1 .or. do_Create1 ) then

      call  FastV2p7_CreateProps_Plus_Master &
       ( E_nlayers, E_ndat, E_nmoms_all, E_ngases, E_npars, Geocape_Path,               & ! Input
         AIRSOMI_DataPath, Select_index, OMI_Channel, OMI_Window, AIRSOMI_SceneIndex,   & ! Input
         do_PGE_TES, do_aerosols, use_hitran, do_GasAbs_WF, do_AOD_WF,                  & ! Input
         ndat, nmoms_2ns, ngases, which_gases, lambdas, wavnums, nlayers,               & ! input
         ncloud_layers, nmoms_all, npars, taug, tausp, taudp, omega, depol, RayXsec,    & ! Main  Output
         solar_flux, aergreeks, albedo, heights, fr, fa, Temperatures, AirColumns,      & ! Main  Output
         L_taudp, L_omega, L_fr, L_fa, dtaug_dLogv,                                     & ! Main  Output
         EarthRadius, Angles, CldFraction, CldTopPressure, CldTopAlbedo,                & ! Other Output
         Ring_Scaling, Albedo_0, Albedo_1, Lambda_Ref_MW2,                              & ! Other Output
         b_fail, b_nmessages, b_messages(1) )                                             ! output

      if ( b_fail ) go to 69

   endif

!  Dump file header. New 2/18/16

   ndatr = ( (ndat-1)/nstep) + 1
   call mkchar(4,c1d,c2d,c3d,c4,ndatr) ; if (c4(1:1).eq.'0')c4(1:1) = '_'
   call mkchar(1,c1,c2d,c3d,c4d,OMI_Window)
   Dumpfile = 'LinDump_FastV2p7_Properties_MW'//C1//'_'//C4//'_'//Transept

!  Dump (linearized, just one Profile Jacobian)
!      Ncloud_layers and Cloud information added, 2/17/16)
!mick fix 2/26/2016 -  added tausp to dump file ; expanded number of L_taudp & L_omega to dump file
!   29 March 2016, Added derivative of Gas absorption w.r.t. Log(VMR) to DUMP
!      --- temperatures and Aircolumns added, 3/30/16. Rayleigh Xsec added 4/1/16
!      --- RingSF+Albedo output added 5/4/16

    if ( do_Dump1 ) then
       open(1,file=Trim(Dumpfile)//'.dat',status='unknown')
       write(1,*)nlayers, ndatr, nmoms_all, npars, nspars, sza_boa(1), vza_boa(1), azm_boa(1)
       write(1,*)ncloud_layers, CldFraction, CldTopPressure, CldTopAlbedo, Ring_Scaling, Albedo_0, Albedo_1
       write(1,*)heights(0)
       do n = 1, nlayers
          write(1,*)heights(n),temperatures(n),Aircolumns(n)
       enddo
       if ( do_aerosols ) then
          DO L = 0, E_nmoms_all
            write(1,'(1p60e12.4)')((aergreeks(1,L,k,m),k=1,5),m=1,6),((aergreeks(2,L,k,m),k=1,5),m=1,6)
          enddo
       endif
       do i = 1, ndatr, nstep
          if ( mod(i,100).eq. 0 ) write(*,*)'Dumping original data # ',i
          write(1,*)i,lambdas(i),wavnums(i),albedo(i),depol(i),RayXsec(i), solar_flux(i)
!          write(1000,*)i,lambdas(i),sum(taudp(1:nlayers,i))  ! check
          do n = 1, nlayers
             !write(1,*)n,taug(n,i),taudp(n,i),omega(n,i),L_taudp(1,n,i),L_omega(1,n,i),fr(n,i),fa(1:5,n,i),&
             !          (tausp(n,i,j),j=1,ngases)
             write(1,*)n,taug(n,i),taudp(n,i),omega(n,i),fr(n,i),fa(1:5,n,i)
             write(1,*)(tausp(n,i,j),j=1,ngases)
             write(1,*)(L_taudp(j,n,i),L_omega(j,n,i),j=1,npars)
             write(1,*)(dtaug_dLogv(n,1,i,j),j=1,ngases),(dtaug_dLogv(n,2,i,j),j=1,ngases)
          enddo
       enddo
       close(1)
       stop'Finished Writing LinDump'
    else if ( .not. do_Create1 ) then
       open(1,file=Trim(Dumpfile)//'.dat',status='old')
       write(*,'(1x,2a)') 'Reading Dump file: ',Trim(Dumpfile)//'.dat'
       read(1,*)nlayers, ndat, nmoms_all, npars, nspars, sza_boa(1), vza_boa(1), azm_boa(1)
       read(1,*)ncloud_layers, CldFraction, CldTopPressure, CldTopAlbedo, Ring_Scaling, Albedo_0, Albedo_1
       read(1,*)heights(0)
       do n = 1, nlayers
          read(1,*)heights(n),temperatures(n),Aircolumns(n)
       enddo
       if ( do_aerosols ) then
          DO L = 0, E_nmoms_all
            read(1,'(1p60e12.4)')((aergreeks(1,L,k,m),k=1,5),m=1,6),((aergreeks(2,L,k,m),k=1,5),m=1,6)
          enddo
       endif
       do i = 1, ndat, nstep
          read(1,*)idum,lambdas(i),wavnums(i),albedo(i),depol(i),RayXsec(i), solar_flux(i)
          if ( mod(i,100).eq. 0 ) write(*,*)'Dump File: Done reading line #',i
          do n = 1, nlayers
            !read(1,*)ndum,taug(n,i),taudp(n,i),omega(n,i),L_taudp(1,n,i),L_omega(1,n,i),fr(n,i),fa(1:5,n,i),&
            !         (tausp(n,i,j),j=1,ngases)
            read(1,*)ndum,taug(n,i),taudp(n,i),omega(n,i),fr(n,i),fa(1:5,n,i)
            read(1,*)(tausp(n,i,j),j=1,ngases)
            read(1,*)(L_taudp(j,n,i),L_omega(j,n,i),j=1,npars)
            read(1,*)(dtaug_dLogv(n,1,i,j),j=1,ngases),(dtaug_dLogv(n,2,i,j),j=1,ngases)
          enddo
       enddo
       close(1)
       write(*,*)'Finished Reading Dump'
    endif

!  5 gases, but only want ozone as first gas (Special case)

   do i = 1, ndat
      do n = 1, nlayers
         tauO3(n,i) = tausp(n,i,1)
      enddo
   enddo

!   write(*,*)'Just after creation'

   if ( Monitor_CPU ) then
      call cpu_time(e3) ; createproptime = e3-e2
   endif

!  Check NLAYERS against given dimensions

   if ( nlayers .gt. E_NLAYERS ) then
      b_fail = .true. ;  b_nmessages =  b_nmessages + 1
      b_messages(b_nmessages) = 'Input nlayers > Allowed dimension' 
   endif

!  Check NLAYERS against VLIDORT dimensions

   if ( nlayers .gt. MAXLAYERS ) then
      b_fail = .true. ;  b_nmessages =  b_nmessages + 1
      b_messages(b_nmessages) = 'Input nlayers > VLIDORT MAXLAYERS dimension' 
   endif

!  Define NLEVELS

   nlevels = nlayers + 1

!  Exception handling for Creation routine

   if ( b_fail ) go to 69

!  Total property output file
!     TAUDP, TAUG, TAURAY, TAUAER, 

   if ( do_debug_ExactTool ) then
      open(1,file='Debug_stuff/FastV2p7_Properties_debug.dat',status='unknown')
      do i = 1, ndat
         taugas = 0.0d0 ; tauray = 0.0d0 ; tautot = 0.0d0
         do n = 1, nlayers
            tauray = tauray + fr(n,i) * omega(n,i) * taudp(n,i)
            taugas = taugas + taug(n,i)
            tautot = tautot + taudp(n,i)
         enddo
         tauaer = tautot - taugas - tauray
         if (.not. use_hitran) then
            write(1,'(1p5e17.7)')lambdas(i),tautot,taugas,tauray,tauaer
         else
            write(1,'(1p5e17.7)')lambdas(ndat-i+1),tautot,taugas,tauray,tauaer
         endif
      enddo
      close(1)
   endif

!  Set L_fr and L_fa = zero, as only have trace-gas absorption Jacobians

   L_fr = 0.0d0 ; L_fa = 0.0d0

!  Check

   if ( do_Dump .or. do_Create ) then
      write(*,*)'..Finished Creation in the Linearized (LPS) Exact tool....'
!   pause     '..Finished Creation in the Exact tool....'
   else
      write(*,*)'..Finished Reading Creation-Dump in the Linearized (LPS) Exact tool....'
   endif

!  Additional linearization inputs

   if ( do_profile_wfs ) then
      LFvary(1:nlayers) = .true. ; LNVary(1:nlayers) = npars
   endif

!  New for 2/17/16. IPA loop
!  =========================

!  fake the Cloud fraction for debugging.....
!   CldFraction = 0.01

!  FD testing special

   FDindex = 20 ; FDindex = 0

!  Scene bookkeeping
!  -----------------

   nscenes = 1 ; do_clearsky = .true. ; Output_index = 1 ; Charscene = 'Clear' 
   if ( CldFraction.lt.CldFraction_Threshold_Upper ) then
      if ( CldFraction.gt.CldFraction_Threshold_Lower ) then
        nscenes = 2 ; do_clearsky(2) = .false. ; output_index(2) = 2 ; Charscene(2) = 'Cloud'
      endif
   else
      nscenes = 1 ; do_clearsky(1) = .false. ; output_index(1) = 2 ; Charscene(1) = 'Cloud'
   endif
   ClrFraction = 1.0_dpm - CldFraction
 
!  Set File headers
!  ----------------

!  Timing for Write-up

   if (monitor_CPU) call cpu_time(e2)

!  general

   Fileheader0 = '_D00_Aer_PP_Latg_S00V00A000_MW0_CF000_Tr000'
!                 1234567890123456789012345678901234567890123
   if ( do_obsgeoms ) Fileheader0(13:16) = 'Obsg' 
   if ( nstreams  .gt.9 ) write( Fileheader0(3:4),'(I2)' ) nstreams
   if ( nstreams  .le.9 ) write( Fileheader0(4:4),'(I1)' ) nstreams
   if (      do_enhanced_ps ) Fileheader0(10:11) = 'ES' 
   if ( .not.do_enhanced_ps ) Fileheader0(10:11) = 'RS'
   if ( .not.do_aerosols ) Fileheader0(6:8) = 'Ray'
   if ( do_Continuum ) Fileheader0(6:8)     = 'NoG'
   Cstokes = 'Scalar' ;  if ( nstokes.gt.1) Cstokes = 'Vector' 

!  OMI-specific Labels

   write(Fileheader0(31:31),'(I1)')OMI_Window
   call mkchar(3,c1d,c2d,Fileheader0(35:37),c4d,NINT(100.0_dpm*CldFraction))
   Fileheader0(39:43) = Transept

!  geometry Labels

   if ( do_obsgeoms ) then
     do v = 1, ngeoms
        Fileheader = Fileheader0 ; c2d = '  ' ; c3d = '   '
        call mkchar(2,c1d,Fileheader(19:20),c3d,c4d,INT(sza_boa(v)))
        call mkchar(2,c1d,Fileheader(22:23),c3d,c4d,INT(vza_boa(v)))
        call mkchar(3,c1d,c2d,Fileheader(25:27),c4d,INT(azm_boa(v)))
        FileHeader_save(v) = 'Exact'//Cstokes//Trim(Fileheader)
     enddo
   else
     do l = 1, nszas
      do j = 1, nvzas
       do k = 1, nazms
        Fileheader = Fileheader0 ; c2d = '  ' ; c3d = '   '
        v = nazms * nvzas * ( l - 1 ) + nazms * ( j - 1 ) + k
        call mkchar(2,c1d,Fileheader(19:20),c3d,c4d,INT(sza_boa(l)))
        call mkchar(2,c1d,Fileheader(22:23),c3d,c4d,INT(vza_boa(j)))
        call mkchar(3,c1d,c2d,Fileheader(25:27),c4d,INT(azm_boa(k)))
        FileHeader_save(v) = 'Exact'//Cstokes//Trim(Fileheader)
       enddo
      enddo
     enddo
   endif

!  File names. Revised 5/6/16

   do v = 1, ngeoms
!  OMI (Coarse-grid) output. Now Includes Timing.
      OMIOutput_AllResults(v)  = Trim(Results_path)//'OMI_Retrieval_Output_'//Trim(Fileheader_save(v))//'.OUT'
      OMIOutput_Timing(v)      = Trim(Results_path)//'OMI_Timing____Output_'//Trim(Fileheader_save(v))//'.TIME'
      if ( do_LPJ_ContourPlot ) then
         OMIOutput_LPJContour(v)  = Trim(Results_path)//'OMI_LPContour_Output_'//Trim(Fileheader_save(v))//'.PLOT'
      endif
!  Fine-grid output and timing. Optional
      if ( do_Output_FineGrid ) then
         fileOutput1_FineGrid(v) = Trim(Results_path)//'RD_'//Trim(Fileheader_save(v))//'.Out'
         fileOutput2_FineGrid(v) = Trim(Results_path)//'LP_'//Trim(Fileheader_save(v))//'.Out'
         fileOutput3_FineGrid(v) = Trim(Results_path)//'LS_'//Trim(Fileheader_save(v))//'.Out'
         fileOutput4_FineGrid(v) = Trim(Results_path)//'SW_'//Trim(Fileheader_save(v))//'.Out'
      endif  
   enddo

!  timing for filename preparation

   if (monitor_CPU) then
      call cpu_time(e3) ; Writetime = Writetime + e3-e2
   endif

!  Prepare for main loop
!  ---------------------

!  Local wavelengths

!   if ( use_hitran ) then 
!     do i = 1, ndat
!       local_lambdas(i) = lambdas(ndat+1-i)
!     enddo
!   else
!     do i = 1, ndat
!       local_lambdas(i) = lambdas(i)
!     enddo
!   enddo

!  Zero the Fine-Grid output

   Intensity_Exact    = 0.0_dpm
   Intensity_LD_Exact = 0.0_dpm
   Intensity_2S_Exact = 0.0_dpm
   Intensity_FO_Exact = 0.0_dpm
   LP_Jacobians_Exact    = 0.0_dpm ; ScatWts_Exact    = 0.0_dpm ; LS_Jacobians_Exact    = 0.0_dpm
   LP_Jacobians_LD_Exact = 0.0_dpm ; ScatWts_LD_Exact = 0.0_dpm ; LS_Jacobians_LD_Exact = 0.0_dpm
   LP_Jacobians_FO_Exact = 0.0_dpm ; ScatWts_FO_Exact = 0.0_dpm ; LS_Jacobians_FO_Exact = 0.0_dpm
   LP_Jacobians_2S_Exact = 0.0_dpm ; ScatWts_2S_Exact = 0.0_dpm ; LS_Jacobians_2S_Exact = 0.0_dpm

!Pushkar's edit, get total gas optical depth for each line to check for jumps
!mick fix 2/26/2016 - limit sum to 1:nlayers

   outgastau = 0.0_dpm
   do i = 1, ndat
      !outgastau(i) = sum(taug(:,i))
      outgastau(i) = sum(taug(1:nlayers,i))
   enddo

!  Skip fine calculation

   if ( do_Skip_Fine ) then
       intensity_Exact(1:ndat,1,1) = solar_flux(1:ndat) 
       intensity_Exact(1:ndat,1,2) = 2.0*solar_flux(1:ndat)
       write(*,*)'WARNING!!!!    Skipping fine calculation to test Raman' 
       go to 67
   endif

!  IPA loop, RT calculations
!  =========================

   do iscene = 1, nscenes
!   do iscene = 1, 1

!  write index

     Iout = Output_index(iscene)
     write(*,'(/A/)')' ** Doing Fine-grid calculation for '//Charscene(iscene)//' Scenario ---'

!  Re-set VLIDORT initialization flag, each scene. Important !!!!!!!

     DO_VLidort_initialize = .false.

!  Local albedo and layers

     local_albedo = 0.0_dpm
     if ( do_clearsky(iscene) ) then
       local_nlayers = nlayers
       local_albedo(1:ndat) = albedo(1:ndat)
     else
       local_nlayers = ncloud_layers
       local_albedo(1:ndat) = CldTopAlbedo
     endif

!  No surface Wfs for the Cloudy scene. IMPORTANT !!!!!!!

     local_do_surf_wfs = do_surface_wfs
     if ( .not.do_clearsky(iscene) ) local_do_surf_wfs = .false.
     
!  Perform Calculation

     call FastV2p7_Exact_LPS_Master ( Monitor_CPU, &
        E_nlayers, E_ndat, E_nmoms_all, E_npars, E_nspars,                   & ! Input dimensions
        E_ngeoms, E_nszas, E_nvzas, E_nazms,                                 & ! Input dimensions
        DO_VLidort_initialize, do_obsgeoms, do_aerosols, use_hitran,         & ! input flags
        do_optional_2stream, bvpindex, do_pdinverse, do_enhanced_ps,         & ! input flags (RT)
        local_do_surf_wfs, do_profile_wfs, npars, nspars, lfvary,lnvary,     & ! Linearization inputs
        nstreams, local_nlayers, ndat, nmoms_all,                            & ! Input numbers
        ngeoms, nstr2, nstokes, nszas, nvzas, nazms,                         & ! Input numbers (Geometry)
        sza_boa, vza_boa, azm_boa, eradius, heights, solar_flux,             & ! input geometry/heights
        lambdas, taug, taudp, omega, depol, aergreeks, local_albedo, fr, fa, & ! input optical
        L_taudp, L_omega, L_fr, L_fa,                                        & ! input optical
        intensity_Exact(:,:,Iout),    LP_Jacobians_exact(:,:,:,:,Iout),    & ! TOTAL   results
        ScatWts_Exact(:,:,:,Iout),    LS_Jacobians_exact(:,:,:,Iout),      & ! TOTAL   results
        intensity_LD_Exact(:,:,Iout), LP_Jacobians_LD_Exact(:,:,:,:,Iout), & ! VLIDORT MS results
        ScatWts_LD_Exact(:,:,:,Iout), LS_Jacobians_LD_Exact(:,:,:,Iout),   & ! VLIDORT MS results
        intensity_FO_Exact(:,:,Iout), LP_Jacobians_FO_Exact(:,:,:,:,Iout), & ! FO      SS results
        ScatWts_FO_Exact(:,:,:,Iout), LS_Jacobians_FO_Exact(:,:,:,Iout),   & ! FO      SS results
        intensity_2S_Exact(:,:,Iout), LP_Jacobians_2S_Exact(:,:,:,:,Iout), & ! 2STREAM MS results
        ScatWts_2S_Exact(:,:,:,Iout), LS_Jacobians_2S_Exact(:,:,:,Iout),   & ! 2STREAM MS results
        Exacttimes(:,Iout), ExactRuntime(Iout), b_fail, b_nmessages, b_messages )  ! Timings and status

!  Exception handling for Master module

     if ( b_fail ) then
       b_messages(b_nmessages+1) = 'LPS Exact Master failed for '//Charscene(iscene)//' Scenario'
       b_nmessages = b_nmessages + 1
       go to 69
     endif

!  Working Only for 1 jacobian (q = 1, qs = 1)

     q = 1 ; qs = 1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   Apply the IPA at this stage. CODE REMOVED 5/6/16
!      --  To Recover, See the file FineGrid_IPACalc_Remnant.f90
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  convert to Level Jacobians ( w.r.t. LOG(Vmr) ), using the chain rule formula
!    Added 3/29/16. Modeled after the MUSICQA code.

!    9/22/16 Modified form of Mapping transformations.
!     -- Designed to work with both TES and Original systems.

!mick fix 10/11/2016 - moved defining of NLEVELS following dump file read
     !nlevels = nlayers + 1
     do g = 1, ngases
        do v = 1, ngeoms
           do n = 0, nlayers
              n1 = n + 1
              if ( n.eq.0 ) then
                 do w = 1, ndat
                    if ( taug(n1,w).gt.0.0_dpm ) then
                       kn1 = LP_Jacobians_exact(q,n1,w,v,iout) ; kn1_gasabs = kn1 ! / taug(n1,w)
                       LP_LogV_Jacobians_Exact(g,n1,w,v,iout) = kn1_gasabs * dtaug_dLogv(n1,1,w,g)
                    endif
                 enddo
              else if ( n.eq.nlayers ) then
                 do w = 1, ndat
                    if ( taug(n,w).gt.0.0_dpm ) then
                       kn2 = LP_Jacobians_exact(q,n,w,v,iout) ; kn2_gasabs = kn2 ! / taug(n,w)
                       LP_LogV_Jacobians_Exact(g,n1,w,v,iout) = kn2_gasabs * dtaug_dLogv(n,2,w,g)
                    endif
                 enddo
              else
                 do w = 1, ndat
                    if ( taug(n1,w).gt.0.0_dpm .and. taug(n,w).gt.0.0_dpm ) then
                       kn1 = LP_Jacobians_exact(q,n1,w,v,iout) ; kn1_gasabs = kn1 ! / taug(n1,w)
                       kn2 = LP_Jacobians_exact(q,n,w,v,iout)  ; kn2_gasabs = kn2 ! / taug(n,w)
                       LP_LogV_Jacobians_Exact(g,n1,w,v,iout) = kn1_gasabs * dtaug_dLogv(n1,1,w,g) &
                                                              + kn2_gasabs * dtaug_dLogv(n,2,w,g)
                    endif
                 enddo
              endif
           enddo
        enddo
     enddo

!  End IPA loop

   enddo

!  debug FD on point 45
!   do n1 = 1, nlevels
!      write(555,*)n1,Intensity_Exact(45,1,1:2),LP_LogV_Jacobians_Exact(1,n1,45,1,1:2)
!   enddo

!  tempo skip

67 continue

!  short test - avoid Convolution and Raman
!  ----------------------------------------

!  Code added 4/1/16.

   if ( do_shorttest ) go to 44

!  Convolution Section
!  ===================

 !  Start timing

   if (monitor_CPU) call cpu_time(e2)

!  Auxiliary Jacobians for Albedo coefficients (5/6/16). Need to be initialized here.

   JExact_AlbCoeffs_C = 0.0_dpm 

!  Introduce the OMI ILS Convolution, 5/4/16
!  -----------------------------------------

!  Three stages : (1) Extract ILS (2) Spline ILS to Fine grid (3) Convolve

   if ( OMIILS_Convolution_Method ) then

!  (1) Get the slit functions
!        ** E_ndat_ILS (dimension) and ndat_ILS (number of ILS points) should both be 505.
!        ** Perform excpetion handling after the call

      Call Extract_OMIILS &
         ( E_ndat_c, E_ndat_ILS, ndat_c, ndat_ILS, lambdas_c, &
           OMI_PixelIndex, OMI_Channel, OMI_Window,           &
           ILS_wavs, ILS_Values, fail_ILS, message_ILS )

      if ( fail_ILS ) then
         b_messages(b_nmessages+1) = Trim(Message_ILS)
         b_messages(b_nmessages+2) = 'Slit function not found : Extract_OMIILS failed'
         b_nmessages = b_nmessages + 2
         b_Fail = .true. ; go to 69
      endif

!  (2) Spline the ILS

      Call Spline_OMIILS &
         ( E_ndat, E_ndat_c, E_ndat_ILS, ndat, lambdas, ndat_c,      &
           ndat_ILS, ILS_wavs, ILS_Values, ILS_offsets, ILS_npoints, &
           fail_ILS, message_ILS )

      if ( fail_ILS ) then
         b_messages(b_nmessages+1) = Trim(Message_ILS)
         b_messages(b_nmessages+2) = 'Slit function Splining to fine grid : Spline_OMIILS failed'
         b_nmessages = b_nmessages + 2
         b_Fail = .true. ; go to 69
      endif

!  (3) Convolve Sun-normalized intensities

      do iscene = 1, nscenes
         Iout = Output_index(iscene) ; qs = 1
         Call Convolve_OMIILS ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, ngeoms, &
                                ILS_Values, ILS_offsets, ILS_npoints,         &
                                lambdas, intensity_Exact(:,1:ngeoms,iout), IExact_C(:,1:ngeoms,iout) )
         do g = 1, ngases
            do n = 1, nlevels
               Call Convolve_OMIILS ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, ngeoms, &
                                      ILS_Values, ILS_offsets, ILS_npoints,         &
                                      lambdas, LP_LogV_Jacobians_Exact(g,n,:,1:ngeoms,iout), LogV_JExact_C(g,n,:,1:ngeoms,iout) )
            enddo
         enddo
         Call Convolve_OMIILS ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, ngeoms, &
                                ILS_Values, ILS_offsets, ILS_npoints,         &
                                lambdas, LS_Jacobians_Exact(qs,:,1:ngeoms,iout), LS_JExact_C(qs,:,1:ngeoms,iout) )
      enddo

!   ---> Albedo-coefficient Jacobians. 5/5/16
   
      do iscene = 1, nscenes
         Iout = Output_index(iscene) ; qs = 1
         JExact_AlbCoeffs_C(1,1:ndat_c,1:ngeoms,iout) = LS_JExact_C(qs,1:ndat_c,1:ngeoms,iout)
         if ( OMI_Window .eq. 2 ) then
            do w = 1, ndat_c
               albfac = ( lambdas_c(w) / Lambda_Ref_MW2 ) - 1.0_dpm
               JExact_AlbCoeffs_C(2,w,1:ngeoms,iout) = albfac * LS_JExact_C(qs,w,1:ngeoms,iout)
            enddo
         endif
      enddo

!  finish OMI ILS convolutione

   endif

!  NON OMI-ILS Convolution, using Asymmetric Gaussians
!  ---------------------------------------------------

   if ( .not. OMIILS_Convolution_Method ) then

!  Either : Convolve fine-grid sun-normalized intensities (Original method, also used by OMI-ILS)
!  Or     : Use Alternative more physical method withfine-grid solar spectrum from SAO

!  Compute Half-widths. This is a temporary fiddle-calculation

      hwhm_c(1) = lambdas_c(2) - lambdas_c(1)
      do w = 2, ndat_c - 1
         hwhm_c(w) = 0.5d0* ( (lambdas_c(w+1)-lambdas_c(w)) + (lambdas_c(w)-lambdas_c(w-1)) )
      enddo
      hwhm_c(ndat_c) = lambdas_c(ndat_c) - lambdas_c(ndat_c-1)
      asym_C = 0.0d0

!  Original method - Convolve Sun-normalized variables

      if ( Original_Convolution_Method ) then

!  ---> Convolve Sun-normalized intensities/Jacobians

         do iscene = 1, nscenes
            Iout = Output_index(iscene); qs = 1
            Call asym_gauss_f2c (lambdas(1:ndat), intensity_Exact(1:ndat,1:ngeoms,iout), ndat, ngeoms, ndat_c, &
                                 hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), IExact_C(1:ndat_c,1:ngeoms,iout) )
            do g = 1, ngases
               do n = 1, nlevels
                  Call asym_gauss_f2c (lambdas(1:ndat), LP_LogV_Jacobians_Exact(g,n,1:ndat,1:ngeoms,iout), ndat, ngeoms, ndat_c, &
                              hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), LogV_JExact_C(g,n,1:ndat_c,1:ngeoms,iout) )
               enddo
            enddo
            Call asym_gauss_f2c (lambdas(1:ndat), LS_Jacobians_Exact(qs,1:ndat,1:ngeoms,iout), ndat, ngeoms, ndat_c, &
                                 hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), LS_JExact_C(qs,1:ndat_c,1:ngeoms,iout) )
         enddo

!   ---> Albedo-coefficient Jacobians. 5/5/16
   
         do iscene = 1, nscenes
            Iout = Output_index(iscene) ; qs = 1
            JExact_AlbCoeffs_C(1,1:ndat_c,1:ngeoms,iout) = LS_JExact_C(qs,1:ndat_c,1:ngeoms,iout)
            if ( OMI_Window .eq. 2 ) then
               do w = 1, ndat_c
                  albfac = ( lambdas_c(w) / Lambda_Ref_MW2 ) - 1.0_dpm
                  JExact_AlbCoeffs_C(2,w,1:ngeoms,iout) = albfac * LS_JExact_C(qs,w,1:ngeoms,iout)
               enddo
            endif
         enddo

      endif

!  Alternative Method - convolution using fine-grid radiances and solar spectrum

      if ( .not. Original_Convolution_Method ) then

!  ---> Convolution of fine-grid solar spectrum

         Call asym_gauss_f2c_1 (lambdas(1:ndat), Solar_flux(1:ndat), ndat, ndat_c, &
                                hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), Solarflux_c(1:ndat_c) )

!  ---> Form fine-grid radiances/Jacobians, by multiplying sun-normalized RT output with fine grid solar spectrum

         do iscene = 1, nscenes
            Iout = Output_index(iscene) ; qs = 1
            do i = 1, ndat
               intensity_Exact(i,1:ngeoms,iout) = intensity_Exact(i,1:ngeoms,iout) * Solar_flux(i)
               do g = 1, ngases
                  do n = 1, nlevels
                     LP_LogV_Jacobians_Exact(g,n,i,1:ngeoms,iout) = LP_LogV_Jacobians_Exact(g,n,i,1:ngeoms,iout) * Solar_flux(i)
                  enddo
               enddo
               LS_Jacobians_Exact(qs,i,1:ngeoms,iout) = LS_Jacobians_Exact(qs,i,1:ngeoms,iout) * Solar_flux(i)
            enddo
         enddo

!  ---> Convolution of fine-grid radiances and Jacobians

         do iscene = 1, nscenes
            Iout = Output_index(iscene) ; qs = 1
            Call asym_gauss_f2c (lambdas(1:ndat), intensity_Exact(1:ndat,1:ngeoms,iout), ndat, ngeoms, ndat_c, &
                                 hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), IExact_C(1:ndat_c,1:ngeoms,iout) )
            do g = 1, ngases
               do n = 1, nlevels
                  Call asym_gauss_f2c (lambdas(1:ndat), LP_LogV_Jacobians_Exact(g,n,1:ndat,1:ngeoms,iout), ndat, ngeoms, ndat_c, &
                            hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), LogV_JExact_C(g,n,1:ndat_c,1:ngeoms,iout) )
               enddo
            enddo
            Call asym_gauss_f2c (lambdas(1:ndat), LS_Jacobians_Exact(qs,1:ndat,1:ngeoms,iout), ndat, ngeoms, ndat_c, &
                              hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), LS_JExact_C(qs,1:ndat_c,1:ngeoms,iout) )
         enddo

!  ---> Restoration of convolved sun-normalized RT output

         do iscene = 1, nscenes
            Iout = Output_index(iscene) ; qs = 1
            do w = 1, ndat_c
               IExact_c(w,1:ngeoms,iout) = IExact_C(w,1:ngeoms,iout) / Solarflux_c(w)
               do g = 1, ngases
                  do n = 1, nlevels
                     LogV_JExact_C(g,n,w,1:ngeoms,iout) = LogV_JExact_C(g,n,w,1:ngeoms,iout) / Solarflux_c(w)
                  enddo
               enddo
               LS_JExact_C(qs,w,1:ngeoms,iout) = LS_JExact_C(qs,w,1:ngeoms,iout) / Solarflux_c(w)
            enddo
         enddo

!   ---> Albedo-coefficient Jacobians. 5/5/16
   
         do iscene = 1, nscenes
            Iout = Output_index(iscene) ; qs = 1
            JExact_Albcoeffs_C(1,1:ndat_c,1:ngeoms,iout) = LS_JExact_C(qs,1:ndat_c,1:ngeoms,iout)
            if ( OMI_Window .eq. 2 ) then
               do w = 1, ndat_c
                  albfac = ( lambdas_c(w) / Lambda_Ref_MW2 ) - 1.0_dpm
                  JExact_AlbCoeffs_C(2,w,1:ngeoms,iout) = albfac * LS_JExact_C(qs,w,1:ngeoms,iout)
               enddo
            endif
         enddo

!  End original vs. Alterative method...

      endif

!  End Non-OMI-ILS method

   endif

!  Timing on Convolution

   if (monitor_CPU) then
      call cpu_time(e3); Convolutiontime = e3-e2
   endif

!  Raman Filling in
!  ================

   if (monitor_CPU) call cpu_time(e2)

!  Sioris Raman calculation
!  ------------------------

   if ( Sioris_Raman_Method ) then

!  get the spectroscopy

      CALL Sioris_Raman_Getspec (RamanData)

!  Convolve layer optical depths (use reverse convolution subroutine)

      if ( OMIILS_Convolution_Method ) then
         Call Convolve_OMIILS_TR ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, nlayers, &
                                   ILS_Values, ILS_offsets, ILS_npoints,          &
                                   lambdas, taudp(1:nlayers,:), taudp_c(:,1:nlayers) )
      else
         Call asym_gauss_f2c_TR (lambdas(1:ndat), taudp(1:nlayers, 1:ndat), ndat, nlayers, ndat_c, &
                                  hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), taudp_c(1:ndat_c,1:nlayers) )
      endif

!  Call Sioris-Raman filliing factors, for both scenes
!   ---------- ONLY WORKS FOR  OBSERVATIONAL GEOMETRY

      do iscene = 1, nscenes

!  Upwelling flag

        SR_do_upwelling = .true.

!  write index

         Iout = Output_index(iscene)
         write(*,'(/A/)')' ** Doing Sioris-Raman calculation for '//Charscene(iscene)//' Scenario ---'

!  Local settings (scene-dependent)

         if ( do_clearsky(iscene) ) then
            SR_nlayers = nlayers       ; SR_albedo = sum(albedo(1:ndat))/dble(ndat)
         else
            SR_nlayers = ncloud_layers ; SR_albedo = CldTopAlbedo
         endif

!  Call for each geometry + exception handling

         do v = 1, ngeoms
            Call Sioris_Raman_Compute &
            ( SR_nlayers, ndat_c, sza_boa(v), vza_boa(v), scatang_boa(v), SR_albedo, SR_do_upwelling, &
              RamanData, temperatures(1:SR_nlayers), Aircolumns(1:SR_nlayers),                        &
              lambdas_c, OMI_solarflux_c(1:ndat_c), taudp_c(1:ndat_c,1:SR_nlayers),                   &
              OMI_filling_c(1:ndat_c,v,iout), Fail_SRaman, Message_SRaman)
            if ( Fail_SRaman ) then
               b_messages(b_nmessages+1) = Trim(Message_SRaman)
               b_messages(b_nmessages+2) = 'Sioris-Raman Call failed for '//Charscene(iscene)//' Scenario'
               b_nmessages = b_nmessages + 2
               b_Fail = .true. ; go to 69
            endif
         enddo

!  End Scene loop

      enddo

!  Debug
!      do w = 1, ndat_c
!         write(888,*)lambdas_c(w),OMI_solarflux_c(w),OMI_filling_c(w,1,1:2)
!      enddo

   endif

!  Alternative Raman calculation, Full calculation with LRRS 2p3
!  -------------------------------------------------------------

!  Programmed up, 4/1/16.

   if ( .not.Sioris_Raman_Method ) then

!  First, develop additional convolutions
!  ======================================

!  Convolve layer optical depths (USE TRANSPOSE, DO NOT USE REVERSE.......)

      if ( OMIILS_Convolution_Method ) then
         Call Convolve_OMIILS_T ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, nlayers, &
                                  ILS_Values, ILS_offsets, ILS_npoints,          &
                                  lambdas, taudp(1:nlayers,:), deltau_c(1:nlayers,:) )
      else
         Call asym_gauss_f2c_T (lambdas(1:ndat), taudp(1:nlayers, 1:ndat), ndat, nlayers, ndat_c, &
                                hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), deltau_c(1:nlayers,1:ndat_c) )
      endif

!  Convolution of depolarization and rayleigh cross-sections

      if ( OMIILS_Convolution_Method ) then
         Call Convolve_OMIILS_1 ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, &
                                  ILS_Values, ILS_offsets, ILS_npoints, &
                                  lambdas, depol, depol_c )
         Call Convolve_OMIILS_1 ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, &
                                  ILS_Values, ILS_offsets, ILS_npoints, &
                                  lambdas, RayXsec, RayXsec_c )
      else
         Call asym_gauss_f2c_1 (lambdas(1:ndat), depol(1:ndat), ndat, ndat_c, &
                                hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), depol_c(1:ndat_c) )
         Call asym_gauss_f2c_1 (lambdas(1:ndat), RayXsec(1:ndat), ndat, ndat_c, &
                                hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), RayXsec_c(1:ndat_c) )
      endif

!  Convolution of Omegamoms. Only developed so far for rayleigh

      omegamoms_c = 0.0_dpm
      if ( OMIILS_Convolution_Method ) then
         Call Convolve_OMIILS_T ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, nlayers, &
                                  ILS_Values, ILS_offsets, ILS_npoints,          &
                                  lambdas, omega(1:nlayers,:), omegamoms_c(1:nlayers,0,:) )
      else
         Call asym_gauss_f2c_T (lambdas(1:ndat), omega(1:nlayers, 1:ndat), ndat, nlayers, ndat_c, &
                                hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), omegamoms_c(1:nlayers,0,1:ndat_c) )
      endif

      if ( do_aerosols ) then
         Nmoments_LRRS = nmoms_2ns                  ! PLACEHOLDER, needs development
      else
         nmoments_LRRS = 2
         do i = 1, ndat
            beta2 = (1.0_dpm - depol(i) ) / (2.0_dpm + depol(i) )
            omegamoms_f(1:nlayers,i) = omega(1:nlayers,i) * beta2
         enddo
         if ( OMIILS_Convolution_Method ) then
            Call Convolve_OMIILS_T ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, nlayers, &
                                     ILS_Values, ILS_offsets, ILS_npoints,          &
                                     lambdas, omegamoms_f(1:nlayers,:), omegamoms_c(1:nlayers,2,:) )
         else
            Call asym_gauss_f2c_T (lambdas(1:ndat), omegamoms_f(1:nlayers, 1:ndat), ndat, nlayers, ndat_c, &
                                   hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), omegamoms_c(1:nlayers,2,1:ndat_c) )
         endif
      endif

!  Call to the Master LRRS Module
!  ==============================

      do iscene = 1, nscenes

!  write index

         Iout = Output_index(iscene)
         write(*,'(/A/)')' ** Doing LRRS-Raman calculation for '//Charscene(iscene)//' Scenario ---'

!  Local settings (scene-dependent)

         LRRS_albedo = 0.0_dpm
         if ( do_clearsky(iscene) ) then
            LRRS_nlayers = nlayers
            if ( OMIILS_Convolution_Method ) then
               Call Convolve_OMIILS_1 ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, &
                                        ILS_Values, ILS_offsets, ILS_npoints, &
                                        lambdas, albedo, LRRS_Albedo )
            else
               Call asym_gauss_f2c_1 (lambdas(1:ndat), albedo(1:ndat), ndat, ndat_c, &
                                      hwhm_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), LRRS_Albedo(1:ndat_c) )
            endif
         else
            LRRS_nlayers = ncloud_layers ; LRRS_albedo(1:ndat_c) = CldTopAlbedo
         endif

!  Call to LRRS master for each geometry, with exception handling

         do v = 1, ngeoms
            call FastV2p7_LRRS2p3_Master &
              ( E_ndat_c, E_nlayers, E_nmoms_all,                 &
                ndat_c, LRRS_nlayers, nstreams, nmoments_LRRS, EarthRadius,     &    
                sza_boa(v), vza_boa(v), azm_boa(v), heights, temperatures, Aircolumns,     &
                lambdas_c, OMI_solarflux_c, LRRS_Albedo, RayXsec_c, Depol_c, deltau_c, omegamoms_c, &
                elastic_c, raman_c, OMI_filling_c(:,v,iout), Filling_SS_c,      &
                fail_LRRS_Raman, n_LRRS_messages, LRRS_messages )
            if ( Fail_LRRS_Raman ) then
               do m = 1, n_LRRS_messages
                  b_messages(b_nmessages+m) = Trim(LRRS_messages(m))
               enddo
               b_nmessages = b_nmessages + n_LRRS_messages
               b_messages(b_nmessages+1) = 'LRRS Raman Call failed for '//Charscene(iscene)//' Scenario'
               b_nmessages = b_nmessages + 1
               b_Fail = .true. ; go to 69
            endif
         enddo

!  End scene loop

      enddo

   endif

!  timing on Raman calculation

   if (monitor_CPU) then
      call cpu_time(e3); RamanFilltime = e3-e2
   endif

!  Final Calculation of OMI Intensity/Jacobians
!  ============================================

!   ---> apply Ring Scaling and assign Ring-scaling Jacobian
!   ---> IPA calculation on normalized quantities

   if (monitor_CPU) call cpu_time(e2)

   do v = 1, ngeoms
      do w = 1, ndat_c
         do iscene = 1, nscenes
            Iout = Output_index(iscene)

!            Filling = ( 1.0_dpm + OMI_filling_c(w,v,iout) ) * OMI_solarflux_c(w)   !  Physical units
            Filling = 1.0_dpm + Ring_Scaling * OMI_filling_c(w,v,iout)
            Rad_Component(Iout)                    = Filling * IExact_C(w,v,iout)
            LPJ_Component(1:ngases,1:nlevels,Iout) = Filling * LP_LogV_Jacobians_Exact(1:ngases,1:nlevels,w,v,Iout)
            LAJ_Component(1:2,Iout)                = Filling * JExact_AlbCoeffs_C(1:2,w,v,iout)
            LRJ_Component(Iout)                    = OMI_filling_c(w,v,iout) * IExact_C(w,v,iout)

!   ---> Cloud-fraction Jacobian OMI_LFJExact, added 9/21/16

            if ( nscenes.eq.2.and.iscene.eq.nscenes ) then
               OMI_RadExact_c(w,v)                    = ClrFraction * Rad_Component(1) + &
                                                        CldFraction * Rad_Component(2)
               OMI_LPJExact_c(1:ngases,1:nlevels,w,v) = ClrFraction * LPJ_Component(1:ngases,1:nlevels,1) + &
                                                        CldFraction * LPJ_Component(1:ngases,1:nlevels,2)
               OMI_LAJExact_c(1:2,w,v)                = ClrFraction * LAJ_Component(1:2,1) + &
                                                        CldFraction * LAJ_Component(1:2,2)
               OMI_LRJExact_c(w,v)                    = ClrFraction * LRJ_Component(1) + &
                                                        CldFraction * LRJ_Component(2)
               OMI_LFJExact_c(w,v)                    = LRJ_Component(2) - LRJ_Component(1) 
            else if ( nscenes .eq. 1 ) then
               OMI_RadExact_c(w,v)                    = Rad_Component(Iout)
               OMI_LPJExact_c(1:ngases,1:nlevels,w,v) = LPJ_Component(1:ngases,1:nlevels,Iout)
               OMI_LAJExact_c(1:2,w,v)                = LAJ_Component(1:2,Iout)
               OMI_LRJExact_c(w,v)                    = LRJ_Component(Iout) 
               OMI_LFJExact_c(w,v)                    = 0.0_dpm
            endif

         enddo
!         write(999,'(i5,f12.6,1p5e17.7)')w,lambdas_c(w), OMI_solarflux_c(w), IExact_C(w,v,1:2), OMI_filling_c(w,v,1:2)
      enddo
   enddo

   if (monitor_CPU) then
      call cpu_time(e3); IPACalctime = e3-e2
      PostProcTime = Convolutiontime + RamanFilltime + IPACalctime
   endif

!  Continuation point for Shorttest option. 4/1/16.

44 continue

   write(*,*)'Finished Calculation, Cloud Fraction = ',CldFraction

!  ===================================
!  write results (for NON_HITRAN only)
!  ===================================

!  timing

   if (monitor_CPU) call cpu_time(e2)

!  Only 1 gas is output for OMI
!   -- changed, 9/21/16. All gases.

!   g = 1  !  First gas

!  OMI Results with coarse-grid calculation. 4/1/16, revised 5/6/16,
!       9/21/16: Cloud-fraction Jacobian, All trace gas Jacobians
!     EVERYTHING OUTPUT to ONE FILE

   if ( .not. do_shorttest ) then
      do v = 1, ngeoms
         open(20,file=ADJUSTL(TRIM(OMIOutput_AllResults(v))),status='replace')    ! Radiances/Jacobians
         do w = 1, ndat_c
            write(20,74)w,lambdas_c(w), OMI_solarflux_c(w), OMI_RadExact_c(w,v), &
                                        OMI_LAJexact_c(1:2,w,v), OMI_LRJExact_c(w,v), OMI_LFJExact_c(w,v)
            do n = 1, nlevels
               write(20,77) heights(n-1), (OMI_LPJexact_c(g,n,w,v), g = 1, ngases)
            enddo
         enddo
         close(20)
      enddo
   endif

!  OPTIONAL Contour output for Ozone profile weighting functions

   if ( .not. do_shorttest .and.do_LPJ_ContourPlot ) then
      g = 1  !  First gas (Ozone)
      do v = 1, ngeoms
         open(24,file=ADJUSTL(TRIM(OMIOutput_LPJContour(v))),status='replace')    ! O3 Profile Jacobians
         do w = 1, ndat_c
            do n = 1, nlevels
               write(24,78) lambdas_c(w), heights(n-1), OMI_LPJexact_c(g,n,w,v)
            enddo
            write(24,'(A)')' '
         enddo
         close(24)
      enddo
   endif

!  OPTIONAL fine-grid calculation. Only if flagged, 5/6/16. 
!     ------------- Radiances/Jacobians, all 3 scenes. Only O3 profile Jacobian

   if ( do_Output_Finegrid ) then
     do v = 1, ngeoms
       open(20,file=ADJUSTL(TRIM(fileOutput1_FineGrid(v))),status='replace')
       open(22,file=ADJUSTL(TRIM(fileOutput2_FineGrid(v))),status='replace')
       open(23,file=ADJUSTL(TRIM(fileOutput3_FineGrid(v))),status='replace')
       open(24,file=ADJUSTL(TRIM(fileOutput4_FineGrid(v))),status='replace')
       do i = 1, ndat
         write(20,57)i,lambdas(i),outgastau(i), &
                       intensity_Exact(i,v,1:3), intensity_LD_Exact(i,v,1:3),  &
                       intensity_FO_Exact(i,v,1:3), intensity_2S_Exact(i,v,1:3)
         write(22,59)i,lambdas(i),((LP_Jacobians_Exact(1,n,i,v,iscene),n=1,nlayers),iscene=1,3)
         write(23,58)i,lambdas(i),LS_Jacobians_Exact(1,i,v,1:3)
         write(24,59)i,lambdas(i),((Scatwts_Exact(n,i,v,iscene),n=1,nlayers),iscene=1,3)
       enddo
       close(22) ; close(23) ; close(24) ; close(20)
     enddo
   endif

!  Monitoring

   if (monitor_CPU) then
      call cpu_time(e3); Writetime = Writetime + e3-e2
   endif

!  Write Format statements

57 format(i5,f10.4,1pe12.4,1p12e19.10)
59 format(i5,f10.4,1p200e19.10)
58 format(i5,f10.4,1p3e19.10)

74 format(i5,f12.6,1p6e17.7)
77 format(i5,f12.6,1p4e17.7)
78 format(2f12.6,1pe17.7)

!  =============
!  FINAL SECTION
!  =============

!  Overall timing

   if ( monitor_CPU ) then
      call cpu_time(e2) ; Overalltime = e2-e1
   endif

!  write Timings. ALWAYS DO THIS

   if (monitor_CPU) then
     do v = 1, ngeoms
       call WriteTiming_Exact & 
          ( OMIOutput_Timing(v), nscenes, charscene, output_index, do_Create, &
            CreatePropTime, Exacttimes, ExactRunTime, ConvolutionTime,        &
            RamanFillTime, IPACalcTime, PostProcTime, WriteTime, OverallTime )
     enddo
   endif

!  normal Finish

   go to 70

!  Error Finish

69 continue

   write(*,*)'FastV2p7_Exact_LPS_Driver: aborted, here are messages ---'
   write(*,'(a,i3)')' * Number of error messages = ', b_nmessages
   do m = 1, b_nmessages
      write(*,'(a,i3,a,a)')' -- Message # ', m,' : ',b_messages(m)
   enddo

70 continue

!  End 

   STOP
END PROGRAM FastV2p7_Exact_LPS_Driver
