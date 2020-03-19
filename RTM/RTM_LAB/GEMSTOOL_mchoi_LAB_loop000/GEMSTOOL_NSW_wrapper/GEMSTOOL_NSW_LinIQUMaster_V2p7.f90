program GEMSTOOL_NSW_Linearized_Main_IQU

!   Target: GEMSTOOL_mchoi_t14_00
!   main test objective: convolution of optical properties for LRRS calculation (0.04nm res?)      

!   Version 1. NSW tool. Rob Spurr, Ukkyo Jeong, Yeonjin Jung.
!   12-16 August 2013.

!   Version 2A. Yonsei Team and R. Spurr, 20-24 October 2013
!   Version 2B. R. Spurr, June 2014 - Bulk Property Aerosols

!  @@ Rob fix 7/23/14, add H2O scaling flag and value

!   Version 3. Rob Spurr, 10/17/16 - 10/28/16.
!              - 10/17/16. Upgrade to VLIDORT 2.7. No real changes to this Wrapper.
!              - 10/19/16. Introduction of SIF facility.        NSW only !!!!!
!              - 10/21/16. User-defined aerosol options.        Both UVN/NSW
!              - 10/25/16. Introduction of Land BRDF facility.  Both UVN/NSW
!              - 10/25/16. Introduction of WATERLEAVING code.   Not HERE, UVN ONLY !!!
!              - 10/25/16. RTCALC_PLUS_FINAL --> separated and renamed as NSW_SetJacobians

!  Upgrade 11/30/16. Fluorescence Treatment (linear/Exact-only)

!  @@ Rob Fix 12/19/16. for Woogjung Kim
!  Separate Version Created for generating Stokes-Q and Stokes-U Jacobian Output
!    Designed to work with the Module GEMSTOOL_NSW_SetJacobians_IQU_m
!    Only need to define the new Stokes-Q and Stokes-U output arrays

!  Tool Modules
!  ------------

!  Common Type Structures

   use GEMSTOOL_Input_Types_m

!  Configuration file-read, only want the main subroutine, 10/21/16

   use GEMSTOOL_Read_Configfiles_m, Only : GEMSTOOL_Read_Configfiles

!  Initialization routine

   use GEMSTOOL_RTInitialize_m

!  Both aerosol modules have been extensively rewritten, 10/21/16.

   use GEMSTOOL_AerProperties_m
   use GEMSTOOL_AerProperties_Plus_m             ! New, June 2014

!  Clouds and albedo modules

   use GEMSTOOL_CldProperties_m
   use GEMSTOOL_Albedo_m

!  New BRDF module 10/25/16

   use GEMSTOOL_BRDFSurface_m

!  RT modules

   use GEMSTOOL_RTCALC_m
   use GEMSTOOL_RTCALC_PLUS_m

!  NSW-specific Modules. 12 August 2013
!     - New SIF module 10/18/16
!     - SetJacobians module added, 10/25/16.

   use GEMSTOOL_NSW_pthgas_profiles_m
   use GEMSTOOL_NSW_cross_sections_m
   use GEMSTOOL_NSW_cross_sections_PLUS_m
   use GEMSTOOL_NSW_regular_iops_m
   use GEMSTOOL_NSW_Fluorescence_m
   use GEMSTOOL_NSW_linearized_iops_m

!  Rob upgrade 12/19/16. Use new module for getting additional Stokes-Q and Stokes-U Jacobians

!   use GEMSTOOL_NSW_SetJacobians_m
   use GEMSTOOL_NSW_SetJacobians_IQU_m

!  add LUT                                      ! New, Feb 2015
   use GEMSTOOL_NSW_cross_sections_LUT_m
   use GEMSTOOL_NSW_cross_sections_PLUS_LUT_m

!  Use VLIDORT I/O type structures

   USE VLIDORT_PARS
   USE VLIDORT_IO_DEFS
   USE VLIDORT_LIN_IO_DEFS

   USE VLIDORT_AUX

   USE VLIDORT_INPUTS
   USE VLIDORT_MASTERS

   USE VLIDORT_L_INPUTS
   USE VLIDORT_LPS_MASTERS
   USE VLIDORT_LCS_MASTERS

!  Add the BRDF and SLEAVE type structures, 10/25/16

   USE VBRDF_SUP_MOD
   USE VBRDF_LinSup_Mod !added by mchoi

   USE VSLEAVE_SUP_MOD

! !  Add LRRS Parts
   USE FastV2p7_LRRS2p3_Master_m, Only : FastV2p7_LRRS2p3_Master
! ! Add Convolution from FastV2p7_Exact_IQU_Driver
   use AsymGauss_Convolution_Tool_m

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   implicit none

!  GEMSTOOL inputs, structure
!  --------------------------

   TYPE(GEMSTOOL_Config_Inputs) :: GEMSTOOL_INPUTS

!  VLIDORT structures
!  ------------------

!  VLIDORT input structures

   TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn
   TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn

!  VLIDORT supplements i/o structure

   TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup

!  VLIDORT output structure

   TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

!  VLIDORT linearized input structures

   TYPE(VLIDORT_Fixed_LinInputs)          :: VLIDORT_LinFixIn
   TYPE(VLIDORT_Modified_LinInputs)       :: VLIDORT_LinModIn

!  VLIDORT linearized supplements i/o structure

   TYPE(VLIDORT_LinSup_InOut)             :: VLIDORT_LinSup

!  VLIDORT linearized output structure

   TYPE(VLIDORT_LinOutputs)               :: VLIDORT_LinOut

!  VBRDF output, added 10/25/16

   TYPE(VBRDF_Sup_Outputs)                :: GEMSTOOL_BRDF_Results
   TYPE(VBRDF_LinSup_Outputs)             :: GEMSTOOL_BRDF_Results_Lin !mchoi

!  VSLEAVE output, added 10/26/16

   TYPE(VSLEAVE_Sup_Outputs)              :: GEMSTOOL_SLEAVE_Results

!  GEMSTOOL dimensioning
!  =====================

!  geophysical and IOP

!   integer, parameter :: maxwavnums   = 101

   ! integer, parameter :: maxwavnums   = 40000 ! if do_LRRS_cal is T 
   integer, parameter :: maxwavnums   = 30000
   integer, parameter :: maxgases     = 2
   integer, parameter :: maxAlbCoeffs = 5  ! This is hard-wired in the "Closure" type structure                

!  Aerosol and Cloud Moment Dimensions
!   ---> These should be the same as each other
!   ---> These should be the same as MAXMOMENTS_INPUT in Vlidort_pars.f90

   integer, parameter :: maxaermoms  = 500
   integer, parameter :: maxcldmoms  = 500

!  Messages

   integer, parameter :: maxmessages = 25

!  Weighting functions. Addition, June 2014
!   ---> Maxwfs should be the same as MAX_ATMOSWFS in Vlidort_pars.f90, should be at least Maxaerwfs + 2

   integer, parameter :: maxaerwfs  = 16
!   integer, parameter :: maxwfs     = 5     ! suggested value for profile linearizations
   integer, parameter :: maxwfs     = 30  ! maximum possible value for Column Jacobians with aerosol-bulk

!  Miscellaneous
!  -------------

!  Proxy Flux-output flag

   logical                         :: do_SphericalAlbedo

!  Proxy aerosol flag

   logical                         :: do_aerosols

!  Proxy cloud flag

   logical                         :: do_clouds

!  Proxy for gases

   integer                         :: ngases
   character*4                     :: which_gases(maxgases)
   logical                         :: do_gases(maxgases)

!  @@ Rob fix 7/23/14, Updated to include H2O scaling flag and value

   Logical                         :: do_H2OScaling
   real(kind=fpk)                  :: H2OScaling

   Logical                         :: do_CH4Scaling
   real(kind=fpk)                  :: CH4Scaling

!  Rob, 10/18/16. Local variable for SIF scaling
!   real(kind=fpk)                  :: SIFScaling

!  Jacobian control
!  ----------------

!  Surface. 10/18/16. SIFScaling Jacobian flag introduced, renamed 11/30/16

   logical                         :: do_Surface_Jacobians
!   logical                         :: do_SIFScaling_Jacobian   ! Old
   logical                         :: do_SIF_Jacobians

!  gases

   logical                         :: do_GasProfile_Jacobians
   logical                         :: do_gas_wfs(maxgases)

!  T-shift and surface-pressure and H2O scaling
!  @@ Rob fix 7/23/14, Updated to include H2O scaling Jacobian flag

   logical                         :: do_Tshift_Jacobian
   logical                         :: do_Surfpress_Jacobian
   logical                         :: do_H2OScaling_Jacobian
   logical                         :: do_CH4Scaling_Jacobian

!  Aerosol optical depth profile WF

   logical                         :: do_AerOpdepProfile_Jacobians

!  Aerosol Bulk Property WFs. New June 2014. PARSLIST added, 10/21/16
!    Bookkeeping variables
!mick fix 1/24/2017 - added Do_LinearRef, Do_LinearPSD, nlin(2), & npsd(2)

   logical                         :: do_AerBulk_Jacobians
   INTEGER                         :: n_AerBulk_Total_wfs, n_AerBulk_Group_wfs(5), AerBulk_Mask_wfs(maxaerwfs)
   REAL(fpk)                       :: AerBulk_pars(maxaerwfs)
   LOGICAL                         :: Do_LinearRef, Do_LinearPSD, Do_LinearEps
   INTEGER                         :: nlin(2), npsd(2)
   REAL(fpk)                       :: PARSLIST(10,2)

!  Flag for normalized WF output

   logical                         :: do_normalized_wfoutput

!!  Flag for doing Hitran
!
   logical                         :: do_Hitran

!  Wavenumbers
!  -----------

   real(fpk)   , dimension(maxwavnums) :: wav, wav_shift_wl,wav_shift_wn,wav_ext
   integer                             :: nwav, nwav_ext

!!!For LRRS   
!  Wavenumbers for LRRS mono options   
   logical, parameter      :: do_LRRS_cal = .false.
   INTEGER, PARAMETER :: N_RRS_N2REF = 48, &
   N_RRS_O2REF = 185, N_RRS_ALL = 234, nwav_shift=234

   ! real(fpk)   , dimension(maxwavnums,N_RRS_ALL) :: wav_lrrs_WN, wav_lrrs_WL
   ! integer                             :: nwav_lrrs  
   !  Excitation wavelength
   REAL(FPK) :: LAMBDA_EXCIT
   Real(fpk), dimension(N_RRS_ALL) :: LAMBDAS_RANKED, WAVENUM_RANKED
   INTEGER ::   NPOINTS_MONO, W_EXCIT_WN, W_EXCIT_WL, N_LOUTPUT_read
   
   !  Positions
   REAL(FPK) :: O2POS (N_RRS_O2REF)
   REAL(FPK) :: N2POS (N_RRS_N2REF)
   REAL(FPK) :: SIGMA, SIGMA_0, LAMBDAS_ALL(N_RRS_ALL), WAVENUM_ALL (N_RRS_ALL)
   INTEGER ::   J, S, NSHIFTS, INDEX_ALL ( N_RRS_ALL ), WW, ia,ib
   integer      :: i_scene_lrrs, n_scenes_lrrs

!  SolarFlux for LRRS. dataset
   Integer, Parameter :: N_solar_read = 800001 !12000-15000 0.01 cm-1 resolution. ;-- GeffToon!!
   ! Integer, Parameter :: N_solar_read = 300001 !12000-15000 0.01 cm-1 resolution. ;-- GeffToon!!
   ! Integer, Parameter :: N_solar_read = 80093 !12000-15000 0.01 cm-1 resolution. ;-- SAO!!


   REAL(FPK), dimension(N_solar_read):: WN_solar_read, WL_solar_read, &
         Solar_read_nm, Solar_read_cm, Solar_read_photon


!  Interpolate - same var name with 'GEMSTOOL_AerProperties.f90 line 1403 ;mchoi
   Integer :: wastart, wa, wa1, wa2, wb, wc
   Real(FPK) :: fa1, fa2
   REAL(FPK), dimension(nwav_shift) :: Solar_photon_shift
   logical :: trawl

!  for convolution 
   logical :: do_gauss_conv_before_LRRS_calculation


!!!Endo of For LRRS


!  Surface pressure

!   real(fpk)                    :: surface_pressure

!  Critical size for aerosol moment

   real(fpk)   , parameter      :: momsize_cutoff = 1.0d-03

!  critical SSA limit

   real(fpk)   , PARAMETER      :: OMEGA_LIM = 0.999999d0

!  CO2 mixing ratio. @@@@@@ Updated to 395, 2/28/13

!   real(fpk)   , PARAMETER      :: CO2_PPMV_MIXRATIO = 385.0d0
   real(fpk)   , PARAMETER      :: CO2_PPMV_MIXRATIO = 395.0d0

!  File path names
!   Rob Fix, 10/18/16 SIF Datapath added.

   character*256 :: SIF_Datapath
   character*256 :: PhysicsDataPath
   character*256 :: PTH_profile_name
   character*256 :: GAS_profile_names(maxgases)
   character*256 :: Hitran_path
   character*256 :: Lut_path
   character*256 :: ConfigPath
   character*256 :: alb_filename

!  Number of stokes components (Proxy) 
!  Number of Muller matrix coefficients ( = 1, for NSTOKES = 1, = 6 for NSTOKES = 4 )
!  Number of Geometries, solar zenith angles (Proxies)

   integer             :: nstokes, nmuller, n_geometries, n_szas, ii, i_layer, i_moment


!  ATMOSPHERIC PROPERTY variables
!  ==============================

!  Trace gas PTH and Layer column amounts
!  -----------------------------------------

    !   Units should be consistent, e.g:
    !   gas columns are in mol/cm^2  or inverse [DU]
!  @@ Rob fix 7/23/14, add H2O scaling derivative

    integer  :: nlayers, nlevels

    real(fpk)   , DIMENSION( 0:maxlayers )               :: LEVEL_TEMPS
    real(fpk)   , DIMENSION( 0:maxlayers )               :: LEVEL_HEIGHTS
    real(fpk)   , DIMENSION( 0:maxlayers )               :: LEVEL_PRESS
    REAL(fpk)   , DIMENSION( 0:maxlayers, maxgases )     :: LEVEL_VMRS

    REAL(fpk)   , DIMENSION ( 0:maxlayers )              :: TSHAPE

    REAL(fpk)   , DIMENSION( maxlayers )                 :: AIRCOLUMNS
    REAL(fpk)   , DIMENSION( maxlayers )                 :: dAIRCOLUMNS_dS
    REAL(fpk)                                            :: dAIRCOLUMN_dP

    REAL(fpk)   , DIMENSION( maxlayers, maxgases, 2 )    :: GASCOLUMNS
    REAL(fpk)   , DIMENSION( maxlayers, maxgases, 2 )    :: dGASCOLUMNS_dS
    REAL(fpk)   , DIMENSION( maxlayers, maxgases, 2 )    :: dGASCOLUMNS_dv
    REAL(fpk)   , DIMENSION( maxgases )                  :: dGASCOLUMN_dP

!  @@ Rob fix 7/23/14, add H2O scaling derivative

    REAL(fpk)   , DIMENSION( maxlayers, 2 )              :: dH2OCOLUMN_dF
    REAL(fpk)   , DIMENSION( maxlayers, 2 )              :: dCH4COLUMN_dF

!  Cross-sections
!  --------------

    !   Units should be consistent, e.g:
    !   X-sections  are in cm^2/mol  or [DU]

!  Rayleigh cross-sections and depolarization output

      real(fpk),    dimension ( maxwavnums ) :: RAYLEIGH_XSEC
      real(fpk),    dimension ( maxwavnums ) :: RAYLEIGH_DEPOL

!  Gas cross-sections, all levels, all gases. Also, P/T derivatives.

      REAL(fpk),    DIMENSION( maxwavnums, 0:maxlayers, maxgases ) :: GAS_XSECS
      REAL(fpk),    DIMENSION( maxwavnums, 0:maxlayers, maxgases ) :: dGASXSECS_dT
      REAL(fpk),    DIMENSION( maxwavnums, 0:maxlayers, maxgases ) :: dGASXSECS_dP

!  Aerosol/Cloud optical properties
!  --------------------------------

!  aerosol Loading output

    real(fpk)   , dimension ( maxlayers )      :: Loading_aerosol
    real(fpk)   , dimension ( maxlayers )      :: DLoading_Dtau
    real(fpk)   , dimension ( maxlayers, 2 )   :: Dloading_Dpars

!   Aerosol optical properties (AOPS)
!     10/19/16.  aerosol_nscatmoms now with wavenumber dependence

    INTEGER     , DIMENSION( maxwavnums )                   :: aerosol_nscatmoms
    LOGICAL,      DIMENSION ( maxlayers )                   :: AEROSOL_LAYERFLAGS 
    real(fpk)   , DIMENSION ( maxlayers, maxwavnums )       :: AEROSOL_DELTAU

    real(fpk)   , DIMENSION ( maxlayers )                   :: aertau_unscaled
    real(fpk)   , DIMENSION ( maxwavnums )                  :: AOD_SCALING
    real(fpk)   , DIMENSION ( maxwavnums )                  :: AEROSOL_SSALBS 
    real(fpk)   , DIMENSION ( 6, 0:maxaermoms, maxwavnums ) :: AEROSOL_SCATMOMS

!  Aerosol linearized optical properties

    real(fpk)   , DIMENSION ( maxlayers, maxaerwfs )                   :: L_aertau_unscaled
    real(fpk)   , DIMENSION ( maxlayers, maxwavnums, maxaerwfs )       :: L_AEROSOL_DELTAU
    real(fpk)   , DIMENSION ( maxwavnums, maxaerwfs )                  :: L_AOD_SCALING
    real(fpk)   , DIMENSION ( maxwavnums, maxaerwfs )                  :: L_AEROSOL_SSALBS 
    real(fpk)   , DIMENSION ( 6, 0:maxaermoms, maxwavnums, maxaerwfs ) :: L_AEROSOL_SCATMOMS

!  Aerosol distribution characterisstics.

    real(fpk),     DIMENSION ( 5, 2 )             :: AEROSOL_DISTCHARS

!  Cloud loading

    real(fpk)   , dimension ( maxlayers)                   :: loading_cloud

!  Cloud optical properties

    INTEGER                                                :: cloud_nscatmoms
    LOGICAL,      DIMENSION( maxlayers )                   :: CLOUD_LAYERFLAGS
    real(fpk)   , DIMENSION( maxlayers, maxwavnums )       :: CLOUD_DELTAU

    real(fpk)   , DIMENSION ( maxlayers )                  :: cldtau_unscaled
    real(fpk)   , DIMENSION( maxwavnums )                  :: COD_SCALING
    real(fpk)   , DIMENSION( maxwavnums )                  :: CLOUD_SSALBS 
    real(fpk)   , DIMENSION( 6, 0:maxcldmoms, maxwavnums ) :: CLOUD_SCATMOMS 

!  Surface properties
!  ------------------

    real(fpk)    :: LAMBERTIAN_ALBEDO (maxwavnums)
    real(fpk)    :: ClosureTerms(maxwavnums,MaxAlbCoeffs)



!  Rob 10/18/16. Added GEMSTOOL_SIF arrays
!   real(fpk) :: Gemstool_SIF          (Maxwavnums,Max_Geometries)
!   real(fpk) :: dGemstoolSIF_dScaling (Maxwavnums,Max_Geometries)

!  VLIDORT PROXIES
!  ===============

!  Regular optical properties

    INTEGER      :: NGREEK_MOMENTS_INPUT
    real(fpk)    :: DELTAU_VERT_INPUT     ( MAXLAYERS )
    real(fpk)    :: OMEGA_TOTAL_INPUT     ( MAXLAYERS )
    real(fpk)    :: GREEKMAT_TOTAL_INPUT  ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Linearized option properties

    REAL(fpk)    :: L_DELTAU_VERT_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
    REAL(fpk)    :: L_OMEGA_TOTAL_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
    REAL(fpk)    :: L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Extended for LLR calculation: include wavnumber dimension
    real(fpk)    :: DELTAU_VERT_INPUT_ext     ( MAXLAYERS, maxwavnums)
    real(fpk)    :: OMEGA_TOTAL_INPUT_ext     ( MAXLAYERS, maxwavnums)
    real(fpk)    :: GREEKMAT_TOTAL_INPUT_ext  ( 0:MAXMOMENTS_INPUT, MAXLAYERS, maxwavnums) ! only for Intensity
   !  real(fpk)    :: temp_moments_ext     ((1+MAXMOMENTS_INPUT)*MAXLAYERS, maxwavnums)
    real(fpk)    :: temp_moments_ext     ( 0:MAXMOMENTS_INPUT, maxwavnums)
   !  real(fpk)    :: RAYLEIGH_XSEC_ext  (maxwavnums)
   !  real(fpk)    :: RAYLEIGH_DEPOL_ext (maxwavnums)

!  shifted for LLR calculation: include wavnumber dimension (nwav_shift= 234)
    real(fpk)    :: DELTAU_VERT_INPUT_shift     ( MAXLAYERS, nwav_shift)
    real(fpk)    :: OMEGA_TOTAL_INPUT_shift     ( MAXLAYERS, nwav_shift)
    real(fpk)    :: GREEKMAT_TOTAL_INPUT_shift  ( MAXLAYERS, 0:MAXMOMENTS_INPUT, nwav_shift) ! only for Intensity
    real(fpk)    :: RAYLEIGH_XSEC_shift  (nwav_shift)
    real(fpk)    :: RAYLEIGH_DEPOL_shift (nwav_shift)
    real(fpk)    :: Lambertian_albedo_shift (nwav_shift)   

!  LRRS Variables. Not using the Linearized LRRS code (YET!!!). 2/4/16

    real(fpk) :: EarthRadius
    real(fpk) :: sza_boa(1), vza_boa(1), azm_boa(1)
 
    ! integer        :: LRRS_nlayers, Nmoments_LRRS
    ! real(kind=dpm) :: LRRS_Albedo   (E_ndat_c)
    real(kind=fpk) :: wn_temp    (maxwavnums)
    real(kind=fpk) :: fluxes_c    (maxwavnums)
    real(kind=fpk) :: lambda_c    (maxwavnums)
 
    real(kind=fpk) :: Elastic_c    (nwav_shift) !(maxwavnums)
    real(kind=fpk) :: Raman_c      (nwav_shift) !(maxwavnums)
    real(kind=fpk) :: Filling_SS_c (nwav_shift) !(maxwavnums)
    real(kind=fpk) :: OMI_filling_c   (nwav_shift) !(maxwavnums)
 
    real(kind=fpk) :: Solar_flux_all (maxwavnums)
    real(kind=fpk) :: Elastic_all    (maxwavnums)
    real(kind=fpk) :: Raman_all      (maxwavnums)
    real(kind=fpk) :: Filling_SS_all (maxwavnums)
    real(kind=fpk) :: OMI_filling_all   (maxwavnums)

    ! using convolution!
    
    real(kind=fpk), DIMENSION (maxwavnums) :: FWHM_c, hw1e_c, asym_c
    real(kind=fpk), DIMENSION (nwav_shift) :: Solar_flux_c, RAYLEIGH_XSEC_c, RAYLEIGH_DEPOL_c
    real(kind=fpk), DIMENSION (nwav_shift) :: Lambertian_albedo_c
    real(fpk)    :: DELTAU_VERT_INPUT_ext_c     ( MAXLAYERS, nwav_shift)
    real(fpk)    :: OMEGA_TOTAL_INPUT_ext_c     ( MAXLAYERS, nwav_shift)
    real(fpk)    :: GREEKMAT_TOTAL_INPUT_ext_c  ( 0:MAXMOMENTS_INPUT, MAXLAYERS, nwav_shift) 
    real(fpk)    :: temp_moments_ext_c     ( 0:MAXMOMENTS_INPUT, nwav_shift)
   !  real(fpk)    :: temp_moments_ext_c     ((1+MAXMOMENTS_INPUT)*MAXLAYERS, nwav_shift)
    integer :: ndat_c, nspec
 
    logical        :: Fail_LRRS_Raman
    integer        :: n_LRRS_messages
    character*100  :: LRRS_Messages(100)
    integer, parameter :: b_max_messages   = 25
    logical       :: b_fail
    integer       :: b_nmessages
    character*100 :: b_messages (b_max_messages)
 
    ! for reading flux files
    character*100 temp0, temp1
    character*300 temp2,file
 
    !  For LRRS preparation
    real(fpk)    :: DELTAU_INPUT_UNSCALED  ( MAXLAYERS, maxwavnums )
    real(fpk)    :: OMEGAMOMS_ELASTIC_UNSCALED  ( MAXLAYERS, 0:MAXMOMENTS_INPUT, maxwavnums )
    

    ! For BRDF Jacobian mchoi
    integer :: BRDF_N_SURFACE_WFS, N_SURFACE_WFS_Fin
    real(fpk), dimension(16) :: BRDF_K_FP_values !4+12
    logical :: do_BRDF

!  Additional output (for Layers-to-Levels Transformation)
!  =================

!  Layer gas optical depths

    REAL(fpk)    :: LAYER_GASOPDEPS       ( maxlayers, maxgases )

!  Chain-Rule Transformation from Layer-GASOPDEP to Level-VMR Jacobians

    REAL(fpk)    :: dGASOPDEPS_dV       ( maxlayers, maxgases, 2 )

!  Output arrays
!  =============

!  Overall amounts
!  ---------------

!  Radiances, polarization, fluxes

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAVNUMS)                     :: STOKES
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,  MAXWAVNUMS)                     :: DOLP
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAVNUMS)                     :: ACTINIC
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAVNUMS)                     :: REGFLUX

!  Jacobians
!  @@ Rob    fix 7/23/14 , add H2O scaling Jacobians (RADIANCE_WSCALEWFS)
!  @@ Y.Jung fix 2/1/15  , add CH4 scaling Jacobians (RADIANCE_MSCALEWFS)

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,0:MAXLAYERS,MAXGASES)  :: RADIANCE_GASWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXLAYERS)             :: RADIANCE_AERPROFWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXAERWFS)             :: RADIANCE_AERBULKWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MaxAlbCoeffs)          :: RADIANCE_ALBEDOWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: RADIANCE_TSHIFTWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: RADIANCE_SURFPWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: RADIANCE_WSCALEWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: RADIANCE_MSCALEWFS

!  Derivatives of Stokes Radiance quantities w.r.t SIF parameters. New, 10/18/16, 11/30/16
!       Gaussian = Amplitude-Scaling, Linear = Scaling Constant and Scaling Gradient

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,2)                     :: RADIANCE_SIFPARAMWFS

!  Stokes-Q and Stokes-U Jacobians
!  @@ Rob added 12/19/16

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,0:MAXLAYERS,MAXGASES)  :: STOKES_Q_GASWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXLAYERS)             :: STOKES_Q_AERPROFWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXAERWFS)             :: STOKES_Q_AERBULKWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MaxAlbCoeffs)          :: STOKES_Q_ALBEDOWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: STOKES_Q_TSHIFTWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: STOKES_Q_SURFPWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: STOKES_Q_WSCALEWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: STOKES_Q_MSCALEWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,2)                     :: STOKES_Q_SIFPARAMWFS

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,0:MAXLAYERS,MAXGASES)  :: STOKES_U_GASWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXLAYERS)             :: STOKES_U_AERPROFWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXAERWFS)             :: STOKES_U_AERBULKWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MaxAlbCoeffs)          :: STOKES_U_ALBEDOWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: STOKES_U_TSHIFTWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: STOKES_U_SURFPWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: STOKES_U_WSCALEWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                       :: STOKES_U_MSCALEWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,2)                     :: STOKES_U_SIFPARAMWFS

! BRDF Jacobian ; mchoi


   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,12) :: RADIANCE_BRDFWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,12) :: STOKES_Q_BRDFWFS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,12) :: STOKES_U_BRDFWFS



!  Air mass factors

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXLAYERS,MAXGASES)  :: AMFSCATWTS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXGASES)            :: AMFTOTAL

!  Scene quantities
!  ----------------

!  Radiances (I) + other stokes components (Q,U,V)

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAVNUMS)                    :: STOKES_ic

!  Degree of linear/circular polarization (DOLP/DOCP)

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                      :: DOLP_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)                      :: DOCP_ic

!  Actinic and Regular Flux output

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAVNUMS)                     :: ACTINIC_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAVNUMS)                     :: REGFLUX_ic

!  Derivatives of Stokes Radiance quantities w.r.t GAS VMRS
!           [ Gas weighting function w.r.t. LEVEL VMRS ] 

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,0:MAXLAYERS,MAXGASES)  :: RADIANCE_GASWFS_ic

!  AMF scattering weights, total AMF

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXLAYERS,MAXGASES)  :: AMFSCATWTS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXGASES)            :: AMFTOTAL_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXGASES)            :: TOTALWF_SAVED_ic

!  Derivatives of Stokes Radiance quantities w.r.t AEROSOL OPTICAL DEPTH PROFILE AT REFERENCE WAVELENGTH
!           [ weighting function w.r.t. LAYER optical depth amounts ] 

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXLAYERS)  :: RADIANCE_AERPROFWFS_ic

!  Derivatives of Stokes Radiance quantities w.r.t AEROSOL BULK PROPERTIES

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXAERWFS)  :: RADIANCE_AERBULKWFS_ic

!  Derivatives of Stokes Radiance quantities w.r.t Surface albedo Coefficients

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MaxAlbCoeffs)  :: RADIANCE_ALBEDOWFS_ic

!  Derivatives of Stokes Radiance quantities w.r.t SIF parameters. New, 10/18/16, 11/30/16
!       Gaussian = Amplitude-Scaling, Linear = Scaling Constant and Scaling Gradient

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,2)  :: RADIANCE_SIFPARAMWFS_ic

!  Derivatives of Stokes Radiance quantities w.r.t. T-Shift

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)  :: RADIANCE_TSHIFTWFS_ic

!  Derivatives of Stokes Radiance quantities w.r.t. Surface pressure

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)  :: RADIANCE_SURFPWFS_ic

!  Derivatives of Stokes Radiance quantities w.r.t. H2O Scaling

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)  :: RADIANCE_WSCALEWFS_ic

!  Derivatives of Stokes Radiance quantities w.r.t. CH4 Scaling

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)  :: RADIANCE_MSCALEWFS_ic

!  STOKES-Q and STOKES-U Jacobian output arrays
!  --------------------------------------------

!  Introduced by R. Spurr, 12/19/16. Names are same as above

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,0:MAXLAYERS,MAXGASES) :: STOKES_Q_GASWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXLAYERS)            :: STOKES_Q_AERPROFWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXAERWFS)            :: STOKES_Q_AERBULKWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MaxAlbCoeffs)         :: STOKES_Q_ALBEDOWFS_ic

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,2) :: STOKES_Q_SIFPARAMWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)   :: STOKES_Q_TSHIFTWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)   :: STOKES_Q_SURFPWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)   :: STOKES_Q_WSCALEWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)   :: STOKES_Q_MSCALEWFS_ic

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,0:MAXLAYERS,MAXGASES) :: STOKES_U_GASWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXLAYERS)            :: STOKES_U_AERPROFWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MAXAERWFS)            :: STOKES_U_AERBULKWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,MaxAlbCoeffs)         :: STOKES_U_ALBEDOWFS_ic

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,2) :: STOKES_U_SIFPARAMWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)   :: STOKES_U_TSHIFTWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)   :: STOKES_U_SURFPWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)   :: STOKES_U_WSCALEWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS)   :: STOKES_U_MSCALEWFS_ic

   !BRDF jacobian by mchoi
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,12) :: RADIANCE_BRDFWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,12) :: STOKES_Q_BRDFWFS_ic
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAVNUMS,12) :: STOKES_U_BRDFWFS_ic

!  Exception handling
!  ==================

!  Highest level

   integer             ::  nmessages
   LOGICAL             ::  FAIL_Main
   CHARACTER(LEN=256)  ::  messages_Main(maxmessages)

!  Subroutine level.
!     * Bookkeep stuff added 10/21/16, BRDF variables introduced 10/25/16

   integer             ::  nmessages_BRDF
   LOGICAL             ::  Fail1_Cld, Fail1_Aer, Fail2_Aer, Fail_cfg, Fail_pthg, &
                           Fail_xsec, Fail_RTinit, Fail_bookkeep, Fail_BRDF
   CHARACTER(LEN=100)  ::  message_cfg,  message_pthg, message_Xsec, message_RTinit, &
                           message_Bookkeep(3),  messages_BRDF(3),                   &
                           message_AerOptical(5), message_AerLoading(3) , message_CldOptical(4)

!  local variables
!  ===============

!  SIF variables. Added 11/30/16

   Logical   :: SIF_ExactOnly
   Integer   :: n_SIFPars
   real(fpk) :: Iso_Value, Polynomial, dPolynomial_dSL(2), factor, SIFPars(2)

!  Control for using wavenumbers

   logical                 :: do_wavnums

!  Temperature shift (shoudl be part of the configuration)

   ! real(fpk), parameter    :: Tshift = 20.0d0
   real(fpk)               :: Tshift

!  Surface press

   real(fpk)               :: Surfpress

!  Interpolate aerosols and clouds flags

!   logical, parameter      :: interpolate_aerosols = .false.
   logical, parameter      :: interpolate_aerosols = .true.

!   logical, parameter      :: interpolate_clouds   = .false.
   logical, parameter      :: interpolate_clouds   = .true.

!  Cloud flag (temporary)

!mick fix 6/18/2013 - commented out now
   logical, parameter      :: do_CLOUDS_V2 = .true.

!  Local cloud scene

   real(fpk)    :: cfrac, scene_weight(2)
   logical      :: local_cloud(2)
   integer      :: i_scene, n_scenes
   character*1  :: c1

!  Miscellaneous

   real(fpk)    :: range, step, interval
   integer      :: m, n, w, k, q, qs, g, gh2o, gch4, Errorstatus, UpDn_Index, dir, ll, o
   character*5  :: COUTPUT(2), CLEVEL
   data COUTPUT / 'TOAUP','BOADN' /

!  Indices for Finite Difference testing, Perturbation FDeps

!    FDAer Must be True of false
!    FDGas Must be 1, 2, 3.... up to ngases
!    FDLev Must be 0, 1, 2, 3 ... up to nlayers
!    FDLay Must be 1, 2, 3 ... up to nlayers

   logical   :: FDAer, FDsfp
   integer   :: FDGas, FDLev, FDLay, FDBul
   real(fpk) :: FDeps

!  Debug variables (you can comment these out when finished with them)

   logical, parameter :: XSEC_debug_flag = .true.
   ! logical, parameter :: IOP_debug_flag = .true.
  logical, parameter :: IOP_debug_flag = .false.
   character*60        :: IOP_debug_filename

!  Start program
  !@@@@@@@@@@@@@

!  Initialise exception handling

   messages_Main  = ' '
   nmessages = 0
   fail_Main = .false.
   IOP_debug_filename = '14August13_IOP_debug.out'

!  the FD test values 
!  ------------------

   FDeps = 0.0d0          ! Should be zero for baseline run

!  Gas profile FD testing

   FDGas = -9999
   FDLev = -9999

!  Aerosol profile WF testing

   FDAer = .false.      ! Y.Jung
   !FDAer = .true.
   !FDLay = 18

!  Aerosol Bulk WF testing

   FDBul = 0

!  Surface pressure

   FDSfp = .false.













   
!  Configuration file setup
!  ========================

!  set the Configuration path,  CALL THE READ_CONFIGURATION FILES
!  12 August 2013. Use NSW-specific path here

   do_wavnums = .true.
   ConfigPath = 'GEMSTOOL_NSW_Configfiles/'
   CALL GEMSTOOL_Read_Configfiles &
    ( do_wavnums, ConfigPath, GEMSTOOL_INPUTS, fail_cfg, message_cfg  )

!  Exception handling

   if ( fail_cfg ) then
      fail_Main = .true.
      nmessages = nmessages + 1
      messages_Main(nmessages) = adjustl(trim(message_cfg))
      go to 4555
   endif

!  No water-leaving (10/25/16)

   if ( GEMSTOOL_INPUTS%Sleave%DO_WaterLeaving ) then
      fail_Main = .true.
      messages_main(nmessages+1) = 'This is the NSW tool: WATER-LEAVING not allowed, turn off Flag'
      nmessages = nmessages + 1
      go to 4555
   endif

!  Proxies (flags and Numbers)

   UpDn_index = 2 ; if ( GEMSTOOL_INPUTS%RTControl%do_upwelling ) UpDn_index = 1
   do_SphericalAlbedo  = GEMSTOOL_INPUTS%RTControl%do_SphericalAlbedo

   do_aerosols  = GEMSTOOL_INPUTS%Atmosph%DO_AEROSOLS
   do_clouds    = GEMSTOOL_INPUTS%Atmosph%DO_CLOUDS
   nstokes      = GEMSTOOL_INPUTS%RTControl%NVlidort_nstokes
   nmuller = 1 ; if ( nstokes.gt.1 ) nmuller = 6

   n_szas       = GEMSTOOL_INPUTS%Geometry%n_GEMS_geometries        ! Only for Observational Geometry Mode.
   n_geometries = GEMSTOOL_INPUTS%Geometry%n_GEMS_geometries

!  Copy gas information to proxies

   ngases = GEMSTOOL_INPUTS%TraceGas%ngases
   which_gases = '    ' ; do_gases = .false.
   which_gases(1:ngases) = GEMSTOOL_INPUTS%TraceGas%which_gases(1:ngases)
   do_gases(1:ngases)    = GEMSTOOL_INPUTS%TraceGas%do_gases(1:ngases)

!  @@ Rob fix 7/23/14, add H2O scaling flag and value

   do_H2OScaling =  GEMSTOOL_INPUTS%Tracegas%do_H2OScaling
   H2OScaling    =  GEMSTOOL_INPUTS%Tracegas%H2OScaling
   print*, 'H2OScaling = ', H2OScaling

   gh2o = 0
   if ( do_H2OScaling ) then
      do g = 1, ngases
         if ( which_gases(g).eq.'H2O ' ) gh2o = g
      enddo
   endif

!  @@ Y.Jung fix 1/2/15, add CH4 scaling flag and value

   do_CH4Scaling =  GEMSTOOL_INPUTS%Tracegas%do_CH4Scaling
   CH4Scaling    =  GEMSTOOL_INPUTS%Tracegas%CH4Scaling
   print*, 'CH4Scaling = ', CH4Scaling

   gch4 = 0
   if ( do_CH4Scaling ) then
      do g = 1, ngases
         if ( which_gases(g).eq.'CH4 ' ) gch4 = g
      enddo
   endif

! @@ Rob Fix 10/18/16. Proxy for the SifScaling value. REMOVED 11/30/16.
!   SIFScaling = 0.0d0 ; if ( GEMSTOOL_INPUTS%Sleave%do_SIF ) SIFScaling =  GEMSTOOL_INPUTS%Sleave%SIF755_Amplitude

!  Proxies for Linearization control ( Copy type structure variables to local variables )
!     add SIFScaling 10/18/16, renamed 11/30/16

   do_GasProfile_Jacobians        = GEMSTOOL_INPUTS%LinControl%do_GasProfile_Jacobians
   do_AerOpdepProfile_Jacobians   = GEMSTOOL_INPUTS%LinControl%do_AerOpdepProfile_Jacobians
   do_AerBulk_Jacobians           = GEMSTOOL_INPUTS%LinControl%do_AerBulk_Jacobians
   do_Surface_Jacobians           = GEMSTOOL_INPUTS%LinControl%do_Surface_Jacobians
   do_Tshift_Jacobian             = GEMSTOOL_INPUTS%LinControl%do_Tshift_Jacobian
   do_SurfPress_Jacobian          = GEMSTOOL_INPUTS%LinControl%do_SurfPress_Jacobian
   do_H2OScaling_Jacobian         = GEMSTOOL_INPUTS%LinControl%do_H2OScaling_Jacobian
   do_CH4Scaling_Jacobian         = GEMSTOOL_INPUTS%LinControl%do_CH4Scaling_Jacobian
   do_SIF_Jacobians               = GEMSTOOL_INPUTS%LinControl%do_SIF_Jacobians
   do_normalized_wfoutput         = GEMSTOOL_INPUTS%LinControl%do_normalized_wfoutput
   do_Hitran                      = GEMSTOOL_INPUTS%LinControl%do_Hitran


   Tshift = GEMSTOOL_INPUTS%LinControl%Tshift_read
   ! write(*,*) 'Tshift = ', Tshift
!  Set local gas wfs flags

   if ( do_GasProfile_Jacobians ) then
      do_gas_wfs(1:ngases) = GEMSTOOL_INPUTS%TraceGas%do_gas_wfs(1:ngases)
   else
      do_gas_wfs = .false.
   endif

!  Proxy names

   GAS_profile_names = ' '
   GAS_profile_names(1:ngases) = GEMSTOOL_INPUTS%TraceGas%GAS_profile_names(1:ngases) 

!  PTH and GAS profile, from Initial profile 12 august 2013.

   PTH_profile_name            = 'atmos_mchoi.dat'    !Y.Jung
   GAS_profile_names(1:ngases) = 'atmos_mchoi.dat'    !Y.Jung
   ! PTH_profile_name            = 'atmos_mchoi_t12.dat'    !Y.Jung
   ! GAS_profile_names(1:ngases) = 'atmos_mchoi_t12.dat'    !Y.Jung
   ! alb_filename                = 'alb.dat'      !Y.Jung
!  PTH_profile_name            = 'Initial_NSW_atmosphere_12aug13.dat'
!   GAS_profile_names(1:ngases) = 'Initial_NSW_atmosphere_12aug13.dat'

  !  PTH_profile_name            = 'Profile_PTH_US76_GAS_discoveraq_fin_mchoi_alt0.dat'    !mchoi
  !  GAS_profile_names(1:ngases) = 'Profile_PTH_US76_GAS_discoveraq_fin_mchoi_alt0.dat'    !mchoi
  !  alb_filename                = 'alb.dat'      !???


!  Data paths

   ! PhysicsDataPath = '../GEMSTOOL_physicsdata/'
   ! Hitran_path     = trim(PhysicsDataPath)//'HITRAN_Xsections_data/HITRAN_2016_all_mchoi/'
   PhysicsDataPath = '../GEMSTOOL_physicsdata/' !'../../GEMSTOOL_physicsdata/'
   Hitran_path     = '../../../HITRAN_Xsections_data_for_all/HITRAN_2016_all_mchoi/'
   Lut_path     = trim(PhysicsDataPath)//'LUT_Xsections_data/results_files_V500/'
   ! Lut_path     = trim(PhysicsDataPath)//'LUT_Xsections_data/results_files_V500/'
   ! Lut_path     = trim(PhysicsDataPath)//'LUT_Xsections_data/results_files_V420/'
   ! Lut_path     = trim(PhysicsDataPath)//'LUT_Xsections_data/results_files/'


!  Add SIF data path 10/18/16

   SIF_datapath = trim(PhysicsDataPath)//'Fluorescence_data/'

!  Checking
!  --------

   if ( .not. do_aerosols .and. do_AerOpdepProfile_Jacobians ) then
      fail_Main = .true.
      nmessages = nmessages + 1
      messages_Main(nmessages) = 'No aerosols!. Cannot have aerosol profile Jacobians - check input flags!!!'
      go to 4555
   endif

   if ( .not. do_aerosols .and. do_AerBulk_Jacobians ) then
      fail_Main = .true.
      nmessages = nmessages + 1
      messages_Main(nmessages) = 'No aerosols!. Cannot have aerosol Bulk-property Jacobians - check input flags!!!'
      go to 4555
   endif

!   PUT IN SOME MORE CHECKS PLEASE .......................!!!!!

!  Initialize the RT models and wavelengths
!  ========================================

!  (VLIDORT, FO), based on these Configuration file settings

   call GEMSTOOL_VLIDORT_Initialize &
   ( GEMSTOOL_INPUTS, VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, fail_RTinit, message_RTinit )

   if ( fail_RTinit ) then
      fail_Main = .true.
      nmessages = nmessages + 1
      messages_Main(nmessages) = adjustl(trim(message_RTinit))
      go to 4555
   endif

!  setting the Linearized VLIDORT control input here
!  -------------------------------------------------

!  Initialize - not used

!  10/17/16. These 2 variables no longer in Version 2.7
!   VLIDORT_LinFixIn%Cont%TS_DO_SURFBB_LINEARIZATION = .false.
!   VLIDORT_LinFixIn%Cont%TS_DO_LTE_LINEARIZATION    = .false.

   VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES     = ' '
   VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES    = ' '

   VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS     = .false.
   VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS       = 0

!  Initialize - main control variables

!   VLIDORT_LinFixIn%Cont%TS_DO_SIMULATION_ONLY        =  .false. !  Moved 10/17/16
   VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY        =  .false.
   VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  =  .false.
   VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION =  .false.
   VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   =  .false.
   VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION =  .false.
   VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         =  .false.

! 10/17/16. 2 new variables in V2p7, not required here....

   VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF   =  .false.
   VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF =  .false.

   VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS         = 0
   VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = 0
   VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 0

   VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG   = .false.
   VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER = 0

!  Initialize - supplemental BRDF variables

   VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC = ZERO
   VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0        = ZERO
   VLIDORT_LinSup%BRDF%TS_LS_BRDF_F          = ZERO
   VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0   = ZERO
   VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F     = ZERO

   VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY      = ZERO
   VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY = ZERO

!mick fix 1/24/2017 - added linearized SLEAVE input initializations
!  Initialize - supplemental SLEAVE variables

   VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC  = ZERO
   VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES = ZERO
   VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0        = ZERO
   VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0   = ZERO

!  Profile Jacobian flag

   if ( do_GasProfile_Jacobians .or. do_AerOpdepProfile_Jacobians ) then
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION =  .true.
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   =  .true.
   endif

!  COLUMN Jacobians, 15 august 2013. Tshift and/or Surface-pressure
!  @@ Rob fix 7/23/14, Update for H2O scaling flag and value

   Q = 0
!   if ( do_Tshift_Jacobian )     Q = Q + 1
!   if ( do_Surfpress_Jacobian )  Q = Q + 1
!   if ( do_H2OScaling_Jacobian ) Q = Q + 1

   if ( do_H2OScaling_Jacobian ) Q = Q + 1
   if ( do_CH4Scaling_Jacobian ) Q = Q + 1
   if ( do_Tshift_Jacobian )     Q = Q + 1
   if ( do_Surfpress_Jacobian )  Q = Q + 1

!  BULK (COLUMN) AEROSOL Jacobians, New section June 2014
!    10/21/16. PARSLIST Added, Exception handling introduced.
!mick fix 1/24/2017 - added Do_LinearRef, Do_LinearPSD, Do_LinearEps, nlin(2), & npsd(2) to output

   if ( do_AerBulk_Jacobians ) then
      call GEMSTOOL_AERWFS_BOOKKEEP &
         ( maxaerwfs, GEMSTOOL_INPUTS,                & ! Inputs
           n_AerBulk_Total_wfs, n_AerBulk_Group_wfs,  & ! OUTPUT, aerosol Jacobian Bookkeeping
           AerBulk_Mask_wfs, AerBulk_pars,            & ! OUTPUT, aerosol Jacobian Bookkeeping
           do_LinearRef, do_LinearPSD, Do_LinearEps,  & ! OUTPUT, aerosol Jacobian Bookkeeping
           NLIN, NPSD, Parslist,                      & ! OUTPUT, aerosol Jacobian Bookkeeping
           Fail_Bookkeep, Message_Bookkeep )            ! Exception handling
      if ( Fail_Bookkeep ) then
         do m = 1, 3
            Messages_Main(m+nmessages) = trim(Message_Bookkeep(m))
         enddo
         nmessages = nmessages + 3
         Fail_Main = .true. ; go to 4555
      endif
      ! write(*,*) 'Q,n_AerBulk_Total_wfs in Master: ', Q,n_AerBulk_Total_wfs
      Q = Q + n_AerBulk_Total_wfs
      ! write(*,*) Q
      ! pause
   endif
    
!  debug (OK)
!   write(*,*) 'aa', n_AerBulk_Total_wfs
!   write(*,*) 'aa2=',n_AerBulk_Group_wfs
!   do n = 1, 16
!      write(*,*),'bb-', AerBulk_Mask_wfs(n),AerBulk_pars(n)
!   enddo
!   pause'book'

!  Check on Number of weighting functions

   if ( do_AerBulk_Jacobians .and. Q.gt.MAX_ATMOSWFS ) then
      fail_Main = .true.
      messages_main(nmessages+1) = 'Number of Bulk Property Jacobians > Dimension MAXW_ATMOSWFS ==> Increase VLIDORT parameter'
      nmessages = nmessages + 1
      go to 4555
   endif

!  Final Column property tally

   if ( Q.gt.0 ) then
      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  =  .true.
      VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS         = Q
   endif

!  Surface jacobians - Lambertian albedo only

   if ( do_Surface_Jacobians ) then
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .true.
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 1  ! changed after GEMSTOOL_BRDF module
       !mchoi
   endif

!  If no linearization, then set the "simulation-only" flag
!  @@ Rob fix 7/23/14, Update for H2O scaling flag and value

   VLIDORT_LinModIn%MCont%TS_do_atmos_linearization   =  do_GasProfile_Jacobians.or.do_AerOpdepProfile_Jacobians.or.&
                                                         do_Tshift_Jacobian.or.do_Surfpress_Jacobian            .or.&
                                                         do_H2OScaling_Jacobian.or.do_CH4Scaling_Jacobian       .or.&
                                                         do_AerBulk_Jacobians
   VLIDORT_LinModIn%MCont%TS_do_linearization         =  VLIDORT_LinModIn%MCont%TS_do_atmos_linearization .or. &
                                                         do_Surface_Jacobians

   VLIDORT_LinModIn%MCont%TS_do_simulation_only        = .not.VLIDORT_LinModIn%MCont%TS_do_linearization
!   VLIDORT_LinFixIn%Cont%TS_do_simulation_only        = .not.VLIDORT_LinModIn%MCont%TS_do_linearization  ! Moved 10/17/16

!  Rob, 10/18/16, set VLIDORT linearization control for Flourescence Jacobians.
!               -->  These are pre-initialized Above.
!       11/30/16, Upgraded to include linear parameterization (2 Jacobians, not 1)

   n_SIFPars = 0 ; SIFPars = 0.0d0
   if ( GEMSTOOL_INPUTS%Sleave%DO_SIF ) then
      if ( GEMSTOOL_INPUTS%LinControl%do_SIF_Jacobians ) then
         VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS  = .true.
         if ( GEMSTOOL_INPUTS%Sleave%DO_SIF_LinearParam ) then
            SIFPars(1) = GEMSTOOL_INPUTS%Sleave%SIF755_Scaling_Constant
            SIFPars(2) = GEMSTOOL_INPUTS%Sleave%SIF755_Scaling_Gradient
            n_SIFPars  = 2
         else
            SIFPars(1) = GEMSTOOL_INPUTS%Sleave%SIF755_Amplitude
            n_SIFpars  = 1
         endif
         VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS = n_SIFPars
      endif
   endif

!  Bookkeeping to count the number of Profile weighting functions
!  Dont forget the aerosol profile !!!!!

   q = 0
   if (do_GasProfile_Jacobians ) then
      do g = 1, ngases
         if ( do_gases(g) .and. do_gas_wfs(g) ) then
            q = q + 1
         endif
      enddo
   endif
   if (do_AerOpdepProfile_Jacobians ) q = q + 1
   VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = q

!  Set up wavenumbers. Check Number of wavenumbers does not exceed MAXWAVNUMS

   range = GEMSTOOL_INPUTS%Wavenums%wavnum_finish - GEMSTOOL_INPUTS%Wavenums%wavnum_start 
   step  = range / GEMSTOOL_INPUTS%Wavenums%wavnum_res
   nwav  = INT(step) + 1
   if ( nwav .gt. maxwavnums) then
      fail_Main = .true.
      messages_main(nmessages+1) = 'Number of wavenumbers > Dimension MAXWAVNUMS ==> Increase MAXWAVNUMS'
      nmessages = nmessages + 1
      go to 4555
   endif


   do n = 1, nwav
      wav(n) = GEMSTOOL_INPUTS%Wavenums%wavnum_start + GEMSTOOL_INPUTS%Wavenums%wavnum_res * real ( n - 1, fpk )
   enddo
   ! write(*,*) wav(1), wav(nwav), step, nwav







   
   !  write(*,*) wav_ext(1), wav_ext(nwav_ext), step, nwav_ext
   ! pause

!   print*,'NWAV = ', nwav

!  =============================================================================
!                   ATMOSPHERIC and SURFACE PROPERTY CREATION
!  =============================================================================

!  1A. get GAS and PTH profiles
!  *****************************

!  12 August 2013. Use NSW-specific Profiles here

!  Assumes that the GRIDDING is fine enough for aerosols

!  @@ Rob fix 7/23/14, Update for H2O scaling flag and value

   call GEMSTOOL_NSW_PTHGAS_PROFILES &
           ( Maxlayers, MaxGases, ngases, which_gases, do_H2OScaling, do_CH4Scaling,   &  ! Input (control)
             gh2o, H2OScaling, gch4, CH4Scaling, Tshift, FDGas, FDLev, FDsfp, FDeps,   &  ! Input (numbers/Fd test)
             PhysicsDataPath, PTH_profile_name, GAS_profile_names,                     &  ! Input
             nlevels, level_heights, level_temps, level_press, level_vmrs,             &  ! Output (LEVELS)
             nlayers, aircolumns, daircolumns_dS, daircolumn_dP, Tshape,               &  ! Output (LAYERS). Air Density
             gascolumns, dgascolumns_dS, dgascolumn_dP, dgascolumns_dv,                &  ! Output (LAYERS). Gas densities
             dh2ocolumn_df, dch4column_df,                                             &  ! Output (LAYERS). Gas densities
             fail_Pthg, message_pthg )                                                    ! output

!  Surface pressure

   Surfpress = level_press(nlayers)
    write(*,*)'Surfpress  = ',  Surfpress
! pause
   
!  debug
!   do n = 0, nlayers
!      write(*,'(i4,1p7e17.7)')n,level_heights(n), level_temps(n), level_press(n), ( level_vmrs(n,k), k=1,4 )
!   enddo
!   pause'Debug PTHG'

!  Exception handling

   if ( fail_Pthg ) then
      fail_Main = .true.
      messages_main(nmessages+1) = adjustl(trim(message_pthg))
      nmessages = nmessages + 1
      go to 4555
   endif

!  set the cloud fraction proxy

!mick fix 6/18/2013 - added if condition
   if ( do_clouds ) cfrac = GEMSTOOL_INPUTS%Clouds%Cloud_fraction

!  Set the Cloud-top heights and pressures

   if ( do_clouds .and. GEMSTOOL_INPUTS%Clouds%do_Scattering_Clouds ) then
      if ( GEMSTOOL_INPUTS%Clouds%Cloud_boundaries ) then
         GEMSTOOL_INPUTS%Clouds%ctop_height   = level_heights(GEMSTOOL_INPUTS%Clouds%ctop_level)
         GEMSTOOL_INPUTS%Clouds%cbot_height   = level_heights(GEMSTOOL_INPUTS%Clouds%cbot_level)
         GEMSTOOL_INPUTS%Clouds%ctop_pressure = level_press  (GEMSTOOL_INPUTS%Clouds%ctop_level)
         GEMSTOOL_INPUTS%Clouds%cbot_pressure = level_press  (GEMSTOOL_INPUTS%Clouds%cbot_level)
      endif
   endif

!  set the profile linearization control

   if ( VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION ) then
      do n = 1, nlayers
         VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(n)   = .true.
         VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(n) = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
      enddo
   endif

!  Set the Column Linearization control

   if ( VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION ) then
      do n = 1, nlayers
         VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(n)   = .true.
         VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(n) = VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS
      enddo
   endif

!  1B. Get Gas and Rayleigh Cross-sections. Strictly level output for Gas Cross-sections.
!  ****************************************

!  Special for GEMS:
!     -- ALL CROSS_SECTIONS are on LEVELS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     -- O2/H2O/CO2/CH4 Cross-sections will be P-T dependent, derived from HITRAN

!  21 October 2013. develop handles for the LUT alternative.

   if ( do_Tshift_Jacobian .or.do_Surfpress_Jacobian ) then

      if ( do_Hitran ) then
         call GEMSTOOL_NSW_Cross_SECTIONS_PLUS &
          ( maxwavnums, maxlayers, maxgases, CO2_PPMV_MIXRATIO,     & ! INPUT
            nlayers, nwav, level_temps, level_press, wav,           & ! INPUT
            ngases, which_gases, hitran_path,                       & ! Input   !! Y.Jung
            do_Tshift_Jacobian, do_Surfpress_Jacobian,              & ! Input
            do_gases,                                               & ! In/Out  !! Y.Jung
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL,                          & ! OUTPUT
            GAS_XSECS, dGASXSECS_dT, dGASXSECS_dP,                  & ! OUTPUT
            fail_xsec, message_xsec )                                 ! output
      else
         call GEMSTOOL_NSW_Cross_SECTIONS_PLUS_LUT &
          ( maxwavnums, maxlayers, maxgases, CO2_PPMV_MIXRATIO,     & ! INPUT
            nlayers, nwav, level_temps, level_press, wav,           & ! INPUT
            ngases, which_gases, lut_path,                          & ! Input   !! Y.Jung
            do_Tshift_Jacobian, do_Surfpress_Jacobian,              & ! Input
            do_gases,                                               & ! In/Out  !! Y.Jung
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL,                          & ! OUTPUT
            GAS_XSECS, dGASXSECS_dT, dGASXSECS_dP,                  & ! OUTPUT
            fail_xsec, message_xsec )                                 ! output
      endif

      !If necessary, reduce the number of gas derivatives to be done
      !if a gas is not present in the spectral range
      do g = 1, ngases
         if ( do_GasProfile_Jacobians .and. .not.do_gases(g) &
              .and. do_gas_wfs(g) ) then
            do_gas_wfs(g) = .false.
            if ( VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION ) then
              VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = &
                VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS - 1
              VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:nlayers) = &
                VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
            endif
         end if
      end do
   else

     if ( do_Hitran ) then
        call GEMSTOOL_NSW_Cross_SECTIONS &
          ( maxwavnums, maxlayers, maxgases, CO2_PPMV_MIXRATIO,   & ! INPUT
            nlayers, nwav, level_temps, level_press, wav,         & ! INPUT
            ngases, which_gases, hitran_path,                     & ! Input     !! Y.Jung
            do_gases,                                             & ! In/Out    !! Y.Jung
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL, GAS_XSECS,             & ! OUTPUT
            fail_xsec, message_xsec )                               ! output
            !write(*,*)do_gases
            ! write(*,*) GAS_XSECS(1:2,1,2)
      else
        call GEMSTOOL_NSW_Cross_SECTIONS_LUT &
          ( maxwavnums, maxlayers, maxgases, CO2_PPMV_MIXRATIO,   & ! INPUT
            nlayers, nwav, level_temps, level_press, wav,         & ! INPUT
            ngases, which_gases, lut_path,                        & ! Input     !! Y.Jung
            do_gases,                                             & ! In/Out    !! Y.Jung
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL, GAS_XSECS,             & ! OUTPUT
            fail_xsec, message_xsec )                               ! output
      endif

   endif

!  Multiply Tshift derivatives by Shape factor

   if ( do_Tshift_Jacobian ) then
      do w = 1, nwav
         do n = 0, nlayers
            dGASXSECS_dT(w,n,1:ngases) = dGASXSECS_dT(w,n,1:ngases) * Tshape(n)
         enddo
      enddo
   endif

   if ( do_Surfpress_Jacobian ) then
      do w = 1, nwav
         dGASXSECS_dP(w,nlayers,1:ngases) = dGASXSECS_dP(w,nlayers,1:ngases) * 1013.25d0
      enddo
   endif

!  Debug.
!    Print out Rayleigh Xsec + depol, 4 gas sections for the lowest level.

  if ( XSEC_debug_flag ) then
     open(57,file='Debug_CrossSections_original.out',status = 'replace')
     do w = 1, nwav
        write(57,*)w,wav(w),RAYLEIGH_XSEC(w), RAYLEIGH_DEPOL(w)!, (GAS_XSECS(w,nlayers,g),g=1,4)
     enddo
     close(57) 
!   pause'Debug Cross sections'
  endif

!  Exception handling

   if ( fail_Xsec ) then
      fail_Main = .true.
      messages_main(nmessages+1) = adjustl(trim(message_xsec))
      nmessages = nmessages + 1
      go to 4555
   endif

!  2. get aerosol and cloud properties
!  ***********************************

!  if aerosols not present, zero the aerosol data

   if ( .not. do_aerosols ) then
       aerosol_layerflags = .false.
       Loading_aerosol    = 0.0d0
       aerosol_nscatmoms  = 0
       aertau_unscaled    = 0.0d0
       aod_scaling        = 0.0d0
       aerosol_deltau     = 0.0d0
       aerosol_ssalbs     = 0.0d0
       aerosol_scatmoms   = 0.0d0
       aerosol_distchars  = 0.0d0
       go to 654
    endif

!  Generation of Aerosol optical properties
!    @@@ User-defined aerosol option is now a part of the code.       10/21/16.
!    @@@ BOTH THESE AEROSOL ROUTINES HAVE BEEN EXTENSIVELY REWRITTEN. 10/21/16.

!     1.  aerosol Loading ( optical depth profile )
!     2a. (Option) Call to the Regular Mie  or Tmatrix programs
!     2b. (Option) Optical properties from User-defined Aerosol Read files. 
!     3.  Convert Mie/Tmatrix/User output (Microscopic) to IOP output (macroscopic)
!mick fix 1/24/2017 - added do_LinearRef, do_LinearPSD, Do_LinearEps, NLIN, & NPSD to
!                     "GEMSTOOL_AER_PROPERTIES_PLUS" input

   write(*,*) 'Doing aerosol'

   if ( do_AerBulk_Jacobians ) then
      call GEMSTOOL_AER_PROPERTIES_PLUS &
      ( maxlayers, maxwavnums, maxaermoms, maxaerwfs, interpolate_aerosols, do_wavnums, & ! Dimensions
        nlayers, nmuller, nwav, wav, level_heights, GEMSTOOL_INPUTS, momsize_cutoff,    & ! Inputs
        n_AerBulk_Total_wfs, n_AerBulk_Group_wfs, AerBulk_pars,                         & ! Inputs, aerosol Jacobian Bookkeeping
        do_LinearRef, do_LinearPSD, Do_LinearEps, NLIN, NPSD, ParsList,                 & ! Inputs, aerosol Jacobian Bookkeeping
        aerosol_layerflags, Loading_aerosol, Dloading_Dtau, Dloading_Dpars,             & ! OUTPUT, aerosol Loading
        aerosol_nscatmoms, aertau_unscaled, aod_scaling, aerosol_distchars,             & ! OUTPUT, aerosol optical properties
        aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms,                               & ! OUTPUT, aerosol optical properties
        L_aertau_unscaled, L_aod_scaling, L_aerosol_deltau,                             & ! OUTPUT, aerosol Linearized OPs
        L_aerosol_ssalbs, L_aerosol_scatmoms,                                           & ! OUTPUT, aerosol Linearized OPs
        fail1_Aer, fail2_Aer, Message_AerLoading, Message_AerOptical )                    ! Exception handling
   else
      call GEMSTOOL_AER_PROPERTIES &
      ( maxlayers, maxwavnums, maxaermoms,                                                & ! Dimensions
        interpolate_aerosols, do_wavnums, FDAer, FDLay, FDBul, FDeps,                     & ! Flags and FD control
        nlayers, nmuller, nwav, wav, level_heights, GEMSTOOL_INPUTS, momsize_cutoff,      & ! Inputs
        aerosol_layerflags, Loading_aerosol, aerosol_nscatmoms, aertau_unscaled,          & ! OUTPUT, aerosol control
        aod_scaling, aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms, aerosol_distchars, & ! OUTPUT, aerosol optical
        fail1_Aer, fail2_Aer, Message_AerLoading, Message_AerOptical )                      ! Exception handling
   endif

!  Debug aerosols
!   open(5,file='aod_profile.out',status='replace')
!   do n = 1, nlayers
!       write(5,55) n, aertau_unscaled(n),loading_aerosol(n)
!   enddo
!   close(5)
!55 format(i3,1p1000e15.5)

   open(10,file='./GEMSTOOL_NSW_Results/Aer_SSA.dat', status='replace')
   do ii = 1, nwav
       write(10,'(f11.5,2x,f15.8)')wav(ii), aerosol_ssalbs(ii)
   enddo
   close(10)
!stop'tempo'
   open(11,file='./GEMSTOOL_NSW_Results/Aer_prof.dat', status='replace')
!mick fix 6/19/2014 - modified extent of do loop for sake of aerosol_deltau
   !do ii = 0, nlayers
   do ii = 1, nlayers
       write(11,'(f11.5,1p30000e15.8)')level_heights(ii), (aerosol_deltau(ii,w), w=1, nwav)
   enddo
   close(11)

   !resave of Aer_prof with different form
   open(11,file='./GEMSTOOL_NSW_Results/Aer_prof_edit.dat', status='replace')
   do ii = 1, nwav
     write(11,'(f11.4,1p100e15.8)') wav(ii), (aerosol_deltau(w,ii), w=1,nlayers)
   enddo  
   close(11)


!  Exception handling

   if ( fail1_Aer ) then
      messages_Main(nmessages+1:nmessages+3) = message_AerLoading(1:3)
      nmessages = nmessages + 4
      messages_Main(nmessages) = 'Aerosol Loading failed (regridding first)- look at preceding 3 messages'
      go to 4555
   endif
   if ( fail2_Aer ) then
      messages_Main(nmessages+1:nmessages+5) = message_AerOptical(1:5)
      nmessages = nmessages + 6
      if ( do_AerBulk_Jacobians ) then
         messages_Main(nmessages) = 'Aerosol-Linearized Mie/Tmatrix failed (regridding first) - look at preceding 5 messages'
      else
         messages_Main(nmessages) = 'Aerosol-Regular Mie/Tmatrix failed (regridding first) - look at preceding 5 messages'
      endif
      go to 4555
   endif

!  Dimensioning check on number of moments
!    10/19/16, take into account wave number dependence

   if ( MAXVAL(aerosol_nscatmoms(1:nwav)) .gt. MAXMOMENTS_INPUT ) then
      fail_Main = .true.
      messages_main(nmessages+1) = 'Number of Aerosol Expansion Coefficients > MAXMOMENTS_INPUT'
      nmessages = nmessages + 1
      messages_main(nmessages+1) = 'Action: Increase MAXMOMENTS_INPUT, preferably to the same as maxaermoms'
      nmessages = nmessages + 1
      go to 4555
   endif

!  Continuation point for skipping aerosols

654 continue

!  if clouds not present, zero the cloud data and move on.

   if ( .not. do_clouds ) then
       cloud_layerflags = .false.
       Loading_cloud    = 0.0d0
       cldtau_unscaled  = 0.0d0
       cod_scaling      = 0.0d0      
       cloud_nscatmoms  = 0
       cloud_deltau     = 0.0d0
       cloud_ssalbs     = 0.0d0
       cloud_scatmoms   = 0.0d0
       go to 664
    endif

!  Generation of Cloud optical properties
!     1. Call to the Regular Mie program
!     2. Convert Mie output (Microsopic) to IOP output (macroscopic)

   write(*,*) 'Doing cloud'

   call GEMSTOOL_CLD_PROPERTIES &
      ( maxlayers, maxwavnums, maxcldmoms, interpolate_clouds, do_wavnums,             & ! Dimensions and Flags 
        nlayers, nmuller, nwav, wav, level_heights, GEMSTOOL_INPUTS, momsize_cutoff,   & ! Inputs3x
        cloud_layerflags, Loading_cloud, cloud_nscatmoms, cldtau_unscaled,             & ! OUTPUT, aerosol control
        cod_scaling, cloud_deltau, cloud_ssalbs, cloud_scatmoms,                       & ! OUTPUT, aerosol optical
        fail1_Cld, Message_CLdOptical )

!  Exception handling

   if ( fail1_Cld ) then
      messages_Main(nmessages+1:nmessages+4) = message_CldOptical(1:4)
      nmessages = nmessages + 5
      messages_Main(nmessages) = 'Cloud Mie calculation failed - look at preceding 5 messages'
      go to 4555
   endif

!  Dimensioning check on number of moments

   if ( cloud_nscatmoms .gt. MAXMOMENTS_INPUT ) then
      fail_Main = .true.
      messages_main(nmessages+1) = 'Number of Cloud Expansion Coefficients > MAXMOMENTS_INPUT'
      nmessages = nmessages + 1
      messages_main(nmessages+1) = 'Action: Increase MAXMOMENTS_INPUT, preferably to the same as maxcldmoms'
      nmessages = nmessages + 1
      go to 4555
   endif

!  Continuation point for skipping clouds

664 continue

!  3. Prepare surface property
!  ***************************

!  Zero the Lambertian albedo

   lambertian_albedo = 0.0d0  

!  Zero the Supplement results
    
!   GEMSTOOL_BRDF_Results   = 0.0d0
   GEMSTOOL_SLEAVE_Results%SL_SLTERM_ISOTROPIC  = 0.0d0
   GEMSTOOL_SLEAVE_Results%SL_SLTERM_USERANGLES = 0.0d0
   GEMSTOOL_SLEAVE_Results%SL_SLTERM_F_0        = 0.0d0
   GEMSTOOL_SLEAVE_Results%SL_USER_SLTERM_F_0   = 0.0d0

!  BRDF for Land Surfaces
!  ----------------------

!  10/25/16. Initial Version for Surfaces that are NOT OCEAN
!    -- Same BRDF for all wavelengths, done outside the main wavelength loop.
!    -- Up to 3 error messages allowed.

   do_BRDF = GEMSTOOL_INPUTS%BRDF%do_Gemstool_BRDF

   if ( GEMSTOOL_INPUTS%BRDF%do_Gemstool_BRDF ) then

      CALL GEMSTOOL_BRDFSurface &
      ( GEMSTOOL_INPUTS, BRDF_N_SURFACE_WFS, BRDF_K_FP_values, &
        GEMSTOOL_BRDF_Results, GEMSTOOL_BRDF_Results_Lin, & ! Input and output ! added by mchoi
        fail_BRDF, nmessages_BRDF, messages_BRDF ) ! Exception handling
      !   write(*,*) 'LinIQUMaster: BRDF_N_SURFACE_WFS ',BRDF_N_SURFACE_WFS   
      !   pause !mchoi
      if ( fail_BRDF ) then
         fail_Main = .true.
         do m = 1, nmessages_brdf
            messages_main(nmessages+m) = Trim(messages_BRDF(m))
         enddo
         nmessages = nmessages + nmessages_brdf
         messages_main(nmessages+1) = 'GEMSTOOL Land-surface BRDF calculation failed - look at preceding messages'
         nmessages = nmessages + 1
         go to 4555
      endif
      ! write(*,*) 'LinIQUMaster: BRDF_N_SURFACE_WFS ',BRDF_N_SURFACE_WFS
      ! write(*,*) BRDF_K_FP_values  
      ! pause !mchoi

      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS  = BRDF_N_SURFACE_WFS
        
   endif

   N_SURFACE_WFS_fin = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS


   !  Lambertian albedo
!  -----------------

   if ( .not. GEMSTOOL_INPUTS%BRDF%do_Gemstool_BRDF ) then

     call GEMSTOOL_Albedo &
        ( Maxwavnums, MaxAlbCoeffs, nwav, wav, do_wavnums, GEMSTOOL_INPUTS, & ! inputs
          PhysicsDataPath, alb_filename,                                    & ! inputs
          lambertian_albedo, ClosureTerms )                                   ! output

!  Debug
!  do w = 1, nwav
!     write(*,*)w,wav_shift_WN(ww),lambertian_albedo(w)
!  enddo
!  pause

   endif

!  ===========================================================================
!          CALCULATION SECTION : MAIN WAVELENGTH LOOP WITH VLIDORT
!  ===========================================================================

!  TOA or BOA output level

  if ( VLIDORT_FixIn%Bool%TS_do_upwelling ) then
     VLIDORT_ModIn%MUserVal%TS_user_levels(1) = GEMSTOOL_INPUTS%RTControl%observation_height !zero
     !VLIDORT_ModIn%MUserVal%TS_user_levels(2) = 3 !GEMSTOOL_INPUTS%RTControl%observation_height !zero
  else if ( VLIDORT_FixIn%Bool%TS_do_dnwelling ) then
     VLIDORT_ModIn%MUserVal%TS_user_levels(1) = dble(nlayers)
  endif

!  Set up IPA. Number of Scenes.

   local_cloud(1:2)  = .false.
   scene_weight(1:2) = 1.0_fpk
!mick fix 6/18/2014 - switched to permanent variable "do_clouds" now
   !if ( .not. do_clouds_V2 ) then
   if ( .not. do_clouds ) then
      n_scenes = 1   ! Clear sky calculation, no clouds
   else
      if ( cfrac .eq. 0.0_fpk ) then
         n_scenes = 1   ! Clear sky calculation, cloud fraction = 0
      else if ( cfrac .eq. 1.0_fpk ) then
         n_scenes = 1 ;  local_cloud(1) = .true. ! Fully cloudy calculation, cloudfraction = 1
      else
         n_scenes = 2 ;  local_cloud(2) = .true.  ! Partially cloudy
         scene_weight(1) = 1.0_fpk - cfrac
         scene_weight(2) = cfrac
      endif
   endif

!  Zero the Final output

   STOKES  = 0.0_fpk
   DOLP    = 0.0_fpk
   ACTINIC = 0.0_fpk
   REGFLUX = 0.0_fpk

   RADIANCE_GASWFS      = 0.0_fpk
   RADIANCE_AERPROFWFS  = 0.0_fpk
   RADIANCE_AERBULKWFS  = 0.0_fpk ! New, June 2014
   RADIANCE_ALBEDOWFS   = 0.0_fpk
   RADIANCE_TSHIFTWFS   = 0.0_fpk
   RADIANCE_SURFPWFS    = 0.0_fpk
   RADIANCE_WSCALEWFS   = 0.0_fpk ! New 23 July 2014
   RADIANCE_MSCALEWFS   = 0.0_fpk ! New 01 Feb  2015
   RADIANCE_SIFPARAMWFS = 0.0_fpk ! New, 10/18/16, 11/30/16    ! SIF introduction

   AMFSCATWTS          = 0.0_fpk
   AMFTOTAL            = 0.0_fpk

!  Stokes Q/U Jacobians, initialization added 12/19/16 by Rob

   STOKES_Q_GASWFS      = 0.0_fpk ; STOKES_U_GASWFS      = 0.0_fpk
   STOKES_Q_ALBEDOWFS   = 0.0_fpk ; STOKES_U_ALBEDOWFS   = 0.0_fpk
   STOKES_Q_SIFPARAMWFS = 0.0_fpk ; STOKES_U_SIFPARAMWFS = 0.0_fpk
   STOKES_Q_TSHIFTWFS   = 0.0_fpk ; STOKES_U_TSHIFTWFS   = 0.0_fpk
   STOKES_Q_SURFPWFS    = 0.0_fpk ; STOKES_U_SURFPWFS    = 0.0_fpk
   STOKES_Q_WSCALEWFS   = 0.0_fpk ; STOKES_U_WSCALEWFS   = 0.0_fpk
   STOKES_Q_MSCALEWFS   = 0.0_fpk ; STOKES_U_MSCALEWFS   = 0.0_fpk
   STOKES_Q_AERPROFWFS  = 0.0_fpk ; STOKES_U_AERPROFWFS  = 0.0_fpk
   STOKES_Q_AERBULKWFS  = 0.0_fpk ; STOKES_U_AERBULKWFS  = 0.0_fpk

   dir = UpDn_Index

!  Start the CLoud Scene loop. 1 or 2 Scenes (Clear-sky, Cloudy)

   do i_scene = 1, n_scenes

!mick mod 1/24/2017 - added more descriptive feedback to user
      write(*,*)
      if ( .not.do_clouds ) then
         write(*,*)'Doing clear scene'
      else
         if ( cfrac .eq. 0.0_fpk ) then
            write(*,*)'Doing clear scene'
         else if ( cfrac .eq. 1.0_fpk ) then
            write(*,*)'Doing cloudy scene'
         else
            if (i_scene .eq. 1) then
               write(*,*)'Doing partly cloudy scene: clear portion'
            else
               write(*,*)'Doing partly cloudy scene: cloudy portion'
            endif
         endif
      endif



!  Start MAIN WAVELENGTH/WAVENUMBER LOOP WITH VLIDORT
!  --------------------------------------------------
      ! open(340,file='./GEMSTOOL_NSW_Results/GASOPDEPS.dat',status='replace')

!      nwav = 1
      DO w = 1, nwav

!  progress

        if (w == 1) write(*,*)
        !if (mod(w,10000).eq.0) write(*,*)' Linearized VLIDORT Calculation,  doing wavelength/wavenumber # ',w
        if (mod(w,1).eq.0) write(*,*)' Linearized VLIDORT Calculation,  doing wavelength/wavenumber # ',w

!  Fluorescence 
!  ------------

        Polynomial = 1.0d0 ; dPolynomial_dSL = 0.0d0 ; SIF_ExactOnly = .false.   ! Initialized for safety

!  New code to implement SIF, 10/18/16. tested 10/21/16. Updated 10/26/15

        if ( GEMSTOOL_INPUTS%Sleave%do_SIF ) then

!  ** Upgrade 11/30/16. Only one call if you are using the linear parameterization
!                       In this case you are just getting SIF 755 value.

            if ( ( GEMSTOOL_INPUTS%Sleave%DO_SIF_LinearParam.and.w.eq.1) .or. &
                 .not. GEMSTOOL_INPUTS%Sleave%DO_SIF_LinearParam ) then
               Call GEMSTOOL_Fluorescence &
                ( do_wavnums, wav(w), GEMSTOOL_INPUTS, SIF_datapath, &
                  GEMSTOOL_SLEAVE_Results )
            endif

!   Upgrade 11/30/16. calculate polynomial for the SIF linear parameterization

            SIF_ExactOnly = GEMSTOOL_INPUTS%Sleave%DO_SIF_ExactOnly
            if ( GEMSTOOL_INPUTS%Sleave%DO_SIF_LinearParam ) then
               factor =  wav(w) ; if ( do_wavnums ) factor = 1.0d+07/wav(w)
               factor = 1.0d0 - ( factor / 755.0d0 )
               Polynomial = SIFPars(1) + SIFPars(2) * Factor
               dPolynomial_dSL(1) = 1.0d0
               dPolynomial_dSL(2) = Factor
            else
               Polynomial = 1.0d0
               dPolynomial_dSL(1) = 1.0d0 / SIFpars(1)
            endif

         endif

!  Develop the VLIDORT-specific optical properties
!  -----------------------------------------------

!  @@ Rob fix 7/23/14, add H2O scaling derivative

         if ( (do_GasProfile_Jacobians.or.do_AerOpdepProfile_Jacobians) .or. &
              (do_Tshift_Jacobian.or.do_Surfpress_Jacobian.or.do_H2OScaling_Jacobian.or.do_CH4Scaling_Jacobian .or. &
               do_AerBulk_Jacobians) ) then
            call GEMSTOOL_NSW_linearized_iops &
        ( maxlayers, maxmoments_input, maxstokes_sq, max_atmoswfs,                 & ! VLIDORT dimensions
          maxwavnums, maxaermoms, maxcldmoms, maxgases, maxaerwfs,                 & ! GEMSTOOL dimensions
          IOP_debug_filename, IOP_Debug_flag,                                      & ! Debug IOP check  
          do_GasProfile_Jacobians, do_AerOpdepProfile_Jacobians,                   & ! Profile Jacobian control flags
          do_AerBulk_Jacobians, do_Tshift_Jacobian,                                & ! Column/Bulk Jacobian control flags
          do_Surfpress_Jacobian, do_H2OScaling_Jacobian, do_CH4Scaling_Jacobian,   & ! Column/Bulk Jacobian control flags
          w, wav(w), omega_lim, do_gases, do_gas_wfs, ngases,                      & ! GEMSTOOL Inputs (Control)
          nstokes, nlayers, nwav, n_AerBulk_Total_wfs, AerBulk_pars,               & ! GEMSTOOL Inputs (Control)
          do_aerosols, aerosol_nscatmoms, local_cloud(i_scene), cloud_nscatmoms,   & ! GEMSTOOL Inputs (Aer/clouds)
          gascolumns, dGASCOLUMNS_dS, dGASCOLUMN_dP, dGASCOLUMNS_dV,               & ! GEMSTOOL Inputs (GasColumnsderivs)
          dH2OCOLUMN_dF, dCH4COLUMN_dF,                                            & ! GEMSTOOL Inputs (GasColumnsderivs)
          aircolumns, dAIRCOLUMNS_dS, dAIRCOLUMN_dP,                               & ! GEMSTOOL Inputs (AirColumns + derivs)
          Tshift, SurfPress, H2OScaling, gh2o, CH4Scaling, gch4,                   & ! GEMSTOOL Inputs (T-shift,SurfP,Scaling)
          gas_xsecs, dgasxsecs_dT, dgasxsecs_dP, Rayleigh_xsec, Rayleigh_depol,    & ! GEMSTOOL Inputs (Cross-sections, depol)
          aerosol_layerflags, aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms,    & ! GEMSTOOL Inputs (Aerosol)
          L_aerosol_deltau, L_aerosol_ssalbs, L_aerosol_scatmoms,                  & ! GEMSTOOL Inputs (Aerosol linearized)
          cloud_layerflags, cloud_deltau, cloud_ssalbs, cloud_scatmoms,            & ! GEMSTOOL Inputs (Clouds)
          DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT,                                    & ! VLIDORT Outputs (Regular iops)
          NGREEK_MOMENTS_INPUT, GREEKMAT_TOTAL_INPUT,                              & ! VLIDORT Outputs (Regular iops)
          L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,        & ! VLIDORT Outputs (linearized)
          LAYER_GASOPDEPS, dGASOPDEPS_dV )                                           ! Layer-to-Level transformation
         else
            call GEMSTOOL_NSW_regular_iops &
        ( maxlayers, maxmoments_input, maxstokes_sq,                            & ! VLIDORT dimensions
          maxwavnums, maxaermoms, maxcldmoms, maxgases,                         & ! GEMSTOOL dimensions
          IOP_debug_filename, IOP_debug_flag,                                   & ! Debug IOP check  
          w, wav(w), omega_lim, do_gases, ngases, nstokes, nlayers, nwav,               & ! GEMSTOOL Inputs
          do_aerosols, aerosol_nscatmoms, local_cloud(i_scene), cloud_nscatmoms,        & ! GEMSTOOL Inputs
          gas_xsecs,                                                            & ! GEMSTOOL Inputs
          gascolumns, aircolumns, Rayleigh_xsec, Rayleigh_depol,                & ! GEMSTOOL Inputs
          aerosol_layerflags, aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms, & ! GEMSTOOL Inputs
          cloud_layerflags, cloud_deltau, cloud_ssalbs, cloud_scatmoms,         & ! GEMSTOOL Inputs
          DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT,                                 & ! VLIDORT Outputs
          NGREEK_MOMENTS_INPUT, GREEKMAT_TOTAL_INPUT )                            ! VLIDORT Outputs
         endif
         ! ! write(*,*) w, wav(w), nwav
         ! write(*,*) w, wav(w), omega_lim, do_gases, ngases, nstokes, nlayers, nwav
         ! write(*,*) GASCOLUMNS(:,1,:) ! GAS_XSECS(1,:,1)
         ! write(*,*)DELTAU_VERT_INPUT
   ! write(*,*) NGREEK_MOMENTS_INPUT
   ! write(*,*) GREEKMAT_TOTAL_INPUT(0:30,20, 1)
   
         ! write(*,*) NGREEK_MOMENTS_INPUT
         ! pause

   ! !put some results to the LRRS preparation ----not here
   ! do i_layer = 1, maxlayers
   !    DELTAU_INPUT_UNSCALED(i_layer,w) = DELTAU_VERT_INPUT(i_layer)
   !    do i_moment = 0, MAXMOMENTS_INPUT
   !       OMEGAMOMS_ELASTIC_UNSCALED(i_layer,i_moment,w)=GREEKMAT_TOTAL_INPUT(i_moment,i_layer,1)
   !    enddo
   !    ! write(*,*) DELTAU_INPUT_UNSCALED(i_layer,w)
   !    ! write(*,*) OMEGAMOMS_ELASTIC_UNSCALED(i_layer, 0:20, w)
   ! enddo




!  Do the RT Calculation (2 parts: Actual and Final)
!  -------------------------------------------------

!  @@ Rob   fix 7/23/14, add H2O scaling derivative
!  @@ YJung fix 2/1/15 , add CH4 scaling derivative
!  @@ Rob, 10/18/16, Introduce SIF arguments (prepared SIF data) --> replaced by SLEAVE structure
!                  --> Include flag for using SIF Jacobians (renamed 11/30/16)
!  @@ Rob, 10/25/16, Introduce BRDF and SLEAVE Type structures
!  @@ Rob, 10/25/16, Renamed module GEMSTOOL_NSW_SetJacobians
!  @@ Rob, 11/30/16, Introduce SIF inputs (for SIF linear parameterization). revised naming.
!  @@ Rob, 12/19/16, Introduce Stokes-Q and Stokes-U output via subroutine GEMSTOOL_NSW_SetJacobians_IQU

         if ( ( do_GasProfile_Jacobians.or.do_AerOpdepProfile_Jacobians ) .or. &
              ( do_Tshift_Jacobian.or.do_Surfpress_Jacobian.or.do_H2OScaling_Jacobian.or.do_CH4Scaling_Jacobian .or. &
                do_AerBulk_Jacobians ) .or. &
              ( do_Surface_Jacobians .or. do_SIF_Jacobians ) ) then

            call GEMSTOOL_RTCALC_PLUS_Actual &
           (  do_GasProfile_Jacobians, do_AerOpdepProfile_Jacobians,                    & ! Main control
              GEMSTOOL_INPUTS%Sleave%do_SIF,Polynomial, dPolynomial_dSL, SIF_ExactOnly, & ! SIF control
              VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                            & ! VLIDORT regular inputs
              VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup,                   & ! VLIDORT linearized inputs
              nlayers, NGREEK_MOMENTS_INPUT, level_heights, lambertian_albedo(w),   & ! Control/Surface Proxies
              DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,           & ! Atmos-Optical Proxies
              L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,     & ! Atmos-Optical Proxies (linearized)
              GEMSTOOL_SLEAVE_Results, GEMSTOOL_BRDF_Results,GEMSTOOL_BRDF_Results_Lin, & ! Surface Supplement Results !mchoi
              VLIDORT_Out, VLIDORT_LinOut )                                           ! VLIDORT Output results

            call GEMSTOOL_NSW_SetJacobians_IQU &
         (  MAXWAVNUMS, MaxAlbCoeffs, maxgases, maxmessages, maxaerwfs, W, dir, do_SphericalAlbedo,      & ! GEMSTOOL Dimensions and W
            do_normalized_wfoutput, do_GasProfile_Jacobians, do_AerOpdepProfile_Jacobians,               & ! GEMSTOOL Jacobian control
            do_Surface_Jacobians, do_SIF_Jacobians, do_AerBulk_Jacobians,                                & ! GEMSTOOL Jacobian control
            do_Tshift_Jacobian, do_Surfpress_Jacobian, do_H2OScaling_Jacobian, do_CH4Scaling_Jacobian,   & ! GEMSTOOL Jacobian control
            do_gases, do_gas_wfs, ngases, nlayers, nstokes, n_geometries, n_AerBulk_Total_wfs, n_SIFPars,& ! Control numbers/gases
            Lambertian_albedo, SIFPars, Tshift, SurfPress, H2OScaling, CH4Scaling,                       & ! necessary for normaliz
            Level_VMRs, AERTAU_UNSCALED, AerBulk_pars,                                                   & ! necessary for normaliz
            VLIDORT_Out, VLIDORT_LinOut,                                              & ! VLIDORT Results (includes linearized)
            LAYER_GASOPDEPS, dGASOPDEPS_dV,                                           & ! Layer-to-Level transformations
            STOKES_ic, ACTINIC_ic, REGFLUX_ic, DOLP_ic, DOCP_ic,                      & ! Main program, GEMSTOOL Stokes Output
            RADIANCE_GASWFS_ic,    RADIANCE_AERPROFWFS_ic,  RADIANCE_AERBULKWFS_ic,   & ! Main program, GEMSTOOL I-Jacobians Output
            RADIANCE_ALBEDOWFS_ic, RADIANCE_SIFPARAMWFS_ic, RADIANCE_TSHIFTWFS_ic,    & ! Main program, GEMSTOOL I-Jacobians Output
            RADIANCE_SURFPWFS_ic,  RADIANCE_WSCALEWFS_ic,   RADIANCE_MSCALEWFS_ic,    & ! Main program, GEMSTOOL I-Jacobians Output
            AMFSCATWTS_ic, AMFTOTAL_ic, TOTALWF_SAVED_ic,                             & ! Main program, GEMSTOOL AMF Output
            STOKES_Q_GASWFS_ic,    STOKES_Q_AERPROFWFS_ic,  STOKES_Q_AERBULKWFS_ic,   & ! Main program, GEMSTOOL Q-Jacobians Output
            STOKES_Q_ALBEDOWFS_ic, STOKES_Q_SIFPARAMWFS_ic, STOKES_Q_TSHIFTWFS_ic,    & ! Main program, GEMSTOOL Q-Jacobians Output
            STOKES_Q_SURFPWFS_ic,  STOKES_Q_WSCALEWFS_ic,   STOKES_Q_MSCALEWFS_ic,    & ! Main program, GEMSTOOL Q-Jacobians Output
            STOKES_U_GASWFS_ic,    STOKES_U_AERPROFWFS_ic,  STOKES_U_AERBULKWFS_ic,   & ! Main program, GEMSTOOL U-Jacobians Output
            STOKES_U_ALBEDOWFS_ic, STOKES_U_SIFPARAMWFS_ic, STOKES_U_TSHIFTWFS_ic,    & ! Main program, GEMSTOOL U-Jacobians Output
            STOKES_U_SURFPWFS_ic,  STOKES_U_WSCALEWFS_ic,   STOKES_U_MSCALEWFS_ic,    & ! Main program, GEMSTOOL U-Jacobians Output
            Errorstatus, nmessages, messages_Main, &                                     ! Main program, exception handling
            do_BRDF, N_SURFACE_WFS_fin, BRDF_K_FP_values, &   !) !mchoi
            RADIANCE_BRDFWFS_ic, STOKES_Q_BRDFWFS_ic,STOKES_U_BRDFWFS_ic) !mchoi
         else

            call GEMSTOOL_RTCALC_Actual &
            ( VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                          & ! VLIDORT Actual inputs
              Polynomial, SIF_ExactOnly,                                          & ! SIF control
              nlayers, NGREEK_MOMENTS_INPUT, level_heights, lambertian_albedo(w), & ! Control Proxies
              DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,         & ! Atmos-Optical Proxies
              GEMSTOOL_SLEAVE_Results, GEMSTOOL_BRDF_Results,                     & ! Surface Supplement Results
              VLIDORT_Out )                                                         ! VLIDORT output results
              
            call GEMSTOOL_RTCALC_Final &
            ( MAXWAVNUMS, maxmessages, W, do_SphericalAlbedo,      & ! GEMSTOOL control
              nstokes, n_geometries, dir, VLIDORT_Out,             & ! VLIDORT Results
              STOKES_ic, ACTINIC_ic, REGFLUX_ic, DOLP_ic, DOCP_ic, & ! Main program, GEMSTOOL Stokes Output
              Errorstatus, nmessages, messages_Main )                ! Main program, exception handling

         endif

!      if ( w.eq.1 ) then
!         write(*,*)aerosol_deltau(18,w)
!         write(*,*) VLIDORT_LinOut%Prof%TS_PROFILEWF(5,18,1,1,1,1), VLIDORT_Out%Main%TS_STOKES(1,1,1,1)
!      endif

!  exit if serious error

         if ( Errorstatus .eq. 1 ) then
            write(c1,'(I1)')i_scene
            messages_Main(nmessages+1) = 'Failure for scene # '//c1
            nmessages = nmessages+1
            go to 4555
         endif

!  Add to the total (Apply IPA)
!  ============================

!  This section updated 12/19/16 by Rob, to include Stokes-Q and Stokes-U Jacobians

!  Stokes vectors
         STOKES(1:n_geometries,1:nstokes,w) = STOKES(1:n_geometries,1:nstokes,w) + &
         scene_weight(i_scene) * STOKES_ic(1:n_geometries,1:nstokes,w)
         if ( nstokes .gt. 1 )  DOLP(1:n_geometries,w) = DOLP(1:n_geometries,w) + scene_weight(i_scene) * DOLP_ic(1:n_geometries,w)

!  Gas Profile Jacobians

         if ( do_GasProfile_Jacobians ) then
            do g = 1, ngases
               do n = 0, nlayers
                  radiance_gaswfs(1,w,n,g) = radiance_gaswfs(1,w,n,g)  + scene_weight(i_scene) * radiance_gaswfs_ic(1,w,n,g)
                  if (n.gt.0)AMFSCATWTS(1,w,n,g) = AMFSCATWTS(1,w,n,g) + scene_weight(i_scene) * AMFSCATWTS_ic(1,w,n,g)
               enddo
               if ( nstokes .gt. 1 ) then
                 do n = 0, nlayers
                   Stokes_Q_gaswfs(1,w,n,g) = Stokes_Q_gaswfs(1,w,n,g)  + scene_weight(i_scene) * Stokes_Q_gaswfs_ic(1,w,n,g)
                   Stokes_U_gaswfs(1,w,n,g) = Stokes_U_gaswfs(1,w,n,g)  + scene_weight(i_scene) * Stokes_U_gaswfs_ic(1,w,n,g)
                 enddo
               endif 
               AMFTOTAL(1,w,g)      = AMFTOTAL(1,w,g)      + scene_weight(i_scene) * AMFTOTAL_ic(1,w,g)
            enddo
         endif

!  Aerosol profile Jacobians

         if ( do_AerOpdepProfile_Jacobians ) then
            do n = 1, nlayers
               radiance_aerprofwfs(1,w,n) = radiance_aerprofwfs(1,w,n) + scene_weight(i_scene) * radiance_aerprofwfs_ic(1,w,n)
               if ( nstokes .gt. 1 ) then
                 Stokes_Q_aerprofwfs(1,w,n) = Stokes_Q_aerprofwfs(1,w,n)  + scene_weight(i_scene) * Stokes_Q_aerprofwfs_ic(1,w,n)
                 Stokes_U_aerprofwfs(1,w,n) = Stokes_U_aerprofwfs(1,w,n)  + scene_weight(i_scene) * Stokes_U_aerprofwfs_ic(1,w,n)
               endif 
            enddo
         endif

!  Surface Jacobians !mchoi

         if ( do_Surface_Jacobians ) then
            if (.not. do_BRDF) then

               radiance_Albedowfs(1:n_geometries,w,1) = radiance_Albedowfs(1:n_geometries,w,1) + &
                                    scene_weight(i_scene) * radiance_Albedowfs_ic(1:n_geometries,w,1)  
               if ( nstokes .gt. 1 ) then
                  Stokes_Q_Albedowfs(1:n_geometries,w,1) = Stokes_Q_Albedowfs(1:n_geometries,w,1) + &
                                    scene_weight(i_scene) * Stokes_Q_Albedowfs_ic(1:n_geometries,w,1)  
                  Stokes_U_Albedowfs(1:n_geometries,w,1) = Stokes_U_Albedowfs(1:n_geometries,w,1) + &
                                    scene_weight(i_scene) * Stokes_U_Albedowfs_ic(1:n_geometries,w,1)  
               endif
            else  !if (do_BRDF) then
               do qs = 1, N_SURFACE_WFS_fin
                  radiance_BRDFwfs(1:n_geometries,w,qs) = radiance_BRDFwfs(1:n_geometries,w,qs) + &
                                       scene_weight(i_scene) * radiance_BRDFwfs_ic(1:n_geometries,w,qs)  
                  if ( nstokes .gt. 1 ) then
                     Stokes_Q_BRDFwfs(1:n_geometries,w,qs) = Stokes_Q_BRDFwfs(1:n_geometries,w,qs) + &
                                       scene_weight(i_scene) * Stokes_Q_BRDFwfs_ic(1:n_geometries,w,qs)  
                     Stokes_U_BRDFwfs(1:n_geometries,w,qs) = Stokes_U_BRDFwfs(1:n_geometries,w,qs) + &
                                       scene_weight(i_scene) * Stokes_U_BRDFwfs_ic(1:n_geometries,w,qs)  
                  endif
               enddo !do qs = 1, n_SIFPars

               ! write(*,*)'FINAL'
               ! write(*,*)radiance_BRDFwfs(1:n_geometries,w,:)
               ! pause
            endif !if (.not. do_BRDF) then

          endif

!  T-shift Jacobians

         if ( do_Tshift_Jacobian ) then
            radiance_Tshiftwfs(1:n_geometries,w) = radiance_Tshiftwfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * radiance_Tshiftwfs_ic(1:n_geometries,w)  
            if ( nstokes .gt. 1 ) then
               Stokes_Q_Tshiftwfs(1:n_geometries,w) = Stokes_Q_Tshiftwfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * Stokes_Q_Tshiftwfs_ic(1:n_geometries,w)  
               Stokes_U_Tshiftwfs(1:n_geometries,w) = Stokes_U_Tshiftwfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * Stokes_U_Tshiftwfs_ic(1:n_geometries,w)  
            endif
         endif

!  Surface pressure Jacobians

         if ( do_Surfpress_Jacobian ) then
            radiance_Surfpwfs(1:n_geometries,w) = radiance_Surfpwfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * radiance_Surfpwfs_ic(1:n_geometries,w)  
            if ( nstokes .gt. 1 ) then
               Stokes_Q_Surfpwfs(1:n_geometries,w) = Stokes_Q_Surfpwfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * Stokes_Q_Surfpwfs_ic(1:n_geometries,w)  
               Stokes_U_Surfpwfs(1:n_geometries,w) = Stokes_U_Surfpwfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * Stokes_U_Surfpwfs_ic(1:n_geometries,w)  
            endif
         endif

!  @@ Rob fix 7/23/14, add H2O scaling derivative

         if ( do_H2OScaling_Jacobian ) then
            radiance_WScalewfs(1:n_geometries,w) = radiance_WScalewfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * radiance_WScalewfs_ic(1:n_geometries,w)  
            if ( nstokes .gt. 1 ) then
               Stokes_Q_WScalewfs(1:n_geometries,w) = Stokes_Q_WScalewfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * Stokes_Q_WScalewfs_ic(1:n_geometries,w)  
               Stokes_U_WScalewfs(1:n_geometries,w) = Stokes_U_WScalewfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * Stokes_U_WScalewfs_ic(1:n_geometries,w)  
            endif
         endif

!  @@ Y.Jung 2/1/15, add CH4 scaling derivative

         if ( do_CH4Scaling_Jacobian ) then
            radiance_MScalewfs(1:n_geometries,w) = radiance_MScalewfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * radiance_MScalewfs_ic(1:n_geometries,w)  
            if ( nstokes .gt. 1 ) then
               Stokes_Q_MScalewfs(1:n_geometries,w) = Stokes_Q_MScalewfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * Stokes_Q_MScalewfs_ic(1:n_geometries,w)  
               Stokes_U_MScalewfs(1:n_geometries,w) = Stokes_U_MScalewfs(1:n_geometries,w) + &
                                  scene_weight(i_scene) * Stokes_U_MScalewfs_ic(1:n_geometries,w)  
            endif
         endif

!  @@ Rob Fix 10/18/16, add SIF derivative
!             11/30/16, update for linear parameteriation, rename output arrays --> SIFParamWFs

         if ( do_SIF_Jacobians ) then
            do qs = 1, n_SIFPars
               radiance_SIFParamwfs(1:n_geometries,w,qs) = radiance_SIFParamwfs(1:n_geometries,w,qs) + &
                                  scene_weight(i_scene) * radiance_SIFParamwfs_ic(1:n_geometries,w,qs)
               if ( nstokes .gt. 1 ) then
                  Stokes_Q_SIFParamwfs(1:n_geometries,w,qs) = Stokes_Q_SIFParamwfs(1:n_geometries,w,qs) + &
                                  scene_weight(i_scene) * Stokes_Q_SIFParamwfs_ic(1:n_geometries,w,qs)
                  Stokes_U_SIFParamwfs(1:n_geometries,w,qs) = Stokes_U_SIFParamwfs(1:n_geometries,w,qs) + &
                                  scene_weight(i_scene) * Stokes_U_SIFParamwfs_ic(1:n_geometries,w,qs)
              endif
            enddo
         endif

!  Bulk Aerosol Jacobians

         if ( do_AerBulk_Jacobians ) then
            do q = 1, n_AerBulk_Total_wfs
               radiance_aerbulkwfs(1,w,q) = radiance_aerbulkwfs(1,w,q) + scene_weight(i_scene) * radiance_aerbulkwfs_ic(1,w,q)
               if ( nstokes .gt. 1 ) then
                  Stokes_Q_aerbulkwfs(1,w,q) = Stokes_Q_aerbulkwfs(1,w,q) + scene_weight(i_scene) * Stokes_Q_aerbulkwfs_ic(1,w,q)
                  Stokes_U_aerbulkwfs(1,w,q) = Stokes_U_aerbulkwfs(1,w,q) + scene_weight(i_scene) * Stokes_U_aerbulkwfs_ic(1,w,q)
               endif
            enddo
         endif

!  fluxes

         if ( do_SphericalAlbedo ) then
            ACTINIC(1:n_szas,1:nstokes,w) = ACTINIC(1:n_szas,1:nstokes,w) + &
                                                            scene_weight(i_scene) * ACTINIC_ic(1:n_szas,1:nstokes,w)
            REGFLUX(1:n_szas,1:nstokes,w) = REGFLUX(1:n_szas,1:nstokes,w) + &
                                                            scene_weight(i_scene) * REGFLUX_ic(1:n_szas,1:nstokes,w)
         endif

!  Finish wavelength/wavenumber loop

         ! !write total column gas optical depth
         ! write(340,567)w,wav_shift_WN(ww),sum(LAYER_GASOPDEPS(1:NLAYERS,1)),sum(LAYER_GASOPDEPS(1:NLAYERS,2)),&
         !                       sum(LAYER_GASOPDEPS(1:NLAYERS,3)),&
         !                       sum(LAYER_GASOPDEPS(1:NLAYERS,4))



      
         ! pause
      
      enddo !DO w = 1, nwav !loop for wavelength


if (do_LRRS_cal) then      
! for LRRS mono


! !  Set up wavnumber for shift (all points at once.)
   range = ((wav(nwav)+200.00) - (wav(1)-200.00)) ! to cover Raman scattering range (-200 cm-1 to 200 cm-1)
   ! step = range / GEMSTOOL_INPUTS%Wavenums%wavnum_res ;---trying to change this res.
   interval = 0.01d0
   step = range / interval
   nwav_ext = INT(step) + 1
   do n = 1, nwav_ext
      ! wav_ext(n) = wav(1)-200.00 + GEMSTOOL_INPUTS%Wavenums%wavnum_res * real ( n - 1, fpk )
      wav_ext(n) = wav(1)-200.00 + interval * real ( n - 1, fpk ) !0.01d0 * real ( n - 1, fpk )
   enddo



!  N2 Spectroscopy Transitions (wavenumber shifts)

!  setting the Linearized VLIDORT control input here
!  -------------------------------------------------

!  Initialize - not used

!  10/17/16. These 2 variables no longer in Version 2.7
!   VLIDORT_LinFixIn%Cont%TS_DO_SURFBB_LINEARIZATION = .false.
!   VLIDORT_LinFixIn%Cont%TS_DO_LTE_LINEARIZATION    = .false.

VLIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES     = ' '
VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES    = ' '

VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS     = .false.
VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS       = 0

!  Initialize - main control variables

!   VLIDORT_LinFixIn%Cont%TS_DO_SIMULATION_ONLY        =  .false. !  Moved 10/17/16
VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY        =  .false.
VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  =  .false.
VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION =  .false.
VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   =  .false.
VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION =  .false.
VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         =  .false.

! 10/17/16. 2 new variables in V2p7, not required here....

VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF   =  .false.
VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF =  .false.

VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS         = 0
VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = 0
VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 0

VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG   = .false.
VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER = 0

!  Initialize - supplemental BRDF variables

VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC = ZERO
VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0        = ZERO
VLIDORT_LinSup%BRDF%TS_LS_BRDF_F          = ZERO
VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0   = ZERO
VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F     = ZERO

VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY      = ZERO
VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY = ZERO

!mick fix 1/24/2017 - added linearized SLEAVE input initializations
!  Initialize - supplemental SLEAVE variables

VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC  = ZERO
VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES = ZERO
VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0        = ZERO
VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0   = ZERO

!  Profile Jacobian flag

if ( do_GasProfile_Jacobians .or. do_AerOpdepProfile_Jacobians ) then
   VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION =  .true.
   VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   =  .true.
endif

!  COLUMN Jacobians, 15 august 2013. Tshift and/or Surface-pressure
!  @@ Rob fix 7/23/14, Update for H2O scaling flag and value

Q = 0
!   if ( do_Tshift_Jacobian )     Q = Q + 1
!   if ( do_Surfpress_Jacobian )  Q = Q + 1
!   if ( do_H2OScaling_Jacobian ) Q = Q + 1

if ( do_H2OScaling_Jacobian ) Q = Q + 1
if ( do_CH4Scaling_Jacobian ) Q = Q + 1
if ( do_Tshift_Jacobian )     Q = Q + 1
if ( do_Surfpress_Jacobian )  Q = Q + 1

!  BULK (COLUMN) AEROSOL Jacobians, New section June 2014
!    10/21/16. PARSLIST Added, Exception handling introduced.
!mick fix 1/24/2017 - added Do_LinearRef, Do_LinearPSD, Do_LinearEps, nlin(2), & npsd(2) to output

if ( do_AerBulk_Jacobians ) then
   call GEMSTOOL_AERWFS_BOOKKEEP &
      ( maxaerwfs, GEMSTOOL_INPUTS,                & ! Inputs
      n_AerBulk_Total_wfs, n_AerBulk_Group_wfs,  & ! OUTPUT, aerosol Jacobian Bookkeeping
      AerBulk_Mask_wfs, AerBulk_pars,            & ! OUTPUT, aerosol Jacobian Bookkeeping
      do_LinearRef, do_LinearPSD, Do_LinearEps,  & ! OUTPUT, aerosol Jacobian Bookkeeping
      NLIN, NPSD, Parslist,                      & ! OUTPUT, aerosol Jacobian Bookkeeping
      Fail_Bookkeep, Message_Bookkeep )            ! Exception handling
   if ( Fail_Bookkeep ) then
      do m = 1, 3
         Messages_Main(m+nmessages) = trim(Message_Bookkeep(m))
      enddo
      nmessages = nmessages + 3
      Fail_Main = .true. ; go to 4555
   endif
   Q = Q + n_AerBulk_Total_wfs
endif

!  debug (OK)
!   write(*,*) 'aa', n_AerBulk_Total_wfs
!   write(*,*) 'aa2=',n_AerBulk_Group_wfs
!   do n = 1, 16
!      write(*,*),'bb-', AerBulk_Mask_wfs(n),AerBulk_pars(n)
!   enddo
!   pause'book'

!  Check on Number of weighting functions

if ( do_AerBulk_Jacobians .and. Q.gt.MAX_ATMOSWFS ) then
   fail_Main = .true.
   messages_main(nmessages+1) = 'Number of Bulk Property Jacobians > Dimension MAXW_ATMOSWFS ==> Increase VLIDORT parameter'
   nmessages = nmessages + 1
   go to 4555
endif

!  Final Column property tally

if ( Q.gt.0 ) then
   VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  =  .true.
   VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS         = Q
endif

!  Surface jacobians - Lambertian albedo only

if ( do_Surface_Jacobians ) then
   VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .true.
   VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 1
endif


!  If no linearization, then set the "simulation-only" flag
!  @@ Rob fix 7/23/14, Update for H2O scaling flag and value

VLIDORT_LinModIn%MCont%TS_do_atmos_linearization   =  do_GasProfile_Jacobians.or.do_AerOpdepProfile_Jacobians.or.&
                                                      do_Tshift_Jacobian.or.do_Surfpress_Jacobian            .or.&
                                                      do_H2OScaling_Jacobian.or.do_CH4Scaling_Jacobian       .or.&
                                                      do_AerBulk_Jacobians
VLIDORT_LinModIn%MCont%TS_do_linearization         =  VLIDORT_LinModIn%MCont%TS_do_atmos_linearization .or. &
                                                      do_Surface_Jacobians

VLIDORT_LinModIn%MCont%TS_do_simulation_only        = .not.VLIDORT_LinModIn%MCont%TS_do_linearization
!   VLIDORT_LinFixIn%Cont%TS_do_simulation_only        = .not.VLIDORT_LinModIn%MCont%TS_do_linearization  ! Moved 10/17/16

!  Rob, 10/18/16, set VLIDORT linearization control for Flourescence Jacobians.
!               -->  These are pre-initialized Above.
!       11/30/16, Upgraded to include linear parameterization (2 Jacobians, not 1)

n_SIFPars = 0 ; SIFPars = 0.0d0
if ( GEMSTOOL_INPUTS%Sleave%DO_SIF ) then
   if ( GEMSTOOL_INPUTS%LinControl%do_SIF_Jacobians ) then
      VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS  = .true.
      if ( GEMSTOOL_INPUTS%Sleave%DO_SIF_LinearParam ) then
         SIFPars(1) = GEMSTOOL_INPUTS%Sleave%SIF755_Scaling_Constant
         SIFPars(2) = GEMSTOOL_INPUTS%Sleave%SIF755_Scaling_Gradient
         n_SIFPars  = 2
      else
         SIFPars(1) = GEMSTOOL_INPUTS%Sleave%SIF755_Amplitude
         n_SIFpars  = 1
      endif
      VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS = n_SIFPars
   endif
endif

!  Bookkeeping to count the number of Profile weighting functions
!  Dont forget the aerosol profile !!!!!

q = 0
if (do_GasProfile_Jacobians ) then
   do g = 1, ngases
      if ( do_gases(g) .and. do_gas_wfs(g) ) then
         q = q + 1
      endif
   enddo
endif
if (do_AerOpdepProfile_Jacobians ) q = q + 1
VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = q

!  Set up wavenumbers. Check Number of wavenumbers does not exceed MAXWAVNUMS

! range = GEMSTOOL_INPUTS%Wavenums%wavnum_finish - GEMSTOOL_INPUTS%Wavenums%wavnum_start 
! step  = range / GEMSTOOL_INPUTS%Wavenums%wavnum_res
! nwav  = INT(step) + 1
! if ( nwav .gt. maxwavnums) then
!    fail_Main = .true.
!    messages_main(nmessages+1) = 'Number of wavenumbers > Dimension MAXWAVNUMS ==> Increase MAXWAVNUMS'
!    nmessages = nmessages + 1
!    go to 4555
! endif

!  =============================================================================
!                   ATMOSPHERIC and SURFACE PROPERTY CREATION
!  =============================================================================

!  1A. get GAS and PTH profiles
!  *****************************

!  12 August 2013. Use NSW-specific Profiles here

!  Assumes that the GRIDDING is fine enough for aerosols

!  @@ Rob fix 7/23/14, Update for H2O scaling flag and value

call GEMSTOOL_NSW_PTHGAS_PROFILES &
( Maxlayers, MaxGases, ngases, which_gases, do_H2OScaling, do_CH4Scaling,   &  ! Input (control)
  gh2o, H2OScaling, gch4, CH4Scaling, Tshift, FDGas, FDLev, FDsfp, FDeps,   &  ! Input (numbers/Fd test)
  PhysicsDataPath, PTH_profile_name, GAS_profile_names,                     &  ! Input
  nlevels, level_heights, level_temps, level_press, level_vmrs,             &  ! Output (LEVELS)
  nlayers, aircolumns, daircolumns_dS, daircolumn_dP, Tshape,               &  ! Output (LAYERS). Air Density
  gascolumns, dgascolumns_dS, dgascolumn_dP, dgascolumns_dv,                &  ! Output (LAYERS). Gas densities
  dh2ocolumn_df, dch4column_df,                                             &  ! Output (LAYERS). Gas densities
  fail_Pthg, message_pthg )                                                    ! output

!  Surface pressure

  Surfpress = level_press(nlayers)
!   write(*,*)'Surfpress  = ',  Surfpress
  if ( fail_Pthg ) then
   fail_Main = .true.
   messages_main(nmessages+1) = adjustl(trim(message_pthg))
   nmessages = nmessages + 1
   go to 4555
  endif

!  set the cloud fraction proxy

!mick fix 6/18/2013 - added if condition
  if ( do_clouds ) cfrac = GEMSTOOL_INPUTS%Clouds%Cloud_fraction

  !  Set the Cloud-top heights and pressures
  
     if ( do_clouds .and. GEMSTOOL_INPUTS%Clouds%do_Scattering_Clouds ) then
        if ( GEMSTOOL_INPUTS%Clouds%Cloud_boundaries ) then
           GEMSTOOL_INPUTS%Clouds%ctop_height   = level_heights(GEMSTOOL_INPUTS%Clouds%ctop_level)
           GEMSTOOL_INPUTS%Clouds%cbot_height   = level_heights(GEMSTOOL_INPUTS%Clouds%cbot_level)
           GEMSTOOL_INPUTS%Clouds%ctop_pressure = level_press  (GEMSTOOL_INPUTS%Clouds%ctop_level)
           GEMSTOOL_INPUTS%Clouds%cbot_pressure = level_press  (GEMSTOOL_INPUTS%Clouds%cbot_level)
        endif
     endif
  
  !  set the profile linearization control
  
     if ( VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION ) then
        do n = 1, nlayers
           VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(n)   = .true.
           VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(n) = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
        enddo
     endif
  
  !  Set the Column Linearization control
  
     if ( VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION ) then
        do n = 1, nlayers
           VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(n)   = .true.
           VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(n) = VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS
        enddo
     endif
     
!  1B. Get Gas and Rayleigh Cross-sections. Strictly level output for Gas Cross-sections.
!  ****************************************

!  Special for GEMS:
!     -- ALL CROSS_SECTIONS are on LEVELS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     -- O2/H2O/CO2/CH4 Cross-sections will be P-T dependent, derived from HITRAN

!  21 October 2013. develop handles for the LUT alternative.

     if ( do_Tshift_Jacobian .or.do_Surfpress_Jacobian ) then

      if ( do_Hitran ) then
         call GEMSTOOL_NSW_Cross_SECTIONS_PLUS &
          ( maxwavnums, maxlayers, maxgases, CO2_PPMV_MIXRATIO,     & ! INPUT
            nlayers, nwav_ext, level_temps, level_press, wav_ext,  & ! INPUT
            ngases, which_gases, hitran_path,                       & ! Input   !! Y.Jung
            do_Tshift_Jacobian, do_Surfpress_Jacobian,              & ! Input
            do_gases,                                               & ! In/Out  !! Y.Jung
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL,                          & ! OUTPUT
            GAS_XSECS, dGASXSECS_dT, dGASXSECS_dP,                  & ! OUTPUT
            fail_xsec, message_xsec )                                 ! output
      else
         call GEMSTOOL_NSW_Cross_SECTIONS_PLUS_LUT &
          ( maxwavnums, maxlayers, maxgases, CO2_PPMV_MIXRATIO,     & ! INPUT
            nlayers, nwav_ext, level_temps, level_press, wav_ext,  & ! INPUT
            ngases, which_gases, lut_path,                          & ! Input   !! Y.Jung
            do_Tshift_Jacobian, do_Surfpress_Jacobian,              & ! Input
            do_gases,                                               & ! In/Out  !! Y.Jung
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL,                          & ! OUTPUT
            GAS_XSECS, dGASXSECS_dT, dGASXSECS_dP,                  & ! OUTPUT
            fail_xsec, message_xsec )                                 ! output
      endif

      !If necessary, reduce the number of gas derivatives to be done
      !if a gas is not present in the spectral range
      do g = 1, ngases
         if ( do_GasProfile_Jacobians .and. .not.do_gases(g) &
              .and. do_gas_wfs(g) ) then
            do_gas_wfs(g) = .false.
            if ( VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION ) then
              VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = &
                VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS - 1
              VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:nlayers) = &
                VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
            endif
         end if
      end do
      else

     if ( do_Hitran ) then
        call GEMSTOOL_NSW_Cross_SECTIONS &
          ( maxwavnums, maxlayers, maxgases, CO2_PPMV_MIXRATIO,   & ! INPUT
            nlayers, nwav_ext, level_temps, level_press, wav_ext,         & ! INPUT
            ngases, which_gases, hitran_path,                     & ! Input     !! Y.Jung
            do_gases,                                             & ! In/Out    !! Y.Jung
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL, GAS_XSECS,             & ! OUTPUT
            fail_xsec, message_xsec )                               ! output
            ! write(*,*)do_gases
            ! write(*,*) GAS_XSECS(1:nwav_shift,1,2)
            ! pause
      else
        call GEMSTOOL_NSW_Cross_SECTIONS_LUT &
          ( maxwavnums, maxlayers, maxgases, CO2_PPMV_MIXRATIO,   & ! INPUT
            nlayers, nwav_ext, level_temps, level_press, wav_ext,         & ! INPUT
            ngases, which_gases, lut_path,                        & ! Input     !! Y.Jung
            do_gases,                                             & ! In/Out    !! Y.Jung
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL, GAS_XSECS,             & ! OUTPUT
            fail_xsec, message_xsec )                               ! output
      endif

   endif

!  Multiply Tshift derivatives by Shape factor

   if ( do_Tshift_Jacobian ) then
      do ww = 1, nwav_ext
         do n = 0, nlayers
            dGASXSECS_dT(ww,n,1:ngases) = dGASXSECS_dT(ww,n,1:ngases) * Tshape(n)
         enddo
      enddo
   endif

   if ( do_Surfpress_Jacobian ) then
      do ww = 1, nwav_ext
         dGASXSECS_dP(ww,nlayers,1:ngases) = dGASXSECS_dP(ww,nlayers,1:ngases) * Surfpress!1013.25d0
      enddo
   endif

!  Debug.
!    Print out Rayleigh Xsec + depol, 4 gas sections for the lowest level.

!   if ( XSEC_debug_flag ) then
!      open(57,file='Debug_CrossSections.out',status = 'replace')
!      do ww = 1, nwav_ext
!         write(57,*)ww,wav_ext(ww),RAYLEIGH_XSEC(ww), RAYLEIGH_DEPOL(ww)!, (GAS_XSECS(w,nlayers,g),g=1,4)
!      enddo
!      close(57) 
! !   pause'Debug Cross sections'
!   endif
! pause

!  Exception handling

   if ( fail_Xsec ) then
      fail_Main = .true.
      messages_main(nmessages+1) = adjustl(trim(message_xsec))
      nmessages = nmessages + 1
      go to 4555
   endif


!  2. get aerosol and cloud properties
!  ***********************************

!  if aerosols not present, zero the aerosol data

   if ( .not. do_aerosols ) then
      aerosol_layerflags = .false.
      Loading_aerosol    = 0.0d0
      aerosol_nscatmoms  = 0
      aertau_unscaled    = 0.0d0
      aod_scaling        = 0.0d0
      aerosol_deltau     = 0.0d0
      aerosol_ssalbs     = 0.0d0
      aerosol_scatmoms   = 0.0d0
      aerosol_distchars  = 0.0d0
      go to 9654
   endif

   if ( do_AerBulk_Jacobians ) then
      call GEMSTOOL_AER_PROPERTIES_PLUS &
      ( maxlayers, maxwavnums, maxaermoms, maxaerwfs, interpolate_aerosols, do_wavnums, & ! Dimensions
        nlayers, nmuller, nwav_ext, wav_ext, level_heights, GEMSTOOL_INPUTS, momsize_cutoff,    & ! Inputs
        n_AerBulk_Total_wfs, n_AerBulk_Group_wfs, AerBulk_pars,                         & ! Inputs, aerosol Jacobian Bookkeeping
        do_LinearRef, do_LinearPSD, Do_LinearEps, NLIN, NPSD, ParsList,                 & ! Inputs, aerosol Jacobian Bookkeeping
        aerosol_layerflags, Loading_aerosol, Dloading_Dtau, Dloading_Dpars,             & ! OUTPUT, aerosol Loading
        aerosol_nscatmoms, aertau_unscaled, aod_scaling, aerosol_distchars,             & ! OUTPUT, aerosol optical properties
        aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms,                               & ! OUTPUT, aerosol optical properties
        L_aertau_unscaled, L_aod_scaling, L_aerosol_deltau,                             & ! OUTPUT, aerosol Linearized OPs
        L_aerosol_ssalbs, L_aerosol_scatmoms,                                           & ! OUTPUT, aerosol Linearized OPs
        fail1_Aer, fail2_Aer, Message_AerLoading, Message_AerOptical )                    ! Exception handling
   else
      call GEMSTOOL_AER_PROPERTIES &
      ( maxlayers, maxwavnums, maxaermoms,                                                & ! Dimensions
        interpolate_aerosols, do_wavnums, FDAer, FDLay, FDBul, FDeps,                     & ! Flags and FD control
        nlayers, nmuller, nwav_ext, wav_ext, level_heights, GEMSTOOL_INPUTS, momsize_cutoff,      & ! Inputs
        aerosol_layerflags, Loading_aerosol, aerosol_nscatmoms, aertau_unscaled,          & ! OUTPUT, aerosol control
        aod_scaling, aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms, aerosol_distchars, & ! OUTPUT, aerosol optical
        fail1_Aer, fail2_Aer, Message_AerLoading, Message_AerOptical )                      ! Exception handling
   endif

!  Exception handling

   if ( fail1_Aer ) then
      messages_Main(nmessages+1:nmessages+3) = message_AerLoading(1:3)
      nmessages = nmessages + 4
      messages_Main(nmessages) = 'Aerosol Loading failed (regridding first)- look at preceding 3 messages'
      go to 4555
   endif
   if ( fail2_Aer ) then
      messages_Main(nmessages+1:nmessages+5) = message_AerOptical(1:5)
      nmessages = nmessages + 6
      if ( do_AerBulk_Jacobians ) then
         messages_Main(nmessages) = 'Aerosol-Linearized Mie/Tmatrix failed (regridding first) - look at preceding 5 messages'
      else
         messages_Main(nmessages) = 'Aerosol-Regular Mie/Tmatrix failed (regridding first) - look at preceding 5 messages'
      endif
      go to 4555
   endif

!  Dimensioning check on number of moments
!    10/19/16, take into account wave number dependence

   if ( MAXVAL(aerosol_nscatmoms(1:nwav_ext)) .gt. MAXMOMENTS_INPUT ) then
      fail_Main = .true.
      messages_main(nmessages+1) = 'Number of Aerosol Expansion Coefficients > MAXMOMENTS_INPUT'
      nmessages = nmessages + 1
      messages_main(nmessages+1) = 'Action: Increase MAXMOMENTS_INPUT, preferably to the same as maxaermoms'
      nmessages = nmessages + 1
      go to 4555
   endif

!  Continuation point for skipping aerosols

9654 continue

   !  if clouds not present, zero the cloud data and move on.
   
      if ( .not. do_clouds ) then
          cloud_layerflags = .false.
          Loading_cloud    = 0.0d0
          cldtau_unscaled  = 0.0d0
          cod_scaling      = 0.0d0      
          cloud_nscatmoms  = 0
          cloud_deltau     = 0.0d0
          cloud_ssalbs     = 0.0d0
          cloud_scatmoms   = 0.0d0
          go to 9664
       endif
   
   !  Generation of Cloud optical properties
   !     1. Call to the Regular Mie program
   !     2. Convert Mie output (Microsopic) to IOP output (macroscopic)
   
      write(*,*) 'Doing cloud'
   
      call GEMSTOOL_CLD_PROPERTIES &
         ( maxlayers, maxwavnums, maxcldmoms, interpolate_clouds, do_wavnums,             & ! Dimensions and Flags 
           nlayers, nmuller, nwav_ext, wav_ext, level_heights, GEMSTOOL_INPUTS, momsize_cutoff,   & ! Inputs3x
           cloud_layerflags, Loading_cloud, cloud_nscatmoms, cldtau_unscaled,             & ! OUTPUT, aerosol control
           cod_scaling, cloud_deltau, cloud_ssalbs, cloud_scatmoms,                       & ! OUTPUT, aerosol optical
           fail1_Cld, Message_CLdOptical )
   
   !  Exception handling
   
      if ( fail1_Cld ) then
         messages_Main(nmessages+1:nmessages+4) = message_CldOptical(1:4)
         nmessages = nmessages + 5
         messages_Main(nmessages) = 'Cloud Mie calculation failed - look at preceding 5 messages'
         go to 4555
      endif
   
   !  Dimensioning check on number of moments
   
      if ( cloud_nscatmoms .gt. MAXMOMENTS_INPUT ) then
         fail_Main = .true.
         messages_main(nmessages+1) = 'Number of Cloud Expansion Coefficients > MAXMOMENTS_INPUT'
         nmessages = nmessages + 1
         messages_main(nmessages+1) = 'Action: Increase MAXMOMENTS_INPUT, preferably to the same as maxcldmoms'
         nmessages = nmessages + 1
         go to 4555
      endif
   
   !  Continuation point for skipping clouds
   
9664 continue

!  3. Prepare surface property
!  ***************************

!  Zero the Lambertian albedo

lambertian_albedo = 0.0d0  

!  Zero the Supplement results
    
!   GEMSTOOL_BRDF_Results   = 0.0d0
   GEMSTOOL_SLEAVE_Results%SL_SLTERM_ISOTROPIC  = 0.0d0
   GEMSTOOL_SLEAVE_Results%SL_SLTERM_USERANGLES = 0.0d0
   GEMSTOOL_SLEAVE_Results%SL_SLTERM_F_0        = 0.0d0
   GEMSTOOL_SLEAVE_Results%SL_USER_SLTERM_F_0   = 0.0d0

!  BRDF for Land Surfaces
!  ----------------------

!  10/25/16. Initial Version for Surfaces that are NOT OCEAN
!    -- Same BRDF for all wavelengths, done outside the main wavelength loop.
!    -- Up to 3 error messages allowed.
   do_BRDF = GEMSTOOL_INPUTS%BRDF%do_Gemstool_BRDF

   if ( GEMSTOOL_INPUTS%BRDF%do_Gemstool_BRDF ) then

      CALL GEMSTOOL_BRDFSurface &
      ( GEMSTOOL_INPUTS, BRDF_N_SURFACE_WFS, BRDF_K_FP_values, &
        GEMSTOOL_BRDF_Results,  GEMSTOOL_BRDF_REsults_Lin, & ! Input and output ! editted by mchoi
        fail_BRDF, nmessages_BRDF, messages_BRDF ) ! Exception handling

      if ( fail_BRDF ) then
         fail_Main = .true.
         do m = 1, nmessages_brdf
            messages_main(nmessages+m) = Trim(messages_BRDF(m))
         enddo
         nmessages = nmessages + nmessages_brdf
         messages_main(nmessages+1) = 'GEMSTOOL Land-surface BRDF calculation failed - look at preceding messages'
         nmessages = nmessages + 1
         go to 4555
      endif
      VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS  = BRDF_N_SURFACE_WFS
   endif
  
   ! N_SURFACE_WFS_fin = VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS

!  Lambertian albedo
!  -----------------

   if ( .not. GEMSTOOL_INPUTS%BRDF%do_Gemstool_BRDF ) then

     call GEMSTOOL_Albedo &
        ( Maxwavnums, MaxAlbCoeffs, nwav_ext, wav_ext, do_wavnums, GEMSTOOL_INPUTS, & ! inputs
          PhysicsDataPath, alb_filename,                                    & ! inputs
          lambertian_albedo, ClosureTerms )                                   ! output

!  Debug
!  do w = 1, nwav
!     write(*,*)w,wav_shift_WN(ww),lambertian_albedo(w)
!  enddo
!  pause

   endif
  
!  ===========================================================================
!          CALCULATION SECTION : MAIN WAVELENGTH LOOP WITH VLIDORT
!  ===========================================================================

!  TOA or BOA output level

   if ( VLIDORT_FixIn%Bool%TS_do_upwelling ) then
      VLIDORT_ModIn%MUserVal%TS_user_levels(1) = GEMSTOOL_INPUTS%RTControl%observation_height !zero
      !VLIDORT_ModIn%MUserVal%TS_user_levels(2) = 3 !GEMSTOOL_INPUTS%RTControl%observation_height !zero
   else if ( VLIDORT_FixIn%Bool%TS_do_dnwelling ) then
      VLIDORT_ModIn%MUserVal%TS_user_levels(1) = dble(nlayers)
   endif
 
 !  Set up IPA. Number of Scenes.
 
    local_cloud(1:2)  = .false.
    scene_weight(1:2) = 1.0_fpk
 !mick fix 6/18/2014 - switched to permanent variable "do_clouds" now
    !if ( .not. do_clouds_V2 ) then
    if ( .not. do_clouds ) then
       n_scenes = 1   ! Clear sky calculation, no clouds
    else
       if ( cfrac .eq. 0.0_fpk ) then
          n_scenes = 1   ! Clear sky calculation, cloud fraction = 0
       else if ( cfrac .eq. 1.0_fpk ) then
          n_scenes = 1 ;  local_cloud(1) = .true. ! Fully cloudy calculation, cloudfraction = 1
       else
          n_scenes = 2 ;  local_cloud(2) = .true.  ! Partially cloudy
          scene_weight(1) = 1.0_fpk - cfrac
          scene_weight(2) = cfrac
       endif
    endif
 
! !  !  Zero the Final output
 
!     STOKES  = 0.0_fpk
!     DOLP    = 0.0_fpk
!     ACTINIC = 0.0_fpk
!     REGFLUX = 0.0_fpk
 
!     RADIANCE_GASWFS      = 0.0_fpk
!     RADIANCE_AERPROFWFS  = 0.0_fpk
!     RADIANCE_AERBULKWFS  = 0.0_fpk ! New, June 2014
!     RADIANCE_ALBEDOWFS   = 0.0_fpk
!     RADIANCE_TSHIFTWFS   = 0.0_fpk
!     RADIANCE_SURFPWFS    = 0.0_fpk
!     RADIANCE_WSCALEWFS   = 0.0_fpk ! New 23 July 2014
!     RADIANCE_MSCALEWFS   = 0.0_fpk ! New 01 Feb  2015
!     RADIANCE_SIFPARAMWFS = 0.0_fpk ! New, 10/18/16, 11/30/16    ! SIF introduction
 
!     AMFSCATWTS          = 0.0_fpk
!     AMFTOTAL            = 0.0_fpk
 
!  !  Stokes Q/U Jacobians, initialization added 12/19/16 by Rob
 
!     STOKES_Q_GASWFS      = 0.0_fpk ; STOKES_U_GASWFS      = 0.0_fpk
!     STOKES_Q_ALBEDOWFS   = 0.0_fpk ; STOKES_U_ALBEDOWFS   = 0.0_fpk
!     STOKES_Q_SIFPARAMWFS = 0.0_fpk ; STOKES_U_SIFPARAMWFS = 0.0_fpk
!     STOKES_Q_TSHIFTWFS   = 0.0_fpk ; STOKES_U_TSHIFTWFS   = 0.0_fpk
!     STOKES_Q_SURFPWFS    = 0.0_fpk ; STOKES_U_SURFPWFS    = 0.0_fpk
!     STOKES_Q_WSCALEWFS   = 0.0_fpk ; STOKES_U_WSCALEWFS   = 0.0_fpk
!     STOKES_Q_MSCALEWFS   = 0.0_fpk ; STOKES_U_MSCALEWFS   = 0.0_fpk
!     STOKES_Q_AERPROFWFS  = 0.0_fpk ; STOKES_U_AERPROFWFS  = 0.0_fpk
!     STOKES_Q_AERBULKWFS  = 0.0_fpk ; STOKES_U_AERBULKWFS  = 0.0_fpk
 
    dir = UpDn_Index
 
 !  Start the CLoud Scene loop. 1 or 2 Scenes (Clear-sky, Cloudy)
    
    n_scenes_lrrs = n_scenes
    do i_scene_lrrs = 1, n_scenes_lrrs
            
      !mick mod 1/24/2017 - added more descriptive feedback to user
            ! write(*,*)
            ! if ( .not.do_clouds ) then
            !    write(*,*)'Doing clear scene'
            ! else
            !    if ( cfrac .eq. 0.0_fpk ) then
            !       write(*,*)'Doing clear scene'
            !    else if ( cfrac .eq. 1.0_fpk ) then
            !       write(*,*)'Doing cloudy scene'
            !    else
            !       if (i_scene_lrrs .eq. 1) then
            !          write(*,*)'Doing partly cloudy scene: clear portion'
            !       else
            !          write(*,*)'Doing partly cloudy scene: cloudy portion'
            !       endif
            !    endif
            ! endif

!  Start MAIN WAVELENGTH/WAVENUMBER LOOP WITH VLIDORT - LRRS
!  --------------------------------------------------
      ! open(340,file='./GEMSTOOL_NSW_Results/GASOPDEPS.dat',status='replace')

!      nwav = 1

      
   DO ww = 1, nwav_ext

   !  progress
   
            if (ww == 1) write(*,*)
            !if (mod(w,10000).eq.0) write(*,*)' Linearized VLIDORT Calculation,  doing wavelength/wavenumber # ',w
            ! if (mod(ww,1).eq.0) write(*,*)' Linearized VLIDORT Calculation,  doing wavelength/wavenumber # ',w
            if (mod(ww,1000).eq.0) write(*,*)' LRRS shifting,  doing wavelength/wavenumber # ',ww
   
   !  Fluorescence 
   !  ------------
   
            Polynomial = 1.0d0 ; dPolynomial_dSL = 0.0d0 ; SIF_ExactOnly = .false.   ! Initialized for safety
   
   !  New code to implement SIF, 10/18/16. tested 10/21/16. Updated 10/26/15
   
            if ( GEMSTOOL_INPUTS%Sleave%do_SIF ) then
   
   !  ** Upgrade 11/30/16. Only one call if you are using the linear parameterization
   !                       In this case you are just getting SIF 755 value.
   
               if ( ( GEMSTOOL_INPUTS%Sleave%DO_SIF_LinearParam.and.ww.eq.1) .or. &
                     .not. GEMSTOOL_INPUTS%Sleave%DO_SIF_LinearParam ) then
                  Call GEMSTOOL_Fluorescence &
                     ( do_wavnums, wav_ext(ww), GEMSTOOL_INPUTS, SIF_datapath, &
                     GEMSTOOL_SLEAVE_Results )
               endif
   
   !   Upgrade 11/30/16. calculate polynomial for the SIF linear parameterization
   
               SIF_ExactOnly = GEMSTOOL_INPUTS%Sleave%DO_SIF_ExactOnly
               if ( GEMSTOOL_INPUTS%Sleave%DO_SIF_LinearParam ) then
                  factor =  wav_ext(ww) ; if ( do_wavnums ) factor = 1.0d+07/wav_ext(ww)
                  factor = 1.0d0 - ( factor / 755.0d0 )
                  Polynomial = SIFPars(1) + SIFPars(2) * Factor
                  dPolynomial_dSL(1) = 1.0d0
                  dPolynomial_dSL(2) = Factor
               else
                  Polynomial = 1.0d0
                  dPolynomial_dSL(1) = 1.0d0 / SIFpars(1)
               endif
   
            endif
   
   !  Develop the VLIDORT-specific optical properties
   !  -----------------------------------------------
   
   !  @@ Rob fix 7/23/14, add H2O scaling derivative
   
            if ( (do_GasProfile_Jacobians.or.do_AerOpdepProfile_Jacobians) .or. &
                  (do_Tshift_Jacobian.or.do_Surfpress_Jacobian.or.do_H2OScaling_Jacobian.or.do_CH4Scaling_Jacobian .or. &
                  do_AerBulk_Jacobians) ) then
               call GEMSTOOL_NSW_linearized_iops &
            ( maxlayers, maxmoments_input, maxstokes_sq, max_atmoswfs,                 & ! VLIDORT dimensions
               maxwavnums, maxaermoms, maxcldmoms, maxgases, maxaerwfs,                 & ! GEMSTOOL dimensions
               IOP_debug_filename, IOP_Debug_flag,                                      & ! Debug IOP check  
               do_GasProfile_Jacobians, do_AerOpdepProfile_Jacobians,                   & ! Profile Jacobian control flags
               do_AerBulk_Jacobians, do_Tshift_Jacobian,                                & ! Column/Bulk Jacobian control flags
               do_Surfpress_Jacobian, do_H2OScaling_Jacobian, do_CH4Scaling_Jacobian,   & ! Column/Bulk Jacobian control flags
               ww, wav_ext(ww), omega_lim, do_gases, do_gas_wfs, ngases,           & ! GEMSTOOL Inputs (Control)
               nstokes, nlayers, nwav_ext, n_AerBulk_Total_wfs, AerBulk_pars,       & ! GEMSTOOL Inputs (Control)
               do_aerosols, aerosol_nscatmoms, local_cloud(i_scene_lrrs), cloud_nscatmoms,   & ! GEMSTOOL Inputs (Aer/clouds)
               gascolumns, dGASCOLUMNS_dS, dGASCOLUMN_dP, dGASCOLUMNS_dV,               & ! GEMSTOOL Inputs (GasColumnsderivs)
               dH2OCOLUMN_dF, dCH4COLUMN_dF,                                            & ! GEMSTOOL Inputs (GasColumnsderivs)
               aircolumns, dAIRCOLUMNS_dS, dAIRCOLUMN_dP,                               & ! GEMSTOOL Inputs (AirColumns + derivs)
               Tshift, SurfPress, H2OScaling, gh2o, CH4Scaling, gch4,                   & ! GEMSTOOL Inputs (T-shift,SurfP,Scaling)
               gas_xsecs, dgasxsecs_dT, dgasxsecs_dP, Rayleigh_xsec, Rayleigh_depol,    & ! GEMSTOOL Inputs (Cross-sections, depol)
               aerosol_layerflags, aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms,    & ! GEMSTOOL Inputs (Aerosol)
               L_aerosol_deltau, L_aerosol_ssalbs, L_aerosol_scatmoms,                  & ! GEMSTOOL Inputs (Aerosol linearized)
               cloud_layerflags, cloud_deltau, cloud_ssalbs, cloud_scatmoms,            & ! GEMSTOOL Inputs (Clouds)
               DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT,                                    & ! VLIDORT Outputs (Regular iops)
               NGREEK_MOMENTS_INPUT, GREEKMAT_TOTAL_INPUT,                              & ! VLIDORT Outputs (Regular iops)
               L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,        & ! VLIDORT Outputs (linearized)
               LAYER_GASOPDEPS, dGASOPDEPS_dV )                                           ! Layer-to-Level transformation
            else
               call GEMSTOOL_NSW_regular_iops &
            ( maxlayers, maxmoments_input, maxstokes_sq,                            & ! VLIDORT dimensions
               maxwavnums, maxaermoms, maxcldmoms, maxgases,                         & ! GEMSTOOL dimensions
               IOP_debug_filename, IOP_debug_flag,                                   & ! Debug IOP check  
               ww, wav_ext(ww), omega_lim, do_gases, ngases, nstokes, nlayers, nwav_ext, & ! GEMSTOOL Inputs
               do_aerosols, aerosol_nscatmoms, local_cloud(i_scene_lrrs), cloud_nscatmoms,        & ! GEMSTOOL Inputs
               gas_xsecs,                                                            & ! GEMSTOOL Inputs
               gascolumns, aircolumns, Rayleigh_xsec, Rayleigh_depol,                & ! GEMSTOOL Inputs
               aerosol_layerflags, aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms, & ! GEMSTOOL Inputs
               cloud_layerflags, cloud_deltau, cloud_ssalbs, cloud_scatmoms,         & ! GEMSTOOL Inputs
               DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT,                                 & ! VLIDORT Outputs
               NGREEK_MOMENTS_INPUT, GREEKMAT_TOTAL_INPUT )                            ! VLIDORT Outputs
            endif
         
         ! save parameter for LRRS calculation      
         ! 1. GREEKMAT_TOTAL_INPUT -> GREEKMAT_TOTAL_INPUT_ext
         ! 2. OMEGA_TOTAL_INPUT -> OMEGA_TOTAL_INPUT_ext
         ! 3. DELTAU_VERT_INPUT -> DELTAU_VERT_INPUT_ext

            do ia=1, MAXLAYERS
            DELTAU_VERT_INPUT_ext(ia,ww)=DELTAU_VERT_INPUT(ia)
            OMEGA_TOTAL_INPUT_ext(ia,ww)=OMEGA_TOTAL_INPUT(ia)
               do ib=0, MAXMOMENTS_INPUT
               GREEKMAT_TOTAL_INPUT_ext(ib,ia,ww)= &
                                       GREEKMAT_TOTAL_INPUT(ib,ia,1) * OMEGA_TOTAL_INPUT(ia)
               ENDDO
            enddo
            
                                     ! write(*,*)DELTAU_VERT_INPUT_ext(1,ww),DELTAU_VERT_INPUT_ext(nlayers,ww)
            ! write(*,*)OMEGA_TOTAL_INPUT_ext(1,ww),OMEGA_TOTAL_INPUT_ext(nlayers,ww)
            ! ! write(*,*)GREEKMAT_TOTAL_INPUT_ext(:,1,ww)
            ! ! write(*,*)GREEKMAT_TOTAL_INPUT_ext(:,nlayers,ww)
            ! write(*,*)GREEKMAT_TOTAL_INPUT_ext(0:30,1,ww)
            ! write(*,*)GREEKMAT_TOTAL_INPUT_ext(0:30,20,ww)
            ! write(*,*)
            ! pause


         !write(*,*) ww, wav_shift_WN(ww), nwav_shift
         !write(*,*) ww, wav_shift_WN(ww), omega_lim, do_gases, ngases, nstokes, nlayers, nwav_shift
         !write(*,*) GASCOLUMNS(:,1,:) !GAS_XSECS(1,:,1)


         ! write(*,*)DELTAU_VERT_INPUT
         ! write(*,*)
         ! write(*,*)OMEGA_TOTAL_INPUT
         ! write(*,*)
         ! write(*,*)GREEKMAT_TOTAL_INPUT(:,2,1)

         ! pause

      ! write(*,*) NGREEK_MOMENTS_INPUT
      ! write(*,*) GREEKMAT_TOTAL_INPUT(0:30,20, 1)
      
   
      ! !put some results to the LRRS preparation  ---LRRS; is this here right?
      ! do i_layer = 1, maxlayers
      !    DELTAU_INPUT_UNSCALED(i_layer,ww) = DELTAU_VERT_INPUT(i_layer)
      !    do i_moment = 0, MAXMOMENTS_INPUT
      !       OMEGAMOMS_ELASTIC_UNSCALED(i_layer,i_moment,ww)=GREEKMAT_TOTAL_INPUT(i_moment,i_layer,1)
      !    enddo
      !    ! write(*,*) DELTAU_INPUT_UNSCALED(i_layer,w)
      !    ! write(*,*) OMEGAMOMS_ELASTIC_UNSCALED(i_layer, 0:20, w)
      ! enddo
   
   
   
   !  Do the RT Calculation (2 parts: Actual and Final)
   !  -------------------------------------------------
   
   !  @@ Rob   fix 7/23/14, add H2O scaling derivative
   !  @@ YJung fix 2/1/15 , add CH4 scaling derivative
   !  @@ Rob, 10/18/16, Introduce SIF arguments (prepared SIF data) --> replaced by SLEAVE structure
   !                  --> Include flag for using SIF Jacobians (renamed 11/30/16)
   !  @@ Rob, 10/25/16, Introduce BRDF and SLEAVE Type structures
   !  @@ Rob, 10/25/16, Renamed module GEMSTOOL_NSW_SetJacobians
   !  @@ Rob, 11/30/16, Introduce SIF inputs (for SIF linear parameterization). revised naming.
   !  @@ Rob, 12/19/16, Introduce Stokes-Q and Stokes-U output via subroutine GEMSTOOL_NSW_SetJacobians_IQU
   
   !          if ( ( do_GasProfile_Jacobians.or.do_AerOpdepProfile_Jacobians ) .or. &
   !                ( do_Tshift_Jacobian.or.do_Surfpress_Jacobian.or.do_H2OScaling_Jacobian.or.do_CH4Scaling_Jacobian .or. &
   !                   do_AerBulk_Jacobians ) .or. &
   !                ( do_Surface_Jacobians .or. do_SIF_Jacobians ) ) then
                         
   !             call GEMSTOOL_RTCALC_PLUS_Actual &
   !             (  do_GasProfile_Jacobians, do_AerOpdepProfile_Jacobians,                    & ! Main control
   !                GEMSTOOL_INPUTS%Sleave%do_SIF,Polynomial, dPolynomial_dSL, SIF_ExactOnly, & ! SIF control
   !                VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                            & ! VLIDORT regular inputs
   !                VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup,                   & ! VLIDORT linearized inputs
   !                nlayers, NGREEK_MOMENTS_INPUT, level_heights, lambertian_albedo(ww),   & ! Control/Surface Proxies
   !                DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,           & ! Atmos-Optical Proxies
   !                L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,     & ! Atmos-Optical Proxies (linearized)
   !                GEMSTOOL_SLEAVE_Results, GEMSTOOL_BRDF_Results,                       & ! Surface Supplement Results
   !                VLIDORT_Out, VLIDORT_LinOut )                                           ! VLIDORT Output results
                  
   !             call GEMSTOOL_NSW_SetJacobians_IQU &
   !          (  MAXWAVNUMS, MaxAlbCoeffs, maxgases, maxmessages, maxaerwfs, WW, dir, do_SphericalAlbedo,      & ! GEMSTOOL Dimensions and W
   !             do_normalized_wfoutput, do_GasProfile_Jacobians, do_AerOpdepProfile_Jacobians,               & ! GEMSTOOL Jacobian control
   !             do_Surface_Jacobians, do_SIF_Jacobians, do_AerBulk_Jacobians,                                & ! GEMSTOOL Jacobian control
   !             do_Tshift_Jacobian, do_Surfpress_Jacobian, do_H2OScaling_Jacobian, do_CH4Scaling_Jacobian,   & ! GEMSTOOL Jacobian control
   !             do_gases, do_gas_wfs, ngases, nlayers, nstokes, n_geometries, n_AerBulk_Total_wfs, n_SIFPars,& ! Control numbers/gases
   !             Lambertian_albedo, SIFPars, Tshift, SurfPress, H2OScaling, CH4Scaling,                       & ! necessary for normaliz
   !             Level_VMRs, AERTAU_UNSCALED, AerBulk_pars,                                                   & ! necessary for normaliz
   !             VLIDORT_Out, VLIDORT_LinOut,                                              & ! VLIDORT Results (includes linearized)
   !             LAYER_GASOPDEPS, dGASOPDEPS_dV,                                           & ! Layer-to-Level transformations
   !             STOKES_ic, ACTINIC_ic, REGFLUX_ic, DOLP_ic, DOCP_ic,                      & ! Main program, GEMSTOOL Stokes Output
   !             RADIANCE_GASWFS_ic,    RADIANCE_AERPROFWFS_ic,  RADIANCE_AERBULKWFS_ic,   & ! Main program, GEMSTOOL I-Jacobians Output
   !             RADIANCE_ALBEDOWFS_ic, RADIANCE_SIFPARAMWFS_ic, RADIANCE_TSHIFTWFS_ic,    & ! Main program, GEMSTOOL I-Jacobians Output
   !             RADIANCE_SURFPWFS_ic,  RADIANCE_WSCALEWFS_ic,   RADIANCE_MSCALEWFS_ic,    & ! Main program, GEMSTOOL I-Jacobians Output
   !             AMFSCATWTS_ic, AMFTOTAL_ic, TOTALWF_SAVED_ic,                             & ! Main program, GEMSTOOL AMF Output
   !             STOKES_Q_GASWFS_ic,    STOKES_Q_AERPROFWFS_ic,  STOKES_Q_AERBULKWFS_ic,   & ! Main program, GEMSTOOL Q-Jacobians Output
   !             STOKES_Q_ALBEDOWFS_ic, STOKES_Q_SIFPARAMWFS_ic, STOKES_Q_TSHIFTWFS_ic,    & ! Main program, GEMSTOOL Q-Jacobians Output
   !             STOKES_Q_SURFPWFS_ic,  STOKES_Q_WSCALEWFS_ic,   STOKES_Q_MSCALEWFS_ic,    & ! Main program, GEMSTOOL Q-Jacobians Output
   !             STOKES_U_GASWFS_ic,    STOKES_U_AERPROFWFS_ic,  STOKES_U_AERBULKWFS_ic,   & ! Main program, GEMSTOOL U-Jacobians Output
   !             STOKES_U_ALBEDOWFS_ic, STOKES_U_SIFPARAMWFS_ic, STOKES_U_TSHIFTWFS_ic,    & ! Main program, GEMSTOOL U-Jacobians Output
   !             STOKES_U_SURFPWFS_ic,  STOKES_U_WSCALEWFS_ic,   STOKES_U_MSCALEWFS_ic,    & ! Main program, GEMSTOOL U-Jacobians Output
   !             Errorstatus, nmessages, messages_Main )                                     ! Main program, exception handling
   
   !          else
               
   !             call GEMSTOOL_RTCALC_Actual &
   !             ( VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                          & ! VLIDORT Actual inputs
   !                Polynomial, SIF_ExactOnly,                                          & ! SIF control
   !                nlayers, NGREEK_MOMENTS_INPUT, level_heights, lambertian_albedo(ww), & ! Control Proxies
   !                DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,         & ! Atmos-Optical Proxies
   !                GEMSTOOL_SLEAVE_Results, GEMSTOOL_BRDF_Results,                     & ! Surface Supplement Results
   !                VLIDORT_Out )                                                         ! VLIDORT output results
                 
   !             call GEMSTOOL_RTCALC_Final &
   !             ( MAXWAVNUMS, maxmessages, WW, do_SphericalAlbedo,      & ! GEMSTOOL control
   !                nstokes, n_geometries, dir, VLIDORT_Out,             & ! VLIDORT Results
   !                STOKES_ic, ACTINIC_ic, REGFLUX_ic, DOLP_ic, DOCP_ic, & ! Main program, GEMSTOOL Stokes Output
   !                Errorstatus, nmessages, messages_Main )                ! Main program, exception handling
   
   !          endif
   
   ! !      if ( w.eq.1 ) then
   ! !         write(*,*)aerosol_deltau(18,w)
   ! !         write(*,*) VLIDORT_LinOut%Prof%TS_PROFILEWF(5,18,1,1,1,1), VLIDORT_Out%Main%TS_STOKES(1,1,1,1)
   ! !      endif
   
   ! !  exit if serious error
   
   !          if ( Errorstatus .eq. 1 ) then
   !             write(c1,'(I1)')i_scene
   !             messages_Main(nmessages+1) = 'Failure for scene # '//c1
   !             nmessages = nmessages+1
   !             go to 4555
   !          endif

   !          print*,999
   !          pause




enddo  !loop for wavelength !DO ww = 1, nwav_ext

enddo  !loop for scene ! do i_scene_lrrs = 1, n_scenes_lrrs


! open(99,file='./GEMSTOOL_NSW_Results/temp_GREEKMAT.out', status='replace')
! do ww=1, nwav_ext  !loop for wavelength !DO ww = 1, nwav_ext
! write(99,990) ww, wav_ext(ww), GREEKMAT_TOTAL_INPUT_ext(0,1,ww), GREEKMAT_TOTAL_INPUT_ext(2,1,ww)
! ENDDO
! close(99)

! 990 format(i5,f10.2,1p2e20.10)


! Solar flux read and interpol to wav_ext lines
! ;----GeffToon database
! N_solar_read = 300001 !12000-15000 0.01 cm-1 resolution. ;-- GeffToon!!
temp0 = '/home/mchoi/OSSE_tool/data_solar_spectrum/SOLSPEC_Toon_convolved'
file = trim(temp0)//'/Solar_SOLSPEC_TOON_for_O2AB_07000_15000_mchoi.dat'
OPEN(1,file=trim(file),status='old')
do w=1, N_solar_read
read(1,*) WN_solar_read(w),WL_solar_read(w),Solar_read_cm(w),Solar_read_nm(w),Solar_read_photon(w)
enddo
close(1)


! ! ;----SAO database
! ! N_solar_read = 80093 !12000-15000 0.01 cm-1 resolution. ;-- SAO!!
! temp0 = '/home/mchoi/HIMAP/data_solar_spectrum/SAO2010'
! file = trim(temp0)//'/Solar_SOLSPEC_TOON_WNorder_mchoi.dat' !re-arrange for WN ascending order
! OPEN(1,file=trim(file),status='old')
! do w=1, N_solar_read
! read(1,*) WN_solar_read(w),WL_solar_read(w),Solar_read_cm(w),Solar_read_nm(w),Solar_read_photon(w)
! enddo
! close(1)



! ! After all wavenumber calculation, do LRRS

N2POS = (/ &
-194.3015, -186.4226, -178.5374, -170.6459, -162.7484, &
-154.8453, -146.9368, -139.0233, -131.1049, -123.1819, &
-115.2547, -107.3235, -99.3886, -91.4502, -83.5086, &
-75.5642, -67.6171, -59.6676, -51.7162, -43.7629, &
-35.8080, -27.8519, -19.8950, -11.9373, 11.9373, &
19.8950, 27.8519, 35.8080, 43.7629, 51.7162, &
59.6676, 67.6171, 75.5642, 83.5086, 91.4502, &
99.3886, 107.3235, 115.2547, 123.1819, 131.1049, &
139.0233, 146.9368, 154.8453, 162.7484, 170.6459, &
178.5374, 186.4226, 194.3015 /)

!  O2 spectroscopy Transitions

O2POS = (/ &
-185.5861, -185.5690, -185.5512, -174.3154, -174.2980, &
-174.2802, -163.0162, -162.9980, -162.9809, -151.6906, &
-151.6730, -151.6551, -140.3405, -140.3230, -140.3046, &
-128.9677, -128.9490, -128.9314, -117.5740, -117.5560, &
-117.5372, -108.2632, -106.1613, -106.1430, -106.1240, &
-104.3008, -96.8136, -94.7318, -94.7120, -94.6934, &
-92.8515, -85.3489, -83.2871, -83.2671, -83.2473, &
-81.3868, -73.8694, -71.8294, -71.8080, -71.7876, &
-69.9077, -62.3771, -60.3611, -60.3373, -60.3157, &
-58.4156, -50.8730, -48.8851, -48.8571, -48.8332, &
-46.9116, -40.1158, -39.3565, -37.4069, -37.3688, &
-37.3406, -35.3953, -27.8241, -25.9472, -25.8745, &
-25.8364, -25.8125, -23.8629, -16.2529, -16.2529, &
-14.3761, -14.3033, -14.1686, -12.2918, -2.1392, &
-2.1202, -2.1016, -2.0843, -2.0843, -2.0818, &
-2.0614, -2.0398, -2.0159, -2.0116, -1.9877, &
-1.9735, -1.9496, -1.9455, -1.9217, -1.9003, &
-1.8803, -1.8768, -1.8605, -1.8422, -0.1347, &
0.1347, 1.8422, 1.8605, 1.8768, 1.8803, &
1.9003, 1.9217, 1.9455, 1.9496, 1.9735, &
1.9877, 2.0116, 2.0159, 2.0398, 2.0614, &
2.0818, 2.0843, 2.0843, 2.1016, 2.1202, &
2.1392, 12.2918, 14.1686, 14.3033, 14.3761, &
16.2529, 16.2529, 23.8629, 25.8125, 25.8364, &
25.8745, 25.9472, 27.8241, 35.3953, 37.3406, &
37.3688, 37.4069, 39.3565, 40.1158, 46.9116, &
48.8332, 48.8571, 48.8851, 50.8730, 58.4156, &
60.3157, 60.3373, 60.3611, 62.3771, 69.9077, &
71.7876, 71.8080, 71.8294, 73.8694, 81.3868, &
83.2473, 83.2671, 83.2871, 85.3489, 92.8515, &
94.6934, 94.7120, 94.7318, 96.8136, 104.3008, &
106.1240, 106.1430, 106.1613, 108.2632, 115.7318, &
117.5372, 117.5560, 117.5740, 119.6952, 128.9314, &
128.9490, 128.9677, 140.3046, 140.3230, 140.3405, &
151.6551, 151.6730, 151.6906, 162.9809, 162.9980, &
163.0162, 174.2802, 174.2980, 174.3154, 185.5512, &
185.5690, 185.5861, 196.7919, 196.8100, 196.8269 /)

DO W = 1, nwav

   ! Initialization
   ! nwav_shift = 234 ! defined in header
   ! wav_shift_wl = 0.0d0
   ! wav_shift_wn = 0.0d0
   
   !! To calculate shifted wn, wl, w_excit
   !  Initialize
   LAMBDAS_RANKED = 0.0d0
   NPOINTS_MONO   = 0
   W_EXCIT_WL     = 0
   W_EXCIT_WN     = 0
   
   LAMBDA_EXCIT = 1.0D+07 /wav(w)
   
   !  baseline values of GAMMA
   SIGMA_0 = 1.0D+07 / LAMBDA_EXCIT
   !  initialize counts
   S  = 0
   !  Shifted wavelengths (N2)
   DO J = 1, N_RRS_N2REF
   S = S + 1
   SIGMA = SIGMA_0 + N2POS(J)
   LAMBDAS_ALL(S) = 1.0D+7 / SIGMA
   ENDDO
   !  Shifted wavelengths (O2)
   DO J = 1, N_RRS_O2REF
   S = S + 1
   SIGMA = SIGMA_0 + O2POS(J)
   LAMBDAS_ALL(S) = 1.0D+7 / SIGMA
   ENDDO
   
   !  Total number of shifts
   NSHIFTS = N_RRS_N2REF + N_RRS_O2REF
   NPOINTS_MONO = NSHIFTS + 1
   LAMBDAS_ALL(NPOINTS_MONO) = LAMBDA_EXCIT
   !  Use INDEXX to rank the wavelengths
   !  Find the position of the excitation wavelength in the rank
   
   CALL INDEXX ( NPOINTS_MONO, LAMBDAS_ALL, INDEX_ALL)
   
   DO WW = 1, nwav_shift
      LAMBDAS_RANKED(WW) = LAMBDAS_ALL(INDEX_ALL(WW))
      IF ( LAMBDAS_RANKED(WW) .EQ. LAMBDA_EXCIT) W_EXCIT_WL = WW ! always 120
   enddo
   DO WW = 1, nwav_shift
      WAVENUM_RANKED(WW) = 1.0D+7 / LAMBDAS_RANKED(nwav_shift+1-WW) 
      IF ( WAVENUM_RANKED(WW) .EQ. SIGMA_0) W_EXCIT_WN = WW ! always 115
   ENDDO


   wav_shift_WN(1:nwav_shift)=WAVENUM_RANKED(1:nwav_shift)
   wav_shift_WL(1:nwav_shift)=LAMBDAS_RANKED(1:nwav_shift)
   
   
   !!!!-------- from here! RRS calculation start!!!
   ! To do list
   ! 1. Preparation of Solar Flux (unit check)
   ! 2. LRRS mono calculation...
   


   asym_c(1:maxwavnums) = 0.0d0
   ! FWHM_c(1:maxwavnums) = 0.51953125 ! cm-1 ~ 0.03 nm at 760 nm
   ! FWHM_c(1:maxwavnums) = 1.73 ! cm-1 ~ 0.1 nm at 760 nm
   FWHM_c(1:maxwavnums) = 8.6503906 ! cm-1 ~ 0.5 nm at 760 nm
   ! FWHM_c(1:maxwavnums) = 17.290039 ! cm-1 ~ 1.0 nm at 760 nm
   
   ! FWHM_c(1:maxwavnums) = 0.69 ! nm if....0.04nm res.
   ! FWHM_c(1:maxwavnums) = 0.1d0 ! nm original wav_ext res.
   
   hw1e_c(1:maxwavnums) = FWHM_c / (2.0 * sqrt(log(2.0)))
   
   ! w=1
   ! write(*,*)WN_solar_read(w),WL_solar_read(w),Solar_read_cm(w),Solar_read_nm(w),Solar_read_photon(w)
   ! w=N_solar_read
   ! write(*,*)WN_solar_read(w),WL_solar_read(w),Solar_read_cm(w),Solar_read_nm(w),Solar_read_photon(w)
   ! pause
   
   
   
   do_gauss_conv_before_LRRS_calculation = .false.

   if (do_gauss_conv_before_LRRS_calculation) then 


   ! first convolution test: solar flux
   ! for convolution: (fwave, fspec, nf, nc, hw1e, asym, cwave, cspec)
   ! fwave: wn_solar_read
   ! fspec: solar_read_photon
   ! nf: N_solar_read
   ! nc: nwav_shift or nwav
   ! hw1e: half-width at 1/e intensity
   ! ------FWHM = 2.0 * sqrt(ln(2.0)) * hw1e = 1.66551 * hw1e
   ! asym: should be 0 for gaussian
   ! cwave: wav_shift_WN, or wav
   ! cspec: solar_flux_c(1:nwav_shift)

   !---convolution of solar spectrum for wav(w)
   Call asym_gauss_f2c_1 (wn_solar_read(1:N_solar_read), solar_read_photon(1:N_solar_read), N_solar_read, &
   nwav, hw1e_c(1:nwav), asym_c(1:nwav), wav(1:nwav), Solar_Flux_all(1:nwav) )

!    open(999,file='./GEMSTOOL_NSW_Results/test_solar_convolved_fortran_FWHM8.6503906.out',status='replace')
!    do ww = 1, nwav
!       write(999,700)ww, wav(ww), Solar_Flux_all(ww)
!    enddo
!    close(999)

!    open(999,file='./GEMSTOOL_NSW_Results/test_solar_original.out',status='replace')
!    do ww = 1, N_solar_read
!       write(999,700)ww, WN_solar_read (ww), Solar_read_photon(ww)
!    enddo
!    close(999)

! 700 format(i5,f10.2,1p4e20.10)

!    pause



   !-------Convolution for other parameter for nwav_shift
   !----1) convolution based on WN grid/order -> 2) reverse them to the wl order
   !----Solar_flux_c                        
   Call asym_gauss_f2c_1 (wn_solar_read(1:N_solar_read), solar_read_photon(1:N_solar_read), N_solar_read, &
                          nwav_shift, hw1e_c(1:nwav_shift), asym_c(1:nwav_shift), &
                          wav_shift_WN(1:nwav_shift), Solar_flux_c(1:nwav_shift) )
   ! do ww = 1, nwav_shift-1
   ! write(*,*) wav_shift_WN(ww+1) - wav_shift_WN(ww)
   ! enddo
   ! pause
   
   !---RAYLEIGH_XSEC_c, RAYLEIGH_DEPOL_c, Lambertian_albedo_c
   Call asym_gauss_f2c_1 (wav_ext(1:nwav_ext), RAYLEIGH_XSEC(1:nwav_ext), nwav_ext, &
                          nwav_shift, hw1e_c(1:nwav_shift), asym_c(1:nwav_shift), &
                          wav_shift_WN(1:nwav_shift), RAYLEIGH_XSEC_c(1:nwav_shift) )
   Call asym_gauss_f2c_1 (wav_ext(1:nwav_ext), RAYLEIGH_DEPOL(1:nwav_ext), nwav_ext, &
                          nwav_shift, hw1e_c(1:nwav_shift), asym_c(1:nwav_shift), &
                          wav_shift_WN(1:nwav_shift), RAYLEIGH_DEPOL_c(1:nwav_shift) )
   Call asym_gauss_f2c_1 (wav_ext(1:nwav_ext), Lambertian_albedo(1:nwav_ext), nwav_ext, &
                          nwav_shift, hw1e_c(1:nwav_shift), asym_c(1:nwav_shift), &
                          wav_shift_WN(1:nwav_shift), Lambertian_albedo_c(1:nwav_shift) )
   
   
   Call asym_gauss_f2c_T (wav_ext(1:nwav_ext), DELTAU_VERT_INPUT_ext(1:nlayers, 1:nwav_ext),nwav_ext,nlayers, &
                        nwav_shift, hw1e_c(1:nwav_shift), asym_c(1:nwav_shift), &
                        wav_shift_WN(1:nwav_shift), DELTAU_VERT_INPUT_ext_c(1:nlayers, 1:nwav_shift))
   
! currently best
   do wb = 1, nlayers
      temp_moments_ext(0:NGREEK_MOMENTS_INPUT,1:nwav_ext) = GREEKMAT_TOTAL_INPUT_ext(0:NGREEK_MOMENTS_INPUT,wb,1:nwav_ext)
      Call asym_gauss_f2c_T (wav_ext(1:nwav_ext), temp_moments_ext(0:NGREEK_MOMENTS_INPUT,1:nwav_ext), &
                           nwav_ext,NGREEK_MOMENTS_INPUT+1, &
                           nwav_shift, hw1e_c(1:nwav_shift), asym_c(1:nwav_shift), &
                           wav_shift_WN(1:nwav_shift), temp_moments_ext_c(0:NGREEK_MOMENTS_INPUT,1:nwav_shift))
      GREEKMAT_TOTAL_INPUT_ext_c(0:NGREEK_MOMENTS_INPUT,wb,1:nwav_shift) = temp_moments_ext_c(0:NGREEK_MOMENTS_INPUT,1:nwav_shift)           
   enddo



   !---reverse
   Solar_photon_shift = 0.0d0
   RAYLEIGH_XSEC_shift = 0.0d0
   RAYLEIGH_DEPOL_shift = 0.0d0
   RAYLEIGH_DEPOL_shift = 0.0d0
   Lambertian_albedo_shift = 0.0d0
   DELTAU_VERT_INPUT_shift = 0.0d0
   GREEKMAT_TOTAL_INPUT_shift = 0.0d0

   do ww=1, nwav_shift
      Solar_photon_shift(1+nwav_shift-ww) = Solar_flux_c(ww)
      RAYLEIGH_XSEC_shift(1+nwav_shift-ww) = RAYLEIGH_XSEC_c(ww)
      RAYLEIGH_DEPOL_shift(1+nwav_shift-ww) = RAYLEIGH_DEPOL_c(ww)
      Lambertian_albedo_shift(1+nwav_shift-ww) = Lambertian_albedo_c(ww)

      do wb = 1, nlayers 
         DELTAU_VERT_INPUT_shift(wb,1+nwav_shift-ww) = DELTAU_VERT_INPUT_ext_c(wb,ww)
         ! OMEGA_TOTAL_INPUT_shift(wb,1+nwav_shift-ww) = OMEGA_TOTAL_INPUT_ext_c(wb,ww)
         do wc = 0, NGREEK_MOMENTS_INPUT
            GREEKMAT_TOTAL_INPUT_shift(wb,wc,1+nwav_shift-ww) = GREEKMAT_TOTAL_INPUT_ext_c(wc,wb,ww)
         enddo !do wc = 0, maxmoments_input
      enddo !do wb = 1, nlayers    
   
   enddo ! do ww=1, nwav_shift



   else   !if (do_gauss_conv_before_LRRS_calculation) then 
      !---else -> monochromatic calculation. resolution is up to 0.01.

 ! save Solar Flux for wav(W)
      wastart = 1
      wa = wastart ; trawl = .true.
      do while (trawl)
         if ( wav(w) .ge. WN_solar_read(wa) .and. wav(w) .lt. WN_solar_read(wa+1)) then
            trawl = .false.
         else
            wa = wa + 1
            ! write(*,*)  wav_shift_wn(ww), WN_solar_read(wa)
         endif
      enddo
      wa1 = wa ; wa2 = wa+1
      fa1 = ( WN_solar_read(wa2) - wav(w) ) / ( WN_solar_read(wa2) -  WN_solar_read(wa1) )
      fa2 = 1.0d0 - fa1
      ! ww -> 1+nwav_shift-ww ; then cm-1 order to nm order
      Solar_Flux_all(w) = fa1 * Solar_read_photon(wa1) + fa2 * Solar_read_photon(wa2)
   
   
   
      !!!!-------- from here! RRS calculation start!!!
      ! To do list
      ! 1. Preparation of Solar Flux (unit check)
      ! 2. LRRS mono calculation...
      
   
      !--interpolateion of Solar flux to wav_shift_WN points
      wastart = 1
      do ww = 1, nwav_shift
         ! for Solar flux
         wa = wastart ; trawl = .true.
         do while (trawl)
            if ( wav_shift_wn(ww) .ge. WN_solar_read(wa) .and. wav_shift_wn(ww) .lt. WN_solar_read(wa+1)) then
               trawl = .false.
            else
               wa = wa + 1
               ! write(*,*)  wav_shift_wn(ww), WN_solar_read(wa)
            endif
         enddo
         wa1 = wa ; wa2 = wa+1
         fa1 = ( WN_solar_read(wa2) - wav_shift_wn(ww) ) / ( WN_solar_read(wa2) -  WN_solar_read(wa1) )
         fa2 = 1.0d0 - fa1
         ! ww -> 1+nwav_shift-ww ; then cm-1 order to nm order
         Solar_photon_shift(1+nwav_shift-ww) = fa1 * Solar_read_photon(wa1) + fa2 * Solar_read_photon(wa2)
         ! write(*,*) wa, wav_shift_wn(ww), WN_solar_read(wa), WN_solar_read(wa+1), fa1, fa2
         ! pause
   
         ! for ext properties
         ! Extended Variables
         ! --DELTAU_VERT_INPUT_ext, OMEGA_TOTAL_INPUT_ext: ( MAXLAYERS, maxwavnums)
         ! --GREEKMAT_TOTAL_INPUT_ext:  (0:MAXMOMENTS_INPUT, MAXLAYERS, maxwavnums)
         ! Target Variables
         ! --DELTAU_VERT_INPUT_shift, OMEGA_TOTAL_INPUT_shift: ( MAXLAYERS, nwav_shift)
         ! --GREEKMAT_TOTAL_INPUT_shift:  (0:MAXMOMENTS_INPUT, MAXLAYERS, nwav_shift)
   
         wa = wastart ; trawl = .true.
         do while (trawl)
            if ( wav_shift_wn(ww) .ge. wav_ext(wa) .and. wav_shift_wn(ww) .lt. wav_ext(wa+1)) then
               trawl = .false.
               ! write(*,*)  wav_shift_wn(ww), wav_ext(wa),wav_ext(wa+1)
            else
               wa = wa + 1
            endif
         enddo !do while (trawl)
         wa1 = wa ; wa2 = wa+1
         fa1 = ( wav_ext(wa2) - wav_shift_wn(ww) ) / ( wav_ext(wa2) -  wav_ext(wa1) )
         fa2 = 1.0d0 - fa1
         ! write(*,*)Fa1,Fa2
   
         ! ww -> 1+nwav_shift-ww ; then cm-1 order to nm order
         RAYLEIGH_XSEC_shift(1+nwav_shift-ww) = fa1 * RAYLEIGH_XSEC(wa1) + fa2 * RAYLEIGH_XSEC(wa2)
         RAYLEIGH_DEPOL_shift(1+nwav_shift-ww) = fa1 * RAYLEIGH_DEPOL(wa1) + fa2 * RAYLEIGH_DEPOL(wa2)
         Lambertian_albedo_shift(1+nwav_shift-ww) = fa1 * Lambertian_albedo(wa1) + fa2 * Lambertian_albedo(wa2)
   
         do wb = 1, nlayers
            DELTAU_VERT_INPUT_shift(wb,1+nwav_shift-ww) = fa1 * DELTAU_VERT_INPUT_ext(wb,wa1) + &
                                                          fa2 * DELTAU_VERT_INPUT_ext(wb,wa2)
            OMEGA_TOTAL_INPUT_shift(wb,1+nwav_shift-ww) = fa1 * OMEGA_TOTAL_INPUT_ext(wb,wa1) + &
                                                          fa2 * OMEGA_TOTAL_INPUT_ext(wb,wa2)
            do wc = 0, maxmoments_input !GREEKMAT_TOTAL_INPUT = Greekmat * omega!! mchoi (ref. from o3bin_m_iopsetup.f90)
               GREEKMAT_TOTAL_INPUT_shift(wb,wc,1+nwav_shift-ww) = &
                      fa1 * GREEKMAT_TOTAL_INPUT_ext(wc,wb,wa1) + &
                      fa2 * GREEKMAT_TOTAL_INPUT_ext(wc,wb,wa2) 
                           ! change dimension order wc,wb -> wb,wc
            enddo !do wc = 0, maxmoments_input
         enddo !do wb = 1, nlayers
      enddo !do ww = 1, nwav_shift
   endif !!if (do_gauss_conv_before_LRRS_calculation) then // else

   !link to FastV2p7_LRRS2p3_Master module!! 

    !Call to LRRS master for each geometry, with exception handling
      EarthRadius = GEMSTOOL_INPUTS%TimePos%earthradius
      sza_boa(1) = GEMSTOOL_INPUTS%Geometry%GEMS_szas(1)
      vza_boa(1) = GEMSTOOL_INPUTS%Geometry%GEMS_vzas(1)
      azm_boa(1) = GEMSTOOL_INPUTS%Geometry%GEMS_azms(1)

      N_LOUTPUT_read = Nint(GEMSTOOL_INPUTS%RTControl%observation_height)
      
      ! ! fluxes_c(w) = 0.07/10**7
      ! ! do w=1, nwav
      ! ! lambda_c(w) =  wav_shift_WN(ww)!10.0**7/wav_shift_WN(ww)
      ! ! enddo
      ! ! Lambertian_albedo(w) = 0.05d0
      ! ! write(*,*) sza_boa, vza_boa, azm_boa, EarthRadius, 10.0**7/wav_shift_WN(ww)

      ! ! do v = 1, ngeoms !<- should be 1 to 1

      ! ! write(*,*) lambda_c(w), fluxes_c(w)

      call FastV2p7_LRRS2p3_Master &
      ( nwav_shift, maxlayers, MAXMOMENTS_INPUT,                 &
      nwav_shift, nlayers, 8, 500, EarthRadius,     &    
      sza_boa(1), vza_boa(1), azm_boa(1), LEVEL_HEIGHTS, LEVEL_TEMPS, AIRCOLUMNS,     &
      wav_shift_WL(1:nwav_shift), Solar_photon_shift, Lambertian_albedo_shift, & 
      RAYLEIGH_XSEC_shift, Rayleigh_depol_shift, DELTAU_VERT_INPUT_shift, GREEKMAT_TOTAL_INPUT_shift, &
      elastic_c, raman_c, OMI_filling_c, Filling_SS_c,      &
      NPOINTS_MONO, W_EXCIT_WL, LAMBDA_EXCIT, N_LOUTPUT_read, & !by mchoi; reading
      fail_LRRS_Raman, n_LRRS_messages, LRRS_messages )
      ! !! --for reference
      ! !   ( E_ndat_c, E_nlayers, E_nmoms_all,                 &
      ! !     ndat_c, LRRS_nlayers, nstreams, nmoments_LRRS, EarthRadius,     &    
      ! !     sza_boa(v), vza_boa(v), azm_boa(v), heights, temperatures, Aircolumns,     &
      ! !     lambdas_c, OMI_solarflux_c, LRRS_Albedo, RayXsec_c, Depol_c, deltau_c, omegamoms_c, &
      ! !     elastic_c, raman_c, OMI_filling_c(:,v,iout), Filling_SS_c,      &
      ! !     fail_LRRS_Raman, n_LRRS_messages, LRRS_messages )

      if ( Fail_LRRS_Raman ) then
      do m = 1, n_LRRS_messages
      b_messages(b_nmessages+m) = Trim(LRRS_messages(m))
      enddo
      b_nmessages = b_nmessages + n_LRRS_messages
      b_messages(b_nmessages+1) = 'LRRS Raman Call failed for Scenario'
      b_nmessages = b_nmessages + 1
      b_Fail = .true. ; go to 69
      !endif

      else

      ! write(*,*) elastic_c(:)
      ! write(*,*)
      ! write(*,*) Raman_c(:)
      ! pause
      ! write(*,*) 
      ! write(*,*) raman_c(W_EXCIT_WL)
      ! write(*,*) 
      ! write(*,*) OMI_filling_c(W_EXCIT_WL)

      elastic_all(w)=elastic_c(1)
      raman_all(w)=raman_c(1)
      OMI_filling_all(w)=OMI_filling_c(1) ! Def: (one - E/R)
      ! write(*,*) OMI_filling_c(1)
      ! write(*,*) Filling_SS_c(1)

      write(*,*)w,wav(w),OMI_filling_all(w),Solar_Flux_all(w)
      endif !if ( Fail_LRRS_Raman ) then

      ! enddo

               
enddo  !loop for wavelength !DO w = 1, nwav

! write LRRS results   
   

open(999,file='./GEMSTOOL_NSW_Results/LRRS_results_mchoi.out',status='replace')
   do w = 1, nwav
      write(999,900)w, wav(w), elastic_all(w), raman_all(w), OMI_filling_all(w), Solar_Flux_all(w)
   enddo
close(999)

900 format(i5,f10.2,1p4e20.10)

!      close(340)
! 567 format(i5,f10.2,1p8e17.7)

!  Debug
!   if ( FDBUL .eq. 0 ) then
!      do w = 1, nwav
!         if (w.eq.1)write(443,'(12F11.7)')(AerBulk_pars(q),q=1,n_AerBulk_Total_wfs)
!         write(443,'(i3,1pe20.12,1p12e15.6)')w,STOKES_ic(1,1,w),(radiance_aerbulkwfs_ic(1,w,q),q=1,n_AerBulk_Total_wfs)
!      enddo
!   endif

!  End Scene loop



   endif ! if (do_LRRS_cal) then ;--should moved to later
   enddo !loop for scene ! do i_scene = 1, n_scenes




!  ============================
!       Write output
!  ===========================

!  12 August 2013. Use NSW-specific path for the Results.   TO DO !!!!!!!!!!!!!!!!!!!!!

   clevel = COUTPUT(UpDn_index)

!  Radiance and Stokes-vector output
!  ---------------------------------
!  Original before changing output altitude
   !    if ( nstokes .eq. 1 ) then
   !       open(1,file='./GEMSTOOL_NSW_Results/Stokes_I_Scalar_'//clevel//'_linearized.out',status='replace')
   !       do w = 1, nwav
   !          write(1,45)w, wav_shift_WN(ww),(stokes(k,1,w),k=1,n_geometries)
   !       enddo
   !    else
   !       open(1,file='./GEMSTOOL_NSW_Results/Stokes_IQU_DoLP_'//clevel//'.out',status='replace')
   !       do w = 1, nwav
   !          write(1,46)w, wav_shift_WN(ww),((stokes(k,q,w),k=1,n_geometries),q=1,nstokes),(dolp(k,w),k=1,n_geometries)
   !       enddo
   !    endif
   !    close(1)
   ! !45 format(1pe20.10, 1p3e20.10e3)
   ! !46 format(1pe20.10, 1p5e20.10e3)
   ! 45 format(i5,f10.2,1p3e20.10)
   ! 46 format(i5,f10.2,1p12e20.10)

!  After changing output height   
   if ( nstokes .eq. 1 ) then
      open(1,file='./GEMSTOOL_NSW_Results/Stokes_I_Scalar_'//clevel//'_linearized.out',status='replace')
      do w = 1, nwav
         write(1,35)w, wav(w),GEMSTOOL_INPUTS%RTControl%observation_height,&
                              (stokes(k,1,w),k=1,n_geometries)
      enddo
   else
      open(1,file='./GEMSTOOL_NSW_Results/Stokes_IQU_DoLP_'//clevel//'.out',status='replace')
      do w = 1, nwav
         write(1,36)w, wav(w),GEMSTOOL_INPUTS%RTControl%observation_height,&
                              ((stokes(k,q,w),k=1,n_geometries),q=1,nstokes),(dolp(k,w),k=1,n_geometries)
      enddo
   endif
   close(1)
!45 format(1pe20.10, 1p3e20.10e3)
!46 format(1pe20.10, 1p5e20.10e3)
   
! 45 format(i5,f10.2,1p3e20.10)
! 46 format(i5,f10.2,1p12e20.10)

! 35 format(i5,f10.2,f8.3,1p3e20.10)
! 36 format(i5,f10.2,f8.3,1p12e20.10)

!editted in t12
45 format(i5,f10.2,1p3e20.10e3)
46 format(i5,f10.2,1p12e20.10e3)

35 format(i5,f10.2,f8.3,1p3e20.10e3)
36 format(i5,f10.2,f8.3,1p12e20.10e3)


! open(1,file='./GEMSTOOL_NSW_Results/Raman_'//clevel//'.out',status='replace')
!       do w = 1, nwav
!          write(1,37)w, wav_shift_WN(ww), elastic_c(w), raman_c(w)
!       enddo
! close(1)
! 37 format(i5,f10.2,1p12e20.10)



!  Flux output (currently no Jacobians for these quantities)
!  -----------

!  Only the intensity fluxes 

   !if ( do_SphericalAlbedo ) then
   !   open(1,file='./GEMSTOOL_NSW_Results/Stokes_ActinicFlux_'//clevel//'.out',status='replace')
   !   open(2,file='./GEMSTOOL_NSW_Results/Stokes_RegularFlux_'//clevel//'.out',status='replace')
   !   do w = 1, nwav
   !      write(1,45)w, wav_shift_WN(ww),(Actinic(k,1,w),k=1,n_szas)
   !      write(2,45)w, wav_shift_WN(ww),(Regflux(k,1,w),k=1,n_szas) 
   !   enddo
   !endif

!  Jacobian output
!  ---------------

!  Rob Upgrade, 12/19/16. Write-up of Stokes-Q and Stokes-U results

!  Trace Gas Profile Jacobians, ONE GEOMETRY ONLY
!     Output also the AMF scattering weights

   if ( do_GasProfile_Jacobians ) then
        !print*, do_gases !! TFTF
        !print*, do_gas_wfs !! TTTT
     do g = 1, ngases
       !if ( do_gases(g) .and. do_gas_wfs(g) ) then
         open(1,file=TRIM('./GEMSTOOL_NSW_Results/Stokes_I_GasWFs_'//clevel//'.out'//which_gases(g)),status='replace')
!         open(2,file=TRIM('./GEMSTOOL_NSW_Results/Stokes_I_AMFScatWts_'//clevel//'.out'//which_gases(g)),status='replace')
!         open(3,file=TRIM('./GEMSTOOL_NSW_Results/Stokes_I_AMFTotal_'//clevel//'.out'//which_gases(g)),status='replace')

         do k = 0, nlayers
!           write(1,47)k,level_heights(k),(radiance_gaswfs(1,w,k,g),w = 1,nwav)
!           if (k.gt.0)write(2,47)k,half*(level_heights(k)+level_heights(k-1)),(AMFSCATWTS(1,w,k,g),w = 1,nwav)
           do w = 1, nwav
              if ( FDGas.eq.g.and.FDLev.eq.k) write(23,*)w,Level_vmrs(k,g),stokes(1,1,w),radiance_gaswfs(1,w,k,g)
           enddo
         enddo

         do w = 1, nwav
             write(1,48)w, wav(w), (radiance_gaswfs(1,w,k,g),k=0,nlayers)
         enddo

         close(1) !; close(2) ; close(3)

         if (nstokes.gt.1) then
           open(1,file=TRIM('./GEMSTOOL_NSW_Results/Stokes_Q_GasWFs_'//clevel//'.out'//which_gases(g)),status='replace')
           open(2,file=TRIM('./GEMSTOOL_NSW_Results/Stokes_U_GasWFs_'//clevel//'.out'//which_gases(g)),status='replace')
           do w = 1, nwav
             write(1,48)w, wav(w), (Stokes_Q_gaswfs(1,w,k,g),k=0,nlayers)
             write(2,48)w, wav(w), (Stokes_U_gaswfs(1,w,k,g),k=0,nlayers)
           enddo
           close(1) ; close(2)
         endif

!       endif
     enddo
   endif
!47 format(i5,f9.3,1p1000e20.10e3)
!48 format(1pe20.10,1p1000e20.10e3)
48 format(i5,f10.2,1p30000e20.10e3)

!  Aerosol profile Jacobians (LAYER JACOBIAN). ONE GEOMETRY ONLY

   if ( do_AerOpdepProfile_Jacobians ) then
      open(1,file=TRIM('./GEMSTOOL_NSW_Results/Stokes_I_AerProfWFs_'//clevel//'.out'),status='replace')
      do k = 1, nlayers
!         write(1,47)k,level_heights(k),(radiance_aerprofwfs(1,w,k),w = 1,nwav)
         do w = 1, nwav
            if ( FDAer.and.FDLay.eq.k ) write(73,*)w,aertau_unscaled(k), stokes(1,1,w),radiance_aerprofwfs(1,w,k)
         enddo
      enddo
      do w = 1, nwav
        write(1,48)w, wav(w), (radiance_aerprofwfs(1,w,k),k=1,nlayers)
      enddo
      close(1)

      if (nstokes.gt.1) then
        open(1,file=TRIM('./GEMSTOOL_NSW_Results/Stokes_Q_AerProfWFs_'//clevel//'.out'),status='replace')
        open(2,file=TRIM('./GEMSTOOL_NSW_Results/Stokes_U_AerProfWFs_'//clevel//'.out'),status='replace')
        do w = 1, nwav
          write(1,48)w, wav(w), (Stokes_Q_aerprofwfs(1,w,k),k=1,nlayers)
          write(2,48)w, wav(w), (Stokes_U_aerprofwfs(1,w,k),k=1,nlayers)
        enddo
        close(1) ; close(2)
      endif

   endif

!  Surface Jacobians. ALL GEOMETRIES

   if ( do_Surface_Jacobians ) then
      if (.not. do_BRDF) then
         open(1,file='./GEMSTOOL_NSW_Results/Stokes_I_AlbWFs_'//clevel//'.out',status='replace')
         do w = 1, nwav
            write(1,45)w, wav(w),(radiance_Albedowfs(k,w,1),k=1,n_geometries)
         enddo
         close(1)

         if (nstokes.gt.1) then
         open(1,file='./GEMSTOOL_NSW_Results/Stokes_Q_AlbWFs_'//clevel//'.out',status='replace')
         open(2,file='./GEMSTOOL_NSW_Results/Stokes_U_AlbWFs_'//clevel//'.out',status='replace')
         do w = 1, nwav
            write(1,45)w, wav(w), (Stokes_Q_Albedowfs(k,w,1),k=1,n_geometries)
            write(2,45)w, wav(w), (Stokes_U_Albedowfs(k,w,1),k=1,n_geometries)
         enddo
         close(1) ; close(2)
         endif

      else !if (do_BRDF) then
         
         open(1,file='./GEMSTOOL_NSW_Results/Stokes_I_BRDFWFs_'//clevel//'.out',status='replace')
         do w = 1, nwav
            write(1,41)w, wav(w),((radiance_BRDFwfs(k,w,qs),k=1,n_geometries),qs = 1, N_SURFACE_WFS_fin)
         enddo
         close(1)

         if (nstokes.gt.1) then
         open(1,file='./GEMSTOOL_NSW_Results/Stokes_Q_BRDFWFs_'//clevel//'.out',status='replace')
         open(2,file='./GEMSTOOL_NSW_Results/Stokes_U_BRDFWFs_'//clevel//'.out',status='replace')
         do w = 1, nwav
            write(1,41)w, wav(w), ((Stokes_Q_BRDFwfs(k,w,qs),k=1,n_geometries),qs = 1, N_SURFACE_WFS_fin)
            write(2,41)w, wav(w), ((Stokes_U_BRDFwfs(k,w,qs),k=1,n_geometries),qs = 1, N_SURFACE_WFS_fin)
         enddo
         close(1) ; close(2)
         endif
         

      endif !if (do_BRDF) then
   endif

   if (do_BRDF) then
      !---write input BRDF Kernel_Factor_Parameters
      open(11,file='./GEMSTOOL_NSW_Results/Input_BRDF_Kernel_Factor_Parameter.dat', status='replace')
      do ii = 1, N_SURFACE_WFS_fin
         write(11,'(f11.5)') BRDF_K_FP_values(ii)
      enddo
      close(11)
   endif !if (do BRDF) then

!  Rob Fix, 10/18/16, SIF Surface-leaving Jacobians. ALL GEOMETRIES.
!           11/30/16, update for linear parameteriation, rename output arrays --> SIFParamWFs

   if ( do_SIF_Jacobians ) then

      open(1,file='./GEMSTOOL_NSW_Results/Stokes_I_SIFPARAMWFs_'//clevel//'.out',status='replace')
      do w = 1, nwav
         write(1,41)w, wav(w),((radiance_SIFPARAMwfs(k,w,qs),k=1,n_geometries),qs = 1, n_SIFPars)
      enddo
      close(1)

      if (nstokes.gt.1) then
        open(1,file='./GEMSTOOL_NSW_Results/Stokes_Q_SIFPARAMWFs_'//clevel//'.out',status='replace')
        open(2,file='./GEMSTOOL_NSW_Results/Stokes_U_SIFPARAMWFs_'//clevel//'.out',status='replace')
        do w = 1, nwav
          write(1,41)w, wav(w),((Stokes_Q_SIFPARAMwfs(k,w,qs),k=1,n_geometries),qs = 1, n_SIFPars)
          write(2,41)w, wav(w),((Stokes_U_SIFPARAMwfs(k,w,qs),k=1,n_geometries),qs = 1, n_SIFPars)
        enddo
        close(1) ; close(2)
      endif

   endif
!41 format(1pe20.10,1p6e20.10e3)
41 format(i5,f10.2,1p20e20.10e3)

!  T-shift Jacobians. ALL GEOMETRIES

   if ( do_Tshift_Jacobian ) then

      open(1,file='./GEMSTOOL_NSW_Results/Stokes_I_TSHIFTWFs_'//clevel//'.out',status='replace')
      do w = 1, nwav
         write(1,45)w, wav(w),(radiance_Tshiftwfs(k,w),k=1,n_geometries)
      enddo
      close(1)

      if (nstokes.gt.1) then
        open(1,file='./GEMSTOOL_NSW_Results/Stokes_Q_TSHIFTWFs_'//clevel//'.out',status='replace')
        open(2,file='./GEMSTOOL_NSW_Results/Stokes_U_TSHIFTWFs_'//clevel//'.out',status='replace')
        do w = 1, nwav
          write(1,45)w, wav(w), (Stokes_Q_Tshiftwfs(k,w),k=1,n_geometries)
          write(2,45)w, wav(w), (Stokes_U_Tshiftwfs(k,w),k=1,n_geometries)
        enddo
        close(1) ; close(2)
      endif

   endif

!  Surface pressure Jacobians. ALL GEOMETRIES

   if ( do_Surfpress_Jacobian ) then

      open(1,file='./GEMSTOOL_NSW_Results/Stokes_I_SURFPWFs_'//clevel//'.out',status='replace')
      do w = 1, nwav
         write(1,45)w, wav(w),(radiance_Surfpwfs(k,w),k=1,n_geometries)
      enddo
      close(1)

      if (nstokes.gt.1) then
        open(1,file='./GEMSTOOL_NSW_Results/Stokes_Q_SURFPWFs_'//clevel//'.out',status='replace')
        open(2,file='./GEMSTOOL_NSW_Results/Stokes_U_SURFPWFs_'//clevel//'.out',status='replace')
        do w = 1, nwav
          write(1,45)w, wav(w), (Stokes_Q_Surfpwfs(k,w),k=1,n_geometries)
          write(2,45)w, wav(w), (Stokes_U_Surfpwfs(k,w),k=1,n_geometries)
        enddo
        close(1) ; close(2)
      endif

   endif

!  H2O Scaling Jacobians. ALL GEOMETRIES
!  @@ Rob fix 7/23/14

   if ( do_H2OScaling_Jacobian ) then

      open(1,file='./GEMSTOOL_NSW_Results/Stokes_I_WSCALEWFs_'//clevel//'.out',status='replace')
      do w = 1, nwav
         write(1,45)w, wav(w),(radiance_WScalewfs(k,w),k=1,n_geometries)
      enddo
      close(1)

      if (nstokes.gt.1) then
        open(1,file='./GEMSTOOL_NSW_Results/Stokes_Q_WSCALEWFs_'//clevel//'.out',status='replace')
        open(2,file='./GEMSTOOL_NSW_Results/Stokes_U_WSCALEWFs_'//clevel//'.out',status='replace')
        do w = 1, nwav
          write(1,45)w, wav(w), (Stokes_Q_WScalewfs(k,w),k=1,n_geometries)
          write(2,45)w, wav(w), (Stokes_U_WScalewfs(k,w),k=1,n_geometries)
        enddo
        close(1) ; close(2)
      endif

   endif

!  CH4 Scaling Jacobians. ALL GEOMETRIES
!  @@ Y.Jung fix 2/1/15

   if ( do_CH4Scaling_Jacobian ) then

      open(1,file='./GEMSTOOL_NSW_Results/Stokes_I_MSCALEWFs_'//clevel//'.out',status='replace')
      do w = 1, nwav
         write(1,45)w, wav(w),(radiance_MScalewfs(k,w),k=1,n_geometries)
      enddo
      close(1)

      if (nstokes.gt.1) then
        open(1,file='./GEMSTOOL_NSW_Results/Stokes_Q_MSCALEWFs_'//clevel//'.out',status='replace')
        open(2,file='./GEMSTOOL_NSW_Results/Stokes_U_MSCALEWFs_'//clevel//'.out',status='replace')
        do w = 1, nwav
          write(1,45)w, wav(w), (Stokes_Q_MScalewfs(k,w),k=1,n_geometries)
          write(2,45)w, wav(w), (Stokes_U_MScalewfs(k,w),k=1,n_geometries)
        enddo
        close(1) ; close(2)
      endif

   endif

!  Aerosol bulk Jacobians

    print*, 'total number of bulk aerosol parameters = ',n_AerBulk_Total_wfs
   if ( do_AerBulk_Jacobians ) then

      open(1,file=TRIM('./GEMSTOOL_NSW_Results/Stokes_I_AerBulkWFs_'//clevel//'.out'),status='replace')
        do w = 1, nwav
            write(1,48) w, wav(w), (radiance_aerbulkwfs(1,w,k),k=1,n_AerBulk_Total_wfs)
        enddo
!      do k = 1, n_AerBulk_Total_wfs
!         write(1,47)k,AerBulk_pars(k),(radiance_aerbulkwfs(1,w,k),w = 1,nwav)
!         do w = 1, nwav
!            if ( FDAer.and.FDLay.eq.k ) write(73,*) w,AerBulk_pars(k), stokes(1,1,w),radiance_aerbulkwfs(1,w,k)
!         enddo
!      enddo
      close(1)

      if (nstokes.gt.1) then
        open(1,file=TRIM('./GEMSTOOL_NSW_Results/Stokes_Q_AerBulkWFs_'//clevel//'.out'),status='replace')
        open(2,file=TRIM('./GEMSTOOL_NSW_Results/Stokes_U_AerBulkWFs_'//clevel//'.out'),status='replace')
        do w = 1, nwav
            write(1,48) w, wav(w), (Stokes_Q_aerbulkwfs(1,w,k),k=1,n_AerBulk_Total_wfs)
            write(2,48) w, wav(w), (Stokes_U_aerbulkwfs(1,w,k),k=1,n_AerBulk_Total_wfs)
        enddo
        close(1) ; close(2)
      endif

   endif

!  normal finish

   go to 4556

!  =====================
!     Error section
!  =====================

4555  continue

   write(*,*)
   write(*,*)'GEMSTOOL_NSW_Linearized_Main: Errors encountered, Look at file: GEMSTOOL_NSW_Errors.LOG'
   write(*,*)

   open(1,file='GEMSTOOL_NSW_Errors.LOG',status='unknown')
   write(1,'(a/a,i5)')'Errors in Setups/Inputs, GEMSTOOL_NSW_Linearized_Main WRAPPER **',&
                      '  * Number of error messages = ',nmessages
   do k = 1, nmessages
     write(1,'(a,i4,a,a)')'  * Message # ',k, ' : ',Adjustl(trim(messages_Main(k)))
   enddo
   close(1)
   stop

!  Finish
!  ======

4556 continue


   open(11,file='./GEMSTOOL_NSW_Results/Input_alt_pres.out', status='replace')
   do ii = 0, nlayers
   !do ii = 1, nlayers
      write(11,'(f11.5,2x,f11.5)')level_heights(ii), level_press(ii)
   enddo
   close(11)

   open(11,file='./GEMSTOOL_NSW_Results/Input_wav.out', status='replace')
      write(11,'(f11.5,2x,f11.5,2x,f11.5)') wav(1),wav(nwav),wav(2)-wav(1)
   close(11)



   write(*,*)
   write(*,*)'GEMSTOOL_NSW_Linearized_Main: Successful Run !!!!!!!!!'
   write(*,*)

   69 continue

   write(*,*)'FastV2p7_Exact_LPS_Driver: aborted, here are messages ---'
   write(*,'(a,i3)')' * Number of error messages = ', b_nmessages
   do m = 1, b_nmessages
      write(*,'(a,i3,a,a)')' -- Message # ', m,' : ',b_messages(m)
   enddo



   stop
end program GEMSTOOL_NSW_Linearized_Main_IQU
