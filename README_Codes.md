# OSSE_aerosol/README_Codes.md

Manual to run the simulation of OSSE_aerosol code for information contents analysis of aerosol vertical profiling

## Step 1: High spectral resolution spectrum calculation

### 1. Overall guidance
This program handles all processes to calculate high-spectral-resolution spectra over O2 absorption bands. 
- Inputs
  - Set 1: from the practical CLARS-FTS measurement condition over LA basin.  
    - Solar zenith/azimuth angle
    - target longitude/latitude
    - measurement date/time
  - Set 2: that can be changed  
    - Measurement geometry: viewing zenith/azimuth angle
    - Aerosol microphysical properties: size distribution (effective radius/variance, fine-mode fraction), refractive indicies
    - Aerosol vertical distribution: total AOD, peak height and half width of half maximum of Gaussian shape distribution
    - Surface reflectance: Lambertian surface reflectance, BRDF, or BRDF+BPDF
    - Atmospheric vertical layer/resolution: pressure-temperature-height profile, gases profile
    - Number of streams for DISORT calculation 
    - Spectral range/resolution
    - etc..  
  - The parameters are controlled by the master code and others following codes.
- Outputs
  - Location (e.g.): ./Results/All_Results_justcal/Case000_cl201907052231West/ii0000_AOD0.30_Peak01.0_Width0.3_SZA35.50_VZA00.01_RAA030.00_SCA144.50/LAB_step00/
  - Stokes parameters (I, Q, and U) with a unit of normalized radiance: Stokes_IQU_DoLP_TOAUP.out
  - Jacobians for 
    - Aerosol bulk parameters: Stokes_I_AerBulkWFs_TOAUP.out (and Q/U)
    - Surface BRDF parameters: Stokes_I_BRDFWFs_TOAUP.out (and Q/U)
    - H2O scaling: Stokes_I_WSCALEWFs_TOAUP.out (and Q/U)
    - Temperature shift: Stokes_I_TSHIFTWFs_TOAUP.out (and Q/U)
    - Surface pressure: Stokes_I_SURFPWFs_TOAUP.out (and Q/U)
    - SIF parameters (not used in this study)
    - AOD profile (not used in this study)
  - Unitless normalized radiance (def. radiance when incident solar radiance at top-of-atmosphere is assumed 1; unity)  

      
### 2. Code structures
- **osse_aerosol_0_master.pro**   
  - osse_aerosol_1_cal_pthg.pro
    - sub_module_setup_pthg_clars_new_multiple.pro
    - sub_make_xsec_lut_using_absco_v5_mchoi_o2ab_new_multiple.pro
  - osse_aerosol_2_cal_aodprof.pro
    - sub_integrate_result_svo_lab_loop.pro
  - osse_aerosol_3_cal_rtm.pro
    - sub_integrate_result_svo_lab_loop.pro   

  
## Step 2. Estimation of information contents and retrieval uncertainty  

### 1. Overall guidance
This process calculates information contents (i.e., degree of freedom for signal; DOFS) of geophyeical parameters and its retrieval uncertainties for aerosol profiling parameters and other parameters.  
Instrument characteristics are considered. e.g., Specific signal-to-noise ratio (SNR) and spectral resolution (i.e., full-width-half-maximum; FWHM)
- Inputs
  - High-resolution (e.g., 0.02 cm-1 in this study) spectra calculated from Step 1.  
    - Stokes parameters
    - Jacobians
  - A priori uncertainty
    - Aerosol/Surface/other ancilliary parameters uncertainty
    - They can be analyzed from the standard deviation of climatological database
    - Or, they also can be defined as specific number (e.g., 100%)
  - Instrument model
    - Instrument line shape function: SINC function (for FTS-type), Gaussian function (for grating-type)
    - Signal-to-noise ratio (SNR) for radiance (I): e.g., 100-1000
    - SNR for Q and U (i.e., assumed as SNR of radiance / sqrt(2) in this study)
    - Measurement noise can be defined as signal / SNR, or other definition can be used.
- Outputs
  - Averaging Kernels
  - Degree of freedom for signal (DOFS) for each geophysical parameters
  - A posterior error
    - Total retrieval error
      - Smoothing error
      - Measurement noise error
      - Cross-state error (or, interference error)
  - Estimated retrieval bias due to measurement bias
  
### 2. Code structures
- Location: ./Analysis_DOFS
- analysis_dofs_satellite_nadir_ab_tropomi.pro 
  - output location: (e.g.) ./Analysis_DOFS/result_satellite_nadir_ab_tropomi/save_AOD0.10_Peak00.2_Width0.2/DFS_RelE.xdr
  - 100 aerosol loading cases: AOD [0.1, 0.3, 0.5, 1.0]; APH [0.2, 0.6, 1.0, 1.5, 2.0 km]; ALT [0.2, 0.6, 1.0, 1.5, 2.0 km]  
  - 
  - 
  ...
- 
