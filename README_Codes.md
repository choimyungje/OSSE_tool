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

  
## Step 2. 
