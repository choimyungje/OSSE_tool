# OSSE_aerosol/README_Codes.md

Manual to run the simulation of OSSE_aerosol code for information contents analysis of aerosol vertical profiling

## 1. Overall guidance  
  
- **osse_aerosol_0_master.pro**  
This "master" code handles all processes to calculate high-spectral-resolution spectra over O2 absorption bands. 
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
    - Some parameters are controlled from the master code and others are controlled from following codes.
  - Outputs
    - High-spectral-resolution Stokes parameters (i.g., I, Q, and U) with a unit of normalized radiance
    - Jacobians for 
      - aerosol bulk parameters/or AOD profile
      - surface parameters
      - H2O scaling
      - Temperature shift
      - Surface pressure
      - SIF parameters
    
    
Outputs of this simulations are as below.

  - osse_aerosol_1_cal_pthg.pro
    - sub_module_setup_pthg_clars_new_multiple.pro
    - sub_make_xsec_lut_using_absco_v5_mchoi_o2ab_new_multiple.pro
  - osse_aerosol_2_cal_aodprof.pro
    - sub_integrate_result_svo_lab_loop.pro
  - osse_aerosol_3_cal_rtm.pro
    - sub_integrate_result_svo_lab_loop.pro  
    
  
