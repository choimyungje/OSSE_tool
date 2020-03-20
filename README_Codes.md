# OSSE_aerosol/README_Codes.md

Manual to run the simulation of OSSE_aerosol code for information contents analysis of aerosol vertical profiling

## 1. Overall guidance  
  
- **osse_aerosol_0_master.pro**  
This "master" code is to handle all processes for calculating high-resolution spectra over O2 absorption bands. The simulation is focused on the real CLARS-FTS measurement over LA basin area, e.g., Solar zenith/azimuth angle, target longitude/latitude, measurement date/time, etc. Other independent variables are as below.  
  - Measurement geometry: viewing zenith/azimuth angle
  - Aerosol microphysical properties: size distribution (effective radius/variance, fine-mode fraction), refractive indicies
  - Aerosol vertical distribution: total AOD, peak height and half width of half maximum of Gaussian shape distribution
  - Surface reflectance: Lambertian surface reflectance, BRDF, or BRDF+BPDF
  - Atmospheric vertical layer/resolution: pressure-temperature-height profile, gases profile
  - Number of streams for DISORT calculation 
  - Spectral range/resolution
  - 
    
Outputs of this simulations are as below.

  - osse_aerosol_1_cal_pthg.pro
    - sub_module_setup_pthg_clars_new_multiple.pro
    - sub_make_xsec_lut_using_absco_v5_mchoi_o2ab_new_multiple.pro
  - osse_aerosol_2_cal_aodprof.pro
    - sub_integrate_result_svo_lab_loop.pro
  - osse_aerosol_3_cal_rtm.pro
    - sub_integrate_result_svo_lab_loop.pro  
    
  
