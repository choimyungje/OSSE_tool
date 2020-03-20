# OSSE_aerosol/README_Codes.md

Manual to run the simulation of OSSE_aerosol code for information contents analysis of aerosol vertical profiling

## 1. Order of codes
- osse_aerosol_0_master.pro  
This code summons below "1-3" codes.  
  - osse_aerosol_1_cal_pthg.pro
    - sub_module_setup_pthg_clars_new_multiple.pro
    - sub_make_xsec_lut_using_absco_v5_mchoi_o2ab_new_multiple.pro
  - osse_aerosol_2_cal_aodprof.pro
    - sub_integrate_result_svo_lab_loop.pro
  - osse_aerosol_3_cal_rtm.pro
    - sub_integrate_result_svo_lab_loop.pro  
    
  
