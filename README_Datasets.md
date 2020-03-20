# OSSE_aerosol/REAEME_Datasets.md

Required datasets to simulate OSSE for aerosol vertical profiling

### 1. AERONET
* Ground-based AERONET: Aerosol microphysical and optical seasonal climatology over Los Angeles (LA) basin area
* Location: ./Datasets/data_AERONET_LA_for_model/save_L20_results_20191023/save_AeroBulk_season.xdr
* Temporal period of data: 1998-2018
* The name of 6 sites: "MISR-JPL", "UCLA", "El_Segundo", "CalTech", "Santa_Monica_Colg", "Mount_Wilson"
* The version of data: AERONET Version 3 Level 2 all-points inversion data
* Resource/reference: https://aeronet.gsfc.nasa.gov/  
  
### 2. MiniMPL
* Ground-based MiniMPL: Aerosol vertical distribution seasonal climatology
* Location: ./Datasets/data_megacities/save_results_20191023/save_AeroBulk_peak_halfwidth_minimpl_season.xdr
* Temporal period of raw data: 2012-2016
* The name of site: "CalTech"
* Resource/reference: https://megacities.jpl.nasa.gov/  

### 3. MODIS MAIAC
* Satellite-based MODIS MAIAC: Surface reflectance BRDF parameters monthly climatology
* Location: ./Datasets/data_MAIAC_BRDF/save_MODIS_BRDF_Kernel_20190529_O2ABD.txt
* Temporal period of the original L2 data: 2000-2019
* Spatial resolution of the original L2 data: 1 km × 1 km
* Monthly BRDF parameters (isotropic, geometric, and volumetric) from multi-year period are calculated for 34 CLARS-FTS targets.
* Interpolation is conducted from MODIS channels to O2 A, B, and 1∆ bands using SCIAMACHY surface LER products
* Resource/reference (MODIS MAIAC): https://earthdata.nasa.gov/ 
* Resource/reference (SCIAMACHY surface LER): http://www.temis.nl/surface/scia_ler.html  

### 4. SCIAMACHY
* Satellite-based SCIAMACHY surface LER: Lambertian surface reflectance monthly climatology (back-up)
* Location: ./Datasets/data_SCIAMACHY/save_SCIAMACHY_LER_min_O2A_O2B_O2delta.txt  
* Spatial resolution of the original L2 data: 0.5° × 0.5°
* This can be used when lambertian surface is only assumed.
* Resource/reference: http://www.temis.nl/surface/scia_ler.html  

### 5. POLDER-GRASP
* Satellite-based POLDER/GRASP: Surface BRDF/BPDF parameters monthly climatology (back-up)
* Location: ./Datasets/data_GRASP_POLDER/POLDER_BRDF_BPDF_Kernels_for_O2_ABD.txt
* Tempoarl period of the original L3 data: 2008-2013
* Spatial resolution of the original L2 data: 7 km × 6 km
* Spatial average is conducted within CLARS-FTS target area
* Interpolation is conducted from POLDER channels to O2 A, B, and 1∆ bands using SCIAMACHY surface LER products
* Resource/reference: https://www.grasp-open.com/products/polder-data-release/  

### 6. MERRA2
* Meteorological field reanalysis dataset: Pressure, Temparature, H2O vertical profiles
* Location: ./Datasets/data_PTHG/data_MERRA2_dn_20191003
* Temporal resolution: 3-hour (instantaneous data are interpolated to CLARS-FTS measurement time)
* Spatial resolution: 0.625° × 0.5°
* Spatial average is conducted within CLARS-FTS target area
* Resource/reference: https://earthdata.nasa.gov/ 
* MERRA-2 inst3_3d_asm_Nv: 3d,3-Hourly,Instantaneous,Model-Level,Assimilation,Assimilated Meteorological Fields 0.625 x 0.5 degree V5.12.4 (M2I3NVASM) at GES DISC  

### 7. CLARS-FTS
* Target infomation file: longitude, latitude, viewing zenith angle, viewing azimuth angle
* Location: ./Datasets/data_CLARS_FTS/save_target_info/save_target_info.xdr (and ~.txt)  
* Log-file: measurement time, sun geometry, pressure and temperature at Mt. Wilson(CLARS-FTS; 1.67 km a.s.l)
* Location: ./Datasets/data_CLARS_FTS/log4spectra_XB_4_Tom/2019-07-01.txt  
* Resource/reference: contact to authors or https://megacities.jpl.nasa.gov/

### 8. Other PTHG profiles
* US standard 1976 Pressure-Temperature-Height profile
* Location: ./Datasets/data_PTHG/data_originally_used/PTH_us_standard_1976_res100m_max.dat
* Discover-AQ in-situ profile: Pressure, temperature, H2O, O3, NO2, HCHO, SO2, and CO
* Location: ./Datasets/data_PTHG/data_originally_used/DiscoverAQ_profile_1.dat
* Temporary profile for CO2 retrieval: Pressure, temperature, H2O, CO2, O2, CH4
* Location: ./Datasets/data_PTHG/data_originally_used/Initial_NSW_atmosphere_12aug13.dat

### 9. Solar spectrum
* SOLSPEC solar continuum spectrum + Geoffrey Toon's solar absorption spectrum
* Location: ./Datasets/data_solar_spectrum/SOLSPEC_Toon_convolved/Solar_SOLSPEC_TOON_for_O2AB_07000_15000_mchoi.xdr (and ~.txt)
* Resource/reference (Geoff Toon): https://mark4sun.jpl.nasa.gov/
* Resource/reference (SOLSPEC): http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+A/611/A1

### 10. Instrument line shape function
* (1) Gaussian function, (2) Sinc + boxcar function
* 0.02 nm spectral sampling interval
* Location: ./Datasets/data_ILS_0p02/n3001/  


