# README: osse_aerosol_0_master.pro

## 1. 

### 1. 
- l0: 'i_loop_start' will be 0, 1, 2, .., etc. The code aims to calculate multiple aerosol scenarios by using multiple nodes of machine (i.e., isoprene tb1 to tb 24). the number of loop will be defined later. 

- l12: i_case = 0 to 3. There are each CLARS measurement cases of four different seasons and sites. (It will determine measurement time, solar zenith/azimuth angle; viewing zenith/azimuth angle if we want to simulate CLARS-FTS LABS measurement mode, not satellite).  

- l27-29: inputs of aerosol optical depth (AOD), peak height (PH), and half width (HW; half width half maximum). These variables can contain multiple values.  

- l32-33: inputs of viewing zenith angle (VZA) and Relative azimuth angle (RAA) for satellite measurement scenarios.  

- l52: the number of scenarios we want to calculate.  

- l60: 'input_loop_begin' = 0; default setting of 'zero' if we want to put the first scenario set (i.e., 'i_loop_start' of 0) to RTM000.  

- l66: 'ii_start' is an initial case to simulate among n_ii cases. Generally it is '0'.

- l67-68: 'n_core': the number of nodes. 'nwn_cal_bin':  the number of cases for each nodes. The number of scenarios (or cases) are divided by 'n_core'. E.g., when the total cases are 100 and the 'n_core' is 5, then five codes will be run simultaneously from different nodes to calculate each 20 cases.

- l71-78: 'ii0_set', 'ii1_set': the beginning and final cases for each nodes.

- l91: 'i_loop_start' is multipled by 10. Therefore, the 2nd set of cases will be calculated by using RTM010~, and the 3rd set of cases will be run by using RTM020~.

- l150: 'n_loop': the number of runs within each node (from 1 to 10). If 'n_loop' is 5, 'i_loop_start' is 1 (-> RTM starts from 010), then the RTMs 010, 011, 012, 013, 014 are used.

- l151: 'w_interval_read': spectral resolution of high resolution calculation. This is set as 0.02 cm-1 for the OSSEs, but can be changed. It changes the calculation time.

- l152: 'i_prof': '0' is for aerosol bulk Jacobian (including Gaussian shape vertical distribution parameters; peak height and half width half maximum); '1' is for AOD profile Jacobian (i.e., Jaocibians for layer by layer AOD).

- l178: 'workspace_pwd': location of the highest directory for this work ("OSSE_tool" folder)

- l201-: Read CLARS-FTS measurement log file and extract location information 
  - i.e., SNR, Solar zenith angle (SZA), Solar azimuth angle (SAA), Time, surface pressure and temperature at CLARS location
  - VZA and viewing azimuth angle are provided, but we don't use them for satellite but use the values defined at l32-33.
  
- l288: A default CLARS target is 'WestPasadena'. Other site can be used. The name of sites are based on the CLARS-FTS measurement log file.

- l295-308: defining location for output results

- l323-362: save setting information of geometry, time, and location.

- l369-374: declare 'Temperature shift' value as 0.01 K 

- l


