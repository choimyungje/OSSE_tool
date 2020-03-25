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

- l150: 


