Pro osse_aerosol_0_master, i_loop_start;, i_case, i_loop_start;, n_loop, wn_interval_read
 
switch_justcal = 2 ;---keep as "2"
switch_constraint_aerosol = 0 ;--keep as "0"
 ; 0=starting from a priori (all pars are different w/ true values) ;--used for OSSE information contents analysis
 ; 1=starting from true value for aerosol (identical aerosol pars except for AOD, PH, HW; all different other pars)
 ; 2=starting from all true values (true-in-ture-out test for all pars!)

switch_perturb = 0 ;--keep as "0"
;--
;----------------------case selection: default is i_case of 0
i_case = 0 ;--2019/07/05 "WestPasadena" case (cl201907052231West.143)

if i_case eq 0 then begin
yy_read = 2019 & mm_read = 07 & dd_read = 05 & loc_select_idx_no = 3 ;cl201907052231West.143
endif
if i_case eq 1 then begin
yy_read = 2019 & mm_read = 04 & dd_read = 07 & loc_select_idx_no = 1 ;cl201904072232West.047
endif
if i_case eq 2 then begin
yy_read = 2019 & mm_read = 01 & dd_read = 04 & loc_select_idx_no = 33 ;cl201901042236West.141
endif
if i_case eq 3 then begin
yy_read = 2018 & mm_read = 10 & dd_read = 01 & loc_select_idx_no = 3 ; cl201810012156West.147
endif

AOD_set = [0.5];[0.1, 0.3, 0.5, 1.0] ;[0.3, 1.0];[0.1, 0.3, 0.5, 1.0]
PH_set  = [1.0];[0.2, 0.6, 1.0, 1.5, 2.0];[1.0, 4.0, 7.0, 10.0]
HW_set  = [0.3];[0.2, 0.6, 1.0, 1.5, 2.0];[1.0, 2.0, 3.0, 4.0]

;;--for nadir
VZA_set2 = [0.01]
RAA_set2 = [30.0]

;--SPEXone angles except for nadir
; VZA_set2 =  [20, 57];[20.0,+29.0,+48.0,+59.0,+66.0];[57];[0.01];[20.0];[0.01,+29.0,+48.0,+59.0,+66.0] ;total 9 angles
; RAA_set2 =  [30.0, 150.0];[0.01, 179.99, 30.0, 150.0, 45.0, 135.0, 60.0, 120.0]
; VZA_set2 =  [0.01];[20.0,+29.0,+48.0,+59.0,+66.0];[57];[0.01];[20.0];[0.01,+29.0,+48.0,+59.0,+66.0] ;total 9 angles
; RAA_set2 =  [30.0];[30.0, 150.0];[0.01, 179.99, 30.0, 150.0, 45.0, 135.0, 60.0, 120.0]

; ;--MAIA angles except for nadir
; VZA_set2 =  [20.0,+29.0,+48.0,+59.0,+66.0]
; RAA_set2 =  [30.0, 150.0]

;---when other angles/refractive indices are fixed.
; SZA_set2 = [30.0]
; ; REFR_set2 = ;range = 1.35-1.55
; REFI_set2 = [0.0005, 0.001, 0.005, 0.01];;range = 0.0005 - 0.0225 
; nFMF_set2 = [0.996, 0.997, 0.998, 0.999];0.99956 - 0.9998
; ; vFMF_set2 = ;0.287 - 0.951

n_ii = n_elements(AOD_set)*n_elements(PH_set)*n_elements(HW_set)$
      *n_elements(VZA_set2)*n_elements(RAA_set2);*n_elements(SZA_set2)$
      ;*n_elements(REFI_set2)*n_elements(nFMF_set2)

; n_ii = 50

;--for checking 

input_loop_begin = 0 ;---starting point of RTM number
;;i_loop_start = starting point of RTM for each loop (ii0~ii1)
;;i_loop = each RTM node for each ii loop

;--final RTM location: i_loop+i_loop_start+input_loop_begin

ii_start = 0  ;initial of n_ii list
n_core = 1;20 ;-the number of simultaneously calculatable scenarios. using tb1-tb24 (max is 24 on isoprene machine)
nwn_cal_bin = round((n_ii-ii_start)/n_core)


ii0_set = fltarr(n_core)*!values.f_nan
ii1_set = fltarr(n_core)*!values.f_nan 
for i_core=0L, n_core-1 do begin 
ii0_set(i_core) = ii_start + (i_core)*nwn_cal_bin*1
ii1_set(i_core) = ii_start + double(i_core+1)*double(nwn_cal_bin)*1 - 1
if i_core eq n_core-1 then ii1_set(i_core)= n_ii
print, i_core, ii0_set(i_core), ii1_set(i_core), ii1_set(i_core)-ii0_set(i_core)+1
endfor

;;----original
ii0 = ii0_set(i_loop_start)
; ; ii0 = 1000
ii1 = ii1_set(i_loop_start)

; ;--only for AOD0.0 test
; ii0 = 0
; ii1 = 0

print, ii0, ii1

i_loop_start = i_loop_start*10 ;--for RTM code allocation (e.g., )

ii = 0 ;---sub loop of AOD, PH, HW

for i_n_ii = 0L, n_ii-1 do begin

; for i_SZA = 0L, 0 do begin;  n_elements(SZA_set2)-1 do begin
for i_RAA = 0L, n_elements(RAA_set2)-1 do begin
for i_VZA = 0L, n_elements(VZA_set2)-1 do begin
for i_HW = 0L, n_elements(HW_set)-1 do begin
for i_PH = 0L, n_elements(PH_set)-1 do begin
for i_AOD = 0L, n_elements(AOD_set)-1 do begin
; for i_REFI = 0L, 0 do begin;  n_elements(REFI_set2)-1 do begin
; for i_FMF = 0L, 0 do begin;  n_elements(nFMF_set2)-1 do begin


if ii ge ii0 and ii le ii1 then begin ;--ii0/ii1 are different for each node (tb1-tb24)

i_combo_ini = 0
i_combo_fin = 0

for i_combo = i_combo_ini, i_combo_fin do begin


SysB_I = 1.00d ; meaningless here; - 0.04d ;1.00+0.04; or 1.00-0.04

; if i_combo eq i_combo_ini then switch_Syn_noiseadd = 1 ;or 2
; if i_combo ne i_combo_ini then switch_Syn_noiseadd = 0
; ; switch_Syn_noiseadd = 1 ; 1= random noise added, 2= random+sysnoise added, 0 = no noise
switch_Syn_cal = 0 ;--default = 1

if i_combo eq 0 then begin ; A-D, I-QU
switch_band_ABD = 2 ;--1=A-only, 2=A+B, 3=A+B+1∆ 
switch_Pol = 1; 0 = Radiance-only; 1 = I + pol; 2 = pol-only
switch_DoLP_Ipol = 2 ;fixed;;;;--0=DoLP, 1=Ipol, 2=Q and U
switch_inversion_first = 0 ;--default = 0
endif

; if i_combo eq 1 then begin ; A-D, I
; switch_band_ABD = 1
; switch_Pol = 0; 0 = Radiance-only; 1 = I + pol; 2 = pol-only
; switch_DoLP_Ipol = 2 ;fixed;;;;--0=DoLP, 1=Ipol, 2=Q and U
; switch_inversion_first = 1 ;---should be 1 for other three combonation tests
; endif

; if i_combo eq 2 then begin ; A, I-QU
; switch_band_ABD = 0
; switch_Pol = 1; 0 = Radiance-only; 1 = I + pol; 2 = pol-only
; switch_DoLP_Ipol = 2 ;fixed;;;;--0=DoLP, 1=Ipol, 2=Q and U
; switch_inversion_first = 1 ;---should be 1 for other three combonation tests
; endif

; if i_combo eq 3 then begin ; A, I
; switch_band_ABD = 0
; switch_Pol = 0; 0 = Radiance-only; 1 = I + pol; 2 = pol-only
; switch_DoLP_Ipol = 2 ;fixed;;;;--0=DoLP, 1=Ipol, 2=Q and U
; switch_inversion_first = 1 ;---should be 1 for other three combonation tests
; endif

n_loop = 10;6 ;8 ;--number of cpu to use within each node (e.g., tb1); physical no. of each node's cpu is 8; and calculatable no. of cpu is 16.
wn_interval_read = 0.1 ;--interval of spectrum (monochromatic spectra)
i_prof = 0 ;---Jacobian calculation for Aerosol Bulk (0) vs AOD profile (1)

nn = 0 ;meaningless, --number of iteration for convergence test..; default setting is "zero(0)" for 
sign_next_iter = 0 ;meaningless, ;--initialization of for convergence flag using "retrieval" mode; not used in DOFS calculation mode


; ;--target of synthetic data
AOD_syn = AOD_set(i_AOD)
PH_syn = PH_set(i_PH)
HW_syn = HW_set(i_HW)

; SNR_I  = 200 ;relative total radiance noise
; Noise_D  = 0.0050 ;absolute DoLP noise ;---not used here
; SNR_QU = SNR_I/sqrt(2.0)  ;relative to Q and U

i_step_ini = 0 ;meaningless in justcal mode
i_step_fin = 19 ;meaningless in justcal mode
 
str_AOD = strtrim(string(AOD_syn,f='(f4.2)'),2)
str_PH = strtrim(string(PH_syn,f='(f04.1)'),2)
str_HW = strtrim(string(HW_syn,f='(f3.1)'),2)

str_case = 'Case'+string(i_case,f='(i3.3)');'AOD'+str_AOD+'_PH'+str_PH+'_HW'+str_HW
str_ii = 'ii'+string(ii,f='(i4.4)')

i_step = 0 ; fixed
workspace_pwd = '/home/mchoi/OSSE_tool'

;--required additional parameters for CLARS-FTS measurement 
;--for each measurement's 
;--PointingAzi(->VAA), PointingDep(=90.0-VZA)
;--SZA, SAA
;--surface pressure and surface temperature
;--H2O profile

;--surface temperature should be from MERRA-2 also......
;--select specific data from log file first.
; loc_long_set = ['Spectralon','Fontana','RanchoCucamunga','Riverside','LakeMatt', $
;                 'Norco','Pomona','210Bend','Corona','NorthOC', $
;                 '60Industry','SantaFeDam','OCairPort','AngelsStadium','HuntinghtonBeach', $
;                 'LaMirada','605and60','SantaAnitaRace','LongBeach405','Downey', $
;                 'ELAwater','PalosVerdes','MarinaDelRey','DownTownFar','DownTownNear', $
;                 'WestPasadena','SantaMonicaMountains','Glendale','UniversalCityNew','ShermanOaks', $
;                 'WoodlandHills','VanNuysAirPort','CanogaPark', 'Northridge']
;---readcol site information
restore, workspace_pwd + '/Datasets/data_CLARS_FTS/save_target_info/save_target_info.xdr'
; save, loc_set, lat_set, lon_set, alt_set, VZA_set, VAA_set, $
;       filename = save_pwd +'/save_taraget_info.xdr'

CLARSFTS_pwd = workspace_pwd + '/Datasets/data_CLARS_FTS'

str_yy_read = string(yy_read,f='(i4.4)')
str_mm_read = string(mm_read,f='(i2.2)')
str_dd_read = string(dd_read,f='(i2.2)')

log_files = file_search(CLARSFTS_pwd+'/log4spectra_XB_4_Tom'+ $
                       '/'+str_yy_read+'-'+str_mm_read+'-'+str_dd_read+'.txt', count=n_log_files) ;read real CLARS-FTS measurement log

lines = file_lines(log_files(0))
log=strarr(1)
n_data=lines-1
data=strarr(n_data)

openr, 1, log_files
readf, 1, log
readf, 1, data
close, 1

log_dum = strsplit(log(0), string(9B), /extract) ;--Tab = String(9B)
filename_idx = where(log_Dum eq 'File', n_filename_idx)
loc_idx = where(log_dum eq 'TaNa', n_loc_idx)
; lat_idx = where(log_dum eq 'TLat', n_lat_idx)
; lon_idx = where(log_dum eq 'TLon', n_lon_idx)
; alt_idx = where(log_dum eq 'TAlt', n_alt_idx)
; PointingAzi_idx = where(log_dum eq 'PointingAzi', n_PointingAzi_idx)
; PointingDep_idx = where(log_dum eq 'PointingDep', n_PointingDep_idx)
SNR_ch1_idx = where(log_dum eq 'SNR_Chan1',n_SNR_ch1_idx)
SNR_ch2_idx = where(log_dum eq 'SNR_Chan2',n_SNR_ch1_idx)
SZA_idx = where(log_dum eq 'SunZenith', n_SZA_idx)
SAA_idx = where(log_dum eq 'SunAzimuth', n_SAA_idx)
hour_utc_idx = where(log_dum eq 'Time', n_hour_UTC_idx)

Pres_CLARS_idx = where(log_dum eq 'OuPr', n_Pres_CLARS_idx)
Temp_CLARS_idx = where(log_dum eq 'OuTe', n_Temp_CLARS_idx) ; fahrenheit

files_all = strarr(n_data)
loc_all = strarr(n_data)
SNR_ch1_all = fltarr(n_data) * !values.f_nan
SNR_ch2_all = fltarr(n_data) * !values.f_nan
SZA_all = fltarr(n_data) * !values.f_nan
SAA_all = fltarr(n_data) * !values.f_nan
RAA_all = fltarr(n_data) * !values.f_nan
VAA_all = fltarr(n_data) * !values.f_nan
VZA_all = fltarr(n_data) * !values.f_nan
loc_idx_all = fltarr(n_data) * !values.f_nan
Pres_CLARS_all = fltarr(n_data) * !values.f_nan
Temp_CLARS_all = fltarr(n_data) * !values.f_nan
hour_utc_all = fltarr(n_data) * !values.f_nan


for l=0L, n_data-1 do begin
  temp = strsplit(data(l), string(9B), /extract) ;---separated by tab =string(9B)
  
  files_all(l) = temp(filename_idx(0))
  loc_all(l) = temp(loc_idx(0))
  loc_idx_all(l) = where(loc_set eq loc_all(l))
  
  ; ;;;--exact CLARS geometry
  ; SZA_all(l) = temp(SZA_idx(0))
  ; SAA_all(l) = temp(SAA_idx(0))
  ; VZA_all(l) = VZA_set(loc_idx_all(l))
  ; VAA_all(l) = VAA_set(loc_idx_all(l))
  ; RAA_all(l) = abs(180.0 - abs(SAA_all(l) - VAA_all(l)))
  
  ;---SZA and SAA from CLARS, VZA and RAA/VAA from satellite simulation
  SZA_all(l) = temp(SZA_idx(0));SZA_set2(i_SZA)
  SAA_all(l) = temp(SAA_idx(0))
  VZA_all(l) = VZA_set2(i_VZA) ;VZA_set(loc_idx_all(l))
  RAA_all(l) = RAA_set2(i_RAA)  ;abs(180.0 - abs(SAA_all(l) - VAA_all(l)))
  VAA_all(l) = SAA_all(l) - RAA_all(l);VAA_set(loc_idx_all(l))
  

  hour_utc_all(l) = temp(hour_utc_idx(0)) ;---local time = UTC - 8

  Pres_CLARS_all(l) = temp(Pres_CLARS_idx(0)) ; hPa
  ; Temp_CLARS_all(l) = temp(Temp_CLARS_idx(0)) ; F
  Temp_CLARS_all(l) = (temp(Temp_CLARS_idx(0))-32.0)/1.8 + 273.15 ;--F -> C -> Kelvin
  ; Surf_temp_F=62.2
  ; Surf_temp_C=(Surf_temp_F-32.0)/1.8
  ; Surf_temp_K=Surf_temp_C+273.15

  SNR_ch1_all(l) = temp(SNR_ch1_idx(0))
  SNR_ch2_all(l) = temp(SNR_ch2_idx(0))
endfor

;--select first West_Pasadena data
loc_select_idx = where(loc_all eq 'WestPasadena', n_loc_Select_idx)
loc_set_idx = where(loc_set eq 'WestPasadena', n_loc_set_idx)

idx_dat = loc_select_idx(loc_select_idx_no)
idx_loc = loc_set_idx(0)
temp = strsplit(files_all(idx_dat),'.',/extract)

if switch_justcal eq 0 then savefolder_fin = string('Results_ret_justcal_'+temp(0))+'_'+str_ii
if switch_justcal eq 2 then savefolder_fin = string('Results_justcal_'+temp(0))+'_'+str_ii

savefolder_num = strmid(temp(0),2,12)

Jump96: ;;-for final AOD profile calculation when sign_next_iter eq 2

Result_pwd0 = workspace_pwd+'/Results/'+savefolder_fin

if (i_step eq 0) and (i_combo eq i_combo_ini) and (switch_inversion_first eq 0) then spawn, 'rm -rf '+ Result_pwd0+' '+Result_pwd0+'_temp'

file_mkdir, Result_pwd0

if switch_constraint_aerosol le 1 then savefolder_fin1 = string('All_Results_justcal') ;always 0 in here.
; if switch_constraint_aerosol le 1 then savefolder_fin1 = string('All_Results_justcal'+temp(0)) ;always 0 in here.
; if switch_constraint_aerosol eq 2 then savefolder_fin1 = string('All_Results_TITO_'+temp(0))
Result_pwd1 = workspace_pwd+'/Results/'+savefolder_fin1

; str_sub = str_case+'_'+temp(0)+'/'+str_ii;+'_AOD'+str_AOD+'_Peak'+str_PH+'_Width'+str_HW
str_sub = str_case+'_'+temp(0)+'/'+str_ii+'_AOD'+str_AOD+'_Peak'+str_PH+'_Width'+str_HW
str_sub1 = '_AOD'+str_AOD+'_Peak'+str_PH+'_Width'+str_HW
str_step = string(i_step, f='(i2.2)')

;--step 1) PTHG profile first
switch_PTHG_cal = 1
if switch_PTHG_cal eq 1 and sign_next_iter ne 2 then begin


  if i_step eq 0 then begin
  state_vectors_pwd = Result_pwd0+'/state_vectors_step'
  file_mkdir, state_vectors_pwd
  
  SZA_fin = SZA_all(idx_dat)
  VZA_fin = VZA_all(idx_dat)
  SAA_fin = SAA_all(idx_dat)
  VAA_fin = VAA_all(idx_dat)
  RAA_fin = RAA_all(idx_dat) ;abs(180.0 - abs(SAA_fin - VAA_fin))
  SCA_fin = !radeg*acos(-cos(!dtor*SZA_fin)*cos(!dtor*VZA_fin)+sin(!dtor*SZA_fin)*sin(!dtor*VZA_fin)*cos(!dtor*RAA_fin))

  str_SZA = strtrim(string(SZA_fin,f='(f05.2)'),2)
  str_VZA = strtrim(string(VZA_fin,f='(f05.2)'),2)
  str_RAA = strtrim(string(RAA_fin,f='(f06.2)'),2)
  str_SCA = strtrim(string(SCA_fin,f='(f06.2)'),2)

  str_sub2 = 'SZA'+str_SZA+'_VZA'+str_VZA+'_RAA'+str_RAA+'_SCA'+str_SCA

  openw, 1, state_vectors_pwd+'/Info_geometry_SZA_VZA_RAA_SCA.txt'
  printf, 1, SZA_fin
  printf, 1, VZA_fin
  printf, 1, SAA_fin
  printf, 1, VAA_fin
  printf, 1, RAA_fin
  printf, 1, SCA_fin
  close, 1

  openw, 1, state_vectors_pwd+'/Info_date_time_YYYY_MM_DD_HH_utc.txt'
  printf, 1, yy_read
  printf, 1, mm_read
  printf, 1, dd_read
  printf, 1, hour_utc_all(idx_dat)
  close, 1

  openw, 1, state_vectors_pwd+'/Info_location_Lon_Lat_Alt.txt'
  printf, 1, lon_set(idx_loc)
  printf, 1, lat_set(idx_loc)
  printf, 1, alt_set(idx_loc)
  close, 1
  endif;if i_step eq 0 then begin


; endif;if n_file eq 0 then begin
; endif
endif ;if switch_PTHG_cal eq 1 then begin

if i_step eq 0 then begin ;--use a-priori value as input
  ;---Tshift
  openw, 1, state_vectors_pwd+'/state_vectors_tshift_step'+str_step+'.txt'
  printf, 1, 0.01d ;--not set as zero in any condition; forced to 0.01 or -0.01
  close, 1
endif

files = file_search(Result_pwd1+'/*'+str_case+'_'+temp(0)+'/ii*'+str_sub1+'_'+str_sub2, count=n_files)

if n_files ge 1 then begin
  spawn, 'rm -rf '+Result_pwd0 ;-delete and skip
  goto, Jump80
endif

if n_files eq 0 then begin

;----Step 1) Pressure-Tempareture-Height-Gases (PTHG) profile preparation 
if ii eq 0 then begin ;---calculation once at first scenario (and use it after following calculation)
  ; if ii eq ii0 then begin ;---temporary
  ; if i_loop_start eq 0 and ii eq 0 then begin ;---original
  print, 'step 1) PTHG profile calculation -- start'
  osse_aerosol_1_cal_pthg, workspace_pwd, i_step, i_loop_start, input_loop_begin, $
                                            lon_set(idx_loc), lat_set(idx_loc), alt_set(idx_loc), $
                                            Pres_CLARS_all(idx_dat), Temp_CLARS_all(idx_dat), $
                                            yy_read, mm_read, dd_read, hour_utc_all(idx_dat), $
                                            SZA_all(idx_dat), VZA_all(idx_dat), SAA_all(idx_dat), VAA_all(idx_dat), $
                                            i_case, savefolder_fin, savefolder_num, $
                                            switch_perturb, switch_constraint_aerosol,ii, wn_interval_read
  print, 'step 1) PTHG profile calculation -- finish'
endif ;if ii eq ii0 then begin


if ii ne 0 then begin ;---copy some results of first scenario and use it after following calculation
; if ii ne ii0 then begin ;---temporary
; if i_loop_start ne 0 or ii ne 0 then begin
  file_0 = file_search(Result_pwd1+'/*/*/state_vectors_step/state_vectors_surfpress_step00.txt')
  file_copy, file_0(0), state_vectors_pwd+'/state_vectors_surfpress_step00.txt', /overwrite

  file_0 = file_search(Result_pwd1+'/*/*/state_vectors_step/NSW_heightgrid_1.txt'); NSW_heightgrid_1.txt
  file_copy, file_0(0), state_vectors_pwd+'/NSW_heightgrid_1.txt', /overwrite

  file_0 = file_search(Result_pwd1+'/*/*/state_vectors_step/atmos_mchoi_step00.txt'); NSW_heightgrid_1.txt
  file_copy, file_0(0), state_vectors_pwd+'/atmos_mchoi_step00.txt', /overwrite

  ;---copy of NSW_heightgrid_1.dat file to each node
  ; workspace_pwd = workspace_pwd
  original_file = workspace_pwd+'/RTM/GEMSTOOL_physicsdata/PTH_PROFILES'+'/NSW_heightgrid_1_ii0000.dat'
  for input_loop = 0L, 9 do begin
      rtm_pwd = workspace_pwd+'/RTM/RTM_LAB/'+'GEMSTOOL_mchoi_LAB_loop'+string(input_loop+i_loop_start+input_loop_begin,f='(i3.3)')+'/GEMSTOOL_NSW_wrapper'
      copy_target_file = rtm_pwd + '/NSW_heightgrid_1.dat'
      FILE_COPY, original_file, copy_target_file, /overwrite
  endfor
  ;---copy to state vectors folder 
  copy_target_file = state_vectors_pwd+'/NSW_heightgrid_1.txt'
  FILE_COPY, original_file, copy_target_file, /overwrite

  ;---copy of Xsec data to each RTM nodes
  if i_loop_start ge 2 then begin
  for input_loop = 0L, 9 do begin
      RTM_pwd = workspace_pwd+'/RTM/RTM_LAB/'+'GEMSTOOL_mchoi_LAB_loop'+string(0+0+0,f='(i3.3)')
      original_folder = RTM_pwd+'/GEMSTOOL_physicsdata/LUT_Xsections_data/results_files_V500'
      
      ; input_loop=0
      RTM1_pwd = workspace_pwd+'/RTM/RTM_LAB/'+'GEMSTOOL_mchoi_LAB_loop'+string(input_loop+i_loop_start+input_loop_begin,f='(i3.3)')
      target_folder = RTM1_pwd+'/GEMSTOOL_physicsdata/LUT_Xsections_data/results_files_V500'
      ; spawn, 'rm -rf '+target_folder
      
      target_location = RTM1_pwd+'/GEMSTOOL_physicsdata/LUT_Xsections_data' ;/results_files_V500'
      FILE_COPY, original_folder, target_location, /overwrite, /RECURSIVE
  endfor
endif ;if i_loop_start ge 1 then begin


for input_loop = 0L, 9 do begin ;--copy PTHG informations to 2nd-10th RTM folders
;NSW_heightgrid_1 in PTH_PROFILES folder
file_title = 'NSW_heightgrid_1'+'_ii'+string(0,f='(i4.4)')+'.dat' 
original_file = workspace_pwd+'/RTM/GEMSTOOL_physicsdata/PTH_PROFILES'+'/'+file_title;'NSW_heightgrid_1.dat'
rtm_pwd = workspace_pwd+'/RTM/RTM_LAB/'+'GEMSTOOL_mchoi_LAB_loop'+string(input_loop+i_loop_start+input_loop_begin,f='(i3.3)')+'/GEMSTOOL_NSW_wrapper'
copy_target_file = rtm_pwd + '/NSW_heightgrid_1.dat'
FILE_COPY, original_file, copy_target_file, /overwrite

;atmos_mchoi in PTH_PROFILES folder and GAS_PROFILES folder
file_title = 'atmos_mchoi'+'_ii'+string(0,f='(i4.4)')+'.dat' 
original_file = workspace_pwd+'/RTM/GEMSTOOL_physicsdata/PTH_PROFILES'+'/'+file_title;'NSW_heightgrid_1.dat'
rtm1_pwd = workspace_pwd+'/RTM/RTM_LAB/'+'GEMSTOOL_mchoi_LAB_loop'+string(input_loop+i_loop_start+input_loop_begin,f='(i3.3)')+'/GEMSTOOL_physicsdata/PTH_PROFILES'
rtm2_pwd = workspace_pwd+'/RTM/RTM_LAB/'+'GEMSTOOL_mchoi_LAB_loop'+string(input_loop+i_loop_start+input_loop_begin,f='(i3.3)')+'/GEMSTOOL_physicsdata/GAS_PROFILES'

copy_target_file = rtm1_pwd + '/atmos_mchoi.dat'
FILE_COPY, original_file, copy_target_file, /overwrite

copy_target_file = rtm2_pwd + '/atmos_mchoi.dat'
FILE_COPY, original_file, copy_target_file, /overwrite
endfor

      
endif ;if i_loop_start ne 0 or ii ne 0 then begin



for i_band = 0, switch_band_ABD, 1 do begin; 

  ;--calculate n_loop
  if i_band eq 0 then begin ;--O2A
    ; original
    wn0_fin0 = 12930.00d ; 775.193
    wn1_fin0 = 13250.00d ; 754.716 

    ; ;continuum for DSCOVR
    ; wn0_fin0 = 12750.00d ; 775.193
    ; wn1_fin0 = 12930.00d ; 754.716 
  endif

  if i_band eq 1 then begin ;--O2B
    ; original
    wn0_fin0 = 14300.00d ; 699.301
    wn1_fin0 = 14600.00d ; 684.931

    ; ;continuum for DSCOVR
    ; wn0_fin0 = 14600.00d ; 699.301
    ; wn1_fin0 = 14800.00d ; 684.931
  endif
  if i_band eq 2 then begin ;--O2B
    wn0_fin0 = 7720.00d ; 1298.70 ;---current test
    wn1_fin0 = 8020.00d ; 1250.000 ;---current test
  endif


  ;--preparation of calculation
  str_band = string(i_band, f='(i2.2)')
  wn_interval=wn_interval_read;0.02d0
  wn0_fin = wn0_fin0
  wn1_fin = wn1_fin0
  str_wnrange_fin  = string(wn0_fin,f='(f08.2)')+'_'+string(wn1_fin,f='(f08.2)')
  str_wn0_fin = string(wn0_fin,f='(f09.3)')
  str_wn1_fin = string(wn1_fin,f='(f09.3)')
  str_wn1_fin_aodprof = string(wn0_fin+wn_interval,f='(f09.3)')
  str_wnrange_fin_aodprof = string(wn0_fin,f='(f08.2)')+'_'+string(wn0_fin+wn_interval,f='(f08.2)')
  str_wn_interval=string(wn_interval,f='(f07.3)')

  str_lon_target = strtrim(lon_set(idx_loc),2)
  str_lat_target = strtrim(lat_set(idx_loc),2)
  str_alt_target = strtrim(alt_set(idx_loc),2) ;but not used
  str_yy = strtrim(string(yy_read,f='(i4.4)'), 2)
  str_mm = strtrim(string(mm_read,f='(i2.2)'), 2)
  str_dd = strtrim(string(dd_read,f='(i2.2)'), 2)
  str_hh_utc = strtrim(hour_utc_all(idx_dat), 2)
  str_SZA = strtrim(SZA_all(idx_dat), 2)
  str_VZA = strtrim(VZA_all(idx_dat), 2)
  str_SAA = strtrim(SAA_all(idx_dat), 2)
  str_VAA = strtrim(VAA_all(idx_dat), 2)
  str_idx_loc = strtrim(idx_loc,2)

  ;---Dividing number of channels to calculate within each node (multi-cpu core calculation)
  ;---fixed n_loop and change nwn_cal_bin
  n_wn_total=(wn1_fin-wn0_fin)/wn_interval+1
  ;n_loop = 16 ->>>read as input
  if n_wn_total le 50 then n_loop = 1
  nwn_cal_bin = round(n_wn_total/n_loop)

  ;--for checking 
  wn0_loop = fltarr(n_loop)*!values.f_nan
  wn1_loop = fltarr(n_loop)*!values.f_nan 
  for i_loop=0L, n_loop-1 do begin 
  wn0_loop(i_loop) = wn0_fin + (i_loop)*nwn_cal_bin*wn_interval
  wn1_loop(i_loop) = wn0_fin + double(i_loop+1)*double(nwn_cal_bin)*wn_interval-0.5*wn_interval ;--safely
  if i_loop eq n_loop-1 then wn1_loop(i_loop)= wn1_fin+0.05*wn_interval
  print, i_loop, wn0_loop(i_loop), wn1_loop(i_loop),  fix(((wn1_loop(i_loop)-wn0_loop(i_loop)))/wn_interval) + 1 ,f='(3(f10.3),1x,1(f7.2))'
  endfor

  str_n_loop = string(n_loop,f='(i3.3)')
  str_sign_next_iter = string(sign_next_iter,f='(i1.1)')


  ;--step 2) LAB calculation but only for AOD profile setting; Calculated RTM results are meaningless. We only take the AOD profile results
  if i_band eq 0 then begin ;--only for A band; other bands will use the A band results
  print, 'step 2) LAB calculation only for AOD profile setting -- start'
  result_pwd = result_pwd0 + '/LAB_step'+str_step+'_aodprof_setup/'+str_wnrange_fin_aodprof

  str_nn = string(nn,f='(i2.2)')
  spawn, 'rm -rf '+result_pwd+'/*'

  osse_aerosol_2_cal_aodprof, $
                            workspace_pwd, i_band, i_Step, 0, i_loop_start, input_loop_begin, 1, $
                            wn0_fin, wn0_fin+wn_interval, wn_interval, $
                            str_lon_target,str_lat_target,$
                            str_yy, str_mm, str_dd, str_hh_utc, $
                            str_SZA, str_VZA, str_SAA, str_VAA, idx_loc, sign_next_iter, $
                            nn, i_case, switch_constraint_aerosol, switch_perturb, $
                            savefolder_num, switch_justcal, ii, I_STEP_FIN, $
                            AOD_syn, PH_syn, HW_syn;, REFI_set2(i_REFI), nFMF_set2(i_FMF)

  cd, workspace_pwd+'/Codes'

  ;---integration
  sub_integrate_result_svo_lab_loop, result_pwd
  print, 'step 2) LAB calculation only for AOD profile setting -- finish'

  check_file = file_search(Result_pwd0+'/LAB_step00_aodprof_setup/*/Stokes_IQU_DoLP_TOAUP.out', count = n_check_file)
  if n_check_file eq 0 then begin
    print, 'Problem in state_vectors. -> STOP' 
    stop
  endif

  endif ;if i_band eq 0 then begin


  ;--step 3) RTM calculation
  print, 'step 3) RTM calculation -- start'
  result_pwd = Result_pwd0+ '/LAB_step'+str_step+'/'+str_wnrange_fin
  spawn, 'rm -rf '+result_pwd+'/*'

  for i_loop=0L, n_loop-1 do begin
    str_loop = string(i_loop, f='(i3.3)')
    str_loop_start = string(i_loop_start, f='(i3.3)')
    str_i_prof = string(i_prof, f='(i1.1)')
    str_switch_justcal = string(switch_justcal,f='(i1.1)')
    str_rtm = string(i_loop+i_loop_start+input_loop_begin, f='(i3.3)')
    str_input_loop_begin = string(input_loop_begin,f='(i3.3)')

    file_mkdir, workspace_pwd+'/log_run_files'
    spawn, 'rm -rf '+workspace_pwd+'/log_run_files/log_lab_rtm'+str_rtm+'*'
    
    spawn, 'nohup nice idl -e osse_aerosol_3_cal_rtm,' + $
                            ; workspace_pwd+','+$
                            str_band+','+str_step+','+str_loop+','+str_loop_start+','+str_input_loop_begin+','+str_n_loop+','+$
                            str_wn0_fin+','+str_wn1_fin+','+str_wn_interval+','+str_i_prof+','+$
                            str_lon_target+','+str_lat_target+','+$
                            str_yy+','+str_mm+','+str_dd+','+str_hh_utc+','+$
                            str_SZA+','+str_VZA+','+str_SAA+','+str_VAA+','+str_idx_loc+','+$
                            savefolder_num+','+str_switch_justcal+','+string(ii,f='(i4.4)')+ $
                            '>& ../log_run_files/log_lab'+'_ii'+string(ii,f='(i4.4)')+'_rtm'+str_rtm+'.txt &' ;+'_b'+str_band+'_s'+str_step+


  endfor

  runtime_m = 0.0d0
  rumtime_m_add = 1.0d0/6

  JUMP1: 

  result_files = file_search(result_pwd+'/Loop_*/Stokes_IQU_DoLP_TOAUP.out', count = n_result_files)

  if n_result_files ne n_loop then begin

  str_runtime_m = string(runtime_m,f='(f5.1)')   
  print, string(n_result_files, f='(i3.3)'), ' of ', $
        string(n_loop, f='(i3.3)'), $
        '   run time: '+str_runtime_m
  runtime_m = runtime_m+rumtime_m_add
  wait, 10.0       

  goto, JUMP1
  endif

  if n_result_files eq n_loop then print, '       All loop calculations finish (LAB) !!'

  cd, workspace_pwd+'/Codes'
  ;---integration
  sub_integrate_result_svo_lab_loop, result_pwd

  print, 'step 3) RTM calculation -- finish'

  spawn, 'rm -rf '+workspace_pwd+'/log_run_files/log_lab'+'_ii'+string(ii,f='(i4.4)')+'_*.txt'

endfor ; for i_band = 0, switch_band_ABD, 1 do begin; 

file_mkdir, result_pwd1+'/'+str_sub+'_'+str_sub2;+'_'+str_sub3
spawn, 'rm -rf '+ result_pwd1+'/'+str_sub+'_'+str_sub2;+'_'+str_sub3
spawn, 'mv '+result_pwd0 + ' ' + result_pwd1+'/'+str_sub+'_'+str_sub2;+'_'+str_sub3

;----re-move from _temp
; file_copy, result_pwd0+'_temp',result_pwd0, /overwrite, /recursive

stop

endif ;if n_files eq 0 then begin

endfor ;for i_combo = 0L, 3 do begin

endif ;if ii ge ii0 and ii le ii1 then begin


Jump80:

ii = ii + 1 ;for next scenario calculation

; endfor ;for i_FMF = 0L, n_elements(nFMF_set2)-1 do begin
; endfor ;for i_REFI = 0L, n_elements(REFI_set2)-1 do begin

endfor ;for i_AOD = 0L, n_elements(AOD_set)-1 do begin
endfor ;for i_PH = 0L, n_elements(PH_set)-1 do begin
endfor ;for i_HW = 0L, n_elements(HW_set)-1 do begin

endfor ;for i_VZA = 0L, n_elements(VZA_set)-1 do begin
endfor ;for i_RAA = 0L, n_elements(RAA_set)-1 do begin
; endfor ;for i_SZA = 0L, n_elements(SZA_set2)-1 do begin

endfor ;for i_n_ii = 0L, n_ii-1 do begin


stop
end