pro osse_aerosol_3_cal_rtm, i_band, i_step, input_loop, i_loop_start, input_loop_begin, n_loop, $
                                       wn0_read, wn1_read,wn_interval_read, i_prof, $
                                       lon_target, lat_target, $
                                       yy_read, mm_read, dd_read, hh_utc_read, $
                                       SZA_read, VZA_read, SAA_read, VAA_read, idx_loc, $
                                       savefolder_fin, switch_justcal, ii

str_ii = 'ii'+string(ii,f='(i4.4)')

;---0.01cm-1 resolution calculation
workspace_pwd = '/home/mchoi/OSSE_tool' ;should be revised for own pwd
rtm_pwd = workspace_pwd+'/RTM/RTM_LAB/'+'GEMSTOOL_mchoi_LAB_loop'+string(input_loop+i_loop_start+input_loop_begin,f='(i3.3)')+'/GEMSTOOL_NSW_wrapper'

if switch_justcal eq 0 then begin
result_pwd = file_search(workspace_pwd+'/Results/'+'Results_ret_justcal_cl'+strtrim(string(savefolder_fin),2)+'*_'+str_ii)
endif
; if switch_justcal eq 1 then begin
; result_pwd = file_search(workspace_pwd+'/Results/'+'Results_justcal_cl'+strtrim(string(savefolder_fin),2)+'*')
; endif
if switch_justcal eq 2 then begin
result_pwd = file_search(workspace_pwd+'/Results/'+'Results_justcal_cl'+strtrim(string(savefolder_fin),2)+'*_'+str_ii)
endif

result_pwd = result_pwd(0)

if i_prof eq 0 then begin ;--bulk option
; i_prof = 0 & loading_case_input = 3 & str_title = 'AerosolBulk';---bulk & gaussian
loading_case_input = 3 & str_title = 'AerosolBulk';---bulk & gaussian
endif

if i_prof eq 1 then begin ;--Prof option
; i_prof = 1 & loading_Case_input = 4 & str_title = 'AerosolProf' ;---AOD profile from file
loading_Case_input = 4 & str_title = 'AerosolProf' ;---AOD profile from file
endif


; ;---Define height node based on Pressure-Height grid
; ;---read PH file
; readcol, PH_file, level, Pres_read, Height_read, /silent
PH_file = workspace_pwd+'/RTM/GEMSTOOL_physicsdata/PTH_PROFILES'+'/NSW_heightgrid_1.dat'
readcol, PH_file, level, altitude_read, f=('D,D'),/silent
n_level=n_elements(level)

;---layer km setting
layer_thick_km = altitude_read(0:n_level-2)-altitude_read(1:n_level-1)
n_layer = n_elements(layer_thick_km)

;---read state vectors: aerosol
str_step = 'step'+string(i_step,f='(i2.2)')
str_step00 = 'step'+string(0,f='(i2.2)')
state_vectors_pwd = result_pwd+ '/state_vectors_step'
file_mkdir, state_vectors_pwd

          
;-----read state vectors
;--1) write first; afterthat
readcol, state_vectors_pwd+'/state_vectors_aerosol_'+str_step+'.txt', SV_aerosols, format='d',/silent
          REFR_1 = SV_aerosols(0)
          REFI_1 = SV_aerosols(1)
          PSD_rad_1 = SV_aerosols(2)
          PSD_sigma_1 = SV_aerosols(3)
          REFR_2 = SV_aerosols(4)
          REFI_2 = SV_aerosols(5)
          PSD_rad_2 = SV_aerosols(6)
          PSD_sigma_2 = SV_aerosols(7)
          FMF = SV_aerosols(8)
          AOD_fin = SV_aerosols(9)
          peak_height_read = SV_aerosols(10)
          half_width_read = SV_aerosols(11)

readcol, state_vectors_pwd+'/state_vectors_albedo_'+str_step00+'_O2A.txt', SV_albedo, /silent
         albedo_S_c0_O2A = SV_albedo(0)
         albedo_S_c1_O2A = SV_albedo(1)

readcol, state_vectors_pwd+'/state_vectors_albedo_'+str_step00+'_O2D.txt', SV_albedo, /silent
         albedo_S_c0_O2D = SV_albedo(0)
         albedo_S_c1_O2D = SV_albedo(1)

readcol, state_vectors_pwd+'/state_vectors_albedo_'+str_step00+'_O2B.txt', SV_albedo, /silent
         albedo_S_c0_O2B = SV_albedo(0)
         albedo_S_c1_O2B = SV_albedo(1)

readcol, state_vectors_pwd+'/state_vectors_surfpress_'+str_step+'.txt', SV_surfpress, /silent
         surfpress=SV_surfpress(0)

readcol, state_vectors_pwd+'/state_vectors_tshift_'+str_step+'.txt', SV_tshift, /silent
         tshift_read = SV_tshift(0)

readcol, state_vectors_pwd+'/state_vectors_h2oscaling_'+str_step+'.txt', SV_H2Oscaling, /silent
         H2Oscaling_read = SV_H2Oscaling(0)


readcol, state_vectors_pwd+'/state_vectors_brdf_'+str_step+'_O2A.txt', SV_brdf, /silent
K_iso_O2A = SV_brdf(0)
K_vol_O2A = SV_brdf(1)
K_geo_O2A = SV_brdf(2)
K_pol_O2A = SV_brdf(3)
Par1_geo_O2A = SV_brdf(4)
Par2_geo_O2A = SV_brdf(5)
Par1_pol_O2A = SV_brdf(6)
Par2_pol_O2A = SV_brdf(7)
Par3_pol_O2A = SV_brdf(8)

readcol, state_vectors_pwd+'/state_vectors_brdf_'+str_step+'_O2B.txt', SV_brdf, /silent
K_iso_O2B = SV_brdf(0)
K_vol_O2B = SV_brdf(1)
K_geo_O2B = SV_brdf(2)
K_pol_O2B = SV_brdf(3)
Par1_geo_O2B = SV_brdf(4)
Par2_geo_O2B = SV_brdf(5)
Par1_pol_O2B = SV_brdf(6)
Par2_pol_O2B = SV_brdf(7)
Par3_pol_O2B = SV_brdf(8)

readcol, state_vectors_pwd+'/state_vectors_brdf_'+str_step+'_O2D.txt', SV_brdf, /silent
K_iso_O2D = SV_brdf(0)
K_vol_O2D = SV_brdf(1)
K_geo_O2D = SV_brdf(2)
K_pol_O2D = SV_brdf(3)
Par1_geo_O2D = SV_brdf(4)
Par2_geo_O2D = SV_brdf(5)
Par1_pol_O2D = SV_brdf(6)
Par2_pol_O2D = SV_brdf(7)
Par3_pol_O2D = SV_brdf(8)

;-------aodprof is updated as copying file itself around line # 542.

;---save SV_all txt file 
openw, 1, state_vectors_pwd+'/state_vectors_all_ini_apriori_O2A.txt'
printf, 1, K_iso_O2A
printf, 1, K_vol_O2A
printf, 1, K_geo_O2A
for i=0L, n_elements(SV_aerosols)-1 do begin
printf, 1, SV_aerosols(i)
endfor
printf, 1, surfpress
printf, 1, H2Oscaling_read
printf, 1, tshift_read
close, 1

;---save SV_all txt file 
openw, 1, state_vectors_pwd+'/state_vectors_all_ini_apriori_O2B.txt'
printf, 1, K_iso_O2B
printf, 1, K_vol_O2B
printf, 1, K_geo_O2B
for i=0L, n_elements(SV_aerosols)-1 do begin
printf, 1, SV_aerosols(i)
endfor
printf, 1, surfpress
printf, 1, H2Oscaling_read
printf, 1, tshift_read
close, 1

;---save SV_all txt file 
openw, 1, state_vectors_pwd+'/state_vectors_all_ini_apriori_O2D.txt'
printf, 1, K_iso_O2D
printf, 1, K_vol_O2D
printf, 1, K_geo_O2D
for i=0L, n_elements(SV_aerosols)-1 do begin
printf, 1, SV_aerosols(i)
endfor
printf, 1, surfpress
printf, 1, H2Oscaling_read
printf, 1, tshift_read
close, 1

openw, 1, state_vectors_pwd+'/state_vectors_all_ini_apriori_O2AB.txt'
printf, 1, K_iso_O2A
printf, 1, K_vol_O2A
printf, 1, K_geo_O2A
printf, 1, K_iso_O2B
printf, 1, K_vol_O2B
printf, 1, K_geo_O2B
for i=0L, n_elements(SV_aerosols)-1 do begin
printf, 1, SV_aerosols(i)
endfor
printf, 1, surfpress
printf, 1, H2Oscaling_read
printf, 1, tshift_read
close, 1

openw, 1, state_vectors_pwd+'/state_vectors_all_ini_apriori_O2ABD.txt'
printf, 1, K_iso_O2A
printf, 1, K_vol_O2A
printf, 1, K_geo_O2A
printf, 1, K_iso_O2B
printf, 1, K_vol_O2B
printf, 1, K_geo_O2B
printf, 1, K_iso_O2D
printf, 1, K_vol_O2D
printf, 1, K_geo_O2D
for i=0L, n_elements(SV_aerosols)-1 do begin
printf, 1, SV_aerosols(i)
endfor
printf, 1, surfpress
printf, 1, H2Oscaling_read
printf, 1, tshift_read
close, 1

;----other fixed conditions
Lat_dum=[  34.221d, lat_target]
Lon_dum=[-118.0571, lon_target]
Year_dum=yy_read
Month_dum=mm_read
DayofMonth_dum=dd_read
Hour_dum = hh_utc_read
Minute_dum = 00
Second_dum = 00

SZA_meas = [SZA_read, SZA_read]
SAA_meas = [SAA_read, SAA_read]
VZA_meas = [68.633659, VZA_read]
VAA_meas = [8.8333292, VAA_read]
RAA_meas = abs(180.0 - abs(SAA_meas - VAA_meas))
Scat_meas = !radeg*acos(-cos(!dtor*SZA_meas)*cos(!dtor*VZA_meas)+sin(!dtor*SZA_meas)*sin(!dtor*VZA_meas)*cos(!dtor*RAA_meas))

level_node_km = [1.673+0.0001, 1.673]
idx_target = 1 ;--0:CLARS   1:LABS
solz_fin = SZA_meas(idx_target)
senz_fin = VZA_meas(idx_target)
rela_fin = RAA_meas(idx_target)
scat_fin = Scat_meas(idx_target)

lat_dum_fin = lat_Dum(idx_target)
lon_dum_fin = lon_Dum(idx_target)

level_node_km_fin = level_node_km(idx_target)

wl_AOD = 760.00 

if i_band eq 0 then begin ;--O2A
  wn0_fin0 = double(wn0_read)
  wn1_fin0 = double(wn1_read)
  K_iso_fin = K_iso_O2A
  K_geo_fin = K_geo_O2A
  K_vol_fin = K_vol_O2A
  K_pol_fin = K_pol_O2A
  Par1_geo_fin = Par1_geo_O2A
  Par2_geo_fin = Par2_geo_O2A 
  Par1_pol_fin = Par1_pol_O2A 
  Par2_pol_fin = Par2_pol_O2A 
  Par3_pol_fin = Par3_pol_O2A 
  albedo_c0 = albedo_S_c0_O2A 
  albedo_c1 = albedo_S_c1_O2A 
  i_Xsec = 0 ; ABSCO
endif

if i_band eq 1 then begin ;--O2B
  wn0_fin0 = double(wn0_read)
  wn1_fin0 = double(wn1_read)
  K_iso_fin = K_iso_O2B 
  K_geo_fin = K_geo_O2B 
  K_vol_fin = K_vol_O2B 
  K_pol_fin = K_pol_O2B 
  Par1_geo_fin = Par1_geo_O2B 
  Par2_geo_fin = Par2_geo_O2B 
  Par1_pol_fin = Par1_pol_O2B 
  Par2_pol_fin = Par2_pol_O2B 
  Par3_pol_fin = Par3_pol_O2B 
  albedo_c0 = albedo_S_c0_O2B 
  albedo_c1 = albedo_S_c1_O2B 
  i_Xsec = 0 ; ABSCO
endif

if i_band eq 2 then begin ;--O2-1∆
  wn0_fin0 = double(wn0_read)
  wn1_fin0 = double(wn1_read)
  K_iso_fin = K_iso_O2D  
  K_geo_fin = K_geo_O2D 
  K_vol_fin = K_vol_O2D  
  K_pol_fin = K_pol_O2D  
  Par1_geo_fin = Par1_geo_O2D 
  Par2_geo_fin = Par2_geo_O2D 
  Par1_pol_fin = Par1_pol_O2D 
  Par2_pol_fin = Par2_pol_O2D 
  Par3_pol_fin = Par3_pol_O2D 
  albedo_c0 = albedo_S_c0_O2D 
  albedo_c1 = albedo_S_c1_O2D 
  i_Xsec = 1 ; HITRAN
endif


wn0_fin = double(wn0_read)
wn1_fin = double(wn1_read)
wn_start_albedo = wn0_read 
wn_end_albedo = wn1_read 
wn_interval=wn_interval_read;0.02d0

n_wn_total=(wn1_fin-wn0_fin)/wn_interval+1
nwn_cal_bin = round(n_wn_total/n_loop)

if n_loop eq 0 then n_loop = 1

wn0_loop = fltarr(n_loop)*!values.f_nan
wn1_loop = fltarr(n_loop)*!values.f_nan 

for i_loop=0L, n_loop-1 do begin 
  wn0_loop(i_loop) = wn0_fin + (i_loop)*nwn_cal_bin*wn_interval
  wn1_loop(i_loop) = wn0_fin + double(i_loop+1)*double(nwn_cal_bin)*wn_interval-0.5*wn_interval ;--safely
  if i_loop eq n_loop-1 then wn1_loop(i_loop)= wn1_fin+0.05*wn_interval
  print, i_loop, wn0_loop(i_loop), wn1_loop(i_loop),  fix(((wn1_loop(i_loop)-wn0_loop(i_loop)))/wn_interval) + 1 ,f='(3(f10.3),1x,1(f7.2))'
endfor



AOD_node = [AOD_fin] 
Xsec_node = ['ABSCOv5','HITRAN2016']

level_node_km_set = level_node_km_fin 
level_node_idx = [interpol(Level, altitude_read, level_node_km_set,/quadratic)]

SZA = [solz_fin] 
VZA = [senz_Fin] 
RAA = [rela_fin] 
 
albedo = albedo_c0 
 
Peak_height_node = [peak_height_read] 
Half_width_node  = [half_width_read]

idx_Single_AOD_layer = [20] 
 
for i_level = 0, 0 do begin ;n_elements(level_node_idx)-1 do begin
for i_AOD    = 0, n_elements(AOD_node)-1 do begin
for i_peak   = 0L, n_elements(Peak_height_node)-1 do begin
for i_Half   = 0L, n_elements(Half_width_node)-1 do begin

 
for i_single_AOD_layer = 0L, 0 do begin ; new 06/17/19

;--------------- Aerosol Loading --------------------------
Top_aerosol_layer  = 15.00 
Bot_aerosol_layer  = altitude_read(n_level-1)+0.01 
Loading_case       = loading_Case_input    ;(1 = uniform, 2 = exponential, 3 = GDF, 4 = User-defined)
AOD                = AOD_node(i_AOD)
reference_wl       = wl_AOD  
relaxation_para    = 0.2  
peak_height        = Peak_height_node(i_peak) 
half_width         = Half_width_node(i_Half) 

 ;------------- Atmosphere Configure file control----------
if AOD_node(i_AOD) eq 0.0 then begin
;GEOMSTOOL_Atmosphere.cfg file
Trace_gases_switch = 'T'
Aerosol_cal_switch = 'F'
Mie_cal_switch = 'F'
Tmat_cal_swith = 'F'
Cloud_cal_swith = 'F'
;GEMSTOOL_LinControl.cfg
do_GasProfile_Jacobians      = 'F   ! GasProfile_Jacobians'
do_AerOpdepProfile_Jacobians = 'F   ! AODProfile_Jacobian'
do_AerBulk_Jacobians         = 'F   ! AerBulk_Jacobian'
do_Surface_Jacobians         = 'T   ! Surface_Jacobians'
do_Tshift_Jacobian           = 'F   ! Tshift_Jacobian'
do_SurfPress_Jacobian        = 'T   ! SurfPress_Jacobian'
;do_SIF_Jacobians             = 'F   ! SIF_Jacobians ! New 10/18/16, 11/30/16'
do_normalized_wfoutput       = 'F   ! normalized_wfoutput'
if i_Xsec eq 0 then do_hitran= 'F   ! hitran'
if i_Xsec eq 1 then do_hitran= 'T   ! hitran'
  ;--------------H2O scaling------------------------------------------
  do_H2OScaling_Jacobian       = 'F   ! H2OScaling_Jacobian'
  do_H2OScaling = 'T'
endif

if AOD_node(i_AOD) gt 0.0 then begin
if i_prof eq 0 then begin ;---Bulk ;--jacobian!
  ;GEOMSTOOL_Atmosphere.cfg file
  Trace_gases_switch = 'T'
  Aerosol_cal_switch = 'T'
  Mie_cal_switch = 'T'
  Tmat_cal_swith = 'F'
  Cloud_cal_swith = 'F'
  ;GEMSTOOL_LinControl.cfg
  do_GasProfile_Jacobians      = 'F   ! GasProfile_Jacobians'
  do_AerOpdepProfile_Jacobians = 'F   ! AODProfile_Jacobian'
  do_AerBulk_Jacobians         = 'T   ! AerBulk_Jacobian'
  do_Surface_Jacobians         = 'T   ! Surface_Jacobians'
  do_Tshift_Jacobian           = 'T   ! Tshift_Jacobian'
  do_SurfPress_Jacobian        = 'T   ! SurfPress_Jacobian'
  ;do_SIF_Jacobians             = 'F   ! SIF_Jacobians !  New 10/18/16, 11/30/16'
  do_normalized_wfoutput       = 'F   ! normalized_wfoutput'
  if i_Xsec eq 0 then do_hitran= 'F   ! hitran'
  if i_Xsec eq 1 then do_hitran= 'T   ! hitran'
  ;--------------H2O scaling------------------------------------------
  if i_band eq 0 or i_band eq 1 then begin
  do_H2OScaling_Jacobian       = 'T   ! H2OScaling_Jacobian'
  do_H2OScaling = 'T'
  endif
  if i_band eq 2 then begin ;---1∆-band
  do_H2OScaling_Jacobian       = 'F   ! H2OScaling_Jacobian'
  do_H2OScaling = 'T'
  endif
endif ;if i_prof eq 0 then begin

if i_prof eq 1 then begin
  ;GEOMSTOOL_Atmosphere.cfg file
  Trace_gases_switch = 'T'
  Aerosol_cal_switch = 'T'
  Mie_cal_switch = 'T'
  Tmat_cal_swith = 'F'
  Cloud_cal_swith = 'F'
  ;GEMSTOOL_LinControl.cfg
  do_GasProfile_Jacobians      = 'T   ! GasProfile_Jacobians'
  do_AerOpdepProfile_Jacobians = 'T   ! AODProfile_Jacobian'
  do_AerBulk_Jacobians         = 'F   ! AerBulk_Jacobian'
  do_Surface_Jacobians         = 'F   ! Surface_Jacobians'
  do_Tshift_Jacobian           = 'F   ! Tshift_Jacobian'
  do_SurfPress_Jacobian        = 'F   ! SurfPress_Jacobian'
  ;do_SIF_Jacobians             = 'F   ! SIF_Jacobians ! New 10/18/16, 11/30/16'
  do_normalized_wfoutput       = 'F   ! normalized_wfoutput'
  if i_Xsec eq 0 then do_hitran= 'F   ! hitran'
  if i_Xsec eq 1 then do_hitran= 'T   ! hitran'
  ;--------------H2O scaling------------------------------------------
  do_H2OScaling_Jacobian       = 'F   ! H2OScaling_Jacobian'
  do_H2OScaling = 'F'
endif ;if i_prof eq 1 then begin
endif ;if AOD_node(i_AOD) gt 0.0 then begin


;--------------LinControl_AerBulk control
do_AerBulk_LoadPars_Jacobians = 'T     ! Flag for AerBulk_LoadPars_Jacobians'
do_AerBulk_RefIndex_Jacobians = 'T     ! Flag for AerBulk_RefIndex_Jacobians'
do_AerBulk_ShapeFac_Jacobians = 'T     ! Flag for AerBulk_ShapeFac_Jacobians'
do_AerBulk_SizeDist_Jacobians = 'T     ! Flag for AerBulk_SizeDist_Jacobians'
do_AerBulk_BmodFrac_Jacobians = 'T     ! Flag for AerBulk_BmodFrac_Jacobians'


;--------------H2O scaling------------------------------------------
value_H2OScaling = H2Oscaling_read ;1.0d

;--------------CH4 scaling (not used);------------------------------
do_CH4Scaling_Jacobian       = 'F   ! CH4Scaling_Jacobian'
do_CH4Scaling = 'F' ;----not used
value_CH4Scaling = 1.0d ;---not used
;-------------Tshift setting read
Tshift_value = tshift_read;1.00d0 ;---Jacobian is not change almost..

;-------------- RT control --------------------------------
n_stream  = 4 ;should be revised
n_stokes  = 3
Upwelling = 'T'
;--here: output_level?
do_Path_Radiance           = 'T        ! turn on flag for "Path Radiance"      (TBD)'
do_Spherical_Albedo        = 'T        ! turn on flag for "Spherical Albedo"   (TBD)'
do_2_wayTransmittance      = 'T        ! turn on flag for "2-wayTransmittance" (TBD)'
do_new_FO_code             = 'F        ! turn on, using new FO code'
do_FO_code_regular_PS_mode = 'T        ! turn on, using FO code regular PS mode'

;-------------- SIF control --------------------------------
;--turn off for CLARS FTS test
DO_SIF = 'F'
DO_SIF_ExactOnly = 'T' ;not changed
DO_SIF_LinearParam = 'T' ;not changed
SIF755_Scaling_Constant = 1.0 ;default
SIF755_Scaling_Gradient = 26.97 ;default
SIF755_Amplitude = 1.0 ;default ;--not used when LinearParm is 'T' 
do_SIF_Jacobians = 'F'


for r = 0L, n_elements(albedo)-1 do begin
for l = 0L, n_elements(SZA)-1 do begin
for j = 0L, n_elements(VZA)-1 do begin
for k = 0L, n_elements(RAA)-1 do begin

for i_loop=input_loop, input_loop do begin ;--wn loop; ;originally from 0L

wn0 = wn0_loop(i_loop)
wn1 = wn1_loop(i_loop)

str_wn0 = string(wn0,f='(f9.3)')+'0000'
str_wn1 = string(wn1,f='(f9.3)')+'0000'

;============== Make Configuration files ========================================
  
  cd, rtm_pwd

file_copy, state_vectors_pwd+'/state_vectors_aodprof_'+str_step+'.txt', './aod_loading_user.input', /overwrite



;---other cfg files are in below pwd.  
cd, 'GEMSTOOL_NSW_Configfiles'

spawn, 'cp -f ./original_cfg/*.cfg ./'

;---------- TimePosition configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_TimePosition.cfg temp_TimePosition.cfg'
openr, 20, 'temp_TimePosition.cfg'
a=0L
openw, 10, 'GEMSTOOL_TimePosition.cfg'
while (not EOF(20)) do begin
  readf, 20, header
  if (a eq 0) then printf, 10, Lat_dum_fin,   format = '(f11.6)' else $
  if (a eq 1) then printf, 10, Lon_dum_fin,   format = '(f11.6)' else $
  if (a eq 3) then printf, 10, Year_dum,   format = '(i4.4)' else $
  if (a eq 4) then printf, 10, Month_dum,   format = '(i2.2)' else $
  if (a eq 5) then printf, 10, DayofMonth_dum,   format = '(i2.2)' else $
  if (a eq 6) then printf, 10, Hour_dum,   format = '(i2.2)' else $
  if (a eq 7) then printf, 10, Minute_dum,   format = '(i2.2)' else $
  if (a eq 8) then printf, 10, Second_dum,   format = '(i2.2)' else $

  printf, 10, header
  a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_TimePosition.cfg'

;---------- Wavenumber configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_Wavenums.cfg temp_Wavenums.cfg'
openr, 20, 'temp_Wavenums.cfg'
a=0L
openw, 10, 'GEMSTOOL_Wavenums.cfg'
while (not EOF(20)) do begin
  readf, 20, header
  ; if (a eq 0) then printf, 10, wn0,   format = '(f9.3)' else $
  ; if (a eq 1) then printf, 10, wn1,   format = '(f9.3)' else $
  if (a eq 0) then printf, 10, str_wn0,   format = '(a)' else $
  if (a eq 1) then printf, 10, str_wn1,   format = '(a)' else $
  if (a eq 2) then printf, 10, wn_interval,   format = '(f11.5)' else $
  
  printf, 10, header
  a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_Wavenums.cfg'

;---------- Atmosphere configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_Atmosphere.cfg temp_Atmosphere.cfg'
openr, 20, 'temp_Atmosphere.cfg'
a=0L
openw, 10, 'GEMSTOOL_Atmosphere.cfg'
while (not EOF(20)) do begin
  readf, 20, header
  if (a eq 0) then printf, 10, Trace_gases_switch,   format = '(a1)' else $
  if (a eq 1) then printf, 10, Aerosol_cal_switch,   format = '(a1)' else $
  if (a eq 2) then printf, 10, Mie_cal_switch,   format = '(a1)' else $
  if (a eq 3) then printf, 10, Tmat_cal_swith,   format = '(a1)' else $
  if (a eq 4) then printf, 10, Cloud_cal_swith,   format = '(a1)' else $

  printf, 10, header
  a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_Atmosphere.cfg'

;---------- LinControl configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_LinControl.cfg temp_LinControl.cfg'
openr, 20, 'temp_LinControl.cfg'
a=0L
openw, 10, 'GEMSTOOL_LinControl.cfg'
while (not EOF(20)) do begin
  readf, 20, header
  if (a eq 00)  then printf, 10, do_GasProfile_Jacobians,   format = '(a1)' else $
  if (a eq 01)  then printf, 10, do_AerOpdepProfile_Jacobians,   format = '(a1)' else $
  if (a eq 02)  then printf, 10, do_AerBulk_Jacobians,   format = '(a1)' else $
  if (a eq 03)  then printf, 10, do_Surface_Jacobians,   format = '(a1)' else $
  if (a eq 04)  then printf, 10, do_Tshift_Jacobian,   format = '(a1)' else $
  if (a eq 05)  then printf, 10, do_SurfPress_Jacobian,   format = '(a1)' else $
  if (a eq 06)  then printf, 10, do_H2OScaling_Jacobian,   format = '(a1)' else $
  if (a eq 07)  then printf, 10, do_CH4Scaling_Jacobian,   format = '(a1)' else $
  if (a eq 08)  then printf, 10, do_SIF_Jacobians,   format = '(a1)' else $
  if (a eq 09)  then printf, 10, do_normalized_wfoutput,   format = '(a1)' else $
  if (a eq 10)  then printf, 10, do_hitran,   format = '(a1)' else $

  printf, 10, header
  a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_LinControl.cfg'

;---------- LinControl_AerBulk configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_LinControl_AerBulk.cfg temp_LinControl_AerBulk.cfg'
openr, 20, 'temp_LinControl_AerBulk.cfg'
a=0L
openw, 10, 'GEMSTOOL_LinControl_AerBulk.cfg'
while (not EOF(20)) do begin
  readf, 20, header
  if (a eq 0) then printf, 10, do_AerBulk_LoadPars_Jacobians,   format = '(a)' else $
  if (a eq 1) then printf, 10, do_AerBulk_RefIndex_Jacobians,   format = '(a)' else $
  if (a eq 2) then printf, 10, do_AerBulk_ShapeFac_Jacobians,   format = '(a)' else $
  if (a eq 3) then printf, 10, do_AerBulk_SizeDist_Jacobians,   format = '(a)' else $
  if (a eq 4) then printf, 10, do_AerBulk_BmodFrac_Jacobians,   format = '(a)' else $

  printf, 10, header
  a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_LinControl_AerBulk.cfg'

;---------- Aerosol Loading configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_AerosolLoading.cfg temp_AerosolLoading.cfg'
openr, 20, 'temp_AerosolLoading.cfg'
a=0L
openw, 10, 'GEMSTOOL_AerosolLoading.cfg'
while (not EOF(20)) do begin
	readf, 20, header
	if (a eq 1) then printf, 10, Top_aerosol_layer,  format = '(f9.4)' else $
	if (a eq 2) then printf, 10, Bot_aerosol_layer,  format = '(f9.4)' else $
	if (a eq 5) then printf, 10, loading_case,       format = '(i1)'   else $
  if (a eq 6) then printf, 10, AOD,                format = '(f7.4)' else $
	if (a eq 7) then printf, 10, reference_wl,       format = '(f9.4)' else $
	if (a eq 8) then printf, 10, relaxation_para,    format = '(f3.1)' else $
	if (a eq 9) then printf, 10, peak_height,        format = '(f7.4)' else $
	if (a eq 10)then printf, 10, half_width,         format = '(f7.4)' else $ 
	printf, 10, header
	a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_AerosolLoading.cfg'

;---------- Trace gases configuration file ----------------------------------

if i_band eq 0 or i_band eq 1 then begin ;---O2 and H2O for O2A band

header=strarr(1)
spawn, 'cp GEMSTOOL_TraceGases.cfg temp_TraceGases.cfg'
openr, 20, 'temp_TraceGases.cfg'
a=0L
openw, 10, 'GEMSTOOL_TraceGases.cfg'
while (not EOF(20)) do begin
	readf, 20, header
	if (a eq 0) then printf, 10, 2,  format = '(i1)' else $   ;fixed
	if (a eq 1) then printf, 10, 'O2  T T uu.prf',  format = '(a)' else $
	if (a eq 2) then printf, 10, 'H2O T T uu.prf',  format = '(a)' else $
  if (a eq 3) then printf, 10, do_H2OScaling, value_H2OScaling, format = '(a1, 1x, f10.5)' else $
	if (a eq 4) then printf, 10, do_CH4Scaling, value_CH4Scaling, format = '(a1, 1x, f10.5)' else $
  if (a eq 5) then printf, 10, ' ' else $
	printf, 10, header
	a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_TraceGases.cfg'
endif 

;---------- Tshift configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_Tshift.cfg temp_Tshift.cfg'
openr, 20, 'temp_Tshift.cfg'
a=0L
openw, 10, 'GEMSTOOL_Tshift.cfg'
while (not EOF(20)) do begin
	readf, 20, header
	if (a eq 0) then printf, 10, Tshift_value,  format = '(f8.4)' else $   ;fixed
	printf, 10, header
	a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_Tshift.cfg'

;---------- Geometries configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_Geometries.cfg temp_Geometries.cfg'
openr, 20, 'temp_Geometries.cfg'
a=0L
openw, 10, 'GEMSTOOL_Geometries.cfg'
while (not EOF(20)) do begin
	readf, 20, header
	if (a eq 1)then printf, 10, SZA(l),VZA(j),RAA(k),        format = '(3(f6.2,3x))' else $
	printf, 10, header
	a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_Geometries.cfg'

;---------- RTControl configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_RTControl.cfg temp_RTControl.cfg'
openr, 20, 'temp_RTControl.cfg'
a=0L
openw, 10, 'GEMSTOOL_RTControl.cfg'
while (not EOF(20)) do begin
    readf, 20, header
    
; ;-------------- RT control --------------------------------
    if (a eq 0)then printf, 10, n_stream,        format = '(i2)' else $
    if (a eq 1)then printf, 10, n_stokes,        format = '(i1)' else $
    if (a eq 2)then printf, 10, upwelling,       format = '(a1)' else $
    ; if (a eq 3)then printf, 10, level_node_idx(i_level), format = '(f6.3)' else $ ;format = '(i2)' else $ ;
    if (a eq 3)then printf, 10, 0, format = '(f6.3)' else $ ;format = '(i2)' else $ ;
    if (a eq 4)then printf, 10, do_Path_Radiance,       format = '(a)' else $
    if (a eq 5)then printf, 10, do_Spherical_Albedo,       format = '(a)' else $
    if (a eq 6)then printf, 10, do_2_wayTransmittance,       format = '(a)' else $
    if (a eq 7)then printf, 10, do_new_FO_code,       format = '(a)' else $
    if (a eq 8)then printf, 10, do_FO_code_regular_PS_mode,       format = '(a)' else $
    
    printf, 10, header
    a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_RTControl.cfg'

;---------- Albedo configuration file ---------------------------------- ;---modified with slope
header=strarr(1)
spawn, 'cp GEMSTOOL_AlbedoClosure.cfg temp_AlbedoClosure.cfg'
openr, 20, 'temp_AlbedoClosure.cfg'
a=0L
openw, 10, 'GEMSTOOL_AlbedoClosure.cfg'
while (not EOF(20)) do begin
  readf, 20, header
  if (a eq 1)then printf, 10, wn_start_albedo, wn_end_albedo, 2,        format = '(f11.5,1x,f11.5,1x,i3)' else $
  if (a eq 2)then printf, 10, 1, albedo(r),        format = '(i3,f15.10)' else $
  ; if (a eq 3)then printf, 10, 1, albedo_c1(r),        format = '(i3,f15.10)' else $
  if (a eq 3)then printf, 10, 2, albedo_c1,        format = '(i3,f15.10)' else $
  printf, 10, header
  a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_AlbedoClosure.cfg'

;---------- BRDF configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_BRDF.cfg temp_BRDF.cfg'
openr, 20, 'temp_BRDF.cfg'
a=0L
openw, 10, 'GEMSTOOL_BRDF.cfg'
while (not EOF(20)) do begin
  readf, 20, header
  if (a eq 0)then printf, 10, 'T', format='(a)' else $
  if (a eq 1)then printf, 10,  3, format='(i1)' else $ ;--should be "4" when BPDF is used
  if (a eq 2)then printf, 10,  1, 0, K_iso_fin, 0.0d, 0.0d, 0.0d,       format='(i2,",",i2,",",f7.3,",",f4.1,",",f4.1,",",f4.1)' else $
  if (a eq 3)then printf, 10,  3, 0, K_vol_fin, 0.0d, 0.0d, 0.0d,       format='(i2,",",i2,",",f7.3,",",f4.1,",",f4.1,",",f4.1)' else $
  if (a eq 4)then printf, 10,  4, 2, K_geo_fin, Par1_geo_fin, Par2_geo_fin, 0.0d,       format='(i2,",",i2,",",f7.3,",",f4.1,",",f4.1,",",f4.1)' else $
  if (a eq 5)then printf, 10, 14, 3, K_pol_fin, Par1_pol_fin, Par2_pol_fin, 1.0d, format='(i2,",",i2,",",f7.3,",",f7.3,",",f7.3,",",f7.3)' else $
  
  printf, 10, header
  a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_BRDF.cfg'

;---------- SIF configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_SurfaceLeaving.cfg temp_SurfaceLeaving.cfg'
openr, 20, 'temp_SurfaceLeaving.cfg'
a=0L
openw, 10, 'GEMSTOOL_SurfaceLeaving.cfg'
while (not EOF(20)) do begin
  readf, 20, header
  if (a eq 3)then printf, 10, DO_SIF, format = '(a)' else $
  if (a eq 13)then printf, 10, DO_SIF_ExactOnly, format = '(a)' else $
  if (a eq 14)then printf, 10, DO_SIF_LinearParam, format = '(a)' else $
  if (a eq 15)then printf, 10, SIF755_Scaling_Constant, format = '(f6.3)' else $ 
  if (a eq 16)then printf, 10, SIF755_Scaling_Gradient, format = '(f6.3)' else $ 
  if (a eq 17)then printf, 10, SIF755_Amplitude, format = '(f6.3)' else $ 
  printf, 10, header
  a=a+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_SurfaceLeaving.cfg'

;--for Bi-modal case
;---------- Mie aerosol configuration file ----------------------------------
header=strarr(1)
spawn, 'cp GEMSTOOL_MieAerosols.cfg temp_MieAerosols.cfg'
openr, 20, 'temp_MieAerosols.cfg'
i=0L
openw, 10, 'GEMSTOOL_MieAerosols.cfg'
while (not EOF(20)) do begin
	readf, 20, header
  if (i eq 0) then printf, 10, 'T',                           format = '(a)'       else $
	if (i eq 1) then printf, 10, FMF,                           format = '(f8.6)'       else $
	if (i eq 3) then printf, 10, PSD_rad_1,    PSD_rad_2,       format = '(2(f5.3,2x))' else $
	if (i eq 4) then printf, 10, PSD_sigma_1,  PSD_sigma_2,     format = '(2(f5.3,2x))' else $
  if (i eq 13)then printf, 10, REFR_1,       REFR_2,          format = '(2(f5.3,2x))' else $
	if (i eq 14)then printf, 10, REFI_1,       REFI_2,          format = '(2(f7.5,2x))' else $
	printf, 10, header
	i=i+1L
endwhile
close, 10
close, 20
spawn, 'rm temp_MieAerosols.cfg'


cd, '..'

;----move output files to specific folder
str_AOD='AOD'+strtrim(string(AOD_node(i_AOD),f='(f7.4)'),2)
str_SZA='SZA'+strtrim(string(SZA(l),f='(i3.3)'),2)
str_VZA='VZA'+strtrim(string(VZA(j),f='(i3.3)'),2)
str_RAA='RAA'+strtrim(string(RAA(k),f='(i3.3)'),2)
str_alb='ALB'+strtrim(string(albedo(r),f='(f8.6)'),2)
str_alb_s='ALB_SCIAMACHY';+strtrim(string(albedo_c0,f='(f15.10)'),2)

str_level_idx='Level'+string(i_level,f='(i1.1)')+'_'+strtrim(string(level_node_idx(i_level),f='(f6.3)'),2)+'idx'
str_level_km='Level'+string(i_level,f='(i1.1)')+'_'+strtrim(string(level_node_km(i_level),f='(f6.3)'),2)+'km'

str_peak='Peak'+strtrim(string(peak_height,f='(f7.4)'),2)
str_Half='Half'+strtrim(string(half_width,f='(f7.4)'),2)

str_wnrange_fin  = string(wn0_fin,f='(f08.2)')+'_'+string(wn1_fin,f='(f08.2)')
str_wnrange = 'Loop_'+string(wn0,f='(f08.2)')+'_'+string(wn1,f='(f08.2)')

str_Xsec = Xsec_node(i_Xsec)

save_pwd=result_pwd+'/LAB_'+str_step+'/'+str_wnrange_fin+'/'+str_wnrange
file_mkdir, save_pwd
spawn, 'rm -rf '+save_pwd+'/*'



;--write input aerosol option as .dat txt file
openw, 1, save_pwd +'/Input_aerosol.dat'
printf, 1, REFR_1 ; 00_n_real_mode1
printf, 1, REFI_1 ; 01_n_imag_mode1
printf, 1, PSD_rad_1 ; 02_radius_mode1
printf, 1, PSD_sigma_1 ; 03_stddev_mode1
printf, 1, REFR_2 ; 04_n_real_mode2
printf, 1, REFI_2 ; 05_n_imag_mode2
printf, 1, PSD_rad_2 ; 06_radius_mode2
printf, 1, PSD_sigma_2 ; 07_stddev_mode2
printf, 1, FMF ; 08_fine-mode_fraction
printf, 1, AOD ; 09_AOD at reference wl
printf, 1, peak_height ; 10_Profile_peak
printf, 1, half_width ; 11_Profile_half_width
;printf, 1, relaxation_para ; 11_Profile_half_width
close, 1

;--write input surface option as .dat txt file
;--only if surface is lambertian assumption
;--can be changed as reading wl0,wl1, n_coeffi, 1, albedo(r) as in Albdeoclosure.cfg file
;--but it is one constant as here
openw, 1, save_pwd +'/Input_surface_lambertian.dat'
printf, 1, albedo_c0, albedo_c1 ;albedo(r) ; 00_n_real_mode1
close, 1

;printf, 1, 'FMF, PSD_rad_1, PSD_rad_2, PSD_sigma_1, PSD_sigma_2, REFR, REFI, '

;--run vlidort
spawn, './Gemstool_NSW_IQULinearized.exe'

spawn, 'cp fort.33 ./GEMSTOOL_NSW_Results/Optical_depth_all.dat'
spawn, 'cp fort.34 ./GEMSTOOL_NSW_Results/Optical_depth_gases.dat'

spawn, 'mv -f ./GEMSTOOL_NSW_Results/*.out* '+save_pwd
spawn, 'mv -f ./GEMSTOOL_NSW_Results/*.dat '+save_pwd


;---cp modified (if gaussian option) aod profile to save_pwd and state vectors pwd
file_copy, './aod_loading_user.input', save_pwd+'/aod_loading_user.input', /overwrite
file_copy, './aod_loading_user.input', state_vectors_pwd+'/state_vectors_aodprof_'+str_step+'.txt',  /overwrite

cd, workspace_pwd;'/home/mchoi/OSSE_tool'


endfor ; for i_wnloop


;---just integrate individual loop out files to one final integrated file















endfor ; for k = 0L, n_elements(RAA)-1 do begin
endfor ; for j = 0L, n_elements(VZA)-1 do begin
endfor ; for i = 0L, n_elements(SZA)-1 do begin
endfor ; for r = 0L, n_elements(albedo)-1 do begin
; endfor ; for m=10L, 25, 15 do begin ;10=M2 (Coarse and moderate absorbing), 25=N8 (small and non-absorbing)

endfor ;for i_single_AOD_layer = 0L, n_elements(idx_Single_AOD_layer)-1 do begin ; new 06/17/19

endfor; for i_Half   = 0L, n_elements(Half_width_node)-1 do begin
endfor ;for i_peak   = 0L, n_elements(Peak_height_node)-1 do begin
endfor ; for i_AOD=0L, n_elements(AOD_node)-1 do begin
endfor ; for i_level=0, 3 do begin
; endfor ; for i_Xsec=0, 1 do begin 

; endfor ; for i_prof=0, 1 do begin

; endfor ;for i_test = 0L, 0 do begin ; BRDF FD test

; endfor ;for i_band = 0L, 2 do begin


END
