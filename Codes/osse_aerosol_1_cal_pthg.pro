pro osse_aerosol_1_cal_pthg, workspace_pwd, i_step, i_loop_start, input_loop_begin, $
                                                lon_target, lat_target, alt_target, $
                                                Pres_CLARS, Temp_CLARS, $
                                                yy_read, mm_read, dd_read, hh_utc_read, $
                                                SZA_read, VZA_read, SAA_read, VAA_read, $
                                                i_case, savefolder_fin, savefolder_num, $
                                                switch_perturb, switch_constraint_aerosol, ii, wn_interval_read
                                                

result_pwd = workspace_pwd+'/Results/'+savefolder_fin
state_vectors_pwd = result_pwd+ '/state_vectors_step'
file_mkdir, state_vectors_pwd


;---read state vectors
str_step = 'step'+string(i_step,f='(i2.2)')

;-----read state vectors

if i_step ge 1 then begin
readcol, state_vectors_pwd+'/state_vectors_surfpress_'+str_step+'.txt', SV_surfpress, /silent
         surfpress=SV_surfpress(0)
endif

readcol, state_vectors_pwd+'/state_vectors_tshift_'+str_step+'.txt', SV_tshift, /silent
         tshift_read = SV_tshift(0)

;----other fixed conditions
Lat_dum=[  34.221d, lat_target] ;--[CLARS-FTS SVO mode, CLARS-FTS LABS mode]
Lon_dum=[-118.0571, lon_target]
Year_dum=yy_read
Month_dum=mm_read
DayofMonth_dum=dd_read
Hour_dum = hh_utc_read; local time + 8  = UTC
Minute_dum = 00
Second_dum = 00

SZA_meas = [SZA_read, SZA_read]
SAA_meas = [SAA_read, SAA_read]
VZA_meas = [68.633659, VZA_read]
VAA_meas = [8.8333292, VAA_read]
RAA_meas = abs(180.0 - abs(SAA_meas - VAA_meas))

Scat_meas = !radeg*acos(-cos(!dtor*SZA_meas)*cos(!dtor*VZA_meas)+sin(!dtor*SZA_meas)*sin(!dtor*VZA_meas)*cos(!dtor*RAA_meas))
sfc_alt_km = [1.673d, alt_target];[1.673d, 0.134d]

level_node_km = [1.673+0.0001, 1.673]

idx_target = 1 ;--0:CLARS   1:LABS

solz_fin = SZA_meas(idx_target)
senz_fin = VZA_meas(idx_target)
rela_fin = RAA_meas(idx_target)
scat_fin = Scat_meas(idx_target)

CLARS_alt_km_fin = sfc_alt_km(0)
CLARS_pres_hpa_fin = Pres_CLARS;--hPa; sfc_pres_hPa(0)
CLARS_temp_K_fin = Temp_CLARS;--K; sfc_temp_C(0)

LABS_alt_km_fin = sfc_alt_km(1)

;---LABS surface pressure/temperature will be defined after MERRA-2 calculation
if i_step ge 1 then LABS_pres_hpa_fin=surfpress;(as updated state vector);-- sfc_pres_hPa(1) (previous fixed)

lat_dum_fin = lat_Dum(idx_target);mean(lat_Dum) ;--mean of LABS and SVO
lon_dum_fin = lon_Dum(idx_target);mean(lon_Dum) ;--mean of LABS and SVO


level_node_km_fin = level_node_km(idx_target)

switch_cal_pthg = 1
if switch_cal_pthg eq 1  then begin
switch_synthetic = 1
sub_module_setup_pthg_clars_new_multiple, i_step, $
                    idx_target, lat_dum_fin,lon_dum_fin,Year_dum,Month_dum,DayofMonth_dum,Hour_dum,Minute_dum, $
                    CLARS_alt_km_fin, CLARS_pres_hPa_fin, CLARS_temp_K_fin, $
                    LABS_alt_km_fin, LABS_pres_hpa_fin, workspace_pwd, tshift_read, $
                    switch_synthetic, switch_perturb, switch_constraint_aerosol, $ ;--
                    i_case, savefolder_num, $ ;--for true value
                    result_pwd, ii,i_loop_start,input_loop_begin ;--deleted RTM_pwd and insert workspace_pwd


endif ;if switch_cal_pthg eq 1 then begin
; endif ;if i_loop_start eq 0 then begin


; if i_step eq 0 then begin
;---copy

for input_loop = 0L, 9 do begin
;NSW_heightgrid_1 in PTH_PROFILES folder
file_title = 'NSW_heightgrid_1'+'_ii'+string(ii,f='(i4.4)')+'.dat' 
original_file = workspace_pwd+'/RTM/GEMSTOOL_physicsdata/PTH_PROFILES'+'/'+file_title;'NSW_heightgrid_1.dat'
rtm_pwd = workspace_pwd+'/RTM/RTM_LAB/'+'GEMSTOOL_mchoi_LAB_loop'+string(input_loop+i_loop_start,f='(i3.3)')+'/GEMSTOOL_NSW_wrapper'
copy_target_file = rtm_pwd + '/NSW_heightgrid_1.dat'
FILE_COPY, original_file, copy_target_file, /overwrite

;atmos_mchoi in PTH_PROFILES folder and GAS_PROFILES folder
file_title = 'atmos_mchoi'+'_ii'+string(ii,f='(i4.4)')+'.dat' 
original_file = workspace_pwd+'/RTM/GEMSTOOL_physicsdata/PTH_PROFILES'+'/'+file_title;'NSW_heightgrid_1.dat'
rtm1_pwd = workspace_pwd+'/RTM/RTM_LAB/'+'GEMSTOOL_mchoi_LAB_loop'+string(input_loop+i_loop_start,f='(i3.3)')+'/GEMSTOOL_physicsdata/PTH_PROFILES'
rtm2_pwd = workspace_pwd+'/RTM/RTM_LAB/'+'GEMSTOOL_mchoi_LAB_loop'+string(input_loop+i_loop_start,f='(i3.3)')+'/GEMSTOOL_physicsdata/GAS_PROFILES'

copy_target_file = rtm1_pwd + '/atmos_mchoi.dat'
FILE_COPY, original_file, copy_target_file, /overwrite

copy_target_file = rtm2_pwd + '/atmos_mchoi.dat'
FILE_COPY, original_file, copy_target_file, /overwrite

endfor



;--calculate Xsec LUT from each i_loop_start folder and be copied to 1-9 next folders.
sub_make_xsec_lut_using_absco_v5_mchoi_o2ab_new_multiple, workspace_pwd, i_loop_start,input_loop_begin, wn_interval_read




;---copy to state vectors folder 
file_title = 'NSW_heightgrid_1'+'_ii'+string(ii,f='(i4.4)')+'.dat' 
original_file = workspace_pwd+'/RTM/GEMSTOOL_physicsdata/PTH_PROFILES'+'/'+file_title;'NSW_heightgrid_1.dat'
copy_target_file = state_vectors_pwd+'/NSW_heightgrid_1.txt'
FILE_COPY, original_file, copy_target_file, /overwrite
; endif ;if i_step eq 0 then begin


;---copy of PTHG file according to state vector folder
file_title = 'atmos_mchoi'+'_ii'+string(ii,f='(i4.4)')+'.dat'
original_file = workspace_pwd+'/RTM/GEMSTOOL_physicsdata/PTH_PROFILES/'+file_title;'atmos_mchoi.dat'
copy_target_file = state_vectors_pwd+'/atmos_mchoi_'+str_step+'.txt'
FILE_COPY, original_file, copy_target_file, /overwrite

;---print
; print, 'PTHG  and Xsection calculation is finished!!'

;---just delete below processes

END
