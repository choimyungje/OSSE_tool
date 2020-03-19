Pro runtest8_4_save_as_synthetic_data, i_step, ii

switch_O2D = 1

switch_Pol = 1

workspace_pwd = '/home/mchoi/HIMAP/workspace_t16_12_LABS_lamb'

; yy_read = 2019 & str_yy_read = string(yy_read,f='(i4.4)')
; mm_read = 07   & str_mm_read = string(mm_read,f='(i2.2)')
; dd_read = 05   & str_dd_read = string(dd_read,f='(i2.2)')
; str_yymmdd_folder = str_yy_read+'-'+str_mm_read+'-'+str_dd_read
; savefolder_fin = '201907052231'
; str_savefolder_fin = strtrim(string(savefolder_fin),2)


; AOD_set = [0.6];[0.2, 0.6, 1.0]
; Peak_set = [1.0];[0.2, 0.6, 1.0, 1.4, 1.8]
; width_set = [0.3];[0.1, 0.3, 0.5]


; for i_AOD = 0L, n_elements(AOD_set)-1 do begin
; for i_peak = 0L, n_elements(Peak_set)-1 do begin
; for i_width = 0L, n_elements(width_set)-1 do begin

; str_AOD = string(AOD_set(i_AOD),f='(f3.1)')
; str_Peak = string(Peak_set(i_peak),f='(f3.1)')
; str_Width = string(width_set(i_width),f='(f3.1)')

; str_sub = 'AOD'+str_AOD+'_Peak'+str_Peak+'_Width'+str_Width


str_ii = 'ii'+string(ii,f='(i4.4)')


; data_pwd1 = file_search(workspace_pwd+'/Results/Synthetic_cl'+str_savefolder_fin+'*'+str_sub)
data_pwd1 = file_search(workspace_pwd+'/Results/Results_justcal8_*'+str_ii)
data_pwd1 = data_pwd1(0)
SV_pwd1 = data_pwd1+'/state_vectors_step'
save_pwd0 = data_pwd1+'/Analysis'
file_mkdir, save_pwd0

; ;---measurement data
; meas_file = file_search('/home/mchoi/HIMAP/data_CLARS_FTS/CoaddedSPC_XB/'+str_yymmdd_folder+'/','*'+str_savefolder_fin+'*.*', count=n_meas_pwd)
; readcol, meas_file, wn_meas, S000_read
; wn_meas = wn_meas(1:n_elements(wn_meas)-1)
; S000_read = S000_read(1:n_elements(S000_read)-1)
; ; n_meas_read = n_elements(wn_meas_read)

; ;---calibration data
; mcal_file = '/home/mchoi/HIMAP/data_CLARS_FTS/calibration_20190715/201907152039_sun_CLARS_2Chan_V2.cal'
; readcol, mcal_file, wn_mcal, S_mcal

; S000_read = S000_read * S_mcal ;--then unit should be W/m^2/cm^-1/sr


;---grid plot
readcol, SV_pwd1+'/NSW_heightgrid_1.txt', level, altitude, /silent

str_step0 = 'step'+string(i_step-1,f='(i2.2)')
str_step = 'step'+string(i_step,f='(i2.2)')
str_step2 = 'step'+string(i_step+1,f='(i2.2)')

i_prof = 0 & str_prof = 'AerosolBulk';---then bulk option and gaussian
; i_prof = 1 & str_prof = 'AerosolProf';---then AOD profile option 
; i_brdf = 0 ;---then albedo jacobian
i_brdf = 1 ; then brdf jacobian

for i_band = 0L, switch_O2D do begin

if i_Band eq 0 then begin 
; wn_range = [12900.,13250.]
wn_range = [12930.,13250.]
wn_range_ret = [12960.-10, 13220.+10]
endif
if i_Band eq 1 then begin 
; wn_range = [07680.,08020.]   
wn_range = [07720.,08020.]   
wn_range_ret = [7750.-10, 7990.+10]
endif

str_subfolder_wnrange = string(wn_range(0), f='(f08.2)')+'_'+string(wn_range(1), f='(f08.2)')
save_pwd = save_pwd0
; save_pwd = save_pwd0 + '/'+str_subfolder_wnrange
file_mkdir, save_pwd


; data_pwd_fin = $
;       [data_pwd0 + '/SVO_lamb_ALB0.999000_AOD0.0000_ABSCOv5'+'/'+str_subfolder_wnrange, $ ; SVO
;        data_pwd1 + '/LABS_BRDF_BPDF_'+str_prof+'_'+str_step+'_ABSCOv5'+'/'+str_subfolder_wnrange]
; data_pwd_fin = $
;       [data_pwd1 + '/SVO_'+str_prof+'_'+str_step+'/'+str_subfolder_wnrange, $ ; SVO
;        data_pwd1 + '/LAB_'+str_prof+'_'+str_step+'/'+str_subfolder_wnrange]
data_pwd_fin = $
      [data_pwd1 + '/SVO_'+str_step+'/'+str_subfolder_wnrange, $ ; SVO
       data_pwd1 + '/LAB_'+str_step+'/'+str_subfolder_wnrange] ;LAB

n_sin = n_elements(data_pwd_fin)

wn_range_ini = wn_range_ret
sampling_res = 0.10d; 0.02d ;wn_cal(1)-wn_cal(0)
wn_range_fin=[wn_range_ini(0)+10.0, wn_range_ini(1)-10.0] ;should check
n_wn_inner = round(((wn_range_fin(1)-wn_range_fin(0))/sampling_res) + 1)
wn_inner = wn_range_fin(0)+sampling_Res*findgen(n_wn_inner) ; identical with original data 
  
; pwd_meas = '/home/mchoi/HIMAP/data_CLARS_FTS/20190529/cl20190529pol_3cm'
; refm_file = pwd_meas + '/cl201905292119000_.001' ;Spectralon target; revised as 3cm; just for wn points settting.

;-for operational products
; pwd_meas = '/home/mchoi/HIMAP/data_CLARS_FTS/CoaddedSPC_XB'
; refm_file = pwd_meas + '/2019-07-15/cl201907151503Ange.001' ;Spectralon target; revised as 3cm; just for wn points settting.

; readcol, refm_file, wn_refm, S_refm, format='D,D', /silent

; ;--2019-11-05; original
; idx_select = where(wn_meas ge wn_range_fin(0) and wn_meas le wn_range_fin(1), n_idx_select)
; wn_conv = wn_meas(idx_select)
; n_conv = n_elements(wn_conv) 
;--2019-11-05; modified
; idx_select = where(wn_meas ge wn_range_fin(0) and wn_meas le wn_range_fin(1), n_idx_select)
meas_res = 0.50d/2.0
n_conv = round(((wn_range_fin(1)-wn_range_fin(0))/meas_res) + 1)
wn_conv = wn_range_fin(0)+meas_res*findgen(n_conv)
; n_conv = n_elements(wn_conv) 

; ;---extended wavenumber for SNR calculation ;---just use logged SNR!
; idx_select_ext = where(wn_meas ge 12500.00 and wn_meas le 14000.00, n_idx_select)
; wn_conv_ext = wn_meas(idx_select_ext)
; n_conv_ext = n_elements(wn_conv_ext)

n_stokes = 4 ;I, Q, U, D


;----0)SVO
array_fin0 = dblarr(n_conv)*!values.d_nan

Solar_cal_conv = array_Fin0
I_cal_SVO = array_fin0
Q_cal_SVO = array_fin0
U_cal_SVO = array_fin0
D_cal_SVO = array_fin0

I2_cal_SVO = array_fin0
Q2_cal_SVO = array_fin0
U2_cal_SVO = array_fin0
D2_cal_SVO = array_fin0

;----1)LAB - SCIAMACHY_lambertian
I_cal_LAB = array_fin0
Q_cal_LAB = array_fin0
U_cal_LAB = array_fin0
D_cal_LAB = array_fin0
DoLP_cal_LAB = array_fin0

I2_cal_LAB = array_fin0
Q2_cal_LAB = array_fin0
U2_cal_LAB = array_fin0
D2_cal_LAB = array_fin0

;----read state vectors

readcol, SV_pwd1 + '/state_vectors_aerosol_'+str_step+'.txt' , prop_aerosol_LAB, format='d' ,/silent
; prop_aerosol_LAB = prop_aerosol_LAB[9:11];---only three (AOD, peak height, half width) for Jacobian.. 
; n_prop_aerosol_LAB = n_elements(prop_aerosol_LAB) 
n_prop_aerosol_LAB = n_elements(prop_aerosol_LAB) ;--for all data

; readcol, SV_pwd1 + '/state_vectors_albedo_'+str_step+'.txt', prop_sfcref_LAB, format='d'  ,/silent
;--not used but just put dump or temporary values
prop_sfcref_LAB = [0.985, 0.000001]
n_prop_sfcref_LAB = n_elements(prop_sfcref_LAB)

readcol, SV_pwd1 + '/state_vectors_aodprof_'+str_step+'.txt' , dum, prop_AODprof_LAB, format='i, d' ,/silent
n_prop_aodprof_LAB = n_elements(prop_aodprof_LAB)

idx_aodprof2 = where(prop_AODprof_LAB ne 0.0, n_idx_aodprof2);[11:25] ;--specific layers
prop_AODprof2_LAB = prop_AODprof_LAB(idx_aodprof2)
n_prop_aodprof2_LAB = n_elements(prop_AODprof2_LAB)

readcol, SV_pwd1 + '/state_vectors_brdf_'+str_step+'_O2A.txt', prop_brdf_LAB, format='d', /silent
prop_brdf_LAB = prop_brdf_LAB[0:2] ;--only for kernel scailing factors
n_prop_brdf_LAB = n_elements(prop_brdf_LAB)
;--O2A O2D seperated
readcol, SV_pwd1 + '/state_vectors_brdf_'+str_step+'_O2A.txt', prop_brdf_LAB, format='d', /silent
prop_brdf_LAB_O2A = prop_brdf_LAB[0:2] ;--only for kernel scailing factors
readcol, SV_pwd1 + '/state_vectors_brdf_'+str_step+'_O2D.txt', prop_brdf_LAB, format='d', /silent
prop_brdf_LAB_O2D = prop_brdf_LAB[0:2] ;--only for kernel scailing factors


readcol, SV_pwd1 + '/state_vectors_surfpress_'+str_step+'.txt' , prop_sfcpres_LAB, format='d' ,/silent
n_prop_sfcpres_LAB = n_elements(prop_sfcpres_LAB)

readcol, SV_pwd1 + '/state_vectors_h2oscaling_'+str_step+'.txt' , prop_h2oscaling_LAB, format='d' ,/silent
n_prop_h2oscaling_LAB = n_elements(prop_h2oscaling_LAB)

readcol, SV_pwd1 + '/state_vectors_tshift_'+str_step+'.txt' , prop_tshift_LAB, format='d' ,/silent
n_prop_tshift_LAB = n_elements(prop_tshift_LAB)



;---layer km setting
n_level=n_elements(level)
layer_thick_km = altitude(0:n_level-2)-altitude(1:n_level-1)
n_layer = n_elements(layer_thick_km)
altitude_layer = (altitude(0:n_level-2)+altitude(1:n_level-1))/2.0
; print, layer_thick_km[idx_aodprof2]


; JS_aerosol_LAB = dblarr(n_conv, n_prop_aerosol_LAB, 4) * !values.d_nan
; JS_sfcref_LAB = dblarr(n_conv, n_prop_sfcref_LAB, 4) * !values.d_nan
; JS_sfcpres_LAB = dblarr(n_conv, n_prop_sfcpres_LAB, 4) * !values.d_nan
; JS_h2oscaling_LAB = dblarr(n_conv, n_prop_h2oscaling_LAB, 4) * !values.d_nan
; JS_tshift_LAB = dblarr(n_conv, n_prop_tshift_LAB, 4) * !values.d_nan
; JS_aodprof2_LAB = dblarr(n_conv, n_prop_aodprof2_LAB, 4) * !values.d_nan
; JS_brdf_LAB = dblarr(n_conv, n_prop_brdf_LAB, 4) * !values.d_nan


for i_sin = 1, 1 do begin;i_sin = 0L, n_sin-1 do begin ;--i_sin 0: SVO ; i_sin 1: LAB

;--1) read IQU
readcol, data_pwd_fin(i_sin) + '/Stokes_IQU_DoLP_TOAUP.out', dum, wn_cal, dum, I_cal, Q_cal, U_cal, D_cal, format='D,D,D,D,D,D,D', /silent

D_cal = sqrt(Q_cal^2.0d + U_cal^2.0d)
DoLP_cal = sqrt((Q_cal^2.0d + U_cal^2.0d) / I_cal^2.0d)


; idx_sampling=where((round(wn_cal*100) mod 2) eq 0,n_idx_sampling)
; wn_cal=wn_Cal(idx_sampling)
; I_cal=I_cal(idx_sampling)
; Q_cal=Q_cal(idx_sampling)
; U_cal=U_cal(idx_sampling)
; D_cal=D_cal(idx_sampling)


;---2. read SOLSPEC+GeffToon solar irradiance for o2b/a + 1.27um range also
temp_pwd='/home/mchoi/HIMAP/data_solar_spectrum/SOLSPEC_Toon_convolved'
restore, temp_pwd+'/Solar_SOLSPEC_TOON_for_O2AB_07000_15000_mchoi.xdr'
;----save, wn_toon_select, SOLSPEC_Toon_select_cm, /variables, filename=temp_pwd+'/Solar_SOLSPEC_TOON_'+string_range+'_mchoi.xdr'
wn_sao=wn_toon_select
solar_sao=SOLSPEC_Toon_select_cm

;----3. doppler effect application to solar ;---ignored because very little change
; temp_pwd='/home/mchoi/HIMAP/for_doppler/doppler_mchoi/bin'
; readcol, temp_pwd+'/osds.out', osds, format='D',/silent
wn_sao2=wn_sao   ;wn_sao*(1.-osds(0)*10.^(-6.))

;----solar spectrum convolution to the simulation
;----should be changed to the convolution shape, but just interpolation in this calculation
solar_sao1=interpol(solar_sao,wn_sao,wn_sao2)

;----interpol to the calculation wn grid
solar_cal=interpol(solar_sao1, wn_sao2, wn_cal)
; I2_cal=I_cal * solar_cal ; units=W/m2/cm-1 ; ODD number!
; Q2_cal=Q_cal * solar_cal ; units=W/m2/cm-1 ; ODD number!
; U2_cal=U_cal * solar_cal ; units=W/m2/cm-1 ; ODD number!

; ;--for wavenumber shift calculation (later)
; if i_sin eq 0  then begin
; wn_refc = wn_cal
; I_refc  = I_cal * solar_cal
; Q_refc  = Q_cal * solar_cal
; U_refc  = U_cal * solar_cal
; endif

; ;--2) read Jacobian (calculated) for i_sin=1 LAB
; if i_sin eq 1 then begin ;---LAB

; JS_Aerosol_LAB_read = dblarr(n_elements(I_cal), n_prop_aerosol_LAB, 4) * !values.d_nan
; JS_sfcref_LAB_read = dblarr(n_elements(I_cal), n_prop_sfcref_LAB, 4) * !values.d_nan
; JS_sfcpres_LAB_read = dblarr(n_elements(I_cal), n_prop_sfcpres_LAB, 4) * !values.d_nan
; JS_h2oscaling_LAB_read = dblarr(n_elements(I_cal), n_prop_h2oscaling_LAB, 4) * !values.d_nan
; JS_tshift_LAB_read = dblarr(n_elements(I_cal), n_prop_tshift_LAB, 4)*!values.d_nan
; JS_aodprof2_LAB_read = dblarr(n_elements(I_cal), n_prop_aodprof2_LAB, 4) * !values.d_nan
; JS_brdf_LAB_read = dblarr(n_elements(I_Cal), n_prop_brdf_LAB, 4) * !values.d_nan


; for i_stokes=0, n_stokes-2 do begin

;   if i_brdf eq 0 then begin
;     if i_stokes eq 0 then file_read = data_pwd_fin(i_sin) + '/Stokes_I_AlbWFs_TOAUP_slope.out'
;     if i_stokes eq 1 then file_read = data_pwd_fin(i_sin) + '/Stokes_Q_AlbWFs_TOAUP_slope.out'
;     if i_stokes eq 2 then file_read = data_pwd_fin(i_sin) + '/Stokes_U_AlbWFs_TOAUP_slope.out'
;     openr,1,file_read
;     n_lines = file_lines(file_read)
;     data=strarr(n_lines)
;     readf,1,data
;     close,1
;     for l=0L, n_lines-1 do begin
;     dum = strsplit(data(l), ' ',/extract)
;     JS_sfcref_LAB_read(l,0:n_prop_sfcref_LAB-1,i_stokes) = dum(2:2+n_prop_sfcref_LAB-1)
;     endfor
;     if i_stokes eq 2 then begin ;---calculate DoLP jacobian
;       Q2U2 = (Q_cal * Q_cal + U_cal * U_cal)
;       for i=0L, n_prop_sfcref_LAB-1 do begin
;       JS_I = JS_sfcref_LAB_read(*,i,0)
;       JS_Q = JS_sfcref_LAB_read(*,i,1)
;       JS_U = JS_sfcref_LAB_read(*,i,2)
;       JS_sfcref_LAB_read(*,i,3) = $
;                     1./(2.*I_cal)*(Q2U2^(-0.5))*(2.*Q_cal*JS_Q+2.*U_cal*JS_U) + $
;                     (-1.0)*(Q2U2^(0.5))/(I_cal^2.0)*JS_I
;       endfor
;     endif
;   endif ;if i_brdf eq 0 then begin

;   if i_brdf eq 1 then begin ;--brdf jacobian
  
;     if i_stokes eq 0 then file_read = data_pwd_fin(i_sin) + '/Stokes_I_BRDFWFs_TOAUP.out'
;     if i_stokes eq 1 then file_read = data_pwd_fin(i_sin) + '/Stokes_Q_BRDFWFs_TOAUP.out'
;     if i_stokes eq 2 then file_read = data_pwd_fin(i_sin) + '/Stokes_U_BRDFWFs_TOAUP.out'
;     openr, 1, file_read
;     n_lines = file_lines(file_read)
;     data = strarr(n_lines)
;     readf, 1, data
;     close, 1
;     for l=0L, n_lines-1 do begin
;     dum = strsplit(data(l), ' ',/extract)
;     JS_brdf_LAB_read(l,0:n_prop_brdf_LAB-1,i_stokes) = dum(2:2+n_prop_brdf_LAB-1)
;     endfor
;     if i_stokes eq 2 then begin ;---calculate DoLP jacobian
;       Q2U2 = (Q_cal * Q_cal + U_cal * U_cal)
;       for i=0L, n_prop_brdf_LAB-1 do begin
;       JS_I = JS_brdf_LAB_read(*,i,0)
;       JS_Q = JS_brdf_LAB_read(*,i,1)
;       JS_U = JS_brdf_LAB_read(*,i,2)
;       JS_brdf_LAB_read(*,i,3) = $
;                     1./(2.*I_cal)*(Q2U2^(-0.5))*(2.*Q_cal*JS_Q+2.*U_cal*JS_U) + $
;                     (-1.0)*(Q2U2^(0.5))/(I_cal^2.0)*JS_I
;       endfor
;     endif
    
;   endif


;   if i_prof eq 1 then begin
;     ;---AODprof
;     if i_stokes eq 0 then file_read = data_pwd_fin(i_sin) + '/Stokes_I_AerProfWFs_TOAUP.out'
;     if i_stokes eq 1 then file_read = data_pwd_fin(i_sin) + '/Stokes_Q_AerProfWFs_TOAUP.out'
;     if i_stokes eq 2 then file_read = data_pwd_fin(i_sin) + '/Stokes_U_AerProfWFs_TOAUP.out'
;     openr,1,file_read
;     n_lines = file_lines(file_read)
;     data=strarr(n_lines)
;     readf,1,data
;     close,1
;     for l=0L, n_lines-1 do begin
;     dum = strsplit(data(l), ' ',/extract)
;     JS_aodprof2_LAB_read(l,0:n_prop_aodprof2_LAB-1,i_stokes) = dum(2+idx_aodprof2(0):2+idx_aodprof2(0)+n_prop_aodprof2_LAB-1)
;     endfor
;     if i_stokes eq 2 then begin ;---calculate DoLP jacobian
;         Q2U2 = (Q_cal * Q_cal + U_cal * U_cal)
;         for i=0L, n_prop_aodprof2_LAB-1 do begin
;         JS_I = JS_aodprof2_LAB_read(*,i,0)
;         JS_Q = JS_aodprof2_LAB_read(*,i,1)
;         JS_U = JS_aodprof2_LAB_read(*,i,2)
;         JS_aodprof2_LAB_read(*,i,3) = $
;                       1./(2.*I_cal)*(Q2U2^(-0.5))*(2.*Q_cal*JS_Q+2.*U_cal*JS_U) + $
;                       (-1.0)*(Q2U2^(0.5))/(I_cal^2.0)*JS_I
;         endfor
;     endif
;   endif

;   if i_prof eq 0 then begin ;---bulk
;     ;----Aerosol Bulk
;     if i_stokes eq 0 then file_read = data_pwd_fin(i_sin) + '/Stokes_I_AerBulkWFs_TOAUP.out'
;     if i_stokes eq 1 then file_read = data_pwd_fin(i_sin) + '/Stokes_Q_AerBulkWFs_TOAUP.out'
;     if i_stokes eq 2 then file_read = data_pwd_fin(i_sin) + '/Stokes_U_AerBulkWFs_TOAUP.out'
;     openr,1,file_read
;     n_lines = file_lines(file_read)
;     data=strarr(n_lines)
;     readf,1,data
;     close,1
;     for l=0L, n_lines-1 do begin
;     dum = strsplit(data(l), ' ',/extract)
;     JS_Aerosol_LAB_read(l,0:n_prop_aerosol_LAB-1,i_stokes) = dum(2:2+n_prop_aerosol_LAB-1)
;     endfor
;     if i_stokes eq 2 then begin ;---calculate DoLP jacobian
;         Q2U2 = (Q_cal * Q_cal + U_cal * U_cal)
;         for i=0L, n_prop_aerosol_LAB-1 do begin
;         JS_I = JS_Aerosol_LAB_read(*,i,0)
;         JS_Q = JS_Aerosol_LAB_read(*,i,1)
;         JS_U = JS_Aerosol_LAB_read(*,i,2)
;         JS_Aerosol_LAB_read(*,i,3) = $
;                       1./(2.*I_cal)*(Q2U2^(-0.5))*(2.*Q_cal*JS_Q+2.*U_cal*JS_U) + $
;                       (-1.0)*(Q2U2^(0.5))/(I_cal^2.0)*JS_I
;         endfor
;     endif ;;----Aerosol Bulk
;     ;----Surface pressure
;     if i_stokes eq 0 then file_read = data_pwd_fin(i_sin) + '/Stokes_I_SURFPWFs_TOAUP.out'
;     if i_stokes eq 1 then file_read = data_pwd_fin(i_sin) + '/Stokes_Q_SURFPWFs_TOAUP.out'
;     if i_stokes eq 2 then file_read = data_pwd_fin(i_sin) + '/Stokes_U_SURFPWFs_TOAUP.out'
;     openr,1,file_read
;     n_lines = file_lines(file_read)
;     data=strarr(n_lines)
;     readf,1,data
;     close,1
;     for l=0L, n_lines-1 do begin
;     dum = strsplit(data(l), ' ',/extract)
;     JS_sfcpres_LAB_read(l,0:n_prop_sfcpres_LAB-1,i_stokes) = dum(2:2+n_prop_sfcpres_LAB-1)
;     endfor
;     if i_stokes eq 2 then begin ;---calculate DoLP jacobian
;         Q2U2 = (Q_cal * Q_cal + U_cal * U_cal)
;         for i=0L, n_prop_sfcpres_LAB-1 do begin
;         JS_I = JS_sfcpres_LAB_read(*,i,0)
;         JS_Q = JS_sfcpres_LAB_read(*,i,1)
;         JS_U = JS_sfcpres_LAB_read(*,i,2)
;         JS_sfcpres_LAB_read(*,i,3) = $
;                       1./(2.*I_cal)*(Q2U2^(-0.5))*(2.*Q_cal*JS_Q+2.*U_cal*JS_U) + $
;                       (-1.0)*(Q2U2^(0.5))/(I_cal^2.0)*JS_I
;         endfor
;     endif ;----Surface pressure
;     ;----H2O scaling
;     if i_stokes eq 0 then file_read = file_search(data_pwd_fin(i_sin) + '/Stokes_I_WSCALEWFs_TOAUP.out', count=n_file_read)
;     if i_stokes eq 1 then file_read = file_search(data_pwd_fin(i_sin) + '/Stokes_Q_WSCALEWFs_TOAUP.out', count=n_file_read)
;     if i_stokes eq 2 then file_read = file_search(data_pwd_fin(i_sin) + '/Stokes_U_WSCALEWFs_TOAUP.out', count=n_file_read)
;     ; JS_h2oscaling_LAB_read(*)=0.0d
;     if n_file_read eq 1 then begin
;     openr,1,file_read(0)
;     n_lines = file_lines(file_read(0))
;     data=strarr(n_lines)
;     readf,1,data
;     close,1
;     for l=0L, n_lines-1 do begin
;     dum = strsplit(data(l), ' ',/extract)
;     JS_h2oscaling_LAB_read(l,0:n_prop_h2oscaling_LAB-1,i_stokes) = dum(2:2+n_prop_h2oscaling_LAB-1)
;     endfor
;     if i_stokes eq 2 then begin ;---calculate DoLP jacobian
;         Q2U2 = (Q_cal * Q_cal + U_cal * U_cal)
;         for i=0L, n_prop_H2Oscaling_LAB-1 do begin
;         JS_I = JS_h2oscaling_LAB_read(*,i,0)
;         JS_Q = JS_h2oscaling_LAB_read(*,i,1)
;         JS_U = JS_h2oscaling_LAB_read(*,i,2)
;         JS_h2oscaling_LAB_read(*,i,3) = $
;                       1./(2.*I_cal)*(Q2U2^(-0.5))*(2.*Q_cal*JS_Q+2.*U_cal*JS_U) + $
;                       (-1.0)*(Q2U2^(0.5))/(I_cal^2.0)*JS_I
;         endfor
;     endif ;if i_stokes eq 2 then begin ;---calculate DoLP jacobian
;     endif ;if n_file_read eq 1 then begin
;     ;----H2O scaling
;     ;----T shift
;     if i_stokes eq 0 then file_read = file_search(data_pwd_fin(i_sin) + '/Stokes_I_TSHIFTWFs_TOAUP.out', count=n_file_read)
;     if i_stokes eq 1 then file_read = file_search(data_pwd_fin(i_sin) + '/Stokes_Q_TSHIFTWFs_TOAUP.out', count=n_file_read)
;     if i_stokes eq 2 then file_read = file_search(data_pwd_fin(i_sin) + '/Stokes_U_TSHIFTWFs_TOAUP.out', count=n_file_read)
;     ; JS_h2oscaling_LAB_read(*)=0.0d
;     if n_file_read eq 1 then begin
;     openr,1,file_read(0)
;     n_lines = file_lines(file_read(0))
;     data=strarr(n_lines)
;     readf,1,data
;     close,1
;     for l=0L, n_lines-1 do begin
;     dum = strsplit(data(l), ' ',/extract)
;     JS_tshift_LAB_read(l,0:n_prop_tshift_LAB-1,i_stokes) = dum(2:2+n_prop_tshift_LAB-1)
;     endfor
;     if i_stokes eq 2 then begin ;---calculate DoLP jacobian
;         Q2U2 = (Q_cal * Q_cal + U_cal * U_cal)
;         for i=0L, n_prop_H2Oscaling_LAB-1 do begin
;         JS_I = JS_tshift_LAB_read(*,i,0)
;         JS_Q = JS_tshift_LAB_read(*,i,1)
;         JS_U = JS_tshift_LAB_read(*,i,2)
;         JS_tshift_LAB_read(*,i,3) = $
;                       1./(2.*I_cal)*(Q2U2^(-0.5))*(2.*Q_cal*JS_Q+2.*U_cal*JS_U) + $
;                       (-1.0)*(Q2U2^(0.5))/(I_cal^2.0)*JS_I
;         endfor
;     endif ;if i_stokes eq 2 then begin ;---calculate DoLP jacobian
;     endif ;if n_file_read eq 1 then begin
;     ;----T shift
;   endif ;if i_prof eq 0 then begin ;---bulk
; endfor ;for i_stokes=0, n_stokes-2 do begin
; endif ;if i_sin eq 1 then begin

; ;---filling nan_idx as 'zero' for JS_sfcpres_LAB_read, JS_h2oscaling_LAB_read, JS_tshift_LAB_read
; nan_idx = where(finite(JS_sfcpres_LAB_read) eq 0, n_nan_idx)
; if n_nan_idx ge 1 then JS_sfcpres_LAB_read(nan_idx) = 0.0d
; nan_idx = where(finite(JS_h2oscaling_LAB_read) eq 0, n_nan_idx)
; if n_nan_idx ge 1 then JS_h2oscaling_LAB_read(nan_idx) = 0.0d
; nan_idx = where(finite(JS_tshift_LAB_read) eq 0, n_nan_idx)
; if n_nan_idx ge 1 then JS_tshift_LAB_read(nan_idx) = 0.0d




; ;----Jacobian smoothing test
; if i_sin eq 1 then begin 


; testvar0 = JS_aodprof2_LAB_read(*,0,0)
; testvar1 = JS_aodprof2_LAB_read(*,n_prop_aodprof2_LAB-1,0)
; n_wn = n_elements(I_cal)
; mean_of_sides = [testvar0(0), 0.5*(testvar0(0:n_wn-3)+testvar0(2:n_wn-1)), testvar0(n_wn-1)]

; median_temp0 = testvar0 * !Values.f_nan
; median_temp1 = testvar0 * !Values.f_nan
; median_temp2 = testvar0 * !Values.f_nan

; n_points = 7 ;!odd number

; for i_median=0+(n_points-1)/2, n_wn-1-(n_points-1)/2 do $
;   median_temp0(i_median) = median(testvar0[i_median-2:i_median+2],/even)

; for i_median=0+(n_points-1)/2, n_wn-1-(n_points-1)/2 do $
;   median_temp1(i_median) = median(median_temp0[i_median-2:i_median+2],/even)

; for i_median=0+(n_points-1)/2, n_wn-1-(n_points-1)/2 do $
;   median_temp2(i_median) = median(median_temp1[i_median-2:i_median+2],/even)


; temp = mean_of_sides/testvar0


; cgplot, wn_Cal, temp, yrange=[-10,1000]

; cgplot, wn_cal, testvar0, yrange=[-0.1,0.1]
; cgplot, wn_cal, smooth(testvar0,5), psym=1, xrange=[13050,13070], yrange=[-0.1,0.1]

; xrange_dum = [13050,13070]
; cgplot, wn_cal, testvar0, psym=1, xrange=xrange_dum, yrange=[-0.1,0.1]
; cgplot, wn_cal, median_temp0, psym=1, xrange=xrange_dum, /overplot, color='red'
; cgplot, wn_cal, median_temp2, psym=1, xrange=xrange_dum, /overplot, color='blue'

; var3 = median_temp0(sort(median_temp0))

; cgplot, var3, psym=1, xrange=[17000,18000]

; ;---smoothing method thinking....
; ;-option 1) original jacobian order itself -> cons) hard to fit in local several points outliers
; ;-option 2) sorting according to Jacobian 
; ;-option 3) sorting accoridng to Intensity or DoLP <- this might be better....


; stop
; endif ;if i_sin eq 1 then

usable_idx = where(wn_cal ge wn_range(0) and wn_cal le wn_range(1), n_usable_idx)

if i_sin eq 0 then S_ref_all = [[I_cal], [Q_cal], [U_cal], [solar_cal]]
if i_sin eq 1 then begin
; S_ref_all_I = [[I_cal], [JS_aerosol_LAB_read(*,*,0)], [JS_sfcpres_LAB_read(*,*,0)], [JS_sfcref_LAB_read(*,*,0)]]
; S_ref_all_Q = [[Q_cal], [JS_aerosol_LAB_read(*,*,1)], [JS_sfcpres_LAB_read(*,*,1)], [JS_sfcref_LAB_read(*,*,1)]]
; S_ref_all_U = [[U_cal], [JS_aerosol_LAB_read(*,*,2)], [JS_sfcpres_LAB_read(*,*,2)], [JS_sfcref_LAB_read(*,*,2)]]
; S_ref_all = [[S_ref_all_I], [S_ref_all_Q], [S_ref_all_U]]
; node = [0, 15, 30]; 0_0=I_cal; 1_12=aerosol; 13_13=sfcpres; 14_15=sfcref ; and same for Q and U


if i_prof eq 0 and i_brdf eq 0 then begin
;---State vector: sfcref (2) and AOD Bulk (12)
S_ref_all_I = [[I_cal], [JS_sfcref_LAB_read(*,*,0)],[JS_Aerosol_LAB_read(*,*,0)]]
S_ref_all_Q = [[Q_cal], [JS_sfcref_LAB_read(*,*,1)],[JS_Aerosol_LAB_read(*,*,1)]]
S_ref_all_U = [[U_cal], [JS_sfcref_LAB_read(*,*,2)],[JS_Aerosol_LAB_read(*,*,2)]]
S_ref_all = [[S_ref_all_I], [S_ref_all_Q], [S_ref_all_U]]
n_all_each = (size(s_ref_all_I,/dimension))[1]
node = [0, n_all_each, n_all_each*2] 

endif 

if i_prof eq 1 and i_brdf eq 0 then begin
;---State vector: sfcref (2) and AODprof (16)
S_ref_all_I = [[I_cal], [JS_sfcref_LAB_read(*,*,0)],[JS_aodprof2_LAB_read(*,*,0)]]
S_ref_all_Q = [[Q_cal], [JS_sfcref_LAB_read(*,*,1)],[JS_aodprof2_LAB_read(*,*,1)]]
S_ref_all_U = [[U_cal], [JS_sfcref_LAB_read(*,*,2)],[JS_aodprof2_LAB_read(*,*,2)]]
S_ref_all = [[S_ref_all_I], [S_ref_all_Q], [S_ref_all_U]]
n_all_each = (size(s_ref_all_I,/dimension))[1]
node = [0, n_all_each, n_all_each*2] 
endif 

if i_prof eq 0 and i_brdf eq 1 then begin
;---State vector: BRDF(3 + 0+0+2) + AOD Bulk (12) + sfcpres (1) + H2O scaling (1)
; S_ref_all_I = [[I_cal], [JS_brdf_LAB_read(*,*,0)],[JS_Aerosol_LAB_read(*,*,0)],[JS_sfcpres_LAB_read(*,*,0)],[JS_h2oscaling_LAB_read(*,*,0)],[JS_tshift_LAB_read(*,*,0)]]
; S_ref_all_Q = [[Q_cal], [JS_brdf_LAB_read(*,*,1)],[JS_Aerosol_LAB_read(*,*,1)],[JS_sfcpres_LAB_read(*,*,1)],[JS_h2oscaling_LAB_read(*,*,1)],[JS_tshift_LAB_read(*,*,1)]]
; S_ref_all_U = [[U_cal], [JS_brdf_LAB_read(*,*,2)],[JS_Aerosol_LAB_read(*,*,2)],[JS_sfcpres_LAB_read(*,*,2)],[JS_h2oscaling_LAB_read(*,*,2)],[JS_tshift_LAB_read(*,*,2)]]
S_ref_all_I = [[I_cal]]
S_ref_all_Q = [[Q_cal]]
S_ref_all_U = [[U_cal]]
S_ref_all_D = [[D_cal]]
S_ref_all_DoLP = [[DoLP_cal]]

S_ref_all = [[S_ref_all_I], [S_ref_all_Q], [S_ref_all_U],[S_ref_all_D],[S_ref_all_DoLP]] ;--solar_cal is added 2019-10-18

n_all_each = 1;(size(s_ref_all_I,/dimension))[1]
node = [0, n_all_each, n_all_each*2, n_all_each*3, n_all_each*4]
endif 

if i_prof eq 1 and i_brdf eq 1 then begin ;---target
;---State vector: sfcref (2) and AODprof (6)
S_ref_all_I = [[I_cal], [JS_brdf_LAB_read(*,*,0)],[JS_aodprof2_LAB_read(*,*,0)]]
S_ref_all_Q = [[Q_cal], [JS_brdf_LAB_read(*,*,1)],[JS_aodprof2_LAB_read(*,*,1)]]
S_ref_all_U = [[U_cal], [JS_brdf_LAB_read(*,*,2)],[JS_aodprof2_LAB_read(*,*,2)]]
S_ref_all = [[S_ref_all_I], [S_ref_all_Q], [S_ref_all_U]]
n_all_each = (size(s_ref_all_I,/dimension))[1]
node = [0, n_all_each, n_all_each*2] ; 0_0 (n=1)=I_cal; 1_9 (n=9)=brdf; 10_25 (n=16)=AODprof2; and same for Q and U
endif 

endif

n_var = (size(S_ref_all,/dimension))[1]

for i_var = 0L, n_var-1 do begin

;--original
; S_ref = S_ref_all(usable_idx,i_var) 
;--solar irradiance multiply to all Spectrum and jacobian
S_ref = S_ref_all(usable_idx,i_var) ;* solar_cal(usable_idx)

if i_var eq 0 then begin


; ;----modified convolution using ILS!! ---0.01 res
; readcol, '/home/mchoi/HIMAP/for_ILS/ILS_sinc_box_opd3_idl.txt', ILS_read, /silent
; ;----sampling resolution again.
; wn_ILS = -5.00d + findgen(1001)*0.01d
; ; idx_sampling=where((round(wn_ILS*100d) mod 2) eq 0,n_idx_sampling)
; ; wn_ILS=wn_ILS(idx_sampling)
; ; ILS_read = ILS_read(idx_sampling)

;----modified convolution using ILS!! ---0.02 res (old)
; readcol, '/home/mchoi/HIMAP/for_ILS/ILS_sinc_box_opd3_idl_res0p02_n501.txt', ILS_read, /silent
; ;----sampling resolution again.
; wn_ILS = -5.00d + findgen(501)*0.02d



;----modified convolution using ILS!! ---0.02 res (2019-10-17 updated)
; ILS_pwd = '/home/mchoi/HIMAP/for_ILS_0p02/n3001_operational'
; if i_band eq 0 then ILS_file = file_search(ILS_pwd+'/ILS_sinc_box_idl_res0p02_n3001_operational_opd05.200_FWHM00.12_O2A.txt')
; if i_band eq 1 then ILS_file = file_search(ILS_pwd+'/ILS_sinc_box_idl_res0p02_n3001_operational_opd05.200_FWHM00.12_O2D.txt')

; ;----modified convolution using ILS!! ---0.10 res (2019-11-05 updated)
ILS_pwd = '/home/mchoi/HIMAP/for_ILS_0p10/n301_operational'
if i_band eq 0 then ILS_file = file_search(ILS_pwd+'/ILS_sinc_box_idl_res0p10_n301_operational_opd01.200_FWHM00.50_O2A.txt')
if i_band eq 1 then ILS_file = file_search(ILS_pwd+'/ILS_sinc_box_idl_res0p10_n301_operational_opd01.200_FWHM00.50_O2D.txt')

readcol, ILS_file(0), ILS_read, /silent
n_point_ILS = n_Elements(ILS_read)


; ILS_file = '/home/mchoi/HIMAP/for_ILS_0p02/ILS_sinc_box_idl_res0p02_n501_'+str_opd_FWHM_band+'.txt'
; readcol, ILS_file, ILS_read, /silent

;----sampling resolution again.
; wn_ILS = -5.00d + findgen(501)*0.02d
; wn_ILS = findgen(n_point_ILS)*0.02d - max(findgen(n_point_ILS)*0.02d)/2.0
; wn_within = max(findgen(n_point_ILS)*0.02d)/2.0
wn_ILS = findgen(n_point_ILS)*0.10d - max(findgen(n_point_ILS)*0.10d)/2.0
wn_within = max(findgen(n_point_ILS)*0.10d)/2.0



; idx_sampling=where((round(wn_ILS*100d) mod 2) eq 0,n_idx_sampling)
; wn_ILS=wn_ILS(idx_sampling)
; ILS_read = ILS_read(idx_sampling)


ILS_norm = ILS_read / total(ILS_read) ;normalization
n_point_ILS = n_elements(ILS_norm)



;----reduced wn_ILS and n_point_ILS
; wn_ILS = wn_ILS(500:500+100)






ILS_fn2 = I_cal*0.0d 
ILS_fn2(0:n_point_ILS-1) = ILS_norm(0:n_point_ILS-1)

IFG_ILS = FFT(ILS_fn2, -1, /center)
;---convolution in the level of interferogram
idx_inner = where(wn_cal ge wn_range_fin(0) and wn_cal le wn_range_fin(1))
endif ;if i_var eq 0 then begin

switch_FFT = 1

if switch_FFT eq 1 then begin
IFG_orig = FFT(S_ref,-1,/center)
IFG_conv_ILS = n_elements(IFG_orig) * IFG_orig * IFG_ILS
var_conv3 = FFT(IFG_conv_ILS, +1, /center)
idx_usable = idx_inner + (n_point_ILS-1)/2;shifting n/2 points of gaussian
var_conv3 = var_conv3(idx_usable)
;---interpolation to the measurment grid
S_conv = interpol(var_conv3,wn_inner,wn_conv)
endif

; if switch_FFT eq 0 then begin
; ;----Just convolution; not using FFT
idx_inner = where(wn_cal ge wn_range_fin(0) and wn_cal le wn_range_fin(1), n_idx_inner)
var_conv4 = fltarr(n_idx_inner) * !values.f_nan
for i=0L, n_idx_inner-1 do begin
var_conv4(i) = total(ILS_norm(*) * S_ref(idx_inner(i)-(n_elements(ILS_norm)-1)/2.:idx_inner(i)+(n_elements(ILS_norm)-1)/2.),/nan)
endfor
; S_conv = interpol(var_conv4,wn_inner,wn_conv)
; endif


if i_sin eq 0 then begin
  if i_var eq 0 then I_cal_SVO(*)=S_conv
  if i_var eq 1 then Q_cal_SVO(*)=S_conv
  if i_var eq 2 then U_cal_SVO(*)=S_conv
  if i_var eq 3 then Solar_cal_conv(*)=S_conv
endif 

if i_sin eq 1 then begin
  ; ;----I
  ; if i_var eq node(0) then I_cal_LAB(*)=S_conv
  ; if i_var ge node(0)+1 and i_var le node(0)+1+(n_prop_aerosol_LAB-1) then JS_aerosol_LAB(*,i_var-1,0)=S_conv
  ; if i_var ge node(0)+1+(n_prop_aerosol_LAB) and $
  ;    i_var le node(0)+1+(n_prop_aerosol_LAB+n_prop_sfcpres_LAB-1) then JS_sfcpres_LAB(*,i_var-(1+n_prop_aerosol_LAB),0)=S_conv
  ; if i_var ge node(0)+1+(n_prop_aerosol_LAB+n_prop_sfcpres_LAB) and $
  ;    i_var le node(0)+1+(n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_sfcref_LAB-1) then $
  ;    JS_sfcref_LAB(*,i_var-(1+n_prop_aerosol_LAB+n_prop_sfcpres_LAB),0)=S_conv
  ; ;----Q
  ; if i_var eq node(1) then Q_cal_LAB(*)=S_conv
  ; if i_var ge node(1)+1 and i_var le node(1)+1+(n_prop_aerosol_LAB-1) then JS_aerosol_LAB(*,i_var-1-node(1),1)=S_conv
  ; if i_var ge node(1)+1+(n_prop_aerosol_LAB) and $
  ;    i_var le node(1)+1+(n_prop_aerosol_LAB+n_prop_sfcpres_LAB-1) then JS_sfcpres_LAB(*,i_var-(1+n_prop_aerosol_LAB)-node(1),1)=S_conv
  ; if i_var ge node(1)+1+(n_prop_aerosol_LAB+n_prop_sfcpres_LAB) and $
  ;    i_var le node(1)+1+(n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_sfcref_LAB-1) then $
  ;    JS_sfcref_LAB(*,i_var-(1+n_prop_aerosol_LAB+n_prop_sfcpres_LAB)-node(1),1)=S_conv
  ; ;----U
  ; if i_var eq node(2) then U_cal_LAB(*)=S_conv
  ; if i_var ge node(2)+1 and i_var le node(2)+1+(n_prop_aerosol_LAB-1) then JS_aerosol_LAB(*,i_var-1-node(2),2)=S_conv
  ; if i_var ge node(2)+1+(n_prop_aerosol_LAB) and $
  ;    i_var le node(2)+1+(n_prop_aerosol_LAB+n_prop_sfcpres_LAB-1) then JS_sfcpres_LAB(*,i_var-(1+n_prop_aerosol_LAB)-node(2),2)=S_conv
  ; if i_var ge node(2)+1+(n_prop_aerosol_LAB+n_prop_sfcpres_LAB) and $
  ;    i_var le node(2)+1+(n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_sfcref_LAB-1) then $
  ;    JS_sfcref_LAB(*,i_var-(1+n_prop_aerosol_LAB+n_prop_sfcpres_LAB)-node(2),2)=S_conv

  if i_prof eq 0 and i_brdf eq 0  then begin ;---Aerosol Bulk
    ;----I
    if i_var eq node(0) then I_cal_LAB(*)=S_conv
    if i_var ge node(0)+1 and i_var le node(0)+1+(n_prop_sfcref_LAB-1) then JS_sfcref_LAB(*,i_var-1,0)=S_conv
    if i_var ge node(0)+1+(n_prop_sfcref_LAB) and $
      i_var le node(0)+1+(n_prop_sfcref_LAB+n_prop_aerosol_LAB-1) then JS_aerosol_LAB(*,i_var-(1+n_prop_sfcref_LAB),0)=S_conv
    ;----Q
    if i_var eq node(1) then Q_cal_LAB(*)=S_conv
    if i_var ge node(1)+1 and i_var le node(1)+1+(n_prop_sfcref_LAB-1) then JS_sfcref_LAB(*,i_var-1-node(1),1)=S_conv
    if i_var ge node(1)+1+(n_prop_sfcref_LAB) and $
      i_var le node(1)+1+(n_prop_sfcref_LAB+n_prop_aerosol_LAB-1) then JS_aerosol_LAB(*,i_var-(1+n_prop_sfcref_LAB)-node(1),1)=S_conv
    ;----U
    if i_var eq node(2) then U_cal_LAB(*)=S_conv
    if i_var ge node(2)+1 and i_var le node(2)+1+(n_prop_sfcref_LAB-1) then JS_sfcref_LAB(*,i_var-1-node(2),2)=S_conv
    if i_var ge node(2)+1+(n_prop_sfcref_LAB) and $
      i_var le node(2)+1+(n_prop_sfcref_LAB+n_prop_aerosol_LAB-1) then JS_aerosol_LAB(*,i_var-(1+n_prop_sfcref_LAB)-node(2),2)=S_conv
  endif ;if i_prof eq 1 then begin
  
  if i_prof eq 1 and i_brdf eq 0 then begin ;---AOD profile
    ;----I
    if i_var eq node(0) then I_cal_LAB(*)=S_conv
    if i_var ge node(0)+1 and i_var le node(0)+1+(n_prop_sfcref_LAB-1) then JS_sfcref_LAB(*,i_var-1,0)=S_conv
    if i_var ge node(0)+1+(n_prop_sfcref_LAB) and $
      i_var le node(0)+1+(n_prop_sfcref_LAB+n_prop_aodprof2_LAB-1) then JS_aodprof2_LAB(*,i_var-(1+n_prop_sfcref_LAB),0)=S_conv
    ;----Q
    if i_var eq node(1) then Q_cal_LAB(*)=S_conv
    if i_var ge node(1)+1 and i_var le node(1)+1+(n_prop_sfcref_LAB-1) then JS_sfcref_LAB(*,i_var-1-node(1),1)=S_conv
    if i_var ge node(1)+1+(n_prop_sfcref_LAB) and $
      i_var le node(1)+1+(n_prop_sfcref_LAB+n_prop_aodprof2_LAB-1) then JS_aodprof2_LAB(*,i_var-(1+n_prop_sfcref_LAB)-node(1),1)=S_conv
    ;----U
    if i_var eq node(2) then U_cal_LAB(*)=S_conv
    if i_var ge node(2)+1 and i_var le node(2)+1+(n_prop_sfcref_LAB-1) then JS_sfcref_LAB(*,i_var-1-node(2),2)=S_conv
    if i_var ge node(2)+1+(n_prop_sfcref_LAB) and $
      i_var le node(2)+1+(n_prop_sfcref_LAB+n_prop_aodprof2_LAB-1) then JS_aodprof2_LAB(*,i_var-(1+n_prop_sfcref_LAB)-node(2),2)=S_conv
  endif 

  if i_prof eq 0 and i_brdf eq 1  then begin ;---Aerosol Bulk; BRDF ;---second target
    ;----I
    if i_var eq node(0) then I_cal_LAB(*)=S_conv
    ; if i_var ge node(0)+1 and i_var le node(0)+1+(n_prop_BRDF_LAB-1) then $
    ;    JS_brdf_LAB(*,i_var-1,0)=S_conv
    ; if i_var ge node(0)+1+(n_prop_BRDF_LAB) and $
    ;    i_var le node(0)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB-1) then $
    ;    JS_aerosol_LAB(*,i_var-(1+n_prop_BRDF_LAB),0)=S_conv
    ; if i_var ge node(0)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB) and $
    ;    i_var le node(0)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB-1) then $
    ;    JS_sfcpres_LAB(*,i_var-(1+n_prop_BRDF_LAB+n_prop_aerosol_LAB),0)=S_conv
    ; if i_var ge node(0)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB) and $
    ;    i_var le node(0)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB-1) then $
    ;    JS_h2oscaling_LAB(*,i_var-(1+n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB),0)=S_conv
    ; if i_var ge node(0)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB) and $
    ;    i_var le node(0)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB+n_prop_tshift_LAB-1) then $
    ;    JS_tshift_LAB(*,i_var-(1+n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB),0)=S_conv            
    ;----Q
    if i_var eq node(1) then Q_cal_LAB(*)=S_conv
    ; if i_var ge node(1)+1 and i_var le node(1)+1+(n_prop_BRDF_LAB-1) then $
    ;    JS_brdf_LAB(*,i_var-1-node(1),1)=S_conv
    ; if i_var ge node(1)+1+(n_prop_BRDF_LAB) and $
    ;   i_var le node(1)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB-1) then $
    ;   JS_aerosol_LAB(*,i_var-(1+n_prop_BRDF_LAB)-node(1),1)=S_conv
    ; if i_var ge node(1)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB) and $
    ;    i_var le node(1)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB-1) then $
    ;    JS_sfcpres_LAB(*,i_var-(1+n_prop_BRDF_LAB+n_prop_aerosol_LAB)-node(1),1)=S_conv
    ; if i_var ge node(1)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB) and $
    ;    i_var le node(1)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB-1) then $
    ;    JS_h2oscaling_LAB(*,i_var-(1+n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB)-node(1),1)=S_conv
    ; if i_var ge node(1)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB) and $
    ;    i_var le node(1)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB+n_prop_tshift_LAB-1) then $
    ;    JS_tshift_LAB(*,i_var-(1+n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB)-node(1),1)=S_conv                  
    ;----U
    if i_var eq node(2) then U_cal_LAB(*)=S_conv
    ; if i_var ge node(2)+1 and i_var le node(2)+1+(n_prop_BRDF_LAB-1) then $
    ;    JS_brdf_LAB(*,i_var-1-node(2),2)=S_conv
    ; if i_var ge node(2)+1+(n_prop_BRDF_LAB) and $
    ;    i_var le node(2)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB-1) then $
    ;    JS_aerosol_LAB(*,i_var-(1+n_prop_BRDF_LAB)-node(2),2)=S_conv
    ; if i_var ge node(2)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB) and $
    ;    i_var le node(2)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB-1) then $
    ;    JS_sfcpres_LAB(*,i_var-(1+n_prop_BRDF_LAB+n_prop_aerosol_LAB)-node(2),2)=S_conv
    ; if i_var ge node(2)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB) and $
    ;    i_var le node(2)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB-1) then $
    ;    JS_h2oscaling_LAB(*,i_var-(1+n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB)-node(2),2)=S_conv             
    ; if i_var ge node(2)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB) and $
    ;    i_var le node(2)+1+(n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB+n_prop_tshift_LAB-1) then $
    ;    JS_tshift_LAB(*,i_var-(1+n_prop_BRDF_LAB+n_prop_aerosol_LAB+n_prop_sfcpres_LAB+n_prop_h2oscaling_LAB)-node(2),2)=S_conv             
    if i_var eq node(3) then D_cal_LAB(*)=S_conv
    if i_var eq node(4) then DoLP_cal_LAB(*)=S_conv
    ; if i_var eq node(4) then stop

 endif ;if i_prof eq 1 then begin

  if i_prof eq 1 and i_brdf eq 1 then begin ;---AOD profile; BRDF ->>> Target
    ;----I
    if i_var eq node(0) then I_cal_LAB(*)=S_conv
    if i_var ge node(0)+1 and i_var le node(0)+1+(n_prop_BRDF_LAB-1) then JS_BRDF_LAB(*,i_var-1,0)=S_conv
    if i_var ge node(0)+1+(n_prop_BRDF_LAB) and $
      i_var le node(0)+1+(n_prop_BRDF_LAB+n_prop_aodprof2_LAB-1) then JS_aodprof2_LAB(*,i_var-(1+n_prop_BRDF_LAB),0)=S_conv
    ;----Q
    if i_var eq node(1) then Q_cal_LAB(*)=S_conv
    if i_var ge node(1)+1 and i_var le node(1)+1+(n_prop_BRDF_LAB-1) then JS_brdf_LAB(*,i_var-1-node(1),1)=S_conv
    if i_var ge node(1)+1+(n_prop_BRDF_LAB) and $
      i_var le node(1)+1+(n_prop_BRDF_LAB+n_prop_aodprof2_LAB-1) then JS_aodprof2_LAB(*,i_var-(1+n_prop_BRDF_LAB)-node(1),1)=S_conv
    ;----U
    if i_var eq node(2) then U_cal_LAB(*)=S_conv
    if i_var ge node(2)+1 and i_var le node(2)+1+(n_prop_BRDF_LAB-1) then JS_brdf_LAB(*,i_var-1-node(2),2)=S_conv
    if i_var ge node(2)+1+(n_prop_BRDF_LAB) and $
      i_var le node(2)+1+(n_prop_BRDF_LAB+n_prop_aodprof2_LAB-1) then JS_aodprof2_LAB(*,i_var-(1+n_prop_BRDF_LAB)-node(2),2)=S_conv
  endif 

endif   

endfor ;for i_var = 0L, n_var-1 do begin
endfor ; for i=0L, n_sin-1 do begin



;---shifted wn is set as final points
wn_conv_fin = wn_conv ; - wn_shift_avg ;---always confusing about direction...

;---calculated Values are interpolated to the final wn points (shifting involved)
; Solar_cal_fin = interpol(solar_cal_conv,wn_conv,wn_conv_fin)
; I_cal_SVO_fin = interpol(I_cal_SVO,wn_conv,wn_conv_fin) 
I_cal_LAB_fin = interpol(I_cal_LAB,wn_conv,wn_conv_fin)

; Q_cal_SVO_fin  = interpol(Q_cal_SVO,wn_conv,wn_conv_fin) 
Q_cal_LAB_fin = interpol(Q_cal_LAB,wn_conv,wn_conv_fin)

; U_cal_SVO_fin  = interpol(U_cal_SVO,wn_conv,wn_conv_fin) 
U_cal_LAB_fin = interpol(U_cal_LAB,wn_conv,wn_conv_fin)

; D = DoLP 
; D_cal_LAB_fin = sqrt((Q_cal_LAB_fin*Q_cal_LAB_fin+U_cal_LAB_fin*U_cal_LAB_fin)/(I_cal_LAB_fin*I_cal_LAB_fin))

;--altenativley, D = Ipol = sqrt(Q^2+U^2), not divided by I
; D_cal_LAB_fin = sqrt(Q_cal_LAB_fin*Q_cal_LAB_fin+U_cal_LAB_fin*U_cal_LAB_fin)
D_cal_LAB_fin = interpol(D_cal_LAB,wn_conv,wn_conv_fin)
; DoLP_cal_LAB_fin = interpol(DoLP_cal_LAB,wn_conv,wn_conv_fin)



; JS_aerosol_LAB_fin = JS_aerosol_LAB*!values.d_nan
; JS_sfcpres_LAB_fin = JS_sfcpres_LAB*!values.d_nan
; JS_H2Oscaling_LAB_fin = JS_H2Oscaling_LAB*!values.d_nan
; JS_sfcref_LAB_fin    = JS_sfcref_LAB*!values.d_nan
; JS_aodprof2_LAB_fin  = JS_aodprof2_LAB*!values.d_nan
; JS_brdf_LAB_fin  = JS_brdf_LAB*!values.d_nan
; JS_tshift_LAB_fin = JS_tshift_LAB*!values.d_nan


; Q2U2 = (Q_cal_LAB_fin * Q_cal_LAB_fin + U_cal_LAB_fin * U_cal_LAB_fin)

; for i_stokes=0, 2 do begin

; if i_brdf eq 0 then begin
;   for i_var = 0L, n_prop_sfcref_LAB-1 do begin
;     JS_sfcref_LAB_fin(*,i_var,i_stokes)=interpol(JS_sfcref_LAB(*,i_var,i_stokes),wn_conv,wn_conv_fin)
;     if i_stokes eq 2 then begin
;       JS_I = JS_sfcref_LAB_fin(*,i_var,0)
;       JS_Q = JS_sfcref_LAB_fin(*,i_var,1)
;       JS_U = JS_sfcref_LAB_fin(*,i_var,2)
;       ; JS_sfcref_LAB_fin(*,i_var,3) = $
;       ;               1./(2.*I_cal_LAB_fin)*(Q2U2^(-0.5))*(2.*Q_cal_LAB_fin*JS_Q+2.*U_cal_LAB_fin*JS_U) + $
;       ;               (-1.0)*(Q2U2^(0.5))/(I_cal_LAB_fin^2.0)*JS_I
;       JS_sfcref_LAB_fin(*,i_var,3) = (Q_cal_LAB_fin * JS_Q + U_cal_LAB_fin * JS_U)/sqrt(Q2U2)
;     endif
;   endfor
; endif

; if i_brdf eq 1 then begin
;   for i_var = 0L, n_prop_brdf_LAB-1 do begin
;     JS_brdf_LAB_fin(*,i_var,i_stokes)=interpol(JS_brdf_LAB(*,i_var,i_stokes),wn_conv,wn_conv_fin)
;     if i_stokes eq 2 then begin
;       JS_I = JS_brdf_LAB_fin(*,i_var,0)
;       JS_Q = JS_brdf_LAB_fin(*,i_var,1)
;       JS_U = JS_brdf_LAB_fin(*,i_var,2)
;       ; JS_brdf_LAB_fin(*,i_var,3) = $
;       ;               1./(2.*I_cal_LAB_fin)*(Q2U2^(-0.5))*(2.*Q_cal_LAB_fin*JS_Q+2.*U_cal_LAB_fin*JS_U) + $
;       ;               (-1.0)*(Q2U2^(0.5))/(I_cal_LAB_fin^2.0)*JS_I
;       JS_brdf_LAB_fin(*,i_var,3) = (Q_cal_LAB_fin * JS_Q + U_cal_LAB_fin * JS_U)/sqrt(Q2U2)
;     endif
;   endfor
; endif

; if i_prof eq 0 then begin
;   for i_var = 0L, n_prop_aerosol_LAB-1 do begin
;     ;JS_aerosol_LAB_fin(*,i_var,i_stokes)=interpol(JS_aerosol_LAB(*,i_var,i_stokes),wn_conv,wn_conv_fin) 
;     JS_aerosol_LAB_fin(*,i_var,i_stokes)=interpol(JS_aerosol_LAB(*,i_var,i_stokes),wn_conv,wn_conv_fin) ;* prop_aerosol_LAB(i_var) ;--alog
;     if i_stokes eq 2 then begin
;       JS_I = JS_aerosol_LAB_fin(*,i_var,0)
;       JS_Q = JS_aerosol_LAB_fin(*,i_var,1)
;       JS_U = JS_aerosol_LAB_fin(*,i_var,2)
;       ; JS_aerosol_LAB_fin(*,i_var,3) = $
;       ;               1./(2.*I_cal_LAB_fin)*(Q2U2^(-0.5))*(2.*Q_cal_LAB_fin*JS_Q+2.*U_cal_LAB_fin*JS_U) + $
;       ;               (-1.0)*(Q2U2^(0.5))/(I_cal_LAB_fin^2.0)*JS_I
;       JS_aerosol_LAB_fin(*,i_var,3) = (Q_cal_LAB_fin * JS_Q + U_cal_LAB_fin * JS_U)/sqrt(Q2U2)              
;     endif
;   endfor
;   for i_var = 0L, n_prop_sfcpres_LAB-1 do begin
;     JS_sfcpres_LAB_fin(*,i_var,i_stokes)=interpol(JS_sfcpres_LAB(*,i_var,i_stokes),wn_conv,wn_conv_fin)
;     if i_stokes eq 2 then begin
;       JS_I = JS_sfcpres_LAB_fin(*,i_var,0)
;       JS_Q = JS_sfcpres_LAB_fin(*,i_var,1)
;       JS_U = JS_sfcpres_LAB_fin(*,i_var,2)
;       ; JS_sfcpres_LAB_fin(*,i_var,3) = $
;       ;               1./(2.*I_cal_LAB_fin)*(Q2U2^(-0.5))*(2.*Q_cal_LAB_fin*JS_Q+2.*U_cal_LAB_fin*JS_U) + $
;       ;               (-1.0)*(Q2U2^(0.5))/(I_cal_LAB_fin^2.0)*JS_I
;       JS_sfcpres_LAB_fin(*,i_var,3) = (Q_cal_LAB_fin * JS_Q + U_cal_LAB_fin * JS_U)/sqrt(Q2U2)
                    
;     endif
;   endfor
;   for i_var = 0L, n_prop_H2Oscaling_LAB-1 do begin
;     JS_H2Oscaling_LAB_fin(*,i_var,i_stokes)=interpol(JS_H2Oscaling_LAB(*,i_var,i_stokes),wn_conv,wn_conv_fin)
;     if i_stokes eq 2 then begin
;       JS_I = JS_H2Oscaling_LAB_fin(*,i_var,0)
;       JS_Q = JS_H2Oscaling_LAB_fin(*,i_var,1)
;       JS_U = JS_H2Oscaling_LAB_fin(*,i_var,2)
;       ; JS_H2Oscaling_LAB_fin(*,i_var,3) = $
;       ;               1./(2.*I_cal_LAB_fin)*(Q2U2^(-0.5))*(2.*Q_cal_LAB_fin*JS_Q+2.*U_cal_LAB_fin*JS_U) + $
;       ;               (-1.0)*(Q2U2^(0.5))/(I_cal_LAB_fin^2.0)*JS_I
;       JS_H2Oscaling_LAB_fin(*,i_var,3) = (Q_cal_LAB_fin * JS_Q + U_cal_LAB_fin * JS_U)/sqrt(Q2U2)
;     endif
;   endfor
;   for i_var = 0L, n_prop_tshift_LAB-1 do begin
;     JS_tshift_LAB_fin(*,i_var,i_stokes)=interpol(JS_tshift_LAB(*,i_var,i_stokes),wn_conv,wn_conv_fin)
;     if i_stokes eq 2 then begin
;       JS_I = JS_tshift_LAB_fin(*,i_var,0)
;       JS_Q = JS_tshift_LAB_fin(*,i_var,1)
;       JS_U = JS_tshift_LAB_fin(*,i_var,2)
;       ; JS_tshift_LAB_fin(*,i_var,3) = $
;       ;               1./(2.*I_cal_LAB_fin)*(Q2U2^(-0.5))*(2.*Q_cal_LAB_fin*JS_Q+2.*U_cal_LAB_fin*JS_U) + $
;       ;               (-1.0)*(Q2U2^(0.5))/(I_cal_LAB_fin^2.0)*JS_I
;       JS_tshift_LAB_fin(*,i_var,3) = (Q_cal_LAB_fin * JS_Q + U_cal_LAB_fin * JS_U)/sqrt(Q2U2)
;     endif
;   endfor
; endif

; if i_prof eq 1 then begin
;   for i_var = 0L, n_prop_aodprof2_LAB-1 do begin
;     JS_aodprof2_LAB_fin(*,i_var,i_stokes)=interpol(JS_aodprof2_LAB(*,i_var,i_stokes),wn_conv,wn_conv_fin)
;     if i_stokes eq 2 then begin
;       JS_I = JS_aodprof2_LAB_fin(*,i_var,0)
;       JS_Q = JS_aodprof2_LAB_fin(*,i_var,1)
;       JS_U = JS_aodprof2_LAB_fin(*,i_var,2)
;       ; JS_aodprof2_LAB_fin(*,i_var,3) = $
;       ;               1./(2.*I_cal_LAB_fin)*(Q2U2^(-0.5))*(2.*Q_cal_LAB_fin*JS_Q+2.*U_cal_LAB_fin*JS_U) + $
;       ;               (-1.0)*(Q2U2^(0.5))/(I_cal_LAB_fin^2.0)*JS_I
;       JS_aodprof2_LAB_fin(*,i_var,3) = (Q_cal_LAB_fin * JS_Q + U_cal_LAB_fin * JS_U)/sqrt(Q2U2)
                    
;     endif
;   endfor
; endif


; endfor ;for i_stokes=0, 2 do begin

; ;--wavenumber shift calculation
; wn_Ref1 = wn_conv_fin
; I_ref1 = I_cal_LAB_fin
; wn_meas1 = wn_ref1
; I_meas1 = S000_read(idx_Select)

; dwn=abs(wn_meas1(1)-wn_meas1(0))*0.002
; n_ref2 = round(((max(wn_meas1)-min(wn_meas1))/dwn))+1
; wn_ref2 = min(wn_meas1)+findgen(n_ref2)*dwn
; I_ref2 = interpol(I_ref1,wn_Ref1,wn_Ref2)
; I_meas2 = interpol(I_meas1,wn_ref1,wn_ref2)
; lag= -1000.0+findgen(2001)
; C_Corr=c_correlate(I_ref2,I_meas2,lag,/double)
; max_idx=where(c_corr eq max(c_corr))
; wn_shift = lag(max_idx(0))*dwn


; ;---if no shift in wavenumber
; if i_band eq 0 then begin
;   wn_conv_O2A = wn_conv_fin
;   I_cal_LAB_O2A = I_cal_LAB_fin
;   Q_cal_LAB_O2A = Q_cal_LAB_fin
;   U_cal_LAB_O2A = U_cal_LAB_fin
;   D_cal_LAB_O2A = D_cal_LAB_fin
;   JS_aerosol_LAB_O2A = JS_aerosol_LAB_fin
;   JS_sfcref_LAB_O2A = JS_sfcref_LAB_fin
;   JS_aodprof2_LAB_O2A = JS_aodprof2_LAB_fin
;   JS_brdf_LAB_O2A = JS_brdf_LAB_fin
;   JS_sfcpres_LAB_O2A = JS_sfcpres_LAB_fin
;   JS_H2Oscaling_LAB_O2A = JS_H2Oscaling_LAB_fin
;   ;--measurement info
;   I_meas_LAB_O2A = S000_read(idx_Select)
; endif
; if i_band eq 1 then begin
;   wn_conv_O2D = wn_conv_fin
;   I_cal_LAB_O2D = I_cal_LAB_fin
;   Q_cal_LAB_O2D = Q_cal_LAB_fin
;   U_cal_LAB_O2D = U_cal_LAB_fin
;   D_cal_LAB_O2D = D_cal_LAB_fin
;   JS_aerosol_LAB_O2D = JS_aerosol_LAB_fin
;   JS_sfcref_LAB_O2D = JS_sfcref_LAB_fin
;   JS_aodprof2_LAB_O2D = JS_aodprof2_LAB_fin
;   JS_brdf_LAB_O2D = JS_brdf_LAB_fin
;   JS_sfcpres_LAB_O2D = JS_sfcpres_LAB_fin
;   JS_H2Oscaling_LAB_O2D = JS_H2Oscaling_LAB_fin
;   ;--measurement info
;   I_meas_LAB_O2D = S000_read(idx_Select)
; endif

;---reflecting shift
if i_band eq 0 then begin
  wn_conv_O2A = wn_conv_fin
  I_meas_LAB_O2A = I_cal_LAB_fin
  Q_meas_LAB_O2A = Q_cal_LAB_fin
  U_meas_LAB_O2A = U_cal_LAB_fin
  D_meas_LAB_O2A = D_cal_LAB_fin
  DoLP_meas_LAB_O2A = sqrt((Q_cal_LAB_fin^2.0+U_cal_LAB_fin^2.0)/I_cal_LAB_fin^2.0);DoLP_cal_LAB_fin
  ; I_cal_LAB_O2A = I_cal_LAB_fin
  ; Q_cal_LAB_O2A = Q_cal_LAB_fin
  ; U_cal_LAB_O2A = U_cal_LAB_fin
  ; D_cal_LAB_O2A = D_cal_LAB_fin
  ; JS_aerosol_LAB_O2A = JS_aerosol_LAB_fin
  ; JS_sfcref_LAB_O2A = JS_sfcref_LAB_fin
  ; JS_aodprof2_LAB_O2A = JS_aodprof2_LAB_fin
  ; JS_brdf_LAB_O2A = JS_brdf_LAB_fin
  ; JS_sfcpres_LAB_O2A = JS_sfcpres_LAB_fin
  ; JS_H2Oscaling_LAB_O2A = JS_H2Oscaling_LAB_fin
  ; JS_tshift_LAB_O2A = JS_tshift_LAB_fin
  ; ;--measurement info
  ; I_meas_LAB_O2A = interpol(S000_read(idx_Select),wn_conv_fin-wn_shift, wn_Conv_fin) ;--direction...
  ; D_meas_LAB_O2A = I_meas_LAB_O2A * !values.f_nan ;--should be updated later
endif
if i_band eq 1 then begin
  wn_conv_O2D = wn_conv_fin
  I_meas_LAB_O2D = I_cal_LAB_fin
  Q_meas_LAB_O2D = Q_cal_LAB_fin
  U_meas_LAB_O2D = U_cal_LAB_fin
  D_meas_LAB_O2D = D_cal_LAB_fin
  DoLP_meas_LAB_O2D = sqrt((Q_cal_LAB_fin^2.0+U_cal_LAB_fin^2.0)/I_cal_LAB_fin^2.0);DoLP_cal_LAB_fin
  ; I_cal_LAB_O2D = I_cal_LAB_fin
  ; Q_cal_LAB_O2D = Q_cal_LAB_fin
  ; U_cal_LAB_O2D = U_cal_LAB_fin
  ; D_cal_LAB_O2D = D_cal_LAB_fin
  ; JS_aerosol_LAB_O2D = JS_aerosol_LAB_fin
  ; JS_sfcref_LAB_O2D = JS_sfcref_LAB_fin
  ; JS_aodprof2_LAB_O2D = JS_aodprof2_LAB_fin
  ; JS_brdf_LAB_O2D = JS_brdf_LAB_fin
  ; JS_sfcpres_LAB_O2D = JS_sfcpres_LAB_fin
  ; JS_H2Oscaling_LAB_O2D = JS_H2Oscaling_LAB_fin
  ; JS_tshift_LAB_O2D = JS_tshift_LAB_fin
  ; ;--measurement info
  ; I_meas_LAB_O2D = interpol(S000_read(idx_Select),wn_conv_fin-wn_shift, wn_Conv_fin)
  ; D_meas_LAB_O2D = I_meas_LAB_O2D * !values.f_nan ;--should be updated later
endif


endfor ;for i_band = 0L, 1 do begin

;---Synthetic measurment database
;---Adding random error to I_cal_LAB_O2A, I_cal_LAB_O2D

;---save synthetic measurement data to xdr file

save, wn_conv_O2A, I_meas_LAB_O2A, Q_meas_LAB_O2A, U_meas_LAB_O2A, D_meas_LAB_O2A, DoLP_meas_LAB_O2A, $
      wn_conv_O2D, I_meas_LAB_O2D, Q_meas_LAB_O2D, U_meas_LAB_O2D, D_meas_LAB_O2D, DoLP_meas_LAB_O2D, $
      filename = SV_pwd1+'/synthetic_meas_O2AD.xdr'





; endfor ;for i_width = 0L, n_elements(width_set)-1 do begin
; endfor ;for i_peak = 0L, n_elements(Peak_set)-1 do begin
; endfor ;for i_AOD = 0L, n_elements(AOD_set)-1 do begin







end