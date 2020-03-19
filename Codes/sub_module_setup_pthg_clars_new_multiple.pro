pro sub_module_setup_pthg_clars_new_multiple, i_step, $ 
					idx_target, lat_r,lon_r,yy_r,mm_r,dd_r,hhh_r,mmm_r, $
					CLARS_alt_km_r, CLARS_pres_hPa_r, CLARS_temp_K_r, $
					LABS_alt_km_r, LABS_pres_hpa_r, workspace_pwd, tshift, $
					switch_synthetic, switch_perturb, switch_constraint_aerosol, $
					i_case, savefolder_num, $
					result_pwd,ii,i_loop_start,input_loop_begin

;--original
; pro module_setup_pthg_clars, idx_target, lat_r,lon_r,yy_r,mm_r,dd_r,hhh_r,mmm_r, $
; 					CLARS_alt_km_r, CLARS_pres_hPa_r, CLARS_temp_K_r, $
; 					LABS_alt_km_r, LABS_pres_hpa_r, LABS_temp_C_r, workspace_pwd, tshift

print, 'Process: making PTHG profiles...'


;Myungje Choi, 11-15-2018
;recalculation for NSW
;adopting "Initial_NSW_atmospehre_12aug13" file but more grid.

;Myungje Choi, 11-08-2018
;original PTH is from US standard atmosphere 1976; https://www.digitaldutch.com/atmoscalc/tableatmosphere.htm
;original Trace gases profile: "DiscoverAQ_profile_1.dat" which is a default option of GEMSTOOL.

;-modified from interpolation_pthg_for_nsw_using_measurement_o2ab.pro

; ;--save Pressure-Tempereture-Gas VMR
; file_title0 = 'atmos_mchoi'+string() ;--same with atmos_mchoi_t08
; file_title = file_title0+'.dat'


;----optaining pressure level info. from measurement
Lat_dum = lat_r ;34.201691
Lon_dum = lon_r; -118.1717
Year_dum = yy_r; 2018
Month_dum = mm_r ;12
DayofMonth_dum = dd_r ;??
Hour_dum = hhh_r;12+8 ;UTC
Min_dum = mmm_r;20.0

jd_obs=julday(Month_dum,DayofMonth_dum,Year_dum,Hour_dum,Min_dum,00)
jd0=julday(Month_dum,DayofMonth_dum,Year_dum) ;---real time matched
; jd_obs=julday(12,03,2018,20,00,00)
; jd0 = julday(12,03,2018) ;---just one sample -> applied to other days

;--Historical Surf_press_HG range of JPL = 992.11853hPa - 1027.9132hPa (29.38in - 30.44in)
; Surf_press_HG=30.44  ;29.38 ;30.44 ;30.06
; Surf_press_hPa=3376.85 * Surf_press_HG / 100. ;--unit conversion from inches of Hg to hPa
; Surf_temp_F=62.2
; Surf_temp_C=(Surf_temp_F-32.0)/1.8
; Surf_temp_K=Surf_temp_C+273.15
; Surf_elevation_km=0.100
; Surf_elevation_feet=1138.45

Surf_elevation_km = LABS_alt_km_r
; Surf_pres_hPa = LABS_pres_hPa_r ;--read from MERRA2 directly
; Surf_temp_C = LABS_temp_C_r ;--should be changed directly from MERRA2
; Surf_temp_K = Surf_temp_C+273.15 ;--should be changed directly from MERRA2

CLARS_temp_K = CLARS_temp_K_r ;CLARS_temp_C_r+273.15



; H_alt = (1.0-(Surf_press_hPa/1013.25)^0.190284)*145366.45

US1976_file = '/home/mchoi/OSSE_tool/Datasets/data_PTHG/data_originally_used/PTH_us_standard_1976_res100m_max.dat'
discoveraq_file = '/home/mchoi/OSSE_tool/Datasets/data_PTHG/data_originally_used/DiscoverAQ_profile_1.dat'

readcol, US1976_file, num_us, alt_us, Temp_us, Pres_us, /silent
readcol, discoveraq_file, num_daq, pres_daq, temp_daq, H2O_daq, O3_daq, NO2_daq, HCHO_daq, SO2_daq, CO_daq, /silent


;---finding sea-level pressure first from inteprolated MERRA-2 data

;---read MERRA-2 data 
;--location
; pwd_MERRA2='/home/mchoi/OSSE_tool/for_PTHG_setting/data_MERRA2'
pwd_MERRA2='/home/mchoi/OSSE_tool/Datasets/data_PTHG/data_MERRA2_dn_20191003'

caldat, jd0, mm0,dd0,yy0
caldat, jd0+1, mm1,dd1,yy1
str_yymmdd0=string(yy0,f='(i4.4)')+string(mm0,f='(i2.2)')+string(dd0,f='(i2.2)')
str_yymmdd1=string(yy1,f='(i4.4)')+string(mm1,f='(i2.2)')+string(dd1,f='(i2.2)')

file0=pwd_MERRA2+'/MERRA2_400.inst3_3d_asm_Nv.'+str_yymmdd0+'.nc4'
file1=pwd_MERRA2+'/MERRA2_400.inst3_3d_asm_Nv.'+str_yymmdd1+'.nc4'



;--read lat lon lev info (identical)
ncid=NCDF_open(file0)
ncdf_varget, ncid, 'lat',lat
ncdf_varget, ncid, 'lon',lon
ncdf_varget, ncid, 'lev',lev
ncdf_varget, ncid, 'time',time
NCDF_close, ncid
n_lat=n_elements(lat)
n_lon=n_elements(lon)
n_lev=n_elements(lev)
n_time=n_elements(time)

;--finding value_locate for lat, lon, time

;--read 0 info
ncid=NCDF_open(file0)
ncdf_varget, ncid, 'T', T0 ;temperature
ncdf_varget, ncid, 'H', H0 ;mid_layer_heights
ncdf_varget, ncid, 'PL', PL0 ;mid-level pressure
ncdf_varget, ncid, 'QV', QV0 ;Specific Humidity
ncdf_varget, ncid, 'PS', PS0 ;Surface pressure
ncdf_varget, ncid, 'SLP', SLP0 ;Sea-level pressure
NCDF_close, ncid

;--read 0 info
ncid=NCDF_open(file1)
ncdf_varget, ncid, 'T', T1 ;temperature
ncdf_varget, ncid, 'H', H1 ;mid_layer_heights
ncdf_varget, ncid, 'PL', PL1 ;mid-level pressure
ncdf_varget, ncid, 'QV', QV1 ;Specific Humidity
ncdf_varget, ncid, 'PS', PS1 ;Surface pressure
ncdf_varget, ncid, 'SLP', SLP1 ;Sea-level pressure
NCDF_close, ncid

;--two data integrated
n_time2=2*n_time
T2=fltarr(n_lon,n_lat,n_lev,n_time2)*!values.f_nan
H2=fltarr(n_lon,n_lat,n_lev,n_time2)*!values.f_nan
PL2=fltarr(n_lon,n_lat,n_lev,n_time2)*!values.f_nan
QV2=fltarr(n_lon,n_lat,n_lev,n_time2)*!values.f_nan
PS2=fltarr(n_lon,n_lat,n_time2)*!values.f_nan
SLP2=fltarr(n_lon,n_lat,n_time2)*!values.f_nan

;--spatio-temporal interpolated
T3=fltarr(n_lev)*!values.f_nan
H3=fltarr(n_lev)*!values.f_nan
PL3=fltarr(n_lev)*!values.f_nan
QV3=fltarr(n_lev)*!values.f_nan

T2(*,*,*,0:7)=T0(*,*,*,*)
H2(*,*,*,0:7)=H0(*,*,*,*)
PL2(*,*,*,0:7)=PL0(*,*,*,*)
QV2(*,*,*,0:7)=QV0(*,*,*,*)

T2(*,*,*,8:15)=T1(*,*,*,*)
H2(*,*,*,8:15)=H1(*,*,*,*)
PL2(*,*,*,8:15)=PL1(*,*,*,*)
QV2(*,*,*,8:15)=QV1(*,*,*,*)

PS2(*,*,0:7) =PS0(*,*,*)
PS2(*,*,8:15)=PS1(*,*,*)
SLP2(*,*,0:7) =SLP0(*,*,*)
SLP2(*,*,8:15)=SLP1(*,*,*)

jd_3hours=julday(1,1,1,03,00,00)-julday(1,1,1,00,00,00)
jd0=julday(mm0,dd0,yy0,00,00,00)+findgen(n_time2)*jd_3hours

;--finding spatiotemporal location
loc_lon=value_locate(lon,lon_dum)
loc_lat=value_locate(lat,lat_dum)
loc_jd=value_locate(jd0,jd_obs)

inv00=1./(map_2points(lon(loc_lon),lat(loc_lat),lon_dum,lat_dum,/meters)/1000.)^2.
inv01=1./(map_2points(lon(loc_lon),lat(loc_lat+1),lon_dum,lat_dum,/meters)/1000.)^2.
inv10=1./(map_2points(lon(loc_lon+1),lat(loc_lat),lon_dum,lat_dum,/meters)/1000.)^2.
inv11=1./(map_2points(lon(loc_lon+1),lat(loc_lat+1),lon_dum,lat_dum,/meters)/1000.)^2.

total_inv=inv00+inv01+inv10+inv11
wgt00=inv00/total_inv
wgt01=inv01/total_inv
wgt10=inv10/total_inv
wgt11=inv11/total_inv


for i_lev=0L, n_lev-1 do begin
;--Temperature: T2 to T3
temp0=wgt00*T2(loc_lon,loc_lat,i_lev,loc_jd)+wgt01*T2(loc_lon,loc_lat+1,i_lev,loc_jd)+ $
		wgt10*T2(loc_lon+1,loc_lat,i_lev,loc_jd)+wgt11*T2(loc_lon+1,loc_lat+1,i_lev,loc_jd)
temp1=wgt00*T2(loc_lon,loc_lat,i_lev,loc_jd+1)+wgt01*T2(loc_lon,loc_lat+1,i_lev,loc_jd+1)+ $
		wgt10*T2(loc_lon+1,loc_lat,i_lev,loc_jd+1)+wgt11*T2(loc_lon+1,loc_lat+1,i_lev,loc_jd+1)		
T3(i_lev)=interpol([temp0,temp1],[jd0(loc_jd),jd0(loc_jd+1)],jd_obs)
;--Height: H2 to H3
temp0=wgt00*H2(loc_lon,loc_lat,i_lev,loc_jd)+wgt01*H2(loc_lon,loc_lat+1,i_lev,loc_jd)+ $
		wgt10*H2(loc_lon+1,loc_lat,i_lev,loc_jd)+wgt11*H2(loc_lon+1,loc_lat+1,i_lev,loc_jd)
temp1=wgt00*H2(loc_lon,loc_lat,i_lev,loc_jd+1)+wgt01*H2(loc_lon,loc_lat+1,i_lev,loc_jd+1)+ $
		wgt10*H2(loc_lon+1,loc_lat,i_lev,loc_jd+1)+wgt11*H2(loc_lon+1,loc_lat+1,i_lev,loc_jd+1)		
H3(i_lev)=interpol([temp0,temp1],[jd0(loc_jd),jd0(loc_jd+1)],jd_obs)
;--Pressure: PL2 to PL3
temp0=wgt00*PL2(loc_lon,loc_lat,i_lev,loc_jd)+wgt01*PL2(loc_lon,loc_lat+1,i_lev,loc_jd)+ $
		wgt10*PL2(loc_lon+1,loc_lat,i_lev,loc_jd)+wgt11*PL2(loc_lon+1,loc_lat+1,i_lev,loc_jd)
temp1=wgt00*PL2(loc_lon,loc_lat,i_lev,loc_jd+1)+wgt01*PL2(loc_lon,loc_lat+1,i_lev,loc_jd+1)+ $
		wgt10*PL2(loc_lon+1,loc_lat,i_lev,loc_jd+1)+wgt11*PL2(loc_lon+1,loc_lat+1,i_lev,loc_jd+1)		
PL3(i_lev)=interpol([temp0,temp1],[jd0(loc_jd),jd0(loc_jd+1)],jd_obs)
;--Specific Humidity: QV2 to QV3
temp0=wgt00*QV2(loc_lon,loc_lat,i_lev,loc_jd)+wgt01*QV2(loc_lon,loc_lat+1,i_lev,loc_jd)+ $
		wgt10*QV2(loc_lon+1,loc_lat,i_lev,loc_jd)+wgt11*QV2(loc_lon+1,loc_lat+1,i_lev,loc_jd)
temp1=wgt00*QV2(loc_lon,loc_lat,i_lev,loc_jd+1)+wgt01*QV2(loc_lon,loc_lat+1,i_lev,loc_jd+1)+ $
		wgt10*QV2(loc_lon+1,loc_lat,i_lev,loc_jd+1)+wgt11*QV2(loc_lon+1,loc_lat+1,i_lev,loc_jd+1)		
QV3(i_lev)=interpol([temp0,temp1],[jd0(loc_jd),jd0(loc_jd+1)],jd_obs)
endfor

;--Surface pressure: PS2 to PS3
temp0=wgt00*PS2(loc_lon,loc_lat,loc_jd)+wgt01*PS2(loc_lon,loc_lat+1,loc_jd)+ $
		wgt10*PS2(loc_lon+1,loc_lat,loc_jd)+wgt11*PS2(loc_lon+1,loc_lat+1,loc_jd)
temp1=wgt00*PS2(loc_lon,loc_lat,loc_jd+1)+wgt01*PS2(loc_lon,loc_lat+1,loc_jd+1)+ $
		wgt10*PS2(loc_lon+1,loc_lat,loc_jd+1)+wgt11*PS2(loc_lon+1,loc_lat+1,loc_jd+1)		
PS3=interpol([temp0,temp1],[jd0(loc_jd),jd0(loc_jd+1)],jd_obs)

;--Sea-level pressure: SLP2 to SLP3
temp0=wgt00*SLP2(loc_lon,loc_lat,loc_jd)+wgt01*SLP2(loc_lon,loc_lat+1,loc_jd)+ $
		wgt10*SLP2(loc_lon+1,loc_lat,loc_jd)+wgt11*SLP2(loc_lon+1,loc_lat+1,loc_jd)
temp1=wgt00*SLP2(loc_lon,loc_lat,loc_jd+1)+wgt01*SLP2(loc_lon,loc_lat+1,loc_jd+1)+ $
		wgt10*SLP2(loc_lon+1,loc_lat,loc_jd+1)+wgt11*SLP2(loc_lon+1,loc_lat+1,loc_jd+1)		
SLP3=interpol([temp0,temp1],[jd0(loc_jd),jd0(loc_jd+1)],jd_obs)

;---Surface pressure and temperature directly from MERRA-2
;-----for synthetic data!!!!  -> switch_synthetic=1
;-----switch_perturb=0; no random noise for surface pressure.
;-----switch_perturb=1; random noise added to surface pressure.
if i_step eq 0 and switch_synthetic eq 1 and switch_perturb eq 0 then Surf_pres_hPa = PS3/100.
if i_step eq 0 and switch_synthetic eq 1 and switch_perturb eq 1 then Surf_pres_hPa = PS3/100. + randomu(seed, 1,/normal) * 3.0 ;--pa -> hPa ;--1st step: a priori from MERRA-2;  ;+-3K from TROPOMI ATBD
;-----for retrieval process!!!-> switch_synthetic=0
;-----switch_constraint_aerosol = 0 or 1;  -> just starting from a priori MERRA-2 values
;-----switch_constraint_aerosol = 2;  -> from true value of synthetic database 
if i_step eq 0 and switch_synthetic eq 0 and switch_constraint_aerosol le 1 then Surf_pres_hPa = PS3/100. 
if i_step eq 0 and switch_synthetic eq 0 and switch_constraint_aerosol eq 2 then begin
	str_savefolder_fin = strtrim(string(savefolder_num),2)
      str_sub = 'Case'+string(i_case,f='(i3.3)')
      if switch_perturb eq 0 then begin ;--fixed AOD, PH, HW
      Synthetic_pwd = file_search(workspace_pwd+'/Results/All_Synthetic_fixed_cl'+str_savefolder_fin+'*/'+str_sub+'*/state_vectors_step')
      endif
      if switch_perturb eq 1 then begin ;--random all pars
      Synthetic_pwd = file_search(workspace_pwd+'/Results/All_Synthetic_random_cl'+str_savefolder_fin+'*/'+str_sub+'*/state_vectors_step')
      endif
      Synthetic_pwd = Synthetic_pwd(0)
      file_SV_t_sfcpres = file_search(Synthetic_pwd+'/state_vectors_surfpress_step00.txt')
      readcol, file_SV_t_sfcpres(0), SV_t_sfcpres, format='d',/silent
	Surf_pres_hPa = SV_t_sfcpres
endif


if i_step gt 0 then Surf_pres_hPa = LABS_pres_hpa_r

;---surface temperature from intepolation using MERRA-2
Surf_temp_K = interpol(T3, PL3, PS3)


;--setting pres_fin profile ;----------previous 25 layers / 26 levels
; alt_us = alt_us * 0.001
; 	n_interval0 = 10
; 	n_interval1 = 90
; 	interval0 = (Surf_pres_hPa-clars_pres_hPa_r)/(n_interval0)
; 	interval1 = (clars_pres_hPa_r-10)/(n_interval1)
; 	Pres_fin = [Surf_pres_hPa-findgen(n_interval0)*interval0, clars_pres_hPa_r-findgen(n_interval1)*interval1, $
; 			10, 5, 1, 0.5, 0.1]
			

; 	node = [Surf_elevation_km, CLARS_alt_km_r, 15, 60]
; 	n_points = [6,13,7]
; 	interval = (node(1:3)-node(0:2))/n_points

; 	if idx_target eq 0 then begin ;SVO
; 	alt_fin_target= [node(1)+findgen(n_points(1))*interval(1), $
; 				node(2)+findgen(n_points(2))*interval(2), $
; 				node(3)]	
; 	endif
; 	if idx_target eq 1 then begin ;LABS
; 	alt_fin_target= [node(0)+findgen(n_points(0))*interval(0), $
; 				node(1)+findgen(n_points(1))*interval(1), $
; 				node(2)+findgen(n_points(2))*interval(2), $
; 				node(3)]
; 	endif

;--modified version (2019/11/13);-----
alt_us = alt_us * 0.001
	n_interval0 = 10
	n_interval1 = 90
	interval0 = (Surf_pres_hPa-clars_pres_hPa_r)/(n_interval0)
	interval1 = (clars_pres_hPa_r-10)/(n_interval1)
	Pres_fin = [Surf_pres_hPa-findgen(n_interval0)*interval0, clars_pres_hPa_r-findgen(n_interval1)*interval1, $
			10, 5, 1, 0.5, 0.1]
	; n_interval0 = 100.
	; interval2 = (Surf_pres_hPa-10.)/(n_interval0)
	; Pres_fin = [Surf_pres_hPa-findgen(n_interval0)*interval0, $
	; 		10, 5, 1, 0.5, 0.1]
			

	node = [Surf_elevation_km, CLARS_alt_km_r, CLARS_alt_km_r+(CLARS_alt_km_r-Surf_elevation_km), 15, 60]
	n_points = [6,6,6,6]
	interval = (node(1:4)-node(0:3))/n_points

	if idx_target eq 0 then begin ;SVO
	alt_fin_target= [node(1)+findgen(n_points(1))*interval(1), $
				node(2)+findgen(n_points(2))*interval(2), $
				node(3)+findgen(n_points(3))*interval(3), $
				node(4)]	
	endif
	if idx_target eq 1 then begin ;LABS
	alt_fin_target= [node(0)+findgen(n_points(0))*interval(0), $
				node(1)+findgen(n_points(1))*interval(1), $
				node(2)+findgen(n_points(2))*interval(2), $
				node(3)+findgen(n_points(3))*interval(3), $
				node(4)]
	endif



n_level=n_elements(Pres_fin)




;---read another profile file: 
Initial_NSW_file = '/home/mchoi/OSSE_tool/Datasets/data_PTHG/data_originally_used/Initial_NSW_atmosphere_12aug13.dat'
readcol, Initial_NSW_file, pres_nswi, temp_nswi, h2o_nswi, co2_nswi, o2_nswi, ch4_nswi, co21_nswi, co22_nswi, /silent



H2O_fin=interpol(h2o_nswi,pres_nswi/100.,pres_fin) ; not spline
CO2_fin=interpol(co2_nswi,pres_nswi/100.,pres_fin) ; not spline
O2_fin =interpol(o2_nswi, pres_nswi/100.,pres_fin) ; not spline
CH4_fin=interpol(ch4_nswi,pres_nswi/100.,pres_fin) ; not spline
CO21_fin=interpol(co21_nswi,pres_nswi/100.,pres_fin) ; not spline
CO22_fin=interpol(co22_nswi,pres_nswi/100.,pres_fin) ; not spline


; ---"fin" setting is
; ---pres_fin, temp_fin, alt_fin
; ---H2O_fin, CO2_fin, CH4_fin, CO21_fin, CO22_fin
; ---to fin 3 directly!

Pres_fin3=Pres_fin
H2O_fin3=interpol(QV3/0.622, PL3/100, Pres_fin3)
Temp_fin2=interpol(T3, PL3/100., Pres_fin3)
; Temp_fin3 = Temp_fin2*!values.f_nan
Temp_fin3 = Temp_fin2

; ;--original ; 0-CLARS: linear based on sfc andn CLARS pressure; CLARS-upper: interpolation
; ; ;Temp_fin3(0:n_interval0-1)=Temp_fin2(0:n_interval0-1) - (Temp_fin2(0)-Surf_temp_K)
; Temp_fin3(0:n_interval0-1)=interpol([Surf_temp_K,CLARS_temp_K],[0,n_interval0-1],[findgen(n_interval0)])
; Temp_fin3(n_interval0:n_elements(temp_fin3)-1)=Temp_fin2(n_interval0:n_elements(temp_fin3)-1) - (Temp_fin2(n_interval0)-CLARS_temp_K)

;--modified one ; 0-top: interpolation (not using CLARS location pressure/temperature)
; Temp_fin3(0:n_interval0-1)=Temp_fin2(0:n_interval0-1) - (Temp_fin2(0)-Surf_temp_K)
; Temp_fin3(0:n_interval0-1)=interpol([Surf_temp_K,CLARS_temp_K],[0,n_interval0-1],[findgen(n_interval0)])
; Temp_fin3(n_interval0:n_elements(temp_fin3)-1)=Temp_fin2(n_interval0:n_elements(temp_fin3)-1) - (Temp_fin2(n_interval0)-CLARS_temp_K)


; Alt_fin3=interpol(H3/1000.,PL3/100.,Pres_fin3,/spline)


;--again from NSW initial dataset
CO2_fin3=interpol(co2_nswi,pres_nswi/100.,pres_fin3)
CO21_fin3=interpol(co21_nswi,pres_nswi/100.,pres_fin3)
CO22_fin3=interpol(co22_nswi,pres_nswi/100.,pres_fin3)
O2_fin3=interpol(O2_nswi,pres_nswi/100.,pres_fin3)
CH4_fin3=interpol(CH4_nswi,pres_nswi/100.,pres_fin3)




;--from GEMSTOOL_NSW_pthg_profiles.f90 fortran code
alt_fin = dblarr(n_level)*!values.d_nan
alt_fin(n_level-1)=Surf_elevation_km
alt_fin(n_level-1-n_interval0)=CLARS_alt_km_r ;--two points

;--reverse of temp_fin3 -> temp_fin4
;--reverse of pres_fin3 -> pres_fin4
temp_fin4 = double(reverse(temp_fin3))
pres_fin4 = double(reverse(pres_fin3))



ccon = - 9.81d0 * 28.9d0 / 8314.0d0 * 500.0d0
; for n = n_level-1, 1, -1 do begin
; avit = (1.0d0/temp_fin4(n-1))+(1.0d0/temp_fin4(n))
; alt_fin(n-1) = alt_fin(n) - alog(pres_fin4(n)/pres_fin4(n-1))/avit/ccon
; endfor


;----should be changed!!!

for n = n_level-1, n_level-1-8, -1 do begin
avit = (1.0d0/temp_fin4(n-1))+(1.0d0/temp_fin4(n))
alt_fin(n-1) = alt_fin(n) - alog(pres_fin4(n)/pres_fin4(n-1))/avit/ccon
endfor

for n = n_level-1-10, 1, -1 do begin
avit = (1.0d0/temp_fin4(n-1))+(1.0d0/temp_fin4(n))
alt_fin(n-1) = alt_fin(n) - alog(pres_fin4(n)/pres_fin4(n-1))/avit/ccon
endfor

;----original

; openw, 9, rtm_pwd+'/../GEMSTOOL_physicsdata/PTH_PROFILES/'+file_title
; for i = n_elements(pres_fin)-1, 0L, -1 do begin
; 	; printf, 9, pres_fin(i)*100., temp_fin(i), H2O_fin(i), CO2_fin(i), O2_fin(i), CH4_fin(i), CO21_fin(i), CO22_fin(i), $
; 	printf, 9, pres_fin3(i)*100., temp_fin3(i), H2O_fin3(i), O2_fin3(i), $
; 	          f='((1(f13.2),x),(3(e14.4),x))'
; endfor
; close, 9

; openw, 9, rtm_pwd+'/../GEMSTOOL_physicsdata/GAS_PROFILES/'+file_title
; for i = n_elements(pres_fin)-1, 0L, -1 do begin
; 	; printf, 9, pres_fin(i)*100., temp_fin(i), H2O_fin(i), CO2_fin(i), O2_fin(i), CH4_fin(i), CO21_fin(i), CO22_fin(i), $
; 	printf, 9, pres_fin3(i)*100., temp_fin3(i), H2O_fin3(i), O2_fin3(i), $
; 	          f='((1(f13.2),x),(3(e14.4),x))'
; endfor
; close, 9

; openw, 9, rtm_pwd+'/NSW_heightgrid_1.dat'
; for i = 0L, n_level-1 do begin
; 	printf, 9, i, alt_fin(i), $
; 	          f='((1(i2.2),x),(1(f18.15),x))'
; endfor
; close, 9





;-----Then, setted products
;-----pres_fin3*100, temp_fin3, H2O_fin3, O2_fin3
;-----Pres_fin3->reverse to Pres_fin4, and alt_fin

;-----strategy
;-----interpolated to alt_fin_target ->Pres_fin5
;-----reverse Pres_fin5 -> Pres_fin6, and interpol 
;Pres_fin5 = interpol(Pres_fin3,alt_fin,alt_fin_target,/spline)

log_Pres_fin5 = interpol(alog(Pres_Fin3),reverse(alt_fin),alt_fin_target)
Pres_fin5 = exp(log_Pres_fin5)

Temp_fin5 = interpol(temp_fin3,pres_fin3,pres_fin5) + tshift ;--added 
H2O_fin5 = interpol(H2O_fin3,pres_fin3,pres_fin5)
O2_fin5  = interpol(O2_fin3,pres_fin3,pres_fin5)



;---revised
; ;--save Pressure-Tempereture-Gas VMR
; file_title0 = 'atmos_mchoi'+string() ;--same with atmos_mchoi_t08
; file_title = file_title0+'.dat'

file_title = 'atmos_mchoi'+'_ii'+string(ii,f='(i4.4)')+'.dat'
; openw, 9, rtm_pwd+'/../GEMSTOOL_physicsdata/PTH_PROFILES/'+file_title
; openw, 9, rtm_pwd+'/../../GEMSTOOL_physicsdata/PTH_PROFILES/'+file_title
openw, 9, workspace_pwd+'/RTM/GEMSTOOL_physicsdata/PTH_PROFILES/'+file_title
for i = n_elements(alt_fin_target)-1, 0L, -1 do begin
	printf, 9, pres_fin5(i)*100., temp_fin5(i), H2O_fin5(i), O2_fin5(i), $
	          f='((1(f13.2),x),(3(e14.4),x))'
endfor
close, 9
; openw, 9, rtm_pwd+'/../GEMSTOOL_physicsdata/GAS_PROFILES/'+file_title
; openw, 9, rtm_pwd+'/../../GEMSTOOL_physicsdata/GAS_PROFILES/'+file_title
openw, 9, workspace_pwd+'/RTM/GEMSTOOL_physicsdata/GAS_PROFILES/'+file_title
for i = n_elements(alt_fin_target)-1, 0L, -1 do begin
	printf, 9, pres_fin5(i)*100., temp_fin5(i), H2O_fin5(i), O2_fin5(i), $
	          f='((1(f13.2),x),(3(e14.4),x))'
endfor
close, 9

file_title = 'NSW_heightgrid_1'+'_ii'+string(ii,f='(i4.4)')+'.dat'
openw, 9, workspace_pwd+'/RTM/GEMSTOOL_physicsdata/PTH_PROFILES/'+file_title;'/NSW_heightgrid_1.dat'
for i = n_elements(alt_fin_target)-1, 0L, -1 do begin
	printf, 9, n_elements(alt_fin_target)-1-i, alt_fin_target(i), $
	          f='((1(i2.2),x),(1(f18.15),x))'
endfor
close, 9






; spawn, 'cp ./'+file_title + rtm_pwd +' ../GEMSTOOL_physicsdata/PTH_PROFILES/'
; spawn, 'cp ./'+file_title + rtm_pwd +' ../GEMSTOOL_physicsdata/GAS_PROFILES/'

; ;---calculate LUT dataset
; cd, rtm_pwd+'/../GEMSTOOL_physicsdata/LUT_Xsections_data/'
; cd, rtm_pwd+'/../../GEMSTOOL_physicsdata/LUT_Xsections_data/'
; cd, workspace_pwd+'/RTM/GEMSTOOL_physicsdata/LUT_Xsections_data/'
; spawn, 'idl -e "make_xsec_lut_using_absco_v5_mchoi_o2ab" ' 


;---original
; cd, workspace_pwd+'/RTM/GEMSTOOL_physicsdata/LUT_Xsections_data/'
; spawn, 'idl -e "make_xsec_lut_using_absco_v5_mchoi_o2ab" ' 
; cd, workspace_pwd

;---updated 2019-11-23 ;--this is calculated in the runtest7_1_svo_lab_pthg_only.pro file
; make_xsec_lut_using_absco_v5_mchoi_o2ab_new_multiple, i_loop_start, input_loop_begin


; switch_plot_profile = 0

; if switch_plot_profile eq 1 then begin

; cd, workspace_pwd+'/RTM/GEMSTOOL_physicsdata/PTH_PROFILES';'/home/mchoi/OSSE_tool/for_PTHG_setting'
; ; spawn, 'cd /home/mchoi/OSSE_tool/for_PTHG_setting'

; ; ;--save layer-Altitude corresponding with layer-Pressure
; ; openw, 9, './NSW_heightgrid_'+file_title0+'.dat'
; ; for i = n_elements(pres_fin)-1, 0, -1 do begin
; ; 	printf, 9, n_elements(pres_fin)-1-i, pres_fin3(i), alt_fin3(i), $
; ; 	          f='(1(i12,3x), 1(f10.4,3x), 1(f7.4))'
; ; endfor
; ; close, 9

; ; spawn, 'cp ./NSW_heightgrid_'+file_title0+'.dat ../GEMSTOOL_mchoi_t12_00/GEMSTOOL_NSW_wrapper'

; ; spawn, 'cp ./NSW_heightgrid_'+file_title0+'.dat ../GEMSTOOL_mchoi_t08_00/GEMSTOOL_physicsdata/PTH_PROFILES/'

; ; --plot each data
; ; --Y-axis=alt_fin; X-axis= Temp_fin, Pres_fin, H2O_fin, O3_fin, NO2_fin, HCHO_fin, SO2_fin, CO_fin
; For i = 0L, 3 do begin

; if i eq 0 then begin
; file_dum = 'Profile_00_temp.eps'
; title_dum = 'Temperature (K)'
; xvar=Temp_fin5
; endif
; if i eq 1 then begin
; file_dum = 'Profile_01_pres.eps'
; title_dum = 'Pressure (hPa)'
; xvar=Pres_fin5
; endif
; if i eq 2 then begin
; file_dum = 'Profile_02_H2O.eps'
; title_dum = 'H2O (VMR)'
; xvar=H2O_fin5
; endif
; ; if i eq 3 then begin
; ; file_dum = '03_CO2.eps'
; ; title_dum = 'CO2 (VMR)'
; ; xvar=CO2_fin3
; ; endif
; if i eq 3 then begin
; file_dum = 'Profile_04_O2.eps'
; title_dum = 'O2 (VMR)'
; xvar=O2_fin5
; endif
; ; if i eq 5 then begin
; ; file_dum = '05_CH4.eps'
; ; title_dum = 'CH4 (VMR)'
; ; xvar=CH4_fin3
; ; endif
; ; if i eq 6 then begin
; ; file_dum = '06_CO2_1.eps'
; ; title_dum = 'CO2_1 (VMR)'
; ; xvar=CO21_fin3
; ; endif
; ; if i eq 7 then begin
; ; file_dum = '07_CO2_2.eps'
; ; title_dum = 'CO2_2 (VMR)'
; ; xvar=CO22_fin3
; ; endif

; yvar=(alt_fin_target) ;alt_fin

; xrange_dum=[min(xvar,/nan),max(xvar,/nan)]
; yrange_dum=[min(yvar,/nan),max(yvar,/nan)]

; ; usersym for circle filled
; A = FINDGEN(17) * (!PI*2/16.)
; USERSYM, COS(A) , SIN(A) , THICK=0.5, /FILL

; loadct,39,/silent
; !p.font=0
; set_plot,'ps'

; device, /encapsulated, /color, bits_per_pixel = 8, filename = file_dum, $
;         /helvetica, /bold, font_size = 15, xs = 30, ys = 30

; plot, xvar, yvar, psym = -8, xrange=xrange_dum, yrange=yrange_dum, $
;       background=0, color=0, $
; 	  xtitle=title_dum, ytitle='Altitude (km)', charsize=1,charthick=1.5, $
; 	  xstyle=0, ystyle=0, xmargin=[8,4], ymargin=[4,4], /nodata

; oplot, xvar, yvar, color=50, psym = -8, symsize=2, thick=5

; device,/close


; endfor

; endif;if switch_plot_profile eq 1 then begin

;--save surface pressure for state vector txt file ;--case of i_step = 0
if i_step eq 0 then begin
; result_pwd = Result_pwd0   ;workspace_pwd+'/Results'
state_vectors_pwd = result_pwd+ '/state_vectors_step'
str_step = 'step'+string(i_step,f='(i2.2)')
openw, 1, state_vectors_pwd+'/state_vectors_surfpress_'+str_step+'.txt'
printf, 1, Surf_pres_hPa;1008.327410d0 ;---initial value
close, 1
endif


; cd, workspace_pwd



end