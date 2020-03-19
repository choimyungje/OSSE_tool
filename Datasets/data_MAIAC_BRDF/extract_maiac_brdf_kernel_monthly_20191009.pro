Pro extract_maiac_brdf_kernel_monthly_20191009

wn_start = [12900.00d, 14300.00d, 7700.00d]
wn_end   = [13250.00d, 14600.00d, 8000.00d]
wl_fin = 10.^7/ ((wn_start+wn_end)/2.0) ;764.818   692.042  1273.885
n_wl_Fin = n_elements(wl_Fin)

;----differently for each site
restore, '/home/mchoi/HIMAP/data_CLARS_FTS/save_target_info/save_target_info.xdr'
; save, loc_set, lat_set, lon_set, alt_set, VZA_set, VAA_set, $
;       filename = save_pwd +'/save_target_info.xdr'


work_pwd = '/home/mchoi/HIMAP/data_MAIAC_BRDF'
save_pwd = work_pwd+'/save_monthly_site_results'
file_mkdir, save_pwd

n_set = n_elements(loc_set)
n_mm = 12
K_iso_avg = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
K_iso_std = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
K_vol_avg = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
K_vol_std = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
K_geo_avg = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
K_geo_std = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
Albedo_avg = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
Albedo_std = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
;_alog
K_iso_avg_alog = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
K_iso_std_alog = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
K_vol_avg_alog = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
K_vol_std_alog = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
K_geo_avg_alog = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
K_geo_std_alog = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
Albedo_avg_alog = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
Albedo_std_alog = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan

;---how to set lon-lat?

nx=1200L ;for 1 km 
ny=1200L ;for 1 km 
wx=43200L ;for 1km  
wy=21600L ;for 1km 

h_idx = dblarr(nx,ny)*!values.d_nan & h_idx(*)= 08;h_idx_read ;--08 for LA basin
v_idx = dblarr(nx,ny)*!values.d_nan & v_idx(*)= 05;v_idx_read ;--05 for LA basin
tx = (h_idx(*,*)) 
ty = (v_idx(*,*))
lat=replicate(1.d,nx) # ((90.d)-double(ty)/(18.d)*(180.d)-(dindgen(ny)+.5d)/double(wy)*(180.d))
lon=((dindgen(nx)+tx*nx)/double(wx)*(360.d)-(180.d))#((1.d)/(sin((dindgen(ny)+ty*ny)*!dpi/double(wy))>1E-8))

;--finding location idx first

wl_modis = [644.9, 855.6 ,465.5, 553.5, 1241.9 ,1629.0 ,2113.1, 412.]
n_wl_modis = n_elements(wl_modis)
sort_idx = sort(wl_modis)

wl_modis_sort = wl_modis(sort_idx)

sz = size(lon,/dimension)
;lon_target, lat_target

; ;---reduced area
; lon_range = [min(lon_set)-0.1,max(lon_set)+0.1] ;[-118.330, -117.390]
; lat_range = [min(lat_set)-0.1,max(lat_set)+0.1] ;[  33.794,   34.170]
;------33.794, -118.330     34.170, -117.390
; ; ;--finding within lon-lat box
; idx = where(lon ge lon_range(0) and lon le lon_range(1) and lat ge lat_range(0) and lat le lat_range(1), n_idx)
; idx_x = idx mod sz(0)
; idx_y = idx / sz(0)



;--setup for wl interpolation using SCIAMACHY 
;---interpolated to the O2ABD channels
restore, '/home/mchoi/HIMAP/data_sfcref_for_tropomi/save_LER_SCIAMACHY.xdr'
;save, WL_SCIAMACHY, LER_SCIAMACHY_all, filename='save_LER_SCIAMACHY.xdr'

;;--wl_modis_sort
;;--412.000      465.500      553.500      644.900      855.600      1241.90      1629.00      2113.10
wl_test = [wl_modis_sort, wl_fin]
n_wl_test = n_elements(wl_test)
LER_SCIAMACHY_all_test = fltarr(n_wl_test, 12) * !values.f_nan
for i=0L, 11 do LER_SCIAMACHY_all_test(*,i) = interpol(LER_SCIAMACHY_all(*,i),WL_SCIAMACHY,wl_test)

ratio_o2b = fltarr(12)*!values.f_nan
ratio_o2a = fltarr(12)*!values.f_nan
ratio_o2d = fltarr(12)*!values.f_nan

for i=0L, 11 do ratio_o2b(i) = (LER_SCIAMACHY_all_test(9,i)-LER_SCIAMACHY_all_test(3,i))/(LER_SCIAMACHY_all_test(4,i)-LER_SCIAMACHY_all_test(3,i))
for i=0L, 11 do ratio_o2a(i) = (LER_SCIAMACHY_all_test(8,i)-LER_SCIAMACHY_all_test(3,i))/(LER_SCIAMACHY_all_test(4,i)-LER_SCIAMACHY_all_test(3,i))
for i=0L, 11 do ratio_o2d(i) = (LER_SCIAMACHY_all_test(10,i)-LER_SCIAMACHY_all_test(3,i))/(LER_SCIAMACHY_all_test(4,i)-LER_SCIAMACHY_all_test(3,i))




;--read each month MAIAC data and calculate!


for mm = 0, n_mm-1 do begin

str_yyyydoyset = strarr(1)
for yyyy = 2000, 2019 do begin
  doy_set = [julday(mm+1,01,yyyy) - julday(01,00,yyyy): julday(mm+2,00,yyyy) - julday(01,00,yyyy)] + 4. ;--8-day composite -> 4 day shift for center
  negative_idx = where(doy_set gt julday(12,31,yyyy) - julday(01,00,yyyy), n_negative_idx)
  str_yyyydoy = string(yyyy,f='(i4.4)')+string(doy_set,f='(i3.3)')
  if n_negative_idx ge 1 then begin
  doy_set(negative_idx) = doy_set(negative_idx) - (julday(12,31,yyyy)-julday(1,0,yyyy))
  str_yyyydoy(negative_idx) = string(yyyy+1,f='(i4.4)')+string(doy_set(negative_idx),f='(i3.3)')
  endif

  str_yyyydoyset = [str_yyyydoyset, str_yyyydoy]
  
endfor
n_str_yyyydoyset = n_elements(str_yyyydoyset)
str_yyyydoyset = str_yyyydoyset(1:n_str_yyyydoyset-1)

L2_pwd = work_pwd+'/data_hdf_dn_20191008'
files = file_search(l2_pwd + '/MCD19A3.A'+str_yyyydoyset+'*.hdf', count=n_files)


;-----------------------------

K_iso_modis = dblarr(n_files, n_set, n_wl_modis)*!values.d_nan
K_vol_modis = dblarr(n_files, n_set, n_wl_modis)*!values.d_nan
K_geo_modis = dblarr(n_files, n_set, n_wl_modis)*!values.d_nan
Albedo_modis = dblarr(n_files, n_set, n_wl_modis)*!values.d_nan

K_iso_O2ABD = dblarr(n_files, n_set, n_wl_fin)*!values.d_nan
K_vol_O2ABD = dblarr(n_files, n_set, n_wl_fin)*!values.d_nan
K_geo_O2ABD = dblarr(n_files, n_set, n_wl_fin)*!values.d_nan
Albedo_O2ABD = dblarr(n_files, n_set, n_wl_fin)*!values.d_nan


for i_file = 0L, n_files-1 do begin

  print, 'Month: ', mm+1, '    files: ', i_file, n_files 

  file=files(i_file)

  temp_idx1 = strsplit(file,'MCD19A3.A',/extract, /regex)
  temp_idx2 = strsplit(temp_idx1(1), '.', /extract)
  temp_idx3 = strsplit(temp_idx2(1), 'hv', /extract)
  h_idx_read = temp_idx3(0)
  v_idx_read = temp_idx3(1)

  hdfid=hdf_sd_start(file)

  ;---Kiso -> K_iso
  valid_range = hdf_sd_attinfo(hdfid, 'Kiso', 'valid_range')
  scale_factor = hdf_sd_attinfo(hdfid, 'Kiso', 'scale_factor')
  add_offset = hdf_sd_attinfo(hdfid, 'Kiso', 'add_offset')
  nan_value = hdf_sd_attinfo(hdfid, 'Kiso', '_FillValue')
  hdf_sd_varread, hdfid, 'Kiso', var_read
  nan_idx=where(var_read eq nan_value.data(0) or var_read lt valid_range.data(0) or var_read gt valid_range.data(1) , n_nan_idx)
  if n_nan_idx ge 1 then var_read(nan_idx)=!values.d_nan
  var_read = (var_read - add_offset.data(0))*scale_factor.data(0)
  K_iso = var_read

  ;---Kvol -> K_vol
  valid_range = hdf_sd_attinfo(hdfid, 'Kvol', 'valid_range')
  scale_factor = hdf_sd_attinfo(hdfid, 'Kvol', 'scale_factor')
  add_offset = hdf_sd_attinfo(hdfid, 'Kvol', 'add_offset')
  nan_value = hdf_sd_attinfo(hdfid, 'Kvol', '_FillValue')
  hdf_sd_varread, hdfid, 'Kvol', var_read
  nan_idx=where(var_read eq nan_value.data(0) or var_read lt valid_range.data(0) or var_read gt valid_range.data(1) , n_nan_idx)
  if n_nan_idx ge 1 then var_read(nan_idx)=!values.d_nan
  var_read = (var_read - add_offset.data(0))*scale_factor.data(0)
  K_vol = var_read

  ;---Kgeo -> K_geo
  valid_range = hdf_sd_attinfo(hdfid, 'Kgeo', 'valid_range')
  scale_factor = hdf_sd_attinfo(hdfid, 'Kgeo', 'scale_factor')
  add_offset = hdf_sd_attinfo(hdfid, 'Kgeo', 'add_offset')
  nan_value = hdf_sd_attinfo(hdfid, 'Kgeo', '_FillValue')
  hdf_sd_varread, hdfid, 'Kgeo', var_read
  nan_idx=where(var_read eq nan_value.data(0) or var_read lt valid_range.data(0) or var_read gt valid_range.data(1) , n_nan_idx)
  if n_nan_idx ge 1 then var_read(nan_idx)=!values.d_nan
  var_read = (var_read - add_offset.data(0))*scale_factor.data(0)
  K_geo = var_read

  ;---Sur_albedo -> Albedo
  valid_range = hdf_sd_attinfo(hdfid, 'Sur_albedo', 'valid_range')
  scale_factor = hdf_sd_attinfo(hdfid, 'Sur_albedo', 'scale_factor')
  add_offset = hdf_sd_attinfo(hdfid, 'Sur_albedo', 'add_offset')
  nan_value = hdf_sd_attinfo(hdfid, 'Sur_albedo', '_FillValue')
  hdf_sd_varread, hdfid, 'Sur_albedo', var_read
  nan_idx=where(var_read eq nan_value.data(0) or var_read lt valid_range.data(0) or var_read gt valid_range.data(1) , n_nan_idx)
  if n_nan_idx ge 1 then var_read(nan_idx)=!values.d_nan
  var_read = (var_read - add_offset.data(0))*scale_factor.data(0)
  Albedo = var_read

  ;nan_idx filtering 
  nan_idx = where(K_iso eq 0.0d or K_vol eq 0.0d or K_geo eq 0.0d or Albedo eq 0.0d, n_nan_idx)
  if n_nan_idx ge 1 then begin
    K_iso(nan_idx)=!values.d_nan
    K_vol(nan_idx)=!values.d_nan
    K_geo(nan_idx)=!values.d_nan
    Albedo(nan_idx)=!values.d_nan
  endif
  

  hdf_sd_end, hdfid

  ;-----------------------------------------------------------------

  for i_wl = 0L, n_wl_modis-1 do begin

  K_iso_select = K_iso(*,*,sort_idx(i_wl))
  ;--Nan-value remove for griddata
  usable_idx = where(finite(K_iso_select) eq 1, n_usable_idx)
  if n_usable_idx ge 1 then begin
  K_iso_modis(i_file,*,i_wl) = griddata(lon(usable_idx),lat(usable_idx),K_iso_select(usable_idx), xout=lon_set, yout=lat_set)
  endif

  K_vol_select = K_vol(*,*,sort_idx(i_wl))
  ;--Nan-value remove for griddata
  usable_idx = where(finite(K_vol_select) eq 1, n_usable_idx)
  if n_usable_idx ge 1 then begin
  K_vol_modis(i_file,*,i_wl) = griddata(lon(usable_idx),lat(usable_idx),K_vol_select(usable_idx), xout=lon_set, yout=lat_set)
  endif

  K_geo_select = K_geo(*,*,sort_idx(i_wl))
  ;--Nan-value remove for griddata
  usable_idx = where(finite(K_geo_select) eq 1, n_usable_idx)
  if n_usable_idx ge 1 then begin
  K_geo_modis(i_file,*,i_wl) = griddata(lon(usable_idx),lat(usable_idx),K_geo_select(usable_idx), xout=lon_set, yout=lat_set)
  endif

  Albedo_select = Albedo(*,*,sort_idx(i_wl))
  ;--Nan-value remove for griddata
  usable_idx = where(finite(Albedo_select) eq 1, n_usable_idx)
  if n_usable_idx ge 1 then begin
  Albedo_modis(i_file,*,i_wl) = griddata(lon(usable_idx),lat(usable_idx),Albedo_select(usable_idx), xout=lon_set, yout=lat_set)
  endif


  endfor ;for i_wl = 0L, n_wl_modis-1 do begin

  ;---intepolate modis channels to the O2ABD channels using each month SCIAMACHY measurement
  ;---from K_iso_modis -> K_iso_O2ABD and others also
  
  ;---interpolate to O2A band
  ;   412.000      465.500      553.500      644.900      855.600      1241.90      1629.00      2113.10
  K_iso_O2ABD(i_file,*,0)  = K_iso_modis(i_file,*,3) * (1. - ratio_o2a(mm))+ K_iso_modis(i_file,*,4) * ratio_o2a(mm)
  K_vol_O2ABD(i_file,*,0)  = K_vol_modis(i_file,*,3) * (1. - ratio_o2a(mm))+ K_vol_modis(i_file,*,4) * ratio_o2a(mm)
  K_geo_O2ABD(i_file,*,0)  = K_geo_modis(i_file,*,3) * (1. - ratio_o2a(mm))+ K_geo_modis(i_file,*,4) * ratio_o2a(mm)
  Albedo_O2ABD(i_file,*,0) = Albedo_modis(i_file,*,3) * (1. - ratio_o2a(mm))+ Albedo_modis(i_file,*,4) * ratio_o2a(mm)
  
  K_iso_O2ABD(i_file,*,1)  = K_iso_modis(i_file,*,3) * (1. - ratio_o2b(mm))+ K_iso_modis(i_file,*,4) * ratio_o2b(mm)
  K_vol_O2ABD(i_file,*,1)  = K_vol_modis(i_file,*,3) * (1. - ratio_o2b(mm))+ K_vol_modis(i_file,*,4) * ratio_o2b(mm)
  K_geo_O2ABD(i_file,*,1)  = K_geo_modis(i_file,*,3) * (1. - ratio_o2b(mm))+ K_geo_modis(i_file,*,4) * ratio_o2b(mm)
  Albedo_O2ABD(i_file,*,1) = Albedo_modis(i_file,*,3) * (1. - ratio_o2b(mm))+ Albedo_modis(i_file,*,4) * ratio_o2b(mm)
  
  K_iso_O2ABD(i_file,*,2)  = K_iso_modis(i_file,*,3) * (1. - ratio_o2d(mm))+ K_iso_modis(i_file,*,4) * ratio_o2d(mm)
  K_vol_O2ABD(i_file,*,2)  = K_vol_modis(i_file,*,3) * (1. - ratio_o2d(mm))+ K_vol_modis(i_file,*,4) * ratio_o2d(mm)
  K_geo_O2ABD(i_file,*,2)  = K_geo_modis(i_file,*,3) * (1. - ratio_o2d(mm))+ K_geo_modis(i_file,*,4) * ratio_o2d(mm)
  Albedo_O2ABD(i_file,*,2) = Albedo_modis(i_file,*,3) * (1. - ratio_o2d(mm))+ Albedo_modis(i_file,*,4) * ratio_o2d(mm)
  
endfor ;for i_file = 0L, n_files-1 do begin

;---calculate avg and std of monthly kernel info of each site

for i_ch=0L, n_wl_fin-1 do begin
for i_set=0L, n_set-1 do begin  
  K_iso_avg(i_ch,mm,i_set) =   mean(K_iso_O2ABD(*,i_set,i_ch),/nan)
  K_iso_std(i_ch,mm,i_set) = stddev(K_iso_O2ABD(*,i_set,i_ch),/nan)
  K_vol_avg(i_ch,mm,i_set) =   mean(K_vol_O2ABD(*,i_set,i_ch),/nan)
  K_vol_std(i_ch,mm,i_set) = stddev(K_vol_O2ABD(*,i_set,i_ch),/nan)
  K_geo_avg(i_ch,mm,i_set) =   mean(K_geo_O2ABD(*,i_set,i_ch),/nan)
  K_geo_std(i_ch,mm,i_set) = stddev(K_geo_O2ABD(*,i_set,i_ch),/nan)
  Albedo_avg(i_ch,mm,i_set) =   mean(Albedo_O2ABD(*,i_set,i_ch),/nan)
  Albedo_std(i_ch,mm,i_set) = stddev(Albedo_O2ABD(*,i_set,i_ch),/nan)

  K_iso_avg_alog(i_ch,mm,i_set) =   mean(alog(K_iso_O2ABD(*,i_set,i_ch)),/nan)
  K_iso_std_alog(i_ch,mm,i_set) = stddev(alog(K_iso_O2ABD(*,i_set,i_ch)),/nan)
  K_vol_avg_alog(i_ch,mm,i_set) =   mean(alog(K_vol_O2ABD(*,i_set,i_ch)),/nan)
  K_vol_std_alog(i_ch,mm,i_set) = stddev(alog(K_vol_O2ABD(*,i_set,i_ch)),/nan)
  K_geo_avg_alog(i_ch,mm,i_set) =   mean(alog(K_geo_O2ABD(*,i_set,i_ch)),/nan)
  K_geo_std_alog(i_ch,mm,i_set) = stddev(alog(K_geo_O2ABD(*,i_set,i_ch)),/nan)
  Albedo_avg_alog(i_ch,mm,i_set) =   mean(alog(Albedo_O2ABD(*,i_set,i_ch)),/nan)
  Albedo_std_alog(i_ch,mm,i_set) = stddev(alog(Albedo_O2ABD(*,i_set,i_ch)),/nan)


endfor 
endfor



endfor ;for mm = 0L, n_mm-1 do begin


save, K_iso_avg, K_iso_std, K_vol_avg, K_vol_std, K_geo_avg, K_geo_std, Albedo_avg, Albedo_std, $
      filename = save_pwd + '/save_MAIAC_BRDF_month_site.xdr'


save, K_iso_avg_alog, K_iso_std_alog, K_vol_avg_alog, K_vol_std_alog, $
      K_geo_avg_alog, K_geo_std_alog, Albedo_avg_alog, Albedo_std_alog, $
      filename = save_pwd + '/save_MAIAC_BRDF_month_site_alog.xdr'


stop
end
