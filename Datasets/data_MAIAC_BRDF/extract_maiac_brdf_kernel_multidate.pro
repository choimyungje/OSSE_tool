Pro extract_maiac_brdf_kernel_multidate

lat_target = 34.1101d
lon_target = -117.9695d

; ;---including all LABS site -- but this range includes ocean pixels a little.
; lon_range = [-118.470, -117.390]
; lat_range = [  33.678,   34.170]

;---reduced area
lon_range = [-118.330, -117.390]
lat_range = [  33.794,   34.170]


yyyy = 2019


; target_jd = julday(05,31,yyyy)
; doy = target_jd - julday(01,00,yyyy)
; doy_set = doy + indgen(8)

target_month_jd0 = julday(05,01,yyyy)
doy_set = target_month_jd0 - julday(01,00,yyyy) +indgen(8 + 31 - 1)
str_yyyy_doy_set = string(yyyy,f='(i4.4)')+string(doy_set,f='(i3.3)')

L2_pwd = '/home/mchoi/HIMAP/data_MAIAC_BRDF/data_hdf'
files = file_search(l2_pwd + '/MCD19A3.A'+str_yyyy_doy_set+'*.hdf', count=n_files)


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

; ;--finding within lon-lat box
idx = where(lon ge lon_range(0) and lon le lon_range(1) and lat ge lat_range(0) and lat le lat_range(1), n_idx)
idx_x = idx mod sz(0)
idx_y = idx / sz(0)

;-----------------------------

K_iso_all = dblarr(n_files, n_idx, n_wl_modis)*!values.d_nan
K_vol_all = dblarr(n_files, n_idx, n_wl_modis)*!values.d_nan
K_geo_all = dblarr(n_files, n_idx, n_wl_modis)*!values.d_nan

Albedo_all = dblarr(n_files, n_idx, n_wl_modis)*!values.d_nan

for i_file = 0L, n_files-1 do begin

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
  K_iso_all(i_file,*,i_wl) = K_iso_select(idx)

  K_vol_select = K_vol(*,*,sort_idx(i_wl))
  K_vol_all(i_file,*,i_wl) = K_vol_select(idx)

  K_geo_select = K_geo(*,*,sort_idx(i_wl))
  K_geo_all(i_file,*,i_wl) = K_geo_select(idx)

  Albedo_select = Albedo(*,*,sort_idx(i_wl))
  Albedo_all(i_file,*,i_wl) = Albedo_select(idx)

  endfor ;for i_wl = 0L, n_wl_modis-1 do begin

  


endfor ;for i_file = 0L, n_files-1 do begin

;---nan_idx in specific cases





;---mean and stddev for K_iso_all, K_vol_all, K_geo_all for modis channel first.
K_iso_LABS_mean = fltarr(n_wl_modis) * !values.f_nan
K_vol_LABS_mean = fltarr(n_wl_modis) * !values.f_nan
K_geo_LABS_mean = fltarr(n_wl_modis) * !values.f_nan

K_iso_LABS_stddev = fltarr(n_wl_modis) * !values.f_nan
K_vol_LABS_stddev = fltarr(n_wl_modis) * !values.f_nan
K_geo_LABS_stddev = fltarr(n_wl_modis) * !values.f_nan

Albedo_LABS_mean = fltarr(n_wl_modis) * !values.f_nan
Albedo_LABS_stddev = fltarr(n_wl_modis) * !values.f_nan

for i_wl = 0L, n_wl_modis -1 do begin

K_iso_LABS_mean(i_wl)   =   mean(K_iso_all(*,*,i_wl),/nan)
K_iso_LABS_stddev(i_wl) = stddev(K_iso_all(*,*,i_wl),/nan)

K_vol_LABS_mean(i_wl)   =   mean(K_vol_all(*,*,i_wl),/nan)
K_vol_LABS_stddev(i_wl) = stddev(K_vol_all(*,*,i_wl),/nan)

K_geo_LABS_mean(i_wl)   =   mean(K_geo_all(*,*,i_wl),/nan)
K_geo_LABS_stddev(i_wl) = stddev(K_geo_all(*,*,i_wl),/nan)

Albedo_LABS_mean(i_wl)   =   mean(Albedo_all(*,*,i_wl),/nan)
Albedo_LABS_stddev(i_wl) = stddev(Albedo_all(*,*,i_wl),/nan)

endfor


openw, 1, './save_MODIS_BRDF_Kernel_201305.txt'
printf, 1, 'Wavelength Wavenumber K_iso K_vol K_geo Albedo stddev_iso stddev_vol stddev_geo stddev_albedo'
for i=0L, n_wl_modis-1 do begin
printf, 1, wl_modis_sort(i), K_iso_LABS_mean(i), K_vol_LABS_mean(i), K_geo_LABS_mean(i), Albedo_LABS_mean(i), $
                             K_iso_LABS_stddev(i), K_vol_LABS_stddev(i), K_geo_LABS_stddev(i), Albedo_LABS_stddev(i), format='(9f)'
endfor
close, 1

save, wl_modis_sort, K_iso_LABS_mean, K_vol_LABS_mean, K_geo_LABS_mean, Albedo_LABS_mean, $
      K_iso_LABS_stddev, K_vol_LABS_stddev, K_geo_LABS_stddev, Albedo_LABS_stddev, filename = 'save_MAIAC_K_201305.xdr'


;---interpolated to the O2ABD channels
restore, '/home/mchoi/HIMAP/data_sfcref_for_tropomi/save_LER_GOME2.xdr'
;save, WL_GOME2, LER_GOME2_all, filename='save_LER_GOME2.xdr'

restore, '/home/mchoi/HIMAP/data_sfcref_for_tropomi/save_LER_SCIAMACHY.xdr'
;save, WL_SCIAMACHY, LER_SCIAMACHY_all, filename='save_LER_SCIAMACHY.xdr'


wn_start = [12900.00d, 14300.00d, 7700.00d]
wn_end   = [13250.00d, 14600.00d, 8000.00d]
wl_fin = 10.^7/ ((wn_start+wn_end)/2.0) ;764.818   692.042  1273.885
n_wl_Fin = n_elements(wl_Fin)
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

ratio_o2b_avg =   mean(ratio_o2b,/nan)
ratio_o2b_std = stddev(ratio_o2b,/nan)
ratio_o2a_avg =   mean(ratio_o2a,/nan)
ratio_o2a_std = stddev(ratio_o2a,/nan)
ratio_o2d_avg =   mean(ratio_o2d,/nan)
ratio_o2d_std = stddev(ratio_o2d,/nan)


;---calculating K_iso, K_vol, K_geo, NDVI, K_pol for O2A, O2B, and O2D bands


K_iso_fin = fltarr(n_wl_Fin)*!values.f_nan
K_vol_fin = fltarr(n_Wl_fin)*!values.f_nan
K_geo_fin = fltarr(n_wl_fin)*!values.f_nan
; K_pol_fin = fltarr(n_wl_fin)*!values.f_nan

K_iso_fin_stddev = fltarr(n_wl_Fin)*!values.f_nan
K_vol_fin_stddev = fltarr(n_Wl_fin)*!values.f_nan
K_geo_fin_stddev = fltarr(n_wl_fin)*!values.f_nan
; K_pol_fin_stddev = fltarr(n_wl_fin)*!values.f_nan

;   412.000      465.500      553.500      644.900      855.600      1241.90      1629.00      2113.10
K_iso_fin(0)= K_iso_LABS_mean(3) * (1. - ratio_o2a_avg) + K_iso_LABS_mean(4) * ratio_o2a_avg
K_vol_fin(0)= K_vol_LABS_mean(3) * (1. - ratio_o2a_avg) + K_vol_LABS_mean(4) * ratio_o2a_avg
K_geo_fin(0)= K_geo_LABS_mean(3) * (1. - ratio_o2a_avg) + K_geo_LABS_mean(4) * ratio_o2a_avg
; K_pol_fin(0)= K_pol_LABS_mean(3) * (1. - ratio_o2a_avg) + K_pol_LABS_mean(4) * ratio_o2a_avg

K_iso_fin(1)= K_iso_LABS_mean(3) * (1. - ratio_o2b_avg) + K_iso_LABS_mean(4) * ratio_o2b_avg
K_vol_fin(1)= K_vol_LABS_mean(3) * (1. - ratio_o2b_avg) + K_vol_LABS_mean(4) * ratio_o2b_avg
K_geo_fin(1)= K_geo_LABS_mean(3) * (1. - ratio_o2b_avg) + K_geo_LABS_mean(4) * ratio_o2b_avg
; K_pol_fin(1)= K_pol_LABS_mean(3) * (1. - ratio_o2b_avg) + K_pol_LABS_mean(4) * ratio_o2b_avg

K_iso_fin(2)= K_iso_LABS_mean(3) * (1. - ratio_o2d_avg) + K_iso_LABS_mean(4) * ratio_o2d_avg
K_vol_fin(2)= K_vol_LABS_mean(3) * (1. - ratio_o2d_avg) + K_vol_LABS_mean(4) * ratio_o2d_avg
K_geo_fin(2)= K_geo_LABS_mean(3) * (1. - ratio_o2d_avg) + K_geo_LABS_mean(4) * ratio_o2d_avg
; K_pol_fin(2)= K_pol_LABS_mean(3) * (1. - ratio_o2d_avg) + K_pol_LABS_mean(4) * ratio_o2d_avg

;---uncertainty together
K_iso_fin_stddev(0)= K_iso_LABS_stddev(3) * (1. - ratio_o2a_avg) + K_iso_LABS_stddev(4) * ratio_o2a_avg
K_vol_fin_stddev(0)= K_vol_LABS_stddev(3) * (1. - ratio_o2a_avg) + K_vol_LABS_stddev(4) * ratio_o2a_avg
K_geo_fin_stddev(0)= K_geo_LABS_stddev(3) * (1. - ratio_o2a_avg) + K_geo_LABS_stddev(4) * ratio_o2a_avg
; K_pol_fin_stddev(0)= K_pol_LABS_stddev(3) * (1. - ratio_o2a_avg) + K_pol_LABS_stddev(4) * ratio_o2a_avg

K_iso_fin_stddev(1)= K_iso_LABS_stddev(3) * (1. - ratio_o2b_avg) + K_iso_LABS_stddev(4) * ratio_o2b_avg
K_vol_fin_stddev(1)= K_vol_LABS_stddev(3) * (1. - ratio_o2b_avg) + K_vol_LABS_stddev(4) * ratio_o2b_avg
K_geo_fin_stddev(1)= K_geo_LABS_stddev(3) * (1. - ratio_o2b_avg) + K_geo_LABS_stddev(4) * ratio_o2b_avg
; K_pol_fin_stddev(1)= K_pol_LABS_stddev(3) * (1. - ratio_o2b_avg) + K_pol_LABS_stddev(4) * ratio_o2b_avg

K_iso_fin_stddev(2)= K_iso_LABS_stddev(3) * (1. - ratio_o2d_avg) + K_iso_LABS_stddev(4) * ratio_o2d_avg
K_vol_fin_stddev(2)= K_vol_LABS_stddev(3) * (1. - ratio_o2d_avg) + K_vol_LABS_stddev(4) * ratio_o2d_avg
K_geo_fin_stddev(2)= K_geo_LABS_stddev(3) * (1. - ratio_o2d_avg) + K_geo_LABS_stddev(4) * ratio_o2d_avg
; K_pol_fin_stddev(2)= K_pol_LABS_stddev(3) * (1. - ratio_o2d_avg) + K_pol_LABS_stddev(4) * ratio_o2d_avg

; NDVI_fin = NDVI_LABS_mean
; NDVI_fin_stddev = NDVI_LABS_stddev

;---write results as txt file
openw, 1, './MAIAC_BRDF_Kernels_for_O2_ABD.txt'
for i=0L, 2 do begin
printf, 1, wl_fin(i), K_iso_fin(i),K_vol_fin(i),K_geo_fin(i), $ 
	     K_iso_fin_stddev(i),K_vol_fin_stddev(i),K_geo_fin_stddev(i), f='(11(f10.5,1x))'
endfor
close, 1



stop


;-----figure
cgplot, 1, 1, psym = 1, charsize = 2, symsize=2, xtitle = 'wavelength(nm)', ytitle = 'K_iso or LER', /nodata, $
	  xrange = [300., 2200.], yrange = [0.0, 0.4]

cgplot, wl_modis_sort, K_iso_LABS_mean, psym = 2, color='black', /overplot, symsize=2
cgplot, wl_modis_sort, Albedo_LABS_mean, psym = 2, color='purple', /overplot, symsize=2
; cgplot, wl_POLDER, K_iso_GRASP_LABS_mean, psym = 2, color='red', /overplot, symsize=2
; cgplot, WL_GOME2, LER_GOME2_all(*,5-1), psym = -4, color='green', /overplot, symsize=2
cgplot, WL_SCIAMACHY, LER_SCIAMACHY_all(*,5-1), psym = -4, color='blue', /overplot, symsize=2

cgtext, 1000, 0.2, 'MODIS MAIAC', color='black', charsize=2
cgtext, 1000, 0.1, 'POLDER GRASP', color='red', charsize = 2

;-----figure

cgplot, 1, 1, psym = 1, charsize = 2, symsize=2, xtitle = 'wavelength(nm)', ytitle = 'LER', /nodata, $
	  xrange = [600., 1300.], yrange = [0.0, 0.32], charthick = 2
cgloadct,39,/silent

for i=0L, n_elements(wl_modis_sort)-1 do cgplot, [wl_modis_sort(i),wl_modis_sort(i)], [-1.0, 2.0], color='Blue', /overplot, thick=5
for i=0L, n_elements(wl_fin)-1 do cgplot, [wl_fin(i),wl_fin(i)], [-1.0, 2.0], color='Orange', /overplot, thick=5
for i=0L, 11 do cgplot, WL_SCIAMACHY, LER_SCIAMACHY_all(*,i), psym = -1, color='grey', /overplot, symsize=2

for i=0L, 11 do cgplot, wl_test, LER_SCIAMACHY_all_test(*,i), psym = 16, color='black', /overplot, symsize=2

cgtext, 900, 0.10, 'MAIAC bands center', color='Blue', charsize=2, charthick = 2
cgtext, 900, 0.08, 'Target O!d2!n absorption bands center', color='Orange', charsize = 2, charthick = 2
cgtext, 900, 0.06, 'SCIAMACHY climatological minimum LER', color='Grey', charsize=2, charthick = 2

cgplot, wl_modis_sort, K_vol_LABS_mean, psym = 1
cgplot, wl_POLDER, K_vol_GRASP_LABS_mean, psym = 2, color='red', /overplot

cgplot, wl_modis_sort, K_geo_LABS_mean, psym = 1
cgplot, wl_POLDER, K_geo_GRASP_LABS_mean, psym = 2, color='red', /overplot
stop


stop
end
