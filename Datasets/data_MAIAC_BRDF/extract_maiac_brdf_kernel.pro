Pro extract_maiac_brdf_kernel

lat_target = 34.1101d
lon_target = -117.9695d

yyyy = 2019

target_jd = julday(05,29,yyyy)
doy = target_jd - julday(01,00,yyyy)
doy_set = doy-indgen(8)
str_yyyy_doy_set = string(yyyy,f='(i4.4)')+string(doy_set,f='(i3.3)')

L2_pwd = '/home/mchoi/HIMAP/data_MAIAC_BRDF/data_hdf'
files = file_search(l2_pwd + '/MCD19A3.A'+str_yyyy_doy_set+'*.hdf', count=n_files)

file=files(0)

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

;---how to set lon-lat?

nx=1200L ;for 1 km 
ny=1200L ;for 1 km 
wx=43200L ;for 1km  
wy=21600L ;for 1km 

h_idx = dblarr(nx,ny)*!values.d_nan & h_idx(*)= h_idx_read
v_idx = dblarr(nx,ny)*!values.d_nan & v_idx(*)= v_idx_read
tx = (h_idx(*,*))
ty = (v_idx(*,*))
lat=replicate(1.d,nx) # ((90.d)-double(ty)/(18.d)*(180.d)-(dindgen(ny)+.5d)/double(wy)*(180.d))
lon=((dindgen(nx)+tx*nx)/double(wx)*(360.d)-(180.d))#((1.d)/(sin((dindgen(ny)+ty*ny)*!dpi/double(wy))>1E-8))

hdf_sd_end, hdfid

;-----------------------------------------------------------------

wl_modis = [644.9, 855.6 ,465.5, 553.5, 1241.9 ,1629.0 ,2113.1, 412.]

sz = size(lon,/dimension)
;lon_target, lat_target

; ;--finding 0.1x0.1 degree range first

idx = where(abs(lat_target - lat) lt 0.3 and abs(lon_target-lon) lt 0.3, n_idx)
if n_idx ge 1 then begin
  distance=fltarr(n_idx)*!values.f_nan
  for i=0L, n_idx-1 do begin
      distance(i)=map_2points(lon_target,lat_target,lon(idx(i)),lat(idx(i)),/meters)
  endfor
idx_closest=where(distance eq min(distance,/nan), n_idx_closest)
idx2=idx(idx_closest(0))
idx2_x = idx2 mod sz(1)
idx2_y = idx2 / sz(1)

;---calculate within 10 km
idx_10km=where(distance le 10000. , n_idx_10km)
idx3=idx(idx_10km)
idx_10km_x = idx3 mod sz(1)
idx_10km_y = idx3 / sz(1)

endif




; lon(idx(idx_closest(0)))
; lat(idx(idx_closest(0)))


;--extracting K_iso, K_vol, K_geo

K_iso_target = K_iso(idx2_x,idx2_y,*)
K_vol_target = K_vol(idx2_x,idx2_y,*)
K_geo_target = K_geo(idx2_x,idx2_y,*)

Albedo_target = Albedo(idx2_x,idx2_y,*)


idx_sort= sort(wl_modis)
wl_modis_sort = wl_modis(idx_sort)
K_iso_target_sort = K_iso_target(idx_sort)
K_geo_target_sort = K_geo_target(idx_sort)
K_vol_target_sort = K_vol_target(idx_sort)
Albedo_target_sort = Albedo_target(idx_sort)

save, wl_modis_sort, K_iso_target_sort, K_geo_target_sort, K_vol_target_sort, filename = 'save_K_modis.xdr'


wl_CLARS =  [770.0, 690.0, 1270.0]; O2-B, O2-A, O2- 1270nm
wn_CLARS = 10.^7/wl_CLARS ;--14492.754       12987.013       7874.0156

K_iso_CLARS = interpol(K_iso_target_sort,wl_modis_sort,wl_CLARS)
K_geo_CLARS = interpol(K_geo_target_sort,wl_modis_sort,wl_CLARS)
K_vol_CLARS = interpol(K_vol_target_sort,wl_modis_sort,wl_CLARS)
Albedo_CLARS = interpol(Albedo_target_sort,wl_modis_sort,wl_CLARS)

;---average within 10 km 

K_iso_10km = fltarr(n_idx_10km, 8) * !values.f_nan
K_vol_10km = fltarr(n_idx_10km, 8) * !values.f_nan
K_geo_10km = fltarr(n_idx_10km, 8) * !values.f_nan
Albedo_10km = fltarr(n_idx_10km, 8) * !values.f_nan

for i=0L, n_idx_10km-1 do begin
K_iso_10km(i,*) = K_iso(idx_10km_x(i),idx_10km_y(i),*)
K_vol_10km(i,*) = K_vol(idx_10km_x(i),idx_10km_y(i),*)
K_geo_10km(i,*) = K_geo(idx_10km_x(i),idx_10km_y(i),*)
Albedo_10km(i,*) = Albedo(idx_10km_x(i),idx_10km_y(i),*)
endfor

K_iso_10km_mean = fltarr(8)*!values.f_nan
K_vol_10km_mean = fltarr(8)*!values.f_nan
K_geo_10km_mean = fltarr(8)*!values.f_nan
Albedo_10km_mean = fltarr(8)*!values.f_nan
for i=0L, 7 do begin ;about channel
K_iso_10km_mean(i)=mean(K_iso_10km(*,i),/nan)
K_vol_10km_mean(i)=mean(K_vol_10km(*,i),/nan)
K_geo_10km_mean(i)=mean(K_geo_10km(*,i),/nan)
Albedo_10km_mean(i)=mean(Albedo_10km(*,i),/nan)
endfor

idx_sort= sort(wl_modis)
wl_modis_sort = wl_modis(idx_sort)
K_iso_10km_sort = K_iso_10km_mean(idx_sort)
K_geo_10km_sort = K_geo_10km_mean(idx_sort)
K_vol_10km_sort = K_vol_10km_mean(idx_sort)
Albedo_10km_sort = Albedo_10km_mean(idx_sort)


K_iso_CLARS_10km = interpol(K_iso_10km_sort,wl_modis_sort,wl_CLARS)
K_geo_CLARS_10km = interpol(K_geo_10km_sort,wl_modis_sort,wl_CLARS)
K_vol_CLARS_10km = interpol(K_vol_10km_sort,wl_modis_sort,wl_CLARS)
Albedo_CLARS_10km = interpol(Albedo_10km_sort,wl_modis_sort,wl_CLARS)






; cgplot, wl_modis, K_iso_target, psym=1, color='blue', yrange = [0.0, 1.0]
; cgplot, wl_CLARS, K_iso_CLARS, psym=1, color='red', /overplot

; cgplot, wl_modis, K_geo_target, psym=1, color='blue'
; cgplot, wl_CLARS, K_geo_CLARS, psym=1, color='red', /overplot

; cgplot, wl_modis, K_vol_target, psym=1, color='blue'
; cgplot, wl_CLARS, K_vol_CLARS, psym=1, color='red', /overplot

openw, 1, './save_MODIS_BRDF_Kernel_20190529_O2ABD.txt'
printf, 1, 'Wavelength Wavenumber K_iso K_geo K_vol Albedo'
for i=0L, 2 do begin
printf, 1, wl_clars(i), wn_clars(i), K_iso_CLARS(i), K_geo_CLARS(i), K_vol_CLARS(i), Albedo_CLARS(i), format='(6f)'
endfor
close, 1

; openw, 1, './save_MODIS_BRDF_Kernel_20190529_modis_Ch.txt'
; printf, 1, 'Wavelength Wavenumber K_iso K_geo K_vol Albedo'
; for i=0L, 2 do begin
; printf, 1, wl_clars(i), wn_clars(i), K_iso_CLARS(i), K_geo_CLARS(i), K_vol_CLARS(i), Albedo_CLARS(i), format='(6f)'
; endfor
; close, 1
;calculate BRF



stop
end
