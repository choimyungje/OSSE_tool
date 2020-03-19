pro interpolation_bpdf

lat_target = 34.1101d
lon_target = -117.9695d

yyyy = 2013
mm = 05
str_yyyymm = string(yyyy,f='(i4.4)')+string(mm,f='(i2.2)')

L3_pwd = '/home/mchoi/HIMAP/data_GRASP_POLDER/data_L3'
files = file_search(L3_pwd + '/GRASP_POLDER_L3_'+str_yyyymm+'.01degree.nc', count = n_files)
;or climatological value


file = files(0)

;--read lat lon lev info (identical)
ncid=NCDF_open(file)
ncdf_varget, ncid, 'Latitude',lat
ncdf_varget, ncid, 'Longitude',lon
ncdf_varget, ncid, 'NDVI',NDVI_read
ncdf_varget, ncid, 'LandBPDFMaignanBreon670',K_BPDF_read

ncdf_varget, ncid, 'Ross_Li_BRDF_443_isotropic_parameter',K_iso_443_read
ncdf_varget, ncid, 'Ross_Li_BRDF_490_isotropic_parameter',K_iso_490_read
ncdf_varget, ncid, 'Ross_Li_BRDF_565_isotropic_parameter',K_iso_565_read
ncdf_varget, ncid, 'Ross_Li_BRDF_670_isotropic_parameter',K_iso_670_read
ncdf_varget, ncid, 'Ross_Li_BRDF_865_isotropic_parameter',K_iso_865_read
ncdf_varget, ncid, 'Ross_Li_BRDF_1020_isotropic_parameter',K_iso_1020_read

ncdf_varget, ncid, 'Ross_Li_BRDF_670_geometric_parameter',K_geo_670_read
ncdf_varget, ncid, 'Ross_Li_BRDF_670_volumetric_parameter',K_vol_670_read


NCDF_close, ncid
n_lat=n_elements(lat)
n_lon=n_elements(lon)
; n_NDVI=n_elements(NDVI)
; n_time=n_elements(time)

sz = size(lat, /dimension)

idx = where(abs(lat_target - lat) lt 0.3 and abs(lon_target-lon) lt 0.3, n_idx)
if n_idx ge 1 then begin
  distance=fltarr(n_idx)*!values.f_nan
  for i=0L, n_idx-1 do begin
      distance(i)=map_2points(lon_target,lat_target,lon(idx(i)),lat(idx(i)),/meters)
  endfor
idx_closest=where(distance eq min(distance,/nan), n_idx_closest)
idx2=idx(idx_closest(0))
idx2_x = idx2 mod sz(0)
idx2_y = idx2 / sz(0)

;---calculate within 10 km
idx_10km=where(distance le 10000. , n_idx_10km)
idx3=idx(idx_10km)
idx_10km_x = idx3 mod sz(0)
idx_10km_y = idx3 / sz(0)

endif


wl_POLDER = [443., 490., 565., 670., 865., 1020.]
n_wl_POLDER = n_elements(wl_POLDER)

NDVI = replicate(NDVI_read(idx2_x,idx2_y), n_wl_POLDER)
K_BPDF = replicate(K_BPDF_read(idx2_x,idx2_y), n_wl_POLDER)

K_iso = [K_iso_443_read(idx2_x,idx2_y), $
	   K_iso_490_read(idx2_x,idx2_y), $
	   K_iso_565_read(idx2_x,idx2_y), $
	   K_iso_670_read(idx2_x,idx2_y), $
	   K_iso_865_read(idx2_x,idx2_y), $
	   K_iso_1020_read(idx2_x,idx2_y)]

K_geo = replicate(K_geo_670_read(idx2_x,idx2_y), n_wl_POLDER)
K_vol = replicate(K_vol_670_read(idx2_x,idx2_y), n_wl_POLDER)


restore, '/home/mchoi/HIMAP/data_MAIAC_BRDF/save_K_modis.xdr'

stop


end