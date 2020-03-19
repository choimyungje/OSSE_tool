pro extract_kernel_from_grasp_d3_monthly_site

wn_start = [12900.00d, 14300.00d, 7700.00d]
wn_end   = [13250.00d, 14600.00d, 8000.00d]
wl_fin = 10.^7/ ((wn_start+wn_end)/2.0) ;764.818   692.042  1273.885
n_wl_Fin = n_elements(wl_Fin)

;----differently for each site
restore, '/home/mchoi/HIMAP/data_CLARS_FTS/save_target_info/save_target_info.xdr'
; save, loc_set, lat_set, lon_set, alt_set, VZA_set, VAA_set, $
;       filename = save_pwd +'/save_target_info.xdr'

work_pwd = '/home/mchoi/HIMAP/data_GRASP_POLDER'
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
K_pol_avg = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
K_pol_std = fltarr(n_wl_fin, n_mm, n_set) * !values.f_nan
NDVI_avg = fltarr(n_mm, n_set) * !values.f_nan
NDVI_std = fltarr(n_mm, n_set) * !values.f_nan


;--setup for wl interpolation using SCIAMACHY 
;---interpolated to the O2ABD channels
restore, '/home/mchoi/HIMAP/data_sfcref_for_tropomi/save_LER_SCIAMACHY.xdr'
;save, WL_SCIAMACHY, LER_SCIAMACHY_all, filename='save_LER_SCIAMACHY.xdr'

;;--wl_modis_sort
wl_POLDER = [443., 490., 565., 670., 865., 1020.]
n_wl_POLDER = n_elements(wl_POLDER)

wl_test = [wl_polder, wl_fin]
n_wl_test = n_elements(wl_test)
LER_SCIAMACHY_all_test = fltarr(n_wl_test, 12) * !values.f_nan
for i=0L, 11 do LER_SCIAMACHY_all_test(*,i) = interpol(LER_SCIAMACHY_all(*,i),WL_SCIAMACHY,wl_test)

ratio_o2b = fltarr(12)*!values.f_nan
ratio_o2a = fltarr(12)*!values.f_nan
ratio_o2d = fltarr(12)*!values.f_nan

for i=0L, 11 do ratio_o2a(i) = (LER_SCIAMACHY_all_test(6,i)-LER_SCIAMACHY_all_test(3,i))/(LER_SCIAMACHY_all_test(4,i)-LER_SCIAMACHY_all_test(3,i))
for i=0L, 11 do ratio_o2b(i) = (LER_SCIAMACHY_all_test(7,i)-LER_SCIAMACHY_all_test(3,i))/(LER_SCIAMACHY_all_test(4,i)-LER_SCIAMACHY_all_test(3,i))
for i=0L, 11 do ratio_o2d(i) = (LER_SCIAMACHY_all_test(8,i)-LER_SCIAMACHY_all_test(3,i))/(LER_SCIAMACHY_all_test(4,i)-LER_SCIAMACHY_all_test(3,i))



;--read each month MAIAC data and calculate!



for mm = 0, n_mm-1 do begin

data_pwd = work_pwd + '/data_wget_dn_20191008'

for yyyy=2008, 2013 do begin
files_yy = file_search(data_pwd, '/PARASOL_GRASP-D3_'+string(yyyy,f='(i4.4)')+'-'+string(mm+1,f='(i2.2)')+'*_sin-proj-6480x3240_V2-06.h5', count = n_file_yy)
if n_file_yy ge 1 then begin
	if yyyy eq 2008 then files = files_yy
	if yyyy gt 2008 then files = [files, files_yy] 
endif
endfor

n_files = n_elements(files)



;---extract lon,lat at first time only

if mm eq 0 then begin
	file = files(0)
	h5fid=h5f_open(file) ;---hdf5 format..
		;---Lat
		dataset_id = h5d_open(h5fid, '/L2-GRASP/Latitude')
		var_read = h5d_read(dataset_id)
		nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
		if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
		lat = var_read
		h5d_close, dataset_id
		;---Lon
		dataset_id = h5d_open(h5fid, '/L2-GRASP/Longitude')
		var_read = h5d_read(dataset_id)
		nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
		if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
		lon = var_read
		h5d_close, dataset_id
	h5f_close, h5fid
endif



; idx_close = where(abs(lon - lon_target) le 0.1 and abs(lat - lat_target) le 0.1, n_idx_close)
; ;--finding within lon-lat box
; idx = where(lon ge lon_range(0) and lon le lon_range(1) and lat ge lat_range(0) and lat le lat_range(1), n_idx)

; sz = size(lon,/dimension)
; idx_x = idx mod sz(0)
; idx_y = idx / sz(0)




; wl_POLDER = [443., 490., 565., 670., 865., 1020.]
; n_wl_POLDER = n_elements(wl_POLDER)

K_iso_polder = fltarr(n_files, n_set, n_wl_POLDER)*!values.f_nan
K_geo_polder = fltarr(n_files, n_set, n_wl_POLDER)*!values.f_nan
K_vol_polder = fltarr(n_files, n_set, n_wl_POLDER)*!values.f_nan
K_pol_polder = fltarr(n_files, n_set, n_wl_POLDER)*!values.f_nan
; NDVI_polder = fltarr(n_files, n_set)*!values.f_nan

K_iso_O2ABD = fltarr(n_files, n_set, n_wl_POLDER)*!values.f_nan
K_geo_O2ABD = fltarr(n_files, n_set, n_wl_POLDER)*!values.f_nan
K_vol_O2ABD = fltarr(n_files, n_set, n_wl_POLDER)*!values.f_nan
K_pol_O2ABD = fltarr(n_files, n_set, n_wl_POLDER)*!values.f_nan
NDVI_O2ABD = fltarr(n_files, n_set)*!values.f_nan


for i_file = 0L, n_files-1 do begin

  print, 'Month: ', mm+1, '    files: ', i_file, n_files 

	file = files(i_file)


	h5fid=h5f_open(file) ;---hdf5 format..

	; ;---Lat
	; dataset_id = h5d_open(h5fid, '/L2-GRASP/Latitude')
	; var_read = h5d_read(dataset_id)
	; nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	; if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	; lat = var_read
	; h5d_close, dataset_id

	; ;---Lon
	; dataset_id = h5d_open(h5fid, '/L2-GRASP/Longitude')
	; var_read = h5d_read(dataset_id)
	; nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	; if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	; lon = var_read
	; h5d_close, dataset_id

	;---RossLi_443_1
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi443_1')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_443_1 = var_read
	h5d_close, dataset_id
	;---RossLi_443_2
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi443_2')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_443_2 = var_read
	h5d_close, dataset_id
	;---RossLi_443_3
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi443_3')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_443_3 = var_read
	h5d_close, dataset_id

	;---RossLi_490_1
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi490_1')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_490_1 = var_read
	h5d_close, dataset_id
	;---RossLi_490_2
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi490_2')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_490_2 = var_read
	h5d_close, dataset_id
	;---RossLi_490_3
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi490_3')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_490_3 = var_read
	h5d_close, dataset_id

	;---RossLi_565_1
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi565_1')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_565_1 = var_read
	h5d_close, dataset_id
	;---RossLi_565_2
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi565_2')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_565_2 = var_read
	h5d_close, dataset_id
	;---RossLi_565_3
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi565_3')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_565_3 = var_read
	h5d_close, dataset_id

	;---RossLi_670_1
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi670_1')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_670_1 = var_read
	h5d_close, dataset_id
	;---RossLi_670_2
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi670_2')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_670_2 = var_read
	h5d_close, dataset_id
	;---RossLi_670_3
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi670_3')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_670_3 = var_read
	h5d_close, dataset_id

	;---RossLi_865_1
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi865_1')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_865_1 = var_read
	h5d_close, dataset_id
	;---RossLi_865_2
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi865_2')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_865_2 = var_read
	h5d_close, dataset_id
	;---RossLi_865_3
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi865_3')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_865_3 = var_read
	h5d_close, dataset_id

	;---RossLi_1020_1
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi1020_1')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_1020_1 = var_read
	h5d_close, dataset_id
	;---RossLi_1020_2
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi1020_2')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_1020_2 = var_read
	h5d_close, dataset_id
	;---RossLi_1020_3
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBRDFRossLi1020_3')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	RossLi_1020_3 = var_read
	h5d_close, dataset_id

	;---BPDF_443
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBPDFMaignanBreon443')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	BPDF_443 = var_read
	h5d_close, dataset_id
	;---BPDF_490
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBPDFMaignanBreon490')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	BPDF_490 = var_read
	h5d_close, dataset_id
	;---BPDF_565
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBPDFMaignanBreon565')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	BPDF_565 = var_read
	h5d_close, dataset_id
	;---BPDF_670
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBPDFMaignanBreon670')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	BPDF_670 = var_read
	h5d_close, dataset_id
	;---BPDF_865
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBPDFMaignanBreon865')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	BPDF_865 = var_read
	h5d_close, dataset_id
	;---BPDF_1020
	dataset_id = h5d_open(h5fid, '/L2-GRASP/LandBPDFMaignanBreon1020')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	BPDF_1020 = var_read
	h5d_close, dataset_id

	;---NDVI
	dataset_id = h5d_open(h5fid, '/L2-GRASP/NDVI')
	var_read = h5d_read(dataset_id)
	nan_idx = where(var_Read eq -9999.0d, n_nan_idx)
	if n_nan_idx ge 1 then var_read(nan_idx) = !values.d_nan
	NDVI = var_read
	h5d_close, dataset_id

	h5f_close, h5fid


K_iso_read = [[[RossLi_443_1]],[[RossLi_490_1]],[[RossLi_565_1]],[[RossLi_670_1]],[[RossLi_865_1]],[[RossLi_1020_1]]]
K_vol_read = [[[RossLi_443_2]],[[RossLi_490_2]],[[RossLi_565_2]],[[RossLi_670_2]],[[RossLi_865_2]],[[RossLi_1020_2]]]
K_geo_read = [[[RossLi_443_3]],[[RossLi_490_3]],[[RossLi_565_3]],[[RossLi_670_3]],[[RossLi_865_3]],[[RossLi_1020_3]]]
K_pol_read = [[[BPDF_443]],[[BPDF_490]],[[BPDF_565]],[[BPDF_670]],[[BPDF_865]],[[BPDF_1020]]]

;---interpolate each date polder wl data to site

for i_wl = 0L, n_wl_POLDER-1 do begin

K_iso_select = K_iso_read(*,*,i_wl)
usable_idx = where(finite(K_iso_select) eq 1 and finite(lon) eq 1 and finite(lat) eq 1, n_usable_idx)
if n_usable_idx ge 1 then begin
K_iso_polder(i_file,*,i_wl) = griddata(lon(usable_idx),lat(usable_idx),K_iso_select(usable_idx), xout=lon_set, yout=lat_set)
endif

K_vol_select = K_vol_read(*,*,i_wl)
usable_idx = where(finite(K_vol_select) eq 1 and finite(lon) eq 1 and finite(lat) eq 1, n_usable_idx)
if n_usable_idx ge 1 then begin
K_vol_polder(i_file,*,i_wl) = griddata(lon(usable_idx),lat(usable_idx),K_vol_select(usable_idx), xout=lon_set, yout=lat_set)
endif

K_geo_select = K_geo_read(*,*,i_wl)
usable_idx = where(finite(K_geo_select) eq 1 and finite(lon) eq 1 and finite(lat) eq 1, n_usable_idx)
if n_usable_idx ge 1 then begin
K_geo_polder(i_file,*,i_wl) = griddata(lon(usable_idx),lat(usable_idx),K_geo_select(usable_idx), xout=lon_set, yout=lat_set)
endif

K_pol_select = K_pol_read(*,*,i_wl)
usable_idx = where(finite(K_pol_select) eq 1 and finite(lon) eq 1 and finite(lat) eq 1, n_usable_idx)
if n_usable_idx ge 1 then begin
K_pol_polder(i_file,*,i_wl) = griddata(lon(usable_idx),lat(usable_idx),K_pol_select(usable_idx), xout=lon_set, yout=lat_set)
endif

if i_wl eq 0 then begin
	NDVI_select = NDVI(*,*)
	usable_idx = where(finite(NDVI_select) eq 1 and finite(lon) eq 1 and finite(lat) eq 1, n_usable_idx)
	if n_usable_idx ge 1 then begin
	NDVI_O2ABD(i_file,*) = griddata(lon(usable_idx),lat(usable_idx),NDVI_select(usable_idx), xout=lon_set, yout=lat_set)
	endif
endif

endfor ; for i_wl = 0L, n_wl_POLDER-1 do begin

;---interpolate polder channels to O2ABD channels using each month SCIAMACHY measurement

  K_iso_O2ABD(i_file,*,0)  = K_iso_polder(i_file,*,3) * (1. - ratio_o2a(mm))+ K_iso_polder(i_file,*,4) * ratio_o2a(mm)
  K_vol_O2ABD(i_file,*,0)  = K_vol_polder(i_file,*,3) * (1. - ratio_o2a(mm))+ K_vol_polder(i_file,*,4) * ratio_o2a(mm)
  K_geo_O2ABD(i_file,*,0)  = K_geo_polder(i_file,*,3) * (1. - ratio_o2a(mm))+ K_geo_polder(i_file,*,4) * ratio_o2a(mm)
  K_pol_O2ABD(i_file,*,0)  = K_pol_polder(i_file,*,3) * (1. - ratio_o2a(mm))+ K_pol_polder(i_file,*,4) * ratio_o2a(mm)
  
  K_iso_O2ABD(i_file,*,1)  = K_iso_polder(i_file,*,3) * (1. - ratio_o2b(mm))+ K_iso_polder(i_file,*,4) * ratio_o2b(mm)
  K_vol_O2ABD(i_file,*,1)  = K_vol_polder(i_file,*,3) * (1. - ratio_o2b(mm))+ K_vol_polder(i_file,*,4) * ratio_o2b(mm)
  K_geo_O2ABD(i_file,*,1)  = K_geo_polder(i_file,*,3) * (1. - ratio_o2b(mm))+ K_geo_polder(i_file,*,4) * ratio_o2b(mm)
  K_pol_O2ABD(i_file,*,1)  = K_pol_polder(i_file,*,3) * (1. - ratio_o2b(mm))+ K_pol_polder(i_file,*,4) * ratio_o2d(mm)
  
  K_iso_O2ABD(i_file,*,2)  = K_iso_polder(i_file,*,3) * (1. - ratio_o2d(mm))+ K_iso_polder(i_file,*,4) * ratio_o2d(mm)
  K_vol_O2ABD(i_file,*,2)  = K_vol_polder(i_file,*,3) * (1. - ratio_o2d(mm))+ K_vol_polder(i_file,*,4) * ratio_o2d(mm)
  K_geo_O2ABD(i_file,*,2)  = K_geo_polder(i_file,*,3) * (1. - ratio_o2d(mm))+ K_geo_polder(i_file,*,4) * ratio_o2d(mm)
  K_pol_O2ABD(i_file,*,2)  = K_pol_polder(i_file,*,3) * (1. - ratio_o2d(mm))+ K_pol_polder(i_file,*,4) * ratio_o2d(mm)
  
  ;NDVI_O2ABD is already saved before.

endfor ; for i_file = 0L, n_files-1 do begin

;---calculate avg and std of monthly kernel info of each site


for i_set=0L, n_set-1 do begin
for i_ch=0L, n_wl_fin-1 do begin
  K_iso_avg(i_ch,mm,i_set) =   mean(K_iso_O2ABD(*,i_set,i_ch),/nan)
  K_iso_std(i_ch,mm,i_set) = stddev(K_iso_O2ABD(*,i_set,i_ch),/nan)
  K_vol_avg(i_ch,mm,i_set) =   mean(K_vol_O2ABD(*,i_set,i_ch),/nan)
  K_vol_std(i_ch,mm,i_set) = stddev(K_vol_O2ABD(*,i_set,i_ch),/nan)
  K_geo_avg(i_ch,mm,i_set) =   mean(K_geo_O2ABD(*,i_set,i_ch),/nan)
  K_geo_std(i_ch,mm,i_set) = stddev(K_geo_O2ABD(*,i_set,i_ch),/nan)
  K_pol_avg(i_ch,mm,i_set) =   mean(K_pol_O2ABD(*,i_set,i_ch),/nan)
  K_pol_std(i_ch,mm,i_set) = stddev(K_pol_O2ABD(*,i_set,i_ch),/nan)
endfor
NDVI_avg(mm,i_set) =   mean(NDVI_O2ABD(*,i_set),/nan)
NDVI_std(mm,i_set) = stddev(NDVI_O2ABD(*,i_set),/nan)
endfor



endfor ;for mm = 0, n_mm-1 do begin


save, K_iso_avg, K_iso_std, K_vol_avg, K_vol_std, K_geo_avg, K_geo_std, $
	K_pol_avg, K_pol_std, NDVI_avg, NDVI_std, $
      filename = save_pwd + '/save_POLDER_BRDF_month_site.xdr'



stop
end