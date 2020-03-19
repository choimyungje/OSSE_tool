Pro remake_lut_crs_v5_mchoi_o2ab

; loadct,39
; device,decomposed=0

gas = [ 'h2o', 'o2']
index = [ '01', '07']

for ngas = 0, 1 do begin 

;	goto, goto_test
	hdf_pwd = '/home/mchoi/HIMAP/GEMSTOOL_mchoi_t11_00/GEMSTOOL_physicsdata/LUT_Xsections_data/ABSCO_table_hdf'
	files_name = ['h2o_hitran12.h5', $
			  'o2_v151005_cia_mlawer_v151005r1_narrow.h5'] ;     V500,  ---no available data for ch4


	;---first_Read = O2A band range data 
	file = hdf_pwd + '/'+files_name(ngas)

	print, '+++++++++++++++++++++++++++++++++'
	print, '	 READ ' + gas(ngas) + ' HDF'
	print, '+++++++++++++++++++++++++++++++++'

	file_id=H5F_OPEN(file)

	data_vmr 	= 'Broadener_01_VMR'
	dataset_id 	= H5D_OPEN(file_id, data_vmr)
	im_vmr0 		= H5D_READ(dataset_id)				;; H2O VMR [0, 0.03, 0.06]  

	data_crs	= 'Gas_'+index(ngas)+'_Absorption'
	dataset_id	= H5D_OPEN(file_id, data_crs)
	crs0			= H5D_READ(dataset_id)				;; [n_elements(wn), 3, 17, 71] = [wavenumber, temperature, pressure]
													;; crs [cm2molecule-1]
	data_index	= 'Gas_Index'
	dataset_id	= H5D_OPEN(file_id, data_index)
	gas_index0	= H5D_READ(dataset_id)      
 
	data_press	= 'Pressure'
	dataset_id	= H5D_OPEN(file_id, data_press)
	im_press	= H5D_READ(dataset_id)					;; 71 pressure levels [Pa] 
	press0		= im_press*0.01
 
	data_temp	= 'Temperature'
	dataset_id	= H5D_OPEN(file_id, data_temp)
	temp0		= H5D_READ(dataset_id)	;; [17, 71] [K]
 
	data_wn		= 'Wavenumber'
	dataset_id	= H5D_OPEN(file_id, data_wn)
	wn0			= H5D_READ(dataset_id)					;; n_elements(wn) = [31001, 16001, 19001]

	;---Second_Read = O2B band range data ; hdf format, not h5
	hdf_pwd = '/home/mchoi/HIMAP/GEMSTOOL_mchoi_t11_00/GEMSTOOL_physicsdata/LUT_Xsections_data/ABSCO_table_hdf/O2B_bands'
	files_name = ['h2o_vhitran16_14340_14590.hdf', $
		        'o2_vhitran16_14340_14590_nocia.hdf'] ;     V500,  ---no available data for ch4


	file = hdf_pwd + '/'+files_name(ngas)

	print, '+++++++++++++++++++++++++++++++++'
	print, '	 READ ' + gas(ngas) + ' HDF'
	print, '+++++++++++++++++++++++++++++++++'

	file_id=H5F_OPEN(file)
	data_vmr 	= 'Broadener_01_VMR'
	dataset_id 	= H5D_OPEN(file_id, data_vmr)
	im_vmr1 		= H5D_READ(dataset_id)				;; H2O VMR [0, 0.03, 0.06]  
      
	data_crs	= 'Gas_'+index(ngas)+'_Absorption'
	dataset_id	= H5D_OPEN(file_id, data_crs)
	crs1			= H5D_READ(dataset_id)				;; [n_elements(wn), 3, 17, 71] = [wavenumber, temperature, pressure]
													;; crs [cm2molecule-1]
	data_index	= 'Gas_Index'
	dataset_id	= H5D_OPEN(file_id, data_index)
	gas_index1	= H5D_READ(dataset_id)      
 
	data_press	= 'Pressure'
	dataset_id	= H5D_OPEN(file_id, data_press)
	im_press	= H5D_READ(dataset_id)					;; 71 pressure levels [Pa] 
	press1		= im_press*0.01
 
	data_temp	= 'Temperature'
	dataset_id	= H5D_OPEN(file_id, data_temp)
	temp1		= H5D_READ(dataset_id)	;; [17, 71] [K]
 
	data_wn		= 'Wavenumber'
	dataset_id	= H5D_OPEN(file_id, data_wn)
	wn1			= H5D_READ(dataset_id)					;; n_elements(wn) = [31001, 16001, 19001]

	;----integration to one variable
	crs=[crs0,crs1];
	wn=[wn0,wn1]
	;-----others are identical
	im_vmr=im_vmr0
	gas_idx=gas_index0
	press=press0
	temp=temp0
	

	print, '-------------------------------'
	print, '   REMAKE '+gas(ngas) + ' LUT'
	print, '-------------------------------'

	help,wn,temp,press,crs

	;	goto, skip
	crs0 = dblarr(n_elements(wn), 17, 64) & crs0(*,*,*)=!values.f_nan
	crs1 = dblarr(n_elements(wn), 17, 64) & crs1(*,*,*)=!values.f_nan
	crs2 = dblarr(n_elements(wn), 17, 64) & crs2(*,*,*)=!values.f_nan
	crs0(*,*,*) = crs(*,0,*,*)
	crs1(*,*,*) = crs(*,1,*,*)
	crs2(*,*,*) = crs(*,2,*,*)

	path = 'ABSCO_V500_O2AB_binary_files/'
	
	file_mkdir, path

	openw, 11, path + gas(ngas)+'_wn_binary.out'
	writeu, 11, wn;, format=
	close, 11
	print, wn(0:10)

	openw,  12, path + gas(ngas)+'_press_binary.out'
	writeu, 12, press
	close, 12
	print, press(10:15)

	openw,  13, path + gas(ngas)+'_temp_binary.out'
	writeu, 13, temp
	close, 13
	print, temp(0:10, 15)

	openw, 14, path + gas(ngas)+'_crs0_binary.out'
	writeu, 14, crs0
	close, 14

	openw, 14, path + gas(ngas)+'_crs1_binary.out'
	writeu, 14, crs1
	close, 14

	openw, 14, path + gas(ngas)+'_crs2_binary.out'
	writeu, 14, crs2
	close, 14

;print, crs(1000)

skip:
endfor


stop
end
