
loadct,39
device,decomposed=0


;goto, skip


;gas = [ 'H2O', 'CO2', 'O2', 'CH4']
gas = [ 'h2o', 'co2', 'o2', 'ch4']
index = [ '01', '02', '07', '06']

for ngas = 0, 3 do begin 

;	goto, goto_test

	nfile = findfile('files/' + gas(ngas) + '_v4.2.0*.hdf', count = c)
	file = nfile(0)

	print, '+++++++++++++++++++++++++++++++++'
	print, '	 READ ' + gas(ngas) + ' HDF'
	print, '+++++++++++++++++++++++++++++++++'

	file_id=H5F_OPEN(file)

	data_vmr 	= 'Broadener_01_VMR'
	dataset_id 	= H5D_OPEN(file_id, data_vmr)
	im_vmr 		= H5D_READ(dataset_id)				;; H2O VMR [0, 0.03, 0.06]  

	data_crs	= 'Gas_'+index(ngas)+'_Absorption'
	dataset_id	= H5D_OPEN(file_id, data_crs)
	crs			= H5D_READ(dataset_id)				;; [n_elements(wn), 3, 17, 71] = [wavenumber, temperature, pressure]
													;; crs [cm2molecule-1]
	data_index	= 'Gas_Index'
	dataset_id	= H5D_OPEN(file_id, data_index)
	gas_index	= H5D_READ(dataset_id)      
 
	data_press	= 'Pressure'
	dataset_id	= H5D_OPEN(file_id, data_press)
	im_press	= H5D_READ(dataset_id)					;; 71 pressure levels [Pa] 
	press		= im_press*0.01
 
	data_temp	= 'Temperature'
	dataset_id	= H5D_OPEN(file_id, data_temp)
	temp		= H5D_READ(dataset_id)	;; [17, 71] [K]
 
	data_wn		= 'Wavenumber'
	dataset_id	= H5D_OPEN(file_id, data_wn)
	wn			= H5D_READ(dataset_id)					;; n_elements(wn) = [31001, 16001, 19001]

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

	path = 'ABSCO_binary_files/'
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
