PRO MAKE_XSEC_LUT_V300, 	 						 $
		INPUT_PATH, atmos_filename, wn_range, wn_res			   ;; INPUT



;; ========================================
;;  Test Inputs
;input_path = '../../input/'
;atmos_filename = 'test_profile.dat'
;wn_range = [[ 12949.161, 13179.192], [6161.2438, 6279.4460], [4798.1846, 4888.9297]]
;wn_range = [[ 12953.161, 13179.192], [6161.2438, 6279.4460], [4798.1846, 4888.9297]]

;loadct, 39
;device, decomposed = 0
;!p.background = -1
;!p.charsize = 2.

;;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;  	READ atmospheric profile
;;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	atmos_file = INPUT_PATH+ atmos_filename
print, atmos_file

	nlevels = float(file_lines(atmos_file))

	im_press = dblarr(nlevels) & im_press(*) = !values.f_nan
	im_temp  = dblarr(nlevels) & im_temp(*)  = !values.f_nan
	im_h2o   = dblarr(nlevels) & im_h2o(*)   = !values.f_nan
	im_co2   = dblarr(nlevels) & im_co2(*)   = !values.f_nan
	im_o2    = dblarr(nlevels) & im_o2(*)    = !values.f_nan
	im_ch4   = dblarr(nlevels) & im_ch4(*)   = !values.f_nan

	openr,1, atmos_file
	for i = 0, nlevels-1 do begin 
		readf, 1, dum1, dum2, dum3, dum4, dum5, dum6, 	$
			format = '(f13.2, 5e16.7)'
		im_press(i) = dum1
		im_temp(i)  = dum2
		im_h2o(i)   = dum3
		im_co2(i)   = dum4
		im_o2(i)    = dum5
		im_ch4(i)   = dum6
		print, im_press(i), im_temp(i), im_h2o(i), im_co2(i), im_o2(i), im_ch4(i)
	endfor 
	close,1

;;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;; 		Calculated Xsec using LUT 
;;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

gas = [ 'h2o', 'co2', 'o2', 'ch4']
gas2 = [ 'H2O', 'CO2', 'O2', 'CH4']
index = [ '01', '02', '07', '06']
nxsecs = [ 66003, 35002, 31001, 35002 ]

crs_wn_range = [ [12920, 13230], [6140, 6300], [4760, 4950] ]

wn_range = [ [12929.8400, 13220.6600], [6140.9400, 6317.5600], [4779.9600, 4922.3800] ]

for ngas = 0, n_elements(gas)-1 do begin

;	path = 'GEMSTOOL_07mar2017/GEMSTOOL_physicsdata/LUT_Xsections_data/'
	path = '/data2/OCO2/'
	path1 = path + 'ABSCO_binary_files/'

	crs0_file = path1 + gas(ngas)+'_crs0_binary.out'
	crs1_file = path1 + gas(ngas)+'_crs1_binary.out'
	crs2_file = path1 + gas(ngas)+'_crs2_binary.out'

	p_file = path1 + gas(ngas)+'_press_binary.out'
	t_file = path1 + gas(ngas)+'_temp_binary.out'
	w_file = path1 + gas(ngas)+'_wn_binary.out'

	print, '+++++++++++++++++++++++++++++++++++++'
	print, '     READ ' + gas(ngas) + ' HDF'
	print, '+++++++++++++++++++++++++++++++++++++'

	press = fltarr(64) & press(*) = !values.f_nan
    ;; 64 height layers  &  17 temperature levels
	temp = fltarr(17,64) & temp(*,*) = !values.F_nan
	wavenumber = dblarr(nxsecs(ngas)) & wavenumber(*) = !values.F_nan
	crs0 = dblarr(nxsecs(ngas), 17, 64) & crs0(*,*,*) = !values.F_nan
    ;; crs_nband = dblarr(wn, temp, pressure)
	crs1 = dblarr(nxsecs(ngas), 17, 64) & crs1(*,*,*) = !values.F_nan
	crs2 = dblarr(nxsecs(ngas), 17, 64) & crs2(*,*,*) = !values.F_nan

	openr, 2, w_file  
	readu, 2, wavenumber
	close, 2 

	openr, 3, t_file  
	readu, 3, temp
	close, 3 

	openr, 4, p_file  
	readu, 4, press
	close, 4 

	openr, 5, crs0_file  
	readu, 5, crs0
	close, 5 

	openr, 6, crs1_file  
	readu, 6, crs1
	close, 6 

	openr, 7, crs2_file  
	readu, 7, crs2
	close, 7 

;	crs = dblarr(nxsecs(ngas), 3, 17, 64) & crs(*,*,*,*) = !values.f_nan
;
;	crs(*,0,*,*) = crs0(*,*,*)
;	crs(*,1,*,*) = crs1(*,*,*)
;	crs(*,2,*,*) = crs2(*,*,*)

	wn_res = 0.01
	nwn = fltarr(3) & nwn(*)=!values.f_nan

	for nband = 0, 2 do begin 					;; 1) FOR

		nband_str = 'band' + strcompress ( string(nband+1, format='(i)'), /remove_all)

;		c1 = (double(fix(wn_range[1,nband]/wn_res,type=3))*wn_res + wn_res) - double(fix(wn_range[0,nband]/wn_res,type=3))*wn_res
;		c2 = fix(c1/wn_res, type=3)
;		nwn(nband) = c2+1


		c1 = (wn_range(1,nband)-wn_range(0,nband))
		c2 = string(c1,format='(f10.3)')
		c3 = float(c2)/wn_res
		nwn(nband) = fix(c3)+1
		print, 'Band '+string(nband+1,format='(i0.0)')+'   No. of wavenumbers '+string(nwn(nband), format='(i0.0)')

		gas_xsecs 		= dblarr(nwn(nband), nlevels) ;& gas_xsecs(*,*)=!values.F_nan
		level_dxsecs_dT = dblarr(nwn(nband), nlevels) 
		level_dxsecs_dP = dblarr(nwn(nband), nlevels)
		
		pressure = fltarr(nlevels) & pressure(*) = !values.f_nan
		temperature = fltarr(nlevels) & temperature(*) = !values.f_nan

		for n = 0, nlevels-1 do begin 			;; 2) FOR
			str_level = strcompress(string(n+1,format='(i)'),/remove_all)

	
			pressure(n) = im_press(n)*0.01
			temperature(n) = im_temp(n)
	
			result 		= size(temp)  
            ntemp       = result[1]
			new_temp	= fltarr(ntemp) & new_temp(*)  =!values.F_nan  ; 17 temperatures 
		
			for t = 0, ntemp-1 do begin 		;; 3) FOR

				x = fltarr(n_elements(press)) & x(*)=!values.F_nan
				y = fltarr(n_elements(press)) & y(*)=!values.F_nan
				count1 = 0

				for p = 0,n_elements(press)-1 do begin 
					x(count1) = press(p) 
					y(count1) = temp(t,p) 
					count1	  = count1+1 
				endfor 
				new_temp(t)=interpol(y,x,pressure(n),/spline)
                ;; new temperature interpolated by input pressure levels.. an array with size of n_elements(temp[*,0])
                ;new_temp(t)=interpol(temp[t,*], press, pressure[n])

			endfor 								;; 3) ENDFOR

			temp_interval=10.			
			temp_index=(temperature(n)-new_temp(0))/temp_interval

			x=fltarr(n_elements(press)) & x(*)=!values.F_nan
			y=fltarr(n_elements(press)) & y(*)=!values.F_nan

			count2=0

			for p=0,n_elements(press)-1	do begin
				x(count2)=press(p)
				y(count2)=count2
				count2=count2+1 
			endfor
			
			press_index=interpol(y,x,pressure(n),/spline)

			wn = wn_range(0,nband)+wn_res*indgen(nwn(nband),/double)

			index1 = where(wavenumber ge crs_wn_range(0,nband) and wavenumber le crs_wn_range(1,nband),  wn_count1)

			if wn_count1 gt 0 then begin 

				im_crs0 = crs0(index1,*,*)
				im_crs1 = crs1(index1,*,*)
				im_crs2 = crs2(index1,*,*)

				im_wavenumber = wavenumber(index1)
	
;				wn_index1 = where(im_wavenumber(*) eq wn(0), count_wn_index)
;				wn_index1 = value_locate(im_wavenumber, wn[0])
				wn_index1 = where(fix(im_wavenumber[*]/wn_res, type=3)*wn_res eq fix(wn[0]/wn_Res, type=3)*wn_Res)
				wn_index  = wn_res/0.01*indgen(nwn(nband))+wn_index1(0)

				h2o_range  	= [0.00, 0.03, 0.06]
				num_h2o 	= findgen(3)
				h2o_index 	= interpol(num_h2o, h2o_range, im_h2o(n))

				fac = 1.
				if ngas eq 1 then begin 
					if nband eq 1 then begin 
						fac = 1.0038
					endif 
					if nband eq 2 then begin 
						fac = 0.9946
					endif
				endif
				if ngas eq 2 then begin 
					fac = 1.0125
				endif

				test = dblarr(nwn(nband), 3) & test(*,*) = !values.f_nan

				test(*, 0) = interpolate(im_crs0, wn_index, temp_index, press_index(0), /grid)
				test(*, 1) = interpolate(im_crs1, wn_index, temp_index, press_index(0), /grid)
				test(*, 2) = interpolate(im_crs2, wn_index, temp_index, press_index(0), /grid)

				test_index = findgen(n_elements(wn_index))
				gas_xsecs(*,n) = interpolate(test, test_index, h2o_index, /grid) *fac 


				d = dblarr(nwn(nband), nlevels) 

				result = size(temp)
				nt = ntemp
				np = n_elements(press)
				dummy_T = dblarr(nwn(nband), nt) & dummy_T(*,*)=!values.f_nan
				dummy_P = dblarr(nwn(nband), np) & dummy_P(*,*)=!values.f_nan
				
				test_dummy_T = dblarr(nwn(nband),3) & test_dummy_T(*,*)=!values.f_nan
				test_dummy_P = dblarr(nwn(nband),3) & test_dummy_P(*,*)=!values.f_nan

				for tindex = 0, nt-1 do begin  
					test_dummy_T(*,0)	= interpolate(im_crs0, wn_index, tindex, press_index(0), /grid) 
					test_dummy_T(*,1)	= interpolate(im_crs1, wn_index, tindex, press_index(0), /grid) 
					test_dummy_T(*,2)	= interpolate(im_crs2, wn_index, tindex, press_index(0), /grid) 
					dummy_T(*,tindex)	= interpolate(test_dummy_T, test_index, h2o_index, /grid) 
				endfor
				for pindex = 0, np-1 do begin  
					test_dummy_P(*,0)	= interpolate(im_crs0, wn_index, temp_index, pindex, /grid)
					test_dummy_P(*,1)	= interpolate(im_crs1, wn_index, temp_index, pindex, /grid)
					test_dummy_P(*,2)	= interpolate(im_crs2, wn_index, temp_index, pindex, /grid)
					dummy_P(*,pindex)	= interpolate(test_dummy_P, test_index, h2o_index, /grid) 
				endfor

				dT = fltarr(nt) & dT(*)=!values.f_nan
				dxsecs_T 	= dblarr(nwn(nband),nt) & dxsecs_T(*,*)	 = !values.F_nan
				dxsecs_dT 	= dblarr(nwn(nband),nt) & dxsecs_dT(*,*) = !values.F_nan

				for tindex = 0, nt-1 do begin  
					dxsecs_T(*,tindex) 	= gas_xsecs(*,n) - dummy_T(*,tindex) 
					dT(tindex) 		 	= temperature(n) - new_temp(tindex) 
					dxsecs_dT(*,tindex) = dxsecs_T(*,tindex) / dT(tindex)
				endfor
nn=0
				sort_index = sort(abs(dT))
				dT_min_index = sort_index(nn)
				
				goto, skip1
				window, 0, xsize = 800, ysize = 500
				loadct, 39
				device, decomposed = 0
				!p.font = 1
				device,set_font='Helvetica Bold',/tt_font
				!P.background = -1  & !P.charsize = 2.
				!x.margin = [13, 4] & !y.margin = [4, 2]

				col_index = findgen(nt)*(254/nt)
				sort_index = sort(abs(dT))

				max1 = max(dxsecs_dT(*,*),/nan)
				min1 = min(dxsecs_dT(*,*),/nan)
				yint = (max1-min1)/(nt+10)
				ypos = max1-yint

				for tt = nt-1, 0, -1 do begin  
					tindex = sort_index(tt)
					tt2 = nt-1 - tt
;					print, tindex, col_index(tt2)
					
					if tt eq 6 then begin 
						plot, wn(*), dxsecs_dT(*,tindex), background=-1, col=0, $
							xtickformat = '(i)', 								$
							xstyle = 1, ystyle = 1, 							$
							ytickformat = '(e10.2)',							$
							yrange =[min1-yint, max1+yint],						$
							title = gas(ngas) +' Level '+str_level,				$ 
							xtitle = 'Wavenumber [cm!N!E-1!N]', ytitle = 'dXsecs/dT', /nodata
						oplot, wn(*), dxsecs_dT(*,tindex), col = col_index(tt2)
						str_dT = strcompress(string(dT(tindex), format='(f8.2)'),/remove_all)
						xyouts, wn(1000), ypos-tt*yint, str_dT, col = col_index(tt2)
					endif else begin
						str_dT = strcompress(string(dT(tindex), format='(f8.2)'),/remove_all)
if tt lt 6 then begin 
						oplot, wn(*), dxsecs_dT(*,tindex), col = col_index(tt2)
						xyouts, wn(1000), ypos-tt*yint, str_dT, col = col_index(tt2)
endif

					endelse
				endfor 

				device, decompose = 1
				saveimage, 'plot_dXsecs_dT_'+gas(ngas)+'_band'+nband_str+'_'+str_level+'.jpg',/jpeg
				device, decompose = 0
				skip1:


				dP = fltarr(np) & dP(*)=!values.f_nan
				dxsecs_P  = dblarr(nwn(nband), np) & dxsecs_P(*,*)  = !values.F_nan
				dxsecs_dP = dblarr(nwn(nband), np) & dxsecs_dP(*,*) = !values.F_nan

				for pindex = 0, np-1 do begin
					dxsecs_P(*,pindex) 	= gas_xsecs(*,n) - dummy_P(*,pindex) 
					dP(pindex) 			= Pressure(n) - press(pindex) 
					dxsecs_dP(*,pindex) = dxsecs_P(*,pindex)/dP(pindex) 
				endfor

				sort_index = sort(abs(dP))
				dP_min_index = sort_index(nn)
				
				goto, skip2
				window, 1, xsize = 800, ysize = 500
				col_index = findgen(np)*(256./np)
				col_index2 = findgen(20)*13
				sort_index = sort(abs(dP))

				max1 = max(dxsecs_dP(*,*),/nan)
				min1 = min(dxsecs_dP(*,*),/nan)
				yint = (max1-min1)/(30)
				ypos = max1-yint*12

				;for tt = np-1, 5, -1 do begin  
				for tt = 19, 0, -1 do begin  
					pindex = sort_index(tt)
					tt2 = 19 - tt
					;tt2 = np-1 - tt
					print, tt, dP(pindex),pindex, col_index(tt2)

					;if tt eq np-1 then begin 
					if tt eq 19 then begin 

						plot, wn(*), dxsecs_dP(*,pindex), background=-1, col=0, $
							xtickformat = '(i)', 								$
							xstyle = 1, ystyle = 1, 							$
							ytickformat = '(e10.2)',							$
							yrange =[min1-yint, max1+yint],						$
							title = gas(ngas) +' Level '+str_level,				$ 
							xtitle = 'Wavenumber [cm!N!E-1!N]', ytitle = 'dXsecs/dP', /nodata
						oplot, wn(*), dxsecs_dP(*,pindex), col = col_index2(tt2)
						str_dP = strcompress(string(dP(pindex), format='(f8.2)'),/remove_all)
						xyouts, wn(1000), ypos-tt*yint, str_dP, col = col_index2(tt2)

					endif else begin
						oplot, wn(*), dxsecs_dP(*,pindex), col = col_index2(tt2)
						str_dP = strcompress(string(dP(pindex), format='(f8.2)'),/remove_all)
						xyouts, wn(1000), ypos-tt*yint, str_dP, col = col_index2(tt2)
					endelse
				endfor

				device, decompose = 1
				saveimage, 'plot_dXsecs_dP_'+gas(ngas)+'_band'+nband_str+'_'+str_level+'.jpg',/jpeg
				device, decompose = 0
				skip2:

				level_dxsecs_dT(*,n) = dxsecs_dT(*,dT_min_index)
				level_dxsecs_dP(*,n) = dxsecs_dP(*,dP_min_index)
				;print, dT_min_index

			endif else begin
				;               print, ' no lines at band '+nband_str+' for '+gas(ngas)+' !!! '
			endelse

		endfor									;; 2) ENDFOR

;help, gas_xsecs, level_dxsecs_dT, level_dxsecs_dP
		path  = 'GEMSTOOL_07mar2017/GEMSTOOL_physicsdata/LUT_Xsections_data/' 
		path2 = path + 'results_files/'

;		crs_filename = path2 + gas2(ngas)+'_crs_'+nband_str+'_ascii.out'
;		openw, 5, crs_filename
;			for n = 0, nlevels-1 do begin           
;				printf, 5, gas_xsecs(*,n), format='(1000000e20.8)'
;			endfor
;		close, 5
	
;		dxsecs_dT_filename = path2 + gas2(ngas)+'_dXsecs_dT_'+nband_str+'_ascii.out'
;		openw, 5, dxsecs_dT_filename
;			for n = 0, nlevels-1 do begin           
;				printf, 5, level_dxsecs_dT(*,n), format='(1000000e20.8)'
;			endfor
;		close, 5

;		dxsecs_dP_filename = path2 + gas2(ngas)+'_dXsecs_dP_'+nband_str+'_ascii.out'
;		openw, 5, dxsecs_dP_filename
;			for n = 0, nlevels-1 do begin           
;				printf, 5, level_dxsecs_dP(*,n), format='(1000000e20.8)'
;			endfor
;		close, 5

		crs_filename = path2 + gas2(ngas)+'_crs_'+nband_str+'_binary.out'
		openw, 5, crs_filename
		writeu, 5, gas_xsecs
		close, 5

		dxsecs_dT_filename = path2 + gas2(ngas)+'_dXsecs_dT_'+nband_str+'_binary.out'
		openw, 5, dxsecs_dT_filename
		writeu, 5, level_dxsecs_dT
		close, 5

		dxsecs_dP_filename = path2 + gas2(ngas)+'_dXsecs_dP_'+nband_str+'_binary.out'
		openw, 5, dxsecs_dP_filename
		writeu, 5, level_dxsecs_dP
		close, 5

	DELVARX, level_dxsecs_dP, level_dxsecs_dT, $
		dxsecs_P, dP, dxsecs_dP, dxsecs_T, dT, dxsecs_dT, gas_xsecs

	endfor										;; 1) ENDFOR(gas)

	print, '!! COPY '+ gas(ngas) + 'crs file !! '

	DELVARX, level_dxsecs_dP, level_dxsecs_dT, $
		dxsecs_P, dP, dxsecs_dP, dxsecs_T, dT, dxsecs_dT, gas_xsecs

endfor

openw, 5, './GEMSTOOL_07mar2017/GEMSTOOL_NSW_wrapper/user_nwn.input'
for k=0, 2 do printf, 5, nwn[k], format='(i0.0)'
close, 5


;	DELVARX, level_dxsecs_dP, level_dxsecs_dT, $
;		dxsecs_P, dP, dxsecs_dP, dxsecs_T, dT, dxsecs_dT, gas_xsecs

END
