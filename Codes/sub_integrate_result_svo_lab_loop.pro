Pro sub_integrate_result_svo_lab_loop, result_pwd

; if i_band eq 0 then str_range = '12900.00_13250.00'
; if i_band eq 1 then str_range = '07800.00_08000.00'

; ;-----case) all data are calculated well without missing
; result_pwd = '/home/mchoi/HIMAP/workspace_t16_12_LABS_lamb/Results'+$
;              '/SVO_AerosolProf_step'+string(i_step,f='(i2.2)')+'/'+str_range

; ; result_pwd = '/home/mchoi/HIMAP/workspace_t16_12_LABS_lamb/Results'+$
; ;              '/LABS_BRDF_BPDF_AerosolProf_step'+string(i_step,f='(i2.2)')+'_ABSCOv5/12900.00_13250.00'
; ; result_pwd = '/home/mchoi/HIMAP/workspace_t16_12_LABS_lamb/Results'+$
; ;              '/LAB_AerosolProf_step'+string(i_step,f='(i2.2)')+'_ABSCOv5/07700.00_08000.00'

file_str = ['Stokes_IQU_DoLP_TOAUP.out', $
            'Aer_prof_edit.dat', 'Aer_SSA.dat', $
            'Optical_depth_all.dat', 'Optical_depth_gases.dat', $
            'Stokes_I_AerBulkWFs_TOAUP.out', 'Stokes_Q_AerBulkWFs_TOAUP.out', 'Stokes_U_AerBulkWFs_TOAUP.out', $
            'Stokes_I_AlbWFs_TOAUP.out',   'Stokes_Q_AlbWFs_TOAUP.out',   'Stokes_U_AlbWFs_TOAUP.out', $
            'Stokes_I_SURFPWFs_TOAUP.out', 'Stokes_Q_SURFPWFs_TOAUP.out', 'Stokes_U_SURFPWFs_TOAUP.out', $
            'Stokes_I_TSHIFTWFs_TOAUP.out', 'Stokes_Q_TSHIFTWFs_TOAUP.out', 'Stokes_U_TSHIFTWFs_TOAUP.out', $
            'Stokes_I_BRDFWFs_TOAUP.out', 'Stokes_Q_BRDFWFs_TOAUP.out', 'Stokes_U_BRDFWFs_TOAUP.out', $ ;
            'Stokes_I_WSCALEWFs_TOAUP.out', 'Stokes_Q_WSCALEWFs_TOAUP.out', 'Stokes_U_WSCALEWFs_TOAUP.out', $ ;--H2O scaling
            'Stokes_I_AerProfWFs_TOAUP.out', 'Stokes_Q_AerProfWFs_TOAUP.out', 'Stokes_U_AerProfWFs_TOAUP.out', $ ;--AOD profile
            'Stokes_I_GasWFs_TOAUP.outH2O', 'Stokes_I_GasWFs_TOAUP.outH2O', 'Stokes_I_GasWFs_TOAUP.outH2O', $ ;--H2O profile
            'Stokes_I_GasWFs_TOAUP.outO2', 'Stokes_Q_GasWFs_TOAUP.outO2', 'Stokes_U_GasWFs_TOAUP.outO2', $ ;--O2 profile
            'dump']

for i = 0L, n_elements(file_str)-1 do begin

files = file_search(result_pwd+'/Loop*/'+file_str(i), count=n_files)
if n_files ge 1 then begin
for f=0L, n_files-1 do begin
  openr,1,files(f)
  n_lines = file_lines(files(f))
  ; print, n_lines
  data=strarr(n_lines)
  readf,1,data
  close,1
  ; print, data(n_lines-1)
  if f eq 0 then data_sum = data
  if f ge 1 then data_Sum=[data_sum,data]
endfor

;--print, 
openw, 1, result_pwd+'/'+file_str(i)
for l=0L, n_elements(data_sum)-1 do begin
printf, 1, data_sum(l)
endfor
close, 1
endif ;if n_files ge 1 then begin
endfor ;for i = 0L, n_elements(file_str)-1 do begin

;---copy other input conditions
file_str0 = ['Input_aerosol.dat','Input_alt_pres.out','Input_surface_lambertian.dat','Input_BRDF_Kernel_Factor_Parameter.dat']
for i = 0L, n_elements(file_str0)-1 do begin

files = file_search(result_pwd+'/Loop*/'+file_str0(i), count=n_files)
if n_files ge 1 then begin
file_copy, files(0), result_pwd+'/'+file_str0(i), /OVERWRITE

endif ;if n_files ge 1 then begin

endfor ;for i = 0L, n_elements(file_str)-1 do begin

file = file_Search(result_pwd+'/Stokes_I_AlbWFs_TOAUP.out', count=n_file)
if n_file ge 1 then begin
;---additional calculation for albedo slope jacobian
readcol, result_pwd+'/Stokes_I_AlbWFs_TOAUP.out', temp, wn, JI_albedo,/silent
JI_slope = JI_albedo*!values.d_nan
n_wn = n_elements(wn)
openw, 1, result_pwd+'/Stokes_I_AlbWFs_TOAUP_slope.out'
for i_wn=0L, n_wn-1 do begin
fac = wn(i_wn) - mean([wn(0),wn(n_wn-1)])
JI_slope(i_wn)=JI_albedo(i_wn)*fac
printf, 1, i_wn+1, wn(i_wn), JI_albedo(i_wn), JI_slope(i_wn), format='(i5, f10.2, 1000E20.10)'
endfor
close, 1
endif

file = file_Search(result_pwd+'/Stokes_Q_AlbWFs_TOAUP.out', count=n_file)
if n_file ge 1 then begin
readcol, result_pwd+'/Stokes_Q_AlbWFs_TOAUP.out', temp, wn, JQ_albedo,/silent
JQ_slope = JQ_albedo*!values.d_nan
n_wn = n_elements(wn)
openw, 1, result_pwd+'/Stokes_Q_AlbWFs_TOAUP_slope.out'
for i_wn=0L, n_wn-1 do begin
fac = wn(i_wn) - mean([wn(0),wn(n_wn-1)])
JQ_slope(i_wn)=JQ_albedo(i_wn)*fac
printf, 1, i_wn+1, wn(i_wn), JQ_albedo(i_wn), JQ_slope(i_wn), format='(i5, f12.5, 100E20.10)'
endfor
close, 1
endif

file = file_Search(result_pwd+'/Stokes_U_AlbWFs_TOAUP.out', count=n_file)
if n_file ge 1 then begin
readcol, result_pwd+'/Stokes_U_AlbWFs_TOAUP.out', temp, wn, JU_albedo,/silent
JU_slope = JU_albedo*!values.d_nan
n_wn = n_elements(wn)
openw, 1, result_pwd+'/Stokes_U_AlbWFs_TOAUP_slope.out'
for i_wn=0L, n_wn-1 do begin
fac = wn(i_wn) - mean([wn(0),wn(n_wn-1)])
JU_slope(i_wn)=JU_albedo(i_wn)*fac
printf, 1, i_wn+1, wn(i_wn), JU_albedo(i_wn), JU_slope(i_wn), format='(i5, f12.2, 100E20.10)'
endfor
close, 1
endif


; ;-----case) if there are missing data... unforunately. -> just one case
; result_pwd = '/home/mchoi/HIMAP/GEMSTOOL_mchoi_t16_02_LABS_lamb/GEMSTOOL_NSW_wrapper/GEMSTOOL_NSW_Results'+ $
;              '/LABS_lamb_MAIAC_AOD0.1142_HITRAN2016/07700.00_08000.00'
; result_backup_pwd = '/home/mchoi/HIMAP/GEMSTOOL_mchoi_t16_02_LABS_lamb/GEMSTOOL_NSW_wrapper/GEMSTOOL_NSW_Results'+ $
;              '/LABS_lamb_MAIAC_AOD0.1142_HITRAN2016/07719.98_07979.98'

; file_str = ['Stokes_IQU_DoLP_TOAUP.out', $
;             'Aer_prof_edit.dat', 'Aer_SSA.dat', $
;             'Optical_depth_all.dat', 'Optical_depth_gases.dat', $
;             'Stokes_I_AerBulkWFs_TOAUP.out', 'Stokes_Q_AerBulkWFs_TOAUP.out', 'Stokes_U_AerBulkWFs_TOAUP.out', $
;             'Stokes_I_AlbWFs_TOAUP.out',   'Stokes_Q_AlbWFs_TOAUP.out',   'Stokes_U_AlbWFs_TOAUP.out', $
;             'Stokes_I_SURFPWFs_TOAUP.out', 'Stokes_Q_SURFPWFs_TOAUP.out', 'Stokes_U_SURFPWFs_TOAUP.out', $
;             'Stokes_I_TSHIFTWFs_TOAUP.out', 'Stokes_Q_TSHIFTWFs_TOAUP.out', 'Stokes_U_TSHIFTWFs_TOAUP.out', $
;             'Stokes_I_BRDFWFs_TOAUP.out', 'Stokes_Q_BRDFWFs_TOAUP.out', 'Stokes_U_BRDFWFs_TOAUP.out']


; for i = 0L, n_elements(file_str)-1 do begin

; backup_files = file_search(result_backup_pwd+'/Loop*/'+file_str(i), count=n_backup_files)
; if n_backup_files ge 1 then begin
; openr,1,backup_files(0)
; n_lines = file_lines(backup_files(0))
; ; print, n_lines
; data_backup=strarr(n_lines)
; readf,1,data_backup
; close,1
; endif 

; files = file_search(result_pwd+'/Loop*/'+file_str(i), count=n_files)
; if n_files ge 1 then begin
; for f=0L, n_files-1 do begin
;   openr,1,files(f)
;   n_lines = file_lines(files(f))
;   ; print, n_lines
;   data=strarr(n_lines)
;   readf,1,data
;   close,1
;   ; print, data(n_lines-1)
;   if f eq 0 then begin
;     data_sum = data
;     data_sum = [data_sum,data_backup(f)]
;   endif
;   if f ge 1 and f le n_files-2 then begin
;   data_Sum=[data_sum,data]
;   data_sum = [data_sum,data_backup(f)]
;   endif
;   if f eq n_files-1 then begin
;   data_Sum=[data_sum,data]
;   endif
; endfor


; ;--print, 
; openw, 1, result_pwd+'/'+file_str(i)
; for l=0L, n_elements(data_sum)-1 do begin
; printf, 1, data_sum(l)
; endfor
; close, 1
; endif ;if n_files ge 1 then begin
; endfor ;for i = 0L, n_elements(file_str)-1 do begin

; ;---copy other input conditions
; file_str0 = ['Input_aerosol.dat','Input_alt_pres.out','Input_surface_lambertian.dat','Input_BRDF_Kernel_Factor_Parameter.dat']
; for i = 0L, n_elements(file_str0)-1 do begin

; files = file_search(result_pwd+'/Loop*/'+file_str0(i), count=n_files)
; if n_files ge 1 then begin
; file_copy, files(0), result_pwd+'/'+file_str0(i)

; endif ;if n_files ge 1 then begin

; endfor ;for i = 0L, n_elements(file_str)-1 do begin


end