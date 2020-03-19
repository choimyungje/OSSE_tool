pro interpolation_pthg_for_nsw

;Myungje Choi, 11-15-2018
;recalculation for NSW
;adopting "Initial_NSW_atmospehre_12aug13" file but more grid.

;Myungje Choi, 11-08-2018
;original PTH is from US standard atmosphere 1976; https://www.digitaldutch.com/atmoscalc/tableatmosphere.htm
;original Trace gases profile: "DiscoverAQ_profile_1.dat" which is a default option of GEMSTOOL.


US1976_file = './PTH_us_standard_1976_res100m.dat'
discoveraq_file = './DiscoverAQ_profile_1.dat'

readcol, US1976_file, num_us, alt_us, Pres_us, Temp_us
readcol, discoveraq_file, num_daq, pres_daq, temp_daq, H2O_daq, O3_daq, NO2_daq, HCHO_daq, SO2_daq, CO_daq

alt_us = alt_us * 0.001
; alt_fin = [findgen(25)*0.2, $ ;0.0~4.8, 0.2 interval (level: 25)
;            5.0+findgen(11)*0.5] ; 5.0~10.0, 0.5 interval (level: 11)

alt_fin = [findgen(51)*0.2] ;200 m resolution, 51 levels.

;---read another profile file: 
Initial_NSW_file = 'Initial_NSW_atmosphere_12aug13.dat'
readcol, Initial_NSW_file, pres_nswi, temp_nswi, h2o_nswi, co2_nswi, o2_nswi, ch4_nswi, co21_nswi, co22_nswi




;---recalculation approach
;---1) P, T, H = US standard
;---2) H2O, CO2, O2, CH4, CO21, CO22 --> linear interpolation using Initial_NSW_file


pres_fin = interpol(pres_us, alt_us, alt_fin)
temp_fin = interpol(temp_us, alt_us, alt_fin) ;
;temp_fin=interpol(temp_us, alog(pres_us), alog(pres_fin)); almost identical with above interpolation < 


n_level=n_elements(alt_fin)

; H2O_fin=interpol(H2O_daq,pres_daq,pres_fin) ; not spline
; O3_fin=interpol(O3_daq,pres_daq,pres_fin) 
; NO2_fin=interpol(NO2_daq,pres_daq,pres_fin) 
; HCHO_fin=interpol(HCHO_daq,pres_daq,pres_fin) 
; SO2_fin=interpol(SO2_daq,pres_daq,pres_fin) 
; CO_fin=interpol(CO_daq,pres_daq,pres_fin) 

H2O_fin=interpol(h2o_nswi,pres_nswi,pres_fin) ; not spline
CO2_fin=interpol(co2_nswi,pres_nswi,pres_fin) ; not spline
O2_fin =interpol(o2_nswi, pres_nswi,pres_fin) ; not spline
CH4_fin=interpol(h2o_nswi,pres_nswi,pres_fin) ; not spline
CO21_fin=interpol(co21_nswi,pres_nswi,pres_fin) ; not spline
CO22_fin=interpol(co22_nswi,pres_nswi,pres_fin) ; not spline



;--save Pressure-Tempereture-Gas VMR
openw, 9, './atmos_mchoi.dat'
; for i = 0L, n_elements(pres_fin)-1 do begin
; 	printf, 9, i+1, pres_fin(i), temp_fin(i), H2O_fin(i), O3_fin(i), NO2_fin(i), HCHO_fin(i), SO2_fin(i), CO_fin(i), $
; 	          f='(1(i3,3x), (1(f12.7),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x))'
; endfor

for i = n_elements(pres_fin)-1, 0L, -1 do begin
	printf, 9, pres_fin(i)*100., temp_fin(i), H2O_fin(i), CO2_fin(i), O2_fin(i), CH4_fin(i), CO21_fin(i), CO22_fin(i), $
	          f='((1(f13.2),x),(7(e14.4),x))'
endfor
close, 9



;--save layer-Altitude corresponding with layer-Pressure
openw, 9, './NSW_heightgrid_0km_10km_mchoi.dat'
for i = n_elements(pres_fin)-1, 0, -1 do begin
	printf, 9, 0, alt_fin(i), $
	          f='(1(i12,3x), (1(f7.4)))'
endfor
close, 9

;--plot each data
;--Y-axis=alt_fin; X-axis= Temp_fin, Pres_fin, H2O_fin, O3_fin, NO2_fin, HCHO_fin, SO2_fin, CO_fin
For i = 0L, 7 do begin

if i eq 0 then begin
file_dum = '00_temp.eps'
title_dum = 'Temperature (K)'
xvar=Temp_fin
endif
if i eq 1 then begin
file_dum = '01_pres.eps'
title_dum = 'Pressure (hPa)'
xvar=Pres_fin
endif
if i eq 2 then begin
file_dum = '02_H2O.eps'
title_dum = 'H2O (VMR)'
xvar=H2O_fin
endif
if i eq 3 then begin
file_dum = '03_CO2.eps'
title_dum = 'CO2 (VMR)'
xvar=CO2_fin
endif
if i eq 4 then begin
file_dum = '04_O2.eps'
title_dum = 'O2 (VMR)'
xvar=O2_fin
endif
if i eq 5 then begin
file_dum = '05_CH4.eps'
title_dum = 'CH4 (VMR)'
xvar=CH4_fin
endif
if i eq 6 then begin
file_dum = '06_CO2_1.eps'
title_dum = 'CO2_1 (VMR)'
xvar=CO21_fin
endif
if i eq 7 then begin
file_dum = '07_CO2_2.eps'
title_dum = 'CO2_2 (VMR)'
xvar=CO22_fin
endif

yvar=alt_fin

xrange_dum=[min(xvar,/nan),max(xvar,/nan)]
yrange_dum=[min(yvar,/nan),max(yvar,/nan)]

; usersym for circle filled
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A) , SIN(A) , THICK=0.5, /FILL

loadct,39,/silent
!p.font=0
set_plot,'ps'

device, /encapsulated, /color, bits_per_pixel = 8, filename = file_dum, $
        /helvetica, /bold, font_size = 15, xs = 20, ys = 30

plot, xvar, yvar, psym = -8, xrange=xrange_dum, yrange=yrange_dum, $
      background=0, color=0, $
	  xtitle=title_dum, ytitle='Altitude (km)', charsize=2,charthick=1.5, $
	  xstyle=0, ystyle=1, xmargin=[8,4], ymargin=[4,4], /nodata

oplot, xvar, yvar, color=50, psym = -8, symsize=2, thick=5

device,/close




endfor




stop


end