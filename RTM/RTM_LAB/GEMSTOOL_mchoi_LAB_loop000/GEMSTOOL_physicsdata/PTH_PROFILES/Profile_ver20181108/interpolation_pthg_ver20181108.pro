pro interpolation_pthg_ver20181108

;Myungje Choi, 11-08-2018
;original PTH is from US standard atmosphere 1976; https://www.digitaldutch.com/atmoscalc/tableatmosphere.htm
;original Trace gases profile: "DiscoverAQ_profile_1.dat" which is a default option of GEMSTOOL.


US1976_file = './PTH_us_standard_1976_res100m.dat'
discoveraq_file = './DiscoverAQ_profile_1.dat'

readcol, US1976_file, num_us, alt_us, Pres_us, Temp_us
readcol, discoveraq_file, num_daq, pres_daq, temp_daq, H2O_daq, O3_daq, NO2_daq, HCHO_daq, SO2_daq, CO_daq

alt_us = alt_us * 0.001

alt_fin = [findgen(25)*0.2, $ ;0.0~4.8, 0.2 interval (level: 25)
           5.0+findgen(11)*0.5] ; 5.0~10.0, 0.5 interval (level: 11)
pres_fin = interpol(pres_us, alt_us, alt_fin)
temp_fin = interpol(temp_us, alt_us, alt_fin) ;
;temp_fin=interpol(temp_us, alog(pres_us), alog(pres_fin)); almost identical with above interpolation < 



n_level=n_elements(alt_fin)

H2O_fin=interpol(H2O_daq,pres_daq,pres_fin) ; not spline
O3_fin=interpol(O3_daq,pres_daq,pres_fin) 
NO2_fin=interpol(NO2_daq,pres_daq,pres_fin) 
HCHO_fin=interpol(HCHO_daq,pres_daq,pres_fin) 
SO2_fin=interpol(SO2_daq,pres_daq,pres_fin) 
CO_fin=interpol(CO_daq,pres_daq,pres_fin) 

;--save Pressure-Tempereture-Gas VMR
openw, 9, './Profile_PTH_US76_GAS_discoveraq_0km_10km_ver20181108_mchoi.dat'

for i = 0L, n_elements(pres_fin)-1 do begin

	printf, 9, i+1, pres_fin(i), temp_fin(i), H2O_fin(i), O3_fin(i), NO2_fin(i), HCHO_fin(i), SO2_fin(i), CO_fin(i), $
	          f='(1(i3,3x), (1(f12.7),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x))'

endfor

close, 9

;--save layer-Altitude corresponding with layer-Pressure

openw, 9, './UVN_heightgrid_0km_10km_US76_ver20181108_mchoi.dat'

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
file_dum = '03_O3.eps'
title_dum = 'O3 (VMR)'
xvar=O3_fin
endif
if i eq 4 then begin
file_dum = '04_NO2.eps'
title_dum = 'NO2 (VMR)'
xvar=NO2_fin
endif
if i eq 5 then begin
file_dum = '05_HCHO.eps'
title_dum = 'HCHO (VMR)'
xvar=HCHO_fin
endif
if i eq 6 then begin
file_dum = '06_SO2.eps'
title_dum = 'SO2 (VMR)'
xvar=SO2_fin
endif
if i eq 7 then begin
file_dum = '07_CO.eps'
title_dum = 'CO (VMR)'
xvar=CO_fin
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