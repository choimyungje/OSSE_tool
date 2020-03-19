pro Profile_PT_from_US1976_and_G_from_discoverAQ_and_H_from_vlidort_0km


;alt = 0.0   ;km unit
;pressure = -127.515 + 1141.65*exp(-0.107084*alt)                 ; hPa unit


fixed_ph_file = 'D:\GEMSTOOL\GEMSTOOL_ver1\GEMSTOOL_ver1\GEMSTOOL_UVN_wrapper\UVN_heightgrid_1.dat'
discoveraq_ph_file = 'D:\GEMSTOOL\GEMSTOOL_ver1\GEMSTOOL_ver1\GEMSTOOL_physicsdata\PTH_PROFILES\DiscoverAQ_profile_1_editted_to_US76_T_number_only.dat'

readcol, fixed_ph_file, num_fixed, alti_fixed
readcol, discoveraq_ph_file, num_daq, pres_daq, temp_daq, H2O_daq, O3_daq, NO2_daq, HCHO_daq, SO2_daq, CO_daq

alti_fixed_inv = fltarr(n_elements(alti_fixed)) & alti_fixed_inv(*) = !values.f_nan

for i=0L, n_elements(alti_fixed)-1 do begin
alti_fixed_inv(i) = alti_fixed(n_elements(alti_fixed)-1  - i)
endfor


pres_mchoi = pres_daq
pres_mchoi(0) = 1013.25
pres_mchoi(32) = 540.88



pth_us1976_file = 'D:\GOCI_OP\LUT_vlidort\htp_us1976.txt'
readcol, pth_us1976_file, alti_us76, temp_us76, pres_us76


temp_mchoi = interpol(temp_us76, alog(pres_us76), alog(pres_mchoi))    ; 실제로 해보니 fitting이 잘맞음


n_layers = n_elements(pres_mchoi)

alti_mchoi = fltarr(n_layers) & alti_mchoi(*) = !values.f_nan

zsur = 0.0
alti_mchoi(0) = zsur / 1000.0
ccon = - 9.81 * 28.9 / 8314.0 * 500.0


;fix this equation
for n = 1L, n_layers-1 do begin
    avit = (1.0/temp_mchoi(n))+(1.0/temp_mchoi(n-1))
    alti_mchoi(n) = alti_mchoi(n-1) - alog(pres_mchoi(n-1)/pres_mchoi(n))/avit/ccon
endfor

print, pres_mchoi(32), alti_mchoi(32)
stop

H2O_mchoi  = interpol(H2O_daq, pres_daq, pres_mchoi)
O3_mchoi   = interpol(O3_daq, pres_daq, pres_mchoi)
NO2_mchoi  = interpol(NO2_daq, pres_daq, pres_mchoi)
HCHO_mchoi = interpol(HCHO_daq, pres_daq, pres_mchoi)
SO2_mchoi  = interpol(SO2_daq, pres_daq, pres_mchoi)
CO_mchoi   = interpol(CO_daq, pres_daq, pres_mchoi)


;--save Pressure-Tempereture-Gas VMR
openw, 9, 'D:\GOCI_OP\LUT_vlidort\Profile_PTH_US76_GAS_discoveraq_mchoi.dat'

for i = 0L, n_elements(pres_mchoi)-1 do begin

	printf, 9, i+1, pres_mchoi(i), temp_mchoi(i), H2O_mchoi(i), O3_mchoi(i), NO2_mchoi(i), HCHO_mchoi(i), SO2_mchoi(i), CO_mchoi(i), $
	          f='(1(i3,3x), (1(f12.7),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x), (1(E11.4),x))'

endfor

close, 9

;--save layer-Altitude corresponding with layer-Pressure

openw, 9, 'D:\GOCI_OP\LUT_vlidort\UVN_heightgrid_0km_US76_fin_mchoi.dat'

for i = n_elements(pres_mchoi)-1, 0, -1 do begin

	printf, 9, 0, alti_mchoi(i), $
	          f='(1(i12,3x), (1(f7.4)))'

endfor

close, 9

stop


end