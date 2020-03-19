pro gauss_asyconv2, solref_wl,solrefspec,ap_lambda,wave_range,convol_spec,n_sample, FWHM, asym, con_rag,c_shift,c_squeeze

hw1e = FWHM / (2*sqrt(alog(2)))
avg_lambda=(wave_range[1]+wave_range[0])/2

ccd_lambda =  avg_lambda + c_shift + (ap_lambda - avg_lambda) * (1.0 + c_squeeze )

convol_spec = dblarr(n_sample)  ;convolved measured spectrum array

for i_lambda = 0d, n_sample-1 do begin
   index=where((solref_wl gt (ccd_lambda[i_lambda]-con_rag)) and $
               (solref_wl lt (ccd_lambda[i_lambda]+con_rag)))
   n_index = size(index,/n_elements)

   tmpsum = 0.0d
   for i_index = 0d, n_index-1 do begin
      d_lambda = ccd_lambda[i_lambda] - solref_wl[index[i_index]]
      gauss_ft=exp(-(d_lambda)^2.0/ (hw1e * ( 1.0 + sign(d_lamda) * asym)) ^2.0)
      tmpsum = tmpsum + gauss_ft
      convol_spec[i_lambda] = convol_spec[i_lambda] + $
                              solrefspec[index[i_index]] * gauss_ft 
      
   endfor
   convol_spec[i_lambda] = convol_spec[i_lambda] / tmpsum
   
endfor

end

