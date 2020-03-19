rm *.o *.mod


gfortran -c -Wall ../../lrrs_def/lrrs_pars.f90
gfortran -c -Wall ../../lrrs_def/lrrs_inputs_def.f90
gfortran -c -Wall ../../lrrs_def/lrrs_outputs_def.f90
gfortran -c -Wall ../../lrrs_def/lrrs_sup_brdf_def.f90
gfortran -c -Wall ../../lrrs_def/lrrs_sup_sleave_def.f90
gfortran -c -Wall ../../lrrs_def/lrrs_sup_ss_def.f90
gfortran -c -Wall ../../lrrs_def/lrrs_sup_def.f90
gfortran -c -Wall ../../lrrs_def/lrrs_io_defs.f90
gfortran -c -Wall ../../lrrs_L_def/lrrs_lin_inputs_def.f90
gfortran -c -Wall ../../lrrs_L_def/lrrs_lin_outputs_def.f90
gfortran -c -Wall ../../lrrs_L_def/lrrs_lin_sup_brdf_def.f90
gfortran -c -Wall ../../lrrs_L_def/lrrs_lin_sup_sleave_def.f90
gfortran -c -Wall ../../lrrs_L_def/lrrs_lin_sup_ss_def.f90
gfortran -c -Wall ../../lrrs_L_def/lrrs_lin_sup_def.f90
gfortran -c -Wall ../../lrrs_L_def/lrrs_lin_io_defs.f90
gfortran -c -Wall ../../lrrs_main/lrrs_aux2.f90
gfortran -c -Wall ../../lrrs_main/lrrs_raman_spectroscopy.f90
gfortran -c -Wall ../lrrs_Phasfunc_Sup_Accessory.f90
gfortran -c -Wall ../lrrs_Phasfunc_LinSup_Accessory.f90
gfortran -c -Wall ../../lrrs_sup/brdf/brdf_sup_inputs_def.f90  
gfortran -c -Wall ../../lrrs_sup/brdf/brdf_sup_outputs_def.f90 
gfortran -c -Wall ../../lrrs_sup/brdf/brdf_findpar.f90
gfortran -c -Wall ../../lrrs_sup/brdf/brdf_sup_aux.f90  
gfortran -c -Wall ../../lrrs_sup/brdf/brdf_sup_kernels.f90 
gfortran -c -Wall ../../lrrs_sup/brdf/brdf_sup_routines.f90
gfortran -c -Wall ../../lrrs_sup/brdf/brdf_sup_masters.f90 
gfortran -c -Wall ../../lrrs_sup/brdf/brdf_sup_mod.f90 
gfortran -c -Wall lrrs_iopsetup_Aux.f90
gfortran -c -Wall lrrs_iopsetup_WorkRoutines.f90
gfortran -c -Wall lrrs_LPbin_M_iopsetup.f90
exit

gfortran -c -Wall ../lrrs_def/lrrs_inputs_def.f90
gfortran -c -Wall ../lrrs_main/lrrs_geometry.f90
gfortran -c -Wall ../lrrs_def/lrrs_outputs_def.f90
gfortran -c -Wall ../lrrs_main/lrrs_aux2.f90
gfortran -c -Wall ../lrrs_main/lrrs_deltamscaling.f90
gfortran -c -Wall ../lrrs_main/lrrs_raman_spectroscopy.f90
gfortran -c -Wall ../lrrs_main/lrrs_generate_ramanops.f90
gfortran -c -Wall lrrs_L_raman_spectroscopy.f90
gfortran -c -Wall lrrs_L_generate_ramanops.f90
gfortran -c -Wall lrrs_L_setup_master.f90
