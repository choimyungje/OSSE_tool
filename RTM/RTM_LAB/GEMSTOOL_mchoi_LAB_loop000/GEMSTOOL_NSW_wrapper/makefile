# @@@@@@@@@@@@@@@ PATH VARIABLES

#  NSW-specific MAKEFILE. 12 August 2013

#  upgraded to VLIDORT Version 2.7, 18 October 2016.
#  upgrades for Fluorescence and BRDFs, 18-25 October 2016.

#  Parent path

GEMSTOOL_PATH = ..

# Source code paths. VLIDORT Upgraded 10/17/16.

VLIDORT_V2p7_PATH   = $(GEMSTOOL_PATH)/vlidort_2p7_sourcecode
VLIDORT_STRC_PATH   = $(VLIDORT_V2p7_PATH)/vlidort_def
VLIDORT_MAIN_PATH   = $(VLIDORT_V2p7_PATH)/vlidort_main_FO_1p4
FO_MAIN_PATH        = $(VLIDORT_V2p7_PATH)/fo_main_1p4
VSUP_PATH           = $(VLIDORT_V2p7_PATH)/vsup

GEMSTOOL_MAIN_PATH  = $(GEMSTOOL_PATH)/GEMSTOOL_sourcecode
RTSMIE_MAIN_PATH    = $(GEMSTOOL_PATH)/RTSMie_sourcecode
TMATRIX_MAIN_PATH   = $(GEMSTOOL_PATH)/Tmatrix_sourcecode
XSECHITRAN_PATH     = $(GEMSTOOL_PATH)/XsecHitran_sourcecode
XSECLUT_PATH        = $(GEMSTOOL_PATH)/XsecLut_sourcecode

# # Source code paths for LRRS V2p3 mchoi 03/25/19
LRRS_DEF_PATH  = $(GEMSTOOL_PATH)/RTS_LRRS2p3_sourcecode/lrrs_def
LRRS_MAIN_PATH = $(GEMSTOOL_PATH)/RTS_LRRS2p3_sourcecode/lrrs_main
# # Source code path for Sioris Raman code mchoi 03/25/19
# SIORIS_PATH = $(GEMSTOOL_PATH)/SiorisRaman_sourcecode

# # Convolution code path
CONV_PATH = $(GEMSTOOL_PATH)/Convolution_sourcecode


#  Wrapper path (perhaps do not need this). Upgraded 10/17/16.
GEMSTOOL_TEST_PATH   = $(GEMSTOOL_PATH)/GEMSTOOL_NSW_wrapper

#  Modulus and Object files
MOD_PATH = $(GEMSTOOL_PATH)/Modules_NSW
OBJ_PATH = $(GEMSTOOL_PATH)/Objects_NSW
MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o

# @@@@@@@@@@@@@@@ Shell

SHELL = /bin/bash

# @@@@@@@@@@@@@@@ Define FORTRAN90 compiler to use, can be defined here
#   -- Alternative, defined on makefile command line

#  Intel
#FC = ifort
#FC = /mnt/home/hiren.jethva/intel/fc/9.1.043/bin/ifort

#  gfortran
FC = gfortran

#  g95
#FC = g95
#FC = /mnt/home/hiren.jethva/Download/g95-install/bin/i686-unknown-linux-gnu-g95

#  NAG
#FC = f95

# pgf90
#FC = pgf90

# @@@@@@@@@@@@@@@ Define FORTRAN90 compiler flags

# Additional flags for Intel
ifeq ($(FC), ifort)
	FFLAGS = -c -warn all -check all
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -module $(MOD_PATH)
#	FFLAGS_DEBUG = -g -warn all -check all -traceback
endif

# Additional flags for gfortran
ifeq ($(FC), gfortran)
	FFLAGS = -c -mcmodel=large
#	FFLAGS = -c -mcmodel=medium
#	FFLAGS = -c -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
#	FFLAGS = -c -Wall -fbounds-check -frange-check -mcmodel=medium
#	FFLAGS = -c -O3 -mcmodel=medium
# Rob setting
#	FFLAGS = -c -fbounds-check -frange-check -mcmodel=medium
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -J$(MOD_PATH)
	FFLAGS_DEBUG = -Wall -fbounds-check -fbacktrace
endif

# Additional flags for g95
ifeq ($(FC), g95)
	FFLAGS = -c
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
#      older g95
#	FFLAGS_DEBUG = -g -Wall -fbounds-check
#      g95 v0.92
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -ftrace=full
endif

# Additional flags for NAG
#ifeq ($(FC), f95)
ifeq ($(FC), $(filter $(FC), f95 nagfor))
	FFLAGS = -c
	FFLAGS := $(FFLAGS) -mdir $(MOD_PATH) -I$(MOD_PATH)
#	FFLAGS = -c -w=obs -w=unused -C=array -gline
#	FFLAGS_DEBUG = -w=obs -w=unused -C=array -gline
#	FFLAGS_DEBUG = -w=obs -w=unused -C=array -C=undefined -gline
	FFLAGS_DEBUG = -w=obs -w=unused -C=all -C=undefined -gline
#	FFLAGS_DEBUG = -g -C=all -C=undefined -gline -mtrace=all -nan
#	FFLAGS_DEBUG = -w=obs -w=unused -g -C=all -C=undefined -gline -mtrace=all -nan
endif

# Additional flags for pgf90
ifeq ($(FC), pgf90)
        FFLAGS = -c
        FFLAGS := $(FFLAGS) -I$(MOD_PATH) -module $(MOD_PATH)
#       FFLAGS_DEBUG = -g -warn all -check all -traceback
endif

# For debug build, use "make DEBUG=t"
ifeq ($(DEBUG), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_DEBUG)
endif


# Define list of source files
# (Note: ordering is important because of dependencies)

V2p7_SOURCES =
V2p7_SOURCES +=   \
   $(VLIDORT_STRC_PATH)/vlidort_pars.f90

V2p7_SOURCES +=   \
   $(VLIDORT_MAIN_PATH)/lapack_tools.f90              \
   $(VLIDORT_STRC_PATH)/vlidort_inputs_def.f90        \
   $(VLIDORT_STRC_PATH)/vlidort_sup_brdf_def.f90      \
   $(VLIDORT_STRC_PATH)/vlidort_sup_ss_def.f90        \
   $(VLIDORT_STRC_PATH)/vlidort_sup_sleave_def.f90    \
   $(VLIDORT_STRC_PATH)/vlidort_sup_def.f90           \
   $(VLIDORT_STRC_PATH)/vlidort_outputs_def.f90	      \
   $(VLIDORT_STRC_PATH)/vlidort_io_defs.f90	      \
   $(VLIDORT_STRC_PATH)/vlidort_work_def.f90   	      \
   $(VLIDORT_MAIN_PATH)/vlidort_aux.f90               \
   $(VLIDORT_MAIN_PATH)/vlidort_geometry.f90          \
   $(VLIDORT_MAIN_PATH)/vlidort_Taylor.f90            \
   $(VLIDORT_MAIN_PATH)/vlidort_inputs.f90            \
   $(VLIDORT_MAIN_PATH)/vlidort_miscsetups.f90        \
   $(VLIDORT_MAIN_PATH)/vlidort_multipliers.f90       \
   $(VLIDORT_MAIN_PATH)/vlidort_corrections.f90       \
   $(VLIDORT_MAIN_PATH)/vlidort_thermalsup.f90        \
   $(VLIDORT_MAIN_PATH)/vlidort_solutions.f90         \
   $(VLIDORT_MAIN_PATH)/vlidort_bvproblem.f90         \
   $(VLIDORT_MAIN_PATH)/vlidort_intensity.f90         \
   $(VLIDORT_MAIN_PATH)/vlidort_pack.f90	      \
   $(VLIDORT_MAIN_PATH)/vlidort_unpack.f90	      \
   $(VLIDORT_MAIN_PATH)/vlidort_writemodules.f90      \
   $(VLIDORT_MAIN_PATH)/vlidort_masters.f90

V2p7_L_SOURCES =
V2p7_L_SOURCES += \
   $(VLIDORT_STRC_PATH)/vlidort_lin_inputs_def.f90    \
   $(VLIDORT_STRC_PATH)/vlidort_lin_sup_brdf_def.f90  \
   $(VLIDORT_STRC_PATH)/vlidort_lin_sup_ss_def.f90    \
   $(VLIDORT_STRC_PATH)/vlidort_lin_sup_sleave_def.f90\
   $(VLIDORT_STRC_PATH)/vlidort_lin_sup_def.f90       \
   $(VLIDORT_STRC_PATH)/vlidort_lin_outputs_def.f90   \
   $(VLIDORT_STRC_PATH)/vlidort_lin_io_defs.f90       \
   $(VLIDORT_STRC_PATH)/vlidort_lin_work_def.f90      \
   $(VLIDORT_MAIN_PATH)/vlidort_l_inputs.f90         \
   $(VLIDORT_MAIN_PATH)/vlidort_la_miscsetups.f90    \
   $(VLIDORT_MAIN_PATH)/vlidort_la_corrections.f90   \
   $(VLIDORT_MAIN_PATH)/vlidort_ls_corrections.f90   \
   $(VLIDORT_MAIN_PATH)/vlidort_l_thermalsup.f90     \
   $(VLIDORT_MAIN_PATH)/vlidort_lpc_solutions.f90    \
   $(VLIDORT_MAIN_PATH)/vlidort_lpc_bvproblem.f90    \
   $(VLIDORT_MAIN_PATH)/vlidort_lbbf_jacobians.f90   \
   $(VLIDORT_MAIN_PATH)/vlidort_ls_wfsurface.f90     \
   $(VLIDORT_MAIN_PATH)/vlidort_ls_wfsleave.f90      \
   $(VLIDORT_MAIN_PATH)/vlidort_l_writemodules.f90   \
   $(VLIDORT_MAIN_PATH)/vlidort_l_pack.f90	     \
   $(VLIDORT_MAIN_PATH)/vlidort_l_unpack.f90

V2p7_LPS_SOURCES =
V2p7_LPS_SOURCES += \
   $(VLIDORT_MAIN_PATH)/vlidort_lp_miscsetups.f90    \
   $(VLIDORT_MAIN_PATH)/vlidort_lp_corrections.f90   \
   $(VLIDORT_MAIN_PATH)/vlidort_lp_solutions.f90     \
   $(VLIDORT_MAIN_PATH)/vlidort_lp_bvproblem.f90     \
   $(VLIDORT_MAIN_PATH)/vlidort_lp_wfatmos.f90       \
   $(VLIDORT_MAIN_PATH)/vlidort_lp_pack.f90	     \
   $(VLIDORT_MAIN_PATH)/vlidort_lp_unpack.f90	     \
   $(VLIDORT_MAIN_PATH)/vlidort_lps_masters.f90

V2p7_LCS_SOURCES =
V2p7_LCS_SOURCES += \
   $(VLIDORT_MAIN_PATH)/vlidort_lc_miscsetups.f90    \
   $(VLIDORT_MAIN_PATH)/vlidort_lc_corrections.f90   \
   $(VLIDORT_MAIN_PATH)/vlidort_lc_solutions.f90     \
   $(VLIDORT_MAIN_PATH)/vlidort_lc_bvproblem.f90     \
   $(VLIDORT_MAIN_PATH)/vlidort_lc_wfatmos.f90       \
   $(VLIDORT_MAIN_PATH)/vlidort_lc_pack.f90	     \
   $(VLIDORT_MAIN_PATH)/vlidort_lc_unpack.f90	     \
   $(VLIDORT_MAIN_PATH)/vlidort_lcs_masters.f90

V2p7_SOURCES_PLUS =
V2p7_SOURCES_PLUS +=   \
   $(V2p7_SOURCES)     \
   $(V2p7_L_SOURCES)   \
   $(V2p7_LCS_SOURCES) \
   $(V2p7_LPS_SOURCES)

# (Include vector supplement source files). Upgrade 10/17/16
include $(VSUP_PATH)/makefile.vsup

# (Include first-order source files). Upgrade 10/17/16
include $(FO_MAIN_PATH)/makefile.fo

#  First-Order sources. Retired, 10/17/16.
#FO_SOURCES   = 
#FO_SOURCES +=   \
#   $(FO_PATH_Generics)/FO_geometry_pool_01Jun12.f90     \
#   $(FO_PATH_Generics)/FO_geometry_SSonly_01Jun12.f90   \
#   $(FO_PATH_Generics)/FO_Vector_Spherfuncs.f90         \
#   $(FO_PATH_Generics)/FO_Vector_RTCalcs.f90            \
#   $(FO_PATH_Masters)/FO_Vector_SS_Master_Mk2.f90

# HITRAN CRS source

CRS_SOURCES   = 
CRS_SOURCES +=   \
   $(XSECHITRAN_PATH)/Linespec_aux.f90          \
   $(XSECHITRAN_PATH)/Linespec_basic.f90        \
   $(XSECHITRAN_PATH)/get_hitran_crs.f90

CRS_SOURCES_PLUS   = 
CRS_SOURCES_PLUS +=   \
   $(XSECHITRAN_PATH)/Linespec_aux.f90          \
   $(XSECHITRAN_PATH)/Linespec_basic.f90        \
   $(XSECHITRAN_PATH)/Linespec_plus.f90         \
   $(XSECHITRAN_PATH)/get_hitran_crs.f90        \
   $(XSECHITRAN_PATH)/get_hitran_crs_plus.f90

# LUT CRS source

CRS2_SOURCES   =
CRS2_SOURCES +=   \
	$(XSECLUT_PATH)/get_lut_crs.f90

CRS2_SOURCES_PLUS   = 
CRS2_SOURCES_PLUS +=   \
    $(XSECLUT_PATH)/get_lut_crs.f90        \
    $(XSECLUT_PATH)/get_lut_crs_plus.f90

# RTS Mie and Tmatrix sources
# 10/21/16. Bimodal masters retired, new aerosol scheme.

MIE_SOURCES   = 
MIE_SOURCES +=   \
   $(RTSMIE_MAIN_PATH)/RTSMie_parameters.f90          \
   $(RTSMIE_MAIN_PATH)/RTSMie_distribution.f90        \
   $(RTSMIE_MAIN_PATH)/RTSMie_sourcecode_2.f90

#   $(RTSMIE_MAIN_PATH)/RTSMie_master_bimodal.f90

MIE_SOURCES_PLUS   = 
MIE_SOURCES_PLUS +=   \
   $(MIE_SOURCES)     \
   $(RTSMIE_MAIN_PATH)/RTSMie_sourcecode_plus_2.f90

#   $(RTSMIE_MAIN_PATH)/RTSMie_master_bimodal_PLUS.f90

TMAT_SOURCES   = 
TMAT_SOURCES +=   \
   $(TMATRIX_MAIN_PATH)/tmat_parameters.f90     \
   $(TMATRIX_MAIN_PATH)/tmat_distributions.f90  \
   $(TMATRIX_MAIN_PATH)/tmat_functions.f90      \
   $(TMATRIX_MAIN_PATH)/tmat_makers.f90         \
   $(TMATRIX_MAIN_PATH)/tmat_scattering.f90     \
   $(TMATRIX_MAIN_PATH)/tmat_master.f90
#   $(TMATRIX_MAIN_PATH)/tmat_master_bimodal.f90

TMAT_SOURCES_PLUS   = 
TMAT_SOURCES_PLUS +=   \
   $(TMAT_SOURCES)     \
   $(TMATRIX_MAIN_PATH)/tmat_functions_PLUS.f90      \
   $(TMATRIX_MAIN_PATH)/tmat_makers_PLUS.f90         \
   $(TMATRIX_MAIN_PATH)/tmat_scattering_PLUS.f90     \
   $(TMATRIX_MAIN_PATH)/tmat_master_PLUS.f90
#   $(TMATRIX_MAIN_PATH)/tmat_master_bimodal_PLUS.f90

#  GEMSTOOL SOURCECODE
#    NSW-specific PTHgas, Cross-section, IOPS Modules. 12 August 2013
#  Added Fluorescence module, 10/18/16
#  Added BRDFSurface  module, 10/24/16

GEMSTOOL_SOURCES   = 
GEMSTOOL_SOURCES +=   \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_Type_Structures.f90         \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_numerical.f90               \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_Read_Configfiles.f90        \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_RTInitialize.f90            \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_pthg_profiles.f90       \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_Rayleigh.f90                \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_cross_sections.f90      \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_cross_sections_LUT.f90  \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_Albedo.f90                  \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_Fluorescence.f90        \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_BRDFSurface.f90             \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_AerLoading.f90              \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_AerProperties.f90           \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_CldProperties.f90           \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_Regular_iops.f90        \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_Regular_RTCALC.f90

GEMSTOOL_SOURCES_PLUS   = 
GEMSTOOL_SOURCES_PLUS +=   \
   $(GEMSTOOL_SOURCES)                                           \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_cross_sections_PLUS.f90    \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_cross_sections_PLUS_LUT.f90    \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_AerProperties_PLUS.f90         \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_Linearized_iops.f90        \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_Linearized_RTCALC.f90          \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_SetJacobians.f90

#  This set by R. Spurr, upgrade for Stokes Q/U Jacobians 12/19/16

GEMSTOOL_SOURCES_PLUS_ALT   = 
GEMSTOOL_SOURCES_PLUS_ALT +=   \
   $(GEMSTOOL_SOURCES)                                           \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_cross_sections_PLUS.f90    \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_cross_sections_PLUS_LUT.f90    \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_AerProperties_PLUS.f90         \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_Linearized_iops.f90        \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_Linearized_RTCALC.f90          \
   $(GEMSTOOL_MAIN_PATH)/GEMSTOOL_NSW_SetJacobians_IQU.f90



#  LRRS Raman Model. 03/25/19 mchoi. LRRS V2p3 from R. Spurr
SOURCES_LRRS =
SOURCES_LRRS  +=  \
   $(LRRS_DEF_PATH)/lrrs_pars.f90			\
   $(LRRS_DEF_PATH)/lrrs_inputs_def.f90			\
   $(LRRS_DEF_PATH)/lrrs_outputs_def.f90		\
   $(LRRS_MAIN_PATH)/lrrs_aux2.f90			\
   $(LRRS_MAIN_PATH)/lrrs_io_check.f90			\
   $(LRRS_MAIN_PATH)/lrrs_deltamscaling.f90		\
   $(LRRS_MAIN_PATH)/lrrs_raman_spectroscopy.f90	\
   $(LRRS_MAIN_PATH)/lrrs_generate_ramanops.f90		\
   $(LRRS_MAIN_PATH)/lrrs_setup_master.f90		\
   $(LRRS_MAIN_PATH)/lrrs_brdf.f90			\
   $(LRRS_MAIN_PATH)/lrrs_smallnos.f90			\
   $(LRRS_MAIN_PATH)/lrrs_corrections_1.f90		\
   $(LRRS_MAIN_PATH)/lrrs_rtsolutions_1.f90		\
   $(LRRS_MAIN_PATH)/lrrs_bvproblem.f90			\
   $(LRRS_MAIN_PATH)/lrrs_postprocessing_1.f90		\
   $(LRRS_MAIN_PATH)/lrrs_raman_intensity_1.f90		\
   $(LRRS_MAIN_PATH)/lrrs_elastic_intensity_1.f90	\
   $(LRRS_MAIN_PATH)/lrrs_elastic_master_1.f90		\
   $(LRRS_MAIN_PATH)/lrrs_sources_master_1.f90		\
   $(LRRS_MAIN_PATH)/lrrs_fourier_master_1.f90		\
   $(LRRS_MAIN_PATH)/lrrs_converge.f90			\
   $(LRRS_MAIN_PATH)/lrrs_write.f90			\
   $(LRRS_MAIN_PATH)/lrrs_main_master.f90


# Convolution code
SOURCES_CONV = 
SOURCES_CONV  +=  \
   $(CONV_PATH)/AsymGaus_Convolution_Tool.f90   \
   $(CONV_PATH)/OMIILS_Convolution_Tool.f90

#  Main definitions
#     $(FO_SOURCES)         \   Add when needed, Not required for Version 2.7
#     $(TMAT_SOURCES)       \   Add when needed.

# 10/17/16. Add FO sources for Version 2.7
# 10/18/16. Added SLEAVE sources for SIF
# 10/24/16. Added BRDF   sources 

SOURCES = 
SOURCES =  $(FO_SOURCES_Vector) + \
   $(V2p7_SOURCES)      \
   $(VBRDF_SUP_SOURCES)   \
   $(VSLEAVE_SUP_SOURCES) \
   $(CRS_SOURCES)       \
   $(CRS2_SOURCES)      \
   $(MIE_SOURCES)       \
   $(TMAT_SOURCES)      \
   $(GEMSTOOL_SOURCES)  \
   $(GEMSTOOL_TEST_PATH)/GEMSTOOL_NSW_RadMaster_V2p7.f90

SOURCES_PLUS = 
SOURCES_PLUS = $(FO_SOURCES_L_Vector) + \
   $(V2p7_SOURCES_PLUS)      \
   $(VBRDF_SUP_SOURCES)      \
   $(VSLEAVE_SUP_SOURCES)    \
   $(CRS_SOURCES_PLUS)       \
   $(CRS2_SOURCES_PLUS)      \
   $(MIE_SOURCES_PLUS)       \
   $(TMAT_SOURCES_PLUS)      \
   $(GEMSTOOL_SOURCES_PLUS)  \
   $(GEMSTOOL_TEST_PATH)/GEMSTOOL_NSW_LinMaster_V2p7.f90

#  This set by R. Spurr, upgrade for Stokes Q/U Jacobians 12/19/16

SOURCES_PLUS_ALT = 
SOURCES_PLUS_ALT = $(FO_SOURCES_L_Vector) + \
   $(V2p7_SOURCES_PLUS)      \
   $(VBRDF_SUP_SOURCES)      \
   $(VBRDF_LINSUP_SOURCES)   \
   $(VSLEAVE_SUP_SOURCES)    \
   $(CRS_SOURCES_PLUS)       \
   $(CRS2_SOURCES_PLUS)      \
   $(MIE_SOURCES_PLUS)       \
   $(TMAT_SOURCES_PLUS)      \
   $(SOURCES_LRRS)   \
   $(SOURCES_CONV)   \
   $(GEMSTOOL_SOURCES_PLUS_ALT)  \
   $(GEMSTOOL_TEST_PATH)/FastV2p7_LRRS2p3_Master.f90 \
   $(GEMSTOOL_TEST_PATH)/GEMSTOOL_NSW_LinIQUMaster_V2p7.f90
#  $(SOURCES_LRRS)           
#  $(SOURCES_SIORIS)
## VBRDF_LINSUP_SOURCES is added by mchoi 05/29/19      
   

# Define pattern rules for creating object files:

.SUFFIXES:

$(OBJ_PATH)/%.o : $(VLIDORT_STRC_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VLIDORT_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(GEMSTOOL_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(XSECHITRAN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(XSECLUT_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(RTSMIE_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(TMATRIX_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(GEMSTOOL_TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# #  LRRS Version 2p3 sourcecode
$(OBJ_PATH)/%.o : $(LRRS_DEF_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(LRRS_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# # Convolution sourcecode
$(OBJ_PATH)/%.o : $(CONV_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@



# #  LRRS Version 2p5a sourcecode
# $(OBJ_PATH)/%.o : $(LRRS_DEF_PATH)/%.f90
# 	$(FC) $(FFLAGS) $< -o $@
# $(OBJ_PATH)/%.o : $(LRRS_MAIN_PATH)/%.f90
# 	$(FC) $(FFLAGS) $< -o $@
# $(OBJ_PATH)/%.o : $(LRRS_SUP_PATH)/%.f90
# 	$(FC) $(FFLAGS) $< -o $@
# $(OBJ_PATH)/%.o : $(LRRS_TEST_PATH)/%.f90
# 	$(FC) $(FFLAGS) $< -o $@
# $(OBJ_PATH)/%.o : $(LRRS_SETUP_PATH)/%.f90
# 	$(FC) $(FFLAGS) $< -o $@   
# # For utility source files
# $(OBJ_PATH)/%.o : $(UTIL_PATH)/%.f90
# 	$(FC) $(FFLAGS) $< -o $@


# #  Sioris Raman sourcecode
# $(OBJ_PATH)/%.o : $(SIORIS_PATH)/%.f90
# 	$(FC) $(FFLAGS) $< -o $@

#  Add these later.........

#$(OBJ_PATH)/%.o : $(FO_PATH_Generics)/%.f90
#	$(FC) $(FFLAGS) $< -o $@
#$(OBJ_PATH)/%.o : $(FO_PATH_Masters)/%.f90
#	$(FC) $(FFLAGS) $< -o $@

# Define object files

# F90SOURCES   := $(notdir $(filter %.f90, $(SOURCES)))
# F90OBJECTS   := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES)))

# F90SOURCES_PLUS   := $(notdir $(filter %.f90, $(SOURCES_PLUS)))
# F90OBJECTS_PLUS   := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_PLUS)))

#  This set by R. Spurr, upgrade for Stokes Q/U Jacobians 12/19/16

F90SOURCES_PLUS_ALT  := $(notdir $(filter %.f90, $(SOURCES_PLUS_ALT)))
F90OBJECTS_PLUS_ALT   := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_PLUS_ALT)))

# Define desired targets

# ALL: Gemstool_NSW_Regular.exe Gemstool_NSW_Linearized.exe Gemstool_NSW_IQULinearized.exe
ALL: Gemstool_NSW_IQULinearized.exe

# Gemstool_NSW_Regular.exe: $(F90OBJECTS)
# 	$(FC) $^ -o $@

# Gemstool_NSW_Linearized.exe: $(F90OBJECTS_PLUS)
# 	$(FC) $^ -o $@

Gemstool_NSW_IQULinearized.exe: $(F90OBJECTS_PLUS_ALT)
	$(FC) $^ -o $@

#  Make clean

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.LOG *.exe
