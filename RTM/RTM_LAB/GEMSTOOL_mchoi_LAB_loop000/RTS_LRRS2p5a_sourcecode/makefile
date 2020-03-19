#
# This file serves as the head Makefile for LIDORT-RRS and its "tester" programs
# UPGRADE TO LRRS 2.5, September 2015 to February 2017
#

#  Compilation with OMP flags. 10.10.18
#      o3bin REENGINEERED test, done September 2017
#      LPbin REENGINEERED test, new  October   2018

# Definition of Variables
# =======================

# Sourcecode paths

LRRS_DEF_PATH    = lrrs_def
LRRS_MAIN_PATH   = lrrs_main
LRRS_SUP_PATH    = lrrs_sup
LRRS_L_DEF_PATH  = lrrs_L_def
LRRS_L_MAIN_PATH = lrrs_L_main
LRRS_L_SUP_PATH  = lrrs_L_sup

# Test drectory and IOPsetup paths

LRRS_TEST_PATH   = lrrs_test
LRRS_SETUP_PATH  = $(LRRS_TEST_PATH)/iopsetup

# Utilities path

UTIL_PATH        = util

# Modules and Objects

MOD_PATH = mod
OBJ_PATH = obj
MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o

# Define default shell make will use
# ==================================

SHELL = /bin/bash

# Define FORTRAN90 compiler to use
# ================================

#    (can be defined here, but usually defined on makefile command line)

#  Intel
#FC = ifort
#  gfortran
FC = gfortran
#  g95
#FC = g95
#  NAG
#FC = f95
#FC = nagfor

# Define FORTRAN90 compiler flags
# ###############################

FFLAGS = -c

# Additional flags for Intel
ifeq ($(FC), ifort)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -module $(MOD_PATH)
	FFLAGS_DEBUG = -g -warn all -check all -traceback
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -openmp
endif

# Additional flags for gfortran
ifeq ($(FC), gfortran)
	FFLAGS := $(FFLAGS) -mcmodel=large -I$(MOD_PATH) -J$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -fbacktrace
#	FFLAGS_DEBUG = -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
#	FFLAGS_DEBUG = -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow -fbacktrace
#	FFLAGS_DEBUG = -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-local-zero
	FFLAGS_DEBUG = -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-integer=-100 -finit-real=nan 
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -fopenmp -frecursive
endif

# Additional flags for NAG
#ifeq ($(FC), f95)
ifeq ($(FC), $(filter $(FC), f95 nagfor))
	FFLAGS := $(FFLAGS) -mdir $(MOD_PATH) -I$(MOD_PATH)
#	FFLAGS_DEBUG = -g -C=all -C=undefined -gline -mtrace=all -nan
#	FFLAGS_DEBUG = -w=obs -w=unused -g -C=all -C=undefined -gline -mtrace=all -nan
#	FFLAGS_DEBUG = -w=obs -w=unused -C=all -C=undefined -gline -mtrace=all
	FFLAGS_DEBUG = -w=obs -w=unused -C=all -C=undefined -gline
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -openmp
endif

# For debug build, use "make DEBUG=t"
ifeq ($(DEBUG), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_DEBUG)
endif

# For optimized build, use "make OPT=t"
ifeq ($(OPT), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPT)
endif

# For parallel build using OpenMP, use "make OPENMP=t"
ifeq ($(OPENMP), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
endif

# Define list of source files
# ===========================

# (Note: ordering is important because of dependencies)

BASE_SOURCES =
SOURCES =
L_SOURCES =

BASE_SOURCES +=   \
   $(LRRS_DEF_PATH)/lrrs_pars.f90

# For RTS debug build
ifeq ($(RTS_DEBUG), t)
   BASE_SOURCES +=   \
      $(UTIL_PATH)/rts_debug.f90
endif

ifeq ($(OPENMP), t)
   SOURCES +=   \
      $(BASE_SOURCES) \
      $(LRRS_DEF_PATH)/lrrs_inputs_def.f90		\
      $(LRRS_DEF_PATH)/lrrs_sup_brdf_def.f90		\
      $(LRRS_DEF_PATH)/lrrs_sup_sleave_def.f90		\
      $(LRRS_DEF_PATH)/lrrs_sup_ss_def.f90		\
      $(LRRS_DEF_PATH)/lrrs_sup_def.f90			\
      $(LRRS_DEF_PATH)/lrrs_outputs_def.f90		\
      $(LRRS_DEF_PATH)/lrrs_io_defs.f90			\
      $(LRRS_MAIN_PATH)/lrrs_aux2.f90			\
      $(LRRS_MAIN_PATH)/lrrs_io_check.f90		\
      $(LRRS_MAIN_PATH)/lrrs_geometry.f90		\
      $(LRRS_MAIN_PATH)/lrrs_Taylor.f90			\
      $(LRRS_MAIN_PATH)/lrrs_deltamscaling.f90		\
      $(LRRS_MAIN_PATH)/lrrs_raman_spectroscopy.f90	\
      $(LRRS_MAIN_PATH)/lrrs_generate_ramanops.f90	\
      $(LRRS_MAIN_PATH)/lrrs_rrsoptical_master.f90	\
      $(LRRS_MAIN_PATH)/lrrs_postprocessing_1.f90	\
      $(LRRS_MAIN_PATH)/lrrs_corrections_1.f90		\
      $(LRRS_MAIN_PATH)/lrrs_miscsetups_1.f90		\
      $(LRRS_MAIN_PATH)/lrrs_rtsolutions_1.f90		\
      $(LRRS_MAIN_PATH)/lrrs_bvproblem.f90		\
      $(LRRS_MAIN_PATH)/lrrs_raman_intensity_1.f90	\
      $(LRRS_MAIN_PATH)/lrrs_elastic_intensity_1.f90	\
      $(LRRS_MAIN_PATH)/lrrs_elastic_master_1_REENG.f90	\
      $(LRRS_MAIN_PATH)/lrrs_sources_master_1.f90	\
      $(LRRS_MAIN_PATH)/lrrs_fourier_master_1_REENG.f90	\
      $(LRRS_MAIN_PATH)/lrrs_converge.f90		\
      $(LRRS_MAIN_PATH)/lrrs_writemodules.f90		\
      $(LRRS_MAIN_PATH)/lrrs_main_master_REENG.f90
else
   SOURCES +=   \
      $(BASE_SOURCES) \
      $(LRRS_DEF_PATH)/lrrs_inputs_def.f90		\
      $(LRRS_DEF_PATH)/lrrs_sup_brdf_def.f90		\
      $(LRRS_DEF_PATH)/lrrs_sup_sleave_def.f90		\
      $(LRRS_DEF_PATH)/lrrs_sup_ss_def.f90		\
      $(LRRS_DEF_PATH)/lrrs_sup_def.f90			\
      $(LRRS_DEF_PATH)/lrrs_outputs_def.f90		\
      $(LRRS_DEF_PATH)/lrrs_io_defs.f90			\
      $(LRRS_MAIN_PATH)/lrrs_aux2.f90			\
      $(LRRS_MAIN_PATH)/lrrs_io_check.f90		\
      $(LRRS_MAIN_PATH)/lrrs_geometry.f90		\
      $(LRRS_MAIN_PATH)/lrrs_Taylor.f90			\
      $(LRRS_MAIN_PATH)/lrrs_deltamscaling.f90		\
      $(LRRS_MAIN_PATH)/lrrs_raman_spectroscopy.f90	\
      $(LRRS_MAIN_PATH)/lrrs_generate_ramanops.f90	\
      $(LRRS_MAIN_PATH)/lrrs_rrsoptical_master.f90	\
      $(LRRS_MAIN_PATH)/lrrs_postprocessing_1.f90	\
      $(LRRS_MAIN_PATH)/lrrs_corrections_1.f90		\
      $(LRRS_MAIN_PATH)/lrrs_miscsetups_1.f90		\
      $(LRRS_MAIN_PATH)/lrrs_rtsolutions_1.f90		\
      $(LRRS_MAIN_PATH)/lrrs_bvproblem.f90		\
      $(LRRS_MAIN_PATH)/lrrs_raman_intensity_1.f90	\
      $(LRRS_MAIN_PATH)/lrrs_elastic_intensity_1.f90	\
      $(LRRS_MAIN_PATH)/lrrs_elastic_master_1.f90	\
      $(LRRS_MAIN_PATH)/lrrs_sources_master_1.f90	\
      $(LRRS_MAIN_PATH)/lrrs_fourier_master_1.f90	\
      $(LRRS_MAIN_PATH)/lrrs_converge.f90		\
      $(LRRS_MAIN_PATH)/lrrs_writemodules.f90		\
      $(LRRS_MAIN_PATH)/lrrs_main_master.f90
endif

#  add linearized sources

ifeq ($(OPENMP), t)
L_SOURCES += \
   $(SOURCES) \
   $(LRRS_L_DEF_PATH)/lrrs_lin_inputs_def.f90		\
   $(LRRS_L_DEF_PATH)/lrrs_lin_sup_brdf_def.f90		\
   $(LRRS_L_DEF_PATH)/lrrs_lin_sup_sleave_def.f90	\
   $(LRRS_L_DEF_PATH)/lrrs_lin_sup_ss_def.f90		\
   $(LRRS_L_DEF_PATH)/lrrs_lin_sup_def.f90		\
   $(LRRS_L_DEF_PATH)/lrrs_lin_outputs_def.f90		\
   $(LRRS_L_DEF_PATH)/lrrs_lin_io_defs.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_io_check.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_raman_spectroscopy.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_generate_ramanops.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_rrsoptical_master.f90     \
   $(LRRS_L_MAIN_PATH)/lrrs_L_postprocessing_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_corrections_plus_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_miscsetups_1.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_rtsolutions_1.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_bvproblem.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_sleave.f90                \
   $(LRRS_L_MAIN_PATH)/lrrs_L_raman_jacobians_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_elastic_jacobians_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_elastic_master_1_REENG.f90\
   $(LRRS_L_MAIN_PATH)/lrrs_L_sources_master_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_fourier_master_1_REENG.f90\
   $(LRRS_L_MAIN_PATH)/lrrs_L_converge.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_writemodules.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_main_master_REENG.f90
else
L_SOURCES += \
   $(SOURCES) \
   $(LRRS_L_DEF_PATH)/lrrs_lin_inputs_def.f90		\
   $(LRRS_L_DEF_PATH)/lrrs_lin_sup_brdf_def.f90		\
   $(LRRS_L_DEF_PATH)/lrrs_lin_sup_sleave_def.f90	\
   $(LRRS_L_DEF_PATH)/lrrs_lin_sup_ss_def.f90		\
   $(LRRS_L_DEF_PATH)/lrrs_lin_sup_def.f90		\
   $(LRRS_L_DEF_PATH)/lrrs_lin_outputs_def.f90		\
   $(LRRS_L_DEF_PATH)/lrrs_lin_io_defs.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_io_check.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_raman_spectroscopy.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_generate_ramanops.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_rrsoptical_master.f90     \
   $(LRRS_L_MAIN_PATH)/lrrs_L_postprocessing_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_corrections_plus_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_miscsetups_1.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_rtsolutions_1.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_bvproblem.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_sleave.f90                \
   $(LRRS_L_MAIN_PATH)/lrrs_L_raman_jacobians_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_elastic_jacobians_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_elastic_master_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_sources_master_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_fourier_master_1.f90	\
   $(LRRS_L_MAIN_PATH)/lrrs_L_converge.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_writemodules.f90		\
   $(LRRS_L_MAIN_PATH)/lrrs_L_main_master.f90
endif

# (Include supplement source files)
include $(LRRS_SUP_PATH)/makefile.sup
include $(LRRS_L_SUP_PATH)/makefile.L_sup


# Main tests
# ----------

# Bin tests (Use bin PARS file)
#   -- Rob mod 5/12/17 for 2p5a, Accessory module for Phase functions

ifeq ($(OPENMP), t)
   SOURCES_BIN = $(SOURCES) +  \
      $(BRDF_SUP_SOURCES)                         	\
      $(SLEAVE_SUP_SOURCES)                       	\
      $(LRRS_SETUP_PATH)/lrrs_iopsetup_Aux.f90          \
      $(LRRS_SETUP_PATH)/lrrs_iopsetup_WorkRoutines.f90 \
      $(LRRS_TEST_PATH)/lrrs_Phasfunc_Sup_Accessory.f90 \
      $(LRRS_SETUP_PATH)/lrrs_o3bin_M_iopsetup.f90	\
      $(LRRS_TEST_PATH)/lrrs_2p5_sup_accessories_1.f90  \
      $(LRRS_TEST_PATH)/lrrs_o3bin_M_REENG_tester.f90
else
   SOURCES_BIN = $(SOURCES) +  \
      $(BRDF_SUP_SOURCES)                         	\
      $(SLEAVE_SUP_SOURCES)                       	\
      $(LRRS_SETUP_PATH)/lrrs_iopsetup_Aux.f90          \
      $(LRRS_SETUP_PATH)/lrrs_iopsetup_WorkRoutines.f90 \
      $(LRRS_TEST_PATH)/lrrs_Phasfunc_Sup_Accessory.f90 \
      $(LRRS_SETUP_PATH)/lrrs_o3bin_M_iopsetup.f90	\
      $(LRRS_TEST_PATH)/lrrs_2p5_sup_accessories_1.f90 	\
      $(LRRS_TEST_PATH)/lrrs_o3bin_M_tester.f90
endif

# LPBin tests (Use bin PARS file)
#   -- Rob mod 5/12/17 for 2p5a, Accessory module for Phase functions

ifeq ($(OPENMP), t)
  SOURCES_BIN_LP = $(L_SOURCES) +  \
     $(BRDF_LINSUP_SOURCES)                         	  \
     $(SLEAVE_LINSUP_SOURCES)                       	  \
     $(LRRS_SETUP_PATH)/lrrs_iopsetup_Aux.f90             \
     $(LRRS_SETUP_PATH)/lrrs_iopsetup_WorkRoutines.f90    \
     $(LRRS_TEST_PATH)/lrrs_Phasfunc_Sup_Accessory.f90    \
     $(LRRS_TEST_PATH)/lrrs_Phasfunc_LinSup_Accessory.f90 \
     $(LRRS_SETUP_PATH)/lrrs_LPbin_M_iopsetup.f90         \
     $(LRRS_TEST_PATH)/lrrs_2p5_sup_accessories_1.f90     \
     $(LRRS_TEST_PATH)/lrrs_LPbin_M_REENG_tester.f90
else
  SOURCES_BIN_LP = $(L_SOURCES) +  \
     $(BRDF_LINSUP_SOURCES)                         	  \
     $(SLEAVE_LINSUP_SOURCES)                       	  \
     $(LRRS_SETUP_PATH)/lrrs_iopsetup_Aux.f90             \
     $(LRRS_SETUP_PATH)/lrrs_iopsetup_WorkRoutines.f90    \
     $(LRRS_TEST_PATH)/lrrs_Phasfunc_Sup_Accessory.f90    \
     $(LRRS_TEST_PATH)/lrrs_Phasfunc_LinSup_Accessory.f90 \
     $(LRRS_SETUP_PATH)/lrrs_LPbin_M_iopsetup.f90         \
     $(LRRS_TEST_PATH)/lrrs_2p5_sup_accessories_1.f90     \
     $(LRRS_TEST_PATH)/lrrs_LPbin_M_tester.f90
endif


# Utilities

SOURCES_UTIL_LRRS =
SOURCES_UTIL_LRRS +=   \
   $(UTIL_PATH)/lrrs_diff.f90

# Define pattern rules for creating object files:
# ###############################################

.SUFFIXES:

# For f90 source files
$(OBJ_PATH)/%.o : $(LRRS_DEF_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(LRRS_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(LRRS_SUP_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(LRRS_L_DEF_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(LRRS_L_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(LRRS_L_SUP_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(LRRS_TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(LRRS_SETUP_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# For utility source files
$(OBJ_PATH)/%.o : $(UTIL_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# Define object files
# ###################

# Bin tests

F90SOURCES_BIN := $(notdir $(filter %.f90, $(SOURCES_BIN)))
F90OBJECTS_BIN := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_BIN)))

F90SOURCES_BIN_LP := $(notdir $(filter %.f90, $(SOURCES_BIN_LP)))
F90OBJECTS_BIN_LP := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_BIN_LP)))

#  Utilities

F90SOURCES_UTIL_LRRS := $(notdir $(filter %.f90, $(SOURCES_UTIL_LRRS)))
F90OBJECTS_UTIL_LRRS := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_UTIL_LRRS)))

# Define desired targets
# ######################

# Test sets

serial_bin:   lrrs_o3bin_M_tester.exe lrrs_LPbin_M_tester.exe

REENG_bin:   lrrs_o3bin_M_REENG_tester.exe lrrs_LPbin_M_REENG_tester.exe 

# Bin tests (serial)

lrrs_o3bin_M_tester.exe: $(F90OBJECTS_BIN)
	$(FC) $^ -o $@
lrrs_LPbin_M_tester.exe: $(F90OBJECTS_BIN_LP)
	$(FC) $^ -o $@

# Bin tests (re-engineered with internal OpenMP)

lrrs_o3bin_M_REENG_tester.exe: $(F90OBJECTS_BIN)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@
lrrs_LPbin_M_REENG_tester.exe: $(F90OBJECTS_BIN_LP)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@

#lrrs_o3bin_M_REENG_tester.exe: $(F90OBJECTS_BIN)
#	$(FC) $^ -o $@
#lrrs_LPbin_M_REENG_tester.exe: $(F90OBJECTS_BIN_LP)
#	$(FC) $^ -o $@

# Diff utility

lrrs_diff: $(F90OBJECTS_UTIL_LRRS)
	$(FC) $^ -o $@

# Make clean

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

