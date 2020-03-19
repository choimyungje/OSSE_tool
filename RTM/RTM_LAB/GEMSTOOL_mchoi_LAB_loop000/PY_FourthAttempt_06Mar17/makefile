# Define some variables

MOD_PATH  = mod
OBJ_PATH  = obj

MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o

# Define default shell make will use

#SHELL = /bin/bash

# Define FORTRAN90 compiler to use 
#  gfortran or Ifort

FC = gfortran
#FC = ifort

# Define FORTRAN90 compiler flags

FFLAGS = -c

# Additional flags for gfortran
ifeq ($(FC), gfortran)
	FFLAGS_TRACE = -Wall -C -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
	FFLAGS := $(FFLAGS) $(FFLAGS_TRACE)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -J$(MOD_PATH)
endif

# Additional flags for ifort
ifeq ($(FC), ifort)
	FFLAGS_TRACE = -warn all -check all -traceback
	FFLAGS := $(FFLAGS) $(FFLAGS_TRACE)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -module $(MOD_PATH)
endif

#  Additional Optimization and Debug flags (Uncomment as required)
#	FFLAGS_OPTIM :=  -O3
#	FFLAGS_DEBUG := -g $(FFLAGS_TRACE)
#	FFLAGS := $(FFLAGS) $(FFLAGS_DEBUG) 
#	FFLAGS := $(FFLAGS) $(FFLAGS_OPTIM) 

# Define list of source files
# (Note: ordering is important because of dependencies)

SOURCES_2 =
SOURCES_2 +=   \
   SMCConv_Interp.f90    \
   SMCConv_Develop.f90   \
   SMCConv_Expand.f90    \
   PY_distribution.f90   \
   Extract_TriAxData_Mk2.f90   \
   Create_TriAx_for_GEMSTOOL_Mk4.f90  \
   WriteGT.f90             \
   Create_TriAx_Driver_Mk4.f90

# Define pattern rules for creating object files:

.SUFFIXES:
$(OBJ_PATH)/%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

# Define object files

F90SOURCES_2 := $(notdir $(filter %.f90, $(SOURCES_2)))
F90OBJECTS_2 := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_2)))

ALL: Create_TriAx_Driver_Mk4.exe

Create_TriAx_Driver_Mk4.exe: $(F90OBJECTS_2)
	$(FC) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

