#
#  MAKEFILE 
#  Compiles code required to run set of GCE microphysics/radiation diagnostics
#
#  Derek Posselt
#  23 July 2008
#  University of Michigan
# 
#  LINUX COMPILE OPTIONS

# Intel FORTRAN compiler
#F90 = ifort
#LD = ifort

# F90FLAGS = -r8 -O0 -g -traceback -check all -align dcommons -fp-model precise -c 
# LDFLAGS  = -r8 -O0 -g -traceback -check all -align dcommons -fp-model precise
#F90FLAGS = -r8 -O2 -align dcommons -fp-model source -c 
#LDFLAGS  = -r8 -O2 -align dcommons -fp-model source

#F90FLAGS1 = -r8 -O2 -fp-model source -c 
#LDFLAGS1 =  -r8 -O2 -fp-model source
# F90FLAGS1 = -r8 -O0 -fp-model precise -g -traceback -check all -c 
# LDFLAGS1 =  -r8 -O0 -fp-model precise -g -traceback -check all

# GNU Fortran Compiler
F90 = gfortran
LD = gfortran

F90FLAGS = -O2 -finit-local-zero -ffree-form -ffree-line-length-0 -fPIC -c
LDFLAGS = -O2 -finit-local-zero -ffree-form -ffree-line-length-0 -fPIC

#F90FLAGS = -O0 -g -fdefault-real-8 -finit-local-zero -Wall -fbacktrace -Wunderflow -Warray-bounds -ffree-form -ffree-line-length-0 -c
#LDFLAGS = -O0 -g -fdefault-real-8 -finit-local-zero -Wall -fbacktrace -Wunderflow -Warray-bounds -ffree-form -ffree-line-length-0 

F90FLAGS1 = -O2 -c
LDFLAGS1 = -O2

# DRIVER PROGRAM
PROGMAIN = cloud_column_model.f90

# EXECUTABLE NAMES
CMD  = cloud_column_model.x
RAD_LIB1 = radiation.a
RAD_LIB2 = radiation_2009.a
MPHYS_LIB = microphysics.a
DYN_LIB = dynamics.a

#-----------------------------------------------------------------------


SOURCES = drive_micro_rad.f90 ideal_setup.f90 check_file.f90 simobs.f90 \
	  read_forcing.f90 ensemble_fcst.f90 

OBJECTS = drive_micro_rad.o   ideal_setup.o   check_file.o   simobs.o   \
	  read_forcing.o   ensemble_fcst.o  

SOURCES1 = check_file.f90

OBJECTS1 = check_file.o

#all:	$(CMD) rmobj
all:	$(CMD)

$(CMD): $(OBJECTS) $(MPHYS_OBJ) $(RAD_OBJ) $(MPHYS_LIB) $(RAD_LIB1) $(RAD_LIB2) $(DYN_LIB)
	$(LD) $(LDFLAGS) $(OBJECTS) $(MPHYS_LIB) $(RAD_LIB1) $(RAD_LIB2) $(DYN_LIB) -o $(@) $(PROGMAIN)

# Code used in main

drive_micro_rad.o: drive_micro_rad.f90 $(MPHYS_LIB) $(RAD_LIB) $(DYN_LIB)
	$(F90) $(F90FLAGS) drive_micro_rad.f90

ideal_setup.o: ideal_setup.f90
	$(F90) $(F90FLAGS) ideal_setup.f90

read_forcing.o: read_forcing.f90
	$(F90) $(F90FLAGS) read_forcing.f90

simobs.o: simobs.f90
	$(F90) $(F90FLAGS) simobs.f90

ensemble_fcst.o: ensemble_fcst.f90 drive_micro_rad.f90
	$(F90) $(F90FLAGS) ensemble_fcst.f90

check_file.o: check_file.f90
	$(F90) $(F90FLAGS1) check_file.f90

clean: rmobj
	rm -f $(CMD) $(CMD1)

rmobj:
	rm -f core *.mod *.o
