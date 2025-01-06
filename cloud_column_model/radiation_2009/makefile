#
#  MAKEFILE 
#  Compiles radiation code required to run GCE column model
#
#  Derek Posselt
#  23 July 2008
#  University of Michigan
# 
#  LINUX COMPILE OPTIONS

# Intel FORTRAN compiler
#F90 = ifort
#LD = ifort

# F90FLAGS = -i4 -r8 -O0 -g -check all -traceback -align commons -fp-model precise -c 
# LDFLAGS = -i4 -r8 -O0 -g -check all -traceback -align commons -fp-model precise

#F90FLAGS = -i4 -r8 -O2 -align dcommons -fp-model precise -c 
#LDFLAGS =  -i4 -r8 -O2 -align dcommons -fp-model precise

# GNU fortran compiler
F90 = gfortran
LD = gfortran

F90FLAGS = -O2 -finit-local-zero -fPIC -c
LDFLAGS = -O2 -finit-local-zero -fPIC

CMD  = radiation_2009.a

#-----------------------------------------------------------------------

RAD_SRC   = module_ra_goddard_gce.f90 update_radiation.f90 zang.f90

RAD_OBJ   = module_ra_goddard_gce.o   update_radiation.o   zang.o

all:	$(CMD)

$(CMD): $(RAD_OBJ)
	$(AR) -r $(CMD) $(RAD_OBJ)

# Code used in radiation

module_ra_goddard_gce.o:  module_ra_goddard_gce.f90
	$(F90) $(F90FLAGS) module_ra_goddard_gce.f90

update_radiation.o:  update_radiation.f90
	$(F90) $(F90FLAGS) update_radiation.f90

zang.o:  zang.f90
	$(F90) $(F90FLAGS) zang.f90

clean:
	rm -f *.o *.mod *.a
