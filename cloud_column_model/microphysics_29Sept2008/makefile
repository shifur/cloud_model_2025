#
#  MAKEFILE 
#  Compiles microphysics code required to run GCE column model
#
#  Derek Posselt
#  23 July 2008
#  University of Michigan
# 
#  LINUX COMPILE OPTIONS

# Intel FORTRAN compiler
# F90 = ifort
# LD = ifort

# F90FLAGS = -r8 -O2 -align dcommons -fp-model precise -c
# LDFLAGS =  -r8 -O2 -align dcommons -fp-model precise

# GNU Fortran Compiler
F90 = gfortran
LD = gfortran

F90FLAGS = -O2 -finit-local-zero -fPIC -c
LDFLAGS = -O2 -finit-local-zero -fPIC

CMD  = microphysics.a

#-----------------------------------------------------------------------

MPHYS_SRC   = consatrh.for saticerh.for tervrh.for gamma1.for

MPHYS_OBJ   = consatrh.o   saticerh.o   tervrh.o   gamma1.o

all:	$(CMD)

$(CMD): $(MPHYS_OBJ)
	$(AR) -r $(CMD) $(MPHYS_OBJ)

# Code used in microphysics

consatrh.o:  consatrh.for
	$(F90) $(F90FLAGS) consatrh.for

tervrh.o:  tervrh.for
	$(F90) $(F90FLAGS) tervrh.for

gamma1.o:  gamma1.for
	$(F90) $(F90FLAGS) gamma1.for

saticerh.o:  saticerh.for
	$(F90) $(F90FLAGS) saticerh.for

clean:
	rm -f *.o *.a
