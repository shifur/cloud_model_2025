#
#  MAKEFILE 
#  Compiles advection code required to run GCE column model
#
#  Derek Posselt
#  02 October 2008
#  University of Michigan
# 
#  LINUX COMPILE OPTIONS

# Intel FORTRAN compiler
#F90 = ifort
#LD = ifort

# F90FLAGS = -i4 -r8 -O0 -g -traceback -check all -align dcommons -fp-model precise -c 
# LDFLAGS = -i4 -r8 -O0 -g -traceback -check all -align dcommons -fp-model precise

# F90FLAGS = -r8 -O2 -align dcommons -fp-model precise -c
# LDFLAGS =  -r8 -O2 -align dcommons -fp-model precise

# GNU Fortran Compiler
F90 = gfortran
LD = gfortran

F90FLAGS = -O2 -finit-local-zero -fPIC -c
LDFLAGS = -O2 -finit-local-zero -fPIC

CMD  = dynamics.a

#-----------------------------------------------------------------------

DYN_SRC   = fadv.for fadvuw.for fadvect.for advectn.for advect.for \
            advectak.for akcoef.for satdt_ak.for two_ak.for

DYN_OBJ   = fadv.o   fadvuw.o   fadvect.o   advectn.o   advect.o   \
            advectak.o   akcoef.o   satdt_ak.o   two_ak.o

all:	$(CMD)

$(CMD): $(DYN_OBJ)
	$(AR) -r $(CMD) $(DYN_OBJ)

# Code used in advection

fadv.o:  fadv.for
	$(F90) $(F90FLAGS) fadv.for

fadvuw.o:  fadvuw.for
	$(F90) $(F90FLAGS) fadvuw.for

fadvect.o:  fadvect.for
	$(F90) $(F90FLAGS) fadvect.for

advectn.o:  advectn.for
	$(F90) $(F90FLAGS) advectn.for

advect.o:  advect.for
	$(F90) $(F90FLAGS) advect.for

advectak.o:  advectak.for
	$(F90) $(F90FLAGS) advectak.for

akcoef.o:  akcoef.for
	$(F90) $(F90FLAGS) akcoef.for

satdt_ak.o:  satdt_ak.for
	$(F90) $(F90FLAGS) satdt_ak.for

two_ak.o:  two_ak.for
	$(F90) $(F90FLAGS) two_ak.for

clean:
	rm -f *.o *.a
