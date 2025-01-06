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
# F90 = ifort
# LD = ifort

# F90FLAGS = -i4 -r8 -O0 -g -check all -traceback -align commons -fp-model precise -c 
# LDFLAGS = -i4 -r8 -O0 -g -check all -traceback -align commons -fp-model precise

# F90FLAGS = -r8 -O2 -align dcommons -fp-model precise -c
# LDFLAGS =  -r8 -O2 -align dcommons -fp-model precise

# GNU Fortran Compiler
F90 = gfortran
LD = gfortran

F90FLAGS = -O2 -finit-local-zero -fPIC -c
LDFLAGS = -O2 -finit-local-zero -fPIC

CMD  = radiation.a

#-----------------------------------------------------------------------

RAD_SRC   = b10kdis.for ch4exps.for co2kdis.for comexps.for \
	h2okdis.for opt4.for    rflx.for    terp1.for   \
	cfcexps.for ch4kdis.for coeff.for   comkdis.for \
	irrad.for   pradrat.for sorad.for   trip.for    \
	cfckdis.for cldovlp.for column.for  fito3.for   \
	n2okdis.for radrat.for  tablup.for  b10exps.for \
	solir.for   soluv.for   zang.for    cldflx.for  \
	co2exps.for conexps.for deledd.for  h2oexps.for \
	interp.for  search.for  n2oexps.for cldscale.for

RAD_OBJ   = b10kdis.o   ch4exps.o   co2kdis.o   comexps.o   \
	h2okdis.o   opt4.o      rflx.o      terp1.o     \
	cfcexps.o   ch4kdis.o   coeff.o     comkdis.o   \
	irrad.o     pradrat.o   sorad.o     trip.o      \
	cfckdis.o   cldovlp.o   column.o    fito3.o     \
	n2okdis.o   radrat.o    tablup.o    b10exps.o   \
	solir.o     soluv.o     zang.o      cldflx.o    \
	co2exps.o   conexps.o   deledd.o    h2oexps.o   \
	interp.o    search.o    n2oexps.o   cldscale.o

all:	$(CMD)

$(CMD): $(RAD_OBJ)
	$(AR) -r $(CMD) $(RAD_OBJ)

# Code used in radiation

b10kdis.o:  b10kdis.for
	$(F90) $(F90FLAGS) b10kdis.for

ch4exps.o:  ch4exps.for
	$(F90) $(F90FLAGS) ch4exps.for

co2kdis.o:  co2kdis.for
	$(F90) $(F90FLAGS) co2kdis.for

comexps.o:  comexps.for
	$(F90) $(F90FLAGS) comexps.for

h2okdis.o:  h2okdis.for
	$(F90) $(F90FLAGS) h2okdis.for

opt4.o:     opt4.for
	$(F90) $(F90FLAGS) opt4.for

rflx.o:     rflx.for
	$(F90) $(F90FLAGS) rflx.for

terp1.o:    terp1.for
	$(F90) $(F90FLAGS) terp1.for

cfcexps.o:  cfcexps.for
	$(F90) $(F90FLAGS) cfcexps.for

ch4kdis.o:  ch4kdis.for
	$(F90) $(F90FLAGS) ch4kdis.for

coeff.o:    coeff.for
	$(F90) $(F90FLAGS) coeff.for

comkdis.o:  comkdis.for
	$(F90) $(F90FLAGS) comkdis.for

irrad.o:    irrad.for
	$(F90) $(F90FLAGS) irrad.for

pradrat.o:  pradrat.for
	$(F90) $(F90FLAGS) pradrat.for

sorad.o:    sorad.for
	$(F90) $(F90FLAGS) sorad.for

trip.o:     trip.for
	$(F90) $(F90FLAGS) trip.for

cfckdis.o:  cfckdis.for
	$(F90) $(F90FLAGS) cfckdis.for

cldovlp.o:  cldovlp.for
	$(F90) $(F90FLAGS) cldovlp.for

column.o:   column.for
	$(F90) $(F90FLAGS) column.for

fito3.o:    fito3.for
	$(F90) $(F90FLAGS) fito3.for

n2okdis.o:  n2okdis.for
	$(F90) $(F90FLAGS) n2okdis.for

radrat.o:   radrat.for
	$(F90) $(F90FLAGS) radrat.for

tablup.o:   tablup.for
	$(F90) $(F90FLAGS) tablup.for

b10exps.o:  b10exps.for
	$(F90) $(F90FLAGS) b10exps.for

solir.o:    solir.for
	$(F90) $(F90FLAGS) solir.for

soluv.o:    soluv.for
	$(F90) $(F90FLAGS) soluv.for

zang.o:     zang.for
	$(F90) $(F90FLAGS) zang.for

cldflx.o:   cldflx.for
	$(F90) $(F90FLAGS) cldflx.for

co2exps.o:  co2exps.for
	$(F90) $(F90FLAGS) co2exps.for

conexps.o:  conexps.for
	$(F90) $(F90FLAGS) conexps.for

deledd.o:   deledd.for
	$(F90) $(F90FLAGS) deledd.for

h2oexps.o:  h2oexps.for
	$(F90) $(F90FLAGS) h2oexps.for

n2oexps.o:  n2oexps.for
	$(F90) $(F90FLAGS) n2oexps.for

cldscale.o:  cldscale.for
	$(F90) $(F90FLAGS) cldscale.for

interp.o:  interp.for
	$(F90) $(F90FLAGS) interp.for

search.o:  search.for
	$(F90) $(F90FLAGS) search.for

clean:
	rm -f *.o *.a
