#
# Makefile for G-value calcuation
# requires Fortran 2003
#

include	../make.conf

.SUFFIXES:
.SUFFIXES:	.f90 .o

.f90.o:
	$(FC) $(INC) $(F90FLAGS) -c -o $@ $<

all:	src
src:	gvalue

gvalue:	gvalue.o mod_constants.o mod_file_utils.o
gvalue: class_grid.o class_medium.o class_work.o class_mixture.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

gvalue.o: mod_constants.o mod_file_utils.o
gvalue.o: class_grid.o class_medium.o
class_medium.o: mod_file_utils.o mod_constants.o class_grid.o class_mixture.o
class_work.o: mod_file_utils.o class_medium.o class_grid.o class_mixture.o

clean:
	$(RM) -f *~
	$(RM) -f *.o
	$(RM) -f *.mod
	$(RM) -f gvalue
