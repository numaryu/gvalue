#
# Makefile for G-value calcuation
# requires Fortran 2003
#

include	../make.conf

.SUFFIXES:
.SUFFIXES:	.f90 .o

.f90.o:
	$(FC) $(INC) $(F90FLAGS) -c -o $@ $<

all:	gvalue

gvalue:	gvalue.o mod_constants.o mod_file_utils.o
gvalue: class_grid.o class_orbital.o class_work.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

gvalue.o: mod_constants.o mod_file_utils.o
gvalue.o: class_grid.o class_orbital.o
class_orbital.o: mod_file_utils.o mod_constants.o class_grid.o

clean:
	rm -f *~
	rm -f *.o
	rm -f *.mod
	rm -f gvalue
