#
# Makefile for G-value calcuation
# requires Fortran 2003
#

include make.conf

.PHONY:	all src clean

SUBDIRS	= src

all clean:
	for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir -f Makefile $@ ;\
	done
