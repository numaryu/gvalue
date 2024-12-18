#
# Makefile for G-value calcuation
# requires Fortran 2003
#

include make.conf

.PHONY:	all src tests clean

SUBDIRS	= src

all:	src

src clean:
	for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir -f Makefile $@ ;\
	done

tests: src
	$(MAKE) -C tests -f Makefile $@
