#
# Makefile for tests
#

include ../make.conf

.PHONY: tests

SUBDIRS = test_helium

tests:
	for dir in $(SUBDIRS); do \
	  echo "Testing in $$dir" ;\
	  $(MAKE) -C $$dir -f Makefile ;\
	done

