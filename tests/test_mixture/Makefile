include ../../make.conf
include ../make_test.conf

.PHONY: run clean

TEST_NAME = test_mixture

BIN	= gvalue
BINDIR	= ../../src/
PARAMDIR = ../../parameters
PARAM	= helium_nml.in neon_nml.in
INPUT	= test_mixture.in

DATAFILE = $(shell echo $(INPUT) | sed "s/.in/\.dat/")
DATAFILE += $(shell echo $(INPUT) | sed "s/.in/\_degradation*.dat/")
DATAFILE += $(shell echo $(INPUT) | sed "s/.in/\_gvalue.dat/")

all:	run clean

run:
	if [ -L $(BIN) ]; then \
	  $(RM) -f $(BIN) ;\
	fi
	$(LN) -s $(BINDIR)/$(BIN)
	for pfile in $(PARAM); do \
	  if [ -L $${pfile} ]; then \
	    $(RM) -f $${pfile} ;\
	  fi ;\
	  $(LN) -s $(PARAMDIR)/$${pfile} ;\
	done

	./$(BIN) $(INPUT)
	pytest -s -v --junit-xml $(TEST_NAME).xml

clean:
	$(RM) -f *~
	$(RM) -f $(BIN)
	$(RM) -f $(PARAM)
	$(RM) -f $(DATAFILE)
