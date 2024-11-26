include ../../make.conf
include ../make_test.conf

.PHONY: run clean

TEST_NAME = test_neon

BIN	= gvalue
BINDIR	= ../../src/
PARAMDIR = ../../parameters
PARAM	= neon_nml.in
INPUT	= test_neon.in

DATAFILE = $(shell echo $(INPUT) | sed "s/.in/\.dat/")
DATAFILE += $(shell echo $(INPUT) | sed "s/.in/\_degradation*.dat/")
DATAFILE += $(shell echo $(INPUT) | sed "s/.in/\_gvalue.dat/")

all:	run clean

run:
	if [ -L $(BIN) ]; then \
	  $(RM) -f $(BIN) ;\
	fi
	$(LN) -s $(BINDIR)/$(BIN)
	if [ -L $(PARAM) ]; then \
	  $(RM) -f $(PARAM) ;\
	fi
	$(LN) -s $(PARAMDIR)/$(PARAM)

	./$(BIN) $(INPUT)
	pytest -s -v --junit-xml $(TEST_NAME).xml

clean:
	$(RM) -f *~
	$(RM) -f $(BIN)
	$(RM) -f $(PARAM)
	$(RM) -f $(DATAFILE)
