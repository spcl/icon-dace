.PHONY: clean tests

SRCS = $(wildcard *.pf)
OBJS = $(SRCS:.pf=$(OBJ_EXT))

tests: $(OBJS)

testSuites.inc: $(SRCS)

FFLAGS += -I$(SRC_DIR) -I$(PFUNIT)/mod

%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $<  $@

%$(OBJ_EXT): %.F90
	$(F90) -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod
