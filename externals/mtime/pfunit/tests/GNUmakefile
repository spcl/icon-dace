# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

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
