.PHONY: clean tests all

# Individual stuff (Christian)
F90=gfortran
PFUNIT=/informatik/home/zoller/scratch/lib/pfunit-serial
LIB_DIR=/informatik/home/zoller/scratch/lib

TOP_DIR := $(shell pwd)
TEST_DIR=$(TOP_DIR)/tests
SRC_DIR=$(TOP_DIR)/../src
INCL_DIR=$(TOP_DIR)/../include

VPATH = . $(SRC_DIR) $(TEST_DIR)

# Decide the file extensions.
include $(PFUNIT)/include/extensions.mk

EXE = tests

PFUNIT_LIBS = -L$(PFUNIT)/lib -lpfunit 
PFUNIT_INCL = -I$(PFUNIT)/mod -I$(PFUNIT)/include
PFUNIT_DRV = $(PFUNIT)/include/driver.F90

SRC_OBJS = $(SRC_DIR)/*$(OBJ_EXT)
TEST_INCL =  -I$(TEST_DIR)
TEST_OBJS = $(TEST_DIR)/*$(OBJ_EXT)

FFLAGS += -I$(INCL_DIR) $(PFUNIT_INCL) $(TEST_INCL) $(PFUNIT_LIBS)

all: clean $(EXE) # Kein Plan, weshalb das clean nötig ist
	./$(EXE)

$(EXE): testSuites.inc from_example_test.pf tests
	$(F90) -o $@ $(PFUNIT_DRV) $(TEST_OBJS) $(SRC_OBJS) $(FFLAGS)

clean:
	make -C $(TEST_DIR) clean
	rm -f $(EXE) *$(OBJ_EXT)

tests: 
	make -C $(TEST_DIR) tests

export SRC_DIR
export F90
export FFLAGS
export PFUNIT

#PFUNITROOT     =/informatik/home/zoller/scratch/lib/pfunit-serial
#PFUNIT_INCLUDE = -I$(PFUNITROOT)/mod -I$(PFUNITROOT)/include
#PFUNIT_LIB     = -L$(PFUNIT)/lib -lpfunit 
#
#VPATH         = ../src tests
#
#SRC_DIR  = ../src
#
#INCLUDE = -I../include
#LIBS    = $(PFUNIT_LIB)
#
#F90     = gfortran
#FFLAGS  = $(INCLUDE) -Wall
#
#%.f90: %.pf
#	$(CHECKROOT)/bin/checkmk $< > $@
#
#all: example_test
#	 ./example_test
#
#example_test: example_test.o
#	$(CC) $(CFLAGS) -o $@ $< $(SRC_DIR)/mtime_*.o $(LIBS)  
#
#example_test.o: example_test.c
#
#example_test.c: example_test.check
#
#
#clean:
#	rm -f example_test *.o *.c
