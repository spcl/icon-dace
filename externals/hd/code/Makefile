# Makefile - Makefile for compiling the HD model
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann and Ha Ho-Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
####################################################
#
#    Makefile for the Fortran 90 program
#
####################################################
#
#
####################################################
#
#    Declaration of Variables
#
####################################################
#
#.SILENT:
#
STDROOT      = ${PWD}
SRCDIR       = $(STDROOT)/src
OBJDIR       = $(STDROOT)/obj

#
#########
#  include file with compiler options and libraries according to your needs
#########
#
include Fopts
#

####################################################
#
#     Declaration of the Object Files
#
####################################################
#
include ObjFiles
#
#
####################################################
#
#     Actions that can be performed
#
####################################################
#
exe:    parexe
#
#SRC=$(SRCDIR)/*
#TAGS: $(SRC)
#	etags $(SRC)
#tags: $(SRC)
#	ctags $(SRC)
#
parexe: $(HDOBJ)
	echo linking parallel program
	( $(LDPAR) -o $(PROGRAM) $(HDOBJ) $(LDFLG) )
#
clean:
	echo cleaning up
	rm -f $(PROGRAM)
	rm -f $(OBJDIR)/*.o
	rm -f $(OBJDIR)/*.T
	rm -f $(OBJDIR)/*.mod
	rm -f $(OBJDIR)/*.lst
#+++++
	rm -f $(OBJDIR)/*.L
#
####################################################
#
#     Dependencies of the Data Modules
#
####################################################
#
include ObjDependencies
#
