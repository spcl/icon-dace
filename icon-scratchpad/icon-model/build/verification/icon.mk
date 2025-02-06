# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------
################################# DIRECTORIES ##################################

# Path to the root source directory:
srcdir:= /home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model

# Relative paths to the directories with the ICON source code:
subdirs:= src support

# Relative (to $(srcdir)) paths to the root source directories of all bundled
# packages:
bundled_srcdirs =  externals/ecrad externals/cdi externals/mtime externals/fortran-support externals/math-support externals/math-interpolation

# Relative paths to the root build directories of all bundled packages:
bundled_builddirs= externals/ecrad externals/cdi externals/mtime externals/fortran-support/build externals/math-support/build externals/math-interpolation/build

# Relative paths to the root build directories of the CMake-based bundled
# packages:
bundled_cmake_builddirs= externals/fortran-support/build externals/math-support/build externals/math-interpolation/build

# Relative paths to the root build directories of the Autoconf-based bundled
# packages:
bundled_config_builddirs= $(filter-out $(bundled_cmake_builddirs),$(bundled_builddirs))

# Relative path to the directory with the Fortran module files:
moddir:= mod

# When the two-step compilation is enabled, Fortran module files that are
# generated together with the object file become an unwanted side effect. More
# often than not, Fortran compilers do not have special flags to suppress the
# generation of the module files. Therefore, we use the following directory as
# a trash bin:
#moddir_null:= $(moddir)/null

# Relative path to the directory with the generated executables:
bindir:= bin

# Relative path to the directory with the generated archives:
libdir:= lib

# Paths to the installation directories:
prefix:= /usr/local
exec_prefix:= ${prefix}

########################### BUILDING TOOLS AND FLAGS ###########################

# Compilers and utilities:
SHELL= /bin/bash
AR= ar
ARFLAGS= cr
RANLIB= ranlib
CC= mpicc
FC= mpif90
CXX= g++
CUDACXX= 
HIPCXX= 
PYTHON= python3
CMAKE= cmake
PVCS= ${PYTHON} $(srcdir)/utils/pvcs.py
DEPLIST= ${PYTHON} $(srcdir)/utils/mkhelper/deplist.py
DEPGEN= ${PYTHON} $(srcdir)/utils/mkhelper/depgen.py
MODCMP= ${PYTHON} $(srcdir)/utils/mkhelper/fortmodcmp.py
DSL4JSB= ${PYTHON} $(srcdir)/externals/jsbach/scripts/dsl4jsb/dsl4jsb.py
SB2PP= 
DACE_SUBST_PP= ${PYTHON} $(srcdir)/sdfgs/utils/dace-substf.py
DACE_GENFI= ${PYTHON} $(srcdir)/sdfgs/utils/dace-genfi.py
DACE_PRINT_INTEGRATIONS= ${PYTHON} $(srcdir)/sdfgs/utils/dace-print-integrations.py
FPP= mpif90 -E -cpp
INSTALL= /usr/bin/install -c
MKDIR_P= /usr/bin/mkdir -p

# Fortran compiler flags:
FCFLAGS= -Iexternals/ecrad/mod -Iexternals/cdi/src -Iexternals/mtime/src -Iexternals/fortran-support/build/src/mod -Iexternals/math-support/build/src/mod -Iexternals/math-interpolation/build/src/mod -g -O2 -I/usr/include/ -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-missing-include-dirs -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE -DHAVE_FC_ATTRIBUTE_CONTIGUOUS -D__ICON__ -D__NO_ICON_LES__ -D__NO_ICON_OCEAN__ -D__NO_JSBACH__ -D__NO_JSBACH_HD__ -D__NO_QUINCY__ -D__NO_ICON_WAVES__ -D__NO_AES__ -D__ECRAD -D__NO_RTE_RRTMGP__ -D__NO_ICON_COMIN__ -DHAVE_ACM_LICENSE -DNOMPI -D__NO_ICON_TESTBED__
make_FCFLAGS= -I$(moddir)
ICON_FCFLAGS= -std=legacy

DEPGEN_FCFLAGS= 

# C compiler and preprocessor flags:
CFLAGS= -g -march=native -O3
CPPFLAGS= -I. -DHAVE_CONFIG_H -I/usr/include -I/usr/include/libxml2

# Secondary GPU compiler flags:
CUDAFLAGS= 
HIPFLAGS= 

# Linker flags and libraries:
LDFLAGS= -L/usr/lib/x86_64-linux-gnu/ -L/home/primrose/Work/IconGrounds/icon-dace/.dacecache/radiation/build -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/ -Wl,-rpath -Wl,/home/primrose/Work/IconGrounds/icon-dace/.dacecache/radiation/build
BUNDLED_LIBFILES= externals/math-interpolation/build/src/libmath-interpolation.a externals/math-support/build/src/libmath-support.a externals/fortran-support/build/src/libfortran-support.a externals/mtime/src/.libs/libmtime.a externals/cdi/src/.libs/libcdi_f2003.a externals/cdi/src/.libs/libcdi.a externals/ecrad/libradiation.a externals/ecrad/libifsrrtm.a externals/ecrad/libutilities.a externals/ecrad/libifsaux.a 
BUNDLED_EXTRA_LIBS= $(shell $(SHELL) ./collect.extra-libs)
# We use `-rpath` to make sure the correct dace generated libraries are found
#BUNDLED_EXTRA_LIBS+= -Wl,-rpath -Wl,$(shell pwd)/sdfgs/build -Lsdfgs/build
#BUNDLED_EXTRA_LIBS+= $(foreach sdfg,$(shell $(DACE_PRINT_INTEGRATIONS) $(srcdir)/sdfgs/integrations.yaml "sdfg_names"), -l$(sdfg))
LIBS= -leccodes -lnetcdff -lnetcdf -lopenblas -lradiation

############################ SOURCE FILE COLLECTION ############################

# Common infix for all preprocessed files:
pp_infix:= .pp-

# Auxiliary function. Receives a space-separated list of subdirectories inside
# $(srcdir) and a space-separated list of shell-like wildcards. Returns a
# space-separated list of relative (with respect to $(srcdir)) paths to files
# residing in the given subdirectories and matching at least one of the given
# patterns. Files, as well as the second and higher order subdirectories, with
# names starting with a dot are ignored. Files containing the preprocessor
# infix '*$(pp_infix)*' are ignored. Duplicates are eliminated.
find_files= $(sort $(patsubst $(srcdir)/%,%, \
                $(shell find $(addprefix $(srcdir)/,$(1)) \
                    -name '.' -o -name '..' -o -name '.*' -prune \
                    -o \! \( -type d -o -name '*$(pp_infix)*' \) \
                    -a \( $(foreach p,$(2),-name $(p) -o) -false \) -print)))

# All source files:
all_files:= $(call find_files,$(subdirs),'*.f90' '*.inc' '*.incf' '*.c' '*.cu' '*.hip.cc')

# Relative (to $(srcdir)) paths to the root source directories of all extra
# packages to be passed to the version provenance collection tool:
PVCS_subdirs= $(bundled_srcdirs)

# JSBACH sources (must be preprocessed with DSL4JSB):
#JSBACH_subdir:= externals/jsbach/src
#JSBACH_dsl_files:= $(call find_files,$(JSBACH_subdir),'*.f90')
#PVCS_subdirs+= externals/jsbach

# DACE sources:
#DACE_subdir:= externals/dace_icon/src_for_icon
#all_files+= $(call find_files,$(DACE_subdir),'*.f90' '*.inc' '*.incf' '*.c')
#PVCS_subdirs+= externals/dace_icon

# EMVORADO sources:
#EMVORADO_subdir:= externals/emvorado
#all_files+= $(call find_files,$(EMVORADO_subdir),'*.f90' '*.incf')
#PVCS_subdirs+= $(EMVORADO_subdir)

# ART sources:
#ART_subdir:= externals/art
#all_files+= $(call find_files,$(ART_subdir),'*.f90' '*.inc' '*.incf')
#PVCS_subdirs+= $(ART_subdir)

# Fortran '.f90' source files:
f90_files:= $(filter %.f90,$(all_files))
# C '.c' source files:
c_files:= $(filter %.c,$(all_files))
# CUDA '.cu' source files:
#cu_files:= $(filter %.cu,$(all_files))
# HIP '.hip.cc' source files:
#hip_files:= $(filter %.hip.cc,$(all_files))

# Relative paths to the directories with Fortran include files:
inc_subdirs:= $(patsubst %/,%,$(sort $(dir $(filter %.inc %.incf,$(all_files)))))

# Extend Fortran compiler flags enabling the location of include files:
make_FCFLAGS+= $(foreach d,$(inc_subdirs),-I$(srcdir)/$(d) -I$(srcdir)/$(d))

############################# PREPROCESSING CHAIN ##############################

# List of all preprocessed files:
pp_files:=

######## DSL4JSB PREPROCESSING #########
#DSL4JSB_in_f90_files:= $(JSBACH_dsl_files)
########################################
#DSL4JSB_infix:= $(pp_infix)jsb
#DSL4JSB_out_f90_files:= $(DSL4JSB_in_f90_files:.f90=$(DSL4JSB_infix).f90)
#f90_files:= $(filter-out $(DSL4JSB_in_f90_files),$(f90_files)) $(DSL4JSB_out_f90_files)
#pp_files+= $(DSL4JSB_out_f90_files)

####### DACE SUBST PREPROCESSING #######
#DACE_SUBST_in_f90_files:= $(shell $(DACE_PRINT_INTEGRATIONS) $(srcdir)/sdfgs/integrations.yaml "icon_f90_files")
########################################
#DACE_SUBST_infix:= $(pp_infix)dsub

#DACE_SUBST_out_f90_files:= $(DACE_SUBST_in_f90_files:.f90=$(DACE_SUBST_infix).f90)
# Add interface files
#DACE_SUBST_out_f90_files+= $(shell $(DACE_PRINT_INTEGRATIONS) $(srcdir)/sdfgs/integrations.yaml "f90_interfaces")

#f90_files:= $(filter-out $(DACE_SUBST_in_f90_files),$(f90_files)) $(DACE_SUBST_out_f90_files)
#pp_files+= $(DACE_SUBST_out_f90_files)

#### EXPLICIT FORTRAN PREPROCESSING ####
FPP_in_f90_files:= $(f90_files)
########################################
FPP_infix:= $(pp_infix)fpp
FPP_out_f90_files:= $(FPP_in_f90_files:.f90=$(FPP_infix).f90)
f90_files:= $(filter-out $(FPP_in_f90_files),$(f90_files)) $(FPP_out_f90_files)
pp_files+= $(FPP_out_f90_files)

####### SERIALBOX2 PREPROCESSING #######
#SB2_in_f90_files:= $(f90_files)
########################################
#SB2_infix:= $(pp_infix)sb2
#SB2_out_f90_files:= $(SB2_in_f90_files:.f90=$(SB2_infix).f90)
#f90_files:= $(filter-out $(SB2_in_f90_files),$(f90_files)) $(SB2_out_f90_files)
#pp_files+= $(SB2_out_f90_files)

# Dependency files
# (this finalizes the list of files eligible for compilation):
dep_files:= $(addsuffix .d,$(f90_files) $(c_files) $(cu_files) $(hip_files)) c_binding.d

############################ GENERATED EXECUTABLES #############################

# Generated executable files:
exec_files:=

# Auxiliary function. Receives a space-separated list of relative paths to the
# original Fortran '.f90' files and returns a space-separated list of paths to
# the respective object files (not necessarily in the respective order though),
# which change depending on the preprocessing steps applied to the original
# source files.
get_obj_names= $(patsubst %.f90,%.o,$(filter $(foreach f,$(1),$(f:f90=)%f90),$(f90_files)))

# ICON executable:
ICON_exec_file:= $(bindir)/icon
ICON_prog_obj_files:= $(call get_obj_names,src/drivers/icon.f90) version.o
#ICON_prog_obj_files+= support/dummy_hip_events.o
exec_files+= $(ICON_exec_file)

############################## GENERATED ARCHIVES ##############################

# Generated archive (library) files:
lib_files:=

# Auxiliary function. Receives a space-separated list of paths to objects and
# returns a space-separated list of paths to their prerequisite objects. The
# input paths are excluded from the output.
get_prereq_objs= $(filter-out $(1),$(filter %.o, \
                     $(shell $(DEPLIST) @deplist.config -t $(1) -f $(dep_files))))

# Common convenience library (holds prerequisite objects of the generated
# executables):
common_lib_file:= $(libdir)/libicon.a
common_lib_prog_obj_files:= $(ICON_prog_obj_files)
common_lib_obj_files:= $(call get_prereq_objs,$(common_lib_prog_obj_files))
lib_files+= $(common_lib_file)

# Stamp files of the building subdirectories
# (this finalizes the list of generated files):
dir_files:= $(filter-out ./.dirstamp,$(addsuffix .dirstamp,$(sort $(dir $(dep_files) $(pp_files) $(exec_files) $(lib_files))) $(moddir)/))

########################## PGI/NVIDIA INLINE LIBRARY ###########################

# Name of the inline library (a directory):
#pgi_inlib_name:= $(libdir)/icon.il

# List of procedures to be inlined:
#pgi_inlib_procedures:=

# List of original Fortran '.f90' files containing $(pgi_inlib_procedures):
#pgi_inlib_f90_files:=

# List of files to be compiled with the GPU register limit flag:
##pgi_inlib_reglim_f90_files:=

# Lengthy space-separated lists $(pgi_inlib_procedures), $(pgi_inlib_f90_files)
# and $(pgi_inlib_reglim_f90_files) are set in a separate file for the sake of
# readability:
#include inlib.mk

# Auxiliary function. Turns a space-separated list into a comma-separated one.
join_comma_symbol:= ,
join_comma= $(subst $() ,$(join_comma_symbol),$(strip $(1)))

# Additional PGI/NVIDIA Fortran compiler flags enabling the inline library
# generation and usage:
#pgi_inlib_common_arg:= lib:$(pgi_inlib_name),name:$(call join_comma,$(pgi_inlib_procedures))
#pgi_inlib_ex_FCFLAGS= -Mextract=$(pgi_inlib_common_arg)
##dir_files+= $(moddir_null)/icon.il/.dirstamp
##pgi_inlib_ex_FCFLAGS+= -J$(moddir_null)/icon.il
#pgi_inlib_ex_FCFLAGS+= -J$(moddir)
#pgi_inlib_in_FCFLAGS= -Minline=$(pgi_inlib_common_arg),reshape

# Additional flag limiting the GPU register usage (valid only if the inline
# library is enabled):
##pgi_inlib_reglim_FCFLAGS= -gpu=maxregcount:96

########################################

# List of object files that correspond to the source files from
# $(pgi_inlib_f90_files) and are needed for the current configuration:
#pgi_inlib_obj_files:= \
#  $(filter $(common_lib_obj_files) $(common_lib_prog_obj_files), \
#    $(call get_obj_names,$(pgi_inlib_f90_files)))

# List of the respective dependency files:
#pgi_inlib_dep_files:= $(pgi_inlib_obj_files:.o=.f90.d)

# List of Fortran module files declared (or extended via submodules) that
# correspond to $(pgi_inlib_obj_files):
#pgi_inlib_mod_files:= \
#  $(filter %.mod.proxy, \
#    $(shell $(DEPLIST) -p $(pgi_inlib_obj_files) -f $(pgi_inlib_dep_files)))

# List of objects that are generated using the inline library:
#pgi_inlib_target_obj_files:= \
#  $(filter %.o, \
#    $(shell $(DEPLIST) -p $(pgi_inlib_obj_files) -f $(dep_files)))

# If a source file from $(pgi_inlib_f90_files) declares a submodule, the object
# file of the respective parent module should be generated using the inline
# library as well:
#pgi_inlib_target_obj_files+= \
#  $(patsubst %$(pp_infix)modstamp.f90,%.o, \
#    $(filter %.o %$(pp_infix)modstamp.f90, \
#      $(shell $(DEPLIST) --no-hints --max-depth 1 -t $(pgi_inlib_mod_files) -f $(dep_files))))

# To avoid circular dependencies, we need to account for situations when we
# have a dependency A -> B -> C, where A and C belong to $(pgi_inlib_obj_files)
# but B does not. In order to get all such B files and include them into
# $(pgi_inlib_obj_files), we find the intersection of two sets: dependencies
# (prerequisites) and dependents (targets) of $(pgi_inlib_obj_files):
#pgi_inlib_obj_files:= \
#  $(filter $(pgi_inlib_target_obj_files), \
#    $(shell $(DEPLIST) --no-hints -t $(pgi_inlib_obj_files) -f $(dep_files)))
#pgi_inlib_dep_files:= $(pgi_inlib_obj_files:.o=.f90.d)

# List of Fortran module and submodule files that must be created before the
# inline library:
#pgi_inlib_prereq_mod_files:= \
#  $(filter %.mod.proxy %.smod.sproxy, \
#    $(filter-out \
#      $(shell $(DEPLIST) --no-hints -p $(pgi_inlib_obj_files) -f $(pgi_inlib_dep_files)), \
#      $(shell $(DEPLIST) --no-hints -t $(pgi_inlib_obj_files) $(pgi_inlib_obj_files:.o=$(pp_infix)modstamp.f90) -f $(pgi_inlib_dep_files))))

############################### MAKEFILE TWEAKS ################################

# Silent rule prefixes:
V= 0
ifeq ($(V),0)
silent_AR=      @echo "  AR      " $@;
silent_CC=      @echo "  CC      " $@;
silent_CMAKE=   @echo "  CMAKE   " $(@D);
silent_CONFIG=  @echo "  CONFIG  " $(@D);
silent_CUDACXX= @echo "  CUDACXX " $@;
silent_DEPGEN=  @echo "  DEPGEN  " $@;
silent_DSL4JSB= @echo "  DSL4JSB " $@;
silent_FC=      @echo "  FC      " $@;
silent_FCEX=    @echo "  FC (EX) " $(@D);
silent_FCIL=    @echo "  FC (IL) " $@;
silent_FCLD=    @echo "  FCLD    " $@;
silent_FPP=     @echo "  FPP     " $@;
silent_GEN=     @echo "  GEN     " $@;
silent_HIPCXX=  @echo "  HIPCXX  " $@;
silent_MKDIR=   @echo "  MKDIR   " $(@D);
silent_MOD=     @echo "  MOD    <" $<;
silent_SB2=     @echo "  SB2PP   " $@;
silent_DACE_SUBST=@echo " DACE_SUBST_PP " $@;
silent_none=    @
endif

# Disable built-in suffix rules:
.SUFFIXES:
# Delete partially updated files:
.DELETE_ON_ERROR:
# Targets not associated with files:
.PHONY: all preprocess depend check install mostlyclean clean distclean \
        config-bundled check-bundled dummy-depend $(bundled_builddirs) \
        force-create-version sanitize-mod-proxies

# Makefiles of the bundled packages:
bundled_makefiles= $(addsuffix /Makefile,$(bundled_builddirs))

# Targets that do not need the dependency files:
nodep_targets:= preprocess depend mostlyclean clean distclean config-bundled \
                check-bundled dummy-depend $(bundled_builddirs) \
                $(bundled_makefiles) version.c version.o $(dep_files) \
                $(pp_files)

# Selective search paths:
vpath %.f90 $(srcdir)
vpath %.c $(srcdir)
vpath %.cu $(srcdir)
vpath %.hip.cc $(srcdir)
vpath sdfgs/%.sdfgz $(srcdir)
vpath sdfgs/%_module_definitions.yaml $(srcdir)
vpath sdfgs/%_associations.yaml $(srcdir)

################### PATTERN- AND TARGET-SPECIFIC ASSIGNMENTS ###################

# Pattern- and target-specific assignments can lead to unexpected side-effects.
# To have better control over them, all such assignments should be done in this
# section, once all global variables have their final values.

# Pattern- and target-specific assignments are propagated to the prerequisites
# and override the global assignments. Therefore, we introduce the following
# match-anything pattern assignments to prevent that:
#%: make_FCFLAGS:= $(make_FCFLAGS)
%: ICON_FCFLAGS:= $(ICON_FCFLAGS)
#   The following variables have to be expanded recursively (i.e. they are
#   supposed to be lazily expanded when running the recipes). The eval/value
#   trick below prevents propagation of the target-specific values to the
#   prerequisites while keeping the recursive flavor of the variables:
#$(eval %: silent_FC= $(value silent_FC))

# Target-specific variables enabling the PGI/NVIDIA inlining during compilation:
#$(pgi_inlib_target_obj_files): silent_FC= $(silent_FCIL)
#$(pgi_inlib_target_obj_files): make_FCFLAGS+= $(pgi_inlib_in_FCFLAGS)

# Compile $(pgi_inlib_reglim_f90_files) with the GPU register limit flag:
##$(call get_obj_names,$(pgi_inlib_reglim_f90_files)): make_FCFLAGS+= $(pgi_inlib_reglim_FCFLAGS)

# Pattern-specific Fortran compiler flags:
#   Note the difference in how older versions of GNU Make (at least up to 3.81)
#   treat pattern-specific variable ASSIGNMENTS and pattern RULES: the last
#   matching ASSIGNMENT and the first matching RULE win. Newer versions of GNU
#   Make (starting at least 4.1) give preference to the most specialized pattern
#   in both cases. This means that, as long as we support older version of GNU
#   Make, we have to declare the most general ASSIGNMENTS first and the most
#   general RULES last.


############################# USER INTERFACE RULES #############################

# Default rule:
all: $(exec_files)

# Explicit preprocessing rule:
preprocess: $(pp_files)

# Explicit dependency generation rule:
depend: $(dep_files)

# Generate a list of files that are going to be compiled in the current
# configuration:
srclist: $(dep_files)
	$(silent_GEN)$(DEPLIST) -t $(ICON_prog_obj_files) -f $^ 2>/dev/null | sed '/\(\.o\|\.mod\.proxy\|.smod\.sproxy\)$$/d' | sort > $@

# Check rule:
check: all check-bundled

# A general rule for clean recipes is to delete only those files that have been
# created explicitly, abandoning the idea of trying to delete all possible side
# effect files, e.g. Cray's *.i and *.lst files, MacOS's *.dSYM directories,
# Python's cache files, etc.

# Delete the results of compilation and linking:
mostlyclean: $(bundled_builddirs)
	rm -f $(exec_files) $(lib_files) $(addsuffix .o,$(basename $(f90_files) $(c_files) $(cu_files) $(hip_files))) version.o
	rm -f $(moddir)/*.mod $(moddir)/*.mod.proxy
	rm -f $(moddir)/*.smod $(moddir)/*.smod.sproxy
#	rm -rf $(moddir_null) $(f90_files:.f90=$(pp_infix)modstamp.f90)
#	rm -rf $(pgi_inlib_name)

# Delete files generated at the building stage (strictly speaking, we should
# also translate 'clean' into 'distclean' for the bundled packages when the
# delayed configuration is enabled but that would be counterproductive,
# therefore, we do not do that):
clean_files= srclist version.c $(pp_files)

# With all components and all preprocessing steps enabled, the list above
# contains ~3k filenames, which might result into the following error if we
# provide it on the command line as is:
#   make[1]: /bin/bash: Argument list too long
# To circumvent the issue, we have to split $(clean_files) into two lists.

# Auxiliary function. Receives a name of a variable containing a space-separated
# list. Splits the list into two ones that are twice as short (in terms of the
# number of elements). Each element of the input list can be found in one of the
# two output lists. In other words, if the input list L is a set, the result is
# two disjoint sets L1 and L2 such that |L1| - |L2| <= 1 and L = L1 U L2.
# Variables that contain the results have names <NAME>_1 and <NAME>_2, where
# <NAME> is the input name.
split_list= $(foreach e,$($(1)), \
                $(if $(filter $(words $($(1)_1)),$(words $($(1)_2))), \
                    $(eval $(1)_1+=$(e)),$(eval $(1)_2+=$(e))))

$(call split_list,clean_files)
clean: mostlyclean
	$(if $(clean_files_1),rm -f $(clean_files_1))
	$(if $(clean_files_2),rm -f $(clean_files_2))

# Delete everything generated at the configure stage (and clean the created
# directories if they are empty):
distclean: clean
	rm -f config.h config.log config.status depgen.c.config depgen.f90.config deplist.config collect.extra-libs # inlib.mk
	rm -f run/collect.set-up.info run/set-up.info run/SETUP.config
	rm -f $(dep_files)
	rm -f $(dir_files)
	test -z "$(bundled_cmake_builddirs)" || rm -rf $(bundled_cmake_builddirs)
	@for dir in $(bindir) $(libdir) $(moddir); do \
          if test -d "$$dir"; then \
            echo "find '$$dir' -type d -empty -delete"; \
            find "$$dir" -type d -empty -delete; \
          fi; \
	done
	@if test '.' != '$(srcdir)'; then \
	  for dir in $(subdirs) externals run; do \
	    if test -d "$$dir"; then \
	      echo "find '$$dir' -type d -empty -delete"; \
	      find "$$dir" -type d -empty -delete; \
	    fi; \
	  done; \
	fi
	rm -f icon.mk

# Installation rule:
install: all
	$(INSTALL) -d $(DESTDIR)${exec_prefix}/bin && $(INSTALL) $(exec_files) $(DESTDIR)${exec_prefix}/bin
	$(silent_none)for d in $(filter externals/hd externals/comin/build,$(bundled_subdirs)); do \
	  $(MAKE) -C "$$d" $@ V=$(V); \
	done
#	$(silent_none)$(MAKE) -C externals/mtime/python $@ V=$(V)
#	$(silent_none)$(MAKE) -C externals/yac/python $@ V=$(V)

# Explicit configuration rule for the bundled packages:
config-bundled: $(bundled_makefiles)

# Check rule for the bundled packages (run the tests serially to make the output
# readable and prevent from possible overloading due to multiple simultaneous
# MPI runs):
check-bundled: config-bundled $(bundled_cmake_builddirs) $(if $(filter check,$(MAKECMDGOALS)),all) # run the tests after 'all' when 'make check'
	@test -z '$(bundled_builddirs)' && echo "The list of bundled packages is empty: nothing to check." && exit 0; \
	fail=; pass=; \
	for d in $(bundled_config_builddirs); do \
	  if $(MAKE) -C "$$d" check V=$(V); then pass="$$pass$$d "; \
	  else fail="$$fail$$d "; fi; \
	done; \
	for d in $(bundled_cmake_builddirs); do \
	  if $(MAKE) -C "$$d" test VERBOSE=$(filter-out 0,$(V)); then pass="$$pass$$d "; \
	  else fail="$$fail$$d "; fi; \
	done; \
	test -n "$$pass" && echo "PASS: $$pass"; \
	if test -n "$$fail"; then echo "FAIL: $$fail" && false; fi # exit code of the last command must be zero if $fail is empty

######################## COMPILATION AND LINKING RULES #########################

# Dependencies of the ICON executable:
$(ICON_exec_file): $(ICON_prog_obj_files) $(common_lib_file) $(BUNDLED_LIBFILES)

# Common linking rule for the executables:
$(exec_files): | $(dir_files)
	$(silent_FCLD)$(FC) -o $@ $(make_FCFLAGS) $(FCFLAGS) $(ICON_FCFLAGS) $(LDFLAGS) $+ $(BUNDLED_EXTRA_LIBS) $(LIBS)

# Dependencies of the common convenience library:
$(common_lib_file): $(common_lib_obj_files)

# Common archive generation rule:
$(lib_files): | $(dir_files)
	$(silent_AR)rm -f $@ && $(AR) $(ARFLAGS) $@ $^ && $(RANLIB) $@

# All objects that can be built with the inline library depend on it:
#   The compiler documentation suggests that we use $(pgi_inlib_name) as a
#   target and as a prerequisite but $(pgi_inlib_name) is a directory and using
#   directories as makefile targets is generally a bad idea. In particular, a
#   target directory is not removed on error, which means that if 'make' is
#   interrupted during the generation of the library, the subsequent 'make'
#   command will consider the library successfully created and erroneously
#   proceed with the next step. Therefore we use the table of contents of the
#   library, which is a regular file, as a target and as a dependency:
#$(pgi_inlib_target_obj_files): $(pgi_inlib_name)/TOC

# PGI/NVIDIA inline library generation rule (note that the source files are
# provided to the compiler in the topological order):
#$(pgi_inlib_name)/TOC: $(pgi_inlib_prereq_mod_files) $(pgi_inlib_obj_files:.o=.f90) inlib.mk | $(dir_files)
#	$(silent_FCEX)rm -rf $(@D) && $(FC) $(pgi_inlib_ex_FCFLAGS) $(make_FCFLAGS) $(FCFLAGS) $(ICON_FCFLAGS) -cpp $(filter %.f90,$^)

# Fortran compilation rule:
#%.o: %.f90 | $(dir_files) $(bundled_builddirs) sanitize-mod-proxies
#	$(silent_FC)/usr/bin/mkdir -p $(moddir_null)/$@ && $(FC) -o $@ -c -J$(moddir_null)/$@ $(make_FCFLAGS) $(FCFLAGS) $(ICON_FCFLAGS) -cpp $<

%.o: %.f90 | $(dir_files) $(bundled_builddirs) sanitize-mod-proxies
	$(silent_FC)$(FC) -o $@ -c -J$(moddir) $(make_FCFLAGS) $(FCFLAGS) $(ICON_FCFLAGS) -cpp $<

# Fortran module file generation rule:
#%$(pp_infix)modstamp.f90: %.f90 | $(dir_files) $(bundled_builddirs) sanitize-mod-proxies
#	$(silent_MOD)$(FC) -c -J$(moddir) $(make_FCFLAGS) $(FCFLAGS) $(ICON_FCFLAGS)  -cpp $< && touch $@

# Fortran module file tracking and recovery rule:
$(moddir)/%.mod.proxy:| sanitize-mod-proxies
	@if test -n '$<'; then \
	  if test ! -f '$(@:.proxy=)'; then rm -f '$<'; $(MAKE) -f icon.mk '$<'; fi; \
	  if cmp '$@' '$(@:.proxy=)' >/dev/null 2>&1 || $(MODCMP) '$@' '$(@:.proxy=)' gnu 2>/dev/null; then :; \
	  else cp '$(@:.proxy=)' '$@' 2>/dev/null; fi; \
	fi

# Fortran submodule file tracking and recovery rule:
$(moddir)/%.smod.sproxy:| sanitize-mod-proxies
	@if test -z '$<'; then \
	  echo "Cannot find Fortran source file providing submodule '$(basename $(@F:.sproxy=))'." >&2; \
	else \
	  if test ! -f '$(@:.sproxy=)'; then rm -f '$<'; $(MAKE) '$<'; fi; \
	  if cmp '$@' '$(@:.sproxy=)' >/dev/null 2>&1 || $(MODCMP) '$@' '$(@:.sproxy=)' gnu 2>/dev/null; then :; \
	  else cp '$(@:.sproxy=)' '$@' 2>/dev/null; fi; \
	fi

# Delete all Fortran (sub)module proxy files that do not have an existing
# (sub)module to be a proxy of, i.e. if <filename>.proxy exists but <filename>
# does not, delete <filename>.proxy:
sanitize-mod-proxies:
	@rm -f $(filter-out $(addsuffix .proxy,$(wildcard $(moddir)/*.mod)),$(wildcard $(moddir)/*.mod.proxy)) $(filter-out $(addsuffix .sproxy,$(wildcard $(moddir)/*.smod)),$(wildcard $(moddir)/*.smod.sproxy))

# C compilation rule:
%.o: %.c | $(dir_files)
	$(silent_CC)$(CC) -o $@ -c $(CFLAGS) $(CPPFLAGS) $<

# CUDA compilation rule:
#%.o: %.cu | $(dir_file)
#	$(silent_CUDACXX)$(CUDACXX) -o $@ -c $(CUDAFLAGS) $<

# HIP compilation rule:
#%.o: %.hip.cc | $(dir_file)
#	$(silent_HIPCXX)$(HIPCXX) -o $@ -c $(HIPFLAGS) $<

######################### DSL4JSB PREPROCESSING RULES ##########################

#%$(DSL4JSB_infix).f90: %.f90 | $(dir_files)
#	$(silent_DSL4JSB)$(DSL4JSB) -i $< -o $@

######################## DACE SUBST PREPROCESSING RULES ########################
.SECONDEXPANSION:
#%$(DACE_SUBST_infix).f90: \
  %.f90 \
  $(srcdir)/sdfgs/integrations.yaml \
  $(srcdir)/sdfgs/utils/dace-substf.py \
  $$(shell $$(DACE_PRINT_INTEGRATIONS) $$(srcdir)/sdfgs/integrations.yaml "association_files" $$*.f90) \
  | $(dir_files)
#	$(silent_DACE_SUBST)$(DACE_SUBST_PP) $(word 2,$^) $(wordlist 4,100000,$^) $< $@

# Don't delete intermediate files
# (`.NOTINTERMEDIATE` doesn't support pattern rules, so compute intermediate files explicitly)
#.NOTINTERMEDIATE: \
  $(foreach sdfg,$(shell $(DACE_PRINT_INTEGRATIONS) $(srcdir)/sdfgs/integrations.yaml "sdfg_names"),\
    sdfgs/build/$(sdfg)_bindings.f90 \
    sdfgs/build/lib$(sdfg).so \
    sdfgs/$(sdfg)_optimized.sdfgz \
    sdfgs/$(sdfg)_simplified.sdfgz \
    sdfgs/$(sdfg)_associations.yaml \
    sdfgs/$(sdfg)_module_definitions.yaml \
    sdfgs/$(sdfg)_unsimplified.sdfgz \
  )

# F90 binddings for SDFGs
# TODO: add unsimplified if it exists
#sdfgs/build/%_bindings.f90: \
  sdfgs/%_simplified.sdfgz \
  $(srcdir)/sdfgs/meta_data.yaml \
  sdfgs/%_module_definitions.yaml \
  sdfgs/build/lib%.so \
  $(srcdir)/sdfgs/utils/dace-genfi.py \
  | $(dir_files)
#	$(DACE_GENFI) $< $(word 2,$^) $(word 3,$^) $@

# Compiled SDFG
# We use symlinks instead of copies (e.g., `-o` option of `sdfgcc`) to avoid stale copies
# Need to touch the `.so` file, because DaCe won't regenerate it if source hasn't changed
#sdfgs/build/lib%.so: \
  sdfgs/%_optimized.sdfgz \
  .dace.conf \
  | $(dir_files)
#	sdfgcc $< && \
#	  touch $@ && \
#	  ln -sf ../../.dacecache/$*/build/lib$*.so ../../.dacecache/$*/include/$*.h sdfgs/build

# symlink `.dace.conf`
#.dace.conf: $(srcdir)/.dace.conf
#	ln -s $< $@

# Optimized SDFG
#sdfgs/%_optimized.sdfgz: \
  sdfgs/%_simplified.sdfgz \
  $(srcdir)/sdfgs/optimize_%.py \
  | $(dir_files)
#	$(srcdir)/sdfgs/optimize_$*.py $< $@

# Simplify
#sdfgs/%_simplified.sdfgz: \
  sdfgs/%_unsimplified.sdfgz \
  $(srcdir)/sdfgs/utils/dace-simplify.py \
  | $(dir_files)
#	$(srcdir)/sdfgs/utils/dace-simplify.py $< $@

# Generate unsimplified
#sdfgs/%_unsimplified.sdfgz \
  sdfgs/%_module_definitions.yaml \
  sdfgs/%_associations.yaml \
  &: \
    $(srcdir)/sdfgs/integrations.yaml \
    $(srcdir)/sdfgs/generate.py \
    | $(dir_files)
#	$(srcdir)/sdfgs/generate.py $* $< sdfgs

##################### EXPLICIT FORTRAN PREPROCESSING RULES #####################

%$(FPP_infix).f90: %.f90 | $(dir_files)
	$(silent_FPP)$(FPP) $(make_FCFLAGS) $(FCFLAGS) $(ICON_FCFLAGS) $< >$@

######################## SERIALBOX2 PREPROCESSING RULES ########################

#%$(SB2_infix).f90: %.f90  | $(dir_files)
#	$(silent_SB2)$(SB2PP) -o $@ $<

############################ BUNDLED PACKAGE RULES #############################

# We need to tolerate the absence of makefiles in the bundled packages when
# running clean rules. If the delayed configuration is disabled, we tolerate the
# absence of the makefiles for the 'distclean' target only. Otherwise, we
# additionally enable the tolerance for the 'mostlyclean' and 'clean' targets
# because they can be requested before the bundled packages are configured:
nomakefile_bundled_targets:= distclean
nomakefile_bundled_targets+= mostlyclean clean

# Make Autoconf-based bundled packages:
$(bundled_config_builddirs):
	@if test -f '$@/Makefile'; then \
	  $(MAKE) -C $@ $(filter all mostlyclean clean distclean,$(MAKECMDGOALS)) V=$(V); \
	else \
	  test -z '$(filter-out $(nomakefile_bundled_targets),$(or $(MAKECMDGOALS),all))'; \
	fi

# Make CMake-based bundled packages:
$(bundled_cmake_builddirs):
	@if test -f '$@/Makefile'; then \
	  $(MAKE) -C $@ $(filter all clean,$(MAKECMDGOALS) $(if $(filter mostlyclean distclean,$(MAKECMDGOALS)),clean)) VERBOSE=$(filter-out 0,$(V)); \
	else \
	  test -z '$(filter-out $(nomakefile_bundled_targets),$(or $(MAKECMDGOALS),all))'; \
	fi

# Build the bundled yaxt before configuring the bundled yac:
externals/yac/Makefile:|  $(filter externals/yaxt,$(bundled_builddirs))

# Build the bundled mtime and yaxt before building the bundled yac:
externals/yac: $(filter externals/mtime externals/yaxt,$(bundled_builddirs))

# Build the bundled yaxt and ppm before building the bundled cdi:
externals/cdi: $(filter externals/yaxt externals/ppm,$(bundled_builddirs))

# Build the bundled yac before the bundled hd:
externals/hd: $(filter externals/yac,$(bundled_builddirs))

# Configure the bundled fortran-support before configuring the bundled
# math-support:
externals/math-support/build/Makefile: $(filter externals/fortran-support/build/Makefile,$(bundled_makefiles))

# Build the bundled fortran-support before building the bundled-math support:
externals/math-support/build: $(filter externals/fortran-support/build,$(bundled_builddirs))

# Configure the bundled fortran-support and math-support before configuring the
# bundled math-interpolation:
externals/math-interpolation/build/Makefile: $(filter externals/fortran-support/build/Makefile externals/math-support/build/Makefile,$(bundled_makefiles))

# Build the bundled fortran-support and math-support before building the bundled
# math-interpolation:
externals/math-interpolation/build: $(filter externals/fortran-support/build externals/math-support/build,$(bundled_builddirs))

# Configure delayed bundled packages:
export FC CC CXX

externals/ecrad: $(if $(filter-out mostlyclean clean distclean,$(or $(MAKECMDGOALS),all)),externals/ecrad/Makefile)
externals/ecrad/Makefile: icon.mk
	$(silent_CONFIG)$(MKDIR_P) $(@D) && cd $(@D) && '/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/ecrad/configure'  'CPPFLAGS=-I/usr/include -I/usr/include/libxml2' 'CXX=g++' 'CXXFLAGS=-g -march=native' 'LDFLAGS=-L/usr/lib64 -L/usr/lib' 'LIBS=-lxml2 -lfyaml -leccodes -llapack -lblas -lnetcdff -lnetcdf -lstdc++' 'MPI_LAUNCH=mpiexec' '--enable-openmp' '--enable-grib2' '--enable-mixed-precision' '--enable-rte-rrtmgp' 'CC=mpicc' 'FC=mpif90' 'LDFLAGS=-L/usr/lib/x86_64-linux-gnu/ -L/home/primrose/Work/IconGrounds/icon-dace/.dacecache/radiation/build' 'LIBS=-leccodes -lnetcdff -lnetcdf -lopenblas -lradiation' '--enable-acm-license' '--disable-mixed-precision' '--disable-edmf' '--disable-les' '--disable-ocean' '--disable-jsbach' '--disable-coupling' '--disable-aes' '--disable-rte-rrtmgp' '--enable-ecrad' '--disable-mpi' '--disable-mpi-checks' '--disable-openmp' '--disable-loop-exchange' '--enable-dace-subst=no' '--enable-explicit-fpp' '--disable-option-checking' '--prefix=/usr/local' '--cache-file=/dev/null' '--srcdir=/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/ecrad' 'FCFLAGS=-g -O2 -I/usr/include/ -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-missing-include-dirs -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE -std=legacy' 'NETCDF_FCFLAGS=' 'NETCDF_FCLIBS=' '--enable-silent-rules=yes' '--enable-pgi-inlib=no' $(if $(silent_CONFIG),--silent)

externals/cdi: $(if $(filter-out mostlyclean clean distclean,$(or $(MAKECMDGOALS),all)),externals/cdi/Makefile)
externals/cdi/Makefile: icon.mk
	$(silent_CONFIG)$(MKDIR_P) $(@D) && cd $(@D) && '/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/cdi/configure'  'CPPFLAGS=-I/usr/include -I/usr/include/libxml2' 'CXXFLAGS=-g -march=native' 'LDFLAGS=-L/usr/lib64 -L/usr/lib' 'LIBS=-lxml2 -lfyaml -leccodes -llapack -lblas -lnetcdff -lnetcdf -lstdc++' '--enable-grib2' '--enable-mixed-precision' '--enable-rte-rrtmgp' 'CC=mpicc' 'FC=mpif90' 'LDFLAGS=-L/usr/lib/x86_64-linux-gnu/ -L/home/primrose/Work/IconGrounds/icon-dace/.dacecache/radiation/build' 'LIBS=-leccodes -lnetcdff -lnetcdf -lopenblas -lradiation' '--enable-acm-license' '--disable-mixed-precision' '--disable-edmf' '--disable-les' '--disable-ocean' '--disable-jsbach' '--disable-coupling' '--disable-aes' '--disable-rte-rrtmgp' '--enable-ecrad' '--disable-mpi-checks' '--disable-loop-exchange' '--enable-dace-subst=no' '--enable-explicit-fpp' '--disable-option-checking' '--prefix=/usr/local' '--cache-file=/dev/null' '--srcdir=/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/cdi' 'CFLAGS=-g -march=native -O2' 'FCFLAGS=-g -O2 -I/usr/include/ -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-missing-include-dirs -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE -std=legacy' 'F77=no' 'CXX=no' 'UTIL_LINUX_UUID_C_INCLUDE=' 'UTIL_LINUX_UUID_C_LIB=' 'OSSP_UUID_C_INCLUDE=' 'OSSP_UUID_C_LIB=' 'DCE_UUID_C_INCLUDE=' 'DCE_UUID_C_LIB=' 'MPIROOT=' 'MPI_C_INCLUDE=' 'MPI_C_LIB=' 'MPI_LAUNCH=/usr/bin//mpiexec' 'MPI_FC_MOD=' 'MPI_FC_LIB=' 'PKG_CONFIG=' 'YAXT_C_INCLUDE=' 'YAXT_C_LIB=' 'YAXT_FC_MOD=' 'YAXT_FC_LIB=' 'PPM_CORE_C_INCLUDE=' 'PPM_CORE_C_LIB=' 'BUILD_CFLAGS=' 'BUILD_FCFLAGS=' 'BUILD_LDFLAGS=' 'BUILD_LIBS=' 'BUILD_MPI_C_LIB=' 'BUILD_MPI_FC_LIB=' 'BUILD_CC=' 'BUILD_CXX=' 'BUILD_FC=' 'BUILD_F77=' '--enable-silent-rules=yes' '--disable-maintainer-mode' '--enable-mpi=no' '--enable-iso-c-interface' '--enable-cf-interface=no' '--disable-ruby-interface' '--disable-python-interface' '--disable-openmp' '--disable-shared' '--enable-static' '--enable-grib' '--disable-across' '--enable-cgribex' '--disable-cdi-app' '--enable-ppm-dist-array=no' '--without-util-linux-uuid' '--without-ossp-uuid' '--without-dce-uuid' '--without-threads' '--without-fdb5' '--without-szlib' '--with-netcdf' '--with-eccodes=yes' '--without-grib_api' '--with-on-demand-check-programs' '--without-example-programs' 'acx_cv_have_netcdf2=yes' 'acx_cv_have_netcdf4=yes' 'acx_cv_have_pnetcdf=no' 'acx_cv_have_libnc_dap=no' 'acx_cv_have_nc4hdf5=no' $(if $(silent_CONFIG),--silent)

externals/mtime: $(if $(filter-out mostlyclean clean distclean,$(or $(MAKECMDGOALS),all)),externals/mtime/Makefile)
externals/mtime/Makefile: icon.mk
	$(silent_CONFIG)$(MKDIR_P) $(@D) && cd $(@D) && '/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/mtime/configure'  'CPPFLAGS=-I/usr/include -I/usr/include/libxml2' 'CXX=g++' 'CXXFLAGS=-g -march=native' 'LDFLAGS=-L/usr/lib64 -L/usr/lib' 'LIBS=-lxml2 -lfyaml -leccodes -llapack -lblas -lnetcdff -lnetcdf -lstdc++' 'MPI_LAUNCH=mpiexec' '--enable-openmp' '--enable-grib2' '--enable-mixed-precision' '--enable-rte-rrtmgp' 'CC=mpicc' 'FC=mpif90' 'LDFLAGS=-L/usr/lib/x86_64-linux-gnu/ -L/home/primrose/Work/IconGrounds/icon-dace/.dacecache/radiation/build' 'LIBS=-leccodes -lnetcdff -lnetcdf -lopenblas -lradiation' '--enable-acm-license' '--disable-mixed-precision' '--disable-edmf' '--disable-les' '--disable-ocean' '--disable-jsbach' '--disable-coupling' '--disable-aes' '--disable-rte-rrtmgp' '--enable-ecrad' '--disable-mpi' '--disable-mpi-checks' '--disable-openmp' '--disable-loop-exchange' '--enable-dace-subst=no' '--enable-explicit-fpp' '--disable-option-checking' '--prefix=/usr/local' '--cache-file=/dev/null' '--srcdir=/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/mtime' 'CFLAGS=-g -march=native -O2' 'FCFLAGS=-g -O2 -I/usr/include/ -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-missing-include-dirs -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE -std=legacy' '--enable-silent-rules=yes' '--disable-maintainer-mode' '--enable-shared=no' '--enable-static' '--disable-examples' '--disable-fortran-hl' '--enable-python=no' '--with-pic=no' $(if $(silent_CONFIG),--silent)

externals/fortran-support/build: $(if $(filter-out mostlyclean clean distclean,$(or $(MAKECMDGOALS),all)),externals/fortran-support/build/Makefile)
externals/fortran-support/build/Makefile: icon.mk
	$(silent_CMAKE)$(MKDIR_P) $(@D) && cd $(@D) && "${CMAKE}"  '-Wno-dev' '--no-warn-unused-cli' '-GUnix Makefiles' '-DCMAKE_EXE_LINKER_FLAGS=-L/usr/lib/x86_64-linux-gnu/ -L/home/primrose/Work/IconGrounds/icon-dace/.dacecache/radiation/build' '-DCMAKE_C_STANDARD_LIBRARIES=-leccodes -lnetcdff -lnetcdf -lopenblas -lradiation' '-DCMAKE_Fortran_STANDARD_LIBRARIES=-leccodes -lnetcdff -lnetcdf -lopenblas -lradiation' '-DCMAKE_EXE_LINKER_FLAGS=' '-DCMAKE_C_STANDARD_LIBRARIES=' '-DCMAKE_Fortran_STANDARD_LIBRARIES=' '-DCMAKE_INSTALL_PREFIX=/usr/local' '/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/fortran-support' '-DCMAKE_C_FLAGS=-g -march=native -O2 -I/usr/include -I/usr/include/libxml2' '-DCMAKE_Fortran_FLAGS=-g -O2 -I/usr/include/ -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-missing-include-dirs -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE -std=legacy' '-DFS_ENABLE_OMP=no' '-DFS_ENABLE_MIXED_PRECISION=no' '-DFS_ENABLE_OPENACC=OFF' '-DBUILD_SHARED_LIBS:BOOL=OFF' '-DBUILD_TESTING:BOOL=OFF' '-DCMAKE_BUILD_TYPE:STRING=NOEXTRAFLAGS' $(if $(silent_CMAKE),>/dev/null)

externals/math-support/build: $(if $(filter-out mostlyclean clean distclean,$(or $(MAKECMDGOALS),all)),externals/math-support/build/Makefile)
externals/math-support/build/Makefile: icon.mk
	$(silent_CMAKE)$(MKDIR_P) $(@D) && cd $(@D) && "${CMAKE}"  '-Wno-dev' '--no-warn-unused-cli' '-GUnix Makefiles' '-DCMAKE_EXE_LINKER_FLAGS=-L/usr/lib/x86_64-linux-gnu/ -L/home/primrose/Work/IconGrounds/icon-dace/.dacecache/radiation/build' '-DCMAKE_Fortran_STANDARD_LIBRARIES=-leccodes -lnetcdff -lnetcdf -lopenblas -lradiation' '-DCMAKE_EXE_LINKER_FLAGS=' '-DCMAKE_Fortran_STANDARD_LIBRARIES=' '-DCMAKE_INSTALL_PREFIX=/usr/local' '/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/math-support' '-DCMAKE_Fortran_FLAGS=-g -O2 -I/usr/include/ -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-missing-include-dirs -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE -std=legacy' '-DMS_ENABLE_LOOP_EXCHANGE=no' '-DMS_ENABLE_MIXED_PRECISION=no' '-DMS_ENABLE_DIM_SWAP=no' '-DMS_ENABLE_OPENACC=OFF' '-DCMAKE_PREFIX_PATH:PATH=/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/build/verification/externals/fortran-support/build' '-DBUILD_SHARED_LIBS:BOOL=OFF' '-DBUILD_TESTING:BOOL=OFF' '-DBUILD_STANDALONE:BOOL=OFF' '-DCMAKE_BUILD_TYPE:STRING=NOEXTRAFLAGS' $(if $(silent_CMAKE),>/dev/null)

externals/math-interpolation/build: $(if $(filter-out mostlyclean clean distclean,$(or $(MAKECMDGOALS),all)),externals/math-interpolation/build/Makefile)
externals/math-interpolation/build/Makefile: icon.mk
	$(silent_CMAKE)$(MKDIR_P) $(@D) && cd $(@D) && "${CMAKE}"  '-Wno-dev' '--no-warn-unused-cli' '-GUnix Makefiles' '-DCMAKE_EXE_LINKER_FLAGS=-L/usr/lib/x86_64-linux-gnu/ -L/home/primrose/Work/IconGrounds/icon-dace/.dacecache/radiation/build' '-DCMAKE_Fortran_STANDARD_LIBRARIES=-leccodes -lnetcdff -lnetcdf -lopenblas -lradiation' '-DCMAKE_EXE_LINKER_FLAGS=' '-DCMAKE_Fortran_STANDARD_LIBRARIES=' '-DCMAKE_INSTALL_PREFIX=/usr/local' '/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/math-interpolation' '-DCMAKE_Fortran_FLAGS=-g -O2 -I/usr/include/ -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-missing-include-dirs -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE -std=legacy' '-DMI_ENABLE_LOOP_EXCHANGE=no' '-DMI_ENABLE_MIXED_PRECISION=no' '-DMI_ENABLE_DIM_SWAP=no' '-DMI_ENABLE_OPENACC=OFF' '-DCMAKE_PREFIX_PATH:PATH=/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/build/verification/externals/fortran-support/build;/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/build/verification/externals/math-support/build' '-DBUILD_SHARED_LIBS:BOOL=OFF' '-DBUILD_TESTING:BOOL=OFF' '-DBUILD_STANDALONE:BOOL=OFF' '-DCMAKE_BUILD_TYPE:STRING=NOEXTRAFLAGS' $(if $(silent_CMAKE),>/dev/null)

# Relink executables if any of the bundled libraries is updated (the semicolon
# is required to support parallel rebuild):
$(BUNDLED_LIBFILES): $(bundled_builddirs);

############################### AUXILIARY RULES ################################

# Version file generation rules:
version.o: config.h
version.c: force-create-version
	$(silent_GEN)$(PVCS) --srcdir $(srcdir) --output $@ --subdirs $(PVCS_subdirs)

# Directory creation rule:
%/.dirstamp:
	$(silent_MKDIR)$(MKDIR_P) $(@D) && touch $@

# Keep directory stamps:
.PRECIOUS: $(dir_files)

######################### DEPENDENCY GENERATION RULES ##########################

# Fortran dependency generation rule:
%.f90.d: %.f90 icon.mk depgen.f90.config | $(dir_files)
	$(silent_DEPGEN)$(DEPGEN) @depgen.f90.config -o $@ --obj-name $(@:.f90.d=.o) $(if ,--fc-mod-stamp-name $(@:.f90.d=$(pp_infix)modstamp.f90)) -i $< -- $(DEPGEN_FCFLAGS) -J$(moddir) $(make_FCFLAGS) $(FCFLAGS) $(ICON_FCFLAGS)

# C dependency generation rule:
%.c.d: %.c icon.mk depgen.c.config | $(dir_files)
	$(silent_DEPGEN)$(DEPGEN) @depgen.c.config -o $@ --obj-name $(@:.c.d=.o) -i $< -- $(CFLAGS) $(CPPFLAGS)

# CUDA dependency generation rule:
#%.cu.d: %.cu icon.mk depgen.c.config | $(dir_files)
#	$(silent_DEPGEN)$(DEPGEN) @depgen.c.config --pp-inc-sys -o $@ --obj-name $(@:.cu.d=.o) -i $< -- $(CUDAFLAGS)

# HIP dependency generation rule:
#%.hip.cc.d: %.hip.cc icon.mk depgen.c.config | $(dir_files)
#	$(silent_DEPGEN)$(DEPGEN) @depgen.c.config --pp-inc-sys -o $@ --obj-name $(@:.hip.cc.d=.o) -i $< -- $(HIPFLAGS)

# Dependency generation rule for Fortran-to-C bindings:
c_binding.d: icon.mk
	$(silent_DEPGEN):;{ \
	  echo '$(call get_obj_names,src/io/restart/mo_c_restart_util.f90): #-hint support/util_multifile_restart.o'; \
	  echo '$(call get_obj_names,support/mo_util_uuid.f90): #-hint support/util_uuid.o'; \
	  echo '#$(call get_obj_names,src/shared/mo_index_list.f90): #-hint support/index_list.o'; \
	  echo '#$(call get_obj_names,$(DACE_subdir)/mo_wigos.f90): #-hint $(DACE_subdir)/dace_md5_ifc.o $(DACE_subdir)/md5.o'; \
	} >$@

# Dummy dependency file generation rule (called by config.status):
dummy-depend: | $(dir_files)
	@for file in $(dep_files); do \
	  test -e "$$file" || touch "$$file"; \
	done

# Include dependencies if required:
ifneq (,$(filter-out $(nodep_targets),$(or $(MAKECMDGOALS),all)))
include $(dep_files)
endif
