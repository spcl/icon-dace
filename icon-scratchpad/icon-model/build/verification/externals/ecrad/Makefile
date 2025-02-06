# Path to the directory with the source files:
srcdir:= /home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/ecrad

# Path to the directory with the Fortran module files:
moddir:= mod

# Paths to the installation directories:
prefix= /usr/local
exec_prefix= ${prefix}

# Package tarname:
PACKAGE_TARNAME= ecrad-1.1.0

# Utilities:
SHELL= /bin/bash
FC= mpif90
PYTHON= python
DEPLIST= ${PYTHON} $(srcdir)/mkhelper/deplist.py
DEPGEN= ${PYTHON} $(srcdir)/mkhelper/depgen.py --src-root='/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/ecrad' --pp-enable --pp-eval-expr \
        --pp-inc-flag='-I' --pp-inc-order='inc,flg' --pp-macro-flag='-D' \
        --fc-enable --fc-mod-ext='mod.proxy' --fc-mod-upper='no' --fc-inc-flag='-I' \
        --fc-inc-order='src,flg' --fc-mod-dir-flag='-J' --fc-external-mods='netcdf'
MODCMP= ${PYTHON} $(srcdir)/mkhelper/fortmodcmp.py
AR= ar
ARFLAGS= cr
INSTALL= /usr/bin/install -c
INSTALL_DATA= ${INSTALL} -m 644
GIT= git
TAR= tar
BZIP2 = bzip2

# Fortran compiler and flags:
FCFLAGS= -I$(moddir) -J$(moddir) -I$(srcdir)/include  -g -O2 -I/usr/include/ -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-missing-include-dirs -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE -std=legacy -D__ECRAD_LITTLE_ENDIAN
LDFLAGS= -L/usr/lib/x86_64-linux-gnu/ -L/home/primrose/Work/IconGrounds/icon-dace/.dacecache/radiation/build
LIBS=  -leccodes -lnetcdff -lnetcdf -lopenblas -lradiation

# Silent rule prefixes:
V= 0
ifeq ($(V), 0)
silent_FC=      @echo "  FC      " $@;
silent_MKDIR=   @echo "  MKDIR   " $(@D);
silent_DEPGEN=  @echo "  DEPGEN  " $@;
silent_AR=      @echo "  AR      " $@;
endif

lib_LIBRARIES:= libifsaux.a libifsrrtm.a libradiation.a libutilities.a

ignored_SOURCES:=                                \
    ifsrrtm/rrtm_ecrt_140gp_mcica.F90            \
    ifsrrtm/rrtm_rrtm_140gp_mcica.F90            \
    ifsrrtm/srtm_gas_optical_depth_test.F90      \
    ifsrrtm/srtm_spcvrt_mcica.F90                \
    ifsrrtm/srtm_srtm_224gp_mcica.F90            \
    radiation/radiation_adding_ica_sw_test.F90   \
    radiation/radiation_adding_ica_sw_test2.F90  \
    radiation/radiation_ice_optics_baran2016.F90 \
    radiation/radiation_psrad.F90                \
    radiation/radiation_psrad_rrtm.F90

libifsaux_a_SOURCES:= $(filter-out $(ignored_SOURCES),$(patsubst $(srcdir)/%,%,$(shell find $(srcdir)/ifsaux -name '*.F90')))
libifsrrtm_a_SOURCES:= $(filter-out $(ignored_SOURCES),$(patsubst $(srcdir)/%,%,$(shell find $(srcdir)/ifsrrtm -name '*.F90')))
libradiation_a_SOURCES:= $(filter-out $(ignored_SOURCES),$(patsubst $(srcdir)/%,%,$(shell find $(srcdir)/radiation -name '*.F90')))
libutilities_a_SOURCES:= $(filter-out $(ignored_SOURCES),$(patsubst $(srcdir)/%,%,$(shell find $(srcdir)/utilities -name '*.F90')))

all_SOURCES:= $(libifsaux_a_SOURCES) $(libifsrrtm_a_SOURCES) $(libradiation_a_SOURCES) $(libutilities_a_SOURCES)

# Dependency files:
deps:= $(addsuffix .d,$(all_SOURCES))

# Stamp files of the building subdirectories:
dirstamps= $(addsuffix .dirstamp, $(sort $(dir $(deps)))) $(moddir)/.dirstamp

# PGI cross-file function inlining via an inline library:
INLIB= 0
ifeq ($(INLIB), 1)

# Additional silence prefixes:
ifeq ($(V), 0)
silent_FCEX= @echo "  FC (EX) " $@;
silent_FCIN= @echo "  FC (IL) " $@;
endif

# Name of the inline library (a directory):
inlib_name:= ecrad.il

# Additional compiler flags enabling the inline library generation:
inlib_extract_FCFLAGS:= -Mextract=lib:$(inlib_name),reshape,name:adding_ica_sw,name:beta2alpha,name:calc_fluxes_no_scattering_lw,name:calc_ice_optics_baran,name:calc_ice_optics_baran2016,name:calc_ice_optics_baran2017,name:calc_ice_optics_fu_lw,name:calc_ice_optics_fu_sw,name:calc_ice_optics_yi_lw,name:calc_ice_optics_yi_sw,name:calc_liq_optics_lindner_li,name:calc_liq_optics_slingo,name:calc_liq_optics_socrates,name:calc_no_scattering_transmittance_lw,name:calc_ref_trans_lw,name:calc_ref_trans_sw,name:calc_two_stream_gammas_lw,name:calc_two_stream_gammas_sw,name:cloud_generator_acc,name:delta_eddington_extensive,name:delta_eddington_scat_od,name:fast_adding_ica_lw,name:initialize_acc,name:radiation_liquid_optics_socrates,name:uniform_distribution_acc


# Additional compiler flags enabling the inline library usage:
inlib_inline_FCFLAGS:= -Minline=lib:$(inlib_name),reshape

# List of source files containing functions that need to be inlined:
inlib_SOURCES:=                    \
		radiation/radiation_adding_ica_lw.F90 \
		radiation/radiation_adding_ica_sw.F90 \
		radiation/radiation_aerosol_optics.F90 \
		radiation/radiation_cloud_generator_acc.F90 \
		radiation/radiation_cloud_optics.F90 \
		radiation/radiation_ice_optics_baran.F90 \
		radiation/radiation_ice_optics_baran2017.F90 \
		radiation/radiation_ice_optics_fu.F90 \
		radiation/radiation_ice_optics_yi.F90 \
		radiation/radiation_liquid_optics_slingo.F90 \
		radiation/radiation_liquid_optics_socrates.F90 \
		radiation/radiation_two_stream.F90 \
		radiation/radiation_random_numbers.F90

		#radiation/radiation_cloud_cover.F90 \ # leads to cyclic dependency, missing inlinling is acceptable for now, affected: beta2alpha

# List of files where we need to manually limit GPU register usage
inlib_limit_registers_OBJECTS:=            \
		radiation/radiation_mcica_acc_lw.o \
		radiation/radiation_mcica_acc_sw.o

# And limiting flags
inlib_limit_registers_FCFLAGS:= -gpu=maxregcount:96

########################################

# Object files that correspond to the source files in $(inlib_SOURCES):
inlib_objs:= $(inlib_SOURCES:.F90=.o)

# List of objects that are generated using the inline library:
inlib_target_objs:= $(filter %.o,$(shell $(DEPLIST) -p $(inlib_objs) -f $(deps)))

# To avoid circular dependencies, we need to account for situations when we
# have a dependency A -> B -> C, where A and C belong to $(inlib_objs) but B
# does not. In order to get all such B files and include them into
# $(inlib_objs), we find the intersection of two sets: dependencies
# (prerequisites) and dependents (targets) of $(inlib_objs):
inlib_objs:= $(filter $(inlib_target_objs),$(shell $(DEPLIST) -t $(inlib_objs) -f $(deps)))

# The commented code below generates a warning message if the list of source
# files to be included into the inline library is automatically extended:
#
# inlib_extra_objs = $(filter-out $(inlib_SOURCES:.F90=.o),$(inlib_objs))
# ifneq ($(inlib_extra_objs),)
# $(warning Additional files to be included into the inline library: $(inlib_extra_objs:.o=.F90))
# endif

# Now we need to create a list of Fortran module files that must be created
# before any of the object files $(inlib_objs) gets generated. First, we need a
# subset of dependency files that correspond to $(inlib_objs):
inlib_deps:= $(inlib_objs:.o=.F90.d)

# List of Fortran module files that at least on element of $(inlib_objs)
# depends on:
inlib_mods:= $(filter %.mod.proxy,$(shell $(DEPLIST) -f $(inlib_deps)))

# List of Fortran modules declared in source files that correspond to
# $(inlib_objs):
inlib_internal_mods:= $(filter %.mod.proxy,$(shell $(DEPLIST) -p $(inlib_objs) -f $(inlib_deps)))

# List of Fortran modules that must be created before the inline library:
inlib_prereq_mods:= $(filter-out $(inlib_internal_mods),$(inlib_mods))

endif # INLIB == 1

# Selective search path:
vpath %.F90 $(srcdir)

# Disable built-in suffix rules:
.SUFFIXES:
# Targets not associated with files:
.PHONY: all depend dummy-depend mostlyclean clean distclean \
        maintainer-clean install install-libs install-mods dist \
        check sanitize-mod-proxies
# Targets that do not need the inclusion of the dependency files:
NO_INC_TARGETS:= depend dummy-depend mostlyclean clean distclean \
                 maintainer-clean dist
# Keep directory stamps:
.PRECIOUS: $(dirstamps)

# Default rule:
all: $(lib_LIBRARIES)

# Explicit dependency generation rule:
depend: $(deps)

# Delete the results of compilation and linking:
mostlyclean: $(bundled_subdirs)
	rm -f $(lib_LIBRARIES) $(all_SOURCES:.F90=.o)
	rm -f $(moddir)/*.mod $(moddir)/*.mod.proxy
	rm -f $(notdir $(all_SOURCES:.F90=.i))
	test x'$(INLIB)' != x1 || rm -rf $(inlib_name)

# Delete files generated at the building stage:
clean: mostlyclean

# Delete everything generated at the configure stage (and clean the created directories if they are empty):
distclean: clean
	find . -name '*.pyc' -delete -o -name '*.pyo' -delete -o -name '__pycache__' -delete
	rm -f $(deps) $(dirstamps)
	@for dir in ifsaux ifsrrtm radiation utilities $(moddir); do \
	  if test -d "$$dir"; then \
	    echo "find '$$dir' -type d -empty -delete"; \
	    find "$$dir" -type d -empty -delete; \
	  fi; \
	done
	rm -f config.log config.status Makefile

# Delete files generated at the autoreconf stage:
maintainer-clean: distclean
	rm -rf autom4te.cache

# Installation rules:
install: all install-libs install-mods

# Check rule:
check: all

# Tarball creation rule:
dist:
	@if test ! -e /home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/ecrad/.git; then echo "'/home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/ecrad' is not a git repository" >&2; exit 1; fi
	$(GIT) -C /home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/externals/ecrad archive --prefix=$(PACKAGE_TARNAME)/ --format tar -o /home/primrose/Work/IconGrounds/icon-dace/icon-scratchpad/icon-model/build/verification/externals/ecrad/$(PACKAGE_TARNAME).tar HEAD
	rm -f $(PACKAGE_TARNAME).tar.bz2 && BZIP2=$${BZIP2--9} $(BZIP2) $(PACKAGE_TARNAME).tar

ifeq ($(INLIB), 1)

# Pattern- and target-specific assignments are propagated to the prerequisites
# and override the global assignments. Therefore, we introduce the following
# match-anything pattern assignments to prevent that (the eval/value combination
# is required to keep the original global values of the variables without
# changing their flavors, i.e. keep them recursively expanded):
${eval %: silent_FC= $(value silent_FC)}
${eval %: FCFLAGS= $(value FCFLAGS)}

# Target-specific variables for objects that use the inline library:
$(inlib_target_objs): silent_FC= $(silent_FCIN)
$(inlib_target_objs): FCFLAGS+= $(inlib_inline_FCFLAGS)
$(inlib_limit_registers_OBJECTS): FCFLAGS+= $(inlib_limit_registers_FCFLAGS)

# All object that can be built with the inline library depend on it:
$(inlib_target_objs): $(inlib_name)

# The inline library generation rule. Note that the source files are provided
# to the compiler in the topological order. We also have to delete the partially
# generated library if the compiler fails:
$(inlib_name): $(inlib_prereq_mods) $(inlib_objs:.o=.F90)
	$(silent_FCEX)rm -rf $@ && $(FC) $(FCFLAGS) $(inlib_extract_FCFLAGS)  $(filter-out $(inlib_prereq_mods),$^) || (ec=$$?; rm -rf $@; exit $$ec)

endif # INLIB == 1

# Fortran compilation rule:
%.o: %.F90 | $(dirstamps)
	$(silent_FC)$(FC) -o $@ -c $(FCFLAGS)  $<

# Static library generation rule:
%.a:
	$(silent_AR)rm -f $@ && $(AR) $(ARFLAGS) $@ $^

# Fortran module file rule:
$(moddir)/%.mod.proxy:| sanitize-mod-proxies
	@if test -z '$<'; then \
	  echo "Cannot find Fortran source file providing module '$(basename $(@F:.proxy=))'." >&2; \
	else \
	  if test ! -f '$(@:.proxy=)'; then rm -f '$<'; $(MAKE) '$<'; fi; \
	  if cmp '$@' '$(@:.proxy=)' >/dev/null 2>&1 || $(MODCMP) '$@' '$(@:.proxy=)' gnu 2>/dev/null; then :; \
	  else cp '$(@:.proxy=)' '$@' 2>/dev/null; fi; \
	fi

# Delete all Fortran module proxy files that do not have an existing module to
# be a proxy of, i.e. if <filename>.proxy exists but <filename> does not,
# delete <filename>.proxy:
sanitize-mod-proxies:
	@rm -f $(filter-out $(addsuffix .proxy,$(wildcard $(moddir)/*.mod)),$(wildcard $(moddir)/*.mod.proxy))

# Directory creation rule:
%/.dirstamp:
	$(silent_MKDIR)/usr/bin/mkdir -p $(@D) && touch $@

# Fortran dependency generation rule:
%.F90.d: %.F90 Makefile | $(dirstamps)
	$(silent_DEPGEN)$(DEPGEN) -o $@ --obj-name $(@:.F90.d=.o) -i $< -- $(FCFLAGS)

# Dummy dependency file generation rule (called by config.status): 
dummy-depend: | $(dirstamps)
	@for file in $(deps); do \
	  test -e "$$file" || touch "$$file"; \
	done

# Library installation rule:
install-libs: $(lib_LIBRARIES)
	$(INSTALL) -d $(DESTDIR)${exec_prefix}/lib && $(INSTALL_DATA) $^ $(DESTDIR)${exec_prefix}/lib

# Fortran module files installation rule:
install-mods: $(filter %.mod.proxy,$(shell $(DEPLIST) -f $(deps)))
	$(INSTALL) -d $(DESTDIR)${prefix}/include && $(INSTALL_DATA) $(basename $^) $(DESTDIR)${prefix}/include

libifsaux.a: $(libifsaux_a_SOURCES:.F90=.o)
libifsrrtm.a: $(libifsrrtm_a_SOURCES:.F90=.o)
libradiation.a: $(libradiation_a_SOURCES:.F90=.o)
libutilities.a: $(libutilities_a_SOURCES:.F90=.o)

current_targets:= $(strip $(MAKECMDGOALS))
ifeq (,$(current_targets))
current_targets:= all
endif

ifneq (,$(filter-out $(NO_INC_TARGETS),$(current_targets)))
include $(deps)
endif
