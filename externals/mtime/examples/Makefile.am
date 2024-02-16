# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
if ENABLE_EXAMPLES
noinst_PROGRAMS =     \
    callback_test     \
    comp_weights      \
    duration          \
    example           \
    example_hl        \
    iconatm           \
    iconoce           \
    iconoce_hl        \
    model_integration \
    modulo            \
    output_control    \
    recurrence        \
    repetitor         \
    simulate_iau      \
    tas               \
    test_dace         \
    time_calculus     \
    test_jd_logic     \
    uniq
endif ENABLE_EXAMPLES

AM_CPPFLAGS = -I$(top_srcdir)/include
AM_FCFLAGS = $(FCMODINC)$(top_builddir)/src
LDADD = $(top_builddir)/src/libmtime.la

AM_DEFAULT_SOURCE_EXT = .f90
duration_SOURCES = duration.f90 mo_kind.f90
iconatm_SOURCES = iconatm.f90 mo_event_manager.f90
iconoce_SOURCES = iconoce.f90 mo_event_manager.f90
model_integration_SOURCES = model_integration.c
modulo_SOURCES = modulo.c
recurrence_SOURCES = recurrence.c
test_dace_SOURCES = test_dace.f90 mo_event_manager.f90
time_calculus_SOURCES = time_calculus.c
test_jd_logic_SOURCES = test_jd_logic.f90

if FCMODUC
mo_event_manager_mod = MO_EVENT_MANAGER.$(FCMODEXT)
mo_exception_mod = MO_EXCEPTION.$(FCMODEXT)
mo_kind_mod = MO_KIND.$(FCMODEXT)
mtime_error_handling_mod = MTIME_ERROR_HANDLING.$(FCMODEXT)
mtime_hl_mod = MTIME_HL.$(FCMODEXT)
mtime_mod = MTIME.$(FCMODEXT)
else !FCMODUC
mo_event_manager_mod = mo_event_manager.$(FCMODEXT)
mo_exception_mod = mo_exception.$(FCMODEXT)
mo_kind_mod = mo_kind.$(FCMODEXT)
mtime_error_handling_mod = mtime_error_handling.$(FCMODEXT)
mtime_hl_mod = mtime_hl.$(FCMODEXT)
mtime_mod = mtime.$(FCMODEXT)
endif !FCMODUC

$(mo_event_manager_mod): mo_event_manager.$(OBJEXT)
$(mo_kind_mod): mo_kind.$(OBJEXT)

$(mo_event_manager_mod) $(mo_kind_mod):
	@if test -z $<; then \
	  echo "Cannot find Fortran source file providing module '$(basename $(@F))'." >&2; \
	elif test ! -f $@; then \
	  rm -f $<; $(MAKE) $(AM_MAKEFLAGS) $<; \
	fi

callback_test.$(OBJEXT)    \
comp_weights.$(OBJEXT)     \
duration.$(OBJEXT)         \
example.$(OBJEXT)          \
iconatm.$(OBJEXT)          \
iconoce.$(OBJEXT)          \
iconoce_hl.$(OBJEXT)       \
mo_event_manager.$(OBJEXT) \
output_control.$(OBJEXT)   \
repetitor.$(OBJEXT)        \
simulate_iau.$(OBJEXT)     \
tas.$(OBJEXT)              \
test_dace.$(OBJEXT)        \
test_jd_logic.$(OBJEXT)    \
uniq.$(OBJEXT): $(top_builddir)/src/$(mtime_mod)

example_hl.$(OBJEXT)   \
iconoce_hl.$(OBJEXT)   \
simulate_iau.$(OBJEXT): $(top_builddir)/src/$(mtime_hl_mod)

duration.$(OBJEXT):                                 \
    $(mo_kind_mod)                                  \
    $(top_builddir)/src/$(mtime_error_handling_mod)

iconatm.$(OBJEXT)   \
iconoce.$(OBJEXT)   \
test_dace.$(OBJEXT): $(mo_event_manager_mod)

clean-local:
	rm -rf *.dSYM

MOSTLYCLEANFILES = $(mo_event_manager_mod) $(mo_exception_mod) $(mo_kind_mod)

mostlyclean-local:
	rm -f *.i
