#
# Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
#
# Redistribution and use in source and binary forms, with or without
# modification, are  permitted provided that the following conditions are
# met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# Neither the name of the DKRZ GmbH nor the names of its contributors
# may be used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
EXTRA_PROGRAMS = \
	check_prng

if with_fortran
EXTRA_PROGRAMS += \
	show_ieee_emulation \
	check_prng_f \
	graph_build
endif

noinst_PROGRAMS =
if with_example_programs
noinst_PROGRAMS += $(EXTRA_PROGRAMS)
endif
examples: $(EXTRA_PROGRAMS)

if with_fortran
makef90depend_ignore_args = \
	--ignore-module=ieee_arithmetic \
	--ignore-module=iso_c_binding \
	--ignore-module=omp_lib

if USE_METIS
EXTRA_PROGRAMS += \
	graph_partition
else
makef90depend_ignore_args += --ignore-module=ppm_graph_partition_serial
endif
endif

AM_FCFLAGS = $(FPP_INCOPT)$(top_srcdir)/include/f90 \
	$(FPP_INCOPT)../include/f77 \
	$(FC_MOD_FLAG)../include/f90 $(MPI_FC_INCLUDE)

AM_CFLAGS = -I$(top_srcdir)/src $(MPI_C_INCLUDE)

# more programs in this directory are Fortran programs, so this is default...
LDADD = ../src/libtestdata.la ../src/libscalesppm.la

show_ieee_emulation_SOURCES = show_ieee_emulation.f90

graph_build_SOURCES = graph_build.f90

graph_partition_SOURCES = graph_partition.f90

# ...and C programs have to be treated specially
check_prng_SOURCES = check_prng.c
check_prng_LDADD = ../src/libscalesppmcore.la

check_prng_f_SOURCES = check_prng_f.f90

# overrides to make regular compilations also use libtool
COMPILE = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(CC) -static $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)
FCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(FC) -static $(AM_FCFLAGS) $(FCFLAGS)

if with_fortran
./$(DEPDIR)/FC.deps: $(SOURCES) Makefile
	$(AM_V_GEN)unique=`for i in $(SOURCES); do \
	  if test -f "$$i"; then echo "$$i"; else echo "$(srcdir)/$$i"; fi; \
	  done | grep '\.f90$$' ` ; \
	  FC="$(FC)" FPP="$(FPP)" \
	  $(PERL) -I$(top_srcdir)/scripts \
	  $(top_srcdir)/scripts/makef90depends \
	  $(makef90depend_ignore_args) \
	  --strip-obj-dirname --parse-cache=./$(DEPDIR)/FC.depdb \
	  --fc-mod-flag=$(FC_MOD_FLAG) --fc-def-opt=$(FPP_DEFOPT) \
	  --fc-inc-flag=$(FPP_INCOPT) --fc-mod-suf=.$(FCMODEXT) \
	  --fc-mod-case=$(FCMODCASE) --prog-fpp='$(FPP)' \
	  --fpp-sizeof-int=$(C_SIZEOF_INT) \
	  --fpp-sizeof-long=$(C_SIZEOF_LONG) \
	  --fpp-sizeof-long-long=$(C_SIZEOF_LONG_LONG) \
	  --fpp-char-is-`test $(C_CHAR_IS_UNSIGNED) = 0 || echo un`signed \
	  -- \
	  $(AM_FCFLAGS) $(FCFLAGS) $(FPPFLAGS) -- $$unique >./$(DEPDIR)/FC.deps

@am__include@ @am__quote@./$(DEPDIR)/FC.deps@am__quote@ # am--include-marker
endif

clean-local:
	-rm -f *.$(FCMODEXT) i.*.L *.i ??????????.dbg
	-rm -rf *.dSYM

# When an out-of-source build takes place, the 'distclean' recipe of the
# top-level Makefile expects the directory containing this Makefile to be
# empty as soon as its 'distclean' recipe is executed. This is not the case
# when the Makefile is generated with Automake 1.16 or later because the
# $(DEPDIR) subdirectory is not deleted. Therefore, we delete it explicitly
# as the older versions of Automake do.
distclean-local:
	-rm -rf '$(DEPDIR)'
