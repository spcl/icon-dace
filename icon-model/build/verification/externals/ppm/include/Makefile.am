#
# Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
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
SUBDIRS=f90 f77

EXTRA_DIST=xpfpa.h xpfpa_func.h cfortran.h

coreincludedir = $(includedir)/core
ppmincludedir = $(includedir)/ppm

nodist_coreinclude_HEADERS = core/core.h \
	core/combinatorics.h \
	core/ppm_std_type_kinds_mp.h \
	core/ppm_math_extensions.h \
	core/ppm_combinatorics.h \
	core/ppm_extents.h \
	core/ppm_random.h \
	core/ppm_rectilinear.h \
	core/ppm_xfuncs.h

nodist_ppminclude_HEADERS = \
	ppm/ppm.h \
	ppm/ppm_set_partition_base.h \
	ppm/ppm_uniform_partition.h

if ENABLE_DISTRIBUTED_ARRAY_C
nodist_ppminclude_HEADERS += \
	ppm/dist_array.h
endif

if USE_MPI
nodist_coreinclude_HEADERS += \
	core/ppm_extents_mp.h
endif


distclean-local:
	-rm -f $(nodist_coreinclude_HEADERS) \
		$(nodist_ppminclude_HEADERS)
	@for d in $(SUBDIRS) core ppm; do \
		dd=$$d/$(DEPDIR) ; test ! -d $$dd \
			|| (rm -f $$dd/FC.dep* ; rmdir $$dd) ; \
		test "$(abs_builddir)" = "$(abs_srcdir)" \
			-o $$d = . -o ! -d $$d || rmdir $$d ; done

dist-hook:
	-rmdir $(distdir)/core $(distdir)/ppm

REWRITE_HEADER = $(AM_V_GEN)$(PERL) -I$(top_srcdir)/scripts \
	  $(top_srcdir)/scripts/header2installheader.pl \
	  --config-header=config.h \
	  --srcdir=$(top_srcdir)/src \
	  --c-sizeof-int=$(C_SIZEOF_INT) \
	  --c-sizeof-long=$(C_SIZEOF_LONG) \
	  --c-sizeof-long-long=$(C_SIZEOF_LONG_LONG) \
	  --c-char-is-`test $(C_CHAR_IS_UNSIGNED) = 0 || echo un`signed \
	  -D PPM_DSO_API_EXPORT -D PPM_DSO_INTERNAL \
	  --dstdir=.

core/core.h: $(top_srcdir)/src/core/core.h config.h
	$(REWRITE_HEADER) $<

core/combinatorics.h: $(top_srcdir)/src/core/combinatorics.h \
  config.h
	$(REWRITE_HEADER) $<

core/ppm_combinatorics.h: $(top_srcdir)/src/core/ppm_combinatorics.h \
  config.h
	$(REWRITE_HEADER) $<

core/ppm_std_type_kinds_mp.h: $(top_srcdir)/src/core/ppm_std_type_kinds_mp.h \
  config.h
	$(REWRITE_HEADER) $<

core/ppm_math_extensions.h: $(top_srcdir)/src/core/ppm_math_extensions.h \
  config.h
	$(REWRITE_HEADER) $<

core/ppm_extents.h: $(top_srcdir)/src/core/ppm_extents.h config.h
	$(REWRITE_HEADER) $<

core/ppm_random.h: $(top_srcdir)/src/core/ppm_random.h config.h
	$(REWRITE_HEADER) $<

core/ppm_rectilinear.h: $(top_srcdir)/src/core/ppm_rectilinear.h config.h
	$(REWRITE_HEADER) $<

core/ppm_xfuncs.h: $(top_srcdir)/src/core/ppm_xfuncs.h config.h
	$(REWRITE_HEADER) $<

core/ppm_extents_mp.h: $(top_srcdir)/src/core/ppm_extents_mp.h config.h
	$(REWRITE_HEADER) $<

ppm/ppm_set_partition_base.h: $(top_srcdir)/src/ppm/ppm_set_partition_base.h \
  config.h
	$(REWRITE_HEADER) $<

ppm/ppm_uniform_partition.h: $(top_srcdir)/src/ppm/ppm_uniform_partition.h \
  config.h
	$(REWRITE_HEADER) $<

ppm/dist_array.h: $(top_srcdir)/src/ppm/dist_array.h config.h
	$(REWRITE_HEADER) $<

ppm/ppm.h: $(top_srcdir)/src/ppm/ppm.h config.h
	$(REWRITE_HEADER) $<


