#
# Copyright  (C)  2022 Jörg Behrens <behrens@dkrz.de>
#                      Moritz Hanke <hanke@dkrz.de>
#                      Thomas Jahns <jahns@dkrz.de>
#
# Author: Jörg Behrens <behrens@dkrz.de>
#         Moritz Hanke <hanke@dkrz.de>
#         Thomas Jahns <jahns@dkrz.de>
#
# Maintainer: Jörg Behrens <behrens@dkrz.de>
#             Moritz Hanke <hanke@dkrz.de>
#             Thomas Jahns <jahns@dkrz.de>
# URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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

AM_CPPFLAGS = -I../include -I$(srcdir)/../include $(MPI_C_INCLUDE)
AM_FCFLAGS  = $(MPI_FC_MOD) $(FPP_INCOPT)../include
lib_LTLIBRARIES = libyaxt_c.la
if with_fortran
lib_LTLIBRARIES += libyaxt.la
else
AM_CPPFLAGS += -DWITHOUT_FORTRAN
endif

SUBDIRS = mpi-workarounds .
DIST_SUBDIRS = mpi-workarounds .

# headers from subdir xt will be installed (after substituting
# definitions from config.h). These must be listed in a sequence in
# which they can be combined in yaxt.h.
libyaxt_c_la_SOURCES =						\
	xt/xt_core.h						\
	xt/xt_config.h						\
	xt/sort_common.h					\
	xt/mergesort.h						\
	xt/quicksort.h						\
	xt/xt_sort.h						\
	xt_sort_internal.h					\
	xt/xt_idxlist.h						\
	xt/xt_idxempty.h					\
	xt/xt_idxlist_collection.h				\
	xt/xt_idxsection.h					\
	xt/xt_idxstripes.h					\
	xt/xt_idxvec.h 						\
	xt/xt_mpi.h						\
	xt/xt_redist.h						\
	xt/xt_redist_collection.h				\
	xt/xt_redist_collection_static.h			\
	xt/xt_redist_p2p.h					\
	xt/xt_redist_repeat.h					\
	xt/xt_stripe.h						\
	xt/xt_xmap.h						\
	xt/xt_xmap_all2all.h					\
	xt/xt_xmap_dist_dir.h					\
	xt/xt_xmap_intersection.h				\
	xt/xt_idxmod.h						\
	xt/xt_request.h						\
	xt/xt_request_msgs.h					\
	xt/xt_request_msgs_ebuf.h				\
	xt/xt_request_msgs_packed.h				\
	xt_request_msgs_ddt_packed.h				\
	xt_arithmetic_util.h					\
	xt_arithmetic_long.h					\
	xt_config_internal.h					\
	xt_config.c						\
	xt_cover.h						\
	xt_cover.c						\
	xt_idxempty_internal.h					\
	xt_idxempty.c						\
	xt_idxlist_internal.h					\
	xt_idxlist.c						\
	xt_idxlist_intersection.c				\
	xt_idxlist_collection.c					\
	xt_idxlist_collection_internal.h			\
	xt_idxlist_unpack.h					\
	xt_idxlist_unpack.c 					\
	xt_idxmod.c						\
	xt_idxsection_internal.h				\
	xt_idxsection.c						\
	xt_idxstripes_internal.h				\
	xt_idxstripes.c						\
	xt_idxvec.c 						\
	xt_idxvec_internal.h					\
	xt_mpi_internal.h					\
	xt_mpi.c						\
	xt_redist_internal.h					\
	xt_redist.c						\
	xt_redist_msg.c						\
	xt_redist_collection.c					\
	xt_redist_collection_static.c				\
	xt_redist_p2p.c						\
	xt_redist_repeat.c					\
	xt_redist_single_array_base.c				\
	xt/xt_redist_single_array_base.h			\
	xt_request.c						\
	xt_request_msgs.c					\
	xt_request_msgs_ebuf.c					\
	xt_request_msgs_packed.c				\
	xt_request_msgs_ddt_packed.c				\
	xt_request_internal.h					\
	xt_request_msgs_internal.h				\
	xt_request_msgs_ebuf_internal.h				\
	xt_stripe.c						\
	xt_stripe_util.h					\
	xt_xmap_internal.h					\
	xt_xmap.c						\
	xt_xmap_all2all.c					\
	xt_xmap_dist_dir.c					\
	xt_xmap_dist_dir_common.h				\
	xt_xmap_dist_dir_common.c				\
	xt/xt_xmap_dist_dir_intercomm.h				\
	xt_xmap_dist_dir_intercomm.c				\
	xt_xmap_intersection.c					\
	xt_xmap_intersection_common.h				\
	xt_xmap_intersection_ext.c				\
	xt_init_internal.h					\
	xt_init.c						\
	core/core.c						\
	core/core.h						\
	core/ppm_xfuncs.h					\
	core/symprefix.h					\
	core/ppm_visibility.h					\
	core/xmalloc.c						\
	core/xstdio.c						\
	ensure_array_size.c					\
	ensure_array_size.h					\
	xt_quicksort_base.h					\
	quicksort.c						\
	mergesort.c						\
	instr.h							\
	xt_sort.c						\
	xt_exchanger.c						\
	xt_exchanger_team.c					\
	xt_exchanger_vtable.c					\
	xt_exchanger.h						\
	xt_exchanger_simple_base.c				\
	xt_exchanger_simple_base.h				\
	xt_exchanger_irecv_send.c				\
	xt_exchanger_irecv_send.h				\
	xt_exchanger_irecv_isend.c				\
	xt_exchanger_irecv_isend.h				\
	xt_exchanger_irecv_isend_packed.c			\
	xt_exchanger_irecv_isend_packed.h			\
	xt_exchanger_irecv_isend_ddt_packed.c			\
	xt_exchanger_irecv_isend_ddt_packed.h			\
	xt_exchanger_neigh_alltoall.h				\
	xt_exchanger_mix_isend_irecv.c				\
	xt_exchanger_mix_isend_irecv.h				\
	xt_gpu.h						\
	xt_gpu.c						\
	xt_ddt.h						\
	xt_ddt_internal.h					\
	xt_ddt.c

if USE_NB_A2A
libyaxt_c_la_SOURCES +=						\
	xt_exchanger_neigh_alltoall.c
endif
if HAVE_CUDA
libyaxt_c_la_SOURCES += xt_cuda.h				\
	xt_cuda.c
endif
libyaxt_la_SOURCES =						\
	xt_idxlist_collection_f.f90				\
	xt_xmap_intersection_f.f90				\
	yaxt_f2c.c						\
	core/core_cf.c						\
	xt_slice_c_loc.inc					\
	xt_mpi_f.f90						\
	xt_core_f.f90						\
	xt_config_f.f90						\
	xt_sort_f.f90						\
	xt_idxlist_f.f90					\
	xt_idxvec_f.f90						\
	xt_idxstripes_f.f90					\
	xt_idxsection_f.f90					\
	xt_xmap_f.f90						\
	xt_redist_f.f90						\
	xt_redist_int_i2.f90					\
	xt_redist_int_i4.f90					\
	xt_redist_int_i8.f90					\
	xt_redist_real_sp.f90					\
	xt_redist_real_dp.f90					\
	xt_redist_logical.f90					\
	xt_request_f.f90					\
	yaxt.f90

libyaxt_la_LIBADD = libyaxt_c.la $(MPI_FC_LIB)
libyaxt_la_LDFLAGS = -version-info 11:1:10
libyaxt_c_la_LIBADD = $(MPI_C_LIB)
libyaxt_c_la_LDFLAGS = -version-info 10:4:9

if XT_ADD_MPI_WORKAROUND
libyaxt_c_la_LIBADD += mpi-workarounds/libmpi_workaround.la
endif

EXTRA_DIST = m4-templates/xt_redist_template_module.f90.m4	\
	m4-templates/f03-redist-gen.m4				\
	m4-templates/int-prologue-i2.m4				\
	m4-templates/int-prologue-i4.m4				\
	m4-templates/int-prologue-i8.m4				\
	m4-templates/real-prologue-dp.m4			\
	m4-templates/real-prologue-sp.m4			\
	xt_redist_p2p_ext.c					\
	xt_mpi_stripe_parse.c

f03_redist_module_template=`i=m4-templates/xt_redist_template_module.f90.m4 ; \
	if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; `
update_orig=true

xt_redist_real_sp.f90: \
	m4-templates/xt_redist_template_module.f90.m4 \
	m4-templates/f03-redist-gen.m4
	$(M4) -I$(srcdir)/m4-templates -I$(top_srcdir)/m4 \
		-Dfname=$@ -Dmodname=xt_redist_real_sp \
		-Dprologue='{real-prologue-sp.m4}' \
		-Dtypetag=sp -Dtypedecl='REAL(sp)' \
		$< >"$@.tmp"
	if $(update_orig) ; then \
		cmp "$@.tmp" "$(srcdir)/$@" >/dev/null 2>&1 \
		|| mv "$@.tmp" "$(srcdir)/$@" ; else \
		cp "$@.tmp" "$@" ; fi
	-rm -f "$@.tmp"

xt_redist_real_dp.f90: \
	m4-templates/xt_redist_template_module.f90.m4 \
	m4-templates/f03-redist-gen.m4
	$(M4) -I$(srcdir)/m4-templates -I$(top_srcdir)/m4 \
		-Dfname=$@ -Dmodname=xt_redist_real_dp \
		-Dprologue='{real-prologue-dp.m4}' \
		-Dtypetag=dp -Dtypedecl='REAL(dp)' \
		$< >"$@.tmp"
	if $(update_orig) ; then \
		cmp "$@.tmp" "$(srcdir)/$@" >/dev/null 2>&1 \
		|| mv "$@.tmp" "$(srcdir)/$@" ; else \
		cp "$@.tmp" "$@" ; fi
	-rm -f "$@.tmp"

xt_redist_logical.f90: \
	m4-templates/xt_redist_template_module.f90.m4 \
	m4-templates/f03-redist-gen.m4
	$(M4) -I$(srcdir)/m4-templates -I$(top_srcdir)/m4 \
		-Dfname=$@ -Dmodname=xt_redist_logical \
		-Dtypetag=l -Dtypedecl='LOGICAL' \
		$< >"$@.tmp"
	if $(update_orig) ; then \
		cmp "$@.tmp" "$(srcdir)/$@" >/dev/null 2>&1 \
		|| mv "$@.tmp" "$(srcdir)/$@" ; else \
		cp "$@.tmp" "$@" ; fi
	-rm -f "$@.tmp"

xt_redist_int_i4.f90: \
	m4-templates/xt_redist_template_module.f90.m4 \
	m4-templates/f03-redist-gen.m4
	$(M4) -I$(srcdir)/m4-templates -I$(top_srcdir)/m4 \
		-Dfname=$@ -Dmodname=xt_redist_int_i4 \
		-Dprologue='{int-prologue-i4.m4}' \
		-Dtypetag=i4 -Dtypedecl='INTEGER(i4)' \
		$< >"$@.tmp"
	if $(update_orig) ; then \
		cmp "$@.tmp" "$(srcdir)/$@" >/dev/null 2>&1 \
		|| mv "$@.tmp" "$(srcdir)/$@" ; else \
		cp "$@.tmp" "$@" ; fi
	-rm -f "$@.tmp"

xt_redist_int_i2.f90: \
	m4-templates/xt_redist_template_module.f90.m4 \
	m4-templates/f03-redist-gen.m4
	$(M4) -I$(srcdir)/m4-templates -I$(top_srcdir)/m4 \
		-Dfname=$@ -Dmodname=xt_redist_int_i2 \
		-Dprologue='{int-prologue-i2.m4}' \
		-Dtypetag=i2 -Dtypedecl='INTEGER(i2)' \
		$< >"$@.tmp"
	if $(update_orig) ; then \
		cmp "$@.tmp" "$(srcdir)/$@" >/dev/null 2>&1 \
		|| mv "$@.tmp" "$(srcdir)/$@" ; else \
		cp "$@.tmp" "$@" ; fi
	-rm -f "$@.tmp"

xt_redist_int_i8.f90: \
	m4-templates/xt_redist_template_module.f90.m4 \
	m4-templates/f03-redist-gen.m4
	$(M4) -I$(srcdir)/m4-templates -I$(top_srcdir)/m4 \
		-Dfname=$@ -Dmodname=xt_redist_int_i8 \
		-Dprologue='{int-prologue-i8.m4}' \
		-Dtypetag=i8 -Dtypedecl='INTEGER(i8)' \
		$< >"$@.tmp"
	if $(update_orig) ; then \
		cmp "$@.tmp" "$(srcdir)/$@" >/dev/null 2>&1 \
		|| mv "$@.tmp" "$(srcdir)/$@" ; else \
		cp "$@.tmp" "$@" ; fi
	-rm -f "$@.tmp"


xtinclude_Makefile.inc: Makefile
	$(AM_V_GEN)headers=`echo $(libyaxt_c_la_SOURCES) | sed 's/[ 	][ 	]*/ /g' | tr ' ' '\n' | grep 'xt/.*\.h$$'` ; \
	  echo "xtinclude_HEADERS = \\" >$@ ; \
	  for header in $$headers ; do echo "	$$header \\" ; done \
	    | sed '$$s/\\$$//' >>$@

all-local: xtinclude_Makefile.inc

clean-local:
	-rm -f *.$(FCMODEXT) i.*.L *.i *.dbg


# When an out-of-source build takes place, the recipe of the rule expects that
# the object files are already deleted from the 'core' subdirectory, which is
# not always the case for a parallel run, i.e. when 'make distclean -j2'.
# Therefore, we explicitly specify a dependency on the 'mostlyclean' rule.
# Moreover, the 'distclean' recipe of the top-level Makefile expects the
# directory containing this Makefile to be empty as soon as its 'distclean'
# recipe is executed. This is not the case when the Makefile is generated with
# Automake 1.16 or later because the $(DEPDIR) subdirectories are not deleted.
# Therefore, we delete them explicitly as the older versions of Automake do.
# We also need to clean-up the 'pkgconfig' and 'core' subdirectories before
# trying to delete them in the out-of-source build case. First, we make sure that
# the directory stamps are deleted by introducing a dependency on
# 'distclean-generic'. Second, we delete the pkgconfig files, which, otherwise,
# are deleted later by the top-level Makefile.
PKG_CFG_FILES = pkgconfig/yaxt_c.pc
if with_fortran
PKG_CFG_FILES += pkgconfig/yaxt.pc
endif
distclean-local: mostlyclean-compile distclean-generic
	-rm -f xtinclude_Makefile.inc
	-rm -rf '$(DEPDIR)' 'core/$(DEPDIR)' 'mpi-profile/$(DEPDIR)'
	-if test "$(abs_builddir)" != "$(abs_srcdir)" ; then \
	  rm -f pkgconfig/yaxt.pc pkgconfig/yaxt_c.pc; rmdir core pkgconfig mpi-workarounds; fi

install-exec-local:
	$(mkinstalldirs) "$(DESTDIR)$(libdir)/pkgconfig"
	for pc in $(PKG_CFG_FILES) ; do \
		$(install_sh_DATA) "$$pc" "$(DESTDIR)$(libdir)/$$pc" ; done

uninstall-local:
	rm -f "$(DESTDIR)$(libdir)/pkgconfig/yaxt.pc" \
		"$(DESTDIR)$(libdir)/pkgconfig/yaxt_c.pc"
	-rmdir "$(DESTDIR)$(libdir)/pkgconfig"

if with_fortran
./$(DEPDIR)/FC.deps: $(libyaxt_la_SOURCES) Makefile
	@unique=`for i in $(SOURCES); do \
	  if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | grep '\.f90$$' ` ; \
	  FC="$(FC)" FPP="$(FPP)" \
	  $(PERL) -I$(top_srcdir)/scripts \
	  $(top_srcdir)/scripts/makef90depends \
	  --ignore-module=iso_c_binding --ignore-module=mpi \
	  --fc-mod-flag=$(FC_MOD_FLAG) --fc-def-opt=$(FPP_DEFOPT) \
	  --fc-inc-flag=$(FPP_INCOPT) --fc-mod-suf='.$(FCMODEXT)' \
	  --fc-mod-case=$(FCMODCASE) --prog-fpp='$(FPP)' \
	  --fpp-sizeof-int=$(C_SIZEOF_INT) \
	  --fpp-sizeof-long=$(C_SIZEOF_LONG) \
	  --fpp-sizeof-long-long=$(C_SIZEOF_LONG_LONG) \
	  --fpp-char-is-`test $(C_CHAR_IS_UNSIGNED) = 0 || echo un`signed \
	  --strip-obj-dirname --parse-cache=./$(DEPDIR)/FC.depdb \
	  --obj-suf='.lo' -- \
	  $(AM_FCFLAGS) $(FCFLAGS) $(FPPFLAGS) -- $$unique \
	  >./$(DEPDIR)/FC.deps.temp && \
	cmp ./$(DEPDIR)/FC.deps.temp ./$(DEPDIR)/FC.deps >/dev/null 2>&1 \
	  || ln -f ./$(DEPDIR)/FC.deps.temp ./$(DEPDIR)/FC.deps \
	  && touch ./$(DEPDIR)/FC.deps

@am__include@ @am__quote@./$(DEPDIR)/FC.deps@am__quote@ # am--include-marker
endif

#
# Local Variables:
# coding: utf-8
# show-trailing-whitespace: t
# End:
#
