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

check_PROGRAMS = 						\
	test_mpi_generate_datatype				\
	test_mpi_smartdedup					\
	test_exchanger_parallel					\
	test_arithmetic_long					\
	test_idxempty						\
	test_idxvec						\
	test_idxlist_collection					\
	test_idxmod						\
	test_idxsection						\
	test_idxstripes						\
	test_redist_single_array_base				\
	test_redist_single_array_base_parallel			\
	test_redist_collection					\
	test_redist_collection_parallel				\
	test_redist_collection_static				\
	test_redist_collection_static_parallel			\
	test_redist_p2p						\
	test_redist_p2p_parallel				\
	test_redist_repeat					\
	test_redist_repeat_parallel				\
	test_request_parallel					\
	test_xmap_all2all					\
	test_xmap_all2all_fail					\
	test_xmap_all2all_parallel				\
	test_xmap_all2all_intercomm_parallel			\
	test_xmap_dist_dir					\
	test_xmap_dist_dir_parallel				\
	test_xmap_dist_dir_intercomm_parallel			\
	test_xmap_intersection_parallel				\
	test_initialized_finalized				\
	test_sort						\
	test_uid						\
	test_ddt

check_LTLIBRARIES = libtestutil.la

if with_fortran
check_PROGRAMS += 						\
	test_misc_f						\
	test_idxempty_f						\
	test_idxvec_f						\
	test_idxlist_collection_f				\
	test_idxmod_f						\
	test_idxsection_f					\
	test_idxstripes_f					\
	test_redist_single_array_base_f				\
	test_redist_single_array_base_parallel_f		\
	test_redist_collection_f				\
	test_redist_collection_parallel_f			\
	test_redist_collection_static_f				\
	test_redist_collection_static_parallel_f		\
	test_redist_p2p_f					\
	test_redist_p2p_parallel_f				\
	test_redist_repeat_f					\
	test_redist_repeat_parallel_f				\
	test_xmap_all2all_f					\
	test_xmap_all2all_fail_f				\
	test_xmap_all2all_parallel_f				\
	test_xmap_all2all_intercomm_parallel_f			\
	test_xmap_dist_dir_f					\
	test_xmap_dist_dir_parallel_f				\
	test_xmap_dist_dir_intercomm_parallel_f			\
	test_xmap_intersection_parallel_f			\
	test_initialized_finalized_f				\
	test_ut							\
	test_perf						\
	test_perf_stripes					\
	test_sort_f						\
	test_yaxt

check_LTLIBRARIES += libtestutil_f.la
endif


if ! with_on_demand_check_programs
all: $(check_LTLIBRARIES) $(check_PROGRAMS)
endif


libtestutil_la_SOURCES = tests.c \
	cksum.c cksum.h \
	ctest_common.c ctest_common.h \
	test_idxlist_utils.c \
	test_xmap_common.c test_xmap_common.h \
	test_xmap_common_parallel.c \
	test_xmap_common_intercomm_parallel.c \
	test_redist_single_array_base_common.c \
	test_redist_common.h test_redist_common.c

libtestutil_f_la_SOURCES = ftest_common.f90 \
	ftest_common_c.c ftest_common.h \
	test_idxlist_utils_f.f90 \
	test_xmap_common_f.f90 \
	test_xmap_common_parallel_f.f90 \
	test_xmap_common_intercomm_parallel_f.f90 \
	test_redist_common_f.f90

libtestutil_f_la_LDFLAGS = -static
libtestutil_la_LDFLAGS = -static
# overrides to make regular compilations also use libtool
# static because the libtool libraries in this directory are meant to
# be statically linked and there is no use for the PIC object
LTCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(CC) -static $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)
COMPILE = $(LTCOMPILE)
LTFCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(FC) -static $(AM_FCFLAGS) $(FCFLAGS)
FCCOMPILE = $(LTFCCOMPILE)
LINK = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=link $(CCLD) $(AM_CFLAGS) $(CFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@
FCLINK = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=link $(FCLD) $(AM_FCFLAGS) $(FCFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@

test_arithmetic_long_SOURCES = test_arithmetic_long.c
test_arithmetic_long_LDADD = libtestutil.la
test_exchanger_parallel_SOURCES = test_exchanger_parallel.c tests.h
test_exchanger_parallel_LDADD =  libtestutil.la \
	../src/xt_gpu.lo \
	../src/xt_ddt.lo \
	../src/xt_config.lo \
	../src/xt_exchanger_irecv_isend.lo \
	../src/xt_exchanger_irecv_send.lo \
	../src/xt_exchanger_irecv_isend_packed.lo \
	../src/xt_exchanger_irecv_isend_ddt_packed.lo \
	../src/xt_exchanger_mix_isend_irecv.lo \
	../src/xt_exchanger_simple_base.lo
if USE_NB_A2A
test_exchanger_parallel_LDADD += ../src/xt_exchanger_neigh_alltoall.lo
endif
if HAVE_CUDA
test_exchanger_parallel_LDADD += ../src/xt_cuda.lo
endif
test_exchanger_parallel_LDADD += ../src/xt_exchanger.lo \
	../src/xt_redist_msg.lo \
	../src/xt_request_msgs.lo \
	../src/xt_request_msgs_ebuf.lo \
	../src/xt_request_msgs_packed.lo \
	../src/xt_request_msgs_ddt_packed.lo \
	../src/xt_request.lo \
	../src/xt_mpi.lo \
	../src/core/xmalloc.lo \
	../src/core/core.lo $(MPI_C_LIB)

@ifGNUmake@../src/%.lo: ../src/%.c
@ifGNUmake@	$(MAKE) -C ../src $(notdir $@)

test_idxempty_SOURCES = test_idxempty.c tests.h test_idxlist_utils.h
test_idxempty_f_SOURCES = test_idxempty_f.f90
test_idxempty_f_LDADD = $(XT_FC_LDADD)
test_idxvec_SOURCES = test_idxvec.c tests.h test_idxlist_utils.h
test_idxvec_f_SOURCES = test_idxvec_f.f90
test_idxvec_f_LDADD = $(XT_FC_LDADD)
test_idxlist_collection_SOURCES = test_idxlist_collection.c tests.h \
	test_idxlist_utils.h
test_idxlist_collection_f_SOURCES = test_idxlist_collection_f.f90
test_idxlist_collection_f_LDADD = $(XT_FC_LDADD)
test_idxmod_SOURCES = test_idxmod.c tests.h
test_idxmod_f_SOURCES = test_idxmod_f.f90
test_idxmod_f_LDADD = $(XT_FC_LDADD)
test_idxsection_SOURCES = test_idxsection.c tests.h test_idxlist_utils.h
test_idxsection_f_SOURCES = test_idxsection_f.f90
test_idxsection_f_LDADD = $(XT_FC_LDADD)
test_idxstripes_SOURCES = test_idxstripes.c tests.h test_idxlist_utils.h
test_idxstripes_f_SOURCES = test_idxstripes_f.f90
test_idxstripes_f_LDADD = $(XT_FC_LDADD)
test_redist_single_array_base_SOURCES = test_redist_single_array_base.c tests.h
test_redist_single_array_base_f_SOURCES = test_redist_single_array_base_f.f90
test_redist_single_array_base_f_LDADD = $(XT_FC_LDADD)
test_redist_single_array_base_parallel_SOURCES = test_redist_single_array_base_parallel.c tests.h
test_redist_single_array_base_parallel_f_SOURCES = test_redist_single_array_base_parallel_f.f90
test_redist_single_array_base_parallel_f_LDADD = $(XT_FC_LDADD)
test_redist_collection_SOURCES = test_redist_collection.c tests.h
test_redist_collection_f_SOURCES = test_redist_collection_f.f90 \
	test_redist_collection_displace_f.f90 xt_slice_c_loc.inc
test_redist_collection_f_LDADD = $(XT_FC_LDADD)
test_redist_collection_parallel_SOURCES = test_redist_collection_parallel.c \
	tests.h
test_redist_collection_parallel_f_SOURCES = \
	test_redist_collection_parallel_f.f90 xt_slice_c_loc.inc
test_redist_collection_parallel_f_LDADD = $(XT_FC_LDADD)
test_redist_collection_static_SOURCES = test_redist_collection_static.c tests.h
test_redist_collection_static_f_SOURCES = test_redist_collection_static_f.f90
test_redist_collection_static_f_LDADD = $(XT_FC_LDADD)
test_redist_collection_static_parallel_SOURCES = \
	test_redist_collection_static_parallel.c tests.h
test_redist_collection_static_parallel_f_SOURCES = \
	test_redist_collection_static_parallel_f.f90
test_redist_collection_static_parallel_f_LDADD = $(XT_FC_LDADD)
test_redist_p2p_SOURCES = test_redist_p2p.c tests.h
test_redist_p2p_f_SOURCES = test_redist_p2p_f.f90
test_redist_p2p_f_LDADD = $(XT_FC_LDADD)
test_redist_p2p_parallel_SOURCES = test_redist_p2p_parallel.c tests.h
test_redist_p2p_parallel_f_SOURCES = test_redist_p2p_parallel_f.f90
test_redist_p2p_parallel_f_LDADD = $(XT_FC_LDADD)
test_redist_repeat_SOURCES = test_redist_repeat.c tests.h
test_redist_repeat_f_SOURCES = test_redist_repeat_f.f90
test_redist_repeat_f_LDADD = $(XT_FC_LDADD)
test_redist_repeat_parallel_SOURCES = test_redist_repeat_parallel.c tests.h
test_redist_repeat_parallel_f_SOURCES = test_redist_repeat_parallel_f.f90
test_redist_repeat_parallel_f_LDADD = $(XT_FC_LDADD)
test_request_parallel_SOURCES = test_request_parallel.c tests.h
test_xmap_all2all_SOURCES = test_xmap_all2all.c test_xmap_common.h
test_xmap_all2all_f_SOURCES = test_xmap_all2all_f.f90
test_xmap_all2all_f_LDADD = $(XT_FC_LDADD)
test_xmap_all2all_fail_SOURCES = test_xmap_all2all_fail.c tests.h
test_xmap_all2all_fail_f_SOURCES = test_xmap_all2all_fail_f.f90
test_xmap_all2all_fail_f_LDADD = $(XT_FC_LDADD)
test_xmap_all2all_parallel_SOURCES = test_xmap_all2all_parallel.c tests.h
test_xmap_all2all_parallel_f_SOURCES = test_xmap_all2all_parallel_f.f90
test_xmap_all2all_parallel_f_LDADD = $(XT_FC_LDADD)
test_xmap_all2all_intercomm_parallel_SOURCES = \
	test_xmap_all2all_intercomm_parallel.c tests.h
test_xmap_all2all_intercomm_parallel_f_SOURCES = \
	test_xmap_all2all_intercomm_parallel_f.f90
test_xmap_all2all_intercomm_parallel_f_LDADD = $(XT_FC_LDADD)
test_xmap_dist_dir_SOURCES = test_xmap_dist_dir.c test_xmap_common.h
test_xmap_dist_dir_f_SOURCES = test_xmap_dist_dir_f.f90
test_xmap_dist_dir_f_LDADD = $(XT_FC_LDADD)
test_xmap_dist_dir_parallel_SOURCES = test_xmap_dist_dir_parallel.c tests.h
test_xmap_dist_dir_parallel_f_SOURCES = test_xmap_dist_dir_parallel_f.f90
test_xmap_dist_dir_parallel_f_LDADD = $(XT_FC_LDADD)
test_xmap_dist_dir_intercomm_parallel_SOURCES = \
	test_xmap_dist_dir_intercomm_parallel.c test_xmap_common.h tests.h
test_xmap_dist_dir_intercomm_parallel_f_SOURCES = \
	test_xmap_dist_dir_intercomm_parallel_f.f90
test_xmap_dist_dir_intercomm_parallel_f_LDADD = $(XT_FC_LDADD)
test_xmap_intersection_parallel_SOURCES = test_xmap_intersection_parallel.c tests.h
test_xmap_intersection_parallel_f_SOURCES = test_xmap_intersection_parallel_f.f90
test_xmap_intersection_parallel_f_LDADD = $(XT_FC_LDADD)
test_initialized_finalized_SOURCES = test_initialized_finalized.c tests.h
test_initialized_finalized_f_SOURCES = test_initialized_finalized_f.f90
test_initialized_finalized_f_LDADD = $(XT_FC_LDADD)
test_ut_SOURCES = test_ut.f90
test_ut_LDADD = $(XT_FC_LDADD)
test_perf_SOURCES = test_perf.f90
test_perf_LDADD = $(XT_FC_LDADD)
test_perf_stripes_SOURCES = test_perf_stripes.f90
test_perf_stripes_LDADD = $(XT_FC_LDADD)
test_sort_SOURCES = test_sort.c tests.h
test_sort_f_SOURCES = test_sort_f.f90
test_sort_f_LDADD = $(XT_FC_LDADD)
test_uid_SOURCES = test_uid.c
test_yaxt_SOURCES = test_yaxt.f90
test_yaxt_LDADD = $(XT_FC_LDADD)
test_misc_f_SOURCES = test_misc_f.f90
test_misc_f_LDADD = $(XT_FC_LDADD)
test_mpi_generate_datatype_SOURCES = test_mpi_generate_datatype.c tests.h
test_mpi_smartdedup_SOURCES = test_mpi_smartdedup.c tests.h
test_mpi_smartdedup_LDADD =  libtestutil.la ../src/xt_mpi.lo \
	../src/core/xmalloc.lo ../src/core/core.lo $(MPI_C_LIB)
test_ddt_SOURCES = test_ddt.c tests.h
test_ddt_LDADD = libtestutil.la ../src/xt_ddt.lo ../src/xt_gpu.lo
if HAVE_CUDA
test_ddt_LDADD += ../src/xt_cuda.lo
endif
test_ddt_LDADD += ../src/xt_mpi.lo \
	../src/core/xmalloc.lo ../src/core/core.lo $(MPI_C_LIB)
AM_CPPFLAGS = -I../inst_headers -I$(top_srcdir)/include $(MPI_C_INCLUDE)
AM_FCFLAGS = $(FC_MOD_FLAG)../inst_headers/f90 \
	$(FPP_INCOPT)../include $(MPI_FC_MOD)
AM_LDFLAGS = -L../src
LDADD = libtestutil.la ../src/libyaxt_c.la $(MPI_C_LIB)
XT_FC_LDADD = libtestutil_f.la libtestutil.la ../src/libyaxt.la \
	../src/libyaxt_c.la $(MPI_FC_LIB)

TESTS = 							\
	test_misc_run						\
	test_mpi_generate_datatype_run				\
	test_mpi_smartdedup_run					\
	test_exchanger_parallel_run				\
	test_arithmetic_long					\
	test_idxempty_run					\
	test_idxvec_run						\
	test_idxlist_collection_run				\
	test_idxmod_run						\
	test_idxsection_run					\
	test_idxstripes_run					\
	test_redist_single_array_base_run			\
	test_redist_single_array_base_parallel_run		\
	test_redist_p2p_run					\
	test_redist_collection_run				\
	test_redist_collection_parallel_run			\
	test_redist_collection_static_run			\
	test_redist_collection_static_parallel_run		\
	test_redist_p2p_parallel_run				\
	test_redist_repeat_run					\
	test_redist_repeat_parallel_run				\
	test_request_parallel_run				\
	test_xmap_all2all_run					\
	test_xmap_all2all_fail_run				\
	test_xmap_all2all_parallel_run				\
	test_xmap_all2all_intercomm_parallel_run		\
	test_xmap_dist_dir_run					\
	test_xmap_dist_dir_parallel_run				\
	test_xmap_dist_dir_intercomm_parallel_run		\
	test_xmap_intersection_parallel_run			\
	test_initialized_finalized_run				\
	test_exported_symbols					\
	test_ut_run						\
	test_perf_run						\
	test_perf_stripes_run					\
	test_sort_run						\
	test_uid_run						\
	test_yaxt_run						\
	test_xmap_all2all_fail_run				\
	test_ddt_run

AUTOMAKE_OPTIONS = color-tests

EXTRA_DIST = core

if with_fortran
./$(DEPDIR)/FC.deps: $(SOURCES) Makefile
	@unique=`for i in $(SOURCES); do \
	  if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | grep '\.f90$$' ` ; \
	  libnamerewriterules=`for i in $(libtestutil_f_la_OBJECTS); do \
	  echo "$$i" | sed -e 'h;s/^/--obj-transform-name=/' \
	  -e 's/\.lo$$/\.f90$$/;s/$$/=/;G;s/\n//' ; done` ; \
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
	  $$libnamerewriterules \
	  --obj-transform-name='\.f90$$=.$(OBJEXT)' -- \
	  $(AM_FCFLAGS) $(FCFLAGS) $(FPPFLAGS) -- $$unique \
	  >./$(DEPDIR)/FC.deps.temp && \
	cmp ./$(DEPDIR)/FC.deps.temp ./$(DEPDIR)/FC.deps >/dev/null 2>&1 \
	  || ln -f ./$(DEPDIR)/FC.deps.temp ./$(DEPDIR)/FC.deps

@am__include@ @am__quote@./$(DEPDIR)/FC.deps@am__quote@ # am--include-marker
endif

clean-local:
	-rm -f *.$(FCMODEXT) *.i *.dbg
	-rm -rf *.dSYM

# When an out-of-source build takes place, the 'distclean' recipe of the
# top-level Makefile expects the directory containing this Makefile to be
# empty as soon as its 'distclean' recipe is executed. This is not the case
# when the Makefile is generated with Automake 1.16 or later because the
# $(DEPDIR) subdirectory is not deleted. Therefore, we delete it explicitly
# as the older versions of Automake do.
distclean-local:
	-rm -rf '$(DEPDIR)'

#
# Local Variables:
# coding: utf-8
# show-trailing-whitespace: t
# End:
#
