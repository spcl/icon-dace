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
AM_FCFLAGS = $(MPI_FC_INCLUDE) $(FPP_INCOPT)$(top_srcdir)/include/f90 \
	$(FPP_INCOPT)$(top_srcdir)/include/f77 \
	$(FPP_INCOPT)../include/f77 \
	$(FPP_INCOPT)$(srcdir)/core \
	`if expr x'$<' : '.*\.f90' >/dev/null ; then echo '$<' | sed -e 's&^$(VPATH)/&&' -e 's&/[^/]*$$&/&;t prep' -e 's&.*&.&' -e ':prep' -e 's&^&$(FC_MODOUT)&'; fi` \
	$(FC_MOD_FLAG)core $(FC_MOD_FLAG)crypto $(FC_MOD_FLAG)ppm \
	$(FC_MOD_FLAG)solver

AM_CPPFLAGS = $(MPI_C_INCLUDE) -I$(top_srcdir)/include -I../include \
	$(PARMETIS_C_INCLUDE) $(CRYPTO_C_INCLUDE)

lib_LTLIBRARIES=libscalesppmcore.la
if with_fortran
lib_LTLIBRARIES += libscalesppm.la
else
AM_CPPFLAGS += -DWITHOUT_FORTRAN
endif

libscalesppmcore_la_LDFLAGS = -version-info 2:3:2
libscalesppm_la_LDFLAGS = -version-info 3:2:3

libscalesppm_la_SOURCES=ppm/scales_ppm.f90 \
	ppm/ppm_uniform_partition.f90 \
	ppm/ppm_uniform_partition_c.c \
	ppm/ppm_uniform_partition.h \
	ppm/ppm_general_block_decomposition.f90 \
	ppm/ppm_m1d.f90 \
	ppm/ppm_set_partition.f90 \
	ppm/ppm_set_repartition.f90 \
	ppm/ppm_set_partition_base.f90 \
	ppm/ppm_set_partition_base.h \
	ppm/ppm_set_partition_base_c.c \
	ppm/ppm_set_partition_base_cf.c \
	solver/solver_public.f90 \
	solver/solver_internal.f90 \
	solver/solver_config.f90 \
	solver/linear_algebra.f90 \
	solver/preconditioners.f90 \
	solver/spectral_methods.f90 \
	solver/solvers.f90 \
	solver/solver_all.f90 \
	\
	crypto/ppm_checksum.f90 crypto/ppm_checksum_c.c \
	crypto/ppm_checksum.h crypto/ppm_checksum_fwrap.c \
	\
	core/core.c \
	core/core.h \
	core/core_cf.c \
	core/symprefix.h \
	core/ppm_visibility.h \
	core/ppm_base.f90 \
	core/ppm_statistics.f90 \
	core/combinatorics.h \
	core/ppm_combinatorics.f90 \
	core/ppm_combinatorics.h \
	core/ppm_combinatorics_c.c \
	core/ppm_combinatorics_fwrap.c \
	core/ppm_math_extensions.f90 \
	core/ppm_real_sp_dp_edit_descriptor.f90 \
	core/ppm_math_extensions_c.c \
	core/ppm_math_extensions_cf.c \
	core/ppm_math_extensions_ddp_c.c \
	core/ppm_math_extensions_ddp_cf.c \
	core/ppm_math_extensions.h \
	core/ppm_math_extensions_internal.f90 \
	core/ppm_fpu_underflow.h \
	core/ppm_compare.f90 \
	core/ppm_external.f90 \
	core/yarandom.c \
	core/yafrand.c \
	core/yarandom.h \
	core/ppm_irand_internal.f90 \
	core/ppm_random.f90 \
	core/ppm_random.h \
	core/ppm_random_cf.c \
	core/ppm_random_c.c \
	core/ppm_omp_assign.h \
	core/ppm_random_mt_cf.c \
	core/ppm_random_mt_c.c \
	core/ppm_f90_io_lun.f90 \
	core/ppm_strio.f90 \
	core/ppm_strio_c.c \
	core/ppm_posix.f90 \
	core/ppm_posix_c.c \
	core/ppm_extents.f90 \
	core/ppm_extents_c.c \
	core/ppm_extents_f2c.c \
	core/ppm_extents.h \
	core/ppm_rectilinear.f90 \
	core/ppm_rectilinear_c.c \
	core/ppm_rectilinear.h \
	core/ppm_strided_extents.f90 \
	core/ppm_strided_extents.h \
	core/ppm_strided_extents_c.c \
	core/ppm_sparse_mask_index.f90 \
	core/ppm_compact_mask_index.f90 \
	core/ppm_graph_csr.f90 \
	core/ppm_graph_alist.f90 \
	core/heap.h \
	core/heap.c \
	core/heap_fwrap.c \
	core/bsearch.h \
	core/bsearch.c \
	core/bsearch_fwrap.c \
	core/insertion_sort.h \
	core/insertion_sort.c \
	core/insertion_sort_fwrap.c \
	core/swapmacros.h \
	core/ppm_xfuncs.h \
	core/xmalloc.c \
	core/xstdio.c \
	core/fptr_api.h \
	core/array.h \
	core/array_api.h \
	core/array.c \
	core/bittab.h \
	core/bittab_api.h \
	core/bittab.c \
	core/minmax.h \
	core/msort.h \
	core/msort_api.h \
	core/msort.c \
	core/qsort_r.c \
	core/qsort_fwrap.c \
	core/qsort_r.h \
	core/qsort_r_api.h \
	core/ppm_std_type_kinds.f90

EXTRA_libscalesppm_la_SOURCES = ppm/ppm_graph_partition_mpi.f90 \
	ppm/ppm_graph_partition_serial.f90 \
	ppm/ppm_distributed.f90 \
	ppm/dist_array_f2c.c \
	solver/linear_algebra_multi.f90	\
	solver/preconditioners_multi.f90 \
	solver/solver_all_multi.f90 \
	solver/solver_config_multi.f90 \
	solver/solver_internal_multi.f90 \
	solver/solvers_multi.f90 \
	solver/spectral_methods_multi.f90 \
	core/ieee_arithmetic.f90 \
	core/ppm_extents_mp.f90 \
	core/ppm_ncdf_dump.f90 core/ppm_ncdf_dump_wrappers.f90 \
	core/errhandler.c core/errhandler.h \
	core/errhandler_cf.c \
	crypto/md5.c crypto/md5.h \
	core/ppm_std_type_kinds_mp_c.c \
	core/ppm_std_type_kinds_mp.h \
	core/ppm_std_type_kinds_mp.f90 \
	core/ppm_real_sp_dp.inc \
	core/ppmcommon.inc

libscalesppmcore_la_SOURCES = \
	ppm/ppm_set_partition_base.h \
	ppm/ppm_uniform_partition_c.c \
	ppm/ppm_uniform_partition.h \
	ppm/ppm.h \
	ppm/ppm_c.c \
	crypto/ppm_checksum_c.c \
	crypto/ppm_checksum.h \
	core/core.c \
	core/core.h \
	core/symprefix.h \
	core/combinatorics.h \
	core/ppm_combinatorics.h \
	core/ppm_combinatorics_c.c \
	core/ppm_extents_c.c \
	core/ppm_extents.h \
	core/yarandom.c \
	core/yafrand.c \
	core/yarandom.h \
	core/ppm_random.h \
	core/ppm_omp_assign.h \
	core/ppm_random_c.c \
	core/ppm_random_mt_c.c \
	core/heap.h \
	core/heap.c \
	core/bsearch.h \
	core/bsearch.c \
	core/insertion_sort.h \
	core/insertion_sort.c \
	core/swapmacros.h \
	core/ppm_math_extensions_c.c \
	core/ppm_math_extensions_ddp_c.c \
	core/ppm_math_extensions.h \
	core/ppm_fpu_underflow.h \
	core/ppm_rectilinear_c.c \
	core/ppm_rectilinear.h \
	core/ppm_strided_extents.h \
	core/ppm_strided_extents_c.c \
	core/ppm_xfuncs.h \
	core/xmalloc.c \
	core/xstdio.c \
	core/fptr_api.h \
	core/array.h \
	core/array_api.h \
	core/array.c \
	core/bittab.h \
	core/bittab_api.h \
	core/bittab.c \
	core/minmax.h \
	core/msort.h \
	core/msort_api.h \
	core/msort.c \
	core/qsort_r.c \
	core/qsort_r.h \
	core/qsort_r_api.h

EXTRA_libscalesppmcore_la_SOURCES = \
	core/errhandler.c \
	core/errhandler.h \
	core/ppm_std_type_kinds_mp_c.c \
	core/ppm_std_type_kinds_mp.h \
	core/ppm_math_extensions_ddp_mp_c.c \
	core/ppm_extents_mp_c.c \
	core/ppm_extents_mp.h \
	ppm/dist_array_c.c \
	ppm/dist_array.h \
	crypto/md5.c \
	crypto/md5.h

AM_FPPFLAGS =

GENERATED_F90 = ppm/dist_array.f90 ppm/ppm_ptr_bnds_remap.f90
CLEANFILES =
nodist_libscalesppm_la_SOURCES =

if MAINTAINER_MODE
# extra dependencies on m4 macro files
@am__include@ @am__quote@$(srcdir)/m4.deps@am__quote@ # am--include-marker

$(srcdir)/m4.deps: \
	$(sort $(wildcard $(srcdir)/*/*.f90.m4))
	$(AM_V_GEN){ srcdirpat=`printf '%s' '$(srcdir)/' \
		| $(PERL) -e 'print quotemeta(<>);'` ; \
	  for f in $^ ; do \
	      sed -n -e '/include(/{' \
	        -e 's@include({\([^}]*\)})dnl@'"$$f"': $$(top_srcdir)/m4/\1@' \
		-e 's@'"$$srcdirpat"'\(.*\)\.m4: @\1: @;p;}' \
		"$$f" ; done >$@.new ; \
	  cmp $@ $@.new >/dev/null 2>&1 && rm $@.new || mv $@.new $@ ; }

%.f90: %.f90.m4
	$(AM_V_GEN)$(MKDIR_P) ppm \
		&& (echo '!' && echo "!   AUTOGENERATED FILE, DON'T MODIFY!" &&\
		echo '!' && echo '!    Base file: $(patsubst $(srcdir)/%,%,$<)' && \
		echo '!' && \
		$(M4) -I "$(top_srcdir)/m4" $<) \
			>$@.temp \
		&& mv $@.temp "$(srcdir)/$@"

CLEANFILES += $(GENERATED_F90)
endif

if USE_MPI
AM_FCFLAGS += $(FPP_DEFOPT)USE_MPI
if USE_MPI_MOD
AM_FCFLAGS += $(FPP_DEFOPT)USE_MPI_MOD
endif

libscalesppm_la_SOURCES += \
	core/errhandler.c \
	core/errhandler_cf.c \
	core/errhandler.h \
	core/ppm_xmpi.c \
	core/ppm_extents_mp.f90 \
	core/ppm_extents_mp_c.c \
	core/ppm_extents_mp_cf.c \
	core/ppm_extents_mp.h \
	core/ppm_std_type_kinds_mp.f90 \
	core/ppm_std_type_kinds_mp_c.c \
	core/ppm_std_type_kinds_mp_cf.c \
	core/ppm_std_type_kinds_mp.h \
	core/ppm_math_extensions_ddp_mp_c.c \
	\
	ppm/ppm_distributed.f90

libscalesppmcore_la_SOURCES += \
	core/errhandler.c \
	core/errhandler.h \
	core/ppm_xmpi.c \
	core/ppm_std_type_kinds_mp_c.c \
	core/ppm_std_type_kinds_mp.h \
	core/ppm_math_extensions_ddp_mp_c.c \
	core/ppm_extents_mp_c.c \
	core/ppm_extents_mp.h

if ENABLE_DISTRIBUTED_ARRAY_C
libscalesppmcore_la_SOURCES += \
	ppm/dist_array_c.c \
	ppm/dist_array.h

libscalesppm_la_SOURCES += \
	ppm/dist_array_c.c \
	ppm/dist_array.h \
	ppm/dist_array_f2c.c

libscalesppm_la_SOURCES += \
	ppm/dist_array.f90

endif
endif

if HAVE_POINTER_REMAP
AM_FCFLAGS += $(FPP_DEFOPT)HAVE_POINTER_REMAP
else
libscalesppm_la_SOURCES += \
	ppm/ppm_ptr_bnds_remap.f90
endif

if USE_PARMETIS
AM_FCFLAGS += $(FPP_DEFOPT)USE_PARMETIS
libscalesppm_la_SOURCES += ppm/ppm_graph_partition_mpi.f90
if HAVE_PARMETIS_V3
AM_FCFLAGS += $(FPP_DEFOPT)HAVE_PARMETIS_V3
libscalesppm_la_SOURCES += ppm/parmetis_wrap.c
endif
endif
if USE_METIS
AM_FCFLAGS += $(FPP_DEFOPT)USE_METIS
if HAVE_METIS_V4
AM_FCFLAGS += $(FPP_DEFOPT)HAVE_METIS_V4
endif
libscalesppm_la_SOURCES += ppm/ppm_graph_partition_serial.f90
endif
if USE_FC_NETCDF
AM_FCFLAGS += $(NETCDF_FC_INCLUDE)
libscalesppm_la_SOURCES += core/ppm_ncdf_dump.f90
if USE_FC_NETCDF_MOD
AM_FCFLAGS += $(FPP_DEFOPT)USE_NETCDF_MOD
libscalesppm_la_SOURCES += core/ppm_ncdf_dump_wrappers.f90
endif
endif
if IEEE_ARITHMETIC_EMULATION
libscalesppm_la_SOURCES += core/ieee_arithmetic.f90 \
	core/ieee_arithmetic_emulation.c
endif

if ! USE_CRYPTO

libscalesppm_la_SOURCES += \
	crypto/md5.c crypto/md5.h

libscalesppmcore_la_SOURCES += \
	crypto/md5.c crypto/md5.h

endif

libscalesppm_la_LIBADD = $(PARMETIS_FC_LIB) $(METIS_FC_LIB) \
	$(NETCDF_FC_LIB) $(MPI_FC_LIB) $(CRYPTO_FC_LIB)

libscalesppmcore_la_LIBADD = $(MPI_C_LIB) $(CRYPTO_C_LIB) $(LIBM)

if with_fortran
noinst_LTLIBRARIES = libtestdata.la

libtestdata_la_SOURCES = testdata/random_data.f90
libtestdata_la_LDFLAGS = -static $(AM_LDFLAGS)
libtestdata_la_LIBADD = libscalesppm.la
endif

pkgconfigdir = $(libdir)/pkgconfig

nodist_pkgconfig_DATA = pkgconfig/scales-ppm-core.pc

if with_fortran
nodist_pkgconfig_DATA += pkgconfig/scales-ppm.pc
endif

DISTCLEANFILES = $(nodist_pkgconfig_DATA)

clean-local:
	-rm -f *.$(FCMODEXT) */*.$(FCMODEXT) i.*.L *.i */*.i ??????????.dbg
	-rm -rf *.dSYM */*.dSYM

# When an out-of-source build takes place, the recipe of the rule expects that
# the object files are already deleted from the 'core' subdirectory, which is
# not always the case for a parallel run, i.e. when 'make distclean -j2'.
# Therefore, we explicitly specify a dependency on the 'mostlyclean' rule.
# Moreover, the 'distclean' recipe of the top-level Makefile expects the
# directory containing this Makefile to be empty as soon as its 'distclean'
# recipe is executed. This is not the case when the Makefile is generated with
# Automake 1.16 or later because the $(DEPDIR) subdirectories are not deleted.
# Therefore, we delete them explicitly as the older versions of Automake do.
distclean-local: mostlyclean
	-rm -rf '$(DEPDIR)' *'/$(DEPDIR)' \
		core/'$(DEPDIR)' \
		core/ppm_real_sp_dp_edit_descriptor.inc
	-rmdir core crypto pkgconfig ppm solver testdata || exit 0

EXTRA_DIST = main_doc.dox ppm/howto.dox solver/howto.dox \
	ppm/dist_array.f90.m4 ppm/ppm_ptr_bnds_remap.f90.m4 \
	m4.deps


makef90depend_ignore_args = \
	--ignore-module=netcdf \
        --ignore-module=ieee_arithmetic \
        --ignore-module=iso_c_binding \
        --ignore-module=mpi \
        --ignore-module=omp_lib

if with_fortran
./$(DEPDIR)/FC.deps: $(libscalesppm_la_SOURCES) \
	$(nodist_libscalesppm_la_SOURCES) \
	$(libtestdata_la_SOURCES) \
	core/$(am__dirstamp) crypto/$(am__dirstamp) \
	ppm/$(am__dirstamp) solver/$(am__dirstamp) \
	Makefile ../scripts/makef90depends
	$(AM_V_GEN)unique=`for i in $(libscalesppm_la_SOURCES) \
	  $(nodist_libscalesppm_la_SOURCES) $(libtestdata_la_SOURCES); do \
	  if test -f "$$i"; then echo "$$i"; else echo "$(srcdir)/$$i"; fi; \
	  done | grep '\.f90$$' ` ; \
	  FC="$(FC)" FPP="$(FPP)" \
	  $(PERL) -I$(top_srcdir)/scripts \
	  $(top_srcdir)/scripts/makef90depends \
	  $(makef90depend_ignore_args) \
	  --fc-mod-flag=$(FC_MOD_FLAG) --fc-def-opt=$(FPP_DEFOPT) \
	  --fc-inc-flag=$(FPP_INCOPT) --fc-mod-suf=.$(FCMODEXT) \
	  --fc-mod-case=$(FCMODCASE) --prog-fpp='$(FPP)' \
	  --fpp-sizeof-int=$(C_SIZEOF_INT) \
	  --fpp-sizeof-long=$(C_SIZEOF_LONG) \
	  --fpp-sizeof-long-long=$(C_SIZEOF_LONG_LONG) \
	  --fpp-char-is-`test $(C_CHAR_IS_UNSIGNED) = 0 || echo un`signed \
	  --parse-cache=./$(DEPDIR)/FC.depdb \
	  --obj-transform-name="^`printf '%s' '$(srcdir)/' \
		| $(PERL) -e 'print quotemeta(<>);'`=" \
	  --obj-transform-name="\.(?:f9[05]|for|f|fpp|ftn|f0[38])$$=.lo" \
	  -- \
	  $(AM_FCFLAGS) $(FCFLAGS) $(FPPFLAGS) -- $$unique \
	  >./$(DEPDIR)/FC.deps.temp && \
	cmp ./$(DEPDIR)/FC.deps.temp ./$(DEPDIR)/FC.deps >/dev/null 2>&1 \
	  || ln -f ./$(DEPDIR)/FC.deps.temp ./$(DEPDIR)/FC.deps

@am__include@ @am__quote@./$(DEPDIR)/FC.deps@am__quote@ # am--include-marker
endif

html: ppm_doxygen_deps
all: ppm_doxygen_deps

.PHONY: ppm_doxygen_deps
ppm_doxygen_deps:
	$(AM_V_GEN){ top_srcdir=$(top_srcdir) srcdir=$(srcdir) ; \
		export top_srcdir srcdir ; \
		doxygen_FCFLAGS="$(AM_FCFLAGS)" ; \
		doxygen_FCFLAGS=`echo x" $${doxygen_FCFLAGS}" \
		| $(PERL) -w -p -e 'substr($$_,0,1,"");chomp;' \
		-e 'my ($$tsdir,$$sdir)=(quotemeta($$ENV{'top_srcdir'}),quotemeta($$ENV{'srcdir'}));' \
		-e 's/ $(FPP_INCOPT)$${tsdir}/ \\$$(FPP_INCOPT)\\$$(top_srcdir)/g;s/ $(FPP_INCOPT)$${sdir}/ \\$$(FPP_INCOPT)\\$$(top_srcdir)\\/src/g;' \
		-e 's/ $(FPP_INCOPT)\\.\\.\\// \\$$(FPP_INCOPT)\\$$(top_builddir)\\//g;' \
		-e 's/ $(FPP_INCOPT)([^\\$$\\/])/ \\$$(FPP_INCOPT)\\$$(top_builddir)\\/src\\/$$1/g;' \
		-e 's/^\\s+//;s/\\s+$$//;'` && \
		find $(srcdir) -type f -a ! \( -name 'Makefile.*' \
		-o -name Makefile -o -name '#*' \
		-o -name '*\(\.in\|\.m4\|\.list\|\.o\|\.lo\|~\)' \) \
		| sed -e \
		'1s/^/ppm_doxygen_deps = /;1!s/^/	/;$$q;s/$$/ \\/' \
		&& echo 'doxygen_FCFLAGS = '"$${doxygen_FCFLAGS}" ; } \
		>$@.list.temp \
		&& { cmp $@.list.temp $@.list >/dev/null 2>&1 \
			&& rm -f $@.list.temp || mv $@.list.temp $@.list ; }


#
# Local Variables:
# license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
# license-default: "bsd"
# End:
