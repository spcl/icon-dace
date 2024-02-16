#
# Copyright  (C)  2010-2016  Thomas Jahns <jahns@dkrz.de>
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
AM_FCFLAGS = $(FPP_INCOPT)$(top_srcdir)/include/f90 \
	$(FPP_INCOPT)../include/f77 \
	$(FC_MOD_FLAG)../include/f90 $(MPI_FC_INCLUDE) \
	$(FC_MOD_FLAG)../example

AM_CPPFLAGS = -I$(top_srcdir)/include $(MPI_C_INCLUDE)

check_PROGRAMS = \
	test_uniform_partition_c \
	test_combinatorics_c \
	test_ddp_math_c \
	test_factorization

if with_fortran
check_PROGRAMS += \
	test_qsort \
	test_insertion_sort \
	test_bsearch \
	test_uniform_partition \
	test_irand \
	test_strio test_posix_f test_extents test_sparse_mask_index \
	test_combinatorics \
	test_compact_mask_index \
	test_graph_csr test_ieee_emulation \
	test_ddp_math \
	test_checksum \
	test_abort_ppm \
	test_heap \
	test_set_partition test_set_partition_base \
	test_set_repartition \
	test_m1d
endif

if HAVE_POINTER_REMAP
AM_FCFLAGS += $(FPP_DEFOPT)HAVE_POINTER_REMAP
endif

if USE_MPI
AM_FCFLAGS += $(FPP_DEFOPT)USE_MPI
if USE_MPI_MOD
AM_FCFLAGS += $(FPP_DEFOPT)USE_MPI_MOD
endif
if HAVE_POINTER_REMAP
AM_FCFLAGS += $(FPP_DEFOPT)HAVE_POINTER_REMAP
endif

check_PROGRAMS += \
	test_ddp_math_mp
if with_fortran
check_PROGRAMS+= test_strided_extents test_set_repartition_mp \
	test_graph_dist_mp
endif

if ENABLE_DISTRIBUTED_ARRAY_C
check_PROGRAMS += \
	test_distributed_array_c_mp
if with_fortran
check_PROGRAMS += \
	test_distributed_array_mp
endif
endif
endif

if ! with_on_demand_check_programs
all: $(check_PROGRAMS)
endif

LDADD = ../src/libtestdata.la ../src/libscalesppm.la

test_m1d_SOURCES = test_m1d.f90

test_qsort_SOURCES = test_qsort.f90

test_insertion_sort_SOURCES = test_insertion_sort.f90

test_bsearch_SOURCES = test_bsearch.f90

test_checksum_SOURCES = test_checksum.f90

test_uniform_partition_SOURCES = test_uniform_partition.f90

test_uniform_partition_c_SOURCES = test_uniform_partition_c.c
test_uniform_partition_c_LDADD = ../src/libscalesppmcore.la

test_strided_extents_SOURCES = test_strided_extents.f90
test_strided_extents_LDADD = $(LDADD) $(MPI_FC_LIB)

test_irand_SOURCES = test_irand.f90

test_posix_f_SOURCES = test_posix_f.f90

test_strio_SOURCES = test_strio.f90 test_strio_c.c

test_extents_SOURCES = test_extents.f90

test_sparse_mask_index_SOURCES = test_sparse_mask_index.f90

test_combinatorics_SOURCES = test_combinatorics.f90

test_compact_mask_index_SOURCES = test_compact_mask_index.f90

test_graph_csr_SOURCES = test_graph_csr.f90

test_ieee_emulation_SOURCES = test_ieee_emulation.f90

test_ddp_math_SOURCES = test_ddp_math.f90

test_ddp_math_c_SOURCES = test_ddp_math_c.c
test_ddp_math_c_LDADD = ../src/libscalesppmcore.la

test_set_partition_SOURCES = test_set_partition.f90

test_set_repartition_SOURCES = test_set_repartition.f90

test_set_partition_base_SOURCES = test_set_partition_base.f90

test_set_repartition_mp_SOURCES = test_set_repartition_mp.f90
test_set_repartition_mp_LDADD = $(LDADD) $(MPI_FC_LIB)

test_distributed_array_mp_SOURCES = test_distributed_array_mp.f90
test_distributed_array_mp_LDADD = ../src/libscalesppm.la $(MPI_FC_LIB)

test_distributed_array_c_mp_SOURCES = test_distributed_array_c_mp.c
test_distributed_array_c_mp_LDADD = ../src/libscalesppmcore.la $(MPI_C_LIB)

test_abort_ppm_SOURCES = test_abort_ppm.f90

test_heap_SOURCES = test_heap.f90

test_graph_dist_mp_SOURCES = test_graph_dist_mp.f90
test_graph_dist_mp_LDADD = $(LDADD) $(MPI_FC_LIB)

test_factorization_SOURCES = test_factorization.c
test_factorization_LDADD = ../src/libscalesppmcore.la

test_combinatorics_c_SOURCES = test_combinatorics_c.c
test_combinatorics_c_LDADD = ../src/libscalesppmcore.la

test_ddp_math_mp_SOURCES = test_ddp_math_mp.c
test_ddp_math_mp_LDADD = ../src/libscalesppmcore.la $(MPI_C_LIB)

EXTRA_DIST=test_strided_extents.f90	\
	test_graph_dist_mp.f90	\
	test_ddp_math_mp.c \
	test_set_repartition_mp.f90	\
	test_distributed_array_mp.f90 \
	test_distributed_array_c_mp.c

if USE_FC_NETCDF
AM_FCFLAGS += $(NETCDF_FC_INCLUDE)
endif

# overrides to make regular compilations also use libtool
COMPILE = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(CC) -static $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)
FCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(FC) -static $(AM_FCFLAGS) $(FCFLAGS)

makef90depend_ignore_args = \
	--ignore-module=mpi \
        --ignore-module=ieee_arithmetic \
        --ignore-module=iso_c_binding \
        --ignore-module=omp_lib

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
	  $(AM_FCFLAGS) $(FCFLAGS) $(FPPFLAGS) -- $$unique \
	  >./$(DEPDIR)/FC.deps.temp && \
	ln -f ./$(DEPDIR)/FC.deps.temp ./$(DEPDIR)/FC.deps

@am__include@ @am__quote@./$(DEPDIR)/FC.deps@am__quote@ # am--include-marker
endif

clean-local:
	-rm -f *.$(FCMODEXT) i.*.L *.i atconfig ??????????.dbg
	-rm -rf *.dSYM
	test ! -f '$(TESTSUITE)' || \
	 $(SHELL) '$(TESTSUITE)' --clean

if MAINTAINER_MODE
# The `:;' works around a Bash 3.2 bug when the output is not writable.
$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	$(AM_V_GEN):;{ \
	  echo '# Signature of the current package.' && \
	  echo 'm4_define([AT_PACKAGE_NAME],' && \
	  echo '  [$(PACKAGE_NAME)])' && \
	  echo 'm4_define([AT_PACKAGE_TARNAME],' && \
	  echo '  [$(PACKAGE_TARNAME)])' && \
	  echo 'm4_define([AT_PACKAGE_VERSION],' && \
	  echo '  [$(PACKAGE_VERSION)])' && \
	  echo 'm4_define([AT_PACKAGE_STRING],' && \
	  echo '  [$(PACKAGE_STRING)])' && \
	  echo 'm4_define([AT_PACKAGE_BUGREPORT],' && \
	  echo '  [$(PACKAGE_BUGREPORT)])'; \
	  echo 'm4_define([AT_PACKAGE_URL],' && \
	  echo '  [$(PACKAGE_URL)])'; \
	} >'$(srcdir)/package.m4'
endif

EXTRA_DIST += testsuite.at $(srcdir)/package.m4 $(TESTSUITE) atlocal.in
TESTSUITE = $(srcdir)/testsuite

check-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' $(TESTSUITEFLAGS)

installcheck-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(bindir)' \
	  $(TESTSUITEFLAGS)

# We want to be able to distribute PPM through a git repository. For this
# purpose, we maintain a branch with all automatically generated files that are
# normally ignored in the master branch. By doing this, we minimize the number
# of tools that are required to build the package. For example, the users do
# not need to have the required version of Automake on their computers because
# the aforementioned branch already contains all the Makefile.in files. The
# problem that we have to solve is that git (for good reasons) does not track
# the timestamps of the files, which means that 'make' might try to rebuild the
# generated files even though it is not required, i.e. the contents of the
# prerequisites of the files are not changed but the timestamps of the
# prerequisites might be more recent than the timestamps of the generated
# files. In such case, the users might get error messages about missing
# utilities. In order to prevent this, we expand the AM_MAINTAINER_MODE in the
# configure.ac file, which, by default, disables the dependencies of the
# Makefiles on configure.ac, Makefile.am and other files that implement the
# building system. Those dependencies can be enabled by running the configure
# script with '--enable-maintainer-mode' argument. The same should apply to the
# testsuite files, i.e. to the $(TESTSUITE) and package.m4, which is why we
# modified the example presented in the Autotest documentation:
# https://www.gnu.org/software/autoconf/manual/autoconf-2.69/html_node/Making-testsuite-Scripts.html
# In order to be able build the master branch, where we don't have the
# generated files, as well as the release branch, where we have the generated
# files and don't want them to get updated when calling 'make', we need to
# meet the following requirements:
#   1. If the maintainer mode is disabled (default) and an automatically
#      generated file exists, the recipe of the rule for the file must not be
#      executed.
#   2. If the maintainer mode is enabled, an automatically generated file
#      exists but its prerequisite has a more recent timestamp, the recipe of
#      the rule must update the file.
#   3. If an automatically generated file does not exist, the recipe of the
#      rule must be executed regardless of whether the maintainer mode is
#      enabled or not and, if the required utility is found on the computer,
#      the file must be created.

AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): @MAINTAINER_MODE_TRUE@ $(srcdir)/testsuite.at $(srcdir)/package.m4
	$(AM_V_GEN)cd '$(srcdir)' ; \
		$(AUTOTEST) -I . -o testsuite.tmp testsuite.at \
		&& mv $@.tmp $@

atconfig:
	cd .. ; ./config.status test/atconfig

# When an out-of-source build takes place, the 'distclean' recipe of the
# top-level Makefile expects the directory containing this Makefile to be
# empty as soon as its 'distclean' recipe is executed. This is not the case
# when the Makefile is generated with Automake 1.16 or later because the
# $(DEPDIR) subdirectory is not deleted. Therefore, we delete it explicitly
# as the older versions of Automake do.
distclean-local:
	-rm -rf '$(DEPDIR)'
