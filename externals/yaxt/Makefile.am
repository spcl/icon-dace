#
# Makefile.am --- Makefile for yaxt
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
SUBDIRS = . src inst_headers tests examples perf
DIST_SUBDIRS = src inst_headers tests examples perf .

ACLOCAL_AMFLAGS=-I m4

EXTRA_DIST = util/sunf95preproc-wrapper \
	util/sxpreproc-wrapper util/xlfpreproc-wrapper \
	util/crayftnpreproc-wrapper \
	util/pbs_aprun \
	util/mpi_launch_poe \
	util/mpi_launch_bgrunjob \
	util/icomp15/ld \
	util/nagdynlib/ld \
	config/checksrc/intelmpi_datatype.c \
	config/checkdoc/intelmpi_datatype.txt \
	config/checksrc/openmpi_struct_dt.c \
	config/checkdoc/openmpi_struct_dt.txt \
	config/checksrc/openmpi_datatype.c \
	config/checkdoc/openmpi_datatype.txt \
	config/checksrc/openmpi_dup.c \
	config/checkdoc/openmpi_dup.txt \
	config/checksrc/pe5228_iv32584.c \
	config/checkdoc/pe5228_iv32584.txt \
	config/checksrc/incorrect_mpi_abort_exitcode.c \
	config/checkdoc/incorrect_mpi_abort_exitcode.txt \
	config/checksrc/mpich_3.4_yaksa_hindexed.c \
	config/checkdoc/mpich_3.4_yaksa_hindexed.txt \
	config/checkpatch/mpich_3.4_yaksa_hindexed.patch \
	config/checkpatch/openmpi_datatype.opal_datatype_add.c.patch \
	config/checkpatch/openmpi_datatype.opal_datatype_optimize.c.patch \
	config/checksrc/mpich_unweighted_topo_comm_dup.c \
	config/checkdoc/mpich_unweighted_topo_comm_dup.txt \
	config/checkpatch/mpich_3.x_unweighted.patch \
	config/checkdoc/mpich_neighbor_alltoallw.txt \
	config/checksrc/mpich_neighbor_alltoallw.c \
	config/checkpatch/mpich_neighbor_alltoallw.patch \
	config/checkdoc/mpich_zero_size_contiguous.txt \
	config/checkpatch/mpich_zero_size_contiguous.patch \
	config/checksrc/mpich_zero_size_contiguous.c \
	config/checkdoc/mpich_4.0.x_ddt_transfer.txt \
	config/checkpatch/mpich_4.0.x_ddt_transfer.patch \
	config/checksrc/mpich_4.0.x_ddt_transfer.c \
	config/checkdoc/mpi_fortran_startup.txt \
	config/checksrc/mpi_fortran_startup.f90 \
	contrib/00nagfor-libtool-patch/README.txt \
	contrib/00nagfor-libtool-patch/nagfor-libtool-2.4.2.patch \
	contrib/00nagfor-libtool-patch/nagfor-libtool-2.4.6.patch \
	contrib/00nagfor-libtool-patch/nagfor-libtool-2.4.7.patch \
	contrib/01aix-deplib-rpath-patch/README \
	contrib/01aix-deplib-rpath-patch/aix-deplib-libtool.patch \
	contrib/02nagfor53-shared-patch/README.txt \
	contrib/02nagfor53-shared-patch/nagfor53-shared.patch \
	contrib/03ltmain-ld-groups-patch/README \
	contrib/03ltmain-ld-groups-patch/ltmain-ld-groups-libtool-2.4.6.patch \
	contrib/03ltmain-ld-groups-patch/ltmain-ld-groups-libtool-2.4.2.patch \
	contrib/04ltmain-xlinker-patch/README \
	contrib/04ltmain-xlinker-patch/ltmain-xlinker-patch.patch \
	contrib/05macos-nagfor-patch/macos-nagfor-libtool-2.4.2.patch \
	contrib/05macos-nagfor-patch/macos-nagfor.patch \
	contrib/06ltmain_nag_pthread-patch/ltmain_nag_pthread.patch \
	contrib/06ltmain_nag_pthread-patch/ltmain_nag_pthread-libtool-2.4.2.patch \
	contrib/07ltmain-early-xcompile-patch/README \
	contrib/07ltmain-early-xcompile-patch/ltmain-early-xcompile-libtool-2.4.6.patch \
	contrib/08ltmain-parallel-printf-patch/ltmain-parallel-printf-libtool-2.4.2.patch \
	contrib/08ltmain-parallel-printf-patch/ltmain-parallel-printf-libtool-2.4.6.patch \
	contrib/09debian-no-overlink-patch/README.txt \
	contrib/09debian-no-overlink-patch/deplib_binary-libtool-2.4.2.patch \
	contrib/09debian-no-overlink-patch/deplib_binary.patch \
	contrib/09debian-no-overlink-patch/link_all_deplibs-libtool-2.4.2.patch \
	contrib/09debian-no-overlink-patch/link_all_deplibs.patch \
	doc/Makefile \
	doc/Doxyfile.in doc/README.build doc/cfortran.doc \
	doc/cfortran.html doc/cfortran.txt doc/styleguide.txt \
	doc/images/decomp_idxlist_collection.png \
	doc/images/decomp_idxlist_collection.svg \
	doc/images/decomp_idxmod.png \
	doc/images/decomp_idxsection.png \
	doc/images/decomp_idxsection.svg \
	doc/images/decomp_idxstripes.png \
	doc/images/decomp_idxstripes.svg \
	doc/images/decomp_idxvec.png \
	doc/images/decomp_idxvec.svg \
	doc/images/neg_size_idxsection_nnnn.png \
	doc/images/neg_size_idxsection_nnnp.png \
	doc/images/neg_size_idxsection_nnpn.png \
	doc/images/neg_size_idxsection_nnpp.png \
	doc/images/neg_size_idxsection_npnn.png \
	doc/images/neg_size_idxsection_npnp.png \
	doc/images/neg_size_idxsection_nppn.png \
	doc/images/neg_size_idxsection_nppp.png \
	doc/images/neg_size_idxsection_pnnn.png \
	doc/images/neg_size_idxsection_pnnp.png \
	doc/images/neg_size_idxsection_pnpn.png \
	doc/images/neg_size_idxsection_pnpp.png \
	doc/images/neg_size_idxsection_ppnn.png \
	doc/images/neg_size_idxsection_ppnp.png \
	doc/images/neg_size_idxsection_pppn.png \
	doc/images/neg_size_idxsection_pppp.png \
	doc/images/stride_idxsection_x.png \
	doc/images/stride_idxsection_xy.png \
	doc/src/idxlists.dox \
	doc/src/idxsection.dox \
	doc/src/mainpage.dox \
	doc/src/redist_colls.dox \
	doc/src/rrobin.dox \
	doc/src/xmaps.dox \
	doc/src/yaxt_concept.dox \
	doc/src/correct_use_of_addresses.dox \
	doc/src/multi-threading.dox \
	include/cfortran.h \
	m4/tj_find_type.m4 \
	m4/forloop2.m4 \
	scripts/cpp.pm \
	scripts/header2installheader.pl \
	scripts/generate_doc.pl \
	scripts/makef90depends \
	scripts/reconfigure \
	scripts/reconfigure-mistral \
	scripts/reconfigure-levante \
	LICENSE COPYING


all-local:
	@for d in $(SUBDIRS) ; do \
          test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile \
	    || ./config.status $$d/Makefile ; done

clean-local:
	@for d in $(SUBDIRS) ; do \
          test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile \
	    || ./config.status $$d/Makefile ; done

# The recipe of the rule needs to be called after 'clean-local' because the
# latter regenerates Makefiles. In order to ensure the order in the case of a
# parallel run, i.e. when 'make distclean -j2', we express the dependency
# explicitly:
distclean-local: clean-local
	-for d in $(SUBDIRS) ; do \
		test $$d = . || (rm -f $$d/Makefile ; test "$(abs_builddir)" = "$(abs_srcdir)" || rmdir $$d) ; done
	-rm -rf config/workarounds
	-test "$(abs_builddir)" = "$(abs_srcdir)" || rmdir include util config
install-exec-local:
	@for d in $(SUBDIRS); do \
	  test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile || \
	  ./config.status $$d/Makefile ; done

install-data-local:
	@for d in $(SUBDIRS); do \
	  test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile || \
	  ./config.status $$d/Makefile ; done

check-local:
	@for d in $(SUBDIRS); do \
	  test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile || \
	  ./config.status $$d/Makefile ; done

#
# Local Variables:
# mode: Makefile
# coding: utf-8
# license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
# license-default: "bsd"
# End:
#
