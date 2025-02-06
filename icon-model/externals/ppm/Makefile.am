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
SUBDIRS = . src include example test doc
DIST_SUBDIRS = src include example test doc .

EXTRA_DIST = util/add-f90-module \
	util/naglinkhack.pl \
	util/sunf95preproc-wrapper \
	util/sxpreproc-wrapper \
	util/xlfpreproc-wrapper \
	util/pbs_aprun util/crayftnpreproc-wrapper \
	util/mpi_launch_poe \
	util/mpi_launch_bgrunjob \
	util/icomp15/ld \
	scripts/build-makef90depends-dist.sh \
	scripts/cpp.pm \
	scripts/makef90depends \
	scripts/header2installheader.pl \
	scripts/reconfigure \
	scripts/reconfigure-mistral \
	scripts/recreate-testsuite.sh \
	scripts/generate-testsuite-m4-input \
	m4/forloop2.m4 \
	config/checksrc/basic_mpi_rma.c \
	config/checkdoc/basic_mpi_rma.txt \
	config/checksrc/openmpi-typedup-rma.c \
	config/checkdoc/openmpi-typedup-rma.txt \
	config/checksrc/incorrect_mpi_abort_exitcode.c \
	config/checkdoc/incorrect_mpi_abort_exitcode.txt \
	contrib/01aix-deplib-rpath-patch/aix-deplib-libtool.patch \
	contrib/01aix-deplib-rpath-patch/README.txt \
	contrib/01aix-deplib-rpath-patch/README \
	contrib/00nagfor-libtool-patch/README.txt \
	contrib/00nagfor-libtool-patch/nagfor-libtool-2.4.6.patch \
	contrib/00nagfor-libtool-patch/nagfor-libtool.patch \
	contrib/00nagfor-libtool-patch/nagfor-libtool-2.4.2.patch \
	contrib/libtool-patch-gen.sh \
	contrib/02nagfor53-shared-patch/README.txt \
	contrib/02nagfor53-shared-patch/nagfor53-shared.patch \
	contrib/03ltmain-ld-groups-patch/README \
	contrib/03ltmain-ld-groups-patch/ltmain-ld-groups-libtool-2.4.6.patch \
	contrib/03ltmain-ld-groups-patch/ltmain-ld-groups-libtool-2.4.2.patch \
	contrib/04ltmain-xlinker-patch/README \
	contrib/04ltmain-xlinker-patch/ltmain-xlinker-patch.patch \
	contrib/05macos-nagfor-patch/macos-nagfor.patch \
	contrib/06ltmain_nag_pthread-patch/ltmain_nag_pthread.patch \
	contrib/07ltmain-early-xcompile-patch/README \
	contrib/07ltmain-early-xcompile-patch/ltmain-early-xcompile-libtool-2.4.6.patch \
	LICENSE

ACLOCAL_AMFLAGS=-I m4

clean-local:
	@for d in $(SUBDIRS); do \
	  test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile || \
	  ./config.status $$d/Makefile ; done
	-rm -f *.$(FCMODEXT)

# The recipe of the rule needs to be called after 'clean-local' because the
# latter regenerates Makefiles. In order to ensure the order in the case of a
# parallel run, i.e. when 'make distclean -j2', we express the dependency
# explicitly:
distclean-recursive: clean-local
distclean-local: clean-local
	-rm -f src/ppm_doxygen_deps.list
	@for d in $(SUBDIRS); do \
		dd=$$d/$(DEPDIR) ; test ! -d $$dd \
		|| (rm -f $$dd/FC.dep* ; rmdir $$dd ) ; \
		test $$d = . || { rm -f $$d/Makefile ; \
			[ -r $$d/Makefile.in ] || rmdir $$d ; } ; done

all-local:
	@for d in $(SUBDIRS); do \
	  test $$d = . -o $(srcdir)/$$d/Makefile.in -ot $$d/Makefile || \
	  ./config.status $$d/Makefile ; done

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
