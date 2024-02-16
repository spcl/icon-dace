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

xtincludedir = $(includedir)/xt

nodist_include_HEADERS = yaxt.h

SUBDIRS = . f90
DIST_SUBDIRS = f90 .

../src/xtinclude_Makefile.inc:
	(cd ../src ; if test -r Makefile ; then \
	  $(MAKE) xtinclude_Makefile.inc ; fi)

@am__include@ @am__quote@../src/xtinclude_Makefile.inc@am__quote@ # am--include-marker

yaxt.h: Makefile $(xtinclude_HEADERS)
	$(AM_V_GEN)( echo '/* license */' ; echo '#ifndef YAXT_H' ; echo '#define YAXT_H' ; \
	for header in $(xtinclude_HEADERS) ; \
	  do echo "#include <$$header>" ; done ; \
	echo '#endif' \
	) >$@.$$$$.tmp ; \
	if test ! -f $@ || ! cmp $@ $@.$$$$.tmp >/dev/null ; then \
	mv $@.$$$$.tmp $@ ; else rm $@.$$$$.tmp ; fi

xt/.dir_stamp:
	$(AM_V_GEN)$(MKDIR_P) xt && touch xt/.dir_stamp

$(xtinclude_HEADERS): ../src/$@ ../include/config.h \
	$(top_srcdir)/scripts/header2installheader.pl xt/.dir_stamp
	$(AM_V_GEN) headerbasedir=`test -f '../src/$@' || echo '$(srcdir)/'` ; \
	  $(PERL) -I$(top_srcdir)/scripts \
	  $(top_srcdir)/scripts/header2installheader.pl \
	  --config-header=../include/config.h \
	  --srcdir=$$headerbasedir../src \
	  --c-sizeof-int=$(C_SIZEOF_INT) \
	  --c-sizeof-long=$(C_SIZEOF_LONG) \
	  --c-sizeof-long-long=$(C_SIZEOF_LONG_LONG) \
	  --c-char-is-`test $(C_CHAR_IS_UNSIGNED) = 0 || echo un`signed \
	  --dstdir=. \
	  $$headerbasedir../src/$@

clean-local:
	rm -rf xt yaxt.h

distclean-local:
	-test "$(abs_builddir)" = "$(abs_srcdir)" || rmdir f90

install-data-local: yaxt.h $(xtinclude_HEADERS)
	@echo " $(MKDIR_P) '$(DESTDIR)$(xtincludedir)'"; \
          $(MKDIR_P) "$(DESTDIR)$(xtincludedir)" || exit 1; \
	for p in $(xtinclude_HEADERS); do \
          echo " $(INSTALL_HEADER) $$p '$(DESTDIR)$(xtincludedir)'"; \
          $(INSTALL_HEADER) $$p "$(DESTDIR)$(xtincludedir)" || exit $$?; \
        done

uninstall-local:
	rm -rf "$(DESTDIR)$(xtincludedir)"

#
# Local Variables:
# coding: utf-8
# show-trailing-whitespace: t
# End:
#
