#
# Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
#
# Version: 1.0
# Author: Thomas Jahns <jahns@dkrz.de>
# Maintainer: Thomas Jahns <jahns@dkrz.de>
# URL: https://www.dkrz.de/redmine/projects/scales-ppm
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
EXTRA_DIST = BMBF-Gef-Logo-200.png \
	doxygen_extra.css \
	scales_logo2-web.png \
	footer.htmli

html-local: html/index.html

../../src/ppm_doxygen_deps.list:
	if test -f ../../src/Makefile ; then \
		cd ../../src ; $(MAKE) ppm_doxygen_deps ; fi

@am__include@ @am__quote@../../src/ppm_doxygen_deps.list@am__quote@ # am--include-marker

html/index.html: scales-ppm.doxygen html/BMBF-Gef-Logo-200.png \
	$(ppm_doxygen_deps) ../../src/ppm_doxygen_deps.list
	doxygen_FCFLAGS='$(doxygen_FCFLAGS)' $(DOXYGEN) scales-ppm.doxygen
if ! HAVE_DOXYGEN_HTML_EXTRA_STYLESHEET
	cat $(srcdir)/doxygen_extra.css >>html/doxygen.css
endif

html/BMBF-Gef-Logo-200.png: $(srcdir)/BMBF-Gef-Logo-200.png
	mkdir -p html ; cp $(srcdir)/BMBF-Gef-Logo-200.png html/

clean-local:
	rm -rf html latex
#
# Local Variables:
# license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
# license-default: "bsd"
# End:
