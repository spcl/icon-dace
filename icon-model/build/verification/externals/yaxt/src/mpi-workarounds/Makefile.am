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

CFLAGS = @XT_WORKAROUND_MPI_CFLAGS@
CPPFLAGS = @XT_WORKAROUND_MPI_CPPFLAGS@

if XT_ADD_MPI_WORKAROUND
noinst_LTLIBRARIES = libmpi_workaround.la
libmpi_workaround_la_SOURCES =
libmpi_workaround_la_LIBADD = $(LIBOBJS) $(xt_mpi_workaround_LIBS) $(MPI_C_LIB)
endif

EXTRA_DIST = xt_xmlInitParser.c

distclean-local:
	-rm -rf '$(DEPDIR)'

#
# Local Variables:
# coding: utf-8
# show-trailing-whitespace: t
# End:
#
