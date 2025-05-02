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

EXTRA_PROGRAMS = 						\
	rrobin							\
	idxlist_examples					\
	idxsection_examples					\
	row2col							\
	tr_gp2fs


if with_fortran
EXTRA_PROGRAMS += 						\
	row2col_f						\
	unstruct_halo_exchange_f
endif

noinst_PROGRAMS =
if with_example_programs
noinst_PROGRAMS += $(EXTRA_PROGRAMS)
endif
examples: $(EXTRA_PROGRAMS)

rrobin_SOURCES = rrobin.c

idxlist_examples_SOURCES = idxlist_examples.c print_index_list.h

idxsection_examples_SOURCES = idxsection_examples.c print_index_list.h

row2col_SOURCES = row2col.c

XT_FC_LDADD = ../src/libyaxt.la ../src/libyaxt_c.la $(MPI_FC_LIB)

row2col_f_SOURCES = row2col_f.f90
row2col_f_LDADD = $(XT_FC_LDADD)

tr_gp2fs_SOURCES = tr_gp2fs.c
tr_gp2fs_LDADD = $(LDADD) -lm

unstruct_halo_exchange_f_SOURCES = unstruct_halo_exchange_f.f90
unstruct_halo_exchange_f_LDADD = $(XT_FC_LDADD)

AM_CFLAGS = -I../inst_headers $(MPI_C_INCLUDE)
AM_FCFLAGS = $(FC_MOD_FLAG)../inst_headers/f90 \
	$(FPP_INCOPT)../include $(MPI_FC_MOD)
LDADD = ../src/libyaxt_c.la $(MPI_C_LIB)

# overrides to make regular compilations also use libtool
COMPILE = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(CC) -static $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)
FCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(FC) -static $(AM_FCFLAGS) $(FCFLAGS)

clean-local:
	-rm -f *.i *.dbg
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
# license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
# coding: utf-8
# license-default: "bsd"
# End:
