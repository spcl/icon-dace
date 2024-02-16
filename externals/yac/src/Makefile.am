#
# @file Makefile.am
#
# @copyright Copyright  (C)  2020 Moritz Hanke <hanke@dkrz.de>
#                                 Rene Redler <rene.redler@mpimet.mpg.de>
#                                 Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#                                 Teresa Holfeld <teresa.holfeld@zmaw.de>
#                                 Thomas Jahns <jahns@dkrz.de>
#
# @version 1.0
# @author Moritz Hanke <hanke@dkrz.de>
#         Rene Redler <rene.redler@mpimet.mpg.de>
#         Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#         Thomas Jahns <jahns@dkrz.de>
#
# Keywords:
# Maintainer: Moritz Hanke <hanke@dkrz.de>
#             Rene Redler <rene.redler@mpimet.mpg.de>
#             Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#             Thomas Jahns <jahns@dkrz.de>
# URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
#
# This file is part of YAC.
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

lib_LIBRARIES = libyac.a

include_HEADERS = yac_interface.h
nodist_include_HEADERS = yac_config.h

if FCMODUC
mo_yac_finterface_mod = MO_YAC_FINTERFACE.$(FCMODEXT)
else
mo_yac_finterface_mod = mo_yac_finterface.$(FCMODEXT)
endif

nodist_include_HEADERS += $(mo_yac_finterface_mod)

pkgconfigdir = $(libdir)/pkgconfig
dist_noinst_DATA = pkgconfig/yac.pc.in
nodist_pkgconfig_DATA = pkgconfig/yac.pc

libyac_a_SOURCES = \
        area.c \
        area.h \
        bnd_circle.c \
        bnd_triangle.c \
        check_overlap.c \
        clipping.c \
        clipping.h \
        component.c \
        component.h \
        config_yaml.c \
        config_yaml.h \
        core/core.c \
        core/core.h \
        core/ppm_xfuncs.h \
        core/symprefix.h \
        core/xmalloc.c \
        core/xstdio.c \
        couple_config.c \
        couple_config.h \
        dist_grid.c \
        dist_grid.h \
        dist_grid_utils.c \
        dist_grid_utils.h \
        ensure_array_size.c \
        ensure_array_size.h \
        event.c \
        event.h \
        fields.c \
        fields.h \
        geometry.h \
        grid.c \
        grid.h \
        grid_cell.c \
        grid_cell.h \
        grid_curve2d.c \
        grid_reg2d.c \
        grid_reg2d_common.c \
        grid_reg2d_common.h \
        grid_unstruct.c \
        instance.c \
        instance.h \
        interp_grid.c \
        interp_grid.h \
        interp_method.c \
        interp_method.h \
        interp_method_avg.c \
        interp_method_avg.h \
        interp_method_callback.c \
        interp_method_callback.h \
        interp_method_check.c \
        interp_method_check.h \
        interp_method_conserv.c \
        interp_method_conserv.h \
        interp_method_creep.c \
        interp_method_creep.h \
        interp_method_file.c \
        interp_method_file.h \
        interp_method_fixed.c \
        interp_method_fixed.h \
        interp_method_hcsbb.c \
        interp_method_hcsbb.h \
        interp_method_nnn.c \
        interp_method_nnn.h \
        interp_method_spmap.c \
        interp_method_spmap.h \
        interp_stack_config.c \
        interp_stack_config.h \
        interp_weights.c \
        interp_weights.h \
        interpolation.c \
        interpolation.h \
        interpolation_exchange.c \
        interpolation_exchange.h \
        interpolation_utils.c \
        interpolation_utils.h \
        interpolation_direct.c \
        interpolation_direct.h \
        interpolation_direct_mf.c \
        interpolation_direct_mf.h \
        interpolation_fixed.c \
        interpolation_fixed.h \
        interpolation_sum_mvp_at_src.c \
        interpolation_sum_mvp_at_src.h \
        interpolation_sum_mvp_at_tgt.c \
        interpolation_sum_mvp_at_tgt.h \
        intersection.c \
        interval_tree.c \
        interval_tree.h \
        io_utils.c \
        io_utils.h \
        mergesort.c \
        mo_yac_finterface.F90 \
        proc_sphere_part.c \
        proc_sphere_part.h \
        quicksort.c \
        quicksort_template.h \
        quicksort_template_2.h \
        sphere_part.c \
        sphere_part.h \
        utils.c \
        utils.h \
        version.h \
        yac_finterface.F90 \
        yac_interface.c \
        yac_interface.h \
        yac_lapack_interface.c \
        yac_lapack_interface.h \
        yac_mpi.c \
        yac_mpi.h \
        yac_redirstdout.c \
        feenableexcept.h

AM_CFLAGS = \
        $(PIC_CFLAGS) \
        $(MPI_CFLAGS) \
        $(YAXT_CFLAGS) \
        $(NETCDF_CFLAGS) \
        $(FYAML_CFLAGS) \
        $(MTIME_CFLAGS) \
        $(LAPACK_CFLAGS)

yac_finterface.$(OBJEXT): $(mo_yac_finterface_mod)
$(mo_yac_finterface_mod): mo_yac_finterface.$(OBJEXT)
	@if test ! -f $@; then rm -f $<; $(MAKE) $<; fi

mostlyclean-local:
	-rm -f *.$(FCMODEXT) i.*.L *.i

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
FCDEFS = $(FCPPDEF)HAVE_CONFIG_H
DEFAULT_FCINCLUDES = $(FCPPINC).@yac_fcisrc@
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_FCINCLUDES) $(FCINCLUDES) \
	$(AM_FCFLAGS) $(FCFLAGS)
