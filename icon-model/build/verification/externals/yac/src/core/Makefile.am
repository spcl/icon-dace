# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

lib_LIBRARIES = libyac_core.a

nodist_include_HEADERS = \
    yac_core.h

if FCMODUC
yac_core_mod = YAC_CORE.$(FCMODEXT)
else
yac_core_mod = yac_core.$(FCMODEXT)
endif

nodist_include_HEADERS += \
    $(yac_core_mod)

libyac_core_a_SOURCES = \
    area.c \
    area.h \
    basic_grid.c \
    basic_grid.h \
    basic_grid_data.c \
    basic_grid_data.h \
    bnd_circle.c \
    bnd_triangle.c \
    check_overlap.c \
    clipping.c \
    clipping.h \
    dist_grid.c \
    dist_grid.h \
    dist_grid_internal.h \
    ensure_array_size.c \
    ensure_array_size.h \
    feenableexcept.h \
    field_data.c \
    field_data.h \
    field_data_set.c \
    field_data_set.h \
    geometry.h \
    grid_cell.c \
    grid_cell.h \
    grid_curve2d.c \
    grid_cloud.c \
    grid_reg2d.c \
    grid_reg2d_common.c \
    grid_reg2d_common.h \
    grid_unstruct.c \
    interp_grid.c \
    interp_grid.h \
    interp_grid_internal.h \
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
    interp_method_internal.h \
    interp_method_ncc.c \
    interp_method_ncc.h \
    interp_method_nnn.c \
    interp_method_nnn.h \
    interp_method_spmap.c \
    interp_method_spmap.h \
    interp_stack_config.c \
    interp_stack_config.h \
    interp_weights.c \
    interp_weights.h \
    interp_weights_internal.h \
    interpolation.c \
    interpolation.h \
    interpolation_direct.c \
    interpolation_direct.h \
    interpolation_direct_mf.c \
    interpolation_direct_mf.h \
    interpolation_exchange.c \
    interpolation_exchange.h \
    interpolation_fixed.c \
    interpolation_fixed.h \
    interpolation_internal.h \
    interpolation_sum_mvp_at_src.c \
    interpolation_sum_mvp_at_src.h \
    interpolation_sum_mvp_at_tgt.c \
    interpolation_sum_mvp_at_tgt.h \
    interpolation_utils.c \
    interpolation_utils.h \
    intersection.c \
    interval_tree.c \
    interval_tree.h \
    io_utils.c \
    io_utils.h \
    location.c \
    location.h \
    mergesort.c \
    ppm/core.c \
    ppm/core.h \
    ppm/ppm_xfuncs.h \
    ppm/symprefix.h \
    ppm/xmalloc.c \
    ppm/xstdio.c \
    proc_sphere_part.c \
    proc_sphere_part.h \
    quicksort.c \
    quicksort_template.h \
    quicksort_template_2.h \
    remote_point.c \
    remote_point.h \
    sphere_part.c \
    sphere_part.h \
    utils_common.h \
    utils_core.c \
    utils_core.h \
    yac_assert.h \
    yac_core.F90 \
    yac_lapack_interface.c \
    yac_lapack_interface.h \
    yac_mpi.c \
    yac_mpi.h \
    yac_mpi_common.h \
    yac_mpi_internal.h \
    yac_types.h

yac_core_h_headers = \
    yac_mpi.h \
    yac_types.h \
    grid_cell.h \
    basic_grid_data.h \
    location.h \
    field_data.h \
    basic_grid.h \
    dist_grid.h \
    interp_grid.h \
    interp_method.h \
    interp_method_avg.h \
    interp_method_check.h \
    interp_method_conserv.h \
    interp_method_creep.h \
    interp_method_ncc.h \
    interp_method_nnn.h \
    interp_method_fixed.h \
    interp_method_file.h \
    interp_method_hcsbb.h \
    interp_method_spmap.h \
    interp_method_callback.h \
    interp_stack_config.h \
    interpolation.h \
    interp_weights.h \
    yac_assert.h \
    io_utils.h

AM_CPPFLAGS = \
    -I$(top_builddir)/src

AM_CFLAGS = \
    $(PIC_CFLAGS) \
    $(OPENMP_CFLAGS) \
    $(MPI_CFLAGS) \
    $(YAXT_CFLAGS) \
    $(NETCDF_CFLAGS) \
    $(LAPACK_CFLAGS)

AM_FCFLAGS = $(PIC_FCFLAGS)

$(yac_core_mod): yac_core.$(OBJEXT)
	@if test ! -f $@; then rm -f $<; $(MAKE) $(AM_MAKEFLAGS) $<; fi

all_local =

if ENABLE_DEPRECATED

if FCMODUC
mo_yac_core_mod = MO_YAC_CORE.$(FCMODEXT)
else
mo_yac_core_mod = mo_yac_core.$(FCMODEXT)
endif

nodist_include_HEADERS += \
    $(mo_yac_core_mod)

libyac_core_a_SOURCES += yac_core_deprecated.f90

yac_core_deprecated.$(OBJEXT): $(yac_core_mod)

$(mo_yac_core_mod): yac_core_deprecated.$(OBJEXT)
	@if test ! -f $@; then rm -f $<; $(MAKE) $(AM_MAKEFLAGS) $<; fi

all_local += manifest.txt

manifest.txt: $(libyac_core_a_OBJECTS)
	$(AM_V_GEN):;{ \
	  for obj in $(libyac_core_a_OBJECTS); do \
	    echo "$(abs_builddir)/$$obj"; \
	  done; \
	} >$@

MOSTLYCLEANFILES = manifest.txt

endif ENABLE_DEPRECATED

all_local += yac_core.h

# REUSE-IgnoreStart
yac_core.h: $(yac_core_h_headers)
	$(AM_V_GEN):;{ \
	echo '// Copyright (c) 2024 The YAC Authors'; \
	echo '//'; \
	echo '// SPDX-License-Identifier: BSD-3-Clause'; \
	echo ''; \
	echo '#ifndef YAC_CORE_H'; \
	echo '#define YAC_CORE_H'; \
	echo ''; \
	for f in $(yac_core_h_headers); do \
	  echo "// $$f"; \
	  sed '1,\|^// YAC PUBLIC HEADER START$$|d;\|^// YAC PUBLIC HEADER STOP$$|,$$d' `test -f "$$f" || echo '$(srcdir)/'`"$$f"; \
	  echo ''; \
	done; \
	echo '#endif // YAC_CORE_H'; \
	} >$@
# REUSE-IgnoreEnd

all-local: $(all_local)

mostlyclean-local:
	-rm -f *.$(FCMODEXT) yac_core.h

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
FCDEFS = $(FCPPDEF)HAVE_CONFIG_H
DEFAULT_FCINCLUDES = $(FCPPINC).@yac_fcisrc@ $(FCPPINC)$(top_builddir)/src
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_FCINCLUDES) $(FCINCLUDES) \
	$(AM_FCFLAGS) $(FCFLAGS)
