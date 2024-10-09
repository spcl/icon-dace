# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

AM_CPPFLAGS = -DYAC_CYTHON
AM_CFLAGS = \
    $(PIC_CFLAGS) \
    -I$(top_srcdir)/src/mci \
    $(MPI_CFLAGS) \
    $(PYTHON_CFLAGS)

LDADD = \
    $(top_builddir)/src/mci/libyac_mci.a \
    $(top_builddir)/src/core/libyac_core.a \
    $(OPENMP_CFLAGS) \
    $(LAPACK_CLIBS) \
    $(MTIME_CLIBS) \
    $(FYAML_CLIBS) \
    $(NETCDF_CLIBS) \
    $(YAXT_CLIBS) \
    $(MPI_CLIBS) \
    $(RPATH_CFLAGS)

# A dummy library just to trigger Automake to inject all required rules:
EXTRA_LIBRARIES = libdummypyyac.a
nodist_libdummypyyac_a_SOURCES = yac_pyx.c

all-local: python-bindings

setup_env = \
  CC='$(CC)' \
  LDFLAGS='$(LDFLAGS)' \
  LDADD='$(LDADD)' \
  LIBS='$(LIBS)' \
  YAC_OVERRIDE_PY_BUILD=$${YAC_OVERRIDE_PY_BUILD-yes}

.PHONY: python-bindings
python-bindings: setup.py yac_pyx.$(OBJEXT)
	@cmd="$(setup_env) $(PYTHON) $< `test -z '$(AM_V_GEN)' && echo '-v' || echo '-q'` build --build-lib $(builddir)"; \
	if test -n '$(AM_V_GEN)'; then \
	  echo "  GEN     " $@; \
	else \
	  echo "$$cmd"; \
	fi; \
	eval "$$cmd"

yac_pyx.c: yac.pyx
	@cmd="$(PYTHON) -m cython `test -z '$(AM_V_GEN)' && echo '-v'` -o $@ -3 --cleanup 3 --module-name='yac._yac' $<"; \
	if test -n '$(AM_V_GEN)'; then \
	  echo "  GEN     " $@; \
	else \
	  echo "$$cmd"; \
	fi; \
	eval "$$cmd"

pip_install_args = \
  --disable-pip-version-check \
  --no-cache-dir \
  --no-input \
  install \
  --ignore-installed \
  --no-build-isolation \
  --no-deps \
  --no-index \
  --no-warn-script-location

install-exec-local:
	$(setup_env) $(PYTHON) -m pip $(pip_install_args) . --prefix='$(DESTDIR)$(prefix)'

mostlyclean-local:
	-rm -f yac_pyx.c

clean-local:
	-rm -rf build yac.egg-info yac

EXTRA_DIST = \
    src/__init__.py \
    src/examples/__init__.py \
    src/examples/driver.py \
    src/examples/netcdf_reader.py \
    src/examples/netcdf_writer.py \
    src/examples/noisegenerator.py \
    src/examples/plotter.py \
    src/utils/__init__.py \
    src/utils/__init__.py \
    src/utils/grid_utils.py \
    yac.pyx
