AM_CFLAGS = -I$(top_srcdir)/src $(MPI_CFLAGS)
LDADD = $(top_builddir)/src/libyac.a \
        $(LAPACK_CLIBS) \
        $(MTIME_CLIBS) \
        $(FYAML_CLIBS) \
        $(NETCDF_CLIBS) \
        $(YAXT_CLIBS) \
        $(MPI_CLIBS)

all-local: python-bindings

setup_env = \
  CC='$(CC)' \
  LDSHARED='' \
  CPPFLAGS='$(CPPFLAGS)' \
  CFLAGS='$(AM_CFLAGS) $(CFLAGS)' \
  LDFLAGS='$(LDFLAGS)' \
  LIBS='$(LDADD) $(LIBS)'

.PHONY: python-bindings
python-bindings: setup.py
	$(AM_V_GEN)$(setup_env) $(PYTHON) $< `test -z '$(AM_V_GEN)' && echo '-v' || echo '-q'` build --build-lib $(builddir)

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
	-rm -f yac.c

clean-local:
	-rm -rf tests/seq_noise
	-rm -rf tests/par_noise
	-rm -rf tests/seq_nc_noise
	-rm -rf tests/par_nc_noise
	-rm -f tests/seq_noise.nc
	-rm -f tests/par_noise.nc
	-rm -rf tests/noise_0
	-rm -rf build yac.egg-info yac

TEST_EXTENSIONS = .sh

EXTRA_DIST = \
        __init__.py \
        examples/__init__.py \
        examples/driver.py \
        examples/netcdf_reader.py \
        examples/netcdf_writer.py \
        examples/noisegenerator.py \
        examples/plotter.py \
        utils/__init__.py \
        utils/grid_utils.py

TESTS = \
        tests/seq_plot_noise_test.sh \
        tests/par_plot_noise_test.sh \
        tests/seq_noise2nc_test.sh \
        tests/par_noise2nc_test.sh \
        tests/seq_nc2plot_test.sh \
        tests/par_nc2plot_test.sh \
        tests/multithreading.sh \
        tests/test_interface.sh

tests/seq_nc2plot_test.log: tests/seq_noise2nc_test.log
if !ENABLE_CONCURRENT_MPI_TESTS
# Run MPI tests ony by one:
tests/par_nc2plot_test.log: tests/par_noise2nc_test.log
tests/par_noise2nc_test.log: tests/par_plot_noise_test.log
tests/multithreading.log: tests/par_noise2nc_test.log
endif

# The test script reside in subdirectories that do not have their own makefiles
# and Automake generates the update rules for them in the root Makefile.
# However, the root Makefile does not know that we need the files for the
# testing and does not update them before running. Therefore, we have to
# duplicate the rules here (with some modifications, see the automatically
# generated rules in $(top_srcdir)/tests/Makefile.in):
tests/seq_plot_noise_test.sh: $(top_builddir)/config.status $(srcdir)/tests/seq_plot_noise_test.sh.in
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@
tests/par_plot_noise_test.sh: $(top_builddir)/config.status $(srcdir)/tests/par_plot_noise_test.sh.in
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@
tests/seq_noise2nc_test.sh: $(top_builddir)/config.status $(srcdir)/tests/seq_noise2nc_test.sh.in
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@
tests/par_noise2nc_test.sh: $(top_builddir)/config.status $(srcdir)/tests/par_noise2nc_test.sh.in
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@
tests/seq_nc2plot_test.sh: $(top_builddir)/config.status $(srcdir)/tests/seq_nc2plot_test.sh.in
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@
tests/par_nc2plot_test.sh: $(top_builddir)/config.status $(srcdir)/tests/par_nc2plot_test.sh.in
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@
tests/multithreading.sh: $(top_builddir)/config.status $(srcdir)/tests/multithreading.sh.in
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@
