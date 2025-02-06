# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

AM_TESTSUITE_SUMMARY_HEADER = ' for $(PACKAGE_STRING) ($(subdir))'
if ENABLE_OPENMP
AM_TESTS_ENVIRONMENT = OMP_NUM_THREADS=$${OMP_NUM_THREADS-1}; export OMP_NUM_THREADS;
endif

TEST_EXTENSIONS = .sh

TESTS = \
    test_coroutines.sh \
    test_interface.sh \
    test_multithreading.sh \
    test_nc2plot.sh \
    test_nc2plot_parallel.sh \
    test_noise2nc.sh \
    test_noise2nc_parallel.sh \
    test_noise2plot.sh \
    test_noise2plot_parallel.sh

if TEST_MPI
# Some tests require output files of the others:
test_nc2plot.log: test_noise2nc.log
test_nc2plot_parallel.log: test_noise2nc_parallel.log
if !ENABLE_CONCURRENT_MPI_TESTS
# Run MPI tests one by one:
test_noise2plot_parallel.log: test_noise2plot.log
test_noise2plot.log: test_nc2plot_parallel.log
test_noise2nc_parallel.log: test_nc2plot.log
test_noise2nc.log: test_multithreading.log
test_multithreading.log: test_interface.log
test_interface.log: test_coroutines.log
endif
endif

EXTRA_DIST = \
    test_coroutines.py \
    test_interface.py

clean-local:
	-rm -rf *.nc *.plots
