# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

SUBDIRS = src grids

if ENABLE_MCI
if !EXTERNAL_MTIME
SUBDIRS += mtime
endif
endif ENABLE_MCI

if !EXTERNAL_LAPACK
SUBDIRS += clapack
endif

if ENABLE_PYTHON_BINDINGS
SUBDIRS += python
endif

if ENABLE_EXAMPLES
SUBDIRS += examples
endif

if ENABLE_TOOLS
SUBDIRS += tools
endif

SUBDIRS += tests

EXTRA_DIST = \
    config/checkdoc \
    config/checksrc \
    input/couple_config_test.yaml \
    input/coupling_test.yaml \
    input/coupling_test2.yaml \
    input/coupling_test3_comp_a_comp_b.yaml \
    input/coupling_test3_comp_a_comp_c.yaml \
    input/coupling_test3_comp_a_intra.yaml \
    input/coupling_test3_comp_b_comp_c.yaml \
    input/coupling_test3_comp_b_intra.yaml \
    input/coupling_test3_comp_c_intra.yaml \
    input/coupling_test3_default.yaml \
    input/coupling_test3_default_c.yaml \
    input/coupling_test4.yaml \
    input/coupling_test5.yaml \
    input/coupling_test6.json \
    input/coupling_test6_local.yaml \
    input/coupling_test9.yaml \
    input/cube_10x10.csv \
    input/instance_test_1_1.yaml \
    input/instance_test_1_10.yaml \
    input/instance_test_1_2.yaml \
    input/instance_test_1_3.yaml \
    input/instance_test_1_4.yaml \
    input/instance_test_1_5.yaml \
    input/instance_test_1_6.yaml \
    input/instance_test_1_7.yaml \
    input/instance_test_1_8.yaml \
    input/instance_test_1_9.yaml \
    input/test_def_datetime.yaml

# We do not have a clean way to uninstall the Python bindings yet:
distuninstallcheck_listfiles = find . \( -name site-packages -o -name dist-packages \) -prune -o -type f -print
