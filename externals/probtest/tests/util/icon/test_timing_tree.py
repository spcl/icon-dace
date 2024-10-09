"""
This module contains unit tests to verify the functionality of the TimingTree class,
including reading log files, loading from JSON, intersection and subtraction of trees,
growing trees with new nodes, and adding trees together.
"""

from datetime import datetime

import numpy as np
import pandas as pd

from util.icon.extract_timings import read_logfile
from util.tree import TimingTree

TIMING_FILE_1 = "tests/data/timing_example_1.txt"
TIMING_FILE_2 = "tests/data/timing_example_2.txt"
TIMING_FILE_3 = "tests/data/timing_example_3.txt"

JSON_REFERENCE = "tests/data/reference"
JSON_TREE_REFERENCE = "tests/data/reference_tree.json"
JSON_META_REFERENCE = "tests/data/reference_meta.json"
JSON_DATA_REFERENCE = "tests/data/reference_data.json"

JSON_ADD_REFERENCE = "tests/data/add"

pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)


def assert_trees_equal(t1, t2):
    for i in range(t1.meta_data["n_tables"]):
        diff = t1.data[i].values - t2.data[i].values
        diff_sum = np.sum(np.abs(diff))

        assert t1.meta_data == t2.meta_data, f"meta data does not match for table {i}"
        assert diff_sum < 1e-12, f"data does not match for table {i}"
        assert (
            t1.root[i].to_ancestry_name_list() == t2.root[i].to_ancestry_name_list()
        ), f"tree does not match for table {i}"


def test_read_timing():
    for timing_file in (TIMING_FILE_1, TIMING_FILE_2, TIMING_FILE_3):
        tt = TimingTree.from_logfile(timing_file, read_logfile)

        assert tt.data is not None, "did not properly initialize data"
        assert tt.meta_data is not None, "did not properly initialize meta data"
        assert tt.root is not None, "did not properly initialize tree"


def test_json_load():
    tt_json = TimingTree.from_json(JSON_REFERENCE)
    tt = TimingTree.from_logfile(TIMING_FILE_1, read_logfile)

    assert_trees_equal(tt_json, tt)


def test_intersection():
    tt1 = TimingTree.from_logfile(TIMING_FILE_1, read_logfile)
    tt2 = TimingTree.from_logfile(TIMING_FILE_2, read_logfile)

    intersection_nodes = tt1.root[-1].intersection(tt2.root[-1])

    names = {n.get_ancestry_name() for n in intersection_nodes}

    ref_names = {
        "root>total>integrate_nh>nh_solve>nh_solve.cellcomp",
        "root>upper_atmosphere",
        "root>total>integrate_nh>physics",
        "root>total>integrate_nh>physics>nwp_radiation>nwp_ecrad_store_fluxes",
        "root>exch_data",
        "root>nh_diagnostics",
        "root>total>integrate_nh>nh_solve>nh_solve.veltend",
        "root>total>integrate_nh>physics>nwp_radiation",
        "root>total>integrate_nh>physics>nwp_turbulence>nwp_turbtrans",
        "root>total>integrate_nh>physics>rediag_prog_vars",
        "root>total>integrate_nh>physics>cloud_cover",
        "root>total>integrate_nh>physics>nwp_turbulence>nwp_turbdiff",
        "root>total>write_restart",
        "root>total>integrate_nh>physics>sso",
        "root>total>integrate_nh>physics>nwp_radiation>nwp_radiation_upscale",
        "root>total>integrate_nh>nh_hdiff",
        "root>total>integrate_nh>transport>adv_horiz",
        "root>total>integrate_nh>transport>adv_horiz>adv_hflx",
        "root>total>integrate_nh>transport>adv_vert>adv_vflx",
        "root>total>integrate_nh>nesting>nesting.bdy_interp",
        "root>total>integrate_nh>physics>phys_acc_sync>global_sum",
        "root>total>integrate_nh>nh_solve>nh_solve.exch",
        "root>total>integrate_nh>physics>phys_acc_sync>ordglb_sum",
        "root>exch_data>exch_data.wait",
        "root>total>integrate_nh>nh_solve",
        "root>total>integrate_nh",
        "root>total>integrate_nh>physics>nwp_radiation>nwp_ecrad_utilities",
        "root>total",
        "root>total>integrate_nh>nh_solve>nh_solve.vimpl",
        "root>total>integrate_nh>physics>nwp_turbulence",
        "root>total>write_restart>write_restart_communication",
        "root>total>write_restart>write_restart_io",
        "root>total>integrate_nh>physics>diagnose_pres_temp",
        "root>total>integrate_nh>physics>nwp_microphysics",
        "root>total>integrate_nh>physics>nwp_surface",
        "root>total>integrate_nh>physics>nwp_convection",
        "root>total>integrate_nh>physics>nwp_radiation>preradiaton",
        "root>total>integrate_nh>nh_solve>nh_solve.vnupd",
        "root>total>integrate_nh>physics>nwp_radiation>nwp_radiation_downscale",
        "root>total>integrate_nh>nesting",
        "root>total>integrate_nh>transport>adv_vert",
        "root>total>integrate_nh>transport",
        "root>total>integrate_nh>physics>phys_u_v",
        "root>total>integrate_nh>physics>satad",
        "root>total>integrate_nh>physics>phys_acc_sync",
        "root>total>integrate_nh>physics>radheat",
        "root>total>integrate_nh>physics>nwp_radiation>nwp_ecrad_ecrad",
    }

    assert names == ref_names, "set of intersecting nodes does not match reference"


def test_sub():
    tt1 = TimingTree.from_logfile(TIMING_FILE_1, read_logfile)
    tt2 = TimingTree.from_logfile(TIMING_FILE_2, read_logfile)

    sub_nodes = tt1.root[-1].sub(tt2.root[-1])

    names = {n.get_ancestry_name() for n in sub_nodes}

    ref_names = {"root>total>integrate_nh>nh_solve>nh_solve.edgecomp"}

    assert names == ref_names, "set of differing nodes does not match reference"


def test_grow():
    tt1 = TimingTree.from_logfile(TIMING_FILE_1, read_logfile)
    tt2 = TimingTree.from_logfile(TIMING_FILE_2, read_logfile)

    diff_nodes = tt1.root[-1].sub(tt2.root[-1])

    tt2.grow(diff_nodes, -1)

    new_node = tt2.find("nh_solve.edgecomp", -1)

    assert (
        new_node.name == "nh_solve.edgecomp"
    ), "did not add (non-present) node nh_solve.edgecomp"


def test_add():
    tt1 = TimingTree.from_logfile(TIMING_FILE_1, read_logfile)
    tt2 = TimingTree.from_logfile(TIMING_FILE_2, read_logfile)
    tt_added = TimingTree.from_json(JSON_ADD_REFERENCE)

    tt1.add(tt2)

    assert_trees_equal(tt1, tt_added)


def test_get_sorted_finish_times():
    tt_json = TimingTree.from_json(JSON_REFERENCE)

    dates = tt_json.get_sorted_finish_times()
    assert dates == [
        datetime(2022, 6, 26, 20, 11, 23)
    ], "sorted finish time does not match reference"
