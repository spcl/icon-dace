"""
This module contains test cases for validating the functionality of the
`performance_cli` command-line interface (CLI).
The tests are focused on checking the correctness of timing and performance data
generated from the CLI, including comparing the results with reference data.
"""

import os

from tests.helpers import (
    assert_empty_df,
    assert_empty_list,
    pandas_error,
    run_performance_cli,
    store_timings_as_potential_new_ref,
)
from util.tree import TimingTree


def test_performance_cli_timing(new_ref, timing_logfile, tmp_dir, df_ref_performance):
    timing_database = os.path.join(tmp_dir, "test")
    run_performance_cli(timing_logfile, timing_database)
    df_test = TimingTree.from_json(timing_database)

    store_timings_as_potential_new_ref(tmp_dir, new_ref)

    for perf_data, test_data in zip(df_ref_performance.data, df_test.data):
        err = pandas_error(perf_data, test_data)

        assert_empty_df(err, "Performance datasets are not equal!")


def test_performance_cli_tree(timing_logfile, tmp_dir, df_ref_performance):
    tree_database = os.path.join(tmp_dir, "test")
    run_performance_cli(timing_logfile, tree_database)
    df_test = TimingTree.from_json(tree_database)

    for ref_data, test_data in zip(df_ref_performance.root, df_test.root):
        ref_nodes = set(ref_data.to_ancestry_name_list())
        cur_nodes = set(test_data.to_ancestry_name_list())

        diff_nodes = list(ref_nodes - cur_nodes)

        assert_empty_list(
            diff_nodes,
            "The test tree does not contain the following nodes:",
        )
