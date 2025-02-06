"""
This module contains tests for verifying the correctness of the tolerance CLI
command.
The tests compare the output of the tolerance CLI with reference tolerance data
and check for discrepancies.
"""

import os

from tests.helpers import (
    assert_empty_df,
    load_pandas,
    pandas_error,
    run_tolerance_cli,
    store_as_potential_new_ref,
)


def test_tolerance_cli(ref_data, df_ref_tolerance, tmp_dir, new_ref):
    stats_file_name = os.path.join(ref_data, "stats_{member_id}.csv")
    tolerance_file_name = os.path.join(tmp_dir, "tolerance.csv")
    run_tolerance_cli(stats_file_name, tolerance_file_name, member_type="dp")
    df_test = load_pandas(tolerance_file_name, index_col=[0, 1])
    err = pandas_error(df_ref_tolerance, df_test)

    store_as_potential_new_ref(tolerance_file_name, new_ref)

    assert_empty_df(err, "Tolerance datasets are not equal!")
