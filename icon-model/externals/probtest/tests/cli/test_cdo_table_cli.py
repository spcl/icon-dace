"""
This module contains pytest-based tests for the CDO table CLI functionality in the
PROBTEST suite. It verifies the correctness of the CDO table generation by comparing
the output with reference data and checks the behavior under different perturbation
amplitudes.
"""

import os

import pytest

from tests.helpers import (
    assert_empty_df,
    load_pandas,
    pandas_error,
    run_cdo_table_cli,
    store_as_potential_new_ref,
)


def test_cdo_table_cli(tmp_dir, ensemble, df_ref_cdo_table, new_ref):
    cdo_table_file = os.path.join(tmp_dir, "cdo_table.csv")
    run_cdo_table_cli(tmp_dir, cdo_table_file, ensemble)
    df_test = load_pandas(cdo_table_file, index_col=[0, 1])
    err = pandas_error(df_ref_cdo_table, df_test)

    store_as_potential_new_ref(cdo_table_file, new_ref)

    assert_empty_df(err, "CDO table datasets are not equal!")


@pytest.mark.xfail(
    reason="perturb amplitude 10e-14 < 10e-12",
    strict=True,
)
def test_cdo_table_cli_with_too_small_perturb_amplitude(
    tmp_dir, too_small_ensemble, df_ref_cdo_table
):
    cdo_table_file = os.path.join(tmp_dir, "cdo_table.csv")
    run_cdo_table_cli(tmp_dir, cdo_table_file, too_small_ensemble)
    df_test = load_pandas(cdo_table_file, index_col=[0, 1])
    err = pandas_error(df_ref_cdo_table, df_test)

    assert_empty_df(err, "CDO table datasets are not equal!")
