"""
This module contains test cases for verifying the functionality of the
`stats_cli` command-line interface (CLI).
"""

import os

import pytest

from tests.helpers import (
    assert_empty_df,
    load_pandas,
    pandas_error,
    run_stats_cli,
    store_as_potential_new_ref,
)


def test_stats_cli_no_ensemble(new_ref, nc_with_t_u_v, df_ref_stats):
    tmp_path = os.path.dirname(nc_with_t_u_v)
    stats_file = os.path.join(tmp_path, "stats_{member_id}.csv")
    run_stats_cli(tmp_path, stats_file.format(member_id="ref"), ensemble=False)
    df_test = load_pandas(
        stats_file.format(member_id="ref"), index_col=[0, 1, 2], header=[0, 1]
    )
    err = pandas_error(df_ref_stats, df_test)

    store_as_potential_new_ref(stats_file.format(member_id="ref"), new_ref)

    assert_empty_df(err, "Stats datasets are not equal!")


@pytest.mark.xfail(
    reason="perturb amplitude 10e-14 < 10e-12",
    strict=True,
)
@pytest.mark.parametrize("member", range(1, 11))
def test_stats_cli_ensemble_with_too_small_perturb_amplitude_for_member(
    tmp_dir, too_small_ensemble, df_ref_ensemble_stats, member
):
    stats_file = os.path.join(tmp_dir, "stats_{member_id}.csv")
    run_stats_cli(
        tmp_dir,
        stats_file,
        ensemble=True,
        perturbed_model_output_dir=too_small_ensemble,
    )
    df_test = load_pandas(
        stats_file.format(member_id="dp_" + str(member)),
        index_col=[0, 1, 2],
        header=[0, 1],
    )
    err = pandas_error(df_ref_ensemble_stats[member], df_test)

    assert_empty_df(err, "Stats datasets are not equal!")


@pytest.mark.parametrize("member", range(1, 11))
def test_stats_cli_ensemble_for_member(
    tmp_dir, ensemble, df_ref_ensemble_stats, new_ref, member
):
    stats_file = os.path.join(tmp_dir, "stats_{member_id}.csv")
    run_stats_cli(
        tmp_dir, stats_file, ensemble=True, perturbed_model_output_dir=ensemble
    )
    df_test = load_pandas(
        stats_file.format(member_id="dp_" + str(member)),
        index_col=[0, 1, 2],
        header=[0, 1],
    )
    err = pandas_error(df_ref_ensemble_stats[member], df_test)

    store_as_potential_new_ref(
        stats_file.format(member_id="dp_" + str(member)), new_ref
    )

    assert_empty_df(err, "Stats datasets are not equal!")
