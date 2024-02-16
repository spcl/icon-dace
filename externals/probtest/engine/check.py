import warnings

import click

from util.click_util import cli_help
from util.constants import CHECK_THRESHOLD, compute_statistics
from util.dataframe_ops import (
    compute_div_dataframe,
    compute_rel_diff_dataframe,
    parse_probtest_csv,
)
from util.log_handler import logger


def check_intersection(df_ref, df_cur):
    # Check if variable names in reference and test case have any intersection
    # Check if numbers of time steps agree
    skip_test = 0
    if not set(df_ref.index.intersection(df_cur.index)):
        logger.info(
            "WARNING: No intersection between variables in input and reference file."
        )
        skip_test = 1
        return skip_test, df_ref, df_cur

    # Check if there are variable missing in the reference or test case
    non_common_vars = list(set(df_ref.index) ^ set(df_cur.index))
    missing_in_ref = []
    missing_in_cur = []
    for var in non_common_vars:
        if var in df_cur.index:
            missing_in_ref.append(var[1])
        else:
            missing_in_cur.append(var[1])
    # Remove multiple entries of a variable due to different altitude levels
    missing_in_ref = list(set(missing_in_ref))
    missing_in_cur = list(set(missing_in_cur))
    if missing_in_ref:
        warning_msg = (
            "WARNING: The following variables are in the test case but not in the"
            " reference case and therefore not tested: {}".format(
                ", ".join(missing_in_ref)
            )
        )
        warnings.warn(warning_msg, UserWarning)
    if missing_in_cur:
        warning_msg = (
            "WARNING: The following variables are in the reference case but not in the"
            " test case and therefore not tested: {}".format(", ".join(missing_in_cur))
        )
        warnings.warn(warning_msg, UserWarning)

    # Remove rows without intersection
    df_ref = df_ref[~df_ref.index.isin(non_common_vars)]
    df_cur = df_cur[~df_cur.index.isin(non_common_vars)]

    # Make sure they have the same number of time steps
    if len(df_ref.columns) > len(df_cur.columns):
        logger.info(
            "WARNING: The reference includes more timesteps than the test case. "
            "Only the first {} time step(s) are tested.\n".format(
                len(df_cur.columns) // len(compute_statistics)
            )
        )
        df_ref = df_ref.iloc[:, : len(df_cur.columns)]
    elif len(df_ref.columns) < len(df_cur.columns):
        logger.info(
            "WARNING: The reference includes less timesteps than the test case. "
            "Only the first {} time step(s) are tested.\n".format(
                len(df_ref.columns) // len(compute_statistics)
            )
        )
        df_cur = df_cur.iloc[:, : len(df_ref.columns)]
    return skip_test, df_ref, df_cur


def check_variable(diff_df, df_tol):
    out = diff_df - df_tol

    selector = (out > CHECK_THRESHOLD).any(axis=1)

    return len(out[selector].index) == 0, diff_df[selector], df_tol[selector]


@click.command()
@click.option(
    "--input-file-ref",
    help=cli_help["input_file_ref"],
)
@click.option(
    "--input-file-cur",
    help=cli_help["input_file_cur"],
)
@click.option(
    "--tolerance-file-name",
    help=cli_help["tolerance_file_name"],
)
@click.option(
    "--factor",
    type=float,
    help=cli_help["factor"],
)
def check(input_file_ref, input_file_cur, tolerance_file_name, factor):
    df_tol = parse_probtest_csv(tolerance_file_name, index_col=[0, 1])

    logger.info("applying a factor of {} to the spread".format(factor))
    df_tol *= factor

    df_ref = parse_probtest_csv(input_file_ref, index_col=[0, 1, 2])
    df_cur = parse_probtest_csv(input_file_cur, index_col=[0, 1, 2])

    logger.info(
        "checking {} against {} using tolerances from {}".format(
            input_file_cur, input_file_ref, tolerance_file_name
        )
    )

    # check if variables are available in reference file
    skip_test, df_ref, df_cur = check_intersection(df_ref, df_cur)
    if skip_test:  # No intersection
        logger.info("RESULT: check FAILED")
        exit(1)

    # compute relative difference
    diff_df = compute_rel_diff_dataframe(df_ref, df_cur)
    # take maximum over height
    diff_df = diff_df.groupby(["file_ID", "variable"]).max()

    out, err, tol = check_variable(diff_df, df_tol)

    div = compute_div_dataframe(err, tol)

    if out:
        logger.info("RESULT: check PASSED!")
    else:
        logger.info("RESULT: check FAILED")
        logger.info("Differences")
        print(err)
        logger.info("\nTolerance")
        print(tol)
        logger.info("\nError relative to tolerance")
        print(div)

    exit(0 if out else 1)
