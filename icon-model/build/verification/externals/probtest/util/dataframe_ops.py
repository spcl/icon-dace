"""
This module provides functions for statistical analysis and comparison of
datasets, primarily for performance testing and validation.
It includes utilities to handle data reading, processing, and comparison against
reference datasets with specified tolerances.
"""

import sys
import warnings

import numpy as np
import pandas as pd

from util.constants import CHECK_THRESHOLD, compute_statistics
from util.file_system import file_names_from_pattern
from util.log_handler import logger
from util.model_output_parser import model_output_parser

pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.float_format", lambda x: f"{x:,.2e}")


def force_monotonic(dataframe):
    stats = list(dataframe.columns.levels[1])
    for s in stats:
        dataframe.loc[:, (slice(None), s)] = dataframe.loc[:, (slice(None), s)].cummax(
            axis=1
        )


def compute_rel_diff_dataframe(df1, df2):
    """This implementation is similar to the numpy.isclose function:
    (absolute(a - b) <= (atol + rtol * absolute(b)) ),
    assuming atol==rtol and moving the right hand side to the left."""
    out = (df1 - df2) / (1.0 + df1.abs())
    out = out.abs()
    return out


def compute_div_dataframe(df1, df2):
    # avoid division by 0 and put nan instead
    out = df1 / df2.replace({0: np.nan})
    # put 0 if numerator is 0 as well
    out[df1 == 0] = 0
    return out


def parse_probtest_csv(path, index_col):
    df = pd.read_csv(path, index_col=index_col, header=[0, 1])

    times = df.columns.levels[0].astype(int)
    df.columns = df.columns.set_levels(times, level=0)

    # the dataframe's time column will be read as string,
    # thus ordered like "0", "1", "10", "11", .. "2", ...
    new_cols = pd.MultiIndex.from_product(
        [sorted(df.columns.levels[0]), df.columns.levels[1]],
        names=df.columns.names,
    )

    return pd.DataFrame(df, columns=new_cols)


def read_input_file(label, file_name, specification):
    """Read input file file_name using the specification."""
    try:
        file_parser = model_output_parser[specification["format"].lower()]
    except KeyError:
        logger.error(
            "No parser defined for format `%s` of file `%s`.",
            specification["format"],
            file_name,
        )
        sys.exit(1)

    var_dfs = file_parser(label, file_name, specification)

    if len(var_dfs) == 0:
        logger.error("Could not find any variables in `%s`", file_name)
        logger.error("Wrong file format or specification? Fid: `%s` ", label)
        sys.exit(1)

    # different variables in a file have same timestamps:
    # concatenate along variable axis
    return pd.concat(var_dfs, axis=0)


def df_from_file_ids(file_id, input_dir, file_specification):
    """
    file_id: [[file_type, file_pattern], [file_type, file_pattern], ...]
        List of 2-tuples. The 2-tuple combines two strings. The first sets the
        file_type and must be a key in file_specification. The second string
        is a file name pattern. file_pattern is extended to real file names using
        glob. file_pattern may contain simple shell-style wildcards such as "*".
    file_specification: {file_type: specification, ...}
        file_type: str
            Name of the file type specification
        specification: dict(format: str, **kwargs)
            dictionary that specifies the file format defined by the key.
            format: str
                The format defines the corresponding parser function.
            **kwargs:
                Other keyword arguments can be used by the corresponding parser.
                Examples:
                time_dim: str
                    The time dimension in the model output files
                horizontal_dims: str
                    List possible horizontal dimensions. If multiple horizontal
                    dimensions are used their names must be combined in one string
                    separated by ":".
    """

    # Collect data frames for each combination of file id (fid) and
    # specification (spec). Frames for the same fid and spec represent
    # different timestamps and have to be concatenated along time-axis (axis=1).
    # Time-concatenated frames from different ids and/or specifications will be
    # concatenated along variable-axis (axis=0).
    fid_dfs = []
    for file_type, file_pattern in file_id:
        input_files, err = file_names_from_pattern(input_dir, file_pattern)
        if err > 0:
            logger.info(
                "Can not find any files for file_pattern %s. Continue.", file_pattern
            )
            continue

        try:
            specification = file_specification[file_type]
        except KeyError:
            logger.error(
                "No parser defined for format `%s` of file_pattern `%s`.",
                file_type,
                file_pattern,
            )
            sys.exit(1)

        file_dfs = []
        for f in input_files:
            var_df = read_input_file(
                label=f"{file_type}:{file_pattern}",
                file_name=f"{input_dir}/{f}",
                specification=specification,
            )
            file_dfs.append(var_df)

        # same file IDs and file type specification will have same variables but
        # with different timestamps: concatenate along time axis
        fid_dfs.append(pd.concat(file_dfs, axis=1))

    if len(fid_dfs) == 0:
        logger.error("Could not find any file.")
        sys.exit(2)

    fid_dfs = unify_time_index(fid_dfs)
    # different file IDs will have different variables but with same timestamps:
    # concatenate along variable axis
    df = pd.concat(fid_dfs, axis=0)

    return df


def unify_time_index(fid_dfs):
    """
    Unify the column index of all data frames by replacing the time with a
    continuous range from 0 to N-1

    fid_dfs: Iterable(pd.DataFrame)
    """
    fid_dfs_out = []
    for df in fid_dfs:
        # Find out number of time steps in the column MultiIndex.
        # Is there an easier way to do so without assuming the order of indices?
        time_multiindex_index = df.columns.names.index("time")
        ntime = len(df.columns.levels[time_multiindex_index])

        unique_times = list(df.columns.levels[time_multiindex_index])
        unique_times.sort()
        df = df.reindex(columns=unique_times, level=time_multiindex_index)

        df.columns = df.columns.set_levels(range(ntime), level="time")

        fid_dfs_out.append(df)

    return fid_dfs_out


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
            f" reference case and therefore not tested: {', '.join(missing_in_ref)}"
        )
        warnings.warn(warning_msg, UserWarning)
    if missing_in_cur:
        warning_msg = (
            "WARNING: The following variables are in the reference case but not in the"
            f" test case and therefore not tested: {', '.join(missing_in_cur)}"
        )
        warnings.warn(warning_msg, UserWarning)

    # Remove rows without intersection
    df_ref = df_ref[~df_ref.index.isin(non_common_vars)]
    df_cur = df_cur[~df_cur.index.isin(non_common_vars)]

    # Make sure they have the same number of time steps
    if len(df_ref.columns) > len(df_cur.columns):
        logger.info(
            "WARNING: The reference includes more timesteps than the test case. "
            "Only the first %s time step(s) are tested.\n",
            len(df_cur.columns) // len(compute_statistics),
        )
        df_ref = df_ref.iloc[:, : len(df_cur.columns)]
    elif len(df_ref.columns) < len(df_cur.columns):
        logger.info(
            "WARNING: The reference includes less timesteps than the test case. "
            "Only the first %s time step(s) are tested.\n",
            len(df_ref.columns) // len(compute_statistics),
        )
        df_cur = df_cur.iloc[:, : len(df_ref.columns)]
    return skip_test, df_ref, df_cur


def check_variable(diff_df, df_tol):
    out = diff_df - df_tol

    selector = (out > CHECK_THRESHOLD).any(axis=1)

    return len(out[selector].index) == 0, diff_df[selector], df_tol[selector]


def parse_check(tolerance_file_name, input_file_ref, input_file_cur, factor):
    """
    Parses all necessary data to perform a check from tolerance, reference and
    input files, applying scaling factor to tolerance values.

    Args:
        tolerance_file_name (str): Path to the CSV file containing tolerance values.
        input_file_ref (str): Path to the reference input CSV file.
        input_file_cur (str): Path to the current input CSV file.
        factor (float): Scaling factor to be applied to the tolerance values.

    Returns:
        tuple: A tuple containing three DataFrames:
            - df_tol (pandas.DataFrame): The tolerance DataFrame with values
                                         scaled by the provided factor.
            - df_ref (pandas.DataFrame): The reference DataFrame parsed from the
                                         reference input file.
            - df_cur (pandas.DataFrame): The current DataFrame parsed from the
                                         current input file.
    """
    df_tol = parse_probtest_csv(tolerance_file_name, index_col=[0, 1])

    logger.info("applying a factor of %s to the spread", factor)
    df_tol *= factor

    df_ref = parse_probtest_csv(input_file_ref, index_col=[0, 1, 2])
    df_cur = parse_probtest_csv(input_file_cur, index_col=[0, 1, 2])

    logger.info(
        "checking %s against %s using tolerances from %s",
        input_file_cur,
        input_file_ref,
        tolerance_file_name,
    )

    return df_tol, df_ref, df_cur


def test_stats_file_with_tolerances(
    tolerance_file_name, input_file_ref, input_file_cur, factor
):

    df_tol, df_ref, df_cur = parse_check(
        tolerance_file_name, input_file_ref, input_file_cur, factor
    )

    # check if variables are available in reference file
    skip_test, df_ref, df_cur = check_intersection(df_ref, df_cur)
    if skip_test:  # No intersection
        logger.error("RESULT: check FAILED")
        sys.exit(1)

    # compute relative difference
    diff_df = compute_rel_diff_dataframe(df_ref, df_cur)
    # take maximum over height
    diff_df = diff_df.groupby(["file_ID", "variable"]).max()

    out, err, tol = check_variable(diff_df, df_tol)

    return out, err, tol
