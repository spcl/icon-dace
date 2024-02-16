import sys

import numpy as np
import pandas as pd

from util.constants import CHECK_THRESHOLD
from util.file_system import file_names_from_pattern
from util.log_handler import logger
from util.model_output_parser import model_output_parser

pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.float_format", "{:,.2e}".format)


def force_monotonic(dataframe):
    stats = list(dataframe.columns.levels[1])
    for s in stats:
        dataframe.loc[:, (slice(None), s)] = dataframe.loc[:, (slice(None), s)].cummax(
            axis=1
        )


def compute_rel_diff_dataframe(df1, df2):
    average = (df1 + df2) / 2
    out = (df1 - df2) / average
    out = out.abs()
    # put 0 if both numbers are very small
    zeros = np.logical_and(df1.abs() < CHECK_THRESHOLD, df2.abs() < CHECK_THRESHOLD)
    out[zeros] = 0.0
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
            "No parser defined for format `{}` of file `{}`.".format(
                specification["format"], file_name
            )
        )
        sys.exit(1)

    var_dfs = file_parser(label, file_name, specification)

    if len(var_dfs) == 0:
        logger.error("Could not find any variables in `{}`".format(file_name))
        logger.error("Wrong file format or specification? Fid: `{}` ".format(label))
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
                "Can not find any files for file_pattern {}. Continue.".format(
                    file_pattern
                )
            )
            continue

        try:
            specification = file_specification[file_type]
        except KeyError:
            logger.error(
                "No parser defined for format `{}` of file_pattern `{}`.".format(
                    file_type, file_pattern
                )
            )
            sys.exit(1)

        file_dfs = []
        for f in input_files:
            var_df = read_input_file(
                label="{}:{}".format(file_type, file_pattern),
                file_name="{}/{}".format(input_dir, f),
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
