import sys
from collections.abc import Iterable

import numpy as np
import pandas as pd
import xarray

from util.constants import compute_statistics
from util.log_handler import logger
from util.utils import numbers
from util.xarray_ops import statistics_over_horizontal_dim

# Definition of available model output data parsers
#
# All available parsers are summarized in the dict `model_output_parser` at the
# end of this module.
#
# Each parser expects two arguments:
# file_id: str
#    identifier of the file type
# filename: str
#    full file name
# specification: dict
#    Dictionary with type specific settings
#
# Each parser returns a list of Pandas DataFrame:
#    index = pd.MultiIndex.from_product(
#        [[file_id], [varname], height], names=("file_ID", "variable", "height")
#    )
#    columns = pd.MultiIndex.from_product(
#        [time, compute_statistics], names=("time", "statistic")
#    )
#
#    return [pd.DataFrame(matrix, index=index, columns=columns)]


def parse_netcdf(file_id, filename, specification):
    logger.debug("parse NetCDF file {}".format(filename))
    time_dim = specification["time_dim"]
    horizontal_dims = specification["horizontal_dims"]
    fill_value_key = specification.get("fill_value_key", None)
    ds = xarray.open_dataset(filename, decode_cf=False)

    var_tmp = __get_variables(ds, time_dim, horizontal_dims)

    var_dfs = []

    for v in var_tmp:
        sub_df = dataframe_from_ncfile(
            file_id=file_id,
            filename=filename,
            varname=v,
            time_dim=time_dim,
            horizontal_dims=horizontal_dims,
            xarray_ds=ds,
            fill_value_key=fill_value_key,
        )
        var_dfs.append(sub_df)

    ds.close()
    return var_dfs


def __get_variables(data, time_dim, horizontal_dims):
    # return a list of variable names from the dataset data that have a time dimension
    # and horizontal dimension or in case there is no time dimension just the variables
    #  with horizontal dimension
    all_variables = data.variables.keys()
    all_variables = [
        v for v in all_variables if np.issubdtype(data.variables[v].dtype, np.number)
    ]

    variables = []

    if isinstance(horizontal_dims, Iterable):
        horizontal_dims_unpacked = [x for hd in horizontal_dims for x in hd.split(":")]

    if time_dim is not None:
        if isinstance(horizontal_dims, Iterable):
            for v in all_variables:
                if (
                    any(hd in horizontal_dims_unpacked for hd in data.variables[v].dims)
                    and time_dim in data.variables[v].dims
                ):
                    variables.append(v)
        else:
            variables = [
                v
                for v in all_variables
                if time_dim in data.variables[v].dims and v != time_dim
            ]
    else:
        if isinstance(horizontal_dims, Iterable):
            for v in all_variables:
                if any(hd in horizontal_dims_unpacked for hd in data.variables[v].dims):
                    variables.append(v)
        else:
            logger.error(
                "horizontal_dims have to be specified when no time_dim is given"
            )
            sys.exit(1)

    return variables


def dataframe_from_ncfile(
    file_id, filename, varname, time_dim, horizontal_dims, xarray_ds, fill_value_key
):
    statistics = statistics_over_horizontal_dim(
        xarray_ds[varname],
        horizontal_dims,
        compute_statistics,
        fill_value_key,
    )

    first_stat = statistics[0]
    if len(first_stat.dims) == 2:
        height_name = (
            first_stat.dims[0] if first_stat.dims[0] != time_dim else first_stat.dims[1]
        )  # might be 'height', 'height_2', 'alt', 'plev', ...
        height = xarray_ds[height_name].values
        matrix = np.empty(
            (first_stat.shape[0] * len(statistics), first_stat.shape[1])
        )  # shape: (mean_max_min of each height, time)

        # weave mean max min into time dimension
        for i, stat in enumerate(statistics):
            matrix[i :: len(statistics), :] = stat.values
    elif len(first_stat.dims) == 1:
        if first_stat.dims[0] == time_dim:
            height = np.array([-1])
            # matrix needs to have 2 dimensions for DataFrame constructor
            matrix = np.empty((first_stat.size * len(statistics), 1))
            # weave mean max min into time dimension
            for i, stat in enumerate(statistics):
                matrix[i :: len(statistics), 0] = stat.values
        else:
            height = xarray_ds[first_stat.dims[0]].values
            # matrix needs to have 2 dimensions for DataFrame constructor
            matrix = np.empty((len(statistics), first_stat.size))
            # weave mean max min into time dimension
            for i, stat in enumerate(statistics):
                matrix[i, :] = stat.values
    elif len(first_stat.dims) == 0:
        height = np.array([-1])
        # matrix needs to have 2 dimensions for DataFrame constructor
        matrix = np.empty((len(statistics), 1))
        # weave mean max min into time dimension
        for i, stat in enumerate(statistics):
            matrix[i :: len(statistics), 0] = stat.values
    else:
        logger.error(
            (
                "Unknown number of dimension for first_stat of variable '{}'. "
                + "Dims: {}"
            ).format(varname, first_stat.dims)
        )
        sys.exit(1)

    if time_dim is not None:
        time = xarray_ds[time_dim].values
    else:
        # Derive a pseudo time from filename. This is required to process multiple
        # files of the same file type if the file type has not time dimension.
        time = [numbers(filename)]

    index = pd.MultiIndex.from_product(
        [[file_id], [varname], height], names=("file_ID", "variable", "height")
    )
    columns = pd.MultiIndex.from_product(
        [time, compute_statistics], names=("time", "statistic")
    )

    return pd.DataFrame(matrix.T, index=index, columns=columns)


def parse_csv(file_id, filename, specification):
    """
    file_id: str
        identifier of the file type
    filename: str
        full file name
    specification: dict(parser_args, time_dim, horizontal_dims)
        parser_args: dict
            passed directly to pandas.read_csv

    The (first) index of the read csv (i.e. usually the rows) is expected to
    represent the time dimension.
    Each column contains a individual variable
    """
    logger.debug("parse CSV file {}".format(filename))

    csv = pd.read_csv(filename, **specification["parser_args"])

    # transpose data such that time is along columns
    csv = csv.transpose()

    if csv.columns.nlevels == 1:
        # use the values in csv as dummy for each of the expected compute_statistics
        matrix = np.array(csv).repeat(len(compute_statistics), 1)
        columns = csv.columns
        # regular Index -> "height" information does not apply
        height = [-1]
    else:
        n_time = len(csv.columns.levels[0])
        # use "height" for independent levels of a variable
        height = np.arange(csv.columns.size / n_time, dtype=int)

        # convert to proper multidimensional array
        array = np.array(csv).reshape(csv.index.size, n_time, -1)
        # transpose such that time is in last dimension
        array = array.transpose(0, 2, 1)  # index, "height", time
        # collapse index and "height" dimensions
        array = array.reshape(-1, n_time)
        matrix = array.repeat(len(compute_statistics), 1)

        columns = csv.columns.levels[0]

    index = pd.MultiIndex.from_product(
        [[file_id], csv.index, height], names=("file_ID", "variable", "height")
    )
    columns = pd.MultiIndex.from_product(
        [columns, compute_statistics], names=("time", "statistic")
    )
    return [pd.DataFrame(matrix, index=index, columns=columns)]


model_output_parser = {  # global lookup dict
    "netcdf": parse_netcdf,
    "csv": parse_csv,
}
