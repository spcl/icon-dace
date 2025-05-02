"""
CLI for CDO Table Generation

This module computes and generates a CDO table by comparing model output data
against perturbed model data.
"""

import tempfile
from pathlib import Path

import click
import numpy as np
import pandas as pd
import xarray as xr

from util import model_output_parser
from util.click_util import CommaSeperatedInts, cli_help
from util.constants import cdo_bins
from util.dataframe_ops import df_from_file_ids
from util.file_system import file_names_from_pattern
from util.log_handler import logger


def compute_rel_diff(var1, var2):
    rel_diff = np.zeros_like(var1)

    mask_0 = np.logical_and(np.abs(var1) < 1e-15, np.abs(var2) < 1e-15)
    mask_inf = np.logical_and(np.abs(var1) < 1e-15, np.abs(var2) > 1e-15)
    mask_valid = np.logical_and(np.abs(var1) > 1e-15, np.abs(var2) > 1e-15)

    rel_diff[mask_0] = 0
    rel_diff[mask_inf] = 999
    rel_diff[mask_valid] = np.abs(
        (var1[mask_valid] - var2[mask_valid]) / var1[mask_valid]
    )

    return rel_diff


def rel_diff_stats(
    file_id,
    filename,
    varname,
    time_dim,
    horizontal_dims,
    xarray_ds,
    fill_value_key,
):  # pylint: disable=unused-argument
    dims = xarray_ds[varname].dims
    dataarray = xarray_ds[varname]
    time = xarray_ds[time_dim].values

    # compute maximum over all axis except time
    dims_without_time = [d for d in dims if d != time_dim]
    amax = dataarray.max(dim=dims_without_time)

    # allocate empty matrix
    ncol = 9
    matrix = np.empty((amax.size * ncol))

    # compute histogram of relative differences
    for i in range(amax.size):
        data = dataarray[{time_dim: i}].values.flatten()
        if fill_value_key and fill_value_key in dataarray.attrs:
            mask = data != dataarray.attrs[fill_value_key]
            data = data[mask]

        hist, _ = np.histogram(data, bins=[0] + cdo_bins)
        matrix[i * ncol] = amax[i]
        matrix[i * ncol + 1 : i * ncol + ncol] = hist

    index = pd.MultiIndex.from_product(
        [[file_id], [varname]], names=("file_ID", "variable")
    )
    columns = pd.MultiIndex.from_product(
        [time, ["rel_diff"] + cdo_bins], names=["time", "statistic"]
    )

    return pd.DataFrame(matrix[np.newaxis, :], index=index, columns=columns)


@click.command()
@click.option(
    "--model-output-dir",
    help=cli_help["model_output_dir"],
)
@click.option(
    "--file-id",
    nargs=2,
    type=str,
    multiple=True,
    metavar="FILE_TYPE FILE_PATTERN",
    help=cli_help["file_id"],
)
@click.option(
    "--member-num",
    type=CommaSeperatedInts(),
    default="10",
    help=cli_help["member_num"],
)
@click.option(
    "--member-type",
    type=str,
    default="",
    help=cli_help["member_type"],
)
@click.option(
    "--perturbed-model-output-dir",
    help=cli_help["perturbed_model_output_dir"],
)
@click.option(
    "--cdo-table-file",
    help=cli_help["cdo_table_file"],
)
@click.option(
    "--file-specification",
    type=list,
    help=cli_help["file_specification"],
)
def cdo_table(
    model_output_dir,
    file_id,
    member_num,
    member_type,
    perturbed_model_output_dir,
    cdo_table_file,
    file_specification,
):
    # TODO: A single perturbed run provides enough data to make proper statistics.
    #       refactor cdo_table interface to reflect that
    if len(member_num) == 1:
        member_num = list(range(1, member_num[0] + 1))
    if member_type:
        member_id = member_type + "_" + str(member_num[0])
    else:
        member_id = str(member_num[0])

    file_specification = file_specification[0]  # can't store dicts as defaults in click
    assert isinstance(file_specification, dict), "must be dict"

    # save original method and restore at the end of this module
    dataframe_from_ncfile_orig = model_output_parser.dataframe_from_ncfile
    # modify netcdf parse method:
    model_output_parser.dataframe_from_ncfile = rel_diff_stats

    # step 1: compute rel-diff netcdf files
    with tempfile.TemporaryDirectory() as tmpdir:
        for _, file_pattern in file_id:
            ref_files, err = file_names_from_pattern(model_output_dir, file_pattern)
            if err > 0:
                logger.info(
                    "did not find any files for pattern %s. Continue.", file_pattern
                )
                continue
            ref_files.sort()
            perturb_files, err = file_names_from_pattern(
                perturbed_model_output_dir.format(member_id=member_id), file_pattern
            )
            if err > 0:
                logger.info(
                    "did not find any files for pattern %s. Continue.", file_pattern
                )
                continue
            perturb_files.sort()

            for rf, pf in zip(ref_files, perturb_files):
                if not rf.endswith(".nc") or not pf.endswith(".nc"):
                    continue
                ref_data = xr.open_dataset(f"{model_output_dir}/{rf}")
                perturb_data = xr.open_dataset(
                    f"{perturbed_model_output_dir.format(member_id=member_id)}/{pf}"
                )
                diff_data = ref_data.copy()
                varnames = [
                    v
                    for v in list(ref_data.keys())
                    if "time" in ref_data.variables.get(v).dims
                ]

                for v in varnames:
                    diff_data.variables.get(v).values = compute_rel_diff(
                        ref_data.variables.get(v).values,
                        perturb_data.variables.get(v).values,
                    )

                diff_data.to_netcdf(f"{tmpdir}/{rf}")
                ref_data.close()
                perturb_data.close()
                diff_data.close()

        # step 2: generate dataframe from precomputed relative differences
        df = df_from_file_ids(file_id, tmpdir, file_specification)

        # normalize histogram component of DataFrame
        times = np.array(df.columns.levels[0], dtype=int)
        for t in times:
            df.loc[:, (t, cdo_bins)] = df.loc[:, (t, cdo_bins)].div(
                df.loc[:, (t, cdo_bins)].sum(axis=1), axis=0
            )

        logger.info("writing cdo table to %s.", cdo_table_file)

        Path(cdo_table_file).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(cdo_table_file)
        model_output_parser.dataframe_from_ncfile = dataframe_from_ncfile_orig
