"""
Test helper utilities

This module provides a collection of utility functions and CLI wrappers for
handling NetCDF data, pandas DataFrames, and various performance and
perturbation operations.
It includes functions to load and process data, perform statistical analysis,
and run CLI commands related to performance, tolerance, perturbation, and
selection of ensemble members.
"""

import logging
import os
import random
import shutil

import numpy as np
import pandas as pd
import xarray as xr
from click.testing import CliRunner

from engine.cdo_table import cdo_table
from engine.performance import performance
from engine.perturb import perturb
from engine.select_members import select_members
from engine.stats import stats
from engine.tolerance import tolerance


def load_netcdf(path):
    return xr.load_dataset(path)


def load_pandas(file, index_col=None, header=None):
    if index_col is None:
        index_col = [0, 1, 2]
    if header is None:
        header = [0, 1]
    return pd.read_csv(file, index_col=index_col, header=header)


def store_timings_as_potential_new_ref(timings, dst):
    # get all files ending with *.json
    for root, _, files in os.walk(timings):
        for file in files:
            if file.endswith(".json"):
                store_as_potential_new_ref(os.path.join(root, file), dst)


def store_as_potential_new_ref(src, dst):
    filename = os.path.basename(src)
    shutil.copyfile(src, os.path.join(dst, filename))


def pandas_error(df_ref, df_cur):
    diff = (df_ref - df_cur).abs()
    err_mask = (diff > 0).any(axis=1)

    return diff[err_mask]


def check_netcdf(data_ref, data_cur):
    diff_keys = set(data_ref.keys()) - set(data_cur.keys())
    same_keys = set(data_ref.keys()).intersection(set(data_cur.keys()))

    err = []
    for key in same_keys:
        diff = np.fabs(data_ref[key] - data_cur[key])
        if np.sum(diff) > 0:
            err.append(key)

    return list(diff_keys), err


def assert_empty_list(lst, msg):
    assert lst == [], f"{msg}:\n{lst}"


def assert_empty_df(df, msg):
    assert len(df.values) == 0, f"{msg}:\n{df}"


def run_performance_cli(timing_regex, timing_database):
    args = [
        "--timing-regex",
        timing_regex,
        "--timing-database",
        timing_database,
    ]
    run_cli(performance, args)


def run_tolerance_cli(
    stats_file_name, tolerance_file_name, member_type=None, member_num=10
):
    args = [
        "--stats-file-name",
        stats_file_name,
        "--tolerance-file-name",
        tolerance_file_name,
        "--member-num",
        str(member_num),
    ]
    if member_type is not None:
        args.append("--member-type")
        args.append(member_type)
    run_cli(tolerance, args)


def generate_ensemble(tmp_path, filename, perturb_amplitude):
    return run_perturb_cli(tmp_path, filename, perturb_amplitude)


def run_perturb_cli(model_input_dir, files, perturb_amplitude, member_num=10):
    perturbed_model_input_dir = f"{model_input_dir}/experiments/{{member_id}}"
    args = [
        "--model-input-dir",
        model_input_dir,
        "--perturbed-model-input-dir",
        perturbed_model_input_dir,
        "--files",
        files,
        "--member-num",
        str(member_num),
        "--member-type",
        "dp",
        "--variable-names",
        "U,V",
        "--perturb-amplitude",
        f"{perturb_amplitude}",
        "--no-copy-all-files",
    ]
    run_cli(perturb, args)

    return perturbed_model_input_dir


def run_stats_cli(
    model_output_dir, stats_file_name, ensemble, perturbed_model_output_dir=None
):
    args = [
        "--model-output-dir",
        model_output_dir,
        "--stats-file-name",
        stats_file_name,
        "--member-type",
        "dp",
        "--file-id",
        "NetCDF",
        "*.nc",
        "--file-specification",
        [
            {
                "NetCDF": {
                    "format": "netcdf",
                    "time_dim": "time",
                    "horizontal_dims": ["lat", "lon"],
                }
            }
        ],
    ]
    args += (
        ["--perturbed-model-output-dir", perturbed_model_output_dir]
        if perturbed_model_output_dir
        else []
    )
    args.append("--ensemble" if ensemble else "--no-ensemble")

    run_cli(stats, args)


def run_cdo_table_cli(model_output_dir, cdo_table_file, perturbed_model_output_dir):
    args = [
        "--model-output-dir",
        model_output_dir,
        "--cdo-table-file",
        cdo_table_file,
        "--member-type",
        "dp",
        "--file-id",
        "NetCDF",
        "*.nc",
        "--perturbed-model-output-dir",
        perturbed_model_output_dir,
        "--file-specification",
        [
            {
                "NetCDF": {
                    "format": "netcdf",
                    "time_dim": "time",
                    "horizontal_dims": ["lat", "lon"],
                }
            }
        ],
    ]
    run_cli(cdo_table, args)


def run_select_members_cli(
    stats_file_name,
    selected_members_file_name,
    tolerance_file_name,
    test_tolerance=False,
    max_member_num=15,
    iterations=50,
    max_factor=50.0,
    log=None,
):
    args = [
        "--stats-file-name",
        stats_file_name,
        "--selected-members-file-name",
        selected_members_file_name,
        "--tolerance-file-name",
        tolerance_file_name,
        "--max-member-num",
        str(max_member_num),
        "--iterations",
        str(iterations),
        "--max-factor",
        str(max_factor),
    ]
    if test_tolerance:
        args.append("--test-tolerance")
    return run_cli(select_members, args, log)


def run_cli(command, args, log=None):
    if log:
        log.set_level(logging.INFO)

    runner = CliRunner()
    result = runner.invoke(command, args)
    if not log:
        catch_error(result)
    else:
        return log.text


def catch_error(result):
    if result.exit_code != 0:
        error_message = "Error executing command:\n" + result.output
        if result.exception:
            error_message += "\nException: " + str(result.exception)
        raise Exception(error_message)


def create_random_stats_file(filename, configurations, seed, perturbation):
    random.seed(seed)
    max_time_dim = max(config["time_dim"] for config in configurations)
    time_header = ",".join(f"{t},{t},{t}" for t in range(max_time_dim))
    header = [
        f"time,,,{time_header}",
        "statistic,," + ",max,mean,min" * max_time_dim,
        "file_ID,variable,height,,,,,,,,,",
    ]

    data = []
    for config in configurations:
        time_dim = config["time_dim"]
        height_dim = config["height_dim"]
        variable = config["variable"]
        file_format = config["file_format"]

        for h in range(height_dim):
            row = f"{file_format},{variable},{h}.0"
            for t in range(time_dim):
                base_mean = round(random.uniform(0, 5), 5)
                mean = base_mean + round(random.uniform(-perturbation, perturbation), 5)
                max_val = mean + round(random.uniform(0, perturbation), 5)
                min_val = mean - round(random.uniform(0, perturbation), 5)
                row += f",{max_val},{mean},{min_val}"
            for _ in range(time_dim, max_time_dim):
                row += ",,,"
            data.append(row)

    with open(filename, "w", encoding="utf-8") as f:
        for line in header:
            f.write(line + "\n")
        for row in data:
            f.write(row + "\n")
