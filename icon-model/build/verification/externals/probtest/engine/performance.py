"""
CLI for Performance Timing

This module provides a command-line interface (CLI) for managing performance
timing data.
It allows users to either append new timing data to an existing database or
overwrite it with new data.
"""

import os
import shutil
import sys
from pathlib import Path

import click

from util.click_util import cli_help
from util.file_system import file_names_from_pattern
from util.icon.extract_timings import read_logfile
from util.log_handler import logger
from util.tree import TREEFILE_TEMPLATE, TimingTree


@click.command()
@click.option(
    "--timing-regex",
    help=cli_help["timing_regex"],
)
@click.option(
    "--timing-database",
    help=cli_help["timing_database"],
)
@click.option("--append-time", help=cli_help["append_time"], type=bool, default=False)
def performance(timing_regex, timing_database, append_time):
    timing_file_name_base, timing_regex = os.path.split(timing_regex)
    if timing_file_name_base == "":
        timing_file_name_base = "."

    timing_file_name, err = file_names_from_pattern(timing_file_name_base, timing_regex)
    if err > 0:
        logger.info("Did not find any files for regex %s", timing_regex)
        sys.exit(1)

    if len(timing_file_name) > 1:
        logger.critical(
            "Found too many files for regex '%s' in '%s':",
            timing_regex,
            timing_file_name_base,
        )
        for tf in timing_file_name:
            logger.critical(tf)
        sys.exit(1)
    else:
        timing_file_name = timing_file_name_base + "/" + timing_file_name[0]
        logger.info("Found timing file %s for regex %s", timing_file_name, timing_regex)

    tt = TimingTree.from_logfile(timing_file_name, read_logfile)

    if TimingTree.input_exists(timing_database):
        if append_time:
            base = TimingTree.from_json(timing_database)
            new_time = tt.meta_data["finish_time"]
            new_revision = tt.meta_data["revision"]
            if (
                new_time in base.meta_data["finish_time"]
                and new_revision in base.meta_data["revision"]
            ):
                logger.info(
                    "this log is already present in the database: %s", timing_database
                )
            else:
                if base.meta_data["n_tables"] != tt.meta_data["n_tables"]:
                    logger.critical(
                        "Cannot add timing file %s\n"
                        + "number of tables does not match. Expected %s, got %s",
                        timing_database,
                        base.meta_data["n_tables"],
                        tt.meta_data["n_tables"],
                    )
                    sys.exit(1)
                for i in range(base.meta_data["n_tables"]):
                    # check if base and tt contain the same nodes (ignores structure)
                    if len(base.root[i].intersection(tt.root[i])) != len(
                        base.root[i].to_list()
                    ):
                        backup = (
                            f"{timing_database}_{base.meta_data['finish_time'][-1]}"
                        )
                        logger.info("tree changed, saving backup as: %s", backup)
                        shutil.copy(
                            TREEFILE_TEMPLATE.format(base=f"{timing_database}_{i}"),
                            TREEFILE_TEMPLATE.format(base=f"{backup}_{i}"),
                        )
                base.add(tt)
                base.json_dump(timing_database)
                logger.info(
                    "appended sample from %s to the database: %s",
                    new_time,
                    timing_database,
                )
        else:  # not append_time
            timing_database_base = os.path.dirname(timing_database)
            tt.json_dump(timing_database)
            logger.info("overwrite old database: %s", timing_database)

    else:
        timing_database_base = os.path.dirname(timing_database)
        if timing_database_base:
            logger.info(
                "creating directory %s, if it doesn't already exist",
                timing_database_base,
            )
        else:
            logger.info("creating database locally")
        Path(timing_database_base).mkdir(parents=True, exist_ok=True)
        tt.json_dump(timing_database)
        logger.info("created new database: %s", timing_database)
