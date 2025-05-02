"""
CLI for performance checks

This module provides a command-line interface (CLI) tool for evaluating the
performance of a current experiment by comparing its timing data with a
reference.
It assesses whether the current experiment's runtime is within acceptable limits
based on specified parameters and provides feedback on whether the performance
check has passed or failed.
"""

import sys

import click

from util.click_util import cli_help
from util.log_handler import logger
from util.tree import TimingTree


@click.command()
@click.option(
    "--timing-current",
    help=cli_help["timing_current"],
)
@click.option(
    "--timing-reference",
    help=cli_help["timing_reference"],
)
@click.option(
    "--measurement-uncertainty",
    help=cli_help["measurement_uncertainty"],
    type=float,
    default=2,
)
@click.option(
    "--tolerance-factor", help=cli_help["tolerance_factor"], type=float, default=1.1
)
@click.option(
    "--new-reference-threshold",
    help=cli_help["new_reference_threshold"],
    type=float,
    default=0.95,
)
@click.option("--i-table", type=int, help=cli_help["i_table"], default=-1)
def performance_check(
    timing_current,
    timing_reference,
    i_table,
    measurement_uncertainty,
    tolerance_factor,
    new_reference_threshold,
):
    ttcur = TimingTree.from_json(timing_current)
    ttref = TimingTree.from_json(timing_reference)

    total_time_cur = (
        ttcur.data[i_table].loc[("total", slice(None)), "total max (s)"].values[0]
    )
    total_time_ref = (
        ttref.data[i_table].loc[("total", slice(None)), "total max (s)"].values[0]
    )

    if measurement_uncertainty < 0:
        logger.error("measurement_uncertainty needs to be positive")
    if tolerance_factor < 1:
        logger.error("tolerance_factor needs to be greater than 1")
    if new_reference_threshold < 0 or new_reference_threshold > 1:
        logger.error("new_reference_threshold needs to be between 0 and 1")

    allowed_time = (total_time_ref + measurement_uncertainty) * tolerance_factor

    logger.info("Current runtime")
    logger.info(total_time_cur)
    logger.info("Allowed runtime")
    logger.info(allowed_time)
    logger.info("Reference runtime")
    logger.info(total_time_ref)

    if total_time_cur <= (allowed_time):
        logger.info("RESULT: performance_check PASSED!")
        if total_time_cur < total_time_ref * new_reference_threshold:
            logger.info(
                "The current experiment ran a lot faster than the reference. "
                + "Consider updating the reference."
            )
        sys.exit(0)
    else:
        logger.info("RESULT: performance_check FAILED")
        sys.exit(1)
