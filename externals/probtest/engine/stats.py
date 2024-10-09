"""
CLI for computing stats

This command line tool provides functionality for:
- Creating and saving statistics dataframes from specified model output files.
- Verifying that lists of values are monotonically increasing.
- Generating statistics for both ensemble and reference model runs.
"""

import os
from multiprocessing import Pool
from pathlib import Path

import click

from util.click_util import CommaSeperatedInts, cli_help
from util.dataframe_ops import df_from_file_ids
from util.log_handler import logger


def monotonically_increasing(li):
    return all(x <= y for x, y in zip(li[:-1], li[1:]))


def create_stats_dataframe(input_dir, file_id, stats_file_name, file_specification):
    df = df_from_file_ids(file_id, input_dir, file_specification)

    logger.info("writing stats file to %s", stats_file_name)

    Path(stats_file_name).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(stats_file_name)

    return df


def process_member(
    m_num,
    member_type,
    model_output_dir,
    perturbed_model_output_dir,
    file_id,
    stats_file_name,
    file_specification,
):
    if m_num == 0:
        input_dir = model_output_dir
        m_id = "ref"
    else:
        m_id = str(m_num)
        if member_type:
            m_id = member_type + "_" + m_id
        input_dir = perturbed_model_output_dir.format(member_id=m_id)
    create_stats_dataframe(
        input_dir,
        file_id,
        stats_file_name.format(member_id=m_id),
        file_specification,
    )


@click.command()
@click.option(
    "--ensemble/--no-ensemble",
    is_flag=True,
    help=cli_help["ensemble"],
)
@click.option(
    "--stats-file-name",
    help=cli_help["stats_file_name"],
)
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
    "--file-specification",
    type=list,
    help=cli_help["file_specification"],
)
def stats(
    ensemble,
    stats_file_name,
    model_output_dir,
    file_id,
    member_num,
    member_type,
    perturbed_model_output_dir,
    file_specification,
):
    file_specification = file_specification[0]  # can't store dicts as defaults in click
    assert isinstance(file_specification, dict), "must be dict"

    # compute stats for the ensemble and the reference run
    if ensemble:
        # Add 0 to the list of member numbers to include the reference run
        if len(member_num) == 1:
            member_num = list(range(member_num[0] + 1))
        else:
            member_num.append(0)
        with Pool() as p:
            args = [
                (
                    m_num,
                    member_type,
                    model_output_dir,
                    perturbed_model_output_dir,
                    file_id,
                    stats_file_name,
                    file_specification,
                )
                for m_num in member_num
            ]
            p.starmap(process_member, args)
    else:
        create_stats_dataframe(
            model_output_dir,
            file_id,
            stats_file_name.format(member_id=os.path.basename(model_output_dir)),
            file_specification,
        )
