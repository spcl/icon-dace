from pathlib import Path

import click

from util.click_util import CommaSeperatedInts, cli_help
from util.dataframe_ops import df_from_file_ids
from util.log_handler import logger


def monotonically_increasing(L):
    return all(x <= y for x, y in zip(L[:-1], L[1:]))


def create_stats_dataframe(input_dir, file_id, stats_file_name, file_specification):
    df = df_from_file_ids(file_id, input_dir, file_specification)

    logger.info("writing stats file to {}".format(stats_file_name))

    Path(stats_file_name).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(stats_file_name)

    return df


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

    # compute stats for the ensemble run
    if ensemble:
        if len(member_num) == 1:
            member_num = [i for i in range(1, member_num[0] + 1)]
        for m_num in member_num:
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

    # compute the stats for the reference.
    # For ensembles, this file is named
    # stats_{member_id} -> stats_ref (used again in "tolerance")
    create_stats_dataframe(
        model_output_dir,
        file_id,
        stats_file_name.format(member_id="ref"),
        file_specification,
    )
