import json
import pathlib
from pathlib import Path

import click
from jinja2 import Environment, FileSystemLoader, StrictUndefined

from util.click_util import CommaSeperatedInts, cli_help
from util.log_handler import logger


@click.command()
@click.option("--codebase-install", help=cli_help["codebase_install"], default="")
@click.option(
    "--template-name",
    default=pathlib.Path(__file__).parent.parent.absolute()
    / "templates"
    / "ICON.jinja",
    help=cli_help["template_name"],
)
@click.option("--experiment-name", help=cli_help["experiment_name"], default="")
@click.option(
    "--file-id",
    nargs=2,
    type=str,
    multiple=True,
    metavar="FILE_TYPE FILE_PATTERN",
    help=cli_help["file_id"],
)
@click.option("--reference", help=cli_help["reference"], default="")
@click.option(
    "--config",
    default="probtest.json",
    help=cli_help["config"],
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
    "--perturb-amplitude",
    type=float,
    default=1e-14,
    help=cli_help["perturb_amplitude"],
)
@click.option(
    "--timing-current",
    default="",
    help=cli_help["timing_current"],
)
@click.option(
    "--timing-reference",
    default="",
    help=cli_help["timing_reference"],
)
@click.option(
    "--append-time",
    type=bool,
    default=False,
    help=cli_help["append_time"],
)
def init(
    codebase_install,
    experiment_name,
    file_id,
    reference,
    config,
    template_name,
    member_num,
    member_type,
    perturb_amplitude,
    timing_current,
    timing_reference,
    append_time,
):
    template_partition = str(template_name).rpartition("/")
    env = Environment(
        loader=FileSystemLoader(template_partition[0]), undefined=StrictUndefined
    )
    template = env.get_template(template_partition[2])
    # The template is supposed to be a valid json file so that in can work as
    # a default PROBTEST_CONFIG (even without running init first)

    # emit warnings if variables are not set
    warn_template = "init argument '--{}' not set. default to '{}'"
    if not codebase_install:
        logger.warning(warn_template.format("codebase_install", ""))
    if not experiment_name:
        logger.warning(warn_template.format("experiment_name", ""))
    if not file_id:
        logger.warning(warn_template.format("file_id", ""))
    if not reference:
        logger.warning(warn_template.format("reference", ""))
    if not member_num:
        logger.warning(warn_template.format("member_num", member_num))
    if not member_type:
        logger.warning(warn_template.format("member_type", member_type))
    if not perturb_amplitude:
        logger.warning(warn_template.format("perturb_amplitude", perturb_amplitude))
    if not timing_current:
        logger.warning(warn_template.format("timing_current", timing_current))
    if not timing_reference:
        logger.warning(warn_template.format("timing_reference", timing_reference))
    if not append_time:
        logger.warning(warn_template.format("append_time", append_time))

    # compose render dictionary
    render_dict = {}
    render_dict["experiment_name"] = experiment_name
    render_dict["codebase_install"] = Path(codebase_install).resolve()
    render_dict["reference"] = Path(reference).resolve()
    render_dict["member_num"] = member_num
    render_dict["member_type"] = member_type
    render_dict["perturb_amplitude"] = perturb_amplitude
    render_dict["timing_current"] = timing_current
    render_dict["timing_reference"] = timing_reference
    render_dict["append_time"] = append_time

    # render jinja
    rendered = template.render(render_dict)
    # append file_id via json
    json_dict = json.loads(rendered)
    json_dict["default"]["file_id"] = file_id
    json_dict["default"]["member_num"] = member_num
    rendered = json.dumps(json_dict, indent=2)
    # print file
    with open(config, "w") as probtest_config:
        probtest_config.write(rendered)

    print("Successfully wrote probtest configuration to " + config)
    return 0
