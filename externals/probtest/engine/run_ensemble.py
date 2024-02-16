import os
import re
import subprocess
import time
from pathlib import Path

import click

from util.click_util import CommaSeperatedInts, CommaSeperatedStrings, cli_help
from util.log_handler import logger
from util.utils import get_seed_from_member_num


def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


# replace an assignment statement (left=*right* to left=new)
# "right" is used here to identify multiple assignments of "left"
def replace_assignment(line, left, right_new, right_old, seed):
    if line.strip().startswith("#") or (left not in line):
        return line

    # make sure this is a regular assign statement (X = Y)
    split = line.split("=")
    if len(split) != 2:
        return line

    lhs, rhs = split
    lhs = lhs.strip()
    rhs = rhs.strip()

    if right_old and (right_old not in rhs):
        return line

    if lhs != left:
        return line

    right_new = right_new.format(seed=seed)

    out_line = "{}={}\n".format(left, right_new)
    logger.info("generating new line: {}".format(out_line.replace("\n", "")))
    return out_line


# replace all strings matching "old" substring by "new" substring
def replace_string(line, old, new):
    out_line = re.sub(old, new, line) if old in line else line
    return out_line


def prepare_perturbed_run_script(
    runscript,
    perturbed_runscript,
    experiment_name,
    modified_experiment_name,
    lhs,
    rhs_new,
    rhs_old,
    seed,
):
    in_file = open(runscript, "r")
    out_file = open(perturbed_runscript, "w")

    if rhs_old is None:
        rhs_old = [None] * len(rhs_new)
    if 1 == len(rhs_old) < len(rhs_new) and rhs_old[0] == "None":
        rhs_old = [None] * len(rhs_new)
    rhs_old = [None if r == "None" else r for r in rhs_old]

    for line in in_file:
        # replace input directory with the ones given in config file
        for lh, rh_old, rh_new in zip(lhs, rhs_old, rhs_new):
            out_line = replace_assignment(line, lh, rh_new, rh_old, seed)
            # replace first match
            if out_line != line:
                break

        # rename the experiment name
        if line == out_line:
            out_line = replace_string(line, experiment_name, modified_experiment_name)

        out_file.write(out_line)

    logger.info("writing model run script to: {}".format(perturbed_runscript))
    out_file.close()
    in_file.close()

    return


def append_job(job, job_list, parallel):
    p = subprocess.Popen(job)
    if not parallel:
        p.communicate()
        time.sleep(5)
        test_job_returncode(p)
    else:
        job_list.append(p)


def finalize_jobs(job_list, dry, parallel):
    if parallel and not dry:
        last_exception = None
        for job in job_list:
            job.communicate()
            try:
                test_job_returncode(job)
            except subprocess.CalledProcessError as e:
                logger.error(e)
                last_exception = e
        if last_exception:
            raise last_exception


def test_job_returncode(job):
    """Test job return code."""
    if job.returncode != 0:
        raise subprocess.CalledProcessError(returncode=job.returncode, cmd=job.args)


@click.command()
@click.option(
    "--run-dir",
    help=cli_help["run_dir"],
)
@click.option(
    "--perturbed-run-dir",
    help=cli_help["perturbed_run_dir"],
)
@click.option(
    "--run-script-name",
    help=cli_help["run_script_name"],
)
@click.option(
    "--perturbed-run-script-name",
    help=cli_help["perturbed_run_script_name"],
)
@click.option(
    "--experiment-name",
    help=cli_help["experiment_name"],
)
@click.option(
    "--perturbed-experiment-name",
    help=cli_help["perturbed_run_script_name"],
)
@click.option(
    "--submit-command",
    help=cli_help["submit_command"],
)
@click.option(
    "--member-num",
    default="10",
    type=CommaSeperatedInts(),
    help=cli_help["member_num"],
)
@click.option(
    "--member-type",
    type=str,
    default="",
    help=cli_help["member_type"],
)
@click.option(
    "--parallel/--no-parallel",
    is_flag=True,
    help=cli_help["parallel"],
)
@click.option(
    "--dry/--no-dry",
    is_flag=True,
    help=cli_help["dry"],
)
@click.option(
    "--lhs",
    type=CommaSeperatedStrings(),
    help=cli_help["lhs"],
)
@click.option(
    "--rhs-new",
    type=CommaSeperatedStrings(),
    help=cli_help["rhs_new"],
)
@click.option(
    "--rhs-old",
    type=CommaSeperatedStrings(),
    help=cli_help["rhs_old"],
)
def run_ensemble(
    run_dir,
    perturbed_run_dir,
    run_script_name,
    perturbed_run_script_name,
    experiment_name,
    perturbed_experiment_name,
    submit_command,
    member_num,
    member_type,
    parallel,
    dry,
    lhs,
    rhs_new,
    rhs_old,
):
    perturbed_run_dir = perturbed_run_dir if perturbed_run_dir else run_dir
    os.chdir(run_dir)
    job_list = []

    # run the reference
    if not dry:
        logger.info("running unperturbed reference")
        job = submit_command.split() + [run_script_name]
        append_job(job, job_list, parallel)

    # run the ensemble
    Path(perturbed_run_dir).mkdir(exist_ok=True, parents=True)
    os.chdir(perturbed_run_dir)
    if len(member_num) == 1:
        member_num = [i for i in range(1, member_num[0] + 1)]
    for m_num in member_num:
        m_id = str(m_num)
        if member_type:
            m_id = member_type + "_" + m_id
        runscript = "{}/{}".format(run_dir, run_script_name)
        perturbed_runscript = "{}/{}".format(
            perturbed_run_dir, perturbed_run_script_name.format(member_id=m_id)
        )

        prepare_perturbed_run_script(
            runscript,
            perturbed_runscript,
            experiment_name,
            perturbed_experiment_name.format(member_id=m_id),
            lhs,
            rhs_new,
            rhs_old,
            get_seed_from_member_num(m_num),
        )

        if not dry:
            job = submit_command.split() + [
                perturbed_run_script_name.format(member_id=m_id)
            ]
            logger.info("running the model with '{}'".format(" ".join(job)))
            append_job(job, job_list, parallel)

    finalize_jobs(job_list, dry, parallel)

    logger.info("model finished!")
