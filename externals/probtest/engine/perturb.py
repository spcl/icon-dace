import os
import shutil

import click
import numpy as np

from util.click_util import CommaSeperatedInts, CommaSeperatedStrings, cli_help
from util.log_handler import logger
from util.netcdf_io import nc4_get_copy
from util.utils import get_seed_from_member_num


def create_perturb_files(in_path, in_files, out_path, copy_all_files=False):
    path = os.path.abspath(in_path)
    if not os.path.exists(out_path):
        logger.info("creating new directory: {}".format(out_path))
        os.makedirs(out_path)
    data = [
        nc4_get_copy("{}/{}".format(path, f), "{}/{}".format(out_path, f))
        for f in in_files
    ]

    if copy_all_files:
        all_files = os.listdir(path)
        # disregard the input files which are copied above
        other_files = [f for f in all_files if f not in in_files]
        # copy all other files
        for f in other_files:
            shutil.copy("{}/{}".format(in_path, f), out_path)

    return data


def perturb_array(array, s, a):
    shape = array.shape
    np.random.seed(s)
    p = (
        np.random.rand(*shape) * 2 - 1
    ) * a + 1  # *2-1: map to [-1,1), *a: rescale to amplitude, +1 perturb around 1
    parray = np.copy(array * p)
    return parray


@click.command()
@click.option("--experiment-name", help=cli_help["experiment_name"], default="")
@click.option(
    "--model-input-dir",
    help=cli_help["model_input_dir"],
)
@click.option(
    "--perturbed-model-input-dir",
    help=cli_help["perturbed_model_input_dir"],
)
@click.option(
    "--files",
    type=CommaSeperatedStrings(),
    help=cli_help["files"],
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
    "--variable-names",
    type=CommaSeperatedStrings(),
    help=cli_help["variable_names"],
)
@click.option(
    "--perturb-amplitude",
    type=float,
    help=cli_help["perturb_amplitude"],
)
@click.option(
    "--copy-all-files/--no-copy-all-files",
    is_flag=True,
    help=cli_help["copy_all_files"],
)
def perturb(
    experiment_name,
    model_input_dir,
    perturbed_model_input_dir,
    files,
    member_num,
    member_type,
    variable_names,
    perturb_amplitude,
    copy_all_files,
):
    if len(member_num) == 1:
        member_num = [i for i in range(1, member_num[0] + 1)]
    for m_num in member_num:
        m_id = str(m_num)
        if member_type:
            m_id = member_type + "_" + m_id
        perturbed_model_input_dir_member_id = perturbed_model_input_dir.format(
            member_id=m_id
        )
        data = create_perturb_files(
            model_input_dir,
            files,
            perturbed_model_input_dir_member_id,
            copy_all_files,
        )
        for d in data:
            for vn in variable_names:
                d.variables[vn][:] = perturb_array(
                    d.variables[vn][:],
                    get_seed_from_member_num(m_num, use_64_bits=False),
                    perturb_amplitude,
                )
            d.close()
