"""
This module provides utility functionality for loading configuration defaults for
probtest from a JSON configuration file and defining custom Click parameter
types for handling comma-separated values.
It also contains help messages for various command-line options.
"""

import json
import os
import pathlib

import click

from util import dataframe_ops
from util.log_handler import logger


def load_defaults(sections):
    configfile = os.environ.get("PROBTEST_CONFIG", "probtest.json")
    config = {}

    if not os.path.isfile(configfile):
        logger.warning(
            "careful, configfile %s does not exist in %s. No defaults.",
            configfile,
            os.getcwd(),
        )
        logger.warning(
            "If you need e.g. ICON defaults, try: export PROBTEST_CONFIG=%s",
            pathlib.Path(__file__).parent.parent.absolute() / "templates/ICON.jinja",
        )
    else:
        logger.info("reading config file from %s", configfile)
        with open(configfile, encoding="utf-8") as f:
            try:
                config = json.load(f)
            except json.decoder.JSONDecodeError as err:
                logger.error("failed to load configfile %s", configfile)
                logger.error(err)

    tmp = {}
    if config:
        # set default variables
        if "default" in config.keys():
            tmp.update(config.get("default"))
        # give current sections precedence
        for i in range(len(sections)):
            section = "-".join(sections[: i + 1])
            if section in config.keys():
                tmp.update(config.get(section))

    return tmp


class CommaSeperatedInts(click.ParamType):
    """
    A custom Click parameter type for handling comma-separated integers.

    This class defines a Click parameter type that converts a comma-separated
    string input into a list of integers.
    It is used to parse command-line arguments where multiple integer values are
    provided as a single comma-separated string.
    """

    name = "CommaSeperatedInts"

    def convert(self, value, param, ctx):
        if isinstance(value, list):
            return value
        if not isinstance(value, str):
            self.fail(f"Input must be a string, found {value}", param, ctx)
        return [int(e) for e in filter(lambda x: x != "", value.split(","))]


class CommaSeperatedStrings(click.ParamType):
    """
    A custom Click parameter type for handling comma-separated strings.

    This class defines a Click parameter type that converts a comma-separated
    string input into a list of strings.
    It is used to parse command-line arguments where multiple string values are
    provided as a single comma-separated string.
    """

    name = "CommaSeperatedStrings"

    def convert(self, value, param, ctx):
        if isinstance(value, list):
            return value
        if not isinstance(value, str):
            self.fail(f"Input must be a int, found {value}", param, ctx)
        return list(filter(lambda x: x != "", value.split(",")))


cli_help = {
    "model_input_dir": r"the directory where the model input is stored",
    "perturbed_model_input_dir": r"Template for the directory where the perturbed "
    + r"model input is stored. Must contain '\{member_id\}'.",
    "model_output_dir": r"the directory where the model output is stored",
    "perturbed_model_output_dir": r"Template for the directory where the perturbed "
    + r"model output is stored. Must contain '\{member_id\}''.",
    "experiment_name": r"the name of the experiment to be run",
    "perturbed_experiment_name": r"Template for the name of the perturbed experiments. "
    + r"Must contain '\{member_id\}'.",
    "tolerance_file_name": r"the name of the file containing the tolerances "
    + r"(per time step and variable)",
    "stats_file_name": r"the name of the stats file. No absolute path here, it will "
    + r"always be created in the (perturbed) model_output_dir.",
    "member_num": r"number of ensemble members or a comma separated list of members "
    + r'(e.g. "1,3,14")',
    "member_type": r"precision of experiment (e.g. double or mixed). "
    + r"The type is part of the created member_id, which is equal to "
    + r"(member_type+'_'+str(member_num))",
    "perturb_amplitude": r"the amplitude of the relative perturbation",
    "files": r"the files that need to be perturbed (comma separated list)",
    "variable_names": r"the variables that are perturbed (comma separated list)",
    "copy_all_files": r"copy all files from the model_input_dir directory",
    "file_id": r"A unique identifier and file pattern FILE_PATTERN of the files "
    + r"containing the variables to be analysed and the file specification label "
    + r"FILE_TYPE. FILE_PATTERN may contain simple shell-style wildcards such as "
    + r'"*" and will be expanded internally by glob. Put FILE_PATTERN in quotes to '
    + r"avoid early glob expansion by the calling shell.",
    "ensemble": r"For ensemble stats: the sub-directory where the ensemble outputs are",
    "test_tolerance": r"Test with how many stats files out "
    + r"of x (x=total_member_num) the probtest passes given a specific tolerance file.",
    "total_member_num": r"Number of ensemble members used "
    + r"to select the x members from (min_member_num<=x<=max_member_num).",
    "min_member_num": r"Minimum number of members to select.",
    "max_member_num": r"Maximum number of members to select.",
    "min_factor": r"Minimum tolerance factor used to select members with.",
    "max_factor": r"Maximum tolerance factor used to select members with.",
    "iterations": r"Maximum number of iterations to select a random sample per "
    + r"number of members.",
    "selected_members_file_name": r"The name of the file in which the selected "
    + r"members and the used factor will be written.",
    "file_specification": "Specify how different file types shall be read. This "
    + r"option must be defined in the json config file. See doc string of  "
    + r"df_from_file_ids for the specification.",
    "input_file_ref": r"reference file to check against",
    "input_file_cur": r"current file to be tested",
    "factor": r"relaxation factor for the tolerance values",
    "timing_regex": r"regex for the file that contains the latest log",
    "timing_names": r"the name of the timing entries to be displayed (comma separated)",
    "timing_database": r"a persistent file to keep track of performance history",
    "append_time": r"if true: append to the performance data, if false: overwrite the "
    + r"performance data (default: false)",
    "run_dir": r"directory from where the run is started "
    + r"(with 'submit_command model_run_script_name')",
    "perturbed_run_dir": r"directory from where the perturbed run is started "
    + r"(with 'submit_command model_run_script_name')",
    "run_script_name": r"name of the original experiment runscript",
    "perturbed_run_script_name": r"Template for the perturbed experiment name. Must "
    + r"contain '\{member_id\}'.",
    "lhs": r"replace assignments in the runscript. For multiples, use comma separated "
    + r"list. Note that the new right hand side can depend on \{member_id\} define "
    + r"left hand side",
    "rhs_new": r"define new right hand side",
    "rhs_old": r"define old right hand side (optional, put None if not needed)",
    "submit_command": r"How a model simulation is submitted",
    "parallel": r"can the jobs run in parallel?",
    "dry": r"only generate runscripts, do not run the model",
    "savedir": r"the directory where the plots are stored",
    "cdo_table_file": r"file to store the cdo table into",
    "variables": r"select variables to print CDO diff from ensemble",
    "times": r"select times to print CDO diff from ensemble",
    "histogram": r"print out full histogram of relative differences",
    "codebase_install": r"the directory where the code base is installed",
    "reference": r"the directory where reference files are read from and written to",
    "config": r"the name of the config file that is being generated",
    "template_name": r"path to the template for the config file",
    "i_table": r"which table to plot, must be an int between 0 and n_tables-1",
    "timing_current": r"the database containing the timings for the current run",
    "timing_reference": r"the database containing the reference timings",
    "measurement_uncertainty": r"How much time in [s] the current time can deviate "
    + "from the reference.",
    "tolerance_factor": r"The factor by which the current time can deviate from the "
    + "reference.",
    "new_reference_threshold": r"The factor by which the current time can be faster "
    + r"than the reference before a warning gets printed.",
}

del dataframe_ops
