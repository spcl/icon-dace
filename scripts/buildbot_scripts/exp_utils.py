#! /usr/bin/env python3

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

import sys

from model_paths import *
from buildbot_config import BuildbotConfig
from icon_paths import buildbot_list_path
from pathlib import Path


def rmexp(experiment_names, builders, with_configureflags, without_configureflags, machines, list_name):

    full_list_name = buildbot_list_path / list_name

    if Path(full_list_name).exists():
      thisList = BuildbotConfig.from_pickle(full_list_name)
    else:
      print("did not find experiment list {}".format(full_list_name))
      sys.exit(1)

    thisList.remove_experiments(experiment_names, builders=builders, machines=machines, with_config=with_configureflags, without_config=without_configureflags)

    thisList.to_pickle(full_list_name)

def addexp(experiment_names, builders, with_configureflags, without_configureflags, machines, runflags, list_name):

    full_list_name = buildbot_list_path / list_name

    if Path(full_list_name).exists():
      thisList = BuildbotConfig.from_pickle(full_list_name)
    else:
      print("did not find experiment list {}".format(full_list_name))
      sys.exit(1)

    thisList.add_experiments(experiment_names, builders=builders, machines=machines, with_config=with_configureflags, without_config=without_configureflags, run_flags=runflags)

    thisList.to_pickle(full_list_name)

def adddep(from_builder, from_experiment, to_builder, to_experiment, builders, with_configureflags, without_configureflags, machines, list_name):

    full_list_name = buildbot_list_path / list_name

    if Path(full_list_name).exists():
      thisList = BuildbotConfig.from_pickle(full_list_name)
    else:
      print("did not find experiment list {}".format(full_list_name))
      sys.exit(1)

    thisList.add_dependency_manager(from_experiment, to_experiment, from_builder, to_builder, builders, machines, with_configureflags, without_configureflags)

    thisList.to_pickle(full_list_name)
