# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

'''Provide setup info and load environment for ICON runs'''

import os
import sys

from icon_paths import run_path

def _parse_line(line):
    return line.strip("\n").replace("'","").split("=", 1)

_info = None

def info():
    global _info
    if _info is None:
        with open(str(run_path / "set-up.info"), "r") as f:
            _info = {k: v for k, v in [_parse_line(l) for l in f.readlines()]}
    return _info

_env_loaded = False

def load():
    global _env_loaded
    if not _env_loaded:
        if 'MODULESHOME' in os.environ:
            sys.path.insert(1, os.path.join(os.environ['MODULESHOME'], 'init'))
            try:
                from python import module # Environment Modules
            except ModuleNotFoundError:
                from env_modules_python import module  # lmod
            module('unload', 'python3') # recover module hidden by conda env
            module('load', *info().get('use_load_modules', '').split())
            sys.path.pop(1)
            _env_loaded = True
