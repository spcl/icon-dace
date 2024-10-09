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

import click
from pathlib import Path

from icon_paths import run_path

# Utilities for dictionary comparison
def relate_with_config(this, other):
    if not this: return True               # if no with_config is given, return all builders
    if not other: return False             # discard builders where no information is given
    b = [other.get(key) == item for key, item in this.items()]
    return all(b)

def relate_without_config(this, other):
    if not this: return True               # if no without_config is given, return all builders
    if not other: return True              # builders without config information are fine
    b = [other.get(key) == item for key, item in this.items() if key in other]
    return not any(b)

def config_dict_to_string(d):
    if not d: return ""
    return " ".join(config_dict_to_list(d))

def config_dict_to_list(d):
    if not d: return []
    return ["{}={}".format(key, value) if value else key for key, value in d.items()]

# param type for click: comma seperated lists
class WhitespaceSeperatedList(click.ParamType):
    name = 'whitespace-seperated-list'
    def convert(self, value, param, ctx):
        if not isinstance(value, str):
            self.fail(
                "Input must be a string, found {}".format(value),
                param,
                ctx
            )
        return list(filter(lambda x: x!="", value.split(" ")))

# param type for click: comma seperated lists
class WhitespaceSeperatedFiles(click.ParamType):
    name = 'whitespace-seperated-regex-list'
    def convert(self, value, param, ctx):
        if not isinstance(value, str):
            self.fail(
                "Input must be a string, found {}".format(value),
                param,
                ctx
            )
        args = list(filter(lambda x: x!="", value.split(" ")))

        regex_args = [arg for arg in args if "*" in arg]
        args = [arg for arg in args if arg not in regex_args]

        # check if there are experiments that do not exist
        not_found = [arg for arg in args if not Path("{}/{}".format(run_path, arg)).exists()]

        for f in not_found:
            print("did not find {} in {}".format(f, run_path))

        # not found
        args = [arg for arg in args if arg not in not_found]

        for ra in regex_args:
            parent = Path(ra).parent.name
            parent = parent + "/" if parent != "" else ""
            args += ["{}{}".format(parent, p.name) for p in Path(run_path).glob(ra)]

        return args

# param type for click: comma seperated lists
class WhitespaceSeperatedDict(click.ParamType):
    name = 'whitespace-seperated-dict'
    def convert(self, value, param, ctx):
        if not isinstance(value, str):
            self.fail(
                "Input must be a string, found {}".format(value),
                param,
                ctx
            )
        args = list(filter(lambda x: x!="", value.split(" ")))
        out = {}
        for arg in args:
            if len(arg.split("=")) == 2:
                key, val = arg.split("=")
                out[key] = val
            elif len(arg.split("=")) == 1:
                out[arg] = None
            else:
                self.fail(
                    "Input is not a 1 to 1 assignment: {}".format(arg),
                    param,
                    ctx
                )
        return out