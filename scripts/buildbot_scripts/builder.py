# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

from util import config_dict_to_list, config_dict_to_string
from icon_paths import base_path
import subprocess

class Builder(object):
    def __init__(self, name, machine, script, config, flag):
        self.name = name
        self.machine = machine
        self.script = script
        self.config = config
        self.flag = flag

        # a builders submit command is only available once set-up.info has been written
        # this is after the configure step. The submit-command will be added just in time.
        self.submit = None

    def to_string(self):
        out = "{} ({}):".format(self.name, self.flag)
        if self.config: out += "{}".format(config_dict_to_string(self.config))
        return out

    def configure_make(self):
        # do nothing if builder is inactive
        if self.flag == "Inactive": return 0

        # execute configure wrapper (note that for buildbot that contains the "make" step)
        if not self.script: 
            print("Configure wrapper missing for builder {}. Check your create_list_<list> file.".format(self.name))
            return 1
        cmd = [self.script] + config_dict_to_list(self.config)

        sp = subprocess.run(cmd, shell=False, cwd=base_path, stderr=subprocess.STDOUT, universal_newlines=True)

        return sp.returncode
