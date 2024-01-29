# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

import subprocess
import sys
import pathlib

from batch_job import BatchJob

class CmdLineJob(BatchJob):
    def __init__(self, cmd, cwd):
        super().__init__(cmd, cwd)
        self.system = "Commandline"

    # Command line jobs cannot be submitted in parallel - wait at submission
    def wait(self):
        print("Command line jobs run sequentially. Continue.")

    def poll(self, timeout):
        print("Command line jobs run sequentially. Nothing to poll.")
        return True

    def submit(self, script):
        if len(self.parents) > 0:
            print("Dependencies are not supported for {}-jobs".format(self.system))
            sys.exit(1)

        # store output as LOG-file following the buildbot nameing convention
        full_cmd = "./{}".format(script)
        with open("{}/LOG.{}.o".format(self.cwd, script), "wb") as out:
            self.job = subprocess.Popen(full_cmd, shell=False, stdout=out, stderr=out, cwd=self.cwd, encoding="UTF-8")
            # wait for job to finish before starting next one
            self.returncode = self.job.wait()
        # Add a symlink following the same naming conventions, but for non-generated runscrits
        pathlib.Path("{}/LOG.{}.run.o".format(self.cwd, script)).symlink_to(pathlib.Path("{}/LOG.{}.o".format(self.cwd, script)))
