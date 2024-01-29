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
class BatchJob(object):
    def __init__(self, cmd, cwd):
        self.system = "undefined"
        self.cmd = cmd
        self.cwd = cwd
        self.jobid = None
        self.job = None
        self.returncode = None
        self.parents = []

    def add_parent(self, parent):
        self.parents.append(parent)

    def wait(self):
        self.poll(timeout=None)

    def poll(self, timeout):
        """Check if task is still running.

        Waits up to specified timeout in seconds for job to finish. If job
        finishes in time or has finished before, return True and set returncode.
        """
        try:
            returncode = self.job.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            return False
        else:
            self.returncode = returncode
            return True

    def cancel(self):
        self.job.cancel()
