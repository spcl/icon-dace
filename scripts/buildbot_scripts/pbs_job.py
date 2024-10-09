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

import subprocess
import re
import sys

from batch_job import BatchJob

debugOutput = False

class PBSJob(BatchJob):
    def __init__(self, cmd, cwd):
        super().__init__(cmd, cwd)
        self.system = "PBS"

    def submit(self, script):
        submit_cmd = self.cmd.split()

        if len(self.parents) > 0:
            parent_ids = [p.jobid for p in self.parents]
            submit_cmd += ["--after", ",".join(parent_ids)]

        submit_cmd.append(script)
        print(f"submitting PBS job: '{submit_cmd}'")
        qsub = subprocess.Popen(submit_cmd,
                                shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                cwd=self.cwd,
                                encoding="UTF-8")

        try:
            stdout = qsub.stdout.readlines()[0]
            if debugOutput:
                print(f"|qsub: stdout = {stdout}|")
            self.jobid = re.findall(r"(\d+)\.\w+", stdout)[0]
            if debugOutput:
                print(f"|qsub: jobId = {self.jobid}|")
        except:
            for line in qsub.stderr.readlines():
                print(f"|qsub STDERR: {line}|")

        if not self.jobid:
            print(f"Parsing jobid from pbs job failed, got {self.jobid}")
            sys.exit(1)

        # qsub will return immediately. That's why we pass along qwait as the
        # real job. It's results will be handled by the wait() method
        qwaitCmd = f'qwait {self.jobid}'
        if debugOutput:
            print(f'|qwait call: "{qwaitCmd}"|')
        self.job = subprocess.Popen(qwaitCmd,
                                    shell=True,
                                    stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    encoding="utf-8")

    def poll(self, timeout):
        # PBS does not return the exitcode of the submitted script. instead it prints it to stdout
        try:
            stdout, stderr = self.job.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            return False

        if debugOutput:
            print(f"|qwait:stdout: {stdout}|")
            print(f"|qwait:stderr: {stderr}|")

        try:
            exit_code = int(stderr.split()[-1])
        except:
            print('pbs_job.py: Could not get a proper return value from "qwait"')
            print(f"qwait: stderr = |{stderr}|")
            exit_code = 1

        if debugOutput:
            print(f"qwait: exit_code = |{exit_code}|")

        self.returncode = exit_code

        return True # poll successfull, job finished.

    def cancel(self):
        if None is not self.jobid:
            subprocess.Popen(["qdel",self.jobid])
            qdel = subprocess.Popen(["qdel",self.jobid],
                                    shell=False,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    cwd=self.cwd,
                                    encoding="UTF-8")
            print(qdel.stdout.realines())
            print(qdel.stderr.realines())
        else:
            print('Cannot find jobid to cancel job!')
