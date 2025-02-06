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

import json
import signal
import sys
import time
from pathlib import Path

import icon_env
from buildbot_config import BuildbotConfig
from cmd_line_job import CmdLineJob
from icon_paths import run_path, buildbot_list_path, base_path
from pbs_job import PBSJob
from slurm_job import SlurmJob


def runexp():
    # Get build settings
    setup_dict = icon_env.info()

    builder = "{}_{}".format(setup_dict["use_target"].upper(), setup_dict["use_compiler"])

    with open(str(run_path / "buildbot_config.json"), "r") as file:
        config = json.load(file)

    builder = config["builder"]
    list_name = config["list_name"]

    print(f"initializing experiments for builder {builder} in list {list_name}")

    full_list_name = buildbot_list_path / list_name

    if Path(full_list_name).exists():
        thisList = BuildbotConfig.from_pickle(full_list_name)
    else:
        print(f"did not find experiment list {full_list_name}")
        sys.exit(1)

    bobj = thisList.builder_meta[builder]
    bobj.submit = setup_dict["use_submit"]

    # Set up mkexp run-time environment (currently breaks on DWD and Balfrin)
    if bobj.machine not in ['dwd_nec', 'balfrin']:
        icon_env.load()

    if bobj.submit.startswith("sbatch"):
        BatchJob = SlurmJob
    elif bobj.submit.startswith("qsub"):
        BatchJob = PBSJob
    elif bobj.submit == "":
        BatchJob = CmdLineJob
    else:
        print(f"no batch job helper class implemented for machine {bobj.machine}")
        sys.exit(1)

    # set up the experiments for job submission
    # builder info and batch system information need to be added
    eobjs = thisList.get_experiments_by_builder(builder).flatten()

    for eobj in eobjs:
        # BatchJob keeps track of dependencies and is responsible for submitting to the queue
        eobj.batch_job = BatchJob(bobj.submit, run_path)

    for eobj in eobjs:
        # finally submit the jobs
        eobj.submit(bobj)
        exp_jobid = eobj.batch_job.jobid
        if exp_jobid is not None:
            exp_file = eobj.get_run_name(relative=True)
            print(f"Submitted {exp_file} (jobID:{exp_jobid})")

    def cleanupOnInterupt(signum,frame):
        print("runexp: caught signal",signum, frame)
        print("runexp: Try to cancel all jobs ...")
        for eobj in eobjs:
            eobj.cancel()
    signal.signal(signal.SIGTERM, cleanupOnInterupt)

    # wait for completion of all gathered processes and do some output while waiting
    incomplete_jobs = set(eobjs)
    while len(incomplete_jobs) > 0:
        completed_jobs = set()
        for eobj in incomplete_jobs:
            if eobj.batch_job.poll(timeout=60):
                _status = f"{'OK' if eobj.batch_job.returncode == 0 else 'FAILED'}"
                print(f"job finished: {eobj.get_run_name(relative=True)} (jobID:{eobj.batch_job.jobid}): {_status}")
                completed_jobs.add(eobj)
            else:
                N = len(incomplete_jobs) - len(completed_jobs)
                if N > 1:
                    print(f"Waiting for {N} jobs:")
                    for e in (set(incomplete_jobs) - set(completed_jobs)):
                        print(f"\t{e.get_run_name(relative=True)} ({str(e.batch_job.jobid)})")
                elif N == 1:
                    e = [*incomplete_jobs][0]
                    print(f"Waiting for the last job: {e.get_run_name(relative=True)} ({e.batch_job.jobid})")
        incomplete_jobs.difference_update(completed_jobs)

    # wait a bit longer, because queing or file system might create logfiles
    # with a certain delay
    time.sleep(20)

    # 1) parse returncodes of all experiments and fail if at least one of them did fail
    # 2) write STATUS file which is supposed to be shown by buildbot for better overview
    loop_status_exp_file = Path(base_path) / "LOOP_STATUS_EXP_FILE"
    experimentsFailed    = []
    exit_code            = -1
    with open(str(loop_status_exp_file), "w") as file:
        for eobj in eobjs:
            exp_jobid = eobj.batch_job.jobid
            exp_file = eobj.get_run_name(relative=False)
            if eobj.batch_job.returncode is not None:
                if eobj.batch_job.returncode == 0:
                    experimentsFailed.append(False)
                    exit_msg = f"OK (jobID:{exp_jobid})"
                else:
                    experimentsFailed.append(True)
                    exit_msg = f"FAILED (jobID:{exp_jobid})"
            else:
                experimentsFailed.append(True)
                print(f"internal error, did not get returncode for {eobj.name} (jobID:{exp_jobid})")
                exit_msg = "FAILED"

            file.write("{0:50}: {1}\n".format(exp_file.name, exit_msg))

    if True in experimentsFailed:
        exit_code = 1
    else:
        exit_code = 0
    sys.exit(exit_code)

if __name__ == "__main__":
    runexp()
