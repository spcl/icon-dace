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

from builder import Builder
from experiment import Experiment
from buildbot_config import BuildbotConfig
import pandas as pd

def test():
    bc = BuildbotConfig()
    bc.add_machine("daint")
    bc.add_builders("daint_cpu", "daint", "blabla", config={"a": 2, "b": 3})
    bc.add_builders("daint_gpu", "daint", "blabla", config={"c": 4, "d": 5})
    bc.add_machine("mistral")
    bc.add_builders("mistral_cpu", "mistral", "blabla", config={"a": 2, "b": 3})
    bc.add_builders("mistral_gpu", "mistral", "blabla", config={"c": 4, "d": 5})

    bc.add_experiment("mch_ch_lowres")
    bc.add_experiment("mch_ch_lowres", builders=["daint_cpu"], run_flags="--node 3")
    bc.add_experiment("atm_amip", run_flags="--node 90", with_config={"a": 2})
    bc.add_experiment("atm_amip", run_flags="--node 80", without_config={"b": 3})

    bc.add_experiment("atm_ape", machines=["mistral"])

    bc.add_builders(["mistral_intel", "mistral_gcc"], "mistral", "blabla")

    bc.add_experiment("atm_amip", run_flags="--node 80", without_config={"b": 3})

    print(bc.data)

    return

if __name__ == "__main__":
    test()