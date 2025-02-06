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

from buildbot_config import BuildbotConfig
import pytest

@pytest.fixture
def base_bc():
    bc = BuildbotConfig()

    bc.add_machine("daint")
    bc.add_builders("daint_cpu", "daint", "dummy_config_wrapper", config={"a": 2, "b": 3})
    bc.add_builders("daint_gpu", "daint", "dummy_config_wrapper", config={"c": 4, "d": 5})
    bc.add_machine("mistral")
    bc.add_builders("mistral_cpu", "mistral", "dummy_config_wrapper", config={"a": 2, "b": 3})
    bc.add_builders("mistral_gpu", "mistral", "dummy_config_wrapper", config={"c": 4, "d": 5})

    bc.add_experiment("mch_ch_lowres")
    bc.add_experiment("mch_ch_lowres", builders=["daint_cpu"], run_flags="--node 3")
    bc.add_experiment("atm_amip", run_flags="--node 90", with_config={"a": 2})
    bc.add_experiment("atm_amip", run_flags="--node 80", without_config={"b": 3})

    bc.add_experiment("atm_ape", machines=["mistral"])

    bc.add_builders(["mistral_intel", "mistral_gcc"], "mistral", "dummy_config_wrapper")

    bc.add_experiment("atm_amip", run_flags="--node 80", without_config={"b": 3})
    bc.add_experiment("mch_opr")

    bc.add_builders("dwd_nec", "nec", "blubb", flag="Inactive")

    return bc

def test_rmexp_machines(base_bc):
    base_bc.remove_experiment("mch_opr", machines=["daint"])

    assert "mch_opr" not in base_bc.data.index, "Found Experiment after 'remove_experiments'"

def test_set_builder_flags(base_bc):
    base_bc.set_builder_flag("dwd_nec", "Active")

    assert base_bc.builder_meta["dwd_nec"].flag == "Active", "Builder flag has not been set correctly"
