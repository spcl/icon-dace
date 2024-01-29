# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

from experiment import Experiment
from icon_paths import run_path

from pathlib import Path

def test_make_runscript():
    target = "check.atm_amip.run"
    p = run_path / Path(target)
    if p.exists(): p.unlink()

    exp = Experiment("checksuite.icon-dev/check.atm_amip", "dummy_builder")

    status = exp.make_runscript()

    assert status == 0, "make_target_runscript command failed"
    assert exp.run_name == target, "experiment filename is wrong"
    assert (run_path / Path(exp.run_name)).exists(), "experiment file has not been created"

def test_make_runscript_hardcoded():
    target = "exp.run_ICON_02_R2B13_lam.run"
    p = run_path / Path(target)
    if p.exists(): p.unlink()

    exp = Experiment("checksuite.xce.dwd.de/exp.run_ICON_02_R2B13_lam.run", "dummy_builder")

    status = exp.make_runscript()

    assert status == 0, "make_target_runscript command failed"
    assert exp.run_name == target, "experiment filename is wrong"
    assert (run_path / Path(exp.run_name)).exists(), "experiment file has not been created"
    assert (run_path / Path(exp.run_name)).is_symlink(), "hardcoded experiment file is not symlink"
