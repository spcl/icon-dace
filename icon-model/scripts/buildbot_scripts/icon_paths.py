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

from pathlib import Path

buildbot_script_path = Path(__file__).parent.absolute()
buildbot_list_path = (buildbot_script_path / "experiment_lists").absolute()
base_path = buildbot_script_path.parents[1].absolute()
run_path = (base_path / "run").absolute()
