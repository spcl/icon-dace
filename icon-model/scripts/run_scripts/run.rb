#!/usr/bin/env ruby

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

require './experiments'
include TestCases

exp = Experiment.new('vertival_mixing',
                     TestCase::TC33,
                     'R2B04',
                     [20,30,50,100,100,200,200,300,500,1000,1000,1000])

exp.run('squall')
