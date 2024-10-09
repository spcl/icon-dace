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

# -*- coding: utf-8 -*-

import pkgutil
li = list(pkgutil.iter_modules())
print(li)

import pandas as pd

data = pd.DataFrame()

print("successfully instantiated a pandas DataFrame!")
