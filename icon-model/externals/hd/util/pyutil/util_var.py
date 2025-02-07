# -*- coding: utf-8 -*-
#
# util_var.py - Utilities for cariable definitions and characteristics, e.g. for plotting
# 
# Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#

"""
###############################################################################
## Utility for Variable definitions and characteristics, e.g. for plotting
#
###############################################################################
## Author: Stefan Hagemann (Hereon, Institute of Coastal Systems, GER)
##
"""


class Variable:
    number = ""

    def __init__(self, number):
        self.no = number
        self.longname = "Long variable name"
        self.name = "Variable"
        self.tag = "tag"
        self.unit = "unit"
        self.ivarbar = 0  # Set no. for colorbar setting
        self.ldiff = False  # Absolute (F) or Difference (T) values

        self.vmin = None
        self.vmax = None
        self.ncol = None
        self.cmap = None


